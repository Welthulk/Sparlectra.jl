# Copyright 2023-2026 Udo Schmitz
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# file: src/adapters/cgmes/cgmes_dynamics.jl
# purpose: the CGMES DY profile to DynamicData (#506 part 1, switch
#          cgmes_import.read_dynamics): inertia and damping of the
#          synchronous machine models, droop and governor time constant of
#          the turbine governors whose parameter meaning the specification
#          states. Every DY class gets one report line with the objects
#          mapped, not mapped and without counterpart in the EQ.

# Synchronous machine models. Each inherits RotatingMachineDynamics, whose
# attributes the CGMES 2.4.15 DY RDFS (IEC 61970-302) describes as
#   inertia: "Inertia constant of generator or motor and mechanical load (H)
#            (>0) ... per unit on the generator MVA base, usually expressed
#            as MW*second/MVA or just second"
#   damping: "Damping torque coefficient (D)", a per-unit value
# and links to its machine by SynchronousMachineDynamics.SynchronousMachine.
# Inheritance chain: SynchronousMachineSimplified -> SynchronousMachineDynamics;
# SynchronousMachineTimeConstantReactance and SynchronousMachineEquivalentCircuit
# -> SynchronousMachineDetailed -> SynchronousMachineDynamics ->
# RotatingMachineDynamics. SynchronousMachineUserDefined is left out: its
# parameters are proprietary and it is counted as not mapped.
const _DY_MACHINE_CLASSES = [:SynchronousMachineSimplified, :SynchronousMachineTimeConstantReactance, :SynchronousMachineEquivalentCircuit]

# Turbine governors with a mapped parameter, linked to their machine model by
# TurbineGovernorDynamics.SynchronousMachineDynamics. `droop`: the attribute
# that holds the permanent droop R, or with `inverse` the gain K = 1/R.
# `time`: the attribute the specification names the governor lag, `nothing`
# where the class has none (its time constants are turbine or servo stages
# of a multi-block model, not the first-order governor lag of DynamicData).
# Attribute descriptions from the CGMES 2.4.15 DY RDFS (same in CIM17):
#   GovSteam0.r      "Permanent droop (R)"; t1 "Steam bowl time constant"
#   GovHydro1.rperm  "Permanent droop (R) (>0)"; tg "Gate servo time constant"
#   GovHydro2.rperm  "Permanent droop (Rperm)"; tg "Gate servo time constant"
#   GovSteam1.k      "Governor gain (reciprocal of droop) (K) (>0)";
#            t1      "Governor lag time constant (T1)"
#   mwbase (all four) "Base for power values (MWbase) (>0). Unit = MW"
# The droop is per unit on MWbase. It is converted to the rating of the unit
# (R * S / MWbase), which needs MWbase > 0 and a known rating S; without
# either the droop is not mapped and counted.
const _DY_GOVERNORS = [
  (class = :GovSteam0, droop = :r, inverse = false, time = nothing),
  (class = :GovHydro1, droop = :rperm, inverse = false, time = nothing),
  (class = :GovHydro2, droop = :rperm, inverse = false, time = nothing),
  (class = :GovSteam1, droop = :k, inverse = true, time = :t1),
]

# Whether DY object `o` reaches an object of another profile (its EQ
# counterpart) over references in either direction, across DY objects.
# `incoming` maps a DY object to the objects that refer to it.
function _dyReachesEQ(o::CIMObject, store::CGMESStore, dy_files::Set{String}, incoming::Dict{String,Vector{String}})::Bool
  seen = Set{String}([o.mrid])
  queue = [o.mrid]
  while !isempty(queue)
    m = popfirst!(queue)
    obj = store.objects[m]
    for n in vcat(collect(values(obj.refs)), get(incoming, m, String[]))
      n in seen && continue
      push!(seen, n)
      t = get(store.objects, n, nothing)
      t === nothing && continue
      t.source in dy_files || return true
      push!(queue, n)
    end
  end
  return false
end

"""
    _mapDynamics!(net, store, ctx)

Fill `ProSumer.dynamics` from the DY profile of `store` (switch
`cgmes_import.read_dynamics`). Runs after the injections are mapped, before
the island slack choice, which reads the unit type.

A synchronous machine model (`_DY_MACHINE_CLASSES`) gives its machine the
inertia, the damping, the rating (`SynchronousMachine.ratedS` of the EQ,
not a placeholder of 9999 MVA and above) and the unit type
`SynchronousUnit`. A turbine governor of `_DY_GOVERNORS` on a mapped machine
model adds the droop, converted from MWbase to that rating, and the governor
time constant where its class has one. A function block with
`enabled = false`, a second machine model or governor on one machine, and an
object whose machine is not part of the network are not mapped. Every other
DY class is counted, not mapped. One `notice:` per DY class states the
objects mapped, not mapped and without counterpart in the EQ (no
`SynchronousMachine` behind the reference, or for the classes not mapped no
reference path to an EQ object); further notices count the reasons. A
record whose values break the `DynamicData` rules is dropped with a
`warning:`.
"""
function _mapDynamics!(net, store::CGMESStore, ctx::_MapCtx)
  dy_files = Set{String}(f.name for f in store.files if :DY in f.profiles && !f.skipped)
  if isempty(dy_files)
    push!(ctx.messages, "notice: cgmes_import.read_dynamics is on, the delivery carries no DY profile")
    return nothing
  end
  is_dy(o::CIMObject) = o.source in dy_files
  dy_objects(cls::Symbol) = [o for o in objectsOf(store, cls) if is_dy(o)]
  mapped = Dict{Symbol,Int}()
  not_mapped = Dict{Symbol,Int}()
  no_eq = Dict{Symbol,Int}()
  reasons = Dict{Tuple{Symbol,String},Int}()
  bump!(d, cls) = (d[cls] = get(d, cls, 0) + 1)
  note!(cls, why) = (reasons[(cls, why)] = get(reasons, (cls, why), 0) + 1)
  skip!(cls, why) = (bump!(not_mapped, cls); note!(cls, why))
  # the machine prosumer behind a SynchronousMachine mRID, or why there is none
  function machine_of(sm_mrid)
    sm_mrid === nothing && return (nothing, :no_eq)
    sm = get(store.objects, sm_mrid, nothing)
    (sm === nothing || sm.class != :SynchronousMachine) && return (nothing, :no_eq)
    el = get(ctx.element_by_mrid, sm_mrid, nothing)
    (el === nothing || first(el) != :generator) && return (nothing, :not_in_network)
    return (last(el), sm)
  end

  # per prosumer: the DynamicData keywords and the DY objects they came from
  acc = Dict{Int,Dict{Symbol,Any}}()
  sources = Dict{Int,Vector{Tuple{Symbol,String}}}()
  model_of = Dict{String,Int}()   # machine model mRID -> prosumer
  for cls in _DY_MACHINE_CLASSES, o in dy_objects(cls)
    if boolval(o, :enabled, true) == false
      skip!(cls, "function block disabled")
      continue
    end
    k, sm = machine_of(get(o.refs, :SynchronousMachine, nothing))
    if k === nothing
      sm === :no_eq ? bump!(no_eq, cls) : skip!(cls, "machine not part of the network")
      continue
    end
    if haskey(acc, k)
      skip!(cls, "second machine model on one machine")
      continue
    end
    ps = net.prosumpsVec[k]
    rs = num(sm, :ratedS)
    usable = rs !== nothing && rs > 0.0 && rs < Sparlectra.SC_RATEDS_PLACEHOLDER_MVA
    # the rating H and D refer to; left out when the unit carries it already
    rating = (usable && ps.ratedS != rs) ? rs : nothing
    acc[k] = Dict{Symbol,Any}(:inertia_s => num(o, :inertia), :damping_pu => num(o, :damping, 0.0), :rating_mva => rating, :unit_type => Sparlectra.SynchronousUnit)
    sources[k] = [(cls, o.mrid)]
    model_of[o.mrid] = k
  end

  governed = Set{Int}()
  for g in _DY_GOVERNORS, o in dy_objects(g.class)
    cls = g.class
    if boolval(o, :enabled, true) == false
      skip!(cls, "function block disabled")
      continue
    end
    smd_mrid = get(o.refs, :SynchronousMachineDynamics, nothing)
    smd = smd_mrid === nothing ? nothing : get(store.objects, smd_mrid, nothing)
    if smd === nothing || machine_of(get(smd.refs, :SynchronousMachine, nothing))[2] === :no_eq
      bump!(no_eq, cls)
      continue
    end
    k = get(model_of, smd_mrid, nothing)
    if k === nothing
      skip!(cls, "its machine model is not mapped")
      continue
    end
    if k in governed
      skip!(cls, "second governor on one machine")
      continue
    end
    # droop on MWbase, converted to the rating of the unit
    raw = num(o, g.droop)
    r_base = (raw === nothing || raw <= 0.0) ? nothing : (g.inverse ? 1.0 / raw : raw)
    mwbase = num(o, :mwbase)
    s_unit = something(acc[k][:rating_mva], net.prosumpsVec[k].ratedS, Some(nothing))
    droop = nothing
    if r_base === nothing
      note!(cls, "droop left out: absent or not positive")
    elseif mwbase === nothing || mwbase <= 0.0
      note!(cls, "droop left out: mwbase 0 or absent, the base of the droop is not stated")
    elseif s_unit === nothing
      note!(cls, "droop left out: the machine has no usable ratedS")
    else
      droop = r_base * s_unit / mwbase
    end
    tgov = g.time === nothing ? nothing : num(o, g.time)
    if droop === nothing && tgov === nothing
      skip!(cls, "neither a usable droop nor a governor time constant")
      continue
    end
    droop === nothing || (acc[k][:droop_pu] = droop)
    tgov === nothing || (acc[k][:governor_time_s] = tgov)
    push!(governed, k)
    push!(sources[k], (cls, o.mrid))
  end

  n_units = 0
  for k in sort!(collect(keys(acc)))
    kw = filter(p -> last(p) !== nothing, acc[k])
    dd = try
      Sparlectra.DynamicData(; kw...)
    catch e
      # the record rules of DynamicData (a negative inertia, a droop of 0);
      # the delivery's values are reported and the unit stays without data
      e isa Sparlectra.DynamicDataError || rethrow()
      push!(ctx.messages, "warning: DY data of $(Sparlectra.getCompName(net.prosumpsVec[k].comp)) not taken: $(sprint(showerror, e))")
      for (cls, _) in sources[k]
        skip!(cls, "values break the DynamicData rules")
      end
      continue
    end
    net.prosumpsVec[k].dynamics = dd
    n_units += 1
    for (cls, _) in sources[k]
      bump!(mapped, cls)
    end
  end

  # every other class the DY files created: not mapped, or without any
  # reference path to an object of another profile
  handled = Set{Symbol}(vcat(_DY_MACHINE_CLASSES, [g.class for g in _DY_GOVERNORS]))
  incoming = Dict{String,Vector{String}}()
  for (m, obj) in store.objects, t in values(obj.refs)
    t_obj = get(store.objects, t, nothing)
    (t_obj !== nothing && is_dy(t_obj)) && push!(get!(Vector{String}, incoming, t), m)
  end
  dy_classes = Set{Symbol}()
  for (cls, mrids) in store.byclass, m in mrids
    is_dy(store.objects[m]) || continue
    push!(dy_classes, cls)
    cls in handled && continue
    _dyReachesEQ(store.objects[m], store, dy_files, incoming) ? bump!(not_mapped, cls) : bump!(no_eq, cls)
  end

  push!(ctx.messages, "notice: DY profile read (cgmes_import.read_dynamics): $(n_units) unit(s) with dynamic data")
  for cls in sort!(collect(dy_classes))
    push!(ctx.messages, "notice: DY class $(cls): $(get(mapped, cls, 0)) mapped, $(get(not_mapped, cls, 0)) not mapped, $(get(no_eq, cls, 0)) without counterpart in EQ")
  end
  for ((cls, why), n) in sort!(collect(reasons))
    push!(ctx.messages, "notice: DY $(cls): $(n) object(s), $(why)")
  end
  return nothing
end
