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

# file: src/hvdc_pair_control.jl
# purpose: back-to-back HVDC pairing controller (outer loop, issue #297
#          Draft B): couples the two converter injections of an HVDC link
#          so the transfer P obeys the pairing invariant
#          P_to = P_from_transfer - losses while each terminal holds either
#          a fixed Q or a voltage target (per-side secant like the machine
#          controller). No angle coupling is introduced: two areas joined
#          only by the pair stay separate electrical islands.
#          Two modes: :setpoint (transfer is a given) and :island_feed
#          (grid-forming receiving converter IS the island reference; the
#          sending side mirrors the island draw plus the converter loss).

"""
    HvdcPairControl <: AbstractOuterController

Pairing controller for the two converter injections of a back-to-back (or
point-to-point) HVDC link. Sign convention follows the MATPOWER `dcline`:
`p_transfer_mw` is the active power LEAVING the from-side AC bus into the
link; the from-side injection is `-p_transfer_mw`, the to-side injection is
`p_transfer_mw - loss` with `loss = loss_mw + loss_fraction * |p_transfer|`.
In the default `:setpoint` mode the transfer is a control setpoint (HVDC
has no angle coupling), so the active-power side needs no iteration; only
voltage-target terminals iterate via the per-side secant on the terminal Q
within its reactive range.

In `:island_feed` mode the dependency inverts: the receiving converter is
grid-forming, it IS the reference (slack) of its island, and its output is
the island balance outcome. The controller then mirrors that outcome onto
the sending side each outer iteration, `P_from = -(P_island + loss)`, with
an honest `at_limit` once the island draw exceeds `p_rating_mw` (the power
flow's slack still balances the island, so the model cannot show the real
collapse; the flag marks the violated rating).
State is owned by `run_control!`; construct via [`addHvdcPairControl!`](@ref).
"""
mutable struct HvdcPairControl <: AbstractOuterController
  name::String
  from_bus::String
  to_bus::String
  from_prosumer::Int
  to_prosumer::Int
  p_transfer_mw::Float64
  loss_mw::Float64
  loss_fraction::Float64
  p_rating_mw::Union{Nothing,Float64}
  # per-terminal reactive handling: fixed Q or voltage target (exclusive)
  from_q_mvar::Union{Nothing,Float64}
  to_q_mvar::Union{Nothing,Float64}
  from_vset_pu::Union{Nothing,Float64}
  to_vset_pu::Union{Nothing,Float64}
  from_qmin_mvar::Float64
  from_qmax_mvar::Float64
  to_qmin_mvar::Float64
  to_qmax_mvar::Float64
  deadband_vm_pu::Float64
  max_outer_iters::Int
  enabled::Bool
  status::Symbol
  converged::Bool
  at_limit::Bool
  p_applied::Bool
  # a rating-clamped transfer means the TARGET transfer is not reachable:
  # honest semantics are at_limit = true and converged = false, exactly like
  # a reactance clamped at its range end
  transfer_clamped::Bool
  # live actuator values and secant state per side (vset mode only)
  from_q_now::Float64
  to_q_now::Float64
  from_vm_now::Union{Nothing,Float64}
  to_vm_now::Union{Nothing,Float64}
  from_prev_q::Union{Nothing,Float64}
  from_prev_vm::Union{Nothing,Float64}
  to_prev_q::Union{Nothing,Float64}
  to_prev_vm::Union{Nothing,Float64}
  outer_iters::Int
  # :setpoint (transfer given) or :island_feed (grid-forming to-side slack,
  # transfer mirrored from the island balance each outer iteration)
  mode::Symbol
  # island_feed only: the mirror is settled when the applied transfer and
  # the freshly derived one agree within this band (MW)
  deadband_p_mw::Float64
end

"""
    _hvdc_pair_controllers(net) -> Vector{HvdcPairControl}

Collect the HVDC pair controllers stored on `net` (shared registry
`net.machineControls`, `hasproperty`-guarded like the other collectors).
"""
function _hvdc_pair_controllers(net)::Vector{HvdcPairControl}
  hasproperty(net, :machineControls) || return HvdcPairControl[]
  return HvdcPairControl[c for c in net.machineControls if c isa HvdcPairControl]
end

"""
    clearHvdcPairControllers!(net)

Remove all HVDC pair controllers from `net` (other controllers stay). The
persistent link records in `net.hvdcLinks` are kept, only their
`controller_name` is reset: the link exists independently of how it is
steered.
"""
function clearHvdcPairControllers!(net::Net)
  filter!(c -> !(c isa HvdcPairControl), net.machineControls)
  for (i, l) in enumerate(net.hvdcLinks)
    l.controller_name === nothing || (net.hvdcLinks[i] = _hvdc_link_with_controller(l, nothing))
  end
  return net
end

"""
    addHvdcLink!(net; from_bus, to_bus, name = nothing, kind = :b2b,
                 status = 1, from_prosumer = nothing, to_prosumer = nothing)

Register a Stage-0 HVDC link record on a hand-built net so the result layer
reports it (`HVDC Link Flows` table, `ACPFlowReport.hvdc_links`,
`hvdc_links.csv`). Importers and [`addHvdcPairControl!`](@ref) register
links automatically; this helper covers programmatic nets that keep fixed
converter injections without a controller (including a grid-forming
reference terminal, which a setpoint pair could not carry). Terminals
resolve like the pair controller (unique generator-type injection per bus,
explicit index on ambiguity) but without role restrictions: the record is
bookkeeping, not an actuator. Returns the created [`HvdcLink`](@ref).
"""
function addHvdcLink!(net::Net; from_bus::String, to_bus::String, name::Union{Nothing,String} = nothing, kind::Symbol = :b2b, status::Int = 1, from_prosumer::Union{Nothing,Int} = nothing, to_prosumer::Union{Nothing,Int} = nothing)
  kind in (:b2b, :p2p) || error("addHvdcLink!: kind must be :b2b or :p2p, got $(kind)")
  resolve = function (bus::String, explicit::Union{Nothing,Int}, side::String)
    busIdx = geNetBusIdx(net = net, busName = bus)
    explicit === nothing || begin
      (1 <= explicit <= length(net.prosumpsVec)) || error("addHvdcLink!: $(side)_prosumer $(explicit) out of range")
      getPosumerBusIndex(net.prosumpsVec[explicit]) == busIdx || error("addHvdcLink!: $(side)_prosumer $(explicit) is not connected to bus $(bus)")
      return explicit
    end
    gens = [i for (i, p) in enumerate(net.prosumpsVec) if isGenerator(p) && getPosumerBusIndex(p) == busIdx]
    isempty(gens) && error("addHvdcLink!: no injection found at $(side) bus $(bus)")
    length(gens) > 1 && error("addHvdcLink!: $(length(gens)) injections at $(side) bus $(bus); pass $(side)_prosumer to pick one.")
    return gens[1]
  end
  f_idx = resolve(from_bus, from_prosumer, "from")
  t_idx = resolve(to_bus, to_prosumer, "to")
  f_idx == t_idx && error("addHvdcLink!: from and to resolve to the same prosumer ($(f_idx)).")
  link = HvdcLink(something(name, string("HVDC_", from_bus, "_", to_bus)), geNetBusIdx(net = net, busName = from_bus), geNetBusIdx(net = net, busName = to_bus), f_idx, t_idx, status, :api, kind, nothing)
  push!(net.hvdcLinks, link)
  return link
end

# register the controller on the persistent link record: match on the
# terminal prosumers, update in place when the importer already pushed the
# link, otherwise create one (API-built net; kind defaults to :b2b, the
# controller's own abstraction)
function _hvdc_register_link!(net::Net, ctrl::HvdcPairControl)
  for (i, l) in enumerate(net.hvdcLinks)
    if (l.from_prosumer == ctrl.from_prosumer && l.to_prosumer == ctrl.to_prosumer) || (l.from_prosumer == ctrl.to_prosumer && l.to_prosumer == ctrl.from_prosumer)
      net.hvdcLinks[i] = _hvdc_link_with_controller(l, ctrl.name)
      return nothing
    end
  end
  fbus = geNetBusIdx(net = net, busName = ctrl.from_bus)
  tbus = geNetBusIdx(net = net, busName = ctrl.to_bus)
  push!(net.hvdcLinks, HvdcLink(ctrl.name, fbus, tbus, ctrl.from_prosumer, ctrl.to_prosumer, 1, :api, :b2b, ctrl.name))
  return nothing
end

# Resolve the converter prosumer at a bus: explicit index wins, otherwise
# there must be exactly one generator-type prosumer at the bus (the machine
# controller's resolution rule). A voltage-regulating (PV) prosumer is
# acceptable for a P-only side: its reactive power is a solver outcome and
# the pair controller then never touches Q there. It is rejected only when
# the side also requests a fixed Q or a voltage target, because both would
# fight the PV mechanism.
function _hvdc_resolve_prosumer(net::Net, bus::String, explicit::Union{Nothing,Int}, side::String; wants_reactive::Bool, require_slack::Bool = false)::Int
  busIdx = geNetBusIdx(net = net, busName = bus)
  idx = if explicit !== nothing
    (1 <= explicit <= length(net.prosumpsVec)) || error("HvdcPairControl: $(side)_prosumer $(explicit) out of range")
    getPosumerBusIndex(net.prosumpsVec[explicit]) == busIdx || error("HvdcPairControl: $(side)_prosumer $(explicit) is not connected to bus $(bus)")
    explicit
  else
    gens = [i for (i, p) in enumerate(net.prosumpsVec) if isGenerator(p) && getPosumerBusIndex(p) == busIdx]
    isempty(gens) && error("HvdcPairControl: no injection found at $(side) bus $(bus)")
    length(gens) > 1 && error("HvdcPairControl: $(length(gens)) injections at $(side) bus $(bus); pass $(side)_prosumer to pick one.")
    gens[1]
  end
  (wants_reactive && net.prosumpsVec[idx].isRegulated) && error("HvdcPairControl: the converter at $(bus) is voltage-regulating (PV); its Q is a solver outcome, so the $(side) side cannot carry q_mvar or vset_pu.")
  is_slack = getNodeType(net.nodeVec[busIdx]) == Slack
  if require_slack
    # island_feed: the grid-forming converter IS the island reference; its
    # output is the balance outcome the controller mirrors
    is_slack || error("HvdcPairControl: island_feed needs the to-side converter as the island reference (slack), but bus $(bus) is not a reference. Model it as EXTERNALNETWORKINJECTION with referencePri.")
  else
    # at the reference even P is a solver outcome; a slack-bus injection can
    # never be one half of a setpoint pair
    is_slack && error("HvdcPairControl: bus $(bus) is the reference (slack) of its island; its injection balances the island and cannot follow a transfer setpoint. For a grid-forming receiving converter use mode = :island_feed instead.")
  end
  return idx
end

# Stage-0 imports clamp the converter injections to their snapshot value
# (minP == maxP, minQ == maxQ). Under pair control the device is steerable,
# so equality clamps on the controlled prosumer are lifted; real ranges
# (min < max) are kept untouched.
function _hvdc_unclamp!(ps::ProSumer)
  if ps.minP !== nothing && ps.maxP !== nothing && ps.minP == ps.maxP
    ps.minP = nothing
    ps.maxP = nothing
  end
  if ps.minQ !== nothing && ps.maxQ !== nothing && ps.minQ == ps.maxQ
    ps.minQ = nothing
    ps.maxQ = nothing
  end
  return ps
end

"""
    addHvdcPairControl!(net; from_bus, to_bus, p_transfer_mw, ...)

Add a back-to-back HVDC pairing controller coupling the converter
injections at `from_bus` and `to_bus`.

# Arguments
- `net::Net`: the network.
- `from_bus::String`, `to_bus::String`: AC terminal buses of the two
  converters. The pair fixes the sign convention: `p_transfer_mw` leaves
  the from side.
- `mode::Symbol = :setpoint`: `:setpoint` steers a given transfer;
  `:island_feed` models a grid-forming receiving converter that is the
  reference (slack) of its island. There the transfer is derived from the
  island balance each outer iteration and mirrored onto the sending side,
  `P_from = -(P_island + loss)`; `p_transfer_mw` must be omitted and the
  to side carries neither `q_mvar` nor `vset_pu` (the slack holds its own
  voltage). The island draw is read at the reference bus (net injection
  plus local load; a shunt at the PCC is not supported).
- `p_transfer_mw::Float64`: transfer setpoint in MW (signed; negative
  reverses the link). Required in `:setpoint` mode, forbidden in
  `:island_feed`.
- `deadband_p_mw::Float64 = 1e-3`: island_feed only, the mirror counts as
  settled when applied and derived transfer agree within this band.
- `loss_mw::Float64 = 0.0`, `loss_fraction::Float64 = 0.0`: converter loss
  model, `loss = loss_mw + loss_fraction * |p_transfer_mw|` (the MATPOWER
  `LOSS0`/`LOSS1` pair maps directly).
- `p_rating_mw = nothing`: optional transfer rating; `|p_transfer_mw|` is
  clamped to it with honest `at_limit`.
- `from_q_mvar`, `to_q_mvar`: fixed terminal reactive injection in MVAr.
- `from_vset_pu`, `to_vset_pu`: terminal voltage target instead of fixed Q
  (per side exclusive with `q_mvar`); requires the matching
  `from_qmin_mvar`/`from_qmax_mvar` resp. `to_qmin_mvar`/`to_qmax_mvar`
  range (defaults from the prosumer's `minQ`/`maxQ` when present).
- `deadband_vm_pu::Float64 = 1e-3`, `max_outer_iters::Int = 20`,
  `enabled::Bool = true`, `name = nothing`: as in the other controllers.
- `from_prosumer`, `to_prosumer`: explicit prosumer indices when a bus
  carries more than one injection.

The active-power side is a setpoint (applied in the first outer
iteration); the invariant `P_to = p_transfer_mw - loss` holds exactly
after every apply step. Voltage-target terminals iterate a per-side
secant on their reactive injection, mirroring `addMachineVoltageControl!`.
Fails for unknown buses, ambiguous or regulated prosumers, identical
from/to prosumers, a per-side `q_mvar`/`vset_pu` double specification,
a missing reactive range in vset mode, or a second pair controller on one
of the prosumers.
"""
function addHvdcPairControl!(
  net::Net;
  from_bus::String,
  to_bus::String,
  mode::Symbol = :setpoint,
  p_transfer_mw::Union{Nothing,Float64} = nothing,
  deadband_p_mw::Float64 = 1e-3,
  loss_mw::Float64 = 0.0,
  loss_fraction::Float64 = 0.0,
  p_rating_mw::Union{Nothing,Float64} = nothing,
  from_q_mvar::Union{Nothing,Float64} = nothing,
  to_q_mvar::Union{Nothing,Float64} = nothing,
  from_vset_pu::Union{Nothing,Float64} = nothing,
  to_vset_pu::Union{Nothing,Float64} = nothing,
  from_qmin_mvar::Union{Nothing,Float64} = nothing,
  from_qmax_mvar::Union{Nothing,Float64} = nothing,
  to_qmin_mvar::Union{Nothing,Float64} = nothing,
  to_qmax_mvar::Union{Nothing,Float64} = nothing,
  deadband_vm_pu::Float64 = 1e-3,
  max_outer_iters::Int = 20,
  enabled::Bool = true,
  name::Union{Nothing,String} = nothing,
  from_prosumer::Union{Nothing,Int} = nothing,
  to_prosumer::Union{Nothing,Int} = nothing,
)
  from_bus == to_bus && error("HvdcPairControl: from_bus equals to_bus ($(from_bus)); a pair couples two different terminals.")
  mode in (:setpoint, :island_feed) || error("HvdcPairControl: unknown mode $(mode); use :setpoint or :island_feed.")
  if mode == :setpoint
    p_transfer_mw === nothing && error("HvdcPairControl: setpoint mode needs p_transfer_mw.")
  else
    # island_feed derives the transfer from the island balance; a setpoint
    # would contradict the grid-forming semantics
    p_transfer_mw === nothing || error("HvdcPairControl: island_feed derives the transfer from the island balance; omit p_transfer_mw.")
    (to_q_mvar !== nothing || to_vset_pu !== nothing) && error("HvdcPairControl: island_feed's to side is the island reference and holds its own voltage; omit to_q_mvar/to_vset_pu.")
  end
  deadband_p_mw > 0.0 || error("HvdcPairControl: deadband_p_mw must be positive, got $(deadband_p_mw)")
  loss_mw >= 0.0 || error("HvdcPairControl: loss_mw must be nonnegative, got $(loss_mw)")
  (0.0 <= loss_fraction < 1.0) || error("HvdcPairControl: loss_fraction must be in [0, 1), got $(loss_fraction)")
  p_rating_mw === nothing || p_rating_mw > 0.0 || error("HvdcPairControl: p_rating_mw must be positive, got $(p_rating_mw)")
  (from_q_mvar !== nothing && from_vset_pu !== nothing) && error("HvdcPairControl: from side has both q_mvar and vset_pu; pick one.")
  (to_q_mvar !== nothing && to_vset_pu !== nothing) && error("HvdcPairControl: to side has both q_mvar and vset_pu; pick one.")
  from_vset_pu === nothing || from_vset_pu > 0.0 || error("HvdcPairControl: from_vset_pu must be positive")
  to_vset_pu === nothing || to_vset_pu > 0.0 || error("HvdcPairControl: to_vset_pu must be positive")

  f_idx = _hvdc_resolve_prosumer(net, from_bus, from_prosumer, "from"; wants_reactive = from_q_mvar !== nothing || from_vset_pu !== nothing)
  t_idx = _hvdc_resolve_prosumer(net, to_bus, to_prosumer, "to"; wants_reactive = to_q_mvar !== nothing || to_vset_pu !== nothing, require_slack = mode == :island_feed)
  f_idx == t_idx && error("HvdcPairControl: from and to resolve to the same prosumer ($(f_idx)); a pair needs two converters.")
  for c in _hvdc_pair_controllers(net)
    c.enabled || continue
    (c.from_prosumer in (f_idx, t_idx) || c.to_prosumer in (f_idx, t_idx)) && error("HvdcPairControl: a converter of this pair is already controlled by $(c.name).")
  end

  fps = net.prosumpsVec[f_idx]
  tps = net.prosumpsVec[t_idx]
  # reactive ranges: explicit kwargs win, otherwise the prosumer's own
  # limits; only required for a voltage-target side
  fqlo = something(from_qmin_mvar, fps.minQ === nothing ? nothing : fps.minQ, missing)
  fqhi = something(from_qmax_mvar, fps.maxQ === nothing ? nothing : fps.maxQ, missing)
  tqlo = something(to_qmin_mvar, tps.minQ === nothing ? nothing : tps.minQ, missing)
  tqhi = something(to_qmax_mvar, tps.maxQ === nothing ? nothing : tps.maxQ, missing)
  if from_vset_pu !== nothing
    (fqlo === missing || fqhi === missing) && error("HvdcPairControl: from_vset_pu needs a reactive range; pass from_qmin_mvar/from_qmax_mvar.")
    # Stage-0 equality clamps carry no usable range either
    fqlo < fqhi || error("HvdcPairControl: empty from-side reactive range [$(fqlo), $(fqhi)] MVAr; pass from_qmin_mvar/from_qmax_mvar.")
  end
  if to_vset_pu !== nothing
    (tqlo === missing || tqhi === missing) && error("HvdcPairControl: to_vset_pu needs a reactive range; pass to_qmin_mvar/to_qmax_mvar.")
    tqlo < tqhi || error("HvdcPairControl: empty to-side reactive range [$(tqlo), $(tqhi)] MVAr; pass to_qmin_mvar/to_qmax_mvar.")
  end

  # same cross-controller courtesy warning as machine/shunt controllers:
  # two devices steering one voltage fight each other in the outer loop
  for (side_bus, vset) in ((from_bus, from_vset_pu), (to_bus, to_vset_pu))
    vset === nothing && continue
    for tf in net.trafos
      for w in (tf.side1, tf.side2, tf.side3)
        w === nothing && continue
        for c in w.controls
          if c.target_bus == side_bus && c.target_vm_pu !== nothing
            @warn "HvdcPairControl: a transformer tap controller already regulates the voltage of bus $(side_bus); two controllers on one voltage fight each other."
          end
        end
      end
    end
  end

  _hvdc_unclamp!(fps)
  # island_feed leaves the reference injection untouched: it is the island's
  # balance, not an actuator of this controller
  mode == :setpoint && _hvdc_unclamp!(tps)

  ctrl = HvdcPairControl(
    something(name, string("B2B_", from_bus, "_", to_bus)),
    from_bus,
    to_bus,
    f_idx,
    t_idx,
    # island_feed starts at 0 and mirrors the island balance from the first
    # outer iteration on
    mode == :setpoint ? p_transfer_mw : 0.0,
    loss_mw,
    loss_fraction,
    p_rating_mw,
    from_q_mvar,
    to_q_mvar,
    from_vset_pu,
    to_vset_pu,
    fqlo === missing ? -Inf : Float64(fqlo),
    fqhi === missing ? Inf : Float64(fqhi),
    tqlo === missing ? -Inf : Float64(tqlo),
    tqhi === missing ? Inf : Float64(tqhi),
    deadband_vm_pu,
    max_outer_iters,
    enabled,
    :idle,
    false,
    false,
    false,
    false,
    fps.qVal === nothing ? 0.0 : fps.qVal,
    tps.qVal === nothing ? 0.0 : tps.qVal,
    nothing,
    nothing,
    nothing,
    nothing,
    nothing,
    nothing,
    0,
    mode,
    deadband_p_mw,
  )
  push!(net.machineControls, ctrl)
  _hvdc_register_link!(net, ctrl)
  return ctrl
end

# defensive prosumer lookup at actuation time (index still on the named bus)
function _hvdc_pair_prosumer(net::Net, ctrl::HvdcPairControl, idx::Int, bus::String)::ProSumer
  (1 <= idx <= length(net.prosumpsVec)) || error("HvdcPairControl $(ctrl.name): prosumer index $(idx) out of range")
  ps = net.prosumpsVec[idx]
  getPosumerBusIndex(ps) == geNetBusIdx(net = net, busName = bus) || error("HvdcPairControl $(ctrl.name): prosumer $(idx) is not connected to bus $(bus)")
  return ps
end

# transfer after the optional rating clamp, and the resulting pair of
# injection setpoints (from side, to side) honoring the invariant
function _hvdc_pair_setpoints(ctrl::HvdcPairControl)
  transfer = ctrl.p_transfer_mw
  clamped = false
  if ctrl.p_rating_mw !== nothing && abs(transfer) > ctrl.p_rating_mw
    transfer = sign(transfer) * ctrl.p_rating_mw
    clamped = true
  end
  loss = ctrl.loss_mw + ctrl.loss_fraction * abs(transfer)
  return (from_p = -transfer, to_p = transfer - loss, transfer = transfer, loss = loss, clamped = clamped)
end

# island_feed: what the grid-forming converter currently delivers into its
# island (MW), i.e. the terminal output of the reference injection at the
# PCC. Single-island runs have the slack's NET injection in _pƩGen (add the
# local load back). Island runs write only voltages back to the parent net,
# so the fallback derives the output from the fresh branch flows instead
# (run_control! calls calcNetLosses! after every inner solve): local load
# plus the active power flowing from the PCC into its in-service branches.
# Returns nothing before the first solve. A shunt at the PCC is not
# supported (documented in addHvdcPairControl!).
function _hvdc_island_delivered(net::Net, ctrl::HvdcPairControl)::Union{Nothing,Float64}
  busIdx = geNetBusIdx(net = net, busName = ctrl.to_bus)
  node = net.nodeVec[busIdx]
  local_load = node._pƩLoad === nothing ? 0.0 : node._pƩLoad
  node._pƩGen === nothing || return node._pƩGen + local_load
  total = local_load
  any_flow = false
  for br in net.branchVec
    br.status == 1 || continue
    if br.fromBus == busIdx && br.fBranchFlow !== nothing
      total += br.fBranchFlow.pFlow
      any_flow = true
    elseif br.toBus == busIdx && br.tBranchFlow !== nothing
      total += br.tBranchFlow.pFlow
      any_flow = true
    end
  end
  return any_flow ? total : nothing
end

# invert the loss model: the transfer whose to-side delivery equals the
# island draw, delivered = transfer - loss_mw - loss_fraction * |transfer|
function _hvdc_island_feed_transfer(ctrl::HvdcPairControl, delivered::Float64)::Float64
  num = delivered + ctrl.loss_mw
  return num >= 0.0 ? num / (1.0 - ctrl.loss_fraction) : num / (1.0 + ctrl.loss_fraction)
end

# --- AbstractOuterController protocol ---------------------------------------

control_name(ctrl::HvdcPairControl) = string(ctrl.name, " B2B")
control_enabled(ctrl::HvdcPairControl) = ctrl.enabled
control_initialize!(::HvdcPairControl, ::Net, context) = NoControlState()
control_status(ctrl::HvdcPairControl, ::AbstractControlState)::Symbol = ctrl.status
control_is_converged(ctrl::HvdcPairControl, ::AbstractControlState)::Bool = ctrl.converged
control_is_blocked(ctrl::HvdcPairControl, ::AbstractControlState)::Bool = ctrl.at_limit

# a side is voltage-converged when it has no voltage target (fixed Q or
# untouched) or its measured terminal voltage sits inside the deadband
_hvdc_side_converged(vset::Union{Nothing,Float64}, vm::Union{Nothing,Float64}, deadband::Float64)::Bool = vset === nothing || (vm !== nothing && abs(vm - vset) <= deadband)

function control_evaluate!(ctrl::HvdcPairControl, net::Net, ::AbstractControlState, context)
  ctrl.from_vm_now = get_bus_vm_pu(net, ctrl.from_bus)
  ctrl.to_vm_now = get_bus_vm_pu(net, ctrl.to_bus)
  # island_feed requires the grid-forming to bus to BE the island reference
  # at runtime, not only at registration: a user who demoted it (e.g. after
  # synchronously tying the islands) invalidates the mirror semantics. The
  # controller reports that honestly and stops mirroring; it never switches
  # modes silently.
  if ctrl.mode == :island_feed && getNodeType(net.nodeVec[geNetBusIdx(net = net, busName = ctrl.to_bus)]) != Slack
    ctrl.converged = false
    ctrl.status = :invalid_topology
    ctrl.outer_iters = context.outer_iteration
    return nothing
  end
  # island_feed: the mirror is settled when the applied transfer matches the
  # one derived from the fresh island balance
  mirror_ok = true
  if ctrl.mode == :island_feed
    delivered = _hvdc_island_delivered(net, ctrl)
    mirror_ok = delivered !== nothing && abs(_hvdc_island_feed_transfer(ctrl, delivered) - ctrl.p_transfer_mw) <= ctrl.deadband_p_mw
  end
  ctrl.converged = ctrl.p_applied && !ctrl.transfer_clamped && mirror_ok && _hvdc_side_converged(ctrl.from_vset_pu, ctrl.from_vm_now, ctrl.deadband_vm_pu) && _hvdc_side_converged(ctrl.to_vset_pu, ctrl.to_vm_now, ctrl.deadband_vm_pu)
  ctrl.status = ctrl.converged ? :converged : (ctrl.at_limit ? :at_limit : :active)
  ctrl.outer_iters = context.outer_iteration
  return nothing
end

# per-side secant step toward the voltage target (machine-controller logic)
function _hvdc_side_q_step(vset::Float64, vm::Union{Nothing,Float64}, q0::Float64, prev_q::Union{Nothing,Float64}, prev_vm::Union{Nothing,Float64}, qlo::Float64, qhi::Float64)::Float64
  vm === nothing && return q0
  err = vset - vm
  dq = if prev_q !== nothing && prev_vm !== nothing && abs(q0 - prev_q) > 1e-9
    slope = (vm - prev_vm) / (q0 - prev_q)
    slope > _MACHINE_CTRL_MIN_SLOPE ? err / slope : _MACHINE_CTRL_BOOTSTRAP_FRACTION * (err > 0.0 ? (qhi - q0) : (qlo - q0))
  else
    _MACHINE_CTRL_BOOTSTRAP_FRACTION * (err > 0.0 ? (qhi - q0) : (qlo - q0))
  end
  return clamp(q0 + dq, qlo, qhi)
end

function control_propose_update!(ctrl::HvdcPairControl, net::Net, ::AbstractControlState, context)
  if ctrl.mode == :island_feed && ctrl.status != :invalid_topology
    # refresh the mirrored transfer from the island balance of the last solve
    delivered = _hvdc_island_delivered(net, ctrl)
    delivered === nothing || (ctrl.p_transfer_mw = _hvdc_island_feed_transfer(ctrl, delivered))
  end
  sp = _hvdc_pair_setpoints(ctrl)
  # fixed-Q sides are one-shot setpoints like P; vset sides iterate
  new_from_q = if ctrl.from_vset_pu !== nothing
    ctrl.converged ? ctrl.from_q_now : _hvdc_side_q_step(ctrl.from_vset_pu, ctrl.from_vm_now, ctrl.from_q_now, ctrl.from_prev_q, ctrl.from_prev_vm, ctrl.from_qmin_mvar, ctrl.from_qmax_mvar)
  elseif ctrl.from_q_mvar !== nothing
    ctrl.from_q_mvar
  else
    ctrl.from_q_now
  end
  new_to_q = if ctrl.to_vset_pu !== nothing
    ctrl.converged ? ctrl.to_q_now : _hvdc_side_q_step(ctrl.to_vset_pu, ctrl.to_vm_now, ctrl.to_q_now, ctrl.to_prev_q, ctrl.to_prev_vm, ctrl.to_qmin_mvar, ctrl.to_qmax_mvar)
  elseif ctrl.to_q_mvar !== nothing
    ctrl.to_q_mvar
  else
    ctrl.to_q_now
  end
  return (from_p = sp.from_p, to_p = sp.to_p, transfer_clamped = sp.clamped, from_q = new_from_q, to_q = new_to_q, old_from_q = ctrl.from_q_now, old_to_q = ctrl.to_q_now)
end

function control_apply_update!(ctrl::HvdcPairControl, net::Net, ::AbstractControlState, update::NamedTuple, context)::Bool
  moved = false
  fps = _hvdc_pair_prosumer(net, ctrl, ctrl.from_prosumer, ctrl.from_bus)
  tps = _hvdc_pair_prosumer(net, ctrl, ctrl.to_prosumer, ctrl.to_bus)
  fnode = net.nodeVec[geNetBusIdx(net = net, busName = ctrl.from_bus)]
  tnode = net.nodeVec[geNetBusIdx(net = net, busName = ctrl.to_bus)]
  # active power: applied once (and again only if the setpoints drifted);
  # bus sums and prosumer values move by the same delta, the machine
  # controller's coherence rule
  f_p_now = fps.pVal === nothing ? 0.0 : fps.pVal
  t_p_now = tps.pVal === nothing ? 0.0 : tps.pVal
  if ctrl.mode == :island_feed
    # only the sending side is an actuator; the to side is the island's
    # reference and its injection stays the solver's balance outcome.
    # invalid_topology stops all mirroring until the user fixes the model.
    if ctrl.status != :invalid_topology && (!ctrl.p_applied || abs(f_p_now - update.from_p) > 1e-9)
      addGenPower!(node = fnode, p = update.from_p - f_p_now, q = nothing)
      fps.pVal = update.from_p
      ctrl.p_applied = true
      moved = true
    end
  elseif !ctrl.p_applied || abs(f_p_now - update.from_p) > 1e-9 || abs(t_p_now - update.to_p) > 1e-9
    addGenPower!(node = fnode, p = update.from_p - f_p_now, q = nothing)
    addGenPower!(node = tnode, p = update.to_p - t_p_now, q = nothing)
    fps.pVal = update.from_p
    tps.pVal = update.to_p
    ctrl.p_applied = true
    moved = true
  end
  if abs(update.from_q - update.old_from_q) > 1e-9
    addGenPower!(node = fnode, p = nothing, q = update.from_q - (fps.qVal === nothing ? 0.0 : fps.qVal))
    fps.qVal = update.from_q
    ctrl.from_prev_q = update.old_from_q
    ctrl.from_prev_vm = ctrl.from_vm_now
    ctrl.from_q_now = update.from_q
    moved = true
  end
  if abs(update.to_q - update.old_to_q) > 1e-9
    addGenPower!(node = tnode, p = nothing, q = update.to_q - (tps.qVal === nothing ? 0.0 : tps.qVal))
    tps.qVal = update.to_q
    ctrl.to_prev_q = update.old_to_q
    ctrl.to_prev_vm = ctrl.to_vm_now
    ctrl.to_q_now = update.to_q
    moved = true
  end
  # honest at_limit: transfer clamped by the rating, or a voltage-target
  # side clamped at its reactive bound without having converged
  ctrl.transfer_clamped = update.transfer_clamped
  ctrl.transfer_clamped && (ctrl.converged = false)
  from_q_clamped = ctrl.from_vset_pu !== nothing && (isapprox(ctrl.from_q_now, ctrl.from_qmin_mvar; atol = 1e-9) || isapprox(ctrl.from_q_now, ctrl.from_qmax_mvar; atol = 1e-9))
  to_q_clamped = ctrl.to_vset_pu !== nothing && (isapprox(ctrl.to_q_now, ctrl.to_qmin_mvar; atol = 1e-9) || isapprox(ctrl.to_q_now, ctrl.to_qmax_mvar; atol = 1e-9))
  ctrl.at_limit = ctrl.transfer_clamped || (!ctrl.converged && (from_q_clamped || to_q_clamped))
  ctrl.at_limit && (ctrl.status = :at_limit)
  return moved
end

function control_element_descriptor(ctrl::HvdcPairControl, net::Net)::Union{Nothing,NamedTuple}
  sp = _hvdc_pair_setpoints(ctrl)
  return (
    name = control_name(ctrl),
    element = string("hvdc@", ctrl.from_bus, "-", ctrl.to_bus),
    device = ctrl.mode == :island_feed ? "back-to-back HVDC pair (grid-forming)" : "back-to-back HVDC pair",
    actuator = :hvdc_p_transfer_mw,
    actuator_min = ctrl.p_rating_mw === nothing ? -Inf : -ctrl.p_rating_mw,
    actuator_max = ctrl.p_rating_mw === nothing ? Inf : ctrl.p_rating_mw,
    quantity = :hvdc_transfer,
    target = string(ctrl.from_bus, "->", ctrl.to_bus),
    target_value = ctrl.p_transfer_mw,
    discrete = false,
    enabled = ctrl.enabled,
    status = ctrl.status,
    converged = ctrl.converged,
    at_limit = ctrl.at_limit,
  )
end

function control_report_rows(ctrl::HvdcPairControl, net::Net, ::AbstractControlState, context)
  sp = _hvdc_pair_setpoints(ctrl)
  return [(
    controller_name = control_name(ctrl),
    controller_type = "HvdcPairControl",
    mode = ctrl.mode,
    link = string(ctrl.from_bus, "->", ctrl.to_bus),
    status = ctrl.status,
    converged = ctrl.converged,
    at_limit = ctrl.at_limit,
    p_transfer_mw = sp.transfer,
    loss_mw = sp.loss,
    from_p_mw = sp.from_p,
    to_p_mw = sp.to_p,
    from_q_mvar = ctrl.from_q_now,
    to_q_mvar = ctrl.to_q_now,
    from_vm_pu = ctrl.from_vm_now === nothing ? missing : ctrl.from_vm_now,
    to_vm_pu = ctrl.to_vm_now === nothing ? missing : ctrl.to_vm_now,
  )]
end

function control_trace_rows(ctrl::HvdcPairControl, net::Net, ::AbstractControlState, context)
  sp = _hvdc_pair_setpoints(ctrl)
  return [(
    outer_iteration = context.outer_iteration,
    controller_name = control_name(ctrl),
    controller_type = "HvdcPairControl",
    link = string(ctrl.from_bus, "->", ctrl.to_bus),
    status = ctrl.status,
    converged = ctrl.converged,
    at_limit = ctrl.at_limit,
    p_transfer_mw = sp.transfer,
    from_q_mvar = ctrl.from_q_now,
    to_q_mvar = ctrl.to_q_now,
  )]
end

"""
    printHvdcPairControllerSummary(io::IO, net::Net)
    printHvdcPairControllerSummary(net::Net)

Engineering-style summary of the registered HVDC pair controllers, one
block per pair: link and direction, transfer and loss, per-terminal P/Q
and voltage state, and the honest limit flags. Prints nothing when no pair
controller is registered.
"""
function printHvdcPairControllerSummary(io::IO, net::Net)
  ctrls = _hvdc_pair_controllers(net)
  isempty(ctrls) && return
  println(io, "\nHVDC Pair Control Summary (back-to-back)")
  println(io, "----------------------------------------")
  for c in ctrls
    sp = _hvdc_pair_setpoints(c)
    println(io, control_name(c), " (link ", c.from_bus, " -> ", c.to_bus, ")")
    c.mode == :island_feed && println(io, "  mode               : island_feed (grid-forming to side)")
    println(io, "  transfer           : ", @sprintf("%.3f MW", sp.transfer), c.mode == :island_feed ? " (mirrored from island draw)" : "", c.p_rating_mw === nothing ? "" : @sprintf(" (rating %.1f MW)", c.p_rating_mw))
    println(io, "  loss               : ", @sprintf("%.3f MW", sp.loss))
    println(io, "  from injection     : ", @sprintf("%.3f MW / %.3f MVAr", sp.from_p, c.from_q_now))
    println(io, "  to injection       : ", @sprintf("%.3f MW / %.3f MVAr", sp.to_p, c.to_q_now), c.mode == :island_feed ? " (island balance outcome)" : "")
    c.from_vset_pu === nothing || println(io, "  from voltage       : ", c.from_vm_now === nothing ? "-" : @sprintf("%.4f pu", c.from_vm_now), @sprintf(" (target %.4f pu)", c.from_vset_pu))
    c.to_vset_pu === nothing || println(io, "  to voltage         : ", c.to_vm_now === nothing ? "-" : @sprintf("%.4f pu", c.to_vm_now), @sprintf(" (target %.4f pu)", c.to_vset_pu))
    println(io, "  converged          : ", c.converged)
    println(io, "  at_limit           : ", c.at_limit)
    println(io, "  status             : ", c.status)
  end
end
printHvdcPairControllerSummary(net::Net) = printHvdcPairControllerSummary(stdout, net)
