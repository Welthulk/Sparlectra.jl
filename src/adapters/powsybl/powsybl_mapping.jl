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

# file: src/adapters/powsybl/powsybl_mapping.jl
# purpose: build a Net from PowSyBl tables. Bus model: one bus per
#          bus-breaker bus, retained switches as links. Every element
#          mapper counts what it read, built and skipped, and every skip
#          carries a reason in the report. Units arrive as pypowsybl
#          delivers them (ohm, siemens, kV, degrees, MW, MVar, A) and are
#          converted here, once, to Sparlectra's per-unit convention:
#          branch impedance on the to-side nominal voltage, the complex
#          ratio at the from side (the MATPOWER model calcAdmittance uses).

# Conventions settled against the OpenLoadFlow references of the shipped
# example bundles (the numbers are the worst bus deviations measured):
# - two-winding ratio: Sparlectra's ratio is V_from_pu / V_to_pu at no load
#   (the MATPOWER tap), IIDM gives V2 = rho * V1 in kV with rho already
#   carrying rated_u2 / rated_u1, so ratio = vn_to / (rho * vn_from); the
#   two alternatives considered first (rho on the bus nominal voltages, rho
#   as is) were 0.78 pu off on four_substations and could not converge
#   micro_grid_be, and only this form reproduces the case14 taps (0.978).
# - phase shift: IIDM alpha turns V2 against V1, MATPOWER's shift turns the
#   from side against the internal to side: shift = -alpha (+alpha was 5.4
#   degrees off at the four_substations PST).
# - magnetizing admittance: OpenLoadFlow places the whole g, b on pi side 1,
#   behind the ideal transformer on the side-2 base (1.5e-7 pu on
#   four_substations; side 2 gave 3.4e-4, a split 1.7e-4).
# - dangling line: g, b at the network terminal, nothing at the boundary
#   (3e-11 pu on the eurostag tie lines; split 9e-4, boundary side 1.8e-3).
# - lines across nominal voltages: the physical conductor, which is the
#   ratio branch on the to-side base of pi_branch_pu_between_levels; the
#   geometric base vn1 * vn2 / S of PowSyBl's importers is not a load-flow
#   model (0.93 pu off on ieee14).

# A limit or reactive bound at or beyond this magnitude is "unbounded"
# (pypowsybl writes 1.797693e308 for a missing bound).
const _POWSYBL_UNBOUNDED = 1.0e300

_powsybl_finite(x::Float64) = isfinite(x) && abs(x) < _POWSYBL_UNBOUNDED
_powsybl_optional(x::Float64) = _powsybl_finite(x) ? x : nothing

struct _PowsyblBuildContext
  net::Net
  tables::PowsyblTables
  opts::PowsyblAdapterOptions
  report::PowsyblImportReport
  idx::Dict{String,Int}                 # bus-breaker id to Sparlectra bus index
  vn::Dict{String,Float64}              # bus-breaker id to nominal voltage
  component::Dict{String,Int}           # bus-breaker id to synchronous component (-1 when none)
  vl_nominal::Dict{String,Float64}      # voltage level id to nominal voltage
  bus_view_vn::Dict{String,Float64}     # bus-view id to nominal voltage
  limits::Dict{String,Vector{Int}}      # element id to rows of operational_limits
  machine_plans::Vector{NamedTuple}     # remote voltage controls to attach once every prosumer exists
end

_powsybl_vn(ctx::_PowsyblBuildContext, bus::AbstractString) = ctx.vn[bus]

# The bus of an element, or nothing with the skip recorded.
function _powsybl_bus(ctx::_PowsyblBuildContext, kind::AbstractString, id::AbstractString, bus::AbstractString)
  haskey(ctx.idx, bus) && return bus
  _powsybl_skip!(ctx.report, kind, id, isempty(bus) ? "no bus" : "bus $(bus) not built")
  return nothing
end

# Branch status from the two terminal flags: both closed is in service,
# both open is out of service, one open is a one-sided open branch.
function _powsybl_terminal_status(c1::Bool, c2::Bool)
  status = (c1 || c2) ? 1 : 0
  return (status = status, from_status = c1 ? 1 : 0, to_status = c2 ? 1 : 0)
end

"""
    build_net_from_powsybl(tables::PowsyblTables, opts::PowsyblAdapterOptions; name = "") -> (Net, PowsyblImportReport)

Build a Sparlectra network from PowSyBl tables: buses and links, lines,
two- and three-winding transformers, generators (with the slack per
synchronous component), loads, shunts, static var compensators,
dangling and tie lines, HVDC links. The report says what was built,
what was skipped and why, and which generator carries the reference of
each component. `name` defaults to the manifest's case name.
"""
function build_net_from_powsybl(tables::PowsyblTables, opts::PowsyblAdapterOptions; name::AbstractString = "")::Tuple{Net,PowsyblImportReport}
  net_name = isempty(name) ? String(get(tables.manifest, "case", "powsybl")) : String(name)
  net = Net(name = net_name, baseMVA = opts.base_mva)
  report = PowsyblImportReport()
  vl = tables.voltage_levels
  vl_nominal = Dict{String,Float64}(vl.id[i] => vl.nominal_v[i] for i in eachindex(vl.id))
  bv = tables.buses
  bus_view_vn = Dict{String,Float64}(bv.id[i] => get(vl_nominal, bv.voltage_level_id[i], NaN) for i in eachindex(bv.id))
  ol = tables.operational_limits
  limits = Dict{String,Vector{Int}}()
  for i in eachindex(ol.element_id)
    push!(get!(limits, ol.element_id[i], Int[]), i)
  end
  ctx = _PowsyblBuildContext(net, tables, opts, report, Dict{String,Int}(), Dict{String,Float64}(), Dict{String,Int}(), vl_nominal, bus_view_vn, limits, NamedTuple[])
  _map_powsybl_buses!(ctx)
  _map_powsybl_links!(ctx)
  _map_powsybl_lines!(ctx)
  _map_powsybl_two_winding_transformers!(ctx)
  _map_powsybl_three_winding_transformers!(ctx)
  _map_powsybl_generators!(ctx)
  _map_powsybl_loads!(ctx)
  _map_powsybl_shunts!(ctx)
  _map_powsybl_svcs!(ctx)
  _map_powsybl_dangling_lines!(ctx)
  _map_powsybl_hvdc!(ctx)
  _attach_powsybl_machine_controls!(ctx)
  return (net, report)
end

# remote_regulation = :remote: the unit is a PQ injection with an outer-loop
# MachineVoltageControl that drives the regulated bus onto the target, the
# way cgmes_import.machine_control does it. Attached last, once every
# prosumer exists and the bus types are settled; a control the framework
# refuses (target bus already voltage-held, isolated) stays a notice and
# the unit keeps its target_q.
function _attach_powsybl_machine_controls!(ctx::_PowsyblBuildContext)
  refreshBusTypesFromProsumers!(ctx.net)
  for plan in ctx.machine_plans
    tidx = get(ctx.net.busDict, plan.target_bus, nothing)
    reason = tidx === nothing ? "target bus $(plan.target_bus) not built" : ctx.net.nodeVec[tidx]._nodeType != PQ ? "target bus $(plan.target_bus) is already voltage-held" : nothing
    if reason !== nothing
      _powsybl_notice!(ctx.report, "generator $(plan.name): remote voltage control not attached ($(reason)); stays PQ with target_q")
      continue
    end
    # OpenLoadFlow holds a remotely regulated bus exactly; the controller's
    # default deadband of 1e-3 pu would leave that much slack
    addMachineVoltageControl!(ctx.net; bus = plan.bus, target_bus = plan.target_bus, target_vm_pu = plan.target_vm_pu, qmin_mvar = plan.qmin, qmax_mvar = plan.qmax, deadband_vm_pu = 1.0e-6, prosumer_index = plan.prosumer_index, name = plan.name)
    _powsybl_notice!(ctx.report, "generator $(plan.name): remote voltage control of bus $(plan.target_bus) at $(round(plan.target_vm_pu; digits = 5)) pu attached (outer loop)")
  end
  return nothing
end

# --- 2.1 buses and links ------------------------------------------------------

function _map_powsybl_buses!(ctx::_PowsyblBuildContext)
  bb = ctx.tables.bus_breaker_view_buses
  bv = ctx.tables.buses
  sync = Dict{String,Int}(bv.id[i] => bv.synchronous_component[i] for i in eachindex(bv.id))
  vl_sub = Dict{String,String}(ctx.tables.voltage_levels.id[i] => ctx.tables.voltage_levels.substation_id[i] for i in eachindex(ctx.tables.voltage_levels.id))
  _powsybl_count!(ctx.report, "bus"; read = length(bb.id))
  unconnected = 0
  for i in eachindex(bb.id)
    id = bb.id[i]
    vlid = bb.voltage_level_id[i]
    vn = get(ctx.vl_nominal, vlid, NaN)
    if !(vn > 0.0)
      _powsybl_skip!(ctx.report, "bus", id, haskey(ctx.vl_nominal, vlid) ? "voltage level $(vlid) has nominal_v $(vn)" : "voltage level $(vlid) not in voltage_levels")
      continue
    end
    # the IIDM state seeds the start voltages the way SV does on the CGMES path
    vm = (isfinite(bb.v_mag[i]) && bb.v_mag[i] > 0.0) ? bb.v_mag[i] / vn : 1.0
    va = isfinite(bb.v_angle[i]) ? bb.v_angle[i] : 0.0
    addBus!(net = ctx.net, busName = id, vn_kV = vn, vm_pu = vm, va_deg = va)
    busIdx = ctx.net.busDict[id]
    # the bus-view id travels in the component id, the external-id channel
    # the SCF export writes, so results join back to the OLF bus table
    ctx.net.nodeVec[busIdx].comp.cID = bb.bus_id[i]
    ctx.idx[id] = busIdx
    ctx.vn[id] = vn
    ctx.report.bus_view[id] = bb.bus_id[i]
    ctx.report.bus_substation[id] = get(vl_sub, vlid, "")
    comp = get(sync, bb.bus_id[i], bb.synchronous_component[i])
    ctx.component[id] = comp < 0 ? -1 : comp
    comp < 0 && (unconnected += 1)
    _powsybl_count!(ctx.report, "bus"; built = 1)
  end
  unconnected > 0 && _powsybl_notice!(ctx.report, "$(unconnected) bus(es) belong to no synchronous component (built, unconnected)")
  return nothing
end

function _map_powsybl_links!(ctx::_PowsyblBuildContext)
  sw = ctx.tables.switches
  retained = 0
  contracted = 0
  for i in eachindex(sw.id)
    if !sw.retained[i]
      contracted += 1
      continue
    end
    retained += 1
    id = sw.id[i]
    b1 = sw.bus_breaker_bus1_id[i]
    b2 = sw.bus_breaker_bus2_id[i]
    if b1 == b2
      _powsybl_skip!(ctx.report, "switch", id, "both ends on bus-breaker bus $(b1)")
      continue
    end
    (_powsybl_bus(ctx, "switch", id, b1) === nothing || _powsybl_bus(ctx, "switch", id, b2) === nothing) && continue
    ctx.vn[b1] == ctx.vn[b2] || throw(ArgumentError("powsybl import: retained switch $(id) joins bus-breaker buses $(b1) ($(ctx.vn[b1]) kV) and $(b2) ($(ctx.vn[b2]) kV) of different nominal voltage"))
    addLink!(net = ctx.net, fromBus = b1, toBus = b2, status = sw.open[i] ? 0 : 1)
    _powsybl_count!(ctx.report, "switch"; built = 1)
  end
  _powsybl_count!(ctx.report, "switch"; read = length(sw.id))
  contracted > 0 && _powsybl_notice!(ctx.report, "$(contracted) non-retained switch(es) contracted by PowSyBl")
  return nothing
end

# --- limits -------------------------------------------------------------------

# Permanent current limits become the branch rating in MVA at the nominal
# voltage of the limit's side (the smaller of the two sides wins);
# temporary limits are kept as (duration_s, value_A) pairs in the report.
function _powsybl_branch_rating(ctx::_PowsyblBuildContext, id::AbstractString, vn1::Float64, vn2::Float64)
  rows = get(ctx.limits, id, nothing)
  rows === nothing && return nothing
  ol = ctx.tables.operational_limits
  rating = nothing
  other_types = Set{String}()
  temps = Tuple{Int,Float64}[]
  for i in rows
    if ol.type[i] != "CURRENT"
      push!(other_types, ol.type[i])
      continue
    end
    value = ol.value[i]
    _powsybl_finite(value) || continue
    if ol.acceptable_duration[i] < 0
      vn = ol.side[i] == "TWO" ? vn2 : vn1
      s = sqrt(3.0) * vn * value / 1000.0
      rating = rating === nothing ? s : min(rating, s)
    else
      push!(temps, (ol.acceptable_duration[i], value))
    end
  end
  isempty(temps) || (ctx.report.temporary_limits[String(id)] = temps)
  for t in sort!(collect(other_types))
    _powsybl_notice!(ctx.report, "limit type $(t) of $(id) skipped (only CURRENT limits are mapped)")
  end
  return rating
end

# --- 2.2 lines ------------------------------------------------------------------

function _map_powsybl_lines!(ctx::_PowsyblBuildContext)
  ln = ctx.tables.lines
  _powsybl_count!(ctx.report, "line"; read = length(ln.id))
  for i in eachindex(ln.id)
    id = ln.id[i]
    b1 = _powsybl_bus(ctx, "line", id, ln.bus_breaker_bus1_id[i])
    b1 === nothing && continue
    b2 = _powsybl_bus(ctx, "line", id, ln.bus_breaker_bus2_id[i])
    b2 === nothing && continue
    if b1 == b2
      _powsybl_skip!(ctx.report, "line", id, "both ends on bus $(b1)")
      continue
    end
    vn1 = ctx.vn[b1]
    vn2 = ctx.vn[b2]
    st = _powsybl_terminal_status(ln.connected1[i], ln.connected2[i])
    rating = _powsybl_branch_rating(ctx, id, vn1, vn2)
    # Sparlectra's pi model splits the branch shunt in equal halves; the
    # symmetric part travels on the branch, the excess of the larger side
    # becomes a bus shunt on that side, so unequal g1/g2 and b1/b2 are exact
    g_sym = min(ln.g1[i], ln.g2[i])
    b_sym = min(ln.b1[i], ln.b2[i])
    _powsybl_add_pi_line!(ctx, b1, b2, ln.r[i], ln.x[i], 2.0 * g_sym, 2.0 * b_sym, st, rating)
    _powsybl_bus_shunt!(ctx, b1, ln.g1[i] - g_sym, ln.b1[i] - b_sym)
    _powsybl_bus_shunt!(ctx, b2, ln.g2[i] - g_sym, ln.b2[i] - b_sym)
    _powsybl_count!(ctx.report, "line"; built = 1)
  end
  return nothing
end

# A shunt in siemens at a bus, as MW and MVar at the nominal voltage.
function _powsybl_bus_shunt!(ctx::_PowsyblBuildContext, bus::AbstractString, g::Float64, b::Float64)
  (g == 0.0 && b == 0.0) && return nothing
  (isfinite(g) && isfinite(b)) || return nothing
  vn = ctx.vn[bus]
  addShuntMatpower!(net = ctx.net, busName = String(bus), Gs = g * vn^2, Bs = b * vn^2)
  return nothing
end

# A pi line in ohm and siemens between two buses. Across different nominal
# voltages PowSyBl keeps a plain conductor (no ideal transformer), and
# OpenLoadFlow's default line model (linePerUnitMode = IMPEDANCE) is exactly
# that physical line; on Sparlectra's convention the physical line is the
# ratio branch of the CGMES path: impedance on the to-side base with the
# ratio vn_to / vn_from at the from side (I1 = (V1 - V2) / Z in kV and A
# reproduces term by term).
function _powsybl_add_pi_line!(ctx::_PowsyblBuildContext, b1::AbstractString, b2::AbstractString, r::Float64, x::Float64, g::Float64, b::Float64, st, rating)
  pu = pi_branch_pu_between_levels(r = r, x = x, g = g, b = b, vn_from_kV = ctx.vn[b1], vn_to_kV = ctx.vn[b2], baseMVA = ctx.net.baseMVA)
  if pu.ratio === nothing
    addPIModelACLine!(net = ctx.net, fromBus = String(b1), toBus = String(b2), r_pu = pu.r_pu, x_pu = pu.x_pu, b_pu = pu.b_pu, g_pu = pu.g_pu, status = st.status, ratedS = rating, from_status = st.from_status, to_status = st.to_status)
  else
    _addPIModelTrafo_by_idx!(net = ctx.net, from = ctx.idx[b1], to = ctx.idx[b2], r_pu = pu.r_pu, x_pu = pu.x_pu, b_pu = pu.b_pu, g_pu = pu.g_pu, status = st.status, ratedS = rating, ratio = pu.ratio, shift_deg = 0.0, from_status = st.from_status, to_status = st.to_status)
  end
  return nothing
end

# --- 2.3 two-winding transformers ------------------------------------------------

# Sparlectra's ratio (the MATPOWER tap, V_from_pu / V_to_pu at no load)
# from the IIDM ratio rho at side 1 (V2 = rho * V1 in kV).
_powsybl_2wt_ratio(rho::Float64, vn_from::Float64, vn_to::Float64)::Float64 = vn_to / (rho * vn_from)

function _map_powsybl_two_winding_transformers!(ctx::_PowsyblBuildContext)
  tw = ctx.tables.two_windings_transformers
  _powsybl_count!(ctx.report, "2wt"; read = length(tw.id))
  for i in eachindex(tw.id)
    id = tw.id[i]
    b1 = _powsybl_bus(ctx, "2wt", id, tw.bus_breaker_bus1_id[i])
    b1 === nothing && continue
    b2 = _powsybl_bus(ctx, "2wt", id, tw.bus_breaker_bus2_id[i])
    b2 === nothing && continue
    if b1 == b2
      _powsybl_skip!(ctx.report, "2wt", id, "both ends on bus $(b1)")
      continue
    end
    vn1 = ctx.vn[b1]
    vn2 = ctx.vn[b2]
    # the series impedance is given at side 2 (the to side): pu on the
    # to-side base, the branch carries no shunt of its own
    g_sh = isfinite(tw.g_at_current_tap[i]) ? tw.g_at_current_tap[i] : 0.0
    b_sh = isfinite(tw.b_at_current_tap[i]) ? tw.b_at_current_tap[i] : 0.0
    r_pu, x_pu, b_pu, g_pu = toPU_RXBG(r = tw.r_at_current_tap[i], x = tw.x_at_current_tap[i], g = 0.0, b = 0.0, v_kv = vn2, baseMVA = ctx.net.baseMVA)
    ratio = _powsybl_2wt_ratio(tw.rho[i], vn1, vn2)
    # OpenLoadFlow puts the whole magnetizing admittance on pi side 1, which
    # sits behind the ideal transformer on the side-2 base: seen from bus 1
    # that is y / ratio^2 on the side-2 base, (g, b) * vn2^2 / ratio^2 in MW/MVar
    k = vn2^2 / ratio^2 / vn1^2
    _powsybl_bus_shunt!(ctx, b1, g_sh * k, b_sh * k)
    shift = -(isfinite(tw.alpha[i]) ? tw.alpha[i] : 0.0)
    st = _powsybl_terminal_status(tw.connected1[i], tw.connected2[i])
    rating = _powsybl_branch_rating(ctx, id, vn1, vn2)
    rated_s = _powsybl_optional(tw.rated_s[i])
    _addPIModelTrafo_by_idx!(net = ctx.net, from = ctx.idx[b1], to = ctx.idx[b2], r_pu = r_pu, x_pu = x_pu, b_pu = b_pu, g_pu = g_pu, status = st.status, ratedU = tw.rated_u1[i], ratedS = rating === nothing ? rated_s : rating, ratio = ratio, shift_deg = shift, from_status = st.from_status, to_status = st.to_status)
    _powsybl_count!(ctx.report, "2wt"; built = 1)
  end
  return nothing
end

# --- 2.4 three-winding transformers ----------------------------------------------

function _map_powsybl_three_winding_transformers!(ctx::_PowsyblBuildContext)
  t3 = ctx.tables.three_windings_transformers
  _powsybl_count!(ctx.report, "3wt"; read = length(t3.id))
  for i in eachindex(t3.id)
    id = t3.id[i]
    legs = [(k = k, bus = t3[Symbol("bus_breaker_bus$(k)_id")][i], connected = t3[Symbol("connected$(k)")][i]) for k in 1:3]
    missing = [l.bus for l in legs if !haskey(ctx.idx, l.bus)]
    if !isempty(missing)
      _powsybl_skip!(ctx.report, "3wt", id, "bus $(first(missing)) not built")
      continue
    end
    vn_star = t3.rated_u0[i]
    if !(vn_star > 0.0)
      _powsybl_skip!(ctx.report, "3wt", id, "rated_u0 $(vn_star) is not a voltage")
      continue
    end
    star = id * "_star"
    addBus!(net = ctx.net, busName = star, vn_kV = vn_star, isAux = true)
    ctx.idx[star] = ctx.net.busDict[star]
    ctx.vn[star] = vn_star
    ctx.component[star] = ctx.component[legs[1].bus]
    for l in legs
      k = l.k
      vn_k = ctx.vn[l.bus]
      # every leg is a two-winding transformer from its bus to the star
      # point: impedance at the star side (rated_u0), ratio at the bus side
      rho = t3[Symbol("rho$(k)")][i]
      ratio = _powsybl_2wt_ratio(rho, vn_k, vn_star)
      # the leg's magnetizing admittance sits on the bus side of the leg's
      # pi behind the ideal transformer, like the two-winding case
      g_sh = t3[Symbol("g$(k)_at_current_tap")][i]
      b_sh = t3[Symbol("b$(k)_at_current_tap")][i]
      g_sh = isfinite(g_sh) ? g_sh : 0.0
      b_sh = isfinite(b_sh) ? b_sh : 0.0
      r_pu, x_pu, b_pu, g_pu = toPU_RXBG(r = t3[Symbol("r$(k)_at_current_tap")][i], x = t3[Symbol("x$(k)_at_current_tap")][i], g = 0.0, b = 0.0, v_kv = vn_star, baseMVA = ctx.net.baseMVA)
      kf = vn_star^2 / ratio^2 / vn_k^2
      _powsybl_bus_shunt!(ctx, l.bus, g_sh * kf, b_sh * kf)
      alpha = t3[Symbol("alpha$(k)")][i]
      shift = -(isfinite(alpha) ? alpha : 0.0)
      # a disconnected leg is open at its own terminal, the network side
      st = _powsybl_terminal_status(l.connected, true)
      _addPIModelTrafo_by_idx!(net = ctx.net, from = ctx.idx[l.bus], to = ctx.idx[star], r_pu = r_pu, x_pu = x_pu, b_pu = b_pu, g_pu = g_pu, status = st.status, ratedU = t3[Symbol("rated_u$(k)")][i], ratedS = _powsybl_optional(t3[Symbol("rated_s$(k)")][i]), ratio = ratio, shift_deg = shift, from_status = st.from_status, to_status = st.to_status)
    end
    _powsybl_count!(ctx.report, "3wt"; built = 1)
  end
  return nothing
end

# --- 2.5 generators and the slack -------------------------------------------------

# Reactive limits of a unit: the plain bounds for MIN_MAX, the bounds at the
# active setpoint for CURVE; an unbounded or NaN value is no limit.
function _powsybl_q_limits(kind::AbstractString, min_q::Float64, max_q::Float64, min_q_at::Float64, max_q_at::Float64)
  lo, hi = kind == "CURVE" ? (min_q_at, max_q_at) : (min_q, max_q)
  return (_powsybl_optional(lo), _powsybl_optional(hi))
end

function _powsybl_choose_slacks(ctx::_PowsyblBuildContext)::Dict{Int,Int}
  g = ctx.tables.generators
  # candidates per component: regulating, connected, on a built bus
  best = Dict{Int,Int}()
  order = sortperm(collect(g.id))
  for i in order
    (g.connected[i] && g.voltage_regulator_on[i] && g.target_v[i] > 0.0) || continue
    bus = g.bus_breaker_bus_id[i]
    haskey(ctx.idx, bus) || continue
    comp = ctx.component[bus]
    comp < 0 && continue
    current = get(best, comp, 0)
    if current == 0 || g.max_p[i] > g.max_p[current]
      best[comp] = i
    end
  end
  reasons = Dict{Int,String}(comp => "largest max_p among regulating units" for comp in keys(best))
  for want in ctx.opts.slack_ids
    j = findfirst(==(want), g.id)
    if j === nothing
      _powsybl_notice!(ctx.report, "slack override $(want) is not a generator id")
      continue
    end
    bus = g.bus_breaker_bus_id[j]
    if !(g.connected[j] && haskey(ctx.idx, bus))
      _powsybl_notice!(ctx.report, "slack override $(want) is disconnected or on a bus that was not built")
      continue
    end
    comp = ctx.component[bus]
    best[comp] = j
    reasons[comp] = "override"
  end
  if !ctx.opts.multi_slack && length(best) > 1
    keep = minimum(keys(best))
    for comp in sort!(collect(keys(best)))
      comp == keep && continue
      _powsybl_notice!(ctx.report, "multi_slack is off: component $(comp) keeps no slack (reference stays in component $(keep))")
      delete!(best, comp)
    end
  end
  for comp in sort!(collect(keys(best)))
    i = best[comp]
    push!(ctx.report.slack, (component = comp, generator = g.id[i], bus = g.bus_breaker_bus_id[i], reason = reasons[comp]))
  end
  # components with buses but no candidate
  comps = Set{Int}(c for c in values(ctx.component) if c >= 0)
  for comp in sort!(collect(setdiff(comps, keys(best))))
    ctx.opts.multi_slack || continue
    vsc = ctx.tables.vsc_converter_stations
    has_vsc = any(k -> vsc.voltage_regulator_on[k] && haskey(ctx.idx, vsc.bus_breaker_bus_id[k]) && ctx.component[vsc.bus_breaker_bus_id[k]] == comp, eachindex(vsc.id))
    push!(ctx.report.slack, (component = comp, generator = "", bus = "", reason = has_vsc ? "no regulating generator, a VSC station regulates the voltage" : "no slack candidate"))
  end
  return best
end

function _map_powsybl_generators!(ctx::_PowsyblBuildContext)
  g = ctx.tables.generators
  _powsybl_count!(ctx.report, "generator"; read = length(g.id))
  slack_rows = Set{Int}(values(_powsybl_choose_slacks(ctx)))
  for i in eachindex(g.id)
    id = g.id[i]
    if !g.connected[i]
      _powsybl_skip!(ctx.report, "generator", id, "disconnected")
      continue
    end
    bus = _powsybl_bus(ctx, "generator", id, g.bus_breaker_bus_id[i])
    bus === nothing && continue
    vn = ctx.vn[bus]
    q_lo, q_hi = _powsybl_q_limits(g.reactive_limits_kind[i], g.min_q[i], g.max_q[i], g.min_q_at_target_p[i], g.max_q_at_target_p[i])
    (q_lo === nothing && q_hi === nothing) && _powsybl_notice!(ctx.report, "generator $(id) carries no reactive limits")
    p = isfinite(g.target_p[i]) ? g.target_p[i] : 0.0
    q = isfinite(g.target_q[i]) ? g.target_q[i] : 0.0
    pv = g.voltage_regulator_on[i] && g.target_v[i] > 0.0
    vset = nothing
    if pv
      regulated = g.regulated_bus_id[i]
      own_view = ctx.report.bus_view[bus]
      remote = !isempty(regulated) && regulated != own_view
      if remote && ctx.opts.remote_regulation === :pq
        pv = false
        _powsybl_notice!(ctx.report, "generator $(id) regulates bus $(regulated) (own bus $(own_view)): mapped as PQ with target_q")
      elseif remote && ctx.opts.remote_regulation === :remote
        pv = false
        target_bus = g.regulated_bus_breaker_bus_id[i]
        vn_reg = get(ctx.bus_view_vn, regulated, vn)
        vn_reg > 0.0 || (vn_reg = vn)
        qwide = 10.0 * ctx.net.baseMVA
        push!(ctx.machine_plans, (name = id, bus = String(bus), target_bus = String(target_bus), target_vm_pu = g.target_v[i] / vn_reg, qmin = something(q_lo, -qwide), qmax = something(q_hi, qwide), prosumer_index = length(ctx.net.prosumpsVec) + 1))
      else
        vn_reg = remote ? get(ctx.bus_view_vn, regulated, vn) : vn
        vn_reg > 0.0 || (vn_reg = vn)
        vset = g.target_v[i] / vn_reg
        remote && _powsybl_notice!(ctx.report, "generator $(id) regulates bus $(regulated) (own bus $(own_view)): setpoint $(round(vset; digits = 5)) pu held at the own bus")
      end
    end
    addProsumer!(
      net = ctx.net,
      busName = String(bus),
      type = "GENERATOR",
      p = p,
      q = pv ? nothing : q,
      pMin = _powsybl_optional(g.min_p[i]),
      pMax = _powsybl_optional(g.max_p[i]),
      qMin = q_lo,
      qMax = q_hi,
      referencePri = i in slack_rows ? String(bus) : nothing,
      vm_pu = vset,
      # OpenLoadFlow's default balance (distributed slack proportional to
      # max_p over the generators with a nonzero target; measured on ieee57:
      # four equal shares of 85.5 MW, none for the three units at target_p
      # 0): carried as the participation factor, so
      # power_flow.distributed_slack.p_mode = imported reproduces it
      participationFactor = (p != 0.0 && _powsybl_finite(g.max_p[i]) && g.max_p[i] > 0.0) ? g.max_p[i] : 0.0,
    )
    _powsybl_count!(ctx.report, "generator"; built = 1)
  end
  return nothing
end

# --- 2.6 loads, shunts, static var compensators --------------------------------------

function _map_powsybl_loads!(ctx::_PowsyblBuildContext)
  ld = ctx.tables.loads
  _powsybl_count!(ctx.report, "load"; read = length(ld.id))
  types = Dict{String,Int}()
  for i in eachindex(ld.id)
    id = ld.id[i]
    if !ld.connected[i]
      _powsybl_skip!(ctx.report, "load", id, "disconnected")
      continue
    end
    bus = _powsybl_bus(ctx, "load", id, ld.bus_breaker_bus_id[i])
    bus === nothing && continue
    types[ld.type[i]] = get(types, ld.type[i], 0) + 1
    addProsumer!(net = ctx.net, busName = String(bus), type = "ENERGYCONSUMER", p = isfinite(ld.p0[i]) ? ld.p0[i] : 0.0, q = isfinite(ld.q0[i]) ? ld.q0[i] : 0.0)
    _powsybl_count!(ctx.report, "load"; built = 1)
  end
  isempty(types) || _powsybl_notice!(ctx.report, "load types: " * join(("$(k) $(v)" for (k, v) in sort!(collect(types))), ", "))
  return nothing
end

function _map_powsybl_shunts!(ctx::_PowsyblBuildContext)
  sh = ctx.tables.shunt_compensators
  _powsybl_count!(ctx.report, "shunt"; read = length(sh.id))
  for i in eachindex(sh.id)
    id = sh.id[i]
    if !sh.connected[i]
      _powsybl_skip!(ctx.report, "shunt", id, "disconnected")
      continue
    end
    bus = _powsybl_bus(ctx, "shunt", id, sh.bus_breaker_bus_id[i])
    bus === nothing && continue
    vn = ctx.vn[bus]
    g = isfinite(sh.g[i]) ? sh.g[i] : 0.0
    b = isfinite(sh.b[i]) ? sh.b[i] : 0.0
    sh.model_type[i] == "LINEAR" || _powsybl_notice!(ctx.report, "shunt $(id) has model $(sh.model_type[i]): g and b taken as delivered")
    # siemens at the bus nominal voltage: MW and MVar at 1 pu
    addShuntMatpower!(net = ctx.net, busName = String(bus), Gs = g * vn^2, Bs = b * vn^2)
    _powsybl_count!(ctx.report, "shunt"; built = 1)
  end
  return nothing
end

function _map_powsybl_svcs!(ctx::_PowsyblBuildContext)
  sv = ctx.tables.static_var_compensators
  _powsybl_count!(ctx.report, "svc"; read = length(sv.id))
  qwide = 10.0 * ctx.net.baseMVA
  for i in eachindex(sv.id)
    id = sv.id[i]
    if !sv.connected[i]
      _powsybl_skip!(ctx.report, "svc", id, "disconnected")
      continue
    end
    bus = _powsybl_bus(ctx, "svc", id, sv.bus_breaker_bus_id[i])
    bus === nothing && continue
    vn = ctx.vn[bus]
    # Q = b * U^2 with b in siemens and U the nominal voltage in kV
    q_lo = isfinite(sv.b_min[i]) ? sv.b_min[i] * vn^2 : -qwide
    q_hi = isfinite(sv.b_max[i]) ? sv.b_max[i] * vn^2 : qwide
    if max(abs(q_lo), abs(q_hi)) > qwide
      _powsybl_notice!(ctx.report, "svc $(id): Q range $(round(q_lo; digits = 1))..$(round(q_hi; digits = 1)) MVar from b_min/b_max clamped to ±$(qwide)")
      q_lo = max(q_lo, -qwide)
      q_hi = min(q_hi, qwide)
    end
    voltage = sv.regulation_mode[i] == "VOLTAGE" && sv.regulating[i] && sv.target_v[i] > 0.0
    addProsumer!(net = ctx.net, busName = String(bus), type = "STATICVARCOMPENSATOR", p = 0.0, q = voltage ? nothing : (isfinite(sv.target_q[i]) ? sv.target_q[i] : 0.0), qMin = q_lo, qMax = q_hi, vm_pu = voltage ? sv.target_v[i] / vn : nothing)
    _powsybl_count!(ctx.report, "svc"; built = 1)
  end
  return nothing
end

# --- 2.7 dangling lines and tie lines --------------------------------------------------

# A dangling line as a branch from its bus to an auxiliary bus; PowSyBl's
# model keeps the whole shunt admittance at the network terminal.
function _powsybl_add_dangling_branch!(ctx::_PowsyblBuildContext, bus::AbstractString, aux::AbstractString, r::Float64, x::Float64, g::Float64, b::Float64, connected::Bool, rating)
  st = _powsybl_terminal_status(connected, true)
  _powsybl_add_pi_line!(ctx, bus, aux, r, x, 0.0, 0.0, st, rating)
  _powsybl_bus_shunt!(ctx, bus, g, b)
  return nothing
end

function _map_powsybl_dangling_lines!(ctx::_PowsyblBuildContext)
  dl = ctx.tables.dangling_lines
  tl = ctx.tables.tie_lines
  _powsybl_count!(ctx.report, "dangling_line"; read = length(dl.id))
  _powsybl_count!(ctx.report, "tie_line"; read = length(tl.id))
  row = Dict{String,Int}(dl.id[i] => i for i in eachindex(dl.id))
  # tie lines: the two paired dangling lines in series through the X-node
  for i in eachindex(tl.id)
    id = tl.id[i]
    k1 = get(row, tl.dangling_line1_id[i], nothing)
    k2 = get(row, tl.dangling_line2_id[i], nothing)
    if k1 === nothing || k2 === nothing
      _powsybl_skip!(ctx.report, "tie_line", id, "dangling line $(k1 === nothing ? tl.dangling_line1_id[i] : tl.dangling_line2_id[i]) not in dangling_lines")
      continue
    end
    b1 = _powsybl_bus(ctx, "tie_line", id, dl.bus_breaker_bus_id[k1])
    b1 === nothing && continue
    b2 = _powsybl_bus(ctx, "tie_line", id, dl.bus_breaker_bus_id[k2])
    b2 === nothing && continue
    xnode = (isempty(tl.pairing_key[i]) ? id : tl.pairing_key[i]) * "_xnode"
    if !haskey(ctx.idx, xnode)
      addBus!(net = ctx.net, busName = xnode, vn_kV = ctx.vn[b1], isAux = true)
      ctx.idx[xnode] = ctx.net.busDict[xnode]
      ctx.vn[xnode] = ctx.vn[b1]
      ctx.component[xnode] = ctx.component[b1]
    end
    _powsybl_add_dangling_branch!(ctx, b1, xnode, dl.r[k1], dl.x[k1], dl.g[k1], dl.b[k1], tl.connected1[i] && dl.connected[k1], _powsybl_branch_rating(ctx, dl.id[k1], ctx.vn[b1], ctx.vn[b1]))
    _powsybl_add_dangling_branch!(ctx, b2, xnode, dl.r[k2], dl.x[k2], dl.g[k2], dl.b[k2], tl.connected2[i] && dl.connected[k2], _powsybl_branch_rating(ctx, dl.id[k2], ctx.vn[b2], ctx.vn[b2]))
    p0 = dl.p0[k1] + dl.p0[k2]
    q0 = dl.q0[k1] + dl.q0[k2]
    (isfinite(p0) && isfinite(q0) && (p0 != 0.0 || q0 != 0.0)) && addProsumer!(net = ctx.net, busName = xnode, type = "ENERGYCONSUMER", p = p0, q = q0)
    _powsybl_count!(ctx.report, "tie_line"; built = 1)
  end
  for i in eachindex(dl.id)
    id = dl.id[i]
    if dl.paired[i]
      _powsybl_skip!(ctx.report, "dangling_line", id, "paired into tie line $(dl.tie_line_id[i])")
      continue
    end
    if !dl.connected[i]
      _powsybl_skip!(ctx.report, "dangling_line", id, "disconnected")
      continue
    end
    bus = _powsybl_bus(ctx, "dangling_line", id, dl.bus_breaker_bus_id[i])
    bus === nothing && continue
    aux = id * "_boundary"
    addBus!(net = ctx.net, busName = aux, vn_kV = ctx.vn[bus], isAux = true)
    ctx.idx[aux] = ctx.net.busDict[aux]
    ctx.vn[aux] = ctx.vn[bus]
    ctx.component[aux] = ctx.component[bus]
    _powsybl_add_dangling_branch!(ctx, bus, aux, dl.r[i], dl.x[i], dl.g[i], dl.b[i], true, _powsybl_branch_rating(ctx, id, ctx.vn[bus], ctx.vn[bus]))
    p0 = isfinite(dl.p0[i]) ? dl.p0[i] : 0.0
    q0 = isfinite(dl.q0[i]) ? dl.q0[i] : 0.0
    (p0 != 0.0 || q0 != 0.0) && addProsumer!(net = ctx.net, busName = aux, type = "ENERGYCONSUMER", p = p0, q = q0)
    _powsybl_count!(ctx.report, "dangling_line"; built = 1)
  end
  return nothing
end

# --- 2.8 HVDC -----------------------------------------------------------------------------

function _powsybl_station(ctx::_PowsyblBuildContext, id::AbstractString)
  vsc = ctx.tables.vsc_converter_stations
  k = findfirst(==(id), vsc.id)
  k === nothing || return (kind = :vsc, row = k)
  lcc = ctx.tables.lcc_converter_stations
  k = findfirst(==(id), lcc.id)
  k === nothing || return (kind = :lcc, row = k)
  return nothing
end

# One converter station as a fixed PCC injection. `p_ac` is the active
# power the station exchanges with the AC grid in load convention
# (positive draws from the grid).
function _powsybl_add_station!(ctx::_PowsyblBuildContext, station, id::AbstractString, p_ac::Float64)
  if station.kind === :vsc
    vsc = ctx.tables.vsc_converter_stations
    k = station.row
    bus = _powsybl_bus(ctx, "hvdc", id, vsc.bus_breaker_bus_id[k])
    bus === nothing && return false
    vn = ctx.vn[bus]
    q_lo, q_hi = _powsybl_q_limits(vsc.reactive_limits_kind[k], vsc.min_q[k], vsc.max_q[k], vsc.min_q_at_target_p[k], vsc.max_q_at_target_p[k])
    voltage = vsc.voltage_regulator_on[k] && vsc.target_v[k] > 0.0
    # generator convention on the prosumer: injection positive
    addProsumer!(net = ctx.net, busName = String(bus), type = "GENERATOR", p = -p_ac, q = voltage ? nothing : (isfinite(vsc.target_q[k]) ? vsc.target_q[k] : 0.0), qMin = q_lo, qMax = q_hi, vm_pu = voltage ? vsc.target_v[k] / vn : nothing)
  else
    lcc = ctx.tables.lcc_converter_stations
    k = station.row
    bus = _powsybl_bus(ctx, "hvdc", id, lcc.bus_breaker_bus_id[k])
    bus === nothing && return false
    pf = lcc.power_factor[k]
    # a line-commutated converter always absorbs reactive power
    q_abs = (isfinite(pf) && 0.0 < pf <= 1.0) ? abs(p_ac) * tan(acos(pf)) : 0.0
    addProsumer!(net = ctx.net, busName = String(bus), type = "ENERGYCONSUMER", p = p_ac, q = q_abs)
  end
  return true
end

function _map_powsybl_hvdc!(ctx::_PowsyblBuildContext)
  hv = ctx.tables.hvdc_lines
  _powsybl_count!(ctx.report, "hvdc"; read = length(hv.id))
  isempty(hv.id) && return nothing
  ctx.opts.hvdc_mode === :fixed_injection || throw(ArgumentError("powsybl import: hvdc_mode $(ctx.opts.hvdc_mode) is not implemented for PowSyBl bundles in this release; use fixed_injection"))
  for i in eachindex(hv.id)
    id = hv.id[i]
    if !(hv.connected1[i] && hv.connected2[i])
      _powsybl_skip!(ctx.report, "hvdc", id, "a converter station is disconnected")
      continue
    end
    s1 = _powsybl_station(ctx, hv.converter_station1_id[i])
    s2 = _powsybl_station(ctx, hv.converter_station2_id[i])
    if s1 === nothing || s2 === nothing
      _powsybl_skip!(ctx.report, "hvdc", id, "converter station $(s1 === nothing ? hv.converter_station1_id[i] : hv.converter_station2_id[i]) not found")
      continue
    end
    mode = hv.converters_mode[i]
    rect_first = mode != "SIDE_1_INVERTER_SIDE_2_RECTIFIER"
    rect, inv = rect_first ? (s1, s2) : (s2, s1)
    rect_id, inv_id = rect_first ? (hv.converter_station1_id[i], hv.converter_station2_id[i]) : (hv.converter_station2_id[i], hv.converter_station1_id[i])
    loss(st) = st.kind === :vsc ? ctx.tables.vsc_converter_stations.loss_factor[st.row] : ctx.tables.lcc_converter_stations.loss_factor[st.row]
    lf_rect = isfinite(loss(rect)) ? loss(rect) : 0.0
    lf_inv = isfinite(loss(inv)) ? loss(inv) : 0.0
    target = isfinite(hv.target_p[i]) ? hv.target_p[i] : 0.0
    # measured on four_substations against OpenLoadFlow: the rectifier
    # draws target_p / (1 - lf / 100), the inverter delivers
    # target_p * (1 - lf / 100); the DC line loss is not modelled
    p_rect = target / (1.0 - lf_rect / 100.0)
    p_inv = target * (1.0 - lf_inv / 100.0)
    ok = _powsybl_add_station!(ctx, rect, id, p_rect)
    ok = _powsybl_add_station!(ctx, inv, id, -p_inv) && ok
    ok || continue
    _powsybl_notice!(ctx.report, "hvdc $(id): $(rect_id) draws $(round(p_rect; digits = 3)) MW, $(inv_id) delivers $(round(p_inv; digits = 3)) MW (fixed injections)")
    _powsybl_count!(ctx.report, "hvdc"; built = 1)
  end
  return nothing
end
