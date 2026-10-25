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

# file: src/adapters/cgmes/cgmes_mapping.jl
# purpose: layer 5 of the CGMES importer, map CIM objects onto the
# Sparlectra Net (Stage 1: bus-branch, lines, 2W/3W transformers with fixed
# SSH tap positions, shunts, loads, machines, equivalent injections, slack
# selection). Short-circuit source data is always harvested into
# `CGMESShortCircuitData` (read, not evaluated).

"""
Typed short-circuit source data harvested during import.
Each entry carries the CIM mRID, the object name and, where resolvable,
the Sparlectra bus name. Values stay in CGMES units (Ω, S, A, pu on machine
base); nothing is evaluated in Stage 1. `importCGMES` keeps the records of
the sources it mapped onto the network (`ExternalNetworkInjection`,
`SynchronousMachine`, `EquivalentInjection`, `AsynchronousMachine`); a
skipped one (terminal disconnected, out of service) feeds no fault current.

`power_electronics_connections` holds one record per
`PowerElectronicsConnection` (CGMES 3.0) that the import mapped (only with
`map_power_electronics = true`), with `maxIFault_pu` (the maximum fault
current in per unit of the rated current), `ratedS_MVA` and `ratedU_kV`.
`maxIFault` is a CIM17 attribute that no CGMES 3.0 profile defines; it is
read where a delivery carries it, and is `nothing` otherwise.

A `SynchronousMachine` record carries `converter_unit`, the class name of
its generating unit when that is a `WindGeneratingUnit` or a
`SolarGeneratingUnit` (`nothing` otherwise), and for such a machine
`pmax_MW` = `GeneratingUnit.maxOperatingP`. `runShortCircuit!` takes such a
machine as a converter current source, not as a voltage source behind
`x''d`. `windGenUnitType` holds the delivered `WindGeneratingUnit.windGenUnitType`
(`onshore`, `offshore`, `nothing` when absent or not a wind unit); the CGMES
export writes the unit class and this value back.
"""
struct CGMESShortCircuitData
  external_network_injections::Vector{NamedTuple}
  synchronous_machines::Vector{NamedTuple}
  ac_line_segments::Vector{NamedTuple}
  transformer_ends::Vector{NamedTuple}
  equivalent_injections::Vector{NamedTuple}
  asynchronous_machines::Vector{NamedTuple}
  power_electronics_connections::Vector{NamedTuple}
end

# one row of _svTapStepFits (cgmes_report.jl): a tap changer whose SvTapStep
# fits the SV flows of its transformer clearly worse than another step
const _SvTapFitRow = NamedTuple{(:transformer, :tap_changer, :sv_step, :best_step, :sv_mismatch_mva, :best_mismatch_mva),Tuple{String,String,Int,Int,Float64,Float64}}

# What the import did with the DY profile, kept for the import report
# (dynamicsReportLines): the DY files of the delivery, whether they were
# read (cgmes_import.read_dynamics), the units that got DynamicData and the
# per-class and per-reason counts of _mapDynamics!. `nothing` on the result
# when the delivery carries no DY profile.
struct CGMESDynamicsSummary
  files::Vector{String}
  read::Bool
  units::Vector{String}
  classes::Vector{NamedTuple{(:class, :mapped, :not_mapped, :no_eq),Tuple{Symbol,Int,Int,Int}}}
  reasons::Vector{NamedTuple{(:class, :count, :reason),Tuple{Symbol,Int,String}}}
end

"""
Result of `importCGMES`: the `Net`, the merged `CGMESStore`, the bus-branch
topology, the harvested short-circuit data, the selected slack bus and the
importer's skip/notice messages. `branch_side_of_terminal` (CGMES Terminal
mRID → `(branchIdx, :from|:to)`) and `skipped_equipment` (mRIDs that were not
mapped) provide the provenance `compareWithSV` needs for its `SvPowerFlow`
comparison. `no_sv_buses` lists the created buses whose TopologicalNode has no
usable `SvVoltage` (none, or one with `v <= 0`); they start at the voltage
of a neighbour with an SV voltage, or at `1.0 pu / 0°` where there is none,
which weakens any comparison against the delivery's SV state
(fixed-reference self-check, `compareWithSV`).
`sv_tap_step_fits` lists the tap changers whose `SvPowerFlow` fits another
step clearly better than their `SvTapStep` (transformer, tap changer, both
steps and both mismatches in MVA); the import counts them in one notice.
`dynamics` says what the import did with a DY profile (`nothing` without
one): its files, whether `cgmes_import.read_dynamics` read them, the units
that got `DynamicData` and the counts per DY class and reason;
`dynamicsReportLines` renders it.

`element_by_mrid` maps the mRID of every CGMES object the importer turned
into a network element to `(kind, index)`, recorded where the element is
created:

- `:branch`: `ACLineSegment`, `SeriesCompensator` (also a boundary line built
  as a ratio branch), `EquivalentBranch` (across two voltage levels a branch
  with ratio 1); index into `net.branchVec`.
- `:transformer`: a two-winding `PowerTransformer` under its own mRID; a
  three-winding one under the mRID of each `PowerTransformerEnd`, one star
  leg each (the `PowerTransformer` mRID itself is not a key); index into
  `net.branchVec`.
- `:shunt`: `LinearShuntCompensator`, `NonlinearShuntCompensator`, and the
  charging of a line whose two terminals sit on one TopologicalNode under
  `"<line mRID>#charging"`; index into `net.shuntVec`.
- `:load`: `EnergyConsumer`, `ConformLoad`, `NonConformLoad`,
  `StationSupply`, `AsynchronousMachine`; index into `net.prosumpsVec`.
- `:generator`: `SynchronousMachine`, `StaticVarCompensator`,
  `EnergySource` when `map_energy_source = true`, and
  `PowerElectronicsConnection` (CGMES 3.0) when
  `map_power_electronics = true`; index into `net.prosumpsVec`.
- `:external_grid`: `ExternalNetworkInjection`, `EquivalentInjection`,
  `VsConverter`, `CsConverter`; index into `net.prosumpsVec`.
- `:link`: a closed switch between two buses; index into `net.linkVec`.
- `:tap_controller`: a `RatioTapChanger` or phase tap changer that became a
  controller (`tap_control = true`); index into `_tap_controllers(net)`.
  A tap changer that joins another one's group as follower has no entry.

Skipped objects have no entry; they are in `skipped_equipment`, which also
holds every `EnergySource` of an import with `map_energy_source = false`
and every `PowerElectronicsConnection` of an import with
`map_power_electronics = false`.
Buses are addressed by their TopologicalNode (`topo.bus_name`), not here.
"""
struct CGMESImportResult
  net::Sparlectra.Net
  store::CGMESStore
  topo::CGMESTopology
  shortcircuit::CGMESShortCircuitData
  slack_bus::String
  messages::Vector{String}
  branch_side_of_terminal::Dict{String,Tuple{Int,Symbol}}
  skipped_equipment::Set{String}
  no_sv_buses::Vector{String}
  element_by_mrid::Dict{String,Tuple{Symbol,Int}}
  sv_tap_step_fits::Vector{_SvTapFitRow}
  dynamics::Union{Nothing,CGMESDynamicsSummary}
end

# mutable importer bookkeeping threaded through the mapping helpers.
# vset_min_pu/vset_max_pu ride along here rather than through every mapping
# signature; see _voltageSetpoint for what the band does.
# split_sides: boundary CN → (source file → per-side bus name), filled by
# _detectNonCancellingBoundarySides! before the element mapping runs.
struct _MapCtx
  messages::Vector{String}
  branch_side::Dict{String,Tuple{Int,Symbol}}
  skipped::Set{String}
  tap_control::Bool
  machine_control::Bool
  ignore_connected::Bool
  vset_min_pu::Float64
  vset_max_pu::Float64
  split_sides::Dict{String,Dict{String,String}}
  # machine remote-voltage controllers planned by _mapInjections! and attached
  # by importCGMES only after bus types and isolation are final
  machine_ctrl_plans::Vector{NamedTuple}
  # first-seen counters for structural identity keys, keyed by (kind, bus):
  # shared across all mapping passes so capture and a later export assign the
  # same <k> even when equipment of one class arrives from different loops
  eqp_counters::Dict{Tuple{String,String},Int}
  # cgmes_import.placeholder_guards = strict: a suspected placeholder value
  # aborts the import instead of being skipped with a warning
  strict_guards::Bool
  # prosumer index per source mRID, filled during the pending flush; lets
  # post-build steps (HVDC pair attachment) address the exact injection a
  # CGMES object became instead of guessing by bus and value
  prosumer_by_mrid::Dict{String,Int}
  # tap controller per TapChangerControl mRID (#322): several tap changers
  # referencing ONE RegulatingControl are a regulated parallel-transformer
  # GROUP; the first one becomes the master controller, later ones join as
  # followers instead of spawning a second, fighting controller
  tap_group_masters::Dict{String,Any}
  # CGMES mRID -> (kind, index) of the element it became, filled where the
  # element is created; becomes CGMESImportResult.element_by_mrid
  element_by_mrid::Dict{String,Tuple{Symbol,Int}}
  # TapChanger mRID -> the controller it became. The position in
  # _tap_controllers(net) is read only once all controllers exist: a
  # controller is appended to the FIRST transformer on its bus pair, so a
  # later parallel unit can move the controllers created before it
  tap_controller_by_mrid::Vector{Pair{String,Sparlectra.PowerTransformerControl}}
  # cgmes_import.map_energy_source: EnergySource objects become PQ
  # injections; off, they are counted and listed as not mapped
  map_energy_source::Bool
  # cgmes_import.map_power_electronics: PowerElectronicsConnection objects
  # become PQ or voltage-regulating injections; off, they are counted and
  # listed as not mapped
  map_power_electronics::Bool
  # bus index of every three-winding star point, in creation order; its start
  # voltage is set by _startStarBuses! once every branch exists
  star_buses::Vector{Int}
  # cgmes_import.equivalent_branch_min_z_pu: an EquivalentBranch below this
  # impedance becomes a bus link instead of a branch
  eb_min_z_pu::Float64
  # mapped units per class whose P limit was a placeholder read as no limit
  # (_pLimitsRead!); one notice in _mapInjections! counts them
  p_limit_placeholders::Dict{String,Int}
  # cgmes_import.equivalent_injection_as_source: an EquivalentInjection on
  # the SV angle reference node or on a boundary node becomes a
  # voltage-regulating source at the SV voltage of its node (the slack when
  # it holds the SV angle reference)
  ei_as_source::Bool
end
_MapCtx(tap_control::Bool = false, ignore_connected::Bool = false; machine_control::Bool = false, vset_min_pu::Float64 = CGMES_VSET_MIN_PU, vset_max_pu::Float64 = CGMES_VSET_MAX_PU, strict_guards::Bool = false, map_energy_source::Bool = false, map_power_electronics::Bool = false, eb_min_z_pu::Float64 = 1e-6, ei_as_source::Bool = false) =
  _MapCtx(String[], Dict{String,Tuple{Int,Symbol}}(), Set{String}(), tap_control, machine_control, ignore_connected, vset_min_pu, vset_max_pu, Dict{String,Dict{String,String}}(), NamedTuple[], Dict{Tuple{String,String},Int}(), strict_guards, Dict{String,Int}(), Dict{String,Any}(), Dict{String,Tuple{Symbol,Int}}(), Pair{String,Sparlectra.PowerTransformerControl}[], map_energy_source, map_power_electronics, Int[], eb_min_z_pu, Dict{String,Int}(), ei_as_source)

# SSH `Terminal.connected` through the ctx switch: `ignore_connected = true`
# treats everything as connected (diagnostic mode for snapshots whose SSH
# flags contradict their SV state)
_conn(ctx::_MapCtx, flag::Bool) = ctx.ignore_connected || flag

# SSH `Equipment.inService`, CGMES 3.0 moved the operational status onto the
# equipment itself; 2.4.15 deliveries do not carry it (absent → in service).
# This is NOT redundant with Terminal.connected: ReliCapGrid parks
# out-of-service generators with connected=true, inService=false, no SvVoltage
# and a placeholder regulation target. Ignoring the flag imported 413 MW of
# phantom generation into Svedala alone, the single largest reason the solved
# corridor flows missed the delivery's SV state. Honors ignore_connected for
# consistency with _conn (the override means "treat everything as live").
_inService(ctx::_MapCtx, o::CIMObject)::Bool = ctx.ignore_connected || something(boolval(o, :inService), true)

# --- tap changers (fixed operating point in Stage 1) ------------------------
#
# Step source order: SvTapStep position (the solved state, present in every
# conformity data set and in real exchanges) → SSH `step` (initial operating
# point) → the middle position (_tapMiddleStep: neutralStep, the band middle
# without one), named once per import (_mapTransformers!). The MicroGrid BaseCase SV was produced with active tap
# control, so SvTapStep differs from SSH for regulated units; a fixed-tap
# import must use the solved position to reproduce SV. Stage 2 controllers
# will start from SSH instead.

"""SvTapStep positions: tap-changer mRID → solved position."""
function _svTapPositions(store::CGMESStore)::Dict{String,Float64}
  out = Dict{String,Float64}()
  for sv in objectsOf(store, :SvTapStep)
    tc = get(sv.refs, :TapChanger, nothing)
    tc === nothing && continue
    p = num(sv, :position)
    p === nothing && continue
    out[tc] = p
  end
  return out
end

# Middle position of a changer whose delivery carries no step: its
# neutralStep, the band middle without one, 0 without a band
function _tapMiddleStep(tc::CIMObject)::Float64
  n = num(tc, :neutralStep)
  n === nothing || return n
  lo, hi = num(tc, :lowStep), num(tc, :highStep)
  (lo === nothing || hi === nothing) && return 0.0
  return Float64(Sparlectra._tap_band_middle(lo, hi))
end

# true when the delivery carries a step for `tc` (SvTapStep or SSH step)
_tapStepGiven(tc::CIMObject, svsteps::Dict{String,Float64})::Bool = haskey(svsteps, tc.mrid) || num(tc, :step) !== nothing

_tapStep(tc::CIMObject, svsteps::Dict{String,Float64}) = something(get(svsteps, tc.mrid, nothing), num(tc, :step), _tapMiddleStep(tc))

# linear step formula of a RatioTapChanger: 1 + (step - neutral) * increment
function _ratioTapLinearCorrection(rtc::CIMObject, step::Float64)::Float64
  neutral = num(rtc, :neutralStep, 0.0)
  incr = num(rtc, :stepVoltageIncrement, 0.0)
  return 1.0 + (step - neutral) * incr / 100.0
end

# (step, ratio) points of the RatioTapChangerTable of `rtc`, sorted by step;
# `nothing` when the changer references no table, empty when the table has
# no point with both values
function _ratioTapTable(store::CGMESStore, rtc::CIMObject)::Union{Nothing,Vector{Tuple{Float64,Float64}}}
  tmrid = get(rtc.refs, :RatioTapChangerTable, nothing)
  tmrid === nothing && return nothing
  pts = Tuple{Float64,Float64}[]
  for pt in objectsOf(store, :RatioTapChangerTablePoint)
    get(pt.refs, :RatioTapChangerTable, "") == tmrid || continue
    s = num(pt, :step)
    r = num(pt, :ratio)
    (s === nothing || r === nothing) && continue
    push!(pts, (s, r))
  end
  sort!(pts; by = first)
  return pts
end

# true when a point of the RatioTapChangerTable of `rtc` carries a non-zero
# r, x, g or b percent correction (applied to the impedance by _tapTableImpedanceFactors)
function _ratioTapTableHasImpedance(store::CGMESStore, rtc::CIMObject)::Bool
  tmrid = get(rtc.refs, :RatioTapChangerTable, nothing)
  tmrid === nothing && return false
  for pt in objectsOf(store, :RatioTapChangerTablePoint)
    get(pt.refs, :RatioTapChangerTable, "") == tmrid || continue
    any(!iszero(something(num(pt, attr), 0.0)) for attr in (:r, :x, :g, :b)) && return true
  end
  return false
end

# Ratio correction of `rtc` at `step` and whether it came from the table.
# A table that covers the step gives the ratio of its point: the table is the
# source of truth, its values need not follow the linear formula and need
# not be uniform. A fractional step (continuous SvTapStep) or a step without
# its own point inside the table range is interpolated between the
# neighbouring points. Without a table, or with a step outside it, the
# linear formula applies.
function _ratioTapCorrectionAt(rtc::CIMObject, table::Union{Nothing,Vector{Tuple{Float64,Float64}}}, step::Float64)::Tuple{Float64,Bool}
  (table === nothing || isempty(table) || step < table[1][1] || step > table[end][1]) && return (_ratioTapLinearCorrection(rtc, step), false)
  k = findfirst(p -> p[1] >= step, table)
  s2, r2 = table[k]
  s2 == step && return (r2, true)
  s1, r1 = table[k - 1]
  return (r1 + (r2 - r1) * (step - s1) / (s2 - s1), true)
end

"""
Tabular phase-tap model of `ptc` (#294 point 2): build a `:tabular`
`PhaseTapChangerModel` from the referenced `PhaseTapChangerTable` rows.
Returns `(model, impedance_note)` or `nothing` when the table is
unresolvable or empty. `impedance_note` flags rows with non-zero r/x/g/b
percent corrections; the model carries ratio and angle, the corrections of
the step reach the branch through `_tapTableImpedanceFactors`.
"""
function _phaseTapTableModel(store::CGMESStore, ptc::CIMObject, step::Int)
  tmrid = get(ptc.refs, :PhaseTapChangerTable, nothing)
  tmrid === nothing && return nothing
  rows = Sparlectra.TapTablePoint[]
  impedance_note = false
  for pt in objectsOf(store, :PhaseTapChangerTablePoint)
    get(pt.refs, :PhaseTapChangerTable, "") == tmrid || continue
    s = num(pt, :step)
    s === nothing && continue
    push!(rows, Sparlectra.TapTablePoint(step = Int(round(s)), ratio = num(pt, :ratio, 1.0), angle_deg = num(pt, :angle, 0.0)))
    if any(!iszero(something(num(pt, attr), 0.0)) for attr in (:r, :x, :g, :b))
      impedance_note = true
    end
  end
  isempty(rows) && return nothing
  sort!(rows; by = r -> r.step)
  # a step outside the table (defensive: inconsistent SSH/SV data) is clamped
  # to the nearest table row rather than erroring the whole import
  step_clamped = clamp(step, rows[1].step, rows[end].step)
  haskey_step = any(r -> r.step == step_clamped, rows)
  haskey_step || (step_clamped = rows[argmin([abs(r.step - step_clamped) for r in rows])].step)
  neutral = Int(round(num(ptc, :neutralStep, Float64(rows[1].step))))
  any(r -> r.step == neutral, rows) || (neutral = rows[1].step)
  model = Sparlectra.PhaseTapChangerModel(kind = :tabular, step = step_clamped, neutralStep = neutral, convention = :direct_regulating_vector, table = rows)
  return (model = model, impedance_note = impedance_note)
end

# A CGMES PhaseTapChangerLinear as a :tabular model with generated points
# (0.20.0, no fourth kind): angle (k - neutral) * stepPhaseShiftIncrement,
# ratio 1 at every step. xMin/xMax are given in ohm on the end base, which
# this builder does not know, so the points carry no x_pu (the branch keeps
# its series reactance).
function _phaseTapLinearModel(ptc::CIMObject, step::Int)
  neutral = Int(round(num(ptc, :neutralStep, 0.0)))
  low = Int(round(num(ptc, :lowStep, Float64(min(step, neutral)))))
  high = Int(round(num(ptc, :highStep, Float64(max(step, neutral)))))
  low <= neutral <= high || (low = min(low, neutral); high = max(high, neutral))
  inc = num(ptc, :stepPhaseShiftIncrement, 0.0)
  rows = Sparlectra.TapTablePoint[Sparlectra.TapTablePoint(step = k, ratio = 1.0, angle_deg = (k - neutral) * inc) for k in low:high]
  return Sparlectra.PhaseTapChangerModel(kind = :tabular, step = clamp(step, low, high), neutralStep = neutral, convention = :direct_regulating_vector, table = rows)
end

function _phaseTapRatioShift(store::CGMESStore, ptc::CIMObject, svsteps::Dict{String,Float64})::Tuple{Float64,Float64}
  step = Int(round(_tapStep(ptc, svsteps)))
  neutral = Int(round(num(ptc, :neutralStep, 0.0)))
  low = Int(round(num(ptc, :lowStep, Float64(min(step, neutral)))))
  high = Int(round(num(ptc, :highStep, Float64(max(step, neutral)))))
  # the class name maps to the kind through the one table (tap_changer_kinds.jl)
  kind = Sparlectra.tap_changer_kind(:cgmes, String(ptc.class)).kind
  if kind === :tabular
    built = ptc.class == :PhaseTapChangerLinear ? (model = _phaseTapLinearModel(ptc, step), impedance_note = false) : _phaseTapTableModel(store, ptc, step)
    built === nothing && return (1.0, 0.0)
    res = Sparlectra.calcPhaseTapAngleRatio(built.model)
    return (res.effective_ratio, res.effective_shift_deg)
  end
  vincr = num(ptc, :voltageStepIncrement)
  model = Sparlectra.PhaseTapChangerModel(;
    kind = kind,
    step = step,
    lowStep = low,
    highStep = high,
    neutralStep = neutral,
    # CGMES voltageStepIncrement is in percent; PhaseTapChangerModel expects
    # the per-step fraction (see calcPhaseTapFraction)
    voltage_step_increment = vincr === nothing ? nothing : vincr / 100.0,
    winding_connection_angle_deg = num(ptc, :windingConnectionAngle),
    # the branch's complex ratio t applies the regulating vector directly on
    # the tapped winding (validated against MicroGrid SV, SvTapStep=16)
    convention = :direct_regulating_vector,
  )
  res = Sparlectra.calcPhaseTapAngleRatio(model)
  return (res.effective_ratio, res.effective_shift_deg)
end

# The table point of a tabular tap changer at its step: the row of that
# step, else (inconsistent SSH/SV data) the row nearest to it, the same
# clamping `_phaseTapTableModel` applies. `nothing` without a table or rows.
function _tapTablePointAtStep(store::CGMESStore, tc::CIMObject, step::Float64)::Union{Nothing,CIMObject}
  tabkey, pointcls = tc.class == :RatioTapChanger ? (:RatioTapChangerTable, :RatioTapChangerTablePoint) : (:PhaseTapChangerTable, :PhaseTapChangerTablePoint)
  tmrid = get(tc.refs, tabkey, nothing)
  tmrid === nothing && return nothing
  best = nothing
  bestd = Inf
  for pt in objectsOf(store, pointcls)
    get(pt.refs, tabkey, "") == tmrid || continue
    s = num(pt, :step)
    s === nothing && continue
    d = abs(s - step)
    d < bestd && (best = pt; bestd = d)
  end
  return best
end

# Reactance curve of a linear or non-linear phase tap changer
# (PhaseTapChangerLinear/NonLinear.xMin, xMax): the reactance follows a
# u-shaped curve over the shift angle, xMin at the neutral step and xMax at
# the step with the largest angle,
#   x(a) = xMin + (xMax - xMin) * (sin(a/2) / sin(amax/2))^2,
# with amax the larger angle of lowStep and highStep. Returns
# (xmin, xmax, amax), or nothing when the changer carries no curve (other
# class, attribute missing, xMin <= 0, xMin == xMax, no angle range) or when
# xMax/xMin leaves the plausibility band of the tap corrections
# (_phaseTapXCurveNotice! reports that case).
function _phaseTapXCurve(store::CGMESStore, tc::CIMObject)
  tc.class in (:PhaseTapChangerLinear, :PhaseTapChangerSymmetrical, :PhaseTapChangerAsymmetrical) || return nothing
  xmin, xmax = num(tc, :xMin), num(tc, :xMax)
  (xmin === nothing || xmax === nothing || xmin <= 0.0 || xmax == xmin) && return nothing
  CGMES_TAP_CORR_MIN <= xmax / xmin <= CGMES_TAP_CORR_MAX || return nothing
  amax = 0.0
  for k in (:lowStep, :highStep)
    st = num(tc, k)
    st === nothing && continue
    amax = max(amax, abs(last(_phaseTapRatioShift_atstep(store, tc, Int(round(st))))))
  end
  amax > 0.0 || return nothing
  return (xmin = xmin, xmax = xmax, amax = amax)
end

# factor of the reactance curve at `step` relative to xMin, 1.0 without a
# curve; the series reactance of the end stays as delivered at the neutral
# step
function _phaseTapXCurveFactor(store::CGMESStore, tc::CIMObject, step::Float64)::Float64
  c = _phaseTapXCurve(store, tc)
  c === nothing && return 1.0
  a = last(_phaseTapRatioShift_atstep(store, tc, Int(round(step))))
  return (c.xmin + (c.xmax - c.xmin) * (sind(a / 2) / sind(c.amax / 2))^2) / c.xmin
end

# impedance factors (r, x, g, b) of one tap changer at `step`: the table
# point of a ratio or tabular changer, the reactance curve of a linear or
# non-linear phase changer
function _tapChangerImpedanceFactors(store::CGMESStore, tc::CIMObject, step::Float64)::NTuple{4,Float64}
  tc.class in (:RatioTapChanger, :PhaseTapChangerTabular) && return _tapTablePointFactors(store, tc, step)
  return (1.0, _phaseTapXCurveFactor(store, tc, step), 1.0, 1.0)
end

# Impedance factors of one tabular tap changer at `step`: a
# RatioTapChangerTablePoint or PhaseTapChangerTablePoint gives r, x, g, b as
# deviations in percent of the nominal value, so the factors are (1 + r/100)
# on the series resistance, (1 + x/100) on the series reactance and the same
# for the magnetizing g and b. (1, 1, 1, 1) without a table point.
function _tapTablePointFactors(store::CGMESStore, tc::CIMObject, step::Float64)::NTuple{4,Float64}
  pt = _tapTablePointAtStep(store, tc, step)
  pt === nothing && return (1.0, 1.0, 1.0, 1.0)
  return (1.0 + num(pt, :r, 0.0) / 100.0, 1.0 + num(pt, :x, 0.0) / 100.0, 1.0 + num(pt, :g, 0.0) / 100.0, 1.0 + num(pt, :b, 0.0) / 100.0)
end

# Whether the impedance factors of `tc` apply: a ratio or tabular changer
# whose tap correction at its step lies in the plausibility band (outside it
# `_applyEndTaps` ignores the changer, and so do the impedance factors), or
# a linear or non-linear phase changer with a reactance curve.
function _tapTableChangerApplies(store::CGMESStore, svsteps::Dict{String,Float64}, tc::CIMObject)::Bool
  tc.class in (:RatioTapChanger, :PhaseTapChangerTabular) || _phaseTapXCurve(store, tc) !== nothing || return false
  corr = tc.class == :RatioTapChanger ? first(_ratioTapCorrectionAt(tc, _ratioTapTable(store, tc), _tapStep(tc, svsteps))) : first(_phaseTapRatioShift(store, tc, svsteps))
  return CGMES_TAP_CORR_MIN <= corr <= CGMES_TAP_CORR_MAX
end

# Impedance factors of the tap changers on end `e` at their step. Without
# tap_control the step is the solved SvTapStep, with tap_control the SSH
# start step; a controller that moves the step then sets the impedance of
# the new step through the branch's TapStepImpedance
# (`_attachTapStepImpedance!`). Returns `(r, x, g, b)`, all 1.0 without
# table points.
function _tapTableImpedanceFactors(store::CGMESStore, svsteps::Dict{String,Float64}, e::CIMObject)::NTuple{4,Float64}
  fr, fx, fg, fb = 1.0, 1.0, 1.0, 1.0
  for tc in _tapChangersOfEnd(store, e.mrid)
    _tapTableChangerApplies(store, svsteps, tc) || continue
    f = _tapChangerImpedanceFactors(store, tc, _tapStep(tc, svsteps))
    fr *= f[1]
    fx *= f[2]
    fg *= f[3]
    fb *= f[4]
  end
  return (fr, fx, fg, fb)
end

# --- Stage-2 tap control (importCGMES(..., tap_control = true)) -------------
#
# With `tap_control = true` the importer starts from the SSH tap positions
# and attaches Sparlectra outer-loop controllers (`addPowerTransformerControl!`)
# for every tap changer whose `controlEnabled` flag AND TapChangerControl
# `enabled` flag are set. The reference SvTapStep positions are then found by
# the control loop (run via `run_sparlectra`), not copied from SV.

# addPowerTransformerControl! returns the net, not the controller: the new
# controller is the one entry of _tap_controllers(net) that `before` lacks.
# Its position is resolved at the end of importCGMES (see _MapCtx).
function _recordNewTapController!(ctx::_MapCtx, net, before::Vector{Sparlectra.PowerTransformerControl}, tc_mrid::String)
  known = IdDict{Sparlectra.PowerTransformerControl,Nothing}(c => nothing for c in before)
  for c in Sparlectra._tap_controllers(net)
    haskey(known, c) || push!(ctx.tap_controller_by_mrid, tc_mrid => c)
  end
  return nothing
end

# Tap range of the ratio tap changer `tc` on a branch with the present ratio
# `tap_ratio`: ratio limits, step and, with a RatioTapChangerTable, the
# ratio of every table point (`points`). The range lies around the nominal
# ratio, the present ratio without the correction of the start step; at the
# branch from side the correction multiplies the ratio, at the to side it
# divides it. The linear formula keeps its uniform increment; a table hands
# a controller its ratio per step, and `tap_step` is then the mean step over
# the range (probe size, report).
function _ratioTapRange(store::CGMESStore, svsteps::Dict{String,Float64}, tc::CIMObject, tap_ratio::Float64, on_from_side::Bool)
  low = Int(round(num(tc, :lowStep, 0.0)))
  high = Int(round(num(tc, :highStep, 0.0)))
  table = _ratioTapTable(store, tc)
  corr0, _ = _ratioTapCorrectionAt(tc, table, _tapStep(tc, svsteps))
  rho_nom = on_from_side ? tap_ratio / corr0 : tap_ratio * corr0
  ratios = Float64[(c = first(_ratioTapCorrectionAt(tc, table, Float64(k))); on_from_side ? rho_nom * c : rho_nom / c) for k in min(low, high):max(low, high)]
  tmin, tmax = minimum(ratios), maximum(ratios)
  (table === nothing || isempty(table)) && return (tap_min = tmin, tap_max = tmax, tap_step = abs(rho_nom * num(tc, :stepVoltageIncrement, 0.0) / 100.0), points = nothing)
  return (tap_min = tmin, tap_max = tmax, tap_step = high != low ? (tmax - tmin) / abs(high - low) : 0.0, points = sort!(unique(ratios)))
end

function _setRatioTapRange!(br, rng)
  br.has_ratio_tap = true
  br.tap_min = rng.tap_min
  br.tap_max = rng.tap_max
  br.tap_step = rng.tap_step
  # the step ratios of a CGMES changer are linear in the step (or table
  # points), not the fraction grid of nameplate data; the tap controller
  # steps on this grid
  br.tap_ratio_grid = :linear
  # a linear changer steps on the uniform grid, also where another changer
  # of the branch gave it table points before
  br.tap_ratio_points = rng.points
  return br
end

# Phase range of the phase tap changer `tc` on a branch with the present
# shift `shift_deg`: shift limits and mean step around the base shift, the
# present shift without the contribution of the start step. An end at the
# branch to side enters negated, as in _applyEndTaps.
function _phaseTapRange(store::CGMESStore, svsteps::Dict{String,Float64}, tc::CIMObject, shift_deg::Float64, on_from_side::Bool)
  low = Int(round(num(tc, :lowStep, 0.0)))
  high = Int(round(num(tc, :highStep, 0.0)))
  _, ps0 = _phaseTapRatioShift(store, tc, svsteps)
  base_shift = on_from_side ? shift_deg - ps0 : shift_deg + ps0
  _, ps_low = _phaseTapRatioShift_atstep(store, tc, low)
  _, ps_high = _phaseTapRatioShift_atstep(store, tc, high)
  s1 = on_from_side ? base_shift + ps_low : base_shift - ps_low
  s2 = on_from_side ? base_shift + ps_high : base_shift - ps_high
  return (phase_min_deg = min(s1, s2), phase_max_deg = max(s1, s2), phase_step_deg = high > low ? abs(s2 - s1) / (high - low) : 1.0)
end

function _setPhaseTapRange!(br, rng)
  br.has_phase_tap = true
  br.phase_min_deg = rng.phase_min_deg
  br.phase_max_deg = rng.phase_max_deg
  br.phase_step_deg = rng.phase_step_deg
  return br
end

# (ratio, shift in degrees) of tap changer `tc` at its import step, the
# correction _applyEndTaps folds into the branch (a RatioTapChanger has no
# shift); ratio 1 and shift 0 for a PhaseTapChangerTabular whose table does
# not resolve.
function _tapCorrectionAt(store::CGMESStore, svsteps::Dict{String,Float64}, tc::CIMObject)::Tuple{Float64,Float64}
  tc.class == :RatioTapChanger && return (first(_ratioTapCorrectionAt(tc, _ratioTapTable(store, tc), _tapStep(tc, svsteps))), 0.0)
  return _phaseTapRatioShift(store, tc, svsteps)
end

# Whether _applyEndTaps folds the correction of `tc` into the branch: not
# for a correction outside the plausibility band, nor for a
# PhaseTapChangerTabular whose table does not resolve.
function _tapCorrectionApplied(store::CGMESStore, svsteps::Dict{String,Float64}, tc::CIMObject)::Bool
  tc.class == :PhaseTapChangerTabular && _phaseTapTableModel(store, tc, Int(round(_tapStep(tc, svsteps)))) === nothing && return false
  corr = first(_tapCorrectionAt(store, svsteps, tc))
  return CGMES_TAP_CORR_MIN <= corr <= CGMES_TAP_CORR_MAX
end

# Tap range, table ratios and impedance per step of the tap changers on the
# ends `ends` (end, whether it sits at the branch from side) of branch
# `branch_idx`, regulated or not, so that a controller configured later
# (control section) steps on the delivered range and table. Without this a
# configured controller stepped on the constructor grid (0.9 to 1.1 in
# steps of 0.00625). The changers the import regulates itself
# (`controlled`, from _attachTapControl!) keep the range of their
# controller; of the others the first ratio and the first phase changer of
# the branch give the range, since the branch has one ratio and one shift.
# Left out: a changer whose correction the import ignores
# (_tapCorrectionApplied) and a changer without a range. Ratio, shift and
# impedance of the start step stay as they are, so the solve does not
# change.
function _attachTapRanges!(net, store::CGMESStore, svsteps::Dict{String,Float64}, ctx::_MapCtx, branch_idx::Int, ends::Vector{Tuple{CIMObject,Bool}}, controlled::Vector{Tuple{CIMObject,Bool}}, trafo_name::AbstractString)
  br_pos = findfirst(b -> b.branchIdx == branch_idx, net.branchVec)
  br_pos === nothing && return nothing
  br = net.branchVec[br_pos]
  ratio_done = any(c -> c[1].class == :RatioTapChanger, controlled)
  phase_done = any(c -> c[1].class != :RatioTapChanger, controlled)
  ranged = Tuple{CIMObject,Bool}[]
  for (e, on_from_side) in ends, tc in _tapChangersOfEnd(store, e.mrid)
    any(c -> c[1] === tc, controlled) && continue
    ratio = tc.class == :RatioTapChanger
    (ratio ? ratio_done : phase_done) && continue
    _tapCorrectionApplied(store, svsteps, tc) || continue
    if ratio
      rng = _ratioTapRange(store, svsteps, tc, br.tap_ratio, on_from_side)
      rng.tap_max > rng.tap_min || continue
      _setRatioTapRange!(br, rng)
      ratio_done = true
    else
      rng = _phaseTapRange(store, svsteps, tc, br.phase_shift_deg, on_from_side)
      rng.phase_max_deg > rng.phase_min_deg || continue
      _setPhaseTapRange!(br, rng)
      phase_done = true
    end
    push!(ranged, (tc, on_from_side))
  end
  _attachTapStepImpedance!(net, store, svsteps, ctx, branch_idx, vcat(controlled, ranged), trafo_name)
  return nothing
end

# A controller moves the ratio or the angle of its branch; the r/x/g/b
# deviations of a tap table belong to the step. For every tabular changer
# with a range on the branch (`changers`: changer and whether it sits at
# the branch from side; the regulated ones first) the branch gets a
# TapStepImpedance: the branch values without the factors of these
# changers at their start step, and per step the value the actuator takes
# there with the factors of that step. A tap controller applies it after
# each move. A second table on the same actuator cannot be told apart by
# the branch value; its impedance stays at the start step, with a notice.
function _attachTapStepImpedance!(net, store::CGMESStore, svsteps::Dict{String,Float64}, ctx::_MapCtx, branch_idx::Int, changers::Vector{Tuple{CIMObject,Bool}}, trafo_name::AbstractString)
  isempty(changers) && return nothing
  br_pos = findfirst(b -> b.branchIdx == branch_idx, net.branchVec)
  br_pos === nothing && return nothing
  br = net.branchVec[br_pos]
  ratio_values, ratio_factors = Float64[], NTuple{4,Float64}[]
  phase_values, phase_factors = Float64[], NTuple{4,Float64}[]
  f0 = (1.0, 1.0, 1.0, 1.0)
  for (tc, on_from_side) in changers
    _tapTableChangerApplies(store, svsteps, tc) || continue
    low = Int(round(num(tc, :lowStep, 0.0)))
    high = Int(round(num(tc, :highStep, 0.0)))
    steps = low:high
    facs = [_tapChangerImpedanceFactors(store, tc, Float64(s)) for s in steps]
    all(f -> f == (1.0, 1.0, 1.0, 1.0), facs) && continue
    ratio = tc.class == :RatioTapChanger
    name = something(str(tc, :name), tc.mrid)
    if !isempty(ratio ? ratio_values : phase_values)
      push!(ctx.messages, "notice: tap changer $(trafo_name)/$(name), a second table on the same branch; its r/x/g/b deviations stay at the start step")
      continue
    end
    start = _tapChangerImpedanceFactors(store, tc, _tapStep(tc, svsteps))
    any(iszero, start) && (push!(ctx.messages, "notice: tap changer $(trafo_name)/$(name), the table point of the start step removes the impedance; its r/x/g/b deviations stay at the start step"); continue)
    f0 = f0 .* start
    if ratio
      # the branch tap_ratio holds the end correction at the from side as a
      # factor, at the to side as a divisor (see _attachTapControl!). The
      # corrections come from the RatioTapChangerTable where there is one,
      # like the ratios the controller steps to (tap_ratio_points)
      table = _ratioTapTable(store, tc)
      corr0 = first(_ratioTapCorrectionAt(tc, table, _tapStep(tc, svsteps)))
      vals = [(c = first(_ratioTapCorrectionAt(tc, table, Float64(s))); on_from_side ? br.tap_ratio / corr0 * c : br.tap_ratio * corr0 / c) for s in steps]
      append!(ratio_values, vals)
      append!(ratio_factors, facs)
    else
      ps0 = last(_phaseTapRatioShift(store, tc, svsteps))
      vals = [(ps = last(_phaseTapRatioShift_atstep(store, tc, s)); on_from_side ? br.phase_shift_deg - ps0 + ps : br.phase_shift_deg + ps0 - ps) for s in steps]
      append!(phase_values, vals)
      append!(phase_factors, facs)
    end
  end
  (isempty(ratio_values) && isempty(phase_values)) && return nothing
  br.tap_step_impedance = Sparlectra.TapStepImpedance(br.r_pu / f0[1], br.x_pu / f0[2], br.g_from_pu / f0[3], br.b_from_pu / f0[4], br.g_to_pu / f0[3], br.b_to_pu / f0[4], ratio_values, ratio_factors, phase_values, phase_factors)
  return nothing
end

# --- tap changer sets of the tap estimation ----------------------------------
#
# The import stores on every transformer branch the tap changers the state
# estimation releases (`se_tap_changers`, a Sparlectra.SETapChangerSet with
# origin :cgmes). Power flow, controllers and export do not read it; the
# branch ratio, shift, ranges and impedance stay as the mapping above sets
# them.

# The single-step PhaseTapChangerLinear that the own CGMES export writes to
# carry a fixed branch shift (cgmes_export.jl, _writeLinearPtc): lowStep 0,
# highStep 1, neutralStep 0, normalStep 1, step 1, and a name ending in
# "_PTC". All five values and the suffix must match. Such a changer is the
# shift of the branch, not a regulator: a re-import of an exported nameplate
# Delta-u changer would otherwise release a one-step angle changer.
function _isExportShiftPtc(tc::CIMObject)::Bool
  tc.class == :PhaseTapChangerLinear || return false
  endswith(something(str(tc, :name), ""), "_PTC") || return false
  return (num(tc, :lowStep), num(tc, :highStep), num(tc, :neutralStep), num(tc, :normalStep), num(tc, :step)) == (0.0, 1.0, 0.0, 1.0, 1.0)
end

# Points of the RatioTapChangerTable or PhaseTapChangerTable of `tc`, sorted
# by step: the ratio (and angle) of the point and its r, x, g, b percent
# deviations as factors 1 + p/100. A ratio point needs a ratio (as
# _ratioTapTable reads it), a phase point takes ratio 1 and angle 0 where it
# leaves them out (as _phaseTapTableModel reads it). `nothing` without a
# usable point, or when a step is not a whole number or two points share a
# step (the estimation needs whole, strictly ascending steps).
function _seTapTablePoints(store::CGMESStore, tc::CIMObject)::Union{Nothing,Vector{Sparlectra.SETapPoint}}
  ratio_tc = tc.class == :RatioTapChanger
  tabkey, pointcls = ratio_tc ? (:RatioTapChangerTable, :RatioTapChangerTablePoint) : (:PhaseTapChangerTable, :PhaseTapChangerTablePoint)
  tmrid = get(tc.refs, tabkey, nothing)
  tmrid === nothing && return nothing
  pts = Sparlectra.SETapPoint[]
  for pt in objectsOf(store, pointcls)
    get(pt.refs, tabkey, "") == tmrid || continue
    s = num(pt, :step)
    r = ratio_tc ? num(pt, :ratio) : num(pt, :ratio, 1.0)
    (s === nothing || r === nothing) && continue
    isinteger(s) || return nothing
    push!(pts, Sparlectra.SETapPoint(step = Int(s), ratio = r, angle_deg = ratio_tc ? 0.0 : num(pt, :angle, 0.0), r = 1.0 + num(pt, :r, 0.0) / 100.0, x = 1.0 + num(pt, :x, 0.0) / 100.0, g = 1.0 + num(pt, :g, 0.0) / 100.0, b = 1.0 + num(pt, :b, 0.0) / 100.0))
  end
  isempty(pts) && return nothing
  sort!(pts; by = p -> p.step)
  all(k -> pts[k].step < pts[k + 1].step, 1:length(pts) - 1) || return nothing
  all(p -> p.ratio > 0.0, pts) || return nothing
  return pts
end

_seTapTableHasImpedance(pts::Vector{Sparlectra.SETapPoint})::Bool = any(p -> (p.r, p.x, p.g, p.b) != (1.0, 1.0, 1.0, 1.0), pts)

# Tap changer `tc` (at the branch from side when `on_from_side`) for the tap
# estimation, or the reason why it cannot be one. The numbers are those the
# import computes the branch ratio and shift from:
# - RatioTapChanger: _ratioTapCorrectionAt with the unrounded step; with a
#   RatioTapChangerTable kind :tabular, continued outside the table by the
#   linear formula (outside_du = stepVoltageIncrement);
# - asymmetrical and symmetrical: _phaseTapRatioShift (rounded step, band
#   defaults around step and neutral);
# - PhaseTapChangerLinear: _phaseTapLinearModel (band widened to neutral,
#   step clamped into it);
# - PhaseTapChangerTabular: _phaseTapTableModel (its neutral and the step
#   snapped onto a table row).
# Impedance data only where the import applies it: table points with a
# factor other than 1, and the accepted xMin/xMax curve (_phaseTapXCurve).
function _seTapChanger(store::CGMESStore, svsteps::Dict{String,Float64}, tc::CIMObject, on_from_side::Bool)::Union{Sparlectra.SETapChanger,String}
  name = something(str(tc, :name), tc.mrid)
  side = on_from_side ? :from : :to
  assumed = !_tapStepGiven(tc, svsteps)
  if tc.class == :RatioTapChanger
    neutral = Int(round(num(tc, :neutralStep, 0.0)))
    lo, hi = minmax(Int(round(num(tc, :lowStep, 0.0))), Int(round(num(tc, :highStep, 0.0))))
    du = num(tc, :stepVoltageIncrement, 0.0) / 100.0
    step = _tapStep(tc, svsteps)
    table = _ratioTapTable(store, tc)
    (table === nothing || isempty(table)) && return Sparlectra.SETapChanger(name = name, mrid = tc.mrid, step_assumed = assumed, kind = :ratio, role = :ratio, side = side, du = du, neutral = neutral, low = lo, high = hi, step = step)
    pts = _seTapTablePoints(store, tc)
    pts === nothing && return "its RatioTapChangerTable has steps that are not whole and strictly ascending"
    return Sparlectra.SETapChanger(name = name, mrid = tc.mrid, step_assumed = assumed, kind = :tabular, role = :ratio, side = side, neutral = neutral, low = lo, high = hi, step = step, table = pts, outside_du = du, impedance = _seTapTableHasImpedance(pts) ? :table : :none)
  end
  step = Int(round(_tapStep(tc, svsteps)))
  neutral = Int(round(num(tc, :neutralStep, 0.0)))
  lo = Int(round(num(tc, :lowStep, Float64(min(step, neutral)))))
  hi = Int(round(num(tc, :highStep, Float64(max(step, neutral)))))
  if tc.class == :PhaseTapChangerTabular
    # resolves: _attachSETapChangers! passes only changers the import
    # applied (_tapCorrectionApplied), and an unresolved table is not one
    built = something(_phaseTapTableModel(store, tc, step))
    pts = _seTapTablePoints(store, tc)
    pts === nothing && return "its PhaseTapChangerTable has steps that are not whole and strictly ascending"
    m = built.model
    return Sparlectra.SETapChanger(name = name, mrid = tc.mrid, step_assumed = assumed, kind = :tabular, role = :phase, side = side, neutral = m.neutralStep, low = pts[1].step, high = pts[end].step, step = m.step, table = pts, impedance = _seTapTableHasImpedance(pts) ? :table : :none)
  end
  xc = _phaseTapXCurve(store, tc)
  xcurve = xc === nothing ? nothing : (xc.xmin, xc.xmax, xc.amax)
  ximp = xc === nothing ? :none : :x_curve
  if tc.class == :PhaseTapChangerLinear
    lo <= neutral <= hi || (lo = min(lo, neutral); hi = max(hi, neutral))
    return Sparlectra.SETapChanger(name = name, mrid = tc.mrid, step_assumed = assumed, kind = :linear, role = :phase, side = side, dphi_deg = num(tc, :stepPhaseShiftIncrement, 0.0), neutral = neutral, low = lo, high = hi, step = clamp(step, lo, hi), x_curve = xcurve, x_curve_form = :sin2, impedance = ximp)
  end
  lo, hi = minmax(lo, hi)
  # both values are present: the import computed the correction of this
  # changer from them (_phaseTapRatioShift), which throws without them
  # (calcPhaseTapFraction, PhaseTapChangerModel for :asymmetrical), so a
  # delivery without them never reaches this point
  vinc = something(num(tc, :voltageStepIncrement))
  tc.class == :PhaseTapChangerSymmetrical && return Sparlectra.SETapChanger(name = name, mrid = tc.mrid, step_assumed = assumed, kind = :symmetrical, role = :phase, side = side, du = vinc / 100.0, neutral = neutral, low = lo, high = hi, step = step, x_curve = xcurve, x_curve_form = :sin2, impedance = ximp)
  alpha = something(num(tc, :windingConnectionAngle))
  return Sparlectra.SETapChanger(name = name, mrid = tc.mrid, step_assumed = assumed, kind = :asymmetrical, role = :phase, side = side, du = vinc / 100.0, alpha_deg = alpha, neutral = neutral, low = lo, high = hi, step = step, x_curve = xcurve, x_curve_form = :sin2, impedance = ximp)
end

# Store the tap changer set of transformer branch `branch_idx`: the changers
# of `ends` (end, whether it sits at the branch from side) in
# _tapChangersOfEnd order, on the branch ratio `ratio_n` before the end taps
# (the structural ratio _applyEndTaps starts from). Left out, each with a
# notice: a changer the import ignored (_tapCorrectionApplied), a changer
# whose table the estimation cannot use, and the single-step phase changer
# of the own export (_isExportShiftPtc). The import applied the last two,
# so their ratio and shift at the import step go into ratio_n and
# shift_n_deg as a fixed part: the set then gives the branch ratio and
# shift at the import steps, and the released changers do not absorb the
# left-out factor. A transformer without changers gets an empty set and is
# never released.
function _attachSETapChangers!(net, store::CGMESStore, svsteps::Dict{String,Float64}, ctx::_MapCtx, branch_idx::Int, ends::Vector{Tuple{CIMObject,Bool}}, ratio_n::Float64, trafo_name::AbstractString)
  br_pos = findfirst(b -> b.branchIdx == branch_idx, net.branchVec)
  br_pos === nothing && return nothing
  chs = Sparlectra.SETapChanger[]
  shift_n = 0.0
  for (e, on_from_side) in ends, tc in _tapChangersOfEnd(store, e.mrid)
    name = something(str(tc, :name), tc.mrid)
    if !_tapCorrectionApplied(store, svsteps, tc)
      push!(ctx.messages, "notice: tap changer $(trafo_name)/$(name) is not estimable, the import ignored it")
      continue
    end
    shift_ptc = _isExportShiftPtc(tc)
    ch = shift_ptc ? nothing : _seTapChanger(store, svsteps, tc, on_from_side)
    if ch isa Sparlectra.SETapChanger
      push!(chs, ch)
      continue
    end
    # fixed part, with the sides of _applyEndTaps: from multiplies and adds
    pr, ps = _tapCorrectionAt(store, svsteps, tc)
    ratio_n = on_from_side ? ratio_n * pr : ratio_n / pr
    shift_n += on_from_side ? ps : -ps
    push!(ctx.messages, shift_ptc ? "notice: tap changer $(trafo_name)/$(name) is the single-step phase changer of a CGMES export: its shift of $(round(ps; digits = 6)) deg stays fixed, the tap estimation does not release it" : "notice: tap changer $(trafo_name)/$(name) is not estimable, $(ch): its ratio $(round(pr; digits = 6)) and shift $(round(ps; digits = 6)) deg at step $(_tapStep(tc, svsteps)) stay fixed")
  end
  net.branchVec[br_pos].se_tap_changers = Sparlectra.SETapChangerSet(ratio_n = ratio_n, shift_n_deg = shift_n, changers = chs, origin = :cgmes)
  return nothing
end

# A transformer-type branch the import builds for a line or an
# EquivalentBranch across voltage levels, or for a phase-stamped
# EquivalentBranch, has no tap changer: an empty set keeps it out of the tap
# estimation (the constructor marks every transformer branch with a default
# ratio grid, which is no changer model).
function _attachNoSETapChangers!(br)
  br.se_tap_changers = Sparlectra.SETapChangerSet(ratio_n = br.ratio, shift_n_deg = br.angle, changers = Sparlectra.SETapChanger[], origin = :cgmes)
  return br
end

"""
Branch tap-range fields and controller for one controlled transformer end.
Returns the tap changers that got a controller or joined a controller
group, each with `on_from_side`.
"""
function _attachTapControl!(net, store::CGMESStore, topo::CGMESTopology, ctx::_MapCtx, svsteps::Dict{String,Float64}, e::CIMObject, branch_idx::Int, on_from_side::Bool, trafo_name::String)::Vector{Tuple{CIMObject,Bool}}
  attached = Tuple{CIMObject,Bool}[]
  # loop-invariant: branch_idx is a function argument, resolve it once
  br_pos = findfirst(b -> b.branchIdx == branch_idx, net.branchVec)
  br_pos === nothing && return attached
  for tc in _tapChangersOfEnd(store, e.mrid)
    something(boolval(tc, :controlEnabled), false) || continue
    tcc = ref(store, tc, :TapChangerControl)
    tcc === nothing && continue
    something(boolval(tcc, :enabled), false) || continue
    mode = something(enumval(tcc, :mode), "")
    br = net.branchVec[br_pos]
    name = something(str(tc, :name), tc.mrid)

    if tc.class == :RatioTapChanger && endswith(mode, ".voltage")
      # step range as branch tap-ratio limits around the nominal ratio, from
      # the same corrections as the fixed-tap import (table or linear)
      _setRatioTapRange!(br, _ratioTapRange(store, svsteps, tc, br.tap_ratio, on_from_side))
      # master/slave group dedup (#322): several tap changers referencing
      # ONE TapChangerControl are a regulated parallel-transformer group.
      # The first one becomes the master controller; every later one joins
      # as a FOLLOWER (its branch keeps the tap machinery set above so the
      # master can mirror its moves) instead of spawning a second
      # independent controller on the same target bus, which would fight
      # the first through circulating reactive power.
      if haskey(ctx.tap_group_masters, tcc.mrid)
        master = ctx.tap_group_masters[tcc.mrid]
        push!(master.followers, string(branch_idx))
        push!(ctx.messages, "tap control: $(trafo_name)/$(name) follows the group of $(master.trafo) (shared TapChangerControl)")
        push!(attached, (tc, on_from_side))
        continue
      end
      # regulated bus from the TapChangerControl terminal
      tinfo = _regulatedBus(store, topo, tcc)
      tinfo === nothing && (push!(ctx.messages, "notice: tap controller $(name), control terminal unresolved, skipped"); continue)
      target = num(tcc, :targetValue, 0.0)
      target <= 0.0 && (push!(ctx.messages, "notice: tap controller $(name), no positive voltage target, skipped"); continue)
      deadband = num(tcc, :targetDeadband, 0.0)
      before = Sparlectra._tap_controllers(net)
      Sparlectra.addPowerTransformerControl!(
        net;
        trafo = string(branch_idx),
        mode = :voltage,
        target_bus = tinfo.bus,
        target_vm_pu = target / tinfo.vn_kV,
        deadband_vm_pu = max(deadband / tinfo.vn_kV / 2.0, 1e-3),
        control_ratio = true,
        control_phase = false,
      )
      _recordNewTapController!(ctx, net, before, tc.mrid)
      push!(attached, (tc, on_from_side))
      # remember the master by its unambiguous trafo id (one voltage
      # controller per imported transformer)
      for c in Sparlectra._tap_controllers(net)
        c.trafo == string(branch_idx) && (ctx.tap_group_masters[tcc.mrid] = c)
      end
      push!(ctx.messages, "tap control: $(trafo_name)/$(name) → voltage $(round(target / tinfo.vn_kV; digits = 4)) pu at $(tinfo.bus)")
    elseif tc.class != :RatioTapChanger && endswith(mode, ".activePower")
      # phase tap: step range as branch phase limits around the base shift
      _setPhaseTapRange!(br, _phaseTapRange(store, svsteps, tc, br.phase_shift_deg, on_from_side))
      # monitored branch side from the control terminal
      tccterm = get(tcc.refs, :Terminal, nothing)
      side = tccterm === nothing ? nothing : get(ctx.branch_side, tccterm, nothing)
      from_name = _busDictName(net, br.fromBus)
      to_name = _busDictName(net, br.toBus)
      (from_name === nothing || to_name === nothing) && (push!(ctx.messages, "notice: tap controller $(name), bus name lookup failed, skipped"); continue)
      # the controller monitors the branch in stored (from → to) orientation;
      # a control terminal on the to side flips the target sign (the loss
      # asymmetry between the two ends is far below the deadband)
      target_branch = (from_name, to_name)
      target = num(tcc, :targetValue, 0.0)
      (side !== nothing && side[2] == :to) && (target = -target)
      deadband = max(num(tcc, :targetDeadband, 1.0), 0.5)
      before = Sparlectra._tap_controllers(net)
      Sparlectra.addPowerTransformerControl!(
        net;
        trafo = string(branch_idx),
        mode = :branch_active_power,
        target_branch = target_branch,
        p_target_mw = target,
        deadband_p_mw = deadband / 2.0,
        control_ratio = false,
        control_phase = true,
      )
      _recordNewTapController!(ctx, net, before, tc.mrid)
      push!(attached, (tc, on_from_side))
      push!(ctx.messages, "tap control: $(trafo_name)/$(name) → active power $(target) MW on $(target_branch[1]) → $(target_branch[2])")
    else
      push!(ctx.messages, "notice: tap controller $(name), unsupported combination $(tc.class)/$(mode), skipped")
    end
  end
  return attached
end

"""Phase-tap (ratio, shift) of `ptc` evaluated at an explicit step."""
function _phaseTapRatioShift_atstep(store::CGMESStore, ptc::CIMObject, step::Int)::Tuple{Float64,Float64}
  fake = Dict{String,Float64}(ptc.mrid => Float64(step))
  return _phaseTapRatioShift(store, ptc, fake)
end

# Whether a prosumer energizes its component. A source is a generator-type
# unit (`isGenerator`, i.e. `proSumptionType == Injection`) that either
# injects active power (`pVal > 0`, generator sign) or holds a voltage
# (`isRegulating`: `isRegulated`, or one of `vset_adjust`, `vstep_pu`,
# `tap_steps_down`, `tap_steps_up`). The slack unit is covered by the second
# test, the importer marks it `isRegulated`. Loads never count, and neither
# does a unit without regulation that has P = 0 or draws power (`pVal < 0`):
# real deliveries carry boundary fragments whose only units are
# EquivalentInjections with P = 0 or drawing power, and counting them kept
# those fragments energized without anything that could feed their load or
# hold their voltage. A machine whose regulation targets a remote bus
# (`machine_control = true`) is not `isRegulated` at its own bus;
# `_deEnergizeSourcelessComponents!` counts it through its control plan.
_isCGMESSource(ps)::Bool = Sparlectra.isGenerator(ps) && (something(ps.pVal, 0.0) > 0.0 || Sparlectra.isRegulating(ps))

# Suffix of the de-energizing notices for the P that loads fed into a
# de-energized component (load sign P < 0), empty when there is none.
_droppedLoadInfeed(p::Float64)::String = p > 0.0 ? ", load infeed P = $(round(p; digits = 2)) MW dropped" : ""

# Switch off the buses `buses`: injections zeroed, every branch touching
# them out of service, the buses typed Isolated. Returns the number of
# zeroed injections, the unsupplied load P, the dropped load infeed P and
# the dropped generation P of units that fed in (MW). The caller reports;
# markIsolatedBuses! sets the voltages to 0.
function _deEnergizeBuses!(net, buses::Vector{Int})::Tuple{Int,Float64,Float64,Float64}
  busset = Set(buses)
  # Zero the injections rather than deleting the prosumers: the components
  # stay visible in reports, they just no longer contribute to the balance.
  nload = 0
  p_load = 0.0
  p_fed = 0.0
  p_gen = 0.0
  for ps in net.prosumpsVec
    ps.comp.cFrom_bus in busset || continue
    # unsupplied: what a load consumes, and what a generator-type unit drew
    # (pVal <= 0). A load that feeds in (pVal < 0) is kept apart: netting it
    # against the drawn power would report a smaller or even negative
    # unsupplied load. A unit that fed in (pVal > 0, only in an island
    # switched off for its unit types) is dropped generation.
    p = something(ps.pVal, 0.0)
    if Sparlectra.isGenerator(ps)
      p > 0.0 ? (p_gen += p) : (p_load -= p)
    else
      p >= 0.0 ? (p_load += p) : (p_fed -= p)
    end
    ps.pVal = 0.0
    ps.qVal = 0.0
    nload += 1
  end
  # Any branch touching the buses goes out of service. Testing only one end
  # is enough -- a branch with one end inside cannot have the other end in
  # an energized component, or the component would have been connected.
  for br in net.branchVec
    (br.fromBus in busset || br.toBus in busset) && (br.status = 0)
  end
  for i in buses
    i in net.isoNodes && continue
    Sparlectra.setNodeType!(net.nodeVec[i], "Isolated")
    push!(net.isoNodes, i)
  end
  return nload, p_load, p_fed, p_gen
end

# Real deliveries contain connected component(s) without any source: every
# generator/injection in them is switched off, or the feeding branches are
# open (RealGrid: 13 pockets totalling 22 buses). Solving those is impossible
# (they have no reference), so they are de-energized deliberately: loads
# zeroed, branches out of service, buses isolated, one message per component.
# Without this the solver aborts on "island without reference".

"""
De-energize every connected component that contains no source.

A CGMES snapshot regularly carries network parts that are switched off in
reality: open bays, stub connections, feeders parked out of service. They still
appear as buses with loads attached. Left alone, such a component has load but
no generation and no reference, an unsolvable subsystem that would fail the
power flow for a reason that has nothing to do with the interesting network.

The pass therefore zeroes the injections, opens the branches and marks the
buses `Isolated` (`markIsolatedBuses!` later sets their voltage to 0). This
changes the imported model, so every affected component is reported, and one
notice counts the components, their buses and the load left unsupplied.

A source is a generator-type unit with P > 0 (generator sign) or with voltage
regulation (`_isCGMESSource`). Voltage regulation includes a machine that
regulates a remote bus of the same component with `machine_control = true`.
A unit without regulation that has P = 0 or draws power (P < 0) does not keep
a component energized, so a component whose units only draw power is
de-energized as well. A load is never a source, also not one that feeds power
in (P < 0 in load sign). The unsupplied P of a component is the P its loads
consume plus the P its drawing units drew; the P that feeding loads put in is
reported as dropped load infeed and not netted against it.

Connectivity is evaluated over in-service branches AND closed links: the
electrical view, since a component fed through a closed switch is energized.
"""
function _deEnergizeSourcelessComponents!(net, ctx::_MapCtx)
  n = length(net.nodeVec)
  n == 0 && return nothing
  parent = collect(1:n)
  find(i) = (while parent[i] != i
    parent[i] = parent[parent[i]]
    i = parent[i]
  end; i)
  function unite(a, b)
    ra, rb = find(a), find(b)
    ra == rb && return nothing
    parent[max(ra, rb)] = min(ra, rb)
    return nothing
  end
  for br in net.branchVec
    br.status == 1 && unite(br.fromBus, br.toBus)
  end
  for l in net.linkVec
    l.status == 1 && unite(Int(l.fromBus), Int(l.toBus))
  end

  # Group buses by component root and record whether any source sits inside
  # (_isCGMESSource: P > 0 or voltage regulation; loads and drawing units do
  # not energize anything).
  has_source = Dict{Int,Bool}()
  members = Dict{Int,Vector{Int}}()
  for i = 1:n
    r = find(i)
    push!(get!(Vector{Int}, members, r), i)
    haskey(has_source, r) || (has_source[r] = false)
  end
  for ps in net.prosumpsVec
    _isCGMESSource(ps) || continue
    has_source[find(ps.comp.cFrom_bus)] = true
  end
  # A machine that regulates a remote bus (machine_control = true) is not
  # isRegulated at its own bus, so _isCGMESSource sees only its P. Its
  # voltage regulation counts all the same, when the target bus lies in the
  # machine's own component; a target elsewhere holds nothing here.
  for plan in ctx.machine_ctrl_plans
    t = get(net.busDict, plan.target_bus, nothing)
    o = Int(net.prosumpsVec[plan.prosumer_idx].comp.cFrom_bus)
    t !== nothing && find(t) == find(o) && (has_source[find(o)] = true)
  end

  n_comp = 0
  n_bus = 0
  p_unsupplied = 0.0
  p_dropped = 0.0
  # sorted by root (the smallest bus index of the component), so the
  # messages come in the same order on every run
  for (r, buses) in sort!(collect(members); by = first)
    has_source[r] && continue
    # A single bus with nothing attached is a harmless leftover (open bay,
    # unconnected boundary node) -- skip it silently instead of drowning the
    # message log. Anything larger, or a lone bus that does carry an injection,
    # is reported.
    length(buses) == 1 && isempty(filter(ps -> ps.comp.cFrom_bus == buses[1], net.prosumpsVec)) && continue
    nload, p_load, p_fed, _ = _deEnergizeBuses!(net, buses)
    push!(ctx.messages, "notice: de-energized a component of $(length(buses)) bus(es) around $(net.nodeVec[first(buses)].comp.cName): no source inside (no unit with P > 0 or voltage regulation); $(nload) injection(s) zeroed, load P = $(round(p_load; digits = 2)) MW unsupplied$(_droppedLoadInfeed(p_fed))")
    n_comp += 1
    n_bus += length(buses)
    p_unsupplied += p_load
    p_dropped += p_fed
  end
  n_comp > 0 && push!(ctx.messages, "notice: $(n_comp) component(s) without source de-energized: $(n_bus) bus(es), load P = $(round(p_unsupplied; digits = 2)) MW unsupplied$(_droppedLoadInfeed(p_dropped))")
  isempty(net.isoNodes) || sort!(unique!(net.isoNodes))
  return nothing
end

# A voltage controller whose target bus is held by a generator (slack or PV)
# cannot influence that voltage: the machine fixes it. Such controllers occur
# in real data sets (MicroGrid BE-TR2_3 regulates the slack busbar); leaving
# them active would burn outer iterations and report a "regulation" that is
# none. Bus types are known only after refreshBusTypesFromProsumers!, so this
# runs as a post-mapping pass.
function _disableIneffectiveTapControllers!(net, ctx::_MapCtx)
  for ctrl in Sparlectra._tap_controllers(net)
    ctrl.enabled || continue
    ctrl.mode in (:voltage, :voltage_and_branch_active_power) || continue
    ctrl.target_bus === nothing && continue
    idx = get(net.busDict, ctrl.target_bus, nothing)
    idx === nothing && continue
    ntype = net.nodeVec[idx]._nodeType
    (ntype == Sparlectra.Slack || ntype == Sparlectra.PV) || continue
    ctrl.enabled = false
    push!(ctx.messages, "notice: tap controller on $(ctrl.target_bus) disabled, target bus is $(ntype) (voltage held by a generator)")
  end
  return nothing
end

"""Bus-dictionary name of node index `idx` (the name `geNetBusIdx` accepts)."""
function _busDictName(net, idx::Int)::Union{Nothing,String}
  for (k, v) in net.busDict
    v == idx && return k
  end
  return nothing
end

"""Regulated bus (name, vn) of a TapChangerControl, or `nothing`."""
function _regulatedBus(store::CGMESStore, topo::CGMESTopology, tcc::CIMObject)
  tmrid = get(tcc.refs, :Terminal, nothing)
  tmrid === nothing && return nothing
  t = get(store.objects, tmrid, nothing)
  t === nothing && return nothing
  tn = tnOfTerminal(t)
  (tn === nothing || !haskey(topo.bus_name, tn)) && return nothing
  return (bus = topo.bus_name[tn], vn_kV = topo.vn_kV[tn])
end

# tap changers attached to a transformer end (RatioTapChanger.TransformerEnd)
function _tapChangersOfEnd(store::CGMESStore, endmrid::String)
  out = CIMObject[]
  for cls in (:RatioTapChanger, :PhaseTapChangerLinear, :PhaseTapChangerSymmetrical, :PhaseTapChangerAsymmetrical, :PhaseTapChangerTabular)
    for tc in objectsOf(store, cls)
      get(tc.refs, :TransformerEnd, "") == endmrid && push!(out, tc)
    end
  end
  return out
end

# fold the fixed tap corrections of transformer end `e` into (ratio, shift).
# `on_from_side`: the end sits on the branch's from bus → ratio multiplies;
# Plausibility band for a single tap correction factor, same philosophy as the
# vset band: real tap ranges stay within a few ten percent of neutral (RealGrid
# tabular rows: 0.9 .. 1.1), while FullGrid's tabular table carries the literal
# placeholder row ratio 9.99 / angle 0.99° (the set's 0.99/99.99 family), a
# 10:1 off-nominal that alone produces a ~400 pu mismatch. Outside the band the
# tap is ignored with a warning; the transformer keeps its nominal ratio.
const CGMES_TAP_CORR_MIN = 0.5
const CGMES_TAP_CORR_MAX = 2.0

function _tapCorrPlausible(kind::String, tc::CIMObject, corr::Float64, messages::Vector{String}; strict::Bool = false)::Bool
  CGMES_TAP_CORR_MIN <= corr <= CGMES_TAP_CORR_MAX && return true
  msg = "$(kind) $(something(str(tc, :name), tc.mrid)) declares an implausible tap correction ($(round(corr; digits = 3)), band $(CGMES_TAP_CORR_MIN)..$(CGMES_TAP_CORR_MAX)), placeholder suspected"
  strict && error("cgmes_import.placeholder_guards = strict: " * msg)
  push!(messages, "warning: " * msg * ", tap ignored")
  return false
end

# notice of the reactance curve of the phase tap changer `tc` at its step;
# a warning when xMax/xMin leave the plausibility band of the tap
# corrections (FullGrid carries the placeholder pair 0.99/9.99), then the
# curve is ignored and the end keeps its reactance
function _phaseTapXCurveNotice!(store::CGMESStore, svsteps::Dict{String,Float64}, tc::CIMObject, messages::Vector{String}; strict::Bool = false)
  tc.class in (:PhaseTapChangerLinear, :PhaseTapChangerSymmetrical, :PhaseTapChangerAsymmetrical) || return nothing
  xmin, xmax = num(tc, :xMin), num(tc, :xMax)
  (xmin === nothing || xmax === nothing || xmin <= 0.0 || xmax == xmin) && return nothing
  name = something(str(tc, :name), tc.mrid)
  if !(CGMES_TAP_CORR_MIN <= xmax / xmin <= CGMES_TAP_CORR_MAX)
    msg = "$(tc.class) $(name) declares an implausible reactance curve (xMin $(xmin), xMax $(xmax), band of xMax/xMin $(CGMES_TAP_CORR_MIN)..$(CGMES_TAP_CORR_MAX)), placeholder suspected"
    strict && error("cgmes_import.placeholder_guards = strict: " * msg)
    push!(messages, "warning: " * msg * ", reactance curve ignored")
    return nothing
  end
  step = _tapStep(tc, svsteps)
  f = _phaseTapXCurveFactor(store, tc, step)
  f == 1.0 || push!(messages, "notice: $(tc.class) $(name): step $(step), series reactance $(round(f; digits = 5)) times the end value (xMin/xMax curve)")
  return nothing
end

# otherwise it divides (referred through the winding).
function _applyEndTaps(store::CGMESStore, svsteps::Dict{String,Float64}, e::CIMObject, ratio::Float64, shift::Float64, on_from_side::Bool, messages::Vector{String}; strict::Bool = false)::Tuple{Float64,Float64}
  for tc in _tapChangersOfEnd(store, e.mrid)
    if tc.class == :RatioTapChanger
      table = _ratioTapTable(store, tc)
      step = _tapStep(tc, svsteps)
      corr, from_table = _ratioTapCorrectionAt(tc, table, step)
      (table === nothing || from_table) || push!(messages, "notice: RatioTapChanger $(something(str(tc, :name), tc.mrid)): RatioTapChangerTable does not cover step $(step), linear step formula used")
      _ratioTapTableHasImpedance(store, tc) && push!(messages, "notice: RatioTapChanger $(something(str(tc, :name), tc.mrid)): RatioTapChangerTable carries per-step r/x/g/b corrections, the point of the step applies to the impedance")
      _tapCorrPlausible("RatioTapChanger", tc, corr, messages; strict = strict) || continue
      ratio = on_from_side ? ratio * corr : ratio / corr
    else
      if tc.class == :PhaseTapChangerTabular
        built = _phaseTapTableModel(store, tc, Int(round(_tapStep(tc, svsteps))))
        if built === nothing
          push!(messages, "notice: PhaseTapChangerTabular $(something(str(tc, :name), tc.mrid)), table unresolved, tap ignored")
          continue
        end
        row = Sparlectra.calcPhaseTapTable(built.model)
        push!(messages, "tabular phase tap $(something(str(tc, :name), tc.mrid)): step $(built.model.step) → ratio $(round(row.effective_ratio; digits = 5)), shift $(round(row.effective_shift_deg; digits = 3))°$(built.impedance_note ? " (table with r/x/g/b corrections, the row of the step applies to the impedance)" : "")")
      end
      pratio, pshift = _phaseTapRatioShift(store, tc, svsteps)
      _tapCorrPlausible(String(tc.class), tc, pratio, messages; strict = strict) || continue
      _phaseTapXCurveNotice!(store, svsteps, tc, messages; strict = strict)
      ratio = on_from_side ? ratio * pratio : ratio / pratio
      # Angle sign: an end-2 (to-side) tap angle enters NEGATED: the branch
      # ratio is t1/t2, so θ_eff = θ1 − θ2 (standard end-referral semantics;
      # PowSyBl does the same when moving an end-2 changer to end 1). The
      # deciding fixture is RealGrid: its four flowing end-2 tabular PSTs
      # reproduce the SV state only with the flip (SV-start max|F| 5 pu vs
      # 395 pu without). KNOWN AMBIGUITY: the ENTSO-E PSEI pair
      # PST_PTChLin_PTE1/PTE2 encodes identical parameters on end 1/end 2 and
      # its PTE2 SV expects the UNflipped angle. Under this (standard-side)
      # convention the PTE2 sets deviate by dva 0.29° / dp 0.79 MW. Real
      # deliveries win over the 2-bus conformity toys.
      shift += on_from_side ? pshift : -pshift
    end
  end
  return (ratio, shift)
end

# --- bus creation (lazy: only TNs referenced by mapped equipment) -----------

# An SvVoltage with v <= 0 is left out, so the node counts as one without SV:
# _ensureBus! gives it the 1.0 pu / 0° start (later replaced by the start of
# a neighbour, _startNoSvBusesFromNeighbours!), the three-winding star start
# does not take it, and it lands in no_sv_buses. compareWithSV and the
# measurement table skip such entries in the same way. Real deliveries carry
# energized nodes with v = 0 in the SV; taken as given, they started the
# solve at zero voltage.
function _svVoltageByTN(store::CGMESStore)::Dict{String,Tuple{Float64,Float64}}
  out = Dict{String,Tuple{Float64,Float64}}()
  for sv in objectsOf(store, :SvVoltage)
    tn = get(sv.refs, :TopologicalNode, nothing)
    tn === nothing && continue
    v = num(sv, :v)
    a = num(sv, :angle, 0.0)
    (v === nothing || v <= 0.0) && continue
    out[tn] = (v, a)
  end
  return out
end

function _ensureBus!(net::Sparlectra.Net, created::Set{String}, topo::CGMESTopology, svmap::Dict{String,Tuple{Float64,Float64}}, tn::String; isAux::Bool = false, busname::Union{Nothing,String} = nothing)
  # busname overrides the TN's bus for per-side split buses (see
  # _detectNonCancellingBoundarySides!): same TN attributes, different identity.
  bus = something(busname, topo.bus_name[tn])
  bus in created && return bus
  vn = topo.vn_kV[tn]
  vm_pu, va_deg = 1.0, 0.0
  if haskey(svmap, tn) && vn > 0.0
    v, a = svmap[tn]
    vm_pu = v / vn
    va_deg = a
  end
  Sparlectra.addBus!(net = net, busName = bus, vn_kV = vn, vm_pu = vm_pu, va_deg = va_deg, isAux = isAux)
  push!(created, bus)
  # the imported state, kept apart from the start value that every solve
  # overwrites: the state-estimation range check substitutes it for an
  # implausible voltage measurement at this bus
  (haskey(svmap, tn) && vn > 0.0) && (net._sv_vm_pu[net.busDict[bus]] = vm_pu)
  # Identity capture for a later CGMES export: remember the original TN mRID
  # under the structural bus key. Aux buses and per-side split buses are
  # synthetic identities without a 1:1 source TN, the exporter mints fresh
  # deterministic ids for them instead.
  if !isAux && busname === nothing
    net.cgmes_ids[cgmesKeyTopologicalNode(bus)] = cgmesCanonicalMrid(tn)
  end
  return bus
end

"""
Detect assembled boundary nodes whose equivalent injections do **not** cancel
and assign each contributing area (source file) its own side bus.

An AC X-node with both sides present carries two declarations of the *same*
exchange: they cancel, the equivalents are discarded, and the node joins the
areas galvanically. When the pair does NOT cancel, the two injections are two
different devices (in practice: the HVDC converters of a DC border crossing;
the residue is the DC loss). Joining the sides galvanically is then wrong
twice over: it creates an AC corridor where none exists, and it violates KCL
against the delivery's own SV state (ReliCapGrid `BP_SD-EH_DC1`: 94.8 MW flow
into the node on one side, 72.6 MW out on the other, impossible through one
node). Each side therefore keeps its line and its converter injection on its
own bus.

The side of a piece of equipment is the file that defined it, the only
criterion that exists, since both sides' terminals reference the same CN (and
here even the same TN). The first source (sorted) keeps the base bus name,
further sides get `<bus>@<k>`.
"""
function _detectNonCancellingBoundarySides!(ctx::_MapCtx, store::CGMESStore, topo::CGMESTopology)
  by_cn = Dict{String,Vector{Tuple{String,Float64}}}()   # cn → [(source, p)]
  for ei in objectsOf(store, :EquivalentInjection)
    _inService(ctx, ei) || continue
    ts = get(topo.terminals, ei.mrid, nothing)
    (ts === nothing || isempty(ts)) && continue
    cn = get(ts[1].refs, :ConnectivityNode, nothing)
    (cn === nothing || !(cn in store.boundary)) && continue
    push!(get!(Vector{Tuple{String,Float64}}, by_cn, cn), (ei.source, num(ei, :p, 0.0)))
  end
  for (cn, members) in by_cn
    length(members) >= 2 || continue
    _eiPairCancels([m[2] for m in members]) && continue
    sources = sort(unique(first.(members)))
    length(sources) >= 2 || continue
    # find the base bus name via any terminal on this CN that resolves to a TN
    base = nothing
    for t in objectsOf(store, :Terminal)
      get(t.refs, :ConnectivityNode, nothing) == cn || continue
      tn = tnOfTerminal(t)
      (tn === nothing || !haskey(topo.bus_name, tn)) && continue
      base = topo.bus_name[tn]
      break
    end
    base === nothing && continue
    sides = Dict{String,String}()
    for (k, src) in enumerate(sources)
      sides[src] = k == 1 ? base : string(base, "@", k)
    end
    ctx.split_sides[cn] = sides
    push!(
      ctx.messages,
      "notice: boundary node $(base), the two sides' equivalents do not cancel (device injections, e.g. HVDC converters); node split per side: $(join([string(basename_short(s), " → ", b) for (s, b) in sort(collect(sides))], ", "))",
    )
  end
  return nothing
end

# short human-readable tag for a source file name in messages
basename_short(s::AbstractString) = first(splitext(basename(String(s))))

"""
Bus for one terminal of `eq`, honoring the per-side split of non-cancelling
boundary nodes: when the terminal's CN is split, the equipment lands on its
own area's side bus; everything else falls back to the TN bus.
"""
function _busForTerminal!(net, created, topo::CGMESTopology, svmap, ctx::_MapCtx, eq::CIMObject, seq::Int, tn::String)::String
  if !isempty(ctx.split_sides)
    ts = get(topo.terminals, eq.mrid, nothing)
    if ts !== nothing && length(ts) >= seq
      cn = get(ts[seq].refs, :ConnectivityNode, nothing)
      if cn !== nothing && haskey(ctx.split_sides, cn)
        side = get(ctx.split_sides[cn], eq.source, nothing)
        side === nothing || return _ensureBus!(net, created, topo, svmap, tn; busname = side)
      end
    end
  end
  return _ensureBus!(net, created, topo, svmap, tn)
end

# --- element mapping --------------------------------------------------------

# Retained switches in the bus-branch view: a Switch-family object whose two
# terminals sit on different TopologicalNodes couples those buses (e.g. the
# NL busbar coupler in MicroGrid). Closed → zero-impedance `addLink!`;
# open → ignored. SSH `open` overrides the EQ default `normalOpen`.
const _SWITCH_CLASSES = (:Switch, :Breaker, :Disconnector, :LoadBreakSwitch, :Jumper, :ProtectedSwitch)

# THE switch-openness rule, shared verbatim by the bus-branch link mapping
# below and the node-breaker topology processor (#314): SSH `open`
# overrides the EQ default `normalOpen`, and an out-of-service switch
# cannot conduct whatever its open flag says.
_switchIsOpen(ctx::_MapCtx, sw::CIMObject)::Bool = something(boolval(sw, :open), boolval(sw, :normalOpen), false) || !_inService(ctx, sw)

function _mapSwitches!(net, store, topo, created, svmap, ctx::_MapCtx)
  for cls in _SWITCH_CLASSES
    for sw in objectsOf(store, cls)
      t1 = busOfEquipment(topo, sw, 1)
      t2 = busOfEquipment(topo, sw, 2)
      (t1 === nothing || t2 === nothing) && continue
      t1.tn == t2.tn && continue
      name = something(str(sw, :name), sw.mrid)
      # an out-of-service switch cannot conduct, whatever its open flag says
      isopen = something(boolval(sw, :open), boolval(sw, :normalOpen), false) || !_inService(ctx, sw)
      if isopen || !(_conn(ctx, t1.connected) && _conn(ctx, t2.connected))
        push!(ctx.messages, "notice: open $(cls) $(name) between $(t1.bus) and $(t2.bus), buses stay separate")
        continue
      end
      b1 = _ensureBus!(net, created, topo, svmap, t1.tn)
      b2 = _ensureBus!(net, created, topo, svmap, t2.tn)
      lidx = Sparlectra.addLink!(net = net, fromBus = b1, toBus = b2, status = 1)
      ctx.element_by_mrid[sw.mrid] = (:link, lidx)
      # identity capture: the link re-exports as a closed Breaker with the
      # original switch mRID
      k = cgmesNextBusEquipmentIndex!(ctx.eqp_counters, "LNK", _cgmesTrafoPairKey(b1, b2))
      lnkKey = cgmesKeyBusEquipment("LNK", _cgmesTrafoPairKey(b1, b2), k)
      net.cgmes_ids[lnkKey] = cgmesCanonicalMrid(sw.mrid)
      ts = get(topo.terminals, sw.mrid, nothing)
      if ts !== nothing
        length(ts) >= 1 && (net.cgmes_ids[cgmesKeyTerminal(lnkKey, 1)] = cgmesCanonicalMrid(ts[1].mrid))
        length(ts) >= 2 && (net.cgmes_ids[cgmesKeyTerminal(lnkKey, 2)] = cgmesCanonicalMrid(ts[2].mrid))
      end
      push!(ctx.messages, "notice: closed $(cls) $(name), bus link $(b1) = $(b2)")
    end
  end
end

# mRID of the first phase tap changer on a transformer end, if any
function _firstPhaseTapChangerMrid(store::CGMESStore, endmrid::String)::Union{Nothing,String}
  for tc in _tapChangersOfEnd(store, endmrid)
    tc.class == :RatioTapChanger && continue
    return tc.mrid
  end
  return nothing
end

# mRID of the equipment's sequence-1 terminal, if the topology knows one
function _firstTerminalMrid(topo::CGMESTopology, obj::CIMObject)::Union{Nothing,String}
  ts = get(topo.terminals, obj.mrid, nothing)
  (ts === nothing || isempty(ts)) && return nothing
  return ts[1].mrid
end

# Structural-identity capture for single-terminal equipment. The (kind, bus)
# counter is consumed for EVERY created unit of the class, also when there is
# no mRID to record (mrid = nothing), so the export, which walks the created
# vectors in the same order, assigns identical <k> values. `children` records
# dependent mRIDs under `<key>|<suffix>` (T1 terminal, GU, RC).
function _captureBusEquipment!(net, ctx::_MapCtx, kind::AbstractString, bus::AbstractString, mrid::Union{Nothing,String}; children::NamedTuple = NamedTuple())
  k = cgmesNextBusEquipmentIndex!(ctx.eqp_counters, kind, bus)
  mrid === nothing && return nothing
  key = cgmesKeyBusEquipment(kind, bus, k)
  net.cgmes_ids[key] = cgmesCanonicalMrid(mrid)
  for (suffix, cmrid) in pairs(children)
    cmrid === nothing && continue
    net.cgmes_ids[string(key, "|", suffix)] = cgmesCanonicalMrid(cmrid)
  end
  return nothing
end

# record the (branchIdx, side) provenance of the just-added branch for the
# two terminals of a line-like equipment (sequence 1 = from, 2 = to)
function _recordBranchTerminals!(ctx::_MapCtx, net, topo::CGMESTopology, eq::CIMObject)
  ts = get(topo.terminals, eq.mrid, nothing)
  ts === nothing && return nothing
  idx = net.branchVec[end].branchIdx
  length(ts) >= 1 && (ctx.branch_side[ts[1].mrid] = (idx, :from))
  length(ts) >= 2 && (ctx.branch_side[ts[2].mrid] = (idx, :to))
  return nothing
end

function _mapLines!(net, store, topo, created, svmap, baseMVA, ctx::_MapCtx)
  # First-seen parallel index per bus pair: must count exactly the lines
  # that reach net.linesAC so the exporter (walking net.linesAC in the same
  # order) reproduces identical structural keys.
  acl_parallel = Dict{Tuple{String,String},Int}()
  for cls in (:ACLineSegment, :SeriesCompensator)
    for line in objectsOf(store, cls)
      t1 = busOfEquipment(topo, line, 1)
      t2 = busOfEquipment(topo, line, 2)
      name = something(str(line, :name), line.mrid)
      if t1 === nothing || t2 === nothing
        push!(ctx.messages, "skip: $(cls) $(name), terminal without TopologicalNode")
        push!(ctx.skipped, line.mrid)
        continue
      end
      # Terminal.connected handling (r0.9.10). Both ends open (or
      # Equipment.inService == false, a de-commissioned line contributes no
      # charging either) → the branch is out of service. Exactly ONE end
      # open → the branch is created with per-terminal flags and reduces to
      # its exact pi Schur complement at the closed bus: it draws its FULL
      # charging there (the earlier substitute shunt with HALF the charging
      # undercounted it by about a factor of two; RealGrid carries 511 such
      # lines against 228 fully open ones).
      in_service = _inService(ctx, line)
      c1 = in_service && _conn(ctx, t1.connected)
      c2 = in_service && _conn(ctx, t2.connected)
      nopen = (c1 ? 0 : 1) + (c2 ? 0 : 1)
      if nopen == 1
        push!(ctx.messages, "notice: $(cls) $(name) has one open terminal, kept as a partially open branch (full charging at the closed bus)")
      end
      status = nopen == 0 ? 1 : 0
      b1 = _busForTerminal!(net, created, topo, svmap, ctx, line, 1, t1.tn)
      b2 = _busForTerminal!(net, created, topo, svmap, ctx, line, 2, t2.tn)
      if b1 == b2
        # Self-loop within one topological node (e.g. a busbar link modeled
        # as a line): there is no longitudinal element to build, but the
        # charging still acts at the bus. Keep it as a shunt instead of
        # tripping the branch constructor's from!=to assertion.
        if cls == :ACLineSegment && status == 1
          g_full = num(line, :gch, 0.0)
          b_full = num(line, :bch, 0.0)
          if g_full != 0.0 || b_full != 0.0
            Sparlectra.addShuntMatpower!(net = net, busName = b1, Gs = g_full * t2.vn_kV^2, Bs = b_full * t2.vn_kV^2)
            # the line itself builds no element; its charging is addressed
            # under a derived key
            ctx.element_by_mrid[string(line.mrid, "#charging")] = (:shunt, length(net.shuntVec))
          end
        end
        push!(ctx.messages, "notice: $(cls) $(name) connects bus $(b1) to itself (both terminals on one topological node), no branch created, charging kept as shunt")
        push!(ctx.skipped, line.mrid)
        continue
      end
      # branch model convention (calcAdmittance): series/shunt admittance on
      # the TO-side voltage base, complex ratio t at the from side
      # a CGMES line carries one charging admittance (bch, gch), so the two
      # terminal arms are its halves: the symmetric pi of every builder
      g_half = cls == :ACLineSegment ? 0.5 * num(line, :gch, 0.0) : 0.0
      b_half = cls == :ACLineSegment ? 0.5 * num(line, :bch, 0.0) : 0.0
      pu = Sparlectra.pi_branch_pu_between_levels(
        r = num(line, :r, 0.0),
        x = num(line, :x, 0.0),
        g_from = g_half,
        b_from = b_half,
        g_to = g_half,
        b_to = b_half,
        vn_from_kV = t1.vn_kV,
        vn_to_kV = t2.vn_kV,
        baseMVA = baseMVA,
      )
      r_pu, x_pu, b_pu, g_pu = pu.r_pu, pu.x_pu, pu.b_pu, pu.g_pu
      if pu.ratio === nothing
        Sparlectra.addPIModelACLine!(net = net, fromBus = b1, toBus = b2, r_pu = r_pu, x_pu = x_pu, b_pu = b_pu, g_pu = g_pu, status = status, from_status = c1 ? 1 : 0, to_status = c2 ? 1 : 0)
        ctx.element_by_mrid[line.mrid] = (:branch, length(net.branchVec))
        _recordBranchTerminals!(ctx, net, topo, line)
        # Identity capture (line + its two terminals) for a later CGMES
        # export. Terminal keys use the equipment sequence (T1 = from side),
        # independent of the lexicographic pair order inside the key.
        k = cgmesNextParallelIndex!(acl_parallel, b1, b2)
        aclKey = cgmesKeyACLineSegment(b1, b2, k)
        net.cgmes_ids[aclKey] = cgmesCanonicalMrid(line.mrid)
        ts = get(topo.terminals, line.mrid, nothing)
        if ts !== nothing
          length(ts) >= 1 && (net.cgmes_ids[cgmesKeyTerminal(aclKey, 1)] = cgmesCanonicalMrid(ts[1].mrid))
          length(ts) >= 2 && (net.cgmes_ids[cgmesKeyTerminal(aclKey, 2)] = cgmesCanonicalMrid(ts[2].mrid))
        end
      else
        # boundary line across different nominal voltages (e.g. 380 kV area
        # bus vs 400 kV boundary node): identical physical conductor, so the
        # 2WT formula with ratedU1 == ratedU2 degenerates to ratio = vn2/vn1
        # on a PI branch (impedance on the to-side base).
        ratio = pu.ratio
        from = Sparlectra.geNetBusIdx(net = net, busName = b1)
        to = Sparlectra.geNetBusIdx(net = net, busName = b2)
        Sparlectra._addPIModelTrafo_by_idx!(net = net, from = from, to = to, r_pu = r_pu, x_pu = x_pu, b_pu = b_pu, g_pu = g_pu, status = status, ratio = ratio, shift_deg = 0.0, from_status = c1 ? 1 : 0, to_status = c2 ? 1 : 0)
        _attachNoSETapChangers!(net.branchVec[end])
        # still a line in CGMES terms, so a branch, not a transformer
        ctx.element_by_mrid[line.mrid] = (:branch, length(net.branchVec))
        _recordBranchTerminals!(ctx, net, topo, line)
        # the ratio branch counts as a transformer branch on this bus pair,
        # consume the PT counter (no capture: the source object is a line)
        # so real transformers on the same pair keep consistent <k> values
        cgmesNextBusEquipmentIndex!(ctx.eqp_counters, "PT", _cgmesTrafoPairKey(b1, b2))
        push!(ctx.messages, "notice: $(cls) $(name) spans $(t1.vn_kV)/$(t2.vn_kV) kV, mapped as boundary-line ratio branch ($(round(ratio; digits=4)))")
      end
    end
  end
  # the equivalent branches continue this counter (_mapEquivalentBranches!)
  return acl_parallel
end

# Phase stamp of a non-reciprocal EquivalentBranch (z12 from terminal 1, z21
# from terminal 2, per unit on one base) on the branch model of
# calcAdmittance: |t| = 1, phi = arg(z21/z12)/2 and the series admittance
# ys = c e^{-j phi} with c = sqrt(|z12|/|z21|)/z12, so |zs| is the geometric
# mean sqrt(|z12| |z21|) and zs has the angle of z12 e^{j phi}. The arms
# y0f = c (1 - e^{-j phi}) and y0t = c e^{-j phi} (e^{-j phi} - 1) make both
# row sums of the two-port zero: Y11 = -Y12 = c, Y22 = -Y21 = c e^{-2j phi}.
# For |z21| = |z12| that is the two-port I1 = (V1 - V2)/z12,
# I2 = (V2 - V1)/z21 exactly. A branch always has |Y12| = |Y21|, so
# otherwise every entry is off in magnitude by sqrt(|z21|/|z12|) or its
# inverse. With zero row sums that error acts on V1 - V2, of the order
# |V1 - V2| delta / (2 |z|) per terminal with delta = ||z21|/|z12| - 1|.
# Putting the whole error on Y21 alone (Y11, Y12, Y22 exact) acts on V1
# instead and is of the order of the short-circuit current.
# _ebFromPhaseStamp is the inverse; it gives back z12 and z21 with the same
# geometric mean of the magnitudes, their ratio is not in the branch. The
# delivered values stay in net.equivalent_branch_metadata for the export.
function _ebPhaseStamp(z12::ComplexF64, z21::ComplexF64)
  phi = angle(z21 / z12) / 2
  c = sqrt(abs(z12) / abs(z21)) / z12
  return (phi = phi, zs = cis(phi) / c, y0f = c * (1 - cis(-phi)), y0t = c * cis(-phi) * (cis(-phi) - 1))
end

# z12 and z21 back from a stamped branch (shift phi in rad, series zs, arms
# y0f and y0t), both with the magnitude sqrt(|z12| |z21|); `nothing` when the
# arms do not fit the stamp of zs and phi, so the branch is no stamped
# EquivalentBranch (any more)
function _ebFromPhaseStamp(phi::Float64, zs::ComplexF64, y0f::ComplexF64, y0t::ComplexF64)
  iszero(zs) && return nothing
  c = cis(phi) / zs
  tol = 1e-9 * abs(c)
  abs(y0f - c * (1 - cis(-phi))) <= tol || return nothing
  abs(y0t - c * cis(-phi) * (cis(-phi) - 1)) <= tol || return nothing
  return (z12 = 1 / c, z21 = cis(2 * phi) / c)
end

# EquivalentBranch: the series impedance of an equivalent network between two
# terminals, r and x in ohm, no charging. Unmapped, the buses of the
# equivalent network that hang on such branches alone end up isolated and
# the solve runs on another network than the delivery's SV state.
# - Mapped only with both terminals connected and the object in service;
#   otherwise it goes to ctx.skipped like any other skipped equipment.
# - Both terminals on one TopologicalNode: no element, one notice each.
# - Same nominal voltage at both terminals: a pi branch without charging, r
#   and x on that voltage base.
# - Different nominal voltages: a branch with per-unit ratio 1 and the
#   impedance on the base of terminal 1, which is what the SvPowerFlow of
#   such objects in real deliveries matches. The boundary-line rule of
#   _mapLines! (ratio vn2/vn1, impedance on the terminal-2 base) misses those
#   flows by an order of magnitude.
# - r21 and x21 (the impedance seen from terminal 2) are optional in CGMES and
#   equal to r and x when absent. With other values the object is the
#   two-port I1 = (V1 - V2)/z12, I2 = (V2 - V1)/z21, which fits the
#   SvPowerFlow of real deliveries best. It is built as a transformer branch
#   with the phase stamp of _ebPhaseStamp on both level kinds: exact for
#   |z21| = |z12|, otherwise built with the geometric mean of the two
#   magnitudes, which the notice counts per branch.
# Identity capture: a line continues the line counter of its bus pair
# (`acl_parallel` from _mapLines!) under `EQB|<busA>|<busB>|<k>`, a
# transformer branch (across levels or stamped) the transformer counter under
# `EQBT|<busA>|<busB>|<k>`, so the exporter, which walks the lines and the
# transformer branches in the same order, finds the key and writes the
# branch back as EquivalentBranch.
function _mapEquivalentBranches!(net, store, topo, created, svmap, baseMVA, ctx::_MapCtx, acl_parallel::Dict{Tuple{String,String},Int})
  ebs = collect(CIMObject, objectsOf(store, :EquivalentBranch))
  isempty(ebs) && return nothing
  n_same = 0
  n_cross = 0
  n_loop = 0
  n_open = 0
  n_noterm = 0
  n_link = 0
  n_stamp = 0
  n_mag = 0
  max_mag = 0.0
  n_asym_zero = 0
  for eb in ebs
    t1 = busOfEquipment(topo, eb, 1)
    t2 = busOfEquipment(topo, eb, 2)
    name = something(str(eb, :name), eb.mrid)
    if t1 === nothing || t2 === nothing
      push!(ctx.messages, "skip: EquivalentBranch $(name), terminal without TopologicalNode")
      push!(ctx.skipped, eb.mrid)
      n_noterm += 1
      continue
    end
    if t1.tn == t2.tn
      push!(ctx.messages, "notice: EquivalentBranch $(name) connects bus $(t1.bus) to itself (both terminals on one topological node), no branch created")
      push!(ctx.skipped, eb.mrid)
      n_loop += 1
      continue
    end
    if !(_inService(ctx, eb) && _conn(ctx, t1.connected) && _conn(ctx, t2.connected))
      push!(ctx.skipped, eb.mrid)
      n_open += 1
      continue
    end
    r = num(eb, :r, 0.0)
    x = num(eb, :x, 0.0)
    z12 = complex(r, x)
    z21 = complex(something(num(eb, :r21), r), something(num(eb, :x21), x))
    # an impedance below cgmes_import.equivalent_branch_min_z_pu (both
    # directions, per unit on the base of terminal 1) builds no branch: the
    # two buses are joined as a link, like a closed switch. The link takes
    # the next LNK counter of its pair without capture, so the switches the
    # export walks after it keep their keys.
    zb1 = t1.vn_kV^2 / baseMVA
    if max(abs(z12), abs(z21)) / zb1 < ctx.eb_min_z_pu
      b1 = _busForTerminal!(net, created, topo, svmap, ctx, eb, 1, t1.tn)
      b2 = _busForTerminal!(net, created, topo, svmap, ctx, eb, 2, t2.tn)
      lidx = Sparlectra.addLink!(net = net, fromBus = b1, toBus = b2, status = 1)
      ctx.element_by_mrid[eb.mrid] = (:link, lidx)
      cgmesNextBusEquipmentIndex!(ctx.eqp_counters, "LNK", _cgmesTrafoPairKey(b1, b2))
      push!(ctx.messages, "notice: EquivalentBranch $(name) has an impedance below $(ctx.eb_min_z_pu) pu, no branch created, bus link $(b1) = $(b2)")
      n_link += 1
      continue
    end
    asym = abs(z21 - z12) > 1e-9 * abs(z12)
    # a zero impedance on one side has no angle to stamp; such an object keeps
    # the symmetric branch with r and x and is counted apart
    stamped = asym && !iszero(z12) && !iszero(z21)
    asym && !stamped && (n_asym_zero += 1)
    b1 = _busForTerminal!(net, created, topo, svmap, ctx, eb, 1, t1.tn)
    b2 = _busForTerminal!(net, created, topo, svmap, ctx, eb, 2, t2.tn)
    # both impedances on the base of terminal 1, also across levels (per-unit
    # ratio 1)
    r_pu, x_pu, _, _ = Sparlectra.toPU_RXBG(r = r, x = x, v_kv = t1.vn_kV, baseMVA = baseMVA)
    if t1.vn_kV == t2.vn_kV && !stamped
      Sparlectra.addPIModelACLine!(net = net, fromBus = b1, toBus = b2, r_pu = r_pu, x_pu = x_pu, b_pu = 0.0, g_pu = 0.0, status = 1)
      k = cgmesNextParallelIndex!(acl_parallel, b1, b2)
      eqbKey = cgmesKeyEquivalentBranch(b1, b2, k)
    else
      from = Sparlectra.geNetBusIdx(net = net, busName = b1)
      to = Sparlectra.geNetBusIdx(net = net, busName = b2)
      if stamped
        r21_pu, x21_pu, _, _ = Sparlectra.toPU_RXBG(r = real(z21), x = imag(z21), v_kv = t1.vn_kV, baseMVA = baseMVA)
        st = _ebPhaseStamp(complex(r_pu, x_pu), complex(r21_pu, x21_pu))
        Sparlectra._addPIModelTrafo_by_idx!(net = net, from = from, to = to, r_pu = real(st.zs), x_pu = imag(st.zs), b_pu = 0.0, g_pu = 0.0, status = 1, ratio = 1.0, shift_deg = rad2deg(st.phi),
          g_from_pu = real(st.y0f), b_from_pu = imag(st.y0f), g_to_pu = real(st.y0t), b_to_pu = imag(st.y0t))
        n_stamp += 1
        dev = abs(abs(z21) / abs(z12) - 1.0)
        if dev > 1e-9
          n_mag += 1
          max_mag = max(max_mag, dev)
        end
      else
        Sparlectra._addPIModelTrafo_by_idx!(net = net, from = from, to = to, r_pu = r_pu, x_pu = x_pu, b_pu = 0.0, g_pu = 0.0, status = 1, ratio = 1.0, shift_deg = 0.0)
      end
      _attachNoSETapChangers!(net.branchVec[end])
      k = cgmesNextBusEquipmentIndex!(ctx.eqp_counters, "PT", _cgmesTrafoPairKey(b1, b2))
      eqbKey = cgmesKeyEquivalentBranchTransformer(b1, b2, k)
    end
    t1.vn_kV == t2.vn_kV ? (n_same += 1) : (n_cross += 1)
    # the delivered values, for the export and as the equivalent-branch
    # marker of the results; the solver reads the branch only
    net.equivalent_branch_metadata[net.branchVec[end].branchIdx] = (r = r, x = x, r21 = real(z21), x21 = imag(z21), vn_kV = t1.vn_kV, stamped = stamped)
    ctx.element_by_mrid[eb.mrid] = (:branch, length(net.branchVec))
    _recordBranchTerminals!(ctx, net, topo, eb)
    net.cgmes_ids[eqbKey] = cgmesCanonicalMrid(eb.mrid)
    ts = get(topo.terminals, eb.mrid, nothing)
    if ts !== nothing
      length(ts) >= 1 && (net.cgmes_ids[cgmesKeyTerminal(eqbKey, 1)] = cgmesCanonicalMrid(ts[1].mrid))
      length(ts) >= 2 && (net.cgmes_ids[cgmesKeyTerminal(eqbKey, 2)] = cgmesCanonicalMrid(ts[2].mrid))
    end
  end
  n_skip = n_loop + n_open + n_noterm
  asym = n_stamp == 0 ? "" : ", $(n_stamp) asymmetric (r21, x21 differ from r, x) phase-stamped (ratio 1, shift arg(z21/z12)/2), $(n_mag) of them with |z21| other than |z12|"
  n_mag == 0 || (asym *= " (max ||z21|/|z12| - 1| = $(round(max_mag; sigdigits = 3)), built with the geometric mean of |z12| and |z21|)")
  n_asym_zero == 0 || (asym *= ", $(n_asym_zero) asymmetric with a zero impedance, mapped symmetric with r, x")
  n_link == 0 || (asym *= ", $(n_link) below $(ctx.eb_min_z_pu) pu as bus links")
  push!(ctx.messages, "notice: EquivalentBranch: $(n_same + n_cross) mapped ($(n_same) same voltage level, $(n_cross) across voltage levels), $(n_skip) skipped ($(n_loop) on one topological node, $(n_open) open or out of service, $(n_noterm) without topological node)$(asym)")
  return nothing
end

function _endsOfTransformer(store::CGMESStore, pt::CIMObject)::Vector{CIMObject}
  ends = [e for e in objectsOf(store, :PowerTransformerEnd) if get(e.refs, :PowerTransformer, "") == pt.mrid]
  sort!(ends; by = e -> something(num(e, :endNumber), 1.0))
  return ends
end

function _terminalInfo(store::CGMESStore, topo::CGMESTopology, e::CIMObject)
  tmrid = get(e.refs, :Terminal, nothing)
  tmrid === nothing && return nothing
  t = get(store.objects, tmrid, nothing)
  t === nothing && return nothing
  tn = tnOfTerminal(t)
  (tn === nothing || !haskey(topo.bus_name, tn)) && return nothing
  return (bus = topo.bus_name[tn], vn_kV = topo.vn_kV[tn], connected = terminalConnected(t), tn = tn)
end

# A tap changer sits at the terminal of its own end and the series impedance
# lies on the inner side of the tap (the convention the SV of real deliveries
# and OpenLoadFlow follow). The branch model has one ideal transformer, at
# the from side, ahead of the impedance on the to-side base, so the end that
# carries the tap changers is the from side. A transformer with tap changers
# only on end 2 is therefore built from end 2 to end 1. With tap changers on
# both ends end 1 stays the from side and the end-2 correction enters the
# same ratio (its impedance then sits on the wrong side of that tap).
_trafoEndHasTaps(store::CGMESStore, e::CIMObject)::Bool = !isempty(_tapChangersOfEnd(store, e.mrid))

function _map2WTrafo!(net, store, topo, created, svmap, svsteps, baseMVA, pt, ends, ctx::_MapCtx)
  e1, e2 = ends[1], ends[2]
  i1 = _terminalInfo(store, topo, e1)
  i2 = _terminalInfo(store, topo, e2)
  name = something(str(pt, :name), pt.mrid)
  if i1 === nothing || i2 === nothing
    push!(ctx.messages, "skip: PowerTransformer $(name), end terminal without TopologicalNode")
    push!(ctx.skipped, pt.mrid)
    return
  end
  taps1 = _trafoEndHasTaps(store, e1)
  taps2 = _trafoEndHasTaps(store, e2)
  taps1 && taps2 && push!(ctx.messages, "notice: PowerTransformer $(name) has tap changers on both ends, end 1 is the branch from side and the end-2 tap acts at the end-1 terminal")
  # (ef, et): the ends on the branch's from and to side
  turned = taps2 && !taps1
  ef, et = turned ? (e2, e1) : (e1, e2)
  i_f, i_t = turned ? (i2, i1) : (i1, i2)
  Uf = something(num(ef, :ratedU), i_f.vn_kV)
  Ut = something(num(et, :ratedU), i_t.vn_kV)
  # refer both end impedances to the to-side end: the branch model
  # (calcAdmittance) expects series/shunt admittance on the TO-side base,
  # ratio at the from side
  z = Sparlectra.calc2WTEndsReferredRXGB(;
    r1 = num(ef, :r, 0.0), x1 = num(ef, :x, 0.0), g1 = num(ef, :g, 0.0), b1 = num(ef, :b, 0.0),
    r2 = num(et, :r, 0.0), x2 = num(et, :x, 0.0), g2 = num(et, :g, 0.0), b2 = num(et, :b, 0.0),
    U1 = Uf, U2 = Ut,
  )
  # tap table r/x/g/b of both ends act on the transformer's series
  # impedance and magnetizing admittance as a whole
  fr1, fx1, fg1, fb1 = _tapTableImpedanceFactors(store, svsteps, ef)
  fr2, fx2, fg2, fb2 = _tapTableImpedanceFactors(store, svsteps, et)
  fr, fx, fg, fb = fr1 * fr2, fx1 * fx2, fg1 * fg2, fb1 * fb2
  r_pu, x_pu, b_pu, g_pu = Sparlectra.toPU_RXBG(r = z.r * fr, x = z.x * fx, g = z.g * fg, b = z.b * fb, v_kv = i_t.vn_kV, baseMVA = baseMVA)
  # magnetizing admittance: half of the total of both ends on each terminal
  # arm (0.41.0), the from arm behind the ratio. The SV of SmallGrid and
  # of MicroGrid BE and NL is computed this way (unit rows at the
  # transformer ends 0.000 MVAr); the per-end placement of 0.20.0 left up to
  # 3.6 MVAr there. The PowSyBl import keeps the per-end placement, the 3W
  # legs put the end admittance at the star node (_map3WTrafo!).
  g_from_pu = g_to_pu = 0.5 * g_pu
  b_from_pu = b_to_pu = 0.5 * b_pu
  ratio = (Uf / Ut) / (i_f.vn_kV / i_t.vn_kV)
  # the structural ratio, the neutral point of the tap estimation's set
  ratio_struct = ratio
  shift = 0.0
  ratio, shift = _applyEndTaps(store, svsteps, ef, ratio, shift, true, ctx.messages; strict = ctx.strict_guards)
  ratio, shift = _applyEndTaps(store, svsteps, et, ratio, shift, false, ctx.messages; strict = ctx.strict_guards)
  b1 = _ensureBus!(net, created, topo, svmap, i1.tn)
  b2 = _ensureBus!(net, created, topo, svmap, i2.tn)
  bf, bt = turned ? (b2, b1) : (b1, b2)
  # per-terminal flags (r0.9.10): a transformer with exactly one open
  # terminal stays in the model as its Schur reduction at the closed bus
  # (magnetizing branch seen through the ratio) instead of being dropped
  pt_in_service = _inService(ctx, pt)
  ptc_f = pt_in_service && _conn(ctx, i_f.connected)
  ptc_t = pt_in_service && _conn(ctx, i_t.connected)
  status = (ptc_f && ptc_t) ? 1 : 0
  (ptc_f != ptc_t) && push!(ctx.messages, "notice: PowerTransformer $(name) has one open terminal, kept as a partially open branch")
  from = Sparlectra.geNetBusIdx(net = net, busName = bf)
  to = Sparlectra.geNetBusIdx(net = net, busName = bt)
  ratedS = num(e1, :ratedS)
  Sparlectra._addPIModelTrafo_by_idx!(net = net, from = from, to = to, r_pu = r_pu, x_pu = x_pu, b_pu = b_pu, g_pu = g_pu, status = status, ratedU = Uf, ratedS = ratedS, ratio = ratio, shift_deg = shift, from_status = ptc_f ? 1 : 0, to_status = ptc_t ? 1 : 0, g_from_pu = g_from_pu, b_from_pu = b_from_pu, g_to_pu = g_to_pu, b_to_pu = b_to_pu)
  ctx.element_by_mrid[pt.mrid] = (:transformer, length(net.branchVec))
  idx = net.branchVec[end].branchIdx
  tmf = get(ef.refs, :Terminal, nothing)
  tmt = get(et.refs, :Terminal, nothing)
  tmf === nothing || (ctx.branch_side[tmf] = (idx, :from))
  tmt === nothing || (ctx.branch_side[tmt] = (idx, :to))
  # identity capture (transformer + ends + terminals): the PT counter runs
  # over every transformer-typed branch on the bus pair, in creation order.
  # The exporter writes the from side as end 1, so the end and terminal keys
  # follow the branch orientation: every mRID stays on its own bus
  k = cgmesNextBusEquipmentIndex!(ctx.eqp_counters, "PT", _cgmesTrafoPairKey(b1, b2))
  ptKey = cgmesKeyPowerTransformer(b1, b2, k)
  net.cgmes_ids[ptKey] = cgmesCanonicalMrid(pt.mrid)
  net.cgmes_ids[cgmesKeyTransformerEnd(ptKey, 1)] = cgmesCanonicalMrid(ef.mrid)
  net.cgmes_ids[cgmesKeyTransformerEnd(ptKey, 2)] = cgmesCanonicalMrid(et.mrid)
  tmf === nothing || (net.cgmes_ids[cgmesKeyTerminal(ptKey, 1)] = cgmesCanonicalMrid(tmf))
  tmt === nothing || (net.cgmes_ids[cgmesKeyTerminal(ptKey, 2)] = cgmesCanonicalMrid(tmt))
  # a phase tap changer keeps its mRID on the synthetic single-step changer
  # the export writes for the flattened shift
  ptc = something(_firstPhaseTapChangerMrid(store, ef.mrid), _firstPhaseTapChangerMrid(store, et.mrid), Some(nothing))
  ptc === nothing || (net.cgmes_ids[string(ptKey, "|PTC")] = cgmesCanonicalMrid(ptc))
  controlled = ctx.tap_control ? vcat(_attachTapControl!(net, store, topo, ctx, svsteps, ef, idx, true, name), _attachTapControl!(net, store, topo, ctx, svsteps, et, idx, false, name)) : Tuple{CIMObject,Bool}[]
  _attachTapRanges!(net, store, svsteps, ctx, idx, [(ef, true), (et, false)], controlled, name)
  _attachSETapChangers!(net, store, svsteps, ctx, idx, [(ef, true), (et, false)], ratio_struct, name)
end

function _map3WTrafo!(net, store, topo, created, svmap, svsteps, baseMVA, pt, ends, ctx::_MapCtx)
  name = something(str(pt, :name), pt.mrid)
  infos = [_terminalInfo(store, topo, e) for e in ends]
  if any(i -> i === nothing, infos)
    push!(ctx.messages, "skip: 3W PowerTransformer $(name), end terminal without TopologicalNode")
    push!(ctx.skipped, pt.mrid)
    return
  end
  e1 = ends[1]
  i1 = infos[1]
  U1 = something(num(e1, :ratedU), i1.vn_kV)
  vn_aux = i1.vn_kV

  # star (AUX) bus at the HV-side voltage base. The star point has no
  # SvVoltage of its own: it starts at the HV bus SV value here, and
  # _startStarBuses! replaces that with the value of its own node balance
  # once every branch exists (the HV value stays where a side has no SV)
  auxname = "AUX3WT_" * replace(name, r"\s+" => "_")
  if !Sparlectra.hasBusInNet(net = net, busName = auxname)
    vm0, va0 = haskey(svmap, i1.tn) ? (svmap[i1.tn][1] / i1.vn_kV, svmap[i1.tn][2]) : (1.0, 0.0)
    Sparlectra.addBus!(net = net, busName = auxname, vn_kV = vn_aux, vm_pu = vm0, va_deg = va0, isAux = true)
  end
  aux = Sparlectra.geNetBusIdx(net = net, busName = auxname)
  aux in ctx.star_buses || push!(ctx.star_buses, aux)

  # A leg whose side terminal is open stays in the model open at its side
  # bus (the per-terminal status of the 2W transformers): as long as another
  # end is connected the star is energized and the winding stays magnetized
  # from it, so the leg keeps its magnetizing arm at the star as the input
  # admittance of the open leg. A unit out of service or with every end open
  # has no live star; its legs are fully open.
  pt_in_service = _inService(ctx, pt)
  star_live = pt_in_service && any(i -> _conn(ctx, i.connected), infos)
  n_open = star_live ? count(i -> !_conn(ctx, i.connected), infos) : 0
  n_open > 0 && push!(ctx.messages, "notice: 3W PowerTransformer $(name) has $(n_open) open end terminal(s), the leg kept open at its side bus with the magnetizing admittance at the star")

  side_buses = String[]
  for (k, e) in enumerate(ends)
    ik = infos[k]
    Uk = something(num(e, :ratedU), ik.vn_kV)
    # end impedances and admittances are given on their own end base; the
    # tap table r/x/g/b of the end act on its leg
    fr, fx, fg, fb = _tapTableImpedanceFactors(store, svsteps, e)
    r_pu, x_pu, b_pu, g_pu = Sparlectra.toPU_RXBG(r = num(e, :r, 0.0) * fr, x = num(e, :x, 0.0) * fx, g = num(e, :g, 0.0) * fg, b = num(e, :b, 0.0) * fb, v_kv = ik.vn_kV, baseMVA = baseMVA)
    # structural ratio of the leg seen from the star (V_star / V_side at no load)
    t_struct = (U1 / Uk) / (vn_aux / ik.vn_kV)
    bk = _ensureBus!(net, created, topo, svmap, ik.tn)
    side = Sparlectra.geNetBusIdx(net = net, busName = bk)
    side_st = (pt_in_service && _conn(ctx, ik.connected)) ? 1 : 0
    star_st = star_live ? 1 : 0
    status = side_st * star_st
    # A leg with tap changers runs side bus → star: the tap sits at the
    # terminal of its end, the leg impedance on the inner side (the
    # arrangement the SV of real deliveries and OpenLoadFlow use, and the
    # one the PowSyBl builder uses for every leg). Ratio (structural and
    # tap) at the from side, impedance and magnetizing admittance on the
    # star base. A tap correction does not scale the impedance, so this is
    # not the same as the star → side leg with z·|c|². A leg without tap
    # changers keeps star → side with the impedance at the side bus, which
    # is the same circuit as the turned leg.
    #
    # The magnetizing admittance of every end sits at the star node, the
    # side arm is zero (0.42.0). The 3W terminal flows then close against
    # the SV at the star, not at the side terminal, and MicroGrid BE solves
    # closer to its SV. g and b are taken as delivered, sign included.
    turned = _trafoEndHasTaps(store, e)
    if turned
      ratio, shift = _applyEndTaps(store, svsteps, e, 1.0 / t_struct, 0.0, true, ctx.messages; strict = ctx.strict_guards)
      t2 = t_struct^2
      # the star is the to side: g, b referred to the star base through the
      # structural ratio only, as the to arm on the star bus
      Sparlectra._addPIModelTrafo_by_idx!(net = net, from = side, to = aux, r_pu = r_pu * t2, x_pu = x_pu * t2, b_pu = b_pu / t2, g_pu = g_pu / t2, status = status, ratedU = Uk, ratedS = num(e, :ratedS), ratio = ratio, shift_deg = shift, from_status = side_st, to_status = star_st, g_from_pu = 0.0, b_from_pu = 0.0, g_to_pu = g_pu / t2, b_to_pu = b_pu / t2)
    else
      # star → side bus: the end impedance on its own end base = the leg's
      # TO side, which is what calcAdmittance expects. The star is the from
      # side, so g, b on the end base are the from arm: it sits behind the
      # ratio and reaches the star bus as y / t_struct², the end admittance
      # referred to the star base
      Sparlectra._addPIModelTrafo_by_idx!(net = net, from = aux, to = side, r_pu = r_pu, x_pu = x_pu, b_pu = b_pu, g_pu = g_pu, status = status, ratedU = Uk, ratedS = num(e, :ratedS), ratio = t_struct, shift_deg = 0.0, from_status = star_st, to_status = side_st, g_from_pu = g_pu, b_from_pu = b_pu, g_to_pu = 0.0, b_to_pu = 0.0)
    end
    # one branch per leg, so the leg is addressed by its end, not by the
    # transformer
    ctx.element_by_mrid[e.mrid] = (:transformer, length(net.branchVec))
    tm = get(e.refs, :Terminal, nothing)
    tm === nothing || (ctx.branch_side[tm] = (net.branchVec[end].branchIdx, turned ? :from : :to))
    # legs also consume the pair counter so a hypothetical 2W transformer on
    # the same pair keeps a consistent <k>
    cgmesNextBusEquipmentIndex!(ctx.eqp_counters, "PT", _cgmesTrafoPairKey(auxname, bk))
    push!(side_buses, bk)
    # #294 point 4: the star-equivalent leg is an ordinary PI-model trafo
    # branch, so the same controller attachment and tap range as for 2W
    # transformers apply per leg; the tapped end is the from side of its leg
    legname = string(name, " (leg ", k, ")")
    legidx = net.branchVec[end].branchIdx
    controlled = ctx.tap_control ? _attachTapControl!(net, store, topo, ctx, svsteps, e, legidx, true, legname) : Tuple{CIMObject,Bool}[]
    _attachTapRanges!(net, store, svsteps, ctx, legidx, [(e, true)], controlled, legname)
    # a tapped leg runs side -> star from the structural ratio 1/t_struct
    # (the start of _applyEndTaps above), an untapped leg has no changer
    _attachSETapChangers!(net, store, svsteps, ctx, legidx, turned ? [(e, true)] : Tuple{CIMObject,Bool}[], turned ? 1.0 / t_struct : t_struct, legname)
  end
  # identity capture under the three-winding key (sorted side buses): the
  # exporter reassembles the star as one PowerTransformer with three ends
  # when the reconstruction is exact, and reuses these mRIDs
  if length(side_buses) == 3
    k3 = cgmesNextBusEquipmentIndex!(ctx.eqp_counters, "PT3", join(sort(side_buses), "|"))
    pt3Key = cgmesKeyPowerTransformer3W(side_buses[1], side_buses[2], side_buses[3], k3)
    net.cgmes_ids[pt3Key] = cgmesCanonicalMrid(pt.mrid)
    for (k, e) in enumerate(ends)
      net.cgmes_ids[cgmesKeyTransformerEnd(pt3Key, k)] = cgmesCanonicalMrid(e.mrid)
      tm = get(e.refs, :Terminal, nothing)
      tm === nothing || (net.cgmes_ids[cgmesKeyTerminal(pt3Key, k)] = cgmesCanonicalMrid(tm))
      ptc = _firstPhaseTapChangerMrid(store, e.mrid)
      ptc === nothing || (net.cgmes_ids[string(pt3Key, "|PTC", k)] = cgmesCanonicalMrid(ptc))
    end
  end
end

function _mapTransformers!(net, store, topo, created, svmap, svsteps, baseMVA, ctx::_MapCtx)
  for pt in objectsOf(store, :PowerTransformer)
    ends = _endsOfTransformer(store, pt)
    if length(ends) == 2
      _map2WTrafo!(net, store, topo, created, svmap, svsteps, baseMVA, pt, ends, ctx)
    elseif length(ends) == 3
      _map3WTrafo!(net, store, topo, created, svmap, svsteps, baseMVA, pt, ends, ctx)
    else
      push!(ctx.messages, "skip: PowerTransformer $(something(str(pt, :name), pt.mrid)) with $(length(ends)) ends")
    end
  end
  # the changers of mapped transformers whose delivery carries no step: one
  # notice each, or one for all of them (Sparlectra._tap_step_assumed_lines)
  items = Tuple{String,Float64}[]
  for pt in objectsOf(store, :PowerTransformer), e in _endsOfTransformer(store, pt)
    (haskey(ctx.element_by_mrid, pt.mrid) || haskey(ctx.element_by_mrid, e.mrid)) || continue
    for tc in _tapChangersOfEnd(store, e.mrid)
      _tapStepGiven(tc, svsteps) && continue
      push!(items, (string(something(str(pt, :name), pt.mrid), "/", something(str(tc, :name), tc.mrid)), _tapMiddleStep(tc)))
    end
  end
  for line in Sparlectra._tap_step_assumed_lines(items)
    push!(ctx.messages, "notice: " * line)
  end
end

"""
    _startStarBuses!(net, ctx, topo, svmap) -> Int

Start voltage of every three-winding star point from its own node balance.

The star point has no SvVoltage. Its provisional start, the SV voltage of
end 1, is off by the tap of that end and by the drop across the winding, so
at the SV state a current flows out of the star that the delivery does not
have, and Newton can end on a second solution at low voltage. With the start
voltages `V_k` of the sides, the star row of the admittance matrix gives the
one voltage at which no current leaves the star:

    V_s = -sum_k Y_sk V_k / Y_ss

The stamps are those `createYBUS` builds, summed over the branches at the
star (a leg with its side end open adds its input admittance to `Y_ss`
only). Applied when the star has at least one closed leg and every closed
leg reaches a side whose TopologicalNode carries an SV voltage above zero;
otherwise the provisional start stays. Runs once every branch is built, so
the stamps hold every ratio and impedance the mapping set, and `Y_ss` holds
the magnetizing admittance of the ends, which the legs carry as their arm at
the star. The importer puts no bus shunt on a star bus. Returns the number
of star buses set.
"""
function _startStarBuses!(net, ctx::_MapCtx, topo::CGMESTopology, svmap::Dict{String,Tuple{Float64,Float64}})::Int
  isempty(ctx.star_buses) && return 0
  # buses whose TopologicalNode carries a usable SV voltage (v <= 0 counts
  # as missing, as in compareWithSV)
  with_sv = Set{Int}()
  for (tn, (v, _)) in svmap
    (v > 0.0 && get(topo.vn_kV, tn, 0.0) > 0.0) || continue
    i = get(net.busDict, get(topo.bus_name, tn, ""), 0)
    i == 0 || push!(with_sv, i)
  end
  legs = Dict{Int,Vector{Int}}(s => Int[] for s in ctx.star_buses)
  for (j, br) in enumerate(net.branchVec)
    haskey(legs, br.fromBus) && push!(legs[br.fromBus], j)
    haskey(legs, br.toBus) && push!(legs[br.toBus], j)
  end
  volt(nd) = something(nd._vm_pu, 1.0) * cis(deg2rad(something(nd._va_deg, 0.0)))
  n_set = 0
  for s in ctx.star_buses
    y_ss = zero(ComplexF64)
    i_sides = zero(ComplexF64)   # sum_k Y_sk V_k over the closed legs
    n_closed = 0
    usable = true
    for j in legs[s]
      br = net.branchVec[j]
      state = Sparlectra._branch_terminal_state(br)
      state == :open && continue
      at_from = br.fromBus == s
      if state != :closed
        # one end open: only a closed star end sees the leg, as the reduced
        # input admittance on the diagonal
        (state == :open_to) == at_from && (y_ss += Sparlectra._open_terminal_yin(br))
        continue
      end
      side = at_from ? br.toBus : br.fromBus
      if !(side in with_sv)
        usable = false
        break
      end
      y11, y12, y21, y22 = Sparlectra.calcAdmittance(br, br.comp.cVN, net.baseMVA)
      y_ss += at_from ? y11 : y22
      i_sides += (at_from ? y12 : y21) * volt(net.nodeVec[side])
      n_closed += 1
    end
    (usable && n_closed > 0 && abs(y_ss) > 0.0) || continue
    vs = -i_sides / y_ss
    Sparlectra.setVmVa!(node = net.nodeVec[s], vm_pu = abs(vs), va_deg = rad2deg(angle(vs)))
    n_set += 1
  end
  return n_set
end

"""
    _startNoSvBusesFromNeighbours!(net, ctx, no_sv_buses) -> Int

Start voltage of the buses without a usable SvVoltage from a neighbour.

Such a bus gets 1.0 pu / 0° from `_ensureBus!`. Where its neighbours sit
far from 0°, that start can lead Newton to a second solution: a leaf bus
without load behind a branch converges at 0 V, a short circuit behind the
branch impedance that pulls its neighbour down (seen on a real delivery
with neighbours at -22° to -26°). Each listed bus therefore takes magnitude
and angle of the neighbour reached over the closed branch or closed link
of the smallest |z| (a link counts as zero impedance, a tie goes to the
lower bus index). The search runs breadth-first from the buses with a
start of their own, so a chain of such buses is covered; a bus picks only
from buses set before its layer, so the result does not depend on the
order within a layer. A bus whose group
touches no bus with a start keeps the default. Three-winding star points
are left alone (`_startStarBuses!` sets them by their node balance) and
are no source either; isolated buses carry no start and are neither.
Returns the number of buses set.
"""
function _startNoSvBusesFromNeighbours!(net, ctx::_MapCtx, no_sv_buses::Vector{String})::Int
  isempty(no_sv_buses) && return 0
  n = length(net.nodeVec)
  blocked = falses(n)   # star points and isolated buses: no source
  for s in ctx.star_buses
    blocked[s] = true
  end
  # a listed star point is skipped like an isolated bus: set from a
  # neighbour it would still be no source, and a listed bus behind it would
  # enter the next layer without a source (index 0)
  pending = falses(n)
  for b in no_sv_buses
    i = get(net.busDict, b, 0)
    (i == 0 || blocked[i] || Sparlectra.isIsolated(net.nodeVec[i])) || (pending[i] = true)
  end
  any(pending) || return 0
  for i in 1:n
    Sparlectra.isIsolated(net.nodeVec[i]) && (blocked[i] = true)
  end
  # (neighbour, |z|) over closed branches and closed links
  adj = [Tuple{Int,Float64}[] for _ in 1:n]
  for br in net.branchVec
    Sparlectra._branch_terminal_state(br) == :closed || continue
    z = abs(complex(br.r_pu, br.x_pu))
    push!(adj[br.fromBus], (br.toBus, z))
    push!(adj[br.toBus], (br.fromBus, z))
  end
  for l in net.linkVec
    l.status == 1 || continue
    push!(adj[l.fromBus], (l.toBus, 0.0))
    push!(adj[l.toBus], (l.fromBus, 0.0))
  end
  known(j) = !pending[j] && !blocked[j]
  layer = [i for i in 1:n if pending[i] && any(e -> known(e[1]), adj[i])]
  n_set = 0
  while !isempty(layer)
    picks = Tuple{Int,Int}[]
    for i in layer
      # smallest |z|, a tie goes to the lower bus index
      best, best_z = 0, Inf
      for (j, z) in adj[i]
        if known(j) && (z < best_z || (z == best_z && j < best))
          best, best_z = j, z
        end
      end
      push!(picks, (i, best))
    end
    for (i, j) in picks
      src = net.nodeVec[j]
      Sparlectra.setVmVa!(node = net.nodeVec[i], vm_pu = something(src._vm_pu, 1.0), va_deg = something(src._va_deg, 0.0))
      pending[i] = false
    end
    n_set += length(picks)
    next = Set{Int}()
    for (i, _) in picks, (j, _) in adj[i]
      pending[j] && push!(next, j)
    end
    layer = sort!(collect(next))
  end
  return n_set
end

function _mapShunts!(net, store, topo, created, svmap, ctx::_MapCtx)
  # shared service/connection gate for both shunt classes; returns the
  # terminal info or `nothing` after pushing the skip message
  function _shunt_terminal(sh, clsname, name)
    t = busOfEquipment(topo, sh, 1)
    if t === nothing
      push!(ctx.messages, "skip: $(clsname) $(name), no TopologicalNode")
      push!(ctx.skipped, sh.mrid)
      return nothing
    end
    if !_inService(ctx, sh)
      push!(ctx.messages, "skip: $(clsname) $(name), out of service (SSH inService=false)")
      push!(ctx.skipped, sh.mrid)
      return nothing
    end
    if !_conn(ctx, t.connected)
      push!(ctx.messages, "skip: $(clsname) $(name), terminal disconnected")
      push!(ctx.skipped, sh.mrid)
      return nothing
    end
    return t
  end

  # Placeholder guard for shunt admittances, same philosophy as the vset
  # plausibility band: FullGrid ships a NonlinearShuntCompensatorPoint with
  # b = g = 0.99 S (the set's 0.99/99.99 placeholder family), at 225 kV that
  # is a 50-GW shunt and no power flow can converge on it. The cap of
  # 10 × baseMVA is the scale the SVC rating clamp already uses; every real
  # compensator in the cached deliveries sits orders of magnitude below it.
  # A shunt above the cap is skipped with a warning, not clamped: clamping
  # would invent a specific size for a value that carries no information.
  shunt_cap_mva = 10.0 * net.baseMVA
  function _shunt_plausible(clsname, name, mrid, B_mvar, G_mw)::Bool
    max(abs(B_mvar), abs(G_mw)) <= shunt_cap_mva && return true
    msg = "$(clsname) $(name) declares an implausible admittance (B $(round(B_mvar; digits = 1)) MVAr / G $(round(G_mw; digits = 1)) MW at nominal voltage, cap ±$(round(Int, shunt_cap_mva))), placeholder suspected"
    ctx.strict_guards && error("cgmes_import.placeholder_guards = strict: " * msg)
    push!(ctx.messages, "warning: " * msg * ", shunt skipped")
    push!(ctx.skipped, mrid)
    return false
  end

  for sh in objectsOf(store, :LinearShuntCompensator)
    name = something(str(sh, :name), sh.mrid)
    t = _shunt_terminal(sh, "LinearShuntCompensator", name)
    t === nothing && continue
    sections = something(num(sh, :sections), num(sh, :normalSections, 0.0))
    b = num(sh, :bPerSection, 0.0) * sections
    g = num(sh, :gPerSection, 0.0) * sections
    _shunt_plausible("LinearShuntCompensator", name, sh.mrid, b * t.vn_kV^2, g * t.vn_kV^2) || continue
    bus = _ensureBus!(net, created, topo, svmap, t.tn)
    # MATPOWER-style shunt: MW/MVAr injected at 1.0 pu of the bus vn
    Sparlectra.addShuntMatpower!(net = net, busName = bus, Gs = g * t.vn_kV^2, Bs = b * t.vn_kV^2)
    ctx.element_by_mrid[sh.mrid] = (:shunt, length(net.shuntVec))
    _captureBusEquipment!(net, ctx, "SH", bus, sh.mrid; children = (T1 = _firstTerminalMrid(topo, sh),))
  end

  # NonlinearShuntCompensator (#294 point 7): the characteristic is a set of
  # per-section points; the effective admittance at `sections = s` is the SUM
  # of the point b/g values with sectionNumber <= s (each point is one
  # switched group, same reading as PowSyBl's conversion).
  for sh in objectsOf(store, :NonlinearShuntCompensator)
    name = something(str(sh, :name), sh.mrid)
    t = _shunt_terminal(sh, "NonlinearShuntCompensator", name)
    t === nothing && continue
    sections = something(num(sh, :sections), num(sh, :normalSections, 0.0))
    pts = [pt for pt in objectsOf(store, :NonlinearShuntCompensatorPoint) if get(pt.refs, :NonlinearShuntCompensator, "") == sh.mrid]
    if isempty(pts)
      push!(ctx.messages, "skip: NonlinearShuntCompensator $(name), no characteristic points")
      push!(ctx.skipped, sh.mrid)
      continue
    end
    active = [pt for pt in pts if something(num(pt, :sectionNumber), 1.0) <= sections]
    b = sum((num(pt, :b, 0.0) for pt in active); init = 0.0)
    g = sum((num(pt, :g, 0.0) for pt in active); init = 0.0)
    _shunt_plausible("NonlinearShuntCompensator", name, sh.mrid, b * t.vn_kV^2, g * t.vn_kV^2) || continue
    bus = _ensureBus!(net, created, topo, svmap, t.tn)
    Sparlectra.addShuntMatpower!(net = net, busName = bus, Gs = g * t.vn_kV^2, Bs = b * t.vn_kV^2)
    ctx.element_by_mrid[sh.mrid] = (:shunt, length(net.shuntVec))
    # the aggregated characteristic re-exports as a LinearShuntCompensator,
    # keeping the original mRID
    _captureBusEquipment!(net, ctx, "SH", bus, sh.mrid; children = (T1 = _firstTerminalMrid(topo, sh),))
    push!(ctx.messages, "nonlinear shunt $(name): sections $(Int(round(sections)))/$(length(pts)) → B $(round(b * t.vn_kV^2; digits = 3)) MVAr, G $(round(g * t.vn_kV^2; digits = 3)) MW at $(bus)")
  end
end

# Bank regulators of the delivery (cgmes_import.shunt_control): a
# LinearShuntCompensator with SSH controlEnabled and an enabled
# voltage RegulatingControl at its own TopologicalNode becomes a switched
# bank (ShuntVoltageControl in discrete mode) that takes over the mapped
# shunt (addShuntVoltageControl! with shunt_idx): blocks of bPerSection at
# the nominal voltage, from 0 to maximumSections, the SSH sections as start,
# the target in pu of the regulated bus, the deadband half of
# targetDeadband (as for the tap controllers), at least 1e-3 pu. Each case
# the controller cannot represent stays a fixed shunt with a notice: a
# remote target (ShuntVoltageControl holds its own bus), a nonlinear
# characteristic (no uniform block), a conductance, a target outside the
# vset band, a bus that is not PQ, a second bank on one bus. Runs after
# the bus types are final.
function _attachShuntControls!(net, store::CGMESStore, topo::CGMESTopology, ctx::_MapCtx)
  n_ctrl = 0
  n_kept = 0
  keep(name, why) = (push!(ctx.messages, "notice: shunt control: $(name) stays a fixed shunt, $(why)"); n_kept += 1)
  for sh in objectsOf(store, :NonlinearShuntCompensator)
    haskey(ctx.element_by_mrid, sh.mrid) || continue
    something(boolval(sh, :controlEnabled), false) || continue
    ref(store, sh, :RegulatingControl) === nothing && continue
    keep(something(str(sh, :name), sh.mrid), "a nonlinear characteristic has no uniform block")
  end
  for sh in objectsOf(store, :LinearShuntCompensator)
    el = get(ctx.element_by_mrid, sh.mrid, nothing)
    (el === nothing || el[1] !== :shunt) && continue
    something(boolval(sh, :controlEnabled), false) || continue
    rc = ref(store, sh, :RegulatingControl)
    rc === nothing && continue
    name = something(str(sh, :name), sh.mrid)
    endswith(something(enumval(rc, :mode), ""), ".voltage") || (keep(name, "its RegulatingControl is no voltage control"); continue)
    something(boolval(rc, :enabled), true) || continue
    t = busOfEquipment(topo, sh, 1)
    t === nothing && continue
    own = topo.bus_name[t.tn]
    tinfo = _regulatedBus(store, topo, rc)
    tinfo === nothing && (keep(name, "control terminal unresolved"); continue)
    tinfo.bus == own || (keep(name, "remote target bus $(tinfo.bus) (a bank regulator holds its own bus)"); continue)
    num(sh, :gPerSection, 0.0) == 0.0 || (keep(name, "it has a conductance per section"); continue)
    target = num(rc, :targetValue, 0.0) / tinfo.vn_kV
    (ctx.vset_min_pu <= target <= ctx.vset_max_pu) || (keep(name, "voltage target $(round(target; digits = 4)) pu outside $(ctx.vset_min_pu) to $(ctx.vset_max_pu)"); continue)
    step = num(sh, :bPerSection, 0.0) * t.vn_kV^2
    nmax = num(sh, :maximumSections, 0.0)
    (step != 0.0 && nmax >= 1.0) || (keep(name, "no block (bPerSection or maximumSections)"); continue)
    busidx = Sparlectra.geNetBusIdx(net = net, busName = own)
    Sparlectra.getNodeType(net.nodeVec[busidx]) == Sparlectra.PQ || (keep(name, "bus $(own) is $(Sparlectra.getNodeType(net.nodeVec[busidx])), its voltage is held otherwise"); continue)
    any(c -> c isa Sparlectra.ShuntVoltageControl && c.bus == own, net.machineControls) && (keep(name, "bus $(own) has a bank regulator already"); continue)
    lo, hi = minmax(0.0, nmax * step)
    deadband = max(num(rc, :targetDeadband, 0.0) / tinfo.vn_kV / 2.0, 1e-3)
    Sparlectra.addShuntVoltageControl!(net; bus = own, target_vm_pu = target, bs_min_mvar = lo, bs_max_mvar = hi, step_mvar = abs(step), deadband_vm_pu = deadband,
                                       name = name, shunt_idx = el[2])
    n_ctrl += 1
    push!(ctx.messages, "shunt control: $(name) → voltage $(round(target; digits = 4)) pu at $(own), $(Int(nmax)) blocks of $(round(abs(step); digits = 3)) MVAr")
  end
  (n_ctrl + n_kept) > 0 && push!(ctx.messages, "notice: cgmes_import.shunt_control = true: $(n_ctrl) bank regulator(s) attached, $(n_kept) controlled bank(s) kept as fixed shunts")
  return nothing
end

function _mapLoads!(net, store, topo, created, svmap, ctx::_MapCtx)
  for cls in (:EnergyConsumer, :ConformLoad, :NonConformLoad, :StationSupply)
    for load in objectsOf(store, cls)
      t = busOfEquipment(topo, load, 1)
      name = something(str(load, :name), load.mrid)
      if t === nothing
        push!(ctx.messages, "skip: $(cls) $(name), no TopologicalNode")
        push!(ctx.skipped, load.mrid)
        continue
      end
      if !_inService(ctx, load)
        push!(ctx.messages, "skip: $(cls) $(name), out of service (SSH inService=false)")
        push!(ctx.skipped, load.mrid)
        continue
      end
      if !_conn(ctx, t.connected)
        push!(ctx.messages, "skip: $(cls) $(name), terminal disconnected")
        push!(ctx.skipped, load.mrid)
        continue
      end
      p = num(load, :p, 0.0)
      q = num(load, :q, 0.0)
      bus = _ensureBus!(net, created, topo, svmap, t.tn)
      Sparlectra.addProsumer!(net = net, busName = bus, type = "ENERGYCONSUMER", p = p, q = q, defer_bus_type_refresh = true)
      ctx.element_by_mrid[load.mrid] = (:load, length(net.prosumpsVec))
      # identity capture: all load classes re-export as EnergyConsumer, the
      # original mRID is what survives the roundtrip
      _captureBusEquipment!(net, ctx, "EC", bus, load.mrid; children = (T1 = _firstTerminalMrid(topo, load),))
    end
  end
end

# regulating-control voltage setpoint for a machine/injection.
# Returns `(vm_pu, voltage_controlled, remote_bus)`:
# - local target → (target/vn, true, nothing)
# - remote target with allow_remote (machine_control) → (target/vn_remote,
#   true, remote bus name), the caller wires an outer-loop controller
# - remote target without allow_remote → (nothing, true, nothing);
#   the caller keeps the unit voltage-controlled at its own bus, holding the
#   bus's SV/start voltage instead of a stale fixed SSH q
# - no active voltage control → (nothing, false, nothing)
# Default plausibility band for a voltage RegulatingControl target, in p.u. of
# the regulated bus's nominal voltage. Deliberately far wider than any real
# setpoint: across every cached ENTSO-E/ReliCapGrid delivery (1187 targets) the
# actual range is 0.92 .. 1.15 p.u. The band exists to catch placeholder values,
# not to enforce an operating range.
#
# Svedala ships five generator nodes with targetValue = 0.001 kV (≈ 5e-5 p.u.):
# the same five that carry no SvVoltage, i.e. units the exporting tool did not
# include in its solved state. Taken literally they produce a PV/reference bus
# at essentially zero volts, and the power flow then "converges" onto a
# collapsed all-zero solution that is formally correct and physically nonsense.
#
# These are DEFAULTS, not fixed law: the numbers come from the deliveries we
# have seen, so they are overridable per import (`vset_min_pu`/`vset_max_pu`,
# also as `cgmes_import.*` config keys) instead of silently imposing a policy on
# a delivery that legitimately regulates outside the band.
const CGMES_VSET_MIN_PU = 0.5
const CGMES_VSET_MAX_PU = 1.5

"""
Voltage setpoint of a regulating unit, as
`(vset_pu, is_voltage_regulating, remote_bus)`.

`vset_pu === nothing` together with `true` means "regulates voltage, but no
usable local setpoint", the caller then holds the unit PV at its own bus
voltage (see `_pvVoltage`). That is the fallback for remote controls (unless
`allow_remote` is set) and for implausible target values.

With `allow_remote = true` (machine control, #294 point 3) a resolvable remote
voltage control returns the target bus name in `remote_bus` and `vset_pu` in
p.u. of the *remote* bus's nominal voltage; the same plausibility band applies
there. `remote_bus === nothing` always means "regulate locally or not at all".
"""
function _voltageSetpoint(
  store::CGMESStore,
  topo::CGMESTopology,
  obj::CIMObject,
  own_tn::String,
  messages::Vector{String};
  vset_min_pu::Float64 = CGMES_VSET_MIN_PU,
  vset_max_pu::Float64 = CGMES_VSET_MAX_PU,
  allow_remote::Bool = false,
)::Tuple{Union{Nothing,Float64},Bool,Union{Nothing,String}}
  rc = ref(store, obj, :RegulatingControl)
  rc === nothing && return (nothing, false, nothing)
  mode = something(enumval(rc, :mode), "")
  endswith(mode, ".voltage") || return (nothing, false, nothing)
  something(boolval(rc, :enabled), true) || return (nothing, false, nothing)
  target = num(rc, :targetValue)
  (target === nothing || target <= 0.0) && return (nothing, false, nothing)
  tmrid = get(rc.refs, :Terminal, nothing)
  tmrid === nothing && return (nothing, false, nothing)
  t = get(store.objects, tmrid, nothing)
  t === nothing && return (nothing, false, nothing)
  tn = tnOfTerminal(t)
  (tn === nothing || !haskey(topo.vn_kV, tn)) && return (nothing, false, nothing)
  is_remote = tn != own_tn
  if is_remote && !allow_remote
    push!(messages, "notice: $(something(str(obj, :name), obj.mrid)) has a remote voltage RegulatingControl (target bus $(get(topo.bus_name, tn, tn))), held PV at its own bus (enable machine_control for remote regulation)")
    return (nothing, true, nothing)
  end
  vn = topo.vn_kV[tn]
  # A zero/absent nominal voltage makes the p.u. conversion meaningless.
  vn > 0.0 || return (nothing, true, nothing)
  vset = target / vn
  if vset < vset_min_pu || vset > vset_max_pu
    push!(
      messages,
      "warning: $(something(str(obj, :name), obj.mrid)) declares an implausible voltage target ($(target) kV = $(round(vset; sigdigits = 3)) pu at $(get(topo.bus_name, tn, tn)), nominal $(vn) kV), ignored, unit held PV at the bus voltage derived from the nominal data",
    )
    return (nothing, true, nothing)
  end
  return (vset, true, is_remote ? get(topo.bus_name, tn, tn) : nothing)
end

# limit hull over both sign-convention readings (see the machine qMin/qMax
# comment in _mapInjections!): [min(a,-b), max(b,-a)]
_limitHullMin(a::Union{Nothing,Float64}, b::Union{Nothing,Float64}) = (a === nothing || b === nothing) ? nothing : min(a, -b)
_limitHullMax(a::Union{Nothing,Float64}, b::Union{Nothing,Float64}) = (a === nothing || b === nothing) ? nothing : max(b, -a)

# An active power limit of 9999 MW and above in magnitude is the way a
# delivery says "no limit stated" (a placeholder such as +-99999.9 MW on
# units without a limit, MATPOWER and IIDM cases 9999 on every unit); no unit is that
# large. Taken as a limit, such a unit wins the slack ranking by size and the
# island slack choice by reserve, and takes almost all of the pmax and
# headroom weights of the distributed slack. The importer reads it as no
# limit (nothing). One constant with the short circuit, which reads a
# maximum of that size as no machine size: Sparlectra.SC_PMAX_PLACEHOLDER_MW.
# The core keeps such values, because MATPOWER cases carry 9999 as their
# ordinary Pmax and the readers there must not change.
_pLimit(x::Union{Nothing,Float64}) = (x === nothing || abs(x) >= Sparlectra.SC_PMAX_PLACEHOLDER_MW) ? nothing : x

# Both P limits of one mapped unit through _pLimit; a unit with a
# placeholder on either side is counted once under `class` for the notice.
function _pLimitsRead!(ctx::_MapCtx, class::String, lo::Union{Nothing,Float64}, hi::Union{Nothing,Float64})
  plo = _pLimit(lo)
  phi = _pLimit(hi)
  ((plo === nothing && lo !== nothing) || (phi === nothing && hi !== nothing)) && (ctx.p_limit_placeholders[class] = get(ctx.p_limit_placeholders, class, 0) + 1)
  return plo, phi
end

"""
Collect `ReactiveCapabilityCurve` points: curve mRID → `CurveData` tuples
`(x = P, y1 = minQ, y2 = maxQ)`, sorted by P. The x axis is the machine's own
active power in the CGMES machine convention (negative = injecting), so it can
legitimately span both signs (MicroGrid BE-G1: −100 … +100 MW).
"""
function _reactiveCapabilityCurvePoints(store::CGMESStore)::Dict{String,Vector{NTuple{3,Float64}}}
  out = Dict{String,Vector{NTuple{3,Float64}}}()
  countOf(store, :ReactiveCapabilityCurve) == 0 && return out
  curves = Set(c.mrid for c in objectsOf(store, :ReactiveCapabilityCurve))
  for cd in objectsOf(store, :CurveData)
    c = get(cd.refs, :Curve, nothing)
    (c === nothing || !(c in curves)) && continue
    x = num(cd, :xvalue)
    y1 = num(cd, :y1value)
    y2 = num(cd, :y2value)
    (x === nothing || y1 === nothing || y2 === nothing) && continue
    push!(get!(Vector{NTuple{3,Float64}}, out, c), (x, y1, y2))
  end
  for pts in values(out)
    sort!(pts; by = first)
  end
  return out
end

"""
Q limits from a capability curve, evaluated at the machine's scheduled active
power `p_machine` (CGMES machine convention, i.e. the raw SSH value): linear
interpolation between curve points, clamped to the curve's P domain outside it.

The interpolated `(y1, y2)` pair then goes through the SAME sign-convention
hull as the scalar `minQ`/`maxQ` path, the ENTSO-E sets are inconsistent
about machine Q signs, and a capability curve is no more trustworthy in that
respect than the scalars next to it. For the symmetric curves the test sets
ship this is the identity; for asymmetric data it is deliberately safe rather
than strictly faithful (the same Stage-2 note as for the scalars applies).
"""
function _curveQHull(pts::Vector{NTuple{3,Float64}}, p_machine::Float64)::Union{Nothing,Tuple{Float64,Float64}}
  isempty(pts) && return nothing
  p = clamp(p_machine, pts[1][1], pts[end][1])
  i = something(findlast(t -> t[1] <= p, pts), 1)
  y1, y2 = if i >= length(pts)
    (pts[end][2], pts[end][3])
  else
    t0, t1 = pts[i], pts[i+1]
    f = t1[1] == t0[1] ? 0.0 : (p - t0[1]) / (t1[1] - t0[1])
    (t0[2] + f * (t1[2] - t0[2]), t0[3] + f * (t1[3] - t0[3]))
  end
  # order the pair first: the hull formula assumes (min, max) roles, and a
  # delivery that swaps y1/y2 wholesale would otherwise produce an inverted,
  # rejected band instead of the intended one
  ylo, yhi = minmax(y1, y2)
  qlo = _limitHullMin(ylo, yhi)
  qhi = _limitHullMax(ylo, yhi)
  (qlo === nothing || qhi === nothing || (qhi - qlo) < 1e-6) && return nothing
  return (qlo, qhi)
end

# PV voltage for a unit: local target if given, else the SV/start voltage of
# its own bus (remote-controlled units, slack fallback)
function _pvVoltage(net, bus::String, vset::Union{Nothing,Float64})::Float64
  vset !== nothing && return vset
  node = net.nodeVec[Sparlectra.geNetBusIdx(net = net, busName = bus)]
  return something(node._vm_pu, 1.0)
end

# The referencePriority of a machine or network injection as the prosumer
# keeps it: CIM states an integer, 0 means "no preference" and so does any
# value below 1 (the importer reads numbers as Float64, a negative or
# fractional value below 1 is no priority in the CGMES sense either)
_cgmesReferencePriority(prio::Float64)::Int = (isfinite(prio) && prio >= 1.0) ? round(Int, prio) : 0

# slack candidate bookkeeping for the slack selection
#
# The tuples are laid out so that plain ascending `sort!` implements the ranking.
# SV angle references: class rank first (0 for a machine or a network
# injection, 1 for an EquivalentInjection or a source made of one, so the
# order among the units the importer took before stays as it was), then
# -rating so the larger unit comes first, then the bus name as a
# deterministic last resort. referencePriority markers: priority ascending
# (lower referencePriority = stronger), then -rating, then bus name.
mutable struct _SlackScan
  sv_refs::Vector{Tuple{Int,Float64,String}}           # (class rank, -rating, bus)
  by_priority::Vector{Tuple{Float64,Float64,String}}   # (priority, -rating, bus)
  enis::Vector{Tuple{Float64,String}}                  # (maxP/rating, bus)
  machines::Vector{Tuple{Float64,String,String}}       # (ratedS, bus, machine name)
end
_SlackScan() = _SlackScan([], [], [], [])

"""
Pick the primary slack bus from the collected candidates.

Fallback chain, strongest evidence first:
 1. the SV angle reference (`TopologicalIsland.AngleRefTopologicalNode`) on
    a node that carries a `SynchronousMachine`, an
    `ExternalNetworkInjection` or an `EquivalentInjection`: it is the
    reference the exporting tool actually used;
 2. SSH `referencePriority` markers;
 3. the only / the largest `ExternalNetworkInjection`: a network equivalent is
    the natural reference;
 4. the largest `SynchronousMachine`.

`island_label(bus)` names the island of a bus for the messages (the
`TopologicalIsland` name from the SV, else the electrical island number);
`has_sv_ref` says whether the SV declares any angle reference at all, so the
message of step 4 does not claim there is none when the declared one
carries no injection.

Errors out when a delivery offers none of them: an angle reference cannot
be invented.
"""
function _chooseSlack(scan::_SlackScan, messages::Vector{String}, island_label::Function, has_sv_ref::Bool)::String
  if !isempty(scan.sv_refs)
    # Ascending sort ranks by the tuple layout documented at _SlackScan.
    sort!(scan.sv_refs)
    bus = scan.sv_refs[1][3]
    push!(messages, "notice: $(island_label(bus)): angle reference $(bus) taken from SV")
    return bus
  end
  if !isempty(scan.by_priority)
    sort!(scan.by_priority)
    bus = scan.by_priority[1][3]
    top = scan.by_priority[1][1]
    ties = count(t -> t[1] == top, scan.by_priority)
    if ties > 1
      # CGMES leaves the winner undefined when several units share the
      # highest referencePriority; we break the tie by rating, then by name,
      # which need not match the exporting tool's choice (RealGrid: two
      # machines at priority 1, its documentation names the other one).
      push!(messages, "notice: $(ties) units share referencePriority $(top), slack chosen by largest rating; the source tool may have picked a different one")
    end
    push!(messages, "slack: referencePriority marker → $(bus)")
    return bus
  end
  # An ExternalNetworkInjection represents the rest of the grid, so it is the
  # natural reference when no explicit priority is given. The single-ENI case is
  # separated only to produce a clearer message.
  if length(scan.enis) == 1
    push!(messages, "slack: single ExternalNetworkInjection → $(scan.enis[1][2])")
    return scan.enis[1][2]
  end
  if length(scan.enis) > 1
    sort!(scan.enis; by = t -> -t[1])
    push!(messages, "slack: largest ExternalNetworkInjection → $(scan.enis[1][2])")
    return scan.enis[1][2]
  end
  isempty(scan.machines) && error("CGMES import: no slack candidate found (no SV angle reference with an injection, no referencePriority, no ExternalNetworkInjection, no SynchronousMachine)")
  sort!(scan.machines; by = t -> -t[1])
  _, bus, mname = scan.machines[1]
  if has_sv_ref
    push!(messages, "slack: largest SynchronousMachine → $(bus)")
  else
    push!(messages, "notice: $(island_label(bus)): no angle reference in SV or SSH, largest machine $(mname) at $(bus) as slack")
  end
  return bus
end

# A pending tuple that can hold the angle reference: a SynchronousMachine,
# an ExternalNetworkInjection, or an EquivalentInjection (also one imported
# as a source). The type string of the last two is EXTERNALNETWORKINJECTION.
# An HVDC converter shares that type string but is a fixed PCC injection
# (`hvdc_converter = true`); EnergySource and PowerElectronicsConnection
# tuples carry `slack_candidate = false`; a StaticVarCompensator has no
# active power. None of these can be the reference.
_isReferenceUnit(inj)::Bool = inj.type in ("SYNCHRONOUSMACHINE", "EXTERNALNETWORKINJECTION") && get(inj, :slack_candidate, true) && !get(inj, :hvdc_converter, false)

"""
Decide the slack bus set. With `multi_slack` every **electrical** island (in
the sense of `electricalIslandComponents`, i.e. closed links count as
connections) receives at most one reference: the island of the primary
slack (`_chooseSlack`) keeps it, and every other island takes its SV angle
reference. A further SV angle reference in an island that already has one
is not promoted, and every injection other than the reference keeps the
bus type its own control gives it (PV for a regulating unit, PQ
otherwise). Two references in one island would be physically wrong. On
clean data the guard is a no-op (Svedala's two references sit in two
electrical islands once `Equipment.inService` is honored), but it protects
against deliveries whose island declarations and switch states disagree.

`island_label(bus)` names the island of a bus for the messages.

An island that has a source but neither an angle reference here nor a
voltage-regulating unit gets its slack later, from
`_slackForIslandsWithoutReference!`, once de-energization and isolation have
decided which islands exist.
"""
function _selectSlackBuses(net, ctx::_MapCtx, primary::String, angle_ref_buses::Set{String}, pending, multi_slack::Bool, island_label::Function)::Set{String}
  chosen = Set{String}([primary])
  # Only angle references that carry a unit able to hold the reference are
  # candidates (_isReferenceUnit): a reference node without one cannot hold
  # a voltage and would produce a slack bus with nothing behind it.
  candidates = [b for b in angle_ref_buses if b != primary && any(inj -> inj.bus == b && _isReferenceUnit(inj), pending)]
  isempty(candidates) && return chosen
  if !multi_slack
    push!(ctx.messages, "notice: the SV declares further angle reference(s) ($(join(sort(candidates), ", "))), not promoted to slack; pass multi_slack = true to give separate electrical islands their own reference")
    return chosen
  end

  # One reference per electrical island. The electrical view is essential here:
  # counting only branches, retained switches would split one island into
  # several and we would hand out a second slack inside a single island.
  island_of = Sparlectra.electricalIslandOfBus(net)
  # electrical island -> the bus holding its reference
  holder = Dict{Int,String}()
  pi = get(island_of, get(net.busDict, primary, 0), 0)
  # Island 0 means "bus unknown/isolated", nothing to reserve in that case.
  pi == 0 || (holder[pi] = primary)
  # Sorted so the outcome does not depend on Set iteration order.
  for b in sort(candidates)
    isl = get(island_of, get(net.busDict, b, 0), 0)
    if isl == 0 || haskey(holder, isl)
      # the line names the holder's island; b may come from another SV
      # TopologicalIsland, then its own island is named as well
      h = get(holder, isl, primary)
      h_isl, b_isl = island_label(h), island_label(b)
      b_of = b_isl == h_isl ? "" : " of $(b_isl)"
      push!(ctx.messages, "notice: $(h_isl): angle reference $(h), other injections keep their own bus type; the SV angle reference $(b)$(b_of) lies in the same electrical island and is not promoted to slack")
      continue
    end
    holder[isl] = b
    push!(chosen, b)
    push!(ctx.messages, "notice: $(island_label(b)): angle reference $(b) taken from SV")
    push!(ctx.messages, "slack: additional reference $(b) for electrical island $(isl)")
  end
  return chosen
end

# The slack unit of an island without reference, among the prosumer indices
# `ks` of its units; `degree` is the number of neighbour buses per bus.
# Returns (index, reason). Candidates are only the units with infeed
# (pVal > 0, generator sign); the caller passes only islands that have one.
# A unit with P = 0 or drawing power is never the slack, also not by its
# referencePriority or its reserve, because it feeds nothing into the
# island. Among the candidates:
#  1. the strongest referencePriority (the smallest positive value), ties to
#     the smallest prosumer index;
#  2. else the largest active power reserve maxP - pVal among the units with
#     a known maxP on buses whose degree is not 1;
#  3. else the largest infeed on a bus whose degree is not 1, and when only
#     degree-1 buses carry infeed (all buses of a two-bus fragment have
#     degree 1), the largest infeed on a degree-1 bus.
# In an island of two or more buses "degree not 1" means degree greater than
# 1; the only bus of a single-bus island has degree 0.
function _islandSlackUnit(net, ks::Vector{Int}, degree::Vector{Int})::Tuple{Int,String}
  pval(k) = something(net.prosumpsVec[k].pVal, 0.0)
  inner(k) = degree[Int(net.prosumpsVec[k].comp.cFrom_bus)] != 1
  feed = [k for k in ks if pval(k) > 0.0]
  prio = [k for k in feed if net.prosumpsVec[k].referencePriority > 0]
  if !isempty(prio)
    k = first(sort(prio; by = k -> (net.prosumpsVec[k].referencePriority, k)))
    return k, "referencePriority $(net.prosumpsVec[k].referencePriority)"
  end
  reserve(k) = net.prosumpsVec[k].maxP - pval(k)
  with_reserve = [k for k in feed if inner(k) && net.prosumpsVec[k].maxP !== nothing && isfinite(net.prosumpsVec[k].maxP)]
  if !isempty(with_reserve)
    k = first(sort(with_reserve; by = k -> (-reserve(k), k)))
    return k, "largest active power reserve off the degree-1 buses, $(round(reserve(k); digits = 2)) MW"
  end
  pool = filter(inner, feed)
  note = isempty(pool) ? "only degree-1 buses carry infeed" : "no known maxP off the degree-1 buses"
  isempty(pool) && (pool = feed)
  k = first(sort(pool; by = k -> (-pval(k), k)))
  return k, "largest infeed, P = $(round(pval(k); digits = 2)) MW, $(note)"
end

"""
Give every energized electrical island that has infeed but no reference a
slack.

Without this, such an island leaves the import without a Slack or PV bus. The
island-wise solve then promotes one of its units to the reference and keeps
that bus at its start voltage, so the result there depends on the start
value only. Real deliveries carry such islands as fragments of a network
equivalent: two buses, one branch, two `EquivalentInjection`s without
`referencePriority` or `maxP`.

Runs after `_deEnergizeSourcelessComponents!` and `markIsolatedBuses!`, when
the islands are final. The islands are those the power flow solves
(`Sparlectra._power_flow_islands`): a bus the marking isolates only because it hangs on
the network through closed switches alone belongs to the island of its
switch group, so the parts it joins get one slack, not one each. An island
qualifies when it has no Slack bus, no PV bus and at
least one generator-type unit with infeed (`pVal > 0`, generator sign). An
island whose units only draw power has no source and is already
de-energized; one whose only source regulates its voltage has a PV bus, and
one whose only source is a machine regulating a remote bus is left to
`_attachMachineControls!`, which holds that machine PV at its own bus.

Slack candidates are only the units with infeed; a unit with P = 0 or
drawing power is never one, and a static var compensator neither. Among
them (`_islandSlackUnit`): the strongest `referencePriority`, else the
largest active power reserve `maxP - pVal` on a bus of degree greater than 1
(degree: neighbour buses over closed in-service branches and closed links;
the only bus of a single-bus island counts as such a bus), else the largest
infeed on a bus of degree greater than 1, and only when no such bus carries
infeed, the largest infeed on a bus of degree 1. The unit becomes the slack
exactly as the importer's own slack choice does it: `referencePri` on its
bus, `isRegulated`, and as voltage setpoint the start voltage of its bus
(the SV voltage, 1.0 pu without one); none of its units regulates, so no
local setpoint exists. Each choice is a `warning:` message naming island,
unit, reason and setpoint, and one notice counts the islands. Called only
with `multi_slack = true`.

The unit type of a unit's `DynamicData` (CGMES DY import, see
`cgmes_import.read_dynamics`) narrows the candidates: a grid-following unit
cannot hold an island and is never the slack. An island whose units are all
grid-following is switched off like a part without source (injections
zeroed, branches out of service, buses isolated) with one `warning:` each
and one counting notice. When grid-following units carry all the infeed but
another unit of unknown or other type sits in the island, no slack is chosen
and a `warning:` says so. Synchronous and grid-forming units, and units
without `DynamicData`, are candidates as described above.
"""
function _slackForIslandsWithoutReference!(net, store::CGMESStore, ctx::_MapCtx)
  nbus = length(net.nodeVec)
  # neighbour buses over closed in-service branches and closed links; a
  # parallel circuit to the same neighbour does not raise the degree
  nbrs = [Set{Int}() for _ = 1:nbus]
  for br in net.branchVec
    Sparlectra._branch_terminal_state(br) === :closed || continue
    f, t = Int(br.fromBus), Int(br.toBus)
    (f != t && 1 <= f <= nbus && 1 <= t <= nbus) || continue
    push!(nbrs[f], t)
    push!(nbrs[t], f)
  end
  for l in net.linkVec
    l.status == 1 || continue
    f, t = Int(l.fromBus), Int(l.toBus)
    (f != t && 1 <= f <= nbus && 1 <= t <= nbus) || continue
    push!(nbrs[f], t)
    push!(nbrs[t], f)
  end
  degree = length.(nbrs)
  units = Dict{Int,Vector{Int}}()   # bus -> prosumer indices of candidate units
  for (k, ps) in enumerate(net.prosumpsVec)
    Sparlectra.isGenerator(ps) || continue
    ps.comp.cTyp == Sparlectra.StaticVarCompensator && continue
    push!(get!(Vector{Int}, units, Int(ps.comp.cFrom_bus)), k)
  end
  # unit names for the messages, from the CGMES object each unit came from
  name_of = Dict{Int,String}()
  for (mrid, (kind, k)) in ctx.element_by_mrid
    kind in (:generator, :external_grid) || continue
    o = get(store.objects, mrid, nothing)
    name_of[k] = o === nothing ? mrid : something(str(o, :name), mrid)
  end
  n_islands = 0
  # the islands of the power flow: a part that hangs on the rest through a
  # switch-only bus is no island of its own, it shares the reference of the
  # rest, and a slack of its own would be a second reference in one island
  islands, held, _ = Sparlectra._power_flow_islands(net)
  # a switch-only bus is at 0 pu after the marking; its group voltage sits
  # on the first member that is not isolated
  fed = Sparlectra._linkFedIsolatedBuses(net)
  n_off = 0
  for (buses, has_ref) in zip(islands, held)
    has_ref && continue
    ks_all = reduce(vcat, (get(units, b, Int[]) for b in buses); init = Int[])
    any(k -> something(net.prosumpsVec[k].pVal, 0.0) > 0.0, ks_all) || continue
    around = something(_busDictName(net, first(buses)), net.nodeVec[first(buses)].comp.cName)
    # a unit stated grid-following (DynamicData.unit_type) cannot hold an
    # island and is no candidate; a unit without DynamicData stays one
    ks = filter(k -> !_isGridFollowingUnit(net.prosumpsVec[k]), ks_all)
    if !any(k -> something(net.prosumpsVec[k].pVal, 0.0) > 0.0, ks)
      if isempty(ks)
        # every unit of the island is grid-following: nothing forms voltage
        # and frequency, the island cannot stand in operation
        nload, p_load, p_fed, p_gen = _deEnergizeBuses!(net, buses)
        n_off += 1
        push!(ctx.messages, "warning: electrical island of $(length(buses)) bus(es) around $(around) has infeed from grid-following units only and cannot stand in operation; switched off, $(nload) injection(s) zeroed, generation P = $(round(p_gen; digits = 2)) MW dropped, load P = $(round(p_load; digits = 2)) MW unsupplied$(_droppedLoadInfeed(p_fed))")
      else
        push!(ctx.messages, "warning: electrical island of $(length(buses)) bus(es) around $(around) has infeed but no reference, and every unit with infeed is grid-following; no slack chosen")
      end
      continue
    end
    k, reason = _islandSlackUnit(net, ks, degree)
    ps = net.prosumpsVec[k]
    bus = Int(ps.comp.cFrom_bus)
    vm = something(net.nodeVec[get(fed, bus, bus)]._vm_pu, 1.0)
    ps.referencePri = bus
    ps.isRegulated = true
    ps.vm_pu = vm
    # the reference angle as for the slack of the flush in _mapInjections!
    ps.va_deg = something(net.nodeVec[get(fed, bus, bus)]._va_deg, 0.0)
    n_islands += 1
    at = something(_busDictName(net, bus), net.nodeVec[bus].comp.cName)
    # a converter can hold an island only when it forms the grid; without
    # DynamicData the delivery does not say which kind it is, so the warning
    # states the assumption
    converter = (ps.comp.cTyp == Sparlectra.PowerElectronicsConnection && ps.dynamics === nothing) ? "; its unit type is unknown and assumed grid-following, an island of grid-following units only cannot stand in operation" : ""
    push!(ctx.messages, "warning: electrical island of $(length(buses)) bus(es) around $(around) has infeed but no reference; $(get(name_of, k, "unit $(k)")) at $(at) made its slack ($(reason)), voltage setpoint $(round(vm; digits = 4)) pu from the start voltage of its bus$(converter)")
  end
  # the switched-off islands lost their branches: mark them isolated (voltage
  # 0) with the rule of the import
  n_off > 0 && Sparlectra.markIsolatedBuses!(net = net)
  (n_islands == 0 && n_off == 0) && return nothing
  Sparlectra.refreshBusTypesFromProsumers!(net)
  n_islands > 0 && push!(ctx.messages, "notice: $(n_islands) electrical island(s) without reference got a slack from the importer (one warning each)")
  n_off > 0 && push!(ctx.messages, "notice: $(n_off) electrical island(s) with grid-following units only switched off (one warning each)")
  return nothing
end

# A unit whose DynamicData states it grid-following. A unit without
# DynamicData is not: its type is unknown and the island rules of the import
# treat it as before.
_isGridFollowingUnit(ps)::Bool = ps.dynamics !== nothing && ps.dynamics.unit_type == Sparlectra.GridFollowingUnit

"""
Angle-reference nodes declared by the SV profile: `TopologicalIsland` names
one `AngleRefTopologicalNode` per island. That is the exporting tool's own
reference choice, the most reliable slack information a delivery can carry,
and the only one that scales to multi-island models (ReliCapGrid/Svedala
ships two islands).

!!! warning "Two different island notions"
    A CIM `TopologicalIsland` is the exporting tool's partitioning and need
    **not** agree with the AC island Sparlectra sees after evaluating the
    switch states, which is why `_selectSlackBuses` re-checks the references
    against `electricalIslandOfBus` instead of trusting the CIM island count.
    Instructive history: Svedala *appeared* to disagree (both references in
    one electrical island), but that was an importer artifact: switches with
    SSH `Equipment.inService = false` were treated as closed and merged the
    two islands. With `inService` honored, the electrical view agrees with the
    delivery's two `TopologicalIsland`s exactly. The cross-check stays: it
    costs little and catches deliveries where the disagreement is real.

Returns `TopologicalNode` mRID => island name (`IdentifiedObject.name` of the
`TopologicalIsland`, its mRID when it has none), for the slack messages.
"""
function _svAngleRefNodes(store::CGMESStore)::Dict{String,String}
  out = Dict{String,String}()
  for isl in objectsOf(store, :TopologicalIsland)
    tn = get(isl.refs, :AngleRefTopologicalNode, nothing)
    tn === nothing || (out[tn] = something(str(isl, :name), isl.mrid))
  end
  return out
end

# Do the EquivalentInjections collected on one boundary node describe the
# same exchange seen from both sides? True only when both sides are present
# AND their p values (approximately) cancel, the defining property of an
# assembled AC X-node pair. A non-cancelling pair is two different devices
# (HVDC converters; the residue is the DC loss) and must be kept. The
# tolerance only absorbs rounding in the delivery (largest observed residue of
# a genuine pair: 1.6 % / 0.025 MW; smallest genuine non-pair: 23.2 % / 22 MW).
function _eiPairCancels(ps::Vector{Float64})::Bool
  length(ps) >= 2 || return false
  pmax = maximum(abs, ps)
  return abs(sum(ps)) <= max(0.1 * pmax, 1.0)
end

# Boundary-node identity of an EquivalentInjection's terminal (see the
# ei_p_bnode comment in _mapInjections!): the terminal's ConnectivityNode if
# that CN comes from a boundary file (CGMES 3.0 style), else its
# TopologicalNode if that TN does (2.4.15 TP_BD style), else `nothing`, the
# EI does not sit on a boundary node at all (e.g. an internal equivalent).
function _eiBoundaryNode(store::CGMESStore, topo::CGMESTopology, ei::CIMObject)::Union{Nothing,String}
  ts = get(topo.terminals, ei.mrid, nothing)
  (ts === nothing || isempty(ts)) && return nothing
  t = ts[1]
  cn = get(t.refs, :ConnectivityNode, nothing)
  (cn !== nothing && cn in store.boundary) && return cn
  tn = tnOfTerminal(t)
  (tn !== nothing && tn in store.boundary) && return tn
  return nothing
end

# Does an EquivalentInjection sit on a boundary node? Its node comes from a
# boundary file (_eiBoundaryNode), or the delivery marks the terminal's
# ConnectivityNode or TopologicalNode with `boundaryPoint = true` (a delivery
# that ships the boundary nodes without the boundary profiles).
function _eiOnBoundary(store::CGMESStore, topo::CGMESTopology, ei::CIMObject)::Bool
  _eiBoundaryNode(store, topo, ei) === nothing || return true
  ts = get(topo.terminals, ei.mrid, nothing)
  (ts === nothing || isempty(ts)) && return false
  t = ts[1]
  for node in (get(t.refs, :ConnectivityNode, nothing), tnOfTerminal(t))
    node === nothing && continue
    o = get(store.objects, node, nothing)
    (o !== nothing && something(boolval(o, :boundaryPoint), false)) && return true
  end
  return false
end

# names for a notice: up to five, then an ellipsis
_shortList(names::Vector{String})::String = length(names) <= 5 ? join(names, ", ") : join(names[1:5], ", ") * ", ..."

"""
Decide which remote machine-control plans survive (#294 point 3).

A plan falls back to the Stage-1 behavior (machine held PV at its own bus)
when its target bus is not part of the built network, is already voltage-held
(PV/slack), or the machine's own bus is voltage-held. Fallbacks themselves
create new PV buses, so the scan runs to a fixpoint before the surviving
plans claim their targets (one controller per target bus, document order).
"""
function _resolveMachineControlPlans!(pending::Vector{NamedTuple}, net, slack_buses::Set{String}, messages::Vector{String})
  vheld = Set{String}(slack_buses)
  for inj in pending
    inj.vm_pu === nothing || push!(vheld, inj.bus)
  end
  fallback = function (i, inj, reason)
    push!(messages, "notice: $(something(get(inj, :name, nothing), inj.bus)), remote voltage control not attachable ($(reason)); held PV at its own bus")
    pending[i] = merge(inj, (vm_pu = _pvVoltage(net, inj.bus, nothing), isRegulated = true, remote_target = nothing))
    push!(vheld, inj.bus)
  end
  changed = true
  while changed
    changed = false
    for (i, inj) in enumerate(pending)
      rt = get(inj, :remote_target, nothing)
      rt === nothing && continue
      if !haskey(net.busDict, rt.bus)
        fallback(i, inj, "target bus $(rt.bus) is not part of the built network")
        changed = true
      elseif rt.bus in vheld
        fallback(i, inj, "target bus $(rt.bus) is already voltage-held (PV/slack)")
        changed = true
      elseif inj.bus in vheld
        fallback(i, inj, "its own bus is voltage-held (PV/slack)")
        changed = true
      end
    end
  end
  claimed = Set{String}()
  for (i, inj) in enumerate(pending)
    rt = get(inj, :remote_target, nothing)
    rt === nothing && continue
    if rt.bus in claimed
      fallback(i, inj, "target bus $(rt.bus) is already claimed by another machine controller")
    else
      push!(claimed, rt.bus)
    end
  end
  return nothing
end

"""
    _detectHvdcPairs(store) -> Vector{@NamedTuple{from::String, to::String, extra::Int, segments::Bool}}

Group the HVDC converters (`VsConverter`/`CsConverter`) by their DC-side
connectivity: every `ACDCConverterDCTerminal`/`DCTerminal` joins its DC
conducting equipment's nodes (`DCBaseTerminal.DCNode`, node-breaker
deliveries use `DCTopologicalNode`), so two converters joined through DC
nodes, segments, or switches land in one component. A component with
exactly two converters is a link candidate (back-to-back when no
`DCLineSegment` participates, point-to-point otherwise); other component
sizes are reported via `extra` so the caller can notice them. Reads the
raw store only, no electrical DC model.
"""
function _detectHvdcPairs(store::CGMESStore)
  # union-find over DC node mrids
  parent = Dict{String,String}()
  function _find(x::String)
    get!(parent, x, x)
    while parent[x] != x
      parent[x] = parent[parent[x]]
      x = parent[x]
    end
    return x
  end
  _union(a::String, b::String) = (parent[_find(a)] = _find(b); nothing)

  eq_nodes = Dict{String,Vector{String}}()
  for cls in (:ACDCConverterDCTerminal, :DCTerminal)
    for t in objectsOf(store, cls)
      eq = get(t.refs, :DCConductingEquipment, nothing)
      node = get(t.refs, :DCNode, get(t.refs, :DCTopologicalNode, nothing))
      (eq === nothing || node === nothing) && continue
      push!(get!(Vector{String}, eq_nodes, String(eq)), String(node))
    end
  end
  segment_ids = Set{String}(o.mrid for o in objectsOf(store, :DCLineSegment))
  seg_roots = Set{String}()
  for (eq, nodes) in eq_nodes
    for n in nodes
      _union(nodes[1], n)
    end
    eq in segment_ids && push!(seg_roots, _find(nodes[1]))
  end

  by_root = Dict{String,Vector{String}}()
  for cls in (:VsConverter, :CsConverter)
    for conv in objectsOf(store, cls)
      nodes = get(eq_nodes, conv.mrid, nothing)
      nodes === nothing && continue
      push!(get!(Vector{String}, by_root, _find(nodes[1])), conv.mrid)
    end
  end
  pairs = NamedTuple[]
  for (root, convs) in by_root
    if length(convs) == 2
      a, b = sort(convs)
      push!(pairs, (from = a, to = b, extra = 0, segments = _find(root) in Set(_find(r) for r in seg_roots)))
    else
      push!(pairs, (from = first(sort(convs)), to = "", extra = length(convs), segments = false))
    end
  end
  return pairs
end

# Attaching the detected HVDC converter pairs as `addHvdcPairControl!`
# controllers (`cgmes_import.hvdc_mode = :paired_control`) happens in the
# mapping step below: the from side is the converter whose PCC injection
# exports into the link (negative injection); the transfer and the loss come
# from the two SSH operating points, each terminal keeps its SSH reactive
# value as a fixed setpoint. A pair that cannot be attached (converter
# skipped at import, isolated bus, inconsistent snapshot signs, negative
# loss) degrades to the Stage-0 fixed injections with a notice; the import
# itself never fails there.
#
# persistent link record for the result layer (r0.9.9): one HvdcLink per
# detected two-converter pair whose terminals were imported; runs in every
# hvdc_mode. Orientation matches the attach rule (exporting side = from).
# Returns the record or nothing when a converter is missing from the net.
function _push_cgmes_hvdc_link!(net, ctx::_MapCtx, pair)
  pair.extra == 0 || return nothing
  a_idx = get(ctx.prosumer_by_mrid, pair.from, nothing)
  b_idx = get(ctx.prosumer_by_mrid, pair.to, nothing)
  (a_idx === nothing || b_idx === nothing) && return nothing
  pa = net.prosumpsVec[a_idx].pVal === nothing ? 0.0 : net.prosumpsVec[a_idx].pVal
  pb = net.prosumpsVec[b_idx].pVal === nothing ? 0.0 : net.prosumpsVec[b_idx].pVal
  f_idx, t_idx = pa <= pb ? (a_idx, b_idx) : (b_idx, a_idx)
  fbus = Sparlectra.getPosumerBusIndex(net.prosumpsVec[f_idx])
  tbus = Sparlectra.getPosumerBusIndex(net.prosumpsVec[t_idx])
  names = Sparlectra._bus_name_by_idx(net)
  link = Sparlectra.HvdcLink(string("HVDC_", get(names, fbus, string(fbus)), "_", get(names, tbus, string(tbus))), fbus, tbus, f_idx, t_idx, 1, :cgmes, pair.segments ? :p2p : :b2b, nothing)
  push!(net.hvdcLinks, link)
  return link
end

function _attachHvdcPairs!(net, ctx::_MapCtx, store::CGMESStore)
  for pair in _detectHvdcPairs(store)
    _push_cgmes_hvdc_link!(net, ctx, pair)
    if pair.extra != 0
      push!(ctx.messages, "notice: HVDC component around $(pair.from) has $(pair.extra) converter(s), only two-converter links are paired, kept as fixed injections")
      continue
    end
    a_idx = get(ctx.prosumer_by_mrid, pair.from, nothing)
    b_idx = get(ctx.prosumer_by_mrid, pair.to, nothing)
    if a_idx === nothing || b_idx === nothing
      push!(ctx.messages, "notice: HVDC pair $(pair.from)/$(pair.to), a converter was not imported (skipped or out of service), kept as fixed injections")
      continue
    end
    pa = net.prosumpsVec[a_idx].pVal === nothing ? 0.0 : net.prosumpsVec[a_idx].pVal
    pb = net.prosumpsVec[b_idx].pVal === nothing ? 0.0 : net.prosumpsVec[b_idx].pVal
    # the exporting side carries the negative injection
    f_idx, t_idx = pa <= pb ? (a_idx, b_idx) : (b_idx, a_idx)
    transfer = -(net.prosumpsVec[f_idx].pVal === nothing ? 0.0 : net.prosumpsVec[f_idx].pVal)
    to_p = net.prosumpsVec[t_idx].pVal === nothing ? 0.0 : net.prosumpsVec[t_idx].pVal
    loss = transfer - to_p
    if transfer < 0.0 || loss < 0.0
      push!(ctx.messages, "notice: HVDC pair $(pair.from)/$(pair.to), SSH operating points are inconsistent (transfer $(round(transfer; digits = 1)) MW, loss $(round(loss; digits = 1)) MW), kept as fixed injections")
      continue
    end
    bus_names = Sparlectra._bus_name_by_idx(net)
    from_bus = Sparlectra._effective_bus_name(bus_names, net, Sparlectra.getPosumerBusIndex(net.prosumpsVec[f_idx]))
    to_bus = Sparlectra._effective_bus_name(bus_names, net, Sparlectra.getPosumerBusIndex(net.prosumpsVec[t_idx]))
    try
      Sparlectra.addHvdcPairControl!(
        net;
        from_bus = from_bus,
        to_bus = to_bus,
        p_transfer_mw = transfer,
        loss_mw = loss,
        from_q_mvar = net.prosumpsVec[f_idx].qVal,
        to_q_mvar = net.prosumpsVec[t_idx].qVal,
        name = string("HVDC_", from_bus, "_", to_bus),
        from_prosumer = f_idx,
        to_prosumer = t_idx,
      )
      kind = pair.segments ? "point-to-point" : "back-to-back"
      push!(ctx.messages, "notice: HVDC pair attached as controller ($(kind)): $(from_bus) -> $(to_bus), transfer $(round(transfer; digits = 1)) MW, loss $(round(loss; digits = 1)) MW")
    catch err
      push!(ctx.messages, "notice: HVDC pair $(from_bus)/$(to_bus) could not be attached ($(sprint(showerror, err))), kept as fixed injections")
    end
  end
  return nothing
end

"""
Attach the surviving machine-control plans as `MachineVoltageControl`
instances. Runs only after `refreshBusTypesFromProsumers!` and isolation
marking, so the final node types and `isoNodes` are authoritative, a plan
whose buses did not end up solvable PQ is dropped with a notice, mirroring
`_disableIneffectiveTapControllers!`.

A machine whose remote regulation is the only source of its island (no
slack, no PV bus, no unit with infeed; `_deEnergizeSourcelessComponents!`
kept the island energized for that regulation) gets no controller either:
the importer gives such an island no slack, and the solve would hold the
machine bus at its start voltage as island reference, where a Q controller
has nothing to act on. The machine is held PV at its own bus at its start
voltage instead, as without `machine_control`, with a notice.
"""
function _attachMachineControls!(net, ctx::_MapCtx)
  # per island of the power flow (_power_flow_islands, switch-only buses
  # with their group): does a node hold a voltage (Slack or PV), and does a
  # generator-type unit inject (pVal > 0)
  islands, held, _ = Sparlectra._power_flow_islands(net)
  island_of = Dict{Int,Int}()
  for (i, buses) in enumerate(islands), b in buses
    island_of[b] = i
  end
  infeed = falses(length(islands))
  for ps in net.prosumpsVec
    i = get(island_of, Int(ps.comp.cFrom_bus), 0)
    i != 0 && Sparlectra.isGenerator(ps) && something(ps.pVal, 0.0) > 0.0 && (infeed[i] = true)
  end
  for plan in ctx.machine_ctrl_plans
    tidx = get(net.busDict, plan.target_bus, nothing)
    oidx = get(net.busDict, plan.bus, nothing)
    reason = if tidx === nothing
      "target bus $(plan.target_bus) not found"
    elseif tidx in net.isoNodes
      "target bus $(plan.target_bus) is isolated"
    elseif net.nodeVec[tidx]._nodeType != Sparlectra.PQ
      "target bus $(plan.target_bus) is $(net.nodeVec[tidx]._nodeType)"
    elseif oidx !== nothing && oidx in net.isoNodes
      "machine bus $(plan.bus) is isolated"
    elseif net.prosumpsVec[plan.prosumer_idx].referencePri !== nothing
      # _slackForIslandsWithoutReference! made the machine the slack of its
      # island; a slack holds its own voltage and takes no Q controller
      "the machine is the slack of its island"
    else
      nothing
    end
    if reason !== nothing
      push!(ctx.messages, "notice: machine control $(plan.name), skipped ($(reason))")
      continue
    end
    isl = get(island_of, oidx, 0)
    if isl != 0 && !held[isl] && !infeed[isl]
      ps = net.prosumpsVec[plan.prosumer_idx]
      ps.isRegulated = true
      ps.vm_pu = _pvVoltage(net, plan.bus, nothing)
      Sparlectra.refreshBusTypesFromProsumers!(net)
      held[isl] = true
      push!(ctx.messages, "notice: machine control $(plan.name) not attached (its island has no slack, no PV bus and no infeed); held PV at its own bus")
      continue
    end
    try
      Sparlectra.addMachineVoltageControl!(
        net;
        bus = plan.bus,
        target_bus = plan.target_bus,
        target_vm_pu = plan.target_vm_pu,
        qmin_mvar = plan.qmin,
        qmax_mvar = plan.qmax,
        prosumer_index = plan.prosumer_idx,
        name = plan.name,
      )
      push!(ctx.messages, "machine control: $(plan.name) → voltage $(round(plan.target_vm_pu; digits = 4)) pu at $(plan.target_bus)")
    catch err
      # the importer degrades gracefully: an unattachable controller becomes a
      # visible notice (the machine stays a plain PQ injection with its SSH q)
      push!(ctx.messages, "notice: machine control $(plan.name), not attached ($(sprint(showerror, err)))")
    end
  end
  return nothing
end

"""
    _hvdcCurrentSourceTransfer(a, b) -> (p_dc, from_mrid, to_mrid, note)

DC power of a current-source (LCC) HVDC link whose two converters state
their operating point as a DC current and a DC voltage instead of an active
power, so `ACDCConverter.p` reads zero at BOTH ends.

The operating point is fully determined all the same: one converter
declares `CsConverter.targetIdc` (its `pPccControl` is `dcCurrent`), the
other declares `ACDCConverter.targetUdc` (`pPccControl` is `dcVoltage`), and
their product is the DC power. `from_mrid` is the end that DRAWS from the AC
grid, taken from `CsConverter.operatingMode` (`rectifier` draws,
`inverter` delivers); without a usable operating mode the current-controlled
end is treated as the rectifier, which is the normal LCC arrangement.

Returns `p_dc === nothing` when the pair states neither a current nor a
voltage target, so the caller can report the gap instead of inventing a
number. The DC line resistance is deliberately NOT applied here: it would
split the transfer into a different value per end (on FullGrid 75.625 MW
drawn against 75.0 MW delivered, matching SV to the converters' pole
losses), which needs the DC-side model this importer does not build. The
`note` states the inputs so the import message can show them.
"""
function _hvdcCurrentSourceTransfer(a::CIMObject, b::CIMObject)
  idc_a = num(a, :targetIdc, 0.0)
  idc_b = num(b, :targetIdc, 0.0)
  udc_a = num(a, :targetUdc, 0.0)
  udc_b = num(b, :targetUdc, 0.0)
  # the current target lives on one converter, the voltage target on the
  # other; either assignment is acceptable, both missing is not
  cur, volt = idc_a > 0.0 ? (a, b) : (b, a)
  idc = max(idc_a, idc_b)
  udc = max(udc_a, udc_b)
  (idc > 0.0 && udc > 0.0) || return (nothing, "", "", "targetIdc=$(round(idc; digits = 1)) A, targetUdc=$(round(udc; digits = 1)) kV")
  # kV times A is kW
  p_dc = udc * idc / 1000.0
  mode_cur = String(something(str(cur, :operatingMode), ""))
  mode_volt = String(something(str(volt, :operatingMode), ""))
  drawing = occursin("inverter", lowercase(mode_cur)) && !occursin("inverter", lowercase(mode_volt)) ? volt : cur
  delivering = drawing === cur ? volt : cur
  note = "targetUdc $(round(udc; digits = 1)) kV times targetIdc $(round(idc; digits = 1)) A"
  return (p_dc, drawing.mrid, delivering.mrid, note)
end

"""
    _hvdcDerivedSlackPower(store, ctx) -> Dict{String,Float64}

Derived SSH active power for HVDC converters that regulate the DC VOLTAGE
instead of their PCC power, keyed by converter mRID (load convention, the
same sign `ACDCConverter.p` uses).

Why this exists: in CIM the active power of a DC-voltage-controlling
converter is a RESULT of the DC power balance, not a setpoint. An SSH
snapshot may therefore legitimately carry `ACDCConverter.p = 0` for that
end while the opposite end carries the scheduled transfer. Reading both
ends literally, as the Stage-0 mapping did, injects the transfer at one end
and nothing at the other: the link swallows its whole throughput. On the
FullGrid test configuration that is 150 MW entering at BE_Busbar_HVDC2 and
0 MW leaving at BE_Busbar_HVDC1, and it produced the largest angle
deviations of the whole delivery against the SV profile (measured
max |dva| 21.9 degrees, halved to 10.8 by this derivation).

The derivation is deliberately narrow and fires only when the current
reading is provably wrong: a detected two-converter link, exactly one end
at (numerically) zero power, that end declaring `targetUdc > 0`, and the
other end carrying a nonzero schedule. The derived value is the opposite
end's power reduced by both converters' declared `idleLoss`, the only loss
term CIM states as a plain active power (`resistiveLoss` is an ohmic value
and `switchingLoss` a power per current, both of which need the DC-side
model this importer deliberately does not build).

A pair with BOTH ends at zero is NOT derived: its transfer lives only in
`targetUdc`/`CsConverter.pPccControl` and cannot be recovered from the SSH
power at all. Those pairs get a notice instead of a guessed number.
"""
function _hvdcDerivedSlackPower(store::CGMESStore, ctx::_MapCtx)
  derived = Dict{String,Float64}()
  # 0.01 MW: below any real schedule, above the noise of a "0e+000" literal
  const_zero_tol = 0.01
  by_mrid = Dict{String,CIMObject}()
  for cls in (:VsConverter, :CsConverter)
    for c in objectsOf(store, cls)
      by_mrid[c.mrid] = c
    end
  end
  for pair in _detectHvdcPairs(store)
    pair.extra == 0 || continue
    a = get(by_mrid, pair.from, nothing)
    b = get(by_mrid, pair.to, nothing)
    (a === nothing || b === nothing) && continue
    pa = num(a, :p, 0.0)
    pb = num(b, :p, 0.0)
    za = abs(pa) < const_zero_tol
    zb = abs(pb) < const_zero_tol
    if za && zb
      # A current-source (LCC) link states its operating point as a DC
      # CURRENT and a DC VOLTAGE rather than as an active power, so both ends
      # legitimately read p = 0. The setpoint is present, just in two
      # attributes on two converters: targetIdc on the current-controlling
      # end and targetUdc on the voltage-controlling one, whose product is
      # the DC power. Measured on FullGrid: 150 kV times 500 A is 75 MW, and
      # the SV profile confirms it (udc 150.0/151.25 kV, idc 500 A on both
      # converters). Same rule as the p-based derivation above, only the
      # source attribute differs.
      p_dc, from_mrid, to_mrid, note = _hvdcCurrentSourceTransfer(a, b)
      if p_dc === nothing
        push!(ctx.messages, "notice: HVDC pair $(something(str(a, :name), a.mrid)) / $(something(str(b, :name), b.mrid)) carries no SSH active power at either end and no usable targetUdc/targetIdc pair either ($(note)); its transfer is not recoverable from the delivery and both ends stay at zero injection")
        continue
      end
      loss = num(a, :idleLoss, 0.0) + num(b, :idleLoss, 0.0)
      # load convention: the rectifier draws the DC power from the AC grid,
      # the inverter delivers it less the declared idle losses
      derived[from_mrid] = p_dc
      derived[to_mrid] = -p_dc + loss
      nm(m) = something(str(get(by_mrid, m, a), :name), m)
      push!(ctx.messages, "notice: HVDC pair $(nm(from_mrid)) / $(nm(to_mrid)) states its operating point as DC current and DC voltage, not as active power; transfer derived as $(round(p_dc; digits = 3)) MW ($(note)) instead of leaving both ends at zero")
      continue
    end
    (za || zb) || continue                       # both scheduled: nothing to derive
    slack, other = za ? (a, b) : (b, a)
    p_other = za ? pb : pa
    # the zero end must actually declare DC-voltage control, otherwise a
    # genuine zero schedule would be overwritten with an invented transfer
    num(slack, :targetUdc, 0.0) > 0.0 || continue
    loss = num(a, :idleLoss, 0.0) + num(b, :idleLoss, 0.0)
    # load convention: the other end draws p_other from the AC grid, so the
    # DC-slack end delivers that much less the losses, i.e. a negative p
    p_slack = -(p_other) + (p_other > 0 ? loss : -loss)
    derived[slack.mrid] = p_slack
    push!(ctx.messages, "notice: $(something(str(slack, :name), slack.mrid)) regulates the DC voltage (targetUdc=$(round(num(slack, :targetUdc, 0.0); digits = 1)) kV), so its SSH active power is a result, not a setpoint; PCC power derived from $(something(str(other, :name), other.mrid)) as $(round(p_slack; digits = 2)) MW (transfer $(round(p_other; digits = 1)) MW less $(round(loss; digits = 2)) MW idle loss) instead of the literal $(round(za ? pa : pb; digits = 1)) MW")
  end
  return derived
end

# EnergySource (cgmes_import.map_energy_source). Real deliveries carry
# generation as EnergySource; without it
# that generation is missing and a comparison with the delivered SV state
# says nothing. CGMES 2.4.15 gives EnergySource.activePower and
# reactivePower in the SSH in load sign convention ("positive sign means
# flow out from a node"), the same convention as SvPowerFlow at its
# terminal, so a generating source has activePower < 0 and the injection is
# the negated value. Mapped as a fixed PQ injection of type GENERATOR (an
# Injection, so it counts as generation in the bus sums): no voltage
# regulation, no Q limits, and slack_candidate = false keeps it out of the
# importer's slack choice from the delivery's reference data. In an island
# without other reference it can still become the slack, from the importer
# (_slackForIslandsWithoutReference!) or, with multi_slack = false, from the
# island-wise solve (see _noticeEnergySourceOnlyIslands!). Switched off, the
# objects go to ctx.skipped (not mapped), so compareWithSV and the
# measurement key index leave their terminals out, and one notice counts
# them.
function _mapEnergySources!(pending::Vector{NamedTuple}, net, store, topo, created, svmap, ctx::_MapCtx)
  sources = collect(CIMObject, objectsOf(store, :EnergySource))
  if !ctx.map_energy_source
    isempty(sources) && return nothing
    foreach(es -> push!(ctx.skipped, es.mrid), sources)
    # count and P sum stand in the table of classes not mapped (one sign
    # convention); this line names the switch
    push!(ctx.messages, "notice: EnergySource not mapped, switch cgmes_import.map_energy_source (count and SSH P sum in the table of classes not mapped)")
    return nothing
  end
  n_mapped = 0
  n_skipped = 0
  p_sum = 0.0
  q_sum = 0.0
  for es in sources
    t = busOfEquipment(topo, es, 1)
    name = something(str(es, :name), es.mrid)
    reason = if t === nothing
      "no TopologicalNode"
    elseif !_conn(ctx, t.connected)
      "terminal disconnected"
    elseif !_inService(ctx, es)
      "out of service, SSH inService=false"
    else
      nothing
    end
    if reason !== nothing
      push!(ctx.messages, "skip: EnergySource $(name) ($(reason))")
      push!(ctx.skipped, es.mrid)
      n_skipped += 1
      continue
    end
    bus = _ensureBus!(net, created, topo, svmap, t.tn)
    p = -num(es, :activePower, 0.0)               # load convention -> injection
    q = -num(es, :reactivePower, 0.0)
    push!(
      pending,
      (
        bus = bus,
        type = "GENERATOR",
        p = p,
        q = q,
        pMin = nothing,
        pMax = nothing,
        qMin = nothing,
        qMax = nothing,
        vm_pu = nothing,
        isRegulated = false,
        slack_candidate = false,
        mrid = es.mrid,
      ),
    )
    n_mapped += 1
    p_sum += p
    q_sum += q
  end
  push!(ctx.messages, "notice: map_energy_source = true: $(n_mapped) EnergySource object(s) mapped as PQ injections (P = $(round(p_sum; digits = 2)) MW, Q = $(round(q_sum; digits = 2)) MVAr injected), $(n_skipped) skipped")
  return nothing
end

# The CGMES 3.0 subclasses of PowerElectronicsUnit. The reader keeps every
# object under its own class, so the abstract base is listed only in case a
# delivery instantiates it.
const _PE_UNIT_CLASSES = [:PowerElectronicsUnit, :BatteryUnit, :PhotoVoltaicUnit, :PowerElectronicsWindUnit]

# Active power limits of each PowerElectronicsConnection from its units
# (PowerElectronicsUnit.minP/maxP in EQ, nameplate values in generator sign
# like GeneratingUnit.min/maxOperatingP): the sum over the units that
# reference the connection, `nothing` for a bound that one of them does not
# state, so a partial sum never reads as the full limit. A placeholder bound
# of a unit (_pLimit) counts as not stated; `placeholder` marks the
# connections with one for the notice.
function _powerElectronicsPLimits(store::CGMESStore)::Dict{String,@NamedTuple{lo::Union{Nothing,Float64}, hi::Union{Nothing,Float64}, placeholder::Bool}}
  units = Dict{String,Vector{CIMObject}}()
  for cls in _PE_UNIT_CLASSES, u in objectsOf(store, cls)
    pec = get(u.refs, :PowerElectronicsConnection, nothing)
    pec === nothing || push!(get!(Vector{CIMObject}, units, pec), u)
  end
  out = Dict{String,@NamedTuple{lo::Union{Nothing,Float64}, hi::Union{Nothing,Float64}, placeholder::Bool}}()
  for (pec, us) in units
    raw_lo = [num(u, :minP) for u in us]
    raw_hi = [num(u, :maxP) for u in us]
    lo = _pLimit.(raw_lo)
    hi = _pLimit.(raw_hi)
    placeholder = count(isnothing, lo) > count(isnothing, raw_lo) || count(isnothing, hi) > count(isnothing, raw_hi)
    out[pec] = (lo = any(isnothing, lo) ? nothing : sum(something, lo), hi = any(isnothing, hi) ? nothing : sum(something, hi), placeholder = placeholder)
  end
  return out
end

# PowerElectronicsConnection (CGMES 3.0: wind and solar parks, batteries
# behind a converter). Without it their
# generation is missing from the network. Behind the switch
# cgmes_import.map_power_electronics (default off like map_energy_source;
# the release decides both): off, every connection goes to ctx.skipped and
# one notice names the switch, while count and P sum stand in the table of
# classes not mapped. The SSH p and q of the connection use the load sign
# like EnergySource.activePower ("positive sign means flow out from a
# node"), so the injection is the negated value. The connection is a RegulatingCondEq:
# with RegulatingCondEq.controlEnabled (SSH, absent counts as on, as for the
# SVC) and a voltage RegulatingControl that is enabled and has a target
# above zero it regulates like a SynchronousMachine (_voltageSetpoint: PV at
# its own bus on the target; a target outside vset_min_pu/vset_max_pu holds
# it PV at its bus voltage with a warning; a remote target holds it PV at its
# own bus or, with machine_control, gives it an outer-loop controller). No
# target, a target <= 0, controlEnabled = false, a disabled control or a
# control in another mode leave a fixed PQ injection. The Q
# limits are the minQ/maxQ hull of the machines; a regulating connection
# without usable minQ/maxQ gets the machines' wide fallback from ratedS, a
# PQ one keeps none. P limits come from its PowerElectronicsUnits. ratedU,
# ratedS and maxIFault go to the short-circuit harvest
# (collectShortCircuitData), where the connection is a current source. An
# injection like EnergySource, and like it no slack candidate of the delivery
# (slack_candidate = false; CGMES gives a connection no referencePriority); in
# an island without other reference it can still become the slack
# (_slackForIslandsWithoutReference!). The prosumer carries the component type
# PowerElectronicsConnection, so the short-circuit engine neither takes its
# maxP as the size of a machine at the same bus nor counts it as a
# synchronous source.
function _mapPowerElectronicsConnections!(pending::Vector{NamedTuple}, net, store, topo, created, svmap, ctx::_MapCtx)
  pecs = collect(CIMObject, objectsOf(store, :PowerElectronicsConnection))
  isempty(pecs) && return nothing
  if !ctx.map_power_electronics
    foreach(pec -> push!(ctx.skipped, pec.mrid), pecs)
    # same form as the EnergySource notice: the count and P sum stand in the
    # table of classes not mapped, this line names the switch
    push!(ctx.messages, "notice: PowerElectronicsConnection not mapped, switch cgmes_import.map_power_electronics (count and SSH P sum in the table of classes not mapped)")
    return nothing
  end
  plimits = _powerElectronicsPLimits(store)
  n_mapped = 0
  n_regulating = 0
  n_skipped = 0
  p_sum = 0.0
  q_sum = 0.0
  for pec in pecs
    t = busOfEquipment(topo, pec, 1)
    name = something(str(pec, :name), pec.mrid)
    reason = if t === nothing
      "no TopologicalNode"
    elseif !_conn(ctx, t.connected)
      "terminal disconnected"
    elseif !_inService(ctx, pec)
      "out of service, SSH inService=false"
    else
      nothing
    end
    if reason !== nothing
      push!(ctx.messages, "skip: PowerElectronicsConnection $(name) ($(reason))")
      push!(ctx.skipped, pec.mrid)
      n_skipped += 1
      continue
    end
    bus = _ensureBus!(net, created, topo, svmap, t.tn)
    p = -num(pec, :p, 0.0)                        # load convention -> injection
    q = -num(pec, :q, 0.0)
    vset, vctrl, remote_bus = if something(boolval(pec, :controlEnabled), true)
      _voltageSetpoint(store, topo, pec, t.tn, ctx.messages; vset_min_pu = ctx.vset_min_pu, vset_max_pu = ctx.vset_max_pu, allow_remote = ctx.machine_control)
    else
      (nothing, false, nothing)
    end
    is_remote = remote_bus !== nothing
    qlo = _limitHullMin(num(pec, :minQ), num(pec, :maxQ))
    qhi = _limitHullMax(num(pec, :minQ), num(pec, :maxQ))
    if qlo === nothing || qhi === nothing || (qhi - qlo) < 1e-6
      qlo = nothing
      qhi = nothing
      if vctrl
        qwide = max(num(pec, :ratedS, 0.0), net.baseMVA)
        qlo = -qwide
        qhi = qwide
        push!(ctx.messages, "notice: PowerElectronicsConnection $(name) regulates voltage and has no usable minQ/maxQ, using ±$(qwide) MVAr")
      end
    end
    lim = get(plimits, pec.mrid, (lo = nothing, hi = nothing, placeholder = false))
    pmin, pmax = lim.lo, lim.hi
    lim.placeholder && (ctx.p_limit_placeholders["PowerElectronicsConnection"] = get(ctx.p_limit_placeholders, "PowerElectronicsConnection", 0) + 1)
    push!(
      pending,
      (
        bus = bus,
        type = "POWERELECTRONICSCONNECTION",
        p = p,
        q = q,
        pMin = pmin,
        pMax = pmax,
        qMin = qlo,
        qMax = qhi,
        vm_pu = (vctrl && !is_remote) ? _pvVoltage(net, bus, vset) : nothing,
        isRegulated = vctrl && !is_remote,
        slack_candidate = false,
        name = name,
        remote_target = is_remote ? (bus = remote_bus, vset_pu = vset) : nothing,
        mrid = pec.mrid,
      ),
    )
    n_mapped += 1
    vctrl && (n_regulating += 1)
    p_sum += p
    q_sum += q
  end
  push!(ctx.messages, "notice: $(n_mapped) PowerElectronicsConnection object(s) mapped ($(n_regulating) voltage-regulating, $(n_mapped - n_regulating) PQ; P = $(round(p_sum; digits = 2)) MW, Q = $(round(q_sum; digits = 2)) MVAr injected), $(n_skipped) skipped")
  return nothing
end

# A mapped EnergySource or PowerElectronicsConnection is no slack candidate
# of the delivery, but with P > 0 it is generation, so
# _deEnergizeSourcelessComponents! keeps the island it feeds energized. With multi_slack = true
# _slackForIslandsWithoutReference! gives such an island a slack before this
# runs, so the notice below only fires with multi_slack = false: then an
# island whose only sources are such objects leaves the import without a
# reference. The island-wise solve (power_flow.islands) gives such an island
# its strongest generating unit as reference (detect_ac_islands, status
# promote_generator_ref), so the P of that source is a result of the island
# balance and no longer the SSH value; with islands disabled the solve
# rejects the island. One notice per such island makes that visible. A Slack
# or PV bus in the island (a regulating SVC or PowerElectronicsConnection)
# means the solver picks another reference, so those islands are left out.
function _noticeEnergySourceOnlyIslands!(net, store, ctx::_MapCtx)
  src = Dict{Int,Tuple{Symbol,String}}()   # prosumer index -> (class, name)
  classes = Symbol[]
  ctx.map_energy_source && push!(classes, :EnergySource)
  ctx.map_power_electronics && push!(classes, :PowerElectronicsConnection)
  for cls in classes
    for o in objectsOf(store, cls)
      k = get(ctx.prosumer_by_mrid, o.mrid, nothing)
      k === nothing || (src[k] = (cls, something(str(o, :name), o.mrid)))
    end
  end
  isempty(src) && return nothing
  # generating units per bus that the solver could promote (an SVC carries
  # no active power and is never promoted)
  units = Dict{Int,Vector{Int}}()
  for (k, ps) in enumerate(net.prosumpsVec)
    Sparlectra.isGenerator(ps.proSumptionType) || continue
    ps.comp.cTyp == Sparlectra.StaticVarCompensator && continue
    push!(get!(Vector{Int}, units, Int(ps.comp.cFrom_bus)), k)
  end
  # the islands of the power flow (switch-only buses with their group)
  islands, held, _ = Sparlectra._power_flow_islands(net)
  for (buses, has_ref) in zip(islands, held)
    has_ref && continue
    ks = reduce(vcat, (get(units, b, Int[]) for b in buses); init = Int[])
    (isempty(ks) || !all(k -> haskey(src, k), ks)) && continue
    names = sort!([last(src[k]) for k in ks])
    shown = length(names) <= 5 ? join(names, ", ") : join(names[1:5], ", ") * ", ..."
    classes = join(sort!(unique([String(first(src[k])) for k in ks])), " and ")
    around = something(_busDictName(net, first(buses)), net.nodeVec[first(buses)].comp.cName)
    which = length(ks) == 1 ? "it" : "the strongest of them"
    push!(ctx.messages, "notice: electrical island of $(length(buses)) bus(es) around $(around) has no slack and no source but $(classes) ($(shown)); the island-wise solve makes $(which) the island reference, so its P follows the island balance instead of the SSH value (with power_flow.islands.enabled = false the solve rejects the island)")
  end
  return nothing
end

# CGMES classes with an active power value in the SSH profile (2.4.15 and
# 3.0), with the attribute the reader stores it under. All of them use the
# load sign: positive means flow out of the node. The reader keeps every
# subclass under its own class, so the EnergyConsumer and RotatingMachine
# subclasses are listed one by one; RotatingMachine.p arrives as :p on
# SynchronousMachine and AsynchronousMachine, ACDCConverter.p on VsConverter
# and CsConverter. PowerElectronicsConnection exists in 3.0 only.
const _SSH_ACTIVE_POWER_ATTR = [
  :EnergyConsumer => :p,
  :ConformLoad => :p,
  :NonConformLoad => :p,
  :StationSupply => :p,
  :EnergySource => :activePower,
  :EquivalentInjection => :p,
  :ExternalNetworkInjection => :p,
  :SynchronousMachine => :p,
  :AsynchronousMachine => :p,
  :VsConverter => :p,
  :CsConverter => :p,
  :PowerElectronicsConnection => :p,
]

# CGMES classes the importer turns into network elements (the class table in
# docs/src/cgmes_import.md, the kinds of CGMESImportResult.element_by_mrid).
# EnergySource and PowerElectronicsConnection join them only with their
# switch on, see _unmappedClasses.
const _MAPPED_EQUIPMENT_CLASSES = [
  :ACLineSegment,
  :SeriesCompensator,
  :EquivalentBranch,
  :PowerTransformer,
  :LinearShuntCompensator,
  :NonlinearShuntCompensator,
  :EnergyConsumer,
  :ConformLoad,
  :NonConformLoad,
  :StationSupply,
  :AsynchronousMachine,
  :SynchronousMachine,
  :StaticVarCompensator,
  :ExternalNetworkInjection,
  :EquivalentInjection,
  :VsConverter,
  :CsConverter,
]

# Node equipment: a busbar section or a junction is a point of its
# TopologicalNode, so the topology consumes it.
const _NODE_EQUIPMENT_CLASSES = [:BusbarSection, :Junction]

# The switch family of CGMES 2.4.15 and 3.0. The importer links the classes
# of _SWITCH_CLASSES (and the topology processor merges them); the others
# (Fuse, GroundDisconnector, Cut, DisconnectingCircuitBreaker) it neither
# links nor merges.
const _SWITCH_FAMILY_CLASSES = [_SWITCH_CLASSES..., :Fuse, :GroundDisconnector, :Cut, :DisconnectingCircuitBreaker]

# A switch the topology consumes: open (the importer's rule: SSH `open` over
# EQ `normalOpen`, out of service counts as open), a terminal disconnected, or
# its terminals on fewer than two TopologicalNodes (merged by the delivery's
# TP profile or by the topology processor). A closed switch between two nodes
# that is no link is not consumed: the coupling it makes is missing.
function _switchConsumed(flags::_MapCtx, sw::CIMObject, ts::Vector{CIMObject})::Bool
  _switchIsOpen(flags, sw) && return true
  any(t -> !_conn(flags, terminalConnected(t)), ts) && return true
  return length(unique(tnOfTerminal(t) for t in ts if tnOfTerminal(t) !== nothing)) < 2
end

# The table of classes not mapped in the import report, so that missing
# generation, load or network elements show up there. The candidates are the
# classes of every object outside the boundary set that a Terminal names as
# its ConductingEquipment, plus the classes of _SSH_ACTIVE_POWER_ATTR (they
# count even without a terminal).
#
# An object is mapped when its mRID is a key of the element map, when a key
# `<mRID>#<part>` derived from it is (the charging shunt of a line on one
# node), or when it is the PowerTransformer of a mapped PowerTransformerEnd
# (a three-winding leg). The parts of a mapped object (transformer ends, tap
# changers, the units behind a generator) and the topology itself
# (ConnectivityNode, TopologicalNode) carry no terminal of their own and are
# never candidates.
#
# A class gets a row when none of its objects is mapped and some object is
# not accounted for otherwise. Accounted for are
# - node equipment (_NODE_EQUIPMENT_CLASSES), consumed by the topology;
# - a switch of _SWITCH_FAMILY_CLASSES the topology consumes (_switchConsumed);
# - an object out of service or without a connected terminal on a
#   TopologicalNode, if its class is one the importer maps
#   (_MAPPED_EQUIPMENT_CLASSES; EnergySource and PowerElectronicsConnection
#   with their switch on): the importer skips it on purpose and says so.
# Objects of a class the importer does not map count in any state.
# The row holds the number of objects not accounted for, the sum of their SSH
# active power in MW in load sign (`nothing` where none of them carries an
# SSH value) and how many of them carry none. Reads only the store and the
# element map, so the importer and the tests share it.
function _unmappedClasses(store::CGMESStore, element_by_mrid::AbstractDict; map_energy_source::Bool = false, map_power_electronics::Bool = false, ignore_connected::Bool = false)
  # carries the importer's in-service and connected rules with its override
  flags = _MapCtx(false, ignore_connected)
  terms = Dict{String,Vector{CIMObject}}()
  for t in objectsOf(store, :Terminal)
    ce = get(t.refs, :ConductingEquipment, nothing)
    (ce === nothing || ce in store.boundary) && continue
    push!(get!(Vector{CIMObject}, terms, ce), t)
  end
  mapped = Set{String}()
  for key in keys(element_by_mrid)
    m = haskey(store.objects, key) ? String(key) : String(first(split(key, '#'; limit = 2)))
    push!(mapped, m)
    o = get(store.objects, m, nothing)
    if o !== nothing && o.class === :PowerTransformerEnd
      pt = get(o.refs, :PowerTransformer, nothing)
      pt === nothing || push!(mapped, pt)
    end
  end
  candidates = Set{Symbol}(first(p) for p in _SSH_ACTIVE_POWER_ATTR)
  for ce in keys(terms)
    o = get(store.objects, ce, nothing)
    o === nothing || push!(candidates, o.class)
  end
  with_path = Set{Symbol}(_MAPPED_EQUIPMENT_CLASSES)
  map_energy_source && push!(with_path, :EnergySource)
  map_power_electronics && push!(with_path, :PowerElectronicsConnection)
  no_terms = CIMObject[]
  rows = @NamedTuple{class::Symbol, count::Int, p_mw::Union{Nothing,Float64}, no_p::Int}[]
  for cls in sort!(collect(candidates); by = String)
    cls in _NODE_EQUIPMENT_CLASSES && continue
    objs = [o for o in objectsOf(store, cls) if !(o.mrid in store.boundary)]
    (isempty(objs) || any(o -> o.mrid in mapped, objs)) && continue
    if cls in _SWITCH_FAMILY_CLASSES
      filter!(o -> !_switchConsumed(flags, o, get(terms, o.mrid, no_terms)), objs)
    elseif cls in with_path
      filter!(o -> _inService(flags, o) && any(t -> _conn(flags, terminalConnected(t)) && tnOfTerminal(t) !== nothing, get(terms, o.mrid, no_terms)), objs)
    end
    isempty(objs) && continue
    k = findfirst(p -> first(p) === cls, _SSH_ACTIVE_POWER_ATTR)
    vals = Union{Nothing,Float64}[k === nothing ? nothing : num(o, last(_SSH_ACTIVE_POWER_ATTR[k])) for o in objs]
    no_p = count(isnothing, vals)
    p_mw = no_p == length(vals) ? nothing : sum(v -> something(v, 0.0), vals)
    push!(rows, (class = cls, count = length(objs), p_mw = p_mw, no_p = no_p))
  end
  return rows
end

# RegulatingCondEq.controlEnabled (SSH): false switches the unit's voltage
# regulation off, whatever its RegulatingControl says. Absent counts as true
# (Sparlectra's own export and older deliveries do not write it). Read for
# every class the importer maps with voltage regulation: SynchronousMachine,
# ExternalNetworkInjection, StaticVarCompensator. An EquivalentInjection is
# no RegulatingCondEq and has no such flag; it regulates only with
# regulationCapability (EQ) and regulationStatus (SSH) both true.
_regulationEnabled(obj::CIMObject)::Bool = something(boolval(obj, :controlEnabled), true)

# one notice per class naming the units with controlEnabled = false
function _noteRegulationOff!(messages::Vector{String}, ctrl_off, cls::Symbol)
  names = [u.name for u in ctrl_off if u.class == cls]
  isempty(names) && return nothing
  shown = length(names) <= 5 ? join(names, ", ") : join(names[1:5], ", ") * ", ..."
  push!(messages, "notice: $(length(names)) $(cls)(s) with controlEnabled = false, no voltage regulation: $(shown)")
  return nothing
end

# IdentifiedObject.description of an EquivalentInjection the 2.4.15 export
# wrote for a PowerElectronicsConnection (2.4.15 has no class for it). The
# import gives such an injection the component type
# PowerElectronicsConnection again: the short circuit then flags its island
# for the fault current the export could not carry, and a new export writes
# it as EquivalentInjection again.
const _PEC_EI_DESCRIPTION = "PowerElectronicsConnection"

function _mapInjections!(net, store, topo, created, svmap, ctx::_MapCtx; multi_slack::Bool = true)
  scan = _SlackScan()
  curve_pts = _reactiveCapabilityCurvePoints(store)
  angle_refs = _svAngleRefNodes(store)
  # bus name -> name of the SV TopologicalIsland whose angle reference it is
  angle_ref_island = Dict{String,String}(topo.bus_name[tn] => isl for (tn, isl) in angle_refs if haskey(topo.bus_name, tn))
  angle_ref_buses = Set{String}(keys(angle_ref_island))
  isempty(angle_ref_buses) || push!(ctx.messages, "notice: SV declares $(length(angle_ref_buses)) angle-reference node(s), used as slack candidate(s): $(join(sort(collect(angle_ref_buses)), ", "))")
  pending = NamedTuple[]

  # mapped units with controlEnabled = false (machines, network injections,
  # SVCs): names for one notice per class, mRIDs for the slack notice in the
  # flush below
  ctrl_off = @NamedTuple{class::Symbol, mrid::String, name::String}[]
  for sm in objectsOf(store, :SynchronousMachine)
    t = busOfEquipment(topo, sm, 1)
    name = something(str(sm, :name), sm.mrid)
    if t === nothing
      push!(ctx.messages, "skip: SynchronousMachine $(name), no TopologicalNode")
      push!(ctx.skipped, sm.mrid)
      continue
    end
    if !_conn(ctx, t.connected)
      push!(ctx.messages, "skip: SynchronousMachine $(name), terminal disconnected")
      push!(ctx.skipped, sm.mrid)
      continue
    end
    # ReliCapGrid parks out-of-service units exactly like this: connected=true,
    # inService=false, no SvVoltage, placeholder regulation target. Importing
    # them added 413 MW of phantom generation in Svedala alone.
    if !_inService(ctx, sm)
      push!(ctx.messages, "skip: SynchronousMachine $(name), out of service (SSH inService=false)")
      push!(ctx.skipped, sm.mrid)
      continue
    end
    bus = _ensureBus!(net, created, topo, svmap, t.tn)
    gu = ref(store, sm, :GeneratingUnit)
    # P limits of the unit, a placeholder of 9999 MW and above read as none
    # (_pLimit). The Q pair is kept as delivered even when it is a
    # placeholder: a machine without minQ/maxQ gets the band +-max(ratedS,
    # baseMVA) below, which can bind, while +-99999.9 is an open band as
    # stated.
    gu_pmin, gu_pmax = gu === nothing ? (nothing, nothing) : _pLimitsRead!(ctx, "GeneratingUnit", num(gu, :minOperatingP), num(gu, :maxOperatingP))
    ratedS = num(sm, :ratedS, 0.0)
    # GeneratingUnit.normalPF → distributed-slack participation weight
    # (issue #192). Same storage rule as the MATPOWER APF path: only positive
    # finite values count, and carrying the value never changes results while
    # power_flow.distributed_slack is disabled.
    npf = gu === nothing ? nothing : num(gu, :normalPF)
    participation = (npf !== nothing && isfinite(npf) && npf > 0.0) ? npf : nothing
    # Q-limit source, most faithful first:
    #  1. InitialReactiveCapabilityCurve evaluated at the machine's scheduled P
    #     (issue #294 point 1, the first data-faithful limits; MicroGrid
    #     BE-G1 writes the degenerate 0/0 scalar pair precisely because its
    #     real limits live in the curve),
    #  2. the scalar minQ/maxQ hull,
    #  3. wide symmetric limits so voltage regulation is not spuriously
    #     Q-pinned when neither source is usable.
    qwide = max(ratedS, net.baseMVA)
    qlo = nothing
    qhi = nothing
    curve = get(sm.refs, :InitialReactiveCapabilityCurve, nothing)
    if curve !== nothing && haskey(curve_pts, curve)
      cq = _curveQHull(curve_pts[curve], num(sm, :p, 0.0))
      if cq !== nothing
        qlo, qhi = cq
        push!(ctx.messages, "notice: SynchronousMachine $(name), Q limits [$(round(qlo; digits = 1)), $(round(qhi; digits = 1))] MVAr from its ReactiveCapabilityCurve at P = $(num(sm, :p, 0.0)) MW")
      end
    end
    if qlo === nothing
      qlo = _limitHullMin(num(sm, :minQ), num(sm, :maxQ))
      qhi = _limitHullMax(num(sm, :minQ), num(sm, :maxQ))
    end
    if qlo === nothing || qhi === nothing || (qhi - qlo) < 1e-6
      qlo = -qwide
      qhi = qwide
      push!(ctx.messages, "notice: SynchronousMachine $(name) has no usable scalar minQ/maxQ, using ±$(qwide) MVAr")
    end
    prio = num(sm, :referencePriority, 0.0)
    # the rating ranks the slack candidates; a placeholder rating (9999 MVA
    # and above, the rule of the short circuit) ranks like no rating, so a
    # unit without data does not become the largest machine of the delivery
    rank_s = ratedS >= Sparlectra.SC_RATEDS_PLACEHOLDER_MVA ? 0.0 : ratedS
    # an SV angle reference outranks referencePriority: it is what the
    # exporting tool actually used
    bus in angle_ref_buses && push!(scan.sv_refs, (0, -rank_s, bus))
    prio >= 1.0 && push!(scan.by_priority, (prio, -rank_s, bus))
    push!(scan.machines, (rank_s, bus, name))
    # RegulatingCondEq.controlEnabled (SSH) switches the machine's voltage
    # regulation off: no PV bus, no remote plan, a PQ injection at the SSH
    # p/q with the Q limits above (_regulationEnabled). The flag concerns the
    # regulating control only, so the slack bookkeeping above stays: a
    # machine chosen as slack holds its bus at the SV/start voltage (slack
    # branch where the pending tuples are flushed).
    ctrl_enabled = _regulationEnabled(sm)
    ctrl_enabled || push!(ctrl_off, (class = :SynchronousMachine, mrid = sm.mrid, name = name))
    vset, vctrl, remote_bus = if ctrl_enabled
      _voltageSetpoint(store, topo, sm, t.tn, ctx.messages; vset_min_pu = ctx.vset_min_pu, vset_max_pu = ctx.vset_max_pu, allow_remote = ctx.machine_control)
    else
      (nothing, false, nothing)
    end
    sm_rc = ref(store, sm, :RegulatingControl)
    # A resolvable remote control (machine_control = true) keeps the machine
    # a PQ injection with the SSH operating point; the outer-loop controller
    # planned here moves its Q until the remote bus hits the target. Whether
    # the plan survives is decided in _resolveMachineControlPlans! once all
    # regulating units and the slack selection are known.
    is_remote = remote_bus !== nothing
    push!(
      pending,
      (
        bus = bus,
        type = "SYNCHRONOUSMACHINE",
        participationFactor = participation,
        # the unit's own referencePriority travels onto its prosumer (0.30.2)
        # whatever the slack selection below makes of it
        referencePriority = _cgmesReferencePriority(prio),
        p = -num(sm, :p, 0.0),                       # CGMES machine sign: p < 0 = injection
        q = -num(sm, :q, 0.0),
        pMin = gu_pmin,
        pMax = gu_pmax,
        # The ENTSO-E test sets are inconsistent about the sign convention of
        # machine minQ/maxQ (MicroGrid NL: load convention like p/q, the SV
        # q lies in [-maxQ, -minQ]; SmallGrid: injection convention). Stage 1
        # uses the hull of both readings so no unit is wrongly Q-pinned;
        # strict limit fidelity is a Stage-2 topic.
        qMin = qlo,
        qMax = qhi,
        vm_pu = (vctrl && !is_remote) ? _pvVoltage(net, bus, vset) : nothing,
        isRegulated = vctrl && !is_remote,
        name = name,
        remote_target = is_remote ? (bus = remote_bus, vset_pu = vset) : nothing,
        mrid = sm.mrid,
        terminal_mrid = _firstTerminalMrid(topo, sm),
        gu_mrid = gu === nothing ? nothing : gu.mrid,
        rc_mrid = sm_rc === nothing ? nothing : sm_rc.mrid,
      ),
    )
  end
  _noteRegulationOff!(ctx.messages, ctrl_off, :SynchronousMachine)

  # AsynchronousMachine (#294 point 6), Stage-0 model: a fixed PQ operating
  # point from SSH (`RotatingMachine.p`/`q`, load convention, a motor
  # consumes with p > 0, a generator-kind machine carries p < 0). No voltage
  # regulation and no slack candidacy: the induction machine's voltage/slip
  # dependence is a dynamics topic, not a steady-state power-flow one.
  # Fixture: FullGrid (ASM_1, ≈1 MW motor).
  for am in objectsOf(store, :AsynchronousMachine)
    t = busOfEquipment(topo, am, 1)
    name = something(str(am, :name), am.mrid)
    if t === nothing
      push!(ctx.messages, "skip: AsynchronousMachine $(name), no TopologicalNode")
      push!(ctx.skipped, am.mrid)
      continue
    end
    if !_conn(ctx, t.connected)
      push!(ctx.messages, "skip: AsynchronousMachine $(name), terminal disconnected")
      push!(ctx.skipped, am.mrid)
      continue
    end
    if !_inService(ctx, am)
      push!(ctx.messages, "skip: AsynchronousMachine $(name), out of service (SSH inService=false)")
      push!(ctx.skipped, am.mrid)
      continue
    end
    bus = _ensureBus!(net, created, topo, svmap, t.tn)
    p = num(am, :p, 0.0)
    q = num(am, :q, 0.0)
    Sparlectra.addProsumer!(net = net, busName = bus, type = "ASYNCHRONOUSMACHINE", p = p, q = q, defer_bus_type_refresh = true)
    # a consumption-type prosumer (motor convention), hence a load
    ctx.element_by_mrid[am.mrid] = (:load, length(net.prosumpsVec))
    # The prosumer constructor collapses the type string to a generic Load
    # component; keep the true component type on the comp so the
    # short-circuit evaluation can recognize motors (skip + lower-bound
    # flag).
    net.prosumpsVec[end].comp.cTyp = Sparlectra.AsynchronousMachine
    _captureBusEquipment!(net, ctx, "ASM", bus, am.mrid; children = (T1 = _firstTerminalMrid(topo, am),))
    # No internal staging vocabulary in user-facing importer messages.
    push!(ctx.messages, "notice: AsynchronousMachine $(name), fixed PQ operating point from SSH (p=$(round(p; digits = 2)) MW, q=$(round(q; digits = 2)) MVAr, load convention) at $(bus)")
  end

  for eni in objectsOf(store, :ExternalNetworkInjection)
    t = busOfEquipment(topo, eni, 1)
    name = something(str(eni, :name), eni.mrid)
    if t === nothing || !_conn(ctx, t.connected)
      push!(ctx.messages, "skip: ExternalNetworkInjection $(name), missing/disconnected terminal")
      push!(ctx.skipped, eni.mrid)
      continue
    end
    # must run before the slack-candidate bookkeeping below: an out-of-service
    # network equivalent is no reference candidate
    if !_inService(ctx, eni)
      push!(ctx.messages, "skip: ExternalNetworkInjection $(name), out of service (SSH inService=false)")
      push!(ctx.skipped, eni.mrid)
      continue
    end
    bus = _ensureBus!(net, created, topo, svmap, t.tn)
    prio = num(eni, :referencePriority, 0.0)
    # the size for the ranking: maxP, a placeholder (_pLimit) read as none
    maxP = abs(something(_pLimit(num(eni, :maxP)), num(eni, :p, 0.0)))
    bus in angle_ref_buses && push!(scan.sv_refs, (0, -maxP, bus))
    prio >= 1.0 && push!(scan.by_priority, (prio, -maxP, bus))
    push!(scan.enis, (maxP, bus))
    # controlEnabled = false: a PQ injection at the SSH p/q, as for the
    # machines; the slack candidacy above stays, as slack it holds the
    # SV/start voltage of its bus instead of its control target
    ctrl_enabled = _regulationEnabled(eni)
    ctrl_enabled || push!(ctrl_off, (class = :ExternalNetworkInjection, mrid = eni.mrid, name = name))
    vset, vctrl, _ = ctrl_enabled ? _voltageSetpoint(store, topo, eni, t.tn, ctx.messages; vset_min_pu = ctx.vset_min_pu, vset_max_pu = ctx.vset_max_pu) : (nothing, false, nothing)
    eni_rc = ref(store, eni, :RegulatingControl)
    # P limits from EQ minP/maxP as delivered, in generator sign, as for the
    # EquivalentInjection below: the CIM names them limits "of the
    # injection", PowSyBl takes them unchanged as generator minP/maxP and
    # writes them back so, and our export writes the prosumer pMin/pMax the
    # same way. The hull would widen an asymmetric range such as (-5, 15)
    # to (-15, 15). A missing, degenerate (0/0) or inverted pair gives none,
    # and so does a pair with a placeholder of 9999 MW and above (_pLimit).
    plo, phi = _pLimitsRead!(ctx, "ExternalNetworkInjection", num(eni, :minP), num(eni, :maxP))
    (plo === nothing || phi === nothing || (phi - plo) < 1e-6) && (plo = phi = nothing)
    push!(
      pending,
      (
        bus = bus,
        type = "EXTERNALNETWORKINJECTION",
        referencePriority = _cgmesReferencePriority(prio),
        p = -num(eni, :p, 0.0),                      # load convention → injection
        q = -num(eni, :q, 0.0),
        pMin = plo,
        pMax = phi,
        qMin = _limitHullMin(num(eni, :minQ), num(eni, :maxQ)),
        qMax = _limitHullMax(num(eni, :minQ), num(eni, :maxQ)),
        vm_pu = vctrl ? _pvVoltage(net, bus, vset) : nothing,
        isRegulated = vctrl,
        mrid = eni.mrid,
        terminal_mrid = _firstTerminalMrid(topo, eni),
        rc_mrid = eni_rc === nothing ? nothing : eni_rc.mrid,
      ),
    )
  end
  _noteRegulationOff!(ctx.messages, ctrl_off, :ExternalNetworkInjection)

  # EquivalentInjections per boundary node: an assembled model carries one EI
  # per SIDE on each X-node (e.g. BE-Inj-X… and NL-Inj-X…). Two or more EIs
  # on one boundary node therefore mean both sides are present and the
  # equivalents must not be applied on top of the real flows. A dangling
  # line keeps its single EI. (A count of real terminals is NOT a valid
  # criterion: SmallGrid has X-nodes with two parallel dangling circuits
  # and one EI.)
  #
  # Node identity is version-dependent. A CGMES 2.4.15 boundary set ships
  # TopologicalNodes (TP_BD), so both sides' EIs meet on one boundary TN. A
  # CGMES 3.0 boundary carries only ConnectivityNodes: each area's own TP
  # assigns its own TN to the shared CN, and the two sides need not even agree
  # on one TN (ReliCapGrid's Espheim-Portheim border: two TNs, one CN). The
  # shared identity is therefore the terminal's boundary CN when there is one,
  # the boundary TN otherwise.
  #
  # The p values are collected (not just counted) because presence of both
  # sides alone is not sufficient to discard: an assembled AC X-node pair is
  # two views of the SAME exchange and cancels exactly (measured: MicroGrid
  # 0.0 %, ReliCapGrid AC borders 0.0-1.6 % of max|p|), while the pair at a
  # DC crossing consists of the two CONVERTER injections whose difference is
  # the DC loss (ReliCapGrid BP_SD-EH_DC1: +94.57 / −72.60, 23.2 %). Those
  # must survive: they are devices, not a double-counted equivalent.
  ei_p_bnode = Dict{String,Vector{Float64}}()
  # StaticVarCompensator: a controllable reactive source with P = 0, the
  # MATPOWER-compatible model (generator with Pg = 0). CGMES capacitiveRating /
  # inductiveRating are REACTANCES in Ohm (RealGrid: ±202.5 Ω @ 225 kV ↔
  # 250/-100 MVAr), so the Q limits come from Q = vn²/X. The droop `slope` is
  # deliberately ignored (standard practice), with a notice.
  for svc in objectsOf(store, :StaticVarCompensator)
    t = busOfEquipment(topo, svc, 1)
    name = something(str(svc, :name), svc.mrid)
    if t === nothing || !_conn(ctx, t.connected)
      push!(ctx.messages, "skip: StaticVarCompensator $(name), missing/disconnected terminal")
      push!(ctx.skipped, svc.mrid)
      continue
    end
    if !_inService(ctx, svc)
      push!(ctx.messages, "skip: StaticVarCompensator $(name), out of service (SSH inService=false)")
      push!(ctx.skipped, svc.mrid)
      continue
    end
    bus = _ensureBus!(net, created, topo, svmap, t.tn)
    qwide = 10.0 * net.baseMVA
    xcap = num(svc, :capacitiveRating)
    xind = num(svc, :inductiveRating)
    q_from_x(x) = (x === nothing || x == 0.0) ? nothing : t.vn_kV^2 / x
    qcap = q_from_x(xcap)
    qind = q_from_x(xind)
    qlo = qcap === nothing || qind === nothing ? -qwide : min(qcap, qind)
    qhi = qcap === nothing || qind === nothing ? qwide : max(qcap, qind)
    if max(abs(qlo), abs(qhi)) > qwide
      push!(ctx.messages, "notice: StaticVarCompensator $(name), implausible Q range $(round(qlo; digits=1))..$(round(qhi; digits=1)) MVAr from ratings, clamped to ±$(qwide)")
      qlo = max(qlo, -qwide)
      qhi = min(qhi, qwide)
    end
    slope = num(svc, :slope, 0.0)
    slope != 0.0 && push!(ctx.messages, "notice: StaticVarCompensator $(name), droop slope $(slope) ignored (constant-voltage model)")
    # controlEnabled = false: PQ at the SSH q, both for a voltage
    # RegulatingControl and for the sVCControlMode/voltageSetPoint fallback
    ctrl_enabled = _regulationEnabled(svc)
    ctrl_enabled || push!(ctrl_off, (class = :StaticVarCompensator, mrid = svc.mrid, name = name))
    vset, vctrl, _ = ctrl_enabled ? _voltageSetpoint(store, topo, svc, t.tn, ctx.messages; vset_min_pu = ctx.vset_min_pu, vset_max_pu = ctx.vset_max_pu) : (nothing, false, nothing)
    if ctrl_enabled && !vctrl && endswith(something(enumval(svc, :sVCControlMode), ""), ".voltage")
      vsp = num(svc, :voltageSetPoint)
      (vsp !== nothing && vsp > 0.0) && ((vset, vctrl) = (vsp / t.vn_kV, true))
    end
    # Consistency guard against snapshots whose control flags contradict their
    # own solved state (RealGrid: RC enabled with target 1.0 pu while the SV
    # voltages sit at 1.008-1.031 and SSH q = 0, the reference clearly did
    # not let the SVCs regulate; importing them as PV at 1.0 fights the whole
    # 225 kV profile and breaks convergence). If the SV start voltage of the
    # bus deviates from the target by more than 0.005 pu, fall back to PQ
    # with the SSH operating point.
    if vctrl && vset !== nothing
      node = net.nodeVec[Sparlectra.geNetBusIdx(net = net, busName = bus)]
      startvm = node._vm_pu
      if startvm !== nothing && abs(startvm - vset) > 0.005
        push!(ctx.messages, "notice: StaticVarCompensator $(name), voltage target $(round(vset; digits=4)) pu contradicts the SV state ($(round(startvm; digits=4)) pu); imported as PQ with the SSH operating point")
        vctrl = false
      end
    end
    svc_rc = ref(store, svc, :RegulatingControl)
    push!(
      pending,
      (
        bus = bus,
        type = "STATICVARCOMPENSATOR",
        p = 0.0,
        q = -num(svc, :q, 0.0),                      # load convention → injection
        pMin = 0.0,
        pMax = 0.0,
        qMin = qlo,
        qMax = qhi,
        vm_pu = vctrl ? _pvVoltage(net, bus, vset) : nothing,
        isRegulated = vctrl,
        mrid = svc.mrid,
        terminal_mrid = _firstTerminalMrid(topo, svc),
        rc_mrid = svc_rc === nothing ? nothing : svc_rc.mrid,
      ),
    )
  end
  _noteRegulationOff!(ctx.messages, ctrl_off, :StaticVarCompensator)

  # HVDC converter stations, Stage-0 model: each converter becomes a fixed
  # injection at its AC connection point with the SSH operating point, the
  # industry-standard load-flow treatment of a point-to-point link (MATPOWER's
  # dcline does the same with two bounded dummy generators). The one property
  # that matters for the power flow is that HVDC has NO angle coupling: the
  # transfer is a control setpoint, not the result of an angle difference. The
  # DC middle (DCLineSegment/DCNode) is deliberately NOT mapped, so two areas
  # joined only through HVDC correctly remain separate electrical islands with
  # their own angle references. The loss shows up as the difference between
  # the two converters' setpoints (ReliCapGrid GA-BH: −100.0 / +109.1 MW).
  hvdc_derived_p = _hvdcDerivedSlackPower(store, ctx)
  for cls in (:VsConverter, :CsConverter)
    for conv in objectsOf(store, cls)
      t = busOfEquipment(topo, conv, 1)
      name = something(str(conv, :name), conv.mrid)
      if t === nothing || !_conn(ctx, t.connected)
        push!(ctx.messages, "skip: $(cls) $(name), missing/disconnected AC terminal")
        push!(ctx.skipped, conv.mrid)
        continue
      end
      if !_inService(ctx, conv)
        push!(ctx.messages, "skip: $(cls) $(name), out of service (SSH inService=false)")
        push!(ctx.skipped, conv.mrid)
        continue
      end
      # a DC-voltage-controlling converter reports its power as a result, not
      # as a schedule; _hvdcDerivedSlackPower supplies the derived value for
      # exactly those ends and is empty for every ordinary converter
      p_ssh = get(hvdc_derived_p, conv.mrid, num(conv, :p, 0.0))
      q_ssh = num(conv, :q, 0.0)
      bus = _ensureBus!(net, created, topo, svmap, t.tn)
      push!(
        pending,
        (
          bus = bus,
          type = "EXTERNALNETWORKINJECTION",
          p = -p_ssh,                              # load convention → injection
          q = -q_ssh,
          pMin = -p_ssh,
          pMax = -p_ssh,
          qMin = -q_ssh,
          qMax = -q_ssh,
          vm_pu = nothing,
          isRegulated = false,
          # the mRID lets the HVDC pair attachment find this exact prosumer
          mrid = conv.mrid,
          # a fixed PCC injection, never the angle reference (_isReferenceUnit)
          hvdc_converter = true,
        ),
      )
      push!(ctx.messages, "notice: $(cls) $(name), HVDC converter: fixed PCC injection at $(bus) (p=$(round(-p_ssh; digits = 1)) MW, q=$(round(-q_ssh; digits = 1)) MVAr from SSH; DC side not mapped, no angle coupling)")
    end
  end

  # Collect on the raw terminals, not via busOfEquipment: an EI whose TN does
  # not resolve (3.0 split-TN border) must still be counted, otherwise its
  # counterpart on the other side sees a count of 1 and survives.
  for ei in objectsOf(store, :EquivalentInjection)
    # an out-of-service EI is no evidence that its side is present
    _inService(ctx, ei) || continue
    bnode = _eiBoundaryNode(store, topo, ei)
    bnode === nothing && continue
    push!(get!(Vector{Float64}, ei_p_bnode, bnode), num(ei, :p, 0.0))
  end

  # mapped EquivalentInjections with regulationCapability = true and no SSH
  # regulationStatus: counted for one notice below
  n_ei_status_missing = 0
  # EquivalentInjections imported as sources (equivalent_injection_as_source),
  # and those of them without a usable internal impedance (x missing or not
  # positive), for one notice each below
  ei_sources = String[]
  ei_sources_no_z = String[]
  for ei in objectsOf(store, :EquivalentInjection)
    t = busOfEquipment(topo, ei, 1)
    name = something(str(ei, :name), ei.mrid)
    if t === nothing
      push!(ctx.messages, "skip: EquivalentInjection $(name), no TopologicalNode")
      push!(ctx.skipped, ei.mrid)
      continue
    end
    if !_conn(ctx, t.connected)
      push!(ctx.messages, "skip: EquivalentInjection $(name), terminal disconnected")
      push!(ctx.skipped, ei.mrid)
      continue
    end
    if !_inService(ctx, ei)
      push!(ctx.messages, "skip: EquivalentInjection $(name), out of service (SSH inService=false)")
      push!(ctx.skipped, ei.mrid)
      continue
    end
    bnode = _eiBoundaryNode(store, topo, ei)
    if bnode !== nothing && _eiPairCancels(get(ei_p_bnode, bnode, Float64[]))
      push!(ctx.messages, "skip: EquivalentInjection $(name), assembled boundary node $(topo.bus_name[t.tn]) has both sides present")
      push!(ctx.skipped, ei.mrid)
      continue
    end
    bus = _busForTerminal!(net, created, topo, svmap, ctx, ei, 1, t.tn)
    # regulationCapability (EQ) and regulationStatus (SSH) both true plus a
    # positive target → PV-like boundary. An EquivalentInjection is no
    # RegulatingCondEq: it has neither controlEnabled nor a
    # RegulatingControl. Either flag absent counts as false; a capable unit
    # without a status (SSH missing or incomplete) is counted for the notice.
    vset = nothing
    ei_capable = something(boolval(ei, :regulationCapability), false)
    ei_status = boolval(ei, :regulationStatus)
    ei_capable && ei_status === nothing && (n_ei_status_missing += 1)
    if ei_capable && something(ei_status, false)
      target = num(ei, :regulationTarget, 0.0)
      target > 0.0 && (vset = target / t.vn_kV)
    end
    # Q limits in the order of the machines: the ReactiveCapabilityCurve at
    # the scheduled p, else the scalar minQ/maxQ hull. Unlike a machine an
    # EquivalentInjection has no rating for a wide fallback, so a missing
    # or degenerate pair (0/0) leaves it without Q limits.
    qlo = nothing
    qhi = nothing
    ei_curve = get(ei.refs, :ReactiveCapabilityCurve, nothing)
    if ei_curve !== nothing && haskey(curve_pts, ei_curve)
      cq = _curveQHull(curve_pts[ei_curve], num(ei, :p, 0.0))
      if cq !== nothing
        qlo, qhi = cq
        push!(ctx.messages, "notice: EquivalentInjection $(name): Q limits [$(round(qlo; digits = 1)), $(round(qhi; digits = 1))] MVAr from its ReactiveCapabilityCurve at P = $(num(ei, :p, 0.0)) MW")
      end
    end
    if qlo === nothing
      qlo = _limitHullMin(num(ei, :minQ), num(ei, :maxQ))
      qhi = _limitHullMax(num(ei, :minQ), num(ei, :maxQ))
      (qlo === nothing || qhi === nothing || (qhi - qlo) < 1e-6) && (qlo = qhi = nothing)
    end
    # P limits from EQ minP/maxP as delivered, in generator sign: the CIM
    # names them limits "of the injection" (p in SSH is load sign), PowSyBl
    # takes them unchanged as generator minP/maxP and writes them back so,
    # and the export of a converter or a network injection writes its
    # prosumer pMin/pMax the same way. Unlike the scalar Q pair (hull, the
    # test sets disagree on its sign) the hull would widen an asymmetric
    # range such as (-5, 15) to (-15, 15). A missing, degenerate (0/0) or
    # inverted pair gives none, and so does a pair with a placeholder of
    # 9999 MW and above (_pLimit). Distributed slack (p_mode pmax/headroom,
    # the P-limit check) and the slack ranking read them.
    plo, phi = _pLimitsRead!(ctx, "EquivalentInjection", num(ei, :minP), num(ei, :maxP))
    (plo === nothing || phi === nothing || (phi - plo) < 1e-6) && (plo = phi = nothing)
    # An EquivalentInjection on the SV angle reference node can hold the
    # reference (MicroGrid BE: the X-node TN_Border_ST23 carries only its
    # EquivalentInjection). It ranks behind a machine or a network injection
    # on an SV angle reference, so the slack of every delivery the importer
    # solved before stays where it was.
    bus in angle_ref_buses && push!(scan.sv_refs, (1, -abs(something(phi, num(ei, :p, 0.0))), bus))
    # equivalent_injection_as_source: the injection stands for the grid
    # behind a boundary node (or the reference node) and becomes a source
    # that holds the SV voltage of its node (PV at its SSH P); the flush
    # gives it the SV angle as reference angle (u_ref_angle in an SCF
    # export), which only the slack among them holds in the solve. Its Q
    # limits stay as read above.
    as_source = ctx.ei_as_source && (bus in angle_ref_buses || _eiOnBoundary(store, topo, ei))
    if as_source
      push!(ei_sources, name)
      x_ei = num(ei, :x)
      (x_ei === nothing || !(x_ei > 0.0)) && push!(ei_sources_no_z, name)
    end
    push!(
      pending,
      (
        bus = bus,
        type = "EXTERNALNETWORKINJECTION",
        p = -num(ei, :p, 0.0),
        q = -num(ei, :q, 0.0),
        pMin = plo,
        pMax = phi,
        qMin = qlo,
        qMax = qhi,
        vm_pu = as_source ? _pvVoltage(net, bus, nothing) : vset,
        isRegulated = as_source || vset !== nothing,
        source = as_source,
        # no `mrid` on purpose (export identity, see the flush below); the
        # source mRID only feeds element_by_mrid
        src_mrid = ei.mrid,
        converter = str(ei, :description) == _PEC_EI_DESCRIPTION,
      ),
    )
  end
  n_ei_status_missing > 0 && push!(ctx.messages, "notice: $(n_ei_status_missing) EquivalentInjection(s) with regulationCapability = true and no SSH regulationStatus, no voltage regulation")
  isempty(ei_sources) || push!(ctx.messages, "notice: $(length(ei_sources)) EquivalentInjection(s) imported as source at the SV voltage of their node (cgmes_import.equivalent_injection_as_source): $(_shortList(ei_sources))")
  isempty(ei_sources_no_z) || push!(ctx.messages, "notice: $(length(ei_sources_no_z)) EquivalentInjection(s) imported as source without internal impedance (r/x missing or x not positive, no short-circuit power): $(_shortList(ei_sources_no_z))")

  _mapEnergySources!(pending, net, store, topo, created, svmap, ctx)
  _mapPowerElectronicsConnections!(pending, net, store, topo, created, svmap, ctx)

  # island names for the slack messages: the SV TopologicalIsland of an
  # angle reference, else the electrical island number of the bus (computed
  # once, only when a message needs it)
  island_of = Ref{Union{Nothing,Dict{Int,Int}}}(nothing)
  island_label = function (bus::String)
    haskey(angle_ref_island, bus) && return "island $(angle_ref_island[bus])"
    island_of[] === nothing && (island_of[] = Sparlectra.electricalIslandOfBus(net))
    return "electrical island $(get(island_of[], get(net.busDict, bus, 0), 0))"
  end
  slack_bus = _chooseSlack(scan, ctx.messages, island_label, !isempty(angle_ref_buses))
  # One slack per island: every SV angle-reference bus that carries a unit
  # becomes a reference, not just the primary one. A multi-island delivery
  # (ReliCapGrid/Svedala: two islands) is otherwise unsolvable, each island
  # needs its own angle reference.
  slack_buses = _selectSlackBuses(net, ctx, slack_bus, angle_ref_buses, pending, multi_slack, island_label)
  ctx.machine_control && _resolveMachineControlPlans!(pending, net, slack_buses, ctx.messages)
  marked = Set{String}()
  for inj in pending
    # only a machine, a network injection or an EquivalentInjection takes
    # the reference role (_isReferenceUnit), never an EnergySource, a
    # PowerElectronicsConnection, an SVC or an HVDC converter on the slack bus
    is_slack = _isReferenceUnit(inj) && inj.bus in slack_buses && !(inj.bus in marked)
    is_slack && push!(marked, inj.bus)
    # the slack and a source hold the start angle of their bus (the SV angle
    # where present) as reference angle, so an SCF export (u_ref_angle) and
    # a conversion to an external grid keep it
    va = (is_slack || get(inj, :source, false)) ? something(net.nodeVec[Sparlectra.geNetBusIdx(net = net, busName = inj.bus)]._va_deg, 0.0) : nothing
    vm = inj.vm_pu
    if is_slack && vm === nothing
      # slack needs a voltage target: use the SV/start value of its bus
      node = net.nodeVec[Sparlectra.geNetBusIdx(net = net, busName = inj.bus)]
      vm = something(node._vm_pu, 1.0)
      # a unit with controlEnabled = false lands here too: its control
      # target does not apply, the start voltage does
      off_i = findfirst(u -> u.mrid == get(inj, :mrid, nothing), ctrl_off)
      off_i === nothing || push!(ctx.messages, "notice: slack $(ctrl_off[off_i].class) $(ctrl_off[off_i].name) has controlEnabled = false and holds its bus at the start voltage $(round(vm; digits = 4)) pu (SV where present), not at its RegulatingControl target")
    end
    Sparlectra.addProsumer!(
      net = net,
      busName = inj.bus,
      type = inj.type,
      p = inj.p,
      q = inj.q,
      pMin = inj.pMin,
      pMax = inj.pMax,
      qMin = inj.qMin,
      qMax = inj.qMax,
      vm_pu = vm,
      va_deg = va,
      referencePri = is_slack ? inj.bus : nothing,
      isRegulated = inj.isRegulated || is_slack,
      # only machine tuples carry a participation factor; every other pending
      # injection type reads out as `nothing` here
      participationFactor = get(inj, :participationFactor, nothing),
      # only machine and network-injection tuples carry a priority; an
      # EquivalentInjection or an SVC states none (0)
      referencePriority = get(inj, :referencePriority, 0),
      defer_bus_type_refresh = true,
    )
    # Record the index per source mRID for post-build steps (HVDC pairs).
    inj_mrid = get(inj, :mrid, nothing)
    inj_mrid === nothing || (ctx.prosumer_by_mrid[String(inj_mrid)] = lastindex(net.prosumpsVec))
    # every pending tuple carries its source mRID, an EquivalentInjection
    # as src_mrid (`something` throws for a tuple with neither); SVCs,
    # EnergySources and PowerElectronicsConnections are injection-type
    # prosumers like the machines, converters and equivalents share the
    # external-grid type string
    addr_mrid = String(something(inj_mrid, get(inj, :src_mrid, nothing)))
    addr_kind = inj.type == "EXTERNALNETWORKINJECTION" ? :external_grid : :generator
    ctx.element_by_mrid[addr_mrid] = (addr_kind, lastindex(net.prosumpsVec))
    # The prosumer constructor collapses the type string to a generic
    # Generator/Load component; keep the concrete class on the comp (same
    # pattern as the AsynchronousMachine marker) so a later CGMES export can
    # write each unit back as what it was.
    if inj.type == "SYNCHRONOUSMACHINE"
      net.prosumpsVec[end].comp.cTyp = Sparlectra.SynchronousMachine
      _captureBusEquipment!(net, ctx, "SM", inj.bus, get(inj, :mrid, nothing); children = (T1 = get(inj, :terminal_mrid, nothing), GU = get(inj, :gu_mrid, nothing), RC = get(inj, :rc_mrid, nothing)))
    elseif inj.type == "EXTERNALNETWORKINJECTION" && get(inj, :converter, false)
      # an EquivalentInjection the 2.4.15 export wrote for a converter (see
      # _PEC_EI_DESCRIPTION): the converter type again, and like a
      # PowerElectronicsConnection no captured identity, so the export
      # counter of external injections at the bus stays as the export saw it
      net.prosumpsVec[end].comp.cTyp = Sparlectra.PowerElectronicsConnection
    elseif inj.type == "EXTERNALNETWORKINJECTION"
      net.prosumpsVec[end].comp.cTyp = Sparlectra.ExternalNetworkInjection
      # EquivalentInjections arrive under the same type string without an
      # mRID: they consume the counter but keep a minted export identity
      _captureBusEquipment!(net, ctx, "ENI", inj.bus, get(inj, :mrid, nothing); children = (T1 = get(inj, :terminal_mrid, nothing), RC = get(inj, :rc_mrid, nothing)))
    elseif inj.type == "STATICVARCOMPENSATOR"
      net.prosumpsVec[end].comp.cTyp = Sparlectra.StaticVarCompensator
      _captureBusEquipment!(net, ctx, "SVC", inj.bus, get(inj, :mrid, nothing); children = (T1 = get(inj, :terminal_mrid, nothing), RC = get(inj, :rc_mrid, nothing)))
    elseif inj.type == "POWERELECTRONICSCONNECTION"
      # no captured identity: the 2.4.15 exporter writes it as an
      # EquivalentInjection under a minted mRID
      net.prosumpsVec[end].comp.cTyp = Sparlectra.PowerElectronicsConnection
    end
    # surviving remote-control plans pick up their prosumer index here (the
    # prosumer just pushed is the machine); the controller itself is attached
    # by _attachMachineControls! after bus types and isolation are final
    rt = get(inj, :remote_target, nothing)
    if rt !== nothing
      push!(ctx.machine_ctrl_plans, (bus = inj.bus, prosumer_idx = length(net.prosumpsVec), name = something(get(inj, :name, nothing), inj.bus), target_bus = rt.bus, target_vm_pu = rt.vset_pu, qmin = inj.qMin, qmax = inj.qMax))
    end
  end
  n_pf = count(inj -> get(inj, :participationFactor, nothing) !== nothing, pending)
  n_pf > 0 && push!(ctx.messages, "participation factors: $(n_pf) machines with normalPF")
  if !isempty(ctx.p_limit_placeholders)
    counts = join(["$(n) $(cls)" for (cls, n) in sort(collect(ctx.p_limit_placeholders))], ", ")
    push!(ctx.messages, "notice: P limits of $(round(Int, Sparlectra.SC_PMAX_PLACEHOLDER_MW)) MW and above read as no limit (placeholder): $(counts)")
  end
  return slack_bus
end

# --- short-circuit harvest (§7.7: read, not evaluated) ---------------------

function _busOrNothing(topo::CGMESTopology, obj::CIMObject)
  ts = get(topo.terminals, obj.mrid, nothing)
  ts === nothing && return nothing
  tn = tnOfTerminal(ts[1])
  tn === nothing && return nothing
  return get(topo.bus_name, tn, nothing)
end

function collectShortCircuitData(store::CGMESStore, topo::CGMESTopology)::CGMESShortCircuitData
  enis = NamedTuple[]
  for o in objectsOf(store, :ExternalNetworkInjection)
    push!(
      enis,
      (
        mrid = o.mrid,
        name = str(o, :name),
        bus = _busOrNothing(topo, o),
        maxInitialSymShCCurrent_A = num(o, :maxInitialSymShCCurrent),
        minInitialSymShCCurrent_A = num(o, :minInitialSymShCCurrent),
        maxR1ToX1Ratio = num(o, :maxR1ToX1Ratio),
        minR1ToX1Ratio = num(o, :minR1ToX1Ratio),
        maxR0ToX0Ratio = num(o, :maxR0ToX0Ratio),
        maxZ0ToZ1Ratio = num(o, :maxZ0ToZ1Ratio),
        ikSecond = boolval(o, :ikSecond),
        governorSCD = num(o, :governorSCD),
      ),
    )
  end
  machines = NamedTuple[]
  for o in objectsOf(store, :SynchronousMachine)
    # A machine on a wind or solar unit is a converter unit in fact: the
    # record names the class of the unit (`converter_unit`), and the engine
    # takes it as a current source. Its
    # rated current comes from ratedS, else from the maximum active power of
    # the unit, so `pmax_MW` carries GeneratingUnit.maxOperatingP for such a
    # machine only: a thermal machine keeps the size search of
    # sc_machine_size through the network unit, unchanged.
    gu = ref(store, o, :GeneratingUnit)
    converter_unit = (gu !== nothing && gu.class in (:WindGeneratingUnit, :SolarGeneratingUnit)) ? String(gu.class) : nothing
    pmax = converter_unit === nothing ? nothing : num(gu, :maxOperatingP)
    # WindGeneratingUnit.windGenUnitType (onshore, offshore) rides along so
    # the CGMES export writes the unit as delivered
    wind_type = gu === nothing || gu.class !== :WindGeneratingUnit ? nothing : str(gu, :windGenUnitType)
    wind_type === nothing || (wind_type = String(last(split(wind_type, '.'))))
    push!(
      machines,
      (
        mrid = o.mrid,
        name = str(o, :name),
        bus = _busOrNothing(topo, o),
        satDirectSubtransX_pu = num(o, :satDirectSubtransX),
        satDirectTransX_pu = num(o, :satDirectTransX),
        r0_pu = num(o, :r0),
        x0_pu = num(o, :x0),
        r2_pu = num(o, :r2),
        x2_pu = num(o, :x2),
        earthing = boolval(o, :earthing),
        ratedS_MVA = num(o, :ratedS),
        ratedU_kV = num(o, :ratedU),
        converter_unit = converter_unit,
        pmax_MW = pmax,
        windGenUnitType = wind_type,
      ),
    )
  end
  lines = NamedTuple[]
  for o in objectsOf(store, :ACLineSegment)
    if num(o, :r0) !== nothing || num(o, :x0) !== nothing || num(o, :shortCircuitEndTemperature) !== nothing
      push!(
        lines,
        (mrid = o.mrid, name = str(o, :name), r0_ohm = num(o, :r0), x0_ohm = num(o, :x0), b0ch_S = num(o, :b0ch), g0ch_S = num(o, :g0ch), shortCircuitEndTemperature_C = num(o, :shortCircuitEndTemperature)),
      )
    end
  end
  ptes = NamedTuple[]
  for o in objectsOf(store, :PowerTransformerEnd)
    if num(o, :r0) !== nothing || num(o, :x0) !== nothing || boolval(o, :grounded) !== nothing
      push!(
        ptes,
        (mrid = o.mrid, name = str(o, :name), r0_ohm = num(o, :r0), x0_ohm = num(o, :x0), grounded = boolval(o, :grounded), rground_ohm = num(o, :rground), xground_ohm = num(o, :xground)),
      )
    end
  end
  eis = NamedTuple[]
  for o in objectsOf(store, :EquivalentInjection)
    if num(o, :r) !== nothing || num(o, :x) !== nothing || num(o, :r0) !== nothing
      push!(
        eis,
        (mrid = o.mrid, name = str(o, :name), bus = _busOrNothing(topo, o), r_ohm = num(o, :r), x_ohm = num(o, :x), r0_ohm = num(o, :r0), x0_ohm = num(o, :x0), r2_ohm = num(o, :r2), x2_ohm = num(o, :x2)),
      )
    end
  end
  # Asynchronous machines: everything the IEC 60909-0 §6.7 motor
  # impedance needs: locked-rotor current ratio, R/X of the locked rotor,
  # and the rated data to form S_rM (directly, or via mechanical power with
  # efficiency and power factor).
  asms = NamedTuple[]
  for o in objectsOf(store, :AsynchronousMachine)
    push!(
      asms,
      (
        mrid = o.mrid,
        name = str(o, :name),
        bus = _busOrNothing(topo, o),
        iaIrRatio = num(o, :iaIrRatio),
        rxLockedRotorRatio = num(o, :rxLockedRotorRatio),
        efficiency_percent = num(o, :efficiency),
        ratedMechanicalPower_MW = num(o, :ratedMechanicalPower),
        polePairNumber = num(o, :polePairNumber),
        ratedS_MVA = num(o, :ratedS),
        ratedU_kV = num(o, :ratedU),
        ratedPowerFactor = num(o, :ratedPowerFactor),
      ),
    )
  end
  # Converter-connected units (CGMES 3.0): PowerElectronicsConnection.maxIFault
  # is the maximum fault current "in per-unit of rated current" (CIM17, type
  # PU), the rated current follows from ratedS and ratedU. No CGMES 3.0
  # profile (EQ, SC, SSH, DY, OP) defines maxIFault, so a conforming
  # delivery has none and the record carries nothing; it is read where a
  # delivery carries it beyond the profile. importCGMES keeps the records
  # of the connections it mapped.
  pecs = NamedTuple[]
  for o in objectsOf(store, :PowerElectronicsConnection)
    push!(pecs, (mrid = o.mrid, name = str(o, :name), bus = _busOrNothing(topo, o), maxIFault_pu = num(o, :maxIFault), ratedS_MVA = num(o, :ratedS), ratedU_kV = num(o, :ratedU)))
  end
  return CGMESShortCircuitData(enis, machines, lines, ptes, eis, asms, pecs)
end

# --- public API -------------------------------------------------------------

"""
    importCGMES(; path, baseMVA=100.0, require_boundary=true, name=nothing) -> CGMESImportResult

Import a CGMES delivery (EQ+SSH+TP, plus boundary set) into a Sparlectra
`Net` (Stage 1: bus-branch topology from the TP profile, fixed SSH tap
positions). Short-circuit source data is harvested into the result's
`shortcircuit` container (read, not evaluated). Start voltages are taken
from the SV profile where present. A three-winding star point has no SV
voltage; it starts at the voltage at which no current leaves it, given the
SV voltages of its sides.

`path` is a folder, a ZIP or a vector of either, or a `CGMESStore` already
loaded with `loadCGMES`. The store is not modified: a delivery without a TP
profile gets its topological nodes from the node-breaker topology
processor, which writes them into the store, so the import then works on a
`deepcopy`. Otherwise `result.store` is the given store itself.

`result.element_by_mrid` names the element every mapped CGMES object
became (see `CGMESImportResult`).

An `EquivalentBranch` with both terminals connected and in service becomes
a branch without charging: on one voltage level a pi branch with `r`, `x`
on that base, across two levels a branch with per-unit ratio 1 and the
impedance on the base of terminal 1. One whose terminals sit on the same
TopologicalNode builds nothing and gets a notice. One whose impedance lies
below `equivalent_branch_min_z_pu` (per unit on the base of terminal 1, both
directions, default `1e-6`) builds no branch either: its two buses are
joined as a bus link, with a notice. With `r21`/`x21` other
than `r`/`x` the object becomes a transformer branch with ratio 1, phase
shift `arg(z21/z12)/2`, series impedance of magnitude `sqrt(|z12| |z21|)`
and unequal shunt arms with zero row sums (phase stamp, on both level
kinds). It reproduces the two-port exactly for `|z21| = |z12|`; otherwise
each admittance entry is off in magnitude by `sqrt(|z21|/|z12|)` or its
inverse. One notice counts the mapped, skipped and stamped objects and the
largest magnitude difference.

With `tap_control = true` the SSH tap positions are the start point and the
CGMES-defined tap changers become outer-loop controllers. With
`machine_control = true` a machine whose voltage `RegulatingControl` points
at a *different* bus becomes a PQ injection with an outer-loop
`MachineVoltageControl` regulating that remote bus (instead of being held PV
at its own bus); plans whose target bus is already voltage-held, isolated,
or claimed by another machine fall back to the held-PV behavior with a
notice in `result.messages`.

A `SynchronousMachine`, `ExternalNetworkInjection` or `StaticVarCompensator`
with SSH `controlEnabled = false` does not regulate: it is a PQ injection at
its SSH `p`/`q` with its Q limits, whatever its `RegulatingControl` says,
and a machine gets no `MachineVoltageControl`; an absent flag counts as
`true`. The flag does not touch the slack choice: `referencePriority` and
the SV angle reference still apply, and as slack the unit holds its bus at
the start voltage (the SV value where present). One notice per class counts
these units. An `EquivalentInjection` has no `controlEnabled`; it regulates
only with EQ `regulationCapability = true`, SSH `regulationStatus = true` and
a positive `regulationTarget` (an absent flag counts as `false`). One notice
counts the units with `regulationCapability = true` and no
`regulationStatus`.

An active power limit of 9999 MW and above in magnitude
(`GeneratingUnit.minOperatingP`/`maxOperatingP`, `minP`/`maxP` of an
`ExternalNetworkInjection`, `EquivalentInjection` or `PowerElectronicsUnit`)
is a placeholder and arrives as no limit (`nothing`): the unit is not ranked
by it in a slack choice and gets no weight from it in the distributed slack.
One notice counts these units per class. Q limits keep their delivered
values.

With `map_energy_source = true` every in-service, connected `EnergySource`
on a TopologicalNode becomes a fixed PQ injection (`GENERATOR` prosumer,
`p = -activePower`, `q = -reactivePower` from the SSH, no regulation, no Q
limits, no slack candidate of the delivery), addressed as `:generator` in
`element_by_mrid`; one notice counts the mapped and skipped objects with
their P and Q sums. In an island whose only infeed is `EnergySource` the
importer makes the largest of them the slack, by the rule for islands
without a reference, with a warning (its P then follows the island
balance). With `multi_slack = false` such an island gets no slack: the
island-wise solve makes the strongest of them its reference, and one notice
per such island says so. With
the default `false` the objects stay out of the network, are listed in
`skipped_equipment`, and one notice names the key that maps them; their
count and P sum stand in the table of classes not mapped.

With `map_power_electronics = true` every in-service, connected
`PowerElectronicsConnection` (CGMES 3.0) on a TopologicalNode becomes an
injection prosumer of component type `PowerElectronicsConnection` with
`p = -p` and `q = -q` from the SSH, addressed as `:generator` in
`element_by_mrid`. With
`controlEnabled` and an enabled voltage `RegulatingControl` with a target
above zero it regulates like a `SynchronousMachine` (local target: PV at
its bus; a target outside `vset_min_pu`/`vset_max_pu`: held PV at its bus
voltage with a warning; remote target: held PV at its bus, or an
outer-loop controller with `machine_control = true`), otherwise it is a
fixed PQ injection. Q limits
come from the `minQ`/`maxQ` hull (a regulating one without them gets
±max(`ratedS`, `baseMVA`)), P limits from the summed `minP`/`maxP` of its
`PowerElectronicsUnit`s. It is no slack candidate of the delivery; an
island whose only sources are such PQ injections is treated like an
`EnergySource` island (an importer slack, or with `multi_slack = false` the
same notice). One notice counts the mapped, regulating and skipped
objects with their P and Q sums. Each mapped connection leaves a record in
`shortcircuit.power_electronics_connections` (`maxIFault`, `ratedS`,
`ratedU`), from which `runShortCircuit!` forms its fault current.
`maxIFault` is a CIM17 attribute outside the CGMES 3.0 profiles, so a
conforming delivery has none and its connections flag their island in the
short circuit instead. With the default `false` the connections
stay out of the network like unmapped `EnergySource` objects: listed in
`skipped_equipment`, one notice naming the key, count and P sum in the
table of classes not mapped. An `EquivalentInjection` with
`IdentifiedObject.description` = `PowerElectronicsConnection` (the CGMES
2.4.15 export writes a connection so) gets the component type
`PowerElectronicsConnection`, whatever the key: it has no fault current,
so the short circuit flags its island.

A `SynchronousMachine` on a `WindGeneratingUnit` or `SolarGeneratingUnit`
is mapped like any machine for the power flow; its short-circuit record
names the unit class (`converter_unit`) and carries
`GeneratingUnit.maxOperatingP`, so `runShortCircuit!` takes it as a
converter current source (`short_circuit.converter_current_factor` times
the rated current from `ratedS`, else from `maxOperatingP`).

With `shunt_control = true` every mapped `LinearShuntCompensator` with SSH
`controlEnabled` and an enabled voltage `RegulatingControl` at its own bus
gets a switched-bank `ShuntVoltageControl` that takes over its shunt: blocks
of `bPerSection` at the nominal voltage from 0 to `maximumSections`, the SSH
`sections` as start, the target in pu of the bus, half the `targetDeadband`
as deadband. A remote target, a target outside the `vset` band, a
conductance per section, a bus that is not PQ, a second bank on one bus and
a `NonlinearShuntCompensator` keep the shunt fixed with a notice.

With `read_dynamics = true` the DY profile is read as well and the machines
get their `DynamicData` (`_mapDynamics!`): inertia, damping and rating from
the synchronous machine models, droop and governor time constant from the
turbine governors whose parameter meaning the specification states. One
`notice:` per DY class counts the objects mapped, not mapped and without
counterpart in the EQ. Without the switch the DY file is skipped and no
result changes; a store passed in with DY objects is not read either.

The SV angle reference (`TopologicalIsland.AngleRefTopologicalNode`) also
takes an `EquivalentInjection` as slack, ranked behind a
`SynchronousMachine` or `ExternalNetworkInjection` on an SV angle
reference (MicroGrid BE: the X-node of the boundary). With
`equivalent_injection_as_source = true` an `EquivalentInjection` on the SV
angle reference node or on a boundary node (boundary set or
`boundaryPoint = true`) becomes a voltage-regulating source at the SV
voltage of its node with its SSH P; the one on the SV angle reference is
the slack and holds the SV angle as well, the others carry it as reference
angle for an SCF export. `r`/`x` with `x > 0` give its short-circuit
record; one notice counts the sources and one names those without
impedance. Off, the injection stays as mapped (PQ at the SSH p and q, PV
with `regulationStatus`, slack when it holds the SV angle reference).

Every class with terminals or with an SSH active power value of which no
object outside the boundary set is mapped gets one row `class not mapped:`
with its object count and the sum of its SSH P in MW, below one header
notice that states the load sign (positive means flow out of the node). A
row whose objects carry no SSH active power has no P sum. An object is
mapped when it became a network element (also as a three-winding leg or as
the charging shunt of a line on one node). The topology consumes
`BusbarSection` and `Junction`, and a switch (`Switch`, `Breaker`,
`Disconnector`, `LoadBreakSwitch`, `Jumper`, `ProtectedSwitch`, `Fuse`,
`GroundDisconnector`, `Cut`, `DisconnectingCircuitBreaker`) that is open,
disconnected or inside one TopologicalNode; a closed switch between two
nodes is a bus link for the first six classes and missing for the others.
Parts of a mapped object (transformer ends, tap changers, generating units)
and the connectivity carry no terminal of their own and get no row. An
object out of service or without a connected terminal on a TopologicalNode
does not count when its class is one the importer maps.
"""
function importCGMES(;
  path,
  baseMVA::Float64 = 100.0,
  # how bus shunts are represented; consumed while the shunts are built, so it
  # has to arrive here rather than being stamped on the finished network
  bus_shunt_model = :admittance,
  require_boundary::Bool = true,
  tap_control::Bool = false,
  machine_control::Bool = false,
  ignore_connected::Bool = false,
  multi_slack::Bool = true,
  equivalent_branch_min_z_pu::Float64 = 1e-6,
  vset_min_pu::Float64 = CGMES_VSET_MIN_PU,
  vset_max_pu::Float64 = CGMES_VSET_MAX_PU,
  strict_placeholder_guards::Bool = false,
  infer_base_voltages::Bool = false,
  hvdc_mode::Symbol = :injections,
  map_energy_source::Bool = false,
  map_power_electronics::Bool = false,
  shunt_control::Bool = false,
  read_dynamics::Bool = false,
  equivalent_injection_as_source::Bool = false,
  name::Union{Nothing,String} = nothing,
)::CGMESImportResult
  vset_min_pu <= vset_max_pu || throw(ArgumentError("importCGMES: vset_min_pu ($(vset_min_pu)) must not exceed vset_max_pu ($(vset_max_pu))"))
  hvdc_mode in (:injections, :paired_control) || throw(ArgumentError("importCGMES: hvdc_mode must be :injections or :paired_control, got :$(hvdc_mode)"))
  # the topology processor below is the one step that writes into the store;
  # a caller's store is copied only when it will run
  store = if path isa CGMESStore
    _nodeBreakerProcessorApplies(path) ? deepcopy(path) : path
  else
    # the DY profile only with cgmes_import.read_dynamics; without it the DY
    # file is skipped as before ("profile not in import set")
    loadCGMES(path; profile_filter = read_dynamics ? IMPORT_PROFILE_TAGS_WITH_DY : IMPORT_PROFILE_TAGS)
  end
  equivalent_branch_min_z_pu >= 0.0 || throw(ArgumentError("importCGMES: equivalent_branch_min_z_pu must not be negative, got $(equivalent_branch_min_z_pu)"))
  ctx = _MapCtx(tap_control, ignore_connected; machine_control = machine_control, vset_min_pu = vset_min_pu, vset_max_pu = vset_max_pu, strict_guards = strict_placeholder_guards, map_energy_source = map_energy_source, map_power_electronics = map_power_electronics, eb_min_z_pu = equivalent_branch_min_z_pu, ei_as_source = equivalent_injection_as_source)
  ignore_connected && push!(ctx.messages, "notice: ignore_connected = true, SSH Terminal.connected flags are overridden, everything is treated as in service")
  messages = ctx.messages
  unresolved = unresolvedReferences(store)
  if !isempty(unresolved)
    topo_refs = count(u -> u.key in (:TopologicalNode, :ConnectivityNode), unresolved)
    if topo_refs > 0 && require_boundary
      # Print the full why-analysis (supplied models, missing declared
      # dependencies, unresolved-reference histogram) for CLI users and carry
      # it in the typed error so service runs can write it into
      # cgmes_import_report.md:
      # a failed customer delivery must name its own gap, wherever it fails.
      analysis = importFailureAnalysis(store)
      print(analysis)
      hint = infer_base_voltages ? " infer_base_voltages is set, additionally set require_boundary=false to import without the boundary." : ""
      throw(CGMESImportError("CGMES import: $(topo_refs) unresolved topology references. See the import analysis (CLI: printed above; Web UI runs: cgmes_import_report.md). Pass the missing boundary/dependency as additional path, or set require_boundary=false.$(hint)", analysis))
    end
    push!(ctx.messages, "notice: $(length(unresolved)) unresolved references (non-fatal)")
  end
  for f in store.files
    f.skipped && push!(ctx.messages, "file skipped: $(f.name), $(f.skip_reason)")
  end
  # a DY profile the switch leaves out is named once with the switch, so the
  # report does not only list it as a skipped file
  dy_off = read_dynamics ? nothing : _dynamicsNotRead(store)
  dy_off === nothing || push!(ctx.messages, "notice: " * only(dynamicsReportLines(dy_off)))

  # Multi-valued references (#294 point 9): protect against silently reading
  # the wrong value: `ref()` returns only the FIRST occurrence of a repeated
  # property; the full list is available via `refsAll`. One notice per
  # (class, property) so a future mapping that touches such a property cannot
  # fall into the last-wins trap unwarned.
  if !isempty(store.multirefs)
    counts = Dict{Tuple{Symbol,Symbol},Int}()
    for ((mrid, key), _) in store.multirefs
      o = get(store.objects, mrid, nothing)
      o === nothing && continue
      counts[(o.class, key)] = get(counts, (o.class, key), 0) + 1
    end
    for ((cls, key), n) in sort!(collect(counts))
      push!(ctx.messages, "notice: multi-valued reference $(cls).$(key) on $(n) object(s): ref() reads the first value, the full list is in refsAll")
    end
  end

  # node-breaker deliveries without a TP profile (#314): derive the bus
  # partition from connectivity nodes and switch states BEFORE the
  # bus-branch topology is built; TP-carrying deliveries skip this
  applyNodeBreakerTopologyProcessor!(store, ctx)

  # with neither a TP profile nor node-breaker connectivity there is no
  # topology information at all; explain instead of mapping into a void
  # (boundary-set TNs alone do not count: TP_BD only describes the seam,
  # so a delivery whose own area has no TN would otherwise fall through
  # to a confusing downstream error, e.g. "no slack candidate")
  if !any(tn -> !(tn.mrid in store.boundary), objectsOf(store, :TopologicalNode))
    analysis = importFailureAnalysis(store)
    print(analysis)
    throw(CGMESImportError("CGMES import: the delivery carries no topology information, neither a TP profile (TopologicalNode) nor usable node-breaker data (ConnectivityNode plus switches). A bus-branch delivery needs its TP file; a node-breaker delivery needs its EQ connectivity. See the import analysis (CLI: printed above; Web UI runs: cgmes_import_report.md).", analysis))
  end

  topo = buildTopology(store; infer_base_voltages = infer_base_voltages, messages = ctx.messages)
  svmap = _svVoltageByTN(store)
  netname = something(name, "cgmes_import")
  net = Sparlectra.Net(name = netname, baseMVA = baseMVA, bus_shunt_model = Sparlectra.normalize_bus_shunt_model(bus_shunt_model))
  created = Set{String}()

  # must run before any element mapping: line endpoints and equivalent
  # injections at non-cancelling boundary nodes land on per-side buses
  _detectNonCancellingBoundarySides!(ctx, store, topo)

  acl_parallel = _mapLines!(net, store, topo, created, svmap, baseMVA, ctx)
  _mapEquivalentBranches!(net, store, topo, created, svmap, baseMVA, ctx, acl_parallel)
  # with tap control the SSH positions are the start point and the outer-loop
  # controllers find the solved positions themselves; without it the SvTapStep
  # positions reproduce the SV state directly (Stage 1)
  svsteps = tap_control ? Dict{String,Float64}() : _svTapPositions(store)
  _mapTransformers!(net, store, topo, created, svmap, svsteps, baseMVA, ctx)
  _mapSwitches!(net, store, topo, created, svmap, ctx)
  _mapShunts!(net, store, topo, created, svmap, ctx)
  _mapLoads!(net, store, topo, created, svmap, ctx)
  slack_bus = _mapInjections!(net, store, topo, created, svmap, ctx; multi_slack = multi_slack)
  # dynamic data onto the mapped machines; before the island slack choice,
  # which reads the unit type
  dynamics = read_dynamics ? _mapDynamics!(net, store, ctx) : dy_off
  # every element is mapped, so the stamps of the star rows are final
  _startStarBuses!(net, ctx, topo, svmap)
  Sparlectra.refreshBusTypesFromProsumers!(net)
  # Real deliveries carry TopologicalNodes whose only equipment is a load or a
  # switch that never becomes a branch (RealGrid: 198 such stubs). Mark them so
  # the solver skips them instead of reporting single-bus islands without a
  # reference; the bus stays in the model for provenance.
  _deEnergizeSourcelessComponents!(net, ctx)
  # The marking counts branches only. A bus that a closed switch joins to a
  # bus with a branch is marked as well, but the solver merges the switch
  # group onto one node and solves it with that bus (_merged_pf_net), so the
  # message counts such buses apart. They stay marked: N-1, the scenario
  # steps and the short circuit apply the same branch-only rule (QSTS mode A
  # marks its own steps again with links = true).
  n_iso = length(Sparlectra.markIsolatedBuses!(net = net))
  n_fed = length(Sparlectra._linkFedIsolatedBuses(net))
  n_iso > n_fed && push!(ctx.messages, "notice: $(n_iso - n_fed) isolated bus(es) without any in-service branch: marked isolated, excluded from the power flow")
  n_fed > 0 && push!(ctx.messages, "notice: $(n_fed) bus(es) without an in-service branch marked isolated, but a closed switch joins each to a bus with one: solved with that bus")
  # islands with a source but without reference get a slack; needs the final
  # isolation, and runs before the EnergySource notice, which then only names
  # islands left without reference (multi_slack = false)
  multi_slack && _slackForIslandsWithoutReference!(net, store, ctx)
  # needs the final isolation: an isolated bus belongs to no island
  _noticeEnergySourceOnlyIslands!(net, store, ctx)
  # every class with terminals or SSH active power of which no object is
  # mapped gets one row (_unmappedClasses), so missing generation, load or
  # network elements show up in the import report; the header states the
  # sign convention once for all rows, and a row without any SSH value has no
  # P sum
  unmapped = _unmappedClasses(store, ctx.element_by_mrid; map_energy_source = ctx.map_energy_source, map_power_electronics = ctx.map_power_electronics, ignore_connected = ctx.ignore_connected)
  isempty(unmapped) || push!(ctx.messages, "notice: classes not mapped, P sum from SSH in load sign (positive means flow out of the node):")
  for r in unmapped
    p_sum = r.p_mw === nothing ? "" : string(r.no_p == 0 ? "" : ", $(r.no_p) without SSH P", ", P sum $(round(r.p_mw; digits = 2)) MW")
    push!(ctx.messages, "notice: class not mapped: $(r.class), $(r.count) object(s)$(p_sum)")
  end
  tap_control && _disableIneffectiveTapControllers!(net, ctx)
  machine_control && _attachMachineControls!(net, ctx)
  shunt_control && _attachShuntControls!(net, store, topo, ctx)
  # HVDC pairing (#297 Draft B): detection runs on the raw store, so even in
  # the default Stage-0 mode the pairs are named in the messages; only
  # :paired_control attaches steerable controllers.
  if hvdc_mode === :paired_control
    _attachHvdcPairs!(net, ctx, store)
  else
    for pair in _detectHvdcPairs(store)
      pair.extra == 0 || continue
      # Stage-0 keeps the fixed injections but still records the link for
      # the result layer (r0.9.9)
      _push_cgmes_hvdc_link!(net, ctx, pair)
      kind = pair.segments ? "point-to-point" : "back-to-back"
      push!(ctx.messages, "notice: HVDC pair detected ($(kind)): $(pair.from) / $(pair.to), fixed injections (Stage 0); set cgmes_import.hvdc_mode = paired_control for a steerable pair")
    end
  end

  sc = collectShortCircuitData(store, topo)
  # a source the import did not map (terminal disconnected, out of service,
  # no node, a boundary side that cancels, a connection behind its switch)
  # is no element of the network and feeds no fault current. Its record
  # still carries the bus of its node, and that node can be live through
  # other equipment, so the record goes here and not in the engine.
  filter!(r -> haskey(ctx.prosumer_by_mrid, r.mrid), sc.power_electronics_connections)
  for field in (:external_network_injections, :synchronous_machines, :equivalent_injections, :asynchronous_machines)
    filter!(r -> haskey(ctx.element_by_mrid, r.mrid), getfield(sc, field))
  end
  isempty(sc.ac_line_segments) || push!(ctx.messages, "short-circuit data: $(length(sc.ac_line_segments)) lines (read, not evaluated)")
  isempty(sc.synchronous_machines) || push!(ctx.messages, "short-circuit data: $(length(sc.synchronous_machines)) machines (read, not evaluated)")
  isempty(sc.external_network_injections) || push!(ctx.messages, "short-circuit data: $(length(sc.external_network_injections)) external injections (read, not evaluated)")

  # SV coverage: buses whose TN carries no usable SvVoltage get the 1.0 pu /
  # 0° fallback from _ensureBus!, and then the start of a neighbour with one
  # where there is such a neighbour (_startNoSvBusesFromNeighbours!; runs
  # after the star start and the isolation marking, so star points and
  # isolated buses keep theirs). Their share decides how meaningful a
  # comparison against the delivery's SV state is (fixed-reference self-check,
  # compareWithSV), so record the list and say so in the messages. Per-side
  # split buses share their TN's SV classification and are covered through the
  # primary bus name; vn <= 0 counts as "no usable SV" because _ensureBus!
  # cannot form a per-unit magnitude there either, and an SvVoltage with
  # v <= 0 never reaches svmap (_svVoltageByTN).
  no_sv_buses = String[]
  for (tn, bus) in topo.bus_name
    bus in created || continue
    (haskey(svmap, tn) && topo.vn_kV[tn] > 0.0) || push!(no_sv_buses, bus)
  end
  sort!(unique!(no_sv_buses))
  isempty(no_sv_buses) || push!(ctx.messages, "notice: $(length(no_sv_buses)) of $(length(created)) buses have no usable SvVoltage (list in result.no_sv_buses)")
  n_nb = _startNoSvBusesFromNeighbours!(net, ctx, no_sv_buses)
  n_nb > 0 && push!(ctx.messages, "notice: $(n_nb) bus(es) without a usable SvVoltage start at the voltage of a neighbour with one")
  # A tap step the SV was not computed with keeps the import off its own SV
  # whatever the model does; one line here, the list per tap changer in
  # result.sv_tap_step_fits (_svTapStepFits, cgmes_report.jl)
  sv_tap_step_fits = _svTapStepFits(net, store, topo, ctx.element_by_mrid)
  isempty(sv_tap_step_fits) || push!(ctx.messages, "notice: $(length(sv_tap_step_fits)) tap changer(s) whose SvTapStep fits the SV flows of its transformer clearly worse than another step (list in result.sv_tap_step_fits)")

  # Everything else is collected silently and inspected via result.messages.
  # "warning:" is reserved for data defects where the importer had to substitute
  # a value to keep the model usable, those must not stay hidden in a vector
  # nobody reads, so they also go to the log.
  _logCGMESWarnings(ctx.messages)

  # the harvest travels with the network as well (the container is
  # field-identical), so an export or a short-circuit run of the network
  # alone has the source data of the delivery
  for field in fieldnames(CGMESShortCircuitData)
    append!(getfield(net.sc_sources, field), getfield(sc, field))
  end
  # every tap controller exists now, so its position in _tap_controllers(net)
  # is final (see _MapCtx); looked up by identity of the recorded controller
  if !isempty(ctx.tap_controller_by_mrid)
    tap_pos = IdDict{Sparlectra.PowerTransformerControl,Int}(c => k for (k, c) in enumerate(Sparlectra._tap_controllers(net)))
    for (mrid, c) in ctx.tap_controller_by_mrid
      ctx.element_by_mrid[mrid] = (:tap_controller, tap_pos[c])
    end
  end
  return CGMESImportResult(net, store, topo, sc, slack_bus, ctx.messages, ctx.branch_side, ctx.skipped, no_sv_buses, ctx.element_by_mrid, sv_tap_step_fits, dynamics)
end

"""
Emit the collected `warning:` messages through the logging system.

The mapping functions only collect strings; a data defect that made the
importer substitute a derived value has to reach the user's log as well, not
just `result.messages`. Grouped into a single `@warn` so a delivery with many
defective objects does not flood the log.
"""
function _logCGMESWarnings(messages::Vector{String})
  warnings = [m for m in messages if startswith(m, "warning:")]
  isempty(warnings) && return nothing
  @warn "CGMES import: $(length(warnings)) data defect(s): substituted values, see result.messages" details = join(warnings, "\n")
  return nothing
end

"""
    createNetFromCGMES(; path, baseMVA=100.0, require_boundary=true, name=nothing) -> Net

Thin wrapper around [`importCGMES`](@ref) returning only the `Net`.
"""
createNetFromCGMES(; kwargs...) = importCGMES(; kwargs...).net
