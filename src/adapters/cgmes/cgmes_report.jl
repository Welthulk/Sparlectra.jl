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

# file: src/adapters/cgmes/cgmes_report.jl
# purpose: diagnostic-first reporting for the CGMES importer (Stage 0).
# `summarizeCGMES` inspects any delivery, including broken ones, and
# reports profiles, version, class histogram and unresolved references
# before any Net mapping exists.

"""
Result of `summarizeCGMES`: everything a user needs to judge a CGMES
delivery before importing it.
"""
struct CGMESSummary
  files::Vector{CGMESFileInfo}
  version::String
  object_count::Int
  class_histogram::Vector{Pair{Symbol,Int}}   # sorted, descending
  unresolved_count::Int
  unresolved_sample::Vector{String}           # up to 5 example target mRIDs
  boundary_missing_hint::Bool
end

"""
    summarizeCGMES(; path, read_dynamics = false) -> CGMESSummary

Read a CGMES delivery (folder, ZIP, or vector of paths, e.g. base case +
boundary) and produce a diagnostic summary without building a `Net`. Works
on incomplete data sets; unresolved `TopologicalNode`/`ConnectivityNode`
references raise the `boundary_missing_hint`. `read_dynamics` reads the DY
profile as `importCGMES` does with the same keyword; without it a DY file is
listed as skipped.
"""
function summarizeCGMES(; path, read_dynamics::Bool = false)::CGMESSummary
  store = loadCGMES(path; profile_filter = read_dynamics ? IMPORT_PROFILE_TAGS_WITH_DY : IMPORT_PROFILE_TAGS)
  histogram = sort([cls => length(mrids) for (cls, mrids) in store.byclass]; by = p -> -p.second)
  unresolved = unresolvedReferences(store)
  # boundary hint: dangling refs whose key points at topology container objects
  boundary_keys = (:TopologicalNode, :ConnectivityNode)
  boundary_missing = any(u -> u.key in boundary_keys, unresolved)
  sample = [u.target for u in Iterators.take(unresolved, 5)]
  return CGMESSummary(store.files, store.version, length(store.objects), histogram, length(unresolved), sample, boundary_missing)
end

"""
    compareWithSV(result::CGMESImportResult) -> NamedTuple

Compare the solved state of `result.net` (after `runpf!`) with the SV
profile of the imported data set, the numeric acceptance oracle of the
importer.

Voltages: per-bus Δvm/Δva vs `SvVoltage` (max/RMS + rows). Isolated buses
are skipped, except a bus that a closed switch joins to a bus with a branch:
the solver merges it with that bus, and it is compared at the voltage of its
switch group. Angles are only
defined up to one constant per island, an IGM cut out of the continental
CGM keeps the CGM's global angle reference while the local solve pins its
own slack, which shows up as a uniform offset of tens of degrees that says
nothing about the state. The comparison therefore removes the median angle
offset PER ISLAND (the islands the power flow solves: branches and closed
bus links both connect, also through a bus isolated only for lack of a
branch of its own) and judges the aligned deltas (`dva_aligned`, feeding
`max_dva`/`rms_dva`); the raw `dva` stays in the rows. `va_ref_offset_deg`
reports the offset of the largest island, the one that describes the
delivery against the CGM reference. Aligning every island with the main
island's offset would measure the island reference instead of the state.

The alignment removes REFERENCE differences, never state differences. A
large `dva_aligned` on a bus inside the main island is therefore a real
disagreement with the SV profile and has to be explained from the model,
not from the reference: on the FullGrid test configuration the largest
angle deltas sit at the HVDC converter buses, whose SSH data declare
150 MW drawn at one end and 0 MW injected at the other while
`cgmes_import.hvdc_mode = injections` maps the two ends as independent
fixed injections without DC coupling.

Flows (`.flows`): per-terminal comparison vs `SvPowerFlow` in the CGMES
sign convention (power flowing *into* the equipment):
- `:branch_from`/`:branch_to`: line/transformer terminals, model flow
  computed from the solved voltages (real exchanges carry these; the
  ENTSO-E conformity sets only ship injection terminals),
- `:shunt`: `LinearShuntCompensator` and `NonlinearShuntCompensator` at
  the solved bus voltage (the nonlinear one with the admittance the import
  summed from its characteristic points up to the SSH section),
- `:load`: SSH load values (data-consistency check),
- `:units`: machines/ENIs/EquivalentInjections aggregated per closed-link
  cluster against the solved balance of the cluster (covers PV and slack
  units; a bus without links is a cluster of one).
Terminals of skipped equipment are excluded.
"""
function compareWithSV(result::CGMESImportResult)
  net = result.net
  store = result.store
  topo = result.topo
  raw_rows = NamedTuple{(:bus, :vm_pu, :sv_vm_pu, :dvm, :va_deg, :sv_va_deg, :dva),Tuple{String,Float64,Float64,Float64,Float64,Float64,Float64}}[]
  iso = Set(net.isoNodes)
  fed = Sparlectra._linkFedIsolatedBuses(net)
  row_bus = Int[]   # bus whose state a row reads (the switch group's for a link-fed bus)
  for sv in objectsOf(store, :SvVoltage)
    tn = get(sv.refs, :TopologicalNode, nothing)
    tn === nothing && continue
    haskey(topo.bus_name, tn) || continue
    bus = topo.bus_name[tn]
    haskey(net.busDict, bus) || continue
    # de-energized/isolated buses carry no solved state: comparing them
    # against the SV profile would only measure the de-energization. A bus
    # that a closed switch joins to a bus with a branch is solved with it and
    # compared at the voltage of its switch group (_linkFedIsolatedBuses).
    idx = net.busDict[bus]
    if idx in iso
      idx = get(fed, idx, 0)
      idx == 0 && continue
    end
    vn = topo.vn_kV[tn]
    v = num(sv, :v)
    (v === nothing || vn <= 0.0 || v <= 0.0) && continue
    node = net.nodeVec[idx]
    vm = something(node._vm_pu, NaN)
    va = something(node._va_deg, NaN)
    sv_vm = v / vn
    sv_va = num(sv, :angle, 0.0)
    push!(raw_rows, (bus = bus, vm_pu = vm, sv_vm_pu = sv_vm, dvm = vm - sv_vm, va_deg = va, sv_va_deg = sv_va, dva = va - sv_va))
    push!(row_bus, idx)
  end
  flows = _compareFlowsWithSV(result)
  aligned_rows = NamedTuple{(:bus, :vm_pu, :sv_vm_pu, :dvm, :va_deg, :sv_va_deg, :dva, :dva_aligned),Tuple{String,Float64,Float64,Float64,Float64,Float64,Float64,Float64}}[]
  isempty(raw_rows) && return (n = 0, max_dvm = NaN, rms_dvm = NaN, max_dva = NaN, rms_dva = NaN, va_ref_offset_deg = NaN, rows = aligned_rows, flows = flows)
  # One offset PER ISLAND, not one for the whole delivery. An angle is
  # defined up to a constant per island, so a single global offset fits the
  # main island and shifts every other island instead of aligning it.
  # Measured on the four-island FullGrid run: the
  # two-bus island around the promoted reference moved by 3.3 degrees when
  # the per-island offset replaced the global one. The large deltas on that
  # run are NOT explained by this, they sit inside the MAIN island (see the
  # docstring). The reported reference offset is the one of the LARGEST
  # island: it is the number that describes the delivery against the CGM
  # reference.
  #
  # The islands are the ELECTRICAL ones, closed links included. The
  # branch-only view of detect_ac_islands made a part that hangs on the rest
  # through a closed switch alone an island of its own, and its median then
  # absorbed its own deviation: on sp_casePST a 2 degree disagreement behind
  # the breaker was reported as 0. The islands are those the power flow
  # solves (_power_flow_island_of_bus): the electrical view still split at a
  # bus the branch-only marking isolates although closed switches join it to
  # both sides, and each side got its own offset. `Sparlectra.` is required:
  # this file is the submodule CGMESImporter and the island map lives in the
  # parent (an unqualified call compiles and throws UndefVarError at run
  # time).
  island_map = Sparlectra._power_flow_island_of_bus(net)
  row_island = Int[get(island_map, b, 0) for b in row_bus]
  alignment = Sparlectra._sv_align_angles(Float64[r.dva for r in raw_rows], row_island)
  va_ref_offset_deg = alignment.offset_deg
  for (r, dva_aligned) in zip(raw_rows, alignment.aligned)
    push!(aligned_rows, (bus = r.bus, vm_pu = r.vm_pu, sv_vm_pu = r.sv_vm_pu, dvm = r.dvm, va_deg = r.va_deg, sv_va_deg = r.sv_va_deg, dva = r.dva, dva_aligned = dva_aligned))
  end
  dvm = [abs(r.dvm) for r in aligned_rows]
  dva = [abs(r.dva_aligned) for r in aligned_rows]
  return (n = length(aligned_rows), max_dvm = maximum(dvm), rms_dvm = sqrt(sum(abs2, dvm) / length(dvm)), max_dva = maximum(dva), rms_dva = sqrt(sum(abs2, dva) / length(dva)), va_ref_offset_deg = va_ref_offset_deg, rows = aligned_rows, flows = flows)
end

# Bus voltages and bus injections S = V conj(Y V) in MVA (load sign) of the
# current state, as compareWithSV and the SV export read them. createYBUS
# leaves isolated buses out and returns the matrix in active-bus space, so it
# is re-embedded at full size with zero rows for them. A link-fed isolated
# bus (`_linkFedIsolatedBuses`, returned as `fed`) takes the voltage of its
# switch group and the power of its own shunts at that voltage, the only
# admittance it has without a branch. Returns `(V, S, fed)`.
function _svBusState(net)
  n_bus = length(net.nodeVec)
  fed = Sparlectra._linkFedIsolatedBuses(net)
  V = ComplexF64[something(nd._vm_pu, 1.0) * cis(deg2rad(something(nd._va_deg, 0.0))) for nd in net.nodeVec]
  for (b, m) in fed
    V[b] = V[m]
  end
  Yred = Sparlectra.createYBUS(net = net, sparse = true)
  Y = size(Yred, 1) == n_bus ? Yred : Sparlectra._expand_ybus_for_isolated_nodes(Yred, n_bus, net.isoNodes)
  S = V .* conj.(Y * V) .* net.baseMVA
  for sh in net.shuntVec
    # the shunts createYBUS stamps (in service, admittance model)
    (sh.status == 0 || sh.model != :Y) && continue
    haskey(fed, sh.busIdx) || continue
    S[sh.busIdx] += abs2(V[sh.busIdx]) * conj(sh.y_pu_shunt) * net.baseMVA
  end
  return V, S, fed
end

# an EnergySource (cgmes_import.map_energy_source) and a CGMES 3.0
# PowerElectronicsConnection (cgmes_import.map_power_electronics) are units
# only when mapped; unmapped ones are in skipped_equipment and never reach
# the class test
const _UNIT_CLASSES = (:SynchronousMachine, :ExternalNetworkInjection, :EquivalentInjection, :EnergySource, :PowerElectronicsConnection)
const _LOAD_CLASSES = (:EnergyConsumer, :ConformLoad, :NonConformLoad, :StationSupply)

function _compareFlowsWithSV(result::CGMESImportResult)
  net = result.net
  store = result.store
  topo = result.topo
  empty_rows = NamedTuple{(:kind, :name, :bus, :sv_p, :sv_q, :p, :q, :dp, :dq),Tuple{Symbol,String,String,Float64,Float64,Float64,Float64,Float64,Float64}}[]

  # solved complex bus voltages and the resulting bus injections. Every real
  # CGMES 3.0 delivery has isolated buses (out-of-service equipment); their
  # rows are filtered below: they carry no solved state, so comparing them
  # would only measure the de-energization. A bus that a closed switch joins
  # to a bus with a branch is solved and stays in (`fed`, see _svBusState).
  V, S, fed = _svBusState(net)
  iso = Set(b for b in net.isoNodes if !haskey(fed, b))

  # scheduled load per bus (consumption positive)
  load_sum = zeros(ComplexF64, length(net.nodeVec))
  for ps in net.prosumpsVec
    Sparlectra.isGenerator(ps.proSumptionType) && continue
    i = ps.comp.cFrom_bus
    load_sum[i] += complex(something(ps.pVal, 0.0), something(ps.qVal, 0.0))
  end

  branch_of = Dict(br.branchIdx => br for br in net.branchVec)
  rows = empty_rows
  # Units are compared per closed-link cluster, not per bus. A closed link
  # carries power that the bus balance S = V conj(Y V) does not see, so a
  # machine behind a retained breaker would be compared with the shunt and
  # load of its own node only. Over the whole cluster the link flows cancel.
  # A bus without links is a cluster of one, its row is the bus balance.
  reps = Sparlectra._active_link_representative_map(net)
  cluster_members = Dict{Int,Vector{Int}}()
  for b in eachindex(reps)
    push!(get!(Vector{Int}, cluster_members, reps[b]), b)
  end
  unit_sv = Dict{Int,ComplexF64}()      # cluster representative → Σ SvPowerFlow of unit terminals
  unit_names = Dict{Int,Vector{String}}()

  for pf in objectsOf(store, :SvPowerFlow)
    tmrid = get(pf.refs, :Terminal, nothing)
    tmrid === nothing && continue
    t = get(store.objects, tmrid, nothing)
    t === nothing && continue
    eqm = get(t.refs, :ConductingEquipment, nothing)
    eqm === nothing && continue
    eqm in result.skipped_equipment && continue
    eq = get(store.objects, eqm, nothing)
    eq === nothing && continue
    p_sv = num(pf, :p)
    q_sv = num(pf, :q)
    (p_sv === nothing || q_sv === nothing) && continue
    tn = tnOfTerminal(t)
    bus = tn === nothing ? nothing : get(topo.bus_name, tn, nothing)
    busidx = bus === nothing ? nothing : get(net.busDict, bus, nothing)
    # de-energized buses carry no solved state, their rows would only
    # measure the de-energization, not the power flow
    busidx !== nothing && busidx in iso && continue
    name = something(str(eq, :name), eqm)

    if haskey(result.branch_side_of_terminal, tmrid)
      idx, side = result.branch_side_of_terminal[tmrid]
      br = get(branch_of, idx, nothing)
      (br === nothing || br.status == 0) && continue
      Smodel = Sparlectra._sv_branch_terminal_flow(br, V, side, net.baseMVA)
      kind = side == :from ? :branch_from : :branch_to
      push!(rows, (kind = kind, name = name, bus = something(bus, "?"), sv_p = p_sv, sv_q = q_sv, p = real(Smodel), q = imag(Smodel), dp = real(Smodel) - p_sv, dq = imag(Smodel) - q_sv))
    elseif eq.class in _LOAD_CLASSES
      p = num(eq, :p, 0.0)
      q = num(eq, :q, 0.0)
      push!(rows, (kind = :load, name = name, bus = something(bus, "?"), sv_p = p_sv, sv_q = q_sv, p = p, q = q, dp = p - p_sv, dq = q - q_sv))
    elseif eq.class == :LinearShuntCompensator
      busidx === nothing && continue
      sections = something(num(eq, :sections), num(eq, :normalSections, 0.0))
      vn = topo.vn_kV[tn]
      vm2 = abs2(V[busidx])
      p = num(eq, :gPerSection, 0.0) * sections * vn^2 * vm2
      q = -num(eq, :bPerSection, 0.0) * sections * vn^2 * vm2
      push!(rows, (kind = :shunt, name = name, bus = something(bus, "?"), sv_p = p_sv, sv_q = q_sv, p = p, q = q, dp = p - p_sv, dq = q - q_sv))
    elseif eq.class == :NonlinearShuntCompensator
      busidx === nothing && continue
      # the admittance is the sum of the characteristic points up to the SSH
      # section; the mapped shunt carries it, so the row reads it there
      # instead of evaluating the points a second time. A compensator the
      # import left out (no points, implausible value) has no row.
      addr = get(result.element_by_mrid, eqm, nothing)
      (addr === nothing || first(addr) != :shunt) && continue
      sh = net.shuntVec[last(addr)]
      Smodel = abs2(V[busidx]) * conj(sh.y_pu_shunt) * net.baseMVA
      push!(rows, (kind = :shunt, name = name, bus = something(bus, "?"), sv_p = p_sv, sv_q = q_sv, p = real(Smodel), q = imag(Smodel), dp = real(Smodel) - p_sv, dq = imag(Smodel) - q_sv))
    elseif eq.class in _UNIT_CLASSES
      busidx === nothing && continue
      rep = reps[busidx]
      unit_sv[rep] = get(unit_sv, rep, zero(ComplexF64)) + complex(p_sv, q_sv)
      push!(get!(Vector{String}, unit_names, rep), name)
    end
  end

  for (rep, sv) in unit_sv
    # model: flow into the units = -(net injection + load) of the cluster
    Smodel = -sum(S[b] + load_sum[b] for b in cluster_members[rep])
    bus = net.nodeVec[rep].comp.cName
    name = join(sort(unit_names[rep]), "+")
    push!(rows, (kind = :units, name = name, bus = bus, sv_p = real(sv), sv_q = imag(sv), p = real(Smodel), q = imag(Smodel), dp = real(Smodel) - real(sv), dq = imag(Smodel) - imag(sv)))
  end

  isempty(rows) && return (n = 0, max_dp = NaN, max_dq = NaN, rms_dp = NaN, rms_dq = NaN, rows = rows)
  dp = [abs(r.dp) for r in rows]
  dq = [abs(r.dq) for r in rows]
  return (n = length(rows), max_dp = maximum(dp), max_dq = maximum(dq), rms_dp = sqrt(sum(abs2, dp) / length(dp)), rms_dq = sqrt(sum(abs2, dq) / length(dq)), rows = rows)
end

# --- tap steps against the SV flows ------------------------------------------
#
# A tapped two-winding transformer whose SvPowerFlow at both terminals is met
# clearly better by another step than by its SvTapStep: the SV was computed
# with a step other than the one it states, so no import that takes the
# SvTapStep can reproduce the SV there. The importer counts such tap
# changers in one notice and lists them in result.sv_tap_step_fits.
#
# Criterion, per tap changer with an SvTapStep, the other changers of the
# transformer held at theirs: the mismatch of step k is |S(k) - S_sv| at
# terminal 1 plus the same at terminal 2, in MVA, where S(k) is the model
# flow into the transformer at the SV voltages of its two buses with the
# transformer built at step k. The best step has the smallest mismatch over
# lowStep..highStep. It is reported when it is not the SvTapStep, its
# mismatch is at most a tenth of the SvTapStep mismatch and the gain is at
# least 1 MVA. The absolute floor keeps an SV the model reproduces (mismatch
# near zero at the SvTapStep) from flagging a neighbour step on noise.
const _SV_TAP_FIT_RATIO = 0.1
const _SV_TAP_FIT_MIN_GAIN_MVA = 1.0
# the row type _SvTapFitRow is defined in cgmes_mapping.jl, ahead of
# CGMESImportResult, which carries the rows

# Mismatch in MVA of the transformer built at `steps`. The transformer is
# mapped again by _map2WTrafo! into `scratch`, a network of its two buses,
# so the check follows the importer's transformer model (tap side, tap
# tables) instead of a copy of it; each call appends one branch.
function _svTapStepMismatch(scratch::Sparlectra.Net, created::Set{String}, ctx::_MapCtx, store::CGMESStore, topo::CGMESTopology, svmap::Dict{String,Tuple{Float64,Float64}}, pt::CIMObject, ends::Vector{CIMObject}, terminals::Vector{String}, svflow::Dict{String,ComplexF64}, steps::Dict{String,Float64})::Float64
  nbr = length(scratch.branchVec)
  _map2WTrafo!(scratch, store, topo, created, svmap, steps, scratch.baseMVA, pt, ends, ctx)
  length(scratch.branchVec) == nbr && return NaN
  br = scratch.branchVec[end]
  V = ComplexF64[something(nd._vm_pu, 1.0) * cis(deg2rad(something(nd._va_deg, 0.0))) for nd in scratch.nodeVec]
  m = 0.0
  for t in terminals
    side = get(ctx.branch_side, t, nothing)
    (side === nothing || side[1] != br.branchIdx) && return NaN
    m += abs(Sparlectra._sv_branch_terminal_flow(br, V, side[2], scratch.baseMVA) - svflow[t])
  end
  return m
end

# The tap changers of the in-service two-winding transformers of an import
# whose SvTapStep fits the SV flows clearly worse than another step (see the
# criterion above). Only transformers with an SvPowerFlow at both terminals
# and an SvVoltage (v > 0) at both buses are judged; three-winding
# transformers are not, their star point has no SV voltage.
function _svTapStepFits(net::Sparlectra.Net, store::CGMESStore, topo::CGMESTopology, element_by_mrid::Dict{String,Tuple{Symbol,Int}})::Vector{_SvTapFitRow}
  rows = _SvTapFitRow[]
  svsteps = _svTapPositions(store)
  isempty(svsteps) && return rows
  svmap = _svVoltageByTN(store)
  svflow = Dict{String,ComplexF64}()
  for pf in objectsOf(store, :SvPowerFlow)
    t = get(pf.refs, :Terminal, nothing)
    p = num(pf, :p)
    q = num(pf, :q)
    (t === nothing || p === nothing || q === nothing) && continue
    svflow[t] = complex(p, q)
  end
  for pt in objectsOf(store, :PowerTransformer)
    entry = get(element_by_mrid, pt.mrid, nothing)
    (entry === nothing || entry[1] != :transformer) && continue
    net.branchVec[entry[2]].status == 1 || continue
    ends = _endsOfTransformer(store, pt)
    length(ends) == 2 || continue
    terminals = String[]
    for e in ends
      t = get(e.refs, :Terminal, nothing)
      (t === nothing || !haskey(svflow, t)) && break
      push!(terminals, t)
    end
    length(terminals) == 2 || continue
    infos = [_terminalInfo(store, topo, e) for e in ends]
    any(i -> i === nothing || i.vn_kV <= 0.0 || !haskey(svmap, i.tn), infos) && continue
    tcs = [tc for e in ends for tc in _tapChangersOfEnd(store, e.mrid) if haskey(svsteps, tc.mrid)]
    isempty(tcs) && continue
    # every terminal counts as connected in the scratch build: the service
    # state is the imported branch's, checked above
    scratch = Sparlectra.Net(name = "sv_tap_fit", baseMVA = net.baseMVA)
    created = Set{String}()
    ctx = _MapCtx(false, true)
    steps = copy(svsteps)
    m_sv = _svTapStepMismatch(scratch, created, ctx, store, topo, svmap, pt, ends, terminals, svflow, steps)
    (isfinite(m_sv) && m_sv >= _SV_TAP_FIT_MIN_GAIN_MVA) || continue
    for tc in tcs
      low = num(tc, :lowStep)
      high = num(tc, :highStep)
      (low === nothing || high === nothing) && continue
      sv_pos = svsteps[tc.mrid]
      best_step, best = Int(round(sv_pos)), m_sv
      for k in Int(round(low)):Int(round(high))
        k == sv_pos && continue
        steps[tc.mrid] = Float64(k)
        m = _svTapStepMismatch(scratch, created, ctx, store, topo, svmap, pt, ends, terminals, svflow, steps)
        isfinite(m) && m < best && ((best_step, best) = (k, m))
      end
      steps[tc.mrid] = sv_pos
      (best_step != sv_pos && best <= _SV_TAP_FIT_RATIO * m_sv && m_sv - best >= _SV_TAP_FIT_MIN_GAIN_MVA) || continue
      push!(rows, (transformer = something(str(pt, :name), pt.mrid), tap_changer = something(str(tc, :name), tc.mrid), sv_step = Int(round(sv_pos)), best_step = best_step, sv_mismatch_mva = m_sv, best_mismatch_mva = best))
    end
  end
  return rows
end

"""
    shortCircuitCoverage(sc::CGMESShortCircuitData) -> Vector{NamedTuple}

Per-class completeness of the harvested short-circuit source data: one row
per element class with the record count and, per attribute, how many records
actually carry a value (`attribute => filled/total`). Identification fields
(`mrid`, `name`, `bus`) are excluded: coverage describes the electrical
attributes a future IEC 60909 evaluation (issue #277) would consume.
"""
function shortCircuitCoverage(sc::CGMESShortCircuitData)::Vector{NamedTuple}
  rows = NamedTuple[]
  for (label, records) in (
    ("ExternalNetworkInjection", sc.external_network_injections),
    ("SynchronousMachine", sc.synchronous_machines),
    ("ACLineSegment", sc.ac_line_segments),
    ("PowerTransformerEnd", sc.transformer_ends),
    ("EquivalentInjection", sc.equivalent_injections),
    ("AsynchronousMachine", sc.asynchronous_machines),
    ("PowerElectronicsConnection", sc.power_electronics_connections),
  )
    n = length(records)
    attrs = NamedTuple[]
    if n > 0
      for key in keys(records[1])
        key in (:mrid, :name, :bus) && continue
        filled = count(r -> r[key] !== nothing, records)
        push!(attrs, (attribute = key, filled = filled, total = n))
      end
    end
    push!(rows, (class = label, records = n, attributes = attrs))
  end
  return rows
end

"""
    printShortCircuitCoverage(io, sc)

Readable rendering of [`shortCircuitCoverage`](@ref): per class the record
count and each attribute's fill rate, `✓` when complete. The short-circuit
service prints it into `run.log`; the CGMES import report of a service run
(`cgmes_import_report.md`) shows the same rows as a table.
"""
function printShortCircuitCoverage(io::IO, sc::CGMESShortCircuitData)
  for row in shortCircuitCoverage(sc)
    if row.records == 0
      # NOTE: the harvest stores line/transformer/EI records only when at
      # least one short-circuit attribute is present, so for those classes
      # "none" means "no records WITH short-circuit data", not necessarily
      # "class absent from the delivery" (machines/ENIs are always recorded).
      println(io, "  ", rpad(row.class, 28), "none (no records with short-circuit attributes)")
      continue
    end
    println(io, "  ", rpad(row.class, 28), row.records, " record(s)")
    for a in row.attributes
      marker = a.filled == a.total ? "✓" : (a.filled == 0 ? "-" : string(a.filled, "/", a.total))
      println(io, "    ", rpad(String(a.attribute), 34), marker)
    end
  end
  return nothing
end

# plain `print(s)`/`string(s)` must render the same multi-line report as the
# REPL display: without this delegation they fall back to the default struct
# show, which dumps every field into one unreadable line (seen in notebooks)
Base.show(io::IO, s::CGMESSummary) = show(io, MIME"text/plain"(), s)

function Base.show(io::IO, ::MIME"text/plain", s::CGMESSummary)
  println(io, "CGMES summary")
  println(io, "  version: ", isempty(s.version) ? "unknown/mixed" : s.version)
  println(io, "  files:   ", length(s.files), " (", count(f -> !f.skipped, s.files), " read, ", count(f -> f.skipped, s.files), " skipped)")
  for f in s.files
    tag = f.skipped ? "skip: $(f.skip_reason)" : "read"
    println(io, "    - ", f.name, "  [", join(sort(collect(f.profiles)), ","), "]  ", tag)
  end
  println(io, "  objects: ", s.object_count, " in ", length(s.class_histogram), " classes")
  shown = min(15, length(s.class_histogram))
  for (cls, n) in Iterators.take(s.class_histogram, shown)
    println(io, "    ", rpad(string(cls), 34), n)
  end
  shown < length(s.class_histogram) && println(io, "    … ", length(s.class_histogram) - shown, " more classes")
  print(io, "  unresolved references: ", s.unresolved_count)
  if s.unresolved_count > 0
    print(io, "  (sample: ", join(s.unresolved_sample, ", "), ")")
    s.boundary_missing_hint && print(io, "\n  ⚠ topology references unresolved, boundary set likely missing")
  end
  println(io)
end
