# Copyright 2023–2026 Udo Schmitz
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

# file: src/cgmes/cgmes_report.jl
# purpose: diagnostic-first reporting for the CGMES importer (Stage 0).
# `summarizeCGMES` inspects any delivery — including broken ones — and
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
    summarizeCGMES(; path) -> CGMESSummary

Read a CGMES delivery (folder, ZIP, or vector of paths — e.g. base case +
boundary) and produce a diagnostic summary without building a `Net`. Works
on incomplete data sets; unresolved `TopologicalNode`/`ConnectivityNode`
references raise the `boundary_missing_hint`.
"""
function summarizeCGMES(; path)::CGMESSummary
  store = loadCGMES(path)
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
profile of the imported data set — the numeric acceptance oracle of the
importer (concept §8).

Voltages: per-bus Δvm/Δva vs `SvVoltage` (max/RMS + rows). Angles are only
defined up to one constant per island — an IGM cut out of the continental
CGM keeps the CGM's global angle reference while the local solve pins its
own slack, which shows up as a uniform offset of tens of degrees that says
nothing about the state. The comparison therefore removes the median angle
offset (`va_ref_offset_deg`, reported separately) and judges the aligned
deltas (`dva_aligned`, feeding `max_dva`/`rms_dva`); the raw `dva` stays in
the rows. Secondary islands with their own reference may keep a residual
offset.

Flows (`.flows`): per-terminal comparison vs `SvPowerFlow` in the CGMES
sign convention (power flowing *into* the equipment):
- `:branch_from`/`:branch_to` — line/transformer terminals, model flow
  computed from the solved voltages (real exchanges carry these; the
  ENTSO-E conformity sets only ship injection terminals),
- `:shunt` — `LinearShuntCompensator` at the solved bus voltage,
- `:load` — SSH load values (data-consistency check),
- `:units` — machines/ENIs/EquivalentInjections aggregated per bus against
  the solved bus balance (covers PV and slack units).
Terminals of skipped equipment are excluded.
"""
function compareWithSV(result::CGMESImportResult)
  net = result.net
  store = result.store
  topo = result.topo
  raw_rows = NamedTuple{(:bus, :vm_pu, :sv_vm_pu, :dvm, :va_deg, :sv_va_deg, :dva),Tuple{String,Float64,Float64,Float64,Float64,Float64,Float64}}[]
  for sv in objectsOf(store, :SvVoltage)
    tn = get(sv.refs, :TopologicalNode, nothing)
    tn === nothing && continue
    haskey(topo.bus_name, tn) || continue
    bus = topo.bus_name[tn]
    haskey(net.busDict, bus) || continue
    # de-energized/isolated buses carry no solved state — comparing them
    # against the SV profile would only measure the de-energization
    net.busDict[bus] in net.isoNodes && continue
    vn = topo.vn_kV[tn]
    v = num(sv, :v)
    (v === nothing || vn <= 0.0 || v <= 0.0) && continue
    node = net.nodeVec[net.busDict[bus]]
    vm = something(node._vm_pu, NaN)
    va = something(node._va_deg, NaN)
    sv_vm = v / vn
    sv_va = num(sv, :angle, 0.0)
    push!(raw_rows, (bus = bus, vm_pu = vm, sv_vm_pu = sv_vm, dvm = vm - sv_vm, va_deg = va, sv_va_deg = sv_va, dva = va - sv_va))
  end
  flows = _compareFlowsWithSV(result)
  aligned_rows = NamedTuple{(:bus, :vm_pu, :sv_vm_pu, :dvm, :va_deg, :sv_va_deg, :dva, :dva_aligned),Tuple{String,Float64,Float64,Float64,Float64,Float64,Float64,Float64}}[]
  isempty(raw_rows) && return (n = 0, max_dvm = NaN, rms_dvm = NaN, max_dva = NaN, rms_dva = NaN, va_ref_offset_deg = NaN, rows = aligned_rows, flows = flows)
  # Median instead of mean: robust against the few buses of secondary
  # islands whose own reference differs from the main island's.
  sorted_dva = sort([r.dva for r in raw_rows if isfinite(r.dva)])
  m = length(sorted_dva)
  va_ref_offset_deg = m == 0 ? 0.0 : (isodd(m) ? sorted_dva[(m + 1) ÷ 2] : (sorted_dva[m ÷ 2] + sorted_dva[m ÷ 2 + 1]) / 2.0)
  for r in raw_rows
    push!(aligned_rows, (bus = r.bus, vm_pu = r.vm_pu, sv_vm_pu = r.sv_vm_pu, dvm = r.dvm, va_deg = r.va_deg, sv_va_deg = r.sv_va_deg, dva = r.dva, dva_aligned = r.dva - va_ref_offset_deg))
  end
  dvm = [abs(r.dvm) for r in aligned_rows]
  dva = [abs(r.dva_aligned) for r in aligned_rows]
  return (n = length(aligned_rows), max_dvm = maximum(dvm), rms_dvm = sqrt(sum(abs2, dvm) / length(dvm)), max_dva = maximum(dva), rms_dva = sqrt(sum(abs2, dva) / length(dva)), va_ref_offset_deg = va_ref_offset_deg, rows = aligned_rows, flows = flows)
end

const _UNIT_CLASSES = (:SynchronousMachine, :ExternalNetworkInjection, :EquivalentInjection)
const _LOAD_CLASSES = (:EnergyConsumer, :ConformLoad, :NonConformLoad, :StationSupply)

function _compareFlowsWithSV(result::CGMESImportResult)
  net = result.net
  store = result.store
  topo = result.topo
  empty_rows = NamedTuple{(:kind, :name, :bus, :sv_p, :sv_q, :p, :q, :dp, :dq),Tuple{Symbol,String,String,Float64,Float64,Float64,Float64,Float64,Float64}}[]

  # solved complex bus voltages and the resulting bus injections. createYBUS
  # returns the matrix in active-bus space when isolated buses exist — every
  # real CGMES 3.0 delivery has some (out-of-service equipment), so re-embed it
  # to full size instead of giving up on the whole comparison. Isolated buses
  # get zero rows, and their per-bus rows are filtered below: they carry no
  # solved state, so comparing them would only measure the de-energization.
  n_bus = length(net.nodeVec)
  iso = Set(net.isoNodes)
  V = ComplexF64[something(nd._vm_pu, 1.0) * cis(deg2rad(something(nd._va_deg, 0.0))) for nd in net.nodeVec]
  Yred = Sparlectra.createYBUS(net = net, sparse = true)
  Y = size(Yred, 1) == n_bus ? Yred : Sparlectra._expand_ybus_for_isolated_nodes(Yred, n_bus, net.isoNodes)
  S = V .* conj.(Y * V) .* net.baseMVA

  # scheduled load per bus (consumption positive)
  load_sum = zeros(ComplexF64, length(net.nodeVec))
  for ps in net.prosumpsVec
    Sparlectra.isGenerator(ps.proSumptionType) && continue
    i = ps.comp.cFrom_bus
    load_sum[i] += complex(something(ps.pVal, 0.0), something(ps.qVal, 0.0))
  end

  branch_of = Dict(br.branchIdx => br for br in net.branchVec)
  rows = empty_rows
  unit_sv = Dict{Int,ComplexF64}()      # bus index → Σ SvPowerFlow of unit terminals
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
    # de-energized buses carry no solved state — their rows would only
    # measure the de-energization, not the power flow
    busidx !== nothing && busidx in iso && continue
    name = something(str(eq, :name), eqm)

    if haskey(result.branch_side_of_terminal, tmrid)
      idx, side = result.branch_side_of_terminal[tmrid]
      br = get(branch_of, idx, nothing)
      (br === nothing || br.status == 0) && continue
      ys = inv(br.r_pu + im * br.x_pu)
      ysh2 = (br.g_pu + im * br.b_pu) / 2
      tr = br.ratio == 0.0 ? 1.0 + 0im : br.tap_ratio * cis(deg2rad(br.phase_shift_deg))
      Vf, Vt = V[br.fromBus], V[br.toBus]
      Smodel = side == :from ? Vf * conj(((ys + ysh2) / abs2(tr)) * Vf - (ys / conj(tr)) * Vt) * net.baseMVA : Vt * conj((ys + ysh2) * Vt - (ys / tr) * Vf) * net.baseMVA
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
    elseif eq.class in _UNIT_CLASSES
      busidx === nothing && continue
      unit_sv[busidx] = get(unit_sv, busidx, zero(ComplexF64)) + complex(p_sv, q_sv)
      push!(get!(Vector{String}, unit_names, busidx), name)
    end
  end

  for (busidx, sv) in unit_sv
    # model: flow into the units = -(net bus injection + local load)
    Smodel = -(S[busidx] + load_sum[busidx])
    bus = net.nodeVec[busidx].comp.cName
    name = join(sort(unit_names[busidx]), "+")
    push!(rows, (kind = :units, name = name, bus = bus, sv_p = real(sv), sv_q = imag(sv), p = real(Smodel), q = imag(Smodel), dp = real(Smodel) - real(sv), dq = imag(Smodel) - imag(sv)))
  end

  isempty(rows) && return (n = 0, max_dp = NaN, max_dq = NaN, rms_dp = NaN, rms_dq = NaN, rows = rows)
  dp = [abs(r.dp) for r in rows]
  dq = [abs(r.dq) for r in rows]
  return (n = length(rows), max_dp = maximum(dp), max_dq = maximum(dq), rms_dp = sqrt(sum(abs2, dp) / length(dp)), rms_dq = sqrt(sum(abs2, dq) / length(dq)), rows = rows)
end

"""
    shortCircuitCoverage(sc::CGMESShortCircuitData) -> Vector{NamedTuple}

Per-class completeness of the harvested short-circuit source data: one row
per element class with the record count and, per attribute, how many records
actually carry a value (`attribute => filled/total`). Identification fields
(`mrid`, `name`, `bus`) are excluded — coverage describes the electrical
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
count and each attribute's fill rate, `✓` when complete. This is what
`cgmes.log` prints under "Short-circuit source data".
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
      marker = a.filled == a.total ? "✓" : (a.filled == 0 ? "—" : string(a.filled, "/", a.total))
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
    s.boundary_missing_hint && print(io, "\n  ⚠ topology references unresolved — boundary set likely missing")
  end
  println(io)
end
