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

# file: src/contingency.jl
# purpose: N-1 contingency batch API (Phase 4 of the multi-core work,
#          issue task_multicore_parallel): branch-outage cases generated
#          from the net or from imported FOR001 lists, solved on isolated
#          deep copies (optionally in parallel on Julia threads), evaluated
#          against voltage and loading limits, reported as a result table
#          and CSV. The base net is never mutated.

"""
    ContingencyCase

One contingency to evaluate: the outage of a single network element.

# Fields
- `name::String`: display name of the case (defaults to the element name).
- `kind::Symbol`: `:branch` (a line/transformer outage) or `:gen` (a generator
  outage).
- `element::String`: for `:branch`, the branch component name as reported by
  `getCompName(branch.comp)`, resolved against `net.branchVec`; for `:gen`, the
  generator prosumer's component name, resolved against `net.prosumpsVec`. Both
  can share a component name (parallel circuits, several units at one bus), so
  `generateN1Branches` / `generateN1Generators` disambiguate duplicates as
  `"<name>#<index>"` and the resolver verifies the index carries that name
  before use.
- `weight::Float64`: a per-case importance/probability weight (default `1.0`),
  carried through to the [`ContingencyResult`](@ref) so a weighted ranking can
  read it. Set it via the 4-argument constructor, or attach outage rates in bulk
  with [`applyContingencyWeights`](@ref) / [`readContingencyWeightsCSV`](@ref).
  A weight of `0` is a pure ranking weight, NOT a skip switch: the case is still
  solved and reported, and only a weighted ranking neutralizes it. To leave an
  outage out of the batch entirely, drop it at generation time with the
  `generateN1Branches` screening filters.
"""
struct ContingencyCase
  name::String
  kind::Symbol
  element::String
  weight::Float64
  function ContingencyCase(name::String, kind::Symbol, element::String, weight::Real = 1.0)
    kind in (:branch, :gen) || throw(ArgumentError("ContingencyCase: kind must be :branch or :gen, got :$(kind)."))
    (isfinite(weight) && weight >= 0.0) || throw(ArgumentError("ContingencyCase: weight must be a finite value >= 0, got $(weight)."))
    return new(name, kind, element, Float64(weight))
  end
end

ContingencyCase(element::String) = ContingencyCase(element, :branch, element)

"""
    OverloadRecord

One overloaded branch under a contingency, collected per case in
[`ContingencyResult`](@ref).`overloads` so a reporting consumer reads the
numbers directly instead of parsing a string.

# Fields
- `name::String`: the branch component name.
- `loading_pct::Float64`: post-contingency loading, `100 · s_MVA / sn_MVA`.
- `loading_base_pct::Float64`: the same branch's loading in the SOLVED base
  case (before the outage); `NaN` when the base did not converge or the branch
  was unrated at base. The reference that makes the post-outage value readable.
- `delta_pct::Float64`: `loading_pct - loading_base_pct` (a jump from 40 % to
  105 % is a very different event than 98 % to 105 %); `NaN` if the base is
  unknown.
- `s_MVA::Float64`: post-contingency apparent power, `max(|S_from|, |S_to|)`.
- `sn_MVA::Float64`: the branch rating the loading is measured against.
"""
struct OverloadRecord
  name::String
  loading_pct::Float64
  loading_base_pct::Float64
  delta_pct::Float64
  s_MVA::Float64
  sn_MVA::Float64
end

"""
    ContingencyResult

Outcome of one [`ContingencyCase`](@ref) evaluated by
[`runContingencies!`](@ref).

# Fields
- `name::String`: the case name.
- `weight::Float64`: the originating case's importance/probability weight
  (default `1.0`), carried through unchanged for weighted ranking.
- `converged::Bool`: whether the post-outage power flow converged.
- `iterations::Int`: Newton iterations of the contingency solve (0 when the
  solve never ran; the sum across ladder stages that were attempted).
- `start_used::Symbol`: the start-value ladder stage that converged
  (`:warm`, `:apslf`, `:dc`, or `:flat`), `:none` for a failed case.
- `max_vm_pu::Float64` / `min_vm_pu::Float64`: voltage envelope over the
  non-isolated buses of the solved case (NaN when not solved).
- `max_branch_loading_pct::Float64`: maximum loading over all rated
  branches (|S| at either end against `sn_MVA`); NaN when the case did not
  solve or no branch carries a finite rating.
- `severity::Float64`: a rankable severity, `weight · max(0, max loading - 100)`
  for a converged case (`0` when nothing is overloaded), and `NaN` for a failed
  case so the failures sort to the TOP. This is what
  [`printContingencyResults`](@ref) ranks by, giving `weight` its consumer.
- `overloads::Vector{OverloadRecord}`: the branches loaded above 100 percent,
  each with its loading and the delta to base, worst first.
- `voltage_violations::Vector{String}`: bus names outside the
  `[vm_min_pu, vm_max_pu]` band.
- `island_count::Int`: AC islands of the post-outage topology (0 when the
  net could not be evaluated).
- `shed_load_mw::Float64`: load disconnected by islanding (the total load in
  reference-less islands); `0.0` when the case solves or does not island.
- `error::Union{Nothing,String}`: `nothing` on success; otherwise the
  failure in one line ("islanded without reference", the solver status, or
  the exception message). Failures are REPORTED, never thrown.
"""
struct ContingencyResult
  name::String
  weight::Float64
  converged::Bool
  iterations::Int
  start_used::Symbol
  max_vm_pu::Float64
  min_vm_pu::Float64
  max_branch_loading_pct::Float64
  severity::Float64
  overloads::Vector{OverloadRecord}
  voltage_violations::Vector{String}
  island_count::Int
  shed_load_mw::Float64
  error::Union{Nothing,String}
end

"""
    generateN1Branches(net::Net; include_transformers = true,
                       min_vn_kV = 0.0, min_sn_MVA = 0.0, name_pattern = nothing)
        -> Vector{ContingencyCase}

Generate one branch-outage [`ContingencyCase`](@ref) per in-service branch
of `net` (aggregate `status == 1`; partially open branches are skipped,
their outage is already half-effective). With
`include_transformers = false` only AC-line branches are listed; a branch
counts as a transformer when its component type is `Trafo` (net-builder
path) or its winding ratio is nonzero (the classical MATPOWER indicator:
line rows carry `ratio = 0`).

Optional screening filters keep only the outages worth simulating on a large
grid (all default to "no filter"):
- `min_vn_kV`: keep a branch only if the HIGHER of its two endpoint voltages is
  at least this value (so a step-down transformer touching the EHV level is
  kept). Screens N-1 down to the transmission grid.
- `min_sn_MVA`: keep a branch only if it carries a finite `sn_MVA` rating of at
  least this value. Unrated branches are dropped when this filter is active.
- `name_pattern`: keep a branch only if its component name matches this
  `Regex` or contains this substring.
"""
function generateN1Branches(net::Net; include_transformers::Bool = true,
                            min_vn_kV::Real = 0.0, min_sn_MVA::Real = 0.0,
                            name_pattern::Union{Nothing,AbstractString,Regex} = nothing)::Vector{ContingencyCase}
  name_count = Dict{String,Int}()
  for br in net.branchVec
    name = getCompName(br.comp)
    name_count[name] = get(name_count, name, 0) + 1
  end
  cases = ContingencyCase[]
  for br in net.branchVec
    br.status == 1 || continue
    _branch_terminal_state(br) === :closed || continue
    is_transformer = br.comp.cTyp === Trafo || br.ratio != 0.0
    include_transformers || !is_transformer || continue
    # voltage-level screen: the higher endpoint voltage is the level the branch
    # belongs to (keeps EHV/HV transformers when screening for the top grid)
    if min_vn_kV > 0.0
      vn = max(getNodeVn(net.nodeVec[Int(br.fromBus)]), getNodeVn(net.nodeVec[Int(br.toBus)]))
      vn >= min_vn_kV || continue
    end
    # rating screen: only branches with a usable rating; an unrated branch
    # cannot satisfy a positive threshold, so it is dropped
    if min_sn_MVA > 0.0
      rating = br.sn_MVA
      (rating !== nothing && isfinite(rating) && rating >= min_sn_MVA) || continue
    end
    name = getCompName(br.comp)
    (name_pattern === nothing || occursin(name_pattern, name)) || continue
    # parallel circuits share a component name; disambiguate by branch index
    # so every circuit gets its OWN outage case (a bare name would always
    # resolve to the first circuit and evaluate the same outage twice)
    element = name_count[name] > 1 ? string(name, "#", br.branchIdx) : name
    push!(cases, ContingencyCase(element, :branch, element))
  end
  return cases
end

"""
    generateN1Generators(net::Net; min_pg_MW = 0.0, name_pattern = nothing)
        -> Vector{ContingencyCase}

Generate one generator-outage [`ContingencyCase`](@ref) (`kind = :gen`) per
in-service generator prosumer of `net` (`isGenerator`, i.e. an injection: a
generator, external-grid feed-in, or synchronous machine). A generator outage
removes ONLY that one unit's injection; the topology is unchanged, so the lost
active power must be picked up by the slack (or, with
`distributed_slack_enabled = true` on [`runContingencies!`](@ref), shared over
the surviving participants). Removing a bus's last voltage-regulating source
demotes it to PQ, and removing the island's only reference is reported as
`islanded without reference` with the stranded generation named.

Optional screening filters (default to "no filter"):
- `min_pg_MW`: keep a generator only if `|Pg|` is at least this value (skips
  tiny or zero-output units).
- `name_pattern`: keep names matching a `Regex` or containing a substring.

Generator component names are not unique (several units can share one name or
bus); duplicates are disambiguated as `"<name>#<prosumerIndex>"`.
"""
function generateN1Generators(net::Net; min_pg_MW::Real = 0.0,
                              name_pattern::Union{Nothing,AbstractString,Regex} = nothing)::Vector{ContingencyCase}
  name_count = Dict{String,Int}()
  for ps in net.prosumpsVec
    isGenerator(ps) || continue
    name = getCompName(ps.comp)
    name_count[name] = get(name_count, name, 0) + 1
  end
  cases = ContingencyCase[]
  for (pidx, ps) in enumerate(net.prosumpsVec)
    isGenerator(ps) || continue
    if min_pg_MW > 0.0
      abs(something(ps.pVal, 0.0)) >= min_pg_MW || continue
    end
    name = getCompName(ps.comp)
    (name_pattern === nothing || occursin(name_pattern, name)) || continue
    # generator names are not unique (unlike branch names); disambiguate by the
    # prosumer index so every unit resolves to its OWN outage
    element = name_count[name] > 1 ? string(name, "#", pidx) : name
    push!(cases, ContingencyCase(element, :gen, element))
  end
  return cases
end

# resolve a case element against branchVec: either a plain component name
# (first match) or the disambiguated "<name>#<branchIdx>" form, where the
# index must carry the expected name (guards against stale case lists)
function _resolve_contingency_branch(cnet::Net, element::String)::Union{Nothing,Int}
  hash_pos = findlast('#', element)
  if hash_pos !== nothing
    idx = tryparse(Int, element[(hash_pos+1):end])
    name = element[1:(hash_pos-1)]
    if idx !== nothing && 1 <= idx <= length(cnet.branchVec) && getCompName(cnet.branchVec[idx].comp) == name
      return idx
    end
    return nothing
  end
  return findfirst(b -> getCompName(b.comp) == element, cnet.branchVec)
end

# resolve a generator case element against prosumpsVec: a plain component name
# (first generator match) or the disambiguated "<name>#<prosumerIndex>" form,
# where the index must be a generator carrying the expected name
function _resolve_contingency_generator(cnet::Net, element::String)::Union{Nothing,Int}
  hash_pos = findlast('#', element)
  if hash_pos !== nothing
    idx = tryparse(Int, element[(hash_pos+1):end])
    name = element[1:(hash_pos-1)]
    if idx !== nothing && 1 <= idx <= length(cnet.prosumpsVec) && isGenerator(cnet.prosumpsVec[idx]) && getCompName(cnet.prosumpsVec[idx].comp) == name
      return idx
    end
    return nothing
  end
  return findfirst(ps -> isGenerator(ps) && getCompName(ps.comp) == element, cnet.prosumpsVec)
end

# take one generator out of service on a case-local net: delete just that
# prosumer. buildComplexSVec reads generation straight from prosumpsVec, so the
# delete already removes this unit's injection from the solve. node._pƩGen is
# what detect_ac_islands reads for the stranded-generation classification, but
# the solved template holds SOLVED dispatch there (a reference bus carries its
# balancing injection, not its schedule), so subtracting the removed unit's
# SCHEDULED P would give a nonsense (even negative) figure. Recompute the bus
# aggregate from the SURVIVING generators' scheduled P/Q instead, so the island
# report reflects the generation actually stranded. Then rebuild bus types and
# Q-limits so the bus drops PV/Slack -> PQ when no regulating or slack unit
# survives (mirrors removeProsumer!'s tail).
#
# CALL AT MOST ONCE PER CASE. `deleteat!` shifts every later prosumer index, and
# `generateN1Generators` disambiguates units by exactly that index
# (`"<name>#<prosumerIndex>"`), which `_resolve_contingency_generator` resolves
# BEFORE the delete. A single removal per case-local deepcopy is therefore safe,
# but a second index-based removal on the same net would hit the wrong unit. Any
# future multi-element case (e.g. a gen-plus-branch combination) must re-resolve
# by name after each removal rather than reuse a pre-computed index.
function _remove_contingency_generator!(cnet::Net, pidx::Int)
  busIdx = getPosumerBusIndex(cnet.prosumpsVec[pidx])
  deleteat!(cnet.prosumpsVec, pidx)
  node = cnet.nodeVec[busIdx]
  node._pƩGen = 0.0
  node._qƩGen = 0.0
  for ps in cnet.prosumpsVec
    (isGenerator(ps) && getPosumerBusIndex(ps) == busIdx) || continue
    node._pƩGen += something(ps.pVal, 0.0)
    node._qƩGen += something(ps.qVal, 0.0)
  end
  refreshBusTypesFromProsumers!(cnet)
  _buildQLimits!(cnet)
  return nothing
end

"""
    generateContingenciesFromFOR001(net::Net) -> Vector{ContingencyCase}

Map the FOR001 contingency branch names imported from MATPOWER metadata
(`mpc.for001_contingencies`, stored on `net.for001Contingencies`) to
[`ContingencyCase`](@ref)s. Names that do not resolve to a branch are still
listed; [`runContingencies!`](@ref) reports them as failed cases with an
actionable error instead of dropping them silently.
"""
generateContingenciesFromFOR001(net::Net)::Vector{ContingencyCase} = [ContingencyCase(name, :branch, name) for name in net.for001Contingencies]

"""
    applyContingencyWeights(cases, weights; default = 1.0) -> Vector{ContingencyCase}

Return a copy of `cases` with each case's `weight` set from `weights`
(an `AbstractDict` mapping case name to a non-negative number), looked up by
`case.name`. Cases whose name is absent from `weights` get `default`. Useful to
attach per-branch outage rates read with [`readContingencyWeightsCSV`](@ref)
before calling [`runContingencies!`](@ref); [`ContingencyCase`](@ref) is
immutable, so this builds new cases rather than mutating in place.
"""
function applyContingencyWeights(cases::AbstractVector{ContingencyCase}, weights::AbstractDict; default::Real = 1.0)::Vector{ContingencyCase}
  return ContingencyCase[ContingencyCase(c.name, c.kind, c.element, Float64(get(weights, c.name, default))) for c in cases]
end

"""
    readContingencyWeightsCSV(path) -> Dict{String,Float64}

Read per-case outage weights from a two-column CSV mapping case name to weight.
The delimiter is auto-detected as `;` or `,` (semicolon wins when both appear).
A header line naming the columns (e.g. `name;weight`) is recognized and skipped;
blank lines and lines beginning with `#` are ignored. Feed the result to
[`applyContingencyWeights`](@ref). Throws an `ArgumentError` on a malformed row
or a negative / non-finite weight.
"""
function readContingencyWeightsCSV(path::AbstractString)::Dict{String,Float64}
  weights = Dict{String,Float64}()
  for (lineno, raw) in enumerate(eachline(path))
    line = strip(raw)
    (isempty(line) || startswith(line, '#')) && continue
    delim = occursin(';', line) ? ';' : ','
    fields = split(line, delim)
    length(fields) >= 2 || throw(ArgumentError("readContingencyWeightsCSV: line $(lineno) needs a name and a weight column, got: $(raw)"))
    name = strip(fields[1])
    wtext = strip(fields[2])
    # skip a header row: a non-numeric second field on the first data line
    w = tryparse(Float64, wtext)
    if w === nothing
      lineno == 1 && continue          # header like "name;weight"
      throw(ArgumentError("readContingencyWeightsCSV: line $(lineno) has a non-numeric weight \"$(wtext)\"."))
    end
    (isfinite(w) && w >= 0.0) || throw(ArgumentError("readContingencyWeightsCSV: line $(lineno) weight must be finite and >= 0, got $(w)."))
    weights[String(name)] = w
  end
  return weights
end

# evaluate the solved contingency net against the limit band; pure reader
function _contingency_metrics(cnet::Net, vm_min_pu::Float64, vm_max_pu::Float64, base_loadings::Dict{String,Float64})
  iso = Set(cnet.isoNodes)
  bus_names = Dict{Int,String}(idx => name for (name, idx) in cnet.busDict)
  vmin = Inf
  vmax = -Inf
  violations = String[]
  for node in cnet.nodeVec
    node.busIdx in iso && continue
    vm = something(node._vm_pu, NaN)
    isfinite(vm) || continue
    vmin = min(vmin, vm)
    vmax = max(vmax, vm)
    if vm < vm_min_pu || vm > vm_max_pu
      push!(violations, get(bus_names, node.busIdx, getCompName(node.comp)))
    end
  end
  max_loading = NaN
  overloads = OverloadRecord[]
  for br in cnet.branchVec
    br.status == 1 || continue
    rating = br.sn_MVA
    (rating === nothing || !isfinite(rating) || rating <= 0.0) && continue
    s_from = br.fBranchFlow === nothing ? 0.0 : hypot(something(br.fBranchFlow.pFlow, 0.0), something(br.fBranchFlow.qFlow, 0.0))
    s_to = br.tBranchFlow === nothing ? 0.0 : hypot(something(br.tBranchFlow.pFlow, 0.0), something(br.tBranchFlow.qFlow, 0.0))
    s_mva = max(s_from, s_to)
    loading = 100.0 * s_mva / rating
    (isnan(max_loading) || loading > max_loading) && (max_loading = loading)
    if loading > 100.0
      name = getCompName(br.comp)
      # base loading is looked up by name; parallel circuits share a name, so
      # their base value is an approximation (their base loadings are close)
      base = get(base_loadings, name, NaN)
      push!(overloads, OverloadRecord(name, loading, base, loading - base, s_mva, rating))
    end
  end
  # worst-loaded branch first, so the result list and the report read top-down
  sort!(overloads; by = o -> -o.loading_pct)
  return (; vmin = isfinite(vmin) ? vmin : NaN, vmax = isfinite(vmax) ? vmax : NaN, violations, max_loading, overloads)
end

# base-case loading per rated branch (component name -> loading_pct), computed
# ONCE on the solved template so every contingency can report the delta to the
# pre-outage state. Empty when the base did not converge (deltas are then NaN).
function _base_branch_loadings(net::Net)::Dict{String,Float64}
  loadings = Dict{String,Float64}()
  for br in net.branchVec
    br.status == 1 || continue
    rating = br.sn_MVA
    (rating === nothing || !isfinite(rating) || rating <= 0.0) && continue
    s_from = br.fBranchFlow === nothing ? 0.0 : hypot(something(br.fBranchFlow.pFlow, 0.0), something(br.fBranchFlow.qFlow, 0.0))
    s_to = br.tBranchFlow === nothing ? 0.0 : hypot(something(br.tBranchFlow.pFlow, 0.0), something(br.tBranchFlow.qFlow, 0.0))
    loadings[getCompName(br.comp)] = 100.0 * max(s_from, s_to) / rating
  end
  return loadings
end

# The per-case start-value ladder (#331 Phase 1). Each stage is ONE bounded
# solve on the case-local net with a distinct start recipe; the first stage
# that converges wins and is reported as `start_used`. Ordered subset of these.
const CONTINGENCY_RESCUE_LADDER_VALUES = (:warm, :apslf, :dc, :flat)

"""
    _validate_contingency_ladder(ladder; context) -> Vector{Symbol}

Validate a contingency start-value ladder: a non-empty, duplicate-free ordered
subset of `(:warm, :apslf, :dc, :flat)`. Throws an `ArgumentError` naming
`context` otherwise. Shared by the `contingency.rescue_ladder` config section
and the `runContingencies!` keyword so both reject the same inputs.
"""
function _validate_contingency_ladder(ladder::AbstractVector{Symbol}; context::AbstractString = "contingency.rescue_ladder")::Vector{Symbol}
  isempty(ladder) && throw(ArgumentError("$(context) must be a non-empty ordered subset of $(collect(CONTINGENCY_RESCUE_LADDER_VALUES)); got an empty list."))
  for s in ladder
    s in CONTINGENCY_RESCUE_LADDER_VALUES || throw(ArgumentError("$(context): unknown stage :$(s); allowed stages are $(collect(CONTINGENCY_RESCUE_LADDER_VALUES))."))
  end
  allunique(ladder) || throw(ArgumentError("$(context) must not contain duplicate stages; got $(ladder)."))
  return Symbol[ladder...]
end

# Build the result row for an outage that strands one or more islands without a
# valid angle reference. The message is deliberately specific so the operator
# reads the cause, not just the symptom: a load-only island is a clean load-shed
# outcome (no generation to reference at all), whereas an island that carries
# generation but no voltage-controlled source is the rarer stranded-generation
# case. Both keep the leading "islanded" token so result filtering by island
# stays simple. The MW figures are the total scheduled load / generation summed
# over the reference-less island(s) (`refless` are the detect_ac_islands rows
# with chosen_ref_bus == 0).
function _islanded_contingency_result(name::AbstractString, weight::Float64, island_count::Int, refless::AbstractVector)
  # fallback: the solver reported islanding the pre-merge snapshot did not see
  # (no refless rows to size), so report the cause generically rather than 0 MW
  isempty(refless) && return ContingencyResult(name, weight, false, 0, :none, NaN, NaN, NaN, NaN, OverloadRecord[], String[], island_count, 0.0, "islanded without reference")
  shed_mw = sum(r.total_load_p_mw for r in refless)
  stranded_gen_mw = sum(r.total_gen_p_mw for r in refless)
  msg = if stranded_gen_mw > 1e-6
    "islanded without reference: $(round(shed_mw; digits = 1)) MW load, $(round(stranded_gen_mw; digits = 1)) MW generation stranded (no voltage-controlled source)"
  else
    "islanded: load-only, $(round(shed_mw; digits = 1)) MW load disconnected"
  end
  return ContingencyResult(name, weight, false, 0, :none, NaN, NaN, NaN, NaN, OverloadRecord[], String[], island_count, shed_mw, msg)
end

# one worker evaluation: template copy -> outage -> ladder of start recipes ->
# metrics. Pure per-case function: touches no shared container, all mutation
# happens on the case-local deepcopy. The outage topology (branch removed,
# islands marked) is fixed ONCE; only the start voltages vary across stages,
# restored from a per-case snapshot before each stage.
function _run_one_contingency(template::Net, case::ContingencyCase, vm_min_pu::Float64, vm_max_pu::Float64, maxIte::Int, tol::Float64, ladder::Vector{Symbol}, pf_kwargs, base_loadings::Dict{String,Float64})
  cnet = deepcopy(template)
  # apply the outage: a branch removal (topology change) or a generator removal
  # (injection change, topology unchanged). Both then go through the same island
  # detection, start-value ladder, and metrics.
  if case.kind === :gen
    pidx = _resolve_contingency_generator(cnet, case.element)
    if pidx === nothing
      return ContingencyResult(case.name, case.weight, false, 0, :none, NaN, NaN, NaN, NaN, OverloadRecord[], String[], 0, 0.0, "unknown generator $(case.element) (not found in net.prosumpsVec)")
    end
    _remove_contingency_generator!(cnet, pidx)
  else
    idx = _resolve_contingency_branch(cnet, case.element)
    if idx === nothing
      return ContingencyResult(case.name, case.weight, false, 0, :none, NaN, NaN, NaN, NaN, OverloadRecord[], String[], 0, 0.0, "unknown branch $(case.element) (not found in net.branchVec)")
    end
    removeBranch!(net = cnet, branchNr = idx)
  end
  markIsolatedBuses!(net = cnet, log = false)
  island_report = detect_ac_islands(cnet)
  island_count = length(island_report.rows)
  # Rows for reference-less islands (chosen_ref_bus == 0), kept so the islanding
  # error below can be reported SPECIFICALLY (load-only vs. stranded generation).
  # We do NOT pre-empt the solve on this: the solver is the authority on whether
  # a case is truly unsolvable, because it may link-merge synchronously-tied
  # islands (e.g. grid-forming HVDC) and reference an AC component that looks
  # reference-less on this pre-merge snapshot. The ladder's first stage throws
  # the validation error immediately, so letting it run costs nothing.
  refless = [r for r in island_report.rows if r.chosen_ref_bus == 0]
  # snapshot the warm (template) start; the :flat template fallback keeps its
  # flatstart flag, so :warm must NOT clear it (bitwise with the pre-#331 path)
  snap = _snapshot_start_voltages(cnet)
  template_flatstart = cnet.flatstart
  total_it = 0
  last_error = nothing
  for stage in ladder
    _restore_start_voltages!(cnet, snap)
    it_stage = 0
    try
      erg = 1
      if stage === :warm
        cnet.flatstart = template_flatstart
        it_stage, erg = runpf!(cnet, maxIte, tol, 0; islands_enabled = true, pf_kwargs...)
      elseif stage === :flat
        cnet.flatstart = true
        it_stage, erg = runpf!(cnet, maxIte, tol, 0; islands_enabled = true, pf_kwargs...)
      elseif stage === :dc
        # flat magnitudes with DC-projected start angles (same recipe the
        # solver rescue ladder's :dc_seed variant reaches)
        cnet.flatstart = false
        _dc_seed_rectangular_angles!(cnet, PowerFlowConfig(max_iter = maxIte, tol = tol))
        it_stage, erg = runpf!(cnet, maxIte, tol, 0; islands_enabled = true, pf_kwargs...)
      elseif stage === :apslf
        # APSLF start via the config-driven solve (rescue OFF: one bounded
        # attempt, not the full solver ladder); pf_kwargs are not forwarded on
        # this path (islands handled by the config)
        cnet.flatstart = false
        it_stage, erg = runpf!(cnet, PowerFlowConfig(rescue = false, islands_enabled = true, max_iter = maxIte, tol = tol, apslf_start = ApslfStartConfig(enabled = true, order = 40)))
      end
      total_it += it_stage
      if erg == 0
        calcNetLosses!(cnet)
        m = _contingency_metrics(cnet, vm_min_pu, vm_max_pu, base_loadings)
        # severity ranks converged cases by overload weighted by outage weight;
        # an un-overloaded case scores 0, and a case without any rating scores 0
        severity = case.weight * max(0.0, isnan(m.max_loading) ? 0.0 : m.max_loading - 100.0)
        return ContingencyResult(case.name, case.weight, true, total_it, stage, m.vmax, m.vmin, m.max_loading, severity, m.overloads, m.violations, island_count, 0.0, nothing)
      end
      last_error = "power flow did not converge (status $(erg))"
    catch err
      err isa InterruptException && rethrow(err)
      total_it += it_stage
      msg = sprint(showerror, err)
      # a removal that splits off a reference-less island surfaces as the
      # island-validation error; no start recipe fixes a topology fact, so
      # stop the ladder and report the specific cause (load-only vs. stranded
      # generation), sized from the pre-computed reference-less island rows
      if occursin("reference", msg) && (occursin("island", msg) || occursin("Island", msg))
        return _islanded_contingency_result(case.name, case.weight, island_count, refless)
      end
      last_error = first(split(msg, '\n'))
    end
  end
  return ContingencyResult(case.name, case.weight, false, total_it, :none, NaN, NaN, NaN, NaN, OverloadRecord[], String[], island_count, 0.0, last_error === nothing ? "power flow did not converge" : last_error)
end

"""
    runContingencies!(net::Net, cases::Vector{ContingencyCase};
                      vm_min_pu = 0.9, vm_max_pu = 1.1,
                      maxIte = 30, tol = 1e-8,
                      rescue_ladder = [:warm],
                      parallel_enabled = nothing, parallel_max_tasks = nothing,
                      parallel_min_work_items = nothing,
                      kwargs...) -> Vector{ContingencyResult}

Evaluate branch-outage contingencies. The base `net` is NEVER mutated: a
solved template copy is created first (the warm start), and each case works on
its own `deepcopy` of that template, removes its branch via
[`removeBranch!`](@ref), solves, and evaluates the `[vm_min_pu, vm_max_pu]`
band plus branch loadings against `sn_MVA` ratings. Cases are returned in
input order; failures (no convergence, islanding without reference, unknown
element) are REPORTED in the result, never thrown.

`rescue_ladder` is the per-case start-value ladder (#331): an ordered,
duplicate-free subset of `(:warm, :apslf, :dc, :flat)`, default `[:warm]`
(the pre-#331 single warm solve). Each stage is ONE bounded `runpf!` on the
case-local net with a distinct start recipe, tried in order until one
converges; the winning stage is reported as `start_used`:
- `:warm`: the template (base-case) voltages;
- `:apslf`: an APSLF start (needs `using AnalyticLoadFlow`; the stage is
  dropped with a warning when the extension is not loaded);
- `:dc`: flat magnitudes with DC-projected start angles;
- `:flat`: `flatstart = true`.
The full solver rescue ladder (`:settled_qlimits`, `:autodamp`) is NOT used
per case, only for the base case (below).

When the base case does not converge, it is retried through the solver rescue
ladder (`runpf!` with `rescue = true`) before falling back to a flat template
with a warning. Remaining `kwargs...` are forwarded to the `:warm`/`:flat`/`:dc`
contingency solves (the `:apslf` config path does not forward them).

`retry_flat_start` is DEPRECATED (kept one minor cycle): `retry_flat_start =
true` now just appends `:flat` to the ladder.

The batch fans out over Julia threads in `runtime.parallel.max_tasks`
chunks (gated by `runtime.parallel.*`; the three `parallel_*` keywords
override the active configuration). Parallel and serial runs produce
identical results.
"""
function runContingencies!(
  net::Net,
  cases::Vector{ContingencyCase};
  vm_min_pu::Float64 = 0.9,
  vm_max_pu::Float64 = 1.1,
  maxIte::Int = 30,
  tol::Float64 = 1e-8,
  rescue_ladder::Vector{Symbol} = [:warm],
  retry_flat_start::Union{Nothing,Bool} = nothing,
  parallel_enabled::Union{Nothing,Bool} = nothing,
  parallel_max_tasks::Union{Nothing,Int} = nothing,
  parallel_min_work_items::Union{Nothing,Int} = nothing,
  kwargs...,
)::Vector{ContingencyResult}
  vm_min_pu < vm_max_pu || throw(ArgumentError("runContingencies!: vm_min_pu must be below vm_max_pu."))
  isempty(cases) && return ContingencyResult[]

  # resolve the per-case start-value ladder (#331 Phase 1)
  ladder = _validate_contingency_ladder(rescue_ladder; context = "runContingencies!: rescue_ladder")
  # retry_flat_start is DEPRECATED: it is now an alias for appending :flat to
  # the ladder. Kept one minor cycle; remove after 0.9.x.
  if retry_flat_start !== nothing
    @warn "runContingencies!: the keyword retry_flat_start is deprecated; pass rescue_ladder instead (retry_flat_start = true is now an alias for appending :flat to the ladder)."
    retry_flat_start && !(:flat in ladder) && push!(ladder, :flat)
  end
  # :apslf needs the AnalyticLoadFlow.jl extension loaded; drop it with one
  # warning rather than hard-erroring per case when it is absent
  if :apslf in ladder && Base.get_extension(Sparlectra, :SparlectraAnalyticLoadFlowExt) === nothing
    @warn "runContingencies!: the :apslf ladder stage needs AnalyticLoadFlow.jl (`using AnalyticLoadFlow`); it is not loaded, so :apslf is dropped from the ladder."
    ladder = filter(!=(:apslf), ladder)
    isempty(ladder) && (ladder = [:warm])
  end

  # warm-start template: the solved base case. The solve happens on a COPY,
  # so the caller's net stays untouched.
  template = deepcopy(net)
  base_converged = try
    _, base_erg = runpf!(template, maxIte, tol, 0; islands_enabled = true, kwargs...)
    base_erg == 0
  catch err
    err isa InterruptException && rethrow(err)
    false
  end
  if !base_converged
    # #331 Phase 1 item 4: give the base case the full solver rescue ladder
    # (alternate start, autodamp, settled Q-limits, DC seed) before deciding it
    # is unsolvable; only then fall back to the flat template.
    rescued = deepcopy(net)
    base_converged = try
      _, r_erg = runpf!(rescued, PowerFlowConfig(rescue = true, islands_enabled = true, max_iter = maxIte, tol = tol))
      r_erg == 0
    catch err
      err isa InterruptException && rethrow(err)
      false
    end
    if base_converged
      template = rescued
    else
      @warn "runContingencies!: the base case did not converge even with the solver rescue ladder; contingencies start FLAT instead of warm."
      template = deepcopy(net)
      template.flatstart = true
    end
  end
  # template hygiene (follow-up item 1): the workers need the solved
  # VOLTAGES, not the base solver status or its Q-limit event logs. Clearing
  # them keeps every per-case deepcopy lean (and avoids the one-time
  # deepcopy compile of the status closure types on the first worker).
  template._rectangular_pf_status = nothing
  template._dc_pf_status = nothing
  empty!(template.qLimitLog)
  empty!(template.qLimitEvents)
  empty!(template.qLimitInitialPVRows)

  # base-case branch loadings, computed ONCE, so every contingency can report
  # the delta to the pre-outage loading. Only meaningful when the base actually
  # solved; a flat-fallback template has no valid flows, so leave deltas NaN.
  base_loadings = Dict{String,Float64}()
  if base_converged
    calcNetLosses!(template)
    base_loadings = _base_branch_loadings(template)
  end

  pf_kwargs = kwargs
  results = Vector{Any}(undef, length(cases))
  parallel_on, parallel_cap, parallel_min_items = _resolve_parallel_runtime(parallel_enabled, parallel_max_tasks, parallel_min_work_items)
  use_parallel = parallel_on && Threads.nthreads() > 1 && parallel_cap > 1 && length(cases) >= parallel_min_items

  if use_parallel
    chunks = collect(Iterators.partition(eachindex(cases), cld(length(cases), parallel_cap)))
    tasks = [
      Threads.@spawn begin
        for idx in chunks[ci]
          results[idx] = _run_one_contingency(template, cases[idx], vm_min_pu, vm_max_pu, maxIte, tol, ladder, pf_kwargs, base_loadings)
        end
      end for ci in eachindex(chunks)
    ]
    foreach(wait, tasks)
  else
    for idx in eachindex(cases)
      results[idx] = _run_one_contingency(template, cases[idx], vm_min_pu, vm_max_pu, maxIte, tol, ladder, pf_kwargs, base_loadings)
    end
  end
  return ContingencyResult[results[i] for i in eachindex(cases)]
end

"""
    printContingencyResults([io::IO], results::Vector{ContingencyResult}; max_rows = 50)

Print the contingency results as a fixed-width table (style of
`printShortCircuitResult`): case, convergence, iterations, voltage
envelope, worst loading, violation counts, shed load, weight, severity, and
the error line for failed cases. Rows beyond `max_rows` are summarized in one
line.

By default (`sort_by = :severity`) the rows are ranked so the failed cases come
first and the converged cases follow by descending `severity`
(`weight · max(0, loading - 100)`); this keeps the worst contingency on the
first page of a long list. Pass `sort_by = :none` to keep the input order (the
order `writeContingencyResultsCSV` always uses).
"""
printContingencyResults(results::Vector{ContingencyResult}; max_rows::Int = 50, sort_by::Symbol = :severity) = printContingencyResults(stdout, results; max_rows = max_rows, sort_by = sort_by)

function printContingencyResults(io::IO, results::Vector{ContingencyResult}; max_rows::Int = 50, sort_by::Symbol = :severity)
  ordered = if sort_by === :severity
    # failed cases (NaN severity) first, then converged by severity descending
    sort(results; by = r -> (isnan(r.severity) ? -Inf : -r.severity))
  elseif sort_by === :none
    results
  else
    throw(ArgumentError("printContingencyResults: sort_by must be :severity or :none, got :$(sort_by)."))
  end
  println(io, "N-1 contingency results (", length(results), " case(s), ", count(r -> r.converged, results), " converged)")
  println(io, "-"^151)
  @printf(io, "%-28s %-9s %5s %-6s %9s %9s %12s %8s %8s %7s %8s %7s %8s  %s\n", "case", "converged", "iter", "start", "Vmin[pu]", "Vmax[pu]", "loading[%]", "overld", "V-viol", "islands", "shed[MW]", "weight", "severity", "error")
  println(io, "-"^151)
  shown = 0
  for r in ordered
    shown >= max_rows && break
    shown += 1
    @printf(
      io,
      "%-28s %-9s %5d %-6s %9s %9s %12s %8d %8d %7d %8s %7s %8s  %s\n",
      _fitColumn(r.name, 28),
      r.converged ? "yes" : "NO",
      r.iterations,
      String(r.start_used),
      isnan(r.min_vm_pu) ? "-" : @sprintf("%.4f", r.min_vm_pu),
      isnan(r.max_vm_pu) ? "-" : @sprintf("%.4f", r.max_vm_pu),
      isnan(r.max_branch_loading_pct) ? "-" : @sprintf("%.1f", r.max_branch_loading_pct),
      length(r.overloads),
      length(r.voltage_violations),
      r.island_count,
      r.shed_load_mw > 0.0 ? @sprintf("%.1f", r.shed_load_mw) : "-",
      @sprintf("%.2f", r.weight),
      isnan(r.severity) ? "-" : @sprintf("%.2f", r.severity),
      r.error === nothing ? "" : r.error,
    )
  end
  shown < length(ordered) && println(io, "... ", length(ordered) - shown, " more row(s) not shown (max_rows = ", max_rows, ")")
  println(io, "-"^151)
  return nothing
end

"""
    writeContingencyResultsCSV(path::AbstractString, results::Vector{ContingencyResult})

Write the contingency results as a semicolon-separated CSV (one row per
case, always in input order): name, weight, converged, iterations, the
start-value ladder stage that converged, voltage envelope, worst loading,
severity, the overloaded branches (each `name@loading%`, the per-branch delta to
base lives in the structured `overloads` field), voltage-violation list, island
count, shed load in MW, and the error text. Returns `path`.
"""
function writeContingencyResultsCSV(path::AbstractString, results::Vector{ContingencyResult})
  open(path, "w") do io
    println(io, "name;weight;converged;iterations;start_used;min_vm_pu;max_vm_pu;max_branch_loading_pct;severity;overloads;voltage_violations;island_count;shed_load_mw;error")
    for r in results
      overloads = join(["$(o.name)@$(round(o.loading_pct; digits = 1))" for o in r.overloads], ",")
      println(io, join([r.name, r.weight, r.converged, r.iterations, r.start_used, r.min_vm_pu, r.max_vm_pu, r.max_branch_loading_pct, r.severity, overloads, join(r.voltage_violations, ","), r.island_count, r.shed_load_mw, r.error === nothing ? "" : r.error], ";"))
    end
  end
  return path
end

"""
    ContingencyReport

Aggregate view over a contingency batch (built by
[`buildContingencyReport`](@ref)): case counts by outcome, the total and worst
load shed, the single worst branch loading, the worst weighted `severity`, and
the branches overloaded by the most contingencies. Print it with
[`printContingencyReport`](@ref).
"""
struct ContingencyReport
  n_cases::Int
  n_converged::Int
  n_islanded::Int
  n_nonconverged::Int
  n_with_overload::Int
  n_with_voltage_violation::Int
  total_shed_load_mw::Float64
  worst_shed_case::String
  worst_shed_mw::Float64
  worst_loading_case::String
  worst_loading_branch::String
  worst_loading_pct::Float64
  worst_severity_case::String
  worst_severity::Float64
  top_overloaded::Vector{Pair{String,Int}}
end

# an islanded-without-reference outcome (branch or generator) carries "island"
# in its error line; used to separate topological load shed from plain divergence
_is_islanded_result(r::ContingencyResult) = r.error !== nothing && occursin("island", r.error)

"""
    buildContingencyReport(results; top = 10) -> ContingencyReport

Summarize a [`runContingencies!`](@ref) batch: count cases by outcome, the total
and worst load shed, the worst single branch loading across all cases, and the
`top` branches overloaded by the most contingencies (ties broken by name).
"""
function buildContingencyReport(results::AbstractVector{ContingencyResult}; top::Int = 10)::ContingencyReport
  n_converged = count(r -> r.converged, results)
  n_islanded = count(_is_islanded_result, results)
  n_nonconverged = count(r -> !r.converged && !_is_islanded_result(r), results)
  n_with_overload = count(r -> !isempty(r.overloads), results)
  n_with_voltage_violation = count(r -> !isempty(r.voltage_violations), results)
  total_shed = sum(r -> r.shed_load_mw, results; init = 0.0)

  worst_shed_case = ""
  worst_shed_mw = 0.0
  for r in results
    if r.shed_load_mw > worst_shed_mw
      worst_shed_mw = r.shed_load_mw
      worst_shed_case = r.name
    end
  end

  worst_loading_case = ""
  worst_loading_branch = ""
  worst_loading_pct = NaN
  for r in results, o in r.overloads
    if isnan(worst_loading_pct) || o.loading_pct > worst_loading_pct
      worst_loading_pct = o.loading_pct
      worst_loading_branch = o.name
      worst_loading_case = r.name
    end
  end

  # worst weighted severity among the converged cases (failed cases have NaN
  # severity; they are counted above, not ranked here)
  worst_severity_case = ""
  worst_severity = 0.0
  for r in results
    if !isnan(r.severity) && r.severity > worst_severity
      worst_severity = r.severity
      worst_severity_case = r.name
    end
  end

  # how many distinct contingencies overload each branch
  counts = Dict{String,Int}()
  for r in results, o in r.overloads
    counts[o.name] = get(counts, o.name, 0) + 1
  end
  ranked = sort!(collect(counts); by = p -> (-p.second, p.first))
  top_overloaded = ranked[1:min(top, length(ranked))]

  return ContingencyReport(length(results), n_converged, n_islanded, n_nonconverged, n_with_overload, n_with_voltage_violation, total_shed, worst_shed_case, worst_shed_mw, worst_loading_case, worst_loading_branch, worst_loading_pct, worst_severity_case, worst_severity, top_overloaded)
end

printContingencyReport(report::ContingencyReport) = printContingencyReport(stdout, report)

"""
    printContingencyReport([io], report::ContingencyReport)

Print the aggregate [`ContingencyReport`](@ref) as a compact summary block:
outcome counts, load shed, worst loading, and the most-overloaded branches.
"""
function printContingencyReport(io::IO, report::ContingencyReport)
  println(io, "N-1 contingency report (", report.n_cases, " case(s))")
  println(io, "-"^60)
  @printf(io, "  converged              : %d\n", report.n_converged)
  @printf(io, "  islanded (load shed)   : %d\n", report.n_islanded)
  @printf(io, "  non-converged          : %d\n", report.n_nonconverged)
  @printf(io, "  with overload          : %d\n", report.n_with_overload)
  @printf(io, "  with voltage violation : %d\n", report.n_with_voltage_violation)
  @printf(io, "  total load shed        : %.1f MW\n", report.total_shed_load_mw)
  if report.worst_shed_mw > 0.0
    @printf(io, "  worst load shed        : %.1f MW  (%s)\n", report.worst_shed_mw, report.worst_shed_case)
  end
  if !isnan(report.worst_loading_pct)
    @printf(io, "  worst branch loading   : %.1f%%  (%s in %s)\n", report.worst_loading_pct, report.worst_loading_branch, report.worst_loading_case)
  end
  if report.worst_severity > 0.0
    @printf(io, "  worst severity         : %.2f  (%s)\n", report.worst_severity, report.worst_severity_case)
  end
  if !isempty(report.top_overloaded)
    println(io, "  most-overloaded branches:")
    for (name, cnt) in report.top_overloaded
      @printf(io, "    %-28s %d contingenc%s\n", _fitColumn(name, 28), cnt, cnt == 1 ? "y" : "ies")
    end
  end
  println(io, "-"^60)
  return nothing
end
