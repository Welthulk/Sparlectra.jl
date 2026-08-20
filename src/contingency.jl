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
- `kind::Symbol`: only `:branch` in this version.
- `element::String`: the branch's component name as reported by
  `getCompName(branch.comp)`; resolved against `net.branchVec` at run time.
  Parallel circuits can share one component name; `generateN1Branches`
  disambiguates duplicates as `"<name>#<branchIdx>"`, and the resolver
  accepts that suffix (the index is verified against the name before use).
"""
struct ContingencyCase
  name::String
  kind::Symbol
  element::String
  function ContingencyCase(name::String, kind::Symbol, element::String)
    kind === :branch || throw(ArgumentError("ContingencyCase: only kind = :branch is supported in this version, got :$(kind)."))
    return new(name, kind, element)
  end
end

ContingencyCase(element::String) = ContingencyCase(element, :branch, element)

"""
    ContingencyResult

Outcome of one [`ContingencyCase`](@ref) evaluated by
[`runContingencies!`](@ref).

# Fields
- `name::String`: the case name.
- `converged::Bool`: whether the post-outage power flow converged.
- `iterations::Int`: Newton iterations of the contingency solve (0 when the
  solve never ran).
- `max_vm_pu::Float64` / `min_vm_pu::Float64`: voltage envelope over the
  non-isolated buses of the solved case (NaN when not solved).
- `max_branch_loading_pct::Float64`: maximum loading over all rated
  branches (|S| at either end against `sn_MVA`); NaN when the case did not
  solve or no branch carries a finite rating.
- `overloaded_branches::Vector{String}`: branch names with loading above
  100 percent.
- `voltage_violations::Vector{String}`: bus names outside the
  `[vm_min_pu, vm_max_pu]` band.
- `island_count::Int`: AC islands of the post-outage topology (0 when the
  net could not be evaluated).
- `error::Union{Nothing,String}`: `nothing` on success; otherwise the
  failure in one line ("islanded without reference", the solver status, or
  the exception message). Failures are REPORTED, never thrown.
"""
struct ContingencyResult
  name::String
  converged::Bool
  iterations::Int
  max_vm_pu::Float64
  min_vm_pu::Float64
  max_branch_loading_pct::Float64
  overloaded_branches::Vector{String}
  voltage_violations::Vector{String}
  island_count::Int
  error::Union{Nothing,String}
end

"""
    generateN1Branches(net::Net; include_transformers::Bool = true) -> Vector{ContingencyCase}

Generate one branch-outage [`ContingencyCase`](@ref) per in-service branch
of `net` (aggregate `status == 1`; partially open branches are skipped,
their outage is already half-effective). With
`include_transformers = false` only AC-line branches are listed; a branch
counts as a transformer when its component type is `Trafo` (net-builder
path) or its winding ratio is nonzero (the classical MATPOWER indicator:
line rows carry `ratio = 0`).
"""
function generateN1Branches(net::Net; include_transformers::Bool = true)::Vector{ContingencyCase}
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
    name = getCompName(br.comp)
    # parallel circuits share a component name; disambiguate by branch index
    # so every circuit gets its OWN outage case (a bare name would always
    # resolve to the first circuit and evaluate the same outage twice)
    element = name_count[name] > 1 ? string(name, "#", br.branchIdx) : name
    push!(cases, ContingencyCase(element, :branch, element))
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

"""
    generateContingenciesFromFOR001(net::Net) -> Vector{ContingencyCase}

Map the FOR001 contingency branch names imported from MATPOWER metadata
(`mpc.for001_contingencies`, stored on `net.for001Contingencies`) to
[`ContingencyCase`](@ref)s. Names that do not resolve to a branch are still
listed; [`runContingencies!`](@ref) reports them as failed cases with an
actionable error instead of dropping them silently.
"""
generateContingenciesFromFOR001(net::Net)::Vector{ContingencyCase} = [ContingencyCase(name, :branch, name) for name in net.for001Contingencies]

# evaluate the solved contingency net against the limit band; pure reader
function _contingency_metrics(cnet::Net, vm_min_pu::Float64, vm_max_pu::Float64)
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
  overloaded = String[]
  for br in cnet.branchVec
    br.status == 1 || continue
    rating = br.sn_MVA
    (rating === nothing || !isfinite(rating) || rating <= 0.0) && continue
    s_from = br.fBranchFlow === nothing ? 0.0 : hypot(something(br.fBranchFlow.pFlow, 0.0), something(br.fBranchFlow.qFlow, 0.0))
    s_to = br.tBranchFlow === nothing ? 0.0 : hypot(something(br.tBranchFlow.pFlow, 0.0), something(br.tBranchFlow.qFlow, 0.0))
    loading = 100.0 * max(s_from, s_to) / rating
    (isnan(max_loading) || loading > max_loading) && (max_loading = loading)
    loading > 100.0 && push!(overloaded, getCompName(br.comp))
  end
  return (; vmin = isfinite(vmin) ? vmin : NaN, vmax = isfinite(vmax) ? vmax : NaN, violations, max_loading, overloaded)
end

# one worker evaluation: template copy -> outage -> solve -> metrics.
# Pure per-case function: touches no shared container, all mutation happens
# on the case-local deepcopy.
function _run_one_contingency(template::Net, case::ContingencyCase, vm_min_pu::Float64, vm_max_pu::Float64, maxIte::Int, tol::Float64, pf_kwargs)
  cnet = deepcopy(template)
  idx = _resolve_contingency_branch(cnet, case.element)
  if idx === nothing
    return ContingencyResult(case.name, false, 0, NaN, NaN, NaN, String[], String[], 0, "unknown branch $(case.element) (not found in net.branchVec)")
  end
  removeBranch!(net = cnet, branchNr = idx)
  markIsolatedBuses!(net = cnet, log = false)
  island_count = length(detect_ac_islands(cnet).rows)
  it = 0
  try
    it, erg = runpf!(cnet, maxIte, tol, 0; islands_enabled = true, pf_kwargs...)
    if erg != 0
      return ContingencyResult(case.name, false, it, NaN, NaN, NaN, String[], String[], island_count, "power flow did not converge (status $(erg))")
    end
    calcNetLosses!(cnet)
    m = _contingency_metrics(cnet, vm_min_pu, vm_max_pu)
    return ContingencyResult(case.name, true, it, m.vmax, m.vmin, m.max_loading, m.overloaded, m.violations, island_count, nothing)
  catch err
    err isa InterruptException && rethrow(err)
    msg = sprint(showerror, err)
    # a removal that splits off a reference-less island surfaces as the
    # island-validation error; report it under the stable short label
    if occursin("reference", msg) && (occursin("island", msg) || occursin("Island", msg))
      return ContingencyResult(case.name, false, it, NaN, NaN, NaN, String[], String[], island_count, "islanded without reference")
    end
    return ContingencyResult(case.name, false, it, NaN, NaN, NaN, String[], String[], island_count, first(split(msg, '\n')))
  end
end

"""
    runContingencies!(net::Net, cases::Vector{ContingencyCase};
                      vm_min_pu = 0.9, vm_max_pu = 1.1,
                      maxIte = 30, tol = 1e-8,
                      parallel_enabled = nothing, parallel_max_tasks = nothing,
                      parallel_min_work_items = nothing,
                      kwargs...) -> Vector{ContingencyResult}

Evaluate branch-outage contingencies. The base `net` is NEVER mutated: a
solved template copy is created first (warm start: every contingency solve
starts from the base-case voltages; when the base case itself does not
converge, the contingencies start flat and a warning is printed once), and
each case works on its own `deepcopy` of that template, removes its branch
via [`removeBranch!`](@ref), solves, and evaluates the
`[vm_min_pu, vm_max_pu]` band plus branch loadings against `sn_MVA`
ratings. Cases are returned in input order; failures (no convergence,
islanding without reference, unknown element) are REPORTED in the result,
never thrown.

Remaining `kwargs...` are forwarded to the contingency `runpf!` solves.
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
  parallel_enabled::Union{Nothing,Bool} = nothing,
  parallel_max_tasks::Union{Nothing,Int} = nothing,
  parallel_min_work_items::Union{Nothing,Int} = nothing,
  kwargs...,
)::Vector{ContingencyResult}
  vm_min_pu < vm_max_pu || throw(ArgumentError("runContingencies!: vm_min_pu must be below vm_max_pu."))
  isempty(cases) && return ContingencyResult[]

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
    @warn "runContingencies!: the base case did not converge — contingencies start FLAT instead of warm (decision recorded in the multi-core analysis)."
    template = deepcopy(net)
    template.flatstart = true
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
          results[idx] = _run_one_contingency(template, cases[idx], vm_min_pu, vm_max_pu, maxIte, tol, pf_kwargs)
        end
      end for ci in eachindex(chunks)
    ]
    foreach(wait, tasks)
  else
    for idx in eachindex(cases)
      results[idx] = _run_one_contingency(template, cases[idx], vm_min_pu, vm_max_pu, maxIte, tol, pf_kwargs)
    end
  end
  return ContingencyResult[results[i] for i in eachindex(cases)]
end

"""
    printContingencyResults([io::IO], results::Vector{ContingencyResult}; max_rows = 50)

Print the contingency results as a fixed-width table (style of
`printShortCircuitResult`): case, convergence, iterations, voltage
envelope, worst loading, violation counts, and the error line for failed
cases. Rows beyond `max_rows` are summarized in one line.
"""
printContingencyResults(results::Vector{ContingencyResult}; max_rows::Int = 50) = printContingencyResults(stdout, results; max_rows = max_rows)

function printContingencyResults(io::IO, results::Vector{ContingencyResult}; max_rows::Int = 50)
  println(io, "N-1 contingency results (", length(results), " case(s), ", count(r -> r.converged, results), " converged)")
  println(io, "-"^116)
  @printf(io, "%-28s %-9s %5s %9s %9s %12s %8s %8s %7s  %s\n", "case", "converged", "iter", "Vmin[pu]", "Vmax[pu]", "loading[%]", "overld", "V-viol", "islands", "error")
  println(io, "-"^116)
  shown = 0
  for r in results
    shown >= max_rows && break
    shown += 1
    @printf(
      io,
      "%-28s %-9s %5d %9s %9s %12s %8d %8d %7d  %s\n",
      _fitColumn(r.name, 28),
      r.converged ? "yes" : "NO",
      r.iterations,
      isnan(r.min_vm_pu) ? "-" : @sprintf("%.4f", r.min_vm_pu),
      isnan(r.max_vm_pu) ? "-" : @sprintf("%.4f", r.max_vm_pu),
      isnan(r.max_branch_loading_pct) ? "-" : @sprintf("%.1f", r.max_branch_loading_pct),
      length(r.overloaded_branches),
      length(r.voltage_violations),
      r.island_count,
      r.error === nothing ? "" : r.error,
    )
  end
  shown < length(results) && println(io, "... ", length(results) - shown, " more row(s) not shown (max_rows = ", max_rows, ")")
  println(io, "-"^116)
  return nothing
end

"""
    writeContingencyResultsCSV(path::AbstractString, results::Vector{ContingencyResult})

Write the contingency results as a semicolon-separated CSV (one row per
case, input order): name, converged, iterations, voltage envelope, worst
loading, overloaded branch list, voltage-violation list, island count, and
the error text. Returns `path`.
"""
function writeContingencyResultsCSV(path::AbstractString, results::Vector{ContingencyResult})
  open(path, "w") do io
    println(io, "name;converged;iterations;min_vm_pu;max_vm_pu;max_branch_loading_pct;overloaded_branches;voltage_violations;island_count;error")
    for r in results
      println(io, join([r.name, r.converged, r.iterations, r.min_vm_pu, r.max_vm_pu, r.max_branch_loading_pct, join(r.overloaded_branches, ","), join(r.voltage_violations, ","), r.island_count, r.error === nothing ? "" : r.error], ";"))
    end
  end
  return path
end
