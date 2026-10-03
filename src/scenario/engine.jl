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

# file: src/scenario/engine.jl
# purpose: the scenario engine: it evaluates patch scenarios and N-1 outages on REUSED
#          per-chunk working copies instead of one deepcopy per case. The
#          full run is the existing contingency logic (ladder, metrics)
#          verbatim; runContingencies! routes through this engine, and the
#          case118 CSV fixture pins the results byte for byte to the
#          pre-engine implementation.

# One resolved work item of a batch. Resolution happens ONCE on the
# unmutated template (the same indices the per-case resolution of the old
# implementation produced), so the workers never run name lookups.
struct _ScenarioOutageItem
  name::String
  weight::Float64
  kind::Symbol        # :branch or :gen
  internal::Int       # branchVec / prosumpsVec position on the template
end

# a case whose element did not resolve: reported, never thrown (old contract)
struct _ScenarioFailedItem
  name::String
  weight::Float64
  message::String
end

# a general patch scenario (anything that is not a single status-0 outage)
struct _ScenarioPatchItem
  scenario::Scenario
end

const _ScenarioItem = Union{_ScenarioOutageItem,_ScenarioFailedItem,_ScenarioPatchItem}

"""
    ScenarioWorker

One reusable working copy of the engine template. After every evaluation
the worker is reset to the template state (patch undo plus solver-state
reset), so the next scenario starts from exactly the state a fresh
deepcopy would have provided.
"""
struct ScenarioWorker
  net::Net
end

"""
    ScenarioEngine

The scenario batch engine: the solved base template
(never mutated after construction), the base branch loadings every result
reports its deltas against, the resolved solve parameters, the optional
[`ScenarioIndex`](@ref) that patch scenarios address components through,
and, with `screening_mode` other than `:off`, the screening state:
the base Jacobian factorization plus everything one Woodbury-corrected
Newton step needs. A distributed-slack batch screens on the AUGMENTED
system: lambda as the extra state, a generator outage zeroes
the lost unit's participation factor and renormalizes the rest (a rank-1
column update the generic correction covers), and a participating unit
whose estimated output leaves its P band flags the scenario. `screen ===
nothing` with screening requested means the state could not be built
(base not converged, residual inconsistent, voltage-dependent
injections, no surviving participant); every scenario then gets the full
run.
"""
struct ScenarioEngine
  template::Net
  base_converged::Bool
  base_loadings::Dict{String,Float64}
  index::Union{Nothing,ScenarioIndex}
  vm_min_pu::Float64
  vm_max_pu::Float64
  maxIte::Int
  tol::Float64
  ladder::Vector{Symbol}
  pf_kwargs::Any
  screening_mode::Symbol
  screening_margin_pct::Float64
  screen::Any
  # the buses the base case solves on as references (slack buses and the
  # bus an island without a slack was given): a case names a reference as
  # taken over only when it is none of these
  base_reference_buses::Set{Int}
  # the cold check of the warm active set (0.30.3): when switched on, a warm
  # result within this margin above vm_min_pu, or with a violation, is solved
  # cold as well
  warm_cold_check::Bool
  warm_cold_check_margin_pu::Float64
end

"""
    ScenarioEngine(net; vm_min_pu, vm_max_pu, maxIte, tol, ladder, index, pf_kwargs)

Build the engine template from `net`: solve the base case (with the full
solver rescue ladder as fallback, then a flat template with a warning,
exactly the old `runContingencies!` preamble), clear the solver status
and Q-limit logs the workers must not inherit, and record the base branch
loadings. `net` itself is never mutated.
"""
function ScenarioEngine(net::Net; vm_min_pu::Float64 = 0.9, vm_max_pu::Float64 = 1.1, maxIte::Int = 30, tol::Float64 = 1e-8, ladder::Vector{Symbol} = Symbol[:warm], index::Union{Nothing,ScenarioIndex} = nothing, pf_kwargs = NamedTuple(), screening_mode::Symbol = :off, screening_margin_pct::Float64 = 10.0, warm_active_set::Bool = true, warm_note::Union{Nothing,Base.RefValue{String}} = nothing, warm_cold_check::Bool = true, warm_cold_check_margin_pu::Float64 = DEFAULT_WARM_COLD_CHECK_MARGIN_PU)
  screening_mode in CONTINGENCY_SCREENING_MODE_VALUES || throw(ArgumentError("ScenarioEngine: screening_mode must be one of $(CONTINGENCY_SCREENING_MODE_VALUES), got :$(screening_mode)."))
  (isfinite(screening_margin_pct) && screening_margin_pct >= 0.0) || throw(ArgumentError("ScenarioEngine: screening_margin_pct must be a finite value >= 0."))
  vm_min_pu < vm_max_pu || throw(ArgumentError("ScenarioEngine: vm_min_pu must be below vm_max_pu."))
  (isfinite(warm_cold_check_margin_pu) && warm_cold_check_margin_pu >= 0.0) || throw(ArgumentError("ScenarioEngine: warm_cold_check_margin_pu must be a finite value >= 0."))
  template = deepcopy(net)
  base_converged = try
    _, base_erg = runpf!(template, maxIte, tol, 0; islands_enabled = true, pf_kwargs...)
    base_erg == 0
  catch err
    _rethrow_unless_solver_failure(err)
    false
  end
  if !base_converged
    # #331 Phase 1 item 4: the base case gets the full solver rescue ladder
    # (alternate start, autodamp, settled Q-limits, DC seed) before being
    # declared unsolvable; only then fall back to the flat template.
    # The ladder runs on the SAME keywords as the base solve above, each
    # strategy only overriding what it names: a default PowerFlowConfig
    # here dropped distributed slack, auto_slack, the Q-limit settings,
    # newton_update and power_mode, and a rescued template then carried a
    # different slack distribution than the batch asked for (#456). The
    # plain attempt is not repeated, it is the base solve that just failed.
    base_kw = (; islands_enabled = true, pf_kwargs...)
    base_flatstart = get(base_kw, :opt_flatstart, net.flatstart) === true
    for (name, variant_kw) in _rescue_kwargs_variants(base_kw, base_flatstart)
      # a fresh copy per strategy: a failed attempt leaves switched bus
      # types and a promoted slack behind, not only voltages
      rescued = deepcopy(net)
      println("rescue: AC solve did not converge, retrying with strategy '", name, "'.")
      rescued_ok = try
        r_it, r_erg = runpf!(rescued, maxIte, tol, 0; variant_kw...)
        r_erg == 0 && println("rescue: strategy '", name, "' converged after ", r_it, " iteration(s).")
        r_erg == 0
      catch err
        _rethrow_unless_solver_failure(err)
        false
      end
      if rescued_ok
        base_converged = true
        template = rescued
        break
      end
    end
    base_converged || println("rescue: no strategy converged.")
    if !base_converged
      @warn "ScenarioEngine: the base case did not converge even with the solver rescue ladder; scenarios start FLAT instead of warm."
      template = deepcopy(net)
      template.flatstart = true
    end
  end
  # Warm active set (0.30.2, `contingency.warm_active_set`): the base
  # case's clamped machines (bus name => reached side), read before the
  # hygiene below empties the logs, start every worker solve as PQ at that
  # limit (`qlimit_warm_clamps`); the release rule stays free to undo them.
  # Without a converged base case or with Q limits off there is nothing to
  # carry over, and one line says so.
  if warm_active_set
    qlimits_on = get(pf_kwargs, :qlimits_enabled, true) === true
    note = if !base_converged
      "Warm active set: not applied, the base case did not converge."
    elseif !qlimits_on
      "Warm active set: not applied, Q limits are off."
    else
      bus_names = Dict{Int,String}(idx => name for (name, idx) in template.busDict)
      warm = Dict{String,Symbol}()
      for (bus, side) in template.qLimitEvents
        side in (:min, :max) || continue
        name = get(bus_names, bus, nothing)
        name === nothing || (warm[name] = side)
      end
      pf_kwargs = (; pf_kwargs..., qlimit_warm_clamps = warm)
      "Warm active set: $(length(warm)) clamped machine bus(es) of the base case carried into every outage."
    end
    # the service passes a Ref and writes the line into run.log; a library
    # call without one sees it on the console
    warm_note === nothing ? println(note) : (warm_note[] = note)
  end
  # template hygiene: the workers need the solved VOLTAGES, not the base
  # solver status or its Q-limit event logs
  template._rectangular_pf_status = nothing
  template._dc_pf_status = nothing
  empty!(template.qLimitLog)
  empty!(template.qLimitEvents)
  empty!(template.qLimitInitialPVRows)
  base_loadings = Dict{String,Float64}()
  if base_converged
    calcNetLosses!(template)
    base_loadings = _base_branch_loadings(template)
  end
  screen = screening_mode === :off ? nothing : _build_screening_state(template, base_converged, tol, pf_kwargs, base_loadings)
  base_refs = Set{Int}(Int(ps.comp.cFrom_bus) for ps in template.prosumpsVec if isSlack(ps))
  # islands always promote a generating unit (the solver does the same)
  for row in detect_ac_islands(template; promote_generators = true).rows
    row.chosen_ref_bus > 0 && push!(base_refs, Int(row.chosen_ref_bus))
  end
  return ScenarioEngine(template, base_converged, base_loadings, index, vm_min_pu, vm_max_pu, maxIte, tol, ladder, pf_kwargs, screening_mode, screening_margin_pct, screen, base_refs, warm_cold_check, warm_cold_check_margin_pu)
end

# --- worker reset -------------------------------------------------------------

# Copy every scalar-valued field (numbers, enums, symbols, strings, and
# `nothing`) from src to dst. Component references, vectors, and dicts are
# deliberately left alone: reassigning them would alias template-owned
# mutable objects into the worker.
function _scenario_copy_scalar_fields!(dst::T, src::T) where {T}
  for f in fieldnames(T)
    # the power-mode cache belongs to the worker, never to the template:
    # a template without one (`nothing`, e.g. a rescued or flat base case)
    # would otherwise wipe the worker's cache after every scenario (#457)
    f === :_power_cache && continue
    v = getfield(src, f)
    if v === nothing || v isa Number || v isa Symbol || v isa AbstractString || v isa Enum
      setfield!(dst, f, v)
    end
  end
  return nothing
end

"""
    _reset_scenario_worker!(work, template) -> Nothing

Reset a working copy to the template state after an evaluation. The patch
itself is already undone (undo log replay, prosumer reinsertion); what
remains is everything a SOLVE mutates: node voltages and types, branch
status/flow/loss fields, prosumer write-backs, the isolated-node registry,
slack rows, Q-limit tables and logs, loss accumulators, and the solver
status. The step-3 acceptance test proves the reset with the SCF
reflection comparison: after evaluations with solves, the worker equals
the template bitwise on every compared field. Transformer tap CONTROLLER
state is not reset here: scenario nets with active tap controllers are
rejected at validation time (step 2), so no engine path moves a tap.
"""
function _reset_scenario_worker!(work::Net, template::Net)
  _scenario_copy_scalar_fields!(work, template)
  for i in eachindex(template.nodeVec)
    _scenario_copy_scalar_fields!(work.nodeVec[i], template.nodeVec[i])
  end
  for i in eachindex(template.branchVec)
    wb = work.branchVec[i]
    tb = template.branchVec[i]
    _scenario_copy_scalar_fields!(wb, tb)
    # flow objects are mutable; fresh copies keep the template unaliased
    wb.fBranchFlow = tb.fBranchFlow === nothing ? nothing : deepcopy(tb.fBranchFlow)
    wb.tBranchFlow = tb.tBranchFlow === nothing ? nothing : deepcopy(tb.tBranchFlow)
  end
  for i in eachindex(template.prosumpsVec)
    _scenario_copy_scalar_fields!(work.prosumpsVec[i], template.prosumpsVec[i])
  end
  for i in eachindex(template.shuntVec)
    _scenario_copy_scalar_fields!(work.shuntVec[i], template.shuntVec[i])
  end
  for i in eachindex(template.linkVec)
    _scenario_copy_scalar_fields!(work.linkVec[i], template.linkVec[i])
  end
  for (wv, tv) in ((work.slackVec, template.slackVec), (work.isoNodes, template.isoNodes), (work.qmin_pu, template.qmin_pu), (work.qmax_pu, template.qmax_pu))
    empty!(wv)
    append!(wv, tv)
  end
  for (wv, tv) in ((work.totalLosses, template.totalLosses), (work.totalBusPower, template.totalBusPower), (work.qLimitLog, template.qLimitLog), (work.qLimitInitialPVRows, template.qLimitInitialPVRows))
    empty!(wv)
    append!(wv, tv)
  end
  empty!(work.qLimitEvents)
  merge!(work.qLimitEvents, template.qLimitEvents)
  work.control_result = template.control_result
  work._rectangular_pf_status = template._rectangular_pf_status
  work._dc_pf_status = template._dc_pf_status
  return nothing
end

# --- screening ----------------------------------------------------------------

# Everything ONE Woodbury-corrected Newton step needs, captured once from the
# solved template: the expanded Y-bus (bus index == matrix index, isolated
# buses as zero rows), the solved state, the residual layout, the factorized
# base Jacobian, and the bridge set (outages that split the graph are never
# screened). Immutable and shared read-only across worker threads.
struct ScreeningState
  Ybus::SparseMatrixCSC{ComplexF64,Int}
  V0::Vector{ComplexF64}
  S0::Vector{ComplexF64}
  bus_types::Vector{Symbol}
  Vset::Vector{Float64}
  slack_idx::Int
  non_slack::Vector{Int}
  J0::SparseMatrixCSC{Float64,Int}
  lu0::Any
  bridges::Set{Int}
  base_island_count::Int
  iso::Set{Int}
  # distributed slack: the base participation state with
  # lambda = 0 (the solved template dispatch already carries the base
  # correction, so the augmented residual vanishes at zero). Shared
  # read-only; every scenario works on its own copy. `nothing` for
  # classical batches.
  ds::Any
  # base-case loading in percent per branch index (NaN without a rating):
  # the reference of the change-based loading flag
  base_loading::Vector{Float64}
end

# Bridges of the closed-branch graph via iterative Tarjan low-link. Closed
# bus links count as connectivity (they cannot be outaged here but they keep
# a doubly-connected pair doubly connected); parallel circuits are handled
# by edge-id skipping, so a twin circuit's partner is never a bridge.
function _scenario_bridge_branches(net::Net)::Set{Int}
  n = length(net.nodeVec)
  adj = [Tuple{Int,Int}[] for _ in 1:n]
  next_link_id = -1
  for (bi, br) in enumerate(net.branchVec)
    (br.status == 1 && _branch_terminal_state(br) === :closed) || continue
    f = Int(br.fromBus)
    t = Int(br.toBus)
    f == t && continue
    push!(adj[f], (t, bi))
    push!(adj[t], (f, bi))
  end
  for lk in net.linkVec
    lk.status == 1 || continue
    f = Int(lk.fromBus)
    t = Int(lk.toBus)
    (1 <= f <= n && 1 <= t <= n && f != t) || continue
    push!(adj[f], (t, next_link_id))
    push!(adj[t], (f, next_link_id))
    next_link_id -= 1
  end
  disc = zeros(Int, n)
  low = zeros(Int, n)
  bridges = Set{Int}()
  timer = 0
  # iterative DFS: (vertex, incoming edge id, next adjacency position)
  stack = Tuple{Int,Int,Int}[]
  for start in 1:n
    disc[start] != 0 && continue
    isempty(adj[start]) && continue
    timer += 1
    disc[start] = low[start] = timer
    push!(stack, (start, 0, 1))
    while !isempty(stack)
      v, in_edge, pos = stack[end]
      if pos <= length(adj[v])
        stack[end] = (v, in_edge, pos + 1)
        w, eid = adj[v][pos]
        eid == in_edge && continue
        if disc[w] == 0
          timer += 1
          disc[w] = low[w] = timer
          push!(stack, (w, eid, 1))
        else
          low[v] = min(low[v], disc[w])
        end
      else
        pop!(stack)
        if !isempty(stack)
          p = stack[end][1]
          low[p] = min(low[p], low[v])
          low[v] > disc[p] && in_edge > 0 && push!(bridges, in_edge)
        end
      end
    end
  end
  return bridges
end

# Build the screening state, or return nothing when screening cannot be
# trusted on this template: base not converged, voltage-dependent
# injections (the specified-S vector is state dependent), a residual that
# is not actually converged at the captured state (Q-limit or controller
# write-back left types and schedule inconsistent), or a singular base
# Jacobian. A distributed-slack batch builds the AUGMENTED
# system: the participation state at lambda = 0 (the solved template
# dispatch already carries the base correction), one extra state and
# residual row, same factorization machinery.
function _build_screening_state(template::Net, base_converged::Bool, tol::Float64, pf_kwargs, base_loadings::AbstractDict = Dict{String,Float64}())
  base_converged || return nothing
  has_voltage_dependent_control(template) && return nothing
  n = length(template.nodeVec)
  Yred = createYBUS(net = template, sparse = true)
  Ybus = size(Yred, 1) == n ? Yred : _expand_ybus_for_isolated_nodes(Yred, n, template.isoNodes)
  V0, slack_idx = initialVrect(template; flatstart = false)
  # The EFFECTIVE specification at the solved state, not the schedule: a
  # Q-limit clamp demotes a PV bus to PQ with Q held at the limit, so the
  # scheduled Q of buildComplexSVec would leave a spurious base residual
  # there. The calculated injections ARE what the solver converged on
  # (identical to the schedule wherever no clamp acted), and the base
  # residual then vanishes by construction on every P/Q row.
  S0 = calc_injections(Ybus, V0)
  Vset = _bus_voltage_setpoints_from_prosumers(template)
  bus_types = Vector{Symbol}(undef, n)
  for (k, node) in enumerate(template.nodeVec)
    bt = getNodeType(node)
    if bt == Slack
      bus_types[k] = :Slack
    elseif bt == PV
      bus_types[k] = :PV
    elseif bt == PQ
      bus_types[k] = :PQ
    elseif bt == Isolated
      bus_types[k] = :PQ
      S0[k] = 0.0 + 0.0im
    else
      return nothing
    end
  end
  # distributed slack (4b): rebuild the participation state on the solved
  # template with the SAME options the batch hands the solver; a build
  # failure or a classical fallback leaves ds = nothing (classical screen
  # for ref_only fallback, full runs when the batch still solves with
  # dslack because the dimensions would not match)
  kw = Dict(pairs(pf_kwargs))
  ds = nothing
  if get(kw, :distributed_slack_enabled, false) == true
    ds = try
      build_distributed_slack_state(template, bus_types; p_mode = get(kw, :distributed_slack_p_mode, :pg_weighted), fallback = get(kw, :distributed_slack_fallback, :error), weights = get(kw, :distributed_slack_weights, Dict{String,Float64}()), respect_p_limits = get(kw, :distributed_slack_respect_p_limits, true), island_label = template.name)
    catch err
      _rethrow_unless_solver_failure(err)
      return nothing
    end
    ds === nothing && return nothing
    ds.lambda = 0.0
    ds.lambda_trial = 0.0
    ds.delta_lambda = 0.0
  end
  F0 = try
    mismatch_rectangular(Ybus, V0, S0, bus_types, Vset, slack_idx; dslack = ds)
  catch err
    _rethrow_unless_solver_failure(err)
    return nothing
  end
  # the captured state must BE the solved state; a large residual means the
  # capture is inconsistent with what the solver converged on
  maximum(abs, F0) <= max(1.0e-6, 100.0 * tol) || return nothing
  non_slack = non_slack_indices(n, slack_idx)
  iso_set = Set(template.isoNodes)
  J0 = try
    _screen_freeze_isolated!(build_rectangular_jacobian_pq_pv(Ybus, V0, bus_types, Vset, slack_idx; dslack = ds), non_slack, iso_set)
  catch err
    _rethrow_unless_solver_failure(err)
    return nothing
  end
  # a singular base Jacobian is the expected failure here: no screening
  lu0 = try
    lu(J0)
  catch err
    _rethrow_unless_solver_failure(err)
    return nothing
  end
  base_loading = Float64[get(base_loadings, getCompName(br.comp), NaN) for br in template.branchVec]
  return ScreeningState(Ybus, V0, S0, bus_types, Vset, slack_idx, non_slack, J0, lu0, _scenario_bridge_branches(template), length(detect_ac_islands(template).rows), iso_set, ds, base_loading)
end

# Per-scenario copy of the shared participation state (lambda_trial is
# mutated during the residual evaluations, and the shared state is read
# concurrently by other chunks). For a generator outage the lost unit's
# factor goes to zero and the surviving factors are renormalized;
# `nothing` when the lost unit was the LAST participant (the augmented
# system loses its lambda state, dimensions no longer match the base
# factorization, the full run owns that).
function _screen_dslack_for_outage(ds, kind::Symbol, internal::Int)
  ds === nothing && return (nothing, true)
  table = ds.gen_table
  if kind === :gen
    surviving = [row for row in table if row.gen_index != internal]
    if length(surviving) != length(table)
      isempty(surviving) && return (nothing, false)
      wsum = sum(row.weight for row in surviving)
      (isfinite(wsum) && wsum > 0.0) || return (nothing, false)
      n_alpha = length(ds.alpha)
      alpha = zeros(Float64, n_alpha)
      table2 = similar(surviving, 0)
      for row in surviving
        a = row.weight / wsum
        push!(table2, (bus = row.bus, gen_index = row.gen_index, weight = row.weight, alpha = a))
        alpha[row.bus] += a
      end
      return (DistributedSlackState(alpha, 0.0, 0.0, 0.0, ds.mode, ds.respect_p_limits, table2, ds.dropped), true)
    end
  end
  return (DistributedSlackState(copy(ds.alpha), 0.0, 0.0, 0.0, ds.mode, ds.respect_p_limits, copy(table), ds.dropped), true)
end

# Freeze the isolated-bus rows of a screening Jacobian: an isolated bus
# carries a zero Y-bus row, zero injections and vm = va = 0, so its P/Q
# residual rows are identically zero and the raw Jacobian is SINGULAR
# (found on RealGrid, 198 de-energized buses). The solver escapes through
# the island decomposition; the full-dimension screening system instead
# pins delta Vr = delta Vi = 0 with two unit entries per isolated bus
# (their rows hold no stored entries, so this is a pure insertion, and
# J1 - J0 cancels the entries again).
function _screen_freeze_isolated!(J::SparseMatrixCSC{Float64,Int}, non_slack::Vector{Int}, iso::Set{Int})
  isempty(iso) && return J
  nn = length(non_slack)
  for (k, bus) in enumerate(non_slack)
    bus in iso || continue
    J[2 * k - 1, k] = 1.0
    J[2 * k, nn + k] = 1.0
  end
  return J
end

# Solve (J0 + ΔJ) x = b through the factorized base: ΔJ = J1 - J0 is confined
# to the rows and columns an outage touches, so the correction is a small
# Woodbury block (column-side U Sel' or row-side Sel W', whichever is
# smaller). Returns nothing when the update is not low-rank or the capacitance
# matrix is singular; the caller then falls back to the full run.
function _screen_woodbury_solve(lu0, J0::SparseMatrixCSC{Float64,Int}, J1::SparseMatrixCSC{Float64,Int}, b::Vector{Float64})
  dJ = J1 - J0
  nz = nnz(dJ)
  nz == 0 && return lu0 \ b
  rows_nz, cols_nz, _ = findnz(dJ)
  rows = sort!(unique(rows_nz))
  cols = sort!(unique(cols_nz))
  k = min(length(rows), length(cols))
  k > 32 && return nothing
  x0 = lu0 \ b
  local corr
  try
    if length(cols) <= length(rows)
      U = Matrix(dJ[:, cols])
      Z = lu0 \ U
      M = Matrix{Float64}(LinearAlgebra.I, length(cols), length(cols)) + Z[cols, :]
      corr = Z * (M \ x0[cols])
    else
      m = size(J0, 1)
      U = zeros(Float64, m, length(rows))
      for (i, r) in enumerate(rows)
        U[r, i] = 1.0
      end
      W = Matrix(dJ[rows, :])
      Z = lu0 \ U
      M = Matrix{Float64}(LinearAlgebra.I, length(rows), length(rows)) + W * Z
      corr = Z * (M \ (W * x0))
    end
  catch err
    _rethrow_unless_solver_failure(err)
    return nothing
  end
  return x0 - corr
end

# An estimate is only trusted when the one-step solution nearly satisfies
# the outage equations: above this max-abs residual (pu power / pu voltage
# rows) the scenario is flagged for the full run regardless of its margins.
# History: 0.02 came from the case118 acceptance (a 0.047 pu residual hid a
# 0.05 pu voltage-drop underestimate); case300 then forced 0.005, its
# outage 57-63 collapses buses 63/64/526 to 0.84 pu while the one-step
# residual is only 0.0104 (a flat mismatch near the local collapse), so
# the gate holds a factor of two below that.
const _SCREENING_TRUST_MISMATCH_PU = 0.005

# the answer of _screen_outage for a scenario that CANNOT be estimated:
# the reason names the category (the calibration report breaks the
# screening share's denominator down by exactly these)
_screen_unscreenable(reason::Symbol) = (flagged = true, estimate = nothing, est_overloads = OverloadRecord[], est_violations = String[], reason = reason)

# Screen one outage: one Woodbury-corrected Newton step from the base state
# under the outaged model, then estimated loadings and voltages from the
# trial state. Always returns a NamedTuple; `estimate === nothing` with a
# `reason` when the scenario is not screenable (:no_screening_state,
# :branch_not_closed, :bridge, :isolated_endpoint, :slack_unit,
# :invalid_bus, :isolated_bus, :last_participant, :numeric), otherwise the
# flag decision plus the estimate (`reason` is :structural_pv_gen when the
# structural flag fired, :ok otherwise).
function _screen_outage(engine::ScenarioEngine, it::_ScenarioOutageItem)
  sc = engine.screen
  sc === nothing && return _screen_unscreenable(:no_screening_state)
  template = engine.template
  outaged_branch = 0
  local Yb, S1, types1
  if it.kind === :branch
    br = template.branchVec[it.internal]
    (br.status == 1 && _branch_terminal_state(br) === :closed) || return _screen_unscreenable(:branch_not_closed)
    it.internal in sc.bridges && return _screen_unscreenable(:bridge)
    f = Int(br.fromBus)
    t = Int(br.toBus)
    (f in sc.iso || t in sc.iso) && return _screen_unscreenable(:isolated_endpoint)
    y11, y12, y21, y22 = calcAdmittance(br, br.comp.cVN, template.baseMVA)
    Yb = copy(sc.Ybus)
    Yb[f, f] -= y11
    Yb[t, t] -= y22
    Yb[f, t] -= y12
    Yb[t, f] -= y21
    S1 = sc.S0
    types1 = sc.bus_types
    outaged_branch = it.internal
  else
    ps = template.prosumpsVec[it.internal]
    # losing a slack unit changes the reference structure; the full run owns that
    isSlack(ps) && return _screen_unscreenable(:slack_unit)
    bus = getPosumerBusIndex(ps)
    (1 <= bus <= length(template.nodeVec)) || return _screen_unscreenable(:invalid_bus)
    bus in sc.iso && return _screen_unscreenable(:isolated_bus)
    Yb = sc.Ybus
    S1 = copy(sc.S0)
    S1[bus] -= ComplexF64(something(ps.pVal, 0.0), something(ps.qVal, 0.0)) / template.baseMVA
    types1 = sc.bus_types
    if sc.bus_types[bus] === :PV
      survivor = any(other !== ps && isGenerator(other) && isRegulating(other) && getPosumerBusIndex(other) == bus for other in template.prosumpsVec)
      if !survivor
        types1 = copy(sc.bus_types)
        types1[bus] = :PQ
      end
    end
  end
  ds1, screenable = _screen_dslack_for_outage(sc.ds, it.kind, it.internal)
  screenable || return _screen_unscreenable(:last_participant)
  # Structural flag (case300 Gen_115#205 finding): the outage of a
  # REGULATING unit shrinks its bus's aggregated Q capability, which is
  # invisible to a residual step whenever the unit's schedule is (near)
  # zero, yet the full solve may clamp at the smaller Q limit and demote
  # the bus (vmin collapsed by 0.062 pu behind a 4e-13 pu residual there,
  # at a bus the BASE solve had already clamped from PV to PQ, so the bus
  # type alone is not the criterion). Any bus that carries regulating
  # generation, whether currently PV or Q-clamped PQ, is in that class; a
  # modeled PV-to-PQ demotion (last regulator lost, types1 changed) shows
  # in the residual and needs no flag. No m2 gate can see this class, so
  # it is flagged structurally, like bridge outages.
  structural_flag = false
  if it.kind === :gen && types1 === sc.bus_types
    ps0 = template.prosumpsVec[it.internal]
    bus0 = getPosumerBusIndex(ps0)
    if 1 <= bus0 <= length(sc.bus_types) && sc.bus_types[bus0] !== :Slack
      regulating_here = any(isGenerator(other) && isRegulating(other) && getPosumerBusIndex(other) == bus0 for other in template.prosumpsVec)
      regulating_here && (structural_flag = true)
    end
  end
  F1 = mismatch_rectangular(Yb, sc.V0, S1, types1, sc.Vset, sc.slack_idx; dslack = ds1)
  m1 = maximum(abs, F1)
  J1 = try
    _screen_freeze_isolated!(build_rectangular_jacobian_pq_pv(Yb, sc.V0, types1, sc.Vset, sc.slack_idx; dslack = ds1), sc.non_slack, sc.iso)
  catch err
    _rethrow_unless_solver_failure(err)
    return _screen_unscreenable(:numeric)
  end
  rhs = Vector{Float64}(undef, length(F1))
  @. rhs = -F1
  dx = _screen_woodbury_solve(sc.lu0, sc.J0, J1, rhs)
  dx === nothing && return _screen_unscreenable(:numeric)
  Vest = _apply_rectangular_delta(sc.V0, dx, sc.slack_idx, sc.non_slack, 1.0)
  all(isfinite, reinterpret(Float64, Vest)) || return _screen_unscreenable(:numeric)
  # the post-step residual is evaluated at the FULL step state, for the
  # augmented system including the estimated lambda
  lambda_est = ds1 === nothing ? 0.0 : dx[end]
  ds1 === nothing || (ds1.lambda_trial = lambda_est)
  F2 = mismatch_rectangular(Yb, Vest, S1, types1, sc.Vset, sc.slack_idx; dslack = ds1)
  m2 = maximum(abs, F2)
  reduced = m2 < m1
  # estimated loadings over the rated, fully closed branches (the estimate
  # mirrors _contingency_metrics but reads the trial state, writes nothing)
  max_loading = NaN
  loading_flag = false
  est_overloads = OverloadRecord[]
  for (bi, br2) in enumerate(template.branchVec)
    bi == outaged_branch && continue
    br2.status == 1 || continue
    _branch_terminal_state(br2) === :closed || continue
    rating = br2.sn_MVA
    (rating === nothing || !isfinite(rating) || rating <= 0.0) && continue
    from = Int(br2.fromBus)
    to = Int(br2.toBus)
    (from in sc.iso || to in sc.iso) && continue
    s_from = abs(_closed_branch_flow_pu(Vest, from, to, br2, 1)) * template.baseMVA
    s_to = abs(_closed_branch_flow_pu(Vest, to, from, br2, 2)) * template.baseMVA
    s_mva = max(s_from, s_to)
    loading = 100.0 * s_mva / rating
    (isnan(max_loading) || loading > max_loading) && (max_loading = loading)
    # Change-based loading flag: the error of the one-step estimate grows
    # with the change the outage causes on a branch, the base value is the
    # converged solution. A branch flags when its estimate plus that change
    # (at least 1 point, at most the margin) reaches 100 percent. Without a
    # base value the full margin applies, as before. A branch the base case
    # already loads near its limit no longer flags every outage of the
    # network (case_ACTIVSg2000: a transformer at 92.7 percent flagged all
    # 2756 screenable outages; zero of 3638 missed with the rule, and the
    # calibration networks case118, case300 and the tiled grids unchanged).
    base_b = sc.base_loading[bi]
    unc = isnan(base_b) ? engine.screening_margin_pct : min(engine.screening_margin_pct, max(1.0, abs(loading - base_b)))
    loading + unc >= 100.0 && (loading_flag = true)
    if loading > 100.0
      cname = getCompName(br2.comp)
      base = get(engine.base_loadings, cname, NaN)
      push!(est_overloads, OverloadRecord(cname, loading, base, loading - base, s_mva, rating))
    end
  end
  sort!(est_overloads; by = o -> -o.loading_pct)
  vmin = Inf
  vmax = -Inf
  est_violations = String[]
  bus_names = Dict{Int,String}(idx => name for (name, idx) in template.busDict)
  for node in template.nodeVec
    node.busIdx in sc.iso && continue
    vm = abs(Vest[node.busIdx])
    isfinite(vm) || continue
    vmin = min(vmin, vm)
    vmax = max(vmax, vm)
    if vm < engine.vm_min_pu || vm > engine.vm_max_pu
      push!(est_violations, get(bus_names, node.busIdx, getCompName(node.comp)))
    end
  end
  vmin = isfinite(vmin) ? vmin : NaN
  vmax = isfinite(vmax) ? vmax : NaN
  # flag rule: near a loading limit, near a voltage band limit
  # (margin_pct of the band width), a screening step that failed to reduce
  # the mismatch, or a post-step residual too large to trust the
  # linearization (the case118 acceptance run showed a 0.047 pu residual
  # hiding a 0.05 pu voltage-drop underestimate; 0.02 pu keeps a safety
  # factor of two below that while releasing the well-conditioned bulk)
  vband = engine.vm_max_pu - engine.vm_min_pu
  vmargin = engine.screening_margin_pct / 100.0 * vband
  flagged = structural_flag
  flagged |= !reduced
  flagged |= !(m2 <= _SCREENING_TRUST_MISMATCH_PU)
  # distributed slack (4b): limits are NOT applied in the estimate; a
  # participating unit whose estimated output leaves its P band flags the
  # scenario, so the water filling of the full solve handles exactly those
  if ds1 !== nothing
    for row in ds1.gen_table
      (1 <= row.gen_index <= length(template.prosumpsVec)) || continue
      ps2 = template.prosumpsVec[row.gen_index]
      p_est = something(ps2.pVal, 0.0) + row.alpha * lambda_est * template.baseMVA
      lo = ps2.minP
      hi = ps2.maxP
      if (lo !== nothing && isfinite(lo) && p_est < lo) || (hi !== nothing && isfinite(hi) && p_est > hi)
        flagged = true
        break
      end
    end
  end
  flagged |= loading_flag
  flagged |= !isnan(vmin) && vmin <= engine.vm_min_pu + vmargin
  flagged |= !isnan(vmax) && vmax >= engine.vm_max_pu - vmargin
  estimate = (max_loading_pct = max_loading, vmin_pu = vmin, vmax_pu = vmax, mismatch_start = m1, mismatch_after = m2)
  return (flagged = flagged, estimate = estimate, est_overloads = est_overloads, est_violations = est_violations, reason = structural_flag ? :structural_pv_gen : :ok)
end

# result row for a scenario that STAYED screened (no full run): the metric
# columns carry the estimates, start_used reads :screen
function _screened_result(engine::ScenarioEngine, it::_ScenarioOutageItem, scr)::ContingencyResult
  est = scr.estimate
  severity = it.weight * max(0.0, isnan(est.max_loading_pct) ? 0.0 : est.max_loading_pct - 100.0)
  sc = engine.screen
  return ContingencyResult(it.name, it.weight, true, 0, :screen, est.vmax_pu, est.vmin_pu, est.max_loading_pct, severity, scr.est_overloads, scr.est_violations, sc === nothing ? 0 : sc.base_island_count, 0.0, nothing)
end

"""
    ScenarioResult

Outcome of one scenario under an engine with screening: the
full [`ContingencyResult`](@ref) surface (all its fields forward)
plus `screened` (true when the full run was skipped and the metric columns
carry the screening estimates) and `screening_estimate` (the estimate
tuple: estimated worst loading, voltage envelope, and the screening-step
mismatch before and after; `nothing` when no estimate was computed).
"""
struct ScenarioResult
  result::ContingencyResult
  screened::Bool
  screening_estimate::Union{Nothing,NamedTuple}
end

function Base.getproperty(r::ScenarioResult, f::Symbol)
  f in (:result, :screened, :screening_estimate) && return getfield(r, f)
  return getfield(getfield(r, :result), f)
end

Base.propertynames(::ScenarioResult) = (fieldnames(ContingencyResult)..., :result, :screened, :screening_estimate)

"""
    writeContingencyResultsCSV(path, results::Vector{ScenarioResult}; format = "technical")

The [`ContingencyResult`](@ref) CSV with the two screening columns appended:
`screened` and the compact estimate
`max_loading_pct|vmin_pu|vmax_pu` (empty when no estimate was computed).
With screening `:off` the engine returns plain `ContingencyResult` rows and
the classic writer runs with the same `format`, so the two stay in the same
shape for one run. `format` (`technical`/`excel_de`/`excel_us`) sets the
delimiter and number formatting like every other CSV artifact a run writes
(issue #376); `excel_de` reproduces the previous hardcoded semicolon
delimiter.
"""
function writeContingencyResultsCSV(path::AbstractString, results::Vector{ScenarioResult}; format = result_csv_format())
  header = ("name", "weight", "converged", "iterations", "start_used", "min_vm_pu", "max_vm_pu", "max_branch_loading_pct", "severity", "overloads", "voltage_violations", "island_count", "shed_load_mw", "error", "screened", "screening_estimate")
  resolved_format = _resolve_detailed_csv_format(format)
  rows = (
    begin
      r = sr.result
      overloads = join(["$(o.name)@$(round(o.loading_pct; digits = 1))" for o in r.overloads], ",")
      # the screening estimate is one cell with its own separator, so its
      # numbers are formatted here and travel as text
      est = sr.screening_estimate === nothing ? "" : string(_format_csv_number(sr.screening_estimate.max_loading_pct, resolved_format), "|", _format_csv_number(sr.screening_estimate.vmin_pu, resolved_format), "|", _format_csv_number(sr.screening_estimate.vmax_pu, resolved_format))
      (
        r.name,
        r.weight,
        r.converged,
        r.iterations,
        String(r.start_used),
        r.min_vm_pu,
        r.max_vm_pu,
        r.max_branch_loading_pct,
        r.severity,
        overloads,
        join(r.voltage_violations, ","),
        r.island_count,
        r.shed_load_mw,
        r.error === nothing ? "" : r.error,
        sr.screened,
        est,
      )
    end for sr in results
  )
  write_result_csv(path, header, rows; format = format)
  return path
end

printContingencyResults(results::Vector{ScenarioResult}; max_rows::Int = 50, sort_by::Symbol = :severity) = printContingencyResults(stdout, results; max_rows = max_rows, sort_by = sort_by)

function printContingencyResults(io::IO, results::Vector{ScenarioResult}; max_rows::Int = 50, sort_by::Symbol = :severity)
  screened = count(r -> r.screened, results)
  screened > 0 && println(io, "screened rows (start = screen): ", screened, " of ", length(results), " carry first-order ESTIMATES in the metric columns, not full solves")
  return printContingencyResults(io, ContingencyResult[r.result for r in results]; max_rows = max_rows, sort_by = sort_by)
end

buildContingencyReport(results::AbstractVector{ScenarioResult}; top::Int = 10) = buildContingencyReport(ContingencyResult[r.result for r in results]; top = top)

# --- evaluation ---------------------------------------------------------------

# What a catch of this file may turn into an outcome is a failure of the
# SOLVE: the island validation, a singular matrix, a participation set that
# cannot be built, a diverged iteration. An interrupt and an abort go on as
# they are, and so does a defect of the calling code (a keyword the solver
# does not know, a missing binding, a wrong type or field): reported as
# "not converged", "islanded" or "screening unavailable" it would look like
# a finding about the network. An unknown keyword did read as "islanded
# without reference", because the error text lists the solver's keywords.
function _rethrow_unless_solver_failure(err)
  (err isa InterruptException || err isa PowerFlowAborted) && rethrow(err)
  (err isa MethodError || err isa UndefKeywordError || err isa UndefVarError || err isa TypeError || err isa FieldError) && rethrow(err)
  return nothing
end

# The note and the disconnected load of the buses an outage cut off
# (`nothing` and 0.0 when it cut off none). The text keeps the leading
# "islanded" token the result filters look for; the scheduled values of
# the prosumers at those buses are summed, a generating unit counts with
# its scheduled active power (0 MW for a synchronous condenser).
function _cut_off_bus_note(net::Net, cut_off::Vector{Int})
  isempty(cut_off) && return (nothing, 0.0)
  names = Dict{Int,String}(idx => bus for (bus, idx) in net.busDict)
  buses = Set(cut_off)
  load_mw = 0.0
  gen_mw = 0.0
  units = 0
  for ps in net.prosumpsVec
    Int(ps.comp.cFrom_bus) in buses || continue
    p = something(ps.pVal, 0.0)
    if isGenerator(ps)
      units += 1
      gen_mw += p
    else
      load_mw += p
    end
  end
  listed = join((get(names, i, string(i)) for i in first(cut_off, 5)), ", ")
  length(cut_off) > 5 && (listed *= ", ...")
  parts = String["$(round(load_mw; digits = 1)) MW load disconnected"]
  units > 0 && push!(parts, "$(units) generating unit(s) with $(round(gen_mw; digits = 1)) MW out of service")
  note = "islanded: $(length(cut_off) == 1 ? "bus" : "buses") $(listed) cut off, $(join(parts, ", ")); the remaining network solved"
  return (note, load_mw)
end

# The note of a case that solved on a reference the outage made necessary.
# It follows the note of cut-off buses where both apply. The wording keeps
# clear of the token the islanding filter looks for: a taken-over
# reference alone is a solved case, not a loss of load.
function _reference_note(net::Net, buses::Vector{Int}, cut_off_note)
  isempty(buses) && return cut_off_note
  names = Dict{Int,String}(idx => bus for (bus, idx) in net.busDict)
  text = "reference taken over by $(length(buses) == 1 ? "bus" : "buses") $(join((get(names, b, string(b)) for b in buses), ", ")) (the outage removed the reference)"
  return cut_off_note === nothing ? text : string(cut_off_note, "; ", text)
end

# The post-outage evaluation, verbatim from the pre-engine
# _run_one_contingency: island marking and precheck, the per-case
# start-value ladder over a warm-start snapshot, metrics on success. The
# snapshot is taken AFTER markIsolatedBuses! on purpose: isolation zeroes
# the isolated bus voltages, and the warm start of every ladder stage must
# include those zeroes (bitwise with the old per-case deepcopy path).
function _evaluate_outaged_net!(engine::ScenarioEngine, work::Net, name::String, weight::Float64; cold::Bool = false)
  # warm active set: the cold retry solves without the base case's clamps
  pf_kw = cold ? Base.structdiff(engine.pf_kwargs, NamedTuple{(:qlimit_warm_clamps,)}) : engine.pf_kwargs
  # buses the outage itself cuts off: isolated after the marking, not before.
  # A single bus without a connection is no island for the detector below
  # (it drops isolated buses), so without this the outage of a radial
  # branch would read as a clean case although its load is gone or its
  # unit is out of service.
  isolated_before = Set{Int}(i for i in eachindex(work.nodeVec) if isIsolated(work.nodeVec[i]))
  markIsolatedBuses!(net = work, log = false)
  cut_off = Int[i for i in eachindex(work.nodeVec) if isIsolated(work.nodeVec[i]) && !(i in isolated_before)]
  cut_off_note, cut_off_load_mw = _cut_off_bus_note(work, cut_off)
  # the same detection the solve runs: an island without a slack or PV bus
  # takes its best generating unit (stated reference priority first, then
  # the strongest) regardless of auto_slack
  island_report = detect_ac_islands(work; promote_generators = true)
  island_count = length(island_report.rows)
  # rows for reference-less islands, kept so an islanding failure is
  # reported specifically (load-only vs. stranded generation). The solver
  # stays the authority on solvability: it may link-merge synchronously
  # tied islands this pre-merge snapshot cannot see.
  refless = [r for r in island_report.rows if r.chosen_ref_bus == 0]
  snap = _snapshot_start_voltages(work)
  template_flatstart = work.flatstart
  total_it = 0
  last_error = nothing
  for stage in engine.ladder
    _restore_start_voltages!(work, snap)
    it_stage = 0
    try
      erg = 1
      if stage === :warm
        work.flatstart = template_flatstart
        it_stage, erg = runpf!(work, engine.maxIte, engine.tol, 0; islands_enabled = true, pf_kw...)
      elseif stage === :flat
        work.flatstart = true
        it_stage, erg = runpf!(work, engine.maxIte, engine.tol, 0; islands_enabled = true, pf_kw...)
      elseif stage === :dc
        # flat magnitudes with DC-projected start angles
        work.flatstart = false
        _dc_seed_rectangular_angles!(work, PowerFlowConfig(max_iter = engine.maxIte, tol = engine.tol))
        it_stage, erg = runpf!(work, engine.maxIte, engine.tol, 0; islands_enabled = true, pf_kw...)
      elseif stage === :apslf
        # APSLF start values, one bounded attempt on the same keywords as
        # every other stage. A default PowerFlowConfig here forwarded only
        # auto_slack and solved this stage on a single slack and default
        # Q-limit settings whatever the batch asked for (#456).
        work.flatstart = false
        it_stage, erg = runpf!(work, engine.maxIte, engine.tol, 0; islands_enabled = true, pf_kw..., apslf_start_enabled = true, apslf_start_order = 40)
      end
      total_it += it_stage
      if erg == 0
        calcNetLosses!(work)
        m = _contingency_metrics(work, engine.vm_min_pu, engine.vm_max_pu, engine.base_loadings)
        severity = weight * max(0.0, isnan(m.max_loading) ? 0.0 : m.max_loading - 100.0)
        # the remaining network solved; what the outage cut off is part of
        # the result (the load it disconnected, the note naming the buses)
        # a reference the outage removed and another unit took over: the
        # unit auto_slack promoted on the whole network, or the bus an
        # island was given; named next to what the outage cut off
        # read AFTER the solve: the unit auto_slack promoted is a slack of
        # the working copy by now, and an island without one is given its
        # bus by the same detection the solve used
        taken = Set{Int}(Int(ps.comp.cFrom_bus) for ps in work.prosumpsVec if isSlack(ps))
        for row in detect_ac_islands(work; promote_generators = true).rows
          (row.n_ref == 0 && row.chosen_ref_bus > 0) && push!(taken, Int(row.chosen_ref_bus))
        end
        setdiff!(taken, engine.base_reference_buses, cut_off, isolated_before)
        note = _reference_note(work, sort!(collect(taken)), cut_off_note)
        return ContingencyResult(name, weight, true, total_it, stage, m.vmax, m.vmin, m.max_loading, severity, m.overloads, m.violations, island_count, cut_off_load_mw, note)
      end
      last_error = "power flow did not converge (status $(erg))"
    catch err
      _rethrow_unless_solver_failure(err)
      total_it += it_stage
      msg = sprint(showerror, err)
      # an outage that splits off a reference-less island surfaces as the
      # island-validation error; no start recipe fixes a topology fact
      if occursin("reference", msg) && (occursin("island", msg) || occursin("Island", msg))
        return _islanded_contingency_result(name, weight, island_count, refless)
      end
      last_error = first(split(msg, '\n'))
    end
  end
  return ContingencyResult(name, weight, false, total_it, :none, NaN, NaN, NaN, NaN, OverloadRecord[], String[], island_count, 0.0, last_error === nothing ? "power flow did not converge" : last_error)
end

"""
    evaluate!(engine, worker, item) -> ContingencyResult

Evaluate one work item on the worker's net and reset the worker
afterwards. An outage item reproduces the pre-engine contingency
semantics exactly: a branch outage opens both terminals (equivalent to
the old structural removal for the Y-bus, isolation marking, and
metrics, all of which treat a fully open branch as absent), and a
generator outage runs `_remove_contingency_generator!` with the deleted
prosumer reinserted afterwards. A patch item applies its operations
through [`apply!`](@ref) with the documented scenario semantics.
"""
function evaluate!(engine::ScenarioEngine, worker::ScenarioWorker, it::_ScenarioOutageItem)::ContingencyResult
  result = _evaluate_outage_item!(engine, worker, it; cold = false)
  # warm active set (0.30.2): an outage that does not converge from the
  # base case's clamps is solved again on a freshly reset worker from the
  # file's PV/PQ state, so the switch never loses an outage the plain start
  # solves (case300: one outage converged only without it)
  if !result.converged && haskey(engine.pf_kwargs, :qlimit_warm_clamps)
    cold = _evaluate_outage_item!(engine, worker, it; cold = true)
    cold.converged && return _warm_cold_result(cold)
  end
  # cold check (0.30.3): a converged warm result near the lower voltage
  # limit or with a violation is solved cold as well, the less favourable
  # result counts (where the active set is not unique the warm start can
  # land on a more favourable valid solution than the cold one)
  if _warm_needs_cold_check(engine, result)
    cold = _evaluate_outage_item!(engine, worker, it; cold = true)
    _contingency_less_favourable(cold, result) && return _warm_checked_result(cold, result)
  end
  return result
end

function _evaluate_outage_item!(engine::ScenarioEngine, worker::ScenarioWorker, it::_ScenarioOutageItem; cold::Bool)::ContingencyResult
  work = worker.net
  removed_prosumer = nothing
  if it.kind === :gen
    removed_prosumer = work.prosumpsVec[it.internal]
    _remove_contingency_generator!(work, it.internal)
  else
    br = work.branchVec[it.internal]
    br.status = 0
    br.from_status = 0
    br.to_status = 0
    # a bus shunt part that stems from this branch (MATPOWER reimport of an
    # asymmetric branch shunt) leaves with it; the worker reset restores it
    _remove_branch_shunt_parts!(work, it.internal)
  end
  result = _evaluate_outaged_net!(engine, work, it.name, it.weight; cold = cold)
  removed_prosumer === nothing || insert!(work.prosumpsVec, it.internal, removed_prosumer)
  _reset_scenario_worker!(work, engine.template)
  return result
end

# a cold-retry result: start_used names the retry, the note says why
# The cold check of the warm active set (0.30.3): only a CONVERGED warm
# result is checked (a failed one already gets the cold retry), and only
# where it matters: the lowest voltage within the margin above the lower
# limit, or any overload or voltage violation.
function _warm_needs_cold_check(engine::ScenarioEngine, r::ContingencyResult)::Bool
  engine.warm_cold_check || return false
  haskey(engine.pf_kwargs, :qlimit_warm_clamps) || return false
  r.converged || return false
  return r.min_vm_pu < engine.vm_min_pu + engine.warm_cold_check_margin_pu || !isempty(r.overloads) || !isempty(r.voltage_violations)
end

# "less favourable", in this order: a converged result against a failed one
# never loses (a failed cold solve leaves the warm result standing); more
# violations (overloads plus voltage violations); then the lower lowest
# voltage; then the higher largest branch loading. Voltages within 1e-6 pu
# and loadings within 0.01 percentage points count as equal (the tolerance
# of the N-1 power-mode comparison, #459): two solves of the same solution
# differ in the last digits, and without it the rounding noise decided
# (sp_case118: 25 of 240 outages took a cold result identical to the warm).
const _WARM_CHECK_VM_TOL_PU = 1.0e-6
const _WARM_CHECK_LOADING_TOL_PCT = 0.01
function _contingency_less_favourable(a::ContingencyResult, b::ContingencyResult)::Bool
  a.converged || return false
  b.converged || return true
  va = length(a.overloads) + length(a.voltage_violations)
  vb = length(b.overloads) + length(b.voltage_violations)
  va != vb && return va > vb
  abs(a.min_vm_pu - b.min_vm_pu) > _WARM_CHECK_VM_TOL_PU && return a.min_vm_pu < b.min_vm_pu
  return a.max_branch_loading_pct > b.max_branch_loading_pct + _WARM_CHECK_LOADING_TOL_PCT
end

function _warm_checked_result(cold::ContingencyResult, warm::ContingencyResult)::ContingencyResult
  why = string("warm active set result (min Vm ", round(warm.min_vm_pu; digits = 5), " pu) replaced by the less favourable cold solve")
  note = (cold.error === nothing || isempty(cold.error)) ? why : string(cold.error, "; ", why)
  return ContingencyResult(cold.name, cold.weight, cold.converged, cold.iterations, :warm_cold_check, cold.max_vm_pu, cold.min_vm_pu, cold.max_branch_loading_pct, cold.severity, cold.overloads, cold.voltage_violations, cold.island_count, cold.shed_load_mw, note)
end

function _warm_cold_result(r::ContingencyResult)::ContingencyResult
  why = "warm active set did not converge, solved from the file's PV/PQ state"
  note = (r.error === nothing || isempty(r.error)) ? why : string(r.error, "; ", why)
  return ContingencyResult(r.name, r.weight, r.converged, r.iterations, :warm_cold, r.max_vm_pu, r.min_vm_pu, r.max_branch_loading_pct, r.severity, r.overloads, r.voltage_violations, r.island_count, r.shed_load_mw, note)
end

evaluate!(::ScenarioEngine, ::ScenarioWorker, it::_ScenarioFailedItem)::ContingencyResult = ContingencyResult(it.name, it.weight, false, 0, :none, NaN, NaN, NaN, NaN, OverloadRecord[], String[], 0, 0.0, it.message)

function evaluate!(engine::ScenarioEngine, worker::ScenarioWorker, it::_ScenarioPatchItem)::ContingencyResult
  engine.index === nothing && return ContingencyResult(it.scenario.name, it.scenario.weight, false, 0, :none, NaN, NaN, NaN, NaN, OverloadRecord[], String[], 0, 0.0, "patch scenarios need an engine built with a ScenarioIndex")
  work = worker.net
  undo = apply!(work, it.scenario.ops, engine.index)
  result = _evaluate_outaged_net!(engine, work, it.scenario.name, it.scenario.weight)
  restore!(work, undo)
  _reset_scenario_worker!(work, engine.template)
  # warm active set: the same cold retry as for an outage item
  if !result.converged && haskey(engine.pf_kwargs, :qlimit_warm_clamps)
    undo = apply!(work, it.scenario.ops, engine.index)
    cold = _evaluate_outaged_net!(engine, work, it.scenario.name, it.scenario.weight; cold = true)
    restore!(work, undo)
    _reset_scenario_worker!(work, engine.template)
    cold.converged && return _warm_cold_result(cold)
  end
  # the same cold check as for an outage item
  if _warm_needs_cold_check(engine, result)
    undo = apply!(work, it.scenario.ops, engine.index)
    cold = _evaluate_outaged_net!(engine, work, it.scenario.name, it.scenario.weight; cold = true)
    restore!(work, undo)
    _reset_scenario_worker!(work, engine.template)
    _contingency_less_favourable(cold, result) && return _warm_checked_result(cold, result)
  end
  return result
end

# --- batch --------------------------------------------------------------------

# One batch item under the active screening mode: an outage item is
# screened first (when a screening state exists); a scenario that stays
# below every flag threshold keeps its estimate as the result row, a
# flagged one gets the full run WITH the estimate attached, and anything
# unscreenable (bridge outage, slack unit, patch scenario, failed screen)
# goes straight to the full run. Only the full-run branch touches the
# worker, so screened scenarios never pay the worker reset.
function _evaluate_item(engine::ScenarioEngine, worker::ScenarioWorker, it)::ScenarioResult
  if it isa _ScenarioOutageItem && engine.screen !== nothing
    scr = _screen_outage(engine, it)
    if scr.estimate !== nothing
      if engine.screening_mode === :only || !scr.flagged
        return ScenarioResult(_screened_result(engine, it, scr), true, scr.estimate)
      end
      return ScenarioResult(evaluate!(engine, worker, it), false, scr.estimate)
    end
  end
  return ScenarioResult(evaluate!(engine, worker, it), false, nothing)
end

# Progress reporting of a batch (Web UI N-1 counter). The
# caller's callable sits in an `Any` field and the keyword is
# @nospecialize, so the batch runner compiles ONCE, not once per closure
# type (compile-time rule: a callback keyword goes behind @nospecialize or
# a concrete struct); the dynamic call costs well under a
# microsecond against a full power-flow solve per item. `done` is an
# atomic counter because the parallel branch finishes items on several
# threads at once: every finished item increments it exactly once and the
# callback receives the value AFTER its own increment, so the values the
# callback sees across a batch are exactly 1:total, each once. On threads
# they may arrive out of order (4 can be reported after 5); a consumer that
# displays a counter keeps the maximum.
struct _BatchProgress
  callback::Any                  # nothing, or callable (done::Int, total::Int)
  done::Threads.Atomic{Int}
  total::Int
end

_BatchProgress(callback, total::Int) = _BatchProgress(callback, Threads.Atomic{Int}(0), total)

# one finished item: count it and tell the caller. Screened items, failed
# items (unknown element) and solved items all count, the base-case solve
# never does (it runs in the ScenarioEngine constructor, before the batch).
function _progress_tick!(p::_BatchProgress)
  done = Threads.atomic_add!(p.done, 1) + 1   # atomic_add! returns the OLD value
  p.callback === nothing || p.callback(done, p.total)
  return nothing
end

# The batch runner: chunked over Julia threads exactly like the old
# runContingencies! fan-out, with the serial path running the SAME
# evaluate! on one worker (threads rule: serial path is the same function).
# The per-chunk workers are created SERIALLY before spawning (threads
# rule: never create shared-resource copies under threads), and worker
# state is indexed by CHUNK, never by threadid. With screening :off the
# batch returns the historical Vector{ContingencyResult}; otherwise every
# row is a ScenarioResult.
function _run_engine_batch(engine::ScenarioEngine, items::Vector{_ScenarioItem}; parallel_enabled::Union{Nothing,Bool} = nothing, parallel_max_tasks::Union{Nothing,Int} = nothing, parallel_min_work_items::Union{Nothing,Int} = nothing, @nospecialize(progress = nothing))
  screening_off = engine.screening_mode === :off
  isempty(items) && return screening_off ? ContingencyResult[] : ScenarioResult[]
  eval_one = screening_off ? evaluate! : _evaluate_item
  results = Vector{Any}(undef, length(items))
  # shared by all chunks; only read and atomically incremented in the tasks
  tracker = _BatchProgress(progress, length(items))
  parallel_on, parallel_cap, parallel_min_items = _resolve_parallel_runtime(parallel_enabled, parallel_max_tasks, parallel_min_work_items)
  use_parallel = parallel_on && Threads.nthreads() > 1 && parallel_cap > 1 && length(items) >= parallel_min_items
  if use_parallel
    chunks = collect(Iterators.partition(eachindex(items), cld(length(items), parallel_cap)))
    chunk_workers = [ScenarioWorker(deepcopy(engine.template)) for _ in chunks]
    tasks = [
      Threads.@spawn begin
        for idx in chunks[ci]
          # an N-1 batch is the longest run the Web UI offers; the abort
          # token is process-wide, so it is visible inside this task
          sparlectra_check_abort()
          results[idx] = eval_one(engine, chunk_workers[ci], items[idx])
          # no assignment here: the task body must not rebind enclosing names
          _progress_tick!(tracker)
        end
      end for ci in eachindex(chunks)
    ]
    foreach(wait, tasks)
  else
    serial_worker = ScenarioWorker(deepcopy(engine.template))
    for idx in eachindex(items)
      sparlectra_check_abort()
      results[idx] = eval_one(engine, serial_worker, items[idx])
      _progress_tick!(tracker)
    end
  end
  screening_off && return ContingencyResult[results[i] for i in eachindex(items)]
  return ScenarioResult[results[i] for i in eachindex(items)]
end

# resolve ContingencyCases against the template ONCE (the old code
# resolved on the fresh per-case deepcopy, which carries the same indices)
function _engine_items_from_cases(template::Net, cases::Vector{ContingencyCase})::Vector{_ScenarioItem}
  items = _ScenarioItem[]
  for case in cases
    if case.kind === :gen
      pidx = _resolve_contingency_generator(template, case.element)
      push!(items, pidx === nothing ? _ScenarioFailedItem(case.name, case.weight, "unknown generator $(case.element) (not found in net.prosumpsVec)") : _ScenarioOutageItem(case.name, case.weight, :gen, pidx))
    else
      idx = _resolve_contingency_branch(template, case.element)
      push!(items, idx === nothing ? _ScenarioFailedItem(case.name, case.weight, "unknown branch $(case.element) (not found in net.branchVec)") : _ScenarioOutageItem(case.name, case.weight, :branch, idx))
    end
  end
  return items
end

# A single status-0 op IS an N-1 outage and gets the contingency-exact
# path (that identity is what makes N-1 the special case of the model);
# everything else is a patch scenario with the documented apply! semantics.
function _engine_item_from_scenario(s::Scenario, index::ScenarioIndex)::_ScenarioItem
  if length(s.ops) == 1
    op = s.ops[1]
    if op.op === :status && op.value == 0.0 && op.target in (:branch, :transformer, :generator, :external_grid)
      internal = get(index.internal_by_id, op.id, nothing)
      internal === nothing && return _ScenarioFailedItem(s.name, s.weight, string("component id ", op.id, " has no internal index"))
      return _ScenarioOutageItem(s.name, s.weight, op.target in (:branch, :transformer) ? :branch : :gen, internal)
    end
  end
  return _ScenarioPatchItem(s)
end

"""
    runScenarios!(net, set::ScenarioSet; index::ScenarioIndex, kwargs...) -> Vector{ContingencyResult}
    runScenarios!(net, scenarios::Vector{Scenario}; index::ScenarioIndex, kwargs...)

Evaluate a scenario set on `net`: validate against
the index and the net (active tap controllers reject a `tap_pos` patch),
expand the N-1 modes through the existing generators, and run every
scenario on the engine. A scenario that is a single status-0 outage
reproduces [`runContingencies!`](@ref) exactly; general patch scenarios
(setpoint changes, load scalings, multi-op outages) run the same ladder
and metrics on the patched working copy. `net` is never mutated. The
keyword surface matches `runContingencies!` (`vm_min_pu`, `vm_max_pu`,
`maxIte`, `tol`, `rescue_ladder`, the `parallel_*` overrides); remaining
keywords reach the per-scenario power-flow solves. Results are returned
in scenario order; failures are reported in the result, never thrown.
`progress` reports the batch progress exactly as in
[`runContingencies!`](@ref) (one call per finished scenario, thread safe
callable required).

Screening: with `screening_mode = :flag` every
non-islanding single outage is estimated first with one Woodbury-corrected
Newton step on the base factorization, and only scenarios whose estimate
comes within `screening_margin_pct` of a limit (or whose screening step
did not reduce the mismatch) get the full solve; `:only` reports the
estimates without full runs where an estimate exists. The programmatic
default is `:off` (every scenario fully solved, bit-identical to the
pre-screening engine); the service path applies the
`contingency.screening.*` configuration, whose default is `:off` as well
(screening is a deliberate opt-in). With
screening active the return type is `Vector{ScenarioResult}` (the
`ContingencyResult` surface plus `screened` and `screening_estimate`).
"""
function runScenarios!(
  net::Net,
  set::ScenarioSet;
  index::ScenarioIndex,
  vm_min_pu::Float64 = 0.9,
  vm_max_pu::Float64 = 1.1,
  maxIte::Int = 30,
  tol::Float64 = 1e-8,
  rescue_ladder::Vector{Symbol} = [:warm],
  screening_mode::Symbol = :off,
  screening_margin_pct::Float64 = 10.0,
  parallel_enabled::Union{Nothing,Bool} = nothing,
  parallel_max_tasks::Union{Nothing,Int} = nothing,
  parallel_min_work_items::Union{Nothing,Int} = nothing,
  auto_slack::Bool = true,
  warm_active_set::Bool = true,
  warm_cold_check::Bool = true,
  warm_cold_check_margin_pu::Float64 = DEFAULT_WARM_COLD_CHECK_MARGIN_PU,
  warm_note::Union{Nothing,Base.RefValue{String}} = nothing,
  @nospecialize(progress = nothing),
  kwargs...,
)
  validate_scenarios(set, index, net)
  scenarios = expand_scenarios(set, net, index)
  isempty(scenarios) && return screening_mode === :off ? ContingencyResult[] : ScenarioResult[]
  ladder = _validate_contingency_ladder(rescue_ladder; context = "runScenarios!: rescue_ladder")
  engine = ScenarioEngine(net; vm_min_pu = vm_min_pu, vm_max_pu = vm_max_pu, maxIte = maxIte, tol = tol, ladder = ladder, index = index, pf_kwargs = (; auto_slack = auto_slack, kwargs...), screening_mode = screening_mode, screening_margin_pct = screening_margin_pct, warm_active_set = warm_active_set, warm_note = warm_note, warm_cold_check = warm_cold_check, warm_cold_check_margin_pu = warm_cold_check_margin_pu)
  items = _ScenarioItem[_engine_item_from_scenario(s, index) for s in scenarios]
  return _run_engine_batch(engine, items; parallel_enabled = parallel_enabled, parallel_max_tasks = parallel_max_tasks, parallel_min_work_items = parallel_min_work_items, progress = progress)
end

runScenarios!(net::Net, scenarios::Vector{Scenario}; kwargs...) = runScenarios!(net, ScenarioSet(scenarios = scenarios); kwargs...)

# Take the branch-derived parts of every bus shunt that stem from branch
# `bidx` out of the shunt admittance of the working copy (0.20.0): the
# MATPOWER importer records the terminal excess of an asymmetric branch
# shunt as a part of the bus shunt, and an outage of the branch removes
# that admittance with the branch. The parts stay recorded, so the worker
# reset (a copy of the template) restores the admittance.
function _remove_branch_shunt_parts!(work::Net, bidx::Int)
  for sh in work.shuntVec
    part = get(sh.branch_parts, bidx, nothing)
    part === nothing && continue
    sh.y_pu_shunt -= part
  end
  return nothing
end
