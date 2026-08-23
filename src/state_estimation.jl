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

# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 1.1.2026
# file: src/state_estimation.jl
# purpose: classical WLS state estimation (runse!) with observability and
#          redundancy analysis and an optional PMU reference-angle offset
#          state

using Printf
"""
    SEResult

Result container for the first classical WLS state-estimation run.

`vaRefOffsetDeg` is the estimated common PMU reference-angle offset α in
degrees (the slack-bus angle expressed in the PMU time base). It is `nothing`
when no offset state was part of the estimation, i.e. when there are no
active `VaMeas` measurements or `pmu_ref_offset = :off`.
"""
struct SEResult
  voltages::Vector{ComplexF64}
  converged::Bool
  iterations::Int
  residualNorm::Float64
  residuals::Vector{Float64}
  objectiveJ::Float64
  dof::Int # degrees of freedom ν=m−n
  jWithin3Sigma::Bool
  vaRefOffsetDeg::Union{Nothing,Float64}
end

"""
    numeric_rank(A; tol=nothing) -> Int

Compute numerical matrix rank using singular values.
"""
function numeric_rank(A::AbstractMatrix{<:Real}; tol = nothing)
  s = svdvals(Matrix(A))
  isempty(s) && return 0
  if tol === nothing
    tol = eps(eltype(s)) * max(size(A)...) * maximum(s)
  end
  return count(>(tol), s)
end

## Internal helper: find slack bus index in the current network.
@inline function _find_slack_idx(net::Net)
  for (i, node) in enumerate(net.nodeVec)
    if getNodeType(node) == Slack
      return i
    end
  end
  error("No slack bus found for state estimation")
end

## Internal helper: true when the active set contains PMU angle measurements.
@inline function _has_active_va(measurements::Vector{Measurement})
  return any(m -> m.active && m.typ == VaMeas, measurements)
end

## Internal helper: decide whether the PMU reference-offset state α is carried.
## :auto adds α as soon as active VaMeas rows exist; :off treats PMU angles as
## already slack-referenced.
@inline function _va_offset_active(measurements::Vector{Measurement}, mode::Symbol)
  mode == :off && return false
  mode == :auto || error("pmu_ref_offset mode must be :auto or :off, got :$(mode)")
  return _has_active_va(measurements)
end

## Internal helper: map estimator state vector x -> complex bus voltages V.
## x = [θ(non-slack); Vm(all buses)] with θ_slack fixed to 0.
## A trailing PMU reference-offset state α (if present) is ignored here and
## extracted separately via _va_offset_from_state.
function _state_to_voltage(x::Vector{Float64}, slackIdx::Int, nbus::Int)
  nθ = nbus - 1
  θ = zeros(Float64, nbus)
  vm = Vector{Float64}(undef, nbus)

  for i = 1:nbus
    vm[i] = x[nθ+i]
  end

  p = 0
  for i = 1:nbus
    if i == slackIdx
      θ[i] = 0.0
    else
      p += 1
      θ[i] = x[p]
    end
  end

  return vm .* cis.(θ)
end

## Internal helper: extract the PMU reference-offset state α (radians) from x.
@inline function _va_offset_from_state(x::Vector{Float64}, withVaOffset::Bool)
  return withVaOffset ? x[end] : 0.0
end

## Internal helper: build initial state vector for WLS iterations.
## flatstart=true uses θ=0, Vm=1 except slack magnitude initialized from net.
## withVaOffset=true appends the PMU reference-offset state α (radians),
## initialized to 0.
function _initial_state_vector(net::Net, slackIdx::Int; flatstart::Bool = true, withVaOffset::Bool = false)
  nbus = length(net.nodeVec)
  nθ = nbus - 1
  x = zeros(Float64, nθ + nbus + (withVaOffset ? 1 : 0))

  p = 0
  for i = 1:nbus
    if i != slackIdx
      p += 1
      x[p] = flatstart ? 0.0 : deg2rad(something(net.nodeVec[i]._va_deg, 0.0))
    end
  end

  for i = 1:nbus
    vm0 = something(net.nodeVec[i]._vm_pu, 1.0)
    x[nθ+i] = flatstart ? 1.0 : vm0
  end

  if flatstart
    x[nθ+slackIdx] = something(net.nodeVec[slackIdx]._vm_pu, 1.0)
  end

  return x
end

## Internal helper: return only currently active measurements.
function _active_measurements(measurements::Vector{Measurement})
  return [m for m in measurements if m.active]
end

## Internal helper: return active measurements and their original indices.
function _active_measurements_with_indices(measurements::Vector{Measurement})
  idx = Int[]
  active = Measurement[]
  for (k, m) in enumerate(measurements)
    if m.active
      push!(idx, k)
      push!(active, m)
    end
  end
  return active, idx
end

@inline function _set_measurement_active(m::Measurement, active::Bool)
  return Measurement(typ = m.typ, value = m.value, sigma = m.sigma, active = active, busIdx = m.busIdx, branchIdx = m.branchIdx, direction = m.direction, id = m.id)
end

## Internal helper: assemble measurement value vector z from active measurements.
function _measurement_vector(measurements::Vector{Measurement})
  z = Vector{Float64}(undef, length(measurements))
  @inbounds for (i, m) in enumerate(measurements)
    z[i] = m.value
  end
  return z
end

## Internal helper: compute normalized residuals and row-level diagnostics.
function _residual_diagnostics(H::Matrix{Float64}, r::Vector{Float64}, w::Vector{Float64})
  m, _ = size(H)
  G = (H' .* reshape(w, 1, :)) * H
  GI = pinv(G)
  WInv = Diagonal(1.0 ./ w)
  Ω = Symmetric(WInv - H * GI * H')

  rn = zeros(Float64, m)
  for i in eachindex(r)
    σr2 = max(real(Ω[i, i]), 0.0)
    if σr2 <= eps(Float64)
      rn[i] = 0.0
    else
      rn[i] = r[i] / sqrt(σr2)
    end
  end
  return rn
end

@inline function _chi_square_zscore(j::Float64, ν::Int)
  ν <= 0 && return Inf
  return (j - ν) / sqrt(2.0 * ν)
end

function _build_measurement_suspicion_report(activeMeas::Vector{Measurement}, activeIdx::Vector{Int}, r::Vector{Float64}, rn::Vector{Float64}; normalizedThreshold::Float64 = 3.0)
  rows = NamedTuple[]
  for i in eachindex(activeMeas)
    m = activeMeas[i]
    push!(rows, (local_index = i, measurement_index = activeIdx[i], id = m.id, typ = Symbol(m.typ), residual = r[i], normalized_residual = rn[i], abs_normalized_residual = abs(rn[i]), suspicious = abs(rn[i]) >= normalizedThreshold))
  end

  sort!(rows; by = x -> x.abs_normalized_residual, rev = true)
  return rows
end

## Internal helper: assemble diagonal weight entries w = 1/σ² from measurements.
function _weight_vector(measurements::Vector{Measurement})
  w = Vector{Float64}(undef, length(measurements))
  @inbounds for (i, m) in enumerate(measurements)
    w[i] = m.weight
  end
  return w
end

## Internal helper: evaluate nonlinear measurement model h(x).
## Uses nodal complex injections S = V .* conj(Ybus*V) and branch formulas in _measurement_prediction.
function _predict_measurements(measurements::Vector{Measurement}, net::Net, V::Vector{ComplexF64}, Ybus::AbstractMatrix{ComplexF64}; vaOffsetRad::Float64 = 0.0)
  Sbus_MVA = calc_injections(Ybus, V) .* net.baseMVA
  h = Vector{Float64}(undef, length(measurements))
  @inbounds for (i, m) in enumerate(measurements)
    h[i] = _measurement_prediction(m, net, V, Sbus_MVA, vaOffsetRad)
  end
  return h
end

## Internal helper: finite-difference Jacobian H = ∂h/∂x.
##
## Inputs and meaning:
## - measurements: active measurement set defining h(x)
## - x: current state vector x = [θ(non-slack); Vm(all)] plus, when
##   withVaOffset=true, a trailing PMU reference-offset state α (radians)
## - slackIdx, nbus: mapping context for x -> complex voltages V
## - Ybus: network admittance used inside measurement predictions
## - eps: perturbation size ε for forward-difference sensitivities
##
## Numerical formula (column-wise perturbation):
##   H[i,k] ≈ (h_i(x + εe_k) - h_i(x)) / ε
## where e_k is the k-th unit vector.
##
## Returns:
## - H: Jacobian matrix of size (m_measurements × n_states)
## - h0: base prediction vector h(x), reused by caller to avoid recomputation.
function _measurement_jacobian_fd(measurements::Vector{Measurement}, net::Net, x::Vector{Float64}, slackIdx::Int, nbus::Int, Ybus::AbstractMatrix{ComplexF64}; eps::Float64 = 1e-6, withVaOffset::Bool = false)
  h0 = _predict_measurements(measurements, net, _state_to_voltage(x, slackIdx, nbus), Ybus; vaOffsetRad = _va_offset_from_state(x, withVaOffset))
  m = length(measurements)
  n = length(x)
  H = zeros(Float64, m, n)

  for k = 1:n
    xk = copy(x)
    xk[k] += eps
    hk = _predict_measurements(measurements, net, _state_to_voltage(xk, slackIdx, nbus), Ybus; vaOffsetRad = _va_offset_from_state(xk, withVaOffset))
    @inbounds for i = 1:m
      H[i, k] = (hk[i] - h0[i]) / eps
    end
  end

  return H, h0
end

## Internal helper: create row-wise bipartite adjacency from Jacobian sparsity.
function _adjacency_from_sparsity(H::AbstractMatrix{<:Real})
  m, n = size(H)
  adj = [Int[] for _ = 1:m]
  for i = 1:m
    for j = 1:n
      H[i, j] != 0.0 && push!(adj[i], j)
    end
  end
  return adj, n
end

## Internal helper: maximum bipartite matching via Hopcroft-Karp.
##
## What this does here:
## - Left side nodes = measurement rows i of Jacobian H.
## - Right side nodes = state columns j of Jacobian H.
## - Edge (i,j) exists if H[i,j] != 0 (measurement i depends on state j).
##
## Why we need it:
## - Structural observability is checked without numeric values, only sparsity pattern.
## - If maximum matching size equals number of state columns, every state can be
##   structurally "covered" by at least one independent measurement relation.
##
## Algorithm note:
## - Hopcroft-Karp alternates BFS (build layers) and DFS (augmenting paths),
##   giving O(E*sqrt(V)) complexity and robust performance for these checks.
function _hopcroft_karp(adj::Vector{Vector{Int}}, nRight::Int)
  nLeft = length(adj)
  pairU = zeros(Int, nLeft)
  pairV = zeros(Int, nRight)
  dist = fill(-1, nLeft)

  function bfs!()
    q = Int[]
    for u = 1:nLeft
      if pairU[u] == 0
        dist[u] = 0
        push!(q, u)
      else
        dist[u] = -1
      end
    end

    found = false
    head = 1
    while head <= length(q)
      u = q[head]
      head += 1
      for v in adj[u]
        u2 = pairV[v]
        if u2 == 0
          found = true
        elseif dist[u2] == -1
          dist[u2] = dist[u] + 1
          push!(q, u2)
        end
      end
    end
    return found
  end

  function dfs!(u::Int)
    for v in adj[u]
      u2 = pairV[v]
      if u2 == 0 || (dist[u2] == dist[u] + 1 && dfs!(u2))
        pairU[u] = v
        pairV[v] = u
        return true
      end
    end
    dist[u] = -1
    return false
  end

  matching = 0
  while bfs!()
    for u = 1:nLeft
      if pairU[u] == 0 && dfs!(u)
        matching += 1
      end
    end
  end

  return matching
end

## Internal helper: test if row i is numerically redundant.
##
## Idea:
## - Remove measurement row i from Jacobian H.
## - If rank(H_without_i) is still full (= number of state columns), the row is
##   not required for numerical observability and is therefore redundant.
## - If not, row i is a numerical critical measurement.
function _numerical_row_redundant(H::AbstractMatrix{<:Real}, i::Int; tol = nothing)
  m, n = size(H)
  @assert 1 <= i <= m
  m <= 1 && return false
  keep = [k for k = 1:m if k != i]
  return numeric_rank(H[keep, :]; tol = tol) == n
end

## Internal helper: test if row i is structurally redundant.
##
## Idea:
## - Build bipartite graph from sparsity pattern of H.
## - Remove all edges of row i.
## - If maximum matching still covers all state columns, the row is structurally
##   redundant; otherwise it is structurally critical.
function _structural_row_redundant(H::AbstractMatrix{<:Real}, i::Int)
  m, n = size(H)
  @assert 1 <= i <= m
  adj, _ = _adjacency_from_sparsity(H)
  adj[i] = Int[]
  return _hopcroft_karp(adj, n) == n
end

## Internal helper: weighted least-squares objective J = r'Wr = Σ_i w_i r_i².
@inline _wls_objective(r::Vector{Float64}, w::Vector{Float64}) = sum(w .* (r .^ 2))

## Internal helper: χ²-like 3σ plausibility check using normal approximation
## z = (J - ν) / sqrt(2ν), with ν = m - n.
function _j_within_3sigma_band(j::Float64, ν::Int)
  ν <= 0 && return false
  z = (j - ν) / sqrt(2.0 * ν)
  return abs(z) <= 3.0
end

## Internal helper: classify redundancy/observability quality.
function _redundancy_quality(observable::Bool, ν::Int, hasCritical::Bool)
  if !observable
    return :not_observable
  elseif ν <= 0 || hasCritical
    return :critical
  else
    return :good
  end
end

## Internal helper: compute numerical + structural observability metrics from Jacobian.
##
## Computed quantities:
## - Numerical observability via matrix rank(H)
## - Structural observability via maximum bipartite matching
## - Redundancy metrics ν = m - n and ρ = m/n
## - Critical measurements by single-row removal tests (numerical + structural)
##
## activeOriginalIdx maps local Jacobian rows back to original measurement indices.
function _evaluate_observability_from_jacobian(H::Matrix{Float64}, activeOriginalIdx::Vector{Int}; tol = nothing)
  m, n = size(H)
  ν = m - n
  ρ = n > 0 ? m / n : Inf

  nrank = numeric_rank(H; tol = tol)
  adj, ncols = _adjacency_from_sparsity(H)
  mm = _hopcroft_karp(adj, ncols)

  numObs = nrank == ncols
  strObs = mm == ncols

  criticalNum = Int[]
  criticalStr = Int[]
  if m > 0
    for i = 1:m
      _numerical_row_redundant(H, i; tol = tol) || push!(criticalNum, activeOriginalIdx[i])
      _structural_row_redundant(H, i) || push!(criticalStr, activeOriginalIdx[i])
    end
  end

  hasCritical = !isempty(criticalNum) || !isempty(criticalStr)

  # Unobservable state columns (the "dark" states): computed only on the
  # not-observable path so observable sets pay nothing extra. The null
  # space of H spans exactly the state directions no measurement sees;
  # a state column with a nonzero component in ANY basis vector cannot be
  # estimated (its value can drift inside the null space without changing
  # a single measurement). `nullspace` returns an ORTHONORMAL basis, so a
  # small absolute threshold on the entries is meaningful. Dense SVD cost:
  # acceptable at workshop and distribution-network sizes (a 2700 x 2700
  # null-space probe is seconds), and it only runs on failures.
  dark = Int[]
  if !numObs && n > 0 && m > 0
    ns = nullspace(H)
    for j = 1:n
      any(k -> abs(ns[j, k]) > 1e-8, axes(ns, 2)) && push!(dark, j)
    end
  end

  return (
    numerical_observable = numObs,
    structural_observable = strObs,
    numerical_rank = nrank,
    structural_matching = mm,
    n_states = ncols,
    n_measurements = m,
    redundancy = ν,
    redundancy_ratio = ρ,
    dof = ν,
    numerical_critical_measurement_indices = criticalNum,
    structural_critical_measurement_indices = criticalStr,
    unobservable_state_columns = dark,
    quality = _redundancy_quality(numObs && strObs, ν, hasCritical),
  )
end

"""
    numerical_observable(H; tol=nothing) -> Bool

Numerical observability test on a Jacobian-like matrix `H`.
Returns `true` when `rank(H) == n` (full column rank).
"""
function numerical_observable(H::AbstractMatrix{<:Real}; tol = nothing)
  _, n = size(H)
  return numeric_rank(H; tol = tol) == n
end

"""
    structural_observable(H) -> Bool

Structural observability test on a Jacobian-like matrix `H`.
Returns `true` when the maximum bipartite matching size equals the number of
state columns `n`.
"""
function structural_observable(H::AbstractMatrix{<:Real})
  _, n = size(H)
  adj, _ = _adjacency_from_sparsity(H)
  return _hopcroft_karp(adj, n) == n
end

"""
    numerical_row_redundant(H, i; tol=nothing) -> Bool

Check if row `i` remains numerically redundant in `H`.
"""
function numerical_row_redundant(H::AbstractMatrix{<:Real}, i::Int; tol = nothing)
  return _numerical_row_redundant(H, i; tol = tol)
end

"""
    structural_row_redundant(H, i) -> Bool

Check if row `i` remains structurally redundant in `H`.
"""
function structural_row_redundant(H::AbstractMatrix{<:Real}, i::Int)
  return _structural_row_redundant(H, i)
end

"""
    evaluate_observability_matrix(H; tol=nothing) -> NamedTuple

Evaluate global observability and single-row criticality directly on a matrix
`H` (without building a network model).

The result includes `unobservable_state_columns::Vector{Int}`: for a
numerically NOT observable `H`, the state columns with a component above
tolerance in any null-space basis vector, i.e. exactly the states no
measurement pins down (their union partitions the network into observable
islands). Empty for observable systems; computed only on the
not-observable path (a dense null-space probe, fine at workshop and
distribution-network sizes).
"""
function evaluate_observability_matrix(H::AbstractMatrix{<:Real}; tol = nothing)
  m, _ = size(H)
  idx = collect(1:m)
  return _evaluate_observability_from_jacobian(Matrix{Float64}(H), idx; tol = tol)
end

"""
    evaluate_local_observability_matrix(H, stateCols; tol=nothing) -> NamedTuple

Evaluate local observability on a matrix `H` restricted to selected
`stateCols`.

Limitation: this column-restricted submatrix test is NECESSARY but not
sufficient; a positive verdict can be wrong when the touching rows couple
the selected states to neighbor states that are themselves undetermined
(the rigorous per-state answer is `unobservable_state_columns` from the
global check).
"""
function evaluate_local_observability_matrix(H::AbstractMatrix{<:Real}, stateCols::Vector{Int}; tol = nothing)
  isempty(stateCols) && error("evaluate_local_observability_matrix: stateCols must not be empty")
  nstates = size(H, 2)
  for c in stateCols
    1 <= c <= nstates || error("evaluate_local_observability_matrix: state column out of bounds: $(c)")
  end

  rows = Int[]
  for i in axes(H, 1)
    for j in stateCols
      if !iszero(H[i, j])
        push!(rows, i)
        break
      end
    end
  end

  isempty(rows) && error("evaluate_local_observability_matrix: no row touches selected stateCols")
  local_idx = [rows[k] for k in eachindex(rows)]
  base = _evaluate_observability_from_jacobian(Matrix{Float64}(H[rows, stateCols]), local_idx; tol = tol)
  return merge(base, (rows = rows, stateCols = copy(stateCols)))
end

"""
    evaluate_global_observability(net, measurements; kwargs...) -> NamedTuple

Evaluate global observability on active measurements using the finite-difference
measurement Jacobian.

Includes global redundancy metrics
- `redundancy = r = m - n`
- `redundancy_ratio = ρ = m / n`
- `dof = ν = m - n` (a COUNT difference: for an observable set it equals
  `m - rank(H)`; for an unobservable set it can be negative and then reads
  as a shortfall, not a redundancy)

Quality classes:
- `:good`: observable and no critical single measurement
- `:critical`: observable, but at least one single critical measurement (or ν <= 0)
- `:not_observable`: not observable

For a not-observable set the result additionally names the dark states in
`unobservable_state_columns` (state columns touched by the null space of
`H`); empty when observable. This is the rigorous per-state answer that
the column-restricted local check cannot give (see
[`evaluate_local_observability`](@ref)).
"""
function evaluate_global_observability(net::Net, measurements::Vector{Measurement}; flatstart::Bool = true, jacEps::Float64 = 1e-6, tol = nothing, pmuRefOffset::Symbol = state_estimation_config().pmu_ref_offset)
  activeMeas, activeIdx = _active_measurements_with_indices(measurements)
  isempty(activeMeas) && error("evaluate_global_observability: no active measurements")

  nbus = length(net.nodeVec)
  slackIdx = _find_slack_idx(net)
  Ybus = createYBUS(net = net)
  withVaOffset = _va_offset_active(activeMeas, pmuRefOffset)
  x = _initial_state_vector(net, slackIdx; flatstart = flatstart, withVaOffset = withVaOffset)
  H, _ = _measurement_jacobian_fd(activeMeas, net, x, slackIdx, nbus, Ybus; eps = jacEps, withVaOffset = withVaOffset)

  return _evaluate_observability_from_jacobian(H, activeIdx; tol = tol)
end

function evaluate_global_observability(net::Net; flatstart::Bool = true, jacEps::Float64 = 1e-6, tol = nothing, pmuRefOffset::Symbol = state_estimation_config().pmu_ref_offset)
  return evaluate_global_observability(net, Measurement[m for m in net.measurements]; flatstart = flatstart, jacEps = jacEps, tol = tol, pmuRefOffset = pmuRefOffset)
end

"""
    measurement_jacobian(net; flatstart=true, jacEps=1e-6, pmuRefOffset=...) -> NamedTuple

Build the measurement Jacobian `H` of the ACTIVE measurements on `net`,
labeled for humans: the matrix the observability checks and the WLS
normal equations run on, with one described row per measurement and one
described column per state.

Returns `(H, rows, cols)`:
- `H::Matrix{Float64}`: m x n finite-difference Jacobian at the flat (or
  stored) start state, the same evaluation `evaluate_global_observability`
  uses.
- `rows`: one NamedTuple per active measurement,
  `(index, type, location, sigma)`; `index` is the position in
  `net.measurements`, `location` names the bus (injections, voltages) or
  the oriented branch (flows).
- `cols`: state-column labels in Jacobian order, `"Va(bus)"` for every
  non-slack bus, then `"Vm(bus)"` for every bus, plus `"alpha"` when PMU
  `Va` measurements activate the reference-offset state.

Errors when no active measurement exists. Intended for measurement-matrix
reports and placement studies; see the state-estimation suite summary and
the workshop's observability deep dive.
"""
function measurement_jacobian(net::Net; flatstart::Bool = true, jacEps::Float64 = 1e-6, pmuRefOffset::Symbol = state_estimation_config().pmu_ref_offset)
  activeMeas, activeIdx = _active_measurements_with_indices(Measurement[m for m in net.measurements])
  isempty(activeMeas) && error("measurement_jacobian: no active measurements")
  nbus = length(net.nodeVec)
  slackIdx = _find_slack_idx(net)
  Ybus = createYBUS(net = net)
  withVaOffset = _va_offset_active(activeMeas, pmuRefOffset)
  x = _initial_state_vector(net, slackIdx; flatstart = flatstart, withVaOffset = withVaOffset)
  H, _ = _measurement_jacobian_fd(activeMeas, net, x, slackIdx, nbus, Ybus; eps = jacEps, withVaOffset = withVaOffset)

  name_by_idx = _bus_name_by_idx(net)
  busname(i) = get(name_by_idx, i, string(i))
  rows = NamedTuple[]
  for (k, m) in enumerate(activeMeas)
    location = if m.typ in (PflowMeas, QflowMeas)
      br = net.branchVec[m.branchIdx]
      ends = (busname(Int(br.fromBus)), busname(Int(br.toBus)))
      m.direction === :to ? string(ends[2], "->", ends[1]) : string(ends[1], "->", ends[2])
    else
      busname(m.busIdx)
    end
    push!(rows, (index = activeIdx[k], type = m.typ, location = location, sigma = m.sigma))
  end
  cols = String[]
  for i = 1:nbus
    i == slackIdx || push!(cols, string("Va(", busname(i), ")"))
  end
  for i = 1:nbus
    push!(cols, string("Vm(", busname(i), ")"))
  end
  withVaOffset && push!(cols, "alpha")
  return (H = H, rows = rows, cols = cols)
end

"""
    evaluate_local_observability(net, measurements, stateCols; kwargs...) -> NamedTuple

Evaluate local observability on selected Jacobian columns (`stateCols`).

Procedure:
1) Build global Jacobian `H` from currently active measurements.
2) Keep only rows that have at least one nonzero entry in the selected columns.
   These rows correspond to measurements that are locally sensitive to the
   requested states.
3) Evaluate observability/redundancy on the reduced matrix `Hlocal`.

Returned NamedTuple extends global metrics with:
- `rows`: selected row indices (within global active-Jacobian row numbering)
- `stateCols`: copied input state-column selection.

Interpretation:
- `:good` means local states are observable with positive redundancy and no
  single critical measurement.
- `:critical` means still observable but vulnerable to a single outage (or ν <= 0).
- `:not_observable` means local states cannot be uniquely reconstructed.

Limitation: this column-restricted test is NECESSARY but not sufficient; a
positive verdict can be wrong when the touching rows couple the selected
states to neighbor states that are themselves undetermined. The rigorous
per-state answer is `unobservable_state_columns` from
[`evaluate_global_observability`](@ref).
"""
function evaluate_local_observability(net::Net, measurements::Vector{Measurement}, stateCols::Vector{Int}; flatstart::Bool = true, jacEps::Float64 = 1e-6, tol = nothing, pmuRefOffset::Symbol = state_estimation_config().pmu_ref_offset)
  isempty(stateCols) && error("evaluate_local_observability: stateCols must not be empty")

  activeMeas, activeIdx = _active_measurements_with_indices(measurements)
  isempty(activeMeas) && error("evaluate_local_observability: no active measurements")

  nbus = length(net.nodeVec)
  slackIdx = _find_slack_idx(net)
  Ybus = createYBUS(net = net)
  withVaOffset = _va_offset_active(activeMeas, pmuRefOffset)
  x = _initial_state_vector(net, slackIdx; flatstart = flatstart, withVaOffset = withVaOffset)
  H, _ = _measurement_jacobian_fd(activeMeas, net, x, slackIdx, nbus, Ybus; eps = jacEps, withVaOffset = withVaOffset)

  nstates = size(H, 2)
  for c in stateCols
    1 <= c <= nstates || error("evaluate_local_observability: state column out of bounds: $(c)")
  end

  localRows = Int[]

  for i in axes(H, 1)
    for j in stateCols
      if !iszero(H[i, j])
        push!(localRows, i)
        break
      end
    end
  end

  isempty(localRows) && error("evaluate_local_observability: no active measurement touches selected stateCols")

  Hlocal = H[localRows, stateCols]
  localOriginalIdx = [activeIdx[i] for i in localRows]
  base = _evaluate_observability_from_jacobian(Hlocal, localOriginalIdx; tol = tol)

  return merge(base, (rows = localRows, stateCols = copy(stateCols)))
end

function evaluate_local_observability(net::Net, stateCols::Vector{Int}; flatstart::Bool = true, jacEps::Float64 = 1e-6, tol = nothing, pmuRefOffset::Symbol = state_estimation_config().pmu_ref_offset)
  return evaluate_local_observability(net, Measurement[m for m in net.measurements], stateCols; flatstart = flatstart, jacEps = jacEps, tol = tol, pmuRefOffset = pmuRefOffset)
end

"""
    runse!(net, measurements; kwargs...) -> SEResult

Run a first classical nonlinear weighted least-squares state estimator.

State representation:
- bus voltage angles for all non-slack buses (radians)
- bus voltage magnitudes for all buses (p.u.)
- optional PMU reference-angle offset α (radians), appended automatically
  when active `VaMeas` measurements exist and `pmu_ref_offset = :auto`

PMU angle measurements (`VaMeas`, degrees) are modeled as
`z = θ_i + α + e`: the network angles stay slack-referenced, α maps them
into the common PMU time base. With `pmu_ref_offset = :off` the offset
state is omitted and PMU angles are assumed to be slack-referenced already.
"""
function runse!(net::Net, measurements::Vector{Measurement}, cfg::StateEstimationConfig)
  return _runse_with_config!(net, measurements, cfg)
end

function _runse_with_config!(net::Net, measurements::Vector{Measurement}, cfg::StateEstimationConfig)
  maxIte = cfg.max_iter
  tol = cfg.tol
  flatstart = cfg.flatstart
  jacEps = cfg.jac_eps
  updateNet = cfg.update_net
  activeMeas = _active_measurements(measurements)
  isempty(activeMeas) && error("runse!: no active measurements")
  # r0.9.10 scope guard: the estimator's measurement model does not carry
  # the one-sided branch reduction yet; a flow measurement on a partially
  # open branch would be evaluated against the wrong model. Rejected with a
  # clear message (follow-up tracked in the terminal-status review notes);
  # measurement generation never creates such measurements.
  for meas in activeMeas
    meas.branchIdx === nothing && continue
    (1 <= meas.branchIdx <= length(net.branchVec)) || continue
    st = _branch_terminal_state(net.branchVec[meas.branchIdx])
    if st == :open_from || st == :open_to
      error("runse!: measurement $(meas.id) references branch $(meas.branchIdx), which is open at one terminal ($(st)). Flow measurements on partially open branches are not supported yet; remove the measurement or close the terminal.")
    end
  end

  nbus = length(net.nodeVec)
  slackIdx = _find_slack_idx(net)
  Ybus = createYBUS(net = net)

  withVaOffset = _va_offset_active(activeMeas, cfg.pmu_ref_offset)
  x = _initial_state_vector(net, slackIdx; flatstart = flatstart, withVaOffset = withVaOffset)
  z = _measurement_vector(activeMeas)
  w = _weight_vector(activeMeas)

  converged = false
  r = zeros(Float64, length(activeMeas))
  iteDone = 0

  for ite = 1:maxIte
    H, h = _measurement_jacobian_fd(activeMeas, net, x, slackIdx, nbus, Ybus; eps = jacEps, withVaOffset = withVaOffset)
    r = z - h

    # Normal equations for WLS:
    #   Δx = (H'WH)^{-1} H'W r
    # with W = diag(w).
    HW = H' .* reshape(w, 1, :)
    G = HW * H
    g = HW * r

    Δx = solve_linear(G, g; allow_pinv = true)
    x .+= Δx

    iteDone = ite
    if norm(Δx, Inf) < tol
      converged = true
      break
    end
  end

  Vest = _state_to_voltage(x, slackIdx, nbus)
  vaOffsetRad = _va_offset_from_state(x, withVaOffset)
  hfinal = _predict_measurements(activeMeas, net, Vest, Ybus; vaOffsetRad = vaOffsetRad)
  r = z - hfinal
  jval = _wls_objective(r, w)
  ν = length(activeMeas) - length(x)

  if updateNet
    for i = 1:nbus
      net.nodeVec[i]._vm_pu = abs(Vest[i])
      if i == slackIdx
        net.nodeVec[i]._va_deg = 0.0
      else
        net.nodeVec[i]._va_deg = rad2deg(angle(Vest[i]))
      end
    end
    calcNetLosses!(net, Vest)
  end

  return SEResult(Vest, converged, iteDone, norm(r), r, jval, ν, _j_within_3sigma_band(jval, ν), withVaOffset ? rad2deg(vaOffsetRad) : nothing)
end

function runse!(
  net::Net,
  measurements::Vector{Measurement};
  maxIte::Int = state_estimation_config().max_iter,
  tol::Float64 = state_estimation_config().tol,
  flatstart::Bool = state_estimation_config().flatstart,
  jacEps::Float64 = state_estimation_config().jac_eps,
  updateNet::Bool = state_estimation_config().update_net,
  pmuRefOffset::Symbol = state_estimation_config().pmu_ref_offset,
)
  cfg = StateEstimationConfig(max_iter = maxIte, tol = tol, flatstart = flatstart, jac_eps = jacEps, update_net = updateNet, pmu_ref_offset = pmuRefOffset)
  return _runse_with_config!(net, measurements, cfg)
end

function runse!(
  net::Net;
  maxIte::Int = state_estimation_config().max_iter,
  tol::Float64 = state_estimation_config().tol,
  flatstart::Bool = state_estimation_config().flatstart,
  jacEps::Float64 = state_estimation_config().jac_eps,
  updateNet::Bool = state_estimation_config().update_net,
  pmuRefOffset::Symbol = state_estimation_config().pmu_ref_offset,
)
  cfg = StateEstimationConfig(max_iter = maxIte, tol = tol, flatstart = flatstart, jac_eps = jacEps, update_net = updateNet, pmu_ref_offset = pmuRefOffset)
  return runse!(net, Measurement[m for m in net.measurements], cfg)
end

"""
    validate_measurements(net, measurements; kwargs...) -> NamedTuple

Run state-estimation diagnostics on currently active measurements and return a
machine-readable report with:
- global bad-data consistency check (`global_consistency`)
- χ²-like objective plausibility summary
- largest-normalized-residual ranking
- suspicious measurement list (threshold-based)
"""
function validate_measurements(net::Net, measurements::Vector{Measurement}; maxIte::Int = 12, tol::Float64 = 1e-6, flatstart::Bool = true, jacEps::Float64 = 1e-6, normalizedThreshold::Float64 = 3.0, pmuRefOffset::Symbol = state_estimation_config().pmu_ref_offset)
  activeMeas, activeIdx = _active_measurements_with_indices(measurements)
  isempty(activeMeas) && error("validate_measurements: no active measurements")

  netLocal = deepcopy(net)
  result = runse!(netLocal, activeMeas; maxIte = maxIte, tol = tol, flatstart = flatstart, jacEps = jacEps, updateNet = false, pmuRefOffset = pmuRefOffset)

  nbus = length(netLocal.nodeVec)
  slackIdx = _find_slack_idx(netLocal)
  withVaOffset = result.vaRefOffsetDeg !== nothing
  x = _initial_state_vector(netLocal, slackIdx; flatstart = flatstart, withVaOffset = withVaOffset)
  nθ = nbus - 1
  p = 0
  for i in eachindex(netLocal.nodeVec)
    if i != slackIdx
      p += 1
      x[p] = angle(result.voltages[i])
    end
    x[nθ+i] = abs(result.voltages[i])
  end
  if withVaOffset
    x[end] = deg2rad(result.vaRefOffsetDeg)
  end

  Ybus = createYBUS(net = netLocal)
  H, h = _measurement_jacobian_fd(activeMeas, netLocal, x, slackIdx, nbus, Ybus; eps = jacEps, withVaOffset = withVaOffset)
  z = _measurement_vector(activeMeas)
  w = _weight_vector(activeMeas)
  r = z - h
  rn = _residual_diagnostics(H, r, w)
  ranking = _build_measurement_suspicion_report(activeMeas, activeIdx, r, rn; normalizedThreshold = normalizedThreshold)

  suspicious = [row for row in ranking if row.suspicious]
  zscore = _chi_square_zscore(result.objectiveJ, result.dof)

  return (
    converged = result.converged,
    global_consistency = result.converged && result.jWithin3Sigma,
    objective = (value = result.objectiveJ, dof = result.dof, zscore = zscore, within_3sigma = result.jWithin3Sigma),
    largest_normalized_residual = isempty(ranking) ? nothing : first(ranking),
    suspicious_measurements = suspicious,
    measurement_ranking = ranking,
    residuals = r,
    normalized_residuals = rn,
    result = result,
  )
end

function validate_measurements(net::Net; kwargs...)
  return validate_measurements(net, Measurement[m for m in net.measurements]; kwargs...)
end

"""
    runse_diagnostics(net, measurements; deactivate_and_rerun=false, kwargs...) -> NamedTuple

Extended diagnostics workflow around `validate_measurements` with optional
deactivate-and-rerun logic for the currently largest suspicious measurement.
"""
function runse_diagnostics(net::Net, measurements::Vector{Measurement}; deactivate_and_rerun::Bool = false, maxIte::Int = 12, tol::Float64 = 1e-6, flatstart::Bool = true, jacEps::Float64 = 1e-6, normalizedThreshold::Float64 = 3.0, pmuRefOffset::Symbol = state_estimation_config().pmu_ref_offset)
  base = validate_measurements(net, measurements; maxIte = maxIte, tol = tol, flatstart = flatstart, jacEps = jacEps, normalizedThreshold = normalizedThreshold, pmuRefOffset = pmuRefOffset)
  rerun = nothing

  if deactivate_and_rerun && !isempty(base.suspicious_measurements)
    target = first(base.suspicious_measurements)
    meas2 = copy(measurements)
    idx = target.measurement_index
    meas2[idx] = _set_measurement_active(meas2[idx], false)
    rerun = validate_measurements(net, meas2; maxIte = maxIte, tol = tol, flatstart = flatstart, jacEps = jacEps, normalizedThreshold = normalizedThreshold, pmuRefOffset = pmuRefOffset)
    rerun = (deactivated_measurement_index = idx, deactivated_measurement_id = measurements[idx].id, diagnostics = rerun)
  end

  return (diagnostics = base, rerun = rerun)
end

function runse_diagnostics(net::Net; kwargs...)
  return runse_diagnostics(net, Measurement[m for m in net.measurements]; kwargs...)
end

@inline function _base_diag_report(diag)
  return hasproperty(diag, :diagnostics) ? diag.diagnostics : diag
end

"""
    summarize_se_diagnostics(diag) -> NamedTuple

Create a compact interpretation summary for a diagnostics object returned by
`validate_measurements` or `runse_diagnostics`.

`global_consistency` is interpreted as:
- `true`: SE converged and objective is inside χ²-like 3σ plausibility band
- `false`: either non-convergence or implausibly large objective
"""
function summarize_se_diagnostics(diag)
  base = _base_diag_report(diag)
  suspicious_count = length(base.suspicious_measurements)
  total = length(base.measurement_ranking)

  reason = if !base.converged
    "State estimation did not converge."
  elseif !base.objective.within_3sigma
    "Objective J is outside the χ²-like 3σ plausibility band (possible bad data/model mismatch)."
  else
    "No global inconsistency detected."
  end

  return (global_consistency = base.global_consistency, converged = base.converged, objective = base.objective, suspicious_count = suspicious_count, total_measurements = total, reason = reason)
end

@inline function _global_consistency_label(summary)
  return summary.global_consistency ? "PASS (globally plausible)" : "FAIL (check objective/residuals)"
end

"""
    print_se_diagnostics(io, diag; topN=10)    

Pretty-print diagnostics from `validate_measurements` or `runse_diagnostics`
including:
- explanation of `global_consistency`
- tabular measurement ranking
- BAD/OK marker per measurement
- optional rerun comparison if present
"""
function print_se_diagnostics(diag; io = stdout, topN::Int = 10, format::Symbol = :plain)
  base = _base_diag_report(diag)
  summary = summarize_se_diagnostics(diag)

  nrows = min(topN, length(base.measurement_ranking))
  if format == :markdown
    println(io, "## State-estimation diagnostics")
    println(io)
    println(io, "- **Converged:** $(summary.converged)")
    println(io, "- **Global consistency:** $(summary.global_consistency) — $(_global_consistency_label(summary))")
    @printf(io, "- **Objective J:** %.6f (dof=%d, z=%.3f, within_3sigma=%s)\n", summary.objective.value, summary.objective.dof, summary.objective.zscore, string(summary.objective.within_3sigma))
    println(io, "- **Interpretation:** $(summary.reason)")
    println(io, "- **Suspicious measurements:** $(summary.suspicious_count) / $(summary.total_measurements)")
    println(io)
    println(io, "### Measurement ranking (largest |normalized residual|)")
    println(io, "| Idx | ID | Type | Residual | Norm.Res | Flag |")
    println(io, "|---:|:---|:---|---:|---:|:---:|")
    for k = 1:nrows
      row = base.measurement_ranking[k]
      flag = row.suspicious ? "BAD" : "OK"
      @printf(io, "| %d | %s | %s | %.5f | %.5f | %s |\n", row.measurement_index, row.id, string(row.typ), row.residual, row.normalized_residual, flag)
    end
  else
    println(io, "State-estimation diagnostics")
    println(io, "--------------------------------------------------------------")
    @printf(io, "Converged: %s\n", string(summary.converged))
    @printf(io, "Global consistency: %s (%s)\n", string(summary.global_consistency), _global_consistency_label(summary))
    @printf(io, "Objective J: %.6f (dof=%d, z=%.3f, within_3sigma=%s)\n", summary.objective.value, summary.objective.dof, summary.objective.zscore, string(summary.objective.within_3sigma))
    @printf(io, "Interpretation: %s\n", summary.reason)
    @printf(io, "Suspicious measurements: %d / %d\n", summary.suspicious_count, summary.total_measurements)

    println(io, "\nMeasurement ranking (largest |normalized residual|)")
    println(io, "----------------------------------------------------------------------------------------------")
    @printf(io, "%5s %-28s %-12s %12s %12s %8s\n", "Idx", "ID", "Type", "Residual", "Norm.Res", "Flag")
    println(io, "----------------------------------------------------------------------------------------------")
    for k = 1:nrows
      row = base.measurement_ranking[k]
      flag = row.suspicious ? "BAD" : "OK"
      @printf(io, "%5d %-28s %-12s %12.5f %12.5f %8s\n", row.measurement_index, row.id, string(row.typ), row.residual, row.normalized_residual, flag)
    end
  end

  if hasproperty(diag, :rerun) && !isnothing(diag.rerun)
    rerun = diag.rerun
    rerun_summary = summarize_se_diagnostics(rerun.diagnostics)
    if format == :markdown
      println(io)
      println(io, "### Deactivate-and-rerun")
      println(io, "- **Deactivated measurement:** idx=$(rerun.deactivated_measurement_index), id=$(rerun.deactivated_measurement_id)")
      @printf(io, "- **Objective before:** %.6f\n", base.objective.value)
      @printf(io, "- **Objective after:** %.6f\n", rerun.diagnostics.objective.value)
      @printf(io, "- **Objective after stats:** dof=%d, z=%.3f, within_3sigma=%s\n", rerun_summary.objective.dof, rerun_summary.objective.zscore, string(rerun_summary.objective.within_3sigma))
      println(io, "- **Global consistency after rerun:** $(rerun_summary.global_consistency) — $(_global_consistency_label(rerun_summary))")
      println(io, "- **Interpretation after rerun:** $(rerun_summary.reason)")
    else
      println(io, "\nDeactivate-and-rerun")
      println(io, "--------------------------------------------------------------")
      @printf(io, "Deactivated measurement: idx=%d id=%s\n", rerun.deactivated_measurement_index, rerun.deactivated_measurement_id)
      @printf(io, "Objective before: %.6f\n", base.objective.value)
      @printf(io, "Objective after : %.6f\n", rerun.diagnostics.objective.value)
      @printf(io, "Objective after stats: dof=%d, z=%.3f, within_3sigma=%s\n", rerun_summary.objective.dof, rerun_summary.objective.zscore, string(rerun_summary.objective.within_3sigma))
      @printf(io, "Global consistency after rerun: %s (%s)\n", string(rerun_summary.global_consistency), _global_consistency_label(rerun_summary))
      @printf(io, "Interpretation after rerun: %s\n", rerun_summary.reason)
    end
  end
end

# Backward-compatible positional-IO method kept for existing callers/tests.
function print_se_diagnostics(io::IO, diag; topN::Int = 10, format::Symbol = :plain)
  return print_se_diagnostics(diag; io = io, topN = topN, format = format)
end
