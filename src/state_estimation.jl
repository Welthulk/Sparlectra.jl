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

# file: src/state_estimation.jl

"""
    SEResult

Result container for the first classical WLS state-estimation run.
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

## Internal helper: map estimator state vector x -> complex bus voltages V.
## x = [θ(non-slack); Vm(all buses)] with θ_slack fixed to 0.
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

## Internal helper: build initial state vector for WLS iterations.
## flatstart=true uses θ=0, Vm=1 except slack magnitude initialized from net.
function _initial_state_vector(net::Net, slackIdx::Int; flatstart::Bool = true)
  nbus = length(net.nodeVec)
  nθ = nbus - 1
  x = zeros(Float64, nθ + nbus)

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

## Internal helper: assemble measurement value vector z from active measurements.
function _measurement_vector(measurements::Vector{Measurement})
  z = Vector{Float64}(undef, length(measurements))
  @inbounds for (i, m) in enumerate(measurements)
    z[i] = m.value
  end
  return z
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
function _predict_measurements(measurements::Vector{Measurement}, net::Net, V::Vector{ComplexF64}, Ybus::AbstractMatrix{ComplexF64})
  Sbus_MVA = calc_injections(Ybus, V) .* net.baseMVA
  h = Vector{Float64}(undef, length(measurements))
  @inbounds for (i, m) in enumerate(measurements)
    h[i] = _measurement_prediction(m, net, V, Sbus_MVA)
  end
  return h
end

## Internal helper: finite-difference Jacobian H = ∂h/∂x.
##
## Inputs and meaning:
## - measurements: active measurement set defining h(x)
## - x: current state vector x = [θ(non-slack); Vm(all)]
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
function _measurement_jacobian_fd(measurements::Vector{Measurement}, net::Net, x::Vector{Float64}, slackIdx::Int, nbus::Int, Ybus::AbstractMatrix{ComplexF64}; eps::Float64 = 1e-6)
  h0 = _predict_measurements(measurements, net, _state_to_voltage(x, slackIdx, nbus), Ybus)
  m = length(measurements)
  n = length(x)
  H = zeros(Float64, m, n)

  for k = 1:n
    xk = copy(x)
    xk[k] += eps
    hk = _predict_measurements(measurements, net, _state_to_voltage(xk, slackIdx, nbus), Ybus)
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
    quality = _redundancy_quality(numObs && strObs, ν, hasCritical),
  )
end

"""
    evaluate_global_observability(net, measurements; kwargs...) -> NamedTuple

Evaluate global observability on active measurements using the finite-difference
measurement Jacobian.

Includes global redundancy metrics
- `redundancy = r = m - n`
- `redundancy_ratio = ρ = m / n`
- `dof = ν = m - n`

Quality classes:
- `:good`: observable and no critical single measurement
- `:critical`: observable, but at least one single critical measurement (or ν <= 0)
- `:not_observable`: not observable
"""
function evaluate_global_observability(net::Net, measurements::Vector{Measurement}; flatstart::Bool = true, jacEps::Float64 = 1e-6, tol = nothing)
  activeMeas, activeIdx = _active_measurements_with_indices(measurements)
  isempty(activeMeas) && error("evaluate_global_observability: no active measurements")

  nbus = length(net.nodeVec)
  slackIdx = _find_slack_idx(net)
  Ybus = createYBUS(net = net)
  x = _initial_state_vector(net, slackIdx; flatstart = flatstart)
  H, _ = _measurement_jacobian_fd(activeMeas, net, x, slackIdx, nbus, Ybus; eps = jacEps)

  return _evaluate_observability_from_jacobian(H, activeIdx; tol = tol)
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
"""
function evaluate_local_observability(net::Net, measurements::Vector{Measurement}, stateCols::Vector{Int}; flatstart::Bool = true, jacEps::Float64 = 1e-6, tol = nothing)
  isempty(stateCols) && error("evaluate_local_observability: stateCols must not be empty")

  activeMeas, activeIdx = _active_measurements_with_indices(measurements)
  isempty(activeMeas) && error("evaluate_local_observability: no active measurements")

  nbus = length(net.nodeVec)
  slackIdx = _find_slack_idx(net)
  Ybus = createYBUS(net = net)
  x = _initial_state_vector(net, slackIdx; flatstart = flatstart)
  H, _ = _measurement_jacobian_fd(activeMeas, net, x, slackIdx, nbus, Ybus; eps = jacEps)

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

"""
    runse!(net, measurements; kwargs...) -> SEResult

Run a first classical nonlinear weighted least-squares state estimator.

State representation:
- bus voltage angles for all non-slack buses (radians)
- bus voltage magnitudes for all buses (p.u.)
"""
function runse!(net::Net, measurements::Vector{Measurement}; maxIte::Int = 12, tol::Float64 = 1e-6, flatstart::Bool = true, jacEps::Float64 = 1e-6, updateNet::Bool = true)
  activeMeas = _active_measurements(measurements)
  isempty(activeMeas) && error("runse!: no active measurements")

  nbus = length(net.nodeVec)
  slackIdx = _find_slack_idx(net)
  Ybus = createYBUS(net = net)

  x = _initial_state_vector(net, slackIdx; flatstart = flatstart)
  z = _measurement_vector(activeMeas)
  w = _weight_vector(activeMeas)

  converged = false
  r = zeros(Float64, length(activeMeas))
  iteDone = 0

  for ite = 1:maxIte
    H, h = _measurement_jacobian_fd(activeMeas, net, x, slackIdx, nbus, Ybus; eps = jacEps)
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
  hfinal = _predict_measurements(activeMeas, net, Vest, Ybus)
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

  return SEResult(Vest, converged, iteDone, norm(r), r, jval, ν, _j_within_3sigma_band(jval, ν))
end
