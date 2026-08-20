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

# file: src/condition_number.jl
# purpose: condition-number diagnostics for Newton-Raphson Jacobians.
#          Hager/Higham 1-norm estimator on the LU factorization (sparse
#          friendly) with an exact dense SVD path for small systems; no
#          dependencies beyond LinearAlgebra and SparseArrays.

"""
    condestJacobian(J::AbstractMatrix{<:Real}; exact::Bool = false) -> Float64

Estimate the 1-norm condition number ``\\kappa_1(J) = \\|J\\|_1 \\cdot \\|J^{-1}\\|_1``
of a power-flow Jacobian.

- default: Hager/Higham estimator for ``\\|J^{-1}\\|_1`` on the LU
  factorization, so it works on the sparse Jacobians the Newton step
  factors anyway. The estimate matches the exact condition number's order
  of magnitude; a factor of 2 to 3 spread is normal.
- `exact = true`: exact 2-norm condition number via dense SVD. This
  materializes the full matrix, so use it only for small systems.

The matrix must be real (the rectangular Newton-Raphson Jacobian is; the
Hager sign-vector test below is not defined for complex entries) and
square. Throws an `ArgumentError` otherwise, and rethrows the factorization
error for a structurally singular matrix.
"""
function condestJacobian(J::AbstractMatrix{<:Real}; exact::Bool = false)
  n, m = size(J)
  n == m || throw(ArgumentError("condestJacobian needs a square matrix, got $(n)x$(m)"))
  n == 0 && throw(ArgumentError("condestJacobian needs a non-empty matrix"))
  exact && return cond(Matrix(J), 2)

  F = lu(J)            # the Newton step needs this factorization anyway
  normJ = opnorm(J, 1) # exact 1-norm (maximum absolute column sum)

  # Hager iteration for the inverse 1-norm: alternate solves with J and J'
  # on sign vectors, following the gradient of x -> norm(J \ x, 1) on the
  # unit simplex. Typically converges within 2 to 3 rounds; 5 is a hard cap.
  x = fill(1.0 / n, n)
  gamma = 0.0
  for _ = 1:5
    y = F \ x
    gamma = norm(y, 1)
    xi = sign.(y)
    z = F' \ xi
    j = argmax(abs.(z))
    # no improving vertex left: the current estimate is (locally) maximal
    abs(z[j]) <= dot(z, x) && break
    fill!(x, 0.0)
    x[j] = 1.0
  end
  return normJ * gamma
end

# Reduce full-bus-order condest inputs to the active (non-isolated)
# subsystem. De-energized or link-merged buses carry zero Ybus rows and
# columns and no Jacobian entries, so the FULL Jacobian is structurally
# singular whenever such buses exist; the operative system the Newton step
# effectively solves is the active-bus submatrix. `Ybus` may arrive in full
# dimension (solver path, expanded with zero rows/cols) or already reduced
# to active-bus order (createYBUS output); both are accepted and checked.
# Returns `(Ybus, V, bus_types, Vset, slack_idx)` in active-bus order.
function _active_condest_inputs(Ybus, V::Vector{ComplexF64}, bus_types::Vector{Symbol}, Vset::Vector{Float64}, slack_idx::Int, iso_nodes::Vector{Int})
  n = length(V)
  iso = [b for b in iso_nodes if 1 <= b <= n]
  if isempty(iso)
    size(Ybus, 1) == n || throw(ArgumentError("condestJacobian: Ybus dimension $(size(Ybus, 1)) does not match bus count $(n)"))
    return (Ybus = Ybus, V = V, bus_types = bus_types, Vset = Vset, slack_idx = slack_idx)
  end
  iso_mask = falses(n)
  for b in iso
    iso_mask[b] = true
  end
  active = [k for k = 1:n if !iso_mask[k]]
  length(active) >= 2 || throw(ArgumentError("condestJacobian needs at least two active (non-isolated) buses, got $(length(active))"))
  slack_active = findfirst(==(slack_idx), active)
  slack_active === nothing && throw(ArgumentError("condestJacobian: slack bus $(slack_idx) is isolated"))
  # Accept both full-dimension and already reduced Ybus. The reduced form
  # from createYBUS walks the active buses in ascending order, which is
  # exactly the order of `active` here.
  Yact = if size(Ybus, 1) == n
    Ybus[active, active]
  elseif size(Ybus, 1) == length(active)
    Ybus
  else
    throw(ArgumentError("condestJacobian: Ybus dimension $(size(Ybus, 1)) matches neither bus count $(n) nor active count $(length(active))"))
  end
  return (Ybus = Yact, V = V[active], bus_types = bus_types[active], Vset = Vset[active], slack_idx = slack_active)
end

# Lazy condition estimate for the system the solver just factored. Called
# by the rectangular solver's finalization; the returned zero-argument
# closure is stored in the solver status (field `jacobian_condest`) and
# computes kappa only when a report actually asks for it, so a solve pays
# nothing for it. The closure snapshots the per-solve vectors (cheap O(n)
# copies; Ybus is shared, it is not mutated after assembly) and caches its
# result, so repeated report calls factorize at most once. Failures are NOT
# caught here: the consumer decides how to log them (see
# `_jacobian_condest`).
function _make_rectangular_condest_thunk(Ybus, V::Vector{ComplexF64}, bus_types::Vector{Symbol}, Vset::Vector{Float64}, slack_idx::Int, iso_nodes::Vector{Int})
  Vc = copy(V)
  btc = copy(bus_types)
  vsc = copy(Vset)
  isoc = copy(iso_nodes)
  cache = Ref{Union{Missing,Float64}}(missing)
  return function ()
    cache[] === missing || return cache[]::Float64
    inp = _active_condest_inputs(Ybus, Vc, btc, vsc, slack_idx, isoc)
    J = build_rectangular_jacobian_pq_pv_sparse(inp.Ybus, inp.V, inp.bus_types, inp.Vset, inp.slack_idx)
    kappa = condestJacobian(J)
    cache[] = kappa
    return kappa
  end
end

"""
    condestJacobian(net::Net; exact::Bool = false) -> Float64

Estimate the condition number of the rectangular Newton-Raphson Jacobian
at the net's CURRENT voltage state (`_vm_pu` / `_va_deg` per node): after a
solve that is the operating point, after a failed solve the last iterate,
which is exactly where conditioning questions arise. Builds the same
sparse PQ/PV Jacobian the solver factors (slack row fixed, PV rows as
magnitude equations) and forwards to the matrix method.

De-energized buses (`net.isoNodes`) are excluded: they carry zero Ybus
rows and no Jacobian equations, so the estimate describes the active
subsystem the solver actually solves. Prefer the value stored by the last
rectangular solve (see `_jacobian_condest`) when one exists: for nets with
closed bus links the solver works on an internally contracted copy, which
this standalone reconstruction cannot reproduce.

Throws an `ArgumentError` when the net has no slack bus among the active
buses or fewer than two active buses.
"""
function condestJacobian(net::Net; exact::Bool = false)
  nodes = net.nodeVec
  n = length(nodes)
  n >= 2 || throw(ArgumentError("condestJacobian needs a net with at least two buses"))
  bus_types = Vector{Symbol}(undef, n)
  for (k, node) in enumerate(nodes)
    t = getNodeType(node)
    # isolated buses map to neutral PQ, mirroring the solver; they are
    # filtered out below anyway
    bus_types[k] = t == Slack ? :Slack : t == PV ? :PV : :PQ
  end
  slack_idx = findfirst(==(:Slack), bus_types)
  slack_idx === nothing && throw(ArgumentError("condestJacobian: net has no slack bus"))
  V = ComplexF64[something(node._vm_pu, 1.0) * cis(deg2rad(something(node._va_deg, 0.0))) for node in nodes]
  Vset = _bus_voltage_setpoints_from_prosumers(net)
  # createYBUS assembles over active buses only, so its dimension is
  # already reduced whenever net.isoNodes is non-empty; the reducer accepts
  # that form directly (this used to raise a DimensionMismatch).
  Ybus = createYBUS(net = net)
  inp = _active_condest_inputs(Ybus, V, bus_types, Vset, slack_idx, net.isoNodes)
  J = build_rectangular_jacobian_pq_pv_sparse(inp.Ybus, inp.V, inp.bus_types, inp.Vset, inp.slack_idx)
  return condestJacobian(J; exact = exact)
end

"""
    _jacobian_condest(net::Net; warn_on_failure = true, context = "result output") -> Union{Nothing,Float64}

Jacobian condition estimate for reporting. Prefers the lazy estimate the
last rectangular solve stored in its status (field `jacobian_condest`, a
cached closure over the exact system the solver factored, including link
contraction and Q-limit active-set state). Falls back to the standalone
reconstruction `condestJacobian(net)` when no solver status carries one
(e.g. APSLF or DC runs). Returns `nothing` when the estimate fails; the
failure is logged as `@warn` (default) or `@debug`
(`warn_on_failure = false`) with `context` as the message prefix, and
never propagates into the report.
"""
function _jacobian_condest(net::Net; warn_on_failure::Bool = true, context::String = "result output", require_thunk::Bool = false)
  rect_status = rectangular_pf_status(net)
  thunk = rect_status !== nothing && hasproperty(rect_status, :jacobian_condest) ? rect_status.jacobian_condest : nothing
  # require_thunk: run-overview metadata only rates the system the solver
  # actually factored; without a thunk (APSLF, DC) it reports nothing
  # instead of reconstructing a Jacobian the run never used.
  thunk === nothing && require_thunk && return nothing
  try
    thunk !== nothing && return thunk()::Float64
    return condestJacobian(net)
  catch err
    # a single active bus (e.g. a two-bus net whose line was opened at one
    # terminal) has no meaningful condition number: expected model state,
    # not a failure, so it never warns in the result output
    degenerate = err isa ArgumentError && occursin("at least two active", err.msg)
    if warn_on_failure && !degenerate
      @warn "$(context): Jacobian condition estimate skipped" exception = err
    else
      @debug "$(context): Jacobian condition estimate skipped" exception = err
    end
    return nothing
  end
end

# Verdict shared by reportCondition and the result/diagnostics writers so
# every surface grades identically.
function _condition_verdict(kappa::Float64)::String
  return kappa < 1e6 ? "well conditioned" : kappa < 1e10 ? "borderline" : kappa < 1e14 ? "poorly conditioned, convergence at risk" : "numerically singular (Float64 exhausted)"
end

# One formatted report line, shared by reportCondition, the classic result
# log (always on since 0.9.7) and diagnose.log. Reports the
# attainable relative accuracy kappa * eps(Float64) instead of an abstract
# digits-lost count; when the solver tolerance is passed, states directly
# whether that tolerance is reachable at this conditioning.
function _condition_report_line(kappa::Float64; tol::Union{Nothing,Float64} = nothing)::String
  floor_rel = kappa * eps(Float64)
  base = "kappa1(J) = $(round(kappa, sigdigits = 3)), attainable accuracy ~ $(round(floor_rel, sigdigits = 2)), $(_condition_verdict(kappa))"
  tol === nothing && return base
  return base * (floor_rel <= tol ? " (tol $(tol) reachable)" : " (tol $(tol) NOT reachable at this conditioning)")
end

"""
    reportCondition(J::AbstractMatrix{<:Real}; iter::Int = -1) -> Float64

Print the estimated condition number of a Newton-Raphson Jacobian together
with the attainable relative accuracy (`kappa * eps(Float64)`) and a
plain-language verdict, then return the estimate. Pass `iter >= 0` to prefix
the line with the Newton iteration number when logging inside a solver loop.
"""
function reportCondition(J::AbstractMatrix{<:Real}; iter::Int = -1)
  kappa = condestJacobian(J)
  tag = iter >= 0 ? "it=$(iter) " : ""
  println(tag, _condition_report_line(kappa))
  return kappa
end