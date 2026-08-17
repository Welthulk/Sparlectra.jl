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

"""
    condestJacobian(net::Net; exact::Bool = false) -> Float64

Estimate the condition number of the rectangular Newton-Raphson Jacobian
at the net's CURRENT voltage state (`_vm_pu` / `_va_deg` per node): after a
solve that is the operating point, after a failed solve the last iterate,
which is exactly where conditioning questions arise. Builds the same
sparse PQ/PV Jacobian the solver factors (slack row fixed, PV rows as
magnitude equations) and forwards to the matrix method.

Throws an `ArgumentError` when the net has no slack bus or fewer than two
buses.
"""
function condestJacobian(net::Net; exact::Bool = false)
  nodes = net.nodeVec
  n = length(nodes)
  n >= 2 || throw(ArgumentError("condestJacobian needs a net with at least two buses"))
  bus_types = Vector{Symbol}(undef, n)
  for (k, node) in enumerate(nodes)
    t = getNodeType(node)
    # isolated buses stay neutral PQ rows, mirroring the solver's mapping
    bus_types[k] = t == Slack ? :Slack : t == PV ? :PV : :PQ
  end
  slack_idx = findfirst(==(:Slack), bus_types)
  slack_idx === nothing && throw(ArgumentError("condestJacobian: net has no slack bus"))
  V = ComplexF64[something(node._vm_pu, 1.0) * cis(deg2rad(something(node._va_deg, 0.0))) for node in nodes]
  Vset = _bus_voltage_setpoints_from_prosumers(net)
  Ybus = createYBUS(net = net)
  J = build_rectangular_jacobian_pq_pv_sparse(Ybus, V, bus_types, Vset, slack_idx)
  return condestJacobian(J; exact = exact)
end

# Verdict shared by reportCondition and the result/diagnostics writers so
# every surface grades identically.
function _condition_verdict(kappa::Float64)::String
  return kappa < 1e6 ? "well conditioned" : kappa < 1e10 ? "borderline" : kappa < 1e14 ? "poorly conditioned, convergence at risk" : "numerically singular (Float64 exhausted)"
end

# One formatted report line, shared by reportCondition, the classic result
# log (output.condition_number = true) and diagnose.log.
function _condition_report_line(kappa::Float64)::String
  digits_lost = round(Int, log10(kappa))
  return "kappa1(J) = $(round(kappa, sigdigits = 3)) (about $(digits_lost) digits lost, $(_condition_verdict(kappa)))"
end

"""
    reportCondition(J::AbstractMatrix{<:Real}; iter::Int = -1) -> Float64

Print the estimated condition number of a Newton-Raphson Jacobian together
with the rough number of significant digits lost and a plain-language
verdict, then return the estimate. Pass `iter >= 0` to prefix the line
with the Newton iteration number when logging inside a solver loop.
"""
function reportCondition(J::AbstractMatrix{<:Real}; iter::Int = -1)
  kappa = condestJacobian(J)
  tag = iter >= 0 ? "it=$(iter) " : ""
  println(tag, _condition_report_line(kappa))
  return kappa
end
