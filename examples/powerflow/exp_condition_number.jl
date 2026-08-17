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

# Date: 2026-08-17
# file: examples/powerflow/exp_condition_number.jl
# purpose: estimates the condition number of a power-flow Jacobian at the
#          solved operating point (condestJacobian / reportCondition) and
#          shows how conditioning degrades on a stressed system

using Sparlectra
using LinearAlgebra
using SparseArrays

include(joinpath(@__DIR__, "..", "internal", "example_header.jl"))

# The rectangular Newton-Raphson solves a REAL system; assemble the complex
# power-flow Jacobian blocks dS/dV, dS/dconj(V) into the equivalent real
# 2n x 2n form so the estimator sees the matrix class the solver factors.
function real_block_jacobian(net)
  Ybus = createYBUS(net = net)
  V = [something(n._vm_pu, 1.0) * cis(deg2rad(something(n._va_deg, 0.0))) for n in net.nodeVec]
  J11, J12, _, _ = Sparlectra.build_complex_jacobian(Ybus, V)
  # dS/dV = J11 + J12 acting on (dV, dconj(V)); in real coordinates
  # (Re V, Im V) the differential dS = J11 dV + J12 dconj(V) becomes
  # [Re; Im] = [Re(J11+J12)  -Im(J11-J12); Im(J11+J12)  Re(J11-J12)] [dRe; dIm]
  A = J11 + J12
  B = J11 - J12
  J = [real(A) -imag(B); imag(A) real(B)]
  # without the slack constraint the full dS/dV is structurally singular
  # (the reference voltage is free); drop the slack bus row/column pair
  # exactly like the solver fixes the reference, leaving the regular
  # all-PQ load-flow Jacobian
  n = length(V)
  slack = findfirst(nd -> Sparlectra.getNodeType(nd) == Sparlectra.Slack, net.nodeVec)
  keep = setdiff(1:(2 * n), (slack, n + slack))
  return J[keep, keep]
end

function main(; casefile::AbstractString = "case14.m")
  print_example_banner("examples/powerflow/exp_condition_number.jl", "estimates the condition number of a power-flow Jacobian at the solved operating point")
  case_path = ensure_casefile(casefile)

  net = createNetFromMatPowerFile(filename = case_path)
  ite, erg = runpf!(net, 40, 1e-8, 0; method = :rectangular)
  erg == 0 || error("power flow did not converge (erg=$(erg))")
  println("solved $(casefile) in $(ite) iteration(s)")

  J = real_block_jacobian(net)
  println()
  println("=== Jacobian conditioning at the solution ===")
  kappa = reportCondition(J)

  # exact 2-norm condition number for comparison (dense SVD; only sensible
  # for small nets like case14)
  println("exact 2-norm condition number: ", round(cond(Matrix(J), 2), sigdigits = 3))

  # a nearly isolated bus (all incident admittances scaled down) drives the
  # Jacobian towards singularity; the estimator makes that visible without
  # waiting for the solver to stall
  println()
  println("=== Stressed variant: one nearly disconnected bus ===")
  Jbad = copy(J)
  weak = size(J, 1) ÷ 2
  Jbad[weak, :] .*= 1e-9
  Jbad[:, weak] .*= 1e-9
  Jbad[weak, weak] = max(abs(Jbad[weak, weak]), 1e-12)
  reportCondition(Jbad)

  return (kappa = kappa,)
end

run_example(main)
