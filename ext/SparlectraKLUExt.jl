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

# file: ext/SparlectraKLUExt.jl
# purpose: the KLU sparse LU for power mode (0.30.1). Loaded by `using KLU`
#          next to Sparlectra; the base package does not depend on KLU.
#          On power-flow Jacobians of 5.7k to 60k unknowns KLU's numeric
#          refactorization is 6 to 20 times faster than UMFPACK's lu!; on
#          very large Jacobians with heavy fill-in (an 82000-bus synthetic
#          case) UMFPACK is faster, which is why KLU is not the default of
#          a single run. One context per solve and per task:
#          a KLU factorization is never shared across threads.

module SparlectraKLUExt

using Sparlectra
using KLU
using SparseArrays
using LinearAlgebra: ldiv!

mutable struct KluReuseNewtonContext <: Sparlectra.AbstractNewtonSolverContext
  fact::Union{Nothing,KLU.KLUFactorization{Float64,Int64}}
  nvar::Int
  colptr::Vector{Int64}
  rowval::Vector{Int64}
  analyze_count::Int
  refactor_count::Int
  fallback_count::Int
  rhs::Vector{Float64}
  sol::Vector{Float64}
  assembly::Sparlectra.RectangularJacobianAssembly
end

KluReuseNewtonContext() = KluReuseNewtonContext(nothing, 0, Int64[], Int64[], 0, 0, 0, Float64[], Float64[], Sparlectra.RectangularJacobianAssembly())

Sparlectra._newton_full_factorization(::KluReuseNewtonContext, J::SparseMatrixCSC{Float64,Int64}) = klu(J)
# klu! refactors numerically on the symbolic analysis of the stored object
Sparlectra._newton_refactorize!(ctx::KluReuseNewtonContext, J::SparseMatrixCSC{Float64,Int64}) = klu!(ctx.fact, J)
Sparlectra._newton_context_backend(::KluReuseNewtonContext) = :klu

function __init__()
  Sparlectra._POWER_MODE_LINEAR_CONTEXT[] = () -> KluReuseNewtonContext()
  return nothing
end

end
