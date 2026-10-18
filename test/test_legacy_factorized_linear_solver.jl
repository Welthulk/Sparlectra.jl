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

# file: test/test_legacy_factorized_linear_solver.jl
# purpose: legacy profile: the testsets of test/test_factorized_linear_solver.jl that pin an earlier
#          change (an issue, a bug, a bit-identity or off-path guarantee).
#          They moved here unchanged from test_factorized_linear_solver.jl; the helpers they use
#          are defined in that file, which the legacy group includes first.
using Sparlectra
using Test

function run_legacy_factorized_linear_solver_tests()
  # cut from testset: Factorized linear-solver backend (umfpack_reuse)
    @testset "Structural guard catches silent pattern drift" begin (function ()
      ctx = UmfpackReuseNewtonContext()
      rhs = [1.0, 2.0]
      J1 = sparse([1.0 2.0; 3.0 4.0])
      x1 = solve_newton_factorized!(ctx, J1, rhs; pattern_changed = false)
      @test ctx.analyze_count == 1
      @test ctx.refactor_count == 0
      @test norm(J1 * x1 - rhs) < 1e-12

      # Same pattern, new values: numeric refactorization only.
      J2 = copy(J1)
      nonzeros(J2)[1] = 10.0
      x2 = solve_newton_factorized!(ctx, J2, rhs; pattern_changed = false)
      @test ctx.analyze_count == 1
      @test ctx.refactor_count == 1
      @test norm(J2 * x2 - rhs) < 1e-12

      # Altered sparsity pattern with pattern_changed = false: the guard must
      # trigger a re-analysis instead of a wrong refactorization.
      J3 = sparse([1.0 0.0; 0.0 4.0])
      x3 = solve_newton_factorized!(ctx, J3, rhs; pattern_changed = false)
      @test ctx.analyze_count == 2
      @test norm(J3 * x3 - rhs) < 1e-12
      @test ctx.fallback_count == 0

      # Explicit pattern_changed forces a fresh analysis even for an
      # unchanged structure (active-set switch semantics).
      x4 = solve_newton_factorized!(ctx, J3, rhs; pattern_changed = true)
      @test ctx.analyze_count == 3
      @test norm(J3 * x4 - rhs) < 1e-12
    end)() end

    @testset "In-place Jacobian assembly matches structural build" begin (function ()
      # 4-bus synthetic system exercising PQ+PV rows and the injection chain
      # terms (duplicate-triplet merging) — issue #292 stage 3.
      y12 = 1.0 / (0.01 + 0.1im)
      y23 = 1.0 / (0.02 + 0.15im)
      y34 = 1.0 / (0.015 + 0.12im)
      y14 = 1.0 / (0.03 + 0.2im)
      Ybus = sparse(
        ComplexF64[
          y12+y14 -y12 0 -y14
          -y12 y12+y23 -y23 0
          0 -y23 y23+y34 -y34
          -y14 0 -y34 y14+y34
        ],
      )
      V = ComplexF64[1.0, 0.98 + 0.02im, 1.01 - 0.01im, 0.99 + 0.03im]
      Vset = abs.(V)
      bus_types = [:PQ, :PQ, :PV, :PQ]
      slack_idx = 1
      dP = [0.0, 0.4, 0.0, 0.2]
      dQ = [0.0, 0.1, 0.0, 0.3]
      build = (Vx, bt; kwargs...) -> Sparlectra.build_rectangular_jacobian_pq_pv(Ybus, Vx, bt, Vset, slack_idx; dPinj_dVm = dP, dQinj_dVm = dQ, structural_pattern = true, kwargs...)

      asm = Sparlectra.RectangularJacobianAssembly()
      J1 = build(V, bus_types; assembly = asm)
      Jref = build(V, bus_types)
      @test J1.colptr == Jref.colptr
      @test J1.rowval == Jref.rowval
      @test J1.nzval ≈ Jref.nzval
      @test asm.valid

      # Second call with a new iterate: in-place refresh of the SAME matrix.
      V2 = V .* (1.0 .+ 0.01im)
      J2 = build(V2, bus_types; assembly = asm)
      @test J2 === J1
      Jref2 = build(V2, bus_types)
      @test J2.colptr == Jref2.colptr
      @test J2.rowval == Jref2.rowval
      @test J2.nzval ≈ Jref2.nzval

      # Active-set change (PV -> PQ): invalidate, structural rebuild matches.
      asm.valid = false
      bus_types_pq = [:PQ, :PQ, :PQ, :PQ]
      J3 = build(V2, bus_types_pq; assembly = asm)
      Jref3 = build(V2, bus_types_pq)
      @test J3.colptr == Jref3.colptr
      @test J3.rowval == Jref3.rowval
      @test J3.nzval ≈ Jref3.nzval
      @test J3 !== J2
    end)() end

    @testset "State estimation normal equations: reuse keeps the old answers" begin (function ()
      # The WLS solve keeps one analysis per pattern of G, but must give what
      # `solve_linear` gave before on every route, the fallbacks for a
      # singular G (least-squares step instead of an abort) included.
      # Each matrix is solved twice with the same pattern and new values,
      # so the second solve runs on the kept analysis.
      spd = sparse([4.0 1.0 0.0; 1.0 3.0 1.0; 0.0 1.0 2.0])
      unsym = sparse([4.0 1.0 0.0; 1.5 3.0 1.0; 0.0 0.5 2.0])
      cases = [
        # (name, G, route counter that must show the reuse)
        ("symmetric (Cholesky)", spd, :chol),
        ("unsymmetric (LU)", unsym, :lu),
        ("singular symmetric", sparse([1.0 1.0 0.0; 1.0 1.0 0.0; 0.0 0.0 0.0]), :none),
        ("singular unsymmetric", sparse([1.0 2.0 0.0; 1.0 2.0 0.0; 0.0 1.0 0.0]), :none),
        ("diagonal", sparse(Diagonal([2.0, 3.0, 4.0])), :none),
      ]
      g = [1.0, 2.0, 3.0]
      for (name, G, route) in cases
        ctx = Sparlectra._SENormalEquationsSolver()
        for scale in (1.0, 1.5)
          Gs = copy(G)
          Gs.nzval .*= scale
          x_old = Sparlectra.solve_linear(Gs, g; allow_pinv = true, svd_max_n = 20_000)
          x_new = Sparlectra._se_solve_normal_equations!(ctx, Gs, g)
          @test isapprox(x_new, x_old; rtol = 1e-12, atol = 1e-12)
        end
        if route === :chol
          @test ctx.chol_analyze_count == 1 && ctx.chol_refactor_count == 1
        elseif route === :lu
          @test ctx.lu_ctx.analyze_count == 1 && ctx.lu_ctx.refactor_count == 1 && ctx.lu_ctx.fallback_count == 0
        end
      end
    end)() end

  return true
end
