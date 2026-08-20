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

# file: test/test_factorized_linear_solver.jl
# purpose: tests the umfpack_reuse linear-solver backend of the rectangular
#          Newton step: UMFPACK equivalence, factorization reuse,
#          pattern-drift guards, singular fallback, config and Web UI wiring.
#          The former klu backend was removed in 0.9.10; the config and Web UI
#          tests assert that "klu" is rejected.

using Sparlectra
using Test
using SparseArrays
using LinearAlgebra

function _reuse_two_island_net()::Net
  island_net = Net(name = "reuse_islands", baseMVA = 100.0)
  for busName in ("A1", "A2", "B1", "B2")
    addBus!(net = island_net, busName = busName, vn_kV = 110.0)
  end
  addPIModelACLine!(net = island_net, fromBus = "A1", toBus = "A2", r_pu = 0.01, x_pu = 0.10, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = island_net, fromBus = "B1", toBus = "B2", r_pu = 0.01, x_pu = 0.10, b_pu = 0.0, status = 1)
  addProsumer!(net = island_net, busName = "A1", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "A1")
  addProsumer!(net = island_net, busName = "A2", type = "ENERGYCONSUMER", p = 10.0, q = 3.0)
  addProsumer!(net = island_net, busName = "B1", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "B1")
  addProsumer!(net = island_net, busName = "B2", type = "ENERGYCONSUMER", p = 8.0, q = 2.0)
  return island_net
end

function run_factorized_linear_solver_tests()
  @testset "Factorized linear-solver backend (umfpack_reuse)" begin
    @testset "umfpack_reuse equivalence and counters" begin
      net_umf = createTest3BusNet()
      _, erg_umf = runpf!(net_umf, 20, 1e-8, 0; method = :rectangular)
      net_reuse = createTest3BusNet()
      _, erg_reuse = runpf!(net_reuse, 20, 1e-8, 0; method = :rectangular, linear_solver = :umfpack_reuse)
      @test erg_umf == 0
      @test erg_reuse == 0
      for i in eachindex(net_umf.nodeVec)
        @test isapprox(net_umf.nodeVec[i]._vm_pu, net_reuse.nodeVec[i]._vm_pu; atol = 1e-10)
        @test isapprox(net_umf.nodeVec[i]._va_deg, net_reuse.nodeVec[i]._va_deg; atol = 1e-10)
      end
      status_umf = Sparlectra.rectangular_pf_status(net_umf)
      status_reuse = Sparlectra.rectangular_pf_status(net_reuse)
      @test status_umf.status == status_reuse.status
      @test status_umf.linear_solver === :umfpack
      @test status_reuse.linear_solver === :umfpack_reuse
      @test status_reuse.linear_solver_analyze_count >= 1
      @test status_reuse.linear_solver_fallback_count == 0
      @test status_umf.linear_solver_analyze_count == 0
      @test status_umf.linear_solver_refactor_count == 0
    end

    @testset "Refactorization reuse across iterations" begin
      net = createTest3BusNet()
      _, erg = runpf!(net, 20, 1e-10, 0; method = :rectangular, linear_solver = :umfpack_reuse, qlimits_enabled = false)
      @test erg == 0
      status = Sparlectra.rectangular_pf_status(net)
      @test status.linear_solver === :umfpack_reuse
      # With the structural (value-independent) Jacobian pattern exactly one
      # symbolic analysis runs; everything else is numeric refactorization.
      @test status.linear_solver_analyze_count == 1
      @test status.linear_solver_refactor_count >= 1
      @test status.linear_solver_fallback_count == 0
    end

    @testset "Q-limit active-set pattern change re-analyzes" begin
      net = createTest3BusNet()
      setQLimits!(net = net, qmin_MVar = -1.0, qmax_MVar = 1.0, busName = "STATION1")
      _, erg = runpf!(net, 20, 1e-6, 0; method = :rectangular, linear_solver = :umfpack_reuse)
      @test erg == 0
      @test getNodeType(net.nodeVec[2]) == Sparlectra.PQ
      status = Sparlectra.rectangular_pf_status(net)
      @test status.linear_solver === :umfpack_reuse
      @test status.linear_solver_analyze_count >= 2
    end

    @testset "Structural guard catches silent pattern drift" begin
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
    end

    @testset "In-place Jacobian assembly matches structural build" begin
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
    end

    @testset "Singular system falls back to the umfpack chain" begin
      ctx = UmfpackReuseNewtonContext()
      J_singular = sparse([1.0 1.0; 1.0 1.0])
      rhs = [1.0, 1.0]
      x = solve_newton_factorized!(ctx, J_singular, rhs; pattern_changed = false)
      @test ctx.fallback_count >= 1
      @test ctx.fact === nothing
      @test all(isfinite, x)
      # The context recovers after a fallback: the next solve re-analyzes.
      J_ok = sparse([2.0 0.0; 0.0 3.0])
      x_ok = solve_newton_factorized!(ctx, J_ok, [2.0, 3.0]; pattern_changed = false)
      @test ctx.analyze_count >= 1
      @test norm(J_ok * x_ok - [2.0, 3.0]) < 1e-12
    end

    @testset "Multi-island run carries per-island counters" begin
      island_net = _reuse_two_island_net()
      profile = Dict{Symbol,Any}()
      _, erg = runpf!(island_net, 20, 1e-8, 0; method = :rectangular, linear_solver = :umfpack_reuse, islands_enabled = true, performance_profile = profile)
      @test erg == 0
      statuses = profile[:ac_island_solver_statuses]
      @test length(statuses) == 2
      for (_, island_status) in statuses
        @test island_status.status == :converged
        @test island_status.linear_solver === :umfpack_reuse
        @test island_status.linear_solver_analyze_count >= 1
      end
    end

    @testset "Configuration validation and defaults (klu rejected)" begin
      @test Sparlectra.PowerFlowConfig(Dict{String,Any}()).linear_solver === :umfpack
      @test powerflow_config().linear_solver === :umfpack
      raw_reuse = Dict{String,Any}("power_flow" => Dict{String,Any}("linear_solver" => "umfpack_reuse"))
      @test Sparlectra.PowerFlowConfig(raw_reuse).linear_solver === :umfpack_reuse
      # the removed klu backend and arbitrary values are rejected alike
      raw_klu = Dict{String,Any}("power_flow" => Dict{String,Any}("linear_solver" => "klu"))
      @test_throws ArgumentError Sparlectra.PowerFlowConfig(raw_klu)
      raw_bad = Dict{String,Any}("power_flow" => Dict{String,Any}("linear_solver" => "nonsense"))
      @test_throws ArgumentError Sparlectra.PowerFlowConfig(raw_bad)
      net = createTest3BusNet()
      @test_throws ArgumentError runpf_rectangular!(net; linear_solver = :klu)
      overrides = Sparlectra.validate_gui_config_overrides(Dict{String,Any}("power_flow.linear_solver" => "umfpack_reuse"))
      @test overrides["power_flow"]["linear_solver"] == "umfpack_reuse"
      @test_throws ArgumentError Sparlectra.validate_gui_config_overrides(Dict{String,Any}("power_flow.linear_solver" => "klu"))
    end

    @testset "Web UI option spec, rendering, and sidecar round-trip" begin
      spec = Sparlectra._webui_option_spec("power_flow_linear_solver")
      @test spec.config_key == "power_flow.linear_solver"
      @test spec.default == "umfpack"
      @test spec.control === :select
      @test spec.section === :expert
      @test spec.save_in_case_sidecar
      @test Tuple(String(v) for v in spec.allowed_values) == ("umfpack", "umfpack_reuse")
      @test "power_flow_linear_solver" in Sparlectra._WEBUI_CASE_PROFILE_FIELDS
      @test Sparlectra._webui_normalize_case_profile_form_value("power_flow_linear_solver", "umfpack_reuse") == "umfpack_reuse"
      @test_throws ArgumentError Sparlectra._webui_normalize_case_profile_form_value("power_flow_linear_solver", "klu")

      form_html = Sparlectra.render_powerflow_form()
      expert_parts = split(form_html, "<summary>Advanced options</summary>")
      @test length(expert_parts) == 2
      expert_html = expert_parts[2]
      @test occursin("name=\"power_flow_linear_solver\"", expert_html)
      @test occursin("<option value=\"umfpack\" selected>", expert_html)
      @test occursin("<option value=\"umfpack_reuse\"", expert_html)
      @test !occursin("<option value=\"klu\"", expert_html)
      @test occursin("href=\"/help/power_flow.linear_solver\"", form_html)
      @test Sparlectra.resolve_webui_help_topic("power_flow.linear_solver") !== nothing
      excerpt = Sparlectra.load_webui_help_excerpt("power_flow.linear_solver")
      @test excerpt !== nothing
      @test occursin("power_flow.linear_solver", excerpt)
    end
  end
end
