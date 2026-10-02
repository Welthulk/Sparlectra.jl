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
#          tests assert that "klu" is rejected. Since 0.30.1 KLU returns as
#          the sparse LU of power mode through a package extension (loaded
#          by `using KLU`); the last testset below runs it when KLU is
#          loadable and says so.

using Sparlectra
const _KLU_AVAILABLE = Base.find_package("KLU") !== nothing
_KLU_AVAILABLE && @eval using KLU
using Test
using SparseArrays
using LinearAlgebra
using Random

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
  @testset "Factorized linear-solver backend (umfpack_reuse)" begin (function ()
    @testset "umfpack_reuse equivalence and counters" begin (function ()
      net_umf = createTest3BusNet()
      _, erg_umf = runpf!(net_umf, 20, 1e-8, 0; method = :rectangular, linear_solver = :umfpack)
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
    end)() end

    @testset "Refactorization reuse across iterations" begin (function ()
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
    end)() end

    @testset "Q-limit active-set pattern change re-analyzes" begin (function ()
      net = createTest3BusNet()
      setQLimits!(net = net, qmin_MVar = -1.0, qmax_MVar = 1.0, busName = "STATION1")
      _, erg = runpf!(net, 20, 1e-6, 0; method = :rectangular, linear_solver = :umfpack_reuse)
      @test erg == 0
      @test getNodeType(net.nodeVec[2]) == Sparlectra.PQ
      status = Sparlectra.rectangular_pf_status(net)
      @test status.linear_solver === :umfpack_reuse
      @test status.linear_solver_analyze_count >= 2
    end)() end

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

    @testset "Singular system falls back to the umfpack chain" begin (function ()
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

    @testset "State estimation linear solver: umfpack default, klu switch" begin (function ()
      @test Sparlectra.StateEstimationConfig(Dict{String,Any}()).linear_solver === :umfpack
      raw_klu = Dict{String,Any}("state_estimation" => Dict{String,Any}("linear_solver" => "klu"))
      @test Sparlectra.StateEstimationConfig(raw_klu).linear_solver === :klu
      raw_bad = Dict{String,Any}("state_estimation" => Dict{String,Any}("linear_solver" => "nonsense"))
      @test_throws ArgumentError Sparlectra.StateEstimationConfig(raw_bad)
      # the default is UMFPACK even in a session that loaded KLU for power
      # mode (this file loads it when it can): the estimator's arithmetic
      # must not depend on what else the session loaded
      @test Sparlectra._newton_context_backend(Sparlectra._SENormalEquationsSolver().lu_ctx) === :umfpack_reuse
      scf_dir = joinpath(dirname(@__DIR__), "data", "scf")
      load14() = begin
        net = importSCF(joinpath(scf_dir, "sp_case14.scf.json"))
        readMeasurementsCSV!(net; file = joinpath(scf_dir, "sp_case14.measurements.csv"))
        net
      end
      ref = runse!(load14())
      # klu without the extension: UMFPACK and one warning naming the
      # extension; driven through the configuration, so the key's way into
      # the solver is what is tested
      saved = Sparlectra._POWER_MODE_LINEAR_CONTEXT[]
      try
        Sparlectra._POWER_MODE_LINEAR_CONTEXT[] = nothing
        res = @test_logs (:warn, r"SparlectraKLUExt") match_mode = :any with_state_estimation_config(() -> runse!(load14()); linear_solver = :klu)
        @test res.converged == ref.converged && res.iterations == ref.iterations
        @test res.voltages == ref.voltages
      finally
        Sparlectra._POWER_MODE_LINEAR_CONTEXT[] = saved
      end
      if _KLU_AVAILABLE
        println("      state estimation klu backend: RAN")
        @test Sparlectra._newton_context_backend(Sparlectra._SENormalEquationsSolver(:klu).lu_ctx) === :klu
        res = with_state_estimation_config(() -> runse!(load14()); linear_solver = :klu)
        @test res.converged == ref.converged && res.iterations == ref.iterations
        @test maximum(abs.(res.voltages .- ref.voltages)) <= 1e-10
      else
        println("      state estimation klu backend: SKIPPED (KLU.jl not loadable in this session)")
      end
    end)() end

    @testset "State estimation gain matrix: symmetric_gain mirrors the upper triangle" begin (function ()
      # the plain product H'WH is symmetric only up to rounding; this H
      # is chosen so it is not (precondition, otherwise the test proves
      # nothing), and then the plain G goes to LU
      rng = Random.Xoshiro(1)
      H = sprand(rng, 60, 20, 0.3)
      w = rand(rng, 60) .* 1e3
      G0 = Sparlectra._se_gain_matrix(H, w, false)
      @test !ishermitian(G0)
      @test G0 == H' * (Diagonal(w) * H)
      G = Sparlectra._se_gain_matrix(H, w, true)
      @test ishermitian(G)
      @test triu(G) == triu(G0)   # the upper triangle as computed, bit for bit
      g = H' * (w .* ones(60))
      ctx = Sparlectra._SENormalEquationsSolver()
      dx = Sparlectra._se_solve_normal_equations!(ctx, G, g)
      @test ctx.chol_analyze_count == 1 && ctx.lu_ctx.analyze_count == 0
      dx0 = Sparlectra._se_solve_normal_equations!(Sparlectra._SENormalEquationsSolver(), G0, g)
      @test norm(dx - dx0) <= 1e-10 * norm(dx0)
      @test Sparlectra.StateEstimationConfig(Dict{String,Any}()).symmetric_gain === false
      @test Sparlectra.StateEstimationConfig(Dict{String,Any}("state_estimation" => Dict{String,Any}("symmetric_gain" => true))).symmetric_gain === true
      # the key reaches the estimator: same iterations, estimates at the
      # rounding level of the other factorization
      scf_dir = joinpath(dirname(@__DIR__), "data", "scf")
      load14() = begin
        net = importSCF(joinpath(scf_dir, "sp_case14.scf.json"))
        readMeasurementsCSV!(net; file = joinpath(scf_dir, "sp_case14.measurements.csv"))
        net
      end
      ref = runse!(load14())
      res = with_state_estimation_config(() -> runse!(load14()); symmetric_gain = true)
      @test res.converged == ref.converged && res.iterations == ref.iterations
      @test maximum(abs.(res.voltages .- ref.voltages)) <= 1e-9
    end)() end

    @testset "Multi-island run carries per-island counters" begin (function ()
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
    end)() end

    @testset "Configuration validation and defaults (klu rejected)" begin (function ()
      @test Sparlectra.PowerFlowConfig(Dict{String,Any}()).linear_solver === :umfpack_reuse
      @test powerflow_config().linear_solver === :umfpack_reuse
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
    end)() end

    @testset "Web UI option spec, rendering, and sidecar round-trip" begin (function ()
      spec = SparlectraApp._webui_option_spec("power_flow_linear_solver")
      @test spec.config_key == "power_flow.linear_solver"
      @test spec.default == "umfpack_reuse"
      @test spec.control === :select
      @test spec.section === :expert
      @test spec.save_in_case_sidecar
      @test Tuple(String(v) for v in spec.allowed_values) == ("umfpack", "umfpack_reuse")
      @test "power_flow_linear_solver" in SparlectraApp._WEBUI_CASE_PROFILE_FIELDS
      @test SparlectraApp._webui_normalize_case_profile_form_value("power_flow_linear_solver", "umfpack_reuse") == "umfpack_reuse"
      @test_throws ArgumentError SparlectraApp._webui_normalize_case_profile_form_value("power_flow_linear_solver", "klu")

      # the expert options render on the Settings page
      form_html = SparlectraApp.render_settings_page()
      expert_parts = split(form_html, "<summary>Advanced options</summary>")
      @test length(expert_parts) == 2
      expert_html = expert_parts[2]
      @test occursin("name=\"power_flow_linear_solver\"", expert_html)
      @test occursin("<option value=\"umfpack_reuse\" selected>", expert_html)
      @test occursin("<option value=\"umfpack\"", expert_html)
      # klu is no value of power_flow.linear_solver (it is one of the
      # power-mode LU select since 0.30.2, so only the linear-solver select
      # is searched)
      linear_select = match(r"<select id=\"power_flow_linear_solver\"[^>]*>.*?</select>"s, expert_html)
      @test linear_select !== nothing && !occursin("<option value=\"klu\"", linear_select.match)
      linear_topic = SparlectraApp.resolve_webui_help_topic("power_flow.linear_solver")
      @test linear_topic !== nothing && !isempty(linear_topic.hint)
      @test occursin("href=\"$(SparlectraApp.webui_help_page_url("power_flow.linear_solver"))\"", form_html)
    end)() end
    @testset "power mode: KLU extension and power_mode_lu" begin (function ()
      path = joinpath(dirname(@__DIR__), "data", "mpower", "sp_case118.m")
      state(n) = [b._vm_pu * cis(deg2rad(b._va_deg)) for b in n.nodeVec]
      lu_status(n) = Sparlectra.rectangular_pf_status(n)
      ref_net = Sparlectra.createNetFromMatPowerFile(filename = path)
      it_ref, erg_ref = runpf!(ref_net, 30, 1e-10, 0; qlimits_enabled = false)
      @test erg_ref == 0
      # without the extension (hook cleared, as in a session without
      # `using KLU`): auto is UMFPACK without a decision, klu falls back to
      # UMFPACK with one warning naming the extension
      saved = Sparlectra._POWER_MODE_LINEAR_CONTEXT[]
      try
        Sparlectra._POWER_MODE_LINEAR_CONTEXT[] = nothing
        @test Sparlectra.power_mode_linear_solver_backend() === :umfpack_reuse
        net = Sparlectra.createNetFromMatPowerFile(filename = path)
        runpf!(net, 30, 1e-10, 0; qlimits_enabled = false, power_mode = true)
        st = lu_status(net)
        @test (st.power_mode_lu, st.power_mode_lu_choice, st.power_mode_lu_source) == (:auto, :umfpack, :klu_not_loaded)
        @test isnan(st.power_mode_lu_klu_est_flops)
        @test maximum(abs.(state(net) .- state(ref_net))) <= 1e-12
        net = Sparlectra.createNetFromMatPowerFile(filename = path)
        @test_logs (:warn, r"SparlectraKLUExt") match_mode = :any runpf!(net, 30, 1e-10, 0; qlimits_enabled = false, power_mode = true, power_mode_lu = :klu)
        @test lu_status(net).power_mode_lu_choice === :umfpack
      finally
        Sparlectra._POWER_MODE_LINEAR_CONTEXT[] = saved
      end
      @test_throws ArgumentError runpf!(Sparlectra.createNetFromMatPowerFile(filename = path), 30, 1e-10, 0; qlimits_enabled = false, power_mode = true, power_mode_lu = :pardiso)
      if !_KLU_AVAILABLE
        println("      power mode KLU extension: SKIPPED (KLU.jl not loadable in this session)")
        return nothing
      end
      println("      power mode KLU extension: RAN")
      @test Sparlectra.power_mode_linear_solver_backend() === :klu
      # the rule of auto: KLU up to the threshold of KLU's symbolic flop
      # estimate, UMFPACK above
      thr = Sparlectra.POWER_MODE_LU_KLU_MAX_FLOPS
      @test Sparlectra._power_mode_lu_choice(thr) === :klu && Sparlectra._power_mode_lu_choice(nextfloat(thr)) === :umfpack
      # every value gives the same solution; klu and umfpack say so, auto
      # decides on the first solve from the symbolic estimate (sp_case118 is
      # far below the threshold: KLU)
      for (lu, backend) in ((:klu, :klu), (:umfpack, :umfpack_reuse), (:auto, :klu))
        net = Sparlectra.createNetFromMatPowerFile(filename = path)
        prof = Dict{Symbol,Any}()
        it_pm, erg_pm = runpf!(net, 30, 1e-10, 0; qlimits_enabled = false, power_mode = true, power_mode_lu = lu, performance_profile = prof)
        @test erg_pm == 0 && it_pm == it_ref
        @test maximum(abs.(state(net) .- state(ref_net))) <= 1e-12
        st = lu_status(net)
        @test st.power_mode_lu === lu && prof[:linear_solver_backend] === backend
        @test Sparlectra._newton_context_backend(net._power_cache.linear_ctx) === backend
        if lu === :auto
          @test st.power_mode_lu_source === :decided && st.power_mode_lu_choice === :klu
          @test 0 < st.power_mode_lu_klu_est_flops <= thr && occursin("KLU symbolic estimate", st.power_mode_lu_line)
        else
          @test st.power_mode_lu_source === :configured && isnan(st.power_mode_lu_klu_est_flops)
        end
      end
      # warm solves on the same net: KLU refactors on the kept analysis; under
      # auto the second solve takes the decision of the first, also after a
      # branch outage (once per network), on a copy of the net (the N-1 and
      # scenario workers) and not after reset_power_mode!
      for lu in (:klu, :auto)
        net = Sparlectra.createNetFromMatPowerFile(filename = path)
        runpf!(net, 30, 1e-10, 0; qlimits_enabled = false, power_mode = true, power_mode_lu = lu)
        first_status = lu_status(net)
        fresh = Sparlectra.createNetFromMatPowerFile(filename = path)
        for (a, b) in zip(net.nodeVec, fresh.nodeVec)
          a._vm_pu = b._vm_pu
          a._va_deg = b._va_deg
        end
        worker = deepcopy(net)
        it2, erg2 = runpf!(net, 30, 1e-10, 0; qlimits_enabled = false, power_mode = true, power_mode_lu = lu)
        @test erg2 == 0 && it2 == it_ref
        @test maximum(abs.(state(net) .- state(ref_net))) <= 1e-12
        ctx = net._power_cache.linear_ctx
        @test ctx.analyze_count == 1 && ctx.refactor_count >= it_ref && ctx.fallback_count == 0
        st = lu_status(net)
        @test st.power_mode_lu_choice === first_status.power_mode_lu_choice
        lu === :klu && continue
        @test st.power_mode_lu_source === :earlier_decision && st.power_mode_lu_klu_est_flops == first_status.power_mode_lu_klu_est_flops
        @test occursin("earlier solve", st.power_mode_lu_line)
        setBranchStatus!(net.branchVec[1], false)
        runpf!(net, 30, 1e-10, 0; qlimits_enabled = false, power_mode = true)
        @test lu_status(net).power_mode_lu_source === :earlier_decision
        runpf!(worker, 30, 1e-10, 0; qlimits_enabled = false, power_mode = true)
        @test lu_status(worker).power_mode_lu_source === :earlier_decision
        # from the file start again (a converged state needs no factorization,
        # and a solve without one decides nothing)
        Sparlectra.reset_power_mode!(net)
        for (a, b) in zip(net.nodeVec, fresh.nodeVec)
          a._vm_pu = b._vm_pu
          a._va_deg = b._va_deg
        end
        runpf!(net, 30, 1e-10, 0; qlimits_enabled = false, power_mode = true)
        @test lu_status(net).power_mode_lu_source === :decided
      end
      # island path: the island nets are fresh in every solve and record
      # into the memo of the caller's network, one decision per island
      island_net = _reuse_two_island_net()
      for round in 1:2
        _, erg = runpf!(island_net, 20, 1e-8, 0; islands_enabled = true, power_mode = true)
        @test erg == 0
        rows = lu_status(island_net).power_mode_lu_islands
        @test length(rows) == 2
        @test all(r -> r.source === (round == 1 ? :decided : :earlier_decision), rows)
        @test occursin("island 1:", lu_status(island_net).power_mode_lu_line) && occursin("island 2:", lu_status(island_net).power_mode_lu_line)
      end
      @test sort(collect(keys(island_net._power_cache.lu_memo.decisions))) == [1, 3]
    end)() end
  end)() end
end
