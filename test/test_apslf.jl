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

# file: test/test_apslf.jl
# purpose: tests the optional APSLF (AnalyticLoadFlow.jl) solver integration:
#          config validation, Web UI form parsing, the not-installed error
#          path, adapter mapping, and case14 agreement with Newton-Raphson
#
# Tests the APSLF (AnalyticLoadFlow.jl) integration: config validation,
# controller rejection, Web UI form parsing, the adapter mapping, the
# standalone solve, and the start-value generator. AnalyticLoadFlow.jl is a
# REQUIRED dependency since 0.10.0, so nothing here is conditional any more:
# the whole surface runs in every profile that includes this file.

using Sparlectra
using Test
using Random
using LinearAlgebra
using SparseArrays

import AnalyticLoadFlow

function run_apslf_tests()
    # the helpers of the no-polish set and the platform probe live at
    # function scope
    quiet = OutputConfig(logfile_results=:off, console_summary=false, startup_latency_hint=false)
    cfg_nr = SparlectraConfig(powerflow=PowerFlowConfig(solver=:rectangular, rescue=false), output=quiet)
    cfg_ap = SparlectraConfig(powerflow=PowerFlowConfig(solver=:apslf), output=quiet)
    function ring3()
        net = Net(name="apslf_ring3", baseMVA=100.0)
        addBus!(net=net, busName="B1", vn_kV=110.0, vm_pu=1.02, va_deg=0.0)
        addBus!(net=net, busName="B2", vn_kV=110.0, vm_pu=1.0, va_deg=0.0)
        addBus!(net=net, busName="B3", vn_kV=110.0, vm_pu=1.0, va_deg=0.0)
        addPIModelACLine!(net=net, fromBus="B1", toBus="B2", r_pu=0.010, x_pu=0.080, b_pu=0.0, status=1)
        addPIModelACLine!(net=net, fromBus="B2", toBus="B3", r_pu=0.011, x_pu=0.085, b_pu=0.0, status=1)
        addPIModelACLine!(net=net, fromBus="B3", toBus="B1", r_pu=0.012, x_pu=0.090, b_pu=0.0, status=1)
        addProsumer!(net=net, busName="B1", type="EXTERNALNETWORKINJECTION", referencePri="B1", vm_pu=1.02, va_deg=0.0)
        addProsumer!(net=net, busName="B2", type="GENERATOR", p=20.0, q=5.0)
        addProsumer!(net=net, busName="B3", type="LOAD", p=30.0, q=10.0)
        ok, msg = validate!(net=net)
        ok || error("ring3 invalid: $msg")
        return net
    end
    scf5 = abspath(joinpath(dirname(@__DIR__), "data", "scf", "sp_case5.scf.json"))
    # on a failure the layers below the solver are probed and printed, so
    # a run on another machine names the first broken layer by itself
    # (complex sparse LU, ldiv! into a column view, AnalyticLoadFlow's own
    # 9-bus PV case) without a separate diagnostic session
    function _apslf_platform_probe()
        println("      APSLF platform probe: Julia ", VERSION, " on ", Sys.KERNEL, " ", Sys.MACHINE, ", BLAS ", BLAS.get_config())
        rng = Random.MersenneTwister(1)
        n = 400
        A = SparseArrays.sprandn(rng, ComplexF64, n, n, 0.02) + LinearAlgebra.I * 10
        b = randn(rng, ComplexF64, n)
        F = LinearAlgebra.lu(A)
        println("      probe 1 complex sparse LU residual: ", LinearAlgebra.norm(A * (F \ b) - b) / LinearAlgebra.norm(b))
        C = zeros(ComplexF64, n, 3)
        LinearAlgebra.ldiv!(@view(C[:, 2]), F, b)
        println("      probe 2 ldiv! into a column view residual: ", LinearAlgebra.norm(A * C[:, 2] - b) / LinearAlgebra.norm(b), ", other columns untouched: ", all(iszero, C[:, 1]) && all(iszero, C[:, 3]))
        # the Padé evaluation solves a small dense complex system with LAPACK
        # (the one BLAS/LAPACK call of the series path); the geometric series
        # of exp(0.7 + 0.3im) must evaluate to that value, and a pure-Julia
        # elimination of the same kind of system must agree with the library solve
        z8 = 0.7 + 0.3im
        c = ComplexF64[z8^n / factorial(big(n)) for n in 0:24]
        println("      probe 8 Pade [12/12] of the exp series (expected ", exp(z8), "): ", AnalyticLoadFlow.pade_eval(c, 12, 12))
        Ad = ComplexF64[c[abs(i - j) + 1] * (1 + 0.1im * (i + j)) for i in 1:12, j in 1:12] + 2I
        bd = ComplexF64[c[i] for i in 1:12]
        x_lib = Ad \ bd
        x_ref = let A2 = copy(Ad), b2 = copy(bd), n = 12
            for k in 1:n, i in (k + 1):n
                f = A2[i, k] / A2[k, k]
                for j in k:n; A2[i, j] -= f * A2[k, j]; end
                b2[i] -= f * b2[k]
            end
            x = zeros(ComplexF64, n)
            for i in n:-1:1
                acc = b2[i]
                for j in (i + 1):n; acc -= A2[i, j] * x[j]; end
                x[i] = acc / A2[i, i]
            end
            x
        end
        println("      probe 8 dense complex 12x12 solve, LAPACK against pure Julia elimination: max diff ", maximum(abs.(x_lib .- x_ref)))
        # AnalyticLoadFlow's own 9-bus case (three PV buses) through the two
        # solver modes Sparlectra can select; the :direct mode is the default
        case9 = AnalyticLoadFlow.demo_case_9bus()
        for mode in (:direct, :outer)
            res9 = AnalyticLoadFlow.solve_pf_apslf(case9; mode=mode, order=24, use_pade=true, nr_polish=false)
            println("      probe 3 AnalyticLoadFlow 9-bus PV case, mode ", mode, ", no polish: converged = ", res9.converged)
        end
        # AnalyticLoadFlow keeps a sparse Y sparse even below its own
        # sparse threshold, so a small Sparlectra model runs the sparse
        # series path while the demo above runs the dense one: the same
        # spec once as delivered (sparse) and once densified
        for (label, build_spec) in (("ring3", ring3), ("sp_case118", () -> Sparlectra.createNetFromMatPowerFile(filename=abspath(joinpath(dirname(@__DIR__), "data", "mpower", "sp_case118.m")), flatstart=false, enable_pq_gen_controllers=true, bus_shunt_model=:admittance, matpower_shift_sign=1.0, matpower_shift_unit=:deg, matpower_ratio=:normal, tap_changer_model=:ideal)))
            spec = Sparlectra._apslf_spec_from_model(Sparlectra.buildPfModel(build_spec()))
            for (form, Y) in (("sparse", SparseArrays.sparse(spec.Y)), ("dense", Matrix(spec.Y)))
                res = AnalyticLoadFlow.solve_pf_apslf(merge(spec, (Y=Y,)); mode=:direct, order=24, use_pade=true, nr_polish=false)
                println("      probe 5 AnalyticLoadFlow on the Sparlectra spec of ", label, ", Y ", form, ": converged = ", res.converged)
            end
            # the adapter asks for the coefficients as well (the radius needs them)
            res_c = AnalyticLoadFlow.solve_pf_apslf(spec; mode=:direct, order=24, use_pade=true, nr_polish=false, return_coeffs=true)
            println("      probe 6 the same with return_coeffs = true (as the adapter calls it): converged = ", res_c.converged)
            # is the voltage behind that flag right? ALF's own mismatch of its
            # result, Sparlectra's mismatch of the same result, and the
            # distance to the Newton solution of the same model
            model_p = Sparlectra.buildPfModel(build_spec())
            net_nr = build_spec()
            runpf!(net_nr, 30, 1e-10, 0)
            vm_nr_p = [net_nr.nodeVec[i]._vm_pu for i in model_p.busIdx_net]
            alf_mis = AnalyticLoadFlow.max_mismatch_on_specY(spec, res_c.V)
            bt_p = Symbol[Sparlectra._apslf_bus_type(b) for b in res_c.bustype]
            S_p = ComplexF64[complex(real(model_p.Sspec[i]), res_c.Q[i]) for i in eachindex(model_p.Sspec)]
            sp_mis = maximum(abs.(Sparlectra.mismatch_rectangular(model_p.Ybus, res_c.V, S_p, bt_p, model_p.Vset, model_p.slack_idx)))
            println("      probe 7 ", label, ": ALF mismatch of its own V ", alf_mis, ", Sparlectra mismatch of that V ", sp_mis, ", max |V| deviation from Newton ", maximum(abs.(abs.(res_c.V) .- vm_nr_p)), ", ALF Q on bus 1..3 ", res_c.Q[1:min(3, end)], ", Sspec ", model_p.Sspec[1:min(3, end)])
        end
        # the same two modes on Sparlectra's ring3 through the adapter, against NR
        ref = ring3()
        runpf!(ref, 30, 1e-10, 0)
        vm_ref = getfield.(ref.nodeVec, :_vm_pu)
        # the adapter path in its three layers, each printed with the
        # numbers a run on another machine can be compared against: the
        # model as runpf_external! builds it (no Q-limits), AnalyticLoadFlow
        # on exactly that spec, solvePf on it, then the full runpf_external!
        model_4 = Sparlectra.buildPfModel(ring3(); include_limits=false)
        spec_4 = Sparlectra._apslf_spec_from_model(model_4)
        println("      probe 4 model: busType ", model_4.busType, ", slack ", model_4.slack_idx, ", Sspec ", model_4.Sspec, ", Vset ", model_4.Vset, ", Qmin ", spec_4.Qmin, ", Qmax ", spec_4.Qmax)
        println("      probe 4 model: Ybus ", Matrix(model_4.Ybus))
        println("      probe 4 Newton |V| in PF order: ", vm_ref[model_4.busIdx_net])
        res_4 = AnalyticLoadFlow.solve_pf_apslf(spec_4; mode=:direct, order=24, use_pade=true, nr_polish=false, return_coeffs=true)
        println("      probe 4 AnalyticLoadFlow on that spec: converged = ", res_4.converged, ", reason ", get(res_4, :reason, :none), ", outer_iters ", res_4.outer_iters, ", |V| ", abs.(res_4.V), ", Q ", res_4.Q, ", bustype ", res_4.bustype, ", ALF mismatch ", AnalyticLoadFlow.max_mismatch_on_specY(spec_4, res_4.V))
        sol_4 = Sparlectra.solvePf(apslf_solver(order=24, use_pade=true, nr_polish=false, mode=:direct), model_4; tol=1e-8)
        println("      probe 4 solvePf on that model: converged = ", sol_4.converged, ", residual_inf ", sol_4.residual_inf, ", |V| ", abs.(sol_4.V), ", bustype_final ", sol_4.meta.bustype_final)
        for mode in (:direct, :outer)
            probe_net = ring3()
            _, st, sol = runpf_external!(probe_net, apslf_solver(order=24, use_pade=true, nr_polish=false, mode=mode); tol=1e-8)
            println("      probe 4 runpf_external! mode ", mode, ": residual_inf ", sol.residual_inf, ", |V| ", abs.(sol.V), ", net |V| ", getfield.(probe_net.nodeVec, :_vm_pu), ", ALF mismatch of that V ", AnalyticLoadFlow.max_mismatch_on_specY(spec_4, sol.V))
            println("      probe 4 Sparlectra ring3 through the adapter, mode ", mode, ": status ", st, ", converged = ", sol.converged, ", max |dVm| vs NR ", maximum(abs.(getfield.(probe_net.nodeVec, :_vm_pu) .- vm_ref)))
        end
    end
    # SPARLECTRA_APSLF_PROBE=1 prints the probe on a passing run too
    get(ENV, "SPARLECTRA_APSLF_PROBE", "") == "1" && _apslf_platform_probe()
    # one outer set, so a standalone call still runs every group after a
    # failing one (a top-level testset throws at its end)
    @testset "APSLF" begin
        @testset "APSLF (AnalyticLoadFlow.jl) integration" begin
            @testset "Config validation: power_flow.solver / apslf / apslf_start" begin
                default_cfg = Sparlectra.SparlectraConfig()
                @test default_cfg.powerflow.solver === :rectangular
                @test default_cfg.powerflow.apslf.order == 24
                @test default_cfg.powerflow.apslf.convergence_radius === true
                @test default_cfg.powerflow.apslf.use_pade === true
                @test default_cfg.powerflow.apslf.nr_polish === false
                @test default_cfg.powerflow.apslf_start.enabled === false
                @test default_cfg.powerflow.apslf_start.order == 40

                @test_throws ArgumentError Sparlectra.SparlectraConfig(Dict("power_flow" => Dict("solver" => "polar")))
                cfg_apslf = Sparlectra.SparlectraConfig(Dict("power_flow" => Dict("solver" => "apslf")))
                @test cfg_apslf.powerflow.solver === :apslf

                @test_throws ArgumentError Sparlectra.SparlectraConfig(Dict("power_flow" => Dict("apslf" => Dict("order" => 0))))
                @test_throws ArgumentError Sparlectra.SparlectraConfig(Dict("power_flow" => Dict("apslf_start" => Dict("order" => -1))))

                # solver=apslf together with apslf_start.enabled=true is rejected (the
                # start-value generator only makes sense ahead of the NR solve).
                @test_throws ArgumentError Sparlectra.SparlectraConfig(Dict("power_flow" => Dict("solver" => "apslf", "apslf_start" => Dict("enabled" => true))))

                # apslf_start.enabled=true alone (solver stays rectangular) is fine.
                cfg_hybrid = Sparlectra.SparlectraConfig(Dict("power_flow" => Dict("apslf_start" => Dict("enabled" => true, "order" => 12))))
                @test cfg_hybrid.powerflow.apslf_start.enabled
                @test cfg_hybrid.powerflow.apslf_start.order == 12
                @test cfg_hybrid.powerflow.solver === :rectangular

                # Unknown keys under the new sections are rejected like any other section.
                bad_key_file = test_scratch_path(".yaml")
                write(bad_key_file, "power_flow:\n  apslf:\n    bogus_key: 1\n")
                # the file carries no config_version on purpose; that once-per-session
                # warning is expected here and must not reach the test output
                run_with_expected_warnings(["declares no config_version"]) do
                    @test_throws ArgumentError Sparlectra.load_sparlectra_config(bad_key_file; reload=true)
                end
            end

            @testset "Controller + APSLF solver rejection" begin
                # load_fixture_net: the shipped sp_case14 carries a REAL declared tap
                # controller (no download, no hand-attached controller)
                net = Sparlectra.importSCF(abspath(joinpath(dirname(@__DIR__), "data", "scf", "sp_case14.scf.json")))
                @test length(Sparlectra.collect_outer_controllers(net)) == 1

                cfg = Sparlectra.SparlectraConfig(powerflow=Sparlectra.PowerFlowConfig(solver=:apslf), output=OutputConfig(logfile_results=:off))
                # No AnalyticLoadFlow.jl is required here: the controller check runs
                # before the solver is even constructed, so this rejects identically
                # regardless of whether the extension is loaded.
                @test_throws ArgumentError run_sparlectra(net=net, config=cfg)
            end

            @testset "WebUI form parsing for APSLF fields -> effective config" begin
                form = Dict{String,Any}(
                    "casefile" => "case14.m",
                    "power_flow_solver" => "apslf",
                    "power_flow_apslf_order" => "25",
                    "power_flow_apslf_use_pade" => "false",
                    "power_flow_apslf_nr_polish" => "true",
                    "power_flow_apslf_start_enabled" => "false",
                    "power_flow_apslf_start_order" => "40",
                )
                request = Sparlectra.powerflow_webui_request(form)
                overrides = request["config_overrides"]
                @test overrides["power_flow.solver"] == "apslf"
                @test overrides["power_flow.apslf.order"] === 25
                @test overrides["power_flow.apslf.use_pade"] === false
                @test overrides["power_flow.apslf.nr_polish"] === true
                @test overrides["power_flow.apslf_start.enabled"] === false
                @test overrides["power_flow.apslf_start.order"] === 40

                nested = Sparlectra.validate_gui_config_overrides(overrides)
                cfg, _ = Sparlectra._load_api_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, nested)
                @test cfg.powerflow.solver === :apslf
                @test cfg.powerflow.apslf.order == 25
                @test cfg.powerflow.apslf.use_pade === false
                @test cfg.powerflow.apslf.nr_polish === true

                # Unchecked checkbox explicitly submits "false" (hidden+checkbox HTML
                # pairing), the same semantics as the Q-limit checkbox pattern: the key
                # is present with an explicit false, not omitted.
                disabled_form = copy(form)
                disabled_form["power_flow_apslf_use_pade"] = "false"
                @test Sparlectra.powerflow_webui_request(disabled_form)["config_overrides"]["power_flow.apslf.use_pade"] === false

                # A field genuinely absent from the submitted form (e.g. a stale
                # client) is skipped rather than defaulted to false.
                missing_field_form = copy(form)
                delete!(missing_field_form, "power_flow_apslf_start_enabled")
                @test !haskey(Sparlectra.powerflow_webui_request(missing_field_form)["config_overrides"], "power_flow.apslf_start.enabled")

                # GUI-layer validation rejects the same invalid values as the core config.
                @test_throws ArgumentError Sparlectra.validate_gui_config_overrides(Dict{String,Any}("power_flow.solver" => "polar"))
                @test_throws ArgumentError Sparlectra.validate_gui_config_overrides(Dict{String,Any}("power_flow.apslf.order" => 0))

                # Downstream config assembly still rejects the solver+start conflict
                # even when both values arrive through the GUI override path.
                conflict_overrides = Dict{String,Any}("power_flow.solver" => "apslf", "power_flow.apslf_start.enabled" => true)
                conflict_nested = Sparlectra.validate_gui_config_overrides(conflict_overrides)
                @test_throws ArgumentError Sparlectra._load_api_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, conflict_nested)
            end

            @testset "the solver is always available" begin
                # AnalyticLoadFlow.jl is a required dependency: `using Sparlectra` is
                # enough, there is no session that has to load anything extra and no
                # not-installed error path left to take
                solver = Sparlectra.apslf_solver()
                @test solver isa Sparlectra.ApslfSolver
                @test solver isa Sparlectra.AbstractExternalSolver
                @test solver.order == 24 && solver.use_pade && !solver.nr_polish && solver.mode === :direct && solver.convergence_radius
                tuned = Sparlectra.apslf_solver(order=40, use_pade=false, nr_polish=false, mode=:outer, convergence_radius=false)
                @test (tuned.order, tuned.use_pade, tuned.nr_polish, tuned.mode, tuned.convergence_radius) == (40, false, false, :outer, false)
            end

            @testset "residual judged against the final active set (0.13.0)" begin
                # a PV machine with a band it cannot hold: AnalyticLoadFlow clamps it
                # at Qmax and reports the bus as PQ. The adapter used to test the
                # voltage setpoint of that bus anyway and labelled a solved case not
                # converged (case118: 19 clamped machines, 0.027 pu "mismatch").
                function _clamp_net()
                    net = Net(name="apslf_clamp", baseMVA=100.0)
                    addBus!(net=net, busName="B1", vn_kV=110.0)
                    addBus!(net=net, busName="B2", vn_kV=110.0)
                    addBus!(net=net, busName="B3", vn_kV=110.0)
                    addProsumer!(net=net, busName="B1", type="EXTERNALNETWORKINJECTION", vm_pu=1.02, va_deg=0.0, referencePri="B1")
                    addProsumer!(net=net, busName="B2", type="SYNCHRONOUSMACHINE", p=10.0, q=0.0, vm_pu=1.05, qMin=-10.0, qMax=10.0)
                    addProsumer!(net=net, busName="B3", type="ENERGYCONSUMER", p=80.0, q=25.0)
                    addPIModelACLine!(net=net, fromBus="B1", toBus="B2", r_pu=0.01, x_pu=0.08, b_pu=0.0, status=1)
                    addPIModelACLine!(net=net, fromBus="B2", toBus="B3", r_pu=0.02, x_pu=0.12, b_pu=0.0, status=1)
                    return net
                end
                net = _clamp_net()
                cfg = SparlectraConfig(powerflow=PowerFlowConfig(solver=:apslf, tol=1e-8), output=OutputConfig(logfile_results=:off))
                r = run_sparlectra(net=net, config=cfg)
                @test r.final_converged
                @test r.final_mismatch < 1e-8
                b2 = net.busDict["B2"]
                @test haskey(net.qLimitEvents, b2)
                @test net.qLimitEvents[b2] === :max
                st = Sparlectra.rectangular_pf_status(net)
                @test hasproperty(st, :apslf_convergence_radius) && isfinite(st.apslf_convergence_radius)
                @test occursin("dmin", String(st.apslf_convergence_line))
                # the header prints the radius next to the condition number
                hdr_path = tempname()
                open(hdr_path, "w") do io
                    redirect_stdout(io) do
                        printACPFlowResults(net, r.elapsed_s, r.iterations, 1e-8, false, ""; converged=r.final_converged, solver=:apslf)
                    end
                end
                @test occursin("APSLF radius   :", read(hdr_path, String))
                # switched off: the solve stays, the line says so
                net = _clamp_net()
                cfg_off = SparlectraConfig(powerflow=PowerFlowConfig(solver=:apslf, tol=1e-8, apslf=Sparlectra.ApslfConfig(convergence_radius=false)), output=OutputConfig(logfile_results=:off))
                r = run_sparlectra(net=net, config=cfg_off)
                @test r.final_converged
                @test occursin("not evaluated", String(Sparlectra.rectangular_pf_status(net).apslf_convergence_line))
            end

            @testset "Adapter mapping (PFModel -> AnalyticLoadFlow spec, PF ordering)" begin
                net = createTest3BusNet()
                model = buildPfModel(net; flatstart=true, include_limits=false)
                n = length(model.busIdx_net)

                spec = Sparlectra._apslf_spec_from_model(model)
                @test spec.Y === model.Ybus
                @test spec.bustype == model.busType
                @test spec.Pspec ≈ real.(model.Sspec)
                @test spec.Qspec ≈ imag.(model.Sspec)
                @test spec.Vm == model.Vset
                @test spec.Qmin == fill(-Inf, n)
                @test spec.Qmax == fill(Inf, n)
                @test spec.slack == model.slack_idx

                model_q = buildPfModel(net; flatstart=true, include_limits=true)
                spec_q = Sparlectra._apslf_spec_from_model(model_q)
                @test spec_q.Qmin == model_q.qmin_pu
                @test spec_q.Qmax == model_q.qmax_pu

                solver = Sparlectra.ApslfSolver(order=20, use_pade=true)
                sol = solvePf(solver, model)
                @test sol isa PFSolution
                @test length(sol.V) == n
                @test sol.meta.solver === :apslf
                @test sol.meta.mode isa Symbol
            end

            @testset "APSLF start-value generator guard and nr_polish=false in start mode" begin
                net = Sparlectra.importSCF(abspath(joinpath(dirname(@__DIR__), "data", "scf", "sp_case5.scf.json")))
                model = buildPfModel(net; flatstart=true, include_limits=false)

                # Disabled: no-op, restores the exact incoming vector.
                Vraw = copy(model.V0)
                V_disabled, summary_disabled = Sparlectra._run_guarded_apslf_start(model.Ybus, Vraw, model.Sspec, model.busType, model.Vset, model.slack_idx; enabled=false, order=40, baseMVA=model.baseMVA)
                @test !summary_disabled.apslf_start_enabled
                @test !summary_disabled.apslf_start_attempted
                @test V_disabled == Vraw

                # Use an already-converged NR solution as the incoming start: its own
                # mismatch sits at machine precision, so the raw (nr_polish=false)
                # APSLF candidate — which does not reach that precision without
                # polish — must always be rejected here, independent of any
                # solver/case-specific magnitude details (unlike comparing against
                # the flat start's own mismatch, which is not a stable invariant to
                # hard-code). This doubles as the nr_polish=false contract check: if
                # the generator used nr_polish=true internally, its candidate would
                # reach comparable machine-precision accuracy and could spuriously
                # "improve" on an already-converged start, defeating this test.
                net_converged = deepcopy(net)
                runpf!(net_converged, 30, 1e-8, 0)
                model_converged = buildPfModel(net_converged; flatstart=false, include_limits=false)
                Vgood = copy(model_converged.V0)
                @test Sparlectra.mismatchInf(model, Vgood) < 1e-6
                V_good, summary_good = Sparlectra._run_guarded_apslf_start(model.Ybus, Vgood, model.Sspec, model.busType, model.Vset, model.slack_idx; enabled=true, order=40, baseMVA=model.baseMVA)
                @test summary_good.apslf_start_enabled
                @test summary_good.apslf_start_attempted
                @test summary_good.apslf_start_accepted === false
                @test summary_good.apslf_start_reason === :not_improved
                @test V_good == Vgood

                # A deliberately perturbed, worse-than-APSLF start must be accepted.
                rng = Random.MersenneTwister(42)
                bad = model.V0 .* cis.(0.6 .* (rand(rng, length(model.V0)) .- 0.5))
                bad[model.slack_idx] = model.V0[model.slack_idx]
                V_bad, summary_bad = Sparlectra._run_guarded_apslf_start(model.Ybus, bad, model.Sspec, model.busType, model.Vset, model.slack_idx; enabled=true, order=40, baseMVA=model.baseMVA)
                @test summary_bad.apslf_start_accepted === true
                @test summary_bad.apslf_start_reason === :improved
                @test V_bad != bad
                @test all(v -> isfinite(real(v)) && isfinite(imag(v)), V_bad)
            end
        end
        @testset "APSLF agrees with NR without polish (ring3, sp_case5)" begin
            # The series alone is the solution; a Newton polish would hide a
            # wrong series result (seen on Windows / Julia 1.13.0, where the
            # series did not converge on cases with PV buses). Every APSLF test
            # runs without polish for that reason. Two shipped cases, the direct
            # adapter call and the framework run, each against NR to 1e-6 pu.
            # Never gate or relax.
            for (label, build) in (("ring3", ring3), ("sp_case5", () -> Sparlectra.importSCF(scf5)))
                r_nr = run_sparlectra(net=build(), config=cfg_nr)
                @test r_nr.final_converged
                vm_nr = getfield.(r_nr.net.nodeVec, :_vm_pu)
                va_nr = getfield.(r_nr.net.nodeVec, :_va_deg)

                net_direct = build()
                _, status_direct, sol_direct = runpf_external!(net_direct, apslf_solver(); tol=1e-8)
                @test status_direct == 0
                @test sol_direct.converged
                @test maximum(abs.(getfield.(net_direct.nodeVec, :_vm_pu) .- vm_nr)) < 1e-6
                @test maximum(abs.(getfield.(net_direct.nodeVec, :_va_deg) .- va_nr)) < 1e-4

                r_ap = run_sparlectra(net=build(), config=cfg_ap)
                @test r_ap.final_converged
                (status_direct == 0 && r_ap.final_converged) || _apslf_platform_probe()
                @test r_ap.diagnostics.solver === :apslf
                vm_ap = getfield.(r_ap.net.nodeVec, :_vm_pu)
                va_ap = getfield.(r_ap.net.nodeVec, :_va_deg)
                @test maximum(abs.(vm_nr .- vm_ap)) < 1e-6
                @test maximum(abs.(va_nr .- va_ap)) < 1e-4
                st = Sparlectra.rectangular_pf_status(r_ap.net)
                @test String(st.apslf_convergence_level) == "GRN"
                println("      APSLF against NR on ", label, ": max |dVm| ", maximum(abs.(vm_nr .- vm_ap)), " pu, ", st.apslf_convergence_line)
            end
        end

        @testset "APSLF workshop runs with its assertions" begin
            # the Literate workshop is executable Julia with an @assert next to every
            # printed number; running it here keeps the notebook from drifting. The
            # Q-limit section runs on the shipped data/mpower/sp_case118.m, so the
            # workshop needs no download and this test is no longer gated.
            workshop = abspath(joinpath(dirname(@__DIR__), "docs", "lit", "workshop_apslf.jl"))
            @test isfile(workshop)
            @test isfile(joinpath(dirname(@__DIR__), "data", "mpower", "sp_case118.m"))
            mod = Module(:WorkshopApslfRun)
            redirect_stdout(devnull) do
                Base.include(mod, workshop)
            end
            @test true
            println("      APSLF workshop: RAN")
        end

        return true
    end
end