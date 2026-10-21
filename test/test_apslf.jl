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
    # on a failure the adapter path is printed layer by layer, so a run on
    # another machine names the first layer whose numbers differ: the model as
    # runpf_external! builds it, AnalyticLoadFlow on exactly that spec, solvePf
    # on that model, runpf_external! on the net; each against Newton. stderr on
    # purpose: the quiet suite runner captures stdout and replays only a few
    # lines per failure block
    function _apslf_platform_probe()
        println(stderr, "      APSLF platform probe: Julia ", VERSION, " on ", Sys.KERNEL, " ", Sys.MACHINE, ", AnalyticLoadFlow ", pkgversion(AnalyticLoadFlow), ", BLAS ", BLAS.get_config())
        for (label, build) in (("ring3", ring3), ("sp_case5", () -> Sparlectra.importSCF(scf5)))
            net_nr = build()
            runpf!(net_nr, 30, 1e-10, 0)
            model = Sparlectra.buildPfModel(build(); include_limits=false)
            spec = Sparlectra._apslf_spec_from_model(model)
            vm_nr = [net_nr.nodeVec[i]._vm_pu for i in model.busIdx_net]
            println(stderr, "      probe ", label, " model: busType ", model.busType, ", slack ", model.slack_idx, ", Sspec ", model.Sspec, ", Vset ", model.Vset)
            println(stderr, "      probe ", label, " model: Ybus ", Matrix(model.Ybus))
            println(stderr, "      probe ", label, " Newton |V| in PF order: ", vm_nr)
            res = AnalyticLoadFlow.solve_pf_apslf(spec; mode=:direct, order=24, nr_polish=false, return_coeffs=true)
            println(stderr, "      probe ", label, " AnalyticLoadFlow on that spec: converged = ", res.converged, ", reason ", get(res, :reason, :none), ", outer_iters ", res.outer_iters, ", |V| ", abs.(res.V), ", Q ", res.Q, ", bustype ", res.bustype, ", ALF mismatch ", AnalyticLoadFlow.max_mismatch_on_specY(spec, res.V), ", max |dVm| vs NR ", maximum(abs.(abs.(res.V) .- vm_nr)))
            sol = Sparlectra.solvePf(apslf_solver(), model; tol=1e-8)
            println(stderr, "      probe ", label, " solvePf on that model: converged = ", sol.converged, ", residual_inf ", sol.residual_inf, ", |V| ", abs.(sol.V), ", max |dVm| vs NR ", maximum(abs.(abs.(sol.V) .- vm_nr)))
            net_ap = build()
            _, st, sol2 = runpf_external!(net_ap, apslf_solver(); tol=1e-8)
            println(stderr, "      probe ", label, " runpf_external!: status ", st, ", converged = ", sol2.converged, ", residual_inf ", sol2.residual_inf, ", net |V| ", getfield.(net_ap.nodeVec, :_vm_pu), ", max |dVm| vs NR ", maximum(abs.(getfield.(net_ap.nodeVec, :_vm_pu) .- getfield.(net_nr.nodeVec, :_vm_pu))))
        end
    end
    # SPARLECTRA_APSLF_PROBE=1 prints the probe on a passing run too
    get(ENV, "SPARLECTRA_APSLF_PROBE", "") == "1" && _apslf_platform_probe()
    # one outer set, so a standalone call still runs every group after a
    # failing one (a top-level testset throws at its end)
    @testset "APSLF" begin
          @testset "Config validation: power_flow.solver / apslf / apslf_start" begin (function ()
              default_cfg = Sparlectra.SparlectraConfig()
              @test default_cfg.powerflow.solver === :rectangular
              @test default_cfg.powerflow.apslf.order == 24
              @test default_cfg.powerflow.apslf.convergence_radius === true
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

              # An unknown key under the new sections is warned about and
              # dropped like any other file key; the file still loads (a
              # stored file must keep loading across releases, so an unknown
              # or removed key is noise, never a failure).
              bad_key_file = test_scratch_path(".yaml")
              write(bad_key_file, "power_flow:\n  apslf:\n    bogus_key: 1\n")
              # the file carries no config_version on purpose; that once-per-session
              # warning may appear as well and must not reach the test output
              cfg_bad = @test_logs (:warn, r"Unknown Sparlectra configuration key power_flow\.apslf\.bogus_key is ignored") match_mode = :any Sparlectra.load_sparlectra_config(bad_key_file; reload=true)
              @test cfg_bad isa Sparlectra.SparlectraConfig
              @test cfg_bad.powerflow.apslf.order == Sparlectra.SparlectraConfig().powerflow.apslf.order
          end)() end

          @testset "Controller + APSLF solver: runs without them, with a warning" begin (function ()
              # the shipped sp_case14 carries a REAL declared tap controller; an
              # APSLF run does not refuse it any more (0.31.0): it solves the
              # static setpoints and one warning names what is not applied
              net = Sparlectra.importSCF(abspath(joinpath(dirname(@__DIR__), "data", "scf", "sp_case14.scf.json")))
              @test length(Sparlectra.collect_outer_controllers(net)) == 1

              cfg = Sparlectra.SparlectraConfig(powerflow=Sparlectra.PowerFlowConfig(solver=:apslf), output=OutputConfig(logfile_results=:off, console_summary=false))
              res = @test_logs (:warn, r"power_flow.solver=apslf runs without controllers: 1 outer-loop controller") match_mode=:any redirect_stdout(devnull) do
                  run_sparlectra(net=net, config=cfg)
              end
              @test res.final_converged
          end)() end

          @testset "WebUI form parsing for APSLF fields -> effective config" begin (function ()
              form = Dict{String,Any}(
                  "casefile" => "case14.m",
                  "power_flow_solver" => "apslf",
                  "power_flow_apslf_order" => "25",
                  "power_flow_apslf_nr_polish" => "true",
                  "power_flow_apslf_start_enabled" => "false",
                  "power_flow_apslf_start_order" => "40",
              )
              request = SparlectraApp.powerflow_webui_request(form)
              overrides = request["config_overrides"]
              @test overrides["power_flow.solver"] == "apslf"
              @test overrides["power_flow.apslf.order"] === 25
              @test overrides["power_flow.apslf.nr_polish"] === true
              @test overrides["power_flow.apslf_start.enabled"] === false
              @test overrides["power_flow.apslf_start.order"] === 40

              nested = Sparlectra.validate_gui_config_overrides(overrides)
              cfg, _ = Sparlectra._load_api_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, nested)
              @test cfg.powerflow.solver === :apslf
              @test cfg.powerflow.apslf.order == 25
              @test cfg.powerflow.apslf.nr_polish === true

              # Unchecked checkbox explicitly submits "false" (hidden+checkbox HTML
              # pairing), the same semantics as the Q-limit checkbox pattern: the key
              # is present with an explicit false, not omitted.
              disabled_form = copy(form)
              disabled_form["power_flow_apslf_nr_polish"] = "false"
              @test SparlectraApp.powerflow_webui_request(disabled_form)["config_overrides"]["power_flow.apslf.nr_polish"] === false

              # A field genuinely absent from the submitted form (e.g. a stale
              # client) is skipped rather than defaulted to false.
              missing_field_form = copy(form)
              delete!(missing_field_form, "power_flow_apslf_start_enabled")
              @test !haskey(SparlectraApp.powerflow_webui_request(missing_field_form)["config_overrides"], "power_flow.apslf_start.enabled")

              # GUI-layer validation rejects the same invalid values as the core config.
              @test_throws ArgumentError Sparlectra.validate_gui_config_overrides(Dict{String,Any}("power_flow.solver" => "polar"))
              @test_throws ArgumentError Sparlectra.validate_gui_config_overrides(Dict{String,Any}("power_flow.apslf.order" => 0))

              # Downstream config assembly still rejects the solver+start conflict
              # even when both values arrive through the GUI override path.
              conflict_overrides = Dict{String,Any}("power_flow.solver" => "apslf", "power_flow.apslf_start.enabled" => true)
              conflict_nested = Sparlectra.validate_gui_config_overrides(conflict_overrides)
              @test_throws ArgumentError Sparlectra._load_api_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, conflict_nested)
          end)() end

          @testset "the solver is always available" begin (function ()
              # AnalyticLoadFlow.jl is a required dependency: `using Sparlectra` is
              # enough, there is no session that has to load anything extra and no
              # not-installed error path left to take
              solver = Sparlectra.apslf_solver()
              @test solver isa Sparlectra.ApslfSolver
              @test solver isa Sparlectra.AbstractExternalSolver
              @test solver.order == 24 && !solver.nr_polish && solver.mode === :direct && solver.convergence_radius
              tuned = Sparlectra.apslf_solver(order=40, nr_polish=false, mode=:outer, convergence_radius=false)
              @test (tuned.order, tuned.nr_polish, tuned.mode, tuned.convergence_radius) == (40, false, :outer, false)
          end)() end

          @testset "Adapter mapping (PFModel -> AnalyticLoadFlow spec, PF ordering)" begin (function ()
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

              solver = Sparlectra.ApslfSolver(order=20)
              sol = solvePf(solver, model)
              @test sol isa PFSolution
              @test length(sol.V) == n
              @test sol.meta.solver === :apslf
              @test sol.meta.mode isa Symbol
          end)() end

        @testset "APSLF agrees with NR without polish (ring3, sp_case5)" begin (function ()
            # The series alone is the solution; a Newton polish would hide a
            # wrong series result (AnalyticLoadFlow 0.9.14 returns one flagged
            # as converged, 2e-2 pu off on ring3). Every APSLF test
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
                # a converged run carries no diagnosis
                @test !hasproperty(st, :apslf_diagnosis_reason)
                println("      APSLF against NR on ", label, ": max |dVm| ", maximum(abs.(vm_nr .- vm_ap)), " pu, ", st.apslf_convergence_line)
            end
            # A non-converged run carries AnalyticLoadFlow's diagnosis (0.10.2)
            # into the status: sp_case5 at order 4 (Q limits off) stops with
            # the Pade value still moving by 1.4e-3 pu between order 2 and 4,
            # converged from order 8 on. The convergence decision is the
            # package's; only the reason is checked.
            cfg_low = SparlectraConfig(powerflow=PowerFlowConfig(solver=:apslf, apslf=Sparlectra.ApslfConfig(order=4), qlimits=Sparlectra.QLimitConfig(ignore_q_limits=true)), output=quiet)
            r_low = run_sparlectra(net=Sparlectra.importSCF(scf5), config=cfg_low)
            @test !r_low.final_converged
            st_low = Sparlectra.rectangular_pf_status(r_low.net)
            @test st_low.apslf_diagnosis_reason === :pade_not_converged
            @test isfinite(st_low.apslf_pole_min_abs) && st_low.apslf_pole_min_bus >= 1
            @test occursin(String(st_low.apslf_diagnosis_line), r_low.reason_text)
            # The same island next to a small second island, run with an output
            # directory as the service does: the island-wise error message and
            # the island's solver log carry the diagnosis of the failed island,
            # and the service reads its keys for the run metadata.
            net_isl = Sparlectra.importSCF(scf5)
            for b in ("X1", "X2")
                addBus!(net=net_isl, busName=b, vn_kV=110.0)
            end
            addProsumer!(net=net_isl, busName="X1", type="EXTERNALNETWORKINJECTION", vm_pu=1.0, va_deg=0.0, referencePri="X1")
            addProsumer!(net=net_isl, busName="X2", type="ENERGYCONSUMER", p=1.0, q=0.2)
            addPIModelACLine!(net=net_isl, fromBus="X1", toBus="X2", r_pu=0.001, x_pu=0.01, b_pu=0.0, status=1)
            profile = Dict{Symbol,Any}(:output_dir => mktempdir())
            err_isl = try
                redirect_stdout(devnull) do
                    run_sparlectra(net=net_isl, config=cfg_low, performance_profile=profile)
                end
                nothing
            catch err
                err isa ErrorException || rethrow()
                err
            end
            st_isl = Sparlectra._islandwise_apslf_diagnosis(profile)
            @test st_isl !== nothing && st_isl.apslf_diagnosis_reason === :pade_not_converged
            @test err_isl !== nothing && occursin(String(st_isl.apslf_diagnosis_line), err_isl.msg)
            @test occursin(String(st_isl.apslf_diagnosis_line), read(joinpath(profile[:output_dir], "ac_island_$(st_isl.island_id)_solver.log"), String))
            # A fixed phase shift is modelled through the unsymmetric Y-bus
            # Sparlectra hands to AnalyticLoadFlow: the shipped sp_casePST
            # with its shifter at -20 degrees agrees with NR with the
            # package's own default germ (no pin since 0.30.3, #470).
            scf_pst = joinpath(dirname(@__DIR__), "data", "scf", "sp_casePST.scf.json")
            pst20() = (n = Sparlectra.importSCF(scf_pst); foreach(b -> b.phase_shift_deg != 0 && (b.phase_shift_deg = -20.0), n.branchVec); n)
            r_nr_pst = run_sparlectra(net=pst20(), config=cfg_nr)
            r_ap_pst = run_sparlectra(net=pst20(), config=cfg_ap)
            @test r_nr_pst.final_converged
            @test r_ap_pst.final_converged
            @test maximum(abs.(getfield.(r_nr_pst.net.nodeVec, :_vm_pu) .- getfield.(r_ap_pst.net.nodeVec, :_vm_pu))) < 1e-6
            @test maximum(abs.(getfield.(r_nr_pst.net.nodeVec, :_va_deg) .- getfield.(r_ap_pst.net.nodeVec, :_va_deg))) < 1e-4
            # 0.30.2: the slack and PV voltage come from the generator, not
            # from the bus. A net built without a bus vm_pu used to solve
            # APSLF with the slack at 1.0 pu while Newton held the
            # generator's 1.02 pu.
            function ring3_gen_vm()
                net = Net(name="apslf_ring3_gen_vm", baseMVA=100.0)
                addBus!(net=net, busName="B1", vn_kV=110.0)
                addBus!(net=net, busName="B2", vn_kV=110.0)
                addBus!(net=net, busName="B3", vn_kV=110.0)
                addPIModelACLine!(net=net, fromBus="B1", toBus="B2", r_pu=0.010, x_pu=0.080, b_pu=0.0, status=1)
                addPIModelACLine!(net=net, fromBus="B2", toBus="B3", r_pu=0.011, x_pu=0.085, b_pu=0.0, status=1)
                addPIModelACLine!(net=net, fromBus="B3", toBus="B1", r_pu=0.012, x_pu=0.090, b_pu=0.0, status=1)
                addProsumer!(net=net, busName="B1", type="EXTERNALNETWORKINJECTION", referencePri="B1", vm_pu=1.02, va_deg=0.0)
                addProsumer!(net=net, busName="B2", type="SYNCHRONOUSMACHINE", p=20.0, vm_pu=1.03)
                addProsumer!(net=net, busName="B3", type="LOAD", p=30.0, q=10.0)
                ok, msg = validate!(net=net)
                ok || error("ring3_gen_vm invalid: $msg")
                return net
            end
            @test Sparlectra.buildPfModel(ring3_gen_vm(); include_limits=false).Vset[1:2] == [1.02, 1.03]
            for cfg in (cfg_nr, cfg_ap)
                r = run_sparlectra(net=ring3_gen_vm(), config=cfg)
                @test r.final_converged
                @test r.net.nodeVec[1]._vm_pu ≈ 1.02 atol = 1e-10
                @test r.net.nodeVec[2]._vm_pu ≈ 1.03 atol = 1e-8
            end
            # 0.60.0: the slack angle reaches AnalyticLoadFlow (Vslack). With
            # the sp_case5 slack at 10 degrees, APSLF fixed the slack at 0 and
            # returned the solution rotated by 10 degrees (0.178 pu in
            # complex V); through the start generator Newton then kept the
            # rotated slack. Solver, framework run and start path are
            # compared with Newton in complex V.
            function case5_va10()
                net = Sparlectra.importSCF(scf5)
                foreach(n -> Sparlectra.isSlack(n) && (n._va_deg = 10.0), net.nodeVec)
                return net
            end
            cvec(net) = [n._vm_pu * cis(deg2rad(n._va_deg)) for n in net.nodeVec]
            slack5 = findfirst(Sparlectra.isSlack, case5_va10().nodeVec)
            r_nr10 = run_sparlectra(net=case5_va10(), config=cfg_nr)
            @test r_nr10.final_converged
            V_nr10 = cvec(r_nr10.net)
            @test rad2deg(angle(V_nr10[slack5])) ≈ 10.0 atol = 1e-12
            net_d10 = case5_va10()
            _, status_d10, _ = runpf_external!(net_d10, apslf_solver(); tol=1e-8)
            @test status_d10 == 0
            @test maximum(abs.(cvec(net_d10) .- V_nr10)) < 1e-5
            r_ap10 = run_sparlectra(net=case5_va10(), config=cfg_ap)
            @test r_ap10.final_converged
            @test maximum(abs.(cvec(r_ap10.net) .- V_nr10)) < 1e-5
            cfg_st = SparlectraConfig(powerflow=PowerFlowConfig(solver=:rectangular, rescue=false, apslf_start=Sparlectra.ApslfStartConfig(enabled=true, order=40)), output=quiet)
            r_st10 = run_sparlectra(net=case5_va10(), config=cfg_st)
            @test r_st10.final_converged
            @test maximum(abs.(cvec(r_st10.net) .- V_nr10)) < 1e-5
            m10 = Sparlectra.buildPfModel(case5_va10(); include_limits=false)
            Vg10, sg10 = Sparlectra._run_guarded_apslf_start(m10.Ybus, m10.V0, m10.Sspec, m10.busType, m10.Vset, m10.slack_idx; enabled=true, order=40, baseMVA=m10.baseMVA)
            @test sg10.apslf_start_accepted
            @test maximum(abs.(Vg10 .- V_nr10[m10.busIdx_net])) < 1e-5
        end)() end


        return true
    end
end