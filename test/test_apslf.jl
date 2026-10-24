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
# purpose: tests the APSLF (AnalyticLoadFlow.jl) solver integration: a run
#          with a declared controller, and agreement with Newton-Raphson
#          without polish on small shipped cases, the reactive limits
#          (clamp set and band check against Newton on sp_case118, both
#          sides) and the injection source (node-only edits move no solver)
#
# AnalyticLoadFlow.jl is a required dependency, so nothing here is
# conditional. The config value domains of power_flow.apslf live in
# test_configuration_coverage.jl, the form parsing in test_webui.jl.

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

        @testset "APSLF reactive limits: same clamps as Newton, band checked (sp_case118)" begin (function ()
            # 0.60.0: AnalyticLoadFlow compares the bus injection with the
            # band, and buildPfModel handed it the generator band. On
            # sp_case118 the PV machines at buses 4 and 85 (Q load 23 and 11
            # MVAr) ended above Qmax, other machines were clamped away from
            # their limit, and the run reported the limits as ok. With the
            # band in injection terms APSLF ends on Newton's clamp set and
            # voltages; the Q-limit stability check of both runs passes.
            case118 = joinpath(dirname(@__DIR__), "data", "mpower", "sp_case118.m")
            r_nr, r_ap = redirect_stdout(devnull) do
                (run_sparlectra(casefile=basename(case118), path=dirname(case118), config=cfg_nr),
                 run_sparlectra(casefile=basename(case118), path=dirname(case118), config=cfg_ap))
            end
            clamps(r) = sort!(collect(keys(r.net.qLimitEvents)))
            @test r_nr.final_converged && r_ap.final_converged
            @test clamps(r_ap) == clamps(r_nr)
            @test 4 in clamps(r_ap) && 85 in clamps(r_ap)
            @test r_ap.limit_validation_status === :ok
            @test Sparlectra.rectangular_pf_status(r_ap.net).q_limit_active_set_ok
            @test r_ap.diagnostics.qlimit_stability_status === :ok
            # same clamp set, same solution: measured 1.8e-15 pu
            @test maximum(abs.(getfield.(r_ap.net.nodeVec, :_vm_pu) .- getfield.(r_nr.net.nodeVec, :_vm_pu))) <= 1e-12
            # The limit status is checked on the solved net, not taken from
            # the convergence flag: the same case solved with the generator
            # band (the fix bypassed in a local model) fails the check with
            # buses 4 and 85 among the machines outside their band.
            net = Sparlectra.createNetFromMatPowerFile(filename=case118)
            m = Sparlectra.buildPfModel(net)
            ql = [something(net.nodeVec[b]._qƩLoad, 0.0) / net.baseMVA for b in m.busIdx_net]
            m_gen = Sparlectra.PFModel(Ybus=m.Ybus, baseMVA=m.baseMVA, busIdx_net=m.busIdx_net, busType=m.busType, slack_idx=m.slack_idx, Vset=m.Vset, Sspec=m.Sspec, V0=m.V0, qmin_pu=m.qmin_pu .+ ql, qmax_pu=m.qmax_pu .+ ql)
            sol = Sparlectra.solvePf(apslf_solver(), m_gen)
            @test sol.converged
            Sparlectra.applyPfSolution!(net, m_gen, sol)
            Sparlectra._apslf_record_clamps!(net, sol, sol.iters)
            band = Sparlectra._apslf_qband_status(net, true, true)
            st_bad = Sparlectra._apslf_pf_status(true, sol.residual_inf, sol.iters; extra=band)
            @test !st_bad.q_limit_active_set_ok
            @test Sparlectra._rectangular_run_status(st_bad).limit_validation_status === :fail
            @test occursin("4 (PV outside its band", st_bad.reason_text) && occursin("85 (PV outside its band", st_bad.reason_text)
            # A fixed-Q generator next to the regulating machine: its Q is in
            # the bus band at both ends (buildQLimits!), the node Q load is
            # the only shift. B2: machine band -5..5 MVAr, fixed 3 MVAr, load
            # 20 MVAr; the voltage set point 1.03 pu needs more than 8 MVAr,
            # so both solvers clamp the bus at 5 + 3 = 8 MVAr generation.
            function ring3_fixed_q()
                net = Net(name="apslf_ring3_fixed_q", baseMVA=100.0)
                foreach(b -> addBus!(net=net, busName=b, vn_kV=110.0), ("B1", "B2", "B3"))
                addPIModelACLine!(net=net, fromBus="B1", toBus="B2", r_pu=0.010, x_pu=0.080, b_pu=0.0, status=1)
                addPIModelACLine!(net=net, fromBus="B2", toBus="B3", r_pu=0.011, x_pu=0.085, b_pu=0.0, status=1)
                addPIModelACLine!(net=net, fromBus="B3", toBus="B1", r_pu=0.012, x_pu=0.090, b_pu=0.0, status=1)
                addProsumer!(net=net, busName="B1", type="EXTERNALNETWORKINJECTION", referencePri="B1", vm_pu=1.0, va_deg=0.0)
                addProsumer!(net=net, busName="B2", type="SYNCHRONOUSMACHINE", p=20.0, vm_pu=1.03, qMin=-5.0, qMax=5.0)
                addProsumer!(net=net, busName="B2", type="GENERATOR", p=5.0, q=3.0)
                addProsumer!(net=net, busName="B2", type="LOAD", p=10.0, q=20.0)
                addProsumer!(net=net, busName="B3", type="LOAD", p=30.0, q=10.0)
                return net
            end
            r3_nr = run_sparlectra(net=ring3_fixed_q(), config=cfg_nr)
            r3_ap = run_sparlectra(net=ring3_fixed_q(), config=cfg_ap)
            for r in (r3_nr, r3_ap)
                @test r.final_converged && r.limit_validation_status === :ok
                @test get(r.net.qLimitEvents, 2, :none) === :max
                @test r.net.nodeVec[2]._qƩGen ≈ 8.0 atol = 1e-6
            end
            @test abs(r3_ap.net.nodeVec[2]._vm_pu - r3_nr.net.nodeVec[2]._vm_pu) < 1e-6
            # The Qmin side: the machine at bus 46 (Q load 21 MVAr) runs at
            # -1.26 MVAr under Newton. With its band narrowed to 5..10 MVAr
            # it falls below Qmin; both solvers clamp it at 5 MVAr, with
            # the same clamp set and voltages and the limit status ok.
            # Measured 8.9e-16 pu.
            name_of(net, bus) = only(k for (k, v) in net.busDict if v == bus)
            function case118_qmin()
                net = Sparlectra.createNetFromMatPowerFile(filename=case118)
                Sparlectra.setQLimits!(net=net, qmin_MVar=5.0, qmax_MVar=10.0, busName=name_of(net, 46))
                return net
            end
            rq_nr, rq_ap = redirect_stdout(devnull) do
                (run_sparlectra(net=case118_qmin(), config=cfg_nr), run_sparlectra(net=case118_qmin(), config=cfg_ap))
            end
            for r in (rq_nr, rq_ap)
                @test r.final_converged && r.limit_validation_status === :ok
                @test get(r.net.qLimitEvents, 46, :none) === :min
                @test r.net.nodeVec[46]._qƩGen ≈ 5.0 atol = 1e-6
            end
            @test clamps(rq_ap) == clamps(rq_nr)
            @test rq_ap.diagnostics.qlimit_stability_status === :ok
            @test maximum(abs.(getfield.(rq_ap.net.nodeVec, :_vm_pu) .- getfield.(rq_nr.net.nodeVec, :_vm_pu))) <= 1e-12
            # One source of the injections: buildPfModel builds them from the
            # prosumers like Newton (buildComplexSVec), not from the node
            # sums. A node-only edit (addBusLoadPower!, read by no solver)
            # moved APSLF and nothing else: ring3 with 4 MVAr added at B2
            # that way gave Vm 0.988093 pu against Newton's 0.990335. Here
            # node-only loads at the clamped PV bus 4 and the PQ bus 2 leave
            # both solvers on the result of the unedited case.
            function case118_nodeonly()
                net = Sparlectra.createNetFromMatPowerFile(filename=case118)
                addBusLoadPower!(net=net, busName=name_of(net, 4), p=10.0, q=10.0)
                addBusLoadPower!(net=net, busName=name_of(net, 2), p=10.0, q=10.0)
                return net
            end
            rn_nr, rn_ap = redirect_stdout(devnull) do
                (run_sparlectra(net=case118_nodeonly(), config=cfg_nr), run_sparlectra(net=case118_nodeonly(), config=cfg_ap))
            end
            @test rn_nr.final_converged && rn_ap.final_converged && rn_ap.limit_validation_status === :ok
            @test clamps(rn_ap) == clamps(rn_nr) == clamps(r_nr)
            @test maximum(abs.(getfield.(rn_ap.net.nodeVec, :_vm_pu) .- getfield.(rn_nr.net.nodeVec, :_vm_pu))) <= 1e-12
            @test maximum(abs.(getfield.(rn_ap.net.nodeVec, :_vm_pu) .- getfield.(r_ap.net.nodeVec, :_vm_pu))) <= 1e-12
            # the write-back adds the load the solve saw: slack generation as Newton's
            slack = findfirst(Sparlectra.isSlack, rn_ap.net.nodeVec)
            @test rn_ap.net.nodeVec[slack]._pƩGen ≈ rn_nr.net.nodeVec[slack]._pƩGen atol = 1e-6
            @test rn_ap.net.nodeVec[slack]._qƩGen ≈ rn_nr.net.nodeVec[slack]._qƩGen atol = 1e-6
        end)() end


        return true
    end
end