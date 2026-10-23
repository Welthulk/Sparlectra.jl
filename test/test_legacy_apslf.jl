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

# file: test/test_legacy_apslf.jl
# purpose: legacy profile: the testsets of test/test_apslf.jl that pin an earlier
#          change (an issue, a bug, a bit-identity or off-path guarantee).
#          They moved here unchanged from test_apslf.jl; the helpers they use
#          are defined in that file, which the legacy group includes first.
using Sparlectra
using Test

function run_legacy_apslf_tests()
  # cut from testset: APSLF (AnalyticLoadFlow.jl) integration
  (function ()
      # the helpers of the no-polish set and the platform probe live at
      # function scope
      quiet = OutputConfig(logfile_results=:off, console_summary=false, startup_latency_hint=false)
      cfg_ap = SparlectraConfig(powerflow=PowerFlowConfig(solver=:apslf), output=quiet)

              @testset "AnalyticLoadFlow outer-fallback bounds" begin (function ()
                  # The size bounds of the outer-mode fallback (_apslf_solve runs
                  # it itself) are mirrored from AnalyticLoadFlow 0.10.0, which
                  # provides no binding for them: they are literal keyword
                  # defaults of solve_pf_apslf in src/solver_core.jl. The loaded
                  # package's source must still carry the same defaults, and a
                  # binding named like a fallback bound must not have appeared;
                  # either change fails here, and the constants then have to be
                  # re-checked or read from the package.
                  mirrored_from = v"0.10.0"
                  @test pkgversion(AnalyticLoadFlow) >= mirrored_from
                  @test (Sparlectra._APSLF_OUTER_FALLBACK_NBUS_MAX, Sparlectra._APSLF_OUTER_FALLBACK_PV_BUS_MAX) == (2000, 250)
                  alf_core = read(joinpath(pkgdir(AnalyticLoadFlow), "src", "solver_core.jl"), String)
                  alf_default(key) = (m = match(Regex("$(key)\\) \\? Int\\(kw_forward\\.$(key)\\) : (\\d+)"), alf_core); m === nothing ? nothing : parse(Int, m[1]))
                  @test alf_default("outer_fallback_nbus_max") == Sparlectra._APSLF_OUTER_FALLBACK_NBUS_MAX
                  @test alf_default("outer_fallback_pv_bus_max") == Sparlectra._APSLF_OUTER_FALLBACK_PV_BUS_MAX
                  @test isempty(filter(n -> occursin("fallback", lowercase(String(n))), names(AnalyticLoadFlow; all = true)))
              end)() end

              @testset "residual judged against the final active set (0.13.0)" begin (function ()
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
                  # #461: the clamp counts as a PV->PQ switching event of the status
                  @test r.diagnostics.pv_pq_switching_events == 1
                  @test isfinite(st.apslf_series_radius) && occursin("series radius R", String(st.apslf_convergence_line))
                  # #461 on the shipped sp_case118 with every case bus number
                  # shifted by 1000: the radius buses are case bus numbers (an
                  # internal index 1..118 cannot pass), and the PV->PQ count of
                  # the status equals the clamps in the Q-limit log
                  mpc = Sparlectra.MatpowerIO.read_case(joinpath(dirname(@__DIR__), "data", "mpower", "sp_case118.m"))
                  mpc.bus[:, 1] .+= 1000
                  mpc.gen[:, 1] .+= 1000
                  mpc.branch[:, 1:2] .+= 1000
                  net118 = Sparlectra.createNetFromMatPowerCase(mpc = mpc)
                  r118 = run_sparlectra(net = net118, config = cfg_ap)
                  @test r118.final_converged
                  st118 = Sparlectra.rectangular_pf_status(net118)
                  @test st118.apslf_convergence_bus in 1001:1118
                  @test st118.apslf_series_radius_bus in 1001:1118
                  @test occursin(string("at bus ", st118.apslf_convergence_bus, ","), String(st118.apslf_convergence_line))
                  @test !isempty(net118.qLimitLog)
                  @test r118.diagnostics.pv_pq_switching_events == length(net118.qLimitLog)
                  # switched off: the solve stays, the line says so
                  net = _clamp_net()
                  cfg_off = SparlectraConfig(powerflow=PowerFlowConfig(solver=:apslf, tol=1e-8, apslf=Sparlectra.ApslfConfig(convergence_radius=false)), output=OutputConfig(logfile_results=:off))
                  r = run_sparlectra(net=net, config=cfg_off)
                  @test r.final_converged
                  @test occursin("not evaluated", String(Sparlectra.rectangular_pf_status(net).apslf_convergence_line))
              end)() end

              @testset "APSLF start-value generator guard and nr_polish=false in start mode" begin (function ()
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
              end)() end

  end)()
  return true
end
