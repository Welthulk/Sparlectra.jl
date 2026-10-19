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

# file: test/test_legacy_testgrid.jl
# purpose: legacy profile: the testsets of test/testgrid.jl that pin an earlier
#          change (an issue, a bug, a bit-identity or off-path guarantee).
#          They moved here unchanged from testgrid.jl; the helpers they use
#          are defined in that file, which the legacy group includes first.
using Sparlectra
using Test

# #375: the PQ->PV release decides on the voltage side. The shipped
# Zeng/Chiang 14-bus case (loads +25 percent, tight Q bands) ended on a
# non-physical solution under active_set: three machines at Qmax with the
# voltage above the setpoint, zero releases, the Q-V check reporting all
# three. With the voltage rule the run reaches the physical solution the
# one-at-a-time mode finds (bus 3 clamped at Qmax at 1.0064 pu, bus 1 at
# Qmin, bus 2 back on PV) and the Q-V check is clean. The feeder
# with a machine clamped at Qmax and its voltage BELOW the setpoint must
# stay clamped.
function test_active_set_voltage_side_release()::Bool
  @testset "Active set re-enable on the voltage side (#375)" begin (function ()
    scf_dir = joinpath(dirname(@__DIR__), "data", "scf")
    zeng = joinpath(scf_dir, "case14_zeng_p306_activeSet_A.scf.json")
    cfg = Sparlectra.resolve_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, zeng).config
    @test cfg.powerflow.qlimits.enforcement_mode === :active_set
    r = run_sparlectra(casefile = zeng, config = cfg)
    @test r.final_converged
    st = Sparlectra.rectangular_pf_status(r.net)
    @test st.qlimit_reenable_events >= 1
    # bus 2 (Qmax, voltage above its setpoint) is released; bus 3 stays at
    # Qmax with the voltage below the setpoint, bus 1 at Qmin
    @test !haskey(r.net.qLimitEvents, 2)
    @test get(r.net.qLimitEvents, 3, nothing) === :max && get(r.net.qLimitEvents, 1, nothing) === :min
    @test isapprox(r.net.nodeVec[3]._vm_pu, 1.0064; atol = 1e-3)
    @test isempty(Sparlectra.qvCharacteristicViolations(r.net))
    # clamped at Qmax with the voltage below the setpoint: no release
    feeder = joinpath(scf_dir, "feeder3_hardpv_A.scf.json")
    fcfg = Sparlectra.resolve_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, feeder).config
    rf = run_sparlectra(casefile = feeder, config = fcfg)
    @test rf.final_converged
    b2 = rf.net.busDict["B2"]
    @test get(rf.net.qLimitEvents, b2, nothing) === :max
    @test Sparlectra.rectangular_pf_status(rf.net).qlimit_reenable_events == 0
    # the margin is a configuration key that reaches the net
    @test rf.net.reenable_v_hyst_pu == fcfg.powerflow.qlimits.reenable_v_hyst_pu

    # The margin is GUI-editable since 0.20.2 and therefore a case-scope key:
    # the Zeng case carries a case configuration file, so a general file's
    # value no longer reaches it (the case is self-contained); another
    # margin comes as an explicit override, like any other form option.
    cfg_a = Sparlectra.resolve_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, zeng, Dict{String,Any}("power_flow.qlimits.reenable_v_hyst_pu" => 1.0)).config
    @test cfg_a.powerflow.qlimits.reenable_v_hyst_pu == 1.0
    ra = run_sparlectra(casefile = zeng, config = cfg_a)
    @test ra.final_converged
    @test Sparlectra.rectangular_pf_status(ra.net).qlimit_reenable_events == 0
    @test any(row -> row.significant, Sparlectra.qvCharacteristicViolations(ra.net))
    # run B (r above) lands on the voltages of the classic one-at-a-time mode
    cfg_c = Sparlectra.resolve_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, zeng, Dict{String,Any}("power_flow.qlimits.enforcement_mode" => "classic_one_at_a_time")).config
    rc = run_sparlectra(casefile = zeng, config = cfg_c)
    @test rc.final_converged
    @test isempty(Sparlectra.qvCharacteristicViolations(rc.net))
    @test isapprox(rc.net.nodeVec[3]._vm_pu, 1.0064; atol = 1e-3)
    # The final Q-limit check judges both modes alike: the classic run keeps
    # bus 6 on PV 0.68 MVAr over Qmax, inside the 1 MVAr hysteresis, and the
    # check says so instead of failing the run; the active-set run is
    # accepted with an ok or within_hysteresis verdict (up to 0.17.2 a
    # released machine inside the band failed it).
    st_c = Sparlectra.rectangular_pf_status(rc.net)
    @test st_c.final_q_check_status === :within_hysteresis
    row6 = only(filter(row -> row.busI == 6, st_c.final_q_check_rows))
    @test row6.side === :high && row6.class === :within_hysteresis
    @test isapprox(row6.dev_pu, 0.0068; atol = 5e-4)
    @test st.final_q_check_status in (:ok, :within_hysteresis)
    @test st.final_q_check_status !== :remaining_pv_q_limit_violations
    # with Q-limit handling off a converged run has nothing to check: the
    # status says so instead of claiming there was no converged solution
    cfg_d = Sparlectra.resolve_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, zeng, Dict{String,Any}("power_flow.qlimits.enabled" => false)).config
    rd = run_sparlectra(casefile = zeng, config = cfg_d)
    @test rd.final_converged
    @test Sparlectra.rectangular_pf_status(rd.net).final_q_check_status === :qlimits_disabled
    @test Sparlectra.final_q_check_line(:qlimits_disabled, NamedTuple[], 100.0) == "Final Q-limit check: not evaluated (Q-limit handling disabled)."
    @test Sparlectra.final_q_check_line(:not_evaluated, NamedTuple[], 100.0) == "Final Q-limit check: not evaluated (no converged solution)."
    # No bound on the voltage difference between the two modes on purpose:
    # switching variants need not end on the same point (measured:
    # 1.3e-3 pu apart at bus 6, the machine on its Qmax edge, clamped again by
    # the active set and kept on PV inside the hysteresis band by the classic
    # mode). What the test guards is that each mode ends on a PHYSICAL point:
    # the Q-V check above is clean for both.
    # regression on the shipped case118: the rule changes neither the clamped
    # set nor the iteration count, and nothing oscillates
    c118 = joinpath(dirname(@__DIR__), "data", "mpower", "sp_case118.m")
    for margin in (1.0, 1e-4)
      cfg118 = Sparlectra.resolve_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, c118, Dict{String,Any}("power_flow.qlimits.reenable_v_hyst_pu" => margin)).config
      r118 = run_sparlectra(casefile = c118, config = cfg118)
      @test r118.final_converged
      @test sort!(collect(keys(r118.net.qLimitEvents))) == [4, 19, 31, 32, 54, 72, 73, 77, 85, 87]
      @test r118.iterations <= 10
      @test Sparlectra.rectangular_pf_status(r118.net).oscillating_buses == 0
      @test Sparlectra.rectangular_pf_status(r118.net).final_q_check_status === :ok
    end
    # the feeder machine sits ON its clamp, so the final check has nothing to say
    @test Sparlectra.rectangular_pf_status(rf.net).final_q_check_status === :ok
  end)() end
  return true
end

# The final Q-limit check every enforcement mode ends with (task
# final-q-check): judged by the SIZE of the overshoot with two
# thresholds, the switching hysteresis and final_q_accept_pu. Synthetic
# injections on the 3-bus net, so every class is hit exactly, plus the
# configuration default and validation and the warning of the bounded class.
function test_final_q_limit_classification()::Bool
  @testset "Final Q-limit check: two thresholds" begin (function ()
    net = createTest3BusNet(cooldown = 1, hyst_pu = 0.01, qlim_min = -15.0, qlim_max = 15.0)
    pv = geNetBusIdx(net = net, busName = "STATION1")
    nb = length(net.nodeVec)
    bus_types = [Sparlectra.getNodeType(node) == Sparlectra.Slack ? :Slack : Sparlectra.getNodeType(node) == Sparlectra.PV ? :PV : :PQ for node in net.nodeVec]
    @test bus_types[pv] == :PV
    qmin_pu = fill(-Inf, nb)
    qmax_pu = fill(Inf, nb)
    qmin_pu[pv] = -0.15
    qmax_pu[pv] = 0.15
    qload = Sparlectra.build_qload_pu(net)
    # the injection that puts the machine at Qmax + dev (or Qmin - dev)
    injections(dev; side = :high) = [ComplexF64(0.0, i == pv ? ((side === :high ? qmax_pu[pv] + dev : qmin_pu[pv] - dev) - qload[pv]) : 0.0) for i in 1:nb]
    classify(dev; hyst = 0.01, accept = 0.02, side = :high) = classify_final_q_limits(net, injections(dev; side), bus_types, qmin_pu, qmax_pu; q_hyst_pu = hyst, final_q_accept_pu = accept, tol = 1e-8)
    # at the limit within the tolerance: no row, ok
    @test classify(0.0).status === :ok
    @test isempty(classify(0.0).rows)
    # inside the hysteresis: a row, accepted, class within_hysteresis
    within = classify(0.005)
    @test within.status === :within_hysteresis
    @test length(within.rows) == 1 && within.rows[1].class === :within_hysteresis && within.rows[1].side === :high
    @test isapprox(within.rows[1].dev_pu, 0.005; atol = 1e-12)
    @test within.max_dev_pu == within.rows[1].dev_pu
    @test occursin("bus $(within.rows[1].busI) over Qmax by 0.0050 pu / 0.50 MVAr (within_hysteresis)", final_q_check_line(within.status, within.rows, net.baseMVA))
    # beyond the hysteresis, inside the bound: bounded
    bounded = classify(0.015)
    @test bounded.status === :bounded_q_limit_violation
    @test bounded.bounded == 1 && bounded.violations == 0
    @test bounded.buses == string(bounded.rows[1].busI, ":bounded")
    # beyond the bound: a remaining violation, the run is not accepted
    @test classify(0.03).status === :remaining_pv_q_limit_violations
    @test classify(0.03).violations == 1
    # the mirror image under Qmin
    low = classify(0.015; side = :low)
    @test low.rows[1].side === :low && low.rows[1].class === :bounded
    # regression: both thresholds at zero is the strict check of 0.17.2
    @test classify(0.005; hyst = 0.0, accept = 0.0).status === :remaining_pv_q_limit_violations
    # a PQ bus is not judged, whatever its injection
    bus_types_pq = copy(bus_types)
    bus_types_pq[pv] = :PQ
    @test classify_final_q_limits(net, injections(0.03), bus_types_pq, qmin_pu, qmax_pu; q_hyst_pu = 0.01, final_q_accept_pu = 0.02, tol = 1e-8).status === :ok

    # the acceptance decision: within and bounded are accepted, bounded with
    # a warning naming bus and overshoot, a violation rejects the run
    finalize(dev) = Sparlectra._finalize_rectangular_qlimit_summary(net, ones(ComplexF64, nb), injections(dev), bus_types, qmin_pu, qmax_pu, true, :none;
      verbose = 0, qlimit_trace_enabled = false, q_hyst_pu = 0.01, final_q_accept_pu = 0.02, tol = 1e-8, pv_table_rows = 30,
      qlimit_guard_accept_bounded_violations = false, qlimit_guard_max_remaining_violations = 0)
    @test finalize(0.005).converged
    @test finalize(0.005).final_q_check.status === :within_hysteresis
    accepted = @test_logs (:warn, r"bounded violation accepted") finalize(0.015)
    @test accepted.converged
    @test accepted.rejection_reason === :bounded_q_limit_violation
    rejected = finalize(0.03)
    @test !rejected.converged
    @test rejected.rejection_reason === :remaining_pv_q_limit_violations

    # the key: default twice the hysteresis, never below it, stamped on the net
    cfg = Sparlectra.QLimitConfig(Dict{String,Any}("hysteresis_pu" => 0.01))
    @test cfg.final_q_accept_pu == 0.02
    @test Sparlectra.QLimitConfig(Dict{String,Any}("hysteresis_pu" => 0.0)).final_q_accept_pu == 0.0
    @test Sparlectra.QLimitConfig(Dict{String,Any}("hysteresis_pu" => 0.01, "final_q_accept_pu" => 0.05)).final_q_accept_pu == 0.05
    @test_throws ArgumentError Sparlectra.QLimitConfig(Dict{String,Any}("hysteresis_pu" => 0.01, "final_q_accept_pu" => 0.005))
    # auto (the template value) follows the hysteresis, so a user file that
    # raises hysteresis_pu alone is not rejected by an inherited literal
    @test Sparlectra.QLimitConfig(Dict{String,Any}("hysteresis_pu" => 0.03, "final_q_accept_pu" => "auto")).final_q_accept_pu == 0.06
    @test "power_flow.qlimits.final_q_accept_pu" in Sparlectra.GUI_EDITABLE_CONFIG_KEYS
    @test Net(name = "n", baseMVA = 100.0, q_hyst_pu = 0.03).final_q_accept_pu == 0.06

    # the same status for every mode at the same operating point: sp_case14
    # under the three enforcement modes
    case = joinpath(dirname(@__DIR__), "data", "scf", "sp_case14.scf.json")
    statuses = Symbol[]
    for mode in ("active_set", "classic_one_at_a_time", "classic_simultaneous")
      cfg14 = Sparlectra.resolve_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, case, Dict{String,Any}("power_flow.qlimits.enforcement_mode" => mode)).config
      r14 = run_sparlectra(casefile = case, config = cfg14)
      @test r14.final_converged
      push!(statuses, Sparlectra.rectangular_pf_status(r14.net).final_q_check_status)
    end
    @test length(unique(statuses)) == 1
    @test statuses[1] in (:ok, :within_hysteresis)
  end)() end
  return true
end

function run_legacy_testgrid_tests()
  # the call that sat in the block "Power flow scenarios", which stays in pf
  @test test_final_q_limit_classification() == true
  @testset "Active set voltage-side release" begin (function ()
    @test test_active_set_voltage_side_release() == true
  end)() end

  return true
end
