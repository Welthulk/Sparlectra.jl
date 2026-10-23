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

# file: test/test_tap_controller.jl
# purpose: tests the control framework and tap controllers: voltage and
#          branch active-power tap control, PST reactance coupling, SVC
#          shunts, machine RVC, YAML instantiation, and run_control!

# A user-defined controller that switches itself off when it is initialized,
# as one whose target cannot be found would; run_control! then ends with
# :no_active_controllers after the first solve.
mutable struct _SelfDisablingTestController <: Sparlectra.AbstractOuterController
  enabled::Bool
end
function Sparlectra.control_initialize!(ctrl::_SelfDisablingTestController, ::Any, context)
  ctrl.enabled = false
  return Sparlectra.NoControlState()
end

function run_tap_controller_tests()
  _runner_cfg(; max_iter::Int = 30, tol::Float64 = 1e-9, control::ControlConfig = ControlConfig()) = SparlectraConfig(powerflow = PowerFlowConfig(max_iter = max_iter, tol = tol), output = OutputConfig(logfile_results = :off), control = control)

  function _build_net()
    net = Net(name = "tap_ctrl", baseMVA = 100.0)
    addBus!(net = net, busName = "Slack", vn_kV = 110.0)
    addBus!(net = net, busName = "Mid", vn_kV = 110.0)
    addBus!(net = net, busName = "Load", vn_kV = 110.0)

    addProsumer!(net = net, busName = "Slack", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.02, va_deg = 0.0, referencePri = "Slack")
    addProsumer!(net = net, busName = "Load", type = "ENERGYCONSUMER", p = -70.0, q = -20.0)

    addPIModelTrafo!(net = net, fromBus = "Slack", toBus = "Mid", r_pu = 0.01, x_pu = 0.08, b_pu = 0.00, ratio = 1.0, shift_deg = 0.0, status = 1)
    addPIModelACLine!(net = net, fromBus = "Mid", toBus = "Load", r_pu = 0.02, x_pu = 0.12, b_pu = 0.01, status = 1)

    tbr = getNetBranch(net = net, fromBus = "Slack", toBus = "Mid")
    tbr.has_ratio_tap = true
    tbr.has_phase_tap = true
    tbr.tap_min = 0.95
    tbr.tap_max = 1.05
    tbr.tap_step = 0.0125
    tbr.phase_min_deg = -10.0
    tbr.phase_max_deg = 10.0
    tbr.phase_step_deg = 1.0

    return net, tbr
  end

  # _build_net plus a second AC island with its own slack, load and line:
  # two islands that both have branches
  function _build_two_island_net()
    net, tbr = _build_net()
    for bus in ("Slack2", "Load2")
      addBus!(net = net, busName = bus, vn_kV = 110.0)
    end
    addProsumer!(net = net, busName = "Slack2", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "Slack2")
    addProsumer!(net = net, busName = "Load2", type = "ENERGYCONSUMER", p = -30.0, q = -10.0)
    addPIModelACLine!(net = net, fromBus = "Slack2", toBus = "Load2", r_pu = 0.02, x_pu = 0.12, b_pu = 0.01, status = 1)
    return net, tbr
  end

  @testset "Tap controller API validation" begin (function ()
    net, _ = _build_net()
    @test_throws ErrorException addPowerTransformerControl!(net; trafo = "does_not_exist", mode = :voltage, target_bus = "Load", target_vm_pu = 1.0)
    @test_throws ErrorException addPowerTransformerControl!(net; trafo = "1", mode = :voltage)
    @test_throws ErrorException addPowerTransformerControl!(net; trafo = "1", mode = :voltage, target_bus = "Load", target_vm_pu = 1.0, voltage_error_metric = :invalid)
    # one controller per transformer, whichever alias names it (branch index or cID)
    anet, atbr = _build_net()
    addPowerTransformerControl!(anet; trafo = string(atbr.branchIdx), mode = :voltage, target_bus = "Load", target_vm_pu = 0.99, control_ratio = true, control_phase = false)
    @test_throws ErrorException addPowerTransformerControl!(anet; trafo = atbr.comp.cID, mode = :voltage, target_bus = "Load", target_vm_pu = 0.99, control_ratio = true, control_phase = false)
  end)() end

  @testset "Split Schraegregelung: two disjoint controllers on one transformer" begin (function ()
    # Voltage controller on the ratio tap plus active-power controller on the
    # phase tap of the SAME transformer; a parallel line provides the loop the
    # phase tap needs to shift flow.
    net = Net(name = "schraeg_split", baseMVA = 100.0)
    for bus in ("Slack", "Mid", "Load")
      addBus!(net = net, busName = bus, vn_kV = 110.0)
    end
    addProsumer!(net = net, busName = "Slack", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.02, va_deg = 0.0, referencePri = "Slack")
    addProsumer!(net = net, busName = "Load", type = "ENERGYCONSUMER", p = -70.0, q = -20.0)
    addPIModelTrafo!(net = net, fromBus = "Slack", toBus = "Mid", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, ratio = 1.0, shift_deg = 0.0, status = 1)
    addPIModelACLine!(net = net, fromBus = "Slack", toBus = "Mid", r_pu = 0.03, x_pu = 0.2, b_pu = 0.0, status = 1)
    addPIModelACLine!(net = net, fromBus = "Mid", toBus = "Load", r_pu = 0.02, x_pu = 0.12, b_pu = 0.01, status = 1)

    tbr = getNetBranch(net = net, fromBus = "Slack", toBus = "Mid")
    tbr.has_ratio_tap = true
    tbr.tap_min = 0.95
    tbr.tap_max = 1.08
    tbr.tap_step = 0.0125
    tbr.has_phase_tap = true
    tbr.phase_min_deg = -10.0
    tbr.phase_max_deg = 10.0
    tbr.phase_step_deg = 0.5

    result0 = run_sparlectra(net = net, config = _runner_cfg())
    @test result0.numerical_converged
    p0 = get_branch_p_from_to_mw(net, "Slack", "Mid")

    addPowerTransformerControl!(net; trafo = string(tbr.branchIdx), mode = :voltage, target_bus = "Load", target_vm_pu = 1.03, control_ratio = true, control_phase = false, deadband_vm_pu = 5e-3)
    addPowerTransformerControl!(net; trafo = string(tbr.branchIdx), mode = :branch_active_power, target_branch = ("Slack", "Mid"), p_target_mw = round(p0) + 8.0, control_ratio = false, control_phase = true, deadband_p_mw = 2.0)
    @test length(Sparlectra._tap_controllers(net)) == 2

    # Per-actuator exclusivity: both actuators are taken now.
    @test_throws ErrorException addPowerTransformerControl!(net; trafo = string(tbr.branchIdx), mode = :voltage, target_bus = "Load", target_vm_pu = 1.0)
    @test_throws ErrorException addPowerTransformerControl!(net; trafo = string(tbr.branchIdx), mode = :branch_active_power, target_branch = ("Slack", "Mid"), p_target_mw = 0.0, control_ratio = false, control_phase = true)

    result = run_sparlectra(net = net, config = _runner_cfg())
    @test result.numerical_converged
    ctrls = Sparlectra._tap_controllers(net)
    @test all(c -> c.converged, ctrls)
    @test abs(get_bus_vm_pu(net, "Load") - 1.03) <= 5e-3
    @test abs(get_branch_p_from_to_mw(net, "Slack", "Mid") - (round(p0) + 8.0)) <= 2.0

    # Report rows: one row per channel, no duplication across controllers.
    cres = latest_control_result(net)
    @test cres !== nothing
    @test length(cres.controllers) == 2
    @test length(unique([row.controller_name for row in cres.controllers])) == 2
  end)() end

  # typed model on the winding (0.20.0): the branch shows the model's shift
  # at its step, the build records the precedence line, and a phase
  # controller moves the model step (5 instead of 3) so the solved flow
  # matches a branch built directly at step 5
  @testset "typed phase model drives the branch and the controller moves its step" begin (function ()
    build_model_net = function (step::Int)
      mnet = Net(name = "pst_model_step$(step)", baseMVA = 100.0)
      for bus in ("Slack", "Mid", "Load")
        addBus!(net = mnet, busName = bus, vn_kV = 110.0)
      end
      addProsumer!(net = mnet, busName = "Slack", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.02, va_deg = 0.0, referencePri = "Slack")
      addProsumer!(net = mnet, busName = "Load", type = "ENERGYCONSUMER", p = -70.0, q = -20.0)
      model = PhaseTapChangerModel(kind = :asymmetrical, step = step, lowStep = -10, highStep = 10, neutralStep = 0, voltage_step_increment = 0.01, winding_connection_angle_deg = 90.0)
      addPIModelTrafo!(net = mnet, fromBus = "Slack", toBus = "Mid", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, ratio = 1.0, shift_deg = 0.0, status = 1, phase_taps = model)
      addPIModelACLine!(net = mnet, fromBus = "Slack", toBus = "Mid", r_pu = 0.03, x_pu = 0.2, b_pu = 0.0, status = 1)
      addPIModelACLine!(net = mnet, fromBus = "Mid", toBus = "Load", r_pu = 0.02, x_pu = 0.12, b_pu = 0.01, status = 1)
      return mnet, getNetBranch(net = mnet, fromBus = "Slack", toBus = "Mid"), model
    end
    net3, br3, m3 = build_model_net(3)
    expect3 = calcPhaseTapAngleRatio(m3)
    @test br3.taps_derived
    @test isapprox(br3.angle, expect3.effective_shift_deg; atol = 1e-12)
    @test isapprox(br3.ratio, expect3.effective_ratio; atol = 1e-12)
    @test br3.phase_shift_deg == br3.angle
    @test any(n -> occursin("replaced by model at phase step 3", n), net3.tapModelNotices)
    # a winding without a model leaves the branch on its explicit values
    @test !getNetBranch(net = net3, fromBus = "Mid", toBus = "Load").taps_derived
    # the flow a direct build at step 5 gives
    net5, br5, m5 = build_model_net(5)
    r5 = run_sparlectra(net = net5, config = _runner_cfg())
    @test r5.numerical_converged
    p5 = get_branch_p_from_to_mw(net5, "Slack", "Mid")
    r3 = run_sparlectra(net = net3, config = _runner_cfg())
    @test r3.numerical_converged
    p3 = get_branch_p_from_to_mw(net3, "Slack", "Mid")
    @test abs(p5 - p3) > 0.5
    # the controller moves the model from step 3 toward the step-5 flow
    addPowerTransformerControl!(net3; trafo = string(br3.branchIdx), mode = :branch_active_power, target_branch = ("Slack", "Mid"), p_target_mw = p5, control_ratio = false, control_phase = true, deadband_p_mw = 0.25 * abs(p5 - p3))
    rc = run_sparlectra(net = net3, config = _runner_cfg())
    @test rc.numerical_converged
    @test m3.step == 5
    @test isapprox(br3.angle, calcPhaseTapAngleRatio(m3).effective_shift_deg; atol = 1e-12)
    @test abs(get_branch_p_from_to_mw(net3, "Slack", "Mid") - p5) <= 1e-6 + 0.25 * abs(p5 - p3)
    rows = latest_control_result(net3).controllers
    @test length(rows) == 1
    @test rows[1].phase_tap_position == 5
  end)() end

  @testset "Voltage deadband is evaluated in pu Vm space" begin (function ()
    @test Sparlectra._voltage_within_deadband(1.2009, 1.200, 1e-3)
    @test !Sparlectra._voltage_within_deadband(1.2025, 1.200, 1e-3)
    @test isapprox(Sparlectra._voltage_control_error(1.201, 1.200, :vm), 1e-3; atol = 1e-12)
    @test isapprox(Sparlectra._voltage_control_error(1.201, 1.200, :vm2), 2.401e-3; atol = 1e-12)
  end)() end

  @testset "Voltage tap controller (discrete ratio)" begin (function ()
    net, tbr = _build_net()
    result0 = run_sparlectra(net = net, config = _runner_cfg())
    erg0 = result0.numerical_converged ? 0 : 1
    @test erg0 == 0
    vm0 = get_bus_vm_pu(net, "Load")
    @test vm0 > 0.98

    addPowerTransformerControl!(net;
      trafo = string(tbr.branchIdx),
      mode = :voltage,
      target_bus = "Load",
      target_vm_pu = 0.98,
      control_ratio = true,
      control_phase = false,
      is_discrete = true,
      deadband_vm_pu = 5e-3,
      max_outer_iters = 8,
    )

    ratio_before = tbr.tap_ratio
    result = run_sparlectra(net = net, config = _runner_cfg())
    erg = result.numerical_converged ? 0 : 1
    @test erg == 0
    @test tbr.tap_ratio != ratio_before
    @test tbr.tap_ratio > ratio_before
    ctrl = Sparlectra._tap_controllers(net)[1]
    if ctrl.converged
      @test abs(get_bus_vm_pu(net, "Load") - 0.98) <= 0.05
    else
      @test ctrl.at_limit
    end

    # Two AC islands that both have branches. The direction probe solved
    # with islands switched off and refused this net with the
    # islands-disabled error, although the run solves it island by island.
    # The controlled island must end exactly as the same island alone.
    pf = PowerFlowConfig(max_iter = 30, tol = 1e-9)
    add_tap_ctrl! = n -> addPowerTransformerControl!(n; trafo = string(getNetBranch(net = n, fromBus = "Slack", toBus = "Mid").branchIdx), mode = :voltage, target_bus = "Load", target_vm_pu = 1.035, control_ratio = true, control_phase = false, is_discrete = true, deadband_vm_pu = 5e-3, max_outer_iters = 8)
    rnet, rtbr = _build_net()
    add_tap_ctrl!(rnet)
    rres = run_control!(rnet; pf_config = pf)
    inet, itbr = _build_two_island_net()
    islands = Sparlectra.detect_ac_islands(inet).rows
    @test length(islands) == 2 && all(row -> row.n_branch > 0, islands)
    _, ierg0 = redirect_stdout(() -> runpf!(inet; config = pf), devnull)
    @test ierg0 == 0
    @test abs(get_bus_vm_pu(inet, "Load") - 1.035) > 5e-3   # starts outside the deadband
    add_tap_ctrl!(inet)
    ires = redirect_stdout(() -> run_control!(inet; pf_config = pf), devnull)
    @test ires.status === :converged
    @test only(Sparlectra._tap_controllers(inet)).converged
    @test itbr.tap_ratio != 1.0
    @test abs(get_bus_vm_pu(inet, "Load") - 1.035) <= 5e-3
    @test rres.status === :converged
    @test itbr.tap_ratio == rtbr.tap_ratio
    @test isapprox(get_bus_vm_pu(inet, "Load"), get_bus_vm_pu(rnet, "Load"); atol = 1e-9)

    # An isolated bus keeps its type through the solve. Every pass after the
    # first solves the same net again; the bus written back as PQ was an
    # island without a reference there and the second solve stopped.
    xnet, xtbr = _build_net()
    addBus!(net = xnet, busName = "Spare", vn_kV = 110.0)
    addProsumer!(net = xnet, busName = "Spare", type = "ENERGYCONSUMER", p = -5.0, q = -1.0)
    @test markIsolatedBuses!(net = xnet, log = false) == [4]
    add_tap_ctrl!(xnet)
    xres = redirect_stdout(() -> run_control!(xnet; pf_config = pf), devnull)
    @test xres.status === :converged && xres.powerflow_solves > 1
    @test Sparlectra.isIsolated(xnet.nodeVec[4])
    @test xtbr.tap_ratio == rtbr.tap_ratio
    @test isapprox(get_bus_vm_pu(xnet, "Load"), get_bus_vm_pu(rnet, "Load"); atol = 1e-9)
    @test redirect_stdout(() -> runpf!(xnet; config = pf), devnull)[2] == 0

    # A failed probe solve on the same net. The island path throws instead
    # of returning a status, the probe must still answer as on one island
    # (0 when the perturbed solve fails) and restore the branch. One Newton
    # iteration keeps the solved net converged but not the probe step.
    tight = PowerFlowConfig(max_iter = 1, tol = 1e-9)
    r_run = itbr.tap_ratio
    @test r_run + itbr.tap_step <= itbr.tap_max   # the probe can step up
    rdir = redirect_stdout(() -> Sparlectra._ratio_probe_direction(inet, itbr, only(Sparlectra._tap_controllers(inet)), tight), devnull)
    @test rdir == 0.0
    @test itbr.tap_ratio == r_run && itbr.ratio == r_run
    Sparlectra.clearTapControllers!(inet)
    addPowerTransformerControl!(inet; trafo = string(itbr.branchIdx), mode = :branch_active_power, target_branch = ("Slack", "Mid"), p_target_mw = 0.0, control_ratio = false, control_phase = true)
    redirect_stdout(() -> runpf!(inet; config = pf), devnull)
    phi_run, x_run = itbr.phase_shift_deg, itbr.x_pu
    pdir = redirect_stdout(() -> Sparlectra._phase_probe_direction(inet, itbr, only(Sparlectra._tap_controllers(inet)), tight), devnull)
    @test pdir == 0.0
    @test itbr.phase_shift_deg == phi_run && itbr.angle == phi_run && itbr.x_pu == x_run

    # an unreachable target: the tap runs to tap_min and the loop ends
    # :blocked in a pass whose probe solved one step above the limit. The
    # closing solve leaves the voltages of the final tap on the net, equal
    # to a single solve at that tap, and its iterations count in the total.
    blocked_run = function (max_outer::Int)
      bnet, bbr = _build_net()
      addPowerTransformerControl!(bnet; trafo = string(bbr.branchIdx), mode = :voltage, target_bus = "Load", target_vm_pu = 1.2, deadband_vm_pu = 5e-3)
      return bnet, bbr, run_control!(bnet; control_config = ControlConfig(max_outer_iterations = max_outer))
    end
    bnet, bbr, bres = blocked_run(20)
    @test bres.status === :blocked
    @test bbr.tap_ratio == bbr.tap_min && only(Sparlectra._tap_controllers(bnet)).at_limit
    @test bres.powerflow_solves == bres.outer_iterations + 1
    snet, sbr = _build_net()
    (sbr.tap_ratio, sbr.ratio) = (bbr.tap_ratio, bbr.ratio)
    @test runpf!(snet, 30, 1e-10, 0)[2] == 0
    @test maximum(abs(bnet.nodeVec[k]._vm_pu - snet.nodeVec[k]._vm_pu) for k in eachindex(snet.nodeVec)) < 1e-8
    _, _, short = blocked_run(bres.outer_iterations - 1)
    @test short.status === :max_outer_iterations
    @test bres.last_pf_iterations > 0 && bres.total_pf_iterations == short.total_pf_iterations + bres.last_pf_iterations
  end)() end

  @testset "Disabled tap controllers are skipped" begin (function ()
    net, tbr = _build_net()
    addPowerTransformerControl!(net;
      trafo = string(tbr.branchIdx),
      mode = :voltage,
      target_bus = "Load",
      target_vm_pu = 0.98,
      control_ratio = true,
      control_phase = false,
      enabled = false,
    )

    ratio_before = tbr.tap_ratio
    result = run_sparlectra(net = net, config = _runner_cfg())
    erg = result.numerical_converged ? 0 : 1
    @test erg == 0
    @test tbr.tap_ratio == ratio_before
    @test result.control_status === :disabled
    @test result.numerical_converged
    @test result.solution_available
    @test result.final_converged
    @test result.outcome === :converged
    @test result.reason === :none
  end)() end

  @testset "Disabled control framework preserves successful baseline PF" begin (function ()
    net, tbr = _build_net()
    addPowerTransformerControl!(net;
      trafo = string(tbr.branchIdx),
      mode = :voltage,
      target_bus = "Load",
      target_vm_pu = 0.98,
      control_ratio = true,
      control_phase = false,
    )

    ratio_before = tbr.tap_ratio
    result = run_sparlectra(net = net, config = _runner_cfg(control = ControlConfig(enabled = false)))
    @test tbr.tap_ratio == ratio_before
    @test result.control_status === :disabled
    @test result.numerical_converged
    @test result.solution_available
    @test result.final_converged
    @test result.outcome === :converged
    @test result.reason === :none
    @test latest_control_result(net).powerflow_solves == 1
  end)() end

  @testset "Controlled framework result composes converged control status" begin (function ()
    net, tbr = _build_net()
    baseline = run_sparlectra(net = net, config = _runner_cfg())
    @test baseline.final_converged
    addPowerTransformerControl!(net;
      trafo = string(tbr.branchIdx),
      mode = :voltage,
      target_bus = "Load",
      target_vm_pu = get_bus_vm_pu(net, "Load"),
      control_ratio = true,
      control_phase = false,
    )
    result = run_sparlectra(net = net, config = _runner_cfg())
    @test result.control_status === :converged
    @test result.numerical_converged
    @test result.solution_available
    @test result.final_converged
    @test result.outcome === :converged
  end)() end

  @testset "Controlled framework result rejects exhausted outer loop" begin (function ()
    net, tbr = _build_net()
    addPowerTransformerControl!(net;
      trafo = string(tbr.branchIdx),
      mode = :voltage,
      target_bus = "Load",
      target_vm_pu = 0.90,
      control_ratio = true,
      control_phase = false,
      is_discrete = true,
      max_outer_iters = 1,
    )
    # the inner solve has max_iter = 30; the outer loop stops at its own cap of 1
    result = run_sparlectra(net = net, config = _runner_cfg(control = ControlConfig(max_outer_iterations = 1)))
    @test result.control_status === :max_outer_iterations
    @test result.numerical_converged
    @test result.solution_available
    @test !result.final_converged
    @test result.outcome === :control_max_outer_iterations
    @test result.reason === :control_max_outer_iterations
    @test occursin("max_outer_iterations", result.reason_text)
  end)() end

  @testset "Direct run_control! returns ControlRunResult" begin (function ()
    net, tbr = _build_net()
    addPowerTransformerControl!(net;
      trafo = string(tbr.branchIdx),
      mode = :voltage,
      target_bus = "Load",
      target_vm_pu = 0.98,
      control_ratio = true,
      control_phase = false,
      is_discrete = true,
      max_outer_iters = 4,
    )
    # once with an explicit PF config, once with the default one resolved
    pf_cfg = PowerFlowConfig(max_iter = 30, tol = 1e-9, method = :rectangular)
    for pf_kw in ((; pf_config = pf_cfg), (;))
      result = run_control!(net; controllers = collect_outer_controllers(net), control_config = ControlConfig(max_outer_iterations = 4), verbose = 0, pf_kw...)
      @test result isa ControlRunResult
      @test result.powerflow_solves >= 1
      @test result.status in (:converged, :blocked, :max_outer_iterations, :pf_failed)
      @test !isempty(result.controllers)
      @test latest_control_result(net) === result
      @test net.control_result === result
    end
  end)() end

  @testset "run_control! stores terminal no-controller result on Net" begin (function ()
    net, _ = _build_net()
    @test net.control_result === nothing
    @test latest_control_result(net) === nothing
    result = run_control!(net; controllers = AbstractOuterController[], pf_config = PowerFlowConfig(max_iter = 2), control_config = ControlConfig(enabled = true), verbose = 0)
    @test result.status == :no_controllers
    @test latest_control_result(net) === result
    @test net.control_result === result
  end)() end

  @testset "Disabled control still runs one baseline PF" begin (function ()
    # switched off on the controller, or for the whole run in ControlConfig
    for (ctrl_enabled, run_enabled) in ((false, true), (true, false))
      net, tbr = _build_net()
      addPowerTransformerControl!(net;
        trafo = string(tbr.branchIdx),
        mode = :voltage,
        target_bus = "Load",
        target_vm_pu = 0.98,
        control_ratio = true,
        control_phase = false,
        enabled = ctrl_enabled,
      )
      result = run_control!(net;
        controllers = collect_outer_controllers(net),
        pf_config = PowerFlowConfig(method = :rectangular, max_iter = 30, tol = 1e-9),
        control_config = ControlConfig(enabled = run_enabled),
        verbose = 0,
      )
      @test result isa ControlRunResult
      @test result.status == :disabled
      @test result.converged == true
      @test result.powerflow_solves == 1
      @test result.last_pf_iterations >= 1
      @test result.total_pf_iterations == result.last_pf_iterations
      @test result.last_pf_status == :ok
      @test latest_control_result(net) === result
      @test net.control_result === result
    end

    # a controller that switches itself off during initialization leaves no
    # active controller after the first solve; that solve still counts
    snet, _ = _build_net()
    sres = run_control!(snet;
      controllers = AbstractOuterController[_SelfDisablingTestController(true)],
      pf_config = PowerFlowConfig(method = :rectangular, max_iter = 30, tol = 1e-9),
      verbose = 0,
    )
    @test sres.status == :no_active_controllers
    @test sres.powerflow_solves == 1
    @test sres.last_pf_iterations >= 1
    @test sres.total_pf_iterations == sres.last_pf_iterations
  end)() end

  @testset "Tap controller report row reads the branch flow" begin (function ()
    add_ctrl! = (net, tbr) -> addPowerTransformerControl!(net;
      trafo = string(tbr.branchIdx),
      mode = :voltage_and_branch_active_power,
      target_bus = "Load",
      target_vm_pu = 0.98,
      target_branch = ("Slack", "Mid"),
      p_target_mw = 10.0,
      control_ratio = true,
      control_phase = true,
      is_discrete = true,
      max_outer_iters = 4,
    )
    # a bare runpf! leaves the branch flows unset (calcNetLosses! writes
    # them): the report builds and the row says missing, then the value
    # after calcNetLosses!; it threw a FieldError on the unset flow before
    bare, bare_tbr = _build_net()
    add_ctrl!(bare, bare_tbr)
    ite, erg = runpf!(bare, 30, 1e-8, 0)
    @test erg == 0
    @test ismissing(only(buildACPFlowReport(bare; ite = ite).transformer_controls).achieved_p_mw)
    calcNetLosses!(bare)
    @test only(buildACPFlowReport(bare; ite = ite).transformer_controls).achieved_p_mw == get_branch_p_from_to_mw(bare, "Slack", "Mid")

    net, tbr = _build_net()
    add_ctrl!(net, tbr)
    result = run_sparlectra(net = net, config = _runner_cfg())
    erg = result.numerical_converged ? 0 : 1
    @test erg == 0

    report = buildACPFlowReport(net; ct = 0.0, ite = 1, tol = 1e-9, converged = true, solver = :rectangular)
    @test length(report.transformer_controls) == 1
    @test !ismissing(only(report.transformer_controls).achieved_p_mw)
    # the classic summary section prints for the controlled net
    @test !isempty(sprint(io -> printTapControllerSummary(io, net)))
  end)() end

  @testset "Controller can be attached during transformer creation" begin (function ()
    net = Net(name = "tap_ctrl_on_create", baseMVA = 100.0)
    addBus!(net = net, busName = "B1", vn_kV = 110.0)
    addBus!(net = net, busName = "B2", vn_kV = 110.0)
    addBus!(net = net, busName = "B3", vn_kV = 110.0)
    addProsumer!(net = net, busName = "B1", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.01, va_deg = 0.0, referencePri = "B1")
    addProsumer!(net = net, busName = "B3", type = "ENERGYCONSUMER", p = -20.0, q = -5.0)
    addPIModelACLine!(net = net, fromBus = "B2", toBus = "B3", r_pu = 0.02, x_pu = 0.12, b_pu = 0.01, status = 1)

    ctrl = PowerTransformerControl(trafo = "", mode = :voltage, target_bus = "B3", target_vm_pu = 0.99, control_ratio = true, control_phase = false)
    addPIModelTrafo!(net = net, fromBus = "B1", toBus = "B2", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, ratio = 1.0, shift_deg = 0.0, status = 1, controls = [ctrl])
    ctrls = Sparlectra._tap_controllers(net)
    @test length(ctrls) == 1
    @test ctrls[1].trafo == string(getNetBranch(net = net, fromBus = "B1", toBus = "B2").branchIdx)
    @test length(net.trafos[1].side1.controls) == 1
    @test net.trafos[1].tapSideNumber == 1
  end)() end

  # --- machine remote voltage control (#294 point 3) ------------------------
  # Shares this file with the tap controllers because both exercise the same
  # AbstractOuterController framework and run_control! loop.

  function _build_rvc_net(; qmin::Float64 = -50.0, qmax::Float64 = 50.0)
    net = Net(name = "machine_rvc", baseMVA = 100.0)
    addBus!(net = net, busName = "Slack", vn_kV = 110.0)
    addBus!(net = net, busName = "GenBus", vn_kV = 110.0)
    addBus!(net = net, busName = "Load", vn_kV = 110.0)
    addProsumer!(net = net, busName = "Slack", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.02, va_deg = 0.0, referencePri = "Slack")
    addProsumer!(net = net, busName = "GenBus", type = "SYNCHRONOUSMACHINE", p = 30.0, q = 0.0, qMin = qmin, qMax = qmax)
    addProsumer!(net = net, busName = "Load", type = "ENERGYCONSUMER", p = -70.0, q = -20.0)
    addPIModelACLine!(net = net, fromBus = "Slack", toBus = "GenBus", r_pu = 0.02, x_pu = 0.12, b_pu = 0.01, status = 1)
    addPIModelACLine!(net = net, fromBus = "GenBus", toBus = "Load", r_pu = 0.02, x_pu = 0.12, b_pu = 0.01, status = 1)
    return net
  end

  @testset "Machine RVC: API validation" begin (function ()
    net = _build_rvc_net()
    @test_throws ErrorException addMachineVoltageControl!(net; bus = "Load", target_bus = "GenBus", target_vm_pu = 1.0)              # no generator at bus
    @test_throws ErrorException addMachineVoltageControl!(net; bus = "GenBus", target_bus = "GenBus", target_vm_pu = 1.0)            # local target
    @test_throws ErrorException addMachineVoltageControl!(net; bus = "GenBus", target_bus = "Slack", target_vm_pu = 1.0)             # voltage-held target
    addMachineVoltageControl!(net; bus = "GenBus", target_bus = "Load", target_vm_pu = 1.05)
    @test_throws ErrorException addMachineVoltageControl!(net; bus = "GenBus", target_bus = "Load", target_vm_pu = 1.04)             # duplicate machine + target
    @test length(collect_outer_controllers(net)) == 1
    clearMachineControllers!(net)
    @test isempty(collect_outer_controllers(net))

    # cross-type check: a tap controller on the same target bus stays PQ, so
    # only the construction-time warning can flag the double regulation
    xnet = _build_rvc_net()
    addBus!(net = xnet, busName = "TapSide", vn_kV = 20.0)
    ctrl = Sparlectra.PowerTransformerControl(trafo = "", mode = :voltage, target_bus = "Load", target_vm_pu = 1.03)
    add2WTrafo!(net = xnet, fromBus = "Load", toBus = "TapSide", sn_mva = 40.0, vk_percent = 12.0, vkr_percent = 0.4, pfe_kw = 20.0, i0_percent = 0.2, controls = [ctrl])
    @test_logs (:warn, r"tap controller already regulates") addMachineVoltageControl!(xnet; bus = "GenBus", target_bus = "Load", target_vm_pu = 1.05)
  end)() end

  @testset "Machine RVC: secant loop reaches a remote target" begin (function ()
    net = _build_rvc_net()
    addMachineVoltageControl!(net; bus = "GenBus", target_bus = "Load", target_vm_pu = 1.05, deadband_vm_pu = 5e-4)
    pf = PowerFlowConfig(max_iter = 30, tol = 1e-9)
    result = run_control!(net; controllers = collect_outer_controllers(net), pf_config = pf, control_config = ControlConfig(max_outer_iterations = 15), verbose = 0)
    @test result.status == :converged
    @test abs(get_bus_vm_pu(net, "Load") - 1.05) <= 5e-4
    ctrl = Sparlectra._machine_controllers(net)[1]
    @test ctrl.converged
    @test !ctrl.at_limit
    @test ctrl.qmin_mvar <= ctrl.q_mvar <= ctrl.qmax_mvar
    # prosumer and bus-level Q stayed coherent with the controller state
    @test isapprox(net.prosumpsVec[2].qVal, ctrl.q_mvar; atol = 1e-9)
    rows = buildMachineControllerReportRows(net)
    @test length(rows) == 1
    @test rows[1].status == "converged"
    @test rows[1].target_bus == "Load"
    # limit row: an unreachable target parks the machine at its reactive
    # bound, at_limit and not converged
    rvc_lim_net = _build_rvc_net(qmin = -2.0, qmax = 2.0)
    addMachineVoltageControl!(rvc_lim_net; bus = "GenBus", target_bus = "Load", target_vm_pu = 1.0, deadband_vm_pu = 5e-4)
    rvc_lim = run_control!(rvc_lim_net; controllers = collect_outer_controllers(rvc_lim_net), pf_config = pf, control_config = ControlConfig(max_outer_iterations = 15), verbose = 0)
    @test rvc_lim.status in (:converged, :blocked)
    rvc_lim_ctrl = Sparlectra._machine_controllers(rvc_lim_net)[1]
    @test (rvc_lim_ctrl.at_limit, rvc_lim_ctrl.converged) == (true, false)
    @test isapprox(rvc_lim_ctrl.q_mvar, -2.0; atol = 1e-6)
    @test abs(get_bus_vm_pu(rvc_lim_net, "Load") - 1.0) > 5e-4

    # The runtime.parallel switch of the run configuration decides the island
    # fan-out of the control solves, while the active configuration says the
    # opposite: a SparlectraConfig handed to run_control! (control loop and
    # disabled path) and the configuration of run_sparlectra, which hands it
    # on to run_control! (with and without the DC seed, which rebuilds it).
    # The fan-out needs a second thread, and gates and CI run with one, so
    # the runs go to a child process with two threads. Its stderr lands in
    # the same text: a warning or error of the child shows in the assertions.
    script = test_scratch_path(".jl")
    write(script, raw"""
    using Sparlectra
    function two_island_net()
      net = Net(name = "machine_rvc_islands", baseMVA = 100.0)
      addBus!(net = net, busName = "Slack", vn_kV = 110.0)
      addBus!(net = net, busName = "GenBus", vn_kV = 110.0)
      addBus!(net = net, busName = "Load", vn_kV = 110.0)
      addProsumer!(net = net, busName = "Slack", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.02, va_deg = 0.0, referencePri = "Slack")
      addProsumer!(net = net, busName = "GenBus", type = "SYNCHRONOUSMACHINE", p = 30.0, q = 0.0, qMin = -50.0, qMax = 50.0)
      addProsumer!(net = net, busName = "Load", type = "ENERGYCONSUMER", p = -70.0, q = -20.0)
      addPIModelACLine!(net = net, fromBus = "Slack", toBus = "GenBus", r_pu = 0.02, x_pu = 0.12, b_pu = 0.01, status = 1)
      addPIModelACLine!(net = net, fromBus = "GenBus", toBus = "Load", r_pu = 0.02, x_pu = 0.12, b_pu = 0.01, status = 1)
      # a separate two-bus island with its own reference
      addBus!(net = net, busName = "Slack2", vn_kV = 110.0)
      addBus!(net = net, busName = "Load2", vn_kV = 110.0)
      addProsumer!(net = net, busName = "Slack2", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "Slack2")
      addProsumer!(net = net, busName = "Load2", type = "ENERGYCONSUMER", p = -20.0, q = -5.0)
      addPIModelACLine!(net = net, fromBus = "Slack2", toBus = "Load2", r_pu = 0.02, x_pu = 0.12, b_pu = 0.01, status = 1)
      addMachineVoltageControl!(net; bus = "GenBus", target_bus = "Load", target_vm_pu = 1.05, deadband_vm_pu = 5e-4)
      return net
    end
    function run_cfg(parallel::Bool; control::Bool = true, dc_seed::Bool = false)
      pf = PowerFlowConfig(max_iter = 30, tol = 1e-9, islands_parallel_min_buses = 1, start_mode = Sparlectra.StartModeConfig(dc_seed_unconditional = dc_seed))
      return SparlectraConfig(
        powerflow = pf,
        control = ControlConfig(enabled = control, max_outer_iterations = 15),
        output = OutputConfig(logfile_results = :off, startup_latency_hint = false),
        runtime = RuntimeConfig(parallel = ParallelRuntimeConfig(enabled = parallel)),
      )
    end
    println("threads ", Threads.nthreads())
    for passed in (false, true), (entry, control, dc_seed) in ((:run_control, true, false), (:run_control, false, false), (:run_sparlectra, true, false), (:run_sparlectra, true, true))
      net = two_island_net()
      profile = Dict{Symbol,Any}(:enabled => true)
      cfg = run_cfg(passed; control = control, dc_seed = dc_seed)
      with_sparlectra_config(run_cfg(!passed)) do
        redirect_stdout(devnull) do
          if entry === :run_control
            run_control!(net; pf_config = cfg, control_config = cfg.control, verbose = 0, performance_profile = profile)
          else
            run_sparlectra(net = net, config = cfg, performance_profile = profile)
          end
        end
      end
      fanout = haskey(get(profile, :timings, Dict{Symbol,Any}()), :parallel_wall_time)
      println("case ", entry, " control=", control, " dc_seed=", dc_seed, " passed=", passed, " status=", latest_control_result(net).status, " fanout=", fanout)
    end
    """)
    library = normpath(joinpath(@__DIR__, ".."))
    buffer = IOBuffer()
    run(pipeline(Cmd(`$(Base.julia_cmd()) --startup-file=no --threads=2 --project=$(library) $(script)`; ignorestatus = true); stdout = buffer, stderr = buffer))
    out = String(take!(buffer))
    cases = collect(eachmatch(r"case (\w+) control=(\w+) dc_seed=(\w+) passed=(\w+) status=(\w+) fanout=(\w+)", out))
    @test occursin("threads 2", out)
    @test !occursin(r"Warning|Error|ERROR", out)
    @test length(cases) == 8
    length(cases) == 8 || println(out)
    for c in cases
      # the passed configuration decides the fan-out, not the active one; the
      # case fields only name a failing case in the test report
      case = (entry = c[1], control = c[2], dc_seed = c[3])
      @test (; case..., fanout = c[6]) == (; case..., fanout = c[4])
      @test (; case..., status = c[5]) == (; case..., status = c[2] == "true" ? "converged" : "disabled")
    end
  end)() end

  @testset "Controller priority: a bank waits for the tap of its bus" begin (function ()
    # the load sags the bus so far that the tap (rank 0, the default) runs
    # to tap_min = 0.95; the bank on the same bus gets the rank under test
    build = function (bank_priority::Int; tap_priority::Int = 0, tap_min = nothing)
      net = Net(name = "prio_ctrl", baseMVA = 100.0)
      for bus in ("Slack", "Mid", "Load")
        addBus!(net = net, busName = bus, vn_kV = 110.0)
      end
      addProsumer!(net = net, busName = "Slack", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.02, va_deg = 0.0, referencePri = "Slack")
      addProsumer!(net = net, busName = "Load", type = "LOAD", p = 90.0, q = 35.0)
      addPIModelTrafo!(net = net, fromBus = "Slack", toBus = "Mid", r_pu = 0.01, x_pu = 0.08, b_pu = 0.00, ratio = 1.0, shift_deg = 0.0, status = 1)
      addPIModelACLine!(net = net, fromBus = "Mid", toBus = "Load", r_pu = 0.02, x_pu = 0.12, b_pu = 0.01, status = 1)
      tbr = getNetBranch(net = net, fromBus = "Slack", toBus = "Mid")
      tbr.has_ratio_tap = true
      (tbr.tap_min, tbr.tap_max, tbr.tap_step) = (0.95, 1.05, 0.0125)
      addPowerTransformerControl!(net; trafo = string(tbr.branchIdx), mode = :voltage, target_bus = "Load", target_vm_pu = 1.0, deadband_vm_pu = 0.005, priority = tap_priority, tap_min = tap_min)
      add_bank = () -> addShuntVoltageControl!(net; bus = "Load", target_vm_pu = 1.0, bs_min_mvar = 0.0, bs_max_mvar = 30.0, step_mvar = 5.0, deadband_vm_pu = 0.02, priority = bank_priority)
      # the cross-type warning (two controllers on one voltage) fires at
      # equal rank; at different ranks with the same target the bank is the
      # intended cascade behind the tap
      bank = bank_priority == tap_priority ? (@test_logs (:warn, r"already regulates the voltage of bus Load \(target 1.0 pu, priority \d+\)") add_bank()) : (@test_logs add_bank())
      return net, tbr, bank
    end

    # one pass: the tap is outside its band and moves, the bank waits; at
    # equal rank both move in the same pass, as without ranks
    for (prio, bank_after_one_pass) in ((1, 0.0), (0, 5.0))
      net, tbr, bank = build(prio)
      run_control!(net; control_config = ControlConfig(max_outer_iterations = 1))
      @test tbr.tap_ratio ≈ 0.9875
      @test bank.bs_mvar == bank_after_one_pass
    end

    # dead time of a sequential QSTS run: the locked tap (rank 0) stays in
    # its rank, is evaluated and moves nothing. Out of its band (Load at
    # 0.895 pu, target 1.0 +/- 0.005) and not at its limit it holds the bank
    # at rank 1; in its band or at tap_min the bank acts. Before the rule a
    # locked tap dropped out of the ranks and the bank went to 30 MVAr.
    for (case, prep, bank_acts) in (("out of band", (_, _) -> nothing, false),
                                    ("in its band", (t, _) -> (t.deadband_vm_pu = 0.15), true),
                                    ("at its limit", (_, br) -> Sparlectra._set_ratio_tap!(br, br.tap_min), true))
      net, tbr, bank = build(1)
      tap = only(Sparlectra._tap_controllers(net))
      prep(tap, tbr)
      ratio_before = tbr.tap_ratio
      tap.locked = true
      res = run_control!(net)
      @test (case, tbr.tap_ratio == ratio_before, bank.bs_mvar > 0.0) == (case, true, bank_acts)
      @test !any(r -> r.controller_type == "PowerTransformerControl", res.trace)
      bank_acts || @test res.status === :blocked && get_bus_vm_pu(net, "Load") < 0.995
    end

    # full run: the bank acts first in a pass in which the tap stands at
    # tap_min, and then lifts the voltage the tap cannot reach
    net, tbr, bank = build(1)
    res = run_control!(net)
    @test tbr.tap_ratio == tbr.tap_min
    @test only(Sparlectra._tap_controllers(net)).at_limit
    @test bank.bs_mvar > 0.0
    first_bank_pass = minimum(r.outer_iteration for r in res.trace if r.controller_type == "ShuntVoltageControl")
    @test first_bank_pass > 1
    @test only(r for r in res.trace if r.controller_type == "PowerTransformerControl" && r.outer_iteration == first_bank_pass).tap_ratio == tbr.tap_min
    # four whole steps leave the ratio a rounding error above tap_min (at
    # its limit already); the snap onto tap_min in the next pass is a move
    # along the limit and lets the bank act in that pass
    first_tap_at_limit = minimum(r.outer_iteration for r in res.trace if r.controller_type == "PowerTransformerControl" && r.at_limit)
    @test first_bank_pass == first_tap_at_limit + 1
    # the bank ends inside its band; read from the trace, which holds the
    # voltage of the last loop solve
    last_bank = filter(r -> r.controller_type == "ShuntVoltageControl", res.trace)[end]
    @test last_bank.converged && last_bank.achieved_vm_pu >= 0.98

    # equal ranks act together whatever the number: both at 3 give the run
    # of both at 0, where the bank switches before the tap reaches its limit
    finals = map((0, 3)) do r
      n, b, k = build(r; tap_priority = r)
      run_control!(n)
      (b.tap_ratio, k.bs_mvar, get_bus_vm_pu(n, "Load"))
    end
    @test finals[1] == finals[2]
    @test finals[1][1] > 0.95

    # with different targets the two walk against each other at any ranks,
    # so the warning stays for a bank at rank 1 with its own target
    wnet, _, _ = build(0)
    clearShuntControllers!(wnet)
    @test_logs (:warn, r"already regulates the voltage of bus Load") addShuntVoltageControl!(wnet; bus = "Load", target_vm_pu = 1.03, bs_min_mvar = 0.0, bs_max_mvar = 30.0, step_mvar = 5.0, deadband_vm_pu = 0.005, priority = 1)
    # a target on the edge of the tap band agrees: a bank at 1.01 pu (band
    # 0.02) at rank 1 behind the tap at 1.0 +/- 0.01 pu gives no warning,
    # although 1.01 - 1.0 lies a little above 0.01 in floats
    enet, _, _ = build(0)
    clearShuntControllers!(enet)
    only(Sparlectra._tap_controllers(enet)).deadband_vm_pu = 0.01
    @test_logs addShuntVoltageControl!(enet; bus = "Load", target_vm_pu = 1.01, bs_min_mvar = 0.0, bs_max_mvar = 30.0, step_mvar = 5.0, deadband_vm_pu = 0.02, priority = 1)

    # config surface: every type takes the key; a negative rank fails at
    # apply time naming the entry
    cnet, _, _ = build(0)
    clearShuntControllers!(cnet)
    @test all(spec -> "priority" in spec.optional, values(Sparlectra._CONTROLLER_TYPE_SPECS))
    bank_entry = Dict{String,Any}("type" => "shunt_voltage", "name" => "bank", "bus" => "Mid", "target_vm_pu" => 1.0, "bs_min_mvar" => 0.0, "bs_max_mvar" => 30.0, "step_mvar" => 5.0, "priority" => 2)
    @test applyConfiguredControllers!(cnet, ControlConfig(enabled = true, controllers = Any[bank_entry])) == 1
    @test Sparlectra.control_priority(only(Sparlectra._shunt_controllers(cnet))) == 2
    clearShuntControllers!(cnet)
    bank_entry["priority"] = -1
    err = try
      applyConfiguredControllers!(cnet, ControlConfig(enabled = true, controllers = Any[bank_entry]))
      nothing
    catch e
      e
    end
    @test err isa ArgumentError && occursin("control.controllers.bank", err.msg) && occursin("priority must be 0 or larger", err.msg)

    # a control range on the tap inside its physical range 0.95 .. 1.05:
    # the tap stops at 0.975, stands at its limit there, and the bank
    # switches in the pass after; the ratio never leaves the range
    rnet, rtbr, rbank = build(1; tap_min = 0.975)
    rres = run_control!(rnet)
    rtaps = filter(r -> r.controller_type == "PowerTransformerControl", rres.trace)
    @test rtbr.tap_ratio == 0.975 && only(Sparlectra._tap_controllers(rnet)).at_limit
    @test all(r -> r.tap_ratio >= 0.975 - 1e-12, rtaps)
    @test rbank.bs_mvar > 0.0
    @test minimum(r.outer_iteration for r in rres.trace if r.controller_type == "ShuntVoltageControl") == minimum(r.outer_iteration for r in rtaps if r.at_limit) + 1
    @test (only(buildTapControllerReportRows(rnet)).ratio_tap_min, only(buildTapControllerReportRows(rnet)).ratio_tap_max) == (0.975, 1.05)
    # the same range on a tabular changer: the tap steps over its table
    # points, stops on 0.98 (the next point 0.9625 lies outside the range),
    # stands at its limit there, and the bank switches in the pass after
    pnet, ptbr, pbank = build(1; tap_min = 0.975)
    ptbr.tap_ratio_points = [0.95, 0.9625, 0.98, 0.99, 1.0, 1.0125, 1.025, 1.05]
    pres = run_control!(pnet)
    ptaps = filter(r -> r.controller_type == "PowerTransformerControl", pres.trace)
    @test ptbr.tap_ratio == 0.98 && only(Sparlectra._tap_controllers(pnet)).at_limit
    @test all(r -> r.tap_ratio in ptbr.tap_ratio_points && r.tap_ratio >= 0.975, ptaps)
    @test pbank.bs_mvar > 0.0
    @test minimum(r.outer_iteration for r in pres.trace if r.controller_type == "ShuntVoltageControl") == minimum(r.outer_iteration for r in ptaps if r.at_limit) + 1

    # B7: a bank at rank 0 on its top block (as carried over from a step
    # before), target 0.972 pu below the tap band 0.975 .. 0.985 pu. The
    # bus starts below the bank band, so the bank stands at its limit while
    # the tap at rank 1 steps up. The last tap step lifts the bus from
    # 0.965 pu over the bank band (0.967 .. 0.977 pu) to 0.981 pu. The bank
    # is then not at its limit any more (the error points inward) and has
    # to switch back. Before the fix the flag of the earlier pass stayed
    # set and the run ended :converged with the bank on its top block.
    b7net = Net(name = "b7", baseMVA = 100.0)
    for bus in ("Slack", "Mid", "Load")
      addBus!(net = b7net, busName = bus, vn_kV = 110.0)
    end
    addProsumer!(net = b7net, busName = "Slack", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.02, va_deg = 0.0, referencePri = "Slack")
    addProsumer!(net = b7net, busName = "Load", type = "LOAD", p = 90.0, q = 35.0)
    addPIModelTrafo!(net = b7net, fromBus = "Slack", toBus = "Mid", r_pu = 0.01, x_pu = 0.08, b_pu = 0.00, ratio = 1.0, shift_deg = 0.0, status = 1)
    addPIModelACLine!(net = b7net, fromBus = "Mid", toBus = "Load", r_pu = 0.02, x_pu = 0.12, b_pu = 0.01, status = 1)
    b7tbr = getNetBranch(net = b7net, fromBus = "Slack", toBus = "Mid")
    b7tbr.has_ratio_tap = true
    (b7tbr.tap_min, b7tbr.tap_max, b7tbr.tap_step) = (0.90, 1.10, 0.0125)
    addPowerTransformerControl!(b7net; trafo = string(b7tbr.branchIdx), mode = :voltage, target_bus = "Load", target_vm_pu = 0.98, deadband_vm_pu = 0.005, priority = 1)
    b7bank = @test_logs (:warn, r"already regulates the voltage of bus Load") addShuntVoltageControl!(b7net; bus = "Load", target_vm_pu = 0.972, bs_min_mvar = 0.0, bs_max_mvar = 10.0, step_mvar = 5.0, deadband_vm_pu = 0.005, bs_start_mvar = 10.0, priority = 0)
    b7res = run_control!(b7net)
    b7rows = filter(r -> r.controller_type == "ShuntVoltageControl", b7res.trace)
    # the bank switched back from its top block
    @test b7bank.bs_mvar < 10.0
    # the flag stands only on an outer block with the error pointing out:
    # on the top block with the bus below the band, on block 0 above it
    @test any(r -> r.at_limit, b7rows)
    @test all(r -> !r.at_limit || (r.bs_mvar == 10.0 && r.achieved_vm_pu < 0.967) || (r.bs_mvar == 0.0 && r.achieved_vm_pu > 0.977), b7rows)
    # no :converged with the bank outside its band
    @test !(b7res.status === :converged && abs(get_bus_vm_pu(b7net, "Load") - 0.972) > 0.005)
    # the same with parked: bank band 0.95 .. 0.97 pu, tap at 1.0 pu. The
    # bank starts at 0 and climbs to its top block itself, or starts there;
    # the tap then lifts the bus far above the bank band. A flag kept from
    # an earlier pass, or a secant over the bank's own step with the tap's
    # moves mixed in, left the bank parked there and the run :converged.
    # The tap swings between two steps afterwards (its band is narrower
    # than one step), so its status is not checked.
    for start in (0.0, 10.0)
      pknet = Net(name = "park", baseMVA = 100.0)
      for bus in ("Slack", "Mid", "Load")
        addBus!(net = pknet, busName = bus, vn_kV = 110.0)
      end
      addProsumer!(net = pknet, busName = "Slack", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.02, va_deg = 0.0, referencePri = "Slack")
      addProsumer!(net = pknet, busName = "Load", type = "LOAD", p = 90.0, q = 35.0)
      addPIModelTrafo!(net = pknet, fromBus = "Slack", toBus = "Mid", r_pu = 0.01, x_pu = 0.08, b_pu = 0.00, ratio = 1.0, shift_deg = 0.0, status = 1)
      addPIModelACLine!(net = pknet, fromBus = "Mid", toBus = "Load", r_pu = 0.02, x_pu = 0.12, b_pu = 0.01, status = 1)
      pktbr = getNetBranch(net = pknet, fromBus = "Slack", toBus = "Mid")
      pktbr.has_ratio_tap = true
      (pktbr.tap_min, pktbr.tap_max, pktbr.tap_step) = (0.90, 1.10, 0.0125)
      addPowerTransformerControl!(pknet; trafo = string(pktbr.branchIdx), mode = :voltage, target_bus = "Load", target_vm_pu = 1.0, deadband_vm_pu = 0.005, priority = 1)
      pkbank = @test_logs (:warn, r"already regulates the voltage of bus Load") addShuntVoltageControl!(pknet; bus = "Load", target_vm_pu = 0.96, bs_min_mvar = 0.0, bs_max_mvar = 10.0, step_mvar = 5.0, deadband_vm_pu = 0.01, bs_start_mvar = start, priority = 0)
      pkres = run_control!(pknet)
      pkrows = filter(r -> r.controller_type == "ShuntVoltageControl", pkres.trace)
      @test any(r -> r.bs_mvar == 10.0, pkrows)
      @test !(pkres.status === :converged && pkbank.parked && abs(get_bus_vm_pu(pknet, "Load") - 0.96) > 0.01)
    end

    # config surface: an entry for a transformer that already carries a tap
    # controller adds nothing but sets that controller's range, a side the
    # entry leaves out keeps its value; a range outside the physical one or
    # an empty range fails naming the entry, on both paths
    apply_tap = (n, extra) -> applyConfiguredControllers!(n, ControlConfig(enabled = true, controllers = Any[merge(Dict{String,Any}("type" => "power_transformer", "name" => "tap", "trafo" => string(getNetBranch(net = n, fromBus = "Slack", toBus = "Mid").branchIdx), "mode" => "voltage"), extra)]))
    apply_err = (n, extra) -> try
      apply_tap(n, extra)
      nothing
    catch e
      e
    end
    tnet, _, _ = build(1)
    tctrl = only(Sparlectra._tap_controllers(tnet))
    @test apply_tap(tnet, Dict{String,Any}("tap_min" => 0.975)) == 0
    @test (tctrl.tap_min, tctrl.tap_max) == (0.975, nothing)
    @test apply_tap(tnet, Dict{String,Any}("tap_max" => 1.0)) == 0
    @test (tctrl.tap_min, tctrl.tap_max) == (0.975, 1.0)
    for (extra, text) in ((Dict{String,Any}("tap_min" => 0.9), "lies outside the physical tap range 0.95 .. 1.05"), (Dict{String,Any}("tap_min" => 1.0, "tap_max" => 0.99), "must be below tap_max"))
      e = apply_err(tnet, extra)
      @test e isa ArgumentError && occursin("control.controllers.tap", e.msg) && occursin(text, e.msg)
    end
    Sparlectra.clearTapControllers!(tnet)
    e = apply_err(tnet, Dict{String,Any}("target_bus" => "Load", "target_vm_pu" => 1.0, "tap_max" => 1.2))
    @test e isa ArgumentError && occursin("control.controllers.tap", e.msg) && occursin("tap_max = 1.2 lies outside", e.msg)
    @test apply_tap(tnet, Dict{String,Any}("target_bus" => "Load", "target_vm_pu" => 1.0, "tap_min" => 0.975)) == 1
    @test only(Sparlectra._tap_controllers(tnet)).tap_min == 0.975

    # a start ratio of 1.0 outside the control range (tap_min above it)
    # goes onto the nearer end of the range in the first pass, whatever the
    # voltage error: in band, and with the voltage too low, which alone asks
    # for a lower ratio. Plain grid (tap_min 1.0125 is that end), a
    # modelled transformer (1.5 kV steps on 110 kV: tap_min 1.02 lies
    # between step 1 at 1.0136 and step 2 at 1.0273, the first inside) and
    # a tabular changer (tap_min 1.02 between the points 1.0125 and 1.03,
    # the tap goes onto 1.03 and stays on its table). The move logs one
    # warning with the transformer, its start ratio and step, the limit and
    # the ratio it goes to.
    start_net_grid = function (modelled::Bool)
      net = Net(name = "range_start", baseMVA = 100.0)
      for bus in ("Slack", "Mid", "Load")
        addBus!(net = net, busName = bus, vn_kV = 110.0)
      end
      addProsumer!(net = net, busName = "Slack", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.02, va_deg = 0.0, referencePri = "Slack")
      addProsumer!(net = net, busName = "Load", type = "LOAD", p = 70.0, q = 20.0)
      taps = modelled ? PowerTransformerTaps(Vn_kV = 110.0, step = 0, lowStep = -8, highStep = 8, neutralStep = 0, voltageIncrement_kV = 1.5) : nothing
      addPIModelTrafo!(net = net, fromBus = "Slack", toBus = "Mid", r_pu = 0.01, x_pu = 0.08, b_pu = 0.00, ratio = 1.0, shift_deg = 0.0, status = 1, taps = taps)
      addPIModelACLine!(net = net, fromBus = "Mid", toBus = "Load", r_pu = 0.02, x_pu = 0.12, b_pu = 0.01, status = 1)
      sbr = getNetBranch(net = net, fromBus = "Slack", toBus = "Mid")
      modelled || ((sbr.has_ratio_tap, sbr.tap_min, sbr.tap_max, sbr.tap_step) = (true, 0.95, 1.05, 0.0125))
      return net, sbr
    end
    start_net = function (kind::Symbol)
      net, sbr = start_net_grid(kind === :modelled)
      kind === :tabular && (sbr.tap_ratio_points = [0.95, 0.975, 1.0, 1.0125, 1.03, 1.05])
      return net, sbr
    end
    v_start = let n = first(start_net(:plain))
      runpf!(n, 30, 1e-10, 0)
      get_bus_vm_pu(n, "Load")
    end
    for (kind, lo, ratio_in) in ((:plain, 1.0125, 1.0125), (:modelled, 1.02, 1.0 + 3.0 / 110.0), (:tabular, 1.02, 1.03)), dv in (0.0, 0.03)
      snet, sbr = start_net(kind)
      addPowerTransformerControl!(snet; trafo = string(sbr.branchIdx), mode = :voltage, target_bus = "Load", target_vm_pu = v_start + dv, deadband_vm_pu = 0.005, tap_min = lo)
      start_label, end_label = Dict(:plain => ("tap position \\+0", "tap position \\+1"), :modelled => ("model step 0", "model step 2"), :tabular => ("table point 3 of 6", "table point 5 of 6"))[kind]
      start_warning = Regex(string("transformer ", sbr.comp.cName, " starts at ratio 1\\.00000 \\(", start_label,
        "\\) below its control range \\(tap_min = ", replace(string(lo), "." => "\\."), "\\) and goes to ", replace(@sprintf("%.5f", ratio_in), "." => "\\."),
        " \\(", end_label, "\\) in the first pass"))
      sres = @test_logs (:warn, start_warning) run_control!(snet)
      @test sbr.tap_ratio ≈ ratio_in atol = 1e-12
      @test first(r.tap_ratio for r in sres.trace if r.controller_type == "PowerTransformerControl") ≈ ratio_in atol = 1e-12
      @test sres.status === :blocked && only(Sparlectra._tap_controllers(snet)).at_limit
    end
  end)() end

end
