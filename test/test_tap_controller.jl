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
    pf_cfg = PowerFlowConfig(max_iter = 30, tol = 1e-9, method = :rectangular)
    result = run_control!(net; controllers = collect_outer_controllers(net), pf_config = pf_cfg, control_config = ControlConfig(max_outer_iterations = 4), verbose = 0)
    @test result isa ControlRunResult
    @test result.powerflow_solves >= 1
    @test result.status in (:converged, :blocked, :max_outer_iterations, :pf_failed)
    @test !isempty(result.controllers)
    @test latest_control_result(net) === result
    @test net.control_result === result
  end)() end

  @testset "Direct run_control! resolves default PF config" begin (function ()
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

    ratio_before = tbr.tap_ratio
    result = run_control!(net; controllers = collect_outer_controllers(net), control_config = ControlConfig(max_outer_iterations = 4), verbose = 0)
    @test result isa ControlRunResult
    @test result.powerflow_solves >= 1
    @test result.status in (:converged, :blocked, :max_outer_iterations, :pf_failed)
    @test latest_control_result(net) === result
    @test net.control_result === result
    @test tbr.tap_ratio != ratio_before || result.status in (:converged, :blocked, :max_outer_iterations, :pf_failed)
  end)() end

  @testset "Net initializes control_result as nothing" begin (function ()
    net = Net(name = "x", baseMVA = 100.0)
    @test net.control_result === nothing
    @test latest_control_result(net) === nothing
  end)() end

  @testset "run_control! stores terminal no-controller result on Net" begin (function ()
    net, _ = _build_net()
    result = run_control!(net; controllers = AbstractOuterController[], pf_config = PowerFlowConfig(max_iter = 2), control_config = ControlConfig(enabled = true), verbose = 0)
    @test result.status == :no_controllers
    @test latest_control_result(net) === result
    @test net.control_result === result
  end)() end

  @testset "Disabled control still runs one baseline PF" begin (function ()
    net, tbr = _build_net()
    addPowerTransformerControl!(net;
      trafo = string(tbr.branchIdx),
      mode = :voltage,
      target_bus = "Load",
      target_vm_pu = 0.98,
      control_ratio = true,
      control_phase = false,
    )
    result = run_control!(net;
      controllers = collect_outer_controllers(net),
      pf_config = PowerFlowConfig(method = :rectangular, max_iter = 30, tol = 1e-9),
      control_config = ControlConfig(enabled = false),
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
  end)() end

  @testset "All disabled controllers still run one baseline PF" begin (function ()
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
    result = run_control!(net;
      controllers = collect_outer_controllers(net),
      pf_config = PowerFlowConfig(method = :rectangular, max_iter = 30, tol = 1e-9),
      control_config = ControlConfig(enabled = true),
      verbose = 0,
    )
    @test result.status == :disabled
    @test result.converged == true
    @test result.powerflow_solves == 1
    @test result.last_pf_status == :ok
    @test latest_control_result(net) === result

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

  @testset "run_sparlectra net path honors per-call cfg.control" begin (function ()
    net, tbr = _build_net()
    addPowerTransformerControl!(net;
      trafo = string(tbr.branchIdx),
      mode = :voltage,
      target_bus = "Load",
      target_vm_pu = 0.98,
      control_ratio = true,
      control_phase = false,
    )
    cfg = SparlectraConfig(
      powerflow = PowerFlowConfig(method = :rectangular, max_iter = 30, tol = 1e-9),
      control = ControlConfig(enabled = false, max_outer_iterations = 3, trace = false),
    )
    result = run_sparlectra(net = net, config = cfg)
    @test result.numerical_converged
    result = latest_control_result(net)
    @test result !== nothing
    @test result.status == :disabled
    @test result.powerflow_solves == 1
    @test result.converged == true
  end)() end

  @testset "Outer-loop limits are separated from inner PF max_iter" begin (function ()
    net, tbr = _build_net()
    addPowerTransformerControl!(net;
      trafo = string(tbr.branchIdx),
      mode = :voltage,
      target_bus = "Load",
      target_vm_pu = 0.90,
      control_ratio = true,
      control_phase = false,
      is_discrete = true,
      max_outer_iters = 5,
    )
    low_inner_pf = PowerFlowConfig(max_iter = 1, tol = 1e-9, method = :rectangular)
    cfg = ControlConfig(max_outer_iterations = 4, trace = true)
    result = run_control!(net; controllers = collect_outer_controllers(net), pf_config = low_inner_pf, control_config = cfg, verbose = 0)
    @test result.outer_iterations <= 4
    @test result.outer_iterations <= 5
    @test result.powerflow_solves >= 1
    @test latest_control_result(net) === result
  end)() end

  @testset "Transformer controller trace rows include stable keys" begin (function ()
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
    pf_cfg = PowerFlowConfig(max_iter = 30, tol = 1e-9, method = :rectangular)
    result = run_control!(net; controllers = collect_outer_controllers(net), pf_config = pf_cfg, control_config = ControlConfig(max_outer_iterations = 4, trace = true), verbose = 0)
    @test !isempty(result.controllers)
    @test !isempty(result.trace)
    row = first(result.trace)
    @test haskey(row, :outer_iteration)
    @test haskey(row, :controller_name)
    @test haskey(row, :controller_type)
    @test haskey(row, :status)
    @test haskey(row, :tap_ratio)
    @test haskey(row, :phase_shift_deg)
    @test row.controller_type == "PowerTransformerControl"
  end)() end


  @testset "Tap controller reporting rows and classic section" begin (function ()
    net, tbr = _build_net()
    addPowerTransformerControl!(net;
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
    result = run_sparlectra(net = net, config = _runner_cfg())
    erg = result.numerical_converged ? 0 : 1
    @test erg == 0

    report = buildACPFlowReport(net; ct = 0.0, ite = 1, tol = 1e-9, converged = true, solver = :rectangular)
    @test length(report.transformer_controls) == 1
    row = report.transformer_controls[1]
    @test haskey(row, :controller_name)
    @test haskey(row, :achieved_p_mw)
    @test !ismissing(row.achieved_p_mw)
    @test row.control_type == "OLTC+PST"

    txt = sprint(io -> printTapControllerSummary(io, net))
    @test occursin("tap position", txt)
    @test occursin("status", txt)
    @test occursin("Power sign convention", txt)
  end)() end

  @testset "Branch derives tap limits from PowerTransformerTaps" begin (function ()
    taps = PowerTransformerTaps(Vn_kV = 110.0, step = 0, lowStep = -4, highStep = 6, neutralStep = 1, voltageIncrement_kV = 1.1)
    w1 = PowerTransformerWinding(110.0, 0.0, 0.12, 0.0, 0.0, 1.0, 0.0, 110.0, 100.0, taps, true, nothing)
    w2 = PowerTransformerWinding(20.0, 0.0, 0.0, 0.0, 0.0, nothing, 0.0, 20.0, 100.0, nothing, true, nothing)
    trafo_comp = Sparlectra.getBranchComp(110.0, 1, 2, 1, "2WT")
    trafo = PowerTransformer(trafo_comp, true, w1, w2, nothing, Sparlectra.Ratio)

    br = Branch(
      branchIdx = 1,
      from = 1,
      to = 2,
      baseMVA = 100.0,
      branch = trafo,
      id = 1,
      ratio = 1.0,
      side = 1,
      vn_kV = 110.0,
      values_are_pu = true,
    )

    pu_per_step = taps.tapStepPercent / 100.0
    expected_min = min(1.0 + (taps.lowStep - taps.neutralStep) * pu_per_step, 1.0 + (taps.highStep - taps.neutralStep) * pu_per_step)
    expected_max = max(1.0 + (taps.lowStep - taps.neutralStep) * pu_per_step, 1.0 + (taps.highStep - taps.neutralStep) * pu_per_step)
    expected_step = abs(pu_per_step)

    @test isapprox(br.tap_min, expected_min; atol = 1e-12)
    @test isapprox(br.tap_max, expected_max; atol = 1e-12)
    @test isapprox(br.tap_step, expected_step; atol = 1e-12)
    @test isapprox(taps.tapStepPercent, 1.0; atol = 1e-12)
    @test isapprox(taps.neutralU, 110.0; atol = 1e-12)
    @test isapprox(taps.neutralU_ratio, 1.0; atol = 1e-12)
  end)() end

  @testset "PowerTransformerTaps neutral voltage derivation" begin (function ()
    taps_ratio = PowerTransformerTaps(Vn_kV = 20.0, step = 0, lowStep = -2, highStep = 2, neutralStep = 0, voltageIncrement_kV = 0.4, neutralU_ratio = 1.05)
    @test isapprox(taps_ratio.neutralU, 21.0; atol = 1e-12)
    @test isapprox(taps_ratio.neutralU_ratio, 1.05; atol = 1e-12)

    taps_neutral = PowerTransformerTaps(Vn_kV = 20.0, step = 0, lowStep = -2, highStep = 2, neutralStep = 0, voltageIncrement_kV = 0.4, neutralU = 19.5)
    @test isapprox(taps_neutral.neutralU_ratio, 0.975; atol = 1e-12)
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

end
