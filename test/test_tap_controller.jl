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

  @testset "Tap controller API validation" begin
    net, _ = _build_net()
    @test_throws ErrorException addPowerTransformerControl!(net; trafo = "does_not_exist", mode = :voltage, target_bus = "Load", target_vm_pu = 1.0)
    @test_throws ErrorException addPowerTransformerControl!(net; trafo = "1", mode = :voltage)
    @test_throws ErrorException addPowerTransformerControl!(net; trafo = "1", mode = :voltage, target_bus = "Load", target_vm_pu = 1.0, voltage_error_metric = :invalid)
  end

  @testset "One-controller guard works across transformer aliases" begin
    net, tbr = _build_net()
    addPowerTransformerControl!(net; trafo = string(tbr.branchIdx), mode = :voltage, target_bus = "Load", target_vm_pu = 0.99, control_ratio = true, control_phase = false)
    @test_throws ErrorException addPowerTransformerControl!(net; trafo = tbr.comp.cID, mode = :voltage, target_bus = "Load", target_vm_pu = 0.99, control_ratio = true, control_phase = false)
  end

  @testset "Split Schraegregelung: two disjoint controllers on one transformer" begin
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
  end

  @testset "PST reactance coupling X(alpha) (#274)" begin
    # PST loop net: trafo Slack→Mid with a parallel line so the phase tap can
    # shift flow; the typed model is attached to the controlled winding
    # closure-local net deliberately NOT named `net`: assigning a name that
    # is also a testset local would write the captured outer variable and
    # alias every built fixture
    build_xalpha = function (; model = nothing)
      pnet = Net(name = "pst_xalpha", baseMVA = 100.0)
      for bus in ("Slack", "Mid", "Load")
        addBus!(net = pnet, busName = bus, vn_kV = 110.0)
      end
      addProsumer!(net = pnet, busName = "Slack", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.02, va_deg = 0.0, referencePri = "Slack")
      addProsumer!(net = pnet, busName = "Load", type = "ENERGYCONSUMER", p = -70.0, q = -20.0)
      addPIModelTrafo!(net = pnet, fromBus = "Slack", toBus = "Mid", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, ratio = 1.0, shift_deg = 0.0, status = 1)
      addPIModelACLine!(net = pnet, fromBus = "Slack", toBus = "Mid", r_pu = 0.03, x_pu = 0.2, b_pu = 0.0, status = 1)
      addPIModelACLine!(net = pnet, fromBus = "Mid", toBus = "Load", r_pu = 0.02, x_pu = 0.12, b_pu = 0.01, status = 1)
      ptbr = getNetBranch(net = pnet, fromBus = "Slack", toBus = "Mid")
      ptbr.has_phase_tap = true
      ptbr.phase_min_deg = -10.0
      ptbr.phase_max_deg = 10.0
      ptbr.phase_step_deg = 0.5
      model === nothing || (pnet.trafos[1].side1.phase_taps = model)
      return pnet, ptbr
    end
    formula_model = () -> Sparlectra.PhaseTapChangerModel(kind = :symmetrical, step = 0, lowStep = -10, highStep = 10, neutralStep = 0, voltage_step_increment = 0.01, x_min = 0.08, x_max = 0.16)

    # mapping unit checks through the controller-resolution path
    net, tbr = build_xalpha(model = formula_model())
    addPowerTransformerControl!(net; trafo = string(tbr.branchIdx), mode = :branch_active_power, target_branch = ("Slack", "Mid"), p_target_mw = 0.0, control_ratio = false, control_phase = true, deadband_p_mw = 2.0)
    ctrl = only(Sparlectra._tap_controllers(net))
    @test Sparlectra._phase_tap_model_of(net, ctrl) === net.trafos[1].side1.phase_taps
    @test Sparlectra._phase_tap_reactance_at(net, ctrl, 0.0) ≈ 0.08 atol = 1e-12
    x5 = Sparlectra._phase_tap_reactance_at(net, ctrl, 5.0)
    @test x5 !== nothing && 0.08 < x5 <= 0.16
    @test Sparlectra._phase_tap_reactance_at(net, ctrl, 5.0) ≈ Sparlectra.calcPhaseTapReactance(formula_model(), 5.0) atol = 1e-12

    # tabular: nearest row by angle, missing x_pu row → nothing
    tab = Sparlectra.PhaseTapChangerModel(
      kind = :tabular,
      step = 1,
      lowStep = 0,
      highStep = 2,
      neutralStep = 1,
      table = [Sparlectra.TapTablePoint(step = 0, ratio = 1.0, angle_deg = -3.0, x_pu = 0.07), Sparlectra.TapTablePoint(step = 1, ratio = 1.0, angle_deg = 0.0, x_pu = 0.08), Sparlectra.TapTablePoint(step = 2, ratio = 1.0, angle_deg = 3.0, x_pu = 0.11)],
    )
    net.trafos[1].side1.phase_taps = tab
    @test Sparlectra._phase_tap_reactance_at(net, ctrl, 2.4) ≈ 0.11 atol = 1e-12   # nearest row: +3°
    @test Sparlectra._phase_tap_reactance_at(net, ctrl, -1.6) ≈ 0.07 atol = 1e-12  # nearest row: -3°
    net.trafos[1].side1.phase_taps = nothing
    @test Sparlectra._phase_tap_reactance_at(net, ctrl, 5.0) === nothing

    # functional: an accepted tap move updates x_pu to the model value at the
    # final angle, and the next solve reflects it
    net1, tbr1 = build_xalpha(model = formula_model())
    r0 = run_sparlectra(net = net1, config = _runner_cfg())
    @test r0.numerical_converged
    p0 = get_branch_p_from_to_mw(net1, "Slack", "Mid")
    addPowerTransformerControl!(net1; trafo = string(tbr1.branchIdx), mode = :branch_active_power, target_branch = ("Slack", "Mid"), p_target_mw = round(p0) + 8.0, control_ratio = false, control_phase = true, deadband_p_mw = 2.0)
    r1 = run_sparlectra(net = net1, config = _runner_cfg())
    @test r1.numerical_converged
    @test abs(tbr1.phase_shift_deg) > 0.0
    mctrl1 = only(Sparlectra._tap_controllers(net1))
    xexp = Sparlectra._phase_tap_reactance_at(net1, mctrl1, tbr1.phase_shift_deg)
    @test xexp !== nothing
    @test tbr1.x_pu ≈ xexp atol = 1e-12
    @test tbr1.x_pu > 0.08

    # without a typed model the run is unchanged: static reactance
    net2, tbr2 = build_xalpha()
    addPowerTransformerControl!(net2; trafo = string(tbr2.branchIdx), mode = :branch_active_power, target_branch = ("Slack", "Mid"), p_target_mw = round(p0) + 8.0, control_ratio = false, control_phase = true, deadband_p_mw = 2.0)
    r2 = run_sparlectra(net = net2, config = _runner_cfg())
    @test r2.numerical_converged
    @test tbr2.x_pu == 0.08

    # probe consistency: tracking active vs. static picks the same direction,
    # and the probe restores angle and reactance
    net3, tbr3 = build_xalpha(model = formula_model())
    run_sparlectra(net = net3, config = _runner_cfg())
    addPowerTransformerControl!(net3; trafo = string(tbr3.branchIdx), mode = :branch_active_power, target_branch = ("Slack", "Mid"), p_target_mw = 0.0, control_ratio = false, control_phase = true)
    ctrl3 = only(Sparlectra._tap_controllers(net3))
    x_before = tbr3.x_pu
    phi_before = tbr3.phase_shift_deg
    dir_model = Sparlectra._phase_probe_direction(net3, tbr3, ctrl3, 30, 1e-8, 0, :rectangular)
    @test tbr3.x_pu == x_before
    @test tbr3.phase_shift_deg == phi_before
    net4, tbr4 = build_xalpha()
    run_sparlectra(net = net4, config = _runner_cfg())
    addPowerTransformerControl!(net4; trafo = string(tbr4.branchIdx), mode = :branch_active_power, target_branch = ("Slack", "Mid"), p_target_mw = 0.0, control_ratio = false, control_phase = true)
    ctrl4 = only(Sparlectra._tap_controllers(net4))
    dir_static = Sparlectra._phase_probe_direction(net4, tbr4, ctrl4, 30, 1e-8, 0, :rectangular)
    @test dir_model == dir_static
    @test abs(dir_model) == 1.0
  end

  @testset "SVC shunt voltage controller + controllable elements (#227)" begin
    # NOTE: the closure-local net must NOT be named `net` — an assignment to
    # a name that is also a testset local would write the captured outer
    # variable (Julia closure capture), aliasing every built fixture.
    build_svc_net = function ()
      svcnet = Net(name = "svc_ctrl", baseMVA = 100.0)
      for bus in ("Slack", "Mid", "Load")
        addBus!(net = svcnet, busName = bus, vn_kV = 110.0)
      end
      addProsumer!(net = svcnet, busName = "Slack", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.02, va_deg = 0.0, referencePri = "Slack")
      addProsumer!(net = svcnet, busName = "Load", type = "ENERGYCONSUMER", p = -80.0, q = -30.0)
      addPIModelACLine!(net = svcnet, fromBus = "Slack", toBus = "Mid", r_pu = 0.02, x_pu = 0.12, b_pu = 0.01, status = 1)
      addPIModelACLine!(net = svcnet, fromBus = "Mid", toBus = "Load", r_pu = 0.02, x_pu = 0.12, b_pu = 0.01, status = 1)
      return svcnet
    end

    # API validation
    vnet = build_svc_net()
    @test_throws ErrorException addShuntVoltageControl!(vnet; bus = "Slack", target_vm_pu = 1.0, bs_min_mvar = -60.0, bs_max_mvar = 60.0)  # voltage-held bus
    @test_throws ErrorException addShuntVoltageControl!(vnet; bus = "Load", target_vm_pu = 1.0, bs_min_mvar = 60.0, bs_max_mvar = -60.0)   # inverted range
    addShuntVoltageControl!(vnet; bus = "Load", target_vm_pu = 1.0, bs_min_mvar = -60.0, bs_max_mvar = 60.0)
    @test_throws ErrorException addShuntVoltageControl!(vnet; bus = "Load", target_vm_pu = 1.01, bs_min_mvar = -60.0, bs_max_mvar = 60.0)  # duplicate bus
    @test length(Sparlectra._shunt_controllers(vnet)) == 1
    clearShuntControllers!(vnet)
    @test isempty(Sparlectra._shunt_controllers(vnet))

    # in-range regulation: the secant loop pulls the overvoltage down to the
    # setpoint by moving into the inductive range
    net = build_svc_net()
    baseline = run_sparlectra(net = net)
    @test baseline.numerical_converged
    vm0 = get_bus_vm_pu(net, "Load")
    @test vm0 > 1.05
    addShuntVoltageControl!(net; bus = "Load", target_vm_pu = 1.0, bs_min_mvar = -60.0, bs_max_mvar = 60.0, deadband_vm_pu = 1e-3)
    result = run_sparlectra(net = net)
    @test result.numerical_converged
    ctrl = only(Sparlectra._shunt_controllers(net))
    @test ctrl.converged
    @test !ctrl.at_limit
    @test abs(get_bus_vm_pu(net, "Load") - 1.0) <= 1e-3
    @test -60.0 < ctrl.bs_mvar < 0.0
    # the actuated shunt carries the same susceptance (MATPOWER stamping)
    sh = net.shuntVec[ctrl.shunt_idx]
    @test imag(sh.y_pu_shunt) * net.baseMVA ≈ ctrl.bs_mvar atol = 1e-9

    # honest limit region: an unreachable target clamps the susceptance and
    # reports at_limit — constant-B, the reactive output then follows V²
    lnet = build_svc_net()
    run_sparlectra(net = lnet)
    addShuntVoltageControl!(lnet; bus = "Load", target_vm_pu = 1.0, bs_min_mvar = -10.0, bs_max_mvar = 10.0, deadband_vm_pu = 1e-3)
    run_sparlectra(net = lnet)
    lctrl = only(Sparlectra._shunt_controllers(lnet))
    @test !lctrl.converged
    @test lctrl.at_limit
    @test lctrl.status == :at_limit
    @test lctrl.bs_mvar == -10.0
    @test get_bus_vm_pu(lnet, "Load") > 1.0 + 1e-3

    # generic controllable-element view: uniform records for all three
    # controller families
    els = controllableElements(net)
    @test length(els) == 1
    @test els[1].device == "SVC (variable shunt)"
    @test els[1].actuator == :shunt_bs_mvar
    @test (els[1].actuator_min, els[1].actuator_max) == (-60.0, 60.0)
    @test els[1].quantity == :bus_voltage
    @test els[1].converged === true
    cres = latest_control_result(net)
    @test cres !== nothing
    @test length(cres.elements) == 1

    mixed = Net(name = "mixed_els", baseMVA = 100.0)
    for bus in ("Slack", "GenBus", "Load")
      addBus!(net = mixed, busName = bus, vn_kV = 110.0)
    end
    addProsumer!(net = mixed, busName = "Slack", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.02, va_deg = 0.0, referencePri = "Slack")
    addProsumer!(net = mixed, busName = "GenBus", type = "SYNCHRONOUSMACHINE", p = 30.0, q = 0.0, qMin = -50.0, qMax = 50.0)
    addProsumer!(net = mixed, busName = "Load", type = "ENERGYCONSUMER", p = -70.0, q = -20.0)
    addPIModelACLine!(net = mixed, fromBus = "Slack", toBus = "GenBus", r_pu = 0.02, x_pu = 0.12, b_pu = 0.01, status = 1)
    addPIModelACLine!(net = mixed, fromBus = "GenBus", toBus = "Load", r_pu = 0.02, x_pu = 0.12, b_pu = 0.01, status = 1)
    addMachineVoltageControl!(mixed; bus = "GenBus", target_bus = "Load", target_vm_pu = 1.03)
    tnet, ttbr = _build_net()
    addPowerTransformerControl!(tnet; trafo = string(ttbr.branchIdx), mode = :voltage, target_bus = "Load", target_vm_pu = 0.99, control_ratio = true, control_phase = false)
    mrows = controllableElements(mixed)
    trows = controllableElements(tnet)
    @test length(mrows) == 1
    @test mrows[1].actuator == :machine_q_mvar
    @test mrows[1].target == "Load"
    @test length(trows) == 1
    @test trows[1].actuator == :tap_ratio
    @test trows[1].quantity == :bus_voltage
    @test trows[1].discrete === true
  end

  @testset "STATCOM current-limit mode (#297 Draft A)" begin
    # weak two-line corridor: the load bus sags visibly, the compensator at
    # Mid works against the sag. NOTE: closure-local nets must not shadow a
    # testset-local `net` (Julia closure capture, see the SVC testset).
    build_weak_net = function (p_load = 60.0, q_load = 25.0)
      wnet = Net(name = "statcom_ctrl", baseMVA = 100.0)
      for bus in ("Slack", "Mid", "Load")
        addBus!(net = wnet, busName = bus, vn_kV = 110.0)
      end
      addProsumer!(net = wnet, busName = "Slack", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "Slack")
      # NOTE the sign convention: type LOAD with POSITIVE p/q consumes (the
      # SVC testset's ENERGYCONSUMER uses negative values for the same thing).
      # Loads are PROSUMER objects: the AC solver builds its injections from
      # prosumpsVec (buildComplexSVec), so load variants are parameterized
      # here instead of edited via the node-sum helpers.
      addProsumer!(net = wnet, busName = "Load", type = "LOAD", p = p_load, q = q_load)
      addProsumer!(net = wnet, busName = "Mid", type = "GENERATOR", p = 0.0, q = 0.0)
      addPIModelACLine!(net = wnet, fromBus = "Slack", toBus = "Mid", r_pu = 0.02, x_pu = 0.20, b_pu = 0.0, status = 1)
      addPIModelACLine!(net = wnet, fromBus = "Mid", toBus = "Load", r_pu = 0.02, x_pu = 0.20, b_pu = 0.0, status = 1)
      ok, msg = validate!(net = wnet)
      ok || error("test net invalid: $msg")
      return wnet
    end

    # registration validation
    vnet = build_weak_net()
    @test_throws ErrorException addMachineVoltageControl!(vnet; bus = "Mid", target_bus = "Load", target_vm_pu = 1.0, s_max_mva = 12.0, i_max_ka = 0.1)      # both ratings
    @test_throws ErrorException addMachineVoltageControl!(vnet; bus = "Mid", target_bus = "Load", target_vm_pu = 1.0, s_max_mva = 12.0, qmin_mvar = -5.0)    # rating + box
    @test_throws ErrorException addMachineVoltageControl!(vnet; bus = "Mid", target_bus = "Load", target_vm_pu = 1.0, s_max_mva = -3.0)                       # negative rating
    @test_throws ErrorException addMachineVoltageControl!(vnet; bus = "Mid", target_bus = "Load", target_vm_pu = 1.0, i_max_ka = -0.1)                        # negative current

    # i_max_ka converts via sqrt(3) * Un * Imax at registration
    cnet = build_weak_net()
    addMachineVoltageControl!(cnet; bus = "Mid", target_bus = "Load", target_vm_pu = 1.0, i_max_ka = 0.05)
    cctrl = only([c for c in cnet.machineControls if c isa MachineVoltageControl])
    @test cctrl.limit_mode === :current
    @test cctrl.s_max_mva ≈ sqrt(3.0) * 110.0 * 0.05 atol = 1e-12

    # at the limit the delivered Q tracks V * S_max: the defining STATCOM
    # acceptance (issue #297 Draft A) — a small rating cannot lift the sag,
    # the controller parks at_limit and q sits ON the live bound
    net = build_weak_net()
    addMachineVoltageControl!(net; bus = "Mid", target_bus = "Load", target_vm_pu = 1.0, s_max_mva = 12.0)
    run_control!(net)
    ctrl = only([c for c in net.machineControls if c isa MachineVoltageControl])
    vt = get_bus_vm_pu(net, "Mid")
    @test !ctrl.converged
    @test ctrl.at_limit
    @test ctrl.status == :at_limit
    @test ctrl.q_mvar ≈ vt * 12.0 atol = 2e-2
    @test ctrl.qmax_mvar ≈ vt * 12.0 atol = 1e-9      # live bound from the solved terminal voltage
    @test ctrl.qmin_mvar ≈ -vt * 12.0 atol = 1e-9     # symmetric converter limit
    @test vt < 1.0                                     # depressed voltage: bound BELOW the 1.0-pu rating
    # element view: STATCOM vocabulary with live bounds
    els = controllableElements(net)
    @test length(els) == 1
    @test els[1].device == "STATCOM (VSC)"
    @test els[1].actuator == :machine_q_mvar
    @test els[1].actuator_max ≈ vt * 12.0 atol = 1e-9
    @test els[1].at_limit === true
    # report row carries the mode fields
    row = only(buildMachineControllerReportRows(net))
    @test row.limit_mode === :current
    @test row.s_max_mva ≈ 12.0
    @test !ismissing(row.machine_vm_pu)

    # deeper sag moves the bound: heavier load, fresh run — the delivered Q
    # must FOLLOW the lower terminal voltage (linear V * I_max tracking)
    hnet = build_weak_net(75.0, 30.0)
    addMachineVoltageControl!(hnet; bus = "Mid", target_bus = "Load", target_vm_pu = 1.0, s_max_mva = 12.0)
    run_control!(hnet)
    hctrl = only([c for c in hnet.machineControls if c isa MachineVoltageControl])
    hvt = get_bus_vm_pu(hnet, "Mid")
    @test hctrl.at_limit
    @test hvt < vt
    @test hctrl.q_mvar ≈ hvt * 12.0 atol = 2e-2
    @test hctrl.q_mvar < ctrl.q_mvar                   # less Q available at the deeper sag

    # in range the STATCOM behaves like the constant-Q controller: with a
    # rating that never binds, both modes converge into the same deadband
    # around the target. (The trajectories differ slightly — the bootstrap
    # step is a fraction of the AVAILABLE headroom, which the live bound
    # changes — so the comparison is deadband-level, not bitwise.)
    anet = build_weak_net()
    addMachineVoltageControl!(anet; bus = "Mid", target_bus = "Load", target_vm_pu = 0.95, s_max_mva = 200.0)
    run_control!(anet)
    bnet = build_weak_net()
    addMachineVoltageControl!(bnet; bus = "Mid", target_bus = "Load", target_vm_pu = 0.95, qmin_mvar = -200.0, qmax_mvar = 200.0)
    run_control!(bnet)
    actrl = only([c for c in anet.machineControls if c isa MachineVoltageControl])
    bctrl = only([c for c in bnet.machineControls if c isa MachineVoltageControl])
    @test actrl.converged
    @test bctrl.converged
    @test abs(get_bus_vm_pu(anet, "Load") - 0.95) <= actrl.deadband_vm_pu
    @test abs(get_bus_vm_pu(bnet, "Load") - 0.95) <= bctrl.deadband_vm_pu
    @test actrl.q_mvar ≈ bctrl.q_mvar atol = 1.0

    # constant-Q mode is untouched by the STATCOM addition: mode fields are
    # inert defaults and the descriptor keeps the RVC vocabulary
    @test bctrl.limit_mode === :constant_q
    @test bctrl.s_max_mva === nothing
    bels = controllableElements(bnet)
    @test bels[1].device == "machine remote voltage control"
  end

  @testset "SVC vs STATCOM limit characteristic (#297 Draft E)" begin
    # The Draft E acceptance: on the SAME depressed-voltage case the SVC's
    # delivered Q collapses with V^2 (constant clamped susceptance through
    # the Y-bus) while the STATCOM's falls only linearly with V (constant
    # current). Both devices are rated identically at 1.0 pu, so at V < 1
    # the STATCOM must deliver MORE reactive support by the factor 1/V.
    rating_mvar = 10.0
    build_sag_net = function ()
      snet = Net(name = "limit_contrast", baseMVA = 100.0)
      for bus in ("Slack", "Mid", "Load")
        addBus!(net = snet, busName = bus, vn_kV = 110.0)
      end
      addProsumer!(net = snet, busName = "Slack", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "Slack")
      # positive LOAD values consume: the corridor sags below 1.0 pu
      addProsumer!(net = snet, busName = "Load", type = "LOAD", p = 60.0, q = 25.0)
      addPIModelACLine!(net = snet, fromBus = "Slack", toBus = "Mid", r_pu = 0.02, x_pu = 0.20, b_pu = 0.0, status = 1)
      addPIModelACLine!(net = snet, fromBus = "Mid", toBus = "Load", r_pu = 0.02, x_pu = 0.20, b_pu = 0.0, status = 1)
      return snet
    end

    # SVC at Mid, clamped capacitive: Q_delivered = V^2 * B_max
    svcnet = build_sag_net()
    addShuntVoltageControl!(svcnet; bus = "Mid", target_vm_pu = 1.0, bs_min_mvar = -rating_mvar, bs_max_mvar = rating_mvar)
    run_sparlectra(net = svcnet)
    svc = only(Sparlectra._shunt_controllers(svcnet))
    @test svc.at_limit
    @test svc.bs_mvar == rating_mvar
    v_svc = get_bus_vm_pu(svcnet, "Mid")
    q_svc = v_svc^2 * rating_mvar                     # constant-B limit region
    sh = svcnet.shuntVec[svc.shunt_idx]
    @test imag(sh.y_pu_shunt) * svcnet.baseMVA ≈ rating_mvar atol = 1e-9

    # STATCOM at Mid (needs a machine there), same 1.0-pu rating
    stnet = build_sag_net()
    addProsumer!(net = stnet, busName = "Mid", type = "GENERATOR", p = 0.0, q = 0.0)
    ok, msg = validate!(net = stnet)
    ok || error(msg)
    addMachineVoltageControl!(stnet; bus = "Mid", target_bus = "Load", target_vm_pu = 1.0, s_max_mva = rating_mvar)
    run_control!(stnet)
    st = only([c for c in stnet.machineControls if c isa MachineVoltageControl])
    @test st.at_limit
    v_st = get_bus_vm_pu(stnet, "Mid")
    q_st = st.q_mvar

    # the characteristic split: quadratic vs linear collapse under sag
    @test v_svc < 1.0
    @test v_st < 1.0
    @test q_st ≈ v_st * rating_mvar atol = 2e-2       # linear: V * I_max
    @test q_svc < v_svc * rating_mvar                  # quadratic sits BELOW the linear law at V < 1
    @test q_st / rating_mvar ≈ v_st atol = 2e-3        # normalized delivered Q equals V ...
    @test q_svc / rating_mvar ≈ v_svc^2 atol = 2e-3    # ... versus V^2
  end

  @testset "Voltage deadband is evaluated in pu Vm space" begin
    @test Sparlectra._voltage_within_deadband(1.2009, 1.200, 1e-3)
    @test !Sparlectra._voltage_within_deadband(1.2025, 1.200, 1e-3)
    @test isapprox(Sparlectra._voltage_control_error(1.201, 1.200, :vm), 1e-3; atol = 1e-12)
    @test isapprox(Sparlectra._voltage_control_error(1.201, 1.200, :vm2), 2.401e-3; atol = 1e-12)
  end

  @testset "Voltage tap controller (discrete ratio)" begin
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
  end

  @testset "Branch active power controller (phase)" begin
    # A radial path gives a phase shifter no lever on P (the load dictates
    # the flow); the parallel line provides the loop the controller needs —
    # the honest probe direction (post-#274 flow refresh) exposes that.
    net, tbr = _build_net()
    addPIModelACLine!(net = net, fromBus = "Slack", toBus = "Mid", r_pu = 0.03, x_pu = 0.2, b_pu = 0.0, status = 1)
    result0 = run_sparlectra(net = net, config = _runner_cfg())
    erg0 = result0.numerical_converged ? 0 : 1
    @test erg0 == 0
    p0 = get_branch_p_from_to_mw(net, "Slack", "Mid")

    # one discrete step of 1.25° moves ≈ 8 MW on this loop — the target must
    # sit on the discrete grid within the deadband, otherwise a discrete
    # controller can only oscillate around it
    p_target = p0 - 8.0
    addPowerTransformerControl!(net;
      trafo = string(tbr.branchIdx),
      mode = :branch_active_power,
      target_branch = ("Slack", "Mid"),
      p_target_mw = p_target,
      control_ratio = false,
      control_phase = true,
      is_discrete = true,
      deadband_p_mw = 4.0,
      max_outer_iters = 12,
    )

    result = run_sparlectra(net = net, config = _runner_cfg())
    erg = result.numerical_converged ? 0 : 1
    @test erg == 0
    @test tbr.phase_shift_deg != 0.0
    ctrl = only(Sparlectra._tap_controllers(net))
    @test ctrl.converged
    @test abs(get_branch_p_from_to_mw(net, "Slack", "Mid") - p_target) <= 4.0
  end

  @testset "Disabled tap controllers are skipped" begin
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
  end

  @testset "Disabled control framework preserves successful baseline PF" begin
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
  end

  @testset "Controlled framework result composes converged control status" begin
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
  end

  @testset "Controlled framework result rejects exhausted outer loop" begin
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
  end

  @testset "Direct run_control! returns ControlRunResult" begin
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
  end

  @testset "Direct run_control! resolves default PF config" begin
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
  end

  @testset "Net initializes control_result as nothing" begin
    net = Net(name = "x", baseMVA = 100.0)
    @test net.control_result === nothing
    @test latest_control_result(net) === nothing
  end

  @testset "run_control! stores terminal no-controller result on Net" begin
    net, _ = _build_net()
    result = run_control!(net; controllers = AbstractOuterController[], pf_config = PowerFlowConfig(max_iter = 2), control_config = ControlConfig(enabled = true), verbose = 0)
    @test result.status == :no_controllers
    @test latest_control_result(net) === result
    @test net.control_result === result
  end

  @testset "Disabled control still runs one baseline PF" begin
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
    @test result.last_pf_status == :ok
    @test latest_control_result(net) === result
    @test net.control_result === result
  end

  @testset "All disabled controllers still run one baseline PF" begin
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
  end

  @testset "run_sparlectra net path honors per-call cfg.control" begin
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
  end

  @testset "YAML controller instantiation (#305)" begin
    # Round trip through the real YAML reader: a declared power_transformer
    # controller must be active and reproduce the programmatic setup.
    net_prog, tbr_prog = _build_net()
    addPowerTransformerControl!(net_prog; trafo = string(tbr_prog.branchIdx), mode = :voltage, target_bus = "Load", target_vm_pu = 0.98, control_ratio = true, control_phase = false)
    cfg_prog = SparlectraConfig(powerflow = PowerFlowConfig(method = :rectangular, max_iter = 30, tol = 1e-9), control = ControlConfig(enabled = true, max_outer_iterations = 8, trace = false))
    res_prog = run_sparlectra(net = net_prog, config = cfg_prog)
    @test res_prog.numerical_converged

    net_yaml, tbr_yaml = _build_net()
    cfg_yaml = mktemp() do path, io
      write(
        io,
        """
power_flow:
  max_iter: 30
  tol: 1.0e-9
control:
  enabled: true
  max_outer_iterations: 8
  trace: false
  controllers:
    tap_main:
      type: power_transformer
      trafo: $(tbr_yaml.branchIdx)
      mode: voltage
      target_bus: Load
      target_vm_pu: 0.98
""",
      )
      close(io)
      load_sparlectra_config(path; reload = true)
    end
    @test length(cfg_yaml.control.controllers) == 1
    res_yaml = run_sparlectra(net = net_yaml, config = cfg_yaml)
    @test res_yaml.numerical_converged
    yaml_controllers = collect_outer_controllers(net_yaml)
    @test length(yaml_controllers) == 1
    @test first(yaml_controllers) isa PowerTransformerControl
    # identical outcome: same controller, same solved state
    vm_prog = [something(n._vm_pu, NaN) for n in net_prog.nodeVec]
    vm_yaml = [something(n._vm_pu, NaN) for n in net_yaml.nodeVec]
    @test isapprox(vm_prog, vm_yaml; atol = 1e-9)
    # idempotency: applying the same declaration again adds nothing
    @test applyConfiguredControllers!(net_yaml, cfg_yaml.control) == 0

    # FACTS limit modes through the config surface (#297 Drafts A/F): a
    # STATCOM machine controller (s_max_mva) and an SSSC series controller
    # (v_inj_max_pu) instantiate from declarations, and re-applying skips
    # both. The series re-apply is a REGRESSION test: the idempotency check
    # used the wrong field names (from_bus/to_bus instead of fromBus/toBus)
    # and crashed with a FieldError instead of skipping.
    facts_net = Net(name = "yaml_facts", baseMVA = 100.0)
    for b in ("A", "M1", "M2", "B")
      addBus!(net = facts_net, busName = b, vn_kV = 110.0)
    end
    addProsumer!(net = facts_net, busName = "A", type = "EXTERNALNETWORKINJECTION", referencePri = "A", vm_pu = 1.0, va_deg = 0.0)
    addProsumer!(net = facts_net, busName = "B", type = "ENERGYCONSUMER", p = -60.0, q = -20.0)
    addProsumer!(net = facts_net, busName = "M1", type = "GENERATOR", p = 0.0, q = 0.0)
    addPIModelACLine!(net = facts_net, fromBus = "A", toBus = "M1", r_pu = 0.01, x_pu = 0.10, b_pu = 0.0, status = 1)
    addPIModelACLine!(net = facts_net, fromBus = "M1", toBus = "B", r_pu = 0.01, x_pu = 0.10, b_pu = 0.0, status = 1)
    addPIModelACLine!(net = facts_net, fromBus = "A", toBus = "M2", r_pu = 0.02, x_pu = 0.20, b_pu = 0.0, status = 1)
    addPIModelACLine!(net = facts_net, fromBus = "M2", toBus = "B", r_pu = 0.02, x_pu = 0.20, b_pu = 0.0, status = 1)
    okf, msgf = validate!(net = facts_net)
    okf || error("facts net invalid: $msgf")
    facts_cfg = ControlConfig(
      enabled = true,
      controllers = Any[
        Dict{String,Any}("type" => "machine_voltage", "name" => "statcom_m1", "bus" => "M1", "target_bus" => "B", "target_vm_pu" => 0.97, "s_max_mva" => 25.0),
        Dict{String,Any}("type" => "series_reactance", "name" => "sssc_a_m2", "from_bus" => "A", "to_bus" => "M2", "p_target_mw" => 30.0, "v_inj_max_pu" => 0.05),
      ],
    )
    @test applyConfiguredControllers!(facts_net, facts_cfg) == 2
    fac_machine = only([c for c in facts_net.machineControls if c isa MachineVoltageControl])
    fac_series = only([c for c in facts_net.machineControls if c isa Sparlectra.SeriesReactanceControl])
    @test fac_machine.limit_mode === :current
    @test fac_machine.s_max_mva ≈ 25.0
    @test fac_series.limit_mode === :injected_voltage
    @test fac_series.v_inj_max_pu ≈ 0.05
    @test applyConfiguredControllers!(facts_net, facts_cfg) == 0

    # structural validation fails at config-load time with named entries
    @test_throws ArgumentError ControlConfig(Dict("control" => Dict("controllers" => Dict("x" => Dict("type" => "bogus_type")))))
    @test_throws ArgumentError ControlConfig(Dict("control" => Dict("controllers" => Dict("x" => Dict("type" => "series_reactance")))))
    @test_throws ArgumentError ControlConfig(Dict("control" => Dict("controllers" => Dict("x" => Dict("type" => "shunt_voltage", "bus" => "B", "target_vm_pu" => 1.0, "bs_min_mvar" => -1.0, "bs_max_mvar" => 1.0, "typo_key" => 1)))))
    # reference validation fails at apply time naming the entry
    net_bad, _ = _build_net()
    bad_cfg = ControlConfig(controllers = Any[Dict{String,Any}("type" => "shunt_voltage", "name" => "s1", "bus" => "NoSuchBus", "target_vm_pu" => 1.0, "bs_min_mvar" => -5.0, "bs_max_mvar" => 5.0)])
    err = try
      applyConfiguredControllers!(net_bad, bad_cfg)
      nothing
    catch e
      e
    end
    @test err isa ArgumentError
    @test occursin("control.controllers.s1", sprint(showerror, err))
  end

  @testset "Outer-loop limits are separated from inner PF max_iter" begin
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
  end

  @testset "Transformer controller trace rows include stable keys" begin
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
  end


  @testset "run_sparlectra net entry point" begin
    net, _ = _build_net()
    result_a = run_sparlectra(net = net, config = _runner_cfg(max_iter = 20))
    @test result_a.numerical_converged
    @test result_a.iterations >= 1

  end
  @testset "Tap controller reporting rows and classic section" begin
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
  end

  @testset "Branch derives tap limits from PowerTransformerTaps" begin
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
  end

  @testset "PowerTransformerTaps neutral voltage derivation" begin
    taps_ratio = PowerTransformerTaps(Vn_kV = 20.0, step = 0, lowStep = -2, highStep = 2, neutralStep = 0, voltageIncrement_kV = 0.4, neutralU_ratio = 1.05)
    @test isapprox(taps_ratio.neutralU, 21.0; atol = 1e-12)
    @test isapprox(taps_ratio.neutralU_ratio, 1.05; atol = 1e-12)

    taps_neutral = PowerTransformerTaps(Vn_kV = 20.0, step = 0, lowStep = -2, highStep = 2, neutralStep = 0, voltageIncrement_kV = 0.4, neutralU = 19.5)
    @test isapprox(taps_neutral.neutralU_ratio, 0.975; atol = 1e-12)
  end

  @testset "Controller can be attached during transformer creation" begin
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
  end

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

  @testset "Machine RVC: API validation" begin
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
  end

  @testset "Machine RVC: secant loop reaches a remote target" begin
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
  end

  @testset "Machine RVC: honest at_limit when the target is unreachable" begin
    net = _build_rvc_net(qmin = -2.0, qmax = 2.0)
    addMachineVoltageControl!(net; bus = "GenBus", target_bus = "Load", target_vm_pu = 1.0, deadband_vm_pu = 5e-4)
    pf = PowerFlowConfig(max_iter = 30, tol = 1e-9)
    result = run_control!(net; controllers = collect_outer_controllers(net), pf_config = pf, control_config = ControlConfig(max_outer_iterations = 15), verbose = 0)
    # the loop settles (converged-or-blocked), but the controller itself is
    # parked at its reactive bound with the target still out of reach
    @test result.status in (:converged, :blocked)
    ctrl = Sparlectra._machine_controllers(net)[1]
    @test ctrl.at_limit
    @test !ctrl.converged
    @test isapprox(ctrl.q_mvar, -2.0; atol = 1e-6)
    @test abs(get_bus_vm_pu(net, "Load") - 1.0) > 5e-4
  end
end
