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

# file: test/test_legacy_tap_controller.jl
# purpose: legacy profile: the testsets of test/test_tap_controller.jl that pin an earlier
#          change (an issue, a bug, a bit-identity or off-path guarantee).
#          They moved here unchanged from test_tap_controller.jl; the helpers they use
#          are defined in that file, which the legacy group includes first.
using Sparlectra
using Test

function run_legacy_tap_controller_tests()
  # cut from the runner function
  (function ()
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
    _runner_cfg(; max_iter::Int = 30, tol::Float64 = 1e-9, control::ControlConfig = ControlConfig()) = SparlectraConfig(powerflow = PowerFlowConfig(max_iter = max_iter, tol = tol), output = OutputConfig(logfile_results = :off), control = control)
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

    @testset "PST reactance coupling X(alpha) (#274)" begin (function ()
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
      probe_cfg = PowerFlowConfig(max_iter = 30, tol = 1e-8, start_mode = Sparlectra.StartModeConfig(flatstart = true))
      net3, tbr3 = build_xalpha(model = formula_model())
      run_sparlectra(net = net3, config = _runner_cfg())
      addPowerTransformerControl!(net3; trafo = string(tbr3.branchIdx), mode = :branch_active_power, target_branch = ("Slack", "Mid"), p_target_mw = 0.0, control_ratio = false, control_phase = true)
      ctrl3 = only(Sparlectra._tap_controllers(net3))
      x_before = tbr3.x_pu
      phi_before = tbr3.phase_shift_deg
      dir_model = Sparlectra._phase_probe_direction(net3, tbr3, ctrl3, probe_cfg)
      @test tbr3.x_pu == x_before
      @test tbr3.phase_shift_deg == phi_before
      net4, tbr4 = build_xalpha()
      run_sparlectra(net = net4, config = _runner_cfg())
      addPowerTransformerControl!(net4; trafo = string(tbr4.branchIdx), mode = :branch_active_power, target_branch = ("Slack", "Mid"), p_target_mw = 0.0, control_ratio = false, control_phase = true)
      ctrl4 = only(Sparlectra._tap_controllers(net4))
      dir_static = Sparlectra._phase_probe_direction(net4, tbr4, ctrl4, probe_cfg)
      @test dir_model == dir_static
      @test abs(dir_model) == 1.0
    end)() end

    @testset "SVC shunt voltage controller + controllable elements (#227)" begin (function ()
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
    end)() end

    @testset "Master/slave transformer group (#322)" begin (function ()
      # substation fixture: two identical parallel transformers HV -> LV
      build_substation = function ()
        gnet = Net(name = "tap_group", baseMVA = 100.0)
        addBus!(net = gnet, busName = "HV", vn_kV = 110.0)
        addBus!(net = gnet, busName = "LV", vn_kV = 20.0)
        addBus!(net = gnet, busName = "Feeder", vn_kV = 20.0)
        addProsumer!(net = gnet, busName = "HV", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.02, va_deg = 0.0, referencePri = "HV")
        addProsumer!(net = gnet, busName = "Feeder", type = "ENERGYCONSUMER", p = 40.0, q = 12.0)
        addPIModelTrafo!(net = gnet, fromBus = "HV", toBus = "LV", r_pu = 0.006, x_pu = 0.12, b_pu = 0.0, ratio = 1.0, shift_deg = 0.0, status = 1)
        addPIModelTrafo!(net = gnet, fromBus = "HV", toBus = "LV", r_pu = 0.006, x_pu = 0.12, b_pu = 0.0, ratio = 1.0, shift_deg = 0.0, status = 1)
        addPIModelACLine!(net = gnet, fromBus = "LV", toBus = "Feeder", r_pu = 0.01, x_pu = 0.05, b_pu = 0.0, status = 1)
        for br in gnet.branchVec
          if br.ratio != 0.0
            br.has_ratio_tap = true
            br.tap_min = 0.9
            br.tap_max = 1.1
            br.tap_step = 0.0125
          end
        end
        ok, msg = validate!(net = gnet)
        ok || error("test net invalid: $msg")
        ids = [string(br.branchIdx) for br in gnet.branchVec if br.ratio != 0.0]
        return gnet, ids
      end

      # the physics this feature prevents: misaligned FIXED taps on parallel
      # units drive a circulating reactive power through the loop
      cnet, _ = build_substation()
      ct = [br for br in cnet.branchVec if br.ratio != 0.0]
      ct[1].tap_ratio = 1.05
      ct[1].ratio = 1.05
      ct[2].tap_ratio = 0.95
      ct[2].ratio = 0.95
      _, cerg = runpf!(cnet, 30, 1e-8, 0)
      @test cerg == 0
      calcNetLosses!(cnet)
      @test abs(ct[1].fBranchFlow.qFlow - ct[2].fBranchFlow.qFlow) > 30.0   # massive Q split

      # group control: master steps toward a target that REQUIRES movement,
      # the follower mirrors step-synchronously, the loop stays circulation-free
      # NOTE the deadband: synchronized steps double the voltage effect per
      # master step, so a group needs a deadband of at least half its
      # AGGREGATED step effect (documented in control_framework.md); the
      # single-unit default of 1e-3 pu is finer than this group can resolve
      gnet, ids = build_substation()
      addPowerTransformerControl!(gnet; trafo = ids[1], followers = [ids[2]], mode = :voltage, target_bus = "Feeder", target_vm_pu = 1.03, deadband_vm_pu = 5e-3)
      gres = run_control!(gnet)
      @test gres.converged
      b1 = Sparlectra._find_trafo_branch(gnet, ids[1])
      b2 = Sparlectra._find_trafo_branch(gnet, ids[2])
      @test b1.tap_ratio != 1.0                                  # the master actually moved
      @test isapprox(b1.tap_ratio, b2.tap_ratio; atol = 1e-9)    # taps stay aligned
      @test abs(get_bus_vm_pu(gnet, "Feeder") - 1.03) <= 5e-3 + 1e-9
      calcNetLosses!(gnet)
      tq1 = b1.fBranchFlow.qFlow
      tq2 = b2.fBranchFlow.qFlow
      @test abs(tq1 - tq2) < 1.0                                 # near-equal sharing, no circulation
      # element view names the group
      gels = controllableElements(gnet)
      @test occursin("+1 follower", only(filter(e -> e.actuator in (:tap_ratio, :tap_ratio_and_phase_shift), gels)).element)

      # validation: self-follow, occupied follower, double-follow, and a
      # controller on a bus that already follows
      vnet, vids = build_substation()
      @test_throws ErrorException addPowerTransformerControl!(vnet; trafo = vids[1], followers = [vids[1]], mode = :voltage, target_bus = "Feeder", target_vm_pu = 1.0)
      addPowerTransformerControl!(vnet; trafo = vids[1], followers = [vids[2]], mode = :voltage, target_bus = "Feeder", target_vm_pu = 1.0)
      @test_throws ErrorException addPowerTransformerControl!(vnet; trafo = vids[2], mode = :voltage, target_bus = "Feeder", target_vm_pu = 1.0)
      wnet2, wids = build_substation()
      @test_throws ErrorException addPowerTransformerControl!(wnet2; trafo = wids[1], followers = [wids[2], wids[2]], mode = :voltage, target_bus = "Feeder", target_vm_pu = 1.0)
      @test_throws ErrorException addPowerTransformerControl!(wnet2; trafo = wids[1], followers = [wids[2]], mode = :branch_active_power, target_branch = ("HV", "LV"), p_target_mw = 10.0, control_ratio = false, control_phase = true)

      # config surface: followers list instantiates the group
      ynet, yids = build_substation()
      ycfg = ControlConfig(enabled = true, controllers = Any[Dict{String,Any}("type" => "power_transformer", "trafo" => yids[1], "followers" => Any[yids[2]], "mode" => "voltage", "target_bus" => "Feeder", "target_vm_pu" => 1.0)])
      @test applyConfiguredControllers!(ynet, ycfg) == 1
      yc = only(Sparlectra._tap_controllers(ynet))
      @test yc.followers == [yids[2]]
    end)() end

    @testset "MSC/MSR discrete shunt bank (#324)" begin (function ()
      build_sag = function ()
        snet = Net(name = "msc_ctrl", baseMVA = 100.0)
        for bus in ("Slack", "Mid", "Load")
          addBus!(net = snet, busName = bus, vn_kV = 110.0)
        end
        addProsumer!(net = snet, busName = "Slack", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "Slack")
        addProsumer!(net = snet, busName = "Load", type = "LOAD", p = 60.0, q = 25.0)
        addPIModelACLine!(net = snet, fromBus = "Slack", toBus = "Mid", r_pu = 0.02, x_pu = 0.20, b_pu = 0.0, status = 1)
        addPIModelACLine!(net = snet, fromBus = "Mid", toBus = "Load", r_pu = 0.02, x_pu = 0.20, b_pu = 0.0, status = 1)
        ok, msg = validate!(net = snet)
        ok || error("test net invalid: $msg")
        return snet
      end

      # registration validation: bad step, no admissible block in the range,
      # start value snapped onto the grid
      vnet = build_sag()
      @test_throws ErrorException addShuntVoltageControl!(vnet; bus = "Mid", target_vm_pu = 1.0, bs_min_mvar = -40.0, bs_max_mvar = 40.0, step_mvar = -5.0)
      @test_throws ErrorException addShuntVoltageControl!(vnet; bus = "Mid", target_vm_pu = 1.0, bs_min_mvar = 3.0, bs_max_mvar = 7.0, step_mvar = 10.0)
      snap = addShuntVoltageControl!(vnet; bus = "Mid", target_vm_pu = 1.0, bs_min_mvar = -40.0, bs_max_mvar = 40.0, step_mvar = 10.0, bs_start_mvar = 13.0)
      @test snap.bs_mvar == 10.0                       # 13 snapped to the nearest block
      @test snap.name == "MSC_Mid"
      clearShuntControllers!(vnet)

      # park behavior: the bank steps toward the target in whole blocks,
      # never overshoots (trunc semantics), and parks on the last step
      # BEFORE crossing when no whole block improves further — without
      # exhausting the outer-iteration budget (anti-hunting)
      pnet = build_sag()
      ctrl = addShuntVoltageControl!(pnet; bus = "Mid", target_vm_pu = 0.95, bs_min_mvar = -40.0, bs_max_mvar = 40.0, step_mvar = 10.0)
      pres = run_control!(pnet)
      @test ctrl.bs_mvar % 10.0 == 0.0
      @test !ctrl.converged
      @test ctrl.parked
      @test ctrl.status == :parked
      @test !ctrl.at_limit
      @test pres.outer_iterations < ctrl.max_outer_iters       # parked, not exhausted
      @test get_bus_vm_pu(pnet, "Mid") < 0.95                  # approached from below, never crossed

      # outermost block: an unreachable target clamps on the last block and
      # reports at_limit (the constant-B region applies unchanged)
      lnet = build_sag()
      lctrl = addShuntVoltageControl!(lnet; bus = "Mid", target_vm_pu = 1.0, bs_min_mvar = -20.0, bs_max_mvar = 20.0, step_mvar = 10.0)
      run_control!(lnet)
      @test lctrl.at_limit
      @test lctrl.status == :at_limit
      @test lctrl.bs_mvar == 20.0

      # element view and report rows carry the bank vocabulary
      els = controllableElements(pnet)
      @test length(els) == 1
      @test els[1].device == "MSC/MSR (switched shunt bank)"
      @test els[1].discrete === true
      prow = only(Sparlectra.control_report_rows(ctrl, pnet, Sparlectra.NoControlState(), (outer_iteration = 0,)))
      @test prow.step_mvar ≈ 10.0
      @test prow.step_position == round(Int, ctrl.bs_mvar / 10.0)
      @test prow.parked === true

      # continuous SVC untouched: inert mode field, unchanged vocabulary
      cnet = build_sag()
      cctrl = addShuntVoltageControl!(cnet; bus = "Mid", target_vm_pu = 0.95, bs_min_mvar = -40.0, bs_max_mvar = 40.0)
      run_control!(cnet)
      @test cctrl.step_mvar === nothing
      @test cctrl.converged
      cels = controllableElements(cnet)
      @test cels[1].device == "SVC (variable shunt)"
      @test cels[1].discrete === false

      # config surface: step_mvar instantiates the discrete bank, re-apply
      # skips (idempotency)
      ynet = build_sag()
      ycfg = ControlConfig(enabled = true, controllers = Any[Dict{String,Any}("type" => "shunt_voltage", "name" => "bank_mid", "bus" => "Mid", "target_vm_pu" => 0.95, "bs_min_mvar" => -40.0, "bs_max_mvar" => 40.0, "step_mvar" => 10.0)])
      @test applyConfiguredControllers!(ynet, ycfg) == 1
      yctrl = only(Sparlectra._shunt_controllers(ynet))
      @test yctrl.step_mvar ≈ 10.0
      @test applyConfiguredControllers!(ynet, ycfg) == 0

      # A bank on its top block next to a tap controller on the same bus,
      # both acting in every pass: bank target 0.972 pu below the tap band
      # 0.975 .. 0.985 pu. The bus starts below the bank band, so the bank
      # stands at its limit while the tap steps up. The last tap step lifts
      # the bus over the bank band (0.967 .. 0.977 pu). The bank is then not
      # at its limit any more (the error points inward) and has to switch
      # back. Before the fix the flag of the earlier pass stayed set and the
      # run ended :converged with the bank on its top block.
      function bank_tap_net(name)
        bnet = Net(name = name, baseMVA = 100.0)
        for bus in ("Slack", "Mid", "Load")
          addBus!(net = bnet, busName = bus, vn_kV = 110.0)
        end
        addProsumer!(net = bnet, busName = "Slack", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.02, va_deg = 0.0, referencePri = "Slack")
        addProsumer!(net = bnet, busName = "Load", type = "LOAD", p = 90.0, q = 35.0)
        addPIModelTrafo!(net = bnet, fromBus = "Slack", toBus = "Mid", r_pu = 0.01, x_pu = 0.08, b_pu = 0.00, ratio = 1.0, shift_deg = 0.0, status = 1)
        addPIModelACLine!(net = bnet, fromBus = "Mid", toBus = "Load", r_pu = 0.02, x_pu = 0.12, b_pu = 0.01, status = 1)
        tbr = getNetBranch(net = bnet, fromBus = "Slack", toBus = "Mid")
        tbr.has_ratio_tap = true
        (tbr.tap_min, tbr.tap_max, tbr.tap_step) = (0.90, 1.10, 0.0125)
        return bnet, tbr
      end
      b7net, b7tbr = bank_tap_net("bank_limit")
      addPowerTransformerControl!(b7net; trafo = string(b7tbr.branchIdx), mode = :voltage, target_bus = "Load", target_vm_pu = 0.98, deadband_vm_pu = 0.005)
      b7bank = @test_logs (:warn, r"already regulates the voltage of bus Load") addShuntVoltageControl!(b7net; bus = "Load", target_vm_pu = 0.972, bs_min_mvar = 0.0, bs_max_mvar = 10.0, step_mvar = 5.0, deadband_vm_pu = 0.005, bs_start_mvar = 10.0)
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
        pknet, pktbr = bank_tap_net("bank_park")
        addPowerTransformerControl!(pknet; trafo = string(pktbr.branchIdx), mode = :voltage, target_bus = "Load", target_vm_pu = 1.0, deadband_vm_pu = 0.005)
        pkbank = @test_logs (:warn, r"already regulates the voltage of bus Load") addShuntVoltageControl!(pknet; bus = "Load", target_vm_pu = 0.96, bs_min_mvar = 0.0, bs_max_mvar = 10.0, step_mvar = 5.0, deadband_vm_pu = 0.01, bs_start_mvar = start)
        pkres = run_control!(pknet)
        pkrows = filter(r -> r.controller_type == "ShuntVoltageControl", pkres.trace)
        @test any(r -> r.bs_mvar == 10.0, pkrows)
        @test !(pkres.status === :converged && pkbank.parked && abs(get_bus_vm_pu(pknet, "Load") - 0.96) > 0.01)
      end
    end)() end

    @testset "STATCOM current-limit mode (#297 Draft A)" begin (function ()
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

      # priority: a bank at rank 1 waits until the STATCOM stands at its
      # current limit and switches in the next pass. At its limit the STATCOM
      # keeps following V * S_max; those moves do not hold the bank, so the
      # run converges within the default budget of 20 passes
      pnet = build_weak_net()
      addMachineVoltageControl!(pnet; bus = "Mid", target_bus = "Load", target_vm_pu = 1.0, s_max_mva = 12.0)
      pbank = addShuntVoltageControl!(pnet; bus = "Load", target_vm_pu = 1.0, bs_min_mvar = 0.0, bs_max_mvar = 60.0, step_mvar = 5.0, deadband_vm_pu = 0.02, priority = 1)
      pres = run_control!(pnet)
      @test pres.status === :converged
      @test pbank.bs_mvar > 0.0
      first_bank = minimum(r.outer_iteration for r in pres.trace if r.controller_type == "ShuntVoltageControl")
      first_statcom_at_limit = minimum(r.outer_iteration for r in pres.trace if r.controller_type == "MachineVoltageControl" && r.at_limit)
      @test first_bank == first_statcom_at_limit + 1
    end)() end

    @testset "SVC vs STATCOM limit characteristic (#297 Draft E)" begin (function ()
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
    end)() end

    @testset "Control loop console output once per run (#387)" begin (function ()
      # the inner solver printed its full diagnostic set on every control
      # pass; now the first pass prints it, later passes get one line each,
      # the wrong-branch block stays silent with the detection off, and the
      # result header names the passes next to the last-pass iteration count
      function _controlled_net()
        net, tbr = _build_net()
        addPowerTransformerControl!(net; trafo = string(tbr.branchIdx), mode = :voltage, target_bus = "Load", target_vm_pu = 0.98, control_ratio = true, control_phase = false, is_discrete = true, deadband_vm_pu = 5e-3, max_outer_iters = 8)
        return net
      end
      function _run_captured(net, control::ControlConfig; wrong_branch::Symbol = :warn)
        pf_cfg = PowerFlowConfig(max_iter = 30, tol = 1e-9, wrong_branch_detection = wrong_branch)
        # redirect_stdout needs a real stream, not an IOBuffer
        out_path = tempname()
        logger = Test.TestLogger()
        res = open(out_path, "w") do io
          redirect_stdout(io) do
            Base.with_logger(logger) do
              run_control!(net; pf_config = pf_cfg, control_config = control, verbose = 1)
            end
          end
        end
        return res, read(out_path, String), logger.logs
      end
      net = _controlled_net()
      res, out, logs = _run_captured(net, ControlConfig())
      # the loop may end at the tap limit; what matters here is that it ran
      # several passes with the power flow solving each time
      @test res.last_pf_status === :ok
      @test res.powerflow_solves > 1
      @test res.total_pf_iterations >= res.last_pf_iterations
      @test count("Q-Limit Active-Set Summary", out) == 1
      @test count("rectangular convergence:", out) == 1
      @test count("Wrong-branch check:", out) == 1
      # one line per solve after the first; a run that ends :blocked closes
      # with one closing solve and its own line
      closing_lines = count("control closing solve after pass ", out)
      @test closing_lines == (res.status === :blocked ? 1 : 0)
      @test count("control pass ", out) + closing_lines == res.powerflow_solves - 1
      @test occursin("no active-set change", out)
      @test count(l -> occursin("start projection selected", string(l.message)), logs) <= 1
      # detection off: no wrong-branch block at all, not even on the first pass
      net = _controlled_net()
      res, out, logs = _run_captured(net, ControlConfig(); wrong_branch = :off)
      @test res.last_pf_status === :ok
      @test !occursin("Wrong-branch check:", out)
      # the full per-pass output stays available for debugging the loop
      net = _controlled_net()
      res, out, logs = _run_captured(net, ControlConfig(verbose_passes = true))
      @test count("Q-Limit Active-Set Summary", out) == res.powerflow_solves
      # log_iterations off: quiet passes without the summary line
      net = _controlled_net()
      res, out, logs = _run_captured(net, ControlConfig(log_iterations = false))
      @test !occursin("control pass ", out)
      # result header: passes and the total next to the last-pass count
      net = _controlled_net()
      result = run_sparlectra(net = net, config = _runner_cfg())
      hdr_path = tempname()
      open(hdr_path, "w") do io
        redirect_stdout(io) do
          printACPFlowResults(net, result.elapsed_s, result.iterations, 1e-9, false, ""; converged = result.final_converged)
        end
      end
      hdr_text = read(hdr_path, String)
      @test occursin("Control passes :", hdr_text)
      @test occursin("last pass", hdr_text)
      @test occursin(string(net.control_result.total_pf_iterations, " inner iterations in total"), hdr_text)

      # the AC island table once per run on a net with two AC islands, also
      # when control passes and tap probes solve it several times (every solve
      # printed it); with the full diagnostics the run's first solve prints
      # it. A direct solve prints it only with verbose.
      function _island_output(f)
        out_path = tempname()
        # the failing auto stages log their start projection; kept off the console
        res = open(io -> redirect_stdout(() -> Base.with_logger(f, Test.TestLogger()), io), out_path, "w")
        return res, read(out_path, String)
      end
      function _island_tables(f)
        res, out = _island_output(f)
        return res, count("AC island detection:", out)
      end
      for diag in (:compact, :full)
        inet, itbr = _build_two_island_net()
        addPowerTransformerControl!(inet; trafo = string(itbr.branchIdx), mode = :voltage, target_bus = "Load", target_vm_pu = 1.035, control_ratio = true, control_phase = false, is_discrete = true, deadband_vm_pu = 5e-3, max_outer_iters = 8)
        icfg = SparlectraConfig(powerflow = PowerFlowConfig(max_iter = 30, tol = 1e-9), output = OutputConfig(logfile_results = :off, console_diagnostics = diag))
        ires, ntab = _island_tables(() -> run_sparlectra(net = inet, config = icfg))
        @test ires.final_converged
        @test inet.control_result.powerflow_solves > 1
        @test ntab == 1
      end
      for (v, expected) in ((0, 0), (1, 1))
        inet, _ = _build_two_island_net()
        (_, ierg), ntab = _island_tables(() -> runpf!(inet, 30, 1e-9, v; islands_enabled = true))
        @test ierg == 0
        @test ntab == expected
      end
      # an auto run calls the execution once per escalation stage on the same
      # net and still prints one table (max_iter = 2 and tol = 1e-12 make
      # every stage fail)
      inet, _ = _build_two_island_net()
      acfg = SparlectraConfig(powerflow = PowerFlowConfig(mode = :auto, max_iter = 2, tol = 1e-12), output = OutputConfig(logfile_results = :off))
      _, ntab = _island_tables(() -> run_sparlectra(net = inet, config = acfg))
      @test count(a -> !isempty(a.solver), Sparlectra.auto_pf_record(inet).attempts) > 1
      @test ntab == 1
      # auto_slack on two islands without a slack: the run promotes the unit
      # with the stated reference priority on bus 1, and the table names that
      # bus as the reference of island 1, not the voltage-regulating unit on
      # bus 2 that the island detection would pick without the promotion
      for diag in (:compact, :full)
        snet = Net(name = "two_islands_no_slack", baseMVA = 100.0)
        for bus in ("A", "B", "C", "D", "E")
          addBus!(net = snet, busName = bus, vn_kV = 110.0)
        end
        addProsumer!(net = snet, busName = "A", type = "SYNCHRONOUSMACHINE", p = 50.0, q = 0.0, referencePriority = 2)
        addProsumer!(net = snet, busName = "B", type = "SYNCHRONOUSMACHINE", p = 10.0, vm_pu = 1.01)
        addProsumer!(net = snet, busName = "C", type = "ENERGYCONSUMER", p = -55.0, q = -10.0)
        addPIModelACLine!(net = snet, fromBus = "A", toBus = "C", r_pu = 0.01, x_pu = 0.05, b_pu = 0.0, status = 1)
        addPIModelACLine!(net = snet, fromBus = "B", toBus = "C", r_pu = 0.01, x_pu = 0.05, b_pu = 0.0, status = 1)
        addProsumer!(net = snet, busName = "D", type = "SYNCHRONOUSMACHINE", p = 20.0, vm_pu = 1.0)
        addProsumer!(net = snet, busName = "E", type = "ENERGYCONSUMER", p = -20.0, q = -5.0)
        addPIModelACLine!(net = snet, fromBus = "D", toBus = "E", r_pu = 0.01, x_pu = 0.05, b_pu = 0.0, status = 1)
        scfg = SparlectraConfig(powerflow = PowerFlowConfig(max_iter = 30, tol = 1e-9, auto_slack = true), output = OutputConfig(logfile_results = :off, console_diagnostics = diag))
        sres, out = _island_output(() -> run_sparlectra(net = snet, config = scfg))
        @test sres.final_converged
        @test snet.slackVec == [1]
        @test count("AC island detection:", out) == 1
        @test occursin("island 1: buses=3, branches=2, ref=1, status=has_ref", out)
        @test findfirst("auto_slack:", out).start < findfirst("AC island detection:", out).start
      end
    end)() end

    @testset "YAML controller instantiation (#305)" begin (function ()
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
  config_version: 1
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
    end)() end

  end)()
  return true
end
