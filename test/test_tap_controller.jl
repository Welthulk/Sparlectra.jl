# file: test/test_tap_controller.jl

function run_tap_controller_tests()
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

  @testset "Voltage deadband is evaluated in pu Vm space" begin
    @test Sparlectra._voltage_within_deadband(1.2009, 1.200, 1e-3)
    @test !Sparlectra._voltage_within_deadband(1.2025, 1.200, 1e-3)
    @test isapprox(Sparlectra._voltage_control_error(1.201, 1.200, :vm), 1e-3; atol = 1e-12)
    @test isapprox(Sparlectra._voltage_control_error(1.201, 1.200, :vm2), 2.401e-3; atol = 1e-12)
  end

  @testset "Voltage tap controller (discrete ratio)" begin
    net, tbr = _build_net()
    _, erg0, _ = run_net_acpflow(net = net, max_ite = 30, tol = 1e-9, verbose = 0, method = :rectangular, show_results = false)
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
    _, erg, _ = run_net_acpflow(net = net, max_ite = 30, tol = 1e-9, verbose = 0, method = :rectangular, show_results = false)
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
    net, tbr = _build_net()
    _, erg0, _ = run_net_acpflow(net = net, max_ite = 30, tol = 1e-9, verbose = 0, method = :rectangular, show_results = false)
    @test erg0 == 0
    p0 = get_branch_p_from_to_mw(net, "Slack", "Mid")

    addPowerTransformerControl!(net;
      trafo = string(tbr.branchIdx),
      mode = :branch_active_power,
      target_branch = ("Slack", "Mid"),
      p_target_mw = p0 - 5.0,
      control_ratio = false,
      control_phase = true,
      is_discrete = true,
      deadband_p_mw = 0.8,
      max_outer_iters = 6,
    )

    _, erg, _ = run_net_acpflow(net = net, max_ite = 30, tol = 1e-9, verbose = 0, method = :rectangular, show_results = false)
    @test erg == 0
    @test tbr.phase_shift_deg != 0.0
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
    _, erg, _ = run_net_acpflow(net = net, max_ite = 30, tol = 1e-9, verbose = 0, method = :rectangular, show_results = false)
    @test erg == 0
    @test tbr.tap_ratio == ratio_before
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
    _, erg, _ = run_net_acpflow(net = net, max_ite = 30, tol = 1e-9, verbose = 0, method = :rectangular, show_results = false)
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
end
