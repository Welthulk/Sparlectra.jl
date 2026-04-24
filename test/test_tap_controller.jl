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
    @test_throws ErrorException addTapController!(net; trafo = "does_not_exist", mode = :voltage, target_bus = "Load", target_vm_pu = 1.0)
    @test_throws ErrorException addTapController!(net; trafo = "1", mode = :voltage)
  end

  @testset "Voltage tap controller (discrete ratio)" begin
    net, tbr = _build_net()
    addTapController!(net;
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
    if net.tapControllers[1].converged
      @test abs(get_bus_vm_pu(net, "Load") - 0.98) <= 0.05
    else
      @test net.tapControllers[1].at_limit
    end
  end

  @testset "Branch active power controller (phase)" begin
    net, tbr = _build_net()
    _, erg0, _ = run_net_acpflow(net = net, max_ite = 30, tol = 1e-9, verbose = 0, method = :rectangular, show_results = false)
    @test erg0 == 0
    p0 = get_branch_p_from_to_mw(net, "Slack", "Mid")

    addTapController!(net;
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

  @testset "Tap controller reporting rows and classic section" begin
    net, tbr = _build_net()
    addTapController!(net;
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
end
