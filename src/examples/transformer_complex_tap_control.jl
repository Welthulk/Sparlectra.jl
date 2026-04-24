using Sparlectra

function run_transformer_complex_tap_control_example(; verbose::Int = 1)
  net = Net(name = "trafo_complex_tap_example", baseMVA = 100.0)
  addBus!(net = net, busName = "Slack", vn_kV = 110.0)
  addBus!(net = net, busName = "Mid", vn_kV = 110.0)
  addBus!(net = net, busName = "Load", vn_kV = 110.0)

  addProsumer!(net = net, busName = "Slack", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.02, va_deg = 0.0, referencePri = "Slack")
  addProsumer!(net = net, busName = "Load", type = "ENERGYCONSUMER", p = -60.0, q = -20.0)

  addPIModelTrafo!(net = net, fromBus = "Slack", toBus = "Mid", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, ratio = 1.0, shift_deg = 0.0)
  addPIModelACLine!(net = net, fromBus = "Mid", toBus = "Load", r_pu = 0.02, x_pu = 0.12, b_pu = 0.01)

  tbr = getNetBranch(net = net, fromBus = "Slack", toBus = "Mid")
  tbr.tap_min = 0.95
  tbr.tap_max = 1.05
  tbr.tap_step = 0.0125
  tbr.phase_min_deg = -10.0
  tbr.phase_max_deg = 10.0
  tbr.phase_step_deg = 1.0

  println("Base case:")
  run_net_acpflow(net = net, max_ite = 30, tol = 1e-9, verbose = verbose, method = :rectangular, show_results = false)
  println("  Vm(Load) = ", get_bus_vm_pu(net, "Load"))
  println("  P(Slack->Mid) = ", get_branch_p_from_to_mw(net, "Slack", "Mid"), " MW")

  addTapController!(net; trafo = string(tbr.branchIdx), mode = :voltage, target_bus = "Load", target_vm_pu = 0.99, control_ratio = true, control_phase = false, is_discrete = true)
  run_net_acpflow(net = net, max_ite = 30, tol = 1e-9, verbose = verbose, method = :rectangular, show_results = false)

  empty!(net.tapControllers)
  addTapController!(net; trafo = string(tbr.branchIdx), mode = :branch_active_power, target_branch = ("Slack", "Mid"), p_target_mw = get_branch_p_from_to_mw(net, "Slack", "Mid") - 3.0, control_ratio = false, control_phase = true, is_discrete = true)
  run_net_acpflow(net = net, max_ite = 30, tol = 1e-9, verbose = verbose, method = :rectangular, show_results = false)

  empty!(net.tapControllers)
  addTapController!(net; trafo = string(tbr.branchIdx), mode = :voltage_and_branch_active_power, target_bus = "Load", target_vm_pu = 0.99, target_branch = ("Slack", "Mid"), p_target_mw = get_branch_p_from_to_mw(net, "Slack", "Mid") - 2.0, control_ratio = true, control_phase = true, is_discrete = true)
  run_net_acpflow(net = net, max_ite = 30, tol = 1e-9, verbose = verbose, method = :rectangular, show_results = true)

  return net
end

Base.invokelatest(() -> run_transformer_complex_tap_control_example())
