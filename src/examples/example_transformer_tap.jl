using Sparlectra

"""
    build_transformer_control_demo_net(; ratio=1.0, shift_deg=0.0)

Create a compact 3-bus demonstration network with one transformer branch and one
line branch.
"""
function build_transformer_control_demo_net(; ratio::Float64 = 1.0, shift_deg::Float64 = 0.0)
  net = Net(name = "trafo_complex_tap_demo", baseMVA = 100.0)
  addBus!(net = net, busName = "Slack", vn_kV = 110.0)
  addBus!(net = net, busName = "Mid", vn_kV = 110.0)
  addBus!(net = net, busName = "Load", vn_kV = 110.0)

  addProsumer!(net = net, busName = "Slack", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.01, va_deg = 0.0, referencePri = "Slack")
  addProsumer!(net = net, busName = "Load", type = "ENERGYCONSUMER", p = -70.0, q = -22.0)

  addPIModelTrafo!(net = net, fromBus = "Slack", toBus = "Mid", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, ratio = ratio, shift_deg = shift_deg, status = 1)
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

function _print_snapshot(net::Net, title::String)
  println("\n", title)
  println("  Vm(Load) [pu]          = ", round(get_bus_vm_pu(net, "Load"), digits = 5))
  println("  P(Slack->Mid) [MW]     = ", round(get_branch_p_from_to_mw(net, "Slack", "Mid"), digits = 4))
  println("  Q(Slack->Mid) [MVAr]   = ", round(get_branch_q_from_to_mvar(net, "Slack", "Mid"), digits = 4))
  tbr = getNetBranch(net = net, fromBus = "Slack", toBus = "Mid")
  println("  tap_ratio / phase_deg  = ", round(tbr.tap_ratio, digits = 5), " / ", round(tbr.phase_shift_deg, digits = 5))
end

"""
    run_example_transformer_complex_tap_control(; verbose=0)

Demonstrates:
1. Regular PF with fixed transformer step settings.
2. Same network without fixed steps (controllable start point).
3. Voltage regulator mode (`:voltage`).
4. Phase/Querregler mode (`:branch_active_power`).
5. Coupled Schrägregler mode (`:voltage_and_branch_active_power`).
"""
function run_example_transformer_complex_tap_control(; verbose::Int = 0)
  println("=== Example: fixed-step PF (no controller) ===")
  net_fixed, _ = build_transformer_control_demo_net(ratio = 1.025, shift_deg = 3.0)
  run_net_acpflow(net = net_fixed, max_ite = 30, tol = 1e-9, verbose = verbose, method = :rectangular, show_results = false)
  _print_snapshot(net_fixed, "Fixed-step operating point")

  println("\n=== Example: same network without fixed steps (controller-enabled) ===")
  net_ctrl, tbr = build_transformer_control_demo_net(ratio = 1.0, shift_deg = 0.0)
  run_net_acpflow(net = net_ctrl, max_ite = 30, tol = 1e-9, verbose = verbose, method = :rectangular, show_results = false)
  _print_snapshot(net_ctrl, "Base point before control")

  println("\n--- Voltage regulator (ratio tap) ---")
  addTapController!(net_ctrl; trafo = string(tbr.branchIdx), mode = :voltage, target_bus = "Load", target_vm_pu = 0.985, control_ratio = true, control_phase = false, is_discrete = true, deadband_vm_pu = 2e-3, max_outer_iters = 12)
  run_net_acpflow(net = net_ctrl, max_ite = 30, tol = 1e-9, verbose = verbose, method = :rectangular, show_results = false)
  _print_snapshot(net_ctrl, "After voltage control")

  println("\n--- Querregler / phase shifter (branch active power) ---")
  empty!(net_ctrl.tapControllers)
  p_ref = get_branch_p_from_to_mw(net_ctrl, "Slack", "Mid")
  addTapController!(net_ctrl; trafo = string(tbr.branchIdx), mode = :branch_active_power, target_branch = ("Slack", "Mid"), p_target_mw = p_ref - 3.0, control_ratio = false, control_phase = true, is_discrete = true, deadband_p_mw = 0.5, max_outer_iters = 12)
  run_net_acpflow(net = net_ctrl, max_ite = 30, tol = 1e-9, verbose = verbose, method = :rectangular, show_results = false)
  _print_snapshot(net_ctrl, "After branch active-power control")

  println("\n--- Schrägregler (coupled voltage + branch active power) ---")
  empty!(net_ctrl.tapControllers)
  p_ref2 = get_branch_p_from_to_mw(net_ctrl, "Slack", "Mid")
  addTapController!(
    net_ctrl;
    trafo = string(tbr.branchIdx),
    mode = :voltage_and_branch_active_power,
    target_bus = "Load",
    target_vm_pu = 0.99,
    target_branch = ("Slack", "Mid"),
    p_target_mw = p_ref2 - 2.0,
    control_ratio = true,
    control_phase = true,
    is_discrete = true,
    deadband_vm_pu = 2e-3,
    deadband_p_mw = 0.5,
    max_outer_iters = 16,
  )
  run_net_acpflow(net = net_ctrl, max_ite = 30, tol = 1e-9, verbose = verbose, method = :rectangular, show_results = false)
  _print_snapshot(net_ctrl, "After coupled Schrägregler control")

  println("\nDone.")
  return net_ctrl
end

Base.invokelatest(() -> run_example_transformer_complex_tap_control())
