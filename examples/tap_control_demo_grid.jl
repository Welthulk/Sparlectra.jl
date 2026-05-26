using Sparlectra
using Printf

function build_demo_net()
  net = Net(name = "tap_control_demo_grid", baseMVA = 100.0)
  addBus!(net = net, busName = "Slack", vn_kV = 110.0)
  addBus!(net = net, busName = "Mid", vn_kV = 110.0)
  addBus!(net = net, busName = "Load", vn_kV = 110.0)

  addProsumer!(net = net, busName = "Slack", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.02, va_deg = 0.0, referencePri = "Slack")
  addProsumer!(net = net, busName = "Load", type = "ENERGYCONSUMER", p = -70.0, q = -20.0)

  addPIModelTrafo!(net = net, fromBus = "Slack", toBus = "Mid", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, ratio = 1.0, shift_deg = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "Mid", toBus = "Load", r_pu = 0.02, x_pu = 0.12, b_pu = 0.01, status = 1)

  trafo = getNetBranch(net = net, fromBus = "Slack", toBus = "Mid")
  trafo.has_ratio_tap = true
  trafo.tap_min = 0.95
  trafo.tap_max = 1.05
  trafo.tap_step = 0.0125
  return net, trafo
end

function add_demo_controllers!(net::Net, trafo)
  addTapController!(
    net;
    trafo = string(trafo.branchIdx),
    mode = :voltage,
    target_bus = "Load",
    target_vm_pu = 0.98,
    control_ratio = true,
    control_phase = false,
    is_discrete = true,
    deadband_vm_pu = 5e-3,
    max_outer_iters = 8,
  )
end

function print_control_result(result::ControlRunResult)
  println("Control status: ", result.status)
  println("Outer iterations: ", result.outer_iterations)
  println("PF solves: ", result.powerflow_solves)
  println("Controller rows:")
  for row in result.controllers
    println("  ", row)
  end
  println("Trace rows:")
  for row in result.trace
    println("  ", row)
  end
end

function main()
  net, trafo = build_demo_net()

  _, erg0, _ = run_acpflow(net = net, max_ite = 30, tol = 1e-9, method = :rectangular, show_results = false)
  erg0 == 0 || error("Uncontrolled PF failed with erg=$erg0")
  vm_uncontrolled = get_bus_vm_pu(net, "Load")

  add_demo_controllers!(net, trafo)

  _, erg, _ = run_acpflow(net = net, max_ite = 30, tol = 1e-9, method = :rectangular, show_results = false)
  erg == 0 || error("Controlled PF failed with erg=$erg")
  vm_controlled = get_bus_vm_pu(net, "Load")

  result = latest_control_result(net)
  result === nothing && error("Expected control result on net.control_result")

  @printf("Uncontrolled Load Vm: %.6f pu\n", vm_uncontrolled)
  @printf("Controlled Load Vm:   %.6f pu\n", vm_controlled)
  print_control_result(result)
end

Base.invokelatest(main)
