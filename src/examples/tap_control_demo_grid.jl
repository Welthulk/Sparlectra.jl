using Sparlectra
using Printf

const MAX_ITE = 40
const PF_TOL = 1e-9
const PF_METHOD = :rectangular

# -----------------------------------------------------------------------------
# Tap controller setpoints
# -----------------------------------------------------------------------------
const T2_TARGET_BUS = "B9"
const T2_TARGET_VM_PU = 1.23

const PST1_TARGET_BRANCH = ("B1", "B2")
const PST1_TARGET_P_MW = 7.0

const ST1_TARGET_BUS = "B5"
const ST1_TARGET_VM_PU = 1.07
const ST1_TARGET_BRANCH = ("B4", "B5")
const ST1_TARGET_P_MW = -60.0

"""
    build_tap_control_demo_grid()

Builds a three-voltage-level network with a meshed 220 kV level, a meshed 110 kV
level, and a radial 20 kV level. The grid contains one OLTC transformer, one PST,
and one coupled ratio+phase regulating transformer.
"""
function build_tap_control_demo_grid()
  net = Net(name = "tap_control_demo_grid", baseMVA = 100.0)

  # 220 kV mesh
  addBus!(net = net, busName = "B1", vn_kV = 220.0)
  addBus!(net = net, busName = "B2", vn_kV = 220.0)
  addBus!(net = net, busName = "B3", vn_kV = 220.0)

  # 110 kV mesh
  addBus!(net = net, busName = "B4", vn_kV = 110.0)
  addBus!(net = net, busName = "B5", vn_kV = 110.0)
  addBus!(net = net, busName = "B6", vn_kV = 110.0)

  # 20 kV radial feeder
  addBus!(net = net, busName = "B7", vn_kV = 20.0)
  addBus!(net = net, busName = "B8", vn_kV = 20.0)
  addBus!(net = net, busName = "B9", vn_kV = 20.0)

  # Injections and loads
  addProsumer!(net = net, busName = "B1", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.02, va_deg = 0.0, referencePri = "B1")
  addProsumer!(net = net, busName = "B3", type = "GENERATOR", p = 95.0, q = 12.0)
  addProsumer!(net = net, busName = "B5", type = "ENERGYCONSUMER", p = -90.0, q = -28.0)
  addProsumer!(net = net, busName = "B8", type = "ENERGYCONSUMER", p = -24.0, q = -8.0)
  addProsumer!(net = net, busName = "B9", type = "ENERGYCONSUMER", p = -42.0, q = -14.0)

  # 220 kV mesh branches (PST on B1->B2)
  addPIModelTrafo!(net = net, fromBus = "B1", toBus = "B2", r_pu = 0.004, x_pu = 0.05, b_pu = 0.0, ratio = 1.0, shift_deg = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "B1", toBus = "B3", r_pu = 0.012, x_pu = 0.10, b_pu = 0.003, status = 1)
  addPIModelACLine!(net = net, fromBus = "B2", toBus = "B3", r_pu = 0.010, x_pu = 0.09, b_pu = 0.003, status = 1)

  # 220/110 kV interconnection
  addPIModelTrafo!(net = net, fromBus = "B3", toBus = "B4", r_pu = 0.006, x_pu = 0.07, b_pu = 0.0, ratio = 1.0, shift_deg = 0.0, status = 1)

  # 110 kV mesh with complex tap transformer on B4->B5
  addPIModelTrafo!(net = net, fromBus = "B4", toBus = "B5", r_pu = 0.005, x_pu = 0.06, b_pu = 0.0, ratio = 1.0, shift_deg = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "B4", toBus = "B6", r_pu = 0.014, x_pu = 0.11, b_pu = 0.004, status = 1)
  addPIModelACLine!(net = net, fromBus = "B6", toBus = "B5", r_pu = 0.013, x_pu = 0.10, b_pu = 0.004, status = 1)

  # 110/20 kV interconnection (OLTC on B6->B7)
  addPIModelTrafo!(net = net, fromBus = "B6", toBus = "B7", r_pu = 0.008, x_pu = 0.08, b_pu = 0.0, ratio = 1.0, shift_deg = 0.0, status = 1)

  # 20 kV radial feeder
  addPIModelACLine!(net = net, fromBus = "B7", toBus = "B8", r_pu = 0.018, x_pu = 0.10, b_pu = 0.001, status = 1)
  addPIModelACLine!(net = net, fromBus = "B8", toBus = "B9", r_pu = 0.025, x_pu = 0.12, b_pu = 0.001, status = 1)

  pst1 = getNetBranch(net = net, fromBus = "B1", toBus = "B2")
  st1 = getNetBranch(net = net, fromBus = "B4", toBus = "B5")
  t2 = getNetBranch(net = net, fromBus = "B6", toBus = "B7")

  # T2 OLTC (ratio only)
  t2.has_ratio_tap = true
  t2.has_phase_tap = false
  t2.tap_min = 0.90
  t2.tap_max = 1.10
  t2.tap_step = 0.0125
  t2.phase_shift_deg = 0.0
  t2.angle = 0.0

  # PST1 (phase only)
  pst1.has_ratio_tap = false
  pst1.has_phase_tap = true
  pst1.tap_ratio = 1.0
  pst1.ratio = 1.0
  pst1.phase_shift_deg = 0.0
  pst1.angle = 0.0
  pst1.phase_min_deg = -15.0
  pst1.phase_max_deg = 15.0
  pst1.phase_step_deg = 1.0

  # ST1 Schrägregler (ratio + phase)
  st1.has_ratio_tap = true
  st1.has_phase_tap = true
  st1.tap_ratio = 1.0
  st1.ratio = 1.0
  st1.tap_min = 0.95
  st1.tap_max = 1.05
  st1.tap_step = 0.00625
  st1.phase_shift_deg = 0.0
  st1.angle = 0.0
  st1.phase_min_deg = -10.0
  st1.phase_max_deg = 10.0
  st1.phase_step_deg = 0.5

  return net
end

function _controller_branches(net::Net)
  return (
    pst1 = getNetBranch(net = net, fromBus = "B1", toBus = "B2"),
    st1 = getNetBranch(net = net, fromBus = "B4", toBus = "B5"),
    t2 = getNetBranch(net = net, fromBus = "B6", toBus = "B7"),
  )
end

"""
    configure_tap_controllers!(net)

Configures one voltage, one phase-power, and one coupled complex tap controller.
All controller setpoints are explicitly wired from top-level constants.
"""
function configure_tap_controllers!(net::Net)
  if !isdefined(Sparlectra, :addTapController!)
    error("Sparlectra.addTapController! is not available. Please update this example to the current API name.")
  end

  empty!(net.tapControllers)
  c = _controller_branches(net)

  addTapController!(net;
    trafo = string(c.t2.branchIdx),
    mode = :voltage,
    target_bus = T2_TARGET_BUS,
    target_vm_pu = T2_TARGET_VM_PU,
    control_ratio = true,
    control_phase = false,
    is_discrete = true,
    deadband_vm_pu = 1e-3,
    max_outer_iters = 20,
  )

  addTapController!(net;
    trafo = string(c.pst1.branchIdx),
    mode = :branch_active_power,
    target_branch = PST1_TARGET_BRANCH,
    p_target_mw = PST1_TARGET_P_MW,
    control_ratio = false,
    control_phase = true,
    is_discrete = true,
    deadband_p_mw = 0.5,
    max_outer_iters = 20,
  )

  addTapController!(net;
    trafo = string(c.st1.branchIdx),
    mode = :voltage_and_branch_active_power,
    target_bus = ST1_TARGET_BUS,
    target_vm_pu = ST1_TARGET_VM_PU,
    target_branch = ST1_TARGET_BRANCH,
    p_target_mw = ST1_TARGET_P_MW,
    control_ratio = true,
    control_phase = true,
    is_discrete = true,
    deadband_vm_pu = 1e-3,
    deadband_p_mw = 0.5,
    max_outer_iters = 24,
  )

  return nothing
end

"""
    print_tap_control_targets(net)

Prints configured target values together with achieved branch/bus values and
current tap states for T2, PST1, and ST1.
"""
function print_tap_control_targets(net::Net)
  c = _controller_branches(net)

  println("Tap Control Targets and Results")
  println("-------------------------------")
  println("T2 voltage controller:")
  @printf("  target bus       : %s\n", T2_TARGET_BUS)
  @printf("  target Vm        : %.4f pu\n", T2_TARGET_VM_PU)
  @printf("  achieved Vm      : %.4f pu\n", get_bus_vm_pu(net, T2_TARGET_BUS))
  @printf("  tap ratio        : %.5f\n", c.t2.tap_ratio)
  @printf("  phase shift      : %.5f deg\n", c.t2.phase_shift_deg)

  println("\nPST1 active power controller:")
  @printf("  target branch    : %s -> %s\n", PST1_TARGET_BRANCH[1], PST1_TARGET_BRANCH[2])
  @printf("  target P         : %.3f MW\n", PST1_TARGET_P_MW)
  @printf("  achieved P       : %.3f MW\n", get_branch_p_from_to_mw(net, PST1_TARGET_BRANCH[1], PST1_TARGET_BRANCH[2]))
  @printf("  tap ratio        : %.5f\n", c.pst1.tap_ratio)
  @printf("  phase shift      : %.5f deg\n", c.pst1.phase_shift_deg)

  println("\nST1 Schraegregler:")
  @printf("  target bus       : %s\n", ST1_TARGET_BUS)
  @printf("  target Vm        : %.4f pu\n", ST1_TARGET_VM_PU)
  @printf("  achieved Vm      : %.4f pu\n", get_bus_vm_pu(net, ST1_TARGET_BUS))
  @printf("  target branch    : %s -> %s\n", ST1_TARGET_BRANCH[1], ST1_TARGET_BRANCH[2])
  @printf("  target P         : %.3f MW\n", ST1_TARGET_P_MW)
  @printf("  achieved P       : %.3f MW\n", get_branch_p_from_to_mw(net, ST1_TARGET_BRANCH[1], ST1_TARGET_BRANCH[2]))
  @printf("  tap ratio        : %.5f\n", c.st1.tap_ratio)
  @printf("  phase shift      : %.5f deg\n", c.st1.phase_shift_deg)

  if isdefined(Sparlectra, :printTapControllerSummary)
    printTapControllerSummary(stdout, net)
  end

  return nothing
end

function _print_section(title::String)
  println("\n======================================================================")
  println(title)
  println("======================================================================")
end

function main()
  println("Tap Control Demo Grid")
  net = build_tap_control_demo_grid()

  _print_section("Uncontrolled power flow")
  _, erg0, _ = run_net_acpflow(net = net, max_ite = MAX_ITE, tol = PF_TOL, verbose = 0, method = PF_METHOD, show_results = true)
  erg0 == 0 || error("Uncontrolled power flow failed with erg=$(erg0)")
  print_tap_control_targets(net)

  configure_tap_controllers!(net)

  _print_section("Configured tap controller setpoints")
  print_tap_control_targets(net)

  _print_section("Controlled power flow with outer-loop tap control")
  _, erg1, _ = run_net_acpflow(net = net, max_ite = MAX_ITE, tol = PF_TOL, verbose = 0, method = PF_METHOD, show_results = true)
  erg1 == 0 || error("Controlled power flow failed with erg=$(erg1)")
  print_tap_control_targets(net)

  return net
end

Base.invokelatest(main)
