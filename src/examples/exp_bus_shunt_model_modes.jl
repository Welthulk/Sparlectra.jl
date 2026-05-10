# Demonstrates the two bus-shunt modeling modes.
using Sparlectra

"""
    build_demo_net(; bus_shunt_model="admittance") -> Net

Build a minimal two-bus system with one MATPOWER-style bus shunt. The shunt
values are MW/MVAr at 1.0 pu on the system base. In `"admittance"` mode they
are stamped into Ybus; in `"voltage_dependent_injection"` mode they are kept in
voltage-dependent injection terms.
"""
function build_demo_net(; bus_shunt_model = "admittance")
  net = Net(name = "bus_shunt_model_demo", baseMVA = 100.0, bus_shunt_model = bus_shunt_model)
  addBus!(net = net, busName = "Slack", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
  addBus!(net = net, busName = "Load", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
  addACLine!(net = net, fromBus = "Slack", toBus = "Load", length = 1.0, r = 0.02, x = 0.06)
  addProsumer!(net = net, busName = "Slack", type = "EXTERNALNETWORKINJECTION", referencePri = "Slack")
  addProsumer!(net = net, busName = "Load", type = "ENERGYCONSUMER", p = 50.0, q = 15.0)
  addShunt!(net = net, busName = "Load", pShunt = 2.0, qShunt = 6.0)
  return net
end

"""
    main() -> nothing

Solve the demo network in both bus-shunt modes and print the Y-bus diagonal and
reported shunt power. The default package behavior is `"admittance"`.
"""
function main()
  for mode in ("admittance", "voltage_dependent_injection")
    net = build_demo_net(bus_shunt_model = mode)
    Y = createYBUS(net = net, sparse = false, printYBUS = false)
    iterations, status, _ = run_net_acpflow(net = net, max_ite = 30, tol = 1e-9, verbose = 0, method = :rectangular, show_results = false)
    sh = net.shuntVec[1]
    println("mode = ", mode)
    println("  converged = ", status == 0, ", iterations = ", iterations)
    println("  Ybus[Load,Load] before solve = ", Y[2, 2])
    println("  shunt power after solve = ", sh.p_shunt, " MW, ", sh.q_shunt, " MVAr")
  end
  return nothing
end

if get(ENV, "SPARLECTRA_SUITE_NO_AUTORUN", "0") != "1"
  Base.invokelatest(main)
end
