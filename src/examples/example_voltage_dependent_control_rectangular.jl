# Example: Voltage-dependent Q(U) and P(U) control in the rectangular solver.

using Sparlectra
using Printf

function run_example_voltage_dependent_control_rectangular(; verbose::Int = 1)
  net = Net(name = "voltage_dependent_control", baseMVA = 100.0)

  addBus!(net = net, busName = "Slack", vn_kV = 110.0)
  addBus!(net = net, busName = "Prosumer", vn_kV = 110.0)

  addACLine!(net = net, fromBus = "Slack", toBus = "Prosumer", length = 30.0, r = 0.05, x = 0.45)

  addProsumer!(net = net, busName = "Slack", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "Slack")

  qu_curve = make_characteristic([(103.4, 35.0), (110.0, 0.0), (116.6, -25.0)]; voltage_unit = :kV, value_unit = :MVAr, vn_kV = 110.0, sbase_MVA = net.baseMVA)
  pu_curve = make_characteristic([(103.4, 25.0), (110.0, 10.0), (116.6, 0.0)]; voltage_unit = :kV, value_unit = :MW, vn_kV = 110.0, sbase_MVA = net.baseMVA)

  addProsumer!(
    net = net,
    busName = "Prosumer",
    type = "SYNCHRONOUSMACHINE",
    p = 10.0,
    q = 0.0,
    qu_controller = QUController(qu_curve; qmin_MVAr = -50.0, qmax_MVAr = 50.0, sbase_MVA = net.baseMVA),
    pu_controller = PUController(pu_curve; pmin_MW = 0.0, pmax_MW = 50.0, sbase_MVA = net.baseMVA),
  )

  addProsumer!(net = net, busName = "Prosumer", type = "ENERGYCONSUMER", p = 45.0, q = 18.0)

  ite, erg = runpf!(net, 40, 1e-9, verbose; method = :rectangular, opt_sparse = true)
  converged = (erg == 0)

  if converged
    V = buildVoltageVector(net)
    Sspec, _, _ = buildControlledSVec(net, V)

    vm = abs(V[2])
    p_ctrl_pu, _ = evaluate_controller(net.prosumpsVec[2].puController, vm)
    q_ctrl_pu, _ = evaluate_controller(net.prosumpsVec[2].quController, vm)

    println("\nConverged in $ite iterations")
    @printf("Bus %-10s |V| = %.5f pu, angle = %.3f deg\n", "Slack", abs(V[1]), rad2deg(angle(V[1])))
    @printf("Bus %-10s |V| = %.5f pu, angle = %.3f deg\n", "Prosumer", abs(V[2]), rad2deg(angle(V[2])))
    @printf("Controlled prosumer setpoints at Vm=%.5f pu: P=%.3f MW, Q=%.3f MVar\n", vm, p_ctrl_pu * net.baseMVA, q_ctrl_pu * net.baseMVA)
    @printf("Net specified injection at Prosumer bus: P=%.3f MW, Q=%.3f MVar\n", real(Sspec[2]) * net.baseMVA, imag(Sspec[2]) * net.baseMVA)
    println("\nDetailed results:")
    printACPFlowResults(net, 0.0, ite, 1e-9; converged = true, solver = :rectangular)
  else
    println("Power flow did not converge.")
  end

  return net
end

Base.invokelatest(run_example_voltage_dependent_control_rectangular)
