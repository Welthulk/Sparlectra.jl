# Example: Voltage-dependent Q(U) and P(U) control in the rectangular solver.

using Sparlectra
using Printf

function _maybe_plot_characteristic(ch::PiecewiseLinearCharacteristic; enable_plot::Bool = false, out_file::String = "kink_characteristic.png")
  enable_plot || return

  if isnothing(Base.find_package("Plots"))
    @info "Skipping curve plot: package Plots.jl is not available in the current/global Julia environment."
    @info "Install globally if desired: julia -e 'using Pkg; Pkg.add(\"Plots\")'"
    return
  end

  @eval using Plots
  xs = range(ch.points[1][1], ch.points[end][1]; length = 200)
  ys = [first(evaluate_characteristic(ch, x)) for x in xs]
  xk = [p[1] for p in ch.points]
  yk = [p[2] for p in ch.points]

  plt = plot(xs, ys; lw = 2, xlabel = "Voltage [pu]", ylabel = "P [pu]", label = "Characteristic", title = "Piecewise-linear kink characteristic")
  scatter!(plt, xk, yk; ms = 4, label = "Kink points")
  savefig(plt, out_file)
  @info "Characteristic curve plot written to $(abspath(out_file))"
end

function _run_kink_check()
  kink_curve = make_characteristic([(104.5, 20.0), (110.0, 10.0), (115.5, 20.0)]; voltage_unit = :kV, value_unit = :MW, vn_kV = 110.0, sbase_MVA = 100.0)

  _, slope_left = evaluate_characteristic(kink_curve, 1.0 - 1e-8)
  value_kink, slope_kink = evaluate_characteristic(kink_curve, 1.0)
  _, slope_right = evaluate_characteristic(kink_curve, 1.0 + 1e-8)

  @assert isapprox(value_kink, 0.1; atol = 1e-12)
  @assert slope_left < 0.0
  @assert slope_right > 0.0
  @assert isapprox(slope_kink, slope_left; atol = 1e-12)

  return (curve = kink_curve, slope_left = slope_left, slope_kink = slope_kink, slope_right = slope_right)
end

function run_example_voltage_dependent_control_rectangular(; verbose::Int = 1, plot_curve::Bool = false)
  kink_check = _run_kink_check()
  _maybe_plot_characteristic(kink_check.curve; enable_plot = plot_curve)

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

  ite, erg, etime = run_net_acpflow(
    net = net,
    max_ite = 40,
    tol = 1e-9,
    verbose = verbose,
    opt_sparse = true,
    method = :rectangular,
    show_results = false,
  )
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
    @printf("Kink check (dP/dU): left=%.6f, at kink=%.6f, right=%.6f\n", kink_check.slope_left, kink_check.slope_kink, kink_check.slope_right)
    println("\nDetailed results:")
    printACPFlowResults(net, etime, ite, 1e-9; converged = true, solver = :rectangular)
  else
    println("Power flow did not converge.")
  end

  return net
end

Base.invokelatest(run_example_voltage_dependent_control_rectangular)
