```@meta
EditURL = "../../lit/workshop_tour.jl"
```

# The Sparlectra workshop tour

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_tour.ipynb)

All workshop examples in **one session**: install
[Sparlectra.jl](https://github.com/Welthulk/Sparlectra.jl) once, warm the
compiler up once, and then walk through the chapters without waiting
again. The focused single-topic notebooks cover the same ground with more
narrative; this tour is the fast lane.

After the warm-up (compilation happens there, everything after is fast)
the chapters are:

1. A first power flow
2. Slack types and short-circuit currents
3. Transformer tap control (OLTC)
4. Voltage-dependent reactive power, Q(U)
5. Remote voltage control by a machine
6. State estimation

> **Note:** On Google Colab the install cell takes a few minutes on a
> fresh session (package download and precompilation). Colab's Julia
> version may change over time; this notebook targets Julia ≥ 1.12.

## Warm-up

Julia compiles each function on first use. The tiny network below triggers
that compilation for the solver path once, so every later chapter runs at
full speed. The two timings make the effect visible.

````@example workshop_tour
using Sparlectra

wnet = Net(name = "warmup", baseMVA = 100.0)
addBus!(net = wnet, busName = "A", vn_kV = 110.0)
addBus!(net = wnet, busName = "B", vn_kV = 110.0)
addProsumer!(net = wnet, busName = "A", type = "EXTERNALNETWORKINJECTION", referencePri = "A", vm_pu = 1.0, va_deg = 0.0)
addProsumer!(net = wnet, busName = "B", type = "ENERGYCONSUMER", p = 10.0, q = 3.0)
addPIModelACLine!(net = wnet, fromBus = "A", toBus = "B", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)

t_first = @elapsed runpf!(wnet, 10, 1e-8, 0)
t_second = @elapsed runpf!(wnet, 10, 1e-8, 0)
println("first solve (includes compilation): ", round(t_first; digits = 2), " s")
println("second solve (compiled):            ", round(t_second * 1000; digits = 2), " ms")
````

A solve helper used by all chapters (25 iterations maximum, tolerance
$10^{-8}$):

````@example workshop_tour
function solve!(net; kwargs...)
  etime = @elapsed begin
    ite, erg = runpf!(net, 25, 1e-8, 0; kwargs...)
  end
  erg == 0 || error("Power flow did not converge (status = $erg)")
  calcNetLosses!(net)
  return etime, ite
end
````

## Chapter 1: a first power flow

Seven 110 kV buses in a ring with two cross-connections, an external
network injection as slack at `B1`, a generator at `B3`, loads elsewhere.
The builder is a function because chapter 7 reuses the same network.
The full guided version of this chapter is the
[introduction notebook](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_intro.ipynb).

````@example workshop_tour
function build_ring7(name::String)
  net = Net(name = name, baseMVA = 100.0)
  addBus!(net = net, busName = "B1", vn_kV = 110.0, vm_pu = 1.02, va_deg = 0.0)
  for i in 2:7
    addBus!(net = net, busName = "B$(i)", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
  end
  addPIModelACLine!(net = net, fromBus = "B1", toBus = "B2", r_pu = 0.010, x_pu = 0.080, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "B2", toBus = "B3", r_pu = 0.011, x_pu = 0.085, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "B3", toBus = "B4", r_pu = 0.012, x_pu = 0.090, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "B4", toBus = "B5", r_pu = 0.010, x_pu = 0.080, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "B5", toBus = "B6", r_pu = 0.011, x_pu = 0.085, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "B6", toBus = "B7", r_pu = 0.012, x_pu = 0.090, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "B7", toBus = "B1", r_pu = 0.010, x_pu = 0.080, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "B2", toBus = "B5", r_pu = 0.009, x_pu = 0.070, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "B3", toBus = "B6", r_pu = 0.009, x_pu = 0.070, b_pu = 0.0, status = 1)
  addProsumer!(net = net, busName = "B1", type = "EXTERNALNETWORKINJECTION", referencePri = "B1", vm_pu = 1.02, va_deg = 0.0)
  addProsumer!(net = net, busName = "B3", type = "GENERATOR", p = 60.0, q = 10.0)
  addProsumer!(net = net, busName = "B2", type = "LOAD", p = 35.0, q = 10.0)
  addProsumer!(net = net, busName = "B4", type = "LOAD", p = 45.0, q = 15.0)
  addProsumer!(net = net, busName = "B5", type = "LOAD", p = 25.0, q = 8.0)
  addProsumer!(net = net, busName = "B6", type = "LOAD", p = 30.0, q = 10.0)
  addProsumer!(net = net, busName = "B7", type = "LOAD", p = 20.0, q = 6.0)
  ok, msg = validate!(net = net)
  ok || error("Network validation failed: $msg")
  return net
end

net1 = build_ring7("tour_first_pf")
etime, ite = solve!(net1)
printACPFlowResults(net1, etime, ite, 1e-8)
````

Reading aid: the slack at `B1` covers the difference between 155 MW of
load, 60 MW of scheduled generation, and the network losses. All bus
voltages stay near 1.0 pu.

## Chapter 2: slack types and short-circuit currents

One grid connection, modeled three ways, plus an IEC 60909-0 fault-current
sweep from the declared feeder data. The detailed walk-through with full
reading aids is the
[slack-types notebook](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_slack_short_circuit.ipynb).

````@example workshop_tour
function build_grid(mode::Symbol)
  net = Net(name = "tour_eg8_$(mode)", baseMVA = 100.0)
  for b in ("B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8")
    addBus!(net = net, busName = b, vn_kV = 110.0)
  end
  addPIModelACLine!(net = net, fromBus = "B1", toBus = "B2", r_pu = 0.010, x_pu = 0.060, b_pu = 0.02, status = 1)
  addPIModelACLine!(net = net, fromBus = "B2", toBus = "B3", r_pu = 0.015, x_pu = 0.080, b_pu = 0.02, status = 1)
  addPIModelACLine!(net = net, fromBus = "B3", toBus = "B4", r_pu = 0.020, x_pu = 0.090, b_pu = 0.02, status = 1)
  addPIModelACLine!(net = net, fromBus = "B4", toBus = "B5", r_pu = 0.012, x_pu = 0.070, b_pu = 0.02, status = 1)
  addPIModelACLine!(net = net, fromBus = "B5", toBus = "B6", r_pu = 0.015, x_pu = 0.075, b_pu = 0.02, status = 1)
  addPIModelACLine!(net = net, fromBus = "B6", toBus = "B7", r_pu = 0.018, x_pu = 0.085, b_pu = 0.02, status = 1)
  addPIModelACLine!(net = net, fromBus = "B7", toBus = "B8", r_pu = 0.010, x_pu = 0.055, b_pu = 0.02, status = 1)
  addPIModelACLine!(net = net, fromBus = "B8", toBus = "B1", r_pu = 0.011, x_pu = 0.065, b_pu = 0.02, status = 1)
  addPIModelACLine!(net = net, fromBus = "B2", toBus = "B7", r_pu = 0.020, x_pu = 0.100, b_pu = 0.02, status = 1)
  addPIModelACLine!(net = net, fromBus = "B3", toBus = "B6", r_pu = 0.022, x_pu = 0.110, b_pu = 0.02, status = 1)
  addProsumer!(net = net, busName = "B3", type = "GENERATOR", p = 60.0, vm_pu = 1.01, qMin = -60.0, qMax = 60.0)
  addProsumer!(net = net, busName = "B6", type = "GENERATOR", p = 40.0, vm_pu = 1.00, qMin = -40.0, qMax = 40.0)
  addProsumer!(net = net, busName = "B2", type = "ENERGYCONSUMER", p = 45.0, q = 12.0)
  addProsumer!(net = net, busName = "B4", type = "ENERGYCONSUMER", p = 50.0, q = 15.0)
  addProsumer!(net = net, busName = "B7", type = "ENERGYCONSUMER", p = 40.0, q = 10.0)
  addProsumer!(net = net, busName = "B8", type = "ENERGYCONSUMER", p = 25.0, q = 8.0)
  addExternalGrid!(net = net, busName = "B1", vm_pu = 1.02, sk_max_MVA = 2000.0, sk_min_MVA = 1500.0, rx_max = 0.1, internal_impedance = (mode === :source))
  ok, msg = validate!(net = net)
  ok || error("Network validation failed: $msg")
  return net
end
````

**Ideal slack**: `B1` is pinned at exactly 1.02 pu / 0° and absorbs the
whole imbalance (the `SLACK` row).

````@example workshop_tour
net_slack = build_grid(:slack)
etime, ite = solve!(net_slack)
printACPFlowResults(net_slack, etime, ite, 1e-8)
````

**Non-ideal source**: with `internal_impedance = true` the setpoint moves
to the hidden internal bus (last row, type `SOURCE`); the terminal `B1`
in the first row droops below 1.02 pu.

````@example workshop_tour
net_source = build_grid(:source)
etime, ite = solve!(net_source)
printACPFlowResults(net_source, etime, ite, 1e-8)
````

**Distributed slack**: the generators pick up the imbalance according to
their scheduled output (0.6/0.4); the slack row keeps only the reactive
balance.

````@example workshop_tour
net_dist = build_grid(:slack)
etime, ite = solve!(net_dist; distributed_slack_enabled = true, distributed_slack_p_mode = :pg_weighted)
printACPFlowResults(net_dist, etime, ite, 1e-8)
````

The three connection models side by side. The losses differ because the
flow pattern differs; a negative Q loss means the line charging produces
more reactive power than the flows consume.

````@example workshop_tour
println(rpad("scenario", 20), lpad("Vm(B1) pu", 11), lpad("P loss MW", 11), lpad("Q loss MVAr", 13), "   balanced by")
for (label, net, by) in (
  ("ideal slack", net_slack, "slack bus B1"),
  ("non-ideal source", net_source, "hidden source bus"),
  ("distributed slack", net_dist, "B3 (α=0.6) + B6 (α=0.4)"),
)
  pl, ql = getTotalLosses(net = net)
  println(rpad(label, 20), lpad(string(round(get_bus_vm_pu(net, "B1"); digits = 4)), 11), lpad(string(round(pl; digits = 3)), 11), lpad(string(round(ql; digits = 3)), 13), "   ", by)
end
````

**Short circuit**: the feeder data declared in `addExternalGrid!` feeds
`runShortCircuit!` directly. $I_k''$ is largest at the connection bus and
decays with electrical distance.

````@example workshop_tour
printShortCircuitResult(runShortCircuit!(net_slack; case = :max))
printShortCircuitResult(runShortCircuit!(net_slack; case = :min))
````

## Chapter 3: transformer tap control (OLTC)

A transformer with a ratio tap changer holds the voltage at a remote load
bus. The outer control loop moves the discrete tap until the target is
inside the deadband; the power flow itself stays untouched. Details:
[Control Framework](https://welthulk.github.io/Sparlectra.jl/control_framework/).

````@example workshop_tour
function build_oltc()
  net = Net(name = "tour_oltc", baseMVA = 100.0)
  addBus!(net = net, busName = "Slack", vn_kV = 110.0)
  addBus!(net = net, busName = "MV", vn_kV = 110.0)
  addBus!(net = net, busName = "Load", vn_kV = 110.0)
  addProsumer!(net = net, busName = "Slack", type = "EXTERNALNETWORKINJECTION", referencePri = "Slack", vm_pu = 1.0, va_deg = 0.0)
  addProsumer!(net = net, busName = "Load", type = "ENERGYCONSUMER", p = 60.0, q = 20.0)
  addPIModelTrafo!(net = net, fromBus = "Slack", toBus = "MV", r_pu = 0.004, x_pu = 0.06, b_pu = 0.0, ratio = 1.0, shift_deg = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "MV", toBus = "Load", r_pu = 0.02, x_pu = 0.10, b_pu = 0.01, status = 1)
  # enable the ratio tap on the transformer branch and give it a name the
  # controller can address
  t = getNetBranch(net = net, fromBus = "Slack", toBus = "MV")
  t.comp.cName = "T1"
  t.has_ratio_tap = true
  t.tap_min = 0.90
  t.tap_max = 1.10
  t.tap_step = 0.0125
  ok, msg = validate!(net = net)
  ok || error("Network validation failed: $msg")
  return net
end

net_oltc = build_oltc()
run_sparlectra(net = net_oltc)
println("uncontrolled: Vm(Load) = ", round(get_bus_vm_pu(net_oltc, "Load"); digits = 4), " pu")
````

Now attach the controller (voltage mode, discrete steps) and rerun.
`run_sparlectra` executes the outer control loop automatically when
controllers are present.

````@example workshop_tour
addTapController!(
  net_oltc;
  trafo = "T1",
  mode = :voltage,
  target_bus = "Load",
  target_vm_pu = 1.0,
  control_ratio = true,
  control_phase = false,
  is_discrete = true,
  deadband_vm_pu = 5e-3,
)
run_sparlectra(net = net_oltc)
println("controlled:   Vm(Load) = ", round(get_bus_vm_pu(net_oltc, "Load"); digits = 4), " pu")
printTapControllerSummary(stdout, net_oltc)
````

Reading aid: the summary shows the chosen tap position and the achieved
voltage. With a discrete 0.0125 step the controller stops as soon as the
target is inside the deadband, not at the exact setpoint.

## Chapter 4: voltage-dependent reactive power, Q(U)

A machine can follow a Q(U) droop characteristic: absorb reactive power
when its voltage is high, inject when it is low. Unlike the outer-loop
controllers above, Q(U) is solved **inside** Newton-Raphson. Details:
[Voltage Dependent Control](https://welthulk.github.io/Sparlectra.jl/voltage_dependent_control/).

````@example workshop_tour
function build_qu(p_load::Float64, q_load::Float64)
  net = Net(name = "tour_qu", baseMVA = 100.0)
  addBus!(net = net, busName = "B1", vn_kV = 110.0)
  addBus!(net = net, busName = "B2", vn_kV = 110.0)
  addBus!(net = net, busName = "B3", vn_kV = 110.0)
  addProsumer!(net = net, busName = "B1", type = "EXTERNALNETWORKINJECTION", referencePri = "B1", vm_pu = 1.02, va_deg = 0.0)
  addPIModelACLine!(net = net, fromBus = "B1", toBus = "B2", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "B2", toBus = "B3", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
  qu = QUController(
    make_characteristic(
      [(104.5, 30.0), (107.0, 20.0), (110.0, 0.0), (112.0, -10.0), (115.5, -20.0)];
      voltage_unit = :kV,
      value_unit = :MVAr,
      vn_kV = 110.0,
      sbase_MVA = 100.0,
      interpolation = :polynomial,
    );
    qmin_MVAr = -50.0,
    qmax_MVAr = 50.0,
    sbase_MVA = 100.0,
  )
  addProsumer!(net = net, busName = "B2", type = "SYNCHRONOUSMACHINE", p = 10.0, q = 0.0, qu_controller = qu)
  addProsumer!(net = net, busName = "B3", type = "ENERGYCONSUMER", p = p_load, q = q_load)
  ok, msg = validate!(net = net)
  ok || error("Network validation failed: $msg")
  return net
end
````

Light load: the machine bus sits above 110 kV, so the characteristic asks
the machine to **absorb** reactive power (negative Q).

````@example workshop_tour
net_qu_light = build_qu(5.0, 1.0)
etime, ite = solve!(net_qu_light)
printACPFlowResults(net_qu_light, etime, ite, 1e-8)
````

Heavy load: the voltage sags below 110 kV and the same characteristic
turns the machine into a reactive power **injector** (positive Q).

````@example workshop_tour
net_qu_heavy = build_qu(80.0, 25.0)
etime, ite = solve!(net_qu_heavy)
printACPFlowResults(net_qu_heavy, etime, ite, 1e-8)
````

Reading aid: compare the `Qg` value and the `Control` column (`Q(U)`) of
bus `B2` between the two tables; the sign flips with the voltage level,
exactly along the declared characteristic.

## Chapter 5: remote voltage control by a machine

A machine regulates the voltage at a **different** bus via its reactive
output, the counterpart of a CGMES `RegulatingControl` at a foreign
terminal. The outer loop drives the machine's Q until the remote target
is met, and parks honestly `at_limit` when the reactive range is too
small. Details:
[Remote Voltage Control](https://welthulk.github.io/Sparlectra.jl/remote_voltage_control/).

````@example workshop_tour
function build_rvc(qmin_mvar::Float64, qmax_mvar::Float64)
  net = Net(name = "tour_rvc", baseMVA = 100.0)
  addBus!(net = net, busName = "Slack", vn_kV = 110.0)
  addBus!(net = net, busName = "GenBus", vn_kV = 110.0)
  addBus!(net = net, busName = "Load", vn_kV = 110.0)
  addProsumer!(net = net, busName = "Slack", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.02, va_deg = 0.0, referencePri = "Slack")
  addProsumer!(net = net, busName = "GenBus", type = "SYNCHRONOUSMACHINE", p = 30.0, q = 0.0, qMin = qmin_mvar, qMax = qmax_mvar)
  addProsumer!(net = net, busName = "Load", type = "ENERGYCONSUMER", p = -70.0, q = -20.0)
  addPIModelACLine!(net = net, fromBus = "Slack", toBus = "GenBus", r_pu = 0.02, x_pu = 0.12, b_pu = 0.01, status = 1)
  addPIModelACLine!(net = net, fromBus = "GenBus", toBus = "Load", r_pu = 0.02, x_pu = 0.12, b_pu = 0.01, status = 1)
  return net
end

pf = PowerFlowConfig(max_iter = 30, tol = 1e-9)

net_rvc = build_rvc(-50.0, 50.0)
runpf!(net_rvc; config = pf, verbose = 0)
println("uncontrolled: Vm(Load) = ", round(get_bus_vm_pu(net_rvc, "Load"); digits = 4), " pu")

addMachineVoltageControl!(net_rvc; bus = "GenBus", target_bus = "Load", target_vm_pu = 1.05, deadband_vm_pu = 5e-4)
result = run_control!(net_rvc; controllers = collect_outer_controllers(net_rvc), pf_config = pf, control_config = ControlConfig(max_outer_iterations = 15), verbose = 0)
println("controlled:   Vm(Load) = ", round(get_bus_vm_pu(net_rvc, "Load"); digits = 4), " pu (target 1.05)")
println("loop: status = ", result.status, ", outer iterations = ", result.outer_iterations, ", pf solves = ", result.powerflow_solves)
printMachineControllerSummary(stdout, net_rvc)
````

The honest failure mode: cut the reactive range to ±2 MVAr and the same
target is out of reach. The controller parks at its limit and says so
instead of pretending convergence.

````@example workshop_tour
net_rvc2 = build_rvc(-2.0, 2.0)
runpf!(net_rvc2; config = pf, verbose = 0)
addMachineVoltageControl!(net_rvc2; bus = "GenBus", target_bus = "Load", target_vm_pu = 1.05, deadband_vm_pu = 5e-4)
result2 = run_control!(net_rvc2; controllers = collect_outer_controllers(net_rvc2), pf_config = pf, control_config = ControlConfig(max_outer_iterations = 15), verbose = 0)
println("limited:      Vm(Load) = ", round(get_bus_vm_pu(net_rvc2, "Load"); digits = 4), " pu (target 1.05)")
printMachineControllerSummary(stdout, net_rvc2)
````

## Chapter 6: state estimation

Close the loop: solve a reference power flow on the chapter-1 network,
derive a noisy synthetic measurement set from it, check observability,
and let the WLS estimator reconstruct the state. The full narrative is
the
[state-estimation notebook](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_state_estimation.ipynb).

````@example workshop_tour
using Random

net_se = build_ring7("tour_se")
ite_pf, status_pf = runpf!(net_se, 40, 1e-10, 0)
status_pf == 0 || error("Power flow did not converge")

std = measurementStdDevs(vm = 1e-3, pinj = 1.0, qinj = 1.0, pflow = 0.7, qflow = 0.7)
setMeasurementsFromPF!(
  net_se;
  includeVm = true,
  includePinj = true,
  includeQinj = true,
  includePflow = true,
  includeQflow = true,
  noise = true,
  stddev = std,
  rng = MersenneTwister(42),
)

gobs = evaluate_global_observability(net_se; flatstart = true, jacEps = 1e-6)
println("observability: ", gobs.quality, " (", gobs.n_measurements, " measurements, ", gobs.n_states, " states)")

se = runse!(net_se; maxIte = 12, tol = 1e-6, flatstart = true, jacEps = 1e-6, updateNet = true)
println("SE converged: ", se.converged, " in ", se.iterations, " iterations")
println("objective J:  ", round(se.objectiveJ; digits = 2), " (dof ", se.dof, ", within 3σ: ", se.jWithin3Sigma, ")")
for (name, idx) in sort(collect(net_se.busDict); by = last)
  v = se.voltages[idx]
  println(rpad(name, 4), "  Vm = ", round(abs(v); digits = 4), " pu   Va = ", round(rad2deg(angle(v)); digits = 3), "°")
end
````

Reading aid: with mild noise the estimate reproduces the reference state
to a few 1e-3 pu, and $J$ lands near the degrees of freedom, the
textbook health check for a WLS estimator.

## Where to go next

The focused notebooks with the full narrative:

- [Introduction](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_intro.ipynb)
- [Slack types and short circuit](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_slack_short_circuit.ipynb)
- [State estimation](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_state_estimation.ipynb)

And the documentation for the chapters that have no own notebook yet:

- [Control Framework](https://welthulk.github.io/Sparlectra.jl/control_framework/) and
  the [Workshop](https://welthulk.github.io/Sparlectra.jl/workshop/) tap-control section
- [Voltage Dependent Control](https://welthulk.github.io/Sparlectra.jl/voltage_dependent_control/)
- [Remote Voltage Control](https://welthulk.github.io/Sparlectra.jl/remote_voltage_control/)
- [Feature Matrix](https://welthulk.github.io/Sparlectra.jl/feature_matrix/)

