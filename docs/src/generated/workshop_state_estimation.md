```@meta
EditURL = "../../lit/workshop_state_estimation.jl"
```

# State estimation from noisy measurements

> **Level: Advanced to Expert**, companion of the advanced tour's state-estimation chapter; the observability section is the advanced deep dive.

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_state_estimation.ipynb)

A power flow computes the network state from an exact specification. A
real control room has the opposite problem: it receives many **redundant,
noisy measurements** (voltage magnitudes, injections, branch flows) and
must reconstruct the most likely state from them. That is *state
estimation*: a weighted-least-squares (WLS) fit of the bus voltages to
the measurement set, where each measurement counts with the inverse of
its variance.

In this notebook you build a 7-bus network with
[Sparlectra.jl](https://github.com/Welthulk/Sparlectra.jl), solve a
reference power flow, derive a noisy synthetic measurement set from it,
check observability, and run the estimator: the closed loop that lets
you judge estimation quality against a known truth.

The network is the workshop ring (B1 carries the grid connection, the
diagonals are the cross-ties B2-B5 and B3-B6):

```text
 (slack)
   B1 ---- B2 ---- B3 ---- B4
   |         \    /         |
   |          \  /          |
   |           \/           |
   |           /\           |
   B7 ---- B6 ---- B5 ------+
```

> **Note:** On Google Colab the install cell takes a few minutes on a
> fresh session (package download and precompilation). Colab's Julia
> version may change over time; this notebook targets Julia ≥ 1.12.

## Load the packages

`Random` (standard library) provides the seeded generator that makes the
synthetic measurement noise reproducible.

````@example workshop_state_estimation
using Sparlectra
using Random
````

## Warm-up

Julia compiles each function on first use. This cell warms the two paths
the notebook exercises, the power-flow solver (which produces the
reference state) and the WLS estimator, on a tiny throwaway network, so
the real study runs at full speed.

````@example workshop_state_estimation
wnet = Net(name = "warmup", baseMVA = 100.0)
addBus!(net = wnet, busName = "A", vn_kV = 110.0)
addBus!(net = wnet, busName = "B", vn_kV = 110.0)
addProsumer!(net = wnet, busName = "A", type = "EXTERNALNETWORKINJECTION", referencePri = "A", vm_pu = 1.0, va_deg = 0.0)
addProsumer!(net = wnet, busName = "B", type = "ENERGYCONSUMER", p = 10.0, q = 3.0)
addPIModelACLine!(net = wnet, fromBus = "A", toBus = "B", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
t_pf = @elapsed runpf!(wnet, 10, 1e-8, 0)
setMeasurementsFromPF!(wnet; includeVm = true, includePinj = true, includeQinj = true, includePflow = true, includeQflow = true, noise = false)
t_se = @elapsed runse!(wnet; maxIte = 8, tol = 1e-6, flatstart = true, jacEps = 1e-6, updateNet = false)
println("warm: power flow ", round(t_pf; digits = 2), " s, estimator ", round(t_se; digits = 2), " s (first calls compile)")
````

## Build the study network

Seven 110 kV buses in a ring with two cross-connections, enough meshing
that branch-flow measurements carry real information. An external network
injection at `B1` is the slack, a generator feeds at `B3`, the other
buses carry loads.

````@example workshop_state_estimation
net = Net(name = "workshop_se_7bus", baseMVA = 100.0)

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
````

## Solve the reference power flow

This solved state plays the role of the (in reality unknown) *true*
system state: the measurements are derived from it, and the estimator
never sees it directly.

````@example workshop_state_estimation
ite_pf, status_pf = runpf!(net, 40, 1e-10, 0)
status_pf == 0 || error("Power flow did not converge")
println("reference power flow converged in $ite_pf iterations")
````

## Derive a noisy measurement set

`measurementStdDevs` defines one standard deviation per measurement kind
(pu for voltages, MW/MVar for powers); `setMeasurementsFromPF!` then
reads the solved state and creates the measurements, with Gaussian noise
of exactly those standard deviations, seeded for reproducibility. This is
the standard trick for exercising an estimator: the truth is known, so
the residuals are meaningful.

````@example workshop_state_estimation
std = measurementStdDevs(vm = 1e-3, pinj = 1.0, qinj = 1.0, pflow = 0.7, qflow = 0.7)
setMeasurementsFromPF!(
  net;
  includeVm = true,
  includePinj = true,
  includeQinj = true,
  includePflow = true,
  includeQflow = true,
  noise = true,
  stddev = std,
  rng = MersenneTwister(42),
)
println(length(net.measurements), " measurements created")
````

## Check observability

Estimation only works where the measurement set actually determines the
state. The global check compares measurement count against state count
and probes the numerical rank of the measurement Jacobian.

````@example workshop_state_estimation
gobs = evaluate_global_observability(net; flatstart = true, jacEps = 1e-6)
println("global observability quality: ", gobs.quality)
println("measurements: ", gobs.n_measurements, ", states: ", gobs.n_states)
````

## Run the estimator

`runse!` iterates the WLS normal equations from a flat start. With
`updateNet = true` the estimated voltages are written back into the
network, so the usual result printers show the *estimated* state. The
objective $J$ is the weighted sum of squared residuals; for healthy
Gaussian noise it should land near the number of redundant measurements.

````@example workshop_state_estimation
se = runse!(
  net;
  maxIte = 12,
  tol = 1e-6,
  flatstart = true,
  jacEps = 1e-6,
  updateNet = true,
)

println("SE converged: ", se.converged, " in ", se.iterations, " iterations")
println("objective J:  ", round(se.objectiveJ; digits = 2))
````

## Inspect the estimated state

`SEResult.voltages` holds the estimated complex bus voltage per bus
index. Because the measurements carried only mild noise, the estimate
reproduces the reference power flow closely (magnitudes to a few
1e-3 pu), and the slack reference `B1` comes back at 1.02 pu / 0°. The
degrees of freedom (`dof`) and the 3σ check on $J$ summarize whether the
residuals are consistent with the declared measurement accuracies.

````@example workshop_state_estimation
println("dof (redundant measurements): ", se.dof)
println("J within 3σ band:             ", se.jWithin3Sigma)
println()
for (name, idx) in sort(collect(net.busDict); by = last)
  v = se.voltages[idx]
  println(rpad(name, 4), "  Vm = ", round(abs(v); digits = 4), " pu   Va = ", round(rad2deg(angle(v)); digits = 3), "°")
end
````

## Building measurement sets manually

Real measurement sets are not derived from a solved power flow; they
arrive from SCADA one by one. The `add*Measurement!` helpers resolve bus
and branch references for you, so assembling a set works just like
building the network. A deliberately sparse set like this one is a good
way to explore *when observability breaks down*:

````@example workshop_state_estimation
empty!(net.measurements)

addVmMeasurement!(net; busName = "B1", value = 1.02, sigma = 0.002)
addPinjMeasurement!(net; busName = "B2", value = -35.0, sigma = 1.0)
addQinjMeasurement!(net; busName = "B2", value = -10.0, sigma = 1.0)
addPflowMeasurement!(net; fromBus = "B1", toBus = "B2", value = 22.0, sigma = 0.8, direction = :from)
addQflowMeasurement!(net; branchNr = 1, value = 7.0, sigma = 0.8, direction = :to)

obs = evaluate_global_observability(net; flatstart = true, jacEps = 1e-6)
println("sparse manual set quality: ", obs.quality)
println("measurements: ", obs.n_measurements, ", states: ", obs.n_states)
````

Five measurements against 13 states: the check reports the shortfall
instead of letting the estimator run into a rank-deficient system.

## Observability deep dive: where measurements must sit (advanced)

The checks above said WHETHER the state is observable. This section is
about the WHY and the WHERE: what the numbers in the observability
result mean, which placements make a network observable, and how the
local check points at the corner that a missing meter leaves dark.

The study network again, as a picture (ring plus two chords):

```text
     B1 ──── B2 ──── B3
    /        │        │ ╲
  B7         │        │  B4
    │        │        │  │
    │        └─ B5 ───┼──┘
    │           │     │      chords: B2─B5 and B3─B6
    └─── B6 ────┴─────┘      ring:   B1-B2-B3-B4-B5-B6-B7-B1
```

The STATE the estimator solves for has 13 entries: one voltage angle per
non-slack bus (columns 1..6 for B2..B7; the slack B1 provides the angle
reference) and one voltage magnitude per bus (columns 7..13 for
B1..B7). Every measurement is one ROW of the measurement Jacobian $H$;
the state is observable exactly when those rows span all 13 columns,
i.e. $\mathrm{rank}(H) = 13$. `evaluate_global_observability` answers
that twice: STRUCTURALLY (a bipartite matching on the sparsity pattern,
`structural_matching` of 13 means every state column can be paired with
its own measurement row) and NUMERICALLY (the actual rank; a placement
can be structurally fine but numerically degenerate). The `quality`
verdict compresses it: `:good` needs observability AND redundancy
without critical single points; `:critical` is observable but one lost
measurement would blind it; `:not_observable` speaks for itself.

### A minimal placement, built by hand

The classical minimal recipe: a P/Q flow pair on every branch of a
SPANNING TREE (six branches for seven buses), plus one voltage-magnitude
anchor. The flow pairs fix all relative angles and magnitude ratios
along the tree, the anchor pins the absolute magnitude level, the slack
pins the absolute angle. That is 6 x 2 + 1 = 13 rows for 13 states,
observability with ZERO redundancy:

````@example workshop_state_estimation
empty!(net.measurements)
# re-solve the reference so the branch flows are available as true values
ite_ref, status_ref = runpf!(net, 40, 1e-10, 0)
status_ref == 0 || error("reference power flow did not converge")
calcNetLosses!(net)

tree = [("B1", "B2"), ("B2", "B3"), ("B3", "B4"), ("B4", "B5"), ("B5", "B6"), ("B6", "B7")]
for (f, t) in tree
  addPflowMeasurement!(net; fromBus = f, toBus = t, value = get_branch_p_from_to_mw(net, f, t), sigma = 0.8, direction = :from)
  addQflowMeasurement!(net; fromBus = f, toBus = t, value = get_branch_q_from_to_mvar(net, f, t), sigma = 0.8, direction = :from)
end
addVmMeasurement!(net; busName = "B1", value = net.nodeVec[net.busDict["B1"]]._vm_pu, sigma = 0.002)

gmin = evaluate_global_observability(net; flatstart = true, jacEps = 1e-6)
println("minimal tree set: quality = ", gmin.quality, " (", gmin.n_measurements, " rows, ", gmin.n_states, " states)")
println("  structurally observable: ", gmin.structural_observable, ", numerically observable: ", gmin.numerical_observable, " (rank ", gmin.numerical_rank, ")")
println("  redundancy dof = ", gmin.dof, ", critical measurements: ", length(gmin.numerical_critical_measurement_indices), " of ", gmin.n_measurements)
````

Reading aid: observable, but `quality = :critical` and EVERY row is on
the critical list: with zero redundancy, losing any single meter blinds
some part of the state. Real placements add the ring-closing and chord
flows (or injections) precisely to buy redundancy; the `dof` count is
what the bad-data machinery later feeds on.

### Breaking a corner, and finding it with the local check

Drop the B6-B7 flow pair: the spur toward B7 loses its only meters.
Globally the verdict flips to `:not_observable`; LOCALLY the check can
say which states went dark. `evaluate_local_observability(net, cols)`
restricts $H$ to selected state columns; for B7 those are angle column 6
and magnitude column 13:

````@example workshop_state_estimation
empty!(net.measurements)
for (f, t) in tree[1:5]   ## tree WITHOUT B6-B7
  addPflowMeasurement!(net; fromBus = f, toBus = t, value = get_branch_p_from_to_mw(net, f, t), sigma = 0.8, direction = :from)
  addQflowMeasurement!(net; fromBus = f, toBus = t, value = get_branch_q_from_to_mvar(net, f, t), sigma = 0.8, direction = :from)
end
addVmMeasurement!(net; busName = "B1", value = net.nodeVec[net.busDict["B1"]]._vm_pu, sigma = 0.002)

gbroken = evaluate_global_observability(net; flatstart = true, jacEps = 1e-6)
println("without B6-B7: quality = ", gbroken.quality, " (rank ", gbroken.numerical_rank, " of ", gbroken.n_states, ")")

# local check on B2 (angle col 1, magnitude col 8): still fully covered
lb2 = evaluate_local_observability(net, [1, 8]; flatstart = true, jacEps = 1e-6)
println("local B2: numerically observable = ", lb2.numerical_observable)
# local check on B7 (angle col 6, magnitude col 13): NO measurement even
# touches these states anymore, which the check reports as an error
try
  evaluate_local_observability(net, [6, 13]; flatstart = true, jacEps = 1e-6)
catch err
  println("local B7: ", sprint(showerror, err))
end
````

### Repairing with an injection

A meter does not have to sit ON the dark bus pair. An INJECTION
measurement at B7 couples B7 to every neighbor (B6 and B1), because the
injected power is the sum over its incident branches; its Jacobian row
touches all their states. P and Q injection at B7 restore rank 13:

````@example workshop_state_estimation
p7 = net.nodeVec[net.busDict["B7"]]._pƩGen === nothing ? 0.0 : net.nodeVec[net.busDict["B7"]]._pƩGen
addPinjMeasurement!(net; busName = "B7", value = p7 - 20.0, sigma = 1.0)   ## net injection: 20 MW load
addQinjMeasurement!(net; busName = "B7", value = -6.0, sigma = 1.0)
grepaired = evaluate_global_observability(net; flatstart = true, jacEps = 1e-6)
println("with P/Q injection at B7: quality = ", grepaired.quality, " (rank ", grepaired.numerical_rank, " of ", grepaired.n_states, ")")
````

Placement rules, condensed: flow pairs see the two ends of their branch,
injections see the bus AND all its neighbors, `Vm` anchors the magnitude
level (at least one required), and the slack anchors the angles. A bus
is dark exactly when no row of any of these reaches its columns.

### PMU phasors and the reference-angle offset

A PMU measures the voltage PHASOR: magnitude plus absolute angle,
GPS-synchronized. In the estimator that is a tightly weighted `VmMeas`
and a `VaMeas` in degrees (`addPmuPhasorMeasurement!` adds the pair).
One subtlety makes PMU angles interesting: the PMU time base rarely
coincides with the slack reference, so all PMU angles share a common
unknown OFFSET. With `state_estimation.pmu_ref_offset = auto` (the
default) the estimator adds that offset as an extra state and solves for
it; watch `n_states` grow from 13 to 14. We simulate two PMUs whose time
base is shifted by exactly 2 degrees:

````@example workshop_state_estimation
pmu_shift_deg = 2.0
for bus in ("B4", "B6")
  idx = net.busDict[bus]
  addPmuPhasorMeasurement!(net; busName = bus, vm_pu = net.nodeVec[idx]._vm_pu, va_deg = net.nodeVec[idx]._va_deg + pmu_shift_deg, sigmaVm = 0.002, sigmaVa = 0.02)
end
gpmu = evaluate_global_observability(net; flatstart = true, jacEps = 1e-6)
println("with 2 PMUs: quality = ", gpmu.quality, " (", gpmu.n_measurements, " rows, ", gpmu.n_states, " states, offset state included)")

se_pmu = runse!(net; maxIte = 15, tol = 1e-6, flatstart = true, jacEps = 1e-6, updateNet = false)
println("SE converged: ", se_pmu.converged, ", estimated PMU reference offset: ", round(se_pmu.vaRefOffsetDeg; digits = 3), "° (true shift ", pmu_shift_deg, "°)")
````

Reading aid: the estimator recovers the 2° time-base shift as the extra
state instead of bending the bus angles toward it, and the estimated
network state stays anchored to the slack reference. Without the offset
state (`pmu_ref_offset = off`) the shifted PMU angles would fight the
flow measurements and inflate $J$. Theory and the measurement model:
[State Estimation](https://welthulk.github.io/Sparlectra.jl/state_estimation/).

## Where to go next

- New to Sparlectra? The
  [workshop tour](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_tour.ipynb)
  builds a network from scratch step by step, directly in Colab.
- [Slack types and short-circuit currents](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_slack_short_circuit.ipynb):
  the grid-connection notebook: slack representations and IEC 60909-0
  fault currents.
- [State Estimation](https://welthulk.github.io/Sparlectra.jl/state_estimation/):
  the full documentation: observability analysis, PMU angle measurements,
  bad-data handling, and the measurement model.
- [State-Estimation Configuration](https://welthulk.github.io/Sparlectra.jl/state_estimation_configuration/):
  every `state_estimation.*` configuration key.

