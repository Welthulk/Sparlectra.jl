```@meta
EditURL = "../../lit/workshop_state_estimation.jl"
```

# State estimation from noisy measurements

> **Level: Advanced to Expert**, companion of the advanced tour's state-estimation chapter; the observability section is the advanced deep dive.

> **Note:** This workshop was created with AI assistance and is reviewed
> and curated by the maintainer; it is not a fully machine-generated text.

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

## Warm-up

Julia compiles each function on first use. This cell loads the packages
the chapters share. `Random` (standard library) provides the seeded
generator that makes the synthetic measurement noise reproducible;
`LinearAlgebra` (standard library) contributes `nullspace` for the
observability deep dive; `Printf` (standard library) aligns the printed
values. The notebook needs no helper functions; the
chapters share the study network of Example 1 and build on each other's
measurement sets, and a code cell that uses a name from another cell
says so in a comment. Its last line, `warmup()`, runs every path the
chapters exercise once on tiny throwaway networks (the power-flow solver
that produces the reference state, the WLS estimator, the observability
checks (global, local, on a matrix), the hand-placed and zero-injection
measurements, the measurement Jacobian, the PMU phasors, the observable
islands, the estimate per island, the restoration and the placement),
so the real study runs at full speed; its code is in
`docs/lit/warmup/workshop_state_estimation.jl` of the package. How long
the warm-up takes depends on the machine; every later cell then runs without compile
pauses.

````@example workshop_state_estimation
using Sparlectra
using Random
using LinearAlgebra
using Printf

# compile every path the chapters use once; the warm-up code is in
# docs/lit/warmup/workshop_state_estimation.jl of the installed package
include(pkgdir(Sparlectra, "docs", "lit", "warmup", "workshop_state_estimation.jl"))
warmup()
````

## Build the study network

**Example 1: the closed estimation loop, from network build to estimated
state.** Seven 110 kV buses in a ring with two cross-connections, the
ring drawn in the introduction above, enough meshing
that branch-flow measurements carry real information. An external network
injection at `B1` is the slack, a generator feeds at `B3`, the other
buses carry loads. The sections up to the state inspection all belong to
this example: reference power flow, noisy measurement set, observability
check, estimation.

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
# uses: net (Example 1)
@isdefined(net) || error("Run Example 1 first: it sets up net.")
ite_pf, status_pf = runpf!(net, 40, 1e-10, 0)
status_pf == 0 || error("Power flow did not converge")
println("reference power flow (true state): ", ite_pf, " iterations")
# the true state, kept for the comparison with the estimate below
vm_true = [net.nodeVec[i]._vm_pu for i in eachindex(net.nodeVec)]
va_true = [net.nodeVec[i]._va_deg for i in eachindex(net.nodeVec)]
````

## Derive a noisy measurement set

`measurementStdDevs` defines one standard deviation per measurement kind
(pu for voltages, MW/MVar for powers); `setMeasurementsFromPF!` then
reads the solved state and creates the measurements, with Gaussian noise
of exactly those standard deviations, seeded for reproducibility. This is
the standard trick for exercising an estimator: the truth is known, so
the residuals are meaningful.

````@example workshop_state_estimation
# uses: net (Example 1, Solve the reference power flow)
@isdefined(net) || error("Run Example 1, the section \"Solve the reference power flow\" first: they set up net.")
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
println("measurement set from the true state: ", length(net.measurements), " rows, seeded noise")
````

## Check observability

Estimation only works where the measurement set actually determines the
state. Why that is a RANK question follows from what the estimator
minimizes. With $m$ measurements $z_i$, their models $h_i(x)$, and
accuracies $\sigma_i$, WLS minimizes the weighted squared residuals

$$J(x) = \sum_{i=1}^{m} \frac{(z_i - h_i(x))^2}{\sigma_i^2}$$

and every Newton step solves the normal equations

$$G\,\Delta x = H^{\top} W r, \qquad G = H^{\top} W H, \qquad W = \mathrm{diag}(1/\sigma_i^2)$$

with $H$ the measurement Jacobian. $G$ is invertible exactly when
$\mathrm{rank}(H) = n$, the number of states: observability is not an
extra condition on top of the estimator, it IS the solvability of these
equations. For this workshop network $n = 2\,n_{bus} - 1 = 13$ (one
magnitude per bus, one angle per non-slack bus), so the 13 appearing
below is a formula, not a coincidence. The global check compares $m$
against $n$ and probes the numerical rank of $H$.

````@example workshop_state_estimation
# uses: net (Example 1, Solve the reference power flow, Derive a noisy measurement set)
@isdefined(net) || error("Run Example 1, the section \"Solve the reference power flow\", the section \"Derive a noisy measurement set\" first: they set up net.")
gobs = with_state_estimation_config(flatstart = true, jac_eps = 1e-6) do
  evaluate_global_observability(net)
end
println("global observability of the set: quality = ", gobs.quality)
println("  measurement rows m = ", gobs.n_measurements, ", state variables n = ", gobs.n_states)
````

## Run the estimator

`runse!` iterates the WLS normal equations from a flat start. With
`state_estimation.update_net = true` the estimated voltages are written
back into the network, so the usual result printers show the *estimated*
state. The estimator takes its settings from the active configuration;
`with_state_estimation_config` installs other values for one run. The
objective $J$ is the weighted sum of squared residuals; for healthy
Gaussian noise it should land near the number of redundant measurements.

````@example workshop_state_estimation
# uses: net (Example 1, Solve the reference power flow, Derive a noisy measurement set)
@isdefined(net) || error("Run Example 1, the section \"Solve the reference power flow\", the section \"Derive a noisy measurement set\" first: they set up net.")
se = with_state_estimation_config(max_iter = 12, tol = 1e-6, flatstart = true, jac_eps = 1e-6, update_net = true) do
  runse!(net)
end

println("SE converged: ", se.converged, " in ", se.iterations, " iterations")
@printf("objective J of the estimate [-]: %.2f\n", se.objectiveJ)
````

## Inspect the estimated state

`SEResult.voltages` holds the estimated complex bus voltage per bus
index. Because the measurements carried only mild noise, the estimate
reproduces the reference power flow closely: magnitudes within 1e-3 pu,
angles within 0.02°. The slack `B1` keeps its angle of 0° (the angle
reference); its magnitude is estimated like every other one and lands
near the true 1.02 pu. The degrees of freedom (`dof`) and the 3σ check on $J$ summarize whether the
residuals are consistent with the declared measurement accuracies.

$J$ is more than a convergence number, it is a QUALITY measure: under
the assumed noise model, $J(\hat{x})$ is approximately chi-square
distributed with $\mathrm{dof} = m - n$ degrees of freedom, so its
expected value is about `dof`. A $J$ far ABOVE `dof` signals bad data or
a wrong model; far BELOW signals overfitted (too pessimistic) sigmas.
The ratio is the one number worth glancing at after every run:

````@example workshop_state_estimation
# uses: net (Example 1, Solve the reference power flow, Derive a noisy measurement set, Run the estimator); se (Run the estimator); va_true, vm_true (Solve the reference power flow)
@isdefined(net) || error("Run Example 1, the section \"Solve the reference power flow\", the section \"Derive a noisy measurement set\", the section \"Run the estimator\" first: they set up net.")
@isdefined(se) || error("Run the section \"Run the estimator\" first: it sets up se.")
(@isdefined(va_true) && @isdefined(vm_true)) || error("Run the section \"Solve the reference power flow\" first: it sets up va_true, vm_true.")
@printf("J/dof of the estimate [-]:       %.3f   (healthy noise: approximately 1)\n", se.objectiveJ / se.dof)

println("dof = m - n (redundant rows):    ", se.dof)
println("J within 3σ band:                ", se.jWithin3Sigma)
println()
println("Bus voltages, true (reference power flow) and estimated:")
println("  bus   Vm true [pu]   Vm est [pu]   Va true [deg]   Va est [deg]")
for (name, idx) in sort(collect(net.busDict); by = last)
  v = se.voltages[idx]
  @printf("  %-4s %12.4f %13.4f %15.3f %14.3f\n", name, vm_true[idx], abs(v), va_true[idx], rad2deg(angle(v)))
end
````

## Building measurement sets manually

**Example 2: a deliberately sparse manual set.** Real measurement sets
are not derived from a solved power flow; they
arrive from SCADA one by one. The `add*Measurement!` helpers resolve bus
and branch references for you, so assembling a set works just like
building the network. A deliberately sparse set like this one is a good
way to explore *when observability breaks down*:

````@example workshop_state_estimation
# uses: net (Example 1, Solve the reference power flow, Derive a noisy measurement set, Run the estimator)
@isdefined(net) || error("Run Example 1, the section \"Solve the reference power flow\", the section \"Derive a noisy measurement set\", the section \"Run the estimator\" first: they set up net.")
empty!(net.measurements)

addVmMeasurement!(net; busName = "B1", value = 1.02, sigma = 0.002)
addPinjMeasurement!(net; busName = "B2", value = -35.0, sigma = 1.0)
addQinjMeasurement!(net; busName = "B2", value = -10.0, sigma = 1.0)
addPflowMeasurement!(net; fromBus = "B1", toBus = "B2", value = 22.0, sigma = 0.8, direction = :from)
addQflowMeasurement!(net; branchNr = 1, value = 7.0, sigma = 0.8, direction = :to)

obs = with_state_estimation_config(flatstart = true, jac_eps = 1e-6) do
  evaluate_global_observability(net)
end
println("sparse manual set: quality = ", obs.quality)
println("  measurement rows m = ", obs.n_measurements, ", state variables n = ", obs.n_states)
````

Reading aid (Example 2): five measurements against 13 states: the check
reports the shortfall instead of letting the estimator run into a
rank-deficient system.

## Observability on paper: three small matrices

Before the network-sized deep dive, the whole theory on matrices small
enough to read: three numbered examples (Examples 3 to 5), each comment
names the matrix it talks about.

**Example 3 (`H_B`), the minimal case:** three measurements, three
states, the identity. Observable (structurally: every column finds its
own row in the matching; numerically: rank 3), but with $m = n$ every
single row is CRITICAL, losing any one loses a state:

````@example workshop_state_estimation
H_B = [
  1.0 0.0 0.0
  0.0 1.0 0.0
  0.0 0.0 1.0
]
obs_B = evaluate_observability_matrix(H_B)
println("H_B: numerically observable = ", obs_B.numerical_observable, ", dof = m - n = ", obs_B.dof)
println("  critical rows of H_B: ", obs_B.numerical_critical_measurement_indices)
````

**Example 4 (`H_A`), duplicated information:** the same identity plus
two extra rows, but BOTH duplicating the same information about
states 1 and 2. The lesson: redundancy is PER STATE, not global. The
set has dof = 2, yet row 3 stays critical, because it is still the only
row that sees state 3; the two spare rows protect the wrong place:

````@example workshop_state_estimation
H_A = [
  1.0 0.0 0.0
  0.0 1.0 0.0
  0.0 0.0 1.0
  1.0 1.0 0.0
  1.0 1.0 0.0
]
obs_A = evaluate_observability_matrix(H_A)
println("H_A: numerically observable = ", obs_A.numerical_observable, ", dof = m - n = ", obs_A.dof)
println("  critical rows of H_A: ", obs_A.numerical_critical_measurement_indices)
for i in axes(H_A, 1)
  println("  row ", i, " of H_A: numerical check ", numerical_row_redundant(H_A, i) ? "redundant" : "CRITICAL", ", structural check ", structural_row_redundant(H_A, i) ? "redundant" : "CRITICAL")
end
````

**Example 5 (`H_E`), the bridge to networks:** an incidence-like matrix
for a 4-bus chain, one row per measurement (the listing below says
which). In this
idealized form you can READ the measurement type off the numbers: a
flow between buses a and b sees only their DIFFERENCE, so its row is
the pair +1 in column a and -1 in column b; an injection at bus k
anchors one state, so its row has a single +1 in column k. Check it
row by row: rows 1, 2, 3, 5 each carry a +1/-1 pair (flows), row 4 is
the only single-entry row (the injection at bus 2).

```text
  row 1: flow 1-2    row 2: flow 2-3    row 3: flow 3-4
  row 4: injection at 2                 row 5: flow 1-3
```

One honest caveat before transferring this: the legibility is a
property of the TOY. In a real AC measurement Jacobian an injection row
touches its bus AND every neighbor, a flow row touches three or four
states, and the entries are sensitivities, not ones (the
`measurement_jacobian` output below, in Example 7, says "PflowMeas at
B1->B2 touches 3 states"). There the row DESCRIPTIONS carry the
what-and-where; the numbers only carry the coupling.

````@example workshop_state_estimation
H_E = [
  1.0 -1.0 0.0 0.0
  0.0 1.0 -1.0 0.0
  0.0 0.0 1.0 -1.0
  0.0 1.0 0.0 0.0
  1.0 0.0 -1.0 0.0
]
obs_E = evaluate_observability_matrix(H_E)
println("H_E: numerically observable = ", obs_E.numerical_observable, ", dof = m - n = ", obs_E.dof)
println("  critical rows of H_E: ", obs_E.numerical_critical_measurement_indices)
````

Reading aid (Example 5): row 3 (the only path to bus 4) is critical;
row 4 (the injection at bus 2) is critical too, it is the only row that
fixes an absolute level, the flow rows only see differences; rows 1, 2, 5
form a triangle of alternatives around buses 1..3 and are redundant.
Exactly this pattern, at network scale, is what the deep dive below
builds and breaks (Examples 7 and 8). The extended version of these
examples (with a toy spanning-tree game and tolerance experiments)
lives in `examples/state_estimation/h_matrix_observability_demo.jl`.

## Observability deep dive: where measurements must sit (advanced)

The checks above said WHETHER the state is observable. This section is
about the WHY and the WHERE: what the numbers in the observability
result mean, which placements make a network observable, and how the
local check points at the corner that a missing meter leaves dark.

First the two concepts in plain words:

**Global observability** asks: does the measurement set determine the
COMPLETE network state, every voltage magnitude and every angle, as one
consistent picture? If yes, the WLS estimator has a unique solution; if
no, whole regions of the state can drift without any measurement
noticing, and the normal equations are singular. It is a property of the
WHOLE set against the WHOLE state.

**Local observability** asks the same question for a chosen PART of the
state, typically the voltage and angle of one bus or one area: do the
measurements that touch these states pin them down? A network can be
globally unobservable while most of it is locally fine, and the local
check is the tool that finds the dark corner. Think of it as zooming the
rank question into a neighborhood.

Both are decided by the same object, the measurement Jacobian $H$: one
row per measurement, one column per state, entry = sensitivity of that
measurement to that state. Global observability is full column rank of
$H$; the sharpest way to say it uses the NULL SPACE:

$$H \nu = 0,\ \nu \neq 0 \quad\Longrightarrow\quad \text{state } i \text{ is unobservable when } \nu_i \neq 0$$

The null-space vectors are the dark regions: directions the state can
drift without any measurement changing. Their nonzero components are the
DARK STATES, which Sparlectra reports as `unobservable_state_columns` on
a not-observable global result. Dark states and observable islands are
two answers: the islands are the connected parts whose states are fixed
once one of their buses is the reference (Monticelli and Wu), found by a
separate analysis, see Example 13.

**Example 6: why a local verdict can mislead.** One warning before
using the LOCAL check: restricting $H$ to the
columns you zoomed into (with the rows that touch them) is a NECESSARY
test, not a sufficient one. The smallest counterexample: a single flow
measurement between buses 1 and 2,

$$H = \begin{pmatrix} 1 & -1 \end{pmatrix}$$

The submatrix for state 1 is just $(1)$, full rank, so the local test
says observable. But only the DIFFERENCE $x_1 - x_2$ is determined; the
null-space vector $\nu = (1, 1)$ has a nonzero first component, so
state 1 is not estimable. Both verdicts side by side:

````@example workshop_state_estimation
H_flow = [1.0 -1.0]
loc = evaluate_local_observability_matrix(H_flow, [1])
println("local submatrix test on state 1: observable = ", loc.numerical_observable, "  (misleading)")
println("null-space vector of H (normalized, states 1 and 2): ", vec(round.(nullspace(H_flow); digits = 4)))
let ns = vec(nullspace(H_flow))
end
glob_flow = evaluate_observability_matrix(H_flow)
println("global check, unobservable state columns: ", glob_flow.unobservable_state_columns, "  (the rigorous answer)")
````

Reading aid (Example 6): the local check remains useful, it localizes
candidate regions fast, but a positive local verdict must be confirmed
globally; the docstrings
of both local helpers carry the same warning.

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

**Example 7: a minimal placement, built by hand.** The classical
minimal recipe: a P/Q flow pair on every branch of a
SPANNING TREE (six branches for seven buses), plus one voltage-magnitude
anchor. The flow pairs fix all relative angles and magnitude ratios
along the tree, the anchor pins the absolute magnitude level, the slack
pins the absolute angle. The counting identity says why this is exactly
minimal:

$$m = 2\,(n_{bus} - 1) + 1 = 2\,n_{bus} - 1 = n$$

a spanning tree has $n_{bus} - 1$ branches, each contributing a P and a
Q row, plus the one anchor; $m = n$ means observability with ZERO
redundancy, every equation is needed:

````@example workshop_state_estimation
# uses: net (Example 1, Solve the reference power flow, Derive a noisy measurement set, Run the estimator, Example 2)
@isdefined(net) || error("Run Example 1, the section \"Solve the reference power flow\", the section \"Derive a noisy measurement set\", the section \"Run the estimator\", Example 2 first: they set up net.")
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
````

The same placement, drawn into the line diagram. Legend: `x` = voltage
meter (Vm), `o` = P/Q flow pair on the branch, `(o)` = P/Q injection
pair at the bus (used in the repair below, Example 9):

```text
    [x]B1 ──o── B2 ──o── B3
      /          │        │ ╲
    B7           │        │  o
     │           │        │   ╲
     o           └─ B5 ───┼─── B4
     │              o     │
     └─── B6 ───────┴─────┘
            (chords B2─B5, B3─B6 and ring closure B7─B1: UNMEASURED)

    x  voltage meter        o  P/Q flow pair on the branch
```

The picture already tells the placement story: the six `o` pairs walk a
spanning tree (every bus is reached), the single `x` anchors the
magnitude level, and the three unmeasured branches are exactly where
redundancy would come from.

````@example workshop_state_estimation
# uses: net (Example 1, Solve the reference power flow, Derive a noisy measurement set, Run the estimator, Example 2, Example 7)
@isdefined(net) || error("Run Example 1, the section \"Solve the reference power flow\", the section \"Derive a noisy measurement set\", the section \"Run the estimator\", Example 2, Example 7 first: they set up net.")
gmin = with_state_estimation_config(flatstart = true, jac_eps = 1e-6) do
  evaluate_global_observability(net)
end
println("minimal tree set: quality = ", gmin.quality, " (measurement rows m = ", gmin.n_measurements, ", state variables n = ", gmin.n_states, ")")
println("  structurally observable: ", gmin.structural_observable, ", numerically observable: ", gmin.numerical_observable, " (numerical rank of H = ", gmin.numerical_rank, ")")
println("  redundancy dof = m - n = ", gmin.dof, ", critical measurement rows: ", length(gmin.numerical_critical_measurement_indices), " of ", gmin.n_measurements)
````

Reading aid (Example 7): observable, but `quality = :critical` and EVERY row is on
the critical list: with zero redundancy, losing any single meter blinds
some part of the state. There is a beautiful residual-side view of the
same fact:

$$\mathrm{dof} = m - n, \qquad S = I - H\,G^{-1} H^{\top} W, \qquad \text{measurement } i \text{ critical} \iff S_{ii} = 0$$

$S$ maps measurements to their residuals; a critical measurement has
residual EXACTLY zero, whatever its value, so a gross error in it is
invisible to bad-data detection. That is why `dof` is what the bad-data
machinery feeds on: only redundant measurements leave residual traces.
One precision on Sparlectra's `dof` field: it is the COUNT difference
$m - n$. For an observable set that equals $m - \mathrm{rank}(H)$; for
an unobservable set it can be negative and then reads as a shortfall,
not a redundancy. Real placements add the ring-closing and chord flows
(or injections) precisely to buy this redundancy.

The matrix behind all of this is one call away: `measurement_jacobian`
returns $H$, here for the Example 7 placement, with described rows and
labeled state columns, ready for a
placement report (the state-estimation example suite writes exactly such
a matrix-plus-stability page on every run):

````@example workshop_state_estimation
# uses: net (Example 1, Solve the reference power flow, Derive a noisy measurement set, Run the estimator, Example 2, Example 7)
@isdefined(net) || error("Run Example 1, the section \"Solve the reference power flow\", the section \"Derive a noisy measurement set\", the section \"Run the estimator\", Example 2, Example 7 first: they set up net.")
mj = measurement_jacobian(net)
println("measurement Jacobian H: ", size(mj.H, 1), " measurement rows x ", size(mj.H, 2), " state columns; first columns: ", join(mj.cols[1:4], ", "), ", ...")
r1 = mj.rows[1]
println("row 1 of H: ", r1.type, " at ", r1.location, " touches ", count(j -> abs(mj.H[1, j]) > 1e-9, eachindex(mj.cols)), " states (nonzero entries)")
````

### Breaking a corner, and finding it with the local check

**Example 8: breaking a corner.** Drop the B6-B7 flow pair from the
Example 7 tree: the spur toward B7 loses its only meters.
Globally the verdict flips to `:not_observable`; LOCALLY the check can
say which states went dark. `evaluate_local_observability(net, cols)`
restricts $H$ to selected state columns; for B7 those are angle column 6
and magnitude column 13:

````@example workshop_state_estimation
# uses: net (Example 1, Solve the reference power flow, Derive a noisy measurement set, Run the estimator, Example 2, Example 7); tree (Example 7)
@isdefined(net) || error("Run Example 1, the section \"Solve the reference power flow\", the section \"Derive a noisy measurement set\", the section \"Run the estimator\", Example 2, Example 7 first: they set up net.")
@isdefined(tree) || error("Run Example 7 first: it sets up tree.")
empty!(net.measurements)
for (f, t) in tree[1:5]   ## tree WITHOUT B6-B7
  addPflowMeasurement!(net; fromBus = f, toBus = t, value = get_branch_p_from_to_mw(net, f, t), sigma = 0.8, direction = :from)
  addQflowMeasurement!(net; fromBus = f, toBus = t, value = get_branch_q_from_to_mvar(net, f, t), sigma = 0.8, direction = :from)
end
addVmMeasurement!(net; busName = "B1", value = net.nodeVec[net.busDict["B1"]]._vm_pu, sigma = 0.002)

gbroken = with_state_estimation_config(flatstart = true, jac_eps = 1e-6) do
  evaluate_global_observability(net)
end
println("without the B6-B7 flow pair: quality = ", gbroken.quality, " (numerical rank of H = ", gbroken.numerical_rank, " of n = ", gbroken.n_states, " states)")
# the rigorous per-state answer, straight from the null space: exactly
# B7's states are dark (angle column 6, magnitude column 13)
println("dark states (unobservable_state_columns): ", gbroken.unobservable_state_columns)

# local check on B2 (angle col 1, magnitude col 8): still fully covered
lb2 = with_state_estimation_config(flatstart = true, jac_eps = 1e-6) do
  evaluate_local_observability(net, [1, 8])
end
println("local B2: numerically observable = ", lb2.numerical_observable)
# local check on B7 (angle col 6, magnitude col 13): NO measurement even
# touches these states anymore, which the check reports as an error;
# this is the EASY case for the local test, the hard case (touched but
# still dark) is the Example 6 counterexample, which only the global
# null-space answer catches
try
  with_state_estimation_config(flatstart = true, jac_eps = 1e-6) do
    evaluate_local_observability(net, [6, 13])
  end
catch err
  println("local B7: ", sprint(showerror, err))
end
````

### Repairing with an injection

**Example 9: repairing with an injection.** A meter does not have to
sit ON the dark bus pair. Starting from the broken set of Example 8, an
INJECTION measurement at B7 couples B7 to every neighbor (B6 and B1),
because the
injected power is the sum over its incident branches; its Jacobian row
touches all their states. P and Q injection at B7 restore rank 13:

````@example workshop_state_estimation
# uses: net (Example 1, Solve the reference power flow, Derive a noisy measurement set, Run the estimator, Example 2, Example 7, Example 8)
@isdefined(net) || error("Run Example 1, the section \"Solve the reference power flow\", the section \"Derive a noisy measurement set\", the section \"Run the estimator\", Example 2, Example 7, Example 8 first: they set up net.")
# B7 carries exactly its 20 MW / 6 MVAr load, so the true net injection
# is -20 MW / -6 MVAr (consumption counts negative)
addPinjMeasurement!(net; busName = "B7", value = -20.0, sigma = 1.0)
addQinjMeasurement!(net; busName = "B7", value = -6.0, sigma = 1.0)
grepaired = with_state_estimation_config(flatstart = true, jac_eps = 1e-6) do
  evaluate_global_observability(net)
end
println("with P/Q injection at B7: quality = ", grepaired.quality, " (numerical rank of H = ", grepaired.numerical_rank, " of n = ", grepaired.n_states, " states)")
````

Placement rules, condensed: flow pairs see the two ends of their branch,
injections see the bus AND all its neighbors, `Vm` anchors the magnitude
level (at least one required), and the slack anchors the angles. A bus
is dark exactly when no row of any of these reaches its columns.

### Passive buses: zero-injection pseudo-measurements

A PASSIVE bus carries no load and no generation, so its injection is
exactly zero, and that knowledge is itself a measurement pair
$P = 0,\ Q = 0$, free of charge and free of noise. Sparlectra models it
as PSEUDO-measurements with a very small sigma
(`addZeroInjectionMeasurements!`), not as hard constraints. With
`state_estimation.zero_injection = constraint` the same rows become hard
equality constraints instead. The trade-off of the default is honest:
the smaller the sigma, the closer to a true constraint, but the weight
$1/\sigma^2$ explodes and degrades the conditioning of $G$. The ratio
between a ZIB row at $\sigma = 10^{-6}$ and a normal injection row at
$\sigma = 1$ is

$$\frac{1/(10^{-6})^2}{1/1^2} = 10^{12}$$

twelve orders of magnitude between rows of one matrix, which is why the
documentation flags a conditioning risk for tiny sigmas.

**Example 10: the zero-injection repair.** The demonstration: demote B7
to a passive bus (drop its load in a local copy), leave it dark with
the broken-tree set of Example 8, and watch the FREE zero-injection
knowledge repair observability:

````@example workshop_state_estimation
# uses: net (Example 1, Solve the reference power flow, Derive a noisy measurement set, Run the estimator, Example 2, Example 7, Example 8, Example 9); tree (Example 7)
@isdefined(net) || error("Run Example 1, the section \"Solve the reference power flow\", the section \"Derive a noisy measurement set\", the section \"Run the estimator\", Example 2, Example 7, Example 8, Example 9 first: they set up net.")
@isdefined(tree) || error("Run Example 7 first: it sets up tree.")
net_p = deepcopy(net)
removeProsumer!(net = net_p, busName = "B7", type = "LOAD")
empty!(net_p.measurements)
_, erg_p = runpf!(net_p, 40, 1e-10, 0)
erg_p == 0 || error("passive-bus reference power flow did not converge")
calcNetLosses!(net_p)
for (f, t) in tree[1:5]   ## broken tree again: B7 unmetered
  addPflowMeasurement!(net_p; fromBus = f, toBus = t, value = get_branch_p_from_to_mw(net_p, f, t), sigma = 0.8, direction = :from)
  addQflowMeasurement!(net_p; fromBus = f, toBus = t, value = get_branch_q_from_to_mvar(net_p, f, t), sigma = 0.8, direction = :from)
end
addVmMeasurement!(net_p; busName = "B1", value = net_p.nodeVec[net_p.busDict["B1"]]._vm_pu, sigma = 0.002)

g_no_zib = with_state_estimation_config(flatstart = true, jac_eps = 1e-6) do
  evaluate_global_observability(net_p)
end
println("without ZIB rows: quality = ", g_no_zib.quality, ", dof = m - n = ", g_no_zib.dof, ", unobservable state columns = ", g_no_zib.unobservable_state_columns)

added = addZeroInjectionMeasurements!(net_p; sigma = 1e-6)   ## auto-detects the passive B7
println("zero-injection rows added at the passive bus: ", length(added))
g_zib = with_state_estimation_config(flatstart = true, jac_eps = 1e-6) do
  evaluate_global_observability(net_p)
end
println("with ZIB rows:    quality = ", g_zib.quality, ", dof = m - n = ", g_zib.dof)
````

Reading aid (Example 10): the ZIB pair acts exactly like the injection
repair of Example 9, it couples B7 to its neighbors, but it costs
nothing: the information
was in the model all along. The dof moves from -2 (a SHORTFALL of two
equations) to 0.

**Example 11: the sigma trade-off.** How much does the sigma matter?
The same set as Example 10, estimated twice: once with
the near-constraint sigma and once with a soft one:

````@example workshop_state_estimation
# uses: net_p (Example 10)
@isdefined(net_p) || error("Run Example 10 first: it sets up net_p.")
for zib_sigma in (1e-6, 1e-2)
  netv = deepcopy(net_p)
  empty!(netv.measurements)
  for m in net_p.measurements
    m.id !== nothing && startswith(something(m.id, ""), "ZI") && continue
    push!(netv.measurements, m)
  end
  addZeroInjectionMeasurements!(netv; sigma = zib_sigma)
  sev = with_state_estimation_config(max_iter = 15, tol = 1e-6, flatstart = true, jac_eps = 1e-6, update_net = false) do
    runse!(netv)
  end
  vm_b7 = abs(sev.voltages[netv.busDict["B7"]])
  println("ZIB rows with sigma = ", zib_sigma, " MW / MVar:")
  println("  converged:                                            ", sev.converged)
  @printf("  objective J of the estimate [-]:                      %.1e\n", sev.objectiveJ)
  @printf("  Vm at B7, estimated [pu]:                             %.5f\n", vm_b7)
  @printf("  weight ratio ZIB row / flow row (sigma 0.8 MW) [-]:   %.1e\n", (0.8 / zib_sigma)^2)
end
````

Reading aid (Example 11): with clean synthetic data both sigmas land on the same
state (the constraint is consistent with the flows), and $J$ stays near
zero. The difference is the CONDITIONING headroom: the weight ratio of
the normal equations grows from about 6.4e3 at sigma 1e-2 to 6.4e11 at
sigma 1e-6, and every order of magnitude eats digits in $G$'s
factorization. Rule of thumb: as small as the constraint needs, as
large as the conditioning allows; with noisy real data a too-soft sigma
lets the passive bus drift, a too-hard one amplifies rounding.

### PMU phasors and the reference-angle offset

**Example 12: PMU phasors and the reference-angle offset.** A PMU
measures the voltage PHASOR: magnitude plus absolute angle,
GPS-synchronized. In the estimator that is a tightly weighted `VmMeas`
and a `VaMeas` in degrees (`addPmuPhasorMeasurement!` adds the pair).
One subtlety makes PMU angles interesting: the PMU time base rarely
coincides with the slack reference, so all PMU angles share a common
unknown OFFSET. With `state_estimation.pmu_ref_offset = auto` (the
default) the estimator adds that offset as an extra state and solves for
it; watch `n_states` grow from 13 to 14. We simulate two PMUs whose time
base is shifted by exactly 2 degrees, added on top of the repaired
measurement set from Example 9:

````@example workshop_state_estimation
# uses: net (Example 1, Solve the reference power flow, Derive a noisy measurement set, Run the estimator, Example 2, Example 7, Example 8, Example 9)
@isdefined(net) || error("Run Example 1, the section \"Solve the reference power flow\", the section \"Derive a noisy measurement set\", the section \"Run the estimator\", Example 2, Example 7, Example 8, Example 9 first: they set up net.")
pmu_shift_deg = 2.0
for bus in ("B4", "B6")
  idx = net.busDict[bus]
  addPmuPhasorMeasurement!(net; busName = bus, vm_pu = net.nodeVec[idx]._vm_pu, va_deg = net.nodeVec[idx]._va_deg + pmu_shift_deg, sigmaVm = 0.002, sigmaVa = 0.02)
end
gpmu = with_state_estimation_config(flatstart = true, jac_eps = 1e-6) do
  evaluate_global_observability(net)
end
println("with 2 PMUs: quality = ", gpmu.quality, " (measurement rows m = ", gpmu.n_measurements, ", state variables n = ", gpmu.n_states, ", offset state included)")

se_pmu = with_state_estimation_config(max_iter = 15, tol = 1e-6, flatstart = true, jac_eps = 1e-6, update_net = false) do
  runse!(net)
end
println("SE converged: ", se_pmu.converged)
println("PMU reference-angle offset [deg]:")
@printf("  true (time-base shift put on the PMU angles):   %.3f\n", pmu_shift_deg)
@printf("  estimated by SE:                                %.3f\n", se_pmu.vaRefOffsetDeg)
````

Reading aid (Example 12): the estimator recovers the 2° time-base shift as the extra
state instead of bending the bus angles toward it, and the estimated
network state stays anchored to the slack reference. Without the offset
state (`pmu_ref_offset = off`) the shifted PMU angles would fight the
flow measurements and inflate $J$. Theory and the measurement model:
[State Estimation](https://welthulk.github.io/Sparlectra.jl/state_estimation/).

## When the set is not observable: islands, closing, placement

Examples 8 to 10 repaired the broken corner by hand. The library can do
the three steps itself: find the observable islands, estimate them one
by one or close the gap with scheduled values, and suggest where a meter
is missing. All four examples start from the broken tree of Example 8.

**Example 13: observable islands.** A NETWORK island is a topologically
separate part with its own reference; the workshop net is one. An
OBSERVABLE island lies inside a network island: the buses whose states
the measurements fix once one of them is the reference, labeled "k.j"
(observable island j of network island k). Not observable is decided on
the decoupled linear models: P injections, P flows and PMU angles
against the angles (P-theta), Q rows and Vm against the magnitudes
(Q-V). The deficit counts free directions, the island count counts
pieces; the two are reported separately.

````@example workshop_state_estimation
# uses: net (Example 1, Solve the reference power flow, Derive a noisy measurement set, Run the estimator, Example 2, Example 7, Example 8, Example 9, Example 12); tree (Example 7)
@isdefined(net) || error("Run Example 1, the section \"Solve the reference power flow\", the section \"Derive a noisy measurement set\", the section \"Run the estimator\", Example 2, Example 7, Example 8, Example 9, Example 12 first: they set up net.")
@isdefined(tree) || error("Run Example 7 first: it sets up tree.")
empty!(net.measurements)
for (f, t) in tree[1:5]   ## the broken tree of Example 8
  addPflowMeasurement!(net; fromBus = f, toBus = t, value = get_branch_p_from_to_mw(net, f, t), sigma = 0.8, direction = :from)
  addQflowMeasurement!(net; fromBus = f, toBus = t, value = get_branch_q_from_to_mvar(net, f, t), sigma = 0.8, direction = :from)
end
addVmMeasurement!(net; busName = "B1", value = net.nodeVec[net.busDict["B1"]]._vm_pu, sigma = 0.002)
meas8 = Measurement[deepcopy(m) for m in net.measurements]   ## kept for the placement in Example 16
isl = with_state_estimation_config(flatstart = true, jac_eps = 1e-6) do
  identify_observable_islands(net)
end
print_observable_islands(isl)
println("observable islands: ", isl.n_islands, " (", join(isl.island_labels, ", "), "), B7 alone: ", isl.island_buses[2] == [net.busDict["B7"]])
println("rank deficit: ", isl.rank_deficit, " (angles ", isl.angle_rank_deficit, ", magnitudes ", isl.magnitude_rank_deficit, ")")
println("unobservable branches: ", join(isl.unobservable_branch_names, ", "))
# the chords carry no meter, yet their flows follow from the angles of
# island 1.1
metered = Set(m.branchIdx for m in net.measurements if m.branchIdx !== nothing)
unmetered = [k for k in eachindex(net.branchVec) if !(k in metered)]
observable_unmetered = [k for k in unmetered if !(k in isl.unobservable_branches)]
println("unmeasured, but observable: ", join(["B$(net.branchVec[k].fromBus) - B$(net.branchVec[k].toBus)" for k in observable_unmetered], ", "))
````

Reading aid (Example 13): one network island, two observable islands.
B7 is cut off on both sides, so B6-B7 and B7-B1 are unobservable; the
chords B2-B5 and B3-B6 carry no meter but are observable, because both
their ends lie in island 1.1. "Unmeasured" is not "unobservable". The
deficit is 2 (the angle and the magnitude of B7), the island count is 2
as well, but the two numbers answer different questions: a larger gap
can leave many single-bus islands with a small deficit, or the
reverse.

**Example 14: estimating the observable islands one by one.** With
`state_estimation.observability.on_unobservable = observable_islands`
`runse!` estimates every observable island on its own subnet. An island
without the slack takes its own reference, and its angles are relative
to it. A bus whose magnitude is dark is no estimated island, so B7 first
gets its Vm meter:

````@example workshop_state_estimation
# uses: net (Example 1, Solve the reference power flow, Derive a noisy measurement set, Run the estimator, Example 2, Example 7, Example 8, Example 9, Example 12, Example 13); va_true, vm_true (Solve the reference power flow)
@isdefined(net) || error("Run Example 1, the section \"Solve the reference power flow\", the section \"Derive a noisy measurement set\", the section \"Run the estimator\", Example 2, Example 7, Example 8, Example 9, Example 12, Example 13 first: they set up net.")
(@isdefined(va_true) && @isdefined(vm_true)) || error("Run the section \"Solve the reference power flow\" first: it sets up va_true, vm_true.")
addVmMeasurement!(net; busName = "B7", value = net.nodeVec[net.busDict["B7"]]._vm_pu, sigma = 0.002)
net_i = deepcopy(net)
se_i = with_state_estimation_config(max_iter = 15, tol = 1e-6, flatstart = true, jac_eps = 1e-6, update_net = true, observability = Sparlectra.ObservabilityConfig(on_unobservable = :observable_islands)) do
  runse!(net_i)
end
print_observable_islands(se_i)
oi = se_i.observableIslands
dvm = maximum(abs(abs(se_i.voltages[b]) - vm_true[b]) for b in 1:6)
dva = maximum(abs(rad2deg(angle(se_i.voltages[b])) - va_true[b]) for b in 1:6)
@printf("B1 to B6 against the truth: max |dVm| %.1e pu, max |dVa| %.1e deg\n", dvm, dva)
b7 = net.busDict["B7"]
println("B7: island ", oi.label_of_bus[b7], ", relative angle: ", oi.relative_angle[b7])
p67 = get_branch_p_from_to_mw(net_i, "B6", "B7")
p71 = get_branch_p_from_to_mw(net_i, "B7", "B1")
println("flows on the unobservable branches B6-B7 and B7-B1 [MW]: ", p67, ", ", p71)
````

Reading aid (Example 14): island 1.1 reproduces the truth, island 1.2 is
B7 alone with B7 as its own reference: its magnitude comes from the
meter, its angle says nothing about B1. The flows between the islands
stay empty (NaN), not zero: no measurement fixes them, and a zero would
look like a result.

**Example 15: closing with scheduled values.** `on_unobservable =
scheduled` first fills the deficit with replacement values: zero
injections at passive buses, then the scheduled P and Q of the
prosumers, then U from the stored state, a set point, the median of the
level or the nominal voltage. Each value gets 30 % of its size as sigma
(`restore_sigma_rel`), at least the sigma of a meter. A candidate is
taken only when it lowers the deficit; the rows are marked as
replacement values (`origin = :pseudo_restore`) and never eliminated as
bad data:

````@example workshop_state_estimation
# uses: net (Example 1, Solve the reference power flow, Derive a noisy measurement set, Run the estimator, Example 2, Example 7, Example 8, Example 9, Example 12, Example 13, Example 14); va_true (Solve the reference power flow)
@isdefined(net) || error("Run Example 1, the section \"Solve the reference power flow\", the section \"Derive a noisy measurement set\", the section \"Run the estimator\", Example 2, Example 7, Example 8, Example 9, Example 12, Example 13, Example 14 first: they set up net.")
@isdefined(va_true) || error("Run the section \"Solve the reference power flow\" first: it sets up va_true.")
schedule_run(n) = with_state_estimation_config(max_iter = 15, tol = 1e-6, flatstart = true, jac_eps = 1e-6, update_net = true, observability = Sparlectra.ObservabilityConfig(on_unobservable = :scheduled)) do
  runse!(n)
end
net_s = deepcopy(net)
se_s = schedule_run(net_s)
for r in se_s.restoredRows
  @printf("restored row: %s at %s, value %.1f MW, sigma %.1f MW, source %s\n", r.quantity, r.bus_name, r.value, r.sigma, r.source)
end
rrow = only(se_s.restoredRows)
println("estimated as one island: ", se_s.observableIslands === nothing, ", converged: ", se_s.converged)
dva_s = maximum(abs(rad2deg(angle(se_s.voltages[b])) - va_true[b]) for b in eachindex(va_true))
@printf("all buses against the truth: max |dVa| %.1e deg\n", dva_s)
````

The taken value is the only row that sees the flow into B7, so it is
critical: whatever the schedule says, the estimate follows. Shift the
schedule of the load at that bus by 10 MW and estimate again:

````@example workshop_state_estimation
# uses: net (Example 1, Solve the reference power flow, Derive a noisy measurement set, Run the estimator, Example 2, Example 7, Example 8, Example 9, Example 12, Example 13, Example 14); net_s, rrow, schedule_run (Example 15)
@isdefined(net) || error("Run Example 1, the section \"Solve the reference power flow\", the section \"Derive a noisy measurement set\", the section \"Run the estimator\", Example 2, Example 7, Example 8, Example 9, Example 12, Example 13, Example 14 first: they set up net.")
(@isdefined(net_s) && @isdefined(rrow) && @isdefined(schedule_run)) || error("Run Example 15 first: it sets up net_s, rrow, schedule_run.")
net_w = deepcopy(net)
load = only(ps for ps in net_w.prosumpsVec if getPosumerBusIndex(ps) == rrow.bus)
updatePQ!(load, load.pVal + 10.0, nothing)
se_w = schedule_run(net_w)
@printf("replacement value at %s: %.1f MW (schedule shifted by 10 MW)\n", rrow.bus_name, only(se_w.restoredRows).value)
dp = get_branch_p_from_to_mw(net_w, "B6", "B7") - get_branch_p_from_to_mw(net_s, "B6", "B7")
@printf("estimated flow B6-B7 moves by %.2f MW\n", dp)
````

Reading aid (Example 15): the deficit was one angle, so one value closed
it, the P schedule at B6 (the first bus by number whose injection sees
the cut flow). The restoration picks by rank, not by place, so the
replacement value need not sit at the dark bus B7: an injection at a
neighbor reaches its states as well (Example 9). The Web UI shows such a run yellow, "observable through 1
restored row". A replacement value carries no redundancy: the 10 MW
error of the schedule goes fully into the flow B6-B7, and no residual
shows it.

**Example 16: suggesting meter places.** `suggest_measurement_placement`
adds one meter per step until the set is observable, from candidates of
one kind: `injection` (P, Q or U at a bus), `flow` (a P and Q pair at a
branch end) or `pmu` (a voltage angle). For the one-row kinds the number
of steps equals the deficit (plus one step for the offset of a first
PMU), the minimum; for `flow` the choice is
greedy, without a guarantee of the minimum. Back to the set of
Example 8:

````@example workshop_state_estimation
# uses: meas8 (Example 13); net (Example 1, Solve the reference power flow, Derive a noisy measurement set, Run the estimator, Example 2, Example 7, Example 8, Example 9, Example 12, Example 13, Example 14)
@isdefined(meas8) || error("Run Example 13 first: it sets up meas8.")
@isdefined(net) || error("Run Example 1, the section \"Solve the reference power flow\", the section \"Derive a noisy measurement set\", the section \"Run the estimator\", Example 2, Example 7, Example 8, Example 9, Example 12, Example 13, Example 14 first: they set up net.")
plc = Dict{Symbol,Any}()
for kind in (:injection, :flow, :pmu)
  plc[kind] = with_state_estimation_config(flatstart = true, jac_eps = 1e-6) do
    suggest_measurement_placement(net, meas8; kind = kind)
  end
  r = plc[kind]
  println(kind, ": ", length(r.steps), " step(s), observable: ", r.reached, " (stop: ", r.stop_reason, ")")
  for s in r.steps
    println("  step ", s.step, ": ", s.quantity, " at ", s.location, ", deficit angles ", s.angle_rank_deficit_before, " -> ", s.angle_rank_deficit_after, ", magnitudes ", s.magnitude_rank_deficit_before, " -> ", s.magnitude_rank_deficit_after)
  end
end
steps_of(kind) = [(s.quantity, s.location, s.angle_rank_deficit_after, s.magnitude_rank_deficit_after) for s in plc[kind].steps]
# the suggestion of the first kind built in: rows at the suggested places
# (placeholder values, enough for the observability check)
g16 = with_state_estimation_config(flatstart = true, jac_eps = 1e-6) do
  evaluate_global_observability(net, vcat(meas8, plc[:injection].rows))
end
println("Example 8 set with the injection suggestion: quality = ", g16.quality, ", rank deficit ", g16.rank_deficit)
````

Reading aid (Example 16): two injection rows close the deficit of two,
one flow pair on B6-B7 closes both at once. The PMU list needs two steps
for one angle: the first PMU of a set only fixes the time-base offset of
Example 12, the second fixes the angle of B7. Its magnitude stays dark,
because a PMU here counts as an angle only; the placement stops and says
so. See `examples/state_estimation/state_estimation_observable_islands.jl`
for the same chain on sp_case118 and sp_case1354.

## Where to go next

- The sequel notebook,
  [bad data and parameter estimation](https://github.com/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_se_diagnostics.ipynb):
  what happens when $J/\mathrm{dof}$ is NOT near one: normalized
  residuals and localizability, sequential elimination, robust
  estimation, current measurements, and shunt-parameter estimation.
- [`state_estimation_observable_islands.jl`](https://github.com/Welthulk/Sparlectra.jl/blob/main/examples/state_estimation/state_estimation_observable_islands.jl):
  Examples 13 to 16 on sp_case118 and sp_case1354, with the run time per
  step.
- New to Sparlectra? The
  [workshop tour](https://github.com/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_tour.ipynb)
  builds a network from scratch step by step.
- [Slack types and short-circuit currents](https://github.com/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_slack_short_circuit.ipynb):
  the grid-connection notebook: slack representations and IEC 60909-0
  fault currents.
- [State Estimation](https://welthulk.github.io/Sparlectra.jl/state_estimation/):
  the full documentation: observability analysis, PMU angle measurements,
  bad-data handling, and the measurement model.
- [State-Estimation Configuration](https://welthulk.github.io/Sparlectra.jl/state_estimation_configuration/):
  every `state_estimation.*` configuration key.

