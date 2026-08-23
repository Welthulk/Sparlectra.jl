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
synthetic measurement noise reproducible; `LinearAlgebra` (standard
library) contributes `nullspace` for the observability deep dive.

````@example workshop_state_estimation
using Sparlectra
using Random
using LinearAlgebra
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

$J$ is more than a convergence number, it is a QUALITY measure: under
the assumed noise model, $J(\hat{x})$ is approximately chi-square
distributed with $\mathrm{dof} = m - n$ degrees of freedom, so its
expected value is about `dof`. A $J$ far ABOVE `dof` signals bad data or
a wrong model; far BELOW signals overfitted (too pessimistic) sigmas.
The ratio is the one number worth glancing at after every run:

````@example workshop_state_estimation
println("J / dof = ", round(se.objectiveJ / se.dof; digits = 3), "  (healthy noise: approximately 1)")

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

## Observability on paper: three small matrices

Before the network-sized deep dive, the whole theory on matrices small
enough to read. First the MINIMAL case: three measurements, three
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
println("H_B: observable = ", obs_B.numerical_observable, ", dof = ", obs_B.dof)
println("  critical rows: ", obs_B.numerical_critical_measurement_indices)
````

Now two extra rows, but BOTH duplicating the same information about
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
println("H_A: observable = ", obs_A.numerical_observable, ", dof = ", obs_A.dof)
println("  critical rows: ", obs_A.numerical_critical_measurement_indices)
for i in axes(H_A, 1)
  println("  row ", i, ": ", numerical_row_redundant(H_A, i) ? "redundant" : "CRITICAL", " (structural: ", structural_row_redundant(H_A, i) ? "redundant" : "CRITICAL", ")")
end
````

Finally the bridge to NETWORKS: an incidence-like matrix. Read each row
as a measurement on a 4-bus chain: a flow between buses a and b becomes
the pair +1/-1 in columns a and b (a flow sees a DIFFERENCE), an
injection at bus k becomes a single +1 (it anchors one state). This is
how measurement placement turns into sparsity structure:

```text
  row 1: flow 1-2    row 2: flow 2-3    row 3: flow 3-4
  row 4: injection at 2                 row 5: flow 1-3
```

````@example workshop_state_estimation
H_E = [
  1.0 -1.0 0.0 0.0
  0.0 1.0 -1.0 0.0
  0.0 0.0 1.0 -1.0
  0.0 1.0 0.0 0.0
  1.0 0.0 -1.0 0.0
]
obs_E = evaluate_observability_matrix(H_E)
println("H_E: observable = ", obs_E.numerical_observable, ", dof = ", obs_E.dof)
println("  critical rows: ", obs_E.numerical_critical_measurement_indices)
````

Reading aid: row 3 (the only path to bus 4) is critical; rows 1, 2, 5
form a triangle of alternatives around buses 1..3 and are redundant.
Exactly this pattern, at network scale, is what the deep dive below
builds and breaks. The extended version of these examples (with a toy
spanning-tree game and tolerance experiments) lives in
`examples/state_estimation/h_matrix_observability_demo.jl`.

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
drift without any measurement changing. The set of their nonzero
components partitions the network into OBSERVABLE ISLANDS (the term to
search the literature for is Monticelli's observable-island analysis).
Sparlectra reports exactly these components as
`unobservable_state_columns` on a not-observable global result.

One warning before using the LOCAL check: restricting $H$ to the
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
println("submatrix test on state 1: observable = ", loc.numerical_observable, "  (misleading)")
println("null space of H:           ", vec(round.(nullspace(H_flow); digits = 4)))
glob_flow = evaluate_observability_matrix(H_flow)
println("global dark states:        ", glob_flow.unobservable_state_columns, "  (the rigorous answer)")
````

The local check remains useful, it localizes candidate regions fast,
but a positive local verdict must be confirmed globally; the docstrings
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

The classical minimal recipe: a P/Q flow pair on every branch of a
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
pair at the bus (used in the repair below):

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
gmin = evaluate_global_observability(net; flatstart = true, jacEps = 1e-6)
println("minimal tree set: quality = ", gmin.quality, " (", gmin.n_measurements, " rows, ", gmin.n_states, " states)")
println("  structurally observable: ", gmin.structural_observable, ", numerically observable: ", gmin.numerical_observable, " (rank ", gmin.numerical_rank, ")")
println("  redundancy dof = ", gmin.dof, ", critical measurements: ", length(gmin.numerical_critical_measurement_indices), " of ", gmin.n_measurements)
````

Reading aid: observable, but `quality = :critical` and EVERY row is on
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
returns $H$ with described rows and labeled state columns, ready for a
placement report (the state-estimation example suite writes exactly such
a matrix-plus-stability page on every run):

````@example workshop_state_estimation
mj = measurement_jacobian(net)
println("H is ", size(mj.H, 1), " x ", size(mj.H, 2), "; columns: ", join(mj.cols[1:4], ", "), ", ...")
r1 = mj.rows[1]
println("row 1: ", r1.type, " at ", r1.location, " touches ", count(j -> abs(mj.H[1, j]) > 1e-9, eachindex(mj.cols)), " states")
````

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
# the rigorous per-state answer, straight from the null space: exactly
# B7's states are dark (angle column 6, magnitude column 13)
println("dark states (unobservable_state_columns): ", gbroken.unobservable_state_columns)

# local check on B2 (angle col 1, magnitude col 8): still fully covered
lb2 = evaluate_local_observability(net, [1, 8]; flatstart = true, jacEps = 1e-6)
println("local B2: numerically observable = ", lb2.numerical_observable)
# local check on B7 (angle col 6, magnitude col 13): NO measurement even
# touches these states anymore, which the check reports as an error;
# this is the EASY case for the local test, the hard case (touched but
# still dark) is the counterexample above, which only the global
# null-space answer catches
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

### Passive buses: zero-injection pseudo-measurements

A PASSIVE bus carries no load and no generation, so its injection is
exactly zero, and that knowledge is itself a measurement pair
$P = 0,\ Q = 0$, free of charge and free of noise. Sparlectra models it
as PSEUDO-measurements with a very small sigma
(`addZeroInjectionMeasurements!`), not as hard constraints. The
trade-off is honest: the smaller the sigma, the closer to a true
constraint, but the weight $1/\sigma^2$ explodes and degrades the
conditioning of $G$. The ratio between a ZIB row at $\sigma = 10^{-6}$
and a normal injection row at $\sigma = 1$ is

$$\frac{1/(10^{-6})^2}{1/1^2} = 10^{12}$$

twelve orders of magnitude between rows of one matrix, which is why the
documentation flags a conditioning risk for tiny sigmas.

Demonstration: demote B7 to a passive bus (drop its load in a local
copy), leave it dark with the broken-tree set, and watch the FREE
zero-injection knowledge repair observability:

````@example workshop_state_estimation
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

g_no_zib = evaluate_global_observability(net_p; flatstart = true, jacEps = 1e-6)
println("without ZIB: quality = ", g_no_zib.quality, ", dof = ", g_no_zib.dof, ", dark states = ", g_no_zib.unobservable_state_columns)

added = addZeroInjectionMeasurements!(net_p; sigma = 1e-6)   ## auto-detects the passive B7
println(length(added), " zero-injection rows added at the passive bus")
g_zib = evaluate_global_observability(net_p; flatstart = true, jacEps = 1e-6)
println("with ZIB:    quality = ", g_zib.quality, ", dof = ", g_zib.dof)
````

Reading aid: the ZIB pair acts exactly like the injection repair above,
it couples B7 to its neighbors, but it costs nothing: the information
was in the model all along. The dof moves from -2 (a SHORTFALL of two
equations) to 0.

How much does the sigma matter? The same set estimated twice, once with
the near-constraint sigma and once with a soft one:

````@example workshop_state_estimation
for zib_sigma in (1e-6, 1e-2)
  netv = deepcopy(net_p)
  empty!(netv.measurements)
  for m in net_p.measurements
    m.id !== nothing && startswith(something(m.id, ""), "ZI") && continue
    push!(netv.measurements, m)
  end
  addZeroInjectionMeasurements!(netv; sigma = zib_sigma)
  sev = runse!(netv; maxIte = 15, tol = 1e-6, flatstart = true, jacEps = 1e-6, updateNet = false)
  vm_b7 = abs(sev.voltages[netv.busDict["B7"]])
  println("ZIB sigma = ", zib_sigma, ": converged = ", sev.converged, ", J = ", round(sev.objectiveJ; digits = 6), ", Vm(B7) = ", round(vm_b7; digits = 5), " pu, max/min weight ratio = ", round((0.8 / zib_sigma)^2; digits = 1))
end
````

Reading aid: with clean synthetic data both sigmas land on the same
state (the constraint is consistent with the flows), and $J$ stays near
zero. The difference is the CONDITIONING headroom: the weight ratio of
the normal equations grows from about 6.4e3 at sigma 1e-2 to 6.4e11 at
sigma 1e-6, and every order of magnitude eats digits in $G$'s
factorization. Rule of thumb: as small as the constraint needs, as
large as the conditioning allows; with noisy real data a too-soft sigma
lets the passive bus drift, a too-hard one amplifies rounding.

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

