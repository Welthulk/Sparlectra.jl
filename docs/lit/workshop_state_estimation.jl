# Copyright 2023-2026 Udo Schmitz                                             #src
#                                                                             #src
# Licensed under the Apache License, Version 2.0 (the "License");             #src
# you may not use this file except in compliance with the License.            #src
# You may obtain a copy of the License at                                     #src
#                                                                             #src
#     http://www.apache.org/licenses/LICENSE-2.0                              #src
#                                                                             #src
# Unless required by applicable law or agreed to in writing, software         #src
# distributed under the License is distributed on an "AS IS" BASIS,           #src
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.    #src
# See the License for the specific language governing permissions and         #src
# limitations under the License.                                              #src
#                                                                             #src
# file: docs/lit/workshop_state_estimation.jl                                 #src
# purpose: Literate.jl source of the state-estimation workshop notebook and   #src
#          its Documenter page. Deep dive on the tour's estimation chapter;   #src
#          regenerate                                                          #src
#          the committed outputs with                                         #src
#          `julia --project=docs docs/generate_notebooks.jl`.                 #src

# # State estimation from noisy measurements
#
# > **Level: Advanced to Expert**, companion of the advanced tour's state-estimation chapter; the observability section is the advanced deep dive.
#
# > **Note:** This workshop was created with AI assistance and is reviewed
# > and curated by the maintainer; it is not a fully machine-generated text.
#
# A power flow computes the network state from an exact specification. A
# real control room has the opposite problem: it receives many **redundant,
# noisy measurements** (voltage magnitudes, injections, branch flows) and
# must reconstruct the most likely state from them. That is *state
# estimation*: a weighted-least-squares (WLS) fit of the bus voltages to
# the measurement set, where each measurement counts with the inverse of
# its variance.
#
# In this notebook you build a 7-bus network with
# [Sparlectra.jl](https://github.com/Welthulk/Sparlectra.jl), solve a
# reference power flow, derive a noisy synthetic measurement set from it,
# check observability, and run the estimator: the closed loop that lets
# you judge estimation quality against a known truth.
#
# The network is the workshop ring (B1 carries the grid connection, the
# diagonals are the cross-ties B2-B5 and B3-B6):
#
# ```text
#  (slack)
#    B1 ---- B2 ---- B3 ---- B4
#    |         \    /         |
#    |          \  /          |
#    |           \/           |
#    |           /\           |
#    B7 ---- B6 ---- B5 ------+
# ```

#nb # ## Setup
#nb # This notebook runs in the workshop folder that the installer sets up
#nb # (Julia, Sparlectra, IJulia and the notebooks in one command, see
#nb # [Workshops](https://welthulk.github.io/Sparlectra.jl/#Workshops)).
#nb # This cell checks that Sparlectra is installed and shows its version.
#nb using Pkg
#nb Base.find_package("Sparlectra") === nothing && error("Sparlectra is not installed in this environment; run the workshop installer, see https://welthulk.github.io/Sparlectra.jl/#Workshops")
#nb Pkg.status("Sparlectra")

# ## Warm-up
#
# Every printed number sits next to an @assert (the scenarios-workshop     #src
# pattern): the notebook cannot drift silently.        #src
# Julia compiles each function on first use. This cell loads the packages
# the chapters share. `Random` (standard library) provides the seeded
# generator that makes the synthetic measurement noise reproducible;
# `LinearAlgebra` (standard library) contributes `nullspace` for the
# observability deep dive; `Printf` (standard library) aligns the printed
# values. The notebook needs no helper functions; the
# chapters share the study network of Example 1 and build on each other's
# measurement sets, and a code cell that uses a name from another cell
# says so in a comment. Its last line, `warmup()`, runs every path the
# chapters exercise once on tiny throwaway networks (the power-flow solver
# that produces the reference state, the WLS estimator, the observability
# checks (global, local, on a matrix), the hand-placed and zero-injection
# measurements, the measurement Jacobian, the PMU phasors, the observable
# islands, the estimate per island, the restoration and the placement),
# so the real study runs at full speed; its code is in
# `docs/lit/warmup/workshop_state_estimation.jl` of the package. How long
# the warm-up takes depends on the machine; every later cell then runs without compile
# pauses.

using Sparlectra
using Random
using LinearAlgebra
using Printf

## compile every path the chapters use once; the warm-up code is in
## docs/lit/warmup/workshop_state_estimation.jl of the installed package
include(pkgdir(Sparlectra, "docs", "lit", "warmup", "workshop_state_estimation.jl"))
warmup()

# ## Build the study network
#
# **Example 1: the closed estimation loop, from network build to estimated
# state.** Seven 110 kV buses in a ring with two cross-connections, the
# ring drawn in the introduction above, enough meshing
# that branch-flow measurements carry real information. An external network
# injection at `B1` is the slack, a generator feeds at `B3`, the other
# buses carry loads. The sections up to the state inspection all belong to
# this example: reference power flow, noisy measurement set, observability
# check, estimation.

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

# ## Solve the reference power flow
#
# This solved state plays the role of the (in reality unknown) *true*
# system state: the measurements are derived from it, and the estimator
# never sees it directly.

## uses: net (Example 1)
@isdefined(net) || error("Run Example 1 first: it sets up net.")
ite_pf, status_pf = runpf!(net, 40, 1e-10, 0)
status_pf == 0 || error("Power flow did not converge")
println("reference power flow (true state): ", ite_pf, " iterations")
@assert ite_pf == 5                                                         #src
## the true state, kept for the comparison with the estimate below
vm_true = [net.nodeVec[i]._vm_pu for i in eachindex(net.nodeVec)]
va_true = [net.nodeVec[i]._va_deg for i in eachindex(net.nodeVec)]

# ## Derive a noisy measurement set
#
# `measurementStdDevs` defines one standard deviation per measurement kind
# (pu for voltages, MW/MVar for powers); `setMeasurementsFromPF!` then
# reads the solved state and creates the measurements, with Gaussian noise
# of exactly those standard deviations, seeded for reproducibility. This is
# the standard trick for exercising an estimator: the truth is known, so
# the residuals are meaningful.

## uses: net (Example 1, Solve the reference power flow)
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
@assert length(net.measurements) == 57                                      #src

# ## Check observability
#
# Estimation only works where the measurement set actually determines the
# state. Why that is a RANK question follows from what the estimator
# minimizes. With $m$ measurements $z_i$, their models $h_i(x)$, and
# accuracies $\sigma_i$, WLS minimizes the weighted squared residuals
#
# $$J(x) = \sum_{i=1}^{m} \frac{(z_i - h_i(x))^2}{\sigma_i^2}$$
#
# and every Newton step solves the normal equations
#
# $$G\,\Delta x = H^{\top} W r, \qquad G = H^{\top} W H, \qquad W = \mathrm{diag}(1/\sigma_i^2)$$
#
# with $H$ the measurement Jacobian. $G$ is invertible exactly when
# $\mathrm{rank}(H) = n$, the number of states: observability is not an
# extra condition on top of the estimator, it IS the solvability of these
# equations. For this workshop network $n = 2\,n_{bus} - 1 = 13$ (one
# magnitude per bus, one angle per non-slack bus), so the 13 appearing
# below is a formula, not a coincidence. The global check compares $m$
# against $n$ and probes the numerical rank of $H$.

## uses: net (Example 1, Solve the reference power flow, Derive a noisy measurement set)
@isdefined(net) || error("Run Example 1, the section \"Solve the reference power flow\", the section \"Derive a noisy measurement set\" first: they set up net.")
gobs = with_state_estimation_config(flatstart = true, jac_eps = 1e-6) do
  evaluate_global_observability(net)
end
println("global observability of the set: quality = ", gobs.quality)
@assert gobs.quality == :good                                               #src
println("  measurement rows m = ", gobs.n_measurements, ", state variables n = ", gobs.n_states)
@assert gobs.n_measurements == 57 && gobs.n_states == 13                    #src

# ## Run the estimator
#
# `runse!` iterates the WLS normal equations from a flat start. With
# `state_estimation.update_net = true` the estimated voltages are written
# back into the network, so the usual result printers show the *estimated*
# state. The estimator takes its settings from the active configuration;
# `with_state_estimation_config` installs other values for one run. The
# objective $J$ is the weighted sum of squared residuals; for healthy
# Gaussian noise it should land near the number of redundant measurements.

## uses: net (Example 1, Solve the reference power flow, Derive a noisy measurement set)
@isdefined(net) || error("Run Example 1, the section \"Solve the reference power flow\", the section \"Derive a noisy measurement set\" first: they set up net.")
se = with_state_estimation_config(max_iter = 12, tol = 1e-6, flatstart = true, jac_eps = 1e-6, update_net = true) do
  runse!(net)
end

println("SE converged: ", se.converged, " in ", se.iterations, " iterations")
@assert se.converged && se.iterations == 3                                  #src
@printf("objective J of the estimate [-]: %.2f\n", se.objectiveJ)
@assert 30 < se.objectiveJ < 70                                             #src

# ## Inspect the estimated state
#
# `SEResult.voltages` holds the estimated complex bus voltage per bus
# index. Because the measurements carried only mild noise, the estimate
# reproduces the reference power flow closely: magnitudes within 1e-3 pu,
# angles within 0.02°. The slack `B1` keeps its angle of 0° (the angle
# reference); its magnitude is estimated like every other one and lands
# near the true 1.02 pu. The degrees of freedom (`dof`) and the 3σ check on $J$ summarize whether the
# residuals are consistent with the declared measurement accuracies.
#
# $J$ is more than a convergence number, it is a QUALITY measure: under
# the assumed noise model, $J(\hat{x})$ is approximately chi-square
# distributed with $\mathrm{dof} = m - n$ degrees of freedom, so its
# expected value is about `dof`. A $J$ far ABOVE `dof` signals bad data or
# a wrong model; far BELOW signals overfitted (too pessimistic) sigmas.
# The ratio is the one number worth glancing at after every run:

## uses: net (Example 1, Solve the reference power flow, Derive a noisy measurement set, Run the estimator); se (Run the estimator); va_true, vm_true (Solve the reference power flow)
@isdefined(net) || error("Run Example 1, the section \"Solve the reference power flow\", the section \"Derive a noisy measurement set\", the section \"Run the estimator\" first: they set up net.")
@isdefined(se) || error("Run the section \"Run the estimator\" first: it sets up se.")
(@isdefined(va_true) && @isdefined(vm_true)) || error("Run the section \"Solve the reference power flow\" first: it sets up va_true, vm_true.")
@printf("J/dof of the estimate [-]:       %.3f   (healthy noise: approximately 1)\n", se.objectiveJ / se.dof)
@assert 0.5 < se.objectiveJ / se.dof < 2.0                                  #src

println("dof = m - n (redundant rows):    ", se.dof)
@assert se.dof == 44                                                        #src
println("J within 3σ band:                ", se.jWithin3Sigma)
@assert se.jWithin3Sigma === true                                           #src
println()
println("Bus voltages, true (reference power flow) and estimated:")
println("  bus   Vm true [pu]   Vm est [pu]   Va true [deg]   Va est [deg]")
for (name, idx) in sort(collect(net.busDict); by = last)
  v = se.voltages[idx]
  @printf("  %-4s %12.4f %13.4f %15.3f %14.3f\n", name, vm_true[idx], abs(v), va_true[idx], rad2deg(angle(v)))
  @assert 0.9 < abs(v) < 1.1 && abs(rad2deg(angle(v))) < 10                 #src
end

# ## Building measurement sets manually
#
# **Example 2: a deliberately sparse manual set.** Real measurement sets
# are not derived from a solved power flow; they
# arrive from SCADA one by one. The `add*Measurement!` helpers resolve bus
# and branch references for you, so assembling a set works just like
# building the network. A deliberately sparse set like this one is a good
# way to explore *when observability breaks down*:

## uses: net (Example 1, Solve the reference power flow, Derive a noisy measurement set, Run the estimator)
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
@assert obs.quality == :not_observable                                      #src
println("  measurement rows m = ", obs.n_measurements, ", state variables n = ", obs.n_states)
@assert obs.n_measurements == 5 && obs.n_states == 13                       #src

# Reading aid (Example 2): five measurements against 13 states: the check
# reports the shortfall instead of letting the estimator run into a
# rank-deficient system.
#
# ## Observability on paper: three small matrices
#
# Before the network-sized deep dive, the whole theory on matrices small
# enough to read: three numbered examples (Examples 3 to 5), each comment
# names the matrix it talks about.
#
# **Example 3 (`H_B`), the minimal case:** three measurements, three
# states, the identity. Observable (structurally: every column finds its
# own row in the matching; numerically: rank 3), but with $m = n$ every
# single row is CRITICAL, losing any one loses a state:

H_B = [
  1.0 0.0 0.0
  0.0 1.0 0.0
  0.0 0.0 1.0
]
obs_B = evaluate_observability_matrix(H_B)
println("H_B: numerically observable = ", obs_B.numerical_observable, ", dof = m - n = ", obs_B.dof)
@assert obs_B.numerical_observable === true && obs_B.dof == 0               #src
println("  critical rows of H_B: ", obs_B.numerical_critical_measurement_indices)
@assert obs_B.numerical_critical_measurement_indices == [1, 2, 3]           #src

# **Example 4 (`H_A`), duplicated information:** the same identity plus
# two extra rows, but BOTH duplicating the same information about
# states 1 and 2. The lesson: redundancy is PER STATE, not global. The
# set has dof = 2, yet row 3 stays critical, because it is still the only
# row that sees state 3; the two spare rows protect the wrong place:

H_A = [
  1.0 0.0 0.0
  0.0 1.0 0.0
  0.0 0.0 1.0
  1.0 1.0 0.0
  1.0 1.0 0.0
]
obs_A = evaluate_observability_matrix(H_A)
println("H_A: numerically observable = ", obs_A.numerical_observable, ", dof = m - n = ", obs_A.dof)
@assert obs_A.numerical_observable === true && obs_A.dof == 2               #src
println("  critical rows of H_A: ", obs_A.numerical_critical_measurement_indices)
@assert obs_A.numerical_critical_measurement_indices == [3]                 #src
for i in axes(H_A, 1)
  println("  row ", i, " of H_A: numerical check ", numerical_row_redundant(H_A, i) ? "redundant" : "CRITICAL", ", structural check ", structural_row_redundant(H_A, i) ? "redundant" : "CRITICAL")
  @assert numerical_row_redundant(H_A, i) == (i != 3)                       #src
end

# **Example 5 (`H_E`), the bridge to networks:** an incidence-like matrix
# for a 4-bus chain, one row per measurement (the listing below says
# which). In this
# idealized form you can READ the measurement type off the numbers: a
# flow between buses a and b sees only their DIFFERENCE, so its row is
# the pair +1 in column a and -1 in column b; an injection at bus k
# anchors one state, so its row has a single +1 in column k. Check it
# row by row: rows 1, 2, 3, 5 each carry a +1/-1 pair (flows), row 4 is
# the only single-entry row (the injection at bus 2).
#
# ```text
#   row 1: flow 1-2    row 2: flow 2-3    row 3: flow 3-4
#   row 4: injection at 2                 row 5: flow 1-3
# ```
#
# One honest caveat before transferring this: the legibility is a
# property of the TOY. In a real AC measurement Jacobian an injection row
# touches its bus AND every neighbor, a flow row touches three or four
# states, and the entries are sensitivities, not ones (the
# `measurement_jacobian` output below, in Example 7, says "PflowMeas at
# B1->B2 touches 3 states"). There the row DESCRIPTIONS carry the
# what-and-where; the numbers only carry the coupling.

H_E = [
  1.0 -1.0 0.0 0.0
  0.0 1.0 -1.0 0.0
  0.0 0.0 1.0 -1.0
  0.0 1.0 0.0 0.0
  1.0 0.0 -1.0 0.0
]
obs_E = evaluate_observability_matrix(H_E)
println("H_E: numerically observable = ", obs_E.numerical_observable, ", dof = m - n = ", obs_E.dof)
@assert obs_E.numerical_observable === true && obs_E.dof == 1               #src
println("  critical rows of H_E: ", obs_E.numerical_critical_measurement_indices)
@assert obs_E.numerical_critical_measurement_indices == [3, 4]              #src

# Reading aid (Example 5): row 3 (the only path to bus 4) is critical;
# row 4 (the injection at bus 2) is critical too, it is the only row that
# fixes an absolute level, the flow rows only see differences; rows 1, 2, 5
# form a triangle of alternatives around buses 1..3 and are redundant.
# Exactly this pattern, at network scale, is what the deep dive below
# builds and breaks (Examples 7 and 8). The extended version of these
# examples (with a toy spanning-tree game and tolerance experiments)
# lives in `examples/state_estimation/h_matrix_observability_demo.jl`.
#
# ## Observability deep dive: where measurements must sit (advanced)
#
# The checks above said WHETHER the state is observable. This section is
# about the WHY and the WHERE: what the numbers in the observability
# result mean, which placements make a network observable, and how the
# local check points at the corner that a missing meter leaves dark.
#
# First the two concepts in plain words:
#
# **Global observability** asks: does the measurement set determine the
# COMPLETE network state, every voltage magnitude and every angle, as one
# consistent picture? If yes, the WLS estimator has a unique solution; if
# no, whole regions of the state can drift without any measurement
# noticing, and the normal equations are singular. It is a property of the
# WHOLE set against the WHOLE state.
#
# **Local observability** asks the same question for a chosen PART of the
# state, typically the voltage and angle of one bus or one area: do the
# measurements that touch these states pin them down? A network can be
# globally unobservable while most of it is locally fine, and the local
# check is the tool that finds the dark corner. Think of it as zooming the
# rank question into a neighborhood.
#
# Both are decided by the same object, the measurement Jacobian $H$: one
# row per measurement, one column per state, entry = sensitivity of that
# measurement to that state. Global observability is full column rank of
# $H$; the sharpest way to say it uses the NULL SPACE:
#
# $$H \nu = 0,\ \nu \neq 0 \quad\Longrightarrow\quad \text{state } i \text{ is unobservable when } \nu_i \neq 0$$
#
# The null-space vectors are the dark regions: directions the state can
# drift without any measurement changing. Their nonzero components are the
# DARK STATES, which Sparlectra reports as `unobservable_state_columns` on
# a not-observable global result. Dark states and observable islands are
# two answers: the islands are the connected parts whose states are fixed
# once one of their buses is the reference (Monticelli and Wu), found by a
# separate analysis, see Example 13.
#
# **Example 6: why a local verdict can mislead.** One warning before
# using the LOCAL check: restricting $H$ to the
# columns you zoomed into (with the rows that touch them) is a NECESSARY
# test, not a sufficient one. The smallest counterexample: a single flow
# measurement between buses 1 and 2,
#
# $$H = \begin{pmatrix} 1 & -1 \end{pmatrix}$$
#
# The submatrix for state 1 is just $(1)$, full rank, so the local test
# says observable. But only the DIFFERENCE $x_1 - x_2$ is determined; the
# null-space vector $\nu = (1, 1)$ has a nonzero first component, so
# state 1 is not estimable. Both verdicts side by side:

H_flow = [1.0 -1.0]
loc = evaluate_local_observability_matrix(H_flow, [1])
println("local submatrix test on state 1: observable = ", loc.numerical_observable, "  (misleading)")
@assert loc.numerical_observable === true                                   #src
println("null-space vector of H (normalized, states 1 and 2): ", vec(round.(nullspace(H_flow); digits = 4)))
let ns = vec(nullspace(H_flow))
  @assert length(ns) == 2 && isapprox(abs(ns[1]), abs(ns[2]); atol = 1e-9)  #src
end
glob_flow = evaluate_observability_matrix(H_flow)
println("global check, unobservable state columns: ", glob_flow.unobservable_state_columns, "  (the rigorous answer)")
@assert glob_flow.unobservable_state_columns == [1, 2]                      #src

# Reading aid (Example 6): the local check remains useful, it localizes
# candidate regions fast, but a positive local verdict must be confirmed
# globally; the docstrings
# of both local helpers carry the same warning.
#
# The study network again, as a picture (ring plus two chords):
#
# ```text
#      B1 ──── B2 ──── B3
#     /        │        │ ╲
#   B7         │        │  B4
#     │        │        │  │
#     │        └─ B5 ───┼──┘
#     │           │     │      chords: B2─B5 and B3─B6
#     └─── B6 ────┴─────┘      ring:   B1-B2-B3-B4-B5-B6-B7-B1
# ```
#
# The STATE the estimator solves for has 13 entries: one voltage angle per
# non-slack bus (columns 1..6 for B2..B7; the slack B1 provides the angle
# reference) and one voltage magnitude per bus (columns 7..13 for
# B1..B7). Every measurement is one ROW of the measurement Jacobian $H$;
# the state is observable exactly when those rows span all 13 columns,
# i.e. $\mathrm{rank}(H) = 13$. `evaluate_global_observability` answers
# that twice: STRUCTURALLY (a bipartite matching on the sparsity pattern,
# `structural_matching` of 13 means every state column can be paired with
# its own measurement row) and NUMERICALLY (the actual rank; a placement
# can be structurally fine but numerically degenerate). The `quality`
# verdict compresses it: `:good` needs observability AND redundancy
# without critical single points; `:critical` is observable but one lost
# measurement would blind it; `:not_observable` speaks for itself.
#
# ### A minimal placement, built by hand
#
# **Example 7: a minimal placement, built by hand.** The classical
# minimal recipe: a P/Q flow pair on every branch of a
# SPANNING TREE (six branches for seven buses), plus one voltage-magnitude
# anchor. The flow pairs fix all relative angles and magnitude ratios
# along the tree, the anchor pins the absolute magnitude level, the slack
# pins the absolute angle. The counting identity says why this is exactly
# minimal:
#
# $$m = 2\,(n_{bus} - 1) + 1 = 2\,n_{bus} - 1 = n$$
#
# a spanning tree has $n_{bus} - 1$ branches, each contributing a P and a
# Q row, plus the one anchor; $m = n$ means observability with ZERO
# redundancy, every equation is needed:

## uses: net (Example 1, Solve the reference power flow, Derive a noisy measurement set, Run the estimator, Example 2)
@isdefined(net) || error("Run Example 1, the section \"Solve the reference power flow\", the section \"Derive a noisy measurement set\", the section \"Run the estimator\", Example 2 first: they set up net.")
empty!(net.measurements)
## re-solve the reference so the branch flows are available as true values
ite_ref, status_ref = runpf!(net, 40, 1e-10, 0)
status_ref == 0 || error("reference power flow did not converge")
calcNetLosses!(net)

tree = [("B1", "B2"), ("B2", "B3"), ("B3", "B4"), ("B4", "B5"), ("B5", "B6"), ("B6", "B7")]
for (f, t) in tree
  addPflowMeasurement!(net; fromBus = f, toBus = t, value = get_branch_p_from_to_mw(net, f, t), sigma = 0.8, direction = :from)
  addQflowMeasurement!(net; fromBus = f, toBus = t, value = get_branch_q_from_to_mvar(net, f, t), sigma = 0.8, direction = :from)
end
addVmMeasurement!(net; busName = "B1", value = net.nodeVec[net.busDict["B1"]]._vm_pu, sigma = 0.002)

# The same placement, drawn into the line diagram. Legend: `x` = voltage
# meter (Vm), `o` = P/Q flow pair on the branch, `(o)` = P/Q injection
# pair at the bus (used in the repair below, Example 9):
#
# ```text
#     [x]B1 ──o── B2 ──o── B3
#       /          │        │ ╲
#     B7           │        │  o
#      │           │        │   ╲
#      o           └─ B5 ───┼─── B4
#      │              o     │
#      └─── B6 ───────┴─────┘
#             (chords B2─B5, B3─B6 and ring closure B7─B1: UNMEASURED)
#
#     x  voltage meter        o  P/Q flow pair on the branch
# ```
#
# The picture already tells the placement story: the six `o` pairs walk a
# spanning tree (every bus is reached), the single `x` anchors the
# magnitude level, and the three unmeasured branches are exactly where
# redundancy would come from.

## uses: net (Example 1, Solve the reference power flow, Derive a noisy measurement set, Run the estimator, Example 2, Example 7)
@isdefined(net) || error("Run Example 1, the section \"Solve the reference power flow\", the section \"Derive a noisy measurement set\", the section \"Run the estimator\", Example 2, Example 7 first: they set up net.")
gmin = with_state_estimation_config(flatstart = true, jac_eps = 1e-6) do
  evaluate_global_observability(net)
end
println("minimal tree set: quality = ", gmin.quality, " (measurement rows m = ", gmin.n_measurements, ", state variables n = ", gmin.n_states, ")")
@assert gmin.quality == :critical && gmin.n_measurements == 13 && gmin.n_states == 13  #src
println("  structurally observable: ", gmin.structural_observable, ", numerically observable: ", gmin.numerical_observable, " (numerical rank of H = ", gmin.numerical_rank, ")")
@assert gmin.structural_observable && gmin.numerical_observable && gmin.numerical_rank == 13  #src
println("  redundancy dof = m - n = ", gmin.dof, ", critical measurement rows: ", length(gmin.numerical_critical_measurement_indices), " of ", gmin.n_measurements)
@assert gmin.dof == 0 && length(gmin.numerical_critical_measurement_indices) == 13  #src

# Reading aid (Example 7): observable, but `quality = :critical` and EVERY row is on
# the critical list: with zero redundancy, losing any single meter blinds
# some part of the state. There is a beautiful residual-side view of the
# same fact:
#
# $$\mathrm{dof} = m - n, \qquad S = I - H\,G^{-1} H^{\top} W, \qquad \text{measurement } i \text{ critical} \iff S_{ii} = 0$$
#
# $S$ maps measurements to their residuals; a critical measurement has
# residual EXACTLY zero, whatever its value, so a gross error in it is
# invisible to bad-data detection. That is why `dof` is what the bad-data
# machinery feeds on: only redundant measurements leave residual traces.
# One precision on Sparlectra's `dof` field: it is the COUNT difference
# $m - n$. For an observable set that equals $m - \mathrm{rank}(H)$; for
# an unobservable set it can be negative and then reads as a shortfall,
# not a redundancy. Real placements add the ring-closing and chord flows
# (or injections) precisely to buy this redundancy.
#
# The matrix behind all of this is one call away: `measurement_jacobian`
# returns $H$, here for the Example 7 placement, with described rows and
# labeled state columns, ready for a
# placement report (the state-estimation example suite writes exactly such
# a matrix-plus-stability page on every run):

## uses: net (Example 1, Solve the reference power flow, Derive a noisy measurement set, Run the estimator, Example 2, Example 7)
@isdefined(net) || error("Run Example 1, the section \"Solve the reference power flow\", the section \"Derive a noisy measurement set\", the section \"Run the estimator\", Example 2, Example 7 first: they set up net.")
mj = measurement_jacobian(net)
println("measurement Jacobian H: ", size(mj.H, 1), " measurement rows x ", size(mj.H, 2), " state columns; first columns: ", join(mj.cols[1:4], ", "), ", ...")
@assert size(mj.H) == (13, 13)                                              #src
r1 = mj.rows[1]
println("row 1 of H: ", r1.type, " at ", r1.location, " touches ", count(j -> abs(mj.H[1, j]) > 1e-9, eachindex(mj.cols)), " states (nonzero entries)")
@assert count(j -> abs(mj.H[1, j]) > 1e-9, eachindex(mj.cols)) == 3         #src
#
# ### Breaking a corner, and finding it with the local check
#
# **Example 8: breaking a corner.** Drop the B6-B7 flow pair from the
# Example 7 tree: the spur toward B7 loses its only meters.
# Globally the verdict flips to `:not_observable`; LOCALLY the check can
# say which states went dark. `evaluate_local_observability(net, cols)`
# restricts $H$ to selected state columns; for B7 those are angle column 6
# and magnitude column 13:

## uses: net (Example 1, Solve the reference power flow, Derive a noisy measurement set, Run the estimator, Example 2, Example 7); tree (Example 7)
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
@assert gbroken.quality == :not_observable && gbroken.numerical_rank == 11 && gbroken.n_states == 13  #src
## the rigorous per-state answer, straight from the null space: exactly
## B7's states are dark (angle column 6, magnitude column 13)
println("dark states (unobservable_state_columns): ", gbroken.unobservable_state_columns)
@assert gbroken.unobservable_state_columns == [6, 13]                       #src

## local check on B2 (angle col 1, magnitude col 8): still fully covered
lb2 = with_state_estimation_config(flatstart = true, jac_eps = 1e-6) do
  evaluate_local_observability(net, [1, 8])
end
println("local B2: numerically observable = ", lb2.numerical_observable)
@assert lb2.numerical_observable === true                                   #src
## local check on B7 (angle col 6, magnitude col 13): NO measurement even
## touches these states anymore, which the check reports as an error;
## this is the EASY case for the local test, the hard case (touched but
## still dark) is the Example 6 counterexample, which only the global
## null-space answer catches
try
  with_state_estimation_config(flatstart = true, jac_eps = 1e-6) do
    evaluate_local_observability(net, [6, 13])
  end
catch err
  println("local B7: ", sprint(showerror, err))
  @assert occursin("no active measurement", sprint(showerror, err))         #src
end

# ### Repairing with an injection
#
# **Example 9: repairing with an injection.** A meter does not have to
# sit ON the dark bus pair. Starting from the broken set of Example 8, an
# INJECTION measurement at B7 couples B7 to every neighbor (B6 and B1),
# because the
# injected power is the sum over its incident branches; its Jacobian row
# touches all their states. P and Q injection at B7 restore rank 13:

## uses: net (Example 1, Solve the reference power flow, Derive a noisy measurement set, Run the estimator, Example 2, Example 7, Example 8)
@isdefined(net) || error("Run Example 1, the section \"Solve the reference power flow\", the section \"Derive a noisy measurement set\", the section \"Run the estimator\", Example 2, Example 7, Example 8 first: they set up net.")
## B7 carries exactly its 20 MW / 6 MVAr load, so the true net injection
## is -20 MW / -6 MVAr (consumption counts negative)
addPinjMeasurement!(net; busName = "B7", value = -20.0, sigma = 1.0)
addQinjMeasurement!(net; busName = "B7", value = -6.0, sigma = 1.0)
grepaired = with_state_estimation_config(flatstart = true, jac_eps = 1e-6) do
  evaluate_global_observability(net)
end
println("with P/Q injection at B7: quality = ", grepaired.quality, " (numerical rank of H = ", grepaired.numerical_rank, " of n = ", grepaired.n_states, " states)")
@assert grepaired.quality == :critical && grepaired.numerical_rank == 13    #src

# Placement rules, condensed: flow pairs see the two ends of their branch,
# injections see the bus AND all its neighbors, `Vm` anchors the magnitude
# level (at least one required), and the slack anchors the angles. A bus
# is dark exactly when no row of any of these reaches its columns.
#
# ### Passive buses: zero-injection pseudo-measurements
#
# A PASSIVE bus carries no load and no generation, so its injection is
# exactly zero, and that knowledge is itself a measurement pair
# $P = 0,\ Q = 0$, free of charge and free of noise. Sparlectra models it
# as PSEUDO-measurements with a very small sigma
# (`addZeroInjectionMeasurements!`), not as hard constraints. With
# `state_estimation.zero_injection = constraint` the same rows become hard
# equality constraints instead. The trade-off of the default is honest:
# the smaller the sigma, the closer to a true constraint, but the weight
# $1/\sigma^2$ explodes and degrades the conditioning of $G$. The ratio
# between a ZIB row at $\sigma = 10^{-6}$ and a normal injection row at
# $\sigma = 1$ is
#
# $$\frac{1/(10^{-6})^2}{1/1^2} = 10^{12}$$
#
# twelve orders of magnitude between rows of one matrix, which is why the
# documentation flags a conditioning risk for tiny sigmas.
#
# **Example 10: the zero-injection repair.** The demonstration: demote B7
# to a passive bus (drop its load in a local copy), leave it dark with
# the broken-tree set of Example 8, and watch the FREE zero-injection
# knowledge repair observability:

## uses: net (Example 1, Solve the reference power flow, Derive a noisy measurement set, Run the estimator, Example 2, Example 7, Example 8, Example 9); tree (Example 7)
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
@assert g_no_zib.quality == :not_observable && g_no_zib.dof == -2 && g_no_zib.unobservable_state_columns == [6, 13]  #src

added = addZeroInjectionMeasurements!(net_p; sigma = 1e-6)   ## auto-detects the passive B7
println("zero-injection rows added at the passive bus: ", length(added))
@assert length(added) == 2                                                  #src
g_zib = with_state_estimation_config(flatstart = true, jac_eps = 1e-6) do
  evaluate_global_observability(net_p)
end
println("with ZIB rows:    quality = ", g_zib.quality, ", dof = m - n = ", g_zib.dof)
@assert g_zib.quality == :critical && g_zib.dof == 0                        #src

# Reading aid (Example 10): the ZIB pair acts exactly like the injection
# repair of Example 9, it couples B7 to its neighbors, but it costs
# nothing: the information
# was in the model all along. The dof moves from -2 (a SHORTFALL of two
# equations) to 0.
#
# **Example 11: the sigma trade-off.** How much does the sigma matter?
# The same set as Example 10, estimated twice: once with
# the near-constraint sigma and once with a soft one:

## uses: net_p (Example 10)
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
  @assert sev.converged && sev.objectiveJ < 1e-8 && isapprox(vm_b7, 1.00672; atol = 1e-3)  #src
end

# Reading aid (Example 11): with clean synthetic data both sigmas land on the same
# state (the constraint is consistent with the flows), and $J$ stays near
# zero. The difference is the CONDITIONING headroom: the weight ratio of
# the normal equations grows from about 6.4e3 at sigma 1e-2 to 6.4e11 at
# sigma 1e-6, and every order of magnitude eats digits in $G$'s
# factorization. Rule of thumb: as small as the constraint needs, as
# large as the conditioning allows; with noisy real data a too-soft sigma
# lets the passive bus drift, a too-hard one amplifies rounding.
#
# ### PMU phasors and the reference-angle offset
#
# **Example 12: PMU phasors and the reference-angle offset.** A PMU
# measures the voltage PHASOR: magnitude plus absolute angle,
# GPS-synchronized. In the estimator that is a tightly weighted `VmMeas`
# and a `VaMeas` in degrees (`addPmuPhasorMeasurement!` adds the pair).
# One subtlety makes PMU angles interesting: the PMU time base rarely
# coincides with the slack reference, so all PMU angles share a common
# unknown OFFSET. With `state_estimation.pmu_ref_offset = auto` (the
# default) the estimator adds that offset as an extra state and solves for
# it; watch `n_states` grow from 13 to 14. We simulate two PMUs whose time
# base is shifted by exactly 2 degrees, added on top of the repaired
# measurement set from Example 9:

## uses: net (Example 1, Solve the reference power flow, Derive a noisy measurement set, Run the estimator, Example 2, Example 7, Example 8, Example 9)
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
@assert gpmu.quality == :critical && gpmu.n_measurements == 17 && gpmu.n_states == 14  #src

se_pmu = with_state_estimation_config(max_iter = 15, tol = 1e-6, flatstart = true, jac_eps = 1e-6, update_net = false) do
  runse!(net)
end
println("SE converged: ", se_pmu.converged)
println("PMU reference-angle offset [deg]:")
@printf("  true (time-base shift put on the PMU angles):   %.3f\n", pmu_shift_deg)
@printf("  estimated by SE:                                %.3f\n", se_pmu.vaRefOffsetDeg)
@assert se_pmu.converged && isapprox(se_pmu.vaRefOffsetDeg, pmu_shift_deg; atol = 0.05)  #src

# Reading aid (Example 12): the estimator recovers the 2° time-base shift as the extra
# state instead of bending the bus angles toward it, and the estimated
# network state stays anchored to the slack reference. Without the offset
# state (`pmu_ref_offset = off`) the shifted PMU angles would fight the
# flow measurements and inflate $J$. Theory and the measurement model:
# [State Estimation](https://welthulk.github.io/Sparlectra.jl/state_estimation/).
#
# ## When the set is not observable: islands, closing, placement
#
# Examples 8 to 10 repaired the broken corner by hand. The library can do
# the three steps itself: find the observable islands, estimate them one
# by one or close the gap with scheduled values, and suggest where a meter
# is missing. All four examples start from the broken tree of Example 8.
#
# **Example 13: observable islands.** A NETWORK island is a topologically
# separate part with its own reference; the workshop net is one. An
# OBSERVABLE island lies inside a network island: the buses whose states
# the measurements fix once one of them is the reference, labeled "k.j"
# (observable island j of network island k). Not observable is decided on
# the decoupled linear models: P injections, P flows and PMU angles
# against the angles (P-theta), Q rows and Vm against the magnitudes
# (Q-V). The deficit counts free directions, the island count counts
# pieces; the two are reported separately.

## uses: net (Example 1, Solve the reference power flow, Derive a noisy measurement set, Run the estimator, Example 2, Example 7, Example 8, Example 9, Example 12); tree (Example 7)
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
@assert length(isl.island_buses[1]) == 6 && isl.passes == 1 && isl.region_buses == isl.dark_magnitude_buses == [net.busDict["B7"]]  #src
println("observable islands: ", isl.n_islands, " (", join(isl.island_labels, ", "), "), B7 alone: ", isl.island_buses[2] == [net.busDict["B7"]])
@assert isl.n_islands == 2 && isl.island_labels == ["1.1", "1.2"] && isl.island_buses[2] == [net.busDict["B7"]]  #src
println("rank deficit: ", isl.rank_deficit, " (angles ", isl.angle_rank_deficit, ", magnitudes ", isl.magnitude_rank_deficit, ")")
@assert (isl.rank_deficit, isl.angle_rank_deficit, isl.magnitude_rank_deficit) == (2, 1, 1)  #src
println("unobservable branches: ", join(isl.unobservable_branch_names, ", "))
@assert isl.unobservable_branch_names == ["B6 - B7", "B7 - B1"]             #src
## the chords carry no meter, yet their flows follow from the angles of
## island 1.1
metered = Set(m.branchIdx for m in net.measurements if m.branchIdx !== nothing)
unmetered = [k for k in eachindex(net.branchVec) if !(k in metered)]
observable_unmetered = [k for k in unmetered if !(k in isl.unobservable_branches)]
println("unmeasured, but observable: ", join(["B$(net.branchVec[k].fromBus) - B$(net.branchVec[k].toBus)" for k in observable_unmetered], ", "))
@assert [(net.branchVec[k].fromBus, net.branchVec[k].toBus) for k in observable_unmetered] == [(2, 5), (3, 6)]  #src

# Reading aid (Example 13): one network island, two observable islands.
# B7 is cut off on both sides, so B6-B7 and B7-B1 are unobservable; the
# chords B2-B5 and B3-B6 carry no meter but are observable, because both
# their ends lie in island 1.1. "Unmeasured" is not "unobservable". The
# deficit is 2 (the angle and the magnitude of B7), the island count is 2
# as well, but the two numbers answer different questions: a larger gap
# can leave many single-bus islands with a small deficit, or the
# reverse.
#
# **Example 14: estimating the observable islands one by one.** With
# `state_estimation.observability.on_unobservable = observable_islands`
# `runse!` estimates every observable island on its own subnet. An island
# without the slack takes its own reference, and its angles are relative
# to it. A bus whose magnitude is dark is no estimated island, so B7 first
# gets its Vm meter:

## uses: net (Example 1, Solve the reference power flow, Derive a noisy measurement set, Run the estimator, Example 2, Example 7, Example 8, Example 9, Example 12, Example 13); va_true, vm_true (Solve the reference power flow)
@isdefined(net) || error("Run Example 1, the section \"Solve the reference power flow\", the section \"Derive a noisy measurement set\", the section \"Run the estimator\", Example 2, Example 7, Example 8, Example 9, Example 12, Example 13 first: they set up net.")
(@isdefined(va_true) && @isdefined(vm_true)) || error("Run the section \"Solve the reference power flow\" first: it sets up va_true, vm_true.")
addVmMeasurement!(net; busName = "B7", value = net.nodeVec[net.busDict["B7"]]._vm_pu, sigma = 0.002)
net_i = deepcopy(net)
se_i = with_state_estimation_config(max_iter = 15, tol = 1e-6, flatstart = true, jac_eps = 1e-6, update_net = true, observability = Sparlectra.ObservabilityConfig(on_unobservable = :observable_islands)) do
  runse!(net_i)
end
print_observable_islands(se_i)
@assert se_i.converged && length(se_i.islands) == 2 && all(r -> r.objectiveJ < 1e-12 && r.dof == 0, se_i.islands)  #src
oi = se_i.observableIslands
dvm = maximum(abs(abs(se_i.voltages[b]) - vm_true[b]) for b in 1:6)
dva = maximum(abs(rad2deg(angle(se_i.voltages[b])) - va_true[b]) for b in 1:6)
@printf("B1 to B6 against the truth: max |dVm| %.1e pu, max |dVa| %.1e deg\n", dvm, dva)
@assert dvm < 1e-6 && dva < 1e-6                                            #src
b7 = net.busDict["B7"]
println("B7: island ", oi.label_of_bus[b7], ", relative angle: ", oi.relative_angle[b7])
@assert oi.label_of_bus[b7] == "1.2" && oi.relative_angle[b7] && !any(oi.relative_angle[1:6])  #src
p67 = get_branch_p_from_to_mw(net_i, "B6", "B7")
p71 = get_branch_p_from_to_mw(net_i, "B7", "B1")
println("flows on the unobservable branches B6-B7 and B7-B1 [MW]: ", p67, ", ", p71)
@assert isnan(p67) && isnan(p71)                                            #src

# Reading aid (Example 14): island 1.1 reproduces the truth, island 1.2 is
# B7 alone with B7 as its own reference: its magnitude comes from the
# meter, its angle says nothing about B1. The flows between the islands
# stay empty (NaN), not zero: no measurement fixes them, and a zero would
# look like a result.
#
# **Example 15: closing with scheduled values.** `on_unobservable =
# scheduled` first fills the deficit with replacement values: zero
# injections at passive buses, then the scheduled P and Q of the
# prosumers, then U from the stored state, a set point, the median of the
# level or the nominal voltage. Each value gets 30 % of its size as sigma
# (`restore_sigma_rel`), at least the sigma of a meter. A candidate is
# taken only when it lowers the deficit; the rows are marked as
# replacement values (`origin = :pseudo_restore`) and never eliminated as
# bad data:

## uses: net (Example 1, Solve the reference power flow, Derive a noisy measurement set, Run the estimator, Example 2, Example 7, Example 8, Example 9, Example 12, Example 13, Example 14); va_true (Solve the reference power flow)
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
@assert length(se_s.restoredRows) == 1                                      #src
rrow = only(se_s.restoredRows)
@assert rrow.bus_name == "B6" && rrow.quantity === :P && rrow.value == -30.0 && rrow.sigma == 9.0  #src
println("estimated as one island: ", se_s.observableIslands === nothing, ", converged: ", se_s.converged)
@assert se_s.observableIslands === nothing && se_s.converged                #src
dva_s = maximum(abs(rad2deg(angle(se_s.voltages[b])) - va_true[b]) for b in eachindex(va_true))
@printf("all buses against the truth: max |dVa| %.1e deg\n", dva_s)
@assert dva_s < 1e-6                                                        #src

# The taken value is the only row that sees the flow into B7, so it is
# critical: whatever the schedule says, the estimate follows. Shift the
# schedule of the load at that bus by 10 MW and estimate again:

## uses: net (Example 1, Solve the reference power flow, Derive a noisy measurement set, Run the estimator, Example 2, Example 7, Example 8, Example 9, Example 12, Example 13, Example 14); net_s, rrow, schedule_run (Example 15)
@isdefined(net) || error("Run Example 1, the section \"Solve the reference power flow\", the section \"Derive a noisy measurement set\", the section \"Run the estimator\", Example 2, Example 7, Example 8, Example 9, Example 12, Example 13, Example 14 first: they set up net.")
(@isdefined(net_s) && @isdefined(rrow) && @isdefined(schedule_run)) || error("Run Example 15 first: it sets up net_s, rrow, schedule_run.")
net_w = deepcopy(net)
load = only(ps for ps in net_w.prosumpsVec if getPosumerBusIndex(ps) == rrow.bus)
updatePQ!(load, load.pVal + 10.0, nothing)
se_w = schedule_run(net_w)
@printf("replacement value at %s: %.1f MW (schedule shifted by 10 MW)\n", rrow.bus_name, only(se_w.restoredRows).value)
@assert only(se_w.restoredRows).value == -40.0                              #src
dp = get_branch_p_from_to_mw(net_w, "B6", "B7") - get_branch_p_from_to_mw(net_s, "B6", "B7")
@printf("estimated flow B6-B7 moves by %.2f MW\n", dp)
@assert isapprox(dp, -10.0; atol = 0.1)                                     #src

# Reading aid (Example 15): the deficit was one angle, so one value closed
# it, the P schedule at B6 (the first bus by number whose injection sees
# the cut flow). The restoration picks by rank, not by place, so the
# replacement value need not sit at the dark bus B7: an injection at a
# neighbor reaches its states as well (Example 9). The Web UI shows such a run yellow, "observable through 1
# restored row". A replacement value carries no redundancy: the 10 MW
# error of the schedule goes fully into the flow B6-B7, and no residual
# shows it.
#
# **Example 16: suggesting meter places.** `suggest_measurement_placement`
# adds one meter per step until the set is observable, from candidates of
# one kind: `injection` (P, Q or U at a bus), `flow` (a P and Q pair at a
# branch end) or `pmu` (a voltage angle). For the one-row kinds the number
# of steps equals the deficit (plus one step for the offset of a first
# PMU), the minimum; for `flow` the choice is
# greedy, without a guarantee of the minimum. Back to the set of
# Example 8:

## uses: meas8 (Example 13); net (Example 1, Solve the reference power flow, Derive a noisy measurement set, Run the estimator, Example 2, Example 7, Example 8, Example 9, Example 12, Example 13, Example 14)
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
@assert plc[:injection].reached && steps_of(:injection) == [(:P, "B1", 0, 1), (:Q, "B1", 0, 0)]  #src
@assert plc[:flow].reached && steps_of(:flow) == [(:PQ, "B6 - B7 (from end)", 0, 0)]  #src
@assert !plc[:pmu].reached && plc[:pmu].stop_reason === :no_candidate && steps_of(:pmu) == [(:Va, "B1", 1, 1), (:Va, "B7", 0, 1)]  #src
@assert all(s.angle_rank_deficit_before == 1 && s.magnitude_rank_deficit_before == 1 for s in plc[:pmu].steps)  #src
@assert (plc[:injection].steps[2].angle_rank_deficit_before, plc[:injection].steps[2].magnitude_rank_deficit_before) == (0, 1)  #src
## the suggestion of the first kind built in: rows at the suggested places
## (placeholder values, enough for the observability check)
g16 = with_state_estimation_config(flatstart = true, jac_eps = 1e-6) do
  evaluate_global_observability(net, vcat(meas8, plc[:injection].rows))
end
println("Example 8 set with the injection suggestion: quality = ", g16.quality, ", rank deficit ", g16.rank_deficit)
@assert g16.quality != :not_observable && g16.rank_deficit == 0             #src

# Reading aid (Example 16): two injection rows close the deficit of two,
# one flow pair on B6-B7 closes both at once. The PMU list needs two steps
# for one angle: the first PMU of a set only fixes the time-base offset of
# Example 12, the second fixes the angle of B7. Its magnitude stays dark,
# because a PMU here counts as an angle only; the placement stops and says
# so. See `examples/state_estimation/state_estimation_observable_islands.jl`
# for the same chain on sp_case118 and sp_case1354.
#
# ## Several areas, one estimate
#
# **Example 17: multi-area estimation.** A large net is often run by
# several operators, each with its own meters. `runse_multi_area!` lets
# every area estimate with its own rows: a bus row belongs to the area of
# its bus, a flow row to the area of its measured end. A coordinator joins
# the areas at the border buses (the ends of the tie branches between two
# areas), pass by pass. The result is the integrated estimate on the same
# rows, which the result carries for comparison. Here area 1 is B1, B2 and
# B7, area 2 is B3 to B6, on a fresh copy of the net with a fresh noisy set:

## uses: net (Example 1, Solve the reference power flow, Derive a noisy measurement set, Run the estimator, Example 2, Example 7, Example 8, Example 9, Example 12, Example 13, Example 14)
@isdefined(net) || error("Run Example 1, the section \"Solve the reference power flow\", the section \"Derive a noisy measurement set\", the section \"Run the estimator\", Example 2, Example 7, Example 8, Example 9, Example 12, Example 13, Example 14 first: they set up net.")
net17 = deepcopy(net)
runpf!(net17, 40, 1e-10, 0)
meas17 = generateMeasurementsFromPF(net17; noise = true, rng = MersenneTwister(17))
part17 = area_partition(net17; areas = Dict("B1" => 1, "B2" => 1, "B7" => 1, "B3" => 2, "B4" => 2, "B5" => 2, "B6" => 2))
println("tie branches: ", join(["$(getCompName(net17.nodeVec[net17.branchVec[k].fromBus].comp)) - $(getCompName(net17.nodeVec[net17.branchVec[k].toBus].comp))" for k in part17.tie_branches], ", "))
r17 = runse_multi_area!(net17, meas17; partition = part17)
@printf("converged %s after %d passes, J = %.3f, dof = %d\n", r17.converged, r17.passes, r17.objectiveJ, r17.dof)
@printf("largest difference to the integrated estimate: %.1e pu, %.1e deg\n", r17.comparison.max_dvm_pu, r17.comparison.max_dva_deg)
@assert r17.converged && length(part17.tie_branches) == 3                  #src
@assert r17.comparison.max_dvm_pu < 1e-5 && r17.comparison.max_dva_deg < 1e-3  #src

# Reading aid (Example 17): three lines cross the cut (B2-B3, B6-B7, B2-B5),
# so B2, B3, B5, B6 and B7 are border buses; B1 and B4 are inside. Each pass, both
# areas estimate their interior with the border buses held, and the
# coordinator corrects the border buses from what both areas report. After
# a few passes nothing moves and the estimate equals the integrated one
# within the solver tolerance. `runse_multi_area_diagnostics` does the same
# for bad data: a gross error on a tie branch is found and removed in the
# same order as by `runse_diagnostics`.
#
# **Example 17, step by step on a shipped case.** `sp_case10_areas` has ten
# buses in three areas, stored in the case file. Area 1 holds the slack, area
# 2 two busbar sections on a closed coupler, area 3 a transformer whose ratio
# changer is released for estimation. The case file has the changer at step
# 0; the measurement table was taken with it at step +2, as if the position
# in the model were stale.
#
# Step 1: the case, its measurement table and the partition from the bus
# areas.

case17 = joinpath(pkgdir(Sparlectra), "data", "scf", "sp_case10_areas.scf.json")
table17 = joinpath(pkgdir(Sparlectra), "data", "scf", "sp_case10_areas.measurements.csv")
net10 = importSCF(case17)
readMeasurementTable!(net10; file = table17)
meas10 = Measurement[m for m in net10.measurements]
part10 = area_partition(net10)
busname(n, b) = getCompName(n.nodeVec[b].comp)
for a in part10.areas
  own = [busname(net10, b) for b in eachindex(net10.nodeVec) if part10.area_of_bus[b] == a]
  border = [busname(net10, b) for b in part10.border_buses if part10.area_of_bus[b] == a]
  println("area ", a, ": ", join(own, ", "), "; border buses ", join(border, ", "))
end
println("tie lines: ", join([getCompName(net10.branchVec[k].comp) for k in part10.tie_branches], ", "))
trafo10 = only(k for (k, br) in enumerate(net10.branchVec) if getCompName(br.comp) == "Oberried_110_20")
println("Oberried_110_20 released as ", net10.branchVec[trafo10].tap_est_mode, ", ", length(meas10), " measurement rows")
@assert part10.areas == [1, 2, 3] && length(part10.tie_branches) == 3  #src
@assert net10.branchVec[trafo10].tap_est_mode === :ratio              #src

# Reading aid (step 1): the three tie lines close a ring over the areas, so
# every area has two border buses. `Marbach_110b` sits inside area 2 although
# the coupler joins it to the border bus `Marbach_110`: a coupled group is
# never cut. Every row has one owner, a flow row the area of its measured
# end.
#
# Step 2: the reference, `runse!` on all rows at once.

## uses: meas10, net10 (Example 17)
(@isdefined(meas10) && @isdefined(net10)) || error("Run Example 17 first: it sets up meas10, net10.")
ri10 = runse!(deepcopy(net10), meas10)
tap(rows) = only(t for t in rows if t.name == "Oberried_110_20")
@printf("runse!: converged %s, J = %.3f, dof = %d, changer estimate %.2f, fixed at step %.0f\n", ri10.converged, ri10.objectiveJ, ri10.dof, tap(ri10.tapEstimates).s_est, tap(ri10.tapEstimates).step_fixed)
@assert ri10.converged && tap(ri10.tapEstimates).step_fixed == 2  #src

# Step 3: the multi-area estimate. Each area sees only its own rows.

## uses: meas10, net10, tap (Example 17)
(@isdefined(meas10) && @isdefined(net10) && @isdefined(tap)) || error("Run Example 17 first: it sets up meas10, net10, tap.")
r10 = runse_multi_area!(deepcopy(net10), meas10)
for (k, d) in enumerate(r10.border_change)
  @printf("pass %d: largest border change %.1e\n", k, d)
end
for a in r10.areas
  @printf("area %d: %d own rows, J %.2f\n", a.area, a.rows, a.j)
end
@printf("converged %s after %d passes, J = %.3f, dof = %d, changer fixed at step %.0f\n", r10.converged, r10.passes, r10.objectiveJ, r10.dof, tap(r10.tapEstimates).step_fixed)
@assert r10.converged && tap(r10.tapEstimates).step_fixed == 2  #src

# Reading aid (step 3): in every pass each area estimates its interior with
# all its border buses held, and returns what its rows say about the border:
# the Schur complement of its gain matrix and its reduced gradient. The
# coordinator adds them up and takes one Gauss-Newton step on the six border
# buses. The first passes run with the changer as a continuous state; once
# the border stops moving, area 3 rounds the changer to the nearest step
# (+2), and the jump in the border change after that is the effect of the
# rounding. Two more passes settle the border with the step fixed. The area
# values of J add up to the J of the whole set.
#
# Step 4: the comparison with `runse!`, per border bus and for the changer.

## uses: part10, r10, ri10 (Example 17)
(@isdefined(part10) && @isdefined(r10) && @isdefined(ri10)) || error("Run Example 17 first: it sets up part10, r10, ri10.")
cmp10 = r10.comparison
for row in cmp10.rows
  row.bus in part10.border_buses || continue
  @printf("%-14s %s  multi-area %10.6f  runse! %10.6f  difference %.1e\n", row.name, row.quantity, row.multi_area, row.integrated, row.difference)
end
@printf("largest difference: %.1e pu, %.1e deg; J difference %.1e\n", cmp10.max_dvm_pu, cmp10.max_dva_deg, r10.objectiveJ - ri10.objectiveJ)
println("tap release as in runse!: ", cmp10.tap_release_differences == 0)
@assert cmp10.max_dvm_pu < 1e-6 && cmp10.max_dva_deg < 1e-4                  #src
@assert isapprox(r10.objectiveJ, ri10.objectiveJ; rtol = 1e-6) && r10.dof == ri10.dof  #src

# Reading aid (step 4): the two estimates differ only in the last digits the
# solver tolerance leaves open. That is the point of the method: the areas
# never exchanged a single row, and the result is the integrated estimate
# on the same set, not an approximation of it.
#
# ## Where to go next
#
# - The sequel notebook,
#   [bad data and parameter estimation](https://github.com/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_se_diagnostics.ipynb):
#   what happens when $J/\mathrm{dof}$ is NOT near one: normalized
#   residuals and localizability, sequential elimination, robust
#   estimation, current measurements, and shunt-parameter estimation.
# - [`state_estimation_observable_islands.jl`](https://github.com/Welthulk/Sparlectra.jl/blob/main/examples/state_estimation/state_estimation_observable_islands.jl):
#   Examples 13 to 16 on sp_case118 and sp_case1354, with the run time per
#   step.
# - [`state_estimation_multi_area.jl`](https://github.com/Welthulk/Sparlectra.jl/blob/main/examples/state_estimation/state_estimation_multi_area.jl):
#   Example 17 on sp_case118 with its three MATPOWER areas, with bad data on
#   a tie branch.
# - New to Sparlectra? The
#   [workshop tour](https://github.com/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_tour.ipynb)
#   builds a network from scratch step by step.
# - [Slack types and short-circuit currents](https://github.com/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_slack_short_circuit.ipynb):
#   the grid-connection notebook: slack representations and IEC 60909-0
#   fault currents.
# - [State Estimation](https://welthulk.github.io/Sparlectra.jl/state_estimation/):
#   the full documentation: observability analysis, PMU angle measurements,
#   bad-data handling, and the measurement model.
# - [State-Estimation Configuration](https://welthulk.github.io/Sparlectra.jl/state_estimation_configuration/):
#   every `state_estimation.*` configuration key.
