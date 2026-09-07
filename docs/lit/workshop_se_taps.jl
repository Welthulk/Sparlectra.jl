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
# file: docs/lit/workshop_se_taps.jl                                          #src
# purpose: Literate.jl source of the third state-estimation workshop and its  #src
#          Documenter page: what J means next to its degrees of freedom, tap  #src
#          estimation (ratio, phase, cascade, and the case that is not        #src
#          identifiable), current measurements, PMU phasors, the topology     #src
#          advisory, island-wise estimation, and J versus J_active.           #src
#          Regenerate the committed outputs with                              #src
#          `julia --project=docs docs/generate_notebooks.jl`.                 #src

# # Taps, phasors and topology
#
# > **Level: Expert**, third part of the state-estimation series after the
# > [basics notebook](https://welthulk.github.io/Sparlectra.jl/generated/workshop_state_estimation/)
# > and [bad data and parameter estimation](https://welthulk.github.io/Sparlectra.jl/generated/workshop_se_diagnostics/).
#
# [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_se_taps.ipynb)
#
# > **Note:** This workshop was created with AI assistance and is reviewed
# > and curated by the maintainer; it is not a fully machine-generated text.
#
# The first two chapters treated the network model as given and asked what
# the measurements say about the STATE. This one asks the harder question:
# what if the model itself is wrong? A transformer that is three steps away
# from where the database thinks it is, a breaker whose status never got
# updated, a station behind a transformer that nobody measures. Each of
# those shows up as a large $J$, and each needs a different answer.
#
# Everything below runs on ONE study network. Each section adds a
# capability to the same case, so the numbers stay comparable: you can see
# $J$, the degrees of freedom, and the error against the truth move as the
# measurement set and the released parameters change. The chapter closes
# with the decision tree that ties the sections together: **$J$ is too
# large, what now?**
#
# What you will see, in order:
#
# 1. why $J$ alone means nothing without its degrees of freedom
# 2. what each added measurement actually buys
# 3. estimating an on-load tap changer, a phase shifter, and both at once
# 4. the tap that CANNOT be estimated, and why the estimator freezes it
# 5. what current magnitudes can and cannot do
# 6. PMU phasors and the reference offset the estimator solves for
# 7. the topology advisory: when the model, not the data, is wrong
# 8. island-wise estimation with one reference per island
# 9. $J$ versus $J_\mathrm{active}$ when rows are suppressed

#nb # ## Setup (Colab)
#nb # This cell installs Sparlectra from GitHub (branch `main`) into a fresh
#nb # temporary environment. The isolation matters on Colab: the shared
#nb # default environment ships many preinstalled packages, and installing
#nb # anything there triggers precompilation of that whole stack. Run this
#nb # cell first, once per session; it takes a few minutes.
#nb using Pkg
#nb Pkg.activate(temp = true)
#nb Pkg.add(url = "https://github.com/Welthulk/Sparlectra.jl", rev = "main")
#nb ## To test another branch, set rev to its name, e.g. rev = "dev/r0.10.0".
#nb ## For the latest registered release use: Pkg.add("Sparlectra")
#nb ## Switching versions in a running session? A "[loaded: ...]" note means
#nb ## the old version is still active; restart the runtime, then rerun
#nb ## this cell.

# ## Load the packages
#
# `Random` (standard library) seeds the synthetic measurement noise, and
# `Printf` formats the comparison tables.

using Sparlectra
using Random
using Printf

# ## Warm-up
#
# Julia compiles each function on first use. This cell warms the paths the
# notebook exercises on a tiny throwaway network, so the real study runs at
# full speed.

wnet = Net(name = "warmup", baseMVA = 100.0)
addBus!(net = wnet, busName = "A", vn_kV = 110.0)
addBus!(net = wnet, busName = "B", vn_kV = 110.0)
addProsumer!(net = wnet, busName = "A", type = "EXTERNALNETWORKINJECTION", referencePri = "A", vm_pu = 1.0, va_deg = 0.0)
addProsumer!(net = wnet, busName = "B", type = "ENERGYCONSUMER", p = 10.0, q = 3.0)
addPIModelACLine!(net = wnet, fromBus = "A", toBus = "B", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
t_pf = @elapsed runpf!(wnet, 10, 1e-8, 0)
setMeasurementsFromPF!(wnet; includeVm = true, includePinj = true, includeQinj = true, includePflow = true, includeQflow = true, noise = false)
t_se = @elapsed runse!(wnet)
println("warm: power flow ", round(t_pf; digits = 2), " s, estimator ", round(t_se; digits = 2), " s (first calls compile)")

# ## The study network
#
# A substation fed at 110 kV, feeding a 20 kV bus bar, with a second path
# that closes the loop. The loop matters: **a transformer tap is only
# observable inside a loop or with a voltage measurement behind it**, and
# this chapter is largely about that sentence.
#
# ```text
#   H1 ============== H2         110 kV (coupler line)
#    |                 |
#   T1 (OLTC + PST)   T2 (fixed)
#    |                 |
#   L1 ------------- L2 --- T3 --- L3    20 kV
#                             (radial, unmeasured behind)
# ```
#
# * **T1** carries both changers on one winding: a ratio tap (0.625 % per
#   step, band ±16) and a phase tap (1.25° per step, band ±20). That is the
#   cascade this chapter estimates.
# * **T2** is a plain transformer. It closes the loop that makes T1's tap
#   observable.
# * **T3** feeds a radial station, and later on nothing will measure the
#   voltage behind it. That is the case the estimator refuses to guess.
#
# `applyTapNameplate!` is the one path that declares tap-changer data:
# step size, position band, and the position the changer currently stands
# at. Importers use the same function, so a hand-built network and an
# imported one behave identically.

function build_net(; ratio_step = 0.0, phase_step = 0.0)
  net = Net(name = "workshop_se_taps", baseMVA = 100.0)
  for (b, vn) in (("H1", 110.0), ("H2", 110.0), ("L1", 20.0), ("L2", 20.0), ("L3", 20.0))
    addBus!(net = net, busName = b, vn_kV = vn)
  end
  addProsumer!(net = net, busName = "H1", type = "EXTERNALNETWORKINJECTION", referencePri = "H1", vm_pu = 1.02, va_deg = 0.0)
  addProsumer!(net = net, busName = "L1", type = "ENERGYCONSUMER", p = 25.0, q = 8.0)
  addProsumer!(net = net, busName = "L2", type = "ENERGYCONSUMER", p = 15.0, q = 5.0)
  addProsumer!(net = net, busName = "L3", type = "ENERGYCONSUMER", p = 8.0, q = 3.0)
  addPIModelACLine!(net = net, fromBus = "H1", toBus = "H2", r_pu = 0.010, x_pu = 0.080, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "L1", toBus = "L2", r_pu = 0.020, x_pu = 0.100, b_pu = 0.0, status = 1)
  addPIModelTrafo!(net = net, fromBus = "H1", toBus = "L1", r_pu = 0.002, x_pu = 0.060, b_pu = 0.0, ratio = 1.0, shift_deg = 0.0, status = 1)
  addPIModelTrafo!(net = net, fromBus = "H2", toBus = "L2", r_pu = 0.002, x_pu = 0.060, b_pu = 0.0, ratio = 1.0, shift_deg = 0.0, status = 1)
  addPIModelTrafo!(net = net, fromBus = "L2", toBus = "L3", r_pu = 0.004, x_pu = 0.070, b_pu = 0.0, ratio = 1.0, shift_deg = 0.0, status = 1)
  ## T1 (branch 3) carries BOTH changers; `*_current_step` is where they stand
  applyTapNameplate!(
    net.branchVec[3];
    tap_step = 0.00625, tap_min_step = -16.0, tap_max_step = 16.0, tap_current_step = ratio_step,
    phase_step_deg = 1.25, phase_min_step = -20.0, phase_max_step = 20.0, phase_current_step = phase_step,
    psi_deg = 90.0,
  )
  ## T3 (branch 5): the radial transformer, at neutral unless stated otherwise
  applyTapNameplate!(net.branchVec[5]; tap_step = 0.0125, tap_min_step = -8.0, tap_max_step = 8.0)
  ok, msg = validate!(net = net)
  ok || error("Network validation failed: $msg")
  return net
end

# The truth is a solved power flow the estimator never sees; the
# measurement set is derived from it with seeded Gaussian noise of exactly
# the declared standard deviations. `ratio_step` and `phase_step` put the
# changers of the TRUE network where the estimator will have to find them.

function truth_measurements(; ratio_step = 0.0, phase_step = 0.0, seed = 7, kwargs...)
  truth = build_net(; ratio_step = ratio_step, phase_step = phase_step)
  _, status = runpf!(truth, 40, 1e-10, 0)
  status == 0 || error("truth power flow did not converge")
  calcNetLosses!(truth)
  std = measurementStdDevs(vm = 1e-3, pinj = 1.0, qinj = 1.0, pflow = 0.7, qflow = 0.7)
  meas = generateMeasurementsFromPF(
    truth;
    includeVm = true, includePinj = true, includeQinj = true, includePflow = true, includeQflow = true,
    noise = true, stddev = std, rng = MersenneTwister(seed), kwargs...,
  )
  return truth, meas
end

truth, full_set = truth_measurements(seed = 3)
vm_true = [getNodeVm(n) for n in truth.nodeVec]
println("study network: ", length(truth.nodeVec), " buses, ", length(truth.branchVec), " branches, ", length(full_set), " measurements")

# ## A number that proves nothing: J at zero degrees of freedom
#
# **Example 1.** The state of an $n$-bus network has $2n-1$ unknowns: one
# magnitude per bus, and one angle per bus except the reference. Five buses
# means nine states. Feed the estimator exactly nine measurements and it
# will reproduce them exactly, because it CAN: nine equations, nine
# unknowns, no room left for a residual.
#
# The values below are deliberately biased (every voltage is 4 mV too high,
# every injection 0.9 MW off). Watch $J$ anyway.

names = ["H1", "H2", "L1", "L2", "L3"]
minimal = Measurement[]
for (i, b) in enumerate(names)
  addVmMeasurement!(minimal; net = truth, busName = b, value = getNodeVm(truth.nodeVec[i]) + 0.004, sigma = 1e-3)
end
for (i, b) in enumerate(names)
  i == 1 && continue                          ## the slack injection is free
  load_p = truth.nodeVec[i]._pƩLoad === nothing ? 0.0 : truth.nodeVec[i]._pƩLoad
  addPinjMeasurement!(minimal; net = truth, busName = b, value = -load_p + 0.9, sigma = 1.0)
end

res_min = runse!(build_net(), minimal; maxIte = 40, tol = 1e-10)
@printf("exactly determined: m = %d, states = %d, dof = %d, J = %.6f\n",
        length(minimal), 2 * length(truth.nodeVec) - 1, res_min.dof, res_min.objectiveJ)

# $J = 0$ with garbage input. This is the single most misread number in
# state estimation: **$J$ measures the CONTRADICTION between measurements,
# and without redundancy there is nothing to contradict.** A set with
# `dof = 0` always fits perfectly, no matter how wrong it is.
#
# Adding the reactive injections gives the estimator something to compare:

redundant = deepcopy(minimal)
for (i, b) in enumerate(names)
  i == 1 && continue
  load_q = truth.nodeVec[i]._qƩLoad === nothing ? 0.0 : truth.nodeVec[i]._qƩLoad
  addQinjMeasurement!(redundant; net = truth, busName = b, value = -load_q + 0.7, sigma = 1.0)
end
res_red = runse!(build_net(), redundant; maxIte = 40, tol = 1e-10)
@printf("redundant:          m = %d, dof = %d, J = %.3f, J/dof = %.3f\n",
        length(redundant), res_red.dof, res_red.objectiveJ, res_red.objectiveJ / res_red.dof)

# Reading aid: only now is $J$ a statement about the data. The rule of
# thumb $J/\mathrm{dof} \approx 1$ needs a denominator, and the denominator
# is redundancy you have to pay for with measurements.

# ## What each measurement buys
#
# **Example 2.** Starting from the exactly determined core, one row is
# added at a time. Two columns matter: `J/dof`, the statistical health, and
# `max_err`, the largest deviation of the estimated voltage from the truth
# (a number only this notebook can compute, because only here the truth is
# known).

core = filter(m -> startswith(m.id, "Vm_") || startswith(m.id, "Pinj_"), deepcopy(full_set))
extra = filter(m -> !(startswith(m.id, "Vm_") || startswith(m.id, "Pinj_")), deepcopy(full_set))
growing = deepcopy(core)
for step in 0:6
  step > 0 && push!(growing, extra[step])
  net = build_net()
  res = runse!(net, deepcopy(growing); maxIte = 40, tol = 1e-10)
  err = maximum(abs.([getNodeVm(n) for n in net.nodeVec] .- vm_true))
  @printf("m = %2d  dof = %2d  J = %7.3f  J/dof = %6.3f  max_err = %.2e  added: %s\n",
          length(growing), res.dof, res.objectiveJ, res.dof > 0 ? res.objectiveJ / res.dof : NaN, err,
          step == 0 ? "(core)" : extra[step].id)
end

# Reading aid: the accuracy jump comes with the THIRD reactive injection,
# where the error against the truth drops by an order of magnitude. After
# that the additional rows buy statistics, not accuracy. Note also how
# wildly `J/dof` swings while `dof` is small: with three degrees of freedom
# a factor of two is normal, which is why the band test uses the chi-square
# distribution instead of comparing against 1.

# ## Example 3: the on-load tap changer
#
# Now the model is wrong on purpose. In the true network T1's ratio changer
# stands four steps up; the estimator's model believes it is at neutral.
# Nothing else changes.
#
# First, what that costs if the tap is NOT released:

truth_oltc, meas_oltc = truth_measurements(ratio_step = 4.0)
res_blind = runse!(build_net(), deepcopy(meas_oltc); maxIte = 40, tol = 1e-10)
@printf("tap not released: J = %.1f, dof = %d, inside band: %s\n",
        res_blind.objectiveJ, res_blind.dof, res_blind.jWithin3Sigma)

# A four-step error on one transformer inflates $J$ by a factor of fifty.
# The band test fires, the residuals point everywhere at once, and no
# amount of bad-data elimination will help: **the data is fine, the model
# is wrong.**
#
# `setTapEstimation!` releases the tap position as an additional state.
# Mode `:ratio` releases the ratio changer only:

model_oltc = build_net()
release = setTapEstimation!(model_oltc; trafo = 3, mode = :ratio)
res_oltc = runse!(model_oltc, deepcopy(meas_oltc); maxIte = 40, tol = 1e-10)
est = res_oltc.tapEstimates[1]
@printf("released: J = %.2f, dof = %d | electrical step = %+.3f -> fixed to %+d (truth %+d), out of range: %s\n",
        res_oltc.objectiveJ, res_oltc.dof, est.electrical_step_1, est.fixed_step_1, 4, est.out_of_range)

# The `runse!` info line above the result reports what the fixation cost:
# the free estimate fits marginally better than the grid position (it has
# one more degree of freedom to spend), and snapping to the mechanical step
# gives that degree of freedom back. A large jump there would be the signal
# that the data does NOT agree with any real tap position.
#
# Reading aid: the estimator finds 4.03 steps and FIXES that to the nearest
# mechanical position, +4. The fixation is the point. A tap changer is a
# discrete device; an estimate of "4.03 steps" is not a physical answer,
# and reporting it as one would invite an operator to chase a position that
# does not exist. The continuous value stays visible next to it, because
# its distance from the grid is the honest uncertainty statement: 0.03
# steps means the data agrees with position 4, and 0.5 would mean the data
# cannot tell 4 from 5.

# ## Example 4: the phase shifter
#
# The same transformer, the same procedure, the other changer. In the true
# network the phase tap stands three steps out; the model believes neutral.
# Mode `:pst` needs `alpha_deg`, the nameplate direction of the regulating
# winding. It is a specification, never estimated: a symmetric phase
# shifter injects its voltage at 90° to the line voltage, and the machine
# cannot rotate that angle.

truth_pst, meas_pst = truth_measurements(phase_step = 3.0, seed = 11)
model_pst = build_net()
setTapEstimation!(model_pst; trafo = 3, mode = :pst, alpha_deg = 90.0)
res_pst = runse!(model_pst, deepcopy(meas_pst); maxIte = 40, tol = 1e-10)
est_pst = res_pst.tapEstimates[1]
@printf("PST: J = %.2f | electrical step = %+.3f -> fixed to %+d (truth %+d)\n",
        res_pst.objectiveJ, est_pst.electrical_step_2, est_pst.fixed_step_2, 3)

# ## Example 5: both changers at once
#
# A transformer that regulates voltage AND angle carries two changers on
# one winding, and their effects multiply rather than add: the combined
# ratio is a cascade $t(r_1, r_2)$, not a sum of two independent knobs.
# Mode `:both` releases both positions as states of that cascade.
#
# The true network now has the ratio changer two steps DOWN and the phase
# changer two steps up, an ambiguous-looking combination the estimator has
# to separate:

truth_both, meas_both = truth_measurements(ratio_step = -2.0, phase_step = 2.0, seed = 13)
model_both = build_net()
setTapEstimation!(model_both; trafo = 3, mode = :both, alpha_deg = 90.0)
res_both = runse!(model_both, deepcopy(meas_both); maxIte = 40, tol = 1e-10)
est_both = res_both.tapEstimates[1]
@printf("cascade: J = %.2f | ratio %+.3f -> %+d (truth %+d), phase %+.3f -> %+d (truth %+d)\n",
        res_both.objectiveJ, est_both.electrical_step_1, est_both.fixed_step_1, -2,
        est_both.electrical_step_2, est_both.fixed_step_2, 2)

# Reading aid: both positions come back within 0.1 steps of the truth. They
# are separable because they act differently: the ratio changer moves the
# voltage magnitude across the transformer, the phase changer moves the
# active flow through the loop. Take the loop away (open T2's path) and the
# separation collapses, which is the subject of the next example.

# ## Example 6: the tap that cannot be estimated
#
# T3 feeds a radial station. Its tap changes the voltage behind it, and
# nothing else. If no measurement sits behind that transformer, the tap and
# the unknown downstream voltage are the SAME unknown: any tap position can
# be compensated by a voltage, and the estimator would happily return one
# of infinitely many solutions.
#
# Here the true T3 stands three steps up, and the voltage measurement at L3
# is removed from the set (the station has no voltage transducer):

truth_bridge = build_net()
applyTapNameplate!(truth_bridge.branchVec[5]; tap_step = 0.0125, tap_min_step = -8.0, tap_max_step = 8.0, tap_current_step = 3.0)
runpf!(truth_bridge, 40, 1e-10, 0)
calcNetLosses!(truth_bridge)
std_bridge = measurementStdDevs(vm = 1e-3, pinj = 1.0, qinj = 1.0, pflow = 0.7, qflow = 0.7)
meas_bridge = generateMeasurementsFromPF(
  truth_bridge;
  includeVm = true, includePinj = true, includeQinj = true, includePflow = true, includeQflow = true,
  noise = true, stddev = std_bridge, rng = MersenneTwister(5),
)
filter!(m -> m.id != "Vm_bus_5", meas_bridge)   ## L3 is bus 5: no voltage behind T3

model_bridge = build_net()
setTapEstimation!(model_bridge; trafo = 5, mode = :ratio)
res_bridge = runse!(model_bridge, meas_bridge; maxIte = 40, tol = 1e-10)
est_bridge = res_bridge.tapEstimates[1]
@printf("radial T3: J = %.2f | electrical step = %+.3f -> fixed to %+d (truth %+d)\n",
        res_bridge.objectiveJ, est_bridge.electrical_step_1, est_bridge.fixed_step_1, 3)

# Reading aid: the estimator does NOT return 3, and it does not return a
# random number either. It FREEZES the released tap at its declared
# position and says so in a warning (a bridge transformer without a voltage
# measurement on its cut-off side, where the tap would absorb the
# downstream voltage). $J$ stays healthy, because nothing in the data
# contradicts the model: the missing measurement is exactly the information
# that would have exposed the error.
#
# This is the most important slide of the chapter. **A released tap that
# comes back at its declared position is not a confirmation.** Check
# whether it was estimated at all, and whether the loop or the voltage
# measurement that makes it observable actually exists.

# ## Example 7: what current magnitudes can and cannot do
#
# A current magnitude is the cheapest measurement in the substation and the
# most common one in older telemetry. It also carries no direction and
# linearizes badly near zero. Sparlectra treats it as an AUXILIARY
# measurement, and the observability check says so outright: it evaluates
# the system WITHOUT current magnitudes, because a system that is
# observable only through them is not observable in the sense the estimator
# needs.
#
# First, what they add to a complete set:

_, with_currents = truth_measurements(seed = 3, includeImag = true)
net_plain = build_net(); res_plain = runse!(net_plain, deepcopy(full_set); maxIte = 40, tol = 1e-10)
net_curr = build_net(); res_curr = runse!(net_curr, deepcopy(with_currents); maxIte = 40, tol = 1e-10)
dvm = maximum(abs.([getNodeVm(n) for n in net_plain.nodeVec] .- [getNodeVm(n) for n in net_curr.nodeVec]))
@printf("complete set:  m = %d -> %d, J/dof = %.3f -> %.3f, state moves by %.1e pu\n",
        length(full_set), length(with_currents), res_plain.objectiveJ / res_plain.dof,
        res_curr.objectiveJ / res_curr.dof, dvm)

# Ten extra rows, and the state barely moves. Where the power flows are
# already measured, currents are redundancy, not information.
#
# Now take the flow measurements off one line and see what is left. The
# sparse set has no P/Q on the 20 kV tie; the second set replaces them with
# current magnitudes on the same bay:

sparse_set = filter(m -> !occursin("branch_2", m.id), deepcopy(full_set))
net_sparse = build_net(); append!(net_sparse.measurements, deepcopy(sparse_set))
obs_sparse = evaluate_global_observability(net_sparse; flatstart = true, jacEps = 1e-6)
res_sparse = runse!(net_sparse; maxIte = 40, tol = 1e-10)
err_sparse = maximum(abs.([getNodeVm(n) for n in net_sparse.nodeVec] .- vm_true))

sparse_curr = filter(m -> !occursin("branch_2", m.id) || startswith(m.id, "Imag"), deepcopy(with_currents))
net_sc = build_net(); append!(net_sc.measurements, deepcopy(sparse_curr))
obs_sc = evaluate_global_observability(net_sc; flatstart = true, jacEps = 1e-6)
res_sc = runse!(net_sc; maxIte = 40, tol = 1e-10)
err_sc = maximum(abs.([getNodeVm(n) for n in net_sc.nodeVec] .- vm_true))

@printf("no flows on the tie:      m = %d, observability = %s, dof = %d, max_err = %.2e\n",
        length(sparse_set), obs_sparse.quality, res_sparse.dof, err_sparse)
@printf("currents instead:         m = %d, observability = %s, dof = %d, max_err = %.2e\n",
        length(sparse_curr), obs_sc.quality, res_sc.dof, err_sc)

# Reading aid, and this is the "for and against" in two numbers: the
# currents do NOT change the observability verdict (they are excluded from
# it by construction), but they nearly halve the error against the truth.
# That is the honest case for cheap current telemetry: it improves an
# estimate that already stands on its own, and it cannot carry one that
# does not.
#
# The missing direction is what a PMU adds. `addCurrentPhasorMeasurement!`
# contributes the magnitude AND the angle of a branch current, and the
# angle turns an amount into a flow with a sign.

# ## Example 8: PMU phasors and the reference offset
#
# A PMU reports voltage angles against ITS time base, and that base is not
# the estimator's slack reference. The difference is a constant offset over
# all PMU channels. Subtracting a guess would bias every angle; ignoring it
# would corrupt the state. Sparlectra solves for it: with `pmu_ref_offset:
# auto` (the default) the offset becomes an extra state whenever voltage
# angles are present.
#
# The set below carries PMU angles at H2 and L1, deliberately shifted by
# 4 degrees against the slack reference:

_, pmu_set = truth_measurements(
  seed = 21, includeVa = true, vaBusIdxs = [2, 3], vaRefOffsetDeg = 4.0,
)
res_pmu = runse!(build_net(), pmu_set; maxIte = 40, tol = 1e-10)
@printf("PMU set: m = %d, J = %.2f, dof = %d, estimated reference offset = %+.3f deg (injected %+.1f)\n",
        length(pmu_set), res_pmu.objectiveJ, res_pmu.dof,
        res_pmu.vaRefOffsetDeg === nothing ? NaN : res_pmu.vaRefOffsetDeg, 4.0)

# Reading aid: the offset comes back as 4.0°, and the state is unaffected
# by it. The cost is one degree of freedom, which is the cheapest possible
# price for not having to trust a time synchronisation you cannot verify.

# ## Example 9: when the model, not the data, is wrong
#
# A breaker was opened in the field and the status never reached the model.
# The estimator now believes in a line that carries nothing. The bay's
# transducers report what they see, about zero, on all four rows.

truth_topo = build_net()
for br in (truth_topo.branchVec[2],)                 ## the 20 kV tie is really open
  br.status = 0
  br.from_status = 0
  br.to_status = 0
end
runpf!(truth_topo, 40, 1e-10, 0)
calcNetLosses!(truth_topo)
std_topo = measurementStdDevs(vm = 1e-3, pinj = 1.0, qinj = 1.0, pflow = 0.7, qflow = 0.7)
meas_topo = generateMeasurementsFromPF(
  truth_topo;
  includeVm = true, includePinj = true, includeQinj = true, includePflow = true, includeQflow = true,
  noise = true, stddev = std_topo, rng = MersenneTwister(31),
)
model_topo = build_net()                             ## the model still has it closed
for (dir, val) in ((:from, 0.05), (:to, -0.04))
  addPflowMeasurement!(meas_topo; net = model_topo, fromBus = "L1", toBus = "L2", value = val, sigma = 0.7, direction = dir)
  addQflowMeasurement!(meas_topo; net = model_topo, fromBus = "L1", toBus = "L2", value = val / 3, sigma = 0.7, direction = dir)
end

pre = validate_topology(model_topo, meas_topo)
for f in pre.findings
  println(f.kind, " at ", f.location, " (", f.severity, "): ", f.evidence)
end
res_topo = runse!(model_topo, meas_topo; maxIte = 40, tol = 1e-10)
@printf("wrong status:  J = %.1f, dof = %d, inside band: %s\n", res_topo.objectiveJ, res_topo.dof, res_topo.jWithin3Sigma)

model_fixed = build_net()
for br in (model_fixed.branchVec[2],)
  br.status = 0
  br.from_status = 0
  br.to_status = 0
end
res_fixed = runse!(model_fixed, deepcopy(meas_topo); maxIte = 40, tol = 1e-10)
@printf("status corrected: J = %.1f, dof = %d, inside band: %s\n", res_fixed.objectiveJ, res_fixed.dof, res_fixed.jWithin3Sigma)

# Reading aid: the pre-check names the element and the evidence BEFORE the
# estimation runs, and the estimation still runs. That is deliberate: the
# findings are advisory, and **Sparlectra never switches an element by
# itself.** A topology error looks like bad data in the residuals, so a
# sequential elimination would happily throw away four healthy
# measurements, one after another, and never fix the cause. The correction
# belongs to the operator, and the evidence is there to support it.

# ## Example 10: island-wise estimation
#
# Open the coupler and the tie, and the substation becomes two networks
# with one feed each. There is no such thing as a common angle reference
# across them: the estimator partitions the measurement set by AC island
# and gives each island its own reference.

iso = build_net()
addProsumer!(net = iso, busName = "H2", type = "EXTERNALNETWORKINJECTION", referencePri = "H2", vm_pu = 1.015, va_deg = 0.0)
for b in (1, 2)                                      ## coupler and tie open
  iso.branchVec[b].status = 0
  iso.branchVec[b].from_status = 0
  iso.branchVec[b].to_status = 0
end
validate!(net = iso)
runpf!(iso, 40, 1e-10, 0; islands_enabled = true)
calcNetLosses!(iso)
std_iso = measurementStdDevs(vm = 1e-3, pinj = 1.0, qinj = 1.0, pflow = 0.7, qflow = 0.7)
meas_iso = generateMeasurementsFromPF(
  iso;
  includeVm = true, includePinj = true, includeQinj = true, includePflow = true, includeQflow = true,
  noise = true, stddev = std_iso, rng = MersenneTwister(41),
)

model_iso = build_net()
addProsumer!(net = model_iso, busName = "H2", type = "EXTERNALNETWORKINJECTION", referencePri = "H2", vm_pu = 1.015, va_deg = 0.0)
for b in (1, 2)
  model_iso.branchVec[b].status = 0
  model_iso.branchVec[b].from_status = 0
  model_iso.branchVec[b].to_status = 0
end
validate!(net = model_iso)
res_iso = runse!(model_iso, meas_iso; maxIte = 40, tol = 1e-10)
@printf("total: J = %.2f, dof = %d, islands = %d\n", res_iso.objectiveJ, res_iso.dof, length(res_iso.islands))
for isl in res_iso.islands
  @printf("  island %d: %d buses, converged = %s, J = %.2f, dof = %d, band = %s\n",
          isl.island, isl.n_bus, isl.converged, isl.objectiveJ, isl.dof, isl.band_reason)
end

# Reading aid: the totals are sums over the islands, and each island
# carries its own verdict. That matters when one island is well measured
# and the other is not: a single global $J$ would average the two and hide
# the weak one.

# ## Example 11: J and J_active
#
# The last number of the chapter. One flow measurement is off by 25 MW, far
# outside anything noise explains.

bad_set = deepcopy(full_set)
k = findfirst(m -> m.id == "Pflow_branch_1_from", bad_set)
bad_set[k] = Measurement(
  typ = bad_set[k].typ, value = bad_set[k].value + 25.0, sigma = bad_set[k].sigma,
  branchIdx = bad_set[k].branchIdx, direction = bad_set[k].direction, id = bad_set[k].id,
)
res_bad = runse!(build_net(), deepcopy(bad_set); maxIte = 40, tol = 1e-10)
@printf("plain WLS:    J = %.1f, dof = %d, inside band: %s\n", res_bad.objectiveJ, res_bad.dof, res_bad.jWithin3Sigma)

# In replacement mode the estimator keeps the row in the set but pins it to
# a large sigma during the solve, the way an EMS suppression list does. The
# statistics stay on the ORIGINAL sigma, so the alarm does not go away:

se_cfg = StateEstimationConfig(robust_mode = :replacement, k_suppress = 4.0, suppression_sigma = 2000.0)
res_repl = runse!(build_net(), deepcopy(bad_set), se_cfg)
act = res_repl.activeObjective
@printf("replacement:  J = %.1f (honest), J_active = %.2f over dof = %d, suppressed rows = %d\n",
        res_repl.objectiveJ, act.j, act.dof, act.suppressed)

# Note the warning above both runs: the topology pre-check flags a KCL
# violation at H1, because a 25 MW error on one branch flow does not
# balance against the injections either. A gross measurement error and a
# wrong breaker status look similar from the outside, which is precisely
# why the pre-check reports and never switches.
#
# Reading aid: two numbers, two jobs. $J$ answers "is my data healthy?"
# and must stay loud, because a suppressed row is still a broken
# transducer somebody has to fix. $J_\mathrm{active}$ answers "how well
# does the state fit what I actually trusted?" and is the number to judge
# the ESTIMATE by. Reporting only the second one would turn a suppression
# list into a way of making alarms disappear, which is exactly how EMS
# suppression lists go bad in practice.

# ## J is too large. What now?
#
# The sections above are the four branches of one decision:
#
# | Symptom | Likely cause | What to do |
# |---|---|---|
# | A few rows carry large normalized residuals, the rest are clean | **Measurement error** | Bad-data localization and elimination, see the [diagnostics chapter](https://welthulk.github.io/Sparlectra.jl/generated/workshop_se_diagnostics/) |
# | Residuals are large across a whole area, elimination removes healthy rows without curing $J$ | **Model error**: a tap or a shunt | Release the parameter (`setTapEstimation!`, `setShuntEstimation!`) and check whether $J$ drops |
# | Flow rows near zero on an element the model believes is in service, or flow on one it believes is open | **Topology error** | Read `validate_topology`, correct the status, never let the estimator switch it |
# | $J$ is small or zero and everything looks perfect | **Too few measurements** | Check `dof`: at zero there is nothing to check, and a released tap may have been frozen |
#
# The order matters. Topology first (it corrupts everything downstream),
# then model parameters, then individual measurements. Eliminating rows
# before checking the model is the most common way to arrive at a
# beautifully converged estimate of a network that does not exist.

# ## Where to go next
#
# - Every knob used here, with its defaults and its interactions:
#   [State-Estimation Configuration](https://welthulk.github.io/Sparlectra.jl/state_estimation_configuration/).
# - The measurement model itself: which types exist, what they mean, and
#   how sets travel as CSV files:
#   [State Estimation](https://welthulk.github.io/Sparlectra.jl/state_estimation/).
# - Bad-data localization, robust estimation and shunt estimation in
#   depth: the [diagnostics notebook](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_se_diagnostics.ipynb).
# - New to the estimator? The
#   [basics notebook](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_state_estimation.ipynb)
#   walks from network build to estimated state.
