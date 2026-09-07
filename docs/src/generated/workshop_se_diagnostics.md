```@meta
EditURL = "../../lit/workshop_se_diagnostics.jl"
```

# Bad data and parameter estimation

> **Level: Expert**, sequel to the [state-estimation basics notebook](https://welthulk.github.io/Sparlectra.jl/generated/workshop_state_estimation/); this chapter assumes you know what $J$, `dof`, and observability mean.

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_se_diagnostics.ipynb)

> **Note:** This workshop was created with AI assistance and is reviewed
> and curated by the maintainer; it is not a fully machine-generated text.

The basics notebook ended with one number worth glancing at after every
run: $J/\mathrm{dof} \approx 1$ for healthy Gaussian noise. This chapter
is about the days when that number is NOT close to one. Real telemetry
contains gross errors: a stuck transducer, a sign flip in a RTU mapping,
a CT ratio entered wrong. The estimator then faces three questions, in
exactly this order:

1. **Detection**: is something wrong at all? (the $J$ band test)
2. **Localization**: WHICH measurement is wrong? (normalized residuals
   and the localizability measure $w_{ii}$)
3. **Treatment**: remove it (sequential elimination) or ride through it
   (robust weight modification)?

After the bad-data arc, the same machinery bends to a different purpose:
estimating a MODEL PARAMETER (a shunt reactor's susceptance) as an extra
state, which turns the estimator into a parameter-correction tool.

The study network is the 7-bus workshop ring from the basics chapter
(B1 carries the grid connection, the diagonals are the cross-ties B2-B5
and B3-B6):

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
> version may change over time; this notebook targets Julia >= 1.12.

## Load the packages

`Random` (standard library) seeds the synthetic measurement noise so
every run of this notebook produces the same numbers.

````@example workshop_se_diagnostics
using Sparlectra
using Random
````

## Warm-up

Julia compiles each function on first use. This cell warms the paths the
notebook exercises (power flow, estimator, diagnostics) on a tiny
throwaway network, so the real study runs at full speed.

````@example workshop_se_diagnostics
wnet = Net(name = "warmup", baseMVA = 100.0)
addBus!(net = wnet, busName = "A", vn_kV = 110.0)
addBus!(net = wnet, busName = "B", vn_kV = 110.0)
addProsumer!(net = wnet, busName = "A", type = "EXTERNALNETWORKINJECTION", referencePri = "A", vm_pu = 1.0, va_deg = 0.0)
addProsumer!(net = wnet, busName = "B", type = "ENERGYCONSUMER", p = 10.0, q = 3.0)
addPIModelACLine!(net = wnet, fromBus = "A", toBus = "B", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
t_pf = @elapsed runpf!(wnet, 10, 1e-8, 0)
setMeasurementsFromPF!(wnet; includeVm = true, includePinj = true, includeQinj = true, includePflow = true, includeQflow = true, noise = false)
t_dg = @elapsed validate_measurements(wnet)
println("warm: power flow ", round(t_pf; digits = 2), " s, diagnostics ", round(t_dg; digits = 2), " s (first calls compile)")
````

## The study network and a healthy measurement set

Same ring as in the basics chapter: seven 110 kV buses, ring plus two
chords, an external network injection at `B1` as slack, a generator at
`B3`, loads elsewhere. The solved power flow plays the role of the (in
reality unknown) true state; the measurement set is derived from it with
seeded Gaussian noise of exactly the declared standard deviations.

````@example workshop_se_diagnostics
net = Net(name = "workshop_se_diag_7bus", baseMVA = 100.0)

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

ite_pf, status_pf = runpf!(net, 40, 1e-10, 0)
status_pf == 0 || error("Power flow did not converge")
calcNetLosses!(net)   ## populates the branch flows Example 5 reads

# the truth, kept for later comparisons (the estimator never sees it)
vm_true = [net.nodeVec[i]._vm_pu for i in eachindex(net.nodeVec)]

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
println("reference PF: ", ite_pf, " iterations; ", length(net.measurements), " measurements created")
````

## Detection: the band test on a healthy set

**Example 1: what healthy looks like.** `validate_measurements` runs the
estimator and returns a machine-readable diagnostics report. The
detection layer is the $J$ band test: under the assumed noise model,
$J(\hat{x})$ is approximately chi-square distributed with
$\nu = \mathrm{dof}$ degrees of freedom. Instead of comparing $J$ to
$\nu \pm 3\sqrt{2\nu}$ directly (the raw chi-square is skewed for small
$\nu$, so a symmetric band is dishonest exactly when you need it most),
Sparlectra applies the **Wilson-Hilferty transformation**: the cube root
$(J/\nu)^{1/3}$ is very close to normally distributed even for small
$\nu$, and

```math
z_{WH} = \frac{(J/\nu)^{1/3} - \left(1 - \frac{2}{9\nu}\right)}{\sqrt{\frac{2}{9\nu}}}
```

is the honest z-score. The verdict carries a reason: `:high` means bad
data or a wrong model, `:low` means the residuals are SMALLER than the
declared sigmas promise (overly pessimistic sigmas, or a synthetic
noise-free set), and `small_redundancy` flags verdicts on thin evidence
($\nu < 30$).

````@example workshop_se_diagnostics
rep = validate_measurements(net)
println("converged:          ", rep.converged)
println("global consistency: ", rep.global_consistency)
println("J = ", round(rep.objective.value; digits = 2), ", dof = ", rep.objective.dof, ", J/dof = ", round(rep.objective.value / rep.objective.dof; digits = 3))
println("z_wh = ", round(rep.objective.z_wh; digits = 3), "  (|z| <= 3 is in band; reason = ", rep.objective.reason, ")")
println("suspicious measurements: ", length(rep.suspicious_measurements))
````

Reading aid (Example 1): $J/\mathrm{dof}$ near one, $|z_{WH}| \le 3$, no
suspicious rows: the declared sigmas and the observed scatter agree.
One implementation note for large systems: the report also carries
`omega_path`, which says how the residual covariance diagonal was
computed. On this small network it reads `:dense`; above
`state_estimation.takahashi_min_states` states (default 200) the
diagnostics switch to a sparse Takahashi selected inverse, same numbers,
a fraction of the cost.

````@example workshop_se_diagnostics
println("omega_path: ", rep.omega_path)
````

## Localization: normalized residuals and w_ii

**Example 2: one gross error, and the arithmetic that finds it.** We
corrupt a single active-power flow measurement by +25 MW, about 36 times
its declared sigma of 0.7 MW: the classic stuck-value error. Raw
residuals do NOT point at the culprit, because the WLS fit SPREADS the
error over its neighborhood: the estimator bends the state a little, so
neighboring measurements pick up part of the misfit. Two quantities undo
that spreading:

The **normalized residual** divides each residual by what a healthy
residual at that position would scatter,

```math
r_i^N = \frac{r_i}{\sqrt{\Omega_{ii}}}, \qquad \Omega = \mathrm{diag}(W)^{-1} - H G^{-1} H^{\top}
```

and the largest $|r^N|$ is the classical prime suspect (the
largest-normalized-residual test).

The **residual sensitivity** $w_{ii} = \Omega_{ii} \, w_i \in [0, 1]$
answers a subtler question: how much of an error in measurement $i$
shows up in measurement $i$'s OWN residual? At $w_{ii} = 0$ the row is
critical (its residual is blind, the basics chapter's $S_{ii} = 0$), at
$w_{ii} = 1$ the error lands fully in its own residual. The rule of
thumb from the literature: only rows with $w_{ii} > 0.3$ are
LOCALIZABLE; below that, a large $r^N$ may well be collateral damage
from an error elsewhere. `validate_measurements` reports both per row.

````@example workshop_se_diagnostics
bad_idx = findfirst(m -> m.typ == Sparlectra.PflowMeas, net.measurements)
m0 = net.measurements[bad_idx]
net.measurements[bad_idx] = Measurement(typ = m0.typ, value = m0.value + 25.0, sigma = m0.sigma, busIdx = m0.busIdx, branchIdx = m0.branchIdx, direction = m0.direction, id = m0.id, linkIdx = m0.linkIdx)
println("corrupted ", m0.id, ": ", round(m0.value; digits = 2), " -> ", round(m0.value + 25.0; digits = 2), " MW (sigma ", m0.sigma, ")")

rep_bad = validate_measurements(net)
println()
println("J/dof now: ", round(rep_bad.objective.value / rep_bad.objective.dof; digits = 1), ", z_wh = ", round(rep_bad.objective.z_wh; digits = 1), " (reason = ", rep_bad.objective.reason, ")")
top = rep_bad.largest_normalized_residual
println("prime suspect: ", top.id, "  |r_N| = ", round(top.abs_normalized_residual; digits = 1), ", w_ii = ", round(top.wii; digits = 2), ", localizable = ", top.localizable)
println("suspicious rows (|r_N| >= 3): ", length(rep_bad.suspicious_measurements))
for row in rep_bad.measurement_ranking[1:5]
  println("  ", rpad(row.id, 22), " |r_N| = ", lpad(round(row.abs_normalized_residual; digits = 1), 6), "   w_ii = ", round(row.wii; digits = 2))
end
````

Reading aid (Example 2): detection fired (`:high`, $z_{WH}$ far out of
band), and the corrupted row tops the ranking with a comfortable
localizability margin. Note the collateral damage in the list: several
HEALTHY rows also exceed the suspicion threshold because the bent state
misfits them. This is exactly why bad data is removed one row at a time,
not as a batch.

## Treatment 1: sequential elimination

**Example 3: identify, remove, re-run, repeat.** `runse_diagnostics`
automates the loop: while the band test fails AND a localizable
suspicious row exists AND fewer than `max_eliminations` rows are gone,
it deactivates the worst localizable row and re-runs the estimator. Rows
that encode model knowledge rather than telemetry (zero-injection
pseudo-measurements, derived shunt rows) are protected from elimination.
The trace records every step; the stop reason says why the loop ended.

````@example workshop_se_diagnostics
diag = runse_diagnostics(net; max_eliminations = 3)
println("eliminations: ", length(diag.eliminations), ", stop reason: ", diag.stop_reason)
for e in diag.eliminations
  println("  step ", e.elimination, ": removed ", e.id, " (|r_N| before = ", round(abs(e.normalized_residual_before); digits = 1), ")")
end
fin = diag.final_diagnostics
println("after elimination: J/dof = ", round(fin.objective.value / fin.objective.dof; digits = 3), ", z_wh = ", round(fin.objective.z_wh; digits = 3), ", reason = ", fin.objective.reason)
````

Reading aid (Example 3): one elimination, and it removed exactly the row
we corrupted; the band test returns in band and the loop stops. The
collateral suspects from Example 2 are healthy again without ever being
touched, which is the whole point of the sequential (rather than batch)
strategy.

## Treatment 2: robust estimation (ride through it)

**Example 4: keep the row, dull its influence.** Elimination is the
identification workflow: an operator wants to find and repair the broken
channel. An ONLINE estimator that runs every few seconds wants something
else: converge to a good state even while the gross error is present.
That is the robust weight modification (`robust = true`): per iteration,
each measurement's SOLVE weight is derived from its residual in units of
its original sigma, $t = |r_i| / \sigma_i$:

* $t \le 3$: weight unchanged,
* $3 < t \le 6$: sigma widened to $\sigma (2t/3 - 1)$ (a tangential
  transition zone),
* $t > 6$: sigma set to $|r_i|/3$, which effectively removes the row's
  pull on the state.

Stages are recomputed every iteration, so a row can recover; the
REPORTED statistics ($J$, normalized residuals) stay on the original
sigmas, so the diagnostics remain honest. `SEResult.robustRows` lists
every row that left stage 0.

````@example workshop_se_diagnostics
meas_bad = Measurement[m for m in net.measurements]   ## still contains the +25 MW error

res_plain = runse!(net, meas_bad; maxIte = 20, tol = 1e-8, updateNet = false)
res_robust = runse!(net, meas_bad; maxIte = 20, tol = 1e-8, updateNet = false, robust = true)

err_plain = maximum(abs.(abs.(res_plain.voltages) .- vm_true))
err_robust = maximum(abs.(abs.(res_robust.voltages) .- vm_true))
println("max |Vm_est - Vm_true|: plain = ", round(err_plain; sigdigits = 3), " pu, robust = ", round(err_robust; sigdigits = 3), " pu")
for rw in res_robust.robustRows
  println("  robust row: ", rw.id, "  stage ", rw.stage, ", t = ", round(rw.t; digits = 1), ", sigma widened by ", round(rw.sigma_factor; digits = 1), "x")
end
````

Reading aid (Example 4): the corrupted row sits in stage 2 (its gradient
contribution suppressed), and the robust state error is close to the
healthy-set level, while the plain WLS state is visibly bent. Rule of
thumb: robust for riding through, elimination for finding and fixing;
the two compose (a robust solve makes the residuals of the OTHER rows
cleaner, which sharpens the ranking).

## Current-magnitude measurements

**Example 5: the cheapest extra redundancy.** Substations meter bay
CURRENTS almost everywhere, even where no power transducer is installed.
A current magnitude (`ImagMeas`, in amperes) adds real redundancy, but
it comes with two structural quirks and two safety gates:

* A current magnitude has NO sign: at light flow, $|I|$ barely responds
  to the flow direction, and near zero its derivative blows up. Gate 1:
  `ImagMeas` rows enter the iteration only from
  `state_estimation.imag_activation_iteration` on (default 2), after the
  state has left the flat start. Gate 2: rows whose value is below three
  times their sigma stay out entirely (the near-zero region).
* Because $|I|$ carries no phase information, current measurements are
  EXCLUDED from the observability analysis: the system must be
  observable without them; they only sharpen it.

````@example workshop_se_diagnostics
p12 = get_branch_p_from_to_mw(net, "B1", "B2")
q12 = get_branch_q_from_to_mvar(net, "B1", "B2")
vm1 = net.nodeVec[net.busDict["B1"]]._vm_pu
i12_A = 1000.0 * sqrt(p12^2 + q12^2) / (sqrt(3.0) * 110.0 * vm1)
addImagMeasurement!(net; fromBus = "B1", toBus = "B2", direction = :from, value = i12_A, sigma = 5.0)
println("bay current B1->B2: ", round(i12_A; digits = 1), " A added as ImagMeas")

gobs = evaluate_global_observability(net; flatstart = true, jacEps = 1e-6)
println("observability rows: ", gobs.n_measurements, " of ", length(net.measurements), " active measurements (ImagMeas excluded)")
rep_i = validate_measurements(net)
println("diagnostics rows:   ", length(rep_i.measurement_ranking), " (ImagMeas included, dof = ", rep_i.objective.dof, ")")
````

Reading aid (Example 5): the observability count is one short of the
active measurement count (the current row), while the diagnostics count
includes it: the estimator uses the information, the observability
guarantee never depends on it. With active current rows the residual
statistics become operating-point dependent, one more reason the
activation gate waits for a settled state.

## Parameter estimation: a shunt reactor's susceptance

**Example 6: when the bad data is the MODEL.** A 30 MVar compensation
reactor hangs at `B5`, but the model database carries its susceptance
20 percent too high (a nameplate-vs-asbuilt classic). No telemetry row
is wrong, yet every run misfits around `B5`: bad-data logic would start
eliminating healthy measurements. The honest fix is to ESTIMATE the
parameter: `setShuntEstimation!` releases the shunt's susceptance $B$ as
an additional state, the bay's reactive-power measurement
(`ShuntQMeas`) provides the information, and the estimator recovers the
true value. Write-back into the model is a separate, explicit decision
(`updateShunts`); by default the estimate is only REPORTED.

````@example workshop_se_diagnostics
net_sh = Net(name = "workshop_se_shunt", baseMVA = 100.0)
addBus!(net = net_sh, busName = "S1", vn_kV = 110.0, vm_pu = 1.02, va_deg = 0.0)
addBus!(net = net_sh, busName = "S2", vn_kV = 110.0)
addBus!(net = net_sh, busName = "S3", vn_kV = 110.0)
addPIModelACLine!(net = net_sh, fromBus = "S1", toBus = "S2", r_pu = 0.010, x_pu = 0.080, b_pu = 0.0, status = 1)
addPIModelACLine!(net = net_sh, fromBus = "S2", toBus = "S3", r_pu = 0.011, x_pu = 0.085, b_pu = 0.0, status = 1)
addPIModelACLine!(net = net_sh, fromBus = "S1", toBus = "S3", r_pu = 0.010, x_pu = 0.080, b_pu = 0.0, status = 1)
addProsumer!(net = net_sh, busName = "S1", type = "EXTERNALNETWORKINJECTION", referencePri = "S1", vm_pu = 1.02, va_deg = 0.0)
addProsumer!(net = net_sh, busName = "S2", type = "ENERGYCONSUMER", p = 35.0, q = 12.0)
addProsumer!(net = net_sh, busName = "S3", type = "ENERGYCONSUMER", p = 28.0, q = 9.0)
addShunt!(net = net_sh, busName = "S3", pShunt = 0.0, qShunt = 30.0)   ## the reactor, TRUE value
ok_sh, msg_sh = validate!(net = net_sh)
ok_sh || error("shunt demo net invalid: $msg_sh")
runpf!(net_sh, 40, 1e-10, 0)

# measurements derived from the TRUE state, including the bay Q measurement
std_sh = measurementStdDevs(vm = 1e-3, pinj = 0.1, qinj = 0.1, pflow = 0.1, qflow = 0.1, shuntq = 0.1)
meas_sh = generateMeasurementsFromPF(net_sh; includeShuntQ = true, noise = false, stddev = std_sh)

# now poison the MODEL: +20 percent on the susceptance, then release it
sh = net_sh.shuntVec[1]
b_true = imag(sh.y_pu_shunt)
sh.y_pu_shunt = complex(real(sh.y_pu_shunt), 1.2 * b_true)
setShuntEstimation!(net_sh; busName = "S3")

res_sh = runse!(net_sh, meas_sh; maxIte = 30, tol = 1e-10, updateNet = false)
row = only(res_sh.shuntEstimates)
println("B_model = ", round(row.B_model; digits = 4), " pu (stale)   B_est = ", round(row.B_est; digits = 4), " pu   B_true = ", round(b_true; digits = 4), " pu")
println("model after the run (untouched): ", round(imag(sh.y_pu_shunt); digits = 4), " pu")
runse!(net_sh, meas_sh; maxIte = 30, tol = 1e-10, updateNet = false, updateShunts = true)
println("after updateShunts = true:       ", round(imag(sh.y_pu_shunt); digits = 4), " pu")
````

Reading aid (Example 6): the estimate lands on the true susceptance
while the model keeps its stale value until you explicitly opt into the
write-back. Two practical notes: without a direct bay Q measurement,
`deriveShuntPseudoMeasurements!` can derive one from the bay current and
a LOCAL voltage measurement (the derived row is protected from
elimination, it encodes model knowledge); and a released shunt whose bay
carries no information is FROZEN rather than estimated blindly, the
report row says so (`frozen = true`).

## Where to go next

- The full measurement workflow around these tools: measurement sets
  travel as CSV files (`writeMeasurementsCSV` /
  `readMeasurementsCSV!`), the local Web UI has a state-estimation page
  with synthetic-set generation, and a converged estimate can seed a
  power flow and an N-1 analysis (`runpf_from_se!`). See
  [State Estimation](https://welthulk.github.io/Sparlectra.jl/state_estimation/).
- Busbar couplers and FACTS devices in the measurement model (link
  contraction, flow allocation, the `se_view` report):
  [State Estimation](https://welthulk.github.io/Sparlectra.jl/state_estimation/),
  sections on links and FACTS.
- Every `state_estimation.*` key, including the robust and Takahashi
  thresholds used above:
  [State-Estimation Configuration](https://welthulk.github.io/Sparlectra.jl/state_estimation_configuration/).
- New to the estimator? Start with the
  [state-estimation basics notebook](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_state_estimation.ipynb):
  the closed loop from network build to estimated state, and the
  observability deep dive.

