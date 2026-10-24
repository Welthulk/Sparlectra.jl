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

> **Note:** On Google Colab the install cell takes a minute or two on a
> fresh session (package download and precompilation). Colab's Julia
> version may change over time; this notebook targets Julia >= 1.12.

## Warm-up

Julia compiles each function on first use. This cell loads the packages.
`Random` (standard library) seeds the synthetic measurement noise so
every run of this notebook produces the same numbers; `Printf` (standard
library) aligns the printed values. The notebook needs
no helper functions; every chapter works on the study network of the
next section, so after this cell and that one any chapter can run on
its own. A code cell that uses a name from another cell says so in a
comment. Its last line, `warmup()`, runs every path the notebook
exercises once on tiny throwaway networks (the power flow, synthetic
measurements, the diagnostics report, the elimination loop, the plain
and robust estimator, the observability check with a current
measurement, and the shunt parameter estimation), so the real study
runs at full speed; its code is in
`docs/lit/warmup/workshop_se_diagnostics.jl` of the package. How long
the warm-up takes depends on the machine: a Colab session is several
times slower than a desktop; every later cell then runs without compile
pauses.

````@example workshop_se_diagnostics
using Sparlectra
using Random
using Printf

# compile every path the chapters use once; the warm-up code is in
# docs/lit/warmup/workshop_se_diagnostics.jl of the installed package
include(pkgdir(Sparlectra, "docs", "lit", "warmup", "workshop_se_diagnostics.jl"))
warmup()
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
println("reference power flow (true state):  ", ite_pf, " iterations")
println("measurement set from the true state: ", length(net.measurements), " rows, seeded noise")
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
# uses: net (The study network and a healthy measurement set)
@isdefined(net) || error("Run the section \"The study network and a healthy measurement set\" first: it sets up net.")
rep = validate_measurements(net)
println("Healthy set, band test of the estimate:")
println("  estimator converged:                   ", rep.converged)
println("  global consistency:                    ", rep.global_consistency)
@printf("  objective J [-]:                       %.2f\n", rep.objective.value)
@printf("  degrees of freedom dof [-]:            %d\n", rep.objective.dof)
@printf("  J/dof [-] (healthy: near 1):           %.3f\n", rep.objective.value / rep.objective.dof)
@printf("  z_wh [-] (in band: |z_wh| <= 3):       %.3f   reason = %s\n", rep.objective.z_wh, rep.objective.reason)
println("  suspicious measurements (|r_N| >= 3):  ", length(rep.suspicious_measurements))
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
# uses: rep (Example 1)
@isdefined(rep) || error("Run Example 1 first: it sets up rep.")
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
# uses: net (The study network and a healthy measurement set)
@isdefined(net) || error("Run the section \"The study network and a healthy measurement set\" first: it sets up net.")
bad_idx = findfirst(m -> m.typ == Sparlectra.PflowMeas, net.measurements)
m0 = net.measurements[bad_idx]
net.measurements[bad_idx] = Measurement(typ = m0.typ, value = m0.value + 25.0, sigma = m0.sigma, busIdx = m0.busIdx, branchIdx = m0.branchIdx, direction = m0.direction, id = m0.id, linkIdx = m0.linkIdx)
println("Measurement ", m0.id, " (line B1-B2 at B1), active power flow P [MW]:")
@printf("  measured, healthy (with noise):   %6.2f\n", m0.value)
@printf("  measured after corruption:        %6.2f   (+25 MW)\n", m0.value + 25.0)
@printf("  declared sigma:                   %6.2f\n", m0.sigma)

rep_bad = validate_measurements(net)
println()
println("Band test of the estimate with the corrupted row:")
@printf("  J/dof [-]:                        %6.1f\n", rep_bad.objective.value / rep_bad.objective.dof)
@printf("  z_wh [-]:                         %6.1f   reason = %s\n", rep_bad.objective.z_wh, rep_bad.objective.reason)
@printf("  suspicious rows (|r_N| >= 3):     %6d\n", length(rep_bad.suspicious_measurements))
top = rep_bad.largest_normalized_residual
println("  prime suspect (largest |r_N|):    ", top.id)
println("Ranking, top 5 (localizable if w_ii > 0.3):")
println("  measurement             |r_N| [-]   w_ii [-]   localizable")
for row in rep_bad.measurement_ranking[1:5]
  @printf("  %-22s %9.1f %10.2f   %s\n", row.id, row.abs_normalized_residual, row.wii, row.localizable)
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
# uses: net (The study network and a healthy measurement set)
@isdefined(net) || error("Run the section \"The study network and a healthy measurement set\" first: it sets up net.")
diag = with_state_estimation_config(max_eliminations = 3) do
  runse_diagnostics(net)
end
println("Sequential elimination (max_eliminations = 3): ", length(diag.eliminations), " row(s) removed, stop reason = ", diag.stop_reason)
println("  step   removed measurement      |r_N| before removal [-]")
for e in diag.eliminations
  @printf("  %4d   %-22s %8.1f\n", e.elimination, e.id, abs(e.normalized_residual_before))
end
fin = diag.final_diagnostics
println("Band test of the estimate after elimination:")
@printf("  J/dof [-]:                        %6.3f\n", fin.objective.value / fin.objective.dof)
@printf("  z_wh [-]:                         %6.3f   reason = %s\n", fin.objective.z_wh, fin.objective.reason)
@printf("  suspicious rows (|r_N| >= 3):     %6d\n", length(fin.suspicious_measurements))
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
# uses: bad_idx, m0 (Example 2); net, vm_true (The study network and a healthy measurement set)
(@isdefined(bad_idx) && @isdefined(m0)) || error("Run Example 2 first: it sets up bad_idx, m0.")
(@isdefined(net) && @isdefined(vm_true)) || error("Run the section \"The study network and a healthy measurement set\" first: it sets up net, vm_true.")
meas_bad = Measurement[m for m in net.measurements]   ## still contains the +25 MW error
# the same set with the healthy row restored, as the reference level
meas_ok = copy(meas_bad)
meas_ok[bad_idx] = m0

res_ok = with_state_estimation_config(max_iter = 20, tol = 1e-6, update_net = false) do
  runse!(net, meas_ok)
end
res_plain = with_state_estimation_config(max_iter = 20, tol = 1e-6, update_net = false) do
  runse!(net, meas_bad)
end
res_robust = with_state_estimation_config(max_iter = 20, tol = 1e-6, update_net = false, robust = true) do
  runse!(net, meas_bad)
end

err_ok = maximum(abs.(abs.(res_ok.voltages) .- vm_true))
err_plain = maximum(abs.(abs.(res_plain.voltages) .- vm_true))
err_robust = maximum(abs.(abs.(res_robust.voltages) .- vm_true))
println("Largest voltage magnitude error over B1..B7, max |Vm_est - Vm_true| [pu]:")
@printf("  plain WLS, healthy set (reference):         %.5f\n", err_ok)
@printf("  plain WLS, corrupted row at full weight:    %.5f\n", err_plain)
@printf("  robust WLS, corrupted row down-weighted:    %.5f\n", err_robust)
println("Robust rows (left stage 0), t = |r| / sigma:")
println("  measurement             stage   t [-]   sigma factor [-]")
for rw in res_robust.robustRows
  @printf("  %-22s %5d %7.1f %12.1f\n", rw.id, rw.stage, rw.t, rw.sigma_factor)
end
````

Reading aid (Example 4): the two runs on the corrupted set open with a
warning from the topology precheck, `kcl_violation at bus B1`: the P
balance at `B1` misses by 27 MW, 19 sigma. That is the +25 MW row again,
found before the estimation starts: `B1` is completely measured
(injection and both branch flows), so its balance can be checked on the
raw values. The precheck is advisory, the estimation proceeds, and each
`runse!` call checks its set once, hence two warnings; the reference run
on the healthy set has none. Examples 2 and 3 did not show it because
the diagnostics run their solve without the precheck.

The corrupted row sits in stage 2 (its gradient contribution
suppressed). The robust voltage error (0.00074 pu) is at the level of
the healthy set (0.00076 pu), while the plain WLS error is about
1.5 times as large (0.00113 pu): the bad row bends the plain state. Rule
of thumb: robust for riding through, elimination for finding and fixing;
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
# uses: net (The study network and a healthy measurement set, Example 4)
@isdefined(net) || error("Run the section \"The study network and a healthy measurement set\", Example 4 first: they set up net.")
p12 = get_branch_p_from_to_mw(net, "B1", "B2")
q12 = get_branch_q_from_to_mvar(net, "B1", "B2")
vm1 = net.nodeVec[net.busDict["B1"]]._vm_pu
i12_A = 1000.0 * sqrt(p12^2 + q12^2) / (sqrt(3.0) * 110.0 * vm1)
addImagMeasurement!(net; fromBus = "B1", toBus = "B2", direction = :from, value = i12_A, sigma = 5.0)
println("Line B1-B2 at B1, current magnitude |I| [A]:")
@printf("  true (from the power flow), added as ImagMeas:  %.1f   sigma 5.0\n", i12_A)

gobs = with_state_estimation_config(flatstart = true, jac_eps = 1e-6) do
  evaluate_global_observability(net)
end
rep_i = validate_measurements(net)
println("Measurement rows used, of ", length(net.measurements), " active:")
println("  observability check (ImagMeas excluded):  ", gobs.n_measurements)
println("  diagnostics report (ImagMeas included):   ", length(rep_i.measurement_ranking), ", dof = ", rep_i.objective.dof)
````

Reading aid (Example 5): the observability count is one short of the
active measurement count (the current row), while the diagnostics count
includes it: the estimator uses the information, the observability
guarantee never depends on it. With active current rows the residual
statistics become operating-point dependent, one more reason the
activation gate waits for a settled state.

## Parameter estimation: a shunt reactor's susceptance

**Example 6: when the bad data is the MODEL.** A separate three-bus
network (`S1` to `S3`) has a 30 MVar compensation reactor at `S3`.
`addShunt!` follows the MATPOWER convention: `qShunt` is the reactive
power in MVar at 1 pu voltage, $B > 0$ is a capacitor, $B < 0$ a
reactor. The reactor therefore enters as `qShunt = -30`, which is
$B = -0.30$ pu on 100 MVA. The model database carries its susceptance
20 percent too large (factor 1.2, $B = -0.36$ pu; a nameplate-vs-asbuilt
classic). No telemetry row is wrong, yet every run misfits around `S3`:
bad-data logic would start eliminating healthy measurements. The fix is
to ESTIMATE the parameter: `setShuntEstimation!` releases the shunt's
susceptance $B$ as an additional state, the bay's reactive-power
measurement (`ShuntQMeas`) provides the information, and the estimator
recovers the true value. Whether a converged run writes the estimate
into the model is set by `state_estimation.update_shunts`; the default is
`true`. Run 1 below switches the write-back off to show the estimate on
its own, run 2 switches it on.

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
addShunt!(net = net_sh, busName = "S3", pShunt = 0.0, qShunt = -30.0)   ## the reactor (B < 0), TRUE value
ok_sh, msg_sh = validate!(net = net_sh)
ok_sh || error("shunt demo net invalid: $msg_sh")
runpf!(net_sh, 40, 1e-10, 0)

# measurements derived from the TRUE state, including the bay Q measurement
std_sh = measurementStdDevs(vm = 1e-3, pinj = 0.1, qinj = 0.1, pflow = 0.1, qflow = 0.1, shuntq = 0.1)
meas_sh = generateMeasurementsFromPF(net_sh; includeShuntQ = true, noise = false, stddev = std_sh)

# now falsify the MODEL: factor 1.2 on the susceptance, then release it
sh = net_sh.shuntVec[1]
b_true = imag(sh.y_pu_shunt)
sh.y_pu_shunt = complex(real(sh.y_pu_shunt), 1.2 * b_true)
b_model_before = imag(sh.y_pu_shunt)
setShuntEstimation!(net_sh; busName = "S3")

# run 1: estimate only, no write-back
res_sh = with_state_estimation_config(max_iter = 30, tol = 1e-6, update_net = false, update_shunts = false) do
  runse!(net_sh, meas_sh)
end
row = only(res_sh.shuntEstimates)
b_after_run1 = imag(sh.y_pu_shunt)
# run 2: the same estimation with write-back (the default)
with_state_estimation_config(max_iter = 30, tol = 1e-6, update_net = false, update_shunts = true) do
  runse!(net_sh, meas_sh)
end
b_after_run2 = imag(sh.y_pu_shunt)

println("Shunt at S3 (reactor), susceptance B [pu]:")
@printf("  true value (used to generate measurements):    %.2f\n", b_true)
@printf("  model before estimation (falsified):           %.2f\n", b_model_before)
@printf("  estimated by SE:                               %.2f\n", row.B_est)
@printf("  model after run 1 (update_shunts = false):     %.2f   %s\n", b_after_run1, b_after_run1 == b_model_before ? "unchanged" : "changed")
@printf("  model after run 2 (update_shunts = true):      %.2f   %s\n", b_after_run2, isapprox(b_after_run2, row.B_est; atol = 1e-4) ? "written back" : "not written back")
@printf("SE quality: J = %.1e, dof = %d\n", res_sh.objectiveJ, res_sh.dof)
````

Reading aid (Example 6): both runs open with a warning from the topology
precheck, `kcl_violation at bus S3`: the Q balance at `S3` misses by
5.9 MVar, 34 sigma. This is the wanted result. The measured injection
and bay flows at `S3` do not balance with the model's shunt, and the miss
is the model error itself (0.06 pu on 100 MVA, times $V^2$ at `S3`). The
precheck finds the wrong parameter at the right bus before the
estimation starts. Run 1 leaves the falsified model in place, so run 2
sees the same mismatch, hence two warnings.

The estimate lands on the true susceptance (-0.30 pu). Run 1 leaves the
model at the falsified -0.36 pu, run 2 writes the estimate back. $J$ is
practically zero: the measurements are noise-free, and with $B$ as a
state the model fits them exactly. Two practical notes: without a direct bay Q measurement,
`deriveShuntPseudoMeasurements!` can derive one from the bay current and
a LOCAL voltage measurement (the derived row is protected from
elimination, it encodes model knowledge); and a released shunt whose bay
carries no information is FROZEN rather than estimated blindly, the
report row says so (`frozen = true`).

## Where to go next

- The full measurement workflow around these tools: measurement sets
  travel as measurement tables (`writeMeasurementTable` /
  `readMeasurementTable!`), the local Web UI has a state-estimation page
  with synthetic-set generation, and a converged estimate can seed a
  power flow and an N-1 analysis (`runpf_from_se!`). See
  [State Estimation](https://welthulk.github.io/Sparlectra.jl/state_estimation/).
- Busbar couplers and FACTS devices in the measurement model (link
  contraction, link flow measurements as WLS rows, link flows after the
  estimate, the `se_view` report):
  [State Estimation Extensions](https://welthulk.github.io/Sparlectra.jl/state_estimation_extensions/).
- Every `state_estimation.*` key, including the robust and Takahashi
  thresholds used above:
  [State-Estimation Configuration](https://welthulk.github.io/Sparlectra.jl/state_estimation_configuration/).
- New to the estimator? Start with the
  [state-estimation basics notebook](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_state_estimation.ipynb):
  the closed loop from network build to estimated state, and the
  observability deep dive.

