# Solver Guide

The internal power-flow path is the sparse rectangular complex-state
Newton-Raphson solver. Its layers under `src/powerflow_rectangular/`:

* `rectangular_core_equations.jl`: residual/equation helpers, including `mismatch_rectangular`.
* `rectangular_jacobian_builders.jl`: the analytic rectangular Jacobian builders.
* `rectangular_newton_step.jl`: Newton-step update and damping helpers.
* `rectangular_standalone_solver.jl`: `run_complex_nr_rectangular`, the standalone array-level Newton driver.
* `rectangular_network_solver.jl`: `runpf_rectangular!`, the network-integrated entry point and orchestration.

## Rectangular Complex-State Newton-Raphson

### Motivation and State Vector

The rectangular complex-state NR uses the complex bus voltages

```math
V_i = V_{r,i} + j V_{i,i}
```

as state variables instead of magnitudes and angles. The state vector is

```math
x =
\begin{bmatrix}
V_r(\text{non-slack}) \\
V_i(\text{non-slack})
\end{bmatrix}
\in \mathbb{R}^{2(n-1)}.
```

The complex bus powers are

```math
I = Y_\mathrm{bus} V, \qquad
S = V \odot \overline{I} = V \odot \overline{Y_\mathrm{bus} V}
```

and the specified injections

```math
S_\mathrm{spec} = P_\mathrm{spec} + j Q_\mathrm{spec}.
```

### Bus Types and Residual Definition

```julia
mismatch_rectangular(Ybus, V, S, bus_types, Vset, slack_idx)
    -> Vector{Float64}
```

builds the real-valued residual vector `F(V)`:

* For each **PQ bus** $i \ne \text{slack}$:

  ```math
  \Delta P_i = \Re(S_{\text{calc},i}) - \Re(S_{\text{spec},i})
  \\
  \Delta Q_i = \Im(S_{\text{calc},i}) - \Im(S_{\text{spec},i})
  ```

* For each **PV bus** $i \ne \text{slack}$:

  ```math
  \Delta P_i = \Re(S_{\text{calc},i}) - \Re(S_{\text{spec},i})
  \\
  \Delta V_i = |V_i| - V_\text{set,i}
  ```

The stacked residual vector has size `2(n-1)` and matches the state dimension.

### Analytic Rectangular Newton Step

`complex_newton_step_rectangular` performs one Newton step with an analytic
Jacobian:

```julia
complex_newton_step_rectangular(Ybus, V, S;
    slack_idx::Int = 1,
    damp::Float64  = 1.0,
    autodamp::Bool = false,
    autodamp_min::Float64 = 0.05,
    bus_types::Vector{Symbol},
    Vset::Vector{Float64},
)
```

Steps:

1. Compute currents and powers.
2. Form the residual `F(V)`.
3. Assemble the Jacobian `J = ∂F/∂x`.
4. Solve

   ```math
   J(x_k)\,\Delta x_k = -F(x_k)
   ```

5. Update the state with fixed or automatic damping.

### Automatic Rectangular Newton Damping

For difficult flat-start studies, enable residual-based backtracking:

```julia
runpf!(net, 60, 1e-8, 1;
    opt_flatstart = true,
    autodamp = true,
    autodamp_min = 0.05,
)
```

With `autodamp = true` the solver tests trial step lengths from `damp` down
to `autodamp_min` by halving. The first trial that reduces the maximum
absolute mismatch is accepted; if none improves, the solver continues with
the best finite conservative trial so the next iteration can rebuild the
Jacobian and active-set state. This line-search style damping only chooses
the scalar step length of the computed Newton correction; model, bus
equations, Q-limit logic and Jacobian are untouched.

### Merit-Function Line Search

Autodamp accepts the first trial step that reduces the **maximum absolute
mismatch**, an $\infty$-norm criterion without a sufficient-decrease
guarantee. `power_flow.merit` adds an opt-in acceptance test based on a
scalar merit function inside the same backtracking loop.

| Use | Where |
|---|---|
| Config key | `power_flow.merit.enabled` (default `false`); `armijo_c1`, `fallback_max_mismatch`, `scale_p` / `scale_q` / `scale_v` under the same prefix |
| Log field | `accept_reason` per iteration ([Merit-function line search options](@ref pf-merit)) |

#### Root-finding vs. optimization view

Newton-Raphson solves $F(x) = 0$. The same residual defines the scalar merit
function

```math
f(x) = \frac{1}{2} \lVert W F(x) \rVert_2^2
```

where $W$ is the optional diagonal scaling matrix of the residual-scaling
weights `scale_p`/`scale_q`/`scale_v`. Minimizing $f$
and solving $F(x) = 0$ share their solutions whenever $F(x) = 0$ is
attainable ($f = 0$ there).

The existing Newton direction is a descent direction for $f$ without any
extra factorization, because of the gradient

```math
\nabla f(x) = J(x)^\top W^\top W F(x)
```

With the Newton direction $\Delta x = -J(x)^{-1} F(x)$ and $W = I$:

```math
\nabla f(x)^\top \Delta x = -F(x)^\top J(x)^{-\top} J(x)^\top J(x)^\top J(x)^{-1} F(x)
```

reduces, using $J \Delta x = -F$, to

```math
\nabla f(x)^\top \Delta x = -\lVert W F(x) \rVert_2^2 < 0
```

whenever $F(x) \ne 0$. The directional derivative comes directly from the
already computed mismatch vector.

#### Armijo sufficient-decrease acceptance

For each trial step length $\lambda$ (same order as autodamp: `damp`, then
halved down to `autodamp_min`), the merit line search accepts the first
$\lambda$ satisfying the Armijo condition:

```math
f(x + \lambda \Delta x) \le f(x) + c_1 \, \lambda \, \nabla f(x)^\top \Delta x
```

`power_flow.merit.armijo_c1` is $c_1 \in (0, 0.5)$. Smaller $c_1$ accepts
almost any step that decreases $f$; larger $c_1$ demands a larger fraction
of the predicted linear decrease and backtracks further (Nocedal & Wright,
*Numerical Optimization*, line-search chapter).

Armijo acceptance guarantees descent on the aggregate residual energy $f$;
the max-mismatch criterion only asks whether the single worst bus equation
improved. A step can reduce the worst-bus mismatch while increasing the
overall residual energy elsewhere, or vice versa, so the two criteria can
accept different step lengths in the same iteration.

If no trial satisfies Armijo, `power_flow.merit.fallback_max_mismatch`
selects the fallback:

* `true` (default): the classic max-mismatch criterion for this iteration
  (first improving trial, else the conservative best-finite trial);
* `false`: straight to the conservative best-finite trial.

#### Limits

* Minimizing $f$ is not equivalent to solving $F(x) = 0$: $f$ can have local
  minima with $F(x) \ne 0$, so Armijo at every iteration does not guarantee
  convergence to a root.
* The criterion does not influence *which* solution branch (a high- versus
  a low-voltage operating point) the solver converges to; it only accepts
  or rejects step lengths along Newton's path.
* A PV/PQ active-set switch changes what residual entries *mean* (a
  $\Delta Q$ entry becomes a $\Delta V$ entry, or vice versa). Comparing $f$
  across such a switch is not well-defined, so the merit comparison is
  skipped in that iteration, which falls back to the max-mismatch criterion
  and is logged with `accept_reason = active_set_skip` (log format:
  [Merit-function line search options](@ref pf-merit)).

#### Residual scaling

$P$, $Q$, and voltage-setpoint ($V$) residual entries can differ by orders
of magnitude depending on the per-unit base and grid mix, and one equation
type can then dominate $f$ and weaken the sufficient-decrease guarantee for
the others. `power_flow.merit.scale_p`, `power_flow.merit.scale_q`, and
`power_flow.merit.scale_v` are the diagonal weights $W$ per equation type:
`ΔP_i` entries use `scale_p`; the second equation per non-slack bus uses
`scale_q` for PQ buses and `scale_v` for PV buses. Leave all three at `1.0`
unless diagnostics show one residual type dominating.

### Trust-Region Step Control

`power_flow.trust_region` is an alternative step control to autodamp: a
**scaled-Newton trust region**. It is disabled by default and mutually
exclusive with `power_flow.autodamp`: both decide how far to step along the
Newton direction, so enabling both is a configuration error.

| Use | Where |
|---|---|
| Config key | `power_flow.trust_region.enabled` (default `false`); `step_mode` (`:scaled`, `:dogleg`), `eta_accept`, `expand_threshold`, `expand_factor`, `shrink_factor`, `min_radius` under the same prefix |
| Result field | `accept_reason` (`:dogleg_newton`, `:dogleg_cauchy`, `:dogleg_interp`, `:active_set_skip`); non-convergence reason `:trust_region_collapsed` |

#### Step-construction modes: `scaled` and `dogleg`

Classical trust-region methods restrict the Newton correction to a ball of
radius $\Delta$ around the current iterate, typically choosing the trial
step by a *dogleg* or *Steihaug-CG* interpolation between the
steepest-descent and full Newton directions.
`power_flow.trust_region.step_mode` selects between two constructions; both
share the acceptance rule below.

##### `step_mode = :scaled` (default)

The trial step is the full Newton direction $\Delta x$, rescaled when it
exceeds the radius:

```math
\Delta x_{\text{scaled}} =
\begin{cases}
\Delta x & \lVert \Delta x \rVert \le \Delta \\
\Delta \dfrac{\Delta x}{\lVert \Delta x \rVert} & \lVert \Delta x \rVert > \Delta
\end{cases}
```

The direction is always the analytic Newton direction, as in the
autodamp/fixed-damping paths; only the step *length* is controlled.

##### `step_mode = :dogleg`

`scaled` steps degrade when the Newton direction stops being a useful
descent direction (for example with a very ill-conditioned Jacobian): every
rescale points the same bad way, so the solver either keeps taking uphill
steps or shrinks the radius toward collapse. Dogleg adds a second endpoint,
the steepest-descent (Cauchy) step on the local quadratic model of the merit
function, and interpolates along the Cauchy-to-Newton path.

**Merit-function gradient.** With $f(x) = \frac{1}{2}\lVert F(x) \rVert_2^2$
(unweighted, $W = I$, the convention `m(x)` uses for both step modes), the
gradient is

```math
g(x) = J(x)^{\mathsf{T}} F(x)
```

one transposed sparse matrix-vector product, no additional factorization.

**Cauchy point.** The steepest-descent minimizer of the local quadratic model
$m(p) = f + g^{\mathsf{T}}p + \frac{1}{2}p^{\mathsf{T}}(J^{\mathsf{T}}J)p$
along $-g$ is

```math
\alpha^\ast = \frac{\lVert g \rVert_2^2}{\lVert Jg \rVert_2^2}
\qquad p_C = -\alpha^\ast g
```

one more sparse matrix-vector product ($Jg$), clipped to
$\lVert p_C \rVert \le \Delta$ when used as the trial step.

**Dogleg path.** With the Newton step $p_N = \Delta x$ and radius $\Delta$:

1. $\lVert p_N \rVert \le \Delta$: the full Newton step
   (`accept_reason = :dogleg_newton`, identical to `scaled` where the radius
   does not bind).
2. $\lVert p_C \rVert \ge \Delta$: $p_C$ rescaled to length $\Delta$
   (`accept_reason = :dogleg_cauchy`, pure steepest descent).
3. Otherwise: interpolate $p(\tau) = p_C + \tau (p_N - p_C)$, $\tau \in
   [0,1]$, with $\tau$ such that $\lVert p(\tau) \rVert = \Delta$
   (`accept_reason = :dogleg_interp`). This is the scalar quadratic
   $a\tau^2 + b\tau + c = 0$ with
   $a = \lVert d \rVert_2^2$, $b = 2\, p_C^{\mathsf{T}} d$,
   $c = \lVert p_C \rVert_2^2 - \Delta^2$, $d = p_N - p_C$, solved in closed
   form.

Along the dogleg path $\lVert p(\tau) \rVert$ increases monotonically and
the model value $m(p(\tau))$ decreases monotonically **provided**
$J^{\mathsf{T}}J$ is positive definite (Nocedal & Wright, trust-region
chapter). Near a singular Jacobian, the regime dogleg is meant to help with,
$J^{\mathsf{T}}J$ may only be positive *semidefinite*; the path is still
defined, but its Newton endpoint can be a poor local model. Dogleg therefore
buys graceful degradation (Cauchy-direction progress) for *transiently*
ill-conditioned Jacobians, not a global-convergence guarantee. The two
escalations in the literature, Levenberg-Marquardt
($(J^{\mathsf{T}}J + \lambda I) \Delta x = -J^{\mathsf{T}}F$, a new
factorization per $\lambda$-change) and Steihaug-CG (matrix-free, terminates
on negative curvature or the radius boundary), are not implemented.

**Active-set switches.** A PV/PQ switch changes what the residual entries
mean, so a gradient or Cauchy point computed before the switch is invalid
after it. As in the merit line search, the dogleg comparison is skipped when
a switch happened in this iteration: a single `scaled`-style trial at the
current radius is accepted unconditionally
(`accept_reason = :active_set_skip`) and the radius stays unchanged.

#### Acceptance by merit decrease

Step acceptance reuses the merit function of the
[Merit-Function Line Search](@ref), $f(x) = \frac{1}{2}\lVert F(x) \rVert_2^2$
(unweighted, $W = I$). For a trial step $p$ (the rescaled Newton step or the
dogleg step), the **actual reduction** is

```math
\text{ared} = f(x) - f(x + p)
```

and the **predicted reduction** comes from the local linear model of $F$
around $x$:

```math
\text{pred} = f(x) - \frac{1}{2} \lVert F(x) + J(x)\, p \rVert_2^2
```

The ratio $\rho = \text{ared} / \max(\text{pred}, \varepsilon)$ measures how
well the linear model predicted the improvement. A trial is accepted when
$\rho \ge \eta_{\text{accept}}$ (`power_flow.trust_region.eta_accept`);
`pred` uses a matrix-vector product with the already factored Jacobian, not
a second Jacobian build or linear solve.

#### Radius update

The radius $\Delta$ persists across Newton iterations (autodamp restarts
`damp` per iteration) and adapts from $\rho$:

* **Accepted** ($\rho \ge \eta_{\text{accept}}$): the step is taken. If it
  also hit the radius boundary ($\lVert \Delta x \rVert > \Delta$ in
  `scaled` mode, or `accept_reason` is `:dogleg_cauchy`/`:dogleg_interp` in
  `dogleg` mode) and $\rho \ge$ `expand_threshold`, the radius expands:
  $\Delta \leftarrow \min(\Delta \cdot \texttt{expand\_factor}, \Delta_{\max})$.
  Otherwise the radius is unchanged.
* **Rejected** ($\rho < \eta_{\text{accept}}$): the radius shrinks,
  $\Delta \leftarrow \Delta \cdot \texttt{shrink\_factor}$, and the *same*
  Newton iteration retries with the smaller radius: no new Jacobian build or
  linear solve, only a re-evaluation of the step construction and a fresh
  mismatch evaluation, as autodamp reuses one Newton correction across its
  trials.
* **Collapsed**: if $\Delta$ falls below `power_flow.trust_region.min_radius`
  without an accepted step, the solver declares non-convergence with reason
  `:trust_region_collapsed` instead of looping indefinitely.

#### Limits

* Neither mode is a global-convergence guarantee or a wrong-branch selector;
  both are step-length controls on the existing analytic Newton direction,
  and `dogleg` does not fix a bad start (a start-mode question).
* `dogleg`'s Cauchy (steepest-descent) steps can be slow on badly scaled
  problems, inherent to the method. The dogleg gradient always uses $W = I$;
  weighted descent exists only via `power_flow.merit.scale_p/q/v` on the
  (mutually exclusive) merit-function line search.

General theory: Nocedal & Wright, *Numerical Optimization*, the chapters on
trust-region methods (dogleg, Cauchy point, Steihaug-CG) and line-search
methods (Armijo).

### Start Projection for Difficult Seeds

#### DC-angle flat-start background

A conventional AC flat start sets bus voltages near $1.0\,\mathrm{pu}$ and
all non-slack angles near $0^\circ$. That seed ignores the active-power flow
pattern encoded in topology, branch reactances, transformer phase shifts and
injections; in large or heavily phase-shifted cases the first Newton step
can start far from the physically relevant angle branch.

The DC-angle seed uses the active-power part of the model as a linearized
predictor for the angles: voltage magnitudes near nominal, resistance and
reactive coupling neglected, branch active power approximated by the angle
difference over the branch reactance. This gives the sparse linear system

```math
B'\theta = P
```

where $P$ is the net active-power injection vector and $B'$ is assembled
from the branch susceptance structure. The slack angle fixes the reference;
the angles are clipped by `start_projection_dc_angle_limit_deg` before
seeding Newton.

Blended starts combine the DC-angle predictor with the stored MATPOWER
`VM`/`VA` data or the raw flat start. The projection scans the requested
candidates and keeps the one with the smallest rectangular mismatch, which
helps avoid wrong low-voltage or wrong-angle branches without changing the
convergence criteria:

```julia
runpf!(net, 60, 1e-8, 1;
    opt_flatstart = true,
    start_projection = true,
    start_projection_try_dc_start = true,
    start_projection_try_blend_scan = true,
    start_projection_blend_lambdas = [0.25, 0.5, 0.75],
    start_projection_dc_angle_limit_deg = 60.0,
)
```

The projection sanitizes the raw seed before the scan. The same options are
available through `buildPfModel` and `runpf_external!`, so external solvers
receive the projected `model.V0`.

## Distributed Active-Power Slack

### Why a single slack bus is a modeling artifact

The classical formulation fixes the reference bus's voltage and lets its
computed injection absorb whatever the network needs beyond the scheduled
injections (imbalance plus losses). Physically no single machine does this:
primary control raises many generators together, each by its droop share. The single-slack model concentrates that response on one
bus and distorts the flows around it, visibly so in large imported cases
where the reference machine absorbs tens of MW that real dispatch would
spread over an area.

The distributed slack uses the primary-control picture: a set of
**participants** shares the imbalance according to normalized
**participation factors** $\alpha_i$ with $\sum_i \alpha_i = 1$.

### Formulation

One scalar unknown $\lambda_P$ (the island's total imbalance in p.u.) is
added to the state. Every participant bus $i$ (including the reference bus
if it participates) gets the modified active-power residual

```math
r_{P,i} = P_{\text{calc},i}(V) - P_{\text{spec},i} - \alpha_i\,\lambda_P
```

so its solved injection is its schedule plus its share of the imbalance.
Non-participants keep their residual ($\alpha_i = 0$).

The system stays square through a **role separation at the reference bus**:
its voltage magnitude *and angle* stay fixed (the angle reference is
untouched), but its active power stops being free. The reference bus's P
residual becomes an ordinary equation, **appended as the last row** of the
residual vector (the interleaved `[ΔP_i, ΔQ/ΔV_i]` layout of the `2(n-1)`
voltage rows is untouched). One new unknown ($\lambda_P$), one new equation
(REF-P): the Jacobian gains one column ($-\alpha_i$ at the participant P
rows) and one row (the network derivatives
$\partial P_\text{ref}/\partial V_{r,j}, \partial P_\text{ref}/\partial V_{i,j}$,
plus $-\alpha_\text{ref}$ in the $\lambda_P$ column).

With all weight on the reference bus ($\alpha_\text{ref} = 1$) the appended
equation reads $P_{\text{calc,ref}} - P_{\text{spec,ref}} - \lambda_P = 0$
and decouples from the voltage solution: the classical single-slack result
is reproduced exactly, with $\lambda_P$ reporting the classically
REF-absorbed power. This equivalence is a regression test.

### When and how it applies

| Use | Where |
|---|---|
| Config key | `power_flow.distributed_slack.enabled` (default `false`); `p_mode`, `fallback`, `respect_p_limits` under the same prefix ([Power-Flow Configuration](powerflow_configuration.md)) |
| Result field | $\lambda_P$ in the per-island solver statuses; P-limit warnings counted in the run metadata |

* **Disabled path.** Imported participation factors (MATPOWER `APF`, CGMES
  `GeneratingUnit.normalPF`) are carried on the `ProSumer` but never
  activate the feature; the disabled path is bit-identical to the classical
  solver.
* **Participants** are the generator-type prosumers at the island's REF and
  PV buses at solve start. Fixed injections at PQ buses (Stage-0 HVDC
  converter injections, kept boundary equivalents) never participate.
  Participation is frozen per solve: a participant that hits its Q limit
  and switches PV to PQ keeps its $\alpha_i$; reactive saturation does not
  remove a machine from primary control.
* **Per island.** Every island builds its own participant set and solves for
  its own $\lambda_P$. Islands without a valid participant abort
  (`fallback: error`) or solve classically with a warning
  (`fallback: ref_only`).
* **Weight modes** (`p_mode`): `pg_weighted` (scheduled `Pg`, the default),
  `pmax_weighted`, `headroom_weighted` (`max(maxP - Pg, 0)`), `imported`
  (the `APF`/`normalPF` factor), `explicit` (a config table).
* **Step control.** Autodamp, merit line search and trust region evaluate
  trial states $V + a\,\delta V$; the trial multiplier
  $\lambda_P + a\,\delta\lambda_P$ is staged alongside. In the weighted
  merit norm the REF-P row is scaled with `merit.scale_p` like every other
  active-power residual.
* **P limits are advisory**: with `respect_p_limits = true` each participant
  whose corrected output leaves `[minP, maxP]` produces a warning and is
  counted in the metadata; there is no clamp-and-resolve round.
* The $\lambda$-augmented system requires the sparse Jacobian path (the
  default); the dense builder rejects it.
* **Comparing against a single-slack reference state** (for example a CGMES
  SV profile): the branch flows deviate by the distributed correction while
  the voltages stay put; the comparison then measures the slack-distribution
  difference, not solver error.

## Linear solver backends

The linear solve of the Newton step (`J · δx = -F`) supports two sparse
backends.

| Use | Where |
|---|---|
| Config key | `power_flow.linear_solver`: `umfpack_reuse` (default) or `umfpack`; independent of `power_flow.solver`, which selects the power-flow method (rectangular / apslf / dc) |
| Result field | backend plus the counters `linear_solver_analyze_count`, `linear_solver_refactor_count`, `linear_solver_fallback_count` in the solver status and the performance profile, per island on the island path |

* **`umfpack`**: the standard sparse direct solve (`J \ F` through
  `solve_sparse_system`), with sparse-QR and small-system SVD fallbacks. Each
  iteration pays full symbolic analysis plus numeric factorization.
* **`umfpack_reuse`**: the same UMFPACK multifrontal factorization, but the
  symbolic analysis of the first iteration is kept and reused via
  `lu!(F, J)` (numeric refactorization only). The Jacobian pattern is
  constant across NR iterations unless the active set changes, so the
  analysis is paid once per active set; usually the fastest choice on
  large cases.

A `klu` backend is not offered: KLU is slower than UMFPACK on power-flow
Jacobians, and a KLU factorization shared across threads gives silently
wrong results.

Behavior of the reuse backend:

* A PV/PQ switch changes the sparsity pattern, so that iteration re-runs
  analyze plus factor. Independent of that signal, the stored
  `colptr`/`rowval` fingerprint of the analyzed pattern is compared against
  the current Jacobian before every refactorization, and any drift triggers
  a re-analyze.
* The Jacobian is assembled with a value-independent sparsity pattern, so
  the pattern depends only on the Ybus structure and the bus types.
* Any factorization error (structure mismatch, bad pivots, singularity)
  resets the context and delegates that solve to the `umfpack` chain, so
  the `:singular_newton_step` handling and the QR/SVD fallbacks keep
  working.
* The context covers only the main Newton solve (the dogleg/trust-region
  internals keep their path), lives for one `runpf_rectangular!` invocation
  (one island on the island path), and is never shared across islands or
  threads.

## Solver-Specific Interaction with Power Limits

When a PV bus hits a Q-limit, the solver does not merely clamp a reported
reactive power. It changes the equation type of that bus from `(ΔP, ΔV)`
(PV) to `(ΔP, ΔQ)` (PQ), represented through the `bus_types` vector.

The start of PV to PQ switching can be controlled for difficult cases:

```julia
runpf!(net, 60, 1e-8, 1;
    qlimit_start_iter = 4,
    qlimit_start_mode = :iteration,
)
```

`qlimit_start_mode = :auto` waits until the PV reactive-power requests have
stabilized; the threshold is `qlimit_auto_q_delta_pu` in p.u. The hybrid
`:iteration_or_auto` mode starts when either the configured iteration or the
stabilization criterion is reached. Operational and data-model side:
[Powerlimits Guide](powerlimits.md).

## Voltage-dependent P(U)/Q(U) controls

`PUController` and `QUController` attached to prosumers make the specified
injections state-dependent: the specified power vector is re-evaluated each
Newton iteration as a function of local `|V|`, and local chain-rule terms
are added to the Jacobian. Derivation, API and output semantics (`Type` vs
`Control`): [Voltage-dependent Control](voltage_dependent_control.md).

## Jacobian condition-number diagnostics

When a Newton solve stalls or diverges, the first question is whether the
Jacobian is numerically healthy at the operating point. `condestJacobian(J)`
estimates the 1-norm condition number
``\kappa_1(J) = \|J\|_1 \cdot \|J^{-1}\|_1`` with the Hager/Higham estimator
on the LU factorization, so it works on the sparse Jacobians the Newton step
factors anyway; `exact = true` switches to the exact dense 2-norm via SVD
for small systems. `reportCondition(J)` prints the estimate with the rough
number of Float64 digits lost and a verdict (well conditioned below `1e6`,
borderline below `1e10`, convergence at risk below `1e14`, numerically
singular above).

```julia
kappa = reportCondition(J)          # prints and returns the estimate
kappa = condestJacobian(J)          # estimate only
kappa = condestJacobian(J; exact = true)  # dense SVD, small nets only
kappa = condestJacobian(net)        # solver Jacobian at the net's current state
```

The net method builds the same sparse PQ/PV Jacobian the rectangular solver
factors, at the net's current voltage state: the operating point after a
converged run, the last iterate after a failed one. The diagnostics are
always on, without configuration:

| Use | Where |
|---|---|
| Result field | classic result log: one `Jacobian cond.` line with the verdict (one Hager estimate on the LU the solver factored; non-NR runs reconstruct the Jacobian); run metadata `jacobian_condition_estimate` / `_verdict` / `_line` |
| Web UI | run overview row `Jacobian condition` for rectangular NR runs (`n/a` for APSLF or DC); Diagnose runs (`run_diagnostics = true`, the Diagnose button) add the estimate and a recommendation for an ill-conditioned Jacobian to the `diagnose.log` Diagnosis section |
| Example | `examples/powerflow/exp_condition_number.jl`, including how conditioning collapses when a bus becomes nearly disconnected |
