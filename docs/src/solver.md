# Solver Guide

This page collects numerical notes for the supported internal power-flow path:
the sparse rectangular complex-state Newton–Raphson solver.

The rectangular power-flow implementation is split into focused layers:

* `src/powerflow_rectangular/rectangular_core_equations.jl` contains the core residual/equation helpers, including `mismatch_rectangular`.
* `src/powerflow_rectangular/rectangular_jacobian_builders.jl` contains the analytic rectangular Jacobian builders.
* `src/powerflow_rectangular/rectangular_newton_step.jl` contains Newton-step update and damping helpers.
* `src/powerflow_rectangular/rectangular_standalone_solver.jl` contains `run_complex_nr_rectangular`, the standalone array-level Newton driver.
* `src/powerflow_rectangular/rectangular_network_solver.jl` contains `runpf_rectangular!`, the network-integrated rectangular power-flow entry point and orchestration glue.

## Rectangular Complex-State Newton–Raphson (`powerflow_rectangular/rectangular_network_solver.jl`)

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

The complex bus powers are computed as

```math
I = Y_\mathrm{bus} V, \qquad
S = V \odot \overline{I} = V \odot \overline{Y_\mathrm{bus} V}
```

and the specified injections as

```math
S_\mathrm{spec} = P_\mathrm{spec} + j Q_\mathrm{spec}.
```

### Bus Types and Residual Definition

The function

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

`complex_newton_step_rectangular` performs one Newton step using an analytic
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

Key steps:

1. Compute currents and powers.
2. Form the residual `F(V)`.
3. Assemble the Jacobian `J = ∂F/∂x`.
4. Solve

   ```math
   J(x_k)\,\Delta x_k = -F(x_k)
   ```

5. Update the state with optional fixed damping or automatic damping.

This is the standard Newton linearization of a nonlinear root-finding problem.

### Automatic Rectangular Newton Damping

For difficult flat-start studies, enable residual-based backtracking:

```julia
runpf!(net, 60, 1e-8, 1;
    opt_flatstart = true,
    autodamp = true,
    autodamp_min = 0.05,
)
```

The rectangular solver first computes the ordinary Newton correction from
`J * Δx = -F`. With `autodamp = true`, it tests trial steps from `damp` down to
`autodamp_min` by halving the step length. The first trial that reduces the
maximum absolute mismatch is accepted. If no trial improves the mismatch, the
solver continues with the best finite conservative trial so that the next
iteration can rebuild the Jacobian and active-set state.

This is a line-search style damping strategy. It does not change the power-flow
model, bus equations, Q-limit logic, or Jacobian formulation. It only chooses the
scalar step length for the already computed rectangular Newton correction.

---

### Merit-Function Line Search

Autodamp (above) accepts the first trial step length that reduces the
**maximum absolute mismatch** — an $\infty$-norm criterion with no formal
sufficient-decrease guarantee. `power_flow.merit` adds an optional, opt-in
alternative acceptance test based on a scalar merit function, used *inside*
the same backtracking loop rather than replacing it. It is disabled by
default (`power_flow.merit.enabled = false`), which leaves the historical
autodamp behavior byte-for-byte unchanged.

#### Root-finding vs. optimization view

Newton-Raphson solves the power-flow equations as a root-finding problem,
$F(x) = 0$. The same residual can be viewed through an optimization lens by
defining the scalar merit function

```math
f(x) = \frac{1}{2} \lVert W F(x) \rVert_2^2
```

where $W$ is an optional diagonal scaling matrix (see below). Minimizing
$f$ and solving $F(x) = 0$ share the same solutions whenever $F(x) = 0$ is
attainable ($f = 0$ there), so a solver that drives $f$ down along the way
is consistent with driving $F$ to zero.

The key fact that makes the *existing* Newton direction usable as a merit
descent direction, without any extra factorization, is its gradient:

```math
\nabla f(x) = J(x)^\top W^\top W F(x)
```

With the Newton direction $\Delta x = -J(x)^{-1} F(x)$ (already computed for
the ordinary Newton step) and $W = I$:

```math
\nabla f(x)^\top \Delta x = -F(x)^\top J(x)^{-\top} J(x)^\top J(x)^\top J(x)^{-1} F(x)
```

reduces, using $J \Delta x = -F$, to

```math
\nabla f(x)^\top \Delta x = -\lVert W F(x) \rVert_2^2 < 0
```

whenever $F(x) \ne 0$. So the Newton direction is always a **descent
direction** for $f$, and the directional derivative is available directly
from the already-computed mismatch vector — no extra Jacobian-vector product
or factorization is needed to evaluate it.

#### Armijo sufficient-decrease acceptance

For each backtracking trial step length $\lambda$ (tried in the same order
as autodamp: `damp`, then halved down to `autodamp_min`), the merit line
search accepts the first $\lambda$ satisfying the Armijo condition:

```math
f(x + \lambda \Delta x) \le f(x) + c_1 \, \lambda \, \nabla f(x)^\top \Delta x
```

`power_flow.merit.armijo_c1` is $c_1 \in (0, 0.5)$, the sufficient-decrease
constant. Smaller $c_1$ accepts almost any step that decreases $f$ at all;
larger $c_1$ requires the decrease to be a larger fraction of the predicted
linear decrease, rejecting more trials and backtracking further. This is the
standard Armijo (backtracking) line-search condition — see Nocedal & Wright,
*Numerical Optimization*, the chapter on line-search methods, for the general
theory.

Compared to the existing $\infty$-norm autodamp criterion, Armijo acceptance
gives a **globally aware** descent guarantee on the aggregate residual energy
$f$, whereas the max-mismatch criterion only asks whether the single
worst-offending bus equation improved. A step can reduce the worst-bus
mismatch while increasing the overall residual energy elsewhere, or vice
versa; the two criteria can therefore accept different trial step lengths on
the same iteration.

If no trial step length satisfies Armijo, `power_flow.merit.fallback_max_mismatch`
selects the fallback behavior:

* `true` (default): fall back to the classic max-mismatch criterion for this
  iteration (first improving trial, else the conservative best-finite trial —
  identical to disabling merit for that iteration);
* `false`: skip straight to the conservative best-finite-trial fallback,
  without checking whether any trial merely improved the max mismatch.

#### Limits

* Minimizing $f$ is not equivalent to solving $F(x) = 0$: $f$ can have local
  minima with $F(x) \ne 0$ (e.g. saddle-like residual configurations), so
  satisfying Armijo at every iteration does not by itself guarantee
  convergence to a root.
* The merit criterion does not influence *which* numerical solution branch
  (e.g. a high-voltage vs. a low-voltage operating point) the solver
  converges to; it only accepts or rejects step lengths along the path
  Newton-Raphson already takes.
* A PV/PQ active-set switch (Q-limit handling) changes what the residual
  vector's entries *mean* (a $\Delta Q$ entry becomes a $\Delta V$ entry, or
  vice versa) for the buses that switched. Comparing $f$ across such a switch
  is not well-defined, so the merit comparison is skipped for the Newton
  iteration in which a switch occurred; that iteration falls back to the
  classic max-mismatch criterion and is logged with
  `accept_reason = active_set_skip` (see
  [Merit-function line search options](@ref) in the configuration reference
  for the diagnostic log format).

#### Residual scaling

$P$, $Q$, and voltage-setpoint ($V$) residual entries can differ in
magnitude by orders of magnitude depending on the per-unit base and grid
mix (e.g. a large-$P$/small-$Q$ transmission grid, or PV buses with tight
voltage tolerances). Because $f$ sums squared residuals, an unscaled merit
value can be dominated by whichever equation type happens to have the
largest natural magnitude, which weakens the sufficient-decrease guarantee
for the other equation types. `power_flow.merit.scale_p`,
`power_flow.merit.scale_q`, and `power_flow.merit.scale_v` apply positive
diagonal weights ($W$ above) per equation type — `ΔP_i` entries always use
`scale_p`; the second equation per non-slack bus uses `scale_q` for PQ buses
and `scale_v` for PV buses. Leave all three at `1.0` (the default) unless
diagnostics show one residual type dominating the merit value in a specific
case.

---

### Trust-Region Step Control

`power_flow.trust_region` is a second, alternative Newton step-control
mechanism to autodamp: a **scaled-Newton trust region**. It is disabled by
default and mutually exclusive with `power_flow.autodamp` — both mechanisms
decide the same thing (how far to step along the Newton direction) by
different rules, and layering them is undefined, so enabling both is a
configuration error.

#### Step-construction modes: `scaled` and `dogleg`

Classical trust-region methods restrict the Newton correction to a ball of
radius $\Delta$ around the current iterate, typically choosing the trial step
via a *dogleg* or *Steihaug-CG* interpolation between the steepest-descent
and full Newton directions. `power_flow.trust_region.step_mode` selects
between two trial-step constructions; both reuse the same acceptance rule
(see "Acceptance by merit decrease" below).

##### `step_mode = :scaled` (default)

The trial step is always the full Newton direction $\Delta x$, rescaled down
when it exceeds the radius:

```math
\Delta x_{\text{scaled}} =
\begin{cases}
\Delta x & \lVert \Delta x \rVert \le \Delta \\
\Delta \dfrac{\Delta x}{\lVert \Delta x \rVert} & \lVert \Delta x \rVert > \Delta
\end{cases}
```

This keeps the direction always the analytic Newton direction — identical to
the autodamp/fixed-damping paths — and only the step *length* is controlled,
consistent with how autodamp itself only chooses a scalar step length. This
is the byte-for-byte original implementation and remains the default.

##### `step_mode = :dogleg`

`scaled` steps degrade poorly when the Newton direction itself stops being a
useful descent direction (e.g. a very ill-conditioned Jacobian): every
rescale keeps pointing the same bad way, so the solver either keeps taking
uphill steps (with autodamp's conservative fallback) or repeatedly shrinks
the radius toward collapse. The dogleg step mode adds a second endpoint —
the steepest-descent (Cauchy) step on the local quadratic model of the merit
function — and interpolates along the Cauchy-to-Newton path.

**Merit-function gradient.** With $f(x) = \frac{1}{2}\lVert F(x) \rVert_2^2$
(unweighted, $W = I$, the same convention `m(x)` already uses for both step
modes), the gradient is

```math
g(x) = J(x)^{\mathsf{T}} F(x)
```

a single transposed sparse matrix-vector product; no additional
factorization. This is the first use of a transposed Jacobian product in the
solver.

**Cauchy point.** The steepest-descent minimizer of the local quadratic model
$m(p) = f + g^{\mathsf{T}}p + \frac{1}{2}p^{\mathsf{T}}(J^{\mathsf{T}}J)p$
along $-g$ is

```math
\alpha^\ast = \frac{\lVert g \rVert_2^2}{\lVert Jg \rVert_2^2}
\qquad p_C = -\alpha^\ast g
```

one more sparse matrix-vector product ($Jg$), clipped to $\lVert p_C \rVert
\le \Delta$ when used as the trial step.

**Dogleg path.** With the Newton step $p_N = \Delta x$ (already computed) and
radius $\Delta$:

1. $\lVert p_N \rVert \le \Delta$: take the full Newton step
   (`accept_reason = :dogleg_newton`, identical to today's `scaled` behavior
   in the region where the radius doesn't bind).
2. $\lVert p_C \rVert \ge \Delta$: take $p_C$ rescaled to length $\Delta$
   (`accept_reason = :dogleg_cauchy`, pure steepest descent).
3. Otherwise: interpolate $p(\tau) = p_C + \tau (p_N - p_C)$, $\tau \in
   [0,1]$, choosing $\tau$ so that $\lVert p(\tau) \rVert = \Delta$
   (`accept_reason = :dogleg_interp`). This reduces to a scalar quadratic
   $a\tau^2 + b\tau + c = 0$ with
   $a = \lVert d \rVert_2^2$, $b = 2\, p_C^{\mathsf{T}} d$,
   $c = \lVert p_C \rVert_2^2 - \Delta^2$, $d = p_N - p_C$, solved
   in closed form (no inner solver).

Along the dogleg path, $\lVert p(\tau) \rVert$ is monotonically increasing
and the model value $m(p(\tau))$ is monotonically decreasing **provided**
$J^{\mathsf{T}}J$ is positive definite (Nocedal & Wright, *Numerical
Optimization*, trust-region chapter). Near a singular Jacobian —
precisely the regime dogleg is meant to help with — $J^{\mathsf{T}}J$ may
only be positive *semidefinite*; the interpolated path is still formally
defined, but its Newton endpoint can be a poor local model there. Dogleg is
therefore a partial answer: it buys graceful degradation (Cauchy-direction
progress) for *transiently* ill-conditioned Jacobians, not a global-
convergence guarantee. Two escalation options exist in the literature —
Levenberg–Marquardt ($(J^{\mathsf{T}}J + \lambda I) \Delta x = -J^{\mathsf{T}}F$,
a new factorization per $\lambda$-change) and Steihaug-CG (matrix-free,
terminates on negative curvature or the radius boundary, most robust and
most implementation effort) — both are **out of scope** for this mode and
are not implemented.

**Active-set switches.** A PV/PQ switch changes what the residual entries
mean, so a gradient/Cauchy point computed pre-switch is invalid post-switch.
Mirroring the merit-function line search's `active_set_skip` policy: when a
switch happened this Newton iteration, the dogleg comparison is skipped —
a single `scaled`-style trial is taken at the current radius and accepted
unconditionally (`accept_reason = :active_set_skip`), and the radius is left
unchanged for that iteration.

#### Acceptance by merit decrease

Step acceptance reuses the same merit function as the
[Merit-Function Line Search](@ref), $f(x) = \frac{1}{2}\lVert F(x) \rVert_2^2$
(unweighted, $W = I$) — no second merit computation is implemented. For a
trial step $p$ (the `scaled`-rescaled Newton step, or the dogleg-selected
step), the **actual reduction** is

```math
\text{ared} = f(x) - f(x + p)
```

and the **predicted reduction** comes from the local linear model of $F$
around $x$:

```math
\text{pred} = f(x) - \frac{1}{2} \lVert F(x) + J(x)\, p \rVert_2^2
```

The acceptance ratio $\rho = \text{ared} / \max(\text{pred}, \varepsilon)$
measures how well the linear model predicted the actual improvement. A trial
is accepted when $\rho \ge \eta_{\text{accept}}$
(`power_flow.trust_region.eta_accept`); note this uses the already-factored
Jacobian's matrix-vector product for `pred`, not a second Jacobian build or
linear solve.

#### Radius update

The radius $\Delta$ persists across Newton iterations (unlike autodamp's
per-iteration `damp` restart) and adapts from $\rho$:

* **Accepted** ($\rho \ge \eta_{\text{accept}}$): the step is taken. If the
  step also hit the radius boundary — $\lVert \Delta x \rVert > \Delta$ in
  `scaled` mode, or `accept_reason` is `:dogleg_cauchy`/`:dogleg_interp` in
  `dogleg` mode — and $\rho \ge$ `expand_threshold`, the radius expands:
  $\Delta \leftarrow \min(\Delta \cdot \texttt{expand\_factor}, \Delta_{\max})$.
  Otherwise the radius is unchanged.
* **Rejected** ($\rho < \eta_{\text{accept}}$): the radius shrinks,
  $\Delta \leftarrow \Delta \cdot \texttt{shrink\_factor}$, and the *same*
  Newton iteration retries with the smaller radius — no new Jacobian build or
  linear solve, only a re-evaluation of the step construction (rescale, or
  dogleg branch re-selection at the smaller radius) and a fresh mismatch
  evaluation, mirroring how autodamp reuses one Newton correction across its
  backtracking trials.
* **Collapsed**: if $\Delta$ falls below `power_flow.trust_region.min_radius`
  without an accepted step, the solver declares non-convergence with reason
  `:trust_region_collapsed` rather than looping indefinitely.

#### Limits

* Neither mode is a global-convergence guarantee or a wrong-branch selector;
  both are acceptance/step-length controls layered on the existing analytic
  Newton direction and Jacobian, like the merit-function line search.
* `dogleg` does not fix a bad start. The one documented real-world case that
  motivated this mode (`case_SyntheticUSA.m`, island 1, `classic`/`classic`
  start; see the case matrix) was fully resolved by switching start mode, not
  by step control — dogleg would only have degraded gracefully there instead
  of diverging, not converged. If the start is fundamentally wrong for the
  case, fix the start mode first.
* `dogleg`'s monotonicity guarantee assumes $J^{\mathsf{T}}J \succ 0$; near a
  singular Jacobian this weakens (see above). Levenberg–Marquardt and
  Steihaug-CG, which handle that regime more robustly, are not implemented.
* `dogleg`'s Cauchy (steepest-descent) steps can be slow on badly scaled
  problems — this is inherent to the method, not an implementation defect.
  There is no per-equation weighting knob for the dogleg gradient (it always
  uses $W = I$, matching `m(x)`); if a case needs weighted descent, that is
  currently only available via `power_flow.merit.scale_p/q/v` on the
  (mutually exclusive) merit-function line search.
* Mutually exclusive with `autodamp`; there is currently no combined
  trust-region-with-backtracking-fallback mode.

See Nocedal & Wright, *Numerical Optimization*, the chapters on trust-region
methods (dogleg, Cauchy point, Steihaug-CG) and line-search methods (Armijo),
for the general theory this section adapts, referenced by name only.

---

### Start Projection for Difficult Seeds

#### DC-angle flat-start background

A conventional AC flat start initializes most bus voltages near $1.0\,\mathrm{pu}$
and all non-slack voltage angles near $0^\circ$. That seed is simple and
reproducible, but it ignores the active-power flow pattern that is already
encoded in the network topology, branch reactances, transformer phase shifts,
and specified injections. In large or heavily phase-shifted MATPOWER cases, the
first Newton step may therefore start far away from the physically relevant
angle branch.

The DC-angle seed uses the active-power part of the power-flow model as a
linearized predictor for voltage angles. Under the usual high-voltage,
small-angle assumptions, voltage magnitudes are held near nominal values,
resistance and reactive-power coupling are neglected, and active-power transfer
across a branch is approximated by the angle difference divided by branch
reactance. This gives a sparse linear system of the form

```math
B'\theta = P
```

where $P$ is the net active-power injection vector and $B'$ is assembled from
the branch susceptance structure. The slack angle fixes the reference, and the
resulting angles are clipped by `start_projection_dc_angle_limit_deg` before
being used as a Newton seed. This is not a replacement for the AC solve; it is
only an initialization step that preserves the full AC equations, Q-limit logic,
and rectangular Newton formulation used by the main run.

Blended starts combine this DC-angle predictor with the stored MATPOWER `VM`/`VA`
data or the raw flat start. The projection scans the requested candidates and
keeps the one with the smallest rectangular mismatch, which helps avoid wrong
low-voltage or wrong-angle branches without changing the final convergence
criteria.

The rectangular solver can project the initial voltage before Newton iterations:

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

The projection is inspired by the other workflows. It sanitizes the raw seed,
optionally builds a DC-angle approximation from active-power injections and the
Y-bus off-diagonal susceptances, and optionally scans convex angle/magnitude
blends between the raw seed and the DC start. Sparlectra picks the candidate with
the lowest rectangular mismatch and then runs the ordinary solver.

The same projection options are also available through `buildPfModel` and
`runpf_external!`, so external solvers can receive the projected
`model.V0` without changing their own solver implementation.

---


## Solver-Specific Interaction with Power Limits

When a PV bus hits a Q-limit, the solver does **not** merely clamp a reported
reactive power value. It changes the equation type for that bus:

* from `(ΔP, ΔV)` at a PV bus
* to `(ΔP, ΔQ)` at a PQ bus.

In the rectangular solver this is represented through the `bus_types` vector.

The start of PV→PQ switching can be controlled for difficult cases:

```julia
runpf!(net, 60, 1e-8, 1;
    qlimit_start_iter = 4,
    qlimit_start_mode = :iteration,
)
```

Use `qlimit_start_mode = :auto` to wait until the PV reactive-power
requests have stabilized. The threshold is `qlimit_auto_q_delta_pu` in p.u.
The hybrid `:iteration_or_auto` mode starts when either the configured iteration
or the reactive-power stabilization criterion is reached.

For the operational / data-model side of Q-limit handling, see
[Powerlimits Guide](powerlimits.md).

---

## Voltage-dependent P(U)/Q(U) controls

The rectangular solver also supports state-dependent specified injections via
`PUController` and `QUController` attached to prosumers. In that case, the
specified power vector is re-evaluated each Newton iteration as a function of
local `|V|`, and local chain-rule terms are added to the Jacobian.

For full derivation, formulas, API usage, and output semantics (`Type` vs
`Control`), see [Voltage-dependent Control](voltage_dependent_control.md).
