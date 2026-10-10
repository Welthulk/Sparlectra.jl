# Solver Guide

The internal power-flow path is the sparse rectangular complex-state
Newton-Raphson solver: `runpf_rectangular!` (network level,
`src/powerflow_rectangular/rectangular_network_solver.jl`) over
`run_complex_nr_rectangular` (array level), with the residual helpers
(`mismatch_rectangular`), the analytic Jacobian builders and the step
control (damping, line search, trust region) in the same directory.

## Rectangular Complex-State Newton-Raphson

### Motivation and State Vector

The state variables are the complex bus voltages

```math
V_k = V_{r,k} + j V_{i,k}
```

(real part $V_{r,k}$ and imaginary part $V_{i,k}$ of bus $k$), not magnitudes
and angles. The state vector is

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

builds the residual vector `F(V)`:

* PQ bus $i \ne \text{slack}$:

  ```math
  \Delta P_i = \Re(S_{\text{calc},i}) - \Re(S_{\text{spec},i})
  \\
  \Delta Q_i = \Im(S_{\text{calc},i}) - \Im(S_{\text{spec},i})
  ```

* PV bus $i \ne \text{slack}$:

  ```math
  \Delta P_i = \Re(S_{\text{calc},i}) - \Re(S_{\text{spec},i})
  \\
  \Delta V_i = |V_i| - V_\text{set,i}
  ```

The residual has `2(n-1)` entries, the state dimension.

### Analytic Rectangular Newton Step

`complex_newton_step_rectangular` performs one Newton step with the analytic
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

It computes currents and powers, forms `F(V)`, assembles `J = ∂F/∂x`, solves

```math
J(x_k)\,\Delta x_k = -F(x_k)
```

and applies the step with fixed or automatic damping.

### [Newton Update: Rectangular or Polar](@id newton_update)

`power_flow.newton_update` selects how the step `Δx = (ΔV_r, ΔV_i)` is
applied to every non-slack bus:

| value | update | property |
|---|---|---|
| `rectangular` | `V ← V + ΔV` | A pure rotation `ΔV/V = jb` inflates the magnitude by `sqrt(1 + b²)`. |
| `polar` (default) | `V ← V (1 + a) e^{jb}` with `ΔV/V = a + jb` | Magnitude times `1 + a`, angle plus `b`: exact for rotations and for scalings. |

Same Jacobian, same solution, different iterates; from a flat start on a
large transmission case the difference decides convergence. Damping
(`autodamp`, the merit line search, the trust region) scales the step
before it is applied and works with either update. Example and
measurements: [Newton update](newton_update.md); recipes per case:
[Start strategies by case](start_strategies.md).

### Automatic Rectangular Newton Damping

`autodamp = true` backtracks on the residual: trial step lengths run from
`damp` down to `autodamp_min` by halving, and the first trial that reduces
the maximum absolute mismatch is accepted. If none does, the solver takes
the best finite trial and rebuilds Jacobian and active set in the next
iteration. Only the scalar step length changes.

```julia
runpf!(net, 60, 1e-8, 1;
    opt_flatstart = true,
    autodamp = true,
    autodamp_min = 0.05,
)
```

### Merit-Function Line Search

Autodamp's test (a smaller maximum mismatch) has no sufficient-decrease
guarantee. `power_flow.merit` adds an opt-in Armijo test on a scalar merit
function inside the same backtracking loop.

| Use | Where |
|---|---|
| Config key | `power_flow.merit.enabled` (default `false`); `armijo_c1`, `fallback_max_mismatch`, `scale_p` / `scale_q` / `scale_v` under the same prefix |
| Log field | `accept_reason` per iteration ([Merit-function line search options](@ref pf-merit)) |

#### Root-finding vs. optimization view

The residual defines the merit function

```math
f(x) = \frac{1}{2} \lVert W F(x) \rVert_2^2
```

with $W$ the diagonal matrix of the weights `scale_p`/`scale_q`/`scale_v`;
where $F(x) = 0$ is attainable, the minimizers of $f$ are the roots of $F$.
Its gradient is

```math
\nabla f(x) = J(x)^\top W^\top W F(x)
```

and with $J \Delta x = -F$ the directional derivative

```math
\nabla f(x)^\top \Delta x = F(x)^\top W^\top W \, J(x) \, \Delta x = -\lVert W F(x) \rVert_2^2 < 0
```

is negative whenever $F(x) \ne 0$: the Newton direction descends $f$ at no
extra cost.

#### Armijo sufficient-decrease acceptance

For the trial step lengths $\lambda$ (`damp`, then halved down to
`autodamp_min`) the first one with

```math
f(x + \lambda \Delta x) \le f(x) + c_1 \, \lambda \, \nabla f(x)^\top \Delta x
```

is accepted. `power_flow.merit.armijo_c1` is $c_1 \in (0, 0.5)$: a small
$c_1$ accepts almost any decrease of $f$, a large one demands a larger share
of the predicted decrease and backtracks further (Nocedal and Wright,
*Numerical Optimization*, line-search chapter). Armijo guarantees descent of
the aggregate residual energy, the max-mismatch test only an improvement of
the worst bus equation; a step can do one without the other, so the two
tests can accept different step lengths.

If no trial satisfies Armijo, `power_flow.merit.fallback_max_mismatch`
decides: `true` (default) uses the max-mismatch criterion for this iteration
(first improving trial, else the best finite trial), `false` goes straight
to the best finite trial.

#### Limits

* $f$ can have local minima with $F(x) \ne 0$; Armijo at every iteration
  does not guarantee a root.
* The test chooses step lengths along Newton's path, not the solution branch
  (a high- versus a low-voltage operating point).
* A PV/PQ switch changes what the residual entries mean ($\Delta Q$ becomes
  $\Delta V$ or vice versa), so $f$ is not comparable across it; that
  iteration uses the max-mismatch criterion and logs
  `accept_reason = active_set_skip`.

#### Residual scaling

$P$, $Q$ and $V$ residual entries can differ by orders of magnitude
(per-unit base, grid mix), and one type then dominates $f$.
`power_flow.merit.scale_p`, `scale_q` and `scale_v` are the weights $W$:
`ΔP_i` entries use `scale_p`, the second equation of a non-slack bus
`scale_q` (PQ) or `scale_v` (PV). Leave all three at `1.0` unless
diagnostics show one type dominating.

### Trust-Region Step Control

`power_flow.trust_region` caps the Newton step at an adaptive radius and
accepts trials by merit decrease. It is off by default and mutually
exclusive with `power_flow.autodamp` (both decide the step length; enabling
both is a configuration error).

| Use | Where |
|---|---|
| Config key | `power_flow.trust_region.enabled` (default `false`); `step_mode` (`:scaled`, `:dogleg`), `eta_accept`, `expand_threshold`, `expand_factor`, `shrink_factor`, `min_radius` under the same prefix ([Trust-region step control options](@ref pf-trust-region)) |
| Result field | `accept_reason` (`:dogleg_newton`, `:dogleg_cauchy`, `:dogleg_interp`, `:active_set_skip`); non-convergence reason `:trust_region_collapsed` |

#### Step-construction modes: `scaled` and `dogleg`

`power_flow.trust_region.step_mode` selects the trial step inside the radius
$\Delta$; both modes share the acceptance rule below.

##### `step_mode = :scaled` (default)

The Newton direction $\Delta x$, rescaled when it exceeds the radius:

```math
\Delta x_{\text{scaled}} =
\begin{cases}
\Delta x & \lVert \Delta x \rVert \le \Delta \\
\Delta \dfrac{\Delta x}{\lVert \Delta x \rVert} & \lVert \Delta x \rVert > \Delta
\end{cases}
```

Only the length is controlled, as with autodamp.

##### `step_mode = :dogleg`

`scaled` degrades when the Newton direction is no descent direction (an
ill-conditioned Jacobian): every rescale points the same way, uphill or
into radius collapse. Dogleg interpolates between the steepest-descent
(Cauchy) step of the local quadratic model and the Newton step.

**Merit-function gradient.** With $f(x) = \frac{1}{2}\lVert F(x) \rVert_2^2$
(unweighted, $W = I$, for both step modes):

```math
g(x) = J(x)^{\mathsf{T}} F(x)
```

one transposed sparse matrix-vector product.

**Cauchy point.** The minimizer of the local quadratic model
$m(p) = f + g^{\mathsf{T}}p + \frac{1}{2}p^{\mathsf{T}}(J^{\mathsf{T}}J)p$
along $-g$:

```math
\alpha^\ast = \frac{\lVert g \rVert_2^2}{\lVert Jg \rVert_2^2}
\qquad p_C = -\alpha^\ast g
```

one more sparse product ($Jg$), clipped to $\lVert p_C \rVert \le \Delta$ as
a trial step.

**Dogleg path.** With $p_N = \Delta x$ and radius $\Delta$:

1. $\lVert p_N \rVert \le \Delta$: the full Newton step
   (`accept_reason = :dogleg_newton`, as `scaled` where the radius does not
   bind).
2. $\lVert p_C \rVert \ge \Delta$: $p_C$ rescaled to length $\Delta$
   (`accept_reason = :dogleg_cauchy`, pure steepest descent).
3. Otherwise $p(\tau) = p_C + \tau (p_N - p_C)$, $\tau \in [0,1]$, with
   $\lVert p(\tau) \rVert = \Delta$ (`accept_reason = :dogleg_interp`): the
   scalar quadratic $a\tau^2 + b\tau + c = 0$ with
   $a = \lVert d \rVert_2^2$, $b = 2\, p_C^{\mathsf{T}} d$,
   $c = \lVert p_C \rVert_2^2 - \Delta^2$, $d = p_N - p_C$, solved in closed
   form.

Along the path $\lVert p(\tau) \rVert$ increases and $m(p(\tau))$ decreases
monotonically provided $J^{\mathsf{T}}J$ is positive definite (Nocedal and
Wright, trust-region chapter); near a singular Jacobian it may only be
semidefinite, and the Newton endpoint can be a poor model. Dogleg gives
Cauchy-direction progress through a transiently ill-conditioned Jacobian,
not a global-convergence guarantee. Levenberg-Marquardt and Steihaug-CG are
not implemented.

**Active-set switches.** A PV/PQ switch invalidates a gradient or Cauchy
point computed before it, so the dogleg comparison is skipped in that
iteration: one `scaled`-style trial at the current radius is accepted
unconditionally (`accept_reason = :active_set_skip`) and the radius stays.

#### Acceptance by merit decrease

For a trial step $p$ the actual reduction of
$f(x) = \frac{1}{2}\lVert F(x) \rVert_2^2$ ($W = I$) is

```math
\text{ared} = f(x) - f(x + p)
```

and the predicted reduction from the linear model of $F$ is

```math
\text{pred} = f(x) - \frac{1}{2} \lVert F(x) + J(x)\, p \rVert_2^2
```

A trial is accepted when
$\rho = \text{ared} / \max(\text{pred}, \varepsilon) \ge \eta_{\text{accept}}$
(`power_flow.trust_region.eta_accept`); `pred` is one matrix-vector product
with the factored Jacobian.

#### Radius update

The radius persists across Newton iterations (autodamp restarts `damp`
every iteration):

* Accepted: the step is taken; if it hit the boundary
  ($\lVert \Delta x \rVert > \Delta$ in `scaled` mode, `:dogleg_cauchy` or
  `:dogleg_interp` in `dogleg` mode) and $\rho \ge$ `expand_threshold`, the
  radius becomes $\min(\Delta \cdot \texttt{expand\_factor}, \Delta_{\max})$.
* Rejected: $\Delta \leftarrow \Delta \cdot \texttt{shrink\_factor}$, and the
  same iteration retries with the smaller radius without a new Jacobian or
  linear solve.
* Collapsed: below `power_flow.trust_region.min_radius` without an accepted
  step the solver stops with reason `:trust_region_collapsed`.

#### Limits

Neither mode guarantees global convergence or selects the solution branch;
both control the step length along the Newton direction, and `dogleg` does
not repair a bad start. Cauchy steps can be slow on badly scaled problems.
The dogleg gradient always uses $W = I$; weighted descent exists only
through `power_flow.merit.scale_p/q/v` on the (mutually exclusive) merit
line search.

### Start Projection for Difficult Seeds

#### DC-angle flat-start background

A flat start (magnitudes near $1.0\,\mathrm{pu}$, non-slack angles
$0^\circ$) ignores the active-power pattern of topology, reactances, phase
shifts and injections; on a large or heavily phase-shifted case the first
Newton step can start far from the physical angle branch. The DC-angle seed
takes the angles from the linearized active-power model (magnitudes
nominal, resistance and reactive coupling neglected, branch flow
proportional to the angle difference over the reactance):

```math
B'\theta = P
```

with $P$ the net active-power injections (with a DC model present
$P - P_{\mathrm{inj}}$, the phase-shifter injections the DC power flow
subtracts) and $B'$ assembled from the branch susceptances. The slack angle
is the reference; the angles are clipped by
`start_projection_dc_angle_limit_deg`.

Blended starts combine the DC-angle predictor with the stored MATPOWER
`VM`/`VA` data or the flat start. On a flat start the projection also
measures the ratio profile (`start_projection_ratio_profile`, key
`power_flow.start_mode.ratio_profile`): the flat magnitudes of the PQ buses
multiplied by the product of the off-nominal transformer ratios on the path
from the reference (the ratio sits on the from side: a step from the from
bus to the to bus divides by $\lvert t \rvert$, the reverse step multiplies
by it), so that no transformer carries a circulating current at the start.
The projection keeps the candidate with the smallest residual 2-norm when it
is at least 10 percent below the raw seed's, which helps avoid wrong
low-voltage or wrong-angle branches; the convergence criteria are unchanged.

```julia
runpf!(net, 60, 1e-8, 1;
    opt_flatstart = true,
    start_projection = true,
    start_projection_try_dc_start = true,
    start_projection_try_blend_scan = true,
    start_projection_blend_lambdas = [0.25, 0.5, 0.75],
    start_projection_dc_angle_limit_deg = 60.0,
    start_projection_ratio_profile = true,
)
```

The same keywords exist on `buildPfModel` and `runpf_external!`, so an
external solver receives the projected `model.V0`.

## Distributed Active-Power Slack

### Why a single slack bus is a modeling artifact

The classical formulation lets the reference bus absorb the imbalance plus
the losses. No single machine does this: primary control raises many
generators together, each by its droop share. Concentrating that response
on one bus distorts the flows around it; in large imported cases the
reference machine absorbs tens of MW that real dispatch spreads over an
area. The distributed slack shares the imbalance among **participants** by
normalized **participation factors** $\alpha_i$ with $\sum_i \alpha_i = 1$.

### Formulation

One scalar unknown $\lambda_P$ (the island's total imbalance in p.u.) joins
the state. Every participant bus $i$ (the reference bus included if it
participates) gets the active-power residual

```math
r_{P,i} = P_{\text{calc},i}(V) - P_{\text{spec},i} - \alpha_i\,\lambda_P
```

so its solved injection is its schedule plus its share; non-participants
keep their residual ($\alpha_i = 0$).

The system stays square through a role separation at the reference bus: its
voltage magnitude and angle stay fixed, its active power stops being free.
Its P residual becomes an ordinary equation, appended as the last row (the
interleaved `[ΔP_i, ΔQ/ΔV_i]` layout of the `2(n-1)` voltage rows is
unchanged): one new unknown, one new equation (REF-P), one Jacobian column
($-\alpha_i$ at the participant P rows) and one row (the network derivatives
$\partial P_\text{ref}/\partial V_{r,j}, \partial P_\text{ref}/\partial V_{i,j}$,
plus $-\alpha_\text{ref}$ in the $\lambda_P$ column).

With all weight on the reference bus ($\alpha_\text{ref} = 1$) the appended
equation reads $P_{\text{calc,ref}} - P_{\text{spec,ref}} - \lambda_P = 0$
and decouples from the voltage solution: the classical result is reproduced
exactly, with $\lambda_P$ reporting the classically absorbed power (a
regression test).

### When and how it applies

| Use | Where |
|---|---|
| Config key | `power_flow.distributed_slack.enabled` (default `false`); `p_mode`, `fallback`, `respect_p_limits` under the same prefix ([Distributed active-power slack](@ref pf-distributed-slack)) |
| Result field | $\lambda_P$ in the per-island solver statuses; P-limit warnings counted in the run metadata |

* Imported participation factors (MATPOWER `APF`, CGMES
  `GeneratingUnit.normalPF`) are stored on the `ProSumer` but never activate
  the feature; the disabled path is bit-identical to the classical solver.
* Participants are the generator-type prosumers at the island's REF and PV
  buses at solve start; fixed injections at PQ buses (HVDC converter
  injections, boundary equivalents) never participate. Participation is
  frozen per solve: a participant that hits its Q limit and becomes PQ keeps
  its $\alpha_i$.
* Every island solves for its own $\lambda_P$; an island without a valid
  participant aborts (`fallback: error`) or solves classically with a
  warning (`fallback: ref_only`).
* Autodamp, merit line search and trust region stage the trial multiplier
  $\lambda_P + a\,\delta\lambda_P$ with the trial state $V + a\,\delta V$;
  in the weighted merit norm the REF-P row is scaled with `merit.scale_p`.
* P limits are advisory: with `respect_p_limits = true` a participant whose
  corrected output leaves `[minP, maxP]` produces a warning, counted in the
  metadata; there is no clamp-and-resolve round.
* The augmented system requires the sparse Jacobian path (the default); the
  dense builder rejects it.
* Against a single-slack reference state (a CGMES SV profile) the branch
  flows deviate by the distributed correction while the voltages stay put;
  the comparison measures the slack distribution, not solver error.

## Linear solver backends

The Newton step `J · δx = -F` has two sparse backends.

| Use | Where |
|---|---|
| Config key | `power_flow.linear_solver`: `umfpack_reuse` (default) or `umfpack`; independent of `power_flow.solver`, which selects the power-flow method (rectangular / apslf / dc) |
| Result field | backend plus the counters `linear_solver_analyze_count`, `linear_solver_refactor_count`, `linear_solver_fallback_count` in the solver status and the performance profile, per island on the island path |

* `umfpack`: the sparse direct solve (`J \ F` through `solve_sparse_system`)
  with sparse-QR and small-system SVD fallbacks; every iteration pays
  symbolic analysis plus numeric factorization.
* `umfpack_reuse`: the same factorization with the symbolic analysis of the
  first iteration kept and reused via `lu!(F, J)`. The pattern is constant
  across iterations unless the active set changes, so the analysis is paid
  once per active set; usually the fastest choice on large cases.

KLU is not a value of `power_flow.linear_solver`: for a single solve it is
not faster than `umfpack_reuse`, and a KLU factorization shared across
threads gives silently wrong results. Power mode uses KLU when the extension
is loaded (`using KLU`), one factorization per network, never shared;
`power_flow.power_mode_lu` chooses per network, see
[Power-mode LU](@ref power_mode_lu).

The reuse backend re-runs analyze plus factor after a PV/PQ switch (the
pattern changes) and whenever the stored `colptr`/`rowval` fingerprint of
the analyzed pattern differs from the current Jacobian; the Jacobian is
assembled with a value-independent pattern, so the pattern depends only on
the Ybus structure and the bus types. A factorization error (structure
mismatch, bad pivots, singularity) resets the context and delegates that
solve to the `umfpack` chain, so the `:singular_newton_step` handling and
the QR/SVD fallbacks keep working. The context covers the main Newton solve
only (the trust-region internals keep their path), lives for one
`runpf_rectangular!` call (one island on the island path) and is never
shared across islands or threads.

### [Dishonest Newton](@id dishonest_newton_solver)

`power_flow.jacobian_reuse` (off by default) keeps the LU factors of the
Jacobian for the next step while the iteration converges fast: after a step
that cut the maximum mismatch by at least `jacobian_reuse_min_reduction`
(default 10), the next step solves the new mismatch with the previous
factors, without refilling or refactorising the Jacobian. Otherwise the step
refills and refactorises as usual. The rules:

1. Never in the first two steps.
2. Never in a step in which the Q-limit active set converted or released a
   bus (outer-loop controllers act between solves; each solve starts with
   honest steps).
3. A reused step that does not cut the mismatch by the factor is kept, and
   the next step refactorises.
4. At most `jacobian_reuse_max_steps` (default 3) reused steps in a row.
5. Tolerance and convergence test are unchanged; the final state equals the
   honest Newton solution within the tolerance.
6. Works with and without power mode, on `umfpack_reuse` and on KLU (power
   mode), with the polar and the rectangular update. With
   `linear_solver: umfpack` there is no factorization to keep, and with the
   trust region every step needs the current Jacobian; the run log then
   says that the switch was ignored.

A reused step converges linearly instead of quadratically: more steps, fewer
factorizations. The switch pays where the factorization dominates the step
and many steps are taken; warm-started N-1 and scenario solves get slower
with it. Measurements: [Dishonest Newton](@ref dishonest_newton) on the
performance page. The solver status and the run log carry the reused steps
and the refactorisations.

## Solver-Specific Interaction with Power Limits

A PV bus that hits a Q limit changes its equation type from `(ΔP, ΔV)` to
`(ΔP, ΔQ)` through the `bus_types` vector; the reported reactive power is
not merely clamped. When the switching starts is a keyword:

```julia
runpf!(net, 60, 1e-8, 1;
    qlimit_start_iter = 4,
    qlimit_start_mode = :iteration,
)
```

`qlimit_start_mode = :auto` waits until the PV reactive-power requests have
settled (threshold `qlimit_auto_q_delta_pu` in p.u.); `:iteration_or_auto`
starts at whichever comes first. Mechanics and data model:
[Powerlimits Guide](powerlimits.md); keys:
[Q-limit options and guard](@ref pf-qlimits).

## Voltage-dependent P(U)/Q(U) controls

`PUController` and `QUController` attached to prosumers make the specified
injections state-dependent: the specified power vector is re-evaluated each
Newton iteration as a function of local `|V|`, and local chain-rule terms
are added to the Jacobian. Derivation, API and output semantics (`Type` vs
`Control`): [Voltage-dependent Control](voltage_dependent_control.md).

## Jacobian condition-number diagnostics

`condestJacobian(J)` estimates the 1-norm condition number
``\kappa_1(J) = \|J\|_1 \cdot \|J^{-1}\|_1`` with the Hager/Higham estimator
on the LU factorization the Newton step factors anyway; `exact = true`
computes the dense 2-norm via SVD for small systems. `reportCondition(J)`
prints the estimate, the attainable relative accuracy
$\kappa_1 \cdot \varepsilon_{\mathrm{Float64}}$ and a verdict (well
conditioned below `1e6`, borderline below `1e10`, convergence at risk below
`1e14`, numerically singular above). The net method builds the solver's
Jacobian at the net's current voltage state: the operating point after a
converged run, the last iterate after a failed one.

```julia
kappa = reportCondition(J)          # prints and returns the estimate
kappa = condestJacobian(J)          # estimate only
kappa = condestJacobian(J; exact = true)  # dense SVD, small nets only
kappa = condestJacobian(net)        # solver Jacobian at the net's current state
```

Always on, without configuration:

| Use | Where |
|---|---|
| Result field | classic result log: one `Jacobian cond.` line with the verdict (one Hager estimate on the LU the solver factored; non-NR runs reconstruct the Jacobian); run metadata `jacobian_condition_estimate` / `_verdict` / `_line` |
| Web UI | run overview row `Jacobian condition` for rectangular NR runs (`n/a` for APSLF or DC); Diagnose runs (`run_diagnostics = true`, the Diagnose button) add the estimate and a recommendation for an ill-conditioned Jacobian to the `diagnose.log` Diagnosis section |
| Example | `examples/powerflow/exp_condition_number.jl`, including how conditioning collapses when a bus becomes nearly disconnected |
