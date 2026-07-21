# External Solver Interface

Sparlectra provides a **solver-agnostic interface** that allows external power-flow solvers
(e.g. research prototypes, new implementations, or custom Newton variants)
to consume a canonical network model and return results in a consistent form.

The goal is **clean separation**:

* Sparlectra builds and owns the network model
* External solvers operate purely on exported data
* Results can be compared against Sparlectra’s internal Newton–Raphson solver

---

## Canonical Data Structures

### `PFModel`

`PFModel` represents a compressed, solver-independent power-flow model.

```julia
PFModel(
  Ybus::AbstractMatrix{ComplexF64},
  baseMVA::Float64,
  busIdx_net::Vector{Int},
  busType::Vector{Symbol},   # :Slack | :PV | :PQ
  slack_idx::Int,
  Vset::Vector{Float64},
  Sspec::Vector{ComplexF64},
  V0::Vector{ComplexF64},
  qmin_pu::Vector{Float64},
  qmax_pu::Vector{Float64},
)
```

**Important notes:**

* Only *active buses* are included (isolated buses are removed).
* Bus ordering matches the internal `BusData` ordering used by `createYBUS`.
* `V0` provides a canonical initial state (flat-start or case-based).
* Reactive power limits are optional and may be ignored by external solvers.

A model is constructed using:

```julia
model = buildPfModel(net; flatstart=net.flatstart)
```

For difficult starts, the canonical model builder can apply the same start
projection used by the internal rectangular solver:

```julia
model = buildPfModel(net;
    flatstart = true,
    start_projection = true,
    start_projection_try_dc_start = true,
    start_projection_try_blend_scan = true,
    start_projection_blend_lambdas = [0.25, 0.5, 0.75],
    start_projection_dc_angle_limit_deg = 60.0,
)
```

External solvers receive the projected start in `model.V0`. `runpf_external!`
forwards the same keyword arguments to `buildPfModel`.

---

### `PFSolution`

External solvers must return a `PFSolution`:

```julia
PFSolution(
  V::Vector{ComplexF64},   # PF ordering
  converged::Bool,
  iters::Int,
  residual_inf::Float64,
  meta::Any,
)
```

The voltage vector `V` **must use PF ordering**
(i.e. the same ordering as `model.busIdx_net`).

---

## External Solver Contract

External solvers implement the following interface:

```julia
abstract type AbstractExternalSolver end

solvePf(solver::AbstractExternalSolver, model::PFModel; kwargs...) -> PFSolution
```

Sparlectra itself does **not** impose any numerical method.
Only the data contract must be respected.

---

## Running an External Solver

For convenience, Sparlectra provides:

```julia
runpf_external!(
  net,
  solver;
  tol=1e-8,
  flatstart=net.flatstart,
  include_limits=false,
  show_model=false,
  show_solution=false,
)
```

This function:

1. Builds a `PFModel` from `net`
2. Calls the external solver
3. Evaluates a canonical mismatch norm
4. Writes the solution back into `net`

---

## Output Control

The external solver interface is **silent by default**.

Optional helpers are provided for inspection and debugging:

```julia
showPfModel(model; verbose=false)
showPfSolution(solution)
```

These functions intentionally do **not** overload `Base.show`
to avoid global display side effects.

---

## APSLF (AnalyticLoadFlow.jl)

`ApslfSolver` is a built-in `AbstractExternalSolver` implementation that bridges
to [AnalyticLoadFlow.jl](https://github.com/Welthulk/AnalyticLoadFlow.jl), an
analytic power-series (holomorphic-embedding-style) load-flow solver. It ships
as a Julia [package extension](https://pkgdocs.julialang.org/v1/creating-packages/#Conditional-loading-of-code-in-packages-(extensions))
in `ext/SparlectraAnalyticLoadFlowExt.jl`: AnalyticLoadFlow.jl is only a
**weak dependency** of Sparlectra, so it never gets installed for users who
don't need it. The adapter, its config keys, and the framework wiring only
become active once the host session loads both packages:

```julia
using Sparlectra
using AnalyticLoadFlow   # activates ext/SparlectraAnalyticLoadFlowExt.jl
```

### Constructing a solver

```julia
solver = apslf_solver(; order = 40, use_pade = true, nr_polish = true, mode = :direct)
```

`apslf_solver` is the public reachability point in the base package: if
AnalyticLoadFlow.jl is not loaded, it raises a clear error
(`"AnalyticLoadFlow.jl nicht installiert — ..."`) instead of silently doing
nothing or falling back to another solver.

- `order::Int` — highest power-series coefficient to compute.
- `use_pade::Bool` — evaluate the voltage series via Padé `[L/M]` approximants
  instead of direct Taylor summation (generally improves the convergence
  radius).
- `nr_polish::Bool` — run a Newton-Raphson polishing step on the series
  result to tighten the final mismatch.
- `mode::Symbol` — `:direct` (native PV handling) or `:outer` (PQ-only series
  plus an outer secant loop for PV enforcement).

### Standalone use via the external-solver bridge

```julia
iters, status, sol = runpf_external!(net, apslf_solver(); tol = 1e-8)
```

This follows the same `buildPfModel` → `solvePf` → write-back contract as any
other `AbstractExternalSolver` (see above). The `PFModel → AnalyticLoadFlow`
spec mapping is: `Ybus → Y`, `busType → bustype`, `real/imag(Sspec) →
Pspec/Qspec`, `Vset → Vm`, `qmin_pu/qmax_pu → Qmin/Qmax` (unconstrained when
`model` carries no Q-limits), `slack_idx → slack`. `sol.meta` carries
solver-specific diagnostics: the series/Padé `order`, an APSLF stability
indicator (`dmin`/`pole`/`bus`/`level`, from the distance of Padé poles to the
physical evaluation point `s = 1`), and NR-polish bookkeeping.

### Framework integration

`power_flow.solver = apslf` routes the central framework run
(`run_sparlectra`) through `ApslfSolver` instead of the internal rectangular
Newton-Raphson solver — including per-island handling for networks with
multiple AC islands. `power_flow.apslf_start.enabled = true` uses APSLF as a
**guarded start-value generator** ahead of the rectangular Newton-Raphson
solve instead (`power_flow.solver` stays `rectangular` in that mode). See
[Solver selection (rectangular vs. APSLF)](powerflow_configuration.md#solver-selection-rectangular-vs-apslf)
for the full configuration reference and the [Web UI guide](webui.md) for the
corresponding form controls.

### Capability limits

APSLF is a genuinely different solution method from the rectangular
Newton-Raphson solver, not a drop-in replacement with identical modeling
depth. Concretely, compared to the internal rectangular path:

- **No selectable start voltage.** The series always starts from the
  canonical analytic germ `V(s=0) = 1∠0`; `model.V0`, `start_mode`, and
  `start_projection` have no effect on an APSLF solve. This is inherent to
  the embedding construction, not a missing feature.
- **No OLTC / tap-changer / phase-shifting-transformer control and no
  Q(U)/P(U) voltage-dependent control.** Runs with `power_flow.solver =
  apslf` and any active outer-loop controller (tap, PST, Q(U), P(U)) are
  rejected up front with a clear error — there is no silent fallback to a
  partially-controlled solve.
- **Q-limits are simple PV→PQ only.** AnalyticLoadFlow.jl performs its own
  internal PV↔PQ switching against `Qmin`/`Qmax` during the series solve; it
  does not reproduce the rectangular solver's active-set guard, hysteresis,
  or classical outer-loop enforcement modes (`power_flow.qlimits.guard`,
  `enforcement_mode`, hysteresis/cooldown settings do not apply to APSLF
  runs).
- **No wrong-branch detection/rescue.** `power_flow.wrong_branch_detection`
  and the guarded current-iteration start pre-solve are rectangular-solver-
  specific and are not part of the APSLF path.

## Example: Exporting a Reference Solution

See:

```
examples/export_solution.jl
```

This example:

* Runs Sparlectra’s internal Newton–Raphson solver
* Exports `PFModel` and `PFSolution`
* Optionally executes an external solver path
* Compares voltage magnitude and angle deviations

It is intended as a **reference and regression test** for external solvers.

---



