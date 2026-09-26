# External Solver Interface

External power-flow solvers (research prototypes, custom Newton variants)
consume a canonical network model from Sparlectra and return results in one
fixed form. Sparlectra builds and owns the network model, the external
solver works on exported data only, and the results are comparable with
the internal Newton-Raphson solver.

## Canonical Data Structures

### `PFModel`

A compressed, solver-independent power-flow model:

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

* Only active buses are included (isolated buses are removed).
* Bus ordering matches the internal `BusData` ordering used by `createYBUS`.
* `V0` is the canonical initial state (flat start or case based).
* Reactive power limits are optional; external solvers may ignore them.

```julia
model = buildPfModel(net; flatstart=net.flatstart)
```

For difficult starts the builder applies the start projection of the
internal rectangular solver:

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

External solvers receive the projected start in `model.V0`;
`runpf_external!` forwards the same keywords to `buildPfModel`.

### `PFSolution`

```julia
PFSolution(
  V::Vector{ComplexF64},   # PF ordering
  converged::Bool,
  iters::Int,
  residual_inf::Float64,
  meta::Any,
)
```

`V` uses PF ordering, the ordering of `model.busIdx_net`.

## External Solver Contract

```julia
abstract type AbstractExternalSolver end

solvePf(solver::AbstractExternalSolver, model::PFModel; kwargs...) -> PFSolution
```

Sparlectra imposes no numerical method, only this data contract.

## Running an External Solver

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

Builds a `PFModel` from `net`, calls the solver, evaluates a canonical
mismatch norm and writes the solution back into `net`.

## Output Control

The interface is silent by default. For inspection:

```julia
showPfModel(model; verbose=false)
showPfSolution(solution)
```

Neither overloads `Base.show`, so there are no global display side effects.

## APSLF (AnalyticLoadFlow.jl)

`ApslfSolver` (`src/acpflow/apslf_solver.jl`) is the built-in
`AbstractExternalSolver` for
[AnalyticLoadFlow.jl](https://github.com/Welthulk/AnalyticLoadFlow.jl), an
analytic power-series (holomorphic-embedding style) load-flow solver.
AnalyticLoadFlow.jl is a required dependency, so the solver, the `:apslf`
stages of the contingency rescue ladder and the automatic power-flow mode
are always available:

```julia
using Sparlectra                 # AnalyticLoadFlow comes with it

solver = apslf_solver(order = 24, use_pade = true, nr_polish = false)
```

- `order::Int`: highest power-series coefficient.
- `use_pade::Bool`: evaluate the voltage series via Padé `[L/M]`
  approximants instead of direct Taylor summation (usually a larger
  convergence radius).
- `nr_polish::Bool`: Newton-Raphson polishing step on the series result
  (off by default).
- `convergence_radius::Bool`: evaluate the Padé-pole margin
  (`stability_from_Vcoeff`) and report it as the APSLF convergence radius
  next to the Jacobian condition (on by default; costs about one solve).
- `mode::Symbol`: `:direct` (native PV handling) or `:outer` (PQ-only
  series plus an outer secant loop for PV enforcement).

The residual of an APSLF solution is judged against the bus types the
solver ended with: a machine clamped at a reactive limit is a PQ bus at
that limit, logged like a rectangular Q-limit event (result table `PQ*`,
Q-V check), and its voltage setpoint does not count as a mismatch.

### Standalone use via the external-solver bridge

```julia
iters, status, sol = runpf_external!(net, apslf_solver(); tol = 1e-8)
```

Spec mapping `PFModel → AnalyticLoadFlow`: `Ybus → Y`, `busType → bustype`,
`real/imag(Sspec) → Pspec/Qspec`, `Vset → Vm`, `qmin_pu/qmax_pu → Qmin/Qmax`
(unconstrained when `model` carries no Q-limits), `slack_idx → slack`.
`sol.meta` carries the series/Padé `order`, an APSLF stability indicator
(`dmin`/`pole`/`bus`/`level`, from the distance of the Padé poles to the
evaluation point `s = 1`) and NR-polish bookkeeping.

### Framework integration

`power_flow.solver = apslf` routes `run_sparlectra` through `ApslfSolver`
instead of the rectangular Newton-Raphson solver, per island for
multi-island networks. `power_flow.apslf_start.enabled = true` instead uses
APSLF as a guarded start-value generator ahead of the rectangular solve
(`power_flow.solver` stays `rectangular`). Configuration reference:
[Solver selection (rectangular vs. APSLF)](powerflow_configuration.md#solver-selection-rectangular-vs-apslf);
form controls: [Web UI guide](webui.md).

### Capability limits

APSLF is a different solution method, not a drop-in replacement of the
rectangular path:

- **No selectable start voltage.** The series starts from the analytic germ
  `V(s=0) = 1∠0`; `model.V0`, `start_mode` and `start_projection` have no
  effect. This is inherent to the embedding.
- **No OLTC, tap-changer or phase-shifting-transformer control and no
  Q(U)/P(U) control.** `power_flow.solver = apslf` with any active
  controller is rejected up front with an error; there is no fallback to a
  partially controlled solve.
- **Q-limits are simple PV→PQ only.** AnalyticLoadFlow.jl switches PV↔PQ
  against `Qmin`/`Qmax` inside the series solve; `power_flow.qlimits.guard`,
  `enforcement_mode`, hysteresis and cooldown settings do not apply.
- **No wrong-branch detection or rescue.** `power_flow.wrong_branch_detection`
  and the guarded current-iteration start pre-solve are rectangular only.

## Example: Exporting a Reference Solution

`examples/others/export_solution.jl` runs the internal Newton-Raphson
solver, exports `PFModel` and `PFSolution`, optionally runs an external
solver and compares voltage magnitudes and angles: a reference and
regression test for external solvers.
