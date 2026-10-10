# External Solver Interface

An external power-flow solver (a research prototype, a custom Newton
variant) consumes a canonical network model from Sparlectra and returns its
result in one fixed form, comparable with the internal Newton-Raphson
solver. Sparlectra owns the network model; the external solver works on
exported data.

## Canonical Data Structures

### `PFModel`

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

Only active buses are included (isolated buses are removed), in the
internal `BusData` ordering of `createYBUS`. `V0` is the initial state
(flat start or case based). The reactive limits are optional; an external
solver may ignore them.

```julia
model = buildPfModel(net; flatstart=net.flatstart)
```

`buildPfModel` and `runpf_external!` accept the `start_projection*`
keywords of `runpf!`, so an external solver receives the projected start in
`model.V0` ([Start Projection for Difficult Seeds](@ref)).

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

builds a `PFModel` from `net`, calls the solver, evaluates a canonical
mismatch norm and writes the solution back into `net`.

## Output Control

The interface is silent by default. `showPfModel(model; verbose=false)` and
`showPfSolution(solution)` print summaries; neither overloads `Base.show`.

## APSLF (AnalyticLoadFlow.jl)

`ApslfSolver` (`src/acpflow/apslf_solver.jl`) is the built-in
`AbstractExternalSolver` for
[AnalyticLoadFlow.jl](https://github.com/SOPTIM/AnalyticLoadFlow.jl), an
analytic power-series (holomorphic-embedding style) solver and a required
dependency, so the solver, the `:apslf` stages of the contingency rescue
ladder and the automatic power-flow mode are always available:

```julia
using Sparlectra                 # AnalyticLoadFlow comes with it

solver = apslf_solver(order = 24, nr_polish = false)
```

- `order::Int`: highest power-series coefficient.
- `nr_polish::Bool`: Newton-Raphson polishing of the series result (off by
  default).
- `convergence_radius::Bool`: evaluate the Padé margin
  (`stability_from_Vcoeff`, the distance of the nearest Padé pole to
  `s = 1`) and report it next to the Jacobian condition (on by default;
  costs about one solve). The series radius, a root test on the growth of
  the voltage coefficients, is a different quantity and is reported in
  every APSLF run (below 1: the plain series does not reach `s = 1`).
- `mode::Symbol`: `:direct` (native PV handling) or `:outer` (PQ-only
  series plus an outer secant loop for PV enforcement).

The series is always evaluated with Padé `[L/M]` approximants
(AnalyticLoadFlow 0.10.0); there is no switch. The residual of an APSLF
solution is judged against the bus types the solver ended with: a machine
clamped at a reactive limit is a PQ bus at that limit, logged like a
rectangular Q-limit event (result table `PQ*`, Q-V check), and its voltage
setpoint does not count as a mismatch.

### Standalone use via the external-solver bridge

```julia
iters, status, sol = runpf_external!(net, apslf_solver(); tol = 1e-8)
```

Spec mapping `PFModel → AnalyticLoadFlow`: `Ybus → Y`, `busType → bustype`,
`real/imag(Sspec) → Pspec/Qspec`, `Vset → Vm`, `qmin_pu/qmax_pu → Qmin/Qmax`
(unconstrained when `model` carries no Q-limits), `slack_idx → slack`.
`sol.meta` carries the series/Padé `order`, the Padé margin (`stability`:
`dmin`/`pole`/`bus`/`level`), the series radius (`series_radius`:
`radius`/`bus`) and NR-polish bookkeeping; both bus fields are node indices
of the solved net, the framework run reports case bus numbers. Without a
finite reactive limit the direct kernel runs one pass: nothing can switch,
so further passes would repeat the same series solve.

### Framework integration

`power_flow.solver = apslf` routes `run_sparlectra` through `ApslfSolver`
instead of the rectangular solver, per island.
`power_flow.apslf_start.enabled = true` instead uses APSLF as a guarded
start-value generator ahead of the rectangular solve (`power_flow.solver`
stays `rectangular`). Keys:
[Solver selection (rectangular vs. APSLF)](@ref pf-solver-selection); form
controls: [Web UI guide](webui.md).

### Capability limits

APSLF is a different solution method, not a drop-in replacement of the
rectangular path:

- **No selectable start voltage.** The series starts from the analytic germ
  `V(s=0) = 1∠0`; `model.V0`, `start_mode` and `start_projection` have no
  effect.
- **No OLTC, tap-changer or phase-shifting-transformer control and no
  Q(U)/P(U) control.** With `power_flow.solver = apslf` the outer-loop tap
  and phase-shifter controllers are not run, Q(U)/P(U) controllers and
  voltage-dependent shunts are not applied, and their static setpoints are
  solved; one warning in the run log names how many of each were left out.
- **Q-limits are simple PV→PQ only.** AnalyticLoadFlow switches PV↔PQ
  against `Qmin`/`Qmax` inside the series solve; `power_flow.qlimits.guard`,
  `enforcement_mode`, hysteresis and cooldown do not apply.
- **No wrong-branch detection or rescue.** `power_flow.wrong_branch_detection`
  and the guarded current-iteration start pre-solve are rectangular only.

## Example: Exporting a Reference Solution

`examples/others/export_solution.jl` runs the internal solver, exports
`PFModel` and `PFSolution`, optionally runs an external solver and compares
voltage magnitudes and angles: a reference and regression test for external
solvers.
