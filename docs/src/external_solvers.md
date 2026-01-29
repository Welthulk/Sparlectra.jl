## External Solver Interface

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
model = buildPfModel(net; opt_sparse=true, flatstart=net.flatstart)
```

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
  opt_sparse=true,
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

## Example: Exporting a Reference Solution

See:

```
src/examples/export_solution.jl
```

This example:

* Runs Sparlectra’s internal Newton–Raphson solver
* Exports `PFModel` and `PFSolution`
* Optionally executes an external solver path
* Compares voltage magnitude and angle deviations

It is intended as a **reference and regression test** for external solvers.

---




