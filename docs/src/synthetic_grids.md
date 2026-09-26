# Synthetic Tiled Grids

The dependency-free tiled-grid builder creates reproducible AC power-flow
benchmark networks of scalable size from Julia code, without any case
file.

## Builder API

```julia
using Sparlectra

net, meta = build_synthetic_tiled_grid_net(1000; aspect_ratio = 1.0)
result = run_sparlectra(net = net)
println(result.outcome)
```

`build_tiled_grid_net` is an alias for `build_synthetic_tiled_grid_net`.
The requested bus count is an upper bound: the builder chooses the largest
grid with `rows * cols <= max_buses` while keeping `cols / rows` close to
`aspect_ratio`.

## Topology

A one-voltage-level rectangular grid with row-major bus numbering:

```julia
synthetic_tiled_grid_bus_index(row, col, cols) = (row - 1) * cols + col
```

Bus names are `B_001_001`, `B_001_002`, and so on. Branches:

* horizontal PI-model AC lines between neighboring columns,
* vertical PI-model AC lines between neighboring rows,
* one diagonal PI-model AC line per tile from the upper-left to the
  lower-right bus.

The branch count is

```math
N_{branch} = rows(cols - 1) + (rows - 1)cols + (rows - 1)(cols - 1).
```

All branches use the Sparlectra AC PI convention: series impedance
`r + im*x` in p.u., total shunt admittance `g + im*b` in p.u., split
half/half in Y-bus and branch-flow calculations.

## Electrical setup

* upper-left bus: slack bus,
* lower-left bus: scheduled generator bus,
* upper-right and lower-right buses: scheduled PQ load buses,
* all non-slack buses start with `vm_flat`, the slack bus with `vm_slack`.

Default parameters:

```julia
aspect_ratio = 1.0
base_mva = 100.0
r = 0.01
x = 0.05
g = 0.0
b = 0.0
load_mw_per_right_corner = 50.0
load_mvar_per_right_corner = 15.0
generation_balance = 0.995
vm_slack = 1.0
vm_flat = 1.0
```

The returned metadata holds requested and actual bus counts, `rows`,
`cols`, `branch_count`, bus role lists and scheduled generation/load
values (MW/MVAr; line parameters and voltages in p.u.).

## YAML configuration utility

The example benchmark uses Sparlectra's YAML subset parser:

```julia
cfg = load_yaml_dict("examples/powerflow/exp_synthetic_tiled_grid_pf_perf.yaml.example")
```

It supports comments beginning with `#`, nested 2-space-indented
dictionaries, scalar key-value pairs, booleans, `null`/`~`, integers,
floating-point numbers, symbols such as `:rectangular`, strings, and
one-line scalar lists such as `[100, 300, 500]`.

## Running the example

```bash
julia --project=. examples/powerflow/exp_synthetic_tiled_grid_pf_perf.jl
julia --project=. examples/powerflow/exp_synthetic_tiled_grid_pf_perf.jl 100 300 1000
julia --project=. examples/powerflow/exp_synthetic_tiled_grid_pf_perf.jl examples/powerflow/exp_synthetic_tiled_grid_pf_perf.yaml
# if the .yaml file is missing, the runner tries .yaml.example automatically
julia --project=. examples/powerflow/exp_synthetic_tiled_grid_pf_perf.jl examples/powerflow/exp_synthetic_tiled_grid_pf_perf.yaml --max-buses=5000
```

Without a configuration path, or when neither the YAML file nor its
`.yaml.example` fallback exists, the example says so and uses built-in
defaults. It prints a summary (grid size, branch count, convergence,
iterations, solve time, mismatch, build timing, allocations, total
runtime), writes a timestamped log under `examples/_out` and plots `nbus`
versus solve time in ASCII.

## Limitations

An artificial one-voltage-level grid for solver scaling, diagnostics and
regression checks: no realistic protection, transformer or operational
constraints, and convergence behavior may differ from real transmission
or distribution cases.
