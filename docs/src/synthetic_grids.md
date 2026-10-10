# Synthetic Tiled Grids

The tiled-grid builder creates reproducible AC power-flow benchmark
networks of scalable size from Julia code, without a case file.

## Builder API

```julia
using Sparlectra

net, meta = build_synthetic_tiled_grid_net(1000; aspect_ratio = 1.0)
result = run_sparlectra(net = net)
println(result.outcome)
```

`build_tiled_grid_net` is an alias. The requested bus count is an upper
bound: the builder chooses the largest grid with `rows * cols <= max_buses`
whose `cols / rows` is closest to `aspect_ratio`.

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

All branches are PI-model lines in the Sparlectra convention (series
impedance `r + im*x` and total shunt admittance `g + im*b` in p.u., shunt
split half/half), see [Branch model](branchmodel.md).

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

The metadata holds requested and actual bus counts, `rows`, `cols`,
`branch_count`, the bus roles and the scheduled generation and load
(MW/MVAr).

## YAML configuration utility

`load_yaml_dict` is Sparlectra's YAML subset parser:

```julia
cfg = load_yaml_dict(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH)
```

It supports `#` comments, nested 2-space-indented dictionaries, scalar
key-value pairs (booleans, `null`/`~`, integers, floats, symbols such as
`:rectangular`, strings) and one-line scalar lists such as
`[100, 300, 500]`.

## Running the example

```bash
julia --project=. examples/powerflow/exp_synthetic_tiled_grid_pf_perf.jl
julia --project=. examples/powerflow/exp_synthetic_tiled_grid_pf_perf.jl 100 300 1000
```

Every argument is a bus limit; without one the runner builds one grid of
at most 20 buses. The solver configuration comes from the file named in
`SPARLECTRA_CONFIGURATION_YAML`, else from `examples/configuration.yaml`
if it exists, else from the built-in defaults. The runner prints one row
per limit (limit, bus count, convergence flags, outcome, reason,
iterations, solve time in milliseconds) and writes the same table to a
timestamped log under `examples/_out`.

## Limitations

An artificial one-voltage-level grid for solver scaling and regression
checks: no protection, transformer or operational constraints, and
convergence behavior may differ from real cases.
