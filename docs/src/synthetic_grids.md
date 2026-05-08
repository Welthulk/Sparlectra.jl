# Synthetic Tiled Grids

Sparlectra includes a dependency-free synthetic tiled-grid builder for creating
simple, reproducible AC power-flow benchmark networks directly from Julia code.
The builder is useful when you want scalable network sizes without reading a
MATPOWER, CGMES, or other external case file.

## Builder API

```julia
using Sparlectra

net, meta = build_synthetic_tiled_grid_net(1000; aspect_ratio = 1.0)
ite, erg, elapsed_s = run_net_acpflow(net = net, show_results = false)
```

`build_tiled_grid_net` is an alias for `build_synthetic_tiled_grid_net`.
The requested bus limit is an upper bound: Sparlectra chooses the largest
rectangular grid with `rows * cols <= max_buses` while keeping `cols / rows`
close to `aspect_ratio`.

## Topology

The generated network is a one-voltage-level rectangular grid with deterministic
row-major bus numbering:

```julia
synthetic_tiled_grid_bus_index(row, col, cols) = (row - 1) * cols + col
```

Bus names use the readable form `B_001_001`, `B_001_002`, and so on. The branch
set contains:

* horizontal PI-model AC lines between neighboring columns,
* vertical PI-model AC lines between neighboring rows,
* one diagonal PI-model AC line per rectangular tile from the upper-left bus to
  the lower-right bus.

The expected branch count is

```math
N_{branch} = rows(cols - 1) + (rows - 1)cols + (rows - 1)(cols - 1).
```

All synthetic branches use the Sparlectra AC PI branch convention. The series
impedance is `r + im*x` in p.u. and the total shunt admittance is `g + im*b` in
p.u.; Sparlectra splits the branch shunt half/half in Y-bus and branch-flow
calculations.

## Electrical setup

The default bus role assignment follows the synthetic benchmark convention:

* upper-left bus: slack bus,
* lower-left bus: scheduled generator bus,
* upper-right and lower-right buses: scheduled PQ load buses,
* all non-slack buses start with `vm_flat`, while the slack bus uses `vm_slack`.

Default parameters are intentionally modest:

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

Metadata returned by the builder includes requested and actual bus counts,
`rows`, `cols`, `branch_count`, bus role lists, and scheduled generation/load
values. Scheduled power metadata is reported in MW/MVAr. Line parameters and
voltage magnitudes are reported in p.u.

## YAML configuration utility

The example benchmark uses Sparlectra's small YAML subset parser:

```julia
cfg = load_yaml_dict("src/examples/exp_synthetic_tiled_grid_pf_perf.yaml.example")
```

This parser is intentionally not a full YAML implementation. It supports only
simple example-configuration files: comments beginning with `#`, nested
2-space-indented dictionaries, scalar key-value pairs, booleans, `null`/`~`,
integers, floating-point numbers, symbols such as `:rectangular`, strings, and
one-line scalar lists such as `[100, 300, 500]`.

## Running the example

```bash
julia --project=. src/examples/exp_synthetic_tiled_grid_pf_perf.jl
julia --project=. src/examples/exp_synthetic_tiled_grid_pf_perf.jl 100 300 1000
julia --project=. src/examples/exp_synthetic_tiled_grid_pf_perf.jl src/examples/exp_synthetic_tiled_grid_pf_perf.yaml
# if the .yaml file is missing, the runner tries .yaml.example automatically
julia --project=. src/examples/exp_synthetic_tiled_grid_pf_perf.jl src/examples/exp_synthetic_tiled_grid_pf_perf.yaml --max-buses=5000
```

When no configuration path is supplied, or when neither the requested YAML file nor its `.yaml.example` fallback is available, the example prints a message before using its built-in defaults. The example prints a compact summary with grid size, branch count, convergence,
iterations, solve time, mismatch diagnostics, system-build timing, allocation
counts, and total case runtime. It also writes a timestamped log file under
`src/examples/_out` and prints a small ASCII plot of `nbus` versus solve time.

## Limitations

* The synthetic grid is artificial and intended for solver scaling,
  diagnostics, and regression checks.
* It is a one-voltage-level grid.
* It does not represent realistic protection, voltage-level, transformer, or
  operational constraints.
* Large grids are useful for stress testing, but practical convergence behavior
  may differ from real transmission or distribution cases.
