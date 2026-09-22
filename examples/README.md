# Local examples configuration

Copy `src/configuration.yaml.example` to `examples/configuration.yaml` and adjust
local settings there. The local `examples/configuration.yaml` override file is
intentionally ignored by Git.

# Running the example suites

The example programs are grouped into topic subfolders (`powerflow/`,
`state_estimation/`, `others/`). The suite runners execute a whole group of
examples in fresh subprocesses and write a summary (CSV + Markdown) plus
per-example logs to `examples/_out/<suite>/`:

```bash
julia --project=. examples/run_powerflow_suite.jl
julia --project=. examples/run_state_estimation_suite.jl
julia --project=. examples/run_others_suite.jl
```

Heavy or optional examples are skipped by default (`--include-heavy` /
`--include-optional` enables them); see `--help` and `--list` for all options.
DTF validation examples are covered by `examples/run_val_dtf_suite.jl`.

Two standalone scripts under `examples/others/` are not part of a suite:
`apslf_vs_nr_timing.jl` (solve-time table, CSV and SVG of the APSLF solver
against the rectangular Newton solver on the shipped `sp_` cases and
synthetic tiled grids; the source of the table in the performance page of
the docs) and `apslf_pv_diagnostic.jl` (four checks that name the layer on
which an APSLF solve of a network with PV buses fails on a given machine).
