# Sparlectra.jl

[![Documentation](https://github.com/Welthulk/Sparlectra.jl/actions/workflows/jekyll-gh-pages.yml/badge.svg)](https://welthulk.github.io/Sparlectra.jl/)
[![License: Apache-2.0](https://img.shields.io/badge/License-Apache--2.0-blue.svg)](LICENSE)
[![Julia](https://img.shields.io/badge/Julia-1.x-9558B2.svg)](https://julialang.org/)
[![GitHub release](https://img.shields.io/github/v/release/Welthulk/Sparlectra.jl)](https://github.com/Welthulk/Sparlectra.jl/releases)

**Sparlectra.jl is a Julia package for transparent power-system calculations.**
<a href="https://github.com/Welthulk/Sparlectra.jl/tree/main/"><img align="left" width="100" src="docs/src/assets/logo.png" style="margin-right: 20px" /></a>

It focuses on understandable and inspectable algorithms for AC power flow, PV/PQ switching, MATPOWER-based studies, transformer-control experiments, and experimental weighted-least-squares state estimation.

Sparlectra is intended for engineers, researchers, and students who want to understand, test, and extend power-flow and state-estimation algorithms without hiding the numerical core behind a black box.

---

## Why Sparlectra?

| Need | Sparlectra focus |
|---|---|
| Understand Newton-Raphson power-flow internals | Transparent rectangular complex-state formulation |
| Study PV/PQ switching | Explicit Q-limit and active-set diagnostics |
| Work with MATPOWER cases | Import helpers, local casefile workflow, and comparison diagnostics |
| Teach or inspect power-system algorithms | Readable Julia implementation with direct model access |
| Experiment with state estimation | Experimental nonlinear WLS workflow |
| Integrate external solvers | `PFModel` / `PFSolution` interface |
| Investigate voltage-control behavior | Transformer tap and voltage-control experiments |

Sparlectra is not primarily an optimization framework it focuses on solver mechanics, model transparency, diagnostics, and engineering-oriented experimentation.

---

## Main features

### AC power flow

- Rectangular complex-state Newton-Raphson power flow as the main solver path.
- Sparse-matrix-oriented implementation for realistic network studies.
- PV/PQ bus handling with Q-limit checks and active-set behavior.
- Slack, PV, and PQ bus typing derived from attached power-system components.
- Loss calculations and power-flow result reporting.

### MATPOWER workflow

- MATPOWER-compatible case import helpers.
- On-demand local casefile workflow for common benchmark cases.
- Diagnostics for comparing imported models against reference voltage profiles.
- Support for branch, transformer, shunt, generator, and load data relevant to AC power-flow studies.

### Network modeling

Sparlectra can represent typical transmission-grid components:

- buses
- transmission lines
- two-winding and three-winding transformers
- generators
- loads
- shunts
- topological impedanceless bus links via `addLink!`
- π-equivalent branch models for lines and transformers

Transmission lines and transformers can be represented using π-equivalent models. This supports direct use of many CIM- and manufacturer-style data sets without forcing an early conversion into overly simplified models.

### State estimation

Sparlectra includes an experimental nonlinear weighted-least-squares state-estimation workflow.

Current focus:

- WLS state estimation
- measurement-building helpers
- observability-related utilities
- residual inspection
- zero-injection pseudo-measurements for passive buses

State estimation should currently be treated as experimental. The implementation is useful for learning, prototyping, and diagnostic experiments, but some production-style workflows such as full bad-data processing are still part of the roadmap.

### Transformer and voltage-control experiments

Sparlectra provides functionality for transformer-control studies, including complex tap behavior and voltage-control experiments. This is useful for investigating how transformer ratios, phase shifts, and voltage-control logic interact with AC power-flow calculations.

---

## Installation

Install Sparlectra from the Julia package manager:

```julia
using Pkg
Pkg.add("Sparlectra")
```

Then load it with:

```julia
using Sparlectra
```

---

## Local Web UI

PowerFlow aborts are cooperative. Active status pages report the current
execution phase and cancellation timestamps. Sparse factorization and other
numerical calls may not be interruptible; after 60 seconds of pending abort,
the page offers **Hard reset Web UI**. This marks the run result invalid
(`aborted_unknown`) and cleanly shuts down the local server. Restart it with
`julia --project=. start_webui.jl`. A future process-isolated solver worker
would permit hard termination while keeping the Web UI process alive.

PowerFlow runs can optionally export detailed Excel-friendly CSV artifacts.
Enable **Export detailed result CSV files** to create
`bus_voltages_complex.csv` (polar and rectangular complex bus voltages plus
bus injections/loads) and `branch_flows.csv` (directional active/reactive
branch flows and losses). This option is off by default. New diagnostic runs
write `diagnose.log`; the artifact viewer still supports `diagnose.txt` files
from older run directories.
For Excel installations that expect locale-style CSV files, enable
**Use Excel CSV format with semicolon delimiter** together with the export
option. It changes the field separator from comma to semicolon and is off by
default.

Installed users can start the local Web UI without knowing the package
installation directory:

```julia
using Sparlectra

server = Sparlectra.start_sparlectra_webui(open_browser = true)
wait(server.task)
```

On first start, the Web UI provisions a user-writable configuration file and a small case cache, and creates its output and log directories automatically. Installed users do not need to locate the package directory. Its default output root is:

- Windows: `%LOCALAPPDATA%\Sparlectra\WebUI\runs`
- Linux: `$XDG_STATE_HOME/sparlectra/webui/runs`, or
  `~/.local/state/sparlectra/webui/runs` when `XDG_STATE_HOME` is unset
- macOS: `~/Library/Application Support/Sparlectra/WebUI/runs`

The operation log is stored in the user Web UI `logs` directory, and downloaded/generated
cases are cached in the sibling user Web UI `data/mpower` directory. To choose another root:

```julia
server = Sparlectra.start_sparlectra_webui(
    output_root = "my_sparlectra_runs",
    open_browser = true,
)
wait(server.task)
```

Repository developers should run `julia --project=. start_webui.jl`;
`start_webui.jl` is the maintained developer launcher.

Stop the server with the **Stop Web UI** button, `close(server)`, or `Ctrl+C`.
Run abort is cooperative: the UI changes the run state to `aborting` immediately, while
an already-blocking solver phase may continue until that phase returns. The rectangular
large-case path checks cancellation around Y-bus construction and start projection,
at Newton/Q-limit iteration boundaries, and immediately after each Newton step. After
terminal `aborted`, the active-run guard is released and another run can be submitted;
active runs cannot be deleted until they are terminal.
Queued, running, and aborting status pages refresh automatically every two seconds;
terminal pages stop refreshing, and the manual **Refresh status** link remains available.
The shared layout shows the running Sparlectra version. Support events in
`webui_operations.jsonl` include that version and an unambiguous UTC timestamp with a
`Z` suffix. Automatic status refreshes are not logged as user actions. The log is
compacted after it exceeds 10,000 valid JSONL entries, preserving the newest 1,000
entries; the existing 10 MiB single-backup size guard remains as an additional limit.

---

## Quick start

`run_sparlectra` is the preferred public framework entry point. It owns MATPOWER
import, configuration, optional control-loop execution, solving, post-processing,
and configured output. For AC power-flow examples, `run_acpflow` is kept as a thin
alias with the same minimal configuration-driven signature. Both names return
`SparlectraRunResult`. Solver and import behavior is controlled through
`SparlectraConfig` or YAML rather than a long list of runner keywords.

This example is runnable from a fresh checkout or package installation. The
`ensure_casefile` helper downloads `case14.m` on demand if it is not already
available locally.

```julia
using Sparlectra

case_path = ensure_casefile("case14.m")

result = run_sparlectra(
    casefile = basename(case_path),
    path = dirname(case_path),
)

println(result.outcome)
println(result.iterations)
println(result.final_mismatch)
```

For deterministic configuration-driven MATPOWER batches, list ordered case
names under `matpower_import.cases` and call `run_sparlectra_cases(config = cfg)`.
The batch helper returns one `SparlectraRunResult` per case in configured order,
uses `matpower_import.case` only when the list is empty, and resolves bare
standard case names through `ensure_casefile` on demand. `run_sparlectra` itself
remains a single-run API.

For a manually constructed network, pass the network directly and read the
solved model from `result.net`:

```julia
using Sparlectra

net = Net(name = "demo", baseMVA = 100.0)
# build network ...

cfg = active_sparlectra_config()
result = run_sparlectra(net = net, config = cfg)

vm = getNodeVm(result.net.nodeVec[1])
```

| Layer | Function | Purpose |
|---|---|---|
| Framework | `run_sparlectra` (`run_acpflow` alias) | Import/config/control/solve/output orchestration for one run |
| Framework batch | `run_sparlectra_cases` | Sequential deterministic execution of configured `matpower_import.cases` |
| Solver | `runpf!` | Solve an already built `Net` using `PowerFlowConfig` |
| Control | `run_control!` | Execute outer-loop controllers |
| Import | `createNetFromMatPowerFile` | Convert a MATPOWER file into a `Net` without running the framework workflow |

---

## Minimal modeling example

A typical Sparlectra workflow is:

1. create or import a network,
2. assign components and operating data,
3. run the solver,
4. inspect voltages, flows, losses, and diagnostics.

For most users, importing a MATPOWER-compatible case is the fastest way to start. For custom models, see the network and branch-model documentation.

---

## Documentation

The full documentation is available here:

<https://welthulk.github.io/Sparlectra.jl/>

Useful entry points:

- [Changelog](docs/src/changelog.md) — version history and release notes
- [Networks](docs/src/networks.md) — creating and manipulating network models
- [Branch Model](docs/src/branchmodel.md) — line and transformer branch modeling
- [Import/Export](docs/src/import.md) — importing and exporting network data
- [Workshop](docs/src/workshop.md) — guided examples and exercises
- [State Estimation](docs/src/state_estimation.md) — WLS state-estimation workflow
- [Feature Matrix](docs/src/feature_matrix.md) — power-flow and state-estimation capability overview
- [Network Reports](docs/src/netreports.md) — machine-readable `ACPFlowReport` output
- [Solver Guide](docs/src/solver.md) — numerical solver formulations
- [Transformer Control](docs/src/transformer_control.md) — transformer tap and voltage-control behavior
- [Function Reference](docs/src/reference.md) — public API reference

---

## Bus typing in power flow

`addBus!` defines the electrical node data, such as bus name, nominal voltage, optional voltage magnitude, and optional voltage angle.

The operational power-flow bus type is resolved from attached components:

- **Slack bus**: at least one slack-defining prosumer, for example an external-network injection with reference priority.
- **PV bus**: at least one regulating generator or controller with a voltage setpoint.
- **PQ bus**: default fallback for loads, non-regulating generation, or mixed non-regulating injections.

The legacy `busType` argument in `addBus!` may still be accepted for compatibility, but operational power-flow typing is derived from the network model.

---

## Test cases and benchmark data

Sparlectra does not ship third-party power-system test cases by default.

Instead, MATPOWER-compatible case files such as `case14.m` or `case118.m` can be downloaded on demand using helper functions and scripts. Downloaded files are stored locally and are not part of the Sparlectra source distribution.

This keeps the repository lightweight while still allowing reproducible experiments and benchmark-style workflows.

---

## Who is this for?

Sparlectra may be useful if you are:

- learning AC power-flow and state-estimation methods,
- teaching Newton-Raphson power-flow mechanics,
- testing PV/PQ switching and Q-limit behavior,
- investigating MATPOWER import behavior,
- experimenting with transformer-control concepts,
- building your own solver and need a transparent reference workflow,
- interested in power-system algorithms implemented in Julia.

---
## Contributing

Contributions, bug reports, test cases, and documentation improvements are welcome.

Before contributing, please read:

- [CONTRIBUTING.md](CONTRIBUTING.md)
- [CODE_OF_CONDUCT.md](CODE_OF_CONDUCT.md)

Good first contributions include:

- improving examples,
- adding small reproducible test networks,
- extending documentation,
- reporting MATPOWER import edge cases,
- improving diagnostics and error messages.

---

## Citing Sparlectra

If you use Sparlectra.jl in teaching, research, presentations, or project reports, please cite the repository.

```bibtex
@software{sparlectra_jl,
  title  = {Sparlectra.jl: Transparent Power-Flow and State-Estimation Experiments in Julia},
  author = {Schmitz, Udo},
  year   = {2026},
  url    = {https://github.com/Welthulk/Sparlectra.jl}
}
```

---

## License

Sparlectra.jl is licensed under the Apache License, Version 2.0.

See [LICENSE](LICENSE) for the full license text.

---
