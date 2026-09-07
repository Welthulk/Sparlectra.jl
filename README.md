# Sparlectra.jl

Power-system analysis in Julia: AC and DC power flow, WLS state estimation, IEC 60909 short circuit, scenarios and N-1 contingency, import and export in various formats.

[![Documentation](https://github.com/Welthulk/Sparlectra.jl/actions/workflows/jekyll-gh-pages.yml/badge.svg)](https://welthulk.github.io/Sparlectra.jl/)
[![Version](https://img.shields.io/badge/dynamic/toml?url=https%3A%2F%2Fraw.githubusercontent.com%2FWelthulk%2FSparlectra.jl%2Fmain%2FProject.toml&query=%24.version&label=version&prefix=v&color=blue)](https://github.com/Welthulk/Sparlectra.jl/blob/main/Project.toml)
[![License: Apache-2.0](https://img.shields.io/badge/License-Apache--2.0-blue.svg)](LICENSE)
[![Julia](https://img.shields.io/badge/Julia-1.x-9558B2.svg)](https://julialang.org/)
[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_tour.ipynb)
[![SBOM: SPDX](https://img.shields.io/badge/SBOM-SPDX-blue.svg)](https://github.com/Welthulk/Sparlectra.jl/releases/latest/download/Sparlectra.spdx.json)

<a href="https://github.com/Welthulk/Sparlectra.jl/tree/main/"><img align="left" width="100" src="docs/src/assets/logo.png" style="margin-right: 20px" /></a>


Nothing in Sparlectra is a black box. The rectangular Newton-Raphson solver exposes model construction, Jacobian assembly, PV/PQ switching and convergence at runtime; a DC power flow and the analytic APSLF backend (AnalyticLoadFlow.jl) are available alongside it. State estimation reports observability and bad data, short-circuit sweeps and SE diagnostics use the Takahashi selected inverse. Every run is deterministic and configuration-driven, results are machine-readable. This suits grid studies as well as solver development and teaching.

The full capability list is in the [feature matrix](docs/src/feature_matrix.md).

---

## Try it in your browser

No installation required, the workshop notebooks run on Google Colab:

| Notebook | Open |
|---|---|
| **Tour, basic**: first network, model editing, Q-limits, slack types, OLTC and Q(U) control | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_tour.ipynb) |
| **Tour, advanced**: remote voltage control, HVDC link, state estimation, FACTS incl. UPFC, N-1, threaded sweeps | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_tour_advanced.ipynb) |
| **Transformer control**: parallel units, master/slave groups, CGMES RegulatingControl | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_tour_control.ipynb) |
| **CGMES**: anatomy of a delivery, bus-branch and node-breaker import, SV validation, export round trip | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_tour_cgmes.ipynb) |
| **Slack types and short circuit**: ideal, external grid, distributed slack, IEC 60909-0 currents | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_slack_short_circuit.ipynb) |
| **Distributed slack**: participation weights, normalization, fallback | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_distributed_slack.ipynb) |
| **Transformer taps**: OLTC, PST, X(α) characteristic, 3WT star equivalent | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_transformers.ipynb) |
| **FACTS flow control**: TCSC to UPFC | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_series_compensation.ipynb) |
| **Scenarios and N-1 screening** | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_scenarios.ipynb) |
| **State estimation**: measurement sets, observability, WLS, bad data | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_state_estimation.ipynb) |

The notebooks are generated from [docs/lit/](docs/lit/) and are also part of the [documentation](https://welthulk.github.io/Sparlectra.jl/generated/workshop_tour/).

---

## Installation

As a Julia package:

```julia
using Pkg
Pkg.add("Sparlectra")
```

### Web UI: one-line install

Installs Julia if missing, downloads the latest tagged release into `Sparlectra/` in the current directory, and starts the Web UI.

Linux/macOS:

```sh
curl -fsSL https://raw.githubusercontent.com/Welthulk/Sparlectra.jl/main/tools/install_webui.sh | sh
```

Windows (PowerShell):

```powershell
iwr -useb https://raw.githubusercontent.com/Welthulk/Sparlectra.jl/main/tools/install_webui.bat -OutFile install_webui.bat; .\install_webui.bat
```

> [!IMPORTANT]
> This downloads and runs a script. If you prefer to read it first: [install_webui.sh](tools/install_webui.sh), [install_webui.bat](tools/install_webui.bat).

Re-running the command updates an existing copy (the old one is kept as `Sparlectra.old`). From a checkout, `./start_webui.sh` or `start_webui.bat` starts the Web UI directly. Environment variables for unattended installs and all other options: [Web UI documentation](https://welthulk.github.io/Sparlectra.jl/webui/).

### Startup time
> [!IMPORTANT]
> Julia compiles on first use, so the first run in a fresh process is slow. The installer offers to build a sysimage once (10 to 20 minutes); afterwards `using Sparlectra` returns at once. How to build and use it in your own scripts: [integration guide](https://welthulk.github.io/Sparlectra.jl/integration/).

### SBOM

Every release ships an SPDX SBOM as a release asset: [Sparlectra.spdx.json](https://github.com/Welthulk/Sparlectra.jl/releases/latest/download/Sparlectra.spdx.json). `julia tools/generate_sbom.jl` builds one locally.

---

## Quick start

`run_sparlectra` runs the whole chain: import, configuration, control loops, solve, report. `ensure_casefile` downloads `case14.m` if it is not present locally.

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

An ENTSO-E CGMES delivery: diagnose, import, solve, compare with the delivered SV profile:

```julia
using Sparlectra

summary = summarizeCGMES(path = ["grid.zip", "boundary.zip"])
result  = importCGMES(path = ["grid.zip", "boundary.zip"])
runpf!(result.net, 30, 1e-8, 0)

cmp = compareWithSV(result)
@show cmp.max_dvm
```

Embedding Sparlectra in your own tooling (scripts, sysimage, API): [integration guide](https://welthulk.github.io/Sparlectra.jl/integration/).

---

## Local Web UI

Case management, run history and results in the browser. Cases can be taken from the local cache, uploaded (MATPOWER, DTF, CGMES ZIP) or fetched by name.

<a href="docs/src/assets/webui_v0.10.0.png"><img src="docs/src/assets/webui_v0.10.0.png" alt="Sparlectra Web UI: case, settings, run and result" /></a>

---

## Documentation

Full documentation: <https://welthulk.github.io/Sparlectra.jl/>

- [Feature matrix](docs/src/feature_matrix.md)
- [Integration guide](https://welthulk.github.io/Sparlectra.jl/integration/)
- [Web UI](https://welthulk.github.io/Sparlectra.jl/webui/)
- [Import/Export](docs/src/import.md), [CGMES import](docs/src/cgmes_import.md), [CGMES export](docs/src/cgmes_export.md)
- [Solver guide](docs/src/solver.md), [External solvers](docs/src/external_solvers.md)
- [State estimation](docs/src/state_estimation.md)
- [Short circuit](docs/src/short_circuit.md)
- [Function reference](docs/src/reference.md)
- [Changelog](docs/src/changelog.md)

---

## Contributing

Bug reports, test networks and import edge cases (CGMES, MATPOWER, DTF) are welcome. Questions: [Discussions](https://github.com/Welthulk/Sparlectra.jl/discussions). Please read [CONTRIBUTING.md](CONTRIBUTING.md) and [CODE_OF_CONDUCT.md](CODE_OF_CONDUCT.md) first.

## Citing

```bibtex
@software{sparlectra_jl,
  title  = {Sparlectra.jl: A Power-Flow and State-Estimation Framework in Julia},
  author = {Schmitz, Udo},
  year   = {2026},
  url    = {https://github.com/Welthulk/Sparlectra.jl}
}
```

## License

Apache License 2.0, see [LICENSE](LICENSE).
