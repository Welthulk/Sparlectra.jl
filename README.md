# Sparlectra.jl

Julia package for AC power flow, DC power flow, state estimation and IEC 60909 short circuit

[![Documentation](https://github.com/Welthulk/Sparlectra.jl/actions/workflows/jekyll-gh-pages.yml/badge.svg)](https://welthulk.github.io/Sparlectra.jl/)
[![Version](https://img.shields.io/badge/dynamic/toml?url=https%3A%2F%2Fraw.githubusercontent.com%2FWelthulk%2FSparlectra.jl%2Fmain%2FProject.toml&query=%24.version&label=version&prefix=v&color=blue)](https://github.com/Welthulk/Sparlectra.jl/blob/main/Project.toml)
[![License: Apache-2.0](https://img.shields.io/badge/License-Apache--2.0-blue.svg)](LICENSE)
[![Julia](https://img.shields.io/badge/Julia-1.x-9558B2.svg)](https://julialang.org/)
[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_tour.ipynb)
[![SBOM: SPDX](https://img.shields.io/badge/SBOM-SPDX-blue.svg)](https://github.com/Welthulk/Sparlectra.jl/releases/latest/download/Sparlectra.spdx.json)

**Sparlectra.jl is a Julia framework for steady-state grid analysis: power flow, state estimation, IEC 60909 short circuit, and N-1 contingency studies, with CGMES and MATPOWER exchange.**

<a href="https://github.com/Welthulk/Sparlectra.jl/tree/main/"><img align="left" width="100" src="docs/src/assets/logo.png" style="margin-right: 20px" /></a>

Sparlectra covers the complete workflow from network import through solving to configurable reporting. Grid data can be read from ENTSO-E CGMES (2.4.15 and 3.0), MATPOWER and native DTF sources, or built programmatically. Three solver backends are available: the built-in rectangular Newton-Raphson solver, a linear DC power flow, and an optional analytic power-series solver (APSLF, via the AnalyticLoadFlow.jl package extension) usable standalone, as the primary solver, or as a guarded start-value generator ahead of Newton-Raphson.

Every stage of the numerical pipeline is documented and accessible at runtime: model construction, Jacobian assembly, PV/PQ active-set handling and convergence behaviour can be inspected and instrumented. Together with deterministic, configuration-driven runs, explicit Q-limit and AC-island handling and machine-readable reporting, this suits production grid studies and planning work as well as algorithm development and solver benchmarking.

---

## Try it in your browser

No installation required, the workshop notebooks run on Google Colab:

| Notebook | Open |
|---|---|
| **Workshop tour, basic**: Newcomer to Advanced: first network step by step, model editing and Q-limits, slack types and short circuit, OLTC tap control, Q(U) control | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_tour.ipynb) |
| **Workshop tour, advanced**: Expert and Beyond: remote voltage control, a steerable HVDC link, state estimation, FACTS devices and their limits, N-1 contingency analysis, parallel sweeps on Julia threads | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_tour_advanced.ipynb) |
| **Workshop tour, part 3** (Advanced and up): coordinated transformer control: the circulating-current trap of parallel units, the master/slave group, and the CGMES RegulatingControl mapping | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_tour_control.ipynb) |
| **Slack types and short circuit** (Advanced): ideal slack vs. external-grid source vs. distributed slack on one 8-bus network, plus IEC 60909-0 fault currents | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_slack_short_circuit.ipynb) |
| **Distributed slack** (Advanced): how the participation weights are determined (schedule, capacity, headroom, imported APF/normalPF, explicit), how they are normalized, and the fallback when no valid participant exists | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_distributed_slack.ipynb) |
| **Transformer taps** (Advanced): ratio taps (OLTC) move voltage, phase taps (PST/Schrägregler) move power; device tap math, the closed control loop with the X(α) characteristic, and the 3WT star equivalent | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_transformers.ipynb) |
| **Series compensation (TCSC)** (Expert): why flow follows reactance, and how the series-reactance controller steers a loop-network flow split onto a branch target, including the at_limit case | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_series_compensation.ipynb) |
| **State estimation** (Advanced to Expert): derive a noisy measurement set from a reference power flow, check observability with a placement deep dive (incl. PMU phasors and the reference-angle offset), run the WLS estimator, and see when observability breaks down | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_state_estimation.ipynb) |

The notebooks are generated from the Literate.jl sources in [docs/lit/](docs/lit/); the same content is on the documentation pages under [Notebooks](https://welthulk.github.io/Sparlectra.jl/generated/workshop_tour/).

---

## Why Sparlectra?

| Requirement | Sparlectra approach |
|---|---|
| Reproducible AC power-flow studies | Deterministic, configuration-driven framework runs |
| Insight into Newton-Raphson internals | Rectangular complex-state formulation, open at every stage |
| Robust PV/PQ handling | Explicit Q-limit enforcement with active-set diagnostics |
| Grid data exchange | ENTSO-E CGMES (2.4.15 and 3.0), MATPOWER and native DTF import with validation against the delivered solution; CGMES export with roundtrip-stable mRIDs |
| Short-circuit analysis | Balanced short-circuit currents (Ik'', Sk'', i_p) per IEC 60909-0 from CGMES short-circuit data, with safety flagging of substituted defaults |
| Custom solver integration | Clean `PFModel` / `PFSolution` interface for external solvers |
| Voltage- and tap-control studies | Outer-loop control framework: transformer regulation (OLTC/PST/combined incl. tap-dependent PST reactance), remote voltage control via machine reactive power, and SVC-style variable-shunt voltage control |
| Realistic slack modeling | Distributed active-power slack with configurable participation factors (incl. imported MATPOWER `APF` / CGMES `normalPF`) |
| Alternative solver backend | Optional analytic power-series solver (APSLF, via AnalyticLoadFlow.jl): standalone, as the primary solver, or as an NR start-value generator |
| State estimation | Nonlinear weighted-least-squares workflow |
| Scalability | Sparse-matrix-oriented implementation for realistic network sizes |

---

## Main features

- Rectangular complex-state Newton-Raphson AC power flow, plus a linear DC power flow.
- Sparse-matrix-oriented implementation for realistic network studies.
- PV/PQ bus handling with Q-limit enforcement (machine capability curves where the data provides them) and active-set diagnostics.
- Distributed active-power slack over configurable participation factors: the primary-control picture instead of a single slack machine.
- Grid import from ENTSO-E CGMES 2.4.15 and 3.0 (EQ/SSH/TP/SV, boundary sets, multi-area assemblies, tap controllers, validation against the delivered SV profile, measured across the full ENTSO-E conformity collection including the 6209-bus RealGrid; per-case results incl. the documented non-converging completeness sets in [docs/dev/cgmes_testset_overview.md](docs/dev/cgmes_testset_overview.md)), MATPOWER cases and native DTF files.
- CGMES export as a complete delivery (EQ + TP + SSH + SV, optionally one re-importable ZIP: buses, lines, 2W/3W transformers incl. tap machinery, loads, machines, injections, SVCs, shunts, links, operating point and voltage state) with roundtrip-stable object identity: an exported and re-imported network solves to the same power flow and reproduces the original short-circuit evaluation; imported mRIDs are preserved, everything else gets deterministic ids.
- Balanced short-circuit analysis per IEC 60909-0 (`runShortCircuit!`): initial symmetrical current Ik'' (max/min case), Sk'' and peak current i_p per fault bus from harvested CGMES short-circuit data, with explicit safety flagging where defaults were substituted; all-bus sweeps run on Julia threads or via the Takahashi selected inverse (`sweep_method = :takahashi`).
- Network modeling: buses, lines, transformers, generators, loads, shunts, links, and π-equivalent branch models.
- Outer-loop control framework: transformer tap/voltage control (OLTC, PST with tap-dependent reactance, combined regulation), remote voltage control via machine reactive power, and SVC-style variable-shunt voltage control, all reported through one generic controllable-element view.
- Configuration-driven batch execution for systematic case studies.
- External-solver integration via the `PFModel` / `PFSolution` interface, including an optional analytic power-series solver (APSLF, via AnalyticLoadFlow.jl) usable standalone, as the primary solver, or as a Newton-Raphson start-value generator.
- Nonlinear weighted-least-squares state estimation.
- Machine-readable reporting (`ACPFlowReport`) and an optional local Web UI.

## Scope

- In scope: steady-state grid analysis: AC and DC power flow, WLS state estimation, IEC 60909 short-circuit currents, N-1 contingency batches, CGMES/MATPOWER/DTF exchange.
- Out of scope: dynamic simulation (RMS/EMT), optimal power flow, protection coordination beyond IEC 60909 currents.

---

## Installation

As a Julia package:

```julia
using Pkg
Pkg.add("Sparlectra")
```

```julia
using Sparlectra
```

### Web UI: one-line install (no git, no GitHub checkout)

The install scripts also run straight from this page. They install Julia
when it is missing (official juliaup installer; on Windows via winget),
download Sparlectra at its **latest tagged release** into a `Sparlectra/`
folder inside the current directory, ask once whether to build the optional
[fast-start sysimage](https://welthulk.github.io/Sparlectra.jl/fast_start/)
(a one-time build of 6 to 20 minutes; later Web UI starts then skip the
JIT warm-up, answering `n` or just pressing Enter skips it), and start the
Web UI. Run the command in the directory where Sparlectra should live.

Linux/macOS:

```sh
curl -fsSL https://raw.githubusercontent.com/Welthulk/Sparlectra.jl/main/install_webui.sh | sh
```

Windows (PowerShell):

```powershell
iwr -useb https://raw.githubusercontent.com/Welthulk/Sparlectra.jl/main/install_webui.bat -OutFile install_webui.bat; .\install_webui.bat
```

For unattended installs, `SPARLECTRA_BUILD_SYSIMAGE=1` (build) or `=0`
(skip) answers the sysimage question without a prompt. The sysimage can
also be built later at any time: `julia tools/build_sysimage.jl` from the
Sparlectra folder, or via the Web UI **Fast start** page.

The installer also offers a desktop shortcut (Windows `.lnk`) or launcher
(`.desktop` file, symlink fallback) for restarting the Web UI later;
`SPARLECTRA_CREATE_SHORTCUT=1` or `=0` answers that question without a
prompt.

**Updating:** re-running the same command is the update path. With an
existing `Sparlectra/` copy the installer compares it against the latest
tagged release and offers the update (default yes; `SPARLECTRA_UPDATE=1`
or `=0` answers without a prompt). The previous copy is kept as
`Sparlectra.old` until you delete it, and a fast-start sysimage from the
old version is detected as stale and can be rebuilt in the same run.

### Software Bill of Materials (SBOM)

Every GitHub release ships an SPDX SBOM (`Sparlectra-<tag>.spdx.json`) as a
release asset, generated by CI; the latest one is always available under
the stable link
[releases/latest/download/Sparlectra.spdx.json](https://github.com/Welthulk/Sparlectra.jl/releases/latest/download/Sparlectra.spdx.json)
(also the SBOM badge at the top of this page). It lists the full
dependency tree with licenses; since no `Manifest.toml` is checked in, it
describes the resolution at generation time (stated in its
`CreatorComment`). Generate
one locally with `julia tools/generate_sbom.jl [outfile]` from a checkout
(the PkgToSoftwareBOM tooling lives in a shared build environment and
never enters the package `Project.toml`);
`SPARLECTRA_SBOM_NO_LICENSESCAN=1` skips the per-package license scan.

### Web UI start scripts

The repository root ships ready-made scripts for the local Web UI:

| | Linux/macOS | Windows |
|---|---|---|
| Start (Julia already installed) | `./start_webui.sh` | `start_webui.bat` |
| Install and start | `./install_webui.sh` | `install_webui.bat` |

`start_webui.sh` / `start_webui.bat` start the Web UI from the checkout
(resolving Julia dependencies once on first start). If Julia is missing they
point at the install script instead of failing cryptically.

`install_webui.sh` / `install_webui.bat` are combined install-and-start
scripts: they install Julia when it is missing (official juliaup installer;
on Windows via winget, no git required there), obtain Sparlectra at its
**latest tagged release** (an existing checkout is used in place; a clean
git tree is moved to the release tag, local changes are never touched;
outside a checkout the release is cloned or downloaded next to the script),
offer the optional fast-start sysimage build (see above), and then start
the Web UI.

---

## Quick start

`run_sparlectra` is the primary framework entry point. It orchestrates import, configuration, optional control-loop execution, solving, post-processing and configured output. For AC power-flow scripts, `run_acpflow` remains available as a thin compatibility alias with the same signature.

The example below runs from a fresh checkout or package installation; `ensure_casefile` downloads `case14.m` on demand if it is not present locally.

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

Reading an ENTSO-E CGMES delivery works the same way: diagnose first, then import and solve:

```julia
using Sparlectra

summary = summarizeCGMES(path = ["grid.zip", "boundary.zip"])   # profiles, classes, dangling references
result  = importCGMES(path = ["grid.zip", "boundary.zip"])
runpf!(result.net, 30, 1e-8, 0)

cmp = compareWithSV(result)     # validate against the delivery's own SV profile
@show cmp.max_dvm
```

For custom network construction, batch execution, solver internals, and the local Web UI, see the documentation linked below.

---

## Local Web UI

Sparlectra ships with an optional browser-based local Web UI for power-flow studies, including run history, artifacts and case management. Cases can be selected from the local cache, uploaded (MATPOWER, DTF, CGMES ZIPs) or fetched by name. See the [Web UI documentation](https://welthulk.github.io/Sparlectra.jl/webui/) for setup and configuration.

**Configuration**: case selection, solver settings, control options and output configuration on a single page:

[![Sparlectra Web UI: PowerFlow run configuration](docs/src/assets/webui_v0.8.15.png)](docs/src/assets/webui_v0.8.15.png)

**Power flow run & history**: result with convergence report (left) and the run history (right):

[![Sparlectra Web UI: PowerFlow result and run history](docs/src/assets/webui_powerflow_history.png)](docs/src/assets/webui_powerflow_history.png)

---

## API entry points

| Layer | Function | Purpose |
|---|---|---|
| Framework | `run_sparlectra` (`run_acpflow` alias) | Import/config/control/solve/output orchestration for one run |
| Framework batch | `run_sparlectra_cases` | Sequential deterministic execution of configured `matpower_import.cases` |
| Solver | `runpf!` | Solve an already built `Net` using `PowerFlowConfig` |
| Alternative solver | `apslf_solver` | Reachability point for the APSLF (AnalyticLoadFlow.jl) external-solver backend |
| Control | `run_control!` | Execute outer-loop controllers |
| Import | `createNetFromMatPowerFile` | Convert a MATPOWER file into a `Net` without the full framework workflow |
| Import | `importCGMES` / `createNetFromCGMES` | Read an ENTSO-E CGMES delivery into a `Net`, with `summarizeCGMES` for diagnosis and `compareWithSV` for validation |
| Export | `writeCGMESFiles` | Write a `Net` as a complete CGMES delivery (EQ+TP+SSH+SV, optional ZIP) with roundtrip-stable mRIDs |
| Short circuit | `runShortCircuit!` | Balanced short-circuit currents (IEC 60909-0) from harvested CGMES short-circuit data |

---

## Documentation

Full documentation: <https://welthulk.github.io/Sparlectra.jl/>

Key entry points:

- [Local Web UI](https://welthulk.github.io/Sparlectra.jl/webui/): browser-based local power-flow workflow
- [Networks](docs/src/networks.md): building and manipulating network models
- [Import/Export](docs/src/import.md) · [CGMES Import](docs/src/cgmes_import.md) · [CGMES Export](docs/src/cgmes_export.md): reading and writing grid data
- [Short-Circuit Analysis](docs/src/short_circuit.md): balanced short-circuit currents per IEC 60909-0
- [Branch Model](docs/src/branchmodel.md) · [Remote Voltage Control](docs/src/remote_voltage_control.md): line/transformer modeling, tap and voltage control
- [Solver Guide](docs/src/solver.md) · [External Solvers](docs/src/external_solvers.md): numerical formulations and the `PFModel`/`PFSolution` interface
- [State Estimation](docs/src/state_estimation.md): WLS state-estimation workflow
- [Feature Matrix](docs/src/feature_matrix.md): capability overview
- [Function Reference](docs/src/reference.md) · [Workshop tour](docs/src/generated/workshop_tour.md): API reference and guided examples
- [Changelog](docs/src/changelog.md): version history

---

## Contributing

Contributions, bug reports, test cases and documentation improvements are welcome, particularly reproducible test networks, import edge cases (CGMES, MATPOWER, DTF), and improved diagnostics.

Questions and feedback: see [Discussions](https://github.com/Welthulk/Sparlectra.jl/discussions).

Please read [CONTRIBUTING.md](CONTRIBUTING.md) and [CODE_OF_CONDUCT.md](CODE_OF_CONDUCT.md) before contributing.

---

## Citing Sparlectra

If you use Sparlectra.jl in research, engineering studies, presentations, or reports, please cite the repository:

```bibtex
@software{sparlectra_jl,
  title  = {Sparlectra.jl: A Power-Flow and State-Estimation Framework in Julia},
  author = {Schmitz, Udo},
  year   = {2026},
  url    = {https://github.com/Welthulk/Sparlectra.jl}
}
```

---

## License

Sparlectra.jl is licensed under the Apache License, Version 2.0.

See [LICENSE](LICENSE) for the full license text.
