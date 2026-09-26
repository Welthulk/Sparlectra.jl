# Integration Guide

Sparlectra in your own scripts, jobs and tools: one API call, one
configuration file, one result contract. Details live on the linked pages.

## Quickstart

```julia
using Sparlectra
using SparlectraApp   # service API: environment app/, or app/ on the load path

result = run_sparlectra_api(
  casefile = "case118.m",
  output_dir = "runs/first",
  config_overrides = Dict("power_flow.mode" => "auto"),
)

println(result.status)          # :succeeded or :not_converged
println(result.metadata["auto_profile"])
println(result.output_dir)      # artifacts live here
```

The case file (MATPOWER, CGMES ZIP or DTF) is imported, auto mode picks the
solver strategy, and the run directory carries every artifact. Full call
surface: [Programmatic API](programmatic_api.md).

## Loading both packages

`using Sparlectra, SparlectraApp` works in a Julia started in the checkout
with `--project=app`. Elsewhere:

* Plain REPL (prompt `(@v1.13) pkg>`, the default environment): neither
  package is found, and Julia offers to install `Sparlectra` from the
  registry, which is the library alone. Activate the application
  environment first:

  ```julia
  using Pkg
  Pkg.activate("path/to/Sparlectra/app")   # the checkout's app/ directory
  Pkg.instantiate()                        # first time only
  using Sparlectra, SparlectraApp
  ```

* REPL on the library project (`julia --project=.`, the editor default for
  the checkout): only `SparlectraApp` is missing. Put `app/` on the load
  path, as `start_webui.jl` does, or activate `app/` as above:

  ```julia
  push!(LOAD_PATH, joinpath(pwd(), "app"))   # from the checkout root
  using Sparlectra, SparlectraApp
  ```

## The configuration file

Every option is a key in one YAML file, and `effective_config.yaml` in the
run directory states what took effect. Precedence, highest first: an
explicit override in the call, the case configuration file next to the case
(`<stem>.config.yaml`), the machine's YAML configuration.

**Creating one.** Copy the packaged template; it carries every key with its
default and a comment:

```julia
using Sparlectra
using SparlectraApp   # service API: environment app/, or app/ on the load path
cp(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "my_configuration.yaml")
```

In the Web UI the same file is edited on the Settings page (Configuration
Editor with Check and Refresh), which validates before saving.

**Using one.** Per run, or once per session:

```julia
result = run_sparlectra_api(casefile = "case118.m", output_dir = "runs/a",
                            config_file = "my_configuration.yaml")

Sparlectra.load_sparlectra_config!("my_configuration.yaml"; reload = true)
runpf!(net)                    # every later call uses it
```

Single values need no file: `config_overrides = Dict("power_flow.tol" =>
1e-6)` in the API call, or `--set power_flow.tol=1e-6` on the command line.
Unknown keys and invalid values are refused by name before anything runs.

**A case can bring its own settings.** A `<stem>.config.yaml` next to a
case file wins over the machine configuration, so the case computes the
same numbers on every installation. `write_case_config(case, Dict(...))`
writes one. Only case-scope keys are accepted there (import conventions,
solver surface, state estimation, short circuit); logging, benchmarking,
parallelism and Web UI settings stay with the machine.

Key reference, merge precedence and compatibility policy:
[Configuration](configuration.md). Power-flow options and their trade-offs:
[PowerFlow Configuration](powerflow_configuration.md).

## Startup latency

The first run in a fresh Julia process compiles and is slow, later runs are
fast: keep the process running, or build a [sysimage](sysimage.md) once
(that page also covers the standalone executable, a checkout tool that is
not part of the package). The one-line hint at the first run of a native
process is silenced with `output.startup_latency_hint: false`.

## Auto mode

`power_flow.mode = auto` inspects the imported network (size, AC islands,
R/X profile, phase shifters, start profile, generator Q limits) and fills
the start-value, step-control and Q-limit strategy. On non-convergence a
bounded escalation ladder retries with stronger strategies, up to an
APSLF-seeded Newton-Raphson run and an optional DC fallback
(`power_flow.dc.fallback`, off by default, API-only).

It never changes the tolerance, mutates the network, overrides an option
you set explicitly, or reports a DC result as an AC solution. Every
decision, precedence conflict, escalation attempt and hint lands in the
`auto_mode_decision.log` artifact of the run.

## Q limits

When a generator's reactive limit binds, its bus switches from voltage
control (PV) to fixed injection (PQ): physics, not an error. Bad or missing
limit data is the most common cause of non-convergence and strange
voltages. Auto mode handles switching, guards and enforcement mode, and
reports what it did as one-line hints (`auto_hints` in the result metadata,
the decision log, the Web UI result page). The hints name the buses and
generators to check in the source data.

## Result contract

| Field / artifact | Meaning |
|---|---|
| `result.run_id` | unique id, also the run directory name |
| `result.status` | `:succeeded` or `:not_converged` |
| `result.converged` | AC convergence; never true for a DC fallback |
| `result.metadata["dc_fallback_solution"]` | true when the net carries a DC approximation after AC failure |
| `result.metadata["auto_hints"]` | advisory hint texts; `auto_hint_ids` carries stable keys for programmatic matching |
| `result.metadata["auto_profile"]`, `auto_final_stage`, `auto_final_solver`, `auto_escalation_stages` | what auto mode decided and how far it escalated |
| `result.json` | the machine-readable result in the run directory |
| `effective_config.yaml` | the exact configuration the run used |
| `auto_mode_decision.log` | features, profile, reasons, conflicts, attempts, hints |

## When auto mode fails

Read the hints and `auto_mode_decision.log` first. If the case still does
not converge, open an issue and attach `result.json` and the case file.
