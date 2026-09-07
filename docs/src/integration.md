# Integration Guide

One page for integrators who read nothing else. Everything here links to
the detailed pages instead of repeating them.

## Quickstart

```julia
using Sparlectra

result = run_sparlectra_api(
  casefile = "case118.m",
  output_dir = "runs/first",
  config_overrides = Dict("power_flow.mode" => "auto"),
)

println(result.status)          # :succeeded or :not_converged
println(result.metadata["auto_profile"])
println(result.output_dir)      # artifacts live here
```

That is the whole integration: import happens from the case file (MATPOWER,
CGMES ZIP, or DTF), the auto mode picks the solver strategy, and the run
directory carries every artifact. See
[Programmatic API](programmatic_api.md) for the full call surface.

## The configuration file

Every option Sparlectra has is a key in one YAML file, and a run states in
`effective_config.yaml` what actually took effect. Three ways to set an
option, highest precedence first: an explicit override in the call, the
case configuration file next to the case (`<stem>.config.yaml`), and the
machine's YAML configuration.

**Creating one.** The packaged template is the reference, so start from a
copy of it:

```julia
using Sparlectra
cp(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "my_configuration.yaml")
```

Edit it in any text editor; it carries every key with its default and a
comment. In the Web UI the same file is reachable through the Settings
page (Configuration Editor, plus Check and Refresh), which writes the
values back and validates before saving.

**Using one.** Pass it per run, or load it once for the session:

```julia
result = run_sparlectra_api(casefile = "case118.m", output_dir = "runs/a",
                            config_file = "my_configuration.yaml")

Sparlectra.load_sparlectra_config!("my_configuration.yaml"; reload = true)
runpf!(net)                    # every later call uses it
```

Single values need no file at all: `config_overrides = Dict("power_flow.tol"
=> 1e-6)` in the API call, or `--set power_flow.tol=1e-6` on the command
line. Unknown keys and invalid values are refused by name before anything
runs.

**A case can bring its own settings.** A `<stem>.config.yaml` next to a
case file makes the case self-contained: its keys win over the machine
configuration, so the same case computes the same numbers on every
installation. `write_case_config(case, Dict(...))` writes one, and only
case-scope keys are accepted there (import conventions, the solver
surface, state estimation, short circuit); logging, benchmarking,
parallelism and Web UI settings stay with the machine.

The full key reference, the merge precedence in detail, and the
compatibility policy are in [Configuration](configuration.md); the
power-flow options with their trade-offs are in
[PowerFlow Configuration](powerflow_configuration.md).

## Startup latency

The first run in a fresh Julia process includes compilation and is slow;
later runs in the same process are fast. In order of preference: keep the
process running, or build a sysimage once with `buildSysimage()`. See
[Sysimage](sysimage.md) for build times and details. A one-line hint
reminds you of this at the first run of a native process; silence it with
`output.startup_latency_hint: false`.

Sparlectra needs an installed Julia; the sysimage shortens the start of
that Julia, it does not replace it.

### Using the sysimage without the Web UI

The sysimage is a plain Julia sysimage, not a Web UI feature. `buildSysimage()`
prints the path it wrote (below the Web UI state directory by default,
`target_dir` chooses another one), and any Julia call takes it with
`--sysimage`:

```bash
julia --sysimage /path/to/sparlectra.so --project=/path/to/project my_script.jl
julia --sysimage /path/to/sparlectra.so --project=/path/to/project -e 'using Sparlectra; runpf!(importSCF("case.scf.json"))'
```

The session starts with Sparlectra already compiled: `using Sparlectra`
returns at once, and the first `runpf!`, `runse!` or `runShortCircuit!`
runs at full speed. That holds for batch scripts, an interactive REPL, a
notebook kernel, and a scheduled job alike; nothing about it requires the
Web UI. Rebuild the image after updating Sparlectra, and note that a
sysimage cannot be moved between operating systems.

A relocatable executable with an embedded Julia runtime is not part of the
package. The build script for it lives in the checkout at
`tools/build_app.jl` and can be run directly
(`julia --project=. tools/build_app.jl --flavor=full`).

## Auto mode

`power_flow.mode = auto` inspects the imported network (size, AC islands,
R/X profile, phase shifters, start profile, generator Q limits) and fills
the start-value, step-control, and Q-limit strategy. On non-convergence a
bounded escalation ladder retries with stronger strategies, up to an
APSLF-seeded Newton-Raphson run and an optional DC fallback
(`power_flow.dc.fallback`, off by default and API-only).

What it never does: change the tolerance, mutate the network, override an
option you set explicitly, or report a DC result as an AC solution. Every
decision, precedence conflict, escalation attempt, and hint lands in the
`auto_mode_decision.log` artifact of the run.

## Q limits, in five sentences

Generators have reactive limits; when a limit binds, the bus switches from
voltage control (PV) to fixed injection (PQ), which is physics, not an
error. Bad or missing limit data is the most common cause of
non-convergence and strange voltages. Auto mode handles the switching,
guards, and enforcement-mode selection automatically. It reports what it
did as one-line hints (`auto_hints` in the result metadata, the decision
log, and the Web UI result page). Read the hints; they name the exact
buses and generators to check in the source data.

## Result contract

| Field / artifact | Meaning |
|---|---|
| `result.run_id` | unique id, also the run directory name |
| `result.status` | `:succeeded` or `:not_converged` |
| `result.converged` | AC convergence; NEVER true for a DC fallback |
| `result.metadata["dc_fallback_solution"]` | true when the net carries a DC approximation after AC failure |
| `result.metadata["auto_hints"]` | advisory hint texts; `auto_hint_ids` carries stable keys for programmatic matching |
| `result.metadata["auto_profile"]`, `auto_final_stage`, `auto_final_solver`, `auto_escalation_stages` | what auto mode decided and how far it escalated |
| `result.json` | the machine-readable result in the run directory |
| `effective_config.yaml` | the exact configuration the run used |
| `auto_mode_decision.log` | features, profile, reasons, conflicts, attempts, hints |

## When auto mode fails

Read the hints and `auto_mode_decision.log` first; they usually name the
generators or branches to check in the source data. If the case still does
not converge, open an issue and attach `result.json` and the case file.
