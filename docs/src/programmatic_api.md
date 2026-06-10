# Programmatic Power-Flow API

`run_sparlectra_api` is the stable, non-interactive backend contract intended
for local applications and future GUI integrations. It accepts a MATPOWER case,
a configuration template, an output directory, and controlled dotted-key
overrides. The template is never modified.

```julia
using Sparlectra

casefile = ensure_casefile("case5.m")
result = run_sparlectra_api(
    casefile = casefile,
    config_file = Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH,
    output_dir = "results/example_api_run",
    config_overrides = Dict(
        "power_flow.tol" => 1.0e-8,
        "power_flow.max_iter" => 80,
        "power_flow.autodamp" => true,
        "output.logfile_results" => "compact",
        "benchmark.enabled" => false,
    ),
)

println(result.status)
println(result.artifacts)
```

The API always uses the existing `run_sparlectra` framework for numerical work.
It does not duplicate solver logic and does not read from stdin or request
manual confirmation.

## Result contract

`SparlectraApiResult` separates framework status from transport metadata:

- `status`, `success`, `converged`, and `solution_available` describe run state.
- `iterations`, `final_mismatch`, `reason`, and `message` describe the outcome.
- `casefile`, `config_file`, and `output_dir` record effective paths.
- `logfile`, `result_file`, and `artifacts` expose generated files explicitly.
- `raw_result` retains the underlying `SparlectraRunResult` for Julia callers.

Input and execution failures return `status == :failed` with a stable reason
such as `"casefile_not_found"`, `"invalid_configuration"`,
`"invalid_config_override"`, or `"execution_error"`.

## GUI-editable overrides

Only keys in `GUI_EDITABLE_CONFIG_KEYS` are accepted. The initial allowlist is:

- `power_flow.method`
- `power_flow.tol`
- `power_flow.max_iter`
- `power_flow.autodamp`
- `power_flow.autodamp_min`
- `power_flow.qlimits.enabled`
- `power_flow.wrong_branch_detection`
- `power_flow.start_mode.angle_mode`
- `power_flow.start_mode.voltage_mode`
- `output.logfile_results`
- `benchmark.enabled`
- `benchmark.samples`
- `benchmark.seconds`

Unknown keys, known but non-editable keys, invalid types, unsupported enum
values, and invalid ranges are rejected before execution.

## Effective configuration and artifacts

Every run with a valid configuration writes
`output_dir/effective_config.yaml`. Successful and failed calls also write
`run.log` and `result.json`. Artifact discovery recursively classifies these and
any generated CSV or report files, so clients never need to guess filenames.
Each `SparlectraApiArtifact` includes an absolute path, MIME type, existence
flag, byte size, kind, and description.

## Serialization

Use these dependency-free helpers:

```julia
dict_value = to_dict(result)
named_value = to_namedtuple(result)
json_text = to_json(result)
yaml_text = to_yaml(result)
```

The transport forms omit `raw_result` by default because a solved `Net` is not a
stable JSON/YAML representation. Pass `include_raw_result=true` only for custom
Julia-side inspection.

See `examples/exp_programmatic_api.jl` for a runnable example.
