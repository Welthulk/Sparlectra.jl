# Programmatic API

The stable programmatic surface of Sparlectra is a small set of entry
points; everything else on the [reference pages](reference.md) is internal
and may change between versions.

## Stable entry points

- [`run_sparlectra`](@ref): the configured power flow on a network or
  case file, the core of every service path.
- [`importCGMES`](@ref): import a CGMES delivery;
  `importCGMES(config; path, name)` applies the configuration's import
  settings in one call.
- [`load_sparlectra_config`](@ref): load a Sparlectra YAML configuration
  file into a [`SparlectraConfig`](@ref).
- [`createNetFromMatPowerFile`](@ref): parse a MATPOWER `.m`/`.jl` case
  file and build the network;
  [`Sparlectra.MatpowerIO.read_case`](@ref) and
  [`Sparlectra.createNetFromMatPowerCase`](@ref) separate parsing from
  construction.
- [`run_sparlectra_api`](@ref): the non-interactive service contract
  below, for local applications and GUI integrations.
- [`runContingencies!`](@ref) / [`runScenarios!`](@ref): the N-1 and
  scenario batch entries. Screening defaults to `:off`; the service and
  the Web UI follow
  `contingency.screening.mode`, and `:flag` is an opt-in; see
  [N-1 Contingency Analysis](contingency.md).

## Service API: run_sparlectra_api

`run_sparlectra_api` takes a MATPOWER case, a configuration template
(never modified), an output directory and dotted-key overrides.

```julia
using Sparlectra
using SparlectraApp   # service API: environment app/, or app/ on the load path

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

The API never reads from stdin or asks for confirmation.

## Result contract

`SparlectraApiResult` separates framework status from transport metadata:

- `run_id`: a globally unique UUID string, stable throughout one run.
- `schema_version`: the serialized result contract, currently `"1.0"`.
- `status`, `success`, `converged`, `solution_available`: run state.
- `iterations`, `final_mismatch`, `reason`, `message`: the outcome.
- `casefile`, `config_file`, `output_dir`: effective paths.
- `logfile`, `result_file`, `artifacts`: generated files.
- `raw_result`: the underlying `SparlectraRunResult` for Julia callers.

Input and execution failures return `status == :failed` with a stable
reason such as `"casefile_not_found"`, `"invalid_configuration"`,
`"invalid_config_override"` or `"execution_error"`.

## GUI-editable overrides

Only keys in `GUI_EDITABLE_CONFIG_KEYS` are accepted:

- `power_flow.method`
- `power_flow.tol`
- `power_flow.max_iter`
- `power_flow.autodamp`
- `power_flow.autodamp_min`
- `power_flow.qlimits.enabled`
- `power_flow.qlimits.enforcement_mode`
- `power_flow.wrong_branch_detection`
- `power_flow.start_mode.angle_mode`
- `power_flow.start_mode.voltage_mode`
- `power_flow.start_current_iteration.enabled`
- `power_flow.start_current_iteration.max_iter`
- `power_flow.start_current_iteration.tol`
- `power_flow.start_current_iteration.damping`
- `power_flow.start_current_iteration.accept_only_if_improved`
- `power_flow.start_current_iteration.min_improvement_factor`
- `power_flow.start_current_iteration.vm_min_pu`
- `power_flow.start_current_iteration.vm_max_pu`
- `power_flow.start_current_iteration.max_angle_step_deg`
- `power_flow.start_current_iteration.only_for_large_cases`
- `power_flow.merit.enabled`
- `power_flow.merit.armijo_c1`
- `power_flow.merit.fallback_max_mismatch`
- `output.logfile_results`
- `output.detailed_result_csv_write_mode`
- `output.detailed_result_csv_exporter`
- `output.detailed_result_csv_direct_threshold_buses`
- `benchmark.enabled`
- `benchmark.samples`
- `benchmark.seconds`

Unknown or non-editable keys, invalid types, unsupported enum values and
invalid ranges are rejected before execution.

## Effective configuration and artifacts

Every run with a valid configuration writes `output_dir/effective_config.yaml`
and `output_dir/run_metadata.yaml` (request/lifecycle metadata);
successful and failed calls also write `run.log` and `result.json`.
Artifact discovery classifies these and any generated CSV or report files;
each `SparlectraApiArtifact` carries an absolute path, MIME type,
existence flag, byte size, kind and description.

After a guarded current-iteration start pre-solve the metadata carries
`current_iteration_enabled`, `current_iteration_attempted`,
`current_iteration_accepted`, `current_iteration_iterations`,
`current_iteration_initial_mismatch`, `current_iteration_final_mismatch`,
`current_iteration_reason` and `current_iteration_artifact`
(`current_iteration_start.log`, a start-value diagnostic). With
`power_flow.merit.enabled = true` it also carries `merit_enabled`,
`merit_used_iterations`, `merit_fallback_count`,
`merit_active_set_skip_count`, `merit_initial`, `merit_final` and
`merit_linesearch_artifact` (`merit_linesearch.log`, a line-search
diagnostic). Newton-Raphson remains the power-flow solver in both cases.

## Serialization

```julia
dict_value = to_dict(result)
named_value = to_namedtuple(result)
json_text = to_json(result)
yaml_text = to_yaml(result)
```

All transport forms, including `result.json`, carry `run_id` (the lookup
key) and `schema_version` (check it before decoding fields added later).
`raw_result` is omitted by default (a solved `Net` is not a stable
JSON/YAML representation); pass `include_raw_result=true` for Julia-side
inspection.

Runnable example: `examples/powerflow/exp_programmatic_api.jl`.
