# Programmatic API

The stable programmatic surface of Sparlectra is a small set of entry
points; everything else on the [reference pages](reference.md) is internal
and may change between versions.

## Stable entry points

- [`run_sparlectra`](@ref): the configured power flow on a network or
  case file.
- [`importCGMES`](@ref): import a CGMES delivery;
  `importCGMES(config; path, name)` applies the configuration's import
  settings.
- [`load_sparlectra_config`](@ref): load a YAML configuration file into a
  [`SparlectraConfig`](@ref).
- [`createNetFromMatPowerFile`](@ref): the network from a MATPOWER
  `.m`/`.jl` case file; [`Sparlectra.MatpowerIO.read_case`](@ref) and
  [`Sparlectra.createNetFromMatPowerCase`](@ref) split parsing and
  construction.
- [`run_sparlectra_api`](@ref): the service contract below.
- [`runContingencies!`](@ref) / [`runScenarios!`](@ref): the N-1 and
  scenario batches ([N-1 Contingency Analysis](contingency.md)).

## Service API: run_sparlectra_api

`run_sparlectra_api` takes a MATPOWER case, a configuration template
(never modified), an output directory and dotted-key overrides; it never
reads from stdin or asks for confirmation.

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

## Result contract

`SparlectraApiResult` fields:

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

Only keys in `GUI_EDITABLE_CONFIG_KEYS` (`src/config/config_overrides.jl`)
are accepted: the `power_flow.*` solver, start, Q-limit, merit,
trust-region, rescue and slack keys, `contingency.*`, `output.*`,
`benchmark.*`, `runtime.parallel.enabled`, the import keys of `model.*`,
`matpower_import.*`, `matpower_export.*`, `cgmes_import.*` and
`powsybl_import.*`, `state_estimation.*` and `short_circuit.sweep_method`.
Anything else, and invalid types, enum values or ranges, is rejected
before execution.

## Effective configuration and artifacts

Every run with a valid configuration writes `effective_config.yaml` and
`run_metadata.yaml` into `output_dir`, successful and failed calls also
`run.log` and `result.json` (contents:
[Local PowerFlow Service](powerflow_service.md)). Each
`SparlectraApiArtifact` carries path, MIME type, existence flag, byte
size, kind and description. The `current_iteration_*` and `merit_*`
metadata keys (artifacts `current_iteration_start.log`,
`merit_linesearch.log`) are listed on
[Power-Flow Configuration](powerflow_configuration.md).

## Serialization

```julia
dict_value = to_dict(result)
named_value = to_namedtuple(result)
json_text = to_json(result)
yaml_text = to_yaml(result)
```

Every transport form, `result.json` included, carries `run_id` and
`schema_version` (check it before decoding fields added later);
`raw_result` is omitted unless `include_raw_result=true`.

Runnable example: `examples/powerflow/exp_programmatic_api.jl`.
