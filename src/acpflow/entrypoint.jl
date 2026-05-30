"""
    run_sparlectra(; net=nothing, casefile=nothing, path=nothing, config=nothing, performance_profile=nothing) -> SparlectraRunResult

Run the configuration-driven Sparlectra framework workflow. Pass exactly one of
`net` or `casefile`. File-based runs import a MATPOWER case using
`config.matpower`; in-memory runs solve the supplied network directly. Solver,
control-loop, and output behavior are taken from one `SparlectraConfig`.
"""
function run_sparlectra(; net::Union{Nothing,Net} = nothing, casefile::Union{Nothing,String} = nothing, path::Union{Nothing,String} = nothing, config::Union{Nothing,SparlectraConfig} = nothing, performance_profile = nothing)::SparlectraRunResult
  return _run_sparlectra(; net = net, casefile = casefile, path = path, config = config, performance_profile = performance_profile)
end

"""
    run_acpflow(; net=nothing, casefile=nothing, path=nothing, config=nothing, performance_profile=nothing) -> SparlectraRunResult

Compatibility alias for [`run_sparlectra`](@ref).

This function uses the same configuration-driven framework workflow as
`run_sparlectra`. It does not support the former keyword-heavy runner API.
Solver, import, output, Q-limit, start-mode, and diagnostic behavior must be
configured through `SparlectraConfig` or YAML.
"""
function run_acpflow(;
  net::Union{Nothing,Net} = nothing,
  casefile::Union{Nothing,String} = nothing,
  path::Union{Nothing,String} = nothing,
  config::Union{Nothing,SparlectraConfig} = nothing,
  performance_profile = nothing,
)::SparlectraRunResult
  return run_sparlectra(;
    net = net,
    casefile = casefile,
    path = path,
    config = config,
    performance_profile = performance_profile,
  )
end

function _run_sparlectra(; net::Union{Nothing,Net} = nothing, casefile::Union{Nothing,String} = nothing, path::Union{Nothing,String} = nothing, config::Union{Nothing,SparlectraConfig} = nothing, performance_profile = nothing, emit_output::Bool = true)::SparlectraRunResult
  (net === nothing) == (casefile === nothing) && throw(ArgumentError("run_sparlectra: pass exactly one of net or casefile."))
  net !== nothing && path !== nothing && throw(ArgumentError("run_sparlectra: path is only valid with casefile."))
  cfg = config === nothing ? active_sparlectra_config() : config
  if net === nothing
    run_net = _import_sparlectra_net(casefile::String, path, cfg; performance_profile = performance_profile)
    run_cfg = _resolve_matpower_powerflow_ids_after_import(run_net, cfg; verbose = _runner_verbose(cfg))
  else
    run_net = net
    run_cfg = cfg
  end
  execution = _execute_sparlectra_powerflow!(run_net, run_cfg; performance_profile = performance_profile)
  result = _build_sparlectra_result(run_net, run_cfg, execution, performance_profile)
  return _postprocess_sparlectra_result!(result, run_cfg; emit_output = emit_output)
end
