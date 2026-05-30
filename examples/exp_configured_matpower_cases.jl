#!/usr/bin/env julia
using Sparlectra

function main(args = ARGS)
  config_file = Sparlectra.configuration_path_from_inputs(
    env_var = "SPARLECTRA_CONFIGURATION_YAML",
    fallback_paths = [Sparlectra.USER_SPARLECTRA_CONFIG_PATH, Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH],
  )
  cfg = Sparlectra.load_sparlectra_config(config_file)
  results = Sparlectra.run_sparlectra_cases(config = cfg)
  for result in results
    println(result.net.name, ": ", result.outcome)
  end
  return results
end

if get(ENV, "SPARLECTRA_SUITE_NO_AUTORUN", "0") != "1"
  Base.invokelatest(getfield(@__MODULE__, :main), ARGS)
end
