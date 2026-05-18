#!/usr/bin/env julia
using Sparlectra

function main(args = ARGS)
  config_file = Sparlectra.configuration_path_from_inputs(
    env_var = "SPARLECTRA_CONFIGURATION_YAML",
    fallback_paths = [Sparlectra.USER_SPARLECTRA_CONFIG_PATH],
  )
  return Sparlectra.run_synthetic_tiled_grid_pf_perf(; config_file = config_file, args = collect(String, args))
end

if get(ENV, "SPARLECTRA_SUITE_NO_AUTORUN", "0") != "1"
  Base.invokelatest(getfield(@__MODULE__, :main), ARGS)
end
