using Sparlectra

function main(args = ARGS)
  config_file = Sparlectra.configuration_path_from_inputs(
    env_var = "SPARLECTRA_CONFIGURATION_YAML",
    fallback_paths = [Sparlectra.USER_SPARLECTRA_CONFIG_PATH],
  )
  plot_curve = any(==("--plot"), args) ? true : (any(==("--no-plot"), args) ? false : nothing)
  return Sparlectra.run_voltage_dependent_control_demo(; config_file = config_file, plot_curve = plot_curve)
end

if get(ENV, "SPARLECTRA_SUITE_NO_AUTORUN", "0") != "1"
  Base.invokelatest(main, ARGS)
end
