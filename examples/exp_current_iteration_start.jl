using Sparlectra

function main(; casefile::AbstractString = "case5.m", output_dir::AbstractString = joinpath(@__DIR__, "_out", "current_iteration_start"))
  case_path = ensure_casefile(casefile)
  result = run_sparlectra_api(
    casefile = case_path,
    config_file = Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH,
    output_dir = output_dir,
    config_overrides = Dict(
      "power_flow.start_current_iteration.enabled" => true,
      "power_flow.start_current_iteration.max_iter" => 10,
      "power_flow.start_current_iteration.tol" => 1.0e-3,
      "power_flow.start_current_iteration.damping" => 0.5,
      "power_flow.start_current_iteration.accept_only_if_improved" => true,
      "power_flow.start_current_iteration.min_improvement_factor" => 0.98,
      "power_flow.start_current_iteration.vm_min_pu" => 0.5,
      "power_flow.start_current_iteration.vm_max_pu" => 1.5,
      "power_flow.start_current_iteration.max_angle_step_deg" => 30.0,
      "power_flow.start_current_iteration.only_for_large_cases" => false,
      "output.logfile_results" => "compact",
      "benchmark.enabled" => false,
    ),
  )

  println("run ID: ", result.run_id)
  println("status: ", result.status)
  println("current iteration attempted: ", get(result.metadata, "current_iteration_attempted", false))
  println("current iteration accepted: ", get(result.metadata, "current_iteration_accepted", false))
  println("current iteration reason: ", get(result.metadata, "current_iteration_reason", ""))
  for artifact in result.artifacts
    artifact.name == "current_iteration_start.log" && println("diagnostic artifact: ", artifact.path)
  end
  return result
end

Base.invokelatest(main)
