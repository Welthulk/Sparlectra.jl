using Sparlectra

function main(; casefile::AbstractString = "case5.m", output_dir::AbstractString = joinpath(@__DIR__, "_out", "api_run"))
  case_path = ensure_casefile(casefile)
  result = run_sparlectra_api(
    casefile = case_path,
    config_file = Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH,
    output_dir = output_dir,
    config_overrides = Dict(
      "power_flow.tol" => 1.0e-8,
      "power_flow.max_iter" => 80,
      "power_flow.autodamp" => true,
      "output.logfile_results" => "compact",
      "benchmark.enabled" => false,
    ),
  )

  println("status: ", result.status)
  println("solution available: ", result.solution_available)
  for artifact in result.artifacts
    println(artifact.kind, ": ", artifact.path)
  end
  return result
end

Base.invokelatest(main)
