using Sparlectra

function main(;
  casefile::AbstractString = "case5.m",
  output_root::AbstractString = joinpath(@__DIR__, "_out", "powerflow_service"),
)
  request = Dict(
    "casefile" => ensure_casefile(casefile),
    "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH,
    "output_root" => output_root,
    "config_overrides" => Dict(
      "power_flow.tol" => 1.0e-8,
      "power_flow.max_iter" => 80,
      "benchmark.enabled" => false,
    ),
  )

  result = start_powerflow_run(request)
  result["success"] || return result

  run_id = result["run_id"]
  stored_result = get_powerflow_result(run_id)
  artifacts = list_powerflow_artifacts(run_id)

  println("run ID: ", run_id)
  println("stored status: ", stored_result["status"])
  for artifact in artifacts
    println(artifact["name"], " (", artifact["mime_type"], ")")
  end
  return result
end

Base.invokelatest(main)
