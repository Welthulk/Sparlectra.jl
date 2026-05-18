using Sparlectra
using Test

function _flatten_yaml_paths(dict::AbstractDict, prefix::String = "")
  paths = String[]
  for (k, v) in dict
    key = String(k)
    path = isempty(prefix) ? key : string(prefix, ".", key)
    push!(paths, path)
    if v isa AbstractDict
      append!(paths, _flatten_yaml_paths(v, path))
    end
  end
  return paths
end

function run_configuration_docs_tests()
  @testset "Configuration docs coverage" begin
    yaml = Sparlectra.load_yaml_dict(joinpath(@__DIR__, "..", "src", "configuration.yaml.example"))
    paths = _flatten_yaml_paths(yaml)
    docs_files = [
      joinpath(@__DIR__, "..", "docs", "src", "configuration.md"),
      joinpath(@__DIR__, "..", "docs", "src", "powerflow_configuration.md"),
      joinpath(@__DIR__, "..", "docs", "src", "matpower_import.md"),
      joinpath(@__DIR__, "..", "docs", "src", "state_estimation_configuration.md"),
      joinpath(@__DIR__, "..", "docs", "src", "performance_profiling.md"),
    ]
    docs = join(read.(docs_files, String), "\n")

    for path in paths
      @test occursin(path, docs)
    end

    for key in [
      "output.console_auto_profile", "output.console_diagnostics", "output.console_q_limit_events",
      "output.logfile_results", "output.logfile_diagnostics", "output.logfile_performance", "output.logfile_warnings",
      "performance.level", "benchmark.show_once_output", "matpower_import.auto_profile", "power_flow.qlimits.start_mode",
    ]
      @test occursin(key, docs)
      @test occursin("Allowed values", docs)
    end

    @test occursin("matpower_import.benchmark", docs)
    @test occursin("benchmark.enabled", docs)
  end
end

