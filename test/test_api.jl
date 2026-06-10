using Sparlectra
using Test

function _write_api_test_case(path::AbstractString)
  write(path, """
function mpc = case_api
mpc.version = '2';
mpc.baseMVA = 100;
mpc.bus = [
1 3 0 0 0 0 1 1.0 0 110 1 1.1 0.9;
2 1 40 15 0 0 1 1.0 0 110 1 1.1 0.9;
];
mpc.gen = [
1 100 0 300 -300 1.02 100 1 300 0;
];
mpc.branch = [
1 2 0.01 0.05 0.0 999 999 999 0 0 1 -360 360;
];
""")
  return path
end

function run_api_tests()
  @testset "GUI-ready Sparlectra API" begin
    mktempdir() do tmpdir
      casefile = _write_api_test_case(joinpath(tmpdir, "case_api.m"))
      template = joinpath(tmpdir, "config_template.yaml")
      cp(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, template)
      template_before = read(template, String)
      output_dir = joinpath(tmpdir, "success")

      result = run_sparlectra_api(
        casefile = casefile,
        config_file = template,
        output_dir = output_dir,
        config_overrides = Dict(
          "power_flow.tol" => 1.0e-9,
          "power_flow.max_iter" => 40,
          "power_flow.autodamp" => true,
          "output.logfile_results" => "compact",
          "benchmark.enabled" => false,
        ),
      )

      @test result isa SparlectraApiResult
      @test result.status === :succeeded
      @test result.success
      @test result.converged === true
      @test result.solution_available
      @test result.iterations isa Int
      @test result.final_mismatch isa Float64
      @test result.raw_result isa SparlectraRunResult
      @test read(template, String) == template_before
      @test isfile(joinpath(output_dir, "effective_config.yaml"))
      @test occursin("tol: 1.0e-9", read(joinpath(output_dir, "effective_config.yaml"), String))
      effective_cfg = Sparlectra.load_sparlectra_config(joinpath(output_dir, "effective_config.yaml"); reload = true)
      @test effective_cfg.powerflow.tol == 1.0e-9
      @test effective_cfg.powerflow.max_iter == 40

      kinds = Set(artifact.kind for artifact in result.artifacts)
      @test :log in kinds
      @test :result_json in kinds
      @test :effective_config in kinds
      @test all(artifact -> artifact.exists && artifact.size_bytes !== nothing, result.artifacts)
      @test result.logfile == only(artifact.path for artifact in result.artifacts if artifact.kind === :log)
      @test result.result_file == only(artifact.path for artifact in result.artifacts if artifact.kind === :result_json)

      dict_result = to_dict(result)
      named_result = to_namedtuple(result)
      json_result = to_json(result)
      yaml_result = to_yaml(result)
      @test dict_result["status"] == "succeeded"
      @test !haskey(dict_result, "raw_result")
      @test named_result.status == "succeeded"
      @test occursin("\"status\":\"succeeded\"", json_result)
      @test occursin("status: succeeded", yaml_result)
      @test occursin("\"artifacts\"", read(result.result_file, String))

      write(joinpath(output_dir, "buses.csv"), "bus,vm\n1,1.0\n")
      write(joinpath(output_dir, "report.txt"), "report\n")
      discovered = collect_sparlectra_api_artifacts(output_dir)
      @test :csv in Set(artifact.kind for artifact in discovered)
      @test :report in Set(artifact.kind for artifact in discovered)

      missing = run_sparlectra_api(casefile = joinpath(tmpdir, "missing.m"), config_file = template, output_dir = joinpath(tmpdir, "missing"))
      @test !missing.success
      @test missing.status === :failed
      @test missing.reason == "casefile_not_found"
      @test any(artifact -> artifact.kind === :effective_config, missing.artifacts)

      invalid_config = joinpath(tmpdir, "invalid.yaml")
      write(invalid_config, "power_flow:\n  typo_tol: 1.0e-8\n")
      invalid = run_sparlectra_api(casefile = casefile, config_file = invalid_config, output_dir = joinpath(tmpdir, "invalid"))
      @test !invalid.success
      @test invalid.reason == "invalid_configuration"
      @test invalid.message !== nothing

      invalid_override = run_sparlectra_api(
        casefile = casefile,
        config_file = template,
        output_dir = joinpath(tmpdir, "invalid_override"),
        config_overrides = Dict("power_flow.typo_tol" => 1.0e-8),
      )
      @test !invalid_override.success
      @test invalid_override.reason == "invalid_config_override"
      @test occursin("Unknown", something(invalid_override.message, ""))

      disallowed_override = run_sparlectra_api(
        casefile = casefile,
        config_file = template,
        output_dir = joinpath(tmpdir, "disallowed_override"),
        config_overrides = Dict("runtime.julia_threads" => "auto"),
      )
      @test !disallowed_override.success
      @test disallowed_override.reason == "invalid_config_override"
      @test occursin("not allowed", something(disallowed_override.message, ""))

      invalid_value = run_sparlectra_api(
        casefile = casefile,
        config_file = template,
        output_dir = joinpath(tmpdir, "invalid_value"),
        config_overrides = Dict("power_flow.max_iter" => 0),
      )
      @test !invalid_value.success
      @test invalid_value.reason == "invalid_config_override"
      @test occursin("positive", something(invalid_value.message, ""))

      invalid_type = run_sparlectra_api(
        casefile = casefile,
        config_file = template,
        output_dir = joinpath(tmpdir, "invalid_type"),
        config_overrides = Dict("power_flow.autodamp" => "true"),
      )
      @test !invalid_type.success
      @test invalid_type.reason == "invalid_config_override"
      @test occursin("invalid type", something(invalid_type.message, ""))

      invalid_enum = run_sparlectra_api(
        casefile = casefile,
        config_file = template,
        output_dir = joinpath(tmpdir, "invalid_enum"),
        config_overrides = Dict("power_flow.wrong_branch_detection" => "unknown"),
      )
      @test !invalid_enum.success
      @test invalid_enum.reason == "invalid_config_override"
      @test occursin("must be one of", something(invalid_enum.message, ""))

      @test validate_gui_config_overrides(Dict("power_flow.qlimits.enabled" => false))["power_flow"]["qlimits"]["enabled"] === false
    end
  end
  return nothing
end
