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
      @test !isempty(result.run_id)
      @test result.schema_version == "1.0"
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
      @test dict_result["run_id"] == result.run_id
      @test dict_result["schema_version"] == "1.0"
      @test dict_result["status"] == "succeeded"
      @test !haskey(dict_result, "raw_result")
      @test named_result.run_id == result.run_id
      @test named_result.schema_version == "1.0"
      @test named_result.status == "succeeded"
      @test occursin("\"run_id\":\"$(result.run_id)\"", json_result)
      @test occursin("\"schema_version\":\"1.0\"", json_result)
      @test occursin("run_id: $(result.run_id)", yaml_result)
      @test occursin("schema_version: 1.0", yaml_result)
      @test occursin("\"status\":\"succeeded\"", json_result)
      @test occursin("status: succeeded", yaml_result)
      result_file_text = read(result.result_file, String)
      @test occursin("\"run_id\":\"$(result.run_id)\"", result_file_text)
      @test occursin("\"schema_version\":\"1.0\"", result_file_text)
      @test occursin("\"artifacts\"", result_file_text)

      write(joinpath(output_dir, "buses.csv"), "bus,vm\n1,1.0\n")
      write(joinpath(output_dir, "report.txt"), "report\n")
      discovered = collect_sparlectra_api_artifacts(output_dir)
      @test :csv in Set(artifact.kind for artifact in discovered)
      @test :report in Set(artifact.kind for artifact in discovered)

      missing = run_sparlectra_api(casefile = joinpath(tmpdir, "missing.m"), config_file = template, output_dir = joinpath(tmpdir, "missing"))
      @test !missing.success
      @test missing.status === :failed
      @test missing.reason == "casefile_not_found"
      @test !isempty(missing.run_id)
      @test missing.run_id != result.run_id
      @test missing.schema_version == "1.0"
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
  @testset "Local PowerFlow service" begin
    mktempdir() do tmpdir
      casefile = _write_api_test_case(joinpath(tmpdir, "case_service.m"))
      config_file = joinpath(tmpdir, "service_config.yaml")
      cp(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, config_file)
      output_root = joinpath(tmpdir, "powerflow_service")

      started = start_powerflow_run(Dict(
        "casefile" => casefile,
        "config_file" => config_file,
        "output_root" => output_root,
        "config_overrides" => Dict(
          "power_flow.tol" => 1.0e-8,
          "power_flow.max_iter" => 80,
          "benchmark.enabled" => false,
        ),
      ))

      @test started["success"]
      required_result_fields = ("run_id", "schema_version", "status", "success", "converged", "solution_available", "iterations", "final_mismatch", "reason", "message", "artifacts")
      @test all(field -> haskey(started, field), required_result_fields)
      @test !isempty(started["run_id"])
      @test started["output_dir"] == joinpath(abspath(output_root), started["run_id"])
      @test isfile(joinpath(started["output_dir"], "result.json"))
      @test isfile(joinpath(started["output_dir"], "effective_config.yaml"))

      index_path = joinpath(output_root, POWERFLOW_RUN_INDEX_FILENAME)
      @test isfile(index_path)
      index = load_powerflow_run_index(output_root)
      @test index["schema_version"] == "1.0"
      @test any(run -> run["run_id"] == started["run_id"], index["runs"])
      listed_runs = list_powerflow_runs(output_root)
      @test any(run -> run["run_id"] == started["run_id"] && run["available"], listed_runs)
      @test isempty(list_powerflow_runs(joinpath(tmpdir, "empty_service")))

      failed = start_powerflow_run(Dict(
        "casefile" => casefile,
        "config_file" => config_file,
        "output_root" => output_root,
        "config_overrides" => Dict("power_flow.tol" => -1.0),
      ))
      @test !failed["success"]
      @test isfile(joinpath(failed["output_dir"], "result.json"))
      @test any(run -> run["run_id"] == failed["run_id"], load_powerflow_run_index(output_root)["runs"])

      stored = get_powerflow_result(started["run_id"])
      @test stored["run_id"] == started["run_id"]
      @test stored["schema_version"] == started["schema_version"]

      artifacts = list_powerflow_artifacts(started["run_id"])
      artifact_names = Set(artifact["name"] for artifact in artifacts)
      @test Set(["result.json", "run.log", "effective_config.yaml"]) ⊆ artifact_names

      resolved = resolve_powerflow_artifact(started["run_id"], "result.json")
      @test resolved isa SparlectraApiArtifact
      @test isfile(resolved.path)
      @test !isempty(resolved.mime_type)

      lock(Sparlectra._POWERFLOW_SERVICE_LOCK) do
        empty!(Sparlectra._POWERFLOW_SERVICE_RUNS)
      end
      @test get_powerflow_result(started["run_id"])["reason"] == "run_not_found"
      refresh = refresh_powerflow_run_registry!(output_root)
      @test refresh["status"] == "succeeded"
      @test Set(refresh["loaded_runs"]) == Set([started["run_id"], failed["run_id"]])
      @test get_powerflow_result(started["run_id"])["run_id"] == started["run_id"]
      recovered_artifacts = list_powerflow_artifacts(started["run_id"])
      recovered_names = Set(artifact["name"] for artifact in recovered_artifacts)
      @test Set(["result.json", "run.log", "effective_config.yaml"]) ⊆ recovered_names
      recovered_result = resolve_powerflow_artifact(started["run_id"], "result.json")
      @test recovered_result isa SparlectraApiArtifact
      @test recovered_result.path == joinpath(started["output_dir"], "result.json")

      for unsafe_name in ("../result.json", "../../Project.toml", "/etc/passwd", raw"C:\Users\x\file.txt")
        rejected = resolve_powerflow_artifact(started["run_id"], unsafe_name)
        @test rejected["success"] === false
        @test rejected["reason"] == "unsafe_artifact_name"
      end

      missing_artifact = resolve_powerflow_artifact(started["run_id"], "missing.txt")
      @test missing_artifact["reason"] == "artifact_not_found"

      missing_run = get_powerflow_result("unknown-run-id")
      @test missing_run["success"] === false
      @test missing_run["reason"] == "run_not_found"
      @test list_powerflow_artifacts("unknown-run-id")["reason"] == "run_not_found"
      @test resolve_powerflow_artifact("unknown-run-id", "result.json")["reason"] == "run_not_found"

      corrupt_id = "corrupt-run"
      corrupt_dir = joinpath(output_root, corrupt_id)
      mkpath(corrupt_dir)
      write(joinpath(corrupt_dir, "result.json"), "{not valid json")
      unsafe_id = "unsafe-run"
      missing_id = "missing-run"
      modified_index = load_powerflow_run_index(output_root)
      push!(modified_index["runs"], Dict(
        "run_id" => corrupt_id,
        "output_dir" => corrupt_dir,
        "result_file" => joinpath(corrupt_dir, "result.json"),
      ))
      push!(modified_index["runs"], Dict(
        "run_id" => unsafe_id,
        "output_dir" => tmpdir,
        "result_file" => joinpath(tmpdir, "result.json"),
      ))
      push!(modified_index["runs"], Dict(
        "run_id" => missing_id,
        "output_dir" => joinpath(output_root, missing_id),
        "result_file" => joinpath(output_root, missing_id, "result.json"),
      ))
      open(index_path, "w") do io
        Sparlectra._write_json(io, modified_index)
        println(io)
      end

      lock(Sparlectra._POWERFLOW_SERVICE_LOCK) do
        empty!(Sparlectra._POWERFLOW_SERVICE_RUNS)
      end
      partial_refresh = refresh_powerflow_run_registry!(output_root)
      @test partial_refresh["status"] == "partial"
      @test get_powerflow_result(started["run_id"])["run_id"] == started["run_id"]
      unavailable_reasons = Dict(run["run_id"] => run["reason"] for run in partial_refresh["unavailable_runs"])
      @test unavailable_reasons[corrupt_id] == "invalid_result_file"
      @test unavailable_reasons[unsafe_id] == "unsafe_output_dir"
      @test unavailable_reasons[missing_id] == "result_file_not_found"
      listed_by_id = Dict(run["run_id"] => run for run in list_powerflow_runs(output_root) if haskey(run, "run_id"))
      @test listed_by_id[corrupt_id]["available"]
      @test !listed_by_id[unsafe_id]["available"]
      @test listed_by_id[unsafe_id]["reason"] == "unsafe_output_dir"
      @test !listed_by_id[missing_id]["available"]
      @test listed_by_id[missing_id]["reason"] == "output_dir_not_found"
    end
  end
  return nothing
end
