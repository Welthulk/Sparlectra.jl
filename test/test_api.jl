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
      csv_writer_path = joinpath(tmpdir, "writer.csv")
      Sparlectra._write_namedtuple_csv(
        csv_writer_path,
        [(name = "quoted, \"value\"", empty_missing = Base.missing, empty_nothing = nothing)],
        (:name, :empty_missing, :empty_nothing),
      )
      @test read(csv_writer_path, String) == "name,empty_missing,empty_nothing\n\"quoted, \"\"value\"\"\",,\n"
      semicolon_writer_path = joinpath(tmpdir, "writer_semicolon.csv")
      Sparlectra._write_namedtuple_csv(
        semicolon_writer_path,
        [(name = "quoted; \"value\"", empty_missing = Base.missing, empty_nothing = nothing)],
        (:name, :empty_missing, :empty_nothing);
        delimiter = ';',
      )
      @test read(semicolon_writer_path, String) == "name;empty_missing;empty_nothing\n\"quoted; \"\"value\"\"\";;\n"
      number_rows = [(integer = 1000000, decimal = 1234.56, negative = -1234.5, nan = NaN, inf = Inf, neginf = -Inf)]
      technical_path = joinpath(tmpdir, "writer_technical.csv")
      Sparlectra._write_namedtuple_csv(technical_path, number_rows, propertynames(only(number_rows)); format = "technical")
      @test read(technical_path, String) == "integer,decimal,negative,nan,inf,neginf\n1000000,1234.56,-1234.5,NaN,Inf,-Inf\n"
      german_path = joinpath(tmpdir, "writer_german.csv")
      Sparlectra._write_namedtuple_csv(german_path, number_rows, propertynames(only(number_rows)); delimiter = ';', format = "excel_de")
      @test read(german_path, String) == "integer;decimal;negative;nan;inf;neginf\n1.000.000;1.234,56;-1.234,5;NaN;Inf;-Inf\n"
      us_path = joinpath(tmpdir, "writer_us.csv")
      Sparlectra._write_namedtuple_csv(us_path, number_rows, propertynames(only(number_rows)); format = "excel_us")
      @test read(us_path, String) == "integer,decimal,negative,nan,inf,neginf\n\"1,000,000\",\"1,234.56\",\"-1,234.5\",NaN,Inf,-Inf\n"
      @test_throws ArgumentError Sparlectra._resolve_detailed_csv_format("unknown")

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
          "output.logfile_results" => "full",
          "benchmark.enabled" => false,
        ),
        performance_timing = :compact,
        run_diagnostics = true,
        detailed_result_csv = true,
        detailed_result_csv_format = "excel_de",
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
      run_log = read(joinpath(output_dir, "run.log"), String)
      for marker in ("solver_time:", "representative_time:", "iterations:", "final_mismatch:", "final_status:", "final_outcome:")
        @test occursin(marker, run_log)
      end
      @test occursin("detailed_result_csv_format: excel_de", run_log)
      @test occursin("detailed_result_csv_delimiter: semicolon", run_log)
      @test occursin("detailed_result_csv_decimal_separator: comma", run_log)
      @test occursin("detailed_result_csv_thousands_separator: dot", run_log)
      performance_log = read(joinpath(output_dir, "performance.log"), String)
      @test occursin("Sparlectra single-run phase timing", performance_log)
      @test occursin("api_config_build:", performance_log)
      @test occursin("case_loading_network_solver:", performance_log)
      @test occursin("total:", performance_log)
      diagnostic_log = read(joinpath(output_dir, "diagnose.log"), String)
      @test occursin("Sparlectra PowerFlow diagnostics", diagnostic_log)
      @test occursin("Final limit validation:", diagnostic_log)
      @test !isfile(joinpath(output_dir, "diagnose.txt"))

      bus_csv = read(joinpath(output_dir, "bus_voltages_complex.csv"), String)
      branch_csv = read(joinpath(output_dir, "branch_flows.csv"), String)
      @test startswith(bus_csv, "bus;bus_name;type;vm_pu;va_deg;vn_kV;v_re;v_im;v_complex;v_kV")
      @test startswith(branch_csv, "branch;branch_index;from_bus;to_bus;status;p_from_MW;q_from_MVar;p_to_MW;q_to_MVar")
      @test occursin(r"\d,\d", bus_csv)
      @test length(collect(eachline(IOBuffer(bus_csv)))) > 1
      @test length(collect(eachline(IOBuffer(branch_csv)))) > 1

      kinds = Set(artifact.kind for artifact in result.artifacts)
      @test :log in kinds
      @test :result_json in kinds
      @test :effective_config in kinds
      @test count(artifact -> artifact.kind === :csv && artifact.mime_type == "text/csv", result.artifacts) == 2
      @test all(artifact -> artifact.exists && artifact.size_bytes !== nothing, result.artifacts)
      @test result.logfile == only(artifact.path for artifact in result.artifacts if artifact.name == "run.log")
      @test result.result_file == only(artifact.path for artifact in result.artifacts if artifact.kind === :result_json)

      classic_dir = joinpath(tmpdir, "classic")
      full_dir = joinpath(tmpdir, "full")
      classic = run_sparlectra_api(
        casefile = casefile,
        config_file = template,
        output_dir = classic_dir,
        config_overrides = Dict("output.logfile_results" => "classic", "benchmark.enabled" => false),
      )
      full = run_sparlectra_api(
        casefile = casefile,
        config_file = template,
        output_dir = full_dir,
        config_overrides = Dict("output.logfile_results" => "full", "benchmark.enabled" => false),
        run_diagnostics = true,
      )
      @test classic.success && full.success
      @test !isfile(joinpath(classic_dir, "bus_voltages_complex.csv"))
      @test !isfile(joinpath(classic_dir, "branch_flows.csv"))
      classic_log = read(joinpath(classic_dir, "run.log"), String)
      full_log = read(joinpath(full_dir, "run.log"), String)
      @test !occursin("Full run details", classic_log)
      @test occursin("Full run details", full_log)
      @test occursin("Effective Sparlectra Configuration", full_log)
      @test occursin("diagnostics_artifact: diagnose.log", full_log)
      @test occursin("detailed_result_csv_delimiter: disabled", full_log)
      @test ncodeunits(full_log) > ncodeunits(classic_log)

      legacy_dir = joinpath(tmpdir, "legacy-semicolon")
      legacy = run_sparlectra_api(
        casefile = casefile,
        config_file = template,
        output_dir = legacy_dir,
        config_overrides = Dict("benchmark.enabled" => false),
        detailed_result_csv = true,
        detailed_result_csv_semicolon = true,
      )
      @test legacy.success
      @test startswith(read(joinpath(legacy_dir, "bus_voltages_complex.csv"), String), "bus;")

      failed_diagnostic = joinpath(tmpdir, "failed_diagnose.txt")
      Sparlectra._write_powerflow_diagnostics(failed_diagnostic, result.raw_result; diagnostic_fn = (io, _) -> error("diagnostic test failure"))
      failed_diagnostic_text = read(failed_diagnostic, String)
      @test occursin("Diagnostic generation failed", failed_diagnostic_text)
      @test occursin("diagnostic test failure", failed_diagnostic_text)

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

      @testset "Web UI case resolution" begin
        existing_m = _write_api_test_case(joinpath(tmpdir, "existing_case.m"))
        generated_jl = Sparlectra._resolve_powerflow_casefile(existing_m, joinpath(tmpdir, "cases"))
        @test lowercase(splitext(generated_jl)[2]) == ".jl"
        @test isfile(generated_jl)

        existing_jl = joinpath(tmpdir, "existing_case.jl")
        write(existing_jl, "nothing\n")
        @test Sparlectra._resolve_powerflow_casefile(existing_jl, joinpath(tmpdir, "cases")) == abspath(existing_jl)

        sibling_m = joinpath(tmpdir, "sibling_case.m")
        sibling_jl = joinpath(tmpdir, "sibling_case.jl")
        write(sibling_m, "not parsed when the Julia case exists\n")
        write(sibling_jl, "nothing\n")
        fail_emit = (_, _) -> error("emit should not be called")
        @test Sparlectra._resolve_powerflow_casefile(sibling_m, joinpath(tmpdir, "cases"); emit_julia_case_fn = fail_emit) == abspath(sibling_jl)

        fallback_m = joinpath(tmpdir, "fallback_case.m")
        write(fallback_m, "invalid MATPOWER used to exercise conversion fallback\n")
        fallback_resolved = @test_logs (:warn, r"MATPOWER Julia case generation failed") Sparlectra._resolve_powerflow_casefile(fallback_m, joinpath(tmpdir, "cases"); emit_julia_case_fn = (_, _) -> error("conversion failed"))
        @test fallback_resolved == abspath(fallback_m)

        ensure_calls = NamedTuple[]
        fake_ensure = function(requested; outdir, to_jl)
          push!(ensure_calls, (; requested, outdir, to_jl))
          mfile = joinpath(outdir, first(splitext(requested)) * ".m")
          jlfile = joinpath(outdir, first(splitext(requested)) * ".jl")
          write(mfile, "downloaded MATPOWER placeholder\n")
          write(jlfile, "nothing\n")
          return endswith(lowercase(requested), ".jl") ? jlfile : mfile
        end
        case_directory = joinpath(tmpdir, "downloaded_cases")
        resolved_missing = Sparlectra._resolve_powerflow_casefile("case118.m", case_directory; ensure_casefile_fn = fake_ensure, emit_julia_case_fn = fail_emit)
        @test resolved_missing == abspath(joinpath(case_directory, "case118.jl"))
        @test only(ensure_calls) == (requested = "case118.m", outdir = abspath(case_directory), to_jl = true)

        generation_failure_ensure = function(requested; outdir, to_jl)
          @test to_jl
          mfile = joinpath(outdir, first(splitext(requested)) * ".m")
          write(mfile, "downloaded MATPOWER fallback\n")
          error("Julia case generation failed")
        end
        downloaded_fallback = @test_logs (:warn, r"using the downloaded \.m case") Sparlectra._resolve_powerflow_casefile("case300.jl", case_directory; ensure_casefile_fn = generation_failure_ensure)
        @test downloaded_fallback == abspath(joinpath(case_directory, "case300.m"))

        empty!(ensure_calls)
        resolved_requested_jl = Sparlectra._resolve_powerflow_casefile("case14.jl", case_directory; ensure_casefile_fn = fake_ensure, emit_julia_case_fn = fail_emit)
        @test resolved_requested_jl == abspath(joinpath(case_directory, "case14.jl"))
        @test only(ensure_calls).to_jl

        empty!(ensure_calls)
        @test_throws ArgumentError Sparlectra._resolve_powerflow_casefile(joinpath("missing", "case14.m"), case_directory; ensure_casefile_fn = fake_ensure)
        @test_throws ArgumentError Sparlectra._resolve_powerflow_casefile("https://example.com/case14.m", case_directory; ensure_casefile_fn = fake_ensure)
        @test isempty(ensure_calls)

        rejected_path = start_powerflow_run(Dict("casefile" => joinpath("missing", "case14.m")); case_directory)
        @test rejected_path["reason"] == "invalid_casefile"
        rejected_url = start_powerflow_run(Dict("casefile" => "https://example.com/case14.m"); case_directory)
        @test rejected_url["reason"] == "invalid_casefile"
      end

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

      resolved_by_service = Ref(false)
      service_resolver = function(requested, directory)
        resolved_by_service[] = requested == "case_service.m" && directory == joinpath(tmpdir, "trusted_cases")
        return abspath(casefile)
      end
      resolved_run = start_powerflow_run(Dict(
        "casefile" => "case_service.m",
        "config_file" => config_file,
        "output_root" => output_root,
      ); case_directory = joinpath(tmpdir, "trusted_cases"), case_resolver = service_resolver)
      @test resolved_run["success"]
      @test resolved_by_service[]
      @test resolved_run["casefile"] == abspath(casefile)

      index_path = joinpath(output_root, POWERFLOW_RUN_INDEX_FILENAME)
      @test isfile(index_path)
      index = load_powerflow_run_index(output_root)
      @test index["schema_version"] == "1.0"
      @test any(run -> run["run_id"] == started["run_id"] && haskey(run, "timestamp"), index["runs"])
      listed_runs = list_powerflow_runs(output_root)
      @test any(run -> run["run_id"] == started["run_id"] && run["available"] && haskey(run, "timestamp"), listed_runs)
      @test isempty(list_powerflow_runs(joinpath(tmpdir, "empty_service")))

      entered = Channel{Nothing}(1)
      release = Channel{Nothing}(1)
      controlled_runner = function(request; case_directory = nothing)
        put!(entered, nothing)
        token = request["cancellation_token"]
        while !isready(release)
          Sparlectra._check_powerflow_cancelled!(token)
          yield()
        end
        take!(release)
        return start_powerflow_run(request; case_directory)
      end
      async_request = Dict(
        "casefile" => casefile,
        "config_file" => config_file,
        "output_root" => joinpath(tmpdir, "async-runs"),
      )
      active = Sparlectra.start_webui_powerflow_run(async_request; runner = controlled_runner)
      take!(entered)
      @test Sparlectra.get_webui_powerflow_job(active["run_id"])["status"] == "running"
      @test Sparlectra.get_active_webui_powerflow_job()["run_id"] == active["run_id"]
      rejected_concurrent = Sparlectra.start_webui_powerflow_run(async_request; runner = controlled_runner)
      @test rejected_concurrent["reason"] == "active_run"
      aborted = Sparlectra.abort_webui_powerflow_run(active["run_id"])
      @test aborted["status"] == "aborting"
      @test aborted["abort_status"] == "accepted"
      @test !aborted["success"]
      @test Sparlectra.abort_webui_powerflow_run(active["run_id"])["abort_status"] == "already_aborting"
      @test Sparlectra.start_webui_powerflow_run(async_request; runner = controlled_runner)["reason"] == "active_run"
      @test Sparlectra.abort_webui_powerflow_run("../unsafe")["reason"] == "unsafe_run_id"
      wait(Sparlectra._POWERFLOW_WEBUI_JOBS[active["run_id"]]["task"])
      @test get_powerflow_result(active["run_id"])["status"] == "aborted"
      @test !get_powerflow_result(active["run_id"])["success"]
      @test occursin("Run aborted by user.", read(joinpath(aborted["output_dir"], "run.log"), String))
      @test Sparlectra.abort_webui_powerflow_run(active["run_id"])["abort_status"] == "already_aborted"
      replacement = Sparlectra.start_webui_powerflow_run(async_request; runner = controlled_runner)
      @test replacement["status"] in ("queued", "running")
      put!(release, nothing)
      wait(Sparlectra._POWERFLOW_WEBUI_JOBS[replacement["run_id"]]["task"])
      @test Sparlectra.get_webui_powerflow_job(active["run_id"])["status"] == "aborted"
      @test Sparlectra.get_active_webui_powerflow_job() === nothing

      pre_cancelled = Threads.Atomic{Bool}(true)
      @test_throws Sparlectra.PowerFlowAborted start_powerflow_run(merge(async_request, Dict("cancellation_token" => pre_cancelled)))

      solver_checks = Ref(0)
      solver_profile = Dict{Symbol,Any}(
        :cancellation_check => () -> begin
          solver_checks[] += 1
          solver_checks[] >= 5 && throw(Sparlectra.PowerFlowAborted())
        end,
      )
      solver_net = deepcopy(Sparlectra._registered_powerflow_run(started["run_id"]).raw_result.net)
      @test_throws Sparlectra.PowerFlowAborted runpf!(solver_net; config = powerflow_config(), performance_profile = solver_profile)
      @test solver_checks[] == 5

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
      @test Set(refresh["loaded_runs"]) == Set([started["run_id"], resolved_run["run_id"], failed["run_id"]])
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
  @testset "PowerFlow run deletion safety" begin
    mktempdir() do tmpdir
      output_root = joinpath(tmpdir, "runs")
      outside_dir = joinpath(tmpdir, "outside")
      mkpath(outside_dir)
      outside_file = joinpath(outside_dir, "keep.txt")
      write(outside_file, "keep")

      entries = Any[]
      for run_id in ("run-a", "run-b")
        run_dir = joinpath(output_root, run_id)
        mkpath(run_dir)
        result_file = joinpath(run_dir, "result.json")
        write(result_file, "{}")
        push!(entries, Dict(
          "run_id" => run_id,
          "status" => "succeeded",
          "success" => true,
          "output_dir" => run_dir,
          "result_file" => result_file,
        ))
      end
      push!(entries, Dict(
        "run_id" => "unsafe-index-entry",
        "output_dir" => outside_dir,
        "result_file" => joinpath(outside_dir, "result.json"),
      ))
      Sparlectra._write_powerflow_run_entries!(output_root, entries)

      registered = Sparlectra._api_result(
        run_id = "run-a",
        status = :succeeded,
        success = true,
        output_dir = joinpath(output_root, "run-a"),
        result_file = joinpath(output_root, "run-a", "result.json"),
      )
      Sparlectra._register_powerflow_run!(registered)
      deleted = delete_powerflow_run("run-a"; output_root)
      @test deleted["success"]
      @test !ispath(joinpath(output_root, "run-a"))
      @test get_powerflow_result("run-a")["reason"] == "run_not_found"
      @test all(entry -> get(entry, "run_id", "") != "run-a", load_powerflow_run_index(output_root)["runs"])
      @test isfile(outside_file)

      for unsafe_id in ("../outside", "/tmp/outside", raw"C:\outside", raw"C:\outside\run")
        rejected = delete_powerflow_run(unsafe_id; output_root)
        @test !rejected["success"]
        @test rejected["reason"] == "unsafe_run_id"
      end
      @test delete_powerflow_run("unknown"; output_root)["reason"] == "run_not_found"

      deleted_all = delete_all_powerflow_runs(; output_root)
      @test deleted_all["status"] == "partial"
      @test deleted_all["deleted_runs"] == ["run-b"]
      @test get(only(deleted_all["failed_runs"]), "reason", "") == "unsafe_output_dir"
      @test !ispath(joinpath(output_root, "run-b"))
      @test isfile(outside_file)
      remaining = load_powerflow_run_index(output_root)["runs"]
      @test length(remaining) == 1
      @test remaining[1]["run_id"] == "unsafe-index-entry"
    end
  end
  return nothing
end
