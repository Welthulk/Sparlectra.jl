using Sparlectra
using Test
using Dates

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

function _control_label_test_net()
  net = Net(name = "control_label_cache", baseMVA = 100.0)
  for name in ("Slack", "NoControl", "QControl", "PControl", "BothControl")
    addBus!(net = net, busName = name, vn_kV = 110.0)
  end
  addACLine!(net = net, fromBus = "Slack", toBus = "NoControl", length = 1.0, r = 0.01, x = 0.05)
  addACLine!(net = net, fromBus = "NoControl", toBus = "QControl", length = 1.0, r = 0.01, x = 0.05)
  addACLine!(net = net, fromBus = "QControl", toBus = "PControl", length = 1.0, r = 0.01, x = 0.05)
  addACLine!(net = net, fromBus = "PControl", toBus = "BothControl", length = 1.0, r = 0.01, x = 0.05)
  curve = make_characteristic([(0.9, 0.0), (1.1, 0.0)])
  qu = QUController(curve; qmin_pu = -1.0, qmax_pu = 1.0)
  pu = PUController(curve; pmin_pu = 0.0, pmax_pu = 1.0)
  addProsumer!(net = net, busName = "Slack", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "Slack")
  addProsumer!(net = net, busName = "NoControl", type = "ENERGYCONSUMER", p = 1.0, q = 0.2)
  addProsumer!(net = net, busName = "QControl", type = "SYNCHRONOUSMACHINE", p = 1.0, q = 0.0, qu_controller = qu)
  addProsumer!(net = net, busName = "PControl", type = "SYNCHRONOUSMACHINE", p = 1.0, q = 0.0, pu_controller = pu)
  addProsumer!(net = net, busName = "BothControl", type = "SYNCHRONOUSMACHINE", p = 1.0, q = 0.0, qu_controller = qu)
  addProsumer!(net = net, busName = "BothControl", type = "SYNCHRONOUSMACHINE", p = 1.0, q = 0.0, pu_controller = pu)
  return net
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
      buffered_writer_path = joinpath(tmpdir, "writer_buffered.csv")
      streaming_writer_path = joinpath(tmpdir, "writer_streaming.csv")
      csv_strategy_rows = [(name = "alpha", value = 1.25), (name = "beta, quoted", value = 2.5)]
      buffered_cfg = Sparlectra.OutputConfig(detailed_result_csv_write_mode = :buffered, detailed_result_csv_buffer_initial_bytes = 128)
      streaming_cfg = Sparlectra.OutputConfig(detailed_result_csv_write_mode = :streaming)
      Sparlectra._write_namedtuple_csv(buffered_writer_path, csv_strategy_rows, (:name, :value); config = buffered_cfg)
      Sparlectra._write_namedtuple_csv(streaming_writer_path, csv_strategy_rows, (:name, :value); config = streaming_cfg)
      @test read(buffered_writer_path, String) == "name,value\nalpha,1.25\n\"beta, quoted\",2.5\n"
      @test read(streaming_writer_path, String) == read(buffered_writer_path, String)
      auto_streaming_cfg = Sparlectra.OutputConfig(detailed_result_csv_write_mode = :auto, detailed_result_csv_streaming_threshold_rows = 1)
      @test Sparlectra._select_namedtuple_csv_write_mode(csv_strategy_rows, (:name, :value); config = auto_streaming_cfg, estimated_rows = 2) === :streaming
      safe_default_cfg = Sparlectra.OutputConfig(Dict("output" => Dict("detailed_result_csv_buffer_initial_bytes" => -1, "detailed_result_csv_buffer_max_bytes" => 0, "detailed_result_csv_streaming_threshold_rows" => 0)))
      @test safe_default_cfg.detailed_result_csv_buffer_initial_bytes == 8 * 1024 * 1024
      @test safe_default_cfg.detailed_result_csv_buffer_max_bytes == 64 * 1024 * 1024
      @test safe_default_cfg.detailed_result_csv_streaming_threshold_rows == 100_000
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
      exponent_rows = [(id = "1E5", small = 1.23e-5, large = 1.23e5)]
      exponent_de_path = joinpath(tmpdir, "writer_exponent_de.csv")
      Sparlectra._write_namedtuple_csv(exponent_de_path, exponent_rows, propertynames(only(exponent_rows)); delimiter = ';', format = "excel_de")
      exponent_de = read(exponent_de_path, String)
      @test exponent_de == "id;small;large\n1E5;0,0000123;123.000\n"
      @test !occursin("1,23e", lowercase(exponent_de))
      @test !occursin("=", exponent_de)
      exponent_us_path = joinpath(tmpdir, "writer_exponent_us.csv")
      Sparlectra._write_namedtuple_csv(exponent_us_path, exponent_rows, propertynames(only(exponent_rows)); format = "excel_us")
      exponent_us = read(exponent_us_path, String)
      @test exponent_us == "id,small,large\n1E5,0.0000123,\"123,000\"\n"
      @test !occursin("1.23e", lowercase(exponent_us))
      @test !occursin("=", exponent_us)
      exponent_technical_path = joinpath(tmpdir, "writer_exponent_technical.csv")
      Sparlectra._write_namedtuple_csv(exponent_technical_path, exponent_rows, propertynames(only(exponent_rows)); format = "technical")
      @test read(exponent_technical_path, String) == "id,small,large\n1E5,1.23e-05,123000\n"
      direct_cell_format = Sparlectra.CsvFormatRuntime(Sparlectra._resolve_detailed_csv_format("excel_de"))
      for value in (0.0, -0.0, 1.0, -1.0, 1.234567, -1234567.89, NaN, Inf, -Inf, Base.missing, nothing, "quoted; \"value\"\n")
        old_cell = Sparlectra._csv_field(value, ';', Sparlectra._resolve_detailed_csv_format("excel_de"))
        direct_buffer = IOBuffer()
        Sparlectra.write_csv_cell!(direct_buffer, value, ';', direct_cell_format)
        @test String(take!(direct_buffer)) == old_cell
      end
      @test_throws ArgumentError Sparlectra._resolve_detailed_csv_format("unknown")
      control_net = _control_label_test_net()
      control_cache = Sparlectra._bus_control_flag_cache(control_net)
      for node in control_net.nodeVec
        @test Sparlectra._cached_control_label(control_cache, node.busIdx) == Sparlectra._control_label(control_net, node.busIdx)
      end
      @test Set(Sparlectra._cached_control_label(control_cache, node.busIdx) for node in control_net.nodeVec) == Set(["-", "Q(U)", "P(U)", "Q(U), P(U)"])
      direct_source = read(joinpath(dirname(@__DIR__), "src", "api", "run_api.jl"), String)
      # Keep this as a simple source guard: the behavioral cache equivalence
      # checks above cover correctness, while these checks prevent accidental
      # reintroduction of the old per-row control-label path without brittle
      # regex body extraction.
      @test !occursin("_control_label(net", direct_source)
      @test !occursin("_bus_control_flags", direct_source)

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
      @test result.metadata["qlimits_enabled"] === true
      @test result.metadata["qlimit_guard_enabled"] isa Bool
      @test result.metadata["q_limit_preview_mode"] == "summary"
      @test result.metadata["q_limit_runlog_max_rows"] == 0
      @test read(template, String) == template_before
      @test isfile(joinpath(output_dir, "effective_config.yaml"))
      effective_config_text = read(joinpath(output_dir, "effective_config.yaml"), String)
      @test occursin("tol: 1.0e-9", effective_config_text)
      @test !occursin("runtime_request:", effective_config_text)
      effective_cfg = Sparlectra.load_sparlectra_config(joinpath(output_dir, "effective_config.yaml"); reload = true)
      @test effective_cfg.powerflow.tol == 1.0e-9
      @test effective_cfg.powerflow.max_iter == 40
      @test effective_cfg.matpower.case == "case14.m"
      @test isfile(joinpath(output_dir, "run_metadata.yaml"))
      run_metadata_text = read(joinpath(output_dir, "run_metadata.yaml"), String)
      @test occursin("runtime_request:", run_metadata_text)
      @test occursin("casefile: case_api.m", run_metadata_text)
      @test occursin("solver_status: completed", run_metadata_text)
      @test occursin("final_outcome:", run_metadata_text)
      @test occursin("numerical_converged: true", run_metadata_text)
      @test occursin("solution_available: true", run_metadata_text)
      run_log = read(joinpath(output_dir, "run.log"), String)
      @test occursin("Resolved Q-limit options", run_log)
      @test occursin("Q-limit handling enabled : true", run_log)
      @test occursin("Q-limit guard enabled    : $(effective_cfg.powerflow.qlimits.guard)", run_log)
      @test occursin("Q-limit preview mode     : summary", run_log)
      @test occursin("Q-limit runlog max rows  : 0", run_log)
      @test occursin("Q-limit detail artifact  : q_limit.log", run_log)
      @test occursin("full details        : q_limit.log", run_log)
      @test !occursin("full details     : q_limit_initial_limits.csv", run_log)
      @test occursin("Runtime casefile:", run_log)
      @test occursin("Original MATPOWER import options", run_log)
      @test occursin("Final effective MATPOWER import options", run_log)
      @test occursin("matpower_import.auto_profile: recommend", run_log)
      @test occursin("Wall time   :", run_log)
      @test occursin("Output time :", run_log)
      @test occursin("Solver time :", run_log)
      @test !occursin("representative_time:", run_log)
      @test !occursin("solver_time:          n/a", run_log)
      for marker in ("iterations:", "final_mismatch:", "final_status:", "final_outcome:")
        @test occursin(marker, run_log)
      end
      for phase in ("reading_matpower_case", "building_sparlectra_net", "solving_powerflow", "writing_artifacts", "finalizing_success")
        @test any(timing -> get(timing, "phase", "") == phase, result.service_phase_timings)
        @test occursin(phase, run_log)
      end
      @test all(timing -> get(timing, "elapsed_seconds", 0.0) === nothing || get(timing, "elapsed_seconds", 0.0) >= 0.0, result.service_phase_timings)
      @test occursin("Large case timing summary", run_log)
      result_json = Sparlectra._parse_service_json(read(joinpath(output_dir, "result.json"), String))
      @test haskey(result_json, "service_phase_timings")
      @test result_json["numerical_converged"] === true
      @test result_json["metadata"]["final_outcome"]["iterations"] == result.iterations
      @test result_json["metadata"]["final_outcome"]["final_mismatch"] == result.final_mismatch
      @test any(timing -> get(timing, "phase", "") == "finalizing_success", result_json["service_phase_timings"])
      @test occursin("detailed_result_csv_format: excel_de", run_log)
      @test occursin("detailed_result_csv_exporter: report", run_log)
      @test occursin("detailed_result_csv_delimiter: semicolon", run_log)
      @test occursin("detailed_result_csv_decimal_separator: comma", run_log)
      @test occursin("detailed_result_csv_thousands_separator: dot", run_log)
      performance_log = read(joinpath(output_dir, "performance.log"), String)
      @test occursin("Sparlectra single-run phase timing", performance_log)
      @test occursin("api_config_build:", performance_log)
      @test occursin("case_loading_network_solver:", performance_log)
      @test occursin("total:", performance_log)
      @test occursin("Phase timings", performance_log)
      @test occursin("reading_matpower_case", performance_log)
      diagnostic_log = read(joinpath(output_dir, "diagnose.log"), String)
      @test occursin("Sparlectra PowerFlow diagnostics", diagnostic_log)
      @test occursin("Final limit validation:", diagnostic_log)
      @test !isfile(joinpath(output_dir, "diagnose.txt"))
      q_limit_log = read(joinpath(output_dir, "q_limit.log"), String)
      @test occursin("Resolved Q-limit options", q_limit_log)
      @test occursin("Initial PV Q-limit table", q_limit_log)
      @test occursin("PV->PQ and PQ->PV event details", q_limit_log)
      @test occursin("Final PV Q-limit active-set check", q_limit_log)
      @test occursin("Q-Limit Active-Set Summary", q_limit_log)

      bus_csv = read(joinpath(output_dir, "bus_voltages_complex.csv"), String)
      branch_csv = read(joinpath(output_dir, "branch_flows.csv"), String)
      @test startswith(bus_csv, "bus;bus_name;type;vm_pu;va_deg;vn_kV;v_re;v_im;v_complex;v_kV")
      @test startswith(branch_csv, "branch;branch_index;from_bus;to_bus;status;p_from_MW;q_from_MVar;p_to_MW;q_to_MVar")
      @test occursin(r"\d,\d", bus_csv)
      @test length(collect(eachline(IOBuffer(bus_csv)))) > 1
      @test length(collect(eachline(IOBuffer(branch_csv)))) > 1
      power_cache = Sparlectra._bus_power_component_cache(result.raw_result.net)
      for node in result.raw_result.net.nodeVec
        @test Sparlectra._bus_power_components(power_cache, node.busIdx) == Sparlectra._effective_bus_power_components(result.raw_result.net, node.busIdx)
      end
      direct_csv_dir = joinpath(tmpdir, "direct_csv")
      mkpath(direct_csv_dir)
      direct_cfg = Sparlectra.OutputConfig(detailed_result_csv_exporter = :direct)
      direct_timing = Dict{Symbol,Any}()
      direct_artifacts = Sparlectra._write_detailed_result_csv(direct_csv_dir, result.raw_result; format = "excel_de", config = direct_cfg, timing_metadata = direct_timing)
      @test direct_artifacts == ["bus_voltages_complex.csv", "branch_flows.csv"]
      @test replace(read(joinpath(direct_csv_dir, "bus_voltages_complex.csv"), String), "\r\n" => "\n") == replace(bus_csv, "\r\n" => "\n")
      @test replace(read(joinpath(direct_csv_dir, "branch_flows.csv"), String), "\r\n" => "\n") == replace(branch_csv, "\r\n" => "\n")
      @test direct_timing[:exporter] === :direct
      @test direct_timing[:write_mode] === :streaming
      @test direct_timing[:bus_rows] == length(result.raw_result.net.nodeVec)
      @test direct_timing[:branch_rows] == length(result.raw_result.net.branchVec)
      @test direct_timing[:bus_export_s] >= 0.0
      @test direct_timing[:branch_export_s] >= 0.0
      @test direct_timing[:control_label_cache_s] >= 0.0
      @test Sparlectra._select_detailed_csv_exporter(result.raw_result.net; config = Sparlectra.OutputConfig(detailed_result_csv_exporter = :auto, detailed_result_csv_direct_threshold_buses = 1)) === :direct
      progress_events = String[]
      progress_payloads = Dict{String,Any}[]
      progress_dir = joinpath(tmpdir, "direct_csv_progress")
      mkpath(progress_dir)
      Sparlectra._write_detailed_result_csv(
        progress_dir,
        result.raw_result;
        format = "excel_de",
        config = direct_cfg,
        timing_metadata = Dict{Symbol,Any}(:progress_callback => ((event; fields...) -> (push!(progress_events, String(event)); push!(progress_payloads, Dict{String,Any}(String(k) => v for (k, v) in fields))))),
      )
      @test "detailed_result_csv_export_started" in progress_events
      @test count(==("detailed_result_csv_file_written"), progress_events) == 2
      @test any(p -> get(p, "artifact", "") == "bus_voltages_complex.csv" && get(p, "rows", 0) == length(result.raw_result.net.nodeVec) && haskey(p, "elapsed_s"), progress_payloads)
      @test any(p -> get(p, "artifact", "") == "branch_flows.csv" && get(p, "rows", 0) == length(result.raw_result.net.branchVec) && haskey(p, "elapsed_s"), progress_payloads)
      abort_checks = Ref(0)
      abort_dir = joinpath(tmpdir, "direct_csv_abort_checks")
      mkpath(abort_dir)
      Sparlectra._write_detailed_result_csv(abort_dir, result.raw_result; format = "technical", config = direct_cfg, abort_checker = () -> (abort_checks[] += 1))
      @test abort_checks[] >= 1
      @test length(power_cache) == length(result.raw_result.net.nodeVec)

      kinds = Set(artifact.kind for artifact in result.artifacts)
      @test :log in kinds
      @test :q_limit_log in kinds
      @test :result_json in kinds
      @test :effective_config in kinds
      @test :run_metadata in kinds
      @test count(artifact -> artifact.kind === :csv && artifact.mime_type == "text/csv", result.artifacts) == 2
      @test all(artifact -> artifact.exists && artifact.size_bytes !== nothing, result.artifacts)
      @test "q_limit.log" in Set(artifact.name for artifact in result.artifacts)
      for m in eachmatch(r"full details\s*:\s*(\S+)", run_log)
        artifact_name = m.captures[1]
        @test isfile(joinpath(output_dir, artifact_name))
        @test artifact_name in Set(artifact.name for artifact in result.artifacts)
      end
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
      ignored_csv_format = run_sparlectra_api(
        casefile = casefile,
        config_file = template,
        output_dir = joinpath(tmpdir, "ignored_csv_format"),
        config_overrides = Dict("benchmark.enabled" => false),
        detailed_result_csv = false,
        detailed_result_csv_format = "unknown",
      )
      @test ignored_csv_format.success
      invalid_csv_format = run_sparlectra_api(
        casefile = casefile,
        config_file = template,
        output_dir = joinpath(tmpdir, "invalid_csv_format"),
        config_overrides = Dict("benchmark.enabled" => false),
        detailed_result_csv = true,
        detailed_result_csv_format = "unknown",
      )
      @test !invalid_csv_format.success
      @test invalid_csv_format.reason == "invalid_detailed_result_csv_format"
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

      qlimit_disabled_dir = joinpath(tmpdir, "q_limit_disabled")
      mkpath(qlimit_disabled_dir)
      disabled_metadata = copy(result.metadata)
      disabled_metadata["qlimits_enabled"] = false
      Sparlectra._write_q_limit_log_artifact(qlimit_disabled_dir, result.raw_result, disabled_metadata)
      qlimit_disabled_log = read(joinpath(qlimit_disabled_dir, "q_limit.log"), String)
      @test occursin("Q-limit handling enabled : false", qlimit_disabled_log)
      @test occursin("Q-limit enforcement      : skipped", qlimit_disabled_log)
      @test occursin("Q-limit active set       : not run", qlimit_disabled_log)
      @test occursin("PV->PQ Q-limit events    : 0", qlimit_disabled_log)
      @test occursin("Guarded PV->PQ locks     : 0", qlimit_disabled_log)

      pre_mutation_net = _control_label_test_net()
      for bus in 2:4
        pre_mutation_net.nodeVec[bus]._nodeType = Sparlectra.PV
      end
      pre_mutation_net.qmin_pu = fill(-1.0, length(pre_mutation_net.nodeVec))
      pre_mutation_net.qmax_pu = fill(1.0, length(pre_mutation_net.nodeVec))
      Sparlectra.snapshotPVQLimits!(pre_mutation_net)
      for bus in 2:4
        pre_mutation_net.nodeVec[bus]._nodeType = Sparlectra.PQ
        Sparlectra.logQLimitHit!(pre_mutation_net, 1, bus, :max)
      end
      pre_mutation_io = IOBuffer()
      Sparlectra.printPVQLimitsTable(pre_mutation_net; io = pre_mutation_io, max_rows = 10)
      pre_mutation_text = String(take!(pre_mutation_io))
      @test occursin("rows total       : 3", pre_mutation_text)
      @test !occursin("no PV buses", pre_mutation_text)

      nonconverged_result = SparlectraRunResult(
        pre_mutation_net,
        :not_converged,
        false,
        false,
        :skip,
        false,
        :nr_mismatch_not_converged_active_set_unstable,
        "NR mismatch did not converge; Q-limit active set changed repeatedly",
        80,
        0.0,
        nothing,
        1.0e6,
        :rectangular,
        :not_run,
        nothing,
        (
          q_limit_active_set_ok = false,
          pv_q_limit_violations = 0,
          ref_q_limit_violations = 0,
          pv_pq_switching_events = length(pre_mutation_net.qLimitLog),
          qlimit_active_set_changes = 3,
          qlimit_reenable_events = 0,
          oscillating_buses = 0,
          guarded_narrow_q_pv_buses = 0,
        ),
      )
      nonconverged_dir = joinpath(tmpdir, "q_limit_nonconverged")
      mkpath(nonconverged_dir)
      Sparlectra._write_q_limit_log_artifact(nonconverged_dir, nonconverged_result, result.metadata)
      nonconverged_log = read(joinpath(nonconverged_dir, "q_limit.log"), String)
      @test occursin("Final Q-limit validation skipped: invalid because NR did not converge.", nonconverged_log)
      @test occursin("Final PV/REF Q-limit check: SKIPPED (NR did not converge)", nonconverged_log)

      # Non-converged solver states still carry useful voltages and branch data
      # when a solution/network state is available; detailed CSV export must not
      # be gated only by final convergence.
      nonconverged_solution = SparlectraRunResult(
        result.raw_result.net,
        :not_converged,
        false,
        false,
        :skip,
        false,
        :nr_mismatch_not_converged,
        "NR mismatch did not converge",
        result.raw_result.iterations,
        result.raw_result.elapsed_s,
        result.raw_result.solver_elapsed_s,
        result.raw_result.final_mismatch,
        result.raw_result.method,
        result.raw_result.control_status,
        result.raw_result.performance_profile,
        result.raw_result.diagnostics,
      )
      nonconverged_csv_dir = joinpath(tmpdir, "nonconverged_csv")
      mkpath(nonconverged_csv_dir)
      nonconverged_csv = Sparlectra._write_detailed_result_csv(nonconverged_csv_dir, nonconverged_solution; format = "technical", config = direct_cfg)
      @test nonconverged_csv == ["bus_voltages_complex.csv", "branch_flows.csv"]
      @test isfile(joinpath(nonconverged_csv_dir, "bus_voltages_complex.csv"))
      @test isfile(joinpath(nonconverged_csv_dir, "branch_flows.csv"))
      @test Sparlectra._csv_solution_quality(nonconverged_solution) == "not_converged_last_iterate"
      synthetic_nonconverged_api = Sparlectra._api_result(
        status = :not_converged,
        success = false,
        converged = false,
        solution_available = false,
        iterations = nonconverged_solution.iterations,
        final_mismatch = nonconverged_solution.final_mismatch,
        reason = String(nonconverged_solution.reason),
        message = "PowerFlow run completed, but numerical solver did not converge.",
        output_dir = nonconverged_csv_dir,
        metadata = Dict{String,Any}(
          "service_status" => "completed",
          "numerical_status" => "not_converged",
          "run_status" => "completed_nonconverged",
          "detailed_result_csv_status" => "exported_diagnostic",
          "detailed_result_csv_solution_quality" => "not_converged_last_iterate",
        ),
        raw_result = nonconverged_solution,
      )
      synthetic_nonconverged_dict = Sparlectra.to_dict(synthetic_nonconverged_api)
      @test synthetic_nonconverged_dict["status"] == "not_converged"
      @test synthetic_nonconverged_dict["success"] === false
      @test synthetic_nonconverged_dict["service_status"] == "completed"
      @test synthetic_nonconverged_dict["numerical_status"] == "not_converged"
      @test synthetic_nonconverged_dict["run_status"] == "completed_nonconverged"
      @test !occursin("Final PV/REF Q-limit check: OK", nonconverged_log)

      failed_diagnostic = joinpath(tmpdir, "failed_diagnose.txt")
      Sparlectra._write_powerflow_diagnostics(failed_diagnostic, result.raw_result; diagnostic_fn = (io, _) -> error("diagnostic test failure"))
      failed_diagnostic_text = read(failed_diagnostic, String)
      @test occursin("Diagnostic generation failed", failed_diagnostic_text)
      @test occursin("diagnostic test failure", failed_diagnostic_text)

      for i in 1:5
        Sparlectra.logQLimitHit!(control_net, i, 2 + (i % 3), i % 2 == 0 ? :max : :min)
      end
      qlimit_io = IOBuffer()
      Sparlectra.printQLimitLog(control_net; io = qlimit_io, max_rows = 2, full_details = "q_limit_events.csv")
      qlimit_text = String(take!(qlimit_io))
      @test occursin("PV->PQ switching events:", qlimit_text)
      @test occursin("total events     : 5", qlimit_text)
      @test occursin("rows shown       : 2", qlimit_text)
      @test occursin("rows omitted     : 3", qlimit_text)
      @test occursin("full details     : q_limit_events.csv", qlimit_text)
      qlimit_summary_io = IOBuffer()
      Sparlectra.printQLimitLog(control_net; io = qlimit_summary_io, max_rows = 0, full_details = "q_limit_events.csv")
      qlimit_summary_text = String(take!(qlimit_summary_io))
      @test occursin("rows shown       : 0", qlimit_summary_text)
      @test !occursin("Iteration │ Bus │ Side", qlimit_summary_text)
      pvlimit_net = _control_label_test_net()
      for bus in 2:4
        pvlimit_net.nodeVec[bus]._nodeType = Sparlectra.PV
      end
      pvlimit_net.qmin_pu = fill(-1.0, length(pvlimit_net.nodeVec))
      pvlimit_net.qmax_pu = fill(1.0, length(pvlimit_net.nodeVec))
      pvlimit_io = IOBuffer()
      Sparlectra.printPVQLimitsTable(pvlimit_net; io = pvlimit_io, max_rows = 2, full_details = "q_limit_initial_limits.csv")
      pvlimit_text = String(take!(pvlimit_io))
      @test occursin("rows shown       : 2", pvlimit_text)
      @test occursin("rows omitted", pvlimit_text)
      @test occursin("full details     : q_limit_initial_limits.csv", pvlimit_text)
      pvlimit_summary_io = IOBuffer()
      Sparlectra.printPVQLimitsTable(pvlimit_net; io = pvlimit_summary_io, max_rows = 0, full_details = "q_limit_initial_limits.csv")
      pvlimit_summary_text = String(take!(pvlimit_summary_io))
      @test occursin("rows shown       : 0", pvlimit_summary_text)
      @test occursin("preview          : summary-only", pvlimit_summary_text)
      @test !occursin("Bus │      Qmin", pvlimit_summary_text)
      qlimit_artifact_dir = joinpath(tmpdir, "q_limit_artifacts")
      mkpath(qlimit_artifact_dir)
      qlimit_artifacts = Sparlectra._write_q_limit_detail_artifacts(qlimit_artifact_dir, control_net)
      @test "q_limit_events.csv" in qlimit_artifacts
      @test length(collect(eachline(joinpath(qlimit_artifact_dir, "q_limit_events.csv")))) == 6

      dict_result = to_dict(result)
      named_result = to_namedtuple(result)
      json_result = to_json(result)
      yaml_result = to_yaml(result)
      @test dict_result["run_id"] == result.run_id
      @test dict_result["schema_version"] == "1.0"
      @test dict_result["status"] == "succeeded"
      @test dict_result["qlimits_enabled"] === true
      @test dict_result["q_limit_preview_mode"] == "summary"
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
      @test any(artifact -> artifact.kind === :run_metadata, missing.artifacts)

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

      import_overrides = validate_gui_config_overrides(Dict(
        "matpower_import.auto_profile" => "apply",
        "matpower_import.ratio" => "reciprocal",
        "matpower_import.shift_sign" => -1.0,
        "matpower_import.shift_unit" => "rad",
        "matpower_import.bus_shunt_model" => "voltage_dependent_injection",
      ))
      @test import_overrides["matpower_import"]["auto_profile"] == "apply"
      @test import_overrides["matpower_import"]["ratio"] == "reciprocal"
      @test import_overrides["matpower_import"]["shift_sign"] == -1.0
      @test import_overrides["matpower_import"]["shift_unit"] == "rad"
      @test import_overrides["matpower_import"]["bus_shunt_model"] == "voltage_dependent_injection"
      @test_throws ArgumentError validate_gui_config_overrides(Dict("matpower_import.auto_profile" => "nonsense"))
      @test_throws ArgumentError validate_gui_config_overrides(Dict("matpower_import.ratio" => "sideways"))
      @test_throws ArgumentError validate_gui_config_overrides(Dict("matpower_import.shift_sign" => 0.0))

      manual_import_output = joinpath(tmpdir, "manual_import")
      manual_import_result = run_sparlectra_api(
        casefile = casefile,
        config_file = template,
        output_dir = manual_import_output,
        config_overrides = Dict(
          "matpower_import.auto_profile" => "off",
          "matpower_import.ratio" => "reciprocal",
          "matpower_import.shift_sign" => -1.0,
          "matpower_import.shift_unit" => "rad",
          "benchmark.enabled" => false,
        ),
      )
      @test manual_import_result.success
      manual_effective_cfg = Sparlectra.load_sparlectra_config(joinpath(manual_import_output, "effective_config.yaml"); reload = true)
      @test manual_effective_cfg.matpower.auto_profile === :off
      @test manual_effective_cfg.matpower.ratio === :reciprocal
      @test manual_effective_cfg.matpower.shift_sign == -1.0
      @test manual_effective_cfg.matpower.shift_unit === :rad
      manual_run_log = read(joinpath(manual_import_output, "run.log"), String)
      @test occursin("Original MATPOWER import options", manual_run_log)
      @test occursin("Final effective MATPOWER import options", manual_run_log)
      @test occursin("matpower_import.auto_profile: off", manual_run_log)
      @test occursin("matpower_import.ratio: reciprocal", manual_run_log)
      @test occursin("matpower_import.shift_sign: -1.0", manual_run_log)
      @test occursin("matpower_import.shift_unit: rad", manual_run_log)

      auto_profile_output = joinpath(tmpdir, "auto_profile_recommend")
      auto_profile_result = run_sparlectra_api(
        casefile = casefile,
        config_file = template,
        output_dir = auto_profile_output,
        config_overrides = Dict(
          "matpower_import.auto_profile" => "recommend",
          "benchmark.enabled" => false,
        ),
      )
      @test auto_profile_result.success
      @test isfile(joinpath(auto_profile_output, "matpower_auto_profile.log"))
      @test any(artifact -> artifact.kind === :matpower_auto_profile, auto_profile_result.artifacts)
      auto_profile_log = read(joinpath(auto_profile_output, "matpower_auto_profile.log"), String)
      @test occursin("Runtime casefile:", auto_profile_log)
      @test occursin("Original MATPOWER import options:", auto_profile_log)
      @test occursin("Auto-profile recommendation:", auto_profile_log)
      @test occursin("Final effective MATPOWER import options:", auto_profile_log)
      @test occursin("matpower_import.ratio", auto_profile_log)
      @test occursin("Final effective MATPOWER auto-profile options", auto_profile_log)
      recommend_effective_cfg = Sparlectra.load_sparlectra_config(joinpath(auto_profile_output, "effective_config.yaml"); reload = true)
      @test recommend_effective_cfg.matpower.auto_profile === :recommend
      @test recommend_effective_cfg.matpower.ratio === :normal
      @test recommend_effective_cfg.matpower.shift_sign == 1.0
      @test recommend_effective_cfg.matpower.shift_unit === :deg

      auto_profile_apply_output = joinpath(tmpdir, "auto_profile_apply")
      auto_profile_apply_result = run_sparlectra_api(
        casefile = casefile,
        config_file = template,
        output_dir = auto_profile_apply_output,
        config_overrides = Dict(
          "matpower_import.auto_profile" => "apply",
          "benchmark.enabled" => false,
        ),
      )
      @test auto_profile_apply_result.success
      apply_effective_text = read(joinpath(auto_profile_apply_output, "effective_config.yaml"), String)
      apply_profile_log = read(joinpath(auto_profile_apply_output, "matpower_auto_profile.log"), String)
      @test occursin("auto_profile: apply", apply_effective_text)
      @test occursin("auto_profile = apply", apply_profile_log)
      @test any(artifact -> artifact.kind === :matpower_auto_profile, auto_profile_apply_result.artifacts)

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
        resolved_existing_m = Sparlectra._resolve_powerflow_casefile(existing_m, joinpath(tmpdir, "cases"))
        @test lowercase(splitext(resolved_existing_m)[2]) == ".m"
        @test isfile(resolved_existing_m)

        existing_jl = joinpath(tmpdir, "existing_case.jl")
        write(existing_jl, "nothing\n")
        @test Sparlectra._resolve_powerflow_casefile(existing_jl, joinpath(tmpdir, "cases")) == abspath(existing_jl)

        sibling_m = joinpath(tmpdir, "sibling_case.m")
        sibling_jl = joinpath(tmpdir, "sibling_case.jl")
        write(sibling_m, "not parsed when the Julia case exists\n")
        write(sibling_jl, "nothing\n")
        @test Sparlectra._resolve_powerflow_casefile(sibling_m, joinpath(tmpdir, "cases")) == abspath(sibling_m)
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
        resolved_missing = Sparlectra._resolve_powerflow_casefile("case118.m", case_directory; ensure_casefile_fn = fake_ensure)
        @test resolved_missing == abspath(joinpath(case_directory, "case118.m"))
        @test only(ensure_calls) == (requested = "case118.m", outdir = abspath(case_directory), to_jl = false)
        empty!(ensure_calls)
        write(joinpath(case_directory, "case14.m"), "downloaded MATPOWER placeholder\n")
        write(joinpath(case_directory, "case14.jl"), "nothing\n")
        resolved_requested_jl = Sparlectra._resolve_powerflow_casefile("case14.jl", case_directory; ensure_casefile_fn = fake_ensure)
        @test resolved_requested_jl == abspath(joinpath(case_directory, "case14.m"))
        @test isempty(ensure_calls)

        only_jl_dir = joinpath(tmpdir, "only_jl_cases")
        mkpath(only_jl_dir)
        write(joinpath(only_jl_dir, "case118.jl"), "nothing\n")
        err = @test_throws ArgumentError Sparlectra._resolve_powerflow_casefile("case118.jl", only_jl_dir; ensure_casefile_fn = fake_ensure)
        @test occursin("Generated MATPOWER .jl cache files are not user-selectable", sprint(showerror, err.value))

        empty!(ensure_calls)
        @test_throws ArgumentError Sparlectra._resolve_powerflow_casefile(joinpath("missing", "case14.m"), case_directory; ensure_casefile_fn = fake_ensure)
        @test_throws ArgumentError Sparlectra._resolve_powerflow_casefile("https://example.com/case14.m", case_directory; ensure_casefile_fn = fake_ensure)
        @test isempty(ensure_calls)

        rejected_path = start_powerflow_run(Dict("casefile" => joinpath("missing", "case14.m")); case_directory)
        @test rejected_path["reason"] == "invalid_casefile"
        rejected_url = start_powerflow_run(Dict("casefile" => "https://example.com/case14.m"); case_directory)
        @test rejected_url["reason"] == "invalid_casefile"
      end

      @testset "Broken generated Julia MATPOWER cache fails cleanly" begin
        broken_dir = joinpath(tmpdir, "broken_cache")
        output_broken = joinpath(tmpdir, "broken_output")
        mkpath(broken_dir)
        broken_case = joinpath(broken_dir, "case_broken.jl")
        write(broken_case, "throw(StackOverflowError())\n")
        broken = start_powerflow_run(Dict(
          "casefile" => broken_case,
          "config_file" => config_file,
          "output_root" => output_broken,
          "performance_timing" => "compact",
          "config_overrides" => Dict("benchmark.enabled" => false),
        ))
        @test broken["success"] === false
        @test broken["status"] == "failed"
        @test broken["reason"] == "loading_julia_case_failed"
        @test any(t -> get(t, "phase", "") == "loading_julia_case" && get(t, "status", "") == "failed", broken["service_phase_timings"])
        @test isfile(joinpath(broken["output_dir"], "run.log"))
        @test isfile(joinpath(broken["output_dir"], "result.json"))
        @test isfile(joinpath(broken["output_dir"], "performance.log"))
        @test occursin("StackOverflowError", read(joinpath(broken["output_dir"], "run.log"), String))
        @test occursin("StackOverflowError", read(joinpath(broken["output_dir"], "result.json"), String))
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
      @test haskey(aborted, "abort_phase")
      @test occursin("Current phase:", aborted["message"])
      @test !aborted["success"]
      @test Sparlectra.abort_webui_powerflow_run(active["run_id"])["abort_status"] == "already_aborting"
      @test Sparlectra.start_webui_powerflow_run(async_request; runner = controlled_runner)["reason"] == "active_run"
      @test Sparlectra.abort_webui_powerflow_run("../unsafe")["reason"] == "unsafe_run_id"
      wait(Sparlectra._POWERFLOW_WEBUI_JOBS[active["run_id"]]["task"])
      @test get_powerflow_result(active["run_id"])["status"] == "aborted"
      @test !get_powerflow_result(active["run_id"])["success"]
      @test occursin("Run aborted by user.", read(joinpath(aborted["output_dir"], "run.log"), String))
      @test occursin("Phase active at abort:", read(joinpath(aborted["output_dir"], "run.log"), String))
      @test Sparlectra.abort_webui_powerflow_run(active["run_id"])["abort_status"] == "already_aborted"
      replacement = Sparlectra.start_webui_powerflow_run(async_request; runner = controlled_runner)
      @test replacement["status"] in ("queued", "running")
      put!(release, nothing)
      wait(Sparlectra._POWERFLOW_WEBUI_JOBS[replacement["run_id"]]["task"])
      @test Sparlectra.get_webui_powerflow_job(active["run_id"])["status"] == "aborted"
      @test Sparlectra.get_active_webui_powerflow_job() === nothing

      artifact_error_runner = function(request; case_directory = nothing)
        request["phase_callback"]("solving_powerflow")
        request["phase_callback"]("postprocessing_result")
        request["phase_callback"]("writing_artifacts")
        error("injected artifact failure")
      end
      artifact_error = Sparlectra.start_webui_powerflow_run(async_request; runner = artifact_error_runner)
      wait(Sparlectra._POWERFLOW_WEBUI_JOBS[artifact_error["run_id"]]["task"])
      artifact_error_result = get_powerflow_result(artifact_error["run_id"])
      @test artifact_error_result["status"] == "failed"
      @test artifact_error_result["run_status"] == "failed"
      @test artifact_error_result["solver_status"] == "completed"
      @test artifact_error_result["artifact_status"] == "running"
      @test artifact_error_result["last_phase"] == "finalizing_failed"

      artifact_abort_entered = Channel{Nothing}(1)
      artifact_abort_release = Channel{Nothing}(1)
      artifact_abort_runner = function(request; case_directory = nothing)
        request["phase_callback"]("solving_powerflow")
        request["phase_callback"]("postprocessing_result")
        request["phase_callback"]("writing_artifacts")
        put!(artifact_abort_entered, nothing)
        token = request["cancellation_token"]
        while !isready(artifact_abort_release)
          Sparlectra._check_powerflow_cancelled!(token)
          yield()
        end
        take!(artifact_abort_release)
        return Dict{String,Any}("status" => "succeeded", "success" => true, "run_id" => request["run_id"], "artifacts" => Any[])
      end
      artifact_abort = Sparlectra.start_webui_powerflow_run(async_request; runner = artifact_abort_runner)
      take!(artifact_abort_entered)
      @test Sparlectra.abort_webui_powerflow_run(artifact_abort["run_id"])["status"] == "aborting"
      wait(Sparlectra._POWERFLOW_WEBUI_JOBS[artifact_abort["run_id"]]["task"])
      artifact_abort_result = get_powerflow_result(artifact_abort["run_id"])
      @test artifact_abort_result["status"] == "aborted"
      @test artifact_abort_result["run_status"] == "aborted"
      @test artifact_abort_result["last_phase"] == "finalizing_aborted"

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
      @test Set(["result.json", "run.log", "effective_config.yaml", "run_metadata.yaml", "q_limit.log"]) ⊆ artifact_names

      resolved = resolve_powerflow_artifact(started["run_id"], "result.json")
      @test resolved isa SparlectraApiArtifact
      @test isfile(resolved.path)
      @test !isempty(resolved.mime_type)
      qlimit_resolved = resolve_powerflow_artifact(started["run_id"], "q_limit.log")
      @test qlimit_resolved isa SparlectraApiArtifact
      @test qlimit_resolved.kind === :q_limit_log
      @test qlimit_resolved.mime_type == "text/plain"

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
      @test Set(["result.json", "run.log", "effective_config.yaml", "run_metadata.yaml"]) ⊆ recovered_names
      recovered_result = resolve_powerflow_artifact(started["run_id"], "result.json")
      @test recovered_result isa SparlectraApiArtifact
      @test recovered_result.path == joinpath(started["output_dir"], "result.json")

      stale_csv_id = "stale-csv-artifacts"
      stale_csv_dir = joinpath(output_root, stale_csv_id)
      mkpath(stale_csv_dir)
      stale_csv_result_file = joinpath(stale_csv_dir, "result.json")
      write(stale_csv_result_file, "{}\n")
      write(joinpath(stale_csv_dir, "bus_voltages_complex.csv"), "bus,vm_pu\n1,1.0\n")
      stale_csv_result = Sparlectra._api_result(
        run_id = stale_csv_id,
        status = :succeeded,
        success = true,
        output_dir = stale_csv_dir,
        result_file = stale_csv_result_file,
        artifacts = SparlectraApiArtifact[],
      )
      Sparlectra._register_powerflow_run!(stale_csv_result)
      stale_csv_artifacts = list_powerflow_artifacts(stale_csv_id)
      @test "bus_voltages_complex.csv" in Set(artifact["name"] for artifact in stale_csv_artifacts)
      stale_csv_resolved = resolve_powerflow_artifact(stale_csv_id, "bus_voltages_complex.csv")
      @test stale_csv_resolved isa SparlectraApiArtifact
      @test stale_csv_resolved.kind === :csv
      @test read(stale_csv_resolved.path, String) == "bus,vm_pu\n1,1.0\n"

      stale_root = joinpath(tmpdir, "stale-runs")
      stale_run_id = "stale-run"
      stale_job = Dict{String,Any}(
        "run_id" => stale_run_id,
        "status" => "running",
        "casefile" => casefile,
        "config_file" => config_file,
        "output_root" => stale_root,
        "output_dir" => joinpath(stale_root, stale_run_id),
        "started_at" => Dates.now(Dates.UTC),
        "finished_at" => nothing,
        "message" => "PowerFlow run is active.",
        "abort_requested" => Threads.Atomic{Bool}(false),
        "current_phase" => "writing_artifacts",
        "last_phase" => "writing_artifacts",
        "solver_status" => "completed",
        "artifact_status" => "running",
        "run_status" => "finalizing",
        "last_heartbeat" => "2026-06-18T12:17:07.854Z",
      )
      Sparlectra._write_webui_job_marker!(stale_job, :running, "webui_job_active", "PowerFlow run is active.")
      first_stale_refresh = refresh_powerflow_run_registry!(stale_root)
      @test [result.run_id for result in first_stale_refresh["stale_recovered_runs"]] == [stale_run_id]
      stale_result = get_powerflow_result(stale_run_id)
      @test stale_result["status"] == "interrupted_unknown"
      @test stale_result["final_outcome"] == "webui_restart_lost_live_task"
      @test stale_result["last_phase"] == "writing_artifacts"
      @test stale_result["solver_status"] == "completed"
      @test stale_result["artifact_status"] == "running"
      second_stale_refresh = refresh_powerflow_run_registry!(stale_root)
      @test isempty(second_stale_refresh["stale_recovered_runs"])
      stale_status_html = Sparlectra.render_powerflow_result(stale_result)
      @test occursin("Run state was recovered after Web UI restart.", stale_status_html)
      @test occursin("Last known phase: <code>writing_artifacts</code>", stale_status_html)
      @test occursin("No live solver task is attached anymore.", stale_status_html)

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
