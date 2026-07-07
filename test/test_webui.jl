using Sparlectra
using Dates
using Sockets
using Test

const WEBUI_TRANSCRIPT_MARKERS = (
  "git apply",
  "diff --git",
  "Search evidence",
  "Unfinished items",
  "Created PR metadata",
  "make_pr",
)

const WEBUI_INTERNAL_LIFECYCLE_MARKERS = (
  "webui_start_requested",
  "webui_config_loaded",
  "webui_routes_registered",
  "webui_server_bound",
  "browser_close_monitor_skipped",
)

_assert_no_webui_transcript_markers(text::AbstractString) = foreach(marker -> @test(!occursin(marker, text)), WEBUI_TRANSCRIPT_MARKERS)
_assert_no_webui_internal_lifecycle_markers(text::AbstractString) = foreach(marker -> @test(!occursin(marker, text)), WEBUI_INTERNAL_LIFECYCLE_MARKERS)

function _webui_input_tag(html::AbstractString, name::AbstractString)
  match = Base.match(Regex("<input[^>]*name=\\\"$(name)\\\"[^>]*type=\\\"checkbox\\\"[^>]*>"), html)
  match === nothing && (match = Base.match(Regex("<input[^>]*name=\\\"$(name)\\\"[^>]*>"), html))
  @test match !== nothing
  return match === nothing ? "" : match.match
end

function _webui_select_block(html::AbstractString, name::AbstractString)
  match = Base.match(Regex("<select[^>]*name=\\\"$(name)\\\"[^>]*>.*?</select>"), html)
  @test match !== nothing
  return match === nothing ? "" : match.match
end

function _webui_assert_checked(html::AbstractString, name::AbstractString, expected::Bool)
  tag = _webui_input_tag(html, name)
  @test occursin("type=\"checkbox\"", tag)
  @test occursin(" checked", tag) == expected
end

function _webui_assert_value(html::AbstractString, name::AbstractString, expected::AbstractString)
  @test occursin("name=\"$(name)\"", html)
  @test occursin("value=\"$(expected)\"", _webui_input_tag(html, name))
end

function _webui_assert_selected(html::AbstractString, name::AbstractString, expected::AbstractString)
  block = _webui_select_block(html, name)
  @test occursin("<option value=\"$(expected)\" selected>", block)
end

function _write_webui_test_case(path::AbstractString)
  write(path, """
function mpc = case_webui
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

function _webui_test_form(casefile, config_file, output_root)
  return Dict(
    "casefile" => "",
    "casefile_manual" => casefile,
    "config_file" => config_file,
    "output_root" => output_root,
    "apply_webui_runtime_overrides" => "true",
    "power_flow_tol" => "1e-8",
    "power_flow_max_iter" => "80",
    "power_flow_autodamp" => "on",
    "power_flow_autodamp_min" => "0.05",
    "power_flow_qlimits_enabled" => "on",
    "power_flow_qlimits_enforcement_mode" => "classic_simultaneous",
    "power_flow_wrong_branch_detection" => "warn",
    "power_flow_start_angle_mode" => "dc",
    "power_flow_start_voltage_mode" => "profile_blend",
    "power_flow_start_current_iteration_enabled" => "true",
    "power_flow_start_current_iteration_max_iter" => "7",
    "power_flow_start_current_iteration_tol" => "1e-4",
    "power_flow_start_current_iteration_damping" => "0.4",
    "power_flow_start_current_iteration_accept_only_if_improved" => "true",
    "power_flow_start_current_iteration_min_improvement_factor" => "0.95",
    "power_flow_start_current_iteration_vm_min_pu" => "0.6",
    "power_flow_start_current_iteration_vm_max_pu" => "1.4",
    "power_flow_start_current_iteration_max_angle_step_deg" => "20.0",
    "power_flow_start_current_iteration_only_for_large_cases" => "false",
    "matpower_import_auto_profile" => "recommend",
    "matpower_import_ratio" => "normal",
    "matpower_import_shift_sign" => "1.0",
    "matpower_import_shift_unit" => "deg",
    "matpower_import_bus_shunt_model" => "admittance",
    "matpower_import_pv_voltage_source" => "gen_vg",
    "matpower_import_compare_voltage_reference" => "imported_setpoint",
    "output_logfile_results" => "compact",
    "performance_timing" => "compact",
    "run_diagnostics" => "on",
    "detailed_result_csv" => "on",
    "detailed_result_csv_format" => "excel_de",
    "benchmark_enabled" => "false",
    "benchmark_samples" => "10",
    "benchmark_seconds" => "1.0",
  )
end

function _webui_http_request(port::Integer, method::AbstractString, target::AbstractString; body::AbstractString = "")::String
  socket = connect(ip"127.0.0.1", UInt16(port))
  headers = isempty(body) ? "" : "Content-Type: application/x-www-form-urlencoded\r\nContent-Length: $(ncodeunits(body))\r\n"
  write(socket, "$(method) $(target) HTTP/1.1\r\nHost: 127.0.0.1\r\n$(headers)Connection: close\r\n\r\n$(body)")
  response = read(socket, String)
  close(socket)
  return response
end

function run_webui_tests()
  @testset "Local PowerFlow Web UI" begin
    @test isdefined(Sparlectra, :start_sparlectra_webui)
    @test isdefined(Sparlectra, :default_webui_output_root)
    @test_throws ArgumentError start_sparlectra_webui(host = "0.0.0.0")
    default_root = default_webui_output_root()
    @test isabspath(default_root)
    @test !startswith(normpath(default_root), normpath(pkgdir(Sparlectra)))
    if Sys.islinux()
      withenv("XDG_STATE_HOME" => joinpath(mktempdir(), "state")) do
        @test default_webui_output_root() == joinpath(ENV["XDG_STATE_HOME"], "sparlectra", "webui", "runs")
      end
    end
    launcher_source = read(joinpath(@__DIR__, "..", "start_webui.jl"), String)
    @test occursin("Sparlectra.start_sparlectra_webui", launcher_source)
    @test count("start_sparlectra_webui", launcher_source) == 1
    repository_text = join(read(path, String) for path in (
      joinpath(@__DIR__, "..", "README.md"),
      joinpath(@__DIR__, "..", "docs", "src", "webui.md"),
      joinpath(@__DIR__, "..", "docs", "src", "examples_overview.md"),
    ))
    @test !occursin("exp_webui_powerflow", repository_text)
    @test Sparlectra._format_elapsed_duration(0.4) == "00:00:00.400"
    @test Sparlectra._format_elapsed_duration(3) == "00:00:03.000"
    @test Sparlectra._format_elapsed_duration(83) == "00:01:23.000"
    @test Sparlectra._format_elapsed_duration(3725.125) == "01:02:05.125"
    @test Sparlectra._format_elapsed_duration(nothing) == "—"


    @testset "Configuration refresh controls" begin
      root = mktempdir()
      config_path = joinpath(root, "configuration.yaml")
      write(config_path, "power_flow:\n  start_mode:\n    voltage_mode: bus_vm_va_blend\n  qlimits:\n    enabled: true\n")
      runtime = Sparlectra._SparlectraWebUIRuntime(nothing, root, config_path, Sparlectra.webui_operation_log_path(root), nothing, Sparlectra.start_powerflow_run, false, false, time(), 0, nothing, IOBuffer(), ReentrantLock())
      before_page_text = read(config_path, String)
      page = String(Sparlectra.route_sparlectra_webui("GET", "/powerflow"; output_root = root, runtime).body)
      @test occursin("Check configuration", page)
      @test occursin("Refresh configuration", page)
      @test occursin("Configuration notice:", page)
      @test occursin("#configuration-maintenance", page)
      @test findfirst("Configuration maintenance", page) > findfirst("Advanced / expert options", page)
      @test findfirst("Configuration maintenance", page) > findfirst("Existing case file", page)
      @test read(config_path, String) == before_page_text

      check = Sparlectra.route_sparlectra_webui("POST", "/powerflow/config/check", Dict("config_file" => config_path); output_root = root, runtime)
      check_body = String(check.body)
      @test check.status == 200
      @test occursin("power_flow.qlimits.enforcement_mode", check_body)
      @test occursin("power_flow.start_mode.voltage_mode", check_body)
      @test !occursin(".bak-", check_body)
      @test read(config_path, String) == before_page_text

      refresh = Sparlectra.route_sparlectra_webui("POST", "/powerflow/config/refresh", Dict("config_file" => config_path); output_root = root, runtime)
      refresh_body = String(refresh.body)
      @test refresh.status == 200
      @test occursin("Configuration refreshed and written", refresh_body)
      @test occursin("Restart or reload the Web UI", refresh_body)
      @test occursin(config_path, refresh_body)
      @test !isempty(filter(name -> occursin(r"configuration.yaml\.bak-", name), readdir(root)))
      @test occursin("profile_blend", read(config_path, String))

      duplicate_path = joinpath(root, "duplicate.yaml")
      write(duplicate_path, "power_flow:\n  tol: 1.0e-8\n  tol: 1.0e-7\n")
      duplicate_before = read(duplicate_path, String)
      duplicate = Sparlectra.route_sparlectra_webui("POST", "/powerflow/config/refresh", Dict("config_file" => duplicate_path); output_root = root, runtime)
      @test duplicate.status == 400
      @test occursin("Duplicate YAML key detected", String(duplicate.body))
      @test read(duplicate_path, String) == duplicate_before
      @test isempty(filter(name -> occursin(r"duplicate.yaml\.bak-", name), readdir(root)))

      log = read(runtime.operation_log, String)
      @test occursin("config_refresh_check_started", log)
      @test occursin("config_refresh_check_completed", log)
      @test occursin("config_refresh_write_started", log)
      @test occursin("config_refresh_write_completed", log)
    end

    @testset "Case-specific settings profiles" begin
      root = mktempdir()
      write(joinpath(root, "case145.m"), "% case fixture\n")
      run_id = "case-settings-test"
      metadata = Dict{String,Any}(
        "runtime_casefile" => "case145.m",
        "runtime_casefile_path" => joinpath(root, "case145.m"),
        "webui_request_settings" => Dict{String,Any}(
          "power_flow.tol" => 1.0e-7,
          "power_flow.max_iter" => 42,
          "power_flow.autodamp" => true,
          "power_flow.autodamp_min" => 0.07,
          "power_flow.qlimits.enabled" => true,
          "power_flow.qlimits.enforcement_mode" => :active_set,
          "power_flow.wrong_branch_detection" => :warn,
          "power_flow.start_mode.angle_mode" => :dc,
          "power_flow.start_mode.voltage_mode" => "profile_blend",
          "matpower_import.auto_profile" => "recommend",
          "matpower_import.ratio" => "normal",
          "matpower_import.shift_sign" => 1.0,
          "matpower_import.shift_unit" => "deg",
          "matpower_import.bus_shunt_model" => "admittance",
          "matpower_import.pv_voltage_source" => "gen_vg",
          "matpower_import.compare_voltage_reference" => "imported_setpoint",
          "output.logfile_results" => "compact",
          "benchmark.enabled" => false,
          "benchmark.samples" => 10,
          "benchmark.seconds" => 1.0,
          "performance_timing" => "compact",
          "run_diagnostics" => false,
          "detailed_result_csv" => true,
          "detailed_result_csv_format" => "excel_de",
        ),
      )
      successful = Sparlectra._api_result(
        run_id = run_id,
        status = :succeeded,
        success = true,
        converged = true,
        solution_available = true,
        iterations = 3,
        final_mismatch = 1.0e-9,
        reason = "converged",
        message = "ok",
        casefile = joinpath(root, "case145.m"),
        config_file = "configuration.yaml",
        output_dir = joinpath(root, run_id),
        metadata = metadata,
      )
      Sparlectra._POWERFLOW_SERVICE_RUNS[run_id] = successful
      result_html = Sparlectra.render_powerflow_result(Sparlectra.to_dict(successful))
      @test occursin("Save settings for this case", result_html)
      @test !occursin("Save these settings anyway", result_html)

      save_response = Sparlectra.handle_powerflow_case_settings_save(run_id, Dict{String,String}(); output_root = root, operation_log = root)
      @test save_response.status == 200
      save_html = String(save_response.body)
      profile_path = Sparlectra._webui_case_settings_path(root, joinpath(root, "case145.m"))
      @test isfile(profile_path)
      @test occursin("/powerflow?casefile=$(Sparlectra._webui_urlencode(joinpath(root, "case145.m")))", save_html)
      @test !occursin("/powerflow?casefile=$(Sparlectra._webui_urlencode(profile_path))", save_html)
      profile_text = read(profile_path, String)
      @test occursin("profile_kind: webui_case_settings", profile_text)
      @test occursin("power_flow_tol: 1.0e-7", profile_text)
      @test occursin("power_flow_autodamp: true", profile_text)
      @test occursin("benchmark_enabled: false", profile_text)
      @test occursin("power_flow_qlimits_enforcement_mode: active_set", profile_text)
      @test occursin("detailed_result_csv_format: excel_de", profile_text)
      @test !occursin("effective_config", profile_text)
      @test basename(profile_path) == "case145.sparlectra-webui.yaml"
      profile = Sparlectra.load_yaml_dict(profile_path)
      @test profile["settings"]["power_flow_autodamp"] === true
      @test profile["settings"]["benchmark_enabled"] === false
      @test profile["settings"]["power_flow_max_iter"] == 42
      @test profile["settings"]["benchmark_samples"] == 10
      @test profile["settings"]["benchmark_seconds"] == 1.0
      @test profile["settings"]["power_flow_qlimits_enforcement_mode"] == "active_set"

      loaded_form = String(Sparlectra.route_sparlectra_webui("GET", "/powerflow?casefile=$(Sparlectra._webui_urlencode(joinpath(root, "case145.m")))"; output_root = root).body)
      @test occursin("Case-specific settings loaded from", loaded_form)
      @test occursin("case145.sparlectra-webui.yaml", loaded_form)
      _webui_assert_value(loaded_form, "power_flow_tol", "1.0e-7")
      _webui_assert_checked(loaded_form, "power_flow_autodamp", true)
      _webui_assert_checked(loaded_form, "benchmark_enabled", false)
      _webui_assert_selected(loaded_form, "power_flow_qlimits_enforcement_mode", "active_set")
      _webui_assert_selected(loaded_form, "detailed_result_csv_format", "excel_de")

      case118 = joinpath(root, "case118.m")
      write(case118, "% case fixture\n")
      write(Sparlectra._webui_case_settings_path(root, case118), """
profile_kind: webui_case_settings
schema_version: 1
settings:
  benchmark_enabled: false
  benchmark_samples: 10
  benchmark_seconds: 1.0
  detailed_result_csv: false
  detailed_result_csv_format: excel_de
  matpower_import_auto_profile: apply
  matpower_import_bus_shunt_model: admittance
  matpower_import_compare_voltage_reference: imported_setpoint
  matpower_import_pv_voltage_source: gen_vg
  matpower_import_ratio: normal
  matpower_import_shift_sign: 1.0
  matpower_import_shift_unit: deg
  output_logfile_results: compact
  performance_timing: compact
  power_flow_autodamp: false
  power_flow_autodamp_min: 0.05
  power_flow_max_iter: 80
  power_flow_qlimits_enabled: false
  power_flow_qlimits_enforcement_mode: active_set
  power_flow_start_angle_mode: classic
  power_flow_start_voltage_mode: classic
  power_flow_tol: 1.0e-8
  power_flow_wrong_branch_detection: off
  run_diagnostics: false
""")
      case118_form = String(Sparlectra.route_sparlectra_webui("GET", "/powerflow?casefile=$(Sparlectra._webui_urlencode(case118))"; output_root = root).body)
      @test occursin("Case-specific settings loaded from", case118_form)
      @test occursin("case118.sparlectra-webui.yaml", case118_form)
      _webui_assert_checked(case118_form, "power_flow_autodamp", false)
      _webui_assert_checked(case118_form, "power_flow_qlimits_enabled", false)
      _webui_assert_checked(case118_form, "benchmark_enabled", false)
      _webui_assert_checked(case118_form, "run_diagnostics", false)
      _webui_assert_checked(case118_form, "detailed_result_csv", false)
      _webui_assert_value(case118_form, "power_flow_tol", "1.0e-8")
      _webui_assert_value(case118_form, "power_flow_max_iter", "80")
      _webui_assert_value(case118_form, "power_flow_autodamp_min", "0.05")
      _webui_assert_value(case118_form, "benchmark_samples", "10")
      _webui_assert_value(case118_form, "benchmark_seconds", "1.0")
      _webui_assert_value(case118_form, "matpower_import_shift_sign", "1.0")
      _webui_assert_selected(case118_form, "power_flow_start_angle_mode", "classic")
      _webui_assert_selected(case118_form, "power_flow_start_voltage_mode", "classic")
      _webui_assert_selected(case118_form, "power_flow_wrong_branch_detection", "off")
      _webui_assert_selected(case118_form, "matpower_import_auto_profile", "apply")
      _webui_assert_selected(case118_form, "matpower_import_ratio", "normal")
      _webui_assert_selected(case118_form, "matpower_import_shift_unit", "deg")
      _webui_assert_selected(case118_form, "matpower_import_bus_shunt_model", "admittance")
      _webui_assert_selected(case118_form, "matpower_import_pv_voltage_source", "gen_vg")
      _webui_assert_selected(case118_form, "matpower_import_compare_voltage_reference", "imported_setpoint")
      _webui_assert_selected(case118_form, "output_logfile_results", "compact")
      _webui_assert_selected(case118_form, "performance_timing", "compact")
      _webui_assert_selected(case118_form, "detailed_result_csv_format", "excel_de")

      dropdown_form = String(Sparlectra.route_sparlectra_webui("GET", "/powerflow"; output_root = root, runtime = (; case_directory = root, config_file = "configuration.yaml", operation_log = Sparlectra.webui_operation_log_path(root), startup_config_error = nothing, runner = Sparlectra.start_powerflow_run)).body)
      @test occursin("case118.m ★", dropdown_form)
      @test occursin("<option value=\"case118.m\">case118.m ★</option>", dropdown_form)
      @test !occursin("value=\"case118.m ★\"", dropdown_form)
      @test occursin("data-case-settings-reload=\"true\"", dropdown_form)
      @test occursin("target.searchParams.set('casefile', caseSelect.value)", dropdown_form)
      @test occursin("target.searchParams.set('config_file', configInput.value)", dropdown_form)
      dropdown_loaded_form = String(Sparlectra.route_sparlectra_webui("GET", "/powerflow?casefile=$(Sparlectra._webui_urlencode(case118))&config_file=$(Sparlectra._webui_urlencode("configuration.yaml"))"; output_root = root).body)
      _webui_assert_value(dropdown_loaded_form, "power_flow_max_iter", "80")
      _webui_assert_selected(dropdown_loaded_form, "detailed_result_csv_format", "excel_de")

      case14 = joinpath(root, "case14.m")
      write(case14, "% case fixture\n")
      case14_form = String(Sparlectra.route_sparlectra_webui("GET", "/powerflow?casefile=$(Sparlectra._webui_urlencode(case14))"; output_root = root).body)
      @test !occursin("Case-specific settings loaded from", case14_form)
      _webui_assert_checked(case14_form, "power_flow_autodamp", true)
      _webui_assert_checked(case14_form, "power_flow_qlimits_enabled", true)
      _webui_assert_checked(case14_form, "benchmark_enabled", true)
      _webui_assert_checked(case14_form, "run_diagnostics", false)
      _webui_assert_value(case14_form, "power_flow_tol", "1.0e-5")
      _webui_assert_value(case14_form, "power_flow_max_iter", "80")
      _webui_assert_selected(case14_form, "power_flow_start_angle_mode", "dc")
      _webui_assert_selected(case14_form, "power_flow_start_voltage_mode", "profile_blend")
      start_voltage_select = _webui_select_block(case14_form, "power_flow_start_voltage_mode")
      @test occursin("value=\"profile_blend\"", start_voltage_select)
      @test !occursin("bus_vm_va_blend", start_voltage_select)
      _webui_assert_selected(case14_form, "matpower_import_auto_profile", "recommend")
      _webui_assert_selected(case14_form, "output_logfile_results", "full")

      request_form = _webui_test_form("case145.m", "configuration.yaml", root)
      request_form["power_flow_tol"] = "2e-6"
      request = Sparlectra.powerflow_webui_request(request_form; default_output_root = root)
      @test request["config_overrides"]["power_flow.tol"] == 2.0e-6

      yaml_first_form = copy(request_form)
      delete!(yaml_first_form, "apply_webui_runtime_overrides")
      yaml_first_request = Sparlectra.powerflow_webui_request(yaml_first_form; default_output_root = root)
      @test isempty(yaml_first_request["config_overrides"])
      @test yaml_first_request["config_override_source"] == "user_yaml"

      invalid_case = joinpath(root, "invalid.m")
      write(invalid_case, "% case fixture\n")
      write(Sparlectra._webui_case_settings_path(root, invalid_case), "not: [valid\n")
      invalid_form = String(Sparlectra.route_sparlectra_webui("GET", "/powerflow?casefile=$(Sparlectra._webui_urlencode(invalid_case))"; output_root = root).body)
      @test occursin("PowerFlow run", invalid_form)
      @test !occursin("Case-specific settings loaded from", invalid_form)
      @test occursin("case_settings_load_failed", read(Sparlectra.webui_operation_log_path(root), String))

      unsupported_case = joinpath(root, "unsupported_field.m")
      write(unsupported_case, "% case fixture\n")
      write(Sparlectra._webui_case_settings_path(root, unsupported_case), """
profile_kind: webui_case_settings
schema_version: 1
settings:
  power_flow_autodamp: false
  power_flow_start_angle_mode: [bad]
  detailed_result_csv_format: impossible
  power_flow_tol: 9.0e-7
""")
      unsupported_form = String(Sparlectra.route_sparlectra_webui("GET", "/powerflow?casefile=$(Sparlectra._webui_urlencode(unsupported_case))"; output_root = root).body)
      @test occursin("PowerFlow run", unsupported_form)
      _webui_assert_checked(unsupported_form, "power_flow_autodamp", false)
      _webui_assert_value(unsupported_form, "power_flow_tol", "9.0e-7")
      _webui_assert_selected(unsupported_form, "power_flow_start_angle_mode", "dc")
      _webui_assert_selected(unsupported_form, "detailed_result_csv_format", "excel_us")
      @test occursin("case_settings_field_ignored", read(Sparlectra.webui_operation_log_path(root), String))

      fresh_root = mktempdir()
      mkpath(joinpath(fresh_root, "resolved"))
      write(joinpath(fresh_root, "resolved", "case145.m"), "% case fixture\n")
      fresh_form = _webui_test_form("case145.m", "configuration.yaml", fresh_root)
      fresh_form["power_flow_tol"] = "2.5e-7"
      fresh_form["power_flow_max_iter"] = "37"
      fresh_form["power_flow_qlimits_enforcement_mode"] = "active_set"
      fresh_form["detailed_result_csv"] = "on"
      fresh_form["detailed_result_csv_format"] = "excel_us"
      fresh_request = Sparlectra.powerflow_webui_request(fresh_form; default_output_root = fresh_root)
      fresh_runner = function(worker_request; case_directory = nothing)
        fresh_run_id = worker_request["run_id"]
        output_dir = joinpath(worker_request["output_root"], fresh_run_id)
        mkpath(output_dir)
        result = Sparlectra._api_result(
          run_id = fresh_run_id,
          status = :succeeded,
          success = true,
          converged = true,
          solution_available = true,
          iterations = 2,
          final_mismatch = 1.0e-9,
          reason = "converged",
          message = "ok",
          casefile = joinpath(fresh_root, "resolved", "case145.m"),
          config_file = "configuration.yaml",
          output_dir = output_dir,
          result_file = joinpath(output_dir, "result.json"),
          metadata = Dict{String,Any}(
            "runtime_casefile_path" => joinpath(fresh_root, "resolved", "case145.m"),
          ),
        )
        Sparlectra._POWERFLOW_SERVICE_RUNS[fresh_run_id] = result
        return Sparlectra.to_dict(result)
      end
      fresh_active = Sparlectra.start_webui_powerflow_run(fresh_request; runner = fresh_runner)
      fresh_run_id = fresh_active["run_id"]
      for _ in 1:100
        get(Sparlectra.get_webui_powerflow_job(fresh_run_id), "status", "") in Sparlectra._POWERFLOW_WEBUI_ACTIVE_STATES || break
        sleep(0.05)
      end
      fresh_completed = Sparlectra.get_webui_powerflow_job(fresh_run_id)
      @test get(fresh_completed, "status", "") == "success"
      @test get(get(fresh_completed, "metadata", Dict{String,Any}()), "webui_request_settings", nothing) isa AbstractDict
      fresh_response = Sparlectra.route_sparlectra_webui("POST", "/powerflow/result/$(fresh_run_id)/case-settings/save", Dict{String,String}(); output_root = fresh_root)
      @test fresh_response.status == 200
      fresh_profile_path = Sparlectra._webui_case_settings_path(fresh_root, joinpath(fresh_root, "resolved", "case145.m"))
      @test isfile(fresh_profile_path)
      fresh_profile = Sparlectra.load_yaml_dict(fresh_profile_path)
      fresh_settings = fresh_profile["settings"]
      @test Set(keys(fresh_settings)) == Set(Sparlectra._WEBUI_CASE_PROFILE_FIELDS)
      @test fresh_settings["power_flow_qlimits_enforcement_mode"] == "active_set"
      @test fresh_settings["power_flow_tol"] == 2.5e-7
      @test fresh_settings["power_flow_max_iter"] == 37
      @test fresh_settings["detailed_result_csv"] === true
      @test fresh_settings["detailed_result_csv_format"] == "excel_us"
      @test !haskey(fresh_settings, "effective_config")
      fresh_reloaded_form = String(Sparlectra.route_sparlectra_webui("GET", "/powerflow?casefile=$(Sparlectra._webui_urlencode(joinpath(fresh_root, "resolved", "case145.m")))"; output_root = fresh_root).body)
      _webui_assert_value(fresh_reloaded_form, "power_flow_tol", "2.5e-7")
      _webui_assert_value(fresh_reloaded_form, "power_flow_max_iter", "37")
      _webui_assert_selected(fresh_reloaded_form, "power_flow_qlimits_enforcement_mode", "active_set")
      _webui_assert_checked(fresh_reloaded_form, "detailed_result_csv", true)
      manual_override_form = _webui_test_form(joinpath(fresh_root, "resolved", "case145.m"), "configuration.yaml", fresh_root)
      manual_override_form["power_flow_tol"] = "3e-6"
      manual_override_request = Sparlectra.powerflow_webui_request(manual_override_form; default_output_root = fresh_root)
      @test manual_override_request["config_overrides"]["power_flow.tol"] == 3.0e-6
      manual_override_nested = Sparlectra.validate_gui_config_overrides(manual_override_request["config_overrides"])
      manual_override_cfg, _ = Sparlectra._load_api_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, manual_override_nested)
      @test manual_override_cfg.powerflow.tol == 3.0e-6
      delete!(Sparlectra._POWERFLOW_SERVICE_RUNS, fresh_run_id)
      delete!(Sparlectra._POWERFLOW_WEBUI_JOBS, fresh_run_id)

      failed_run_id = "case-settings-failed"
      failed = Sparlectra._api_result(
        run_id = failed_run_id,
        status = :not_converged,
        success = false,
        converged = false,
        solution_available = true,
        iterations = 80,
        final_mismatch = 1.0,
        reason = "not_converged",
        message = "not converged",
        casefile = joinpath(root, "case145.m"),
        config_file = "configuration.yaml",
        output_dir = joinpath(root, failed_run_id),
        metadata = metadata,
      )
      Sparlectra._POWERFLOW_SERVICE_RUNS[failed_run_id] = failed
      failed_html = Sparlectra.render_powerflow_result(Sparlectra.to_dict(failed))
      @test !occursin("Save settings for this case", failed_html)
      @test occursin("Save these settings anyway", failed_html)
      rejected = Sparlectra.handle_powerflow_case_settings_save(failed_run_id, Dict{String,String}(); output_root = root, operation_log = root)
      @test rejected.status == 400
      accepted = Sparlectra.handle_powerflow_case_settings_save(failed_run_id, Dict("override_non_success" => "true"); output_root = root, operation_log = root)
      @test accepted.status == 200
      log_text = read(Sparlectra.webui_operation_log_path(root), String)
      @test occursin("case_settings_saved", log_text)
      @test occursin("case_settings_save_failed", log_text)
      @test Sparlectra._webui_normalized_case_key("../bad/../../case145.m") == "case145"

      unsupported_run_id = "case-settings-unsupported"
      unsupported_metadata = deepcopy(metadata)
      unsupported_metadata["webui_request_settings"]["power_flow.autodamp"] = Dates.now()
      Sparlectra._POWERFLOW_WEBUI_JOBS[unsupported_run_id] = Dict{String,Any}(
        "run_id" => unsupported_run_id,
        "status" => "success",
        "success" => true,
        "converged" => true,
        "solution_available" => true,
        "iterations" => 3,
        "final_mismatch" => 1.0e-9,
        "reason" => "converged",
        "message" => "ok",
        "casefile" => joinpath(root, "case145.m"),
        "output_root" => root,
        "output_dir" => joinpath(root, unsupported_run_id),
        "metadata" => unsupported_metadata,
      )
      unsupported_response = Sparlectra.handle_powerflow_case_settings_save(unsupported_run_id, Dict{String,String}(); output_root = root, operation_log = root)
      @test unsupported_response.status == 400
      @test occursin("unsupported value type", String(unsupported_response.body))
      traversal_run_id = "case-settings-traversal"
      traversal_metadata = deepcopy(metadata)
      traversal_metadata["runtime_casefile_path"] = "../outside/case145.m"
      Sparlectra._POWERFLOW_WEBUI_JOBS[traversal_run_id] = Dict{String,Any}(
        "run_id" => traversal_run_id,
        "status" => "success",
        "success" => true,
        "converged" => true,
        "solution_available" => true,
        "iterations" => 3,
        "final_mismatch" => 1.0e-9,
        "reason" => "converged",
        "message" => "ok",
        "casefile" => "../outside/case145.m",
        "output_root" => root,
        "output_dir" => joinpath(root, traversal_run_id),
        "metadata" => traversal_metadata,
      )
      traversal_response = Sparlectra.handle_powerflow_case_settings_save(traversal_run_id, Dict{String,String}(); output_root = root, operation_log = root)
      @test traversal_response.status == 400
      delete!(Sparlectra._POWERFLOW_SERVICE_RUNS, run_id)
      delete!(Sparlectra._POWERFLOW_SERVICE_RUNS, failed_run_id)
      delete!(Sparlectra._POWERFLOW_WEBUI_JOBS, unsupported_run_id)
      delete!(Sparlectra._POWERFLOW_WEBUI_JOBS, traversal_run_id)
    end

    registry_before_warmup = Set(keys(Sparlectra._POWERFLOW_SERVICE_RUNS))
    warmup_output = Ref("")
    warmup_runner = function(; output_dir, kwargs...)
      warmup_output[] = output_dir
      write(joinpath(output_dir, "warmup-marker.txt"), "compiled\n")
      return (success = true, reason = nothing, message = nothing)
    end
    warmup_result = Sparlectra._run_sparlectra_webui_warmup(mktempdir(); runner = warmup_runner)
    @test warmup_result.success
    @test !ispath(warmup_output[])
    @test Set(keys(Sparlectra._POWERFLOW_SERVICE_RUNS)) == registry_before_warmup

    @testset "Case and configuration selection" begin
      mktempdir() do start_dir
        application_root = joinpath(start_dir, "Sparlectra")
        case_directory = joinpath(application_root, "data", "mpower")
        config_directory = joinpath(application_root, "examples")
        mkpath(case_directory)
        mkpath(config_directory)
        case14 = joinpath(case_directory, "case14.m")
        case118 = joinpath(case_directory, "case118.m")
        casejl = joinpath(case_directory, "case14.jl")
        for001 = joinpath(case_directory, "FOR001.DAT")
        for002 = joinpath(case_directory, "FOR002.DAT")
        for002_variant = joinpath(case_directory, "FOR002_reference.DAT")
        write(case14, "function mpc = case14\nend\n")
        write(case118, "function mpc = case118\nend\n")
        write(casejl, "nothing\n")
        write(for001, "0 / END\n")
        write(for002, "legacy reference\n")
        write(for002_variant, "legacy reference variant\n")
        for ext in ("csv", "json", "yaml", "log", "md")
          write(joinpath(case_directory, "artifact." * ext), "ignored\n")
        end
        write(joinpath(case_directory, "artifact.sparlectra-webui.yaml"), "power_flow: {}\n")
        primary_config = joinpath(config_directory, "configuration.yaml.example")
        secondary_config = joinpath(config_directory, "study.yaml.example")
        write(primary_config, "power_flow: {}\n")
        write(secondary_config, "power_flow: {}\n")
        write(joinpath(config_directory, "README.md"), "ignored\n")

        @test Sparlectra._webui_application_root(start_dir) == application_root
        @test Sparlectra._webui_casefile_options(application_root) == ["case118.m", "case14.jl", "case14.m", "FOR001.DAT"]
        @test Sparlectra._webui_for002_reference_options_in_directory(case_directory) == ["FOR002.DAT", "FOR002_reference.DAT"]
        @test Sparlectra._webui_config_file_options(application_root) == [primary_config, secondary_config]

        selection_html = Sparlectra.render_powerflow_form(
          application_root = application_root,
          selected_casefile = "case14.m",
          selected_config_file = secondary_config,
        )
        @test occursin("<select id=\"casefile\" name=\"casefile\" data-case-settings-reload=\"true\">", selection_html)
        @test occursin("<option value=\"\">-- choose existing case --</option>", selection_html)
        @test occursin("<option value=\"case14.m\" selected>case14.m</option>", selection_html)
        @test occursin("<option value=\"case14.jl\">case14.jl</option>", selection_html)
        @test occursin("<option value=\"FOR001.DAT\">FOR001.DAT</option>", selection_html)
        primary_selector_html = split(selection_html, "<datalist id=\"for002-reference-candidates\">")[1]
        @test !occursin("<option value=\"FOR002.DAT\">FOR002.DAT</option>", primary_selector_html)
        @test !occursin("<option value=\"FOR002_reference.DAT\">FOR002_reference.DAT</option>", primary_selector_html)
        @test occursin("<option value=\"case118.m\">case118.m</option>", selection_html)
        @test occursin("<input id=\"casefile_manual\" name=\"casefile_manual\" value=\"\" placeholder=\"case14.m or /path/to/FOR001.DAT\">", selection_html)
        @test !occursin("available-casefiles", selection_html)
        @test occursin("MATPOWER citation", selection_html)
        @test occursin("Zimmerman", selection_html)
        @test occursin("Murillo-Sanchez", selection_html)
        @test occursin("IEEE Transactions on Power Systems", selection_html)
        @test occursin("href=\"https://matpower.org\"", selection_html)
        @test occursin("href=\"https://matpower.org/citing/\"", selection_html)
        @test occursin("10.1109/TPWRS.2010.2051168", selection_html)
        @test occursin("ACTIVSg, PEGASE, and RTE", selection_html)
        @test occursin("<form id=\"powerflow-run-form\"", selection_html)
        @test occursin("Start PowerFlow run", selection_html)
        @test occursin("<input type=\"hidden\" name=\"config_file\" value=\"$(secondary_config)\">", selection_html)
        @test occursin("<code>$(secondary_config)</code>", selection_html)
        @test !occursin("README.md", selection_html)
        for ignored in ("artifact.csv", "artifact.json", "artifact.yaml", "artifact.log", "artifact.md", "artifact.sparlectra-webui.yaml")
          @test !occursin(ignored, selection_html)
        end
        @test occursin("Existing case file", selection_html)
        @test occursin("Or type case file path", selection_html)
        @test occursin("Default remains MATPOWER-oriented.", selection_html)
        @test occursin("DTF/FOR001 diagnostics (experimental/internal)", selection_html)
        @test occursin("name=\"for002_reference_file\"", selection_html)
        @test occursin("list=\"for002-reference-candidates\"", selection_html)
        @test occursin("<datalist id=\"for002-reference-candidates\"><option value=\"FOR002.DAT\">FOR002.DAT</option><option value=\"FOR002_reference.DAT\">FOR002_reference.DAT</option></datalist>", selection_html)
        @test occursin("FOR002.DAT is available in the case cache and can be selected here", selection_html)
        @test !occursin("full DTF/FOR001 support", selection_html)
        @test occursin("const updateDatCaseAssistance = function ()", selection_html)
        @test occursin("new RegExp('\\\\.dat\$', 'i').test(effectiveValue)", selection_html)
        @test occursin("caseFormat.value = 'dtf_for001'", selection_html)
        @test occursin("dtfInternalSection.open = true", selection_html)
        @test occursin("target.searchParams.set('casefile', caseSelect.value)", selection_html)

        dat_html = Sparlectra.render_powerflow_form(
          application_root = application_root,
          selected_casefile = "FOR001.DAT",
        )
        @test occursin("<option value=\"FOR001.DAT\" selected>FOR001.DAT</option>", dat_html)
        @test occursin("<input id=\"casefile_manual\" name=\"casefile_manual\" value=\"\"", dat_html)
        @test occursin("<option value=\"auto\">Auto</option>", dat_html)
        @test occursin("<option value=\"dtf_for001\" selected>DTF/FOR001 diagnostics (experimental/internal)</option>", dat_html)
        @test occursin("<details class=\"span-2 dtf-internal-section is-dat-selected\" open>", dat_html)
        @test occursin(".DAT selected:</strong> using internal DTF/FOR001 diagnostics.", dat_html)
        @test occursin("name=\"for002_reference_file\" value=\"\"", dat_html)
        @test !occursin("full DTF/FOR001 support", dat_html)

        for002_reference_html = Sparlectra.render_powerflow_form(
          application_root = application_root,
          selected_casefile = "FOR001.DAT",
          submitted_form = Dict("for002_reference_file" => "/manual/FOR002.DAT"),
        )
        @test occursin("name=\"for002_reference_file\" value=\"/manual/FOR002.DAT\"", for002_reference_html)
        @test !occursin("for002_reference_file\" value=\"FOR002.DAT\"", dat_html)

        submitted_auto_dat_html = Sparlectra.render_powerflow_form(
          application_root = application_root,
          selected_casefile = "FOR001.DAT",
          submitted_form = Dict("case_format" => "auto"),
        )
        @test occursin("<option value=\"auto\" selected>Auto</option>", submitted_auto_dat_html)
        @test occursin("<option value=\"dtf_for001\">DTF/FOR001 diagnostics (experimental/internal)</option>", submitted_auto_dat_html)

        case14_html = Sparlectra.render_powerflow_form(
          application_root = application_root,
          selected_casefile = "case14.m",
        )
        @test occursin("<option value=\"case14.m\" selected>case14.m</option>", case14_html)
        @test occursin("<option value=\"auto\" selected>Auto</option>", case14_html)
        @test occursin("<option value=\"dtf_for001\">DTF/FOR001 diagnostics (experimental/internal)</option>", case14_html)

        casejl_html = Sparlectra.render_powerflow_form(
          application_root = application_root,
          selected_casefile = "case14.jl",
        )
        @test occursin("<option value=\"case14.jl\" selected>case14.jl</option>", casejl_html)
        @test occursin("<option value=\"auto\" selected>Auto</option>", casejl_html)
        @test occursin("<option value=\"dtf_for001\">DTF/FOR001 diagnostics (experimental/internal)</option>", casejl_html)

        manual_dat = joinpath(case_directory, "manual", "FOR001.DAT")
        manual_dat_html = Sparlectra.render_powerflow_form(
          application_root = application_root,
          selected_casefile = manual_dat,
        )
        @test occursin("<input id=\"casefile_manual\" name=\"casefile_manual\" value=\"$(manual_dat)\"", manual_dat_html)
        @test occursin("<option value=\"dtf_for001\" selected>DTF/FOR001 diagnostics (experimental/internal)</option>", manual_dat_html)
        @test occursin("<details class=\"span-2 dtf-internal-section is-dat-selected\" open>", manual_dat_html)

        for manual_matpower in (joinpath(case_directory, "manual", "case14.m"), joinpath(case_directory, "manual", "case14.jl"))
          manual_matpower_html = Sparlectra.render_powerflow_form(
            application_root = application_root,
            selected_casefile = manual_matpower,
          )
          @test occursin("<input id=\"casefile_manual\" name=\"casefile_manual\" value=\"$(manual_matpower)\"", manual_matpower_html)
          @test occursin("<option value=\"auto\" selected>Auto</option>", manual_matpower_html)
          @test occursin("<option value=\"dtf_for001\">DTF/FOR001 diagnostics (experimental/internal)</option>", manual_matpower_html)
        end

        manual_overrides_dropdown_html = Sparlectra.render_powerflow_form(
          application_root = application_root,
          selected_casefile = "case14.m",
          submitted_form = Dict("casefile" => "case14.m", "casefile_manual" => manual_dat),
        )
        @test occursin("<input id=\"casefile_manual\" name=\"casefile_manual\" value=\"$(manual_dat)\"", manual_overrides_dropdown_html)
        @test occursin("<option value=\"dtf_for001\" selected>DTF/FOR001 diagnostics (experimental/internal)</option>", manual_overrides_dropdown_html)

        rm(case14)
        rm(case118)
        manual_case = joinpath(case_directory, "manual_case.m")
        fallback_html = Sparlectra.render_powerflow_form(
          application_root = application_root,
          selected_casefile = manual_case,
        )
        @test occursin("<select id=\"casefile\" name=\"casefile\" data-case-settings-reload=\"true\">", fallback_html)
        @test occursin("<input id=\"casefile_manual\" name=\"casefile_manual\" value=\"$(manual_case)\"", fallback_html)
      end
    end

    @testset "Standalone app-window launcher" begin
      url = "http://127.0.0.1:8080/powerflow"
      linux_lookup = name -> name == "google-chrome" ? "/opt/google/chrome" : nothing
      linux_selected = Sparlectra._webui_browser_open_command(url; platform = :linux, executable_lookup = linux_lookup)
      @test linux_selected !== nothing
      linux_command, linux_strategy = linux_selected
      @test linux_strategy == :app_window
      @test linux_command !== nothing
      @test linux_command.exec == ["/opt/google/chrome", "--app=$(url)", "--window-size=1500,950"]
      @test Sparlectra._webui_app_command(url; platform = :linux, executable_lookup = linux_lookup).exec == linux_command.exec

      xdg_lookup = name -> name == "xdg-open" ? "/usr/bin/xdg-open" : nothing
      xdg_selected = Sparlectra._webui_browser_open_command(url; platform = :linux, executable_lookup = xdg_lookup)
      @test xdg_selected !== nothing
      xdg_command, xdg_strategy = xdg_selected
      @test xdg_strategy == :xdg_open
      @test xdg_command.exec == ["/usr/bin/xdg-open", url]
      @test !any(arg -> startswith(arg, "--app="), xdg_command.exec)

      gio_lookup = name -> name == "gio" ? "/usr/bin/gio" : nothing
      gio_selected = Sparlectra._webui_browser_open_command(url; platform = :linux, executable_lookup = gio_lookup)
      @test gio_selected !== nothing
      gio_command, gio_strategy = gio_selected
      @test gio_strategy == :gio_open
      @test gio_command.exec == ["/usr/bin/gio", "open", url]
      @test !any(arg -> startswith(arg, "--app="), gio_command.exec)

      sensible_lookup = name -> name == "sensible-browser" ? "/usr/bin/sensible-browser" : nothing
      sensible_selected = Sparlectra._webui_browser_open_command(url; platform = :linux, executable_lookup = sensible_lookup)
      @test sensible_selected !== nothing
      sensible_command, sensible_strategy = sensible_selected
      @test sensible_strategy == :sensible_browser
      @test sensible_command.exec == ["/usr/bin/sensible-browser", url]
      @test !any(arg -> startswith(arg, "--app="), sensible_command.exec)

      firefox_lookup = name -> name == "firefox" ? "/usr/bin/firefox" : nothing
      firefox_selected = Sparlectra._webui_browser_open_command(url; platform = :linux, executable_lookup = firefox_lookup)
      @test firefox_selected === nothing

      windows_lookup = name -> name == "msedge.exe" ? raw"C:\Browser\msedge.exe" : nothing
      windows_command = Sparlectra._webui_app_command(url; platform = :windows, executable_lookup = windows_lookup, environment = Dict{String,String}())
      @test windows_command !== nothing
      @test windows_command.exec == [raw"C:\Browser\msedge.exe", "--app=$(url)", "--window-size=1500,950"]

      macos_exists = path -> path == "/Applications/Google Chrome.app"
      macos_command = Sparlectra._webui_app_command(url; platform = :macos, path_exists = macos_exists)
      @test macos_command !== nothing
      @test macos_command.exec == ["open", "-na", "/Applications/Google Chrome.app", "--args", "--app=$(url)", "--window-size=1500,950"]

      missing_lookup = _ -> nothing
      @test Sparlectra._webui_app_command(url; platform = :linux, executable_lookup = missing_lookup) === nothing
      @test Sparlectra._webui_browser_open_command(url; platform = :linux, executable_lookup = missing_lookup) === nothing
    end

    @testset "Markdown-backed contextual help and documentation" begin
      topic = "power_flow.start_mode.voltage_mode"
      @test haskey(Sparlectra.WEBUI_HELP_TOPICS, topic)
      topic_metadata = Sparlectra.resolve_webui_help_topic(topic)
      @test topic_metadata !== nothing
      doc_metadata = Sparlectra.resolve_webui_doc_page(topic_metadata.page)
      @test doc_metadata !== nothing
      doc_path = joinpath(@__DIR__, "..", "docs", "src", doc_metadata.file)
      @test isfile(doc_path)

      markdown_text = Sparlectra.load_webui_markdown_document(topic_metadata.page)
      @test markdown_text !== nothing
      section = Sparlectra.extract_webui_markdown_section(markdown_text, topic_metadata.heading)
      @test section !== nothing
      @test occursin(topic_metadata.selector, section)
      excerpt = Sparlectra.load_webui_help_excerpt(topic)
      @test excerpt !== nothing
      @test occursin("power_flow.start_mode.voltage_mode", excerpt)
      @test occursin("YAML path", excerpt)

      synthetic = "## Target\n\n- item\n\n```julia\nx = 1\n```\n\n### Child\nchild text\n\n## Next\nnot included"
      extracted = Sparlectra.extract_webui_markdown_section(synthetic, "Target")
      @test occursin("- item", extracted)
      @test occursin("```julia", extracted)
      @test occursin("### Child", extracted)
      @test !occursin("not included", extracted)

      link_markdown = "## Local section\n\n[Configuration](powerflow_configuration.md)\n[Start modes](powerflow_configuration.md#start-mode-options)\n[Dot relative](./powerflow_configuration.md)\n[External](https://example.com/reference)\n[Local anchor](#local-section)\n[Unknown](unknown.md)"
      link_html = Sparlectra.render_webui_markdown(link_markdown; current_page = "webui")
      @test occursin("href=\"/docs/powerflow_configuration\"", link_html)
      @test occursin("href=\"/docs/powerflow_configuration#start-mode-options\"", link_html)
      @test count("href=\"/docs/powerflow_configuration\"", link_html) == 2
      @test occursin("href=\"https://example.com/reference\"", link_html)
      @test occursin("href=\"#local-section\"", link_html)
      @test occursin("id=\"local-section\"", link_html)
      @test !occursin("href=\"unknown.md\"", link_html)
      @test occursin("aria-disabled=\"true\"", link_html)

      for unsafe_target in ("../Project.toml", "../../src/Sparlectra.jl", "/etc/passwd", raw"C:\Users\x\file.txt")
        unsafe_html = Sparlectra.render_webui_markdown("[unsafe]($(unsafe_target))"; current_page = "webui")
        @test occursin("aria-disabled=\"true\"", unsafe_html)
        @test !occursin("href=\"$(unsafe_target)\"", unsafe_html)
      end

      @test Set(values(Sparlectra.WEBUI_FORM_HELP_TOPICS)) == Set(keys(Sparlectra.WEBUI_HELP_TOPICS))
      for (configured_topic, metadata) in Sparlectra.WEBUI_HELP_TOPICS
        page_metadata = Sparlectra.resolve_webui_doc_page(metadata.page)
        @test page_metadata !== nothing
        help_doc_path = joinpath(@__DIR__, "..", "docs", "src", page_metadata.file)
        @test isfile(help_doc_path)
        help_markdown = Sparlectra.load_webui_markdown_document(metadata.page)
        @test help_markdown !== nothing
        help_section = Sparlectra.extract_webui_markdown_section(help_markdown, metadata.heading)
        @test help_section !== nothing
        isempty(metadata.selector) || @test occursin(metadata.selector, help_section)
        @test Sparlectra.load_webui_help_excerpt(configured_topic) !== nothing
      end

      current_iteration_help = Dict(
        "power_flow.start_current_iteration.enabled" => ("not a separate power-flow solver", "start-value preconditioner", "before the Newton-Raphson power-flow solver starts", "accepted only if it passes the voltage and angle guards", "current_iteration_start.log"),
        "power_flow.start_current_iteration.max_iter" => ("maximum number of current-iteration pre-solve steps", "extra time", "Default: 10", "pre-solve stops too early"),
        "power_flow.start_current_iteration.tol" => ("stopping tolerance", "not the final Newton-Raphson power-flow tolerance", "Default: 1.0e-3", "improve the starting point"),
        "power_flow.start_current_iteration.damping" => ("damping factor", "Smaller values blend the update", "Default: 0.5", "rejected by voltage or angle guards"),
        "power_flow.start_current_iteration.accept_only_if_improved" => ("improves the existing Sparlectra mismatch metric", "restores the original start values", "Default: enabled", "expert experiments"),
        "power_flow.start_current_iteration.min_improvement_factor" => ("required improvement ratio", "0.98 means", "Default: 0.98", "clearly better starts"),
        "power_flow.start_current_iteration.vm_min_pu" => ("lower voltage-magnitude guard", "candidate is rejected", "Default: 0.5 pu", "candidate voltage minima"),
        "power_flow.start_current_iteration.vm_max_pu" => ("upper voltage-magnitude guard", "unrealistic over-voltage", "Default: 1.5 pu", "candidate voltage maxima"),
        "power_flow.start_current_iteration.max_angle_step_deg" => ("maximum allowed angle change", "angle jump larger than this limit", "Default: 30 degrees", "angle-step guard"),
        "power_flow.start_current_iteration.only_for_large_cases" => ("classifies as large enough", "avoids spending time on small cases", "Default: disabled", "difficult large MATPOWER cases"),
      )
      for (current_iteration_topic, required_fragments) in current_iteration_help
        current_iteration_excerpt = Sparlectra.load_webui_help_excerpt(current_iteration_topic)
        @test current_iteration_excerpt !== nothing
        @test occursin(current_iteration_topic, current_iteration_excerpt)
        @test !startswith(strip(current_iteration_excerpt), "|")
        for fragment in required_fragments
          @test occursin(fragment, current_iteration_excerpt)
        end
      end

      help_response = Sparlectra.route_sparlectra_webui("GET", "/help/$(topic)")
      @test help_response.status == 200
      help_html = String(help_response.body)
      @test occursin("Start voltage mode", help_html)
      @test occursin("power_flow.start_mode.voltage_mode", help_html)
      @test occursin("/docs/powerflow_configuration", help_html)
      @test occursin("class=\"panel help-page help-panel\"", help_html)

      @test occursin("class=\"button back-button\"", help_html)
      @test occursin("document.referrer.startsWith(location.origin)", help_html)
      @test occursin("history.back()", help_html)
      @test occursin("href=\"/powerflow\"", help_html)

      for representative_topic in (
        "power_flow.start_mode.voltage_mode",
        "power_flow.start_mode.angle_mode",
        "power_flow.wrong_branch_detection",
        "benchmark.enabled",
      )
        representative_help = Sparlectra.route_sparlectra_webui("GET", "/help/$(representative_topic)")
        @test representative_help.status == 200
        @test occursin("class=\"panel help-page help-panel\"", String(representative_help.body))
      end

      unknown_help = Sparlectra.route_sparlectra_webui("GET", "/help/power_flow.unknown")
      @test unknown_help.status == 404
      @test occursin("Unknown help topic", String(unknown_help.body))

      docs_index = Sparlectra.route_sparlectra_webui("GET", "/docs")
      @test docs_index.status == 200
      docs_index_html = String(docs_index.body)
      for page in keys(Sparlectra.WEBUI_DOC_PAGES)
        @test occursin("/docs/$(page)", docs_index_html)
      end

      docs_page = Sparlectra.route_sparlectra_webui("GET", "/docs/powerflow_configuration")
      @test docs_page.status == 200
      docs_page_html = String(docs_page.body)
      @test occursin("Power-Flow Configuration", docs_page_html)
      @test occursin("class=\"panel docs-page docs-content\"", docs_page_html)
      @test occursin("id=\"start-mode-options\"", docs_page_html)

      qlimit_strategy_page = Sparlectra.route_sparlectra_webui("GET", "/docs/q_limit_switching_strategy")
      @test qlimit_strategy_page.status == 200
      qlimit_strategy_html = String(qlimit_strategy_page.body)
      for term in ("active_set", "classic_simultaneous", "classic_one_at_a_time")
        @test occursin(term, qlimit_strategy_html)
      end

      matpower_format_page = Sparlectra.route_sparlectra_webui("GET", "/docs/matpower_format")
      @test matpower_format_page.status == 200
      matpower_format_html = String(matpower_format_page.body)
      @test occursin("MATPOWER format", matpower_format_html)
      @test occursin("MATPOWER Manual", matpower_format_html)
      @test occursin("href=\"https://matpower.app/manual/matpower/DataFileFormat.html\" target=\"_blank\" rel=\"noopener noreferrer\"", matpower_format_html)
      @test occursin("href=\"https://matpower.org/documentation/ref-manual/legacy/functions/caseformat.html\" target=\"_blank\" rel=\"noopener noreferrer\"", matpower_format_html)

      webui_docs_page = Sparlectra.route_sparlectra_webui("GET", "/docs/webui")
      @test webui_docs_page.status == 200
      webui_docs_html = String(webui_docs_page.body)
      @test occursin("href=\"/docs/powerflow_configuration\"", webui_docs_html)
      @test occursin("href=\"/docs/powerflow_configuration#start-mode-options\"", webui_docs_html)
      @test occursin("href=\"/docs/performance_profiling\"", webui_docs_html)
      @test !occursin("href=\"/docs/powerflow_configuration\" target=\"_blank\"", webui_docs_html)

      for unsafe_page in ("../Project.toml", "../../src/Sparlectra.jl", "/etc/passwd", raw"C:\Users\x\file.txt")
        unsafe = Sparlectra.handle_webui_doc_page(unsafe_page)
        @test unsafe.status == 404
        @test occursin("Documentation page not found", String(unsafe.body))
        unsafe_route = Sparlectra.route_sparlectra_webui("GET", "/docs/" * Sparlectra._webui_urlencode(unsafe_page))
        @test unsafe_route.status == 404
      end
    end

    mktempdir() do tmpdir
      casefile = _write_webui_test_case(joinpath(tmpdir, "case_webui.m"))
      config_file = joinpath(tmpdir, "config_template.yaml")
      cp(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, config_file)
      output_root = joinpath(tmpdir, "runs")
      form = _webui_test_form(casefile, config_file, output_root)

      request = Sparlectra.powerflow_webui_request(form; default_output_root = output_root)
      @test request["casefile"] == casefile
      existing_form = copy(form)
      existing_form["casefile"] = "  case118.jl  "
      existing_form["casefile_manual"] = "  "
      existing_request = Sparlectra.powerflow_webui_request(existing_form; default_output_root = output_root)
      @test existing_request["casefile"] == "case118.jl"
      @test existing_request["case_format"] == "auto"

      typed_form = copy(form)
      typed_form["casefile"] = "case9.m"
      typed_form["casefile_manual"] = "  case9241pegase.m  "
      typed_request = Sparlectra.powerflow_webui_request(typed_form; default_output_root = output_root)
      @test typed_request["casefile"] == "case9241pegase.m"
      @test typed_request["case_format"] == "auto"

      dtf_form = copy(form)
      dtf_form["casefile"] = "FOR001.DAT"
      dtf_form["casefile_manual"] = "  "
      dtf_form["case_format"] = "dtf_for001"
      dtf_request = Sparlectra.powerflow_webui_request(dtf_form; default_output_root = output_root)
      @test dtf_request["casefile"] == "FOR001.DAT"
      @test dtf_request["case_format"] == "dtf_for001"
      @test dtf_request["for002_reference_file"] === nothing

      auto_dtf_form = copy(dtf_form)
      auto_dtf_form["case_format"] = "auto"
      auto_dtf_request = Sparlectra.powerflow_webui_request(auto_dtf_form; default_output_root = output_root)
      @test auto_dtf_request["casefile"] == "FOR001.DAT"
      @test auto_dtf_request["case_format"] == "auto"

      manual_dtf = joinpath(tmpdir, "FOR001.DAT")
      manual_dtf_form = copy(dtf_form)
      manual_dtf_form["casefile"] = "case9.m"
      manual_dtf_form["casefile_manual"] = "  " * manual_dtf * "  "
      manual_dtf_request = Sparlectra.powerflow_webui_request(manual_dtf_form; default_output_root = output_root)
      @test manual_dtf_request["casefile"] == manual_dtf
      @test manual_dtf_request["case_format"] == "dtf_for001"

      reference_dtf_form = copy(dtf_form)
      reference_dtf_form["for002_reference_file"] = "  /manual/FOR002.DAT  "
      reference_dtf_request = Sparlectra.powerflow_webui_request(reference_dtf_form; default_output_root = output_root)
      @test reference_dtf_request["casefile"] == "FOR001.DAT"
      @test reference_dtf_request["for002_reference_file"] == "/manual/FOR002.DAT"

      cache_for002 = joinpath(tmpdir, "FOR002.DAT")
      write(cache_for002, "legacy reference\n")
      write(manual_dtf, "0 / END\n")
      base_service_request = Dict{String,Any}(
        "config_file" => config_file,
        "output_root" => output_root,
        "case_format" => "dtf_for001",
      )
      for002_primary_response = Sparlectra.start_powerflow_run(merge(copy(base_service_request), Dict{String,Any}("casefile" => "FOR002.DAT")); case_directory = tmpdir)
      @test !for002_primary_response["success"]
      @test for002_primary_response["reason"] == "invalid_casefile"
      @test occursin("FOR002.DAT is a reference/result file and cannot be used as the primary DTF/FOR001 input case.", for002_primary_response["message"])

      manual_for002_primary_response = Sparlectra.start_powerflow_run(merge(copy(base_service_request), Dict{String,Any}("casefile" => cache_for002)); case_directory = tmpdir)
      @test !manual_for002_primary_response["success"]
      @test manual_for002_primary_response["reason"] == "invalid_casefile"
      @test occursin("Use FOR001.DAT as the case and enter FOR002.DAT as optional FOR002 reference file.", manual_for002_primary_response["message"])

      empty_case_form = copy(form)
      empty_case_form["casefile"] = ""
      empty_case_form["casefile_manual"] = ""
      @test_throws ArgumentError Sparlectra.powerflow_webui_request(empty_case_form; default_output_root = output_root)
      empty_case_response = Sparlectra.route_sparlectra_webui("POST", "/powerflow/run", empty_case_form; output_root = output_root)
      @test empty_case_response.status == 400
      @test occursin("Select an existing MATPOWER case or type a case name.", String(empty_case_response.body))
      operation_log_path = Sparlectra.webui_operation_log_path(output_root)
      @test isfile(operation_log_path)
      @test occursin("\"event\":\"validation_error\"", read(operation_log_path, String))
      @test request["config_file"] == config_file
      @test request["output_root"] == output_root
      @test request["detailed_result_csv"] === true
      @test request["detailed_result_csv_format"] == "excel_de"
      @test request["config_overrides"]["power_flow.qlimits.enabled"] === true
      @test request["config_overrides"]["power_flow.qlimits.enforcement_mode"] == "classic_simultaneous"
      @test request["config_overrides"]["power_flow.start_current_iteration.enabled"] === true
      @test request["config_overrides"]["power_flow.start_current_iteration.max_iter"] == 7
      @test request["config_overrides"]["power_flow.start_current_iteration.tol"] == 1.0e-4
      qlimits_disabled_form = copy(form)
      qlimits_disabled_form["power_flow_qlimits_enabled"] = "false"
      @test Sparlectra.powerflow_webui_request(qlimits_disabled_form; default_output_root = output_root)["config_overrides"]["power_flow.qlimits.enabled"] === false
      missing_field_form = copy(form)
      delete!(missing_field_form, "power_flow_start_current_iteration_enabled")
      @test !haskey(Sparlectra.powerflow_webui_request(missing_field_form; default_output_root = output_root)["config_overrides"], "power_flow.start_current_iteration.enabled")
      ci_disabled_form = copy(form)
      ci_disabled_form["power_flow_start_current_iteration_enabled"] = "false"
      @test Sparlectra.powerflow_webui_request(ci_disabled_form; default_output_root = output_root)["config_overrides"]["power_flow.start_current_iteration.enabled"] === false
      csv_disabled_form = copy(form)
      csv_disabled_form["detailed_result_csv"] = "false"
      @test Sparlectra.powerflow_webui_request(csv_disabled_form; default_output_root = output_root)["detailed_result_csv"] === false
      delete!(csv_disabled_form, "detailed_result_csv_format")
      @test Sparlectra.powerflow_webui_request(csv_disabled_form; default_output_root = output_root)["detailed_result_csv_format"] == "excel_us"
      untrusted_form = copy(form)
      untrusted_form["output_root"] = joinpath(tmpdir, "outside")
      @test Sparlectra.powerflow_webui_request(untrusted_form; default_output_root = output_root)["output_root"] == output_root
      overrides = request["config_overrides"]
      @test Set(keys(overrides)) == Set(key for (key, _, _) in Sparlectra._WEBUI_FORM_CONFIG_FIELDS)
      @test overrides["power_flow.tol"] == 1.0e-8
      @test overrides["power_flow.max_iter"] == 80
      @test overrides["power_flow.autodamp"]
      apply_form = copy(form)
      apply_form["matpower_import_auto_profile"] = "apply"
      apply_request = Sparlectra.powerflow_webui_request(apply_form; default_output_root = output_root)
      @test apply_request["config_overrides"]["matpower_import.auto_profile"] == "apply"
      @test overrides["benchmark.enabled"] === false

      invalid_form = copy(form)
      invalid_form["power_flow_max_iter"] = "invalid"
      invalid_response = Sparlectra.route_sparlectra_webui("POST", "/powerflow/run", invalid_form; output_root = output_root)
      invalid_html = String(invalid_response.body)
      @test invalid_response.status == 400
      @test occursin("name=\"casefile_manual\"", invalid_html)
      @test occursin("name=\"config_file\" value=\"$(config_file)\"", invalid_html)
      @test occursin("data-dismissible-alert", invalid_html)
      @test occursin("class=\"alert alert-error error\"", invalid_html)
      @test occursin("class=\"alert-close\"", invalid_html)
      @test occursin("data-powerflow-form", invalid_html)
      @test occursin("form.addEventListener('input', clearAlert, {once: true})", invalid_html)
      @test occursin("form.addEventListener('change', clearAlert, {once: true})", invalid_html)
      @test !occursin("<details class=\"last-errors span-2\">", invalid_html)
      @test occursin("Last errors", invalid_html)
      @test occursin("href=\"/webui/last-errors\"", invalid_html)
      last_errors_response = Sparlectra.route_sparlectra_webui("GET", "/webui/last-errors"; output_root = output_root)
      last_errors_html = String(last_errors_response.body)
      @test last_errors_response.status == 200
      @test occursin("Last errors", last_errors_html)
      @test occursin("validation_error", read(operation_log_path, String))
      @test occursin("/powerflow/run", last_errors_html)
      empty_errors_html = Sparlectra.render_webui_last_errors(joinpath(tmpdir, "missing-operation-log.jsonl"))
      @test occursin("No recent errors.", empty_errors_html)

      api_failure_root = joinpath(tmpdir, "api-failure-runs")
      api_failure_case = _write_webui_test_case(joinpath(tmpdir, "case_with_dcline.m"))
      api_failure_form = _webui_test_form(api_failure_case, config_file, api_failure_root)
      api_failure_message = "Unsupported MATPOWER dcline data: active mpc.dcline entries detected. DC lines are currently not supported; the run was aborted before power-flow solution."
      api_failure_runner = function(worker_request; case_directory = nothing)
        api_failure_run_id = worker_request["run_id"]
        api_failure_output_dir = joinpath(worker_request["output_root"], api_failure_run_id)
        mkpath(api_failure_output_dir)
        result = Sparlectra._api_result(
          run_id = api_failure_run_id,
          status = :failed,
          success = false,
          converged = false,
          solution_available = false,
          iterations = nothing,
          final_mismatch = nothing,
          reason = "unsupported_matpower_dcline",
          message = api_failure_message,
          casefile = api_failure_case,
          config_file = config_file,
          output_dir = api_failure_output_dir,
          result_file = joinpath(api_failure_output_dir, "result.json"),
          metadata = Dict{String,Any}("failure_reason" => "unsupported_matpower_dcline"),
        )
        Sparlectra._POWERFLOW_SERVICE_RUNS[api_failure_run_id] = result
        return Sparlectra.to_dict(result)
      end
      api_failure_response = Sparlectra.handle_powerflow_run(api_failure_form; default_output_root = api_failure_root, runner = api_failure_runner, operation_log = api_failure_root)
      api_failure_run_id = api_failure_response["run_id"]
      wait(Sparlectra._POWERFLOW_WEBUI_JOBS[api_failure_run_id]["task"])
      api_failure_job = Sparlectra.get_webui_powerflow_job(api_failure_run_id)
      @test api_failure_job["status"] == "failed"
      @test api_failure_job["run_status"] == "failed"
      @test api_failure_job["reason"] == "unsupported_matpower_dcline"
      api_failure_result_html = String(Sparlectra.route_sparlectra_webui("GET", "/powerflow/result/$(api_failure_run_id)"; output_root = api_failure_root).body)
      @test occursin(api_failure_message, api_failure_result_html)
      api_failure_last_errors_html = String(Sparlectra.route_sparlectra_webui("GET", "/webui/last-errors"; output_root = api_failure_root).body)
      @test occursin("Last errors", api_failure_last_errors_html)
      @test occursin("unsupported_matpower_dcline", api_failure_last_errors_html)
      @test occursin("case_with_dcline.m", api_failure_last_errors_html)
      @test occursin("Unsupported MATPOWER dcline data", api_failure_last_errors_html)
      @test occursin("/powerflow/result/$(api_failure_run_id)", api_failure_last_errors_html)

      api_success_root = joinpath(tmpdir, "api-success-runs")
      api_success_form = _webui_test_form(casefile, config_file, api_success_root)
      api_success_runner = function(worker_request; case_directory = nothing)
        api_success_run_id = worker_request["run_id"]
        api_success_output_dir = joinpath(worker_request["output_root"], api_success_run_id)
        mkpath(api_success_output_dir)
        result = Sparlectra._api_result(
          run_id = api_success_run_id,
          status = :succeeded,
          success = true,
          converged = true,
          solution_available = true,
          iterations = 2,
          final_mismatch = 1.0e-9,
          reason = "converged",
          message = "PowerFlow run completed.",
          casefile = casefile,
          config_file = config_file,
          output_dir = api_success_output_dir,
          result_file = joinpath(api_success_output_dir, "result.json"),
        )
        Sparlectra._POWERFLOW_SERVICE_RUNS[api_success_run_id] = result
        return Sparlectra.to_dict(result)
      end
      api_success_response = Sparlectra.handle_powerflow_run(api_success_form; default_output_root = api_success_root, runner = api_success_runner, operation_log = api_success_root)
      api_success_run_id = api_success_response["run_id"]
      wait(Sparlectra._POWERFLOW_WEBUI_JOBS[api_success_run_id]["task"])
      api_success_last_errors_html = String(Sparlectra.route_sparlectra_webui("GET", "/webui/last-errors"; output_root = api_success_root).body)
      @test occursin("No recent errors.", api_success_last_errors_html)

      form_html = Sparlectra.render_powerflow_form(output_root = output_root)
      @test occursin("<option value=\"off\">off</option>", form_html)
      @test occursin("<option value=\"recommend\" selected>recommend</option>", form_html)
      @test occursin("<option value=\"apply\">apply</option>", form_html)
      @test occursin("<legend>MATPOWER import conventions</legend>", form_html)
      @test !occursin("matpower_simultaneous", form_html)
      @test !occursin("matpower_one_at_a_time", form_html)
      @test findfirst("Advanced / expert options", form_html) < findfirst("MATPOWER import conventions", form_html)
      @test findfirst("<details class=\"span-2 expert-section\">", form_html) < findfirst("MATPOWER import conventions", form_html)
      @test occursin("Use <strong>off</strong> for manual overrides", form_html)
      @test findfirst("Advanced / expert options", form_html) < findfirst("name=\"run_diagnostics\"", form_html)
      routed_powerflow_html = String(Sparlectra.route_sparlectra_webui("GET", "/powerflow"; output_root = output_root).body)
      @test occursin("MATPOWER import conventions", routed_powerflow_html)
      @test occursin("name=\"matpower_import_auto_profile\"", routed_powerflow_html)
      @test occursin("name=\"matpower_import_ratio\"", routed_powerflow_html)
      @test occursin("name=\"matpower_import_shift_sign\"", routed_powerflow_html)
      @test occursin("name=\"matpower_import_shift_unit\"", routed_powerflow_html)
      @test occursin("Package path:", routed_powerflow_html)
      @test occursin("commit", routed_powerflow_html)
      @test occursin("name=\"power_flow_qlimits_enforcement_mode\"", form_html)
      for mode in ("active_set", "classic_simultaneous", "classic_one_at_a_time")
        @test occursin("<option value=\"$(mode)\"", form_html)
      end
      @test !occursin("<option value=\"matpower_simultaneous\"", form_html)
      @test !occursin("<option value=\"matpower_one_at_a_time\"", form_html)
      @test Sparlectra.resolve_webui_help_topic("power_flow.qlimits.enforcement_mode") !== nothing
      expected_help_topics = Dict(
        "casefile" => "webui.casefile",
        "power_flow_tol" => "power_flow.tol",
        "power_flow_max_iter" => "power_flow.max_iter",
        "power_flow_autodamp" => "power_flow.autodamp",
        "power_flow_autodamp_min" => "power_flow.autodamp_min",
        "power_flow_qlimits_enabled" => "power_flow.qlimits.enabled",
        "power_flow_qlimits_enforcement_mode" => "power_flow.qlimits.enforcement_mode",
        "power_flow_wrong_branch_detection" => "power_flow.wrong_branch_detection",
        "power_flow_start_angle_mode" => "power_flow.start_mode.angle_mode",
        "power_flow_start_voltage_mode" => "power_flow.start_mode.voltage_mode",
        "power_flow_start_current_iteration_enabled" => "power_flow.start_current_iteration.enabled",
        "power_flow_start_current_iteration_max_iter" => "power_flow.start_current_iteration.max_iter",
        "power_flow_start_current_iteration_tol" => "power_flow.start_current_iteration.tol",
        "power_flow_start_current_iteration_damping" => "power_flow.start_current_iteration.damping",
        "power_flow_start_current_iteration_accept_only_if_improved" => "power_flow.start_current_iteration.accept_only_if_improved",
        "power_flow_start_current_iteration_min_improvement_factor" => "power_flow.start_current_iteration.min_improvement_factor",
        "power_flow_start_current_iteration_vm_min_pu" => "power_flow.start_current_iteration.vm_min_pu",
        "power_flow_start_current_iteration_vm_max_pu" => "power_flow.start_current_iteration.vm_max_pu",
        "power_flow_start_current_iteration_max_angle_step_deg" => "power_flow.start_current_iteration.max_angle_step_deg",
        "power_flow_start_current_iteration_only_for_large_cases" => "power_flow.start_current_iteration.only_for_large_cases",
        "matpower_import_auto_profile" => "matpower_import.auto_profile",
        "matpower_import_ratio" => "matpower_import.ratio",
        "matpower_import_shift_sign" => "matpower_import.shift_sign",
        "matpower_import_shift_unit" => "matpower_import.shift_unit",
        "matpower_import_bus_shunt_model" => "matpower_import.bus_shunt_model",
        "matpower_import_pv_voltage_source" => "matpower_import.pv_voltage_source",
        "matpower_import_compare_voltage_reference" => "matpower_import.compare_voltage_reference",
        "output_logfile_results" => "output.logfile_results",
        "benchmark_enabled" => "benchmark.enabled",
        "benchmark_samples" => "benchmark.samples",
        "benchmark_seconds" => "benchmark.seconds",
        "performance_timing" => "webui.performance_timing",
        "run_diagnostics" => "webui.run_diagnostics",
        "detailed_result_csv" => "webui.detailed_result_csv",
        "detailed_result_csv_format" => "webui.detailed_result_csv_format",
      )
      @test all(Sparlectra.WEBUI_FORM_HELP_TOPICS[field] == help_topic for (field, help_topic) in expected_help_topics)
      for (field, help_topic) in expected_help_topics
        @test occursin("name=\"$(field)\"", form_html)
        @test occursin("href=\"/help/$(help_topic)\"", form_html)
      end
      @test !occursin("name=\"output_root\"", form_html)
      @test occursin("Output root", form_html)
      @test occursin(output_root, form_html)
      @test occursin("action=\"/webui/shutdown\"", form_html)
      @test occursin("Stop Web UI", form_html)
      @test !occursin("navigator.sendBeacon('/webui/shutdown'", form_html)
      @test !occursin("fetch('/webui/shutdown', {method: 'POST', keepalive: true}", form_html)
      @test !occursin("pagehide", form_html)
      @test !occursin("beforeunload", form_html)
      @test !occursin("visibilitychange", form_html)
      @test occursin("href=\"/webui/operation-log\"", form_html)
      @test count("class=\"help-link\"", form_html) == length(expected_help_topics)

      @test !occursin("class=\"button back-button\"", form_html)

      @test occursin("id=\"powerflow-run-form\"", form_html)
      @test occursin("onsubmit=\"this.classList.add('is-submitting')", form_html)
      @test occursin("this.setAttribute('aria-busy', 'true')", form_html)
      @test occursin("this.querySelector('button[type=submit]').disabled = true", form_html)
      @test occursin("window.addEventListener('pageshow'", form_html)
      @test occursin("form.classList.remove('is-submitting')", form_html)
      @test occursin("form.removeAttribute('aria-busy')", form_html)
      @test occursin("submitButton.disabled = false", form_html)
      @test occursin("class=\"powerflow-submit\"", form_html)
      @test occursin("class=\"submit-spinner\" aria-hidden=\"true\"", form_html)
      @test occursin("class=\"submit-progress-label\" role=\"status\" aria-live=\"polite\">Running PowerFlow…", form_html)

      @test occursin("src=\"/assets/logo.png\"", form_html)
      @test occursin("alt=\"Sparlectra.jl logo\"", form_html)
      @test occursin("Sparlectra.jl v", form_html)
      @test occursin("Sparlectra.jl v$(Sparlectra.version())", form_html)
      @test occursin("name=\"performance_timing\"", form_html)
      @test occursin("name=\"run_diagnostics\"", form_html)
      @test occursin("Advanced start values", form_html)
      @test occursin("<details class=\"span-2 start-current-iteration-options advanced-start-values\">", form_html)
      @test occursin("<summary>Advanced start values</summary>", form_html)
      @test occursin("Current-iteration pre-solve", form_html)
      for field in (
        "power_flow_start_current_iteration_enabled",
        "power_flow_start_current_iteration_max_iter",
        "power_flow_start_current_iteration_tol",
        "power_flow_start_current_iteration_damping",
        "power_flow_start_current_iteration_accept_only_if_improved",
        "power_flow_start_current_iteration_min_improvement_factor",
        "power_flow_start_current_iteration_vm_min_pu",
        "power_flow_start_current_iteration_vm_max_pu",
        "power_flow_start_current_iteration_max_angle_step_deg",
        "power_flow_start_current_iteration_only_for_large_cases",
      )
        @test occursin("name=\"$(field)\"", form_html)
      end
      @test occursin("name=\"power_flow_start_current_iteration_enabled\" type=\"hidden\" value=\"false\"", form_html)
      @test occursin("name=\"power_flow_start_current_iteration_accept_only_if_improved\" type=\"hidden\" value=\"false\"", form_html)
      @test occursin("name=\"power_flow_start_current_iteration_only_for_large_cases\" type=\"hidden\" value=\"false\"", form_html)
      @test occursin("name=\"detailed_result_csv\" type=\"checkbox\" value=\"true\" checked", form_html)
      @test occursin("class=\"span-2 detailed-csv-options\"", form_html)
      @test occursin("name=\"detailed_result_csv_format\"", form_html)
      @test occursin("<option value=\"technical\">", form_html)
      @test occursin("<option value=\"excel_de\">", form_html)
      @test occursin("<option value=\"excel_us\" selected>", form_html)
      @test occursin("navigator.languages", form_html)
      @test occursin("startsWith('de')", form_html)
      @test !occursin("Use Excel CSV format with semicolon delimiter", form_html)

      logo_response = Sparlectra.route_sparlectra_webui("GET", "/assets/logo.png")
      @test logo_response.status == 200
      @test ("Content-Type" => "image/png") in logo_response.headers
      @test logo_response.body == read(joinpath(@__DIR__, "..", "docs", "src", "assets", "logo.png"))
      for target in ("/assets/tablestyle.css", "/assets/../Project.toml", "/assets/logo.png/extra")
        @test Sparlectra.route_sparlectra_webui("GET", target).status == 404
      end
      routed_form = copy(form)
routed_form["output_root"] = joinpath(tmpdir, "outside-runs")
run_response = Sparlectra.route_sparlectra_webui("POST", "/powerflow/run", routed_form; output_root)
@test run_response.status == 303
result_location = only(header.second for header in run_response.headers if header.first == "Location")
run_id = basename(result_location)
wait(Sparlectra._POWERFLOW_WEBUI_JOBS[run_id]["task"])
result = get_powerflow_result(run_id)
@test result["success"]
@test result["output_dir"] == joinpath(abspath(output_root), run_id)
@test !ispath(joinpath(tmpdir, "outside-runs"))
      run_log_text = read(joinpath(result["output_dir"], "run.log"), String)
      @test occursin("Wall time", run_log_text)
      @test !occursin("Unknown Sparlectra configuration key: power_flow.qlimits.enforcement_mode", run_log_text)
      @test occursin("Wall time   :", run_log_text)
      @test occursin("Output time :", run_log_text)
      @test occursin("Solver time :", run_log_text)
      @test !occursin("Solver time    :unavailable", run_log_text)
      @test !occursin("Representative wall time", run_log_text)
      @test !occursin("representative_time:", run_log_text)
      @test !occursin("solver_time:          n/a", run_log_text)
      result_json_text = read(joinpath(result["output_dir"], "result.json"), String)
      result_json = Sparlectra._parse_service_json(result_json_text)
      @test result_json["run_id"] == run_id
      @test result_json["schema_version"] == "1.0"
      @test haskey(result_json, "service_phase_timings")
      @test !haskey(result_json, "wall_time")
      @test !haskey(result_json, "representative_time")
      operation_log_text = read(operation_log_path, String)
      @test occursin("\"event\":\"powerflow_submitted\"", operation_log_text)
      @test occursin("\"event\":\"powerflow_started\"", operation_log_text)
      @test occursin("\"event\":\"powerflow_completed\"", operation_log_text)
      @test occursin("\"matpower_auto_profile\":\"recommend\"", operation_log_text)
      @test occursin("\"matpower_ratio\":\"normal\"", operation_log_text)
      @test occursin("\"matpower_shift_sign\":1.0", operation_log_text)
      @test occursin("\"matpower_shift_unit\":\"deg\"", operation_log_text)
      @test occursin("\"sparlectra_version\":\"$(Sparlectra.version())\"", operation_log_text)
      @test any(line -> occursin(r"\"timestamp\":\"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}\.\d{3}Z\"", line), eachline(IOBuffer(operation_log_text)))
      @test occursin("\"event\":\"diagnostics_enabled\"", operation_log_text)
      @test occursin("\"event\":\"performance_timing_enabled\"", operation_log_text)
      @test occursin("\"event\":\"detailed_result_csv_export_enabled\"", operation_log_text)
      @test occursin("\"csv_format\":\"excel_de\"", operation_log_text)
      @test occursin("\"delimiter\":\";\"", operation_log_text)
      @test occursin("\"decimal_separator\":\",\"", operation_log_text)
      @test occursin("\"thousands_separator\":\".\"", operation_log_text)
      @test occursin("\"event\":\"detailed_result_csv_exported\"", operation_log_text)
      @test occursin("\"artifacts\":[\"bus_voltages_complex.csv\",\"branch_flows.csv\"]", operation_log_text)
      performance_log_text = read(joinpath(result["output_dir"], "performance.log"), String)
      q_limit_log_text = read(joinpath(result["output_dir"], "q_limit.log"), String)
      for artifact_text in (run_log_text, result_json_text, operation_log_text, performance_log_text, q_limit_log_text)
        _assert_no_webui_transcript_markers(artifact_text)
      end
      manual_import_form = copy(form)
      manual_import_form["matpower_import_auto_profile"] = "off"
      manual_import_form["matpower_import_ratio"] = "reciprocal"
      manual_import_form["matpower_import_shift_sign"] = "-1.0"
      manual_import_form["matpower_import_shift_unit"] = "rad"
      manual_response = Sparlectra.route_sparlectra_webui("POST", "/powerflow/run", manual_import_form; output_root)
      @test manual_response.status == 303
      manual_location = only(header.second for header in manual_response.headers if header.first == "Location")
      manual_run_id = basename(manual_location)
      wait(Sparlectra._POWERFLOW_WEBUI_JOBS[manual_run_id]["task"])
      manual_result = get_powerflow_result(manual_run_id)
      @test manual_result["success"]
      manual_effective_cfg = Sparlectra.load_sparlectra_config(joinpath(manual_result["output_dir"], "effective_config.yaml"); reload = true)
      @test manual_effective_cfg.matpower.auto_profile === :off
      @test manual_effective_cfg.matpower.ratio === :reciprocal
      @test manual_effective_cfg.matpower.shift_sign == -1.0
      @test manual_effective_cfg.matpower.shift_unit === :rad
      manual_run_log = read(joinpath(manual_result["output_dir"], "run.log"), String)
      @test occursin("Final effective MATPOWER import options", manual_run_log)
      @test occursin("matpower_import.ratio: reciprocal", manual_run_log)
      @test !isempty(run_id)
      result_response = Sparlectra.handle_powerflow_result(run_id)
      result_html = String(result_response.body)
      @test result_response.status == 200
      for field in ("run_id", "status", "iterations", "final_mismatch")
        @test occursin(field, result_html)
      end
      @test occursin(run_id, result_html)
      @test !occursin("action=\"/powerflow/abort/$(run_id)\"", result_html)
      failed_result_html = Sparlectra.render_powerflow_result(Dict("run_id" => "failed-run", "status" => "failed", "success" => false))
      @test !occursin("action=\"/powerflow/abort/failed-run\"", failed_result_html)
      converged_result_html = Sparlectra.render_powerflow_result(Dict(
        "run_id" => "converged-run",
        "status" => "success",
        "success" => true,
        "converged" => true,
        "numerical_converged" => true,
        "solution_available" => true,
        "iterations" => 6,
        "final_mismatch" => 6.3e-11,
        "reason" => "converged",
      ))
      @test !occursin("<span class=\"summary-label\">Solver status</span>", converged_result_html)
      @test !occursin("<span class=\"summary-label\">Iterations</span>", converged_result_html)
      @test !occursin("<span class=\"summary-label\">Final mismatch</span>", converged_result_html)
      @test occursin("<span class=\"summary-label\">Run status</span>", converged_result_html)
      @test occursin("<span class=\"summary-label\">Elapsed time</span>", converged_result_html)
      @test occursin("converged", converged_result_html)
      @test occursin("numerical_converged</th><td>true</td>", converged_result_html)
      @test occursin("solution_available</th><td>true</td>", converged_result_html)
      @test occursin("iterations</th><td>6</td>", converged_result_html)
      @test occursin("final_mismatch</th><td>6.3e-11</td>", converged_result_html)
      nonconverged_result_html = Sparlectra.render_powerflow_result(Dict(
        "run_id" => "nonconverged-run",
        "status" => "not_converged",
        "success" => false,
        "service_status" => "completed",
        "numerical_status" => "not_converged",
        "converged" => false,
        "numerical_converged" => false,
        "solution_available" => false,
        "iterations" => 80,
        "reason" => "nr_mismatch_not_converged",
        "final_outcome" => Dict("converged" => false),
      ))
      @test occursin("status-badge status-error", nonconverged_result_html)
      @test occursin("not_converged", nonconverged_result_html)
      @test occursin("converged</th><td>false</td>", nonconverged_result_html)
      @test occursin("numerical_converged</th><td>false</td>", nonconverged_result_html)
      @test occursin("solution_available</th><td>false</td>", nonconverged_result_html)
      @test occursin("iterations</th><td>80</td>", nonconverged_result_html)
      @test occursin("reason</th><td>nr_mismatch_not_converged</td>", nonconverged_result_html)
      missing_metric_html = Sparlectra.render_powerflow_result(Dict("run_id" => "missing-metrics", "status" => "success", "success" => true))
      for field in ("converged", "numerical_converged", "solution_available", "iterations", "final_mismatch", "reason")
        @test occursin("$(field)</th><td>n/a</td>", missing_metric_html)
      end
      for active_status in ("queued", "running")
        active_status_html = Sparlectra.render_powerflow_result(Dict("run_id" => "$(active_status)-run", "status" => active_status))
        @test occursin("http-equiv=\"refresh\"", active_status_html)
        @test occursin("action=\"/powerflow/abort/$(active_status)-run\"", active_status_html)
      end
      for terminal_status in ("success", "failed", "aborted")
        terminal_status_html = Sparlectra.render_powerflow_result(Dict("run_id" => "$(terminal_status)-run", "status" => terminal_status, "elapsed_seconds" => 83))
        @test !occursin("http-equiv=\"refresh\"", terminal_status_html)
        @test occursin("class=\"runtime-card\"", terminal_status_html)
        @test occursin("Elapsed time", terminal_status_html)
        @test occursin("00:01:23.000", terminal_status_html)
      end

      artifacts = list_powerflow_artifacts(run_id)
      artifact_names = Set(artifact["name"] for artifact in artifacts)
      @test Set(("result.json", "run.log", "effective_config.yaml", "run_metadata.yaml", "performance.log", "diagnose.log", "q_limit.log", "bus_voltages_complex.csv", "branch_flows.csv")) ⊆ artifact_names
      artifact_response = Sparlectra.handle_powerflow_artifacts(run_id)
      artifact_html = String(artifact_response.body)
      for name in ("result.json", "run.log", "effective_config.yaml", "run_metadata.yaml", "performance.log", "diagnose.log", "q_limit.log", "bus_voltages_complex.csv", "branch_flows.csv")
        @test occursin(name, artifact_html)
      end
      qlimit_artifact = only(artifact for artifact in artifacts if artifact["name"] == "q_limit.log")
      @test qlimit_artifact["kind"] == "q_limit_log"
      @test qlimit_artifact["mime_type"] == "text/plain"
      qlimit_download = Sparlectra.handle_powerflow_artifact_download(run_id, "q_limit.log")
      @test qlimit_download.status == 200
      @test occursin("Resolved Q-limit options", String(qlimit_download.body))
      @test Sparlectra.route_sparlectra_webui("GET", "/powerflow/artifact/$(run_id)/q_limit.log"; output_root).status == 200
      for name in ("bus_voltages_complex.csv", "branch_flows.csv")
        csv_artifact = only(artifact for artifact in artifacts if artifact["name"] == name)
        @test csv_artifact["kind"] == "csv"
        @test csv_artifact["mime_type"] == "text/csv"
        csv_download = Sparlectra.handle_powerflow_artifact_download(run_id, name)
        @test csv_download.status == 200
        @test ("Content-Type" => "text/csv") in csv_download.headers
        expected_header = name == "bus_voltages_complex.csv" ? "bus;bus_name;type;vm_pu;va_deg" : "branch;branch_index;from_bus;to_bus;status"
        @test startswith(String(csv_download.body), expected_header)
      end
      legacy_diagnostic_path = joinpath(result["output_dir"], "diagnose.txt")
      write(legacy_diagnostic_path, "legacy diagnostics\n")
      legacy_artifacts = list_powerflow_artifacts(run_id)
      @test any(artifact -> artifact["name"] == "diagnose.txt", legacy_artifacts)
      legacy_download = Sparlectra.handle_powerflow_artifact_download(run_id, "diagnose.txt")
      @test legacy_download.status == 200
      @test String(legacy_download.body) == "legacy diagnostics\n"
      large_artifact_path = joinpath(result["output_dir"], "large_preview.log")
      large_artifact_text = repeat("0123456789abcdef", 5000)
      write(large_artifact_path, large_artifact_text)
      large_preview = Sparlectra.handle_powerflow_artifact(run_id, "large_preview.log")
      @test large_preview.status == 200
      large_preview_html = String(large_preview.body)
      @test occursin("Preview truncated", large_preview_html)
      @test !occursin(large_artifact_text, large_preview_html)
      large_download = Sparlectra.handle_powerflow_artifact_download(run_id, "large_preview.log")
      @test String(large_download.body) == large_artifact_text
      result_artifact = Sparlectra.handle_powerflow_artifact(run_id, "result.json")
      @test result_artifact.status == 200
      result_artifact_html = String(result_artifact.body)
      @test occursin("Artifact: result.json", result_artifact_html)
      @test occursin("<main class=\"page artifact-page\">", result_artifact_html)
      @test occursin("class=\"artifact-text-page\"", result_artifact_html)
      @test occursin("class=\"button\" href=\"?download=1\"", result_artifact_html)
      @test occursin("class=\"artifact-text\"", result_artifact_html)

      download = Sparlectra.handle_powerflow_artifact_download(run_id, "result.json")
      @test download.status == 200
      @test any(header -> header.first == "Content-Disposition", download.headers)
      @test Sparlectra.route_sparlectra_webui("GET", "/powerflow/artifact/$(run_id)/result.json"; output_root).status == 200
      @test Sparlectra.route_sparlectra_webui("GET", "/powerflow/artifact/$(run_id)/result.json?download=1"; output_root).status == 200
      zip_download = Sparlectra.handle_powerflow_artifacts_zip(run_id)
      @test zip_download.status == 200
      @test ("Content-Type" => "application/zip") in zip_download.headers
      @test zip_download.body[1:4] == UInt8[0x50, 0x4b, 0x03, 0x04]
      @test occursin("result.json", String(zip_download.body))
      @test !occursin("Project.toml", String(zip_download.body))
      @test Sparlectra.route_sparlectra_webui("GET", "/powerflow/artifact-zip/$(run_id)"; output_root).status == 200
      @test Sparlectra.handle_powerflow_artifacts_zip("unknown-run").status == 404
      operation_log_text = read(operation_log_path, String)
      @test occursin("\"event\":\"artifact_opened\"", operation_log_text)
      @test occursin("\"event\":\"artifact_downloaded\"", operation_log_text)

      for unsafe_name in ("../result.json", "../../Project.toml", "/etc/passwd", raw"C:\Users\x\file.txt")
        unsafe = Sparlectra.handle_powerflow_artifact(run_id, unsafe_name)
        @test unsafe.status == 400
        @test occursin("Unsafe artifact name rejected", String(unsafe.body))
      end

      lock(Sparlectra._POWERFLOW_SERVICE_LOCK) do
        empty!(Sparlectra._POWERFLOW_SERVICE_RUNS)
      end
      refreshed = Sparlectra.handle_powerflow_refresh(output_root)
      @test run_id in refreshed["loaded_runs"]
      @test Sparlectra.route_sparlectra_webui("POST", "/powerflow/refresh"; output_root).status == 303
      @test occursin("\"event\":\"history_refreshed\"", read(operation_log_path, String))
      history_response = Sparlectra.handle_powerflow_history(output_root)
      @test history_response.status == 200
      history_html = String(history_response.body)
      @test occursin(run_id, history_html)
      @test occursin("<th>Date/Time</th>", history_html)
      @test occursin("class=\"status-badge status-success\"", history_html)
      @test occursin("action=\"/powerflow/delete/$(run_id)\"", history_html)
      @test occursin("action=\"/powerflow/delete_all\"", history_html)
      @test occursin("action=\"/powerflow/refresh\"", history_html)
      @test occursin(r"\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}", history_html)
      @test Sparlectra.route_sparlectra_webui("GET", "/powerflow/abort/$(run_id)"; output_root).status == 404
      @test Sparlectra.route_sparlectra_webui("POST", "/powerflow/abort/unknown-run"; output_root).status == 404
      @test Sparlectra.route_sparlectra_webui("POST", "/powerflow/abort/..%2Foutside"; output_root).status == 404
      @test Sparlectra.route_sparlectra_webui("GET", "/powerflow/delete/$(run_id)"; output_root).status == 404
      unsafe_delete = Sparlectra.route_sparlectra_webui("POST", "/powerflow/delete/..%2Foutside"; output_root)
      @test unsafe_delete.status == 400
      @test occursin("Unsafe PowerFlow run ID rejected", String(unsafe_delete.body))
      @test Sparlectra.route_sparlectra_webui("GET", "/webui/shutdown"; output_root).status == 404
      @test Sparlectra.route_sparlectra_webui("POST", "/webui/shutdown"; output_root).status == 503
      @test Sparlectra.route_sparlectra_webui("POST", "/webui/heartbeat"; output_root).status == 204
      disposable_id = "disposable-run"
      disposable_dir = joinpath(output_root, disposable_id)
      mkpath(disposable_dir)
      disposable_result = joinpath(disposable_dir, "result.json")
      write(disposable_result, "{}")
      disposable_entries = load_powerflow_run_index(output_root)["runs"]
      push!(disposable_entries, Dict(
        "run_id" => disposable_id,
        "output_dir" => disposable_dir,
        "result_file" => disposable_result,
      ))
      Sparlectra._write_powerflow_run_entries!(output_root, disposable_entries)
      deleted_response = Sparlectra.route_sparlectra_webui("POST", "/powerflow/delete/$(disposable_id)"; output_root)
      @test deleted_response.status == 303
      @test !ispath(disposable_dir)
      @test all(entry -> get(entry, "run_id", "") != disposable_id, load_powerflow_run_index(output_root)["runs"])
      @test occursin("\"event\":\"run_deleted\"", read(operation_log_path, String))

      entered_webui = Channel{Nothing}(1)
      release_webui = Channel{Nothing}(1)
      slow_runner = function(request; case_directory = nothing)
        request["phase_callback"]("linear_solve")
        put!(entered_webui, nothing)
        take!(release_webui)
        return start_powerflow_run(request; case_directory)
      end
      active_request = Dict(
        "casefile" => casefile,
        "config_file" => config_file,
        "output_root" => output_root,
      )
      lifecycle_event = (event; fields...) -> Sparlectra.record_webui_operation!(output_root, event; route = "/powerflow/run", method = "POST", user_action = false, fields...)
      active = Sparlectra.start_webui_powerflow_run(active_request; runner = slow_runner, event_callback = lifecycle_event)
      take!(entered_webui)
      active_id = active["run_id"]
      active_result_html = String(Sparlectra.route_sparlectra_webui("GET", "/powerflow/result/$(active_id)"; output_root).body)
      active_form_html = String(Sparlectra.route_sparlectra_webui("GET", "/powerflow"; output_root).body)
      active_history_html = String(Sparlectra.route_sparlectra_webui("GET", "/powerflow/history"; output_root).body)
      abort_action = "action=\"/powerflow/abort/$(active_id)\""
      @test occursin("<form method=\"post\" $(abort_action)", active_result_html)
      @test occursin("class=\"runtime-card\"", active_result_html)
      @test occursin("Elapsed time", active_result_html)
      @test occursin("class=\"status-badge status-running\">running</span>", active_result_html)
      @test occursin("http-equiv=\"refresh\"", active_result_html)
      @test occursin("content=\"$(Sparlectra.WEBUI_STATUS_AUTO_REFRESH_SECONDS); url=/powerflow/result/$(active_id)?autorefresh=1\"", active_result_html)
      @test occursin("This page refreshes automatically while the run is active.", active_result_html)
      @test occursin("<td>linear_solve</td>", active_result_html)
      @test occursin(">Refresh status</a>", active_result_html)
      @test occursin("PowerFlow run is running:", active_form_html)
      @test occursin(abort_action, active_form_html)
      @test occursin(abort_action, active_history_html)
      active_log = read(operation_log_path, String)
      status_position = findlast("\"event\":\"powerflow_status_opened\"", active_log)
      started_position = findlast("\"event\":\"powerflow_started\"", active_log)
      @test started_position < status_position
      auto_refresh_count = count(==('\n'), active_log)
      auto_refresh_html = String(Sparlectra.route_sparlectra_webui("GET", "/powerflow/result/$(active_id)?autorefresh=1"; output_root).body)
      @test occursin("http-equiv=\"refresh\"", auto_refresh_html)
      @test count(==('\n'), read(operation_log_path, String)) == auto_refresh_count
      @test !occursin("\"event\":\"powerflow_aborted\",\"run_id\":\"$(active_id)\"", active_log)
      aborted_response = Sparlectra.route_sparlectra_webui("POST", "/powerflow/abort/$(active_id)"; output_root)
      @test aborted_response.status == 303
      aborted_html = String(Sparlectra.handle_powerflow_result(active_id).body)
      @test occursin("aborting", aborted_html)
      @test occursin("class=\"runtime-card\"", aborted_html)
      @test occursin("http-equiv=\"refresh\"", aborted_html)
      @test occursin("Aborting requested. Current phase:", aborted_html)
      @test occursin("This phase may need to finish before cancellation is observed.", aborted_html)
      @test !occursin(abort_action, aborted_html)
      @test occursin("\"event\":\"powerflow_abort_requested\"", read(operation_log_path, String))
      @test occursin("\"current_phase\":\"linear_solve\"", read(operation_log_path, String))
      repeated_abort = Sparlectra.route_sparlectra_webui("POST", "/powerflow/abort/$(active_id)"; output_root)
      @test repeated_abort.status == 303
      @test occursin("\"event\":\"powerflow_abort_already_requested\"", read(operation_log_path, String))
      @test occursin("\"event\":\"powerflow_phase_started\"", read(operation_log_path, String))
      @test !any(
        line -> occursin("\"event\":\"powerflow_phase_started\"", line) && occursin("\"phase\":\"linear_solve\"", line),
        eachline(operation_log_path),
      )
      active_delete = Sparlectra.route_sparlectra_webui("POST", "/powerflow/delete/$(active_id)"; output_root)
      @test active_delete.status == 409
      @test occursin("This run is still active", String(active_delete.body))
      @test occursin("\"event\":\"run_delete_rejected\"", read(operation_log_path, String))
      put!(release_webui, nothing)
      wait(Sparlectra._POWERFLOW_WEBUI_JOBS[active_id]["task"])
      terminal_html = String(Sparlectra.handle_powerflow_result(active_id).body)
      @test occursin("<td>aborted</td>", terminal_html)
      @test occursin("class=\"runtime-card\"", terminal_html)
      @test occursin("Run aborted by user.", terminal_html)
      @test !occursin(abort_action, terminal_html)
      @test !occursin("http-equiv=\"refresh\"", terminal_html)
      terminal_log = read(operation_log_path, String)
      @test any(line -> occursin("\"event\":\"powerflow_aborted\"", line) && occursin("\"run_id\":\"$(active_id)\"", line), eachline(IOBuffer(terminal_log)))
      @test Sparlectra.get_active_webui_powerflow_job() === nothing
      terminal_delete = Sparlectra.route_sparlectra_webui("POST", "/powerflow/delete/$(active_id)"; output_root)
      @test terminal_delete.status == 303

      hard_reset_job = Dict{String,Any}(
        "run_id" => "hard-reset-test",
        "status" => "aborting",
        "casefile" => casefile,
        "config_file" => config_file,
        "output_root" => output_root,
        "output_dir" => joinpath(output_root, "hard-reset-test"),
        "started_at" => Dates.now(Dates.UTC) - Dates.Second(120),
        "finished_at" => nothing,
        "message" => "Abort requested.",
        "abort_requested" => Threads.Atomic{Bool}(true),
        "abort_requested_at" => Dates.now(Dates.UTC) - Dates.Second(Sparlectra.WEBUI_ABORT_HARD_RESET_AFTER_SECONDS + 1),
        "current_phase" => "linear_solve",
        "phase_started_at" => Dates.now(Dates.UTC),
        "last_progress_at" => Dates.now(Dates.UTC),
      )
      Sparlectra._POWERFLOW_WEBUI_JOBS["hard-reset-test"] = hard_reset_job
      hard_reset_html = String(Sparlectra.handle_powerflow_result("hard-reset-test").body)
      @test occursin("Hard reset Web UI", hard_reset_html)
      hard_reset_response = Sparlectra.route_sparlectra_webui("POST", "/powerflow/hard-reset/hard-reset-test"; output_root)
      @test hard_reset_response.status == 200
      @test Sparlectra.get_webui_powerflow_job("hard-reset-test")["status"] == "aborted_unknown"
      @test occursin("\"event\":\"webui_hard_reset_requested\"", read(operation_log_path, String))
      @test !Sparlectra.get_webui_powerflow_job("hard-reset-test")["success"]

      operation_log_page = Sparlectra.route_sparlectra_webui("GET", "/webui/operation-log"; output_root)
      @test operation_log_page.status == 200
      @test occursin("Download operation log", String(operation_log_page.body))
      operation_log_download = Sparlectra.route_sparlectra_webui("GET", "/webui/operation-log/download"; output_root)
      @test operation_log_download.status == 200
      @test any(header -> header.first == "Content-Disposition", operation_log_download.headers)
      static_count_before = count(==('\n'), read(operation_log_path, String))
      @test Sparlectra.route_sparlectra_webui("GET", "/static/sparlectra.css"; output_root).status == 200
      @test count(==('\n'), read(operation_log_path, String)) == static_count_before
      compact_root = joinpath(tmpdir, "compact-log")
      compact_log = Sparlectra.webui_operation_log_path(compact_root)
      mkpath(compact_root)
      open(compact_log, "w") do io
        println(io)
        println(io, "malformed")
        for index in 1:11
          println(io, "{\"event\":\"entry-$(index)\"}")
        end
      end
      @test Sparlectra._compact_webui_operation_log!(compact_log; max_entries = 10, keep_entries = 3)
      compact_lines = readlines(compact_log)
      @test length(compact_lines) == 3
      @test all(line -> Sparlectra._parse_service_json(line) isa AbstractDict, compact_lines)
      @test occursin("\"event\":\"entry-9\"", compact_lines[1])
      @test occursin("\"event\":\"entry-11\"", compact_lines[3])
      @test !Sparlectra._compact_webui_operation_log!(compact_log; max_entries = 10, keep_entries = 3)
      @test !Sparlectra._compact_webui_operation_log!(joinpath(tmpdir, "missing-operation-log.jsonl"); max_entries = 10, keep_entries = 3)
      compact_viewer = Sparlectra.handle_webui_operation_log(compact_log)
      compact_download = Sparlectra.handle_webui_operation_log(compact_log; download = true)
      @test compact_viewer.status == 200
      @test compact_download.status == 200
      @test occursin("entry-11", String(compact_viewer.body))
      @test occursin("entry-11", String(compact_download.body))
      logging_failure_root = joinpath(tmpdir, "logging-is-a-file")
      write(logging_failure_root, "not a directory")
      record_result = @test_logs (:warn, r"Could not record Web UI operation") begin
        Sparlectra.record_webui_operation!(logging_failure_root, "expected_failure")
      end
      @test record_result == false
      route_result = @test_logs (:warn, r"Could not record Web UI operation") begin
        Sparlectra.route_sparlectra_webui("GET", "/powerflow"; output_root = logging_failure_root)
      end
      @test route_result.status == 200
      empty_delete_root = joinpath(tmpdir, "empty-delete-root")
      @test Sparlectra.route_sparlectra_webui("POST", "/powerflow/delete_all"; output_root = empty_delete_root).status == 303
      @test occursin("\"event\":\"all_runs_deleted\"", read(Sparlectra.webui_operation_log_path(empty_delete_root), String))

      @test occursin(run_id, String(Sparlectra.handle_powerflow_result(run_id).body))
      shared_pages = (
        form_html,
        result_html,
        artifact_html,
        result_artifact_html,
        history_html,
      )
      for page_html in shared_pages
        @test occursin("src=\"/assets/logo.png\"", page_html)
        @test occursin("Sparlectra.jl v$(Sparlectra.version())", page_html)
      end

      source = join(read(joinpath(@__DIR__, "..", "src", "webui", file), String) for file in ("docs.jl", "forms.jl", "handlers.jl", "routes.jl", "views.jl", "webui.jl"))
      forbidden = ["run" * "pf!", "run" * "pf_rectangular!", "run_complex_nr_" * "rectangular"]
      @test all(name -> !occursin(name, source), forbidden)
      @test !occursin("Voltage-magnitude/angle blend strategy.", source)

      lock(Sparlectra._POWERFLOW_SERVICE_LOCK) do
        empty!(Sparlectra._POWERFLOW_SERVICE_RUNS)
      end
      probe = listen(ip"127.0.0.1", UInt16(0))
      port = Int(getsockname(probe)[2])
      close(probe)
      lifecycle_io = IOBuffer()
      server = redirect_stderr(devnull) do
        start_sparlectra_webui(port = port, output_root = output_root, auto_shutdown_on_browser_close = false, _lifecycle_io = lifecycle_io)
      end
      startup_output = String(take!(lifecycle_io))
      _assert_no_webui_transcript_markers(startup_output)
      _assert_no_webui_internal_lifecycle_markers(startup_output)
      @test occursin("Sparlectra Web UI is available at http://127.0.0.1:$(port)/powerflow", startup_output)
      @test occursin("Stop: use Stop Web UI in the browser, close(server), or Ctrl+C here.", startup_output)
      @test occursin("Operation log: ", startup_output)
      @test !istaskdone(server.task)
      @test isdir(output_root)
      @test server.runtime.case_directory == Sparlectra.default_webui_case_cache_dir(output_root)
      @test isdir(server.runtime.case_directory)
      @test isfile(server.runtime.config_file)
      @test isfile(joinpath(server.runtime.case_directory, "warmup_case3.jl"))
      @test isfile(server.runtime.operation_log)
      @test !startswith(server.runtime.case_directory, normpath(pkgdir(Sparlectra)))
      @test !startswith(server.runtime.config_file, normpath(pkgdir(Sparlectra)))
      startup_log = read(server.runtime.operation_log, String)
      @test occursin("\"event\":\"webui_start_requested\"", startup_log)
      @test occursin("\"event\":\"webui_config_loaded\"", startup_log)
      @test occursin("\"event\":\"webui_routes_registered\"", startup_log)
      @test occursin("\"event\":\"webui_start\"", startup_log)
      @test occursin("\"event\":\"webui_server_bound\"", startup_log)
      @test occursin("\"event\":\"webui_started\"", startup_log)
      @test !occursin("browser_window_close", startup_log)
      @test occursin("\"sparlectra_version\":\"$(Sparlectra.version())\"", startup_log)
      @test occursin("\"sparlectra_package_path\":", startup_log)
      @test occursin(basename(Sparlectra._sparlectra_package_path()), startup_log)
      @test occursin("\"sparlectra_git_commit\":", startup_log)
      for field in ("output_root", "config_file", "case_cache_dir", "operation_log")
        @test occursin("\"$(field)\":", startup_log)
      end
      @test get_powerflow_result(run_id)["run_id"] == run_id
      try
        @test server.url == "http://127.0.0.1:$(port)/powerflow"
        response_text = _webui_http_request(port, "GET", "/powerflow")
        @test occursin("HTTP/1.1 200 OK", response_text)
        @test occursin("PowerFlow run", response_text)
        @test occursin("Output root", response_text)
        @test occursin(abspath(output_root), response_text)
        @test occursin("Config file", response_text)
        @test occursin(server.runtime.config_file, response_text)
        @test occursin("MATPOWER case cache", response_text)
        @test occursin(server.runtime.case_directory, response_text)
        @test occursin("Operation log", response_text)
        @test occursin(server.runtime.operation_log, response_text)
        @test occursin("MATPOWER import conventions", response_text)
        @test occursin("name=\"matpower_import_auto_profile\"", response_text)
        @test occursin("name=\"matpower_import_ratio\"", response_text)
        @test occursin("name=\"matpower_import_shift_sign\"", response_text)
        @test occursin("name=\"matpower_import_shift_unit\"", response_text)

        logo_response_text = _webui_http_request(port, "GET", "/assets/logo.png")
        @test occursin("HTTP/1.1 200 OK", logo_response_text)
        @test occursin("Content-Type: image/png", logo_response_text)

        heartbeat_response = _webui_http_request(port, "POST", "/webui/heartbeat")
        @test occursin("HTTP/1.1 204", heartbeat_response)
        shutdown_response = _webui_http_request(port, "POST", "/webui/shutdown")
        @test timedwait(() -> istaskdone(server.task), 2.0) == :ok
        shutdown_output = String(take!(lifecycle_io))
        @test occursin("HTTP/1.1 200 OK", shutdown_response)
        @test occursin("Web UI stopped", shutdown_response)
        @test occursin("Sparlectra Web UI stopped by explicit shutdown.", shutdown_output)
        @test !isopen(server.listener)
      finally
        close(server)
      end

      for (name, config_text, expected_marker) in (
        ("new qlimit enforcement mode", "power_flow:\n  qlimits:\n    enforcement_mode: classic_simultaneous\n", "PowerFlow run"),
        ("legacy qlimit enforcement alias", "power_flow:\n  qlimits:\n    enforcement_mode: matpower_one_at_a_time\n", "PowerFlow run"),
        ("removed start-voltage alias", "power_flow:\n  start_mode:\n    voltage_mode: bus_vm_va_blend\n", "Configuration error"),
        ("malformed configuration", "power_flow:\n  qlimits:\n    enforcement_mode: definitely_not_supported\n", "Configuration error"),
      )
        probe = listen(ip"127.0.0.1", UInt16(0))
        config_port = Int(getsockname(probe)[2])
        close(probe)
        config_root = mktempdir()
        config_path = joinpath(config_root, "configuration.yaml")
        write(config_path, config_text)
        config_io = IOBuffer()
        config_server = start_sparlectra_webui(
          port = config_port,
          output_root = joinpath(config_root, "runs"),
          config_file = config_path,
          auto_shutdown_on_browser_close = false,
          _lifecycle_io = config_io,
        )
        try
          config_response = _webui_http_request(config_port, "GET", "/powerflow")
          @test occursin("HTTP/1.1 200 OK", config_response)
          @test !isempty(config_response)
          @test occursin("Sparlectra", config_response)
          @test occursin(expected_marker, config_response)
          config_startup_output = String(take!(config_io))
          @test occursin("Sparlectra Web UI is available at http://127.0.0.1:$(config_port)/powerflow", config_startup_output)
          _assert_no_webui_internal_lifecycle_markers(config_startup_output)
          config_log = read(config_server.runtime.operation_log, String)
          @test occursin("\"event\":\"webui_config_loaded\"", config_log)
          name == "malformed configuration" && @test occursin("error_visible", config_log)
        finally
          close(config_server)
          wait(config_server.task)
        end
      end

      close_io = IOBuffer()
      restarted = start_sparlectra_webui(port = port, output_root = output_root, auto_shutdown_on_browser_close = false, _lifecycle_io = close_io)
      take!(close_io)
      close(restarted)
      @test timedwait(() -> istaskdone(restarted.task), 2.0) == :ok
      close_output = String(take!(close_io))
      @test occursin("Sparlectra Web UI stopped by server close.", close_output)

      interrupt_io = IOBuffer()
      interrupted = start_sparlectra_webui(port = port, output_root = output_root, auto_shutdown_on_browser_close = false, _lifecycle_io = interrupt_io)
      take!(interrupt_io)
      @test Sparlectra._wait_sparlectra_webui(interrupted; wait_for_task = _ -> throw(InterruptException()))
      @test timedwait(() -> istaskdone(interrupted.task), 2.0) == :ok
      interrupt_output = String(take!(interrupt_io))
      @test occursin("Sparlectra Web UI stopped by Ctrl+C.", interrupt_output)
      @test timedwait(() -> istaskdone(interrupted.task), 2.0) == :ok
      @test !isopen(interrupted.listener)

      browser_probe = listen(ip"127.0.0.1", UInt16(0))
      browser_port = Int(getsockname(browser_probe)[2])
      close(browser_probe)
      browser_io = IOBuffer()
      browser_server = start_sparlectra_webui(
        port = browser_port,
        output_root = output_root,
        shutdown_on_browser_close = true,
        open_browser = true,
        _lifecycle_io = browser_io,
        _browser_opener = _ -> run(`$(Base.julia_cmd()) --startup-file=no -e "sleep(0.2)"`; wait = false),
      )
      initial_browser_output = String(take!(browser_io))
      @test isopen(browser_server.listener)
      @test occursin("Sparlectra Web UI is available at http://127.0.0.1:$(browser_port)/powerflow", initial_browser_output)
      _assert_no_webui_internal_lifecycle_markers(initial_browser_output)
      sleep(0.5)
      @test !istaskdone(browser_server.task)
      browser_response = _webui_http_request(browser_port, "GET", "/powerflow")
      @test occursin("HTTP/1.1 200 OK", browser_response)
      @test !isempty(browser_response)
      @test occursin("PowerFlow run", browser_response)
      browser_log = read(browser_server.runtime.operation_log, String)
      @test occursin("\"event\":\"browser_close_monitor_skipped\"", browser_log)
      @test !occursin("browser_window_close", browser_log)
      close(browser_server)
      @test timedwait(() -> istaskdone(browser_server.task), 2.0) == :ok
      @test !isopen(browser_server.listener)

      manual_probe = listen(ip"127.0.0.1", UInt16(0))
      manual_port = Int(getsockname(manual_probe)[2])
      close(manual_probe)
      manual_io = IOBuffer()
      manual_server = start_sparlectra_webui(
        port = manual_port,
        output_root = output_root,
        open_browser = true,
        _lifecycle_io = manual_io,
        _browser_opener = _ -> nothing,
      )
      try
        manual_output = String(take!(manual_io))
        @test isopen(manual_server.listener)
        @test occursin("Sparlectra Web UI is available at http://127.0.0.1:$(manual_port)/powerflow", manual_output)
        manual_response = _webui_http_request(manual_port, "GET", "/powerflow")
        @test occursin("HTTP/1.1 200 OK", manual_response)
        @test occursin("PowerFlow run", manual_response)
      finally
        close(manual_server)
        @test timedwait(() -> istaskdone(manual_server.task), 2.0) == :ok
      end

      heartbeat_probe = listen(ip"127.0.0.1", UInt16(0))
      heartbeat_port = Int(getsockname(heartbeat_probe)[2])
      close(heartbeat_probe)
      heartbeat_io = IOBuffer()
      heartbeat_server = start_sparlectra_webui(
        port = heartbeat_port,
        output_root = output_root,
        auto_shutdown_on_browser_close = true,
        browser_heartbeat_timeout_seconds = 0.2,
        _lifecycle_io = heartbeat_io,
      )
      take!(heartbeat_io)
      sleep(0.3)
      @test isopen(heartbeat_server.listener)
      _webui_http_request(heartbeat_port, "POST", "/webui/heartbeat")
      Sparlectra._webui_begin_request!(heartbeat_server.runtime)
      sleep(0.3)
      @test isopen(heartbeat_server.listener)
      Sparlectra._webui_finish_request!(heartbeat_server.runtime)
      sleep(0.3)
      @test isopen(heartbeat_server.listener)
      close(heartbeat_server)
      @test timedwait(() -> istaskdone(heartbeat_server.task), 2.0) == :ok
      @test !isopen(heartbeat_server.listener)
    end
  end
  return nothing
end
