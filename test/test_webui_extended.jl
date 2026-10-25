# Copyright 2023-2026 Udo Schmitz
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# file: test/test_webui_extended.jl
# purpose: extended Web UI tests against a running local server: config
#          editor and refresh, case selection/import/delete, contextual
#          help, runs through the form, artifacts, the job lifecycle, the
#          server start and stop paths, and the warm-up page
using Sparlectra
using Dates
using Sockets
using Test
import ZipArchives

# _dtf_network_fixture and the form helpers, shared with test_webui.jl
isdefined(@__MODULE__, :_dtf_network_fixture) || include("test_webui_support.jl")

const WEBUI_INTERNAL_LIFECYCLE_MARKERS = (
  "webui_start_requested",
  "webui_config_loaded",
  "webui_routes_registered",
  "webui_server_bound",
  "browser_close_monitor_skipped",
)

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
    "ignore_webui_settings" => "false",
    "power_flow_mode" => "manual",
    "power_flow_tol" => "1e-8",
    "power_flow_max_iter" => "80",
    "power_flow_autodamp" => "on",
    "power_flow_autodamp_min" => "0.05",
    "power_flow_qlimits_enabled" => "on",
    "power_flow_qlimits_enforcement_mode" => "classic_simultaneous",
    "power_flow_solver" => "rectangular",
    "power_flow_apslf_order" => "40",
    "power_flow_apslf_nr_polish" => "true",
    "power_flow_apslf_start_enabled" => "false",
    "power_flow_apslf_start_order" => "40",
    "power_flow_dc_seed_unconditional" => "false",
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
    "power_flow_merit_enabled" => "true",
    "power_flow_merit_armijo_c1" => "1e-4",
    "power_flow_merit_fallback_max_mismatch" => "true",
    "power_flow_trust_region_enabled" => "false",
    "power_flow_trust_region_initial_radius" => "1.0",
    "power_flow_trust_region_eta_accept" => "0.1",
    "power_flow_trust_region_step_mode" => "scaled",
    "power_flow_distributed_slack_enabled" => "false",
    "power_flow_distributed_slack_p_mode" => "pg_weighted",
    "power_flow_external_grid_enabled" => "false",
    "power_flow_external_grid_source" => "auto",
    "power_flow_external_grid_sk_mva" => "2000.0",
    "power_flow_external_grid_rx" => "0.1",
    "power_flow_linear_solver" => "umfpack",
    "power_flow_newton_update" => "polar",
    "power_flow_rescue" => "true",
    "power_flow_auto_slack" => "true",
    "runtime_parallel_enabled" => "true",
    "power_flow_dc_fallback" => "false",
    "cgmes_start_values" => "auto",
    "cgmes_require_boundary" => "true",
    "cgmes_infer_base_voltages" => "false",
    "cgmes_read_dynamics" => "false",
    "cgmes_equivalent_injection_as_source" => "false",
    "matpower_import_auto_profile" => "recommend",
    "matpower_import_ratio" => "normal",
    "matpower_import_shift_sign" => "1.0",
    "matpower_import_shift_unit" => "deg",
    "matpower_import_bus_shunt_model" => "admittance",
    "matpower_import_apply_bus_names" => "false",
    "matpower_import_pv_voltage_source" => "gen_vg",
    "matpower_import_compare_voltage_reference" => "imported_setpoint",
    "transformer_tap_changer_model" => "ideal",
    "matpower_export_write_solution" => "true",
    "matpower_import_dcline_mode" => "pf_injections",
    "cgmes_hvdc_mode" => "injections",
    "output_logfile_results" => "compact",
    "performance_timing" => "compact",
    "detailed_result_csv" => "on",
    "detailed_result_csv_format" => "excel_de",
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

function run_webui_extended_tests()
  @testset "Local PowerFlow Web UI" begin (function ()
    @test isdefined(SparlectraApp, :start_sparlectra_webui)
    @test isdefined(SparlectraApp, :default_webui_output_root)
    @test_throws ArgumentError start_sparlectra_webui(host = "0.0.0.0")
    let
      busy_listener = listen(ip"127.0.0.1", UInt16(0))
      busy_port = Int(getsockname(busy_listener)[2])
      try
        eaddrinuse_io = IOBuffer()
        err = nothing
        try
          redirect_stderr(devnull) do
            start_sparlectra_webui(port = busy_port, output_root = mktempdir(), _lifecycle_io = eaddrinuse_io)
          end
        catch caught
          err = caught
        end
        @test err isa ArgumentError
        @test occursin("already in use", err.msg)
        @test occursin("non-Sparlectra", err.msg)
        @test occursin("port=", err.msg)
        @test occursin("already in use", String(take!(eaddrinuse_io)))
      finally
        close(busy_listener)
      end
    end
    let
      # A RUNNING Sparlectra Web UI on the requested port is not an error:
      # the second start detects it via the /powerflow probe, opens the
      # browser on the existing instance, and returns nothing.
      port_probe = listen(ip"127.0.0.1", UInt16(0))
      running_port = Int(getsockname(port_probe)[2])
      close(port_probe)
      running_io = IOBuffer()
      running = start_sparlectra_webui(port = running_port, output_root = mktempdir(), _lifecycle_io = running_io)
      try
        opened_url = Ref("")
        second = start_sparlectra_webui(
          port = running_port,
          output_root = mktempdir(),
          open_browser = true,
          _lifecycle_io = IOBuffer(),
          _browser_opener = url -> (opened_url[] = url; nothing),
        )
        @test second === nothing
        @test opened_url[] == "http://127.0.0.1:$(running_port)/powerflow"
        # the first instance keeps serving untouched
        @test isopen(running.listener)
        @test occursin("HTTP/1.1 200 OK", _webui_http_request(running_port, "GET", "/powerflow"))
        # 0.32.0: without a sysimage start_webui.jl renders the start page
        # through the server's own routing call before it opens the window
        # (dots in the console meanwhile)
        launcher = Module(:LauncherWarmup)
        Base.include(launcher, joinpath(Sparlectra.SPARLECTRA_ROOT, "tools", "sysimage_launcher.jl"))
        SL = Base.invokelatest(getfield, launcher, :SysimageLauncher)
        sl(name, args...; kw...) = Base.invokelatest(Base.invokelatest(getfield, SL, name), args...; kw...)
        @test sl(:running_on_sysimage; current = "/opt/julia/lib/julia/sys.so") == false
        @test sl(:running_on_sysimage; current = "/home/u/.local/state/sparlectra/webui/sysimage/sparlectra.so") == true
        warm_io = IOBuffer()
        render = () -> SparlectraApp.route_sparlectra_webui("GET", "/powerflow", Dict{String,String}(); output_root = abspath(mktempdir()), runtime = running.runtime).status == 200
        @test sl(:warm_webui_page, render; io = warm_io) == true
        @test occursin(r"ready in \d+ s\.", String(take!(warm_io)))
        # a page that fails to render: the reason is printed, the window opens anyway
        @test sl(:warm_webui_page, () -> error("render broke"); io = warm_io) == false
        @test occursin("not prepared (render broke)", String(take!(warm_io)))
      finally
        close(running)
      end
    end
    default_root = default_webui_output_root()
    @test isabspath(default_root)
    @test !startswith(normpath(default_root), normpath(pkgdir(Sparlectra)))
    if Sys.islinux()
      withenv("XDG_STATE_HOME" => joinpath(mktempdir(), "state")) do
        @test default_webui_output_root() == joinpath(ENV["XDG_STATE_HOME"], "sparlectra", "webui", "runs")
      end
    end
    @test SparlectraApp._format_elapsed_duration(0.4) == "00:00:00.400"
    @test SparlectraApp._format_elapsed_duration(3) == "00:00:03.000"
    @test SparlectraApp._format_elapsed_duration(83) == "00:01:23.000"
    @test SparlectraApp._format_elapsed_duration(3725.125) == "01:02:05.125"
    @test SparlectraApp._format_elapsed_duration(nothing) == "-"


    @testset "Configuration refresh controls" begin (function ()
      root = mktempdir()
      config_path = joinpath(root, "configuration.yaml")
      write(config_path, "power_flow:\n  start_mode:\n    voltage_mode: bus_vm_va_blend\n  qlimits:\n    enabled: true\n")
      runtime = SparlectraApp._SparlectraWebUIRuntime(nothing, root, config_path, SparlectraApp.webui_operation_log_path(root), nothing, SparlectraApp.start_powerflow_run, false, false, time(), 0, nothing, IOBuffer(), ReentrantLock())
      before_page_text = read(config_path, String)
      # rendering the pages that show the maintenance actions and the
      # notice writes nothing
      SparlectraApp.route_sparlectra_webui("GET", "/powerflow/settings"; output_root = root, runtime)
      SparlectraApp.route_sparlectra_webui("GET", "/powerflow"; output_root = root, runtime)
      @test read(config_path, String) == before_page_text

      check = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/config/check", Dict("config_file" => config_path); output_root = root, runtime)
      check_body = String(check.body)
      @test check.status == 200
      @test occursin("power_flow.qlimits.enforcement_mode", check_body)
      @test occursin("power_flow.start_mode.voltage_mode", check_body)
      @test !occursin(".bak-", check_body)
      @test read(config_path, String) == before_page_text

      refresh = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/config/refresh", Dict("config_file" => config_path); output_root = root, runtime)
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
      duplicate = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/config/refresh", Dict("config_file" => duplicate_path); output_root = root, runtime)
      @test duplicate.status == 400
      @test occursin("Duplicate YAML key detected", String(duplicate.body))
      @test read(duplicate_path, String) == duplicate_before
      @test isempty(filter(name -> occursin(r"duplicate.yaml\.bak-", name), readdir(root)))

      log = read(runtime.operation_log, String)
      @test occursin("config_refresh_check_started", log)
      @test occursin("config_refresh_check_completed", log)
      @test occursin("config_refresh_write_started", log)
      @test occursin("config_refresh_write_completed", log)
    end)() end

    @testset "Configuration editor validates before atomic save" begin (function ()
      root = mktempdir()
      config_path = joinpath(root, "configuration.yaml")
      operation_log = SparlectraApp.webui_operation_log_path(root)
      original_text = "matpower_import:\n  matpower_dcline_mode: pf_injections\npower_flow:\n  tol: 1.0e-8\n  autodamp: false\n"
      write(config_path, original_text)

      temp_leftovers() = filter(name -> !(name in ("configuration.yaml", "configuration.yaml.user-keys.txt", "webui_operations.jsonl")) && !occursin(r"configuration\.yaml\.bak-", name), readdir(root))

      for mode in ("ignore_inactive", "pf_injections")
        response = SparlectraApp.handle_powerflow_config_editor_save(
          Dict("config_file" => config_path, "config_text" => "matpower_import:\n  matpower_dcline_mode: $(mode)\n");
          operation_log,
        )
        @test response.status == 200
        @test load_sparlectra_config(config_path; reload = true).matpower.matpower_dcline_mode == Symbol(mode)
        @test isempty(temp_leftovers())
      end

      nested_text = "power_flow:\n  start_current_iteration:\n    enabled: true\n    damping: 0.8\nmatpower_import:\n  matpower_dcline_mode: pf_injections\n"
      nested_response = SparlectraApp.handle_powerflow_config_editor_save(Dict("config_file" => config_path, "config_text" => nested_text); operation_log)
      @test nested_response.status == 200
      nested_cfg = load_sparlectra_config(config_path; reload = true)
      @test nested_cfg.powerflow.start_current_iteration.enabled === true
      @test nested_cfg.powerflow.start_current_iteration.damping == 0.8
      @test isempty(temp_leftovers())

      unchanged_before_invalid = read(config_path, String)
      invalid_enum = SparlectraApp.handle_powerflow_config_editor_save(
        Dict("config_file" => config_path, "config_text" => "matpower_import:\n  matpower_dcline_mode: ignore\n");
        operation_log,
      )
      invalid_enum_body = String(invalid_enum.body)
      @test invalid_enum.status == 400
      @test read(config_path, String) == unchanged_before_invalid
      @test occursin("Configuration could not be saved.", invalid_enum_body)
      @test occursin("matpower_import.matpower_dcline_mode", invalid_enum_body)
      @test occursin("ignore", invalid_enum_body)
      @test occursin("ignore_inactive, pf_injections", invalid_enum_body)
      @test occursin("matpower_dcline_mode: ignore", invalid_enum_body)
      @test isempty(temp_leftovers())

      invalid_enabled = SparlectraApp.handle_powerflow_config_editor_save(
        Dict("config_file" => config_path, "config_text" => "matpower_import:\n  matpower_dcline_mode: enabled\n");
        operation_log,
      )
      invalid_enabled_body = String(invalid_enabled.body)
      @test invalid_enabled.status == 400
      @test read(config_path, String) == unchanged_before_invalid
      @test occursin("matpower_import.matpower_dcline_mode", invalid_enabled_body)
      @test occursin("ignore_inactive, pf_injections", invalid_enabled_body)
      @test isempty(temp_leftovers())

      invalid_range = SparlectraApp.handle_powerflow_config_editor_save(Dict("config_file" => config_path, "config_text" => "power_flow:\n  tol: -1.0\n"); operation_log)
      @test invalid_range.status == 400
      @test read(config_path, String) == unchanged_before_invalid
      @test occursin("powerflow.tol must be positive", String(invalid_range.body))
      @test isempty(temp_leftovers())

      invalid_bool = SparlectraApp.handle_powerflow_config_editor_save(Dict("config_file" => config_path, "config_text" => "output:\n  console_summary: enabled\n"); operation_log)
      invalid_bool_body = String(invalid_bool.body)
      @test invalid_bool.status == 400
      @test read(config_path, String) == unchanged_before_invalid
      @test occursin("Cannot convert", invalid_bool_body)
      @test occursin("Bool", invalid_bool_body)
      @test isempty(temp_leftovers())

      unknown = SparlectraApp.handle_powerflow_config_editor_save(Dict("config_file" => config_path, "config_text" => "matpower_import:\n  not_a_real_option: true\n"); operation_log)
      @test unknown.status == 400
      @test read(config_path, String) == unchanged_before_invalid
      @test occursin("Unknown Sparlectra configuration key: matpower_import.not_a_real_option", String(unknown.body))
      @test isempty(temp_leftovers())

      log = read(operation_log, String)
      @test occursin("config_editor_saved", log)
      @test occursin("config_editor_save_failed", log)
    end)() end

    @testset "Case-specific settings profiles" begin (function ()
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
          "model.auto_profile" => "recommend",
          "matpower_import.ratio" => "normal",
          "matpower_import.shift_sign" => 1.0,
          "matpower_import.shift_unit" => "deg",
          "model.bus_shunt_model" => "admittance",
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
      successful = SparlectraApp._api_result(
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
      SparlectraApp._POWERFLOW_SERVICE_RUNS[run_id] = successful
      result_html = SparlectraApp.render_powerflow_result(SparlectraApp.to_dict(successful))
      @test occursin("Save settings for this case", result_html)
      @test !occursin("Save these settings anyway", result_html)

      # regression: a case-scope key WITHOUT a form control
      # (here the converted sidecar's apply_bus_names) and stored form
      # fields of other pages (the SE generator seed) must SURVIVE the
      # run-based save; replacing the file dropped them and the case's
      # measurement set stopped resolving bus names
      write(Sparlectra.case_config_path(joinpath(root, "case145.m")), "config_version: 1\nscope: case\ncase: case145.m\nmatpower_import:\n  apply_bus_names: false\nform:\n  gen_seed: 7\n")
      save_response = SparlectraApp.handle_powerflow_case_settings_save(run_id, Dict{String,String}(); output_root = root, operation_log = root)
      @test save_response.status == 200
      save_html = String(save_response.body)
      profile_path = SparlectraApp._webui_case_settings_path(root, joinpath(root, "case145.m"))
      @test isfile(profile_path)
      @test occursin("/powerflow?casefile=$(SparlectraApp._webui_urlencode(joinpath(root, "case145.m")))", save_html)
      @test !occursin("/powerflow?casefile=$(SparlectraApp._webui_urlencode(profile_path))", save_html)
      # the case configuration file: case-scope header, nested case-scope sections,
      # request-only fields in the form block, machine scope dropped
      @test basename(profile_path) == "case145.m.config.yaml"
      profile_text = read(profile_path, String)
      @test occursin("scope: case", profile_text)
      @test occursin("case: case145.m", profile_text)
      @test occursin("tol: 1.0e-7", profile_text)
      # autodamp true and active_set are the inherited values: not written
      # (L20), the case still runs with them
      @test !occursin("autodamp: true", profile_text)
      @test !occursin("enforcement_mode", profile_text)
      resolved145 = Sparlectra.resolve_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, joinpath(root, "case145.m")).config.powerflow
      @test (resolved145.max_iter, resolved145.autodamp, resolved145.qlimits.enforcement_mode) == (42, true, :active_set)
      @test !occursin("detailed_result_csv_format", profile_text)   # machine scope since 0.15.1
      @test !occursin("csv_format", profile_text)
      @test !occursin("benchmark", profile_text)   # machine scope stays out
      @test !occursin("logfile_results", profile_text)
      @test !occursin("effective_config", profile_text)
      case_cfg = Sparlectra.load_case_config(joinpath(root, "case145.m"))
      @test case_cfg["power_flow.max_iter"] == 42
      @test case_cfg["model.auto_profile"] == "recommend"
      @test !haskey(case_cfg, "benchmark.enabled")
      @test case_cfg["matpower_import.apply_bus_names"] == false
      form_defaults = SparlectraApp._webui_case_form_defaults(joinpath(root, "case145.m"), root)
      @test !haskey(form_defaults, "detailed_result_csv_format")
      @test !haskey(form_defaults, "performance_timing")   # compact is the spec default
      @test form_defaults["gen_seed"] == 7

      # An SCF case saves through the SAME file; the case file itself is not
      # rewritten.
      scf_case = joinpath(root, "warmup_casePST.scf.json")
      cp(abspath(joinpath(dirname(@__DIR__), "data", "scf", "sp_casePST.scf.json")), scf_case)
      scf_before = read(scf_case, String)
      scf_run_id = "case-settings-scf-test"
      scf_metadata = deepcopy(metadata)
      scf_result = SparlectraApp._api_result(
        run_id = scf_run_id, status = :succeeded, success = true, converged = true,
        solution_available = true, iterations = 3, final_mismatch = 1.0e-9,
        reason = "converged", message = "ok", casefile = scf_case,
        config_file = "configuration.yaml", output_dir = joinpath(root, scf_run_id),
        metadata = scf_metadata,
      )
      SparlectraApp._POWERFLOW_SERVICE_RUNS[scf_run_id] = scf_result
      scf_save = SparlectraApp.handle_powerflow_case_settings_save(scf_run_id, Dict{String,String}(); output_root = root, operation_log = root)
      @test scf_save.status == 200
      @test occursin("Saved case configuration:", String(scf_save.body))
      @test read(scf_case, String) == scf_before        # the case file is untouched
      scf_cc = Sparlectra.load_case_config(scf_case)
      @test scf_cc["power_flow.max_iter"] == 42         # case scope arrives
      @test !haskey(scf_cc, "power_flow.autodamp")      # inherited values stay out
      @test !haskey(scf_cc, "matpower_import.ratio")
      @test !haskey(scf_cc, "output.logfile_results")   # machine scope stays out
      @test !haskey(scf_cc, "benchmark.enabled")
      @test all(Sparlectra.scf_is_case_config_key(k) for k in keys(scf_cc))
      # the file is still a readable revision-1.1 case after the save
      @test length(importSCF(scf_case).nodeVec) == 9
      # reset deletes the case configuration file
      reset_resp = SparlectraApp.handle_powerflow_case_settings_reset(Dict{String,Any}("casefile" => "warmup_casePST.scf.json"); output_root = root, case_directory = root, operation_log = root)
      @test reset_resp.status in (302, 303)
      @test occursin("deleted", String(Dict(reset_resp.headers)["Location"]))
      @test !isfile(Sparlectra.case_config_path(scf_case))

      # the Settings page shows the saved case settings by default;
      # ?case_settings=1 selects that view explicitly and keeps working
      loaded_form = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow/settings?casefile=$(SparlectraApp._webui_urlencode(joinpath(root, "case145.m")))&case_settings=1"; output_root = root).body)
      @test occursin("Case-specific settings loaded from", loaded_form)
      @test occursin("case145.m.config.yaml", loaded_form)
      _webui_assert_value(loaded_form, "power_flow_tol", "1.0e-7")
      _webui_assert_checked(loaded_form, "power_flow_autodamp", true)
      _webui_assert_selected(loaded_form, "power_flow_qlimits_enforcement_mode", "active_set")
      # the CSV format is machine scope: the configuration file's value
      # shows, whatever the run form said before the save
      _webui_assert_selected(loaded_form, "detailed_result_csv_format", "auto")

      notice_off_root = mktempdir()
      notice_off_config = joinpath(notice_off_root, "configuration.yaml")
      write(notice_off_config, "webui:\n  show_case_settings_notice: false\n")
      @test Sparlectra.SparlectraConfig(Sparlectra.load_yaml_dict(notice_off_config)).webui.show_case_settings_notice === false
      @test Sparlectra.SparlectraConfig(Dict()).webui.show_case_settings_notice === true
      @test SparlectraApp._powerflow_show_case_settings_notice(notice_off_config) === false
      @test SparlectraApp._powerflow_show_case_settings_notice("") === true
      @test SparlectraApp._powerflow_show_case_settings_notice(joinpath(notice_off_root, "does_not_exist.yaml")) === true
      notice_off_case = joinpath(notice_off_root, "case145.m")
      write(notice_off_case, "% case fixture\n")
      write(SparlectraApp._webui_case_settings_path(notice_off_root, notice_off_case), """
config_version: 1
scope: case
case: case145.m
power_flow:
  tol: 1.0e-7
""")
      notice_off_form = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow/settings?casefile=$(SparlectraApp._webui_urlencode(notice_off_case))&config_file=$(SparlectraApp._webui_urlencode(notice_off_config))"; output_root = notice_off_root).body)
      @test !occursin("Case-specific settings loaded from", notice_off_form)
      @test occursin("power_flow_tol", notice_off_form) # form itself still prefilled from the profile

      dismiss_root = mktempdir()
      dismiss_config = joinpath(dismiss_root, "configuration.yaml")
      write(dismiss_config, "power_flow:\n  tol: 1.0e-6\n")
      dismiss_case = joinpath(dismiss_root, "case14.m")
      write(dismiss_case, "% case fixture\n")
      write(SparlectraApp._webui_case_settings_path(dismiss_root, dismiss_case), """
config_version: 1
scope: case
case: case14.m
power_flow:
  tol: 1.0e-7
""")
      before_dismiss_form = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow/settings?casefile=$(SparlectraApp._webui_urlencode(dismiss_case))&config_file=$(SparlectraApp._webui_urlencode(dismiss_config))&case_settings=1"; output_root = dismiss_root).body)
      @test occursin("Case-specific settings loaded from", before_dismiss_form)
      @test occursin("action=\"/powerflow/config/dismiss-case-settings-notice\"", before_dismiss_form)
      @test occursin("name=\"config_file\" value=\"$(SparlectraApp._webui_escape(dismiss_config))\"", before_dismiss_form)
      @test occursin("class=\"link-button\"", before_dismiss_form)
      dismiss_response = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/config/dismiss-case-settings-notice", Dict{String,String}("config_file" => dismiss_config, "casefile" => dismiss_case); output_root = dismiss_root)
      @test dismiss_response.status == 303
      @test Dict(dismiss_response.headers)["Location"] == "/powerflow/settings?casefile=$(SparlectraApp._webui_urlencode(dismiss_case))"
      dismiss_config_text = read(dismiss_config, String)
      @test occursin("show_case_settings_notice: false", dismiss_config_text)
      @test occursin("tol: 1.0e-6", dismiss_config_text) # unrelated existing settings preserved
      after_dismiss_form = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow?casefile=$(SparlectraApp._webui_urlencode(dismiss_case))&config_file=$(SparlectraApp._webui_urlencode(dismiss_config))"; output_root = dismiss_root).body)
      @test !occursin("Case-specific settings loaded from", after_dismiss_form)
      missing_config_dismiss = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/config/dismiss-case-settings-notice", Dict{String,String}("config_file" => joinpath(dismiss_root, "does_not_exist.yaml"), "casefile" => dismiss_case); output_root = dismiss_root)
      @test missing_config_dismiss.status == 400
      empty_config_dismiss = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/config/dismiss-case-settings-notice", Dict{String,String}("casefile" => dismiss_case); output_root = dismiss_root)
      @test empty_config_dismiss.status == 400

      case118 = joinpath(root, "case118.m")
      write(case118, "% case fixture\n")
      # the MATPOWER import conventions render on the Case page; the saved
      # case configuration must prefill THAT form
      write(Sparlectra.case_config_path(case118), """
config_version: 1
scope: case
case: case118.m
model:
  auto_profile: apply
  bus_shunt_model: admittance
matpower_import:
  ratio: reciprocal
  shift_sign: -1.0
  shift_unit: rad
  pv_voltage_source: gen_vg
  compare_voltage_reference: imported_setpoint
power_flow:
  max_iter: 81
""")
      case118_case_page = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow/case?casefile=$(SparlectraApp._webui_urlencode(case118))"; output_root = root).body)
      _webui_assert_value(case118_case_page, "matpower_import_shift_sign", "-1.0")
      _webui_assert_selected(case118_case_page, "matpower_import_auto_profile", "apply")
      _webui_assert_selected(case118_case_page, "matpower_import_ratio", "reciprocal")
      _webui_assert_selected(case118_case_page, "matpower_import_shift_unit", "rad")
      _webui_assert_selected(case118_case_page, "matpower_import_bus_shunt_model", "admittance")
      _webui_assert_selected(case118_case_page, "matpower_import_pv_voltage_source", "gen_vg")
      _webui_assert_selected(case118_case_page, "matpower_import_compare_voltage_reference", "imported_setpoint")

      # the chooser marks a case with saved settings; the option value stays
      # the plain file name
      dropdown_case = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow/case"; output_root = root, runtime = (; case_directory = root, config_file = "configuration.yaml", operation_log = SparlectraApp.webui_operation_log_path(root), startup_config_error = nothing, runner = SparlectraApp.start_powerflow_run)).body)
      @test occursin("case118.m ★", dropdown_case)
      @test occursin("data-case-option=\"case118.m\"", dropdown_case)
      @test !occursin("data-case-option=\"case118.m ★\"", dropdown_case)
      dropdown_loaded_form = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow/settings?casefile=$(SparlectraApp._webui_urlencode(case118))&config_file=$(SparlectraApp._webui_urlencode("configuration.yaml"))"; output_root = root).body)
      _webui_assert_value(dropdown_loaded_form, "power_flow_max_iter", "81")
      _webui_assert_selected(dropdown_loaded_form, "detailed_result_csv_format", "auto")

      case14 = joinpath(root, "case14.m")
      write(case14, "% case fixture\n")
      case14_form = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow/settings?casefile=$(SparlectraApp._webui_urlencode(case14))"; output_root = root).body)
      @test !occursin("Case-specific settings loaded from", case14_form)
      _webui_assert_checked(case14_form, "power_flow_autodamp", true)
      _webui_assert_checked(case14_form, "power_flow_qlimits_enabled", true)
      _webui_assert_value(case14_form, "power_flow_tol", "1.0e-5")
      _webui_assert_value(case14_form, "power_flow_max_iter", "80")
      _webui_assert_selected(case14_form, "power_flow_start_angle_mode", "dc")
      _webui_assert_selected(case14_form, "power_flow_start_voltage_mode", "profile_blend")
      start_voltage_select = _webui_select_block(case14_form, "power_flow_start_voltage_mode")
      @test occursin("value=\"profile_blend\"", start_voltage_select)
      @test !occursin("bus_vm_va_blend", start_voltage_select)
      # the auto-profile selector renders on the Case page, in the
      # experimental convention block, at the template default off (0.30.0)
      case14_case_page = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow/case?casefile=$(SparlectraApp._webui_urlencode(case14))"; output_root = root).body)
      _webui_assert_selected(case14_case_page, "matpower_import_auto_profile", "off")
      _webui_assert_selected(case14_form, "output_logfile_results", "full")

      request_form = _webui_test_form("case145.m", "configuration.yaml", root)
      request_form["power_flow_tol"] = "2e-6"
      request = SparlectraApp.powerflow_webui_request(request_form; default_output_root = root)
      @test request["config_overrides"]["power_flow.tol"] == 2.0e-6

      write_solution_form = _webui_test_form("case145.m", "configuration.yaml", root)
      write_solution_form["matpower_export_write_solution"] = "false"
      write_solution_request = SparlectraApp.powerflow_webui_request(write_solution_form; default_output_root = root)
      @test write_solution_request["config_overrides"]["matpower_export.write_solution"] === false

      yaml_first_form = copy(request_form)
      yaml_first_form["ignore_webui_settings"] = "true"
      yaml_first_request = SparlectraApp.powerflow_webui_request(yaml_first_form; default_output_root = root)
      @test isempty(yaml_first_request["config_overrides"])
      @test yaml_first_request["config_override_source"] == "user_yaml"
      yaml_result_html = SparlectraApp.render_powerflow_result(Dict(
        "run_id" => "yaml-first",
        "status" => "not_converged",
        "success" => false,
        "metadata" => Dict("config_override_source" => "user_yaml"),
      ))
      @test occursin("Web UI settings ignored.", yaml_result_html)
      @test occursin("This run used YAML/default configuration values", yaml_result_html)

      invalid_case = joinpath(root, "invalid.m")
      write(invalid_case, "% case fixture\n")
      write(SparlectraApp._webui_case_settings_path(root, invalid_case), "not: [valid\n")
      invalid_form = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow/settings?casefile=$(SparlectraApp._webui_urlencode(invalid_case))"; output_root = root).body)
      @test occursin("Settings", invalid_form)
      @test !occursin("Case-specific settings loaded from", invalid_form)
      @test occursin("case_settings_load_failed", read(SparlectraApp.webui_operation_log_path(root), String))

      unsupported_case = joinpath(root, "unsupported_field.m")
      write(unsupported_case, "% case fixture\n")
      write(SparlectraApp._webui_case_settings_path(root, unsupported_case), """
config_version: 1
scope: case
case: unsupported_field.m
power_flow:
  autodamp: false
  tol: 9.0e-7
form:
  power_flow_start_angle_mode: [bad]
  detailed_result_csv_format: impossible
""")
      unsupported_form = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow/settings?casefile=$(SparlectraApp._webui_urlencode(unsupported_case))&case_settings=1"; output_root = root).body)
      @test occursin("Settings", unsupported_form)
      _webui_assert_checked(unsupported_form, "power_flow_autodamp", false)
      _webui_assert_value(unsupported_form, "power_flow_tol", "9.0e-7")
      _webui_assert_selected(unsupported_form, "power_flow_start_angle_mode", "dc")
      _webui_assert_selected(unsupported_form, "detailed_result_csv_format", "auto")
      @test occursin("case_settings_field_ignored", read(SparlectraApp.webui_operation_log_path(root), String))

      fresh_root = mktempdir()
      mkpath(joinpath(fresh_root, "resolved"))
      write(joinpath(fresh_root, "resolved", "case145.m"), "% case fixture\n")
      fresh_form = _webui_test_form("case145.m", "configuration.yaml", fresh_root)
      fresh_form["power_flow_tol"] = "2.5e-7"
      fresh_form["power_flow_max_iter"] = "37"
      fresh_form["power_flow_qlimits_enforcement_mode"] = "active_set"
      fresh_form["detailed_result_csv"] = "on"
      fresh_form["detailed_result_csv_format"] = "excel_us"
      fresh_request = SparlectraApp.powerflow_webui_request(fresh_form; default_output_root = fresh_root)
      fresh_runner = function(worker_request; case_directory = nothing)
        fresh_run_id = worker_request["run_id"]
        output_dir = joinpath(worker_request["output_root"], fresh_run_id)
        mkpath(output_dir)
        result = SparlectraApp._api_result(
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
        SparlectraApp._POWERFLOW_SERVICE_RUNS[fresh_run_id] = result
        return SparlectraApp.to_dict(result)
      end
      fresh_active = SparlectraApp.start_webui_powerflow_run(fresh_request; runner = fresh_runner)
      fresh_run_id = fresh_active["run_id"]
      for _ in 1:100
        get(SparlectraApp.get_webui_powerflow_job(fresh_run_id), "status", "") in SparlectraApp._POWERFLOW_WEBUI_ACTIVE_STATES || break
        sleep(0.05)
      end
      fresh_completed = SparlectraApp.get_webui_powerflow_job(fresh_run_id)
      @test get(fresh_completed, "status", "") == "success"
      @test get(get(fresh_completed, "metadata", Dict{String,Any}()), "webui_request_settings", nothing) isa AbstractDict
      fresh_response = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/result/$(fresh_run_id)/case-settings/save", Dict{String,String}(); output_root = fresh_root)
      @test fresh_response.status == 200
      fresh_profile_path = SparlectraApp._webui_case_settings_path(fresh_root, joinpath(fresh_root, "resolved", "case145.m"))
      @test isfile(fresh_profile_path)
      fresh_case_path = joinpath(fresh_root, "resolved", "case145.m")
      fresh_cc = Sparlectra.load_case_config(fresh_case_path)
      @test all(Sparlectra.scf_is_case_config_key(k) for k in keys(fresh_cc))
      @test !haskey(fresh_cc, "power_flow.qlimits.enforcement_mode")   # active_set is inherited
      @test fresh_cc["power_flow.tol"] == 2.5e-7
      @test fresh_cc["power_flow.max_iter"] == 37
      fresh_form = SparlectraApp._webui_case_form_defaults(fresh_case_path, joinpath(fresh_root, "resolved"))
      @test !haskey(fresh_form, "detailed_result_csv")   # true is the spec default
      @test !haskey(fresh_form, "detailed_result_csv_format")
      @test !haskey(fresh_cc, "output.csv_format")
      @test !haskey(fresh_form, "effective_config")
      fresh_reloaded_form = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow/settings?casefile=$(SparlectraApp._webui_urlencode(joinpath(fresh_root, "resolved", "case145.m")))&case_settings=1"; output_root = fresh_root).body)
      _webui_assert_value(fresh_reloaded_form, "power_flow_tol", "2.5e-7")
      _webui_assert_value(fresh_reloaded_form, "power_flow_max_iter", "37")
      _webui_assert_selected(fresh_reloaded_form, "power_flow_qlimits_enforcement_mode", "active_set")
      _webui_assert_checked(fresh_reloaded_form, "detailed_result_csv", true)
      manual_override_form = _webui_test_form(joinpath(fresh_root, "resolved", "case145.m"), "configuration.yaml", fresh_root)
      manual_override_form["power_flow_tol"] = "3e-6"
      manual_override_request = SparlectraApp.powerflow_webui_request(manual_override_form; default_output_root = fresh_root)
      @test manual_override_request["config_overrides"]["power_flow.tol"] == 3.0e-6
      manual_override_nested = Sparlectra.validate_gui_config_overrides(manual_override_request["config_overrides"])
      manual_override_cfg, _ = Sparlectra._load_api_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, manual_override_nested)
      @test manual_override_cfg.powerflow.tol == 3.0e-6
      delete!(SparlectraApp._POWERFLOW_SERVICE_RUNS, fresh_run_id)
      delete!(SparlectraApp._POWERFLOW_WEBUI_JOBS, fresh_run_id)

      failed_run_id = "case-settings-failed"
      failed = SparlectraApp._api_result(
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
      SparlectraApp._POWERFLOW_SERVICE_RUNS[failed_run_id] = failed
      failed_html = SparlectraApp.render_powerflow_result(SparlectraApp.to_dict(failed))
      @test !occursin("Save settings for this case", failed_html)
      @test occursin("Save these settings anyway", failed_html)
      rejected = SparlectraApp.handle_powerflow_case_settings_save(failed_run_id, Dict{String,String}(); output_root = root, operation_log = root)
      @test rejected.status == 400
      accepted = SparlectraApp.handle_powerflow_case_settings_save(failed_run_id, Dict("override_non_success" => "true"); output_root = root, operation_log = root)
      @test accepted.status == 200
      log_text = read(SparlectraApp.webui_operation_log_path(root), String)
      @test occursin("case_settings_saved", log_text)
      @test occursin("case_settings_save_failed", log_text)
      @test SparlectraApp._webui_normalized_case_key("../bad/../../case145.m") == "case145"

      unsupported_run_id = "case-settings-unsupported"
      unsupported_metadata = deepcopy(metadata)
      unsupported_metadata["webui_request_settings"]["power_flow.autodamp"] = Dates.now()
      SparlectraApp._POWERFLOW_WEBUI_JOBS[unsupported_run_id] = Dict{String,Any}(
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
      unsupported_response = SparlectraApp.handle_powerflow_case_settings_save(unsupported_run_id, Dict{String,String}(); output_root = root, operation_log = root)
      @test unsupported_response.status == 400
      @test occursin("unsupported value type", String(unsupported_response.body))
      traversal_run_id = "case-settings-traversal"
      traversal_metadata = deepcopy(metadata)
      traversal_metadata["runtime_casefile_path"] = "../outside/case145.m"
      SparlectraApp._POWERFLOW_WEBUI_JOBS[traversal_run_id] = Dict{String,Any}(
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
      traversal_response = SparlectraApp.handle_powerflow_case_settings_save(traversal_run_id, Dict{String,String}(); output_root = root, operation_log = root)
      @test traversal_response.status == 400
      delete!(SparlectraApp._POWERFLOW_SERVICE_RUNS, run_id)
      delete!(SparlectraApp._POWERFLOW_SERVICE_RUNS, failed_run_id)
      delete!(SparlectraApp._POWERFLOW_WEBUI_JOBS, unsupported_run_id)
      delete!(SparlectraApp._POWERFLOW_WEBUI_JOBS, traversal_run_id)
    end)() end

    # the bundled .jl precompile workloads stay out of the case selector
    @test !SparlectraApp._webui_is_user_selectable_case("warmup_case118.jl")

    @testset "Case and configuration selection" begin (function ()
      mktempdir() do start_dir
        application_root = joinpath(start_dir, "Sparlectra")
        case_directory = joinpath(application_root, "data", "mpower")
        config_directory = joinpath(application_root, "examples")
        mkpath(case_directory)
        mkpath(config_directory)
        case14 = joinpath(case_directory, "case14.m")
        case118 = joinpath(case_directory, "case118.m")
        casejl = joinpath(case_directory, "case14.jl")
        warmup_m = joinpath(case_directory, "warmup_internal.m")
        warmup_jl = joinpath(case_directory, "warmup_case3.jl")
        generated_cache = joinpath(case_directory, "generated_cache.jl")
        for001 = joinpath(case_directory, "FOR001.DAT")
        for002 = joinpath(case_directory, "FOR002.DAT")
        for002_variant = joinpath(case_directory, "FOR002_reference.DAT")
        write(case14, "function mpc = case14\nend\n")
        write(case118, "function mpc = case118\nend\n")
        write(casejl, "nothing\n")
        write(warmup_m, "function mpc = warmup_internal\nend\n")
        write(warmup_jl, "nothing\n")
        write(generated_cache, "nothing\n")
        write(for001, _dtf_network_fixture())
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

        @test SparlectraApp._webui_application_root(start_dir) == application_root
        @test SparlectraApp._webui_is_user_selectable_case("case3.m")
        @test SparlectraApp._webui_supported_upload_case_extension("case3.m")
        @test SparlectraApp._webui_supported_upload_case_extension("case3.M")
        @test SparlectraApp._webui_supported_upload_case_extension("FOR001.dat")
        @test SparlectraApp._webui_supported_upload_case_extension("FOR001.DAT")
        @test !SparlectraApp._webui_supported_upload_case_extension("case3.jl")
        @test SparlectraApp._webui_supported_upload_case_extension("case3.zip")   # CGMES deliveries
        @test !SparlectraApp._webui_supported_upload_case_extension("case3.txt")
        @test !SparlectraApp._webui_supported_upload_case_extension("case3")
        @test SparlectraApp._webui_sanitize_upload_filename("case3.m") == ("case3.m", "")
        @test last(SparlectraApp._webui_sanitize_upload_filename("../case3.m")) == "invalid filename"
        @test last(SparlectraApp._webui_sanitize_upload_filename("")) == "empty filename"
        @test SparlectraApp._webui_is_user_selectable_case("case14.m")
        # warmup_ hides only the .jl workloads; a MATPOWER case carrying the
        # prefix stays selectable (warmup_casePST.m regression)
        @test SparlectraApp._webui_is_user_selectable_case("warmup_case3.m")
        @test !SparlectraApp._webui_is_user_selectable_case("warmup_case3.jl")
        @test !SparlectraApp._webui_is_user_selectable_case("warmup_casePST.sparlectra-webui.yaml")
        @test SparlectraApp._webui_is_user_selectable_case("warmup_internal.m")
        @test !SparlectraApp._webui_is_user_selectable_case("generated_cache.jl")
        @test SparlectraApp._webui_casefile_options(application_root) == ["case118.m", "case14.m", "FOR001.DAT", "warmup_internal.m"]
        @test SparlectraApp._webui_for002_reference_options_in_directory(case_directory) == ["FOR002.DAT", "FOR002_reference.DAT"]
        @test SparlectraApp._webui_config_file_options(application_root) == [primary_config, secondary_config]

        # the chooser, upload form, and import-format options live
        # on the Case page; the run page carries the selected case as a
        # hidden field plus the DTF outage run parameters
        selection_html = SparlectraApp.render_case_page(
          application_root = application_root,
          selected_casefile = "case14.m",
          selected_config_file = secondary_config,
        )
        run_selection_html = SparlectraApp.render_powerflow_form(
          application_root = application_root,
          selected_casefile = "case14.m",
          selected_config_file = secondary_config,
        )
        # the value the case chooser shows, and the format the Case page selects
        case_value(html) = (m = match(r"<input id=\"casefile\" name=\"casefile\"[^>]* value=\"([^\"]*)\"", html); m === nothing ? nothing : String(m.captures[1]))
        dtf_open(html) = occursin("dtf-internal-section is-dat-selected\" open>", html)
        @test case_value(selection_html) == "case14.m"
        for name in ("case14.m", "case118.m", "FOR001.DAT", "warmup_internal.m")
          @test occursin("data-case-option=\"$(name)\"", selection_html)
        end
        @test !occursin("data-case-option=\"case14.jl\"", selection_html)
        @test !occursin("warmup_case3.jl", selection_html)
        @test !occursin("generated_cache", selection_html)
        # a FOR002 file is a reference result, never a primary case
        @test !occursin("data-case-option=\"FOR002", selection_html)
        @test !occursin("README.md", selection_html)
        for ignored in ("artifact.csv", "artifact.json", "artifact.yaml", "artifact.log", "artifact.md", "artifact.sparlectra-webui.yaml")
          @test !occursin(ignored, selection_html)
        end
        # .json is a Sparlectra Case Format case (#342); the import validates
        # the content before storing, so a foreign .json is refused by name
        @test occursin("type=\"file\" name=\"casefiles\" accept=\".m,.M,.dat,.DAT,.zip,.ZIP,.json,.yaml,.xml,.XML,.xiidm,.XIIDM\" multiple", selection_html)
        @test occursin("<input type=\"hidden\" name=\"config_file\" value=\"$(secondary_config)\">", selection_html)
        @test !occursin("name=\"for002_reference_file\"", selection_html)
        # the run page submits the SAME selected case and configuration as
        # hidden fields; the case choice and its format live on the Case page
        @test occursin("type=\"hidden\" name=\"casefile\" value=\"case14.m\"", run_selection_html)
        @test occursin("<input type=\"hidden\" name=\"config_file\" value=\"$(secondary_config)\">", run_selection_html)
        @test !occursin("name=\"case_format\"", run_selection_html)
        # the FOR002 reference and DTF outage selection are RUN parameters
        # on the run page, rendered for a DTF case only: a MATPOWER case
        # leaves the five fields out and counts them, the MATPOWER export
        # artifact (every format) stays
        @test !occursin("name=\"for002_reference_file\"", run_selection_html)
        @test !occursin("name=\"dtf_outage_selection_mode\"", run_selection_html)
        @test occursin("data-hidden-options=\"5\"", run_selection_html)
        @test occursin("name=\"matpower_export_requested\"", run_selection_html)

        # a .DAT case selects the DTF format and opens its section
        dat_html = SparlectraApp.render_case_page(
          application_root = application_root,
          selected_casefile = "FOR001.DAT",
        )
        @test case_value(dat_html) == "FOR001.DAT"
        _webui_assert_selected(dat_html, "case_format", "dtf_for001")
        @test dtf_open(dat_html)
        # the run page opens its DTF outage details for a .DAT case, carries
        # the (empty) FOR002 reference field and offers both FOR002 files
        dat_run_html = SparlectraApp.render_powerflow_form(
          application_root = application_root,
          selected_casefile = "FOR001.DAT",
        )
        @test dtf_open(dat_run_html)
        @test occursin("name=\"for002_reference_file\" value=\"\"", dat_run_html)
        @test occursin("<option value=\"FOR002.DAT\">", dat_run_html) && occursin("<option value=\"FOR002_reference.DAT\">", dat_run_html)
        @test occursin("name=\"dtf_outage_selection_mode\"", dat_run_html)
        @test occursin("name=\"matpower_export_requested\"", dat_run_html)
        @test !occursin("data-hidden-options", dat_run_html)

        for002_reference_html = SparlectraApp.render_powerflow_form(
          application_root = application_root,
          selected_casefile = "FOR001.DAT",
          submitted_form = Dict("for002_reference_file" => "/manual/FOR002.DAT"),
        )
        @test occursin("name=\"for002_reference_file\" value=\"/manual/FOR002.DAT\"", for002_reference_html)

        # an explicit format choice wins over the .DAT detection
        submitted_auto_dat_html = SparlectraApp.render_case_page(
          application_root = application_root,
          selected_casefile = "FOR001.DAT",
          submitted_form = Dict("case_format" => "auto"),
        )
        _webui_assert_selected(submitted_auto_dat_html, "case_format", "auto")

        case14_html = SparlectraApp.render_case_page(
          application_root = application_root,
          selected_casefile = "case14.m",
        )
        @test case_value(case14_html) == "case14.m"
        _webui_assert_selected(case14_html, "case_format", "auto")

        casejl_html = SparlectraApp.render_case_page(
          application_root = application_root,
          selected_casefile = "case14.jl",
        )
        @test !occursin("data-case-option=\"case14.jl\"", casejl_html)
        @test case_value(casejl_html) == "case14.jl"
        _webui_assert_selected(casejl_html, "case_format", "auto")

        manual_dat = joinpath(case_directory, "manual", "FOR001.DAT")
        manual_dat_html = SparlectraApp.render_case_page(
          application_root = application_root,
          selected_casefile = manual_dat,
        )
        @test case_value(manual_dat_html) == manual_dat
        _webui_assert_selected(manual_dat_html, "case_format", "dtf_for001")
        @test dtf_open(manual_dat_html)

        for manual_matpower in (joinpath(case_directory, "manual", "case14.m"), joinpath(case_directory, "manual", "case14.jl"))
          manual_matpower_html = SparlectraApp.render_case_page(
            application_root = application_root,
            selected_casefile = manual_matpower,
          )
          @test case_value(manual_matpower_html) == manual_matpower
          _webui_assert_selected(manual_matpower_html, "case_format", "auto")
        end

        # a typed path wins over the chooser value
        manual_overrides_dropdown_html = SparlectraApp.render_case_page(
          application_root = application_root,
          selected_casefile = "case14.m",
          submitted_form = Dict("casefile" => "case14.m", "casefile_manual" => manual_dat),
        )
        @test case_value(manual_overrides_dropdown_html) == manual_dat
        _webui_assert_selected(manual_overrides_dropdown_html, "case_format", "dtf_for001")

        rm(case14)
        rm(case118)
        manual_case = joinpath(case_directory, "manual_case.m")
        fallback_html = SparlectraApp.render_case_page(
          application_root = application_root,
          selected_casefile = manual_case,
        )
        @test case_value(fallback_html) == manual_case
      end
    end)() end

    @testset "Case file import" begin (function ()
      mktempdir() do tmpdir
        case_directory = joinpath(tmpdir, "cases")
        output_root = joinpath(tmpdir, "runs")
        runtime = SparlectraApp._SparlectraWebUIRuntime(nothing, case_directory, "configuration.yaml", SparlectraApp.webui_operation_log_path(output_root), nothing, SparlectraApp.start_powerflow_run, false, false, time(), 0, nothing, IOBuffer(), ReentrantLock())
        upload(name, text) = SparlectraApp.WebUICaseUpload(name, Vector{UInt8}(codeunits(text)))

        response = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/import-cases", Dict("casefiles" => [upload("case_upload.m", "function mpc = case_upload\nend\n")]); output_root, runtime)
        @test response.status == 303
        @test isfile(joinpath(case_directory, "case_upload.m"))
        @test read(joinpath(case_directory, "case_upload.m"), String) == "function mpc = case_upload\nend\n"
        refreshed = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow/case"; output_root, runtime).body)
        @test occursin("data-case-option=\"case_upload.m\"", refreshed)

        dat_response = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/import-cases", Dict("casefiles" => [upload("FOR001.DAT", _dtf_network_fixture())]); output_root, runtime)
        @test dat_response.status == 303
        @test isfile(joinpath(case_directory, "FOR001.DAT"))
        @test SparlectraApp._webui_casefile_options_in_directory(case_directory) == ["case_upload.m", "FOR001.DAT"]

        multi = Dict("casefiles" => [upload("case_a.M", "a"), upload("case_b.dat", "b"), upload("bad.txt", "z")])
        multi_response = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/import-cases", multi; output_root, runtime)
        @test multi_response.status == 303
        @test isfile(joinpath(case_directory, "case_a.M"))
        @test isfile(joinpath(case_directory, "case_b.dat"))
        @test !isfile(joinpath(case_directory, "bad.txt"))
        multi_page = String(SparlectraApp.route_sparlectra_webui("GET", String(Dict(multi_response.headers)["Location"]); output_root, runtime).body)
        @test occursin("rejected 1", multi_page)
        @test occursin("bad.txt", multi_page) && occursin("unsupported extension", multi_page)

        # a stray manifest.json (the index of a table bundle) is no case
        # file and is refused like any other JSON that is none
        manifest = Dict("casefiles" => [upload("manifest.json", "{\"format\": \"powsybl_tables\", \"format_version\": 1, \"case\": \"x\", \"source\": \"test\"}")])
        manifest_response = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/import-cases", manifest; output_root, runtime)
        @test manifest_response.status == 303
        @test !isfile(joinpath(case_directory, "manifest.json"))
        # the .xiidm file itself is a case file: the upload stores it, the case runs by its bare name (it used to be
        # looked up as a MATPOWER case to download, and the run service had
        # no branch for the format), and a measurement set generated for it
        # estimates; a name that is not there is a plain not-found
        iidm_source = joinpath(Sparlectra.SPARLECTRA_ROOT, "data", "powsybl", "ieee14.xiidm")
        iidm_response = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/import-cases", Dict("casefiles" => [SparlectraApp.WebUICaseUpload("ieee14.xiidm", read(iidm_source))]); output_root, runtime)
        @test iidm_response.status == 303
        @test isfile(joinpath(case_directory, "ieee14.xiidm"))
        @test occursin("data-case-option=\"ieee14.xiidm\"", String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow/case"; output_root, runtime).body))
        iidm_run = redirect_stdout(devnull) do
          SparlectraApp.start_powerflow_run(Dict{String,Any}("casefile" => "ieee14.xiidm", "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => joinpath(tmpdir, "runs_iidm")); case_directory = case_directory)
        end
        @test iidm_run["status"] == "succeeded"
        @test iidm_run["converged"] === true
        # the import report is an artifact of the run (the importer prints
        # nothing, so a run from the Web UI never showed the slack choice)
        iidm_report = joinpath(iidm_run["output_dir"], "powsybl_import.log")
        @test isfile(iidm_report)
        @test occursin(r"slack component \d+: B1-G", read(iidm_report, String))
        @test occursin("PowSyBl import report: powsybl_import.log", read(joinpath(iidm_run["output_dir"], "run.log"), String))
        bundle_run = redirect_stdout(devnull) do
          SparlectraApp.start_powerflow_run(Dict{String,Any}("casefile" => joinpath(Sparlectra.SPARLECTRA_ROOT, "test", "fixtures", "powsybl", "ieee14.powsybl"), "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => joinpath(tmpdir, "runs_bundle")); case_directory = case_directory)
        end
        @test bundle_run["status"] == "succeeded"
        missing_run = SparlectraApp.start_powerflow_run(Dict{String,Any}("casefile" => "not_there.xiidm", "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => joinpath(tmpdir, "runs_missing")); case_directory = case_directory)
        @test missing_run["status"] == "failed"
        @test occursin("Case file not found: not_there.xiidm", missing_run["message"])
        # (the runtime of this test set names a configuration file that does
        # not exist; the generator needs a real one)
        iidm_runtime = SparlectraApp._SparlectraWebUIRuntime(nothing, case_directory, Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, SparlectraApp.webui_operation_log_path(output_root), nothing, SparlectraApp.start_powerflow_run, false, false, time(), 0, nothing, IOBuffer(), ReentrantLock())
        # every other run kind goes through on the same uploaded file: the
        # N-1 sweep, the short circuit (the generators are the sources; this
        # file carries no short-circuit data, so every row is flagged and
        # the message says why) and Diagnose (its verdict on the file's
        # start state, as for any other format)
        # (a runtime with a configuration file that exists, for the page and the routes)
        export_runtime = SparlectraApp._SparlectraWebUIRuntime(nothing, case_directory, Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, SparlectraApp.webui_operation_log_path(output_root), nothing, SparlectraApp.start_powerflow_run, false, false, time(), 0, nothing, IOBuffer(), ReentrantLock())
        iidm_kind(tag, extra) = redirect_stdout(devnull) do
          SparlectraApp.start_powerflow_run(merge(Dict{String,Any}("casefile" => "ieee14.xiidm", "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => joinpath(tmpdir, "runs_iidm_" * tag)), extra); case_directory = case_directory)
        end
        iidm_n1 = iidm_kind("n1", Dict{String,Any}("contingency_mode" => true))
        @test iidm_n1["status"] == "succeeded"
        # the outage of line 1-2 takes away the path of the reference unit's
        # output: on a single slack it has no solution with the reactive
        # limits (OpenLoadFlow ends the same way in that setting), with the
        # distributed slack of the run configuration it has one (lowest
        # voltage 0.9256 pu; OpenLoadFlow 0.9243), and the service passes
        # that configuration to every post-outage solve
        @test iidm_n1["metadata"]["contingency_converged"] == 19
        iidm_n1_shared = iidm_kind("n1_shared", Dict{String,Any}("contingency_mode" => true, "config_overrides" => Dict{String,Any}("power_flow.distributed_slack.enabled" => true, "power_flow.distributed_slack.p_mode" => "imported")))
        @test iidm_n1_shared["status"] == "succeeded"
        @test iidm_n1_shared["metadata"]["contingency_converged"] == 20
        # generator outages: the loss of the reference unit hands the
        # reference to the strongest unit that is left, named in the row
        iidm_n1_units = iidm_kind("n1_units", Dict{String,Any}("contingency_mode" => true, "contingency_kind" => "gen", "config_overrides" => Dict{String,Any}("power_flow.distributed_slack.enabled" => true, "power_flow.distributed_slack.p_mode" => "imported")))
        @test iidm_n1_units["status"] == "succeeded"
        @test occursin("reference taken over by bus", read(joinpath(iidm_n1_units["output_dir"], "contingency_n1.csv"), String))
        @test !occursin("no slack bus", read(joinpath(iidm_n1_units["output_dir"], "contingency_n1.csv"), String))
        # short circuit, three states of the data. The shipped file carries
        # neither short-circuit data nor a rated power on any generator:
        # nothing of a source impedance would come from the file, the run
        # refuses. With rated powers the table is complete on the default
        # reactance and the run ends with a WARNING that leads the message.
        # With the extension on every generator it is a plain success.
        iidm_sc = iidm_kind("sc", Dict{String,Any}("short_circuit_mode" => true))
        @test iidm_sc["status"] == "failed"
        @test iidm_sc["reason"] == "short_circuit_data_missing"
        @test occursin("generatorShortCircuit", iidm_sc["message"]) && occursin("ratedS", iidm_sc["message"])
        iidm_text = read(iidm_source, String)
        iidm_ids = [String(m.captures[1]) for m in eachmatch(r"<iidm:generator id=\"([^\"]+)\"", iidm_text)]
        rated_text = replace(iidm_text, r"<iidm:generator id=" => "<iidm:generator ratedS=\"100.0\" id=")
        write(joinpath(case_directory, "ieee14_rated.xiidm"), rated_text)
        sc_extension = join("    <iidm:extension id=\"$(id)\">\n        <gsc:generatorShortCircuit xmlns:gsc=\"http://www.powsybl.org/schema/iidm/ext/generator_short_circuit/1_0\" directSubtransX=\"20.0\"/>\n    </iidm:extension>\n" for id in iidm_ids)
        write(joinpath(case_directory, "ieee14_sc.xiidm"), replace(rated_text, "</iidm:network>" => sc_extension * "</iidm:network>"))
        sc_run(case, tag) = redirect_stdout(devnull) do
          SparlectraApp.start_powerflow_run(Dict{String,Any}("casefile" => case, "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => joinpath(tmpdir, "runs_sc_" * tag), "short_circuit_mode" => true); case_directory = case_directory)
        end
        sc_defaults = sc_run("ieee14_rated.xiidm", "defaults")
        @test sc_defaults["status"] == "warning"
        @test sc_defaults["success"] === true
        @test sc_defaults["reason"] == "short_circuit_defaults_only"
        @test startswith(sc_defaults["message"], "Short circuit completed on defaults: none of the 5 source(s) carries short-circuit data")
        @test sc_defaults["metadata"]["sc_sources"] == 5 && sc_defaults["metadata"]["sc_sources_with_data"] == 0 && sc_defaults["metadata"]["sc_defaults_only"] === true
        @test startswith(read(joinpath(sc_defaults["output_dir"], "run.log"), String), "WARNING: none of the 5 generator(s) carries short-circuit data")
        @test isfile(joinpath(sc_defaults["output_dir"], "short_circuit_max.csv")) && isfile(joinpath(sc_defaults["output_dir"], "short_circuit_min.csv"))
        @test SparlectraApp.webui_status_class(sc_defaults) == "status-warning"
        @test occursin("no short-circuit data in the case", SparlectraApp._webui_short_circuit_summary(sc_defaults))
        sc_data = sc_run("ieee14_sc.xiidm", "data")
        @test sc_data["status"] == "succeeded"
        @test sc_data["metadata"]["sc_sources_with_data"] == 5 && sc_data["metadata"]["sc_flagged_rows"] == 0
        @test sc_data["message"] == "Short circuit completed."
        # the shipped demo file with short-circuit data runs the same way
        sc_shipped = sc_run(joinpath(Sparlectra.SPARLECTRA_ROOT, "data", "powsybl", "ieee14_sc.xiidm"), "shipped")
        @test sc_shipped["status"] == "succeeded"
        @test sc_shipped["metadata"]["sc_sources_with_data"] == 5 && sc_shipped["metadata"]["sc_flagged_rows"] == 0
        # the Short circuit button knows the three states before a run: not
        # offered with its reason for the file without data and rating,
        # offered for the two others
        button_state(case) = match(r"data-sc-state=\"([a-z-]+)\"", String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow?casefile=$(case)"; output_root, runtime = export_runtime).body)).captures[1]
        @test button_state("ieee14.xiidm") == "missing-data"
        @test button_state("ieee14_rated.xiidm") == "ready"
        @test button_state("ieee14_sc.xiidm") == "ready"
        # pageshow also fires on the first load, and the run form's handler
        # sets the disabled flag of every submit button. Evaluated on the
        # buttons of the rendered page (its expression translated to Julia,
        # the suite has no browser), it keeps each button the server
        # disabled (Short circuit without data, QSTS without a series) and
        # enables the others: per button the server's flag => the rule's.
        function pageshow_disabled(case)
          html = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow?casefile=$(case)"; output_root, runtime = export_runtime).body)
          handler = match(r"window\.addEventListener\('pageshow', function \(\) \{\s*const form = document\.getElementById\('powerflow-run-form'\);(.*?)\n\}\);"s, html).captures[1]
          expr = match(r"b\.disabled = ([^;]*);", handler).captures[1]
          jexpr = replace(expr, r"b\.hasAttribute\('([^']*)'\)" => s"has(\"\1\")", r"b\.getAttribute\('([^']*)'\)" => s"attr(\"\1\")", "!==" => "!=", "===" => "==", r"'([^']*)'" => s"\"\1\"")
          rule = eval(Meta.parse("(has, attr) -> " * jexpr))
          states = Pair{Bool,Bool}[]
          for tag in eachmatch(r"<button class=\"powerflow-submit[^>]*>", html)
            attrs = Dict(String(a[1]) => something(a[2], "") for a in eachmatch(r"([a-z-]+)(?:=\"([^\"]*)\")?", tag.match))
            push!(states, haskey(attrs, "disabled") => Base.invokelatest(rule, k -> haskey(attrs, k), k -> get(attrs, k, nothing)))
          end
          return states
        end
        gated = pageshow_disabled("ieee14.xiidm")
        # neither short-circuit data nor a series: Short circuit and QSTS
        @test count(first, gated) == 2 && all(p -> last(p) == first(p), gated)
        # with short-circuit data only QSTS stays disabled
        sc_gated = pageshow_disabled("ieee14_sc.xiidm")
        @test count(first, sc_gated) == 1 && all(p -> last(p) == first(p), sc_gated)
        # the three exports: the SCF file through its route, the MATPOWER
        # and CGMES artifacts with a run; every export runs again through the
        # service and lands on the voltages of the source run
        scf_export = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/export-scf", Dict{String,Any}("casefile" => "ieee14.xiidm", "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH); output_root, runtime = export_runtime)
        @test occursin("Exported ieee14.scf.json", SparlectraApp._webui_urldecode(Dict(scf_export.headers)["Location"]))
        @test isfile(joinpath(case_directory, "ieee14.scf.json"))
        iidm_exports = iidm_kind("exports", Dict{String,Any}("matpower_export_requested" => true, "export_cgmes" => true, "detailed_result_csv" => true))
        @test iidm_exports["status"] == "succeeded"
        @test iidm_exports["metadata"]["matpower_export_status"] == "completed"
        @test iidm_exports["metadata"]["cgmes_export_status"] == "completed"
        cp(joinpath(iidm_exports["output_dir"], String(iidm_exports["metadata"]["matpower_export_file"])), joinpath(case_directory, "ieee14_export.m"))
        cp(joinpath(iidm_exports["output_dir"], String(iidm_exports["metadata"]["cgmes_export_files"])), joinpath(case_directory, "ieee14_export_cgmes.zip"))
        magnitudes(result) = begin
          csv = run_csv(joinpath(result["output_dir"], "bus_voltages_complex.csv"))
          column = findfirst(==("vm_pu"), csv.header)
          sort!([csv.number(row[column]) for row in csv.rows])
        end
        source_vm = magnitudes(iidm_exports)
        for exported in ("ieee14_export.m", "ieee14_export_cgmes.zip", "ieee14.scf.json")
          back = redirect_stdout(devnull) do
            SparlectraApp.start_powerflow_run(Dict{String,Any}("casefile" => exported, "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => joinpath(tmpdir, "runs_back_" * exported), "detailed_result_csv" => true); case_directory = case_directory)
          end
          @test back["status"] == "succeeded"
          back_vm = magnitudes(back)
          @test length(back_vm) == length(source_vm)
          @test maximum(abs.(back_vm .- source_vm)) < 1e-5
        end
        iidm_diag = iidm_kind("diag", Dict{String,Any}("diagnose_mode" => true))
        @test isfile(joinpath(iidm_diag["output_dir"], "diagnose.log"))
        @test iidm_diag["metadata"]["run_mode"] == "diagnose"
        iidm_gen = SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "ieee14.xiidm", "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH); output_root, runtime = iidm_runtime)
        iidm_set = SparlectraApp._webui_urldecode(String(match(r"measurement_file=([^&#]+)", Dict(iidm_gen.headers)["Location"]).captures[1]))
        iidm_se = redirect_stdout(devnull) do
          SparlectraApp.start_powerflow_run(Dict{String,Any}("casefile" => "ieee14.xiidm", "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => joinpath(tmpdir, "runs_iidm_se"), "se_mode" => true, "measurement_file" => iidm_set); case_directory = case_directory)
        end
        @test iidm_se["status"] == "succeeded"

        # An incomplete CGMES delivery triggers the upload-time import
        # analysis: the message names the missing declared dependency and the
        # full report lands next to the case.
        incomplete_eq = """<?xml version="1.0" encoding="UTF-8"?>
<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:cim="http://iec.ch/TC57/2013/CIM-schema-cim16#" xmlns:md="http://iec.ch/TC57/61970-552/ModelDescription/1#">
<md:FullModel rdf:about="urn:uuid:upload-eq-model">
<md:Model.profile>http://entsoe.eu/CIM/EquipmentCore/3/1</md:Model.profile>
<md:Model.DependentOn rdf:resource="urn:uuid:absent-boundary"/>
</md:FullModel>
<cim:TopologicalNode rdf:ID="_upl_tn"><cim:IdentifiedObject.name>UplTN</cim:IdentifiedObject.name><cim:TopologicalNode.BaseVoltage rdf:resource="#_bv_in_boundary"/></cim:TopologicalNode>
</rdf:RDF>
"""
        zipbuf = IOBuffer()
        ZipArchives.ZipWriter(zipbuf) do w
          ZipArchives.zip_newfile(w, "incomplete_EQ.xml")
          write(w, incomplete_eq)
        end
        cgmes_response = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/import-cases", Dict("casefiles" => [SparlectraApp.WebUICaseUpload("incomplete_delivery.zip", take!(zipbuf))]); output_root, runtime)
        @test cgmes_response.status == 303
        @test isfile(joinpath(case_directory, "incomplete_delivery.zip"))
        analysis_file = joinpath(case_directory, "incomplete_delivery.import_analysis.txt")
        @test isfile(analysis_file)
        analysis_text = read(analysis_file, String)
        @test occursin("absent-boundary", analysis_text)
        @test occursin("Verdict:", analysis_text)
        cgmes_page = String(SparlectraApp.route_sparlectra_webui("GET", String(Dict(cgmes_response.headers)["Location"]); output_root, runtime).body)
        @test occursin("missing declared dependencies", cgmes_page)
        @test occursin("absent-boundary", cgmes_page)

        # the profile files of one delivery, selected together in the file
        # picker (the shipped demo folders are exactly that), become one ZIP
        demo = joinpath(dirname(@__DIR__), "data", "cgmes_demo", "sp_case14")
        xml_uploads = [SparlectraApp.WebUICaseUpload(f, read(joinpath(demo, f))) for f in readdir(demo)]
        xml_response = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/import-cases", Dict("casefiles" => xml_uploads); output_root, runtime)
        @test xml_response.status == 303
        @test isfile(joinpath(case_directory, "sp_case14_cgmes.zip"))
        @test occursin("sp_case14_cgmes.zip", SparlectraApp._webui_urldecode(String(Dict(xml_response.headers)["Location"])))
        @test length(importCGMES(path = joinpath(case_directory, "sp_case14_cgmes.zip"), name = "sp_case14_cgmes").net.nodeVec) == 14
        # without the EQ profile the set is refused, nothing is written
        no_eq = [u for u in xml_uploads if !occursin("_EQ", u.filename)]
        no_eq_response = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/import-cases", Dict("casefiles" => no_eq); output_root, runtime)
        @test occursin("EQ profile", SparlectraApp._webui_urldecode(String(Dict(no_eq_response.headers)["Location"])))

        write(joinpath(case_directory, "duplicate.m"), "old")
        duplicate_response = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/import-cases", Dict("casefiles" => [upload("duplicate.m", "new")]); output_root, runtime)
        @test duplicate_response.status == 303
        @test read(joinpath(case_directory, "duplicate.m"), String) == "old"

        traversal_response = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/import-cases", Dict("casefiles" => [upload("../escape.m", "bad")]); output_root, runtime)
        @test traversal_response.status == 303
        @test !isfile(joinpath(tmpdir, "escape.m"))

        oversized = SparlectraApp.handle_powerflow_case_import(Dict("casefiles" => [upload("huge.m", "12345")]); output_root, case_directory, operation_log = output_root, max_file_bytes = 4)
        @test oversized.status == 303
        @test !isfile(joinpath(case_directory, "huge.m"))

        empty_response = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/import-cases", Dict("casefiles" => SparlectraApp.WebUICaseUpload[]); output_root, runtime)
        @test empty_response.status == 303

        for002_response = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/import-cases", Dict("casefiles" => [upload("FOR002.DAT", "reference")]); output_root, runtime)
        @test for002_response.status == 303
        @test isfile(joinpath(case_directory, "FOR002.DAT"))
        # the case selector (Case page) must not offer the FOR002 reference
        # file as a case; the run page's FOR002 datalist is separate
        selector_html = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow/case"; output_root, runtime).body)
        @test !occursin("data-case-option=\"FOR002.DAT\"", selector_html)
        run_selector_html = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow"; output_root, runtime).body)
        primary_selector_html = split(run_selector_html, "<datalist id=\"for002-reference-candidates\">")[1]
        @test !occursin("<option value=\"FOR002.DAT\">FOR002.DAT</option>", primary_selector_html)

        @test !isdir(joinpath(output_root, "runs"))
        @test isempty(SparlectraApp.list_powerflow_runs(output_root))
        log_text = read(SparlectraApp.webui_operation_log_path(output_root), String)
        @test occursin("case_import_completed", log_text)
        @test !occursin("powerflow_run_started", log_text)
      end
    end)() end

    # Enter in the manual "type case file path" field must resolve/copy the
    # case into the case directory without starting a PowerFlow run, so it
    # then appears in the "choose existing case" selector.
    @testset "Manual case path resolve (no PowerFlow run)" begin (function ()
      mktempdir() do tmpdir
        case_directory = joinpath(tmpdir, "cases")
        output_root = joinpath(tmpdir, "runs")
        runtime = SparlectraApp._SparlectraWebUIRuntime(nothing, case_directory, "configuration.yaml", SparlectraApp.webui_operation_log_path(output_root), nothing, SparlectraApp.start_powerflow_run, false, false, time(), 0, nothing, IOBuffer(), ReentrantLock())

        external_dir = joinpath(tmpdir, "external")
        mkpath(external_dir)
        external_dat = joinpath(external_dir, "FOR001_manual.DAT")
        write(external_dat, _dtf_network_fixture())

        response = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/resolve-case", Dict("casefile_manual" => external_dat); output_root, runtime)
        @test response.status == 303
        @test Dict(response.headers)["Location"] == "/powerflow/case?casefile=FOR001_manual.DAT&import_message=Resolved%20case%3A%20FOR001_manual.DAT"
        @test isfile(joinpath(case_directory, "FOR001_manual.DAT"))
        @test !isdir(joinpath(output_root, "runs"))
        @test isempty(SparlectraApp.list_powerflow_runs(output_root))

        refreshed = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow/case"; output_root, runtime).body)
        @test occursin("data-case-option=\"FOR001_manual.DAT\"", refreshed)

        duplicate_response = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/resolve-case", Dict("casefile_manual" => external_dat); output_root, runtime)
        @test duplicate_response.status == 303
        @test occursin("already exists", SparlectraApp._webui_urldecode(Dict(duplicate_response.headers)["Location"]))

        bad_ext = joinpath(external_dir, "notes.txt")
        write(bad_ext, "not a case")
        bad_ext_response = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/resolve-case", Dict("casefile_manual" => bad_ext); output_root, runtime)
        @test occursin("unsupported extension", SparlectraApp._webui_urldecode(Dict(bad_ext_response.headers)["Location"]))

        missing_response = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/resolve-case", Dict("casefile_manual" => joinpath(external_dir, "missing.m")); output_root, runtime)
        @test occursin("case file not found", SparlectraApp._webui_urldecode(Dict(missing_response.headers)["Location"]))

        empty_response = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/resolve-case", Dict("casefile_manual" => ""); output_root, runtime)
        @test occursin("Enter a case name or path", SparlectraApp._webui_urldecode(Dict(empty_response.headers)["Location"]))

        # CGMES aliases typed without "cgmes:" (the five tries of a user's
        # operation log). The packed set is already in the case directory, so
        # fetchCGMESTestSet reuses it and nothing is downloaded; a typo is no
        # alias and gets the list of accepted forms, not the .m/.jl error.
        write(joinpath(case_directory, "cgmes_smallgrid.zip"), "packed set placeholder")
        for (typed, resolved) in [("cgmes smallgrind", nothing), ("cgmes_SmallGrid", "cgmes_smallgrid.zip"), ("cgmes_smallgrid.zip", "cgmes_smallgrid.zip"), ("CGMES SmallGrid.zip", "cgmes_smallgrid.zip"), ("smallgrid", "cgmes_smallgrid.zip")]
          location = SparlectraApp._webui_urldecode(Dict(SparlectraApp.route_sparlectra_webui("POST", "/powerflow/resolve-case", Dict("casefile" => typed); output_root, runtime).headers)["Location"])
          if resolved === nothing
            @test occursin("Accepted: a MATPOWER name such as case118.m", location) && occursin("cgmes:<alias>", location) && occursin("smallgrid_nb", location)
            @test !occursin("Unsupported casefile extension", location)
          else
            @test startswith(location, "/powerflow/case?casefile=$(resolved)&import_message=Resolved case: $(resolved)")
          end
        end

        log_text = read(SparlectraApp.webui_operation_log_path(output_root), String)
        @test occursin("case_resolve_completed", log_text)
        @test occursin("case_resolve_failed", log_text)
        @test !occursin("powerflow_run_started", log_text)
      end
    end)() end

    @testset "Case delete from the combobox (no PowerFlow run)" begin (function ()
      mktempdir() do tmpdir
        case_directory = joinpath(tmpdir, "cases")
        output_root = joinpath(tmpdir, "runs")
        mkpath(case_directory)
        mkpath(output_root)
        runtime = (; case_directory, config_file = "configuration.yaml", operation_log = SparlectraApp.webui_operation_log_path(output_root), startup_config_error = nothing, runner = SparlectraApp.start_powerflow_run)
        victim = joinpath(case_directory, "delete_me.m")
        write(victim, "function mpc = delete_me\nend\n")
        sidecar = SparlectraApp._webui_case_settings_path(output_root, "delete_me.m"; case_directory)
        mkpath(dirname(sidecar))
        write(sidecar, "config_version: 1\nscope: case\ncase: delete_me.m\npower_flow:\n  max_iter: 44\n")

        delete_response = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/delete-case", Dict("casefile" => "delete_me.m"); output_root, runtime)
        @test delete_response.status == 303
        @test occursin("Deleted case: delete_me.m", SparlectraApp._webui_urldecode(Dict(delete_response.headers)["Location"]))
        @test !isfile(victim)
        @test !isfile(sidecar)

        traversal_response = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/delete-case", Dict("casefile" => "../delete_me.m"); output_root, runtime)
        @test occursin("invalid case name", SparlectraApp._webui_urldecode(Dict(traversal_response.headers)["Location"]))

        missing_response = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/delete-case", Dict("casefile" => "missing.m"); output_root, runtime)
        @test occursin("case file not found", SparlectraApp._webui_urldecode(Dict(missing_response.headers)["Location"]))

        empty_response = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/delete-case", Dict{String,Any}(); output_root, runtime)
        @test occursin("Select a case to delete", SparlectraApp._webui_urldecode(Dict(empty_response.headers)["Location"]))

        log_text = read(SparlectraApp.webui_operation_log_path(output_root), String)
        @test occursin("case_delete_completed", log_text)
        @test occursin("case_delete_failed", log_text)
        @test !occursin("powerflow_run_started", log_text)
      end
    end)() end


    @testset "help registry links into the documentation" begin (function ()
      # the same check the docs gate runs: every topic has a hint of at most
      # 160 characters and every doc target is a labelled heading in docs/src
      include(joinpath(@__DIR__, "..", "tools", "check_webui_doc_links.jl"))
      # the binding is defined by the include above: reach it through
      # invokelatest so no world-age warning reaches the test output
      check_fn = Base.invokelatest(getfield, Main, :check_webui_doc_links)
      @test isempty(Base.invokelatest(check_fn; io = devnull))
      @test Set(values(SparlectraApp.WEBUI_FORM_HELP_TOPICS)) == Set(keys(SparlectraApp.WEBUI_HELP_TOPICS))
      # the sections the help pages show are generated source: every doc
      # target has one, and the committed file matches the documentation
      include(joinpath(@__DIR__, "..", "tools", "generate_webui_help_excerpts.jl"))
      excerpt_check = Base.invokelatest(getfield, Main, :check_webui_help_excerpts)
      @test isempty(Base.invokelatest(excerpt_check; io = devnull))
      @test Set(keys(SparlectraApp.WEBUI_HELP_EXCERPTS)) == Set(String(m.doc) for m in values(SparlectraApp.WEBUI_HELP_TOPICS) if !isempty(m.doc))
      @test all(!isempty(e.markdown) for (k, e) in SparlectraApp.WEBUI_HELP_EXCERPTS if startswith(k, "webui/"))
      # sections of the Web UI page ship whole, library sections lead-only
      @test all(!e.truncated for (k, e) in SparlectraApp.WEBUI_HELP_EXCERPTS if startswith(k, "webui/"))
      @test !isempty(SparlectraApp.WEBUI_HELP_MANUAL)
      @test haskey(SparlectraApp.WEBUI_HELP_MANUAL_ANCHORS, "Case formats at a glance")
      @test haskey(SparlectraApp.WEBUI_HELP_REFERENCE_ANCHORS, "Form options")
      # formulas: the Web UI has no math renderer, so every help page it
      # serves shows them converted, without a `$` or a TeX command outside
      # code (the Bad-data topic showed `$r_{n,i} = r_i / \sqrt{\Omega_{ii}}$`
      # raw); a command the converter does not cover stays visible and the
      # excerpt check names it
      @test SparlectraApp._webui_tex_html("r_{n,i} = r_i / \\sqrt{\\Omega_{ii}}") == "r<sub>n,i</sub> = r<sub>i</sub> / √Ω<sub>ii</sub>"
      leftovers = Base.invokelatest(getfield, Main, :_excerpt_math_leftovers)
      @test Base.invokelatest(leftovers, SparlectraApp._webui_render_markdown("a \$\\foo{x} \\le 1\$ b")) == ["\\foo"]
      help_section(path) = match(r"<section class=\"panel help-page rendered-markdown\">.*?</section>"s, String(SparlectraApp.route_sparlectra_webui("GET", path).body)).match
      served = Dict(path => help_section(path) for path in vcat(["/help", "/help/reference"], ["/help/" * t for t in keys(SparlectraApp.WEBUI_HELP_TOPICS)]))
      @test [path for (path, html) in served if !isempty(Base.invokelatest(leftovers, html))] == String[]
      @test occursin("<span class=\"math\">r<sub>n,i</sub> = r<sub>i</sub> / √Ω<sub>ii</sub></span>", served["/help/state_estimation.robust_mode"])
      # every visible control on the Runs, Case and Settings pages carries
      # its help link, and every help link there points at a registered
      # topic (hidden inputs such as config_file carry no link)
      mktempdir() do help_root
        pages_html = string(SparlectraApp.render_powerflow_form(output_root = help_root),
                            SparlectraApp.render_case_page(output_root = help_root),
                            SparlectraApp.render_settings_page(output_root = help_root))
        registered = Set(SparlectraApp.webui_help_page_url(t) for t in values(SparlectraApp.WEBUI_FORM_HELP_TOPICS))
        rendered = Set(String(m.captures[1]) for m in eachmatch(r"class=\"help-link\" href=\"([^\"]+)\"", pages_html))
        @test !isempty(rendered)
        @test issubset(rendered, registered)
        for (field, topic) in SparlectraApp.WEBUI_FORM_HELP_TOPICS
          controls = [m.match for m in eachmatch(Regex("<(?:input|select|textarea)\\b[^>]*\\bname=\"$(field)\"[^>]*>"), pages_html)]
          any(c -> !occursin("type=\"hidden\"", c), controls) || continue
          @test SparlectraApp.webui_help_page_url(topic) in rendered
        end
      end
    end)() end

    mktempdir() do tmpdir
      casefile = _write_webui_test_case(joinpath(tmpdir, "case_webui.m"))
      config_file = joinpath(tmpdir, "config_template.yaml")
      cp(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, config_file)
      output_root = joinpath(tmpdir, "runs")
      form = _webui_test_form(casefile, config_file, output_root)

      request = SparlectraApp.powerflow_webui_request(form; default_output_root = output_root)
      @test request["casefile"] == casefile
      existing_form = copy(form)
      existing_form["casefile"] = "  case118.jl  "
      existing_form["casefile_manual"] = "  "
      existing_request = SparlectraApp.powerflow_webui_request(existing_form; default_output_root = output_root)
      @test existing_request["casefile"] == "case118.jl"
      @test existing_request["case_format"] == "auto"

      typed_form = copy(form)
      typed_form["casefile"] = "case9.m"
      typed_form["casefile_manual"] = "  case9241pegase.m  "
      typed_request = SparlectraApp.powerflow_webui_request(typed_form; default_output_root = output_root)
      @test typed_request["casefile"] == "case9241pegase.m"
      @test typed_request["case_format"] == "auto"

      dtf_form = copy(form)
      dtf_form["casefile"] = "FOR001.DAT"
      dtf_form["casefile_manual"] = "  "
      dtf_form["case_format"] = "dtf_for001"
      dtf_request = SparlectraApp.powerflow_webui_request(dtf_form; default_output_root = output_root)
      @test dtf_request["casefile"] == "FOR001.DAT"
      @test dtf_request["case_format"] == "dtf_for001"
      @test dtf_request["for002_reference_file"] === nothing

      auto_dtf_form = copy(dtf_form)
      auto_dtf_form["case_format"] = "auto"
      auto_dtf_request = SparlectraApp.powerflow_webui_request(auto_dtf_form; default_output_root = output_root)
      @test auto_dtf_request["casefile"] == "FOR001.DAT"
      @test auto_dtf_request["case_format"] == "auto"

      manual_dtf = joinpath(tmpdir, "FOR001.DAT")
      manual_dtf_form = copy(dtf_form)
      manual_dtf_form["casefile"] = "case9.m"
      manual_dtf_form["casefile_manual"] = "  " * manual_dtf * "  "
      manual_dtf_request = SparlectraApp.powerflow_webui_request(manual_dtf_form; default_output_root = output_root)
      @test manual_dtf_request["casefile"] == manual_dtf
      @test manual_dtf_request["case_format"] == "dtf_for001"

      reference_dtf_form = copy(dtf_form)
      reference_dtf_form["for002_reference_file"] = "  /manual/FOR002.DAT  "
      reference_dtf_request = SparlectraApp.powerflow_webui_request(reference_dtf_form; default_output_root = output_root)
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
      for002_primary_response = SparlectraApp.start_powerflow_run(merge(copy(base_service_request), Dict{String,Any}("casefile" => "FOR002.DAT")); case_directory = tmpdir)
      @test !for002_primary_response["success"]
      @test for002_primary_response["reason"] == "invalid_casefile"
      @test occursin("FOR002.DAT is a reference/result file and cannot be used as the primary DTF network input case.", for002_primary_response["message"])

      manual_for002_primary_response = SparlectraApp.start_powerflow_run(merge(copy(base_service_request), Dict{String,Any}("casefile" => cache_for002)); case_directory = tmpdir)
      @test !manual_for002_primary_response["success"]
      @test manual_for002_primary_response["reason"] == "invalid_casefile"
      @test occursin("Use a runnable DTF network case as the case and enter FOR002.DAT as optional FOR002 reference file.", manual_for002_primary_response["message"])

      empty_case_form = copy(form)
      empty_case_form["casefile"] = ""
      empty_case_form["casefile_manual"] = ""
      @test_throws ArgumentError SparlectraApp.powerflow_webui_request(empty_case_form; default_output_root = output_root)
      empty_case_response = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/run", empty_case_form; output_root = output_root)
      @test empty_case_response.status == 400
      # format-neutral wording (MATPOWER is only one of the supported
      # formats), pointing at the Case page
      @test occursin("Select a case first (Case page).", String(empty_case_response.body))
      operation_log_path = SparlectraApp.webui_operation_log_path(output_root)
      @test isfile(operation_log_path)
      @test occursin("\"event\":\"validation_error\"", read(operation_log_path, String))
      @test request["config_file"] == config_file
      @test request["output_root"] == output_root
      @test request["detailed_result_csv"] === true
      # 0.15.1: the CSV format is a configuration override only, no request field
      @test !haskey(request, "detailed_result_csv_format")
      @test request["config_overrides"]["output.csv_format"] == "excel_de"
      @test request["config_overrides"]["power_flow.qlimits.enabled"] === true
      @test request["config_overrides"]["power_flow.qlimits.enforcement_mode"] == "classic_simultaneous"
      @test request["config_overrides"]["power_flow.start_current_iteration.enabled"] === true
      @test request["config_overrides"]["power_flow.start_current_iteration.max_iter"] == 7
      @test request["config_overrides"]["power_flow.start_current_iteration.tol"] == 1.0e-4
      @test request["config_overrides"]["power_flow.merit.enabled"] === true
      @test request["config_overrides"]["power_flow.merit.armijo_c1"] == 1.0e-4
      @test request["config_overrides"]["power_flow.merit.fallback_max_mismatch"] === true
      merit_disabled_form = copy(form)
      merit_disabled_form["power_flow_merit_enabled"] = "false"
      @test SparlectraApp.powerflow_webui_request(merit_disabled_form; default_output_root = output_root)["config_overrides"]["power_flow.merit.enabled"] === false
      merit_missing_form = copy(form)
      delete!(merit_missing_form, "power_flow_merit_enabled")
      @test !haskey(SparlectraApp.powerflow_webui_request(merit_missing_form; default_output_root = output_root)["config_overrides"], "power_flow.merit.enabled")
      @test request["config_overrides"]["power_flow.trust_region.enabled"] === false
      @test request["config_overrides"]["power_flow.trust_region.initial_radius"] == 1.0
      @test request["config_overrides"]["power_flow.trust_region.eta_accept"] == 0.1
      trust_region_enabled_form = copy(form)
      trust_region_enabled_form["power_flow_trust_region_enabled"] = "true"
      @test SparlectraApp.powerflow_webui_request(trust_region_enabled_form; default_output_root = output_root)["config_overrides"]["power_flow.trust_region.enabled"] === true
      trust_region_missing_form = copy(form)
      delete!(trust_region_missing_form, "power_flow_trust_region_enabled")
      @test !haskey(SparlectraApp.powerflow_webui_request(trust_region_missing_form; default_output_root = output_root)["config_overrides"], "power_flow.trust_region.enabled")
      qlimits_disabled_form = copy(form)
      qlimits_disabled_form["power_flow_qlimits_enabled"] = "false"
      @test SparlectraApp.powerflow_webui_request(qlimits_disabled_form; default_output_root = output_root)["config_overrides"]["power_flow.qlimits.enabled"] === false
      missing_field_form = copy(form)
      delete!(missing_field_form, "power_flow_start_current_iteration_enabled")
      @test !haskey(SparlectraApp.powerflow_webui_request(missing_field_form; default_output_root = output_root)["config_overrides"], "power_flow.start_current_iteration.enabled")
      ci_disabled_form = copy(form)
      ci_disabled_form["power_flow_start_current_iteration_enabled"] = "false"
      @test SparlectraApp.powerflow_webui_request(ci_disabled_form; default_output_root = output_root)["config_overrides"]["power_flow.start_current_iteration.enabled"] === false
      csv_disabled_form = copy(form)
      csv_disabled_form["detailed_result_csv"] = "false"
      @test SparlectraApp.powerflow_webui_request(csv_disabled_form; default_output_root = output_root)["detailed_result_csv"] === false
      delete!(csv_disabled_form, "detailed_result_csv_format")
      @test !haskey(SparlectraApp.powerflow_webui_request(csv_disabled_form; default_output_root = output_root)["config_overrides"], "output.csv_format")
      @test SparlectraApp.powerflow_webui_request(form; default_output_root = output_root)["diagnose_mode"] === false
      diagnose_form = copy(form)
      diagnose_form["diagnose_mode"] = "true"
      @test SparlectraApp.powerflow_webui_request(diagnose_form; default_output_root = output_root)["diagnose_mode"] === true
      # browser side of the same value: Diagnose is a second named submit
      # button, and the form disables its buttons while it submits. A
      # disabled submitter is dropped from the form data, so the page copies
      # the submitter's name and value into a hidden input first (without it
      # Diagnose silently ran a normal power flow)
      run_form_html = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow?casefile=$(SparlectraApp._webui_urlencode(casefile))&config_file=$(SparlectraApp._webui_urlencode(config_file))"; output_root = output_root).body)
      @test occursin("type=\"submit\" name=\"diagnose_mode\" value=\"true\"", run_form_html)
      @test !occursin("onsubmit=\"this.classList.add('is-submitting')", run_form_html)
      @test occursin("const copySubmitterValue = function (submitter)", run_form_html)
      @test occursin("event.submitter", run_form_html)
      @test occursin("data-submitter-value", run_form_html)
      @test occursin("staleSubmitterValue", run_form_html)
      # A normal "Start PowerFlow run" never writes diagnose.log — there is no
      # "Run diagnostics" checkbox; only diagnose_mode forces it server-side.
      @test SparlectraApp.powerflow_webui_request(form; default_output_root = output_root)["run_diagnostics"] === false
      run_diagnostics_form = copy(form)
      run_diagnostics_form["run_diagnostics"] = "on"
      @test SparlectraApp.powerflow_webui_request(run_diagnostics_form; default_output_root = output_root)["run_diagnostics"] === false
      untrusted_form = copy(form)
      untrusted_form["output_root"] = joinpath(tmpdir, "outside")
      @test SparlectraApp.powerflow_webui_request(untrusted_form; default_output_root = output_root)["output_root"] == output_root
      overrides = request["config_overrides"]
      # Every GUI-editable field arrives as an override, EXCEPT the ones the
      # form left empty: an empty numeric field means "not stated" and must
      # not send a value (that is how the optional tolerance in MW works).
      stated = Set(key for (key, field, _) in SparlectraApp._WEBUI_FORM_CONFIG_FIELDS
                   if !isempty(strip(String(get(form, field, "")))))
      # the run form's CSV format is the run's central output.csv_format
      # for every run type (0.14.2), carried as an override like any key
      isempty(strip(String(get(form, "detailed_result_csv_format", "")))) || push!(stated, "output.csv_format")
      @test Set(keys(overrides)) == stated
      @test !haskey(overrides, "power_flow.tol_MW")
      @test overrides["power_flow.tol"] == 1.0e-8
      @test overrides["power_flow.max_iter"] == 80
      @test overrides["power_flow.autodamp"]
      # the advanced auto_slack checkbox reaches the run as an override (0.20.4)
      @test overrides["power_flow.auto_slack"] === true
      apply_form = copy(form)
      apply_form["matpower_import_auto_profile"] = "apply"
      apply_request = SparlectraApp.powerflow_webui_request(apply_form; default_output_root = output_root)
      @test apply_request["config_overrides"]["model.auto_profile"] == "apply"
      @test !haskey(overrides, "benchmark.enabled")   # no benchmark option in the Web UI

      invalid_form = copy(form)
      invalid_form["power_flow_max_iter"] = "invalid"
      invalid_response = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/run", invalid_form; output_root = output_root)
      invalid_html = String(invalid_response.body)
      @test invalid_response.status == 400
      # the form comes back with the entered values
      @test occursin("name=\"config_file\" value=\"$(config_file)\"", invalid_html)
      last_errors_response = SparlectraApp.route_sparlectra_webui("GET", "/webui/last-errors"; output_root = output_root)
      last_errors_html = String(last_errors_response.body)
      @test last_errors_response.status == 200
      @test occursin("Last errors", last_errors_html)
      @test occursin("validation_error", read(operation_log_path, String))
      @test occursin("/powerflow/run", last_errors_html)
      empty_errors_html = SparlectraApp.render_webui_last_errors(joinpath(tmpdir, "missing-operation-log.jsonl"))
      @test occursin("No recent errors.", empty_errors_html)

      api_failure_root = joinpath(tmpdir, "api-failure-runs")
      api_failure_case = _write_webui_test_case(joinpath(tmpdir, "case_with_dcline.m"))
      api_failure_form = _webui_test_form(api_failure_case, config_file, api_failure_root)
      api_failure_message = "Unsupported MATPOWER dcline data: active mpc.dcline entries detected. DC lines are currently not supported; the run was aborted before power-flow solution."
      api_failure_runner = function(worker_request; case_directory = nothing)
        api_failure_run_id = worker_request["run_id"]
        api_failure_output_dir = joinpath(worker_request["output_root"], api_failure_run_id)
        mkpath(api_failure_output_dir)
        result = SparlectraApp._api_result(
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
        SparlectraApp._POWERFLOW_SERVICE_RUNS[api_failure_run_id] = result
        return SparlectraApp.to_dict(result)
      end
      api_failure_response = SparlectraApp.handle_powerflow_run(api_failure_form; default_output_root = api_failure_root, runner = api_failure_runner, operation_log = api_failure_root)
      api_failure_run_id = api_failure_response["run_id"]
      wait(SparlectraApp._POWERFLOW_WEBUI_JOBS[api_failure_run_id]["task"])
      api_failure_job = SparlectraApp.get_webui_powerflow_job(api_failure_run_id)
      @test api_failure_job["status"] == "failed"
      @test api_failure_job["run_status"] == "failed"
      @test api_failure_job["reason"] == "unsupported_matpower_dcline"
      api_failure_result_html = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow/result/$(api_failure_run_id)"; output_root = api_failure_root).body)
      @test occursin(api_failure_message, api_failure_result_html)
      api_failure_last_errors_html = String(SparlectraApp.route_sparlectra_webui("GET", "/webui/last-errors"; output_root = api_failure_root).body)
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
        result = SparlectraApp._api_result(
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
        SparlectraApp._POWERFLOW_SERVICE_RUNS[api_success_run_id] = result
        return SparlectraApp.to_dict(result)
      end
      api_success_response = SparlectraApp.handle_powerflow_run(api_success_form; default_output_root = api_success_root, runner = api_success_runner, operation_log = api_success_root)
      api_success_run_id = api_success_response["run_id"]
      wait(SparlectraApp._POWERFLOW_WEBUI_JOBS[api_success_run_id]["task"])
      api_success_last_errors_html = String(SparlectraApp.route_sparlectra_webui("GET", "/webui/last-errors"; output_root = api_success_root).body)
      @test occursin("No recent errors.", api_success_last_errors_html)

      logo_response = SparlectraApp.route_sparlectra_webui("GET", "/assets/logo.png")
      @test logo_response.status == 200
      @test ("Content-Type" => "image/png") in logo_response.headers
      @test logo_response.body == read(joinpath(@__DIR__, "..", "docs", "src", "assets", "logo.png"))
      for target in ("/assets/tablestyle.css", "/assets/../Project.toml", "/assets/logo.png/extra")
        @test SparlectraApp.route_sparlectra_webui("GET", target).status == 404
      end
      routed_form = copy(form)
      routed_form["output_root"] = joinpath(tmpdir, "outside-runs")
      run_response = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/run", routed_form; output_root)
      @test run_response.status == 303
      result_location = only(header.second for header in run_response.headers if header.first == "Location")
      run_id = basename(result_location)
      wait(SparlectraApp._POWERFLOW_WEBUI_JOBS[run_id]["task"])
      result = get_powerflow_result(run_id)
      @test result["success"]
      @test result["output_dir"] == joinpath(abspath(output_root), run_id)
      @test !ispath(joinpath(tmpdir, "outside-runs"))
      run_log_text = read(joinpath(result["output_dir"], "run.log"), String)
      @test !occursin("Unknown Sparlectra configuration key: power_flow.qlimits.enforcement_mode", run_log_text)
      @test occursin("Total time  :", run_log_text)
      @test occursin("Output time :", run_log_text)
      @test occursin("Solver time :", run_log_text)
      result_json = SparlectraApp._parse_service_json(read(joinpath(result["output_dir"], "result.json"), String))
      @test result_json["run_id"] == run_id
      @test result_json["schema_version"] == "1.0"
      @test haskey(result_json, "service_phase_timings")
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
      @test occursin("\"event\":\"performance_timing_enabled\"", operation_log_text)
      @test occursin("\"event\":\"detailed_result_csv_export_enabled\"", operation_log_text)
      @test occursin("\"csv_format\":\"excel_de\"", operation_log_text)
      @test occursin("\"delimiter\":\";\"", operation_log_text)
      @test occursin("\"decimal_separator\":\",\"", operation_log_text)
      @test occursin("\"thousands_separator\":\".\"", operation_log_text)
      @test occursin("\"event\":\"detailed_result_csv_exported\"", operation_log_text)
      @test occursin("\"artifacts\":[\"bus_voltages_complex.csv\",\"branch_flows.csv\"]", operation_log_text)
      @test isfile(joinpath(result["output_dir"], "performance.log"))
      @test isfile(joinpath(result["output_dir"], "q_limit.log"))
      manual_import_form = copy(form)
      manual_import_form["matpower_import_auto_profile"] = "off"
      manual_import_form["matpower_import_ratio"] = "reciprocal"
      manual_import_form["matpower_import_shift_sign"] = "-1.0"
      manual_import_form["matpower_import_shift_unit"] = "rad"
      manual_response = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/run", manual_import_form; output_root)
      @test manual_response.status == 303
      manual_location = only(header.second for header in manual_response.headers if header.first == "Location")
      manual_run_id = basename(manual_location)
      wait(SparlectraApp._POWERFLOW_WEBUI_JOBS[manual_run_id]["task"])
      manual_result = get_powerflow_result(manual_run_id)
      @test manual_result["success"]
      manual_effective_cfg = Sparlectra.load_sparlectra_config(joinpath(manual_result["output_dir"], "effective_config.yaml"); reload = true)
      @test manual_effective_cfg.model.auto_profile === :off
      @test manual_effective_cfg.matpower.ratio === :reciprocal
      @test manual_effective_cfg.matpower.shift_sign == -1.0
      @test manual_effective_cfg.matpower.shift_unit === :rad
      manual_run_log = read(joinpath(manual_result["output_dir"], "run.log"), String)
      # compact console (the default): the conventions that
      # move results are named on one line, the full list stays in
      # effective_config.yaml, which the assertions above read
      @test occursin("Import conventions: ratio=reciprocal", manual_run_log)
      @test occursin("shift=rad", manual_run_log)
      # dishonest Newton (0.30.2): silent when off (the runs above), one
      # settings line and the counts when on; the form values reach
      # effective_config.yaml, survive "save case settings" and come back
      # in the reloaded Settings form
      @test !occursin("Jacobian reuse", manual_run_log)
      @test !haskey(manual_result["metadata"], "jacobian_reuse_reused_steps")
      # power-mode LU (0.30.2): silent without power mode, and the greyed
      # (not submitted) select leaves the configured default in place
      @test !occursin("Power-mode LU", manual_run_log)
      @test !haskey(manual_result["metadata"], "power_mode_lu_choice")
      @test manual_effective_cfg.powerflow.power_mode_lu === :auto
      reuse_form = copy(form)
      reuse_form["power_flow_jacobian_reuse"] = "true"
      reuse_form["power_flow_jacobian_reuse_min_reduction"] = "5"
      reuse_form["power_flow_jacobian_reuse_max_steps"] = "2"
      # the test form runs linear_solver umfpack, which keeps no factorisation
      reuse_form["power_flow_linear_solver"] = "umfpack_reuse"
      # the same run with power mode and a fixed power-mode LU: the choice
      # reaches the solver, run.log and the metadata, and survives the save
      reuse_form["power_flow_power_mode"] = "true"
      reuse_form["power_flow_power_mode_lu"] = "umfpack"
      reuse_response = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/run", reuse_form; output_root)
      reuse_run_id = basename(only(header.second for header in reuse_response.headers if header.first == "Location"))
      wait(SparlectraApp._POWERFLOW_WEBUI_JOBS[reuse_run_id]["task"])
      reuse_result = get_powerflow_result(reuse_run_id)
      @test reuse_result["success"]
      reuse_cfg = Sparlectra.load_sparlectra_config(joinpath(reuse_result["output_dir"], "effective_config.yaml"); reload = true)
      @test (reuse_cfg.powerflow.jacobian_reuse, reuse_cfg.powerflow.jacobian_reuse_min_reduction, reuse_cfg.powerflow.jacobian_reuse_max_steps) == (true, 5.0, 2)
      reuse_run_log = read(joinpath(reuse_result["output_dir"], "run.log"), String)
      @test occursin("Jacobian reuse: on, min reduction 5, max 2 steps", reuse_run_log)
      # on umfpack the service says that the switch cannot act
      umfpack_pf = Sparlectra.SparlectraConfig(Dict{String,Any}("config_version" => 1, "power_flow" => Dict{String,Any}("jacobian_reuse" => true, "linear_solver" => "umfpack"))).powerflow
      @test SparlectraApp._jacobian_reuse_runlog_line(umfpack_pf) == "Jacobian reuse ignored: linear_solver umfpack keeps no factorisation; use umfpack_reuse or KLU."
      @test SparlectraApp._jacobian_reuse_runlog_line(Sparlectra.SparlectraConfig().powerflow) === nothing
      @test occursin(r"Jacobian reuse\s*: \d+ reused steps, \d+ refactorisations", reuse_run_log)
      @test haskey(reuse_result["metadata"], "jacobian_reuse_reused_steps") && haskey(reuse_result["metadata"], "jacobian_reuse_refactorisations")
      @test reuse_cfg.powerflow.power_mode && reuse_cfg.powerflow.power_mode_lu === :umfpack
      @test occursin("Power-mode LU  : umfpack (power_mode_lu = umfpack)", reuse_run_log)
      @test reuse_result["metadata"]["power_mode_lu_choice"] == "umfpack" && reuse_result["metadata"]["power_mode_lu_source"] == "configured"
      reuse_save = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/result/$(reuse_run_id)/case-settings/save", Dict{String,String}(); output_root)
      @test reuse_save.status == 200
      reuse_casefile = form["casefile_manual"]
      reuse_cc = Sparlectra.load_case_config(reuse_casefile)
      @test reuse_cc["power_flow.jacobian_reuse"] == true
      reuse_reloaded = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow/settings?casefile=$(SparlectraApp._webui_urlencode(reuse_casefile))&case_settings=1"; output_root).body)
      _webui_assert_checked(reuse_reloaded, "power_flow_jacobian_reuse", true)
      _webui_assert_value(reuse_reloaded, "power_flow_jacobian_reuse_min_reduction", "5.0")
      _webui_assert_value(reuse_reloaded, "power_flow_jacobian_reuse_max_steps", "2")
      @test reuse_cc["power_flow.power_mode_lu"] == "umfpack"
      _webui_assert_selected(reuse_reloaded, "power_flow_power_mode_lu", "umfpack")
      # N-1 screening from the Settings form reaches the effective configuration
      screen_form = copy(form)
      screen_form["contingency_screening_mode"] = "flag"
      screen_form["contingency_screening_margin_pct"] = "15"
      screen_form["contingency_warm_active_set"] = "true"
      screen_form["contingency_max_iter"] = "12"
      # the checked box (its value after the hidden "false") switches the
      # cold check on against its default off
      screen_form["contingency_warm_cold_check"] = "true"
      screen_request = SparlectraApp.powerflow_webui_request(screen_form; default_output_root = output_root)
      screen_cfg, _ = Sparlectra._load_api_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, Sparlectra.validate_gui_config_overrides(screen_request["config_overrides"]))
      @test (screen_cfg.contingency.screening_mode, screen_cfg.contingency.screening_margin_pct, screen_cfg.contingency.warm_active_set, screen_cfg.contingency.warm_cold_check, screen_cfg.contingency.max_iter) == (:flag, 15.0, true, true, 12)
      # correlated events and pairs (#355) from the same form
      event_form = copy(form)
      event_form["contingency_events_enabled"] = "true"
      event_form["contingency_events_busbars"] = "true"
      event_form["contingency_events_name_groups"] = "true"
      event_form["contingency_events_max_members"] = "4"
      event_form["contingency_pairs_enabled"] = "true"
      event_form["contingency_pairs_max_pairs"] = "3000"
      event_request = SparlectraApp.powerflow_webui_request(event_form; default_output_root = output_root)
      event_cfg, _ = Sparlectra._load_api_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, Sparlectra.validate_gui_config_overrides(event_request["config_overrides"]))
      ev = event_cfg.contingency.events
      @test (ev.enabled, ev.double_circuits, ev.busbars, ev.substations, ev.links, ev.name_groups, ev.max_members) == (true, true, true, false, false, true, 4)
      @test (event_cfg.contingency.pairs.enabled, event_cfg.contingency.pairs.rule, event_cfg.contingency.pairs.max_pairs) == (true, :all, 3000)
      @test ev.name_stem_pattern == Sparlectra.ContingencyEventsConfig().name_stem_pattern
      # the pattern survives a write of the configuration file (a Settings
      # save rewrites the whole file)
      rt_yaml = joinpath(mktempdir(), "rt.yaml")
      SparlectraApp._write_yaml_file(rt_yaml, Sparlectra.load_yaml_dict(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH))
      @test Sparlectra.load_sparlectra_config(rt_yaml; reload = true).contingency.events.name_stem_pattern == Sparlectra.ContingencyEventsConfig().name_stem_pattern
      # stopping the Web UI aborts and awaits a running job before the
      # process may exit (native factorizations in use by worker threads
      # crashed the process on stop): a runner that solves until the abort
      # token is set ends within the wait, the job ends "aborted"
      stop_runner = function (worker_request; case_directory = nothing)
        while true
          Sparlectra.sparlectra_check_abort()
          sleep(0.05)
        end
      end
      stop_job = SparlectraApp.start_webui_powerflow_run(SparlectraApp.powerflow_webui_request(copy(form); default_output_root = output_root); runner = stop_runner)
      sleep(0.3)
      stopped = SparlectraApp._webui_abort_active_jobs_and_wait!(; timeout_s = 10.0)
      @test stopped.active >= 1 && stopped.still_running == 0
      @test get(SparlectraApp.get_webui_powerflow_job(stop_job["run_id"]), "status", "") == "aborted"
      # the saved case settings must not leak into the runs that follow
      rm(SparlectraApp._webui_case_settings_path(output_root, reuse_casefile); force = true)
      @test !isempty(run_id)
      result_response = SparlectraApp.handle_powerflow_result(run_id)
      result_html = String(result_response.body)
      @test result_response.status == 200
      for field in ("run_id", "status", "iterations", "final_mismatch")
        @test occursin(field, result_html)
      end
      @test occursin(run_id, result_html)
      @test !occursin("action=\"/powerflow/abort/$(run_id)\"", result_html)
      failed_result_html = SparlectraApp.render_powerflow_result(Dict("run_id" => "failed-run", "status" => "failed", "success" => false))
      @test !occursin("action=\"/powerflow/abort/failed-run\"", failed_result_html)
      converged_result_html = SparlectraApp.render_powerflow_result(Dict(
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
      @test occursin("converged", converged_result_html)
      @test occursin("numerical_converged</th><td>true</td>", converged_result_html)
      @test occursin("solution_available</th><td>true</td>", converged_result_html)
      @test occursin("iterations</th><td>6</td>", converged_result_html)
      @test occursin("final_mismatch</th><td>6.3e-11</td>", converged_result_html)
      nonconverged_result_html = SparlectraApp.render_powerflow_result(Dict(
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
      missing_metric_html = SparlectraApp.render_powerflow_result(Dict("run_id" => "missing-metrics", "status" => "success", "success" => true))
      for field in ("converged", "numerical_converged", "solution_available", "iterations", "final_mismatch", "reason")
        @test occursin("$(field)</th><td>n/a</td>", missing_metric_html)
      end
      for active_status in ("queued", "running")
        active_status_html = SparlectraApp.render_powerflow_result(Dict("run_id" => "$(active_status)-run", "status" => active_status))
        @test occursin("data-refresh-url=\"/powerflow/result/$(active_status)-run?autorefresh=1\"", active_status_html)
        @test !occursin("http-equiv=\"refresh\"", active_status_html)
        @test occursin("action=\"/powerflow/abort/$(active_status)-run\"", active_status_html)
      end
      for terminal_status in ("success", "failed", "aborted")
        terminal_status_html = SparlectraApp.render_powerflow_result(Dict("run_id" => "$(terminal_status)-run", "status" => terminal_status, "elapsed_seconds" => 83))
        @test !occursin("data-refresh-url=\"", terminal_status_html)
        @test !occursin("http-equiv=\"refresh\"", terminal_status_html)
        @test occursin("class=\"runtime-card\"", terminal_status_html)
        @test occursin("Total time", terminal_status_html)
        @test occursin("00:01:23.000", terminal_status_html)
      end

      artifacts = list_powerflow_artifacts(run_id)
      artifact_names = Set(artifact["name"] for artifact in artifacts)
      @test Set(("result.json", "run.log", "effective_config.yaml", "run_metadata.yaml", "performance.log", "q_limit.log", "bus_voltages_complex.csv", "branch_flows.csv")) ⊆ artifact_names
      @test !("diagnose.log" in artifact_names) # normal "Start PowerFlow run" never writes diagnose.log
      artifact_response = SparlectraApp.handle_powerflow_artifacts(run_id)
      artifact_html = String(artifact_response.body)
      for name in ("result.json", "run.log", "effective_config.yaml", "run_metadata.yaml", "performance.log", "q_limit.log", "bus_voltages_complex.csv", "branch_flows.csv")
        @test occursin(name, artifact_html)
      end
      qlimit_artifact = only(artifact for artifact in artifacts if artifact["name"] == "q_limit.log")
      @test qlimit_artifact["kind"] == "q_limit_log"
      @test qlimit_artifact["mime_type"] == "text/plain"
      qlimit_download = SparlectraApp.handle_powerflow_artifact_download(run_id, "q_limit.log")
      @test qlimit_download.status == 200
      @test occursin("Resolved Q-limit options", String(qlimit_download.body))
      @test SparlectraApp.route_sparlectra_webui("GET", "/powerflow/artifact/$(run_id)/q_limit.log"; output_root).status == 200
      for name in ("bus_voltages_complex.csv", "branch_flows.csv")
        csv_artifact = only(artifact for artifact in artifacts if artifact["name"] == name)
        @test csv_artifact["kind"] == "csv"
        @test csv_artifact["mime_type"] == "text/csv"
        csv_download = SparlectraApp.handle_powerflow_artifact_download(run_id, name)
        @test csv_download.status == 200
        @test ("Content-Type" => "text/csv") in csv_download.headers
        expected_header = name == "bus_voltages_complex.csv" ? "bus;bus_name;type;vm_pu;va_deg" : "branch;branch_index;from_bus;to_bus;status"
        @test startswith(String(csv_download.body), expected_header)
        # the viewer opens a CSV as a table (semicolon file: the delimiter is sniffed), raw=1 as text
        csv_page_html = String(SparlectraApp.handle_powerflow_artifact(run_id, name).body)
        @test occursin("class=\"artifact-csv\"", csv_page_html)
        # a column with numbers has a sort link in its header
        sort_column, sort_index = name == "bus_voltages_complex.csv" ? ("vm_pu", 4) : ("from_bus", 3)
        @test occursin("<th class=\"num\" data-col=\"$(sort_index)\"><a href=\"?sort=$(sort_column)\" title=\"Sort largest first\" data-csv-keep>$(sort_column)</a></th>", csv_page_html)
        @test occursin("<td class=\"num\">", csv_page_html)
        @test !occursin("class=\"artifact-text\"", csv_page_html)
        @test occursin("href=\"$(name)?raw=1\"", csv_page_html)
        csv_raw_html = String(SparlectraApp.handle_powerflow_artifact(run_id, name; raw = true).body)
        @test occursin("class=\"artifact-text\"", csv_raw_html)
        @test occursin(expected_header, csv_raw_html)
        @test !occursin("class=\"artifact-csv\"", csv_raw_html)
        @test occursin("class=\"artifact-text\"", String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow/artifact/$(run_id)/$(name)?raw=1"; output_root).body))
      end
      # a Markdown report renders (heading, table), HTML inside it is escaped, the mime type names the format
      markdown_path = joinpath(result["output_dir"], "note.md")
      write(markdown_path, "# Heading\n\n| bus | vm |\n|---|---|\n| 1 | 1.02 |\n\nText with <script>alert(1)</script>.\n")
      markdown_artifact = only(artifact for artifact in list_powerflow_artifacts(run_id) if artifact["name"] == "note.md")
      @test markdown_artifact["mime_type"] == "text/markdown"
      markdown_html = String(SparlectraApp.handle_powerflow_artifact(run_id, "note.md").body)
      @test occursin("class=\"rendered-markdown\"", markdown_html)
      @test occursin("<h1>Heading</h1>", markdown_html)
      @test occursin("<table>", markdown_html)
      @test occursin("1.02", markdown_html)
      @test !occursin("<script>alert", markdown_html)
      @test occursin("&lt;script&gt;", markdown_html)
      @test occursin("href=\"note.md?raw=1\"", markdown_html)
      markdown_raw_html = String(SparlectraApp.handle_powerflow_artifact(run_id, "note.md"; raw = true).body)
      @test occursin("class=\"artifact-text\"", markdown_raw_html)
      @test !occursin("<h1>Heading</h1>", markdown_raw_html)
      # a comment header (measurement CSV v1) sits above the table, a truncated preview drops its cut row
      comment_csv_path = joinpath(result["output_dir"], "note_rows.csv")
      write(comment_csv_path, "# version: 1\nbus,value\n1,0.5\n2,\"1,5\"\n")
      comment_csv_html = String(SparlectraApp.handle_powerflow_artifact(run_id, "note_rows.csv").body)
      @test occursin("class=\"artifact-csv-comments\"># version: 1</pre>", comment_csv_html)
      @test occursin("<th class=\"num\" data-col=\"1\"><a href=\"?sort=bus\" title=\"Sort largest first\" data-csv-keep>bus</a></th><th class=\"num\" data-col=\"2\"><a href=\"?sort=value\" title=\"Sort largest first\" data-csv-keep>value</a></th>", comment_csv_html)
      @test occursin("<td class=\"num\">1,5</td>", comment_csv_html)
      truncated_table = SparlectraApp._webui_artifact_csv_html("a,b\n1,2\n3,4\n5,"; truncated = true)
      @test occursin("<td class=\"num\">3</td>", truncated_table)
      @test !occursin("<td class=\"num\">5</td>", truncated_table)
      # number display of the CSV table: significant digits (default 5), power
      # columns in whole W/kW/MW, pu and integer columns never converted
      tiny = "1.98394650844228e-16"
      display_csv = "bus,p_mw,vm_pu\n7,$(tiny),$(tiny)\n8,123.456789,123.456789\n"
      default_table = SparlectraApp._webui_artifact_csv_html(display_csv)
      @test occursin("<td class=\"num\">123.46</td><td class=\"num\">123.46</td>", default_table)
      @test occursin("<td class=\"num\">1.9839e-16</td><td class=\"num\">1.9839e-16</td>", default_table)
      @test occursin("<td class=\"num\">7</td>", default_table)
      four_digits = SparlectraApp._webui_artifact_csv_html(display_csv; csv_display = SparlectraApp._WebUICsvDisplay(4, ""))
      @test occursin("<td class=\"num\">123.5</td><td class=\"num\">123.5</td>", four_digits)
      kw_table = SparlectraApp._webui_artifact_csv_html(display_csv; csv_display = SparlectraApp._WebUICsvDisplay(3, "kW"))
      @test occursin("<td class=\"num\">7</td><td class=\"num\">0</td><td class=\"num\">1.98e-16</td>", kw_table)
      @test occursin("<td class=\"num\">123457</td><td class=\"num\">123</td>", kw_table)
      @test occursin("<th class=\"num\" data-col=\"2\"><a href=\"?digits=3&amp;power=kW&amp;sort=p_mw\" title=\"Sort largest first\" data-csv-keep>p_mw [kW]</a></th><th class=\"num\" data-col=\"3\"><a href=\"?digits=3&amp;power=kW&amp;sort=vm_pu\" title=\"Sort largest first\" data-csv-keep>vm_pu</a></th>", kw_table)
      for unit in ("W", "MW")
        @test occursin("<td class=\"num\">7</td><td class=\"num\">0</td>", SparlectraApp._webui_artifact_csv_html(display_csv; csv_display = SparlectraApp._WebUICsvDisplay(3, unit)))
      end
      # a semicolon file reads and shows with the decimal comma, thousands in groups of three
      comma_csv = "bus;q_mvar;vm_pu\n1;0,000000000000000198394650844228;1,0234567\n2;1.234,5678;0,5\n"
      comma_table = SparlectraApp._webui_artifact_csv_html(comma_csv)
      @test occursin("<td class=\"num\">1,9839e-16</td><td class=\"num\">1,0235</td>", comma_table)
      @test occursin("<td class=\"num\">1.234,6</td><td class=\"num\">0,50000</td>", comma_table)
      comma_kw = SparlectraApp._webui_artifact_csv_html(comma_csv; csv_display = SparlectraApp._WebUICsvDisplay(3, "kW"))
      @test occursin("<td class=\"num\">0</td><td class=\"num\">1,02</td>", comma_kw)
      @test occursin("<td class=\"num\">1.234.568</td>", comma_kw)
      @test occursin("<th class=\"num\" data-col=\"2\"><a href=\"?digits=3&amp;power=kW&amp;sort=q_mvar\" title=\"Sort largest first\" data-csv-keep>q_mvar [kvar]</a></th>", comma_kw)
      # the choice comes from the query, stays on the links to the other CSV
      # files and the raw view, and never changes the file or the download
      power_csv_text = "bus;p_mw\n1;123,456789\n"
      write(joinpath(result["output_dir"], "note_power.csv"), power_csv_text)
      routed = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow/artifact/$(run_id)/note_power.csv?digits=4&power=kW"; output_root).body)
      @test occursin("<td class=\"num\">123.457</td>", routed)
      @test occursin("href=\"bus_voltages_complex.csv?digits=4&amp;power=kW\"", routed)
      @test occursin("href=\"note_power.csv?raw=1&amp;digits=4&amp;power=kW\"", routed)
      routed_bad = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow/artifact/$(run_id)/note_power.csv?digits=abc&power=GW"; output_root).body)
      @test occursin("<td class=\"num\">123,46</td>", routed_bad)
      @test occursin("href=\"bus_voltages_complex.csv\"", routed_bad)
      @test SparlectraApp._webui_csv_display(Dict("digits" => "0", "power" => "kw")) == SparlectraApp._WebUICsvDisplay(5, "")
      routed_raw = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow/artifact/$(run_id)/note_power.csv?raw=1&digits=4&power=kW"; output_root).body)
      @test occursin("123,456789", routed_raw)
      @test String(SparlectraApp.handle_powerflow_artifact_download(run_id, "note_power.csv").body) == power_csv_text
      # sorting by a column with numbers: largest first, then smallest
      # first; empty and text cells at the end in file order; an unknown
      # sort column keeps the file order. Every row starts with its row
      # number of the file, which stays with the row through the sort.
      sort_csv = "name,p_mw,vm_pu\nA,-5.5,1.01\nB,12.25,0.99\nC,,1.0\nD,3.0,x\n"
      row_order(h) = [m.captures[1] for m in eachmatch(r"<tr><td class=\"csv-rn\">\d+</td><td>([A-D])</td>", h)]
      row_numbers(h) = [m.captures[1] for m in eachmatch(r"<tr><td class=\"csv-rn\">(\d+)</td>", h)]
      @test row_order(SparlectraApp._webui_artifact_csv_html(sort_csv)) == ["A", "B", "C", "D"]
      sorted_desc = SparlectraApp._webui_artifact_csv_html(sort_csv; csv_display = SparlectraApp._WebUICsvDisplay(5, "", "p_mw", "desc"))
      @test row_order(sorted_desc) == ["B", "D", "A", "C"]
      @test row_numbers(sorted_desc) == ["2", "4", "1", "3"]
      @test occursin("<th class=\"num\" data-col=\"2\"><a href=\"?sort=p_mw&amp;order=asc\" title=\"Sort smallest first\" data-csv-keep>p_mw<span class=\"csv-sort-mark\">&#9660;</span></a></th>", sorted_desc)
      @test occursin("<th class=\"num\" data-col=\"3\"><a href=\"?sort=vm_pu\" title=\"Sort largest first\" data-csv-keep>vm_pu</a></th>", sorted_desc)
      @test occursin("Sorted by p_mw, largest first", sorted_desc)
      @test row_order(SparlectraApp._webui_artifact_csv_html(sort_csv; csv_display = SparlectraApp._WebUICsvDisplay(5, "", "p_mw", "asc"))) == ["A", "D", "B", "C"]
      @test row_order(SparlectraApp._webui_artifact_csv_html(sort_csv; csv_display = SparlectraApp._WebUICsvDisplay(5, "", "vm_pu", "desc"))) == ["A", "C", "B", "D"]
      unsorted = SparlectraApp._webui_artifact_csv_html(sort_csv; csv_display = SparlectraApp._WebUICsvDisplay(5, "", "nope", "desc"))
      @test row_order(unsorted) == ["A", "B", "C", "D"] && !occursin("Sorted by", unsorted)
      # a text column sorts alphabetically ignoring case, A to Z on the first
      # click, Z to A on the second, empty cells at the end
      text_csv = "name,p_mw\nb,1.0\nC,2.0\n,3.0\na,4.0\n"
      text_plain = SparlectraApp._webui_artifact_csv_html(text_csv)
      @test occursin("<th data-col=\"1\"><a href=\"?sort=name&amp;order=asc\" title=\"Sort A to Z\" data-csv-keep>name</a></th>", text_plain)
      text_asc = SparlectraApp._webui_artifact_csv_html(text_csv; csv_display = SparlectraApp._WebUICsvDisplay(5, "", "name", "asc"))
      @test row_numbers(text_asc) == ["4", "1", "2", "3"]
      @test occursin("Sorted by name, A to Z", text_asc)
      @test occursin("<a href=\"?sort=name\" title=\"Sort Z to A\" data-csv-keep>name<span class=\"csv-sort-mark\">&#9650;</span></a>", text_asc)
      @test row_numbers(SparlectraApp._webui_artifact_csv_html(text_csv; csv_display = SparlectraApp._WebUICsvDisplay(5, "", "name", "desc"))) == ["2", "1", "4", "3"]
      # the header of a column is right-aligned like its cells: when most of
      # its non-empty cells are numbers (vm_pu with one NaN), else left
      align_csv = "name,vm_pu,note\nA,1.0,x\nB,NaN,y\nC,0.9,1\n"
      align_html = SparlectraApp._webui_artifact_csv_html(align_csv)
      @test occursin("<th class=\"num\" data-col=\"2\">", align_html)
      @test occursin("<th data-col=\"1\">", align_html) && occursin("<th data-col=\"3\">", align_html)
      @test occursin("<td>A</td><td class=\"num\">1.0000</td><td>x</td>", align_html)
      # the decimal comma of a semicolon file sorts by value
      comma_sorted = SparlectraApp._webui_artifact_csv_html("bus;q_mvar\n1;0,5\n2;1.234,5\n3;-2,0\n"; csv_display = SparlectraApp._WebUICsvDisplay(5, "", "q_mvar", "desc"))
      @test row_numbers(comma_sorted) == ["2", "1", "3"]
      # every row is sorted before the table is cut to its row limit
      limit = SparlectraApp._WEBUI_ARTIFACT_TABLE_ROWS
      long_csv = "id,p_mw\n" * join(("$(i),$(i == limit + 50 ? 999.5 : 1.5)" for i in 1:(limit + 50)), "\n") * "\n"
      long_sorted = SparlectraApp._webui_artifact_csv_html(long_csv; csv_display = SparlectraApp._WebUICsvDisplay(5, "", "p_mw", "desc"))
      @test occursin("<tbody><tr><td class=\"csv-rn\">$(limit + 50)</td><td class=\"num\">$(limit + 50)</td><td class=\"num\">999.50</td>", long_sorted)
      @test occursin("Table shows the first $(limit) of $(limit + 50) rows", long_sorted)
      @test SparlectraApp._webui_csv_display(Dict("sort" => " p_mw ", "order" => "up")) == SparlectraApp._WebUICsvDisplay(5, "", "p_mw", "desc")
      # hidden columns: column numbers from `hide=`, bad entries dropped; a
      # style rule hides them (column j is cell j + 1 after the row number),
      # a number beyond the last column is ignored
      @test SparlectraApp._webui_csv_display(Dict("hide" => "3,x,1,3,0")).hide == "1,3"
      hidden_html = SparlectraApp._webui_artifact_csv_html(sort_csv; csv_display = SparlectraApp._webui_csv_display(Dict("hide" => "1,3,9")))
      @test occursin("<style data-csv-grid-style>.csv-grid th:nth-child(2),.csv-grid td:nth-child(2){display:none}.csv-grid th:nth-child(4),.csv-grid td:nth-child(4){display:none}</style>", hidden_html)
      @test occursin("<style data-csv-grid-style></style>", SparlectraApp._webui_artifact_csv_html(sort_csv))
      # through the route: Apply keeps the sort column and the hidden columns
      # as hidden fields, the sort links and the raw view carry them, the
      # links to the other CSV files do not
      write(joinpath(result["output_dir"], "note_sort.csv"), sort_csv)
      routed_sort = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow/artifact/$(run_id)/note_sort.csv?digits=4&sort=p_mw&order=asc"; output_root).body)
      @test row_order(routed_sort) == ["A", "D", "B", "C"]
      @test occursin("<input type=\"hidden\" name=\"sort\" value=\"p_mw\"><input type=\"hidden\" name=\"order\" value=\"asc\">", routed_sort)
      @test occursin("<input type=\"hidden\" name=\"hide\" value=\"\" data-csv-hide-input disabled>", routed_sort)
      @test occursin("href=\"note_sort.csv?raw=1&amp;digits=4&amp;sort=p_mw&amp;order=asc\"", routed_sort)
      @test occursin("href=\"bus_voltages_complex.csv?digits=4\"", routed_sort)
      @test occursin("<th class=\"num\" data-col=\"2\"><a href=\"?digits=4&amp;sort=p_mw\" title=\"Sort largest first\" data-csv-keep>p_mw<span class=\"csv-sort-mark\">&#9650;</span></a></th>", routed_sort)
      routed_hidden = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow/artifact/$(run_id)/note_sort.csv?sort=p_mw&hide=3%2C1"; output_root).body)
      @test row_order(routed_hidden) == ["B", "D", "A", "C"]
      @test occursin("<a href=\"?sort=p_mw&amp;order=asc&amp;hide=1%2C3\" title=\"Sort smallest first\" data-csv-keep>", routed_hidden)
      @test occursin("<a href=\"?sort=vm_pu&amp;hide=1%2C3\" title=\"Sort largest first\" data-csv-keep>", routed_hidden)
      @test occursin("<input type=\"hidden\" name=\"hide\" value=\"1,3\" data-csv-hide-input>", routed_hidden)
      @test occursin("href=\"note_sort.csv?raw=1&amp;sort=p_mw&amp;hide=1%2C3\"", routed_hidden)
      @test occursin("href=\"bus_voltages_complex.csv\"", routed_hidden)
      @test occursin(".csv-grid th:nth-child(2),.csv-grid td:nth-child(2){display:none}.csv-grid th:nth-child(4),", routed_hidden)
      @test occursin("<input type=\"checkbox\" data-csv-col=\"1\"> name", routed_hidden) && occursin("<input type=\"checkbox\" data-csv-col=\"2\" checked> p_mw", routed_hidden)
      @test occursin("href=\"note_sort.csv?sort=p_mw&amp;hide=1%2C3\"", String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow/artifact/$(run_id)/note_sort.csv?raw=1&sort=p_mw&hide=1,3"; output_root).body))
      large_artifact_path = joinpath(result["output_dir"], "large_preview.log")
      large_artifact_text = repeat("0123456789abcdef", 5000)
      write(large_artifact_path, large_artifact_text)
      large_preview = SparlectraApp.handle_powerflow_artifact(run_id, "large_preview.log")
      @test large_preview.status == 200
      large_preview_html = String(large_preview.body)
      @test occursin("Preview truncated", large_preview_html)
      @test !occursin(large_artifact_text, large_preview_html)
      large_download = SparlectraApp.handle_powerflow_artifact_download(run_id, "large_preview.log")
      @test String(large_download.body) == large_artifact_text
      result_artifact = SparlectraApp.handle_powerflow_artifact(run_id, "result.json")
      @test result_artifact.status == 200
      result_artifact_html = String(result_artifact.body)
      @test occursin("Artifact: result.json", result_artifact_html)
      @test occursin("<main class=\"page artifact-page\">", result_artifact_html)
      @test occursin("class=\"artifact-text-page\"", result_artifact_html)
      @test occursin("class=\"button\" href=\"?download=1\"", result_artifact_html)
      @test occursin("class=\"artifact-text\"", result_artifact_html)

      download = SparlectraApp.handle_powerflow_artifact_download(run_id, "result.json")
      @test download.status == 200
      @test any(header -> header.first == "Content-Disposition", download.headers)
      @test SparlectraApp.route_sparlectra_webui("GET", "/powerflow/artifact/$(run_id)/result.json"; output_root).status == 200
      @test SparlectraApp.route_sparlectra_webui("GET", "/powerflow/artifact/$(run_id)/result.json?download=1"; output_root).status == 200
      zip_download = SparlectraApp.handle_powerflow_artifacts_zip(run_id)
      @test zip_download.status == 200
      @test ("Content-Type" => "application/zip") in zip_download.headers
      @test zip_download.body[1:4] == UInt8[0x50, 0x4b, 0x03, 0x04]
      @test occursin("result.json", String(zip_download.body))
      @test !occursin("Project.toml", String(zip_download.body))
      @test SparlectraApp.route_sparlectra_webui("GET", "/powerflow/artifact-zip/$(run_id)"; output_root).status == 200
      @test SparlectraApp.handle_powerflow_artifacts_zip("unknown-run").status == 404
      operation_log_text = read(operation_log_path, String)
      @test occursin("\"event\":\"artifact_opened\"", operation_log_text)
      @test occursin("\"event\":\"artifact_downloaded\"", operation_log_text)

      for unsafe_name in ("../result.json", "../../Project.toml", "/etc/passwd", raw"C:\Users\x\file.txt")
        unsafe = SparlectraApp.handle_powerflow_artifact(run_id, unsafe_name)
        @test unsafe.status == 400
        @test occursin("Unsafe artifact name rejected", String(unsafe.body))
      end

      lock(SparlectraApp._POWERFLOW_SERVICE_LOCK) do
        empty!(SparlectraApp._POWERFLOW_SERVICE_RUNS)
      end
      refreshed = SparlectraApp.handle_powerflow_refresh(output_root)
      @test run_id in refreshed["loaded_runs"]
      @test SparlectraApp.route_sparlectra_webui("POST", "/powerflow/refresh"; output_root).status == 303
      @test occursin("\"event\":\"history_refreshed\"", read(operation_log_path, String))
      history_response = SparlectraApp.handle_powerflow_history(output_root)
      @test history_response.status == 200
      history_html = String(history_response.body)
      @test occursin(run_id, history_html)
      @test occursin("<th>Date/Time</th>", history_html)
      @test occursin("class=\"status-badge status-success\"", history_html)
      @test occursin("action=\"/powerflow/delete/$(run_id)\"", history_html)
      @test occursin("action=\"/powerflow/delete_all\"", history_html)
      @test occursin("action=\"/powerflow/refresh\"", history_html)
      @test occursin(r"\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}", history_html)
      @test SparlectraApp.route_sparlectra_webui("GET", "/powerflow/abort/$(run_id)"; output_root).status == 404
      @test SparlectraApp.route_sparlectra_webui("POST", "/powerflow/abort/unknown-run"; output_root).status == 404
      @test SparlectraApp.route_sparlectra_webui("POST", "/powerflow/abort/..%2Foutside"; output_root).status == 404
      @test SparlectraApp.route_sparlectra_webui("GET", "/powerflow/delete/$(run_id)"; output_root).status == 404
      unsafe_delete = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/delete/..%2Foutside"; output_root)
      @test unsafe_delete.status == 400
      @test occursin("Unsafe PowerFlow run ID rejected", String(unsafe_delete.body))
      @test SparlectraApp.route_sparlectra_webui("GET", "/webui/shutdown"; output_root).status == 404
      @test SparlectraApp.route_sparlectra_webui("POST", "/webui/shutdown"; output_root).status == 503
      @test SparlectraApp.route_sparlectra_webui("POST", "/webui/heartbeat"; output_root).status == 204
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
      SparlectraApp._write_powerflow_run_entries!(output_root, disposable_entries)
      deleted_response = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/delete/$(disposable_id)"; output_root)
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
      lifecycle_event = (event; fields...) -> SparlectraApp.record_webui_operation!(output_root, event; route = "/powerflow/run", method = "POST", user_action = false, fields...)
      active = SparlectraApp.start_webui_powerflow_run(active_request; runner = slow_runner, event_callback = lifecycle_event)
      take!(entered_webui)
      active_id = active["run_id"]
      active_result_html = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow/result/$(active_id)"; output_root).body)
      active_form_html = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow"; output_root).body)
      active_history_html = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow/history"; output_root).body)
      abort_action = "action=\"/powerflow/abort/$(active_id)\""
      @test occursin("<form method=\"post\" $(abort_action)", active_result_html)
      @test occursin("class=\"runtime-card\"", active_result_html)
      @test occursin("Elapsed time", active_result_html)
      @test occursin("class=\"status-badge status-running\">running</span>", active_result_html)
      @test !occursin("http-equiv=\"refresh\"", active_result_html)
      @test occursin("data-refresh-url=\"/powerflow/result/$(active_id)?autorefresh=1\" data-refresh-seconds=\"$(SparlectraApp.WEBUI_STATUS_AUTO_REFRESH_SECONDS)\"", active_result_html)
      @test occursin("const scheduleAutoRefresh = function ()", active_result_html)
      # Terminal pages must be reached via a REAL reload, never by swapping the
      # fetched <main> in place: scripts inside a DOMParser-imported node are
      # inert per the HTML spec, so a swapped-in final page renders with dead
      # JS (reported as: case-combobox arrow only works after typing + Enter,
      # i.e. after the first forced navigation).
      @test occursin("newMain.hasAttribute('data-refresh-url')", active_result_html)
      @test occursin("window.location.reload()", active_result_html)
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
      auto_refresh_html = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow/result/$(active_id)?autorefresh=1"; output_root).body)
      @test occursin("data-refresh-url=\"/powerflow/result/$(active_id)?autorefresh=1\"", auto_refresh_html)
      @test count(==('\n'), read(operation_log_path, String)) == auto_refresh_count
      @test !occursin("\"event\":\"powerflow_aborted\",\"run_id\":\"$(active_id)\"", active_log)
      aborted_response = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/abort/$(active_id)"; output_root)
      @test aborted_response.status == 303
      aborted_html = String(SparlectraApp.handle_powerflow_result(active_id).body)
      @test occursin("aborting", aborted_html)
      @test occursin("class=\"runtime-card\"", aborted_html)
      @test occursin("data-refresh-url=\"/powerflow/result/$(active_id)?autorefresh=1\"", aborted_html)
      @test !occursin("http-equiv=\"refresh\"", aborted_html)
      @test occursin("Aborting requested. Current phase:", aborted_html)
      @test occursin("This phase may need to finish before cancellation is observed.", aborted_html)
      @test !occursin(abort_action, aborted_html)
      @test occursin("\"event\":\"powerflow_abort_requested\"", read(operation_log_path, String))
      @test occursin("\"current_phase\":\"linear_solve\"", read(operation_log_path, String))
      repeated_abort = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/abort/$(active_id)"; output_root)
      @test repeated_abort.status == 303
      @test occursin("\"event\":\"powerflow_abort_already_requested\"", read(operation_log_path, String))
      @test occursin("\"event\":\"powerflow_phase_started\"", read(operation_log_path, String))
      @test !any(
        line -> occursin("\"event\":\"powerflow_phase_started\"", line) && occursin("\"phase\":\"linear_solve\"", line),
        eachline(operation_log_path),
      )
      active_delete = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/delete/$(active_id)"; output_root)
      @test active_delete.status == 409
      @test occursin("This run is still active", String(active_delete.body))
      @test occursin("\"event\":\"run_delete_rejected\"", read(operation_log_path, String))
      put!(release_webui, nothing)
      wait(SparlectraApp._POWERFLOW_WEBUI_JOBS[active_id]["task"])
      terminal_html = String(SparlectraApp.handle_powerflow_result(active_id).body)
      @test occursin("<td>aborted</td>", terminal_html)
      @test occursin("class=\"runtime-card\"", terminal_html)
      @test occursin("Run aborted by user.", terminal_html)
      @test !occursin(abort_action, terminal_html)
      @test !occursin("data-refresh-url=\"", terminal_html)
      @test !occursin("http-equiv=\"refresh\"", terminal_html)
      terminal_log = read(operation_log_path, String)
      @test any(line -> occursin("\"event\":\"powerflow_aborted\"", line) && occursin("\"run_id\":\"$(active_id)\"", line), eachline(IOBuffer(terminal_log)))
      @test SparlectraApp.get_active_webui_powerflow_job() === nothing
      terminal_delete = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/delete/$(active_id)"; output_root)
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
        "abort_requested_at" => Dates.now(Dates.UTC) - Dates.Second(SparlectraApp.WEBUI_ABORT_HARD_RESET_AFTER_SECONDS + 1),
        "current_phase" => "linear_solve",
        "phase_started_at" => Dates.now(Dates.UTC),
        "last_progress_at" => Dates.now(Dates.UTC),
      )
      SparlectraApp._POWERFLOW_WEBUI_JOBS["hard-reset-test"] = hard_reset_job
      hard_reset_html = String(SparlectraApp.handle_powerflow_result("hard-reset-test").body)
      @test occursin("Hard reset Web UI", hard_reset_html)
      hard_reset_response = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/hard-reset/hard-reset-test"; output_root)
      @test hard_reset_response.status == 200
      @test SparlectraApp.get_webui_powerflow_job("hard-reset-test")["status"] == "aborted_unknown"
      @test occursin("\"event\":\"webui_hard_reset_requested\"", read(operation_log_path, String))
      @test !SparlectraApp.get_webui_powerflow_job("hard-reset-test")["success"]

      operation_log_page = SparlectraApp.route_sparlectra_webui("GET", "/webui/operation-log"; output_root)
      @test operation_log_page.status == 200
      @test occursin("Download operation log", String(operation_log_page.body))
      operation_log_download = SparlectraApp.route_sparlectra_webui("GET", "/webui/operation-log/download"; output_root)
      @test operation_log_download.status == 200
      @test any(header -> header.first == "Content-Disposition", operation_log_download.headers)
      static_count_before = count(==('\n'), read(operation_log_path, String))
      @test SparlectraApp.route_sparlectra_webui("GET", "/static/sparlectra.css"; output_root).status == 200
      @test count(==('\n'), read(operation_log_path, String)) == static_count_before
      compact_root = joinpath(tmpdir, "compact-log")
      compact_log = SparlectraApp.webui_operation_log_path(compact_root)
      mkpath(compact_root)
      open(compact_log, "w") do io
        println(io)
        println(io, "malformed")
        for index in 1:11
          println(io, "{\"event\":\"entry-$(index)\"}")
        end
      end
      @test SparlectraApp._compact_webui_operation_log!(compact_log; max_entries = 10, keep_entries = 3)
      compact_lines = readlines(compact_log)
      @test length(compact_lines) == 3
      @test all(line -> SparlectraApp._parse_service_json(line) isa AbstractDict, compact_lines)
      @test occursin("\"event\":\"entry-9\"", compact_lines[1])
      @test occursin("\"event\":\"entry-11\"", compact_lines[3])
      @test !SparlectraApp._compact_webui_operation_log!(compact_log; max_entries = 10, keep_entries = 3)
      @test !SparlectraApp._compact_webui_operation_log!(joinpath(tmpdir, "missing-operation-log.jsonl"); max_entries = 10, keep_entries = 3)
      compact_viewer = SparlectraApp.handle_webui_operation_log(compact_log)
      compact_download = SparlectraApp.handle_webui_operation_log(compact_log; download = true)
      @test compact_viewer.status == 200
      @test compact_download.status == 200
      @test occursin("entry-11", String(compact_viewer.body))
      @test occursin("entry-11", String(compact_download.body))
      logging_failure_root = joinpath(tmpdir, "logging-is-a-file")
      write(logging_failure_root, "not a directory")
      record_result = @test_logs (:warn, r"Could not record Web UI operation") begin
        SparlectraApp.record_webui_operation!(logging_failure_root, "expected_failure")
      end
      @test record_result == false
      route_result = @test_logs (:warn, r"Could not record Web UI operation") begin
        SparlectraApp.route_sparlectra_webui("GET", "/powerflow"; output_root = logging_failure_root)
      end
      @test route_result.status == 200
      empty_delete_root = joinpath(tmpdir, "empty-delete-root")
      @test SparlectraApp.route_sparlectra_webui("POST", "/powerflow/delete_all"; output_root = empty_delete_root).status == 303
      @test occursin("\"event\":\"all_runs_deleted\"", read(SparlectraApp.webui_operation_log_path(empty_delete_root), String))

      @test occursin(run_id, String(SparlectraApp.handle_powerflow_result(run_id).body))

      lock(SparlectraApp._POWERFLOW_SERVICE_LOCK) do
        empty!(SparlectraApp._POWERFLOW_SERVICE_RUNS)
      end
      probe = listen(ip"127.0.0.1", UInt16(0))
      port = Int(getsockname(probe)[2])
      close(probe)
      lifecycle_io = IOBuffer()
      server = redirect_stderr(devnull) do
        start_sparlectra_webui(port = port, output_root = output_root, auto_shutdown_on_browser_close = false, _lifecycle_io = lifecycle_io)
      end
      startup_output = String(take!(lifecycle_io))
      _assert_no_webui_internal_lifecycle_markers(startup_output)
      @test occursin("Sparlectra Web UI is available at http://127.0.0.1:$(port)/powerflow", startup_output)
      @test !istaskdone(server.task)
      @test isdir(output_root)
      @test server.runtime.case_directory == SparlectraApp.default_webui_case_cache_dir(output_root)
      @test isdir(server.runtime.case_directory)
      @test isfile(server.runtime.config_file)
      @test isfile(joinpath(server.runtime.case_directory, "warmup_case3.jl"))
      @test isfile(joinpath(server.runtime.case_directory, "warmup_case118.jl"))
      @test isfile(server.runtime.operation_log)
      # path boundary, not a string prefix: a sibling directory such as
      # <pkg>-tmp must not count as inside the package
      @test !startswith(server.runtime.case_directory, joinpath(normpath(pkgdir(Sparlectra)), ""))
      @test !startswith(server.runtime.config_file, joinpath(normpath(pkgdir(Sparlectra)), ""))
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
      @test occursin(basename(SparlectraApp._sparlectra_package_path()), startup_log)
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
        @test occursin("Case cache", response_text)
        @test occursin(server.runtime.case_directory, response_text)
        @test occursin("Operation log", response_text)
        @test occursin(server.runtime.operation_log, response_text)
        # the import conventions render on the Case page; smoke
        # the new route over the REAL server socket
        case_response_text = _webui_http_request(port, "GET", "/powerflow/case")
        @test occursin("HTTP/1.1 200 OK", case_response_text)
        @test occursin("MATPOWER import conventions", case_response_text)
        @test occursin("name=\"matpower_import_auto_profile\"", case_response_text)
        @test occursin("name=\"matpower_import_ratio\"", case_response_text)
        @test occursin("name=\"matpower_import_shift_sign\"", case_response_text)
        @test occursin("name=\"matpower_import_shift_unit\"", case_response_text)

        logo_response_text = _webui_http_request(port, "GET", "/assets/logo.png")
        @test occursin("HTTP/1.1 200 OK", logo_response_text)
        @test occursin("Content-Type: image/png", logo_response_text)

        heartbeat_response = _webui_http_request(port, "POST", "/webui/heartbeat")
        @test occursin("HTTP/1.1 204", heartbeat_response)
        # Stop Web UI sends the window geometry; the server saves it into
        # its configuration file before it stops
        shutdown_response = _webui_http_request(port, "POST", "/webui/shutdown"; body = "window_width=1234&window_height=876&window_x=12&window_y=34&window_hidden=0")
        @test timedwait(() -> istaskdone(server.task), 2.0) == :ok
        @test occursin("HTTP/1.1 200 OK", shutdown_response)
        @test !isopen(server.listener)
        stopped_webui = Sparlectra.load_sparlectra_config(server.runtime.config_file; reload = true).webui
        @test (stopped_webui.window_width, stopped_webui.window_height, stopped_webui.window_x, stopped_webui.window_y) == (1234, 876, 12, 34)
        @test occursin("\"event\":\"webui_window_geometry\"", read(server.runtime.operation_log, String))
      finally
        close(server)
      end

      for (name, config_text, expected_marker) in (
        ("new qlimit enforcement mode", "power_flow:\n  qlimits:\n    enforcement_mode: classic_simultaneous\n", "PowerFlow run"),
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
      @test SparlectraApp._wait_sparlectra_webui(interrupted; wait_for_task = _ -> throw(InterruptException()))
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
      SparlectraApp._webui_begin_request!(heartbeat_server.runtime)
      sleep(0.3)
      @test isopen(heartbeat_server.listener)
      SparlectraApp._webui_finish_request!(heartbeat_server.runtime)
      sleep(0.3)
      @test isopen(heartbeat_server.listener)
      close(heartbeat_server)
      @test timedwait(() -> istaskdone(heartbeat_server.task), 2.0) == :ok
      @test !isopen(heartbeat_server.listener)
    end

    @testset "SE settings: one source, and the value arrives" begin (function ()
      # Twice in two days a setting existed, was displayed, and did not act:
      # rankTolFactor was taken and not forwarded, and literal defaults in
      # the service signature outranked the configuration (a 25000-bus run
      # estimated with k_suppress 6.0 while the configuration said 4.0).
      # Both falsified measurements before anyone noticed. This test walks
      # the SE option list mechanically instead of naming single keys, so a
      # third occurrence fails here rather than in someone's results.
      se_cfg = Sparlectra.state_estimation_config()
      # form default == configuration default, for every option that mirrors
      # a configuration field: by its config_key path where it has one
      # (nested keys such as tap_estimation.impedance_correction), else by
      # name (se_<field>). A field whose name meets a configuration BLOCK
      # (se_tap_estimation and the tap_estimation section) mirrors nothing.
      mirrored = 0
      for spec in SparlectraApp.WEBUI_OPTION_SPECS
        startswith(spec.field, "se_") || continue
        path = spec.config_key !== nothing && startswith(spec.config_key, "state_estimation.") ?
               Symbol.(split(spec.config_key, ".")[2:end]) : [Symbol(spec.field[4:end])]
        want = se_cfg
        for p in path
          want = hasproperty(want, p) ? getproperty(want, p) : nothing
          want === nothing && break
        end
        want === nothing && continue
        want isa Union{Symbol,Number,AbstractString} || continue
        mirrored += 1
        got = spec.default
        if want isa Symbol
          @test Symbol(String(got)) === want
        elseif want isa AbstractFloat
          @test Float64(got) ≈ Float64(want)
        else
          @test got == want
        end
      end
      # the list must actually contain the interesting ones, otherwise this
      # loop could pass while covering nothing
      @test mirrored >= 7
      # the tolerance sits AT the finite-difference floor, never below it
      @test se_cfg.tol >= se_cfg.jac_eps

      # the bad-data fields of the inactive mode are grayed out, not dropped:
      # a value saved for the other mode still renders back
      se_args = (; cases = ["sp_case14.scf.json"], selected_case = "sp_case14.scf.json",
                 measurements = ["sp_case14.measurements.csv"])
      staged = SparlectraApp.render_se_form(; se_args...,
        gen_values = Dict{String,String}("se_robust_mode" => "staged", "se_k_suppress" => "7.5", "se_robust_k1" => "2.5"))
      @test occursin("name=\"se_k_suppress\" id=\"se-k-suppress\" value=\"7.5\"", staged)
      @test occursin("name=\"se_robust_k1\" value=\"2.5\"", staged)
      replacement = SparlectraApp.render_se_form(; se_args...,
        gen_values = Dict{String,String}("se_robust_mode" => "replacement", "se_robust_k2" => "9.0"))
      @test occursin("name=\"se_robust_k2\" value=\"9.0\"", replacement)

      # and the other direction: a value that differs from the default has
      # to reach the run. The estimator reports what it used in its metadata,
      # so the run itself is the witness, not an internal.
      mktempdir() do d
        case = joinpath(Sparlectra.SPARLECTRA_ROOT, "data", "scf", "sp_case14.scf.json")
        res = redirect_stdout(devnull) do
          SparlectraApp._run_state_estimation_service(case, Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH,
            joinpath(d, "run"), "se_settings_arrive", joinpath(d, "no_set.csv");
            max_iter = 7, tol = 1.0e-4, k_eliminate = 4.25, robust_mode = :replacement,
            k_suppress = 5.5, suppression_sigma = 1234.0,
            config_overrides = Dict{String,Any}("state_estimation.tap_estimation.impedance_correction" => "off"))
        end
        md = SparlectraApp.to_dict(res)["metadata"]
        @test md["se_tap_impedance_correction"] == "off"
        @test md["se_max_iter"] == 7
        @test Float64(md["se_tol"]) ≈ 1.0e-4
        @test Float64(md["se_k_eliminate"]) ≈ 4.25
        @test String(md["se_robust_mode"]) == "replacement"
        @test Float64(md["se_k_suppress"]) ≈ 5.5
        @test Float64(md["se_suppression_sigma"]) ≈ 1234.0
      end
    end)() end

  end)() end

  @testset "scenario editor and sources" begin (function ()
    # suite migration (no-download rule): the flows run on the shipped
    # sp_case60 bundle instead of a locally cached case118, so the set
    # ALWAYS runs (the old availability gate silently skipped fresh
    # installations)
    begin
      mktempdir() do dir
        root = joinpath(dir, "out")
        mkpath(root)
        casedir = joinpath(dir, "cases")
        mkpath(casedir)
        config_path = Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH
        rt = SparlectraApp._SparlectraWebUIRuntime(nothing, casedir, config_path, SparlectraApp.webui_operation_log_path(root), nothing, SparlectraApp.start_powerflow_run, false, false, time(), 0, nothing, IOBuffer(), ReentrantLock())

        # the shipped SCF bundle in the case directory (it already carries
        # a scenarios block; the flows add their own scenario next to it)
        scf_name = "sp_case60.scf.json"
        scf_path = joinpath(casedir, scf_name)
        cp(joinpath(dirname(@__DIR__), "data", "scf", scf_name), scf_path)
        # the shipped block is mode explicit (the named
        # scenarios are what an editor shows); no staging rewrite needed
        scfcase = Sparlectra.read_scf_json(scf_path)
        idx = Sparlectra.ScenarioIndex(scfcase)
        branch_id = first(r.id for r in scfcase.data.line)
        gen_id = first(id for (id, k) in idx.kind_by_id if k === :generator)
        load_id = first(id for (id, k) in idx.kind_by_id if k === :load)

        # flow 1: create a three-op scenario through the form, save, reload,
        # run with :flag, see the row
        pre_save_bytes = read(scf_path, String)
        form = Dict{String,Any}(
          "casefile" => scf_name,
          "original_name" => "",
          "scenario_name" => "step6 triple",
          "scenario_weight" => "1.5",
          "op_op" => ["status", "set", "scale"],
          "op_target" => ["branch", "generator", "load"],
          "op_component" => [string(branch_id), string(gen_id), string(load_id)],
          "op_field" => ["", "p", ""],
          "op_value" => ["0", "25.0", ""],
          "op_factor" => ["", "", "1.2"],
        )
        SparlectraApp.route_sparlectra_webui("POST", "/powerflow/scenarios/save", form; output_root = root, runtime = rt)
        page = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow/scenarios?case=$(SparlectraApp._webui_urlencode(scf_name))"; output_root = root, runtime = rt).body)
        @test occursin("step6 triple", page)
        @test occursin("<td>3</td>", page)
        run_result = SparlectraApp.start_powerflow_run(Dict{String,Any}("casefile" => scf_path, "config_file" => config_path, "output_root" => root, "contingency_mode" => true, "contingency_kind" => "branch", "scenario_source" => "file_block", "screening_mode" => "flag"))
        @test run_result["status"] == "succeeded"
        @test run_result["metadata"]["contingency_cases_source"] == "scenario_file_block"
        @test run_result["metadata"]["contingency_screening_mode"] == "flag"
        csv_text = read(joinpath(root, run_result["run_id"], "contingency_n1.csv"), String)
        @test occursin("step6 triple", csv_text)
        # the last two columns mark the screened cases
        @test run_csv(joinpath(root, run_result["run_id"], "contingency_n1.csv")).header[(end-1):end] == ["screened", "screening_estimate"]

        # the result page renders the N-1 table (the screened marker must be visible in the browser, not only in the
        # CSV). The file_block run carries ONE multi-op scenario, which is
        # correctly never screened (full solve, no estimate line), so the
        # screened-row assertions run on an n1_all :flag batch below.
        flag_page = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow/result/$(SparlectraApp._webui_urlencode(run_result["run_id"]))"; output_root = root, runtime = rt).body)
        # a run of the file block is titled as scenarios, not as N-1
        @test occursin("Scenarios result", flag_page) && !occursin("N-1 contingency result", flag_page)
        @test occursin("<th>screened</th>", flag_page)
        @test occursin("step6 triple", flag_page)
        n1_flag_run = SparlectraApp.start_powerflow_run(Dict{String,Any}("casefile" => scf_path, "config_file" => config_path, "output_root" => root, "contingency_mode" => true, "contingency_kind" => "branch", "scenario_source" => "n1_all", "screening_mode" => "flag"))
        @test n1_flag_run["status"] == "succeeded"
        @test n1_flag_run["metadata"]["contingency_screened"] > 0
        n1_flag_page = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow/result/$(SparlectraApp._webui_urlencode(n1_flag_run["run_id"]))"; output_root = root, runtime = rt).body)
        @test occursin("one-step estimate (no full solve)", n1_flag_page)
        @test occursin("flagged: estimate before the full run", n1_flag_page)
        @test occursin("class=\"n1-screened\">yes<", n1_flag_page)
        # an :off run renders the table WITHOUT the screened column
        off_run = SparlectraApp.start_powerflow_run(Dict{String,Any}("casefile" => scf_path, "config_file" => config_path, "output_root" => root, "contingency_mode" => true, "contingency_kind" => "branch", "scenario_source" => "file_block", "screening_mode" => "off"))
        @test off_run["status"] == "succeeded"
        off_page = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow/result/$(SparlectraApp._webui_urlencode(off_run["run_id"]))"; output_root = root, runtime = rt).body)
        @test occursin("Scenarios result", off_page)
        @test !occursin("<th>screened</th>", off_page)

        # the editor save touches only the scenarios block: with that block
        # removed from both documents, the re-serialized files are byte
        # identical (stronger than a field comparison; the writer rewrote
        # the whole file and must not have changed anything else)
        doc_before = Sparlectra.scf_json_parse(pre_save_bytes)
        doc_after = Sparlectra.scf_json_parse(read(scf_path, String))
        haskey(doc_before, "sparlectra") && delete!(doc_before["sparlectra"], "scenarios")
        haskey(doc_after, "sparlectra") && delete!(doc_after["sparlectra"], "scenarios")
        @test Sparlectra.scf_json_string(doc_before) == Sparlectra.scf_json_string(doc_after)

        # lost-update guard: a form carrying the digest of an OLDER file
        # state is rejected with the way out, and the file stays untouched
        stale_form = Dict{String,Any}(
          "casefile" => scf_name,
          "original_name" => "",
          "scenario_name" => "stale write",
          "scenario_weight" => "1.0",
          "op_op" => ["status"],
          "op_target" => ["branch"],
          "op_component" => [string(branch_id)],
          "op_field" => [""],
          "op_value" => ["0"],
          "op_factor" => [""],
          "file_digest" => "deadbeef",
        )
        before_bytes = read(scf_path, String)
        stale_body = String(SparlectraApp.route_sparlectra_webui("POST", "/powerflow/scenarios/save", stale_form; output_root = root, runtime = rt).body)
        @test occursin("changed since this form was loaded", stale_body)
        @test read(scf_path, String) == before_bytes
        # with the CURRENT digest the same save goes through
        stale_form["file_digest"] = Sparlectra._config_file_hash(scf_path)
        SparlectraApp.route_sparlectra_webui("POST", "/powerflow/scenarios/save", stale_form; output_root = root, runtime = rt)
        @test occursin("stale write", read(scf_path, String))

        # flow 2: a tap_pos op on a REGULATED transformer is rejected with
        # the controller named (net-aware step-2 validation through the form)
        regnet = Net(name = "step6_reg", baseMVA = 100.0)
        addBus!(net = regnet, busName = "Slack", vn_kV = 110.0)
        addBus!(net = regnet, busName = "Mid", vn_kV = 110.0)
        addBus!(net = regnet, busName = "Load", vn_kV = 110.0)
        addProsumer!(net = regnet, busName = "Slack", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.02, va_deg = 0.0, referencePri = "Slack")
        addProsumer!(net = regnet, busName = "Load", type = "ENERGYCONSUMER", p = -70.0, q = -20.0)
        addPIModelTrafo!(net = regnet, fromBus = "Slack", toBus = "Mid", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, ratio = 1.0, shift_deg = 0.0, status = 1)
        addPIModelACLine!(net = regnet, fromBus = "Mid", toBus = "Load", r_pu = 0.02, x_pu = 0.12, b_pu = 0.01, status = 1)
        tbr = getNetBranch(net = regnet, fromBus = "Slack", toBus = "Mid")
        tbr.has_ratio_tap = true
        tbr.tap_min = 0.95
        tbr.tap_max = 1.05
        tbr.tap_step = 0.0125
        addPowerTransformerControl!(regnet; trafo = string(tbr.branchIdx), mode = :voltage, target_bus = "Load", target_vm_pu = 1.0, control_ratio = true, control_phase = false)
        reg_name = "step6_regulated.scf.json"
        Sparlectra.exportSCF(regnet; file = joinpath(casedir, reg_name), case_name = "step6_regulated")
        regcase = Sparlectra.read_scf_json(joinpath(casedir, reg_name))
        regidx = Sparlectra.ScenarioIndex(regcase)
        trafo_id = first(id for (id, k) in regidx.kind_by_id if k === :transformer)
        reject_form = Dict{String,Any}(
          "casefile" => reg_name,
          "original_name" => "",
          "scenario_name" => "tap fight",
          "scenario_weight" => "1.0",
          "op_op" => ["set"],
          "op_target" => ["transformer"],
          "op_component" => [string(trafo_id)],
          "op_field" => ["tap_pos"],
          "op_value" => ["1"],
          "op_factor" => [""],
        )
        reject_body = String(SparlectraApp.route_sparlectra_webui("POST", "/powerflow/scenarios/save", reject_form; output_root = root, runtime = rt).body)
        @test occursin("regulated by controller", reject_body)
        @test occursin("tap fight", reject_body)
        # the rejected scenario is NOT in the file
        @test !occursin("tap fight", read(joinpath(casedir, reg_name), String))

        # a failed case shows its error text in the result-table row: the
        # generator N-1 on the regulated net takes out its only slack unit
        fail_run = SparlectraApp.start_powerflow_run(Dict{String,Any}("casefile" => joinpath(casedir, reg_name), "config_file" => config_path, "output_root" => root, "contingency_mode" => true, "contingency_kind" => "gen", "screening_mode" => "off"))
        @test fail_run["status"] == "succeeded"
        fail_csv = read(joinpath(root, fail_run["run_id"], "contingency_n1.csv"), String)
        @test occursin("false", fail_csv)
        fail_page = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow/result/$(SparlectraApp._webui_urlencode(fail_run["run_id"]))"; output_root = root, runtime = rt).body)
        @test occursin("N-1 contingency result", fail_page)
        @test occursin("n1-error", fail_page)

        # flow 3: a MATPOWER case gets the export hint and no editor, on the
        # scenarios page and on the main form
        write(joinpath(casedir, "step6_plain.m"), "% case fixture\n")
        plain_page = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow/scenarios?case=step6_plain.m"; output_root = root, runtime = rt).body)
        @test occursin("Scenarios need an SCF case", plain_page)
        @test !occursin("scenario-op-rows", plain_page)
        main_plain = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow?casefile=$(SparlectraApp._webui_urlencode("step6_plain.m"))"; output_root = root, runtime = rt).body)
        @test occursin("scenario-scf-hint", main_plain)
        @test !occursin("data-scenarios-button", main_plain)
        # 0.40.0: an SCF case WITH a scenarios block gets the Scenarios
        # button, which runs the block alone; the N-1 list no longer offers
        # the block, so there is one way to the scenarios
        main_scf = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow?casefile=$(SparlectraApp._webui_urlencode(scf_name))"; output_root = root, runtime = rt).body)
        @test occursin("data-scenarios-button", main_scf)
        @test !occursin("value=\"file_block\"", main_scf)
        @test occursin("scenario-editor-link", main_scf)
        # an SCF case without a block: the editor, but no button
        bare = Sparlectra.read_scf_json(scf_path)
        bare.sparlectra.scenarios = Dict{String,Any}()
        Sparlectra.write_scf_json(bare, joinpath(casedir, "step6_bare.scf.json"))
        main_bare = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow?casefile=step6_bare.scf.json"; output_root = root, runtime = rt).body)
        @test occursin("scenario-editor-link", main_bare)
        @test !occursin("data-scenarios-button", main_bare)
        # the button's request: a contingency run on the file block alone,
        # whatever the N-1 list says
        scen_request = SparlectraApp.powerflow_webui_request(Dict{String,Any}("casefile" => scf_name, "scenario_run" => "true", "scenario_source" => "n1_all"); default_output_root = root)
        @test scen_request["contingency_mode"] === true
        @test scen_request["scenario_source"] == "file_block"
        n1_request = SparlectraApp.powerflow_webui_request(Dict{String,Any}("casefile" => scf_name, "contingency_mode" => "true", "scenario_source" => "n1_all"); default_output_root = root)
        @test n1_request["scenario_source"] == "n1_all"

        # flow 4: the outage counter of a RUNNING N-1 job. The
        # runner holds the batch inside the progress callback of the first
        # finished outage, so the job is observably mid-run: the snapshot
        # carries progress_done/progress_total, and the status page the
        # browser polls (main[data-refresh-url], swapped every
        # WEBUI_STATUS_AUTO_REFRESH_SECONDS) renders the counter element
        # from it. After release the counter ends at done == total.
        ct_entered = Channel{Nothing}(1)
        ct_release = Channel{Symbol}(1)
        ct_runner = function (request; case_directory = nothing)
          inner = request["progress_callback"]
          request["progress_callback"] = function (done, total)
            inner(done, total)
            # exactly one call has done == 1, also on a threaded batch
            if done == 1
              put!(ct_entered, nothing)
              take!(ct_release)
            end
            return nothing
          end
          return SparlectraApp.start_powerflow_run(request; case_directory)
        end
        op_log_path = SparlectraApp.webui_operation_log_path(root)
        ct_event = (event; fields...) -> SparlectraApp.record_webui_operation!(root, event; route = "/powerflow/run", method = "POST", user_action = false, fields...)
        # detailed_result_csv on as the form sends it by default: an N-1
        # batch writes no bus or branch CSV, so the job logs the export
        # skipped, not failed (checked after the released run below)
        ct_request = Dict{String,Any}("casefile" => scf_path, "config_file" => config_path, "output_root" => root, "contingency_mode" => true, "contingency_kind" => "branch", "scenario_source" => "n1_branches", "screening_mode" => "off", "detailed_result_csv" => true)
        n_outages = length(Sparlectra.generateN1Branches(Sparlectra.importSCF(scf_path)))
        for finish in (:release, :abort)
          ct_job = SparlectraApp.start_webui_powerflow_run(copy(ct_request); runner = ct_runner, event_callback = ct_event)
          ct_id = ct_job["run_id"]
          take!(ct_entered)
          mid = SparlectraApp.get_webui_powerflow_job(ct_id)
          @test mid["status"] == "running"
          # L19: the job names the batch phase while the counter runs, not
          # the last phase before the service call (preparing_configuration)
          @test mid["current_phase"] == "running_outages"
          @test mid["progress_total"] == n_outages
          @test 1 <= mid["progress_done"] < n_outages
          mid_html = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow/result/$(ct_id)?autorefresh=1"; output_root = root, runtime = rt).body)
          # the counter element with the snapshot numbers, inside the
          # polled <main>, and the poll script that re-fetches it
          # (on a threaded batch the other chunks keep counting between the
          # snapshot and the render, so the page may show a later count, but
          # never an earlier one and never the end: one chunk is held)
          counter = match(r"id=\"batch-progress\" data-progress-done=\"(\d+)\" data-progress-total=\"(\d+)\"", mid_html)
          @test counter !== nothing
          if counter !== nothing
            @test parse(Int, counter[2]) == n_outages
            @test mid["progress_done"] <= parse(Int, counter[1]) < n_outages
          end
          @test occursin("data-refresh-url=\"/powerflow/result/$(ct_id)?autorefresh=1\"", mid_html)
          @test occursin("const scheduleAutoRefresh = function ()", mid_html)
          if finish === :abort
            # the abort token still reaches the batch: the next outage
            # check after the release ends the run as aborted
            @test SparlectraApp.route_sparlectra_webui("POST", "/powerflow/abort/$(ct_id)"; output_root = root, runtime = rt).status == 303
          end
          put!(ct_release, finish)
          wait(SparlectraApp._POWERFLOW_WEBUI_JOBS[ct_id]["task"])
          final_job = SparlectraApp._POWERFLOW_WEBUI_JOBS[ct_id]
          if finish === :release
            @test final_job["status"] == "success"
            @test final_job["progress_done"] == final_job["progress_total"] == n_outages
            ct_csv_events = [l for l in eachline(op_log_path) if occursin(ct_id, l) && occursin("\"event\":\"detailed_result_csv_export", l)]
            @test any(l -> occursin("\"event\":\"detailed_result_csv_export_skipped\"", l) && occursin("\"reason\":\"contingency_run\"", l), ct_csv_events)
            @test !any(l -> occursin("\"event\":\"detailed_result_csv_export_failed\"", l), ct_csv_events)
            # L19: the phases of an N-1 run reach the operation log in run order
            ct_phases = [m[1] for m in (match(r"\"phase\":\"([a-z_]+)\"", l) for l in eachline(op_log_path) if occursin(ct_id, l) && occursin("\"event\":\"powerflow_phase_started\"", l)) if m !== nothing]
            ct_order = [findfirst(==(p), ct_phases) for p in ("importing_case", "solving_base_case", "running_outages", "writing_artifacts")]
            @test all(!isnothing, ct_order) && issorted(ct_order)
            # a finished job renders no counter (the N-1 table takes over)
            @test !occursin("id=\"batch-progress\"", String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow/result/$(ct_id)"; output_root = root, runtime = rt).body))
          else
            @test final_job["status"] == "aborted"
            @test final_job["progress_done"] < n_outages
          end
          # the operation log stays high level: no line per outage
          @test count(line -> occursin(ct_id, line), eachline(op_log_path)) < n_outages
        end
      end
    end
  end)() end
  _webui_qsts_tests()
  _webui_short_circuit_status_tests()
  return nothing
end

# T13 (QSTS in the Web UI, 0.60.0): the QSTS button of the Runs page on the
# shipped example series sp_case14_day.qsts.zip. The run goes through the
# form route as the browser posts it (casefile, config_file, the select
# qsts_mode_choice and the submitter qsts_mode), the bundled series is
# staged with its companion sp_case14_day.control.yaml, and the result
# carries run_mode "qsts", the files, one SVG per element with a point per
# solved step, and the chart route draws the same bytes. A case without a
# series gets a disabled button, and a request with qsts_mode on it fails
# with qsts_no_series.
function _webui_qsts_tests()
  @testset "QSTS run from the run form (T13)" begin (function ()
    mktempdir() do dir
      root = joinpath(dir, "out")
      casedir = joinpath(dir, "cases")
      mkpath(root)
      mkpath(casedir)
      config_path = Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH
      rt = SparlectraApp._SparlectraWebUIRuntime(nothing, casedir, config_path, SparlectraApp.webui_operation_log_path(root), nothing, SparlectraApp.start_powerflow_run, false, false, time(), 0, nothing, IOBuffer(), ReentrantLock())
      series = "sp_case14_day.qsts.zip"
      # the case chooser offers the shipped series (by content)
      case_page = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow/case"; output_root = root, runtime = rt).body)
      @test occursin("data-case-option=\"$(series)\"", case_page)
      # the Runs page of the series: the button posts qsts_mode, enabled,
      # inside the run form, with the mode select next to it
      page = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow?casefile=$(series)"; output_root = root, runtime = rt).body)
      form_html = page[findfirst("<form id=\"powerflow-run-form\"", page)[1]:findnext("</form>", page, findfirst("<form id=\"powerflow-run-form\"", page)[1])[end]]
      button = match(r"<button[^>]*data-qsts-button[^>]*>", form_html)
      @test button !== nothing
      @test occursin("name=\"qsts_mode\"", button.match)
      @test occursin("data-qsts-state=\"series\"", button.match)
      @test !occursin(" disabled", button.match)
      @test occursin("<select name=\"qsts_mode_choice\"", form_html)
      # a case without a series: the button is disabled, so a browser never
      # sends qsts_mode for it (a disabled control is not in the form data)
      cp(joinpath(dirname(@__DIR__), "data", "scf", "sp_case14.scf.json"), joinpath(casedir, "sp_case14.scf.json"))
      plain_page = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow?casefile=sp_case14.scf.json"; output_root = root, runtime = rt).body)
      plain_button = match(r"<button[^>]*data-qsts-button[^>]*>", plain_page)
      @test plain_button !== nothing && occursin("data-qsts-state=\"no-series\"", plain_button.match) && occursin(" disabled", plain_button.match)
      # the mode select is format-bound: an SCF file reads no QSTS option,
      # so the page leaves it out, and a series leaves out the options of
      # the other formats
      @test !occursin("name=\"qsts_mode_choice\"", plain_page)
      @test !occursin("name=\"matpower_import_ratio\"", String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow/case?casefile=$(series)"; output_root = root, runtime = rt).body))

      # the POST of the form: the select names the mode, the button marks the run
      form = Dict{String,Any}("casefile" => series, "config_file" => config_path, "qsts_mode_choice" => "sequential", "qsts_mode" => "sequential")
      response = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/run", form; output_root = root, runtime = rt)
      @test response.status == 303
      run_id = basename(only(h.second for h in response.headers if h.first == "Location"))
      wait(SparlectraApp._POWERFLOW_WEBUI_JOBS[run_id]["task"])
      result = get_powerflow_result(run_id)
      @test result["success"] === true
      md = result["metadata"]
      @test md["run_mode"] == "qsts"
      @test md["qsts_mode"] == "sequential"
      # staged with its companion, and run.log names it
      @test isfile(joinpath(casedir, "sp_case14_day.control.yaml"))
      od = result["output_dir"]
      @test occursin("control section from sp_case14_day.control.yaml", read(joinpath(od, "run.log"), String))
      # the series states 50 Hz but carries c1 = 0 on every line: the
      # frequency never entered the model, and run.log says so
      @test md["frequency_source"] == "not_used"
      @test occursin("System frequency: not used (no capacitance in farad); case file states 50 Hz\n", read(joinpath(od, "run.log"), String))
      for f in ("qsts_values.csv", "qsts_steps.csv", "qsts_elements.csv", "qsts_switching.csv", "qsts_metrics.csv", "qsts_manifest.json", "qsts_view.json", "qsts_overview.html")
        @test isfile(joinpath(od, f))
      end
      # one SVG per element; the first controller's chart has a point per solved step
      n_solved = count(row -> "solved" in row, run_csv(joinpath(od, "qsts_steps.csv")).rows)
      @test n_solved == md["qsts_solved"] > 0
      default_chart = md["qsts_default_charts"][1]
      view = Sparlectra.scf_json_parse(read(joinpath(od, "qsts_view.json"), String))
      ctl = only(el for el in view["elements"] if el["id"] == default_chart[1])
      @test ctl["kind"] == "controller" && default_chart[2] == "position"
      svg_text = read(joinpath(od, ctl["plot"]), String)
      @test count("<circle", svg_text) == n_solved
      @test count(el -> haskey(el, "plot"), view["elements"]) == md["qsts_plots"] == length(readdir(joinpath(od, "qsts_plots")))
      # the chart route draws the bytes of the service's SVG
      chart = SparlectraApp.route_sparlectra_webui("GET", "/powerflow/result/$(run_id)/qsts-chart?element=$(SparlectraApp._webui_urlencode(ctl["id"]))&quantity=position"; output_root = root, runtime = rt)
      @test chart.status == 200
      @test any(h -> h.first == "Content-Type" && startswith(h.second, "image/svg+xml"), chart.headers)
      @test String(chart.body) == svg_text
      @test SparlectraApp.route_sparlectra_webui("GET", "/powerflow/result/$(run_id)/qsts-chart?element=nobody&quantity=vm_pu"; output_root = root, runtime = rt).status == 400
      # the result page: the summary row and the inline chart
      result_page = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow/result/$(run_id)"; output_root = root, runtime = rt).body)
      @test occursin("QSTS result", result_page)
      @test occursin("data-qsts-chart", result_page) && occursin("<svg", result_page)
      # the overview link of the result page opens the overview itself as
      # text/html, not its source in a text panel, and every chart it embeds
      # resolves relative to it, the way the browser resolves <object data>
      overview_href = match(r"href=\"([^\"]*/qsts_overview\.html)\"", result_page)
      @test overview_href !== nothing
      overview = SparlectraApp.route_sparlectra_webui("GET", overview_href[1]; output_root = root, runtime = rt)
      @test overview.status == 200
      @test any(h -> h.first == "Content-Type" && startswith(h.second, "text/html"), overview.headers)
      # the page policy blocks script, but has to let the charts load:
      # Chromium checks an <object> document against frame-src
      policy = only(h.second for h in overview.headers if h.first == "Content-Security-Policy")
      directives = Dict(String(p[1]) => String.(p[2:end]) for p in (split(strip(d)) for d in split(policy, ';')) if !isempty(p))
      @test "'self'" in get(directives, "frame-src", String[]) && "'self'" in get(directives, "object-src", String[])
      @test !haskey(directives, "script-src") && get(directives, "default-src", String[]) == ["'none'"]
      # String(bytes) takes the vector over, so the body is read once
      overview_text = String(overview.body)
      # the response adds a link bar after <body> (the app window has no
      # back button); without it the document is the file on disk
      @test occursin("href=\"/powerflow/result/$(run_id)\"", overview_text)
      @test occursin("href=\"/powerflow/artifacts/$(run_id)\"", overview_text)
      @test replace(overview_text, r"<nav class=\"artifact-back\".*?</nav>" => "") == read(joinpath(od, "qsts_overview.html"), String)
      plot_refs = [m[1] for m in eachmatch(r"<object type=\"image/svg\+xml\" data=\"([^\"]+)\"", overview_text)]
      @test length(plot_refs) == md["qsts_plots"]
      overview_base = overview_href[1][1:findlast('/', overview_href[1])]
      plot_responses = [SparlectraApp.route_sparlectra_webui("GET", overview_base * ref; output_root = root, runtime = rt) for ref in plot_refs]
      @test all(r -> r.status == 200 && any(h -> h.first == "Content-Type" && startswith(h.second, "image/svg+xml"), r.headers), plot_responses)
      # a chart opened from the artifact list is a Web UI page with a way
      # back; the bare route (no page=1) stays the SVG the overview embeds
      list_html = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow/artifacts/$(run_id)"; output_root = root, runtime = rt).body)
      svg_links = [m[1] for m in eachmatch(r"href=\"(/powerflow/artifact/[^\"]*\.svg[^\"]*)\"", list_html)]
      @test length(svg_links) == md["qsts_plots"]
      @test all(l -> endswith(l, "?page=1"), svg_links)
      chart_page = SparlectraApp.route_sparlectra_webui("GET", first(svg_links); output_root = root, runtime = rt)
      @test chart_page.status == 200
      @test any(h -> h.first == "Content-Type" && startswith(h.second, "text/html"), chart_page.headers)
      chart_html = String(chart_page.body)
      bare_href = first(svg_links)[1:(end - length("?page=1"))]
      @test occursin("href=\"/powerflow/artifacts/$(run_id)\"", chart_html)
      @test occursin("<object type=\"image/svg+xml\" data=\"$(bare_href)\"", chart_html)
      bare = SparlectraApp.route_sparlectra_webui("GET", bare_href; output_root = root, runtime = rt)
      @test bare.status == 200 && any(h -> h.first == "Content-Type" && startswith(h.second, "image/svg+xml"), bare.headers)
      @test startswith(lstrip(String(bare.body)), "<")  && !occursin("<html", String(bare.body))
      # a QSTS run writes no bus or branch CSV of a single power flow: the
      # operation log names the detailed CSV export skipped, never failed
      csv_events = [l for l in eachline(SparlectraApp.webui_operation_log_path(root)) if occursin(run_id, l) && occursin("\"event\":\"detailed_result_csv_export", l)]
      @test any(l -> occursin("\"event\":\"detailed_result_csv_export_skipped\"", l), csv_events)
      @test !any(l -> occursin("\"event\":\"detailed_result_csv_export_failed\"", l), csv_events)
      @test get(md, "qsts_controls_used", nothing) === true
      # L19: the phases of a QSTS run reach the operation log in run order
      qs_phases = [m[1] for m in (match(r"\"phase\":\"([a-z_]+)\"", l) for l in eachline(SparlectraApp.webui_operation_log_path(root)) if occursin(run_id, l) && occursin("\"event\":\"powerflow_phase_started\"", l)) if m !== nothing]
      qs_order = [findfirst(==(p), qs_phases) for p in ("importing_case", "running_steps", "writing_artifacts")]
      @test all(!isnothing, qs_order) && issorted(qs_order)

      # mode A without qsts.snapshot_controls: the companion file is there,
      # but its controllers are not used, and the result page says so
      snap_form = Dict{String,Any}("casefile" => series, "config_file" => config_path, "qsts_mode_choice" => "snapshots", "qsts_mode" => "snapshots")
      snap_response = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/run", snap_form; output_root = root, runtime = rt)
      snap_id = basename(only(h.second for h in snap_response.headers if h.first == "Location"))
      wait(SparlectraApp._POWERFLOW_WEBUI_JOBS[snap_id]["task"])
      snap_md = get_powerflow_result(snap_id)["metadata"]
      @test snap_md["qsts_mode"] == "snapshots" && snap_md["qsts_control_file"] == "sp_case14_day.control.yaml"
      @test get(snap_md, "qsts_controls_used", nothing) === false
      snap_page = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow/result/$(snap_id)"; output_root = root, runtime = rt).body)
      @test occursin("are not used", snap_page) && !occursin("Controllers from", snap_page)

      # a request with qsts_mode on a case without a series is refused by name
      plain_form = Dict{String,Any}("casefile" => "sp_case14.scf.json", "config_file" => config_path, "qsts_mode" => "sequential")
      plain_response = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/run", plain_form; output_root = root, runtime = rt)
      plain_id = basename(only(h.second for h in plain_response.headers if h.first == "Location"))
      wait(SparlectraApp._POWERFLOW_WEBUI_JOBS[plain_id]["task"])
      plain_result = get_powerflow_result(plain_id)
      @test plain_result["success"] === false
      @test plain_result["reason"] == "qsts_no_series"
    end
  end)() end
  return nothing
end

# The three status levels of a short-circuit run on a synthetic case with one
# placeholder machine: island 1 is a feeder with data (A) and a machine with
# a placeholder ratedS and no x''d (B), island 2 a feeder with data (C) and
# one without a current (D, skipped). Rows A and B are flagged on a default,
# rows C and D are flagged and Ik''max there is a lower bound; the run ends
# with a warning because two of four sources carry no data.
function _webui_short_circuit_status_tests()
  @testset "short-circuit status levels on a placeholder machine" begin (function ()
    mktempdir() do dir
      root = joinpath(dir, "out")
      casedir = mkpath(joinpath(dir, "cases"))
      feeder(bus, ik, name) = (mrid = name, name = name, bus = bus, maxInitialSymShCCurrent_A = ik, minInitialSymShCCurrent_A = ik, maxR1ToX1Ratio = 0.1, minR1ToX1Ratio = 0.1, maxR0ToX0Ratio = nothing, maxZ0ToZ1Ratio = nothing, ikSecond = nothing, governorSCD = nothing)
      net = Net(name = "scph", baseMVA = 100.0)
      foreach(b -> addBus!(net = net, busName = b, vn_kV = 110.0), ("A", "B", "C", "D"))
      addPIModelACLine!(net = net, fromBus = "A", toBus = "B", r_pu = 0.01, x_pu = 0.05, b_pu = 0.0, status = 1)
      addPIModelACLine!(net = net, fromBus = "C", toBus = "D", r_pu = 0.01, x_pu = 0.05, b_pu = 0.0, status = 1)
      addProsumer!(net = net, busName = "A", type = "ExternalNetworkInjection", referencePri = "A", vm_pu = 1.0, va_deg = 0.0)
      append!(net.sc_sources.external_network_injections, [feeder("A", 10_000.0, "F_A"), feeder("C", 8_000.0, "F_C"), feeder("D", nothing, "F_D")])
      push!(net.sc_sources.synchronous_machines, (mrid = "G_PH", name = "G_PH", bus = "B", satDirectSubtransX_pu = nothing, satDirectTransX_pu = nothing, r0_pu = nothing, x0_pu = nothing, r2_pu = nothing, x2_pu = nothing, earthing = nothing, ratedS_MVA = 141421.2, ratedU_kV = 110.0))
      case = exportSCF(net; file = joinpath(casedir, "scph.scf.json"), short_circuit = Dict{String,Any}("case" => "max"))
      res = SparlectraApp.start_powerflow_run(Dict{String,Any}("casefile" => case, "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => root, "short_circuit_mode" => true); case_directory = casedir)
      # the row values (status reason, lower bound per bus, placeholder
      # reason) are checked on the library result in test_short_circuit.jl;
      # here the service and the pages carry them
      @test res["status"] == "warning"
      @test res["metadata"]["sc_flagged_rows"] == 4
      @test res["metadata"]["sc_lower_bound_rows"] == 2
      # the CSV carries the row flag, the reasons and the island flag
      lines = readlines(joinpath(res["output_dir"], "short_circuit_max.csv"))
      delim = occursin(';', lines[1]) ? ';' : ','
      @test split(lines[1], delim)[end - 2:end] == ["flagged", "reasons", "lower_bound"]
      # the result page: both counts on the short-circuit card
      page = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow/result/$(SparlectraApp._webui_urlencode(res["run_id"]))"; output_root = root).body)
      card = match(r"<span class=\"summary-label\">Short circuit</span>(.*?)</div>"s, page)
      @test card !== nothing
      @test occursin("4 flagged", card.captures[1]) && occursin("2 lower bound", card.captures[1])
      # the head line of each table: first line of the printed table in
      # run.log, in the metadata, on the result page and above the table of
      # the artifact page
      head_max = "run status: warning (short_circuit_partial_defaults), 4 of 4 rows flagged, 2 rows with a lower bound"
      head_min = "run status: warning (short_circuit_partial_defaults), 4 of 4 rows flagged, 0 rows with a lower bound"
      run_log = readlines(joinpath(res["output_dir"], "run.log"))
      for (head, case) in ((head_max, "max"), (head_min, "min"))
        title = [run_log[i + 1] for i in findall(==(head), run_log)]
        @test length(title) == 1 && startswith(title[1], "Balanced 3-phase short circuit") && endswith(title[1], "case: $(case), c per IEC Table 1")
      end
      @test (res["metadata"]["sc_status_line_max"], res["metadata"]["sc_status_line_min"]) == (head_max, head_min)
      @test occursin(SparlectraApp._webui_escape(head_max), card.captures[1])
      for (name, head) in (("short_circuit_max.csv", head_max), ("short_circuit_min.csv", head_min))
        artifact_page = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow/artifact/$(SparlectraApp._webui_urlencode(res["run_id"]))/$(name)"; output_root = root).body)
        at = findfirst(SparlectraApp._webui_escape(head), artifact_page)
        @test at !== nothing && first(at) < first(findfirst("<table class=\"csv-grid\"", artifact_page))
      end
    end
  end)() end
  return nothing
end
