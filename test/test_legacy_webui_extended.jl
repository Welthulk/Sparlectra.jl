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

# file: test/test_legacy_webui_extended.jl
# purpose: legacy profile: the testsets of test/test_webui_extended.jl that pin an earlier
#          change (an issue, a bug, a bit-identity or off-path guarantee).
#          They moved here unchanged from test_webui_extended.jl; the helpers they use
#          are defined in that file, which the legacy group includes first.
using Sparlectra
using Test

function run_legacy_webui_extended_tests()
  # cut from testset: Local PowerFlow Web UI
    @testset "DC power flow mode" begin (function ()
      mktempdir() do dc_tmpdir
        dc_casefile = _write_webui_test_case(joinpath(dc_tmpdir, "case_webui_dc.m"))
        dc_config_file = joinpath(dc_tmpdir, "config_template.yaml")
        cp(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, dc_config_file)
        dc_output_root = joinpath(dc_tmpdir, "runs")

        # (2a) Rendered form contains the DC selection, with AC selected by default.
        # the solver radio group lives on the Settings page
        rendered_form = String(
          SparlectraApp.route_sparlectra_webui(
            "GET",
            "/powerflow/settings?casefile=$(SparlectraApp._webui_urlencode(dc_casefile))&config_file=$(SparlectraApp._webui_urlencode(dc_config_file))";
            output_root = dc_output_root,
          ).body,
        )
        # AC/APSLF/DC are one unified radio group on the single "power_flow_solver"
        # field -- there is no separate "Berechnungsmodell" field/control anymore,
        # and no dropdown duplicating the "dc" choice.
        @test occursin("<input type=\"radio\" name=\"power_flow_solver\" value=\"rectangular\" data-solver-radio checked>", rendered_form)
        @test occursin("<input type=\"radio\" name=\"power_flow_solver\" value=\"apslf\" data-solver-radio>", rendered_form)
        @test occursin("<input type=\"radio\" name=\"power_flow_solver\" value=\"dc\" data-solver-radio>", rendered_form)
        @test !occursin("<input type=\"radio\" name=\"power_flow_solver\" value=\"apslf\" data-solver-radio checked>", rendered_form)
        @test !occursin("<input type=\"radio\" name=\"power_flow_solver\" value=\"dc\" data-solver-radio checked>", rendered_form)
        @test occursin("DC (linear screening model, replaces Newton-Raphson entirely)", rendered_form)
        @test !occursin("Berechnungsmodell", rendered_form)
        @test !occursin("power_flow_calc_mode", rendered_form)
        @test !occursin("data-calc-mode-radio", rendered_form)
        @test !occursin("data-solver-select", rendered_form)
        @test !occursin("<select name=\"power_flow_solver\"", rendered_form)
        @test occursin("data-ac-only-field", rendered_form)
        # Newton-Raphson's DC-derived start angles are documented inline next to the
        # AC radio, distinct from the standalone DC solver option below it.
        @test occursin("Start angle mode", rendered_form)
        @test occursin("DC pre-pass", rendered_form)

        # (2b) POST with DC mode terminates successfully; result.json carries the DC flag.
        dc_form = _webui_test_form(dc_casefile, dc_config_file, dc_output_root)
        dc_form["power_flow_solver"] = "dc"
        dc_form["detailed_result_csv_format"] = "technical"
        dc_response = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/run", dc_form; output_root = dc_output_root)
        @test dc_response.status == 303
        dc_run_id = basename(only(header.second for header in dc_response.headers if header.first == "Location"))
        wait(SparlectraApp._POWERFLOW_WEBUI_JOBS[dc_run_id]["task"])
        dc_result = get_powerflow_result(dc_run_id)
        @test dc_result["success"]
        dc_result_json = SparlectraApp._parse_service_json(read(joinpath(dc_result["output_dir"], "result.json"), String))
        @test dc_result_json["final_outcome"]["solver"] == "dc"

        # (2c) Result page marks the DC run distinctly.
        dc_result_response = SparlectraApp.handle_powerflow_result(dc_run_id)
        @test dc_result_response.status == 200
        dc_result_html = String(dc_result_response.body)
        @test occursin("DC solution", dc_result_html)
        @test occursin("status-info", dc_result_html)
        @test occursin("<code>dc</code>", dc_result_html)

        # (2d) Artifact listing of the DC run is not empty; result.json downloads with status 200.
        dc_artifacts = list_powerflow_artifacts(dc_run_id)
        @test !isempty(dc_artifacts)
        @test any(artifact -> artifact["name"] == "result.json", dc_artifacts)
        @test any(artifact -> artifact["name"] == "bus_voltages_complex.csv", dc_artifacts)
        @test SparlectraApp.route_sparlectra_webui("GET", "/powerflow/artifact/$(dc_run_id)/result.json?download=1"; output_root = dc_output_root).status == 200

        # (2e) Regression: no solver field, and an explicit AC selection, still
        # produce the unchanged AC report (no "dc" solver in final_outcome).
        ac_implicit_form = _webui_test_form(dc_casefile, dc_config_file, dc_output_root)
        ac_implicit_response = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/run", ac_implicit_form; output_root = dc_output_root)
        @test ac_implicit_response.status == 303
        ac_implicit_run_id = basename(only(header.second for header in ac_implicit_response.headers if header.first == "Location"))
        wait(SparlectraApp._POWERFLOW_WEBUI_JOBS[ac_implicit_run_id]["task"])
        ac_implicit_result = get_powerflow_result(ac_implicit_run_id)
        @test ac_implicit_result["success"]
        ac_implicit_json = SparlectraApp._parse_service_json(read(joinpath(ac_implicit_result["output_dir"], "result.json"), String))
        @test ac_implicit_json["final_outcome"]["solver"] == "rectangular"

        ac_explicit_form = _webui_test_form(dc_casefile, dc_config_file, dc_output_root)
        ac_explicit_form["power_flow_solver"] = "rectangular"
        ac_explicit_response = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/run", ac_explicit_form; output_root = dc_output_root)
        @test ac_explicit_response.status == 303
        ac_explicit_run_id = basename(only(header.second for header in ac_explicit_response.headers if header.first == "Location"))
        wait(SparlectraApp._POWERFLOW_WEBUI_JOBS[ac_explicit_run_id]["task"])
        ac_explicit_result = get_powerflow_result(ac_explicit_run_id)
        @test ac_explicit_result["success"]
        ac_explicit_json = SparlectraApp._parse_service_json(read(joinpath(ac_explicit_result["output_dir"], "result.json"), String))
        @test ac_explicit_json["final_outcome"]["solver"] == "rectangular"
        @test ac_implicit_json["final_outcome"]["converged"] == ac_explicit_json["final_outcome"]["converged"]
        @test ac_implicit_json["final_outcome"]["iterations"] == ac_explicit_json["final_outcome"]["iterations"]

        # History row shows the calculation mode for the DC run, unchanged status columns otherwise.
        dc_runs = SparlectraApp.list_powerflow_runs(dc_output_root)
        dc_history_html = SparlectraApp.render_powerflow_history(dc_runs, dc_output_root)
        @test occursin("<th>Solver</th>", dc_history_html)
        @test occursin("<td>dc</td>", dc_history_html)
        @test occursin("<td>rectangular</td>", dc_history_html)

        # (3) Plausibility: the WebUI DC run's bus angles match a direct rundcpf! call on the
        # same case, proving both paths run the same solver code.
        dc_bus_csv = read(joinpath(dc_result["output_dir"], "bus_voltages_complex.csv"), String)
        dc_bus_rows = split(strip(dc_bus_csv), '\n')[2:end]
        webui_va_deg = Dict{Int,Float64}()
        for row in dc_bus_rows
          fields = split(row, ',')
          webui_va_deg[parse(Int, fields[1])] = parse(Float64, fields[5])
        end
        mpc = Sparlectra.MatpowerIO.read_case(dc_casefile)
        direct_net = Sparlectra.createNetFromMatPowerCase(mpc = mpc, flatstart = true)
        direct_report = rundcpf!(direct_net)
        for row in direct_report.nodes
          @test isapprox(webui_va_deg[row.bus], row.va_deg; atol = 1e-8)
        end
      end
    end)() end

    @testset "Diagnose button submitter preservation" begin (function ()
      mktempdir() do diag_tmpdir
        diag_casefile = _write_webui_test_case(joinpath(diag_tmpdir, "case_webui_diag.m"))
        diag_config_file = joinpath(diag_tmpdir, "config_template.yaml")
        cp(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, diag_config_file)
        diag_output_root = joinpath(diag_tmpdir, "runs")

        # (2a) The rendered form no longer disables submit buttons synchronously
        # before serialization; it preserves the submitter's name/value instead.
        diag_form_html = String(
          SparlectraApp.route_sparlectra_webui(
            "GET",
            "/powerflow?casefile=$(SparlectraApp._webui_urlencode(diag_casefile))&config_file=$(SparlectraApp._webui_urlencode(diag_config_file))";
            output_root = diag_output_root,
          ).body,
        )
        @test !occursin("onsubmit=\"this.classList.add('is-submitting')", diag_form_html)
        @test occursin("<form id=\"powerflow-run-form\" data-powerflow-form method=\"post\" action=\"/powerflow/run\" class=\"panel form-grid powerflow-form-card\">", diag_form_html)
        @test occursin("const copySubmitterValue = function (submitter)", diag_form_html)
        @test occursin("data-submitter-value", diag_form_html)
        @test occursin("event.submitter", diag_form_html)
        @test occursin("powerflowForm.addEventListener('submit', function (event)", diag_form_html)
        @test occursin("button.addEventListener('click', function () { copySubmitterValue(button); })", diag_form_html)
        @test occursin("powerflowForm.classList.add('is-submitting')", diag_form_html)
        @test occursin("powerflowForm.setAttribute('aria-busy', 'true')", diag_form_html)
        @test occursin("powerflowForm.querySelectorAll('button[type=submit]').forEach(function (b) { b.disabled = true; })", diag_form_html)
        @test occursin("staleSubmitterValue", diag_form_html)

        # (1) Repro (still valid post-fix, since it exercises the unchanged server
        # path): a POST body without "diagnose_mode" never writes diagnose.log.
        no_diagnose_form = _webui_test_form(diag_casefile, diag_config_file, diag_output_root)
        haskey(no_diagnose_form, "diagnose_mode") && delete!(no_diagnose_form, "diagnose_mode")
        no_diagnose_response = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/run", no_diagnose_form; output_root = diag_output_root)
        @test no_diagnose_response.status == 303
        no_diagnose_run_id = basename(only(header.second for header in no_diagnose_response.headers if header.first == "Location"))
        wait(SparlectraApp._POWERFLOW_WEBUI_JOBS[no_diagnose_run_id]["task"])
        no_diagnose_result = get_powerflow_result(no_diagnose_run_id)
        @test no_diagnose_result["success"]
        @test !isfile(joinpath(no_diagnose_result["output_dir"], "diagnose.log"))

        # (2b) POST with diagnose_mode=true writes a non-empty diagnose.log containing
        # the "Diagnosis" section from run_diagnostic_artifacts.jl.
        diagnose_form = _webui_test_form(diag_casefile, diag_config_file, diag_output_root)
        diagnose_form["diagnose_mode"] = "true"
        diagnose_response = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/run", diagnose_form; output_root = diag_output_root)
        @test diagnose_response.status == 303
        diagnose_run_id = basename(only(header.second for header in diagnose_response.headers if header.first == "Location"))
        wait(SparlectraApp._POWERFLOW_WEBUI_JOBS[diagnose_run_id]["task"])
        diagnose_result = get_powerflow_result(diagnose_run_id)
        # A diagnose run takes ONE step from the case's own stored VM/VA, so
        # a residual is the normal outcome and `success` (which mirrors
        # convergence) is false: the service completed, the numerics did not
        # converge, and that residual IS the diagnosis. This asserted
        # `success` before this fix, when the forced settings still lost to
        # the case configuration file and the run was secretly an ordinary
        # solve.
        @test diagnose_result["success"] === false
        @test diagnose_result["service_status"] == "completed"
        @test diagnose_result["run_status"] == "completed_nonconverged"
        @test diagnose_result["iterations"] == 1
        diagnose_log_path = joinpath(diagnose_result["output_dir"], "diagnose.log")
        @test isfile(diagnose_log_path)
        @test filesize(diagnose_log_path) > 0
        diagnose_log_text = read(diagnose_log_path, String)
        @test occursin("Diagnosis", diagnose_log_text)

        # (2c) Artifact listing of the diagnostic run includes diagnose.log.
        diagnose_artifacts = list_powerflow_artifacts(diagnose_run_id)
        @test any(artifact -> artifact["name"] == "diagnose.log", diagnose_artifacts)

        # (2d) diagnose.log downloads with status 200 and the same file content.
        diagnose_download = SparlectraApp.route_sparlectra_webui("GET", "/powerflow/artifact/$(diagnose_run_id)/diagnose.log?download=1"; output_root = diag_output_root)
        @test diagnose_download.status == 200
        @test String(diagnose_download.body) == diagnose_log_text

        # Regression: a plain "Start PowerFlow run" (no diagnose_mode field at all,
        # matching what the real submit button sends) still never writes diagnose.log.
        normal_form = _webui_test_form(diag_casefile, diag_config_file, diag_output_root)
        haskey(normal_form, "diagnose_mode") && delete!(normal_form, "diagnose_mode")
        normal_response = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/run", normal_form; output_root = diag_output_root)
        @test normal_response.status == 303
        normal_run_id = basename(only(header.second for header in normal_response.headers if header.first == "Location"))
        wait(SparlectraApp._POWERFLOW_WEBUI_JOBS[normal_run_id]["task"])
        normal_result = get_powerflow_result(normal_run_id)
        @test normal_result["success"]
        @test !isfile(joinpath(normal_result["output_dir"], "diagnose.log"))
      end
    end)() end

  return true
end
