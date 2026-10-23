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

        # the only run of solver=dc through the service: the POST terminates,
        # result.json carries the DC flag
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

        # the artifact listing of the DC run is not empty; result.json downloads
        dc_artifacts = list_powerflow_artifacts(dc_run_id)
        @test !isempty(dc_artifacts)
        @test any(artifact -> artifact["name"] == "result.json", dc_artifacts)
        @test any(artifact -> artifact["name"] == "bus_voltages_complex.csv", dc_artifacts)
        @test SparlectraApp.route_sparlectra_webui("GET", "/powerflow/artifact/$(dc_run_id)/result.json?download=1"; output_root = dc_output_root).status == 200

        # the bus angles of the service run match a direct rundcpf! on the
        # same case, so both paths run the same solver code
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

  return true
end
