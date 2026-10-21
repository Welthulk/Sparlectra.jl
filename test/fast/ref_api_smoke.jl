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

# file: test/fast/ref_api_smoke.jl
# purpose: profile `fast`: the service API smoke (one
#          power flow and one state estimation through the application
#          layer), moved here unchanged from test_api.jl
using Sparlectra
using Test

function run_ref_api_smoke_tests()
    @testset "application layer smoke: one power flow and one state estimation" begin (function ()
      # the service API is the second way to every result; this stays in the
      # pull-request gate so a change to the library cannot break it unseen
      scf5 = joinpath(dirname(dirname(@__DIR__)), "data", "scf", "sp_case5.scf.json")
      meas5 = joinpath(dirname(dirname(@__DIR__)), "data", "scf", "sp_case5.measurements.csv")
      mktempdir() do tmpdir
        pf = run_sparlectra_api(casefile = scf5, config_file = Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, output_dir = joinpath(tmpdir, "pf"))
        @test pf.success
        @test pf.converged
        @test isfile(pf.result_file)
        # the Q-V verdict travels with every run: metadata line, count, run.log
        @test pf.metadata["qv_non_physical_states"] == 0
        @test pf.metadata["qv_characteristic_line"] == "Q-V characteristic: no non-physical generator states."
        @test occursin("Q-V characteristic: no non-physical", read(pf.logfile, String))
        # and names the buses of a non-physical end point: the Zeng/Chiang
        # case with the voltage-side release switched off (margin 1.0) ends
        # with three machines at Qmax above their setpoints
        zeng = joinpath(dirname(dirname(@__DIR__)), "data", "scf", "case14_zeng_p306_activeSet_A.scf.json")
        # (a case-scope key since 0.20.2: the Zeng case has a case
        # configuration file, so the margin comes as an explicit override)
        nonphys = run_sparlectra_api(casefile = zeng, config_file = Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, output_dir = joinpath(tmpdir, "nonphys"), config_overrides = Dict{String,Any}("power_flow.qlimits.reenable_v_hyst_pu" => 1.0))
        @test nonphys.converged
        @test nonphys.metadata["qv_non_physical_states"] == 3
        @test nonphys.metadata["qv_non_physical_buses"] == "2;3;6"
        @test occursin("3 non-physical generator state(s)", read(nonphys.logfile, String))
        se = start_powerflow_run(Dict{String,Any}("casefile" => scf5, "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => joinpath(tmpdir, "runs"), "se_mode" => true, "measurement_file" => meas5))
        @test se["status"] == "succeeded"
        @test se["metadata"]["run_mode"] == "se"
        # the measurements.csv artifact of an SE run is a measurement table
        # (its fixed layout, whatever output.csv_format says) and reads back
        # with readMeasurementTable! to the rows of the input table; the
        # verdict next to it follows the run's CSV format, one line per row
        se_de = start_powerflow_run(Dict{String,Any}("casefile" => scf5, "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => joinpath(tmpdir, "runs"), "se_mode" => true, "measurement_file" => meas5, "config_overrides" => Dict{String,Any}("output.csv_format" => "excel_de")))
        @test se_de["status"] == "succeeded"
        art = joinpath(se_de["output_dir"], "measurements.csv")
        @test isMeasurementTable(art)
        probe = importSCF(scf5)
        n_rows = length(readlines(meas5)) - count(l -> startswith(l, "#") || startswith(l, "timestamp"), readlines(meas5))
        @test length(readMeasurementTable(art)) == n_rows
        @test isempty(readMeasurementTable!(probe; file = art, zero_injection = false).problems)
        verdict = readlines(joinpath(se_de["output_dir"], "measurements_verdict.csv"))
        @test verdict[1] == "line;key;type;measurement_id;residual;normalized_residual;wii;localizable;eliminated;constraint;reason"
        @test length(verdict) == n_rows + 1
        @test any(l -> occursin(r";-?\d+,\d+;", l), verdict)
        # a gross error (the one row of the shipped bad-data set whose value
        # differs from the clean set, scaled by 10): the verdict names that
        # row as the one eliminated, with the normalized residual it was
        # eliminated at
        bad5 = joinpath(dirname(meas5), "sp_case5.measurements.baddata.csv")
        se_bad = start_powerflow_run(Dict{String,Any}("casefile" => scf5, "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => joinpath(tmpdir, "runs"), "se_mode" => true, "measurement_file" => bad5))
        @test se_bad["status"] == "succeeded"
        vcsv = run_csv(joinpath(se_bad["output_dir"], "measurements_verdict.csv"))
        vbad = vcsv.rows
        elim = filter(v -> v[9] == "true", vbad)
        @test length(elim) == se_bad["metadata"]["se_eliminations"] == 1
        art_bad = readlines(joinpath(se_bad["output_dir"], "measurements.csv"))
        gross = findall(((a, b),) -> a.value != b.value, collect(zip(readMeasurementTable(meas5), readMeasurementTable(bad5))))
        @test length(gross) == 1
        @test Int(vcsv.number(elim[1][1])) == findfirst(l -> startswith(l, "timestamp"), art_bad) + only(gross)
        @test abs(vcsv.number(elim[1][6])) >= se_bad["metadata"]["se_k_eliminate"]
        @test elim[1][11] == "eliminated"
        @test occursin(elim[1][4], read(joinpath(se_bad["output_dir"], "se_diagnostics.md"), String))
      end
    end)() end
  return true
end
