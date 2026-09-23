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
# file: app/src/precompile.jl
# purpose: precompile workload of the application layer: the service run
#          behind the Web UI's first click, warmed only when
#          SPARLECTRA_PRECOMPILE_WORKLOAD=full is set (the sysimage build
#          sets it); a library session never pays for it.

using PrecompileTools: @setup_workload, @compile_workload

@setup_workload begin
    _pc_scf = normpath(joinpath(Sparlectra.SPARLECTRA_ROOT, "data", "scf", "sp_casePST.scf.json"))
    _pc_full = get(ENV, "SPARLECTRA_PRECOMPILE_WORKLOAD", "off") == "full"
    @compile_workload begin
        Logging.with_logger(Logging.NullLogger()) do
            redirect_stdout(devnull) do
                if _pc_full
                    # --- service path (the Web UI's first run) --------------------------
                    # Measured 2026-09-22 on a fresh process with the workload above:
                    # the solver paths answered in well under a second, but the first
                    # start_powerflow_run took 38 s (25 s in run_sparlectra_api: the
                    # artifact writers, effective_config.yaml, result.json, the
                    # metadata; 22 s more in the service layer: run id, run index,
                    # lifecycle). That is exactly the first click on the Runs page,
                    # so both layers are warmed here on the tracked SCF fixture.
                    # The run's entry in the process-wide registry is removed again:
                    # nothing of this run may survive in the package image.
                    if isfile(_pc_scf)
                        _pc_out = mktempdir()
                        _pc_api = run_sparlectra_api(casefile=_pc_scf, config_file=DEFAULT_SPARLECTRA_CONFIG_PATH, output_dir=joinpath(_pc_out, "api"))
                        to_dict(_pc_api)
                        _pc_srv = start_powerflow_run(Dict{String,Any}("casefile" => _pc_scf, "config_file" => DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => joinpath(_pc_out, "runs")))
                        _pc_run_id = String(_pc_srv["run_id"])
                        get_powerflow_result(_pc_run_id)
                        list_powerflow_artifacts(_pc_run_id)
                        # the state-estimation run of the same service (the fixture
                        # carries its measurement set): 2.7 s cold without this
                        _pc_se = start_powerflow_run(Dict{String,Any}("casefile" => _pc_scf, "config_file" => DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => joinpath(_pc_out, "runs"), "se_mode" => true, "measurement_file" => ""))
                        # NOT warmed: the Web UI pages (rendering them here cost another
                        # 20 s of precompile time; the Runs page stays cold on its
                        # first open, the sysimage covers that)
                        lock(_POWERFLOW_SERVICE_LOCK) do
                            delete!(_POWERFLOW_SERVICE_RUNS, _pc_run_id)
                            delete!(_POWERFLOW_SERVICE_RUNS, String(_pc_se["run_id"]))
                        end
                        rm(_pc_out; recursive=true, force=true)
                    end
                end
            end
        end
    end
end
