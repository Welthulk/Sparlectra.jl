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
# purpose: precompile workload of the application layer: the Web UI start
#          page in every mode (a start without a sysimage then shows the
#          page in seconds instead of compiling it for a minute or two),
#          and the service run behind the Web UI's first click only when
#          SPARLECTRA_PRECOMPILE_WORKLOAD=full is set (the sysimage build
#          sets it); a library session never pays for either.

using PrecompileTools: @setup_workload, @compile_workload

# a function, not code inside the block: see src/build/precompile.jl of the
# library for why (a guarded block is compiled before it runs)
function _precompile_app_workload()
    _pc_scf = normpath(joinpath(Sparlectra.SPARLECTRA_ROOT, "data", "scf", "sp_casePST.scf.json"))
    Logging.with_logger(Logging.NullLogger()) do
        redirect_stdout(devnull) do
            # --- service path (the Web UI's first run) --------------------------
            # Measured on a fresh process with the workload above:
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
                # the start page is warmed in every mode by
                # _precompile_webui_start_page (the Settings page shares most of
                # its code); the run result pages stay cold on their first open,
                # the sysimage covers those
                lock(_POWERFLOW_SERVICE_LOCK) do
                    delete!(_POWERFLOW_SERVICE_RUNS, _pc_run_id)
                    delete!(_POWERFLOW_SERVICE_RUNS, String(_pc_se["run_id"]))
                end
                rm(_pc_out; recursive=true, force=true)
            end
        end
    end
    return nothing
end

# The start page of the Web UI, rendered the way the server renders it for a
# GET request: a real runtime (a listener on a port the system picks, a case
# directory with one shipped demo case, a copy of the template), the request
# line split into SubStrings like `_webui_read_request` does. Runs in every
# workload mode: without a sysimage the first request of the start page
# otherwise compiles its whole render path on every start of the Web UI;
# here it is paid once, when the package image is built. Everything lives in
# a temporary directory, nothing of the user's state is read or written.
function _precompile_webui_start_page()
    dir = mktempdir()
    listener = nothing
    try
        cases = joinpath(dir, "cases")
        mkpath(cases)
        demo = normpath(joinpath(Sparlectra.SPARLECTRA_ROOT, "data", "scf", "sp_case5.scf.json"))
        isfile(demo) && cp(demo, joinpath(cases, "sp_case5.scf.json"))
        config = joinpath(dir, "configuration.yaml")
        cp(DEFAULT_SPARLECTRA_CONFIG_PATH, config)
        listener = Sockets.listen(Sockets.localhost, 0)
        runtime = _SparlectraWebUIRuntime(listener, cases, config, joinpath(dir, "webui_operations.jsonl"), nothing, start_powerflow_run, false, false, 0.0, 0, nothing, devnull, ReentrantLock())
        method, target = split("GET /powerflow")
        Logging.with_logger(Logging.NullLogger()) do
            redirect_stdout(devnull) do
                route_sparlectra_webui(method, target, Dict{String,String}(); output_root = joinpath(dir, "runs"), runtime)
            end
        end
    finally
        listener === nothing || close(listener)
        rm(dir; recursive = true, force = true)
    end
    return nothing
end

# The workload mode this image was built with, baked in at precompile time
# like `Sparlectra.PRECOMPILE_WORKLOAD_MODE`; the gate runner reads it.
const PRECOMPILE_WORKLOAD_MODE = get(ENV, "SPARLECTRA_PRECOMPILE_WORKLOAD", "off")

@setup_workload begin
    _pc_full = get(ENV, "SPARLECTRA_PRECOMPILE_WORKLOAD", "off") == "full"
    @compile_workload begin
        _precompile_webui_start_page()
        _pc_full && _precompile_app_workload()
    end
end
