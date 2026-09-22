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

# file: src/build/precompile.jl
# purpose: PrecompileTools workload covering the RUN PATHS so first
#          interactive runs skip most JIT compilation (issue #288). Four
#          paths are warmed: the file-based import (MATPOWER, SCF, PGM)
#          with the rectangular solve and loss postprocessing; the
#          programmatic path of the workshop notebooks (builder API,
#          run_sparlectra with the rectangular and the APSLF solver, the
#          status query, the APSLF start ahead of NR); state estimation on
#          the tracked SCF fixture, which carries a measurement set; and
#          the standalone DC power flow.

# Warmed:
# - MATPOWER: data/mpower/warmup_casePST.m (tracked), through import_case,
#   then the config-form rectangular solve and the loss postprocessing.
# - SCF: data/scf/sp_casePST.scf.json (tracked), through import_case, and
#   again through importSCF for the state estimation (the fixture ships
#   59 measurements, so runse! runs on the file alone).
# - PGM: data/scf/pgm_interop.json (tracked), same pipeline as SCF.
# - Programmatic: the 7-bus ring of the APSLF workshop, built with the
#   addBus!/addPIModelACLine!/addProsumer! API and solved through
#   run_sparlectra(net = ..., config = ...) with both solvers and the
#   hybrid start. The APSLF calls also pull the AnalyticLoadFlow methods
#   instantiated with Sparlectra's concrete types into this package image;
#   ALF's own workload only covers its demo types. The ring is known to
#   solve under both solvers (a 3-bus ring did not converge under APSLF).
#   The same ring feeds rundcpf!: the PST warmup case splits into two
#   islands under the DC island analysis (island 2 has no reference) and
#   would abort precompilation.
# - NOT warmed: DTF and CGMES. Neither has a tracked file fixture; their
#   build side is the shared build_net path warmed above, the mappings
#   themselves stay cold by design and are named here so nobody mistakes
#   them for covered.
#
# isfile guards keep precompilation robust if the data directory is
# stripped; stdout is silenced because the import context prints run
# banners, the logger because the SE topology precheck warns on the
# fixture.

using PrecompileTools: @setup_workload, @compile_workload

@setup_workload begin
    _pc_mpower = normpath(joinpath(@__DIR__, "..", "..", "data", "mpower", "warmup_casePST.m"))
    _pc_scf = normpath(joinpath(@__DIR__, "..", "..", "data", "scf", "sp_casePST.scf.json"))
    _pc_pgm = normpath(joinpath(@__DIR__, "..", "..", "data", "scf", "pgm_interop.json"))
    @compile_workload begin
        Logging.with_logger(Logging.NullLogger()) do
            redirect_stdout(devnull) do
                # --- file-based import + rectangular solve ------------------------
                _pc_cfg = load_sparlectra_config(DEFAULT_SPARLECTRA_CONFIG_PATH; reload=true)
                if isfile(_pc_mpower)
                    _pc_imported = import_case(_pc_mpower, _pc_cfg)
                    runpf!(_pc_imported.net; config=_pc_imported.config)
                    calcNetLosses!(_pc_imported.net)
                end
                isfile(_pc_scf) && import_case(_pc_scf, _pc_cfg)
                isfile(_pc_pgm) && import_case(_pc_pgm, _pc_cfg)

                # --- state estimation on the tracked SCF fixture ------------------
                if isfile(_pc_scf)
                    _pc_se_net = importSCF(_pc_scf)
                    runse!(_pc_se_net)
                end

                # --- programmatic net (workshop ring7) + run_sparlectra -----------
                _pc_net = Net(name="precompile_ring7", baseMVA=100.0)
                addBus!(net=_pc_net, busName="B1", vn_kV=110.0, vm_pu=1.02, va_deg=0.0)
                for _pc_i in 2:7
                    addBus!(net=_pc_net, busName="B$(_pc_i)", vn_kV=110.0, vm_pu=1.0, va_deg=0.0)
                end
                addPIModelACLine!(net=_pc_net, fromBus="B1", toBus="B2", r_pu=0.010, x_pu=0.080, b_pu=0.0, status=1)
                addPIModelACLine!(net=_pc_net, fromBus="B2", toBus="B3", r_pu=0.011, x_pu=0.085, b_pu=0.0, status=1)
                addPIModelACLine!(net=_pc_net, fromBus="B3", toBus="B4", r_pu=0.012, x_pu=0.090, b_pu=0.0, status=1)
                addPIModelACLine!(net=_pc_net, fromBus="B4", toBus="B5", r_pu=0.010, x_pu=0.080, b_pu=0.0, status=1)
                addPIModelACLine!(net=_pc_net, fromBus="B5", toBus="B6", r_pu=0.011, x_pu=0.085, b_pu=0.0, status=1)
                addPIModelACLine!(net=_pc_net, fromBus="B6", toBus="B7", r_pu=0.012, x_pu=0.090, b_pu=0.0, status=1)
                addPIModelACLine!(net=_pc_net, fromBus="B7", toBus="B1", r_pu=0.010, x_pu=0.080, b_pu=0.0, status=1)
                addPIModelACLine!(net=_pc_net, fromBus="B2", toBus="B5", r_pu=0.009, x_pu=0.070, b_pu=0.0, status=1)
                addPIModelACLine!(net=_pc_net, fromBus="B3", toBus="B6", r_pu=0.009, x_pu=0.070, b_pu=0.0, status=1)
                addProsumer!(net=_pc_net, busName="B1", type="EXTERNALNETWORKINJECTION", referencePri="B1", vm_pu=1.02, va_deg=0.0)
                addProsumer!(net=_pc_net, busName="B3", type="GENERATOR", p=60.0, q=10.0)
                addProsumer!(net=_pc_net, busName="B2", type="LOAD", p=35.0, q=10.0)
                addProsumer!(net=_pc_net, busName="B4", type="LOAD", p=45.0, q=15.0)
                addProsumer!(net=_pc_net, busName="B5", type="LOAD", p=25.0, q=8.0)
                addProsumer!(net=_pc_net, busName="B6", type="LOAD", p=30.0, q=10.0)
                addProsumer!(net=_pc_net, busName="B7", type="LOAD", p=20.0, q=6.0)
                validate!(net=_pc_net)

                _pc_quiet = OutputConfig(logfile_results=:off, console_summary=false, startup_latency_hint=false)
                _pc_cfg_nr = SparlectraConfig(powerflow=PowerFlowConfig(solver=:rectangular, rescue=false), output=_pc_quiet)
                _pc_cfg_ap = SparlectraConfig(powerflow=PowerFlowConfig(solver=:apslf), output=_pc_quiet)
                _pc_cfg_hyb = SparlectraConfig(powerflow=PowerFlowConfig(solver=:rectangular, apslf_start=ApslfStartConfig(enabled=true)), output=_pc_quiet)

                _pc_r_nr = run_sparlectra(net=deepcopy(_pc_net), config=_pc_cfg_nr)
                _pc_r_ap = run_sparlectra(net=deepcopy(_pc_net), config=_pc_cfg_ap)
                run_sparlectra(net=deepcopy(_pc_net), config=_pc_cfg_hyb)
                rectangular_pf_status(_pc_r_nr.net)
                _pc_st = rectangular_pf_status(_pc_r_ap.net)
                _pc_st.apslf_convergence_line
                _pc_r_ap.final_converged
                _pc_r_ap.final_mismatch
                [n._vm_pu for n in _pc_r_ap.net.nodeVec]

                # --- standalone DC power flow -------------------------------------
                rundcpf!(deepcopy(_pc_net))

                # SPARLECTRA_PRECOMPILE_WORKLOAD=minimal skips the two blocks below:
                # a notebook session (Colab installs and precompiles from scratch
                # every time) never calls the service layer or a tap controller,
                # and those two cost about 40 s of precompile time
                _pc_full = get(ENV, "SPARLECTRA_PRECOMPILE_WORKLOAD", "default") != "minimal"

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
                if _pc_full && isfile(_pc_scf)
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
                    # NOT warmed: the Web UI pages (maintainer decision
                    # 2026-09-22: rendering them here cost another 20 s of
                    # precompile time; the Runs page stays cold on its first
                    # open, the sysimage covers that for the maintainer)
                    lock(_POWERFLOW_SERVICE_LOCK) do
                        delete!(_POWERFLOW_SERVICE_RUNS, _pc_run_id)
                        delete!(_POWERFLOW_SERVICE_RUNS, String(_pc_se["run_id"]))
                    end
                    rm(_pc_out; recursive=true, force=true)
                end

                # --- tap controller path ------------------------------------------
                # sp_case14 carries a declared OLTC controller: the control loop
                # (outer passes, controller write-back) is not on the ring's path
                # and cost 2.3 s on its first run
                _pc_scf14 = normpath(joinpath(@__DIR__, "..", "..", "data", "scf", "sp_case14.scf.json"))
                _pc_full && isfile(_pc_scf14) && run_sparlectra(net=importSCF(_pc_scf14), config=_pc_cfg_nr)
            end
        end
    end
end