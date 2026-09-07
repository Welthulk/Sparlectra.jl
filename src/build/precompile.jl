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
# purpose: PrecompileTools workload covering the RUN PATH so first
#          interactive runs skip most JIT compilation (issue #288). Since
#          adapter stage 3 the run path is import_case -> convert_case ->
#          build_net, so the workload warms exactly that for every
#          adapter with a tracked small case, plus the rectangular solve,
#          loss postprocessing and the standalone DC power flow.

# Warmed per adapter:
# - MATPOWER: data/mpower/warmup_casePST.m (tracked), through import_case,
#   then the config-form rectangular solve and the loss postprocessing.
#   The DC power flow runs on the DTF net below: the PST warmup case
#   splits into two islands under the DC island analysis (island 2 has no
#   reference) and would abort precompilation.
# - SCF: data/scf/sp_casePST.scf.json (tracked), through import_case.
# - PGM: data/scf/pgm_interop.json (tracked), same pipeline as SCF.
# - DTF: a tiny in-memory DTFCase (the synthetic shape the DTF tests use)
#   through convert_case and build_net; the format has no tracked file
#   fixture.
# - CGMES: NOT warmed; the mapping has no tracked fixture. Its build side
#   is the shared build_net path warmed above; the mapping itself stays
#   cold by design and is named here so nobody mistakes it for covered.
#
# isfile guards keep precompilation robust if the data directory is
# stripped; stdout is silenced because the import context prints run
# banners.

using PrecompileTools: @setup_workload, @compile_workload

@setup_workload begin
  _pc_mpower = normpath(joinpath(@__DIR__, "..", "..", "data", "mpower", "warmup_casePST.m"))
  _pc_scf = normpath(joinpath(@__DIR__, "..", "..", "data", "scf", "sp_casePST.scf.json"))
  _pc_pgm = normpath(joinpath(@__DIR__, "..", "..", "data", "scf", "pgm_interop.json"))
  @compile_workload begin
    Logging.with_logger(Logging.NullLogger()) do
      redirect_stdout(devnull) do
        _pc_cfg = load_sparlectra_config(DEFAULT_SPARLECTRA_CONFIG_PATH; reload = true)
        if isfile(_pc_mpower)
          _pc_imported = import_case(_pc_mpower, _pc_cfg)
          runpf!(_pc_imported.net; config = _pc_imported.config)
          calcNetLosses!(_pc_imported.net)
        end
        isfile(_pc_scf) && import_case(_pc_scf, _pc_cfg)
        isfile(_pc_pgm) && import_case(_pc_pgm, _pc_cfg)
        _pc_dtf = DTFImporter.DTFCase(
          "precompile",
          100.0,
          DTFImporter.DTFParams("", Float64[]),
          ["precompile"],
          [110.0],
          DTFImporter.DTFSize("", 2, 1, 0, 0, "SLACK"),
          [DTFImporter.DTFBranch("", 1, 'T', 1, "A", "PV", "SLACK", 1.21, 6.05, 4.0e-5, -1.0e-5, nothing)],
          DTFImporter.DTFCompensation[],
          DTFImporter.DTFTransformerControl[],
          [
            DTFImporter.DTFBus("", 1, 1, 1, "PV", 110.0, 0.0, 0.0, 0.0, 10.0, 2.0, -5.0, 5.0),
            DTFImporter.DTFBus("", 2, 2, 1, "SLACK", 110.0, 0.0, 0.0, 0.0, 20.0, 3.0, -10.0, 10.0),
          ],
          DTFImporter.DTFOutage[],
          DTFImporter.DTFTrailingRecord[],
        )
        # task_import_direct: the run path is the DIRECT importer, so that
        # is what gets warm; the converter stays a product of its own and
        # is warmed separately (the explicit SCF export path)
        _pc_dtf_net = DTFImporter.build_net(_pc_dtf; bus_shunt_model = _pc_cfg.model.bus_shunt_model, tap_changer_model = _pc_cfg.model.tap_changer_model)
        _apply_config_net_parameters!(_pc_dtf_net, _pc_cfg)
        convert_case(DTFAdapter(), _pc_dtf, dtf_adapter_options(_pc_cfg))
        rundcpf!(_pc_dtf_net)
      end
    end
  end
end
