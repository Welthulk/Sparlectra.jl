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

# file: src/adapters/cgmes/converter.jl
# purpose: the CGMES adapter (adapter task stage 3c): convert_case runs the
#          proven configured import and captures the mapped network as the
#          typed SCFCase. The source-system identity travels with the case:
#          the full structural-key mRID registry in the namespaced meta and
#          the per-component mRID under the `cgmes` key of its extra
#          record, so a case built from the file still drives
#          writeCGMESFiles and the mRID-addressed measurement resolution;
#          the CGMESContext is derived at export time from the net and the
#          restored registry, so nothing else has to persist. The SV
#          voltages arrive as the case's start_state (the import applies
#          them per the start-values decision before the capture).

struct CGMESAdapter <: FormatAdapter end

"""
    CGMESAdapterOptions

The adapter-scope options of the CGMES conversion (design decision D4):
the `cgmes_import` configuration surface plus the run kind the start-value
decision may depend on.
"""
Base.@kwdef struct CGMESAdapterOptions
  start_values::Symbol = :auto
  require_boundary::Bool = true
  infer_base_voltages::Bool = false
  hvdc_mode::Symbol = :injections
  run_kind::Symbol = :powerflow
  name::String = ""
end

detect(::Type{CGMESAdapter}, path::AbstractString)::Bool = begin
  isdir(path) && return true
  ext = lowercase(splitext(String(path))[2])
  ext in (".zip", ".xml") || return false
  return isfile(path)
end

options_type(::CGMESAdapter) = CGMESAdapterOptions

"""
    import_net(::CGMESAdapter, path::AbstractString, opts::CGMESAdapterOptions; config = active_sparlectra_config()) -> Net

The CGMES importer of the adapter contract (task_import_direct): the
delivery builds natively through importCGMES. The method maps the
adapter options that the importer consumes; the SV start-value
decision and the run-kind gating live in the config-rich service
layer (_import_cgmes carries the thirteen cgmes_import values), which
stays the path the services use.
"""
function import_net(::CGMESAdapter, path::AbstractString, opts::CGMESAdapterOptions; config::SparlectraConfig = active_sparlectra_config())::Net
  net = createNetFromCGMES(path = String(path), require_boundary = opts.require_boundary, infer_base_voltages = opts.infer_base_voltages, hvdc_mode = opts.hvdc_mode, bus_shunt_model = config.model.bus_shunt_model, name = isempty(opts.name) ? basename(String(path)) : opts.name)
  _apply_config_net_parameters!(net, config)
  return net
end

"""
    cgmes_adapter_options(cfg::SparlectraConfig; run_kind, name) -> CGMESAdapterOptions

The adapter options an effective run configuration implies.
"""
function cgmes_adapter_options(cfg::SparlectraConfig; run_kind::Symbol = :powerflow, name::AbstractString = "")
  cg = cfg.cgmes
  return CGMESAdapterOptions(
    start_values = cg.start_values,
    require_boundary = cg.require_boundary,
    infer_base_voltages = cg.infer_base_voltages,
    hvdc_mode = cg.hvdc_mode,
    run_kind = run_kind,
    name = String(name),
  )
end

"""
    cgmes_enrich_case!(case, net) -> SCFCase

Attach the CGMES source identity of `net` to its typed case: the full
structural-key mRID registry as `meta["cgmes_ids"]`, and the mRID of every
component whose extra record already resolved one under the `cgmes` key of
that record (adapter task stage 3c surface).
"""
function cgmes_enrich_case!(case::SCFCase, net::Net)::SCFCase
  spar = case.sparlectra
  spar === nothing && return case
  isempty(net.cgmes_ids) || (spar.meta["cgmes_ids"] = Dict{String,Any}(String(k) => String(v) for (k, v) in net.cgmes_ids))
  for (_, e) in spar.extra
    e isa AbstractDict || continue
    get(e, "source_id_kind", "") == "cgmes_mrid" || continue
    e["cgmes"] = Dict{String,Any}("mrid" => e["external_id"])
  end
  return case
end

function convert_case(::CGMESAdapter, source::AbstractString, opts::CGMESAdapterOptions)::SCFCase
  overrides = Dict{String,Any}(
    "cgmes_import" => Dict{String,Any}(
      "start_values" => String(opts.start_values),
      "require_boundary" => opts.require_boundary,
      "infer_base_voltages" => opts.infer_base_voltages,
      "hvdc_mode" => String(opts.hvdc_mode),
    ),
  )
  cfg = load_sparlectra_config(DEFAULT_SPARLECTRA_CONFIG_PATH; reload = true, overrides = overrides)
  cg = _import_cgmes(String(source), cfg; name = isempty(opts.name) ? basename(String(source)) : opts.name, run_kind = opts.run_kind)
  case = net_to_scfcase(cg.result.net; source_format = "cgmes", source_reference = basename(String(source)), include_start_state = true)
  return cgmes_enrich_case!(case, cg.result.net)
end
