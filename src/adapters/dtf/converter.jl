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

# file: src/adapters/dtf/converter.jl
# purpose: the DTF adapter (adapter task stage 3b): convert_case turns a
#          parsed native DTFCase into the typed SCFCase. The conversion
#          core stays the proven native builder (DTFImporter.build_net),
#          the built network is captured as its typed case, and the
#          importer detail the SCF data model has no slot for (the DTF
#          branch metadata record, the typed phase-tap model, the outage
#          labels, the original bus names) travels in the namespaced block
#          through the tagged value encoding, so build_net restores the
#          network the native builder produced, field for field.

struct DTFAdapter <: FormatAdapter end

"""
    DTFAdapterOptions

The adapter-scope options of the DTF conversion (design decision D4):
parser strictness and base power, the transformer ratio convention, the
legacy voltage-level collapse, the model settings consumed while the case
is interpreted, and the outage selection (a DTF adapter option per review
point 2 of the adapter task).
"""
Base.@kwdef struct DTFAdapterOptions
  baseMVA::Float64 = 100.0
  strict::Bool = true
  legacy_voltage_level_collapse_230kv::Bool = false
  transformer_ratio_mode::Symbol = :neutral_one
  tap_changer_model::Symbol = :ideal
  bus_shunt_model::Symbol = :admittance
  outage_mode::Symbol = :none
  outage_selection::Vector{String} = String[]
end

detect(::Type{DTFAdapter}, path::AbstractString)::Bool = begin
  ext = lowercase(splitext(String(path))[2])
  ext == ".dat" || return false
  return isfile(path)
end

options_type(::DTFAdapter) = DTFAdapterOptions

"""
    import_net(::DTFAdapter, case, opts::DTFAdapterOptions) -> Net

The DTF importer of the adapter contract (task_import_direct): builds
the network DIRECTLY through `DTFImporter.build_net`; the adapter
options carry every model knob the build takes.
"""
import_net(::DTFAdapter, case, opts::DTFAdapterOptions)::Net = DTFImporter.build_net(case; bus_shunt_model = opts.bus_shunt_model, legacy_voltage_level_collapse_230kv = opts.legacy_voltage_level_collapse_230kv, transformer_ratio_mode = opts.transformer_ratio_mode, tap_changer_model = opts.tap_changer_model)

"""
    dtf_adapter_options(cfg::SparlectraConfig) -> DTFAdapterOptions

The adapter options an effective run configuration implies. The DTF
conventions themselves (`transformer_ratio_mode`, the legacy collapse)
have no configuration keys yet (adapter task stage 3b inventory, open in
issue #1 point 3), so they stay at the importer defaults here.
"""
function dtf_adapter_options(cfg::SparlectraConfig)
  return DTFAdapterOptions(tap_changer_model = cfg.model.tap_changer_model, bus_shunt_model = cfg.model.bus_shunt_model)
end

function convert_case(::DTFAdapter, dtf_case, opts::DTFAdapterOptions)::SCFCase
  net = DTFImporter.build_net(
    dtf_case;
    bus_shunt_model = opts.bus_shunt_model,
    legacy_voltage_level_collapse_230kv = opts.legacy_voltage_level_collapse_230kv,
    transformer_ratio_mode = opts.transformer_ratio_mode,
    tap_changer_model = opts.tap_changer_model,
  )
  case = net_to_scfcase(net; source_format = "dtf", include_start_state = true)
  spar = case.sparlectra
  extra = spar.extra
  # file id per net index, from the recorded build-order indices
  node_id_by_idx = Dict{Int,Int}()
  branch_id_by_idx = Dict{Int,Int}()
  for (idstr, e) in extra
    e isa AbstractDict || continue
    haskey(e, "bus_index") && (node_id_by_idx[_scf_int(e["bus_index"], "extra.bus_index")] = parse(Int, idstr))
    haskey(e, "branch_index") && (branch_id_by_idx[_scf_int(e["branch_index"], "extra.branch_index")] = parse(Int, idstr))
  end
  # original display names (for DTF the reference name itself)
  for (idx, name) in net.busOriginalNameDict
    id = get(node_id_by_idx, idx, nothing)
    id === nothing && continue
    extra[string(id)]["original_name"] = String(name)
  end
  # the DTF branch metadata record, order preserving; and the typed
  # phase-tap model of the winding where one is attached
  trafo_k = 0
  for i in eachindex(net.branchVec)
    id = get(branch_id_by_idx, i, nothing)
    is_trafo = _scf_is_transformer(net.branchVec[i])
    is_trafo && (trafo_k += 1)
    id === nothing && continue
    meta = get(net.matpower_branch_metadata, i, nothing)
    meta === nothing || (extra[string(id)]["dtf_branch"] = scf_encode_value(meta))
    if is_trafo && trafo_k <= length(net.trafos)
      pt = net.trafos[trafo_k].side1.phase_taps
      pt === nothing || (extra[string(id)]["phase_taps"] = scf_encode_value(pt))
    end
  end
  # outage labels: the run path's N-1 selection reads them from the net
  labels = String[DTFImporter.outage_label(o) for o in dtf_case.outages]
  isempty(labels) || (spar.components.for001_contingency = Dict{String,Any}[Dict{String,Any}("label" => l) for l in labels])
  return case
end
