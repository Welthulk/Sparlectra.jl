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

# file: src/import/case_import.jl
# purpose: the one import entry point of the run path (stage 0 of the
#          adapter task): format detection, per-format construction, and the
#          ImportedCase record every service consumes. Format dispatch used
#          to exist five times across the services, and each copy applied a
#          different subset of the configuration.

"""
    ImportedCase

What one case import hands to a service run:

- `net`: the constructed network, configuration parameters stamped.
- `config`: the effective run configuration AFTER the import. It can differ
  from the general configuration that went in: the MATPOWER path may apply
  an auto profile and a projected start, the CGMES path resolves
  `cgmes_import.start_values` against the delivery it actually read.
- `format`: the detected (or requested) case format.
- `provenance`: source path, per-format import records (the CGMES import
  result with its short-circuit data and SV bookkeeping, the parsed DTF
  case, the auto-profile decision), for the artifacts a service writes.
- `studies`: the study definitions a case file carries (`contingencies`,
  `short_circuit`); empty for every foreign format.
- `overrides`: reserved for auto-profile recommendations as override level
  (design decision D11); empty in stage 0, where the MATPOWER context still
  rewrites the run configuration itself.
"""
struct ImportedCase
  net::Net
  config::SparlectraConfig
  format::Symbol
  provenance::Dict{String,Any}
  studies::NamedTuple
  overrides::Dict{String,Any}

  # The imported network carries the configuration it was built with, so a
  # later `runpf!(net; ...)` without an explicit `config` solves it under that
  # configuration instead of the globally active one. One inner constructor
  # covers every format path, so no importer can forget it.
  function ImportedCase(net::Net, config::SparlectraConfig, format::Symbol, provenance::Dict{String,Any}, studies::NamedTuple, overrides::Dict{String,Any})
    net._import_config = config
    return new(net, config, format, provenance, studies, overrides)
  end
end

const _EMPTY_CASE_STUDIES = (contingencies = Dict{String,Any}(), short_circuit = Dict{String,Any}())

function _looks_like_cgmes(case_path::AbstractString)::Bool
  isdir(case_path) && return true
  ext = lowercase(splitext(case_path)[2])
  ext == ".zip" && return true
  ext == ".xml" && return true
  return false
end

function _detect_case_format(case_path::AbstractString; requested::Symbol = :auto)::Symbol
  requested !== :auto && return requested
  _looks_like_cgmes(case_path) && return :cgmes
  ext = lowercase(splitext(case_path)[2])
  # Sparlectra Case Format (#342): self-describing JSON, recognized by its
  # extension (the canonical name is <case>.scf.json)
  ext == ".json" && return :scf
  ext in (".m", ".jl") && return :matpower
  text = read(case_path, String)
  # Native FOR001 test data has explicit section cards and a DTF size card.  Do
  # not infer arbitrary .DAT files unless these FOR001 markers are present.
  has_for001_sections = occursin("##Z", text) && occursin("##L", text) && occursin("##K", text)
  has_size_card = occursin(r"(?m)^##G\s+\d+\s+\d+", text)
  if has_for001_sections && has_size_card
    return :dtf_for001
  end
  ext == ".dat" && throw(ArgumentError("Ambiguous .DAT input; set case_format = :dtf_for001 to use the experimental/internal native DTF path."))
  return :matpower
end

"""
    importCGMES(config; path, name) -> CGMESImportResult

Import a CGMES delivery with everything the configuration says about it, in
one place: the thirteen `cgmes_import` values and the bus shunt model (the
net-parameter stamping lives with the import dispatch since
task_import_direct: once per importer). The four service call sites used
to unpack the same values by hand, and a new option had to be added in four
places or it silently kept its default at three of them (private issue #2;
the review found the accompanying stamp in the wrong function once already).
The keyword form of `importCGMES` stays for programmatic use.
"""
function CGMESImporter.importCGMES(config::SparlectraConfig; path, name::AbstractString, hvdc_mode::Symbol = config.cgmes.hvdc_mode)
  cgmes_cfg = config.cgmes
  result = importCGMES(
    path = path,
    baseMVA = cgmes_cfg.base_mva,
    bus_shunt_model = config.model.bus_shunt_model,
    require_boundary = cgmes_cfg.require_boundary,
    tap_control = cgmes_cfg.tap_control,
    machine_control = cgmes_cfg.machine_control,
    ignore_connected = cgmes_cfg.ignore_connected,
    vset_min_pu = cgmes_cfg.vset_min_pu,
    vset_max_pu = cgmes_cfg.vset_max_pu,
    multi_slack = cgmes_cfg.multi_slack,
    strict_placeholder_guards = cgmes_cfg.placeholder_guards === :strict,
    infer_base_voltages = cgmes_cfg.infer_base_voltages,
    hvdc_mode = hvdc_mode,
    name = String(name),
  )
  return result
end

"""
    _import_cgmes(path, cfg; name, phase_callback) -> NamedTuple

The one CGMES import of the run path: delivery-path resolution, the
configured import, and the `cgmes_import.start_values` decision, which
belongs to the import because `auto` can only resolve against the delivery
that was actually read (does it carry a usable SvVoltage state or not).
Returns the import result, the run configuration with the start decision
applied, and the records the service artifacts need. Import errors
propagate; the caller owns the failure reporting (the power-flow service
writes a diagnostic cgmes.log from what can still be read).
"""
function _import_cgmes(path::AbstractString, cfg::SparlectraConfig; name::AbstractString, run_kind::Symbol = :powerflow, phase_callback = phase -> nothing)
  cgmes_cfg = cfg.cgmes
  paths, boundary_autodetected = _cgmes_delivery_paths(path, cgmes_cfg)
  boundary_autodetected && phase_callback("cgmes_boundary_autodetected")
  # Only the power-flow run models HVDC per the configured mode; every other
  # run kind keeps the plain injection form, exactly the behavior the
  # services had before the shared mapping (review point 0 fix-up).
  hvdc_mode = run_kind === :powerflow ? cgmes_cfg.hvdc_mode : :injections
  hvdc_mode === cgmes_cfg.hvdc_mode || @info "CGMES hvdc_mode $(hvdc_mode) for run kind $(run_kind) (configured $(cgmes_cfg.hvdc_mode) applies to power-flow runs only)"
  result = importCGMES(cfg; path = length(paths) == 1 ? paths[1] : paths, name = name, hvdc_mode = hvdc_mode)
  # cgmes_import.start_values wins over power_flow(.start_mode).flatstart on
  # CGMES runs: :flat = synthetic flat start (default), :sv = the imported
  # SvVoltage state with every competing start-value machine forced off.
  # `auto` resolves only now: whether the delivery actually carries a usable
  # SvVoltage state is an import result, not a configuration property. A
  # delivery is built around its own operating point, so starting there is
  # the honest default; a delivery without SV keeps the flat start.
  sv_bus_count = length(result.net.nodeVec) - length(result.no_sv_buses)
  effective_start_values = if cgmes_cfg.start_values === :auto
    sv_bus_count > 0 ? :sv : :flat
  else
    cgmes_cfg.start_values
  end
  run_powerflow, start_overridden = _cgmes_start_values_powerflow(cfg.powerflow, effective_start_values)
  run_config = _copy_sparlectra_with_powerflow(cfg, run_powerflow)
  start_decision = string(
    "CGMES start values: ", effective_start_values,
    cgmes_cfg.start_values === :auto ? string(" (auto: ", sv_bus_count > 0 ? "delivery carries SvVoltage for $(sv_bus_count) bus(es)" : "no SvVoltage in this delivery", ")") : "",
    effective_start_values === :sv ? " (imported SvVoltage state; start-value machines forced off: start_projection, dc_seed_unconditional, start_current_iteration, apslf_start)" : " (synthetic flat start)",
    isempty(start_overridden) ? "" : string(", overrides: ", join(start_overridden, ", ")),
  )
  return (result = result, run_config = run_config, paths = paths, boundary_autodetected = boundary_autodetected, sv_bus_count = sv_bus_count, effective_start_values = effective_start_values, start_decision = start_decision, start_overridden = start_overridden)
end

"""
    import_case(path, general_config; requested_format, name, performance_profile, phase_callback) -> ImportedCase

Import a case of any format for a service run. Detection follows
`_detect_case_format` unless `requested_format` names the format
explicitly. Construction goes through the same per-format paths the
power-flow service used, so every service gets the identical network for
the identical file: the shared context for MATPOWER and the case format
(auto profile, net cache, projected start, configuration stamping), the
configured CGMES import with the start-values decision, and the DTF build
with its DC-line rejection.

Format policy stays with the caller: a service that does not accept a
format checks `ImportedCase.format` and words its own refusal, because the
wording of those refusals is part of each service's contract.
"""
function import_case(path::AbstractString, general_config::SparlectraConfig; requested_format::Symbol = :auto, run_kind::Symbol = :powerflow, name::AbstractString = basename(path), performance_profile = nothing, phase_callback = phase -> nothing)
  fmt = _detect_case_format(String(path); requested = requested_format)
  provenance = Dict{String,Any}("source_path" => String(path))
  if fmt === :matpower || fmt === :scf
    ctx = _import_sparlectra_context(String(path), nothing, general_config; performance_profile = performance_profile)
    provenance["auto_profile_result"] = ctx.auto_profile_result
    provenance["projected_start_applied"] = ctx.projected_start_applied
    studies = fmt === :scf ? scf_case_studies(String(path)) : _EMPTY_CASE_STUDIES
    # D11: the auto profile's decisions ride on the imported case as dotted
    # overrides (the context still applies them to its effective config;
    # the full merge as an own precedence level is review-point-3 material)
    return ImportedCase(ctx.net, ctx.config, fmt, provenance, studies, ctx.auto_profile_overrides)
  elseif fmt === :cgmes
    phase_callback("reading_cgmes_delivery")
    cg = _import_cgmes(String(path), general_config; name = name, run_kind = run_kind, phase_callback = phase_callback)
    provenance["cgmes_result"] = cg.result
    provenance["cgmes_paths"] = cg.paths
    provenance["cgmes_boundary_autodetected"] = cg.boundary_autodetected
    provenance["cgmes_sv_bus_count"] = cg.sv_bus_count
    provenance["cgmes_effective_start_values"] = cg.effective_start_values
    provenance["cgmes_start_decision"] = cg.start_decision
    provenance["cgmes_start_overridden"] = cg.start_overridden
    # task_import_direct: the CGMES importer built this network natively;
    # it IS the result. The old capture into an SCFCase plus a second
    # build_net doubled the construction (0.9 to 1.8 s and about 190 MB on
    # RealGrid) for nothing anybody asked: converting to SCF is an explicit
    # action, never a step inside an import.
    # D12, corrected by task_import_direct: stamped exactly once PER
    # IMPORTER, and this is the CGMES importer's once (with the run
    # config, exactly what the discarded second build used to apply).
    _apply_config_net_parameters!(cg.result.net, cg.run_config)
    return ImportedCase(cg.result.net, cg.run_config, fmt, provenance, _EMPTY_CASE_STUDIES, Dict{String,Any}())
  elseif fmt === :dtf_for001
    _reject_dtf_dcline_like_content!(String(path))
    dtf_case = DTFImporter.read_dtf(String(path))
    # task_import_direct: the DTF importer builds its network directly; the
    # model options come from the same config keys the old conversion read
    net = DTFImporter.build_net(dtf_case; bus_shunt_model = general_config.model.bus_shunt_model, tap_changer_model = general_config.model.tap_changer_model)
    # D12, corrected by task_import_direct: stamped exactly once PER
    # IMPORTER, and this is the DTF importer's once
    _apply_config_net_parameters!(net, general_config)
    provenance["dtf_case"] = dtf_case
    return ImportedCase(net, general_config, fmt, provenance, _EMPTY_CASE_STUDIES, Dict{String,Any}())
  end
  throw(ArgumentError("import_case: unsupported case format $(fmt) for $(path)."))
end
