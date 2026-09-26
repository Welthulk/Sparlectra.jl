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
# purpose: the one import entry point of the run path: format
#          detection, per-format construction, and the
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
- `overrides`: reserved for auto-profile recommendations as override level;
  empty while the MATPOWER context still rewrites the run configuration itself.
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

"""
    _case_file_lead(case_path; nbytes = 4096) -> String

The first non-blank bytes of a case file, at most `nbytes`, with a UTF-8
byte order mark removed. Empty for a directory, a missing file or an empty
file. This is what format decisions look at: the extension is only a hint,
the content is the evidence. The bytes are not required to be valid UTF-8
(a MATPOWER comment in a Windows code page must not break detection).
"""
function _case_file_lead(case_path::AbstractString; nbytes::Int = 4096)::String
  isfile(case_path) || return ""
  bytes = open(io -> read(io, nbytes), case_path)
  if length(bytes) >= 3 && bytes[1] == 0xEF && bytes[2] == 0xBB && bytes[3] == 0xBF
    bytes = bytes[4:end]
  end
  return String(lstrip(String(bytes)))
end

# a JSON document, hence a Sparlectra or power-grid-model case file, opens
# with the root object
_case_file_looks_like_json(lead::AbstractString)::Bool = startswith(lead, "{")

"""
    _case_format_conflict(case_path, requested) -> Union{Nothing,String}

The reason an explicitly requested case format cannot apply to `case_path`,
or `nothing` when file and format agree. The verdict comes from the CONTENT,
not from the extension: a requested `:scf` needs a file that opens with a
JSON object, a requested `:matpower` one that does not (no MATPOWER reader
takes JSON). `:cgmes` and `:dtf_for001` are left alone on purpose, they
cover directories, ZIPs and the ambiguous `.DAT` files whose whole point is
the explicit choice. A missing file gets no verdict here; the reader names
it.

Found 2026-09-24: a per-case setting `case_format = scf` saved for
`sp_case118.m` sent the MATPOWER file into the SCF reader, which reported
`SCF JSON parse error at byte 1: invalid integer ""` (byte 1 is the `f` of
`function mpc`), a message that blames a JSON file nobody was reading. The
message here names the file, what it actually starts with, and the way out.
"""
function _case_format_conflict(case_path::AbstractString, requested::Symbol)::Union{Nothing,String}
  requested in (:scf, :matpower) || return nothing
  name = repr(basename(case_path))
  if isdir(case_path)
    requested === :scf || return nothing
    return string("The case input format is set to scf (Sparlectra Case Format, a JSON file), but ", name, " is a directory. Set the case input format to auto, or choose a .scf.json file.")
  end
  isfile(case_path) || return nothing
  lead = _case_file_lead(case_path)
  if requested === :scf && !_case_file_looks_like_json(lead)
    what = isempty(lead) ? "is empty" : string("starts with ", repr(first(lead, 24)))
    return string("The case input format is set to scf (Sparlectra Case Format, a JSON file), but ", name, " ", what,
      " and is not a JSON object. Set the case input format to auto so the file decides, or choose a .scf.json file.")
  elseif requested === :matpower && _case_file_looks_like_json(lead)
    return string("The case input format is set to matpower, but ", name, " starts with a JSON object (a Sparlectra case file). ",
      "Set the case input format to auto so the file decides, or choose a MATPOWER .m file.")
  end
  return nothing
end

function _detect_case_format(case_path::AbstractString; requested::Symbol = :auto)::Symbol
  if requested !== :auto
    # an explicit format still has to fit the file, or the wrong reader
    # reports a parse error about a file of the other format
    conflict = _case_format_conflict(case_path, requested)
    conflict === nothing || throw(ArgumentError(conflict))
    return requested
  end
  # a PowSyBl bundle directory or IIDM file comes first: the CGMES probe
  # claims every directory and every .xml
  detect(PowsyblAdapter, case_path) && return :powsybl
  _looks_like_cgmes(case_path) && return :cgmes
  ext = lowercase(splitext(case_path)[2])
  # Sparlectra Case Format (#342): self-describing JSON, recognized by its
  # extension (the canonical name is <case>.scf.json) ...
  ext == ".json" && return :scf
  # ... or by its content: a file that opens with a JSON object is a case
  # file whatever it is called (the extension is only a hint)
  _case_file_looks_like_json(_case_file_lead(case_path)) && return :scf
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
net-parameter stamping lives with the import dispatch:
once per importer). The four service call sites used
to unpack the same values by hand, and a new option had to be added in four
places or it silently kept its default at three of them (the accompanying
stamp was found in the wrong function once already).
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
  # services had before the shared mapping.
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
  # power_flow.flatstart is the one start switch the Settings page offers;
  # under auto it asks for the flat start on a CGMES run as well, an explicit
  # sv or flat here still wins
  effective_start_values = if cgmes_cfg.start_values === :auto
    sv_bus_count > 0 && !cfg.powerflow.start_mode.flatstart ? :sv : :flat
  else
    cgmes_cfg.start_values
  end
  run_powerflow, start_overridden = _cgmes_start_values_powerflow(cfg.powerflow, effective_start_values)
  run_config = _copy_sparlectra_with_powerflow(cfg, run_powerflow)
  start_decision = string(
    "CGMES start values: ", effective_start_values,
    cgmes_cfg.start_values === :auto ? string(" (auto: ", sv_bus_count > 0 ? (cfg.powerflow.start_mode.flatstart ? "power_flow.flatstart asks for the flat start" : "delivery carries SvVoltage for $(sv_bus_count) bus(es)") : "no SvVoltage in this delivery", ")") : "",
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
  # the flat start overrides the other start machines for this run, for
  # every format; the caller reports the list in its run log
  general_config, flatstart_forced = _flatstart_forced_off_config(general_config)
  provenance["flatstart_forced_off"] = flatstart_forced
  if fmt === :matpower || fmt === :scf
    ctx = _import_sparlectra_context(String(path), nothing, general_config; performance_profile = performance_profile)
    provenance["auto_profile_result"] = ctx.auto_profile_result
    provenance["projected_start_applied"] = ctx.projected_start_applied
    studies = fmt === :scf ? scf_case_studies(String(path)) : _EMPTY_CASE_STUDIES
    # the auto profile's decisions ride on the imported case as dotted
    # overrides (the context still applies them to its effective config)
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
    # the CGMES importer built this network natively;
    # it IS the result. The old capture into an SCFCase plus a second
    # build_net doubled the construction (0.9 to 1.8 s and about 190 MB on
    # RealGrid) for nothing anybody asked: converting to SCF is an explicit
    # action, never a step inside an import.
    # stamped exactly once PER IMPORTER, and this is the CGMES
    # importer's once (with the run
    # config, exactly what the discarded second build used to apply).
    _apply_config_net_parameters!(cg.result.net, cg.run_config)
    return ImportedCase(cg.result.net, cg.run_config, fmt, provenance, _EMPTY_CASE_STUDIES, Dict{String,Any}())
  elseif fmt === :powsybl
    phase_callback("reading_powsybl_bundle")
    net, report, tables = _import_powsybl(String(path), powsybl_adapter_options(general_config); name = name)
    # stamped exactly once PER IMPORTER, and this is the PowSyBl
    # importer's once
    _apply_config_net_parameters!(net, general_config)
    provenance["powsybl_report"] = report
    provenance["powsybl_manifest"] = tables.manifest
    return ImportedCase(net, general_config, fmt, provenance, _EMPTY_CASE_STUDIES, Dict{String,Any}())
  elseif fmt === :dtf_for001
    _reject_dtf_dcline_like_content!(String(path))
    dtf_case = DTFImporter.read_dtf(String(path))
    # the DTF importer builds its network directly; the
    # model options come from the same config keys the old conversion read
    net = DTFImporter.build_net(dtf_case; bus_shunt_model = general_config.model.bus_shunt_model, tap_changer_model = general_config.model.tap_changer_model)
    # stamped exactly once PER IMPORTER, and this is the DTF
    # importer's once
    _apply_config_net_parameters!(net, general_config)
    provenance["dtf_case"] = dtf_case
    return ImportedCase(net, general_config, fmt, provenance, _EMPTY_CASE_STUDIES, Dict{String,Any}())
  end
  throw(ArgumentError("import_case: unsupported case format $(fmt) for $(path)."))
end

const DTF_FOR001_UNSUPPORTED_DCLINE_MESSAGE = "DC lines are currently not supported by the native DTF/MATPOWER power-flow path."

# Resolve the delivery parts of one CGMES run: the case itself, any extra
# parts from cgmes_import.path (typically the boundary set), and — for a bare
# base case — a boundary delivery sitting next to it in the same directory
# (ENTSO-E names them "...Boundary..." or "..._BD_..."). Returns the path list
# plus whether the boundary was autodetected. Shared by the power-flow CGMES
# branch and the short-circuit service run.
function _cgmes_delivery_paths(case_path::AbstractString, cgmes_cfg)::Tuple{Vector{String},Bool}
  extra = filter(!isempty, strip.(split(cgmes_cfg.path, ';')))
  paths = String[String(case_path)]
  for p in extra
    p == case_path || push!(paths, String(p))
  end
  length(paths) > 1 && return paths, false
  neighbours = try
    readdir(dirname(abspath(case_path)))
  catch
    String[]
  end
  for n in neighbours
    occursin(r"(?i)boundary|_bd_|_bd\.", n) || continue
    cand = joinpath(dirname(abspath(case_path)), n)
    cand == abspath(case_path) && continue
    (isdir(cand) || endswith(lowercase(n), ".zip")) || continue
    push!(paths, cand)
    return paths, true
  end
  return paths, false
end

function _reject_dtf_dcline_like_content!(case_path::AbstractString)
  # open(...) do closes the handle even when the ArgumentError below aborts
  # the scan mid-file; eachline(path) would keep the descriptor open until GC
  # (blocks file deletion on Windows, EBUSY).
  open(case_path) do io
    for (line_no, line) in enumerate(eachline(io))
      if occursin(r"(?i)\b(HVDC|DCLINE|DC\s*LINE)\b", line)
        throw(ArgumentError("unsupported_dtf_dc_line at line $(line_no): $(DTF_FOR001_UNSUPPORTED_DCLINE_MESSAGE)"))
      end
    end
  end
  return nothing
end
