# Copyright 2023–2026 Udo Schmitz
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

# file: src/webui/forms.jl
# purpose: Web UI input resolution: application root and case directory
#          lookup, case/config file selectors, upload classification, and
#          form value parsing
const _WEBUI_PACKAGE_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

"""
    _webui_application_root([start_dir]) -> String

Resolve the Sparlectra application root used by the local Web UI. The lookup
supports starting Julia from the repository itself, from a parent directory
containing `Sparlectra` or `Sparlectra.jl`, and from the installed package root.
"""
function _webui_application_root(start_dir::AbstractString = pwd())::String
  start_root = abspath(start_dir)
  candidates = unique((
    start_root,
    joinpath(start_root, "Sparlectra"),
    joinpath(start_root, "Sparlectra.jl"),
    _WEBUI_PACKAGE_ROOT,
  ))
  for candidate in candidates
    isdir(joinpath(candidate, "data", "mpower")) && isdir(joinpath(candidate, "examples")) && return normpath(candidate)
  end
  return _WEBUI_PACKAGE_ROOT
end

"""
    _webui_casefile_options(application_root) -> Vector{String}

Return sorted user-selectable case filenames from the Web UI application's
`data/mpower` directory. The selector stays conservative: MATPOWER `.m` files
and copied internal DTF `.DAT` candidates are shown, while generated
Julia cache artifacts, precompile workload cases, result artifacts, and sidecar profiles
stay hidden. Missing or empty directories produce an empty list.
"""
function _webui_casefile_options(application_root::AbstractString)::Vector{String}
  return _webui_casefile_options_in_directory(_webui_case_directory(; application_root))
end

"""
    _webui_bundled_scf_options(application_root) -> Vector{String}

The shipped SCF cases under `data/scf` (`sp_*` only: the demo cases and
the PST format fixture, task_demo_cases_v0100). They are offered in the
case chooser regardless of the cache directory and staged into the cache
on first use ([`_webui_stage_bundled_case!`](@ref)).
"""
function _webui_bundled_scf_options(application_root::AbstractString)::Vector{String}
  dir = joinpath(application_root, "data", "scf")
  isdir(dir) || return String[]
  return sort!([name for name in readdir(dir) if startswith(name, "sp_") && endswith(lowercase(name), ".scf.json")]; by = lowercase)
end

"""
    _webui_stage_bundled_case!(application_root, case_directory, requested) -> Union{Nothing,String}

Resolve a bare requested case NAME against the package's bundled sources
(`data/mpower`, then `data/scf`) and stage it into `case_directory`.
Returns the staged (or already cached) absolute path, or `nothing` when no
bundled source carries the name. An SCF demo case is staged WITH its
sidecars (`<stem>.config.yaml` pins the machine-neutral resolution, the
measurement CSVs belong to the case); existing cache files are never
overwritten, so user-saved settings survive. With no writable
`case_directory` the bundled path itself is returned.
"""
function _webui_stage_bundled_case!(application_root::AbstractString, case_directory::Union{Nothing,AbstractString}, requested::AbstractString)::Union{Nothing,String}
  name = String(requested)
  (isabspath(name) || occursin('/', name) || occursin('\\', name)) && return nothing
  for source_dir in (joinpath(application_root, "data", "mpower"), joinpath(application_root, "data", "scf"))
    bundled = joinpath(source_dir, name)
    isfile(bundled) || continue
    case_directory === nothing && return bundled
    mkpath(case_directory)
    cached = joinpath(String(case_directory), name)
    if !isfile(cached)
      # the companions travel with the case, enumerated by the single
      # definition (stage 4B): config pin, measurement CSVs, weights
      for src in vcat([bundled], case_companion_files(bundled))
        dst = joinpath(String(case_directory), basename(src))
        isfile(dst) || cp(src, dst)
      end
    end
    return cached
  end
  return nothing
end

"""Return the effective Web UI case directory used for selection and imports."""
function _webui_case_directory(; case_directory::Union{Nothing,AbstractString} = nothing, application_root::AbstractString = _webui_application_root(), output_root::AbstractString = default_webui_output_root())::String
  case_directory === nothing || return normpath(String(case_directory))
  package_case_directory = joinpath(application_root, "data", "mpower")
  try
    mkpath(package_case_directory)
    test_path = tempname(package_case_directory)
    open(test_path, "w") do io
      write(io, "")
    end
    rm(test_path; force = true)
    return normpath(package_case_directory)
  catch
    return normpath(default_webui_case_cache_dir(output_root))
  end
end

# .zip covers CGMES deliveries (single ZIP, base case plus boundary, or
# ZIP-in-ZIP) — the importer's container layer opens them in memory.
# .csv covers measurement sets (SE phase 5): the role is decided by a content
# sniff on the version comment, not the extension (see _webui_is_measurement_csv).
# .json is a Sparlectra Case Format case (#342); the import path validates the
# content before it is stored, so a foreign .json is rejected by name
_webui_supported_upload_case_extension(name::AbstractString)::Bool = lowercase(splitext(basename(String(name)))[2]) in (".m", ".dat", ".zip", ".csv", ".json")

"""
    _webui_scf_upload_reason(bytes) -> Union{Nothing,String}

Why an uploaded `.json` is not a usable case file, or `nothing` when it is
one. Parsing and the schema/reference validation run on the bytes BEFORE
anything is stored, so a foreign JSON never lands in the case directory and
the user gets the reader's own message.
"""
function _webui_scf_upload_reason(bytes::Vector{UInt8})::Union{Nothing,String}
  try
    root = scf_json_parse(String(copy(bytes)))
    scf_validate_dataset(root)
    return nothing
  catch err
    return first(split(sprint(showerror, err), '\n'))
  end
end

function _webui_classify_dat_content(path::AbstractString)::Symbol
  lowercase(splitext(path)[2]) == ".dat" || return :not_dat
  text = try
    read(path, String)
  catch
    return :unknown_dat
  end
  try
    case = DTFImporter.read_dtf(path; strict = false)
    if case.size.NGES > 0 && length(case.buses) == case.size.NGES && length(case.branches) == case.size.LGES
      return isempty(case.outages) ? :dtf_network_case : :dtf_network_case_with_outages
    end
  catch
  end
  has_outages = occursin(r"(?im)^\s*AUSFALL\s*$", text) && occursin(r"(?im)^\s*ENDE\s*$", text)
  if has_outages
    return :dtf_outage_file
  elseif _is_for002_reference_dat(path) || occursin(r"(?i)(for002|reference|vergleich|result)", text)
    return :dtf_outage_or_reference
  end
  return :unknown_dat
end

_webui_dat_role_label(role::Symbol)::String = replace(String(role), '_' => ' ')

"""True when a CGMES delivery contains boundary profiles only (no equipment)."""
function _webui_is_cgmes_boundary_only(path::AbstractString)::Bool
  try
    profiles = union((f.profiles for f in CGMESImporter.summarizeCGMES(path = path).files)...)
    return isempty(intersect(profiles, (:EQ, :EQ_OP, :EQ_SC, :TP, :SSH))) && !isempty(intersect(profiles, (:EQ_BD, :EQ_BD_OP, :TP_BD)))
  catch
    return false
  end
end

"""
Pre-analysis of an uploaded CGMES ZIP, shown in the import summary so the
user sees immediately what the delivery contains and whether it is ready to
compute: version, profiles, element counts, and the boundary situation.
"""
function _webui_cgmes_upload_role(path::AbstractString, directory::AbstractString)::String
  summary = try
    CGMESImporter.summarizeCGMES(path = path)
  catch err
    return "⚠️ unreadable CGMES delivery (" * first(sprint(showerror, err), 80) * ")"
  end
  profiles = union((f.profiles for f in summary.files)...)
  counts = Dict(summary.class_histogram)
  cnt(cls) = get(counts, cls, 0)
  version = isempty(summary.version) ? "CGMES" : "CGMES " * summary.version
  proflist = join(sort(String.(collect(profiles))), ",")

  # a boundary-only delivery has no equipment of its own
  if isempty(intersect(profiles, (:EQ, :EQ_OP, :EQ_SC, :TP, :SSH))) && !isempty(intersect(profiles, (:EQ_BD, :EQ_BD_OP, :TP_BD)))
    return "🔗 $(version) boundary set · $(proflist) · $(cnt(:TopologicalNode)) boundary nodes — pair it with a base case"
  end

  parts = String[]
  push!(parts, "$(cnt(:TopologicalNode)) nodes")
  cnt(:ACLineSegment) > 0 && push!(parts, "$(cnt(:ACLineSegment)) lines")
  cnt(:PowerTransformer) > 0 && push!(parts, "$(cnt(:PowerTransformer)) transformers")
  total_loads = cnt(:EnergyConsumer) + cnt(:ConformLoad) + cnt(:NonConformLoad)
  total_loads > 0 && push!(parts, "$(total_loads) loads")
  cnt(:SynchronousMachine) > 0 && push!(parts, "$(cnt(:SynchronousMachine)) machines")
  cnt(:LinearShuntCompensator) > 0 && push!(parts, "$(cnt(:LinearShuntCompensator)) shunts")
  cnt(:SvVoltage) > 0 && push!(parts, "SV reference present")
  content = join(parts, ", ")

  missing_profiles = [p for p in (:EQ, :SSH, :TP) if !(p in profiles)]
  isempty(missing_profiles) || return "❌ $(version) · $(proflist) · $(content) — incomplete: missing profile(s) $(join(String.(missing_profiles), ", "))$(_webui_upload_analysis_note(path))"

  if !summary.boundary_missing_hint
    return "✅ $(version) · $(proflist) · $(content) · self-contained, ready to compute"
  end
  # a boundary delivery next to it completes this one
  neighbours = try
    filter(n -> occursin(r"(?i)boundary|_bd_|_bd\.", n) && joinpath(directory, n) != path, readdir(directory))
  catch
    String[]
  end
  for n in neighbours
    try
      s2 = CGMESImporter.summarizeCGMES(path = [path, joinpath(directory, n)])
      s2.boundary_missing_hint || return "✅ $(version) · $(proflist) · $(content) · boundary found in $(n), ready to compute"
    catch
    end
  end
  # Last resort before bothering the user: the ENTSO-E test configurations are
  # published as base case plus a separate boundary package. If the extracted
  # local test-set cache holds a boundary that completes this delivery, copy
  # it next to the case instead of asking for a second upload.
  fetched = _webui_try_supply_cgmes_boundary(path, directory)
  isempty(fetched) || return "✅ $(version) · $(proflist) · $(content) · boundary supplied automatically ($(fetched)), ready to compute"
  return "⚠️ $(version) · $(proflist) · $(content) · $(summary.unresolved_count) unresolved references — boundary set missing: upload it as well, or type cgmes:&lt;alias&gt;$(_webui_upload_analysis_note(path))"
end

"""
Upload-time import analysis for a CGMES delivery that stayed incomplete
(missing profiles or missing boundary): writes the full
[`importFailureAnalysis`](@ref Sparlectra.CGMESImporter.importFailureAnalysis)
report next to the case as `<case>.import_analysis.txt` and returns a short
message suffix naming the missing declared dependencies. Best-effort — an
unreadable delivery yields an empty suffix and no file.
"""
function _webui_upload_analysis_note(path::AbstractString)::String
  try
    store = CGMESImporter.loadCGMES(path)
    stats = CGMESImporter.importabilityStats(store)
    report_file = string(splitext(path)[1], ".import_analysis.txt")
    open(report_file, "w") do io
      print(io, CGMESImporter.importFailureAnalysis(store))
    end
    if !isempty(stats.missing_dependencies)
      shown = join(Iterators.take(stats.missing_dependencies, 3), ", ")
      length(stats.missing_dependencies) > 3 && (shown *= ", …")
      return " · missing declared dependencies: $(shown) — full analysis in $(basename(report_file))"
    end
    return " — full analysis in $(basename(report_file))"
  catch
    # keep the role summary useful even when the analysis cannot be built
    return ""
  end
end

"""
Look for a boundary delivery in the local ENTSO-E test-set cache that resolves
`path`, copy it next to the case, and return its filename (empty when none
fits). Only uses what is already extracted locally — no download is triggered
from an upload.
"""
function _webui_try_supply_cgmes_boundary(path::AbstractString, directory::AbstractString)::String
  extracted = joinpath(CGMESImporter.cgmesTestSetCacheDir(), "extracted")
  isdir(extracted) || return ""
  candidates = String[]
  for (root, dirs, files) in walkdir(extracted)
    for d in dirs
      occursin(r"(?i)boundary|_bd_v|_bd$", d) && push!(candidates, joinpath(root, d))
    end
  end
  for cand in candidates
    ok = try
      !CGMESImporter.summarizeCGMES(path = [path, cand]).boundary_missing_hint
    catch
      false
    end
    ok || continue
    name = basename(cand) * ".zip"
    dest = joinpath(directory, name)
    ispath(dest) && return name
    try
      open(dest, "w") do io
        CGMESImporter.ZipArchives.ZipWriter(io) do w
          for f in sort(readdir(cand))
            endswith(lowercase(f), ".xml") || continue
            CGMESImporter.ZipArchives.zip_newfile(w, f)
            write(w, read(joinpath(cand, f)))
          end
        end
      end
      return name
    catch
      isfile(dest) && rm(dest; force = true)
      return ""
    end
  end
  return ""
end


_webui_is_runnable_dat_role(role::Symbol)::Bool = role in (:dtf_network_case, :dtf_network_case_with_outages)

"""
    _webui_is_user_selectable_case(name) -> Bool

Return whether `name` should be shown in the normal Web UI case selector.
The bundled precompile workloads are `warmup_*.jl` files and stay out of the
selector; a MATPOWER `.m` case may carry the `warmup_` prefix and stays
selectable. Generated Julia cache
artifacts are also hidden from the selector; users can still enter an
explicit path in the manual case field when they intentionally want to run
such a file.
"""
# Memoization for the per-file content checks of the case selector: the ZIP
# boundary-set detection and the DAT role classification read file content,
# and the form re-scans the whole case directory on every render (measured
# ~10 s per render over a 59-file cache with large CGMES ZIPs). Results are
# keyed by (absolute path, mtime, size); a replaced or edited file is
# re-classified, deleted files age out via the key check on the next hit.
const _WEBUI_CASE_SCAN_LOCK = ReentrantLock()
const _WEBUI_CASE_SCAN_CACHE = Dict{String,Tuple{Float64,Int64,Any}}()

function _webui_file_scan_memo(f::Function, path::AbstractString)
  st = try
    stat(path)
  catch
    return f(path)
  end
  key = string(abspath(path), "|", nameof(f))
  mt = Float64(Base.Filesystem.mtime(st))
  sz = Int64(st.size)
  return lock(_WEBUI_CASE_SCAN_LOCK) do
    hit = get(_WEBUI_CASE_SCAN_CACHE, key, nothing)
    if hit !== nothing && hit[1] == mt && hit[2] == sz
      hit[3]
    else
      value = f(path)
      _WEBUI_CASE_SCAN_CACHE[key] = (mt, sz, value)
      value
    end
  end
end

_webui_classify_dat_content_cached(path::AbstractString)::Symbol = _webui_file_scan_memo(_webui_classify_dat_content, path)
_webui_is_cgmes_boundary_only_cached(path::AbstractString)::Bool = _webui_file_scan_memo(_webui_is_cgmes_boundary_only, path)

"""
    _webui_is_case_json(path) -> Bool

Whether a `.json` file in the case directory is a case file. Decided by
CONTENT, exactly like the DAT role and the CGMES boundary check: the
document is parsed and validated. A run artifact that happens to end in
`.json` therefore stays out of the selector, while both export products
(`.scf.json` and the plain `.pgm.json`) and any validated upload appear.
"""
function _webui_is_case_json(path::AbstractString)::Bool
  try
    scf_validate_dataset(scf_json_parse(read(String(path), String)))
    return true
  catch
    return false
  end
end

_webui_is_case_json_cached(path::AbstractString)::Bool = _webui_file_scan_memo(_webui_is_case_json, path)

function _webui_is_user_selectable_case(name::AbstractString)::Bool
  lowered_name = lowercase(basename(name))
  _, extension = splitext(lowered_name)
  endswith(lowered_name, ".sparlectra-webui.yaml") && return false
  # per-case configuration files travel next to their case and are not cases
  endswith(lowered_name, ".config.yaml") && return false
  endswith(lowered_name, ".contingency-weights.csv") && return false
  # the bundled precompile workloads are .jl files (warmup_case3/case118);
  # a regular MATPOWER case may legitimately carry the prefix
  # (warmup_casePST.m) and must stay selectable
  startswith(lowered_name, "warmup_") && extension == ".jl" && return false
  extension == ".jl" && return false
  # A .json is a case file when it CONTAINS one: the export writes
  # `.scf.json` (full) and `.pgm.json` (plain power-grid-model dataset), an
  # upload is validated before it is stored, and a run artifact that happens
  # to end in .json parses as neither. Hiding everything but `.scf.json` made
  # an exported plain-PGM file unfindable, which is how this was found.
  if extension == ".json"
    isfile(name) || return endswith(lowered_name, ".scf.json") || endswith(lowered_name, ".pgm.json")
    return _webui_is_case_json_cached(name)
  end
  extension in (".m", ".dat", ".zip") || return false
  if extension == ".zip"
    # A boundary set is a companion delivery, not a runnable case — hide it
    # the same way FOR002 reference files are hidden.
    isfile(name) || return true
    return !_webui_is_cgmes_boundary_only_cached(name)
  end
  if extension == ".dat"
    isfile(name) || return !_is_for002_reference_dat(name)
    return _webui_is_runnable_dat_role(_webui_classify_dat_content_cached(name))
  end
  return true
end

function _webui_casefile_options_in_directory(directory::AbstractString)::Vector{String}
  isdir(directory) || return String[]
  files = filter(readdir(directory)) do name
    path = joinpath(directory, name)
    return isfile(path) && _webui_is_user_selectable_case(path)
  end
  return sort!(files; by = lowercase)
end

function _webui_for002_reference_options_in_directory(directory::AbstractString)::Vector{String}
  isdir(directory) || return String[]
  files = filter(readdir(directory)) do name
    path = joinpath(directory, name)
    return isfile(path) && lowercase(splitext(name)[2]) == ".dat" && _webui_classify_dat_content_cached(path) === :dtf_outage_or_reference
  end
  return sort!(files; by = lowercase)
end

function _webui_is_config_file(path::AbstractString)::Bool
  name = lowercase(basename(path))
  return endswith(name, ".yaml") || endswith(name, ".yml") || endswith(name, ".yaml.example") || endswith(name, ".yml.example")
end

"""
    _webui_config_file_options(application_root) -> Vector{String}

Return sorted YAML configuration files and YAML example templates from the Web
UI application's `examples` directory. A local `configuration.yaml` or
`configuration.yaml.example` is ordered first when present.
"""
function _webui_config_file_options(application_root::AbstractString)::Vector{String}
  directory = joinpath(application_root, "examples")
  isdir(directory) || return String[]
  files = filter(path -> isfile(path) && _webui_is_config_file(path), readdir(directory; join = true))
  priority(path) = lowercase(basename(path)) in ("configuration.yaml", "configuration.yaml.example") ? 0 : 1
  return sort!(normpath.(files); by = path -> (priority(path), lowercase(basename(path))))
end

function _webui_form_string(value)::String
  value === nothing && return ""
  value === missing && return ""
  value isa AbstractString && return String(value)
  value isa Symbol && return String(value)
  value isa Bool && return value ? "true" : "false"
  value isa Integer && return string(value)
  value isa AbstractFloat && return isfinite(value) ? string(value) : throw(ArgumentError("Web UI form value must be finite."))
  throw(ArgumentError("Unsupported Web UI form value type $(typeof(value))."))
end

function _webui_form_bool(value)::Bool
  value isa AbstractVector && return any(_webui_form_bool, value)
  value isa Bool && return value
  value === nothing && return false
  value === missing && return false
  value isa AbstractString && return lowercase(strip(value)) in ("1", "true", "yes", "on")
  value isa Integer && return value != 0
  throw(ArgumentError("Unsupported Web UI checkbox value type $(typeof(value))."))
end

function _webui_form_number_string(value)::String
  value isa Bool && throw(ArgumentError("Boolean is not a numeric Web UI form value."))
  value isa Integer && return string(value)
  value isa AbstractFloat && return isfinite(value) ? string(value) : throw(ArgumentError("Web UI numeric form value must be finite."))
  value isa AbstractString && begin
    text = strip(value)
    isempty(text) && throw(ArgumentError("Web UI numeric form value must not be empty."))
    parsed = tryparse(Float64, text)
    parsed === nothing && throw(ArgumentError("Invalid numeric Web UI form value $(repr(text))."))
    isfinite(parsed) || throw(ArgumentError("Web UI numeric form value must be finite."))
    return text
  end
  throw(ArgumentError("Unsupported Web UI numeric value type $(typeof(value))."))
end

# name of the RETIRED settings sidecar; still known so a leftover file can
# be converted once and deleted (D8 of the adapter task)
function _webui_legacy_case_settings_filename(casefile::AbstractString)::String
  stem = splitext(basename(strip(String(casefile))))[1]
  isempty(stem) && throw(ArgumentError("Case-settings profile requires a case filename."))
  return string(stem, ".sparlectra-webui.yaml")
end

function _webui_normalized_case_key(casefile::AbstractString)::String
  stem = replace(lowercase(splitext(basename(strip(String(casefile))))[1]), r"[^a-z0-9]+" => "_")
  stem = strip(stem, '_')
  return isempty(stem) ? "case" : stem
end

# Shared selected-case state (stage 4A harmonization): the Case page is THE
# place to choose a case, but the choice must reach every page that starts
# work on that case (the run form's hidden casefile, the SE page's loader)
# even when the user gets there through the plain nav links, which carry no
# query. One local user, one remembered name, stored under the output root
# so it survives a Web UI restart. A page whose URL names a case explicitly
# still wins and updates the memory.
_webui_selected_case_state_path(output_root::AbstractString)::String = joinpath(String(output_root), "webui_selected_case.txt")

function _webui_remember_selected_case!(output_root::AbstractString, casefile::AbstractString)::String
  path = _webui_selected_case_state_path(output_root)
  name = strip(String(casefile))
  try
    if isempty(name)
      rm(path; force = true)
    else
      mkpath(dirname(path))
      write(path, name)
    end
  catch
    # remembering is a convenience; failing to persist must never break a render
  end
  return name
end

function _webui_recall_selected_case(output_root::AbstractString)::String
  path = _webui_selected_case_state_path(output_root)
  isfile(path) || return ""
  return try
    String(strip(read(path, String)))
  catch
    ""
  end
end

function _webui_resolve_case_profile_source(casefile::AbstractString; case_directory::Union{Nothing,AbstractString} = nothing, application_root::AbstractString = _webui_application_root())::String
  raw = strip(String(casefile))
  isempty(raw) && return ""
  if isabspath(raw)
    return normpath(raw)
  end
  ".." in splitpath(raw) && return ""
  base = case_directory === nothing ? joinpath(application_root, "data", "mpower") : String(case_directory)
  return normpath(joinpath(base, raw))
end

"""
    _webui_case_config_from_settings(settings_raw) -> (keep, dropped)

Split a run's recorded settings into the dotted configuration keys a case
file may carry (`keep`) and the keys that belong to the machine (`dropped`:
output, benchmark, runtime, webui, export scope). Request metadata entries
without a dot (casefile, performance_timing, ...) are neither.
"""
function _webui_case_config_from_settings(settings_raw::AbstractDict)
  keep = Dict{String,Any}()
  dropped = String[]
  for (k, v) in settings_raw
    key = String(k)
    occursin(".", key) || continue
    if scf_is_case_config_key(key)
      keep[key] = v
    else
      push!(dropped, key)
    end
  end
  return keep, sort!(dropped)
end

function _webui_case_settings_path(output_root::AbstractString, casefile::AbstractString; case_directory::Union{Nothing,AbstractString} = nothing)::String
  source = isabspath(strip(String(casefile))) ? normpath(String(casefile)) : _webui_resolve_case_profile_source(casefile; case_directory)
  isempty(source) && throw(ArgumentError("Unsafe MATPOWER case path for case-settings profile."))
  # saved settings live in the case configuration file next to the case
  return case_config_path(source)
end

function _webui_legacy_case_settings_path(output_root::AbstractString, casefile::AbstractString; case_directory::Union{Nothing,AbstractString} = nothing)::String
  source = isabspath(strip(String(casefile))) ? normpath(String(casefile)) : _webui_resolve_case_profile_source(casefile; case_directory)
  isempty(source) && throw(ArgumentError("Unsafe MATPOWER case path for case-settings profile."))
  return joinpath(dirname(source), _webui_legacy_case_settings_filename(source))
end

# Per-case N-1 contingency weight list (issue #331 Phase 5 follow-up): stored
# next to the case as `<stem>.contingency-weights.csv`, the exact two-column
# format `readContingencyWeightsCSV` parses. Mirrors the sidecar-YAML naming so
# it is human-editable and travels with the case (list-excluded, delete-cascaded).
function _webui_case_weights_filename(casefile::AbstractString)::String
  stem = splitext(basename(strip(String(casefile))))[1]
  isempty(stem) && throw(ArgumentError("Contingency weights file requires a case filename."))
  return string(stem, ".contingency-weights.csv")
end

"""
    case_companion_files(case_path) -> Vector{String}

The SINGLE definition of which files belong to one case (stage 4B review
earmark): the case configuration file, the legacy settings sidecar, the
per-case N-1 weights list, and the measurement CSVs including the
baddata and noisy variants (both stem conventions occur in the wild: the
bundled sets use the short stem, generated sets the splitext stem; the
in-file binding comment stays authoritative, the name is only the copy
rule). Staging, the delete cascade, and every future copy path
enumerate companions HERE; a new companion kind is added in this
function and nowhere else. Returns the companions that EXIST next to
the case, as absolute paths.
"""
function case_companion_files(case_path::AbstractString)::Vector{String}
  src = normpath(abspath(String(case_path)))
  dir = dirname(src)
  name = basename(src)
  isdir(dir) || return String[]
  short_stem = endswith(lowercase(name), ".scf.json") ? name[1:(end - 9)] : first(splitext(name))
  companions = String[]
  cfg = case_config_path(src)
  isfile(cfg) && push!(companions, cfg)
  legacy = joinpath(dir, _webui_legacy_case_settings_filename(src))
  isfile(legacy) && push!(companions, legacy)
  weights = joinpath(dir, _webui_case_weights_filename(name))
  isfile(weights) && push!(companions, weights)
  for f in readdir(dir)
    f == name && continue
    startswith(f, string(short_stem, ".")) || continue
    lowered = lowercase(f)
    (endswith(lowered, ".csv") && occursin("measurements", lowered)) || continue
    push!(companions, joinpath(dir, f))
  end
  return unique!(companions)
end

# No output_root parameter (unlike the sidecar helper it mirrors): the file
# lives next to the case, output_root is irrelevant, so it is left out.
function _webui_case_weights_path(casefile::AbstractString; case_directory::Union{Nothing,AbstractString} = nothing)::String
  source = isabspath(strip(String(casefile))) ? normpath(String(casefile)) : _webui_resolve_case_profile_source(casefile; case_directory)
  isempty(source) && throw(ArgumentError("Unsafe MATPOWER case path for contingency weights."))
  return joinpath(dirname(source), _webui_case_weights_filename(source))
end

function _webui_log_case_settings_load(output_root::AbstractString, event::AbstractString; fields...)
  try
    record_webui_operation!(output_root, event; route = "/powerflow", method = "GET", user_action = true, fields...)
  catch
  end
  return nothing
end

function _webui_normalize_case_profile_form_value(field::AbstractString, value)
  value === nothing && return nothing
  value === missing && return nothing
  value isa AbstractVector && throw(ArgumentError("Case-settings field $(field) does not support vector values for form rendering."))
  type = get(_WEBUI_CASE_PROFILE_FIELD_TYPES, String(field), String)
  allowed = get(_WEBUI_CASE_PROFILE_SELECT_VALUES, String(field), nothing)
  if allowed !== nothing && value isa Bool
    !value && "off" in allowed && return "off"
    value && "on" in allowed && return "on"
  end
  normalized = if type === Bool
    _webui_form_bool(value)
  elseif type <: Integer
    value isa Bool && throw(ArgumentError("Case-settings field $(field) must be an integer."))
    if value isa Integer
      Int(value)
    elseif value isa AbstractFloat && isinteger(value) && isfinite(value)
      Int(value)
    else
      text = value isa AbstractString ? strip(value) : _webui_form_string(value)
      parsed = tryparse(type, text)
      parsed === nothing && throw(ArgumentError("Case-settings field $(field) has invalid integer value $(repr(text))."))
      parsed
    end
  elseif type <: AbstractFloat
    text = _webui_form_number_string(value)
    parsed = tryparse(type, text)
    parsed === nothing && throw(ArgumentError("Case-settings field $(field) has invalid numeric value $(repr(text))."))
    parsed
  else
    _webui_form_string(value)
  end
  if allowed !== nothing && !(_webui_form_string(normalized) in allowed)
    throw(ArgumentError("Case-settings field $(field) has unsupported value $(repr(normalized))."))
  end
  return normalized
end


function _webui_config_field_values(config_file::AbstractString)::Dict{String,Any}
  values = Dict{String,Any}()
  path = strip(String(config_file))
  isempty(path) && return values
  isfile(path) || return values
  try
    raw = load_yaml_dict(path)
    for (config_key, field, _) in _WEBUI_FORM_CONFIG_FIELDS
      value = _dotted_config_value(raw, config_key)
      value === nothing && continue
      values[field] = _webui_normalize_case_profile_form_value(field, value)
    end
    return values
  catch
    return Dict{String,Any}()
  end
end

"""
    _webui_case_config_field_values(casefile, case_directory) -> Dict{String,Any}

The form fields a case sets through its configuration levels, low to high:
the deprecated `sparlectra.config` block of an SCF case, then the case
configuration file next to the case. Empty for a case without settings and
for anything unreadable, so the caller can ask unconditionally. The form
has to show these values: it posts a value for EVERY field it renders, and
those count as explicit overrides, so an unseeded form would silently
outrank the very settings it just loaded.
"""
function _webui_case_config_field_values(casefile::AbstractString, case_directory)::Dict{String,Any}
  values = Dict{String,Any}()
  name = strip(String(casefile))
  isempty(name) && return values
  path = isabspath(name) ? name : (case_directory === nothing ? name : joinpath(String(case_directory), name))
  (isfile(path) || isdir(path)) || return values
  merged = Dict{String,Any}()
  if lowercase(splitext(path)[2]) == ".json" && isfile(path)
    try
      merge!(merged, scf_case_config(path))
    catch
    end
  end
  try
    merge!(merged, load_case_config(path))
  catch
    return Dict{String,Any}()
  end
  for (config_key, field, _) in _WEBUI_FORM_CONFIG_FIELDS
    haskey(merged, config_key) || continue
    try
      values[field] = _webui_normalize_case_profile_form_value(field, merged[config_key])
    catch
    end
  end
  return values
end

"""
    _webui_case_form_defaults(casefile, case_directory) -> Dict{String,Any}

The `form` block of the case configuration file: per-case Web UI request
defaults (SE and generator options). Field-named, normalized, unknown
fields ignored. Empty when the case has no configuration file or no block.
"""
function _webui_case_form_defaults(casefile::AbstractString, case_directory)::Dict{String,Any}
  values = Dict{String,Any}()
  name = strip(String(casefile))
  isempty(name) && return values
  path = isabspath(name) ? name : (case_directory === nothing ? name : joinpath(String(case_directory), name))
  cfg_path = try
    case_config_path(path)
  catch
    return values
  end
  isfile(cfg_path) || return values
  block = try
    get(load_yaml_dict(cfg_path), "form", nothing)
  catch
    nothing
  end
  block isa AbstractDict || return values
  for (k, v) in block
    field = String(k)
    haskey(_WEBUI_OPTION_BY_FIELD, field) || continue
    try
      normalized = _webui_normalize_case_profile_form_value(field, v)
      normalized === nothing || (values[field] = normalized)
    catch
    end
  end
  # case_format deliberately has no option spec (it names the input rather
  # than configuring the run); the Case page persists it in the form block,
  # so it is the one non-spec field read back here (stage 4A)
  raw_format = get(block, "case_format", nothing)
  if raw_format !== nothing
    fmt = lowercase(strip(String(raw_format)))
    # scf and pgm name the same reader; both are accepted so the choice
    # made in the form survives the save (see _normalize_case_format)
    fmt in ("auto", "matpower", "dtf_for001", "cgmes", "scf", "pgm") && (values["case_format"] = fmt)
  end
  return values
end

function webui_form_state(; selected_casefile::AbstractString = "", selected_config_file::AbstractString = "", sidecar_profile = nothing, submitted_form = nothing, case_directory = nothing)
  config_path = isempty(selected_config_file) ? DEFAULT_SPARLECTRA_CONFIG_PATH : selected_config_file
  values = Dict{String,Any}(spec.field => spec.default for spec in WEBUI_OPTION_SPECS)
  config_values = _webui_config_field_values(config_path)
  merge!(values, config_values)
  # Case levels seed the form in resolution order (deprecated in-file block,
  # then the case configuration file); no mtime logic anywhere (D5). The
  # controls have to show these values: the form posts a value for EVERY
  # field it renders, and those count as explicit overrides, so an unseeded
  # form would silently outrank the very settings it just loaded (measured:
  # a case asking for distributed slack ran without it because the untouched
  # checkbox posted false).
  case_file_values = _webui_case_config_field_values(selected_casefile, case_directory)
  merge!(values, case_file_values)
  isempty(case_file_values) || (values["_case_file_fields"] = sort!(collect(keys(case_file_values))))
  values["casefile"] = selected_casefile
  values["casefile_manual"] = ""
  values["config_file"] = config_path
  if sidecar_profile isa AbstractDict
    for (field, value) in sidecar_profile
      field == "_profile_path" && continue
      haskey(_WEBUI_OPTION_BY_FIELD, String(field)) || continue
      values[String(field)] = _webui_normalize_case_profile_form_value(String(field), value)
    end
    haskey(sidecar_profile, "_profile_path") && (values["_profile_path"] = sidecar_profile["_profile_path"])
  end
  if submitted_form isa AbstractDict
    for spec in WEBUI_OPTION_SPECS
      if spec.control == :checkbox
        values[spec.field] = _webui_form_bool(_webui_form_value(submitted_form, spec.field, false))
      elseif _webui_form_value(submitted_form, spec.field, nothing) !== nothing
        raw = _webui_form_value(submitted_form, spec.field)
        values[spec.field] = try
          _webui_normalize_case_profile_form_value(spec.field, raw)
        catch
          _webui_form_string(raw)
        end
      end
    end
    for field in ("casefile", "casefile_manual", "config_file")
      raw = _webui_form_value(submitted_form, field, nothing)
      raw === nothing || (values[field] = strip(_webui_form_string(raw)))
    end
  end
  return values
end

"""
    _webui_load_case_settings(output_root, casefile; case_directory) -> Union{Nothing,Dict}

Per-case Web UI form defaults from the case configuration file: the `form`
block (SE and generator options) plus the `_profile_path` marker for the
settings notice. Converts a leftover legacy settings sidecar
(`<stem>.sparlectra-webui.yaml`) ONCE into the case configuration file and
deletes it, with an operation-log line either way (D8 of the adapter task);
when a case configuration file already exists, the stale sidecar is
discarded instead of clobbering the newer file. The config-backed form
fields are seeded separately from the configuration levels
(`_webui_case_config_field_values`).
"""
function _webui_load_case_settings(output_root::AbstractString, casefile::AbstractString; case_directory::Union{Nothing,AbstractString} = nothing)
  path = try
    _webui_case_settings_path(output_root, casefile; case_directory)
  catch err
    _webui_log_case_settings_load(output_root, "case_settings_load_failed"; casefile, status = "rejected", message = sprint(showerror, err))
    return nothing
  end
  source = isabspath(strip(String(casefile))) ? normpath(String(casefile)) : _webui_resolve_case_profile_source(casefile; case_directory)
  legacy = try
    _webui_legacy_case_settings_path(output_root, casefile; case_directory)
  catch
    ""
  end
  if !isempty(legacy) && isfile(legacy)
    if isfile(path)
      rm(legacy; force = true)
      _webui_log_case_settings_load(output_root, "case_settings_sidecar_discarded"; casefile, profile_path = path, status = "converted", message = "legacy sidecar removed; the existing case configuration file wins")
    else
      converted_config = Dict{String,Any}()
      converted_form = Dict{String,Any}()
      try
        data = load_yaml_dict(legacy)
        settings = get(data, "settings", Dict{String,Any}())
        settings isa AbstractDict || (settings = Dict{String,Any}())
        for (k, v) in settings
          field = String(k)
          spec = get(_WEBUI_OPTION_BY_FIELD, field, nothing)
          spec === nothing && continue
          normalized = try
            _webui_normalize_case_profile_form_value(field, v)
          catch
            continue
          end
          normalized === nothing && continue
          if spec.config_key !== nothing
            converted_config[String(spec.config_key)] = normalized
          elseif spec.save_in_case_sidecar
            converted_form[field] = normalized
          end
        end
        write_case_config(source, converted_config; form = converted_form)
        rm(legacy; force = true)
        _webui_log_case_settings_load(output_root, "case_settings_sidecar_converted"; casefile, profile_path = path, status = "converted", config_keys = length(converted_config), form_fields = length(converted_form))
      catch err
        _webui_log_case_settings_load(output_root, "case_settings_load_failed"; casefile, profile_path = legacy, status = "failed", message = string("sidecar conversion failed: ", sprint(showerror, err)))
      end
    end
  end
  if !isfile(path)
    _webui_log_case_settings_load(output_root, "case_settings_not_found"; casefile, profile_path = path, status = "missing")
    return nothing
  end
  # an unreadable or invalid case configuration file must not read as
  # "settings loaded": the same header contract the run path enforces
  # (load_case_config) decides here too
  raw = try
    load_case_config(source)
    load_yaml_dict(path)
  catch err
    _webui_log_case_settings_load(output_root, "case_settings_load_failed"; casefile, profile_path = path, status = "failed", message = sprint(showerror, err))
    return nothing
  end
  profile = Dict{String,Any}()
  block = get(raw, "form", nothing)
  if block isa AbstractDict
    for (k, v) in block
      field = String(k)
      if !haskey(_WEBUI_OPTION_BY_FIELD, field)
        _webui_log_case_settings_load(output_root, "case_settings_field_ignored"; casefile, profile_path = path, status = "ignored", field, message = "unknown field")
        continue
      end
      try
        normalized = _webui_normalize_case_profile_form_value(field, v)
        normalized === nothing && continue
        profile[field] = normalized
      catch err
        _webui_log_case_settings_load(output_root, "case_settings_field_ignored"; casefile, profile_path = path, status = "ignored", field, message = sprint(showerror, err))
      end
    end
  end
  profile["_profile_path"] = path
  _webui_log_case_settings_load(output_root, "case_settings_loaded"; casefile, profile_path = path, status = "loaded", setting_count = length(profile) - 1)
  return profile
end

## Merge a subset of settings into the case's `form` block, creating a
## minimal case configuration file when the case has none yet. Used by the
## SE generator so a case reload restores the generator options (the
## run-based save path never sees the generator form). Only persisted
## fields (WEBUI_OPTION_SPECS) are accepted; existing configuration keys
## stay untouched.
function _webui_merge_case_settings!(output_root::AbstractString, casefile_path::AbstractString, updates::AbstractDict; case_directory::Union{Nothing,AbstractString} = nothing)
  source = isabspath(strip(String(casefile_path))) ? normpath(String(casefile_path)) : _webui_resolve_case_profile_source(casefile_path; case_directory)
  isempty(source) && return nothing
  existing_config = try
    load_case_config(source)
  catch err
    # an UNREADABLE case configuration must not be rewritten from a side
    # path: overwriting here silently threw away whatever the file held.
    # load_case_config returns an empty dict for a missing file, so this
    # branch only fires on real damage; the generator options are simply
    # not persisted this time.
    _webui_log_case_settings_load(output_root, "case_settings_merge_skipped"; casefile = casefile_path, profile_path = case_config_path(source), status = "skipped", message = sprint(showerror, err))
    return nothing
  end
  form = _webui_case_form_defaults(casefile_path, case_directory)
  for (k, v) in updates
    field = String(k)
    field in _WEBUI_CASE_PROFILE_FIELDS || continue
    spec = get(_WEBUI_OPTION_BY_FIELD, field, nothing)
    (spec === nothing || spec.config_key !== nothing) && continue
    form[field] = v
  end
  return write_case_config(source, existing_config; form = form)
end

function _webui_input_value(values::AbstractDict, field::AbstractString, default)::String
  return _webui_escape(_webui_form_string(get(values, field, default)))
end

function _webui_checked(values::AbstractDict, field::AbstractString, default::Bool)::String
  return _webui_parse_bool(get(values, field, default)) ? " checked" : ""
end

function _webui_selected(values::AbstractDict, field::AbstractString, default)
  return get(values, field, default)
end

function _webui_form_value(form::AbstractDict, key::String, default = nothing)
  haskey(form, key) && return form[key]
  symbol_key = Symbol(key)
  haskey(form, symbol_key) && return form[symbol_key]
  return default
end

# Format hint for the case pages (badge preselection, DTF assistance, SC
# button gating). The CONTENT-based detection is `_detect_case_format`, the
# same function `import_case` runs, so the form and the import can never
# disagree about what a resolvable file is (stage 4A review point: no second
# detection path). Only values that do not resolve to an existing path (a
# cgmes: alias not fetched yet, a free-typed name) fall back to the thin
# syntactic pre-stage below, which mirrors the detector's extension rules.
function _webui_case_format_hint(casefile::AbstractString; case_directory::Union{Nothing,AbstractString} = nothing)::Symbol
  value = strip(String(casefile))
  isempty(value) && return :auto
  startswith(lowercase(value), "cgmes:") && return :cgmes
  resolved = try
    _webui_resolve_case_profile_source(value; case_directory)
  catch
    ""
  end
  if !isempty(resolved) && (isfile(resolved) || isdir(resolved))
    return try
      _detect_case_format(resolved)
    catch
      # ambiguous .DAT: the detector refuses to guess, the form still offers
      # the DTF option (choosing it is exactly the way out the error names)
      lowercase(splitext(resolved)[2]) == ".dat" ? :dtf_for001 : :auto
    end
  end
  ext = lowercase(splitext(value)[2])
  ext == ".dat" && return :dtf_for001
  (ext in (".zip", ".xml") || (isempty(ext) && isdir(value))) && return :cgmes
  return :auto
end

# "Short circuit" button gating (the button is only
# selectable when the data is actually there): cheap byte scan over the
# delivery's XML contents for short-circuit source markers — synchronous
# machines, feeder short-circuit currents, or equivalent-injection
# impedances. In-memory via collectCGMESFiles (nested ZIPs included), cached
# per path+mtime so a 6209-bus delivery is scanned once, not on every form
# render. A scan error must NOT lock the button: the service run reports
# missing data with an explicit failure reason either way.
const _WEBUI_SC_DATA_CACHE = Dict{String,Tuple{Float64,Bool}}()
const _WEBUI_SC_DATA_CACHE_LOCK = ReentrantLock()

function _webui_case_has_short_circuit_data(casefile::AbstractString)::Bool
  path = strip(String(casefile))
  (isfile(path) || isdir(path)) || return true
  stamp = try
    mtime(path)
  catch
    0.0
  end
  return lock(_WEBUI_SC_DATA_CACHE_LOCK) do
    hit = get(_WEBUI_SC_DATA_CACHE, path, nothing)
    hit !== nothing && hit[1] == stamp && return hit[2]
    result = try
      # A Sparlectra Case Format case carries its sources in the file itself,
      # so the probe is a key scan, not a delivery walk. This only gates the
      # BUTTON; the service run still reports the real reason.
      if lowercase(splitext(path)[2]) == ".json"
        occursin("\"sc_source\"", read(path, String))
      else
      files = CGMESImporter.collectCGMESFiles(path)
      any(occursin("SynchronousMachine", f.content) || occursin("maxInitialSymShCCurrent", f.content) || occursin("EquivalentInjection.x", f.content) for f in files)
      end
    catch
      true
    end
    _WEBUI_SC_DATA_CACHE[path] = (stamp, result)
    result
  end
end

function _webui_parse_bool(value)::Bool
  return _webui_form_bool(value)
end

function _webui_parse_form_value(value, ::Type{Bool}, field::String)
  return _webui_parse_bool(value)
end

function _webui_parse_form_value(value, type::Type{<:Number}, field::String)
  text = value === nothing ? "" : strip(string(value))
  isempty(text) && throw(ArgumentError("Web UI field $(field) must not be empty."))
  try
    return parse(type, text)
  catch
    throw(ArgumentError("Web UI field $(field) has invalid value $(repr(text))."))
  end
end

_webui_parse_form_value(value, ::Type{String}, field::String) = strip(String(something(value, "")))

"""
    powerflow_webui_request(form; default_output_root="results/powerflow_service")

Convert browser form values into the dictionary accepted by
[`start_powerflow_run`](@ref). Only keys from
[`GUI_EDITABLE_CONFIG_KEYS`](@ref) are emitted as configuration overrides.

Stage 4A precedence: a request-option field the run page no longer renders
falls back to the case configuration file's `form` block before the spec
default, so values saved on the Case and Settings pages keep reaching the
run. A field present in the POST always wins (config-KEY fields need no such
fallback here: omitted overrides fall through to `resolve_config`, which
reads the same case configuration file server-side).
"""
function powerflow_webui_request(form::AbstractDict; default_output_root::AbstractString = "results/powerflow_service", case_directory::Union{Nothing,AbstractString} = nothing)::Dict{String,Any}
  existing_casefile = strip(String(something(_webui_form_value(form, "casefile", ""), "")))
  manual_casefile = strip(String(something(_webui_form_value(form, "casefile_manual", ""), "")))
  casefile = isempty(manual_casefile) ? existing_casefile : manual_casefile
  # any supported case format qualifies here, so the message must not say
  # MATPOWER; choosing happens on the Case page since stage 4A
  isempty(casefile) && throw(ArgumentError("Select a case first (Case page)."))
  stored_form = _webui_case_form_defaults(casefile, case_directory)
  config_file = strip(String(something(_webui_form_value(form, "config_file", ""), "")))
  output_root = String(default_output_root)
  ignore_webui_settings = if _webui_form_value(form, "ignore_webui_settings", nothing) !== nothing
    _webui_parse_form_value(_webui_form_value(form, "ignore_webui_settings", "false"), Bool, "ignore_webui_settings")
  else
    !_webui_parse_form_value(_webui_form_value(form, "apply_webui_runtime_overrides", "true"), Bool, "apply_webui_runtime_overrides")
  end
  apply_runtime_overrides = !ignore_webui_settings
  overrides = Dict{String,Any}()
  if apply_runtime_overrides
    for (config_key, field, type) in _WEBUI_FORM_CONFIG_FIELDS
      config_key in GUI_EDITABLE_CONFIG_KEYS || error("Web UI field $(field) is not GUI-editable.")
      spec = _webui_option_spec(field)
      _webui_form_value(form, field, nothing) === nothing && continue
      raw = _webui_form_value(form, field)
      # An EMPTY field means "not stated", not "zero": the configured value
      # stays. Without this an optional field (the tolerance in MW) could
      # only be offered by inventing a magic number for "unset", and
      # clearing any numeric field would fail the run with a parse error.
      raw isa AbstractString && isempty(strip(raw)) && continue
      overrides[config_key] = _webui_parse_form_value(raw, type, field)
    end
    # "off" in the enforcement-mode control means "no Q-limit handling". The
    # control pair (checkbox + mode) reads as one setting, and a mode picked
    # while the handling is off does nothing while looking like it does
    # (reported from a live session 2026-09-08).
    mode_key = "power_flow.qlimits.enforcement_mode"
    if haskey(overrides, mode_key) && lowercase(strip(String(overrides[mode_key]))) == "off"
      delete!(overrides, mode_key)
      overrides["power_flow.qlimits.enabled"] = false
    end
  end
  # The tolerance is ONE value with a unit, not two competing fields: the
  # unit decides which configuration key the number becomes. MW is the
  # physical statement (converted with the case base when the tolerance
  # meets the net), pu the classical one. Only one of the two keys is ever
  # sent, so a run can never carry both.
  if apply_runtime_overrides
    unit = uppercase(strip(String(something(_webui_form_value(form, "power_flow_tol_unit", "pu"), "pu"))))
    if unit == "MW"
      raw_tol = _webui_form_value(form, "power_flow_tol", nothing)
      if raw_tol !== nothing && !isempty(strip(String(raw_tol)))
        overrides["power_flow.tol_MW"] = _webui_parse_form_value(raw_tol, Float64, "power_flow_tol")
        delete!(overrides, "power_flow.tol")
      end
    end
  end
  request_options = Dict{String,Any}()
  for field in _WEBUI_CASE_PROFILE_EXTRA_FIELDS
    spec = _webui_option_spec(field)
    raw = _webui_form_value(form, field, nothing)
    # absent from the POST (the page no longer renders it; a rendered
    # checkbox always posts via its hidden false input): stored form block,
    # then spec default
    raw === nothing && (raw = get(stored_form, field, spec.default))
    request_options[field] = _webui_parse_form_value(raw, spec.value_type, field)
  end
  case_format = strip(String(something(_webui_form_value(form, "case_format", nothing), "")))
  if isempty(case_format)
    stored_format = strip(String(something(get(stored_form, "case_format", nothing), "")))
    case_format = if !isempty(stored_format)
      stored_format
    else
      # the format hint resolves through _detect_case_format; only the two
      # formats the auto path cannot settle by itself are made explicit
      hint = _webui_case_format_hint(casefile; case_directory)
      hint in (:dtf_for001, :cgmes) ? String(hint) : "auto"
    end
  end
  for002_reference_file = strip(String(something(_webui_form_value(form, "for002_reference_file", ""), "")))
  outage_mode = strip(String(something(_webui_form_value(form, "dtf_outage_selection_mode", "none"), "none")))
  raw_selection = _webui_form_value(form, "dtf_outage_selection", String[])
  selection = raw_selection isa AbstractVector ? String.(raw_selection) : (isempty(strip(String(raw_selection))) ? String[] : [String(raw_selection)])
  return Dict{String,Any}(
    "casefile" => casefile,
    "case_format" => case_format,
    "for002_reference_file" => isempty(for002_reference_file) ? nothing : for002_reference_file,
    "run_dtf_outages" => outage_mode != "none",
    "dtf_outage_selection" => selection,
    "dtf_outage_selection_mode" => outage_mode,
    "compare_for002_outages" => !isempty(for002_reference_file) && outage_mode != "none",
    "write_outage_artifacts" => _webui_parse_form_value(_webui_form_value(form, "write_outage_artifacts", "true"), Bool, "write_outage_artifacts"),
    "write_outage_matpower_exports" => _webui_parse_form_value(_webui_form_value(form, "write_outage_matpower_exports", "false"), Bool, "write_outage_matpower_exports"),
    "matpower_export_requested" => _webui_parse_form_value(_webui_form_value(form, "matpower_export_requested", "false"), Bool, "matpower_export_requested"),
    "config_file" => config_file,
    "output_root" => output_root,
    "config_overrides" => overrides,
    "config_override_source" => apply_runtime_overrides ? "webui_form_runtime" : "user_yaml",
    "performance_timing" => request_options["performance_timing"],
    # Normal "Start PowerFlow run" submissions never write diagnose.log; only
    # the dedicated "Diagnose" self-check does (diagnose_mode forces
    # run_diagnostics = true server-side in start_powerflow_run regardless of
    # this value). There is no "Run diagnostics" checkbox for a normal run.
    "run_diagnostics" => false,
    "diagnose_mode" => _webui_parse_form_value(_webui_form_value(form, "diagnose_mode", "false"), Bool, "diagnose_mode"),
    "short_circuit_mode" => _webui_parse_form_value(_webui_form_value(form, "short_circuit_mode", "false"), Bool, "short_circuit_mode"),
    "import_analysis_mode" => _webui_parse_form_value(_webui_form_value(form, "import_analysis_mode", "false"), Bool, "import_analysis_mode"),
    "contingency_mode" => _webui_parse_form_value(_webui_form_value(form, "contingency_mode", "false"), Bool, "contingency_mode"),
    # N-1 outage kind is a RUN parameter (branch / generator), read as a plain
    # request key like dtf_outage_selection, not a config override
    "contingency_kind" => strip(String(something(_webui_form_value(form, "contingency_kind", "branch"), "branch"))),
    # scenario task step 6: source, screening mode and margin are run
    # parameters; empty form values mean "not set" and keep the service
    # defaults (historical kind list, configured screening)
    "scenario_source" => (v = strip(String(something(_webui_form_value(form, "scenario_source", ""), ""))); isempty(v) ? nothing : v),
    "screening_mode" => (v = strip(String(something(_webui_form_value(form, "screening_mode", ""), ""))); isempty(v) ? nothing : v),
    "screening_margin_pct" => (v = strip(String(something(_webui_form_value(form, "screening_margin_pct", ""), ""))); isempty(v) ? nothing : _webui_parse_form_value(v, Float64, "screening_margin_pct")),
    "detailed_result_csv" => request_options["detailed_result_csv"],
    "detailed_result_csv_format" => request_options["detailed_result_csv_format"],
    "export_cgmes" => request_options["export_cgmes"],
    # state estimation (SE phase 5): its own run kind plus the SE-started
    # chain reference; all read as plain request keys like contingency_kind
    "se_mode" => _webui_parse_form_value(_webui_form_value(form, "se_mode", "false"), Bool, "se_mode"),
    "measurement_file" => strip(String(something(_webui_form_value(form, "measurement_file", ""), ""))),
    "se_max_iter" => _webui_parse_form_value(_webui_form_value(form, "se_max_iter", "30"), Int, "se_max_iter"),
    "se_tol" => _webui_parse_form_value(_webui_form_value(form, "se_tol", "1e-8"), Float64, "se_tol"),
    "se_flatstart" => _webui_parse_form_value(_webui_form_value(form, "se_flatstart", "true"), Bool, "se_flatstart"),
    "se_robust" => _webui_parse_form_value(_webui_form_value(form, "se_robust", "false"), Bool, "se_robust"),
    "se_robust_mode" => strip(String(something(_webui_form_value(form, "se_robust_mode", "off"), "off"))),
    "se_k_eliminate" => _webui_parse_form_value(_webui_form_value(form, "se_k_eliminate", "3.0"), Float64, "se_k_eliminate"),
    "se_robust_k1" => _webui_parse_form_value(_webui_form_value(form, "se_robust_k1", "3.0"), Float64, "se_robust_k1"),
    "se_robust_k2" => _webui_parse_form_value(_webui_form_value(form, "se_robust_k2", "6.0"), Float64, "se_robust_k2"),
    "se_k_suppress" => _webui_parse_form_value(_webui_form_value(form, "se_k_suppress", "6.0"), Float64, "se_k_suppress"),
    "se_suppression_sigma" => _webui_parse_form_value(_webui_form_value(form, "se_suppression_sigma", "2000"), Float64, "se_suppression_sigma"),
    "se_max_eliminations" => _webui_parse_form_value(_webui_form_value(form, "se_max_eliminations", "3"), Int, "se_max_eliminations"),
    "se_update_shunts" => _webui_parse_form_value(_webui_form_value(form, "se_update_shunts", "false"), Bool, "se_update_shunts"),
    "se_report_correlation" => _webui_parse_form_value(_webui_form_value(form, "se_report_correlation", "false"), Bool, "se_report_correlation"),
    "se_tap_estimation" => _webui_parse_form_value(_webui_form_value(form, "se_tap_estimation", "false"), Bool, "se_tap_estimation"),
    "se_start_run_id" => (v = strip(String(something(_webui_form_value(form, "se_start_run_id", ""), ""))); isempty(v) ? nothing : v),
    "se_start_mode" => strip(String(something(_webui_form_value(form, "se_start_mode", "se_state"), "se_state"))),
  )
end

## measurement CSV v1 content sniff (SE phase 5): the role is decided by the
## version comment, never by the extension alone
function _webui_is_measurement_csv(path::AbstractString)::Bool
  isfile(path) || return false
  line = try
    open(io -> readline(io; keep = false), path, "r")
  catch
    return false
  end
  return strip(line) == "# sparlectra-measurements v1"
end

## measurement sets offered on the SE page: .csv files in the case directory
## whose content sniff confirms the v1 version comment (a .csv without it is
## retained but never offered)
function _webui_measurement_options_in_directory(directory::AbstractString)
  isdir(directory) || return String[]
  names = [name for name in readdir(directory) if lowercase(splitext(name)[2]) == ".csv" && _webui_is_measurement_csv(joinpath(directory, name))]
  return sort!(names; by = lowercase)
end

## explicit case binding of a measurement set: the `# case: <name>` comment
## the generator (and the editor) write into the file. The binding lives IN
## the file, not in its name: a renamed or re-uploaded set keeps saying
## which case it belongs to. Empty for sets without one (pre-binding files).
"""
    _webui_measurement_set_provenance(path) -> Dict{String,Any}

The provenance comments of a measurement CSV as a dictionary: `noise` (a
Bool plus the descriptive text), `generator`, `seed`. Empty for a file that
carries none (an older set, or one written by hand).

`noise` is the one that matters beyond bookkeeping: a noise-free set has
J = 0 by construction, and a reader who does not know that will read the
zero as a perfect estimate.
"""
function _webui_measurement_set_provenance(path::AbstractString)::Dict{String,Any}
  out = Dict{String,Any}()
  isfile(path) || return out
  in_taps = false
  try
    for line in eachline(path)
      startswith(line, "#") || break            ## comments are the header block
      text = strip(chopprefix(line, "#"))
      if startswith(text, "noise:")
        detail = strip(chopprefix(text, "noise:"))
        out["noise"] = !startswith(detail, "none")
        out["noise_detail"] = String(detail)
      elseif startswith(text, "generator:")
        out["generator"] = String(strip(chopprefix(text, "generator:")))
      elseif startswith(text, "seed:")
        parsed = tryparse(Int, strip(chopprefix(text, "seed:")))
        parsed === nothing || (out["seed"] = parsed)
      elseif startswith(text, "flow_ends:")
        out["flow_ends"] = String(strip(chopprefix(text, "flow_ends:")))
      elseif startswith(text, "truth:")
        out["truth"] = String(strip(chopprefix(text, "truth:")))
      elseif startswith(text, "truth_value,")
        # per-row generator truth: without it a case file can never produce
        # the measured-vs-truth delta file that a CSV set produces
        parts = split(chopprefix(text, "truth_value,"), ",")
        if length(parts) >= 2
          v = tryparse(Float64, String(parts[end]))
          v === nothing || (get!(out, "truth_values", Dict{String,Any}())[String(join(parts[1:(end-1)], ","))] = v)
        end
      elseif in_taps && !startswith(text, "branch,")
        # the taps table the generator appended: only rows with a non-zero
        # generation deviation matter downstream (the estimator has to absorb
        # exactly those), and losing them is why a case file gave no tap hint
        cols = split(text, ",")
        if length(cols) >= 11
          br = tryparse(Int, String(cols[1]))
          dev = tryparse(Float64, String(cols[11]))
          if br !== nothing && dev !== nothing && dev != 0.0
            push!(get!(out, "tap_deviations", Any[]), Dict{String,Any}("branch" => br, "name" => String(cols[3]), "steps" => dev))
          end
        end
      end
      startswith(text, "sparlectra-taps") && (in_taps = true)
    end
  catch
    return Dict{String,Any}()
  end
  return out
end

function _webui_measurement_set_case(path::AbstractString)::String
  isfile(path) || return ""
  result = ""
  open(path, "r") do io
    for _ = 1:50
      eof(io) && break
      line = strip(readline(io))
      startswith(line, "#") || break
      m = match(r"^#\s*case:\s*(.+)$", line)
      if m !== nothing
        result = String(strip(m.captures[1]))
        break
      end
    end
  end
  return result
end

## documented tap deviations of a measurement set: rows of the
## `sparlectra-taps v1` comment table whose generation_deviation_steps (the
## last column) is nonzero. The SE service warns when such a set runs
## WITHOUT tap estimation (the large J is a model discrepancy, not bad data).
function _measurement_set_tap_deviations(path::AbstractString)::Vector{NamedTuple{(:branch, :steps),Tuple{Int,Float64}}}
  out = NamedTuple{(:branch, :steps),Tuple{Int,Float64}}[]
  intable = false
  # The tap table lists EVERY transformer of the case, so the deviating ones
  # can sit deep inside it: on a 25000-bus set the first deviation was line
  # 949 of the block while the reader stopped after 300, which is why the
  # deviations were never seen and tap estimation never released itself.
  # The limit is therefore large enough for the whole table; the loop still
  # stops at the first non-comment line, so it never reads the data rows.
  for line in _webui_measurement_set_comments(path; limit = 100_000)
    if line == "sparlectra-taps v1"
      intable = true
      continue
    end
    intable || continue
    startswith(line, "branch,") && continue
    parts = split(line, ",")
    length(parts) >= 2 || break
    b = tryparse(Int, String(first(parts)))
    d = tryparse(Float64, String(last(parts)))
    (b === nothing || d === nothing) && break   # end of the table block
    d != 0.0 && push!(out, (branch = b, steps = d))
  end
  return out
end

## structured data rows of a measurement CSV for the table editor: line
## number in the file plus the columns the editor may change (value, sigma,
## active) and the read-only identity (type, location, id). Limited: very
## large sets stay in the text editor / download path.
function _webui_measurement_set_rows(path::AbstractString; limit::Int = 400)
  isfile(path) || return NamedTuple[]
  out = NamedTuple[]
  for (ln, line) in enumerate(eachline(path))
    s = strip(line)
    (isempty(s) || startswith(s, "#") || startswith(s, "type,")) && continue
    parts = split(s, ","; limit = 11)
    length(parts) >= 11 || continue
    push!(out, (line = ln, typ = String(parts[1]), bus = String(parts[2]), from_bus = String(parts[3]), to_bus = String(parts[4]), branch_nr = String(parts[5]), direction = String(parts[7]), value = String(parts[8]), sigma = String(parts[9]), active = String(parts[10]), id = String(parts[11])))
    length(out) >= limit && break
  end
  return out
end

## per-type row counts of a measurement CSV (data rows only), for the
## "what is in this set" summary on the SE page
function _webui_measurement_set_counts(path::AbstractString)::Vector{Tuple{String,Int}}
  isfile(path) || return Tuple{String,Int}[]
  counts = Dict{String,Int}()
  open(path, "r") do io
    for line in eachline(io)
      s = strip(line)
      (isempty(s) || startswith(s, "#") || startswith(s, "type,")) && continue
      typ = String(first(split(s, ","; limit = 2)))
      counts[typ] = get(counts, typ, 0) + 1
    end
  end
  return sort!([(k, v) for (k, v) in counts]; by = first)
end

## leading metadata comments of a measurement CSV (between the version line
## and the header): the generator records the transformer tap positions the
## set was generated from there. Returns at most `limit` lines.
function _webui_measurement_set_comments(path::AbstractString; limit::Int = 40)::Vector{String}
  isfile(path) || return String[]
  out = String[]
  open(path, "r") do io
    first_line = true
    for line in eachline(io)
      s = strip(line)
      if first_line
        first_line = false
        continue   # version comment
      end
      startswith(s, "#") || break
      push!(out, String(strip(lstrip(s, '#'))))
      length(out) >= limit && break
    end
  end
  return out
end
