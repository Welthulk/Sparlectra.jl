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

# file: src/webui/forms.jl
# purpose: Web UI input resolution: application root and case directory
#          lookup, case/config file selectors, upload classification, and
#          form value parsing
const _WEBUI_PACKAGE_ROOT = SPARLECTRA_ROOT   # the library checkout: data, examples and docs live there, not in app/

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

The shipped cases under `data/scf` (every `.scf.json`: the `sp_*` demo
cases, the PST format fixture, and since 0.12.5 the feeder and Q-limit
example cases), every SCF time series under `data/qsts` (since 0.60.0,
recognized by its content), the plain power-grid-model files under `data/PGM`
(`.json`), the CGMES demo deliveries, and since 0.20.0 the PowSyBl IIDM
files under `data/powsybl` (`.xiidm`). They are offered in the case
chooser regardless of the cache directory and staged into the cache on
first use
([`_webui_stage_bundled_case!`](@ref)).
"""
function _webui_bundled_scf_options(application_root::AbstractString)::Vector{String}
  names = String[]
  scf_dir = joinpath(application_root, "data", "scf")
  isdir(scf_dir) && append!(names, [name for name in readdir(scf_dir) if endswith(lowercase(name), ".scf.json")])
  # SCF time series (QSTS, 0.60.0), shipped under data/qsts: the content decides, not the name
  qsts_dir = joinpath(application_root, "data", "qsts")
  isdir(qsts_dir) && append!(names, [name for name in readdir(qsts_dir) if _webui_is_regular_file(joinpath(qsts_dir, name)) && _webui_is_scf_series_cached(joinpath(qsts_dir, name))])
  pgm_dir = joinpath(application_root, "data", "PGM")
  isdir(pgm_dir) && append!(names, [name for name in readdir(pgm_dir) if endswith(lowercase(name), ".json")])
  # the CGMES deliveries exported by Sparlectra itself (data/cgmes_demo/<case>,
  # four profile files each) are offered as one ZIP per case; the ZIP is
  # packed into the case directory on first use
  cgmes_dir = joinpath(application_root, "data", "cgmes_demo")
  isdir(cgmes_dir) && append!(names, [_webui_cgmes_demo_zip_name(name) for name in readdir(cgmes_dir) if isdir(joinpath(cgmes_dir, name))])
  # the PowSyBl demo networks, one IIDM file each
  powsybl_dir = joinpath(application_root, "data", "powsybl")
  isdir(powsybl_dir) && append!(names, [name for name in readdir(powsybl_dir) if isfile(joinpath(powsybl_dir, name)) && endswith(lowercase(name), ".xiidm")])
  return sort!(names; by = lowercase)
end

const _WEBUI_CGMES_DEMO_SUFFIX = "_cgmes.zip"
_webui_cgmes_demo_zip_name(case::AbstractString)::String = string(case, _WEBUI_CGMES_DEMO_SUFFIX)

# the shipped delivery folder behind a bundled CGMES ZIP name, or nothing
function _webui_cgmes_demo_folder(application_root::AbstractString, name::AbstractString)::Union{Nothing,String}
  endswith(name, _WEBUI_CGMES_DEMO_SUFFIX) || return nothing
  case = name[1:(end - length(_WEBUI_CGMES_DEMO_SUFFIX))]
  isempty(case) && return nothing
  folder = joinpath(application_root, "data", "cgmes_demo", case)
  return isdir(folder) ? folder : nothing
end

# pack the four profile files of a shipped delivery into one ZIP (the form
# the Web UI case selector and the CGMES importer take); atomic through a
# temporary file so a half-written ZIP is never picked up
function _webui_pack_cgmes_demo!(folder::AbstractString, dest::AbstractString)::String
  files = sort!([f for f in readdir(folder) if endswith(lowercase(f), ".xml")])
  isempty(files) && error("CGMES demo delivery $(folder) carries no XML profile file")
  tmp = string(dest, ".tmp")
  open(tmp, "w") do io
    CGMESImporter.ZipArchives.ZipWriter(io) do w
      for f in files
        CGMESImporter.ZipArchives.zip_newfile(w, f)
        write(w, read(joinpath(folder, f)))
      end
    end
  end
  mv(tmp, dest; force = true)
  return dest
end

# The hash list of the shipped files: for every file the Web UI copies into a
# user's case directory (a case, its companions, a packed CGMES demo ZIP, a
# warmup workload), the SHA-256 of every committed version of that file. It
# is generated from the version history when the demo cases are built and
# checked in; here it is read as a plain file, nothing comes from git at run
# time.
const _WEBUI_SHIPPED_HASHES_FILE = joinpath("data", "shipped_hashes.toml")

# the directories a shipped file is copied from, in the order staging
# searches them; data/webui holds the warmup workloads of the provisioning
const _WEBUI_SHIPPED_SOURCE_DIRS = ["mpower", "scf", "qsts", "PGM", "powsybl", "webui"]

# A git checkout with core.autocrlf on Windows writes CRLF into the text
# files that .gitattributes does not pin (.yaml, .m, .jl, .xiidm), and a copy
# staged from such a checkout carries CRLF too. The list holds LF hashes, so
# a text file is hashed with every CRLF read as LF, on both sides (list and
# copy). A ZIP is binary and hashed as it is.
function _webui_lf_bytes(bytes::Vector{UInt8})::Vector{UInt8}
  out = Vector{UInt8}(undef, 0)
  sizehint!(out, length(bytes))
  n = length(bytes)
  for i in 1:n
    (bytes[i] == 0x0d && i < n && bytes[i + 1] == 0x0a) && continue
    push!(out, bytes[i])
  end
  return out
end

"""
    _webui_shipped_sha256(path) -> String

SHA-256 (hex) of a shipped file or its copy as the hash list counts it: a
`.zip` byte for byte, every other file with CRLF line ends read as LF, so a
copy from a Windows checkout with `core.autocrlf` matches the list.
"""
function _webui_shipped_sha256(path::AbstractString)::String
  bytes = read(path)
  endswith(lowercase(String(path)), ".zip") || (bytes = _webui_lf_bytes(bytes))
  return bytes2hex(sha256(bytes))
end

# A hash list that cannot be used: unparsable TOML, or a value that is not a
# list of SHA-256 hex strings. Raised by _webui_shipped_hashes and handled
# at start like an unreadable case directory (warning, log entry, no
# refresh).
struct _WebUIShippedHashesDamaged <: Exception
  message::String
end
Base.showerror(io::IO, err::_WebUIShippedHashesDamaged) = print(io, err.message)

"""
    _webui_shipped_hashes(application_root) -> Dict{String,Set{String}}

The hash list `data/shipped_hashes.toml` of the package at
`application_root`: file name in the case directory => SHA-256 (hex, see
[`_webui_shipped_sha256`](@ref)) of every committed version of that file,
the current one included. An absent list gives an empty table, so nothing
is refreshed. A list that does not parse or holds a value of the wrong type
raises `_WebUIShippedHashesDamaged`.
"""
function _webui_shipped_hashes(application_root::AbstractString)::Dict{String,Set{String}}
  path = joinpath(application_root, _WEBUI_SHIPPED_HASHES_FILE)
  isfile(path) || return Dict{String,Set{String}}()
  parsed = try
    TOML.parsefile(path)
  catch err
    err isa TOML.ParserError || rethrow()
    throw(_WebUIShippedHashesDamaged("$(path) does not parse: $(sprint(showerror, err))"))
  end
  table = get(parsed, "files", Dict{String,Any}())
  table isa AbstractDict || throw(_WebUIShippedHashesDamaged("$(path): [files] is not a table"))
  out = Dict{String,Set{String}}()
  for (name, hashes) in table
    ok = hashes isa AbstractVector && all(h -> h isa AbstractString && occursin(r"^[0-9a-f]{64}$", h), hashes)
    ok || throw(_WebUIShippedHashesDamaged("$(path): the entry of $(name) is not a list of SHA-256 hex strings"))
    out[String(name)] = Set{String}(String.(hashes))
  end
  return out
end

# How a shipped file under `name` is written and what its current version
# hashes to, or nothing when the package ships no such file. The copy in the
# case directory is never its own source: a writable checkout uses its own
# data/mpower as the case directory, and there the search goes on to the
# next directory.
function _webui_shipped_source(application_root::AbstractString, name::AbstractString, copy_path::AbstractString)
  folder = _webui_cgmes_demo_folder(application_root, name)
  if folder !== nothing
    # the ZIP is packed, not copied; packing is deterministic (stored
    # entries, no time stamps), so the current version hashes as packed
    write_zip! = dst -> _webui_pack_cgmes_demo!(folder, dst)
    zip_sha256 = () -> mktempdir(dir -> _webui_shipped_sha256(_webui_pack_cgmes_demo!(folder, joinpath(dir, String(name)))))
    return (write! = write_zip!, sha256 = zip_sha256)
  end
  for dir in _WEBUI_SHIPPED_SOURCE_DIRS
    src = joinpath(application_root, "data", dir, String(name))
    isfile(src) || continue
    ispath(copy_path) && samefile(src, copy_path) && continue
    return (write! = dst -> cp(src, dst), sha256 = () -> _webui_shipped_sha256(src))
  end
  return nothing
end

# `<name>.replaced-<date>`, with a counter when that name is taken (a second
# replacement on the same day must not overwrite the first kept copy)
function _webui_replaced_path(path::AbstractString, today::Date)::String
  base = string(path, ".replaced-", Dates.format(today, "yyyy-mm-dd"))
  candidate = base
  n = 1
  while ispath(candidate)
    n += 1
    candidate = string(base, "-", n)
  end
  return candidate
end

# The failures a copy step expects: a full disk or a missing permission
# (IOError, SystemError), and a second Web UI instance that moved or wrote
# the same file in between (ArgumentError from mv onto an existing name or
# from a vanished source).
_webui_expected_copy_failure(err)::Bool = err isa Base.IOError || err isa SystemError || err isa ArgumentError

# write!(path) into a fresh temporary name in the directory of `dst`; the
# caller renames it into place, so a crash or a full disk never leaves a
# truncated file under the real name (a truncated copy would count as
# changed by the user and never be refreshed again)
function _webui_write_temp!(@nospecialize(write!), dst::AbstractString)::String
  tmp = tempname(dirname(dst); cleanup = false)
  try
    write!(tmp)
  catch
    rm(tmp; force = true)
    rethrow()
  end
  return tmp
end

"""
    _webui_stage_shipped_file!(write!, current_sha256, dst; known=nothing, operation_log=nothing, today=Dates.today()) -> NamedTuple

The one copy step for a shipped file. A missing `dst` is written with
`write!` under a temporary name and renamed into place. An existing `dst`
is replaced only when `known` (the hashes of its name from the hash list)
and an `operation_log` are given and its content equals an EARLIER version
(hash as [`_webui_shipped_sha256`](@ref)): equal to the current version
(`current_sha256()`) it is left alone, and content that matches no listed
version (a copy the user changed) is never touched. The replacement is
written to a temporary file first; then the old copy is renamed to
`<name>.replaced-<date>` and the new file into place. The operation log
gets `shipped_case_file_replacing` (naming the kept copy) before the
renames and `shipped_case_file_replaced` after them, both with the old and
the new hash. A failure removes the temporary file, moves the old copy back
when it was already moved, and is logged as `shipped_case_file_replace_failed`;
an expected failure (disk, permission, a concurrent instance) returns
status `:failed`, any other is rethrown. Returns `(status, old_sha256,
new_sha256, kept_as)` with status `:copied`, `:kept`, `:current`,
`:replaced` or `:failed`.
"""
function _webui_stage_shipped_file!(@nospecialize(write!), @nospecialize(current_sha256), dst::AbstractString; known::Union{Nothing,Set{String}} = nothing, operation_log::Union{Nothing,AbstractString} = nothing, today::Date = Dates.today())
  if !isfile(dst)
    tmp = _webui_write_temp!(write!, dst)
    try
      mv(tmp, dst)
    catch err
      rm(tmp; force = true)
      # another instance staged the same file in between: its copy stays
      (err isa ArgumentError && isfile(dst)) || rethrow()
      return (status = :kept, old_sha256 = "", new_sha256 = "", kept_as = "")
    end
    return (status = :copied, old_sha256 = "", new_sha256 = "", kept_as = "")
  end
  (known === nothing || operation_log === nothing) && return (status = :kept, old_sha256 = "", new_sha256 = "", kept_as = "")
  old = _webui_shipped_sha256(dst)
  old in known || return (status = :kept, old_sha256 = old, new_sha256 = "", kept_as = "")
  new = current_sha256()
  old == new && return (status = :current, old_sha256 = old, new_sha256 = new, kept_as = "")
  kept_as = _webui_replaced_path(dst, today)
  log_fields = (route = "/powerflow", method = "START", user_action = false, file = basename(dst), old_sha256 = old, new_sha256 = new)
  tmp = ""
  moved = false
  try
    tmp = _webui_write_temp!(write!, dst)
    record_webui_operation!(operation_log, "shipped_case_file_replacing"; log_fields..., status = "replacing", kept_as = basename(kept_as))
    mv(dst, kept_as)
    moved = true
    mv(tmp, dst)
  catch err
    isempty(tmp) || rm(tmp; force = true)
    restored = false
    if moved && !ispath(dst)
      restored = try
        mv(kept_as, dst)
        true
      catch restore_err
        _webui_expected_copy_failure(restore_err) || rethrow()
        false
      end
    end
    # the kept copy is named whenever it still exists, so nothing goes
    # missing without a trace
    left = moved && !restored ? basename(kept_as) : ""
    record_webui_operation!(operation_log, "shipped_case_file_replace_failed"; log_fields..., status = "failed", kept_as = left, restored, message = sprint(showerror, err))
    _webui_expected_copy_failure(err) || rethrow()
    return (status = :failed, old_sha256 = old, new_sha256 = new, kept_as = left)
  end
  record_webui_operation!(operation_log, "shipped_case_file_replaced"; log_fields..., status = "replaced", kept_as = basename(kept_as))
  return (status = :replaced, old_sha256 = old, new_sha256 = new, kept_as)
end

"""
    _webui_refresh_shipped_copies!(application_root, case_directory, operation_log; known, today) -> Vector

Bring the copies of shipped files in `case_directory` up to the package's
version (0.50.0). Every name of the hash list ([`_webui_shipped_hashes`](@ref))
that exists in the case directory goes through
[`_webui_stage_shipped_file!`](@ref): a copy that equals an earlier
committed version is replaced (old copy kept as `<name>.replaced-<date>`,
operation-log event `shipped_case_file_replaced`), every other copy stays as
it is. Only the listed names are looked at, and only existing files are
hashed. Returns the replacements as `(file, kept_as)`.
"""
function _webui_refresh_shipped_copies!(application_root::AbstractString, case_directory::AbstractString, operation_log::AbstractString; known::Dict{String,Set{String}} = _webui_shipped_hashes(application_root), today::Date = Dates.today())
  replaced = NamedTuple{(:file, :kept_as),Tuple{String,String}}[]
  isdir(case_directory) || return replaced
  for name in sort!(collect(keys(known)))
    # a list entry is a bare file name; anything else would reach outside
    # the case directory
    (basename(name) == name && name != "." && name != "..") || continue
    dst = joinpath(case_directory, name)
    isfile(dst) || continue
    source = _webui_shipped_source(application_root, name, dst)
    source === nothing && continue
    result = _webui_stage_shipped_file!(source.write!, source.sha256, dst; known = known[name], operation_log, today)
    result.status === :replaced && push!(replaced, (file = name, kept_as = basename(result.kept_as)))
  end
  return replaced
end

"""
    _webui_stage_bundled_case!(application_root, case_directory, requested) -> Union{Nothing,String}

Resolve a bare requested case NAME against the package's bundled sources
(`data/mpower`, then `data/scf`, then `data/qsts`, then `data/PGM`, then `data/powsybl`) and
stage it into `case_directory`.
Returns the staged (or already cached) absolute path, or `nothing` when no
bundled source carries the name. An SCF demo case is staged WITH its
sidecars (`<stem>.config.yaml` pins the machine-neutral resolution, the
measurement CSVs belong to the case); existing cache files are not
overwritten here, so user-saved settings survive (a copy that still holds an
earlier shipped version is replaced at Web UI start by
[`_webui_refresh_shipped_copies!`](@ref)). With no writable
`case_directory` the bundled path itself is returned.
"""
function _webui_stage_bundled_case!(application_root::AbstractString, case_directory::Union{Nothing,AbstractString}, requested::AbstractString)::Union{Nothing,String}
  name = String(requested)
  (isabspath(name) || occursin('/', name) || occursin('\\', name)) && return nothing
  demo_folder = _webui_cgmes_demo_folder(application_root, name)
  if demo_folder !== nothing
    target_dir = case_directory === nothing ? mktempdir() : String(case_directory)
    mkpath(target_dir)
    cached = joinpath(target_dir, name)
    _webui_stage_shipped_file!(dst -> _webui_pack_cgmes_demo!(demo_folder, dst), () -> "", cached)
    return cached
  end
  bundled = _webui_bundled_case_source(application_root, name)
  bundled === nothing && return nothing
  case_directory === nothing && return bundled
  mkpath(case_directory)
  cached = joinpath(String(case_directory), name)
  if !isfile(cached)
    # the companions travel with the case, enumerated by the single
    # definition: config pin, measurement CSVs, weights
    for src in vcat([bundled], case_companion_files(bundled))
      dst = joinpath(String(case_directory), basename(src))
      _webui_stage_shipped_file!(d -> cp(src, d), () -> "", dst)
    end
  end
  return cached
end

# The shipped file behind a bare case name (`data/mpower`, then `data/scf`,
# `data/qsts`, `data/PGM`, `data/powsybl`), or nothing. Staging copies it; the format
# hint reads its content while the case is not staged yet. The CGMES demo
# ZIPs are packed on staging and have no single source file.
function _webui_bundled_case_source(application_root::AbstractString, name::AbstractString)::Union{Nothing,String}
  (isabspath(name) || occursin('/', name) || occursin('\\', name)) && return nothing
  for source_dir in (joinpath(application_root, "data", "mpower"), joinpath(application_root, "data", "scf"), joinpath(application_root, "data", "qsts"), joinpath(application_root, "data", "PGM"), joinpath(application_root, "data", "powsybl"))
    bundled = joinpath(source_dir, String(name))
    isfile(bundled) && return bundled
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
# ZIP-in-ZIP), the importer's container layer opens them in memory.
# .csv covers measurement sets (SE phase 5): the role is decided by a content
# sniff on the version comment, not the extension (see _webui_is_measurement_csv).
# .json is a Sparlectra Case Format case (#342); the import path validates the
# content before it is stored, so a foreign .json is rejected by name.
# .yaml is the case configuration sidecar (issue #378 follow-up): a case
# exported with "Save case as"/"Export as SCF" travels as three files
# (case, sidecar, measurements); re-uploading all three together must bring
# the sidecar along, or the settings the export carried are silently lost.
_webui_supported_upload_case_extension(name::AbstractString)::Bool = lowercase(splitext(basename(String(name)))[2]) in (".m", ".dat", ".zip", ".csv", ".json", ".yaml", ".xml", ".xiidm")

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

"""
    _webui_case_config_upload_reason(bytes) -> Union{Nothing,String}

Why an uploaded `.yaml` is not a usable case configuration sidecar, or
`nothing` when it is one: it must parse as YAML and declare `scope: case`.
The `case:` header is deliberately NOT required to name a file already
present - a case, its sidecar, and its measurement sets are uploaded
together in one multi-file selection, and upload order within that
selection is not guaranteed; `load_case_config` still enforces the
`case_config_mismatch` check the moment a run actually reads the file,
which is the point where "does this sidecar belong to that case" actually
matters. Parses on a temporary copy, exactly like the `.json` case-file
check, so a foreign YAML never lands in the case directory under a
`.config.yaml`-shaped name.
"""
function _webui_case_config_upload_reason(bytes::Vector{UInt8})::Union{Nothing,String}
  tmp = tempname()
  try
    write(tmp, bytes)
    raw = load_yaml_dict(tmp)
    scope = lowercase(strip(string(get(raw, "scope", ""))))
    scope == "case" || return "not a case configuration file (scope: case is missing or wrong)"
    return nothing
  catch err
    return first(split(sprint(showerror, err), '\n'))
  finally
    isfile(tmp) && rm(tmp; force = true)
  end
end

function _webui_classify_dat_content(path::AbstractString)::Symbol
  lowercase(splitext(path)[2]) == ".dat" || return :not_dat
  text = try
    read(path, String)
  catch
    return :unknown_dat
  end
  # a network deck is what the DTF reader takes as one, the rule of the
  # format detection of the library (`is_dtf_deck`); the second read only
  # looks for the outage records of a file that has passed
  if DTFImporter.is_dtf_deck(path)
    return isempty(DTFImporter.read_dtf(path; strict = false).outages) ? :dtf_network_case : :dtf_network_case_with_outages
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
    return "🔗 $(version) boundary set · $(proflist) · $(cnt(:TopologicalNode)) boundary nodes - pair it with a base case"
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
  isempty(missing_profiles) || return "❌ $(version) · $(proflist) · $(content) - incomplete: missing profile(s) $(join(String.(missing_profiles), ", "))$(_webui_upload_analysis_note(path))"

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
  return "⚠️ $(version) · $(proflist) · $(content) · $(summary.unresolved_count) unresolved references - boundary set missing: upload it as well, or type cgmes:&lt;alias&gt;$(_webui_upload_analysis_note(path))"
end

"""
Upload-time import analysis for a CGMES delivery that stayed incomplete
(missing profiles or missing boundary): writes the full
[`importFailureAnalysis`](@ref Sparlectra.CGMESImporter.importFailureAnalysis)
report next to the case as `<case>.import_analysis.txt` and returns a short
message suffix naming the missing declared dependencies. Best-effort: an
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
      return " · missing declared dependencies: $(shown) - full analysis in $(basename(report_file))"
    end
    return " - full analysis in $(basename(report_file))"
  catch
    # keep the role summary useful even when the analysis cannot be built
    return ""
  end
end

"""
Look for a boundary delivery in the local ENTSO-E test-set cache that resolves
`path`, copy it next to the case, and return its filename (empty when none
fits). Only uses what is already extracted locally, no download is triggered
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

# Whether a file is an SCF time series (QSTS), by content. is_scf_series
# raises the ArgumentError of the zip reader for a file that ends like a zip
# naming series.json but is none; such a file is no series for the selector
# (as a .zip it is still offered, and a run names the reader's error).
function _webui_is_scf_series(path::AbstractString)::Bool
  try
    return is_scf_series(path)
  catch err
    err isa ArgumentError || rethrow()
    return false
  end
end

_webui_is_scf_series_cached(path::AbstractString)::Bool = _webui_file_scan_memo(_webui_is_scf_series, path)

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
function _webui_is_user_selectable_case(name::AbstractString)::Bool
  lowered_name = lowercase(basename(name))
  _, extension = splitext(lowered_name)
  # PowSyBl: the IIDM file
  endswith(lowered_name, ".xiidm") && return true
  endswith(lowered_name, ".sparlectra-webui.yaml") && return false
  # per-case configuration files travel next to their case and are not cases
  endswith(lowered_name, ".config.yaml") && return false
  endswith(lowered_name, ".contingency-weights.csv") && return false
  # an SCF time series (QSTS) is a case whatever it is called: the content
  # decides (a zip naming series.json)
  (isfile(name) && _webui_is_scf_series_cached(name)) && return true
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
    # A boundary set is a companion delivery, not a runnable case: hide it
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

# `isfile` on a directory entry that is being deleted by another process
# raises EACCES on Windows (seen with a temporary directory of the test
# suite under data/mpower); such an entry is simply not a case file.
function _webui_is_regular_file(path::AbstractString)::Bool
  try
    return isfile(path)
  catch err
    err isa Base.IOError || rethrow(err)
    return false
  end
end

function _webui_casefile_options_in_directory(directory::AbstractString)::Vector{String}
  isdir(directory) || return String[]
  files = filter(readdir(directory)) do name
    path = joinpath(directory, name)
    return _webui_is_regular_file(path) && _webui_is_user_selectable_case(path)
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
# be converted once and deleted
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

# Shared selected-case state: the Case page is THE
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
output, runtime, webui, export scope). Request metadata entries
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

The SINGLE definition of which files belong to one case: the case
configuration file, the legacy settings sidecar, the
per-case N-1 weights list, the controller file `<series>.control.yaml` of
an SCF time series, and the measurement CSVs including the
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
  # the controllers of a QSTS series (<series>.control.yaml)
  control = _qsts_control_file_path(src)
  (control != src && isfile(control) && _webui_is_scf_series_cached(src)) && push!(companions, control)
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
  # the PowSyBl slack list is a YAML list in the configuration and one
  # semicolon-separated text in the form (the reader takes both spellings)
  if String(field) == "powsybl_import_slack_ids" && value isa AbstractVector
    return join((strip(string(v)) for v in value), "; ")
  end
  type = get(_WEBUI_CASE_PROFILE_FIELD_TYPES, String(field), String)
  # `auto` or a number shows as it is stored
  type === Union{Float64,String} && return _webui_form_string(_webui_parse_form_value(value, Union{Float64,String}, String(field)))
  # a bus list shows as comma-separated text; a stored text is checked the
  # way a posted one is
  if type === Vector{Int}
    buses = value isa AbstractVector ? _webui_parse_form_value(join(string.(value), ","), Vector{Int}, String(field)) : _webui_parse_form_value(value, Vector{Int}, String(field))
    return join(string.(buses), ", ")
  end
  value isa AbstractVector && throw(ArgumentError("Case-settings field $(field) does not support vector values for form rendering."))
  allowed = get(_WEBUI_CASE_PROFILE_SELECT_VALUES, String(field), nothing)
  if allowed !== nothing && value isa Bool
    !value && "off" in allowed && return "off"
    value && "on" in allowed && return "on"
  end
  # the zero-injection box stands for state_estimation.zero_injection, whose
  # third value constraint also means "rows present": the box shows checked
  # (the generation handler writes constraint back while it stays checked)
  String(field) == "gen_passive_as_zi" && lowercase(strip(string(value))) == "constraint" && return true
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
      normalized = try
        _webui_normalize_case_profile_form_value(field, value)
      catch err
        # a value the form cannot render (an alias the select does not
        # list, a type the control cannot show) leaves THIS field on its
        # default and the other fields on the file's values; the form
        # fallback is silent by design, the run reads the file itself.
        # Earlier one such value emptied the whole seed.
        err isa ArgumentError || rethrow()
        continue
      end
      values[field] = normalized
    end
    return values
  catch err
    # an unreadable or unparseable configuration file leaves the form on
    # the spec defaults; said aloud, a silent fallback here once hid a
    # list-valued key behind default values on every page
    @warn "Web UI: configuration file $(path) could not seed the form, showing defaults" exception = (err, catch_backtrace())
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
      merge!(merged, scf_case_config(path; dropped = String[]))
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

# --- measurement sigma by voltage level (State Estimation section) ----------
#
# The table of state_estimation.sigma_by_voltage_level on the SE section:
# one row per nominal voltage of the case, plus every configured level the
# case does not have (shown so a save does not drop it). The form posts
# se_sigma_lvl_n and per row i se_sigma_lvl_<i>_<column>; an empty cell is
# "default", not zero.
const _WEBUI_SIGMA_LEVEL_KEY = "state_estimation.sigma_by_voltage_level"
const _WEBUI_SIGMA_LEVEL_SIGMA_COLUMNS = ["p_mw", "q_mvar", "i_a", "u_kv", "rel_pct"]

# Nominal voltages per case, whether it has tap changers and its
# transformers with tap changers, cached by path and mtime: all need the
# imported net, and the sections render on every Runs and Case page load.
# The import runs outside the lock so a slow case does not block other pages.
const _WEBUI_CASE_VN_CACHE = Dict{String,Tuple{Float64,Vector{Float64},String,Bool,Vector{NamedTuple}}}()
const _WEBUI_CASE_VN_CACHE_LOCK = ReentrantLock()

"""
    _webui_case_import_facts(case_path, config_file) -> (levels, error, has_taps, trafos)

What the SE sections need from the imported case, from one cached import:
the distinct bus nominal voltages in kV, highest first, an error text
(empty on success), `has_taps`, whether the case has a transformer with
a declared tap changer (`_case_has_tap_changers`), and `trafos`, the
transformers the SE service releases for the tap estimation
(`_webui_case_tap_transformers`). A case that does not import gives no
levels, the import error, `has_taps = false` and no transformers.
"""
function _webui_case_import_facts(case_path::AbstractString, config_file::AbstractString)
  path = String(case_path)
  (isfile(path) || isdir(path)) || return (Float64[], "case file not found", false, NamedTuple[])
  stamp = mtime(path)
  hit = lock(() -> get(_WEBUI_CASE_VN_CACHE, path, nothing), _WEBUI_CASE_VN_CACHE_LOCK)
  (hit !== nothing && hit[1] == stamp) && return (hit[2], hit[3], hit[4], hit[5])
  levels, err, has_taps, trafos = try
    # the effective configuration the SE run imports the case with
    cfg = resolve_config(isempty(strip(config_file)) ? DEFAULT_SPARLECTRA_CONFIG_PATH : String(config_file), path).config
    net = import_case(path, cfg; run_kind = :state_estimation).net
    (sort!(unique(Float64[nd.comp.cVN for nd in net.nodeVec]); rev = true), "", _case_has_tap_changers(net), _webui_case_tap_transformers(net))
  catch err
    # the expected failure: the case does not import (format, content);
    # the section shows the reason instead of the levels
    err isa InterruptException && rethrow()
    (Float64[], first(split(sprint(showerror, err), '\n')), false, NamedTuple[])
  end
  lock(() -> (_WEBUI_CASE_VN_CACHE[path] = (stamp, levels, err, has_taps, trafos)), _WEBUI_CASE_VN_CACHE_LOCK)
  return (levels, err, has_taps, trafos)
end

## The transformers of a net the SE service releases for the tap estimation
## (in service, a changer of a release mode), one row each: `name`, `mrid`
## (empty outside CGMES), `machine` (a machine transformer by the rule of the
## estimator), `kinds` (the changer kinds, comma separated) and `key`, the
## entry of state_estimation.tap_estimation.exclude that names it: the name
## when no other transformer has it, else the mRID.
function _webui_case_tap_transformers(net)::Vector{NamedTuple}
  ks = Int[k for (k, br) in enumerate(net.branchVec) if br.ratio != 0.0 && _branch_terminal_state(br) == :closed && _tap_release_mode(br) !== nothing]
  mrids = _transformer_mrids(net)
  machine = _tap_exclusions(net, TapEstimationConfig(); ks = ks).reasons
  names = [getCompName(br.comp) for br in net.branchVec if br.ratio != 0.0]
  rows = NamedTuple[]
  for k in ks
    br = net.branchVec[k]
    name = getCompName(br.comp)
    mrid = get(mrids, k, "")
    set = _effective_tap_changers(br)
    kinds = set === nothing ? "" : join(unique(String(ch.kind) for ch in set.changers), ", ")
    key = (count(==(name), names) > 1 && !isempty(mrid)) ? mrid : name
    push!(rows, (branch = k, name = name, mrid = mrid, machine = haskey(machine, k), kinds = kinds, key = key))
  end
  return rows
end

"""
    _webui_case_voltage_levels(case_path, config_file) -> (levels, error)

The distinct bus nominal voltages of a case in kV, highest first, and an
error text (empty on success). A case that does not import gives no levels
and the import error, which the SE section shows next to the table.
"""
function _webui_case_voltage_levels(case_path::AbstractString, config_file::AbstractString)
  levels, err, _, _ = _webui_case_import_facts(case_path, config_file)
  return (levels, err)
end

"""
    _webui_sigma_level_table(case_path, case_directory, config_file) -> NamedTuple

The rows the SE section shows: `(rows, error, editable, stored)` with one
row per nominal voltage of the case and per configured level the case does
not have, each `(vn_kv, values, in_case)`; `values` maps a column to its
text, empty for "default". The values come from the case configuration
file only; without a table there every cell is empty, also when the
configuration file sets the key. `stored` says the case configuration has
rows: the form then posts its table even when all cells were cleared (that
clears the stored table), while an empty table that was empty on the page
posts nothing, so the run resolves the key like an API run. When the case
configuration cannot be read, `editable` is false: the section then posts
no table, and the run stops with the same named error as an API run
instead of replacing the stored value with empty cells.
"""
function _webui_sigma_level_table(case_path::AbstractString, case_directory, config_file::AbstractString)
  path = isabspath(String(case_path)) ? String(case_path) : (case_directory === nothing ? String(case_path) : joinpath(String(case_directory), String(case_path)))
  vns, err = _webui_case_voltage_levels(path, config_file)
  configured = SESigmaLevel[]
  editable = true
  try
    configured = _parse_sigma_by_voltage_level(get(load_case_config(path), _WEBUI_SIGMA_LEVEL_KEY, nothing))
  catch e
    # the expected failures: an unreadable case configuration (header,
    # scope, a key outside the case scope) or an invalid table
    # (ArgumentError), a line the YAML reader refuses (ErrorException); the
    # run names the same error, the section shows it
    (e isa ArgumentError || e isa ErrorException) || rethrow()
    editable = false
    err = isempty(err) ? sprint(showerror, e) : string(err, "; ", sprint(showerror, e))
  end
  text(x) = x === nothing ? "" : string(x)
  rows = NamedTuple[]
  for v in sort!(unique(vcat(vns, [r.vn_kv for r in configured])); rev = true)
    r = findfirst(c -> abs(c.vn_kv - v) < Sparlectra.SE_SIGMA_LEVEL_VN_TOL_KV, configured)
    c = r === nothing ? nothing : configured[r]
    values = Dict{String,String}(
      "p_mw" => c === nothing ? "" : text(c.p_mw), "q_mvar" => c === nothing ? "" : text(c.q_mvar),
      "i_a" => c === nothing ? "" : text(c.i_a), "u_kv" => c === nothing ? "" : text(c.u_kv),
      "rel_pct" => (c === nothing || c.rel_pct == 0.0) ? "" : text(c.rel_pct),
    )
    push!(rows, (vn_kv = v, values = values, in_case = any(x -> abs(x - v) < Sparlectra.SE_SIGMA_LEVEL_VN_TOL_KV, vns)))
  end
  return (rows = rows, error = err, editable = editable, stored = !isempty(configured))
end

"""
    _webui_sigma_levels_from_form(form) -> Union{Nothing,Dict{String,Any}}

The posted sigma table as the block mapping of the configuration (the
nominal voltage as the row key, written by `_sigma_level_key`, => the sigma
columns), or `nothing` when the form carries no
table (`se_sigma_lvl_n` absent: the configured value stays). A row without
any value is left out. An empty table is an empty mapping when the page
showed stored rows (`se_sigma_lvl_stored`: the user cleared them), and
`nothing` otherwise, so an untouched empty table changes no run and writes
nothing. A cell that is not a number, or a row the configuration parser
refuses, throws an `ArgumentError` naming the row.
"""
function _webui_sigma_levels_from_form(form::AbstractDict)::Union{Nothing,Dict{String,Any}}
  n_raw = _webui_form_value(form, "se_sigma_lvl_n", nothing)
  n_raw === nothing && return nothing
  n = something(tryparse(Int, strip(String(n_raw))), 0)
  out = Dict{String,Any}()
  for i in 1:n
    vn_raw = strip(String(something(_webui_form_value(form, "se_sigma_lvl_$(i)_vn_kv", ""), "")))
    isempty(vn_raw) && continue
    vn = tryparse(Float64, vn_raw)
    vn === nothing && throw(ArgumentError("$(_WEBUI_SIGMA_LEVEL_KEY): row $(i): vn_kv '$(vn_raw)' is not a number"))
    row = Dict{String,Any}()
    for col in _WEBUI_SIGMA_LEVEL_SIGMA_COLUMNS
      raw = strip(String(something(_webui_form_value(form, "se_sigma_lvl_$(i)_$(col)", ""), "")))
      isempty(raw) && continue
      x = tryparse(Float64, replace(raw, "," => "."))
      x === nothing && throw(ArgumentError("$(_WEBUI_SIGMA_LEVEL_KEY): row $(i) ($(vn) kV): $(col) '$(raw)' is not a number"))
      row[col] = x
    end
    isempty(row) && continue
    out[_sigma_level_key(vn)] = row
  end
  isempty(out) && strip(String(something(_webui_form_value(form, "se_sigma_lvl_stored", ""), ""))) != "1" && return nothing
  _parse_sigma_by_voltage_level(out)
  return out
end

# --- tap estimation per transformer (Case page) ------------------------------
#
# The table of the transformers with tap changers on the Case page, one
# "estimate" checkbox per transformer, saved into the case configuration as
# state_estimation.tap_estimation.exclude (the unchecked transformers) and
# include_machine_transformers (true when a machine transformer is checked;
# then every unchecked transformer is listed, an unchecked machine
# transformer too). The form posts se_tap_trafo_n, per row i
# se_tap_trafo_<i>_key (the exclude entry), se_tap_trafo_<i>_machine and the
# checkbox se_tap_trafo_<i>_estimate, se_tap_trafo_<i>_excluded (named in the
# configured exclude list), and se_tap_keep_<j> for each configured entry
# without a transformer in this case (kept, so a save does not drop it).
const _WEBUI_TAP_EXCLUDE_KEY = "state_estimation.tap_estimation.exclude"
const _WEBUI_TAP_MACHINE_KEY = "state_estimation.tap_estimation.include_machine_transformers"

"""
    _webui_tap_exclusion_table(case_path, case_directory, config_file) -> NamedTuple

The rows of the tap estimation table of the Case page: `(rows, keep,
include_machine, error, editable)`. `rows` has one entry per transformer the SE service releases
(`_webui_case_tap_transformers`) plus `estimate`, whether the run estimates
it under the effective configuration of the case (not named in
`state_estimation.tap_estimation.exclude`, and for a machine transformer
`include_machine_transformers` set) and `excluded`, whether it is named in
the exclude list. `keep` lists the configured entries
that name no transformer of the case. When the configuration of the case
cannot be read, `editable` is false and `error` names the reason; the page
then posts no table.
"""
function _webui_tap_exclusion_table(case_path::AbstractString, case_directory, config_file::AbstractString)
  path = isabspath(String(case_path)) ? String(case_path) : (case_directory === nothing ? String(case_path) : joinpath(String(case_directory), String(case_path)))
  _, err, _, trafos = _webui_case_import_facts(path, config_file)
  te = TapEstimationConfig()
  editable = true
  try
    te = resolve_config(isempty(strip(config_file)) ? DEFAULT_SPARLECTRA_CONFIG_PATH : String(config_file), path).config.state_estimation.tap_estimation
  catch e
    # the expected failures: a configuration file or case configuration that
    # does not load (ArgumentError from the validators, ErrorException from
    # the YAML reader); the run names the same error, the table shows it
    (e isa ArgumentError || e isa ErrorException) || rethrow()
    editable = false
    err = isempty(err) ? sprint(showerror, e) : string(err, "; ", sprint(showerror, e))
  end
  wanted = Set{String}(te.exclude)
  named(t) = t.name in wanted || (!isempty(t.mrid) && t.mrid in wanted)
  rows = NamedTuple[merge(t, (estimate = !named(t) && (!t.machine || te.include_machine_transformers), excluded = named(t))) for t in trafos]
  known = Set{String}(vcat([t.name for t in trafos], [t.mrid for t in trafos if !isempty(t.mrid)]))
  keep = String[e for e in te.exclude if !(e in known)]
  return (rows = rows, keep = keep, include_machine = te.include_machine_transformers, error = err, editable = editable)
end

"""
    _webui_tap_exclusions_from_form(form) -> Union{Nothing,NamedTuple}

The posted tap estimation table as `(exclude, include_machine)`, the values
of `state_estimation.tap_estimation.exclude` and
`include_machine_transformers`, or `nothing` when the form carries no table
(`se_tap_trafo_n` absent: the configured values stay). A checked machine
transformer sets `include_machine`; every unchecked transformer then goes
into `exclude`. Without one, the unchecked machine transformers are left
out of `exclude` (the default keeps them out and calculates their step),
except those already named in the configured list
(`se_tap_trafo_<i>_excluded`): a user exclusion is not calculated, and a
save must not turn it into the machine default. A table without machine
transformers keeps the configured `include_machine`
(`se_tap_machine_configured`). The kept entries (`se_tap_keep_<j>`) come
first.
"""
function _webui_tap_exclusions_from_form(form::AbstractDict)::Union{Nothing,NamedTuple}
  n_raw = _webui_form_value(form, "se_tap_trafo_n", nothing)
  n_raw === nothing && return nothing
  n = something(tryparse(Int, strip(String(n_raw))), 0)
  text(field) = strip(String(something(_webui_form_value(form, field, ""), "")))
  rows = [(key = text("se_tap_trafo_$(i)_key"), machine = text("se_tap_trafo_$(i)_machine") == "true", estimate = text("se_tap_trafo_$(i)_estimate") == "true", excluded = text("se_tap_trafo_$(i)_excluded") == "true") for i in 1:n]
  filter!(r -> !isempty(r.key), rows)
  include_machine = any(r -> r.machine, rows) ? any(r -> r.machine && r.estimate, rows) : text("se_tap_machine_configured") == "true"
  exclude = String[]
  j = 1
  while (k = text("se_tap_keep_$(j)")) != ""
    push!(exclude, k)
    j += 1
  end
  for r in rows
    r.estimate && continue
    (r.machine && !include_machine && !r.excluded) && continue
    push!(exclude, r.key)
  end
  return (exclude = unique(exclude), include_machine = include_machine)
end

# The case file (or CGMES directory) behind a case name, or nothing (no
# name, no such path). A
# bare name not yet in the case directory is looked up in the bundled
# sources the combobox offers, the way the run resolves it before staging.
function _webui_case_source_file(casefile::AbstractString, case_directory; application_root::AbstractString = _webui_application_root())::Union{Nothing,String}
  name = String(strip(String(casefile)))
  isempty(name) && return nothing
  path = isabspath(name) ? name : (case_directory === nothing ? name : joinpath(String(case_directory), name))
  if !isfile(path) && !isabspath(name) && !occursin('/', name) && !occursin('\\', name)
    bundled = _webui_bundled_case_source(application_root, name)
    bundled === nothing || (path = bundled)
  end
  return ispath(path) ? path : nothing
end

# (stated_hz, can_use) of an SCF/PGM document: its meta.f_nom (nothing when
# absent or invalid, the import reports that), and whether the import turns
# a capacitance into a susceptance with the frequency: an si file with a
# line c1 != 0 (the condition that sets net._frequency_used, see
# _scf_case_frequency); a pu file reads b_pu as it stands. Raw JSON instead
# of an SCFCase keeps the check cheap for every Case page render.
function _webui_scf_frequency_facts(text::AbstractString)::Tuple{Union{Nothing,Float64},Bool}
  root = try
    scf_json_parse(text)
  catch err
    # a malformed file leaves the select as it was; the run names the parse
    # error (Web UI form seeding falls back silently)
    err isa ArgumentError || rethrow()
    return (nothing, true)
  end
  root isa AbstractDict || return (nothing, true)
  spar = get(root, "sparlectra", nothing)
  meta = spar isa AbstractDict ? get(spar, "meta", nothing) : nothing
  stated = nothing
  units = "si"
  if meta isa AbstractDict
    v = get(meta, "f_nom", nothing)
    (v isa Real && !(v isa Bool) && isfinite(v) && v > 0) && (stated = Float64(v))
    u = get(meta, "units", "si")
    u isa AbstractString && (units = lowercase(strip(u)))
  end
  data = get(root, "data", nothing)
  lines = data isa AbstractDict ? get(data, "line", nothing) : nothing
  has_c1(row) = row isa AbstractDict && (c = get(row, "c1", nothing); c isa Real && !(c isa Bool) && c != 0)
  can_use = units == "si" && lines isa AbstractVector && any(has_c1, lines)
  return (stated, can_use)
end

# (stated_hz, can_use) of a case file by content: an SCF/PGM .json file, or
# the base case of an SCF time series zip. Every other format never uses the
# key: MATPOWER gives line charging in per unit, DTF, CGMES and PowSyBl as
# a susceptance; none of their imports converts a capacitance with
# model.nominal_frequency_hz (it only becomes the meta.f_nom of an SCF
# export, which the re-import cancels).
function _webui_case_frequency_facts(path::AbstractString)::Tuple{Union{Nothing,Float64},Bool}
  isfile(path) || return (nothing, false)
  lowercase(splitext(path)[2]) == ".json" && return _webui_scf_frequency_facts(read(path, String))
  if _webui_is_scf_series_cached(path)
    entries = try
      _scf_series_read_entries(path)
    catch err
      # a broken series zip (missing or extra entry): the run names it
      err isa ArgumentError || rethrow()
      return (nothing, true)
    end
    return _webui_scf_frequency_facts(entries["base.scf.json"])
  end
  return (nothing, false)
end

_webui_case_frequency_facts_cached(path::AbstractString) = _webui_file_scan_memo(_webui_case_frequency_facts, path)

# (stated_hz, can_use) for a case name; a name that resolves to no file yet
# says nothing, the select stays enabled.
function _webui_case_frequency_state(casefile::AbstractString, case_directory; application_root::AbstractString = _webui_application_root())::Tuple{Union{Nothing,Float64},Bool}
  path = _webui_case_source_file(casefile, case_directory; application_root)
  path === nothing && return (nothing, true)
  return _webui_case_frequency_facts_cached(path)
end

# Whether the system frequency select is locked for a case: the file states
# its frequency (it wins on import), or no capacitance of the case uses it.
function _webui_case_frequency_locked(casefile::AbstractString, case_directory; application_root::AbstractString = _webui_application_root())::Bool
  stated, can_use = _webui_case_frequency_state(casefile, case_directory; application_root)
  return stated !== nothing || !can_use
end

# meta.f_nom of an SCF/PGM case file or of the base case of an SCF time
# series, or nothing.
_webui_case_file_frequency(casefile::AbstractString, case_directory; application_root::AbstractString = _webui_application_root())::Union{Nothing,Float64} = first(_webui_case_frequency_state(casefile, case_directory; application_root))

# the machine-scope keys an old case file carries in its in-file block:
# ignored by the run and named in the notice so nobody expects them to act
function _webui_case_config_dropped_keys(casefile::AbstractString, case_directory)::Vector{String}
  name = strip(String(casefile))
  isempty(name) && return String[]
  path = isabspath(name) ? name : (case_directory === nothing ? name : joinpath(String(case_directory), name))
  (lowercase(splitext(path)[2]) == ".json" && isfile(path)) || return String[]
  dropped = String[]
  try
    scf_case_config(path; dropped)
  catch
    return String[]   # an unreadable block is reported by the run itself
  end
  return dropped
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
    raw = load_yaml_dict(cfg_path)
    declared = strip(string(get(raw, "case", "")))
    # a shared-stem file of the other case carries the other case's form
    (isempty(declared) || declared == basename(path)) ? get(raw, "form", nothing) : nothing
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
  # so it is the one non-spec field read back here
  raw_format = get(block, "case_format", nothing)
  if raw_format !== nothing
    fmt = lowercase(strip(String(raw_format)))
    # scf and pgm name the same reader; both are accepted so the choice
    # made in the form survives the save (see _normalize_case_format)
    if fmt in ("auto", "matpower", "dtf_for001", "cgmes", "scf", "pgm", "powsybl")
      conflict = _webui_case_format_conflict(path, fmt)
      if conflict === nothing
        values["case_format"] = fmt
      else
        # the stored value contradicts the file it was stored for: it is
        # not loaded (the run falls back to auto), and the reason travels
        # under the underscore keys so the Case page and the run log can say
        # so. The file is left alone; the next save writes what the form
        # shows.
        values["_case_format_conflict"] = conflict
        values["_case_format_stored"] = fmt
      end
    end
  end
  return values
end

# The stored `case_format` judged against the CONTENT of the case it was
# stored for (the library's `_case_format_conflict`, the same verdict the
# import gives). Seen on Windows: `scf` saved for `sp_case118.m`
# sent the MATPOWER file into the JSON reader, and the state estimation
# failed with `invalid integer ""` while the Case page showed nothing wrong.
function _webui_case_format_conflict(path::AbstractString, fmt::AbstractString)::Union{Nothing,String}
  # pgm is a spelling of scf, not a second reader (see _normalize_case_format)
  requested = fmt == "pgm" ? :scf : Symbol(fmt)
  return try
    _case_format_conflict(path, requested)
  catch err
    # the verdict reads the first bytes of the file; a file that cannot be
    # read (permissions, a vanished cache entry) gets no verdict here, the
    # run reports it with its own message
    err isa SystemError || rethrow()
    nothing
  end
end

# Whether the selected case has a case configuration file next to it, the
# test resolve_config uses to take case-scope keys from the defaults.
function _webui_case_has_config_file(casefile::AbstractString, case_directory)::Bool
  name = strip(String(casefile))
  isempty(name) && return false
  path = isabspath(name) ? name : (case_directory === nothing ? name : joinpath(String(case_directory), name))
  (isfile(path) || isdir(path)) || return false
  return isfile(case_config_path(path))
end

function webui_form_state(; selected_casefile::AbstractString = "", selected_config_file::AbstractString = "", sidecar_profile = nothing, submitted_form = nothing, case_directory = nothing, apply_case_levels::Bool = true, case_scope_from_defaults::Bool = false)
  config_path = isempty(selected_config_file) ? DEFAULT_SPARLECTRA_CONFIG_PATH : selected_config_file
  values = Dict{String,Any}(spec.field => spec.default for spec in WEBUI_OPTION_SPECS)
  config_values = _webui_config_field_values(config_path)
  # A case with a case configuration file reads its case-scope keys from
  # that file over the packaged defaults, never from the configuration file
  # (resolve_config). The Settings page asks for that view
  # (case_scope_from_defaults): it showed the configuration file's value
  # (L20: max_iter 25 shown, 80 used), and a save of the unchanged page
  # wrote it into the case file and changed the run.
  if case_scope_from_defaults && apply_case_levels && _webui_case_has_config_file(selected_casefile, case_directory)
    default_values = _webui_config_field_values(DEFAULT_SPARLECTRA_CONFIG_PATH)
    for (config_key, field, _) in _WEBUI_FORM_CONFIG_FIELDS
      scf_is_case_config_key(config_key) || continue
      haskey(default_values, field) ? (config_values[field] = default_values[field]) : delete!(config_values, field)
    end
  end
  merge!(values, config_values)
  # an SCF file states its system frequency (meta.f_nom), and the import
  # gives it precedence over the key: the select shows it, so an untouched
  # form posts the file's value and never triggers the mismatch warning.
  # A case configuration file below still wins (the user saved a value).
  if apply_case_levels
    case_f = _webui_case_file_frequency(selected_casefile, case_directory)
    case_f === nothing || (values["model_nominal_frequency_hz"] = case_f)
  end
  # Case levels seed the form in resolution order (deprecated in-file block,
  # then the case configuration file); no mtime logic anywhere. The
  # controls have to show these values: the form posts a value for EVERY
  # field it renders, and those count as explicit overrides, so an unseeded
  # form would silently outrank the very settings it just loaded (measured:
  # a case asking for distributed slack ran without it because the untouched
  # checkbox posted false).
  # The Settings page shows the configuration file alone unless the case
  # view is requested (apply_case_levels = false): its values reach a run
  # through the configuration precedence, not through a posted form, so an
  # unseeded form cannot outrank the case there.
  case_file_values = apply_case_levels ? _webui_case_config_field_values(selected_casefile, case_directory) : Dict{String,Any}()
  merge!(values, case_file_values)
  isempty(case_file_values) || (values["_case_file_fields"] = sort!(collect(keys(case_file_values))))
  dropped_keys = apply_case_levels ? _webui_case_config_dropped_keys(selected_casefile, case_directory) : String[]
  isempty(dropped_keys) || (values["_case_file_dropped"] = dropped_keys)
  values["casefile"] = selected_casefile
  values["casefile_manual"] = ""
  values["config_file"] = config_path
  if sidecar_profile isa AbstractDict && apply_case_levels
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
    _webui_load_case_settings(output_root, casefile; config_file, case_directory) -> Union{Nothing,Dict}

Per-case Web UI form defaults from the case configuration file: the `form`
block (SE and generator options) plus the `_profile_path` marker for the
settings notice. Converts a leftover legacy settings sidecar
(`<stem>.sparlectra-webui.yaml`) ONCE into the case configuration file and
deletes it, with an operation-log line either way;
when a case configuration file already exists, the stale sidecar is
discarded instead of clobbering the newer file. The conversion takes over
the case-scope keys in which `config_file` differs from the packaged
defaults, as a save that creates the file does. The config-backed form
fields are seeded separately from the configuration levels
(`_webui_case_config_field_values`).
"""
function _webui_load_case_settings(output_root::AbstractString, casefile::AbstractString; config_file::AbstractString, case_directory::Union{Nothing,AbstractString} = nothing)
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
        # the conversion creates the case file: variant A as in a save
        adoption = _webui_case_adoption(config_file, source; skip = collect(keys(converted_config)))
        merge!(converted_config, adoption)
        write_case_config(source, converted_config; form = converted_form)
        rm(legacy; force = true)
        _webui_log_case_settings_load(output_root, "case_settings_sidecar_converted"; casefile, profile_path = path, status = "converted", config_keys = length(converted_config), form_fields = length(converted_form), adopted_keys = length(adoption))
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
      # the measurement-set choice of the SE section has no option spec (it
      # is a file name, not an option); it is kept as a plain name so the
      # section arms that set again
      if field == "measurement_file"
        name = basename(strip(String(v)))
        isempty(name) || (profile[field] = name)
        continue
      end
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
## fields (WEBUI_OPTION_SPECS) are accepted; a config-backed field among
## them is written as its case-scope configuration key, every other
## existing configuration key stays untouched. A merge that creates the
## file takes over the case-scope keys in which `config_file` differs from
## the packaged defaults (variant A, _webui_case_adoption). Returns
## `(path, adopted)`, or `nothing` when nothing was written.
function _webui_merge_case_settings!(output_root::AbstractString, casefile_path::AbstractString, updates::AbstractDict; config_file::AbstractString, case_directory::Union{Nothing,AbstractString} = nothing)
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
  # underscore markers (a stored case_format that contradicts the file) are
  # bookkeeping for the pages, not settings; a merge drops the contradicting
  # value from the file along with them
  filter!(kv -> !startswith(first(kv), "_"), form)
  config = Dict{String,Any}(existing_config)
  for (k, v) in updates
    field = String(k)
    field in _WEBUI_CASE_PROFILE_FIELDS || continue
    spec = get(_WEBUI_OPTION_BY_FIELD, field, nothing)
    spec === nothing && continue
    if spec.config_key !== nothing
      # a config-backed field is a configuration key of the case, never a
      # form-block entry (the generator's zero-injection checkbox)
      key = String(spec.config_key)
      scf_is_case_config_key(key) || continue
      config[key] = v
      delete!(form, field)
      continue
    end
    form[field] = v
  end
  adopted = String[]
  if !isfile(case_config_path(source))
    for (key, installed) in _webui_case_adoption(config_file, source; skip = collect(keys(config)))
      config[key] = installed
      push!(adopted, key)
    end
    sort!(adopted)
  end
  return (; path = write_case_config(source, config; form = form), adopted)
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
# disagree about what a resolvable file is (no second detection
# path). Only values that do not resolve to an existing path (a
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
  # a shipped case that is not staged yet has its content in the package
  # (an SCF series zip would otherwise be taken for CGMES by its extension)
  if !isempty(resolved) && !ispath(resolved)
    bundled = _webui_bundled_case_source(_webui_application_root(), value)
    bundled === nothing || (resolved = bundled)
  end
  if !isempty(resolved) && (isfile(resolved) || isdir(resolved))
    return try
      _detect_case_format(resolved)
    catch err
      # a file the detection refuses by name (a .DAT file that is no
      # network deck): no format is preselected for it, the content has
      # decided and the extension does not overrule it
      err isa ArgumentError || rethrow()
      :auto
    end
  end
  ext = lowercase(splitext(value)[2])
  # a name that resolves to no file has no content to judge; the extension
  # is the only hint there is
  ext == ".dat" && return :dtf_for001
  # PowSyBl before the CGMES rule: a .powsybl bundle is a directory too
  ext in (".powsybl", ".xiidm") && return :powsybl
  (ext in (".zip", ".xml") || (isempty(ext) && isdir(value))) && return :cgmes
  return :auto
end

# "Short circuit" button gating (the button is only
# selectable when the data is actually there): cheap byte scan over the
# delivery's XML contents for short-circuit source markers, synchronous
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
      elseif endswith(lowercase(path), ".xiidm")
        # the sources are the generators: the short-circuit extension, or
        # at least a rated power the default reactance can stand on
        text = read(path, String)
        occursin("generatorShortCircuit", text) || occursin(r"<(?:[A-Za-z0-9]+:)?generator\b[^>]*\bratedS=\"", text)
      elseif endswith(lowercase(path), ".powsybl")
        generators = Sparlectra.read_powsybl_bundle(path).generators
        usable = name -> haskey(generators, name) && any(x -> isfinite(x) && x > 0.0, generators[name])
        usable(:direct_subtrans_x) || usable(:rated_s)
      else
      files = CGMESImporter.collectCGMESFiles(path)
      # the quantities an impedance comes from, the rule of the run: a
      # machine with its reactance or at least its rated power, a feeder
      # with its short-circuit current, an equivalent with its reactance,
      # a motor with its locked-rotor ratio
      any(occursin("SynchronousMachine.satDirectSubtransX", f.content) || occursin("RotatingMachine.ratedS", f.content) || occursin("maxInitialSymShCCurrent", f.content) || occursin("EquivalentInjection.x", f.content) || occursin("AsynchronousMachine.iaIrRatio", f.content) for f in files)
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

# `auto` or a number (power_flow.qlimits.final_q_accept_pu): the word stays
# a word, a number is stored as a number, anything else is refused with the
# field name
function _webui_parse_form_value(value, ::Type{Union{Float64,String}}, field::String)
  text = value === nothing ? "" : strip(string(value))
  (isempty(text) || lowercase(text) == "auto") && return "auto"
  parsed = tryparse(Float64, text)
  parsed === nothing && throw(ArgumentError("Web UI field $(field) takes auto or a number; got $(repr(String(text)))."))
  return parsed
end

# a bus list (power_flow.qlimits.trace_buses, lock_pv_to_pq_buses): comma- or
# space-separated integers >= 1, empty for none; anything else is refused
# with the field name, the way a bad number is
function _webui_parse_form_value(value, ::Type{Vector{Int}}, field::String)
  text = value === nothing ? "" : strip(string(value))
  isempty(text) && return Int[]
  buses = Int[]
  for item in split(text, r"[,\s;]+"; keepempty = false)
    bus = tryparse(Int, item)
    (bus === nothing || bus < 1) && throw(ArgumentError("Web UI field $(field) needs bus numbers (integers >= 1, comma-separated); got $(repr(String(item)))."))
    push!(buses, bus)
  end
  return buses
end

"""
    _webui_apply_qlimits_off!(updates) -> updates

"off" in the enforcement-mode control means "no Q-limit handling". The
control pair (checkbox + mode) reads as one setting, and a mode picked
while the handling is off does nothing while looking like it does
(as reported). The word is not a mode the
configuration knows, so every path that turns form fields into
configuration keys (the run request, the settings page, the case options)
replaces it here by `enabled = false` and drops the mode key; a saved
"off" reached the loader as an `ArgumentError` on the next run
(as reported).
"""
function _webui_apply_qlimits_off!(updates::AbstractDict)
  mode_key = "power_flow.qlimits.enforcement_mode"
  haskey(updates, mode_key) || return updates
  lowercase(strip(String(updates[mode_key]))) == "off" || return updates
  delete!(updates, mode_key)
  updates["power_flow.qlimits.enabled"] = false
  return updates
end

"""
    powerflow_webui_request(form; default_output_root="results/powerflow_service")

Convert browser form values into the dictionary accepted by
[`start_powerflow_run`](@ref). Only keys from
[`GUI_EDITABLE_CONFIG_KEYS`](@ref) are emitted as configuration overrides.

Precedence: a request-option field the run page no longer renders
falls back to the case configuration file's `form` block before the spec
default, so values saved on the Case and Settings pages keep reaching the
run. A field present in the POST always wins (config-KEY fields need no such
fallback here: omitted overrides fall through to `resolve_config`, which
reads the same case configuration file server-side).
"""
function powerflow_webui_request(form::AbstractDict; default_output_root::AbstractString = "results/powerflow_service", case_directory::Union{Nothing,AbstractString} = nothing, operation_log::AbstractString = default_output_root)::Dict{String,Any}
  existing_casefile = strip(String(something(_webui_form_value(form, "casefile", ""), "")))
  manual_casefile = strip(String(something(_webui_form_value(form, "casefile_manual", ""), "")))
  casefile = isempty(manual_casefile) ? existing_casefile : manual_casefile
  # any supported case format qualifies here, so the message must not say
  # MATPOWER; choosing happens on the Case page
  isempty(casefile) && throw(ArgumentError("Select a case first (Case page)."))
  stored_form = _webui_case_form_defaults(casefile, case_directory)
  # a stored format that contradicts the case file is not applied (see
  # _webui_case_form_defaults); the operation log says so, because the
  # user sees a run under `auto` that the case settings do not explain
  if haskey(stored_form, "_case_format_conflict")
    record_webui_operation!(operation_log, "case_format_ignored"; route = "/powerflow/run", method = "POST", user_action = false,
      casefile, stored_format = String(stored_form["_case_format_stored"]), status = "ignored", message = String(stored_form["_case_format_conflict"]))
  end
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
    _webui_apply_qlimits_off!(overrides)
    # the sigma table of the SE section: a table, not a scalar option
    # field, so it is read here; absent from the POST, the configuration
    # (case file first) stays
    sigma_levels = _webui_sigma_levels_from_form(form)
    sigma_levels === nothing || (overrides[_WEBUI_SIGMA_LEVEL_KEY] = sigma_levels)
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
      hint in (:dtf_for001, :cgmes, :powsybl) ? String(hint) : "auto"
    end
  end
  for002_reference_file = strip(String(something(_webui_form_value(form, "for002_reference_file", ""), "")))
  outage_mode = strip(String(something(_webui_form_value(form, "dtf_outage_selection_mode", "none"), "none")))
  raw_selection = _webui_form_value(form, "dtf_outage_selection", String[])
  selection = raw_selection isa AbstractVector ? String.(raw_selection) : (isempty(strip(String(raw_selection))) ? String[] : [String(raw_selection)])
  # the Scenarios button (shown for a case file with a scenarios block) is a
  # contingency run on the file's block and nothing else
  scenario_run = _webui_parse_form_value(_webui_form_value(form, "scenario_run", "false"), Bool, "scenario_run")
  # The QSTS button posts qsts_mode (only when it is the submitter; a
  # disabled button, the case without a series, never does). The mode
  # comes from the select next to it; the button's own value is the
  # fallback when the select is missing. The service checks the value.
  qsts_button = _webui_form_value(form, "qsts_mode", nothing)
  qsts_mode = if qsts_button === nothing || isempty(strip(String(qsts_button)))
    nothing
  else
    choice = strip(String(something(_webui_form_value(form, "qsts_mode_choice", ""), "")))
    isempty(choice) ? strip(String(qsts_button)) : choice
  end
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
    "qsts_mode" => qsts_mode === nothing ? nothing : String(qsts_mode),
    "contingency_mode" => scenario_run || _webui_parse_form_value(_webui_form_value(form, "contingency_mode", "false"), Bool, "contingency_mode"),
    # N-1 outage kind is a RUN parameter (branch / generator), read as a plain
    # request key like dtf_outage_selection, not a config override
    "contingency_kind" => strip(String(something(_webui_form_value(form, "contingency_kind", "branch"), "branch"))),
    # source, screening mode and margin are run
    # parameters; empty form values mean "not set" and keep the service
    # defaults (historical kind list, configured screening)
    "scenario_source" => scenario_run ? "file_block" : (v = strip(String(something(_webui_form_value(form, "scenario_source", ""), ""))); isempty(v) ? nothing : v),
    "screening_mode" => (v = strip(String(something(_webui_form_value(form, "screening_mode", ""), ""))); isempty(v) ? nothing : v),
    "screening_margin_pct" => (v = strip(String(something(_webui_form_value(form, "screening_margin_pct", ""), ""))); isempty(v) ? nothing : _webui_parse_form_value(v, Float64, "screening_margin_pct")),
    "detailed_result_csv" => request_options["detailed_result_csv"],
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
    "se_max_eliminations" => _webui_parse_form_value(_webui_form_value(form, "se_max_eliminations", "5"), Int, "se_max_eliminations"),
    "se_update_shunts" => _webui_parse_form_value(_webui_form_value(form, "se_update_shunts", "false"), Bool, "se_update_shunts"),
    "se_report_correlation" => _webui_parse_form_value(_webui_form_value(form, "se_report_correlation", "false"), Bool, "se_report_correlation"),
    "se_start_run_id" => (v = strip(String(something(_webui_form_value(form, "se_start_run_id", ""), ""))); isempty(v) ? nothing : v),
    "se_start_mode" => strip(String(something(_webui_form_value(form, "se_start_mode", "se_state"), "se_state"))),
  )
end

## measurement CSV v1 content sniff (SE phase 5): the role is decided by the
## version comment, never by the extension alone
# The name of a measurement set the generator or the noise action writes:
# `<stem>.<yyyymmdd-HHMMSS>.csv`, with a counter when two writes fall into
# the same second, so every write is a new file next to the earlier ones.
function _webui_stamped_measurement_name(directory::AbstractString, stem::AbstractString; stamp::AbstractString = Dates.format(Dates.now(), "yyyymmdd-HHMMSS"))::String
  name = string(stem, ".", stamp, ".csv")
  k = 1
  while isfile(joinpath(directory, name))
    k += 1
    name = string(stem, ".", stamp, "-", k, ".csv")
  end
  return name
end

function _webui_is_measurement_csv(path::AbstractString)::Bool
  isfile(path) || return false
  line = try
    open(io -> readline(io; keep = false), path, "r")
  catch
    return false
  end
  return strip(chopprefix(line, "﻿")) == "# sparlectra-measurements v1"
end

## a measurement file the SE page offers, shows and hands out: a measurement
## CSV v1 or a measurement table (version line plus header, isMeasurementTable).
## The editors and the noise action stay CSV v1 only.
_webui_is_measurement_file(path::AbstractString)::Bool = _webui_is_measurement_csv(path) || isMeasurementTable(path)

## measurement sets offered on the SE page: .csv files in the case directory
## whose content sniff confirms the v1 version comment or the measurement
## table (a .csv without either is retained but never offered)
function _webui_measurement_options_in_directory(directory::AbstractString)
  isdir(directory) || return String[]
  names = [name for name in readdir(directory) if lowercase(splitext(name)[2]) == ".csv" && _webui_is_measurement_file(joinpath(directory, name))]
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

"""
    _webui_measurement_set_binding(path) -> NamedTuple{(:case, :models)}

The binding lines in the header comments of a measurement file: `case` from
`# case: <case file name>` (empty without one) and `models` from every
`# cgmes_model: <id>` line, each in the form of `cgmesModelKey` (empty
without any). A leading BOM goes first, blank lines before or between the
header comments are skipped, CRLF is taken; reading stops at the first line
that is no comment, after at most 50 lines.
"""
function _webui_measurement_set_binding(path::AbstractString)::NamedTuple{(:case, :models),Tuple{String,Vector{String}}}
  case = ""
  models = String[]
  isfile(path) || return (case = case, models = models)
  open(path, "r") do io
    for i = 1:50
      eof(io) && break
      # strip keeps U+FEFF, so the BOM is removed explicitly
      line = strip(i == 1 ? lstrip(readline(io), '﻿') : readline(io))
      isempty(line) && continue
      startswith(line, "#") || break
      m = match(r"^#\s*case:\s*(.+)$", line)
      if m !== nothing
        isempty(case) && (case = String(strip(m.captures[1])))
        continue
      end
      m = match(r"^#\s*cgmes_model:\s*(\S+)", line)
      m === nothing || push!(models, CGMESImporter.cgmesModelKey(m.captures[1]))
    end
  end
  return (case = case, models = unique!(models))
end

## the explicit case binding of a measurement file (`# case:` line), empty
## without one; the SE service checks it before a run
_webui_measurement_set_case(path::AbstractString)::String = _webui_measurement_set_binding(path).case

# The EQ model ids of a CGMES case in the case directory, read from the file
# headers. The ArgumentError is the ZIP reader's refusal of a file that ends
# in .zip but is no ZIP (a placeholder, a broken download): such a case has no
# model id, and its run names the reader's error.
function _webui_case_cgmes_model_ids(path::AbstractString)::Vector{String}
  lowercase(splitext(path)[2]) == ".zip" || return String[]
  isfile(path) || return String[]
  try
    return CGMESImporter.cgmesEquipmentModelIds(path)
  catch err
    err isa ArgumentError || rethrow()
    return String[]
  end
end

_webui_case_cgmes_model_ids_cached(path::AbstractString)::Vector{String} = _webui_file_scan_memo(_webui_case_cgmes_model_ids, path)

## the stems a measurement set name may start with to belong to a case: the
## case name without its extension, and for a compound extension
## (`sp_case14.scf.json`) also without the inner one, because the shipped sets
## are named `sp_case14.measurements.csv` and generated ones
## `sp_case14.scf.measurements.<stamp>.csv`
function _webui_case_stems(casefile::AbstractString)::Vector{String}
  stem = String(first(splitext(basename(String(casefile)))))
  stems = [stem]
  inner, ext = splitext(stem)
  lowercase(ext) in (".scf", ".pgm") && !isempty(inner) && push!(stems, String(inner))
  return stems
end

"""
    _webui_measurement_sets_for_case(directory, casefile; also_bound_to) -> Vector{NamedTuple{(:name, :case, :model)}}

The measurement sets in `directory` that belong to `casefile`, sorted by
name. In this order:

1. A set with a `# case:` line belongs to the case when the line names the
   case or one of `also_bound_to` (the source an SCF case was exported
   from, which the SE service accepts as well); `case` is that name. The
   line wins over everything else, so the assign switch of the upload can
   bind any table to any case.
2. A set without a `# case:` line but with `# cgmes_model:` lines belongs to
   a CGMES case (a ZIP in `directory`) whose EQ header carries one of these
   model ids; `model` is the matching id. Several scenarios built on the
   same EQ all list the set.
3. A set without either line belongs to the case when its file name starts
   with a case stem followed by a dot (`sp_case14.measurements.csv` for
   `sp_case14.scf.json`).

Only the header lines of each file and the header of each EQ file are read.
"""
function _webui_measurement_sets_for_case(directory::AbstractString, casefile::AbstractString; also_bound_to::AbstractVector{<:AbstractString} = String[])
  out = NamedTuple{(:name, :case, :model),Tuple{String,String,String}}[]
  isempty(casefile) && return out
  accepted = Set{String}(vcat([basename(String(casefile))], [String(c) for c in also_bound_to if !isempty(c)]))
  stems = _webui_case_stems(casefile)
  # read lazily: only a set with a model line needs the case's EQ ids
  case_models = nothing
  for name in _webui_measurement_options_in_directory(directory)
    binding = _webui_measurement_set_binding(joinpath(directory, name))
    if !isempty(binding.case)
      binding.case in accepted && push!(out, (name = name, case = binding.case, model = ""))
    elseif !isempty(binding.models)
      case_models === nothing && (case_models = _webui_case_cgmes_model_ids_cached(joinpath(directory, basename(String(casefile)))))
      hit = findfirst(in(case_models), binding.models)
      hit === nothing || push!(out, (name = name, case = "", model = binding.models[hit]))
    elseif any(s -> startswith(name, s * "."), stems)
      push!(out, (name = name, case = "", model = ""))
    end
  end
  return out
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
## the delimiter of a measurement CSV, read off its header line (a file
## written under excel_de carries semicolons); comma when no header is found
function _webui_measurement_csv_delimiter(lines)::Char
  for line in lines
    s = strip(line)
    (startswith(s, "type,") || startswith(s, "type;")) && return Sparlectra.result_csv_delimiter(s)
  end
  return ','
end

_webui_measurement_header_line(s::AbstractString)::Bool = startswith(s, "type,") || startswith(s, "type;")

function _webui_measurement_set_rows(path::AbstractString; limit::Int = 400)
  isfile(path) || return NamedTuple[]
  isMeasurementTable(path) && return _webui_measurement_table_rows(path; limit = limit)
  _webui_is_measurement_csv(path) || return NamedTuple[]
  out = NamedTuple[]
  delimiter = _webui_measurement_csv_delimiter(eachline(path))
  for (ln, line) in enumerate(eachline(path))
    s = strip(line)
    (isempty(s) || startswith(s, "#") || _webui_measurement_header_line(s)) && continue
    parts = split(s, delimiter; limit = 11)
    length(parts) >= 11 || continue
    bus = String(parts[2])
    location = isempty(bus) ? string(parts[3], "-", parts[4], " #", parts[5], " ", parts[7]) : string("bus ", bus)
    push!(out, (line = ln, typ = String(parts[1]), location = location, value = String(parts[8]), sigma = String(parts[9]), active = String(parts[10]), id = String(parts[11])))
    length(out) >= limit && break
  end
  return out
end

## the editable rows of a measurement table for the value editor: the file
## line, the type (an Analog row takes its type from the model and shows as
## "Analog"), the key as location, value, sigma and valid as written. Rows
## the reader cannot split (wrong field count) are not offered. The
## ArgumentError is the reader's refusal of the whole file.
function _webui_measurement_table_rows(path::AbstractString; limit::Int = 400)
  rows = try
    readMeasurementTable(path)
  catch err
    err isa ArgumentError || rethrow()
    return NamedTuple[]
  end
  out = NamedTuple[]
  for r in rows
    r.problem === :field_count && continue
    location = !isempty(r.node_mrid) ? string("node ", r.node_mrid) : !isempty(r.terminal_mrid) ? string("terminal ", r.terminal_mrid) : string("Analog ", r.measurement_mrid)
    push!(out, (line = r.line, typ = isempty(r.type) ? "Analog" : r.type, location = location, value = r.value_text, sigma = r.sigma_text, active = r.valid ? "true" : "false", id = ""))
    length(out) >= limit && break
  end
  return out
end

## per-type row counts of a measurement CSV (data rows only), for the
## "what is in this set" summary on the SE page
function _webui_measurement_set_counts(path::AbstractString)::Vector{Tuple{String,Int}}
  isfile(path) || return Tuple{String,Int}[]
  isMeasurementTable(path) && return _webui_measurement_table_counts(path)
  counts = Dict{String,Int}()
  delimiter = _webui_measurement_csv_delimiter(eachline(path))
  open(path, "r") do io
    for line in eachline(io)
      s = strip(line)
      (isempty(s) || startswith(s, "#") || _webui_measurement_header_line(s)) && continue
      typ = String(first(split(s, delimiter; limit = 2)))
      counts[typ] = get(counts, typ, 0) + 1
    end
  end
  return sort!([(k, v) for (k, v) in counts]; by = first)
end

## per-type row counts of a measurement table (column `type`; a row keyed
## by measurement_mrid takes its type from the Analog and counts as
## "Analog"). The ArgumentError is the reader's refusal of a table of
## another version: the page then shows no counts, the run names the reason.
function _webui_measurement_table_counts(path::AbstractString)::Vector{Tuple{String,Int}}
  rows = try
    readMeasurementTable(path)
  catch err
    err isa ArgumentError || rethrow()
    return Tuple{String,Int}[]
  end
  counts = Dict{String,Int}()
  for r in rows
    t = isempty(r.type) ? "Analog" : r.type
    counts[t] = get(counts, t, 0) + 1
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
