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

const WEBUI_OPERATION_LOG_FILENAME = "webui_operations.jsonl"
const WEBUI_OPERATION_LOG_RETENTION_DAYS = 10
const WEBUI_OPERATION_LOG_MAX_ENTRIES = 5_000
const WEBUI_OPERATION_LOG_KEEP_ENTRIES = 1_000
const _WEBUI_OPERATION_LOG_MAX_BYTES = 1 * 1024 * 1024
const _WEBUI_OPERATION_LOG_LOCK = ReentrantLock()

function webui_operation_log_path(output_root::AbstractString)::String
  path = abspath(output_root)
  return basename(path) == WEBUI_OPERATION_LOG_FILENAME ? path : joinpath(path, WEBUI_OPERATION_LOG_FILENAME)
end

function _webui_operation_log_int(name::AbstractString, default::Integer; minimum::Integer = 0)::Int
  raw = get(ENV, name, nothing)
  raw === nothing && return Int(default)
  parsed = tryparse(Int, strip(raw))
  parsed === nothing && return Int(default)
  return max(Int(minimum), parsed)
end

function _webui_operation_log_options(; retention_days = nothing, max_bytes = nothing, max_entries = nothing, keep_entries = nothing)::NamedTuple
  resolved_max_entries = max_entries === nothing ? _webui_operation_log_int("SPARLECTRA_WEBUI_OPERATION_LOG_MAX_ENTRIES", WEBUI_OPERATION_LOG_MAX_ENTRIES) : Int(max_entries)
  resolved_keep_entries = keep_entries === nothing ? _webui_operation_log_int("SPARLECTRA_WEBUI_OPERATION_LOG_KEEP_ENTRIES", WEBUI_OPERATION_LOG_KEEP_ENTRIES) : Int(keep_entries)
  resolved_max_entries = max(0, resolved_max_entries)
  resolved_keep_entries = clamp(resolved_keep_entries, 0, resolved_max_entries)
  return (
    retention_days = retention_days === nothing ? _webui_operation_log_int("SPARLECTRA_WEBUI_OPERATION_LOG_RETENTION_DAYS", WEBUI_OPERATION_LOG_RETENTION_DAYS) : max(0, Int(retention_days)),
    max_bytes = max_bytes === nothing ? _webui_operation_log_int("SPARLECTRA_WEBUI_OPERATION_LOG_MAX_BYTES", _WEBUI_OPERATION_LOG_MAX_BYTES; minimum = 1) : max(1, Int(max_bytes)),
    max_entries = resolved_max_entries,
    keep_entries = resolved_keep_entries,
  )
end

function _webui_operation_entry_time(entry)::Union{Nothing,Dates.DateTime}
  timestamp = entry isa AbstractDict ? get(entry, "timestamp", nothing) : nothing
  timestamp isa AbstractString || return nothing
  text = replace(String(timestamp), "Z" => "")
  for fmt in (dateformat"yyyy-mm-ddTHH:MM:SS.sss", dateformat"yyyy-mm-ddTHH:MM:SS")
    parsed = try
      Dates.DateTime(text, fmt)
    catch
      nothing
    end
    parsed === nothing || return parsed
  end
  return nothing
end

function _write_webui_operation_log_lines!(path::AbstractString, lines::Vector{String})::Bool
  temporary_path = string(path, ".tmp")
  try
    open(temporary_path, "w") do io
      for line in lines
        println(io, line)
      end
    end
    mv(temporary_path, path; force = true)
    return true
  catch
    ispath(temporary_path) && rm(temporary_path; force = true)
    rethrow()
  end
end

function _rotate_webui_operation_log!(path::AbstractString; max_bytes::Integer = _WEBUI_OPERATION_LOG_MAX_BYTES)
  isfile(path) && filesize(path) > max_bytes || return nothing
  rotated = string(path, ".1")
  ispath(rotated) && rm(rotated; force = true)
  mv(path, rotated)
  return nothing
end

function _prune_webui_operation_log!(
  path::AbstractString;
  retention_days::Integer = WEBUI_OPERATION_LOG_RETENTION_DAYS,
  max_bytes::Integer = _WEBUI_OPERATION_LOG_MAX_BYTES,
  max_entries::Integer = WEBUI_OPERATION_LOG_MAX_ENTRIES,
  keep_entries::Integer = WEBUI_OPERATION_LOG_KEEP_ENTRIES,
)::Bool
  max_entries >= 0 || throw(ArgumentError("max_entries must be nonnegative."))
  0 <= keep_entries <= max_entries || throw(ArgumentError("keep_entries must be between zero and max_entries."))
  isfile(path) || return false

  cutoff = Dates.now(Dates.UTC) - Dates.Day(max(0, retention_days))
  valid_lines = String[]
  changed = false
  for line in eachline(path)
    isempty(strip(line)) && (changed = true; continue)
    try
      entry = _parse_service_json(line)
      entry_time = _webui_operation_entry_time(entry)
      if entry_time !== nothing && entry_time >= cutoff
        push!(valid_lines, line)
      else
        changed = true
      end
    catch
      changed = true
    end
  end

  if length(valid_lines) > max_entries
    keep = min(keep_entries, length(valid_lines))
    valid_lines = keep == 0 ? String[] : valid_lines[(end - keep + 1):end]
    changed = true
  end

  if changed
    _write_webui_operation_log_lines!(path, valid_lines)
  end

  if isfile(path) && filesize(path) > max_bytes
    keep = min(keep_entries, length(valid_lines))
    if keep > 0
      _write_webui_operation_log_lines!(path, valid_lines[(end - keep + 1):end])
    else
      _write_webui_operation_log_lines!(path, String[])
    end
    isfile(path) && filesize(path) > max_bytes && _rotate_webui_operation_log!(path; max_bytes)
    changed = true
  end
  return changed
end

function _compact_webui_operation_log!(path::AbstractString; max_entries::Integer = WEBUI_OPERATION_LOG_MAX_ENTRIES, keep_entries::Integer = WEBUI_OPERATION_LOG_KEEP_ENTRIES)::Bool
  max_entries >= 0 || throw(ArgumentError("max_entries must be nonnegative."))
  0 <= keep_entries <= max_entries || throw(ArgumentError("keep_entries must be between zero and max_entries."))
  isfile(path) || return false

  valid_lines = String[]
  for line in eachline(path)
    isempty(strip(line)) && continue
    try
      _parse_service_json(line)
      push!(valid_lines, line)
    catch
    end
  end
  length(valid_lines) > max_entries || return false

  kept_lines = keep_entries == 0 ? String[] : valid_lines[(end - keep_entries + 1):end]
  return _write_webui_operation_log_lines!(path, kept_lines)
end

_webui_operation_timestamp() = string(Dates.format(Dates.now(Dates.UTC), dateformat"yyyy-mm-ddTHH:MM:SS.sss"), 'Z')

"""
    record_webui_operation!(output_root, event; fields...) -> Bool

Append one machine-readable Web UI support event beneath `output_root`. Logging
is best effort: file or serialization failures are reported as warnings and
never escape into request handling.
"""
function record_webui_operation!(output_root::AbstractString, event::AbstractString; fields...)::Bool
  try
    path = webui_operation_log_path(output_root)
    entry = Dict{String,Any}(
      "timestamp" => _webui_operation_timestamp(),
      "sparlectra_version" => string(version()),
      "sparlectra_package_path" => _sparlectra_package_path(),
      "sparlectra_git_commit" => _sparlectra_git_commit_sha(),
      "event" => String(event),
    )
    for (key, value) in fields
      value === nothing || (entry[String(key)] = value)
    end
    lock(_WEBUI_OPERATION_LOG_LOCK) do
      mkpath(dirname(path))
      options = _webui_operation_log_options()
      # Operation logging is best effort and must never break request handling.
      # Retention is time-based first. Size and entry caps are only safety guards
      # for unusually chatty sessions or malformed logs.
      isfile(path) && filesize(path) > options.max_bytes && _prune_webui_operation_log!(path; options...)
      open(path, "a") do io
        _write_json(io, entry)
        println(io)
      end
    end
    return true
  catch err
    @warn "Could not record Web UI operation" event exception = (err, catch_backtrace())
    return false
  end
end

function _webui_log_route!(output_root::AbstractString, event::AbstractString, method::AbstractString, route::AbstractString; status = nothing, fields...)
  return record_webui_operation!(
    output_root,
    event;
    method = uppercase(String(method)),
    route = String(route),
    status,
    user_action = true,
    fields...,
  )
end
