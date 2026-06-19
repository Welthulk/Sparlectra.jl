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
const WEBUI_OPERATION_LOG_MAX_ENTRIES = 10_000
const WEBUI_OPERATION_LOG_KEEP_ENTRIES = 1_000
const _WEBUI_OPERATION_LOG_MAX_BYTES = 10 * 1024 * 1024
const _WEBUI_OPERATION_LOG_LOCK = ReentrantLock()

function webui_operation_log_path(output_root::AbstractString)::String
  path = abspath(output_root)
  return basename(path) == WEBUI_OPERATION_LOG_FILENAME ? path : joinpath(path, WEBUI_OPERATION_LOG_FILENAME)
end

function _rotate_webui_operation_log!(path::AbstractString)
  isfile(path) && filesize(path) >= _WEBUI_OPERATION_LOG_MAX_BYTES || return nothing
  rotated = string(path, ".1")
  ispath(rotated) && rm(rotated; force = true)
  mv(path, rotated)
  return nothing
end

function _compact_webui_operation_log!(
  path::AbstractString;
  max_entries::Integer = WEBUI_OPERATION_LOG_MAX_ENTRIES,
  keep_entries::Integer = WEBUI_OPERATION_LOG_KEEP_ENTRIES,
)::Bool
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
  temporary_path = string(path, ".tmp")
  try
    open(temporary_path, "w") do io
      for line in kept_lines
        println(io, line)
      end
    end
    mv(temporary_path, path; force = true)
  catch
    ispath(temporary_path) && rm(temporary_path; force = true)
    rethrow()
  end
  return true
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
      _rotate_webui_operation_log!(path)
      open(path, "a") do io
        _write_json(io, entry)
        println(io)
      end
      _compact_webui_operation_log!(path)
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
