const WEBUI_OPERATION_LOG_FILENAME = "webui_operations.jsonl"
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
      "timestamp" => Dates.format(Dates.now(), dateformat"yyyy-mm-ddTHH:MM:SS.sss"),
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
