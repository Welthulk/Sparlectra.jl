_api_transport_value(value::Symbol) = String(value)
_api_transport_value(value::Nothing) = nothing
_api_transport_value(value::Missing) = nothing
_api_transport_value(value::AbstractString) = String(value)
_api_transport_value(value::Bool) = value
_api_transport_value(value::Integer) = Int(value)
_api_transport_value(value::AbstractFloat) = isfinite(value) ? Float64(value) : nothing
_api_transport_value(value::NamedTuple) = Dict(String(key) => _api_transport_value(item) for (key, item) in pairs(value))
_api_transport_value(value::AbstractDict) = Dict(String(key) => _api_transport_value(item) for (key, item) in value)
_api_transport_value(value::Tuple) = [_api_transport_value(item) for item in value]
_api_transport_value(value::AbstractVector) = [_api_transport_value(item) for item in value]
_api_transport_value(value) = string(value)

"""Convert artifact metadata to a transport-safe dictionary."""
function to_dict(artifact::SparlectraApiArtifact)::Dict{String,Any}
  return Dict{String,Any}(
    "name" => artifact.name,
    "kind" => String(artifact.kind),
    "path" => artifact.path,
    "mime_type" => artifact.mime_type,
    "exists" => artifact.exists,
    "size_bytes" => artifact.size_bytes,
    "description" => artifact.description,
  )
end

"""
    to_dict(result::SparlectraApiResult; include_raw_result=false) -> Dict{String,Any}

Convert an API result to a transport-safe dictionary. `raw_result` is omitted
by default because the solved network is not a stable JSON/YAML transport type.
"""
function to_dict(result::SparlectraApiResult; include_raw_result::Bool = false)::Dict{String,Any}
  data = Dict{String,Any}(
    "run_id" => result.run_id,
    "schema_version" => result.schema_version,
    "status" => String(result.status),
    "success" => result.success,
    "converged" => result.converged,
    "solution_available" => result.solution_available,
    "iterations" => result.iterations,
    "final_mismatch" => result.final_mismatch,
    "reason" => result.reason,
    "message" => result.message,
    "casefile" => result.casefile,
    "config_file" => result.config_file,
    "output_dir" => result.output_dir,
    "logfile" => result.logfile,
    "result_file" => result.result_file,
    "artifacts" => [to_dict(artifact) for artifact in result.artifacts],
  )
  include_raw_result && (data["raw_result"] = _api_transport_value(result.raw_result))
  return data
end

"""Convert an API result to a transport-safe `NamedTuple`."""
function to_namedtuple(result::SparlectraApiResult; include_raw_result::Bool = false)::NamedTuple
  data = to_dict(result; include_raw_result = include_raw_result)
  ordered = (:run_id, :schema_version, :status, :success, :converged, :solution_available, :iterations, :final_mismatch, :reason, :message, :casefile, :config_file, :output_dir, :logfile, :result_file, :artifacts)
  base = NamedTuple{ordered}(Tuple(data[String(key)] for key in ordered))
  return include_raw_result ? merge(base, (raw_result = data["raw_result"],)) : base
end

function _json_escape(value::AbstractString)::String
  io = IOBuffer()
  for char in value
    char == '"' ? print(io, "\\\"") : char == '\\' ? print(io, "\\\\") : char == '\n' ? print(io, "\\n") : char == '\r' ? print(io, "\\r") : char == '\t' ? print(io, "\\t") : print(io, char)
  end
  return String(take!(io))
end

function _write_json(io::IO, value)
  if value === nothing || value === missing
    print(io, "null")
  elseif value isa Bool
    print(io, value ? "true" : "false")
  elseif value isa Integer
    print(io, value)
  elseif value isa AbstractFloat
    isfinite(value) ? print(io, value) : print(io, "null")
  elseif value isa Symbol || value isa AbstractString
    print(io, '"', _json_escape(string(value)), '"')
  elseif value isa NamedTuple
    _write_json(io, Dict(String(key) => item for (key, item) in pairs(value)))
  elseif value isa AbstractDict
    print(io, '{')
    first_item = true
    for key in sort!(collect(keys(value)); by = string)
      first_item || print(io, ',')
      first_item = false
      _write_json(io, String(key))
      print(io, ':')
      _write_json(io, value[key])
    end
    print(io, '}')
  elseif value isa Tuple || value isa AbstractVector
    print(io, '[')
    for (index, item) in enumerate(value)
      index == 1 || print(io, ',')
      _write_json(io, item)
    end
    print(io, ']')
  else
    _write_json(io, string(value))
  end
  return nothing
end

"""Serialize an API result as JSON without an additional JSON package."""
function to_json(result::SparlectraApiResult; include_raw_result::Bool = false)::String
  io = IOBuffer()
  _write_json(io, to_dict(result; include_raw_result = include_raw_result))
  return String(take!(io))
end

function _yaml_scalar(value)::String
  value === nothing && return "null"
  value isa Bool && return value ? "true" : "false"
  value isa Integer && return string(value)
  value isa AbstractFloat && return isfinite(value) ? string(value) : "null"
  value isa Symbol && return String(value)
  text = string(value)
  isempty(text) && return "\"\""
  occursin(r"[:#\[\]{},&*!|>'\"%@`]|^[-?]|\s$|^\s", text) && return string('"', replace(text, "\\" => "\\\\", "\"" => "\\\""), '"')
  return text
end

function _write_yaml(io::IO, value; indent::Int = 0)
  prefix = " "^indent
  if value isa AbstractDict
    for key in sort!(collect(keys(value)); by = string)
      item = value[key]
      if item isa AbstractDict
        println(io, prefix, key, ":")
        _write_yaml(io, item; indent = indent + 2)
      elseif item isa AbstractVector && all(entry -> !(entry isa AbstractDict || entry isa AbstractVector), item)
        println(io, prefix, key, ": [", join((_yaml_scalar(entry) for entry in item), ", "), "]")
      elseif item isa AbstractVector
        println(io, prefix, key, ":")
        _write_yaml(io, item; indent = indent + 2)
      else
        println(io, prefix, key, ": ", _yaml_scalar(item))
      end
    end
  elseif value isa AbstractVector
    for item in value
      if item isa AbstractDict || item isa AbstractVector
        println(io, prefix, "-")
        _write_yaml(io, item; indent = indent + 2)
      else
        println(io, prefix, "- ", _yaml_scalar(item))
      end
    end
  else
    println(io, prefix, _yaml_scalar(value))
  end
  return nothing
end

"""Serialize an API result as YAML without an additional YAML package."""
function to_yaml(result::SparlectraApiResult; include_raw_result::Bool = false)::String
  io = IOBuffer()
  _write_yaml(io, to_dict(result; include_raw_result = include_raw_result))
  return String(take!(io))
end

function _write_yaml_file(path::AbstractString, value)
  open(path, "w") do io
    _write_yaml(io, value)
  end
  return String(path)
end
