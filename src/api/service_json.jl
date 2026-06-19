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

mutable struct _ServiceJsonParser
  text::String
  index::Int
end

function _service_failure(reason::AbstractString, message::AbstractString; run_id = nothing)::Dict{String,Any}
  failure = Dict{String,Any}(
    "status" => "failed",
    "success" => false,
    "reason" => String(reason),
    "message" => String(message),
  )
  run_id === nothing || (failure["run_id"] = String(run_id))
  return failure
end

function _service_request_value(request::AbstractDict, key::String, default = nothing)
  haskey(request, key) && return request[key]
  symbol_key = Symbol(key)
  haskey(request, symbol_key) && return request[symbol_key]
  return default
end

function _skip_service_json_whitespace!(parser::_ServiceJsonParser)
  while parser.index <= lastindex(parser.text) && isspace(parser.text[parser.index])
    parser.index = nextind(parser.text, parser.index)
  end
  return nothing
end

function _parse_service_json_string!(parser::_ServiceJsonParser)::String
  parser.text[parser.index] == '"' || error("Expected JSON string")
  parser.index = nextind(parser.text, parser.index)
  io = IOBuffer()
  while parser.index <= lastindex(parser.text)
    char = parser.text[parser.index]
    parser.index = nextind(parser.text, parser.index)
    char == '"' && return String(take!(io))
    if char == '\\'
      parser.index <= lastindex(parser.text) || error("Incomplete JSON escape")
      escaped = parser.text[parser.index]
      parser.index = nextind(parser.text, parser.index)
      escaped == '"' ? print(io, '"') :
        escaped == '\\' ? print(io, '\\') :
        escaped == '/' ? print(io, '/') :
        escaped == 'b' ? print(io, '\b') :
        escaped == 'f' ? print(io, '\f') :
        escaped == 'n' ? print(io, '\n') :
        escaped == 'r' ? print(io, '\r') :
        escaped == 't' ? print(io, '\t') :
        escaped == 'u' ? _parse_service_json_unicode_escape!(parser, io) :
        error("Unsupported JSON escape \\$(escaped)")
    else
      print(io, char)
    end
  end
  error("Unterminated JSON string")
end

function _parse_service_json_unicode_escape!(parser::_ServiceJsonParser, io::IO)
  digits = IOBuffer()
  for _ in 1:4
    parser.index <= lastindex(parser.text) || error("Incomplete JSON Unicode escape")
    char = parser.text[parser.index]
    isxdigit(char) || error("Invalid JSON Unicode escape")
    print(digits, char)
    parser.index = nextind(parser.text, parser.index)
  end
  print(io, Char(parse(UInt32, String(take!(digits)); base = 16)))
  return nothing
end

function _parse_service_json_number!(parser::_ServiceJsonParser)
  start = parser.index
  while parser.index <= lastindex(parser.text) && parser.text[parser.index] in ('-', '+', '.', 'e', 'E', '0':'9'...)
    parser.index = nextind(parser.text, parser.index)
  end
  stop = prevind(parser.text, parser.index)
  token = parser.text[start:stop]
  return occursin(r"[.eE]", token) ? parse(Float64, token) : parse(Int, token)
end

function _parse_service_json_literal!(parser::_ServiceJsonParser, literal::String, value)
  stop = nextind(parser.text, parser.index, length(literal))
  parser.text[parser.index:prevind(parser.text, stop)] == literal || error("Invalid JSON literal")
  parser.index = stop
  return value
end

function _parse_service_json_array!(parser::_ServiceJsonParser)::Vector{Any}
  parser.index = nextind(parser.text, parser.index)
  values = Any[]
  _skip_service_json_whitespace!(parser)
  if parser.index <= lastindex(parser.text) && parser.text[parser.index] == ']'
    parser.index = nextind(parser.text, parser.index)
    return values
  end
  while true
    push!(values, _parse_service_json_value!(parser))
    _skip_service_json_whitespace!(parser)
    parser.index <= lastindex(parser.text) || error("Unterminated JSON array")
    separator = parser.text[parser.index]
    parser.index = nextind(parser.text, parser.index)
    separator == ']' && return values
    separator == ',' || error("Expected comma in JSON array")
  end
end

function _parse_service_json_object!(parser::_ServiceJsonParser)::Dict{String,Any}
  parser.index = nextind(parser.text, parser.index)
  values = Dict{String,Any}()
  _skip_service_json_whitespace!(parser)
  if parser.index <= lastindex(parser.text) && parser.text[parser.index] == '}'
    parser.index = nextind(parser.text, parser.index)
    return values
  end
  while true
    _skip_service_json_whitespace!(parser)
    key = _parse_service_json_string!(parser)
    _skip_service_json_whitespace!(parser)
    parser.index <= lastindex(parser.text) && parser.text[parser.index] == ':' || error("Expected colon in JSON object")
    parser.index = nextind(parser.text, parser.index)
    values[key] = _parse_service_json_value!(parser)
    _skip_service_json_whitespace!(parser)
    parser.index <= lastindex(parser.text) || error("Unterminated JSON object")
    separator = parser.text[parser.index]
    parser.index = nextind(parser.text, parser.index)
    separator == '}' && return values
    separator == ',' || error("Expected comma in JSON object")
  end
end

function _parse_service_json_value!(parser::_ServiceJsonParser)
  _skip_service_json_whitespace!(parser)
  parser.index <= lastindex(parser.text) || error("Unexpected end of JSON input")
  char = parser.text[parser.index]
  char == '{' && return _parse_service_json_object!(parser)
  char == '[' && return _parse_service_json_array!(parser)
  char == '"' && return _parse_service_json_string!(parser)
  char == 't' && return _parse_service_json_literal!(parser, "true", true)
  char == 'f' && return _parse_service_json_literal!(parser, "false", false)
  char == 'n' && return _parse_service_json_literal!(parser, "null", nothing)
  (char == '-' || isdigit(char)) && return _parse_service_json_number!(parser)
  error("Unsupported JSON value")
end

function _parse_service_json(text::AbstractString)
  source = String(text)
  isempty(strip(source)) && error("Empty JSON input")
  parser = _ServiceJsonParser(source, firstindex(source))
  value = _parse_service_json_value!(parser)
  _skip_service_json_whitespace!(parser)
  parser.index > lastindex(source) || error("Unexpected trailing JSON content")
  return value
end

function _optional_string(value)
  value === nothing && return nothing
  value isa AbstractString || error("Expected a string or null")
  return String(value)
end

function _optional_bool(value)
  value === nothing && return nothing
  value isa Bool || error("Expected a boolean or null")
  return value
end

function _optional_int(value)
  value === nothing && return nothing
  value isa Integer || error("Expected an integer or null")
  return Int(value)
end

function _optional_float(value)
  value === nothing && return nothing
  value isa Real || error("Expected a number or null")
  return Float64(value)
end

function _reconstruct_powerflow_result(data::AbstractDict, paths)::SparlectraApiResult
  String(get(data, "run_id", "")) == paths.run_id || error("result.json run_id does not match its index entry")
  status = get(data, "status", nothing)
  status isa AbstractString || error("result.json status must be a string")
  success = get(data, "success", nothing)
  success isa Bool || error("result.json success must be a boolean")
  solution_available = get(data, "solution_available", false)
  solution_available isa Bool || error("result.json solution_available must be a boolean")
  return _api_result(
    run_id = paths.run_id,
    schema_version = string(get(data, "schema_version", _SPARLECTRA_API_SCHEMA_VERSION)),
    status = Symbol(status),
    success = success,
    converged = _optional_bool(get(data, "converged", nothing)),
    solution_available = solution_available,
    iterations = _optional_int(get(data, "iterations", nothing)),
    final_mismatch = _optional_float(get(data, "final_mismatch", nothing)),
    reason = _optional_string(get(data, "reason", nothing)),
    message = _optional_string(get(data, "message", nothing)),
    casefile = _optional_string(get(data, "casefile", nothing)),
    config_file = _optional_string(get(data, "config_file", nothing)),
    output_dir = paths.output_dir,
    logfile = _optional_string(get(data, "logfile", nothing)),
    result_file = paths.result_file,
    artifacts = collect_sparlectra_api_artifacts(paths.output_dir),
    service_phase_timings = get(data, "service_phase_timings", Dict{String,Any}[]),
    metadata = _powerflow_lifecycle_metadata(data),
    raw_result = nothing,
  )
end
