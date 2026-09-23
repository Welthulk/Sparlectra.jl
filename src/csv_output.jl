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
# file: src/csv_output.jl
# purpose: CSV formatting and writing shared by every result artifact:
#          the format table (technical, excel_de, excel_us), number and
#          cell formatting, the streaming and buffered NamedTuple writers.
#          Library code (island diagnostics, contingency and scenario
#          results, state estimation) writes CSV through this file, so it
#          lives in the core and not in the service layer.

function _resolve_detailed_csv_format(value)::NamedTuple
  name = String(value)
  name == "technical" && return (name = name, delimiter = ',', decimal_separator = '.', thousands_separator = "")
  name == "excel_de" && return (name = name, delimiter = ';', decimal_separator = ',', thousands_separator = ".")
  name == "excel_us" && return (name = name, delimiter = ',', decimal_separator = '.', thousands_separator = ",")
  throw(ArgumentError("Unsupported detailed_result_csv_format \"$(name)\". Expected technical, excel_de, or excel_us."))
end

function _group_csv_integer(text::AbstractString, separator::AbstractString)::String
  isempty(separator) && return String(text)

  sign = startswith(text, "-") ? "-" : ""
  digits = isempty(sign) ? String(text) : text[2:end]

  # Thousands grouping is only valid for integer digit strings. Leave textual
  # values such as Bool strings unchanged if they accidentally reach this helper.
  (isempty(digits) || !all(isdigit, digits)) && return String(text)

  first_group = mod(length(digits), 3)
  first_group == 0 && (first_group = 3)
  groups = String[digits[1:first_group]]
  for start = (first_group+1):3:length(digits)
    push!(groups, digits[start:(start+2)])
  end
  return sign * join(groups, separator)
end

function _format_csv_number(value::Integer, format)::String
  return _group_csv_integer(string(value), format.thousands_separator)
end

function _format_csv_number(value::AbstractFloat, format)::String
  isnan(value) && return "NaN"
  isinf(value) && return signbit(value) ? "-Inf" : "Inf"
  technical = @sprintf("%.15g", value)
  if format.name != "technical"
    exponent_marker = findfirst(character -> character in ('e', 'E'), technical)
    if exponent_marker !== nothing
      mantissa = technical[begin:prevind(technical, exponent_marker)]
      exponent = parse(Int, technical[nextind(technical, exponent_marker):end])
      sign = startswith(mantissa, "-") ? "-" : ""
      unsigned = isempty(sign) ? mantissa : mantissa[2:end]
      dot_index = findfirst(==('.'), unsigned)
      fractional_digits = dot_index === nothing ? 0 : ncodeunits(unsigned) - dot_index
      digits = replace(unsigned, "." => "")
      decimal_position = ncodeunits(digits) - fractional_digits + exponent
      if decimal_position <= 0
        technical = sign * "0." * repeat("0", -decimal_position) * digits
      elseif decimal_position >= ncodeunits(digits)
        technical = sign * digits * repeat("0", decimal_position - ncodeunits(digits))
      else
        technical = sign * digits[1:decimal_position] * "." * digits[(decimal_position+1):end]
      end
    end
    dot_index = findfirst(==('.'), technical)
    if dot_index !== nothing
      last_nonzero = lastindex(technical)
      while last_nonzero > dot_index && technical[last_nonzero] == '0'
        last_nonzero = prevind(technical, last_nonzero)
      end
      technical = last_nonzero == dot_index ? technical[begin:prevind(technical, dot_index)] : technical[begin:last_nonzero]
      isempty(technical) && (technical = "0")
      technical == "-0" && (technical = "0")
    end
  end
  exponent_marker = format.name == "technical" ? findfirst(character -> character in ('e', 'E'), technical) : nothing
  mantissa = exponent_marker === nothing ? technical : technical[begin:prevind(technical, exponent_marker)]
  exponent = exponent_marker === nothing ? "" : technical[exponent_marker:end]
  dot_index = findfirst(==('.'), mantissa)
  integer_text = dot_index === nothing ? mantissa : mantissa[begin:prevind(mantissa, dot_index)]
  integer_part = _group_csv_integer(integer_text, format.thousands_separator)
  fractional_part = dot_index === nothing ? "" : string(format.decimal_separator, mantissa[nextind(mantissa, dot_index):end])
  return integer_part * fractional_part * exponent
end

struct CsvFormatRuntime
  name::String
  delimiter::Char
  decimal_separator::Char
  thousands_separator::String
end

CsvFormatRuntime(format) = CsvFormatRuntime(String(format.name), format.delimiter, format.decimal_separator, String(format.thousands_separator))

function _csv_needs_quotes(text::AbstractString, delimiter::Char)::Bool
  for character in text
    character in (delimiter, '"', '\r', '\n') && return true
  end
  return false
end

function write_csv_cell!(io::IO, value, delimiter::Char, fmt::CsvFormatRuntime)
  value === missing && return nothing
  value === nothing && return nothing
  text = value isa Bool ? string(value) : value isa Integer ? _format_csv_number(value, fmt) : value isa AbstractFloat ? _format_csv_number(value, fmt) : string(value)
  if _csv_needs_quotes(text, delimiter)
    print(io, '"')
    for character in text
      character == '"' && print(io, '"')
      print(io, character)
    end
    print(io, '"')
  else
    print(io, text)
  end
  return nothing
end

function write_csv_row_direct!(io::IO, delimiter::Char, fmt::CsvFormatRuntime, values...)
  first = true
  for value in values
    first || print(io, delimiter)
    write_csv_cell!(io, value, delimiter, fmt)
    first = false
  end
  println(io)
  return nothing
end

function _format_csv_value(value, format)::String
  value === missing && return ""
  value === nothing && return ""

  # Bool is an Integer subtype in Julia. Format it before Integer values so
  # Excel-oriented thousands grouping cannot turn true/false into t.rue/fa.lse.
  value isa Bool && return string(value)

  value isa Integer && return _format_csv_number(value, format)
  value isa AbstractFloat && return _format_csv_number(value, format)
  return string(value)
end

function _csv_field(value, delimiter::Char, format = _resolve_detailed_csv_format("technical"))::String
  value === missing && return ""
  value === nothing && return ""
  text = _format_csv_value(value, format)
  if any(character -> character in (delimiter, '"', '\r', '\n'), text)
    return "\"" * replace(text, "\"" => "\"\"") * "\""
  end
  return text
end

const _DETAILED_CSV_BUFFER_INITIAL_BYTES_DEFAULT = 8 * 1024 * 1024
const _DETAILED_CSV_BUFFER_MAX_BYTES_DEFAULT = 64 * 1024 * 1024
const _DETAILED_CSV_STREAMING_THRESHOLD_ROWS_DEFAULT = 100_000
const _DETAILED_CSV_DIRECT_THRESHOLD_BUSES_DEFAULT = 10_000

function _csv_write_options(config = nothing)::NamedTuple
  output = config isa SparlectraConfig ? config.output : config isa OutputConfig ? config : nothing
  mode = output === nothing ? :auto : output.detailed_result_csv_write_mode
  mode in OUTPUT_DETAILED_RESULT_CSV_WRITE_MODE_VALUES || throw(ArgumentError("Unsupported detailed_result_csv_write_mode \"$(mode)\". Expected auto, buffered, or streaming."))
  initial_bytes = output === nothing ? _DETAILED_CSV_BUFFER_INITIAL_BYTES_DEFAULT : output.detailed_result_csv_buffer_initial_bytes
  max_bytes = output === nothing ? _DETAILED_CSV_BUFFER_MAX_BYTES_DEFAULT : output.detailed_result_csv_buffer_max_bytes
  threshold_rows = output === nothing ? _DETAILED_CSV_STREAMING_THRESHOLD_ROWS_DEFAULT : output.detailed_result_csv_streaming_threshold_rows
  initial_bytes = initial_bytes < 0 ? _DETAILED_CSV_BUFFER_INITIAL_BYTES_DEFAULT : initial_bytes
  max_bytes = max_bytes <= 0 ? _DETAILED_CSV_BUFFER_MAX_BYTES_DEFAULT : max_bytes
  threshold_rows = threshold_rows <= 0 ? _DETAILED_CSV_STREAMING_THRESHOLD_ROWS_DEFAULT : threshold_rows
  return (mode = mode, initial_bytes = initial_bytes, max_bytes = max_bytes, threshold_rows = threshold_rows)
end

function _estimated_namedtuple_csv_bytes(rows::AbstractVector, columns)::Int
  return length(join(String.(columns), ',')) + 1 + length(rows) * max(1, length(columns)) * 24
end

function _write_namedtuple_csv_row(io::IO, row, columns, delimiter::Char, resolved_format)
  println(io, join((_csv_field(getproperty(row, column), delimiter, resolved_format) for column in columns), delimiter))
  return nothing
end


function _write_namedtuple_csv_buffered(path::AbstractString, rows::AbstractVector, columns, delimiter::Char, resolved_format; initial_bytes::Integer = _DETAILED_CSV_BUFFER_INITIAL_BYTES_DEFAULT)
  buffer = IOBuffer(; sizehint = max(0, Int(initial_bytes)))
  println(buffer, join(String.(columns), delimiter))
  for row in rows
    _write_namedtuple_csv_row(buffer, row, columns, delimiter, resolved_format)
  end
  open(path, "w") do io
    write(io, take!(buffer))
  end
  return path
end

function _write_namedtuple_csv_streaming(path::AbstractString, rows::AbstractVector, columns, delimiter::Char, resolved_format)
  open(path, "w") do io
    println(io, join(String.(columns), delimiter))
    for row in rows
      _write_namedtuple_csv_row(io, row, columns, delimiter, resolved_format)
    end
  end
  return path
end

function _select_namedtuple_csv_write_mode(rows::AbstractVector, columns; config = nothing, estimated_rows::Integer = length(rows))::Symbol
  options = _csv_write_options(config)
  options.mode === :buffered && return :buffered
  options.mode === :streaming && return :streaming
  estimated_rows > options.threshold_rows && return :streaming
  _estimated_namedtuple_csv_bytes(rows, columns) > options.max_bytes && return :streaming
  return :buffered
end

function _write_namedtuple_csv(path::AbstractString, rows::AbstractVector, columns; delimiter::Union{Nothing,Char} = nothing, format = nothing, config = nothing, estimated_rows::Integer = length(rows))
  # `delimiter` only needs to be given when it must override the format's own
  # delimiter (or when there is no format at all); every current call site
  # that passes a resolved `format` name relies on this default, so it must
  # never silently disagree with that format's delimiter (issue #376: this
  # threw for every caller that did not also repeat `delimiter` explicitly
  # once formats other than "technical" started reaching them).
  resolved_format = format === nothing ? (name = "custom", delimiter = delimiter === nothing ? ',' : delimiter, decimal_separator = '.', thousands_separator = "") : _resolve_detailed_csv_format(format)
  delimiter = delimiter === nothing ? resolved_format.delimiter : delimiter
  delimiter in (',', ';') || throw(ArgumentError("CSV delimiter must be ',' or ';'."))
  resolved_format.delimiter == delimiter || throw(ArgumentError("CSV delimiter does not match detailed CSV format $(resolved_format.name)."))
  options = _csv_write_options(config)
  mode = _select_namedtuple_csv_write_mode(rows, columns; config, estimated_rows)
  mode === :streaming && return _write_namedtuple_csv_streaming(path, rows, columns, delimiter, resolved_format)
  return _write_namedtuple_csv_buffered(path, rows, columns, delimiter, resolved_format; initial_bytes = min(options.initial_bytes, options.max_bytes))
end
