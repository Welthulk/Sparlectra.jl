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

# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 10.5.2026
# file: src/config/yamlparams.jl
# purpose: dependency-free YAML subset parser for example and benchmark
#          configuration files (parse_yaml_scalar, load_yaml_dict,
#          merge_yaml_dict!)

"""
    parse_yaml_scalar(raw::AbstractString)

Parse one scalar value from Sparlectra's small dependency-free YAML subset.

This parser is intended for example and benchmark configuration files. It is
not a general-purpose YAML implementation. Supported scalars are booleans,
`null`/`~`, integers, floating-point numbers, symbols written as `:name`,
quoted strings, unquoted strings, and one-line scalar lists such as
`[100, 300, 500]`.

A double-quoted string resolves the two escapes the writers produce, `\\\\`
(one backslash) and `\\"` (a quote); every other backslash stays as it is,
so a hand-written `"C:\\data"` or `"\\d{1,2}"` reads unchanged. In a
single-quoted string `''` is one quote, backslashes are literal.
"""
function parse_yaml_scalar(raw::AbstractString)
  s = strip(raw)
  isempty(s) && return ""

  if startswith(s, "[") && endswith(s, "]")
    inner = strip(s[2:(end-1)])
    isempty(inner) && return Any[]
    return [parse_yaml_scalar(part) for part in _split_yaml_inline_list(inner)]
  end

  if length(s) >= 2 && startswith(s, "\"") && endswith(s, "\"")
    return _yaml_unescape_double(s[2:(end-1)])
  end
  if length(s) >= 2 && startswith(s, "'") && endswith(s, "'")
    return replace(s[2:(end-1)], "''" => "'")
  end

  lower = lowercase(s)
  if lower in ("true", "yes", "on")
    return true
  elseif lower in ("false", "no", "off")
    return false
  elseif lower in ("null", "~")
    return nothing
  elseif startswith(s, ":") && length(s) > 1
    return Symbol(s[2:end])
  end

  if occursin(r"^[+-]?\d+$", s)
    return parse(Int, s)
  end
  if occursin(r"^[+-]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][+-]?\d+)?$", s)
    return parse(Float64, s)
  end

  return s
end

function _split_yaml_inline_list(inner::AbstractString)::Vector{String}
  parts = String[]
  buf = IOBuffer()
  qchar = '\0'
  escaped = false
  for ch in inner
    if qchar != '\0'
      # inside double quotes a backslash escapes the next character, so
      # \" does not end the string
      if escaped
        escaped = false
      elseif qchar == '"' && ch == '\\'
        escaped = true
      elseif ch == qchar
        qchar = '\0'
      end
      print(buf, ch)
    elseif ch == '\'' || ch == '"'
      qchar = ch
      print(buf, ch)
    elseif ch == ','
      push!(parts, strip(String(take!(buf))))
    else
      print(buf, ch)
    end
  end
  qchar == '\0' || error("Unterminated quote in YAML inline list: [$(inner)]")
  push!(parts, strip(String(take!(buf))))
  return parts
end

function _strip_yaml_comment(line::AbstractString)::String
  buf = IOBuffer()
  qchar = '\0'
  escaped = false
  for ch in line
    if qchar != '\0'
      # inside double quotes a backslash escapes the next character, so
      # \" does not end the string and a later # stays text
      if escaped
        escaped = false
      elseif qchar == '"' && ch == '\\'
        escaped = true
      elseif ch == qchar
        qchar = '\0'
      end
      print(buf, ch)
    elseif ch == '\'' || ch == '"'
      qchar = ch
      print(buf, ch)
    elseif ch == '#'
      break
    else
      print(buf, ch)
    end
  end
  return String(take!(buf))
end

# The inverse of _yaml_quoted: in a double-quoted scalar `\\` is one
# backslash and `\"` a quote; any other backslash is kept with the character
# after it, so hand-written Windows paths and regular expressions read as
# written (and a file saved before 0.40.0 with doubled backslashes reads back
# with single ones).
function _yaml_unescape_double(body::AbstractString)::String
  occursin('\\', body) || return String(body)
  out = IOBuffer()
  escaped = false
  for ch in body
    if escaped
      (ch == '\\' || ch == '"') || print(out, '\\')
      print(out, ch)
      escaped = false
    elseif ch == '\\'
      escaped = true
    else
      print(out, ch)
    end
  end
  escaped && print(out, '\\')
  return String(take!(out))
end

"""
    yaml_quoted(text) -> String

The one quoting rule of Sparlectra's YAML writers: `text` in double quotes
with every backslash written as `\\\\` and every quote as `\\"`, the two
escapes [`parse_yaml_scalar`](@ref) resolves, so writing and reading give the
same string for every value. A line break cannot be written in this
line-based subset and throws an `ArgumentError`.
"""
function yaml_quoted(text::AbstractString)::String
  (occursin('\n', text) || occursin('\r', text)) && throw(ArgumentError("YAML subset: a string value cannot contain a line break: $(repr(String(text)))"))
  return string('"', replace(String(text), "\\" => "\\\\", "\"" => "\\\""), '"')
end

# A string is written bare only when the reader gives it back unchanged as a
# string; everything else is quoted. This is the self-inverse criterion, not a
# maintained list: it covers the YAML 1.1 boolean family (live case: the enum
# value "off" was written bare and came back as `false`), null spellings,
# every number form the reader knows including leading-zero integers,
# ":symbol" spellings, inline-list lookalikes and padded values, and it stays
# correct when parse_yaml_scalar learns new forms.
function _yaml_ambiguous_scalar(s::AbstractString)
  v = parse_yaml_scalar(String(s))
  return !(v isa AbstractString && String(v) == String(s))
end

# a string scalar as the writers emit it: bare where that reads back
# unchanged, quoted (yaml_quoted) where a character would end the value early
# (`#` starts a comment, `:` a mapping), the reader would take it for another
# type, or quotes and backslashes would be read differently
function yaml_string_text(text::AbstractString)::String
  s = String(text)
  (isempty(s) || occursin(r"[:#\[\]{},&*!|>'\"%@`\\]|^[-?]|^\s|\s$", s) || _yaml_ambiguous_scalar(s)) && return yaml_quoted(s)
  return s
end

# A mapping key is bare (letters, digits, "_" and "-", or a decimal number
# such as 10.5) or in single or double quotes, which are stripped; the text
# after the colon is the value. The key stays text either way, the reader
# does not type it (the sigma table reads its row keys as kV itself). A
# dotted name such as power_flow.max_iter is no key: the line is refused
# instead of being read as an unknown key. Shared with the line scanner
# _yaml_key_lines of configuration.jl, so both see the same keys.
const _YAML_MAPPING_KEY_RE = r"^(?:\"([^\"]+)\"|'([^']+)'|([A-Za-z0-9_\-]+|[0-9]+\.[0-9]+))\s*:\s*(.*)$"

# the bare key form of _YAML_MAPPING_KEY_RE; the writers quote every other key
const _YAML_BARE_KEY_RE = r"^(?:[A-Za-z0-9_\-]+|[0-9]+\.[0-9]+)$"

# (key, value text) of a "key: value" line without its indentation, or
# nothing when the line is no mapping line
function _yaml_mapping_key(content::AbstractString)::Union{Nothing,Tuple{String,String}}
  m = match(_YAML_MAPPING_KEY_RE, content)
  m === nothing && return nothing
  return (String(something(m.captures[1], m.captures[2], m.captures[3])), String(m.captures[4]))
end

# A mapping key as both writers emit it: bare when it has the bare form,
# otherwise quoted; the reader strips the quotes without unescaping. Double
# quotes when the key has no double quote and does not end in a backslash
# (inside double quotes the comment scanner reads a backslash as escaping
# the next character, so a final one would escape the closing quote), else
# single quotes, where a backslash is plain text. Refused: a key the reader
# never gives back (empty, a line break, both quote kinds), and a single
# quote with a final backslash, which no quote style keeps safe from the
# comment scanner when the value holds a "#".
function _yaml_key_text(key)::String
  s = string(key)
  occursin(_YAML_BARE_KEY_RE, s) && return s
  if !isempty(s) && !occursin(r"[\r\n]", s)
    (occursin('"', s) || endswith(s, '\\')) || return string('"', s, '"')
    occursin('\'', s) || return string('\'', s, '\'')
  end
  throw(ArgumentError("YAML key $(repr(s)) cannot be written: the key is empty, contains a line break, or contains a single quote together with a double quote or a final backslash"))
end

"""
    load_yaml_dict(path::AbstractString)::Dict{String,Any}

Load a simple YAML file into a `Dict{String,Any}`.

Only the subset used by Sparlectra examples is supported: comments beginning
with `#`, 2-space-indented nested dictionaries, scalar key-value pairs, and
one-line scalar lists. A mapping key is bare (`380:`, `10.5:`, `max_iter:`)
or quoted (`"10.5":`, `'380':`, the quotes stripped) and is always read as
text; a dotted name (`power_flow.max_iter:`) is no key and throws, and a
key at the top level must not contain a dot even in quotes. Invalid indentation and unsupported lines throw clear errors. Use a
full YAML package for general YAML documents.
"""
function load_yaml_dict(path::AbstractString)::Dict{String,Any}
  root = Dict{String,Any}()
  stack = Dict{Int,Dict{String,Any}}(0 => root)

  open(path, "r") do io
    for (lineno, rawline) in enumerate(eachline(io))
      line = rstrip(_strip_yaml_comment(rawline))
      isempty(strip(line)) && continue

      indent_count = 0
      for ch in line
        ch == ' ' || break
        indent_count += 1
      end
      indent_count % 2 == 0 || error("Invalid YAML indentation at $(path):$(lineno): indentation must use multiples of 2 spaces")

      content = strip(line)
      m = _yaml_mapping_key(content)
      # a list item or a flow mapping is the usual form of a table in YAML
      # and the one most often copied in; name the block form instead of
      # the bare "unsupported" (the subset itself stays as it is)
      if isnothing(m) && (startswith(content, "- ") || content == "-" || startswith(content, "{"))
        error("Unsupported YAML line at $(path):$(lineno): $(strip(rawline)): list items (\"- ...\") and flow mappings (\"{...}\") are not part of Sparlectra's YAML subset. Write the entries as named mappings in block style, one key per line (see state_estimation.sigma_by_voltage_level and control.controllers in configuration.yaml.example).")
      end
      isnothing(m) && error("Unsupported YAML line at $(path):$(lineno): $(strip(rawline))")

      level = indent_count ÷ 2
      haskey(stack, level) || error("Invalid YAML indentation at $(path):$(lineno): missing parent for indentation level $(level)")
      parent = stack[level]
      key, value_raw = m
      # a top-level key is a section name: a dot there, also in quotes
      # ("power_flow.max_iter"), is a configuration path written flat, and
      # read as a key it would be dropped as unknown or applied past the
      # checks that look at the nested form (sigma table rows)
      level == 0 && occursin('.', key) && error("Unsupported YAML line at $(path):$(lineno): $(strip(rawline)): a top-level key must not contain a dot; write a configuration path as nested keys")

      if isempty(value_raw)
        child = Dict{String,Any}()
        parent[key] = child
        stack[level+1] = child
        for k in collect(keys(stack))
          if k > level + 1
            delete!(stack, k)
          end
        end
      else
        parent[key] = parse_yaml_scalar(value_raw)
        for k in collect(keys(stack))
          if k > level
            delete!(stack, k)
          end
        end
      end
    end
  end

  return root
end

"""
    merge_yaml_dict!(dst::Dict{String,Any}, src::Dict{String,Any})

Recursively merge YAML dictionaries from `src` into `dst` and return `dst`.
Nested dictionaries are merged; non-dictionary values replace existing values.
"""
function merge_yaml_dict!(dst::Dict{String,Any}, src::Dict{String,Any})
  for (key, value) in src
    if haskey(dst, key) && dst[key] isa Dict{String,Any} && value isa Dict{String,Any}
      merge_yaml_dict!(dst[key], value)
    else
      dst[key] = value
    end
  end
  return dst
end

"""
    as_bool(x)::Bool

Convert YAML booleans or boolean-like strings (`true`, `false`, `yes`, `no`,
`on`, `off`) to `Bool`. Invalid values throw an `ArgumentError`.
"""
function as_bool(x)::Bool
  x isa Bool && return x
  if x isa AbstractString
    lower = lowercase(strip(x))
    lower in ("true", "yes", "on") && return true
    lower in ("false", "no", "off") && return false
  end
  throw(ArgumentError("Cannot convert $(repr(x)) to Bool; expected true/false, yes/no, or on/off."))
end

"""
    as_int_vector(x)::Vector{Int}

Convert an integer or a one-dimensional vector of integer values to
`Vector{Int}`. Invalid values throw an `ArgumentError`.
"""
function as_int_vector(x)::Vector{Int}
  if x isa Integer
    return [Int(x)]
  elseif x isa AbstractVector
    values = Int[]
    for item in x
      if item isa Integer
        push!(values, Int(item))
      else
        throw(ArgumentError("Cannot convert $(repr(x)) to Vector{Int}; every element must be an integer."))
      end
    end
    return values
  end
  throw(ArgumentError("Cannot convert $(repr(x)) to Vector{Int}; expected an integer or vector of integers."))
end
