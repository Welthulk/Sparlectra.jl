# Copyright 2023–2026 Udo Schmitz
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.

# file: src/yamlparams.jl

"""
    parse_yaml_scalar(raw::AbstractString)

Parse one scalar value from Sparlectra's small dependency-free YAML subset.

This parser is intended for example and benchmark configuration files. It is
not a general-purpose YAML implementation. Supported scalars are booleans,
`null`/`~`, integers, floating-point numbers, symbols written as `:name`,
quoted strings, unquoted strings, and one-line scalar lists such as
`[100, 300, 500]`.
"""
function parse_yaml_scalar(raw::AbstractString)
  s = strip(raw)
  isempty(s) && return ""

  if startswith(s, "[") && endswith(s, "]")
    inner = strip(s[2:(end-1)])
    isempty(inner) && return Any[]
    return [parse_yaml_scalar(part) for part in _split_yaml_inline_list(inner)]
  end

  if (startswith(s, "\"") && endswith(s, "\"")) || (startswith(s, "'") && endswith(s, "'"))
    return s[2:(end-1)]
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
  for ch in inner
    if qchar != '\0'
      if ch == qchar
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
  for ch in line
    if qchar != '\0'
      if ch == qchar
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

"""
    load_yaml_dict(path::AbstractString)::Dict{String,Any}

Load a simple YAML file into a `Dict{String,Any}`.

Only the subset used by Sparlectra examples is supported: comments beginning
with `#`, 2-space-indented nested dictionaries, scalar key-value pairs, and
one-line scalar lists. Invalid indentation and unsupported lines throw clear
errors. Use a full YAML package for general YAML documents.
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
      m = match(r"^([A-Za-z0-9_\-]+)\s*:\s*(.*)$", content)
      isnothing(m) && error("Unsupported YAML line at $(path):$(lineno): $(strip(rawline))")

      level = indent_count ÷ 2
      haskey(stack, level) || error("Invalid YAML indentation at $(path):$(lineno): missing parent for indentation level $(level)")
      parent = stack[level]
      key = String(m.captures[1])
      value_raw = String(m.captures[2])

      if isempty(value_raw)
        child = Dict{String,Any}()
        parent[key] = child
        stack[level + 1] = child
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
