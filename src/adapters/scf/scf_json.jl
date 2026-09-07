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

# file: src/adapters/scf/scf_json.jl
# purpose: deterministic pretty JSON writer and a minimal parser for the
#          Sparlectra Case Format.
#          Objects print one key per line (sorted, so a rebuild diffs
#          cleanly); a component list prints one object per line, which keeps
#          case files reviewable in Git. Floats use Julia's shortest
#          round-trip representation, so a write/read/write cycle is byte
#          identical (SCF requirement: lossless round trip).

# keys that must keep their given order instead of being sorted (the root
# object reads better in PGM order; everything else is sorted)
const _SCF_ROOT_KEY_ORDER = ("version", "type", "is_batch", "attributes", "data", "sparlectra")

function _scf_json_escape(s::AbstractString)::String
  io = IOBuffer()
  for ch in s
    if ch == '"'
      print(io, "\\\"")
    elseif ch == '\\'
      print(io, "\\\\")
    elseif ch == '\n'
      print(io, "\\n")
    elseif ch == '\r'
      print(io, "\\r")
    elseif ch == '\t'
      print(io, "\\t")
    elseif ch < ' '
      print(io, "\\u", lpad(string(UInt16(ch); base = 16), 4, '0'))
    else
      print(io, ch)
    end
  end
  return String(take!(io))
end

# one scalar; Float64 goes through print (Ryu shortest round trip), and a
# whole-number float keeps its ".0" so the type survives a round trip
function _scf_json_scalar(io::IO, value)
  if value === nothing || value === missing
    print(io, "null")
  elseif value isa Bool
    print(io, value ? "true" : "false")
  elseif value isa Integer
    print(io, value)
  elseif value isa AbstractFloat
    # normalize negative zero: -0.0 and 0.0 are numerically equal but
    # textually different, which would break the byte-identical round trip
    # JSON has no Inf/NaN. Printing null for one silently produced a file the
    # reader then rejected with a misleading message (seen on MATPOWER cases
    # whose branches carry an unlimited rating), so this is a hard error: a
    # writer that can produce a non-finite value must encode it itself, with
    # `scf_infinity_sentinel` in the namespaced block or by omitting the key
    # in the PGM dataset.
    isfinite(value) || throw(ArgumentError("SCF: cannot write the non-finite number $(value); JSON has no representation for it."))
    print(io, iszero(value) ? 0.0 : Float64(value))
  else
    print(io, '"', _scf_json_escape(string(value)), '"')
  end
  return nothing
end

"""
    scf_infinity_sentinel(value) -> value

Encode a possibly non-finite nameplate number for the namespaced
`sparlectra` block. JSON cannot carry `Inf`, and an unlimited rating is real
data (MATPOWER writes `rateA = 0` for "no limit", which arrives as `Inf`), so
it travels as the string `"inf"` / `"-inf"` and comes back as a number on
read. Finite values pass through unchanged. The PGM `data` section never uses
this: a value it cannot express is omitted there.
"""
scf_infinity_sentinel(value) = value isa AbstractFloat && !isfinite(value) ? (isnan(value) ? throw(ArgumentError("SCF: NaN is not a value the format can carry.")) : (value > 0 ? "inf" : "-inf")) : value

"""
    scf_number_or_sentinel(value, context) -> Float64

Read back what [`scf_infinity_sentinel`](@ref) wrote.
"""
function scf_number_or_sentinel(value, context::AbstractString)::Float64
  value == "inf" && return Inf
  value == "-inf" && return -Inf
  value isa Number && return Float64(value)
  throw(ArgumentError("SCF: $(context) must be a number or \"inf\"/\"-inf\", got $(repr(value))."))
end

_scf_is_scalar(v) = v === nothing || v === missing || v isa Bool || v isa Integer || v isa AbstractFloat || v isa AbstractString || v isa Symbol

# name the key when a non-finite value reaches the writer: "cannot write Inf"
# alone leaves the producer unknown
function _scf_json_check_finite(key, v)
  v isa AbstractFloat && !isfinite(v) && throw(ArgumentError("SCF: key \"$(key)\" carries the non-finite number $(v); JSON has no representation for it (use scf_infinity_sentinel in the namespaced block, or omit the key in the PGM dataset)."))
  return nothing
end

function _scf_ordered_keys(d::AbstractDict, key_order)
  ks = String[String(k) for k in keys(d)]
  key_order === nothing && return sort!(ks)
  # listed keys first in the given order, the rest sorted behind them
  head = String[k for k in key_order if k in ks]
  tail = sort!(String[k for k in ks if !(k in key_order)])
  return vcat(head, tail)
end

"""
    _scf_json_write(io, value; indent = 0, key_order = nothing)

Write `value` as pretty JSON: objects and the objects inside a list break
one key per line. Only a vector of plain scalars stays inline (`[1, 2, 3]`),
because breaking a list of node ids over twenty lines helps nobody.

Component rows used to be written as ONE compact line each, which kept
diffs narrow but read like a stream of data rather than a document. The
expanded form is what JSON tooling produces and what a reader expects; a
diff now points at the changed FIELD instead of the changed component.
"""
function _scf_json_write(io::IO, value; indent::Int = 0, key_order = nothing)
  pad = "  "^indent
  if value isa AbstractDict
    if isempty(value)
      print(io, "{}")
      return nothing
    end
    println(io, '{')
    ks = _scf_ordered_keys(value, key_order)
    for (i, key) in enumerate(ks)
      print(io, pad, "  ", '"', _scf_json_escape(String(key)), "\": ")
      _scf_json_check_finite(key, value[key])
      _scf_json_write(io, value[key]; indent = indent + 1)
      println(io, i == length(ks) ? "" : ",")
    end
    print(io, pad, '}')
  elseif value isa AbstractVector
    if isempty(value)
      print(io, "[]")
      return nothing
    end
    if all(_scf_is_scalar, value)
      print(io, '[')
      for (i, item) in enumerate(value)
        i == 1 || print(io, ", ")
        _scf_json_scalar(io, item)
      end
      print(io, ']')
      return nothing
    end
    println(io, '[')
    for (i, item) in enumerate(value)
      print(io, pad, "  ")
      _scf_json_write(io, item; indent = indent + 1)
      println(io, i == length(value) ? "" : ",")
    end
    print(io, pad, ']')
  else
    _scf_json_scalar(io, value)
  end
  return nothing
end

"""
    scf_json_string(root::AbstractDict) -> String

Serialize an SCF root object to its canonical text form (trailing newline
included). Deterministic: the same object always produces the same bytes.
"""
function scf_json_string(root::AbstractDict)::String
  io = IOBuffer()
  _scf_json_write(io, root; key_order = _SCF_ROOT_KEY_ORDER)
  println(io)
  return String(take!(io))
end

## ---------------------------------------------------------------------------
## Reader: a minimal, dependency-free JSON parser in the spirit of the
## repository's minimal YAML reader. Integers stay `Int` (component ids must
## not silently become floats), everything with a fraction or exponent
## becomes `Float64`. Errors name the byte position, so a broken case file
## points at itself instead of at a generic parse failure.
## ---------------------------------------------------------------------------

mutable struct _ScfJsonCursor
  text::String
  pos::Int
end

_scf_json_error(c::_ScfJsonCursor, msg::AbstractString) = throw(ArgumentError(string("SCF JSON parse error at byte ", c.pos, ": ", msg)))

function _scf_skip_ws!(c::_ScfJsonCursor)
  n = ncodeunits(c.text)
  while c.pos <= n
    ch = c.text[c.pos]
    (ch == ' ' || ch == '\n' || ch == '\r' || ch == '\t') || break
    c.pos = nextind(c.text, c.pos)
  end
  return nothing
end

function _scf_expect!(c::_ScfJsonCursor, ch::Char)
  _scf_skip_ws!(c)
  (c.pos <= ncodeunits(c.text) && c.text[c.pos] == ch) || _scf_json_error(c, string("expected '", ch, "'"))
  c.pos = nextind(c.text, c.pos)
  return nothing
end

function _scf_parse_string!(c::_ScfJsonCursor)::String
  _scf_expect!(c, '"')
  io = IOBuffer()
  n = ncodeunits(c.text)
  while true
    c.pos <= n || _scf_json_error(c, "unterminated string")
    ch = c.text[c.pos]
    c.pos = nextind(c.text, c.pos)
    if ch == '"'
      return String(take!(io))
    elseif ch == '\\'
      c.pos <= n || _scf_json_error(c, "unterminated escape")
      esc = c.text[c.pos]
      c.pos = nextind(c.text, c.pos)
      if esc == 'n'
        print(io, '\n')
      elseif esc == 't'
        print(io, '\t')
      elseif esc == 'r'
        print(io, '\r')
      elseif esc == 'b'
        print(io, '\b')
      elseif esc == 'f'
        print(io, '\f')
      elseif esc == 'u'
        hex = c.text[c.pos:(c.pos+3)]
        c.pos += 4
        print(io, Char(parse(UInt16, hex; base = 16)))
      else
        print(io, esc)   # covers '"', '\\', '/'
      end
    else
      print(io, ch)
    end
  end
end

function _scf_parse_number!(c::_ScfJsonCursor)
  start = c.pos
  n = ncodeunits(c.text)
  isfloat = false
  while c.pos <= n
    ch = c.text[c.pos]
    if ch == '-' || ch == '+' || isdigit(ch)
      c.pos = nextind(c.text, c.pos)
    elseif ch == '.' || ch == 'e' || ch == 'E'
      isfloat = true
      c.pos = nextind(c.text, c.pos)
    else
      break
    end
  end
  lit = c.text[start:prevind(c.text, c.pos)]
  if isfloat
    v = tryparse(Float64, lit)
    v === nothing && _scf_json_error(c, string("invalid number ", repr(lit)))
    return v
  end
  v = tryparse(Int, lit)
  v === nothing && _scf_json_error(c, string("invalid integer ", repr(lit)))
  return v
end

function _scf_parse_value!(c::_ScfJsonCursor)
  _scf_skip_ws!(c)
  c.pos <= ncodeunits(c.text) || _scf_json_error(c, "unexpected end of input")
  ch = c.text[c.pos]
  if ch == '{'
    c.pos = nextind(c.text, c.pos)
    out = Dict{String,Any}()
    _scf_skip_ws!(c)
    if c.pos <= ncodeunits(c.text) && c.text[c.pos] == '}'
      c.pos = nextind(c.text, c.pos)
      return out
    end
    while true
      _scf_skip_ws!(c)
      key = _scf_parse_string!(c)
      _scf_expect!(c, ':')
      out[key] = _scf_parse_value!(c)
      _scf_skip_ws!(c)
      c.pos <= ncodeunits(c.text) || _scf_json_error(c, "unterminated object")
      sep = c.text[c.pos]
      c.pos = nextind(c.text, c.pos)
      sep == ',' && continue
      sep == '}' && return out
      _scf_json_error(c, "expected ',' or '}' in object")
    end
  elseif ch == '['
    c.pos = nextind(c.text, c.pos)
    out = Vector{Any}()
    _scf_skip_ws!(c)
    if c.pos <= ncodeunits(c.text) && c.text[c.pos] == ']'
      c.pos = nextind(c.text, c.pos)
      return out
    end
    while true
      push!(out, _scf_parse_value!(c))
      _scf_skip_ws!(c)
      c.pos <= ncodeunits(c.text) || _scf_json_error(c, "unterminated array")
      sep = c.text[c.pos]
      c.pos = nextind(c.text, c.pos)
      sep == ',' && continue
      sep == ']' && return out
      _scf_json_error(c, "expected ',' or ']' in array")
    end
  elseif ch == '"'
    return _scf_parse_string!(c)
  elseif startswith(SubString(c.text, c.pos), "true")
    c.pos += 4
    return true
  elseif startswith(SubString(c.text, c.pos), "false")
    c.pos += 5
    return false
  elseif startswith(SubString(c.text, c.pos), "null")
    c.pos += 4
    return nothing
  else
    return _scf_parse_number!(c)
  end
end

"""
    scf_json_parse(text) -> Dict{String,Any}

Parse an SCF case file. Object keys become `String`, integers stay `Int`
(component ids), numbers with a fraction or exponent become `Float64`.
"""
function scf_json_parse(text::AbstractString)::Dict{String,Any}
  c = _ScfJsonCursor(String(text), 1)
  value = _scf_parse_value!(c)
  _scf_skip_ws!(c)
  c.pos <= ncodeunits(c.text) && _scf_json_error(c, "trailing content after the root object")
  value isa AbstractDict || throw(ArgumentError("SCF: the root of a case file must be a JSON object"))
  return Dict{String,Any}(value)
end
