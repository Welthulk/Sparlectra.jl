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

# file: src/adapters/powsybl/powsybl_csv.jl
# purpose: read and write a PowSyBl table bundle (a directory with
#          manifest.json and one CSV per table) on DelimitedFiles and the
#          JSON reader of the SCF adapter. The conventions both sides keep:
#          comma separated, one header line, strings always double-quoted
#          with embedded quotes doubled, floats in full precision with NaN,
#          Inf and -Inf as literals, booleans true/false, integers plain.
#          The reader also accepts True/False and an empty float cell (as
#          NaN), so a hand-edited fixture still loads, and an empty tap
#          cell (the changer at its middle position); the writer never
#          produces any of these.

const POWSYBL_BUNDLE_FORMAT = "powsybl_tables"
const POWSYBL_BUNDLE_FORMAT_VERSION = 1
const POWSYBL_MANIFEST_FILE = "manifest.json"

_powsybl_bundle_error(dir::AbstractString, what::AbstractString) = ArgumentError("powsybl bundle $(dir): $(what)")

function _powsybl_parse_float(text::AbstractString, dir, table, column, row)::Float64
  s = strip(text)
  isempty(s) && return NaN
  v = tryparse(Float64, s)
  v === nothing && throw(_powsybl_bundle_error(dir, "table $(table) column $(column) row $(row): $(repr(s)) is not a number"))
  return v
end

function _powsybl_parse_int(text::AbstractString, dir, table, column, row)::Int
  s = strip(text)
  v = tryparse(Int, s)
  v === nothing || return v
  # pandas writes an integer column as floats once it holds a NaN; an
  # integral float is accepted, a NaN in an integer column is not
  f = tryparse(Float64, s)
  (f === nothing || !isfinite(f) || f != round(f)) && throw(_powsybl_bundle_error(dir, "table $(table) column $(column) row $(row): $(repr(s)) is not an integer"))
  return Int(f)
end

function _powsybl_parse_bool(text::AbstractString, dir, table, column, row)::Bool
  s = lowercase(strip(text))
  s in ("true", "1") && return true
  s in ("false", "0") && return false
  throw(_powsybl_bundle_error(dir, "table $(table) column $(column) row $(row): $(repr(strip(text))) is not a boolean"))
end

# Tap position cells that may be empty: the changer then sits at its
# middle position (_powsybl_bundle_tap_middle!). An empty cell reads as
# _POWSYBL_NO_TAP until that step replaces it.
const _POWSYBL_NO_TAP = typemin(Int)
const _POWSYBL_TAP_CELLS = Dict{String,Vector{Symbol}}(
  "ratio_tap_changers" => [:tap],
  "phase_tap_changers" => [:tap],
  "3_windings_transformers" => [Symbol(p, "_tap_position", k) for p in ("ratio", "phase") for k in 1:3],
)

function _powsybl_convert_column(cells::AbstractVector{<:AbstractString}, T::DataType, dir, table, column)
  T === String && return String[String(c) for c in cells]
  T === Float64 && return Float64[_powsybl_parse_float(c, dir, table, column, i) for (i, c) in enumerate(cells)]
  if T === Int && column in get(_POWSYBL_TAP_CELLS, String(table), Symbol[])
    return Int[isempty(strip(c)) ? _POWSYBL_NO_TAP : _powsybl_parse_int(c, dir, table, column, i) for (i, c) in enumerate(cells)]
  end
  T === Int && return Int[_powsybl_parse_int(c, dir, table, column, i) for (i, c) in enumerate(cells)]
  return Bool[_powsybl_parse_bool(c, dir, table, column, i) for (i, c) in enumerate(cells)]
end

# The manifest's format_version arrives as whatever number the JSON
# reader produced; both 1 and 1.0 are version 1.
_powsybl_manifest_version(manifest::AbstractDict) = (v = get(manifest, "format_version", nothing); v isa Number && isfinite(v) && v == round(v) ? Int(v) : nothing)

"""
    read_powsybl_manifest(dir) -> Dict{String,Any}

The manifest of a bundle directory, checked for the format marker and the
format version. The first step of [`read_powsybl_bundle`](@ref), also used
by format detection.
"""
function read_powsybl_manifest(dir::AbstractString)::Dict{String,Any}
  path = joinpath(dir, POWSYBL_MANIFEST_FILE)
  isfile(path) || throw(_powsybl_bundle_error(dir, "no $(POWSYBL_MANIFEST_FILE)"))
  manifest = scf_json_parse(read(path, String))
  fmt = get(manifest, "format", nothing)
  fmt == POWSYBL_BUNDLE_FORMAT || throw(_powsybl_bundle_error(dir, "manifest format $(repr(fmt)) is not $(repr(POWSYBL_BUNDLE_FORMAT))"))
  version = _powsybl_manifest_version(manifest)
  version == POWSYBL_BUNDLE_FORMAT_VERSION || throw(_powsybl_bundle_error(dir, "manifest format_version $(repr(get(manifest, "format_version", nothing))) is not $(POWSYBL_BUNDLE_FORMAT_VERSION)"))
  return manifest
end

function _powsybl_read_table(dir::AbstractString, name::AbstractString, entry)::NamedTuple
  entry isa AbstractDict || throw(_powsybl_bundle_error(dir, "manifest entry of table $(name) is not an object"))
  file = get(entry, "file", nothing)
  file isa AbstractString || throw(_powsybl_bundle_error(dir, "table $(name) has no file name in the manifest"))
  path = joinpath(dir, file)
  isfile(path) || throw(_powsybl_bundle_error(dir, "table $(name) has no file $(file)"))
  data, header = readdlm(path, ',', String; quotes = true, header = true)
  names = Symbol[Symbol(strip(String(h))) for h in vec(header)]
  # alias map: an older or newer pypowsybl column name becomes the
  # canonical one before the schema check
  aliases = get(POWSYBL_COLUMN_ALIASES, String(name), nothing)
  if aliases !== nothing
    for (canonical, olds) in aliases, (k, n) in enumerate(names)
      n in olds && (names[k] = canonical)
    end
  end
  schema = POWSYBL_SCHEMA[String(name)]
  types = Dict{Symbol,DataType}(c.name => c.eltype for c in schema)
  for c in schema
    c.required && !(c.name in names) && throw(_powsybl_bundle_error(dir, "table $(name) misses required column $(c.name)"))
  end
  rows = get(entry, "rows", nothing)
  rows isa Number && size(data, 1) != rows && throw(_powsybl_bundle_error(dir, "table $(name) has $(size(data, 1)) rows, the manifest says $(rows)"))
  columns = Any[]
  for (k, n) in enumerate(names)
    cells = size(data, 1) == 0 ? String[] : view(data, :, k)
    push!(columns, _powsybl_convert_column(cells, get(types, n, String), dir, name, n))
  end
  return NamedTuple{Tuple(names)}(Tuple(columns))
end

"""
    read_powsybl_bundle(dir) -> PowsyblTables

Read a PowSyBl table bundle: `manifest.json` (format `powsybl_tables`,
version 1) and every table of [`POWSYBL_TABLE_NAMES`](@ref) from the file
the manifest names. Column types follow [`POWSYBL_SCHEMA`](@ref), never
the data; a required column that is missing, a table that the manifest
does not list or whose file is absent, and a row count that disagrees
with the manifest each fail with an `ArgumentError` naming the bundle,
the table and the column. Columns the schema does not know are kept as
strings. A tap changer row with an empty `tap` cell carries no step: the
changer sits at the middle of `low_tap..high_tap`, the transformer columns
at the current tap are recomputed from the step rows, and the manifest key
`tap_position_assumed` lists the changer for the import notice and the
tap estimation.
"""
function read_powsybl_bundle(dir::AbstractString)::PowsyblTables
  isdir(dir) || throw(_powsybl_bundle_error(dir, "not a directory"))
  manifest = read_powsybl_manifest(dir)
  listed = get(manifest, "tables", nothing)
  listed isa AbstractDict || throw(_powsybl_bundle_error(dir, "manifest has no tables object"))
  kwargs = Dict{Symbol,Any}()
  for name in POWSYBL_TABLE_NAMES
    entry = get(listed, name, nothing)
    entry === nothing && throw(_powsybl_bundle_error(dir, "table $(name) is not listed in the manifest"))
    kwargs[POWSYBL_TABLE_FIELDS[name]] = _powsybl_read_table(dir, name, entry)
  end
  _powsybl_bundle_tap_middle!(dir, kwargs, manifest)
  return PowsyblTables(; manifest = manifest, kwargs...)
end

# A tap changer row with an empty `tap` cell carries no step: the changer
# sits at the middle of low_tap..high_tap (IIDM has no neutral step, the
# rule of the IIDM reader). The transformer columns at the current tap
# (rho, alpha, r/x/g/b at current tap, the 3W tap positions) describe a
# position the file does not name, so they are recomputed from the step
# rows at the taps now set, as the IIDM reader builds them. The changers
# go into the manifest key tap_position_assumed, which the mapping reads
# for the mark and the import notice. An empty 3W tap position without
# such a changer stays an error.
function _powsybl_bundle_tap_middle!(dir::AbstractString, kw::Dict{Symbol,Any}, manifest::Dict{String,Any})
  assumed = Dict{String,Any}[]
  for (field, role) in ((:ratio_tap_changers, "ratio"), (:phase_tap_changers, "phase"))
    t = kw[field]
    for i in eachindex(t.tap)
      t.tap[i] == _POWSYBL_NO_TAP || continue
      t.tap[i] = Sparlectra._tap_band_middle(minmax(t.low_tap[i], t.high_tap[i])...)
      push!(assumed, Dict{String,Any}("id" => t.id[i], "side" => t.side[i], "role" => role, "tap" => t.tap[i]))
    end
  end
  # the step row of the changer (id, side) in `role` at its tap, nothing
  # without a changer
  function step_row(role::String, id::String, side::String)
    t = kw[role == "ratio" ? :ratio_tap_changers : :phase_tap_changers]
    st = kw[role == "ratio" ? :ratio_tap_changer_steps : :phase_tap_changer_steps]
    i = findfirst(k -> t.id[k] == id && t.side[k] == side, eachindex(t.id))
    i === nothing && return nothing
    j = findfirst(k -> st.id[k] == id && st.side[k] == side && st.position[k] == t.tap[i], eachindex(st.id))
    j === nothing && throw(_powsybl_bundle_error(dir, "tap changer $(id)$(isempty(side) ? "" : "/" * side)/$(role): tap $(t.tap[i]) has no step row"))
    return (st = st, j = j, tap = t.tap[i])
  end
  factor(sr, key::Symbol) = sr === nothing ? 1.0 : 1.0 + sr.st[key][sr.j] / 100.0
  rho(sr) = sr === nothing ? 1.0 : sr.st.rho[sr.j]
  tw = kw[:two_windings_transformers]
  t3 = kw[:three_windings_transformers]
  for (id, side) in unique([(d["id"], d["side"]) for d in assumed])
    rr = step_row("ratio", id, side)
    pr = step_row("phase", id, side)
    alpha = pr === nothing ? 0.0 : pr.st.alpha[pr.j]
    if isempty(side)
      i = findfirst(==(id), tw.id)
      i === nothing && throw(_powsybl_bundle_error(dir, "tap changer $(id) has no row in 2_windings_transformers"))
      tw.rho[i] = tw.rated_u2[i] / tw.rated_u1[i] * rho(rr) * rho(pr)
      tw.alpha[i] = alpha
      for key in (:r, :x, :g, :b)
        tw[Symbol(key, "_at_current_tap")][i] = tw[key][i] * factor(rr, key) * factor(pr, key)
      end
    else
      i = findfirst(==(id), t3.id)
      k = findfirst(==(side), ("ONE", "TWO", "THREE"))
      (i === nothing || k === nothing) && throw(_powsybl_bundle_error(dir, "tap changer $(id)/$(side) has no leg in 3_windings_transformers"))
      t3[Symbol("rho", k)][i] = t3.rated_u0[i] / t3[Symbol("rated_u", k)][i] * rho(rr) * rho(pr)
      t3[Symbol("alpha", k)][i] = alpha
      for key in (:r, :x, :g, :b)
        t3[Symbol(key, k, "_at_current_tap")][i] = t3[Symbol(key, k)][i] * factor(rr, key) * factor(pr, key)
      end
      t3[Symbol("ratio_tap_position", k)][i] = rr === nothing ? t3[Symbol("ratio_tap_position", k)][i] : rr.tap
      t3[Symbol("phase_tap_position", k)][i] = pr === nothing ? t3[Symbol("phase_tap_position", k)][i] : pr.tap
    end
  end
  for c in _POWSYBL_TAP_CELLS["3_windings_transformers"], i in eachindex(t3.id)
    t3[c][i] == _POWSYBL_NO_TAP && throw(_powsybl_bundle_error(dir, "table 3_windings_transformers column $(c) row $(i): \"\" is not an integer"))
  end
  isempty(assumed) || (manifest["tap_position_assumed"] = vcat(Any[get(manifest, "tap_position_assumed", Any[])...], assumed))
  return nothing
end

# --- writer -----------------------------------------------------------------

function _powsybl_csv_string(s::AbstractString)::String
  return string('"', replace(String(s), "\"" => "\"\""), '"')
end

_powsybl_csv_cell(v::AbstractString) = _powsybl_csv_string(v)
_powsybl_csv_cell(v::Bool) = v ? "true" : "false"
_powsybl_csv_cell(v::Integer) = string(v)
_powsybl_csv_cell(v::AbstractFloat) = repr(Float64(v))
_powsybl_csv_cell(v) = _powsybl_csv_string(string(v))

function _powsybl_write_table(io::IO, table::NamedTuple)
  names = collect(keys(table))
  println(io, join((_powsybl_csv_string(String(n)) for n in names), ","))
  for row in 1:powsybl_table_rows(table)
    println(io, join((_powsybl_csv_cell(table[n][row]) for n in names), ","))
  end
  return nothing
end

"""
    write_powsybl_bundle(tables::PowsyblTables, dir; manifest_extra = Dict()) -> String

Write a bundle directory: one CSV per table of
[`POWSYBL_TABLE_NAMES`](@ref) (`<table>.csv`, index columns first as the
NamedTuple carries them) and `manifest.json`. The manifest starts from
`tables.manifest`, so a read bundle round-trips with its case, source and
reference fields, then sets `format`, `format_version`, `written` (UTC
now) and the `tables` entries (file, rows, index), and finally takes every
key of `manifest_extra`. Returns the manifest path.
"""
function write_powsybl_bundle(tables::PowsyblTables, dir::AbstractString; manifest_extra::AbstractDict = Dict{String,Any}())::String
  mkpath(dir)
  manifest = Dict{String,Any}(String(k) => v for (k, v) in tables.manifest)
  manifest["format"] = POWSYBL_BUNDLE_FORMAT
  manifest["format_version"] = POWSYBL_BUNDLE_FORMAT_VERSION
  manifest["written"] = Dates.format(Dates.now(Dates.UTC), "yyyy-mm-ddTHH:MM:SSZ")
  entries = Dict{String,Any}()
  for name in POWSYBL_TABLE_NAMES
    table = powsybl_table(tables, name)
    file = name * ".csv"
    open(joinpath(dir, file), "w") do io
      _powsybl_write_table(io, table)
    end
    index = String[String(c.name) for c in POWSYBL_SCHEMA[name] if c.index]
    entries[name] = Dict{String,Any}("file" => file, "rows" => powsybl_table_rows(table), "index" => index)
  end
  manifest["tables"] = entries
  for (k, v) in manifest_extra
    manifest[String(k)] = v
  end
  path = joinpath(dir, POWSYBL_MANIFEST_FILE)
  open(path, "w") do io
    write(io, scf_json_string(manifest))
  end
  return path
end
