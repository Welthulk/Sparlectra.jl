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
#          NaN), so a hand-edited fixture still loads; the writer never
#          produces either.

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

function _powsybl_convert_column(cells::AbstractVector{<:AbstractString}, T::DataType, dir, table, column)
  T === String && return String[String(c) for c in cells]
  T === Float64 && return Float64[_powsybl_parse_float(c, dir, table, column, i) for (i, c) in enumerate(cells)]
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
strings.
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
  return PowsyblTables(; manifest = manifest, kwargs...)
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
