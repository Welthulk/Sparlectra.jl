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

# file: ext/SparlectraPythonCallExt.jl
# purpose: the optional live path of the PowSyBl adapter: with PythonCall
#          and pypowsybl in the session, an IIDM file is read straight into
#          PowsyblTables (every table with all attributes, the same column
#          names and types as a bundle written by tools/powsybl_dump.py),
#          OpenLoadFlow provides the reference solution, and a bundle can be
#          written from Julia. Loaded automatically when PythonCall is in
#          the environment; without it the core's generic functions have no
#          method and import_net says how to produce a bundle instead.

module SparlectraPythonCallExt

using Sparlectra
using PythonCall

const _EXT_TABLE_NAMES = Sparlectra.POWSYBL_TABLE_NAMES

# The Python executable hint of the powsybl_import scope: honoured only
# before PythonCall initialised its interpreter, which happened when this
# extension loaded. The function exists so a caller can state the hint and
# learn whether it took effect.
function _apply_python_exe_hint(python_exe::AbstractString)::Bool
  isempty(python_exe) && return true
  current = try
    pyconvert(String, pyimport("sys").executable)
  catch err
    # the interpreter is up or not; either way the hint cannot change it
    @warn "powsybl_import.python_exe: the running Python could not be queried; the hint $(python_exe) has no effect" exception = err
    return false
  end
  current == python_exe && return true
  @warn "powsybl_import.python_exe: PythonCall is already initialised with $(current); the hint $(python_exe) has no effect in this session (set JULIA_PYTHONCALL_EXE before Julia starts)"
  return false
end

_pyimport_powsybl() = pyimport("pypowsybl")

# A string cell the way tools/powsybl_dump.py writes it: None and a pandas
# NA become the empty string, a float in a string column keeps its repr
# (NaN stays "NaN"), everything else is str(value).
function _string_cell(v)::String
  pyis(v, pybuiltins.None) && return ""
  if pyisinstance(v, pybuiltins.float)
    f = pyconvert(Float64, v)
    return isnan(f) ? "NaN" : (isinf(f) ? (f > 0 ? "Inf" : "-Inf") : repr(f))
  end
  pyconvert(Bool, pyimport("pandas").isna(v)) && return ""
  return pyconvert(String, pystr(v))
end

# One pandas column (a Python list, so every cell is still a Py object) to
# the typed Julia vector of the schema: NaN stays NaN in Float64 columns, a
# missing integer becomes -1 (the dump script's rule), strings follow
# _string_cell, booleans become Bool; unknown columns stay String.
function _column_vector(values, T::DataType)
  if T === Float64
    return Float64[pyis(v, pybuiltins.None) ? NaN : pyconvert(Float64, v) for v in values]
  elseif T === Int
    out = Int[]
    for v in values
      f = pyis(v, pybuiltins.None) ? NaN : pyconvert(Float64, v)
      push!(out, isfinite(f) ? Int(round(f)) : -1)
    end
    return out
  elseif T === Bool
    return Bool[pyconvert(Bool, v) for v in values]
  end
  return String[_string_cell(v) for v in values]
end

# A pypowsybl table (after reset_index) into the NamedTuple of a bundle
# table, index columns first, then the value columns in pypowsybl order.
function _table_from_dataframe(name::AbstractString, df)
  df = df.reset_index()
  names = [pyconvert(String, c) for c in df.columns]
  types = Dict{Symbol,DataType}(c.name => c.eltype for c in Sparlectra.POWSYBL_SCHEMA[String(name)])
  columns = Any[]
  for col in names
    T = get(types, Symbol(col), String)
    values = df[col].tolist()
    push!(columns, _column_vector(values, T))
  end
  return NamedTuple{Tuple(Symbol.(names))}(Tuple(columns))
end

function _get_table(network, name::AbstractString)
  getter = pygetattr(network, "get_" * name)
  df = try
    getter(all_attributes = true)
  catch err
    # a getter without the keyword (older pypowsybl) still answers without it
    err isa PyException && occursin("all_attributes", sprint(showerror, err)) || rethrow()
    getter()
  end
  return _table_from_dataframe(name, df)
end

function _iidm_version(pp, path::AbstractString)::String
  endswith(lowercase(path), ".bz2") && return ""
  head = try
    open(io -> String(read(io, 4096)), path, "r")
  catch
    return ""
  end
  m = match(r"schema/iidm/([0-9_]+)", head)
  return m === nothing ? "" : String(m.captures[1])
end

function _load_network(pp, path::AbstractString)
  return pp.network.load(path)
end

function _tables_of(pp, network, path::AbstractString, case::AbstractString)::Sparlectra.PowsyblTables
  kwargs = Dict{Symbol,Any}()
  for name in _EXT_TABLE_NAMES
    kwargs[Sparlectra.POWSYBL_TABLE_FIELDS[name]] = _get_table(network, name)
  end
  manifest = Dict{String,Any}(
    "case" => case,
    "source" => String(path),
    "source_file" => basename(String(path)),
    "pypowsybl_version" => pyconvert(String, pp.__version__),
    "iidm_version" => _iidm_version(pp, path),
    "all_attributes" => true,
    "format" => "powsybl_tables",
    "format_version" => 1,
  )
  return Sparlectra.PowsyblTables(; manifest = manifest, kwargs...)
end

# The reference load flow with the tolerances tools/powsybl_dump.py uses.
function _run_reference(pp, network)
  parameters = pp.loadflow.Parameters(provider_parameters = pydict(Dict("slackBusPMaxMismatch" => "0.0001", "newtonRaphsonConvEpsPerEq" => "1e-9")))
  results = pp.loadflow.run_ac(network, parameters)
  components = Dict{String,Any}[]
  for r in results
    push!(components, Dict{String,Any}(
      "num" => pyconvert(Int, r.connected_component_num),
      "synchronous_component_num" => pyconvert(Int, r.synchronous_component_num),
      "status" => last(split(pyconvert(String, pystr(r.status)), '.')),
      "iterations" => pyconvert(Int, r.iteration_count),
      "reference_bus_id" => pyconvert(String, pystr(r.reference_bus_id)),
    ))
  end
  bb = network.get_bus_breaker_view_buses().reset_index()
  vls = network.get_voltage_levels()
  rows = NamedTuple{(:id, :bus_breaker_id, :v_mag_kv, :v_angle_deg, :v_pu, :synchronous_component, :nominal_v),Tuple{String,String,Float64,Float64,Float64,Int,Float64}}[]
  n = pyconvert(Int, pylen(bb))
  for i in 0:(n - 1)
    row = bb.iloc[i]
    vl = pyconvert(String, row["voltage_level_id"])
    nominal = pyconvert(Float64, vls.loc[vl, "nominal_v"])
    v_mag = pyconvert(Float64, row["v_mag"])
    push!(rows, (id = pyconvert(String, row["bus_id"]), bus_breaker_id = pyconvert(String, row["id"]), v_mag_kv = v_mag, v_angle_deg = pyconvert(Float64, row["v_angle"]), v_pu = nominal > 0 ? v_mag / nominal : NaN, synchronous_component = pyconvert(Int, row["synchronous_component"]), nominal_v = nominal))
  end
  return (components = components, reference_buses = rows, loadflow = "OpenLoadFlow via pypowsybl.loadflow.run_ac, default parameters except slackBusPMaxMismatch=0.0001 and newtonRaphsonConvEpsPerEq=1e-9")
end

_case_name(path::AbstractString) = begin
  stem = basename(String(path))
  for ext in (".xiidm.bz2", ".xiidm", ".xml")
    endswith(lowercase(stem), ext) && (stem = stem[1:(end - length(ext))]; break)
  end
  stem
end

"""
    read_powsybl_network(path) -> PowsyblTables

Read an IIDM file through pypowsybl: the network is loaded, OpenLoadFlow
runs (so the tables carry the solved state like a bundle written by the
dump script), and every table is converted with all attributes.
"""
function Sparlectra.read_powsybl_network(path::String)::Sparlectra.PowsyblTables
  pp = _pyimport_powsybl()
  network = _load_network(pp, path)
  reference = _run_reference(pp, network)
  tables = _tables_of(pp, network, path, _case_name(path))
  tables.manifest["reference"] = Dict{String,Any}("loadflow" => reference.loadflow, "components" => reference.components)
  return tables
end

"""
    run_powsybl_reference(path) -> NamedTuple

The OpenLoadFlow reference of an IIDM file: `components` (status,
iterations, reference bus per component) and `reference_buses` (one row
per bus-breaker bus).
"""
function Sparlectra.run_powsybl_reference(path::String)
  pp = _pyimport_powsybl()
  network = _load_network(pp, path)
  return _run_reference(pp, network)
end

"""
    dump_powsybl_bundle(path, outdir) -> String

Write the table bundle of an IIDM file from Julia: the tables, the
manifest and `reference_buses.csv`, the layout of `tools/powsybl_dump.py`.
Returns the manifest path.
"""
function Sparlectra.dump_powsybl_bundle(path::String, outdir::String)::String
  pp = _pyimport_powsybl()
  network = _load_network(pp, path)
  reference = _run_reference(pp, network)
  tables = _tables_of(pp, network, path, _case_name(path))
  mkpath(outdir)
  open(joinpath(outdir, "reference_buses.csv"), "w") do io
    println(io, "\"id\",\"bus_breaker_id\",\"v_mag_kv\",\"v_angle_deg\",\"v_pu\",\"synchronous_component\",\"nominal_v\"")
    for r in reference.reference_buses
      println(io, Sparlectra._powsybl_csv_string(r.id), ",", Sparlectra._powsybl_csv_string(r.bus_breaker_id), ",", repr(r.v_mag_kv), ",", repr(r.v_angle_deg), ",", repr(r.v_pu), ",", r.synchronous_component, ",", repr(r.nominal_v))
    end
  end
  extra = Dict{String,Any}("reference" => Dict{String,Any}("file" => "reference_buses.csv", "loadflow" => reference.loadflow, "components" => reference.components))
  return Sparlectra.write_powsybl_bundle(tables, outdir; manifest_extra = extra)
end

end # module
