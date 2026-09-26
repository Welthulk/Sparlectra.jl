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

# file: src/adapters/powsybl/powsybl_tables.jl
# purpose: the in-memory form of a PowSyBl table bundle: one column-vector
#          table per pypowsybl getter, the schema that fixes the column
#          types (never inferred from data), and the empty-table
#          constructors. PowSyBl resolves the IIDM file, the node-breaker
#          topology and the tap positions; Sparlectra reads only these
#          tables. No DataFrames, no Tables.jl: a table is a NamedTuple of
#          vectors, the keys are the pypowsybl column names.

"""
    PowsyblColumn

One column of the bundle schema: its name, its element type (`String`,
`Float64`, `Int` or `Bool`), whether it is an index column of the
pypowsybl table, and whether the network builder reads it (`required`).
A column the builder never reads is read when present and ignored when
absent; a required column that is missing fails the bundle read.
"""
struct PowsyblColumn
  name::Symbol
  eltype::DataType
  index::Bool
  required::Bool
end

_pc(name::Symbol, T::DataType; index::Bool = false, required::Bool = true) = PowsyblColumn(name, T, index, required)

# Column groups that repeat across tables. The pypowsybl order is kept
# inside every table definition, so a written bundle lists the columns as
# the getter would.
const _POWSYBL_BUS_ATTACHMENT = [_pc(:voltage_level_id, String), _pc(:bus_id, String), _pc(:bus_breaker_bus_id, String), _pc(:node, Int; required = false), _pc(:connected, Bool)]
const _POWSYBL_RESULT_PQI = [_pc(:p, Float64; required = false), _pc(:q, Float64; required = false), _pc(:i, Float64; required = false)]
const _POWSYBL_FICTITIOUS = [_pc(:fictitious, Bool; required = false)]

function _powsybl_branch_side_columns(k::Int)
  return [
    _pc(Symbol("voltage_level$(k)_id"), String),
    _pc(Symbol("bus$(k)_id"), String),
    _pc(Symbol("bus_breaker_bus$(k)_id"), String),
    _pc(Symbol("node$(k)"), Int; required = false),
  ]
end

function _powsybl_3wt_leg_columns(k::Int)
  return [
    _pc(Symbol("r$(k)"), Float64), _pc(Symbol("x$(k)"), Float64), _pc(Symbol("g$(k)"), Float64), _pc(Symbol("b$(k)"), Float64),
    _pc(Symbol("rated_u$(k)"), Float64), _pc(Symbol("rated_s$(k)"), Float64),
    _pc(Symbol("ratio_tap_position$(k)"), Int), _pc(Symbol("phase_tap_position$(k)"), Int),
    _pc(Symbol("p$(k)"), Float64; required = false), _pc(Symbol("q$(k)"), Float64; required = false), _pc(Symbol("i$(k)"), Float64; required = false),
    _powsybl_branch_side_columns(k)...,
    _pc(Symbol("connected$(k)"), Bool),
    _pc(Symbol("selected_limits_group_$(k)"), String; required = false),
    _pc(Symbol("rho$(k)"), Float64), _pc(Symbol("alpha$(k)"), Float64),
    _pc(Symbol("r$(k)_at_current_tap"), Float64), _pc(Symbol("x$(k)_at_current_tap"), Float64),
    _pc(Symbol("g$(k)_at_current_tap"), Float64), _pc(Symbol("b$(k)_at_current_tap"), Float64),
  ]
end

"""
    POWSYBL_SCHEMA

The bundle schema, keyed by the manifest table name (the pypowsybl getter
name without `get_`): the columns the reader knows, in pypowsybl order,
with type, index flag and required flag. Every table listed here must be
present in a bundle, with at least its header line. Columns of a bundle
that the schema does not list are read as `String` and carried along
unchanged, so a newer pypowsybl adds nothing that breaks a read.
"""
const POWSYBL_SCHEMA = Dict{String,Vector{PowsyblColumn}}(
  "substations" => [_pc(:id, String; index = true), _pc(:name, String; required = false), _pc(:TSO, String; required = false), _pc(:geo_tags, String; required = false), _pc(:country, String; required = false), _POWSYBL_FICTITIOUS...],
  "voltage_levels" => [_pc(:id, String; index = true), _pc(:name, String; required = false), _pc(:substation_id, String), _pc(:nominal_v, Float64), _pc(:high_voltage_limit, Float64; required = false), _pc(:low_voltage_limit, Float64; required = false), _POWSYBL_FICTITIOUS..., _pc(:topology_kind, String; required = false)],
  "buses" => [_pc(:id, String; index = true), _pc(:name, String; required = false), _pc(:v_mag, Float64), _pc(:v_angle, Float64), _pc(:connected_component, Int), _pc(:synchronous_component, Int), _pc(:voltage_level_id, String), _POWSYBL_FICTITIOUS..., _pc(:fictitious_p0, Float64; required = false), _pc(:fictitious_q0, Float64; required = false)],
  "bus_breaker_view_buses" => [_pc(:id, String; index = true), _pc(:name, String; required = false), _pc(:v_mag, Float64), _pc(:v_angle, Float64), _pc(:connected_component, Int), _pc(:synchronous_component, Int), _pc(:voltage_level_id, String), _pc(:bus_id, String), _POWSYBL_FICTITIOUS..., _pc(:fictitious_p0, Float64; required = false), _pc(:fictitious_q0, Float64; required = false)],
  "switches" => [_pc(:id, String; index = true), _pc(:name, String; required = false), _pc(:kind, String), _pc(:open, Bool), _pc(:retained, Bool), _pc(:voltage_level_id, String), _pc(:bus_breaker_bus1_id, String), _pc(:bus_breaker_bus2_id, String), _pc(:node1, Int; required = false), _pc(:node2, Int; required = false), _POWSYBL_FICTITIOUS...],
  "busbar_sections" => [_pc(:id, String; index = true), _pc(:name, String; required = false), _pc(:v, Float64; required = false), _pc(:angle, Float64; required = false), _pc(:voltage_level_id, String; required = false), _pc(:bus_id, String; required = false), _pc(:bus_breaker_bus_id, String; required = false), _pc(:node, Int; required = false), _pc(:connected, Bool; required = false), _POWSYBL_FICTITIOUS...],
  "lines" => [
    _pc(:id, String; index = true), _pc(:name, String; required = false),
    _pc(:r, Float64), _pc(:x, Float64), _pc(:g1, Float64), _pc(:b1, Float64), _pc(:g2, Float64), _pc(:b2, Float64),
    _pc(:p1, Float64; required = false), _pc(:q1, Float64; required = false), _pc(:i1, Float64; required = false),
    _pc(:p2, Float64; required = false), _pc(:q2, Float64; required = false), _pc(:i2, Float64; required = false),
    _powsybl_branch_side_columns(1)..., _powsybl_branch_side_columns(2)...,
    _pc(:connected1, Bool), _pc(:connected2, Bool), _POWSYBL_FICTITIOUS...,
    _pc(:selected_limits_group_1, String; required = false), _pc(:selected_limits_group_2, String; required = false),
  ],
  "2_windings_transformers" => [
    _pc(:id, String; index = true), _pc(:name, String; required = false),
    _pc(:r, Float64), _pc(:x, Float64), _pc(:g, Float64), _pc(:b, Float64),
    _pc(:rated_u1, Float64), _pc(:rated_u2, Float64), _pc(:rated_s, Float64),
    _pc(:p1, Float64; required = false), _pc(:q1, Float64; required = false), _pc(:i1, Float64; required = false),
    _pc(:p2, Float64; required = false), _pc(:q2, Float64; required = false), _pc(:i2, Float64; required = false),
    _powsybl_branch_side_columns(1)..., _powsybl_branch_side_columns(2)...,
    _pc(:connected1, Bool), _pc(:connected2, Bool), _POWSYBL_FICTITIOUS...,
    _pc(:selected_limits_group_1, String; required = false), _pc(:selected_limits_group_2, String; required = false),
    _pc(:rho, Float64), _pc(:alpha, Float64),
    _pc(:r_at_current_tap, Float64), _pc(:x_at_current_tap, Float64), _pc(:g_at_current_tap, Float64), _pc(:b_at_current_tap, Float64),
  ],
  "3_windings_transformers" => [
    _pc(:id, String; index = true), _pc(:name, String; required = false), _pc(:rated_u0, Float64),
    _powsybl_3wt_leg_columns(1)..., _powsybl_3wt_leg_columns(2)..., _powsybl_3wt_leg_columns(3)...,
    _POWSYBL_FICTITIOUS...,
  ],
  "ratio_tap_changers" => [_pc(:id, String; index = true), _pc(:side, String; required = false), _pc(:tap, Int; required = false), _pc(:solved_tap_position, Float64; required = false), _pc(:low_tap, Int; required = false), _pc(:high_tap, Int; required = false), _pc(:step_count, Int; required = false), _pc(:oltc, Bool; required = false), _pc(:regulating, Bool; required = false), _pc(:target_v, Float64; required = false), _pc(:target_deadband, Float64; required = false), _pc(:regulating_bus_id, String; required = false), _pc(:regulated_side, String; required = false)],
  "ratio_tap_changer_steps" => [_pc(:id, String; index = true), _pc(:position, Int; index = true), _pc(:side, String; required = false), _pc(:rho, Float64; required = false), _pc(:r, Float64; required = false), _pc(:x, Float64; required = false), _pc(:g, Float64; required = false), _pc(:b, Float64; required = false)],
  "phase_tap_changers" => [_pc(:id, String; index = true), _pc(:side, String; required = false), _pc(:tap, Int; required = false), _pc(:solved_tap_position, Float64; required = false), _pc(:low_tap, Int; required = false), _pc(:high_tap, Int; required = false), _pc(:step_count, Int; required = false), _pc(:oltc, Bool; required = false), _pc(:regulating, Bool; required = false), _pc(:regulation_mode, String; required = false), _pc(:regulation_value, Float64; required = false), _pc(:target_deadband, Float64; required = false), _pc(:regulating_bus_id, String; required = false), _pc(:regulated_side, String; required = false)],
  "phase_tap_changer_steps" => [_pc(:id, String; index = true), _pc(:position, Int; index = true), _pc(:side, String; required = false), _pc(:rho, Float64; required = false), _pc(:alpha, Float64; required = false), _pc(:r, Float64; required = false), _pc(:x, Float64; required = false), _pc(:g, Float64; required = false), _pc(:b, Float64; required = false)],
  "generators" => [
    _pc(:id, String; index = true), _pc(:name, String; required = false), _pc(:energy_source, String; required = false),
    _pc(:target_p, Float64), _pc(:min_p, Float64), _pc(:max_p, Float64), _pc(:min_q, Float64), _pc(:max_q, Float64),
    _pc(:min_q_at_target_p, Float64), _pc(:max_q_at_target_p, Float64), _pc(:min_q_at_p, Float64; required = false), _pc(:max_q_at_p, Float64; required = false),
    _pc(:rated_s, Float64; required = false), _pc(:reactive_limits_kind, String), _pc(:target_v, Float64), _pc(:equivalent_local_target_v, Float64; required = false), _pc(:target_q, Float64),
    _pc(:voltage_regulator_on, Bool), _pc(:regulated_element_id, String; required = false), _pc(:regulated_bus_id, String), _pc(:regulated_bus_breaker_bus_id, String),
    _POWSYBL_RESULT_PQI..., _POWSYBL_BUS_ATTACHMENT..., _POWSYBL_FICTITIOUS..., _pc(:condenser, Bool; required = false),
  ],
  "reactive_capability_curve_points" => [_pc(:id, String; index = true), _pc(:num, Int; index = true), _pc(:p, Float64; required = false), _pc(:min_q, Float64; required = false), _pc(:max_q, Float64; required = false)],
  "loads" => [_pc(:id, String; index = true), _pc(:name, String; required = false), _pc(:type, String), _pc(:p0, Float64), _pc(:q0, Float64), _POWSYBL_RESULT_PQI..., _POWSYBL_BUS_ATTACHMENT..., _POWSYBL_FICTITIOUS...],
  "shunt_compensators" => [
    _pc(:id, String; index = true), _pc(:name, String; required = false), _pc(:g, Float64), _pc(:b, Float64), _pc(:model_type, String),
    _pc(:max_section_count, Int), _pc(:section_count, Int), _pc(:solved_section_count, Float64; required = false),
    _pc(:voltage_regulation_on, Bool; required = false), _pc(:target_v, Float64; required = false), _pc(:target_deadband, Float64; required = false), _pc(:regulating_bus_id, String; required = false),
    _POWSYBL_RESULT_PQI..., _POWSYBL_BUS_ATTACHMENT..., _POWSYBL_FICTITIOUS...,
  ],
  "linear_shunt_compensator_sections" => [_pc(:id, String; index = true), _pc(:g_per_section, Float64; required = false), _pc(:b_per_section, Float64; required = false), _pc(:max_section_count, Int; required = false)],
  "dangling_lines" => [
    _pc(:id, String; index = true), _pc(:name, String; required = false),
    _pc(:r, Float64), _pc(:x, Float64), _pc(:g, Float64), _pc(:b, Float64), _pc(:p0, Float64), _pc(:q0, Float64),
    _POWSYBL_RESULT_PQI...,
    _pc(:boundary_p, Float64; required = false), _pc(:boundary_q, Float64; required = false), _pc(:boundary_i, Float64; required = false),
    _pc(:boundary_v_mag, Float64; required = false), _pc(:boundary_v_angle, Float64; required = false),
    _POWSYBL_BUS_ATTACHMENT...,
    _pc(:pairing_key, String), _pc(:ucte_xnode_code, String; required = false), _pc(:paired, Bool), _POWSYBL_FICTITIOUS..., _pc(:tie_line_id, String), _pc(:selected_limits_group, String; required = false),
  ],
  "tie_lines" => [_pc(:id, String; index = true), _pc(:name, String; required = false), _pc(:boundary_line1_id, String; required = false), _pc(:dangling_line1_id, String), _pc(:boundary_line2_id, String; required = false), _pc(:dangling_line2_id, String), _pc(:pairing_key, String), _pc(:ucte_xnode_code, String; required = false), _pc(:connected1, Bool), _pc(:connected2, Bool), _POWSYBL_FICTITIOUS...],
  "hvdc_lines" => [_pc(:id, String; index = true), _pc(:name, String; required = false), _pc(:converters_mode, String), _pc(:target_p, Float64), _pc(:max_p, Float64), _pc(:nominal_v, Float64), _pc(:r, Float64), _pc(:converter_station1_id, String), _pc(:converter_station2_id, String), _pc(:connected1, Bool), _pc(:connected2, Bool), _POWSYBL_FICTITIOUS...],
  "vsc_converter_stations" => [
    _pc(:id, String; index = true), _pc(:name, String; required = false), _pc(:loss_factor, Float64),
    _pc(:min_q, Float64), _pc(:max_q, Float64), _pc(:min_q_at_target_p, Float64), _pc(:max_q_at_target_p, Float64), _pc(:min_q_at_p, Float64; required = false), _pc(:max_q_at_p, Float64; required = false),
    _pc(:reactive_limits_kind, String), _pc(:target_v, Float64), _pc(:target_q, Float64), _pc(:voltage_regulator_on, Bool),
    _pc(:regulated_element_id, String; required = false), _pc(:regulated_bus_id, String), _pc(:regulated_bus_breaker_bus_id, String; required = false),
    _POWSYBL_RESULT_PQI..., _POWSYBL_BUS_ATTACHMENT..., _POWSYBL_FICTITIOUS..., _pc(:hvdc_line_id, String),
  ],
  "lcc_converter_stations" => [_pc(:id, String; index = true), _pc(:name, String; required = false), _pc(:power_factor, Float64), _pc(:loss_factor, Float64), _POWSYBL_RESULT_PQI..., _POWSYBL_BUS_ATTACHMENT..., _POWSYBL_FICTITIOUS..., _pc(:hvdc_line_id, String)],
  "static_var_compensators" => [
    _pc(:id, String; index = true), _pc(:name, String; required = false), _pc(:b_min, Float64), _pc(:b_max, Float64), _pc(:target_v, Float64), _pc(:target_q, Float64),
    _pc(:regulation_mode, String), _pc(:regulating, Bool), _pc(:regulated_element_id, String; required = false), _pc(:regulated_bus_id, String; required = false), _pc(:regulated_bus_breaker_bus_id, String; required = false),
    _POWSYBL_RESULT_PQI..., _POWSYBL_BUS_ATTACHMENT..., _POWSYBL_FICTITIOUS...,
  ],
  "operational_limits" => [_pc(:element_id, String; index = true), _pc(:side, String; index = true), _pc(:type, String; index = true), _pc(:acceptable_duration, Int; index = true), _pc(:group_name, String; index = true), _pc(:element_type, String), _pc(:name, String), _pc(:value, Float64), _POWSYBL_FICTITIOUS..., _pc(:selected, Bool; required = false)],
)

"""
    POWSYBL_TABLE_NAMES

The manifest table names in the order of the pypowsybl getters (and of
the bundle files). `tie_lines` closes the list: the getter exists in
pypowsybl 1.16.1 and pairs the dangling lines of a tie line.
"""
const POWSYBL_TABLE_NAMES = String[
  "substations", "voltage_levels", "buses", "bus_breaker_view_buses", "switches", "busbar_sections",
  "lines", "2_windings_transformers", "3_windings_transformers",
  "ratio_tap_changers", "ratio_tap_changer_steps", "phase_tap_changers", "phase_tap_changer_steps",
  "generators", "reactive_capability_curve_points", "loads", "shunt_compensators", "linear_shunt_compensator_sections",
  "dangling_lines", "tie_lines", "hvdc_lines", "vsc_converter_stations", "lcc_converter_stations", "static_var_compensators",
  "operational_limits",
]

# Two manifest names start with a digit and cannot be Julia field names;
# every other field name equals its manifest name.
const POWSYBL_TABLE_FIELDS = Dict{String,Symbol}(name => (name == "2_windings_transformers" ? :two_windings_transformers : name == "3_windings_transformers" ? :three_windings_transformers : Symbol(name)) for name in POWSYBL_TABLE_NAMES)

"""
    POWSYBL_COLUMN_ALIASES

Column renames between pypowsybl versions, per table: canonical name to
the names older or newer releases used. Empty: pypowsybl 1.16.1 (the
version the schema was verified on) renamed nothing. The reader maps an
alias to its canonical name before the schema check.
"""
const POWSYBL_COLUMN_ALIASES = Dict{String,Dict{Symbol,Vector{Symbol}}}()

_powsybl_vector_type(T::DataType) = T === String ? Vector{String} : T === Float64 ? Vector{Float64} : T === Int ? Vector{Int} : Vector{Bool}

"""
    powsybl_empty_table(name) -> NamedTuple

An empty table with every schema column of `name` (required and
optional), in schema order, with typed empty vectors.
"""
function powsybl_empty_table(name::AbstractString)::NamedTuple
  columns = POWSYBL_SCHEMA[String(name)]
  return NamedTuple{Tuple(c.name for c in columns)}(Tuple(_powsybl_vector_type(c.eltype)() for c in columns))
end

"""
    powsybl_table_rows(table::NamedTuple) -> Int

The row count of a table: the length of its first column, 0 for a table
without columns.
"""
powsybl_table_rows(table::NamedTuple)::Int = isempty(table) ? 0 : length(first(table))

"""
    PowsyblTables

A PowSyBl table bundle in memory: the manifest as a `Dict` and one
`NamedTuple` of column vectors per pypowsybl table (field names as in
`POWSYBL_TABLE_FIELDS`). The keyword constructor defaults every
table to its empty schema form, so a table that a source does not carry
is present with zero rows, never absent.
"""
struct PowsyblTables
  manifest::Dict{String,Any}
  substations::NamedTuple
  voltage_levels::NamedTuple
  buses::NamedTuple
  bus_breaker_view_buses::NamedTuple
  switches::NamedTuple
  busbar_sections::NamedTuple
  lines::NamedTuple
  two_windings_transformers::NamedTuple
  three_windings_transformers::NamedTuple
  ratio_tap_changers::NamedTuple
  ratio_tap_changer_steps::NamedTuple
  phase_tap_changers::NamedTuple
  phase_tap_changer_steps::NamedTuple
  generators::NamedTuple
  reactive_capability_curve_points::NamedTuple
  loads::NamedTuple
  shunt_compensators::NamedTuple
  linear_shunt_compensator_sections::NamedTuple
  dangling_lines::NamedTuple
  tie_lines::NamedTuple
  hvdc_lines::NamedTuple
  vsc_converter_stations::NamedTuple
  lcc_converter_stations::NamedTuple
  static_var_compensators::NamedTuple
  operational_limits::NamedTuple
end

function PowsyblTables(; manifest::AbstractDict = Dict{String,Any}(), kwargs...)
  fields = Any[Dict{String,Any}(String(k) => v for (k, v) in manifest)]
  for name in POWSYBL_TABLE_NAMES
    field = POWSYBL_TABLE_FIELDS[name]
    push!(fields, get(kwargs, field, nothing) === nothing ? powsybl_empty_table(name) : kwargs[field])
  end
  return PowsyblTables(fields...)
end

"""
    powsybl_table(tables::PowsyblTables, name) -> NamedTuple

The table of the manifest name `name` (`"2_windings_transformers"` and the
like), so callers can loop over [`POWSYBL_TABLE_NAMES`](@ref).
"""
powsybl_table(tables::PowsyblTables, name::AbstractString)::NamedTuple = getfield(tables, POWSYBL_TABLE_FIELDS[String(name)])

function Base.show(io::IO, ::MIME"text/plain", tables::PowsyblTables)
  case = get(tables.manifest, "case", "")
  println(io, "PowsyblTables", isempty(case) ? "" : " ($(case))", ":")
  width = maximum(length, POWSYBL_TABLE_NAMES)
  for name in POWSYBL_TABLE_NAMES
    println(io, "  ", rpad(name, width), " ", powsybl_table_rows(powsybl_table(tables, name)), " rows")
  end
  return nothing
end
