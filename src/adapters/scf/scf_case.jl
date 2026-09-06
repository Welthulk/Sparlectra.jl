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

# file: src/adapters/scf/scf_case.jl
# purpose: SCFCase, the typed in-memory form of one SCF document (adapter
#          task stage 2, design decision D1). File reading and writing are
#          serialization of this struct; build_net constructs the Net from
#          it. The `data` section is typed component vectors; the
#          namespaced `sparlectra` block keeps its document sub-objects as
#          dictionaries, because they are free-form by design (extra,
#          meta, controllers) and the reader addresses them by key.

# Optional-field conventions of the row structs:
# - `nothing` on a plain optional means "key absent in the document"; the
#   writer never emits JSON null for these, so two states suffice.
# - Sensor readings are TRI-state: `missing` = key absent, `nothing` = JSON
#   null (the writer emits null for the unused half of a paired sensor),
#   a value = the reading. Serialization reproduces exactly that, which is
#   what keeps the byte-identical round trip through the struct.

const _ScfOptNum = Union{Nothing,Float64}
const _ScfOptInt = Union{Nothing,Int}
const _ScfTriNum = Union{Missing,Nothing,Float64}
const _ScfTriStr = Union{Missing,Nothing,String}

Base.@kwdef mutable struct SCFNodeRow
  id::Int
  u_rated::Float64
end

Base.@kwdef mutable struct SCFLineRow
  id::Int
  from_node::Int
  to_node::Int
  from_status::_ScfOptInt = nothing
  to_status::_ScfOptInt = nothing
  r1::Float64
  x1::Float64
  c1::_ScfOptNum = nothing
  tan1::_ScfOptNum = nothing
  # direct shunt conductance in Siemens, OUTSIDE the PGM vocabulary: tan1
  # (= g/b) cannot spell a conductance without capacitance, so a line with
  # g != 0 and b == 0 (a CGMES gch with no bch, say) carries it here; when
  # present it wins over tan1, and the strict-PGM writer must strip it
  g1::_ScfOptNum = nothing
  i_n::_ScfOptNum = nothing
end

Base.@kwdef mutable struct SCFGenericBranchRow
  id::Int
  from_node::Int
  to_node::Int
  from_status::_ScfOptInt = nothing
  to_status::_ScfOptInt = nothing
  r1::Float64
  x1::Float64
  g1::_ScfOptNum = nothing
  b1::_ScfOptNum = nothing
  k::_ScfOptNum = nothing
  theta::_ScfOptNum = nothing
  sn::_ScfOptNum = nothing
end

Base.@kwdef mutable struct SCFLinkRow
  id::Int
  from_node::Int
  to_node::Int
  from_status::_ScfOptInt = nothing
  to_status::_ScfOptInt = nothing
end

Base.@kwdef mutable struct SCFSourceRow
  id::Int
  node::Int
  status::_ScfOptInt = nothing
  u_ref::_ScfOptNum = nothing
  u_ref_angle::_ScfOptNum = nothing
  sk::_ScfOptNum = nothing
  rx_ratio::_ScfOptNum = nothing
  z01_ratio::_ScfOptNum = nothing
end

Base.@kwdef mutable struct SCFApplianceRow  # sym_load and sym_gen
  id::Int
  node::Int
  status::_ScfOptInt = nothing
  type::_ScfOptInt = nothing
  p_specified::_ScfOptNum = nothing
  q_specified::_ScfOptNum = nothing
  # only a foreign PGM sym_gen carries a setpoint here; ours live in the
  # voltage_regulator row or in extra.vm_pu
  u_ref::_ScfOptNum = nothing
end

Base.@kwdef mutable struct SCFShuntRow
  id::Int
  node::Int
  status::_ScfOptInt = nothing
  g1::_ScfOptNum = nothing
  b1::_ScfOptNum = nothing
end

Base.@kwdef mutable struct SCFVoltageRegulatorRow
  id::Int
  regulated_object::Int
  status::_ScfOptInt = nothing
  u_ref::_ScfOptNum = nothing
  q_min::_ScfOptNum = nothing
  q_max::_ScfOptNum = nothing
end

Base.@kwdef mutable struct SCFVoltageSensorRow
  id::Int
  measured_object::Int
  u_measured::_ScfTriNum = missing
  u_sigma::_ScfTriNum = missing
  u_angle_measured::_ScfTriNum = missing
end

Base.@kwdef mutable struct SCFPowerSensorRow
  id::Int
  measured_object::Int
  measured_terminal_type::_ScfTriStr = missing
  p_measured::_ScfTriNum = missing
  q_measured::_ScfTriNum = missing
  power_sigma::_ScfTriNum = missing
  p_sigma::_ScfTriNum = missing
  q_sigma::_ScfTriNum = missing
end

Base.@kwdef mutable struct SCFCurrentSensorRow
  id::Int
  measured_object::Int
  measured_terminal_type::_ScfTriStr = missing
  angle_measurement_type::_ScfTriStr = missing
  i_measured::_ScfTriNum = missing
  i_sigma::_ScfTriNum = missing
  i_angle_measured::_ScfTriNum = missing
  i_angle_sigma::_ScfTriNum = missing
end

Base.@kwdef mutable struct SCFFaultRow
  id::Int
  fault_object::Int
  status::_ScfOptInt = nothing
  fault_type::Union{Nothing,String} = nothing
  r_f::_ScfOptNum = nothing
  x_f::_ScfOptNum = nothing
end

"""
    SCFData

The PGM-compatible `data` section of one SCF document as typed component
vectors (design decision D1: no stringly-typed row access on the build
path). A component key is written to the file exactly when its vector is
non-empty, which is the writer's behavior for every revision so far.
"""
Base.@kwdef mutable struct SCFData
  node::Vector{SCFNodeRow} = SCFNodeRow[]
  line::Vector{SCFLineRow} = SCFLineRow[]
  generic_branch::Vector{SCFGenericBranchRow} = SCFGenericBranchRow[]
  link::Vector{SCFLinkRow} = SCFLinkRow[]
  source::Vector{SCFSourceRow} = SCFSourceRow[]
  sym_load::Vector{SCFApplianceRow} = SCFApplianceRow[]
  sym_gen::Vector{SCFApplianceRow} = SCFApplianceRow[]
  shunt::Vector{SCFShuntRow} = SCFShuntRow[]
  voltage_regulator::Vector{SCFVoltageRegulatorRow} = SCFVoltageRegulatorRow[]
  sym_voltage_sensor::Vector{SCFVoltageSensorRow} = SCFVoltageSensorRow[]
  sym_power_sensor::Vector{SCFPowerSensorRow} = SCFPowerSensorRow[]
  sym_current_sensor::Vector{SCFCurrentSensorRow} = SCFCurrentSensorRow[]
  fault::Vector{SCFFaultRow} = SCFFaultRow[]
end

## --- namespaced block: typed sub-blocks -------------------------------------
## roles, components and start_state are typed (stage-3 preparation): the
## adapters CONSTRUCT these, and a typo in a hand-built dict key must fail
## on write, not on the next read. The genuinely free-form sub-objects
## (per-id `extra` records, controller declarations, study definitions,
## measurement rows) stay dictionaries by design.

Base.@kwdef mutable struct SCFParticipationEntry
  object::Int
  factor::Float64
end

Base.@kwdef mutable struct SCFSlackRole
  mode::String = "single"
  nodes::Vector{Int} = Int[]
  participation::Vector{SCFParticipationEntry} = SCFParticipationEntry[]
end

Base.@kwdef mutable struct SCFRoles
  slack::Union{Nothing,SCFSlackRole} = nothing
  aux_nodes::Vector{Int} = Int[]
  isolated_nodes::Vector{Int} = Int[]
end

# positions live on an integer grid and the writer emits them as Int;
# a hand-written file may state them as floats, and the number TYPE is
# part of the byte-identical round trip, so both are kept as read
const _ScfOptGridPos = Union{Nothing,Int,Float64}

Base.@kwdef mutable struct SCFTapControllerRow
  index::Int
  alpha_deg::_ScfOptNum = nothing
  step::_ScfOptNum = nothing
  step_deg::_ScfOptNum = nothing
  pos::_ScfOptGridPos = nothing
  pos_min::_ScfOptGridPos = nothing
  pos_max::_ScfOptGridPos = nothing
end

Base.@kwdef mutable struct SCFTapChangerRow
  id::Int
  branch::Int
  side::Union{Nothing,String} = nothing
  ratio_base::_ScfOptNum = nothing
  angle_base_deg::_ScfOptNum = nothing
  controllers::Vector{SCFTapControllerRow} = SCFTapControllerRow[]
  # {"mode" => "fixed"} today; kept open for the control modes the format
  # may grow, exactly as the document carries it
  control::Dict{String,Any} = Dict{String,Any}()
  tap_est_mode::Union{Nothing,String} = nothing
end

Base.@kwdef mutable struct SCFT3WEndRow
  branch::Int
  terminal_node::Int
  u_rated::_ScfOptNum = nothing
  sn::_ScfOptNum = nothing
  role::Union{Nothing,String} = nothing
end

Base.@kwdef mutable struct SCFT3WNodeRef
  role::String
  node::Int
end

Base.@kwdef mutable struct SCFTransformer3WRow
  id::_ScfOptInt = nothing
  star_node::_ScfOptInt = nothing
  ends::Vector{SCFT3WEndRow} = SCFT3WEndRow[]
  leg_direction::Union{Nothing,String} = nothing
  tap_changer::_ScfOptInt = nothing
  # the CGMES-shaped INPUT shorthand: one PowerTransformerEnd-like object
  # per winding; the writer never emits this form, so it stays a document
  # object (validated by _scf_validate_references before conversion)
  nameplate::Union{Nothing,Dict{String,Any}} = nothing
  type::Union{Nothing,String} = nothing
  nodes::Vector{SCFT3WNodeRef} = SCFT3WNodeRef[]
  status::_ScfOptInt = nothing
end

Base.@kwdef mutable struct SCFShuntStateRow
  shunt::Int
  status::_ScfOptInt = nothing
  model::Union{Nothing,String} = nothing
  estimate::Union{Nothing,Bool} = nothing
end

# IEC 60909 source record: `kind` routes it, `node` is the writer's
# resolved file id, and the FIELD VOCABULARY deliberately stays with the
# record (net.sc_sources NamedTuples), so the free fields travel as a
# dictionary, sentinel strings ("inf"/"-inf") included.
Base.@kwdef mutable struct SCFScSourceRow
  kind::String
  node::_ScfOptInt = nothing
  fields::Dict{String,Any} = Dict{String,Any}()
end

Base.@kwdef mutable struct SCFComponents
  tap_changer::Vector{SCFTapChangerRow} = SCFTapChangerRow[]
  transformer3w::Vector{SCFTransformer3WRow} = SCFTransformer3WRow[]
  sc_source::Vector{SCFScSourceRow} = SCFScSourceRow[]
  shunt_state::Vector{SCFShuntStateRow} = SCFShuntStateRow[]
  # the declarative control.controllers schema, applied verbatim through
  # applyConfiguredControllers!; its vocabulary belongs to the control
  # framework, not to the case format
  controllers::Dict{String,Any} = Dict{String,Any}()
  # MATPOWER dcline records (stage 3a): the full per-row record of an
  # active dcline, so the run path restores the dcline metadata and HVDC
  # link registry the direct import kept on the net; free records like
  # sc_source, the vocabulary belongs to the importer
  matpower_dcline::Vector{Dict{String,Any}} = Dict{String,Any}[]
  # FOR001 contingency definitions a converted DTF case carries (stage 3a)
  for001_contingency::Vector{Dict{String,Any}} = Dict{String,Any}[]
end

Base.@kwdef mutable struct SCFStartNode
  vm_pu::_ScfOptNum = nothing
  va_deg::_ScfOptNum = nothing
end

Base.@kwdef mutable struct SCFStartState
  source::Union{Nothing,String} = nothing
  nodes::Dict{Int,SCFStartNode} = Dict{Int,SCFStartNode}()
end

"""
    SCFSparlectra

The namespaced `sparlectra` block of one SCF document. `roles`,
`components` and `start_state` are typed (the adapters construct them);
the genuinely free-form sub-objects (per-id `extra` records, controller
declarations, study definitions, measurement rows) stay dictionaries and
every consumer addresses them by key. `meta`, `roles` and `extra` are
always written; every other block only when non-empty, exactly like the
writer.
"""
Base.@kwdef mutable struct SCFSparlectra
  format_version::String = SCF_FORMAT_VERSION
  meta::Dict{String,Any} = Dict{String,Any}()
  roles::SCFRoles = SCFRoles()
  components::SCFComponents = SCFComponents()
  extra::Dict{String,Any} = Dict{String,Any}()
  measurements::Dict{String,Any} = Dict{String,Any}()
  start_state::SCFStartState = SCFStartState()
  contingencies::Dict{String,Any} = Dict{String,Any}()
  short_circuit::Dict{String,Any} = Dict{String,Any}()
  config::Dict{String,Any} = Dict{String,Any}()
  transformer_types::Dict{String,Any} = Dict{String,Any}()
  # the scenario block (scenario task D3); free document form, typed
  # accessors live in src/scenario/patch.jl
  scenarios::Dict{String,Any} = Dict{String,Any}()
end

"""
    SCFCase

One SCF document in memory (design decision D1): the PGM root scalars,
the typed `data` section, and the namespaced `sparlectra` block
(`nothing` for a strict-PGM file, which carries none). File reading
([`read_scf_json`](@ref)) and writing ([`write_scf_json`](@ref)) are
serialization of this struct; [`build_net`](@ref) is the one network
constructor over it.
"""
Base.@kwdef mutable struct SCFCase
  version::String = SCF_PGM_VERSION
  type::String = "input"
  is_batch::Bool = false
  attributes::Dict{String,Any} = Dict{String,Any}()
  data::SCFData = SCFData()
  sparlectra::Union{Nothing,SCFSparlectra} = SCFSparlectra()
end

## --- dict -> struct ---------------------------------------------------------
## The converters keep the reader's error vocabulary: every numeric check
## uses _scf_num/_scf_int with the same context strings scf_to_net used, so
## a broken field fails with the message the tests and the docs know.

_scf_opt_num(row, key::String, context::String)::_ScfOptNum = haskey(row, key) ? _scf_num(row[key], context) : nothing
_scf_opt_int(row, key::String, context::String)::_ScfOptInt = haskey(row, key) ? _scf_int(row[key], context) : nothing
_scf_tri_num(row, key::String, context::String)::_ScfTriNum = haskey(row, key) ? (row[key] === nothing ? nothing : _scf_num(row[key], context)) : missing
_scf_tri_str(row, key::String)::_ScfTriStr = haskey(row, key) ? (row[key] === nothing ? nothing : String(row[key])) : missing

_scfcase_node(row) = SCFNodeRow(id = _scf_int(_scf_require(row, "id", "a node entry"), "node.id"), u_rated = _scf_num(_scf_require(row, "u_rated", "a node entry"), "node.u_rated"))

function _scfcase_line(row)
  return SCFLineRow(
    id = _scf_int(row["id"], "line.id"),
    from_node = _scf_int(row["from_node"], "line.from_node"),
    to_node = _scf_int(row["to_node"], "to_node"),
    from_status = _scf_opt_int(row, "from_status", "line.from_status"),
    to_status = _scf_opt_int(row, "to_status", "line.to_status"),
    r1 = _scf_num(_scf_require(row, "r1", "a line entry"), "line.r1"),
    x1 = _scf_num(_scf_require(row, "x1", "a line entry"), "line.x1"),
    c1 = _scf_opt_num(row, "c1", "line.c1"),
    tan1 = _scf_opt_num(row, "tan1", "line.tan1"),
    g1 = _scf_opt_num(row, "g1", "line.g1"),
    i_n = _scf_opt_num(row, "i_n", "line.i_n"),
  )
end

function _scfcase_generic_branch(row)
  return SCFGenericBranchRow(
    id = _scf_int(row["id"], "generic_branch.id"),
    from_node = _scf_int(row["from_node"], "generic_branch.from_node"),
    to_node = _scf_int(row["to_node"], "to_node"),
    from_status = _scf_opt_int(row, "from_status", "generic_branch.from_status"),
    to_status = _scf_opt_int(row, "to_status", "generic_branch.to_status"),
    r1 = _scf_num(_scf_require(row, "r1", "a generic_branch entry"), "generic_branch.r1"),
    x1 = _scf_num(_scf_require(row, "x1", "a generic_branch entry"), "generic_branch.x1"),
    g1 = _scf_opt_num(row, "g1", "generic_branch.g1"),
    b1 = _scf_opt_num(row, "b1", "generic_branch.b1"),
    k = _scf_opt_num(row, "k", "generic_branch.k"),
    theta = _scf_opt_num(row, "theta", "generic_branch.theta"),
    sn = _scf_opt_num(row, "sn", "generic_branch.sn"),
  )
end

function _scfcase_link(row)
  return SCFLinkRow(
    id = _scf_int(row["id"], "link.id"),
    from_node = _scf_int(row["from_node"], "link.from_node"),
    to_node = _scf_int(row["to_node"], "link.to_node"),
    from_status = _scf_opt_int(row, "from_status", "link.from_status"),
    to_status = _scf_opt_int(row, "to_status", "link.to_status"),
  )
end

function _scfcase_source(row)
  return SCFSourceRow(
    id = _scf_int(row["id"], "source.id"),
    node = _scf_int(_scf_require(row, "node", "a source entry"), "source.node"),
    status = _scf_opt_int(row, "status", "source.status"),
    u_ref = _scf_opt_num(row, "u_ref", "source.u_ref"),
    u_ref_angle = _scf_opt_num(row, "u_ref_angle", "source.u_ref_angle"),
    sk = _scf_opt_num(row, "sk", "source.sk"),
    rx_ratio = _scf_opt_num(row, "rx_ratio", "source.rx_ratio"),
    z01_ratio = _scf_opt_num(row, "z01_ratio", "source.z01_ratio"),
  )
end

function _scfcase_appliance(row, kind::String)
  return SCFApplianceRow(
    id = _scf_int(row["id"], "$(kind).id"),
    node = _scf_int(row["node"], "$(kind).node"),
    status = _scf_opt_int(row, "status", "$(kind).status"),
    type = _scf_opt_int(row, "type", "$(kind).type"),
    p_specified = _scf_opt_num(row, "p_specified", "$(kind).p_specified"),
    q_specified = _scf_opt_num(row, "q_specified", "$(kind).q_specified"),
    u_ref = _scf_opt_num(row, "u_ref", "sym_gen.u_ref"),
  )
end

function _scfcase_shunt(row)
  return SCFShuntRow(
    id = _scf_int(row["id"], "shunt.id"),
    node = _scf_int(row["node"], "shunt.node"),
    status = _scf_opt_int(row, "status", "shunt.status"),
    g1 = _scf_opt_num(row, "g1", "shunt.g1"),
    b1 = _scf_opt_num(row, "b1", "shunt.b1"),
  )
end

function _scfcase_regulator(row)
  return SCFVoltageRegulatorRow(
    id = _scf_int(row["id"], "voltage_regulator.id"),
    regulated_object = _scf_int(row["regulated_object"], "voltage_regulator.regulated_object"),
    status = _scf_opt_int(row, "status", "voltage_regulator.status"),
    u_ref = _scf_opt_num(row, "u_ref", "voltage_regulator.u_ref"),
    q_min = _scf_opt_num(row, "q_min", "voltage_regulator.q_min"),
    q_max = _scf_opt_num(row, "q_max", "voltage_regulator.q_max"),
  )
end

function _scfcase_voltage_sensor(row)
  return SCFVoltageSensorRow(
    id = _scf_int(row["id"], "sym_voltage_sensor.id"),
    measured_object = _scf_int(row["measured_object"], "sym_voltage_sensor.measured_object"),
    u_measured = _scf_tri_num(row, "u_measured", "u_measured"),
    u_sigma = _scf_tri_num(row, "u_sigma", "u_sigma"),
    u_angle_measured = _scf_tri_num(row, "u_angle_measured", "u_angle_measured"),
  )
end

function _scfcase_power_sensor(row)
  return SCFPowerSensorRow(
    id = _scf_int(row["id"], "sym_power_sensor.id"),
    measured_object = _scf_int(row["measured_object"], "sym_power_sensor.measured_object"),
    measured_terminal_type = _scf_tri_str(row, "measured_terminal_type"),
    p_measured = _scf_tri_num(row, "p_measured", "p_measured"),
    q_measured = _scf_tri_num(row, "q_measured", "q_measured"),
    power_sigma = _scf_tri_num(row, "power_sigma", "power_sigma"),
    p_sigma = _scf_tri_num(row, "p_sigma", "p_measured"),
    q_sigma = _scf_tri_num(row, "q_sigma", "q_measured"),
  )
end

function _scfcase_current_sensor(row)
  return SCFCurrentSensorRow(
    id = _scf_int(row["id"], "sym_current_sensor.id"),
    measured_object = _scf_int(row["measured_object"], "sym_current_sensor.measured_object"),
    measured_terminal_type = _scf_tri_str(row, "measured_terminal_type"),
    angle_measurement_type = _scf_tri_str(row, "angle_measurement_type"),
    i_measured = _scf_tri_num(row, "i_measured", "i_measured"),
    i_sigma = _scf_tri_num(row, "i_sigma", "i_sigma"),
    i_angle_measured = _scf_tri_num(row, "i_angle_measured", "i_angle_measured"),
    i_angle_sigma = _scf_tri_num(row, "i_angle_sigma", "i_angle_sigma"),
  )
end

function _scfcase_fault(row)
  return SCFFaultRow(
    id = _scf_int(_scf_get(row, "id", 0), "fault.id"),
    fault_object = _scf_int(_scf_require(row, "fault_object", "a fault entry"), "fault.fault_object"),
    status = _scf_opt_int(row, "status", "fault.status"),
    fault_type = haskey(row, "fault_type") ? String(row["fault_type"]) : nothing,
    r_f = _scf_opt_num(row, "r_f", "fault.r_f"),
    x_f = _scf_opt_num(row, "x_f", "fault.x_f"),
  )
end

_scfcase_dict(v) = v isa AbstractDict ? Dict{String,Any}(String(k) => x for (k, x) in v) : Dict{String,Any}()

_scfcase_opt_str(row, key::String) = haskey(row, key) ? (row[key] === nothing ? nothing : String(row[key])) : nothing
_scfcase_opt_bool(row, key::String) = haskey(row, key) ? (row[key] === nothing ? nothing : row[key] === true) : nothing

function _scfcase_roles(raw)::SCFRoles
  raw isa AbstractDict || return SCFRoles()
  slack_raw = _scf_get(raw, "slack", nothing)
  slack = if slack_raw isa AbstractDict
    SCFSlackRole(
      mode = String(_scf_get(slack_raw, "mode", "single")),
      nodes = Int[_scf_int(x, "roles.slack.nodes") for x in _scf_get(slack_raw, "nodes", [])],
      participation = SCFParticipationEntry[SCFParticipationEntry(object = _scf_int(_scf_require(p, "object", "a participation entry"), "participation.object"), factor = _scf_num(_scf_require(p, "factor", "a participation entry"), "participation.factor")) for p in _scf_get(slack_raw, "participation", [])],
    )
  else
    nothing
  end
  return SCFRoles(
    slack = slack,
    aux_nodes = Int[_scf_int(x, "roles.aux_nodes") for x in _scf_get(raw, "aux_nodes", [])],
    isolated_nodes = Int[_scf_int(x, "roles.isolated_nodes") for x in _scf_get(raw, "isolated_nodes", [])],
  )
end

_scf_opt_gridpos(row, key::String, context::String)::_ScfOptGridPos = haskey(row, key) ? (row[key] isa Integer ? Int(row[key]) : _scf_num(row[key], context)) : nothing

function _scfcase_tap_controller(row)::SCFTapControllerRow
  return SCFTapControllerRow(
    index = _scf_int(_scf_require(row, "index", "a tap_changer controller"), "controller.index"),
    alpha_deg = _scf_opt_num(row, "alpha_deg", "controller.alpha_deg"),
    step = _scf_opt_num(row, "step", "controller.step"),
    step_deg = _scf_opt_num(row, "step_deg", "controller.step_deg"),
    pos = _scf_opt_gridpos(row, "pos", "controller.pos"),
    pos_min = _scf_opt_gridpos(row, "pos_min", "controller.pos_min"),
    pos_max = _scf_opt_gridpos(row, "pos_max", "controller.pos_max"),
  )
end

function _scfcase_tap_changer(row)::SCFTapChangerRow
  return SCFTapChangerRow(
    id = _scf_int(_scf_get(row, "id", 0), "tap_changer.id"),
    branch = _scf_int(_scf_require(row, "branch", "a tap_changer entry"), "tap_changer.branch"),
    side = _scfcase_opt_str(row, "side"),
    ratio_base = _scf_opt_num(row, "ratio_base", "tap_changer.ratio_base"),
    angle_base_deg = _scf_opt_num(row, "angle_base_deg", "tap_changer.angle_base_deg"),
    controllers = SCFTapControllerRow[_scfcase_tap_controller(c) for c in _scf_get(row, "controllers", [])],
    control = _scfcase_dict(_scf_get(row, "control", nothing)),
    tap_est_mode = _scfcase_opt_str(row, "tap_est_mode"),
  )
end

function _scfcase_transformer3w(row)::SCFTransformer3WRow
  plate = _scf_get(row, "nameplate", nothing)
  return SCFTransformer3WRow(
    id = _scf_opt_int(row, "id", "transformer3w.id"),
    star_node = _scf_opt_int(row, "star_node", "transformer3w.star_node"),
    ends = SCFT3WEndRow[SCFT3WEndRow(branch = _scf_int(_scf_require(e, "branch", "a transformer3w end"), "transformer3w end.branch"), terminal_node = _scf_int(_scf_require(e, "terminal_node", "a transformer3w end"), "transformer3w end.terminal_node"), u_rated = _scf_opt_num(e, "u_rated", "transformer3w end.u_rated"), sn = _scf_opt_num(e, "sn", "transformer3w end.sn"), role = _scfcase_opt_str(e, "role")) for e in (haskey(row, "nameplate") || haskey(row, "type") ? Any[] : _scf_get(row, "ends", []))],
    leg_direction = _scfcase_opt_str(row, "leg_direction"),
    tap_changer = _scf_opt_int(row, "tap_changer", "transformer3w.tap_changer"),
    nameplate = plate isa AbstractDict ? _scfcase_dict(plate) : nothing,
    type = _scfcase_opt_str(row, "type"),
    nodes = SCFT3WNodeRef[SCFT3WNodeRef(role = String(_scf_require(n, "role", "a transformer3w node reference")), node = _scf_int(_scf_require(n, "node", "a transformer3w node reference"), "transformer3w.nodes.node")) for n in _scf_get(row, "nodes", [])],
    status = _scf_opt_int(row, "status", "transformer3w.status"),
  )
end

function _scfcase_sc_source(row)::SCFScSourceRow
  fields = Dict{String,Any}()
  for (k, v) in row
    key = String(k)
    key in ("kind", "node") && continue
    fields[key] = v
  end
  return SCFScSourceRow(kind = String(_scf_get(row, "kind", "")), node = _scf_opt_int(row, "node", "sc_source.node"), fields = fields)
end

function _scfcase_shunt_state(row)::SCFShuntStateRow
  return SCFShuntStateRow(
    shunt = _scf_int(_scf_require(row, "shunt", "a shunt_state entry"), "shunt_state.shunt"),
    status = _scf_opt_int(row, "status", "shunt_state.status"),
    model = _scfcase_opt_str(row, "model"),
    estimate = _scfcase_opt_bool(row, "estimate"),
  )
end

function _scfcase_components(raw)::SCFComponents
  raw isa AbstractDict || return SCFComponents()
  return SCFComponents(
    tap_changer = SCFTapChangerRow[_scfcase_tap_changer(r) for r in _scf_get(raw, "tap_changer", [])],
    transformer3w = SCFTransformer3WRow[_scfcase_transformer3w(r) for r in _scf_get(raw, "transformer3w", []) if r isa AbstractDict],
    sc_source = SCFScSourceRow[_scfcase_sc_source(r) for r in _scf_get(raw, "sc_source", []) if r isa AbstractDict],
    shunt_state = SCFShuntStateRow[_scfcase_shunt_state(r) for r in _scf_get(raw, "shunt_state", []) if r isa AbstractDict],
    controllers = _scfcase_dict(_scf_get(raw, "controllers", nothing)),
    matpower_dcline = Dict{String,Any}[_scfcase_dict(r) for r in _scf_get(raw, "matpower_dcline", []) if r isa AbstractDict],
    for001_contingency = Dict{String,Any}[_scfcase_dict(r) for r in _scf_get(raw, "for001_contingency", []) if r isa AbstractDict],
  )
end

function _scfcase_start_state(raw)::SCFStartState
  raw isa AbstractDict || return SCFStartState()
  nodes = Dict{Int,SCFStartNode}()
  for (idstr, v) in _scf_get(raw, "nodes", Dict{String,Any}())
    id = tryparse(Int, String(idstr))
    id === nothing && continue
    v isa AbstractDict || continue
    nodes[id] = SCFStartNode(vm_pu = _scf_opt_num(v, "vm_pu", "start_state.vm_pu"), va_deg = _scf_opt_num(v, "va_deg", "start_state.va_deg"))
  end
  return SCFStartState(source = _scfcase_opt_str(raw, "source"), nodes = nodes)
end

"""
    scfcase_from_root(root::AbstractDict) -> SCFCase

Convert a parsed (and validated) SCF root object into its typed form. The
converter mirrors the reader's field vocabulary exactly; unknown row
fields of a foreign PGM file are dropped here, which is the documented
scope of the round-trip contract (it covers files Sparlectra wrote).
"""
function scfcase_from_root(root::AbstractDict)::SCFCase
  data_raw = _scf_get(root, "data", Dict{String,Any}())
  data = SCFData(
    node = SCFNodeRow[_scfcase_node(r) for r in _scf_get(data_raw, "node", [])],
    line = SCFLineRow[_scfcase_line(r) for r in _scf_get(data_raw, "line", [])],
    generic_branch = SCFGenericBranchRow[_scfcase_generic_branch(r) for r in _scf_get(data_raw, "generic_branch", [])],
    link = SCFLinkRow[_scfcase_link(r) for r in _scf_get(data_raw, "link", [])],
    source = SCFSourceRow[_scfcase_source(r) for r in _scf_get(data_raw, "source", [])],
    sym_load = SCFApplianceRow[_scfcase_appliance(r, "sym_load") for r in _scf_get(data_raw, "sym_load", [])],
    sym_gen = SCFApplianceRow[_scfcase_appliance(r, "sym_gen") for r in _scf_get(data_raw, "sym_gen", [])],
    shunt = SCFShuntRow[_scfcase_shunt(r) for r in _scf_get(data_raw, "shunt", [])],
    voltage_regulator = SCFVoltageRegulatorRow[_scfcase_regulator(r) for r in _scf_get(data_raw, "voltage_regulator", [])],
    sym_voltage_sensor = SCFVoltageSensorRow[_scfcase_voltage_sensor(r) for r in _scf_get(data_raw, "sym_voltage_sensor", [])],
    sym_power_sensor = SCFPowerSensorRow[_scfcase_power_sensor(r) for r in _scf_get(data_raw, "sym_power_sensor", [])],
    sym_current_sensor = SCFCurrentSensorRow[_scfcase_current_sensor(r) for r in _scf_get(data_raw, "sym_current_sensor", [])],
    fault = SCFFaultRow[_scfcase_fault(r) for r in _scf_get(data_raw, "fault", [])],
  )
  spar_raw = _scf_get(root, "sparlectra", nothing)
  spar = if spar_raw isa AbstractDict
    SCFSparlectra(
      format_version = String(_scf_get(spar_raw, "format_version", SCF_FORMAT_VERSION)),
      meta = _scfcase_dict(_scf_get(spar_raw, "meta", nothing)),
      roles = _scfcase_roles(_scf_get(spar_raw, "roles", nothing)),
      components = _scfcase_components(_scf_get(spar_raw, "components", nothing)),
      extra = _scfcase_dict(_scf_get(spar_raw, "extra", nothing)),
      measurements = _scfcase_dict(_scf_get(spar_raw, "measurements", nothing)),
      start_state = _scfcase_start_state(_scf_get(spar_raw, "start_state", nothing)),
      contingencies = _scfcase_dict(_scf_get(spar_raw, "contingencies", nothing)),
      short_circuit = _scfcase_dict(_scf_get(spar_raw, "short_circuit", nothing)),
      config = _scfcase_dict(_scf_get(spar_raw, "config", nothing)),
      transformer_types = _scfcase_dict(_scf_get(spar_raw, "transformer_types", nothing)),
      scenarios = _scfcase_dict(_scf_get(spar_raw, "scenarios", nothing)),
    )
  else
    nothing
  end
  version_raw = _scf_get(root, "version", SCF_PGM_VERSION)
  return SCFCase(version = string(version_raw), type = String(_scf_get(root, "type", "input")), is_batch = _scf_get(root, "is_batch", false) === true, attributes = _scfcase_dict(_scf_get(root, "attributes", nothing)), data = data, sparlectra = spar)
end

## --- struct -> dict ---------------------------------------------------------
## Emission writes exactly the keys the struct holds: a plain optional is
## written when it is not `nothing`, a tri-state field is written unless it
## is `missing` (with `nothing` as JSON null). The serializer sorts object
## keys, so the bytes do not depend on emission order.

function _scfrow_dict(row)::Dict{String,Any}
  d = Dict{String,Any}()
  for key in propertynames(row)
    v = getproperty(row, key)
    v === nothing && continue
    d[String(key)] = v
  end
  return d
end

# the tri-state structs need the null-emitting variant: `nothing` IS data
function _scfrow_dict_tri(row)::Dict{String,Any}
  d = Dict{String,Any}()
  for key in propertynames(row)
    v = getproperty(row, key)
    v === missing && continue
    d[String(key)] = v
  end
  return d
end

_scfcase_rows(rows::AbstractVector, tri::Bool) = Any[tri ? _scfrow_dict_tri(r) : _scfrow_dict(r) for r in rows]

function _scfdict_roles(roles::SCFRoles)::Dict{String,Any}
  d = Dict{String,Any}()
  slack = roles.slack
  if slack !== nothing
    sdict = Dict{String,Any}("mode" => slack.mode)
    isempty(slack.nodes) || (sdict["nodes"] = copy(slack.nodes))
    isempty(slack.participation) || (sdict["participation"] = Any[Dict{String,Any}("object" => p.object, "factor" => p.factor) for p in slack.participation])
    d["slack"] = sdict
  end
  isempty(roles.aux_nodes) || (d["aux_nodes"] = copy(roles.aux_nodes))
  isempty(roles.isolated_nodes) || (d["isolated_nodes"] = copy(roles.isolated_nodes))
  return d
end

function _scfdict_tap_changer(row::SCFTapChangerRow)::Dict{String,Any}
  d = Dict{String,Any}("id" => row.id, "branch" => row.branch)
  row.side === nothing || (d["side"] = row.side)
  row.ratio_base === nothing || (d["ratio_base"] = row.ratio_base)
  row.angle_base_deg === nothing || (d["angle_base_deg"] = row.angle_base_deg)
  isempty(row.controllers) || (d["controllers"] = Any[_scfrow_dict(c) for c in row.controllers])
  isempty(row.control) || (d["control"] = row.control)
  row.tap_est_mode === nothing || (d["tap_est_mode"] = row.tap_est_mode)
  return d
end

function _scfdict_transformer3w(row::SCFTransformer3WRow)::Dict{String,Any}
  d = Dict{String,Any}()
  row.id === nothing || (d["id"] = row.id)
  row.star_node === nothing || (d["star_node"] = row.star_node)
  isempty(row.ends) || (d["ends"] = Any[_scfrow_dict(e) for e in row.ends])
  row.leg_direction === nothing || (d["leg_direction"] = row.leg_direction)
  row.tap_changer === nothing || (d["tap_changer"] = row.tap_changer)
  row.nameplate === nothing || (d["nameplate"] = row.nameplate)
  row.type === nothing || (d["type"] = row.type)
  isempty(row.nodes) || (d["nodes"] = Any[Dict{String,Any}("role" => n.role, "node" => n.node) for n in row.nodes])
  row.status === nothing || (d["status"] = row.status)
  return d
end

function _scfdict_sc_source(row::SCFScSourceRow)::Dict{String,Any}
  d = Dict{String,Any}(row.fields)
  d["kind"] = row.kind
  row.node === nothing || (d["node"] = row.node)
  return d
end

function _scfdict_components(comps::SCFComponents)::Dict{String,Any}
  d = Dict{String,Any}()
  isempty(comps.tap_changer) || (d["tap_changer"] = Any[_scfdict_tap_changer(r) for r in comps.tap_changer])
  isempty(comps.transformer3w) || (d["transformer3w"] = Any[_scfdict_transformer3w(r) for r in comps.transformer3w])
  isempty(comps.sc_source) || (d["sc_source"] = Any[_scfdict_sc_source(r) for r in comps.sc_source])
  isempty(comps.shunt_state) || (d["shunt_state"] = Any[_scfrow_dict(r) for r in comps.shunt_state])
  isempty(comps.controllers) || (d["controllers"] = comps.controllers)
  isempty(comps.matpower_dcline) || (d["matpower_dcline"] = Any[r for r in comps.matpower_dcline])
  isempty(comps.for001_contingency) || (d["for001_contingency"] = Any[r for r in comps.for001_contingency])
  return d
end

function _scfdict_start_state(start::SCFStartState)::Dict{String,Any}
  d = Dict{String,Any}()
  start.source === nothing || (d["source"] = start.source)
  if !isempty(start.nodes)
    nodes = Dict{String,Any}()
    for (id, n) in start.nodes
      nd = Dict{String,Any}()
      n.vm_pu === nothing || (nd["vm_pu"] = n.vm_pu)
      n.va_deg === nothing || (nd["va_deg"] = n.va_deg)
      nodes[string(id)] = nd
    end
    d["nodes"] = nodes
  end
  return d
end

"""
    scf_root_dict(case::SCFCase) -> Dict{String,Any}

The document form of a typed case, ready for [`scf_json_string`](@ref).
Component keys appear exactly when their vector is non-empty; the
namespaced block writes `meta`, `roles` and `extra` always and every
other sub-block only when non-empty, matching the writer's shape.
"""
function scf_root_dict(case::SCFCase)::Dict{String,Any}
  data = Dict{String,Any}()
  for (key, rows, tri) in (
    ("node", case.data.node, false),
    ("line", case.data.line, false),
    ("generic_branch", case.data.generic_branch, false),
    ("link", case.data.link, false),
    ("source", case.data.source, false),
    ("sym_load", case.data.sym_load, false),
    ("sym_gen", case.data.sym_gen, false),
    ("shunt", case.data.shunt, false),
    ("voltage_regulator", case.data.voltage_regulator, false),
    ("sym_voltage_sensor", case.data.sym_voltage_sensor, true),
    ("sym_power_sensor", case.data.sym_power_sensor, true),
    ("sym_current_sensor", case.data.sym_current_sensor, true),
    ("fault", case.data.fault, false),
  )
    isempty(rows) || (data[key] = _scfcase_rows(rows, tri))
  end
  root = Dict{String,Any}("version" => case.version, "type" => case.type, "is_batch" => case.is_batch, "attributes" => case.attributes, "data" => data)
  spar = case.sparlectra
  if spar !== nothing
    sd = Dict{String,Any}("format_version" => spar.format_version, "meta" => spar.meta, "roles" => _scfdict_roles(spar.roles), "extra" => spar.extra)
    components = _scfdict_components(spar.components)
    isempty(components) || (sd["components"] = components)
    start_state = _scfdict_start_state(spar.start_state)
    isempty(start_state) || (sd["start_state"] = start_state)
    for (key, block) in (("measurements", spar.measurements), ("contingencies", spar.contingencies), ("short_circuit", spar.short_circuit), ("config", spar.config), ("transformer_types", spar.transformer_types), ("scenarios", spar.scenarios))
      isempty(block) || (sd[key] = block)
    end
    root["sparlectra"] = sd
  end
  return root
end

## --- tagged value encoding ---------------------------------------------------
## Importer metadata records (the DTF branch metadata, the typed phase-tap
## model) carry Symbols, Chars, complex numbers and nested typed rows that
## plain JSON cannot spell. The tagged encoding keeps every value round-trip
## exact and ORDER-preserving for named tuples, so a consumer reads back the
## record the importer stored, field for field.

function scf_encode_value(v)
  v === nothing && return nothing
  v isa Bool && return v
  v isa Symbol && return Dict{String,Any}("__sym" => String(v))
  v isa Char && return Dict{String,Any}("__char" => string(v))
  v isa Complex && return Dict{String,Any}("__cplx" => Any[Float64(real(v)), Float64(imag(v))])
  v isa AbstractFloat && return isfinite(v) ? Float64(v) : Dict{String,Any}("__f" => v > 0 ? "inf" : "-inf")
  v isa Integer && return Int(v)
  v isa AbstractString && return String(v)
  v isa TapTablePoint && return Dict{String,Any}("__ttp" => Dict{String,Any}("step" => v.step, "ratio" => v.ratio, "angle_deg" => v.angle_deg, "x_pu" => v.x_pu))
  v isa PhaseTapChangerModel && return Dict{String,Any}("__ptc" => Dict{String,Any}(
    "kind" => String(v.kind), "step" => v.step, "lowStep" => v.lowStep, "highStep" => v.highStep, "neutralStep" => v.neutralStep,
    "voltage_step_increment" => v.voltage_step_increment, "step_phase_shift_increment" => v.step_phase_shift_increment,
    "winding_connection_angle_deg" => v.winding_connection_angle_deg, "x_min" => v.x_min, "x_max" => v.x_max,
    "convention" => String(v.convention), "table" => v.table === nothing ? nothing : Any[scf_encode_value(p) for p in v.table]))
  v isa NamedTuple && return Dict{String,Any}("__nt" => Any[Any[String(k), scf_encode_value(getproperty(v, k))] for k in propertynames(v)])
  v isa AbstractVector && return Any[scf_encode_value(x) for x in v]
  throw(ArgumentError("SCF: no tagged encoding for a value of type $(typeof(v))."))
end

function scf_decode_value(v)
  v === nothing && return nothing
  if v isa AbstractDict
    haskey(v, "__sym") && return Symbol(String(v["__sym"]))
    haskey(v, "__char") && return only(String(v["__char"]))
    haskey(v, "__cplx") && return ComplexF64(Float64(v["__cplx"][1]), Float64(v["__cplx"][2]))
    haskey(v, "__f") && return scf_number_or_sentinel(v["__f"], "tagged float")
    if haskey(v, "__ttp")
      p = v["__ttp"]
      return TapTablePoint(step = Int(p["step"]), ratio = Float64(p["ratio"]), angle_deg = Float64(p["angle_deg"]), x_pu = p["x_pu"] === nothing ? nothing : Float64(p["x_pu"]))
    end
    if haskey(v, "__ptc")
      p = v["__ptc"]
      return PhaseTapChangerModel(
        kind = Symbol(String(p["kind"])), step = Int(p["step"]), lowStep = Int(p["lowStep"]), highStep = Int(p["highStep"]), neutralStep = Int(p["neutralStep"]),
        voltage_step_increment = p["voltage_step_increment"] === nothing ? nothing : Float64(p["voltage_step_increment"]),
        step_phase_shift_increment = p["step_phase_shift_increment"] === nothing ? nothing : Float64(p["step_phase_shift_increment"]),
        winding_connection_angle_deg = p["winding_connection_angle_deg"] === nothing ? nothing : Float64(p["winding_connection_angle_deg"]),
        x_min = p["x_min"] === nothing ? nothing : Float64(p["x_min"]),
        x_max = p["x_max"] === nothing ? nothing : Float64(p["x_max"]),
        convention = Symbol(String(p["convention"])),
        table = p["table"] === nothing ? nothing : TapTablePoint[scf_decode_value(x) for x in p["table"]])
    end
    if haskey(v, "__nt")
      pairs = v["__nt"]
      keys_ = Symbol[Symbol(String(p[1])) for p in pairs]
      vals = Any[scf_decode_value(p[2]) for p in pairs]
      return NamedTuple{Tuple(keys_)}(Tuple(vals))
    end
    return Dict{String,Any}(String(k) => scf_decode_value(x) for (k, x) in v)
  end
  v isa AbstractVector && return Any[scf_decode_value(x) for x in v]
  return v
end

## --- file IO ----------------------------------------------------------------

"""
    read_scf_json(path) -> SCFCase

Read one SCF case file into its typed in-memory form. Validation is the
reader's staged file-level check (schema, study definitions, reference
integrity) on the parsed document, so a broken file fails here with the
same line of reasoning the upload path shows; the model plausibility
stage runs when the network is built ([`build_net`](@ref)).
"""
function read_scf_json(path::AbstractString)::SCFCase
  isfile(path) || throw(ArgumentError("SCF: case file not found: $(path)"))
  root = scf_json_parse(read(String(path), String))
  _scf_validate_schema(root)
  _scf_validate_studies(root)
  _scf_validate_references(root)
  return scfcase_from_root(root)
end

"""
    write_scf_json(case, path) -> String

Serialize a typed case to its canonical file form (deterministic bytes,
see [`scf_json_string`](@ref)) and return the path.
"""
function write_scf_json(case::SCFCase, path::AbstractString)::String
  p = String(path)
  mkpath(dirname(abspath(p)))
  open(p, "w") do io
    print(io, scf_json_string(scf_root_dict(case)))
  end
  return p
end
