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

module DTFImporter

using ..Sparlectra: Net, addBus!, addProsumer!, addShuntMatpower!, _addPIModelACLine_by_idx!, _addPIModelTrafo_by_idx!, geNetBusIdx, validate!, normalize_bus_shunt_model

export DTFCase, DTFParams, DTFSize, DTFBranch, DTFBus, DTFCompensation, DTFTransformerControl, DTFOutage, read_dtf, build_net

"""Native DTF/FOR001 importer parameter card values preserved for auditing."""
struct DTFParams
  raw::String
  values::Vector{Float64}
end

struct DTFSize
  raw::String
  NGES::Int
  LGES::Int
  IKOMP::Int
  IRETRA::Int
  slack::String
end

struct DTFBranch
  raw::String
  index::Int
  kind::Char
  voltage_level_index::Int
  parallel_id::String
  from::String
  to::String
  r_ohm::Float64
  x_ohm::Float64
  g_s::Float64
  b_s::Float64
  imax_ka::Union{Nothing,Float64}
end

struct DTFBus
  raw::String
  index::Int
  bus_type::Int
  voltage_level_index::Int
  name::String
  start_kv::Float64
  angle_deg::Float64
  pd_mw::Float64
  qd_mvar::Float64
  pg_mw::Float64
  qg_mvar::Float64
  qmin_mvar::Union{Nothing,Float64}
  qmax_mvar::Union{Nothing,Float64}
end

struct DTFCompensation
  raw::String
  index::Int
  fields::Vector{String}
end

struct DTFTransformerControl
  raw::String
  index::Int
  regulated_side::String
  unregulated_side::String
  parallel_id::String
  phase_shifter_flag::String
  from::String
  to::String
  nominal_unregulated_kv::Union{Nothing,Float64}
  nominal_regulated_kv::Union{Nothing,Float64}
  longitudinal_range_percent::Union{Nothing,Float64}
  added_voltage_angle_deg::Union{Nothing,Float64}
  max_tap_step::Union{Nothing,Int}
  actual_tap_step::Union{Nothing,Int}
  quadrature_range_percent::Union{Nothing,Float64}
  quadrature_steps::Union{Nothing,Int}
end

struct DTFOutage
  raw::String
  index::Int
  kind::Char
  voltage_level_index::Int
  parallel_id::String
  from::String
  to::String
end

struct DTFCase
  source_path::String
  baseMVA::Float64
  params::DTFParams
  texts::Vector{String}
  nominal_voltages_kv::Vector{Float64}
  size::DTFSize
  branches::Vector{DTFBranch}
  compensations::Vector{DTFCompensation}
  transformer_controls::Vector{DTFTransformerControl}
  buses::Vector{DTFBus}
  outages::Vector{DTFOutage}
end

_slice(s::AbstractString, a::Int, b::Int) = a > lastindex(s) ? "" : s[a:min(b, lastindex(s))]
_parse_float(s::AbstractString) = isempty(strip(s)) ? nothing : parse(Float64, replace(strip(s), 'D' => 'E', 'd' => 'E'))
_parse_int(s::AbstractString) = (x = _parse_float(s); x === nothing ? nothing : Int(round(x)))
_numbers(s::AbstractString) = [parse(Float64, replace(m.match, 'D' => 'E', 'd' => 'E')) for m in eachmatch(r"[-+]?\d*\.?\d+(?:[EeDd][-+]?\d+)?", s)]

function _parse_branch(line::String, index::Int)::DTFBranch
  kind = isempty(line) || line[1] == ' ' ? 'L' : line[1]
  voltage_level_index = line[2] == ' ' ? 1 : parse(Int, string(line[2]))
  parallel_id = strip(_slice(line, 3, 3))
  from = String(strip(_slice(line, 6, 13)))
  to = String(strip(_slice(line, 16, 23)))
  vals = _numbers(_slice(line, 24, lastindex(line)))
  length(vals) >= 4 || error("DTF branch card $(index) has too few numeric fields: $(line)")
  imax = length(vals) >= 5 ? vals[5] : nothing
  return DTFBranch(line, index, kind, voltage_level_index, parallel_id, from, to, vals[1], vals[2], vals[3], vals[4], imax)
end

function _parse_outage(line::String, index::Int)::DTFOutage
  kind = isempty(line) || line[1] == ' ' ? 'L' : line[1]
  voltage_level_index = line[2] == ' ' ? 1 : parse(Int, string(line[2]))
  parallel_id = strip(_slice(line, 3, 3))
  from = String(strip(_slice(line, 6, 13)))
  to = String(strip(_slice(line, 16, 23)))
  return DTFOutage(line, index, kind, voltage_level_index, parallel_id, from, to)
end

function _parse_bus(line::String, index::Int)::DTFBus
  bus_type = line[1] == ' ' ? 0 : parse(Int, string(line[1]))
  voltage_level_index = line[2] == ' ' ? 1 : parse(Int, string(line[2]))
  name = String(strip(_slice(line, 6, 13)))
  vals = _numbers(_slice(line, 14, lastindex(line)))
  length(vals) >= 6 || error("DTF bus card $(index) has too few numeric fields: $(line)")
  qmin = length(vals) >= 7 ? vals[7] : nothing
  qmax = length(vals) >= 8 ? vals[8] : nothing
  return DTFBus(line, index, bus_type, voltage_level_index, name, vals[1], vals[2], vals[3], vals[4], vals[5], vals[6], qmin, qmax)
end

function _parse_transformer_control(line::String, index::Int)::DTFTransformerControl
  parallel_id = String(strip(_slice(line, 3, 3)))
  from = String(strip(_slice(line, 6, 13)))
  to = String(strip(_slice(line, 16, 23)))
  vals = _numbers(_slice(line, 24, lastindex(line)))
  return DTFTransformerControl(line, index, strip(_slice(line, 1, 1)), strip(_slice(line, 2, 2)), parallel_id, "", from, to,
    length(vals) >= 1 ? vals[1] : nothing,
    length(vals) >= 2 ? vals[2] : nothing,
    length(vals) >= 3 ? vals[3] : nothing,
    length(vals) >= 4 ? vals[4] : nothing,
    length(vals) >= 5 ? Int(round(vals[5])) : nothing,
    length(vals) >= 6 ? Int(round(vals[6])) : nothing,
    length(vals) >= 7 ? vals[7] : nothing,
    length(vals) >= 8 ? Int(round(vals[8])) : nothing)
end

function _find_control(branch::DTFBranch, controls::Vector{DTFTransformerControl})
  for c in controls
    if c.from == branch.from && c.to == branch.to && c.parallel_id == branch.parallel_id
      return c
    end
  end
  return nothing
end

"""
    read_dtf(path; baseMVA = 100.0, strict = true) -> DTFCase

Read a legacy DTF/FOR001 fixed-column input file into a typed `DTFCase`.
This MVP parses outage cards into `DTFOutage` records but intentionally does
not execute outage simulations.
"""
function read_dtf(path; baseMVA::Real = 100.0, strict::Bool = true)::DTFCase
  lines = readlines(path)
  length(lines) >= 7 || error("DTF/FOR001 file is too short: $(path)")
  params = DTFParams(lines[1], _numbers(lines[1]))
  texts = String.(lines[2:5])
  nominal = _numbers(lines[6])
  size_vals = _numbers(lines[7])
  length(size_vals) >= 4 || error("DTF/FOR001 network size card is incomplete")
  size_tokens = split(strip(lines[7]))
  size = DTFSize(lines[7], Int(size_vals[1]), Int(size_vals[2]), Int(size_vals[3]), Int(size_vals[4]), length(size_tokens) >= 5 ? String(size_tokens[5]) : "")
  i = 8
  branches = DTFBranch[]
  for idx in 1:size.LGES
    push!(branches, _parse_branch(lines[i], idx)); i += 1
  end
  comps = DTFCompensation[]
  for idx in 1:size.IKOMP
    push!(comps, DTFCompensation(lines[i], idx, split(strip(lines[i])))); i += 1
  end
  controls = DTFTransformerControl[]
  for idx in 1:size.IRETRA
    push!(controls, _parse_transformer_control(lines[i], idx)); i += 1
  end
  buses = DTFBus[]
  for idx in 1:size.NGES
    push!(buses, _parse_bus(lines[i], idx)); i += 1
  end
  outages = DTFOutage[]
  in_outage = false
  outage_index = 0
  while i <= length(lines)
    line = lines[i]
    token = uppercase(strip(line))
    if token == "AUSFALL"
      in_outage = true
    elseif token == "ENDE"
      in_outage = false
    elseif in_outage && !isempty(strip(line))
      outage_index += 1
      push!(outages, _parse_outage(line, outage_index))
    elseif strict && !isempty(strip(line))
      error("Unexpected DTF/FOR001 content after bus cards: $(line)")
    end
    i += 1
  end
  return DTFCase(String(path), Float64(baseMVA), params, texts, Float64.(nominal), size, branches, comps, controls, buses, outages)
end

function _branch_pu(case::DTFCase, branch::DTFBranch)
  idx = branch.voltage_level_index
  1 <= idx <= length(case.nominal_voltages_kv) || error("DTF branch $(branch.index) references missing voltage level $(idx)")
  u_ref_kv = case.nominal_voltages_kv[idx]
  z_base_ohm = u_ref_kv^2 / case.baseMVA
  return (r = branch.r_ohm / z_base_ohm, x = branch.x_ohm / z_base_ohm, g = branch.g_s * z_base_ohm, b = branch.b_s * z_base_ohm, u_ref_kv = u_ref_kv)
end

function _actual_ratio(control::Union{Nothing,DTFTransformerControl})
  control === nothing && return nothing
  control.nominal_unregulated_kv === nothing && return nothing
  control.nominal_regulated_kv === nothing && return nothing
  control.nominal_regulated_kv == 0.0 && return nothing
  # Task-1 MVP: keep the fixed actual ratio used by the existing FOR001 builder.
  return control.nominal_unregulated_kv / control.nominal_regulated_kv
end

"""
    build_net(case; bus_shunt_model = :admittance) -> Net

Build a Sparlectra `Net` from a parsed native DTF/FOR001 case. Branch
R/X/G/B are converted with the branch voltage-level index as reference voltage,
not with the from-side bus voltage. Outages remain parsed metadata and are not
executed by this Task-1 MVP.
"""
function build_net(case::DTFCase; bus_shunt_model = :admittance)::Net
  shunt_model = normalize_bus_shunt_model(bus_shunt_model)
  net = Net(name = isempty(case.texts) ? "DTF" : String(strip(case.texts[1])), baseMVA = case.baseMVA, bus_shunt_model = shunt_model)
  for bus in case.buses
    vn = bus.start_kv > 0 ? bus.start_kv : case.nominal_voltages_kv[bus.voltage_level_index]
    vm = vn / case.nominal_voltages_kv[bus.voltage_level_index]
    addBus!(net = net, busName = bus.name, vn_kV = vn, vm_pu = vm, va_deg = bus.angle_deg, isAux = false, oBusIdx = bus.index)
    if bus.pd_mw != 0.0 || bus.qd_mvar != 0.0
      addProsumer!(net = net, busName = bus.name, type = "ENERGYCONSUMER", p = bus.pd_mw, q = bus.qd_mvar, defer_bus_type_refresh = true)
    end
    if bus.pg_mw != 0.0 || bus.qg_mvar != 0.0 || bus.bus_type == 2
      type = bus.bus_type == 2 ? "SYNCHRONOUSMACHINE" : "GENERATOR"
      ref = bus.bus_type == 2 ? bus.name : nothing
      addProsumer!(net = net, busName = bus.name, type = type, p = bus.pg_mw, q = bus.qg_mvar, qMin = something(bus.qmin_mvar, -case.baseMVA), qMax = something(bus.qmax_mvar, case.baseMVA), referencePri = ref, vm_pu = vm, va_deg = bus.angle_deg)
    end
  end
  for branch in case.branches
    pu = _branch_pu(case, branch)
    from = geNetBusIdx(net = net, busName = branch.from)
    to = geNetBusIdx(net = net, busName = branch.to)
    ratedS = branch.imax_ka === nothing ? nothing : sqrt(3.0) * pu.u_ref_kv * branch.imax_ka
    if branch.kind == 'T'
      control = _find_control(branch, case.transformer_controls)
      _addPIModelTrafo_by_idx!(net = net, from = from, to = to, r_pu = pu.r, x_pu = pu.x, b_pu = pu.b, status = 1, ratedU = pu.u_ref_kv, ratedS = ratedS, ratio = _actual_ratio(control), shift_deg = something(control === nothing ? nothing : control.added_voltage_angle_deg, 0.0))
    else
      _addPIModelACLine_by_idx!(net = net, from = from, to = to, r_pu = pu.r, x_pu = pu.x, b_pu = pu.b, status = 1, ratedS = ratedS)
    end
    net.matpower_branch_metadata[length(net.branchVec)] = (
      orig_name = string(branch.kind, branch.voltage_level_index, branch.parallel_id, " ", branch.from, " -> ", branch.to),
      orig_kind = branch.kind == 'T' ? :transformer : :line,
      orig_index = branch.index,
      dtf_kind = branch.kind,
      parallel_id = branch.parallel_id,
      voltage_level_index = branch.voltage_level_index,
      u_ref_kV = pu.u_ref_kv,
      g_pu = pu.g,
      r_ohm = branch.r_ohm,
      x_ohm = branch.x_ohm,
      g_s = branch.g_s,
      b_s = branch.b_s,
    )
    if pu.g != 0.0
      addShuntMatpower!(net = net, busName = branch.from, Gs = pu.g * case.baseMVA / 2, Bs = 0.0, bus_shunt_model = shunt_model)
      addShuntMatpower!(net = net, busName = branch.to, Gs = pu.g * case.baseMVA / 2, Bs = 0.0, bus_shunt_model = shunt_model)
    end
  end
  ok, msg = validate!(net = net)
  ok || error(msg)
  return net
end

end # module DTFImporter

"""
    createNetFromDTFFile(path; baseMVA = 100.0, strict = true, bus_shunt_model = :admittance) -> Net

Read a legacy DTF/FOR001 file and build a Sparlectra `Net` without routing
through MATPOWER. Outage cards are parsed and preserved in `DTFCase` by
`DTFImporter.read_dtf`, but are not executed by this Task-1 importer MVP.
"""
function createNetFromDTFFile(path; baseMVA::Real = 100.0, strict::Bool = true, bus_shunt_model = :admittance)::Net
  return DTFImporter.build_net(DTFImporter.read_dtf(path; baseMVA = baseMVA, strict = strict); bus_shunt_model = bus_shunt_model)
end
