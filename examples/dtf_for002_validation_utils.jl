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

using Printf
using SparseArrays

_norm_name(s::AbstractString) = uppercase(replace(strip(String(s)), r"\s+" => ""))
_field(s::AbstractString, a::Int, b::Int) = a > lastindex(s) ? "" : s[a:min(b, lastindex(s))]
_numbers(s::AbstractString) = [parse(Float64, replace(m.match, 'D' => 'E', 'd' => 'E')) for m in eachmatch(r"[-+]?\d*\.?\d+(?:[EeDd][-+]?\d+)?", s)]
_normalize_for002_bus_name(s::AbstractString) = replace(strip(String(s)), r"\s+S$" => "")
_csv(x) = x === missing || x === nothing ? "" : x isa AbstractString ? "\"" * replace(x, "\"" => "\"\"") * "\"" : string(x)

mutable struct For002Scenario
  name::String
  buses::Dict{String,NamedTuple}
  branches::Vector{NamedTuple}
  total_p_loss_MW::Union{Nothing,Float64}
  total_q_loss_MVar::Union{Nothing,Float64}
  iterations::Union{Nothing,Int}
end
For002Scenario(name::String) = For002Scenario(name, Dict{String,NamedTuple}(), NamedTuple[], nothing, nothing, nothing)

function _try_parse_for002_bus(line::AbstractString)
  startswith(line, "I ") || return nothing
  name = _normalize_for002_bus_name(_field(line, 3, 13))
  isempty(name) && return nothing
  v = _numbers(_field(line, 15, 28)); gen = _numbers(_field(line, 30, 44)); load = _numbers(_field(line, 46, 60))
  length(v) == 2 && length(gen) == 2 && length(load) == 2 || return nothing
  return (bus_name = name, v_kV = v[1], va_deg = v[2], p_gen_MW = gen[1], q_gen_MVar = gen[2], p_load_MW = load[1], q_load_MVar = load[2])
end

function _try_parse_for002_branch(line::AbstractString, current_bus::Union{Nothing,AbstractString})
  current_bus === nothing && return nothing
  startswith(line, "I ") || return nothing
  strip(_field(line, 3, 13)) == "" || return nothing
  target_field = _field(line, 62, 85)
  to_bus = _normalize_for002_bus_name(_field(target_field, 6, 13))
  isempty(to_bus) && return nothing
  nr = String(strip(_field(target_field, 15, 15)))
  branch_type = String(strip(_field(target_field, 19, 22)))
  flow = _numbers(_field(line, 87, 102)); losses = _numbers(_field(line, 104, 119))
  length(flow) == 2 || return nothing
  return (from_bus = _normalize_for002_bus_name(current_bus), to_bus = to_bus, nr = nr, type = branch_type, p_MW = flow[1], q_MVar = flow[2], p_loss_MW = length(losses) >= 1 ? losses[1] : missing, q_loss_MVar = length(losses) >= 2 ? losses[2] : missing)
end

function parse_for002_ground_load_flow(path::AbstractString)::For002Scenario
  isfile(path) || throw(ArgumentError("FOR002 reference file not found: $path"))
  scenario = For002Scenario("base")
  current_bus::Union{Nothing,String} = nothing
  for raw_line in eachline(path)
    line = replace(raw_line, '\f' => ' ')
    occursin("AUSFALL DES ZWEIGES", uppercase(line)) && !isempty(scenario.buses) && break
    bus = _try_parse_for002_bus(line)
    if bus !== nothing
      scenario.buses[_norm_name(bus.bus_name)] = bus
      current_bus = String(bus.bus_name)
      continue
    end
    branch = _try_parse_for002_branch(line, current_bus)
    branch !== nothing && push!(scenario.branches, branch)
    m = match(r"VERLUSTE\s+([-+]?\d*\.?\d+)\s+MW\s+([-+]?\d*\.?\d+)\s+MVAR", uppercase(line))
    m !== nothing && (scenario.total_p_loss_MW = parse(Float64, m.captures[1]); scenario.total_q_loss_MVar = parse(Float64, m.captures[2]))
    im = match(r"NEWTONITERATION\s+(\d+)\s+ITERATIONEN", line)
    im !== nothing && (scenario.iterations = parse(Int, im.captures[1]))
  end
  return scenario
end

function write_csv(path::AbstractString, rows::Vector{<:NamedTuple}; headers = isempty(rows) ? Symbol[] : collect(keys(rows[1])))
  open(path, "w") do io
    println(io, join(String.(headers), ","))
    for row in rows
      println(io, join((_csv(getproperty(row, h)) for h in headers), ","))
    end
  end
end

function _bus_type_flags(t::Int)
  return (is_slack = t == 2, is_pv = t == 1 || t == 4, is_pq = t == 0 || t == 3)
end

function _find_ps_by_bus(net, bus_idx::Int)
  return [ps for ps in net.prosumpsVec if Sparlectra.getPosumerBusIndex(ps) == bus_idx]
end

function _bus_generation_load(net, bus_idx::Int)
  pg = qg = pl = ql = 0.0
  for ps in _find_ps_by_bus(net, bus_idx)
    if Sparlectra.isGenerator(ps)
      pg += something(ps.pRes, ps.pVal, 0.0); qg += something(ps.qRes, ps.qVal, 0.0)
    else
      pl += something(ps.pVal, 0.0); ql += something(ps.qVal, 0.0)
    end
  end
  return (pg = pg, qg = qg, pl = pl, ql = ql)
end

function branch_ref_dict(s::For002Scenario)
  d = Dict{Tuple{String,String,String},NamedTuple}()
  for b in s.branches
    d[(_norm_name(b.from_bus), _norm_name(b.to_bus), uppercase(strip(b.nr)))] = b
  end
  return d
end

function for002_state_residual_rows(net, case, ref::For002Scenario)
  Y = Sparlectra.createYBUS(net = net, sparse = false)
  V = zeros(ComplexF64, length(net.nodeVec))
  rows = NamedTuple[]
  for bus in case.buses
    idx = bus.index
    node = net.nodeVec[idx]
    key = _norm_name(bus.name)
    haskey(ref.buses, key) || continue
    fb = ref.buses[key]
    vm = fb.v_kV / node.comp.cVN
    va = fb.va_deg
    V[idx] = vm * cis(deg2rad(va))
  end
  S = V .* conj.(Y * V) .* net.baseMVA
  for bus in case.buses
    idx = bus.index
    node = net.nodeVec[idx]
    key = _norm_name(bus.name)
    haskey(ref.buses, key) || continue
    fb = ref.buses[key]
    pnet = fb.p_gen_MW - fb.p_load_MW
    qnet = fb.q_gen_MVar - fb.q_load_MVar
    push!(rows, (bus_name = bus.name, forced_for002_vm_pu = fb.v_kV / node.comp.cVN, forced_for002_va_deg = fb.va_deg, ybus_calc_p_MW = real(S[idx]), ybus_calc_q_MVar = imag(S[idx]), for002_bus_net_p_MW = pnet, for002_bus_net_q_MVar = qnet, d_p_MW = real(S[idx]) - pnet, d_q_MVar = imag(S[idx]) - qnet))
  end
  return rows
end
