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

# Developer notes:
# - Shared helpers validate native Testnetz13 DTF/FOR002 examples.
# - The intended path is DTFImporter.read_dtf -> DTFImporter.build_net -> runpf!;
#   MATPOWER import/export and the generated FOR001 builder are intentionally bypassed.
# - FOR002 is a legacy textual reference report, not a machine-stable schema.
# - State residuals are diagnostic only; branch-flow and solved generator/slack
#   comparisons are the stronger validation signals for now.
# - Console output should stay concise; detailed rows belong in CSV/Markdown and
#   in optional detailed-mode return values.

# FOR002 and DTF names differ in spacing/case, so matching keys must ignore
# report padding while preserving the human-readable names in output rows.
_norm_name(s::AbstractString) = uppercase(replace(strip(String(s)), r"\s+" => ""))
_field(s::AbstractString, a::Int, b::Int) = a > lastindex(s) ? "" : s[a:min(b, lastindex(s))]
_numbers(s::AbstractString) = [parse(Float64, replace(m.match, 'D' => 'E', 'd' => 'E')) for m in eachmatch(r"[-+]?\d*\.?\d+(?:[EeDd][-+]?\d+)?", s)]
_normalize_for002_bus_name(s::AbstractString) = replace(strip(replace(String(s), "*" => "")), r"\s+S$" => "")
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
  # The FOR002 report is fixed-width-ish, but alignment varies around names and
  # numeric fields; parse by tolerant slices plus numeric extraction.
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
  # Branch-flow rows are directional: the current bus is the reported source,
  # and the later comparison matches both source->target and target->source rows.
  return (from_bus = _normalize_for002_bus_name(current_bus), to_bus = to_bus, nr = nr, type = branch_type, p_MW = flow[1], q_MVar = flow[2], p_loss_MW = length(losses) >= 1 ? losses[1] : missing, q_loss_MVar = length(losses) >= 2 ? losses[2] : missing)
end

function parse_for002_ground_load_flow(path::AbstractString)::For002Scenario
  isfile(path) || throw(ArgumentError("FOR002 reference file not found: $path"))
  scenario = For002Scenario("base")
  current_bus::Union{Nothing,String} = nothing
  for raw_line in eachline(path)
    line = replace(raw_line, '\f' => ' ')
    # Stop at the first outage heading; the preceding block is the base case.
    occursin("AUSFALL DES ZWEIGES", uppercase(line)) && !isempty(scenario.buses) && break
    bus = _try_parse_for002_bus(line)
    if bus !== nothing
      scenario.buses[_norm_name(bus.bus_name)] = bus
      current_bus = String(bus.bus_name)
      continue
    end
    branch = _try_parse_for002_branch(line, current_bus)
    branch !== nothing && push!(scenario.branches, branch)
    m = match(r"VERLUSTE.*?P\s*=\s*([-+]?\d*\.?\d+)\s+MW\s+Q\s*=\s*([-+]?\d*\.?\d+)\s+MVAR", uppercase(line))
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
  pg_spec = qg_spec = pg_res = qg_res = pl = ql = 0.0
  has_generator = has_load = is_regulating = false
  for ps in _find_ps_by_bus(net, bus_idx)
    if Sparlectra.isGenerator(ps)
      has_generator = true
      is_regulating |= ps.isRegulated
      p_spec = something(ps.pVal, 0.0); q_spec = something(ps.qVal, 0.0)
      pg_spec += p_spec; qg_spec += q_spec
      pg_res += something(ps.pRes, p_spec); qg_res += something(ps.qRes, q_spec)
    else
      has_load = true
      pl += something(ps.pVal, 0.0); ql += something(ps.qVal, 0.0)
    end
  end
  return (pg_spec = pg_spec, qg_spec = qg_spec, pg_res = pg_res, qg_res = qg_res, pl = pl, ql = ql, has_generator = has_generator, has_load = has_load, is_regulating = is_regulating)
end

function _branch_kcl_arrays(net)
  branch_p = zeros(Float64, length(net.nodeVec)); branch_q = zeros(Float64, length(net.nodeVec))
  # Use solved endpoint flows. They already include branch charging and taps as
  # represented by the branch model; explicit bus shunts are not added here.
  for br in net.branchVec
    if br.fBranchFlow !== nothing
      branch_p[br.fromBus] += something(br.fBranchFlow.pFlow, 0.0)
      branch_q[br.fromBus] += something(br.fBranchFlow.qFlow, 0.0)
    end
    if br.tBranchFlow !== nothing
      branch_p[br.toBus] += something(br.tBranchFlow.pFlow, 0.0)
      branch_q[br.toBus] += something(br.tBranchFlow.qFlow, 0.0)
    end
  end
  return branch_p, branch_q
end

function branch_kcl_rows(net, case, ref::For002Scenario)
  branch_p, branch_q = _branch_kcl_arrays(net)
  rows = NamedTuple[]
  for bus in case.buses
    idx = bus.index; node = net.nodeVec[idx]; key = _norm_name(bus.name)
    haskey(ref.buses, key) || continue
    fb = ref.buses[key]; flags = _bus_type_flags(bus.bus_type); gl = _bus_generation_load(net, idx)
    for002_pnet = fb.p_gen_MW - fb.p_load_MW; for002_qnet = fb.q_gen_MVar - fb.q_load_MVar
    spec_pnet = gl.pg_spec - gl.pl; spec_qnet = gl.qg_spec - gl.ql
    solved_pg = (flags.is_slack || flags.is_pv) ? branch_p[idx] + gl.pl : gl.pg_res
    solved_qg = (flags.is_slack || flags.is_pv) ? branch_q[idx] + gl.ql : gl.qg_res
    res_pnet = solved_pg - gl.pl; res_qnet = solved_qg - gl.ql
    notes = "branch_kcl sums solved branch endpoint flows; branch charging/taps are included in endpoint flows; explicit bus shunts are not added separately"
    (abs(node._pShunt) > 1e-12 || abs(node._qShunt) > 1e-12) && (notes *= "; bus has nonzero explicit shunt fields")
    push!(rows, (bus_name = bus.name, dtf_bus_type = bus.bus_type, is_slack = flags.is_slack, is_pv = flags.is_pv, is_pq = flags.is_pq,
      for002_p_gen_MW = fb.p_gen_MW, for002_q_gen_MVar = fb.q_gen_MVar, for002_p_load_MW = fb.p_load_MW, for002_q_load_MVar = fb.q_load_MVar,
      for002_p_net_MW = for002_pnet, for002_q_net_MVar = for002_qnet,
      model_specified_p_net_MW = spec_pnet, model_specified_q_net_MVar = spec_qnet,
      model_result_p_net_MW = res_pnet, model_result_q_net_MVar = res_qnet,
      branch_kcl_p_net_MW = branch_p[idx], branch_kcl_q_net_MVar = branch_q[idx],
      d_model_result_vs_for002_p_MW = res_pnet - for002_pnet, d_model_result_vs_for002_q_MVar = res_qnet - for002_qnet,
      d_branch_kcl_vs_for002_p_MW = branch_p[idx] - for002_pnet, d_branch_kcl_vs_for002_q_MVar = branch_q[idx] - for002_qnet,
      d_branch_kcl_vs_model_result_p_MW = branch_p[idx] - res_pnet, d_branch_kcl_vs_model_result_q_MVar = branch_q[idx] - res_qnet,
      notes = notes))
  end
  return rows
end

function _likely_q_explanation(genrow, kclrow, residualrow)
  genrow.is_slack && return "slack_q_reporting_difference"
  genrow.is_regulating && abs(genrow.d_qg_result_vs_for002_MVar) > 1e-6 && return "pv_q_result_vs_input_difference"
  genrow.has_generator && !genrow.is_regulating && abs(genrow.model_qg_specified_MVar - genrow.for002_qg_MVar) <= 1e-6 && return "pq_generator_fixed_q_matches_input"
  abs(kclrow.d_branch_kcl_vs_for002_q_MVar) > 10 && abs(kclrow.d_branch_kcl_vs_model_result_q_MVar) <= 1e-6 && return "for002_bus_q_net_not_equal_simple_gen_minus_load"
  abs(kclrow.d_branch_kcl_vs_for002_q_MVar) > 10 && return "branch_flows_match_but_bus_q_table_differs"
  return "needs_manual_review"
end

function q_semantics_diagnostic_rows(gen_rows, kcl_rows, residual_rows)
  gen_by_bus = Dict(_norm_name(r.bus_name) => r for r in gen_rows)
  kcl_by_bus = Dict(_norm_name(r.bus_name) => r for r in kcl_rows)
  res_by_bus = Dict(_norm_name(r.bus_name) => r for r in residual_rows)
  wanted = Set{Tuple{String,String}}()
  for r in first(sort(gen_rows, by = r -> -abs(r.d_qg_result_vs_for002_MVar)), min(10, length(gen_rows)))
    push!(wanted, ("top_abs_d_qg_result_vs_for002", _norm_name(r.bus_name)))
  end
  for r in first(sort(residual_rows, by = r -> -abs(r.d_q_MVar)), min(10, length(residual_rows)))
    push!(wanted, ("top_abs_state_residual_q", _norm_name(r.bus_name)))
  end
  for r in first(sort(kcl_rows, by = r -> -abs(r.d_branch_kcl_vs_for002_q_MVar)), min(10, length(kcl_rows)))
    push!(wanted, ("top_abs_branch_kcl_vs_for002_q", _norm_name(r.bus_name)))
  end
  for r in gen_rows
    (r.is_slack || r.has_generator) && push!(wanted, (r.is_slack ? "slack_bus" : "generator_bus", _norm_name(r.bus_name)))
  end
  for name in ["BETA1 S1", "BETA2 S1", "DELTA1S1", "DELTA2S1"]
    push!(wanted, ("transformer_adjacent_bus", _norm_name(name)))
  end
  rows = NamedTuple[]
  for (category, key) in sort(collect(wanted))
    haskey(gen_by_bus, key) && haskey(kcl_by_bus, key) && haskey(res_by_bus, key) || continue
    g = gen_by_bus[key]; k = kcl_by_bus[key]; r = res_by_bus[key]
    push!(rows, (category = category, bus_name = g.bus_name, dtf_bus_type = g.dtf_bus_type, is_slack = g.is_slack, is_regulating = g.is_regulating,
      for002_q_gen_MVar = g.for002_qg_MVar, model_qg_specified_MVar = g.model_qg_specified_MVar, model_qg_result_MVar = g.model_qg_result_MVar,
      d_qg_result_vs_for002_MVar = g.d_qg_result_vs_for002_MVar, for002_q_net_MVar = g.for002_q_net_MVar,
      model_q_net_result_MVar = g.model_q_net_result_MVar, branch_kcl_q_net_MVar = k.branch_kcl_q_net_MVar,
      state_residual_q_MVar = r.d_q_MVar, likely_explanation = _likely_q_explanation(g, k, r)))
  end
  return rows
end

function branch_ref_dict(s::For002Scenario)
  d = Dict{Tuple{String,String,String},NamedTuple}()
  for b in s.branches
    # Keep the direction in the key; FOR002 prints endpoint flows separately for
    # each direction and the comparison must not collapse them.
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
    # Force the printed FOR002 voltage state into the native Ybus, calculate the
    # implied injection, then compare it with the printed FOR002 bus table.
    # These residuals can be larger than branch-flow deviations because printed
    # magnitudes/angles are rounded and transformer-adjacent nodes are sensitive.
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

mutable struct For002OutageScenario
  scenario_index::Int
  raw_heading::String
  outage_kind::Union{Nothing,Char}
  parallel_id::Union{Nothing,String}
  from_bus::Union{Nothing,String}
  to_bus::Union{Nothing,String}
  buses::Dict{String,NamedTuple}
  branches::Vector{NamedTuple}
  total_p_loss_MW::Union{Nothing,Float64}
  total_q_loss_MVar::Union{Nothing,Float64}
  iterations::Union{Nothing,Int}
end

function For002OutageScenario(i::Int, heading::String)
  kind, pid, from_bus, to_bus = _parse_for002_outage_heading(heading)
  return For002OutageScenario(i, heading, kind, pid, from_bus, to_bus, Dict{String,NamedTuple}(), NamedTuple[], nothing, nothing, nothing)
end

function _parse_for002_outage_heading(line::AbstractString)
  clean = replace(strip(String(line)), r"\s+" => " ")
  # Headings encode parallel id and endpoints; the DTF branch kind is usually
  # unavailable in the legacy FOR002 text and is matched from FOR001 instead.
  m = match(r"AUSFALL DES ZWEIGES\s*([A-Z0-9]?)\s*VON\s+(.+?)\s+NACH\s+(.+?)(?:\s+I)?$", uppercase(clean))
  m === nothing && return (nothing, nothing, nothing, nothing)
  pid = strip(m.captures[1])
  from_bus = _normalize_for002_bus_name(m.captures[2])
  to_bus = _normalize_for002_bus_name(replace(m.captures[3], r"\s+I$" => ""))
  return (nothing, isempty(pid) ? "" : pid, from_bus, to_bus)
end

function _scenario_from_outage(s::For002OutageScenario)
  sc = For002Scenario("outage_$(s.scenario_index)")
  sc.buses = s.buses; sc.branches = s.branches
  sc.total_p_loss_MW = s.total_p_loss_MW; sc.total_q_loss_MVar = s.total_q_loss_MVar; sc.iterations = s.iterations
  return sc
end

function parse_for002_outage_scenarios(path::AbstractString)::Vector{For002OutageScenario}
  isfile(path) || throw(ArgumentError("FOR002 reference file not found: $path"))
  scenarios = For002OutageScenario[]
  current::Union{Nothing,For002OutageScenario} = nothing
  current_bus::Union{Nothing,String} = nothing
  for raw_line in eachline(path)
    line = replace(raw_line, '\f' => ' ')
    if occursin("AUSFALL DES ZWEIGES", uppercase(line))
      # Every outage heading starts a new legacy report block; only blocks that
      # match DTF-listed outage cards are executed by the outage example.
      current = For002OutageScenario(length(scenarios) + 1, String(strip(line)))
      push!(scenarios, current)
      current_bus = nothing
      continue
    end
    current === nothing && continue
    bus = _try_parse_for002_bus(line)
    if bus !== nothing
      current.buses[_norm_name(bus.bus_name)] = bus
      current_bus = String(bus.bus_name)
      continue
    end
    branch = _try_parse_for002_branch(line, current_bus)
    branch !== nothing && push!(current.branches, branch)
    m = match(r"VERLUSTE.*?P\s*=\s*([-+]?\d*\.?\d+)\s+MW\s+Q\s*=\s*([-+]?\d*\.?\d+)\s+MVAR", uppercase(line))
    m !== nothing && (current.total_p_loss_MW = parse(Float64, m.captures[1]); current.total_q_loss_MVar = parse(Float64, m.captures[2]))
    im = match(r"NEWTONITERATION\s+(\d+)\s+ITERATIONEN", line)
    im !== nothing && (current.iterations = parse(Int, im.captures[1]))
  end
  return scenarios
end
