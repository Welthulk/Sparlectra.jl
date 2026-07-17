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

using ..Sparlectra: Net, addBus!, addProsumer!, _addPIModelACLine_by_idx!, _addPIModelTrafo_by_idx!, geNetBusIdx, validate!, normalize_bus_shunt_model, calcSkewAngleTap, calcTapCorrectedRX, transformer_config

export DTFCase, DTFParams, DTFSize, DTFBranch, DTFBus, DTFCompensation, DTFTransformerControl, DTFOutage, DTFTrailingRecord, read_dtf, build_net,
  dtf_branch_key, find_outage_branch_indices, outage_match_diagnostic, apply_single_branch_outage!, case_summary, outage_label

const DTF_TRANSFORMER_RATIO_MODES = (:neutral_one, :winding_over_network)

function _normalize_transformer_ratio_mode(mode)::Symbol
  sym = mode isa Symbol ? mode : Symbol(mode)
  sym in DTF_TRANSFORMER_RATIO_MODES || throw(ArgumentError("transformer_ratio_mode must be one of $(collect(DTF_TRANSFORMER_RATIO_MODES)); got $(sym)."))
  return sym
end

"""Native DTF importer parameter card values preserved for auditing."""
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
  quadrature_max_steps::Union{Nothing,Int}
  quadrature_actual_step::Union{Nothing,Int}
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

struct DTFTrailingRecord
  raw::String
  index::Int
  kind::Symbol
  interpreted_kind::Symbol
  paired_marker_index::Union{Nothing,Int}
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
  trailing_records::Vector{DTFTrailingRecord}
end

_slice(s::AbstractString, a::Int, b::Int) = a > lastindex(s) ? "" : s[a:min(b, lastindex(s))]
_norm_name(s::AbstractString) = uppercase(replace(strip(String(s)), r"\s+" => ""))
_parse_float(s::AbstractString) = isempty(strip(s)) ? nothing : parse(Float64, replace(strip(s), 'D' => 'E', 'd' => 'E'))
_parse_int(s::AbstractString) = (x = _parse_float(s); x === nothing ? nothing : Int(round(x)))
_numbers(s::AbstractString) = [parse(Float64, replace(m.match, 'D' => 'E', 'd' => 'E')) for m in eachmatch(r"[-+]?\d*\.?\d+(?:[EeDd][-+]?\d+)?", s)]
function _nominal_voltages(s::AbstractString)
  vals = _numbers(s)
  length(vals) > 1 && occursin(r"\s", strip(s)) && return vals
  chunks = Float64[]
  for i in 1:5:lastindex(s)
    field = strip(_slice(s, i, i + 4))
    isempty(field) && continue
    push!(chunks, parse(Float64, field))
  end
  return chunks
end

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

"""Return the strict branch-matching key used for DTF parallel branches."""
dtf_branch_key(x) = (kind = x.kind, voltage_level_index = x.voltage_level_index,
  parallel_id = uppercase(strip(x.parallel_id)), from = _norm_name(x.from), to = _norm_name(x.to))

"""Return a concise human-readable DTF outage label."""
outage_label(o::DTFOutage) = string(o.kind, o.voltage_level_index, o.parallel_id, " ", o.from, " -> ", o.to)

"""Find native branch indices matching a DTF outage card with strict parallel-branch keys."""
function find_outage_branch_indices(case::DTFCase, outage::DTFOutage)::Vector{Int}
  key = dtf_branch_key(outage)
  return [i for (i, br) in enumerate(case.branches) if dtf_branch_key(br) == key]
end

"""Build a clear diagnostic for missing or ambiguous DTF outage branch matches."""
function outage_match_diagnostic(case::DTFCase, outage::DTFOutage, matches::Vector{Int}=find_outage_branch_indices(case, outage))
  if isempty(matches)
    return "No native DTF branch matches outage $(outage_label(outage)); key=$(dtf_branch_key(outage))"
  end
  if length(matches) > 1
    candidates = join(matches, ", ")
    return "Ambiguous native DTF branch match for outage $(outage_label(outage)); candidates=$(candidates)"
  end
  return "Matched native DTF branch $(only(matches)) for outage $(outage_label(outage))"
end

"""Safely set exactly one branch out of service and verify no other status changed."""
function apply_single_branch_outage!(net::Net, branch_index::Int)
  1 <= branch_index <= length(net.branchVec) || throw(ArgumentError("branch index $branch_index is outside 1:$(length(net.branchVec))"))
  before = [br.status for br in net.branchVec]
  before[branch_index] == 1 || throw(ArgumentError("matched branch $branch_index was not initially in service"))
  net.branchVec[branch_index].status = 0
  after = [br.status for br in net.branchVec]
  changed = findall(i -> before[i] != after[i], eachindex(before))
  length(changed) == 1 && only(changed) == branch_index || throw(ArgumentError("branch outage changed $(length(changed)) branches, expected exactly one"))
  net.branchVec[branch_index].status == 0 || throw(ArgumentError("matched branch $branch_index is still in service after outage"))
  return nothing
end

"""Return stable summary counts for a parsed DTF case."""
case_summary(case::DTFCase) = (baseMVA = case.baseMVA, bus_count = length(case.buses), branch_count = length(case.branches),
  line_count = count(b -> b.kind != 'T', case.branches), transformer_count = count(b -> b.kind == 'T', case.branches),
  outage_count = length(case.outages), trailing_record_count = length(case.trailing_records), slack_bus = case.size.slack)

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
  control_type = String(strip(_slice(line, 1, 1)))
  variable_flag = String(strip(_slice(line, 2, 2)))
  parallel_id = String(strip(_slice(line, 3, 3)))
  phase_shifter_flag = String(strip(_slice(line, 4, 4)))
  from = String(strip(_slice(line, 6, 13)))
  to = String(strip(_slice(line, 16, 23)))

  # DTF transformer-control cards are fixed-column records. Adjacent numeric
  # fields may appear without separators (for example `400.0231.012.50`), so
  # this parser must not use the generic `_numbers` scanner here.
  return DTFTransformerControl(line, index, control_type, variable_flag, parallel_id, phase_shifter_flag, from, to,
    _parse_float(_slice(line, 26, 30)),
    _parse_float(_slice(line, 31, 35)),
    _parse_float(_slice(line, 36, 40)),
    _parse_float(_slice(line, 41, 45)),
    _parse_int(_slice(line, 46, 50)),
    _parse_int(_slice(line, 51, 55)),
    _parse_float(_slice(line, 56, 60)),
    _parse_int(_slice(line, 61, 65)),
    _parse_int(_slice(line, 66, 70)))
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

Read a legacy DTF fixed-column input file into a typed `DTFCase`.
This MVP parses outage cards into `DTFOutage` records but intentionally does
not execute outage simulations.
"""
function read_dtf(path; baseMVA::Real = 100.0, strict::Bool = true)::DTFCase
  lines = readlines(path)
  length(lines) >= 7 || error("DTF file is too short: $(path)")
  params = DTFParams(lines[1], _numbers(lines[1]))
  texts = String.(lines[2:5])
  nominal = _nominal_voltages(lines[6])
  size_vals = _numbers(lines[7])
  length(size_vals) >= 4 || error("DTF network size card is incomplete")
  size_tokens = split(strip(lines[7]))
  size = DTFSize(lines[7], Int(size_vals[1]), Int(size_vals[2]), Int(size_vals[3]), Int(size_vals[4]), length(size_tokens) >= 5 ? String(size_tokens[end]) : "")
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
  trailing = DTFTrailingRecord[]
  in_outage = false
  outage_index = 0
  trailing_index = 0
  pending_variant_marker = nothing
  while i <= length(lines)
    line = lines[i]
    token = uppercase(strip(line))
    if token == "AUSFALL"
      in_outage = true
      pending_variant_marker = nothing
    elseif token == "ENDE"
      in_outage = false
      pending_variant_marker = nothing
    elseif in_outage && !isempty(strip(line))
      outage_index += 1
      push!(outages, _parse_outage(line, outage_index))
    elseif token == "A"
      trailing_index += 1
      push!(trailing, DTFTrailingRecord(line, trailing_index, :variant_a_marker, :variant_branch_echo_marker, nothing))
      pending_variant_marker = trailing_index
    elseif pending_variant_marker !== nothing && !isempty(strip(line)) && first(strip(line)) in ('L', 'T')
      # legacy validation cases B-E contain post-bus standalone `A` cards followed by a
      # branch-like line for every in-service branch.  The payload duplicates
      # the already counted branch data in a tighter legacy print format; it is
      # preserved as a trailing branch echo for diagnostics, not added as an
      # electrical branch or outage.
      trailing_index += 1
      push!(trailing, DTFTrailingRecord(line, trailing_index, :variant_a_payload, :variant_branch_echo, pending_variant_marker))
      pending_variant_marker = nothing
    elseif strict && !isempty(strip(line))
      error("Unexpected DTF content after bus cards: $(line)")
    end
    i += 1
  end
  return DTFCase(String(path), Float64(baseMVA), params, texts, Float64.(nominal), size, branches, comps, controls, buses, outages, trailing)
end

function _dtf_nominal_voltage_levels(case::DTFCase; legacy_voltage_level_collapse_230kv::Bool = false)
  legacy_voltage_level_collapse_230kv || return case.nominal_voltages_kv
  levels = copy(case.nominal_voltages_kv)
  for i in eachindex(levels)
    (isapprox(levels[i], 231.0; atol = 1e-9) || isapprox(levels[i], 230.0; atol = 1e-9)) && (levels[i] = 230.0)
  end
  return levels
end

function _branch_pu(case::DTFCase, branch::DTFBranch; nominal_voltages_kv = case.nominal_voltages_kv)
  idx = branch.voltage_level_index
  1 <= idx <= length(nominal_voltages_kv) || error("DTF branch $(branch.index) references missing voltage level $(idx)")
  u_ref_kv = nominal_voltages_kv[idx]
  z_base_ohm = u_ref_kv^2 / case.baseMVA
  return (r = branch.r_ohm / z_base_ohm, x = branch.x_ohm / z_base_ohm, g = branch.g_s * z_base_ohm, b = branch.b_s * z_base_ohm, u_ref_kv = u_ref_kv)
end

function _dtf_effective_transformer_tap(case::DTFCase, branch::DTFBranch, control::Union{Nothing,DTFTransformerControl}, from_bus::DTFBus, to_bus::DTFBus; nominal_voltages_kv = case.nominal_voltages_kv, transformer_ratio_mode = :neutral_one)
  ratio_mode = _normalize_transformer_ratio_mode(transformer_ratio_mode)
  no_control = (ratio = nothing, shift_deg = 0.0, model = :fixed_ratio, tap_fraction = 0.0,
    skew_angle_deg = 0.0, effective_complex = 1.0 + 0.0im, convention = :dtf_regulating_vector_reciprocal,
    nominal_unregulated_kv = nothing, nominal_regulated_kv = nothing, from_bus_vn_kV = nominal_voltages_kv[from_bus.voltage_level_index],
    to_bus_vn_kV = nominal_voltages_kv[to_bus.voltage_level_index], winding_over_network_base_ratio = 1.0,
    transformer_ratio_mode = ratio_mode, base_ratio_used = 1.0, effective_ratio = 1.0)
  control === nothing && return no_control
  control.nominal_unregulated_kv === nothing && return nothing
  control.nominal_regulated_kv === nothing && return nothing
  control.nominal_regulated_kv == 0.0 && return nothing
  from_vn = nominal_voltages_kv[from_bus.voltage_level_index]
  to_vn = nominal_voltages_kv[to_bus.voltage_level_index]
  to_vn == 0.0 && return nothing

  winding_over_network_base_ratio = (control.nominal_unregulated_kv / control.nominal_regulated_kv) / (from_vn / to_vn)
  base_ratio = ratio_mode == :neutral_one ? 1.0 : winding_over_network_base_ratio
  tap_fraction = 0.0
  skew_angle_deg = something(control.added_voltage_angle_deg, 0.0)
  model = :fixed_ratio
  effective_complex = 1.0 + 0.0im
  relative_ratio = 1.0
  shift_deg = 0.0
  if control.longitudinal_range_percent !== nothing && control.max_tap_step !== nothing &&
      control.actual_tap_step !== nothing && control.max_tap_step != 0
    tap_fraction = (control.longitudinal_range_percent / 100.0) * control.actual_tap_step / control.max_tap_step
    if tap_fraction != 0.0
      model = skew_angle_deg == 0.0 ? :longitudinal : :skew_angle
      tap_result = calcSkewAngleTap(; tap_fraction = tap_fraction, skew_angle_deg = skew_angle_deg)
      effective_complex = tap_result.regulating_vector
      relative_ratio = tap_result.effective_ratio
      shift_deg = model == :skew_angle ? tap_result.effective_shift_deg : 0.0
    end
  end
  # DTF compatibility treats winding voltages as nameplate/control
  # metadata by default.  The off-nominal network tap is neutral at Stufe 0 and
  # only receives longitudinal/skew deviations relative to that neutral point.
  ratio = base_ratio * relative_ratio
  return (ratio = ratio, shift_deg = shift_deg, model = model, tap_fraction = tap_fraction,
    skew_angle_deg = skew_angle_deg, effective_complex = effective_complex,
    convention = :dtf_regulating_vector_reciprocal,
    nominal_unregulated_kv = control.nominal_unregulated_kv,
    nominal_regulated_kv = control.nominal_regulated_kv,
    from_bus_vn_kV = from_vn,
    to_bus_vn_kV = to_vn,
    winding_over_network_base_ratio = winding_over_network_base_ratio,
    transformer_ratio_mode = ratio_mode,
    base_ratio_used = base_ratio,
    effective_ratio = relative_ratio)
end

function _dtf_transformer_ratio(case::DTFCase, branch::DTFBranch, control::Union{Nothing,DTFTransformerControl}, from_bus::DTFBus, to_bus::DTFBus)
  return _dtf_effective_transformer_tap(case, branch, control, from_bus, to_bus).ratio
end

"""
    build_net(case; bus_shunt_model = :admittance) -> Net

Build a Sparlectra `Net` from a parsed native DTF case. Branch
R/X/G/B are converted with the branch voltage-level index as reference voltage,
not with the from-side bus voltage. Outages remain parsed metadata and are not
executed by this Task-1 MVP.

`tap_changer_model` selects the tap-changer model applied to all transformers
(`transformer.tap_changer_model` in the central configuration; `nothing` reads
the active configuration). `:ideal` keeps the neutral-position series
impedance; `:impedance_correction` re-refers R/X through the tapped winding via
the central `calcTapCorrectedRX` using the parsed regulating vector
`1 + f·e^(jφ)`.
"""
function build_net(case::DTFCase; bus_shunt_model = :admittance, legacy_voltage_level_collapse_230kv::Bool = false, transformer_ratio_mode = :neutral_one, tap_changer_model::Union{Nothing,Symbol} = nothing)::Net
  ratio_mode = _normalize_transformer_ratio_mode(transformer_ratio_mode)
  tap_model = tap_changer_model === nothing ? transformer_config().tap_changer_model : tap_changer_model
  shunt_model = normalize_bus_shunt_model(bus_shunt_model)
  nominal_voltages_kv = _dtf_nominal_voltage_levels(case; legacy_voltage_level_collapse_230kv = legacy_voltage_level_collapse_230kv)
  net = Net(name = isempty(case.source_path) ? "DTF" : basename(case.source_path), baseMVA = case.baseMVA, bus_shunt_model = shunt_model)
  bus_by_name = Dict(bus.name => bus for bus in case.buses)
  for bus in case.buses
    vn = nominal_voltages_kv[bus.voltage_level_index]
    vm = bus.start_kv > 0 ? bus.start_kv / vn : 1.0
    addBus!(net = net, busName = bus.name, vn_kV = vn, vm_pu = vm, va_deg = bus.angle_deg, isAux = false, oBusIdx = bus.index)
    net.busOriginalNameDict[geNetBusIdx(net = net, busName = bus.name)] = bus.name
    if bus.pd_mw != 0.0 || bus.qd_mvar != 0.0
      addProsumer!(net = net, busName = bus.name, type = "ENERGYCONSUMER", p = bus.pd_mw, q = bus.qd_mvar, defer_bus_type_refresh = true)
    end
    if bus.pg_mw != 0.0 || bus.qg_mvar != 0.0 || bus.bus_type == 2
      if bus.bus_type == 2
        addProsumer!(net = net, busName = bus.name, type = "SYNCHRONOUSMACHINE", p = bus.pg_mw, q = bus.qg_mvar, qMin = bus.qmin_mvar, qMax = bus.qmax_mvar, referencePri = bus.name, vm_pu = vm, va_deg = bus.angle_deg)
      elseif bus.bus_type == 1 || bus.bus_type == 4
        addProsumer!(net = net, busName = bus.name, type = "GENERATOR", p = bus.pg_mw, q = bus.qg_mvar, qMin = bus.qmin_mvar, qMax = bus.qmax_mvar, vm_pu = vm, va_deg = bus.angle_deg)
      else
        # DTF bus types 0 and 3 are PQ buses. Do not pass `vm_pu` here:
        # `addProsumer!` treats the presence of a voltage setpoint as voltage
        # regulation, which would incorrectly turn PQ generator injections into
        # PV buses solely because Pg/Qg are present.
        addProsumer!(net = net, busName = bus.name, type = "GENERATOR", p = bus.pg_mw, q = bus.qg_mvar, qMin = bus.qmin_mvar, qMax = bus.qmax_mvar, isRegulated = false)
      end
    end
  end
  for branch in case.branches
    pu = _branch_pu(case, branch; nominal_voltages_kv = nominal_voltages_kv)
    from = geNetBusIdx(net = net, busName = branch.from)
    to = geNetBusIdx(net = net, busName = branch.to)
    # TODO: confirm whether DTF branch current ratings should use the branch
    # voltage-level reference or a physical terminal voltage. Preserve Task-1
    # behavior for now and expose the chosen/reference voltages in metadata.
    ratedS = branch.imax_ka === nothing ? nothing : sqrt(3.0) * pu.u_ref_kv * branch.imax_ka
    control = nothing
    ratio = 1.0
    if branch.kind == 'T'
      control = _find_control(branch, case.transformer_controls)
      tap = _dtf_effective_transformer_tap(case, branch, control, bus_by_name[branch.from], bus_by_name[branch.to]; nominal_voltages_kv = nominal_voltages_kv, transformer_ratio_mode = ratio_mode)
      ratio = tap.ratio
      shift_deg = tap.shift_deg
      # Central tap-changer impedance model (src/equicircuit.jl): with
      # :impedance_correction the series impedance is re-referred through the
      # tapped winding using the parsed regulating vector 1 + f·e^(jφ).
      tap_rx = calcTapCorrectedRX(r_pu = pu.r, x_pu = pu.x, tap_changer_model = tap_model, tap_fraction = tap.tap_fraction, skew_angle_deg = tap.skew_angle_deg)
      _addPIModelTrafo_by_idx!(net = net, from = from, to = to, r_pu = tap_rx.r_pu, x_pu = tap_rx.x_pu, b_pu = pu.b, g_pu = pu.g, status = 1, ratedU = pu.u_ref_kv, ratedS = ratedS, ratio = ratio, shift_deg = shift_deg)
    else
      _addPIModelACLine_by_idx!(net = net, from = from, to = to, r_pu = pu.r, x_pu = pu.x, b_pu = pu.b, status = 1, ratedS = ratedS)
    end
    from_bus = bus_by_name[branch.from]
    to_bus = bus_by_name[branch.to]
    net.matpower_branch_metadata[length(net.branchVec)] = (
      orig_name = string(branch.kind, branch.voltage_level_index, branch.parallel_id, " ", branch.from, " -> ", branch.to),
      source_label = string(branch.kind, branch.voltage_level_index, branch.parallel_id),
      orig_kind = branch.kind == 'T' ? :transformer : :line,
      orig_index = branch.index,
      dtf_kind = branch.kind,
      parallel_id = branch.parallel_id,
      voltage_level_index = branch.voltage_level_index,
      u_ref_kV = pu.u_ref_kv,
      rate_u_ref_kV = pu.u_ref_kv,
      from_bus_vn_kV = nominal_voltages_kv[from_bus.voltage_level_index],
      to_bus_vn_kV = nominal_voltages_kv[to_bus.voltage_level_index],
      g_pu = pu.g,
      b_pu = pu.b,
      transformer_loss_allocation = branch.kind == 'T' ? :native_branch_pi : :not_transformer,
      active_no_load_g_pu = branch.kind == 'T' ? pu.g : 0.0,
      r_ohm = branch.r_ohm,
      x_ohm = branch.x_ohm,
      g_s = branch.g_s,
      b_s = branch.b_s,
      tap_ratio = branch.kind == 'T' ? ratio : 1.0,
      nominal_unregulated_kv = branch.kind == 'T' ? tap.nominal_unregulated_kv : nothing,
      nominal_regulated_kv = branch.kind == 'T' ? tap.nominal_regulated_kv : nothing,
      winding_over_network_base_ratio = branch.kind == 'T' ? tap.winding_over_network_base_ratio : 1.0,
      transformer_ratio_mode = branch.kind == 'T' ? tap.transformer_ratio_mode : ratio_mode,
      base_ratio_used = branch.kind == 'T' ? tap.base_ratio_used : 1.0,
      actual_tap_step = control === nothing ? nothing : control.actual_tap_step,
      max_tap_step = control === nothing ? nothing : control.max_tap_step,
      phase_shift_deg = branch.kind == 'T' ? tap.shift_deg : 0.0,
      dtf_control_model = branch.kind == 'T' ? tap.model : :not_transformer,
      tap_fraction = branch.kind == 'T' ? tap.tap_fraction : 0.0,
      skew_angle_deg = branch.kind == 'T' ? tap.skew_angle_deg : 0.0,
      effective_ratio = branch.kind == 'T' ? tap.effective_ratio : 1.0,
      effective_complex_tap = branch.kind == 'T' ? tap.effective_complex : 1.0 + 0.0im,
      dtf_tap_convention = branch.kind == 'T' ? tap.convention : :not_transformer,
      tap_changer_model = tap_model,
      tap_impedance_correction_factor = branch.kind == 'T' ? tap_rx.factor : 1.0,
    )
  end
  ok, msg = validate!(net = net)
  ok || error(msg)
  append!(net.for001Contingencies, outage_label.(case.outages))
  return net
end

end # module DTFImporter

"""
    createNetFromDTFFile(path; baseMVA = 100.0, strict = true, bus_shunt_model = :admittance, transformer_ratio_mode = :neutral_one) -> Net

Read a legacy DTF file and build a Sparlectra `Net` without routing
through MATPOWER. Outage cards are parsed and preserved in `DTFCase` by
`DTFImporter.read_dtf`, but are not executed by this Task-1 importer MVP.
"""
function createNetFromDTFFile(path; baseMVA::Real = 100.0, strict::Bool = true, bus_shunt_model = :admittance, legacy_voltage_level_collapse_230kv::Bool = false, transformer_ratio_mode = :neutral_one, tap_changer_model::Union{Nothing,Symbol} = nothing)::Net
  return DTFImporter.build_net(DTFImporter.read_dtf(path; baseMVA = baseMVA, strict = strict);
    bus_shunt_model = bus_shunt_model,
    legacy_voltage_level_collapse_230kv = legacy_voltage_level_collapse_230kv,
    transformer_ratio_mode = transformer_ratio_mode,
    tap_changer_model = tap_changer_model)
end
