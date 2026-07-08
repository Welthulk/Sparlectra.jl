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

# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 01.10.2023
# file: src/createnet_powermat.jl
# helper
#! format: off


function _createDict()

  busKeys = ["bus", "type", "Pd", "Qd", "Gs", "Bs", "area", "Vm", "Va", "baseKV", "zone", "Vmax", "Vmin"]
  busDict = Dict{String,Int}()
  for (idx, key) in enumerate(busKeys)
    busDict[key] = idx
  end

  genKeys = ["bus", "Pg", "Qg", "Qmax", "Qmin", "Vg", "mBase", "status", "Pmax", "Pmin", "Pc1", "Pc2", "Qc1min", "Qc1max", "Qc2min", "Qc2max", "ramp_agc", "ramp_10", "ramp_30", "ramp_q", "apf"]
  genDict = Dict{String,Int}()
  for (idx, key) in enumerate(genKeys)    
    genDict[key] = idx
  end

  branchKey = ["fbus", "tbus", "r", "x", "b", "rateA", "rateB", "rateC", "ratio", "angle", "status", "angmin", "angmax"]
  branchDict = Dict{String,Int}()
  for (idx, key) in enumerate(branchKey)    
    branchDict[key] = idx
  end

  return busDict, genDict, branchDict
end

function _normalize_matpower_shift_unit(unit)::Symbol
  value = Symbol(lowercase(String(unit)))
  value in (:deg, :degree, :degrees) && return :deg
  value in (:rad, :radian, :radians) && return :rad
  error(string("matpower_shift_unit must be \"deg\" or \"rad\" (got ", unit, ")."))
end

function _matpower_shift_degrees(raw_shift::Float64; sign::Float64 = 1.0, unit::Symbol = :deg)::Float64
  shift_deg = unit == :rad ? rad2deg(raw_shift) : raw_shift
  return sign * shift_deg
end

function _normalize_matpower_ratio_mode(mode)::Symbol
  value = Symbol(lowercase(String(mode)))
  value === :normal && return :normal
  value === :reciprocal && return :reciprocal
  error(string("matpower_ratio must be \"normal\" or \"reciprocal\" (got ", mode, ")."))
end

function _matpower_import_ratio(raw_ratio::Float64; mode::Symbol = :normal)::Float64
  ratio = raw_ratio == 0.0 ? 1.0 : raw_ratio
  return mode === :reciprocal ? inv(ratio) : ratio
end


#! format: on

"""
    createNetFromMatPowerCase(; mpc, log=false, flatstart=false) -> Net

Builds a Sparlectra `Net` from a MATPOWER-like container `mpc`.

`mpc` can be either:
- a `NamedTuple` with fields `name, baseMVA, bus, gen, branch` (optionally `gencost, bus_name`)
- or a struct with the same field names (e.g. `MatpowerCase`)

All matrices are expected in MATPOWER v2 column conventions.

`bus_shunt_model` controls how MATPOWER bus `Gs`/`Bs` values are represented:
`"admittance"` stamps them into Ybus (default), while
`"voltage_dependent_injection"` keeps them out of Ybus and evaluates their
|V|²-dependent contribution in the rectangular mismatch path.

`matpower_shift_sign` and `matpower_shift_unit` control how MATPOWER branch
`SHIFT` values are converted before they are stored as Sparlectra transformer
phase shifts. Defaults preserve MATPOWER convention: `SHIFT` is in degrees,
positive on the branch from side. PEGASE-style data sets may require
`matpower_shift_unit = "rad"` and/or `matpower_shift_sign = -1`.

`matpower_ratio` controls MATPOWER branch `TAP` import. The default `"normal"`
uses the `TAP` value directly (with MATPOWER `0` treated as `1`). Set
`matpower_ratio = "reciprocal"` when an input data set stores the inverse tap
ratio expected by Sparlectra.
"""
function _apply_matpower_reference_override!(net::Net, slack_orig_idx::Int, bus_idx_by_orig::Dict{Int,Int}; reference_vm_pu::Union{Nothing,Float64} = nothing, reference_va_deg::Union{Nothing,Float64} = nothing)
  slack_orig_idx != 0 || return nothing
  node_idx = get(bus_idx_by_orig, slack_orig_idx, 0)
  node_idx != 0 || return nothing

  node = net.nodeVec[node_idx]
  vm = isnothing(reference_vm_pu) ? (node._vm_pu === nothing ? 1.0 : Float64(node._vm_pu)) : reference_vm_pu
  va = isnothing(reference_va_deg) ? (node._va_deg === nothing ? 0.0 : Float64(node._va_deg)) : reference_va_deg
  setVmVa!(node = node, vm_pu = vm, va_deg = va)
  return nothing
end

function _matpower_metadata_vector(mpc, field::Symbol)
  return hasproperty(mpc, field) ? getproperty(mpc, field) : nothing
end

function _matpower_bus_names(mpc, busData, BUS_I::Int; apply_bus_names::Bool)
  names = Dict{Int,String}()
  for row in eachrow(busData)
    names[Int(row[BUS_I])] = string(Int(row[BUS_I]))
  end
  apply_bus_names || return names
  bus_name = _matpower_metadata_vector(mpc, :bus_name)
  if bus_name === nothing || length(bus_name) != size(busData, 1)
    @warn "MATPOWER bus_name metadata missing or length-mismatched; falling back to numeric BUS_I names." expected = size(busData, 1) actual = bus_name === nothing ? 0 : length(bus_name)
    return names
  end
  seen = Set{String}()
  for (row_index, row) in enumerate(eachrow(busData))
    name = bus_name[row_index]
    !isempty(name) || throw(ArgumentError("MATPOWER bus_name at row $(row_index) is empty."))
    name in seen && throw(ArgumentError("Duplicate MATPOWER bus_name `$(name)`; cannot import ambiguous bus names."))
    push!(seen, name)
    names[Int(row[BUS_I])] = name
  end
  return names
end

function _normalize_matpower_branch_kind(kind)::Symbol
  k = uppercase(strip(String(kind)))
  k in ("L", "LINE", "ACL") && return :line
  k in ("T", "TRAFO", "TRANSFORMER", "2WT") && return :transformer
  return :unknown
end

function _matpower_branch_kind_overrides(mpc, nbranch::Int; apply_branch_kind::Bool)
  kinds = _matpower_metadata_vector(mpc, :branch_kind)
  if !apply_branch_kind
    return nothing
  elseif kinds === nothing || length(kinds) != nbranch
    @warn "MATPOWER branch_kind metadata missing or length-mismatched; using electrical branch classification heuristic." expected = nbranch actual = kinds === nothing ? 0 : length(kinds)
    return nothing
  end
  return [_normalize_matpower_branch_kind(k) for k in kinds]
end

function _matpower_sparlectra_transformer_losses(mpc)
  hasproperty(mpc, :sparlectra) || return NamedTuple[]
  ext = getproperty(mpc, :sparlectra)
  ext === nothing && return NamedTuple[]
  hasproperty(ext, :format_version) && ext.format_version != 1 && throw(ArgumentError("Unsupported Sparlectra MATPOWER extension format_version=$(ext.format_version)."))
  return hasproperty(ext, :transformer_losses) ? collect(ext.transformer_losses) : NamedTuple[]
end

function _matpower_dcline_col(dcline::AbstractMatrix, row::Integer, col::Integer, name::AbstractString; required::Bool = false, default::Float64 = 0.0)
  if size(dcline, 2) >= col
    return Float64(dcline[row, col])
  end
  required && throw(ArgumentError("MATPOWER dcline row $(row) requires column $(col) ($(name)) for :pf_injections mode."))
  return default
end

function _matpower_dcline_terminal_voltage_control(net::Net, bus_idx::Int, bus_type::Int, vm::Union{Nothing,Float64})::Bool
  vm === nothing && return false
  bus_type == 4 && return false
  bus_type == 3 && return false
  (1 <= bus_idx <= length(net.nodeVec)) || return false
  getNodeType(net.nodeVec[bus_idx]) == Isolated && return false
  return true
end

function createNetFromMatPowerCase(; mpc, log::Bool=false, flatstart::Bool=false, cooldown::Int = 0, q_hyst_pu::Float64 = 0.0, enable_pq_gen_controllers::Bool = true, bus_shunt_model = :admittance, matpower_shift_sign::Real = 1.0, matpower_shift_unit = :deg, matpower_ratio = :normal, reference_vm_pu::Union{Nothing,Float64} = nothing, reference_va_deg::Union{Nothing,Float64} = nothing, matpower_pv_voltage_source = :gen_vg, matpower_pv_voltage_mismatch_tol_pu::Float64 = 1e-4, preallocate_network::Symbol = :auto, preallocate_min_buses::Int = 1000, apply_bus_names::Bool = false, apply_branch_names::Bool = false, apply_branch_kind::Bool = false, import_for001_contingencies::Bool = true, matpower_dcline_mode::Symbol = :reject_active, profile::Union{Nothing,AbstractDict}=nothing)::Net
  # Small logger helper
  pInfo(msg::String) = (log ? (@info msg) : nothing)

  # --- Extract fields from NamedTuple / struct uniformly ---
  name    = hasproperty(mpc, :name)    ? getproperty(mpc, :name)    : "mpc"
  baseMVA = hasproperty(mpc, :baseMVA) ? Float64(getproperty(mpc, :baseMVA)) : error("mpc.baseMVA missing")
  busData = hasproperty(mpc, :bus)     ? Matrix{Float64}(getproperty(mpc, :bus)) : error("mpc.bus missing")
  genData = hasproperty(mpc, :gen)     ? Matrix{Float64}(getproperty(mpc, :gen)) : error("mpc.gen missing")
  brData  = hasproperty(mpc, :branch)  ? Matrix{Float64}(getproperty(mpc, :branch)) : error("mpc.branch missing")

  # --- Legacy-compatible dicts (same as your old importer) ---
  busDict, genDict, branchDict = _createDict()
  BUS_I = busDict["bus"]; BUS_TYPE = busDict["type"]; PD = busDict["Pd"]; QD = busDict["Qd"]; GS = busDict["Gs"]; BS = busDict["Bs"]
  BUS_AREA = busDict["area"]; VM = busDict["Vm"]; VA = busDict["Va"]; BASE_KV = busDict["baseKV"]; BUS_ZONE = busDict["zone"]; VMAX = busDict["Vmax"]; VMIN = busDict["Vmin"]
  GEN_BUS = genDict["bus"]; PG = genDict["Pg"]; QG = genDict["Qg"]; QMAX = genDict["Qmax"]; QMIN = genDict["Qmin"]; VG = genDict["Vg"]
  MBASE = genDict["mBase"]; GEN_STATUS = genDict["status"]; PMAX = genDict["Pmax"]; PMIN = genDict["Pmin"]
  F_BUS = branchDict["fbus"]; T_BUS = branchDict["tbus"]; BR_R = branchDict["r"]; BR_X = branchDict["x"]; BR_B = branchDict["b"]
  RATE_A = branchDict["rateA"]; TAP = branchDict["ratio"]; SHIFT = branchDict["angle"]; BR_STATUS = branchDict["status"]

  mp_bus_type = Dict{Int,Int}()
  bus_row_by_i = Dict{Int,Int}()
  gen_rows_by_bus_i = Dict{Int,Vector{Int}}()
  branch_rows_by_from_to = Dict{Tuple{Int,Int},Vector{Int}}()
  nbus = size(busData, 1)
  nbranch = size(brData, 1)
  ngen = size(genData, 1)
  sizehint!(mp_bus_type, nbus)
  sizehint!(bus_row_by_i, nbus)
  @inbounds for r in axes(busData, 1)
    bus_i = Int(busData[r, BUS_I])
    mp_bus_type[bus_i] = Int(busData[r, BUS_TYPE])
    bus_row_by_i[bus_i] = r
  end
  @inbounds for g in axes(genData, 1)
    bus_i = Int(genData[g, GEN_BUS])
    push!(get!(gen_rows_by_bus_i, bus_i, Int[]), g)
  end
  @inbounds for e in axes(brData, 1)
    key = (Int(brData[e, F_BUS]), Int(brData[e, T_BUS]))
    push!(get!(branch_rows_by_from_to, key, Int[]), e)
  end

  shunt_model = normalize_bus_shunt_model(bus_shunt_model)
  shift_unit = _normalize_matpower_shift_unit(matpower_shift_unit)
  shift_sign = Float64(matpower_shift_sign)
  shift_sign in (-1.0, 1.0) || error(string("matpower_shift_sign must be 1 or -1 (got ", matpower_shift_sign, ")."))
  ratio_mode = _normalize_matpower_ratio_mode(matpower_ratio)
  pv_voltage_source = MatpowerIO._normalize_pv_voltage_source(matpower_pv_voltage_source)
  pv_voltage_rows = MatpowerIO.pv_voltage_reference_rows(mpc; matpower_pv_voltage_source = pv_voltage_source, tol = matpower_pv_voltage_mismatch_tol_pu, warn = log)
  pv_vset_by_bus = Dict(row.busI => row.imported_vset for row in pv_voltage_rows)
  matpower_dcline_mode in (:reject_active, :ignore_inactive, :pf_injections) || throw(ArgumentError("matpower_dcline_mode must be one of :reject_active, :ignore_inactive, :pf_injections."))
  matpower_dcline_mode === :pf_injections || MatpowerIO.assert_no_active_dcline(mpc)
  bus_name_by_orig = _matpower_bus_names(mpc, busData, BUS_I; apply_bus_names = apply_bus_names)
  bus_original_name_by_orig = Dict{Int,String}()
  raw_bus_names = _matpower_metadata_vector(mpc, :bus_name)
  if raw_bus_names !== nothing && length(raw_bus_names) == nbus
    for (row_index, row) in enumerate(eachrow(busData))
      name = raw_bus_names[row_index]
      isempty(name) || (bus_original_name_by_orig[Int(row[BUS_I])] = name)
    end
  end
  branch_kind_overrides = _matpower_branch_kind_overrides(mpc, nbranch; apply_branch_kind = apply_branch_kind)
  branch_names = _matpower_metadata_vector(mpc, :branch_name)
  branch_names_valid = branch_names !== nothing && length(branch_names) == nbranch
  apply_branch_names && !branch_names_valid && @warn "MATPOWER branch_name metadata missing or length-mismatched; branch metadata names will not be attached." expected = nbranch actual = branch_names === nothing ? 0 : length(branch_names)
  sparlectra_transformer_losses = _matpower_sparlectra_transformer_losses(mpc)
  transformer_loss_by_row = Dict{Int,NamedTuple}()
  for entry in sparlectra_transformer_losses
    hasproperty(entry, :branch_row) || continue
    transformer_loss_by_row[Int(entry.branch_row)] = entry
  end
  !isempty(transformer_loss_by_row) && @info "Sparlectra MATPOWER transformer-loss metadata found" transformer_loss_records = length(transformer_loss_by_row)

  if log
    @info "Creating new Net: $(name) with baseMVA=$(baseMVA), flatstart=$(flatstart)"
  end
  
  myNet = Net(name = String(name), baseMVA = baseMVA, flatstart = flatstart, cooldown_iters = cooldown, q_hyst_pu = q_hyst_pu, bus_shunt_model = shunt_model)
  do_prealloc = preallocate_network === :on || (preallocate_network === :auto && nbus >= preallocate_min_buses)
  if do_prealloc
    sizehint!(myNet.nodeVec, nbus)
    sizehint!(myNet.busDict, nbus)
    sizehint!(myNet.busOrigIdxDict, nbus)
    sizehint!(myNet.branchVec, nbranch)
    sizehint!(myNet.linesAC, nbranch)
    sizehint!(myNet.trafos, nbranch)
    sizehint!(myNet.prosumpsVec, nbus + ngen)
    sizehint!(myNet.shuntVec, nbus)
    sizehint!(myNet.shuntDict, nbus)
  end
  bus_idx_by_orig = Dict{Int,Int}()
  sizehint!(bus_idx_by_orig, nbus)

  # --- Find slack bus index from BUS_TYPE==3 (MATPOWER) ---
  slackIdx = 0
  for row in eachrow(busData)
    btype = Int(row[BUS_TYPE])
    if btype == 3
      slackIdx = Int(row[BUS_I])
      break
    end
  end

  # --- Buses + Loads + Shunts (same semantics as your existing importer) ---
  mFak = 10.0

  for row in eachrow(busData)
    btype = Int(row[BUS_TYPE])
    kIdx  = Int(row[BUS_I])     # original bus number
    busName = bus_name_by_orig[kIdx]

    raw_vn = Float64(row[BASE_KV])
    vn_kv  = raw_vn <= 0.0 ? 1.0 : raw_vn

    va_deg = Float64(row[VA])
    vm_pu  = Float64(row[VM])
    vm_pu  = vm_pu <= 0.0 ? 1.0 : vm_pu

    pLoad  = Float64(row[PD])
    qLoad  = Float64(row[QD])

    pShunt = Float64(row[GS])
    qShunt = Float64(row[BS])

    zone = Int(row[BUS_ZONE])
    area = Int(row[BUS_AREA])
    vmin = Float64(row[VMIN])
    vmax = Float64(row[VMAX])


    addBus!(
      net = myNet,
      busName = busName,
      vn_kV = vn_kv,
      vm_pu = vm_pu,
      va_deg = va_deg,
      vmin_pu = vmin,
      vmax_pu = vmax,
      isAux = false,
      oBusIdx = kIdx,
      zone = zone,
      area = area,
    )
    # fix 04.02.2026: MATPOWER -> p.u. admittance components:
    #=
    pShunt_pu = pShunt / baseMVA
    qShunt_pu = qShunt / baseMVA
    if pShunt != 0.0 || qShunt != 0.0
      addShunt!(net=myNet, busName=busName, pShunt=pShunt_pu, qShunt=qShunt_pu)
    end
    =#
    # new 04.02.2026: directly add shunt as bus shunt
    pShunt = Float64(row[GS])
    qShunt = Float64(row[BS])

    bus_idx_by_orig[kIdx] = length(myNet.nodeVec)
    if haskey(bus_original_name_by_orig, kIdx)
      myNet.busOriginalNameDict[bus_idx_by_orig[kIdx]] = bus_original_name_by_orig[kIdx]
    end

    if pShunt != 0.0 || qShunt != 0.0
      addShuntMatpower!(net=myNet, busName=busName, Gs=pShunt, Bs=qShunt, bus_shunt_model = shunt_model)
    end
    
    if pLoad != 0.0 || qLoad != 0.0
      qMax = min(abs(mFak * qLoad), baseMVA)
      qMin = -qMax
      pMax = min(abs(mFak * pLoad), baseMVA)
      pMin = -pMax
      addProsumer!(
        net = myNet,
        busName = busName,
        type = "ENERGYCONSUMER",
        p = pLoad,
        q = qLoad,
        pMin = pMin,
        pMax = pMax,
        qMin = qMin,
        qMax = qMax,
        defer_bus_type_refresh = true,
      )
    end
  end

  # --- Branches (line vs transformer decision same as old importer) ---
  for (branch_row_index, row) in enumerate(eachrow(brData))
    fbus_orig = Int(row[F_BUS])
    tbus_orig = Int(row[T_BUS])
    from_idx = get(bus_idx_by_orig, fbus_orig, 0)
    to_idx = get(bus_idx_by_orig, tbus_orig, 0)

    from_idx != 0 || (@warn "bus $(fbus_orig) not found, branch ignored."; continue)
    to_idx != 0 || (@warn "bus $(tbus_orig) not found, branch ignored."; continue)

    r_pu = Float64(row[BR_R])
    x_pu = Float64(row[BR_X])
    b_pu = Float64(row[BR_B])

    ratedS_raw = Float64(row[RATE_A])
    ratedS = ratedS_raw > 0.0 ? ratedS_raw : Inf

    tap_raw = Float64(row[TAP])
    tap_specified = tap_raw != 0.0
    ratio = _matpower_import_ratio(tap_raw; mode = ratio_mode)

    # MATPOWER defines SHIFT on the branch from side. Sparlectra uses the
    # same PI tap placement; changing PST orientation for ratio=1 is therefore
    # equivalent to changing the phase-shift sign. For off-nominal ratios,
    # reversing from/to also moves the ratio tap and is intentionally not done.
    angle  = _matpower_shift_degrees(Float64(row[SHIFT]); sign = shift_sign, unit = shift_unit)
    status = Int(row[BR_STATUS])

    vn_from = getNodeVn(myNet.nodeVec[from_idx])
    vn_to   = getNodeVn(myNet.nodeVec[to_idx])

    # MATPOWER uses TAP == 0 for ordinary lines and a non-zero TAP value for
    # transformers. Preserve explicit nominal-tap transformers (`TAP == 1`) so
    # case300-style branches such as BUS_I 196 -> 2040 are not demoted to lines.
    heuristic_is_line = ((!tap_specified && angle == 0.0) && (vn_from == vn_to))
    kind = branch_kind_overrides === nothing ? :unknown : branch_kind_overrides[branch_row_index]
    isLine = kind === :line ? true : kind === :transformer ? false : heuristic_is_line

    if isLine
      _addPIModelACLine_by_idx!(
        net = myNet,
        from = from_idx,
        to = to_idx,
        r_pu = r_pu,
        x_pu = x_pu,
        b_pu = b_pu,
        status = status,
        ratedS = ratedS,
      )
    else
      _addPIModelTrafo_by_idx!(
        net = myNet,
        from = from_idx,
        to = to_idx,
        r_pu = r_pu,
        x_pu = x_pu,
        b_pu = b_pu,
        status = status,
        ratedS = ratedS,
        ratio = ratio,
        shift_deg = angle,
      )
    end
    imported_branch_index = length(myNet.branchVec)
    loss_meta = get(transformer_loss_by_row, branch_row_index, nothing)
    if loss_meta !== nothing
      g_pu = hasproperty(loss_meta, :g_pu) ? Float64(loss_meta.g_pu) : 0.0
      b_loss_pu = hasproperty(loss_meta, :b_pu) ? Float64(loss_meta.b_pu) : 0.0
      expected_gs_each = g_pu * baseMVA / 2
      expected_bs_each = b_loss_pu * baseMVA / 2
      from_bus_name = bus_name_by_orig[fbus_orig]
      to_bus_name = bus_name_by_orig[tbus_orig]
      if expected_gs_each != 0.0 || expected_bs_each != 0.0
        from_row = bus_row_by_i[fbus_orig]
        to_row = bus_row_by_i[tbus_orig]
        represented =
          abs(Float64(busData[from_row, GS])) + 1e-9 >= abs(expected_gs_each) &&
          abs(Float64(busData[to_row, GS])) + 1e-9 >= abs(expected_gs_each)
        if represented
          @info "Sparlectra transformer-loss metadata restored without adding shunts; standard MATPOWER Gs/Bs already carries the equivalent approximation." branch_row = branch_row_index g_pu = g_pu
        else
          addShuntMatpower!(net = myNet, busName = from_bus_name, Gs = expected_gs_each, Bs = expected_bs_each, bus_shunt_model = shunt_model)
          addShuntMatpower!(net = myNet, busName = to_bus_name, Gs = expected_gs_each, Bs = expected_bs_each, bus_shunt_model = shunt_model)
          @info "Applied Sparlectra transformer-loss metadata as equivalent terminal bus shunts." branch_row = branch_row_index g_pu = g_pu added_Gs_each = expected_gs_each
        end
      end
    end
    if branch_names_valid || branch_kind_overrides !== nothing || loss_meta !== nothing
      myNet.matpower_branch_metadata[imported_branch_index] = (
        orig_name = branch_names_valid ? branch_names[branch_row_index] : nothing,
        orig_kind = kind,
        orig_index = branch_row_index,
        transformer_loss = loss_meta,
      )
    end
  end

  # --- Generators ---
  pq_gen_controller_count = 0
  for row in eachrow(genData)
    status = Int(row[GEN_STATUS])
    status < 1 && continue

    bus = bus_name_by_orig[Int(row[GEN_BUS])]

    pGen = Float64(row[PG])
    qGen = Float64(row[QG])

    qMax = Float64(row[QMAX])
    qMin = Float64(row[QMIN])
    pMax = Float64(row[PMAX])
    pMin = Float64(row[PMIN])
    vm_pu_gen = Float64(row[VG])
    vm_pu = get(pv_vset_by_bus, Int(row[GEN_BUS]), vm_pu_gen)
    mBase = Float64(row[MBASE])
    btype = get(mp_bus_type, Int(row[GEN_BUS]), 1)

    referencePri = (slackIdx == Int(row[GEN_BUS])) ? bus : nothing
    (mBase != baseMVA) && @debug "generator $(bus) has different mBase than network baseMVA (allowed in MATPOWER)" bus mBase baseMVA

    pu_controller = nothing
    qu_controller = nothing
    if enable_pq_gen_controllers && btype == 1
      p_pu = pGen / baseMVA
      q_pu = qGen / baseMVA
      pu_controller = PUController(make_characteristic([(0.0, p_pu), (2.0, p_pu)]); pmin_MW = pMin, pmax_MW = pMax, sbase_MVA = baseMVA)
      qu_controller = QUController(make_characteristic([(0.0, q_pu), (2.0, q_pu)]); qmin_MVAr = qMin, qmax_MVAr = qMax, sbase_MVA = baseMVA)
      pq_gen_controller_count += 1
    end

    if btype == 2 || btype == 3
      # MATPOWER PV/slack voltage setpoints are configurable because some
      # solved cases store BUS.VM values that differ from GEN.VG setpoints.
      node_idx = get(bus_idx_by_orig, Int(row[GEN_BUS]), 0)
      node_idx != 0 && setVmVa!(node = myNet.nodeVec[node_idx], vm_pu = vm_pu)
    end

    addProsumer!(
      net = myNet,
      busName = bus,
      type = "GENERATOR",
      vm_pu = (btype == 2 || btype == 3) ? vm_pu : nothing,
      referencePri = referencePri,
      p = pGen,
      q = qGen,
      pMax = pMax,
      pMin = pMin,
      qMax = qMax,
      qMin = qMin,
      qu_controller = qu_controller,
      pu_controller = pu_controller,
      isRegulated = (btype == 2),
      defer_bus_type_refresh = true,
    )
  end

  if matpower_dcline_mode === :pf_injections
    dcline = _matpower_metadata_vector(mpc, :dcline)
    if dcline !== nothing
      @inbounds for r in axes(dcline, 1)
        size(dcline, 2) >= 3 && dcline[r, 3] == 0.0 && continue
        size(dcline, 2) >= 4 || throw(ArgumentError("MATPOWER dcline row $(r) must have at least 4 columns (F_BUS, T_BUS, BR_STATUS, PF) for :pf_injections mode."))
        fbus = Int(dcline[r, 1])
        tbus = Int(dcline[r, 2])
        f_name = get(bus_name_by_orig, fbus, nothing)
        t_name = get(bus_name_by_orig, tbus, nothing)
        f_name !== nothing || throw(ArgumentError("MATPOWER dcline row $(r) references unknown F_BUS $(fbus)."))
        t_name !== nothing || throw(ArgumentError("MATPOWER dcline row $(r) references unknown T_BUS $(tbus)."))
        f_idx = bus_idx_by_orig[fbus]
        t_idx = bus_idx_by_orig[tbus]
        f_bus_type = get(mp_bus_type, fbus, 1)
        t_bus_type = get(mp_bus_type, tbus, 1)
        pf = Float64(dcline[r, 4])
        input_pt = _matpower_dcline_col(dcline, r, 5, "PT"; required = size(dcline, 2) < 17)
        loss0 = _matpower_dcline_col(dcline, r, 16, "LOSS0"; default = 0.0)
        loss1 = _matpower_dcline_col(dcline, r, 17, "LOSS1"; default = 0.0)
        pt = input_pt
        if size(dcline, 2) >= 17
          pt = pf - (loss0 + loss1 * pf)
        end
        qf = _matpower_dcline_col(dcline, r, 6, "QF")
        qt = _matpower_dcline_col(dcline, r, 7, "QT")
        vf = size(dcline, 2) >= 8 ? Float64(dcline[r, 8]) : nothing
        vt = size(dcline, 2) >= 9 ? Float64(dcline[r, 9]) : nothing
        qminf = size(dcline, 2) >= 12 ? Float64(dcline[r, 12]) : nothing
        qmaxf = size(dcline, 2) >= 13 ? Float64(dcline[r, 13]) : nothing
        qmint = size(dcline, 2) >= 14 ? Float64(dcline[r, 14]) : nothing
        qmaxt = size(dcline, 2) >= 15 ? Float64(dcline[r, 15]) : nothing
        from_voltage_controlled = _matpower_dcline_terminal_voltage_control(myNet, f_idx, f_bus_type, vf)
        to_voltage_controlled = _matpower_dcline_terminal_voltage_control(myNet, t_idx, t_bus_type, vt)
        from_vm = f_bus_type == 4 ? nothing : vf
        to_vm = t_bus_type == 4 ? nothing : vt
        addProsumer!(net = myNet, busName = f_name, type = "GENERATOR", p = -pf, q = qf, qMin = qminf, qMax = qmaxf, vm_pu = from_vm, isRegulated = from_voltage_controlled, defer_bus_type_refresh = true)
        from_prosumer = length(myNet.prosumpsVec)
        addProsumer!(net = myNet, busName = t_name, type = "GENERATOR", p = pt, q = qt, qMin = qmint, qMax = qmaxt, vm_pu = to_vm, isRegulated = to_voltage_controlled, defer_bus_type_refresh = true)
        to_prosumer = length(myNet.prosumpsVec)
        push!(myNet.matpowerDclineMetadata, (
          orig_index = r,
          from_bus = fbus,
          to_bus = tbus,
          from_bus_name = f_name,
          to_bus_name = t_name,
          status = size(dcline, 2) >= 3 ? Float64(dcline[r, 3]) : 1.0,
          pf_mw = pf,
          input_pt_mw = input_pt,
          effective_pt_mw = pt,
          pt_mw = pt,
          loss0_mw = loss0,
          loss1 = loss1,
          loss_mw = pf - pt,
          qf_mvar = qf,
          qt_mvar = qt,
          vf_pu = vf,
          vt_pu = vt,
          qminf_mvar = qminf,
          qmaxf_mvar = qmaxf,
          qmint_mvar = qmint,
          qmaxt_mvar = qmaxt,
          from_prosumer = from_prosumer,
          to_prosumer = to_prosumer,
          from_voltage_controlled = from_voltage_controlled,
          to_voltage_controlled = to_voltage_controlled,
        ))
      end
      !isempty(myNet.matpowerDclineMetadata) && @info "MATPOWER active mpc.dcline rows imported using toggle_dcline-compatible PF injections" active_dcline_count = length(myNet.matpowerDclineMetadata)
    end
  end

  if import_for001_contingencies && hasproperty(mpc, :for001_contingencies) && getproperty(mpc, :for001_contingencies) !== nothing
    append!(myNet.for001Contingencies, getproperty(mpc, :for001_contingencies))
  end

  refreshBusTypesFromProsumers!(myNet)
  _buildQLimits!(myNet)

  if enable_pq_gen_controllers && pq_gen_controller_count > 0
    @info "MATPOWER import: PQ generator limits mapped to constant P(U)/Q(U) controllers" count = pq_gen_controller_count
  end

  _apply_matpower_reference_override!(myNet, slackIdx, bus_idx_by_orig; reference_vm_pu = reference_vm_pu, reference_va_deg = reference_va_deg)

  log && log_bus_shunt_model(myNet)

  ok, msg = validate!(net = myNet)
  ok || @debug "network is invalid: $msg"
  if profile !== nothing
    profile[:network_construction_preallocated] = do_prealloc
    profile[:network_construction_nbus] = nbus
    profile[:network_construction_nbranch] = nbranch
    profile[:network_construction_ngen] = ngen
    profile[:network_construction_subtimings] = Dict(
      :matpower_data_normalization => 0.0,
      :net_object_creation => 0.0,
      :bus_import => 0.0,
      :branch_import => 0.0,
      :generator_prosumer_import => 0.0,
      :shunt_import => 0.0,
      :pq_generator_controller_setup => 0.0,
      :bus_branch_dictionary_construction => 0.0,
      :validation_post_import_consistency => 0.0,
    )
  end

  return myNet
end

function createNetFromMatPowerFile(; filename::String,
    log::Bool=false,
    flatstart::Union{Nothing,Bool}=nothing,
    cooldown::Int = 0,
    q_hyst_pu::Float64 = 0.0,
    enable_pq_gen_controllers::Union{Nothing,Bool}=nothing,
    bus_shunt_model = nothing,
    matpower_shift_sign::Union{Nothing,Real} = nothing,
    matpower_shift_unit = nothing,
    matpower_ratio = nothing,
    reference_vm_pu::Union{Nothing,Float64} = nothing,
    reference_va_deg::Union{Nothing,Float64} = nothing,
    matpower_pv_voltage_source = nothing,
    matpower_pv_voltage_mismatch_tol_pu::Union{Nothing,Float64} = nothing,
    apply_bus_names::Union{Nothing,Bool} = nothing,
    apply_branch_names::Union{Nothing,Bool} = nothing,
    apply_branch_kind::Union{Nothing,Bool} = nothing,
    import_for001_contingencies::Union{Nothing,Bool} = nothing,
    matpower_dcline_mode = nothing,
    verbose::Int = 0,
    profile::Union{Nothing,AbstractDict}=nothing)::Net

  mat_cfg = matpower_import_config()
  pf_cfg = powerflow_config()
  flatstart = isnothing(flatstart) ? pf_cfg.start_mode.flatstart : flatstart
  enable_pq_gen_controllers = isnothing(enable_pq_gen_controllers) ? mat_cfg.enable_pq_gen_controllers : enable_pq_gen_controllers
  bus_shunt_model = isnothing(bus_shunt_model) ? mat_cfg.bus_shunt_model : bus_shunt_model
  matpower_shift_sign = isnothing(matpower_shift_sign) ? mat_cfg.shift_sign : matpower_shift_sign
  matpower_shift_unit = isnothing(matpower_shift_unit) ? mat_cfg.shift_unit : matpower_shift_unit
  matpower_ratio = isnothing(matpower_ratio) ? mat_cfg.ratio : matpower_ratio
  matpower_pv_voltage_source = isnothing(matpower_pv_voltage_source) ? mat_cfg.pv_voltage_source : matpower_pv_voltage_source
  matpower_pv_voltage_mismatch_tol_pu = isnothing(matpower_pv_voltage_mismatch_tol_pu) ? mat_cfg.pv_voltage_mismatch_tol_pu : matpower_pv_voltage_mismatch_tol_pu
  apply_bus_names = isnothing(apply_bus_names) ? mat_cfg.apply_bus_names : apply_bus_names
  apply_branch_names = isnothing(apply_branch_names) ? mat_cfg.apply_branch_names : apply_branch_names
  apply_branch_kind = isnothing(apply_branch_kind) ? mat_cfg.apply_branch_kind : apply_branch_kind
  import_for001_contingencies = isnothing(import_for001_contingencies) ? mat_cfg.import_for001_contingencies : import_for001_contingencies
  matpower_dcline_mode = isnothing(matpower_dcline_mode) ? mat_cfg.matpower_dcline_mode : matpower_dcline_mode
  preallocate_network = mat_cfg.preallocate_network
  preallocate_min_buses = mat_cfg.preallocate_min_buses

  mpc = MatpowerIO.read_case(filename; legacy_compat=true)

  # Build the network first
  net = createNetFromMatPowerCase(; mpc=mpc, log=log, flatstart=flatstart,
                                  cooldown=cooldown, q_hyst_pu=q_hyst_pu, enable_pq_gen_controllers=enable_pq_gen_controllers,
                                  bus_shunt_model=bus_shunt_model, matpower_shift_sign=matpower_shift_sign,
                                  matpower_shift_unit=matpower_shift_unit, matpower_ratio=matpower_ratio,
                                  reference_vm_pu=reference_vm_pu, reference_va_deg=reference_va_deg,
                                  matpower_pv_voltage_source=matpower_pv_voltage_source,
                                  matpower_pv_voltage_mismatch_tol_pu=matpower_pv_voltage_mismatch_tol_pu,
                                  apply_bus_names=apply_bus_names, apply_branch_names=apply_branch_names,
                                  apply_branch_kind=apply_branch_kind, import_for001_contingencies=import_for001_contingencies,
                                  matpower_dcline_mode=matpower_dcline_mode,
                                  preallocate_network=preallocate_network, preallocate_min_buses=preallocate_min_buses, profile=profile)

  # Always apply MATPOWER isolated flags (BUS_TYPE==4) onto net
  MatpowerIO.apply_mp_isolated_buses!(net, mpc; verbose=verbose)

  # Only if not flatstart: take VM/VA from mpc.bus as initial guess
  MatpowerIO.apply_mp_bus_vmva_init!(net, mpc; flatstart=flatstart, verbose=verbose)

  return net
end
