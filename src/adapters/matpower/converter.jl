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

# file: src/adapters/matpower/converter.jl
# purpose: the MATPOWER adapter (adapter task stage 3a): convert_case turns
#          a parsed MatpowerCase into the typed SCFCase, so the network is
#          constructed by build_net like every other format. The conversion
#          mirrors createNetFromMatPowerCase's interpretation rules field
#          by field (classification, conventions, setpoint sources, order),
#          and every SI value is emitted through _scf_stable so the
#          file-unit round trip reproduces the per-unit value bit for bit.

struct MatpowerAdapter <: FormatAdapter end

"""
    MatpowerAdapterOptions

The adapter-scope options of the MATPOWER conversion (design decision D4):
the import conventions plus the model settings that are consumed WHILE the
case is interpreted (the tap-changer impedance correction rewrites R/X at
conversion time). Field names follow the configuration keys.
"""
Base.@kwdef struct MatpowerAdapterOptions
  shift_sign::Float64 = 1.0
  shift_unit::Symbol = :deg
  ratio::Symbol = :normal
  pv_voltage_source::Symbol = :gen_vg
  pv_voltage_mismatch_tol_pu::Float64 = 1e-4
  apply_bus_names::Bool = false
  apply_branch_names::Bool = false
  apply_branch_kind::Bool = false
  matpower_dcline_mode::Symbol = :pf_injections
  import_for001_contingencies::Bool = true
  enable_pq_gen_controllers::Bool = true
  tap_changer_model::Symbol = :ideal
  f_nom::Float64 = 50.0
  log::Bool = false
end

detect(::Type{MatpowerAdapter}, path::AbstractString)::Bool = begin
  ext = lowercase(splitext(String(path))[2])
  ext in (".m", ".jl") || return false
  isfile(path) || return false
  head = open(io -> String(read(io, 4096)), String(path), "r")
  return occursin("mpc.", head) || occursin("MatpowerCase", head)
end

options_type(::MatpowerAdapter) = MatpowerAdapterOptions

"""
    import_net(::MatpowerAdapter, mpc, opts::MatpowerAdapterOptions; kwargs...) -> Net

The MATPOWER importer of the adapter contract (task_import_direct):
builds the network DIRECTLY from the parsed case;
`createNetFromMatPowerCase` is its implementation and keeps its exported
name and signature. The opts cover the adapter scope; run-scope knobs
(bus shunt model, flatstart, cooldown, preallocation) ride through as
keywords exactly as the import dispatch passes them.
"""
function import_net(::MatpowerAdapter, mpc, opts::MatpowerAdapterOptions; kwargs...)::Net
  return createNetFromMatPowerCase(
    mpc = mpc,
    log = opts.log,
    enable_pq_gen_controllers = opts.enable_pq_gen_controllers,
    matpower_shift_sign = opts.shift_sign,
    matpower_shift_unit = opts.shift_unit,
    matpower_ratio = opts.ratio,
    tap_changer_model = opts.tap_changer_model,
    matpower_pv_voltage_source = opts.pv_voltage_source,
    matpower_pv_voltage_mismatch_tol_pu = opts.pv_voltage_mismatch_tol_pu,
    apply_bus_names = opts.apply_bus_names,
    apply_branch_names = opts.apply_branch_names,
    apply_branch_kind = opts.apply_branch_kind,
    import_for001_contingencies = opts.import_for001_contingencies,
    matpower_dcline_mode = opts.matpower_dcline_mode,
    kwargs...,
  )
end

"""
    matpower_adapter_options(cfg::SparlectraConfig) -> MatpowerAdapterOptions

The adapter options an effective run configuration implies; the run path
builds them here so the conversion and the historical direct import read
the same keys.
"""
function matpower_adapter_options(cfg::SparlectraConfig; log::Bool = false)
  mat = cfg.matpower
  return MatpowerAdapterOptions(
    shift_sign = Float64(mat.shift_sign),
    shift_unit = _normalize_matpower_shift_unit(mat.shift_unit),
    ratio = _normalize_matpower_ratio_mode(mat.ratio),
    pv_voltage_source = MatpowerIO._normalize_pv_voltage_source(mat.pv_voltage_source),
    pv_voltage_mismatch_tol_pu = mat.pv_voltage_mismatch_tol_pu,
    apply_bus_names = mat.apply_bus_names,
    apply_branch_names = mat.apply_branch_names,
    apply_branch_kind = mat.apply_branch_kind,
    matpower_dcline_mode = mat.matpower_dcline_mode,
    import_for001_contingencies = mat.import_for001_contingencies,
    enable_pq_gen_controllers = mat.enable_pq_gen_controllers,
    tap_changer_model = cfg.model.tap_changer_model,
    log = log,
  )
end

# stable SI emission helpers: the value written is the one the READER turns
# back into exactly the per-unit/internal value the direct import produced
_mpsi_watt(v::Float64) = _scf_stable(v, x -> x * 1.0e6, x -> x / 1.0e6)
_mpsi_volt(vn_kv::Float64) = _scf_stable(vn_kv, x -> x * 1.0e3, x -> x / 1.0e3)
_mpsi_deg_as_rad(deg::Float64) = _scf_stable(deg, deg2rad, rad2deg)

function convert_case(::MatpowerAdapter, mpc, opts::MatpowerAdapterOptions)::SCFCase
  busDict, genDict, branchDict = _createDict()
  BUS_I = busDict["bus"]; BUS_TYPE = busDict["type"]; PD = busDict["Pd"]; QD = busDict["Qd"]; GS = busDict["Gs"]; BS = busDict["Bs"]
  BUS_AREA = busDict["area"]; VM = busDict["Vm"]; VA = busDict["Va"]; BASE_KV = busDict["baseKV"]; BUS_ZONE = busDict["zone"]; VMAX = busDict["Vmax"]; VMIN = busDict["Vmin"]
  GEN_BUS = genDict["bus"]; PG = genDict["Pg"]; QG = genDict["Qg"]; QMAX = genDict["Qmax"]; QMIN = genDict["Qmin"]; VG = genDict["Vg"]
  GEN_STATUS = genDict["status"]; PMAX = genDict["Pmax"]; PMIN = genDict["Pmin"]; APF = genDict["apf"]
  F_BUS = branchDict["fbus"]; T_BUS = branchDict["tbus"]; BR_R = branchDict["r"]; BR_X = branchDict["x"]; BR_B = branchDict["b"]
  RATE_A = branchDict["rateA"]; TAP = branchDict["ratio"]; SHIFT = branchDict["angle"]; BR_STATUS = branchDict["status"]

  name = hasproperty(mpc, :name) ? String(getproperty(mpc, :name)) : "mpc"
  baseMVA = Float64(getproperty(mpc, :baseMVA))
  busData = convert(Matrix{Float64}, getproperty(mpc, :bus))
  genData = convert(Matrix{Float64}, getproperty(mpc, :gen))
  brData = convert(Matrix{Float64}, getproperty(mpc, :branch))
  nbus = size(busData, 1)
  nbranch = size(brData, 1)
  s_base_va = baseMVA * 1.0e6
  w_nom = 2.0 * pi * opts.f_nom

  pv_voltage_rows = MatpowerIO.pv_voltage_reference_rows(mpc; matpower_pv_voltage_source = opts.pv_voltage_source, tol = opts.pv_voltage_mismatch_tol_pu, warn = opts.log)
  pv_vset_by_bus = Dict(row.busI => row.imported_vset for row in pv_voltage_rows)
  opts.matpower_dcline_mode in (:reject_active, :ignore_inactive, :pf_injections, :paired_control) || throw(ArgumentError("matpower_dcline_mode must be one of :reject_active, :ignore_inactive, :pf_injections, :paired_control."))
  opts.matpower_dcline_mode in (:pf_injections, :paired_control) || MatpowerIO.assert_no_active_dcline(mpc)
  bus_name_by_orig = _matpower_bus_names(mpc, busData, BUS_I; apply_bus_names = opts.apply_bus_names)
  bus_original_name_by_orig = Dict{Int,String}()
  raw_bus_names = _matpower_metadata_vector(mpc, :bus_name)
  if raw_bus_names !== nothing && length(raw_bus_names) == nbus
    for (row_index, row) in enumerate(eachrow(busData))
      orig = raw_bus_names[row_index]
      isempty(orig) || (bus_original_name_by_orig[Int(row[BUS_I])] = orig)
    end
  end
  branch_kind_overrides = _matpower_branch_kind_overrides(mpc, nbranch; apply_branch_kind = opts.apply_branch_kind)
  branch_names = _matpower_metadata_vector(mpc, :branch_name)
  branch_names_valid = branch_names !== nothing && length(branch_names) == nbranch
  transformer_loss_by_row = Dict{Int,NamedTuple}()
  for entry in _matpower_sparlectra_transformer_losses(mpc)
    hasproperty(entry, :branch_row) || continue
    transformer_loss_by_row[Int(entry.branch_row)] = entry
  end
  skip_tap_rx_correction = _matpower_sparlectra_tap_changer_model_marker(mpc) == "impedance_correction"

  case = SCFCase()
  spar = case.sparlectra
  spar.meta["case_name"] = name
  spar.meta["s_base"] = _scf_stable(baseMVA, x -> x * 1.0e6, x -> x / 1.0e6)
  spar.meta["f_nom"] = opts.f_nom
  spar.meta["source_format"] = "matpower"
  extra = spar.extra
  data = case.data

  next_id = Ref(0)
  nid() = (next_id[] += 1)

  mp_bus_type = Dict{Int,Int}()
  node_id_by_orig = Dict{Int,Int}()
  u_rated_by_orig = Dict{Int,Float64}()
  zb_by_orig = Dict{Int,Float64}()
  slack_orig = 0
  isolated_ids = Int[]
  start_nodes = spar.start_state.nodes
  spar.start_state.source = "matpower_case"

  # deferred load emission: the direct import creates the load prosumer of a
  # bus INSIDE the bus loop, so the appliance order is loads (bus order),
  # then generators, then dcline injections; element_index records it
  load_rows = Tuple{Int,Int,Float64,Float64}[]  # (node_id, orig_bus, p, q)

  for (row_index, row) in enumerate(eachrow(busData))
    btype = Int(row[BUS_TYPE])
    orig = Int(row[BUS_I])
    mp_bus_type[orig] = btype
    btype == 3 && slack_orig == 0 && (slack_orig = orig)
    raw_vn = Float64(row[BASE_KV])
    vn_kv = raw_vn <= 0.0 ? 1.0 : raw_vn
    id = nid()
    node_id_by_orig[orig] = id
    u_rated = _mpsi_volt(vn_kv)
    u_rated_by_orig[orig] = u_rated
    zb_by_orig[orig] = u_rated^2 / s_base_va
    push!(data.node, SCFNodeRow(id = id, u_rated = u_rated))
    e = Dict{String,Any}("name" => bus_name_by_orig[orig], "bus_index" => row_index, "source_index" => orig, "external_id" => bus_name_by_orig[orig])
    e["area"] = Int(row[BUS_AREA])
    e["zone"] = Int(row[BUS_ZONE])
    vmin = Float64(row[VMIN])
    vmax = Float64(row[VMAX])
    vmin == 0.9 || (e["vmin_pu"] = vmin)
    vmax == 1.1 || (e["vmax_pu"] = vmax)
    haskey(bus_original_name_by_orig, orig) && (e["original_name"] = bus_original_name_by_orig[orig])
    extra[string(id)] = e
    # start values (D10): the case's own operating point, with the PV and
    # slack setpoint override the direct import applies via setVmVa!
    vm = Float64(row[VM])
    vm = vm <= 0.0 ? 1.0 : vm
    va = Float64(row[VA])
    if btype == 4
      push!(isolated_ids, id)
      vm = 0.0
      va = 0.0
    end
    start_nodes[id] = SCFStartNode(vm_pu = vm, va_deg = va)
    pLoad = Float64(row[PD])
    qLoad = Float64(row[QD])
    (pLoad != 0.0 || qLoad != 0.0) && push!(load_rows, (id, orig, pLoad, qLoad))
    pShunt = Float64(row[GS])
    qShunt = Float64(row[BS])
    if pShunt != 0.0 || qShunt != 0.0
      sid = nid()
      zb_sh = zb_by_orig[orig]
      # the reader's arithmetic forms MW/MVar and divides by baseMVA again;
      # stabilize against exactly that composite, like the SCF writer does
      to_pu = v -> (v * zb_sh * baseMVA) / baseMVA
      g_pu = pShunt / baseMVA
      b_pu = qShunt / baseMVA
      push!(data.shunt, SCFShuntRow(id = sid, node = id, status = 1, g1 = _scf_stable(g_pu, v -> v / zb_sh, to_pu), b1 = _scf_stable(b_pu, v -> v / zb_sh, to_pu)))
      extra[string(sid)] = Dict{String,Any}("name" => string("SHUNT_", orig), "shunt_index" => length(data.shunt))
    end
  end

  # --- branches: same classification rules as the direct import ----------
  branch_index = 0
  for (branch_row_index, row) in enumerate(eachrow(brData))
    fbus = Int(row[F_BUS])
    tbus = Int(row[T_BUS])
    haskey(node_id_by_orig, fbus) || (@warn "bus $(fbus) not found, branch ignored."; continue)
    haskey(node_id_by_orig, tbus) || (@warn "bus $(tbus) not found, branch ignored."; continue)
    r_pu = Float64(row[BR_R])
    x_pu = Float64(row[BR_X])
    b_pu = Float64(row[BR_B])
    ratedS_raw = Float64(row[RATE_A])
    ratedS = ratedS_raw > 0.0 ? ratedS_raw : Inf
    tap_raw = Float64(row[TAP])
    tap_specified = tap_raw != 0.0
    ratio = _matpower_import_ratio(tap_raw; mode = opts.ratio)
    angle = _matpower_shift_degrees(Float64(row[SHIFT]); sign = opts.shift_sign, unit = opts.shift_unit)
    status = Int(row[BR_STATUS])
    vn_from = u_rated_by_orig[fbus] / 1.0e3
    vn_to = u_rated_by_orig[tbus] / 1.0e3
    heuristic_is_line = ((!tap_specified && angle == 0.0) && (vn_from == vn_to))
    loss_meta = get(transformer_loss_by_row, branch_row_index, nothing)
    kind = branch_kind_overrides === nothing ? :unknown : branch_kind_overrides[branch_row_index]
    kind = loss_meta !== nothing ? :transformer : kind
    isLine = kind === :line ? true : kind === :transformer ? false : heuristic_is_line

    id = nid()
    branch_index += 1
    zb = zb_by_orig[tbus]
    e = Dict{String,Any}("name" => string(isLine ? "ACL_" : "T2W_", branch_row_index), "branch_index" => branch_index)
    tap_correction_model = nothing
    tap_correction_factor = 1.0
    if isLine
      push!(data.line, SCFLineRow(
        id = id,
        from_node = node_id_by_orig[fbus],
        to_node = node_id_by_orig[tbus],
        from_status = status,
        to_status = status,
        r1 = _scf_stable(r_pu, v -> v * zb, v -> v / zb),
        x1 = _scf_stable(x_pu, v -> v * zb, v -> v / zb),
        c1 = _scf_stable(b_pu, v -> (v / zb) / w_nom, v -> (v * w_nom) * zb),
        tan1 = 0.0,
        i_n = isfinite(ratedS) ? ratedS * 1.0e6 / (sqrt(3.0) * u_rated_by_orig[tbus]) : nothing,
      ))
    else
      tap_rx = skip_tap_rx_correction ? (r_pu = r_pu, x_pu = x_pu, factor = 1.0) : calcTapCorrectedRX(r_pu = r_pu, x_pu = x_pu, tap_changer_model = opts.tap_changer_model, ratio = ratio)
      tap_correction_model = skip_tap_rx_correction ? :impedance_correction : opts.tap_changer_model
      tap_correction_factor = tap_rx.factor
      g_loss = loss_meta !== nothing && hasproperty(loss_meta, :g_pu) ? Float64(loss_meta.g_pu) : 0.0
      b_eff = b_pu
      if loss_meta !== nothing
        b_loss = hasproperty(loss_meta, :b_pu) ? Float64(loss_meta.b_pu) : 0.0
        b_loss != 0.0 && b_eff == 0.0 && (b_eff = b_loss)
      end
      push!(data.generic_branch, SCFGenericBranchRow(
        id = id,
        from_node = node_id_by_orig[fbus],
        to_node = node_id_by_orig[tbus],
        from_status = status,
        to_status = status,
        r1 = _scf_stable(tap_rx.r_pu, v -> v * zb, v -> v / zb),
        x1 = _scf_stable(tap_rx.x_pu, v -> v * zb, v -> v / zb),
        g1 = _scf_stable(g_loss, v -> v / zb, v -> v * zb),
        b1 = _scf_stable(b_eff, v -> v / zb, v -> v * zb),
        k = ratio,
        theta = _mpsi_deg_as_rad(angle),
        sn = isfinite(ratedS) ? _mpsi_watt(ratedS) : nothing,
      ))
      loss_meta !== nothing && (e["transformer_loss_g_pu"] = g_loss)
    end
    # an unlimited rating (MATPOWER rateA = 0 arrives as Inf) has no PGM
    # slot; the namespaced meta sentinel restores Inf exactly like a
    # Sparlectra-written file
    isfinite(ratedS) || (e["meta"] = Dict{String,Any}("sn_mva" => "inf"))
    # the branch metadata record of the direct import, for report parity
    if branch_names_valid || branch_kind_overrides !== nothing || loss_meta !== nothing || !isLine
      meta = Dict{String,Any}("orig_index" => branch_row_index, "orig_kind" => String(kind), "tap_impedance_correction_factor" => tap_correction_factor)
      branch_names_valid && branch_names[branch_row_index] !== nothing && (meta["orig_name"] = branch_names[branch_row_index])
      tap_correction_model === nothing || (meta["tap_changer_model"] = String(tap_correction_model))
      e["matpower_branch"] = meta
    end
    extra[string(id)] = e
  end

  # --- appliances in the direct import's element order --------------------
  element_index = 0
  participation = SCFParticipationEntry[]
  mFak = 10.0
  for (node_id, orig, pLoad, qLoad) in load_rows
    id = nid()
    element_index += 1
    qMax = min(abs(mFak * qLoad), baseMVA)
    pMax = min(abs(mFak * pLoad), baseMVA)
    push!(data.sym_load, SCFApplianceRow(id = id, node = node_id, status = 1, type = 0, p_specified = _mpsi_watt(pLoad), q_specified = _mpsi_watt(qLoad)))
    extra[string(id)] = Dict{String,Any}(
      "name" => string("LOAD_", orig),
      "element_index" => element_index,
      # the NET's own component vocabulary (what the direct import stores
      # as cTyp), because the build restores exactly this string
      "prosumption_type" => "Load",
      "max_p_mw" => pMax,
      "min_p_mw" => -pMax,
      "max_q_mvar" => qMax,
      "min_q_mvar" => -qMax,
    )
  end

  regulator_rows = Tuple{Int,Int,Float64,Float64,Float64}[]  # (gen_id, node_id, u_ref, qmin, qmax)
  for row in eachrow(genData)
    Int(row[GEN_STATUS]) < 1 && continue
    orig = Int(row[GEN_BUS])
    node_id = node_id_by_orig[orig]
    btype = get(mp_bus_type, orig, 1)
    id = nid()
    element_index += 1
    pGen = Float64(row[PG])
    qGen = Float64(row[QG])
    qMax = Float64(row[QMAX])
    qMin = Float64(row[QMIN])
    pMax = Float64(row[PMAX])
    pMin = Float64(row[PMIN])
    vm_pu = get(pv_vset_by_bus, orig, Float64(row[VG]))
    push!(data.sym_gen, SCFApplianceRow(id = id, node = node_id, status = 1, type = 0, p_specified = _mpsi_watt(pGen), q_specified = _mpsi_watt(qGen)))
    e = Dict{String,Any}("name" => string("GEN_", orig), "element_index" => element_index, "prosumption_type" => "Generator")
    isfinite(pMax) && (e["max_p_mw"] = pMax)
    isfinite(pMin) && (e["min_p_mw"] = pMin)
    if btype == 2
      push!(regulator_rows, (id, node_id, vm_pu, qMin, qMax))
      e["regulated"] = true
    elseif btype == 3
      # the slack machine: the direct import's constructor auto-regulates
      # any machine that carries a voltage setpoint, so the case states
      # exactly that (setpoint plus flag, no regulator row: the bus is the
      # reference, not a PV-regulated node)
      e["vm_pu"] = vm_pu
      e["regulated"] = true
      isfinite(qMax) && (e["max_q_mvar"] = qMax)
      isfinite(qMin) && (e["min_q_mvar"] = qMin)
    else
      isfinite(qMax) && (e["max_q_mvar"] = qMax)
      isfinite(qMin) && (e["min_q_mvar"] = qMin)
    end
    slack_orig == orig && (e["reference_pri"] = true)
    # PQ generator limits become constant P(U)/Q(U) controllers on build,
    # exactly like the direct import (stage-3a carrier; the flag makes the
    # derivation explicit so existing SCF files stay untouched)
    opts.enable_pq_gen_controllers && btype == 1 && (e["pq_gen_controller"] = true)
    apf_raw = length(row) >= APF ? Float64(row[APF]) : 0.0
    (isfinite(apf_raw) && apf_raw > 0.0) && push!(participation, SCFParticipationEntry(object = id, factor = apf_raw))
    extra[string(id)] = e
  end

  # --- dclines as fixed terminal injections (pf_injections semantics) ----
  dcline_records = Dict{String,Any}[]
  if opts.matpower_dcline_mode in (:pf_injections, :paired_control)
    dcline = _matpower_metadata_vector(mpc, :dcline)
    if dcline !== nothing
      for r in axes(dcline, 1)
        size(dcline, 2) >= 3 && dcline[r, 3] == 0.0 && continue
        size(dcline, 2) >= 4 || throw(ArgumentError("MATPOWER dcline row $(r) must have at least 4 columns (F_BUS, T_BUS, BR_STATUS, PF) for :pf_injections mode."))
        fbus = Int(dcline[r, 1])
        tbus = Int(dcline[r, 2])
        haskey(node_id_by_orig, fbus) || throw(ArgumentError("MATPOWER dcline row $(r) references unknown F_BUS $(fbus)."))
        haskey(node_id_by_orig, tbus) || throw(ArgumentError("MATPOWER dcline row $(r) references unknown T_BUS $(tbus)."))
        f_bus_type = get(mp_bus_type, fbus, 1)
        t_bus_type = get(mp_bus_type, tbus, 1)
        pf = Float64(dcline[r, 4])
        input_pt = _matpower_dcline_col(dcline, r, 5, "PT"; required = size(dcline, 2) < 17)
        loss0 = _matpower_dcline_col(dcline, r, 16, "LOSS0"; default = 0.0)
        loss1 = _matpower_dcline_col(dcline, r, 17, "LOSS1"; default = 0.0)
        pt = size(dcline, 2) >= 17 ? pf - (loss0 + loss1 * pf) : input_pt
        qf = _matpower_dcline_col(dcline, r, 6, "QF")
        qt = _matpower_dcline_col(dcline, r, 7, "QT")
        vf = size(dcline, 2) >= 8 ? Float64(dcline[r, 8]) : nothing
        vt = size(dcline, 2) >= 9 ? Float64(dcline[r, 9]) : nothing
        qminf = size(dcline, 2) >= 12 ? Float64(dcline[r, 12]) : nothing
        qmaxf = size(dcline, 2) >= 13 ? Float64(dcline[r, 13]) : nothing
        qmint = size(dcline, 2) >= 14 ? Float64(dcline[r, 14]) : nothing
        qmaxt = size(dcline, 2) >= 15 ? Float64(dcline[r, 15]) : nothing
        # terminal voltage control mirrors the direct import's rule; the
        # Isolated check there can only fire for bus type 4, which is what
        # f_bus_type/t_bus_type already encode at conversion time
        from_ctrl = vf !== nothing && f_bus_type != 4 && f_bus_type != 3
        to_ctrl = vt !== nothing && t_bus_type != 4 && t_bus_type != 3
        entry = Dict{String,Any}("orig_index" => r, "from_bus" => fbus, "to_bus" => tbus, "status" => size(dcline, 2) >= 3 ? Float64(dcline[r, 3]) : 1.0, "pf_mw" => pf, "input_pt_mw" => input_pt, "effective_pt_mw" => pt, "loss0_mw" => loss0, "loss1" => loss1, "qf_mvar" => qf, "qt_mvar" => qt)
        vf === nothing || (entry["vf_pu"] = vf)
        vt === nothing || (entry["vt_pu"] = vt)
        qminf === nothing || (entry["qminf_mvar"] = qminf)
        qmaxf === nothing || (entry["qmaxf_mvar"] = qmaxf)
        qmint === nothing || (entry["qmint_mvar"] = qmint)
        qmaxt === nothing || (entry["qmaxt_mvar"] = qmaxt)
        entry["from_voltage_controlled"] = from_ctrl
        entry["to_voltage_controlled"] = to_ctrl
        for (side, bus_orig, p, q, qmin, qmax, vmv, ctrl) in (("from", fbus, -pf, qf, qminf, qmaxf, f_bus_type == 4 ? nothing : vf, from_ctrl), ("to", tbus, pt, qt, qmint, qmaxt, t_bus_type == 4 ? nothing : vt, to_ctrl))
          id = nid()
          element_index += 1
          push!(data.sym_gen, SCFApplianceRow(id = id, node = node_id_by_orig[bus_orig], status = 1, type = 0, p_specified = _mpsi_watt(p), q_specified = _mpsi_watt(q)))
          e = Dict{String,Any}("name" => string("DCLINE_", r, "_", side), "element_index" => element_index, "prosumption_type" => "Generator")
          qmax === nothing || !isfinite(qmax) || (e["max_q_mvar"] = qmax)
          qmin === nothing || !isfinite(qmin) || (e["min_q_mvar"] = qmin)
          if vmv !== nothing
            # the direct import's constructor auto-regulates any terminal
            # that carries a voltage value, whether or not the dcline row
            # declared it voltage controlled
            e["vm_pu"] = Float64(vmv)
            e["regulated"] = true
          end
          extra[string(id)] = e
          entry[string(side, "_element_index")] = element_index
        end
        push!(dcline_records, entry)
      end
    end
  end
  isempty(dcline_records) || (spar.components.matpower_dcline = dcline_records)
  isempty(dcline_records) || opts.matpower_dcline_mode !== :paired_control || (spar.meta["matpower_dcline_paired_control"] = true)

  # regulator rows in appliance order, ids after every appliance
  for (gen_id, node_id, u_ref, qmin, qmax) in regulator_rows
    rid = nid()
    row = SCFVoltageRegulatorRow(id = rid, regulated_object = gen_id, status = 1, u_ref = u_ref)
    isfinite(qmin) && (row.q_min = _mpsi_watt(qmin))
    isfinite(qmax) && (row.q_max = _mpsi_watt(qmax))
    push!(data.voltage_regulator, row)
  end

  # --- links and tap nameplates from the Sparlectra extension -------------
  sp_links = _matpower_sparlectra_links(mpc)
  if sp_links !== nothing
    for row in eachrow(sp_links)
      fb = Int(row[1])
      tb = Int(row[2])
      haskey(node_id_by_orig, fb) || throw(ArgumentError("mpc.sparlectra.links references unknown fbus $(fb)."))
      haskey(node_id_by_orig, tb) || throw(ArgumentError("mpc.sparlectra.links references unknown tbus $(tb)."))
      id = nid()
      push!(data.link, SCFLinkRow(id = id, from_node = node_id_by_orig[fb], to_node = node_id_by_orig[tb], from_status = Int(row[3]), to_status = Int(row[3])))
    end
  end
  sp_taps = _matpower_sparlectra_tap_changers(mpc)
  if sp_taps !== nothing
    by_index = Dict{Int,Int}()
    for r in data.line
      by_index[_scf_int(extra[string(r.id)]["branch_index"], "branch_index")] = r.id
    end
    for r in data.generic_branch
      by_index[_scf_int(extra[string(r.id)]["branch_index"], "branch_index")] = r.id
    end
    for row in eachrow(sp_taps)
      k = Int(row[1])
      haskey(by_index, k) || throw(ArgumentError("mpc.sparlectra.tap_changers references unknown branch $(k)."))
      bid = by_index[k]
      controllers = SCFTapControllerRow[]
      Float64(row[2]) > 0.0 && push!(controllers, SCFTapControllerRow(index = 1, alpha_deg = 0.0, step = Float64(row[2]), pos = Float64(row[5]), pos_min = Float64(row[3]), pos_max = Float64(row[4])))
      psi = size(sp_taps, 2) >= 10 ? Float64(row[10]) : 0.0
      dustep = size(sp_taps, 2) >= 11 ? Float64(row[11]) : 0.0
      # same nameplate contradiction check as applyTapNameplate! on the
      # direct path (found by the oracle-endpoint switch: the adapter
      # silently preferred phase_du_step when both grids were declared)
      (Float64(row[6]) > 0.0 && dustep > 0.0) && throw(ArgumentError("mpc.sparlectra.tap_changers branch $(k): phase_step_deg and phase_du_step are exclusive (a changer has ONE mechanical grid)."))
      if Float64(row[6]) > 0.0 || dustep > 0.0
        ctrl = SCFTapControllerRow(index = 2, alpha_deg = psi, pos = Float64(row[9]), pos_min = Float64(row[7]), pos_max = Float64(row[8]))
        if dustep > 0.0
          ctrl.step = dustep
        else
          ctrl.step_deg = Float64(row[6])
        end
        push!(controllers, ctrl)
      end
      isempty(controllers) && continue
      tid = nid()
      # the branch row keeps the NEUTRAL position (mpc TAP/SHIFT); the
      # nameplate application moves the live tap off it, exactly like the
      # direct import's applyTapNameplate! call
      push!(spar.components.tap_changer, SCFTapChangerRow(id = tid, branch = bid, side = "from", controllers = controllers, control = Dict{String,Any}("mode" => "fixed")))
    end
  end

  # --- roles, contingencies, provenance -----------------------------------
  if slack_orig != 0
    spar.roles.slack = SCFSlackRole(mode = isempty(participation) ? "single" : "distributed", nodes = Int[node_id_by_orig[slack_orig]], participation = participation)
  end
  isempty(isolated_ids) || (spar.roles.isolated_nodes = isolated_ids)
  if opts.import_for001_contingencies && hasproperty(mpc, :for001_contingencies) && getproperty(mpc, :for001_contingencies) !== nothing
    # the metadata carries outage LABELS (strings), one row per label
    rows = Dict{String,Any}[Dict{String,Any}("label" => String(c)) for c in getproperty(mpc, :for001_contingencies)]
    isempty(rows) || (spar.components.for001_contingency = rows)
  end
  return case
end
