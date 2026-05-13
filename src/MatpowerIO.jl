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

# file: src/MatpowerIO.jl

module MatpowerIO
using LinearAlgebra, Printf

export MatpowerCase, read_case, read_case_m, read_case_julia, build_ybus_matpower, vmva_power_mismatch_stats

using LinearAlgebra

@inline function _angle_delta_deg(calc_deg::Real, ref_deg::Real)::Float64
  return mod(Float64(calc_deg) - Float64(ref_deg) + 180.0, 360.0) - 180.0
end

function _normalize_matpower_shift_unit(unit)::Symbol
  value = Symbol(lowercase(String(unit)))
  value in (:deg, :degree, :degrees) && return :deg
  value in (:rad, :radian, :radians) && return :rad
  error(string("matpower_shift_unit must be \"deg\" or \"rad\" (got ", unit, ")."))
end

function _matpower_shift_radians(raw_shift::Float64; sign::Float64 = 1.0, unit::Symbol = :deg)::Float64
  shift_rad = unit == :rad ? raw_shift : deg2rad(raw_shift)
  return sign * shift_rad
end

import ..Sparlectra: setVmVa!, setNodeType!

"""
Container for a MATPOWER case (case format v2/v2-ish).
All matrices are stored as Float64 matrices, names as Vector{String} when available.
"""
struct MatpowerCase
  name::String
  baseMVA::Float64
  bus::Matrix{Float64}
  gen::Matrix{Float64}
  branch::Matrix{Float64}
  gencost::Union{Nothing,Matrix{Float64}}
  bus_name::Union{Nothing,Vector{String}}
end

"""
	legacy_sort_bus(mpc::MatpowerCase) -> MatpowerCase

Mimics the legacy `casefileparser` behavior:
- Sorts `bus` rows by BUS_I (column 1) ascending.
- Does NOT reorder `gen` or `branch`.
This keeps results compatible with the old importer.
"""
function legacy_sort_bus(mpc::MatpowerCase)::MatpowerCase
  bus = mpc.bus
  perm = sortperm(bus[:, 1])  # sort by BUS_I
  bus_sorted = bus[perm, :]
  return MatpowerCase(mpc.name, mpc.baseMVA, bus_sorted, mpc.gen, mpc.branch, mpc.gencost, mpc.bus_name)
end

"""
	normalize_branch_tap!(branch)

MATPOWER semantics:
- column 9 (TAP) == 0 means "no transformer tap specified" -> treat as 1.0.

We normalize data early so downstream code sees a consistent ratio.
"""
function normalize_branch_tap!(branch::AbstractMatrix{<:Real})
  ncol = size(branch, 2)
  ncol < 9 && return branch
  anz = size(branch, 1)
  @inbounds for e ∈ 1:anz
    if branch[e, 9] == 0.0
      branch[e, 9] = 1.0
    end
  end
  return branch
end

"""
	read_case(path::AbstractString; legacy_compat::Bool=true) -> MatpowerCase

Dispatches based on file extension:
- `.m`  => `read_case_m`
- `.jl` => `read_case_julia` (expects file returns `MatpowerCase` or `NamedTuple`)
"""
function read_case(path::AbstractString; legacy_compat::Bool = true)
  lpath = lowercase(path)
  mpc = if endswith(lpath, ".m")
    read_case_m(path; legacy_compat = legacy_compat)
  elseif endswith(lpath, ".jl")
    read_case_julia(path; legacy_compat = legacy_compat)
  else
    error("Unsupported case file extension. Expected .m or .jl, got: $path")
  end
  return mpc
end

"""
	read_case_julia(path::AbstractString; legacy_compat::Bool=true) -> MatpowerCase

The `.jl` file should either:
(A) return a `MatpowerCase`, OR
(B) return a `NamedTuple` with fields: `baseMVA, bus, gen, branch` and optional `gencost, bus_name, name`.
"""
function read_case_julia(path::AbstractString; legacy_compat::Bool = true)
  obj = include(abspath(path))

  mpc = if obj isa MatpowerCase
    obj
  elseif obj isa NamedTuple
    name = get(obj, :name, splitext(basename(path))[1])
    baseMVA = Float64(obj.baseMVA)
    bus = Matrix{Float64}(obj.bus)
    gen = Matrix{Float64}(obj.gen)
    branch = Matrix{Float64}(obj.branch)
    normalize_branch_tap!(branch)

    gencost = haskey(obj, :gencost) ? (obj.gencost === nothing ? nothing : Matrix{Float64}(obj.gencost)) : nothing
    bus_name = haskey(obj, :bus_name) ? (obj.bus_name === nothing ? nothing : Vector{String}(obj.bus_name)) : nothing
    MatpowerCase(String(name), baseMVA, bus, gen, branch, gencost, bus_name)
  else
    error("Julia case file must return MatpowerCase or NamedTuple, got: $(typeof(obj))")
  end

  return legacy_compat ? legacy_sort_bus(mpc) : mpc
end

"""
	read_case_m(path::AbstractString; legacy_compat::Bool=true) -> MatpowerCase

Best-effort parser for typical MATPOWER case files (*.m).
Parses:
- `mpc.baseMVA`
- `mpc.bus`
- `mpc.gen`
- `mpc.branch`
- optional: `mpc.gencost`
- optional: `mpc.bus_name`
"""
function read_case_m(path::AbstractString; legacy_compat::Bool = true)
  txt = read(path, String)
  name = splitext(basename(path))[1]

  # Remove MATLAB comments (%) but keep line structure
  lines = split(txt, '\n')
  lines_nc = String[]
  for ln in lines
    i = findfirst('%', ln)
    push!(lines_nc, i === nothing ? ln : ln[1:prevind(ln, i)])
  end
  txt = join(lines_nc, "\n")

  baseMVA = parse_baseMVA(txt)

  # fixed widths to match MATPOWER v2 tables
  bus    = parse_matrix_block(txt, "mpc.bus"; ncols = 13)
  gen    = parse_matrix_block(txt, "mpc.gen"; ncols = 21)
  branch = parse_matrix_block(txt, "mpc.branch"; ncols = 13)
  normalize_branch_tap!(branch)

  gencost = try_parse_matrix_block(txt, "mpc.gencost")
  bus_name = try_parse_bus_name(txt)

  mpc = MatpowerCase(name, baseMVA, bus, gen, branch, gencost, bus_name)
  apply_supported_postprocessing!(mpc, txt)
  return legacy_compat ? legacy_sort_bus(mpc) : mpc
end

function apply_supported_postprocessing!(mpc::MatpowerCase, txt::AbstractString)
  bus = mpc.bus
  branch = mpc.branch

  # Support common MATPOWER scripts that store branch r/x in Ohm
  # and convert to per-unit at the end of the case file.
  has_ohm_to_pu =
    occursin("mpc.branch(:, [BR_R BR_X])", txt) &&
    occursin("Vbase^2 / Sbase", txt)
  if has_ohm_to_pu
    base_kV = bus[1, 10]  # BASE_KV column in MATPOWER bus table
    vbase = base_kV * 1e3
    sbase = mpc.baseMVA * 1e6
    zbase = vbase^2 / sbase
    branch[:, 3] ./= zbase  # BR_R
    branch[:, 4] ./= zbase  # BR_X
  end

  # Support common MATPOWER scripts that declare loads in kW/kVAr (or kVA)
  # and convert to MW/MVAr by dividing by 1e3.
  has_pdqd_div_1e3 = occursin("mpc.bus(:, [PD, QD]) = mpc.bus(:, [PD, QD]) / 1e3", txt)
  if has_pdqd_div_1e3
    bus[:, 3] ./= 1e3  # PD
    bus[:, 4] ./= 1e3  # QD
  end

  # Support common conversion from apparent power at fixed pf to PD/QD.
  # Equivalent MATLAB pattern:
  #   pf = ...
  #   mpc.bus(:, QD) = mpc.bus(:, PD) * sin(acos(pf));
  #   mpc.bus(:, PD) = mpc.bus(:, PD) * pf;
  m_pf = match(r"\bpf\s*=\s*([0-9]+(?:\.[0-9]+)?)\s*;", txt)
  has_pf_conversion =
    m_pf !== nothing &&
    occursin("mpc.bus(:, QD) = mpc.bus(:, PD) * sin(acos(pf));", txt) &&
    occursin("mpc.bus(:, PD) = mpc.bus(:, PD) * pf;", txt)
  if has_pf_conversion
    pf = parse(Float64, m_pf.captures[1])
    bus[:, 4] .= bus[:, 3] .* sin(acos(pf))  # QD from PD and pf
    bus[:, 3] .*= pf                          # PD scaled by pf
  end

  return mpc
end

function parse_baseMVA(txt::String)
  m = match(r"mpc\.baseMVA\s*=\s*([0-9]+(?:\.[0-9]+)?)\s*;", txt)
  m === nothing && error("Could not find `mpc.baseMVA = ...;` in MATPOWER file.")
  return parse(Float64, m.captures[1])
end

function try_parse_matrix_block(txt::String, key::String)
  try
    return parse_matrix_block(txt, key)
  catch
    return nothing
  end
end

function parse_matrix_block(txt::String, key::String; ncols::Int = 0)
  re = Regex(replace(key, "." => "\\.") * raw"\s*=\s*\[\s*(.*?)\s*\]\s*;", "s")
  m = match(re, txt)
  m === nothing && error("Could not find matrix block `$key = [ ... ];`")
  body = m.captures[1]

  if ncols > 0
    return parse_numeric_matrix_ncols(body, key; ncols = ncols)
  else
    return parse_numeric_matrix(body, key)
  end
end

@inline function _matpower_row_tokens(row::AbstractString)
  s = strip(row)
  isempty(s) && return nothing
  return split(s)
end

# best-effort: auto width; IMPORTANT: fill with 0.0 (not NaN)
function parse_numeric_matrix(body::AbstractString, key::AbstractString)
  rows = Vector{Vector{Float64}}()
  maxcols = 0

  for rr in eachsplit(body, ';')
    toks = _matpower_row_tokens(rr)
    toks === nothing && continue

    vals = Vector{Float64}(undef, length(toks))
    @inbounds for j in eachindex(toks)
      vals[j] = parse(Float64, toks[j])
    end

    maxcols = max(maxcols, length(vals))
    push!(rows, vals)
  end

  isempty(rows) && error("Matrix `$key` appears empty or could not be parsed.")

  M = zeros(Float64, length(rows), maxcols)
  @inbounds for (i, r) in enumerate(rows)
    for j in eachindex(r)
      M[i, j] = r[j]
    end
  end
  return M
end

# fixed width path: pad/truncate to exactly ncols with zeros
function parse_numeric_matrix_ncols(body::AbstractString, key::AbstractString; ncols::Int)
  row_count = 0
  for rr in eachsplit(body, ';')
    _matpower_row_tokens(rr) === nothing || (row_count += 1)
  end
  row_count == 0 && error("Matrix `$key` appears empty or could not be parsed.")

  M = zeros(Float64, row_count, ncols)
  i = 0
  for rr in eachsplit(body, ';')
    toks = _matpower_row_tokens(rr)
    toks === nothing && continue
    i += 1
    ntok = min(length(toks), ncols)
    @inbounds for j = 1:ntok
      M[i, j] = parse(Float64, toks[j])
    end
  end
  return M
end

function try_parse_bus_name(txt::String)
  re = r"mpc\.bus_name\s*=\s*\{\s*(.*?)\s*\}\s*;"s
  m = match(re, txt)
  m === nothing && return nothing
  body = m.captures[1]

  names = String[]
  for mm in eachmatch(r"'([^']*)'", body)
    push!(names, mm.captures[1])
  end
  isempty(names) ? nothing : names
end

# -----------------------------------------------------------------------------
# MATPOWER Vm/Va comparison helpers
# -----------------------------------------------------------------------------

"""
	_mp_bus_row_index(mpc) -> Dict{Int,Int}

Map MATPOWER BUS_I -> row index in mpc.bus
"""
function _mp_bus_row_index(mpc)
  idx = Dict{Int,Int}()
  @inbounds for r = 1:size(mpc.bus, 1)
    idx[Int(mpc.bus[r, 1])] = r
  end
  return idx
end

"""
	compare_vm_va(net, mpc; show_diff=false, tol_vm=1e-6, tol_va=1e-4, maxlines=20)

Compare net results (node._vm_pu/_va_deg) against MATPOWER reference (mpc.bus VM/VA).
Skips isolated buses (BUS_TYPE==4) by default.
Returns (ok::Bool, stats::NamedTuple).
"""
function compare_vm_va(net, mpc; show_diff::Bool = false, tol_vm::Float64 = 1e-6, tol_va::Float64 = 1e-4, maxlines::Int = 20)
  busrow = _mp_bus_row_index(mpc)

  # --- angle alignment on slack to remove global offset ---
  slack_row = findfirst(r -> Int(mpc.bus[r, 2]) == 3, 1:size(mpc.bus, 1))
  slack_busI = slack_row === nothing ? nothing : Int(mpc.bus[slack_row, 1])

  delta_va = 0.0
  have_delta = false
  if slack_busI !== nothing
    net_slack_k = nothing
    @inbounds for k in eachindex(net.nodeVec)
      busI_k = haskey(net.busOrigIdxDict, k) ? net.busOrigIdxDict[k] : k
      if busI_k == slack_busI
        net_slack_k = k
        break
      end
    end
    if net_slack_k !== nothing
      va_ref_slack  = Float64(mpc.bus[slack_row, 9])
      node_slack    = net.nodeVec[net_slack_k]
      va_calc_slack = (node_slack._va_deg === nothing) ? NaN : Float64(node_slack._va_deg)
      if isfinite(va_ref_slack) && isfinite(va_calc_slack)
        delta_va = _angle_delta_deg(va_calc_slack, va_ref_slack)
        have_delta = true
      end
    end
  end

  # Collect diffs
  diffs = Vector{NamedTuple{(:busIdx, :busI, :vm_ref, :vm_calc, :dvm, :va_ref, :va_calc, :dva, :type_ref),Tuple{Int,Int,Float64,Float64,Float64,Float64,Float64,Float64,Int}}}()

  @inbounds for k in eachindex(net.nodeVec)
    node = net.nodeVec[k]

    # Map to MATPOWER BUS_I if available; else assume k
    busI = haskey(net.busOrigIdxDict, k) ? net.busOrigIdxDict[k] : k
    haskey(busrow, busI) || continue
    r = busrow[busI]

    btype = Int(mpc.bus[r, 2])
    btype == 4 && continue  # skip isolated

    vm_ref = Float64(mpc.bus[r, 8])
    va_ref = Float64(mpc.bus[r, 9])

    vm_calc = (node._vm_pu === nothing) ? NaN : Float64(node._vm_pu)
    va_calc = (node._va_deg === nothing) ? NaN : Float64(node._va_deg)

    va_calc_aligned = have_delta ? (va_calc - delta_va) : va_calc

    dvm = vm_calc - vm_ref
    dva = _angle_delta_deg(va_calc_aligned, va_ref)

    push!(diffs, (busIdx = k, busI = busI, vm_ref = vm_ref, vm_calc = vm_calc, dvm = dvm, va_ref = va_ref, va_calc = va_calc_aligned, dva = dva, type_ref = btype))
  end

  isempty(diffs) && return (false, (msg = "No comparable buses found.",))

  abs_dvm = [abs(d.dvm) for d in diffs]
  abs_dva = [abs(d.dva) for d in diffs]

  max_dvm = maximum(abs_dvm)
  max_dva = maximum(abs_dva)

  # Pass/fail: all within tolerance (ignoring NaNs)
  ok_vm = all(x -> (isnan(x) ? true : x <= tol_vm), abs_dvm)
  ok_va = all(x -> (isnan(x) ? true : x <= tol_va), abs_dva)
  ok = ok_vm && ok_va

  println("\n==================== MATPOWER Vm/Va compare ====================")
  println("Compared buses   : ", length(diffs))
  println("tol_vm / tol_va  : ", tol_vm, " pu / ", tol_va, " deg")
  println("angle alignment  : ", have_delta ? "slack (Δva=$(delta_va) deg)" : "none")
  println("max |dVm|        : ", max_dvm, " pu")
  println("max |dVa|        : ", max_dva, " deg")
  println("status           : ", ok ? "OK" : "FAIL")
  println("===============================================================\n")

  if show_diff
    ord = sortperm(1:length(diffs); by = i -> max(abs(diffs[i].dvm), abs(diffs[i].dva)), rev = true)

    println("Top diffs (up to $maxlines lines):")
    println(" busIdx  BUS_I  type   Vm_ref   Vm_calc    dVm      Va_ref   Va_calc    dVa")
    nshow = min(maxlines, length(ord))
    @inbounds for ii = 1:nshow
      d = diffs[ord[ii]]
      @printf(" %6d  %5d   %2d   %7.4f  %7.4f  %+8.5f   %7.3f  %7.3f  %+8.5f\n", d.busIdx, d.busI, d.type_ref, d.vm_ref, d.vm_calc, d.dvm, d.va_ref, d.va_calc, d.dva)
    end
    println()
  end

  return ok, (max_dvm = max_dvm, max_dva = max_dva, n = length(diffs))
end

"""
	build_ybus_matpower(bus, branch, baseMVA) -> Matrix{ComplexF64}

MATPOWER-style Ybus stamping (π-model + tap/shift + bus shunts):
- series r/x
- line charging b split as b/2 on each end
- off-nominal tap ratio and phase shift on from-side
- bus shunts GS/BS added on diagonal: (GS + j*BS)/baseMVA
Optional `matpower_shift_sign` and `matpower_shift_unit` convert branch `SHIFT` before stamping. Defaults preserve MATPOWER degrees/from-side semantics.
"""
function build_ybus_matpower(bus::AbstractMatrix{<:Real}, branch::AbstractMatrix{<:Real}, baseMVA::Float64; matpower_shift_sign::Real = 1.0, matpower_shift_unit = :deg)
  shift_unit = _normalize_matpower_shift_unit(matpower_shift_unit)
  shift_sign = Float64(matpower_shift_sign)
  shift_sign in (-1.0, 1.0) || error(string("matpower_shift_sign must be 1 or -1 (got ", matpower_shift_sign, ")."))

  nbus = size(bus, 1)
  Y = zeros(ComplexF64, nbus, nbus)

  # MATPOWER BUS_I may be non-contiguous (e.g. 1, 2, 9001).
  # Build a mapping BUS_I -> row index so all stamping uses matrix rows.
  busrow = Dict{Int,Int}()
  for i ∈ 1:nbus
    busrow[Int(bus[i, 1])] = i
  end

  # MATPOWER BUS_I may be non-contiguous (e.g. 1, 2, 9001).
  # Build a mapping BUS_I -> row index so all stamping uses matrix rows.
  busrow = Dict{Int,Int}()
  for i ∈ 1:nbus
    busrow[Int(bus[i, 1])] = i
  end

  # MATPOWER BUS_I may be non-contiguous (e.g. 1, 2, 9001).
  # Build a mapping BUS_I -> row index so all stamping uses matrix rows.
  busrow = Dict{Int,Int}()
  for i ∈ 1:nbus
    busrow[Int(bus[i, 1])] = i
  end
  # bus columns (MATPOWER):
  # BUS_I=1, BUS_TYPE=2, PD=3, QD=4, GS=5, BS=6, ...
  # add shunts: (GS + j*BS)/baseMVA
  for i ∈ 1:nbus
    Gs = bus[i, 5]
    Bs = bus[i, 6]
    if Gs != 0.0 || Bs != 0.0
      Y[i, i] += (Gs + 1im * Bs) / baseMVA
    end
  end

  # branch columns (MATPOWER):
  # F_BUS=1, T_BUS=2, BR_R=3, BR_X=4, BR_B=5, RATEA=6, RATEB=7, RATEC=8,
  # TAP=9, SHIFT=10, BR_STATUS=11, ...
  nbranch = size(branch, 1)
  for e ∈ 1:nbranch
    status = (size(branch, 2) >= 11) ? branch[e, 11] : 1.0
    status == 0 && continue

    f_busI = Int(branch[e, 1])
    t_busI = Int(branch[e, 2])
    haskey(busrow, f_busI) || continue
    haskey(busrow, t_busI) || continue
    f = busrow[f_busI]
    t = busrow[t_busI]
    r = branch[e, 3]
    x = branch[e, 4]
    b = branch[e, 5]

    z = r + 1im * x
    y = inv(z)

    # line charging b -> j*b/2 at both ends
    ysh = 1im * b / 2

    tap = (size(branch, 2) >= 9) ? branch[e, 9] : 1.0
    shift_raw = (size(branch, 2) >= 10) ? Float64(branch[e, 10]) : 0.0

    # enforce MATPOWER semantics again (safe)
    if tap == 0.0
      tap = 1.0
    end

    shift = _matpower_shift_radians(shift_raw; sign = shift_sign, unit = shift_unit)
    a = tap * cis(shift)  # complex tap on from side

    # stamping with off-nominal tap
    Y[f, f] += (y + ysh) / (a * conj(a))
    Y[t, t] += (y + ysh)
    Y[f, t] += -y / conj(a)
    Y[t, f] += -y / a
  end

  return Y
end

"""
    vmva_power_mismatch_stats(mpc::MatpowerCase) -> NamedTuple

Checks whether `mpc.bus[:, 8:9]` (VM/VA) are internally consistent with MATPOWER
power balance equations for this case data.

Important PF semantics:
- Active-power equation is enforced for PQ/PV buses (not slack).
- Reactive-power equation is enforced for PQ buses only.

Returns maxima of residuals in p.u. and MW/MVar for those enforced equation sets,
plus equation counts.
"""
function vmva_power_mismatch_stats(mpc::MatpowerCase; matpower_shift_sign::Real = 1.0, matpower_shift_unit = :deg)
  size(mpc.bus, 2) >= 9 || return (ok = false, msg = "mpc.bus has no VM/VA columns")

  busrow = bus_row_index(mpc)
  nbus = size(mpc.bus, 1)

  # Build Sspec from bus loads + in-service generator setpoints.
  Pinj = zeros(Float64, nbus)
  Qinj = zeros(Float64, nbus)

  @inbounds for r = 1:nbus
    Pinj[r] -= mpc.bus[r, 3] / mpc.baseMVA
    Qinj[r] -= mpc.bus[r, 4] / mpc.baseMVA
  end

  ngen = size(mpc.gen, 1)
  @inbounds for g = 1:ngen
    # GEN_STATUS column 8 in MATPOWER
    if size(mpc.gen, 2) >= 8 && mpc.gen[g, 8] <= 0.0
      continue
    end
    busI = Int(mpc.gen[g, 1])
    haskey(busrow, busI) || continue
    r = busrow[busI]
    Pinj[r] += mpc.gen[g, 2] / mpc.baseMVA
    Qinj[r] += mpc.gen[g, 3] / mpc.baseMVA
  end

  Sspec = ComplexF64.(Pinj, Qinj)
  Y = build_ybus_matpower(mpc.bus, mpc.branch, mpc.baseMVA; matpower_shift_sign = matpower_shift_sign, matpower_shift_unit = matpower_shift_unit)

  Vm = mpc.bus[:, 8]
  Va = mpc.bus[:, 9] .* (pi / 180)
  V = Vm .* cis.(Va)
  Scalc = V .* conj.(Y * V)
  mis = Scalc .- Sspec

  p_rows = Int[]  # PQ + PV
  q_rows = Int[]  # PQ only
  @inbounds for r = 1:nbus
    btype = Int(mpc.bus[r, 2])
    if btype == 1
      push!(p_rows, r)
      push!(q_rows, r)
    elseif btype == 2
      push!(p_rows, r)
    elseif btype == 3 || btype == 4
      # slack / isolated: no enforced equation here for this check
    end
  end

  isempty(p_rows) && return (ok = false, msg = "No enforceable P equations (PQ/PV) to check")
  isempty(q_rows) && return (ok = false, msg = "No enforceable Q equations (PQ) to check")

  p_abs = abs.(real.(mis[p_rows]))
  q_abs = abs.(imag.(mis[q_rows]))

  max_p_pu = maximum(p_abs)
  max_q_pu = maximum(q_abs)
  return (ok = true, n_p = length(p_rows), n_q = length(q_rows), max_p_mis_pu = max_p_pu, max_q_mis_pu = max_q_pu, max_p_mis_MW = max_p_pu * mpc.baseMVA, max_q_mis_MVar = max_q_pu * mpc.baseMVA)
end

# -----------------------------------------------------------------------------
# Voltage helpers (MATPOWER semantics)
# -----------------------------------------------------------------------------

"""
    has_vm_va(mpc::MatpowerCase) -> Bool

True if mpc.bus has VM/VA columns (8/9) and at least one finite entry.
"""
function has_vm_va(mpc::MatpowerCase)::Bool
  size(mpc.bus, 2) >= 9 || return false
  vm = mpc.bus[:, 8]
  va = mpc.bus[:, 9]
  return any(isfinite, vm) && any(isfinite, va)
end

"""
    bus_row_index(mpc::MatpowerCase) -> Dict{Int,Int}

Map MATPOWER BUS_I -> row index in mpc.bus.
"""
function bus_row_index(mpc::MatpowerCase)
  idx = Dict{Int,Int}()
  @inbounds for r = 1:size(mpc.bus, 1)
    idx[Int(mpc.bus[r, 1])] = r
  end
  return idx
end

# -----------------------------------------------------------------------------
# MATPOWER bus voltage import policy
# -----------------------------------------------------------------------------

"""
    apply_matpower_bus_voltage!(net, mpc; flatstart, verbose=0)

Policy:
- flatstart=true  -> do NOT write bus VM/VA into net (true flat start)
- flatstart=false -> write bus VM/VA as *initial guess* (only if node has no meaningful voltage yet)
Never warn; debug-only on conflicts.
"""
function apply_matpower_bus_voltage!(net, mpc::MatpowerIO.MatpowerCase; flatstart::Bool, verbose::Int = 0)
  flatstart && return nothing
  MatpowerIO.has_vm_va(mpc) || return nothing

  busrow = MatpowerIO.bus_row_index(mpc)

  @inbounds for k = 1:length(net.nodeVec)
    node = net.nodeVec[k]

    # Map to MATPOWER BUS_I if available (preferred), else assume same index
    busI = haskey(net.busOrigIdxDict, k) ? net.busOrigIdxDict[k] : k
    haskey(busrow, busI) || continue
    r = busrow[busI]

    btype = Int(mpc.bus[r, 2])
    btype == 4 && continue  # MATPOWER isolated

    vm = Float64(mpc.bus[r, 8])
    va = Float64(mpc.bus[r, 9])

    # Only set if missing/unset
    old_vm = node._vm_pu
    old_va = node._va_deg

    if old_vm === nothing || old_vm <= 0.0
      Sparlectra.setVmVa!(node = node, vm_pu = vm, va_deg = va)
    else
      # debug-only (no warnings)
      if verbose > 1 && (abs(Float64(old_vm) - vm) > 1e-6 || (old_va !== nothing && abs(Float64(old_va) - va) > 1e-6))
        @debug "MATPOWER bus VM/VA ignored at busIdx=$k (BUS_I=$busI): already have vm=$(old_vm),va=$(old_va); mpc has vm=$vm,va=$va"
      end
    end
  end

  return nothing
end

function apply_mp_isolated_buses!(net, mpc; verbose::Int = 0)
  size(mpc.bus, 2) >= 2 || return nothing
  busrow = bus_row_index(mpc)

  for k = 1:length(net.nodeVec)
    busI = haskey(net.busOrigIdxDict, k) ? net.busOrigIdxDict[k] : k
    haskey(busrow, busI) || continue
    r = busrow[busI]
    btype = Int(mpc.bus[r, 2])
    if btype == 4
      # mark isolated: keep consistent with your markIsolatedBuses! semantics
      setNodeType!(net.nodeVec[k], "Isolated")
      setVmVa!(node = net.nodeVec[k], vm_pu = 0.0, va_deg = 0.0)
      if !(k in net.isoNodes)
        push!(net.isoNodes, k)
      end
      verbose > 1 && @debug "MATPOWER marks BUS_I=$busI as isolated -> busIdx=$k"
    end
  end
  return nothing
end

function apply_mp_bus_vmva_init!(net, mpc; flatstart::Bool, verbose::Int = 0)
  flatstart && return nothing
  size(mpc.bus, 2) >= 9 || return nothing

  busrow = bus_row_index(mpc)
  for k = 1:length(net.nodeVec)
    node = net.nodeVec[k]

    # skip isolated (already handled)
    (k in net.isoNodes) && continue

    busI = haskey(net.busOrigIdxDict, k) ? net.busOrigIdxDict[k] : k
    haskey(busrow, busI) || continue
    r = busrow[busI]

    vm = Float64(mpc.bus[r, 8])
    va = Float64(mpc.bus[r, 9])

    # only write if missing/unset
    if node._vm_pu === nothing || node._vm_pu <= 0.0
      setVmVa!(node = node, vm_pu = vm, va_deg = va)
    elseif verbose > 1 && (abs(Float64(node._vm_pu) - vm) > 1e-6)
      @debug "Bus init Vm differs at busIdx=$k: keep $(node._vm_pu), ignore mpc $vm"
    end
  end
  return nothing
end

end # module MatpowerIO
