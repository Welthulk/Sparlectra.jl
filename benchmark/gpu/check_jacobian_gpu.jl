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
#
# file: benchmark/gpu/check_jacobian_gpu.jl
# purpose: proof for the GPU study, stage 1b: the Jacobian, the mismatch, the
#          current I = Ybus V and abs(V) of the device kernels against the
#          CPU code, entry by entry, on the cases of the study at flat start,
#          at the converged state and at perturbed states, single and as a
#          batch. Exit code 0 only if every check passed.
#
# Run: julia --project=benchmark/gpu benchmark/gpu/check_jacobian_gpu.jl [case ...]
# Output CSV: $GPU_STUDY_RUNS/check_jacobian_gpu.csv (default
# benchmark/gpu/results).
#
# The kernels compute I = Ybus V with the unfused complex muladd of
# SparseArrays._spmatmul!. The script first checks which association
# Sparlectra's CPU SpMV (calc_currents) gives on this machine and stops with
# a failure if it is the fused one: the comparison below would then measure
# against another association.
#
# Acceptance per entry, fixed in advance:
#   abs(gpu - cpu) <= 1e-14 * sigma
# with sigma the sum of the absolute values of the terms the CPU adds for
# that entry (sigma = 0 requires gpu == cpu). Reported as well: the plain
# relative error abs(gpu - cpu) / abs(cpu), the absolute error where
# cpu == 0, and the count of bitwise equal entries.

using Sparlectra, CUDA, CUDSS, SparseArrays, LinearAlgebra, Random, Printf, Dates

const EXT = Base.get_extension(Sparlectra, :SparlectraCUDAExt)
EXT === nothing && error("SparlectraCUDAExt not loaded")

const CASES = isempty(ARGS) ? ["case14", "case118", "case1354pegase", "case2869pegase", "case6495rte", "case9241pegase", "case_ACTIVSg10k", "case_ACTIVSg25k", "case_ACTIVSg70k"] : ARGS
const TOL = 1e-14
const ROOT = normpath(joinpath(@__DIR__, "..", ".."))

function out_dir()
  d = get(ENV, "GPU_STUDY_RUNS", "")
  isempty(d) || return d
  d = joinpath(@__DIR__, "results")
  mkpath(d)
  return d
end

bus_type_symbol(t) = t == Sparlectra.Slack ? :Slack : t == Sparlectra.PV ? :PV : :PQ
bus_types_of(net) = [bus_type_symbol(Sparlectra.getNodeType(nd)) for nd in net.nodeVec]

# --- CPU references -----------------------------------------------------------

# I = Ybus * V by rows in ascending column order, the two associations of
# the complex muladd of SparseArrays._spmatmul!
function spmv_unfused(Yt::SparseMatrixCSC{ComplexF64}, V::Vector{ComplexF64})
  n = length(V)
  out = zeros(ComplexF64, n)
  for i = 1:n
    re = 0.0
    im_ = 0.0
    for k in nzrange(Yt, i)
      y = Yt.nzval[k]
      x = V[Yt.rowval[k]]
      re = real(y) * real(x) - (imag(y) * imag(x) - re)
      im_ = real(y) * imag(x) + (imag(y) * real(x) + im_)
    end
    out[i] = ComplexF64(re, im_)
  end
  return out
end

function spmv_fused(Yt::SparseMatrixCSC{ComplexF64}, V::Vector{ComplexF64})
  n = length(V)
  out = zeros(ComplexF64, n)
  for i = 1:n
    re = 0.0
    im_ = 0.0
    for k in nzrange(Yt, i)
      y = Yt.nzval[k]
      x = V[Yt.rowval[k]]
      re = fma(real(y), real(x), -fma(imag(y), imag(x), -re))
      im_ = fma(real(y), imag(x), fma(imag(y), real(x), im_))
    end
    out[i] = ComplexF64(re, im_)
  end
  return out
end

# R_i and M_i: the sums of the absolute values of the terms of real(I_i)
# and imag(I_i)
function row_magnitudes(Yt::SparseMatrixCSC{ComplexF64}, V::Vector{ComplexF64})
  n = length(V)
  R = zeros(n)
  M = zeros(n)
  for i = 1:n, k in nzrange(Yt, i)
    g, h = reim(Yt.nzval[k])
    c, d = reim(V[Yt.rowval[k]])
    R[i] += abs(g * c) + abs(h * d)
    M[i] += abs(g * d) + abs(h * c)
  end
  return R, M
end

# --- entry comparison -----------------------------------------------------------

mutable struct Cmp
  entries::Int
  bitwise::Int
  fails::Int
  max_sigma_err::Float64
  max_rel::Float64
  max_abs_cpu0::Float64
  worst::String
end
Cmp() = Cmp(0, 0, 0, 0.0, 0.0, 0.0, "")

function add!(c::Cmp, gpu::Float64, cpu::Float64, sigma::Float64, label)
  c.entries += 1
  d = abs(gpu - cpu)
  (gpu === cpu || (gpu == cpu && !iszero(cpu))) && (c.bitwise += 1)
  if cpu == 0.0
    c.max_abs_cpu0 = max(c.max_abs_cpu0, d)
  else
    c.max_rel = max(c.max_rel, d / abs(cpu))
  end
  ok = sigma == 0.0 ? gpu == cpu : d <= TOL * sigma
  ok = ok && (isnan(gpu) == isnan(cpu))
  if !ok
    c.fails += 1
    c.fails <= 3 && (c.worst *= string(label(), " gpu=", gpu, " cpu=", cpu, " sigma=", sigma, "; "))
  end
  sigma > 0.0 && (c.max_sigma_err = max(c.max_sigma_err, d / sigma))
  return nothing
end

# sigma of Jacobian entry (row r, column c) and whether it is a PV Q row
function jac_sigma(r, c, nonslack, n, Ybus, V, R, M, types, dP, dQ, cpu, vm_eps)
  i = nonslack[(r + 1) >> 1]
  isP = isodd(r)
  isVr = c <= n - 1
  j = nonslack[isVr ? c : c - (n - 1)]
  a, b = reim(V[i])
  if !isP && types[i] === :PV
    return abs(cpu)
  end
  if j != i
    g, h = reim(Ybus[i, j])
    return (isP == isVr) ? abs(a * g) + abs(b * h) : abs(a * h) + abs(b * g)
  end
  vm = abs(V[i])
  vm_safe = vm > vm_eps ? vm : vm_eps
  dvm = isVr ? a / vm_safe : b / vm_safe
  chain = abs((isP ? dP[i] : dQ[i]) * dvm)
  insorted(i, view(Ybus.rowval, nzrange(Ybus, i))) || return chain
  g, h = reim(Ybus[i, i])
  if isP
    return isVr ? abs(a * g) + abs(b * h) + R[i] + chain : abs(a * h) + abs(b * g) + M[i] + chain
  else
    return isVr ? abs(a * h) + abs(b * g) + M[i] + chain : abs(a * g) + abs(b * h) + R[i] + chain
  end
end

# GPU values (CSR rowptr/colval/vals) against a CPU CSC matrix: pattern
# containment, GPU-only entries exactly zero, every entry by sigma
function compare_jacobian(rowptr, colval, gvals::AbstractVector{Float64}, Jcpu::SparseMatrixCSC{Float64}, nonslack, n, Ybus, V, R, M, types, dP, dQ, vm_eps)
  Ct = copy(transpose(Jcpu))   # CSR of the CPU Jacobian, values moved, not computed
  c = Cmp()
  missing_cpu = 0
  gpu_only_nonzero = 0
  nvar = length(rowptr) - 1
  for r = 1:nvar
    gk = Int(rowptr[r])
    ge = Int(rowptr[r+1]) - 1
    ck = Ct.colptr[r]
    ce = Ct.colptr[r+1] - 1
    while gk <= ge
      col = Int(colval[gk])
      while ck <= ce && Ct.rowval[ck] < col
        missing_cpu += 1
        ck += 1
      end
      if ck <= ce && Ct.rowval[ck] == col
        cpu = Ct.nzval[ck]
        ck += 1
        s = jac_sigma(r, col, nonslack, n, Ybus, V, R, M, types, dP, dQ, cpu, vm_eps)
        add!(c, gvals[gk], cpu, s, () -> "J[$r,$col]")
      else
        gvals[gk] === 0.0 || (gpu_only_nonzero += 1)
      end
      gk += 1
    end
    missing_cpu += ce - ck + 1
  end
  return c, missing_cpu, gpu_only_nonzero
end

function compare_mismatch(Fg::AbstractVector{Float64}, Fc::Vector{Float64}, nonslack, V, S, Vset, types, R, M)
  c = Cmp()
  for (t, i) in enumerate(nonslack)
    a, b = reim(V[i])
    sP = abs(a) * R[i] + abs(b) * M[i] + abs(real(S[i]))
    s2 = types[i] === :PV ? abs(V[i]) + Vset[i] : abs(b) * R[i] + abs(a) * M[i] + abs(imag(S[i]))
    add!(c, Fg[2t-1], Fc[2t-1], sP, () -> "F[$(2t-1)]")
    add!(c, Fg[2t], Fc[2t], s2, () -> "F[$(2t)]")
  end
  return c
end

bitwise_count(a, b) = count(k -> a[k] === b[k], eachindex(a, b))

# which association the CPU SpMV of Sparlectra gives: "unfused", "fused" or
# "neither", by bitwise comparison with the two references
function spmv_association(name)
  p = Sparlectra.ensure_casefile(name * ".m"; to_jl = false)
  net = Sparlectra.createNetFromMatPowerFile(filename = p)
  Ybus = Sparlectra.createYBUS(net = net, sparse = true)
  n = size(Ybus, 1)
  rng = MersenneTwister(118)
  V = ComplexF64.(0.9 .+ 0.2 .* rand(rng, n), 0.3 .* (2 .* rand(rng, n) .- 1))
  Yt = copy(transpose(Ybus))
  Icpu = Sparlectra.calc_currents(Ybus, V)
  nu = bitwise_count(Icpu, spmv_unfused(Yt, V))
  nf = bitwise_count(Icpu, spmv_fused(Yt, V))
  variant = nu == n ? "unfused" : nf == n ? "fused" : "neither"
  return (n = n, unfused = nu, fused = nf, variant = variant)
end

# --- static FMA check -------------------------------------------------------------

# own kernels: plain mul/add/sub.f64 must not occur; fma.rn.f64 only from
# _hypot_cr, 4 per inlined call (the muladd of Base written as fma plus the
# three fma of the correction step). _k_update! is exempt (Base complex
# division and cis).
const EXPECTED_HYPOT_CALLS = Dict(:_k_current_mismatch! => 1, :_k_mismatch_from_current! => 1, :_k_reduce_max! => 0, :_k_jacobian_offdiag! => 0, :_k_jacobian_diag! => 1, :_k_vm! => 1, :_k_qlimit_sample! => 1, :_k_switch_scatter! => 0)

function static_fma_check(Ybus, slack)
  rows = String[]
  ok = true
  for T in (Float64, Float32)
    for (name, f, tt) in EXT._device_kernel_signatures(Ybus, slack; T = T)
      ptx = sprint(io -> CUDA.code_ptx(io, f, tt; kernel = true))
      plain = length(collect(eachmatch(r"\b(mul|add|sub)\.f64\b", ptx)))
      plain32 = length(collect(eachmatch(r"\b(mul|add|sub)\.f32\b", ptx)))
      fmas = length(collect(eachmatch(r"\bfma\.rn\.f64\b", ptx)))
      exempt = name === :_k_update!
      expected = exempt ? -1 : 4 * EXPECTED_HYPOT_CALLS[name]
      pass = exempt || (plain == 0 && plain32 == 0 && fmas == expected)
      ok &= pass
      push!(rows, @sprintf("%-26s T=%-7s plain f64 %3d  plain f32 %3d  fma.rn.f64 %3d (expected %s)  %s", name, T, plain, plain32, fmas, exempt ? "exempt" : string(expected), exempt ? "EXEMPT" : pass ? "ok" : "FAIL"))
    end
  end
  return ok, rows
end

# --- one case ---------------------------------------------------------------------

const QUIET = OutputConfig(logfile_results = :off, console_summary = false, startup_latency_hint = false)

function case_inputs(name)
  p = Sparlectra.ensure_casefile(name * ".m"; to_jl = false)
  start_net = Sparlectra.createNetFromMatPowerFile(filename = p)
  types_start = bus_types_of(start_net)
  r = run_sparlectra(casefile = p, config = SparlectraConfig(output = QUIET))
  r.final_converged || @warn "C1 did not converge" name
  net = r.net
  length(net.nodeVec) == length(start_net.nodeVec) || error("bus count differs between import and solved net")
  [nd.comp.cName for nd in net.nodeVec] == [nd.comp.cName for nd in start_net.nodeVec] || error("bus order differs between import and solved net")
  n = length(net.nodeVec)
  Ybus = Sparlectra.createYBUS(net = net, sparse = true)
  size(Ybus, 1) == n || error("isolated buses: not in the proof")
  Vflat, slack = Sparlectra.initialVrect(net; flatstart = true)
  S = Sparlectra.buildComplexSVec(net)
  Vset = Sparlectra._bus_voltage_setpoints_from_prosumers(net)
  Vflat[slack] = Sparlectra._apply_voltage_magnitude_preserving_angle(Vflat[slack], Vset[slack])
  Vconv = Sparlectra.buildVoltageVector(net)
  types_final = bus_types_of(net)
  return (; p, net, n, Ybus, slack, Vflat, Vconv, S, Vset, types_start, types_final, iterations = length(Sparlectra.rectangular_pf_status(net).mismatch_history))
end

function run_case(name, csv_rows, summary)
  t0 = now()
  ci = case_inputs(name)
  n = ci.n
  slack = ci.slack
  Ybus = ci.Ybus
  Yt = copy(transpose(Ybus))
  nonslack = [i for i = 1:n if i != slack]
  zeros_n = zeros(n)
  rng = MersenneTwister(20261112)
  allpq = [i == slack ? :Slack : :PQ for i = 1:n]
  Vpert = ci.Vconv .+ ComplexF64.(1e-2 .* (2 .* rand(rng, n) .- 1), 1e-2 .* (2 .* rand(rng, n) .- 1))
  Vpert[slack] = ci.Vconv[slack]
  dPr = 0.05 .+ 0.5 .* rand(rng, n)
  dQr = -(0.05 .+ 0.5 .* rand(rng, n))
  states = Any[
    ("flat_start", ci.Vflat, ci.types_start, ci.S, zeros_n, zeros_n),
    ("converged_final_types", ci.Vconv, ci.types_final, ci.S, zeros_n, zeros_n),
    ("converged_start_types", ci.Vconv, ci.types_start, ci.S, zeros_n, zeros_n),
    ("perturbed_all_pq_dinj", Vpert, allpq, ci.S, dPr, dQr),
  ]
  if Sparlectra.has_voltage_dependent_control(ci.net)
    Sc, dPc, dQc = Sparlectra.buildControlledSVec(ci.net, ci.Vconv)
    push!(states, ("converged_controlled_dinj", ci.Vconv, ci.types_final, Sc, dPc, dQc))
  end
  @printf("\n%s: n %d, nnz(Ybus) %d, max Ybus row length %d, C1 iterations %d, PV start %d, PV final %d\n", name, n, nnz(Ybus), maximum(diff(Yt.colptr)), ci.iterations, count(==(:PV), ci.types_start), count(==(:PV), ci.types_final))
  case_ok = true
  single = Dict{String,Any}()
  for (label, V, types, S, dP, dQ) in states
    g = EXT._device_jacobian_mismatch(Ybus, slack, reshape(V, n, 1), [types], reshape(S, n, 1), reshape(ci.Vset, n, 1), reshape(dP, n, 1), reshape(dQ, n, 1); T = Float64)
    g32 = EXT._device_jacobian_mismatch(Ybus, slack, reshape(V, n, 1), [types], reshape(S, n, 1), reshape(ci.Vset, n, 1), reshape(dP, n, 1), reshape(dQ, n, 1); T = Float32)
    single[label] = g
    # CPU references
    # Sparlectra's CPU SpMV, the one mismatch_rectangular and calc_injections use
    Icpu = Sparlectra.calc_currents(Ybus, V)
    I_unf = spmv_unfused(Yt, V)
    I_fus = spmv_fused(Yt, V)
    n_unf = bitwise_count(Icpu, I_unf)
    n_fus = bitwise_count(Icpu, I_fus)
    spmv_variant = n_unf == n ? "unfused" : n_fus == n ? "fused" : "neither"
    Jcpu = Sparlectra.build_rectangular_jacobian_pq_pv_sparse(Ybus, V, types, ci.Vset, slack; dPinj_dVm = dP, dQinj_dVm = dQ, structural_pattern = true)
    # the in-place path of C1 (recorded positions, fill!(nz, 0.0), nz[pos] += val)
    asm = Sparlectra.RectangularJacobianAssembly()
    Sparlectra.build_rectangular_jacobian_pq_pv_sparse(Ybus, V, types, ci.Vset, slack; dPinj_dVm = dP, dQinj_dVm = dQ, structural_pattern = true, assembly = asm)
    Jinp = copy(Sparlectra.build_rectangular_jacobian_pq_pv_sparse(Ybus, V, types, ci.Vset, slack; dPinj_dVm = dP, dQinj_dVm = dQ, structural_pattern = true, assembly = asm))
    Jallpq = Sparlectra.build_rectangular_jacobian_pq_pv_sparse(Ybus, V, allpq, ci.Vset, slack; structural_pattern = true)
    Fcpu = Sparlectra.mismatch_rectangular(Ybus, V, S, types, ci.Vset, slack)
    R, M = row_magnitudes(Yt, V)
    # pattern: GPU pattern == CPU structural pattern with every bus PQ
    Gt = SparseMatrixCSC(2(n - 1), 2(n - 1), Int.(g.jrowptr), Int.(g.jcol), ones(length(g.jcol)))
    Gp = copy(transpose(Gt))
    pattern_ok = Gp.colptr == Jallpq.colptr && Gp.rowval == Jallpq.rowval
    gv = vec(g.J[:, 1])
    cj, missing_cpu, gpu_only_nonzero = compare_jacobian(g.jrowptr, g.jcol, gv, Jcpu, nonslack, n, Ybus, V, R, M, types, dP, dQ, 1e-9)
    cji, missing_inp, _ = compare_jacobian(g.jrowptr, g.jcol, gv, Jinp, nonslack, n, Ybus, V, R, M, types, dP, dQ, 1e-9)
    # the two CPU paths differ at most in the sign of an exact zero
    signzero_cpu = count(k -> Jcpu.nzval[k] !== Jinp.nzval[k], eachindex(Jcpu.nzval))
    value_cpu = count(k -> Jcpu.nzval[k] != Jinp.nzval[k], eachindex(Jcpu.nzval))
    cf = compare_mismatch(vec(g.F[:, 1]), Fcpu, nonslack, V, S, ci.Vset, types, R, M)
    I_bit = bitwise_count(vec(g.I[:, 1]), Icpu)
    vm_bit = bitwise_count(vec(g.vm[:, 1]), abs.(V))
    norm_cpu = maximum(abs.(Fcpu))
    norm_ok = g.norm[1] === norm_cpu
    # the Q-limit sample kernel: Q = imag(S_calc) and abs(V) of every bus
    # against calc_injections and abs, bitwise, and the norm it carries
    req = Int32[collect(Int32(1):Int32(n)); -collect(Int32(1):Int32(n))]
    qs = EXT._device_qlimit_sample(Ybus, slack, V, types, S, ci.Vset, req)
    Scalc = Sparlectra.calc_injections(Ybus, V)
    q_bit = count(i -> qs.values[i] === imag(Scalc[i]), 1:n)
    qvm_bit = count(i -> qs.values[n+i] === abs(V[i]), 1:n)
    qs_ok = q_bit == n && qvm_bit == n && qs.norm === norm_cpu
    # Float32 instance: the kernel converts its Float64 value, so the G4
    # values equal Float32 of the G3 values exactly
    g32v = vec(g32.J[:, 1])
    f32_from_f64 = count(k -> g32v[k] === Float32(gv[k]), eachindex(gv))
    f32_rhs_ok = vec(g32.F[:, 1]) == vec(g.F[:, 1])
    passed = pattern_ok && missing_cpu == 0 && gpu_only_nonzero == 0 && cj.fails == 0 && cf.fails == 0 && spmv_variant == "unfused" && I_bit == n && qs_ok && vm_bit == n && norm_ok && f32_from_f64 == length(gv) && f32_rhs_ok && value_cpu == 0
    case_ok &= passed
    @printf("  %-26s J %8d entries, bitwise %8d, fails %d, max err/sigma %.2e, max rel %.2e, max abs(cpu=0) %.2e | F %7d, bitwise %7d, fails %d, max err/sigma %.2e, max rel %.2e | I %d/%d vm %d/%d bitwise, norm %s, SpMV %s (unfused %d, fused %d of %d), pattern %s, CPU-only %d, GPU-only nonzero %d, in-place: sign-of-zero diffs %d, bitwise vs GPU %d | f32 %d/%d | Q-limit sample Q %d/%d vm %d/%d norm %s | %s\n",
      label, cj.entries, cj.bitwise, cj.fails, cj.max_sigma_err, cj.max_rel, cj.max_abs_cpu0, cf.entries, cf.bitwise, cf.fails, cf.max_sigma_err, cf.max_rel, I_bit, n, vm_bit, n, norm_ok ? "equal" : "DIFF $(g.norm[1]) vs $(norm_cpu)", spmv_variant, n_unf, n_fus, n, pattern_ok, missing_cpu, gpu_only_nonzero, signzero_cpu, cji.bitwise, f32_from_f64, length(gv), q_bit, n, qvm_bit, n, qs.norm === norm_cpu ? "equal" : "DIFF", passed ? "PASS" : "FAIL")
    isempty(cj.worst) || println("    J worst: ", cj.worst)
    isempty(cf.worst) || println("    F worst: ", cf.worst)
    push!(csv_rows, (case = name, state = label, n = n, nnzJ = length(gv), N = 1, J_entries = cj.entries, J_bitwise = cj.bitwise, J_fails = cj.fails, J_max_err_over_sigma = cj.max_sigma_err, J_max_rel = cj.max_rel, J_max_abs_cpu0 = cj.max_abs_cpu0, F_entries = cf.entries, F_bitwise = cf.bitwise, F_fails = cf.fails, F_max_err_over_sigma = cf.max_sigma_err, F_max_rel = cf.max_rel, I_bitwise = I_bit, vm_bitwise = vm_bit, norm_equal = norm_ok, spmv_variant = spmv_variant, pattern_ok = pattern_ok, cpu_only = missing_cpu, gpu_only_nonzero = gpu_only_nonzero, inplace_signzero_diffs = signzero_cpu, inplace_value_diffs = value_cpu, inplace_bitwise_vs_gpu = cji.bitwise, f32_equal = f32_from_f64, qlimit_sample_q_bitwise = q_bit, qlimit_sample_vm_bitwise = qvm_bit, passed = passed))
  end
  # batch: the first four states in one launch with N = 4; each column must
  # equal its single run bit for bit
  bl = [s[1] for s in states[1:4]]
  Vm = reduce(hcat, [s[2] for s in states[1:4]])
  Sm = reduce(hcat, [s[4] for s in states[1:4]])
  dPm = reduce(hcat, [s[5] for s in states[1:4]])
  dQm = reduce(hcat, [s[6] for s in states[1:4]])
  gb = EXT._device_jacobian_mismatch(Ybus, slack, Vm, [s[3] for s in states[1:4]], Sm, repeat(ci.Vset, 1, 4), dPm, dQm; T = Float64)
  batch_ok = true
  for (s, label) in enumerate(bl)
    g1 = single[label]
    eq = gb.J[:, s] == g1.J[:, 1] && all(k -> gb.J[k, s] === g1.J[k, 1], axes(gb.J, 1)) && all(k -> gb.F[k, s] === g1.F[k, 1], axes(gb.F, 1)) && all(k -> gb.I[k, s] === g1.I[k, 1], axes(gb.I, 1)) && gb.norm[s] === g1.norm[1]
    batch_ok &= eq
  end
  @printf("  %-26s N = 4, every column bitwise equal to its single run: %s\n", "batch", batch_ok ? "PASS" : "FAIL")
  push!(csv_rows, (case = name, state = "batch_N4", n = n, nnzJ = size(gb.J, 1), N = 4, J_entries = length(gb.J), J_bitwise = batch_ok ? length(gb.J) : -1, J_fails = batch_ok ? 0 : 1, J_max_err_over_sigma = 0.0, J_max_rel = 0.0, J_max_abs_cpu0 = 0.0, F_entries = length(gb.F), F_bitwise = batch_ok ? length(gb.F) : -1, F_fails = batch_ok ? 0 : 1, F_max_err_over_sigma = 0.0, F_max_rel = 0.0, I_bitwise = -1, vm_bitwise = -1, norm_equal = batch_ok, spmv_variant = "", pattern_ok = true, cpu_only = 0, gpu_only_nonzero = 0, inplace_signzero_diffs = 0, inplace_value_diffs = 0, inplace_bitwise_vs_gpu = 0, f32_equal = 0, qlimit_sample_q_bitwise = -1, qlimit_sample_vm_bitwise = -1, passed = batch_ok))
  case_ok &= batch_ok
  # non-finite injection: a NaN and an Inf in S of the last bus must reach
  # the norm as on the CPU (maximum(abs.(F)) is NaN and Inf)
  k = nonslack[end]
  Snan = copy(ci.S)
  Snan[k] = ComplexF64(NaN, 0.0)
  Sinf = copy(ci.S)
  Sinf[k] = ComplexF64(Inf, 0.0)
  gn = EXT._device_jacobian_mismatch(Ybus, slack, repeat(ci.Vconv, 1, 2), [ci.types_final, ci.types_final], hcat(Snan, Sinf), repeat(ci.Vset, 1, 2), zeros(n, 2), zeros(n, 2); T = Float64)
  cn = (maximum(abs.(Sparlectra.mismatch_rectangular(Ybus, ci.Vconv, Snan, ci.types_final, ci.Vset, slack))), maximum(abs.(Sparlectra.mismatch_rectangular(Ybus, ci.Vconv, Sinf, ci.types_final, ci.Vset, slack))))
  nonfinite_ok = isnan(gn.norm[1]) && isnan(cn[1]) && gn.norm[2] == Inf && cn[2] == Inf
  @printf("  %-26s device norm %s and %s, CPU %s and %s: %s\n", "nonfinite_norm", gn.norm[1], gn.norm[2], cn[1], cn[2], nonfinite_ok ? "PASS" : "FAIL")
  case_ok &= nonfinite_ok
  # the update kernel against _apply_rectangular_delta! (not bitwise by
  # design: complex division and sincos on the device), reported only
  dx = 1e-3 .* (2 .* rand(rng, 2(n - 1)) .- 1)
  upd = Dict{Symbol,Float64}()
  for polar in (true, false)
    Vg = EXT._device_update(Ybus, slack, Vpert, dx, 1.0, polar)
    Vc = Sparlectra._apply_rectangular_delta(Vpert, dx, slack, nonslack, 1.0; update = polar ? :polar : :rectangular)
    upd[polar ? :polar : :rectangular] = maximum(abs.(Vg .- Vc))
    upd[polar ? :polar_bitwise : :rect_bitwise] = bitwise_count(Vg, Vc)
  end
  @printf("  %-26s polar max|dV| %.2e (bitwise %d/%d), rectangular max|dV| %.2e (bitwise %d/%d)\n", "update_vs_cpu", upd[:polar], Int(upd[:polar_bitwise]), n, upd[:rectangular], Int(upd[:rect_bitwise]), n)
  summary[name] = (ok = case_ok, update_polar = upd[:polar], max_row = maximum(diff(Yt.colptr)), start = t0, stop = now())
  return case_ok
end

function main()
  t_start = now()
  commit = strip(read(`git -C $ROOT rev-parse --short HEAD`, String))
  println("check_jacobian_gpu.jl  started ", Dates.format(t_start, "yyyy-mm-dd HH:MM:SS"))
  println("commit ", commit, ", Julia ", VERSION, ", CUDA.jl ", pkgversion(CUDA), ", CUDSS.jl ", pkgversion(CUDSS), ", device ", CUDA.name(CUDA.device()), " (", CUDA.capability(CUDA.device()), "), Sparlectra ", pkgversion(Sparlectra))
  println("CPU have_fma(Float64) = ", Core.Intrinsics.have_fma(Float64), ", acceptance |gpu - cpu| <= ", TOL, " * sigma")
  # which association Sparlectra's CPU SpMV gives here, on case118 with a
  # random iterate (every entry has several summands)
  spmv = spmv_association("case118")
  @printf("CPU SpMV (Sparlectra.calc_currents) on case118, random V: unfused bitwise %d/%d, fused bitwise %d/%d -> %s\n", spmv.unfused, spmv.n, spmv.fused, spmv.n, spmv.variant)
  if spmv.variant != "unfused"
    println("FAIL: the CPU SpMV is ", spmv.variant, ", the kernel _k_current_mismatch! computes the unfused association; the 1e-14 criterion would compare against a different association. Adapt the kernel to the CPU before running the proof.")
    return false
  end
  # static FMA check on the kernels compiled for case14
  ci14 = case_inputs("case14")
  fma_ok, rows = static_fma_check(ci14.Ybus, ci14.slack)
  println("\nStatic FMA check (PTX of the own kernels):")
  foreach(r -> println("  ", r), rows)
  # abs(V) at the edges of Base._hypot (scaling branches, zero, subnormal,
  # Inf, NaN), the branches the network states do not reach
  edge = ComplexF64[0.0, 1e-310 + 3e-310im, 1e-300 + 2e-300im, 1e300 + 1e300im, 1.7e308 + 1e308im, 3.0 + 4.0im, 1.0 + 1e-9im, complex(Inf, 1.0), complex(NaN, 1.0), complex(-2.5, NaN), -1e-200 + 0.0im, 0.6 - 0.8im, 5e-324 + 0.0im, 1.0 + 1.0im]
  n14 = ci14.n
  Ve = [edge; fill(1.0 + 0.0im, max(n14 - length(edge), 0))][1:n14]
  ge = EXT._device_jacobian_mismatch(ci14.Ybus, ci14.slack, reshape(Ve, n14, 1), [ci14.types_start], reshape(ci14.S, n14, 1), reshape(ci14.Vset, n14, 1), zeros(n14, 1), zeros(n14, 1))
  hypot_ok = all(k -> isequal(ge.vm[k, 1], abs(Ve[k])), 1:n14)
  println("abs(V) at the hypot edge cases (zero, subnormal, scaled, Inf, NaN): ", count(k -> isequal(ge.vm[k, 1], abs(Ve[k])), 1:n14), "/", n14, " bitwise, ", hypot_ok ? "PASS" : "FAIL")
  csv_rows = NamedTuple[]
  summary = Dict{String,Any}()
  all_ok = fma_ok && hypot_ok
  for name in CASES
    ok = try
      run_case(name, csv_rows, summary)
    catch e
      e isa CUDA.OutOfGPUMemoryError || rethrow(e)
      println("  ", name, ": out of GPU memory, skipped: ", sprint(showerror, e))
      false
    end
    all_ok &= ok
  end
  dir = out_dir()
  file = joinpath(dir, "check_jacobian_gpu.csv")
  open(file, "w") do io
    ks = keys(first(csv_rows))
    println(io, "# check_jacobian_gpu.jl commit ", commit, " started ", t_start, " finished ", now(), " static_fma_check ", fma_ok ? "pass" : "FAIL", " hypot_edges ", hypot_ok ? "pass" : "FAIL")
    println(io, join(ks, ','))
    for r in csv_rows
      println(io, join((string(getfield(r, k)) for k in ks), ','))
    end
  end
  t_stop = now()
  println("\nsummary:")
  for name in CASES
    haskey(summary, name) || (println("  ", name, ": not run"); continue)
    s = summary[name]
    @printf("  %-16s %s  max Ybus row %3d  update polar max|dV| %.2e\n", name, s.ok ? "PASS" : "FAIL", s.max_row, s.update_polar)
  end
  println("static FMA check: ", fma_ok ? "PASS" : "FAIL")
  println("CSV: ", file)
  println("finished ", Dates.format(t_stop, "yyyy-mm-dd HH:MM:SS"), " (", round((t_stop - t_start).value / 1000; digits = 1), " s), overall ", all_ok ? "PASS" : "FAIL")
  return all_ok
end

exit(Base.invokelatest(main) ? 0 : 1)
