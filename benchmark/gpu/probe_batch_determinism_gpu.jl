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
# file: benchmark/gpu/probe_batch_determinism_gpu.jl
# purpose: GPU study, before stage 2: does the cuDSS deterministic mode hold
#          in the uniform batch (ubatch_size > 1), with and without
#          ubatch_mask? N Jacobians of one pattern, built by the device
#          kernels of stage 1b at N perturbed states of the converged case,
#          solved as one batch REPEATS times with a fresh solver each
#          (analysis, factorization, solve) and once more by refactorization
#          on the last solver. Compared bit for bit per column: the repeats
#          among themselves, the refactorization against the factorization,
#          and each batch column against a single solve of that Jacobian in
#          the same mode. The fast mode (deterministic_mode off) runs the
#          same way for contrast. Float64 (G3) and Float32 (G4).
#          cuDSS 1.0.0 refuses deterministic_mode with ubatch_size > 1
#          (CUDSS_STATUS_NOT_SUPPORTED in the analysis, for every reordering,
#          matching and factorization algorithm tried); the probe records
#          that and runs the non-uniform batch (CudssBatchedSolver, one
#          matrix object per scenario) as the deterministic alternative.
#          Every batch result carries its residual against its own Jacobian:
#          the batch modes run without matching (matching_alg other than
#          none is refused there), and the stage 1b Jacobian needs it (the
#          single solve without matching perturbs pivots as well). So the
#          probe runs a second time on the Jacobians with their columns
#          permuted by the matching cuDSS computes for the single solve of
#          scenario 1, every solver without matching.
#
# Run: julia --project=benchmark/gpu benchmark/gpu/probe_batch_determinism_gpu.jl [case ...]
#   default cases case2869pegase case_ACTIVSg25k; N (default 8), REPEATS
#   (default 10) and SHORT (1: fewer variants on the as-built Jacobians)
#   from the environment.
#   Exit code 0 only if every deterministic-mode comparison that cuDSS
#   accepts is bitwise.

using Sparlectra, CUDA, CUDSS, SparseArrays, LinearAlgebra, Random, Printf, Dates
using CUDA.CUSPARSE: CuSparseMatrixCSR

const EXT = Base.get_extension(Sparlectra, :SparlectraCUDAExt)
EXT === nothing && error("SparlectraCUDAExt not loaded")

const QUIET = OutputConfig(logfile_results = :off, console_summary = false, startup_latency_hint = false)
const NB = parse(Int, get(ENV, "N", "8"))
const REPEATS = parse(Int, get(ENV, "REPEATS", "10"))
const SHORT = get(ENV, "SHORT", "0") == "1"
# a batch result is wrong if its relative residual is above this and above
# 10 times the residual of the single solves of the same Jacobians (these
# reach about 1e-14 in Float64 and 1e-6 in Float32, but 4e-4 in Float32 on
# case_ACTIVSg25k)
const RES_OK = 1e-4
batch_solved(maxres, maxres_single) = isfinite(maxres) && maxres <= max(RES_OK, 10 * maxres_single)

bus_type_symbol(t) = t == Sparlectra.Slack ? :Slack : t == Sparlectra.PV ? :PV : :PQ

# N Jacobians and right-hand sides of one pattern from the device kernels:
# column 1 the converged state, columns 2..N the converged state with a
# perturbation of 1e-3 pu per bus (own seed per column)
function batch_inputs(name, T)
  p = Sparlectra.ensure_casefile(name * ".m"; to_jl = false)
  r = run_sparlectra(casefile = p, config = SparlectraConfig(output = QUIET))
  net = r.net
  n = length(net.nodeVec)
  Ybus = Sparlectra.createYBUS(net = net, sparse = true)
  _, slack = Sparlectra.initialVrect(net; flatstart = true)
  S = Sparlectra.buildComplexSVec(net)
  Vset = Sparlectra._bus_voltage_setpoints_from_prosumers(net)
  Vconv = Sparlectra.buildVoltageVector(net)
  types = [bus_type_symbol(Sparlectra.getNodeType(nd)) for nd in net.nodeVec]
  Vm = repeat(Vconv, 1, NB)
  for s = 2:NB
    rng = MersenneTwister(1000 + s)
    Vm[:, s] .+= ComplexF64.(1e-3 .* (2 .* rand(rng, n) .- 1), 1e-3 .* (2 .* rand(rng, n) .- 1))
    Vm[slack, s] = Vconv[slack]
  end
  g = EXT._device_jacobian_mismatch(Ybus, slack, Vm, fill(types, NB), repeat(S, 1, NB), repeat(Vset, 1, NB), zeros(n, NB), zeros(n, NB); T = T)
  nvar = size(g.F, 1)
  # the right-hand side as the loop has it: -F in the factorization type
  return (n = n, nvar = nvar, nnzJ = size(g.J, 1), rowptr = Int32.(g.jrowptr), colval = Int32.(g.jcol), J = g.J, rhs = T.(-g.F), converged = r.final_converged, matching = true)
end

# one batch solve with a fresh solver; with `refactor` also a refactorization
# and solve on the same solver. Returns the solutions as host matrices.
function batch_solve(bi, T; det::Bool, mask::Union{Nothing,Vector{Cint}} = nothing, refactor::Bool = false, sentinel = T(-7), reordering::String = "default")
  nz = CuVector{T}(vec(bi.J))
  solver = CudssSolver(CuVector{Int32}(bi.rowptr), CuVector{Int32}(bi.colval), nz, "G", 'F')
  cudss_set(solver, "ubatch_size", NB)
  cudss_set(solver, "reordering_alg", reordering)
  bi.matching || cudss_set(solver, "matching_alg", "none")
  det && cudss_set(solver, "deterministic_mode", 1)
  mask === nothing || cudss_set(solver, "ubatch_mask", mask)
  xg = CuVector{T}(fill(sentinel, bi.nvar * NB))
  bg = CuVector{T}(vec(bi.rhs))
  X = CudssMatrix(T, bi.nvar; nbatch = NB)
  cudss_update(X, xg)
  B = CudssMatrix(T, bi.nvar; nbatch = NB)
  cudss_update(B, bg)
  cudss("analysis", solver, X, B)
  cudss("factorization", solver, X, B)
  cudss("solve", solver, X, B; asynchronous = false)
  x1 = reshape(Array(xg), bi.nvar, NB)
  x2 = nothing
  if refactor
    fill!(xg, sentinel)
    cudss_update(solver, CuVector{Int32}(bi.rowptr), CuVector{Int32}(bi.colval), nz)
    cudss("refactorization", solver, X, B)
    cudss("solve", solver, X, B; asynchronous = false)
    x2 = reshape(Array(xg), bi.nvar, NB)
  end
  CUDA.synchronize()
  return x1, x2, Int(cudss_get(solver, "npivots"))
end

# the non-uniform batch: one CSR matrix per scenario in one batched solver
function nonuniform_solve(bi, T; det::Bool)
  A = [CuSparseMatrixCSR{T,Int32}(CuVector{Int32}(bi.rowptr), CuVector{Int32}(bi.colval), CuVector{T}(bi.J[:, s]), (bi.nvar, bi.nvar)) for s = 1:NB]
  x = [CUDA.zeros(T, bi.nvar) for _ = 1:NB]
  b = [CuVector{T}(bi.rhs[:, s]) for s = 1:NB]
  solver = CudssBatchedSolver(A, "G", 'F')
  det && cudss_set(solver, "deterministic_mode", 1)
  bi.matching || cudss_set(solver, "matching_alg", "none")
  cudss("analysis", solver, x, b)
  cudss("factorization", solver, x, b)
  cudss("solve", solver, x, b)
  CUDA.synchronize()
  return reduce(hcat, Array.(x))
end

# cuDSS rejects the combination (the deterministic uniform batch in cuDSS
# 1.0.0), or runs out of device memory in the analysis (the AMD reordering
# of the uniform batch on case_ACTIVSg25k without matching asks for more
# than the card has)
batch_refused(e) = (e isa CUDSS.CUDSSError && occursin("NOT_SUPPORTED", sprint(showerror, e))) || e isa CUDA.OutOfGPUMemoryError

# column s alone, a single-matrix solver in the same mode
function single_solve(bi, T, s; det::Bool)
  nz = CuVector{T}(bi.J[:, s])
  solver = CudssSolver(CuVector{Int32}(bi.rowptr), CuVector{Int32}(bi.colval), nz, "G", 'F')
  det && cudss_set(solver, "deterministic_mode", 1)
  bi.matching || cudss_set(solver, "matching_alg", "none")
  xg = CUDA.zeros(T, bi.nvar)
  bg = CuVector{T}(bi.rhs[:, s])
  cudss("analysis", solver, xg, bg)
  cudss("factorization", solver, xg, bg)
  cudss("solve", solver, xg, bg)
  CUDA.synchronize()
  return Array(xg)
end

bits(x) = reinterpret(UInt8, collect(x))
reldiff(a, b) = maximum(abs.(Float64.(a) .- Float64.(b))) / max(maximum(abs.(Float64.(b))), floatmin(Float64))

# relative residual of the solution x of scenario s in Float64, against the
# Jacobian of that scenario (the CSR values read as CSC give the transpose)
function residual(bi, x::AbstractVector, s)
  A = SparseMatrixCSC(transpose(SparseMatrixCSC(bi.nvar, bi.nvar, Int.(bi.rowptr), Int.(bi.colval), Float64.(bi.J[:, s]))))
  return norm(A * Float64.(x) - Float64.(bi.rhs[:, s]), Inf) / max(norm(Float64.(bi.rhs[:, s]), Inf), floatmin(Float64))
end

const FAILS = String[]

# The batch variants on one set of Jacobians: the uniform batch (with and
# without mask, several reorderings) and the non-uniform batch, each in the
# deterministic and the fast mode where cuDSS accepts it
function run_variants(bi, T, name, label, variants)
  mask = Cint[isodd(s) ? 1 : 0 for s = 1:NB]
  for (det, use_mask, reo) in variants
    msk = use_mask ? mask : nothing
    runs = Matrix{T}[]
    xr = nothing
    refused = ""
    npiv = -1
    for k = 1:REPEATS
      x1, x2, npiv = try
        batch_solve(bi, T; det = det, mask = msk, refactor = k == REPEATS, reordering = reo)
      catch e
        batch_refused(e) || rethrow(e)
        refused = sprint(showerror, e)
        break
      end
      push!(runs, x1)
      k == REPEATS && (xr = x2)
    end
    if !isempty(refused)
      @printf("  %-13s %-7s reordering %-7s uniform batch failed in cuDSS: %s\n", det ? "deterministic" : "fast", msk === nothing ? "no mask" : "mask", reo, first(refused, 90))
      continue
    end
    cols = msk === nothing ? collect(1:NB) : findall(==(1), msk)
    skipped = msk === nothing ? Int[] : findall(==(0), msk)
    # distinct bit patterns per processed column over the repeats
    distinct = [length(unique([bits(x[:, s]) for x in runs])) for s in cols]
    refac_same = all(s -> bits(xr[:, s]) == bits(runs[end][:, s]), cols)
    sentinel_ok = all(s -> all(==(T(-7)), runs[end][:, s]) && all(==(T(-7)), xr[:, s]), skipped)
    single = [single_solve(bi, T, s; det = det) for s in cols]
    vs_single = [bits(runs[1][:, s]) == bits(single[k]) for (k, s) in enumerate(cols)]
    maxrel_single = maximum(reldiff(runs[1][:, s], single[k]) for (k, s) in enumerate(cols))
    maxres = maximum(residual(bi, runs[1][:, s], s) for s in cols)
    maxres_single = maximum(residual(bi, single[k], s) for (k, s) in enumerate(cols))
    # equal bits of wrong (or non-finite) solutions say nothing about
    # determinism: judged only where the batch solves the systems
    solved = batch_solved(maxres, maxres_single)
    ok = all(==(1), distinct) && refac_same && sentinel_ok
    (det && solved) && check(ok, "$name $T $label deterministic $(msk === nothing ? "no mask" : "mask")")
    @printf("  %-13s %-7s reordering %-7s distinct results per column %s | refactorization bitwise %s | masked columns untouched %s | batch column = single solve bitwise %d of %d, max rel diff %.2e | max residual batch %.1e (npivots %d), single %.1e  %s\n",
      det ? "deterministic" : "fast", msk === nothing ? "no mask" : "mask", reo, string(distinct), refac_same, isempty(skipped) ? "-" : string(sentinel_ok), count(vs_single), length(cols), maxrel_single, maxres, npiv, maxres_single, !solved ? "WRONG SOLUTIONS, repeatability not judged" : det ? (ok ? "ok" : "FAIL") : "(reported)")
    flush(stdout)
  end
  # the non-uniform batch, deterministic and fast
  for det in (true, false)
    runs = try
      [nonuniform_solve(bi, T; det = det) for _ = 1:REPEATS]
    catch e
      batch_refused(e) || rethrow(e)
      @printf("  %-13s non-uniform batch failed in cuDSS: %s\n", det ? "deterministic" : "fast", first(sprint(showerror, e), 90))
      continue
    end
    distinct = [length(unique([bits(x[:, s]) for x in runs])) for s = 1:NB]
    single = [single_solve(bi, T, s; det = det) for s = 1:NB]
    vs_single = [bits(runs[1][:, s]) == bits(single[s]) for s = 1:NB]
    maxres = maximum(residual(bi, runs[1][:, s], s) for s = 1:NB)
    maxres_single = maximum(residual(bi, single[s], s) for s = 1:NB)
    solved = batch_solved(maxres, maxres_single)
    ok = all(==(1), distinct)
    (det && solved) && check(ok, "$name $T $label deterministic non-uniform batch")
    @printf("  %-13s non-uniform batch: distinct results per column %s | batch column = single solve bitwise %d of %d, max rel diff %.2e | max residual %.1e, single %.1e  %s\n",
      det ? "deterministic" : "fast", string(distinct), count(vs_single), NB, maximum(reldiff(runs[1][:, s], single[s]) for s = 1:NB), maxres, maxres_single, !solved ? "WRONG SOLUTIONS, repeatability not judged" : det ? (ok ? "ok" : "FAIL") : "(reported)")
    flush(stdout)
  end
end

# The Jacobians with their columns permuted by the matching cuDSS computes
# for the single solve of scenario 1 (one permutation for all scenarios,
# the pattern is shared). The permuted system needs no matching: its
# solvers run with matching_alg none, batch and single alike.
function matched(bi, T)
  n = bi.nvar
  d = CuSparseMatrixCSR{T,Int32}(CuVector{Int32}(bi.rowptr), CuVector{Int32}(bi.colval), CuVector{T}(bi.J[:, 1]), (n, n))
  s = CudssSolver(d, "G", 'F')
  x = CUDA.zeros(T, n)
  b = CuVector{T}(bi.rhs[:, 1])
  cudss("analysis", s, x, b)
  pm = CUDA.zeros(Int32, n)
  cudss_set(s, "perm_matching", pm)
  cudss_get(s, "perm_matching")
  q = Int.(Array(pm))
  isperm(q) || error("perm_matching is not a permutation")
  # CSR of A[:, q] per scenario: as CSC the CSR arrays describe the
  # transpose, whose rows are permuted by q
  cols = map(1:NB) do sc
    At = SparseMatrixCSC(n, n, Int.(bi.rowptr), Int.(bi.colval), bi.J[:, sc])
    copy(transpose(copy(transpose(At))[:, q]))
  end
  rowptr = Int32.(cols[1].colptr)
  colval = Int32.(cols[1].rowval)
  all(c -> c.colptr == cols[1].colptr && c.rowval == cols[1].rowval, cols) || error("the permuted patterns differ")
  J = reduce(hcat, [c.nzval for c in cols])
  zero_diag_before = count(i -> SparseMatrixCSC(n, n, Int.(bi.rowptr), Int.(bi.colval), bi.J[:, 1])[i, i] == 0, 1:n)
  zero_diag_after = count(i -> copy(transpose(cols[1]))[i, i] == 0, 1:n)
  return merge(bi, (rowptr = rowptr, colval = colval, J = J, matching = false, zero_diag = (zero_diag_before, zero_diag_after)))
end

function probe(name)
  for T in (Float64, Float32)
    bi = batch_inputs(name, T)
    @printf("\n%s %s: n %d, nvar %d, nnzJ %d, N %d, repeats %d, C1 converged %s\n", name, T, bi.n, bi.nvar, bi.nnzJ, NB, REPEATS, bi.converged)
    println(" Jacobians as the device map builds them (single solves with the default matching):")
    # deterministic and fast with the default reordering, and the fast mode
    # with the COLAMD-family and AMD reorderings (wrong uniform batch
    # results with the default were reported for cuDSS 0.8.0)
    asbuilt = [(true, false, "default"), (true, true, "default"), (false, false, "default"), (false, true, "default"), (false, false, "algo1"), (false, false, "algo2"), (false, false, "algo3")]
    # SHORT=1: the as-built Jacobians only in the deterministic mode and the
    # fast default (the wrong batch solutions take long on large cases)
    SHORT && (asbuilt = asbuilt[1:3])
    run_variants(bi, T, name, "as built", asbuilt)
    bm = matched(bi, T)
    @printf(" Jacobians with the columns permuted by the cuDSS matching of scenario 1 (zero diagonal entries %d before, %d after), every solver without matching:\n", bm.zero_diag...)
    run_variants(bm, T, name, "matched", [(true, false, "default"), (true, true, "default"), (false, false, "default"), (false, true, "default")])
  end
end

check(ok, what) = (ok || push!(FAILS, what); ok)

function main()
  cases = isempty(ARGS) ? ["case2869pegase", "case_ACTIVSg25k"] : ARGS
  commit = try
    readchomp(`git -C $(@__DIR__) rev-parse --short HEAD`)
  catch e
    e isa ProcessFailedException || e isa Base.IOError || rethrow(e)
    "unknown"
  end
  println("probe_batch_determinism_gpu  started ", now(), "  commit ", commit, "  Julia ", VERSION, "  CUDA.jl ", pkgversion(CUDA), "  CUDSS.jl ", pkgversion(CUDSS), "  device ", CUDA.name(CUDA.device()))
  for name in cases
    probe(name)
  end
  println("\nfinished ", now(), "  failed checks: ", length(FAILS))
  foreach(f -> println("  FAIL ", f), FAILS)
  exit(isempty(FAILS) ? 0 : 1)
end

Base.invokelatest(main)
