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
# file: benchmark/gpu/trace_host_copies_gpu.jl
# purpose: GPU study, stage 1b: an independent trace of the host copies of
#          the device Newton loop, per iteration. CUPTI (CUDA.@profile)
#          records every memcpy; a marker kernel launched at every fault
#          point of the loop attributes each copy to the point it follows.
#          Per run: the copies of each iteration with their sizes, the
#          composition of the Q-limit download, the loop counters for
#          comparison, the host API allocations and synchronizations per
#          segment (pre-flight, each hook point of each iteration), and
#          which records the end of the solve downloads. Second of two runs
#          (the first compiles).
#
# Run: julia --project=benchmark/gpu benchmark/gpu/trace_host_copies_gpu.jl [case ...]
#   default cases case118 case_ACTIVSg25k; Q limits on (default settings)
#   and off (ignore_q_limits), Float64 and float32_mixed.
#   The sample lines count the Q-limit requests; the download after a
#   switch iteration carries 8 B more (the norm after the switch).

using Sparlectra, CUDA, CUDSS, Printf, Dates
const EXT = Base.get_extension(Sparlectra, :SparlectraCUDAExt)
const QUIET = OutputConfig(logfile_results = :off, console_summary = false, startup_latency_hint = false)
cfg(cb; qlimits = QLimitConfig()) = SparlectraConfig(powerflow = PowerFlowConfig(compute_backend = cb, qlimits = qlimits), output = QUIET)
gpu_cb(prec) = Sparlectra.ComputeBackendConfig(mode = :gpu, gpu_precision = prec, gpu_fallback_to_cpu = true)
casefile(name) = Sparlectra.ensure_casefile(name * ".m"; to_jl = false)

# marker kernels: one launch per hook point, so every memcpy can be attributed
# by time and every host API call by the start time of the marker's launch
function mk_marker(x::Int32)
  return nothing
end
const HOOKLOG = Tuple{Symbol,Int}[]
function hook(point, it, st)
  push!(HOOKLOG, (point, it))
  @cuda threads = 1 mk_marker(Int32(it))
  return false
end

# which records the end of the solve downloads (monkeypatched copy of the
# original body with a log line): the last finite iterate (Vlast) unless the
# iteration that measured the device V is the last one, the best (Vbest)
# unless it is the last one
const RECLOG = NamedTuple[]
@eval EXT function _host_records!(L::_LoopLocals, st::DeviceNewtonState, Vh::Vector{ComplexF64}, vh_it::Int, d2h!::F) where {F}
  _fault(:host_records, L.iters, st)
  push!(Main.RECLOG, (last_it = L.last_it, best_it = L.best_it, vh_it = vh_it, vlast_download = L.last_it > 0 && L.last_it != vh_it, vbest_download = L.best_it > 0 && L.best_it != L.last_it, best_mis = L.best_mis))
  n = length(Vh)
  if L.last_it > 0
    L.vlast = L.last_it == vh_it ? copy(Vh) : d2h!(Vector{ComplexF64}(undef, n), st.Vlast, n)
  end
  if L.best_it > 0
    L.vbest = L.best_it == L.last_it ? copy(L.vlast) : d2h!(Vector{ComplexF64}(undef, n), st.Vbest, n)
  end
  return nothing
end

# host API calls reported per segment: allocations and synchronizations
const API_OF_INTEREST = r"alloc|free|ynchronize|StreamQuery|EventQuery"i

# request composition per download (logged in the unpack, monkeypatched copy of the original body)
const REQLOG = NamedTuple[]
@eval EXT function _qlimit_unpack!(link::_QLimitLink)
  q_pv = 0; q_other = 0; v_pv = 0; v_other = 0
  for r in link.req
    b = abs(Int(r))
    pv = link.codes[b] == _BT_PV
    if r > 0
      pv ? (q_pv += 1) : (q_other += 1)
    else
      pv ? (v_pv += 1) : (v_other += 1)
    end
  end
  npv = count(b -> link.codes[b] == _BT_PV, eachindex(link.codes))
  push!(Main.REQLOG, (n_pv = npv, q_pv = q_pv, q_other = q_other, v_pv = v_pv, v_other = v_other))
  q = link.sample.q
  vm = link.sample.vm
  fill!(q, NaN)
  fill!(vm, NaN)
  @inbounds for (k, r) in enumerate(link.req)
    r > 0 ? (q[r] = link.vals[1+k]) : (vm[-r] = link.vals[1+k])
  end
  return nothing
end

function one_run(p, prec, qlimits)
  empty!(HOOKLOG); empty!(REQLOG); empty!(RECLOG)
  EXT._FAULT_HOOK[] = hook
  local r
  prof = CUDA.@profile trace = true begin
    r = Sparlectra._with_device_newton(() -> run_sparlectra(casefile = p, config = cfg(gpu_cb(prec); qlimits = qlimits)))
  end
  EXT._FAULT_HOOK[] = nothing
  return r, prof
end

function analyse(name, prec, qlabel, qlimits)
  p = casefile(name)
  one_run(p, prec, qlimits)        # warm-up (compilation)
  r, prof = one_run(p, prec, qlimits)
  st = Sparlectra.rectangular_pf_status(r.net)
  nhist = length(st.mismatch_history)
  n = length(r.net.nodeVec)
  d = prof.device
  # events in time order
  idx = sortperm(d.start)
  seg = (:none, 0)
  hk = 0
  rows = NamedTuple[]
  # host start time of each marker launch (correlation id of the kernel record)
  h = prof.host
  hstart = Dict(zip(h.id, h.start))
  mk_host = Float64[]
  for i in idx
    nm = d.name[i]
    if occursin("mk_marker", nm)
      hk += 1
      seg = HOOKLOG[hk]
      push!(mk_host, hstart[d.id[i]])
      continue
    end
    kind = occursin("device to pageable", nm) || occursin("device to pinned", nm) ? :d2h :
           occursin("pageable to device", nm) || occursin("pinned to device", nm) ? :h2d :
           occursin("device to device", nm) ? :d2d : occursin("[set", nm) || occursin("memset", lowercase(nm)) ? :memset : :kernel
    kind in (:d2h, :h2d) || continue
    push!(rows, (seg = seg, kind = kind, bytes = coalesce(d.size[i], -1)))
  end
  @assert hk == length(HOOKLOG) "marker count $(hk) vs hooks $(length(HOOKLOG))"
  npre = count(s -> s[1] === :preflight, HOOKLOG)
  its = maximum((s[2] for s in HOOKLOG if s[1] === :mismatch); init = 0)
  d2h = filter(x -> x.kind === :d2h, rows)
  h2d = filter(x -> x.kind === :h2d, rows)
  inloop(x) = x.seg[1] ∉ (:none, :preflight, :preflight_kernels)
  d2h_loop = filter(inloop, d2h); h2d_loop = filter(inloop, h2d)
  vbytes = 16n
  vdl = filter(x -> x.bytes == vbytes, d2h_loop)
  small = filter(x -> x.bytes != vbytes, d2h_loop)
  @printf("\n%s %s Q limits %s: n %d, device solves %d, iterations %d (history %d), conv %s, gpu_newton %s, fallback %d\n", name, prec, qlabel, n, npre, its, nhist, r.final_converged, get(st, :gpu_newton, false), st.linear_solver_fallback_count)
  @printf("  CUPTI d2h: preflight %d, loop %d (%d B); of these V-size (%d B) %d; h2d preflight %d, loop %d (%d B)\n", length(d2h) - length(d2h_loop), length(d2h_loop), sum(x -> x.bytes, d2h_loop; init = 0), vbytes, length(vdl), length(h2d) - length(h2d_loop), length(h2d_loop), sum(x -> x.bytes, h2d_loop; init = 0))
  @printf("  loop counters: d2h_count %d d2h_bytes %d h2d_count %d h2d_bytes %d\n", get(st, :gpu_newton_d2h_count, -1), get(st, :gpu_newton_d2h_bytes, -1), get(st, :gpu_newton_h2d_count, -1), get(st, :gpu_newton_h2d_bytes, -1))
  # per iteration: copies after the mismatch marker of iteration k until the next marker of k+1
  for k = 1:its
    dk = filter(x -> x.seg[2] == k && x.seg[1] !== :none, d2h_loop)
    hk2 = filter(x -> x.seg[2] == k && x.seg[1] !== :none, h2d_loop)
    @printf("  it %2d d2h %s | h2d %s\n", k, join([string(x.seg[1], ":", x.bytes) for x in dk], " "), join([string(x.seg[1], ":", x.bytes) for x in hk2], " "))
  end
  for (k, q) in enumerate(REQLOG)
    @printf("  sample %2d: PV buses %d, Q of PV %d, Q of non-PV %d, |V| of PV %d, |V| of non-PV %d -> %d B\n", k, q.n_pv, q.q_pv, q.q_other, q.v_pv, q.v_other, 8 * (1 + q.q_pv + q.q_other + q.v_pv + q.v_other))
  end
  # host API side: memcpy calls and who issued them is not separable; report names
  names = Dict{String,Int}()
  for nm in h.name
    occursin("emcpy", nm) && (names[nm] = get(names, nm, 0) + 1)
  end
  println("  host API memcpy calls (whole run incl. preflight): ", names)
  # host API calls by segment: a call belongs to the last marker launched
  # before it started on the host (segment :none before the first marker)
  @assert issorted(mk_host) "marker launches out of order"
  hseg(t) = (k = searchsortedlast(mk_host, t); k == 0 ? (:none, 0) : HOOKLOG[k])
  apiseg = Dict{Tuple{Symbol,Int},Dict{String,Int}}()
  malloc_us = Float64[]
  for i in eachindex(h.name)
    nm = h.name[i]
    occursin(API_OF_INTEREST, nm) || continue
    sg = hseg(h.start[i])
    dd = get!(Dict{String,Int}, apiseg, sg)
    dd[nm] = get(dd, nm, 0) + 1
    occursin("cudaMalloc", nm) && sg[1] in (:step, :update) && push!(malloc_us, (h.stop[i] - h.start[i]) * 1e6)
  end
  fmt(dd) = isempty(dd) ? "-" : join(sort([string(k, " ", v) for (k, v) in dd]), ", ")
  # the pre-flight: from the start of the run to the first loop marker; its
  # cuDSS analysis and first factorization follow the :preflight_kernels marker
  pre = Dict{String,Int}()
  for sg in ((:none, 0), (:preflight, 0), (:preflight_kernels, 0)), (k, v) in get(apiseg, sg, Dict{String,Int}())
    pre[k] = get(pre, k, 0) + v
  end
  pre_d2h = count(x -> x.seg[1] === :preflight_kernels, d2h)
  println("  host API, pre-flight (once per solve): d2h copies after the pre-flight kernels (cuDSS analysis and first factorization) ", pre_d2h, "; ", fmt(pre))
  for k = 1:its
    parts = String[]
    for pt in (:vdep, :mismatch, :records, :switch, :step, :update, :check, :final, :host_records, :outcome)
      dd = get(apiseg, (pt, k), nothing)
      dd === nothing || isempty(dd) || push!(parts, string(pt, ": ", fmt(dd)))
    end
    @printf("  it %2d host API %s\n", k, isempty(parts) ? "-" : join(parts, " | "))
  end
  isempty(malloc_us) || @printf("  cudaMalloc in the Newton steps: %d calls, host duration median %.1f us, max %.1f us\n", length(malloc_us), sort(malloc_us)[cld(length(malloc_us), 2)], maximum(malloc_us))
  # the end of the solve: which records were downloaded
  for rl in RECLOG
    @printf("  end-of-solve records: measured iterate %d, last finite %d, best %d (norm %.3e); downloads: Vlast %s, Vbest %s\n", rl.vh_it, rl.last_it, rl.best_it, rl.best_mis, rl.vlast_download, rl.vbest_download)
  end
  @printf("  history: %s\n", join([@sprintf("%.3e", x) for x in st.mismatch_history], " "))
  println("  q-limit log entries: ", length(r.net.qLimitLog))
end

function main()
  cases = isempty(ARGS) ? ["case118", "case_ACTIVSg25k"] : ARGS
  commit = readchomp(`git -C $(@__DIR__) rev-parse --short HEAD`)
  println("trace_host_copies_gpu  started ", Dates.now(), "  commit ", commit, "  CUDA.jl ", pkgversion(CUDA), "  CUDSS.jl ", pkgversion(CUDSS), "  device ", CUDA.name(CUDA.device()))
  for name in cases, prec in (:float64, :float32_mixed)
    for (ql, qc) in (("on (default)", QLimitConfig()), ("off", QLimitConfig(ignore_q_limits = true)))
      analyse(name, prec, ql, qc)
    end
  end
end

Base.invokelatest(main)
