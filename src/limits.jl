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

# file: src/limits.jl
# Helper function: ensures that the vector is filled to at least 'bus'
# and fills with a default value (via append!) up to the required length.
function _ensure_bus_index!(v::Vector{T}, bus::Int, default::T) where {T}
  if length(v) < bus
    append!(v, fill(default, bus - length(v)))
  end
  return v
end

"""
    printQLimitLog(net::Net; sort_by=:iter, io::IO=stdout)

Pretty-prints the structured Q-limit log (`net.qLimitLog`) as a small table.

# Keyword arguments
- `sort_by`:  `:iter` (default) or `:bus` — controls sorting.
- `io`: optional output stream (default = `stdout`).

Each line shows:
where `Side` is `:min` or `:max`.
"""
function printQLimitLog(net::Net; sort_by::Symbol = :iter, io::IO = stdout)
  logdata = net.qLimitLog
  if isempty(logdata)
    println(io, "No Q-limit events recorded.")
    return
  end

  # sort by chosen field
  if sort_by == :bus
    sort!(logdata, by = x -> x.bus)
  else
    sort!(logdata, by = x -> x.iter)
  end

  println(io, "─"^34)
  println(io, " Iteration │ Bus │ Side")
  println(io, "─"^34)

  for ev in logdata
    @printf(io, " %9d │ %3d │ %-4s\n", ev.iter, ev.bus, String(ev.side))
  end

  println(io, "─"^34)
  println(io, "Total events: ", length(logdata))
end

Base.@kwdef struct QLimitEvent
  iter::Int
  bus::Int
  side::Symbol   # :min | :max
end

function logQLimitHit!(net::Net, iter::Int, bus::Int, side::Symbol)
  push!(net.qLimitLog, QLimitEvent(iter = iter, bus = bus, side = side))
  net.qLimitEvents[bus] = side
end

"Returns the last iteration number where `bus` hit a Q-limit, or `nothing`."
function lastQLimitIter(net::Net, bus::Int)
  for i = length(net.qLimitLog):-1:1
    if net.qLimitLog[i].bus == bus
      return net.qLimitLog[i].iter
    end
  end
  return nothing
end

function resetQLimitLog!(net::Net)
  empty!(net.qLimitEvents)
  empty!(net.qLimitLog)
  return nothing
end

"""
    getQLimits_pu(net::Net) -> (qmin_pu, qmax_pu)

Return per-bus Q limits in p.u. (build once if empty).
"""
function getQLimits_pu(net::Net)
  if isempty(net.qmin_pu) || isempty(net.qmax_pu)
    buildQLimits!(net)
  end
  _sanitize_q_limits!(net.qmin_pu, net.qmax_pu, length(net.nodeVec))
  return net.qmin_pu, net.qmax_pu
end

"""
    printPVQLimitsTable(net::Net; io::IO=stdout, max_rows::Int=30)

Print a compact table of PV-bus reactive limits before the PF iteration starts.
Values are shown in **MVAr**.
"""
function printPVQLimitsTable(net::Net; io::IO = stdout, max_rows::Int = 30)
  qmin_pu, qmax_pu = getQLimits_pu(net)
  rows = Tuple{Int,Float64,Float64}[]

  for (bus, node) in enumerate(net.nodeVec)
    getNodeType(node) == PV || continue
    qmin = ((bus <= length(qmin_pu)) && isfinite(qmin_pu[bus])) ? (qmin_pu[bus] * net.baseMVA) : -Inf
    qmax = ((bus <= length(qmax_pu)) && isfinite(qmax_pu[bus])) ? (qmax_pu[bus] * net.baseMVA) : Inf
    push!(rows, (bus, qmin, qmax))
  end

  if isempty(rows)
    println(io, "PV Q-limits (MVAr): no PV buses.")
    return nothing
  end

  println(io, "PV Q-limits before PF run (MVAr):")
  println(io, "──────────────────────────────────────────────")
  println(io, " Bus │      Qmin [MVAr] │      Qmax [MVAr]")
  println(io, "──────────────────────────────────────────────")
  shown = max_rows < 0 ? length(rows) : min(max_rows, length(rows))
  for i in 1:shown
    bus, qmin, qmax = rows[i]
    @printf(io, " %3d │ %15.6f │ %15.6f\n", bus, qmin, qmax)
  end
  if shown < length(rows)
    println(io, " ...")
    @printf(io, " (%d more PV rows omitted; increase max_rows for full output)\n", length(rows) - shown)
  end
  println(io, "──────────────────────────────────────────────")
  return nothing
end
"""
    pv_hit_q_limit(net, pv_names)

Returns `true` if any of the PV buses from `pv_names`
appears in `net.qLimitEvents`.
`pv_names` is a list of bus names (strings).
"""
function pv_hit_q_limit(net, pv_names)
  pv_idx = map(name -> geNetBusIdx(net = net, busName = name), pv_names)
  return any(haskey(net.qLimitEvents, i) for i in pv_idx)
end

# ------------------------------
# Q-limit helpers (local arrays)
# ------------------------------
@inline function has_q_limits(qmin_pu::AbstractVector, qmax_pu::AbstractVector, i::Int)::Bool
  has_lo = (i <= length(qmin_pu)) && isfinite(qmin_pu[i])
  has_hi = (i <= length(qmax_pu)) && isfinite(qmax_pu[i])
  return has_lo || has_hi
end

@inline function q_limit_band(qmin_pu::AbstractVector, qmax_pu::AbstractVector, i::Int, q_hyst_pu::Float64)::Tuple{Float64,Float64}
  lo = ((i <= length(qmin_pu)) && isfinite(qmin_pu[i])) ? (qmin_pu[i] + q_hyst_pu) : -Inf
  hi = ((i <= length(qmax_pu)) && isfinite(qmax_pu[i])) ? (qmax_pu[i] - q_hyst_pu) : Inf
  return lo, hi
end

# ------------------------------
# Q-limit sanitizing (drop-in)
# ------------------------------
function _sanitize_q_limits!(qmin_pu::Vector{Float64}, qmax_pu::Vector{Float64}, nb::Int)
  _ensure_bus_index!(qmin_pu, nb, -Inf)
  _ensure_bus_index!(qmax_pu, nb, Inf)

  @inbounds for i = 1:nb
    # Treat NaN as "no limit"
    if isnan(qmin_pu[i])
      qmin_pu[i] = -Inf
    end
    if isnan(qmax_pu[i])
      qmax_pu[i] = Inf
    end
  end
  return qmin_pu, qmax_pu
end

"""
    validate_q_limit_signs!(qmin_pu, qmax_pu; io::IO=stdout, autocorrect::Bool=false, warn::Bool=true)

Validate Q-limit sign conventions per bus:
- expect `qmin ≤ 0`
- expect `qmax ≥ 0`
- expect `qmin ≤ qmax`

If `autocorrect=true`, suspicious sign-only limits are flipped and inverted ranges are swapped.
"""
function validate_q_limit_signs!(qmin_pu::AbstractVector{Float64}, qmax_pu::AbstractVector{Float64}; io::IO = stdout, autocorrect::Bool = false, warn::Bool = true)
  n = min(length(qmin_pu), length(qmax_pu))
  corrected = 0
  flagged = 0

  for bus in 1:n
    qmin = qmin_pu[bus]
    qmax = qmax_pu[bus]
    (isfinite(qmin) || isfinite(qmax)) || continue

    changed = false
    msgs = String[]

    if isfinite(qmin) && isfinite(qmax) && (qmin > qmax)
      push!(msgs, "qmin > qmax")
      if autocorrect
        qmin_pu[bus], qmax_pu[bus] = qmax, qmin
        qmin, qmax = qmin_pu[bus], qmax_pu[bus]
        changed = true
      end
    end

    if isfinite(qmin) && (qmin > 0.0)
      push!(msgs, "qmin positive")
      if autocorrect
        qmin_pu[bus] = -abs(qmin)
        qmin = qmin_pu[bus]
        changed = true
      end
    end

    if isfinite(qmax) && (qmax < 0.0)
      push!(msgs, "qmax negative")
      if autocorrect
        qmax_pu[bus] = abs(qmax)
        qmax = qmax_pu[bus]
        changed = true
      end
    end

    if !isempty(msgs)
      flagged += 1
      if warn
        action = changed ? "corrected" : "detected"
        @printf(io, "Warning: Q-limit sign issue at bus %d (%s) -> %s [qmin=%.6f, qmax=%.6f]\n", bus, join(msgs, ", "), action, qmin, qmax)
      end
      corrected += changed ? 1 : 0
    end
  end

  return (flagged = flagged, corrected = corrected)
end

function _qlimit_headroom(limit::Float64, headroom::Float64)
  if !isfinite(limit)
    return Inf
  end
  return abs(limit) * max(headroom, 0.0)
end

"""
    printFinalLimitValidation(net::Net; q_headroom::Float64=0.20, io::IO=stdout)

Prints post-PF validation tables for violated Q limits and voltage limits.
Returns `(q_violations, v_violations)`.
"""
function printFinalLimitValidation(net::Net; q_headroom::Float64 = 0.20, io::IO = stdout)
  qmin_pu, qmax_pu = getQLimits_pu(net)
  qrows = Tuple{Int,Float64,Float64,Float64,Float64}[]
  vrows = Tuple{Int,Float64,Float64,Float64}[]
  qgen_missing = 0

  for (bus, node) in enumerate(net.nodeVec)
    qgen_pu = isnothing(node._qƩGen) ? nothing : (Float64(node._qƩGen) / net.baseMVA)
    has_qgen = !isnothing(qgen_pu)
    if !has_qgen && (((bus <= length(qmin_pu)) && isfinite(qmin_pu[bus])) || ((bus <= length(qmax_pu)) && isfinite(qmax_pu[bus])))
      qgen_missing += 1
    end

    if has_qgen && (bus <= length(qmin_pu)) && isfinite(qmin_pu[bus])
      lo = qmin_pu[bus]
      if qgen_pu < (lo - _qlimit_headroom(lo, q_headroom))
        push!(qrows, (bus, qgen_pu * net.baseMVA, lo * net.baseMVA, qmax_pu[bus] * net.baseMVA, (lo - qgen_pu) * net.baseMVA))
      end
    end
    if has_qgen && (bus <= length(qmax_pu)) && isfinite(qmax_pu[bus])
      hi = qmax_pu[bus]
      if qgen_pu > (hi + _qlimit_headroom(hi, q_headroom))
        push!(qrows, (bus, qgen_pu * net.baseMVA, qmin_pu[bus] * net.baseMVA, hi * net.baseMVA, (qgen_pu - hi) * net.baseMVA))
      end
    end

    vm = node._vm_pu
    isnothing(vm) && continue
    vmin = isnothing(node._vmin_pu) ? net.vmin_pu : node._vmin_pu
    vmax = isnothing(node._vmax_pu) ? net.vmax_pu : node._vmax_pu
    vmf = Float64(vm)
    if vmf < vmin || vmf > vmax
      push!(vrows, (bus, vmf, vmin, vmax))
    end
  end

  println(io, "Final limit validation:")
  if isempty(qrows)
    @printf(io, "  Q-limits: no violations beyond %.1f%% headroom.\n", 100 * max(q_headroom, 0.0))
  else
    println(io, "  Q-limit violations (MVAr):")
    println(io, "  Bus │      Qgen │      Qmin │      Qmax │  Violation")
    for (bus, qgen, qmin, qmax, vio) in qrows
      @printf(io, " %4d │ %9.3f │ %9.3f │ %9.3f │ %10.3f\n", bus, qgen, qmin, qmax, vio)
    end
  end
  if qgen_missing > 0
    @printf(io, "  Q-limits: skipped %d bus(es) with finite Q-limits but missing _qƩGen.\n", qgen_missing)
  end

  if isempty(vrows)
    println(io, "  Voltage limits: no violations.")
  else
    println(io, "  Voltage limit violations (p.u.):")
    println(io, "  Bus │        Vm │      Vmin │      Vmax")
    for (bus, vm, vmin, vmax) in vrows
      @printf(io, " %4d │ %9.5f │ %9.5f │ %9.5f\n", bus, vm, vmin, vmax)
    end
  end

  return (q_violations = length(qrows), v_violations = length(vrows))
end

"""
    active_set_q_limits!(
        net, it, nb;
        get_qreq_pu,
        is_pv,
        make_pq!,
        make_pv!,
        qmin_pu,
        qmax_pu,
        pv_orig_mask,
        allow_reenable::Bool,
        q_hyst_pu::Float64,
        cooldown_iters::Int,
        verbose::Int=0,
    ) -> (changed::Bool, reenabled::Bool)

Core PV/Q-limit active-set logic shared by solvers.

Callbacks:
- get_qreq_pu(bus) -> Float64
- is_pv(bus) -> Bool
- make_pq!(bus, q_clamp_pu::Float64, side::Symbol)  # side = :min/:max
- make_pv!(bus)
"""
function active_set_q_limits!(
  net::Net,
  it::Int,
  nb::Int;
  get_qreq_pu,
  is_pv,
  make_pq!,
  make_pv!,
  qmin_pu::AbstractVector,
  qmax_pu::AbstractVector,
  pv_orig_mask,                      # AbstractVector{Bool} length nb
  allow_reenable::Bool,
  q_hyst_pu::Float64,
  cooldown_iters::Int,
  pv_to_pq_blocked::AbstractVector{Bool} = falses(nb),
  verbose::Int = 0,
)
  changed   = false
  reenabled = false

  # --- PV -> PQ ------------------------------------------------------------
  @inbounds for bus = 1:nb
    is_pv(bus) || continue
    pv_to_pq_blocked[bus] && continue
    qreq = get_qreq_pu(bus)

    has_q_limits(qmin_pu, qmax_pu, bus) || continue

    has_hi = (bus <= length(qmax_pu)) && isfinite(qmax_pu[bus])
    has_lo = (bus <= length(qmin_pu)) && isfinite(qmin_pu[bus])

    if has_hi && (qreq > qmax_pu[bus])
      make_pq!(bus, qmax_pu[bus], :max)
      logQLimitHit!(net, it, bus, :max)
      changed = true
      (verbose > 0) && @printf "PV->PQ Bus %d: Q=%.6f > Qmax=%.6f (it=%d)\n" bus qreq qmax_pu[bus] it

    elseif has_lo && (qreq < qmin_pu[bus])
      make_pq!(bus, qmin_pu[bus], :min)
      logQLimitHit!(net, it, bus, :min)
      changed = true
      (verbose > 0) && @printf "PV->PQ Bus %d: Q=%.6f < Qmin=%.6f (it=%d)\n" bus qreq qmin_pu[bus] it
    end
  end

  # --- Optional PQ -> PV ---------------------------------------------------
  if allow_reenable
    @inbounds for bus = 1:nb
      # Guards: was PV originally + has event + currently NOT PV
      pv_orig_mask[bus] || continue
      haskey(net.qLimitEvents, bus) || continue
      is_pv(bus) && continue

      has_q_limits(qmin_pu, qmax_pu, bus) || continue

      qreq = get_qreq_pu(bus)
      lo, hi = q_limit_band(qmin_pu, qmax_pu, bus, q_hyst_pu)

      ready = (qreq > lo) && (qreq < hi)

      if ready && (cooldown_iters > 0)
        last_it = lastQLimitIter(net, bus)
        if !isnothing(last_it) && (it - last_it) < cooldown_iters
          ready = false
        end
      end

      if ready
        make_pv!(bus)
        delete!(net.qLimitEvents, bus)   # clear event after re-enable
        reenabled = true
        (verbose > 0) && @printf "PQ->PV Bus %d: Q=%.6f within (%.6f, %.6f)\n" bus qreq lo hi
      end
    end
  end

  return changed, reenabled
end
