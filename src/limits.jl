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

# Helper function: ensures that the vector is filled to at least 'bus'
# and fills with a default value (via append!) up to the required length.

# file: src/limits.jl
# purpose: reactive-limit bookkeeping for PV buses: the Q-limit event log and
#          its printers, PV limit snapshots, limit sanitation and validation,
#          and the active_set_q_limits! switching helper

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
- `sort_by`:  `:iter` (default) or `:bus`, controls sorting.
- `io`: optional output stream (default = `stdout`).

Each row shows the iteration, the bus, the side (`:min` or `:max`) and the
change: `PV->PQ` for a limit hit, `PV->PQ (group)` for the revert of a
jointly released group, `PQ->PV` for a release (the side is the limit the
bus was clamped at).
"""
function printQLimitLog(net::Net; sort_by::Symbol = :iter, io::IO = stdout, max_rows::Int = 200, full_details::Union{Nothing,AbstractString} = nothing)
  logdata = collect(net.qLimitLog)
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

  println(io, "Q-limit switching events:")
  println(io, "  total events     : ", length(logdata))
  println(io, "  PV->PQ           : ", count(ev -> _qlimit_event_kind(ev) !== :release, logdata))
  println(io, "  PQ->PV           : ", count(ev -> _qlimit_event_kind(ev) === :release, logdata))
  shown = max_rows < 0 ? length(logdata) : min(max_rows, length(logdata))
  println(io, "  rows shown       : ", shown)
  println(io, "  rows omitted     : ", length(logdata) - shown)
  isnothing(full_details) || println(io, "  full details     : ", full_details)
  shown == 0 && return
  println(io, "─"^50)
  println(io, " Iteration │ Bus │ Side │ Change")
  println(io, "─"^50)

  for i in 1:shown
    ev = logdata[i]
    kind = _qlimit_event_kind(ev)
    change = kind === :release ? "PQ->PV" : kind === :group_revert ? "PV->PQ (group)" : "PV->PQ"
    @printf(io, " %9d │ %3d │ %-4s │ %s\n", ev.iter, ev.bus, String(ev.side), change)
  end

  println(io, "─"^50)
end

Base.@kwdef struct QLimitEvent
  iter::Int
  bus::Int
  side::Symbol   # :min | :max
  # :hit for a limit hit; :group_revert when a jointly released group of
  # twice-clamped buses goes back to its limits (#475). A revert is a switch
  # but not a hit: the per-bus hit counts (one-retry guard, max_switches)
  # leave it out. :release for a PQ->PV release, `side` is the limit the
  # bus was clamped at; up to 0.50 the log held the PV->PQ switches only
  kind::Symbol = :hit
end

"""
    logQLimitHit!(net, iter, bus, side)

Record a Q-limit hit (bus, iteration, min or max side) in the network's
Q-limit log.
"""
function logQLimitHit!(net::Net, iter::Int, bus::Int, side::Symbol; kind::Symbol = :hit)
  push!(net.qLimitLog, QLimitEvent(iter = iter, bus = bus, side = side, kind = kind))
  net.qLimitEvents[bus] = side
end

# The kind of a log entry; net.qLimitLog is a Vector{Any}, and an entry
# written by other code without a `kind` field (a test pushes a plain
# NamedTuple) counts as a hit, as every entry did before the field existed.
_qlimit_event_kind(ev)::Symbol = hasproperty(ev, :kind) ? ev.kind : :hit

# The PV->PQ switches of the log (hits and group reverts) without the
# releases; the run summaries report this number as "PV->PQ events".
qlimit_pv_to_pq_event_count(net::Net)::Int = count(ev -> _qlimit_event_kind(ev) !== :release, net.qLimitLog)

function qlimit_switch_count(net::Net, bus::Int)::Int
  return count(ev -> ev.bus == bus && ev.kind === :hit, net.qLimitLog)
end

function qlimit_switch_counts(net::Net)
  counts = Dict{Int,Int}()
  for ev in net.qLimitLog
    ev.kind === :hit || continue
    counts[ev.bus] = get(counts, ev.bus, 0) + 1
  end
  return counts
end

"Returns the last iteration number where `bus` hit a Q-limit, or `nothing`."
function lastQLimitIter(net::Net, bus::Int)
  for i = length(net.qLimitLog):-1:1
    # a release is not a hit; the cooldown measures from the last clamp
    if net.qLimitLog[i].bus == bus && net.qLimitLog[i].kind !== :release
      return net.qLimitLog[i].iter
    end
  end
  return nothing
end

"""
    resetQLimitLog!(net)

Clear the Q-limit event log of the network.
"""
function resetQLimitLog!(net::Net)
  empty!(net.qLimitEvents)
  empty!(net.qLimitLog)
  return nothing
end

function snapshotPVQLimits!(net::Net)
  qmin_pu, qmax_pu = getQLimits_pu(net)
  empty!(net.qLimitInitialPVRows)
  for (bus, node) in enumerate(net.nodeVec)
    getNodeType(node) == PV || continue
    qmin = ((bus <= length(qmin_pu)) && isfinite(qmin_pu[bus])) ? (qmin_pu[bus] * net.baseMVA) : -Inf
    qmax = ((bus <= length(qmax_pu)) && isfinite(qmax_pu[bus])) ? (qmax_pu[bus] * net.baseMVA) : Inf
    push!(net.qLimitInitialPVRows, (bus = bus, qmin_MVAr = qmin, qmax_MVAr = qmax))
  end
  return net.qLimitInitialPVRows
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
function printPVQLimitsTable(net::Net; io::IO = stdout, max_rows::Int = 30, full_details::Union{Nothing,AbstractString} = nothing)
  rows = isempty(net.qLimitInitialPVRows) ? snapshotPVQLimits!(net) : net.qLimitInitialPVRows

  if isempty(rows)
    println(io, "PV Q-limits (MVAr): no PV buses in pre-solve snapshot.")
    return nothing
  end

  shown = max_rows < 0 ? length(rows) : min(max_rows, length(rows))
  println(io, "PV Q-limits before PF run (MVAr):")
  println(io, "  rows total       : ", length(rows))
  println(io, "  rows shown       : ", shown)
  println(io, "  rows omitted     : ", length(rows) - shown)
  isnothing(full_details) || println(io, "  full details     : ", full_details)
  if shown == 0
    println(io, "  preview          : summary-only")
    return nothing
  end
  println(io, "PV->PQ switch comparisons use per-unit internally; switch diagnostics report both p.u. and MVAr.")
  println(io, "──────────────────────────────────────────────")
  println(io, " Bus │      Qmin [MVAr] │      Qmax [MVAr]")
  println(io, "──────────────────────────────────────────────")
  for i = 1:shown
    row = rows[i]
    bus = row.bus
    qmin = row.qmin_MVAr
    qmax = row.qmax_MVAr
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
"""
    has_q_limits(qmin_pu, qmax_pu, i) -> Bool

Whether bus `i` carries a finite reactive band.
"""
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

  for bus = 1:n
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
function printFinalLimitValidation(net::Net; q_headroom::Float64 = 0.20, io::IO = stdout, converged::Bool = true)
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

  if converged
    println(io, "Final Q-limit validation")
  else
    println(io, "Last-iteration Q-limit diagnostic (NR did not converge; values are not a valid final solution)")
  end
  if isempty(qrows)
    @printf(io, "  Q-limits: no violations beyond %.1f%% headroom.\n", 100 * max(q_headroom, 0.0))
  else
    println(io, "  Q-limit violations (MVAr):")
    println(io, "  Displayed values are MVAr. Switching comparisons use internal p.u. values.")
    println(io, "  Bus │ Qgen [MVAr] │ Qmin [MVAr] │ Qmax [MVAr] │ Violation [MVAr]")
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

  qv = printQVCharacteristicCheck(net; io = io, converged = converged)

  return (q_violations = length(qrows), v_violations = length(vrows), qv_violations = length(qv))
end

"""
    _qv_effective_qgen(net, bus, limit_mvar = NaN) -> Float64

The reactive power a clamped machine at `bus` really injects, in MVAr. Two
places hold it and each enforcement mode writes only one of them: the active
set writes the clamp back to the node (`_qƩGen`) and leaves the specification
alone, the classical outer loop applies the clamp in the specification (what
`_effective_bus_power_components` reads) and leaves the node alone. Reading
one of them alone therefore misses every machine the OTHER mode clamped: on
the Zeng case the node said 54.17 MVAr at a 50 MVAr limit under the classical
loop and the check reported nothing, and the specification said 42.4 under the
active set (measured). With a `limit_mvar` given, the candidate
closer to that limit is the one the solver enforced; a re-enabled machine
sits on neither, and stays invisible to the check as it should. A
Q(U)-controlled machine has its own source, the characteristic.
"""
function _qv_effective_qgen(net::Net, bus::Int, limit_mvar::Float64 = NaN)::Float64
  node = net.nodeVec[bus]
  controlled = any(ps -> getPosumerBusIndex(ps) == bus && has_qu_controller(ps), net.prosumpsVec)
  _, q_spec, _, _ = _effective_bus_power_components(net, bus)
  q_spec = Float64(q_spec)
  # a bus with a Q(U) unit: the node is not an independent candidate there,
  # a converged solve writes the same characteristic value to _qƩGen of an
  # unclamped PQ bus (_write_controlled_generation!)
  controlled && return q_spec
  q_node = isnothing(node._qƩGen) ? NaN : Float64(node._qƩGen)
  isfinite(q_node) || return q_spec
  isfinite(q_spec) || return q_node
  isfinite(limit_mvar) || return q_node
  return abs(q_node - limit_mvar) <= abs(q_spec - limit_mvar) ? q_node : q_spec
end

"""
    qvCharacteristicViolations(net; band_pu = 1e-4) -> Vector{NamedTuple}

Generators that ended a solve AT a reactive limit while their voltage sits on
the WRONG side of their own setpoint:

- `Q = Qmax` with `Vm > Vset` (condition (7) of [`qlimitStabilityReport`](@ref) fails), or
- `Q = Qmin` with `Vm < Vset` (condition (6) fails).

Such a point satisfies the power-flow equations and still contradicts the
physical Q-V characteristic: a machine that is already at its upper reactive
limit cannot hold a voltage that is above its setpoint. The literature calls
this a *non-physical* solution (Zeng, Chiang, Neves, Alberto, IJEPES 147
(2023) 108905) or an *anomalous* one (Sundaresh, Rao 2015), and it is a
property of the PV/PQ switching strategy, not of the network.

The rows are a view of [`qlimitStabilityReport`](@ref), not a check of their
own: the buses with a Q-limit event whose report row sits on the limit of
that event and fails the condition of that limit. Q, the set point and the
tolerances (1e-3 MVAr on the limit, 1e-6 pu off the set point) are those of
the stability report, so the two verdicts cannot disagree.

Each row carries `(bus, side, vm_pu, vset_pu, dv_pu, q_MVAr, limit_MVAr,
significant)`. `dv_pu` is `Vm - Vset` and therefore signed; `significant` is
`abs(dv_pu) > band_pu`. The band exists because a converged solve leaves
voltages within its tolerance of the setpoint, so a bus close to its
setpoint can fall to either side by rounding. Report BOTH counts: the raw one
is the literature's definition, the significant one is what a reader can act
on.

Returns an empty vector when no bus hit a limit. The caller decides whether
the state comes from a converged solve (`converged = false` at the call site).
"""
function qvCharacteristicViolations(net::Net; band_pu::Float64 = 1.0e-4)
  rows = NamedTuple[]
  isempty(net.qLimitEvents) && return rows
  for r in qlimitStabilityReport(net)
    side = get(net.qLimitEvents, r.bus, nothing)
    (side === :max || side === :min) || continue
    # only a bus that really sits ON the limit of its event can be
    # non-physical; a stale event whose machine moved back inside its band
    # is not a violation
    r.reached_limit === side || continue
    wrong_side = side === :max ? !r.check_7 : !r.check_6
    wrong_side || continue
    dv = r.u_pu - r.uset_pu
    push!(rows, (bus = r.bus, side = side, vm_pu = r.u_pu, vset_pu = r.uset_pu, dv_pu = dv,
      q_MVAr = r.q_MVar, limit_MVAr = side === :max ? r.qmax_MVar : r.qmin_MVar, significant = abs(dv) > band_pu))
  end
  return rows
end

"""
    printQVCharacteristicCheck(net; io=stdout, band_pu=1e-4, converged=true) -> Vector{NamedTuple}

Print the [`qvCharacteristicViolations`](@ref) of the current solution and
return them. On a non-converged solve the check is skipped with a note: the
state is not a solution, so the question does not apply.
"""
function printQVCharacteristicCheck(net::Net; io::IO = stdout, band_pu::Float64 = 1.0e-4, converged::Bool = true)
  if !converged
    println(io, "  Q-V characteristic: skipped (no converged solution).")
    return NamedTuple[]
  end
  rows = qvCharacteristicViolations(net; band_pu = band_pu)
  if isempty(rows)
    println(io, "  Q-V characteristic: no non-physical generator states.")
    return rows
  end
  strong = count(r -> r.significant, rows)
  @printf(io, "  Q-V characteristic: %d non-physical generator state(s), %d beyond %.0e pu.\n", length(rows), strong, band_pu)
  println(io, "  A machine at its reactive limit whose voltage sits on the wrong side of its")
  println(io, "  setpoint solves the equations but contradicts the physical Q-V curve.")
  println(io, "  Bus │ Side │        Vm │      Vset │        dV │  Q [MVAr] │ limit [MVAr]")
  for r in rows
    @printf(io, " %4d │ %-4s │ %9.5f │ %9.5f │ %+9.5f │ %9.3f │ %12.3f%s\n",
      r.bus, String(r.side), r.vm_pu, r.vset_pu, r.dv_pu, r.q_MVAr, r.limit_MVAr, r.significant ? "" : "  (within band)")
  end
  return rows
end

# --- Q-limit stability check of the final solution ---------------------------
#
# Zhao, Chiang, Li, Ju, "On PV-PQ bus type switching logic in power flow
# computation", PSCC 2008, conditions (6) and (7) of a stable solution:
#   (6) (U - U_set)(Q_max - Q) >= 0,   (7) (U - U_set)(Q - Q_min) <= 0.
# A machine at Q_max must not sit above its set point ((7) fails), one at
# Q_min not below it ((6) fails), and a machine inside its band must hold the
# set point. qvCharacteristicViolations (machines with a Q-limit event, on the
# limit, on the wrong side) is the subset of this check for the event buses,
# read from its rows; there is one check, not two. Set points come from
# _bus_voltage_setpoints_from_prosumers, bus bands from getQLimits_pu. Q is the value of the result report
# (_bus_power_component_cache: the solver result of a PV bus, the clamp of a
# PQ bus through _qv_effective_qgen), and the value from the solved voltages
# (_qlimit_stability_qgen) on a net with a closed link and at a bus that
# carries a Q(U) generator: there the report gives the regulating machine
# its scheduled qVal, not its solved Q.

"""
    QLIMIT_STABILITY_TOL_Q_MVAR

Reactive power tolerance of the Q-limit stability check (1e-3 MVAr): a
machine within this distance of a limit is on that limit, and a band is
exceeded only beyond it. Fixed, not a configuration key: the check asks
whether the solved machine sits on its limit, which `power_flow.tol` (a
mismatch tolerance) does not answer.
"""
const QLIMIT_STABILITY_TOL_Q_MVAR = 1.0e-3

"""
    QLIMIT_STABILITY_TOL_U_PU

Voltage tolerance of the Q-limit stability check (1e-6 pu): a bus within
this distance of its set point holds it. Fixed like
`QLIMIT_STABILITY_TOL_Q_MVAR`.
"""
const QLIMIT_STABILITY_TOL_U_PU = 1.0e-6

# sign of x with a dead band: 0 within tol, +1 / -1 outside; an infinite
# margin (an open limit) counts as positive, so no Inf * 0 can appear
_qlimit_stability_sign(x::Float64, tol::Float64)::Int = abs(x) <= tol ? 0 : (x > 0.0 ? 1 : -1)

# The buses the check judges: a generator with a voltage regulator and a set
# point, in service, solved as PV or PQ, with no reference unit on the bus
# or on a bus joined to it by a closed link. The classic outer loop may hand
# the reference to such a bus; it is Slack then and out of the check like
# every reference bus. A link cluster with a reference is one node whose Q
# the reference balances: the Q of a regulated member is only its share of
# the cluster Q, which no limit binds (the 3-bus test net with a machine of
# 1 MVAr band on a bus coupled to the slack bus: reported as a PV machine at
# 10.8 MVAr outside its band before this exclusion). The island references
# of a net without a slack in some island are left out in
# qlimitStabilityReport by the same cluster rule.
function _qlimit_stability_buses(net::Net)::Vector{Int}
  n = length(net.nodeVec)
  reps = _active_link_representative_map(net)
  regulated = falses(n)
  # indexed by the cluster representative: a reference unit, or a bus the
  # solver runs as Slack, anywhere in the cluster
  reference = falses(n)
  for ps in net.prosumpsVec
    isGenerator(ps) || continue
    b = getPosumerBusIndex(ps)
    1 <= b <= n || continue
    isSlack(ps) && (reference[reps[b]] = true)
    (isRegulating(ps) && !isnothing(ps.vm_pu)) && (regulated[b] = true)
  end
  for b = 1:n
    getNodeType(net.nodeVec[b]) == Slack && (reference[reps[b]] = true)
  end
  iso = Set{Int}(net.isoNodes)
  return [b for b = 1:n if regulated[b] && !reference[reps[b]] && !(b in iso) && getNodeType(net.nodeVec[b]) in (PV, PQ)]
end

# Generator Q (MVAr) of the regulated buses `vbuses` from the solved
# voltages, on the network the solver works on (`wnet`, contracted on a net
# with a closed link): the injection of the bus, or of its link cluster,
# plus what the cluster draws itself (_bus_consumption_MVA at the
# representative: the loads of its members and its voltage-dependent
# shunts, which the Y-bus leaves out) minus the generation the report gives
# the other members. A cluster with a reference never gets here
# (_qlimit_stability_buses), so no member's Q is a share of a slack.
#
# Two cases need it:
#  - a closed link: the solve runs on a contracted copy and its results are
#    split back onto the buses, so the voltages are the value that does not
#    depend on that write-back. Before 0.60.0 the write-back left a clamped
#    bus PV under the classic modes and the report took a wrong value:
#    sp_case60 with a 4.0 .. 5.5 MVAr band at Birkkamp_110 (one closed
#    coupler elsewhere in the net) reported 7.389 MVAr, its voltages give
#    5.500 (_sync_merged_qlimit_outcome!); the report agrees since;
#  - a Q(U) generator on the bus: the report keeps the scheduled qVal of the
#    regulating machine (issue #374 keeps the Q(U) value from being
#    overwritten), so the 3-bus test net with a Q(U) unit next to the
#    machine at STATION1 read 34.379 MVAr at a bus clamped at 8.909.
# Elsewhere the report is the solver result itself and agreed with the
# voltages to 4e-6 MVAr on case14_zeng in all three enforcement modes, so
# the Y-bus, the main cost of the check, is built only for these buses.
#
# The value is in the terms of the bus band (_buildQLimits!): a Q(U) unit
# that follows its characteristic on the regulated bus is in the band with
# its Q at the held voltage, so its Q at the solved voltage is replaced by
# that one. At a bus that holds its set point both are the same.
# `wnet` and `wreps` are the result of _merged_pf_net(net), `vset` the set
# points of _bus_voltage_setpoints_from_prosumers(net). Returns an empty
# dictionary when `vbuses` is empty.
function _qlimit_stability_qgen(net::Net, vbuses::Vector{Int}, cache, wnet::Net, wreps::Vector{Int}, vset::Vector{Float64})::Dict{Int,Float64}
  out = Dict{Int,Float64}()
  isempty(vbuses) && return out
  n = length(net.nodeVec)
  base = net.baseMVA
  Yred = createYBUS(net = wnet, sparse = true, printYBUS = false)
  Y = size(Yred, 1) == n ? Yred : _expand_ybus_for_isolated_nodes(Yred, n, wnet.isoNodes)
  iso = Set{Int}(wnet.isoNodes)
  V = zeros(ComplexF64, n)
  for k = 1:n
    nd = wnet.nodeVec[k]
    (k in iso || isnothing(nd._vm_pu) || isnothing(nd._va_deg)) && continue
    V[k] = Float64(nd._vm_pu) * cis(deg2rad(Float64(nd._va_deg)))
  end
  I = Y * V
  members = Dict{Int,Vector{Int}}()
  for b = 1:n
    push!(get!(members, wreps[b], Int[]), b)
  end
  # the prosumers and shunts of a cluster sit at its representative in wnet
  _, q_cons = _bus_consumption_MVA(wnet, n)
  # Q(U) followers of each bus: their Q at the held voltage minus their Q at
  # the solved voltage, the shift from the solved total to the band terms
  follow = Dict{Int,Float64}()
  wanted = Set{Int}(vbuses)
  for ps in net.prosumpsVec
    isGenerator(ps) || continue
    b = getPosumerBusIndex(ps)
    (b in wanted && _qband_follows_characteristic(ps, true)) || continue
    vm = _safe_vm(abs(V[wreps[b]]))
    follow[b] = get(follow, b, 0.0) + (first(evaluate_controller(ps.quController, vset[b])) - first(evaluate_controller(ps.quController, vm))) * base
  end
  for b in vbuses
    r = wreps[b]
    M = members[r]
    q_inj = imag(V[r] * conj(I[r])) * base
    q_load = q_cons[r]
    q_other = sum((cache[m].q_gen for m in M if m != b && haskey(cache, m)); init = 0.0)
    out[b] = q_inj + q_load - q_other + get(follow, b, 0.0)
  end
  return out
end

"""
    qlimitStabilityReport(net; locked_by = Dict{Int,Symbol}()) -> Vector{NamedTuple}

The stability conditions of Zhao, Chiang, Li, Ju (PSCC 2008) checked on the
solved network, one row per voltage-regulated generator bus (a generator with
a voltage regulator and a set point, no reference unit on the bus or on a bus
joined to it by a closed link, not the reference the solver gives an island
without a slack), whether it ended PV or PQ, switched or locked:

- (6) `(U - U_set)(Q_max - Q) >= 0`: holds at `Q_max` (the factor is zero)
  and fails for a machine at `Q_min` below its set point;
- (7) `(U - U_set)(Q - Q_min) <= 0`: holds at `Q_min` and fails for a
  machine at `Q_max` above its set point.

A bus inside its band (neither at a limit) meets both only at its set point.
A bus that ended PV is `ok` when `Q_min <= Q <= Q_max`; a bus that ended PQ
when (6), (7) and the band hold. Tolerances: `QLIMIT_STABILITY_TOL_Q_MVAR`
(1e-3 MVAr) for Q against a limit, `QLIMIT_STABILITY_TOL_U_PU` (1e-6 pu) for
U against U_set.

Q is the generator Q of the bus (MVAr) as the result report shows it
(`q_source = :report`). On a net with a closed link, and at a bus that
carries a Q(U) generator, it comes from the solved voltages instead
(`q_source = :voltages`): the injection of the bus on the contracted network
the solver works on plus what it draws itself (loads and voltage-dependent
shunts), minus the generation of the other buses of its link cluster, with
a Q(U) unit of the bus at its Q at the held voltage, as the bus band counts
it. The report value is kept as `q_report_MVar`; it differs at a Q(U) bus,
where the report keeps the scheduled Q of the regulating machine. Q_min and Q_max are the
bus bands of `getQLimits_pu`. U_set is the regulator set point
(`_bus_voltage_setpoints_from_prosumers`), not the start voltage. `end_type`
is `:PQ` for a bus that ended PQ or carries a Q-limit event.

Row fields: `bus`, `bus_name`, `end_type`, `q_MVar`, `qmin_MVar`,
`qmax_MVar`, `u_pu`, `uset_pu`, `reached_limit` (`:max`, `:min`, `:none`),
`check_6`, `check_7`, `in_band`, `ok`, `locked_by` (from the `locked_by`
lookup, `:none` for a bus not in it), `q_source` and `q_report_MVar`. The rows are in bus
order. The caller decides whether the state is a solution; the function
reads it as it is.
"""
function qlimitStabilityReport(net::Net; locked_by::AbstractDict{Int,Symbol} = Dict{Int,Symbol}())
  rows = NamedTuple[]
  buses = _qlimit_stability_buses(net)
  isempty(buses) && return rows
  wnet, wreps, merged = _merged_pf_net(net)
  # an island without a slack runs with the bus of its best unit as its
  # reference (detect_ac_islands, the choice the solver makes); that machine
  # balances the island and its Q is free like that of a slack
  # (two_islands_prio: Fennwiese_110 at 41.8 MVAr in a +-30 MVAr band)
  if length(_active_ac_island_components(wnet)) > 1
    # detect_ac_islands refreshes the bus types from the units first, which
    # undid the clamps of the classic modes on the solved net (measured on
    # case14_zeng), so it reads a copy; one island needs no lookup
    island_refs = Set{Int}(row.chosen_ref_bus for row in detect_ac_islands(deepcopy(wnet); promote_generators = true).rows if row.chosen_ref_bus > 0)
    filter!(b -> !(wreps[b] in island_refs), buses)
  end
  isempty(buses) && return rows
  base = net.baseMVA
  vset = _bus_voltage_setpoints_from_prosumers(net)
  qmin_pu, qmax_pu = getQLimits_pu(net)
  cache = _bus_power_component_cache(net)
  # Q from the voltages: every bus on a net with a closed link, else the
  # buses with a Q(U) generator (see _qlimit_stability_qgen)
  vbuses = if merged
    buses
  else
    qu_buses = Set{Int}(getPosumerBusIndex(ps) for ps in net.prosumpsVec if isGenerator(ps) && has_qu_controller(ps))
    [b for b in buses if b in qu_buses]
  end
  q_solved = _qlimit_stability_qgen(net, vbuses, cache, wnet, wreps, vset)
  names = _bus_name_by_idx(net)
  tol_q = QLIMIT_STABILITY_TOL_Q_MVAR
  tol_u = QLIMIT_STABILITY_TOL_U_PU
  for b in buses
    node = net.nodeVec[b]
    isnothing(node._vm_pu) && continue
    u = Float64(node._vm_pu)
    us = vset[b]
    q = get(q_solved, b, Float64(cache[b].q_gen))
    qmin = b <= length(qmin_pu) ? qmin_pu[b] * base : -Inf
    qmax = b <= length(qmax_pu) ? qmax_pu[b] * base : Inf
    end_type = (haskey(net.qLimitEvents, b) || getNodeType(node) == PQ) ? :PQ : :PV
    s_u = _qlimit_stability_sign(u - us, tol_u)
    s_max = _qlimit_stability_sign(qmax - q, tol_q)
    s_min = _qlimit_stability_sign(q - qmin, tol_q)
    check_6 = s_u * s_max >= 0
    check_7 = s_u * s_min <= 0
    in_band = s_max >= 0 && s_min >= 0
    # on both limits (a zero-width band) the side of the event names it
    reached = if s_max == 0 && s_min == 0
      get(net.qLimitEvents, b, :max)
    elseif s_max == 0
      :max
    elseif s_min == 0
      :min
    else
      :none
    end
    ok = end_type === :PV ? in_band : (in_band && check_6 && check_7)
    push!(rows, (bus = b, bus_name = _effective_bus_name(names, net, b), end_type = end_type, q_MVar = q, qmin_MVar = qmin, qmax_MVar = qmax,
      u_pu = u, uset_pu = us, reached_limit = reached, check_6 = check_6, check_7 = check_7, in_band = in_band, ok = ok,
      locked_by = get(locked_by, b, :none), q_source = haskey(q_solved, b) ? :voltages : :report, q_report_MVar = Float64(cache[b].q_gen)))
  end
  return rows
end

# Why a regulated bus could not switch freely at the end of the run, from
# the run configuration and the Q-limit log (no solver state survives the
# solve). Precedence: the configured lock list, the guard, then the active-set
# locks. Under the classic modes only the list locks a bus; their clamp is
# the mode itself (one-sided, no release) and carries :none.
#  - :lock_pv_to_pq_buses: listed in power_flow.qlimits.lock_pv_to_pq_buses;
#  - :guard_lock_pq: active set with the guard, a band narrower than
#    guard_min_q_range_pu whose mode runs the bus PQ (zero range or
#    prefer_pq) or locks it at the first limit it reaches (lock_pq), and the
#    bus ended clamped;
#  - :max_switches: active set, guard_freeze_after_repeated_switching, at
#    least guard_max_switches hits, the bus ended PV (the freeze keeps it);
#  - :oscillation: active set, two or more hits and the bus ended clamped:
#    the one-retry lock of active_set_q_limits! releases it no more.
function _qlimit_stability_locks(net::Net, qlim::QLimitConfig, mode::Symbol, buses::Vector{Int})::Dict{Int,Symbol}
  locks = Dict{Int,Symbol}()
  listed = Set{Int}(qlim.lock_pv_to_pq_buses)
  active_set = mode === :active_set
  qmin_pu, qmax_pu = active_set ? getQLimits_pu(net) : (Float64[], Float64[])
  hits = active_set ? qlimit_switch_counts(net) : Dict{Int,Int}()
  for b in buses
    if b in listed
      locks[b] = :lock_pv_to_pq_buses
      continue
    end
    active_set || continue
    clamped = haskey(net.qLimitEvents, b) || getNodeType(net.nodeVec[b]) == PQ
    if qlim.guard && clamped && b <= length(qmin_pu) && b <= length(qmax_pu) && isfinite(qmin_pu[b]) && isfinite(qmax_pu[b])
      width = abs(qmax_pu[b] - qmin_pu[b])
      guard_mode = width <= eps(Float64) ? qlim.guard_zero_range_mode : qlim.guard_narrow_range_mode
      if width < qlim.guard_min_q_range_pu && guard_mode in (:lock_pq, :prefer_pq)
        locks[b] = :guard_lock_pq
        continue
      end
    end
    h = get(hits, b, 0)
    if !clamped && qlim.guard_freeze_after_repeated_switching && h >= max(qlim.guard_max_switches, 1)
      locks[b] = :max_switches
    elseif clamped && h >= 2
      locks[b] = :oscillation
    end
  end
  return locks
end

"""
    qlimit_stability_line(status, rows) -> String

The Q-limit stability check of a run in one line, for the result print,
`run.log`, the result metadata and the Web UI result page:
`"Q-limit stability violations: n"`, followed by the buses that fail when
`n > 0`. `status` is `:ok` or `:violations` for a checked run, `:not_checked`
for a run without a converged solution, `:last_iterate` for a run stopped at
`power_flow.max_iter` (the check of its last iterate, which is no solution;
after `power_flow.dc.fallback` the net no longer holds that iterate and the
run is `:not_checked`) and `:qlimits_disabled` for a run without Q-limit handling.
"""
function qlimit_stability_line(status::Symbol, rows)::String
  status === :not_checked && return "Q-limit stability violations: not checked (no converged solution)"
  status === :qlimits_disabled && return "Q-limit stability violations: not checked (Q-limit handling disabled)"
  bad = [r for r in rows if !r.ok]
  buses = isempty(bad) ? "" : string(" (bus ", join((r.bus_name for r in bad), ", "), ")")
  status === :last_iterate && return string("Q-limit stability violations: ", length(bad), buses, " at the last iterate (no converged solution)")
  return string("Q-limit stability violations: ", length(bad), buses)
end

# one failing row as the warning names it
function _qlimit_stability_row_text(r)::String
  where_ = r.reached_limit === :max ? "at Qmax" : r.reached_limit === :min ? "at Qmin" : r.in_band ? "inside its band" : "outside its band"
  q_text = r.in_band ? "" : @sprintf(", Q %.3f MVAr outside %.3f .. %.3f", r.q_MVar, r.qmin_MVar, r.qmax_MVar)
  lock_text = r.locked_by === :none ? "" : string(", locked by ", r.locked_by)
  return @sprintf("%s (%s %s, U %.5f pu, U_set %.5f pu%s%s)", r.bus_name, String(r.end_type), where_, r.u_pu, r.uset_pu, q_text, lock_text)
end

# The check after a converged run (_build_sparlectra_result): diagnostics
# fields, the line, and the warning text ("" without a violation). A run
# without a converged solution and a run without Q-limit handling are not
# checked: there is no final state, or no switching to judge. A run stopped
# at max_iter (reason :iteration_limit, :type_identification_divergence or
# :numerical_divergence) gets the rows of its last iterate with the status
# :last_iterate and no warning: the run has failed already, the rows show
# where the bus types stood when it stopped. A run that went on to the DC
# fallback is not checked. `qlim` is the
# configuration of the attempt that produced the result (in auto mode the
# last stage, which may have locked buses or switched the handling off);
# `st` is the solver status of that run (`nothing` off the rectangular
# path), whose final_q_check_status also says when the solver ran without
# Q-limit handling.
function _qlimit_stability_diagnostics(net::Net, qlim::QLimitConfig, numerical_converged::Bool, st)
  disabled = qlim.ignore_q_limits || (st !== nothing && hasproperty(st, :final_q_check_status) && st.final_q_check_status === :qlimits_disabled)
  # after power_flow.dc.fallback the net holds the restored start voltages
  # or the DC state (vm = 1 pu), not the last iterate: not checked
  dc_fallback = st !== nothing && hasproperty(st, :dc_fallback_applied)
  max_iter_stop = !numerical_converged && !dc_fallback && st !== nothing && hasproperty(st, :reason) && Symbol(st.reason) in QLIMIT_MAX_ITER_STOP_REASONS
  status = !numerical_converged ? (max_iter_stop && !disabled ? :last_iterate : :not_checked) : disabled ? :qlimits_disabled : :checked
  if !(status in (:checked, :last_iterate))
    return (qlimit_stability_status = status, qlimit_stability_violations = 0, qlimit_stability_report = NamedTuple[], qlimit_stability_line = qlimit_stability_line(status, NamedTuple[]), qlimit_stability_warning = "")
  end
  # the classic outer loop records its mode; the auto mode may end on another
  # strategy than the configuration names
  mode = st !== nothing && hasproperty(st, :qlimit_enforcement_mode) ? Symbol(st.qlimit_enforcement_mode) : qlim.enforcement_mode
  buses = _qlimit_stability_buses(net)
  rows = qlimitStabilityReport(net; locked_by = _qlimit_stability_locks(net, qlim, mode, buses))
  bad = [r for r in rows if !r.ok]
  status === :last_iterate && return (qlimit_stability_status = status, qlimit_stability_violations = length(bad), qlimit_stability_report = rows, qlimit_stability_line = qlimit_stability_line(status, rows), qlimit_stability_warning = "")
  status = isempty(bad) ? :ok : :violations
  warning = isempty(bad) ? "" :
            string("Q-limit stability check: ", length(bad), " voltage-regulated bus(es) violate the stability conditions of the final solution ",
    "(a machine at Qmax above its set point, at Qmin below it, inside its band off it, or a PV machine outside its band): ",
    join((_qlimit_stability_row_text(r) for r in bad), "; "))
  return (qlimit_stability_status = status, qlimit_stability_violations = length(bad), qlimit_stability_report = rows, qlimit_stability_line = qlimit_stability_line(status, rows), qlimit_stability_warning = warning)
end

"""
    printQLimitStabilityLine(net; io=stdout, converged=true) -> String

Print the line of the Q-limit stability check ([`qlimit_stability_line`](@ref))
and return it. A run through `run_sparlectra` stores its line in the solver
status, with the lock reasons of its configuration; on a network solved by
`runpf!` alone the check runs here on the solved state. Part of the result
print.
"""
function printQLimitStabilityLine(net::Net; io::IO = stdout, converged::Bool = true)
  st = rectangular_pf_status(net)
  line = if st !== nothing && hasproperty(st, :qlimit_stability_line)
    String(st.qlimit_stability_line)
  elseif !converged
    qlimit_stability_line(:not_checked, NamedTuple[])
  elseif st !== nothing && hasproperty(st, :final_q_check_status) && st.final_q_check_status === :qlimits_disabled
    qlimit_stability_line(:qlimits_disabled, NamedTuple[])
  else
    rows = qlimitStabilityReport(net)
    qlimit_stability_line(any(r -> !r.ok, rows) ? :violations : :ok, rows)
  end
  println(io, "  ", line)
  return line
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
        get_vm_pu=nothing,
        get_vset_pu=nothing,
        v_hyst_pu::Float64=1e-4,
        converged_iterate::Bool=false,
        verbose::Int=0,
        io::IO=stdout,
    ) -> (changed::Bool, reenabled::Bool)

Core PV/Q-limit active-set logic shared by solvers.

The PQ->PV release of a clamped bus is decided on the VOLTAGE side (#375,
the "back off" rule of Sundaresh and Rao): a bus clamped at Qmax is
released when `Vm > Vset + v_hyst_pu`, one clamped at Qmin when
`Vm < Vset - v_hyst_pu`. The reactive injection of a clamped bus converges
to the limit itself, so the former test "Q strictly inside the band"
could never be true at a converged point and non-physical solutions were
kept. The Q test remains the fallback for a bus without a stored side or
without the voltage callbacks. Cooldown and the one-retry guard apply to
both forms.

The one-retry guard keeps a bus clamped once it has hit its limit twice,
against chattering in unsettled iterates. On a converged iterate
(`converged_iterate = true`) on which nothing else switched, all such buses
whose voltage sits on the release side are released together, once: a
clamp taken in an early iterate must not outlive the state that caused it
(#475: seven machines stayed at Qmin with the voltage below their
setpoint). When one of the group violates its limit again, the whole group
goes back to its limits (logged as `kind = :group_revert`, not counted as a
hit) and the solve makes no further extra release; released one by one or
all without the revert, a group of mutually supporting machines diverged
(sp_case118, outage of branch 68-81).

Callbacks:
- get_qreq_pu(bus) -> Float64
- get_vm_pu(bus) -> Float64, get_vset_pu(bus) -> Float64 (voltage-side release)
- is_pv(bus) -> Bool
- make_pq!(bus, q_clamp_pu::Float64, side::Symbol)  # side = :min/:max
- make_pv!(bus)
- on_violation!(bus, qreq_pu::Float64, side::Symbol, q_clamp_pu::Float64) -> Bool
  (optional; return true if violation was handled without PV->PQ fallback)
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
  get_vm_pu = nothing,
  get_vset_pu = nothing,
  v_hyst_pu::Float64 = 1e-4,
  converged_iterate::Bool = false,
  lock_pv_to_pq_buses::AbstractVector{Int} = Int[],
  on_violation! = nothing,
  verbose::Int = 0,
  io::IO = stdout,
  max_console_rows::Int = 30,
  qlimit_guard_max_switches::Int = typemax(Int),
  qlimit_guard_freeze_after_repeated_switching::Bool = true,
  qlimit_guard_violation_mode::Symbol = :delayed_switch,
  qlimit_guard_violation_threshold_pu::Float64 = 1e-4,
  lock_mask_buffer::Union{Nothing,BitVector} = nothing,
)
  qlimit_guard_violation_mode in (:delayed_switch, :lock_pq, :ignore) || error("Unsupported qlimit_guard_violation_mode=$(qlimit_guard_violation_mode). Supported: :delayed_switch, :lock_pq, :ignore.")
  qlimit_guard_violation_threshold_pu >= 0.0 || error("qlimit_guard_violation_threshold_pu must be >= 0 (got $(qlimit_guard_violation_threshold_pu)).")
  changed = false
  reenabled = false
  printed_events = 0
  omitted_events = 0
  lock_mask = isnothing(lock_mask_buffer) ? falses(nb) : lock_mask_buffer
  length(lock_mask) == nb || throw(ArgumentError("lock_mask_buffer length mismatch: expected $nb, got $(length(lock_mask))"))
  fill!(lock_mask, false)
  @inbounds for bus in lock_pv_to_pq_buses
    if 1 <= bus <= nb
      lock_mask[bus] = true
    end
  end

  # Per-bus Q-limit event count. The guard below and the re-enable branch both
  # ask "how often has this bus switched so far", which used to be a `count`
  # over the whole event log PER BUS PER ITERATION: on a 13659-bus network
  # that made the Q-limit handling 70 percent of the power flow (7731 of
  # 10957 profile samples, 60 s for 8 iterations). The log GROWS while this
  # function runs (every switch appends), so the counts cannot be built once
  # up front - they are carried forward from where the last lookup stopped,
  # which is exactly what the repeated scan computed, amortized to O(1).
  switch_counts = Dict{Int,Int}()
  counted_upto = 0
  group_reverted = false
  switch_count = function (bus::Int)
    @inbounds while counted_upto < length(net.qLimitLog)
      counted_upto += 1
      ev = net.qLimitLog[counted_upto]
      if ev.kind === :hit
        switch_counts[ev.bus] = get(switch_counts, ev.bus, 0) + 1
      elseif ev.kind === :group_revert
        group_reverted = true
      end
    end
    return get(switch_counts, bus, 0)
  end

  # --- PV -> PQ ------------------------------------------------------------
  extra_group_failed = false
  @inbounds for bus = 1:nb
    is_pv(bus) || continue
    lock_mask[bus] && continue
    qreq = get_qreq_pu(bus)

    has_q_limits(qmin_pu, qmax_pu, bus) || continue

    has_hi = (bus <= length(qmax_pu)) && isfinite(qmax_pu[bus])
    has_lo = (bus <= length(qmin_pu)) && isfinite(qmin_pu[bus])

    side = :none
    qclamp = 0.0
    physical_side = :none
    physical_violation = 0.0
    if has_hi && (qreq > qmax_pu[bus])
      physical_side = :max
      physical_violation = qreq - qmax_pu[bus]
    elseif has_lo && (qreq < qmin_pu[bus])
      physical_side = :min
      physical_violation = qmin_pu[bus] - qreq
    end
    # Use the same hysteresis margin for the initial PV->PQ decision as for
    # later PQ->PV re-enable checks.  With q_hyst_pu = 0.0 this preserves the
    # historical exact-limit behavior; positive margins avoid switching on
    # transient or roundoff-sized Q overshoots during difficult flat starts.
    hi_switch = has_hi ? qmax_pu[bus] + q_hyst_pu : Inf
    lo_switch = has_lo ? qmin_pu[bus] - q_hyst_pu : -Inf

    if has_hi && (qreq > hi_switch)
      side = :max
      qclamp = qmax_pu[bus]
    elseif has_lo && (qreq < lo_switch)
      side = :min
      qclamp = qmin_pu[bus]
    end
    if side == :none && qlimit_guard_violation_mode == :lock_pq && physical_violation > qlimit_guard_violation_threshold_pu
      side = physical_side
      qclamp = side == :max ? qmax_pu[bus] : qmin_pu[bus]
    end
    side == :none && continue

    if qlimit_guard_freeze_after_repeated_switching && switch_count(bus) >= qlimit_guard_max_switches
      if verbose > 0
        if max_console_rows < 0 || printed_events < max_console_rows
          @printf(io, "PV->PQ Bus %d: switching suppressed after %d Q-limit event(s); bus remains in its current active-set state (it=%d)\n", bus, switch_count(bus), it)
          printed_events += 1
        else
          omitted_events += 1
        end
      end
      continue
    end

    # a PV bus with exactly two hits can only come from the extra release of
    # a converged iterate (#475; the normal release stops at one hit)
    in_group = switch_count(bus) == 2
    in_group && (extra_group_failed = true)

    handled = false
    if !isnothing(on_violation!)
      handled = on_violation!(bus, qreq, side, qclamp)
    end

    if handled
      changed = true
    elseif side == :max
      make_pq!(bus, qclamp, side)
      logQLimitHit!(net, it, bus, side; kind = in_group ? :group_revert : :hit)
      changed = true
      if verbose > 0
        if max_console_rows < 0 || printed_events < max_console_rows
          @printf(io, "PV->PQ Bus %d: Q=%.6f pu (%.6f MVAr) > Qmax=%.6f pu (%.6f MVAr) (it=%d)\n", bus, qreq, qreq * net.baseMVA, qmax_pu[bus], qmax_pu[bus] * net.baseMVA, it)
          printed_events += 1
        else
          omitted_events += 1
        end
      end
    else
      make_pq!(bus, qclamp, side)
      logQLimitHit!(net, it, bus, side; kind = in_group ? :group_revert : :hit)
      changed = true
      if verbose > 0
        if max_console_rows < 0 || printed_events < max_console_rows
          @printf(io, "PV->PQ Bus %d: Q=%.6f pu (%.6f MVAr) < Qmin=%.6f pu (%.6f MVAr) (it=%d)\n", bus, qreq, qreq * net.baseMVA, qmin_pu[bus], qmin_pu[bus] * net.baseMVA, it)
          printed_events += 1
        else
          omitted_events += 1
        end
      end
    end
  end

  # The extra group of #475 is released together; when one of its machines
  # violates again, the release did not hold for the group (the machines
  # support each other's voltage), so every machine still PV in it goes back
  # to the limit it had before. The revert is logged as such, not as a hit,
  # and ends the extra releases of this solve; the state is one Newton step
  # from the converged point it came from.
  if extra_group_failed
    last_side = Dict{Int,Symbol}()
    for ev in net.qLimitLog
      ev.kind === :hit && (last_side[ev.bus] = ev.side)
    end
    @inbounds for bus = 1:nb
      (is_pv(bus) && switch_count(bus) == 2) || continue
      side = get(last_side, bus, :none)
      side in (:min, :max) || continue
      make_pq!(bus, side == :max ? qmax_pu[bus] : qmin_pu[bus], side)
      logQLimitHit!(net, it, bus, side; kind = :group_revert)
      changed = true
      if verbose > 0
        @printf(io, "PV->PQ Bus %d: back to Q%s with its group, the joint release did not hold (it=%d)\n", bus, String(side), it)
      end
    end
  end

  # --- Optional PQ -> PV ---------------------------------------------------
  # the extra release of a twice-clamped bus on a converged iterate (#475):
  # every candidate together (the machines of the five ACTIVSg2000 outages
  # hold their voltage only as a group), and only when nothing else switched
  # in this iterate; a group that does not hold is clamped back above
  extra_buses = Int[]
  # one failed group ends the extra releases of this solve (its revert is in
  # the log; switch_count(0) brings the scan up to date)
  switch_count(0)
  extra_allowed = converged_iterate && !group_reverted
  release! = function (bus::Int, side::Symbol, voltage_rule::Bool, qreq::Float64, lo::Float64, hi::Float64)
    make_pv!(bus)
    delete!(net.qLimitEvents, bus)   # clear event after re-enable
    # the release goes into the log with the side the bus was clamped at;
    # not through logQLimitHit!, which records a clamp in qLimitEvents
    push!(net.qLimitLog, QLimitEvent(iter = it, bus = bus, side = side, kind = :release))
    reenabled = true
    if verbose > 0
      if voltage_rule
        @printf(io, "PQ->PV Bus %d: clamped at %s, Vm=%.5f pu %s Vset=%.5f pu (margin %.1e), released\n", bus, String(side), get_vm_pu(bus), side == :max ? "above" : "below", get_vset_pu(bus), v_hyst_pu)
      else
        @printf(io, "PQ->PV Bus %d: Q=%.6f pu (%.6f MVAr) within (%.6f, %.6f) pu ((%.6f, %.6f) MVAr)\n", bus, qreq, qreq * net.baseMVA, lo, hi, lo * net.baseMVA, hi * net.baseMVA)
      end
    end
  end
  if allow_reenable
    @inbounds for bus = 1:nb
      # Guards: was PV originally + has event + currently NOT PV
      pv_orig_mask[bus] || continue
      haskey(net.qLimitEvents, bus) || continue
      is_pv(bus) && continue

      has_q_limits(qmin_pu, qmax_pu, bus) || continue

      # A re-enabled bus that violates its Q limit again is likely chattering
      # between the PQ clamp and the PV voltage constraint. Keep it clamped
      # after the first retry instead of re-enabling it indefinitely. On a
      # converged iterate the twice-clamped buses are candidates for the
      # joint extra release (#475): both clamps may stem from transient
      # iterates, and the converged voltage decides.
      prior_hits = switch_count(bus)
      prior_hits > (extra_allowed ? 2 : 1) && continue
      extra = prior_hits == 2

      qreq = get_qreq_pu(bus)
      lo, hi = q_limit_band(qmin_pu, qmax_pu, bus, q_hyst_pu)

      side = net.qLimitEvents[bus]
      voltage_rule = get_vm_pu !== nothing && get_vset_pu !== nothing && side in (:min, :max)
      if voltage_rule
        # at Qmax the machine would need LESS than Qmax to hold its
        # setpoint when the voltage sits above it (dQ/dV < 0), so it can
        # go back to PV; mirror image at Qmin
        vm = get_vm_pu(bus)
        vset = get_vset_pu(bus)
        ready = side == :max ? (vm > vset + v_hyst_pu) : (vm < vset - v_hyst_pu)
      else
        ready = (qreq > lo) && (qreq < hi)
      end

      if ready && (cooldown_iters > 0)
        last_it = lastQLimitIter(net, bus)
        if !isnothing(last_it) && (it - last_it) < cooldown_iters
          ready = false
        end
      end

      if ready && extra
        # collected, released together after the loop
        push!(extra_buses, bus)
      elseif ready
        release!(bus, side, voltage_rule, qreq, lo, hi)
      end
    end
    if !isempty(extra_buses) && !changed && !reenabled
      for bus in extra_buses
        side = net.qLimitEvents[bus]
        voltage_rule = get_vm_pu !== nothing && get_vset_pu !== nothing && side in (:min, :max)
        lo, hi = q_limit_band(qmin_pu, qmax_pu, bus, q_hyst_pu)
        release!(bus, side, voltage_rule, get_qreq_pu(bus), lo, hi)
      end
    end
  end

  if verbose > 0 && omitted_events > 0
    @printf(io, "PV->PQ event details: %d additional row(s) omitted; increase max_console_rows for full output.\n", omitted_events)
  end

  return changed, reenabled
end
