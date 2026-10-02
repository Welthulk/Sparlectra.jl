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

# file: src/powerflow_rectangular/rectangular_power_mode.jl
# purpose: power mode (0.30.1, `power_flow.power_mode`): the state a repeated
#          solve on one network keeps between solves, so the warm solve
#          neither rebuilds the Ybus nor repeats the symbolic analysis of
#          the sparse LU nor allocates its work arrays again; the choice
#          between KLU and UMFPACK (0.30.2, `power_flow.power_mode_lu`,
#          made once per network from KLU's symbolic flop estimate under
#          `auto`); the hooks through which the KLU package extension
#          supplies its factorization.

"""
    POWER_MODE_LU_KLU_MAX_FLOPS

The threshold of `power_flow.power_mode_lu = auto`: KLU's symbolic estimate
of the factorization flops (`klu_analyze`, BTF plus AMD, no numeric work)
of the first power-mode Jacobian of a network. Up to it power mode
factorizes with KLU, above it with UMFPACK.

The value is fitted on measured Jacobians (structural Jacobian of the
first solve, refactorization time of both backends): KLU was
faster on case2869pegase (estimate 1.95e6), case9241pegase (1.29e7),
mvlv29840 (1.76e6), case_ACTIVSg2000 (1.13e7), case13659pegase (1.52e7)
and the islands 2 and 3 of case_SyntheticUSA (2.71e7, 1.05e7); UMFPACK was
faster on case_ACTIVSg25k (1.26e8) and island 1 of case_SyntheticUSA
(5.60e8). 5e7 lies between the two groups. The estimate does not see
KLU's numeric pivoting, which is what makes KLU slow on the large cases
(5762 off-diagonal pivots on case_SyntheticUSA island 1); the threshold
uses the size of the problem as its stand-in.
"""
const POWER_MODE_LU_KLU_MAX_FLOPS = 5.0e7

"""
    PowerModeLuDecision

The sparse LU `power_flow.power_mode_lu = auto` chose for a network (or one
island of it): `choice` (`:klu` or `:umfpack`) and `klu_est_flops`, KLU's
symbolic flop estimate of the Jacobian the choice was made on, compared
with [`POWER_MODE_LU_KLU_MAX_FLOPS`](@ref).
"""
struct PowerModeLuDecision
  choice::Symbol
  klu_est_flops::Float64
end

# The `auto` decisions of one network: key 0 for the network solved as a
# whole, the smallest bus index of the parent network for an island solved
# on its own. Guarded by a lock: the island nets of one solve record into
# the memo of their parent network, and islands may run on their own tasks.
struct PowerModeLuMemo
  lock::ReentrantLock
  decisions::Dict{Int,PowerModeLuDecision}
end

PowerModeLuMemo() = PowerModeLuMemo(ReentrantLock(), Dict{Int,PowerModeLuDecision}())

_power_mode_lu_lookup(memo::PowerModeLuMemo, key::Int) = lock(() -> get(memo.decisions, key, nothing), memo.lock)
_power_mode_lu_record!(memo::PowerModeLuMemo, key::Int, decision::PowerModeLuDecision) = lock(() -> (memo.decisions[key] = decision), memo.lock)

"""
    PowerModeCache

What a power-mode solve keeps on the network (`net._power_cache`) for the
next solve: the Ybus with the fingerprint of the data it was built from,
the linear-solver context ([`PowerModeLuContext`](@ref): symbolic
analysis, factorization object, the Jacobian assembly buffers), the Newton
iteration workspace, and the `power_flow.power_mode_lu = auto` decisions of
the network (`lu_memo`, under `lu_key`; island nets of a solve record into
the memo of their parent network). The Ybus is reused only when the
fingerprint matches, the context re-analyses by itself when the Jacobian
pattern differs (bus types, bus count), so a changed topology or a changed
branch parameter never reuses stale state; the LU decision is kept across
such changes (once per network). `lu_source` is where the LU choice of the
last solve came from (`:decided`, `:earlier_decision`, `:configured`,
`:klu_not_loaded`).
"""
mutable struct PowerModeCache
  n::Int
  ybus::Any
  ybus_fingerprint::UInt64
  linear_ctx::Any
  workspace::Any
  bus_type_fingerprint::UInt64
  ybus_reuse_count::Int
  context_reuse_count::Int
  lu_memo::PowerModeLuMemo
  lu_key::Int
  lu_source::Symbol
end

PowerModeCache() = PowerModeCache(0, nothing, UInt64(0), nothing, nothing, UInt64(0), 0, 0, PowerModeLuMemo(), 0, :configured)

# The KLU package extension (ext/SparlectraKLUExt.jl) registers its context
# constructor here when `using KLU` loads it; nothing means no KLU in the
# session. Power mode reads it as the availability signal only (it
# factorizes through the `_klu_*` functions below); the estimator's
# `state_estimation.linear_solver = klu` builds its LU context from it.
const _POWER_MODE_LINEAR_CONTEXT = Ref{Any}(nothing)

_klu_available() = _POWER_MODE_LINEAR_CONTEXT[] !== nothing

# Defined by the KLU package extension; the base package calls them only
# after `_klu_available()` said it is loaded.
# _klu_analyze(J) -> (object, est_flops): symbolic analysis only
# _klu_factor_analyzed!(object) -> factorization on that analysis
# _klu_factorize(J), _klu_refactorize!(F, J): analysis plus numeric, and the
#   numeric refactorization on the stored analysis
function _klu_analyze end
function _klu_factor_analyzed! end
function _klu_factorize end
function _klu_refactorize! end

"""
    power_mode_linear_solver_backend() -> Symbol

`:klu` when the KLU extension is loaded (`using KLU` in the session), so
that `power_flow.power_mode_lu` can choose KLU; `:umfpack_reuse` otherwise,
the only sparse LU a power-mode solve then has. Which backend a solve
used is in its status (`power_mode_lu_choice`).
"""
power_mode_linear_solver_backend() = _klu_available() ? :klu : :umfpack_reuse

"""
    PowerModeLuContext

The linear-solver context of a power-mode solve: the reuse scheme of
[`UmfpackReuseNewtonContext`](@ref) (one analysis per Jacobian pattern,
a numeric refactorization per Newton step) with the backend chosen by
`power_flow.power_mode_lu`. `backend` is `:klu` or `:umfpack` and says what
`fact` is. With `decision_pending` set (`auto`, no decision for this
network yet), the first full factorization decides
([`PowerModeLuDecision`](@ref)) and records the decision in `memo` under
`key`. One context per network (per island net on the island path), never
shared across tasks.
"""
mutable struct PowerModeLuContext <: AbstractNewtonSolverContext
  fact::Any
  nvar::Int
  colptr::Vector{Int64}
  rowval::Vector{Int64}
  analyze_count::Int
  refactor_count::Int
  fallback_count::Int
  rhs::Vector{Float64}
  sol::Vector{Float64}
  assembly::RectangularJacobianAssembly
  backend::Symbol
  decision_pending::Bool
  key::Int
  memo::Union{Nothing,PowerModeLuMemo}
  # the decision this solve made; nothing when the choice came from the
  # memo or the configuration (reset at the start of every solve)
  decided::Union{Nothing,PowerModeLuDecision}
end

PowerModeLuContext() = PowerModeLuContext(nothing, 0, Int64[], Int64[], 0, 0, 0, Float64[], Float64[], RectangularJacobianAssembly(), :umfpack, false, 0, nothing, nothing)

function _newton_full_factorization(ctx::PowerModeLuContext, J::SparseMatrixCSC{Float64,Int64})
  ctx.decision_pending && return _power_mode_lu_decide!(ctx, J)
  return ctx.backend === :klu ? _klu_factorize(J) : lu(J)
end
_newton_refactorize!(ctx::PowerModeLuContext, J::SparseMatrixCSC{Float64,Int64}) = ctx.backend === :klu ? _klu_refactorize!(ctx.fact, J) : lu!(ctx.fact, J)
_newton_context_backend(ctx::PowerModeLuContext) = ctx.backend === :klu ? :klu : :umfpack_reuse

# The choice of `auto` from KLU's symbolic analysis of J alone (no numeric
# factorization): KLU up to POWER_MODE_LU_KLU_MAX_FLOPS estimated flops,
# UMFPACK above. KLU's analysis is kept when KLU wins, so the decision
# costs one extra symbolic analysis only when UMFPACK wins.
_power_mode_lu_choice(est_flops::Float64) = est_flops <= POWER_MODE_LU_KLU_MAX_FLOPS ? :klu : :umfpack

function _power_mode_lu_decide!(ctx::PowerModeLuContext, J::SparseMatrixCSC{Float64,Int64})
  K, est_flops = _klu_analyze(J)
  decision = PowerModeLuDecision(_power_mode_lu_choice(est_flops), est_flops)
  ctx.backend = decision.choice
  ctx.decision_pending = false
  ctx.decided = decision
  ctx.memo === nothing || _power_mode_lu_record!(ctx.memo, ctx.key, decision)
  return decision.choice === :klu ? _klu_factor_analyzed!(K) : lu(J)
end

"""
    _power_mode_lu_select!(cache, power_mode_lu) -> PowerModeLuContext

Sets up the linear-solver context of a power-mode solve for
`power_flow.power_mode_lu`: `:umfpack` and `:klu` fix the backend (`:klu`
without the loaded extension falls back to UMFPACK with a warning);
`:auto` without the extension is UMFPACK; with it, the decision the memo
holds for this network (`cache.lu_key`, else the whole network, key 0),
or a pending decision at the first factorization when there is none. The
context keeps its analysis while the backend stays; a change of backend
drops the factorization, so the next Newton step factorizes afresh.
"""
function _power_mode_lu_select!(cache::PowerModeCache, power_mode_lu::Symbol)
  ctx = cache.linear_ctx
  if !(ctx isa PowerModeLuContext)
    ctx = PowerModeLuContext()
    cache.linear_ctx = ctx
  end
  ctx.memo = cache.lu_memo
  ctx.key = cache.lu_key
  ctx.decided = nothing
  ctx.decision_pending = false
  source = :configured
  backend = :umfpack
  if power_mode_lu === :klu
    if _klu_available()
      backend = :klu
    else
      @warn "power_flow.power_mode_lu = klu needs the KLU package extension (SparlectraKLUExt, loaded by `using KLU`), which is not loaded; power mode factorizes with UMFPACK" maxlog = 1
      source = :klu_not_loaded
    end
  elseif power_mode_lu === :auto
    if !_klu_available()
      source = :klu_not_loaded
    else
      # an island that has no decision of its own takes the one of the
      # network solved as a whole (an outage that splits a network keeps
      # the network's choice)
      known = _power_mode_lu_lookup(cache.lu_memo, cache.lu_key)
      known === nothing && cache.lu_key != 0 && (known = _power_mode_lu_lookup(cache.lu_memo, 0))
      if known === nothing
        ctx.decision_pending = true
        source = :decided
        # the decision factorizes afresh
        ctx.fact = nothing
        backend = ctx.backend
      else
        backend = known.choice
        source = :earlier_decision
      end
    end
  elseif power_mode_lu !== :umfpack
    throw(ArgumentError("power_mode_lu must be one of $(POWER_MODE_LU_VALUES), got :$(power_mode_lu)"))
  end
  backend !== ctx.backend && (ctx.fact = nothing)
  ctx.backend = backend
  cache.lu_source = source
  return ctx
end

# The power-mode LU of one solve for its status, the run log and the
# metadata: what was asked for, what was used, where the choice came from,
# KLU's symbolic flop estimate it was made on (NaN for a configured choice).
function _power_mode_lu_status(cache::PowerModeCache, ctx::PowerModeLuContext, power_mode_lu::Symbol)
  source = cache.lu_source
  # a pending decision that never ran: the solve left before its first
  # factorization (converged at the start, or stopped)
  source === :decided && ctx.decided === nothing && (source = :not_decided)
  decision = ctx.decided
  if decision === nothing && source === :earlier_decision
    decision = _power_mode_lu_lookup(cache.lu_memo, cache.lu_key)
    decision === nothing && (decision = _power_mode_lu_lookup(cache.lu_memo, 0))
  end
  est = decision === nothing ? NaN : decision.klu_est_flops
  choice = ctx.backend
  return (
    power_mode_lu = power_mode_lu,
    power_mode_lu_choice = choice,
    power_mode_lu_source = source,
    power_mode_lu_klu_est_flops = est,
    power_mode_lu_threshold_flops = POWER_MODE_LU_KLU_MAX_FLOPS,
    power_mode_lu_line = _power_mode_lu_line(choice, source, est, power_mode_lu),
  )
end

# One line with the choice and the measure it was made on (run.log result
# header, the run metadata).
function _power_mode_lu_line(choice::Symbol, source::Symbol, est::Float64, requested::Symbol)
  if source === :decided || source === :earlier_decision
    side = choice === :klu ? "at most" : "above"
    when = source === :decided ? "decided in this solve" : "decided in an earlier solve of this network"
    return @sprintf("%s (KLU symbolic estimate %.3g flops, %s the threshold %.0e; %s)", choice, est, side, POWER_MODE_LU_KLU_MAX_FLOPS, when)
  end
  source === :klu_not_loaded && return "umfpack (KLU extension not loaded; power_mode_lu = $(requested))"
  source === :not_decided && return "$(choice) (no factorization in this solve, not decided)"
  return "$(choice) (power_mode_lu = $(requested))"
end

# Fingerprint of everything createYBUS reads: the stamped admittances of
# every branch with its terminals and switching state, the shunt
# admittances, the isolated buses, the base. Cheaper than the assembly by
# two orders of magnitude (arithmetic per branch against sparse insertion),
# and a tap or a switching change moves it, so the cached Ybus is never
# stale.
function _ybus_fingerprint(net::Net)::UInt64
  h = hash(length(net.nodeVec), hash(net.baseMVA, UInt64(0x5b1ec7a0)))
  h = hash(net.isoNodes, h)
  for branch in net.branchVec
    h = hash((branch.fromBus, branch.toBus, branch.status, branch.from_status, branch.to_status), h)
    h = hash(calcAdmittance(branch, branch.comp.cVN, net.baseMVA), h)
  end
  for sh in net.shuntVec
    h = hash((sh.busIdx, sh.status, sh.model, sh.y_pu_shunt), h)
  end
  return h
end

function _power_cache!(net::Net)::PowerModeCache
  cache = net._power_cache
  if !(cache isa PowerModeCache)
    cache = PowerModeCache()
    net._power_cache = cache
  end
  return cache
end

# The working nets the solver builds per solve (island nets, the merged
# copy of an active-link merge) start with a copy of the cache (deepcopy,
# see below); they record their `auto` decision into the memo of the
# network the caller holds, under `key` (0 for the merged copy, the
# smallest parent bus index for an island), so the decision is made once
# per network and island and not in every solve. Called by the
# orchestrator before any island task starts.
function _share_power_mode_lu_memo!(working::Net, memo::PowerModeLuMemo, key::Int)
  cache = _power_cache!(working)
  cache.lu_memo = memo
  cache.lu_key = key
  return nothing
end

"""
    reset_power_mode!(net)

Drops the power-mode state of `net` (Ybus, factorization, workspace, the
LU decisions). The next power-mode solve rebuilds it and decides the LU
again. Not needed after a topology or parameter change (the fingerprint
and the pattern guard cover that); for callers that want the memory back.
"""
function reset_power_mode!(net::Net)
  net._power_cache = nothing
  return nothing
end

# A copied network starts without the numeric power-mode state: the
# factorization objects wrap native memory (UMFPACK, KLU) that must not be
# shared or double-freed, and the Ybus fingerprint is recomputed on the
# first solve anyway. The `auto` LU decisions are plain values and are
# copied: a copy is the same network (the workers of N-1 and of the
# scenario engine are copies of the solved template and keep its choice).
function Base.deepcopy_internal(cache::PowerModeCache, ::IdDict)
  copied = PowerModeCache()
  lock(cache.lu_memo.lock) do
    merge!(copied.lu_memo.decisions, cache.lu_memo.decisions)
  end
  copied.lu_key = cache.lu_key
  return copied
end
