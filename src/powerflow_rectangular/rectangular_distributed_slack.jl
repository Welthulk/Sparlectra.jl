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
#
# This file is included inside module Sparlectra. Do not add a module wrapper here.
#
# file: src/powerflow_rectangular/rectangular_distributed_slack.jl
# purpose: distributed active-power slack for the rectangular power flow
# (issue #192). Participant discovery, per-island alpha normalization, and the
# solver-side state threaded through mismatch/Jacobian/Newton step. The REF
# bus keeps the angle reference; the island's P imbalance is absorbed by the
# participating generators via one scalar lambda_P:
#
#     rP_i = P_calc_i(V) - P_spec_i - alpha_i * lambda_P,   sum(alpha_i) = 1
#
# The system gains one unknown (lambda_P) and one equation (the REF bus's P
# residual, which is no longer free). The REF voltage stays fixed — the angle
# reference is untouched; only its P becomes an equation.

"""
Solver-side distributed-slack state for one island solve.

- `alpha`: per-bus participation, `length(alpha) == n`, `sum(alpha) == 1`
  (aggregated over the generators of a bus; zero for non-participants). Fixed
  for the whole solve — a PV→PQ switch during Q-limit handling does NOT remove
  a participant.
- `lambda`: the accepted `lambda_P` in p.u. (persists across the NR loop of
  one `runpf_rectangular!` call, starting at 0.0).
- `lambda_trial`: the value `mismatch_rectangular` actually reads. Step-search
  helpers (autodamp, trust region) evaluate trial states `V + a*δV`; the
  matching trial lambda is `lambda + a*δlambda`, set via
  [`_dslack_set_trial!`](@ref) before each trial mismatch evaluation.
- `gen_table`: generator-level rows for reporting (bus, prosumer index, raw
  weight, alpha share).
"""
mutable struct DistributedSlackState
  alpha::Vector{Float64}
  lambda::Float64
  lambda_trial::Float64
  delta_lambda::Float64
  mode::Symbol
  respect_p_limits::Bool
  gen_table::Vector{NamedTuple{(:bus, :gen_index, :weight, :alpha),Tuple{Int,Int,Float64,Float64}}}
  dropped::Int
end

DistributedSlackState(alpha::Vector{Float64}, mode::Symbol, respect_p_limits::Bool, gen_table, dropped::Int) = DistributedSlackState(alpha, 0.0, 0.0, 0.0, mode, respect_p_limits, gen_table, dropped)

"""Stage the trial lambda for a step of relative length `a` (see `DistributedSlackState`)."""
@inline function _dslack_set_trial!(ds::Union{Nothing,DistributedSlackState}, a::Float64)
  ds === nothing && return nothing
  ds.lambda_trial = ds.lambda + a * ds.delta_lambda
  return nothing
end

"""Accept the staged trial: the last staged `lambda_trial` becomes the iterate's lambda."""
@inline function _dslack_accept!(ds::Union{Nothing,DistributedSlackState})
  ds === nothing && return nothing
  ds.lambda = ds.lambda_trial
  return nothing
end

"""
    build_distributed_slack_state(net, bus_types;
                                  p_mode, fallback, weights,
                                  respect_p_limits, island_label) -> state | nothing

Discover participants and build the per-island alpha vector (Task 3 of #192).

Candidates are the generator-type prosumers at the island's REF or PV buses.
Raw weights per mode: `:pg_weighted` scheduled Pg, `:pmax_weighted` maxP,
`:headroom_weighted` `max(maxP − Pg, 0)`, `:imported`
`ProSumer.participationFactor` (MATPOWER `APF` / CGMES `normalPF`),
`:explicit` the config table (bus name or bus index as key). Invalid
candidates (missing data, weight ≤ 0, NaN/Inf) are dropped with a debug log
and counted. Surviving weights are normalized to `sum(alpha) = 1`, aggregated
per bus.

With no surviving participant, `fallback = :error` throws an `ArgumentError`
naming island and mode; `:ref_only` returns `nothing` with a warning — the
island then solves classically.
"""
function build_distributed_slack_state(
  net::Net,
  bus_types::Vector{Symbol};
  p_mode::Symbol,
  fallback::Symbol,
  weights::AbstractDict{String,Float64} = Dict{String,Float64}(),
  respect_p_limits::Bool = true,
  island_label::AbstractString = "",
)::Union{Nothing,DistributedSlackState}
  n = length(net.nodeVec)
  # explicit weights are keyed by bus name or by bus index written as string
  explicit_by_bus = Dict{Int,Float64}()
  if p_mode === :explicit
    for (k, w) in weights
      idx = get(net.busDict, k, nothing)
      if idx === nothing
        parsed = tryparse(Int, k)
        idx = (parsed !== nothing && 1 <= parsed <= n) ? parsed : nothing
      end
      if idx === nothing
        @debug "distributed slack: explicit weight key does not resolve to a bus" key = k
        continue
      end
      explicit_by_bus[idx] = get(explicit_by_bus, idx, 0.0) + w
    end
  end

  rows = NamedTuple{(:bus, :gen_index, :weight, :alpha),Tuple{Int,Int,Float64,Float64}}[]
  dropped = 0
  for (gen_index, ps) in enumerate(net.prosumpsVec)
    isGenerator(ps) || continue
    bus = getPosumerBusIndex(ps)
    (1 <= bus <= n) || continue
    # Candidates live at REF or PV buses — evaluated on the bus types at
    # discovery time. Fixed injections at PQ buses (Stage-0 HVDC converters,
    # boundary equivalents) are excluded naturally by the bus-type gate.
    bus_types[bus] in (:Slack, :PV) || continue
    w = if p_mode === :pg_weighted
      something(ps.pVal, 0.0)
    elseif p_mode === :pmax_weighted
      something(ps.maxP, 0.0)
    elseif p_mode === :headroom_weighted
      max(something(ps.maxP, 0.0) - something(ps.pVal, 0.0), 0.0)
    elseif p_mode === :imported
      something(ps.participationFactor, 0.0)
    elseif p_mode === :explicit
      get(explicit_by_bus, bus, 0.0)
    else
      throw(ArgumentError("distributed slack: unsupported p_mode $(p_mode)"))
    end
    if !isfinite(w) || w <= 0.0
      dropped += 1
      @debug "distributed slack: candidate dropped" bus gen_index mode = p_mode weight = w
      continue
    end
    push!(rows, (bus = bus, gen_index = gen_index, weight = w, alpha = 0.0))
  end

  if isempty(rows)
    label = isempty(island_label) ? "the island" : island_label
    if fallback === :error
      throw(ArgumentError("distributed slack: no valid participant on $(label) with p_mode=$(p_mode) ($(dropped) candidate(s) dropped). Set power_flow.distributed_slack.fallback=ref_only to fall back to the classical slack."))
    end
    @warn "Distributed slack: no valid participant — falling back to the classical slack for this island." island = label mode = p_mode dropped
    return nothing
  end

  total = sum(r.weight for r in rows)
  alpha = zeros(Float64, n)
  for r in rows
    alpha[r.bus] += r.weight / total
  end
  gen_table = [(bus = r.bus, gen_index = r.gen_index, weight = r.weight, alpha = r.weight / total) for r in rows]
  return DistributedSlackState(alpha, p_mode, respect_p_limits, gen_table, dropped)
end

"""
    _merge_distributed_slack_diagnostics(status_build, performance_profile, net, dslack, Sbase, verbose) -> status_build

Task 5 of #192: attach distributed-slack result metadata to the solver status
NamedTuple (flat `distributed_slack_*` keys, following the `wrong_branch_*`
naming of the existing metadata) and print a compact console summary at
`verbose > 0`. The full participant table goes to the debug log only — the
console stays at a few lines per the repository logging rule.

When `respect_p_limits` is set, each participant whose corrected output
`Pg + alpha_gen * lambda_P * Sbase` leaves `[minP, maxP]` produces one WARN;
the count is reported as `distributed_slack_p_limit_violations`. Stage 1
warns only — it does not clamp and re-solve.
"""
function _merge_distributed_slack_diagnostics(status_build, performance_profile, net::Net, dslack::Union{Nothing,DistributedSlackState}, Sbase::Float64, verbose::Int)
  active = dslack !== nothing
  lambda_pu = active ? dslack.lambda : 0.0
  lambda_mw = lambda_pu * Sbase
  participants = active ? length(dslack.gen_table) : 0
  alpha_sum = active ? sum(dslack.alpha) : 0.0
  dropped = active ? dslack.dropped : 0

  p_limit_violations = 0
  if active && dslack.respect_p_limits
    for r in dslack.gen_table
      ps = net.prosumpsVec[r.gen_index]
      pg = something(ps.pVal, 0.0)
      corrected = pg + r.alpha * lambda_mw
      pmax = ps.maxP
      pmin = ps.minP
      # Relative tolerance: an absolute 1e-6 MW fired ~20 warnings per RealGrid
      # run for pure epsilon overshoots of participants scheduled exactly at
      # Pg = maxP. A violation below 0.01 % of the limit (floor 1e-3 MW) is
      # numerical noise, not a reportable limit breach.
      tol_hi = pmax === nothing ? 0.0 : max(1e-3, 1e-4 * abs(pmax))
      tol_lo = pmin === nothing ? 0.0 : max(1e-3, 1e-4 * abs(pmin))
      hi = pmax !== nothing && corrected > pmax + tol_hi
      lo = pmin !== nothing && corrected < pmin - tol_lo
      if hi || lo
        p_limit_violations += 1
        @warn "Distributed slack: corrected participant output leaves its P limits (stage 1 warns only, no clamping)." bus = getCompName(net.nodeVec[r.bus].comp) gen_index = r.gen_index pg_mw = pg correction_mw = r.alpha * lambda_mw corrected_mw = corrected minP = pmin maxP = pmax
      end
    end
  end

  if active && verbose > 0
    # Top participants by alpha; the full table is debug-only.
    order = sortperm(dslack.gen_table; by = r -> -r.alpha)
    top = first(order, min(3, length(order)))
    println(stdout, "Distributed slack summary:")
    @printf(stdout, "  mode              = %s\n", String(dslack.mode))
    @printf(stdout, "  lambda_P          = %+.3f MW (%.6f p.u.)\n", lambda_mw, lambda_pu)
    @printf(stdout, "  participants      = %d (sum alpha = %.3f, dropped = %d)\n", participants, alpha_sum, dropped)
    for k in top
      r = dslack.gen_table[k]
      @printf(stdout, "  top participant   : %s  alpha=%.4f  dP=%+.3f MW\n", getCompName(net.nodeVec[r.bus].comp), r.alpha, r.alpha * lambda_mw)
    end
    p_limit_violations > 0 && @printf(stdout, "  P-limit violations = %d (see warnings)\n", p_limit_violations)
  end
  active && @debug "Distributed slack participant table" gen_table = dslack.gen_table lambda_pu lambda_mw

  if performance_profile isa AbstractDict
    performance_profile[:distributed_slack_active] = active
    performance_profile[:distributed_slack_lambda_p_mw] = lambda_mw
  end

  # Persist the participant table (bus name resolved here, dP in MW) so
  # reporting code such as printACPFlowResults can show per-bus participation
  # after the solve; the transient DistributedSlackState is gone by then.
  participation_row = NamedTuple{(:bus, :bus_idx, :gen_index, :alpha, :dp_mw, :pg_mw),Tuple{String,Int,Int,Float64,Float64,Float64}}
  participation = participation_row[]
  if active
    for r in dslack.gen_table
      ps = net.prosumpsVec[r.gen_index]
      push!(participation, (bus = getCompName(net.nodeVec[r.bus].comp), bus_idx = r.bus, gen_index = r.gen_index, alpha = r.alpha, dp_mw = r.alpha * lambda_mw, pg_mw = something(ps.pVal, 0.0)))
    end
  end

  extras = (
    distributed_slack_active = active,
    distributed_slack_mode = active ? dslack.mode : :none,
    distributed_slack_lambda_p_pu = lambda_pu,
    distributed_slack_lambda_p_mw = lambda_mw,
    distributed_slack_participants = participants,
    distributed_slack_alpha_sum = alpha_sum,
    distributed_slack_dropped = dropped,
    distributed_slack_p_limit_violations = p_limit_violations,
    distributed_slack_participation = participation,
  )
  return merge(status_build, (status = (; status_build.status..., extras...),))
end
