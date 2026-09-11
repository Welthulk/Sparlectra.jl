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
# Date: 1.1.2026
# file: src/stateestimation/state_estimation.jl
# purpose: classical WLS state estimation (runse!) with observability and
#          redundancy analysis and an optional PMU reference-angle offset
#          state

using Printf
"""
    SEResult

Result container for the first classical WLS state-estimation run.

`vaRefOffsetDeg` is the estimated common PMU reference-angle offset α in
degrees (the slack-bus angle expressed in the PMU time base). It is `nothing`
when no offset state was part of the estimation, i.e. when there are no
active `VaMeas` measurements or `pmu_ref_offset = :off`.

`shuntEstimates` (SE phase 2, case A) carries one row per shunt released via
`setShuntEstimation!`: `(busIdx, busName, B_model, B_est, delta, frozen)`,
susceptances in pu. Frozen rows (no direct measurement, or B column not
observable) keep `B_est = B_model`. `nothing` when no shunt is released.
Each estimated B consumes one degree of freedom (it counts as a state in
ν = m − n).

`robustRows` (SE phase 4, `robust = true`) lists every measurement that left
stage 0 of the two-stage R modification in the FINAL iteration:
`(measurement_index, id, stage, t, sigma_factor)` with `t = |z - h|/sigma`
(original sigma) and `sigma_factor = sigma_mod/sigma`. `nothing` when the
robust mode is off. The modification changes only the solve weights;
`objectiveJ`, the residuals, and every diagnostic statistic stay on the
original sigmas.
"""
struct SEResult
  voltages::Vector{ComplexF64}
  converged::Bool
  iterations::Int
  residualNorm::Float64
  residuals::Vector{Float64}
  objectiveJ::Float64
  dof::Int # degrees of freedom ν=m−n
  jWithin3Sigma::Bool
  vaRefOffsetDeg::Union{Nothing,Float64}
  shuntEstimates::Union{Nothing,Vector{NamedTuple}}
  robustRows::Union{Nothing,Vector{NamedTuple}}
  # island-wise SE: one row per AC island (estimated or skipped); `nothing`
  # on a single-island net. On a multi-island net `voltages` is indexed by
  # the ORIGINAL bus numbering (unestimated islands keep their start
  # values), `objectiveJ`/`dof` are summed over the estimated islands, and
  # `vaRefOffsetDeg` is set only when exactly one island estimated it (the
  # per-island values live in these rows).
  islands::Union{Nothing,Vector{NamedTuple}}
  # tap estimation: one row per released transformer with the
  # continuous (electrical) and rounded (fixed/mechanical) step per
  # regulator; `nothing` without releases. After the estimation the taps
  # are FIXED to the mechanical step and one more run is solved in which
  # the tap is no state any more: `objectiveJ`/`dof`/`voltages` above are
  # that FIXED final run, `tapFixation` carries J and dof before versus
  # after the fixation.
  tapEstimates::Union{Nothing,Vector{NamedTuple}}
  tapFixation::Union{Nothing,NamedTuple}
  # topology validation: the ADVISORY stage-1 precheck findings
  # (status contradictions, dead closed branches, link voltage mismatches,
  # node-balance violations); nothing when the precheck is off or clean.
  # Findings never block or mutate anything.
  topologyFindings::Union{Nothing,Vector{NamedTuple}}
  # J_active: `(j, dof, suppressed)` when replacement-mode suppression
  # removed rows from the STATE (frozen suppression-sigma weights). j/dof
  # exclude those rows, answering "how well does the model fit the data
  # the estimator actually trusted". The honest objectiveJ above keeps
  # them at original sigmas as the alarm signal, and the band test stays
  # on objectiveJ. nothing when no row was suppressed.
  activeObjective::Union{Nothing,NamedTuple}
end

## Two-stage robust weight modification (SE phase 4): per
## measurement `t = |residual| / sigma` with the ORIGINAL sigma.
## - stage 0, t <= k1: weight unchanged
## - stage 1 (tangential), k1 < t <= k2: sigma_mod = sigma (2t/k1 - 1)
## - stage 2 (constant), t > k2: sigma_mod = |residual| / k1 (the gradient
##   contribution is effectively suppressed)
## The knees default to the historical 3/6 and are bitwise identical to the
## pre-0.10 fixed formulas at those values. Stages are recomputed every
## iteration, so a measurement can recover. Virtual measurements
## (sigma <= 1e-6, including ZI pseudo-measurements) are exempt; SHDERIV and
## LINKAGG rows carry real sigmas and remain subject to modification.
## Returns (stage, sigma_mod).
@inline function _robust_stage(absr::Float64, σ::Float64, k1::Float64 = 3.0, k2::Float64 = 6.0)
  if σ <= 1e-6
    return (0, σ)
  end
  t = absr / σ
  if t <= k1
    return (0, σ)
  end
  if t <= k2
    return (1, σ * (2.0 * t / k1 - 1.0))
  end
  return (2, absr / k1)
end

## Effective robust mode of a config: the legacy Bool `robust` selects the
## staged modification only while robust_mode itself is :off, so every
## pre-0.10 call site (robust = true) keeps its exact behavior.
@inline _se_effective_robust_mode(cfg::StateEstimationConfig)::Symbol = (cfg.robust_mode === :off && cfg.robust) ? :staged : cfg.robust_mode

"""
    numeric_rank(A; tol=nothing) -> Int

Compute numerical matrix rank using singular values.
"""
## Sparse/dense crossover, measured: the dense SVD/pinv machinery is
## sub-second up to roughly n = 2000 states and cubic beyond, so 2000 is
## where the dense conveniences stop paying for themselves; the sparse
## routes (SPQR rank, Takahashi diagnostics) take over above.
const _SE_DENSE_LINALG_MAX_N = 2_000
## budget for the per-row criticality classification (m rank tests): the
## step-0 baseline measured 7.6 s at m*n = 210k and non-termination at
## m*n = 11M, so 300k keeps the check within roughly a minute
const _SE_ROW_CRITICALITY_MAX_MN = 300_000
## FD column coloring engages only where it is the lever: below this
## state count the per-column assembly costs milliseconds anyway, and the
## small test networks keep their exact warning contracts (the pattern of
## a tiny system is also the least stable across iterations)
const _SE_FD_COLORING_MIN_STATES = 1_000
## the full residual-covariance matrix (K-matrix report) is m×m dense by
## definition; above this row count the report must be refused by name
## instead of allocating gigabytes
const _SE_FULL_OMEGA_MAX_M = 20_000

## Localizability bound for the replacement suppression. `wii` is the share
## of a row's own error that reaches its residual; below the literature
## guideline 0.3 (Abur/Gomez Exposito, the same bound the bad-data report
## uses for its `localizable` flag) the row is nearly critical: its residual
## is structurally small, so its NORMALIZED residual is inflated and says
## little about the row itself. Suppressing such a row throws away the
## little information it carries and pushes the error onto the neighbour
## that shared its redundancy. Measured on warmup_casePST: suppressing
## Qinj_5 (wii 0.05, raw residual 1.4 sigma, rn 5.8) drove its partner to
## rn 58 and had that healthy row eliminated (task_se_bad_data_v0100).
const _SE_SUPPRESSION_MIN_WII = 0.3

function numeric_rank(A::AbstractMatrix{<:Real}; tol = nothing)
  s = svdvals(Matrix(A))
  isempty(s) && return 0
  if tol === nothing
    tol = eps(eltype(s)) * max(size(A)...) * maximum(s)
  end
  return count(>(tol), s)
end

## Sparse rank via SPQR (task_se_sparse): below the dense threshold the
## exact dense SVD path keeps its historical semantics bit for bit; above
## it the rank comes from the sparse QR factorization, whose column-norm
## tolerance is fed from the SAME effective tolerance the callers derive
## (rank_tol_factor * jacEps * sigma_max, or the SVD default formula).
function numeric_rank(A::SparseMatrixCSC{Float64}; tol = nothing)
  m, n = size(A)
  (m == 0 || n == 0) && return 0
  max(m, n) <= _SE_DENSE_LINALG_MAX_N && return numeric_rank(Matrix(A); tol = tol)
  effTol = tol === nothing ? eps(Float64) * max(m, n) * _sigma_max_estimate(A) : Float64(tol)
  F = qr(A; tol = effTol)
  return rank(F)
end

## Largest singular value dispatcher: exact dense svdvals below the
## threshold (bit-identical to the historical behavior for every existing
## caller), power-iteration estimate above it.
_sigma_max(H::AbstractMatrix{<:Real})::Float64 = maximum(svdvals(Matrix(H)))
_sigma_max(H::SparseMatrixCSC{Float64})::Float64 = max(size(H)...) <= _SE_DENSE_LINALG_MAX_N ? maximum(svdvals(Matrix(H))) : _sigma_max_estimate(H)

## Largest singular value of a sparse matrix without densifying: power
## iteration on AᵀA with a fixed deterministic start vector. 30 iterations
## put the estimate well inside one percent, which is plenty for a rank
## TOLERANCE that carries a factor-10 safety margin anyway; the dense
## paths below the threshold keep using exact svdvals.
function _sigma_max_estimate(A::SparseMatrixCSC{Float64})::Float64
  n = size(A, 2)
  n == 0 && return 0.0
  v = fill(1.0 / sqrt(n), n)
  s = 0.0
  for _ = 1:30
    u = A * v
    v2 = A' * u
    nv = norm(v2)
    nv == 0.0 && return 0.0
    v = v2 ./ nv
    s = sqrt(nv)
  end
  return s
end

## Internal helper: find slack bus index in the current network.
@inline function _find_slack_idx(net::Net)
  for (i, node) in enumerate(net.nodeVec)
    if getNodeType(node) == Slack
      return i
    end
  end
  error("No slack bus found for state estimation")
end

## Internal helper: true when the active set contains PMU angle measurements.
@inline function _has_active_va(measurements::Vector{Measurement})
  # PMU voltage AND current angles share the same time base: either kind
  # activates the common reference-angle offset state alpha
  return any(m -> m.active && (m.typ == VaMeas || m.typ == IaMeas), measurements)
end

## exact residual wrap for angle-type rows (VaMeas, IaMeas): a measured and
## a predicted angle can sit on opposite sides of the +-180 degree seam.
## Conditional shift only, so in-range residuals stay BITWISE untouched
## (mod-based wrapping would cost precision on tiny residuals).
@inline function _wrap_angle_residuals!(r::Vector{Float64}, measurements::Vector{Measurement})
  @inbounds for i in eachindex(measurements)
    typ = measurements[i].typ
    (typ == VaMeas || typ == IaMeas) || continue
    if r[i] > 180.0
      r[i] -= 360.0
    elseif r[i] < -180.0
      r[i] += 360.0
    end
  end
  return r
end

## Internal helper: decide whether the PMU reference-offset state α is carried.
## :auto adds α as soon as active VaMeas rows exist; :off treats PMU angles as
## already slack-referenced.
@inline function _va_offset_active(measurements::Vector{Measurement}, mode::Symbol)
  mode == :off && return false
  mode == :auto || error("pmu_ref_offset mode must be :auto or :off, got :$(mode)")
  return _has_active_va(measurements)
end

## Internal helper: map estimator state vector x -> complex bus voltages V.
## x = [θ(non-slack); Vm(all buses)] with θ_slack fixed to 0.
## A trailing PMU reference-offset state α (if present) is ignored here and
## extracted separately via _va_offset_from_state.
function _state_to_voltage(x::Vector{Float64}, slackIdx::Int, nbus::Int)
  nθ = nbus - 1
  θ = zeros(Float64, nbus)
  vm = Vector{Float64}(undef, nbus)

  for i = 1:nbus
    vm[i] = x[nθ+i]
  end

  p = 0
  for i = 1:nbus
    if i == slackIdx
      θ[i] = 0.0
    else
      p += 1
      θ[i] = x[p]
    end
  end

  return vm .* cis.(θ)
end

## Internal helper: extract the PMU reference-offset state α (radians) from x.
@inline function _va_offset_from_state(x::Vector{Float64}, withVaOffset::Bool)
  return withVaOffset ? x[end] : 0.0
end

## Internal helper: build initial state vector for WLS iterations.
## flatstart=true uses θ=0, Vm=1 except slack magnitude initialized from net.
## State layout (SE phase 2): x = [θ(non-slack); Vm(all); B_1..B_k; α].
## shuntBInit carries the initial susceptances of the k released shunts
## (always the model values, never flat). withVaOffset=true appends the PMU
## reference-offset state α (radians, initialized to 0) as the LAST entry,
## so _va_offset_from_state stays x[end].
function _initial_state_vector(net::Net, slackIdx::Int; flatstart::Bool = true, withVaOffset::Bool = false, shuntBInit::Vector{Float64} = Float64[], tapRInit::Vector{Float64} = Float64[])
  nbus = length(net.nodeVec)
  nθ = nbus - 1
  k = length(shuntBInit)
  nt = length(tapRInit)
  x = zeros(Float64, nθ + nbus + k + nt + (withVaOffset ? 1 : 0))

  p = 0
  for i = 1:nbus
    if i != slackIdx
      p += 1
      x[p] = flatstart ? 0.0 : deg2rad(something(net.nodeVec[i]._va_deg, 0.0))
    end
  end

  for i = 1:nbus
    vm0 = something(net.nodeVec[i]._vm_pu, 1.0)
    x[nθ+i] = flatstart ? 1.0 : vm0
  end

  if flatstart
    x[nθ+slackIdx] = something(net.nodeVec[slackIdx]._vm_pu, 1.0)
  end

  for j = 1:k
    x[nθ+nbus+j] = shuntBInit[j]
  end

  # released transformer tap states: after the shunt B states,
  # before the trailing PMU offset
  for j = 1:nt
    x[nθ+nbus+k+j] = tapRInit[j]
  end

  return x
end

## Internal mapping for released shunt states (SE phase 2, case A):
## - shunt_idxs: indices into net.shuntVec, ascending
## - bus_idxs: the shunts' bus indices
## - col0: state-vector column of B_1 (B_j lives at col0 + j - 1)
## The B states sit between the Vm block and the trailing α, so
## _state_to_voltage (prefix) and _va_offset_from_state (x[end]) stay valid.
struct ShuntStateMap
  shunt_idxs::Vector{Int}
  bus_idxs::Vector{Int}
  col0::Int
end

## Extract the per-bus effective susceptances (busIdx -> B in pu) from the
## state vector for the measurement predictions.
function _shunt_b_override(x::Vector{Float64}, map::ShuntStateMap)
  d = Dict{Int,Float64}()
  for j in eachindex(map.shunt_idxs)
    d[map.bus_idxs[j]] = x[map.col0+j-1]
  end
  return d
end

## Internal helper: return only currently active measurements.
function _active_measurements(measurements::Vector{Measurement})
  return [m for m in measurements if m.active]
end

## Internal helper: return active measurements and their original indices.
function _active_measurements_with_indices(measurements::Vector{Measurement})
  idx = Int[]
  active = Measurement[]
  for (k, m) in enumerate(measurements)
    if m.active
      push!(idx, k)
      push!(active, m)
    end
  end
  return active, idx
end

@inline function _set_measurement_active(m::Measurement, active::Bool)
  return Measurement(typ = m.typ, value = m.value, sigma = m.sigma, active = active, busIdx = m.busIdx, branchIdx = m.branchIdx, direction = m.direction, id = m.id, linkIdx = m.linkIdx)
end

## SE working context after link contraction (SE phase 3).
##
## The estimator runs on `snet`: the net with every closed-link cluster
## contracted onto its representative and densely renumbered (no isolated
## buses), built from the PF contraction (`_merged_pf_net`) plus the island
## subnet builder (`_prepare_island_net`). `buses[new] = original bus`,
## `busmap[original rep] = new`, `reps` is the original-space cluster map.
## `meas`/`origIdx` are the remapped measurement rows and their positions in
## the caller's measurement vector (0 marks a LINKAGG aggregate row).
## `shuntOrig[snet shunt pos] = net.shuntVec position` for write-back.
## `notes` records everything that was excluded or aggregated.
struct SEWorkingNet
  snet::Net
  meas::Vector{Measurement}
  origIdx::Vector{Int}
  reps::Vector{Int}
  buses::Vector{Int}
  busmap::Dict{Int,Int}
  shuntOrig::Vector{Int}
  has_merges::Bool
  notes::Vector{NamedTuple}
end

## Build the SE working context: contract closed-link clusters, remap the
## active measurements per the phase-3 design rules, and always strip
## link-referenced flow measurements (allocation inputs, never WLS rows).
function _se_prepare(net::Net, measurements::Vector{Measurement})
  activeMeas, activeIdx = _active_measurements_with_indices(measurements)
  notes = NamedTuple[]
  # collected instead of warned per row: a station with a closed coupler can
  # collapse dozens of measured branches, and dozens of identical lines hide
  # everything else in a run log
  collapsed_ids = String[]

  # link flow measurements never enter the WLS: the link is not in the Ybus,
  # they constrain only the allocation split (calcLinkFlowsSE!)
  keep = Int[]
  for k in eachindex(activeMeas)
    if activeMeas[k].linkIdx !== nothing
      push!(notes, (reason = :link_measurement, id = activeMeas[k].id, detail = "allocation input, excluded from the WLS"))
      @debug "SE: link measurement $(activeMeas[k].id) excluded from the WLS (allocation input)"
    else
      push!(keep, k)
    end
  end
  activeMeas = activeMeas[keep]
  activeIdx = activeIdx[keep]

  wnet, reps, has_merges = _merged_pf_net(net)
  # fast path only for a pristine net: link merges AND isolated buses both
  # force the dense renumbering below, otherwise isolated buses would enter
  # the state vector and leave the FD Jacobian rank-deficient (their columns
  # have no measurement); found via CGMES imports, which routinely carry
  # isolated nodes
  if !has_merges && isempty(net.isoNodes)
    n = length(net.nodeVec)
    busmap = Dict{Int,Int}(i => i for i = 1:n)
    return SEWorkingNet(net, activeMeas, activeIdx, reps, collect(1:n), busmap, collect(eachindex(net.shuntVec)), false, notes)
  end

  # dense renumbering: reuse the island subnet builder on the contracted net
  # (the neutralized cluster members are exactly wnet.isoNodes)
  isoset = Set(wnet.isoNodes)
  buses = [i for i in eachindex(wnet.nodeVec) if !(i in isoset)]
  busset = Set(buses)
  branches = Int[]
  for (k, br) in enumerate(wnet.branchVec)
    st = _branch_terminal_state(br)
    st == :open && continue
    closedFrom = st != :open_from
    closedTo = st != :open_to
    (closedFrom && !(Int(br.fromBus) in busset)) && continue
    (closedTo && !(Int(br.toBus) in busset)) && continue
    push!(branches, k)
  end
  n_ref = count(b -> getNodeType(wnet.nodeVec[b]) == Slack, buses)
  row = (buses = buses, branches = branches, n_ref = n_ref, chosen_ref_bus = 0)
  snet = _prepare_island_net(wnet, row)
  # the contracted net carries no links any more: closed ones are fused into
  # their representative, open ones are real separations. Emptying the copy
  # also makes a nested _se_prepare (diagnostics reruns) a no-op fast path.
  empty!(snet.linkVec)
  busmap = Dict{Int,Int}(old => new for (new, old) in enumerate(buses))
  branchmap = Dict{Int,Int}(old => new for (new, old) in enumerate(branches))

  # the island builder replaces shuntVec but not shuntDict: rebuild it, and
  # keep the original shunt positions for the updateShunts write-back (the
  # subnet keeps net.shuntVec order, filtered to surviving buses)
  empty!(snet.shuntDict)
  for (k, sh) in enumerate(snet.shuntVec)
    snet.shuntDict[sh.busIdx] = k
  end
  shuntOrig = [k for (k, sh) in enumerate(net.shuntVec) if reps[Int(sh.busIdx)] in busset]
  length(shuntOrig) == length(snet.shuntVec) || error("SE link contraction: shunt bookkeeping out of sync (this is a bug)")

  # measurement remapping onto the contracted, renumbered net
  members_by_rep = Dict{Int,Vector{Int}}()
  for b in eachindex(reps)
    push!(get!(members_by_rep, reps[b], Int[]), b)
  end
  name_by_idx = _bus_name_by_idx(net)
  busname(b) = get(name_by_idx, b, string(b))

  remapped = Measurement[]
  origIdx = Int[]
  # cluster injection bookkeeping: (rep, typ) -> measurement positions
  injPositions = Dict{Tuple{Int,MeasurementType},Vector{Int}}()

  for k in eachindex(activeMeas)
    m = activeMeas[k]
    if m.busIdx !== nothing && (m.typ == PinjMeas || m.typ == QinjMeas)
      rep = reps[m.busIdx]
      if length(members_by_rep[rep]) > 1
        push!(get!(injPositions, (rep, m.typ), Int[]), k)
        continue   # handled cluster-wise below
      end
    end
    if m.typ in (PflowMeas, QflowMeas, ImagMeas) && m.branchIdx !== nothing
      newBr = get(branchmap, m.branchIdx, 0)
      if newBr == 0
        # branch collapsed inside a cluster (shorted by the closed coupler):
        # its flow is not defined by the fused state
        push!(notes, (reason = :collapsed_branch, id = m.id, detail = "branch $(m.branchIdx) collapsed by link contraction"))
        push!(collapsed_ids, String(m.id))
        continue
      end
      push!(remapped, Measurement(typ = m.typ, value = m.value, sigma = m.sigma, active = true, branchIdx = newBr, direction = m.direction, id = m.id))
      push!(origIdx, activeIdx[k])
      continue
    end
    if m.busIdx !== nothing
      # Vm/Va (voltage equality makes the remap exact), single-bus Pinj/Qinj,
      # and the device-referenced ShuntQMeas / bus ImagMeas. A bus outside
      # the working state (isolated) cannot carry a measurement row.
      if !haskey(busmap, reps[m.busIdx])
        push!(notes, (reason = :isolated_bus, id = m.id, detail = "bus $(m.busIdx) is isolated"))
        @warn "SE: measurement $(m.id) excluded (its bus is isolated)"
        continue
      end
      newBus = busmap[reps[m.busIdx]]
      push!(remapped, Measurement(typ = m.typ, value = m.value, sigma = m.sigma, active = true, busIdx = newBus, branchIdx = m.branchIdx, direction = m.direction, id = m.id))
      push!(origIdx, activeIdx[k])
      continue
    end
    push!(remapped, m)
    push!(origIdx, activeIdx[k])
  end

  # cluster injections: aggregable only as a whole; a partial sum
  # would bias the fused balance, so partial coverage excludes all of them
  for ((rep, typ), positions) in sort(collect(injPositions); by = x -> (x[1][1], Int(x[1][2])))
    members = members_by_rep[rep]
    measuredBuses = Set(activeMeas[k].busIdx for k in positions)
    if all(b -> b in measuredBuses, members)
      vals = [activeMeas[k].value for k in positions]
      sigs = [activeMeas[k].sigma for k in positions]
      ids = [activeMeas[k].id for k in positions]
      aggId = "LINKAGG_$(typ == PinjMeas ? "Pinj" : "Qinj")_rep_$(rep)"
      push!(remapped, Measurement(typ = typ, value = sum(vals), sigma = sqrt(sum(sigs .^ 2)), active = true, busIdx = busmap[rep], id = aggId))
      push!(origIdx, 0)
      push!(notes, (reason = :aggregated_cluster_injection, id = aggId, detail = "members: " * join(ids, ", ")))
    else
      for k in positions
        push!(notes, (reason = :partial_cluster_injection, id = activeMeas[k].id, detail = "cluster at $(busname(rep)) not fully measured"))
      end
      @warn "SE: $(length(positions)) $(typ) measurement(s) on the link cluster at $(busname(rep)) excluded (not every cluster member is measured; a partial sum would bias the fused balance)"
    end
  end

  if !isempty(collapsed_ids)
    shown = join(first(collapsed_ids, 5), ", ")
    @warn "SE: $(length(collapsed_ids)) measurement(s) excluded: their branch collapsed inside a closed-link cluster" rows = shown * (length(collapsed_ids) > 5 ? ", ..." : "")
  end
  return SEWorkingNet(snet, remapped, origIdx, reps, buses, busmap, shuntOrig, true, notes)
end

## Internal helper: assemble measurement value vector z from active measurements.
function _measurement_vector(measurements::Vector{Measurement})
  z = Vector{Float64}(undef, length(measurements))
  @inbounds for (i, m) in enumerate(measurements)
    z[i] = m.value
  end
  return z
end

"""
    _suppression_normalized_residuals(H, r, w) -> Union{Nothing,Vector{Float64}}

Normalized residuals `rn = r_i / sqrt(Omega_ii)` for one suppression round,
through the same `_residual_diagnostics` the bad-data elimination uses, so
both decisions read one scale (task_se_bad_data_v0100).

Returns `nothing` when no verdict is possible: no Jacobian was kept, or the
diagnostics refuse this problem size (they error rather than allocate an
unaffordable dense inverse). The caller then suppresses NOTHING. Falling
back to the raw ratio `|r|/sigma` is deliberately not offered: it is the
measure this task removed, and a suppression decided on the wrong scale is
worse than none, because nothing downstream shows which scale was used.

Returns `rn` together with the residual sensitivity `wii`, because the
caller needs both: `rn` says HOW far the row is off, `wii` says whether
that verdict can be trusted for this row at all.

The `wii` guard of `_residual_diagnostics` carries over unchanged: a row
whose `Omega_ii` sits at the numerical floor is not localizable and gets
`rn = 0`, so it can never be suppressed on the strength of a large raw
residual.
"""
function _suppression_normalized_residuals(H, r::Vector{Float64}, w::Vector{Float64})
  H === nothing && return nothing
  return try
    d = _residual_diagnostics(H, r, w)
    (rn = d.rn, wii = d.wii)
  catch err
    @warn "SE suppression: no residual covariance available, no row is suppressed this round" exception = err
    nothing
  end
end

## Internal helper: compute normalized residuals and row-level diagnostics.
##
## Returns a NamedTuple:
## - rn: normalized residuals r_i / sqrt(Ω_ii)
## - wii: residual-sensitivity diagonal w_ii = Ω_ii * w_i (dimensionless, the
##   diagonal of the residual sensitivity matrix W = Ω W_meas).
##   w_ii near 1 means the residual carries almost the full measurement error
##   (well localizable), w_ii near 0 means a nearly critical measurement whose
##   error hides in the state (not localizable). Literature threshold for a
##   localizable measurement: w_ii > 0.3 (literature threshold).
## - omega: residual covariance Ω = W_meas⁻¹ − H G⁻¹ Hᵀ (needed by the
##   optional residual-correlation K-matrix report).
##
## Two paths (Takahashi task, 0.10.0):
## - :takahashi (default above `state_estimation.takahashi_min_states`
##   states): sparse G = H' W H, UMFPACK lu, and the shared selected inverse
##   on the factor pattern. Ω_ii = 1/w_i - h_i G⁻¹ h_iᵀ needs G⁻¹ only at
##   index pairs (j, k) where row i is nonzero in both: every such pair is
##   a structural nonzero of G itself and therefore inside the factor
##   pattern, so no lookup can miss structurally. A miss anyway (or any
##   Takahashi guard) falls back to the dense path with one warning, never
##   a wrong number. `omega = nothing` on this path (Ω_ii only).
## - :dense (small systems, K-matrix requests via `need_full_omega`, and the
##   fallback): the original dense pinv, full Ω returned.
## Both paths also return `state_variances = diag(G⁻¹)` (confidence
## intervals) and the chosen `omega_path`.
function _residual_diagnostics(H::AbstractMatrix{<:Real}, r::Vector{Float64}, w::Vector{Float64}; need_full_omega::Bool = false, minStates::Int = state_estimation_config().takahashi_min_states)
  m, n = size(H)
  # the full residual covariance is m×m dense by definition; refuse by
  # name instead of allocating (task_se_sparse step 3)
  need_full_omega && m > _SE_FULL_OMEGA_MAX_M && error("SE diagnostics: the residual-correlation K report needs the full m×m residual covariance and m=$(m) exceeds $(_SE_FULL_OMEGA_MAX_M); disable state_estimation.report_residual_correlation for sets this large")
  takahashi_reason = ""
  if !need_full_omega && n >= minStates
    sp, takahashi_reason = _residual_diagnostics_takahashi(H, r, w)
    sp !== nothing && return sp
    @warn "SE diagnostics: Takahashi selected inverse unavailable ($(takahashi_reason)); falling back to the dense path"
  end
  # the dense pinv fallback is cubic in n: below the threshold it is the
  # small-system default and the K-report path; above it no affordable
  # route remains and the estimator says so instead of allocating
  n > _SE_DENSE_LINALG_MAX_N && error("SE diagnostics: no affordable route for n=$(n) states (dense fallback capped at $(_SE_DENSE_LINALG_MAX_N)$(isempty(takahashi_reason) ? "" : "; Takahashi refused: " * takahashi_reason))")

  # Matrix(H) restores the exact zeros the sparse assembly dropped, so this
  # small dense path computes bit for bit what the historical dense H did
  Hd = Matrix{Float64}(H)
  G = (Hd' .* reshape(w, 1, :)) * Hd
  GI = pinv(G)
  WInv = Diagonal(1.0 ./ w)
  Ω = Symmetric(WInv - Hd * GI * Hd')

  rn = zeros(Float64, m)
  wii = zeros(Float64, m)
  for i in eachindex(r)
    σr2 = max(real(Ω[i, i]), 0.0)
    wii[i] = σr2 * w[i]
    if σr2 <= eps(Float64)
      rn[i] = 0.0
    else
      rn[i] = r[i] / sqrt(σr2)
    end
  end
  return (rn = rn, wii = wii, omega = Ω, omega_path = :dense, state_variances = Float64[real(GI[j, j]) for j = 1:n])
end

## LDLt-shaped input for the shared Takahashi pattern pass: CHOLMOD
## factorizes the symmetric G with one symmetric permutation, so the
## p == q guard that rejects UMFPACK's automatic unsymmetric pivoting
## never triggers. Returns nothing when the factorization fails (G
## indefinite or singular; the caller then reports the original reason).
function _takahashi_ldlt_input(G::SparseMatrixCSC{Float64})
  F = try
    ldlt(Symmetric(G); check = false)
  catch
    return nothing
  end
  issuccess(F) || return nothing
  LD = sparse(F.LD)
  n = size(LD, 1)
  d = Vector{Float64}(diag(LD))
  any(iszero, d) && return nothing
  L = tril(LD, -1) + sparse(1.0I, n, n)
  U = sparse(Diagonal(d) * SparseMatrixCSC(transpose(L)))
  return (p = F.p, q = F.p, L = L, U = U, Rs = ones(Float64, n))
end

## Sparse Ω_ii via the shared Takahashi selected inverse. Returns
## (result, reason): result === nothing means the caller must use the dense
## path; reason names why (guard text or the counted structural misses).
function _residual_diagnostics_takahashi(H::AbstractMatrix{<:Real}, r::Vector{Float64}, w::Vector{Float64})
  m, n = size(H)
  # forward differences produce EXACT zeros for measurements not depending
  # on a state, so sparse(H) recovers the true structure without a drop
  # tolerance (verified by test, not assumed); since task_se_sparse H
  # already arrives as a SparseMatrixCSC and passes through untouched
  Hs = H isa SparseMatrixCSC{Float64} ? H : sparse(Matrix{Float64}(H))
  G = Hs' * (Diagonal(w) * Hs)
  F = try
    lu(G)
  catch err
    return nothing, "factorization failed: $(sprint(showerror, err))"
  end
  S, ok, info = takahashi_selected_inverse(F)
  if !ok && occursin("unsymmetric", info)
    # UMFPACK's automatic strategy pivots larger G unsymmetrically
    # (p != q) even though G = H'WH is symmetric, which threw away the
    # whole Takahashi path (step-0 baseline: case1354pegase fell back to
    # the dense pinv for exactly this reason). CHOLMOD's LDLt always uses
    # ONE symmetric permutation and its factors are exactly the
    # (p, q, L, U, Rs) shape the shared pattern pass takes: G[p, p] =
    # L * (D * L') with unit-lower L and no row scaling. Tried only as
    # the second attempt, so every case the UMFPACK route already served
    # keeps its bit-identical factors.
    alt = _takahashi_ldlt_input(G)
    alt === nothing || ((S, ok, info) = takahashi_selected_inverse(alt))
  end
  if !ok
    return (nothing, info)
  end

  Ht = SparseMatrixCSC(transpose(Hs))   # column i = row i of H
  rn = zeros(Float64, m)
  wii = zeros(Float64, m)
  misses = 0
  for i = 1:m
    rng = Ht.colptr[i]:(Ht.colptr[i+1]-1)
    acc = 0.0
    for a in rng
      ja = Ht.rowval[a]
      va = Ht.nzval[a]
      for b in rng
        z = S[ja, Ht.rowval[b]]
        if z === nothing
          misses += 1
        else
          acc += va * z * Ht.nzval[b]
        end
      end
    end
    if misses != 0
      return (nothing, "structural-impossibility: $(misses) in-pattern lookups missed")
    end
    σr2 = max(1.0 / w[i] - acc, 0.0)
    wii[i] = σr2 * w[i]
    rn[i] = σr2 <= eps(Float64) ? 0.0 : r[i] / sqrt(σr2)
  end
  sv = Vector{Float64}(undef, n)
  for j = 1:n
    z = S[j, j]
    if z === nothing
      return (nothing, "structural-impossibility: diagonal lookup missed")
    end
    sv[j] = z
  end
  return (rn = rn, wii = wii, omega = nothing, omega_path = :takahashi, state_variances = sv), ""
end

## Internal helper: per-row maximum absolute residual-correlation coefficient
## max_{j != i} |k_ij| with K = D^(-1/2) Ω D^(-1/2), D = diag(Ω). |k_ij| above
## 1/sqrt(2) marks a simple-redundant group: a bad measurement there is not
## reliably distinguishable from its partner.
function _residual_correlation_max(Ω::AbstractMatrix{Float64})
  m = size(Ω, 1)
  d = [max(real(Ω[i, i]), 0.0) for i = 1:m]
  kmax = zeros(Float64, m)
  for i = 1:m
    d[i] <= eps(Float64) && continue
    for j = 1:m
      (j == i || d[j] <= eps(Float64)) && continue
      k = abs(real(Ω[i, j])) / sqrt(d[i] * d[j])
      k > kmax[i] && (kmax[i] = k)
    end
  end
  return kmax
end

@inline function _chi_square_zscore(j::Float64, ν::Int)
  ν <= 0 && return Inf
  return (j - ν) / sqrt(2.0 * ν)
end

## Ranking rows for the bad-data report. `wii` is the residual-sensitivity
## diagonal from `_residual_diagnostics`; `localizable = wii > wiiThreshold`
## (literature threshold 0.3). `kmax` is the optional per-row maximum
## |correlation coefficient| of the normalized residuals; rows get
## `max_abs_correlation = NaN` and `correlation_warning = false` when the
## K-matrix report is disabled. The warning bound is 1/sqrt(2) (the
## result 7: above it, bad data in a simple-redundant group cannot be told
## apart from its partner).
##
## `measurement_index` is the position in the CALLER's measurement vector.
## SENTINEL: 0 marks a LINKAGG cluster aggregate (SE phase 3), which has no
## single source row. Never index `measurements[row.measurement_index]`
## without checking `>= 1` first; the sequential elimination skips these
## rows for exactly that reason.
function _build_measurement_suspicion_report(
  activeMeas::Vector{Measurement},
  activeIdx::Vector{Int},
  r::Vector{Float64},
  rn::Vector{Float64};
  normalizedThreshold::Float64 = 3.0,
  wii::Union{Nothing,Vector{Float64}} = nothing,
  wiiThreshold::Float64 = 0.3,
  kmax::Union{Nothing,Vector{Float64}} = nothing,
)
  rows = NamedTuple[]
  for i in eachindex(activeMeas)
    m = activeMeas[i]
    wv = wii === nothing ? NaN : wii[i]
    kv = kmax === nothing ? NaN : kmax[i]
    push!(
      rows,
      (
        local_index = i,
        measurement_index = activeIdx[i],
        id = m.id,
        typ = Symbol(m.typ),
        residual = r[i],
        normalized_residual = rn[i],
        abs_normalized_residual = abs(rn[i]),
        suspicious = abs(rn[i]) >= normalizedThreshold,
        wii = wv,
        localizable = !isnan(wv) && wv > wiiThreshold,
        max_abs_correlation = kv,
        correlation_warning = !isnan(kv) && kv > 1.0 / sqrt(2.0),
      ),
    )
  end

  sort!(rows; by = x -> x.abs_normalized_residual, rev = true)
  return rows
end

## Internal helper: assemble diagonal weight entries w = 1/σ² from measurements.
function _weight_vector(measurements::Vector{Measurement})
  w = Vector{Float64}(undef, length(measurements))
  @inbounds for (i, m) in enumerate(measurements)
    w[i] = m.weight
  end
  return w
end

## Internal helper: evaluate nonlinear measurement model h(x).
## Uses nodal complex injections S = V .* conj(Ybus*V) and branch formulas in
## _measurement_prediction. With released shunt states (SE phase 2) the passed
## Ybus has the released shunts' model admittances removed from its diagonal;
## `shuntB` (busIdx -> B state in pu) re-adds the injection analytically:
## exactly the removed term, with the state susceptance in place of the model
## value. Ybus itself stays constant across FD perturbations.
## Released transformer taps work the same way: the four stamped
## admittance terms of a released trafo are removed from the Ybus once
## (`_tap_unstamp!`), and `tapOverlay` re-adds its terminal injections
## analytically at the CASCADE tap position of the current state; branch
## rows on such a trafo evaluate against its scratch branch.
function _predict_measurements(measurements::Vector{Measurement}, net::Net, V::Vector{ComplexF64}, Ybus::AbstractMatrix{ComplexF64}; vaOffsetRad::Float64 = 0.0, shuntB::Union{Nothing,Dict{Int,Float64}} = nothing, tapOverlay = nothing)
  Sbus_MVA = calc_injections(Ybus, V) .* net.baseMVA
  if shuntB !== nothing
    for (busIdx, b) in shuntB
      sh = net.shuntVec[net.shuntDict[busIdx]]
      Sbus_MVA[busIdx] += abs2(V[busIdx]) * conj(complex(real(sh.y_pu_shunt), b)) * net.baseMVA
    end
  end
  if tapOverlay !== nothing
    for (k, ov) in tapOverlay
      br = net.branchVec[k]
      f = Int(br.fromBus)
      t = Int(br.toBus)
      Sbus_MVA[f] += V[f] * conj(ov.Y11 * V[f] + ov.Y12 * V[t]) * net.baseMVA
      Sbus_MVA[t] += V[t] * conj(ov.Y21 * V[f] + ov.Y22 * V[t]) * net.baseMVA
    end
  end
  h = Vector{Float64}(undef, length(measurements))
  @inbounds for (i, m) in enumerate(measurements)
    h[i] = _measurement_prediction(m, net, V, Sbus_MVA, vaOffsetRad, shuntB, tapOverlay)
  end
  return h
end

## Internal helper: finite-difference Jacobian H = ∂h/∂x.
##
## Inputs and meaning:
## - measurements: active measurement set defining h(x)
## - x: current state vector x = [θ(non-slack); Vm(all)] plus, when
##   withVaOffset=true, a trailing PMU reference-offset state α (radians)
## - slackIdx, nbus: mapping context for x -> complex voltages V
## - Ybus: network admittance used inside measurement predictions
## - eps: perturbation size ε for forward-difference sensitivities
##
## Numerical formula (column-wise perturbation):
##   H[i,k] ≈ (h_i(x + εe_k) - h_i(x)) / ε
## where e_k is the k-th unit vector.
##
## Column coloring for the finite-difference assembly (task_se_fd_coloring):
## states that share no measurement row can be perturbed TOGETHER, so the
## number of full prediction evaluations drops from the state count to the
## color count (a few dozen on network structures). Derived from the
## sparsity of an assembled H, deterministic (greedy in column order).
struct _FdColoring
  colors::Vector{Int}              # color of each state column
  groups::Vector{Vector{Int}}      # columns per color, ascending
  rows_per_col::Vector{Vector{Int}}  # the column's measurement-row pattern
  row_covered::Vector{BitVector}   # per color: union of the group's rows
end

function _fd_coloring_from_pattern(H::SparseMatrixCSC{Float64})::_FdColoring
  m, n = size(H)
  rv = rowvals(H)
  rows_per_col = [Int[rv[p] for p in nzrange(H, k)] for k = 1:n]
  return _fd_coloring_from_rows(rows_per_col, m, n)
end

## STRUCTURAL coloring input (task_se_fd_coloring, second iteration of the
## design): the pattern of an assembled H is a POINT property, and on
## case13659pegase micro-derivatives at machine precision appear and
## vanish bitwise between iterations, so no point union ever settled and
## the coloring never engaged. The structural pattern is the superset of
## every position that CAN carry a derivative, taken from the topology:
## a bus measurement couples to its bus and (for injections) its Ybus
## neighbors, a branch measurement to both terminals plus the branch's
## released tap columns, PMU-angle rows additionally to the trailing
## reference-offset column, shunt rows to their released B column. The
## in-assembly detector stays armed: a forgotten coupling shows up as a
## difference outside the pattern and falls back loudly, never silently.
function _fd_structural_coloring(measurements::Vector{Measurement}, x_len::Int, slackIdx::Int, nbus::Int, Ybus::AbstractMatrix{ComplexF64}, withVaOffset::Bool, shuntMap::Union{Nothing,ShuntStateMap}, tapMap::Union{Nothing,TapStateMap}, net::Net)
  # link-flow allocation couples a row across the whole cluster; rare, and
  # not worth modeling here
  any(mm -> mm.linkIdx !== nothing, measurements) && return nothing
  m = length(measurements)
  n = x_len
  vac(b) = b == slackIdx ? 0 : (b < slackIdx ? b : b - 1)
  vmc(b) = (nbus - 1) + b
  shunt_col = Dict{Int,Int}()
  if shuntMap !== nothing
    for j in eachindex(shuntMap.shunt_idxs)
      shunt_col[shuntMap.bus_idxs[j]] = shuntMap.col0 + j - 1
    end
  end
  tap_cols_of_branch = Dict{Int,Vector{Int}}()
  tap_branches_at_bus = Dict{Int,Vector{Int}}()
  if tapMap !== nothing
    for j in eachindex(tapMap.branch_idxs)
      k = tapMap.branch_idxs[j]
      cols = Int[]
      tapMap.r1cols[j] == 0 || push!(cols, tapMap.r1cols[j])
      tapMap.r2cols[j] == 0 || push!(cols, tapMap.r2cols[j])
      tap_cols_of_branch[k] = cols
      br = net.branchVec[k]
      for b in (Int(br.fromBus), Int(br.toBus))
        push!(get!(tap_branches_at_bus, b, Int[]), k)
      end
    end
  end
  offset_col = withVaOffset ? n : 0
  rows_per_col = [Int[] for _ = 1:n]
  cols = Int[]
  for (i, mm) in enumerate(measurements)
    empty!(cols)
    t = mm.typ
    if t == VmMeas
      # BOTH own columns: |V| = Vm * |exp(j theta)| is theta-free only in
      # exact arithmetic; in floats |exp(j theta)| is 1 up to rounding, so
      # the FD picks up a last-bit theta coupling (found on
      # case13659pegase, the detector fired). Same reasoning for Va (the
      # atan2 does not cancel the magnitude bit-exactly) and ShuntQ below.
      push!(cols, vmc(mm.busIdx))
      c = vac(mm.busIdx)
      c == 0 || push!(cols, c)
    elseif t == VaMeas
      c = vac(mm.busIdx)
      c == 0 || push!(cols, c)
      push!(cols, vmc(mm.busIdx))
      offset_col == 0 || push!(cols, offset_col)
    elseif t == PinjMeas || t == QinjMeas
      b = mm.busIdx
      for ptr in nzrange(Ybus, b)
        nb = rowvals(Ybus)[ptr]
        c = vac(nb)
        c == 0 || push!(cols, c)
        push!(cols, vmc(nb))
        for k in get(tap_branches_at_bus, nb, Int[])
          append!(cols, tap_cols_of_branch[k])
        end
      end
      haskey(shunt_col, b) && push!(cols, shunt_col[b])
    elseif t == PflowMeas || t == QflowMeas || t == ImagMeas || t == IaMeas
      br = net.branchVec[mm.branchIdx]
      for b in (Int(br.fromBus), Int(br.toBus))
        c = vac(b)
        c == 0 || push!(cols, c)
        push!(cols, vmc(b))
      end
      append!(cols, get(tap_cols_of_branch, mm.branchIdx, Int[]))
      t == IaMeas && offset_col != 0 && push!(cols, offset_col)
    elseif t == ShuntQMeas
      b = mm.busIdx
      push!(cols, vmc(b))
      c = vac(b)
      c == 0 || push!(cols, c)
      haskey(shunt_col, b) && push!(cols, shunt_col[b])
    else
      # unknown coupling: refuse to color rather than guess
      return nothing
    end
    unique!(sort!(cols))
    for c in cols
      (1 <= c <= n) || return nothing
      push!(rows_per_col[c], i)
    end
  end
  # rows were appended in ascending measurement order per column
  return _fd_coloring_from_rows(rows_per_col, m, n)
end

function _fd_coloring_from_rows(rows_per_col::Vector{Vector{Int}}, m::Int, n::Int)::_FdColoring
  # rows -> columns adjacency for the conflict test
  cols_per_row = [Int[] for _ = 1:m]
  for k = 1:n
    for i in rows_per_col[k]
      push!(cols_per_row[i], k)
    end
  end
  cols_of_row = i -> cols_per_row[i]
  colors = zeros(Int, n)
  ncolors = 0
  forbidden = Int[]   # forbidden[c] == k marks color c as taken for column k
  for k = 1:n
    length(forbidden) < ncolors && append!(forbidden, zeros(Int, ncolors - length(forbidden)))
    for i in rows_per_col[k]
      for k2 in cols_of_row(i)
        c = colors[k2]
        c == 0 && continue
        forbidden[c] = k
      end
    end
    chosen = 0
    for c = 1:ncolors
      if forbidden[c] != k
        chosen = c
        break
      end
    end
    if chosen == 0
      ncolors += 1
      push!(forbidden, 0)
      chosen = ncolors
    end
    colors[k] = chosen
  end
  groups = [Int[] for _ = 1:ncolors]
  for k = 1:n
    push!(groups[colors[k]], k)
  end
  row_covered = [falses(m) for _ = 1:ncolors]
  for c = 1:ncolors
    for k in groups[c]
      for i in rows_per_col[k]
        row_covered[c][i] = true
      end
    end
  end
  return _FdColoring(colors, groups, rows_per_col, row_covered)
end

## Returns:
## - H: Jacobian as a SparseMatrixCSC of size (m_measurements × n_states)
##   (task_se_sparse: the dense m×n allocation was the large-network
##   blocker, tens of gigabytes at 25k buses for a matrix with a handful
##   of entries per row)
## - h0: base prediction vector h(x), reused by caller to avoid recomputation.
##
## With a `coloring` (task_se_fd_coloring) the assembly evaluates ONE
## perturbed prediction per color instead of one per state; the values are
## then written in the SAME column-ascending order with the SAME per-entry
## formula, so the resulting CSC is bit-identical to the per-column
## assembly as long as the coloring's pattern holds. A difference showing
## up OUTSIDE the pattern union of a color's group means an exact zero of
## the base point came alive; the assembly then falls back to the
## per-column path for this call with one warning, correct over fast.
##
## Sparsity contract: the forward difference is computed column-wise
## exactly as the dense version did (same perturbation, same per-entry
## formula, bit-identical values). An entry whose difference is an EXACT
## zero means the prediction does not depend on that state; those exact
## zeros are what the Takahashi diagnostics note below already relies on
## (sparse(H) there recovered the true structure without a drop
## tolerance), so dropping them structurally changes no number downstream.
function _measurement_jacobian_fd(measurements::Vector{Measurement}, net::Net, x::Vector{Float64}, slackIdx::Int, nbus::Int, Ybus::AbstractMatrix{ComplexF64}; eps::Float64 = 1e-6, withVaOffset::Bool = false, shuntMap::Union{Nothing,ShuntStateMap} = nothing, tapMap::Union{Nothing,TapStateMap} = nothing, coloring::Union{Nothing,_FdColoring} = nothing)
  _sb(xv) = shuntMap === nothing ? nothing : _shunt_b_override(xv, shuntMap)
  _tov(xv) = tapMap === nothing ? nothing : _tap_overlay(xv, tapMap, net)
  h0 = _predict_measurements(measurements, net, _state_to_voltage(x, slackIdx, nbus), Ybus; vaOffsetRad = _va_offset_from_state(x, withVaOffset), shuntB = _sb(x), tapOverlay = _tov(x))
  m = length(measurements)
  n = length(x)
  predict_at = xk -> _predict_measurements(measurements, net, _state_to_voltage(xk, slackIdx, nbus), Ybus; vaOffsetRad = _va_offset_from_state(xk, withVaOffset), shuntB = _sb(xk), tapOverlay = _tov(xk))

  # Every caller gets the coloring, not just the WLS loop: the diagnostics
  # and the observability entry points build their own Jacobian, and on a
  # 25k-bus network an uncolored assembly is 50000 full prediction sweeps
  # (found in the maintainer's browser pass, where the estimator itself
  # was already fast). Derived here when the caller passed none.
  if coloring === nothing && n >= _SE_FD_COLORING_MIN_STATES
    coloring = _fd_structural_coloring(measurements, n, slackIdx, nbus, Ybus, withVaOffset, shuntMap, tapMap, net)
  end

  if coloring !== nothing && length(coloring.colors) == n
    # one perturbed prediction PER COLOR; every touched difference must
    # stay inside the color group's row pattern, otherwise an exact zero
    # of the base point came alive and the pattern no longer describes
    # the function: fall back to the per-column path for this call
    ncolors = length(coloring.groups)
    hk_by_color = Vector{Vector{Float64}}(undef, ncolors)
    pattern_ok = true
    for c = 1:ncolors
      # The Jacobian assembly is where a large estimation spends its time,
      # diagnostics included: without a check here an abort waits out the
      # whole 38 s diagnostics phase and looks broken.
      sparlectra_check_abort()
      xk = copy(x)
      for k in coloring.groups[c]
        xk[k] += eps
      end
      hk = predict_at(xk)
      covered = coloring.row_covered[c]
      @inbounds for i = 1:m
        if hk[i] != h0[i] && !covered[i]
          pattern_ok = false
          break
        end
      end
      pattern_ok || break
      hk_by_color[c] = hk
    end
    if pattern_ok
      # SAME column-ascending write order and SAME per-entry formula as
      # the per-column path below: the CSC comes out bit-identical
      rowIdx = Int[]
      colIdx = Int[]
      vals = Float64[]
      sizehint!(rowIdx, 12 * m)
      sizehint!(colIdx, 12 * m)
      sizehint!(vals, 12 * m)
      @inbounds for k = 1:n
        hk = hk_by_color[coloring.colors[k]]
        for i in coloring.rows_per_col[k]
          d = hk[i] - h0[i]
          d == 0.0 && continue
          push!(rowIdx, i)
          push!(colIdx, k)
          push!(vals, d / eps)
        end
      end
      return sparse(rowIdx, colIdx, vals, m, n), h0, true
    end
    @warn "SE Jacobian coloring: a difference outside the color pattern (an exact zero of the base point came alive); falling back to the per-column assembly for this call" maxlog = 1
  end

  # triplet assembly; columns are visited in order and rows ascend inside
  # a column, so `sparse` neither reorders nor combines anything
  rowIdx = Int[]
  colIdx = Int[]
  vals = Float64[]
  sizehint!(rowIdx, 12 * m)
  sizehint!(colIdx, 12 * m)
  sizehint!(vals, 12 * m)

  for k = 1:n
    # same reason as in the colored path: this loop is n full prediction
    # sweeps, so an abort must be able to land inside it
    k % 64 == 0 && sparlectra_check_abort()
    xk = copy(x)
    xk[k] += eps
    hk = predict_at(xk)
    @inbounds for i = 1:m
      d = hk[i] - h0[i]
      d == 0.0 && continue
      push!(rowIdx, i)
      push!(colIdx, k)
      push!(vals, d / eps)
    end
  end

  return sparse(rowIdx, colIdx, vals, m, n), h0, false
end

## Internal helper: create row-wise bipartite adjacency from Jacobian sparsity.
function _adjacency_from_sparsity(H::AbstractMatrix{<:Real})
  m, n = size(H)
  adj = [Int[] for _ = 1:m]
  for i = 1:m
    for j = 1:n
      H[i, j] != 0.0 && push!(adj[i], j)
    end
  end
  return adj, n
end

## Sparse specialization (task_se_sparse): O(nnz) over the stored entries
## instead of the m*n scan; same result, an explicitly stored zero (which
## the FD assembly never produces) is skipped like the dense scan does.
function _adjacency_from_sparsity(H::SparseMatrixCSC{Float64})
  m, n = size(H)
  adj = [Int[] for _ = 1:m]
  rv = rowvals(H)
  nz = nonzeros(H)
  for j = 1:n
    for ptr in nzrange(H, j)
      nz[ptr] != 0.0 && push!(adj[rv[ptr]], j)
    end
  end
  # columns are visited in order, so each row's list is already sorted
  return adj, n
end

## Internal helper: maximum bipartite matching via Hopcroft-Karp.
##
## What this does here:
## - Left side nodes = measurement rows i of Jacobian H.
## - Right side nodes = state columns j of Jacobian H.
## - Edge (i,j) exists if H[i,j] != 0 (measurement i depends on state j).
##
## Why we need it:
## - Structural observability is checked without numeric values, only sparsity pattern.
## - If maximum matching size equals number of state columns, every state can be
##   structurally "covered" by at least one independent measurement relation.
##
## Algorithm note:
## - Hopcroft-Karp alternates BFS (build layers) and DFS (augmenting paths),
##   giving O(E*sqrt(V)) complexity and robust performance for these checks.
function _hopcroft_karp(adj::Vector{Vector{Int}}, nRight::Int)
  nLeft = length(adj)
  pairU = zeros(Int, nLeft)
  pairV = zeros(Int, nRight)
  dist = fill(-1, nLeft)

  function bfs!()
    q = Int[]
    for u = 1:nLeft
      if pairU[u] == 0
        dist[u] = 0
        push!(q, u)
      else
        dist[u] = -1
      end
    end

    found = false
    head = 1
    while head <= length(q)
      u = q[head]
      head += 1
      for v in adj[u]
        u2 = pairV[v]
        if u2 == 0
          found = true
        elseif dist[u2] == -1
          dist[u2] = dist[u] + 1
          push!(q, u2)
        end
      end
    end
    return found
  end

  function dfs!(u::Int)
    for v in adj[u]
      u2 = pairV[v]
      if u2 == 0 || (dist[u2] == dist[u] + 1 && dfs!(u2))
        pairU[u] = v
        pairV[v] = u
        return true
      end
    end
    dist[u] = -1
    return false
  end

  matching = 0
  while bfs!()
    for u = 1:nLeft
      if pairU[u] == 0 && dfs!(u)
        matching += 1
      end
    end
  end

  return matching
end

## Internal helper: test if row i is numerically redundant.
##
## Idea:
## - Remove measurement row i from Jacobian H.
## - If rank(H_without_i) is still full (= number of state columns), the row is
##   not required for numerical observability and is therefore redundant.
## - If not, row i is a numerical critical measurement.
function _numerical_row_redundant(H::AbstractMatrix{<:Real}, i::Int; tol = nothing)
  m, n = size(H)
  @assert 1 <= i <= m
  m <= 1 && return false
  keep = [k for k = 1:m if k != i]
  return numeric_rank(H[keep, :]; tol = tol) == n
end

## Internal helper: test if row i is structurally redundant.
##
## Idea:
## - Build bipartite graph from sparsity pattern of H.
## - Remove all edges of row i.
## - If maximum matching still covers all state columns, the row is structurally
##   redundant; otherwise it is structurally critical.
function _structural_row_redundant(H::AbstractMatrix{<:Real}, i::Int)
  m, n = size(H)
  @assert 1 <= i <= m
  adj, _ = _adjacency_from_sparsity(H)
  adj[i] = Int[]
  return _hopcroft_karp(adj, n) == n
end

## Internal helper: weighted least-squares objective J = r'Wr = Σ_i w_i r_i².
@inline _wls_objective(r::Vector{Float64}, w::Vector{Float64}) = sum(w .* (r .^ 2))

## Internal helper: Wilson-Hilferty z-score of the χ² band test (SE phase 4).
## (J/ν)^(1/3) is approximately normal with mean 1 - 2/(9ν) and variance
## 2/(9ν), usable from ν ≈ 3 and consistent for all ν >= 1; unlike the plain
## normal approximation (valid only for ν > 30) it respects the right skew
## of the χ² distribution at small redundancy.
@inline function _wilson_hilferty_zscore(j::Float64, ν::Int)
  ν <= 0 && return Inf
  μ = 1.0 - 2.0 / (9.0 * ν)
  σ = sqrt(2.0 / (9.0 * ν))
  return (cbrt(max(j, 0.0) / ν) - μ) / σ
end

## Full band-test verdict (SE phase 4): pass/fail via Wilson-Hilferty with a
## two-sided reason. `:high` = J too large (bad data, model error); `:low` =
## J implausibly small (sigmas overestimated, the typical signature of
## noise-free synthetic data); `:no_redundancy` = ν = 0, J is structurally 0
## and bad data is invisible (test skipped, reported as not passed);
## `:small_redundancy` flags ν < 30 informationally (the WH test stays
## valid, but the relative spread sqrt(2/ν) reduces its power).
function _band_test_verdict(j::Float64, ν::Int)
  if ν <= 0
    return (passed = false, reason = :no_redundancy, z_wh = Inf, z_legacy = Inf, small_redundancy = true)
  end
  z_wh = _wilson_hilferty_zscore(j, ν)
  z_legacy = (j - ν) / sqrt(2.0 * ν)
  passed = abs(z_wh) <= 3.0
  reason = passed ? :ok : (z_wh > 0.0 ? :high : :low)
  return (passed = passed, reason = reason, z_wh = z_wh, z_legacy = z_legacy, small_redundancy = ν < 30)
end

## Backward-compatible name: the band test now DECIDES via Wilson-Hilferty
## for all ν; the old symmetric z-score is only reported for continuity.
function _j_within_3sigma_band(j::Float64, ν::Int)
  return _band_test_verdict(j, ν).passed
end

## Internal helper: classify redundancy/observability quality.
function _redundancy_quality(observable::Bool, ν::Int, hasCritical::Bool)
  if !observable
    return :not_observable
  elseif ν <= 0 || hasCritical
    return :critical
  else
    return :good
  end
end

## Internal helper: compute numerical + structural observability metrics from Jacobian.
##
## Computed quantities:
## - Numerical observability via matrix rank(H)
## - Structural observability via maximum bipartite matching
## - Redundancy metrics ν = m - n and ρ = m/n
## - Critical measurements by single-row removal tests (numerical + structural)
##
## activeOriginalIdx maps local Jacobian rows back to original measurement indices.
"""
    _column_normalized(H) -> SparseMatrixCSC

Normalize every column of `H` to unit norm for the rank test.

Rank does not change under column scaling, but its NUMERICAL determination
does. Unscaled, a voltage column (entries around 1) sits next to an angle
column whose entries carry the branch admittances of a 765 kV line, so the
columns span orders of magnitude; `sigma_max` then belongs to the largest
scale and genuine small singular values fall below a relative cut. That
read as a rank deficit on a 25000-bus network whose set was four times
overdetermined.

Rows are deliberately NOT scaled by sigma. The rank is invariant under
positive row scaling, so sigma carries no information for this question,
and dividing by it is actively harmful: zero-injection pseudo-measurements
carry sigma 1e-6 against 0.01 to 24 for real ones, so the division lifts
those rows by six orders of magnitude, they take over `sigma_max`, and
everything else drops below the cut. Measured on that network: raw deficit
191, with row scaling 12308, with column normalization alone 0.

A column of norm zero carries no measurement at all and stays zero, so it
keeps producing the rank deficit it should.
"""
function _column_normalized(H::AbstractMatrix{<:Real})
  Hs = sparse(H)
  vals = nonzeros(Hs)
  for j = 1:size(Hs, 2)
    r = nzrange(Hs, j)
    isempty(r) && continue
    nrm = sqrt(sum(abs2, view(vals, r)))
    (isfinite(nrm) && nrm > 0.0) || continue
    for k in r
      vals[k] /= nrm
    end
  end
  return Hs
end

function _evaluate_observability_from_jacobian(H::AbstractMatrix{<:Real}, activeOriginalIdx::Vector{Int}; tol = nothing)
  m, n = size(H)
  ν = m - n
  ρ = n > 0 ? m / n : Inf

  nrank = numeric_rank(H; tol = tol)
  adj, ncols = _adjacency_from_sparsity(H)
  mm = _hopcroft_karp(adj, ncols)

  numObs = nrank == ncols
  strObs = mm == ncols

  # single-row criticality is one rank test per row, so its cost is
  # m times a full decomposition. Step-0 baseline: m*n = 210k (sp_case188,
  # 562 rows) took 7.6 s and 2.9 GB; m*n = 11M (case1354) did not
  # terminate in reasonable time. The m*n budget below keeps the check
  # within roughly a minute; above it the classification is SKIPPED and
  # says so (criticality_skipped), instead of allocating for hours, and
  # quality then reflects observability and redundancy only
  # (task_se_sparse step 3).
  criticalNum = Int[]
  criticalStr = Int[]
  criticalitySkipped = m * n > _SE_ROW_CRITICALITY_MAX_MN
  if m > 0 && !criticalitySkipped
    for i = 1:m
      _numerical_row_redundant(H, i; tol = tol) || push!(criticalNum, activeOriginalIdx[i])
      _structural_row_redundant(H, i) || push!(criticalStr, activeOriginalIdx[i])
    end
  elseif criticalitySkipped
    @warn "observability: single-row criticality skipped (m=$(m) times n=$(n) exceeds the $(_SE_ROW_CRITICALITY_MAX_MN) budget); the per-row rank tests do not pay for themselves at this size"
  end

  hasCritical = !isempty(criticalNum) || !isempty(criticalStr)

  # Unobservable state columns (the "dark" states): computed only on the
  # not-observable path so observable sets pay nothing extra. The (numerical)
  # null space of H spans exactly the state directions no measurement sees;
  # a state column with a nonzero component in ANY basis vector cannot be
  # estimated (its value can drift inside the null space without changing
  # a single measurement). SE phase 4: the basis is taken from the SVD right
  # singular vectors of the singular values BELOW the effective rank
  # tolerance, so an FD-aware `tol` (see the observability entry points)
  # finds the same near-null directions the rank decision saw; a plain
  # `nullspace` with LAPACK's eps-scale default would miss them on
  # forward-difference Jacobians. Dense SVD cost: acceptable at workshop and
  # distribution-network sizes, and it only runs on failures.
  dark = Int[]
  if !numObs && n > 0 && m > 0
    if max(m, n) > _SE_DENSE_LINALG_MAX_N
      @warn "observability: dark-state identification skipped above $(_SE_DENSE_LINALG_MAX_N) states/measurements (m=$(m), n=$(n)); the dense null-space SVD does not pay for itself at this size"
    else
      F = svd(Matrix(H); full = true)
      for k = (nrank+1):n
        for j = 1:n
          abs(F.V[j, k]) > 1e-8 && push!(dark, j)
        end
      end
      unique!(sort!(dark))
    end
  end

  return (
    numerical_observable = numObs,
    structural_observable = strObs,
    numerical_rank = nrank,
    structural_matching = mm,
    n_states = ncols,
    n_measurements = m,
    redundancy = ν,
    redundancy_ratio = ρ,
    dof = ν,
    numerical_critical_measurement_indices = criticalNum,
    structural_critical_measurement_indices = criticalStr,
    criticality_skipped = criticalitySkipped,
    unobservable_state_columns = dark,
    quality = _redundancy_quality(numObs && strObs, ν, hasCritical),
  )
end

"""
    numerical_observable(H; tol=nothing) -> Bool

Numerical observability test on a Jacobian-like matrix `H`.
Returns `true` when `rank(H) == n` (full column rank).
"""
function numerical_observable(H::AbstractMatrix{<:Real}; tol = nothing)
  _, n = size(H)
  return numeric_rank(H; tol = tol) == n
end

"""
    structural_observable(H) -> Bool

Structural observability test on a Jacobian-like matrix `H`.
Returns `true` when the maximum bipartite matching size equals the number of
state columns `n`.
"""
function structural_observable(H::AbstractMatrix{<:Real})
  _, n = size(H)
  adj, _ = _adjacency_from_sparsity(H)
  return _hopcroft_karp(adj, n) == n
end

"""
    numerical_row_redundant(H, i; tol=nothing) -> Bool

Check if row `i` remains numerically redundant in `H`.
"""
function numerical_row_redundant(H::AbstractMatrix{<:Real}, i::Int; tol = nothing)
  return _numerical_row_redundant(H, i; tol = tol)
end

"""
    structural_row_redundant(H, i) -> Bool

Check if row `i` remains structurally redundant in `H`.
"""
function structural_row_redundant(H::AbstractMatrix{<:Real}, i::Int)
  return _structural_row_redundant(H, i)
end

"""
    evaluate_observability_matrix(H; tol=nothing) -> NamedTuple

Evaluate global observability and single-row criticality directly on a matrix
`H` (without building a network model).

The result includes `unobservable_state_columns::Vector{Int}`: for a
numerically NOT observable `H`, the state columns with a component above
tolerance in any null-space basis vector, i.e. exactly the states no
measurement pins down (their union partitions the network into observable
islands). Empty for observable systems; computed only on the
not-observable path (a dense null-space probe, fine at workshop and
distribution-network sizes).
"""
function evaluate_observability_matrix(H::AbstractMatrix{<:Real}; tol = nothing)
  m, _ = size(H)
  idx = collect(1:m)
  return _evaluate_observability_from_jacobian(H, idx; tol = tol)
end

"""
    evaluate_local_observability_matrix(H, stateCols; tol=nothing) -> NamedTuple

Evaluate local observability on a matrix `H` restricted to selected
`stateCols`.

Limitation: this column-restricted submatrix test is NECESSARY but not
sufficient; a positive verdict can be wrong when the touching rows couple
the selected states to neighbor states that are themselves undetermined
(the rigorous per-state answer is `unobservable_state_columns` from the
global check).
"""
function evaluate_local_observability_matrix(H::AbstractMatrix{<:Real}, stateCols::Vector{Int}; tol = nothing)
  isempty(stateCols) && error("evaluate_local_observability_matrix: stateCols must not be empty")
  nstates = size(H, 2)
  for c in stateCols
    1 <= c <= nstates || error("evaluate_local_observability_matrix: state column out of bounds: $(c)")
  end

  rows = Int[]
  for i in axes(H, 1)
    for j in stateCols
      if !iszero(H[i, j])
        push!(rows, i)
        break
      end
    end
  end

  isempty(rows) && error("evaluate_local_observability_matrix: no row touches selected stateCols")
  local_idx = [rows[k] for k in eachindex(rows)]
  base = _evaluate_observability_from_jacobian(H[rows, stateCols], local_idx; tol = tol)
  return merge(base, (rows = rows, stateCols = copy(stateCols)))
end

"""
    evaluate_global_observability(net, measurements; kwargs...) -> NamedTuple

Evaluate global observability on active measurements using the finite-difference
measurement Jacobian.

Includes global redundancy metrics
- `redundancy = r = m - n`
- `redundancy_ratio = ρ = m / n`
- `dof = ν = m - n` (a COUNT difference: for an observable set it equals
  `m - rank(H)`; for an unobservable set it can be negative and then reads
  as a shortfall, not a redundancy)

Quality classes:
- `:good`: observable and no critical single measurement
- `:critical`: observable, but at least one single critical measurement (or ν <= 0)
- `:not_observable`: not observable

For a not-observable set the result additionally names the dark states in
`unobservable_state_columns` (state columns touched by the null space of
`H`); empty when observable. This is the rigorous per-state answer that
the column-restricted local check cannot give (see
[`evaluate_local_observability`](@ref)).

Current-magnitude measurements (`ImagMeas`) are excluded before the Jacobian
is built: currents are auxiliary measurements that must never carry
observability, so the verdict states whether the system is
observable from the power and voltage measurements alone.

The verdict is a two-stage check (SE phase 4):

1. **Structural**: `detect_ac_islands` on the contracted SE net. More than
   one synchronous island containing measured buses yields
   `quality = :not_observable` with `:structural_islands` in `notes`,
   regardless of numeric rank (no measurement ties the islands' angle
   references together).
2. **FD-aware numeric rank**: the Jacobian is built by forward differences,
   so its error floor is O(`jacEps`), not machine epsilon; an eps-scale SVD
   tolerance would count FD noise as rank. With `tol = nothing` the rank
   tolerance therefore defaults to
   `state_estimation.rank_tol_factor * jacEps * sigma_max` (factor default
   10.0). An explicitly passed `tol` wins.
"""
function evaluate_global_observability(net::Net, measurements::Vector{Measurement}; flatstart::Bool = true, jacEps::Float64 = 1e-6, tol = nothing, pmuRefOffset::Symbol = state_estimation_config().pmu_ref_offset, rankTolFactor::Float64 = state_estimation_config().rank_tol_factor)
  # island-wise nets: judge every measured island on its own subnet with its
  # own reference (the estimator solves them the same way) and aggregate;
  # unmeasured islands are excluded from the estimation and reported
  part = _se_island_partition(net, measurements)
  if length(part.rows) > 1
    return _evaluate_global_observability_islands(net, measurements, part; flatstart = flatstart, jacEps = jacEps, tol = tol, pmuRefOffset = pmuRefOffset, rankTolFactor = rankTolFactor)
  end

  # SE phase 3: observability is judged on the contracted net (links fused),
  # with the same measurement remapping the estimator uses
  prep = _se_prepare(net, measurements)
  net = prep.snet
  activeMeas = prep.meas
  activeIdx = prep.origIdx
  # currents must not carry observability: drop ImagMeas rows up front
  keep = [k for k in eachindex(activeMeas) if !(activeMeas[k].typ in (ImagMeas, IaMeas))]
  activeMeas = activeMeas[keep]
  activeIdx = activeIdx[keep]
  isempty(activeMeas) && error("evaluate_global_observability: no active measurements (current magnitudes are excluded; the system must be observable without them)")

  nbus = length(net.nodeVec)
  slackIdx = _find_slack_idx(net)
  Ybus = createYBUS(net = net)
  withVaOffset = _va_offset_active(activeMeas, pmuRefOffset)
  x = _initial_state_vector(net, slackIdx; flatstart = flatstart, withVaOffset = withVaOffset)
  H, _ = _measurement_jacobian_fd(activeMeas, net, x, slackIdx, nbus, Ybus; eps = jacEps, withVaOffset = withVaOffset)

  # Normalize columns before the rank test: one tolerance then means the
  # same thing at 14 buses and at 25000 (see _column_normalized).
  Hs = _column_normalized(H)
  effTol = tol
  if effTol === nothing && !isempty(Hs)
    effTol = rankTolFactor * jacEps * _sigma_max(Hs)
  end
  base = _evaluate_observability_from_jacobian(Hs, activeIdx; tol = effTol)
  cond_notes = Symbol[]

  # stage 1: structural island check on the contracted net. More than one
  # synchronous island with measured buses means no measurement chain ties
  # the angle references together; the set is not observable no matter what
  # the FD rank says.
  notes = copy(cond_notes)
  islands = detect_ac_islands(net)
  if length(islands.rows) > 1
    measuredBuses = Set{Int}()
    for m in activeMeas
      m.busIdx !== nothing && push!(measuredBuses, m.busIdx)
      if m.branchIdx !== nothing && 1 <= m.branchIdx <= length(net.branchVec)
        br = net.branchVec[m.branchIdx]
        push!(measuredBuses, Int(br.fromBus))
        push!(measuredBuses, Int(br.toBus))
      end
    end
    nMeasuredIslands = count(row -> any(b -> b in measuredBuses, row.buses), islands.rows)
    if nMeasuredIslands > 1
      push!(notes, :structural_islands)
      return merge(base, (quality = :not_observable, numerical_observable = false, notes = notes))
    end
  end

  return merge(base, (notes = notes,))
end

## per-island observability merge (multi-island nets): quality is the WORST
## verdict over the measured islands, counts are summed, critical-measurement
## indices are remapped into the caller's measurement indexing. Unmeasured
## islands do not enter the verdict (the estimator skips them); their count
## is reported. `unobservable_state_columns` stays per island (inside the
## `islands` rows), the merged field is left empty.
function _evaluate_global_observability_islands(net::Net, measurements::Vector{Measurement}, part; kwargs...)
  sev(q) = q == :not_observable ? 3 : (q == :critical ? 2 : 1)
  islandsOut = NamedTuple[]
  worst = :good
  nMeas = 0
  nStates = 0
  nrank = 0
  mm = 0
  critNum = Int[]
  critStr = Int[]
  measured = 0
  unmeasured = 0
  for (iid, row) in enumerate(part.rows)
    idxs = part.perIsland[iid]
    if isempty(idxs)
      unmeasured += 1
      push!(islandsOut, (island = iid, n_bus = row.n_bus, measured = false, result = nothing))
      continue
    end
    sub = _se_island_subnet(net, row)
    meas_i, keptIdx = _se_remap_island_measurements(sub, measurements, idxs)
    if isempty(meas_i)
      unmeasured += 1
      push!(islandsOut, (island = iid, n_bus = row.n_bus, measured = false, result = nothing))
      continue
    end
    r = evaluate_global_observability(sub.inet, meas_i; kwargs...)
    measured += 1
    sev(r.quality) > sev(worst) && (worst = r.quality)
    nMeas += r.n_measurements
    nStates += r.n_states
    nrank += r.numerical_rank
    mm += r.structural_matching
    append!(critNum, Int[keptIdx[i] for i in r.numerical_critical_measurement_indices if 1 <= i <= length(keptIdx)])
    append!(critStr, Int[keptIdx[i] for i in r.structural_critical_measurement_indices if 1 <= i <= length(keptIdx)])
    push!(islandsOut, (island = iid, n_bus = row.n_bus, measured = true, result = r))
  end
  measured == 0 && error("evaluate_global_observability: no active measurements (current magnitudes are excluded; the system must be observable without them)")
  ν = nMeas - nStates
  return (
    numerical_observable = worst != :not_observable,
    structural_observable = worst != :not_observable,
    numerical_rank = nrank,
    structural_matching = mm,
    n_states = nStates,
    n_measurements = nMeas,
    redundancy = ν,
    redundancy_ratio = nStates > 0 ? nMeas / nStates : Inf,
    dof = ν,
    numerical_critical_measurement_indices = critNum,
    structural_critical_measurement_indices = critStr,
    unobservable_state_columns = Int[],
    quality = worst,
    notes = Symbol[:island_partition],
    islands = islandsOut,
    n_measured_islands = measured,
    n_unmeasured_islands = unmeasured,
  )
end

function evaluate_global_observability(net::Net; flatstart::Bool = true, jacEps::Float64 = 1e-6, tol = nothing, pmuRefOffset::Symbol = state_estimation_config().pmu_ref_offset, rankTolFactor::Float64 = state_estimation_config().rank_tol_factor)
  return evaluate_global_observability(net, Measurement[m for m in net.measurements]; flatstart = flatstart, jacEps = jacEps, tol = tol, pmuRefOffset = pmuRefOffset, rankTolFactor = rankTolFactor)
end

"""
    measurement_jacobian(net; flatstart=true, jacEps=1e-6, pmuRefOffset=...) -> NamedTuple

Build the measurement Jacobian `H` of the ACTIVE measurements on `net`,
labeled for humans: the matrix the observability checks and the WLS
normal equations run on, with one described row per measurement and one
described column per state.

Returns `(H, rows, cols)`:
- `H::Matrix{Float64}`: m x n finite-difference Jacobian at the flat (or
  stored) start state, the same evaluation `evaluate_global_observability`
  uses.
- `rows`: one NamedTuple per active measurement,
  `(index, type, location, sigma)`; `index` is the position in
  `net.measurements`, `location` names the bus (injections, voltages) or
  the oriented branch (flows).
- `cols`: state-column labels in Jacobian order, `"Va(bus)"` for every
  non-slack bus, then `"Vm(bus)"` for every bus, plus `"alpha"` when PMU
  `Va` measurements activate the reference-offset state.

Errors when no active measurement exists. Intended for measurement-matrix
reports and placement studies; see the state-estimation suite summary and
the workshop's observability deep dive.
"""
function measurement_jacobian(net::Net; flatstart::Bool = true, jacEps::Float64 = 1e-6, pmuRefOffset::Symbol = state_estimation_config().pmu_ref_offset)
  # SE phase 3: the labeled H describes the contracted net (links fused)
  prep = _se_prepare(net, Measurement[m for m in net.measurements])
  net = prep.snet
  activeMeas = prep.meas
  activeIdx = prep.origIdx
  isempty(activeMeas) && error("measurement_jacobian: no active measurements")
  nbus = length(net.nodeVec)
  slackIdx = _find_slack_idx(net)
  Ybus = createYBUS(net = net)
  withVaOffset = _va_offset_active(activeMeas, pmuRefOffset)
  x = _initial_state_vector(net, slackIdx; flatstart = flatstart, withVaOffset = withVaOffset)
  H, _ = _measurement_jacobian_fd(activeMeas, net, x, slackIdx, nbus, Ybus; eps = jacEps, withVaOffset = withVaOffset)

  name_by_idx = _bus_name_by_idx(net)
  busname(i) = get(name_by_idx, i, string(i))
  rows = NamedTuple[]
  for (k, m) in enumerate(activeMeas)
    location = if m.typ in (PflowMeas, QflowMeas)
      br = net.branchVec[m.branchIdx]
      ends = (busname(Int(br.fromBus)), busname(Int(br.toBus)))
      m.direction === :to ? string(ends[2], "->", ends[1]) : string(ends[1], "->", ends[2])
    else
      busname(m.busIdx)
    end
    push!(rows, (index = activeIdx[k], type = m.typ, location = location, sigma = m.sigma))
  end
  cols = String[]
  for i = 1:nbus
    i == slackIdx || push!(cols, string("Va(", busname(i), ")"))
  end
  for i = 1:nbus
    push!(cols, string("Vm(", busname(i), ")"))
  end
  withVaOffset && push!(cols, "alpha")
  return (H = H, rows = rows, cols = cols)
end

"""
    evaluate_local_observability(net, measurements, stateCols; kwargs...) -> NamedTuple

Evaluate local observability on selected Jacobian columns (`stateCols`).

Procedure:
1) Build global Jacobian `H` from currently active measurements.
2) Keep only rows that have at least one nonzero entry in the selected columns.
   These rows correspond to measurements that are locally sensitive to the
   requested states.
3) Evaluate observability/redundancy on the reduced matrix `Hlocal`.

Returned NamedTuple extends global metrics with:
- `rows`: selected row indices (within global active-Jacobian row numbering)
- `stateCols`: copied input state-column selection.

Interpretation:
- `:good` means local states are observable with positive redundancy and no
  single critical measurement.
- `:critical` means still observable but vulnerable to a single outage (or ν <= 0).
- `:not_observable` means local states cannot be uniquely reconstructed.

Limitation: this column-restricted test is NECESSARY but not sufficient; a
positive verdict can be wrong when the touching rows couple the selected
states to neighbor states that are themselves undetermined. The rigorous
per-state answer is `unobservable_state_columns` from
[`evaluate_global_observability`](@ref).

Like the global check, current-magnitude measurements (`ImagMeas`) are
excluded before the Jacobian is built (currents must never carry
observability), and with `tol = nothing` the rank tolerance is FD-aware
(`rank_tol_factor * jacEps * sigma_max` of the local submatrix; see
[`evaluate_global_observability`](@ref)). The structural island stage is
deliberately global-only: a column selection inside ONE island is perfectly
observable even when the net contains further measured islands.
"""
function evaluate_local_observability(net::Net, measurements::Vector{Measurement}, stateCols::Vector{Int}; flatstart::Bool = true, jacEps::Float64 = 1e-6, tol = nothing, pmuRefOffset::Symbol = state_estimation_config().pmu_ref_offset, rankTolFactor::Float64 = state_estimation_config().rank_tol_factor)
  isempty(stateCols) && error("evaluate_local_observability: stateCols must not be empty")

  # SE phase 3: judged on the contracted net, like the global check
  prep = _se_prepare(net, measurements)
  net = prep.snet
  activeMeas = prep.meas
  activeIdx = prep.origIdx
  # currents must not carry observability: drop ImagMeas rows up front
  keep = [k for k in eachindex(activeMeas) if !(activeMeas[k].typ in (ImagMeas, IaMeas))]
  activeMeas = activeMeas[keep]
  activeIdx = activeIdx[keep]
  isempty(activeMeas) && error("evaluate_local_observability: no active measurements (current magnitudes are excluded)")

  nbus = length(net.nodeVec)
  slackIdx = _find_slack_idx(net)
  Ybus = createYBUS(net = net)
  withVaOffset = _va_offset_active(activeMeas, pmuRefOffset)
  x = _initial_state_vector(net, slackIdx; flatstart = flatstart, withVaOffset = withVaOffset)
  H, _ = _measurement_jacobian_fd(activeMeas, net, x, slackIdx, nbus, Ybus; eps = jacEps, withVaOffset = withVaOffset)

  nstates = size(H, 2)
  for c in stateCols
    1 <= c <= nstates || error("evaluate_local_observability: state column out of bounds: $(c)")
  end

  localRows = Int[]

  for i in axes(H, 1)
    for j in stateCols
      if !iszero(H[i, j])
        push!(localRows, i)
        break
      end
    end
  end

  isempty(localRows) && error("evaluate_local_observability: no active measurement touches selected stateCols")

  Hlocal = H[localRows, stateCols]
  localOriginalIdx = [activeIdx[i] for i in localRows]
  # FD-aware rank tolerance on the local submatrix (SE phase 4); an explicit
  # tol wins
  effTol = tol
  if effTol === nothing && !isempty(Hlocal)
    effTol = rankTolFactor * jacEps * _sigma_max(Hlocal)
  end
  base = _evaluate_observability_from_jacobian(Hlocal, localOriginalIdx; tol = effTol)

  return merge(base, (rows = localRows, stateCols = copy(stateCols)))
end

function evaluate_local_observability(net::Net, stateCols::Vector{Int}; flatstart::Bool = true, jacEps::Float64 = 1e-6, tol = nothing, pmuRefOffset::Symbol = state_estimation_config().pmu_ref_offset)
  return evaluate_local_observability(net, Measurement[m for m in net.measurements], stateCols; flatstart = flatstart, jacEps = jacEps, tol = tol, pmuRefOffset = pmuRefOffset)
end

"""
    runse!(net, measurements; kwargs...) -> SEResult

Run a first classical nonlinear weighted least-squares state estimator.

State representation:
- bus voltage angles for all non-slack buses (radians)
- bus voltage magnitudes for all buses (p.u.)
- optional PMU reference-angle offset α (radians), appended automatically
  when active `VaMeas` measurements exist and `pmu_ref_offset = :auto`

PMU angle measurements (`VaMeas`, degrees) are modeled as
`z = θ_i + α + e`: the network angles stay slack-referenced, α maps them
into the common PMU time base. With `pmu_ref_offset = :off` the offset
state is omitted and PMU angles are assumed to be slack-referenced already.
"""
# ---------------------------------------------------------------------------
# Island-wise state estimation (per-island reference)
# ---------------------------------------------------------------------------
# A net can split into several synchronous AC islands (multi-area MATPOWER
# cases, CGMES deliveries with stub islands). Mirroring the power flow's
# island-wise solving, the estimator partitions the measurement set onto the
# islands, estimates every island that carries measurements with its OWN
# angle reference (an island without a slack gets the reference bus
# detect_ac_islands already chose, same MATPOWER-like promotion as the PF),
# and merges the results. Unmeasured islands are skipped and reported.

## kwdef-config copy with overrides (no Setfield dependency)
function _copy_se_config_with(cfg::StateEstimationConfig; kwargs...)
  vals = Dict{Symbol,Any}(f => getfield(cfg, f) for f in fieldnames(StateEstimationConfig))
  for (k, v) in kwargs
    vals[k] = v
  end
  return StateEstimationConfig(; vals...)
end

## partition measurement positions onto the AC islands; `unassigned` collects
## rows whose bus/branch/link belongs to no island (isolated buses).
## detect_ac_islands sees only branch connectivity, but a CLOSED link fuses
## its two buses into one electrical node, so islands connected by closed
## links are merged into one estimation group here (the estimator's own
## contraction then handles the cluster); tearing such a cluster apart
## would break the LINKAGG semantics and drop the link-flow measurements.
function _se_island_partition(net::Net, measurements::Vector{Measurement})
  rep = detect_ac_islands(net)
  raw = rep.rows
  busIsland0 = Dict{Int,Int}()
  for (iid, row) in enumerate(raw)
    for b in row.buses
      busIsland0[b] = iid
    end
  end
  # union-find over island ids via closed links
  parent = collect(1:length(raw))
  function findroot(i)
    while parent[i] != i
      parent[i] = parent[parent[i]]
      i = parent[i]
    end
    return i
  end
  for l in net.linkVec
    l.status == 1 || continue
    a = get(busIsland0, Int(l.fromBus), 0)
    b = get(busIsland0, Int(l.toBus), 0)
    (a == 0 || b == 0) && continue
    ra = findroot(a)
    rb = findroot(b)
    ra != rb && (parent[rb] = ra)
  end
  groups = Dict{Int,Vector{Int}}()
  for i in eachindex(raw)
    push!(get!(groups, findroot(i), Int[]), i)
  end
  rows = NamedTuple[]
  for members in sort!(collect(values(groups)); by = first)
    if length(members) == 1
      push!(rows, raw[members[1]])
    else
      buses = reduce(vcat, (raw[m].buses for m in members))
      branches = reduce(vcat, (raw[m].branches for m in members))
      n_ref = sum(raw[m].n_ref for m in members)
      chosen = 0
      for m in members
        raw[m].chosen_ref_bus > 0 && (chosen = raw[m].chosen_ref_bus; break)
      end
      push!(rows, (buses = buses, branches = branches, n_bus = length(buses), n_ref = n_ref, chosen_ref_bus = chosen))
    end
  end
  busIsland = Dict{Int,Int}()
  for (iid, row) in enumerate(rows)
    for b in row.buses
      busIsland[b] = iid
    end
  end
  perIsland = [Int[] for _ in rows]
  unassigned = Int[]
  for (k, m) in enumerate(measurements)
    m.active || continue
    iid = 0
    if m.linkIdx !== nothing && 1 <= m.linkIdx <= length(net.linkVec)
      iid = get(busIsland, Int(net.linkVec[m.linkIdx].fromBus), 0)
    elseif m.branchIdx !== nothing && 1 <= m.branchIdx <= length(net.branchVec)
      iid = get(busIsland, Int(net.branchVec[m.branchIdx].fromBus), 0)
    elseif m.busIdx !== nothing
      iid = get(busIsland, m.busIdx, 0)
    end
    iid == 0 ? push!(unassigned, k) : push!(perIsland[iid], k)
  end
  return (rows = rows, perIsland = perIsland, unassigned = unassigned)
end

## island subnet plus the index maps needed to remap measurements and to
## mirror shunt write-backs onto the parent net
function _se_island_subnet(net::Net, row)
  inet = _prepare_island_net(net, row)
  # _prepare_island_net promotes a reference only when detect_ac_islands
  # found a regulating generator (the PF's slack must balance power). The
  # estimator only needs an ANGLE PIN, so an island without any candidate
  # (an unpowered stub) still gets one: its first bus. The pin must be
  # prosumer-backed, because refreshBusTypesFromProsumers! (run by every
  # detect_ac_islands call) re-derives node types and would wipe a bare
  # type assignment. The reference prosumer carries no injection and the
  # bus magnitude stays an ordinary state.
  if !any(getNodeType(nd) == Slack for nd in inet.nodeVec)
    refname = getCompName(inet.nodeVec[1].comp)
    addProsumer!(net = inet, busName = refname, type = "EXTERNALNETWORKINJECTION", referencePri = refname, vm_pu = 1.0, va_deg = 0.0)
    refreshBusTypesFromProsumers!(inet)
    push!(inet.slackVec, 1)
  end
  busmap = Dict{Int,Int}(old => new for (new, old) in enumerate(row.buses))
  branchmap = Dict{Int,Int}(old => new for (new, old) in enumerate(row.branches))
  busset = Set(row.buses)
  # _prepare_island_net deepcopies linkVec with stale net-wide indices (the
  # PF contracts links before splitting and never reads them); the SE runs
  # its own contraction on the subnet, so rebuild the vector properly
  linkmap = Dict{Int,Int}()
  links = eltype(net.linkVec)[]
  for (k, l) in enumerate(net.linkVec)
    (Int(l.fromBus) in busset && Int(l.toBus) in busset) || continue
    nl = deepcopy(l)
    nl.fromBus = busmap[Int(l.fromBus)]
    nl.toBus = busmap[Int(l.toBus)]
    nl.linkIdx = length(links) + 1
    push!(links, nl)
    linkmap[k] = nl.linkIdx
  end
  inet.linkVec = links
  shuntpos = [k for (k, sh) in enumerate(net.shuntVec) if Int(sh.busIdx) in busset]
  return (inet = inet, busmap = busmap, branchmap = branchmap, linkmap = linkmap, shuntpos = shuntpos)
end

## remap the island's measurement rows into subnet indexing; returns the
## remapped vector plus, per row, the position in the caller's vector
function _se_remap_island_measurements(sub, measurements::Vector{Measurement}, idxs::Vector{Int})
  meas = Measurement[]
  keptIdx = Int[]
  for k in idxs
    m = measurements[k]
    remapped = if m.linkIdx !== nothing
      nl = get(sub.linkmap, m.linkIdx, 0)
      if nl == 0
        @warn "SE: measurement $(m.id) excluded (its link is not inside the island)"
        nothing
      else
        Measurement(typ = m.typ, value = m.value, sigma = m.sigma, active = m.active, direction = m.direction, id = m.id, linkIdx = nl)
      end
    elseif m.branchIdx !== nothing
      nb = get(sub.branchmap, m.branchIdx, 0)
      if nb == 0
        @warn "SE: measurement $(m.id) excluded (its branch is not inside the island)"
        nothing
      else
        Measurement(typ = m.typ, value = m.value, sigma = m.sigma, active = m.active, branchIdx = nb, direction = m.direction, id = m.id)
      end
    elseif m.busIdx !== nothing
      Measurement(typ = m.typ, value = m.value, sigma = m.sigma, active = m.active, busIdx = sub.busmap[m.busIdx], direction = m.direction, id = m.id)
    else
      nothing
    end
    remapped === nothing && continue
    push!(meas, remapped)
    push!(keptIdx, k)
  end
  return meas, keptIdx
end

## island-wise estimation: per-island solve on the subnet, merged result
function _runse_islands!(net::Net, measurements::Vector{Measurement}, cfg::StateEstimationConfig, part)
  isempty(part.unassigned) || @warn "SE: $(length(part.unassigned)) measurement(s) reference isolated buses outside every island and are ignored"
  mergedV = buildVoltageVector(net)   # start values for unestimated islands
  nOrig = length(net.nodeVec)
  seVm = Float64[something(nd._vm_pu, 1.0) for nd in net.nodeVec]
  seVa = Float64[something(nd._va_deg, 0.0) for nd in net.nodeVec]
  sePinj = zeros(Float64, nOrig)
  seQinj = zeros(Float64, nOrig)
  islandsInfo = NamedTuple[]
  residuals = Float64[]
  Jsum = 0.0
  νsum = 0
  activeJsum = 0.0
  activeνsum = 0
  suppressedSum = 0
  anyActiveObj = false
  iters = 0
  conv = true
  shuntRows = NamedTuple[]
  robustRows = NamedTuple[]
  anyShunt = false
  anyRobust = false
  # tap estimates aggregate across islands: rows are remapped to parent-net
  # branch indices; the fixation J/dof sums are legitimate the same way the
  # island chi-squares add (independent runs)
  tapRows = NamedTuple[]
  anyTap = false
  tapJBefore = 0.0
  tapDofBefore = 0
  tapJAfter = 0.0
  tapDofAfter = 0
  tapAllFixed = true
  tapAnyOffgrid = false
  αdeg = nothing
  nα = 0
  estimated = 0
  for (iid, row) in enumerate(part.rows)
    idxs = part.perIsland[iid]
    refBus = row.n_ref > 0 ? 0 : row.chosen_ref_bus
    if isempty(idxs)
      push!(islandsInfo, (island = iid, n_bus = row.n_bus, estimated = false, reason = :no_measurements, converged = false, iterations = 0, objectiveJ = 0.0, dof = 0, band_reason = :none, z_wh = NaN, promoted_ref_bus = refBus, va_ref_offset_deg = nothing))
      continue
    end
    sub = _se_island_subnet(net, row)
    meas_i, keptIdx = _se_remap_island_measurements(sub, measurements, idxs)
    if isempty(meas_i)
      push!(islandsInfo, (island = iid, n_bus = row.n_bus, estimated = false, reason = :no_measurements, converged = false, iterations = 0, objectiveJ = 0.0, dof = 0, band_reason = :none, z_wh = NaN, promoted_ref_bus = refBus, va_ref_offset_deg = nothing))
      continue
    end
    # update_net on the subnet is the read-back channel for the merged
    # voltages (write-back into the parent happens below, gated by the
    # caller's update_net); shunt write-back is mirrored explicitly
    icfg = _copy_se_config_with(cfg; update_net = true)
    res = _runse_with_config!(sub.inet, meas_i, icfg)
    # the island run registered its chain start state under the subnet copy;
    # collect it for the merged parent-net registration below
    ist = _se_start_state(sub.inet)
    if ist !== nothing
      for (newb, oldb) in enumerate(row.buses)
        seVm[oldb] = ist.vm[newb]
        seVa[oldb] = ist.va[newb]
        sePinj[oldb] = ist.pinj[newb]
        seQinj[oldb] = ist.qinj[newb]
      end
    end
    estimated += 1
    conv &= res.converged
    iters = max(iters, res.iterations)
    Jsum += res.objectiveJ
    νsum += max(res.dof, 0)
    # J_active aggregates like the chi-squares: islands without suppression
    # contribute their plain J/dof, so the merged pair stays comparable
    if res.activeObjective !== nothing
      anyActiveObj = true
      activeJsum += res.activeObjective.j
      activeνsum += res.activeObjective.dof
      suppressedSum += res.activeObjective.suppressed
    else
      activeJsum += res.objectiveJ
      activeνsum += max(res.dof, 0)
    end
    append!(residuals, res.residuals)
    for (newb, oldb) in enumerate(row.buses)
      nd = sub.inet.nodeVec[newb]
      mergedV[oldb] = nd._vm_pu * cis(deg2rad(nd._va_deg))
    end
    if res.vaRefOffsetDeg !== nothing
      nα += 1
      αdeg = res.vaRefOffsetDeg
    end
    if res.shuntEstimates !== nothing
      anyShunt = true
      for r in res.shuntEstimates
        push!(shuntRows, (busIdx = row.buses[r.busIdx], busName = r.busName, B_model = r.B_model, B_est = r.B_est, delta = r.delta, frozen = r.frozen))
      end
      if cfg.update_shunts
        # the island solve wrote into the subnet's shunt copies; mirror the
        # susceptances back onto the parent net (order-preserving filter)
        for (pos, opos) in enumerate(sub.shuntpos)
          net.shuntVec[opos].y_pu_shunt = sub.inet.shuntVec[pos].y_pu_shunt
          net.shuntVec[opos].B_shunt = sub.inet.shuntVec[pos].B_shunt
        end
      end
    end
    if res.robustRows !== nothing
      anyRobust = true
      for r in res.robustRows
        push!(robustRows, (measurement_index = r.measurement_index >= 1 ? keptIdx[r.measurement_index] : 0, id = r.id, stage = r.stage, t = r.t, sigma_factor = r.sigma_factor))
      end
    end
    if res.tapEstimates !== nothing
      anyTap = true
      for tr in res.tapEstimates
        # subnet branch position k maps back through the island's ordered
        # branch list (row.branches[k] is the parent-net index)
        pb = 1 <= tr.branch <= length(row.branches) ? row.branches[tr.branch] : tr.branch
        push!(tapRows, (branch = pb, name = tr.name, mrid = tr.mrid, mode = tr.mode, alpha_deg = tr.alpha_deg, r1_est = tr.r1_est, r2_est = tr.r2_est, electrical_step_1 = tr.electrical_step_1, fixed_step_1 = tr.fixed_step_1, electrical_step_2 = tr.electrical_step_2, fixed_step_2 = tr.fixed_step_2, out_of_range = tr.out_of_range, fixed = tr.fixed, frozen_reason = tr.frozen_reason, island = iid))
      end
      if res.tapFixation !== nothing
        tapJBefore += res.tapFixation.j_before
        tapDofBefore += res.tapFixation.dof_before
        tapJAfter += res.tapFixation.j_after
        tapDofAfter += res.tapFixation.dof_after
        tapAllFixed &= res.tapFixation.fixed
        tapAnyOffgrid |= res.tapFixation.offgrid_residual
      end
      if cfg.update_taps
        # the island solve wrote the fixed positions into its subnet branch
        # copies; mirror them back onto the parent net (like update_shunts)
        for tr in res.tapEstimates
          (1 <= tr.branch <= length(row.branches)) || continue
          pb = row.branches[tr.branch]
          net.branchVec[pb].tap_ratio = sub.inet.branchVec[tr.branch].tap_ratio
          net.branchVec[pb].phase_shift_deg = sub.inet.branchVec[tr.branch].phase_shift_deg
        end
      end
    end
    # per-island band verdict: the summed chi-square test is legitimate
    # (independent chi-squares add), but a single bad island can hide in the
    # sum, so every island carries its own Wilson-Hilferty verdict here
    iv = _band_test_verdict(res.objectiveJ, res.dof)
    push!(islandsInfo, (island = iid, n_bus = row.n_bus, estimated = true, reason = :estimated, converged = res.converged, iterations = res.iterations, objectiveJ = res.objectiveJ, dof = res.dof, band_reason = iv.reason, z_wh = iv.z_wh, promoted_ref_bus = refBus, va_ref_offset_deg = res.vaRefOffsetDeg))
  end
  estimated == 0 && error("state estimation: no island carries any usable measurement")
  if cfg.update_net
    for info in islandsInfo
      info.estimated || continue
      for b in part.rows[info.island].buses
        nd = net.nodeVec[b]
        nd._vm_pu = abs(mergedV[b])
        nd._va_deg = rad2deg(angle(mergedV[b]))
      end
    end
    # merged chain start registration (runpf_from_se!/writeSEStateCSV read
    # it from the parent net); unestimated islands keep their start voltages
    # and a zero balance
    _register_se_start!(net, seVm, seVa, sePinj, seQinj, conv)
  end
  return SEResult(mergedV, conv, iters, norm(residuals), residuals, Jsum, νsum, _j_within_3sigma_band(Jsum, νsum), nα == 1 ? αdeg : nothing, anyShunt ? shuntRows : nothing, anyRobust ? robustRows : nothing, islandsInfo, anyTap ? tapRows : nothing, anyTap ? (fixed = tapAllFixed, j_before = tapJBefore, dof_before = tapDofBefore, j_after = tapJAfter, dof_after = tapDofAfter, offgrid_residual = tapAnyOffgrid) : nothing, nothing, anyActiveObj ? (j = activeJsum, dof = activeνsum, suppressed = suppressedSum) : nothing)
end

"""
    runse!(net; kwargs...) -> StateEstimationResult

Run the weighted-least-squares state estimation on the network's
measurements, island-wise, and write the estimated state back when
`updateNet` is set. The keyword form reads its defaults from the active
configuration.
"""
function runse!(net::Net, measurements::Vector{Measurement}, cfg::StateEstimationConfig)
  # stage-1 topology precheck (advisory): run BEFORE the partition,
  # warn per finding, attach the findings to the result. The estimation
  # itself always proceeds; a finding is information, not a gate.
  topo = nothing
  if cfg.topology_precheck
    pre = validate_topology(net, measurements; k_open = cfg.topology_open_flow_k, k_dead = cfg.topology_dead_flow_k, k_v = cfg.topology_voltage_k, k_kcl = cfg.topology_kcl_k)
    if !isempty(pre.findings)
      topo = pre.findings
      # one line per RUN, not per finding: a weak measurement set produces a
      # dozen findings, and the result object carries them all anyway
      shown = [string(f.kind, " at ", f.location, " (", f.evidence, ", ", f.severity, ")") for f in first(topo, 3)]
      @warn "runse!: topology precheck reported $(length(topo)) finding(s) (advisory, the estimation proceeds)" findings = join(shown, "; ") * (length(topo) > 3 ? "; ..." : "")
    end
  end
  part = _se_island_partition(net, measurements)
  res = length(part.rows) <= 1 ? _runse_with_config!(net, measurements, cfg) : _runse_islands!(net, measurements, cfg, part)
  topo === nothing && return res
  return SEResult(res.voltages, res.converged, res.iterations, res.residualNorm, res.residuals, res.objectiveJ, res.dof, res.jWithin3Sigma, res.vaRefOffsetDeg, res.shuntEstimates, res.robustRows, res.islands, res.tapEstimates, res.tapFixation, topo, res.activeObjective)
end

## per-island diagnostics merge (multi-island nets): every measured island is
## validated on its own subnet with its own reference; the merged report
## keeps the single-island shape (summed chi-square objective, rankings
## remapped into the caller's measurement indexing plus an `island` tag) and
## adds an `islands` vector with the per-island reports.
function _validate_measurements_islands(net::Net, measurements::Vector{Measurement}, part, vkw)
  islandReports = NamedTuple[]
  ranking = NamedTuple[]
  residuals = Float64[]
  rn = Float64[]
  wii = Float64[]
  stateVars = Float64[]
  robustRows = NamedTuple[]
  omegaPaths = Symbol[]
  gatingNotes = NamedTuple[]
  anyRobust = false
  Jsum = 0.0
  νsum = 0
  conv = true
  measured = 0
  for (iid, row) in enumerate(part.rows)
    idxs = part.perIsland[iid]
    if isempty(idxs)
      push!(islandReports, (island = iid, n_bus = row.n_bus, measured = false, report = nothing))
      continue
    end
    sub = _se_island_subnet(net, row)
    meas_i, keptIdx = _se_remap_island_measurements(sub, measurements, idxs)
    if isempty(meas_i)
      push!(islandReports, (island = iid, n_bus = row.n_bus, measured = false, report = nothing))
      continue
    end
    r = validate_measurements(sub.inet, meas_i; vkw...)
    measured += 1
    conv &= r.converged
    Jsum += r.objective.value
    νsum += max(r.objective.dof, 0)
    for rowr in r.measurement_ranking
      push!(ranking, merge(rowr, (measurement_index = rowr.measurement_index >= 1 ? keptIdx[rowr.measurement_index] : 0, island = iid)))
    end
    append!(residuals, r.residuals)
    append!(rn, r.normalized_residuals)
    append!(wii, r.residual_sensitivities)
    r.state_variances !== nothing && append!(stateVars, r.state_variances)
    push!(omegaPaths, r.omega_path)
    hasproperty(r, :gating_notes) && append!(gatingNotes, r.gating_notes)
    if r.robust_rows !== nothing
      anyRobust = true
      for rr in r.robust_rows
        push!(robustRows, merge(rr, (measurement_index = rr.measurement_index >= 1 ? keptIdx[rr.measurement_index] : 0,)))
      end
    end
    push!(islandReports, (island = iid, n_bus = row.n_bus, measured = true, report = r))
  end
  measured == 0 && error("validate_measurements: no active measurements")
  sort!(ranking; by = x -> x.abs_normalized_residual, rev = true)
  verdict = _band_test_verdict(Jsum, νsum)
  within = _j_within_3sigma_band(Jsum, νsum)
  suspicious = [row for row in ranking if row.suspicious]
  return (
    converged = conv,
    global_consistency = conv && within,
    objective = (value = Jsum, dof = νsum, zscore = verdict.z_legacy, z_wh = verdict.z_wh, within_3sigma = within, reason = verdict.reason, small_redundancy = verdict.small_redundancy),
    largest_normalized_residual = isempty(ranking) ? nothing : first(ranking),
    suspicious_measurements = suspicious,
    measurement_ranking = ranking,
    residuals = residuals,
    normalized_residuals = rn,
    residual_sensitivities = wii,
    correlation_enabled = vkw.reportResidualCorrelation,
    robust_rows = anyRobust ? robustRows : nothing,
    omega_path = length(unique(omegaPaths)) == 1 ? first(omegaPaths) : :mixed,
    state_variances = stateVars,
    gating_notes = gatingNotes,
    result = nothing,
    islands = islandReports,
  )
end

function _runse_with_config!(net::Net, measurements::Vector{Measurement}, cfg::StateEstimationConfig)
  maxIte = cfg.max_iter
  flatstart = cfg.flatstart
  jacEps = cfg.jac_eps
  # The Jacobian is built by forward differences with step `jacEps`, so its
  # error floor is O(jacEps): the state step cannot become smaller than that
  # noise, no matter how well the estimate has converged. Asking for a
  # tighter tolerance therefore cannot be met and the run burns its whole
  # iteration budget on a solution that already stands. Measured on a
  # 25000-bus set: tol 1e-8 failed after 30 iterations, tol 1e-6 succeeded
  # after 3 with J/dof 1.00 - the same estimate, only one of them can say so.
  tol = cfg.tol
  if tol < jacEps
    @info "state estimation: tolerance $(tol) is below the finite-difference noise floor (jac_eps = $(jacEps)) and cannot be reached; using $(jacEps)"
    tol = jacEps
  end
  updateNet = cfg.update_net
  # SE phase 3: contract closed-link clusters and remap the measurements; the
  # solve below runs entirely on prep.snet (== net when nothing is merged).
  prep = _se_prepare(net, measurements)
  onet = net
  net = prep.snet
  activeMeas = prep.meas
  activeOrig = prep.origIdx   # caller-vector positions (0 = LINKAGG aggregate)
  isempty(activeMeas) && error("runse!: no active measurements")

  # ImagMeas value gate: near zero current the derivative of
  # |I| is discontinuous, so a current measurement below 3 sigma is excluded
  # for the whole run, not just early iterations. An IaMeas row whose paired
  # magnitude row falls to this gate goes with it (its angle is just as
  # meaningless there); unpaired IaMeas rows are gated dynamically at the
  # predicted current (see the solve loop).
  if any(m -> m.typ == ImagMeas, activeMeas)
    droppedImagLoc = Set{Tuple{Union{Nothing,Int},Union{Nothing,Int},Symbol}}()
    keepGate = Int[]
    for (k, m) in enumerate(activeMeas)
      if m.typ == ImagMeas && m.value < 3.0 * m.sigma
        @info "runse!: current measurement $(m.id) excluded (value $(m.value) A below 3 sigma = $(3.0 * m.sigma) A)"
        push!(droppedImagLoc, (m.branchIdx, m.busIdx, m.direction))
      else
        push!(keepGate, k)
      end
    end
    activeMeas = activeMeas[keepGate]
    activeOrig = activeOrig[keepGate]
    if !isempty(droppedImagLoc)
      keepIa = Int[]
      for (k, m) in enumerate(activeMeas)
        if m.typ == IaMeas && (m.branchIdx, m.busIdx, m.direction) in droppedImagLoc
          @info "runse!: current-angle measurement $(m.id) excluded (:ia_below_current_floor, its paired magnitude failed the 3 sigma value gate)"
        else
          push!(keepIa, k)
        end
      end
      activeMeas = activeMeas[keepIa]
      activeOrig = activeOrig[keepIa]
    end
    isempty(activeMeas) && error("runse!: no active measurements left after the ImagMeas 3-sigma value gate")
  end
  # r0.9.10 scope guard: the estimator's measurement model does not carry
  # the one-sided branch reduction yet; a flow measurement on a partially
  # open branch would be evaluated against the wrong model. Rejected with a
  # clear message (follow-up tracked in the terminal-status review notes);
  # measurement generation never creates such measurements.
  for meas in activeMeas
    meas.branchIdx === nothing && continue
    (1 <= meas.branchIdx <= length(net.branchVec)) || continue
    st = _branch_terminal_state(net.branchVec[meas.branchIdx])
    if st == :open_from || st == :open_to
      error("runse!: measurement $(meas.id) references branch $(meas.branchIdx), which is open at one terminal ($(st)). Flow measurements on partially open branches are not supported yet; remove the measurement or close the terminal.")
    end
  end

  nbus = length(net.nodeVec)
  slackIdx = _find_slack_idx(net)
  Ybus = createYBUS(net = net)

  withVaOffset = _va_offset_active(activeMeas, cfg.pmu_ref_offset)

  # --- released shunt selection (SE phase 2, case A) ---------------------
  # A shunt marked via setShuntEstimation! becomes a B state only when at
  # least one active DIRECT shunt measurement references its bus; otherwise
  # it stays frozen at the model value (warning). Gated measurements (the
  # ImagMeas 3-sigma value gate above) no longer count as direct support.
  shuntCandidates = Int[]
  shuntFrozen = Int[]   # candidate shunts frozen (no measurement / not observable)
  for (si, sh) in enumerate(net.shuntVec)
    sh.estimate || continue
    if sh.status != 1
      @warn "runse!: shunt at bus $(sh.busIdx) is released for estimation but out of service; frozen at the model value"
      push!(shuntFrozen, si)
      continue
    end
    if sh.model != :Y
      # the injection-mode shunt is not stamped into Ybus, so the analytic
      # replacement below would double-count it; rejected in this phase
      error("runse!: shunt at bus $(sh.busIdx) uses the voltage-dependent injection mode ($(sh.model)); shunt estimation supports the default :Y (admittance) model only")
    end
    hasDirect = any(m -> (m.typ == ShuntQMeas && m.busIdx == sh.busIdx) || (m.typ == ImagMeas && m.branchIdx === nothing && m.busIdx == sh.busIdx), activeMeas)
    if !hasDirect
      @warn "runse!: shunt at bus $(sh.busIdx) is released for estimation but has no active direct shunt measurement (ShuntQMeas or bus-referenced ImagMeas); frozen at the model value"
      push!(shuntFrozen, si)
      continue
    end
    push!(shuntCandidates, si)
  end
  # link contraction can land several released shunts on one fused bus; one
  # B state per fused bus in this phase (their susceptances are not
  # separable from bus measurements alone)
  if length(shuntCandidates) > 1
    seen = Dict{Int,Int}()
    for si in shuntCandidates
      b = net.shuntVec[si].busIdx
      if haskey(seen, b)
        error("runse!: two shunts released for estimation resolve to the same (fused) bus $(b); estimate at most one shunt per link cluster in this phase")
      end
      seen[b] = si
    end
  end

  # Ybus handling: remove the released shunts' model admittance from the
  # diagonal ONCE; the prediction re-adds the injection analytically with the
  # state susceptance (see _predict_measurements). Ybus stays constant across
  # FD perturbations, only the analytic term varies.
  released = shuntCandidates
  shuntMap = nothing
  _shunt_col0 = (nbus - 1) + nbus + 1
  _build_shunt_state = function (rel::Vector{Int})
    for si in rel
      sh = net.shuntVec[si]
      Ybus[sh.busIdx, sh.busIdx] -= sh.y_pu_shunt
    end
    bInit = [imag(net.shuntVec[si].y_pu_shunt) for si in rel]
    xr = _initial_state_vector(net, slackIdx; flatstart = flatstart, withVaOffset = withVaOffset, shuntBInit = bInit)
    mapr = isempty(rel) ? nothing : ShuntStateMap(copy(rel), [net.shuntVec[si].busIdx for si in rel], _shunt_col0)
    return xr, mapr
  end
  x, shuntMap = _build_shunt_state(released)

  # Observability guard on the B columns: a released susceptance whose state
  # column no active measurement pins is frozen (state removed, model stamp
  # restored) and the run continues.
  if shuntMap !== nothing
    Htest, _ = _measurement_jacobian_fd(activeMeas, net, x, slackIdx, nbus, Ybus; eps = jacEps, withVaOffset = withVaOffset, shuntMap = shuntMap)
    keep = Int[]
    for (j, si) in enumerate(released)
      col = _shunt_col0 + j - 1
      # No try/catch here: an unmeasured column comes back with
      # numerical_observable = false, it does not raise. The ONLY exception
      # this call has is "state column out of bounds", which means a wrong
      # column index, and absorbing that would silently freeze every
      # released shunt while the run reported "not observable" (measured
      # 2026-09-06, task_silent_catch_v0101). Same class as the island map
      # that stayed empty in compareWithSV.
      obs = evaluate_local_observability_matrix(Htest, [col])
      if obs.numerical_observable
        push!(keep, si)
      else
        sh = net.shuntVec[si]
        Ybus[sh.busIdx, sh.busIdx] += sh.y_pu_shunt   # restore the model stamp
        push!(shuntFrozen, si)
        @warn "runse!: shunt B state at bus $(sh.busIdx) is not observable from the active measurements; frozen at the model value"
      end
    end
    if length(keep) != length(released)
      # rebuild without the frozen columns (the kept shunts' stamps are
      # already removed; remove them again after restoring? no: restore all
      # kept stamps first, then rebuild once from a clean Ybus state)
      for si in keep
        sh = net.shuntVec[si]
        Ybus[sh.busIdx, sh.busIdx] += sh.y_pu_shunt
      end
      released = keep
      x, shuntMap = _build_shunt_state(released)
    end
  end

  # Released transformer taps: state columns after the shunt B
  # states, before the trailing PMU offset. The stamped admittance terms
  # leave the Ybus once; the predictions re-add them at the cascade
  # position of the current state (see _predict_measurements).
  tapMap = _build_tap_state_map(net, (nbus - 1) + nbus + length(released) + 1)
  # frozen-regulator bookkeeping for the report: snet branch idx => reason
  # (:radial_no_voltage_pin structural, :not_observable numerical); fully
  # frozen transformers leave the map and get their model stamp back
  tapFrozen = Dict{Int,Symbol}()
  if tapMap !== nothing
    _tap_unstamp!(Ybus, net, tapMap)
    # structural guard: a bridge transformer whose cut-off side carries no
    # active voltage measurement makes the tap and the downstream voltage
    # indistinguishable; both regulators freeze at the current position
    gf1, gf2, greasons = _tap_guard_freezes(net, tapMap, activeMeas)
    for j in eachindex(greasons)
      greasons[j] == :ok && continue
      tapFrozen[tapMap.branch_idxs[j]] = greasons[j]
      @warn "runse!: released tap on transformer branch $(tapMap.branch_idxs[j]) ($(getCompName(net.branchVec[tapMap.branch_idxs[j]].comp))) is frozen: bridge transformer without a voltage measurement on its cut-off side (the tap would absorb the downstream voltage)"
    end
    (any(gf1) || any(gf2)) && (tapMap = _tap_apply_freezes(tapMap, gf1, gf2, net, Ybus))
  end
  if tapMap !== nothing
    # numerical observability guard on the remaining r columns, the same
    # per-column test the released shunt B states pass; a :both pair can
    # freeze partially (one regulator stays a state, the other at its r0)
    rinitG = _tap_r_init(tapMap)
    xg = withVaOffset ? vcat(x[1:(end-1)], rinitG, x[end]) : vcat(x, rinitG)
    Htest, _ = _measurement_jacobian_fd(activeMeas, net, xg, slackIdx, nbus, Ybus; eps = jacEps, withVaOffset = withVaOffset, shuntMap = shuntMap, tapMap = tapMap)
    nf1 = falses(length(tapMap.branch_idxs))
    nf2 = falses(length(tapMap.branch_idxs))
    _tap_col_dead = function (col::Int)
      # see the shunt guard above: an unmeasured column returns
      # numerical_observable = false, so no catch is needed; a wrong column
      # index must NOT be absorbed into "tap frozen"
      return !evaluate_local_observability_matrix(Htest, [col]).numerical_observable
    end
    for j in eachindex(tapMap.branch_idxs)
      tapMap.r1cols[j] != 0 && _tap_col_dead(tapMap.r1cols[j]) && (nf1[j] = true)
      tapMap.r2cols[j] != 0 && _tap_col_dead(tapMap.r2cols[j]) && (nf2[j] = true)
      if nf1[j] || nf2[j]
        tapFrozen[tapMap.branch_idxs[j]] = :not_observable
        @warn "runse!: tap state on transformer branch $(tapMap.branch_idxs[j]) ($(getCompName(net.branchVec[tapMap.branch_idxs[j]].comp))) is not observable from the active measurements; regulator frozen at the current position"
      end
    end
    (any(nf1) || any(nf2)) && (tapMap = _tap_apply_freezes(tapMap, nf1, nf2, net, Ybus))
  end
  if tapMap !== nothing
    rinit = _tap_r_init(tapMap)
    x = withVaOffset ? vcat(x[1:(end-1)], rinit, x[end]) : vcat(x, rinit)
  end
  _tap_ov(xv) = tapMap === nothing ? nothing : _tap_overlay(xv, tapMap, net)

  z = _measurement_vector(activeMeas)
  w = _weight_vector(activeMeas)

  converged = false
  r = zeros(Float64, length(activeMeas))
  iteDone = 0

  # Current iteration gate: current magnitudes AND current
  # angles stay out of the update step while ite < imag_activation_iteration,
  # so the first linearization(s) from a flat start run on power/voltage
  # measurements only. Row masking keeps measurement indices stable.
  imagGate = cfg.imag_activation_iteration
  nonImagRows = [i for (i, m) in enumerate(activeMeas) if !(m.typ in (ImagMeas, IaMeas))]
  hasImag = length(nonImagRows) < length(activeMeas)

  # IaMeas activity gate: a current ANGLE is meaningless near zero current.
  # Per iteration an Ia row participates only while the predicted current
  # magnitude passes 3 * sigma_I_ref, where sigma_I_ref is the sigma of the
  # paired ImagMeas at the same end (else the configured ampere floor).
  iaRows = [i for (i, m) in enumerate(activeMeas) if m.typ == IaMeas]
  iaThreshold = Dict{Int,Float64}()
  for i in iaRows
    mi = activeMeas[i]
    pair = findfirst(m -> m.typ == ImagMeas && m.branchIdx == mi.branchIdx && m.busIdx == mi.busIdx && m.direction == mi.direction, activeMeas)
    iaThreshold[i] = 3.0 * (pair === nothing ? cfg.ia_current_floor_A : activeMeas[pair].sigma)
  end
  # rows gated at the last evaluated state (drives the final J/dof exclusion)
  iaGatedNow = Set{Int}()
  function _ia_gated!(Vnow::Vector{ComplexF64})
    empty!(iaGatedNow)
    for i in iaRows
      ipred = abs(_end_current_A(net, activeMeas[i], Vnow, shuntMap === nothing ? nothing : _shunt_b_override(x, shuntMap), _tap_ov(x)))
      ipred < iaThreshold[i] && push!(iaGatedNow, i)
    end
    return iaGatedNow
  end

  # robust R modification (SE phase 4): active from robust_start_iteration,
  # modifies ONLY the solve weights of the current linearization. Every
  # statistic (objective, residuals, diagnostics) stays on the original w.
  # robustMode :staged is the two-stage k1/k2 modification (legacy Bool
  # `robust` maps here), :replacement the fixed-sigma suppression above
  # k_suppress.
  robustMode = _se_effective_robust_mode(cfg)
  robustOn = robustMode !== :off
  robustStart = cfg.robust_start_iteration
  robustK1 = cfg.robust_k1
  robustK2 = cfg.robust_k2
  robustKsup = cfg.k_suppress
  robustSupSigma = cfg.suppression_sigma

  # The WLS loop runs as a callable so the tap fixation stage below can
  # re-solve after the taps left the state vector. `gateImag = false` on the
  # fixation run: it starts from an already converged voltage state, so the
  # flat-start current gate (ite < imagGate) must not re-mask current rows.
  # replacement-mode weights: decided OUTSIDE the WLS loop on a converged
  # state (see the suppression rounds below the first solve) and constant
  # during one solve. The closure reassigns this enclosing local on purpose.
  # Released regulator states are bounded by the changer's DECLARED
  # mechanical band. Without it the Gauss-Newton step can drive r1 toward
  # -1, where the cascade t = t_base/((1+r1)(1+r2 e^{j alpha})) is singular:
  # measured 2026-09-06 on a CGMES delivery, one weakly determined
  # transformer ran to -329 electrical steps on a band of about -14 to +18
  # and took the whole estimation with it. The clamp is not a tuning knob,
  # a position outside the band does not exist on the mechanical grid.
  # Hits are counted per column: a regulator that keeps pressing against its
  # bound is a finding (the measurements want a position the changer cannot
  # reach), not something to hide.
  tapStepLimits = tapMap === nothing ? Dict{Int,Float64}() : _tap_step_limits(tapMap, net)
  tapClampHits = Dict{Int,Int}()
  local replacementFrozenW = nothing
  # The Jacobian of the LAST iteration, kept for the suppression rounds:
  # they judge on the normalized residual and therefore need the residual
  # covariance diagonal, which needs H. Re-assembling it would pay a second
  # finite-difference pass (the long pole on large systems, 58 s of a 127 s
  # run on case13659pegase); the last iteration's H sits at the same state
  # up to one step below tol, which the smooth Omega does not notice.
  # The closure reassigns this enclosing local on purpose.
  local lastJacobian = nothing

  function _wls_solve!(gateImag::Bool)
    converged = false
    # robust weight freeze: the stage assignment can toggle a row across a
    # stage boundary every iteration and keep the state jittering above
    # tol (a weight limit cycle, seen with many released taps). Once the
    # step is small the stages are frozen, and the now-constant weights
    # let the smooth fixed point converge.
    local robustFrozenW = nothing
    local lastStep = Inf
    # FD column coloring (task_se_fd_coloring): the coloring comes from
    # the STRUCTURAL pattern (topological couplings), not from an
    # assembled H, because point patterns never settled on large systems
    # (micro-derivatives at machine precision appear and vanish bitwise
    # between iterations; measured on case13659pegase, where both a
    # two-point union and a convergence-gated union kept falling back).
    # Engaged from iteration 1 on large systems only (threshold comment
    # at the constant); the in-assembly detector remains the correctness
    # guard and falls back loudly if a coupling were missing.
    local fdColoring = length(x) >= _SE_FD_COLORING_MIN_STATES ? _fd_structural_coloring(activeMeas, length(x), slackIdx, nbus, Ybus, withVaOffset, shuntMap, tapMap, net) : nothing
    for ite = 1:maxIte
      # one abort check per iteration: a Web UI run that the user gave up
      # on stops here instead of finishing its 20 iterations (the Jacobian
      # assembly below is the long pole, so this is the useful place)
      sparlectra_check_abort()
      H, h, _ = _measurement_jacobian_fd(activeMeas, net, x, slackIdx, nbus, Ybus; eps = jacEps, withVaOffset = withVaOffset, shuntMap = shuntMap, tapMap = tapMap, coloring = fdColoring)
      lastJacobian = H
      r = z - h
      _wrap_angle_residuals!(r, activeMeas)

      wSolve = w
      if robustMode === :replacement
        # replacement weights are judged on a CONVERGED state in the
        # suppression rounds outside this loop (a per-iteration threshold
        # suppressed healthy tight rows during the flat-start transient
        # and the estimate drifted); during one solve they are constant
        replacementFrozenW !== nothing && (wSolve = replacementFrozenW)
      elseif robustOn && ite >= robustStart && robustFrozenW !== nothing
        wSolve = robustFrozenW
      elseif robustOn && ite >= robustStart
        wSolve = copy(w)
        for i in eachindex(activeMeas)
          stage, σmod = _robust_stage(abs(r[i]), activeMeas[i].sigma, robustK1, robustK2)
          stage == 0 && continue
          wSolve[i] = inv(σmod * σmod)
        end
        # freeze when the state has settled, or unconditionally after 15
        # robust iterations: a weight limit cycle with a larger amplitude
        # would otherwise never converge, and constant weights turn the
        # remaining iterations into a plain, smoothly converging WLS
        (lastStep < 1e-4 || ite >= robustStart + 15) && (robustFrozenW = wSolve)
      end

      # Normal equations for WLS:
      #   Δx = (H'WH)^{-1} H'W r
      # with W = diag(w). Sparse throughout (task_se_sparse): H is a
      # SparseMatrixCSC, so G = H' * (W * H) stays sparse and the solve
      # below takes the sparse factorization path; the former dense
      # HW-broadcast form would have densified through the adjoint
      # broadcast.
      if gateImag && hasImag && ite < imagGate
        Hs = H[nonImagRows, :]
        ws = wSolve[nonImagRows]
        G = Hs' * (Diagonal(ws) * Hs)
        g = Hs' * (ws .* r[nonImagRows])
      elseif !isempty(iaRows) && !isempty(_ia_gated!(_state_to_voltage(x, slackIdx, nbus)))
        useRows = [i for i in eachindex(activeMeas) if !(i in iaGatedNow)]
        Hs = H[useRows, :]
        ws = wSolve[useRows]
        G = Hs' * (Diagonal(ws) * Hs)
        g = Hs' * (ws .* r[useRows])
      else
        G = H' * (Diagonal(wSolve) * H)
        g = H' * (wSolve .* r)
      end

      # svd_max_n raised for the WLS normal equations: released extra
      # states (shunt B, tap r) can turn out COLLECTIVELY dependent even
      # when each column passes its individual observability test (e.g.
      # parallel released taps between the same stations); a singular G
      # then degrades to a least-squares step instead of aborting the run
      # weak prior on the tap states of PARALLEL released transformers
      # (full Tikhonov term, gradient included): parallel taps between the
      # same buses share an unobservable difference direction (only their
      # SUM is observable). Without the prior that direction explodes under
      # iteration-varying robust weights (reproduced on a 17-trafo case) or
      # creeps forever when damped on the matrix alone. The prior "stay
      # near the initial position" gives it a real fixed point. It is
      # applied ONLY where the degeneracy exists: a solitary released tap
      # stays prior-free, so the exact-recovery roundtrips are untouched.
      if tapMap !== nothing && tapMap.ncols > 0
        pairCount = Dict{Tuple{Int,Int},Int}()
        for k in tapMap.branch_idxs
          br = net.branchVec[k]
          key = minmax(Int(br.fromBus), Int(br.toBus))
          pairCount[key] = get(pairCount, key, 0) + 1
        end
        # a pair member with its own active flow/current row is individually
        # identifiable (branch flows follow the LIVE tap since the
        # branchFlow_pu fix), and the prior would swamp that small
        # difference-direction information with a penalty scaled by the
        # large sum-direction diagonal, biasing the split away from the
        # measurement optimum (seen as J stuck at 0.6 instead of 0). Keep
        # the prior ONLY for pairs no flow/current measurement can tell
        # apart.
        measuredPairs = Set{Tuple{Int,Int}}()
        for k in tapMap.branch_idxs
          any(m -> m.branchIdx == k && m.typ in (PflowMeas, QflowMeas, ImagMeas), activeMeas) || continue
          br = net.branchVec[k]
          push!(measuredPairs, minmax(Int(br.fromBus), Int(br.toBus)))
        end
        for j in eachindex(tapMap.branch_idxs)
          br = net.branchVec[tapMap.branch_idxs[j]]
          key = minmax(Int(br.fromBus), Int(br.toBus))
          pairCount[key] > 1 || continue
          key in measuredPairs && continue
          for c0 in (tapMap.r1cols[j], tapMap.r2cols[j])
            c0 == 0 && continue
            λ = 1e-3 * abs(G[c0, c0])
            G[c0, c0] += λ
            g[c0] += λ * (rinit[c0 - tapMap.col0 + 1] - x[c0])
          end
        end
      end
      Δx = solve_linear(G, g; allow_pinv = true, svd_max_n = 20_000)
      # Released regulator states get a STEP LIMIT, not a hard band.
      # Measured 2026-09-06 on a CGMES delivery with 11 transformers off by
      # 4 mechanical steps: clamping the state into its declared band (both
      # per component and as an active-set projection) rescued the one
      # runaway transformer but made five that had converged stick to the
      # lower bound and diverge. The transient excursion is part of how the
      # tap state finds its position; blocking it breaks the path. What has
      # to be bounded is how far one iteration may move, so a weakly
      # determined regulator cannot reach the singularity of
      # t = t_base/((1+r1)(1+r2 e^{j alpha})) in a single step.
      # The whole step is scaled by one factor, so voltage and tap states
      # keep moving along one consistent direction.
      if !isempty(tapStepLimits)
        α = 1.0
        for (col, lim) in tapStepLimits
          col <= length(Δx) || continue
          d = abs(Δx[col])
          d <= lim && continue
          f = lim / d
          f < α && (α = f)
          tapClampHits[col] = get(tapClampHits, col, 0) + 1
        end
        α < 1.0 && (Δx .*= α)
      end
      x .+= Δx

      iteDone = ite
      lastStep = norm(Δx, Inf)
      if lastStep < tol
        converged = true
        break
      end
    end
  end
  _wls_solve!(true)

  # --- replacement suppression rounds (robust_mode :replacement) -----------
  # Judged on the CONVERGED residuals, never on flat-start transients: solve
  # plain first, put every real row whose NORMALIZED residual reaches
  # k_suppress on the fixed suppression sigma, re-solve with those weights
  # frozen, and repeat until the suppression set is stable (bounded rounds).
  # Virtual rows (sigma <= 1e-6) are never suppressed; statistics keep the
  # original w.
  #
  # One scale for both decisions (task_se_bad_data_v0100). Until 0.10.0 this
  # used the raw ratio |r_i|/sigma_i while the elimination used
  # rn = r_i/sqrt(Omega_ii). Those differ exactly where it matters: with
  # wii = Omega_ii * w_i, a row at the localizability guideline wii = 0.3 has
  # sqrt(Omega_ii) = 0.55*sigma, so its true normalized residual is about 1.8
  # times the raw ratio, and for a nearly critical row the raw ratio stays
  # small while rn grows without bound. The suppression was therefore weakest
  # precisely where a gross error does the most damage.
  #
  # The Omega diagonal comes from the SAME diagnostics path the elimination
  # uses: one `_residual_diagnostics` call per round on the Jacobian of the
  # last WLS iteration (Takahashi selected inverse above the state
  # threshold, dense below). Its guard carries over for free: a row whose
  # Omega_ii sits at the numerical floor is not localizable, gets rn = 0
  # there, and can therefore never be suppressed on the strength of a large
  # raw residual.
  if robustMode === :replacement && converged
    prevSet = Set{Int}()
    for _ = 1:3
      Vr = _state_to_voltage(x, slackIdx, nbus)
      hr = _predict_measurements(activeMeas, net, Vr, Ybus; vaOffsetRad = _va_offset_from_state(x, withVaOffset), shuntB = shuntMap === nothing ? nothing : _shunt_b_override(x, shuntMap), tapOverlay = _tap_ov(x))
      rr = z - hr
      _wrap_angle_residuals!(rr, activeMeas)
      supDiag = _suppression_normalized_residuals(lastJacobian, rr, w)
      # No usable Omega (no Jacobian kept, or the diagnostics refused this
      # size): suppress NOTHING rather than fall back to the raw ratio the
      # task removed. A missed suppression is a weaker estimate; a wrong one
      # on the wrong scale is a silently wrong estimate.
      supDiag === nothing && break
      rnSup = supDiag.rn
      wiiSup = supDiag.wii
      # Only rows the residual can actually point at. `wii` is the share of
      # the row's own error that reaches its residual; below the literature
      # guideline 0.3 the row is nearly critical, its residual is
      # structurally small and its rn correspondingly inflated. Suppressing
      # such a row does not remove bad data, it removes the little
      # information the row still carries, and the neighbour that shared its
      # redundancy then looks like gross error: measured on the PST warm-up
      # case, suppressing Qinj_5 (wii 0.05, raw residual 1.4 sigma, rn 5.8)
      # drove its partner row to rn 58 and got that HEALTHY row eliminated.
      # The elimination already reports this quantity as `localizable`; the
      # suppression now honours the same bound.
      supSet = Set{Int}(i for i in eachindex(activeMeas) if !_elimination_protected(activeMeas[i]) && abs(rnSup[i]) >= robustKsup && wiiSup[i] > _SE_SUPPRESSION_MIN_WII)
      # say how many rows this round takes out and how far the worst one is
      # past the limit: a suppression set that suddenly counts thousands is
      # a scale or observability problem, not bad data, and nothing else in
      # the run would show it
      if !isempty(supSet)
        worst = maximum(abs(rnSup[i]) for i in supSet)
        @info "state estimation: down-weighting $(length(supSet)) of $(count(i -> !_elimination_protected(activeMeas[i]), eachindex(activeMeas))) eligible row(s) (normalized residual >= $(robustKsup), largest $(round(worst; digits = 2)))"
      end
      supSet == prevSet && break
      prevSet = supSet
      wRepl = copy(w)
      for i in supSet
        wRepl[i] = inv(robustSupSigma * robustSupSigma)
      end
      replacementFrozenW = wRepl
      _wls_solve!(false)
      converged || break
    end
  end

  # --- tap fixation run (always follows the released estimate) --------------
  # After the estimation the released taps are FIXED to the nearest
  # MECHANICAL step and one more run is solved in which the tap is no state
  # any more (it leaves J and the dof). J before versus after the fixation
  # is reported in `tapFixation`; the SEResult voltages/objectiveJ/dof are
  # those of the FIXED final run.
  tapRows = nothing
  tapJBefore = NaN
  tapDofBefore = 0
  tapIteEst = 0
  tapConvEst = false
  tapFixed = false
  if tapMap !== nothing
    # J of the CONTINUOUS solution, same statistics discipline as the final
    # block below (Ia rows gated at this state leave J and dof)
    Vc = _state_to_voltage(x, slackIdx, nbus)
    hC = _predict_measurements(activeMeas, net, Vc, Ybus; vaOffsetRad = _va_offset_from_state(x, withVaOffset), shuntB = shuntMap === nothing ? nothing : _shunt_b_override(x, shuntMap), tapOverlay = _tap_ov(x))
    rC = z - hC
    _wrap_angle_residuals!(rC, activeMeas)
    gatedC = isempty(iaRows) ? Set{Int}() : copy(_ia_gated!(Vc))
    tapJBefore = isempty(gatedC) ? _wls_objective(rC, w) : _wls_objective(rC[[i for i in eachindex(rC) if !(i in gatedC)]], w[[i for i in eachindex(w) if !(i in gatedC)]])
    tapDofBefore = length(activeMeas) - length(gatedC) - length(x)
    tapIteEst = iteDone
    tapConvEst = converged
    tapRows, tapMapFixed, xFixed = _fixate_taps(x, tapMap, net)
    if converged
      # the fixed positions replace the tap states; re-solve from the
      # converged voltages (never from a diverged state: the rounded steps
      # would be meaningless, so a diverged run keeps the continuous result)
      x = xFixed
      tapMap = tapMapFixed
      tapFixed = true
      _wls_solve!(false)
      converged = converged && tapConvEst
    end
  end

  Vest = _state_to_voltage(x, slackIdx, nbus)
  vaOffsetRad = _va_offset_from_state(x, withVaOffset)
  hfinal = _predict_measurements(activeMeas, net, Vest, Ybus; vaOffsetRad = vaOffsetRad, shuntB = shuntMap === nothing ? nothing : _shunt_b_override(x, shuntMap), tapOverlay = _tap_ov(x))
  r = z - hfinal
  _wrap_angle_residuals!(r, activeMeas)
  # Ia rows gated at the FINAL state never entered the solve there and must
  # not enter the statistics either: excluded from J and from the redundancy
  finalGated = isempty(iaRows) ? Set{Int}() : _ia_gated!(Vest)
  for i in finalGated
    @info "runse!: current-angle measurement $(activeMeas[i].id) excluded (:ia_below_current_floor, predicted current below $(round(iaThreshold[i]; digits = 2)) A)"
  end
  jval = isempty(finalGated) ? _wls_objective(r, w) : _wls_objective(r[[i for i in eachindex(r) if !(i in finalGated)]], w[[i for i in eachindex(w) if !(i in finalGated)]])
  ν = length(activeMeas) - length(finalGated) - length(x)
  # J_active: the objective WITHOUT the replacement-suppressed rows. Their
  # frozen solve weight (suppression sigma) already removed them from the
  # STATE; the honest jval above keeps them at original sigmas as the alarm
  # signal, J_active answers "how well does the model fit the data the
  # estimator actually trusted". The band test stays on the honest jval.
  activeObjective = nothing
  if replacementFrozenW !== nothing
    supIdx = Set{Int}(i for i in eachindex(w) if replacementFrozenW[i] != w[i] && !(i in finalGated))
    if !isempty(supIdx)
      keep = [i for i in eachindex(r) if !(i in supIdx) && !(i in finalGated)]
      activeObjective = (j = _wls_objective(r[keep], w[keep]), dof = ν - length(supIdx), suppressed = length(supIdx))
    end
  end

  if updateNet
    # write back to the ORIGINAL net: every original bus reads the estimated
    # voltage of its cluster representative (reps), translated through the
    # dense renumbering (busmap). Identity mappings when nothing was merged.
    Vorig = Vector{ComplexF64}(undef, length(onet.nodeVec))
    for i in eachindex(onet.nodeVec)
      s = get(prep.busmap, prep.reps[i], 0)
      # a bus whose representative did not survive contraction (branchless
      # pocket) keeps its stored state; nothing was estimated for it
      Vorig[i] = s == 0 ? complex(something(onet.nodeVec[i]._vm_pu, 1.0), 0.0) : Vest[s]
      s == 0 && continue
      onet.nodeVec[i]._vm_pu = abs(Vest[s])
      onet.nodeVec[i]._va_deg = rad2deg(angle(Vest[s]))
    end
    calcNetLosses!(onet, Vorig)
    # refresh shunt powers so the post-SE link allocation (and reporting)
    # reads consistent node._pShunt/_qShunt values (phase-3 gap closure)
    updateShuntPowers!(net = onet)

    # SE chain start registration (phase 5): keep the estimated voltages and
    # nodal balances for runpf_from_se!/writeSEStateCSV. The balances are the
    # h(x) bus injections of the estimated state (network Ybus view, released
    # shunt states included); on a contracted net the representative carries
    # the cluster sum and the members carry 0, so a re-contracting power flow
    # sees exactly the estimated balance.
    SbusEst = calc_injections(Ybus, Vest) .* net.baseMVA
    if shuntMap !== nothing
      for (j, si) in enumerate(shuntMap.shunt_idxs)
        sh = net.shuntVec[si]
        bset = x[shuntMap.col0+j-1]
        SbusEst[sh.busIdx] += abs2(Vest[sh.busIdx]) * conj(complex(real(sh.y_pu_shunt), bset)) * net.baseMVA
      end
    end
    if tapMap !== nothing
      # released trafos are unstamped from the Ybus; re-add their terminal
      # injections at the estimated cascade position (same as the prediction)
      for (k, ov) in _tap_ov(x)
        brk = net.branchVec[k]
        fb = Int(brk.fromBus)
        tb = Int(brk.toBus)
        SbusEst[fb] += Vest[fb] * conj(ov.Y11 * Vest[fb] + ov.Y12 * Vest[tb]) * net.baseMVA
        SbusEst[tb] += Vest[tb] * conj(ov.Y21 * Vest[fb] + ov.Y22 * Vest[tb]) * net.baseMVA
      end
    end
    nOrig = length(onet.nodeVec)
    seVm = Vector{Float64}(undef, nOrig)
    seVa = Vector{Float64}(undef, nOrig)
    sePinj = zeros(Float64, nOrig)
    seQinj = zeros(Float64, nOrig)
    for i = 1:nOrig
      seVm[i] = something(onet.nodeVec[i]._vm_pu, 1.0)
      seVa[i] = something(onet.nodeVec[i]._va_deg, 0.0)
      s = get(prep.busmap, prep.reps[i], 0)
      if s != 0 && prep.reps[i] == i
        sePinj[i] = real(SbusEst[s])
        seQinj[i] = imag(SbusEst[s])
      end
    end
    _register_se_start!(onet, seVm, seVa, sePinj, seQinj, converged)
  end

  # --- shunt estimate report and optional write-back (SE phase 2) --------
  # bus indices in the report and the write-back target both live in the
  # ORIGINAL net: prep.buses translates snet indices back, prep.shuntOrig
  # maps snet shunt positions onto net.shuntVec (identity without merges).
  shuntEstimates = nothing
  if !isempty(released) || !isempty(shuntFrozen)
    name_by_idx = _bus_name_by_idx(net)
    rows = NamedTuple[]
    for (j, si) in enumerate(released)
      sh = net.shuntVec[si]
      origBus = prep.buses[sh.busIdx]
      bModel = imag(sh.y_pu_shunt)
      bEst = x[_shunt_col0+j-1]
      push!(rows, (busIdx = origBus, busName = get(name_by_idx, sh.busIdx, string(origBus)), B_model = bModel, B_est = bEst, delta = bEst - bModel, frozen = false))
    end
    for si in sort(shuntFrozen)
      sh = net.shuntVec[si]
      origBus = prep.buses[sh.busIdx]
      bModel = imag(sh.y_pu_shunt)
      push!(rows, (busIdx = origBus, busName = get(name_by_idx, sh.busIdx, string(origBus)), B_model = bModel, B_est = bModel, delta = 0.0, frozen = true))
    end
    shuntEstimates = rows
    # estimation must never silently overwrite model data: write-back only on
    # explicit request, and only from a converged run
    if cfg.update_shunts && converged
      for (j, si) in enumerate(released)
        osh = onet.shuntVec[prep.shuntOrig[si]]
        bEst = x[_shunt_col0+j-1]
        osh.y_pu_shunt = complex(real(osh.y_pu_shunt), bEst)
        osh.B_shunt = bEst
      end
    end
  end

  # robust reporting (SE phase 4): stages at the FINAL residuals; only rows
  # that left stage 0 are listed. origIdx 0 marks a LINKAGG aggregate.
  robustRows = nothing
  if robustOn
    rows = NamedTuple[]
    if robustMode === :replacement
      # report the set that was ACTUALLY suppressed (the frozen replacement
      # weights), not a recomputation at the final state
      if replacementFrozenW !== nothing
        for i in eachindex(activeMeas)
          replacementFrozenW[i] == w[i] && continue
          σ = activeMeas[i].sigma
          push!(rows, (measurement_index = activeOrig[i], id = activeMeas[i].id, stage = 3, t = abs(r[i]) / σ, sigma_factor = robustSupSigma / σ))
        end
      end
    else
      for i in eachindex(activeMeas)
        σ = activeMeas[i].sigma
        stage, σmod = _robust_stage(abs(r[i]), σ, robustK1, robustK2)
        stage == 0 && continue
        push!(rows, (measurement_index = activeOrig[i], id = activeMeas[i].id, stage = stage, t = abs(r[i]) / σ, sigma_factor = σmod / σ))
      end
    end
    robustRows = rows
  end

  # --- tap estimate report -----------------------------------------------
  # branch indices in the report live in the ORIGINAL net: component
  # identity (cID) survives the contraction copies, so match by it
  # (identity mapping without merges). mRIDs are only populated for CGMES
  # cases (net.cgmes_ids); MATPOWER/DTF rows carry an empty mrid and are
  # addressed by branch index and bus numbers.
  tapEstimates = nothing
  tapFixation = nothing
  if tapRows !== nothing || !isempty(tapFrozen)
    mrids = _transformer_mrids(onet)
    cidToOrig = Dict{String,Int}(br.comp.cID => k for (k, br) in enumerate(onet.branchVec))
    rows = NamedTuple[]
    liveIdxs = Set{Int}()
    if tapRows !== nothing
      for tr in tapRows
        push!(liveIdxs, tr.branch)
        ob = get(cidToOrig, net.branchVec[tr.branch].comp.cID, tr.branch)
        push!(rows, (branch = ob, name = tr.name, mrid = get(mrids, ob, ""), mode = tr.mode, alpha_deg = tr.alpha_deg, r1_est = tr.r1_est, r2_est = tr.r2_est, electrical_step_1 = tr.electrical_step_1, fixed_step_1 = tr.fixed_step_1, electrical_step_2 = tr.electrical_step_2, fixed_step_2 = tr.fixed_step_2, out_of_range = tr.out_of_range, fixed = tapFixed, frozen_reason = get(tapFrozen, tr.branch, :none)))
      end
    end
    # transformers whose regulators are ALL frozen left the map; report them
    # at their current position so a released-but-guarded tap never vanishes
    # silently from the result
    for (k, reason) in sort(collect(tapFrozen); by = first)
      k in liveIdxs && continue
      br = net.branchVec[k]
      r10, r20 = _tap_r0_split(br, br.tap_est_mode, br.tap_est_alpha_deg)
      step1 = br.tap_step > 0.0 ? br.tap_step : 0.00625
      s2 = calcSkewAngleTap(tap_fraction = r20, skew_angle_deg = br.tap_est_alpha_deg)
      step2 = br.phase_step_deg > 0.0 ? br.phase_step_deg : 1.25
      ob = get(cidToOrig, br.comp.cID, k)
      push!(rows, (branch = ob, name = getCompName(br.comp), mrid = get(mrids, ob, ""), mode = br.tap_est_mode, alpha_deg = br.tap_est_alpha_deg, r1_est = r10, r2_est = r20, electrical_step_1 = r10 / step1, fixed_step_1 = Int(round(r10 / step1, RoundNearestTiesAway)), electrical_step_2 = s2.effective_shift_deg / step2, fixed_step_2 = Int(round(s2.effective_shift_deg / step2, RoundNearestTiesAway)), out_of_range = false, fixed = false, frozen_reason = reason))
    end
    tapEstimates = rows
    # write-back of the FIXED mechanical positions: only on explicit request
    # (update_taps, default false), only from a converged fixed run, and a
    # frozen regulator keeps its exact model value inside the cascade (the
    # same never-silently-overwrite protection as update_shunts)
    if cfg.update_taps && tapFixed && tapMap !== nothing
      for j in eachindex(tapMap.branch_idxs)
        br = net.branchVec[tapMap.branch_idxs[j]]
        ct = _cascade_tap(br, tapMap.r1_0[j], tapMap.r2_0[j], tapMap.alphas[j])
        ob = get(cidToOrig, br.comp.cID, 0)
        ob == 0 && continue
        onet.branchVec[ob].tap_ratio = ct.tap_ratio
        onet.branchVec[ob].phase_shift_deg = ct.phase_shift_deg
      end
    end
    if tapRows !== nothing
      # :offgrid_tap_residual: a band failure that only
      # APPEARS through the fixation is not bad data. The continuous fit
      # (taps still states) passed or undershot the chi-square band, the
      # fixed run fails it :high, so the J jump comes from the mechanical
      # rounding: the true position sits between steps, or the step table
      # (tap_step/neutral position) is wrong.
      offgrid = tapFixed && _band_test_verdict(tapJBefore, tapDofBefore).reason in (:ok, :low) && _band_test_verdict(jval, ν).reason == :high
      tapFixation = (fixed = tapFixed, j_before = tapJBefore, dof_before = tapDofBefore, j_after = jval, dof_after = ν, iterations_estimation = tapIteEst, iterations_fixation = tapFixed ? iteDone : 0, offgrid_residual = offgrid)
      if tapFixed
        @info "runse!: tap fixation J $(round(tapJBefore; sigdigits = 4)) (dof $(tapDofBefore)) -> J $(round(jval; sigdigits = 4)) (dof $(ν))"
        offgrid && @info "runse!: :offgrid_tap_residual, the J jump comes from the fixation, not from bad data (true position between mechanical steps, or a wrong step table)"
      else
        @warn "runse!: tap estimation did not converge; the taps were NOT fixed and the continuous result is returned"
      end
    end
  end

  return SEResult(Vest, converged, iteDone, norm(r), r, jval, ν, _j_within_3sigma_band(jval, ν), withVaOffset ? rad2deg(vaOffsetRad) : nothing, shuntEstimates, robustRows, nothing, tapEstimates, tapFixation, nothing, activeObjective)
end

function runse!(
  net::Net,
  measurements::Vector{Measurement};
  maxIte::Int = state_estimation_config().max_iter,
  tol::Float64 = state_estimation_config().tol,
  flatstart::Bool = state_estimation_config().flatstart,
  jacEps::Float64 = state_estimation_config().jac_eps,
  updateNet::Bool = state_estimation_config().update_net,
  pmuRefOffset::Symbol = state_estimation_config().pmu_ref_offset,
  imagActivationIteration::Int = state_estimation_config().imag_activation_iteration,
  updateShunts::Bool = state_estimation_config().update_shunts,
  updateTaps::Bool = state_estimation_config().update_taps,
  topologyPrecheck::Bool = state_estimation_config().topology_precheck,
  robust::Bool = state_estimation_config().robust,
  robustStartIteration::Int = state_estimation_config().robust_start_iteration,
  robustMode::Symbol = state_estimation_config().robust_mode,
  robustK1::Float64 = state_estimation_config().robust_k1,
  robustK2::Float64 = state_estimation_config().robust_k2,
  kSuppress::Float64 = state_estimation_config().k_suppress,
  suppressionSigma::Float64 = state_estimation_config().suppression_sigma,
)
  cfg = StateEstimationConfig(max_iter = maxIte, tol = tol, flatstart = flatstart, jac_eps = jacEps, update_net = updateNet, pmu_ref_offset = pmuRefOffset, imag_activation_iteration = imagActivationIteration, update_shunts = updateShunts, update_taps = updateTaps, topology_precheck = topologyPrecheck, robust = robust, robust_start_iteration = robustStartIteration, robust_mode = robustMode, robust_k1 = robustK1, robust_k2 = robustK2, k_suppress = kSuppress, suppression_sigma = suppressionSigma)
  return runse!(net, measurements, cfg)
end

function runse!(
  net::Net;
  maxIte::Int = state_estimation_config().max_iter,
  tol::Float64 = state_estimation_config().tol,
  flatstart::Bool = state_estimation_config().flatstart,
  jacEps::Float64 = state_estimation_config().jac_eps,
  updateNet::Bool = state_estimation_config().update_net,
  pmuRefOffset::Symbol = state_estimation_config().pmu_ref_offset,
  imagActivationIteration::Int = state_estimation_config().imag_activation_iteration,
  updateShunts::Bool = state_estimation_config().update_shunts,
  updateTaps::Bool = state_estimation_config().update_taps,
  topologyPrecheck::Bool = state_estimation_config().topology_precheck,
  robust::Bool = state_estimation_config().robust,
  robustStartIteration::Int = state_estimation_config().robust_start_iteration,
  robustMode::Symbol = state_estimation_config().robust_mode,
  robustK1::Float64 = state_estimation_config().robust_k1,
  robustK2::Float64 = state_estimation_config().robust_k2,
  kSuppress::Float64 = state_estimation_config().k_suppress,
  suppressionSigma::Float64 = state_estimation_config().suppression_sigma,
)
  cfg = StateEstimationConfig(max_iter = maxIte, tol = tol, flatstart = flatstart, jac_eps = jacEps, update_net = updateNet, pmu_ref_offset = pmuRefOffset, imag_activation_iteration = imagActivationIteration, update_shunts = updateShunts, update_taps = updateTaps, topology_precheck = topologyPrecheck, robust = robust, robust_start_iteration = robustStartIteration, robust_mode = robustMode, robust_k1 = robustK1, robust_k2 = robustK2, k_suppress = kSuppress, suppression_sigma = suppressionSigma)
  return runse!(net, Measurement[m for m in net.measurements], cfg)
end

"""
    validate_measurements(net, measurements; kwargs...) -> NamedTuple

Run state-estimation diagnostics on currently active measurements and return a
machine-readable report with:
- global bad-data consistency check (`global_consistency`)
- χ²-like objective plausibility summary
- largest-normalized-residual ranking with the residual-sensitivity diagonal
  `wii` and a `localizable` flag (`wii > wiiThreshold`, default 0.3, the
  established literature threshold)
- suspicious measurement list (threshold-based)
- optional residual-correlation columns (`reportResidualCorrelation`,
  default from `state_estimation.report_residual_correlation`): per row the
  maximum |k_ij| over all partners plus a warning above 1/sqrt(2) (the
  result 7, simple-redundant group)

Current-magnitude measurements below their 3 sigma value gate are excluded,
matching the estimator's own activation rule, so the report describes the
measurement set the estimator actually used. With active `ImagMeas` rows the
residual-sensitivity and correlation figures become load-flow dependent
they hold for the estimated operating point.

`measurement_index` in the ranking rows is the position in the caller's
measurement vector, with one sentinel: `0` marks a `LINKAGG` cluster
aggregate (SE phase 3) that has no single source row. Check `>= 1` before
indexing the measurement vector with it.
"""
function validate_measurements(
  net::Net,
  measurements::Vector{Measurement};
  maxIte::Int = 12,
  tol::Float64 = 1e-6,
  flatstart::Bool = true,
  jacEps::Float64 = 1e-6,
  normalizedThreshold::Float64 = state_estimation_config().k_eliminate,
  pmuRefOffset::Symbol = state_estimation_config().pmu_ref_offset,
  wiiThreshold::Float64 = 0.3,
  reportResidualCorrelation::Bool = state_estimation_config().report_residual_correlation,
  robust::Bool = state_estimation_config().robust,
  robustStartIteration::Int = state_estimation_config().robust_start_iteration,
  robustMode::Symbol = state_estimation_config().robust_mode,
  robustK1::Float64 = state_estimation_config().robust_k1,
  robustK2::Float64 = state_estimation_config().robust_k2,
  kSuppress::Float64 = state_estimation_config().k_suppress,
  suppressionSigma::Float64 = state_estimation_config().suppression_sigma,
)
  # island-wise nets: run the full diagnostics per measured island and merge
  # (measurement_index stays the position in the CALLER's vector, so the
  # sequential elimination in runse_diagnostics works unchanged)
  part = _se_island_partition(net, measurements)
  if length(part.rows) > 1
    vkw = (maxIte = maxIte, tol = tol, flatstart = flatstart, jacEps = jacEps, normalizedThreshold = normalizedThreshold, pmuRefOffset = pmuRefOffset, wiiThreshold = wiiThreshold, reportResidualCorrelation = reportResidualCorrelation, robust = robust, robustStartIteration = robustStartIteration, robustMode = robustMode, robustK1 = robustK1, robustK2 = robustK2, kSuppress = kSuppress, suppressionSigma = suppressionSigma)
    return _validate_measurements_islands(net, measurements, part, vkw)
  end

  # SE phase 3: the diagnostics run on the same contracted net and remapped
  # measurement set as the estimator (measurement_index in the rankings stays
  # the position in the CALLER's vector; 0 marks a LINKAGG aggregate row)
  prep = _se_prepare(net, measurements)
  net = prep.snet
  activeMeas = prep.meas
  activeIdx = prep.origIdx
  # mirror the runse! ImagMeas 3-sigma value gate so the diagnostics see the
  # same measurement set the estimator solves on; IaMeas rows follow their
  # dropped magnitude pair (same rule as the estimator)
  droppedImagLoc = Set{Tuple{Union{Nothing,Int},Union{Nothing,Int},Symbol}}()
  for m in activeMeas
    (m.typ == ImagMeas && m.value < 3.0 * m.sigma) && push!(droppedImagLoc, (m.branchIdx, m.busIdx, m.direction))
  end
  gatingNotes = NamedTuple[]
  keep = Int[]
  for k in eachindex(activeMeas)
    m = activeMeas[k]
    if m.typ == ImagMeas && m.value < 3.0 * m.sigma
      push!(gatingNotes, (reason = :imag_below_value_gate, id = m.id))
    elseif m.typ == IaMeas && (m.branchIdx, m.busIdx, m.direction) in droppedImagLoc
      push!(gatingNotes, (reason = :ia_below_current_floor, id = m.id))
    else
      push!(keep, k)
    end
  end
  activeMeas = activeMeas[keep]
  activeIdx = activeIdx[keep]
  isempty(activeMeas) && error("validate_measurements: no active measurements")

  netLocal = deepcopy(net)
  # topologyPrecheck off: the diagnostics re-run this solve per elimination
  # step, and the caller (or the service) runs the precheck exactly once
  result = runse!(netLocal, activeMeas; maxIte = maxIte, tol = tol, flatstart = flatstart, jacEps = jacEps, updateNet = false, pmuRefOffset = pmuRefOffset, robust = robust, robustStartIteration = robustStartIteration, robustMode = robustMode, robustK1 = robustK1, robustK2 = robustK2, kSuppress = kSuppress, suppressionSigma = suppressionSigma, topologyPrecheck = false)

  # dynamic part of the IaMeas gate: rows below the predicted-current
  # threshold at the ESTIMATED state leave the diagnostic set too
  if any(m -> m.typ == IaMeas, activeMeas)
    cfg0 = state_estimation_config()
    keep2 = Int[]
    for k in eachindex(activeMeas)
      m = activeMeas[k]
      if m.typ == IaMeas
        pair = findfirst(mm -> mm.typ == ImagMeas && mm.branchIdx == m.branchIdx && mm.busIdx == m.busIdx && mm.direction == m.direction, activeMeas)
        thr = 3.0 * (pair === nothing ? cfg0.ia_current_floor_A : activeMeas[pair].sigma)
        ipred = abs(_end_current_A(netLocal, m, result.voltages, nothing))
        if ipred < thr
          push!(gatingNotes, (reason = :ia_below_current_floor, id = m.id))
          continue
        end
      end
      push!(keep2, k)
    end
    activeMeas = activeMeas[keep2]
    activeIdx = activeIdx[keep2]
    isempty(activeMeas) && error("validate_measurements: no active measurements left after the IaMeas gate")
  end

  nbus = length(netLocal.nodeVec)
  slackIdx = _find_slack_idx(netLocal)
  withVaOffset = result.vaRefOffsetDeg !== nothing

  # SE phase 4 (diagnostics unification): the diagnostics Jacobian uses the
  # SAME extended state definition as the estimator. Released shunts (the
  # non-frozen shuntEstimates rows, order = the estimator's B-column order)
  # re-enter as B states at their estimated values; frozen shunts stay at the
  # model value, exactly as in the estimator.
  releasedRows = result.shuntEstimates === nothing ? NamedTuple[] : NamedTuple[rw for rw in result.shuntEstimates if !rw.frozen]
  shuntBInit = Float64[rw.B_est for rw in releasedRows]
  x = _initial_state_vector(netLocal, slackIdx; flatstart = flatstart, withVaOffset = withVaOffset, shuntBInit = shuntBInit)
  nθ = nbus - 1
  p = 0
  for i in eachindex(netLocal.nodeVec)
    if i != slackIdx
      p += 1
      x[p] = angle(result.voltages[i])
    end
    x[nθ+i] = abs(result.voltages[i])
  end
  if withVaOffset
    x[end] = deg2rad(result.vaRefOffsetDeg)
  end

  # tap unification: the estimator's final voltages belong to the
  # FIXED tap positions, but write-back is usually off, so the diagnostics
  # copy still carries the model taps; stamping the fixed cascade here keeps
  # the residuals free of a fake tap discrepancy at the released trafos.
  # Frozen regulators (and unfixed runs) stay at the model position; a
  # partially frozen :both pair conservatively stays at the model position
  # too (the row cannot say which regulator froze).
  if result.tapEstimates !== nothing
    for tr in result.tapEstimates
      (1 <= tr.branch <= length(netLocal.branchVec)) || continue
      (tr.fixed && tr.frozen_reason == :none) || continue
      br = netLocal.branchVec[tr.branch]
      br.ratio == 0.0 && continue
      step1 = br.tap_step > 0.0 ? br.tap_step : 0.00625
      step2 = br.phase_step_deg > 0.0 ? br.phase_step_deg : 1.25
      r1 = tr.mode in (:ratio, :both) ? tr.fixed_step_1 * step1 : 0.0
      r2 = 0.0
      if tr.mode in (:pst, :both)
        if br.phase_du_step > 0.0
          # Delta-u PST: the fixed step counts additional-voltage steps,
          # r2 is linear on that grid (the degree formula below would stamp
          # a wrong angle and fake an 800-MW discrepancy in the ranking)
          r2 = tr.fixed_step_2 * br.phase_du_step
        else
          ϕ = -deg2rad(tr.fixed_step_2 * step2)
          r2 = abs(sin(deg2rad(tr.alpha_deg) - ϕ)) < 1e-12 ? 0.0 : sin(ϕ) / sin(deg2rad(tr.alpha_deg) - ϕ)
        end
      end
      ct = _cascade_tap(br, r1, r2, tr.alpha_deg)
      br.tap_ratio = ct.tap_ratio
      br.phase_shift_deg = ct.phase_shift_deg
    end
  end

  Ybus = createYBUS(net = netLocal)
  shuntMap = nothing
  if !isempty(releasedRows)
    sidxs = Int[]
    sbuses = Int[]
    for rw in releasedRows
      # shuntEstimates carry ORIGINAL bus indices; translate through the
      # contraction (identity when nothing is merged)
      snetBus = prep.busmap[prep.reps[rw.busIdx]]
      si = netLocal.shuntDict[snetBus]
      push!(sidxs, si)
      push!(sbuses, snetBus)
      # remove the model stamp; the prediction re-adds the injection with the
      # state susceptance (same contract as the estimator)
      sh = netLocal.shuntVec[si]
      Ybus[snetBus, snetBus] -= sh.y_pu_shunt
    end
    shuntMap = ShuntStateMap(sidxs, sbuses, nθ + nbus + 1)
  end

  H, h = _measurement_jacobian_fd(activeMeas, netLocal, x, slackIdx, nbus, Ybus; eps = jacEps, withVaOffset = withVaOffset, shuntMap = shuntMap)
  z = _measurement_vector(activeMeas)
  w = _weight_vector(activeMeas)
  r = z - h
  _wrap_angle_residuals!(r, activeMeas)
  # K needs the full Omega, so a correlation request pins the dense path;
  # otherwise the Takahashi selected inverse takes over above the size
  # threshold (with an automatic dense fallback on any guard failure)
  rd = _residual_diagnostics(H, r, w; need_full_omega = reportResidualCorrelation)
  kmax = reportResidualCorrelation && rd.omega !== nothing ? _residual_correlation_max(rd.omega) : nothing
  ranking = _build_measurement_suspicion_report(activeMeas, activeIdx, r, rd.rn; normalizedThreshold = normalizedThreshold, wii = rd.wii, wiiThreshold = wiiThreshold, kmax = kmax)

  suspicious = [row for row in ranking if row.suspicious]
  verdict = _band_test_verdict(result.objectiveJ, result.dof)

  return (
    converged = result.converged,
    global_consistency = result.converged && result.jWithin3Sigma,
    objective = (value = result.objectiveJ, dof = result.dof, zscore = verdict.z_legacy, z_wh = verdict.z_wh, within_3sigma = result.jWithin3Sigma, reason = verdict.reason, small_redundancy = verdict.small_redundancy),
    largest_normalized_residual = isempty(ranking) ? nothing : first(ranking),
    suspicious_measurements = suspicious,
    measurement_ranking = ranking,
    residuals = r,
    normalized_residuals = rd.rn,
    residual_sensitivities = rd.wii,
    correlation_enabled = reportResidualCorrelation,
    robust_rows = result.robustRows,
    omega_path = rd.omega_path,
    state_variances = rd.state_variances,
    gating_notes = gatingNotes,
    result = result,
  )
end

function validate_measurements(net::Net; kwargs...)
  return validate_measurements(net, Measurement[m for m in net.measurements]; kwargs...)
end

## Internal helper: measurements protected from bad-data elimination.
## Zero-injection pseudo-measurements (id prefix "ZI") encode network
## structure, not telemetry, and any near-exact measurement (sigma <= 1e-6)
## is a constraint in disguise; eliminating either would remove model
## knowledge instead of a faulty meter. Derived shunt pseudo-measurements
## (id prefix "SHDERIV", SE phase 2 case B) are protected too: eliminating
## one would silently mask its source measurements (bay current, Vm).
@inline function _elimination_protected(m::Measurement)::Bool
  return startswith(m.id, "ZI") || startswith(m.id, "SHDERIV") || m.sigma <= 1e-6
end

"""
    runse_diagnostics(net, measurements; max_eliminations=3, kwargs...) -> NamedTuple

Extended diagnostics workflow around `validate_measurements` with sequential
bad-data elimination (identification stage 1):

While the χ²-like 3σ band test fails, a suspicious measurement
(`|rn| >= normalizedThreshold`) exists, and fewer than `max_eliminations`
measurements have been removed, the top suspicious measurement is
deactivated and the diagnostics rerun. Zero-injection pseudo-measurements
(id prefix `ZI`) and near-exact measurements (`sigma <= 1e-6`) are never
eliminated.

Returns `(diagnostics, rerun, eliminations, stop_reason, final_diagnostics)`:
- `diagnostics`: the initial `validate_measurements` report
- `eliminations`: one trace row per elimination
  `(elimination, measurement_index, id, typ, normalized_residual_before,
  wii, skipped_unlocalizable, objective_before, objective_after)`.
  `wii` is the localizability of the removed row and
  `skipped_unlocalizable` counts the suspects that were passed over in that
  round BECAUSE they are not localizable, so the trace shows why the
  largest normalized residual was not the one removed.
- `stop_reason`: `:consistent` (band test passed), `:no_suspicious_left`,
  `:no_localizable_suspect` (suspects remain, but none of them is
  localizable: the data cannot say which row is wrong, which points at a
  measurement gap around those rows rather than at bad data),
  `:max_eliminations`, or `:not_converged`
- `final_diagnostics`: the report after the last elimination (equal to
  `diagnostics` when nothing was eliminated)
- `rerun`: backward-compatible view of the FIRST elimination
  (`deactivated_measurement_index`, `deactivated_measurement_id`,
  `diagnostics`), `nothing` when no elimination happened

Backward compatibility: `deactivate_and_rerun = true` is kept as an alias for
`max_eliminations = 1` (exactly one elimination pass, result also in `rerun`).
Unlike the old single pass, the loop only eliminates while the band test
fails; a suspicious measurement in an already consistent set is reported but
not removed.
"""
function runse_diagnostics(
  net::Net,
  measurements::Vector{Measurement};
  deactivate_and_rerun::Bool = false,
  max_eliminations::Union{Nothing,Int} = nothing,
  maxIte::Int = 12,
  tol::Float64 = 1e-6,
  flatstart::Bool = true,
  jacEps::Float64 = 1e-6,
  normalizedThreshold::Float64 = state_estimation_config().k_eliminate,
  pmuRefOffset::Symbol = state_estimation_config().pmu_ref_offset,
  wiiThreshold::Float64 = 0.3,
  reportResidualCorrelation::Bool = state_estimation_config().report_residual_correlation,
  robust::Bool = state_estimation_config().robust,
  robustStartIteration::Int = state_estimation_config().robust_start_iteration,
  robustMode::Symbol = state_estimation_config().robust_mode,
  robustK1::Float64 = state_estimation_config().robust_k1,
  robustK2::Float64 = state_estimation_config().robust_k2,
  kSuppress::Float64 = state_estimation_config().k_suppress,
  suppressionSigma::Float64 = state_estimation_config().suppression_sigma,
)
  # alias handling: an explicit max_eliminations wins; deactivate_and_rerun
  # = true alone maps to exactly one elimination (the pre-0.10 behavior)
  effMax = something(max_eliminations, deactivate_and_rerun ? 1 : 3)
  effMax >= 0 || error("runse_diagnostics: max_eliminations must be >= 0")

  vkw = (maxIte = maxIte, tol = tol, flatstart = flatstart, jacEps = jacEps, normalizedThreshold = normalizedThreshold, pmuRefOffset = pmuRefOffset, wiiThreshold = wiiThreshold, reportResidualCorrelation = reportResidualCorrelation, robust = robust, robustStartIteration = robustStartIteration, robustMode = robustMode, robustK1 = robustK1, robustK2 = robustK2, kSuppress = kSuppress, suppressionSigma = suppressionSigma)
  base = validate_measurements(net, measurements; vkw...)

  meas2 = copy(measurements)
  current = base
  trace = NamedTuple[]
  firstElimination = nothing
  stopReason = :consistent

  while true
    if !current.converged
      stopReason = :not_converged
      break
    end
    if current.objective.within_3sigma
      stopReason = :consistent
      break
    end
    # candidates: suspicious, not elimination-protected (ranking is sorted
    # by |rn| descending, so the first candidate is the top suspect).
    # measurement_index 0 marks a LINKAGG cluster aggregate (SE phase 3):
    # it has no single source row to deactivate and is never eliminated.
    # Localizability is a PRECONDITION for elimination, the same one and the
    # same threshold the replacement suppression uses (`wii > wiiThreshold`,
    # the value the report already prints as `localizable`). Removing a row
    # the residual cannot point at costs observability and buys nothing: a
    # critical row has a structurally zero residual, so a gross error there
    # is invisible and the row can never earn its removal. Until now the
    # SHARPER action (removal) had the WEAKER precondition than the milder
    # one (down-weighting), which is backwards.
    suspects = [row for row in current.measurement_ranking if row.suspicious && row.measurement_index >= 1 && !_elimination_protected(meas2[row.measurement_index])]
    candidates = [row for row in suspects if !isnan(row.wii) && row.wii > wiiThreshold]
    skipped_unlocalizable = length(suspects) - length(candidates)
    if isempty(candidates)
      # Not the same thing as "nothing suspicious left": if suspects remain
      # but none of them is localizable, the data cannot say WHICH row is
      # wrong. That points at a measurement gap around those rows, not at
      # bad data, and it gets its own stop reason instead of a silent stop.
      stopReason = isempty(suspects) ? :no_suspicious_left : :no_localizable_suspect
      break
    end
    if length(trace) >= effMax
      stopReason = :max_eliminations
      break
    end

    target = first(candidates)
    idx = target.measurement_index
    objBefore = current.objective.value
    meas2[idx] = _set_measurement_active(meas2[idx], false)
    current = validate_measurements(net, meas2; vkw...)
    push!(trace, (elimination = length(trace) + 1, measurement_index = idx, id = target.id, typ = target.typ,
                  normalized_residual_before = target.normalized_residual, wii = target.wii,
                  skipped_unlocalizable = skipped_unlocalizable,
                  objective_before = objBefore, objective_after = current.objective.value))
    if firstElimination === nothing
      firstElimination = (deactivated_measurement_index = idx, deactivated_measurement_id = measurements[idx].id, diagnostics = current)
    end
  end

  # stage-2 topology classification: the FINGERPRINT of a
  # topology error is an EXHAUSTED elimination that still fails the band
  # test :high with the surviving suspects clustered at one station. A
  # curable gross error never shows it (the elimination succeeds), so a
  # finding requires the full fingerprint, never a part of it.
  topoFindings = nothing
  if stopReason == :max_eliminations && current.converged && hasproperty(current.objective, :reason) && current.objective.reason == :high
    reps, mems = _topology_station_map(net)
    clusterMin = state_estimation_config().topology_cluster_min
    # station -> surviving suspects (a branch row counts for both terminals)
    byStation = Dict{Int,Vector{NamedTuple}}()
    for row in current.measurement_ranking
      (row.suspicious && row.measurement_index >= 1) || continue
      m = meas2[row.measurement_index]
      stations = Int[]
      m.busIdx !== nothing && push!(stations, reps[m.busIdx])
      if m.branchIdx !== nothing && 1 <= m.branchIdx <= length(net.branchVec)
        br = net.branchVec[m.branchIdx]
        push!(stations, reps[Int(br.fromBus)], reps[Int(br.toBus)])
      end
      for st in unique(stations)
        push!(get!(byStation, st, NamedTuple[]), row)
      end
    end
    # stage-1 agreement (open design decision 2): a precheck finding at the
    # same station is extra evidence, noted but not severity-changing
    pre = validate_topology(net, measurements)
    preStations = Set{Int}()
    for f in pre.findings
      mloc = match(r"^(branch|link) (\d+)", f.location)
      if mloc !== nothing
        idx = parse(Int, mloc.captures[2])
        if mloc.captures[1] == "branch" && 1 <= idx <= length(net.branchVec)
          push!(preStations, reps[Int(net.branchVec[idx].fromBus)], reps[Int(net.branchVec[idx].toBus)])
        elseif mloc.captures[1] == "link" && 1 <= idx <= length(net.linkVec)
          push!(preStations, reps[Int(net.linkVec[idx].fromBus)], reps[Int(net.linkVec[idx].toBus)])
        end
      end
    end
    rows = NamedTuple[]
    # A total order, so the strongest station comes first on every Julia
    # version: most suspects, then the largest normalized residual among
    # them, then the label. Sorting by count alone left ties in Dict
    # iteration order, which changed between Julia 1.12 and 1.13 and moved a
    # different station to the front (found on the CI run of 2026-09-11).
    station_rank(x) = (-length(last(x)), -maximum(r.abs_normalized_residual for r in last(x)), _topology_station_label(net, first(x), mems))
    for (st, suspects) in sort(collect(byStation); by = station_rank)
      length(suspects) >= clusterMin || continue
      notes = Symbol[]
      st in preStations && push!(notes, :precheck_agreement)
      # K-matrix side signal: with correlation columns present, a cluster
      # whose members all correlate above 1/sqrt(2) is one correlated group
      if current.correlation_enabled && all(!isnan(r.max_abs_correlation) && r.max_abs_correlation > 1.0 / sqrt(2.0) for r in suspects)
        push!(notes, :correlated_group)
      end
      push!(rows, (stage = :classification, kind = :topology_error_suspected_at_station, station = st, location = _topology_station_label(net, st, mems), evidence = join(("$(r.id) (|rn| $(round(r.abs_normalized_residual; digits = 1)))" for r in suspects), ", "), severity = :strong, notes = notes))
    end
    isempty(rows) || (topoFindings = rows)
  end

  return (diagnostics = base, rerun = firstElimination, eliminations = trace, stop_reason = stopReason, final_diagnostics = current, topology_findings = topoFindings)
end

function runse_diagnostics(net::Net; kwargs...)
  return runse_diagnostics(net, Measurement[m for m in net.measurements]; kwargs...)
end

@inline function _base_diag_report(diag)
  return hasproperty(diag, :diagnostics) ? diag.diagnostics : diag
end

"""
    summarize_se_diagnostics(diag) -> NamedTuple

Create a compact interpretation summary for a diagnostics object returned by
`validate_measurements` or `runse_diagnostics`.

`global_consistency` is interpreted as:
- `true`: SE converged and objective is inside χ²-like 3σ plausibility band
- `false`: either non-convergence or implausibly large objective
"""
function summarize_se_diagnostics(diag)
  base = _base_diag_report(diag)
  suspicious_count = length(base.suspicious_measurements)
  total = length(base.measurement_ranking)

  # band-test reason (SE phase 4, Wilson-Hilferty): :high and :low are the
  # two failure directions, :no_redundancy the nu = 0 special case. Reports
  # from before the verdict fields default to the legacy :high wording.
  bandReason = hasproperty(base.objective, :reason) ? base.objective.reason : :high
  reason = if !base.converged
    "State estimation did not converge."
  elseif base.objective.within_3sigma
    "No global inconsistency detected."
  elseif bandReason == :low
    "Objective J is implausibly small (Wilson-Hilferty z below -3): the measurement sigmas are likely overestimated, the typical signature of noise-free synthetic data."
  elseif bandReason == :no_redundancy
    "No redundancy (nu = 0): the residuals are structurally zero and bad data is invisible."
  else
    "Objective J is outside the χ²-like 3σ plausibility band (possible bad data/model mismatch)."
  end

  small_redundancy = hasproperty(base.objective, :small_redundancy) ? base.objective.small_redundancy : false

  return (global_consistency = base.global_consistency, converged = base.converged, objective = base.objective, suspicious_count = suspicious_count, total_measurements = total, reason = reason, small_redundancy = small_redundancy)
end

@inline function _global_consistency_label(summary)
  return summary.global_consistency ? "PASS (globally plausible)" : "FAIL (check objective/residuals)"
end

"""
    print_se_diagnostics(io, diag; topN=10)    

Pretty-print diagnostics from `validate_measurements` or `runse_diagnostics`
including:
- explanation of `global_consistency`
- tabular measurement ranking
- BAD/OK marker per measurement
- optional rerun comparison if present
"""
function print_se_diagnostics(diag; io = stdout, topN::Int = 10, format::Symbol = :plain)
  base = _base_diag_report(diag)
  summary = summarize_se_diagnostics(diag)

  nrows = min(topN, length(base.measurement_ranking))
  # wii/localizable columns exist on 0.10 ranking rows; correlation columns
  # carry values only when the K-matrix report was enabled
  hasWii = !isempty(base.measurement_ranking) && hasproperty(first(base.measurement_ranking), :wii)
  showK = hasproperty(base, :correlation_enabled) && base.correlation_enabled
  if format == :markdown
    println(io, "## State-estimation diagnostics")
    println(io)
    println(io, "- **Converged:** $(summary.converged)")
    println(io, "- **Global consistency:** $(summary.global_consistency) — $(_global_consistency_label(summary))")
    @printf(io, "- **Objective J:** %.6f (dof=%d, z=%.3f, within_3sigma=%s)\n", summary.objective.value, summary.objective.dof, summary.objective.zscore, string(summary.objective.within_3sigma))
    if hasproperty(summary.objective, :z_wh)
      @printf(io, "- **Wilson-Hilferty:** z_wh=%.3f, reason=%s%s\n", summary.objective.z_wh, string(summary.objective.reason), summary.small_redundancy ? " (small redundancy: nu < 30, reduced test power)" : "")
    end
    println(io, "- **Interpretation:** $(summary.reason)")
    println(io, "- **Suspicious measurements:** $(summary.suspicious_count) / $(summary.total_measurements)")
    # island-wise runs: per-island J/dof plus Wilson-Hilferty, so a single
    # bad island cannot hide inside the summed chi-square
    if hasproperty(base, :islands) && base.islands !== nothing
      println(io)
      println(io, "### Islands")
      println(io, "| Island | Buses | Measured | J | dof | z_wh | Band |")
      println(io, "|---:|---:|:---:|---:|---:|---:|:---|")
      for irow in base.islands
        if irow.measured && irow.report !== nothing
          o = irow.report.objective
          @printf(io, "| %d | %d | yes | %.4f | %d | %.3f | %s |\n", irow.island, irow.n_bus, o.value, o.dof, o.z_wh, string(o.reason))
        else
          println(io, "| $(irow.island) | $(irow.n_bus) | no | | | | skipped (no measurements) |")
        end
      end
    end
    println(io)
    println(io, "### Measurement ranking (largest |normalized residual|)")
    if hasWii
      if showK
        println(io, "| Idx | ID | Type | Residual | Norm.Res | wii | Loc | MaxK | Flag |")
        println(io, "|---:|:---|:---|---:|---:|---:|:---:|---:|:---:|")
      else
        println(io, "| Idx | ID | Type | Residual | Norm.Res | wii | Loc | Flag |")
        println(io, "|---:|:---|:---|---:|---:|---:|:---:|:---:|")
      end
    else
      println(io, "| Idx | ID | Type | Residual | Norm.Res | Flag |")
      println(io, "|---:|:---|:---|---:|---:|:---:|")
    end
    for k = 1:nrows
      row = base.measurement_ranking[k]
      flag = row.suspicious ? "BAD" : "OK"
      if hasWii
        loc = row.localizable ? "yes" : "no"
        if showK
          kwarn = row.correlation_warning ? " (!)" : ""
          @printf(io, "| %d | %s | %s | %.5f | %.5f | %.3f | %s | %.3f%s | %s |\n", row.measurement_index, row.id, string(row.typ), row.residual, row.normalized_residual, row.wii, loc, row.max_abs_correlation, kwarn, flag)
        else
          @printf(io, "| %d | %s | %s | %.5f | %.5f | %.3f | %s | %s |\n", row.measurement_index, row.id, string(row.typ), row.residual, row.normalized_residual, row.wii, loc, flag)
        end
      else
        @printf(io, "| %d | %s | %s | %.5f | %.5f | %s |\n", row.measurement_index, row.id, string(row.typ), row.residual, row.normalized_residual, flag)
      end
    end
  else
    println(io, "State-estimation diagnostics")
    println(io, "--------------------------------------------------------------")
    @printf(io, "Converged: %s\n", string(summary.converged))
    @printf(io, "Global consistency: %s (%s)\n", string(summary.global_consistency), _global_consistency_label(summary))
    @printf(io, "Objective J: %.6f (dof=%d, z=%.3f, within_3sigma=%s)\n", summary.objective.value, summary.objective.dof, summary.objective.zscore, string(summary.objective.within_3sigma))
    if hasproperty(summary.objective, :z_wh)
      @printf(io, "Wilson-Hilferty: z_wh=%.3f, reason=%s%s\n", summary.objective.z_wh, string(summary.objective.reason), summary.small_redundancy ? " (small redundancy: nu < 30, reduced test power)" : "")
    end
    hasproperty(base, :omega_path) && @printf(io, "Omega path: %s\n", string(base.omega_path))
    @printf(io, "Interpretation: %s\n", summary.reason)
    @printf(io, "Suspicious measurements: %d / %d\n", summary.suspicious_count, summary.total_measurements)
    if hasproperty(base, :gating_notes) && !isempty(base.gating_notes)
      println(io, "Gated measurements (excluded from the solve and the statistics):")
      for gn in base.gating_notes
        println(io, "  ", gn.id, " (", gn.reason, ")")
      end
    end
    # island-wise runs: per-island band test so no island hides in the sum
    if hasproperty(base, :islands) && base.islands !== nothing
      println(io, "\nIslands")
      for irow in base.islands
        if irow.measured && irow.report !== nothing
          o = irow.report.objective
          @printf(io, "  island %d (%d buses): J=%.4f dof=%d z_wh=%.3f reason=%s\n", irow.island, irow.n_bus, o.value, o.dof, o.z_wh, string(o.reason))
        else
          @printf(io, "  island %d (%d buses): skipped (no measurements)\n", irow.island, irow.n_bus)
        end
      end
    end

    println(io, "\nMeasurement ranking (largest |normalized residual|)")
    println(io, "----------------------------------------------------------------------------------------------")
    if hasWii
      if showK
        @printf(io, "%5s %-28s %-12s %12s %12s %8s %5s %8s %8s\n", "Idx", "ID", "Type", "Residual", "Norm.Res", "wii", "Loc", "MaxK", "Flag")
      else
        @printf(io, "%5s %-28s %-12s %12s %12s %8s %5s %8s\n", "Idx", "ID", "Type", "Residual", "Norm.Res", "wii", "Loc", "Flag")
      end
    else
      @printf(io, "%5s %-28s %-12s %12s %12s %8s\n", "Idx", "ID", "Type", "Residual", "Norm.Res", "Flag")
    end
    println(io, "----------------------------------------------------------------------------------------------")
    for k = 1:nrows
      row = base.measurement_ranking[k]
      flag = row.suspicious ? "BAD" : "OK"
      if hasWii
        loc = row.localizable ? "yes" : "no"
        if showK
          kwarn = row.correlation_warning ? "!" : " "
          @printf(io, "%5d %-28s %-12s %12.5f %12.5f %8.3f %5s %7.3f%s %8s\n", row.measurement_index, row.id, string(row.typ), row.residual, row.normalized_residual, row.wii, loc, row.max_abs_correlation, kwarn, flag)
        else
          @printf(io, "%5d %-28s %-12s %12.5f %12.5f %8.3f %5s %8s\n", row.measurement_index, row.id, string(row.typ), row.residual, row.normalized_residual, row.wii, loc, flag)
        end
      else
        @printf(io, "%5d %-28s %-12s %12.5f %12.5f %8s\n", row.measurement_index, row.id, string(row.typ), row.residual, row.normalized_residual, flag)
      end
    end
  end

  # robust R-modification rows (SE phase 4): measurements the FINAL iteration
  # weighted down; the diagnosis itself always runs on original sigmas
  if hasproperty(base, :robust_rows) && base.robust_rows !== nothing && !isempty(base.robust_rows)
    if format == :markdown
      println(io)
      println(io, "### Robust R modification (final iteration, solve weights only)")
      println(io, "| Idx | ID | Stage | t | sigma factor |")
      println(io, "|---:|:---|---:|---:|---:|")
      for rr in base.robust_rows
        @printf(io, "| %d | %s | %d | %.2f | %.3f |\n", rr.measurement_index, rr.id, rr.stage, rr.t, rr.sigma_factor)
      end
    else
      println(io, "\nRobust R modification (final iteration, solve weights only)")
      println(io, "----------------------------------------------------------------------------------------------")
      @printf(io, "%5s %-28s %6s %8s %13s\n", "Idx", "ID", "Stage", "t", "sigma factor")
      println(io, "----------------------------------------------------------------------------------------------")
      for rr in base.robust_rows
        @printf(io, "%5d %-28s %6d %8.2f %13.3f\n", rr.measurement_index, rr.id, rr.stage, rr.t, rr.sigma_factor)
      end
    end
  end

  # sequential-elimination trace (runse_diagnostics results only)
  if hasproperty(diag, :eliminations) && !isempty(diag.eliminations)
    if format == :markdown
      println(io)
      println(io, "### Sequential elimination (stop: $(diag.stop_reason))")
      println(io, "| # | Idx | ID | Type | Norm.Res before | J before | J after |")
      println(io, "|---:|---:|:---|:---|---:|---:|---:|")
      for t in diag.eliminations
        @printf(io, "| %d | %d | %s | %s | %.5f | %.6f | %.6f |\n", t.elimination, t.measurement_index, t.id, string(t.typ), t.normalized_residual_before, t.objective_before, t.objective_after)
      end
    else
      println(io, "\nSequential elimination (stop: $(diag.stop_reason))")
      println(io, "----------------------------------------------------------------------------------------------")
      @printf(io, "%3s %5s %-28s %-12s %14s %12s %12s\n", "#", "Idx", "ID", "Type", "Norm.Res bef.", "J before", "J after")
      println(io, "----------------------------------------------------------------------------------------------")
      for t in diag.eliminations
        @printf(io, "%3d %5d %-28s %-12s %14.5f %12.6f %12.6f\n", t.elimination, t.measurement_index, t.id, string(t.typ), t.normalized_residual_before, t.objective_before, t.objective_after)
      end
    end
  elseif hasproperty(diag, :stop_reason)
    println(io, format == :markdown ? "\n### Sequential elimination\n- none performed (stop: $(diag.stop_reason))" : "\nSequential elimination: none performed (stop: $(diag.stop_reason))")
  end

  if hasproperty(diag, :rerun) && !isnothing(diag.rerun)
    rerun = diag.rerun
    rerun_summary = summarize_se_diagnostics(rerun.diagnostics)
    if format == :markdown
      println(io)
      println(io, "### Deactivate-and-rerun")
      println(io, "- **Deactivated measurement:** idx=$(rerun.deactivated_measurement_index), id=$(rerun.deactivated_measurement_id)")
      @printf(io, "- **Objective before:** %.6f\n", base.objective.value)
      @printf(io, "- **Objective after:** %.6f\n", rerun.diagnostics.objective.value)
      @printf(io, "- **Objective after stats:** dof=%d, z=%.3f, within_3sigma=%s\n", rerun_summary.objective.dof, rerun_summary.objective.zscore, string(rerun_summary.objective.within_3sigma))
      println(io, "- **Global consistency after rerun:** $(rerun_summary.global_consistency) — $(_global_consistency_label(rerun_summary))")
      println(io, "- **Interpretation after rerun:** $(rerun_summary.reason)")
    else
      println(io, "\nDeactivate-and-rerun")
      println(io, "--------------------------------------------------------------")
      @printf(io, "Deactivated measurement: idx=%d id=%s\n", rerun.deactivated_measurement_index, rerun.deactivated_measurement_id)
      @printf(io, "Objective before: %.6f\n", base.objective.value)
      @printf(io, "Objective after : %.6f\n", rerun.diagnostics.objective.value)
      @printf(io, "Objective after stats: dof=%d, z=%.3f, within_3sigma=%s\n", rerun_summary.objective.dof, rerun_summary.objective.zscore, string(rerun_summary.objective.within_3sigma))
      @printf(io, "Global consistency after rerun: %s (%s)\n", string(rerun_summary.global_consistency), _global_consistency_label(rerun_summary))
      @printf(io, "Interpretation after rerun: %s\n", rerun_summary.reason)
    end
  end
end

# Backward-compatible positional-IO method kept for existing callers/tests.
function print_se_diagnostics(io::IO, diag; topN::Int = 10, format::Symbol = :plain)
  return print_se_diagnostics(diag; io = io, topN = topN, format = format)
end

"""
    se_view(net) -> NamedTuple

Static SE view of the network (SE phase 3): reports how the estimator will
treat controllers, shunts, links, and measurements, WITHOUT mutating `net`.

Guarantee: `runse!` never invokes the outer control loop. Every controller is
frozen at its current operating point; the state estimation is a snapshot.
The concept transformation table applies: OLTC/PST/Schraegregler act as
transformers with fixed taps, SVC/STATCOM and reactors/capacitor banks act
as shunts (their susceptance is estimable per phase 2, case A), and series
compensation acts as a fixed branch impedance.

Returned fields:
- `controllers`: one row per registered outer controller
  `(name, type, enabled)`, all frozen during SE
- `shunts`: one row per shunt released via `setShuntEstimation!`
  `(busIdx, busName, status)` with `status = :released` (has an active direct
  bay measurement), `:frozen_no_measurement`, or `:injection_mode_rejected`
- `links`: `(closed_clusters, open_links)`; each closed cluster is
  `(representative, members)` in original bus indices, and the SE runs on the
  contracted net that fuses each cluster onto its representative
- `excluded_measurements`: the notes of the measurement remapping
  (link-referenced flow measurements, partial-cluster injections, collapsed
  branches, aggregations), reasons as in `runse!`

See `print_se_view` for a formatted report.
"""
function se_view(net::Net)
  # controllers: all frozen during SE, listed for transparency
  controllers = NamedTuple[]
  for c in collect_outer_controllers(net)
    push!(controllers, (name = control_name(c), type = Symbol(typeof(c)), enabled = control_enabled(c)))
  end

  # shunts: classification of the estimation releases
  name_by_idx = _bus_name_by_idx(net)
  busname(b) = get(name_by_idx, b, string(b))
  activeMeas = _active_measurements(Measurement[m for m in net.measurements])
  shunts = NamedTuple[]
  for sh in net.shuntVec
    sh.estimate || continue
    status = if sh.model != :Y
      :injection_mode_rejected
    elseif any(m -> (m.typ == ShuntQMeas && m.busIdx == sh.busIdx) || (m.typ == ImagMeas && m.branchIdx === nothing && m.busIdx == sh.busIdx), activeMeas)
      :released
    else
      :frozen_no_measurement
    end
    push!(shunts, (busIdx = sh.busIdx, busName = busname(sh.busIdx), status = status))
  end

  # links: closed clusters (fused by the SE) and open links (real separations)
  reps = _active_link_representative_map(net)
  members_by_rep = Dict{Int,Vector{Int}}()
  for b in eachindex(reps)
    push!(get!(members_by_rep, reps[b], Int[]), b)
  end
  closed_clusters = NamedTuple[]
  for rep in sort(collect(keys(members_by_rep)))
    members = members_by_rep[rep]
    length(members) > 1 || continue
    push!(closed_clusters, (representative = rep, members = members))
  end
  open_links = [(linkIdx = l.linkIdx, fromBus = l.fromBus, toBus = l.toBus) for l in net.linkVec if l.status != 1]

  # excluded / transformed measurements: exactly what the estimator will do
  prep = _se_prepare(net, Measurement[m for m in net.measurements])

  # released transformer taps: mode and fixed direction alpha per
  # released trafo, for the same at-a-glance transparency as the shunts
  taps = NamedTuple[]
  for (k, br) in enumerate(net.branchVec)
    br.tap_est_mode == :none && continue
    push!(taps, (branch = k, name = getCompName(br.comp), mode = br.tap_est_mode, alpha_deg = br.tap_est_alpha_deg, closed = _branch_terminal_state(br) == :closed))
  end

  return (controllers = controllers, shunts = shunts, taps = taps, links = (closed_clusters = closed_clusters, open_links = open_links), excluded_measurements = prep.notes)
end

"""
    print_se_view(view; io=stdout, format=:plain)

Pretty-print an `se_view` report (`format = :plain` or `:markdown`).
"""
function print_se_view(view; io = stdout, format::Symbol = :plain)
  md = format == :markdown
  println(io, md ? "## SE view (frozen operating point)" : "SE view (frozen operating point)")
  md || println(io, "--------------------------------------------------------------")
  println(io, md ? "- **Controllers (frozen):** $(length(view.controllers))" : "Controllers (frozen): $(length(view.controllers))")
  for c in view.controllers
    println(io, md ? "  - $(c.name) [$(c.type)] enabled=$(c.enabled)" : "  $(c.name) [$(c.type)] enabled=$(c.enabled)")
  end
  println(io, md ? "- **Shunt estimation releases:** $(length(view.shunts))" : "Shunt estimation releases: $(length(view.shunts))")
  for s in view.shunts
    println(io, md ? "  - bus $(s.busName): $(s.status)" : "  bus $(s.busName): $(s.status)")
  end
  ntaps = hasproperty(view, :taps) ? length(view.taps) : 0
  println(io, md ? "- **Tap estimation releases:** $(ntaps)" : "Tap estimation releases: $(ntaps)")
  for t in (ntaps > 0 ? view.taps : NamedTuple[])
    println(io, md ? "  - branch $(t.branch) ($(t.name)): mode $(t.mode)$(t.mode in (:pst, :both) ? ", alpha $(t.alpha_deg) deg" : "")$(t.closed ? "" : " (open, ignored)")" : "  branch $(t.branch) ($(t.name)): mode $(t.mode)$(t.mode in (:pst, :both) ? ", alpha $(t.alpha_deg) deg" : "")$(t.closed ? "" : " (open, ignored)")")
  end
  ncl = length(view.links.closed_clusters)
  nol = length(view.links.open_links)
  println(io, md ? "- **Links:** $(ncl) closed cluster(s) fused, $(nol) open link(s) kept as separations" : "Links: $(ncl) closed cluster(s) fused, $(nol) open link(s) kept as separations")
  for cl in view.links.closed_clusters
    println(io, md ? "  - representative bus $(cl.representative), members $(cl.members)" : "  representative bus $(cl.representative), members $(cl.members)")
  end
  println(io, md ? "- **Excluded/transformed measurements:** $(length(view.excluded_measurements))" : "Excluded/transformed measurements: $(length(view.excluded_measurements))")
  for n in view.excluded_measurements
    println(io, md ? "  - $(n.id): $(n.reason) ($(n.detail))" : "  $(n.id): $(n.reason) ($(n.detail))")
  end
  return nothing
end

function print_se_view(io::IO, view; format::Symbol = :plain)
  return print_se_view(view; io = io, format = format)
end

# ---------------------------------------------------------------------------
# SE chain start (SE phase 5): PF from the estimated state
# ---------------------------------------------------------------------------

## Last SE start state per net (weak keys: entries die with their net). One
## entry = (vm, va per original bus; pinj/qinj in MW/MVar per original bus,
## cluster sums on the representative, members 0; converged flag).
const _SE_START_STATES = WeakKeyDict{Any,Any}()

function _register_se_start!(net, vm::Vector{Float64}, va::Vector{Float64}, pinj::Vector{Float64}, qinj::Vector{Float64}, converged::Bool)
  _SE_START_STATES[net] = (vm = vm, va = va, pinj = pinj, qinj = qinj, converged = converged)
  return nothing
end

@inline _se_start_state(net) = get(_SE_START_STATES, net, nothing)

const _SE_STATE_CSV_VERSION = "sparlectra-se-state v1"

"""
    writeSEStateCSV(net; file) -> NamedTuple

Write the last SE start state of `net` (registered by
`runse!(...; updateNet = true)` or `readSEStateCSV!`) as a CSV artifact:

    # sparlectra-se-state v1
    bus,vm_pu,va_deg,pinj_MW,qinj_MVar

One row per bus (bus by name). This is the persistence half of the SE -> PF
chain: a later power-flow run (possibly in another process) restores the
state with `readSEStateCSV!` and starts via `runpf_from_se!`. Errors when no
SE state is registered for `net`.
"""
function writeSEStateCSV(net::Net; file::AbstractString)
  st = _se_start_state(net)
  st === nothing && error("writeSEStateCSV: no state-estimation result registered for this net; run runse!(...; updateNet = true) first")
  name_by_idx = _bus_name_by_idx(net)
  open(file, "w") do io
    println(io, "# ", _SE_STATE_CSV_VERSION)
    println(io, "bus,vm_pu,va_deg,pinj_MW,qinj_MVar")
    for i in eachindex(net.nodeVec)
      println(io, get(name_by_idx, i, string(i)), ",", repr(st.vm[i]), ",", repr(st.va[i]), ",", repr(st.pinj[i]), ",", repr(st.qinj[i]))
    end
  end
  return (count = length(net.nodeVec),)
end

"""
    readSEStateCSV!(net; file) -> NamedTuple

Read an SE state CSV (see `writeSEStateCSV`) written for the SAME case,
write the estimated voltages into `net`, and register the state so
`runpf_from_se!` can start from it. Buses are matched by name; a bus missing
from the file or an unknown bus name aborts with a line-precise error.
"""
function readSEStateCSV!(net::Net; file::AbstractString)
  lines = readlines(file)
  isempty(lines) && error("$(file): empty file")
  expected = "# " * _SE_STATE_CSV_VERSION
  strip(lines[1]) == expected || error("$(file):1: unknown SE state version ('$(strip(lines[1]))'; expected '$(expected)')")
  n = length(net.nodeVec)
  vm = fill(NaN, n)
  va = fill(NaN, n)
  pinj = zeros(Float64, n)
  qinj = zeros(Float64, n)
  for (ln, raw) in enumerate(lines)
    ln <= 2 && continue   # version + header
    line = strip(raw)
    isempty(line) && continue
    fields = split(line, ",")
    length(fields) == 5 || error("$(file):$(ln): expected 5 fields, got $(length(fields))")
    busIdx = try
      geNetBusIdx(net = net, busName = String(strip(fields[1])))
    catch
      error("$(file):$(ln): unknown bus '$(strip(fields[1]))'")
    end
    vals = map(f -> tryparse(Float64, strip(f)), fields[2:5])
    any(v -> v === nothing, vals) && error("$(file):$(ln): non-numeric value")
    vm[busIdx], va[busIdx], pinj[busIdx], qinj[busIdx] = vals
  end
  missingBuses = [i for i = 1:n if isnan(vm[i])]
  isempty(missingBuses) || error("$(file): no state row for bus(es) $(missingBuses); the file does not match this case")
  for i = 1:n
    setVmVa!(node = net.nodeVec[i], vm_pu = vm[i], va_deg = va[i])
  end
  _register_se_start!(net, vm, va, pinj, qinj, true)
  return (count = n,)
end

"""
    runpf_from_se!(net, maxIte, tolerance=1e-8, verbose=0; mode=:se_state, kwargs...) -> NamedTuple

Run a power flow that starts from the last state-estimation result of `net`
(the SE-to-PF chain). Requires a preceding `runse!(...; updateNet =
true)` on this net or a state restored via `readSEStateCSV!`; a clear error
otherwise.

Modes:
- `:se_state` (default): the PF starts from the ESTIMATED VOLTAGES; the
  model injections stay authoritative, so the measurement/model difference
  goes into the (possibly distributed) slack. This is the semantics of the
  `power_flow.start_mode.profile_source = state_estimation` value.
- `:se_snapshot`: the PF additionally takes the NODAL BALANCES from the
  estimation: the bus load aggregates of the WORKING net are replaced by
  `gen_static - inj_est` for the run and restored afterwards; the persistent
  model is never mutated (same protection philosophy as `updateShunts`).
  With consistent PV setpoints the PF then converges in 0 or 1 iterations
  and the slack pickup stays below tolerance. Config value
  `profile_source = se_snapshot`.

`kwargs` are forwarded to `runpf!` (the flat start is forced off; the
estimated voltages ARE the start). Returns `(iterations, erg, converged,
mode, slack_pickup_mw, slack_pickup_mvar)`; the slack pickup is the PF slack
injection minus the estimated injection at the slack bus.
"""
function runpf_from_se!(net::Net, maxIte::Int, tolerance::Float64 = 1e-8, verbose::Int = 0; mode::Symbol = :se_state, kwargs...)
  mode in (:se_state, :se_snapshot) || error("runpf_from_se!: mode must be :se_state or :se_snapshot, got :$(mode)")
  st = _se_start_state(net)
  st === nothing && error("runpf_from_se!: no preceding state-estimation result for this net; run runse!(...; updateNet = true) or readSEStateCSV! first")
  st.converged || error("runpf_from_se!: the registered state-estimation result did not converge; refusing to start a power flow from it")

  # start voltages := estimated state (re-applied defensively; a PF or other
  # write may have overwritten the node metadata since the SE ran)
  for i in eachindex(net.nodeVec)
    setVmVa!(node = net.nodeVec[i], vm_pu = st.vm[i], va_deg = st.va[i])
  end

  savedLoads = nothing
  nProsBefore = length(net.prosumpsVec)
  if mode == :se_snapshot
    # balance takeover on the WORKING net only: per bus, the effective load
    # becomes gen_static - inj_est. The solver assembles injections from the
    # PROSUMERS (node aggregates are derived bookkeeping; a bare aggregate
    # write never reaches the mismatch, which the original phase-5 takeover
    # missed because its noise-free test had takeover == model), so the
    # difference is carried by one temporary delta load prosumer per bus,
    # appended at the tail and removed in the finally block below.
    savedLoads = [(n._pƩLoad, n._qƩLoad) for n in net.nodeVec]
    name_by_idx = _bus_name_by_idx(net)
    iso = Set(net.isoNodes)
    for (i, n) in enumerate(net.nodeVec)
      i in iso && continue
      # the reference injection is the power flow's free balancing variable;
      # pinning it through a load would double-count the slack power
      getNodeType(n) == Slack && continue
      dp = (something(n._pƩGen, 0.0) - st.pinj[i]) - something(n._pƩLoad, 0.0)
      dq = (something(n._qƩGen, 0.0) - st.qinj[i]) - something(n._qƩLoad, 0.0)
      (abs(dp) < 1e-12 && abs(dq) < 1e-12) && continue
      addProsumer!(net = net, busName = name_by_idx[i], type = "ENERGYCONSUMER", p = dp, q = dq, defer_bus_type_refresh = true)
    end
  end

  ite = 0
  erg = 1
  # island-wise nets (the estimator handles them per island) need the
  # island-wise PF solve too; an explicit caller kwarg wins
  pfkw = Dict{Symbol,Any}(kwargs)
  if !haskey(pfkw, :islands_enabled) && length(detect_ac_islands(net).rows) > 1
    pfkw[:islands_enabled] = true
  end
  try
    ite, erg = runpf!(net, maxIte, tolerance, verbose; opt_flatstart = false, pfkw...)
  finally
    # remove the delta prosumers (tail-appended, so pre-existing prosumer
    # indices stay valid) and restore the aggregates their addition bumped
    if length(net.prosumpsVec) > nProsBefore
      deleteat!(net.prosumpsVec, (nProsBefore+1):length(net.prosumpsVec))
    end
    if savedLoads !== nothing
      for (i, n) in enumerate(net.nodeVec)
        net.nodeVec[i]._pƩLoad = savedLoads[i][1]
        net.nodeVec[i]._qƩLoad = savedLoads[i][2]
      end
    end
  end

  # slack pickup: PF slack injection minus the estimated injection there
  # (links never attach to a slack bus, so the uncontracted view is exact)
  # createYBUS drops isolated nodes, so the injections live in iso-shifted
  # indexing; on a multi-island net the pickup sums over every island
  # reference (each island balances its own measurement/model difference)
  iso_sorted = sort(net.isoNodes)
  Vpf = buildVoltageVector(net)
  Vc = isempty(iso_sorted) ? Vpf : Vpf[[i for i in eachindex(Vpf) if !insorted(i, iso_sorted)]]
  Spf = calc_injections(createYBUS(net = net), Vc) .* net.baseMVA
  pickup_mw = 0.0
  pickup_mvar = 0.0
  for (i, nd) in enumerate(net.nodeVec)
    (getNodeType(nd) == Slack && !insorted(i, iso_sorted)) || continue
    s = i - (searchsortedfirst(iso_sorted, i) - 1)
    pickup_mw += real(Spf[s]) - st.pinj[i]
    pickup_mvar += imag(Spf[s]) - st.qinj[i]
  end

  return (iterations = ite, erg = erg, converged = erg == 0, mode = mode, slack_pickup_mw = pickup_mw, slack_pickup_mvar = pickup_mvar)
end
