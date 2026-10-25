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

`voltages` are the estimated bus voltages in the angle reference of the
model: the slack keeps its stored angle (for example the SV angle of a CGMES
import), although the estimate itself runs with the slack at 0.

`vaRefOffsetDeg` is the estimated common PMU reference-angle offset α in
degrees, the angle of the PMU time base against the model angle reference
(near 0 when the PMU angles are in the model reference). It is `nothing`
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

`linkCirculations` holds the estimated circulation states of coupler rings
whose links carry flow measurements (MW for P, Mvar for Q); each counts as
a state in ν. `nothing` when there are none, and always on a net with more
than one AC island (the islands estimate their circulations, the merged
result does not carry them; `calcLinkFlowsSE!` adds them to the link
flows).

`conditioning` describes the linear system of the last Gauss-Newton step:
`(zero_injection, matrix, condition_number, method)`. `zero_injection` is
the mode of the run (`:on`, `:off`, `:constraint`), `matrix` is `:gain`
(the gain matrix `H' W H` with the zero-injection rows weighted) or `:kkt`
(the scaled KKT matrix of `zero_injection = constraint`). The condition
number is the 1-norm estimate of `condestJacobian` on that sparse matrix
(`method = :hager_1norm`, one extra sparse LU), `Inf` when the matrix is
singular. On a net with several AC islands it is the largest island
value. J and the degrees of freedom are `objectiveJ` and `dof`.

`zeroInjectionMultipliers` lists, under `zero_injection = constraint`, the
Lagrange multipliers of the constraint rows at the last step, one entry
per bus: `(bus, kind, lambda_p, lambda_q)` with the bus name, `kind`
`:zero_injection` for a passive bus or `:link_aggregate` for the fused bus
of a closed-link cluster whose members are all passive (its `LINKAGG`
rows), and the multipliers in
1/MW and 1/Mvar (`NaN` where the bus has no row of that kind). A
multiplier is the pull of the measurements against the balance: the
measurement part of the gradient, `H_M' W r_M = C' lambda`. A large value
points at a measurement error near the bus. `nothing` under `on` and
`off`; on a multi-island net the rows of all islands.

`start` names the start of the result run: `:flat` (`flatstart = true`),
`:stored` (`flatstart = false`, the voltages stored in the net) or
`:stored_after_flat`: the run from the flat start did not converge and
the start ladder (`state_estimation.start_ladder`, default on) solved
once more from the stored voltages; `converged` says whether that second
try converged.
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
  # tap estimation: one row per released regulator (tap changer) with its
  # continuous step estimate s, the fixed physical step, the model step and
  # the deviation, `status` (:estimated, :out_of_range or :not_estimated) with its
  # `reason`, note and notice; `nothing` without releases. After the estimation the regulators are FIXED to their
  # physical step and one more run is solved in which they are no states any
  # more: `objectiveJ`/`dof`/`voltages` above are that FIXED final run,
  # `tapFixation` carries J and dof before versus after the fixation, the
  # Delta-J notice and the flags of a rise above the 95 % quantile of
  # chi2(k) (`j_rise_exceeds_quantile`, then `sequential` when the steps
  # fixed one at a time gave the smaller J, `j_simultaneous`). A regulator
  # only relative to a frozen phase regulator stays a state in the fixed run.
  # With state_estimation.tap_estimation.enabled = false no regulator is a
  # state (reason :disabled). Transformers in tap_estimation.exclude have
  # reason :excluded, machine transformers :machine_transformer (status
  # :calculated after a converged run with update_net, unless
  # include_machine_transformers; :q_not_measured_in_group when Q is not
  # measured on every bus of the machine node, note :p_from_schedule when
  # P is not), also after a fallback. `reactance_source` names where the
  # impedance of the transformer comes from (:table, :analytic, :nominal).
  # A run with tap estimation that does not
  # converge is solved once more without it (tap_estimation.fallback):
  # that run is the result, `tapFixation.fallback_used` is set, every row
  # has reason :fallback_not_converged, and `tapFixation.fallback_reason`
  # with `iterations_estimation` names how the first run failed.
  # `tapFixation.failure` `(reason, iterations)` is that failure of a run
  # with tap states (reason :none when it did not fail).
  tapEstimates::Union{Nothing,Vector{NamedTuple}}
  tapFixation::Union{Nothing,NamedTuple}
  # topology validation: the ADVISORY stage-1 precheck findings
  # (status contradictions, dead closed branches, link voltage mismatches,
  # node-balance violations, meshes without an active power row:
  # kind :mesh_without_p_measurement, the solution may not be unique; the
  # tap rows of a changer in such a mesh carry that in their notice);
  # nothing when the precheck is off or clean.
  # Findings never block or mutate anything.
  topologyFindings::Union{Nothing,Vector{NamedTuple}}
  # J_active: `(j, dof, suppressed)` when replacement-mode suppression
  # removed rows from the STATE (frozen suppression-sigma weights). j/dof
  # exclude those rows, answering "how well does the model fit the data
  # the estimator actually trusted". The honest objectiveJ above keeps
  # them at original sigmas as the alarm signal, and the band test stays
  # on objectiveJ. nothing when no row was suppressed.
  activeObjective::Union{Nothing,NamedTuple}
  # link measurements (0.60.0): the estimated circulation states of the
  # coupler rings with measured links, in MW (P) and Mvar (Q), in the order
  # of the working net's link table; `nothing` when there are none or on a
  # merged multi-island result. The link flows themselves come from
  # `calcLinkFlowsSE!` after a run with update_net.
  linkCirculations::Union{Nothing,Vector{Float64}}
  # zero injection (0.60.0): the condition number of the last step's
  # matrix and, under zero_injection = constraint, the multipliers per bus
  conditioning::Union{Nothing,NamedTuple}
  zeroInjectionMultipliers::Union{Nothing,Vector{NamedTuple}}
  # start of the result run (0.60.0): `:flat`, `:stored` (flatstart =
  # false, the voltages stored in the net) or `:stored_after_flat` (the
  # start ladder: the run from the flat start did not converge, this result
  # starts from the stored voltages)
  start::Symbol
end

## Two-stage robust weight modification (SE phase 4): per
## measurement `t = |residual| / sigma` with the ORIGINAL sigma.
## - stage 0, t <= k1: weight unchanged
## - stage 1 (tangential), k1 < t <= k2: sigma_mod = sigma (2t/k1 - 1)
## - stage 2 (constant), t > k2: sigma_mod = |residual| / k1 (the gradient
##   contribution is effectively suppressed)
## The knees default to the historical 3/6 and are bitwise identical to the
## pre-0.10 fixed formulas at those values. Stages are recomputed every
## iteration, so a measurement can recover. Rows with sigma <= 1e-6 are
## exempt here; the row-level method below also exempts every row that
## bad-data elimination and replacement protect.
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

## Stage of one measurement row: rows protected from bad-data elimination
## (`_elimination_protected`: zero-injection rows, also on a LINKAGG
## aggregate of them, SHDERIV rows, sigma <= 1e-6) keep their weight, the
## same test elimination and replacement use. A zero-injection row with
## a real sigma (`zero_injection = on`) encodes network structure; widening
## it would let the estimate leave the balance it stands for.
@inline _robust_stage(m::Measurement, absr::Float64, k1::Float64, k2::Float64) = _elimination_protected(m) ? (0, m.sigma) : _robust_stage(absr, m.sigma, k1, k2)

## Effective robust mode of a config: the legacy Bool `robust` selects the
## staged modification only while robust_mode itself is :off, so every
## pre-0.10 call site (robust = true) keeps its exact behavior.
@inline _se_effective_robust_mode(cfg::StateEstimationConfig)::Symbol = (cfg.robust_mode === :off && cfg.robust) ? :staged : cfg.robust_mode

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
## rn 58 and had that healthy row eliminated.
const _SE_SUPPRESSION_MIN_WII = 0.3

"""
    numeric_rank(A; tol=nothing) -> Int

Compute numerical matrix rank using singular values.
"""
function numeric_rank(A::AbstractMatrix{<:Real}; tol = nothing)
  s = svdvals(Matrix(A))
  isempty(s) && return 0
  if tol === nothing
    tol = eps(eltype(s)) * max(size(A)...) * maximum(s)
  end
  return count(>(tol), s)
end

## Sparse rank via SPQR: below the dense threshold the
## exact dense SVD path keeps its historical semantics bit for bit; above
## it the rank comes from the sparse QR factorization, whose column-norm
## tolerance is fed from the SAME effective tolerance the callers derive
## (rank_tol_factor * jacEps * sigma_max, or the SVD default formula).
"""
    numeric_rank_pivots(H; tol = nothing) -> Int

Numerical rank of the Jacobian read from the LDLt factorization of the
gain matrix `G = H' H` (unit weights, no column scaling: a column that
carries only forward-difference noise must stay below the tolerance, and a
normalization would lift it to unit norm), instead of a second
decomposition of `H`: the number of pivots above `tol^2`, since the pivots
of `H' H` scale like squared singular values and `tol` is the
singular-value bound of `numeric_rank`. The default `tol` is the same
formula as there. The pivot count is exact for a full-rank Jacobian; on a
deficit it can undercount (measured on sp_case1354 with eleven unmeasured
buses: 2700 pivots against a rank of 2702 from the sparse QR, since the
pivots of an LDLt are not the eigenvalues of `G`), and a factorization
CHOLMOD refuses (a singular `G`) states no size at all. So any deficit is
handed to `numeric_rank`, and the method is reported as `:decomposition`;
`:pivots` is reported only for the full-rank answer. Selected with
`state_estimation.rank_method = pivots`; `decomposition` (SVD below
$(_SE_DENSE_LINALG_MAX_N) states, sparse QR above) is the default, the
exact answer at the same cost on every shipped case.
"""
function numeric_rank_pivots(H::AbstractMatrix{<:Real}; tol = nothing)::Tuple{Int,Symbol}
  m, n = size(H)
  (m == 0 || n == 0) && return 0, :pivots
  Hs = H isa SparseMatrixCSC{Float64} ? H : sparse(Matrix{Float64}(H))
  if tol === nothing
    tol = eps(Float64) * max(m, n) * _sigma_max(Hs)
  end
  G = Hs' * Hs
  F = try
    ldlt(Symmetric(G); check = false)
  catch
    nothing
  end
  if F === nothing || !issuccess(F)
    # singular G: the count of the deficit needs the decomposition
    return numeric_rank(H; tol = tol), :decomposition
  end
  d = Vector{Float64}(diag(sparse(F.LD)))
  all(isfinite, d) || return numeric_rank(H; tol = tol), :decomposition
  nrank = count(>(Float64(tol)^2), d)
  # a deficit is confirmed and sized by the decomposition: the pivot count
  # undercounts there (see the docstring), the full-rank verdict is exact
  nrank < n && return numeric_rank(H; tol = tol), :decomposition
  return nrank, :pivots
end

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

## Internal helper: the model angle of the slack in degrees, the angle
## reference of the net. Inside the estimate the slack sits at 0; the warm
## start subtracts this angle from the stored angles, the write-back and the
## result add it again, and the angle predictions (VaMeas, IaMeas) add it, so
## measurements, results and the net stay in the model reference
## (SE start value). A slack without a stored angle has 0.
@inline _se_slack_ref_deg(net::Net, slackIdx::Int)::Float64 = Float64(something(net.nodeVec[slackIdx]._va_deg, 0.0))

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

## Internal helper: the angle the predictions of VaMeas and IaMeas add to
## the estimated angles (radians): the model slack angle `refRad`
## (`_se_slack_ref_deg`) plus the PMU reference-offset state α when it is
## carried. α is therefore the offset of the PMU time base against the
## model angle reference, near 0 when both coincide; under
## pmu_ref_offset = :off the PMU angles are taken in the model reference.
@inline function _va_offset_from_state(x::Vector{Float64}, withVaOffset::Bool, refRad::Float64)
  return withVaOffset ? x[end] + refRad : refRad
end

## Whether the voltages stored in `net` differ from the flat start, the
## condition of the start ladder (`state_estimation.start_ladder`): a net
## without a stored state (every |V| 1 pu or unset, every angle at the slack
## angle) would start the second try where the first one started.
function _se_has_stored_state(net::Net)::Bool
  slackIdx = _find_slack_idx(net)
  vaRef = _se_slack_ref_deg(net, slackIdx)
  for (i, nd) in enumerate(net.nodeVec)
    i == slackIdx && continue
    something(nd._vm_pu, 1.0) != 1.0 && return true
    something(nd._va_deg, 0.0) != vaRef && return true
  end
  return false
end

## the log line of the start ladder (runse! and the SE service write it)
## (`ite` the iterations of the try from the flat start, `nothing` when the
## caller does not know them)
_se_start_ladder_line(ite::Union{Nothing,Int})::String = string("the WLS from the flat start did not converge", ite === nothing ? "" : string(" (", ite, " iteration(s))"), "; second try from the voltages stored in the case (state_estimation.start_ladder)")

## Internal helper: build initial state vector for WLS iterations.
## flatstart=true uses θ=0, Vm=1 except slack magnitude initialized from net;
## with rotate_behind_pst (state_estimation.start_mode.rotate_behind_pst) a
## bus the slack reaches only across phase shifters starts at their shift
## sum, the angles of the power flow flat start relative to the slack
## (_flat_start_shift_angles in network.jl, one function for both).
## flatstart=false takes the stored angles minus the model slack angle
## (_se_slack_ref_deg), so a converged state is the solution of the
## estimate, which keeps the slack at 0.
## State layout: x = [θ(non-slack); Vm(all); B_1..B_k; s; c; α].
## shuntBInit carries the initial susceptances of the k released shunts
## (always the model values, never flat), tapSInit the tap states s. The
## circulation states c of the link rows (`_se_link_ncirc(net)`, 0.60.0)
## start at 0, the minimum-norm split. withVaOffset=true appends the PMU
## reference-offset state α (radians, initialized to 0) as the LAST entry,
## so _va_offset_from_state stays x[end].
function _initial_state_vector(net::Net, slackIdx::Int; flatstart::Bool = true, withVaOffset::Bool = false, shuntBInit::Vector{Float64} = Float64[], tapSInit::Vector{Float64} = Float64[], rotate_behind_pst::Bool = false)
  nbus = length(net.nodeVec)
  nθ = nbus - 1
  k = length(shuntBInit)
  nt = length(tapSInit)
  x = zeros(Float64, nθ + nbus + k + nt + _se_link_ncirc(net) + (withVaOffset ? 1 : 0))
  # nothing on a net without phase shift, so its flat start stays all zero
  theta = flatstart && rotate_behind_pst ? _flat_start_shift_angles(net, slackIdx) : nothing
  vaSlackDeg = _se_slack_ref_deg(net, slackIdx)

  p = 0
  for i = 1:nbus
    if i != slackIdx
      p += 1
      x[p] = flatstart ? (theta === nothing ? 0.0 : theta[i]) : deg2rad(something(net.nodeVec[i]._va_deg, 0.0) - vaSlackDeg)
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

  # released tap regulator states (steps): after the shunt B states,
  # before the trailing PMU offset
  for j = 1:nt
    x[nθ+nbus+k+j] = tapSInit[j]
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
  return Measurement(typ = m.typ, value = m.value, sigma = m.sigma, active = active, busIdx = m.busIdx, branchIdx = m.branchIdx, direction = m.direction, id = m.id, linkIdx = m.linkIdx, zero_injection = m.zero_injection)
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

## Link measurement rows (0.60.0). A closed link has no Ybus entry, so its
## flow is no branch formula. Per link cluster (the buses a closed-link
## component fuses) the link flows f solve the nodal balances A f = b with
## the incidence matrix A of the cluster's links and, per member k,
## b_k = injection spec of k (generation minus load) - shunt draw of k -
## sum of the branch flows leaving k: the computation the allocation after
## the estimate does (`_link_flow_allocation!`). Then f = pinv(A) b + N c.
## pinv(A) b is a function of the fused state alone (unique on a tree, the
## minimum-norm split on a ring); N spans the ring circulations that the
## measured links of the cluster see, and c are extra states, one per
## column of N and per P/Q (state layout [θ; Vm; B; s; c; α]).
##
## One row model per link measurement row; the measurement carries its
## number in `linkIdx` on the working net. `coef` values are entries of
## pinv(A) (or of N for `circ`).
struct _SELinkRow
  typ::MeasurementType
  sign::Float64                                    # -1 for a row read at the to end
  bus::Int                                         # cluster bus in the working net
  spec::ComplexF64                                 # sum of coef * (Pgen - Pload + j (Qgen - Qload)), MVA
  branches::Vector{Tuple{Int,Int,Float64}}         # (working-net branch, end 1 from / 2 to, coef)
  internal::Vector{Tuple{Branch,Float64,Float64}}  # branch inside the cluster: copy, coef of its from end, of its to end
  shunts::Vector{Tuple{Int,Float64}}               # (working-net shunt position, coef)
  circ::Vector{Tuple{Int,Float64}}                 # (circulation state, coef)
  link::Int                                        # linkIdx in the caller's net
end

## The link rows of one working net plus its circulation states: per state
## its type (P or Q) and its column of N as (caller linkIdx, coef).
struct _SELinkTable
  rows::Vector{_SELinkRow}
  ncirc::Int
  circ_typ::Vector{MeasurementType}
  circ_links::Vector{Vector{Tuple{Int,Float64}}}
end

## number of circulation states of a working net (0 without link rows)
@inline _se_link_ncirc(net::Net)::Int = net._se_link_rows isa _SELinkTable ? net._se_link_rows.ncirc : 0

## the circulation states of `x` (they sit in front of the PMU offset), or
## `nothing` when the working net has none
function _se_link_circ(x::Vector{Float64}, net::Net, withVaOffset::Bool)
  nc = _se_link_ncirc(net)
  nc == 0 && return nothing
  last = length(x) - (withVaOffset ? 1 : 0)
  return x[(last-nc+1):last]
end

## write the circulation values `c` into `x` (no-op without states)
function _se_link_circ_set!(x::Vector{Float64}, net::Net, withVaOffset::Bool, c)
  nc = _se_link_ncirc(net)
  (nc == 0 || c === nothing) && return x
  length(c) == nc || error("SE link rows: $(length(c)) circulation values for $(nc) states (this is a bug)")
  last = length(x) - (withVaOffset ? 1 : 0)
  x[(last-nc+1):last] .= c
  return x
end

## insert the tap states `s` into `x` in front of the circulation states and
## the PMU offset
function _se_insert_tap_states(x::Vector{Float64}, s::Vector{Float64}, net::Net, withVaOffset::Bool)
  tail = _se_link_ncirc(net) + (withVaOffset ? 1 : 0)
  return vcat(x[1:(end-tail)], s, x[(end-tail+1):end])
end

## complex end flows (pu) of a branch whose two ends sit on one fused bus at
## voltage v: what calcNetLosses! gives such a branch on the caller's net
function _se_internal_branch_flows(br::Branch, v::ComplexF64)::Tuple{ComplexF64,ComplexF64}
  st = _branch_terminal_state(br)
  st == :open && return (0.0 + 0.0im, 0.0 + 0.0im)
  st == :open_to && return (abs2(v) * conj(_open_terminal_yin(br)), 0.0 + 0.0im)
  st == :open_from && return (0.0 + 0.0im, abs2(v) * conj(_open_terminal_yin(br)))
  Vb = ComplexF64[v, v]
  return (_closed_branch_flow_pu(Vb, 1, 2, br, 1), _closed_branch_flow_pu(Vb, 2, 1, br, 2))
end

## prediction of one link row in MW or Mvar
function _link_row_prediction(row::_SELinkRow, net::Net, V::Vector{ComplexF64}, shuntB::Union{Nothing,Dict{Int,Float64}}, tapOverlay, circ)::Float64
  s = row.spec
  base = net.baseMVA
  for (k, side, a) in row.branches
    br = _prediction_branch(net, k, tapOverlay)
    sb = side == 1 ? branchFlow_pu(br, br.fromBus, br.toBus, 1, V) : branchFlow_pu(br, br.toBus, br.fromBus, 2, V)
    s -= a * sb * base
  end
  vb = V[row.bus]
  for (br, af, at) in row.internal
    sf, st = _se_internal_branch_flows(br, vb)
    s -= (af * sf + at * st) * base
  end
  if !isempty(row.shunts)
    reg = get(net.shuntDict, row.bus, 0)
    for (p, a) in row.shunts
      sh = net.shuntVec[p]
      y = (p == reg && shuntB !== nothing && haskey(shuntB, row.bus)) ? complex(real(sh.y_pu_shunt), shuntB[row.bus]) : sh.y_pu_shunt
      s -= a * abs2(vb) * conj(y) * base
    end
  end
  h = row.typ == PflowMeas ? real(s) : imag(s)
  if circ !== nothing
    for (j, a) in row.circ
      h += a * circ[j]
    end
  end
  return row.sign * h
end

## Build the link rows of the measurement rows `linkMeas` (their positions
## in `meas`) on the dense working net: `reps` the original-space cluster
## map, `busmap`/`branchmap` original -> working net, `shuntOrig` working
## shunt position -> original. Returns the table and the row number per
## entry of `linkMeas` (0 for a row that cannot be built, with its note).
function _se_link_table(net::Net, meas::Vector{Measurement}, linkMeas::Vector{Int}, reps::Vector{Int}, busmap::Dict{Int,Int}, branchmap::Dict{Int,Int}, shuntOrig::Vector{Int}, notes::Vector{NamedTuple})
  rows = _SELinkRow[]
  rowOf = zeros(Int, length(linkMeas))
  circ_typ = MeasurementType[]
  circ_links = Vector{Tuple{Int,Float64}}[]
  # cluster -> its rows (positions in linkMeas)
  byRep = Dict{Int,Vector{Int}}()
  for (q, k) in enumerate(linkMeas)
    l = net.linkVec[meas[k].linkIdx]
    rep = reps[Int(l.fromBus)]
    if !haskey(busmap, rep)
      push!(notes, (reason = :isolated_bus, id = meas[k].id, detail = "link $(meas[k].linkIdx) sits on an isolated bus"))
      @warn "SE: measurement $(meas[k].id) excluded (its link sits on an isolated bus)"
      continue
    end
    push!(get!(byRep, rep, Int[]), q)
  end
  shuntsAt = Dict{Int,Vector{Int}}()   # original bus -> working shunt positions
  for (p, so) in enumerate(shuntOrig)
    sh = net.shuntVec[so]
    (sh.status == 1 && sh.model in (:Y, :VoltageDependentInjection)) || continue
    push!(get!(shuntsAt, Int(sh.busIdx), Int[]), p)
  end
  # one pass over buses, links and branches for all clusters with rows:
  # members, closed links, and the branches with an end in the cluster (the
  # external ones by working-net number and end, the internal ones,
  # collapsed by the contraction, by their original number)
  membersOf = Dict{Int,Vector{Int}}(rep => Int[] for rep in keys(byRep))
  for b in eachindex(reps)
    haskey(membersOf, reps[b]) && push!(membersOf[reps[b]], b)
  end
  linksOf = Dict{Int,Vector{Int}}(rep => Int[] for rep in keys(byRep))
  for (j, l) in enumerate(net.linkVec)
    l.status == 1 || continue
    rep = reps[Int(l.fromBus)]
    haskey(linksOf, rep) && push!(linksOf[rep], j)
  end
  extOf = Dict{Int,Vector{Tuple{Int,Int,Int}}}(rep => Tuple{Int,Int,Int}[] for rep in keys(byRep))     # (original member, working branch, end)
  innerOf = Dict{Int,Vector{Tuple{Int,Int,Int}}}(rep => Tuple{Int,Int,Int}[] for rep in keys(byRep))   # (original branch, from member, to member)
  for (k, br) in enumerate(net.branchVec)
    f = Int(br.fromBus)
    t = Int(br.toBus)
    rf = reps[f]
    rt = reps[t]
    if rf == rt
      haskey(innerOf, rf) && push!(innerOf[rf], (k, f, t))
      continue
    end
    kb = get(branchmap, k, 0)
    kb == 0 && continue   # open or outside the working net: no flow
    haskey(extOf, rf) && push!(extOf[rf], (f, kb, 1))
    haskey(extOf, rt) && push!(extOf[rt], (t, kb, 2))
  end
  for rep in sort!(collect(keys(byRep)))
    members = membersOf[rep]
    pos = Dict{Int,Int}(b => i for (i, b) in enumerate(members))
    clinks = linksOf[rep]
    lpos = Dict{Int,Int}(j => i for (i, j) in enumerate(clinks))
    A = zeros(Float64, length(members), length(clinks))
    for (i, j) in enumerate(clinks)
      A[pos[Int(net.linkVec[j].fromBus)], i] += 1.0
      A[pos[Int(net.linkVec[j].toBus)], i] -= 1.0
    end
    T = pinv(A)
    N = nullspace(A)
    ext = extOf[rep]
    inner = innerOf[rep]
    # circulation states per type, from the measured links of the cluster
    circOf = Dict{MeasurementType,Tuple{Int,Matrix{Float64}}}()   # typ -> (first state, columns)
    if size(N, 2) > 0
      for typ in (PflowMeas, QflowMeas)
        measured = sort!(unique([lpos[meas[linkMeas[q]].linkIdx] for q in byRep[rep] if meas[linkMeas[q]].typ == typ]))
        isempty(measured) && continue
        # the loop directions the measured links see (absolute rank cut,
        # shared with the split after the estimate, losses.jl)
        Nc, _ = _link_measured_cycles(N, measured)
        r = size(Nc, 2)
        r == 0 && continue
        circOf[typ] = (length(circ_typ) + 1, Nc)
        for c = 1:r
          push!(circ_typ, typ)
          push!(circ_links, [(net.linkVec[clinks[i]].linkIdx, Nc[i, c]) for i in eachindex(clinks) if Nc[i, c] != 0.0])
        end
      end
    end
    for q in byRep[rep]
      m = meas[linkMeas[q]]
      jl = lpos[m.linkIdx]
      coef(b) = (a = T[jl, pos[b]]; abs(a) <= 1e-12 ? 0.0 : a)
      spec = 0.0 + 0.0im
      for b in members
        a = coef(b)
        a == 0.0 && continue
        nd = net.nodeVec[b]
        spec += a * complex(something(nd._pƩGen, 0.0) - something(nd._pƩLoad, 0.0), something(nd._qƩGen, 0.0) - something(nd._qƩLoad, 0.0))
      end
      brs = Tuple{Int,Int,Float64}[(kb, side, coef(b)) for (b, kb, side) in ext if coef(b) != 0.0]
      ints = Tuple{Branch,Float64,Float64}[]
      for (k, f, t) in inner
        af = coef(f)
        at = coef(t)
        (af == 0.0 && at == 0.0) && continue
        push!(ints, (deepcopy(net.branchVec[k]), af, at))
      end
      shs = Tuple{Int,Float64}[]
      for b in members
        a = coef(b)
        a == 0.0 && continue
        for p in get(shuntsAt, b, Int[])
          push!(shs, (p, a))
        end
      end
      circ = Tuple{Int,Float64}[]
      if haskey(circOf, m.typ)
        c0, Nc = circOf[m.typ]
        for c in axes(Nc, 2)
          Nc[jl, c] == 0.0 || push!(circ, (c0 + c - 1, Nc[jl, c]))
        end
      end
      push!(rows, _SELinkRow(m.typ, m.direction === :to ? -1.0 : 1.0, busmap[rep], spec, brs, ints, shs, circ, net.linkVec[m.linkIdx].linkIdx))
      rowOf[q] = length(rows)
    end
  end
  return _SELinkTable(rows, length(circ_typ), circ_typ, circ_links), rowOf
end

## Build the SE working context: contract closed-link clusters, remap the
## active measurements per the phase-3 design rules, and turn the link flow
## measurements on closed links into link rows (`_se_link_table`, stored on
## the working net). A link row on an open link is excluded (an open link
## carries no flow).
function _se_prepare(net::Net, measurements::Vector{Measurement})
  activeMeas, activeIdx = _active_measurements_with_indices(measurements)
  notes = NamedTuple[]
  # collected instead of warned per row: a station with a closed coupler can
  # collapse dozens of measured branches, and dozens of identical lines hide
  # everything else in a run log
  collapsed_ids = String[]
  # a working net built by an earlier _se_prepare (diagnostics reruns): its
  # link rows already point into its table
  nested = net._se_link_rows isa _SELinkTable

  # link rows on open links leave here; on a closed link a row stays (it
  # becomes a WLS row below)
  keep = Int[]
  openIds = String[]
  for k in eachindex(activeMeas)
    m = activeMeas[k]
    if m.linkIdx !== nothing && !nested && !(1 <= m.linkIdx <= length(net.linkVec) && net.linkVec[m.linkIdx].status == 1)
      push!(notes, (reason = :open_link, id = m.id, detail = "link $(m.linkIdx) is open or missing"))
      push!(openIds, m.id)
    else
      push!(keep, k)
    end
  end
  isempty(openIds) || @warn "SE: $(length(openIds)) measurement(s) excluded: their link is open (an open link carries no flow)" rows = join(first(openIds, 5), ", ") * (length(openIds) > 5 ? ", ..." : "")
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
  nested && error("SE link rows on a working net that needs a second contraction (this is a bug)")

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

  # the island builder rebuilds shuntVec and shuntDict in the subnet's
  # numbering; keep the original shunt positions for the updateShunts
  # write-back (the subnet keeps net.shuntVec order, filtered to surviving
  # buses)
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

  # link rows: built after the loop, once per cluster
  linkMeas = Int[]

  for k in eachindex(activeMeas)
    m = activeMeas[k]
    if m.linkIdx !== nothing
      push!(linkMeas, k)
      continue
    end
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
      push!(remapped, Measurement(typ = m.typ, value = m.value, sigma = m.sigma, active = true, busIdx = newBus, branchIdx = m.branchIdx, direction = m.direction, id = m.id, zero_injection = m.zero_injection))
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
      # a cluster whose members are all zero-injection rows is the zero
      # injection of the fused bus: the aggregate carries the flag (id
      # unchanged), so it is protected like its members and
      # zero_injection = constraint solves it as a constraint as well
      allZI = all(k -> activeMeas[k].zero_injection, positions)
      aggId = "LINKAGG_$(typ == PinjMeas ? "Pinj" : "Qinj")_rep_$(rep)"
      push!(remapped, Measurement(typ = typ, value = sum(vals), sigma = sqrt(sum(sigs .^ 2)), active = true, busIdx = busmap[rep], id = aggId, zero_injection = allZI))
      push!(origIdx, 0)
      push!(notes, (reason = :aggregated_cluster_injection, id = aggId, detail = "members: " * join(ids, ", ")))
    else
      for k in positions
        push!(notes, (reason = :partial_cluster_injection, id = activeMeas[k].id, detail = "cluster at $(busname(rep)) not fully measured"))
      end
      @warn "SE: $(length(positions)) $(typ) measurement(s) on the link cluster at $(busname(rep)) excluded (not every cluster member is measured; a partial sum would bias the fused balance)"
    end
  end

  if !isempty(linkMeas)
    table, rowOf = _se_link_table(net, activeMeas, linkMeas, reps, busmap, branchmap, shuntOrig, notes)
    for (q, k) in enumerate(linkMeas)
      rowOf[q] == 0 && continue
      m = activeMeas[k]
      push!(remapped, Measurement(typ = m.typ, value = m.value, sigma = m.sigma, active = true, direction = m.direction, id = m.id, linkIdx = rowOf[q]))
      push!(origIdx, activeIdx[k])
    end
    isempty(table.rows) || (snet._se_link_rows = table)
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
both decisions read one scale.

Returns `nothing` when no verdict is possible: no Jacobian was kept, or the
diagnostics refuse this problem size (they error rather than allocate an
unaffordable dense inverse). The caller then suppresses NOTHING. Falling
back to the raw ratio `|r|/sigma` is deliberately not offered: it is the
measure the normalized residual replaced, and a suppression decided on the
wrong scale is worse than none, because nothing downstream shows which
scale was used.

Returns `rn` together with the residual sensitivity `wii`, because the
caller needs both: `rn` says HOW far the row is off, `wii` says whether
that verdict can be trusted for this row at all.

The `wii` guard of `_residual_diagnostics` carries over unchanged: a row
whose `Omega_ii` sits at the numerical floor is not localizable and gets
`rn = 0`, so it can never be suppressed on the strength of a large raw
residual.
"""
function _suppression_normalized_residuals(H, r::Vector{Float64}, w::Vector{Float64}; minStates::Int = state_estimation_config().takahashi_min_states)
  H === nothing && return nothing
  return try
    d = _residual_diagnostics(H, r, w; minStates = minStates)
    (rn = d.rn, wii = d.wii)
  catch err
    # an abort of the run is no missing covariance
    err isa PowerFlowAborted && rethrow()
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
  # the full residual covariance is m×m dense by definition; a set beyond
  # the limit gets the diagnostics WITHOUT the K columns (the caller sees
  # omega = nothing), named in a warning instead of failing the run: the
  # report is on by default since 0.17.4, and a large network must still run
  if need_full_omega && m > _SE_FULL_OMEGA_MAX_M
    @warn "SE diagnostics: the residual-correlation K report needs the full m×m residual covariance and m=$(m) exceeds $(_SE_FULL_OMEGA_MAX_M); the report is skipped for this run"
    need_full_omega = false
  end
  # above the dense state cap the full Omega is never formed: the dense pinv
  # of G is cubic in n. Omega_ii comes from the Takahashi route (or the
  # blockwise sparse solve below), and a K request is served by
  # _residual_correlation_max_blockwise in the caller, column block by
  # column block. Before this, a K request pinned the dense path and a
  # 1354-bus case (2707 states) failed the whole run with "no affordable
  # route", although no part of the diagnostics needs a dense G inverse.
  need_full_omega && n > _SE_DENSE_LINALG_MAX_N && (need_full_omega = false)
  if !need_full_omega && n >= minStates
    sp, takahashi_reason = _residual_diagnostics_takahashi(H, r, w)
    sp !== nothing && return sp
    @warn "SE diagnostics: Takahashi selected inverse unavailable ($(takahashi_reason)); falling back to the $(n > _SE_DENSE_LINALG_MAX_N ? "blockwise sparse" : "dense") path"
  end
  # the dense pinv is cubic in n: below the cap it is the small-system
  # default and the K-report path; above it the blockwise sparse solve
  # (one factorization of G, block right-hand sides) gives the same
  # diagonal without ever holding a dense n×n or m×m matrix
  n > _SE_DENSE_LINALG_MAX_N && return _residual_diagnostics_blockwise(H, r, w)

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
  # tolerance (verified by test, not assumed); H already
  # arrives as a SparseMatrixCSC and passes through untouched
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

## Blockwise residual covariance for systems above the dense cap.
##
## Ω = W⁻¹ − H G⁻¹ Hᵀ with G = Hᵀ W H. Column j of Ω needs one solve
## G x = h_jᵀ (h_j = row j of H) and one sparse product H x, so a block of
## columns costs one multi-RHS solve on the sparse factorization and never
## more than an n×b and an m×b dense block in memory. G is factorized ONCE
## per call (UMFPACK; G is nonsingular whenever the estimate converged).
const _SE_OMEGA_BLOCK = 256

## Sparse G and its factorization. The one expected failure is a singular
## gain matrix (an unobservable set reached the diagnostics); it is turned
## into an error naming that cause, which the callers already report
## (validate_measurements fails the run with the text, the suppression
## round warns and suppresses nothing).
function _omega_gain_factor(Hs::SparseMatrixCSC{Float64,Int}, w::Vector{Float64})
  G = Hs' * (Diagonal(w) * Hs)
  return try
    lu(G)
  catch err
    err isa SingularException || rethrow()
    error("SE diagnostics: the gain matrix G = H'WH is singular (state $(err.info) not observable); no residual covariance exists")
  end
end

## Ω_jj for every row j: per block, X = G⁻¹ Hᵀ[:, block], then
## Ω_jj = 1/w_j − h_j X[:, j] over the nonzeros of row j only.
function _omega_diag_blockwise(Ht::SparseMatrixCSC{Float64,Int}, w::Vector{Float64}, F)
  m = size(Ht, 2)
  d = zeros(Float64, m)
  for c0 = 1:_SE_OMEGA_BLOCK:m
    # one abort point per block: on a large set the blocks add up to minutes
    sparlectra_check_abort()
    cols = c0:min(c0 + _SE_OMEGA_BLOCK - 1, m)
    X = F \ Matrix{Float64}(Ht[:, cols])
    for (jj, j) in enumerate(cols)
      acc = 0.0
      for a = Ht.colptr[j]:(Ht.colptr[j+1]-1)
        acc += Ht.nzval[a] * X[Ht.rowval[a], jj]
      end
      # same clamp as the dense and Takahashi paths: roundoff can push a
      # critical row (Ω_jj = 0 exactly) a hair below zero
      d[j] = max(1.0 / w[j] - acc, 0.0)
    end
  end
  return d
end

## Diagnostics on the blockwise route: the same (rn, wii, omega,
## omega_path, state_variances) tuple as the other two paths, with
## omega = nothing (Ω_ii only, as on the Takahashi path). diag(G⁻¹) comes
## from identity blocks on the same factorization.
function _residual_diagnostics_blockwise(H::AbstractMatrix{<:Real}, r::Vector{Float64}, w::Vector{Float64})
  m, n = size(H)
  Hs = H isa SparseMatrixCSC{Float64,Int} ? H : sparse(Matrix{Float64}(H))
  F = _omega_gain_factor(Hs, w)
  Ht = SparseMatrixCSC(transpose(Hs))
  d = _omega_diag_blockwise(Ht, w, F)
  rn = zeros(Float64, m)
  wii = zeros(Float64, m)
  for i = 1:m
    wii[i] = d[i] * w[i]
    rn[i] = d[i] <= eps(Float64) ? 0.0 : r[i] / sqrt(d[i])
  end
  sv = zeros(Float64, n)
  for c0 = 1:_SE_OMEGA_BLOCK:n
    sparlectra_check_abort()
    cols = c0:min(c0 + _SE_OMEGA_BLOCK - 1, n)
    E = zeros(Float64, n, length(cols))
    for (jj, j) in enumerate(cols)
      E[j, jj] = 1.0
    end
    X = F \ E
    for (jj, j) in enumerate(cols)
      sv[j] = X[j, jj]
    end
  end
  return (rn = rn, wii = wii, omega = nothing, omega_path = :blockwise, state_variances = sv)
end

## K report without the full Ω: the same per-row maximum as
## _residual_correlation_max, computed column block by column block. Two
## passes because every coefficient needs BOTH diagonal entries: the first
## pass collects diag(Ω), the second forms Ω[:, block] = W⁻¹[:, block] −
## H X and reads the column maximum. Ω is symmetric, so the maximum over
## column j is the maximum over row j. Memory: one m×b block.
function _residual_correlation_max_blockwise(H::AbstractMatrix{<:Real}, w::Vector{Float64})
  m = size(H, 1)
  Hs = H isa SparseMatrixCSC{Float64,Int} ? H : sparse(Matrix{Float64}(H))
  F = _omega_gain_factor(Hs, w)
  Ht = SparseMatrixCSC(transpose(Hs))
  d = _omega_diag_blockwise(Ht, w, F)
  kmax = zeros(Float64, m)
  for c0 = 1:_SE_OMEGA_BLOCK:m
    sparlectra_check_abort()
    cols = c0:min(c0 + _SE_OMEGA_BLOCK - 1, m)
    Y = Hs * (F \ Matrix{Float64}(Ht[:, cols]))   # (H G⁻¹ Hᵀ)[:, cols]
    for (jj, j) in enumerate(cols)
      d[j] <= eps(Float64) && continue
      best = 0.0
      for i = 1:m
        (i == j || d[i] <= eps(Float64)) && continue
        # off the diagonal W⁻¹ is zero, so Ω_ij = −Y_ij
        k = abs(Y[i, jj]) / sqrt(d[i] * d[j])
        k > best && (best = k)
      end
      kmax[j] = best
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
##
## `constraintRows` (local indices) are equality constraints
## (`state_estimation.zero_injection = constraint`): `constraint = true`,
## never `suspicious`; the caller has set their normalized residual to zero.
function _build_measurement_suspicion_report(
  activeMeas::Vector{Measurement},
  activeIdx::Vector{Int},
  r::Vector{Float64},
  rn::Vector{Float64};
  normalizedThreshold::Float64 = 3.0,
  wii::Union{Nothing,Vector{Float64}} = nothing,
  wiiThreshold::Float64 = 0.3,
  kmax::Union{Nothing,Vector{Float64}} = nothing,
  constraintRows::Vector{Int} = Int[],
)
  rows = NamedTuple[]
  isCon = falses(length(activeMeas))
  isCon[constraintRows] .= true
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
        suspicious = !isCon[i] && abs(rn[i]) >= normalizedThreshold,
        constraint = isCon[i],
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
## analytically at the step positions of the current state (the changer set
## of the branch, tap_changer_data.jl); branch rows on such a trafo evaluate
## against its scratch branch.
## Link rows (measurements with `linkIdx` on a working net) evaluate their
## row model (`_link_row_prediction`); `linkCirc` carries the circulation
## states of the current state vector (`_se_link_circ`).
function _predict_measurements(measurements::Vector{Measurement}, net::Net, V::Vector{ComplexF64}, Ybus::AbstractMatrix{ComplexF64}; vaOffsetRad::Float64 = 0.0, shuntB::Union{Nothing,Dict{Int,Float64}} = nothing, tapOverlay = nothing, linkCirc = nothing)
  Sbus_MVA = _predicted_injections_MVA(net, V, Ybus, shuntB, tapOverlay)
  h = Vector{Float64}(undef, length(measurements))
  @inbounds for (i, m) in enumerate(measurements)
    h[i] = _predict_row(m, net, V, Sbus_MVA, vaOffsetRad, shuntB, tapOverlay, linkCirc)
  end
  return h
end

@inline function _predict_row(m::Measurement, net::Net, V::Vector{ComplexF64}, Sbus_MVA::Vector{ComplexF64}, vaOffsetRad::Float64, shuntB, tapOverlay, linkCirc)::Float64
  m.linkIdx === nothing && return _measurement_prediction(m, net, V, Sbus_MVA, vaOffsetRad, shuntB, tapOverlay)
  net._se_link_rows isa _SELinkTable || error("link measurement $(m.id) reached the WLS without its link row (this is a bug)")
  return _link_row_prediction(net._se_link_rows.rows[m.linkIdx], net, V, shuntB, tapOverlay, linkCirc)
end

## Row-restricted prediction for the colored FD assembly: re-evaluates only
## `rows` into `h` (preset by the caller), with the SAME injections and the
## SAME per-row function as `_predict_measurements`, so every evaluated
## entry is bitwise the one the full prediction gives.
function _predict_measurement_rows!(h::Vector{Float64}, rows::Vector{Int}, measurements::Vector{Measurement}, net::Net, V::Vector{ComplexF64}, Ybus::AbstractMatrix{ComplexF64}; vaOffsetRad::Float64 = 0.0, shuntB::Union{Nothing,Dict{Int,Float64}} = nothing, tapOverlay = nothing, linkCirc = nothing)
  Sbus_MVA = _predicted_injections_MVA(net, V, Ybus, shuntB, tapOverlay)
  @inbounds for i in rows
    h[i] = _predict_row(measurements[i], net, V, Sbus_MVA, vaOffsetRad, shuntB, tapOverlay, linkCirc)
  end
  return h
end

function _predicted_injections_MVA(net::Net, V::Vector{ComplexF64}, Ybus::AbstractMatrix{ComplexF64}, shuntB::Union{Nothing,Dict{Int,Float64}}, tapOverlay)
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
  return Sbus_MVA
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
## Column coloring for the finite-difference assembly:
## states that share no measurement row can be perturbed TOGETHER, so the
## number of full prediction evaluations drops from the state count to the
## color count (a few dozen on network structures). Derived from the
## sparsity of an assembled H, deterministic (greedy in column order).
struct _FdColoring
  colors::Vector{Int}              # color of each state column
  groups::Vector{Vector{Int}}      # columns per color, ascending
  rows_per_col::Vector{Vector{Int}}  # the column's measurement-row pattern
  row_covered::Vector{BitVector}   # per color: union of the group's rows
  rows_of_color::Vector{Vector{Int}}  # the same union as an ascending row list
end

function _fd_coloring_from_pattern(H::SparseMatrixCSC{Float64})::_FdColoring
  m, n = size(H)
  rv = rowvals(H)
  rows_per_col = [Int[rv[p] for p in nzrange(H, k)] for k = 1:n]
  return _fd_coloring_from_rows(rows_per_col, m, n)
end

## STRUCTURAL coloring input (second iteration of the design): the
## pattern of an assembled H is a POINT property, and on
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
  m = length(measurements)
  n = x_len
  # circulation columns of the link rows sit in front of the PMU offset
  circ0 = n - _se_link_ncirc(net) - (withVaOffset ? 1 : 0) + 1
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
    for ts in tapMap.trafos
      k = ts.branch
      cols = Int[c for c in ts.cols if c != 0]
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
    if mm.linkIdx !== nothing
      # link row: the cluster bus, both terminals of every external branch
      # of the cluster with their tap columns, a released shunt at the
      # cluster bus, and the circulation columns of the row
      net._se_link_rows isa _SELinkTable || return nothing
      lr = net._se_link_rows.rows[mm.linkIdx]
      c = vac(lr.bus)
      c == 0 || push!(cols, c)
      push!(cols, vmc(lr.bus))
      for (k, _, _) in lr.branches
        br = net.branchVec[k]
        for b in (Int(br.fromBus), Int(br.toBus))
          cb = vac(b)
          cb == 0 || push!(cols, cb)
          push!(cols, vmc(b))
        end
        append!(cols, get(tap_cols_of_branch, k, Int[]))
      end
      haskey(shunt_col, lr.bus) && push!(cols, shunt_col[lr.bus])
      for (j, _) in lr.circ
        push!(cols, circ0 + j - 1)
      end
    elseif t == VmMeas
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
  return _FdColoring(colors, groups, rows_per_col, row_covered, [findall(rc) for rc in row_covered])
end

## Returns:
## - H: Jacobian as a SparseMatrixCSC of size (m_measurements × n_states)
##   (the dense m×n allocation was the large-network
##   blocker, tens of gigabytes at 25k buses for a matrix with a handful
##   of entries per row)
## - h0: base prediction vector h(x), reused by caller to avoid recomputation.
##
## With a `coloring` the assembly evaluates ONE
## perturbed prediction per color instead of one per state; the values are
## then written in the SAME column-ascending order with the SAME per-entry
## formula, so the resulting CSC is bit-identical to the per-column
## assembly as long as the coloring's pattern holds. A difference showing
## up OUTSIDE the pattern union of a color's group means an exact zero of
## the base point came alive; the assembly then falls back to the
## per-column path for this call with one warning, correct over fast.
##
## `check_pattern = false` evaluates per color only the rows inside the
## group's pattern (`rows_of_color`) instead of all m rows, about a quarter
## of the work on network structures; the evaluated entries are bitwise
## those of the full prediction (same injections, same per-row function),
## the others are the structural zeros. It skips the detector, which needs
## every row, so the WLS loop arms the detector on the first assembly of
## every solve and runs the restricted path after it (measured on
## case_ACTIVSg25k, 203916 rows, 36 colors: 0.37 s to 0.10 s per assembly,
## H bitwise unchanged).
##
## Sparsity contract: the forward difference is computed column-wise
## exactly as the dense version did (same perturbation, same per-entry
## formula, bit-identical values). An entry whose difference is an EXACT
## zero means the prediction does not depend on that state; those exact
## zeros are what the Takahashi diagnostics note below already relies on
## (sparse(H) there recovered the true structure without a drop
## tolerance), so dropping them structurally changes no number downstream.
function _measurement_jacobian_fd(measurements::Vector{Measurement}, net::Net, x::Vector{Float64}, slackIdx::Int, nbus::Int, Ybus::AbstractMatrix{ComplexF64}; eps::Float64 = 1e-6, withVaOffset::Bool = false, shuntMap::Union{Nothing,ShuntStateMap} = nothing, tapMap::Union{Nothing,TapStateMap} = nothing, coloring::Union{Nothing,_FdColoring} = nothing, check_pattern::Bool = true)
  _sb(xv) = shuntMap === nothing ? nothing : _shunt_b_override(xv, shuntMap)
  # tap columns: forward difference along the analytic tangent dt/ds at the
  # base point x (_tap_overlay_tangent), with the per-column step
  # eps / sigma, so the tap ratio moves by about eps in pu like the other
  # states; every other column keeps eps and is bitwise as before
  _tov(xv) = tapMap === nothing ? nothing : _tap_overlay_tangent(xv, x, tapMap, net)
  # circulation states of the link rows (0.60.0), nothing without
  _lc(xv) = _se_link_circ(xv, net, withVaOffset)
  # the angle rows predict in the model reference (_va_offset_from_state)
  refRad = deg2rad(_se_slack_ref_deg(net, slackIdx))
  h0 =_predict_measurements(measurements, net, _state_to_voltage(x, slackIdx, nbus), Ybus; vaOffsetRad = _va_offset_from_state(x, withVaOffset, refRad), shuntB = _sb(x), tapOverlay = _tov(x), linkCirc = _lc(x))
  m = length(measurements)
  n = length(x)
  epsv = fill(eps, n)
  if tapMap !== nothing
    for ts in tapMap.trafos, c in eachindex(ts.cols)
      ts.cols[c] == 0 || (epsv[ts.cols[c]] = eps / ts.sigma[c])
    end
  end
  predict_at = xk -> _predict_measurements(measurements, net, _state_to_voltage(xk, slackIdx, nbus), Ybus; vaOffsetRad = _va_offset_from_state(xk, withVaOffset, refRad), shuntB = _sb(xk), tapOverlay = _tov(xk), linkCirc = _lc(xk))

  # Every caller gets the coloring, not just the WLS loop: the diagnostics
  # and the observability entry points build their own Jacobian, and on a
  # 25k-bus network an uncolored assembly is 50000 full prediction sweeps
  # (found in a browser pass, where the estimator itself was already
  # fast). Derived here when the caller passed none.
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
        xk[k] += epsv[k]
      end
      if check_pattern
        hk = predict_at(xk)
        covered = coloring.row_covered[c]
        @inbounds for i = 1:m
          if hk[i] != h0[i] && !covered[i]
            pattern_ok = false
            break
          end
        end
        pattern_ok || break
      else
        # only the rows the group's columns touch are re-evaluated; every
        # other row keeps its base value, so its difference is the exact
        # zero the structural pattern says it is
        hk = copy(h0)
        _predict_measurement_rows!(hk, coloring.rows_of_color[c], measurements, net, _state_to_voltage(xk, slackIdx, nbus), Ybus; vaOffsetRad = _va_offset_from_state(xk, withVaOffset, refRad), shuntB = _sb(xk), tapOverlay = _tov(xk), linkCirc = _lc(xk))
      end
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
          push!(vals, d / epsv[k])
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
    xk[k] += epsv[k]
    hk = predict_at(xk)
    @inbounds for i = 1:m
      d = hk[i] - h0[i]
      d == 0.0 && continue
      push!(rowIdx, i)
      push!(colIdx, k)
      push!(vals, d / epsv[k])
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

## Sparse specialization: O(nnz) over the stored entries
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

Normalize every row of `H` to unit norm, then every column, for the rank
test.

Rank does not change under positive row or column scaling, but its
NUMERICAL determination does. Unscaled, a voltage column (entries around 1)
sits next to an angle column whose entries carry the branch admittances of
a 765 kV line, so the columns span orders of magnitude; `sigma_max` then
belongs to the largest scale and genuine small singular values fall below a
relative cut. That read as a rank deficit on a 25000-bus network whose set
was four times overdetermined; the column normalization fixes it.

Columns alone are not enough when one branch is orders of magnitude
stiffer than the rest, a bus coupler entered as a branch with near-zero
impedance (case141: buses 86 and 87, x = 1e-5 pu, admittance 1e5 against
about 250 for the line that feeds them). The coupler's flow and injection
rows dominate both bus columns, the column normalization shrinks the
feeding line's entries to about 1e-5, and the common mode of the two buses
fell below the cut: rank 279 of 281, "not observable", although the
feeding line's flows measure both buses. Normalizing every ROW to unit norm
first takes the coupler's magnitude out of its rows, so the column
normalization afterwards sees comparable entries: rank 281 there, and
unchanged full rank and critical sets on case1354pegase and
case13659pegase.

Rows are scaled to unit norm, deliberately NOT by sigma: dividing by sigma
lifts zero-injection pseudo-measurements (sigma 1e-6 against 0.01 to 24 for
real ones) by six orders of magnitude, they take over `sigma_max`, and
everything else drops below the cut (measured on the 25000-bus network:
raw deficit 191, with sigma row scaling 12308, with column normalization
alone 0). A unit row norm gives every measurement the same weight in the
rank test whatever its sigma or its electrical magnitude.

A row or column of norm zero carries nothing and stays zero, so an
unmeasured state keeps producing the rank deficit it should.
"""
function _column_normalized(H::AbstractMatrix{<:Real})
  Hs = sparse(H)
  vals = nonzeros(Hs)
  # rows first: one pass over the stored entries for the row norms, a
  # second one to divide (CSC has no cheap row access)
  rows = rowvals(Hs)
  rowsq = zeros(Float64, size(Hs, 1))
  for k in eachindex(vals)
    rowsq[rows[k]] += abs2(vals[k])
  end
  for k in eachindex(vals)
    rn = sqrt(rowsq[rows[k]])
    (isfinite(rn) && rn > 0.0) && (vals[k] /= rn)
  end
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

## Critical rows from the diagonal of the residual covariance (issue #394).
## With unit weights Omega = I - H (H'H)^-1 H' is the projector onto the
## residual space, and Omega_ii = 0 exactly when row i is critical (its
## information exists nowhere else, so the estimate reproduces it and its
## residual is structurally zero); weights do not change WHERE the zeros
## are, only the scale in between, so the classification needs no sigmas.
## One factorization of G = H'H and one selected-inverse pass (Takahashi
## above `takahashi_min_states` states, dense pinv below) answer the
## question for every row at once, where the rank test paid one
## decomposition PER row. `wii` is the dimensionless share of a row's own
## error that reaches its residual (0 critical, 1 fully redundant), so the
## threshold is dimensionless too: `tol_wii` (the FD-aware rank tolerance
## squared and made relative to sigma_max, floored at 1e-8, far above the
## rounding of the selected inverse and far below any real redundancy).
## Returns (critical rows, wii per row, method actually used).
function _criticality_from_omega(H::AbstractMatrix{<:Real}, tol_wii::Float64)
  m, n = size(H)
  d = _residual_diagnostics(H, zeros(Float64, m), ones(Float64, m); need_full_omega = false)
  crit = Int[i for i = 1:m if d.wii[i] <= tol_wii]
  return crit, d.wii, d.omega_path
end

## the dimensionless criticality threshold matching the rank tolerance
## `tol` (an absolute singular-value bound): a row whose removal leaves a
## singular value at tol has leverage of order (tol / sigma_max)^2
function _criticality_wii_tolerance(H::AbstractMatrix{<:Real}, tol)::Float64
  floor_wii = 1.0e-8
  tol === nothing && return floor_wii
  smax = _sigma_max(H)
  (isfinite(smax) && smax > 0.0) || return floor_wii
  return max(floor_wii, (Float64(tol) / smax)^2)
end

function _evaluate_observability_from_jacobian(H::AbstractMatrix{<:Real}, activeOriginalIdx::Vector{Int}; tol = nothing, criticality_method::Symbol = :omega, rank_method::Symbol = :decomposition)
  m, n = size(H)
  ν = m - n
  ρ = n > 0 ? m / n : Inf

  rank_method in (:decomposition, :pivots) || error("observability: rank_method must be :decomposition or :pivots (got $(rank_method))")
  rank_method_used = rank_method
  nrank = if rank_method == :pivots
    r, used = numeric_rank_pivots(H; tol = tol)
    rank_method_used = used
    r
  else
    numeric_rank(H; tol = tol)
  end
  adj, ncols = _adjacency_from_sparsity(H)
  mm = _hopcroft_karp(adj, ncols)

  numObs = nrank == ncols
  strObs = mm == ncols

  criticality_method in (:omega, :rank) || error("observability: criticality_method must be :omega or :rank (got $(criticality_method))")
  criticalNum = Int[]
  criticalStr = Int[]
  criticalitySkipped = false
  structuralSkipped = false
  wii = Float64[]
  method_used = criticality_method
  # The structural test (one bipartite matching per row) keeps the m*n
  # budget: it is the smaller cost but still per row, and under :omega the
  # numerical set already carries every structurally critical row (a
  # structurally critical row is numerically critical as well).
  structural_affordable = m * n <= _SE_ROW_CRITICALITY_MAX_MN
  if m > 0 && criticality_method == :omega
    # issue #394: one factorization for all rows, no budget; a system the
    # estimator itself cannot factorize falls back to the budgeted rank
    # tests with one warning
    tol_wii = _criticality_wii_tolerance(H, tol)
    crit_rows = try
      c, wii_all, _ = _criticality_from_omega(H, tol_wii)
      wii = wii_all
      c
    catch err
      # an abort of the run is no unavailable covariance
      err isa PowerFlowAborted && rethrow()
      @warn "observability: criticality via diag(Omega) unavailable ($(sprint(showerror, err))); falling back to the per-row rank tests"
      method_used = :rank
      nothing
    end
    if crit_rows !== nothing
      criticalNum = Int[activeOriginalIdx[i] for i in crit_rows]
      if structural_affordable
        for i = 1:m
          _structural_row_redundant(H, i) || push!(criticalStr, activeOriginalIdx[i])
        end
      else
        structuralSkipped = true
      end
    end
  end
  if m > 0 && method_used == :rank
    # single-row criticality is one rank test per row, so its cost is
    # m times a full decomposition. Step-0 baseline: m*n = 210k (sp_case188,
    # 562 rows) took 7.6 s and 2.9 GB; m*n = 11M (case1354) did not
    # terminate in reasonable time. The m*n budget below keeps the check
    # within roughly a minute; above it the classification is SKIPPED and
    # says so (criticality_skipped), instead of allocating for hours, and
    # quality then reflects observability and redundancy only.
    # Since issue #394 this is the cross-check
    # path (criticality_method = :rank), no longer the default.
    criticalitySkipped = !structural_affordable
    if !criticalitySkipped
      for i = 1:m
        # one rank test per row: an abort point every 64 rows
        i % 64 == 0 && sparlectra_check_abort()
        _numerical_row_redundant(H, i; tol = tol) || push!(criticalNum, activeOriginalIdx[i])
        _structural_row_redundant(H, i) || push!(criticalStr, activeOriginalIdx[i])
      end
    else
      structuralSkipped = true
      @warn "observability: single-row criticality skipped (m=$(m) times n=$(n) exceeds the $(_SE_ROW_CRITICALITY_MAX_MN) budget); the per-row rank tests do not pay for themselves at this size"
    end
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
    structural_criticality_skipped = structuralSkipped,
    criticality_method = method_used,
    rank_method = rank_method_used,
    criticality_wii = wii,
    active_measurement_indices = activeOriginalIdx,
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
    evaluate_observability_matrix(H; tol=nothing, criticality_method=state_estimation_config().criticality_method, rank_method=state_estimation_config().rank_method) -> NamedTuple

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
function evaluate_observability_matrix(H::AbstractMatrix{<:Real}; tol = nothing, criticality_method::Symbol = state_estimation_config().criticality_method, rank_method::Symbol = state_estimation_config().rank_method)
  m, _ = size(H)
  idx = collect(1:m)
  return _evaluate_observability_from_jacobian(H, idx; tol = tol, criticality_method = criticality_method, rank_method = rank_method)
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

  # no row touches the selected states: that IS the verdict (not
  # observable, rank 0), never an error. The tap and shunt release guards
  # of the estimator rely on it: a thinned set can leave a released tap
  # without any measurement, and the guard then freezes the tap instead of
  # aborting the run (a Web UI run).
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
function evaluate_global_observability(net::Net, measurements::Vector{Measurement}; tol = nothing)
  # settings from the registry (#381); `tol` is the numerical rank
  # tolerance of the SVD, not a configuration key
  base = state_estimation_config()
  flatstart = base.flatstart
  jacEps = base.jac_eps
  pmuRefOffset = base.pmu_ref_offset
  rankTolFactor = base.rank_tol_factor
  base_cfg_criticality_method = base.criticality_method
  # island-wise nets: judge every measured island on its own subnet with its
  # own reference (the estimator solves them the same way) and aggregate;
  # unmeasured islands are excluded from the estimation and reported
  part = _se_island_partition(net, measurements)
  if length(part.rows) > 1
    return _evaluate_global_observability_islands(net, measurements, part; tol = tol)
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
  x = _initial_state_vector(net, slackIdx; flatstart = flatstart, withVaOffset = withVaOffset, rotate_behind_pst = base.start_mode.rotate_behind_pst)
  H, _ = _measurement_jacobian_fd(activeMeas, net, x, slackIdx, nbus, Ybus; eps = jacEps, withVaOffset = withVaOffset)

  # Normalize rows, then columns, before the rank test: one tolerance then means the
  # same thing at 14 buses and at 25000 (see _column_normalized).
  Hs = _column_normalized(H)
  effTol = tol
  if effTol === nothing && !isempty(Hs)
    effTol = rankTolFactor * jacEps * _sigma_max(Hs)
  end
  base = _evaluate_observability_from_jacobian(Hs, activeIdx; tol = effTol, criticality_method = base_cfg_criticality_method, rank_method = base.rank_method)
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

function evaluate_global_observability(net::Net; tol = nothing)
  return evaluate_global_observability(net, Measurement[m for m in net.measurements]; tol = tol)
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
  `net.measurements`, `location` names the bus (injections, voltages),
  the oriented branch (flows) or the link (`"link <k>"`, link flow
  measurements, which are WLS rows of their own).
- `cols`: state-column labels in Jacobian order, `"Va(bus)"` for every
  non-slack bus, then `"Vm(bus)"` for every bus, then `"circP(j)"` and
  `"circQ(j)"` for the circulation states of measured coupler rings, plus
  `"alpha"` when PMU `Va` measurements activate the reference-offset state.

Errors when no active measurement exists. Intended for measurement-matrix
reports and placement studies; see the state-estimation suite summary and
the workshop's observability deep dive.
"""
function measurement_jacobian(net::Net)
  base = state_estimation_config()
  flatstart = base.flatstart
  jacEps = base.jac_eps
  pmuRefOffset = base.pmu_ref_offset
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
  x = _initial_state_vector(net, slackIdx; flatstart = flatstart, withVaOffset = withVaOffset, rotate_behind_pst = base.start_mode.rotate_behind_pst)
  H, _ = _measurement_jacobian_fd(activeMeas, net, x, slackIdx, nbus, Ybus; eps = jacEps, withVaOffset = withVaOffset)

  name_by_idx = _bus_name_by_idx(net)
  busname(i) = get(name_by_idx, i, string(i))
  rows = NamedTuple[]
  for (k, m) in enumerate(activeMeas)
    location = if m.linkIdx !== nothing
      string("link ", net._se_link_rows.rows[m.linkIdx].link)
    elseif m.typ in (PflowMeas, QflowMeas)
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
  if net._se_link_rows isa _SELinkTable
    for (j, t) in enumerate(net._se_link_rows.circ_typ)
      push!(cols, string(t == PflowMeas ? "circP(" : "circQ(", j, ")"))
    end
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
function evaluate_local_observability(net::Net, measurements::Vector{Measurement}, stateCols::Vector{Int}; tol = nothing)
  base = state_estimation_config()
  flatstart = base.flatstart
  jacEps = base.jac_eps
  pmuRefOffset = base.pmu_ref_offset
  rankTolFactor = base.rank_tol_factor
  base_cfg_criticality_method = base.criticality_method
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
  x = _initial_state_vector(net, slackIdx; flatstart = flatstart, withVaOffset = withVaOffset, rotate_behind_pst = base.start_mode.rotate_behind_pst)
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
  base = _evaluate_observability_from_jacobian(Hlocal, localOriginalIdx; tol = effTol, criticality_method = base_cfg_criticality_method, rank_method = base.rank_method)

  return merge(base, (rows = localRows, stateCols = copy(stateCols)))
end

function evaluate_local_observability(net::Net, stateCols::Vector{Int}; tol = nothing)
  return evaluate_local_observability(net, Measurement[m for m in net.measurements], stateCols; tol = tol)
end

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
## The join runs over the link groups (`_active_link_representative_map`),
## so it also passes a bus no branch reaches (a busbar with couplers only,
## typed Isolated, in no island of detect_ac_islands); such a bus goes into
## the island of its group with its rows.
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
  reps = _active_link_representative_map(net)
  groupIsland = Dict{Int,Int}()   # link group representative => an island of it
  for b in eachindex(reps)
    a = get(busIsland0, b, 0)
    a == 0 && continue
    g = get!(groupIsland, reps[b], a)
    ra = findroot(g)
    rb = findroot(a)
    ra != rb && (parent[rb] = ra)
  end
  # buses in no island whose link group has one: island root => buses
  linkFed = Dict{Int,Vector{Int}}()
  for b in eachindex(reps)
    (haskey(busIsland0, b) || !haskey(groupIsland, reps[b])) && continue
    push!(get!(linkFed, findroot(groupIsland[reps[b]]), Int[]), b)
  end
  groups = Dict{Int,Vector{Int}}()
  for i in eachindex(raw)
    push!(get!(groups, findroot(i), Int[]), i)
  end
  rows = NamedTuple[]
  for members in sort!(collect(values(groups)); by = first)
    fedBuses = get(linkFed, findroot(members[1]), Int[])
    if length(members) == 1 && isempty(fedBuses)
      push!(rows, raw[members[1]])
    else
      buses = vcat(reduce(vcat, (raw[m].buses for m in members)), fedBuses)
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
    if m.linkIdx !== nothing && net._se_link_rows isa _SELinkTable
      # a link row on a working net: the island of its cluster bus
      iid = get(busIsland, net._se_link_rows.rows[m.linkIdx].bus, 0)
    elseif m.linkIdx !== nothing && 1 <= m.linkIdx <= length(net.linkVec)
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
        Measurement(typ = m.typ, value = m.value, sigma = m.sigma, active = m.active, direction = m.direction, id = m.id, linkIdx = nl, zero_injection = m.zero_injection)
      end
    elseif m.branchIdx !== nothing
      nb = get(sub.branchmap, m.branchIdx, 0)
      if nb == 0
        @warn "SE: measurement $(m.id) excluded (its branch is not inside the island)"
        nothing
      else
        Measurement(typ = m.typ, value = m.value, sigma = m.sigma, active = m.active, branchIdx = nb, direction = m.direction, id = m.id, zero_injection = m.zero_injection)
      end
    elseif m.busIdx !== nothing
      Measurement(typ = m.typ, value = m.value, sigma = m.sigma, active = m.active, busIdx = sub.busmap[m.busIdx], direction = m.direction, id = m.id, zero_injection = m.zero_injection)
    else
      nothing
    end
    remapped === nothing && continue
    push!(meas, remapped)
    push!(keptIdx, k)
  end
  return meas, keptIdx
end

## The per-island estimations (0.30.2): islands of at least
## `power_flow.islands.parallel_min_buses` buses run on their own tasks when
## two or more of them carry measurements, Julia has more than one thread
## and `runtime.parallel.enabled` is on (the master switch, as for the power
## flow); smaller islands follow on this task, in island order. Every island
## works on its own subnet copy and measurement list; nothing mutable is
## shared, and the caller merges the results in island order.
## `machine_rows` is the caller's update_net: `icfg` forces update_net on
## the subnets (the read-back channel), but whether the machine transformers
## are calculated follows the caller, as on a single-island net.
function _se_solve_islands(part, prepared::Vector{Any}, icfg::StateEstimationConfig; machine_rows::Bool, tap_excl::Dict{String,Symbol})::Vector{Any}
  n = length(prepared)
  results = Vector{Any}(nothing, n)
  active = active_sparlectra_config()
  pr = active.runtime.parallel
  min_buses = active.powerflow.islands_parallel_min_buses
  big = [i for i in 1:n if prepared[i] !== nothing && part.rows[i].n_bus >= min_buses]
  cap = parallel_max_tasks(pr)
  if pr.enabled && Threads.nthreads() > 1 && cap > 1 && length(big) >= 2
    tasks = Vector{Union{Nothing,Task}}(nothing, n)
    sem = Base.Semaphore(max(1, cap))
    for i in sort(big; by = i -> part.rows[i].n_bus, rev = true)
      let p = prepared[i]
        tasks[i] = Threads.@spawn begin
          Base.acquire(sem)
          try
            _runse_with_config!(p.sub.inet, p.meas, icfg; machine_rows = machine_rows, tap_excl = tap_excl)
          finally
            Base.release(sem)
          end
        end
      end
    end
    for i in big
      results[i] = fetch(tasks[i])
    end
    for i in 1:n
      (prepared[i] === nothing || tasks[i] !== nothing) && continue
      results[i] = _runse_with_config!(prepared[i].sub.inet, prepared[i].meas, icfg; machine_rows = machine_rows, tap_excl = tap_excl)
    end
  else
    for i in 1:n
      prepared[i] === nothing && continue
      results[i] = _runse_with_config!(prepared[i].sub.inet, prepared[i].meas, icfg; machine_rows = machine_rows, tap_excl = tap_excl)
    end
  end
  return results
end

## island-wise estimation: per-island solve on the subnet, merged result
function _runse_islands!(net::Net, measurements::Vector{Measurement}, cfg::StateEstimationConfig, part; tap_excl::Dict{String,Symbol} = _tap_exclusion_map(net, cfg.tap_estimation))
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
  # conditioning: the island with the largest condition number; the
  # multipliers of every island
  worstCond = nothing
  ziRows = NamedTuple[]
  anyZI = false
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
  # an island that carried tap states (a tapFixation); islands whose rows
  # are all frozen or disabled have rows but no fixation
  anyTapFix = false
  # the failure of the first island WITH tap states that failed; an island
  # without tap states that does not converge is no failure of the tap
  # estimation, and the fallback would solve it the same way again
  tapFailure = (reason = :none, iterations = 0)
  tapJBefore = 0.0
  tapDofBefore = 0
  tapJAfter = 0.0
  tapDofAfter = 0
  tapAllFixed = true
  tapAnyOffgrid = false
  tapIteEst = 0
  tapIteFix = 0
  tapK = 0
  tapAnyExceeds = false
  tapAnySequential = false
  tapJSimSum = 0.0
  tapNotices = String[]
  αdeg = nothing
  nα = 0
  estimated = 0
  # estimated ring circulations of the islands, parent linkIdx => (P, Q)
  linkCircMerged = Dict{Int,Tuple{Float64,Float64}}()
  # 1. every island prepared serially: its subnet and its remapped
  #    measurements (`nothing` when the island carries no usable one)
  prepared = Vector{Any}(nothing, length(part.rows))
  for (iid, row) in enumerate(part.rows)
    idxs = part.perIsland[iid]
    isempty(idxs) && continue
    sub = _se_island_subnet(net, row)
    meas_i, keptIdx = _se_remap_island_measurements(sub, measurements, idxs)
    isempty(meas_i) || (prepared[iid] = (sub = sub, meas = meas_i, keptIdx = keptIdx))
  end
  # update_net on the subnet is the read-back channel for the merged
  # voltages (write-back into the parent happens below, gated by the
  # caller's update_net); shunt write-back is mirrored explicitly
  icfg = _copy_se_config_with(cfg; update_net = true)
  # 2. the estimations, large islands on their own tasks (0.30.2)
  island_results = _se_solve_islands(part, prepared, icfg; machine_rows = cfg.update_net, tap_excl = tap_excl)
  # 3. merged in island order, exactly as the serial loop did
  for (iid, row) in enumerate(part.rows)
    refBus = row.n_ref > 0 ? 0 : row.chosen_ref_bus
    if prepared[iid] === nothing
      push!(islandsInfo, (island = iid, n_bus = row.n_bus, estimated = false, reason = :no_measurements, converged = false, iterations = 0, objectiveJ = 0.0, dof = 0, band_reason = :none, z_wh = NaN, promoted_ref_bus = refBus, va_ref_offset_deg = nothing))
      continue
    end
    sub = prepared[iid].sub
    keptIdx = prepared[iid].keptIdx
    res = island_results[iid]
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
      # link circulations: subnet linkIdx back to the parent's link
      if ist.link_circulation !== nothing
        parentOf = Dict{Int,Int}(v => k for (k, v) in sub.linkmap)
        for (lk, pq) in ist.link_circulation
          k = get(parentOf, lk, 0)
          k == 0 || (linkCircMerged[net.linkVec[k].linkIdx] = pq)
        end
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
    if res.conditioning !== nothing && (worstCond === nothing || res.conditioning.condition_number > worstCond.condition_number)
      worstCond = res.conditioning
    end
    if res.zeroInjectionMultipliers !== nothing
      anyZI = true
      append!(ziRows, res.zeroInjectionMultipliers)
    end
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
        # every field of the row travels (the report fields of the
        # regulator), only the branch index is remapped
        push!(tapRows, merge(tr, (branch = pb, island = iid)))
      end
      if res.tapFixation !== nothing
        tf = res.tapFixation
        anyTapFix = true
        (tapFailure.reason === :none && tf.failure.reason !== :none) && (tapFailure = tf.failure)
        tapJBefore += tf.j_before
        tapDofBefore += tf.dof_before
        tapJAfter += tf.j_after
        tapDofAfter += tf.dof_after
        tapAllFixed &= tf.fixed
        tapAnyOffgrid |= tf.offgrid_residual
        tapIteEst = max(tapIteEst, tf.iterations_estimation)
        tapIteFix = max(tapIteFix, tf.iterations_fixation)
        tapK += tf.k_removed
        tapAnyExceeds |= tf.j_rise_exceeds_quantile
        tapAnySequential |= tf.sequential
        tapJSimSum += tf.j_simultaneous
        isempty(tf.delta_j_notice) || push!(tapNotices, string("island ", iid, ": ", tf.delta_j_notice))
      end
      if cfg.update_taps
        # the island solve wrote the fixed positions into its subnet branch
        # copies; mirror them back onto the parent net (like update_shunts)
        for b in unique(tr.branch for tr in res.tapEstimates)
          (1 <= b <= length(row.branches)) || continue
          _copy_tap_position!(net.branchVec[row.branches[b]], sub.inet.branchVec[b])
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
    # the island runs compute flows and losses on their subnet copies only,
    # the parent net still holds the flows and shunt powers of its last
    # power flow. Recompute both on the merged voltages (and written-back
    # taps), as the single-island path does: calcLinkFlowsSE! reads them and
    # adds the estimated ring circulations. Unestimated islands keep their
    # start voltages in mergedV, so their flows follow those voltages.
    calcNetLosses!(net, mergedV)
    updateShuntPowers!(net = net)
    # merged chain start registration (runpf_from_se!/writeSEStateCSV read
    # it from the parent net); unestimated islands keep their start voltages
    # and a zero balance
    _register_se_start!(net, seVm, seVa, sePinj, seQinj, conv; link_circulation = any(l -> l.status == 1, net.linkVec) ? linkCircMerged : nothing)
  end
  return SEResult(mergedV, conv, iters, norm(residuals), residuals, Jsum, νsum, _j_within_3sigma_band(Jsum, νsum), nα == 1 ? αdeg : nothing, anyShunt ? shuntRows : nothing, anyRobust ? robustRows : nothing, islandsInfo, anyTap ? tapRows : nothing, anyTapFix ? (fixed = tapAllFixed, j_before = tapJBefore, dof_before = tapDofBefore, j_after = tapJAfter, dof_after = tapDofAfter, iterations_estimation = tapIteEst, iterations_fixation = tapIteFix, offgrid_residual = tapAnyOffgrid, k_removed = tapK, delta_j_notice = join(tapNotices, "; "), j_rise_exceeds_quantile = tapAnyExceeds, sequential = tapAnySequential, j_simultaneous = tapJSimSum, fallback_used = false, fallback_reason = :none, failure = tapFailure) : nothing, nothing, anyActiveObj ? (j = activeJsum, dof = activeνsum, suppressed = suppressedSum) : nothing, nothing, worstCond, anyZI ? ziRows : nothing, cfg.flatstart ? :flat : :stored)
end

## `tap_excl`: the transformers the tap estimation keeps out, resolved on
## the net of the caller (`_tap_exclusion_map`); the diagnostics pass the
## map of the uncontracted net into their solve on the contracted copy.
function _runse_configured!(net::Net, measurements::Vector{Measurement}, cfg::StateEstimationConfig; tap_excl::Dict{String,Symbol} = _tap_exclusion_map(net, cfg.tap_estimation))
  # stage-1 topology precheck (advisory): run BEFORE the partition,
  # warn per finding, attach the findings to the result. The estimation
  # itself always proceeds; a finding is information, not a gate.
  topo = nothing
  if cfg.topology_precheck
    pre = validate_topology(net, measurements)
    if !isempty(pre.findings)
      topo = pre.findings
      # one line per RUN, not per finding: a weak measurement set produces a
      # dozen findings, and the result object carries them all anyway. The
      # meshes without P get their own line, a plausibility finding listed
      # first must not hide them.
      mesh = [f for f in topo if f.kind === :mesh_without_p_measurement]
      other = [f for f in topo if f.kind !== :mesh_without_p_measurement]
      isempty(mesh) || @warn "runse!: mesh without P measurement: solution may not be unique ($(length(mesh)) mesh group(s), the estimation proceeds)" meshes = join((f.location for f in first(mesh, 3)), "; ") * (length(mesh) > 3 ? "; ..." : "")
      if !isempty(other)
        shown = [string(f.kind, " at ", f.location, " (", f.evidence, ", ", f.severity, ")") for f in first(other, 3)]
        @warn "runse!: topology precheck reported $(length(other)) finding(s) (advisory, the estimation proceeds)" findings = join(shown, "; ") * (length(other) > 3 ? "; ..." : "")
      end
    end
  end
  part = _se_island_partition(net, measurements)
  solve(c::StateEstimationConfig) = length(part.rows) <= 1 ? _runse_with_config!(net, measurements, c; tap_excl = tap_excl) : _runse_islands!(net, measurements, c, part; tap_excl = tap_excl)
  # the start voltages of the run: a run that does not converge still writes
  # its last iterate into the net (update_net), the fallback and the second
  # try of the start ladder start from these
  vStart = [(nd._vm_pu, nd._va_deg) for nd in net.nodeVec]
  te = cfg.tap_estimation
  # and the model: on a multi-island net an island that converged writes
  # its steps and shunts back (update_taps, update_shunts) even when another
  # island fails, so run 2 would otherwise see the run-1 steps
  withFallback = te.enabled && te.fallback
  # start ladder (state_estimation.start_ladder): only from the flat start
  # and only when the net stores voltages that differ from it; the second
  # try would otherwise repeat the first
  withLadder = cfg.flatstart && cfg.start_ladder && _se_has_stored_state(net)
  modelStart = (withFallback || withLadder) ? _se_writeback_snapshot(net, cfg) : nothing
  restoreStart! = function ()
    for (i, nd) in enumerate(net.nodeVec)
      nd._vm_pu, nd._va_deg = vStart[i]
    end
    _se_writeback_restore!(net, modelStart)
  end
  # one try from the start of `c`: the run, and its fallback without tap
  # estimation when the run with tap states failed
  attempt = function (c::StateEstimationConfig)
    local r = solve(c)
    if withFallback && _se_tap_run_failed(r)
      local reason, ite1 = r.tapFixation.failure
      @warn _se_tap_fallback_line1(reason, ite1)
      restoreStart!()
      local cfgOff = _copy_se_config_with(c; tap_estimation = _tap_config_off(te))
      r = _se_tap_fallback_result(solve(cfgOff), reason, ite1; net = net, measurements = measurements, excl = tap_excl, update_net = c.update_net)
      local line2 = _se_tap_fallback_line2(r)
      r.converged ? (@info line2) : (@warn line2)
    end
    return r
  end
  res = attempt(cfg)
  if withLadder && !res.converged
    @info string("runse!: ", _se_start_ladder_line(res.iterations))
    restoreStart!()
    res = _se_result_with(attempt(_copy_se_config_with(cfg; flatstart = false)); start = :stored_after_flat)
  end
  topo === nothing && return res
  return _se_result_with(res; topologyFindings = topo, tapEstimates = _se_tap_rows_mesh_notice(res.tapEstimates, topo))
end

# the row reasons of a changer that was no state of the run (not released,
# switched off, kept out by the configuration, or the fallback run): its
# step is the transferred one, the mesh notice does not apply
const _TAP_ROW_NOT_RELEASED = (:no_changer, :disabled, :excluded, :machine_transformer, :q_not_measured_in_group, :fallback_not_converged, :fallback_worse_fit)

## the tap rows with the mesh notice appended to every row of a released
## changer whose branch lies in a mesh of a :mesh_without_p_measurement
## finding: Q, |I| and Vm rows alone can admit a second state there, and a
## released changer in such a mesh can end at a wrong step with J near
## zero. Status, reason and note stay; the notice is the only trace on the
## row.
function _se_tap_rows_mesh_notice(rows, topo::Vector{NamedTuple})
  rows === nothing && return nothing
  meshBr = Set{Int}(k for f in topo if f.kind === :mesh_without_p_measurement for k in f.branches)
  isempty(meshBr) && return rows
  text = _tap_notice_text(:mesh_without_p_measurement)
  return NamedTuple[(!(t.reason in _TAP_ROW_NOT_RELEASED) && t.branch in meshBr) ? merge(t, (notice = join(filter(!isempty, [t.notice, text]), "; "),)) : t for t in rows]
end

## the tap estimation settings `te` with the switch off, everything else
## kept (the run without tap estimation of the fallback)
_tap_config_off(te::TapEstimationConfig)::TapEstimationConfig = TapEstimationConfig(enabled = false, fallback = te.fallback, impedance_correction = te.impedance_correction, exclude = te.exclude, include_machine_transformers = te.include_machine_transformers)

## whether a run with tap states failed (the fallback condition, any
## reason): `tapFixation.failure`, set per island run and on a merged
## multi-island result by the first island with tap states that failed
_se_tap_run_failed(res::SEResult)::Bool = res.tapFixation !== nothing && res.tapFixation.failure.reason !== :none

## how a run with tap states failed and after how many iterations of the
## failing solve: `:non_finite` (J or a voltage not finite),
## `:fixation_not_converged` (the estimate converged, the run with the fixed
## steps did not), `:estimation_not_converged`, or `:none` when it did not
## fail
function _se_tap_failure_of(converged::Bool, j::Float64, V::AbstractVector, fixed::Bool, ite::Int, iteEst::Int, iteFix::Int)
  (!isfinite(j) || !all(isfinite, V)) && return (reason = :non_finite, iterations = ite)
  converged && return (reason = :none, iterations = 0)
  fixed && return (reason = :fixation_not_converged, iterations = iteFix)
  return (reason = :estimation_not_converged, iterations = iteEst)
end

## the two run.log lines of the fallback (the service writes the same text)
_se_tap_fallback_line1(reason::Symbol, ite::Int)::String = string("tap estimation: run 1 not converged (", reason, ", ", ite, " iterations), falling back")
## why a failed run 2 failed: `:non_finite` (J or a voltage not finite) or
## `:not_converged`; run.log and the service status use the same reason
_se_tap_run2_reason(res::SEResult)::Symbol = (!isfinite(res.objectiveJ) || !all(isfinite, res.voltages)) ? :non_finite : :not_converged
function _se_tap_fallback_line2(res::SEResult)::String
  res.converged && return string("run 2 without tap estimation: converged, J = ", round(res.objectiveJ; sigdigits = 6), ", dof = ", res.dof)
  return string("run 2 without tap estimation: not converged (", _se_tap_run2_reason(res), ", ", res.iterations, " iterations)")
end

## The result of a run without tap estimation, marked as the fallback of
## run 1 (`reason`, `ite1`): every tap row gets reason
## :fallback_not_converged instead of :disabled (:fallback_worse_fit when
## the SE service rejected a converged run 1 by its acceptance rule
## `worse_fit`), and `tapFixation` says that the fallback was used, with J
## and dof of this run before and after (no fixation) and the iterations of
## run 1. Shared by `runse!` and the SE service, which repeats the
## diagnostics on the run without tap estimation. With `net` and the
## exclusions `excl` of run 1 (`_tap_exclusion_map` on `net` with the
## settings of run 1), the rows of transformers the configuration keeps
## out get that reason back (`:excluded`, `:machine_transformer`), and with
## `update_net` on a converged run the machine transformers are calculated
## from the state written into `net` and the `measurements`.
function _se_tap_fallback_result(res::SEResult, reason::Symbol, ite1::Int; net::Union{Nothing,Net} = nothing, measurements::Vector{Measurement} = Measurement[], excl::Dict{String,Symbol} = Dict{String,Symbol}(), update_net::Bool = false)::SEResult
  rowReason = reason === :worse_fit ? :fallback_worse_fit : :fallback_not_converged
  exclOf(t) = net === nothing ? rowReason : get(excl, net.branchVec[t.branch].comp.cID, rowReason)
  rows = res.tapEstimates === nothing ? nothing : NamedTuple[t.reason === :disabled ? (r = exclOf(t); merge(t, (reason = r, notice = _tap_assumed_notice(t, _tap_notice_text(r))))) : t for t in res.tapEstimates]
  (rows !== nothing && net !== nothing && update_net && res.converged) && (rows = _tap_target_rows(_tap_machine_rows(rows, net, measurements), net))
  notice = string(reason === :worse_fit ? "tap estimation rejected (" : "tap estimation did not converge (", reason, ", ", ite1, " iterations); the result is the run without tap estimation, every tap keeps its transferred step")
  # failure :none: this run has no tap states, run 1 is named by
  # fallback_reason and iterations_estimation
  tf = (fixed = false, j_before = res.objectiveJ, dof_before = res.dof, j_after = res.objectiveJ, dof_after = res.dof, iterations_estimation = ite1, iterations_fixation = 0, offgrid_residual = false, k_removed = 0, delta_j_notice = notice, j_rise_exceeds_quantile = false, sequential = false, j_simultaneous = NaN, fallback_used = true, fallback_reason = reason, failure = (reason = :none, iterations = 0))
  return _se_result_with(res; tapEstimates = rows, tapFixation = tf)
end

## a copy of `res` with some fields replaced (keywords named like the fields)
function _se_result_with(res::SEResult; kwargs...)::SEResult
  vals = Any[get(kwargs, f, getfield(res, f)) for f in fieldnames(SEResult)]
  return SEResult(vals...)
end

## per-island diagnostics merge (multi-island nets): every measured island is
## validated on its own subnet with its own reference; the merged report
## keeps the single-island shape (summed chi-square objective, rankings
## remapped into the caller's measurement indexing plus an `island` tag) and
## adds an `islands` vector with the per-island reports.
# leverage threshold of the suspicion ranking (formerly the keyword
# wiiThreshold = 0.3 of validate_measurements)
const _SE_WII_THRESHOLD = 0.3

function _validate_measurements_islands(net::Net, measurements::Vector{Measurement}, part, tapExcl::Dict{String,Symbol})
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
    r = _validate_measurements(sub.inet, meas_i, tapExcl)
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
    correlation_enabled = state_estimation_config().report_residual_correlation,
    robust_rows = anyRobust ? robustRows : nothing,
    omega_path = length(unique(omegaPaths)) == 1 ? first(omegaPaths) : :mixed,
    state_variances = stateVars,
    gating_notes = gatingNotes,
    result = nothing,
    islands = islandReports,
  )
end

# _se_gain_matrix(H, w, symmetric) -> SparseMatrixCSC: the WLS gain matrix
# G = H' W H. The product is symmetric only up to rounding: entry (i, j)
# sums H[k,i] * (w_k * H[k,j]) and entry (j, i) sums H[k,j] * (w_k * H[k,i]),
# which round differently, so the bitwise symmetry test of the solve sends
# such a G to LU instead of to CHOLMOD on the kept analysis. On
# case_ACTIVSg25k 346k entries of G differ from G' (max |G - G'| 2048 on
# max |G| 8e19), and the LU solve is about ten times slower per iteration
# than the Cholesky one.
# `symmetric = true` (`state_estimation.symmetric_gain`) keeps the upper
# triangle as computed and mirrors it onto the lower one: a G that already
# came out symmetric is unchanged in value, any other one now takes the
# Cholesky route. Off by default: on the se test cases the estimates moved
# by up to 1.6e-10 pu (a deliberately wrong topology that needs 43
# iterations), above the 1e-10 bound the default has to keep.
function _se_gain_matrix(H::SparseMatrixCSC{Float64}, w::AbstractVector{Float64}, symmetric::Bool)
  G = H' * (Diagonal(w) * H)
  return symmetric ? sparse(Symmetric(G, :U)) : G
end

## Standard deviation of the states `cols` from the covariance of the
## estimate (the step states of the tap
## estimation): the diagonal of (H' W H)^-1 at the Jacobian `H` of the last
## WLS iteration on the rows `rows` with the weights `w`; with
## zero-injection constraint rows (`conIdx`, weight 0 in `w`) the upper-left
## block of the inverse of the KKT matrix [G C'; C 0], the covariance of
## the constrained estimate. Returns column => sigma; every column gets
## `nothing` when the matrix is singular (no covariance exists).
function _se_state_sigmas(H::AbstractMatrix{Float64}, w::AbstractVector{Float64}, conIdx::AbstractVector{Int}, rows::Vector{Int}, cols::Vector{Int})::Dict{Int,Union{Nothing,Float64}}
  out = Dict{Int,Union{Nothing,Float64}}(c => nothing for c in cols)
  isempty(cols) && return out
  Hs = sparse(H[rows, :])
  G = Hs' * (Diagonal(w[rows]) * Hs)
  pos = Dict{Int,Int}(r => k for (k, r) in enumerate(rows))
  con = [pos[i] for i in conIdx if haskey(pos, i)]
  K = isempty(con) ? G : (C = Hs[con, :]; [G C'; C spzeros(length(con), length(con))])
  F = lu(K; check = false)
  issuccess(F) || return out
  for c in cols
    e = zeros(size(K, 1))
    e[c] = 1.0
    v = (F \ e)[c]
    out[c] = (isfinite(v) && v > 0.0) ? sqrt(v) : nothing
  end
  return out
end

# Factorization reuse for the WLS normal equations `G * dx = g`. Each
# Gauss-Newton iteration used to hand G to `solve_linear`, whose sparse
# backslash ran a full symbolic analysis plus numeric factorization every
# time; the pattern of G only changes when the row set does (current gate,
# Ia gate, tap fixation), so the analysis is kept per pattern and later
# iterations refactor numerically.
#
# The ROUTE is the one `G \ g` takes, decided on the same tests in the same
# order, so the arithmetic stays that of the old path:
# - non-square, triangular or diagonal G: the old path itself, untouched;
# - bitwise symmetric G (`ishermitian`, which is what backslash checks; the
#   product H'WH is symmetric only up to rounding, some sets come out exact
#   and some do not): CHOLMOD `cholesky(...; check = false)` as backslash
#   does, with the symbolic analysis kept and `cholesky!` redoing only the
#   numeric part (the same factor, the analysis is structural);
# - any other G: LU through `solve_newton_factorized!`, UMFPACK `lu!` on the
#   kept analysis, or KLU with `state_estimation.linear_solver = klu` (see
#   `_se_linear_context`).
# Every failure hands the system to the old path: a Cholesky that is not
# positive definite goes to `solve_linear` (which retries the same Cholesky
# and continues with LU), and `solve_newton_factorized!` falls back to
# `solve_sparse_system` on any LU failure. The fallback chain for a
# singular or ill-conditioned G (LU, sparse QR, small dense SVD) is
# therefore exactly the one the estimator relied on before.
#
# One instance per `_runse_with_config!` call (one island), never shared
# across islands or tasks: the contexts hold mutable factorizations.
mutable struct _SENormalEquationsSolver
  lu_ctx::AbstractNewtonSolverContext
  chol::Union{Nothing,SparseArrays.CHOLMOD.Factor{Float64,Int64}}
  chol_colptr::Vector{Int64}
  chol_rowval::Vector{Int64}
  chol_analyze_count::Int
  chol_refactor_count::Int
end

_SENormalEquationsSolver(backend::Symbol = :umfpack) = _SENormalEquationsSolver(_se_linear_context(backend), nothing, Int64[], Int64[], 0, 0)

# The LU context for `state_estimation.linear_solver`. UMFPACK is chosen
# explicitly, independent of whether some code in the session loaded KLU for
# power mode: the estimator's arithmetic must not depend on what else was
# loaded. :klu takes
# the KLU extension's constructor when it is registered; without it the
# estimator still runs, on UMFPACK, and says once why.
function _se_linear_context(backend::Symbol)
  backend === :umfpack && return UmfpackReuseNewtonContext()
  backend === :klu || throw(ArgumentError("state_estimation.linear_solver must be one of $(STATE_ESTIMATION_LINEAR_SOLVER_VALUES), got :$(backend)"))
  ctor = _POWER_MODE_LINEAR_CONTEXT[]
  if ctor === nothing
    @warn "state_estimation.linear_solver = klu needs the KLU package extension (SparlectraKLUExt, loaded by `using KLU`), which is not loaded; the estimator uses UMFPACK" maxlog = 1
    return UmfpackReuseNewtonContext()
  end
  return ctor()
end

# _se_solve_normal_equations!(ctx, G, g) -> Vector{Float64}: solve the WLS
# normal equations `G * dx = g` with the factorization reuse of `ctx`. Same
# matrix, same right-hand side and same route as `solve_linear(G, g)`; only
# the repeated symbolic analysis is saved. Returns a vector the caller may
# modify.
function _se_solve_normal_equations!(ctx::_SENormalEquationsSolver, G::SparseMatrixCSC{Float64,Int64}, g::AbstractVector{Float64})
  # the tests and their order are those of SparseArrays' `\` for a sparse
  # matrix; the cases it does not factor by Cholesky or LU keep the old path
  if size(G, 1) != size(G, 2) || istril(G) || istriu(G)
    return solve_linear(G, g; allow_pinv = true, svd_max_n = 20_000)
  end
  if ishermitian(G)
    hG = Hermitian(G)
    if ctx.chol === nothing || ctx.chol_colptr != G.colptr || ctx.chol_rowval != G.rowval
      ctx.chol = cholesky(hG; check = false)
      ctx.chol_colptr = copy(G.colptr)
      ctx.chol_rowval = copy(G.rowval)
      ctx.chol_analyze_count += 1
    else
      cholesky!(ctx.chol, hG; check = false)
      ctx.chol_refactor_count += 1
    end
    issuccess(ctx.chol) && return ctx.chol \ g
    # not positive definite: drop the factor (the next iteration analyzes
    # afresh, as the old path did every time) and take the old path, which
    # repeats this Cholesky, fails the same way and continues with LU
    ctx.chol = nothing
    return solve_linear(G, g; allow_pinv = true, svd_max_n = 20_000)
  end
  # the context's solution buffer is overwritten by the next solve, and the
  # caller scales the step in place: hand out a copy
  return copy(solve_newton_factorized!(ctx.lu_ctx, G, g))
end

# _zero_injection_constraint(m) -> Bool: the row is a zero-injection row by
# origin (the `zero_injection` flag the creating code sets; the LINKAGG
# aggregate of a closed-link cluster whose members all carry it) and an
# injection type: the row set that state_estimation.zero_injection =
# constraint solves as equality constraints. Never decided by the id: a
# measurement table uses the meter mRID as id.
@inline _zero_injection_constraint(m::Measurement)::Bool = m.zero_injection && (m.typ == PinjMeas || m.typ == QinjMeas)

# Indices of the rows solved as equality constraints in `meas` under `mode`
# (`state_estimation.zero_injection`): every ZI row with :constraint, none
# otherwise.
_zero_injection_constraint_rows(meas::Vector{Measurement}, mode::Symbol)::Vector{Int} = mode === :constraint ? Int[i for (i, m) in enumerate(meas) if _zero_injection_constraint(m)] : Int[]

# Weights for the residual covariance with constraint rows: the diagnostics
# keep the gain matrix H'WH and take a constraint row as a row of weight
# 1/ZERO_INJECTION_SIGMA^2 (or its own, if tighter). That is the limit the
# constraint stands for, with an error of the order w_measurement / w_ZI on
# Omega, and it keeps the dense, Takahashi and blockwise paths unchanged.
# The fixed ZERO_INJECTION_SIGMA, not state_estimation.zero_injection_sigma_mva:
# that key is the sigma of the weighted form (on) and a constraint row has
# none.
function _constraint_diagnostic_weights(w::Vector{Float64}, conIdx::Vector{Int})::Vector{Float64}
  isempty(conIdx) && return w
  wd = copy(w)
  for i in conIdx
    wd[i] = max(wd[i], inv(ZERO_INJECTION_SIGMA^2))
  end
  return wd
end

# _se_kkt_matrix(G, C) -> (K, a): the scaled KKT matrix of the constrained
# step below, [G a*C'; a*C 0], and its scale a = max diag(G) / max |C|
# (1 when either is zero). One builder for the solve and for the condition
# number of the run.
function _se_kkt_matrix(G::SparseMatrixCSC{Float64,Int64}, C::SparseMatrixCSC{Float64,Int64})
  p = size(C, 1)
  gmax = maximum(abs, diag(G); init = 0.0)
  cmax = maximum(abs, nonzeros(C); init = 0.0)
  a = (gmax > 0.0 && cmax > 0.0) ? gmax / cmax : 1.0
  return [G a*C'; a*C spzeros(Float64, p, p)], a
end

# _se_solve_constrained!(ctx, G, g, C, rc) -> (dx, lambda, fell): one
# Gauss-Newton step of the equality-constrained WLS (Lagrange, KKT of the
# normal equations)
#
#   [ G    a C' ] [ dx ]   [ g    ]
#   [ a C  0    ] [ mu ] = [ a rc ]
#
# G = H_M' W H_M of the measurement rows (constraint rows enter with weight
# zero), C the Jacobian rows of the constraints, rc = 0 - c(x) their
# residuals. The matrix is symmetric indefinite, so it always takes the LU
# route of the context (UMFPACK with the kept analysis, or KLU). The scale
# a = max diag(G) / max |C| brings the constraint rows to the size of G:
# dx does not depend on a, the condition number does. Measured on
# sp_case14/60/188, SmallGrid and MiniGrid: unscaled (a = 1) the KKT matrix
# was worse conditioned than G with sigma rows (up to 8e17), with this a it
# sits at cond(G_M). A singular KKT
# matrix takes the fallback chain of `solve_newton_factorized!` (sparse
# QR, least squares); the second return value says that this happened, so
# the caller can name it (see _kkt_singular_reason). The multipliers of
# the unscaled system are lambda = a * mu (G dx + C' lambda = g).
function _se_solve_constrained!(ctx::_SENormalEquationsSolver, G::SparseMatrixCSC{Float64,Int64}, g::AbstractVector{Float64}, C::SparseMatrixCSC{Float64,Int64}, rc::AbstractVector{Float64})::Tuple{Vector{Float64},Vector{Float64},Bool}
  n = size(G, 1)
  K, a = _se_kkt_matrix(G, C)
  rhs = vcat(g, a .* rc)
  fallbacks = ctx.lu_ctx.fallback_count
  sol = solve_newton_factorized!(ctx.lu_ctx, K, rhs)
  return sol[1:n], a .* sol[(n + 1):end], ctx.lu_ctx.fallback_count > fallbacks
end

# _se_conditioning(G, C, mode) -> NamedTuple: the condition number of the
# matrix of the last Gauss-Newton step, for SEResult.conditioning. With
# constraint rows (C not nothing) the scaled KKT matrix, else the gain
# matrix G as solved (zero-injection rows weighted, the tap prior in).
# condestJacobian: the Hager/Higham 1-norm estimate on one sparse LU, never
# a dense cond. A singular matrix has no finite condition number: the LU
# throws SingularException (UMFPACK) or ZeroPivotException, and the run
# reports Inf; any other error is a bug and propagates.
function _se_conditioning(G::Union{Nothing,SparseMatrixCSC{Float64,Int64}}, C::Union{Nothing,SparseMatrixCSC{Float64,Int64}}, mode::Symbol)
  G === nothing && return nothing
  matrix = C === nothing ? :gain : :kkt
  M = C === nothing ? G : first(_se_kkt_matrix(G, C))
  kappa = try
    condestJacobian(M)
  catch err
    (err isa SingularException || err isa ZeroPivotException) || rethrow()
    Inf
  end
  return (zero_injection = mode, matrix = matrix, condition_number = kappa, method = :hager_1norm)
end

# _se_pin_steps(G, g, pins) -> (G2, g2): the normal equations with the steps
# of some states held (`pins`, pairs column => step). Row and column of a
# held state are cleared except the diagonal, its right-hand side becomes
# diagonal * step, and its coupling moves to the right-hand side of the
# other states, so the solve gives the held steps exactly and the best
# steps of the others given them. Cleared entries stay stored: the pattern
# is that of G and the solver context keeps its analysis. A zero diagonal
# (a dead column) is replaced by 1.
function _se_pin_steps(G::SparseMatrixCSC{Float64,Int64}, g::AbstractVector{Float64}, pins::Vector{Tuple{Int,Float64}})
  G2 = copy(G)
  g2 = collect(Float64, g)
  held = Set{Int}(first.(pins))
  for (c, d) in pins
    g2 .-= Vector(G[:, c]) .* d
  end
  rows = rowvals(G2)
  vals = nonzeros(G2)
  for j in axes(G2, 2), p in nzrange(G2, j)
    i = rows[p]
    (i != j && (i in held || j in held)) && (vals[p] = 0.0)
  end
  for (c, d) in pins
    dc = G[c, c]
    if dc == 0.0
      dc = 1.0
      G2[c, c] = 1.0
    end
    g2[c] = dc * d
  end
  return G2, g2
end

# _se_pin_constraint(C, rc, pins) -> (C2, rc2): the constraint rows of the
# KKT step with the same held steps: their columns are cleared (stored
# zeros) and their contribution moves to the residual.
function _se_pin_constraint(C::SparseMatrixCSC{Float64,Int64}, rc::AbstractVector{Float64}, pins::Vector{Tuple{Int,Float64}})
  C2 = copy(C)
  rc2 = collect(Float64, rc)
  vals = nonzeros(C2)
  for (c, d) in pins
    rc2 .-= Vector(C[:, c]) .* d
    for p in nzrange(C2, c)
      vals[p] = 0.0
    end
  end
  return C2, rc2
end

# Why the KKT matrix of the constrained step was singular: the constraint
# rows C are linearly dependent (two rows state the same balance), or C has
# full row rank and the measurements plus constraints do not fix every
# state (an unobservable part, for example an island that only zero
# injections touch). Sparse QR rank of C', computed once per run and only
# after a fallback.
function _kkt_singular_reason(C::SparseMatrixCSC{Float64,Int64})::String
  p = size(C, 1)
  r = rank(qr(sparse(C')))
  r < p && return "dependent constraints ($(p - r) of $(p) zero-injection rows are combinations of the others)"
  return "the measurements and the constraints do not fix every state (an unobservable part, such as an island only zero injections touch)"
end

## `machine_rows` (default `cfg.update_net`): calculate the steps of the
## machine transformers after a converged run; the island runs pass the
## caller's update_net, their own cfg always has it on.
function _runse_with_config!(net::Net, measurements::Vector{Measurement}, cfg::StateEstimationConfig; machine_rows::Bool = cfg.update_net, tap_excl::Dict{String,Symbol} = _tap_exclusion_map(net, cfg.tap_estimation))
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
  # the model slack angle: the warm start subtracts it, the angle
  # predictions, the write-back and the result add it (_se_slack_ref_deg)
  vaRefDeg = _se_slack_ref_deg(net, slackIdx)
  refRad = deg2rad(vaRefDeg)

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
    xr = _initial_state_vector(net, slackIdx; flatstart = flatstart, withVaOffset = withVaOffset, shuntBInit = bInit, rotate_behind_pst = cfg.start_mode.rotate_behind_pst)
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
      # released shunt while the run reported "not observable" (as
      # measured). Same class as the island map
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

  # Released tap regulators: one step state per released tap changer,
  # after the shunt B states and before the trailing PMU offset. The stamped
  # admittance terms leave the Ybus once; the predictions re-add them at the
  # step positions of the current state (see _predict_measurements).
  tapImpedance = cfg.tap_estimation.impedance_correction
  tapMap, tapOut, tapNoChanger = _build_tap_state_map(net, (nbus - 1) + nbus + length(released) + 1, tapImpedance)
  # transformers whose released regulators all froze: they leave the map,
  # get their model stamp back and are reported at their start position
  tapFrozenStates = copy(tapOut)
  tapPriorCols = Int[]
  # transformers the configuration keeps out (tap_estimation.exclude, machine
  # transformers without include_machine_transformers): `tap_excl`, resolved
  # by the caller on its net before any contraction (_tap_exclusion_map),
  # matched by component id, and frozen with that reason. With the switch
  # off every changer says :disabled (the fallback run gives them their
  # reason back, _se_tap_fallback_result)
  tapExcl = Dict{String,Symbol}()
  if cfg.tap_estimation.enabled
    tapExcl = tap_excl
    _tap_mark_exclusions!(tapFrozenStates, net, tapExcl)
    tapMap === nothing || _tap_mark_exclusions!(tapMap.trafos, net, tapExcl)
  end
  if !cfg.tap_estimation.enabled
    # tap estimation switched off: no regulator becomes a state, every
    # released changer keeps its transferred step and is reported with
    # reason :disabled. The map was never unstamped, so the transformers
    # keep their model stamp in the Ybus.
    for ts in vcat(tapFrozenStates, tapMap === nothing ? TapBranchState[] : tapMap.trafos), c in eachindex(ts.frozen)
      ts.released[c] && (ts.frozen[c] = :disabled)
    end
    tapMap === nothing || append!(tapFrozenStates, TapBranchState[_tap_with_cols(ts, zeros(Int, length(ts.cols)), ts.spos) for ts in tapMap.trafos])
    tapMap = nothing
  end
  if tapMap !== nothing
    _tap_unstamp!(Ybus, net, tapMap)
    # the excluded transformers leave the map before the guards look at it
    if !isempty(tapExcl)
      tapMap, moved = _tap_repack(tapMap, net, Ybus)
      append!(tapFrozenStates, moved)
    end
  end
  if tapMap !== nothing
    # determinability guard: the structural
    # rules (bridge rule, separability, angle anchors with the symmetric
    # start) only write the cause per regulator; what is not estimated
    # decides the rank test on H at the start point below
    xg = _se_insert_tap_states(x, _tap_s_init(tapMap), net, withVaOffset)
    # the rows of the rank test: current rows only above 3 sigma of their
    # predicted current at the start (_tap_rank_rows)
    rankRows = _tap_rank_rows(net, activeMeas, _state_to_voltage(xg, slackIdx, nbus), shuntMap === nothing ? nothing : _shunt_b_override(xg, shuntMap), _tap_overlay(xg, tapMap, net), cfg.ia_current_floor_A)
    rankMeas = activeMeas[rankRows]
    _tap_guard_causes!(net, tapMap, rankMeas, slackIdx, withVaOffset)
    Htest, _ = _measurement_jacobian_fd(activeMeas, net, xg, slackIdx, nbus, Ybus; eps = jacEps, withVaOffset = withVaOffset, shuntMap = shuntMap, tapMap = tapMap)
    # 1) a tap column zero within the rounding noise of its own forward
    # difference on every row: no measurement sees it
    _tap_column_freezes!(tapMap, Htest, _measurement_vector(activeMeas), xg, nbus, jacEps)
    # 2) rank test on the central difference: the forward difference has a
    # truncation error of about jac_eps times the curvature (1e-6 and more
    # relative to the column), far above _TAP_RANK_TOL, so on it no exact
    # null vector would show (measured: smallest scaled singular
    # value 1.1e-8 forward, 2.4e-16 central, symmetrical PST on a bridge)
    Hback, _ = _measurement_jacobian_fd(activeMeas, net, xg, slackIdx, nbus, Ybus; eps = -jacEps, withVaOffset = withVaOffset, shuntMap = shuntMap, tapMap = tapMap)
    _tap_rank_freezes!(tapMap, (0.5 .* (Htest .+ Hback))[rankRows, :], net, rankMeas, slackIdx, withVaOffset)
    _tap_conditional_notes!(tapMap)
    # 3) hint only: a column far below the strongest voltage column of the
    # unweighted H stays estimated and is named
    for (k, c, factor, b) in _tap_weak_hints!(tapMap, Htest, nbus, jacEps)
      chname = only(ts for ts in tapMap.trafos if ts.branch == k).set.changers[c].name
      @info "runse!: tap regulator $(getCompName(net.branchVec[k].comp))$(isempty(chname) ? "" : "/" * chname): $(_tap_weak_text(factor)) (Vm column of bus $(getCompName(net.nodeVec[b].comp))), kept released"
    end
    tapMap, moved = _tap_repack(tapMap, net, Ybus)
    append!(tapFrozenStates, moved)
  end
  if tapMap !== nothing
    # after every freeze: a phase regulator whose angle is tied only over a
    # frozen phase regulator stays a state and is not fixed (6.5)
    _tap_relative_to_frozen!(net, tapMap, tapFrozenStates, activeMeas, slackIdx, withVaOffset)
    tapPriorCols = _tap_prior_columns!(tapMap, net, activeMeas)
    sinit = _tap_s_init(tapMap)
    x = _se_insert_tap_states(x, sinit, net, withVaOffset)
  end
  # the freezes are reported per row (notice); one line per run here
  nTapFrozen = sum((count(c -> ts.released[c] && !(ts.frozen[c] in (:none, :excluded, :machine_transformer)), eachindex(ts.cols)) for ts in vcat(tapFrozenStates, tapMap === nothing ? TapBranchState[] : tapMap.trafos)); init = 0)
  (nTapFrozen > 0 && cfg.tap_estimation.enabled) && @info "runse!: $(nTapFrozen) tap regulator(s) frozen, see tapEstimates[].notice"
  _tap_ov(xv) = tapMap === nothing ? nothing : _tap_overlay(xv, tapMap, net)

  z = _measurement_vector(activeMeas)
  w = _weight_vector(activeMeas)
  # state_estimation.zero_injection = constraint: the ZI rows are equality
  # constraints (Lagrange, see _se_solve_constrained!). They leave the
  # objective (a constraint has no statistical residual, wJ zeroes them)
  # but stay in the redundancy: dof = m_M - (n - p) = m - n as before.
  conIdx = _zero_injection_constraint_rows(activeMeas, cfg.zero_injection)
  wJ = isempty(conIdx) ? w : (wt = copy(w); wt[conIdx] .= 0.0; wt)

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
  local replacementFrozenW = nothing
  # The Jacobian of the LAST iteration, kept for the suppression rounds:
  # they judge on the normalized residual and therefore need the residual
  # covariance diagonal, which needs H. Re-assembling it would pay a second
  # finite-difference pass (the long pole on large systems, 58 s of a 127 s
  # run on case13659pegase); the last iteration's H sits at the same state
  # up to one step below tol, which the smooth Omega does not notice.
  # The closure reassigns this enclosing local on purpose.
  local lastJacobian = nothing
  # factorization reuse for the normal equations across all WLS runs of
  # this call (suppression rounds and the tap fixation re-solve included);
  # the pattern guard inside re-analyzes when the row set changes G
  local seSolver = _SENormalEquationsSolver(cfg.linear_solver)
  # zero_injection = constraint: the reason of the first singular KKT step
  # (nothing while every step factorized); the closure assigns it on purpose
  local kktSingular = nothing
  local symmetricGain = cfg.symmetric_gain
  # the matrix of the last Gauss-Newton step and, with constraint rows, its
  # constraint Jacobian and multipliers, for SEResult.conditioning and
  # zeroInjectionMultipliers; the closure assigns them on purpose
  local lastGain = nothing
  local lastConJac = nothing
  local lastLambda = nothing

  function _wls_solve!(gateImag::Bool)
    converged = false
    # robust weight freeze: the stage assignment can toggle a row across a
    # stage boundary every iteration and keep the state jittering above
    # tol (a weight limit cycle, seen with many released taps). Once the
    # step is small the stages are frozen, and the now-constant weights
    # let the smooth fixed point converge.
    local robustFrozenW = nothing
    local lastStep = Inf
    # FD column coloring: the coloring comes from
    # the STRUCTURAL pattern (topological couplings), not from an
    # assembled H, because point patterns never settled on large systems
    # (micro-derivatives at machine precision appear and vanish bitwise
    # between iterations; measured on case13659pegase, where both a
    # two-point union and a convergence-gated union kept falling back).
    # Engaged from iteration 1 on large systems only (threshold comment
    # at the constant); the in-assembly detector remains the correctness
    # guard and falls back loudly if a coupling were missing. It checks the
    # first assembly of every solve, the later ones evaluate only the rows
    # inside each color's pattern.
    local fdPatternHeld = false
    local fdColoring = length(x) >= _SE_FD_COLORING_MIN_STATES ? _fd_structural_coloring(activeMeas, length(x), slackIdx, nbus, Ybus, withVaOffset, shuntMap, tapMap, net) : nothing
    # tap states of THIS solve (none in the fixation run): the per-iteration
    # step limit in steps, the step size per column for the stop test, and
    # the columns of the parallel-pair prior
    local tapStepLimits = tapMap === nothing ? Dict{Int,Float64}() : _tap_step_limits(tapMap)
    local tapSide = _tap_side_vertices(tapMap)
    local sideHeld = Dict{Int,Float64}()
    local stepScale = _tap_step_scale(tapMap, length(x))
    local priorCols = (tapMap === nothing || tapMap.ncols == 0) ? Int[] : tapPriorCols
    local priorInit = (tapMap === nothing || tapMap.ncols == 0) ? Float64[] : _tap_s_init(tapMap)
    for ite = 1:maxIte
      # one abort check per iteration: a Web UI run that the user gave up
      # on stops here instead of finishing its 20 iterations (the Jacobian
      # assembly below is the long pole, so this is the useful place)
      sparlectra_check_abort()
      # the first assembly of each solve runs the pattern detector over all
      # rows; when the pattern held there, the later ones take the
      # row-restricted path (see the FD notes). A solve whose first
      # assembly fell back keeps the detector, and with it the per-column
      # fallback, on every iteration, as before.
      H, h, fdColored = _measurement_jacobian_fd(activeMeas, net, x, slackIdx, nbus, Ybus; eps = jacEps, withVaOffset = withVaOffset, shuntMap = shuntMap, tapMap = tapMap, coloring = fdColoring, check_pattern = !fdPatternHeld)
      ite == 1 && (fdPatternHeld = fdColored)
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
          stage, σmod = _robust_stage(activeMeas[i], abs(r[i]), robustK1, robustK2)
          stage == 0 && continue
          wSolve[i] = inv(σmod * σmod)
        end
        # freeze when the state has settled, or unconditionally after 15
        # robust iterations: a weight limit cycle with a larger amplitude
        # would otherwise never converge, and constant weights turn the
        # remaining iterations into a plain, smoothly converging WLS
        (lastStep < 1e-4 || ite >= robustStart + 15) && (robustFrozenW = wSolve)
      end
      # constraint rows never enter G or g as weighted rows (weight zero,
      # after the robust stages, which could otherwise give them one); they
      # enter the KKT block of the solve below
      if !isempty(conIdx)
        wSolve = copy(wSolve)
        wSolve[conIdx] .= 0.0
      end

      # Normal equations for WLS:
      #   Δx = (H'WH)^{-1} H'W r
      # with W = diag(w). Sparse throughout: H is a
      # SparseMatrixCSC, so G = H' * (W * H) stays sparse and the solve
      # below takes the sparse factorization path; the former dense
      # HW-broadcast form would have densified through the adjoint
      # broadcast.
      if gateImag && hasImag && ite < imagGate
        Hs = H[nonImagRows, :]
        ws = wSolve[nonImagRows]
        G = _se_gain_matrix(Hs, ws, symmetricGain)
        g = Hs' * (ws .* r[nonImagRows])
      elseif !isempty(iaRows) && !isempty(_ia_gated!(_state_to_voltage(x, slackIdx, nbus)))
        useRows = [i for i in eachindex(activeMeas) if !(i in iaGatedNow)]
        Hs = H[useRows, :]
        ws = wSolve[useRows]
        G = _se_gain_matrix(Hs, ws, symmetricGain)
        g = Hs' * (ws .* r[useRows])
      else
        G = _se_gain_matrix(H, wSolve, symmetricGain)
        g = H' * (wSolve .* r)
      end

      # svd_max_n raised for the WLS normal equations: released extra
      # states (shunt B, tap steps) can turn out COLLECTIVELY dependent even
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
      # applied ONLY where the degeneracy exists (_tap_prior_columns!: a
      # pair no Pflow, Qflow, Imag or Ia row on a member tells apart; with
      # such a row the member is identifiable and the prior would bias the
      # split away from the measurement optimum, seen as J stuck at 0.6
      # instead of 0): a solitary released tap stays prior-free, so the
      # exact-recovery roundtrips are untouched.
      for c0 in priorCols
        λ = 1e-3 * abs(G[c0, c0])
        G[c0, c0] += λ
        g[c0] += λ * (priorInit[c0 - tapMap.col0 + 1] - x[c0])
      end
      lastGain = G
      if isempty(conIdx)
        Δx = _se_solve_normal_equations!(seSolver, G, g)
      else
        lastConJac = H[conIdx, :]
        Δx, lastLambda, kktFell = _se_solve_constrained!(seSolver, G, g, lastConJac, r[conIdx])
        kktFell && kktSingular === nothing && (kktSingular = _kkt_singular_reason(lastConJac))
      end
      # A state with note :sign_from_start keeps the side of its start: a
      # step that reaches or crosses a turning point of |t| goes half the
      # way to it, or onto it and is held there for the rest of this solve
      # once the remaining magnitude effect is below the FD tolerance, and
      # the other states are solved again with that step
      # (_tap_keep_side_pins!). A side state whose step reaches a turning
      # point in the new solve is added; at most one pass per side column.
      if !isempty(tapSide)
        local sidePins = _tap_keep_side_pins!(sideHeld, tapSide, x, Δx, jacEps)
        local sidePass = 0
        while !isempty(sidePins) && sidePass < length(tapSide)
          sidePass += 1
          local Gp, gp = _se_pin_steps(G, g, sidePins)
          if isempty(conIdx)
            Δx = _se_solve_normal_equations!(seSolver, Gp, gp)
          else
            local Cp, rcp = _se_pin_constraint(H[conIdx, :], r[conIdx], sidePins)
            Δx, lastLambda, kktFell = _se_solve_constrained!(seSolver, Gp, gp, Cp, rcp)
            kktFell && kktSingular === nothing && (kktSingular = _kkt_singular_reason(H[conIdx, :]))
          end
          local morePins = [p for p in _tap_keep_side_pins!(sideHeld, tapSide, x, Δx, jacEps) if !(p[1] in first.(sidePins))]
          isempty(morePins) && break
          append!(sidePins, morePins)
        end
      end
      # Released regulator states get a STEP LIMIT, not a hard band.
      # Measured on a CGMES delivery with 11 transformers off by
      # 4 mechanical steps: clamping the state into its declared band (both
      # per component and as an active-set projection) rescued the one
      # runaway transformer but made five that had converged stick to the
      # lower bound and diverge. The transient excursion is part of how the
      # tap state finds its position; blocking it breaks the path. What has
      # to be bounded is how far one iteration may move (_tap_step_limits),
      # so a weakly determined regulator cannot run off its curves in a
      # single step. The whole step is scaled by one factor, so voltage and
      # tap states keep moving along one consistent direction.
      if !isempty(tapStepLimits)
        α = 1.0
        for (col, lim) in tapStepLimits
          d = abs(Δx[col])
          d <= lim && continue
          f = lim / d
          f < α && (α = f)
        end
        α < 1.0 && (Δx .*= α)
      end
      x .+= Δx

      iteDone = ite
      # a tap step counts with its step size (du, or the step angle), so tol
      # keeps its pu meaning for the step states
      lastStep = maximum(abs(Δx[k]) * stepScale[k] for k in eachindex(Δx); init = 0.0)
      if lastStep < tol
        converged = true
        break
      end
    end
  end

  # --- replacement suppression rounds (robust_mode :replacement) -----------
  # Judged on the CONVERGED residuals, never on flat-start transients: solve
  # plain first, put every real row whose NORMALIZED residual reaches
  # k_suppress on the fixed suppression sigma, re-solve with those weights
  # frozen, and repeat until the suppression set is stable (bounded rounds).
  # Virtual rows (sigma <= 1e-6) are never suppressed; statistics keep the
  # original w.
  #
  # One scale for both decisions. Until 0.10.0 this
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
  # The rounds run as a callable so the re-solve with the taps frozen
  # again (below) gets them too. The closure reassigns replacementFrozenW
  # on purpose; every other name in it is local to it.
  function _replacement_rounds!()
    local Vr, hr, rr, supDiag, rnSup, wiiSup, supSet, worst, wRepl
    (robustMode === :replacement && converged) || return
    local prevSet = Set{Int}()
    for _ = 1:3
      Vr = _state_to_voltage(x, slackIdx, nbus)
      hr = _predict_measurements(activeMeas, net, Vr, Ybus; vaOffsetRad = _va_offset_from_state(x, withVaOffset, refRad), shuntB = shuntMap === nothing ? nothing : _shunt_b_override(x, shuntMap), tapOverlay = _tap_ov(x), linkCirc = _se_link_circ(x, net, withVaOffset))
      rr = z - hr
      _wrap_angle_residuals!(rr, activeMeas)
      supDiag = _suppression_normalized_residuals(lastJacobian, rr, _constraint_diagnostic_weights(w, conIdx); minStates = cfg.takahashi_min_states)
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
    return
  end

  _wls_solve!(true)
  _replacement_rounds!()

  # --- tap fixation run (follows a converged released estimate) ------------
  # After the estimation every estimated regulator is FIXED to its nearest
  # physical step (up or down, no release threshold) and one more run is
  # solved in which it is no state any more (it leaves J and the dof). J
  # before versus after the fixation is reported in `tapFixation`; the
  # SEResult voltages/objectiveJ/dof are those of the FIXED final run.
  # A regulator with note :relative_to_frozen is not fixed and stays a
  # state. When the J rise of the fixation exceeds the 95 % quantile of
  # chi2(k), the steps are fixed one at a time as well (the one closest to
  # its physical step first, the rest estimated again) and the result with
  # the smaller J is kept.
  tapMapEst = tapMap
  tapXEst = Float64[]
  tapPositions = nothing
  tapGatedBefore = Set{Int}()
  tapJBefore = NaN
  tapDofBefore = 0
  tapIteEst = 0
  tapConvEst = false
  tapFixed = false
  tapKRemoved = 0
  tapIteFix = 0
  tapExceeds = false
  tapSequential = false
  tapJSim = NaN
  tapJSeq = NaN
  tapQuantile = 0.0
  # standard deviation of each step state at the continuous estimate
  # (column => sigma, `nothing` without a covariance)
  tapSigma = Dict{Int,Union{Nothing,Float64}}()
  # the text the notice of a changer ends with, by (branch, changer): the J
  # of the mirror check of the recheck at the estimate
  tapExtraNotice = Dict{Tuple{Int,Int},String}()
  # the continuous estimate (relative to neutral) and sigma_s of a weak
  # re-estimate of the recheck, by (branch, changer): its row shows them,
  # its step stays the transferred one
  tapWeakEst = Dict{Tuple{Int,Int},Tuple{Float64,Union{Nothing,Float64}}}()
  # J and the gated Ia rows at the current state (Ia rows gated at this
  # state leave J, the same discipline as the final block below)
  tapJNow = function ()
    local Vc, hC, rC, gatedC, jC
    Vc = _state_to_voltage(x, slackIdx, nbus)
    hC = _predict_measurements(activeMeas, net, Vc, Ybus; vaOffsetRad = _va_offset_from_state(x, withVaOffset, refRad), shuntB = shuntMap === nothing ? nothing : _shunt_b_override(x, shuntMap), tapOverlay = _tap_ov(x), linkCirc = _se_link_circ(x, net, withVaOffset))
    rC = z - hC
    _wrap_angle_residuals!(rC, activeMeas)
    gatedC = isempty(iaRows) ? Set{Int}() : copy(_ia_gated!(Vc))
    jC = isempty(gatedC) ? _wls_objective(rC, wJ) : _wls_objective(rC[[i for i in eachindex(rC) if !(i in gatedC)]], wJ[[i for i in eachindex(wJ) if !(i in gatedC)]])
    return jC, gatedC
  end
  # fix the columns `cols` of the current tap map and re-solve; the prior
  # columns follow the renumbering of the kept states. The closure
  # reassigns the enclosing x, tapMap and tapPriorCols on purpose.
  tapFixAndSolve! = function (cols)
    local colmap
    tapMap, x, colmap = _tap_fix_columns(x, tapMap, Set{Int}(cols))
    tapPriorCols = Int[colmap[c] for c in tapPriorCols if haskey(colmap, c)]
    _wls_solve!(false)
    return converged
  end
  # --- recheck at the estimate --------------------------------------------
  # The rank test above judges at the start point, and its verdict moves
  # with that point: in a test with exact measurements and a start three
  # steps off it released a symmetrical phase changer that is not
  # determinable at the truth (fixed at the mirror step 3 instead of 13,
  # J 2.4e-11) and kept out one that is (its ratio partner compensated at
  # -6 instead of -10, J 14.6). So once, at the converged state with the
  # fixed steps, the noise and rank test run again over every fixed
  # changer and every changer the start test kept out with a rank reason
  # (_tap_recheck_map): a fixed one that adds no rank there is withdrawn
  # (:not_observable_at_estimate, back to its transferred step), a kept
  # out one that adds rank is estimated again together with the other
  # free changers of its transformer and fixed (a weak re-estimate is not
  # taken over and the joint estimate runs again without it, see below).
  # A fixed phase changer on a tree edge without an angle anchor
  # with a mirror step (_tap_mirror_step) is solved once at that step;
  # when J there differs from J at the fixed step by less than
  # _TAP_MIRROR_DJ the measurements cannot tell the two steps apart and
  # it is withdrawn as well (:mirror_step_fits). Whether a phase changer has
  # an angle anchor is judged at the converged state as well (the guard of
  # the start test on the recheck rows: an Ia row gated out at the start
  # anchors once its current flows). The fixed partners of a withdrawn
  # changer (or of one put at its transferred step by the band verdict)
  # are estimated again with it at that step (partnerPass!), so no step
  # of the result comes from a joint estimate with a step that was not
  # taken over. It runs after the fixation, and also
  # on a converged run in which the start test kept every released changer
  # out (no exception for it; the caller below sets up an
  # empty map). Every name below is local to this closure; it reassigns
  # only names the WLS closures already reassign, and tapRecheckEst.
  tapRecheckEst = nothing
  tapRecheck! = function ()
    local mirrorJ, mirrorCheck, withdraw!, regName, noCause, unanch, jA, dofA, iteA
    local rmap, xr, rcand, movedIdx, nfin, withdrawList, newly, keptCols, newlyT, again
    local Yr, Hf, Hb, rankRowsR, live, order, candCols, keptR
    local saved, addIdx, addTs, jOf, againF, finTrafos, estTrafos, reasons, colR, valsR, againCols, newTrafos, okR, gatedA, sigA, xestR, viewTrafos, w2
    local touchedB, weakNow, oorR, detLines
    local stale, fixedEst, markStale!, partnerPass!, pairsP, outP
    # J of the fixed run with changer c of transformer j at step m,
    # solved warm from the current state; state, map and the numerics
    # of the last solve are put back (Inf when it does not converge)
    mirrorJ = function (j, c, m)
      local sv, sp, jm
      sv = (x = copy(x), map = tapMap, prior = copy(tapPriorCols), ite = iteDone, conv = converged, gain = lastGain, conJac = lastConJac, lambda = lastLambda, jac = lastJacobian)
      sp = copy(tapMap.trafos[j].spos)
      sp[c] = m
      tapMap = TapStateMap(TapBranchState[k == j ? _tap_with_cols(t, t.cols, sp) : t for (k, t) in enumerate(tapMap.trafos)], tapMap.col0, tapMap.ncols, tapMap.impedance)
      _wls_solve!(false)
      jm = converged ? first(tapJNow()) : Inf
      x = sv.x
      tapMap = sv.map
      tapPriorCols = sv.prior
      iteDone = sv.ite
      converged = sv.conv
      lastGain = sv.gain
      lastConJac = sv.conJac
      lastLambda = sv.lambda
      lastJacobian = sv.jac
      return jm
    end
    # the mirror check over the fixed changers `pairs` (j, c) of the
    # current map: (j, c, J text) of every one to withdraw
    mirrorCheck = function (pairs)
      local out, jf, ts, ch, m, jm
      out = Tuple{Int,Int,String}[]
      isempty(pairs) && return out
      jf = first(tapJNow())
      for (j, c) in pairs
        ts = tapMap.trafos[j]
        (ts.branch, c) in unanch || continue
        ch = ts.set.changers[c]
        m = _tap_mirror_step(ch, ts.spos[c])
        m === nothing && continue
        jm = mirrorJ(j, c, m)
        abs(jm - jf) < _TAP_MIRROR_DJ && push!(out, (j, c, string("J ", @sprintf("%.4g", jf), " at the fixed step ", _tap_notice_step(ch.neutral + ts.spos[c]), ", J ", @sprintf("%.4g", jm), " at the mirror step ", _tap_notice_step(ch.neutral + m))))
      end
      return out
    end
    regName = (ts, c) -> string(getCompName(net.branchVec[ts.branch].comp), isempty(ts.set.changers[c].name) ? "" : "/" * ts.set.changers[c].name)
    # withdraw (j, c, reason) of the current map: frozen with the reason
    # (the freeze vectors are shared with the estimation map, so the
    # report sees it), back to the transferred step, one solve; returns
    # its iterations
    withdraw! = function (list)
      local sp, ts, v
      isempty(list) && return 0
      sp = Dict{Int,Vector{Float64}}()
      for (j, c, reason) in list
        ts = tapMap.trafos[j]
        _tap_freeze!(ts, c, reason)
        v = get!(() -> copy(ts.spos), sp, j)
        v[c] = ts.s0[c]
        @info "runse!: recheck at the estimate: tap regulator $(regName(ts, c)) withdrawn ($(reason)), the transferred step $(_tap_notice_step(ts.set.changers[c].neutral + ts.s0[c])) stays"
      end
      tapMap = TapStateMap(TapBranchState[haskey(sp, k) ? _tap_with_cols(t, t.cols, sp[k]) : t for (k, t) in enumerate(tapMap.trafos)], tapMap.col0, tapMap.ncols, tapMap.impedance)
      _wls_solve!(false)
      return iteDone
    end
    # The partners of a withdrawn changer. A changer withdrawn here
    # (:not_observable_at_estimate, :mirror_step_fits) or put at its
    # transferred step by the band verdict (:out_of_range) was a state of
    # the joint estimate that fixed the other changers of its transformer;
    # their steps belong to a state that no longer exists. Measured on a
    # ratio and phase pair: the phase changer jointly at -11.84 (truth
    # -13) and withdrawn to its transferred -13, the ratio changer kept the
    # -2 of the joint estimate (truth -5). Such a partner is in `stale`
    # ((j, c) of the current map) and is estimated again with the
    # withdrawn changer at its transferred step (partnerPass! below).
    # fixedEst: changer c of transformer j is fixed from an estimate of
    # this run (a state of the estimation map, now a fixed step)
    stale = Set{Tuple{Int,Int}}()
    fixedEst = function (j, c)
      local t
      t = tapMap.trafos[j]
      return t.released[c] && t.frozen[c] === :none && t.cols[c] == 0 && j <= length(tapMapEst.trafos) && tapMapEst.trafos[j].cols[c] != 0
    end
    markStale! = function (wd)
      for (j, c) in wd, k in eachindex(tapMap.trafos[j].cols)
        (k != c && !((j, k) in wd) && fixedEst(j, k)) && push!(stale, (j, k))
      end
      return nothing
    end
    # One pass over the partners `pairs` ((j, c) of the current map):
    # estimated again together, every other changer at its current step
    # (a withdrawn one at its transferred step), then the band verdict, the
    # fixation and the mirror check, as for the first estimate. The rank
    # test at the converged state has kept each of them already (a
    # withdrawn changer comes after its partners in that test), and the
    # sigma_s bound of a newly determinable changer does not apply: this
    # estimate replaces the first one. When the joint run or its fixation
    # does not converge, the state before it is put back and the partners
    # go to their transferred step (:not_converged_at_estimate). Returns
    # the changers this pass withdrew or put at their transferred step.
    partnerPass! = function (pairs)
      local sv, colP, valsP, colsP, trafosP, okP, gatedP, sigP, oorP, viewX, viewSig, wP, out, jOfB
      sv = (x = copy(x), map = tapMap, prior = copy(tapPriorCols), ite = iteDone, gain = lastGain, conJac = lastConJac, lambda = lastLambda, jac = lastJacobian)
      for (j, c) in pairs
        @info "runse!: recheck at the estimate: tap regulator $(regName(tapMap.trafos[j], c)) estimated again, a changer of its transformer was withdrawn to its transferred step"
      end
      colP = tapMap.col0
      valsP = Float64[]
      colsP = Int[]
      trafosP = TapBranchState[]
      for (j, ts) in enumerate(tapMap.trafos)
        let cols = zeros(Int, length(ts.cols))
          for c in eachindex(cols)
            if ts.cols[c] != 0
              push!(valsP, x[ts.cols[c]])
            elseif (j, c) in pairs
              push!(valsP, ts.spos[c])
              push!(colsP, colP)
            else
              continue
            end
            cols[c] = colP
            colP += 1
          end
          push!(trafosP, _tap_with_cols(ts, cols, copy(ts.spos)))
        end
      end
      x = vcat(x[1:(tapMap.col0 - 1)], valsP, x[(tapMap.col0 + tapMap.ncols):end])
      tapMap = TapStateMap(trafosP, tapMap.col0, colP - tapMap.col0, tapMap.impedance)
      tapPriorCols = _tap_prior_columns!(tapMap, net, activeMeas)
      _wls_solve!(false)
      tapIteFix += iteDone
      okP = converged
      oorP = Tuple{TapBranchState,Int,Symbol,Bool}[]
      if okP
        _, gatedP = tapJNow()
        sigP = lastJacobian === nothing ? Dict{Int,Union{Nothing,Float64}}() : _se_state_sigmas(lastJacobian, wJ, conIdx, [i for i in eachindex(activeMeas) if !(i in gatedP)], colsP)
        oorP = _tap_mark_out_of_range!(tapMap, x, sigP, Set{Int}(colsP))
        # the report view: the new estimate and its sigma_s replace the
        # joint ones (the column of the estimation map stays)
        viewX = copy(tapXEst)
        viewSig = copy(tapSigma)
        for (j, c) in pairs
          let col = tapMap.trafos[j].cols[c], ecol = tapMapEst.trafos[j].cols[c]
            viewX[ecol] = x[col]
            viewSig[ecol] = get(sigP, col, nothing)
          end
        end
        okP = tapFixAndSolve!(colsP)
        tapIteFix += iteDone
        if okP
          tapXEst = viewX
          tapSigma = viewSig
        end
      end
      out = Tuple{Int,Int}[]
      if okP
        jOfB = Dict{Int,Int}(ts.branch => j for (j, ts) in enumerate(tapMap.trafos))
        for (ts, c, _, _) in oorP
          push!(out, (jOfB[ts.branch], c))
        end
        wP = Tuple{Int,Int,Symbol}[]
        for (j, c, txt) in mirrorCheck([(j, c) for (j, c) in pairs if tapMap.trafos[j].frozen[c] === :none])
          push!(wP, (j, c, :mirror_step_fits))
          tapExtraNotice[(tapMap.trafos[j].branch, c)] = txt
          push!(out, (j, c))
        end
        tapIteFix += withdraw!(wP)
      else
        for (ts, c, nt, tp) in oorP
          ts.frozen[c] = :none
          ts.notes[c] = nt
          ts.tree_phase[c] = tp
        end
        x = sv.x
        tapMap = sv.map
        tapPriorCols = sv.prior
        iteDone = sv.ite
        converged = true
        lastGain = sv.gain
        lastConJac = sv.conJac
        lastLambda = sv.lambda
        lastJacobian = sv.jac
        tapIteFix += withdraw!([(j, c, :not_converged_at_estimate) for (j, c) in pairs])
        append!(out, pairs)
      end
      return out
    end

    rmap, xr, rcand, movedIdx = _tap_recheck_map(tapMap, tapMapEst, x, tapFrozenStates)
    nfin = length(tapMap.trafos)
    withdrawList = Tuple{Int,Int,Symbol}[]
    newly = Tuple{Int,Int}[]
    keptCols = Set{Int}()
    # (branch, changer) of the phase changers on a tree edge without an
    # angle anchor at the converged state, for the mirror check
    unanch = Set{Tuple{Int,Int}}()
    if !isempty(rcand)
      Yr = copy(Ybus)
      isempty(movedIdx) || _tap_unstamp!(Yr, net, TapStateMap(rmap.trafos[(nfin + 1):end], rmap.col0, 0, rmap.impedance))
      Hf, _ = _measurement_jacobian_fd(activeMeas, net, xr, slackIdx, nbus, Yr; eps = jacEps, withVaOffset = withVaOffset, shuntMap = shuntMap, tapMap = rmap)
      Hb, _ = _measurement_jacobian_fd(activeMeas, net, xr, slackIdx, nbus, Yr; eps = -jacEps, withVaOffset = withVaOffset, shuntMap = shuntMap, tapMap = rmap)
      rankRowsR = _tap_rank_rows(net, activeMeas, _state_to_voltage(xr, slackIdx, nbus), shuntMap === nothing ? nothing : _shunt_b_override(xr, shuntMap), _tap_overlay(xr, rmap, net), cfg.ia_current_floor_A)
      # structural causes and angle anchors at this state: the guard of the
      # start test, in the same order (guard, noise test), on the copies of
      # rmap with the rank rows of this state; the copies carry the flags of
      # the start, so they are cleared first
      for ts in rmap.trafos
        fill!(ts.cause, :none)
        fill!(ts.tree_phase, false)
      end
      _tap_guard_causes!(net, rmap, activeMeas[rankRowsR], slackIdx, withVaOffset)
      # the noise test writes into the copies of rmap only
      _tap_column_freezes!(rmap, Hf, z, xr, nbus, jacEps)
      unanch = Set{Tuple{Int,Int}}((ts.branch, c) for ts in rmap.trafos for c in eachindex(ts.cols) if ts.tree_phase[c])
      live = [k for k in rcand if rmap.trafos[k.j].frozen[k.c] === :none]
      noCause = k -> rmap.trafos[k.j].cause[k.c] === :none
      # the fixed changers first, so a kept-out one becomes
      # determinable only with rank beyond them; inside each group the
      # order of the start test (no structural cause first)
      order = vcat([k for k in live if k.kind === :est && noCause(k)], [k for k in live if k.kind === :est && !noCause(k)], [k for k in live if k.kind === :frz && noCause(k)], [k for k in live if k.kind === :frz && !noCause(k)])
      candCols = Set{Int}(k.col for k in rcand)
      keptR = _tap_rank_kept((0.5 .* (Hf .+ Hb))[rankRowsR, :], [col for col in axes(Hf, 2) if !(col in candCols)], Int[k.col for k in order])
      keptCols = Set{Int}(order[i].col for i in eachindex(order) if keptR[i])
      for k in rcand
        if k.kind === :est && !(k.col in keptCols)
          push!(withdrawList, (k.j, k.c, :not_observable_at_estimate))
        elseif k.kind === :frz && k.col in keptCols
          push!(newly, (k.j, k.c))
        end
      end
    end
    newlyT = sort!(unique(first.(newly)))
    touchedB = Set{Int}(rmap.trafos[j].branch for j in newlyT)
    # the mirror check over every fixed changer that keeps its rank, also
    # on a transformer estimated again below: if that joint estimate fails,
    # its fixed changers are checked here all the same
    for (j, c, txt) in mirrorCheck([(k.j, k.c) for k in rcand if k.kind === :est && k.col in keptCols])
      push!(withdrawList, (j, c, :mirror_step_fits))
      tapExtraNotice[(tapMap.trafos[j].branch, c)] = txt
    end
    tapIteFix += withdraw!(withdrawList)
    # the partners of these withdrawals and of the changers the band
    # verdict of the first estimate put at their transferred step
    markStale!(vcat(Tuple{Int,Int}[(j, c) for (j, c, _) in withdrawList], Tuple{Int,Int}[(j, c) for (j, ts) in enumerate(tapMap.trafos) for c in eachindex(ts.cols) if ts.released[c] && ts.frozen[c] === :out_of_range]))
    # The joint estimate of the newly determinable changers. A newly
    # determinable changer whose sigma_s there is not below
    # _TAP_RECHECK_SIGMA_MAX is not taken over: the state before the joint
    # estimate is put back and solved again, it keeps its transferred step
    # (:weak_at_estimate, the row and the notice keep the estimate and
    # sigma_s), and
    # the joint estimate runs again without it. Each repeat drops at least
    # one changer, so the loop ends. Only the re-estimate gets this bound,
    # never the first round.
    if !isempty(newlyT) && converged
      saved = (x = copy(x), map = tapMap, prior = copy(tapPriorCols), ite = iteDone, gain = lastGain, conJac = lastConJac, lambda = lastLambda, jac = lastJacobian, est = tapMapEst, xest = copy(tapXEst), frozenStates = copy(tapFrozenStates), sigma = copy(tapSigma))
      # set when a pass withdraws a weak changer: a later pass that ends
      # with the state put back still owes the new solve below
      weakWithdrawn = false
      while true
        newlyT = sort!(unique(first.(newly)))
        # per transformer with a newly determinable changer (rmap index):
        # those changers and its other fixed changers that stay
        again = Dict{Int,Vector{Int}}(j => sort!(unique(vcat([c for (jj, c) in newly if jj == j], [k.c for k in rcand if k.j == j && k.kind === :est && !any(wd -> wd[1] == j && wd[2] == k.c, withdrawList)]))) for j in newlyT)
        # the transformers that left the map at the start join it again:
        # out of the Ybus (the overlay stamps them) and out of the frozen
        # list; jOf maps an rmap index to the index of the maps from here
        addIdx = Int[movedIdx[j - nfin] for j in newlyT if j > nfin]
        addTs = TapBranchState[tapFrozenStates[i] for i in addIdx]
        jOf = Dict{Int,Int}(j => (j <= nfin ? j : nfin + findfirst(==(movedIdx[j - nfin]), addIdx)) for j in newlyT)
        againF = Dict{Int,Vector{Int}}(jOf[j] => again[j] for j in newlyT)
        isempty(addTs) || _tap_unstamp!(Ybus, net, TapStateMap(addTs, tapMap.col0, 0, tapMap.impedance))
        tapFrozenStates = TapBranchState[ts for (i, ts) in enumerate(tapFrozenStates) if !(i in addIdx)]
        finTrafos = vcat(tapMap.trafos, addTs)
        estTrafos = vcat(tapMapEst.trafos, addTs)
        reasons = Dict{Tuple{Int,Int},Symbol}()
        # the log lines of the changers estimated again, emitted once the
        # outcome is known: at info level when the joint estimate changes
        # the result (taken over, or a weak changer with its notice), at
        # debug level when it did not converge and the result before it
        # stays unchanged (seen on radial changers kept out at the start)
        detLines = String[]
        for (j, c) in newly
          let ts = finTrafos[jOf[j]]
            reasons[(jOf[j], c)] = ts.frozen[c]
            ts.frozen[c] = :none
            ts.tree_phase[c] = (ts.branch, c) in unanch
            push!(detLines, "runse!: recheck at the estimate: tap regulator $(regName(ts, c)) determinable at the estimate, estimated again with the other free changers of its transformer")
          end
        end
        # map of the joint estimate: the states of the fixed run and the
        # changers to estimate again, at their positions
        colR = tapMap.col0
        valsR = Float64[]
        againCols = Int[]
        newTrafos = TapBranchState[]
        for (jf, ts) in enumerate(finTrafos)
          let cols = zeros(Int, length(ts.cols)), cs = get(againF, jf, Int[])
            for c in eachindex(cols)
              if ts.cols[c] != 0
                push!(valsR, x[ts.cols[c]])
              elseif c in cs
                push!(valsR, ts.spos[c])
                push!(againCols, colR)
              else
                continue
              end
              cols[c] = colR
              colR += 1
            end
            push!(newTrafos, _tap_with_cols(ts, cols, copy(ts.spos)))
          end
        end
        x = vcat(x[1:(tapMap.col0 - 1)], valsR, x[(tapMap.col0 + tapMap.ncols):end])
        tapMap = TapStateMap(newTrafos, tapMap.col0, colR - tapMap.col0, tapMap.impedance)
        tapPriorCols = _tap_prior_columns!(tapMap, net, activeMeas)
        _wls_solve!(false)
        tapIteFix += iteDone
        okR = converged
        weakNow = Tuple{Int,Int,Float64,Union{Nothing,Float64}}[]
        if okR
          # the continuous estimate and its sigma go into the report view
          # of the estimation map (tapMapEst, tapXEst): a changer with a
          # column there gets the new value, one without a column appended
          # behind the state vector (the view is read by the report only)
          jA, gatedA = tapJNow()
          dofA = length(activeMeas) - length(gatedA) - length(x)
          iteA = iteDone
          sigA = lastJacobian === nothing ? Dict{Int,Union{Nothing,Float64}}() : _se_state_sigmas(lastJacobian, wJ, conIdx, [i for i in eachindex(activeMeas) if !(i in gatedA)], againCols)
          # the newly determinable changers without a sigma_s below the bound
          # (none from the covariance counts as not below it)
          for (j, c) in newly
            let col = tapMap.trafos[jOf[j]].cols[c], sg = get(sigA, col, nothing)
              (sg !== nothing && sg < _TAP_RECHECK_SIGMA_MAX) || push!(weakNow, (j, c, x[col], sg))
            end
          end
          okR = isempty(weakNow)
        end
        oorR = Tuple{TapBranchState,Int,Symbol,Bool}[]
        if okR
          # the band verdict of the first round on the joint estimate
          oorR = _tap_mark_out_of_range!(tapMap, x, sigA, Set{Int}(againCols))
          xestR = copy(tapXEst)
          viewTrafos = TapBranchState[]
          for (jf, ts) in enumerate(tapMap.trafos)
            let e = estTrafos[jf], ecols = copy(estTrafos[jf].cols)
              for c in eachindex(ts.cols)
                (ts.cols[c] != 0 && ts.cols[c] in againCols) || continue
                if ecols[c] == 0
                  push!(xestR, x[ts.cols[c]])
                  ecols[c] = length(xestR)
                else
                  xestR[ecols[c]] = x[ts.cols[c]]
                end
                tapSigma[ecols[c]] = get(sigA, ts.cols[c], nothing)
              end
              push!(viewTrafos, _tap_with_cols(e, ecols, e.spos))
            end
          end
          okR = tapFixAndSolve!(againCols)
          tapIteFix += iteDone
          if okR
            tapMapEst = TapStateMap(viewTrafos, tapMapEst.col0, tapMapEst.ncols, tapMapEst.impedance)
            tapXEst = xestR
            # J and dof of the joint estimate and of its fixation, for the
            # fixation statistics of a run that had no fixation before
            tapRecheckEst = (j = jA, gated = gatedA, dof = dofA, ite = iteA, k = length(againCols), jfix = first(tapJNow()))
          end
        end
        if okR || !isempty(weakNow)
          foreach(l -> @info(l), detLines)
        else
          foreach(l -> @debug(l), detLines)
        end
        if okR
          # the mirror check over the changers fixed again (the joint estimate
          # can move the ones checked above)
          w2 = Tuple{Int,Int,Symbol}[]
          for (j, c, txt) in mirrorCheck([(jf, c) for (jf, cs) in againF for c in cs if tapMap.trafos[jf].frozen[c] === :none])
            push!(w2, (j, c, :mirror_step_fits))
            tapExtraNotice[(tapMap.trafos[j].branch, c)] = txt
          end
          tapIteFix += withdraw!(w2)
          # the changers of this joint estimate were estimated with every
          # earlier withdrawal at its transferred step; the withdrawals and
          # band verdicts of this estimate leave their own partners stale
          for (jf, cs) in againF, c in cs
            delete!(stale, (jf, c))
          end
          let jOfB = Dict{Int,Int}(ts.branch => jf for (jf, ts) in enumerate(tapMap.trafos))
            markStale!(vcat(Tuple{Int,Int}[(jOfB[ts.branch], c) for (ts, c, _, _) in oorR], Tuple{Int,Int}[(j, c) for (j, c, _) in w2]))
          end
        else
          # a weak re-estimate, or the joint estimate or its fixation did not
          # converge: the result before it stays, the changers keep their
          # reason, a weak one gets :weak_at_estimate
          isempty(weakNow) && @debug "runse!: recheck at the estimate: the joint estimate of the newly determinable tap regulator(s) did not converge, they keep the transferred step"
          for (j, c, s, sg) in weakNow
            let ts = finTrafos[jOf[j]], ch = ts.set.changers[c], sgText = sg === nothing ? "sigma_s not available" : string("sigma_s ", @sprintf("%.3g", sg), " steps")
              reasons[(jOf[j], c)] = :weak_at_estimate
              tapWeakEst[(ts.branch, c)] = (s, sg)
              weakWithdrawn = true
              tapExtraNotice[(ts.branch, c)] = string("estimated again at step ", @sprintf("%.4g", ch.neutral + s), " with ", sgText)
              @info "runse!: recheck at the estimate: tap regulator $(regName(ts, c)) estimated again at step $(@sprintf("%.4g", ch.neutral + s)) with $(sgText) (bound $(_TAP_RECHECK_SIGMA_MAX)), not taken over (weak_at_estimate), the transferred step $(_tap_notice_step(ch.neutral + ts.s0[c])) stays"
            end
          end
          # a changer the band verdict froze goes back to the state before
          for (ts, c, nt, tp) in oorR
            ts.frozen[c] = :none
            ts.notes[c] = nt
            ts.tree_phase[c] = tp
          end
          for ((jf, c), r) in reasons
            _tap_freeze!(finTrafos[jf], c, r)
          end
          for ts in addTs
            _tap_restamp!(Ybus, net, ts.branch)
          end
          x = saved.x
          tapMap = saved.map
          tapPriorCols = saved.prior
          iteDone = saved.ite
          converged = true
          lastGain = saved.gain
          lastConJac = saved.conJac
          lastLambda = saved.lambda
          lastJacobian = saved.jac
          tapMapEst = saved.est
          tapXEst = saved.xest
          tapFrozenStates = saved.frozenStates
          # a copy: a later pass writes the sigmas of its joint estimate
          # into tapSigma before its fixation, and that fixation can fail
          tapSigma = copy(saved.sigma)
          # again without the weak ones, when a newly determinable one is left
          filter!(k -> !any(w -> w[1] == k[1] && w[2] == k[2], weakNow), newly)
          (!isempty(weakNow) && !isempty(newly)) && continue
          # the final state of a weak withdrawal (in this pass or an earlier
          # one) comes from a new solve with the transferred step, warm from
          # the state put back, so the state and the reported step agree
          if weakWithdrawn
            _wls_solve!(false)
            tapIteFix += iteDone
          end
        end
        break
      end
    end
    # the stale partners estimated again, one pass per round; a round that
    # withdraws a partner leaves the others of its transformer stale again,
    # each round withdraws at least one changer or ends the loop
    while !isempty(stale) && converged
      pairsP = sort!(collect(stale))
      empty!(stale)
      outP = partnerPass!(pairsP)
      markStale!(outP)
    end
    # the ratio partners of every transformer the recheck touched (a
    # changer withdrawn, or one estimated again)
    for ts in tapMap.trafos
      (ts.branch in touchedB || any(c -> ts.frozen[c] in (:not_observable_at_estimate, :mirror_step_fits, :weak_at_estimate, :not_converged_at_estimate), eachindex(ts.cols))) && _tap_recheck_conditional!(ts)
    end
    return nothing
  end
  if tapMap !== nothing
    # J of the CONTINUOUS solution
    tapJBefore, tapGatedBefore = tapJNow()
    tapDofBefore = length(activeMeas) - length(tapGatedBefore) - length(x)
    tapIteEst = iteDone
    tapConvEst = converged
    tapXEst = copy(x)
    if converged
      # sigma of the step from the covariance,
      # on the Jacobian of the last iteration and the rows of J (the Ia rows
      # gated at the estimate left out)
      if lastJacobian !== nothing
        tapSigma = _se_state_sigmas(lastJacobian, wJ, conIdx, [i for i in eachindex(activeMeas) if !(i in tapGatedBefore)], [c for ts in tapMapEst.trafos for c in ts.cols if c != 0])
      end
      # the fixed positions replace the tap states; re-solve from the
      # converged voltages (never from a diverged state: the rounded steps
      # would be meaningless, so a run that did not converge keeps the
      # continuous result, and runse! falls back to the run without tap
      # estimation, see _runse_configured!)
      fixCols = _tap_fixable_columns(tapMapEst)
      tapKRemoved = length(fixCols)
      priorEst = copy(tapPriorCols)
      tapFixed = true
      # an estimate more than _TAP_BAND_SIGMA_MAX sigma_s outside its band
      # is not taken over: frozen with :out_of_range, the fixation puts it
      # at its transferred step (a nearer one goes to the band end step)
      _tap_mark_out_of_range!(tapMapEst, tapXEst, tapSigma, Set{Int}(fixCols))
      tapFixAndSolve!(fixCols)
      tapIteFix = iteDone
      tapJSim, gatedSim = tapJNow()
      # the quantile holds only when both runs gate the same Ia rows (the
      # Delta-J notice below); then a rise above it is the warning case
      tapQuantile = _chi2_quantile95(tapKRemoved)
      tapExceeds = converged && tapKRemoved >= 1 && gatedSim == tapGatedBefore && tapJSim - tapJBefore > tapQuantile
      if tapExceeds && tapKRemoved >= 2
        # sequential fixation from the continuous estimate; with one
        # column it would repeat the simultaneous run
        # the numerics of the last solve (gain, constraint Jacobian,
        # multipliers) go with the state: a rejected sequential run must
        # not leave its own in the result
        simState = (x = copy(x), map = tapMap, prior = copy(tapPriorCols), ite = iteDone, conv = converged, gain = lastGain, conJac = lastConJac, lambda = lastLambda)
        x = copy(tapXEst)
        tapMap = tapMapEst
        tapPriorCols = priorEst
        seqOk = true
        seqIte = 0
        while true
          seqCol = _tap_next_fix(x, tapMap)
          seqCol == 0 && break
          seqOk = tapFixAndSolve!([seqCol])
          seqIte += iteDone
          seqOk || break
        end
        jSeq, gatedSeq = seqOk ? tapJNow() : (Inf, Set{Int}())
        tapJSeq = jSeq
        if seqOk && gatedSeq == gatedSim && jSeq < tapJSim
          tapSequential = true
          tapIteFix = seqIte
        else
          x = simState.x
          tapMap = simState.map
          tapPriorCols = simState.prior
          iteDone = simState.ite
          converged = simState.conv
          lastGain = simState.gain
          lastConJac = simState.conJac
          lastLambda = simState.lambda
        end
      end
      # the recheck at the estimate (tapRecheck! above): once, at the
      # converged state with the fixed steps
      (converged && tapFixed) && tapRecheck!()
      tapPositions = _tap_positions(tapXEst, tapMapEst, x, tapMap)
      converged = converged && tapConvEst
    else
      tapPositions = _tap_positions(tapXEst, tapMapEst, tapXEst, tapMapEst)
    end
  elseif converged && cfg.tap_estimation.enabled && any(ts -> any(c -> ts.released[c] && ts.frozen[c] in _TAP_START_RANK_REASONS, eachindex(ts.cols)), tapFrozenStates)
    # the start test kept every released changer out and the run solved
    # with all of them frozen: the recheck at the estimate runs on this
    # state too, on an empty map of the fixed run (the tap block sits before
    # the link circulations and the PMU offset, _se_insert_tap_states)
    tapMap = TapStateMap(TapBranchState[], length(x) - _se_link_ncirc(net) - (withVaOffset ? 1 : 0) + 1, 0, tapImpedance)
    tapMapEst = tapMap
    tapXEst = copy(x)
    tapRecheck!()
    if tapRecheckEst === nothing
      # no changer became determinable, or its joint estimate did not
      # converge (state and Ybus put back): the run stays without tap states
      tapMap = nothing
      tapMapEst = nothing
      tapXEst = Float64[]
    else
      # the joint estimate is the continuous estimate of this run, its
      # fixation the fixation
      tapFixed = true
      tapConvEst = true
      tapJBefore = tapRecheckEst.j
      tapGatedBefore = tapRecheckEst.gated
      tapDofBefore = tapRecheckEst.dof
      tapIteEst = tapRecheckEst.ite
      tapIteFix -= tapRecheckEst.ite
      tapKRemoved = tapRecheckEst.k
      tapJSim = tapRecheckEst.jfix
      tapQuantile = _chi2_quantile95(tapKRemoved)
      tapPositions = _tap_positions(tapXEst, tapMapEst, x, tapMap)
    end
  end
  # a singular KKT step was solved by the least-squares fallback: say so
  # with the reason, the constrained estimate is then not exact
  kktSingular === nothing || @warn "state estimation: the KKT matrix of the zero-injection constraints was singular, $(kktSingular); the step was solved by the sparse least-squares fallback"

  Vest = _state_to_voltage(x, slackIdx, nbus)
  # the estimate in the model angle reference, for SEResult.voltages and the
  # write-back; Vest itself keeps the slack at 0 for the predictions
  VestModel = vaRefDeg == 0.0 ? Vest : Vest .* cis(refRad)
  vaOffsetRad = _va_offset_from_state(x, withVaOffset, refRad)
  linkCircEst = _se_link_circ(x, net, withVaOffset)
  hfinal = _predict_measurements(activeMeas, net, Vest, Ybus; vaOffsetRad = vaOffsetRad, shuntB = shuntMap === nothing ? nothing : _shunt_b_override(x, shuntMap), tapOverlay = _tap_ov(x), linkCirc = linkCircEst)
  r = z - hfinal
  _wrap_angle_residuals!(r, activeMeas)
  # Ia rows gated at the FINAL state never entered the solve there and must
  # not enter the statistics either: excluded from J and from the redundancy
  finalGated = isempty(iaRows) ? Set{Int}() : _ia_gated!(Vest)
  for i in finalGated
    @info "runse!: current-angle measurement $(activeMeas[i].id) excluded (:ia_below_current_floor, predicted current below $(round(iaThreshold[i]; digits = 2)) A)"
  end
  jval = isempty(finalGated) ? _wls_objective(r, wJ) : _wls_objective(r[[i for i in eachindex(r) if !(i in finalGated)]], wJ[[i for i in eachindex(wJ) if !(i in finalGated)]])
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
      activeObjective = (j = _wls_objective(r[keep], wJ[keep]), dof = ν - length(supIdx), suppressed = length(supIdx))
    end
  end

  # write-back of the FIXED physical tap positions: only on explicit request
  # (update_taps, default false), only from a converged fixed run; a frozen
  # regulator keeps its exact start position (the same
  # never-silently-overwrite protection as update_shunts). It runs BEFORE
  # calcNetLosses! below: before 0.60.0 the losses
  # used the old taps with the voltages of the fixed-tap solve. Branch
  # indices map to the original net by component identity (cID survives
  # the contraction copies). The prediction overlays below stamp scratch
  # branches, so writing the live branch here does not change SbusEst.
  if cfg.update_taps && converged && tapFixed && tapMap !== nothing
    cidToOrigWb = Dict{String,Int}(br.comp.cID => k for (k, br) in enumerate(onet.branchVec))
    for ts in tapMap.trafos
      ob = get(cidToOrigWb, net.branchVec[ts.branch].comp.cID, 0)
      ob == 0 && continue
      # a regulator that stayed a state (:relative_to_frozen) is no fixed
      # step: it keeps its start position like a frozen one
      s_wb = Float64[ts.cols[c] != 0 ? ts.s0[c] : ts.spos[c] for c in eachindex(ts.cols)]
      _tap_apply_position!(onet.branchVec[ob], ts.set, s_wb; s_from = ts.s0, impedance = tapMap.impedance)
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
      Vorig[i] = s == 0 ? complex(something(onet.nodeVec[i]._vm_pu, 1.0), 0.0) : VestModel[s]
      s == 0 && continue
      onet.nodeVec[i]._vm_pu = abs(Vest[s])
      # the net keeps its model angle reference: the slack angle is added
      # in degrees, so the slack gets its stored angle back exactly
      onet.nodeVec[i]._va_deg = rad2deg(angle(Vest[s])) + vaRefDeg
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
      # injections at the estimated step positions (same as the prediction)
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
    _register_se_start!(onet, seVm, seVa, sePinj, seQinj, converged; link_circulation = _se_link_circulation_flows(onet, net, linkCircEst))
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
    # explicit request, and only from a converged run. It runs after the
    # report rows above (on an unmerged net `net === onet`, B_model must be
    # read before the write).
    if cfg.update_shunts && converged
      for (j, si) in enumerate(released)
        osh = onet.shuntVec[prep.shuntOrig[si]]
        bEst = x[_shunt_col0+j-1]
        osh.y_pu_shunt = complex(real(osh.y_pu_shunt), bEst)
        osh.B_shunt = bEst
      end
      # the net update above computed the shunt draw with the model B;
      # refresh it so calcLinkFlowsSE! sees the estimated B (before 0.60.0
      # the link Q flow kept the model B on this path, the island path
      # already refreshed after its shunt mirror)
      (updateNet && !isempty(released)) && updateShuntPowers!(net = onet)
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
        stage, σmod = _robust_stage(activeMeas[i], abs(r[i]), robustK1, robustK2)
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
  if tapMapEst !== nothing || !isempty(tapFrozenStates) || !isempty(tapNoChanger)
    mrids = _transformer_mrids(onet)
    cidToOrig = Dict{String,Int}(br.comp.cID => k for (k, br) in enumerate(onet.branchVec))
    origOf(k) = get(cidToOrig, net.branchVec[k].comp.cID, k)
    rows = NamedTuple[]
    # one row per released regulator: the estimated ones at their estimate
    # and fixed step, the frozen ones (also of transformers that left the
    # map) at their start position, so a released-but-guarded regulator
    # never vanishes silently from the result
    if tapMapEst !== nothing
      for (j, ts) in enumerate(tapMapEst.trafos)
        s_est, s_fix, oor = tapPositions[j]
        ob = origOf(ts.branch)
        append!(rows, _tap_regulator_rows(ts, net.branchVec[ts.branch], ob, get(mrids, ob, ""), s_est, tapFixed ? s_fix : s_est, oor, tapFixed, tapConvEst; sigma_s = [get(tapSigma, ts.cols[c], nothing) for c in eachindex(ts.cols)], voltage_hint = _tap_voltage_hint(net, ts.branch, activeMeas, r), extra_notice = [get(tapExtraNotice, (ts.branch, c), "") for c in eachindex(ts.cols)], weak = [get(tapWeakEst, (ts.branch, c), nothing) for c in eachindex(ts.cols)]))
      end
    end
    for ts in tapFrozenStates
      ob = origOf(ts.branch)
      append!(rows, _tap_regulator_rows(ts, net.branchVec[ts.branch], ob, get(mrids, ob, ""), ts.s0, ts.s0, falses(length(ts.s0)), false, false; voltage_hint = _tap_voltage_hint(net, ts.branch, activeMeas, r), extra_notice = [get(tapExtraNotice, (ts.branch, c), "") for c in eachindex(ts.cols)], weak = [get(tapWeakEst, (ts.branch, c), nothing) for c in eachindex(ts.cols)]))
    end
    for k in tapNoChanger
      ob = origOf(k)
      push!(rows, _tap_no_changer_row(net.branchVec[k], ob, get(mrids, ob, "")))
    end
    # machine transformers: the step calculated from the written-back state
    (machine_rows && converged) && (rows = _tap_machine_rows(rows, onet, measurements))
    # regulators kept out automatically: the step from the control target
    # as a comparison column, from the state written into the net
    (machine_rows && converged) && (rows = _tap_target_rows(rows, onet))
    sort!(rows; by = r -> (r.branch, r.regulator))
    tapEstimates = rows
    # (the update_taps write-back runs before calcNetLosses!, above)
    if tapMapEst !== nothing
      # :offgrid_tap_residual: a band failure that only
      # APPEARS through the fixation is not bad data. The continuous fit
      # (taps still states) passed or undershot the chi-square band, the
      # fixed run fails it :high, so the J jump comes from the rounding to
      # the physical steps: the true position sits between steps, or the
      # step data (increment, neutral position) is wrong.
      offgrid = tapFixed && _band_test_verdict(tapJBefore, tapDofBefore).reason in (:ok, :low) && _band_test_verdict(jval, ν).reason == :high
      # Delta-J notice: J after minus J
      # before with k = the step states the fixed run removed (a regulator
      # whose estimate was already a whole step counts too) and the 95 %
      # quantile of chi2(k). The quantile holds only when both runs gated
      # the same Ia rows (dof_after - dof_before = k); otherwise the notice
      # names both dof. The notice is no acceptance rule; a rise above the
      # quantile started the sequential fixation above and gives a warning.
      # k counts the removed states only: a :relative_to_frozen regulator
      # stays a state in the fixed run
      kRemoved = tapKRemoved
      deltaNotice = ""
      if tapFixed
        dJ = jval - tapJBefore
        ndiff = length(symdiff(tapGatedBefore, Set{Int}(finalGated)))
        deltaNotice = kRemoved == 0 ? string("tap fixation changed J by ", round(dJ; sigdigits = 4), " (no step state removed, every estimated regulator stays a state)") : ndiff == 0 ? string("tap fixation changed J by ", round(dJ; sigdigits = 4), " (k = ", kRemoved, " step state(s) removed, 95 % quantile of chi2(", kRemoved, ") = ", round(_chi2_quantile95(kRemoved); sigdigits = 4), ")") : string("tap fixation changed J by ", round(dJ; sigdigits = 4), " (dof ", tapDofBefore, " before, ", ν, " after, ", ndiff, " Ia row(s) gated differently, no quantile)")
        # defect 3 of the cross check: the rise of the simultaneous
        # fixation above the quantile is named, with the sequential result
        if tapExceeds
          seqText = isnan(tapJSeq) ? "one step state, no sequential fixation" : string("sequential fixation J ", round(tapJSeq; sigdigits = 4), " against simultaneous J ", round(tapJSim; sigdigits = 4), ", kept the ", tapSequential ? "sequential" : "simultaneous", " result")
          deltaNotice = string(deltaNotice, "; the simultaneous fixation raised J by ", round(tapJSim - tapJBefore; sigdigits = 4), ", above the 95 % quantile ", round(tapQuantile; sigdigits = 4), "; ", seqText)
        end
      end
      # j_rise_exceeds_quantile: the simultaneous fixation raised J above
      # the 95 % quantile of chi2(k); sequential: the sequential result was
      # kept; j_simultaneous: J of the simultaneous fixation
      # failure: how this run with tap states failed, the fallback
      # condition of _runse_configured! (reason :none when it did not fail)
      tapFailure = _se_tap_failure_of(converged, jval, Vest, tapFixed, iteDone, tapIteEst, tapIteFix)
      tapFixation = (fixed = tapFixed, j_before = tapJBefore, dof_before = tapDofBefore, j_after = jval, dof_after = ν, iterations_estimation = tapIteEst, iterations_fixation = tapFixed ? tapIteFix : 0, offgrid_residual = offgrid, k_removed = kRemoved, delta_j_notice = deltaNotice, j_rise_exceeds_quantile = tapExceeds, sequential = tapSequential, j_simultaneous = tapFixed ? tapJSim : NaN, fallback_used = false, fallback_reason = :none, failure = tapFailure)
      # a released run that did not converge has no fixation; the fallback
      # in _runse_configured! names it
      if tapFixed
        fixLine = "runse!: tap fixation J $(round(tapJBefore; sigdigits = 4)) (dof $(tapDofBefore)) -> J $(round(jval; sigdigits = 4)) (dof $(ν)); notice: $(deltaNotice)"
        # a rise above the quantile is the warning case (defect 3)
        tapExceeds ? (@warn fixLine) : (@info fixLine)
      end
      offgrid && @info "runse!: :offgrid_tap_residual, the J jump comes from the fixation, not from bad data (true position between physical steps, or wrong step data)"
    end
  end

  conditioning = _se_conditioning(lastGain, lastConJac, cfg.zero_injection)
  ziMultipliers = (isempty(conIdx) || lastLambda === nothing) ? nothing : _se_zero_injection_multipliers(net, activeMeas, activeOrig, conIdx, lastLambda)
  return SEResult(VestModel, converged, iteDone, norm(r), r, jval, ν, _j_within_3sigma_band(jval, ν), withVaOffset ? rad2deg(x[end]) : nothing, shuntEstimates, robustRows, nothing, tapEstimates, tapFixation, nothing, activeObjective, linkCircEst, conditioning, ziMultipliers, flatstart ? :flat : :stored)
end

# The multipliers of the constraint rows `conIdx` of `meas` (on the working
# net `net`), one entry per bus for SEResult.zeroInjectionMultipliers:
# `lambda[k]` belongs to row conIdx[k]. A row with caller position 0 in
# `orig` is the LINKAGG aggregate of a closed-link cluster (kind
# :link_aggregate), every other one a zero-injection row of a passive bus.
function _se_zero_injection_multipliers(net::Net, meas::Vector{Measurement}, orig::Vector{Int}, conIdx::Vector{Int}, lambda::Vector{Float64})::Vector{NamedTuple}
  rows = Dict{Int,Vector{Float64}}()
  kinds = Dict{Int,Symbol}()
  order = Int[]
  for (k, i) in enumerate(conIdx)
    m = meas[i]
    b = m.busIdx::Int   # a constraint row is a bus injection row
    haskey(rows, b) || (push!(order, b); rows[b] = [NaN, NaN]; kinds[b] = orig[i] == 0 ? :link_aggregate : :zero_injection)
    rows[b][m.typ == PinjMeas ? 1 : 2] = lambda[k]
  end
  return NamedTuple[(bus = getCompName(net.nodeVec[b].comp), kind = kinds[b], lambda_p = rows[b][1], lambda_q = rows[b][2]) for b in order]
end

## --- Range check of the voltage measurements -----------------------------
# The plausible band of a measured bus voltage magnitude in pu. A VmMeas
# outside it, at 0 kV or not finite, is no reading of an energized bus:
# either the switching state leaves the bus dead (then 0 kV is right) or
# the value is a gross error. `_se_range_check` sorts the two apart before
# every WLS of runse!, validate_measurements and runse_diagnostics.
const _SE_VM_RANGE_PU = (0.5, 1.5)

_se_vm_out_of_range(v::Float64)::Bool = !isfinite(v) || v < _SE_VM_RANGE_PU[1] || v > _SE_VM_RANGE_PU[2]

## the buses the switching state energizes: members of an AC island
## (closed branches; closed links join islands and link-fed isolated buses)
## that holds a source, i.e. an island the power flow would give a
## reference (`detect_ac_islands` with `promote_generators`, the rule of the
## solve). A bus typed Isolated without a closed link to such an island and
## every bus of an island without a source are not energized.
function _se_energized_buses(net::Net)::Set{Int}
  rows = detect_ac_islands(net; promote_generators = true).rows
  n = length(net.nodeVec)
  parent = collect(1:n)
  function findroot(i::Int)
    while parent[i] != i
      parent[i] = parent[parent[i]]
      i = parent[i]
    end
    return i
  end
  function unite!(a::Int, b::Int)
    ra = findroot(a)
    rb = findroot(b)
    ra != rb && (parent[rb] = ra)
    return nothing
  end
  for row in rows, b in row.buses
    unite!(first(row.buses), b)
  end
  for l in net.linkVec
    l.status == 1 || continue
    f = Int(l.fromBus)
    t = Int(l.toBus)
    (1 <= f <= n && 1 <= t <= n) && unite!(f, t)
  end
  fed = Set{Int}(findroot(Int(row.chosen_ref_bus)) for row in rows if row.chosen_ref_bus != 0)
  return Set{Int}(b for b in 1:n if findroot(b) in fed)
end

## the substitute of a voltage measurement that fails the range check on an
## energized bus: the voltage of the imported state at the bus
## (`net._sv_vm_pu`, the SvVoltage of a CGMES import) and :imported_state,
## else the nominal voltage 1.0 pu and :nominal. A stored value outside the
## band is no substitute either.
function _se_range_substitute(net::Net, bus::Int)::Tuple{Float64,Symbol}
  v = get(net._sv_vm_pu, bus, NaN)
  _se_vm_out_of_range(v) && return (1.0, :nominal)
  return (v, :imported_state)
end

## whether measurement `m` sits at bus `b`: at the bus itself, at a branch
## or a link with an end at it
function _se_measurement_at_bus(net::Net, m::Measurement, b::Int)::Bool
  m.busIdx == b && return true
  if m.branchIdx !== nothing && 1 <= m.branchIdx <= length(net.branchVec)
    br = net.branchVec[m.branchIdx]
    (Int(br.fromBus) == b || Int(br.toBus) == b) && return true
  end
  if m.linkIdx !== nothing && 1 <= m.linkIdx <= length(net.linkVec)
    l = net.linkVec[m.linkIdx]
    (Int(l.fromBus) == b || Int(l.toBus) == b) && return true
  end
  return false
end

const _SE_RANGE_ROW = NamedTuple{(:measurement_index, :id, :bus, :value, :reason, :substitute, :source, :sigma, :sigma_substitute, :deactivated),Tuple{Int,String,Int,Float64,Symbol,Float64,Symbol,Float64,Float64,Vector{Int}}}

## The range check, the step before every WLS of the library. Per active
## VmMeas outside _SE_VM_RANGE_PU (0 kV and non-finite values included):
## 1. topology first: on a bus the switching state does not energize
##    (`_se_energized_buses`) 0 kV is right; the row and every other active
##    row at the bus, its branches and its links are set inactive, so the
##    bus and its branches stay out of the estimation (reason :deenergized,
##    `deactivated` lists the rows);
## 2. on an energized bus the value is a gross error: it is not eliminated
##    but replaced by a pseudo measurement at `_se_range_substitute` with
##    its sigma times state_estimation.bad_data.substitute_sigma_factor
##    (reason :range_check_substituted).
## Returns the checked copy of `measurements` (same order and length, every
## index stays the caller's) and one row per finding. It logs nothing: the
## public entry points log once per call (`_se_log_range_check`), the SE
## service writes the same lines into run.log.
function _se_range_check(net::Net, measurements::Vector{Measurement})::Tuple{Vector{Measurement},Vector{_SE_RANGE_ROW}}
  rows = _SE_RANGE_ROW[]
  hits = [k for (k, m) in enumerate(measurements) if m.active && m.typ == VmMeas && m.busIdx !== nothing && _se_vm_out_of_range(m.value)]
  isempty(hits) && return measurements, rows
  out = copy(measurements)
  energized = _se_energized_buses(net)
  factor = state_estimation_config().bad_data.substitute_sigma_factor
  for k in hits
    # a row a dead bus found earlier already took out
    out[k].active || continue
    m = measurements[k]
    b = Int(m.busIdx)
    if b in energized
      v, src = _se_range_substitute(net, b)
      out[k] = Measurement(typ = m.typ, value = v, sigma = m.sigma * factor, active = true, busIdx = m.busIdx, branchIdx = m.branchIdx, direction = m.direction, id = m.id, linkIdx = m.linkIdx, zero_injection = m.zero_injection)
      push!(rows, _SE_RANGE_ROW((k, m.id, b, m.value, :range_check_substituted, v, src, m.sigma, m.sigma * factor, Int[])))
    else
      off = Int[]
      for j in eachindex(out)
        (out[j].active && _se_measurement_at_bus(net, out[j], b)) || continue
        out[j] = _set_measurement_active(out[j], false)
        push!(off, j)
      end
      push!(rows, _SE_RANGE_ROW((k, m.id, b, m.value, :deenergized, NaN, :none, m.sigma, NaN, off)))
    end
  end
  return out, rows
end

## one text line per row of the range check; the library log and run.log of
## the SE service carry the same text
function _se_range_check_lines(net::Net, rows::AbstractVector)::Vector{String}
  names = _bus_name_by_idx(net)
  lines = String[]
  for r in rows
    bus = get(names, r.bus, string(r.bus))
    head = string("range check: voltage measurement ", r.id, " at bus ", bus, " is ", round(r.value; sigdigits = 6), " pu, outside ", _SE_VM_RANGE_PU[1], " to ", _SE_VM_RANGE_PU[2], " pu")
    if r.reason === :range_check_substituted
      src = r.source === :imported_state ? "imported state" : "nominal voltage"
      push!(lines, string(head, "; the bus is energized: replaced by ", round(r.substitute; sigdigits = 6), " pu (", src, "), sigma ", round(r.sigma; sigdigits = 4), " -> ", round(r.sigma_substitute; sigdigits = 4), " (range_check_substituted)"))
    else
      push!(lines, string(head, "; the bus is not energized by the switching state: the bus and its branches stay out of the estimation, ", length(r.deactivated), " measurement(s) inactive (deenergized)"))
    end
  end
  return lines
end

## the log of the range check of one public call: a warning per substituted
## row (a gross error), an info line per dead bus (0 kV is right there)
function _se_log_range_check(net::Net, rows::AbstractVector, caller::AbstractString)
  isempty(rows) && return nothing
  for (r, line) in zip(rows, _se_range_check_lines(net, rows))
    r.reason === :range_check_substituted ? (@warn string(caller, ": ", line)) : (@info string(caller, ": ", line))
  end
  return nothing
end

## the counter of the range check, as the diagnostics results carry it
_se_range_check_summary(rows::AbstractVector) = (substituted = count(r -> r.reason === :range_check_substituted, rows), deenergized = count(r -> r.reason === :deenergized, rows), rows = rows)

## a diagnostics report of `_validate_measurements` with the range check of
## its call: the counter as field `range_check`, and every ranking row
## (`measurement_ranking`, `suspicious_measurements`,
## `largest_normalized_residual`) with the field `range_check`,
## :range_check_substituted for a substituted row, :none for the others
function _se_report_with_range_check(report, rows::AbstractVector)
  sub = Set{Int}(r.measurement_index for r in rows if r.reason === :range_check_substituted)
  mark(row) = merge(row, (range_check = row.measurement_index in sub ? :range_check_substituted : :none,))
  lnr = report.largest_normalized_residual
  return merge(report, (measurement_ranking = NamedTuple[mark(r) for r in report.measurement_ranking], suspicious_measurements = NamedTuple[mark(r) for r in report.suspicious_measurements], largest_normalized_residual = lnr === nothing ? nothing : mark(lnr), range_check = _se_range_check_summary(rows)))
end

# The public entry points take NO settings (#381): the estimator fetches
# its configuration from the registry itself. A run with other settings
# installs them for its duration with with_state_estimation_config or
# with_sparlectra_config; the service does that with the resolved
# configuration of the case.
"""
    runse!(net, measurements) -> SEResult
    runse!(net) -> SEResult

Run the nonlinear weighted-least-squares state estimation, island by island;
`runse!(net)` uses the measurements stored on the net (`net.measurements`).
The estimated state is written back when `state_estimation.update_net` is set
(default `true`). The estimator reads every setting from the active
configuration (`state_estimation_config()`, #381); nothing is passed per
call. A run with other settings wraps the call:

    with_state_estimation_config(max_iter = 12, update_net = false) do
      runse!(net)
    end

With `state_estimation.topology_precheck` (default `true`) an advisory
topology check runs first: its findings are warned once per run and attached
to the result, and the estimation always proceeds. The tap row of a
released changer on a branch of a mesh without P measurement gets that in
its notice (the step may not be unique); with the check off the rows carry
no such notice.

Before the WLS every active voltage measurement outside 0.5 to 1.5 pu (0 kV
and non-finite values included) goes through the range check: on a bus the
switching state does not energize (isolated, or in an island without a
source) the row and the rows of the bus, its branches and its links are set
inactive for the run (log reason `deenergized`); on an energized bus the
value is replaced by a pseudo measurement at the imported state of the bus
(the `SvVoltage` of a CGMES import, else 1.0 pu) with its sigma times
`state_estimation.bad_data.substitute_sigma_factor` (default 3; warning
with reason `range_check_substituted`). The caller's vector is not changed.

State representation:
- bus voltage angles for all non-slack buses (radians), the slack at 0
- bus voltage magnitudes for all buses (p.u.)
- optional PMU reference-angle offset α (radians), appended automatically
  when active `VaMeas` measurements exist and `pmu_ref_offset = :auto`

The start is the flat profile (`state_estimation.flatstart`, default
`true`: every |V| 1 pu, every angle 0, the slack magnitude from the net;
with `state_estimation.start_mode.rotate_behind_pst` buses behind phase
shifters at their shift sum, as in the power flow flat start) or the
stored state of the net, taken against the stored slack angle. The result
voltages, the write-back and the angle predictions are in the angle
reference of the model: the slack keeps its stored angle.

PMU angle measurements (`VaMeas`, degrees) are modeled as
`z = θ_i + θ_ref + α + e`, with `θ_ref` the stored slack angle: α is the
offset of the PMU time base against the model angle reference. With
`pmu_ref_offset = :off` the offset state is omitted and PMU angles are
taken in the model angle reference.
"""
function runse!(net::Net, measurements::Vector{Measurement})
  cfg = state_estimation_config()
  _tap_warn_unknown_exclusions(net, cfg.tap_estimation, "runse!")
  # range check of the voltage measurements before the WLS
  checked, rc = _se_range_check(net, measurements)
  _se_log_range_check(net, rc, "runse!")
  res = _runse_configured!(net, checked, cfg)
  # one warning per regulator whose estimate left its band. Within 2 sigma_s
  # of the band end step the fixation sets it to that step and the status
  # stays :estimated, so without the warning the result looks normal;
  # farther out the status is :out_of_range and the transferred step stays.
  # The rows of a fallback come from the run without tap states and carry
  # no flag.
  if res.tapEstimates !== nothing
    for t in res.tapEstimates
      t.out_of_range && @warn string("runse!: ", _tap_out_of_range_line(t))
    end
  end
  return res
end

runse!(net::Net) = runse!(net, Measurement[m for m in net.measurements])

"""
    validate_measurements(net, measurements) -> NamedTuple

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

The range check of `runse!` runs first (voltage measurements outside 0.5
to 1.5 pu). The report carries its counter `range_check = (substituted,
deenergized, rows)`, and every ranking row the field `range_check`
(`:range_check_substituted` for a replaced row, `:none` otherwise).
"""
function validate_measurements(net::Net, measurements::Vector{Measurement})
  checked, rc = _se_range_check(net, measurements)
  _se_log_range_check(net, rc, "validate_measurements")
  return _se_report_with_range_check(_validate_measurements(net, checked, _tap_exclusion_map(net, state_estimation_config().tap_estimation)), rc)
end

## the body of `validate_measurements`, with the transformers the tap
## estimation keeps out resolved on the caller's net (`tapExcl`,
## `_tap_exclusion_map`): the solve below runs on the link-contracted copy,
## whose fused and renamed buses would otherwise keep out other
## transformers (machine test, mRID match) than `runse!` on the same net
function _validate_measurements(net::Net, measurements::Vector{Measurement}, tapExcl::Dict{String,Symbol})
  # every setting from the registry (#381); the residual-correlation report
  # and the bad-data thresholds included
  base = state_estimation_config()
  maxIte = base.max_iter
  tol = base.tol
  flatstart = base.flatstart
  jacEps = base.jac_eps
  normalizedThreshold = base.k_eliminate
  pmuRefOffset = base.pmu_ref_offset
  wiiThreshold = _SE_WII_THRESHOLD
  reportResidualCorrelation = base.report_residual_correlation
  robust = base.robust
  robustStartIteration = base.robust_start_iteration
  robustMode = base.robust_mode
  robustK1 = base.robust_k1
  robustK2 = base.robust_k2
  kSuppress = base.k_suppress
  suppressionSigma = base.suppression_sigma
  # island-wise nets: run the full diagnostics per measured island and merge
  # (measurement_index stays the position in the CALLER's vector, so the
  # sequential elimination in runse_diagnostics works unchanged)
  part = _se_island_partition(net, measurements)
  if length(part.rows) > 1
    return _validate_measurements_islands(net, measurements, part, tapExcl)
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
  result = _runse_configured!(netLocal, activeMeas, _se_config_with(base; update_net = false, topology_precheck = false); tap_excl = tapExcl)

  # tap unification: the estimator's final voltages belong to the FIXED tap
  # positions, but write-back is usually off, so the diagnostics copy still
  # carries the model taps; stamping the fixed positions here keeps the
  # residuals (and the Ia gate below) free of a fake tap discrepancy at the
  # released transformers (`_tap_stamp_end_steps!`)
  _tap_stamp_end_steps!(netLocal, result, base.tap_estimation.impedance_correction)

  # dynamic part of the IaMeas gate: rows below the predicted-current
  # threshold at the ESTIMATED state leave the diagnostic set too
  if any(m -> m.typ == IaMeas, activeMeas)
    cfg0 = base
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
  # result.voltages are in the model angle reference, the state keeps the
  # slack at 0 (_se_slack_ref_deg)
  refRadLocal = deg2rad(_se_slack_ref_deg(netLocal, slackIdx))
  for i in eachindex(netLocal.nodeVec)
    if i != slackIdx
      p += 1
      x[p] = angle(result.voltages[i]) - refRadLocal
    end
    x[nθ+i] = abs(result.voltages[i])
  end
  if withVaOffset
    x[end] = deg2rad(result.vaRefOffsetDeg)
  end
  # the estimated ring circulations of the link rows (same table: netLocal
  # is a copy of the working net the run solved on)
  if _se_link_ncirc(netLocal) > 0 && result.linkCirculations === nothing
    @warn "validate_measurements: the run returned no circulation states for the measured coupler rings; their link rows are judged at zero circulation"
  end
  _se_link_circ_set!(x, netLocal, withVaOffset, result.linkCirculations)


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

  # A regulator that stayed a state in the fixed run (note
  # :relative_to_frozen, not fixed) is a state of the diagnostics Jacobian
  # too, at its continuous value (stamped above), after the B states and
  # before the offset, as in the run. Without its column rn, wii and Omega
  # would describe a state one column short of the run's (sum of wii =
  # dof + 1 per regulator). The
  # map releases every changer of the stamped branches, the others are
  # taken out again at their stamped positions (`_tap_repack`).
  tapMap = nothing
  relRows = result.tapEstimates === nothing || result.tapFixation === nothing || !result.tapFixation.fixed ? NamedTuple[] : NamedTuple[tr for tr in result.tapEstimates if tr.note === :relative_to_frozen && !tr.fixed && tr.status === :estimated]
  if !isempty(relRows)
    keepCols = Set{Tuple{Int,Int}}((tr.branch, tr.regulator) for tr in relRows)
    tapMap, _, _ = _build_tap_state_map(netLocal, nθ + nbus + length(releasedRows) + 1, base.tap_estimation.impedance_correction)
    if tapMap !== nothing
      _tap_unstamp!(Ybus, netLocal, tapMap)
      for ts in tapMap.trafos, c in eachindex(ts.cols)
        # any reason but :none takes the column out; this local map is used
        # for the Jacobian only and never reported
        (ts.cols[c] != 0 && !((ts.branch, c) in keepCols)) && (ts.frozen[c] = :fixed_in_run)
      end
      tapMap, _ = _tap_repack(tapMap, netLocal, Ybus)
    end
    if tapMap !== nothing
      sinit = _tap_s_init(tapMap)
      x = _se_insert_tap_states(x, sinit, netLocal, withVaOffset)
    end
  end

  H, h = _measurement_jacobian_fd(activeMeas, netLocal, x, slackIdx, nbus, Ybus; eps = jacEps, withVaOffset = withVaOffset, shuntMap = shuntMap, tapMap = tapMap)
  z = _measurement_vector(activeMeas)
  # constraint rows (state_estimation.zero_injection = constraint) enter
  # the residual covariance at the zero-injection weight, the limit they
  # stand for (see _constraint_diagnostic_weights)
  conIdx = _zero_injection_constraint_rows(activeMeas, base.zero_injection)
  w = _constraint_diagnostic_weights(_weight_vector(activeMeas), conIdx)
  r = z - h
  _wrap_angle_residuals!(r, activeMeas)
  # K needs the full Omega up to the dense state cap, so there a
  # correlation request pins the dense path; otherwise the Takahashi
  # selected inverse takes over above the size threshold (with an
  # automatic fallback on any guard failure)
  rd = _residual_diagnostics(H, r, w; need_full_omega = reportResidualCorrelation, minStates = base.takahashi_min_states)
  # a constraint has no statistical residual: normalized residual zero,
  # never suspicious, marked as constraint in the ranking
  rd.rn[conIdx] .= 0.0
  kmax = nothing
  if reportResidualCorrelation
    if rd.omega !== nothing
      kmax = _residual_correlation_max(rd.omega)
    elseif size(H, 1) <= _SE_FULL_OMEGA_MAX_M
      # above the dense state cap: K column block by column block on the
      # sparse factorization (the m cap stays the run-time budget; beyond
      # it _residual_diagnostics already warned that the report is skipped)
      kmax = _residual_correlation_max_blockwise(H, w)
    end
  end
  ranking = _build_measurement_suspicion_report(activeMeas, activeIdx, r, rd.rn; normalizedThreshold = normalizedThreshold, wii = rd.wii, wiiThreshold = wiiThreshold, kmax = kmax, constraintRows = conIdx)

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

validate_measurements(net::Net) = validate_measurements(net, Measurement[m for m in net.measurements])

## Internal helper: measurements protected from bad-data elimination.
## Zero-injection rows (the `zero_injection` flag, set by the code that
## creates them, also on the LINKAGG aggregate of a cluster of them) encode network
## structure, not telemetry, and any near-exact measurement (sigma <= 1e-6)
## is a constraint in disguise; eliminating either would remove model
## knowledge instead of a faulty meter. Derived shunt pseudo-measurements
## (id prefix "SHDERIV", SE phase 2 case B) are protected too: eliminating
## one would silently mask its source measurements (bay current, Vm).
@inline function _elimination_protected(m::Measurement)::Bool
  return m.zero_injection || startswith(m.id, "SHDERIV") || m.sigma <= 1e-6
end

## _se_report_tap_fallback(report) -> Vector{Int}: the islands of a
## `validate_measurements` report whose run did not converge with tap
## estimation and was solved again without it (`tapFixation.fallback_used`);
## a single-island report gives [1] or []. runse_diagnostics compares a
## control solve with its base run by this list, the SE service reads it on
## the base run.
function _se_report_tap_fallback(report)::Vector{Int}
  _fellBack(r) = r !== nothing && r.tapFixation !== nothing && r.tapFixation.fallback_used
  if hasproperty(report, :islands)
    return [ir.island for ir in report.islands if ir.report !== nothing && _fellBack(ir.report.result)]
  end
  return _fellBack(report.result) ? [1] : Int[]
end

"""
    runse_diagnostics(net, measurements) -> NamedTuple

Extended diagnostics workflow around `validate_measurements` with sequential
bad-data elimination (identification stage 1):

While the χ²-like 3σ band test fails, a suspicious measurement
(`|rn| >= normalizedThreshold`) exists, and fewer than `max_eliminations`
measurements have been removed, the top suspicious measurement is
deactivated and the diagnostics rerun. Zero-injection pseudo-measurements
(id prefix `ZI`), near-exact measurements (`sigma <= 1e-6`) and the pseudo
measurements of the range check are never eliminated.

The range check of `runse!` runs once before the first pass (voltage
measurements outside 0.5 to 1.5 pu: replaced on an energized bus, reason
`range_check_substituted`; out of the estimation with the bus and its
branches on a bus that is not energized, reason `deenergized`); every pass
works on the checked set. The default budget `max_eliminations = 5` is also
the round limit of the loop (one row per round).

Returns `(diagnostics, rerun, eliminations, reverted, stop_reason, final_diagnostics, topology_findings, range_check)`:
- `diagnostics`: the initial `validate_measurements` report
- `eliminations`: one trace row per elimination
  `(elimination, measurement_index, id, typ, normalized_residual_before,
  residual_before, wii, skipped_unlocalizable, objective_before,
  objective_after)`; `residual_before` is the residual of the row in its
  measurement unit at the state before its removal.
  `wii` is the localizability of the removed row and
  `skipped_unlocalizable` counts the suspects that were passed over in that
  round BECAUSE they are not localizable, so the trace shows why the
  largest normalized residual was not the one removed.
- `stop_reason`: `:consistent` (band test passed), `:band_low` (J below the
  band: the sigmas overstate the errors, not bad data, nothing eliminated),
  `:no_suspicious_left`,
  `:no_localizable_suspect` (suspects remain, but none of them is
  localizable: the data cannot say which row is wrong, which points at a
  measurement gap around those rows rather than at bad data),
  `:max_eliminations`, or `:not_converged` (the base run or a control
  solve did not converge; a control solve that converged only in the
  fallback without tap estimation, where the base run estimated the taps,
  counts as not converged: the row goes back to active and is reported in
  `reverted`)
- `reverted`: the elimination attempt that was taken back because its
  control solve did not converge, `nothing` otherwise
- `final_diagnostics`: the report after the last elimination (equal to
  `diagnostics` when nothing was eliminated)
- `rerun`: backward-compatible view of the FIRST elimination
  (`deactivated_measurement_index`, `deactivated_measurement_id`,
  `diagnostics`), `nothing` when no elimination happened
- `topology_findings`: the station clusters of the topology classification,
  `nothing` unless the loop stopped at `:max_eliminations` with J still
  above the band
- `range_check`: the counter of the range check, `(substituted,
  deenergized, rows)`, one row per finding (`measurement_index`, `id`,
  `bus`, `value`, `reason`, `substitute`, `source`, `sigma`,
  `sigma_substitute`, `deactivated`)

Backward compatibility: `deactivate_and_rerun = true` is kept as an alias for
`max_eliminations = 1` (exactly one elimination pass, result also in `rerun`).
Unlike the old single pass, the loop only eliminates while the band test
fails; a suspicious measurement in an already consistent set is reported but
not removed.
"""
function runse_diagnostics(net::Net, measurements::Vector{Measurement})
  # every setting from the registry (#381): the elimination budget
  # (state_estimation.max_eliminations) included
  se_cfg = state_estimation_config()
  _tap_warn_unknown_exclusions(net, se_cfg.tap_estimation, "runse_diagnostics")
  normalizedThreshold = se_cfg.k_eliminate
  wiiThreshold = _SE_WII_THRESHOLD
  effMax = se_cfg.max_eliminations
  effMax >= 0 || error("runse_diagnostics: state_estimation.max_eliminations must be >= 0")

  # the range check once for the whole pass: every control solve below
  # works on the checked set, and a substituted row is never eliminated
  # (it is the replacement of a gross error, and the bus keeps it)
  checked, rc = _se_range_check(net, measurements)
  _se_log_range_check(net, rc, "runse_diagnostics")
  substituted = Set{Int}(r.measurement_index for r in rc if r.reason === :range_check_substituted)
  tapExcl = _tap_exclusion_map(net, se_cfg.tap_estimation)
  diagnose(ms) = _se_report_with_range_check(_validate_measurements(net, ms, tapExcl), rc)
  base = diagnose(checked)
  baseFellBack = _se_report_tap_fallback(base)

  meas2 = copy(checked)
  current = base
  trace = NamedTuple[]
  firstElimination = nothing
  stopReason = :consistent
  revertedAttempt = nothing

  while true
    # one abort point per elimination pass (every pass solves again)
    sparlectra_check_abort()
    if !current.converged
      stopReason = :not_converged
      break
    end
    if current.objective.within_3sigma
      stopReason = :consistent
      break
    end
    # only a :high band failure (J above the band) says bad data; a :low
    # failure means the sigmas overstate the errors, and removing a row can
    # only lower J further, so the elimination does not start there (until
    # 0.32.1 it ran on every band failure and reported the budget or the
    # missing suspects as the reason)
    bandReason = hasproperty(current.objective, :reason) ? current.objective.reason : :high
    if bandReason == :low
      stopReason = :band_low
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
    suspects = [row for row in current.measurement_ranking if row.suspicious && row.measurement_index >= 1 && !_elimination_protected(meas2[row.measurement_index]) && !(row.measurement_index in substituted)]
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
    previous = current
    meas2[idx] = _set_measurement_active(meas2[idx], false)
    current = diagnose(meas2)
    if !current.converged || !issubset(_se_report_tap_fallback(current), baseFellBack)
      # An elimination is confirmed by its control solve. When that solve
      # does not converge, the removal is only claimed, not shown (seen
      # on the MiniGrid: J 870 to 758897 after the third
      # elimination): the row goes back to active, the diagnostics stand at
      # the state before the attempt, the attempt is reported by name, and
      # the loop stops with not_converged.
      # The same holds when the control solve did not converge with tap
      # estimation in an island where the base run did, and runse! solved
      # it again without tap estimation: it then converges on a different
      # state model than the base run, and the removal is not shown for
      # the model the pass started on.
      meas2[idx] = _set_measurement_active(meas2[idx], true)
      revertedAttempt = (elimination = length(trace) + 1, measurement_index = idx, id = target.id, typ = target.typ,
                         normalized_residual_before = target.normalized_residual, wii = target.wii, objective_before = objBefore)
      current = previous
      stopReason = :not_converged
      break
    end
    push!(trace, (elimination = length(trace) + 1, measurement_index = idx, id = target.id, typ = target.typ,
                  normalized_residual_before = target.normalized_residual, residual_before = target.residual, wii = target.wii,
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
    reps, mems = _station_map(net)
    clusterMin = state_estimation_config().topology_cluster_min
    # station -> surviving suspects (a branch row counts for both terminals)
    byStation = Dict{Int,Vector{NamedTuple}}()
    for row in current.measurement_ranking
      (row.suspicious && row.measurement_index >= 1) || continue
      m = meas2[row.measurement_index]
      stations = Int[]
      m.busIdx !== nothing && push!(stations, reps[m.busIdx])
      if m.linkIdx !== nothing && 1 <= m.linkIdx <= length(net.linkVec)
        lk = net.linkVec[m.linkIdx]
        push!(stations, reps[Int(lk.fromBus)], reps[Int(lk.toBus)])
      end
      if m.branchIdx !== nothing && 1 <= m.branchIdx <= length(net.branchVec)
        br = net.branchVec[m.branchIdx]
        push!(stations, reps[Int(br.fromBus)], reps[Int(br.toBus)])
      end
      for st in unique(stations)
        push!(get!(byStation, st, NamedTuple[]), row)
      end
    end
    # stage-1 agreement: a precheck finding at the
    # same station is extra evidence, noted but not severity-changing
    pre = validate_topology(net, checked)
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
    # different station to the front (found on a CI run).
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

  return (diagnostics = base, rerun = firstElimination, eliminations = trace, reverted = revertedAttempt, stop_reason = stopReason, final_diagnostics = current, topology_findings = topoFindings, range_check = _se_range_check_summary(rc))
end

runse_diagnostics(net::Net) = runse_diagnostics(net, Measurement[m for m in net.measurements])

@inline function _base_diag_report(diag)
  return hasproperty(diag, :diagnostics) ? diag.diagnostics : diag
end

# the range check counter of a diagnostics result or report, `nothing` for
# a result built without one
_se_diag_range_check(diag) = hasproperty(diag, :range_check) ? diag.range_check : (hasproperty(_base_diag_report(diag), :range_check) ? _base_diag_report(diag).range_check : nothing)

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

# Flag of one ranking row. The ranking comes from the pass with every row
# active, so one gross error smears into its neighbours and pushes them over
# the threshold too; only the rows the sequential elimination removed are
# bad data ("BAD"), the others above the threshold are "over k". A row the
# range check replaced is "substituted" (never eliminated).
function _se_ranking_flag(row, eliminated::AbstractSet{Int})::String
  get(row, :constraint, false) && return "constraint"
  get(row, :range_check, :none) === :range_check_substituted && return "substituted"
  row.measurement_index in eliminated && return "BAD"
  return row.suspicious ? "over k" : "OK"
end

"""
    print_se_diagnostics(io, diag; topN=10)

Pretty-print diagnostics from `validate_measurements` or `runse_diagnostics`
including:
- explanation of `global_consistency`
- tabular measurement ranking
- marker per measurement: `BAD` for a row the sequential elimination
  removed, `over k` for a row above the threshold that stayed (smearing of
  a gross error elsewhere), `OK`, or `constraint`
- optional rerun comparison if present
"""
function print_se_diagnostics(diag; io = stdout, topN::Int = 10, format::Symbol = :plain)
  base = _base_diag_report(diag)
  summary = summarize_se_diagnostics(diag)

  nrows = min(topN, length(base.measurement_ranking))
  eliminated = hasproperty(diag, :eliminations) ? Set{Int}(t.measurement_index for t in diag.eliminations) : Set{Int}()
  # wii/localizable columns exist on 0.10 ranking rows; correlation columns
  # carry values only when the K-matrix report was enabled
  hasWii = !isempty(base.measurement_ranking) && hasproperty(first(base.measurement_ranking), :wii)
  showK = hasproperty(base, :correlation_enabled) && base.correlation_enabled
  if format == :markdown
    println(io, "## State-estimation diagnostics")
    println(io)
    println(io, "- **Converged:** $(summary.converged)")
    println(io, "- **Global consistency:** $(summary.global_consistency) ($(_global_consistency_label(summary)))")
    @printf(io, "- **Objective J:** %.6f (dof=%d, z=%.3f, within_3sigma=%s)\n", summary.objective.value, summary.objective.dof, summary.objective.zscore, string(summary.objective.within_3sigma))
    if hasproperty(summary.objective, :z_wh)
      @printf(io, "- **Wilson-Hilferty:** z_wh=%.3f, reason=%s%s\n", summary.objective.z_wh, string(summary.objective.reason), summary.small_redundancy ? " (small redundancy: nu < 30, reduced test power)" : "")
    end
    println(io, "- **Interpretation:** $(summary.reason)")
    println(io, "- **Suspicious measurements:** $(summary.suspicious_count) / $(summary.total_measurements)")
    rcs = _se_diag_range_check(diag)
    (rcs !== nothing && !isempty(rcs.rows)) && println(io, "- **Range check of the voltage measurements:** $(rcs.substituted) replaced by a pseudo measurement (`range_check_substituted`), $(rcs.deenergized) at a bus that is not energized (`deenergized`): ", join(("`$(r.id)` ($(r.reason))" for r in rcs.rows), ", "))
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
        # ids sit in code spans in every Markdown table: a Markdown renderer
        # (the Web UI artifact viewer) otherwise reads the underscores of
        # `ZI_QINJ_bus_15` as emphasis and shows `ZIQINJbus_15`
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
      # a constraint row (zero_injection = constraint) is no measurement
      flag = _se_ranking_flag(row, eliminated)
      if hasWii
        loc = row.localizable ? "yes" : "no"
        if showK
          kwarn = row.correlation_warning ? " (!)" : ""
          @printf(io, "| %d | %s | %s | %.5f | %.5f | %.3f | %s | %.3f%s | %s |\n", row.measurement_index, string("`", row.id, "`"), string(row.typ), row.residual, row.normalized_residual, row.wii, loc, row.max_abs_correlation, kwarn, flag)
        else
          @printf(io, "| %d | %s | %s | %.5f | %.5f | %.3f | %s | %s |\n", row.measurement_index, string("`", row.id, "`"), string(row.typ), row.residual, row.normalized_residual, row.wii, loc, flag)
        end
      else
        @printf(io, "| %d | %s | %s | %.5f | %.5f | %s |\n", row.measurement_index, string("`", row.id, "`"), string(row.typ), row.residual, row.normalized_residual, flag)
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
    rcs = _se_diag_range_check(diag)
    if rcs !== nothing && !isempty(rcs.rows)
      @printf(io, "Range check of the voltage measurements: %d substituted, %d deenergized\n", rcs.substituted, rcs.deenergized)
      for r in rcs.rows
        println(io, "  ", r.id, " (", r.reason, ")")
      end
    end
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
      # a constraint row (zero_injection = constraint) is no measurement
      flag = _se_ranking_flag(row, eliminated)
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

  # constraint rows (zero_injection = constraint) are no measurements and
  # sort last in the ranking (normalized residual zero): one line with their
  # count and the largest residual, which shows that they hold
  conRows = [row for row in base.measurement_ranking if get(row, :constraint, false)]
  if !isempty(conRows)
    worst = maximum(abs(row.residual) for row in conRows)
    println(io)
    @printf(io, "%s %d row(s) solved as equality constraints (not measurements), largest |residual| %.3g\n", format == :markdown ? "**Zero-injection constraints:**" : "Zero-injection constraints:", length(conRows), worst)
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
        @printf(io, "| %d | %s | %d | %.2f | %.3f |\n", rr.measurement_index, string("`", rr.id, "`"), rr.stage, rr.t, rr.sigma_factor)
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
        @printf(io, "| %d | %d | %s | %s | %.5f | %.6f | %.6f |\n", t.elimination, t.measurement_index, string("`", t.id, "`"), string(t.typ), t.normalized_residual_before, t.objective_before, t.objective_after)
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
  if hasproperty(diag, :reverted) && diag.reverted !== nothing
    r = diag.reverted
    if format == :markdown
      println(io, "\n- **Reverted:** elimination $(r.elimination) of `$(r.id)` ($(r.typ), normalized residual $(round(r.normalized_residual_before; digits = 5))) was not confirmed: the control solve without it did not converge, so the row stays active and the elimination stopped (not_converged).")
    else
      println(io, "\nReverted: elimination $(r.elimination) of $(r.id) ($(r.typ), rn $(round(r.normalized_residual_before; digits = 5))) was not confirmed: the control solve without it did not converge; the row stays active (stop: not_converged)")
    end
  end

  if hasproperty(diag, :rerun) && !isnothing(diag.rerun)
    rerun = diag.rerun
    rerun_summary = summarize_se_diagnostics(rerun.diagnostics)
    if format == :markdown
      println(io)
      println(io, "### Deactivate-and-rerun")
      println(io, "- **Deactivated measurement:** idx=$(rerun.deactivated_measurement_index), id=`$(rerun.deactivated_measurement_id)`")
      @printf(io, "- **Objective before:** %.6f\n", base.objective.value)
      @printf(io, "- **Objective after:** %.6f\n", rerun.diagnostics.objective.value)
      @printf(io, "- **Objective after stats:** dof=%d, z=%.3f, within_3sigma=%s\n", rerun_summary.objective.dof, rerun_summary.objective.zscore, string(rerun_summary.objective.within_3sigma))
      println(io, "- **Global consistency after rerun:** $(rerun_summary.global_consistency) ($(_global_consistency_label(rerun_summary)))")
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

  # released transformer taps: mode and the released regulators (tap
  # changers of the branch's changer set) per released trafo, for the same
  # at-a-glance transparency as the shunts
  taps = NamedTuple[]
  for (k, br) in enumerate(net.branchVec)
    br.tap_est_mode == :none && continue
    set = _effective_tap_changers(br)
    regs = NamedTuple[]
    if set !== nothing
      for ch in set.changers
        _tap_role_released(ch.role, br.tap_est_mode) || continue
        push!(regs, (changer = ch.name, kind = ch.kind, role = ch.role, side = ch.side, alpha_deg = ch.alpha_deg, du = ch.du, dphi_deg = ch.dphi_deg, neutral = ch.neutral, low = ch.low, high = ch.high))
      end
    end
    push!(taps, (branch = k, name = getCompName(br.comp), mode = br.tap_est_mode, regulators = regs, closed = _branch_terminal_state(br) == :closed))
  end

  return (controllers = controllers, shunts = shunts, taps = taps, links = (closed_clusters = closed_clusters, open_links = open_links), excluded_measurements = prep.notes)
end

# one regulator of an se_view release: kind with the letter of the changer
# type (a ratio, b asymmetrical, c quadrature, d symmetrical, e table,
# f linear), its step data and side
function _se_view_regulator_text(g)::String
  k = g.kind
  if k === :ratio
    return "ratio (a) du $(g.du) side $(g.side)"
  elseif k === :asymmetrical
    return "asymmetrical ($(abs(g.alpha_deg - 90.0) <= 1e-9 ? "c" : "b")) alpha $(g.alpha_deg) du $(g.du) side $(g.side)"
  elseif k === :symmetrical
    return "symmetrical (d) du $(g.du) side $(g.side)"
  elseif k === :linear
    return "linear (f) dphi $(g.dphi_deg) deg side $(g.side)"
  end
  return "table (e) steps $(g.low)..$(g.high) side $(g.side)"
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
    regs = hasproperty(t, :regulators) ? join((_se_view_regulator_text(g) for g in t.regulators), "; ") : ""
    line = string("branch $(t.branch) ($(t.name)): mode $(t.mode)", isempty(regs) ? "" : "; " * regs, t.closed ? "" : " (open, ignored)")
    println(io, md ? "  - " * line : "  " * line)
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

## `link_circulation` (0.60.0): the estimated ring circulation of every
## link that has one, linkIdx => (P MW, Q Mvar), as `calcLinkFlowsSE!`
## reads it; an empty Dict when the net has closed links but no estimated
## circulation, `nothing` without closed links or from readSEStateCSV!.
## `node_vm`/`node_va` are the node voltages of `net` at registration, read
## from the nodes themselves: `_se_link_circulation_current` hands out the
## circulation only while the net still holds exactly these voltages, so a
## later power flow (or any other state write) does not get the circulation
## of an old estimate.
function _register_se_start!(net, vm::Vector{Float64}, va::Vector{Float64}, pinj::Vector{Float64}, qinj::Vector{Float64}, converged::Bool; link_circulation::Union{Nothing,Dict{Int,Tuple{Float64,Float64}}} = nothing)
  node_vm = Union{Nothing,Float64}[n._vm_pu for n in net.nodeVec]
  node_va = Union{Nothing,Float64}[n._va_deg for n in net.nodeVec]
  _SE_START_STATES[net] = (vm = vm, va = va, pinj = pinj, qinj = qinj, converged = converged, link_circulation = link_circulation, node_vm = node_vm, node_va = node_va)
  return nothing
end

## The registered link circulations of `net` when its node voltages are still
## bitwise those the estimate wrote, else `nothing` (no estimate, no closed
## link, or the state changed since, e.g. by runpf!).
function _se_link_circulation_current(net)
  st = _se_start_state(net)
  (st === nothing || st.link_circulation === nothing) && return nothing
  length(st.node_vm) == length(net.nodeVec) || return nothing
  for (i, n) in enumerate(net.nodeVec)
    (isequal(n._vm_pu, st.node_vm[i]) && isequal(n._va_deg, st.node_va[i])) || return nothing
  end
  return st.link_circulation
end

## The circulation part of the link flows of `onet` from the circulation
## states `c` of its working net `wnet`: linkIdx => (P, Q). `nothing` when
## `onet` has no closed link.
function _se_link_circulation_flows(onet::Net, wnet::Net, c)::Union{Nothing,Dict{Int,Tuple{Float64,Float64}}}
  any(l -> l.status == 1, onet.linkVec) || return nothing
  out = Dict{Int,Tuple{Float64,Float64}}()
  (c === nothing || !(wnet._se_link_rows isa _SELinkTable)) && return out
  tbl = wnet._se_link_rows
  for s in 1:tbl.ncirc
    for (lk, a) in tbl.circ_links[s]
      p, q = get(out, lk, (0.0, 0.0))
      out[lk] = tbl.circ_typ[s] == PflowMeas ? (p + a * c[s], q) : (p, q + a * c[s])
    end
  end
  return out
end

@inline _se_start_state(net) = get(_SE_START_STATES, net, nothing)

const _SE_STATE_CSV_VERSION = "sparlectra-se-state v1"

"""
    writeSEStateCSV(net; file) -> NamedTuple

Write the last SE start state of `net` (registered by
`runse!` with `state_estimation.update_net = true` or `readSEStateCSV!`) as a CSV artifact:

    # sparlectra-se-state v1
    bus,vm_pu,va_deg,pinj_MW,qinj_MVar

(delimiter and decimal separator follow `output.csv_format`, `format`
keyword; the reader accepts every format)

One row per bus (bus by name). This is the persistence half of the SE -> PF
chain: a later power-flow run (possibly in another process) restores the
state with `readSEStateCSV!` and starts via `runpf_from_se!`. Errors when no
SE state is registered for `net`.
"""
function writeSEStateCSV(net::Net; file::AbstractString, format = result_csv_format())
  st = _se_start_state(net)
  st === nothing && error("writeSEStateCSV: no state-estimation result registered for this net; run runse! with state_estimation.update_net = true first")
  name_by_idx = _bus_name_by_idx(net)
  # a result artifact (issue #386): output.csv_format like every other CSV of
  # the run; readSEStateCSV! accepts both delimiters
  rows = ((get(name_by_idx, i, string(i)), st.vm[i], st.va[i], st.pinj[i], st.qinj[i]) for i in eachindex(net.nodeVec))
  write_result_csv(file, ("bus", "vm_pu", "va_deg", "pinj_MW", "qinj_MVar"), rows; format = format, comments = [_SE_STATE_CSV_VERSION])
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
  # the file follows output.csv_format of the run that wrote it: the header
  # line says which delimiter (and thereby which decimal separator) it uses
  delimiter = length(lines) >= 2 ? result_csv_delimiter(lines[2]) : ','
  for (ln, raw) in enumerate(lines)
    ln <= 2 && continue   # version + header
    line = strip(raw)
    isempty(line) && continue
    fields = split(line, delimiter)
    length(fields) == 5 || error("$(file):$(ln): expected 5 fields, got $(length(fields))")
    busIdx = try
      geNetBusIdx(net = net, busName = String(strip(fields[1])))
    catch
      error("$(file):$(ln): unknown bus '$(strip(fields[1]))'")
    end
    vals = map(f -> _parse_result_csv_number(f, delimiter), fields[2:5])
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
(the SE-to-PF chain). Requires a preceding `runse!` (`state_estimation.update_net =
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
  st === nothing && error("runpf_from_se!: no preceding state-estimation result for this net; run runse! with state_estimation.update_net = true or readSEStateCSV! first")
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
