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
# Rectangular power-flow wrong-branch diagnostics.

# Date: 29.5.2026
# file: src/powerflow_rectangular/rectangular_wrong_branch.jl
# purpose: wrong-branch diagnostics for rectangular power flow: angle
#          wrapping, circular angle-spread metrics, and plausibility checks
#          on converged solutions

"""
Normalize an angle (degrees) to the interval `(-180, 180]`.

# Why
Wrong-branch diagnostics compare relative angles from different buses and branches.
A canonical wrap prevents false spread inflation due to `±360°` representation jumps.
"""
@inline function _wrap_to_180_deg(angle_deg::Float64)::Float64
  wrapped = mod(angle_deg + 180.0, 360.0) - 180.0
  # Use +180 instead of -180 as canonical boundary to keep comparisons deterministic.
  return wrapped == -180.0 ? 180.0 : wrapped
end

"""
Compute circular angle spread (degrees) from a collection of angles.

Returns the minimal covering arc on the unit circle after wrap-to-180 normalization.
Non-finite entries are ignored.
"""
function _circular_angle_spread_deg(angles_deg)::Float64
  # Circular spread is measured after wrap-to-180 normalization.
  vals = Float64[_wrap_to_180_deg(Float64(a)) for a in angles_deg if isfinite(a)]
  n = length(vals)
  n == 0 && return NaN
  n == 1 && return 0.0

  vals360 = sort(mod.(vals, 360.0))
  # Largest gap complement gives the minimal circular spread.
  largest_gap = vals360[1] + 360.0 - vals360[end]
  for i = 1:(n-1)
    gap = vals360[i+1] - vals360[i]
    if gap > largest_gap
      largest_gap = gap
    end
  end

  spread = 360.0 - largest_gap
  return max(0.0, min(360.0, spread))
end

"""
Build a normalized wrong-branch diagnostic result payload.

# Why
A single NamedTuple shape keeps caller/reporting logic simple and stable, even when
some diagnostics are not available (e.g., no branch checks).
"""
@inline function _wrong_branch_result(;
  status::Symbol,
  reason::Symbol,
  min_vm_pu::Float64 = NaN,
  max_vm_pu::Float64 = NaN,
  low_vm_count::Int = 0,
  high_vm_count::Int = 0,
  angle_spread_deg::Float64 = NaN,
  max_branch_angle_deg::Float64 = NaN,
  worst_branch_angle_deg::Float64 = NaN,
  branch_angle_violation_count::Int = 0,
  worst_branch = nothing,
  lowest_buses::Vector{Int} = Int[],
)
  return (; status, reason, min_vm_pu, max_vm_pu, low_vm_count, high_vm_count, angle_spread_deg, max_branch_angle_deg, worst_branch_angle_deg, branch_angle_violation_count, worst_branch, lowest_buses)
end

# Explicit "diagnostics disabled" payload to distinguish from a checked-and-ok state.
@inline _wrong_branch_not_checked_result() = _wrong_branch_result(status = :not_checked, reason = :disabled)

"""
Convenience overload that accepts `net` as first positional argument.

Forwards to the keyword-based implementation to keep all decision logic in one place.
"""
function _check_wrong_branch_solution(net::Net, V::Vector{ComplexF64}, bus_types::Vector{Symbol}, Vset::Vector{Float64}, slack_idx::Int; min_vm_pu::Float64, max_vm_pu::Float64, max_angle_spread_deg::Float64, max_branch_angle_deg::Float64 = Inf, min_low_vm_count::Int)
  return _check_wrong_branch_solution(V, bus_types, Vset, slack_idx; net = net, min_vm_pu = min_vm_pu, max_vm_pu = max_vm_pu, max_angle_spread_deg = max_angle_spread_deg, max_branch_angle_deg = max_branch_angle_deg, min_low_vm_count = min_low_vm_count)
end

"""
Evaluate a solved voltage state for wrong-branch indicators.

The check is intentionally conservative:
- `:fail` is reserved for invalid states (e.g. non-finite voltages),
- suspicious but finite states are returned as `:warn`,
- final policy (`:warn`, `:fail`, potential rescue path) is decided by caller logic.

When `net` is available, all heuristics (voltage band, angle spread,
branch-angle) are evaluated only on the network's highest nominal voltage
level; lower levels carry normal operating spread that would produce false
SUSPECT verdicts. The reported `min_vm_pu`/`max_vm_pu`/`lowest_buses` refer
to that level. Without a `net` every bus is in scope.

# Returns
A normalized diagnostic NamedTuple created by `_wrong_branch_result`.
"""
function _check_wrong_branch_solution(V::Vector{ComplexF64}, bus_types::Vector{Symbol}, Vset::Vector{Float64}, slack_idx::Int; net::Union{Nothing,Net} = nothing, min_vm_pu::Float64, max_vm_pu::Float64, max_angle_spread_deg::Float64, max_branch_angle_deg::Float64 = Inf, min_low_vm_count::Int)
  # This helper is intentionally conservative: suspicious states are classified
  # as :warn so caller policy (:warn/:fail/:rescue) decides final acceptance.
  n = length(V)
  n == length(bus_types) || throw(ArgumentError("V and bus_types must have same length."))
  n == length(Vset) || throw(ArgumentError("V and Vset must have same length."))
  1 <= slack_idx <= n || throw(ArgumentError("slack_idx must be inside V."))

  if any(v -> !isfinite(real(v)) || !isfinite(imag(v)), V)
    # Non-finite voltages are not recoverable for mismatch/Jacobian logic.
    return _wrong_branch_result(status = :fail, reason = :nonfinite_voltage)
  end

  # Wrong-branch heuristics judge only the network's HIGHEST voltage level
  # (maintainer decision 2026-07-30): sub-transmission dips are normal
  # operating spread, not a wrong-branch indicator — a real 6209-bus CGMES
  # snapshot flagged healthy 45 kV feeders at 0.94 pu as SUSPECT while its
  # 380 kV level was clean. Without a net (pure vector call) or without
  # nominal-voltage data every bus stays in scope, preserving the previous
  # behavior for the unit-level callers.
  level_mask = trues(n)
  if !isnothing(net) && length(net.nodeVec) == n
    iso = Set(net.isoNodes)
    vn_max = 0.0
    for i in eachindex(net.nodeVec)
      i in iso && continue
      vn = getNodeVn(net.nodeVec[i])
      vn > vn_max && (vn_max = vn)
    end
    if vn_max > 0.0
      for i in eachindex(net.nodeVec)
        level_mask[i] = !(i in iso) && isapprox(getNodeVn(net.nodeVec[i]), vn_max; rtol = 1e-6)
      end
    end
  end
  level_indices = findall(level_mask)
  isempty(level_indices) && (level_mask = trues(n); level_indices = collect(1:n))

  vm = abs.(V)
  low_idx = [i for i in level_indices if vm[i] < min_vm_pu]
  high_idx = [i for i in level_indices if vm[i] > max_vm_pu]
  va_deg = rad2deg.(angle.(V))
  slack_ang = va_deg[slack_idx]
  # Relative-angle diagnostics are measured against slack as reference; the
  # spread is reference-invariant, so a slack below the top level is fine.
  rel = _wrap_to_180_deg.(va_deg .- slack_ang)
  angle_spread_deg = _circular_angle_spread_deg(rel[level_indices])

  branch_angle_violation_count = 0
  max_branch_angle_seen_deg = NaN
  worst_branch = nothing
  if !isnothing(net) && isfinite(max_branch_angle_deg)
    # Non-finite branch-angle diagnostics are treated as suspicious and count
    # against wrong-branch detection by escalating to warning/fail conditions.
    max_branch_angle_seen_deg = 0.0
    # Branch-angle diagnostics intentionally scan net.branchVec only.
    # Non-standard couplers/impedanceless links live in net.linkVec and are excluded.
    @inbounds for br in net.branchVec
      br.status == 1 || continue
      from = br.fromBus
      to = br.toBus
      (1 <= from <= n && 1 <= to <= n) || continue
      # Highest-level scope (see level_mask above): only branches with both
      # ends on the top voltage level count for the branch-angle heuristic.
      (level_mask[from] && level_mask[to]) || continue

      theta_from_deg = rad2deg(angle(V[from]))
      theta_to_deg = rad2deg(angle(V[to]))
      raw_diff_deg = _wrap_to_180_deg(theta_from_deg - theta_to_deg)

      # Effective branch-angle check compensates transformer phase shift:
      # bus-angle difference minus branch phase shift is the physical stress basis.
      phase_shift_deg = hasproperty(br, :phase_shift_deg) ? Float64(getproperty(br, :phase_shift_deg)) : 0.0
      eff_diff_deg = _wrap_to_180_deg(raw_diff_deg - phase_shift_deg)
      eff_diff_abs_deg = abs(eff_diff_deg)

      if eff_diff_abs_deg > max_branch_angle_seen_deg
        max_branch_angle_seen_deg = eff_diff_abs_deg
        worst_branch = (branch_index = br.branchIdx, from_bus = from, to_bus = to, active = true, phase_shift_deg = phase_shift_deg, angle_diff_raw_deg = abs(raw_diff_deg), angle_diff_effective_deg = eff_diff_abs_deg, angle_check_basis = :effective_bus_angle_minus_phase_shift)
      end
      if eff_diff_abs_deg > max_branch_angle_deg
        branch_angle_violation_count += 1
      end
    end
  end

  # Keep the three lowest-voltage buses (within the checked level) for
  # compact diagnostics in logs/tables.
  lowest_order = sort(level_indices; by = i -> vm[i])
  lowest_buses = [Int(i) for i in lowest_order[1:min(3, length(lowest_order))]]

  # :warn => suspicious but usable; :fail => invalid (e.g., non-finite voltages);
  # :not_checked is emitted by callers when diagnostics are disabled.
  status = :ok
  reason = :none
  if !isempty(low_idx) && length(low_idx) >= min_low_vm_count
    status = :warn
    reason = :low_voltage_magnitude
  elseif !isempty(high_idx)
    status = :warn
    reason = :high_voltage_magnitude
  elseif angle_spread_deg > max_angle_spread_deg
    status = :warn
    reason = :angle_spread_exceeded
  elseif branch_angle_violation_count > 0
    status = :warn
    reason = :branch_angle_exceeded
  end

  return _wrong_branch_result(
    status = status,
    reason = reason,
    min_vm_pu = minimum(vm[level_indices]),
    max_vm_pu = maximum(vm[level_indices]),
    low_vm_count = length(low_idx),
    high_vm_count = length(high_idx),
    angle_spread_deg = angle_spread_deg,
    max_branch_angle_deg = max_branch_angle_deg,
    worst_branch_angle_deg = max_branch_angle_seen_deg,
    branch_angle_violation_count = branch_angle_violation_count,
    worst_branch = worst_branch,
    lowest_buses = lowest_buses,
  )
end
