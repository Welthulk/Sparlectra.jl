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
# file: src/jacobian_complex.jl
#

# jacobian_complex.jl — Complex-State Newton-Raphson Power Flow Formulation
#
# This module implements a Newton-Raphson power flow solver using complex voltages
# in rectangular coordinates (Vr + jVi) as state variables, as an alternative to
# the conventional polar formulation (Vm, θ).
#
# Features:
# - Rectangular complex-state Newton-Raphson with PQ and PV bus handling
# - Wirtinger calculus-based Jacobian construction for complex power equations
# - Active-set Q-limit management for PV→PQ switching with optional re-enable
# - Both analytic and finite-difference Jacobian options
# - Hysteresis and cooldown mechanisms for robust PV bus management
# - Direct integration with Sparlectra.jl network data structures
#
# Mathematical Foundation:
# - State vector: x = [Vr(non-slack); Vi(non-slack)] ∈ ℝ^(2(n-1))
# - Complex power: S = V .* conj(Y * V) where Y is the bus admittance matrix  
# - PQ buses: ΔP = Re(S_calc - S_spec), ΔQ = Im(S_calc - S_spec)
# - PV buses: ΔP = Re(S_calc - S_spec), ΔV = |V| - V_set
# - Slack bus voltage is held constant throughout the iteration
#
# Key Functions:
# - run_complex_nr_rectangular_for_net!(): Main solver interface
# - runpf_rectangular!(): Convenience wrapper matching runpf!() signature
# - build_complex_jacobian(): Wirtinger-based Jacobian block construction
# - mismatch_rectangular(): Residual function for PQ/PV bus constraints
#
# Note:
# - The FD Jacobian is mathematically dense; sparse storage does not bring much benefit.
# - The analytic Jacobian currently uses a dense rectangular build, even though the
#   underlying structure is sparse (Ybus-like). A true sparse implementation would
#   require a dedicated builder similar to `calcJacobian(...; sparse=true)`.
#
# References:
# - Wirtinger calculus for complex derivatives

using LinearAlgebra
using SparseArrays
using Printf

"""
    build_complex_jacobian(Ybus, V)

Builds the 2n × 2n Wirtinger-type Jacobian blocks for the complex-state
Newton–Raphson formulation.

Given:
    I = Ybus * V
    S = V .* conj.(I)

We construct the blocks:
    J11 = ∂S/∂V
    J12 = ∂S/∂V*
    J21 = ∂conj(S)/∂V
    J22 = ∂conj(S)/∂V*

Returns:
    J11, J12, J21, J22  (all full matrices, not Diagonal)
"""
function build_complex_jacobian(Ybus, V)
  I = Ybus * V
  n = length(V)

  # J11 = diag(conj(I))
  J11 = Matrix(Diagonal(conj.(I)))

  # J12 = diag(V) * conj(Ybus)
  J12 = Matrix(Diagonal(V)) * conj.(Ybus)

  # J21 = conj(J12)
  J21 = conj.(J12)

  # J22 = conj(J11)
  J22 = conj.(J11)

  return J11, J12, J21, J22
end

"""
    mismatch_rectangular(Ybus, V, S, bus_types, Vset, slack_idx) -> F::Vector{Float64}

Compute the real-valued mismatch vector F(V) for the rectangular
complex-state formulation with PQ and PV buses.

For each non-slack bus i:
- if bus_types[i] == :PQ:
      ΔP_i = Re(S_calc[i]) - Re(S_spec[i])
      ΔQ_i = Im(S_calc[i]) - Im(S_spec[i])

- if bus_types[i] == :PV:
      ΔP_i = Re(S_calc[i]) - Re(S_spec[i])
      ΔV_i = |V[i]| - Vset[i]

F is stacked as [ΔP_2, ΔQ/ΔV_2, ..., ΔP_n, ΔQ/ΔV_n] over all non-slack buses.
"""
function mismatch_rectangular(Ybus, V::Vector{ComplexF64}, S::Vector{ComplexF64}, bus_types::Vector{Symbol}, Vset::Vector{Float64}, slack_idx::Int)
  n = length(V)
  @assert length(S) == n
  @assert length(bus_types) == n
  @assert length(Vset) == n
  # Network-based injections for the current state
  S_calc = calc_injections(Ybus, V)

  # F has 2*(n-1) entries: for each non-slack bus two residuals
  # PQ:  ΔP_i, ΔQ_i
  # PV:  ΔP_i, ΔV_i
  F = zeros(Float64, 2 * (n - 1))

  row = 1
  @inbounds for i = 1:n
    if i == slack_idx
      continue
    end

    S_ci = S_calc[i]
    S_si = S[i]

    if bus_types[i] == :PQ
      ΔP = real(S_ci) - real(S_si)
      ΔQ = imag(S_ci) - imag(S_si)

      F[row]   = ΔP
      F[row+1] = ΔQ

    elseif bus_types[i] == :PV
      ΔP = real(S_ci) - real(S_si)
      ΔV = abs(V[i]) - Vset[i]

      F[row]   = ΔP
      F[row+1] = ΔV

    else
      error("mismatch_rectangular: unsupported bus type $(bus_types[i]) at bus $i")
    end

    row += 2
  end

  return F
end

function run_complex_nr_rectangular(
  Ybus,
  V0,
  S;
  slack_idx::Int = 1,
  maxiter::Int = 20,
  tol::Float64 = 1e-8,
  verbose::Bool = false,
  damp::Float64 = 1.0,
  autodamp::Bool = false,
  autodamp_min::Float64 = 1e-3,
  bus_types::Vector{Symbol},
  Vset::Vector{Float64},
  use_fd::Bool = false,
  use_sparse::Bool = false,
  dPinj_dVm::Vector{Float64} = zeros(Float64, length(V0)),
  dQinj_dVm::Vector{Float64} = zeros(Float64, length(V0)),
)
  V = copy(V0)
  history = Float64[]

  for iter = 1:maxiter
    F = mismatch_rectangular(Ybus, V, S, bus_types, Vset, slack_idx)
    max_mis = maximum(abs.(F))
    push!(history, max_mis)

    if max_mis <= tol
      return V, true, iter, history
    end

    if use_fd
      V = complex_newton_step_rectangular_fd(Ybus, V, S; slack_idx = slack_idx, damp = damp, autodamp = autodamp, autodamp_min = autodamp_min, h = 1e-6, bus_types = bus_types, Vset = Vset, dPinj_dVm = dPinj_dVm, dQinj_dVm = dQinj_dVm)
    else
      V = complex_newton_step_rectangular(Ybus, V, S; slack_idx = slack_idx, damp = damp, autodamp = autodamp, autodamp_min = autodamp_min, bus_types = bus_types, Vset = Vset, use_sparse = use_sparse, dPinj_dVm = dPinj_dVm, dQinj_dVm = dQinj_dVm)
    end
  end

  return V, false, maxiter, history
end

function _sanitize_rectangular_start(V::Vector{ComplexF64}, bus_types::Vector{Symbol}, Vset::Vector{Float64}, slack_idx::Int)
  Vs = copy(V)
  @inbounds for k in eachindex(Vs)
    Vk = Vs[k]
    vm = abs(Vk)
    if !isfinite(real(Vk)) || !isfinite(imag(Vk)) || vm <= 0.0
      vm = (bus_types[k] in (:Slack, :PV) && isfinite(Vset[k]) && Vset[k] > 0.0) ? Vset[k] : 1.0
      Vs[k] = ComplexF64(vm, 0.0)
    end
  end
  Vs[slack_idx] = V[slack_idx]
  return Vs
end

function _voltage_magnitude_for_projection(Vraw::Vector{ComplexF64}, bus_types::Vector{Symbol}, Vset::Vector{Float64}, k::Int)
  vm_raw = abs(Vraw[k])
  if bus_types[k] in (:Slack, :PV) && isfinite(Vset[k]) && Vset[k] > 0.0
    return Vset[k]
  elseif isfinite(vm_raw) && vm_raw > 0.0
    return vm_raw
  else
    return 1.0
  end
end

function _dc_angle_start_rectangular(
  Ybus,
  Vraw::Vector{ComplexF64},
  S::Vector{ComplexF64},
  bus_types::Vector{Symbol},
  Vset::Vector{Float64},
  slack_idx::Int;
  dc_angle_limit_deg::Float64 = 60.0,
)
  n = length(Vraw)
  non_slack = non_slack_indices(n, slack_idx)
  nred = length(non_slack)
  nred == 0 && return copy(Vraw)
  pos = build_pos_map(non_slack, n)

  B = zeros(Float64, nred, nred)
  if Ybus isa SparseMatrixCSC
    rv = rowvals(Ybus)
    nz = nonzeros(Ybus)
    @inbounds for j in 1:n
      for ptr in nzrange(Ybus, j)
        i = rv[ptr]
        i == j && continue
        ri = pos[i]
        ri == 0 && continue
        bij = imag(nz[ptr])
        bij == 0.0 && continue
        if j != slack_idx
          cj = pos[j]
          if cj != 0
            B[ri, cj] -= bij
          end
        end
        B[ri, ri] += bij
      end
    end
  else
    @inbounds for i in non_slack
      ri = pos[i]
      for j in 1:n
        j == i && continue
        bij = imag(Ybus[i, j])
        bij == 0.0 && continue
        if j != slack_idx
          cj = pos[j]
          if cj != 0
            B[ri, cj] -= bij
          end
        end
        B[ri, ri] += bij
      end
    end
  end

  P = real.(S[non_slack])
  θred = solve_linear(B, P; allow_pinv = true)
  limit = deg2rad(dc_angle_limit_deg)
  θslack = angle(Vraw[slack_idx])
  Vdc = similar(Vraw)

  @inbounds for k in 1:n
    vm = _voltage_magnitude_for_projection(Vraw, bus_types, Vset, k)
    if k == slack_idx
      Vdc[k] = Vraw[k]
    else
      θ = θslack + clamp(θred[pos[k]], -limit, limit)
      Vdc[k] = ComplexF64(vm * cos(θ), vm * sin(θ))
    end
  end
  return Vdc
end

function _blend_voltage_starts(Vraw::Vector{ComplexF64}, Vdc::Vector{ComplexF64}, λ::Float64, slack_idx::Int)
  0.0 <= λ <= 1.0 || error("blend lambda must satisfy 0 ≤ λ ≤ 1 (got $(λ)).")
  V = similar(Vraw)
  @inbounds for k in eachindex(Vraw)
    if k == slack_idx
      V[k] = Vraw[k]
      continue
    end
    vm = (1.0 - λ) * abs(Vraw[k]) + λ * abs(Vdc[k])
    θ = (1.0 - λ) * angle(Vraw[k]) + λ * angle(Vdc[k])
    V[k] = ComplexF64(vm * cos(θ), vm * sin(θ))
  end
  return V
end

"""
    project_rectangular_start(Ybus, Vraw, S, bus_types, Vset, slack_idx; ...)

Build a projected initial voltage for the rectangular power-flow solver.

When enabled, the projection sanitizes the raw seed, optionally computes a
DC-angle start from active-power injections and the Y-bus off-diagonal
susceptances, and optionally scans convex blends between the raw and DC starts.
The candidate with the lowest rectangular mismatch is returned.
"""
function project_rectangular_start(
  Ybus,
  Vraw::Vector{ComplexF64},
  S::Vector{ComplexF64},
  bus_types::Vector{Symbol},
  Vset::Vector{Float64},
  slack_idx::Int;
  enabled::Bool = false,
  try_dc_start::Bool = true,
  try_blend_scan::Bool = true,
  blend_lambdas::AbstractVector{<:Real} = [0.25, 0.5, 0.75],
  dc_angle_limit_deg::Float64 = 60.0,
  verbose::Int = 0,
)
  enabled || return Vraw
  dc_angle_limit_deg > 0.0 || error("dc_angle_limit_deg must be > 0 (got $(dc_angle_limit_deg)).")

  raw = _sanitize_rectangular_start(Vraw, bus_types, Vset, slack_idx)
  best = raw
  best_name = :raw
  best_mis = _max_rectangular_mismatch(Ybus, raw, S, bus_types, Vset, slack_idx)
  Vdc = nothing

  if try_dc_start
    Vdc = _dc_angle_start_rectangular(Ybus, raw, S, bus_types, Vset, slack_idx; dc_angle_limit_deg = dc_angle_limit_deg)
    dc_mis = _max_rectangular_mismatch(Ybus, Vdc, S, bus_types, Vset, slack_idx)
    if isfinite(dc_mis) && dc_mis < best_mis
      best = Vdc
      best_name = :dc_start
      best_mis = dc_mis
    end
  end

  if try_blend_scan && Vdc !== nothing
    for λ_raw in blend_lambdas
      λ = Float64(λ_raw)
      Vblend = _blend_voltage_starts(raw, Vdc, λ, slack_idx)
      blend_mis = _max_rectangular_mismatch(Ybus, Vblend, S, bus_types, Vset, slack_idx)
      if isfinite(blend_mis) && blend_mis < best_mis
        best = Vblend
        best_name = Symbol("blend_", λ)
        best_mis = blend_mis
      end
    end
  end

  if verbose > 0
    raw_mis = _max_rectangular_mismatch(Ybus, raw, S, bus_types, Vset, slack_idx)
    @info "start projection selected $(best_name)" raw_mismatch = raw_mis projected_mismatch = best_mis
  end
  return best
end

function _validate_rectangular_damping(damp::Float64, autodamp_min::Float64)
  isfinite(damp) || error("damp must be finite (got $(damp)).")
  isfinite(autodamp_min) || error("autodamp_min must be finite (got $(autodamp_min)).")
  0.0 < damp <= 1.0 || error("damp must satisfy 0 < damp ≤ 1 (got $(damp)).")
  0.0 < autodamp_min <= damp || error("autodamp_min must satisfy 0 < autodamp_min ≤ damp (got autodamp_min=$(autodamp_min), damp=$(damp)).")
  return nothing
end

function _apply_rectangular_delta(V::Vector{ComplexF64}, δx::Vector{Float64}, slack_idx::Int, non_slack::Vector{Int}, alpha::Float64)
  n = length(V)
  Vr = real.(V)
  Vi = imag.(V)
  Vr_new = copy(Vr)
  Vi_new = copy(Vi)

  @inbounds for (idx, bus) in enumerate(non_slack)
    Vr_new[bus] += alpha * δx[idx]
    Vi_new[bus] += alpha * δx[(n - 1) + idx]
  end

  Vr_new[slack_idx] = Vr[slack_idx]
  Vi_new[slack_idx] = Vi[slack_idx]
  return ComplexF64.(Vr_new, Vi_new)
end

function _max_rectangular_mismatch(Ybus, V::Vector{ComplexF64}, S::Vector{ComplexF64}, bus_types::Vector{Symbol}, Vset::Vector{Float64}, slack_idx::Int)
  F = mismatch_rectangular(Ybus, V, S, bus_types, Vset, slack_idx)
  return maximum(abs.(F))
end

"""
    choose_rectangular_autodamp(Ybus, V, S, δx, F0; slack_idx, damp, autodamp_min, bus_types, Vset)

Select a Newton step length for the rectangular power-flow solver by backtracking
from `damp` toward `autodamp_min`. The first trial step that reduces the maximum
absolute mismatch is accepted. If no trial reduces the mismatch, the smallest
finite trial is returned so the solver can continue safely with a conservative
step.

Returns `(alpha, Vtrial, trial_mismatch)`.
"""
function choose_rectangular_autodamp(
  Ybus,
  V::Vector{ComplexF64},
  S::Vector{ComplexF64},
  δx::Vector{Float64},
  F0::Vector{Float64};
  slack_idx::Int,
  damp::Float64 = 1.0,
  autodamp_min::Float64 = 1e-3,
  bus_types::Vector{Symbol},
  Vset::Vector{Float64},
)
  _validate_rectangular_damping(damp, autodamp_min)
  non_slack = non_slack_indices(length(V), slack_idx)
  current_mismatch = maximum(abs.(F0))
  best_alpha = autodamp_min
  best_V = _apply_rectangular_delta(V, δx, slack_idx, non_slack, autodamp_min)
  best_mismatch = _max_rectangular_mismatch(Ybus, best_V, S, bus_types, Vset, slack_idx)

  alpha = damp
  while alpha >= autodamp_min
    Vtrial = _apply_rectangular_delta(V, δx, slack_idx, non_slack, alpha)
    trial_mismatch = _max_rectangular_mismatch(Ybus, Vtrial, S, bus_types, Vset, slack_idx)
    if isfinite(trial_mismatch) && trial_mismatch < current_mismatch
      return alpha, Vtrial, trial_mismatch
    end
    if isfinite(trial_mismatch) && trial_mismatch < best_mismatch
      best_alpha = alpha
      best_V = Vtrial
      best_mismatch = trial_mismatch
    end
    alpha *= 0.5
  end

  return best_alpha, best_V, best_mismatch
end

"""
    update_net_voltages_from_complex!(net, V)

Update the bus voltage magnitudes and angles in the network from the
final complex voltages V (in per-unit).
"""
function update_net_voltages_from_complex!(net::Net, V::Vector{ComplexF64})
  nodes = net.nodeVec
  n = length(nodes)
  @assert length(V) == n

  for (k, node) in enumerate(nodes)
    Vk = V[k]
    vm = abs(Vk)
    va_rad = angle(Vk)
    va_deg = rad2deg(va_rad)
    node._vm_pu = vm
    node._va_deg = va_deg
  end
end

function _expand_ybus_for_isolated_nodes(Yred, n::Int, iso_nodes::Vector{Int})
  isempty(iso_nodes) && return Yred

  iso_mask = falses(n)
  for bus in iso_nodes
    if 1 <= bus <= n
      iso_mask[bus] = true
    end
  end

  active = Int[]
  for bus in eachindex(iso_mask)
    iso_mask[bus] || push!(active, bus)
  end

  size(Yred, 1) == length(active) || error("_expand_ybus_for_isolated_nodes: size mismatch between reduced Ybus and active buses.")
  size(Yred, 2) == length(active) || error("_expand_ybus_for_isolated_nodes: Ybus is not square in active-bus space.")

  Yfull = issparse(Yred) ? spzeros(ComplexF64, n, n) : zeros(ComplexF64, n, n)
  Yfull[active, active] = Yred
  return Yfull
end

@inline function _has_vset_adjust_config(ps::ProSumer)::Bool
  return !isnothing(ps.vset_adjust) || !(isnothing(ps.vstep_pu) && isnothing(ps.tap_steps_down) && isnothing(ps.tap_steps_up))
end

@inline function _bus_label(net::Net, bus::Int)::String
  return getCompName(net.nodeVec[bus].comp)
end

@inline function _resolve_vset_adjust_config(ps::ProSumer)::Union{Nothing,VoltageAdjustConfig}
  if !isnothing(ps.vset_adjust)
    return ps.vset_adjust
  end
  has_any = _has_vset_adjust_config(ps)
  has_any || return nothing
  all_defined = !isnothing(ps.vstep_pu) && !isnothing(ps.tap_steps_down) && !isnothing(ps.tap_steps_up)
  all_defined || return nothing
  return VoltageAdjustConfig(Float64(ps.vstep_pu), Int(ps.tap_steps_down), Int(ps.tap_steps_up))
end

function _build_vset_adjust_controllers(net::Net)
  controllers = Dict{Int,NamedTuple{(:prosumer_idx, :config),Tuple{Int,VoltageAdjustConfig}}}()

  for (ps_idx, ps) in enumerate(net.prosumpsVec)
    isGenerator(ps) || continue
    bus = getPosumerBusIndex(ps)

    has_any = _has_vset_adjust_config(ps)
    if has_any
      cfg = _resolve_vset_adjust_config(ps)
      all_defined = !isnothing(cfg)
      all_defined || error("Bus $(_bus_label(net, bus)): invalid voltage adjustment config at prosumer $ps_idx. vstep_pu, tap_steps_down and tap_steps_up must be provided together.")
      cfg.vstep_pu > 0.0 || error("Bus $(_bus_label(net, bus)): invalid vstep_pu=$(cfg.vstep_pu). Must be > 0.")
      cfg.tap_steps_down >= 0 || error("Bus $(_bus_label(net, bus)): invalid tap_steps_down=$(cfg.tap_steps_down). Must be ≥ 0.")
      cfg.tap_steps_up >= 0 || error("Bus $(_bus_label(net, bus)): invalid tap_steps_up=$(cfg.tap_steps_up). Must be ≥ 0.")
      haskey(controllers, bus) && error("Bus $(_bus_label(net, bus)): multiple prosumers define voltage adjustment data. Only one controller per bus is allowed.")

      controllers[bus] = (prosumer_idx = ps_idx, config = cfg)
    else
      partially_set = !isnothing(ps.vstep_pu) || !isnothing(ps.tap_steps_down) || !isnothing(ps.tap_steps_up)
      partially_set && error("Bus $(_bus_label(net, bus)): incomplete voltage adjustment config at prosumer $ps_idx. vstep_pu, tap_steps_down and tap_steps_up must be provided together.")
    end
  end

  return controllers
end

function _try_adjust_vset_on_q_limit!(
  net::Net,
  bus::Int,
  side::Symbol,
  it::Int,
  controllers::Dict{Int,NamedTuple{(:prosumer_idx, :config),Tuple{Int,VoltageAdjustConfig}}},
  base_vset::Vector{Float64},
  Vset::Vector{Float64},
  adjust_counter::Vector{Int},
  qlimit_max_outer::Int,
  verbose::Int,
)::Bool
  cname = _bus_label(net, bus)
  if !haskey(controllers, bus)
    if verbose > 0
      @info "Bus $cname: no voltage adjustment controller -> fallback PV→PQ (it=$it)"
    end
    return false
  end

  ctrl = controllers[bus].config
  vm_base = base_vset[bus]
  vm_min = vm_base - ctrl.tap_steps_down * ctrl.vstep_pu
  vm_max = vm_base + ctrl.tap_steps_up * ctrl.vstep_pu
  vm_old = Vset[bus]
  vm_new = side == :max ? vm_old - ctrl.vstep_pu : vm_old + ctrl.vstep_pu
  can_step = (vm_new >= vm_min - 1e-12) && (vm_new <= vm_max + 1e-12)

  if can_step && adjust_counter[bus] < qlimit_max_outer
    vm_new = clamp(vm_new, vm_min, vm_max)
    Vset[bus] = vm_new
    net.nodeVec[bus]._vm_pu = vm_new
    adjust_counter[bus] += 1
    if verbose > 0
      event = side == :max ? "Qmax violated" : "Qmin violated"
      @info "Bus $cname: $event -> voltage adjusted from $vm_old to $vm_new (it=$it, step=$(adjust_counter[bus]))"
    end
    return true
  end

  if verbose > 0
    @info "Bus $cname: no further voltage steps (vm_min=$vm_min, vm_max=$vm_max, vm=$vm_old) -> fallback PV→PQ (it=$it)"
  end
  return false
end

"""
    build_rectangular_jacobian_pq_pv_sparse(
        Ybus,
        V,
        bus_types,
        Vset,
        slack_idx,
    ) -> SparseMatrixCSC{Float64}

Builds the analytic rectangular Jacobian corresponding to `mismatch_rectangular`
using the sparsity pattern of `Ybus`.

State vector:
    x = [Vr(non-slack); Vi(non-slack)] ∈ ℝ^(2(n-1))

Residual F(V):
    - PQ buses: ΔP_i, ΔQ_i
    - PV buses: ΔP_i, ΔV_i
    - Slack bus: no equations

Jacobian entries are derived from
    S_i(V) = V_i * conj( (Ybus * V)_i )

Wirtinger-based identities:
    ∂S/∂V   = diag(conj(I)) + diag(V) * conj(Ybus)
    ∂S/∂V*  = diag(V) * conj(Ybus)

Chain rule to rectangular:
    ∂S/∂Vr = ∂S/∂V + ∂S/∂V*
    ∂S/∂Vi = j(∂S/∂V - ∂S/∂V*)

With ΔP_i = Re(ΔS_i), ΔQ_i = Im(ΔS_i), ΔV_i = |V_i| - Vset[i].

Returns:
    J :: SparseMatrixCSC{Float64} with size (2(n-1)) × (2(n-1)).
"""
function build_rectangular_jacobian_pq_pv_sparse(
  Ybus::SparseMatrixCSC{ComplexF64},
  V::Vector{ComplexF64},
  bus_types::Vector{Symbol},
  Vset::Vector{Float64},
  slack_idx::Int;
  dPinj_dVm::Vector{Float64} = zeros(Float64, length(V)),
  dQinj_dVm::Vector{Float64} = zeros(Float64, length(V)),
  vm_eps::Float64 = 1e-9,
)
  n = length(V)
  @assert length(bus_types) == n
  @assert length(Vset) == n
  @assert length(dPinj_dVm) == n
  @assert length(dQinj_dVm) == n

  I = Ybus * V

  non_slack = non_slack_indices(n, slack_idx)
  pos_non_slack = build_pos_map(non_slack, n)
  nvar = 2 * (n - 1)   # [Vr(non-slack); Vi(non-slack)]
  m = nvar          # F has 2 equations per non-slack bus

  # Row blocks: for each non-slack bus i
  #   row_block[i]   = index of ΔP row for bus i
  #   row_block[i]+1 = index of ΔQ / ΔV row for bus i
  row_block = zeros(Int, n)
  row = 0
  for i = 1:n
    if i == slack_idx
      continue
    end
    row += 2
    row_block[i] = row - 1  # ΔP row
  end

  # Triplet storage
  Iidx = Int[]
  Jidx = Int[]
  Vals = Float64[]
  # Rough capacity hint (purely heuristic)
  sizehint!(Iidx, 16 * nnz(Ybus))
  sizehint!(Jidx, 16 * nnz(Ybus))
  sizehint!(Vals, 16 * nnz(Ybus))

  # Helpers for sparse access
  rv    = rowvals(Ybus)
  nzval = nonzeros(Ybus)

  # --- 1) Contributions from complex power equations (ΔP, ΔQ) ---------------
  #
  # For bus i, bus j:
  #   ∂S_i/∂Vr_j = conj(I_i) * δ_ij + V_i * conj(Y_ij)
  #   ∂S_i/∂Vi_j = j*(conj(I_i) * δ_ij - V_i * conj(Y_ij))
  #
  # Then:
  #   ∂P_i/∂Vr_j = Re(∂S_i/∂Vr_j)
  #   ∂P_i/∂Vi_j = Re(∂S_i/∂Vi_j)
  #   ∂Q_i/∂Vr_j = Im(∂S_i/∂Vr_j)
  #   ∂Q_i/∂Vi_j = Im(∂S_i/∂Vi_j)
  #
  # For PV buses, only the ΔP row uses these derivatives; the second row is ΔV.

  for j = 1:n
    col_pos = pos_non_slack[j]
    if col_pos == 0
      # Slack bus column -> no state variable
      continue
    end

    colVr = col_pos
    colVi = (n - 1) + col_pos

    for ptr in nzrange(Ybus, j)
      i = rv[ptr]

      # Slack bus has no equations
      if i == slack_idx
        continue
      end

      rb = row_block[i]
      if rb == 0
        continue
      end

      rowP = rb          # ΔP row for bus i
      rowQ = rb + 1      # ΔQ/ΔV row for bus i

      Yij = nzval[ptr]

      # Base contributions from V_i * conj(Y_ij)
      S_dVr = V[i] * conj(Yij)
      S_dVi = im * (-V[i] * conj(Yij))

      # Diagonal term from conj(I_i) * δ_ij
      if i == j
        Ii = I[i]
        S_dVr += conj(Ii)
        S_dVi += im * conj(Ii)
      end

      # Real / imaginary parts
      ∂P_Vr = real(S_dVr)
      ∂P_Vi = real(S_dVi)
      ∂Q_Vr = imag(S_dVr)
      ∂Q_Vi = imag(S_dVi)

      # First equation: always ΔP_i for PQ and PV
      if abs(∂P_Vr) > 0.0
        push!(Iidx, rowP)
        push!(Jidx, colVr)
        push!(Vals, ∂P_Vr)
      end
      if abs(∂P_Vi) > 0.0
        push!(Iidx, rowP)
        push!(Jidx, colVi)
        push!(Vals, ∂P_Vi)
      end

      # Second equation:
      #   - PQ: ΔQ_i -> uses Q derivatives
      #   - PV: ΔV_i -> no contribution from S, handled separately
      bt = bus_types[i]
      if bt == :PQ
        if abs(∂Q_Vr) > 0.0
          push!(Iidx, rowQ)
          push!(Jidx, colVr)
          push!(Vals, ∂Q_Vr)
        end
        if abs(∂Q_Vi) > 0.0
          push!(Iidx, rowQ)
          push!(Jidx, colVi)
          push!(Vals, ∂Q_Vi)
        end
      elseif bt == :PV
        # nothing here, ΔV row added below
      else
        error("build_rectangular_jacobian_pq_pv_sparse: unsupported bus type $(bt) at bus $i")
      end
    end
  end

  # --- 2) Contributions for ΔV_i = |V_i| - Vset[i] on PV buses --------------
  #
  # Only depends on local Vr_i, Vi_i:
  #   |V_i| = sqrt(Vr_i^2 + Vi_i^2)
  #   ∂|V_i|/∂Vr_i = Vr_i / |V_i|
  #   ∂|V_i|/∂Vi_i = Vi_i / |V_i|

  for i = 1:n
    if i == slack_idx || bus_types[i] != :PV
      continue
    end

    rb   = row_block[i]
    rowV = rb + 1  # second row for that bus

    pos = pos_non_slack[i]
    if pos == 0
      continue
    end

    vm = abs(V[i])
    if vm == 0.0
      continue
    end

    dVr = real(V[i]) / vm
    dVi = imag(V[i]) / vm

    colVr = pos
    colVi = (n - 1) + pos

    if abs(dVr) > 0.0
      push!(Iidx, rowV)
      push!(Jidx, colVr)
      push!(Vals, dVr)
    end
    if abs(dVi) > 0.0
      push!(Iidx, rowV)
      push!(Jidx, colVi)
      push!(Vals, dVi)
    end
  end

  # Local chain-rule terms for voltage-dependent specified injections.
  # For ΔP = Pcalc - Pspec(|V|): subtract dPspec/d|V| * d|V|/dVr and d|V|/dVi.
  # For PQ second row ΔQ = Qcalc - Qspec(|V|): analogous subtraction.
  for i in non_slack
    pos = pos_non_slack[i]
    rb = row_block[i]
    pos == 0 && continue
    rb == 0 && continue

    vm = abs(V[i])
    vm_safe = vm > vm_eps ? vm : vm_eps
    dvm_dvr = real(V[i]) / vm_safe
    dvm_dvi = imag(V[i]) / vm_safe
    colVr = pos
    colVi = (n - 1) + pos

    dP = dPinj_dVm[i]
    if dP != 0.0
      push!(Iidx, rb)
      push!(Jidx, colVr)
      push!(Vals, -dP * dvm_dvr)
      push!(Iidx, rb)
      push!(Jidx, colVi)
      push!(Vals, -dP * dvm_dvi)
    end

    if bus_types[i] == :PQ
      dQ = dQinj_dVm[i]
      if dQ != 0.0
        push!(Iidx, rb + 1)
        push!(Jidx, colVr)
        push!(Vals, -dQ * dvm_dvr)
        push!(Iidx, rb + 1)
        push!(Jidx, colVi)
        push!(Vals, -dQ * dvm_dvi)
      end
    end
  end

  return sparse(Iidx, Jidx, Vals, m, nvar)
end

"""
    build_rectangular_jacobian_pq_pv_dense(
        Ybus, V, bus_types, Vset, slack_idx
    ) -> J::Matrix{Float64}

Build the analytic rectangular Jacobian for the mismatch vector `F(V)`
defined in `mismatch_rectangular`.

- State vector: x = [Vr(non-slack); Vi(non-slack)]
- Rows: for each non-slack bus i
    * PQ: [ΔP_i; ΔQ_i]
    * PV: [ΔP_i; ΔV_i]  with ΔV_i = |V_i| - Vset[i]

`bus_types` and `Vset` must be consistent with `mismatch_rectangular`.
"""
function build_rectangular_jacobian_pq_pv_dense(Ybus, V::Vector{ComplexF64}, bus_types::Vector{Symbol}, Vset::Vector{Float64}, slack_idx::Int; dPinj_dVm::Vector{Float64} = zeros(Float64, length(V)), dQinj_dVm::Vector{Float64} = zeros(Float64, length(V)), vm_eps::Float64 = 1e-9)
  n = length(V)
  @assert length(bus_types) == n
  @assert length(Vset) == n
  @assert length(dPinj_dVm) == n
  @assert length(dQinj_dVm) == n

  # --- 1) Wirtinger blocks for S(V) = V .* conj(Ybus * V)
  J11, J12, J21, J22 = build_complex_jacobian(Ybus, V)

  # --- 2) Full 2n×2n rectangular J for ΔP/ΔQ wrt Vr/Vi (all buses)
  # Rows: [ΔP_1..ΔP_n; ΔQ_1..ΔQ_n]
  # Cols: [Vr_1..Vr_n; Vi_1..Vi_n]
  Jrect_full = zeros(Float64, 2n, 2n)

  @inbounds for j = 1:n
    col_sum  = J11[:, j] .+ J12[:, j]  # corresponds to dS/dVr_j
    col_diff = J11[:, j] .- J12[:, j]  # used for dS/dVi_j

    # dP/dVr_j, dQ/dVr_j
    @views Jrect_full[1:n, j]      .= real.(col_sum)
    @views Jrect_full[(n+1):2n, j] .= imag.(col_sum)

    # dP/dVi_j, dQ/dVi_j
    @views Jrect_full[1:n, n+j]      .= -imag.(col_diff)
    @views Jrect_full[(n+1):2n, n+j] .= real.(col_diff)
  end

  # --- 3) Reduce to non-slack variables and rows matching mismatch_rectangular

  non_slack = collect(1:n)
  deleteat!(non_slack, slack_idx)

  nvar = 2 * (n - 1)
  m    = 2 * (n - 1)
  @assert nvar == m

  pos_map = build_pos_map(non_slack, n)

  # Column indices in the full rectangular J that correspond to
  # [Vr(non_slack); Vi(non_slack)]
  col_idx_full = vcat(non_slack, n .+ non_slack)

  J = zeros(Float64, m, nvar)

  row = 1
  @inbounds for i = 1:n
    if i == slack_idx
      continue
    end

    # First row for this bus: ΔP_i
    rowP_full = i                  # P row index in full J
    @views J[row, :] .= Jrect_full[rowP_full, col_idx_full]

    # Second row: ΔQ_i (PQ) or ΔV_i (PV)
    if bus_types[i] == :PQ
      rowQ_full = n + i          # Q row index in full J
      @views J[row+1, :] .= Jrect_full[rowQ_full, col_idx_full]

    elseif bus_types[i] == :PV
      # ΔV_i = |V_i| - Vset[i]
      J[row+1, :] .= 0.0

      pos = pos_map[i]
      if pos != 0
        vm = abs(V[i])
        if vm > 0.0
          dVr = real(V[i]) / vm
          dVi = imag(V[i]) / vm

          # Columns in reduced J:
          #   Vr_i -> index pos
          #   Vi_i -> index (n-1) + pos
          J[row+1, pos]       = dVr
          J[row+1, (n-1)+pos] = dVi
        end
      end
    else
      error("build_rectangular_jacobian_pq_pv: unsupported bus type $(bus_types[i]) at bus $i")
    end

    row += 2
  end

  for i in non_slack
    rowP = 2 * pos_map[i] - 1
    pos = pos_map[i]
    vm = abs(V[i])
    vm_safe = vm > vm_eps ? vm : vm_eps
    dvm_dvr = real(V[i]) / vm_safe
    dvm_dvi = imag(V[i]) / vm_safe

    J[rowP, pos] -= dPinj_dVm[i] * dvm_dvr
    J[rowP, (n-1)+pos] -= dPinj_dVm[i] * dvm_dvi

    if bus_types[i] == :PQ
      rowQ = rowP + 1
      J[rowQ, pos] -= dQinj_dVm[i] * dvm_dvr
      J[rowQ, (n-1)+pos] -= dQinj_dVm[i] * dvm_dvi
    end
  end

  return J
end

"""
    build_rectangular_jacobian_pq_pv(
        Ybus,
        V,
        bus_types,
        Vset,
        slack_idx;
        use_sparse::Bool = false,
    )

Dispatches to either the dense or sparse rectangular Jacobian builder matching
`mismatch_rectangular`.

- If `use_sparse == true` and `Ybus` is a `SparseMatrixCSC{ComplexF64}`, the
  sparse builder is used.
- Otherwise, the dense builder is used.
"""
function build_rectangular_jacobian_pq_pv(Ybus, V::Vector{ComplexF64}, bus_types::Vector{Symbol}, Vset::Vector{Float64}, slack_idx::Int; use_sparse::Bool = false, dPinj_dVm::Vector{Float64} = zeros(Float64, length(V)), dQinj_dVm::Vector{Float64} = zeros(Float64, length(V)), vm_eps::Float64 = 1e-9)
  if use_sparse && Ybus isa SparseMatrixCSC{ComplexF64}
    return build_rectangular_jacobian_pq_pv_sparse(Ybus, V, bus_types, Vset, slack_idx; dPinj_dVm = dPinj_dVm, dQinj_dVm = dQinj_dVm, vm_eps = vm_eps)
  else
    return build_rectangular_jacobian_pq_pv_dense(Ybus, V, bus_types, Vset, slack_idx; dPinj_dVm = dPinj_dVm, dQinj_dVm = dQinj_dVm, vm_eps = vm_eps)
  end
end

"""
    complex_newton_step_rectangular(
        Ybus,
        V,
        S;
        slack_idx,
        damp,
        bus_types,
        Vset,
        use_sparse=false,
    )

Performs one Newton–Raphson step in rectangular coordinates using the analytic
Jacobian that matches `mismatch_rectangular`.

- State: x = [Vr(non-slack); Vi(non-slack)]
- Residual: F(x) = mismatch_rectangular(...)
"""
function complex_newton_step_rectangular(
  Ybus,
  V::Vector{ComplexF64},
  S::Vector{ComplexF64};
  slack_idx::Int,
  damp::Float64 = 1.0,
  autodamp::Bool = false,
  autodamp_min::Float64 = 1e-3,
  bus_types::Vector{Symbol},
  Vset::Vector{Float64},
  use_sparse::Bool = false,
  dPinj_dVm::Vector{Float64} = zeros(Float64, length(V)),
  dQinj_dVm::Vector{Float64} = zeros(Float64, length(V)),
)
  n = length(V)
  @assert length(S) == n
  @assert length(bus_types) == n
  @assert length(Vset) == n
  @assert length(dPinj_dVm) == n
  @assert length(dQinj_dVm) == n

  # Non-slack indices
  non_slack = non_slack_indices(n, slack_idx)
  # Residual matching the FD variant
  F0 = mismatch_rectangular(Ybus, V, S, bus_types, Vset, slack_idx)
  m = length(F0)
  nvar = 2 * (n - 1)
  @assert m == nvar "complex_newton_step_rectangular: mismatch and state dimension differ"

  # Analytic Jacobian (dense or sparse)
  J = build_rectangular_jacobian_pq_pv(Ybus, V, bus_types, Vset, slack_idx; use_sparse = use_sparse, dPinj_dVm = dPinj_dVm, dQinj_dVm = dQinj_dVm)

  # Solve J * δx = -F
  δx = solve_linear(J, -F0; allow_pinv = true)
  if autodamp
    _, Vtrial, _ = choose_rectangular_autodamp(Ybus, V, S, δx, F0; slack_idx = slack_idx, damp = damp, autodamp_min = autodamp_min, bus_types = bus_types, Vset = Vset)
    return Vtrial
  end

  _validate_rectangular_damping(damp, min(autodamp_min, damp))
  return _apply_rectangular_delta(V, δx, slack_idx, non_slack, damp)
end

"""
    run_complex_nr_rectangular_for_net!(net; maxiter=20, tol=1e-8, damp=0.2, verbose=0, use_fd=false)

Run a complex-state Newton-Raphson power flow in rectangular coordinates on a Sparlectra network.

# Arguments
- `net::Net`: Network object containing bus, branch, and generation data
- `maxiter::Int=20`: Maximum number of Newton-Raphson iterations
- `tol::Float64=1e-8`: Convergence tolerance for maximum mismatch
- `damp::Float64=0.2`: Damping factor for Newton step (0 < damp ≤ 1)
- `verbose::Int=0`: Verbosity level (0=quiet, 1=basic info, 2=detailed)
- `use_fd::Bool=false`: Use finite-difference Jacobian instead of analytic

# Returns
- `Tuple{Int, Int}`: (iterations_used, error_code)
  - error_code: 0=converged, 1=max_iterations_reached

# Details
This function implements a complete power flow solver using complex voltages in rectangular 
coordinates (Vr + jVi) as state variables, providing an alternative to conventional 
polar formulations.

## Mathematical Foundation
- **State vector**: x = [Vr(non-slack); Vi(non-slack)] ∈ ℝ^(2(n-1))
- **Complex power**: S = V .* conj(Y * V) where Y is the bus admittance matrix
- **PQ buses**: ΔP = Re(S_calc - S_spec), ΔQ = Im(S_calc - S_spec)
- **PV buses**: ΔP = Re(S_calc - S_spec), ΔV = |V| - V_set
- **Slack bus**: Voltage held constant throughout iterations

## Active Set Q-Limit Management
- **PV→PQ switching**: When reactive power demand violates generator Q-limits
- **Optional PQ→PV re-enable**: With hysteresis band and cooldown mechanisms
- **Robust handling**: Guards against inappropriate switching of non-generator buses

## Algorithm Steps
1. **Initialization**: Extract voltages, build Y-bus, classify bus types
2. **Power specification**: Build S = P + jQ from network loads/generation
3. **Iterative solution**: Newton-Raphson with mismatch function for PQ/PV constraints
4. **Q-limit enforcement**: Active-set management during iterations
5. **Result update**: Write final voltages and computed powers back to network

## Network Integration
- **Input**: Uses `net.nodeVec` for bus data, `net.baseMVA` for per-unit conversion
- **Output**: Updates `node._vm_pu`, `node._va_deg`, `node._pƩGen`, `node._qƩGen`
- **Q-limits**: Integrates with `net.qLimitEvents` and `net.qLimitLog` for tracking
- **Compatibility**: Maintains same interface as `runpf!()` for easy substitution

# See Also
- `runpf_rectangular!()`: Convenience wrapper matching `runpf!()` signature
- `mismatch_rectangular()`: Core mismatch function for PQ/PV constraints
- `build_rectangular_jacobian_pq_pv()`: Analytic Jacobian construction
"""

function run_complex_nr_rectangular_for_net!(
  net::Net;
  maxiter::Int = 20,
  tol::Float64 = 1e-8,
  damp::Float64 = 0.2,
  verbose::Int = 0,
  use_fd::Bool = false,
  autodamp::Bool = false,
  autodamp_min::Float64 = 1e-3,
  opt_sparse::Bool = true,
  opt_flatstart::Bool = net.flatstart,
  pv_table_rows::Int = 30,
  lock_pv_to_pq_buses::AbstractVector{Int} = Int[],
  qlimit_mode::Symbol = :switch_to_pq,
  qlimit_max_outer::Int = 30,
  start_projection::Bool = false,
  start_projection_try_dc_start::Bool = true,
  start_projection_try_blend_scan::Bool = true,
  start_projection_blend_lambdas::AbstractVector{<:Real} = [0.25, 0.5, 0.75],
  start_projection_dc_angle_limit_deg::Float64 = 60.0,
  qlimit_start_iter::Int = 2,
  qlimit_start_mode::Symbol = :iteration,
  qlimit_auto_q_delta_pu::Float64 = 1e-4,
)
  if verbose > 1
    @info "Running complex rectangular NR power flow... use_fd=$use_fd, opt_sparse=$opt_sparse"
  end

  nodes = net.nodeVec
  n     = length(nodes)
  Sbase = net.baseMVA
  if !opt_sparse
    sparse = n > 1000
  else
    sparse = opt_sparse
  end
  Yred = createYBUS(net = net, sparse = sparse, printYBUS = (verbose > 1))
  Ybus = (size(Yred, 1) == n) ? Yred : _expand_ybus_for_isolated_nodes(Yred, n, net.isoNodes)

  # 1) Initial complex voltages V0 and slack index
  V0, slack_idx = initialVrect(net; flatstart = opt_flatstart)

  # 2) Specified complex power injections S (p.u.). For Q(U)/P(U) controllers
  # this vector becomes state-dependent and is re-evaluated per Newton iteration.
  S = buildComplexSVec(net)
  dPinj_dVm = zeros(Float64, n)
  dQinj_dVm = zeros(Float64, n)
  has_vdep_control = has_voltage_dependent_control(net)

  # 3) Bus types and PV setpoints from Node data
  bus_types = Vector{Symbol}(undef, n)
  Vset      = Vector{Float64}(undef, n)

  @inbounds for (k, node) in enumerate(nodes)
    BusType = getNodeType(node)
    if BusType == Slack
      bus_types[k] = :Slack
    elseif BusType == PV
      bus_types[k] = :PV
    elseif BusType == PQ
      bus_types[k] = :PQ
    elseif BusType == Isolated
      # Keep isolated buses in the rectangular state vector as neutral PQ rows.
      # Their injections are forced to zero so they do not affect the solved grid.
      bus_types[k] = :PQ
      S[k] = 0.0 + 0.0im
    else
      error("run_complex_nr_rectangular_for_net!: unsupported bus type at bus $k, given: $(BusType)")
    end

    Vset[k] = isnothing(node._vm_pu) ? 1.0 : node._vm_pu
  end

  V0 = project_rectangular_start(
    Ybus,
    V0,
    S,
    bus_types,
    Vset,
    slack_idx;
    enabled = start_projection,
    try_dc_start = start_projection_try_dc_start,
    try_blend_scan = start_projection_try_blend_scan,
    blend_lambdas = start_projection_blend_lambdas,
    dc_angle_limit_deg = start_projection_dc_angle_limit_deg,
    verbose = verbose,
  )

  # 4) Q-limit data 
  qmin_pu, qmax_pu = getQLimits_pu(net)
  if verbose > 1
    printPVQLimitsTable(net; max_rows = typemax(Int))
  elseif verbose > 0
    printPVQLimitsTable(net; max_rows = pv_table_rows)
  end

  # --- Active-set bookkeeping (rectangular solver) ------------------------
  nb = n  # number of buses

  # PV origin mask (guards for PQ->PV re-enable)
  pv_orig_mask = falses(nb)
  @inbounds for k = 1:nb
    pv_orig_mask[k] = (bus_types[k] == :PV)
  end

  cooldown_iters = hasfield(typeof(net), :cooldown_iters) ? net.cooldown_iters : 0
  q_hyst_pu      = hasfield(typeof(net), :q_hyst_pu) ? net.q_hyst_pu : 0.0
  allow_reenable = (cooldown_iters > 0) || (q_hyst_pu > 0.0)
  qlimit_mode in (:switch_to_pq, :adjust_vset) || error("Unsupported qlimit_mode=$(qlimit_mode). Supported: :switch_to_pq, :adjust_vset.")
  qlimit_max_outer > 0 || error("qlimit_max_outer must be > 0 (got $(qlimit_max_outer)).")
  qlimit_start_iter > 0 || error("qlimit_start_iter must be > 0 (got $(qlimit_start_iter)).")
  qlimit_start_mode in (:iteration, :auto_q_delta, :iteration_or_auto) || error("Unsupported qlimit_start_mode=$(qlimit_start_mode). Supported: :iteration, :auto_q_delta, :iteration_or_auto.")
  qlimit_auto_q_delta_pu >= 0.0 || error("qlimit_auto_q_delta_pu must be >= 0 (got $(qlimit_auto_q_delta_pu)).")

  controllers = qlimit_mode == :adjust_vset ? _build_vset_adjust_controllers(net) : Dict{Int,NamedTuple{(:prosumer_idx, :config),Tuple{Int,VoltageAdjustConfig}}}()
  base_vset = copy(Vset)
  adjust_counter = zeros(Int, nb)

  # Start fresh each PF run
  resetQLimitLog!(net)

  # 5) NR-Loop
  V         = copy(V0)
  history   = Float64[]
  converged = false
  iters     = 0
  prev_pv_qreq_pu = fill(NaN, nb)

  if verbose > 1
    @info "Starting rectangular complex NR power flow..."
    @info "Initial complex voltages V0:" V0
    @info "Slack bus index:" slack_idx
    @info "maxiter = $maxiter, tol = $tol, damp = $damp, autodamp = $autodamp, autodamp_min = $autodamp_min, start_projection = $start_projection"
  end

  for it = 1:maxiter
    iters = it

    if has_vdep_control
      S, dPinj_dVm, dQinj_dVm = buildControlledSVec(net, V)
    end

    # Mismatch with current bus_types and (possibly) voltage-dependent S.
    F = mismatch_rectangular(Ybus, V, S, bus_types, Vset, slack_idx)
    max_mis = maximum(abs.(F))
    push!(history, max_mis)

    (verbose > 1) && @debug "Rectangular NR iteration" iter = it max_mismatch = max_mis

    if max_mis <= tol
      converged = true
      break
    end

    # --- Q-Limit Active Set: PV -> PQ, optional PQ -> PV (rectangular) ------
    changed   = false
    reenabled = false

    Qload_pu = build_qload_pu(net)

    Scalc_pu = calc_injections(Ybus, V)
    current_pv_qreq_pu = fill(NaN, nb)
    @inbounds for bus in eachindex(current_pv_qreq_pu)
      if bus_types[bus] == :PV
        current_pv_qreq_pu[bus] = imag(Scalc_pu[bus]) + Qload_pu[bus]
      end
    end

    qlimit_iter_ready = it >= qlimit_start_iter
    qlimit_auto_ready = false
    if qlimit_start_mode in (:auto_q_delta, :iteration_or_auto)
      max_q_delta = 0.0
      compared = false
      @inbounds for bus in eachindex(current_pv_qreq_pu)
        if isfinite(current_pv_qreq_pu[bus]) && isfinite(prev_pv_qreq_pu[bus])
          max_q_delta = max(max_q_delta, abs(current_pv_qreq_pu[bus] - prev_pv_qreq_pu[bus]))
          compared = true
        end
      end
      qlimit_auto_ready = compared && (max_q_delta <= qlimit_auto_q_delta_pu)
    end
    qlimit_ready = qlimit_start_mode == :iteration ? qlimit_iter_ready :
                   qlimit_start_mode == :auto_q_delta ? qlimit_auto_ready :
                   (qlimit_iter_ready || qlimit_auto_ready)

    if qlimit_ready
      changed, reenabled = active_set_q_limits!(
        net,
        it,
        nb;
        qmin_pu = qmin_pu,
        qmax_pu = qmax_pu,
        pv_orig_mask = pv_orig_mask,
        allow_reenable = qlimit_mode == :switch_to_pq ? allow_reenable : false,
        q_hyst_pu = q_hyst_pu,
        cooldown_iters = cooldown_iters,
        lock_pv_to_pq_buses = lock_pv_to_pq_buses,
        on_violation! = qlimit_mode == :adjust_vset ? ((bus, qreq, side, qclamp) -> _try_adjust_vset_on_q_limit!(net, bus, side, it, controllers, base_vset, Vset, adjust_counter, qlimit_max_outer, verbose)) : nothing,
        verbose = verbose,
        get_qreq_pu = bus -> begin
          (bus_types[bus] == :Slack) && return 0.0
          return imag(Scalc_pu[bus]) + Qload_pu[bus]
        end,
        is_pv = bus -> (bus_types[bus] == :PV),
        make_pq! = (bus, qclamp_gen_pu, side) -> begin
          bus_types[bus] = :PQ
          qinj_pu = qclamp_gen_pu - Qload_pu[bus]
          S[bus] = ComplexF64(real(S[bus]), qinj_pu)
          net.nodeVec[bus]._qƩGen = qclamp_gen_pu * net.baseMVA
        end,
        make_pv! = (bus) -> begin
          bus_types[bus] = :PV
        end,
      )

      # If bus_types/spec changed, mismatch definition changed (ΔQ ↔ ΔV) => rebuild F
      if changed || reenabled
        F = mismatch_rectangular(Ybus, V, S, bus_types, Vset, slack_idx)
        max_mis = maximum(abs.(F))
        history[end] = max_mis  # optional: overwrite last stored value for this iteration
      end
    end
    prev_pv_qreq_pu = current_pv_qreq_pu
    # --- Newton-Stepp (FD oder analytisch) -----------------------------
    if use_fd
      V = complex_newton_step_rectangular_fd(Ybus, V, S; slack_idx = slack_idx, damp = damp, autodamp = autodamp, autodamp_min = autodamp_min, h = 1e-6, bus_types = bus_types, Vset = Vset, dPinj_dVm = dPinj_dVm, dQinj_dVm = dQinj_dVm)
    else
      V = complex_newton_step_rectangular(Ybus, V, S; slack_idx = slack_idx, damp = damp, autodamp = autodamp, autodamp_min = autodamp_min, bus_types = bus_types, Vset = Vset, use_sparse = sparse, dPinj_dVm = dPinj_dVm, dQinj_dVm = dQinj_dVm)
    end

    # Keep slack voltage fixed (safety belt)
    V[slack_idx] = V0[slack_idx]
  end

  # 6) Update voltages back to network
  # --- mirror bus_types back into Net/node types (PV->PQ switching) ---
  @inbounds for k = 1:n
    if bus_types[k] == :PQ
      setNodeType!(net.nodeVec[k], "PQ")
    elseif bus_types[k] == :PV
      setNodeType!(net.nodeVec[k], "PV")
    end
  end
  update_net_voltages_from_complex!(net, V)

  # 7) Compute bus injections from final voltages  
  Sbus_pu = calc_injections(Ybus, V)
  Sbus_MVA = Sbus_pu .* Sbase

  @debug "Final Voltages Mag = " [abs.(V)...]
  @debug "Final Voltages Ang = " [angle.(V) .* (180.0 / π)...]

  isoNodes = net.isoNodes

  @inbounds for (k, node) in enumerate(nodes)
    # Skip isolated nodes (if any)
    if k in isoNodes
      continue
    end

    Sbus      = Sbus_MVA[k]
    Pbus_MW   = real(Sbus)
    Qbus_MVar = imag(Sbus)

    # Slack bus: always write P and Q generation from the solved state
    if node._nodeType == Sparlectra.Slack
      node._pƩGen = Pbus_MW
      node._qƩGen = Qbus_MVar

    elseif node._nodeType == Sparlectra.PV
      # PV bus: Qgen is a result (unless it later switches)
      node._qƩGen = Qbus_MVar

    elseif node._nodeType == Sparlectra.PQ
      # If this bus was forced PV->PQ by Q-limit, keep the clamped generator Q
      # (already written in make_pq! as MVar).
      if haskey(net.qLimitEvents, k)
        @debug "Bus $(k) is PQ due to Q-limit; keeping clamped Qgen = $(node._qƩGen) MVar."
      end
    end
    # PQ buses / pure loads: do not touch _pƩLoad / _pƩGen here.
    # The original load/generation specification remains intact.
  end

  # 8) Update total bus power (sum of complex injections in p.u.)
  p = (sum(real.(Sbus_pu))) * Sbase
  q = (sum(imag.(Sbus_pu))) * Sbase

  if verbose > 1
    @info "Set total bus power to p = $p MW and q = $q MVar"
  end

  setTotalBusPower!(net = net, p = p, q = q)
  updateShuntPowers!(net = net)

  return iters, converged ? 0 : 1
end

"""
    runpf_rectangular!(net, maxIte, tolerance=1e-6, verbose=0)

Runs a rectangular complex-state Newton–Raphson power flow on `net::Net`.

Returns:
    (iterations::Int, status::Int)
where `status == 0` indicates convergence.
"""
function runpf_rectangular!(
  net::Net,
  maxIte::Int,
  tolerance::Float64 = 1e-6,
  verbose::Int = 0;
  opt_fd::Bool = false,
  opt_sparse::Bool = true,
  damp = 1.0,
  autodamp::Bool = false,
  autodamp_min::Float64 = 1e-3,
  opt_flatstart::Bool = net.flatstart,
  pv_table_rows::Int = 30,
  lock_pv_to_pq_buses::AbstractVector{Int} = Int[],
  qlimit_mode::Symbol = :switch_to_pq,
  qlimit_max_outer::Int = 30,
  start_projection::Bool = false,
  start_projection_try_dc_start::Bool = true,
  start_projection_try_blend_scan::Bool = true,
  start_projection_blend_lambdas::AbstractVector{<:Real} = [0.25, 0.5, 0.75],
  start_projection_dc_angle_limit_deg::Float64 = 60.0,
  qlimit_start_iter::Int = 2,
  qlimit_start_mode::Symbol = :iteration,
  qlimit_auto_q_delta_pu::Float64 = 1e-4,
)
  iters, erg = run_complex_nr_rectangular_for_net!(
    net;
    maxiter = maxIte,
    tol = tolerance,
    damp = damp,
    autodamp = autodamp,
    autodamp_min = autodamp_min,
    verbose = verbose,
    use_fd = opt_fd,
    opt_sparse = opt_sparse,
    opt_flatstart = opt_flatstart,
    pv_table_rows = pv_table_rows,
    lock_pv_to_pq_buses = lock_pv_to_pq_buses,
    qlimit_mode = qlimit_mode,
    qlimit_max_outer = qlimit_max_outer,
    start_projection = start_projection,
    start_projection_try_dc_start = start_projection_try_dc_start,
    start_projection_try_blend_scan = start_projection_try_blend_scan,
    start_projection_blend_lambdas = start_projection_blend_lambdas,
    start_projection_dc_angle_limit_deg = start_projection_dc_angle_limit_deg,
    qlimit_start_iter = qlimit_start_iter,
    qlimit_start_mode = qlimit_start_mode,
    qlimit_auto_q_delta_pu = qlimit_auto_q_delta_pu,
  )
  return iters, erg
end

function _active_link_representative_map(net::Net)
  n = length(net.nodeVec)
  parent = collect(1:n)

  function find_root(i::Int)
    while parent[i] != i
      parent[i] = parent[parent[i]]
      i = parent[i]
    end
    return i
  end

  function union_set(a::Int, b::Int)
    ra = find_root(a)
    rb = find_root(b)
    ra == rb && return
    if ra < rb
      parent[rb] = ra
    else
      parent[ra] = rb
    end
  end

  for l in net.linkVec
    l.status == 1 || continue
    union_set(Int(l.fromBus), Int(l.toBus))
  end

  return [find_root(i) for i = 1:n]
end

function _merged_pf_net(net::Net)
  reps = _active_link_representative_map(net)
  all(reps[i] == i for i in eachindex(reps)) && return net, reps, false

  wnet = deepcopy(net)
  n = length(wnet.nodeVec)
  cluster_members = [Int[] for _ = 1:n]
  for bus = 1:n
    push!(cluster_members[reps[bus]], bus)
  end

  for rep = 1:n
    members = cluster_members[rep]
    isempty(members) && continue

    ref = wnet.nodeVec[rep]
    p_load = 0.0
    q_load = 0.0
    p_gen = 0.0
    q_gen = 0.0
    p_sh = 0.0
    q_sh = 0.0
    has_slack = false
    has_pv = false

    for b in members
      nref = wnet.nodeVec[b]
      p_load += isnothing(nref._pƩLoad) ? 0.0 : nref._pƩLoad
      q_load += isnothing(nref._qƩLoad) ? 0.0 : nref._qƩLoad
      p_gen += isnothing(nref._pƩGen) ? 0.0 : nref._pƩGen
      q_gen += isnothing(nref._qƩGen) ? 0.0 : nref._qƩGen
      p_sh += isnothing(nref._pShunt) ? 0.0 : nref._pShunt
      q_sh += isnothing(nref._qShunt) ? 0.0 : nref._qShunt
      has_slack |= (nref._nodeType == Slack)
      has_pv |= (nref._nodeType == PV)
    end

    ref._pƩLoad = p_load
    ref._qƩLoad = q_load
    ref._pƩGen = p_gen
    ref._qƩGen = q_gen
    ref._pShunt = p_sh
    ref._qShunt = q_sh
    ref._nodeType = has_slack ? Slack : (has_pv ? PV : PQ)

    for b in members
      b == rep && continue
      nref = wnet.nodeVec[b]
      nref._pƩLoad = 0.0
      nref._qƩLoad = 0.0
      nref._pƩGen = 0.0
      nref._qƩGen = 0.0
      nref._pShunt = 0.0
      nref._qShunt = 0.0
      nref._nodeType = Isolated
    end
  end

  for br in wnet.branchVec
    f = reps[Int(br.fromBus)]
    t = reps[Int(br.toBus)]
    br.fromBus = f
    br.toBus = t
    if f == t
      br.status = 0
    end
  end

  for sh in wnet.shuntVec
    sh.busIdx = reps[Int(sh.busIdx)]
    if hasproperty(sh.comp, :cFrom_bus) && getfield(sh.comp, :cFrom_bus) !== nothing
      setfield!(sh.comp, :cFrom_bus, sh.busIdx)
    end
  end

  for ps in wnet.prosumpsVec
    if hasproperty(ps.comp, :cFrom_bus) && getfield(ps.comp, :cFrom_bus) !== nothing
      new_bus = reps[Int(getfield(ps.comp, :cFrom_bus))]
      setfield!(ps.comp, :cFrom_bus, new_bus)
    end
    if hasproperty(ps.comp, :cTo_bus) && getfield(ps.comp, :cTo_bus) !== nothing
      new_bus = reps[Int(getfield(ps.comp, :cTo_bus))]
      setfield!(ps.comp, :cTo_bus, new_bus)
    end
  end

  empty!(wnet.isoNodes)
  markIsolatedBuses!(net = wnet, log = false)
  return wnet, reps, true
end

"""
    runpf!(net, maxIte, tolerance=1e-6, verbose=0; method=:rectangular)

Unified AC power flow interface.

Arguments:
- `net::Net`: network
- `maxIte::Int`: maximum iterations
- `tolerance::Float64`: mismatch tolerance
- `verbose::Int`: verbosity level
- `method::Symbol`: `:rectangular` (recommended), `:polar_full` (deprecated), or `:classic` (deprecated)
- `autodamp::Bool`: enable residual-based backtracking for rectangular Newton steps
- `autodamp_min::Float64`: minimum automatic damping factor when `autodamp = true`
- `qlimit_start_iter::Int`: first Newton iteration where PV→PQ Q-limit switching may run in `:iteration` mode
- `qlimit_start_mode::Symbol`: `:iteration`, `:auto_q_delta`, or `:iteration_or_auto` start criterion for PV→PQ switching
- `qlimit_auto_q_delta_pu::Float64`: PV reactive-power request change threshold for automatic switching start

Notes:
- Link-flow recovery (`calcLinkFlowsKCL!`) is method-agnostic and uses solved PF results.
- If active-link merges create internal isolated buses, `:rectangular` currently falls
  back to `:polar_full` for robustness.

Returns:
    (iterations::Int, status::Int)

where `status == 0` indicates convergence.
"""
function runpf!(
  net::Net,
  maxIte::Int,
  tolerance::Float64 = 1e-6,
  verbose::Int = 0;
  method::Symbol = :rectangular,
  opt_fd::Bool = false,
  opt_sparse::Bool = true,
  opt_flatstart::Bool = net.flatstart,
  damp = 1.0,
  autodamp::Bool = false,
  autodamp_min::Float64 = 1e-3,
  pv_table_rows::Int = 30,
  validate_limits_after_pf::Bool = false,
  q_limit_violation_headroom::Float64 = 0.0,
  lock_pv_to_pq_buses::AbstractVector{Int} = Int[],
  qlimit_mode::Symbol = :switch_to_pq,
  qlimit_max_outer::Int = 30,
  start_projection::Bool = false,
  start_projection_try_dc_start::Bool = true,
  start_projection_try_blend_scan::Bool = true,
  start_projection_blend_lambdas::AbstractVector{<:Real} = [0.25, 0.5, 0.75],
  start_projection_dc_angle_limit_deg::Float64 = 60.0,
  qlimit_start_iter::Int = 2,
  qlimit_start_mode::Symbol = :iteration,
  qlimit_auto_q_delta_pu::Float64 = 1e-4,
)
  wnet, reps, has_merges = _merged_pf_net(net)
  refreshBusTypesFromProsumers!(wnet)
  has_vdep_control = has_voltage_dependent_control(wnet)

  function _sync_merged_results_to_original!()
    for i in eachindex(net.nodeVec)
      src = wnet.nodeVec[reps[i]]
      net.nodeVec[i]._vm_pu = src._vm_pu
      net.nodeVec[i]._va_deg = src._va_deg
    end
    updateShuntPowers!(net = net)
  end

  #@info "Running AC Power Flow using method: $(method)"
  if method === :polar_full
    has_vdep_control && error("runpf!: voltage-dependent injections, including P(U)/Q(U) controllers and bus_shunt_model=voltage_dependent_injection, are currently supported only for method=:rectangular.")
    if qlimit_mode != :switch_to_pq
      @warn "runpf!: qlimit_mode=$(qlimit_mode) is only supported for method=:rectangular. Falling back to :switch_to_pq behavior."
    end
    iters, erg = runpf_full!(wnet, maxIte, tolerance, verbose; opt_sparse = opt_sparse, opt_flatstart = opt_flatstart, pv_table_rows = pv_table_rows, lock_pv_to_pq_buses = lock_pv_to_pq_buses, warn_deprecated = true)
    if erg == 0 && has_merges
      _sync_merged_results_to_original!()
    end
    if validate_limits_after_pf && (verbose > 0)
      printFinalLimitValidation(has_merges ? net : wnet; q_headroom = q_limit_violation_headroom)
    end
    return iters, erg
  elseif method === :rectangular
    if has_merges
      has_vdep_control && error("runpf!: voltage-dependent injections, including P(U)/Q(U) controllers and bus_shunt_model=voltage_dependent_injection, are not supported with active-link merge handling in rectangular mode. Disable merges or use a topology without internal isolated buses.")
      if verbose > 0
        @warn "runpf!: rectangular solver detected internal Isolated buses from active-link merges; using rectangular FD fallback instead of :polar_full"
      end
      iters, erg = runpf_rectangular!(wnet, maxIte, tolerance, verbose; opt_fd = true, opt_sparse = opt_sparse, damp = damp, autodamp = autodamp, autodamp_min = autodamp_min, opt_flatstart = opt_flatstart, pv_table_rows = pv_table_rows, lock_pv_to_pq_buses = lock_pv_to_pq_buses, qlimit_mode = qlimit_mode, qlimit_max_outer = qlimit_max_outer, start_projection = start_projection, start_projection_try_dc_start = start_projection_try_dc_start, start_projection_try_blend_scan = start_projection_try_blend_scan, start_projection_blend_lambdas = start_projection_blend_lambdas, start_projection_dc_angle_limit_deg = start_projection_dc_angle_limit_deg, qlimit_start_iter = qlimit_start_iter, qlimit_start_mode = qlimit_start_mode, qlimit_auto_q_delta_pu = qlimit_auto_q_delta_pu)
    else
      iters, erg = runpf_rectangular!(wnet, maxIte, tolerance, verbose; opt_fd = opt_fd, opt_sparse = opt_sparse, damp = damp, autodamp = autodamp, autodamp_min = autodamp_min, opt_flatstart = opt_flatstart, pv_table_rows = pv_table_rows, lock_pv_to_pq_buses = lock_pv_to_pq_buses, qlimit_mode = qlimit_mode, qlimit_max_outer = qlimit_max_outer, start_projection = start_projection, start_projection_try_dc_start = start_projection_try_dc_start, start_projection_try_blend_scan = start_projection_try_blend_scan, start_projection_blend_lambdas = start_projection_blend_lambdas, start_projection_dc_angle_limit_deg = start_projection_dc_angle_limit_deg, qlimit_start_iter = qlimit_start_iter, qlimit_start_mode = qlimit_start_mode, qlimit_auto_q_delta_pu = qlimit_auto_q_delta_pu)
    end
    if erg == 0 && has_merges
      _sync_merged_results_to_original!()
    end
    if validate_limits_after_pf && (verbose > 0)
      printFinalLimitValidation(has_merges ? net : wnet; q_headroom = q_limit_violation_headroom)
    end
    return iters, erg
  elseif method === :classic
    has_vdep_control && error("runpf!: voltage-dependent injections, including P(U)/Q(U) controllers and bus_shunt_model=voltage_dependent_injection, are currently supported only for method=:rectangular.")
    if qlimit_mode != :switch_to_pq
      @warn "runpf!: qlimit_mode=$(qlimit_mode) is only supported for method=:rectangular. Falling back to :switch_to_pq behavior."
    end
    iters, erg = runpf_classic!(wnet, maxIte, tolerance, verbose, opt_sparse, opt_flatstart)
    if erg == 0 && has_merges
      _sync_merged_results_to_original!()
    end
    if validate_limits_after_pf && (verbose > 0)
      printFinalLimitValidation(has_merges ? net : wnet; q_headroom = q_limit_violation_headroom)
    end
    return iters, erg
  else
    error("runpf!: unknown method $(method). Use :rectangular (recommended), :polar_full (deprecated), or :classic (deprecated).")
  end
end
