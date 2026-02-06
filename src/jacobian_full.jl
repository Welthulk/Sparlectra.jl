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

# file: src/jacobian_full.jl
# jacobian_full.jl — Full-system Newton-Raphson including PV identity rows
#
# Features:
# - Full-state NR: θ and V for all non-slack buses
# - PV buses: second equation is replaced by identity row enforcing rPV := V - Vset = 0
# - Active-set handling of PV Q-limits: PV→PQ when Q demand violates limits
# - Optional re-enable PQ→PV with hysteresis band and cooldown (robust to missing fields)
#
# Robustness:
# - If the Net type does not have fields `q_hyst_pu` or `cooldown_iters`, defaults (0.0, 0) are used.
# - If Q-limit arrays don't match the current bus count, we (re)build once.

using LinearAlgebra
using SparseArrays
using Printf

# ------------------------------
# Power feed vector for the full system
# (P for all buses, Q only for PQ; PV gets placeholder rows)
# ------------------------------
function getPowerFeeds_full(busVec::Vector{BusData}, n_pq::Int, n_pv::Int)::Vector{Float64}
  n = n_pq + n_pv
  pf = zeros(Float64, 2n)
  vindex = 0
  @inbounds for bus in busVec
    if bus.type == Sparlectra.PQ
      vindex += 1;
      pf[vindex] = bus.pƩ   # P
      vindex += 1;
      pf[vindex] = bus.qƩ   # Q
    elseif bus.type == Sparlectra.PV
      vindex += 1;
      pf[vindex] = bus.pƩ   # P
      vindex += 1;
      pf[vindex] = 0.0      # rPV placeholder
    end
  end
  return pf
end

# ------------------------------
# Residual for the full system: PV’s second equation becomes rPV = V − Vset
# ------------------------------
"""
    residuum_state_full_withPV(
        Y, busVec, Vset, n_pq, n_pv, log
    ) -> NamedTuple

Compute the full-system residual Δ for the Newton-Raphson iteration including
PV identity rows, and additionally return the complex bus voltages and bus
power injections.

Returns a `NamedTuple` with:
- `Δ::Vector{Float64}`  : stacked active/reactive residuals (P/Q and rPV)
- `V::Vector{ComplexF64}` : complex bus voltages V = Vm * exp(j*Va)
- `S::Vector{ComplexF64}` : complex bus powers S = P + jQ = V .* conj(Y*V)
"""
# Extended residual function used internally by the NR solver.
# Besides the residual Δ, it returns complex bus voltages (V_nr) and nodal
# powers (S_nr). This enables downstream modules (branch flows, losses,
# diagnostics) to use exactly the same voltage vector as the NR residual
# calculation.
function residuum_state_full_withPV(
  Y::AbstractMatrix{ComplexF64},
  busVec::Vector{BusData},
  Vset::Vector{Float64},           # length = number of buses (only relevant for PV)
  n_pq::Int,
  n_pv::Int,
  log::Bool,
)

  # Complex bus voltages
  V = [bus.vm_pu * exp(im * bus.va_rad) for bus in busVec]

  # Bus power injections from network model  
  S = calc_injections(Y, V)

  n = n_pq + n_pv
  Δ = zeros(Float64, 2n)

  vindex = 0
  @inbounds for k in eachindex(busVec)
    bus = busVec[k]
    if bus.type == Sparlectra.PQ
      vindex += 1;
      Δ[vindex] = bus.pƩ - real(S[k])        # ΔP
      vindex += 1;
      Δ[vindex] = bus.qƩ - imag(S[k])        # ΔQ
    elseif bus.type == Sparlectra.PV
      vindex += 1;
      Δ[vindex] = bus.pƩ - real(S[k])        # ΔP
      vindex += 1;
      Δ[vindex] = bus.vm_pu - Vset[k]        # rPV
    end
    # Store calculated bus powers for later use in the Net object
    bus._pRes = real(S[k])
    bus._qRes = imag(S[k])
  end

  if log
    println("residuum_full_withPV: ‖Δ‖ = ", norm(Δ))
  end

  return (Δ = Δ, V = V, S = S)
end

"""
    residuum_full_withPV(
        Y, busVec, Vset, n_pq, n_pv, log
    ) -> Vector{Float64}

Compatibility wrapper that returns only the residual vector Δ.

Internally calls `residuum_state_full_withPV` which also provides the
complex bus voltages and bus powers for further post-processing.
"""
# Compatibility wrapper that returns only the residual vector Δ.
# Internally calls `residuum_state_full_withPV`, which also provides the
# complex bus voltages and bus powers for further post-processing. This
# simpler API is kept for existing solver code.

function residuum_full_withPV(Y::AbstractMatrix{ComplexF64}, busVec::Vector{BusData}, Vset::Vector{Float64}, n_pq::Int, n_pv::Int, log::Bool)::Vector{Float64}
  state = residuum_state_full_withPV(Y, busVec, Vset, n_pq, n_pv, log)
  return state.Δ
end

# ------------------------------
# Full Jacobian matrix including PV voltage columns and PV identity rows
# ------------------------------
function calcJacobian_withPVIdentity(Y::AbstractMatrix{ComplexF64}, busVec::Vector{BusData}, adjBranch::Vector{Vector{Int}}, slackIdx::Int, n_pq::Int, n_pv::Int; log::Bool = false, sparse::Bool = true)

  # Remove only the slack bus from θ/V as usual; do not eliminate anything else
  shiftIJ(idx::Int) = (idx >= slackIdx) ? idx - 1 : idx

  n = n_pq + n_pv
  dim = 2n

  J = sparse ? spzeros(Float64, dim, dim) : zeros(Float64, dim, dim)

  @inbounds for i in eachindex(busVec)
    # SKIP: never build Jacobian rows for the slack bus
    if busVec[i].type == Sparlectra.Slack
      continue
    end

    vm_i = busVec[i].vm_pu
    va_i = busVec[i].va_rad
    type_i = busVec[i].type

    i1 = 2*shiftIJ(i)-1   # P-row i
    i2 = 2*shiftIJ(i)     # Q/constraint row i

    for j in adjBranch[i]
      # SKIP: never build Jacobian columns for the slack bus
      if busVec[j].type == Sparlectra.Slack
        continue
      end

      vm_j = busVec[j].vm_pu
      va_j = busVec[j].va_rad

      j1 = 2*shiftIJ(j)-1   # θ-column j
      j2 = 2*shiftIJ(j)     # V-column j (now also for PV buses)

      if i == j
        Yii = abs(Y[i, i]);
        αii = angle(Y[i, i]);
        arg = -αii
        # P-θ
        J[i1, j1] = vm_i*(Yii*vm_i*sin(arg)) - vm_i*sum(abs(Y[i, k]) * busVec[k].vm_pu * sin(va_i - busVec[k].va_rad - angle(Y[i, k])) for k in adjBranch[i])
        # Q-θ
        J[i2, j1] = -vm_i*(Yii*vm_i*cos(arg)) + vm_i*sum(abs(Y[i, k]) * busVec[k].vm_pu * cos(va_i - busVec[k].va_rad - angle(Y[i, k])) for k in adjBranch[i])
        # P-V (relative)
        J[i1, j2] = vm_i*(Yii*vm_i*cos(arg)) + vm_i*sum(abs(Y[i, k]) * busVec[k].vm_pu * cos(va_i - busVec[k].va_rad - angle(Y[i, k])) for k in adjBranch[i])
        # Q-V (relative)
        J[i2, j2] = vm_i*(Yii*vm_i*sin(arg)) + vm_i*sum(abs(Y[i, k]) * busVec[k].vm_pu * sin(va_i - busVec[k].va_rad - angle(Y[i, k])) for k in adjBranch[i])
      else
        Yij = abs(Y[i, j]);
        arg = va_i - va_j - angle(Y[i, j])
        # P-θ_j
        J[i1, j1] = vm_i * Yij * vm_j * sin(arg)
        # Q-θ_j
        J[i2, j1] = -vm_i * Yij * vm_j * cos(arg)
        # P-V_j (relative)
        J[i1, j2] = vm_i * Yij * vm_j * cos(arg)
        # Q-V_j (relative)
        J[i2, j2] = vm_i * Yij * vm_j * sin(arg)
      end
    end

    # PV constraint: row i2 as (V - Vset) = 0 => ∂/∂ΔV_rel_i = +V_i
    if type_i == Sparlectra.PV
      if issparse(J)
        zero_row!(J, i2)               # sparse path
      else
        @views J[i2, :] .= 0.0         # dense path
      end
      J[i2, 2*shiftIJ(i)] = vm_i
    end
  end

  if log
    print_jacobian(J)
  end

  return J
end

# ------------------------------
# Full Newton-Raphson including PV identity rows (separate from the reduced version)
# ------------------------------
function calcNewtonRaphson_withPVIdentity!(net::Net, Y::AbstractMatrix{ComplexF64}, maxIte::Int; tolerance::Float64 = 1e-6, verbose::Int = 0, sparse::Bool = false, flatStart::Bool = false, angle_limit::Bool = false, debug::Bool = false)
  if verbose > 0
    @info "Running full-system Newton-Raphson with PV identity rows...sparse=$sparse, flatStart=$flatStart, angle_limit=$angle_limit, debug=$debug"
  end
  
  setJacobianAngleLimit(angle_limit)
  setJacobianDebug(debug)

  # Defaults if these fields do not exist on the Net type
  cooldown_iters = hasfield(typeof(net), :cooldown_iters) ? net.cooldown_iters : 0
  q_hyst_pu      = hasfield(typeof(net), :q_hyst_pu) ? net.q_hyst_pu : 0.0

  nodes     = net.nodeVec
  Sbase_MVA = net.baseMVA

  busVec, slackNum     = getBusData(nodes, Sbase_MVA, flatStart)
  busTypeVec, slackIdx = getBusTypeVec(busVec)

  # --- Active-set bookkeeping (PV origin mask for guarded re-enable) ---
  nb = length(busVec)
  pv_orig_mask = falses(nb)
  @inbounds for k = 1:nb
    pv_orig_mask[k] = (busVec[k].type == Sparlectra.PV)
  end

  # Re-enable only if hysteresis or cooldown configured (otherwise keep behavior: no re-enable).
  allow_reenable = (q_hyst_pu > 0.0) || (cooldown_iters > 0)
  
  n_pv = count(b->b.type==Sparlectra.PV, busVec)
  n_pq = count(b->b.type==Sparlectra.PQ, busVec)
  n    = n_pq + n_pv

  # Vset: actual setpoints for PV, 1.0 (dummy) otherwise
  Vset = [(b.type==Sparlectra.PV) ? (isnothing(net.nodeVec[b.idx]._vm_pu) ? b.vm_pu : net.nodeVec[b.idx]._vm_pu) : 1.0 for b in busVec]

  adjBranch = adjacentBranches(Y, debug)
  qmin_pu, qmax_pu = getQLimits_pu(net)
  if verbose > 1
    nshow = min(30, length(qmin_pu))
    println("Q-limits preview (first $nshow buses):")
    for i = 1:nshow
      qmin = isfinite(qmin_pu[i]) ? qmin_pu[i] : -Inf
      qmax = isfinite(qmax_pu[i]) ? qmax_pu[i] : Inf
      @printf "  bus %d: qmin=%g  qmax=%g\n" i qmin qmax
    end
  end


  it = 0;
  erg = 1
  resetQLimitLog!(net)
  while it <= maxIte
    Δ = residuum_full_withPV(Y, busVec, Vset, n_pq, n_pv, (verbose>2))    
    # --- Active-set PV/Q-limits (shared helper in limits.jl) ----------------
    if it > 1
      changed, reenabled = active_set_q_limits!(
        net, it, nb;
        qmin_pu = qmin_pu,
        qmax_pu = qmax_pu,
        pv_orig_mask = pv_orig_mask,
        allow_reenable = allow_reenable,
        q_hyst_pu = q_hyst_pu,
        cooldown_iters = cooldown_iters,
        verbose = verbose,

        # Q required at bus (p.u.) from last residuum evaluation
        get_qreq_pu = bus -> busVec[bus]._qRes,

        # Current PV status in busVec
        is_pv = bus -> (busVec[bus].type == Sparlectra.PV),

        # Enforce PV->PQ with Q clamp
        make_pq! = (bus, qclamp_pu, side) -> begin
          busVec[bus].type = Sparlectra.PQ
          busVec[bus].qƩ   = qclamp_pu
          busVec[bus].isPVactive = false

          # Mirror clamp into Net/Node (MVar)
          net.nodeVec[bus]._qƩGen = qclamp_pu * net.baseMVA
        end,

        # Optional PQ->PV re-enable
        make_pv! = (bus) -> begin
          busVec[bus].type = Sparlectra.PV
          busVec[bus].isPVactive = true
        end,
      )

      # If the active-set changed bus types/specs, recompute residuum immediately
      if changed || reenabled
        Δ = residuum_full_withPV(Y, busVec, Vset, n_pq, n_pv, (verbose>2))
      end

    end

    # Convergence check
    nrm = norm(Δ)
    if verbose > 0
      @printf "2-norm = %.6e, tol = %.6e, ite = %d\n" nrm tolerance it
    end

    if nrm < tolerance
      (verbose>0) && @info "Convergence after $(it) iterations"
      erg = 0;
      break
    end

    # Solve J * Δx = Δ
    J = calcJacobian_withPVIdentity(Y, busVec, adjBranch, slackIdx, n_pq, n_pv; log = (verbose>=3), sparse = sparse)
    Δx = nothing
    try
      Δx = J \ Δ
      if verbose > 3
        println("max |Δ| = ", maximum(abs.(Δx)))
      end
    catch err
      println("Linear solve failed: ", err)
      erg = 2;
      break
    end

    # State update (PV buses get ΔV_rel from the identity row)
    kSlack = 0
    for i in eachindex(busVec)
      if busVec[i].type == Sparlectra.Slack
        kSlack += 1;
        continue
      end
      # Full index without elimination (only shift for slack)
      idx1 = 2*((i >= slackIdx) ? (i-1) : i) - 1
      idx2 = idx1 + 1

      if angle_limit
        busVec[i].va_rad += min(Δx[idx1], 0.7)
      else
        busVec[i].va_rad += Δx[idx1]
      end
      busVec[i].vm_pu += busVec[i].vm_pu * Δx[idx2]

      # Reproject to polar consistency
      Vc = busVec[i].vm_pu * exp(im*busVec[i].va_rad)
      busVec[i].vm_pu = abs(Vc)
      busVec[i].va_rad = angle(Vc)
    end

    it += 1
  end

  # Write results back to the network object
  isoNodes = net.isoNodes
  for bus in busVec
    idx = bus.idx
    if idx in isoNodes
      continue
    end
    idx += count(i -> i < idx, isoNodes)
    setVmVa!(node = nodes[idx], vm_pu = bus.vm_pu, va_deg = rad2deg(bus.va_rad))

    if bus.isPVactive == false
      nodes[idx]._qƩGen = bus._qRes * Sbase_MVA
      @debug "Bus $(bus.idx) hit Q limit; set as PQ bus."
    end

    if bus.type == Sparlectra.PQ
      setBusType!(net, bus.idx, "PQ")  # oder: nodes[idx]._nodeType = Sparlectra.PQ
    elseif bus.type == Sparlectra.PV
      setBusType!(net, bus.idx, "PV")
    elseif bus.type == Sparlectra.Slack
      setBusType!(net, bus.idx, "Slack")  # optional, falls nötig
    end

    if bus.type == Sparlectra.PV
      nodes[idx]._qƩGen = bus._qRes * Sbase_MVA
    elseif bus.type == Sparlectra.Slack
      nodes[idx]._pƩGen = bus._pRes * Sbase_MVA
      nodes[idx]._qƩGen = bus._qRes * Sbase_MVA
    end

    if bus.isPVactive == false
      nodes[idx]._qƩGen = bus._qRes * Sbase_MVA
    end
  end

  p_pu = sum(b -> b._pRes, busVec)
  q_pu = sum(b -> b._qRes, busVec)

  p_MW   = p_pu * net.baseMVA
  q_MVar = q_pu * net.baseMVA

  if verbose > 1
    @info "Set total bus power to p = $p_MW MW and q = $q_MVar MVar"
  end

  setTotalBusPower!(net = net, p = p_MW, q = q_MVar)
  updateShuntPowers!(net = net)
  return it, erg
end

# ------------------------------
# Convenience wrapper similar to runpf!, but for the full system
# ------------------------------
function runpf_full!(net::Net, maxIte::Int, tolerance::Float64 = 1e-6, verbose::Int = 0; opt_sparse::Bool = false, opt_flatstart::Bool = net.flatstart)
  printYBus = (length(net.nodeVec) < 20) && (verbose > 1)
  Y = createYBUS(net = net, sparse = opt_sparse, printYBUS = printYBus)
  if verbose > 1
    Yabs_max = maximum(abs.(Y))
    Ydiag_max = maximum(abs.(diag(Y)))
    Ydiag_imag_max = maximum(abs.(imag.(diag(Y))))
    @printf "Ybus summary: Yabs_max=%.6e, Ydiag_max=%.6e, Ydiag_imag_max=%.6e\n" Yabs_max Ydiag_max Ydiag_imag_max
  end
  return calcNewtonRaphson_withPVIdentity!(net, Y, maxIte; tolerance = tolerance, verbose = verbose, sparse = opt_sparse, flatStart = opt_flatstart, angle_limit = false, debug = false)  
end

