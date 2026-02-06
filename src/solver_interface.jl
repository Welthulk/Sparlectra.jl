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

# file: src/solver_interface.jl
# solver_interface.jl — external solver bridge
#
# Goal:
# - Export a canonical, solver-agnostic PF model (Ybus, specs, start state, bus types).
# - Allow external solvers (HELM, custom NR variants, etc.) to consume the model.
# - Apply a solver solution back into `Net` safely and consistently.

using LinearAlgebra
using SparseArrays

"""
    PFModel

Solver-agnostic power-flow model extracted from a `Net`.

Important indexing note:
- The model uses the same "active-bus ordering" as `BusData` returned by `getBusData(...)`:
  sorted by original bus index and excluding isolated buses.
- `Ybus` must be consistent with that ordering. We build `Ybus` using `createYBUS(...)`
  which internally compresses indices by removing isolated buses.

Fields:
- `Ybus`: bus admittance matrix (n×n), active buses only
- `baseMVA`: system base power
- `busIdx_net`: original bus indices (length n), increasing
- `busType`: :Slack/:PV/:PQ per bus (length n)
- `slack_idx`: slack index in PF-ordering (1..n)
- `Vset`: voltage magnitude setpoints for PV/Slack buses (length n; dummy 1.0 for PQ)
- `Sspec`: specified complex injections in per-unit (length n), P + jQ (net injection)
- `V0`: initial complex voltages in per-unit (length n)
- `qmin_pu`, `qmax_pu`: optional PV reactive limits in per-unit (length n, may be empty)
"""
Base.@kwdef struct PFModel{TY<:AbstractMatrix{ComplexF64}}
  Ybus::TY
  baseMVA::Float64
  busIdx_net::Vector{Int}
  busType::Vector{Symbol}
  slack_idx::Int
  Vset::Vector{Float64}
  Sspec::Vector{ComplexF64}
  V0::Vector{ComplexF64}
  qmin_pu::Vector{Float64} = Float64[]
  qmax_pu::Vector{Float64} = Float64[]
end

"""
    PFSolution

Common result object returned by any external solver.

Fields:
- `V`: final complex bus voltages in PF ordering (length n)
- `converged`: true if solver reached tolerance
- `iters`: iteration/order count (solver-defined)
- `residual_inf`: infinity-norm of canonical mismatch (if available)
- `meta`: optional solver-specific payload (NamedTuple/Dict/etc.)
"""
Base.@kwdef struct PFSolution
  V::Vector{ComplexF64}
  converged::Bool
  iters::Int
  residual_inf::Float64
  meta::Any = nothing
end


"""
    showPfModel(model; io=stdout, verbose=false)

Print a compact human-readable summary of the PFModel.
This is intentionally separate from Base.show(model) to avoid global display side effects.
"""
function showPfModel(model::PFModel; io::IO=stdout, verbose::Bool=false)
  n = length(model.busIdx_net)
  n_slack = count(==( :Slack), model.busType)
  n_pv    = count(==( :PV),    model.busType)
  n_pq    = count(==( :PQ),    model.busType)

  println(io, "==================== PFModel ====================")
  println(io, "Active buses:   ", n)
  println(io, "Slack (pf idx): ", model.slack_idx, "   (net bus idx: ", model.busIdx_net[model.slack_idx], ")")
  println(io, "Bus types:      PQ=", n_pq, "  PV=", n_pv, "  Slack=", n_slack)
  println(io, "Ybus:           ", size(model.Ybus,1), "×", size(model.Ybus,2),
              "   sparse=", isa(model.Ybus, SparseMatrixCSC))
  println(io, "Has Q-limits:   ", (!isempty(model.qmin_pu) && !isempty(model.qmax_pu)))

  if verbose
    println(io, "Net bus indices (active): ", model.busIdx_net)
  end
  return nothing
end

"""
    showPfSolution(sol; io=stdout)

Print a compact human-readable summary of a PFSolution.
"""
function showPfSolution(sol::PFSolution; io::IO=stdout)
  println(io, "==================== PFSolution =================")
  println(io, "Converged:      ", sol.converged)
  println(io, "Iters/Order:    ", sol.iters)
  println(io, "Residual ||F||∞ ", sol.residual_inf)
  if sol.meta !== nothing
    println(io, "Meta:           ", sol.meta)
  end
  return nothing
end



"""
    AbstractExternalSolver

External solvers should subtype this and implement `solvePf(solver, model; kwargs...)`.
"""
abstract type AbstractExternalSolver end

"""
    solvePf(solver, model; kwargs...) -> PFSolution

External solver hook. Must be implemented by the external solver package.

The returned voltage vector `sol.V` must be in PF ordering (`model.busIdx_net` order).
"""
function solvePf(::AbstractExternalSolver, ::PFModel; kwargs...)
  error("solvePf(::AbstractExternalSolver, ::PFModel) is not implemented for this solver.")
end

# -----------------------------------------------------------------------------
# Model builder
# -----------------------------------------------------------------------------

"""
    buildPfModel(net; opt_sparse=true, flatstart=net.flatstart, include_limits=true, verbose=0) -> PFModel

Build a canonical PF model from `net`:
- active buses are defined via `getBusData(net.nodeVec, net.baseMVA, flatstart)`,
  which excludes isolated buses and sorts by bus index. :contentReference[oaicite:5]{index=5}
- Ybus is built via `createYBUS(net=net, sparse=opt_sparse, ...)` and is consistent
  with the same isolated-bus compression. :contentReference[oaicite:6]{index=6}

`include_limits=true` attaches qmin/qmax in per-unit as used by your NR-with-limits code. :contentReference[oaicite:7]{index=7}
"""
function buildPfModel(
  net::Net;
  opt_sparse::Bool = true,
  flatstart::Bool = net.flatstart,
  include_limits::Bool = true,
  verbose::Int = 0,
)
  # 1) BusData provides canonical active-bus ordering (skip iso, sort by idx)
  busVec, slackNum = getBusData(net.nodeVec, net.baseMVA, flatstart)
  n = length(busVec)
  n == 0 && error("buildPfModel: empty busVec (no active buses).")

  # 2) Build Ybus (internally skips isoNodes and compresses indices)
  Ybus = createYBUS(net = net, sparse = opt_sparse, printYBUS = (verbose > 1))
  @assert size(Ybus, 1) == n "buildPfModel: Ybus size $(size(Ybus)) does not match busVec length $n. Check isoNodes / numbering assumptions."
  @assert size(Ybus, 2) == n

  # 3) Types + slack index (PF ordering)
  busTypeVec, slack_idx = getBusTypeVec(busVec)
  @assert slack_idx == slackNum "buildPfModel: slack index mismatch (getBusTypeVec=$slack_idx, getBusData=$slackNum)."

  busType = Vector{Symbol}(undef, n)
  busIdx_net = Vector{Int}(undef, n)
  V0 = Vector{ComplexF64}(undef, n)
  Sspec = Vector{ComplexF64}(undef, n)
  Vset = ones(Float64, n)

  @inbounds for k in 1:n
    b = busVec[k]
    busIdx_net[k] = b.idx

    # NodeType -> Symbol
    if b.type == Slack
      busType[k] = :Slack
    elseif b.type == PV
      busType[k] = :PV
    elseif b.type == PQ
      busType[k] = :PQ
    else
      error("buildPfModel: unsupported BusData type $(b.type) at bus $(b.idx)")
    end

    # Initial V0 from BusData (already respects flatstart)
    V0[k] = ComplexF64(b.vm_pu * cos(b.va_rad), b.vm_pu * sin(b.va_rad))

    # Specified injections in p.u. from BusData
    Sspec[k] = ComplexF64(b.pƩ, b.qƩ)

    # PV/Slack voltage magnitude setpoint (best effort: node._vm_pu if present, else BusData vm)
    if b.type == PV || b.type == Slack
      vm_node = net.nodeVec[b.idx]._vm_pu
      Vset[k] = (vm_node === nothing) ? b.vm_pu : Float64(vm_node)
    else
      Vset[k] = 1.0
    end
  end

  # 4) Optional limits aligned with busVec indexing (as in jacobian_full.jl) :contentReference[oaicite:8]{index=8}
  qmin_pu = Float64[]
  qmax_pu = Float64[]
  if include_limits
    qmin_pu, qmax_pu = getQLimits_pu(net)
    if !(length(qmin_pu) == n && length(qmax_pu) == n)
      # Fallback: do not fail hard; limits are optional.
      qmin_pu = Float64[]
      qmax_pu = Float64[]
      (verbose > 0) && @warn "buildPfModel: q-limit vectors size mismatch; dropping q-limits for PFModel."
    end
  end

  return PFModel(
    Ybus      = Ybus,
    baseMVA   = net.baseMVA,
    busIdx_net = busIdx_net,
    busType   = busType,
    slack_idx = slack_idx,
    Vset      = Vset,
    Sspec     = Sspec,
    V0        = V0,
    qmin_pu   = qmin_pu,
    qmax_pu   = qmax_pu,
  )
end

# -----------------------------------------------------------------------------
# Canonical mismatch / residual
# -----------------------------------------------------------------------------

"""
    mismatchInf(model, V) -> Float64

Compute `||F(V)||_∞` using Sparlectra's rectangular PQ/PV mismatch definition:

    F = mismatch_rectangular(model.Ybus, V, model.Sspec, model.busType, model.Vset, model.slack_idx)

This is the canonical comparison metric for external solvers. :contentReference[oaicite:9]{index=9}
"""
function mismatchInf(model::PFModel, V::Vector{ComplexF64})::Float64
  F = mismatch_rectangular(model.Ybus, V, model.Sspec, model.busType, model.Vset, model.slack_idx)
  return maximum(abs.(F))
end

# -----------------------------------------------------------------------------
# Apply solution back into Net
# -----------------------------------------------------------------------------

"""
    applyPfSolution!(net, model, sol; write_pq_results=true, verbose=0)

Write the external solution back into `net`:
- updates node Vm/Va on all active buses
- optionally writes Slack P/Q and PV Q generation from solved injections (same policy
  as your rectangular NR integration) :contentReference[oaicite:10]{index=10}

Notes:
- Isolated buses (not in model) are left untouched.
- `sol.V` must be PF-ordered (length n == length(model.busIdx_net)).
"""
function applyPfSolution!(
  net::Net,
  model::PFModel,
  sol::PFSolution;
  write_pq_results::Bool = true,
  verbose::Int = 0,
)
  V = sol.V
  n = length(model.busIdx_net)
  @assert length(V) == n "applyPfSolution!: length(sol.V) != model size."

  # 1) Update Vm/Va on active buses
  @inbounds for k in 1:n
    busIdx = model.busIdx_net[k]
    node = net.nodeVec[busIdx]
    Vk = V[k]
    node._vm_pu = abs(Vk)
    node._va_deg = rad2deg(angle(Vk))
  end

  if !write_pq_results
    return nothing
  end

  # 2) Compute nodal injections from solved state (per-unit and MVA)
  Ibus = model.Ybus * V
  Sbus_pu = V .* conj.(Ibus)
  Sbus_MVA = Sbus_pu .* model.baseMVA

  # 3) Write back Slack P/Q and PV Q generation (policy as in jacobian_complex.jl) :contentReference[oaicite:11]{index=11}
  @inbounds for k in 1:n
    busIdx = model.busIdx_net[k]
    node = net.nodeVec[busIdx]
    S_MVA = Sbus_MVA[k]
    Pbus = real(S_MVA)
    Qbus = imag(S_MVA)

    if getNodeType(node) == Slack
      node._pƩGen = Pbus
      node._qƩGen = Qbus
    elseif getNodeType(node) == PV
      node._qƩGen = Qbus
    end
  end

  # 4) Update total bus power (sum of complex injections in p.u., converted to system units)
  p_total = sum(real.(Sbus_pu)) * model.baseMVA
  q_total = sum(imag.(Sbus_pu)) * model.baseMVA
  setTotalBusPower!(net = net, p = p_total, q = q_total)

  (verbose > 1) && @info "applyPfSolution!: wrote Vm/Va and Slack/PV injections back to Net."

  return nothing
end

# -----------------------------------------------------------------------------
# Convenience runner for external solvers
# -----------------------------------------------------------------------------

"""
    runpf_external!(net, solver; tol=1e-8, opt_sparse=true, flatstart=net.flatstart, include_limits=false, verbose=0, solver_kwargs...) -> (iters, status, sol)

Convenience wrapper:
1) build model from net
2) call external `solvePf(solver, model; solver_kwargs...)`
3) compute canonical mismatch infinity norm (unless solver provided it)
4) apply solution back to net
5) return `(iters, status, sol)` where `status == 0` indicates convergence

`include_limits=false` by default: keep Sparlectra "clean" and let external solver decide. :contentReference[oaicite:12]{index=12}
"""
function runpf_external!(
  net::Net,
  solver::AbstractExternalSolver;
  tol::Float64 = 1e-8,
  opt_sparse::Bool = true,
  flatstart::Bool = net.flatstart,
  include_limits::Bool = false,
  verbose::Int = 0,
  show_model::Bool = false,
  show_solution::Bool = false,
  io::IO = stdout,
  solver_kwargs...,
)
  model = buildPfModel(net; opt_sparse=opt_sparse, flatstart=flatstart, include_limits=include_limits, verbose=verbose)

  show_model && showPfModel(model; io=io, verbose=(verbose > 0))

  sol = solvePf(solver, model; tol=tol, solver_kwargs...)

  res = (sol.residual_inf <= 0.0 || !isfinite(sol.residual_inf)) ? mismatchInf(model, sol.V) : sol.residual_inf

  sol2 = PFSolution(
    V = sol.V,
    converged = (sol.converged && (res <= tol)),
    iters = sol.iters,
    residual_inf = res,
    meta = sol.meta,
  )

  show_solution && showPfSolution(sol2; io=io)

  applyPfSolution!(net, model, sol2; write_pq_results=true, verbose=verbose)

  status = sol2.converged ? 0 : 1
  return sol2.iters, status, sol2
end
