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

# file: src/upfc_full_control.jl
# purpose: the FULL UPFC (issue #326): one controller owning the series
#          converter (arbitrary-phase voltage injection, steering line P AND
#          Q independently) and the shunt converter, coupled through the
#          DC-link active-power balance P_se + P_sh = 0.
#
#   Electrical model:
#
#      bus i (from)      series converter          bus j (to)
#        V_i  o----+---[ + V_se - ]---[ z_line ]---o  V_j
#                  |                                    line i->j, current I_s
#            [ shunt conv ]  injects P_sh = -P_se (DC balance) + Q_sh
#                  |
#                 === DC link ===  (couples the two converters)
#
#   The series converter injects a voltage V_se of ARBITRARY phase; it is
#   realised as an equivalent series impedance z_add = V_se / I_s added to the
#   branch (Re(z_add) < 0 when the converter injects active power), reusing the
#   SSSC branch re-stamp mechanism, so no injection carrier is needed. Its
#   active power P_se = Re(V_se·conj(I_s)) = |I_s|^2·Re(z_add) flows through the
#   DC link and is supplied by the shunt converter at bus i (P_sh = -P_se). The
#   phase-shifter degree of freedom (independent P AND Q on one line) is exactly
#   what the #325 quadrature composite lacks: there V_se is orthogonal to I_s,
#   P_se = 0, the DC link idles, and the device collapses to SSSC + STATCOM.
#
#   First-cut scope: the shunt converter provides the DC-link balance plus a
#   reactive SETPOINT (q_shunt_mvar), inside the current-based rating whose
#   reactive headroom is coupled to P_sh by Q_max = sqrt((V_i·s_max)^2 - P_sh^2).
#   Closed-loop shunt VOLTAGE regulation is deferred: coupling a shunt-voltage
#   secant with the line reactive-flow control does not converge in the
#   sequential outer loop (a known behaviour of injection-model UPFCs), and a
#   robust version needs the AC power-flow sensitivity framework (issue #217) or
#   an augmented in-solver state (which would touch the protected solver).

"""
    UpfcFullControl <: AbstractOuterController

Full unified power flow controller (issue #326), the DC-link-coupled model.

The series converter injects a voltage `V_se` of ARBITRARY phase into the line
`fromBus`->`toBus`, steering the from-end line flow to `(p_target_mw,
q_target_mvar)` INDEPENDENTLY (the phase-shifter degree of freedom the #325
quadrature composite lacks). The active power the series converter exchanges
with the line, `P_se = Re(V_se·conj(I_s))`, flows through the DC link and is
supplied by the shunt converter (`P_sh = -P_se`) at `shunt_bus`; the shunt also
delivers the reactive setpoint `q_shunt_mvar`, clamped to the current-based
rating whose reactive headroom is coupled to `P_sh` by
`Q_max = sqrt((V_shunt·s_max)^2 - P_sh^2)`.

`series_phase = :quadrature` constrains `V_se ⟂ I_s` (`P_se = 0`), reducing the
series converter to the SSSC (the regression bridge to #325).

Runtime fields are owned by `run_control!`; construct via
[`addUpfcControl!`](@ref) with `model = :full`.
"""
mutable struct UpfcFullControl <: AbstractOuterController
  name::String
  fromBus::String
  toBus::String
  branch_idx::Int
  shunt_bus::String            # bus carrying the shunt converter (an end of the corridor)
  shunt_prosumer_idx::Int      # the shunt converter's generator prosumer
  p_target_mw::Float64         # from-end line active-power target
  q_target_mvar::Float64       # from-end line reactive-power target (the new DOF)
  q_shunt_mvar::Float64        # shunt reactive setpoint
  v_inj_max_pu::Float64        # series injected-voltage limit
  s_max_mva::Float64           # shunt converter rating at 1.0 pu
  r_base_pu::Float64           # natural branch resistance at registration
  x_base_pu::Float64           # natural branch reactance at registration
  series_phase::Symbol         # :free (full) or :quadrature (SSSC reduction)
  deadband_p_mw::Float64
  deadband_q_mvar::Float64
  max_outer_iters::Int
  enabled::Bool
  # runtime state, owned by the control loop
  status::Symbol
  converged::Bool
  at_limit::Bool
  v_se_pu::ComplexF64          # current series injected voltage
  i_s_pu::ComplexF64           # last measured series current
  p_se_mw::Float64             # series converter active power (into the line)
  q_sh_mvar::Float64           # shunt reactive injection (after the coupled clamp)
  p_sh_mw::Float64             # shunt active injection (= -p_se)
  qmin_mvar::Float64           # live coupled shunt bound
  qmax_mvar::Float64
  achieved_p_mw::Union{Nothing,Float64}
  achieved_q_mvar::Union{Nothing,Float64}
  shunt_vm_pu::Union{Nothing,Float64}     # shunt-bus voltage (reported, not regulated)
  dc_residual_mw::Union{Nothing,Float64}  # |P_se + P_sh| on the re-solved state
  prev_x_pu::Union{Nothing,Float64}       # quadrature secant memory (reactance path)
  prev_p_mw::Union{Nothing,Float64}
  damp::Float64                           # adaptive series step damping (globalisation)
  err_prev::Union{Nothing,Float64}        # previous flow-error norm, for the line search
  upfc_group::Union{Nothing,String}
  outer_iters::Int
end

"""
    _upfc_full_controllers(net) -> Vector{UpfcFullControl}

Collect the full-UPFC controllers stored on `net` (shared registry
`net.machineControls`).
"""
function _upfc_full_controllers(net)::Vector{UpfcFullControl}
  hasproperty(net, :machineControls) || return UpfcFullControl[]
  return UpfcFullControl[c for c in net.machineControls if c isa UpfcFullControl]
end

# complex bus voltage in pu from the solved state
function _upfc_vc(net::Net, busName::String)::ComplexF64
  idx = geNetBusIdx(net = net, busName = busName)
  node = net.nodeVec[idx]
  vm = node._vm_pu === nothing ? 1.0 : node._vm_pu
  va = node._va_deg === nothing ? 0.0 : node._va_deg
  return vm * cis(deg2rad(va))
end

# series admittance from the natural (base) branch impedance
_upfc_y_base(ctrl::UpfcFullControl)::ComplexF64 = 1.0 / complex(ctrl.r_base_pu, ctrl.x_base_pu)

# current-floor so z_add stays finite on a (numerically) currentless branch
const _UPFC_MIN_CURRENT_PU = 1e-4
# series impedance-magnitude guard, mirrors the SSSC resonance exclusion
const _UPFC_EPS_Z = 1e-4
# quadrature reactance-secant min slope / bootstrap nudge
const _UPFC_MIN_SLOPE = 1e-6
const _UPFC_X_BOOTSTRAP = 0.02
# damping of the free-mode series step: the 2x2 solve is EXACT at the frozen
# terminal voltages, but the shunt's DC-balance injection (P_sh = -P_se) moves
# those voltages each re-solve, so a full step can overshoot when the injection
# is large. The step is ADAPTIVELY damped (a fixed-point line search): the
# factor shrinks when the flow error grows and relaxes back when it falls, so
# the coupled iteration is globalised without changing its fixed point.
const _UPFC_SERIES_DAMP = 0.6      # nominal / starting damping
const _UPFC_DAMP_MIN = 0.1         # floor when the iteration is stiff
const _UPFC_DAMP_SHRINK = 0.5      # cut on a growing error
const _UPFC_DAMP_GROW = 1.3        # relax on a shrinking error

"""
    addUpfcFullControl!(net; ...) -> UpfcFullControl

Internal constructor for the full UPFC. Users call [`addUpfcControl!`](@ref)
with `model = :full`, which validates the composite preconditions and forwards
here. Registers one `UpfcFullControl` on `net.machineControls` and returns it.
The shunt converter is the existing generator at `shunt_bus` (resolved and
required, as in the #325 composite); no injection carrier is created (the
series source is realised as an equivalent branch impedance).
"""
function addUpfcFullControl!(
  net::Net;
  fromBus::String,
  toBus::String,
  shunt_bus::String,
  p_target_mw::Float64,
  q_target_mvar::Float64,
  v_inj_max_pu::Float64,
  q_shunt_mvar::Float64 = 0.0,
  s_max_mva::Union{Nothing,Float64} = nothing,
  i_max_ka::Union{Nothing,Float64} = nothing,
  series_phase::Symbol = :free,
  deadband_p_mw::Float64 = 0.5,
  deadband_q_mvar::Float64 = 0.5,
  prosumer_index::Union{Nothing,Int} = nothing,
  name::Union{Nothing,String} = nothing,
  max_outer_iters::Int = 40,
  enabled::Bool = true,
)
  series_phase in (:free, :quadrature) || error("UpfcControl(full): series_phase must be :free or :quadrature, got $(series_phase).")
  # branch existence and line-kind, mirroring the series controller
  br = getNetBranch(net = net, fromBus = fromBus, toBus = toBus)
  if br === nothing
    rev = getNetBranch(net = net, fromBus = toBus, toBus = fromBus)
    rev === nothing && error("UpfcControl(full): no branch between $(fromBus) and $(toBus).")
    error("UpfcControl(full): the branch is stored as $(toBus) to $(fromBus); register with that bus order (the pair fixes the flow measurement direction).")
  end
  br.ratio == 0.0 || error("UpfcControl(full): branch $(fromBus) to $(toBus) is a transformer branch; attach the UPFC to a line branch.")
  any(c -> c isa UpfcFullControl && c.branch_idx == br.branchIdx, net.machineControls) && error("UpfcControl(full): branch $(fromBus) to $(toBus) already has a full UPFC.")
  any(c -> c isa SeriesReactanceControl && c.branch_idx == br.branchIdx, net.machineControls) && error("UpfcControl(full): branch $(fromBus) to $(toBus) already has a series controller.")

  # shunt converter rating (STATCOM current limit), resolved as in addMachineVoltageControl!
  (s_max_mva === nothing) == (i_max_ka === nothing) && error("UpfcControl(full): pass exactly one shunt converter rating, s_max_mva or i_max_ka.")
  shuntIdx = geNetBusIdx(net = net, busName = shunt_bus)
  s_rating = s_max_mva
  if i_max_ka !== nothing
    i_max_ka > 0.0 || error("UpfcControl(full): i_max_ka must be positive, got $(i_max_ka)")
    s_rating = sqrt(3.0) * getNodeVn(net.nodeVec[shuntIdx]) * i_max_ka
  end
  s_rating > 0.0 || error("UpfcControl(full): s_max_mva must be positive, got $(s_rating)")
  v_inj_max_pu > 0.0 || error("UpfcControl(full): v_inj_max_pu must be positive, got $(v_inj_max_pu)")

  # resolve the shunt converter generator (an existing PQ generator at shunt_bus)
  ps_idx = if prosumer_index !== nothing
    (1 <= prosumer_index <= length(net.prosumpsVec)) || error("UpfcControl(full): prosumer_index $(prosumer_index) out of range")
    isGenerator(net.prosumpsVec[prosumer_index]) || error("UpfcControl(full): prosumer $(prosumer_index) is not a generator")
    getPosumerBusIndex(net.prosumpsVec[prosumer_index]) == shuntIdx || error("UpfcControl(full): prosumer $(prosumer_index) is not at bus $(shunt_bus)")
    prosumer_index
  else
    gens = [i for (i, p) in enumerate(net.prosumpsVec) if isGenerator(p) && getPosumerBusIndex(p) == shuntIdx]
    isempty(gens) && error("UpfcControl(full): no generator at shunt bus $(shunt_bus) to carry the shunt converter.")
    length(gens) > 1 && error("UpfcControl(full): $(length(gens)) generators at $(shunt_bus) - pass prosumer_index to pick the shunt converter.")
    gens[1]
  end
  net.prosumpsVec[ps_idx].isRegulated && error("UpfcControl(full): the generator at $(shunt_bus) is voltage-regulating (PV); the shunt converter needs a PQ machine.")

  cname = something(name, string("UPFC_", fromBus, "_", toBus))
  ctrl = UpfcFullControl(
    cname, fromBus, toBus, br.branchIdx, shunt_bus, ps_idx,
    p_target_mw, q_target_mvar, q_shunt_mvar, v_inj_max_pu, s_rating,
    br.r_pu, br.x_pu, series_phase,
    deadband_p_mw, deadband_q_mvar, max_outer_iters, enabled,
    :active, false, false,
    complex(0.0, 0.0), complex(0.0, 0.0), 0.0, 0.0, 0.0,
    -s_rating, s_rating,
    nothing, nothing, nothing, nothing,
    nothing, nothing,
    _UPFC_SERIES_DAMP, nothing,
    nothing, 0,
  )
  push!(net.machineControls, ctrl)
  return ctrl
end

# resolve the controlled branch defensively (as the series controller does)
function _upfc_branch(net::Net, ctrl::UpfcFullControl)::Branch
  (1 <= ctrl.branch_idx <= length(net.branchVec)) || error("UpfcControl(full) $(ctrl.name): branch index $(ctrl.branch_idx) out of range")
  return net.branchVec[ctrl.branch_idx]
end

function _upfc_shunt_prosumer(net::Net, ctrl::UpfcFullControl)::ProSumer
  (1 <= ctrl.shunt_prosumer_idx <= length(net.prosumpsVec)) || error("UpfcControl(full) $(ctrl.name): shunt prosumer index out of range")
  return net.prosumpsVec[ctrl.shunt_prosumer_idx]
end

# --- AbstractOuterController protocol ---------------------------------------

control_name(ctrl::UpfcFullControl) = string(ctrl.name, " UPFC")
control_enabled(ctrl::UpfcFullControl) = ctrl.enabled
control_initialize!(::UpfcFullControl, ::Net, context) = NoControlState()
control_status(ctrl::UpfcFullControl, ::AbstractControlState)::Symbol = ctrl.status
control_is_converged(ctrl::UpfcFullControl, ::AbstractControlState)::Bool = ctrl.converged
control_is_blocked(ctrl::UpfcFullControl, ::AbstractControlState)::Bool = ctrl.at_limit

function control_evaluate!(ctrl::UpfcFullControl, net::Net, ::AbstractControlState, context)
  # measured from-end line flow (reflects the modified branch after re-solve)
  p = get_branch_p_from_to_mw(net, ctrl.fromBus, ctrl.toBus)
  q = get_branch_q_from_to_mvar(net, ctrl.fromBus, ctrl.toBus)
  ctrl.achieved_p_mw = p
  ctrl.achieved_q_mvar = q
  ctrl.shunt_vm_pu = get_bus_vm_pu(net, ctrl.shunt_bus)
  # series current on the modified branch, from the solved terminal voltages
  Vi = _upfc_vc(net, ctrl.fromBus)
  Vj = _upfc_vc(net, ctrl.toBus)
  br = _upfc_branch(net, ctrl)
  y_eff = 1.0 / complex(br.r_pu, br.x_pu)
  ctrl.i_s_pu = y_eff * (Vi - Vj)
  # DC-link residual on the RE-SOLVED state: the balance holds by construction
  # only at the frozen state each iteration; after the re-solve the voltages
  # move, so this is a genuine convergence quantity
  p_se_now = ctrl.series_phase === :quadrature ? 0.0 : real(ctrl.v_se_pu * conj(ctrl.i_s_pu)) * net.baseMVA
  ctrl.dc_residual_mw = abs(p_se_now + ctrl.p_sh_mw)
  # convergence: from-end flow targets (P always; Q only in :free)
  p_ok = abs(p - ctrl.p_target_mw) <= ctrl.deadband_p_mw
  q_ok = ctrl.series_phase === :quadrature ? true : abs(q - ctrl.q_target_mvar) <= ctrl.deadband_q_mvar
  ctrl.converged = p_ok && q_ok
  ctrl.status = ctrl.converged ? :converged : (ctrl.at_limit ? :at_limit : :active)
  ctrl.outer_iters = context.outer_iteration
  return nothing
end

function control_propose_update!(ctrl::UpfcFullControl, net::Net, ::AbstractControlState, context)
  br = _upfc_branch(net, ctrl)
  ctrl.converged && return (r_new = br.r_pu, x_new = br.x_pu, v_se = ctrl.v_se_pu, p_sh = ctrl.p_sh_mw, q_sh = ctrl.q_sh_mvar, p_se = ctrl.p_se_mw)
  base = net.baseMVA
  Vi = _upfc_vc(net, ctrl.fromBus)
  Vj = _upfc_vc(net, ctrl.toBus)
  y_s = _upfc_y_base(ctrl)
  ΔV = Vi - Vj

  local v_se::ComplexF64
  local r_new::Float64
  local x_new::Float64
  if ctrl.series_phase === :quadrature
    # SSSC reduction: 1-DOF reactance secant on the P error (Q ignored), z_add
    # pure imaginary so P_se stays zero
    p_now = ctrl.achieved_p_mw === nothing ? 0.0 : ctrl.achieved_p_mw
    err_p = ctrl.p_target_mw - p_now
    x_now = br.x_pu
    Δx = if ctrl.prev_x_pu !== nothing && ctrl.prev_p_mw !== nothing && abs(x_now - ctrl.prev_x_pu) > 1e-12
      slope = (p_now - ctrl.prev_p_mw) / (x_now - ctrl.prev_x_pu)
      abs(slope) > _UPFC_MIN_SLOPE ? err_p / slope : -sign(err_p) * _UPFC_X_BOOTSTRAP
    else
      # bootstrap: raising reactance lowers the through-flow
      -sign(err_p) * _UPFC_X_BOOTSTRAP
    end
    x_target = x_now + Δx
    is_mag = max(abs(ctrl.i_s_pu), _UPFC_MIN_CURRENT_PU)
    dx_max = ctrl.v_inj_max_pu / is_mag
    Δx_total = clamp(x_target - ctrl.x_base_pu, -dx_max, dx_max)
    x_new = ctrl.x_base_pu + Δx_total
    r_new = ctrl.r_base_pu
    v_se = complex(0.0, Δx_total) * ctrl.i_s_pu   # quadrature by construction
    ctrl.prev_x_pu = x_now
    ctrl.prev_p_mw = p_now
  else
    # full model: with frozen terminal voltages S_from = V_i·conj(I_s),
    # I_s = y_s(ΔV - V_se), so S_from is affine in conj(V_se) with sensitivity
    # dS_from/dconj(V_se) = -V_i·conj(y_s). Solving for V_se to hit the target
    # is therefore an EXACT 2x2 linear solve at the frozen voltages; the outer
    # re-solve chases the true fixed point.
    p_now = (ctrl.achieved_p_mw === nothing ? 0.0 : ctrl.achieved_p_mw) / base
    q_now = (ctrl.achieved_q_mvar === nothing ? 0.0 : ctrl.achieved_q_mvar) / base
    ΔS = complex(ctrl.p_target_mw / base - p_now, ctrl.q_target_mvar / base - q_now)
    # adaptive damping (fixed-point line search): shrink the step when the flow
    # error grew since the last iteration, relax it back when the error fell
    err = abs(ΔS) * base
    if ctrl.err_prev !== nothing
      ctrl.damp = err > ctrl.err_prev ? max(_UPFC_DAMP_MIN, ctrl.damp * _UPFC_DAMP_SHRINK) : min(_UPFC_SERIES_DAMP, ctrl.damp * _UPFC_DAMP_GROW)
    end
    ctrl.err_prev = err
    c1 = Vi * conj(y_s)
    ΔV_se = abs(c1) > 1e-12 ? -conj(ΔS / c1) : complex(0.0, 0.0)
    v_se = ctrl.v_se_pu + ctrl.damp * ΔV_se
    if abs(v_se) > ctrl.v_inj_max_pu
      v_se = v_se * (ctrl.v_inj_max_pu / abs(v_se))
    end
    # realise as an equivalent series impedance z_add = V_se / I_s
    i_s = y_s * (ΔV - v_se)
    z_add = abs(i_s) > _UPFC_MIN_CURRENT_PU ? v_se / i_s : complex(0.0, 0.0)
    r_cand = ctrl.r_base_pu + real(z_add)
    x_cand = ctrl.x_base_pu + imag(z_add)
    # resonance guard: keep the series impedance magnitude out of the eps_z
    # neighbourhood (never error mid-run, clamp instead)
    if abs(complex(r_cand, x_cand)) < _UPFC_EPS_Z
      r_new = ctrl.r_base_pu
      x_new = ctrl.x_base_pu
      v_se = complex(0.0, 0.0)
    else
      r_new = r_cand
      x_new = x_cand
    end
  end

  # series converter active power P_se = Re(V_se·conj(I_s)) at the proposed
  # impedance (0 by construction in quadrature)
  y_new = 1.0 / complex(r_new, x_new)
  i_s_new = y_new * ΔV
  p_se_mw = ctrl.series_phase === :quadrature ? 0.0 : real(v_se * conj(i_s_new)) * base

  # shunt: DC-link balance P_sh = -P_se, reactive setpoint clamped to the
  # coupled bound Q_max = sqrt((V_shunt·s_max)^2 - P_sh^2)
  p_sh_mw = -p_se_mw
  vsh = get_bus_vm_pu(net, ctrl.shunt_bus)
  s_lim = (isfinite(vsh) && vsh > 0.0 ? vsh : 1.0) * ctrl.s_max_mva
  q_head = s_lim^2 - p_sh_mw^2
  qmax = q_head > 0.0 ? sqrt(q_head) : 0.0
  q_sh = clamp(ctrl.q_shunt_mvar, -qmax, qmax)

  return (r_new = r_new, x_new = x_new, v_se = v_se, p_sh = p_sh_mw, q_sh = q_sh, p_se = p_se_mw)
end

function control_apply_update!(ctrl::UpfcFullControl, net::Net, ::AbstractControlState, update::NamedTuple, context)::Bool
  br = _upfc_branch(net, ctrl)
  ps = _upfc_shunt_prosumer(net, ctrl)
  node = net.nodeVec[geNetBusIdx(net = net, busName = ctrl.shunt_bus)]

  moved_branch = abs(update.r_new - br.r_pu) > 1e-12 || abs(update.x_new - br.x_pu) > 1e-12
  br.r_pu = update.r_new
  br.x_pu = update.x_new

  # shunt injection: set the generator P/Q, keeping the node sums coherent
  # (the solve reads ps.pVal/qVal, losses/results read the node sums)
  p_old = ps.pVal === nothing ? 0.0 : ps.pVal
  q_old = ps.qVal === nothing ? 0.0 : ps.qVal
  Δp = update.p_sh - p_old
  Δq = update.q_sh - q_old
  ps.pVal = update.p_sh
  ps.qVal = update.q_sh
  addGenPower!(node = node, p = Δp, q = Δq)
  moved_shunt = abs(Δp) > 1e-9 || abs(Δq) > 1e-9

  ctrl.v_se_pu = update.v_se
  ctrl.p_se_mw = update.p_se
  ctrl.p_sh_mw = update.p_sh
  ctrl.q_sh_mvar = update.q_sh

  # honest at_limit: series voltage clamp OR shunt reactive clamp (the reactive
  # setpoint exceeds the coupled bound), with the flow target still missed
  series_clamped = abs(update.v_se) >= ctrl.v_inj_max_pu * (1 - 1e-9)
  shunt_clamped = abs(ctrl.q_shunt_mvar) > abs(update.q_sh) + 1e-9
  ctrl.at_limit = !ctrl.converged && (series_clamped || shunt_clamped)
  ctrl.at_limit && (ctrl.status = :at_limit)
  return moved_branch || moved_shunt
end

function control_element_descriptor(ctrl::UpfcFullControl, net::Net)::Union{Nothing,NamedTuple}
  return (
    name = control_name(ctrl),
    element = string("branch@", ctrl.fromBus, "-", ctrl.toBus),
    device = "UPFC (full, DC-link coupled)",
    actuator = :upfc_series_voltage,
    actuator_min = 0.0,
    actuator_max = ctrl.v_inj_max_pu,
    quantity = :branch_pq,
    target = string(ctrl.fromBus, "->", ctrl.toBus),
    target_value = ctrl.p_target_mw,
    discrete = false,
    enabled = ctrl.enabled,
    status = ctrl.status,
    converged = ctrl.converged,
    at_limit = ctrl.at_limit,
  )
end

function control_report_rows(ctrl::UpfcFullControl, net::Net, ::AbstractControlState, context)
  return [(
    controller_name = control_name(ctrl),
    controller_type = "UpfcFullControl",
    branch = string(ctrl.fromBus, "->", ctrl.toBus),
    status = ctrl.status,
    converged = ctrl.converged,
    at_limit = ctrl.at_limit,
    achieved_p_mw = ctrl.achieved_p_mw === nothing ? missing : ctrl.achieved_p_mw,
    p_target_mw = ctrl.p_target_mw,
    achieved_q_mvar = ctrl.achieved_q_mvar === nothing ? missing : ctrl.achieved_q_mvar,
    q_target_mvar = ctrl.q_target_mvar,
    v_se_mag_pu = abs(ctrl.v_se_pu),
    v_se_ang_deg = rad2deg(angle(ctrl.v_se_pu)),
    p_se_mw = ctrl.p_se_mw,
    p_sh_mw = ctrl.p_sh_mw,
    q_sh_mvar = ctrl.q_sh_mvar,
    shunt_vm_pu = ctrl.shunt_vm_pu === nothing ? missing : ctrl.shunt_vm_pu,
    dc_residual_mw = ctrl.dc_residual_mw === nothing ? missing : ctrl.dc_residual_mw,
    series_phase = ctrl.series_phase,
  )]
end

control_trace_rows(::UpfcFullControl, ::Net, ::AbstractControlState, context) = NamedTuple[]
