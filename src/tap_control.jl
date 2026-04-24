# Copyright 2023–2026 Udo Schmitz

"""
    TapController

Outer-loop transformer controller definition.

# Fields
- `trafo`: Transformer identifier (branch component name/id or branch index as `String`).
- `mode`: `:voltage`, `:branch_active_power`, or `:voltage_and_branch_active_power`.
- `target_bus`: Controlled bus name for voltage mode.
- `target_branch`: Controlled branch tuple `(fromBus, toBus)` for active-power mode.
- `target_vm_pu`: Voltage target in p.u.
- `p_target_mw`: Active-power target in MW.
- `q_target_mvar`: Optional reactive-power target (reserved for future use).
- `control_ratio`: Enable ratio (`τ`) updates.
- `control_phase`: Enable phase (`φ`) updates.
- `is_discrete`: Use discrete tap-step updates.
- `deadband_vm_pu`: Voltage deadband.
- `deadband_p_mw`: Active-power deadband.
- `max_outer_iters`: Maximum outer-loop iterations.
- `enabled`: Enable/disable controller.
- Runtime fields (`converged`, `at_limit`, `status`, `achieved_*`, `outer_iters`) are
  updated by `run_tap_controllers_outer!`.
"""
mutable struct TapController
  trafo::String
  mode::Symbol
  target_bus::Union{Nothing,String}
  target_branch::Union{Nothing,Tuple{String,String}}
  target_vm_pu::Union{Nothing,Float64}
  p_target_mw::Union{Nothing,Float64}
  q_target_mvar::Union{Nothing,Float64}
  control_ratio::Bool
  control_phase::Bool
  is_discrete::Bool
  deadband_vm_pu::Float64
  deadband_p_mw::Float64
  max_outer_iters::Int
  enabled::Bool
  converged::Bool
  at_limit::Bool
  status::Symbol
  achieved_vm_pu::Union{Nothing,Float64}
  achieved_p_mw::Union{Nothing,Float64}
  outer_iters::Int
end

"""
    TapController(; ...)

Convenience constructor for `TapController` with sane defaults for outer-loop
execution and status initialization.
"""
function TapController(; trafo::String, mode::Symbol, target_bus::Union{Nothing,String}=nothing, target_branch::Union{Nothing,Tuple{String,String}}=nothing,
  target_vm_pu::Union{Nothing,Float64}=nothing, p_target_mw::Union{Nothing,Float64}=nothing, q_target_mvar::Union{Nothing,Float64}=nothing,
  control_ratio::Bool=true, control_phase::Bool=false, is_discrete::Bool=true, deadband_vm_pu::Float64=1e-3, deadband_p_mw::Float64=0.5,
  max_outer_iters::Int=20, enabled::Bool=true)
  return TapController(trafo, mode, target_bus, target_branch, target_vm_pu, p_target_mw, q_target_mvar, control_ratio, control_phase, is_discrete,
    deadband_vm_pu, deadband_p_mw, max_outer_iters, enabled, false, false, :idle, nothing, nothing, 0)
end

"""
    _find_trafo_branch(net, name)

Internal helper to resolve a transformer branch by component name, component id,
or branch index represented as string.
"""
function _find_trafo_branch(net::Net, name::String)::Branch
  for br in net.branchVec
    if br.ratio != 0.0 && (br.comp.cName == name || br.comp.cID == name || string(br.branchIdx) == name)
      return br
    end
  end
  error("TapController: transformer $(name) not found. Use branch component name/id or branch index as String.")
end

"""
    get_bus_vm_pu(net, bus_name)

Return solved voltage magnitude in p.u. for bus `bus_name`.
"""
function get_bus_vm_pu(net::Net, bus_name::String)
  bus = geNetBusIdx(net = net, busName = bus_name)
  return net.nodeVec[bus]._vm_pu
end

"""
    get_branch_p_from_to_mw(net, from_bus, to_bus)

Return active power flow in MW for the oriented branch direction
`from_bus -> to_bus`.
"""
function get_branch_p_from_to_mw(net::Net, from_bus::String, to_bus::String)
  br = getNetBranch(net = net, fromBus = from_bus, toBus = to_bus)
  from = geNetBusIdx(net = net, busName = from_bus)
  return br.fromBus == from ? br.fBranchFlow.pFlow : br.tBranchFlow.pFlow
end

"""
    get_branch_q_from_to_mvar(net, from_bus, to_bus)

Return reactive power flow in MVAr for the oriented branch direction
`from_bus -> to_bus`.
"""
function get_branch_q_from_to_mvar(net::Net, from_bus::String, to_bus::String)
  br = getNetBranch(net = net, fromBus = from_bus, toBus = to_bus)
  from = geNetBusIdx(net = net, busName = from_bus)
  return br.fromBus == from ? br.fBranchFlow.qFlow : br.tBranchFlow.qFlow
end

"""
    addTapController!(net; ...)

Add and validate a transformer tap controller.

Validation rules:
- One active controller per transformer.
- Required target fields must be set according to `mode`.
- `control_ratio` / `control_phase` must match the chosen mode.
"""
function addTapController!(net::Net; trafo::String, mode::Symbol, target_bus::Union{Nothing,String}=nothing, target_branch::Union{Nothing,Tuple{String,String}}=nothing,
  target_vm_pu::Union{Nothing,Float64}=nothing, p_target_mw::Union{Nothing,Float64}=nothing, q_target_mvar::Union{Nothing,Float64}=nothing,
  control_ratio::Bool=true, control_phase::Bool=false, is_discrete::Bool=true, deadband_vm_pu::Float64=1e-3, deadband_p_mw::Float64=0.5,
  max_outer_iters::Int=20, enabled::Bool=true)

  _ = _find_trafo_branch(net, trafo)
  any(c -> c.trafo == trafo && c.enabled, net.tapControllers) && error("TapController: only one active controller per transformer is allowed.")
  mode in (:voltage, :branch_active_power, :voltage_and_branch_active_power) || error("TapController: unsupported mode=$(mode)")

  if mode in (:voltage, :voltage_and_branch_active_power)
    isnothing(target_bus) && error("TapController: target_bus is required for mode=$(mode)")
    isnothing(target_vm_pu) && error("TapController: target_vm_pu is required for mode=$(mode)")
    !control_ratio && error("TapController: control_ratio must be true for voltage control")
  end
  if mode in (:branch_active_power, :voltage_and_branch_active_power)
    isnothing(target_branch) && error("TapController: target_branch is required for mode=$(mode)")
    isnothing(p_target_mw) && error("TapController: p_target_mw is required for mode=$(mode)")
    !control_phase && error("TapController: control_phase must be true for branch active power control")
  end

  push!(net.tapControllers, TapController(;
    trafo = trafo,
    mode = mode,
    target_bus = target_bus,
    target_branch = target_branch,
    target_vm_pu = target_vm_pu,
    p_target_mw = p_target_mw,
    q_target_mvar = q_target_mvar,
    control_ratio = control_ratio,
    control_phase = control_phase,
    is_discrete = is_discrete,
    deadband_vm_pu = deadband_vm_pu,
    deadband_p_mw = deadband_p_mw,
    max_outer_iters = max_outer_iters,
    enabled = enabled,
  ))
  return net
end

"""
    _phase_probe_direction(...)

Internal helper: determines the empirical sign of `ΔP_from_to` for a positive
phase increment (`+phase_step_deg`) on the controlled transformer and branch.
Returns `-1`, `0`, or `+1`.
"""
function _phase_probe_direction(net::Net, br::Branch, ctrl::TapController, max_ite::Int, tol::Float64, verbose::Int, opt_fd::Bool, opt_sparse::Bool, method::Symbol)
  oldphi = br.phase_shift_deg
  step = br.phase_step_deg
  _, erg = runpf!(net, max_ite, tol, verbose; opt_fd = opt_fd, opt_sparse = opt_sparse, method = method)
  erg != 0 && return -1.0
  p0 = get_branch_p_from_to_mw(net, ctrl.target_branch[1], ctrl.target_branch[2])
  br.phase_shift_deg = clamp(oldphi + step, br.phase_min_deg, br.phase_max_deg)
  br.angle = br.phase_shift_deg
  _, erg2 = runpf!(net, max_ite, tol, verbose; opt_fd = opt_fd, opt_sparse = opt_sparse, method = method)
  p1 = erg2 == 0 ? get_branch_p_from_to_mw(net, ctrl.target_branch[1], ctrl.target_branch[2]) : p0
  br.phase_shift_deg = oldphi
  br.angle = oldphi
  return sign(p1 - p0)
end

"""
    run_tap_controllers_outer!(net; ...)

Run outer-loop transformer control around `runpf!`.

Loop per iteration:
1. Solve PF.
2. Evaluate controller errors.
3. Update ratio and/or phase (discrete or continuous-like step).
4. Re-solve PF.
5. Stop when all controllers are converged, blocked at limits, or no movement is
   possible; otherwise stop at `max_outer_iters`.

Returns `(iterations, erg)` where `erg == 0` means successful PF termination.
"""
function run_tap_controllers_outer!(net::Net; max_ite::Int=30, tol::Float64=1e-6, verbose::Int=0, opt_fd::Bool=false, opt_sparse::Bool=false, method::Symbol=:rectangular)
  isempty(net.tapControllers) && return (0, 0)
  _, erg = runpf!(net, max_ite, tol, verbose; opt_fd = opt_fd, opt_sparse = opt_sparse, method = method)
  erg != 0 && return (0, erg)
  calcNetLosses!(net)
  calcLinkFlowsKCL!(net)

  max_outer = maximum(c.max_outer_iters for c in net.tapControllers if c.enabled)
  for it in 1:max_outer
    all_done = true
    any_moved = false

    for ctrl in net.tapControllers
      ctrl.enabled || continue
      br = _find_trafo_branch(net, ctrl.trafo)
      moved = false
      converged_v = true
      converged_p = true

      if ctrl.mode in (:voltage, :voltage_and_branch_active_power)
        vm = get_bus_vm_pu(net, ctrl.target_bus)
        ctrl.achieved_vm_pu = vm
        e_v = vm^2 - ctrl.target_vm_pu^2
        converged_v = abs(e_v) <= ctrl.deadband_vm_pu
        if !converged_v && ctrl.control_ratio
          Δ = ctrl.is_discrete ? br.tap_step : 0.25 * br.tap_step
          new_ratio = clamp(br.tap_ratio + (e_v < 0.0 ? Δ : -Δ), br.tap_min, br.tap_max)
          any_moved = any_moved || (new_ratio != br.tap_ratio)
          br.tap_ratio = new_ratio
          br.ratio = new_ratio
          ctrl.at_limit = ctrl.at_limit || isapprox(new_ratio, br.tap_min; atol = 1e-12) || isapprox(new_ratio, br.tap_max; atol = 1e-12)
        end
      end

      if ctrl.mode in (:branch_active_power, :voltage_and_branch_active_power)
        p = get_branch_p_from_to_mw(net, ctrl.target_branch[1], ctrl.target_branch[2])
        ctrl.achieved_p_mw = p
        e_p = p - ctrl.p_target_mw
        converged_p = abs(e_p) <= ctrl.deadband_p_mw
        if !converged_p && ctrl.control_phase
          direction = _phase_probe_direction(net, br, ctrl, max_ite, tol, 0, opt_fd, opt_sparse, method)
          direction == 0.0 && (direction = -1.0)
          Δ = ctrl.is_discrete ? br.phase_step_deg : 0.25 * br.phase_step_deg
          step = (e_p < 0.0) ? direction * Δ : -direction * Δ
          new_phase = clamp(br.phase_shift_deg + step, br.phase_min_deg, br.phase_max_deg)
          any_moved = any_moved || (new_phase != br.phase_shift_deg)
          br.phase_shift_deg = new_phase
          br.angle = new_phase
          ctrl.at_limit = ctrl.at_limit || isapprox(new_phase, br.phase_min_deg; atol = 1e-12) || isapprox(new_phase, br.phase_max_deg; atol = 1e-12)
        end
      end

      ctrl.converged = converged_v && converged_p
      ctrl.outer_iters = it
      ctrl.status = ctrl.converged ? :converged : (ctrl.at_limit ? :at_limit : :active)
      all_done = all_done && (ctrl.converged || ctrl.at_limit)
      all_done || (all_done = all_done && !moved)
    end

    if all_done || !any_moved
      return (it, 0)
    end

    _, erg2 = runpf!(net, max_ite, tol, verbose; opt_fd = opt_fd, opt_sparse = opt_sparse, method = method)
    erg2 != 0 && return (it, erg2)
    calcNetLosses!(net)
    calcLinkFlowsKCL!(net)
  end

  for ctrl in net.tapControllers
    ctrl.enabled || continue
    ctrl.status = :max_outer_iters
  end
  return (max_outer, 0)
end

"""
    printTapControllerSummary(io, net)

Print a compact controller summary block for all configured tap controllers.
"""
function printTapControllerSummary(io::IO, net::Net)
  isempty(net.tapControllers) && return
  println(io, "\nTransformer Control Summary")
  println(io, "---------------------------")
  for ctrl in net.tapControllers
    br = _find_trafo_branch(net, ctrl.trafo)
    println(io, "Transformer ", ctrl.trafo)
    println(io, "  mode               : ", ctrl.mode)
    println(io, "  target_bus         : ", isnothing(ctrl.target_bus) ? "-" : ctrl.target_bus)
    println(io, "  target_vm_pu       : ", isnothing(ctrl.target_vm_pu) ? "-" : @sprintf("%.4f", ctrl.target_vm_pu))
    println(io, "  achieved_vm_pu     : ", isnothing(ctrl.achieved_vm_pu) ? "-" : @sprintf("%.4f", ctrl.achieved_vm_pu))
    println(io, "  target_branch      : ", isnothing(ctrl.target_branch) ? "-" : string(ctrl.target_branch))
    println(io, "  p_target_mw        : ", isnothing(ctrl.p_target_mw) ? "-" : @sprintf("%.2f", ctrl.p_target_mw))
    println(io, "  achieved_p_mw      : ", isnothing(ctrl.achieved_p_mw) ? "-" : @sprintf("%.2f", ctrl.achieved_p_mw))
    println(io, "  tap_ratio          : ", @sprintf("%.5f", br.tap_ratio))
    println(io, "  phase_shift_deg    : ", @sprintf("%.5f", br.phase_shift_deg))
    println(io, "  discrete           : ", ctrl.is_discrete)
    println(io, "  converged          : ", ctrl.converged)
    println(io, "  at_limit           : ", ctrl.at_limit)
    println(io, "  status             : ", ctrl.status)
  end
end
