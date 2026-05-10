# Copyright 2023–2026 Udo Schmitz

@inline _voltage_within_deadband(vm::Float64, target_vm_pu::Float64, deadband_vm_pu::Float64)::Bool = abs(vm - target_vm_pu) <= deadband_vm_pu
@inline _voltage_control_error(vm::Float64, target_vm_pu::Float64, metric::Symbol)::Float64 = metric == :vm2 ? vm^2 - target_vm_pu^2 : vm - target_vm_pu

"""
    _tap_controllers(net)

Collect all tap controllers configured in `net` from transformer windings.
"""
function _tap_controllers(net::Net)::Vector{PowerTransformerControl}
  controllers = PowerTransformerControl[]
  seen = IdDict{PowerTransformerControl,Bool}()

  for trafo in net.trafos
    for winding in (trafo.side1, trafo.side2, trafo.side3)
      isnothing(winding) && continue
      for ctrl in winding.controls
        haskey(seen, ctrl) && continue
        seen[ctrl] = true
        push!(controllers, ctrl)
      end
    end
  end
  return controllers
end

"""
    clearTapControllers!(net)

Remove all tap controllers from transformer windings in `net`.
"""
function clearTapControllers!(net::Net)
  for trafo in net.trafos
    empty!(trafo.side1.controls)
    empty!(trafo.side2.controls)
    if !isnothing(trafo.side3)
      empty!(trafo.side3.controls)
    end
  end
  return net
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
  error("PowerTransformerControl: transformer $(name) not found. Use branch component name/id or branch index as String.")
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

@inline _controller_mode_label(mode::Symbol) = String(mode)

function _controller_type_label(ctrl::PowerTransformerControl, br::Branch)::String
  ratio = ctrl.control_ratio && br.has_ratio_tap
  phase = ctrl.control_phase && br.has_phase_tap
  if ratio && phase
    return "OLTC+PST"
  elseif ratio
    return "OLTC"
  elseif phase
    return "PST"
  end
  return "fixed transformer"
end

function _controller_status_label(ctrl::PowerTransformerControl)::String
  !ctrl.enabled && return "inactive"
  ctrl.converged && return "converged"
  ctrl.at_limit && return "at_limit"
  ctrl.status == :max_outer_iters && return "not_converged"
  return "active"
end

@inline function _tap_position(value::Float64, neutral::Float64, step::Float64)
  step <= 0.0 && return missing
  return round(Int, (value - neutral) / step)
end

"""
    buildTapControllerReportRows(net::Net)

Build typed, machine-readable tap-controller report rows with engineering fields
for textual and DataFrame-based reporting.
"""
function buildTapControllerReportRows(net::Net)::Vector{NamedTuple}
  rows = NamedTuple[]
  for ctrl in _tap_controllers(net)
    br = _find_trafo_branch(net, ctrl.trafo)
    target_bus = ctrl.target_bus
    target_branch = ctrl.target_branch
    achieved_vm = isnothing(target_bus) ? missing : get_bus_vm_pu(net, target_bus)
    achieved_p = isnothing(target_branch) ? missing : get_branch_p_from_to_mw(net, target_branch[1], target_branch[2])
    push!(rows, (
      controller_name = string(br.comp.cName, " ", _controller_type_label(ctrl, br)),
      transformer_id = br.comp.cName,
      transformer_branch_index = br.branchIdx,
      from_bus = br.fromBus,
      to_bus = br.toBus,
      control_type = _controller_type_label(ctrl, br),
      mode = _controller_mode_label(ctrl.mode),
      target_bus = isnothing(target_bus) ? missing : target_bus,
      target_vm_pu = isnothing(ctrl.target_vm_pu) ? missing : ctrl.target_vm_pu,
      achieved_vm_pu = achieved_vm,
      target_branch_from = isnothing(target_branch) ? missing : target_branch[1],
      target_branch_to = isnothing(target_branch) ? missing : target_branch[2],
      p_target_mw = isnothing(ctrl.p_target_mw) ? missing : ctrl.p_target_mw,
      achieved_p_mw = achieved_p,
      tap_ratio = br.tap_ratio,
      phase_shift_deg = br.phase_shift_deg,
      ratio_tap_position = ctrl.is_discrete && br.has_ratio_tap ? _tap_position(br.tap_ratio, 1.0, br.tap_step) : missing,
      phase_tap_position = ctrl.is_discrete && br.has_phase_tap ? _tap_position(br.phase_shift_deg, 0.0, br.phase_step_deg) : missing,
      ratio_tap_min = br.has_ratio_tap ? br.tap_min : missing,
      ratio_tap_max = br.has_ratio_tap ? br.tap_max : missing,
      ratio_tap_step = br.has_ratio_tap ? br.tap_step : missing,
      phase_tap_min_deg = br.has_phase_tap ? br.phase_min_deg : missing,
      phase_tap_max_deg = br.has_phase_tap ? br.phase_max_deg : missing,
      phase_tap_step_deg = br.has_phase_tap ? br.phase_step_deg : missing,
      discrete = ctrl.is_discrete,
      converged = ctrl.converged,
      at_limit = ctrl.at_limit,
      status = _controller_status_label(ctrl),
      power_direction = isnothing(target_branch) ? missing : string(target_branch[1], " -> ", target_branch[2]),
    ))
  end
  return rows
end

"""
    addPowerTransformerControl!(net; ...)

Add and validate a transformer tap controller.

Validation rules:
- One active controller per transformer.
- Required target fields must be set according to `mode`.
- `control_ratio` / `control_phase` must match the chosen mode.
"""
function addPowerTransformerControl!(net::Net; trafo::String, mode::Symbol, target_bus::Union{Nothing,String}=nothing, target_branch::Union{Nothing,Tuple{String,String}}=nothing,
  target_vm_pu::Union{Nothing,Float64}=nothing, p_target_mw::Union{Nothing,Float64}=nothing, q_target_mvar::Union{Nothing,Float64}=nothing,
  control_ratio::Bool=true, control_phase::Bool=false, is_discrete::Bool=true, deadband_vm_pu::Float64=1e-3, deadband_p_mw::Float64=0.5,
  voltage_error_metric::Symbol=:vm, max_outer_iters::Int=20, enabled::Bool=true)

  br = _find_trafo_branch(net, trafo)
  trafo_obj = findfirst(t -> t.comp.cFrom_bus == br.comp.cFrom_bus && t.comp.cTo_bus == br.comp.cTo_bus, net.trafos)
  max_controls = isnothing(trafo_obj) ? 1 : max(1, net.trafos[trafo_obj].nController)
  controllers = _tap_controllers(net)
  n_active = 0
  for c in controllers
    c.enabled || continue
    cbr = _find_trafo_branch(net, c.trafo)
    cbr.branchIdx == br.branchIdx && (n_active += 1)
  end
  n_active >= max_controls && error("PowerTransformerControl: transformer $(trafo) allows at most $(max_controls) active controller(s).")
  mode in (:voltage, :branch_active_power, :voltage_and_branch_active_power) || error("PowerTransformerControl: unsupported mode=$(mode)")

  if mode in (:voltage, :voltage_and_branch_active_power)
    isnothing(target_bus) && error("PowerTransformerControl: target_bus is required for mode=$(mode)")
    isnothing(target_vm_pu) && error("PowerTransformerControl: target_vm_pu is required for mode=$(mode)")
    !control_ratio && error("PowerTransformerControl: control_ratio must be true for voltage control")
  end
  if mode in (:branch_active_power, :voltage_and_branch_active_power)
    isnothing(target_branch) && error("PowerTransformerControl: target_branch is required for mode=$(mode)")
    isnothing(p_target_mw) && error("PowerTransformerControl: p_target_mw is required for mode=$(mode)")
    !control_phase && error("PowerTransformerControl: control_phase must be true for branch active power control")
  end

  ctrl = PowerTransformerControl(;
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
    voltage_error_metric = voltage_error_metric,
    max_outer_iters = max_outer_iters,
    enabled = enabled)
  if isnothing(trafo_obj)
    error("PowerTransformerControl: no transformer object found for branch $(br.branchIdx)")
  end
  push!(net.trafos[trafo_obj].side1.controls, ctrl)
  return net
end

addTapController!(net::Net; kwargs...) = addPowerTransformerControl!(net; kwargs...)

"""
    _phase_probe_direction(...)

Internal helper: determines the empirical sign of `ΔP_from_to` for a positive
phase increment (`+phase_step_deg`) on the controlled transformer and branch.
Returns `-1`, `0`, or `+1`.
"""
function _phase_probe_direction(net::Net, br::Branch, ctrl::PowerTransformerControl, max_ite::Int, tol::Float64, verbose::Int, opt_fd::Bool, opt_sparse::Bool, method::Symbol; opt_flatstart::Bool=true, pv_table_rows::Int=30, validate_limits_after_pf::Bool=false, q_limit_violation_headroom::Float64=0.0, lock_pv_to_pq_buses::AbstractVector{Int}=Int[])
  oldphi = br.phase_shift_deg
  step = br.phase_step_deg
  _, erg = runpf!(net, max_ite, tol, verbose; opt_fd = opt_fd, opt_sparse = opt_sparse, method = method, opt_flatstart = opt_flatstart, pv_table_rows = pv_table_rows, validate_limits_after_pf = validate_limits_after_pf, q_limit_violation_headroom = q_limit_violation_headroom, lock_pv_to_pq_buses = lock_pv_to_pq_buses)
  erg != 0 && return -1.0
  p0 = get_branch_p_from_to_mw(net, ctrl.target_branch[1], ctrl.target_branch[2])
  br.phase_shift_deg = clamp(oldphi + step, br.phase_min_deg, br.phase_max_deg)
  br.angle = br.phase_shift_deg
  _, erg2 = runpf!(net, max_ite, tol, verbose; opt_fd = opt_fd, opt_sparse = opt_sparse, method = method, opt_flatstart = opt_flatstart, pv_table_rows = pv_table_rows, validate_limits_after_pf = validate_limits_after_pf, q_limit_violation_headroom = q_limit_violation_headroom, lock_pv_to_pq_buses = lock_pv_to_pq_buses)
  p1 = erg2 == 0 ? get_branch_p_from_to_mw(net, ctrl.target_branch[1], ctrl.target_branch[2]) : p0
  br.phase_shift_deg = oldphi
  br.angle = oldphi
  return sign(p1 - p0)
end

"""
    _ratio_probe_direction(...)

Internal helper: determines the empirical sign of `ΔVm_target` for a positive
ratio increment (`+tap_step`) on the controlled transformer and target bus.
Returns `-1`, `0`, or `+1`.
"""
function _ratio_probe_direction(net::Net, br::Branch, ctrl::PowerTransformerControl, max_ite::Int, tol::Float64, verbose::Int, opt_fd::Bool, opt_sparse::Bool, method::Symbol; opt_flatstart::Bool=true, pv_table_rows::Int=30, validate_limits_after_pf::Bool=false, q_limit_violation_headroom::Float64=0.0, lock_pv_to_pq_buses::AbstractVector{Int}=Int[])
  oldratio = br.tap_ratio
  step = br.tap_step
  _, erg = runpf!(net, max_ite, tol, verbose; opt_fd = opt_fd, opt_sparse = opt_sparse, method = method, opt_flatstart = opt_flatstart, pv_table_rows = pv_table_rows, validate_limits_after_pf = validate_limits_after_pf, q_limit_violation_headroom = q_limit_violation_headroom, lock_pv_to_pq_buses = lock_pv_to_pq_buses)
  erg != 0 && return -1.0
  vm0 = get_bus_vm_pu(net, ctrl.target_bus)
  newratio = clamp(oldratio + step, br.tap_min, br.tap_max)
  if isapprox(newratio, oldratio; atol = 1e-12)
    return 0.0
  end
  br.tap_ratio = newratio
  br.ratio = newratio
  _, erg2 = runpf!(net, max_ite, tol, verbose; opt_fd = opt_fd, opt_sparse = opt_sparse, method = method, opt_flatstart = opt_flatstart, pv_table_rows = pv_table_rows, validate_limits_after_pf = validate_limits_after_pf, q_limit_violation_headroom = q_limit_violation_headroom, lock_pv_to_pq_buses = lock_pv_to_pq_buses)
  vm1 = erg2 == 0 ? get_bus_vm_pu(net, ctrl.target_bus) : vm0
  br.tap_ratio = oldratio
  br.ratio = oldratio
  return sign(vm1 - vm0)
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
function run_tap_controllers_outer!(net::Net; max_ite::Int=30, tol::Float64=1e-6, verbose::Int=0, opt_fd::Bool=false, opt_sparse::Bool=false, method::Symbol=:rectangular, autodamp::Bool=false, autodamp_min::Float64=1e-3, start_projection::Bool=false, start_projection_try_dc_start::Bool=true, start_projection_try_blend_scan::Bool=true, start_projection_blend_lambdas::AbstractVector{<:Real}=[0.25, 0.5, 0.75], start_projection_dc_angle_limit_deg::Float64=60.0, opt_flatstart::Bool=true, pv_table_rows::Int=30, validate_limits_after_pf::Bool=false, q_limit_violation_headroom::Float64=0.0, lock_pv_to_pq_buses::AbstractVector{Int}=Int[])
  controllers = _tap_controllers(net)
  isempty(controllers) && return (0, 0)
  if !any(c -> c.enabled, controllers)
    return runpf!(net, max_ite, tol, verbose; opt_fd = opt_fd, opt_sparse = opt_sparse, method = method, autodamp = autodamp, autodamp_min = autodamp_min, start_projection = start_projection, start_projection_try_dc_start = start_projection_try_dc_start, start_projection_try_blend_scan = start_projection_try_blend_scan, start_projection_blend_lambdas = start_projection_blend_lambdas, start_projection_dc_angle_limit_deg = start_projection_dc_angle_limit_deg, opt_flatstart = opt_flatstart, pv_table_rows = pv_table_rows, validate_limits_after_pf = validate_limits_after_pf, q_limit_violation_headroom = q_limit_violation_headroom, lock_pv_to_pq_buses = lock_pv_to_pq_buses)
  end
  _, erg = runpf!(net, max_ite, tol, verbose; opt_fd = opt_fd, opt_sparse = opt_sparse, method = method, autodamp = autodamp, autodamp_min = autodamp_min, start_projection = start_projection, start_projection_try_dc_start = start_projection_try_dc_start, start_projection_try_blend_scan = start_projection_try_blend_scan, start_projection_blend_lambdas = start_projection_blend_lambdas, start_projection_dc_angle_limit_deg = start_projection_dc_angle_limit_deg, opt_flatstart = opt_flatstart, pv_table_rows = pv_table_rows, validate_limits_after_pf = validate_limits_after_pf, q_limit_violation_headroom = q_limit_violation_headroom, lock_pv_to_pq_buses = lock_pv_to_pq_buses)
  erg != 0 && return (0, erg)
  calcNetLosses!(net)
  calcLinkFlowsKCL!(net)

  max_outer = maximum(c.max_outer_iters for c in controllers if c.enabled)
  for it in 1:max_outer
    all_done = true
    any_moved = false

    for ctrl in controllers
      ctrl.enabled || continue
      br = _find_trafo_branch(net, ctrl.trafo)
      moved = false
      converged_v = true
      converged_p = true

      if ctrl.mode in (:voltage, :voltage_and_branch_active_power)
        vm = get_bus_vm_pu(net, ctrl.target_bus)
        ctrl.achieved_vm_pu = vm
        converged_v = _voltage_within_deadband(vm, ctrl.target_vm_pu, ctrl.deadband_vm_pu)
        if !converged_v && ctrl.control_ratio
          e_v = _voltage_control_error(vm, ctrl.target_vm_pu, ctrl.voltage_error_metric)
          direction = _ratio_probe_direction(net, br, ctrl, max_ite, tol, 0, opt_fd, opt_sparse, method; opt_flatstart = opt_flatstart, pv_table_rows = pv_table_rows, validate_limits_after_pf = validate_limits_after_pf, q_limit_violation_headroom = q_limit_violation_headroom, lock_pv_to_pq_buses = lock_pv_to_pq_buses)
          direction == 0.0 && (direction = -1.0)
          Δ = ctrl.is_discrete ? br.tap_step : 0.25 * br.tap_step
          step = (e_v < 0.0) ? direction * Δ : -direction * Δ
          new_ratio = clamp(br.tap_ratio + step, br.tap_min, br.tap_max)
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
          direction = _phase_probe_direction(net, br, ctrl, max_ite, tol, 0, opt_fd, opt_sparse, method; opt_flatstart = opt_flatstart, pv_table_rows = pv_table_rows, validate_limits_after_pf = validate_limits_after_pf, q_limit_violation_headroom = q_limit_violation_headroom, lock_pv_to_pq_buses = lock_pv_to_pq_buses)
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

    _, erg2 = runpf!(net, max_ite, tol, verbose; opt_fd = opt_fd, opt_sparse = opt_sparse, method = method, autodamp = autodamp, autodamp_min = autodamp_min, start_projection = start_projection, start_projection_try_dc_start = start_projection_try_dc_start, start_projection_try_blend_scan = start_projection_try_blend_scan, start_projection_blend_lambdas = start_projection_blend_lambdas, start_projection_dc_angle_limit_deg = start_projection_dc_angle_limit_deg, opt_flatstart = opt_flatstart, pv_table_rows = pv_table_rows, validate_limits_after_pf = validate_limits_after_pf, q_limit_violation_headroom = q_limit_violation_headroom, lock_pv_to_pq_buses = lock_pv_to_pq_buses)
    erg2 != 0 && return (it, erg2)
    calcNetLosses!(net)
    calcLinkFlowsKCL!(net)
  end

  for ctrl in controllers
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
  if isempty(_tap_controllers(net))
    println(io, "\nControl")
    println(io, "-------")
    println(io, "Transformer controls: none")
    return
  end
  rows = buildTapControllerReportRows(net)
  println(io, "\nTransformer Control Summary")
  println(io, "---------------------------")
  println(io, "Power sign convention: achieved_p_mw is positive in the configured target branch direction (from -> to).")
  for row in rows
    println(io, row.controller_name, " (", row.transformer_id, ", ", row.from_bus, " -> ", row.to_bus, ")")
    println(io, "  controller type    : ", row.control_type)
    println(io, "  mode               : ", row.mode)
    println(io, "  target bus         : ", ismissing(row.target_bus) ? "-" : row.target_bus)
    println(io, "  target Vm          : ", ismissing(row.target_vm_pu) ? "-" : @sprintf("%.4f pu", row.target_vm_pu))
    println(io, "  achieved Vm        : ", ismissing(row.achieved_vm_pu) ? "-" : @sprintf("%.4f pu", row.achieved_vm_pu))
    println(io, "  target branch      : ", ismissing(row.target_branch_from) ? "-" : string(row.target_branch_from, " -> ", row.target_branch_to))
    println(io, "  target P           : ", ismissing(row.p_target_mw) ? "-" : @sprintf("%.3f MW", row.p_target_mw))
    println(io, "  achieved P         : ", ismissing(row.achieved_p_mw) ? "-" : @sprintf("%.3f MW", row.achieved_p_mw))
    println(io, "  tap ratio          : ", @sprintf("%.5f", row.tap_ratio))
    println(io, "  phase shift        : ", @sprintf("%.5f deg", row.phase_shift_deg))
    println(io, "  tap position       : ", ismissing(row.ratio_tap_position) ? "-" : @sprintf("%+d", row.ratio_tap_position))
    println(io, "  phase position     : ", ismissing(row.phase_tap_position) ? "-" : @sprintf("%+d", row.phase_tap_position))
    println(io, "  ratio range        : ", ismissing(row.ratio_tap_min) ? "-" : @sprintf("%.5f .. %.5f", row.ratio_tap_min, row.ratio_tap_max))
    println(io, "  ratio step         : ", ismissing(row.ratio_tap_step) ? "-" : @sprintf("%.5f", row.ratio_tap_step))
    println(io, "  phase range        : ", ismissing(row.phase_tap_min_deg) ? "-" : @sprintf("%.5f .. %.5f deg", row.phase_tap_min_deg, row.phase_tap_max_deg))
    println(io, "  phase step         : ", ismissing(row.phase_tap_step_deg) ? "-" : @sprintf("%.5f deg", row.phase_tap_step_deg))
    println(io, "  discrete           : ", row.discrete)
    println(io, "  converged          : ", row.converged)
    println(io, "  at_limit           : ", row.at_limit)
    println(io, "  status             : ", row.status)
    if !row.converged && row.at_limit
      println(io, "  status detail      : target not fully reached because tap/phase limit was hit")
    end
  end
end
