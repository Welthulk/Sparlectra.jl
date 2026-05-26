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

"""
Compatibility wrapper to perform AC power flow analysis on an existing Net.

Parameters:
- max_ite: Int, the maximum number of iterations for the power flow algorithm (default: 30).
- tol: Float64, tolerance for convergence criterion (default: 1e-6).
- casefile: String, the name of the case file to load.
- path: Union{Nothing,String}, the path to the case file (default: nothing).
- verbose: Int, verbosity level for output (default: 0).
- printResultToFile: Bool, flag to print results to a file (default: false).
- printResultAnyCase: Bool, flag to print results even if the power flow fails (default: false).
- autodamp: Bool, enable residual-based Newton step backtracking for `method = :rectangular`.
- autodamp_min: Float64, minimum trial step length for automatic damping.
- matpower_ratio: `:normal`/`"normal"` keeps MATPOWER branch `TAP` values as stored; `:reciprocal`/`"reciprocal"` imports reciprocal transformer ratios.
- reference_vm_pu/reference_va_deg: optional MATPOWER slack/reference voltage override used by the imported network and flat-start seed.

Returns:
- Net, the network object.
"""
function _apply_matpower_flatstart_modes!(net::Net, mpc; voltage_mode = :classic, angle_mode = :classic, profile_source = :matpower_reference, matpower_pv_voltage_source = :gen_vg, matpower_pv_voltage_mismatch_tol_pu::Float64 = 1e-4, performance_profile = nothing)
  vmode = Symbol(voltage_mode)
  amode = Symbol(angle_mode)
  psource = Symbol(profile_source)
  vmode in (:classic, :pv_gen_vg, :pv_bus_vm, :all_bus_vm, :profile_blend) || error("flatstart_voltage_mode must be classic, pv_gen_vg, pv_bus_vm, all_bus_vm, or profile_blend")
  amode in (:classic, :dc, :bus_va_blend, :matpower_va) || error("flatstart_angle_mode must be classic, dc, bus_va_blend, or matpower_va")
  busrow = _perf_profile_time!(performance_profile, :start_projection_matpower_bus_map) do
    MatpowerIO.bus_row_index(mpc)
  end
  _perf_profile_time!(performance_profile, :start_projection_matpower_branch_map) do
    # The current start projection does not need a persistent branch map, but
    # recording the online-branch scan keeps large-case profile output explicit.
    count(e -> size(mpc.branch, 2) < 11 || mpc.branch[e, 11] != 0.0, axes(mpc.branch, 1))
  end
  pv_rows = _perf_profile_time!(performance_profile, :start_projection_matpower_reference_lookup) do
    MatpowerIO.pv_voltage_reference_rows(mpc; matpower_pv_voltage_source = matpower_pv_voltage_source, tol = matpower_pv_voltage_mismatch_tol_pu, warn = false)
  end
  imported_vset = Dict(row.busI => row.imported_vset for row in pv_rows)
  gen_vset = Dict(row.busI => (isempty(row.gen_vgs) ? row.bus_vm : row.gen_vgs[1]) for row in pv_rows)

  for k in eachindex(net.nodeVec)
    node = net.nodeVec[k]
    busI = haskey(net.busOrigIdxDict, k) ? net.busOrigIdxDict[k] : k
    haskey(busrow, busI) || continue
    r = busrow[busI]
    btype = Int(mpc.bus[r, 2])
    bus_vm = Float64(mpc.bus[r, 8])
    bus_va = Float64(mpc.bus[r, 9])

    if vmode == :pv_gen_vg && (btype == 2 || btype == 3)
      setVmVa!(node = node, vm_pu = get(gen_vset, busI, bus_vm))
    elseif vmode == :pv_bus_vm && (btype == 2 || btype == 3)
      setVmVa!(node = node, vm_pu = bus_vm)
    elseif vmode == :all_bus_vm
      setVmVa!(node = node, vm_pu = bus_vm)
    elseif vmode == :profile_blend
      psource == :matpower_reference || error("profile_blend currently requires profile_source=:matpower_reference for MATPOWER imports.")
      current_vm = node._vm_pu === nothing ? 1.0 : Float64(node._vm_pu)
      setVmVa!(node = node, vm_pu = 0.5 * (current_vm + bus_vm))
    elseif vmode == :classic && (btype == 2 || btype == 3)
      setVmVa!(node = node, vm_pu = get(imported_vset, busI, bus_vm))
    end

    if amode == :matpower_va
      setVmVa!(node = node, vm_pu = node._vm_pu === nothing ? bus_vm : Float64(node._vm_pu), va_deg = bus_va)
    elseif amode == :bus_va_blend
      current_va = node._va_deg === nothing ? 0.0 : Float64(node._va_deg)
      setVmVa!(node = node, vm_pu = node._vm_pu === nothing ? bus_vm : Float64(node._vm_pu), va_deg = 0.5 * (current_va + bus_va))
    elseif amode == :dc
      # DC starts are handled by start_projection_try_dc_start in the rectangular solver.
      nothing
    end
  end
  return nothing
end

function _print_ac_pf_nonconvergence(method::Symbol, net::Net)
  rect_status = method === :rectangular ? rectangular_pf_status(net) : nothing
  if rect_status !== nothing && getproperty(rect_status, :numerical_converged)
    @printf("Power flow numerical solve converged, but final status is %s (reason=%s).\n", String(rect_status.status), rect_status.reason_text)
  else
    println("Newton-Raphson did not converge")
  end
  return nothing
end

function _rectangular_public_outcome(rect_status)::Symbol
  rect_status === nothing && return :solver_error
  status = hasproperty(rect_status, :status) ? getproperty(rect_status, :status) : :solver_error
  status === :converged && return :converged
  status === :converged_with_limit_warnings && return :converged_with_limit_warnings
  status === :converged_limits_failed && return :converged_limits_failed
  status === :singular_jacobian && return :singular_jacobian
  status === :not_converged && return :not_converged
  return :solver_error
end

function _has_numerical_voltage_solution(method::Symbol, erg::Int, rect_status)::Bool
  if method === :rectangular && rect_status !== nothing
    return Bool(getproperty(rect_status, :numerical_converged))
  end
  return erg == 0
end

function _solution_available_for_outcome(outcome::Symbol)::Bool
  return outcome in (:converged, :converged_with_limit_warnings, :converged_limits_failed)
end

function _rectangular_limit_validation_status(rect_status)::Symbol
  rect_status === nothing && return :skip
  Bool(getproperty(rect_status, :q_limit_active_set_ok)) && return :ok
  Bool(getproperty(rect_status, :numerical_converged)) ? :fail : :skip
end

function _build_pf_outcome_status(method::Symbol, erg::Int, ite::Int, etime::Float64, rect_status)
  if method === :rectangular && rect_status !== nothing
    outcome = _rectangular_public_outcome(rect_status)
    return (
      outcome = outcome,
      numerical_converged = Bool(getproperty(rect_status, :numerical_converged)),
      solution_available = _solution_available_for_outcome(outcome),
      limit_validation_status = _rectangular_limit_validation_status(rect_status),
      final_converged = Bool(getproperty(rect_status, :final_converged)),
      reason = Symbol(getproperty(rect_status, :reason)),
      reason_text = String(getproperty(rect_status, :reason_text)),
      final_mismatch = Float64(getproperty(rect_status, :final_mismatch)),
      iterations = ite,
      elapsed_s = etime,
    )
  end
  outcome = erg == 0 ? :converged : :not_converged
  return (outcome = outcome, numerical_converged = (erg == 0), solution_available = (erg == 0), limit_validation_status = :skip, final_converged = (erg == 0), reason = (erg == 0 ? :none : :nr_mismatch_not_converged), reason_text = (erg == 0 ? "none" : "NR mismatch did not converge"), final_mismatch = NaN, iterations = ite, elapsed_s = etime)
end

_legacy_erg_from_control_status(status::Symbol)::Int = status == :pf_failed ? 1 : 0

function _legacy_powerflow_config(; max_ite::Int, tol::Float64, method::Symbol, autodamp::Bool, autodamp_min::Float64, opt_flatstart::Bool, start_projection::Bool, start_projection_try_dc_start::Bool, start_projection_try_blend_scan::Bool, start_projection_branch_guard::Bool, start_projection_measure_candidates::Bool, start_projection_accept_unmeasured_dc_start::Bool, start_projection_blend_lambdas, start_projection_dc_angle_limit_deg::Float64, qlimit_start_iter::Int, qlimit_start_mode::Symbol, qlimit_auto_q_delta_pu::Float64, lock_pv_to_pq_buses, qlimit_trace_buses, qlimit_guard::Bool, qlimit_guard_min_q_range_pu::Float64, qlimit_guard_zero_range_mode::Symbol, qlimit_guard_narrow_range_mode::Symbol, qlimit_guard_log::Bool, qlimit_guard_max_switches::Int, qlimit_guard_accept_bounded_violations::Bool, qlimit_guard_max_remaining_violations::Int, qlimit_guard_freeze_after_repeated_switching::Bool, qlimit_guard_violation_mode::Symbol, qlimit_guard_violation_threshold_pu::Float64)
  return PowerFlowConfig(
    method = method,
    tol = tol,
    max_iter = max_ite,
    sparse = true,
    autodamp = autodamp,
    autodamp_min = autodamp_min,
    start_mode = StartModeConfig(
      flatstart = opt_flatstart,
      start_projection = start_projection,
      try_dc_start = start_projection_try_dc_start,
      try_blend_scan = start_projection_try_blend_scan,
      branch_guard = start_projection_branch_guard,
      measure_candidates = start_projection_measure_candidates,
      accept_unmeasured_dc_start = start_projection_accept_unmeasured_dc_start,
      blend_lambdas = Float64[x for x in start_projection_blend_lambdas],
      dc_angle_limit_deg = start_projection_dc_angle_limit_deg,
    ),
    qlimits = QLimitConfig(
      start_iter = qlimit_start_iter,
      start_mode = qlimit_start_mode,
      auto_q_delta_pu = qlimit_auto_q_delta_pu,
      guard = qlimit_guard,
      guard_min_q_range_pu = qlimit_guard_min_q_range_pu,
      guard_zero_range_mode = qlimit_guard_zero_range_mode,
      guard_narrow_range_mode = qlimit_guard_narrow_range_mode,
      guard_max_switches = qlimit_guard_max_switches,
      guard_freeze_after_repeated_switching = qlimit_guard_freeze_after_repeated_switching,
      guard_accept_bounded_violations = qlimit_guard_accept_bounded_violations,
      guard_max_remaining_violations = qlimit_guard_max_remaining_violations,
      guard_violation_mode = qlimit_guard_violation_mode,
      guard_violation_threshold_pu = qlimit_guard_violation_threshold_pu,
      guard_log = qlimit_guard_log,
      trace_buses = Int[x for x in qlimit_trace_buses],
      lock_pv_to_pq_buses = Int[x for x in lock_pv_to_pq_buses],
    ),
  )
end

_as_powerflow_config(config::PowerFlowConfig) = config
_as_powerflow_config(config::SparlectraConfig) = config.powerflow
_matpower_config_for_runner(::Union{Nothing,PowerFlowConfig}) = matpower_import_config()
_matpower_config_for_runner(config::SparlectraConfig) = config.matpower
_output_config_for_runner(::Union{Nothing,PowerFlowConfig}) = output_config()
_output_config_for_runner(config::SparlectraConfig) = config.output

"""
Run AC power flow using the public high-level ACP runner.

Use `run_acpflow(; net = net, ...)` for already constructed in-memory networks (preferred),
or `run_acpflow(; casefile = "case.m", path = "...", ...)` for file-based workflows.
Exactly one of `net` or `casefile` must be provided.
"""
function run_acpflow(;
  net::Union{Nothing,Net} = nothing,
  max_ite::Int = 30,
  tol::Float64 = 1e-6,
  casefile::Union{Nothing,String} = nothing,
  path::Union{Nothing,String} = nothing,
  verbose::Int = 0,  
  printResultToFile::Bool = false,
  printResultAnyCase::Bool = false,
  method::Symbol = :rectangular,
  autodamp::Bool = false,
  autodamp_min::Float64 = 1e-3, start_projection::Bool = false, start_projection_try_dc_start::Bool = true, start_projection_try_blend_scan::Bool = true, start_projection_branch_guard::Bool = true, start_projection_measure_candidates::Bool = true, start_projection_accept_unmeasured_dc_start::Bool = false, start_projection_blend_lambdas::AbstractVector{<:Real} = [0.25, 0.5, 0.75], start_projection_dc_angle_limit_deg::Float64 = 60.0, qlimit_start_iter::Int = 2, qlimit_start_mode::Symbol = :iteration, qlimit_auto_q_delta_pu::Float64 = 1e-4,
  qlimit_guard::Bool = false,
  qlimit_guard_min_q_range_pu::Float64 = 1e-4,
  qlimit_guard_zero_range_mode::Symbol = :lock_pq,
  qlimit_guard_narrow_range_mode::Symbol = :prefer_pq,
  qlimit_guard_max_switches::Int = 10,
  qlimit_guard_freeze_after_repeated_switching::Bool = true,
  qlimit_guard_accept_bounded_violations::Bool = false,
  qlimit_guard_max_remaining_violations::Int = 0,
  qlimit_guard_violation_mode::Symbol = :delayed_switch,
  qlimit_guard_violation_threshold_pu::Float64 = 1e-4,
  qlimit_guard_log::Bool = true,
  qlimit_trace_buses::AbstractVector{Int} = Int[],
  qlimit_lock_reason::Symbol = :manual,
  opt_flatstart::Union{Nothing,Bool} = nothing,
  show_results::Bool = true,
  show_compact_result::Bool = false,
  status_ref::Union{Nothing,Base.RefValue{Any}} = nothing,
  cooldown_iters::Int = 0,
  q_hyst_pu::Float64 = 0.0,
  pv_table_rows::Int = 30,
  check_q_limit_signs::Bool = false,
  autocorrect_q_limit_signs::Bool = false,
  validate_limits_after_pf::Bool = false,
  q_limit_violation_headroom::Float64 = 0.20,
  lock_pv_to_pq_buses::AbstractVector{Int} = Int[],
  enable_pq_gen_controllers::Union{Nothing,Bool} = nothing,
  bus_shunt_model = nothing,
  matpower_shift_sign::Union{Nothing,Real} = nothing,
  matpower_shift_unit = nothing,
  matpower_ratio = nothing,
  reference_vm_pu::Union{Nothing,Float64} = nothing,
  reference_va_deg::Union{Nothing,Float64} = nothing,
  matpower_pv_voltage_source = nothing,
  matpower_pv_voltage_mismatch_tol_pu::Union{Nothing,Float64} = nothing,
  flatstart_voltage_mode = :classic,
  flatstart_angle_mode = :classic,
  imported_matpower_case = nothing,
  config::Union{Nothing,PowerFlowConfig,SparlectraConfig} = nothing,
  performance_profile = nothing,
)
  if net !== nothing
    isnothing(casefile) || throw(ArgumentError("run_acpflow: pass either net or casefile, not both."))
    return _run_acpflow_net!(; net = net, max_ite = max_ite, tol = tol, verbose = verbose, printResultToFile = printResultToFile, printResultAnyCase = printResultAnyCase, method = method, autodamp = autodamp, autodamp_min = autodamp_min, start_projection = start_projection, start_projection_try_dc_start = start_projection_try_dc_start, start_projection_try_blend_scan = start_projection_try_blend_scan, start_projection_branch_guard = start_projection_branch_guard, start_projection_measure_candidates = start_projection_measure_candidates, start_projection_accept_unmeasured_dc_start = start_projection_accept_unmeasured_dc_start, start_projection_blend_lambdas = start_projection_blend_lambdas, start_projection_dc_angle_limit_deg = start_projection_dc_angle_limit_deg, qlimit_start_iter = qlimit_start_iter, qlimit_start_mode = qlimit_start_mode, qlimit_auto_q_delta_pu = qlimit_auto_q_delta_pu, show_results = show_results, lock_pv_to_pq_buses = lock_pv_to_pq_buses, opt_flatstart = isnothing(opt_flatstart) ? true : opt_flatstart, pv_table_rows = pv_table_rows, validate_limits_after_pf = validate_limits_after_pf, q_limit_violation_headroom = q_limit_violation_headroom, qlimit_trace_buses = qlimit_trace_buses, qlimit_lock_reason = qlimit_lock_reason, qlimit_guard = qlimit_guard, qlimit_guard_min_q_range_pu = qlimit_guard_min_q_range_pu, qlimit_guard_zero_range_mode = qlimit_guard_zero_range_mode, qlimit_guard_narrow_range_mode = qlimit_guard_narrow_range_mode, qlimit_guard_log = qlimit_guard_log, qlimit_guard_max_switches = qlimit_guard_max_switches, qlimit_guard_accept_bounded_violations = qlimit_guard_accept_bounded_violations, qlimit_guard_max_remaining_violations = qlimit_guard_max_remaining_violations, qlimit_guard_freeze_after_repeated_switching = qlimit_guard_freeze_after_repeated_switching, qlimit_guard_violation_mode = qlimit_guard_violation_mode, qlimit_guard_violation_threshold_pu = qlimit_guard_violation_threshold_pu, bus_shunt_model = isnothing(bus_shunt_model) ? net.bus_shunt_model : bus_shunt_model, config = config, performance_profile = performance_profile)
  end
  isnothing(casefile) && throw(ArgumentError("run_acpflow: casefile is required when net is not provided."))
  pf_config = isnothing(config) ? powerflow_config() : _as_powerflow_config(config)
  mat_cfg = _matpower_config_for_runner(config)
  out_cfg = _output_config_for_runner(config)
  max_ite = pf_config.max_iter
  tol = pf_config.tol
  method = pf_config.method
  autodamp = pf_config.autodamp
  autodamp_min = pf_config.autodamp_min
  opt_flatstart = isnothing(opt_flatstart) ? pf_config.start_mode.flatstart : opt_flatstart
  flatstart_angle_mode = pf_config.start_mode.angle_mode
  flatstart_voltage_mode = pf_config.start_mode.voltage_mode
  flatstart_profile_source = pf_config.start_mode.profile_source
  start_projection = pf_config.start_mode.start_projection
  start_projection_try_dc_start = pf_config.start_mode.try_dc_start
  start_projection_try_blend_scan = pf_config.start_mode.try_blend_scan
  start_projection_branch_guard = pf_config.start_mode.branch_guard
  start_projection_measure_candidates = pf_config.start_mode.measure_candidates
  start_projection_accept_unmeasured_dc_start = pf_config.start_mode.accept_unmeasured_dc_start
  start_projection_blend_lambdas = pf_config.start_mode.blend_lambdas
  start_projection_dc_angle_limit_deg = pf_config.start_mode.dc_angle_limit_deg
  qlimit_start_iter = pf_config.qlimits.start_iter
  qlimit_start_mode = pf_config.qlimits.start_mode
  qlimit_auto_q_delta_pu = pf_config.qlimits.auto_q_delta_pu
  q_hyst_pu = pf_config.qlimits.hysteresis_pu
  cooldown_iters = pf_config.qlimits.cooldown_iters
  qlimit_guard = pf_config.qlimits.guard
  qlimit_guard_min_q_range_pu = pf_config.qlimits.guard_min_q_range_pu
  qlimit_guard_zero_range_mode = pf_config.qlimits.guard_zero_range_mode
  qlimit_guard_narrow_range_mode = pf_config.qlimits.guard_narrow_range_mode
  qlimit_guard_max_switches = pf_config.qlimits.guard_max_switches
  qlimit_guard_freeze_after_repeated_switching = pf_config.qlimits.guard_freeze_after_repeated_switching
  qlimit_guard_accept_bounded_violations = pf_config.qlimits.guard_accept_bounded_violations
  qlimit_guard_max_remaining_violations = pf_config.qlimits.guard_max_remaining_violations
  qlimit_guard_violation_mode = pf_config.qlimits.guard_violation_mode
  qlimit_guard_violation_threshold_pu = pf_config.qlimits.guard_violation_threshold_pu
  qlimit_guard_log = pf_config.qlimits.guard_log
  qlimit_trace_buses = pf_config.qlimits.trace_buses
  lock_pv_to_pq_buses = pf_config.qlimits.lock_pv_to_pq_buses
  bus_shunt_model = isnothing(bus_shunt_model) ? mat_cfg.bus_shunt_model : bus_shunt_model
  matpower_shift_sign = isnothing(matpower_shift_sign) ? mat_cfg.shift_sign : matpower_shift_sign
  matpower_shift_unit = isnothing(matpower_shift_unit) ? mat_cfg.shift_unit : matpower_shift_unit
  matpower_ratio = isnothing(matpower_ratio) ? mat_cfg.ratio : matpower_ratio
  matpower_pv_voltage_source = isnothing(matpower_pv_voltage_source) ? mat_cfg.pv_voltage_source : matpower_pv_voltage_source
  matpower_pv_voltage_mismatch_tol_pu = isnothing(matpower_pv_voltage_mismatch_tol_pu) ? mat_cfg.pv_voltage_mismatch_tol_pu : matpower_pv_voltage_mismatch_tol_pu
  enable_pq_gen_controllers = isnothing(enable_pq_gen_controllers) ? mat_cfg.enable_pq_gen_controllers : enable_pq_gen_controllers
  show_results = show_results && out_cfg.logfile_results !== :off
  row_limit = out_cfg.result_table_max_rows > 0 ? out_cfg.result_table_max_rows : nothing
  result_mode = out_cfg.logfile_results === :compact ? :summary : out_cfg.logfile_results
  ext = lowercase(splitext(casefile)[2])
  myNet = nothing              # Initialize myNet variable
  in_path = nothing
  out_path = nothing
  if ext in (".m", ".jl")
    if path === nothing
      in_path  = joinpath(pwd(), "data", "mpower", strip(casefile))
      out_path = joinpath(pwd(), "data", "mpower")
      else
      in_path  = joinpath(path, strip(casefile))
      out_path = joinpath(path)
    end

    if !isfile(in_path)
      error("File $(in_path) not found")
    end

    myNet = _perf_profile_time!(performance_profile, :network_construction) do
      if imported_matpower_case === nothing
        createNetFromMatPowerFile(filename = in_path, log = (verbose > 0), flatstart = opt_flatstart, cooldown = cooldown_iters, q_hyst_pu = q_hyst_pu, enable_pq_gen_controllers = enable_pq_gen_controllers, bus_shunt_model = bus_shunt_model, matpower_shift_sign = matpower_shift_sign, matpower_shift_unit = matpower_shift_unit, matpower_ratio = matpower_ratio, reference_vm_pu = reference_vm_pu, reference_va_deg = reference_va_deg, matpower_pv_voltage_source = matpower_pv_voltage_source, matpower_pv_voltage_mismatch_tol_pu = matpower_pv_voltage_mismatch_tol_pu, profile = performance_profile)
      else
        createNetFromMatPowerCase(mpc = imported_matpower_case, log = (verbose > 0), flatstart = opt_flatstart, cooldown = cooldown_iters, q_hyst_pu = q_hyst_pu, enable_pq_gen_controllers = enable_pq_gen_controllers, bus_shunt_model = bus_shunt_model, matpower_shift_sign = matpower_shift_sign, matpower_shift_unit = matpower_shift_unit, matpower_ratio = matpower_ratio, reference_vm_pu = reference_vm_pu, reference_va_deg = reference_va_deg, matpower_pv_voltage_source = matpower_pv_voltage_source, matpower_pv_voltage_mismatch_tol_pu = matpower_pv_voltage_mismatch_tol_pu, preallocate_network = mat_cfg.preallocate_network, preallocate_min_buses = mat_cfg.preallocate_min_buses, profile = performance_profile)
      end
    end
    phase = performance_profile === nothing ? nothing : get(performance_profile, :network_construction, nothing)
    if performance_profile !== nothing && phase !== nothing && phase isa NamedTuple && hasproperty(phase, :bytes)
      performance_profile[:network_construction_allocated_bytes] = getproperty(phase, :bytes)
    end
    if opt_flatstart && (Symbol(flatstart_voltage_mode) != :classic || Symbol(flatstart_angle_mode) != :classic)
      mpc_init = _perf_profile_time!(performance_profile, :start_projection_matpower_reference_parse_lookup) do
        imported_matpower_case === nothing ? MatpowerIO.read_case(in_path; legacy_compat = true) : imported_matpower_case
      end
      _apply_matpower_flatstart_modes!(myNet, mpc_init; voltage_mode = flatstart_voltage_mode, angle_mode = flatstart_angle_mode, profile_source = flatstart_profile_source, matpower_pv_voltage_source = matpower_pv_voltage_source, matpower_pv_voltage_mismatch_tol_pu = matpower_pv_voltage_mismatch_tol_pu, performance_profile = performance_profile)
      if Symbol(flatstart_voltage_mode) in (:all_bus_vm, :profile_blend) || Symbol(flatstart_angle_mode) in (:bus_va_blend, :matpower_va)
        opt_flatstart = false
      end
    end
    if length(myNet.nodeVec) >= out_cfg.result_table_large_case_threshold_buses && out_cfg.logfile_results != :full
      result_mode = out_cfg.result_table_large_case_mode
    end
    
    if verbose > 1
      # --- DEBUG START ---
      @info "DEBUG Full net after import:"
      @printf "DEBUG Full net after import:\n"
      showNet(myNet, verbose = true)

      Y = createYBUS(net=myNet, sparse=false, printYBUS=false)  # dense for inspection      
      V0, slack = initialVrect(myNet; flatstart=myNet.flatstart)
      S = buildComplexSVec(myNet)

      @printf "DEBUG Yabs_max: %.6e\n" maximum(abs.(Y))
      @printf "DEBUG Ydiag_max: %.6e\n" maximum(abs.(diag(Y)))
      @printf "DEBUG Ydiag_imag_max: %.6e\n" maximum(abs.(imag.(diag(Y))))

      @printf "DEBUG slack bus index: %d\n" slack
      @printf "DEBUG initial V0: %s\n" string(V0)
      @printf "DEBUG case: name=%s baseMVA=%.1f flatstart=%s slack=%d\n" myNet.name myNet.baseMVA myNet.flatstart slack
      @printf "DEBUG sums (pu): sumP=%.6f sumQ=%.6f\n" sum(real.(S)) sum(imag.(S))
      @printf "DEBUG V0: Vslack=%s Vmin=%.6f Vmax=%.6f\n" string(V0[slack]) minimum(abs.(V0)) maximum(abs.(V0))
      @printf "DEBUG Y: Ydiag_min=%.6e Ydiag_max=%.6e Yabs_max=%.6e\n" minimum(abs.(diag(Y))) maximum(abs.(diag(Y))) maximum(abs.(Y))
      # --- DEBUG END ---
    end

  else
    error("File extension $(ext) not supported! Use .m or .jl")
  end

  # Resolve PV->PQ lock list (prefer original bus IDs; fallback to internal indices).
  lock_pv_to_pq_buses_resolved = Int[]
  if !isempty(lock_pv_to_pq_buses)
    orig_to_net = Dict{Int,Int}()
    for (net_idx, orig_idx) in myNet.busOrigIdxDict
      orig_to_net[orig_idx] = net_idx
    end
    for bus in lock_pv_to_pq_buses
      if haskey(orig_to_net, bus)
        push!(lock_pv_to_pq_buses_resolved, orig_to_net[bus])
      elseif 1 <= bus <= length(myNet.nodeVec)
        push!(lock_pv_to_pq_buses_resolved, bus)
      elseif verbose > 0
        @warn "lock_pv_to_pq_buses entry not found (ignored): $bus"
      end
    end
    unique!(lock_pv_to_pq_buses_resolved)
    sort!(lock_pv_to_pq_buses_resolved)
  end

  solver_config = _legacy_powerflow_config(
    max_ite = max_ite,
    tol = tol,
    method = method,
    autodamp = autodamp,
    autodamp_min = autodamp_min,
    opt_flatstart = opt_flatstart,
    start_projection = start_projection,
    start_projection_try_dc_start = start_projection_try_dc_start,
    start_projection_try_blend_scan = start_projection_try_blend_scan,
    start_projection_branch_guard = start_projection_branch_guard,
    start_projection_measure_candidates = start_projection_measure_candidates,
    start_projection_accept_unmeasured_dc_start = start_projection_accept_unmeasured_dc_start,
    start_projection_blend_lambdas = start_projection_blend_lambdas,
    start_projection_dc_angle_limit_deg = start_projection_dc_angle_limit_deg,
    qlimit_start_iter = qlimit_start_iter,
    qlimit_start_mode = qlimit_start_mode,
    qlimit_auto_q_delta_pu = qlimit_auto_q_delta_pu,
    lock_pv_to_pq_buses = lock_pv_to_pq_buses_resolved,
    qlimit_trace_buses = qlimit_trace_buses,
    qlimit_guard = qlimit_guard,
    qlimit_guard_min_q_range_pu = qlimit_guard_min_q_range_pu,
    qlimit_guard_zero_range_mode = qlimit_guard_zero_range_mode,
    qlimit_guard_narrow_range_mode = qlimit_guard_narrow_range_mode,
    qlimit_guard_log = qlimit_guard_log,
    qlimit_guard_max_switches = qlimit_guard_max_switches,
    qlimit_guard_accept_bounded_violations = qlimit_guard_accept_bounded_violations,
    qlimit_guard_max_remaining_violations = qlimit_guard_max_remaining_violations,
    qlimit_guard_freeze_after_repeated_switching = qlimit_guard_freeze_after_repeated_switching,
    qlimit_guard_violation_mode = qlimit_guard_violation_mode,
    qlimit_guard_violation_threshold_pu = qlimit_guard_violation_threshold_pu,
  )

  # Run power flow / optional tap-controller outer loop
  ite = 0

  erg = 2
  etime = @elapsed begin
    ite, erg = _perf_profile_time!(performance_profile, :solver_total) do
    controllers = collect_outer_controllers(myNet)
    if isempty(controllers)
      runpf!(myNet, solver_config; verbose = verbose, pv_table_rows = pv_table_rows, validate_limits_after_pf = validate_limits_after_pf, q_limit_violation_headroom = q_limit_violation_headroom, qlimit_lock_reason = qlimit_lock_reason, performance_profile = performance_profile)
    else
      control_result = run_control!(myNet; controllers = controllers, pf_config = solver_config, control_config = control_config(), verbose = verbose, performance_profile = performance_profile)
      # The legacy `erg` flag describes numerical PF success/failure.
      # Control-loop terminal states such as :blocked or :max_outer_iterations are
      # reported through `control_result` and `net.control_result`, not as PF failure.
      (control_result.last_pf_iterations, _legacy_erg_from_control_status(control_result.status))
      end
    end
  end
  solver_elapsed_s = begin
    timings = performance_profile isa AbstractDict ? get(performance_profile, :timings, nothing) : nothing
    row = timings isa AbstractDict ? get(timings, :solver_total, nothing) : nothing
    row isa NamedTuple && hasproperty(row, :elapsed_s) ? Float64(getproperty(row, :elapsed_s)) : nothing
  end

  if show_compact_result
    pv_to_pq_events = length(myNet.qLimitLog)
    pv_to_pq_buses = length(myNet.qLimitEvents)
    rect_status = method === :rectangular ? rectangular_pf_status(myNet) : nothing
    if rect_status !== nothing
      pf_status = _build_pf_outcome_status(method, erg, ite, etime, rect_status)
      limit_text = pf_status.limit_validation_status === :ok ? "OK" : pf_status.limit_validation_status === :fail ? "FAIL" : "SKIP"
      @printf("method=%-12s  outcome=%s  numerical_solution=%s  solution_available=%s  limit_validation=%s  final_converged=%s  iterations=%d  pv2pq_events=%d  pv2pq_buses=%d  final_mismatch=%.9g  time=%8.6f s  reason=%s", String(method), String(pf_status.outcome), pf_status.numerical_converged ? "OK" : "FAIL", string(pf_status.solution_available), limit_text, string(pf_status.final_converged), ite, pv_to_pq_events, pv_to_pq_buses, pf_status.final_mismatch, etime, pf_status.reason_text)
    else
      converged_text = erg == 0 ? "yes" : erg == 1 ? "no" : "error"
      @printf("method=%-12s  converged=%s  iterations=%d  pv2pq_events=%d  pv2pq_buses=%d  time=%8.6f s", String(method), converged_text, ite, pv_to_pq_events, pv_to_pq_buses, etime)
    end
    println()
  end

  if status_ref !== nothing
    rect_status = method === :rectangular ? rectangular_pf_status(myNet) : nothing
    outcome = _rectangular_public_outcome(rect_status)
    if rect_status === nothing
      ctrl_result = latest_control_result(myNet)
      status_ref[] = (converged = (erg == 0), erg = erg, iterations = ite, elapsed_s = etime, solver_elapsed_s = solver_elapsed_s, method = method, outcome = (erg == 0 ? :converged : :not_converged), control_status = isnothing(ctrl_result) ? :none : ctrl_result.status)
    else
      pf_status = _build_pf_outcome_status(method, erg, ite, etime, rect_status)
      ctrl_result = latest_control_result(myNet)
      status_ref[] = (converged = (erg == 0), erg = erg, iterations = ite, elapsed_s = etime, solver_elapsed_s = solver_elapsed_s, method = method, outcome = outcome, control_status = isnothing(ctrl_result) ? :none : ctrl_result.status,
                      numerical_solution = (pf_status.numerical_converged ? "OK" : "FAIL"),
                      solution_available = pf_status.solution_available,
                      limit_validation_status = pf_status.limit_validation_status,
                      numerical_converged = rect_status.numerical_converged,
                      q_limit_active_set_ok = rect_status.q_limit_active_set_ok,
                      final_converged = rect_status.final_converged,
                      status = rect_status.status,
                      nr_converged = rect_status.nr_converged,
                      active_set_converged = rect_status.active_set_converged,
                      reason = rect_status.reason,
                      reason_text = rect_status.reason_text,
                      final_mismatch = pf_status.final_mismatch,
                      pv_q_limit_violations = rect_status.pv_q_limit_violations,
                      ref_q_limit_violations = rect_status.ref_q_limit_violations,
                      final_pv_voltage_residual = rect_status.final_pv_voltage_residual)
    end
  end
  rect_status = method === :rectangular ? rectangular_pf_status(myNet) : nothing
  pf_status = _build_pf_outcome_status(method, erg, ite, etime, rect_status)
  outcome = pf_status.outcome
  if pf_status.solution_available || printResultAnyCase
    # Calculate network losses and print results
    _perf_profile_time!(performance_profile, :postprocess_losses_and_flows) do
      calcNetLosses!(myNet)
      calcLinkFlowsKCL!(myNet)
    end
    jpath = printResultToFile ? out_path : ""
    if show_results || printResultAnyCase
      _perf_profile_time!(performance_profile, :result_output) do
        printACPFlowResults(myNet, etime, ite, tol, printResultToFile, jpath; converged = (outcome == :converged), solver = method, solver_time_s = solver_elapsed_s, result_mode = result_mode, max_rows = row_limit)
        if outcome == :not_converged
          println("Diagnostic last-iterate voltages; Newton-Raphson did not converge.")
        elseif outcome == :converged_limits_failed
          println("Numerical voltage solution available, but Q-limit / active-set validation failed.")
          println("  final status: ", outcome)
          println("  reason: ", rect_status.reason_text)
          println("  remaining PV Q-limit violations: ", rect_status.pv_q_limit_violations)
          println("  REF Q-limit diagnostic violations: ", rect_status.ref_q_limit_violations)
          println("  final PV voltage residual: ", rect_status.final_pv_voltage_residual)
        end
      end
    end
  elseif erg == 1
    _print_ac_pf_nonconvergence(method, myNet)
  else
    @error "Errors during calculation of Newton-Raphson"
  end

  return myNet
end


function _run_acpflow_net!(; net::Net, max_ite::Int = 30, tol::Float64 = 1e-6, verbose::Int = 0, printResultToFile::Bool = false, printResultAnyCase::Bool = false, method::Symbol = :rectangular, autodamp::Bool = false, autodamp_min::Float64 = 1e-3, start_projection::Bool = false, start_projection_try_dc_start::Bool = true, start_projection_try_blend_scan::Bool = true, start_projection_branch_guard::Bool = true, start_projection_measure_candidates::Bool = true, start_projection_accept_unmeasured_dc_start::Bool = false, start_projection_blend_lambdas::AbstractVector{<:Real} = [0.25, 0.5, 0.75], start_projection_dc_angle_limit_deg::Float64 = 60.0, qlimit_start_iter::Int = 2, qlimit_start_mode::Symbol = :iteration, qlimit_auto_q_delta_pu::Float64 = 1e-4, show_results::Bool = true, lock_pv_to_pq_buses::AbstractVector{Int} = Int[], opt_flatstart::Bool = true, pv_table_rows::Int = 30, validate_limits_after_pf::Bool = false, q_limit_violation_headroom::Float64 = 0.0, qlimit_trace_buses::AbstractVector{Int} = Int[], qlimit_lock_reason::Symbol = :manual, qlimit_guard::Bool = false, qlimit_guard_min_q_range_pu::Float64 = 1e-4, qlimit_guard_zero_range_mode::Symbol = :lock_pq, qlimit_guard_narrow_range_mode::Symbol = :prefer_pq, qlimit_guard_log::Bool = true, qlimit_guard_max_switches::Int = 10, qlimit_guard_accept_bounded_violations::Bool = false, qlimit_guard_max_remaining_violations::Int = 0, qlimit_guard_freeze_after_repeated_switching::Bool = true, qlimit_guard_violation_mode::Symbol = :delayed_switch, qlimit_guard_violation_threshold_pu::Float64 = 1e-4, bus_shunt_model = net.bus_shunt_model, config::Union{Nothing,PowerFlowConfig,SparlectraConfig} = nothing, performance_profile = nothing)

  use_active_config = config !== nothing || (max_ite == 30 && tol == 1e-6 && method === :rectangular && !autodamp && autodamp_min == 1e-3 && !start_projection && start_projection_try_dc_start && start_projection_try_blend_scan && start_projection_branch_guard && start_projection_measure_candidates && !start_projection_accept_unmeasured_dc_start && collect(start_projection_blend_lambdas) == [0.25, 0.5, 0.75] && start_projection_dc_angle_limit_deg == 60.0 && qlimit_start_iter == 2 && qlimit_start_mode === :iteration && qlimit_auto_q_delta_pu == 1e-4 && isempty(lock_pv_to_pq_buses) && opt_flatstart && isempty(qlimit_trace_buses) && !qlimit_guard)
  pf_config = config === nothing ? powerflow_config() : _as_powerflow_config(config)
  out_cfg = _output_config_for_runner(config)
  if use_active_config
    max_ite = pf_config.max_iter
    tol = pf_config.tol
    method = pf_config.method
    autodamp = pf_config.autodamp
    autodamp_min = pf_config.autodamp_min
    opt_flatstart = pf_config.start_mode.flatstart
    flatstart_angle_mode = pf_config.start_mode.angle_mode
    flatstart_voltage_mode = pf_config.start_mode.voltage_mode
    start_projection = pf_config.start_mode.start_projection
    start_projection_try_dc_start = pf_config.start_mode.try_dc_start
    start_projection_try_blend_scan = pf_config.start_mode.try_blend_scan
    start_projection_branch_guard = pf_config.start_mode.branch_guard
    start_projection_measure_candidates = pf_config.start_mode.measure_candidates
    start_projection_accept_unmeasured_dc_start = pf_config.start_mode.accept_unmeasured_dc_start
    start_projection_blend_lambdas = pf_config.start_mode.blend_lambdas
    start_projection_dc_angle_limit_deg = pf_config.start_mode.dc_angle_limit_deg
    qlimit_start_iter = pf_config.qlimits.start_iter
    qlimit_start_mode = pf_config.qlimits.start_mode
    qlimit_auto_q_delta_pu = pf_config.qlimits.auto_q_delta_pu
    qlimit_guard = pf_config.qlimits.guard
    qlimit_guard_min_q_range_pu = pf_config.qlimits.guard_min_q_range_pu
    qlimit_guard_zero_range_mode = pf_config.qlimits.guard_zero_range_mode
    qlimit_guard_narrow_range_mode = pf_config.qlimits.guard_narrow_range_mode
    qlimit_guard_max_switches = pf_config.qlimits.guard_max_switches
    qlimit_guard_freeze_after_repeated_switching = pf_config.qlimits.guard_freeze_after_repeated_switching
    qlimit_guard_accept_bounded_violations = pf_config.qlimits.guard_accept_bounded_violations
    qlimit_guard_max_remaining_violations = pf_config.qlimits.guard_max_remaining_violations
    qlimit_guard_violation_mode = pf_config.qlimits.guard_violation_mode
    qlimit_guard_violation_threshold_pu = pf_config.qlimits.guard_violation_threshold_pu
    qlimit_guard_log = pf_config.qlimits.guard_log
    qlimit_trace_buses = pf_config.qlimits.trace_buses
    lock_pv_to_pq_buses = pf_config.qlimits.lock_pv_to_pq_buses
  end
  show_results = show_results && out_cfg.logfile_results !== :off
  row_limit = out_cfg.result_table_max_rows > 0 ? out_cfg.result_table_max_rows : nothing
  result_mode = out_cfg.logfile_results === :compact ? :summary : out_cfg.logfile_results
  if length(net.nodeVec) >= out_cfg.result_table_large_case_threshold_buses && out_cfg.logfile_results != :full
    result_mode = out_cfg.result_table_large_case_mode
  end

  requested_shunt_model = normalize_bus_shunt_model(bus_shunt_model)
  requested_shunt_model == net.bus_shunt_model || error("run_acpflow(net=...): bus_shunt_model must be set when constructing/importing the Net; got $(requested_shunt_model) for a net configured as $(net.bus_shunt_model).")

  solver_config = _legacy_powerflow_config(
    max_ite = max_ite,
    tol = tol,
    method = method,
    autodamp = autodamp,
    autodamp_min = autodamp_min,
    opt_flatstart = opt_flatstart,
    start_projection = start_projection,
    start_projection_try_dc_start = start_projection_try_dc_start,
    start_projection_try_blend_scan = start_projection_try_blend_scan,
    start_projection_branch_guard = start_projection_branch_guard,
    start_projection_measure_candidates = start_projection_measure_candidates,
    start_projection_accept_unmeasured_dc_start = start_projection_accept_unmeasured_dc_start,
    start_projection_blend_lambdas = start_projection_blend_lambdas,
    start_projection_dc_angle_limit_deg = start_projection_dc_angle_limit_deg,
    qlimit_start_iter = qlimit_start_iter,
    qlimit_start_mode = qlimit_start_mode,
    qlimit_auto_q_delta_pu = qlimit_auto_q_delta_pu,
    lock_pv_to_pq_buses = lock_pv_to_pq_buses,
    qlimit_trace_buses = qlimit_trace_buses,
    qlimit_guard = qlimit_guard,
    qlimit_guard_min_q_range_pu = qlimit_guard_min_q_range_pu,
    qlimit_guard_zero_range_mode = qlimit_guard_zero_range_mode,
    qlimit_guard_narrow_range_mode = qlimit_guard_narrow_range_mode,
    qlimit_guard_log = qlimit_guard_log,
    qlimit_guard_max_switches = qlimit_guard_max_switches,
    qlimit_guard_accept_bounded_violations = qlimit_guard_accept_bounded_violations,
    qlimit_guard_max_remaining_violations = qlimit_guard_max_remaining_violations,
    qlimit_guard_freeze_after_repeated_switching = qlimit_guard_freeze_after_repeated_switching,
    qlimit_guard_violation_mode = qlimit_guard_violation_mode,
    qlimit_guard_violation_threshold_pu = qlimit_guard_violation_threshold_pu,
  )

  # Run power flow / optional tap-controller outer loop
  ite = 0

  erg = 2
  etime = @elapsed begin
    ite, erg = _perf_profile_time!(performance_profile, :solver_total) do
    controllers = collect_outer_controllers(net)
    if isempty(controllers)
      runpf!(net, solver_config; verbose = verbose, pv_table_rows = pv_table_rows, validate_limits_after_pf = validate_limits_after_pf, q_limit_violation_headroom = q_limit_violation_headroom, qlimit_lock_reason = qlimit_lock_reason, performance_profile = performance_profile)
    else
      control_result = run_control!(net; controllers = controllers, pf_config = solver_config, control_config = control_config(), verbose = verbose, performance_profile = performance_profile)
      # The legacy `erg` flag describes numerical PF success/failure.
      # Control-loop terminal states such as :blocked or :max_outer_iterations are
      # reported through `control_result` and `net.control_result`, not as PF failure.
      (control_result.last_pf_iterations, _legacy_erg_from_control_status(control_result.status))
      end
    end
  end
  solver_elapsed_s = begin
    timings = performance_profile isa AbstractDict ? get(performance_profile, :timings, nothing) : nothing
    row = timings isa AbstractDict ? get(timings, :solver_total, nothing) : nothing
    row isa NamedTuple && hasproperty(row, :elapsed_s) ? Float64(getproperty(row, :elapsed_s)) : nothing
  end

  rect_status = method === :rectangular ? rectangular_pf_status(net) : nothing
  pf_status = _build_pf_outcome_status(method, erg, ite, etime, rect_status)
  outcome = pf_status.outcome
  if pf_status.solution_available || printResultAnyCase
    # Calculate network losses and print results
    _perf_profile_time!(performance_profile, :postprocess_losses_and_flows) do
      calcNetLosses!(net)
      calcLinkFlowsKCL!(net)
    end
    jpath = ""
    if show_results
      _perf_profile_time!(performance_profile, :result_output) do
        printACPFlowResults(net, etime, ite, tol, printResultToFile, jpath; converged = (outcome == :converged), solver = method, solver_time_s = solver_elapsed_s, result_mode = result_mode, max_rows = row_limit)
        if outcome == :not_converged
          println("Diagnostic last-iterate voltages; Newton-Raphson did not converge.")
        elseif outcome == :converged_limits_failed
          println("Numerical voltage solution available, but Q-limit / active-set validation failed.")
          println("  final status: ", outcome)
          println("  reason: ", rect_status.reason_text)
        end
      end
    end
  elseif erg == 1
    _print_ac_pf_nonconvergence(method, net)
  else
    @error "Errors during calculation of Newton-Raphson"
  end

  return (ite, erg, etime)
end
