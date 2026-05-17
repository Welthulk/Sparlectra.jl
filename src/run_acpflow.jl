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
Function to perform AC power flow analysis.

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
function _apply_matpower_flatstart_modes!(net::Net, mpc; voltage_mode = :classic, angle_mode = :classic, matpower_pv_voltage_source = :gen_vg, matpower_pv_voltage_mismatch_tol_pu::Float64 = 1e-4)
  vmode = Symbol(voltage_mode)
  amode = Symbol(angle_mode)
  vmode in (:classic, :pv_gen_vg, :pv_bus_vm, :all_bus_vm, :bus_vm_va_blend) || error("flatstart_voltage_mode must be classic, pv_gen_vg, pv_bus_vm, all_bus_vm, or bus_vm_va_blend")
  amode in (:classic, :dc, :bus_va_blend, :matpower_va) || error("flatstart_angle_mode must be classic, dc, bus_va_blend, or matpower_va")
  busrow = MatpowerIO.bus_row_index(mpc)
  pv_rows = MatpowerIO.pv_voltage_reference_rows(mpc; matpower_pv_voltage_source = matpower_pv_voltage_source, tol = matpower_pv_voltage_mismatch_tol_pu, warn = false)
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
    elseif vmode == :bus_vm_va_blend
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

function run_acpflow(;
  max_ite::Int = 30,
  tol::Float64 = 1e-6,
  casefile::String,
  path::Union{Nothing,String} = nothing,
  verbose::Int = 0,  
  printResultToFile::Bool = false,
  printResultAnyCase::Bool = false,
  opt_fd::Bool = false,
  opt_sparse::Bool = false,
  method::Symbol = :rectangular,
  autodamp::Bool = false,
  autodamp_min::Float64 = 1e-3, start_projection::Bool = false, start_projection_try_dc_start::Bool = true, start_projection_try_blend_scan::Bool = true, start_projection_blend_lambdas::AbstractVector{<:Real} = [0.25, 0.5, 0.75], start_projection_dc_angle_limit_deg::Float64 = 60.0, qlimit_start_iter::Int = 2, qlimit_start_mode::Symbol = :iteration, qlimit_auto_q_delta_pu::Float64 = 1e-4,
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
  opt_flatstart::Bool = false,
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
  enable_pq_gen_controllers::Bool = true,
  bus_shunt_model = :admittance,
  matpower_shift_sign::Real = 1.0,
  matpower_shift_unit = :deg,
  matpower_ratio = :normal,
  reference_vm_pu::Union{Nothing,Float64} = nothing,
  reference_va_deg::Union{Nothing,Float64} = nothing,
  matpower_pv_voltage_source = :gen_vg,
  matpower_pv_voltage_mismatch_tol_pu::Float64 = 1e-4,
  flatstart_voltage_mode = :classic,
  flatstart_angle_mode = :classic,
  performance_profile = nothing,
)
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
      createNetFromMatPowerFile(filename = in_path, log = (verbose > 0), flatstart = opt_flatstart, cooldown = cooldown_iters, q_hyst_pu = q_hyst_pu, enable_pq_gen_controllers = enable_pq_gen_controllers, bus_shunt_model = bus_shunt_model, matpower_shift_sign = matpower_shift_sign, matpower_shift_unit = matpower_shift_unit, matpower_ratio = matpower_ratio, reference_vm_pu = reference_vm_pu, reference_va_deg = reference_va_deg, matpower_pv_voltage_source = matpower_pv_voltage_source, matpower_pv_voltage_mismatch_tol_pu = matpower_pv_voltage_mismatch_tol_pu)
    end
    if opt_flatstart && (Symbol(flatstart_voltage_mode) != :classic || Symbol(flatstart_angle_mode) != :classic)
      mpc_init = _perf_profile_time!(performance_profile, :matpower_parse_for_flatstart) do
        MatpowerIO.read_case(in_path; legacy_compat = true)
      end
      _apply_matpower_flatstart_modes!(myNet, mpc_init; voltage_mode = flatstart_voltage_mode, angle_mode = flatstart_angle_mode, matpower_pv_voltage_source = matpower_pv_voltage_source, matpower_pv_voltage_mismatch_tol_pu = matpower_pv_voltage_mismatch_tol_pu)
      if Symbol(flatstart_voltage_mode) in (:all_bus_vm, :bus_vm_va_blend) || Symbol(flatstart_angle_mode) in (:bus_va_blend, :matpower_va)
        opt_flatstart = false
      end
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

  # Run power flow / optional tap-controller outer loop
  ite = 0

  erg = 2
  etime = @elapsed begin
    ite, erg = _perf_profile_time!(performance_profile, :solver_total) do
    if isempty(_tap_controllers(myNet))
      runpf!(myNet, max_ite, tol, verbose; opt_fd = opt_fd, opt_sparse = opt_sparse, method = method, autodamp = autodamp, autodamp_min = autodamp_min, start_projection = start_projection, start_projection_try_dc_start = start_projection_try_dc_start, start_projection_try_blend_scan = start_projection_try_blend_scan, start_projection_blend_lambdas = start_projection_blend_lambdas, start_projection_dc_angle_limit_deg = start_projection_dc_angle_limit_deg, qlimit_start_iter = qlimit_start_iter, qlimit_start_mode = qlimit_start_mode, qlimit_auto_q_delta_pu = qlimit_auto_q_delta_pu, opt_flatstart = opt_flatstart, pv_table_rows = pv_table_rows, validate_limits_after_pf = validate_limits_after_pf, q_limit_violation_headroom = q_limit_violation_headroom, lock_pv_to_pq_buses = lock_pv_to_pq_buses_resolved, qlimit_trace_buses = qlimit_trace_buses, qlimit_lock_reason = qlimit_lock_reason, qlimit_guard = qlimit_guard, qlimit_guard_min_q_range_pu = qlimit_guard_min_q_range_pu, qlimit_guard_zero_range_mode = qlimit_guard_zero_range_mode, qlimit_guard_narrow_range_mode = qlimit_guard_narrow_range_mode, qlimit_guard_log = qlimit_guard_log, qlimit_guard_max_switches = qlimit_guard_max_switches, qlimit_guard_accept_bounded_violations = qlimit_guard_accept_bounded_violations, qlimit_guard_max_remaining_violations = qlimit_guard_max_remaining_violations, qlimit_guard_freeze_after_repeated_switching = qlimit_guard_freeze_after_repeated_switching, qlimit_guard_violation_mode = qlimit_guard_violation_mode, qlimit_guard_violation_threshold_pu = qlimit_guard_violation_threshold_pu, performance_profile = performance_profile)
      else
      run_tap_controllers_outer!(myNet; max_ite = max_ite, tol = tol, verbose = verbose, opt_fd = opt_fd, opt_sparse = opt_sparse, method = method, autodamp = autodamp, autodamp_min = autodamp_min, start_projection = start_projection, start_projection_try_dc_start = start_projection_try_dc_start, start_projection_try_blend_scan = start_projection_try_blend_scan, start_projection_blend_lambdas = start_projection_blend_lambdas, start_projection_dc_angle_limit_deg = start_projection_dc_angle_limit_deg, qlimit_start_iter = qlimit_start_iter, qlimit_start_mode = qlimit_start_mode, qlimit_auto_q_delta_pu = qlimit_auto_q_delta_pu, opt_flatstart = opt_flatstart, pv_table_rows = pv_table_rows, validate_limits_after_pf = validate_limits_after_pf, q_limit_violation_headroom = q_limit_violation_headroom, lock_pv_to_pq_buses = lock_pv_to_pq_buses_resolved, qlimit_trace_buses = qlimit_trace_buses, qlimit_lock_reason = qlimit_lock_reason, qlimit_guard = qlimit_guard, qlimit_guard_min_q_range_pu = qlimit_guard_min_q_range_pu, qlimit_guard_zero_range_mode = qlimit_guard_zero_range_mode, qlimit_guard_narrow_range_mode = qlimit_guard_narrow_range_mode, qlimit_guard_log = qlimit_guard_log, qlimit_guard_max_switches = qlimit_guard_max_switches, qlimit_guard_accept_bounded_violations = qlimit_guard_accept_bounded_violations, qlimit_guard_max_remaining_violations = qlimit_guard_max_remaining_violations, qlimit_guard_freeze_after_repeated_switching = qlimit_guard_freeze_after_repeated_switching, qlimit_guard_violation_mode = qlimit_guard_violation_mode, qlimit_guard_violation_threshold_pu = qlimit_guard_violation_threshold_pu)
      end
    end
  end

  if show_compact_result
    pv_to_pq_events = length(myNet.qLimitLog)
    pv_to_pq_buses = length(myNet.qLimitEvents)
    rect_status = method === :rectangular ? rectangular_pf_status(myNet) : nothing
    if rect_status !== nothing
      @printf("method=%-12s  numerical_solution=%s  q_limit_active_set=%s  final_converged=%s  status=%s  iterations=%d  pv2pq_events=%d  pv2pq_buses=%d  time=%8.6f s  reason=%s", String(method), rect_status.numerical_converged ? "OK" : "FAIL", rect_status.q_limit_active_set_ok ? "OK" : "FAIL", string(rect_status.final_converged), String(rect_status.status), ite, pv_to_pq_events, pv_to_pq_buses, etime, rect_status.reason_text)
    else
      converged_text = erg == 0 ? "yes" : erg == 1 ? "no" : "error"
      @printf("method=%-12s  converged=%s  iterations=%d  pv2pq_events=%d  pv2pq_buses=%d  time=%8.6f s", String(method), converged_text, ite, pv_to_pq_events, pv_to_pq_buses, etime)
    end
    println()
  end

  if status_ref !== nothing
    rect_status = method === :rectangular ? rectangular_pf_status(myNet) : nothing
    if rect_status === nothing
      status_ref[] = (converged = (erg == 0), erg = erg, iterations = ite, elapsed_s = etime, method = method)
    else
      status_ref[] = (converged = (erg == 0), erg = erg, iterations = ite, elapsed_s = etime, method = method,
                      numerical_converged = rect_status.numerical_converged,
                      q_limit_active_set_ok = rect_status.q_limit_active_set_ok,
                      final_converged = rect_status.final_converged,
                      status = rect_status.status,
                      nr_converged = rect_status.nr_converged,
                      active_set_converged = rect_status.active_set_converged,
                      reason = rect_status.reason,
                      reason_text = rect_status.reason_text,
                      pv_q_limit_violations = rect_status.pv_q_limit_violations,
                      ref_q_limit_violations = rect_status.ref_q_limit_violations,
                      final_pv_voltage_residual = rect_status.final_pv_voltage_residual)
    end
  end
  if erg == 0 || printResultAnyCase
    # Calculate network losses and print results
    _perf_profile_time!(performance_profile, :postprocess_losses_and_flows) do
      calcNetLosses!(myNet)
      calcLinkFlowsKCL!(myNet)
    end
    jpath = printResultToFile ? out_path : ""
    if show_results || printResultAnyCase
      _perf_profile_time!(performance_profile, :result_output) do
        printACPFlowResults(myNet, etime, ite, tol, printResultToFile, jpath; converged = (erg == 0), solver = method)
      end
    end
  elseif erg == 1
    _print_ac_pf_nonconvergence(method, myNet)
  else
    @error "Errors during calculation of Newton-Raphson"
  end

  return myNet
end

"""
Function to perform AC power flow analysis.

Parameters:
- net: Net, the network object.
- max_ite: Int, the maximum number of iterations for the power flow algorithm (default: 30).
- tol: Float64, tolerance for convergence criterion (default: 1e-6).
- verbose: Int, verbosity level for output (default: 0).
- printResultToFile: Bool, flag to print results to a file (default: false).
- printResultAnyCase: Bool, flag to print results even if the power flow fails (default: false).
- autodamp: Bool, enable residual-based Newton step backtracking for `method = :rectangular`.
- autodamp_min: Float64, minimum trial step length for automatic damping.
"""
function run_net_acpflow(; net::Net, max_ite::Int = 30, tol::Float64 = 1e-6, verbose::Int = 0, printResultToFile::Bool = false, printResultAnyCase::Bool = false, opt_fd::Bool = false, opt_sparse::Bool = false, method::Symbol = :rectangular, autodamp::Bool = false, autodamp_min::Float64 = 1e-3, start_projection::Bool = false, start_projection_try_dc_start::Bool = true, start_projection_try_blend_scan::Bool = true, start_projection_blend_lambdas::AbstractVector{<:Real} = [0.25, 0.5, 0.75], start_projection_dc_angle_limit_deg::Float64 = 60.0, qlimit_start_iter::Int = 2, qlimit_start_mode::Symbol = :iteration, qlimit_auto_q_delta_pu::Float64 = 1e-4, show_results::Bool = true, lock_pv_to_pq_buses::AbstractVector{Int} = Int[], opt_flatstart::Bool = true, pv_table_rows::Int = 30, validate_limits_after_pf::Bool = false, q_limit_violation_headroom::Float64 = 0.0, qlimit_trace_buses::AbstractVector{Int} = Int[], qlimit_lock_reason::Symbol = :manual, qlimit_guard::Bool = false, qlimit_guard_min_q_range_pu::Float64 = 1e-4, qlimit_guard_zero_range_mode::Symbol = :lock_pq, qlimit_guard_narrow_range_mode::Symbol = :prefer_pq, qlimit_guard_log::Bool = true, qlimit_guard_max_switches::Int = 10, qlimit_guard_accept_bounded_violations::Bool = false, qlimit_guard_max_remaining_violations::Int = 0, qlimit_guard_freeze_after_repeated_switching::Bool = true, qlimit_guard_violation_mode::Symbol = :delayed_switch, qlimit_guard_violation_threshold_pu::Float64 = 1e-4, bus_shunt_model = net.bus_shunt_model, performance_profile = nothing)

  requested_shunt_model = normalize_bus_shunt_model(bus_shunt_model)
  requested_shunt_model == net.bus_shunt_model || error("run_net_acpflow: bus_shunt_model must be set when constructing/importing the Net; got $(requested_shunt_model) for a net configured as $(net.bus_shunt_model).")

  # Run power flow / optional tap-controller outer loop
  ite = 0

  erg = 2
  etime = @elapsed begin
    ite, erg = _perf_profile_time!(performance_profile, :solver_total) do
    if isempty(_tap_controllers(net))
      runpf!(net, max_ite, tol, verbose; opt_fd = opt_fd, opt_sparse = opt_sparse, method = method, autodamp = autodamp, autodamp_min = autodamp_min, start_projection = start_projection, start_projection_try_dc_start = start_projection_try_dc_start, start_projection_try_blend_scan = start_projection_try_blend_scan, start_projection_blend_lambdas = start_projection_blend_lambdas, start_projection_dc_angle_limit_deg = start_projection_dc_angle_limit_deg, qlimit_start_iter = qlimit_start_iter, qlimit_start_mode = qlimit_start_mode, qlimit_auto_q_delta_pu = qlimit_auto_q_delta_pu, lock_pv_to_pq_buses = lock_pv_to_pq_buses, opt_flatstart = opt_flatstart, pv_table_rows = pv_table_rows, validate_limits_after_pf = validate_limits_after_pf, q_limit_violation_headroom = q_limit_violation_headroom, qlimit_trace_buses = qlimit_trace_buses, qlimit_lock_reason = qlimit_lock_reason, qlimit_guard = qlimit_guard, qlimit_guard_min_q_range_pu = qlimit_guard_min_q_range_pu, qlimit_guard_zero_range_mode = qlimit_guard_zero_range_mode, qlimit_guard_narrow_range_mode = qlimit_guard_narrow_range_mode, qlimit_guard_log = qlimit_guard_log, qlimit_guard_max_switches = qlimit_guard_max_switches, qlimit_guard_accept_bounded_violations = qlimit_guard_accept_bounded_violations, qlimit_guard_max_remaining_violations = qlimit_guard_max_remaining_violations, qlimit_guard_freeze_after_repeated_switching = qlimit_guard_freeze_after_repeated_switching, qlimit_guard_violation_mode = qlimit_guard_violation_mode, qlimit_guard_violation_threshold_pu = qlimit_guard_violation_threshold_pu, performance_profile = performance_profile)
    else
      run_tap_controllers_outer!(net; max_ite = max_ite, tol = tol, verbose = verbose, opt_fd = opt_fd, opt_sparse = opt_sparse, method = method, autodamp = autodamp, autodamp_min = autodamp_min, start_projection = start_projection, start_projection_try_dc_start = start_projection_try_dc_start, start_projection_try_blend_scan = start_projection_try_blend_scan, start_projection_blend_lambdas = start_projection_blend_lambdas, start_projection_dc_angle_limit_deg = start_projection_dc_angle_limit_deg, qlimit_start_iter = qlimit_start_iter, qlimit_start_mode = qlimit_start_mode, qlimit_auto_q_delta_pu = qlimit_auto_q_delta_pu, opt_flatstart = opt_flatstart, pv_table_rows = pv_table_rows, validate_limits_after_pf = validate_limits_after_pf, q_limit_violation_headroom = q_limit_violation_headroom, lock_pv_to_pq_buses = lock_pv_to_pq_buses, qlimit_trace_buses = qlimit_trace_buses, qlimit_lock_reason = qlimit_lock_reason, qlimit_guard = qlimit_guard, qlimit_guard_min_q_range_pu = qlimit_guard_min_q_range_pu, qlimit_guard_zero_range_mode = qlimit_guard_zero_range_mode, qlimit_guard_narrow_range_mode = qlimit_guard_narrow_range_mode, qlimit_guard_log = qlimit_guard_log, qlimit_guard_max_switches = qlimit_guard_max_switches, qlimit_guard_accept_bounded_violations = qlimit_guard_accept_bounded_violations, qlimit_guard_max_remaining_violations = qlimit_guard_max_remaining_violations, qlimit_guard_freeze_after_repeated_switching = qlimit_guard_freeze_after_repeated_switching, qlimit_guard_violation_mode = qlimit_guard_violation_mode, qlimit_guard_violation_threshold_pu = qlimit_guard_violation_threshold_pu)
      end
    end
  end

  if erg == 0 || printResultAnyCase
    # Calculate network losses and print results
    _perf_profile_time!(performance_profile, :postprocess_losses_and_flows) do
      calcNetLosses!(net)
      calcLinkFlowsKCL!(net)
    end
    jpath = ""
    if show_results
      _perf_profile_time!(performance_profile, :result_output) do
        printACPFlowResults(net, etime, ite, tol, printResultToFile, jpath; converged = (erg == 0), solver = method)
      end
    end
  elseif erg == 1
    _print_ac_pf_nonconvergence(method, net)
  else
    @error "Errors during calculation of Newton-Raphson"
  end

  return (ite, erg, etime)
end
