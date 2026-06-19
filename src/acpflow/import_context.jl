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

function _resolve_sparlectra_casefile(casefile::String, path::Union{Nothing,String})::String
  ext = lowercase(splitext(casefile)[2])
  ext in (".m", ".jl") || throw(ArgumentError("run_sparlectra: file extension $(ext) is not supported; use .m or .jl."))
  filename = path === nothing ? joinpath(pwd(), "data", "mpower", strip(casefile)) : joinpath(path, strip(casefile))
  isfile(filename) || error("File $(filename) not found")
  return filename
end

function _resolve_matpower_lock_pv_to_pq_buses(net::Net, buses::AbstractVector{Int}; verbose::Int = 0)::Vector{Int}
  isempty(buses) && return Int[]

  orig_to_net = Dict{Int,Int}()
  sizehint!(orig_to_net, length(net.busOrigIdxDict))
  for (net_idx, orig_idx) in net.busOrigIdxDict
    orig_to_net[Int(orig_idx)] = Int(net_idx)
  end

  resolved = Int[]
  nb = length(net.nodeVec)
  for bus in buses
    b = Int(bus)
    if haskey(orig_to_net, b)
      push!(resolved, orig_to_net[b])
    elseif 1 <= b <= nb
      push!(resolved, b)
    else
      verbose > 0 && @warn "Configured MATPOWER PV/PQ lock bus was not found and will be ignored" bus = b
    end
  end
  unique!(resolved)
  sort!(resolved)
  return resolved
end

function _copy_qlimits_with(ql::QLimitConfig; kwargs...)::QLimitConfig
  fields = NamedTuple{fieldnames(QLimitConfig)}(getfield.(Ref(ql), fieldnames(QLimitConfig)))
  return QLimitConfig(; fields..., kwargs...)
end

function _copy_start_mode_with(start::StartModeConfig; kwargs...)::StartModeConfig
  fields = NamedTuple{fieldnames(StartModeConfig)}(getfield.(Ref(start), fieldnames(StartModeConfig)))
  return StartModeConfig(; fields..., kwargs...)
end

function _copy_powerflow_with(pf::PowerFlowConfig; kwargs...)::PowerFlowConfig
  fields = NamedTuple{fieldnames(PowerFlowConfig)}(getfield.(Ref(pf), fieldnames(PowerFlowConfig)))
  return PowerFlowConfig(; fields..., kwargs...)
end

function _copy_sparlectra_with_powerflow(cfg::SparlectraConfig, powerflow::PowerFlowConfig)::SparlectraConfig
  return SparlectraConfig(; powerflow = powerflow, state_estimation = cfg.state_estimation, matpower = cfg.matpower, performance = cfg.performance, benchmark = cfg.benchmark, runtime = cfg.runtime, diagnostics = cfg.diagnostics, output = cfg.output, control = cfg.control)
end

function _resolve_matpower_powerflow_ids_after_import(net::Net, cfg::SparlectraConfig; verbose::Int = 0)::SparlectraConfig
  qlimits = cfg.powerflow.qlimits
  resolved = _resolve_matpower_lock_pv_to_pq_buses(net, qlimits.lock_pv_to_pq_buses; verbose = verbose)
  resolved == qlimits.lock_pv_to_pq_buses && return cfg
  qlimits2 = _copy_qlimits_with(qlimits; lock_pv_to_pq_buses = resolved)
  return _copy_sparlectra_with_powerflow(cfg, _copy_powerflow_with(cfg.powerflow; qlimits = qlimits2))
end

function _uses_projected_matpower_start(start::StartModeConfig)::Bool
  return start.flatstart && (start.voltage_mode in (:all_bus_vm, :profile_blend) || start.angle_mode in (:bus_va_blend, :matpower_va))
end

function _copy_sparlectra_with_projected_matpower_start(cfg::SparlectraConfig)::SparlectraConfig
  start_mode = _copy_start_mode_with(cfg.powerflow.start_mode; flatstart = false)
  return _copy_sparlectra_with_powerflow(cfg, _copy_powerflow_with(cfg.powerflow; start_mode = start_mode))
end

function _import_sparlectra_context(casefile::String, path::Union{Nothing,String}, cfg::SparlectraConfig; performance_profile = nothing)
  filename = _resolve_sparlectra_casefile(casefile, path)
  pf_cfg = cfg.powerflow
  mat_cfg = cfg.matpower
  phase_callback = performance_profile isa AbstractDict ? get(performance_profile, :phase_callback, phase -> nothing) : phase -> nothing
  extension = lowercase(splitext(filename)[2])
  phase_callback("reading_matpower_case")
  phase_callback(extension == ".jl" ? "loading_julia_case" : "parsing_matpower_file")
  mpc = _perf_profile_time!(performance_profile, :matpower_case_parse) do
    MatpowerIO.read_case(filename; legacy_compat = true)
  end
  auto_profile_result = nothing
  println(stdout, "Runtime casefile: ", filename)
  print_matpower_import_runtime_options(stdout, "Original MATPOWER import options", cfg)
  if cfg.matpower.auto_profile !== :off
    phase_callback("matpower_auto_profile")
    auto_profile_result = _perf_profile_time!(performance_profile, :matpower_auto_profile) do
      run_matpower_import_auto_profile(mpc, cfg)
    end
    cfg = auto_profile_result.config
    pf_cfg = cfg.powerflow
    mat_cfg = cfg.matpower
    if performance_profile isa AbstractDict
      performance_profile[:matpower_auto_profile_result] = auto_profile_result
      performance_profile[:matpower_auto_profile_casefile] = filename
    end
    if mat_cfg.auto_profile_log
      write_matpower_import_auto_profile(stdout, auto_profile_result, cfg; casefile = filename)
    end
  else
    print_matpower_import_runtime_options(stdout, "Final effective MATPOWER import options", cfg)
  end
  phase_callback("building_sparlectra_net")
  net = _perf_profile_time!(performance_profile, :network_construction) do
    createNetFromMatPowerCase(
      mpc = mpc,
      flatstart = pf_cfg.start_mode.flatstart,
      cooldown = pf_cfg.qlimits.cooldown_iters,
      q_hyst_pu = pf_cfg.qlimits.hysteresis_pu,
      enable_pq_gen_controllers = mat_cfg.enable_pq_gen_controllers,
      bus_shunt_model = mat_cfg.bus_shunt_model,
      matpower_shift_sign = mat_cfg.shift_sign,
      matpower_shift_unit = mat_cfg.shift_unit,
      matpower_ratio = mat_cfg.ratio,
      matpower_pv_voltage_source = mat_cfg.pv_voltage_source,
      matpower_pv_voltage_mismatch_tol_pu = mat_cfg.pv_voltage_mismatch_tol_pu,
      preallocate_network = mat_cfg.preallocate_network,
      preallocate_min_buses = mat_cfg.preallocate_min_buses,
      profile = performance_profile,
    )
  end
  phase_callback("applying_import_options")
  MatpowerIO.apply_mp_isolated_buses!(net, mpc; verbose = 0)
  MatpowerIO.apply_mp_bus_vmva_init!(net, mpc; flatstart = pf_cfg.start_mode.flatstart, verbose = 0)
  projected_start_applied = false
  if pf_cfg.start_mode.flatstart && (pf_cfg.start_mode.voltage_mode != :classic || pf_cfg.start_mode.angle_mode != :classic)
    phase_callback("preparing_start_values")
    _apply_matpower_start_modes!(net, mpc, pf_cfg.start_mode, mat_cfg; performance_profile = performance_profile)
    if _uses_projected_matpower_start(pf_cfg.start_mode)
      net.flatstart = false
      projected_start_applied = true
    end
  end
  run_cfg = projected_start_applied ? _copy_sparlectra_with_projected_matpower_start(cfg) : cfg
  projected_start_applied && @debug "MATPOWER projected start applied; effective solver flatstart disabled for this run."
  return (net = net, config = run_cfg, projected_start_applied = projected_start_applied, auto_profile_result = auto_profile_result)
end

function _import_sparlectra_net(casefile::String, path::Union{Nothing,String}, cfg::SparlectraConfig; performance_profile = nothing)::Net
  return _import_sparlectra_context(casefile, path, cfg; performance_profile = performance_profile).net
end
