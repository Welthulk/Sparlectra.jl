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

# file: src/acpflow/import_context.jl
# purpose: casefile resolution and import-context helpers for the framework
#          runner: MATPOWER/CGMES import into a Net, config copy helpers, and
#          PV/PQ lock-bus resolution

function _resolve_sparlectra_casefile(casefile::String, path::Union{Nothing,String})::String
  ext = lowercase(splitext(casefile)[2])
  ext in (".m", ".jl") || throw(ArgumentError("run_sparlectra: file extension $(ext) is not supported; use .m or .jl."))
  c = String(strip(casefile))
  if path !== nothing
    filename = joinpath(path, c)
    isfile(filename) || error("File $(filename) not found")
    return filename
  end
  # explicit or cwd-relative path
  isfile(c) && return abspath(c)
  # legacy dev-checkout convention: <cwd>/data/mpower/<case>
  legacy = joinpath(pwd(), "data", "mpower", c)
  isfile(legacy) && return legacy
  # a path-like argument that does not exist must fail explicitly instead of
  # being treated as a downloadable bare case name
  occursin(r"[\\/]", c) && error("File $(c) not found")
  # bare case name: same resolution as run_sparlectra_cases (case cache,
  # download on demand); removes the pwd() dependency for registered installs
  return ensure_casefile(c)
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
  # Field splat instead of an explicit keyword list: the old hand-written list
  # silently dropped `cgmes` and `webui`, resetting them to defaults on every
  # copy (latent — the previous callers never read those fields afterwards).
  fields = NamedTuple{fieldnames(SparlectraConfig)}(getfield.(Ref(cfg), fieldnames(SparlectraConfig)))
  return SparlectraConfig(; fields..., powerflow = powerflow)
end

function _copy_start_current_iteration_with(cfg::StartCurrentIterationConfig; kwargs...)::StartCurrentIterationConfig
  fields = NamedTuple{fieldnames(StartCurrentIterationConfig)}(getfield.(Ref(cfg), fieldnames(StartCurrentIterationConfig)))
  return StartCurrentIterationConfig(; fields..., kwargs...)
end

"""
    _cgmes_start_values_powerflow(pf_cfg, start_values) -> (PowerFlowConfig, overridden::Vector{String})

Apply `cgmes_import.start_values` to the power-flow configuration of one CGMES
run. `:flat` forces a synthetic flat start; `:sv` keeps the imported SvVoltage
state and force-disables every competing start-value machine — the same set the
fixed-reference self-check forces off in `_self_check_forced_overrides`
(run_self_check.jl); keep the two lists aligned. Returns the adjusted config
plus the effective config keys that were overridden (for the decision line in
run.log/cgmes.log). MATPOWER/DTF runs never call this.
"""
function _cgmes_start_values_powerflow(pf_cfg::PowerFlowConfig, start_values::Symbol)
  overridden = String[]
  if start_values === :flat
    pf_cfg.start_mode.flatstart || push!(overridden, "power_flow.flatstart=false")
    return _copy_powerflow_with(pf_cfg; start_mode = _copy_start_mode_with(pf_cfg.start_mode; flatstart = true)), overridden
  end
  start_values === :sv || throw(ArgumentError("cgmes_import.start_values must be one of $(CGMES_START_VALUES_VALUES); got $(start_values)."))
  pf_cfg.start_mode.flatstart && push!(overridden, "power_flow.flatstart=true")
  pf_cfg.start_mode.start_projection && push!(overridden, "power_flow.start_mode.start_projection=true")
  pf_cfg.start_mode.dc_seed_unconditional && push!(overridden, "power_flow.start_mode.dc_seed_unconditional=true")
  pf_cfg.start_current_iteration.enabled && push!(overridden, "power_flow.start_current_iteration.enabled=true")
  pf_cfg.apslf_start.enabled && push!(overridden, "power_flow.apslf_start.enabled=true")
  start_mode = _copy_start_mode_with(pf_cfg.start_mode; flatstart = false, start_projection = false, dc_seed_unconditional = false)
  return _copy_powerflow_with(pf_cfg; start_mode = start_mode, start_current_iteration = _copy_start_current_iteration_with(pf_cfg.start_current_iteration; enabled = false), apslf_start = ApslfStartConfig(enabled = false, order = pf_cfg.apslf_start.order)), overridden
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

  # Opt-in binary net cache (issue #292). Only active with auto_profile off:
  # auto-profile derives config rewrites from the parsed case, which a cache
  # hit would skip.
  cache_active = mat_cfg.net_cache_enabled && mat_cfg.auto_profile === :off
  if mat_cfg.net_cache_enabled && mat_cfg.auto_profile !== :off
    @info "matpower_import.net_cache is enabled but inactive because matpower_import.auto_profile is not 'off'." auto_profile = mat_cfg.auto_profile
  end
  cache_components = nothing
  cache_path = ""
  cached = nothing
  if cache_active
    cached = _perf_profile_time!(performance_profile, :net_cache_lookup) do
      components = _net_cache_components(filename, cfg)
      cache_components = components
      cache_path = _net_cache_path(filename, _net_cache_key(components))
      _net_cache_load(cache_path, components)
    end
    if performance_profile isa AbstractDict
      performance_profile[:net_cache] = Dict{String,Any}("enabled" => true, "hit" => cached !== nothing, "path" => cache_path)
    end
  end

  phase_callback("reading_matpower_case")
  phase_callback(extension == ".jl" ? "loading_julia_case" : "parsing_matpower_file")
  mpc = if cached !== nothing
    @info "MATPOWER net cache hit" casefile = filename cache = cache_path
    cached
  else
    _perf_profile_time!(performance_profile, :matpower_case_parse) do
      MatpowerIO.read_case(filename; legacy_compat = true)
    end
  end
  if cache_active && cached === nothing && cache_components !== nothing
    _perf_profile_time!(performance_profile, :net_cache_store) do
      _net_cache_store(cache_path, cache_components, mpc)
    end
  end
  if mat_cfg.matpower_dcline_mode !== :pf_injections
    try
      MatpowerIO.assert_no_active_dcline(mpc; casefile = filename)
    catch err
      if err isa MatpowerIO.UnsupportedMatpowerDclineError
        details = err.details
        println(stdout, "matpower_dcline_detected")
        println(stdout, "matpower_dcline_unsupported")
        println(stdout, "powerflow_aborted_unsupported_matpower_dcline")
        println(stdout, err.message)
        if performance_profile isa AbstractDict
          performance_profile[:unsupported_matpower_dcline] = details
        end
      end
      rethrow()
    end
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
      tap_changer_model = cfg.transformer.tap_changer_model,
      matpower_pv_voltage_source = mat_cfg.pv_voltage_source,
      matpower_pv_voltage_mismatch_tol_pu = mat_cfg.pv_voltage_mismatch_tol_pu,
      apply_bus_names = mat_cfg.apply_bus_names,
      apply_branch_names = mat_cfg.apply_branch_names,
      apply_branch_kind = mat_cfg.apply_branch_kind,
      import_for001_contingencies = mat_cfg.import_for001_contingencies,
      matpower_dcline_mode = mat_cfg.matpower_dcline_mode,
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