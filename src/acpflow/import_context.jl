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
  # .json is the Sparlectra Case Format (#342), handled by its own reader
  ext in (".m", ".jl", ".json") || throw(ArgumentError("run_sparlectra: file extension $(ext) is not supported; use .m, .jl, or .json (Sparlectra Case Format)."))
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

# flat dotted difference of two typed configurations over the sections the
# auto profile may rewrite: the D11 override surface and its log lines
function _config_struct_diff!(out::Dict{String,Any}, prefix::String, a, b)
  for f in fieldnames(typeof(a))
    va = getfield(a, f)
    vb = getfield(b, f)
    path = string(prefix, ".", f)
    if va isa Union{Symbol,AbstractString,Number,Bool,Nothing,AbstractVector,AbstractDict,Tuple} || typeof(va) != typeof(vb)
      isequal(va, vb) || (out[path] = vb isa Symbol ? String(vb) : vb)
    else
      _config_struct_diff!(out, path, va, vb)
    end
  end
  return out
end

function _config_flat_difference(a::SparlectraConfig, b::SparlectraConfig)::Dict{String,Any}
  out = Dict{String,Any}()
  _config_struct_diff!(out, "power_flow", a.powerflow, b.powerflow)
  _config_struct_diff!(out, "matpower_import", a.matpower, b.matpower)
  _config_struct_diff!(out, "model", a.model, b.model)
  return out
end

function _import_sparlectra_context(casefile::AbstractString, path::Union{Nothing,AbstractString}, cfg::SparlectraConfig; performance_profile = nothing)
  filename = _resolve_sparlectra_casefile(casefile, path)
  pf_cfg = cfg.powerflow
  mat_cfg = cfg.matpower
  model_cfg = cfg.model
  phase_callback = performance_profile isa AbstractDict ? get(performance_profile, :phase_callback, phase -> nothing) : phase -> nothing
  extension = lowercase(splitext(filename)[2])

  # Sparlectra Case Format (#342): a self-describing case file that carries
  # its own model, so it bypasses the MATPOWER parse chain entirely. The
  # reader builds through the same constructors and validates itself.
  if extension == ".json"
    phase_callback("reading_scf_case")
    # build_net applies the net parameters and the shunt-model default with
    # the RESOLVED run configuration (D12, corrected by task_import_direct:
    # exactly once PER IMPORTER; for SCF that once lives inside build_net)
    scfcase = _perf_profile_time!(performance_profile, :scf_case_read) do
      read_scf_json(filename)
    end
    net = _perf_profile_time!(performance_profile, :scf_case_build) do
      build_net(scfcase; config = cfg)
    end
    return (net = net, config = cfg, projected_start_applied = false, auto_profile_result = nothing, auto_profile_overrides = Dict{String,Any}())
  end

  # task_import_direct (maintainer decision 2026-09-04): the net cache
  # stored the CONVERTED SCFCase, the very side product the direct import
  # no longer creates, so the cache fell away with the conversion. The key
  # stays readable and says so out loud instead of silently doing nothing;
  # the removal follow-up is recorded with issue #292 in the task report.
  if model_cfg.net_cache_enabled
    @warn "model.net_cache_enabled is inert since the direct import (task_import_direct): the import no longer produces the converted case the cache stored." casefile = filename
  end

  phase_callback("reading_matpower_case")
  phase_callback(extension == ".jl" ? "loading_julia_case" : "parsing_matpower_file")
  mpc = _perf_profile_time!(performance_profile, :matpower_case_parse) do
    MatpowerIO.read_case(filename; legacy_compat = true)
  end
  if mpc !== nothing && mat_cfg.matpower_dcline_mode !== :pf_injections
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
  auto_profile_overrides = Dict{String,Any}()
  # The importer's option dump is DIAGNOSTIC output, so it follows
  # output.console_auto_profile like the analysis below: at :full the whole
  # list, otherwise one line naming the case. Printing twelve MATPOWER
  # option lines before every run made the tool read like a MATPOWER
  # front end (maintainer, 2026-09-05).
  if cfg.output.console_auto_profile === :full
    println(stdout, "Runtime casefile: ", filename)
    print_matpower_import_runtime_options(stdout, "Original MATPOWER import options", cfg)
  else
    println(stdout, "Case: ", filename)
  end
  if cfg.model.auto_profile !== :off
    phase_callback("matpower_auto_profile")
    auto_profile_result = _perf_profile_time!(performance_profile, :matpower_auto_profile) do
      run_matpower_import_auto_profile(mpc, cfg)
    end
    # D11 surface: the profile's decisions as dotted overrides, computed as
    # the flat difference between the incoming and the rewritten effective
    # configuration; carried on the imported case for the services
    auto_profile_overrides = _config_flat_difference(cfg, auto_profile_result.config)
    for (k, v) in sort!(collect(auto_profile_overrides); by = first)
      println(stdout, "auto_profile: ", k, " = ", v)
    end
    cfg = auto_profile_result.config
    pf_cfg = cfg.powerflow
    mat_cfg = cfg.matpower
    model_cfg = cfg.model
    if performance_profile isa AbstractDict
      performance_profile[:matpower_auto_profile_result] = auto_profile_result
      performance_profile[:matpower_auto_profile_casefile] = filename
    end
    if model_cfg.auto_profile_log
      write_matpower_import_auto_profile(stdout, auto_profile_result, cfg; casefile = filename)
    end
  elseif cfg.output.console_auto_profile === :full
    print_matpower_import_runtime_options(stdout, "Final effective MATPOWER import options", cfg)
  elseif cfg.output.console_auto_profile !== :off
    # auto-profile is off, so nothing analyses the conventions: name the
    # three that actually move results, on ONE line. The full option list
    # stays behind :full and in the effective_config.yaml artifact.
    println(stdout, "Import conventions: ratio=", cfg.matpower.ratio, ", shift=", cfg.matpower.shift_unit,
            " (sign ", cfg.matpower.shift_sign, "), bus shunt=", cfg.model.bus_shunt_model)
  end
  phase_callback("building_sparlectra_net")
  # task_import_direct: the MATPOWER importer builds its network DIRECTLY;
  # no intermediate SCFCase, no re-encoding, no loss. Converting into SCF
  # is an explicit action elsewhere (convert_case), never a step in here.
  net = _perf_profile_time!(performance_profile, :network_construction) do
    createNetFromMatPowerCase(
      mpc = mpc,
      log = false,
      flatstart = pf_cfg.start_mode.flatstart,
      cooldown = pf_cfg.qlimits.cooldown_iters,
      q_hyst_pu = pf_cfg.qlimits.hysteresis_pu,
      enable_pq_gen_controllers = mat_cfg.enable_pq_gen_controllers,
      bus_shunt_model = model_cfg.bus_shunt_model,
      matpower_shift_sign = mat_cfg.shift_sign,
      matpower_shift_unit = mat_cfg.shift_unit,
      matpower_ratio = mat_cfg.ratio,
      tap_changer_model = model_cfg.tap_changer_model,
      matpower_pv_voltage_source = mat_cfg.pv_voltage_source,
      matpower_pv_voltage_mismatch_tol_pu = mat_cfg.pv_voltage_mismatch_tol_pu,
      apply_bus_names = mat_cfg.apply_bus_names,
      apply_branch_names = mat_cfg.apply_branch_names,
      apply_branch_kind = mat_cfg.apply_branch_kind,
      import_for001_contingencies = mat_cfg.import_for001_contingencies,
      matpower_dcline_mode = mat_cfg.matpower_dcline_mode,
      preallocate_network = model_cfg.preallocate_network,
      preallocate_min_buses = model_cfg.preallocate_min_buses,
    )
  end
  MatpowerIO.apply_mp_isolated_buses!(net, mpc)
  MatpowerIO.apply_mp_bus_vmva_init!(net, mpc; flatstart = pf_cfg.start_mode.flatstart)
  # D12, corrected by task_import_direct: the net-parameter stamping happens
  # exactly once PER IMPORTER, and this is the MATPOWER importer's once. Do
  # not pull the four format call sites back together; each importer owns
  # its stamping where it finishes.
  _apply_config_net_parameters!(net, cfg)
  net.bus_shunt_model = normalize_bus_shunt_model(model_cfg.bus_shunt_model)
  if performance_profile isa AbstractDict
    # the construction-report shape (counts plus the per-stage placeholder
    # block the performance log prints); the real timing is
    # :network_construction of the direct importer
    performance_profile[:network_construction_nbus] = length(net.nodeVec)
    performance_profile[:network_construction_nbranch] = length(net.branchVec)
    performance_profile[:network_construction_ngen] = count(ps -> isGenerator(ps), net.prosumpsVec)
    performance_profile[:network_construction_subtimings] = Dict(
      :matpower_data_normalization => 0.0,
      :net_object_creation => 0.0,
      :bus_import => 0.0,
      :branch_import => 0.0,
      :generator_prosumer_import => 0.0,
      :shunt_import => 0.0,
      :pq_generator_controller_setup => 0.0,
      :bus_branch_dictionary_construction => 0.0,
      :validation_post_import_consistency => 0.0,
    )
  end
  # the flat-start decision is the RUN's, not the case's (D10)
  net.flatstart = pf_cfg.start_mode.flatstart
  phase_callback("applying_import_options")
  projected_start_applied = false
  if pf_cfg.start_mode.flatstart && (pf_cfg.start_mode.voltage_mode != :classic || pf_cfg.start_mode.angle_mode != :classic)
    phase_callback("preparing_start_values")
    _apply_start_modes!(net, _matpower_raw_starts(mpc), pf_cfg.start_mode; performance_profile = performance_profile)
    if _uses_projected_matpower_start(pf_cfg.start_mode)
      net.flatstart = false
      projected_start_applied = true
    end
  end
  run_cfg = projected_start_applied ? _copy_sparlectra_with_projected_matpower_start(cfg) : cfg
  projected_start_applied && @debug "MATPOWER projected start applied; effective solver flatstart disabled for this run."
  return (net = net, config = run_cfg, projected_start_applied = projected_start_applied, auto_profile_result = auto_profile_result, auto_profile_overrides = auto_profile_overrides)
end

# AbstractString on the outside, String on the inside: the service entry
# points hand paths through as AbstractString, and a SubString must not turn
# a run into a MethodError (it did, on a case file carrying its own config).
function _import_sparlectra_net(casefile::AbstractString, path::Union{Nothing,AbstractString}, cfg::SparlectraConfig; performance_profile = nothing)::Net
  return _import_sparlectra_context(String(casefile), path === nothing ? nothing : String(path), cfg; performance_profile = performance_profile).net
end