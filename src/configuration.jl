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

# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 20.5.2026
# file: src/configuration.jl

using SHA

"""
    StartModeConfig

Typed power-flow start-option configuration. These fields collect the
flat-start and rectangular start-projection controls that otherwise tend to be
forwarded as long keyword lists through example and benchmark call chains.
"""
Base.@kwdef struct StartModeConfig
  flatstart::Bool = false
  angle_mode::Symbol = :dc
  voltage_mode::Symbol = :profile_blend
  profile_source::Symbol = :matpower_reference
  start_projection::Bool = false
  try_dc_start::Bool = true
  try_blend_scan::Bool = true
  branch_guard::Bool = true
  measure_candidates::Bool = true
  accept_unmeasured_dc_start::Bool = false
  reuse_import_data::Bool = true
  blend_lambdas::Vector{Float64} = [0.25, 0.5, 0.75]
  dc_angle_limit_deg::Float64 = 60.0
end

"""
    StartCurrentIterationConfig

Guarded fixed-point current-injection pre-solve configuration. The stage is
disabled by default and, when enabled, only prepares the initial voltage profile
before the normal rectangular Newton-Raphson solve.
"""
Base.@kwdef struct StartCurrentIterationConfig
  enabled::Bool = false
  max_iter::Int = 10
  tol::Float64 = 1.0e-3
  damping::Float64 = 0.5
  accept_only_if_improved::Bool = true
  min_improvement_factor::Float64 = 0.98
  vm_min_pu::Float64 = 0.5
  vm_max_pu::Float64 = 1.5
  max_angle_step_deg::Float64 = 30.0
  only_for_large_cases::Bool = false
end

"""
    QLimitConfig

Typed reactive-power limit switching configuration used by power-flow runners.
"""
Base.@kwdef struct QLimitConfig
  start_iter::Int = 2
  start_mode::Symbol = :iteration
  auto_q_delta_pu::Float64 = 1e-4
  hysteresis_pu::Float64 = 0.01
  cooldown_iters::Int = 1
  guard::Bool = false
  guard_min_q_range_pu::Float64 = 1e-4
  guard_zero_range_mode::Symbol = :lock_pq
  guard_narrow_range_mode::Symbol = :prefer_pq
  guard_max_switches::Int = 10
  guard_freeze_after_repeated_switching::Bool = true
  guard_accept_bounded_violations::Bool = false
  guard_max_remaining_violations::Int = 0
  guard_violation_mode::Symbol = :delayed_switch
  guard_violation_threshold_pu::Float64 = 1e-4
  guard_log::Bool = true
  trace_buses::Vector{Int} = Int[]
  lock_pv_to_pq_buses::Vector{Int} = Int[]
  ignore_q_limits::Bool = false
  enforcement_mode::Symbol = :active_set
end

"""
    IslandPowerFlowConfig

Configuration for AC-island-aware power-flow diagnostics.
"""
Base.@kwdef struct IslandPowerFlowConfig
  enabled::Bool = true
  mode::Symbol = :solve_independent
  reference_policy::Symbol = :matpower_like
  diagnostic_continue_after_failure::Bool = true
end

"""
    PowerFlowConfig

Typed power-flow configuration. It owns solver tolerances, sparse execution
settings, automatic damping, start-mode controls, and Q-limit behavior.
"""
Base.@kwdef struct PowerFlowConfig
  method::Symbol = :rectangular
  tol::Float64 = 1.0e-8
  max_iter::Int = 30
  sparse::Bool = true
  autodamp::Bool = false
  autodamp_min::Float64 = 0.05
  wrong_branch_detection::Symbol = :warn
  wrong_branch_rescue::Bool = false
  wrong_branch_min_vm_pu::Float64 = 0.70
  wrong_branch_max_vm_pu::Float64 = 1.30
  wrong_branch_max_angle_spread_deg::Float64 = 180.0
  wrong_branch_max_branch_angle_deg::Float64 = 90.0
  wrong_branch_min_low_vm_count::Int = 1
  wrong_branch_rescue_max_attempts::Int = 2
  rectangular_workspace_reuse::Bool = true
  rectangular_preallocate_workspace::Symbol = :auto
  rectangular_workspace_min_buses::Int = 1000
  islands_enabled::Bool = true
  islands_mode::Symbol = :solve_independent
  islands_reference_policy::Symbol = :matpower_like
  start_mode::StartModeConfig = StartModeConfig()
  start_current_iteration::StartCurrentIterationConfig = StartCurrentIterationConfig()
  qlimits::QLimitConfig = QLimitConfig()
  islands::IslandPowerFlowConfig = IslandPowerFlowConfig()
end

const WRONG_BRANCH_DETECTION_VALUES = [:off, :warn, :fail, :rescue]
const POWERFLOW_ISLAND_MODE_VALUES = [:solve_independent]
const POWERFLOW_ISLAND_REFERENCE_POLICY_VALUES = [:matpower_like]

"""
    ObservabilityConfig

State-estimation observability diagnostic configuration.
"""
Base.@kwdef struct ObservabilityConfig
  enabled::Bool = true
end

"""
    StateEstimationConfig

Typed state-estimation configuration for future SE runners and diagnostics.
"""
Base.@kwdef struct StateEstimationConfig
  enabled::Bool = true
  method::Symbol = :wls
  tol::Float64 = 1.0e-8
  max_iter::Int = 20
  sparse::Bool = true
  flatstart::Bool = true
  jac_eps::Float64 = 1.0e-6
  update_net::Bool = true
  observability::ObservabilityConfig = ObservabilityConfig()
end

"""
    MatpowerImportConfig

Typed MATPOWER import/example configuration for case selection and import
conventions.
"""
Base.@kwdef struct MatpowerImportConfig
  case::String = ""
  cases::Vector{String} = String[]
  pv_voltage_source::Symbol = :gen_vg
  pv_voltage_mismatch_tol_pu::Float64 = 1.0e-4
  compare_voltage_reference::Symbol = :imported_setpoint
  shift_unit::Symbol = :deg
  shift_sign::Float64 = 1.0
  ratio::Symbol = :normal
  bus_shunt_model::Symbol = :admittance
  auto_profile::Symbol = :off
  auto_profile_log::Bool = true
  enable_pq_gen_controllers::Bool = true
  preallocate_network::Symbol = :auto
  preallocate_min_buses::Int = 1000
  apply_bus_names::Bool = false
  apply_branch_names::Bool = false
  apply_branch_kind::Bool = false
  import_for001_contingencies::Bool = true
  matpower_dcline_mode::Symbol = :pf_injections
end

"""
    TransformerConfig

Typed transformer-modeling configuration shared by all network importers.

`tap_changer_model` selects how transformer tap changers act on the equivalent
circuit:
- `:ideal` — tap steps only change the complex winding ratio; the series
  impedance keeps its neutral-position value (no impedance feedback; previous
  Sparlectra behavior).
- `:impedance_correction` — tap steps additionally re-refer the transformer
  series impedance through the tapped winding (R and X scaled with
  `|1 + f·e^(jφ)|²`, implemented centrally in `src/equicircuit.jl`).

The option applies to all transformers of an imported case and is read by both
the MATPOWER and the native DTF importer.
"""
Base.@kwdef struct TransformerConfig
  tap_changer_model::Symbol = :ideal
end

"""
    PerformanceConfig

Typed performance and diagnostic-volume configuration.
"""
Base.@kwdef struct PerformanceConfig
  enabled::Bool = false
  level::Symbol = :summary
  print_to_console::Bool = true
  write_to_logfile::Bool = true
  show_allocations::Bool = false
  show_iteration_table::Bool = false
  compact_logging::Bool = true
  representative_warmup_runs::Int = 0
  compare_cold_warm::Bool = false
  skip_reference_comparison::Bool = false
  skip_expensive_diagnostics::Bool = true
  skip_branch_neighborhood_report::Bool = true
  max_diagnostic_rows::Int = 25
end

Base.@kwdef struct BenchmarkConfig
  enabled::Bool = true
  methods::Vector{Symbol} = [:rectangular]
  seconds::Float64 = 2.0
  samples::Int = 50
  show_once::Bool = false
  show_once_output::Symbol = :classic
  show_once_max_nodes::Int = 0
end

"""
    RuntimeConfig

Typed runtime/threading configuration for example and benchmark entry points.
"""
Base.@kwdef struct RuntimeConfig
  julia_threads::String = "keep"
  blas_threads::String = "keep"
  print_thread_config::Bool = true
  casefile::String = ""
  case_name::String = ""
  case_source::String = ""
  configured_default_casefile::String = ""
end

"""
    DiagnosticsConfig

Typed diagnostic-output configuration shared by examples and future modules.
"""
Base.@kwdef struct DiagnosticsConfig
  log_effective_config::Bool = false
  console_summary::Bool = true
  console_auto_profile::Symbol = :compact
  console_diagnostics::Symbol = :compact
  console_q_limit_events::Symbol = :summary
  console_max_rows::Int = 100
  logfile_diagnostics::Symbol = :compact
end

"""
    OutputConfig

Typed output and logfile-format configuration. Console and logfile result
streams are intentionally independent so example runners can disable classic
logfile tables without suppressing compact console progress.
"""
Base.@kwdef struct OutputConfig
  console_summary::Bool = true
  console_auto_profile::Symbol = :compact
  console_diagnostics::Symbol = :compact
  console_q_limit_events::Symbol = :summary
  console_max_rows::Int = 100
  logfile_results::Symbol = :off
  result_table_max_rows::Int = 200
  result_table_large_case_threshold_buses::Int = 1000
  result_table_large_case_mode::Symbol = :summary
  detailed_result_csv_write_mode::Symbol = :auto
  detailed_result_csv_exporter::Symbol = :auto
  detailed_result_csv_direct_threshold_buses::Int = 10_000
  detailed_result_csv_buffer_initial_bytes::Int = 8 * 1024 * 1024
  detailed_result_csv_buffer_max_bytes::Int = 64 * 1024 * 1024
  detailed_result_csv_streaming_threshold_rows::Int = 100_000
  logfile_diagnostics::Symbol = :compact
  logfile_performance::Symbol = :compact
  logfile_warnings::Symbol = :table
end

"""
    SparlectraConfig

Central typed configuration assembled once at application or example boundaries.
Module-specific sections own their parsing and validation through constructors
such as `PowerFlowConfig(raw)` and `MatpowerImportConfig(raw)`.
"""
Base.@kwdef struct SparlectraConfig
  powerflow::PowerFlowConfig = PowerFlowConfig()
  state_estimation::StateEstimationConfig = StateEstimationConfig()
  matpower::MatpowerImportConfig = MatpowerImportConfig()
  transformer::TransformerConfig = TransformerConfig()
  performance::PerformanceConfig = PerformanceConfig()
  benchmark::BenchmarkConfig = BenchmarkConfig()
  runtime::RuntimeConfig = RuntimeConfig()
  diagnostics::DiagnosticsConfig = DiagnosticsConfig()
  output::OutputConfig = OutputConfig()
  control::ControlConfig = ControlConfig()
end

const _SPARLECTRA_CONFIG_CACHE = Ref{Union{Nothing,NamedTuple}}(nothing)
const ACTIVE_SPARLECTRA_CONFIG = Ref{SparlectraConfig}(SparlectraConfig())
const DEFAULT_SPARLECTRA_CONFIG_PATH = normpath(joinpath(@__DIR__, "configuration.yaml.example"))
const USER_SPARLECTRA_CONFIG_PATH = normpath(joinpath(@__DIR__, "..", "examples", "configuration.yaml"))

const SUPPORTED_POWERFLOW_METHOD = :rectangular
const POWERFLOW_START_ANGLE_MODE_VALUES = (:classic, :dc, :bus_va_blend, :matpower_va)
const POWERFLOW_START_VOLTAGE_MODE_VALUES = (:classic, :pv_gen_vg, :pv_bus_vm, :all_bus_vm, :profile_blend)
const POWERFLOW_START_PROFILE_SOURCE_VALUES = (:flat, :dc, :bus_metadata, :historical_profile, :matpower_reference, :state_estimation, :scada_snapshot)
const QLIMIT_START_MODE_VALUES = (:iteration, :auto, :iteration_or_auto)
const QLIMIT_ENFORCEMENT_MODE_VALUES = (:active_set, :classic_simultaneous, :classic_one_at_a_time)
const QLIMIT_ENFORCEMENT_MODE_LEGACY_ALIASES = Dict(
  :matpower_simultaneous => :classic_simultaneous,
  :matpower_one_at_a_time => :classic_one_at_a_time,
)
const QLIMIT_GUARD_ZERO_RANGE_MODE_VALUES = (:lock_pq,)
const QLIMIT_GUARD_NARROW_RANGE_MODE_VALUES = (:prefer_pq, :lock_pq)
const QLIMIT_GUARD_VIOLATION_MODE_VALUES = (:delayed_switch, :lock_pq)
const POWERFLOW_ISLAND_MODE_VALUES = (:solve_independent,)
const POWERFLOW_ISLAND_REFERENCE_POLICY_VALUES = (:matpower_like,)
const RECTANGULAR_PREALLOCATE_WORKSPACE_VALUES = (:off, :on, :auto)
const STATE_ESTIMATION_METHOD_VALUES = (:wls,)
const MATPOWER_PV_VOLTAGE_SOURCE_VALUES = (:gen_vg, :bus_vm, :auto, :strict_check)
const MATPOWER_COMPARE_VOLTAGE_REFERENCE_VALUES = (:bus_vm, :gen_vg, :imported_setpoint, :hybrid)
const MATPOWER_SHIFT_UNIT_VALUES = (:deg, :rad)
const MATPOWER_RATIO_VALUES = (:normal, :reciprocal)
const MATPOWER_BUS_SHUNT_MODEL_VALUES = (:admittance, :voltage_dependent_injection)
const MATPOWER_AUTO_PROFILE_VALUES = (:off, :recommend, :apply)
const MATPOWER_DCLINE_MODE_VALUES = (:reject_active, :ignore_inactive, :pf_injections)
const TRANSFORMER_TAP_CHANGER_MODEL_VALUES = (:ideal, :impedance_correction)
const PERFORMANCE_LEVEL_VALUES = (:off, :summary, :iteration, :full)
const OUTPUT_CONSOLE_AUTO_PROFILE_VALUES = (:off, :compact, :full)
const OUTPUT_CONSOLE_DIAGNOSTICS_VALUES = (:off, :compact, :summary, :full)
const OUTPUT_CONSOLE_Q_LIMIT_EVENTS_VALUES = (:off, :summary, :full)
const OUTPUT_LOGFILE_RESULTS_VALUES = (:off, :compact, :classic, :full)
const OUTPUT_RESULT_TABLE_LARGE_CASE_MODE_VALUES = (:summary, :classic, :full)
const OUTPUT_DETAILED_RESULT_CSV_WRITE_MODE_VALUES = (:auto, :buffered, :streaming)
const OUTPUT_DETAILED_RESULT_CSV_EXPORTER_VALUES = (:auto, :report, :direct)
const OUTPUT_LOGFILE_DIAGNOSTICS_VALUES = (:off, :compact, :full)
const OUTPUT_LOGFILE_PERFORMANCE_VALUES = (:off, :compact, :full)
const OUTPUT_LOGFILE_WARNINGS_VALUES = (:off, :summary, :table, :full)
const BENCHMARK_SHOW_ONCE_OUTPUT_VALUES = (:classic, :dataframe, :compact)

_powerflow_method_label(method)::Symbol = Symbol(method) === :polar ? :polar_full : Symbol(method)

function unsupported_powerflow_method_message(method)::String
  return "Only the rectangular AC power-flow solver is supported in this version.\nRequested method: $(_powerflow_method_label(method))."
end

_as_symbol_vector_cfg(x)::Vector{Symbol} = x isa AbstractVector ? Symbol[_as_symbol_cfg(item) for item in x] : [_as_symbol_cfg(x)]

function _validate_rectangular_powerflow_options(; method::Symbol = SUPPORTED_POWERFLOW_METHOD, sparse::Bool = true)
  requested = _powerflow_method_label(method)
  requested === SUPPORTED_POWERFLOW_METHOD || throw(ArgumentError(unsupported_powerflow_method_message(requested)))
  sparse || throw(ArgumentError("Only sparse power-flow matrices are supported in this version. Requested sparse=false."))
  return nothing
end

function _validate_powerflow_solver_config!(raw::AbstractDict)
  pf = _merged_section(raw, "powerflow")
  method = _as_symbol_cfg(_raw_get(pf, "method", SUPPORTED_POWERFLOW_METHOD))
  if haskey(pf, "sparse") || haskey(pf, :sparse) || haskey(pf, "opt_sparse") || haskey(pf, :opt_sparse)
    throw(ArgumentError("power_flow.sparse/opt_sparse is obsolete and no longer configurable. Sparse matrices are mandatory."))
  end
  _validate_rectangular_powerflow_options(method = method, sparse = true)
  benchmark = _raw_section(raw, "benchmark")
  if haskey(benchmark, "methods")
    for m in _as_symbol_vector_cfg(benchmark["methods"])
      requested = _powerflow_method_label(m)
      requested === SUPPORTED_POWERFLOW_METHOD || throw(ArgumentError(unsupported_powerflow_method_message(requested)))
    end
  end
  return nothing
end

function set_sparlectra_config!(cfg::SparlectraConfig)
  ACTIVE_SPARLECTRA_CONFIG[] = cfg
  return cfg
end

active_sparlectra_config()::SparlectraConfig = ACTIVE_SPARLECTRA_CONFIG[]
powerflow_config()::PowerFlowConfig = active_sparlectra_config().powerflow
matpower_import_config()::MatpowerImportConfig = active_sparlectra_config().matpower
transformer_config()::TransformerConfig = active_sparlectra_config().transformer
state_estimation_config()::StateEstimationConfig = active_sparlectra_config().state_estimation
diagnostics_config()::DiagnosticsConfig = active_sparlectra_config().diagnostics
output_config()::OutputConfig = active_sparlectra_config().output
performance_config()::PerformanceConfig = active_sparlectra_config().performance
benchmark_config()::BenchmarkConfig = active_sparlectra_config().benchmark
runtime_config()::RuntimeConfig = active_sparlectra_config().runtime
control_config()::ControlConfig = active_sparlectra_config().control

function configuration_path_from_inputs(; env_var::AbstractString = "SPARLECTRA_CONFIGURATION_YAML", fallback_paths::AbstractVector{<:AbstractString} = String[])
  candidate = strip(get(ENV, env_var, ""))
  !isempty(candidate) && return candidate
  for path in fallback_paths
    file = strip(String(path))
    isempty(file) && continue
    isfile(file) && return file
  end
  return ""
end

_config_key(key::Symbol) = String(key)
_config_key(key) = String(key)

function _raw_get(raw::AbstractDict, key::AbstractString, default)
  haskey(raw, key) && return raw[key]
  symkey = Symbol(key)
  haskey(raw, symkey) && return raw[symkey]
  return default
end

function _raw_section(raw::AbstractDict, key::AbstractString, aliases::AbstractVector{<:AbstractString} = String[])
  value = _raw_get(raw, key, nothing)
  if value === nothing
    for alias in aliases
      value = _raw_get(raw, alias, nothing)
      value === nothing || break
    end
  end
  value === nothing && (value = Dict{String,Any}())
  value isa AbstractDict && return value
  throw(ArgumentError("Configuration section $(repr(key)) must be a dictionary."))
end

_as_string_cfg(x)::String = x isa AbstractString ? String(x) : string(x)
function _as_string_vector_cfg(x)::Vector{String}
  x isa AbstractVector || throw(ArgumentError("matpower_import.cases must be a vector of non-empty case names."))
  values = String[strip(_as_string_cfg(item)) for item in x]
  any(isempty, values) && throw(ArgumentError("matpower_import.cases must not contain empty case names."))
  return values
end
_as_symbol_cfg(x)::Symbol = x isa Symbol ? x : (x isa Bool ? (x ? :on : :off) : Symbol(lowercase(String(x))))
function _as_auto_profile_symbol_cfg(x)::Symbol
  x === false && return :off
  x === true && return :apply
  sym = _as_symbol_cfg(x)
  sym in (:false, :no, :none) && return :off
  sym in (:true, :yes) && return :apply
  return sym
end
_as_bool_cfg(x)::Bool = as_bool(x)
_as_int_cfg(x)::Int = x isa Integer ? Int(x) : parse(Int, String(x))
_as_float_cfg(x)::Float64 = x isa Real ? Float64(x) : parse(Float64, String(x))

function _as_float_vector_cfg(x)::Vector{Float64}
  x isa AbstractVector || throw(ArgumentError("Expected a vector of floating-point values; got $(repr(x))."))
  return Float64[_as_float_cfg(item) for item in x]
end

function _as_int_vector_cfg(x)::Vector{Int}
  return as_int_vector(x)
end

function _validate_nonnegative(name::AbstractString, value::Real)
  value >= 0 || throw(ArgumentError("$(name) must be non-negative; got $(value)."))
  return value
end
function _validate_allowed_symbol(name::AbstractString, value::Symbol, allowed::Tuple)
  if !(value in allowed)
    guidance = name == "power_flow.start_mode.voltage_mode" && value === :bus_vm_va_blend ?
      " The former bus_vm_va_blend alias has been removed; use voltage_mode: profile_blend with profile_source: matpower_reference." : ""
    allowed_text = join(string.(allowed), ", ")
    throw(ArgumentError("Invalid value:\n$(name) = \"$(value)\"\n\nAllowed values:\n$(allowed_text)\n\n$(name) must be one of $(collect(allowed)); got $(value).$(guidance)"))
  end
  return value
end

function _validate_allowed_symbol(name::AbstractString, value::Symbol, allowed::AbstractVector{Symbol})
  if !(value in allowed)
    allowed_text = join(string.(allowed), ", ")
    throw(ArgumentError("Invalid value:\n$(name) = \"$(value)\"\n\nAllowed values:\n$(allowed_text)\n\n$(name) must be one of $(allowed); got $(value)."))
  end
  return value
end

function _canonical_qlimit_enforcement_mode(value::Symbol)::Symbol
  canonical = get(QLIMIT_ENFORCEMENT_MODE_LEGACY_ALIASES, value, value)
  if canonical in QLIMIT_ENFORCEMENT_MODE_VALUES
    return canonical
  end
  legacy_aliases = ["$(alias) -> $(target)" for (alias, target) in sort(collect(QLIMIT_ENFORCEMENT_MODE_LEGACY_ALIASES); by = first)]
  throw(ArgumentError("power_flow.qlimits.enforcement_mode must be one of $(collect(QLIMIT_ENFORCEMENT_MODE_VALUES)); got $(value). Legacy aliases accepted: $(legacy_aliases)."))
end

function _validate_positive(name::AbstractString, value::Real)
  value > 0 || throw(ArgumentError("$(name) must be positive; got $(value)."))
  return value
end

function StartModeConfig(raw::AbstractDict)
  angle_mode = _validate_allowed_symbol("power_flow.start_mode.angle_mode", _as_symbol_cfg(_raw_get(raw, "angle_mode", :dc)), POWERFLOW_START_ANGLE_MODE_VALUES)
  voltage_mode = _validate_allowed_symbol("power_flow.start_mode.voltage_mode", _as_symbol_cfg(_raw_get(raw, "voltage_mode", :profile_blend)), POWERFLOW_START_VOLTAGE_MODE_VALUES)
  profile_source = _validate_allowed_symbol("power_flow.start_mode.profile_source", _as_symbol_cfg(_raw_get(raw, "profile_source", :matpower_reference)), POWERFLOW_START_PROFILE_SOURCE_VALUES)
  return StartModeConfig(
    flatstart = _as_bool_cfg(_raw_get(raw, "flatstart", _raw_get(raw, "opt_flatstart", false))),
    angle_mode = angle_mode,
    voltage_mode = voltage_mode,
    profile_source = profile_source,
    start_projection = _as_bool_cfg(_raw_get(raw, "start_projection", false)),
    try_dc_start = _as_bool_cfg(_raw_get(raw, "try_dc_start", _raw_get(raw, "start_projection_try_dc_start", true))),
    try_blend_scan = _as_bool_cfg(_raw_get(raw, "try_blend_scan", _raw_get(raw, "start_projection_try_blend_scan", true))),
    branch_guard = _as_bool_cfg(_raw_get(raw, "branch_guard", _raw_get(raw, "start_projection_branch_guard", true))),
    measure_candidates = _as_bool_cfg(_raw_get(raw, "measure_candidates", _raw_get(raw, "start_projection_measure_candidates", true))),
    accept_unmeasured_dc_start = _as_bool_cfg(_raw_get(raw, "accept_unmeasured_dc_start", _raw_get(raw, "start_projection_accept_unmeasured_dc_start", false))),
    reuse_import_data = _as_bool_cfg(_raw_get(raw, "reuse_import_data", _raw_get(raw, "start_projection_reuse_import_data", true))),
    blend_lambdas = _as_float_vector_cfg(_raw_get(raw, "blend_lambdas", _raw_get(raw, "start_projection_blend_lambdas", [0.25, 0.5, 0.75]))),
    dc_angle_limit_deg = _validate_positive("start_projection_dc_angle_limit_deg", _as_float_cfg(_raw_get(raw, "dc_angle_limit_deg", _raw_get(raw, "start_projection_dc_angle_limit_deg", 60.0)))),
  )
end

function StartCurrentIterationConfig(raw::AbstractDict)
  max_iter = _as_int_cfg(_raw_get(raw, "max_iter", 10))
  max_iter >= 0 || throw(ArgumentError("power_flow.start_current_iteration.max_iter must be non-negative; got $(max_iter)."))
  damping = _as_float_cfg(_raw_get(raw, "damping", 0.5))
  isfinite(damping) && 0.0 < damping <= 1.0 || throw(ArgumentError("power_flow.start_current_iteration.damping must be finite and in (0, 1]; got $(damping)."))
  min_improvement_factor = _as_float_cfg(_raw_get(raw, "min_improvement_factor", 0.98))
  isfinite(min_improvement_factor) && min_improvement_factor >= 0.0 || throw(ArgumentError("power_flow.start_current_iteration.min_improvement_factor must be finite and non-negative; got $(min_improvement_factor)."))
  vm_min_pu = _validate_positive("power_flow.start_current_iteration.vm_min_pu", _as_float_cfg(_raw_get(raw, "vm_min_pu", 0.5)))
  vm_max_pu = _validate_positive("power_flow.start_current_iteration.vm_max_pu", _as_float_cfg(_raw_get(raw, "vm_max_pu", 1.5)))
  vm_min_pu <= vm_max_pu || throw(ArgumentError("power_flow.start_current_iteration.vm_min_pu must be <= vm_max_pu."))
  return StartCurrentIterationConfig(
    enabled = _as_bool_cfg(_raw_get(raw, "enabled", false)),
    max_iter = max_iter,
    tol = _validate_positive("power_flow.start_current_iteration.tol", _as_float_cfg(_raw_get(raw, "tol", 1.0e-3))),
    damping = damping,
    accept_only_if_improved = _as_bool_cfg(_raw_get(raw, "accept_only_if_improved", true)),
    min_improvement_factor = min_improvement_factor,
    vm_min_pu = vm_min_pu,
    vm_max_pu = vm_max_pu,
    max_angle_step_deg = _validate_positive("power_flow.start_current_iteration.max_angle_step_deg", _as_float_cfg(_raw_get(raw, "max_angle_step_deg", 30.0))),
    only_for_large_cases = _as_bool_cfg(_raw_get(raw, "only_for_large_cases", false)),
  )
end

function QLimitConfig(raw::AbstractDict)
  qlimits_enabled = _raw_get(raw, "enabled", true)
  guard_raw = _raw_get(raw, "guard", Dict{String,Any}())
  guard_enabled_default = guard_raw isa AbstractDict ? _as_bool_cfg(_raw_get(guard_raw, "enabled", false)) : _as_bool_cfg(guard_raw)
  guard_cfg = guard_raw isa AbstractDict ? guard_raw : Dict{String,Any}()
  merged = merge(Dict{Any,Any}(raw), Dict{Any,Any}(guard_cfg))
  return QLimitConfig(
    start_iter = _as_int_cfg(_raw_get(merged, "start_iter", _raw_get(merged, "qlimit_start_iter", 2))),
    start_mode = _validate_allowed_symbol(
      "power_flow.qlimits.start_mode",
      _as_symbol_cfg((_raw_get(merged, "start_mode", nothing) isa AbstractDict) ? _raw_get(merged, "qlimit_start_mode", :iteration) : _raw_get(merged, "start_mode", _raw_get(merged, "qlimit_start_mode", :iteration))),
      QLIMIT_START_MODE_VALUES,
    ),
    auto_q_delta_pu = _validate_nonnegative("qlimit_auto_q_delta_pu", _as_float_cfg(_raw_get(merged, "auto_q_delta_pu", _raw_get(merged, "qlimit_auto_q_delta_pu", 1e-4)))),
    hysteresis_pu = _validate_nonnegative("power_flow.qlimits.hysteresis_pu", _as_float_cfg(_raw_get(merged, "hysteresis_pu", _raw_get(merged, "q_hyst_pu", 0.01)))),
    cooldown_iters = _as_int_cfg(_raw_get(merged, "cooldown_iters", 1)),
    guard = _as_bool_cfg(_raw_get(raw, "qlimit_guard", guard_enabled_default)),
    guard_min_q_range_pu = _validate_nonnegative("qlimit_guard_min_q_range_pu", _as_float_cfg(_raw_get(merged, "min_q_range_pu", _raw_get(merged, "guard_min_q_range_pu", _raw_get(merged, "qlimit_guard_min_q_range_pu", 1e-4))))),
    guard_zero_range_mode = _validate_allowed_symbol("power_flow.qlimits.guard.zero_range_mode", _as_symbol_cfg(_raw_get(merged, "zero_range_mode", _raw_get(merged, "guard_zero_range_mode", _raw_get(merged, "qlimit_guard_zero_range_mode", :lock_pq)))), QLIMIT_GUARD_ZERO_RANGE_MODE_VALUES),
    guard_narrow_range_mode = _validate_allowed_symbol(
      "power_flow.qlimits.guard.narrow_range_mode",
      _as_symbol_cfg(_raw_get(merged, "narrow_range_mode", _raw_get(merged, "guard_narrow_range_mode", _raw_get(merged, "qlimit_guard_narrow_range_mode", :prefer_pq)))),
      QLIMIT_GUARD_NARROW_RANGE_MODE_VALUES,
    ),
    guard_max_switches = _as_int_cfg(_raw_get(merged, "max_switches", _raw_get(merged, "guard_max_switches", _raw_get(merged, "qlimit_guard_max_switches", 10)))),
    guard_freeze_after_repeated_switching = _as_bool_cfg(_raw_get(merged, "freeze_after_repeated_switching", _raw_get(merged, "guard_freeze_after_repeated_switching", _raw_get(merged, "qlimit_guard_freeze_after_repeated_switching", true)))),
    guard_accept_bounded_violations = _as_bool_cfg(_raw_get(merged, "accept_bounded_violations", _raw_get(merged, "guard_accept_bounded_violations", _raw_get(merged, "qlimit_guard_accept_bounded_violations", false)))),
    guard_max_remaining_violations = _as_int_cfg(_raw_get(merged, "max_remaining_violations", _raw_get(merged, "guard_max_remaining_violations", _raw_get(merged, "qlimit_guard_max_remaining_violations", 0)))),
    guard_violation_mode = _validate_allowed_symbol("power_flow.qlimits.guard.violation_mode", _as_symbol_cfg(_raw_get(merged, "violation_mode", _raw_get(merged, "guard_violation_mode", _raw_get(merged, "qlimit_guard_violation_mode", :delayed_switch)))), QLIMIT_GUARD_VIOLATION_MODE_VALUES),
    guard_violation_threshold_pu = _validate_nonnegative("qlimit_guard_violation_threshold_pu", _as_float_cfg(_raw_get(merged, "violation_threshold_pu", _raw_get(merged, "guard_violation_threshold_pu", _raw_get(merged, "qlimit_guard_violation_threshold_pu", 1e-4))))),
    guard_log = _as_bool_cfg(_raw_get(merged, "log", _raw_get(merged, "guard_log", _raw_get(merged, "qlimit_guard_log", true)))),
    trace_buses = _as_int_vector_cfg(_raw_get(merged, "trace_buses", _raw_get(merged, "qlimit_trace_buses", Int[]))),
    lock_pv_to_pq_buses = _as_int_vector_cfg(_raw_get(merged, "lock_pv_to_pq_buses", Int[])),
    ignore_q_limits = _as_bool_cfg(_raw_get(raw, "ignore_q_limits", qlimits_enabled == false)),
    enforcement_mode = _canonical_qlimit_enforcement_mode(_as_symbol_cfg(_raw_get(merged, "enforcement_mode", :active_set))),
  )
end

function _merged_section(raw::AbstractDict, section_name::AbstractString)
  aliases = section_name == "powerflow" ? ["power_flow"] : section_name == "matpower" ? ["matpower_import"] : String[]
  canonical = section_name == "power_flow" ? "powerflow" : section_name == "matpower_import" ? "matpower" : section_name
  primary = Dict{Any,Any}(_raw_section(raw, canonical))
  for alias in aliases
    merge!(primary, Dict{Any,Any}(_raw_section(raw, alias)))
  end
  return primary
end

function PowerFlowConfig(raw::AbstractDict)
  merged = _merged_section(raw, "powerflow")
  method = _as_symbol_cfg(_raw_get(merged, "method", SUPPORTED_POWERFLOW_METHOD))
  if haskey(merged, "sparse") || haskey(merged, :sparse) || haskey(merged, "opt_sparse") || haskey(merged, :opt_sparse)
    throw(ArgumentError("power_flow.sparse/opt_sparse is obsolete and no longer configurable. Sparse matrices are mandatory."))
  end
  _validate_rectangular_powerflow_options(method = method, sparse = true)
  start_raw = _raw_get(merged, "start_values", _raw_section(merged, "start_mode"))
  start_current_iteration_raw = _raw_section(merged, "start_current_iteration")
  qlimit_raw = _raw_section(merged, "qlimits")
  islands_raw = _raw_section(merged, "islands")
  wrong_branch_min_vm_pu = _validate_nonnegative("power_flow.wrong_branch_min_vm_pu", _as_float_cfg(_raw_get(merged, "wrong_branch_min_vm_pu", 0.70)))
  wrong_branch_max_vm_pu = _validate_positive("power_flow.wrong_branch_max_vm_pu", _as_float_cfg(_raw_get(merged, "wrong_branch_max_vm_pu", 1.30)))
  wrong_branch_min_vm_pu <= wrong_branch_max_vm_pu || throw(ArgumentError("power_flow.wrong_branch_min_vm_pu must be <= power_flow.wrong_branch_max_vm_pu."))
  wrong_branch_min_low_vm_count = _as_int_cfg(_raw_get(merged, "wrong_branch_min_low_vm_count", 1))
  wrong_branch_min_low_vm_count >= 0 || throw(ArgumentError("power_flow.wrong_branch_min_low_vm_count must be >= 0."))
  wrong_branch_rescue_max_attempts = _as_int_cfg(_raw_get(merged, "wrong_branch_rescue_max_attempts", 2))
  wrong_branch_rescue_max_attempts >= 0 || throw(ArgumentError("power_flow.wrong_branch_rescue_max_attempts must be >= 0."))
  return PowerFlowConfig(
    method = method,
    tol = _validate_positive("powerflow.tol", _as_float_cfg(_raw_get(merged, "tol", 1.0e-8))),
    max_iter = _as_int_cfg(_raw_get(merged, "max_iter", _raw_get(merged, "max_ite", 30))),
    sparse = true,
    autodamp = _as_bool_cfg(_raw_get(merged, "autodamp", false)),
    autodamp_min = _validate_positive("powerflow.autodamp_min", _as_float_cfg(_raw_get(merged, "autodamp_min", 0.05))),
    wrong_branch_detection = _validate_allowed_symbol("power_flow.wrong_branch_detection", _as_symbol_cfg(_raw_get(merged, "wrong_branch_detection", :warn)), WRONG_BRANCH_DETECTION_VALUES),
    wrong_branch_rescue = _as_bool_cfg(_raw_get(merged, "wrong_branch_rescue", false)),
    wrong_branch_min_vm_pu = wrong_branch_min_vm_pu,
    wrong_branch_max_vm_pu = wrong_branch_max_vm_pu,
    wrong_branch_max_angle_spread_deg = _validate_nonnegative("power_flow.wrong_branch_max_angle_spread_deg", _as_float_cfg(_raw_get(merged, "wrong_branch_max_angle_spread_deg", 180.0))),
    wrong_branch_max_branch_angle_deg = _validate_nonnegative("power_flow.wrong_branch_max_branch_angle_deg", _as_float_cfg(_raw_get(merged, "wrong_branch_max_branch_angle_deg", 90.0))),
    wrong_branch_min_low_vm_count = wrong_branch_min_low_vm_count,
    wrong_branch_rescue_max_attempts = wrong_branch_rescue_max_attempts,
    rectangular_workspace_reuse = _as_bool_cfg(_raw_get(merged, "rectangular_workspace_reuse", true)),
    rectangular_preallocate_workspace = _validate_allowed_symbol("power_flow.rectangular_preallocate_workspace", _as_symbol_cfg(_raw_get(merged, "rectangular_preallocate_workspace", :auto)), RECTANGULAR_PREALLOCATE_WORKSPACE_VALUES),
    rectangular_workspace_min_buses = _as_int_cfg(_raw_get(merged, "rectangular_workspace_min_buses", 1000)),
    islands_enabled = _as_bool_cfg(_raw_get(islands_raw, "enabled", true)),
    islands_mode = _validate_allowed_symbol("power_flow.islands.mode", _as_symbol_cfg(_raw_get(islands_raw, "mode", :solve_independent)), POWERFLOW_ISLAND_MODE_VALUES),
    islands_reference_policy = _validate_allowed_symbol("power_flow.islands.reference_policy", _as_symbol_cfg(_raw_get(islands_raw, "reference_policy", :matpower_like)), POWERFLOW_ISLAND_REFERENCE_POLICY_VALUES),
    start_mode = StartModeConfig(merge(Dict{Any,Any}(merged), Dict{Any,Any}(start_raw))),
    start_current_iteration = StartCurrentIterationConfig(start_current_iteration_raw),
    qlimits = QLimitConfig(merge(Dict{Any,Any}(merged), Dict{Any,Any}(qlimit_raw))),
    islands = IslandPowerFlowConfig(islands_raw),
  )
end

function IslandPowerFlowConfig(raw::AbstractDict)
  merged = haskey(raw, "islands") || haskey(raw, :islands) ? _merged_section(raw, "islands") : raw
  return IslandPowerFlowConfig(
    enabled = _as_bool_cfg(_raw_get(merged, "enabled", true)),
    mode = _validate_allowed_symbol("power_flow.islands.mode", _as_symbol_cfg(_raw_get(merged, "mode", :solve_independent)), POWERFLOW_ISLAND_MODE_VALUES),
    reference_policy = _validate_allowed_symbol("power_flow.islands.reference_policy", _as_symbol_cfg(_raw_get(merged, "reference_policy", :matpower_like)), POWERFLOW_ISLAND_REFERENCE_POLICY_VALUES),
    diagnostic_continue_after_failure = _as_bool_cfg(_raw_get(merged, "diagnostic_continue_after_failure", true)),
  )
end

function ObservabilityConfig(raw::AbstractDict)
  merged = _merged_section(raw, "observability")
  return ObservabilityConfig(enabled = _as_bool_cfg(_raw_get(merged, "enabled", true)))
end

function StateEstimationConfig(raw::AbstractDict)
  merged = _merged_section(raw, "state_estimation")
  if haskey(merged, "sparse") || haskey(merged, :sparse)
    throw(ArgumentError("state_estimation.sparse is obsolete and no longer configurable. Sparse matrices are mandatory."))
  end
  return StateEstimationConfig(
    enabled = _as_bool_cfg(_raw_get(merged, "enabled", true)),
    method = _validate_allowed_symbol("state_estimation.method", _as_symbol_cfg(_raw_get(merged, "method", :wls)), STATE_ESTIMATION_METHOD_VALUES),
    tol = _validate_positive("state_estimation.tol", _as_float_cfg(_raw_get(merged, "tol", 1.0e-8))),
    max_iter = _as_int_cfg(_raw_get(merged, "max_iter", 20)),
    sparse = true,
    flatstart = _as_bool_cfg(_raw_get(merged, "flatstart", true)),
    jac_eps = _validate_positive("state_estimation.jac_eps", _as_float_cfg(_raw_get(merged, "jac_eps", 1.0e-6))),
    update_net = _as_bool_cfg(_raw_get(merged, "update_net", true)),
    observability = ObservabilityConfig(merged),
  )
end

function MatpowerImportConfig(raw::AbstractDict)
  merged = _merged_section(raw, "matpower")
  return MatpowerImportConfig(
    case = strip(_as_string_cfg(_raw_get(merged, "case", _raw_get(raw, "case", "")))),
    cases = _as_string_vector_cfg(_raw_get(merged, "cases", String[])),
    pv_voltage_source = _validate_allowed_symbol("matpower_import.pv_voltage_source", _as_symbol_cfg(_raw_get(merged, "pv_voltage_source", _raw_get(merged, "matpower_pv_voltage_source", :gen_vg))), MATPOWER_PV_VOLTAGE_SOURCE_VALUES),
    pv_voltage_mismatch_tol_pu = _validate_nonnegative("matpower.pv_voltage_mismatch_tol_pu", _as_float_cfg(_raw_get(merged, "pv_voltage_mismatch_tol_pu", _raw_get(merged, "matpower_pv_voltage_mismatch_tol_pu", 1.0e-4)))),
    compare_voltage_reference = _validate_allowed_symbol("matpower_import.compare_voltage_reference", _as_symbol_cfg(_raw_get(merged, "compare_voltage_reference", :imported_setpoint)), MATPOWER_COMPARE_VOLTAGE_REFERENCE_VALUES),
    shift_unit = _validate_allowed_symbol("matpower_import.shift_unit", _as_symbol_cfg(_raw_get(merged, "shift_unit", _raw_get(merged, "matpower_shift_unit", :deg))), MATPOWER_SHIFT_UNIT_VALUES),
    shift_sign = _as_float_cfg(_raw_get(merged, "shift_sign", _raw_get(merged, "matpower_shift_sign", 1.0))),
    ratio = _validate_allowed_symbol("matpower_import.ratio", _as_symbol_cfg(_raw_get(merged, "ratio", _raw_get(merged, "matpower_ratio", :normal))), MATPOWER_RATIO_VALUES),
    bus_shunt_model = _validate_allowed_symbol("matpower_import.bus_shunt_model", _as_symbol_cfg(_raw_get(merged, "bus_shunt_model", :admittance)), MATPOWER_BUS_SHUNT_MODEL_VALUES),
    auto_profile = _validate_allowed_symbol("matpower_import.auto_profile", _as_auto_profile_symbol_cfg(_raw_get(merged, "auto_profile", :recommend)), MATPOWER_AUTO_PROFILE_VALUES),
    auto_profile_log = _as_bool_cfg(_raw_get(merged, "auto_profile_log", true)),
    enable_pq_gen_controllers = _as_bool_cfg(_raw_get(merged, "enable_pq_gen_controllers", true)),
    preallocate_network = _validate_allowed_symbol("matpower_import.preallocate_network", _as_symbol_cfg(_raw_get(merged, "preallocate_network", :auto)), [:off, :on, :auto]),
    preallocate_min_buses = _as_int_cfg(_raw_get(merged, "preallocate_min_buses", 1000)),
    apply_bus_names = _as_bool_cfg(_raw_get(merged, "apply_bus_names", false)),
    apply_branch_names = _as_bool_cfg(_raw_get(merged, "apply_branch_names", false)),
    apply_branch_kind = _as_bool_cfg(_raw_get(merged, "apply_branch_kind", false)),
    import_for001_contingencies = _as_bool_cfg(_raw_get(merged, "import_for001_contingencies", true)),
    matpower_dcline_mode = _validate_allowed_symbol("matpower_import.matpower_dcline_mode", _as_symbol_cfg(_raw_get(merged, "matpower_dcline_mode", :pf_injections)), MATPOWER_DCLINE_MODE_VALUES),
  )
end

function TransformerConfig(raw::AbstractDict)
  merged = _merged_section(raw, "transformer")
  return TransformerConfig(
    tap_changer_model = _validate_allowed_symbol("transformer.tap_changer_model", _as_symbol_cfg(_raw_get(merged, "tap_changer_model", :ideal)), TRANSFORMER_TAP_CHANGER_MODEL_VALUES),
  )
end

"""
    configured_matpower_cases(config) -> Vector{String}

Return configured MATPOWER cases in deterministic execution order. A non-empty
`matpower.cases` list takes precedence over the compatible single-case
`matpower.case` setting.
"""
function configured_matpower_cases(mat_cfg::MatpowerImportConfig)::Vector{String}
  !isempty(mat_cfg.cases) && return copy(mat_cfg.cases)
  case_name = strip(mat_cfg.case)
  return isempty(case_name) ? String[] : [case_name]
end

configured_matpower_cases(cfg::SparlectraConfig)::Vector{String} = configured_matpower_cases(cfg.matpower)

function PerformanceConfig(raw::AbstractDict)
  merged = _merged_section(raw, "performance")
  return PerformanceConfig(
    enabled = _as_bool_cfg(_raw_get(merged, "enabled", _raw_get(merged, "performance_enabled", false))),
    level = _validate_allowed_symbol("performance.level", _as_symbol_cfg(_raw_get(merged, "level", _raw_get(merged, "performance_level", :summary))), PERFORMANCE_LEVEL_VALUES),
    print_to_console = _as_bool_cfg(_raw_get(merged, "print_to_console", _raw_get(merged, "performance_print_to_console", true))),
    write_to_logfile = _as_bool_cfg(_raw_get(merged, "write_to_logfile", _raw_get(merged, "performance_write_to_logfile", true))),
    show_allocations = _as_bool_cfg(_raw_get(merged, "show_allocations", _raw_get(merged, "performance_show_allocations", false))),
    show_iteration_table = _as_bool_cfg(_raw_get(merged, "show_iteration_table", _raw_get(merged, "performance_show_iteration_table", false))),
    compact_logging = _as_bool_cfg(_raw_get(merged, "compact_logging", _raw_get(merged, "performance_compact_logging", true))),
    representative_warmup_runs = _as_int_cfg(_raw_get(merged, "representative_warmup_runs", 0)),
    compare_cold_warm = _as_bool_cfg(_raw_get(merged, "compare_cold_warm", false)),
    skip_reference_comparison = _as_bool_cfg(_raw_get(merged, "skip_reference_comparison", _raw_get(merged, "performance_skip_reference_comparison", false))),
    skip_expensive_diagnostics = _as_bool_cfg(_raw_get(merged, "skip_expensive_diagnostics", _raw_get(merged, "performance_skip_expensive_diagnostics", true))),
    skip_branch_neighborhood_report = _as_bool_cfg(_raw_get(merged, "skip_branch_neighborhood_report", _raw_get(merged, "performance_skip_branch_neighborhood_report", true))),
    max_diagnostic_rows = _as_int_cfg(_raw_get(merged, "max_diagnostic_rows", _raw_get(merged, "performance_max_diagnostic_rows", 25))),
  )
end

function RuntimeConfig(raw::AbstractDict)
  merged = _merged_section(raw, "runtime")
  return RuntimeConfig(
    julia_threads = _as_string_cfg(_raw_get(merged, "julia_threads", "keep")),
    blas_threads = _as_string_cfg(_raw_get(merged, "blas_threads", "keep")),
    print_thread_config = _as_bool_cfg(_raw_get(merged, "print_thread_config", true)),
    casefile = _as_string_cfg(_raw_get(merged, "casefile", "")),
    case_name = _as_string_cfg(_raw_get(merged, "case_name", "")),
    case_source = _as_string_cfg(_raw_get(merged, "case_source", "")),
    configured_default_casefile = _as_string_cfg(_raw_get(merged, "configured_default_casefile", "")),
  )
end

function BenchmarkConfig(raw::AbstractDict)
  merged = _merged_section(raw, "benchmark")
  methods = _as_symbol_vector_cfg(_raw_get(merged, "methods", [:rectangular]))
  for method in methods
    _validate_rectangular_powerflow_options(method = method, sparse = true)
  end
  return BenchmarkConfig(
    enabled = _as_bool_cfg(_raw_get(merged, "enabled", true)),
    methods = methods,
    seconds = _validate_positive("benchmark.seconds", _as_float_cfg(_raw_get(merged, "seconds", 2.0))),
    samples = _as_int_cfg(_raw_get(merged, "samples", 50)),
    show_once = _as_bool_cfg(_raw_get(merged, "show_once", false)),
    show_once_output = _validate_allowed_symbol("benchmark.show_once_output", _as_symbol_cfg(_raw_get(merged, "show_once_output", :classic)), BENCHMARK_SHOW_ONCE_OUTPUT_VALUES),
    show_once_max_nodes = _as_int_cfg(_raw_get(merged, "show_once_max_nodes", 0)),
  )
end

function DiagnosticsConfig(raw::AbstractDict)
  merged = _merged_section(raw, "diagnostics")
  return DiagnosticsConfig(
    log_effective_config = _as_bool_cfg(_raw_get(merged, "log_effective_config", false)),
    console_summary = _as_bool_cfg(_raw_get(merged, "console_summary", true)),
    console_auto_profile = _validate_allowed_symbol("diagnostics.console_auto_profile", _as_symbol_cfg(_raw_get(merged, "console_auto_profile", :compact)), OUTPUT_CONSOLE_AUTO_PROFILE_VALUES),
    console_diagnostics = _validate_allowed_symbol("diagnostics.console_diagnostics", _as_symbol_cfg(_raw_get(merged, "console_diagnostics", :compact)), OUTPUT_CONSOLE_DIAGNOSTICS_VALUES),
    console_q_limit_events = _validate_allowed_symbol("diagnostics.console_q_limit_events", _as_symbol_cfg(_raw_get(merged, "console_q_limit_events", :summary)), OUTPUT_CONSOLE_Q_LIMIT_EVENTS_VALUES),
    console_max_rows = _as_int_cfg(_raw_get(merged, "console_max_rows", 100)),
    logfile_diagnostics = _validate_allowed_symbol("diagnostics.logfile_diagnostics", _as_symbol_cfg(_raw_get(merged, "logfile_diagnostics", :compact)), OUTPUT_LOGFILE_DIAGNOSTICS_VALUES),
  )
end

_output_nonnegative_or_default(value::Integer, default::Integer) = value < 0 ? default : value
_output_positive_or_default(value::Integer, default::Integer) = value <= 0 ? default : value

function OutputConfig(raw::AbstractDict)
  merged = _merged_section(raw, "output")
  return OutputConfig(
    console_summary = _as_bool_cfg(_raw_get(merged, "console_summary", true)),
    console_auto_profile = _validate_allowed_symbol("output.console_auto_profile", _as_symbol_cfg(_raw_get(merged, "console_auto_profile", :compact)), OUTPUT_CONSOLE_AUTO_PROFILE_VALUES),
    console_diagnostics = _validate_allowed_symbol("output.console_diagnostics", _as_symbol_cfg(_raw_get(merged, "console_diagnostics", :compact)), OUTPUT_CONSOLE_DIAGNOSTICS_VALUES),
    console_q_limit_events = _validate_allowed_symbol("output.console_q_limit_events", _as_symbol_cfg(_raw_get(merged, "console_q_limit_events", :summary)), OUTPUT_CONSOLE_Q_LIMIT_EVENTS_VALUES),
    console_max_rows = _as_int_cfg(_raw_get(merged, "console_max_rows", 100)),
    logfile_results = _validate_allowed_symbol("output.logfile_results", _as_symbol_cfg(_raw_get(merged, "logfile_results", :off)), OUTPUT_LOGFILE_RESULTS_VALUES),
    result_table_max_rows = _as_int_cfg(_raw_get(merged, "result_table_max_rows", 200)),
    result_table_large_case_threshold_buses = _as_int_cfg(_raw_get(merged, "result_table_large_case_threshold_buses", 1000)),
    result_table_large_case_mode = _validate_allowed_symbol("output.result_table_large_case_mode", _as_symbol_cfg(_raw_get(merged, "result_table_large_case_mode", :summary)), OUTPUT_RESULT_TABLE_LARGE_CASE_MODE_VALUES),
    detailed_result_csv_write_mode = _validate_allowed_symbol("output.detailed_result_csv_write_mode", _as_symbol_cfg(_raw_get(merged, "detailed_result_csv_write_mode", :auto)), OUTPUT_DETAILED_RESULT_CSV_WRITE_MODE_VALUES),
    detailed_result_csv_exporter = _validate_allowed_symbol("output.detailed_result_csv_exporter", _as_symbol_cfg(_raw_get(merged, "detailed_result_csv_exporter", :auto)), OUTPUT_DETAILED_RESULT_CSV_EXPORTER_VALUES),
    detailed_result_csv_direct_threshold_buses = _output_positive_or_default(_as_int_cfg(_raw_get(merged, "detailed_result_csv_direct_threshold_buses", 10_000)), 10_000),
    detailed_result_csv_buffer_initial_bytes = _output_nonnegative_or_default(_as_int_cfg(_raw_get(merged, "detailed_result_csv_buffer_initial_bytes", 8 * 1024 * 1024)), 8 * 1024 * 1024),
    detailed_result_csv_buffer_max_bytes = _output_positive_or_default(_as_int_cfg(_raw_get(merged, "detailed_result_csv_buffer_max_bytes", 64 * 1024 * 1024)), 64 * 1024 * 1024),
    detailed_result_csv_streaming_threshold_rows = _output_positive_or_default(_as_int_cfg(_raw_get(merged, "detailed_result_csv_streaming_threshold_rows", 100_000)), 100_000),
    logfile_diagnostics = _validate_allowed_symbol("output.logfile_diagnostics", _as_symbol_cfg(_raw_get(merged, "logfile_diagnostics", :compact)), OUTPUT_LOGFILE_DIAGNOSTICS_VALUES),
    logfile_performance = _validate_allowed_symbol("output.logfile_performance", _as_symbol_cfg(_raw_get(merged, "logfile_performance", :compact)), OUTPUT_LOGFILE_PERFORMANCE_VALUES),
    logfile_warnings = _validate_allowed_symbol("output.logfile_warnings", _as_symbol_cfg(_raw_get(merged, "logfile_warnings", :table)), OUTPUT_LOGFILE_WARNINGS_VALUES),
  )
end

function ControlConfig(raw::AbstractDict)
  merged = _merged_section(raw, "control")
  raw_controllers = _raw_get(merged, "controllers", Any[])
  raw_controllers isa AbstractVector || throw(ArgumentError("control.controllers must be a vector."))
  return ControlConfig(
    enabled = _as_bool_cfg(_raw_get(merged, "enabled", true)),
    max_outer_iterations = _as_int_cfg(_raw_get(merged, "max_outer_iterations", 20)),
    trace = _as_bool_cfg(_raw_get(merged, "trace", true)),
    log_iterations = _as_bool_cfg(_raw_get(merged, "log_iterations", true)),
    stop_on_pf_failure = _as_bool_cfg(_raw_get(merged, "stop_on_pf_failure", true)),
    controllers = Any[raw_controllers...],
  )
end

function SparlectraConfig(raw::AbstractDict)
  return SparlectraConfig(
    powerflow = PowerFlowConfig(raw),
    state_estimation = StateEstimationConfig(raw),
    matpower = MatpowerImportConfig(raw),
    transformer = TransformerConfig(raw),
    performance = PerformanceConfig(raw),
    benchmark = BenchmarkConfig(raw),
    runtime = RuntimeConfig(raw),
    diagnostics = DiagnosticsConfig(raw),
    output = OutputConfig(raw),
    control = ControlConfig(raw),
  )
end

_canonical_config_key(key::AbstractString)::String = key == "powerflow" ? "power_flow" : key == "matpower" ? "matpower_import" : String(key)

function _merge_config_overrides(raw::Dict{String,Any}, overrides::AbstractDict)
  merged = deepcopy(raw)
  for (key, value) in overrides
    skey = _canonical_config_key(_config_key(key))
    if value isa AbstractDict
      section = get!(merged, skey, Dict{String,Any}())
      section isa AbstractDict || throw(ArgumentError("Cannot merge override section $(repr(skey)) into scalar value."))
      merged[skey] = _merge_config_overrides(Dict{String,Any}(String(k) => v for (k, v) in section), value)
    else
      merged[skey] = value
    end
  end
  return merged
end

function _validate_known_config_keys(user::AbstractDict, defaults::AbstractDict; path::String = "")
  for (key, value) in user
    skey = _canonical_config_key(_config_key(key))
    isempty(path) && skey in ("_config_sources", "_config_metadata") && continue
    current_path = isempty(path) ? skey : string(path, ".", skey)
    if current_path == "matpower_import.benchmark"
      throw(ArgumentError("Removed Sparlectra configuration key: matpower_import.benchmark.\nUse top-level benchmark.enabled instead, e.g.\n\nbenchmark:\n  enabled: true"))
    end
    if isempty(path) && skey == "methods"
      for m in _as_symbol_vector_cfg(value)
        requested = _powerflow_method_label(m)
        requested === SUPPORTED_POWERFLOW_METHOD || throw(ArgumentError(unsupported_powerflow_method_message(requested)))
      end
      continue
    end
    haskey(defaults, skey) || throw(ArgumentError("Unknown Sparlectra configuration key: $(current_path)"))
    if value isa AbstractDict
      defaults[skey] isa AbstractDict || throw(ArgumentError("Configuration key $(current_path) must be a scalar value."))
      _validate_known_config_keys(value, defaults[skey]; path = current_path)
    end
  end
  return nothing
end

_config_file_hash(path::AbstractString) = isfile(path) ? bytes2hex(sha256(read(path))) : ""
_config_file_mtime(path::AbstractString) = isfile(path) ? stat(path).mtime : 0.0

function _load_and_validate_config(default_path::AbstractString, user_path::AbstractString; cli_overrides::AbstractDict, overrides::AbstractDict)
  isfile(default_path) || throw(ArgumentError("Default Sparlectra config file not found: $(default_path)"))
  if abspath(user_path) != abspath(USER_SPARLECTRA_CONFIG_PATH) && !isfile(user_path)
    throw(ArgumentError("Sparlectra user config file not found: $(user_path)"))
  end
  defaults = load_yaml_dict(default_path)
  user = isfile(user_path) ? load_yaml_dict(user_path) : Dict{String,Any}()
  _validate_known_config_keys(user, defaults)
  _validate_known_config_keys(cli_overrides, defaults)
  _validate_known_config_keys(overrides, defaults)
  raw = merge_yaml_dict!(deepcopy(defaults), Dict{String,Any}(String(_canonical_config_key(_config_key(k))) => v for (k, v) in user))
  raw = _merge_config_overrides(raw, cli_overrides)
  raw = _merge_config_overrides(raw, overrides)
  _validate_powerflow_solver_config!(raw)
  return raw, isfile(user_path)
end

"""
    load_sparlectra_config([user_path]; reload=false, cli_overrides=Dict(), overrides=Dict())

Load, validate, merge, and cache the central typed `SparlectraConfig`.
Configuration precedence is `src/configuration.yaml.example`, optional
`examples/configuration.yaml` (or an explicit `user_path`), CLI-style overrides,
then explicit Julia API overrides. Unknown user keys throw an `ArgumentError`.
"""

"""
    refresh_sparlectra_config_file(path; write=false, backup=true, normalize_deprecated=true, default_path=DEFAULT_SPARLECTRA_CONFIG_PATH)

Compare a user YAML configuration with the current Sparlectra template, add
missing keys from the template, and optionally normalize documented deprecated
aliases. Dry-run mode returns the refreshed YAML without writing. Writes are
explicit and create a timestamped backup unless `backup=false` is requested.
"""
function refresh_sparlectra_config_file(path::AbstractString; write::Bool = false, backup::Bool = true, normalize_deprecated::Bool = true, default_path::AbstractString = DEFAULT_SPARLECTRA_CONFIG_PATH)
  isfile(path) || throw(ArgumentError("Sparlectra config file not found: $(path)"))
  isfile(default_path) || throw(ArgumentError("Default Sparlectra config file not found: $(default_path)"))
  duplicate_keys = _detect_yaml_duplicate_keys(path)
  warnings = String[]
  for key in duplicate_keys
    push!(warnings, "Duplicate YAML key detected: $(key). Refresh did not write because duplicate keys require manual review.")
  end
  defaults = load_yaml_dict(default_path)
  user = load_yaml_dict(path)
  refreshed = deepcopy(user)
  missing_keys = String[]
  _add_missing_config_keys!(refreshed, defaults, missing_keys)
  normalized_keys = normalize_deprecated ? _normalize_deprecated_config_aliases!(refreshed) : String[]
  refreshed_text = _yaml_dict_text(refreshed)
  original_text = read(path, String)
  changed = !isempty(missing_keys) || !isempty(normalized_keys) || refreshed_text != original_text
  backup_path = nothing
  written = false
  success = true
  if write && changed
    if !isempty(duplicate_keys)
      success = false
    else
      if backup
        backup_path = string(path, ".bak-", Dates.format(now(), "yyyymmdd-HHMMSS"))
        cp(path, backup_path; force = false)
      end
      Base.write(path, refreshed_text)
      written = true
    end
  end
  return (;
    success,
    changed,
    written,
    backup_path,
    missing_keys = sort!(missing_keys),
    normalized_keys = sort!(normalized_keys),
    duplicate_keys = sort!(duplicate_keys),
    warnings,
    refreshed_text,
  )
end

function refresh_sparlectra_config_text(text::AbstractString; normalize_deprecated::Bool = true, default_path::AbstractString = DEFAULT_SPARLECTRA_CONFIG_PATH)
  mktemp() do path, io
    write(io, text)
    close(io)
    return refresh_sparlectra_config_file(path; write = false, backup = false, normalize_deprecated, default_path)
  end
end

function _detect_yaml_duplicate_keys(path::AbstractString)::Vector{String}
  duplicates = String[]
  stack_names = Dict{Int,Vector{String}}(0 => String[])
  seen = Dict{Int,Set{String}}(0 => Set{String}())
  open(path, "r") do io
    for (lineno, rawline) in enumerate(eachline(io))
      line = rstrip(_strip_yaml_comment(rawline))
      isempty(strip(line)) && continue
      indent_count = count(==(' '), first(line, something(findfirst(!=(' '), line), lastindex(line) + 1) - 1))
      indent_count % 2 == 0 || continue
      content = strip(line)
      m = match(r"^([A-Za-z0-9_\-]+)\s*:\s*(.*)$", content)
      isnothing(m) && continue
      level = indent_count ÷ 2
      key = String(m.captures[1])
      parent_path = get(stack_names, level, String[])
      level_seen = get!(seen, level, Set{String}())
      compound = isempty(parent_path) ? key : string(join(parent_path, "."), ".", key)
      if key in level_seen
        push!(duplicates, compound)
      else
        push!(level_seen, key)
      end
      if isempty(String(m.captures[2]))
        stack_names[level + 1] = [parent_path; key]
        seen[level + 1] = Set{String}()
      end
      for old_level in collect(keys(stack_names))
        old_level > level + 1 && delete!(stack_names, old_level)
      end
      for old_level in collect(keys(seen))
        old_level > level + 1 && delete!(seen, old_level)
      end
    end
  end
  return unique(duplicates)
end

function _add_missing_config_keys!(dst::Dict{String,Any}, defaults::Dict{String,Any}, missing::Vector{String}; prefix::String = "")
  for key in sort(collect(keys(defaults)))
    path = isempty(prefix) ? key : string(prefix, ".", key)
    if !haskey(dst, key)
      dst[key] = deepcopy(defaults[key])
      push!(missing, path)
    elseif dst[key] isa Dict{String,Any} && defaults[key] isa Dict{String,Any}
      _add_missing_config_keys!(dst[key], defaults[key], missing; prefix = path)
    end
  end
  return dst
end

function _ensure_config_section!(root::Dict{String,Any}, keys::AbstractVector{<:AbstractString})::Dict{String,Any}
  current = root
  for key in keys
    if !haskey(current, key) || !(current[key] isa Dict{String,Any})
      current[key] = Dict{String,Any}()
    end
    current = current[key]
  end
  return current
end

function _normalize_deprecated_config_aliases!(root::Dict{String,Any})::Vector{String}
  normalized = String[]
  pf = haskey(root, "power_flow") ? root["power_flow"] : get(root, "powerflow", nothing)
  pf isa Dict{String,Any} || return normalized
  start_mode = get(pf, "start_mode", nothing)
  if start_mode isa Dict{String,Any} && get(start_mode, "voltage_mode", nothing) in ("bus_vm_va_blend", :bus_vm_va_blend)
    start_mode["voltage_mode"] = "profile_blend"
    start_mode["profile_source"] = "matpower_reference"
    push!(normalized, "power_flow.start_mode.voltage_mode")
    push!(normalized, "power_flow.start_mode.profile_source")
  end
  qlimits = get(pf, "qlimits", nothing)
  if qlimits isa Dict{String,Any}
    mode = get(qlimits, "enforcement_mode", nothing)
    if mode in ("matpower_simultaneous", :matpower_simultaneous)
      qlimits["enforcement_mode"] = "classic_simultaneous"
      push!(normalized, "power_flow.qlimits.enforcement_mode")
    elseif mode in ("matpower_one_at_a_time", :matpower_one_at_a_time)
      qlimits["enforcement_mode"] = "classic_one_at_a_time"
      push!(normalized, "power_flow.qlimits.enforcement_mode")
    end
  end
  return unique(normalized)
end

function _yaml_scalar_text(value)
  value === nothing && return "null"
  value isa Bool && return value ? "true" : "false"
  value isa Symbol && return String(value)
  value isa AbstractString && return occursin(r"[:#\[\],{}]|^$", value) ? repr(String(value)) : String(value)
  value isa AbstractVector && return string("[", join((_yaml_scalar_text(v) for v in value), ", "), "]")
  return string(value)
end

function _write_yaml_dict(io::IO, data::Dict{String,Any}; indent::Int = 0)
  pad = repeat(" ", indent)
  for key in sort(collect(keys(data)))
    value = data[key]
    if value isa Dict{String,Any}
      println(io, pad, key, ":")
      _write_yaml_dict(io, value; indent = indent + 2)
    else
      println(io, pad, key, ": ", _yaml_scalar_text(value))
    end
  end
end

function _yaml_dict_text(data::Dict{String,Any})::String
  io = IOBuffer()
  _write_yaml_dict(io, data)
  return String(take!(io))
end
function load_sparlectra_config(user_path::AbstractString = USER_SPARLECTRA_CONFIG_PATH; default_path::AbstractString = DEFAULT_SPARLECTRA_CONFIG_PATH, reload::Bool = false, cli_overrides::AbstractDict = Dict{String,Any}(), overrides::AbstractDict = Dict{String,Any}())
  abspath_default = abspath(default_path)
  abspath_user = abspath(user_path)
  default_mtime = _config_file_mtime(abspath_default)
  user_mtime = _config_file_mtime(abspath_user)
  default_digest = _config_file_hash(abspath_default)
  user_digest = _config_file_hash(abspath_user)
  cached = _SPARLECTRA_CONFIG_CACHE[]
  if !reload &&
     cached !== nothing &&
     cached.default_path == abspath_default &&
     cached.user_path == abspath_user &&
     cached.default_mtime == default_mtime &&
     cached.user_mtime == user_mtime &&
     cached.default_hash == default_digest &&
     cached.user_hash == user_digest &&
     isempty(cli_overrides) &&
     isempty(overrides)
    return cached.config
  end

  raw, user_found = _load_and_validate_config(abspath_default, abspath_user; cli_overrides = cli_overrides, overrides = overrides)
  config = SparlectraConfig(raw)
  if isempty(cli_overrides) && isempty(overrides)
    _SPARLECTRA_CONFIG_CACHE[] = (; default_path = abspath_default, user_path = abspath_user, default_mtime, user_mtime, default_hash = default_digest, user_hash = user_digest, user_found, config)
  end
  return config
end


function load_sparlectra_config!(user_path::AbstractString = USER_SPARLECTRA_CONFIG_PATH; default_path::AbstractString = DEFAULT_SPARLECTRA_CONFIG_PATH, reload::Bool = false, cli_overrides::AbstractDict = Dict{String,Any}(), overrides::AbstractDict = Dict{String,Any}())
  return set_sparlectra_config!(load_sparlectra_config(user_path; default_path = default_path, reload = reload, cli_overrides = cli_overrides, overrides = overrides))
end

function _print_config_section(io::IO, name::AbstractString, value; indent::AbstractString = "")
  println(io, indent, name, ":")
  child_indent = string(indent, "  ")
  for field in fieldnames(typeof(value))
    field_value = getfield(value, field)
    if field_value isa StartModeConfig || field_value isa QLimitConfig || field_value isa ObservabilityConfig
      _print_config_section(io, String(field), field_value; indent = child_indent)
    else
      println(io, child_indent, field, ": ", field_value)
    end
  end
  return nothing
end

"""
    print_effective_config([io], config::SparlectraConfig)

Print an `Effective Sparlectra Configuration` block with typed module sections.
"""
function print_effective_config(io::IO, config::SparlectraConfig)
  println(io, "==================== Effective Sparlectra Configuration ====================")
  println(io, "default_file: ", DEFAULT_SPARLECTRA_CONFIG_PATH)
  println(io, "user_file: ", USER_SPARLECTRA_CONFIG_PATH)
  println(io, "user_file_found: ", isfile(USER_SPARLECTRA_CONFIG_PATH))
  _print_config_section(io, "power_flow", config.powerflow)
  _print_config_section(io, "state_estimation", config.state_estimation)
  _print_config_section(io, "matpower_import", config.matpower)
  _print_config_section(io, "transformer", config.transformer)
  _print_config_section(io, "performance", config.performance)
  _print_config_section(io, "runtime", config.runtime)
  _print_config_section(io, "diagnostics", config.diagnostics)
  _print_config_section(io, "output", config.output)
  println(io, "==========================================================================")
  return nothing
end

print_effective_config(config::SparlectraConfig) = print_effective_config(stdout, config)
