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

# `power_flow.start_mode.start_projection` is not in `GUI_EDITABLE_CONFIG_KEYS`,
# so the fixed-reference self-check settings below are applied via a merged,
# temporary YAML configuration file rather than run_sparlectra_api's
# `config_overrides` dict (which is validated against that allowlist).
#
# Every start-value machine is forced off so the imported bus voltages reach
# the solver verbatim, whatever the base configuration says:
# - flatstart=false: `initialVrect` keeps the imported vm/va instead of 1.0/0°;
# - start_projection=false: no DC angle start, no blend scan;
# - dc_seed_unconditional=false: no standalone DC pre-solve seeding;
# - start_current_iteration.enabled=false: no current-iteration pre-solve;
# - apslf_start.enabled=false: no APSLF start preconditioner.
# angle_mode/voltage_mode matter only on the MATPOWER import path
# (`_apply_matpower_start_modes!`) and only when flatstart is true; they are
# kept for backward compatibility of the written config, and with
# flatstart=false the MATPOWER path reaches the same case VM/VA start through
# `apply_mp_bus_vmva_init!`. On the CGMES path the SV voltages are already the
# imported bus state, so bypassing the start machinery is all that is needed.
function _self_check_forced_overrides()::Dict{String,Any}
  return Dict{String,Any}(
    "power_flow" => Dict{String,Any}(
      "max_iter" => 1,
      # `flatstart` is the legacy power_flow-level key — the only one the
      # config schema accepts — and PowerFlowConfig merges it into the
      # start_mode raw dict (configuration.jl), so this reliably forces
      # flatstart=false even when the base config sets `flatstart: true`
      # (observed in a WebUI configuration.yaml, where it silently wiped the
      # CGMES SV start on every run).
      "flatstart" => false,
      "start_mode" => Dict{String,Any}(
        "angle_mode" => "matpower_va",
        "voltage_mode" => "all_bus_vm",
        "start_projection" => false,
        "dc_seed_unconditional" => false,
      ),
      "start_current_iteration" => Dict{String,Any}("enabled" => false),
      "apslf_start" => Dict{String,Any}("enabled" => false),
      "qlimits" => Dict{String,Any}("enabled" => false),
    ),
  )
end

function _write_self_check_config_file(config_file::AbstractString)::String
  isfile(config_file) || throw(ArgumentError("Sparlectra configuration file not found: $(config_file)"))
  base = load_yaml_dict(config_file)
  merged = _merge_config_overrides(base, _self_check_forced_overrides())
  path = tempname() * "_self_check_config.yaml"
  _write_yaml_file(path, merged)
  return path
end

# The statement below is asserted by tests: it may only be written when every
# start-value machine really was forced off (which _self_check_forced_overrides
# guarantees for any base configuration).
const _SELF_CHECK_VERBATIM_LINE = "start values taken verbatim from import (no flat start, no DC seed, no start projection, no current-iteration pre-solve, no APSLF start)"

# Start-state residual over the captured per-bus rows: the max |P|/|Q| residual
# at the exact state Newton saw, i.e. the self-check answer independent of the
# one corrective step max_iter=1 still takes.
function _self_check_start_residual(rows)::Float64
  rows isa AbstractVector || return NaN
  worst = NaN
  for r in rows
    for v in (r.p_residual, r.q_residual)
      isfinite(v) || continue
      (isnan(worst) || abs(v) > worst) && (worst = abs(v))
    end
  end
  return worst
end

function _write_self_check_summary(result, output_dir::AbstractString)::Nothing
  raw = result.raw_result
  profile = raw === nothing ? nothing : raw.performance_profile
  rows = profile isa AbstractDict ? get(profile, :initial_residual_rows, nothing) : nothing
  path = joinpath(output_dir, "self_check.log")
  open(path, "w") do io
    println(io, "# Fixed-reference self-check")
    println(io, "forced: max_iter=1, qlimits.enabled=false, flatstart=false, start_mode.start_projection=false, start_mode.dc_seed_unconditional=false, start_current_iteration.enabled=false, apslf_start.enabled=false")
    if raw === nothing
      # The run never reached the solver — no start state existed, so the
      # verbatim-start statement would be an empty claim.
      println(io, "run did not reach the solver (", something(result.reason, "unknown"), "); no residual was evaluated")
      return
    end
    println(io, _SELF_CHECK_VERBATIM_LINE)
    println(io, "note: the slack bus magnitude is pinned to its regulating setpoint (vset) by the solver; a vset that deviates from the imported slack voltage shows up as residual at slack-adjacent buses")
    start_residual = _self_check_start_residual(rows)
    isnan(start_residual) || println(io, "start_state_residual_inf: ", start_residual)
    rows isa AbstractVector && println(io, "residual_rows: ", length(rows), " (self_check_residuals.csv)")
    if raw !== nothing
      # GUI-editable overrides (Web UI form: max_iter, qlimits.enabled) are
      # applied on top of the forced self-check config, so the run may have
      # taken more than the forced single iteration — report what actually
      # ran instead of assuming one step. start_state_residual_inf above is
      # captured before the first iteration and stays the fixed-reference
      # answer either way.
      println(io, "iterations_run: ", raw.iterations)
      println(io, "final_mismatch_after_solve: ", raw.final_mismatch)
      println(io, "outcome: ", raw.outcome)
      raw.iterations isa Integer && raw.iterations > 1 && println(io, "note: GUI/API overrides allowed more than the forced single iteration (max_iter/qlimits are form-editable); the start-state residual above is unaffected")
    end
    no_sv = get(result.metadata, "cgmes_no_sv_buses", nothing)
    if no_sv !== nothing
      buses = get(result.metadata, "cgmes_buses", nothing)
      println(io, "cgmes_no_sv_buses: ", no_sv, buses === nothing ? "" : " of $(buses)", " (flat 1.0 pu / 0° start — a large share weakens the SV comparison)")
    end
  end
  return nothing
end

"""
    run_fixed_reference_self_check(; casefile, config_file=DEFAULT_SPARLECTRA_CONFIG_PATH,
                                    output_dir, config_overrides=Dict(), case_format=:auto) -> SparlectraApiResult

Evaluate the Newton-Raphson mismatch at a case's own stored operating point,
without letting any start-value machinery adjust it. Where a normal run
answers "does this case converge with the configured start/step-control
settings?", the self-check answers a narrower, model-focused question: is the
case's own recorded operating point already close to power balance under
Sparlectra's imported network model? A large residual here points at the
imported network model (branch parameters, shifted angles, a wrong per-unit
convention) rather than at the solver's start guess or step control.

Supported inputs:
- **MATPOWER** (`case_format = :matpower`/`:auto`): evaluates at the case's
  stored `VM`/`VA` columns.
- **CGMES** (`case_format = :cgmes`): evaluates at the delivery's `SvVoltage`
  state, which the importer writes into the bus objects at creation time.
  Buses without a usable `SvVoltage` keep the importer's `1.0 pu / 0°`
  fallback; their count is reported in `self_check.log` and per bus in the
  `has_sv` column of `self_check_residuals.csv` — a large no-SV share weakens
  the interpretation of the residual.

# Arguments
- `casefile`: case path (MATPOWER file or CGMES delivery), forwarded to
  [`run_sparlectra_api`](@ref).
- `config_file`: base configuration; the self-check settings are merged on top
  of it (see "How" below), so any other settings in `config_file` (e.g.
  import conventions) still apply.
- `output_dir`: forwarded to `run_sparlectra_api`; all normal run artifacts
  (including `diagnose.log`) are written here, plus `self_check.log` and
  `self_check_residuals.csv`.
- `config_overrides`: optional further `GUI_EDITABLE_CONFIG_KEYS` overrides,
  applied on top of the self-check settings (e.g. to adjust `power_flow.tol`
  for the reported mismatch classification).

# How
Forces `power_flow.max_iter = 1`, `power_flow.qlimits.enabled = false`, and
every start-value machine off (`start_mode.flatstart = false`,
`start_mode.start_projection = false`, `start_mode.dc_seed_unconditional =
false`, `start_current_iteration.enabled = false`, `apslf_start.enabled =
false`), then runs through [`run_sparlectra_api`](@ref) with diagnostics
enabled. The imported bus voltages therefore reach the solver verbatim; the
only remaining adjustment is the solver pinning the slack magnitude to its
regulating setpoint (noted in `self_check.log`).

# Returns
A `SparlectraApiResult`. `success` reflects whether the run completed, not
whether the residual is small — read `self_check.log`'s
`start_state_residual_inf` (the max |P|/|Q| residual at the unmodified start
state) or `raw_result.final_mismatch` (after the single corrective step) for
that. `self_check_residuals.csv` holds the full per-bus residual export for
attribution.
"""
function run_fixed_reference_self_check(;
  casefile::AbstractString,
  config_file::AbstractString = DEFAULT_SPARLECTRA_CONFIG_PATH,
  output_dir::AbstractString,
  config_overrides::AbstractDict = Dict{String,Any}(),
  case_format = :auto,
)::SparlectraApiResult
  self_check_config_file = _write_self_check_config_file(config_file)
  try
    result = run_sparlectra_api(
      casefile = casefile,
      config_file = self_check_config_file,
      output_dir = output_dir,
      case_format = case_format,
      config_overrides = config_overrides,
      run_diagnostics = true,
    )
    try
      _write_self_check_summary(result, result.output_dir)
    catch err
      @warn "could not write self_check.log" exception = err
    end
    return result
  finally
    rm(self_check_config_file; force = true)
  end
end
