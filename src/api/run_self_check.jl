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
function _self_check_forced_overrides()::Dict{String,Any}
  return Dict{String,Any}(
    "power_flow" => Dict{String,Any}(
      "max_iter" => 1,
      "start_mode" => Dict{String,Any}(
        "angle_mode" => "matpower_va",
        "voltage_mode" => "all_bus_vm",
        "start_projection" => false,
      ),
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

"""
    run_fixed_reference_self_check(; casefile, config_file=DEFAULT_SPARLECTRA_CONFIG_PATH,
                                    output_dir, config_overrides=Dict(), case_format=:auto) -> SparlectraApiResult

Evaluate the Newton-Raphson mismatch at a MATPOWER case's own stored `VM`/`VA`
values, without taking any corrective Newton step. Where a normal run answers
"does this case converge with the configured start/step-control settings?",
the self-check answers a narrower, model-focused question: is the case's own
recorded operating point already close to power balance under Sparlectra's
imported network model? A large residual here points at the imported network
model (branch parameters, shifted angles, a wrong per-unit convention) rather
than at the solver's start guess or step control.

# Arguments
- `casefile`: MATPOWER case path, forwarded to [`run_sparlectra_api`](@ref).
- `config_file`: base configuration; the self-check settings are merged on top
  of it (see "How" below), so any other settings in `config_file` (e.g.
  `matpower_import` conventions) still apply.
- `output_dir`: forwarded to `run_sparlectra_api`; all normal run artifacts
  (including `diagnose.log`) are written here.
- `config_overrides`: optional further `GUI_EDITABLE_CONFIG_KEYS` overrides,
  applied on top of the self-check settings (e.g. to adjust `power_flow.tol`
  for the reported mismatch classification).

# How
Forces `power_flow.max_iter = 1`, `power_flow.start_mode.angle_mode =
:matpower_va`, `power_flow.start_mode.voltage_mode = :all_bus_vm`,
`power_flow.start_mode.start_projection = false`, and
`power_flow.qlimits.enabled = false`, then runs through
[`run_sparlectra_api`](@ref) exactly like a normal run. The returned
`SparlectraApiResult`'s `raw_result.final_mismatch` (and
`raw_result.diagnostics.initial_mismatch`, since only one iteration runs) is
the self-check residual.

# Returns
A `SparlectraApiResult`. `success` reflects whether the run completed, not
whether the residual is small — read `raw_result.final_mismatch` for that.
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
    return run_sparlectra_api(
      casefile = casefile,
      config_file = self_check_config_file,
      output_dir = output_dir,
      case_format = case_format,
      config_overrides = config_overrides,
      run_diagnostics = true,
    )
  finally
    rm(self_check_config_file; force = true)
  end
end
