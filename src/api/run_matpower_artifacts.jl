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

# Keep MATPOWER import artifacts outside the API orchestrator so the run path
# only decides when artifacts are needed, not how each artifact is written.

function _update_effective_matpower_raw!(raw::Dict{String,Any}, cfg::SparlectraConfig)
  mat_raw = get!(raw, "matpower_import", Dict{String,Any}())
  mat_raw isa Dict{String,Any} || return raw
  mat_raw["auto_profile"] = String(cfg.matpower.auto_profile)
  mat_raw["ratio"] = String(cfg.matpower.ratio)
  mat_raw["shift_sign"] = cfg.matpower.shift_sign
  mat_raw["shift_unit"] = String(cfg.matpower.shift_unit)
  mat_raw["bus_shunt_model"] = String(cfg.matpower.bus_shunt_model)
  mat_raw["pv_voltage_source"] = String(cfg.matpower.pv_voltage_source)
  mat_raw["compare_voltage_reference"] = String(cfg.matpower.compare_voltage_reference)
  return raw
end

function _write_matpower_auto_profile_artifact(output_path::AbstractString, profile, cfg::SparlectraConfig; casefile::AbstractString)::String
  artifact = joinpath(output_path, "matpower_auto_profile.log")
  open(artifact, "w") do io
    println(io, "MATPOWER import auto-profile artifact")
    println(io, "====================================")
    write_matpower_import_auto_profile(io, profile, cfg; casefile = casefile)
  end
  return artifact
end

