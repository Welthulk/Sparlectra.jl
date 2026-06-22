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
# file: examples/qlimit_large_case_mode_comparison.jl
using Sparlectra

function _parse_args(args)
  cases = String[]
  start_profiles = collect(Sparlectra.QLIMIT_LARGE_CASE_START_PROFILES)
  output_root = nothing
  full_logs = false
  i = firstindex(args)
  while i <= lastindex(args)
    arg = args[i]
    if arg == "--classic-only"
      start_profiles = [:classic_start]
    elseif arg == "--robust-only"
      start_profiles = [:robust_dc_start]
    elseif arg == "--full-logs"
      full_logs = true
    elseif arg == "--output-root"
      i == lastindex(args) && throw(ArgumentError("--output-root requires a path argument"))
      i += 1
      output_root = args[i]
    elseif startswith(arg, "--")
      throw(ArgumentError("unknown option: $(arg)"))
    else
      push!(cases, arg)
    end
    i += 1
  end
  isempty(cases) && append!(cases, Sparlectra.QLIMIT_LARGE_CASE_DEFAULTS)
  return (cases = cases, start_profiles = start_profiles, output_root = output_root, full_logs = full_logs)
end

function main(args = ARGS)
  parsed = _parse_args(args)
  kwargs = Dict{Symbol,Any}(
    :cases => parsed.cases,
    :start_profiles => parsed.start_profiles,
    :verbose_logs => parsed.full_logs,
    :io => stdout,
    :verbose => true,
  )
  parsed.output_root === nothing || (kwargs[:output_root] = parsed.output_root)
  result = Sparlectra.compare_qlimit_large_case_modes(; kwargs...)
  println("Q-limit large-case mode comparison complete")
  println("rows: ", length(result.rows))
  println("csv : ", result.csv_path)
  println("json: ", result.json_path)
  for row in result.rows
    println(row["case"], " | ", row["start_profile"], " | ", row["mode"], " | ", row["run_status"], " | ", row["numerical_status"], " | ", row["run_directory"])
  end
  return result
end

if get(ENV, "SPARLECTRA_QLIMIT_LARGE_CASE_NO_MAIN", "0") != "1"
  Base.invokelatest(main)
end
