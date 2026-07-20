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

# Date: 2026-06-19
# file: examples/qlimit_large_case_mode_comparison.jl
# purpose: compares Q-limit large-case start-profile/enforcement-mode combinations across configured MATPOWER cases (experimental diagnostic tooling)

using Sparlectra
include(joinpath(@__DIR__, "experimental", "qlimit_large_case_comparison.jl"))
include(joinpath(@__DIR__, "internal", "example_header.jl"))

function _parse_args(args)
  cases = String[]
  start_profiles = collect(QLIMIT_LARGE_CASE_START_PROFILES)
  modes = collect(QLIMIT_LARGE_CASE_MODES)
  output_root = nothing
  verbose_logs = false
  profiles_overridden = false
  modes_overridden = false
  i = firstindex(args)
  while i <= lastindex(args)
    arg = args[i]
    if arg == "--classic-only"
      start_profiles = [:classic_start]
    elseif arg == "--robust-only"
      start_profiles = [:robust_dc_start]
    elseif arg == "--case"
      i == lastindex(args) && throw(ArgumentError("--case requires a case name or path"))
      i += 1
      push!(cases, args[i])
    elseif arg == "--profile"
      i == lastindex(args) && throw(ArgumentError("--profile requires a profile name"))
      i += 1
      profiles_overridden || (empty!(start_profiles); profiles_overridden = true)
      push!(start_profiles, Symbol(args[i]))
    elseif arg == "--mode"
      i == lastindex(args) && throw(ArgumentError("--mode requires a Q-limit mode"))
      i += 1
      modes_overridden || (empty!(modes); modes_overridden = true)
      push!(modes, Symbol(args[i]))
    elseif arg == "--full-logs"
      verbose_logs = true
    elseif arg == "--compact-logs"
      verbose_logs = false
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
  isempty(cases) && append!(cases, QLIMIT_LARGE_CASE_DEFAULTS)
  return (cases = cases, start_profiles = start_profiles, modes = modes, output_root = output_root, verbose_logs = verbose_logs)
end

function main(args = ARGS)
  print_example_banner("examples/qlimit_large_case_mode_comparison.jl", "compares Q-limit large-case start-profile/enforcement-mode combinations across configured MATPOWER cases (experimental diagnostic tooling)")
  parsed = _parse_args(args)
  kwargs = Dict{Symbol,Any}(
    :cases => parsed.cases,
    :start_profiles => parsed.start_profiles,
    :modes => parsed.modes,
    :verbose_logs => parsed.verbose_logs,
    :io => stdout,
    :verbose => true,
  )
  parsed.output_root === nothing || (kwargs[:output_root] = parsed.output_root)
  result = compare_qlimit_large_case_modes(; kwargs...)
  return result
end

# SPARLECTRA_QLIMIT_LARGE_CASE_NO_MAIN=1 is a test-only escape hatch: it lets
# test/test_qlimit_large_case_comparison.jl load this file's function
# definitions (e.g. _parse_args) without triggering the real, network/large
# case-dependent comparison run. It has no effect on normal console/
# interactive use — those always run main() unconditionally, exactly like
# every other examples/*.jl program (see
# examples/exp_transformer_tap_changer_model.jl for why no
# `abspath(PROGRAM_FILE) == @__FILE__` guard is used here).
if get(ENV, "SPARLECTRA_QLIMIT_LARGE_CASE_NO_MAIN", "0") != "1"
  run_example(main)
end
nothing

