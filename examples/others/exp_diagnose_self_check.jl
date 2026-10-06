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

# Date: 2026-09-29
# file: examples/others/exp_diagnose_self_check.jl
# purpose: demonstrates run_fixed_reference_self_check and the narrative diagnose.log report on a case with a deliberately shifted branch reactance
#
# how to run
#
#   SparlectraApp is the application package in the folder app/ (the folder
#   is called `app`, not `SparlectraApp`). This example needs that
#   environment to be the active one; with the library environment
#   (--project=.) `using SparlectraApp` fails.
#
#   Shell, in the repository root:
#     julia --project=app examples/others/exp_diagnose_self_check.jl
#   first time only:
#     julia --project=app -e "using Pkg; Pkg.instantiate()"
#
#   Julia REPL, started in the repository root:
#     ] activate app        the answer must read "Activating project at ...",
#                           "Activating new project" means a wrong folder name
#     ] instantiate         first time only
#     (Backspace leaves the package mode)
#     include("examples/others/exp_diagnose_self_check.jl")
#   the same without the package mode:
#     using Pkg; Pkg.activate("app"); Pkg.instantiate()
#
#   VS Code: choose the environment `app` in the status bar (Julia env), then
#   run the file.

using Sparlectra

# run_fixed_reference_self_check and the diagnose.log report live in
# SparlectraApp (app/src/api), not in the library. Without the app
# environment `using SparlectraApp` ends in a stack trace; this says what to do.
if Base.find_package("SparlectraApp") === nothing
    error("""
      SparlectraApp not found in the active environment ($(Base.active_project())).
      SparlectraApp is the application package in the folder app/ (the folder is called app,
      not SparlectraApp). Activate that environment from the repository root:

        Julia REPL:  ] activate app         (the answer must not say "new project")
                     ] instantiate          (first time only)
                     include("examples/others/exp_diagnose_self_check.jl")
        or:          using Pkg; Pkg.activate("app"); Pkg.instantiate()
        shell:       julia --project=app examples/others/exp_diagnose_self_check.jl
        VS Code:     choose the environment app in the status bar
      """)
end
using SparlectraApp

include(joinpath(@__DIR__, "..", "others", "example_header.jl"))

# Deliberately mis-scale one branch reactance so the case's own stored VM/VA
# is no longer close to power balance under Sparlectra's imported model —
# this is exactly the situation the self-check is meant to surface.
function _write_perturbed_case(source_path::AbstractString, target_path::AbstractString)
    lines = readlines(source_path)
    branch_section = false
    perturbed = false
    open(target_path, "w") do io
        for line in lines
            if occursin("mpc.branch", line)
                branch_section = true
            elseif branch_section && !perturbed && !isempty(strip(line)) && !startswith(strip(line), "%") && !startswith(strip(line), "]")
                fields = split(line, ';'; limit=2)
                cols = split(strip(fields[1]))
                if length(cols) >= 4
                    cols[4] = string(50 * parse(Float64, cols[4])) # x_pu on branch 1: 50x its imported value, well past normal line-to-line variance
                    fields[1] = "\t" * join(cols, "\t")
                    line = join(fields, ";")
                    perturbed = true
                end
            elseif branch_section && startswith(strip(line), "]")
                branch_section = false
            end
            println(io, line)
        end
    end
    return perturbed
end

# the shipped sp_case118 (data/mpower); examples never download a case
# (another case works too: fetch one yourself, e.g. with
# Sparlectra.ensure_casefile("case14.m"), and use its path here)
function main(; casefile::AbstractString=joinpath(pkgdir(Sparlectra), "data", "mpower", "sp_case118.m"), output_root::AbstractString=joinpath(@__DIR__, "..", "_out", "diagnose_self_check"))
    print_example_banner("examples/others/exp_diagnose_self_check.jl", "demonstrates run_fixed_reference_self_check and the narrative diagnose.log report on a case with a deliberately shifted branch reactance")
    source_path = casefile
    mkpath(output_root)
    perturbed_path = joinpath(output_root, "sp_case118_perturbed_branch.m")
    _write_perturbed_case(source_path, perturbed_path)

    self_check = run_fixed_reference_self_check(
        casefile=perturbed_path,
        output_dir=joinpath(output_root, "self_check"),
    )
    println("self-check run ID:      ", self_check.run_id)
    println("self-check residual:    ", self_check.raw_result.final_mismatch)
    println("worst-mismatch bus:     ", get(self_check.raw_result.diagnostics, :worst_mismatch_bus_id, "n/a"))
    println("worst-mismatch equation:", get(self_check.raw_result.diagnostics, :worst_mismatch_equation, "n/a"))
    println()
    println("diagnose.log:")
    println("-------------")
    for artifact in self_check.artifacts
        if artifact.name == "diagnose.log"
            println(read(artifact.path, String))
        end
    end
    return self_check
end

run_example(main)
