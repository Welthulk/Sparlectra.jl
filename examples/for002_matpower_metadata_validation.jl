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

using Sparlectra

"""
    main(; casefile, output_dir = "for002_metadata_artifacts")

Import a FOR001/FOR002 MATPOWER fixture with metadata enabled and write compact
CSV/Markdown artifacts that validation tooling can compare against Builder and
FOR002 exports. FOR002 slack labels such as `BSTADTS1 S` should be normalized by
removing only the explicit trailing slack suffix before comparing them with
`mpc.bus_name` values such as `BSTADTS1`.
"""
function main(; casefile::AbstractString, output_dir::AbstractString = "for002_metadata_artifacts")
  mpc = Sparlectra.MatpowerIO.read_case(casefile; legacy_compat = true)
  net = Sparlectra.createNetFromMatPowerCase(
    mpc = mpc,
    apply_bus_names = true,
    apply_branch_names = true,
    apply_branch_kind = true,
    import_for001_contingencies = true,
  )
  mkpath(output_dir)
  contingency_indices = Sparlectra.MatpowerIO.for001_contingency_branch_indices(mpc)

  open(joinpath(output_dir, "matpower_metadata_summary.md"), "w") do io
    println(io, "# FOR002 MATPOWER metadata validation summary")
    println(io)
    println(io, "- Case: `", casefile, "`")
    println(io, "- Imported buses: ", length(net.nodeVec))
    println(io, "- Imported branches: ", length(net.branchVec))
    println(io, "- FOR001 contingencies: ", length(net.for001Contingencies))
    println(io)
    println(io, "Compare Builder vs FOR002, MATPOWER import vs FOR002, and MATPOWER import vs Builder using the bus names from `mpc.bus_name` and branch outage mapping from `mpc.branch_name`/`mpc.for001_contingencies`; do not reconstruct names from numeric IDs.")
    println(io, "Normalize FOR002 slack labels by stripping a final ` S` suffix only when it denotes the slack marker, e.g. `BSTADTS1 S` -> `BSTADTS1`.")
  end

  open(joinpath(output_dir, "for001_contingencies.csv"), "w") do io
    println(io, "contingency,branch_index")
    for (name, idx) in zip(net.for001Contingencies, contingency_indices)
      println(io, "\"", replace(name, "\"" => "\"\""), "\",", idx)
    end
  end
  return (net = net, mpc = mpc, contingency_indices = contingency_indices)
end

if abspath(PROGRAM_FILE) == @__FILE__
  isempty(ARGS) && error("usage: julia --project=. examples/for002_matpower_metadata_validation.jl CASEFILE [OUTPUT_DIR]")
  Base.invokelatest(main; casefile = ARGS[1], output_dir = length(ARGS) >= 2 ? ARGS[2] : "for002_metadata_artifacts")
end
