# Copyright 2023-2026 Udo Schmitz
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

# file: examples/others/exp_powsybl_iidm_import.jl
# purpose: a PowSyBl IIDM file (.xiidm) imported in Julia, without Python
#          and without a table bundle: the shipped ieee14.xiidm through
#          run_sparlectra as any other case file, then micro_grid_be.xiidm
#          through the reader and the builder with the import report
#          (node-breaker views, tap steps, a file state the reader drops),
#          both solved and compared with the OpenLoadFlow reference that
#          ships next to each file. Theory: powsybl_import.md.

using Sparlectra
using DelimitedFiles

include(joinpath(@__DIR__, "..", "others", "example_header.jl"))

const POWSYBL_DIR = joinpath(dirname(dirname(@__DIR__)), "data", "powsybl")

# The OpenLoadFlow reference next to the file: bus-breaker id => v_pu.
function reference_voltages(case::String)
  data, header = readdlm(joinpath(POWSYBL_DIR, case * ".powsybl", "reference_buses.csv"), ',', String; quotes = true, header = true)
  cols = Dict(Symbol(strip(h)) => k for (k, h) in enumerate(vec(header)))
  return Dict(String(data[i, cols[:bus_breaker_id]]) => parse(Float64, data[i, cols[:v_pu]]) for i in axes(data, 1))
end

# The run configuration that reproduces OpenLoadFlow's defaults: every
# synchronous component solved, the mismatch shared through the imported
# participation factors, limit switching without hysteresis.
function olf_like_config()
  pf = PowerFlowConfig(max_iter = 60, tol = 1e-9, islands_enabled = true, qlimits = Sparlectra.QLimitConfig(hysteresis_pu = 1e-6), distributed_slack = Sparlectra.DistributedSlackConfig(enabled = true, p_mode = :imported))
  return SparlectraConfig(powerflow = pf, output = OutputConfig(logfile_results = :off, startup_latency_hint = false), control = ControlConfig())
end

function worst_deviation(net::Net, reference::Dict{String,Float64})
  worst = 0.0
  for (bus, v_ref) in reference
    idx = get(net.busDict, bus, nothing)
    idx === nothing && continue
    worst = max(worst, abs(net.nodeVec[idx]._vm_pu - v_ref))
  end
  return worst
end

"""
    main()

Two shipped IIDM files, read in Julia. `ieee14.xiidm` goes through
`run_sparlectra(casefile = ...)` like a MATPOWER case: the format is
detected from the file's first bytes, the reader builds the tables, the
builder the network. `micro_grid_be.xiidm` goes through the reader and
the builder by hand so the import report can be shown: the bus views of
the CGMES-derived substations, the tap steps of the transformers, the
remote voltage regulation as an outer-loop control, and the notice that
the file's bus state was solved at other tap positions and is not used as
the start. Both results are compared with the OpenLoadFlow voltages that
ship next to the files.
"""
function main()
  print_example_banner("examples/others/exp_powsybl_iidm_import.jl", "PowSyBl IIDM files imported in Julia, no Python: ieee14 and micro_grid_be against the OpenLoadFlow reference")

  # 1. ieee14: the file as a case file
  ieee14 = joinpath(POWSYBL_DIR, "ieee14.powsybl", "ieee14.xiidm")
  println("1. run_sparlectra on ", relpath(ieee14, dirname(dirname(@__DIR__))))
  res = run_sparlectra(casefile = ieee14, config = olf_like_config())
  res.numerical_converged || error("ieee14 did not converge")
  println("   converged in ", res.iterations, " iteration(s), ", length(res.net.nodeVec), " buses, ", length(res.net.branchVec), " branches")
  println("   worst |dV| against OpenLoadFlow: ", round(worst_deviation(res.net, reference_voltages("ieee14")); sigdigits = 3), " pu")
  println()

  # 2. micro_grid_be: reader, builder, report
  micro = joinpath(POWSYBL_DIR, "micro_grid_be.powsybl", "micro_grid_be.xiidm")
  println("2. read_iidm_tables + build_net_from_powsybl on ", relpath(micro, dirname(dirname(@__DIR__))))
  tables = Sparlectra.read_iidm_tables(micro)
  println("   IIDM version ", tables.manifest["iidm_version"], ", ", length(tables.buses.id), " bus-view buses, ", length(tables.bus_breaker_view_buses.id), " bus-breaker buses, ",
    length(tables.two_windings_transformers.id), " two-winding and ", length(tables.three_windings_transformers.id), " three-winding transformer(s), ",
    length(tables.ratio_tap_changers.id), " ratio and ", length(tables.phase_tap_changers.id), " phase tap changer(s)")
  # the generators regulate remote buses; OpenLoadFlow holds those buses,
  # the remote mode attaches the same as outer-loop machine controls
  net, report = build_net_from_powsybl(tables, PowsyblAdapterOptions(remote_regulation = :remote))
  println(format_powsybl_report(report))
  res2 = run_sparlectra(net = net, config = olf_like_config())
  res2.numerical_converged || error("micro_grid_be did not converge")
  println("   converged in ", res2.iterations, " iteration(s)")
  println("   worst |dV| against OpenLoadFlow: ", round(worst_deviation(net, reference_voltages("micro_grid_be")); sigdigits = 3), " pu")
  println()
  println("Reading the numbers:")
  println("  both files land within the bands of the test suite (ieee14 1e-6 pu, micro_grid_be 1e-5 pu) without Python;")
  println("  the micro_grid_be report names the dropped file state: its tap changers carry solvedTapPosition values")
  println("  that differ from tapPosition, so the state belongs to another network and the import starts flat.")
  return (ieee14 = res, micro_grid_be = res2)
end

run_example(main)
