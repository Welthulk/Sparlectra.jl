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

# file: examples/others/exp_island_reference_priority.jl
# purpose: reference priority (ProSumer.referencePriority, CGMES semantics) on
#          the shipped two-island case data/scf/two_islands_prio.scf.json:
#          which bus each AC island takes as its reference and why, with and
#          without stated priorities, in the power flow and in N-1.

using Sparlectra

include(joinpath(@__DIR__, "..", "others", "example_header.jl"))

"""
    build_case() -> Net

The network behind `data/scf/two_islands_prio.scf.json` (the file is
`exportSCF(build_case(); case_name = "two_islands_prio")`). Seven 110 kV
buses; the tie `Kesselau_110`-`Dornberg_110` is open, so the case runs as two
AC islands:

- north: the grid infeed at `Almhof_110` is the slack (reference priority 1),
  the 60 MW machine at `Brunntal_110` carries priority 2, the 200 MW machine
  at `Kesselau_110` none;
- south: no slack at all. `Dornberg_110` (40 MW, priority 2), `Eschwald_110`
  (150 MW, no priority, the largest unit) and `Fennwiese_110` (60 MW,
  priority 1) all regulate their voltage, `Gruenau_110` carries the load.
"""
function build_case()
  net = Net(name = "two_islands_prio", baseMVA = 100.0)
  for bus in ("Almhof_110", "Brunntal_110", "Kesselau_110", "Dornberg_110", "Eschwald_110", "Fennwiese_110", "Gruenau_110")
    addBus!(net = net, busName = bus, vn_kV = 110.0)
  end
  line!(a, b; status = 1) = addPIModelACLine!(net = net, fromBus = a, toBus = b, r_pu = 0.01, x_pu = 0.06, b_pu = 0.02, status = status)
  line!("Almhof_110", "Brunntal_110")
  line!("Brunntal_110", "Kesselau_110")
  line!("Almhof_110", "Kesselau_110")
  line!("Kesselau_110", "Dornberg_110"; status = 0)
  line!("Dornberg_110", "Eschwald_110")
  line!("Eschwald_110", "Fennwiese_110")
  line!("Fennwiese_110", "Gruenau_110")
  line!("Gruenau_110", "Dornberg_110")
  addProsumer!(net = net, busName = "Almhof_110", type = "EXTERNALNETWORKINJECTION", referencePri = "Almhof_110", vm_pu = 1.02, va_deg = 0.0, referencePriority = 1)
  addProsumer!(net = net, busName = "Brunntal_110", type = "SYNCHRONOUSMACHINE", p = 40.0, q = 0.0, pMin = 0.0, pMax = 60.0, qMin = -30.0, qMax = 30.0, vm_pu = 1.01, isRegulated = true, referencePriority = 2)
  addProsumer!(net = net, busName = "Brunntal_110", type = "ENERGYCONSUMER", p = 50.0, q = 10.0)
  addProsumer!(net = net, busName = "Kesselau_110", type = "SYNCHRONOUSMACHINE", p = 80.0, q = 0.0, pMin = 0.0, pMax = 200.0, qMin = -80.0, qMax = 80.0, vm_pu = 1.0, isRegulated = true)
  addProsumer!(net = net, busName = "Kesselau_110", type = "ENERGYCONSUMER", p = 90.0, q = 20.0)
  addProsumer!(net = net, busName = "Dornberg_110", type = "SYNCHRONOUSMACHINE", p = 20.0, q = 0.0, pMin = 0.0, pMax = 40.0, qMin = -20.0, qMax = 20.0, vm_pu = 1.01, isRegulated = true, referencePriority = 2)
  addProsumer!(net = net, busName = "Eschwald_110", type = "SYNCHRONOUSMACHINE", p = 60.0, q = 0.0, pMin = 0.0, pMax = 150.0, qMin = -60.0, qMax = 60.0, vm_pu = 1.01, isRegulated = true)
  addProsumer!(net = net, busName = "Fennwiese_110", type = "SYNCHRONOUSMACHINE", p = 30.0, q = 0.0, pMin = 0.0, pMax = 60.0, qMin = -30.0, qMax = 30.0, vm_pu = 1.02, isRegulated = true, referencePriority = 1)
  addProsumer!(net = net, busName = "Gruenau_110", type = "ENERGYCONSUMER", p = 100.0, q = 25.0)
  refreshBusTypesFromProsumers!(net)
  return net
end

# bus index -> the bus name of the case file
_bus_names(net::Net) = Dict{Int,String}(idx => name for (name, idx) in net.busDict)

# one line per island: the bus it takes as reference and the reason the
# island report gives for it
function _reference_lines(net::Net)
  names = _bus_names(net)
  report = Sparlectra.detect_ac_islands(net; promote_generators = true)
  return [string("island ", row.island_id, ": reference ", get(names, row.chosen_ref_bus, "none"), " (", row.status, isempty(row.note) ? "" : string(", ", row.note), ")") for row in report.rows]
end

# the outage of the grid infeed at Almhof_110 (the north slack) through the
# N-1 engine; its result row names the bus that took the reference over
function _slack_outage_note(net::Net)
  infeed = findfirst(Sparlectra.isSlack, net.prosumpsVec)
  # generateN1Generators names a unit by its component name, with the
  # prosumer index appended where that name is not unique
  name = getCompName(net.prosumpsVec[infeed].comp)
  case = only(c for c in generateN1Generators(net) if c.element in (name, string(name, "#", infeed)))
  result = only(runContingencies!(net, [case]))
  return (converged = result.converged, note = something(result.error, ""))
end

"""
    main()

Load the shipped two-island case and show the reference choice:

1. with the stated priorities the south island takes `Fennwiese_110`
   (priority 1), not its largest unit;
2. with the priorities cleared it takes its strongest unit by
   `reference_candidate_rank`, `Eschwald_110` (before 0.30.2: the smallest
   PV bus index, `Dornberg_110`);
3. the island-wise power flow converges, each island on its own reference;
4. N-1 outage of the north slack: the north island hands the reference to
   `Brunntal_110` (priority 2), without priorities to `Kesselau_110` (the
   largest unit).
"""
function main()
  print_example_banner("examples/others/exp_island_reference_priority.jl", "reference priority: which bus each AC island takes as its reference, and why")
  case = joinpath(pkgdir(Sparlectra), "data", "scf", "two_islands_prio.scf.json")
  net = importSCF(case)
  names = _bus_names(net)
  println("stated reference priorities:")
  for ps in net.prosumpsVec
    ps.referencePriority > 0 && println("  ", names[Int(ps.comp.cFrom_bus)], ": ", ps.referencePriority)
  end
  with_priority = _reference_lines(net)

  plain = deepcopy(net)
  foreach(ps -> ps.referencePriority = 0, plain.prosumpsVec)
  without_priority = _reference_lines(plain)

  solved = deepcopy(net)
  iterations, erg = redirect_stdout(devnull) do
    runpf!(solved, 40, 1e-8, 0; islands_enabled = true)
  end
  # the active power each island's reference balances
  slack_p = [(bus, solved.nodeVec[net.busDict[bus]]._pƩGen) for bus in ("Almhof_110", "Fennwiese_110")]

  outage = redirect_stdout(devnull) do
    (with = _slack_outage_note(net), without = _slack_outage_note(plain))
  end
  return (with_priority = with_priority, without_priority = without_priority, iterations = iterations, converged = erg == 0, slack_p = slack_p, outage = outage)
end

result = run_example(main)
println()
println("with the stated priorities:")
foreach(line -> println("  ", line), result.with_priority)
println("with every priority cleared:")
foreach(line -> println("  ", line), result.without_priority)
println()
println("island-wise power flow: converged=", result.converged, " after ", result.iterations, " iteration(s)")
for (bus, p) in result.slack_p
  println("  reference ", bus, " balances its island with ", round(p; digits = 2), " MW")
end
println()
println("N-1 outage of the north grid infeed (Almhof_110):")
println("  with priorities:    converged=", result.outage.with.converged, ", ", result.outage.with.note)
println("  without priorities: converged=", result.outage.without.converged, ", ", result.outage.without.note)
