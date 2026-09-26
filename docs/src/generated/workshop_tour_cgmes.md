```@meta
EditURL = "../../lit/workshop_tour_cgmes.jl"
```

Copyright 2023-2026 Udo Schmitz

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

file: docs/lit/workshop_tour_cgmes.jl
purpose: Literate.jl source of the foreign-formats tour: the network
         formats of other tools that Sparlectra reads. CGMES (what an
         ENTSO-E delivery is, reading and analyzing it, bus-branch and
         node-breaker imports with and without a TP profile, #314, SV
         validation, the export round trip), PowSyBl IIDM (the same
         network as pypowsybl wrote it, compared with OpenLoadFlow) and
         the legacy DTF deck (FOR001 with its outage records). Runs on
         the official conformity test sets, fetched on demand, plus the
         shipped PowSyBl and DTF files.

# The Sparlectra workshop tour: foreign formats

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_tour_cgmes.ipynb)

> **Note:** This workshop was created with AI assistance and is reviewed
> and curated by the maintainer; it is not a fully machine-generated text.

> **Level: Expert.** You should be comfortable importing and solving
> networks (basic tour); no prior knowledge of the formats is assumed,
> building that knowledge is the point of this tour.

Sparlectra has formats of its own: the Sparlectra Case Format (`.scf.json`,
a power-grid-model dataset with a Sparlectra block), the plain
power-grid-model JSON inside it, and MATPOWER, the lingua franca of the
basic tour. Everything else is a FOREIGN format: the exchange format of
another tool, with that tool's model of a network, its conventions and
its blind spots. Three of them are read natively, and this tour takes
each one apart:

- **CGMES** (Common Grid Model Exchange Standard), how European TSOs
  exchange grid models: a DELIVERY of several RDF/XML files, each
  carrying one PROFILE of the same network. Chapters 1 to 5 work on the
  official ENTSO-E conformity test sets, downloaded on demand (about
  22 MB, once).
- **PowSyBl IIDM** (`.xiidm`), the XML network model of the LF Energy
  grid framework, with node-breaker substations, tap changers and
  limits; chapter 6 reads a file that pypowsybl wrote and compares the
  result with OpenLoadFlow.
- **DTF**, the fixed-column input deck (FOR001) of a legacy load-flow
  program, with its own outage records; chapter 7 reads one and applies
  an outage from the deck.

1. CGMES: anatomy of a delivery, profiles, and what a summary shows
2. CGMES: when an import cannot work, the analysis report
3. CGMES: bus-branch import, solve, and validation against the shipped state
4. CGMES: node-breaker, with a TP profile and without one (the topology processor)
5. CGMES: the export round trip
6. PowSyBl IIDM: the same network as pypowsybl delivers it, compared
   with OpenLoadFlow
7. DTF: a legacy input deck, its cards, and an outage from its records

## Warm-up and the test data

Three data sources, one per format. The ENTSO-E conformity package
bundles reference networks in several variants;
`ensureCGMESTestConfigurations` downloads and extracts it once into a
local cache and returns the extraction root (chapters 1 to 5). The
PowSyBl file and the DTF deck ship with Sparlectra under `data/powsybl`
and `data/DTF` (chapters 6 and 7), nothing to fetch.

````@example workshop_tour_cgmes
using Sparlectra

root = Sparlectra.CGMESImporter.ensureCGMESTestConfigurations()
println("conformity sets under: ", root)

# the two study deliveries of this tour
microgrid_be = joinpath(root, "MicroGrid", "BaseCase_BC", "CGMES_v2.4.15_MicroGridTestConfiguration_BC_BE_v2")
microgrid_bd = joinpath(root, "MicroGrid", "BaseCase_BC", "CGMES_v2.4.15_MicroGridTestConfiguration_BD_v2")
minigrid_nb = joinpath(root, "MiniGrid", "NodeBreaker", "CGMES_v2.4.15_MiniGridTestConfiguration_BaseCase_Complete_v3")
minigrid_bd = joinpath(root, "MiniGrid", "NodeBreaker", "CGMES_v2.4.15_MiniGridTestConfiguration_Boundary_v3")

# the shipped files of chapters 6 and 7
powsybl_dir = joinpath(pkgdir(Sparlectra), "data", "powsybl", "micro_grid_be.powsybl")
dtf_dir = joinpath(pkgdir(Sparlectra), "data", "DTF")
println("PowSyBl file: ", joinpath(powsybl_dir, "micro_grid_be.xiidm"))
println("DTF deck:     ", joinpath(dtf_dir, "FOR001.DAT"))
````

The first power-flow solve of a session compiles the solver (about a
minute); a two-bus warm-up net takes that hit here, so the timing of the
MicroGrid solve in Chapter 3 is the solve, not the compiler:

````@example workshop_tour_cgmes
wnet = Net(name = "warmup", baseMVA = 100.0)
addBus!(net = wnet, busName = "A", vn_kV = 110.0)
addBus!(net = wnet, busName = "B", vn_kV = 110.0)
addProsumer!(net = wnet, busName = "A", type = "EXTERNALNETWORKINJECTION", referencePri = "A", vm_pu = 1.0, va_deg = 0.0)
addProsumer!(net = wnet, busName = "B", type = "ENERGYCONSUMER", p = 10.0, q = 3.0)
addPIModelACLine!(net = wnet, fromBus = "A", toBus = "B", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
t_first = @elapsed runpf!(wnet, 10, 1e-8, 0; islands_enabled = true)
t_second = @elapsed runpf!(wnet, 10, 1e-8, 0; islands_enabled = true)
calcNetLosses!(wnet)
println("warm-up solve: first ", round(t_first; digits = 2), " s (compiles), second ", round(t_second * 1000; digits = 2), " ms")
````

## Chapter 1: anatomy of a delivery

**Example 1.1: what is in the box.** MicroGrid BE is the Belgian half of
a two-country reference model, as a sketch:

```text
    NL area (boundary) ~~~~ X-nodes ~~~~ BE area (this delivery)
                                          380 kV ring: 2 substations
                                          220/110 kV under 2 transformers
                                          gens, loads, shunt, PST
```

One delivery = one folder (or ZIP) of profile files. The important ones:

- **EQ** (equipment): what exists, lines, transformers, machines, their
  physical parameters. The skeleton.
- **SSH** (steady-state hypothesis): the operating point, set points,
  switch states, in-service flags.
- **TP** (topology): the sender's bus partition, `TopologicalNode`s and
  which terminal sits on which node.
- **SV** (state variables): the sender's solved state, voltages and
  flows. Reference material, not input.
- **Boundary set** (EQ_BD/TP_BD): the X-nodes where areas stitch
  together; shipped separately so neighbors agree on the seam.

`summarizeCGMES` reads a delivery without importing it:

````@example workshop_tour_cgmes
s = summarizeCGMES(path = [microgrid_be, microgrid_bd])
print(s)
````

Reading aid (Example 1.1): the summary names the files with their
profiles, the CGMES version, object counts per class, and the declared
model dependencies. Nothing is mapped yet; this is the "what did I just
receive" view.

## Chapter 2: when an import cannot work

**Example 2.1: the missing boundary set.** Deliveries reference objects
they do not define, most prominently the boundary X-nodes. Import the BE
files WITHOUT the boundary set and the import aborts; `analyzeCGMES`
explains the gap in plain language instead of a bare error:

````@example workshop_tour_cgmes
a_report = analyzeCGMES(path = microgrid_be)
print(a_report)
````

Reading aid (Example 2.1): the analysis lists the supplied models, the
declared dependencies (`md:Model.DependentOn`), which of them are
missing from the input, and a histogram of unresolved references. The
missing boundary set is named BY MODEL ID: exactly the file to ask the
sender for. The same report lands in `cgmes.log` when a Web UI or
service import aborts.

## Chapter 3: bus-branch import, solve, validate

**Example 3.1: the working import.** With the boundary set the BE model
(network of Example 1.1) imports into an ordinary `Net`; everything the
other workshops do applies from here:

````@example workshop_tour_cgmes
res = importCGMES(path = [microgrid_be, microgrid_bd], name = "microgrid_be")
net = res.net
println("imported: ", length(net.nodeVec), " buses, ", length(net.branchVec), " branches, ", length(net.prosumpsVec), " injections")
etime = @elapsed ((ite, erg) = runpf!(net, 30, 1e-8, 0; islands_enabled = true))
println("solved: status ", erg, " in ", ite, " iterations")
````

The classical result tables, exactly as in the basic tour:
`calcNetLosses!` derives branch flows and losses from the converged
voltages, `printACPFlowResults` prints bus voltages and branch flows.

````@example workshop_tour_cgmes
calcNetLosses!(net)
printACPFlowResults(net, etime, ite, 1e-8)
````

The delivery also carries the sender's solved state (SV), which makes it
its own reference: `compareWithSV` checks the re-solved voltages and
flows against what the sender shipped:

````@example workshop_tour_cgmes
cmp = compareWithSV(res)
println("compared ", cmp.n, " buses against the shipped SV: max |dVm| = ", round(cmp.max_dvm; sigdigits = 3), " pu, max |dVa| = ", round(cmp.max_dva; sigdigits = 3), " deg")
````

Reading aid (Example 3.1): agreement to a few 1e-4 pu means the model
mapped faithfully, the remaining difference is solver tolerance and
rounding in the shipped SV. A LARGE deviation would point at a mapping
or data problem; the per-bus table in the comparison names where.

## Chapter 4: node-breaker, with and without TP

**Example 4.1: node-breaker with the shipped TP.** Node-breaker
deliveries model every busbar, breaker, and disconnector explicitly;
the bus partition is the result of aggregating connectivity across
closed switches. Sketch of one MiniGrid substation:

```text
   busbar A ═══╤═══════╤═══   every ─[/]─ is a breaker or
              [/]     [/]     disconnector with its own
   busbar B ═══╧══╤════╧═══   SSH switch state; equipment
                 [/]          connects through bays, not
                line          directly to a busbar
```

As long as the sender ships a TP profile, the importer simply consumes
the sender's aggregation result:

````@example workshop_tour_cgmes
res_tp = importCGMES(path = [minigrid_nb, minigrid_bd], name = "minigrid_nb")
println("with TP: ", length(res_tp.net.nodeVec), " buses, ", length(res_tp.net.branchVec), " branches")
````

**Example 4.2: node-breaker WITHOUT a TP profile.** EMS and
substation-level exports often ship EQ+SSH only, topology expressed as
`ConnectivityNode`s plus switch states, and no topology processor ever
ran. Since 0.9.16 Sparlectra derives the partition itself: connectivity
nodes aggregate across closed switches (SSH `open` overriding EQ
`normalOpen`, out-of-service counts as open), retained switches stay
bus couplers, and the derived nodes feed the unchanged import pipeline.
We simulate such a delivery by withholding the TP and SV files of
Example 4.1's set:

````@example workshop_tour_cgmes
files_no_tp = [f for f in readdir(minigrid_nb; join = true) if endswith(f, ".xml") && !occursin("_TP", basename(f)) && !occursin("_SV", basename(f))]
res_notp = importCGMES(path = vcat(files_no_tp, minigrid_bd), name = "minigrid_nb_no_tp")
println("without TP: ", length(res_notp.net.nodeVec), " buses, ", length(res_notp.net.branchVec), " branches")
for m in res_notp.messages
  occursin("topology processor", m) && println("  ", m)
end

# the twins must describe the same network: solve both and compare
net_a = res_notp.net
net_b = res_tp.net
net_a.flatstart = true
net_b.flatstart = true
_, ea = runpf!(net_a, 40, 1e-8, 0; islands_enabled = true)
_, eb = runpf!(net_b, 40, 1e-8, 0; islands_enabled = true)
va = sort([n._vm_pu for n in net_a.nodeVec])
vb = sort([n._vm_pu for n in net_b.nodeVec])
println("both solve (", ea, "/", eb, "); max |Vm| difference across the sorted buses: ", maximum(abs.(va .- vb)))
````

Reading aid (Example 4.2): the derived partition reproduces the
sender's TP bus for bus (up to naming: derived buses take their busbar
names), and the solved voltages agree to numerical precision. The
processor announces itself in the import messages and runs ONLY when no
usable TP is present; TP-carrying deliveries take the unchanged path of
Example 4.1. One honest caveat from the conformity sweep: a TP that
assigns terminals of ONE connectivity node to DIFFERENT topological
nodes (FullGrid's completeness set does this on a load node) cannot be
derived from the connectivity graph; the processor then produces the
graph-consistent partition and the affected injections sit one bus
apart.

## Chapter 5: the export round trip

**Example 5.1: write it back out.** The imported network (Example 3.1)
exports as a complete CGMES 2.4.15 delivery (EQ+TP+SSH+SV) with
roundtrip-stable identities: mRIDs recorded on import are reused, so
renaming nothing changes nothing. The re-import of the export solves to
the same power flow:

````@example workshop_tour_cgmes
outdir = mktempdir()
files = writeCGMESFiles(net; path = outdir)
println(length(files), " profile files written")
res2 = importCGMES(path = outdir, name = "roundtrip")
ite2, erg2 = runpf!(res2.net, 30, 1e-8, 0; islands_enabled = true)
vm1 = sort([n._vm_pu for n in net.nodeVec])
vm2 = sort([n._vm_pu for n in res2.net.nodeVec])
println("re-import solves (status ", erg2, "); max |Vm| difference to the original: ", maximum(abs.(vm1 .- vm2)))
````

Reading aid (Example 5.1): the export writes what the network IS, not
what the original files said; the shipped SV of the export carries the
CURRENT solved state, so a receiving tool starts from it. Details and
the identity rules: [CGMES Export](https://welthulk.github.io/Sparlectra.jl/cgmes_export/).

## Chapter 6: PowSyBl IIDM, the same network compared with OpenLoadFlow

**What IIDM is.** PowSyBl, the LF Energy grid framework, imports CGMES
too and keeps its networks in its own XML model, IIDM (`.xiidm`). Where
CGMES spreads one network over profiles and RDF references, IIDM is one
document: substations with voltage levels, each level either a
bus-breaker topology (named buses) or a node-breaker topology (numbered
nodes joined by switches and internal connections), and the equipment
attached to buses or nodes; tap changers with their step tables,
operational limits, HVDC links and dangling lines sit right on the
elements. What pypowsybl resolves before it hands out tables (the bus
views of a node-breaker substation, the current tap step, the reactive
limits at the setpoint), Sparlectra's reader resolves in Julia, so the
file needs no Python and no conversion.

**Example 6.1: an IIDM file.** The MicroGrid BE of this tour ships with
Sparlectra in that form under `data/powsybl`, as pypowsybl wrote it.
Next to the file lies the solution OpenLoadFlow (PowSyBl's load flow)
computed for it, `reference_buses.csv`, so the import can be checked
against an independent solver:

````@example workshop_tour_cgmes
using DelimitedFiles
tables = Sparlectra.read_iidm_tables(joinpath(powsybl_dir, "micro_grid_be.xiidm"))
println("IIDM ", tables.manifest["iidm_version"], ": ", length(tables.buses.id), " buses, ", length(tables.two_windings_transformers.id), " two-winding and ", length(tables.three_windings_transformers.id), " three-winding transformer(s), ", length(tables.dangling_lines.id), " dangling lines")
````

The generators of the MicroGrid regulate remote buses; OpenLoadFlow holds
those buses at their targets, and `remote_regulation = :remote` attaches
the same as outer-loop machine controls. The import report lists what
was built, the slack decision and every notice:

````@example workshop_tour_cgmes
pnet, preport = build_net_from_powsybl(tables, PowsyblAdapterOptions(remote_regulation = :remote))
println(format_powsybl_report(preport))
````

Reading aid (Example 6.1): the notice about the file's bus state is the
reader being careful: the file was solved at other tap positions
(`solvedTapPosition`) than the ones it carries as the current position,
so that state belongs to another network and the import starts flat.

**Example 6.2: the comparison.** OpenLoadFlow's defaults are: every
synchronous component solved, the active-power mismatch shared over the
generators in proportion to `max_p` (the importer carries that as
participation factors), and reactive limits switched without
hysteresis. The matching Sparlectra configuration, then the bus voltages
against the reference file:

````@example workshop_tour_cgmes
olf_like = SparlectraConfig(
  powerflow = PowerFlowConfig(max_iter = 60, tol = 1e-9, islands_enabled = true, qlimits = Sparlectra.QLimitConfig(hysteresis_pu = 1e-6), distributed_slack = Sparlectra.DistributedSlackConfig(enabled = true, p_mode = :imported)),
  output = OutputConfig(logfile_results = :off, startup_latency_hint = false),
  control = ControlConfig(),
)
pres = run_sparlectra(net = pnet, config = olf_like)
@assert pres.numerical_converged
ref, hdr = readdlm(joinpath(powsybl_dir, "reference_buses.csv"), ',', String; quotes = true, header = true)
cols = Dict(Symbol(strip(h)) => k for (k, h) in enumerate(vec(hdr)))
worst = 0.0
for i in axes(ref, 1)
  bus = String(ref[i, cols[:bus_breaker_id]])
  haskey(pnet.busDict, bus) || continue
  dv = abs(pnet.nodeVec[pnet.busDict[bus]]._vm_pu - parse(Float64, ref[i, cols[:v_pu]]))
  global worst = max(worst, dv)
end
println("worst |dV| against OpenLoadFlow over ", size(ref, 1), " buses: ", round(worst; sigdigits = 3), " pu")
@assert worst < 1e-5
````

Reading aid (Example 6.2): two load-flow programs with different
conventions (ratio at the from side against `rho` at side 1, the
magnetizing admittance placement, the slack distribution) agree to
1e-9 pu once the conventions are mapped, which is what the PowSyBl
import settled. The same comparison runs for any `.xiidm` of your own
next to a reference written with pypowsybl; without one, the file alone
imports just the same. Details: [PowSyBl Import](https://welthulk.github.io/Sparlectra.jl/powsybl_import/).

## Chapter 7: DTF, a legacy input deck with its own outage records

**What DTF is.** A fixed-column text deck of a legacy load-flow program,
the format of the Testnetz13 validation examples: FOR001 is the input
(the network plus run parameters), FOR002 the program's printed result
report. A deck is a sequence of CARDS in a fixed order: parameter and
text cards, the nominal voltages the voltage-level indices refer to, a
size card with the bus and branch counts and the NAMED slack bus, the
branch cards (`L` for a line, `T` for a transformer, with impedances in
per unit of the level), compensation cards, transformer-control cards
(winding voltages, the longitudinal tap range and step, an optional
skew-angle regulator), the bus cards (type, level index, name, start
voltage, load, generation), and after the buses an optional block of
OUTAGE records between `AUSFALL` and `ENDE`, one branch outage per line.
There is no per-terminal switch: a branch is in or out. Sparlectra reads
the deck into typed records, raw lines included, and builds the network
from them.

**Example 7.1: reading the deck.** `read_dtf` parses the cards;
everything the deck says stays on the case object before any network
exists:

````@example workshop_tour_cgmes
case = Sparlectra.DTFImporter.read_dtf(joinpath(dtf_dir, "FOR001.DAT"); strict = false)
println("deck: base ", case.baseMVA, " MVA, ", length(case.buses), " buses, ", length(case.branches), " branches (", count(b -> b.kind == 'T', case.branches), " transformers), ", length(case.transformer_controls), " transformer control(s), ", length(case.outages), " outage record(s)")
println("nominal voltages (kV): ", join(case.nominal_voltages_kv, ", "), "; slack bus of the size card: ", case.size.slack)
@assert length(case.buses) > 0 && length(case.branches) > 0
````

**Example 7.2: the network and its solve.** `build_net` turns the records
into a `Net` (transformer transverse admittance on the branch, the named
slack as the reference, PQ buses kept as fixed injections); from there
it is the ordinary solver:

````@example workshop_tour_cgmes
dnet = Sparlectra.DTFImporter.build_net(case)
println("network: ", length(dnet.nodeVec), " buses, ", length(dnet.branchVec), " branches")
@assert length(dnet.nodeVec) == length(case.buses)
dite, derg = runpf!(dnet, 50, 1e-8, 0)
@assert derg == 0
calcNetLosses!(dnet)
println("solved in ", dite, " iterations; losses ", round(dnet.totalLosses[end][1]; digits = 3), " MW")
````

Reading aid (Example 7.2): the FOR002 report that ships next to the deck
is the legacy program's result for the same deck; the example
`examples/dtf/dtf_validation_base.jl` parses it and compares bus
voltages, branch flows and generator reactive power line by line. That
comparison is how the DTF path was validated, and the place to look when
a deck of your own disagrees with its old report.

**Example 7.3: an outage from the deck.** The records between `AUSFALL`
and `ENDE` name branches by kind, level, parallel identifier and the two
bus names; `find_outage_branch_indices` resolves one to the base
network's branches (exactly one match is required, a miss or an
ambiguity is a diagnostic, not a guess), `apply_single_branch_outage!`
takes the branch out, and the solve repeats:

````@example workshop_tour_cgmes
if isempty(case.outages)
  println("this deck carries no outage records")
else
  outage = first(case.outages)
  matches = Sparlectra.DTFImporter.find_outage_branch_indices(case, outage)
  println("first outage record resolves to branch index(es) ", matches, ": ", Sparlectra.DTFImporter.outage_match_diagnostic(case, outage, matches))
  if length(matches) == 1
    onet = Sparlectra.DTFImporter.build_net(case)
    Sparlectra.DTFImporter.apply_single_branch_outage!(onet, matches[1])
    markIsolatedBuses!(net = onet, log = false)
    oite, oerg = runpf!(onet, 50, 1e-8, 0; islands_enabled = true)
    @assert oerg == 0
    dv = maximum(abs(onet.nodeVec[i]._vm_pu - dnet.nodeVec[i]._vm_pu) for i in eachindex(dnet.nodeVec))
    println("with the outage: solved in ", oite, " iterations, largest |Vm| change ", round(dv; sigdigits = 3), " pu")
  end
end
````

Reading aid (Example 7.3): the outage records are part of the deck, so a
DTF study is reproducible from the file alone; the Web UI and the service
apply them the same way (`for001Contingencies`). What the deck cannot
say (a one-sided open branch, a switch) has no record, and the importer
does not invent one.

## Where to go next

- [PowSyBl Import](https://welthulk.github.io/Sparlectra.jl/powsybl_import/):
  what the IIDM reader resolves, the conventions against OpenLoadFlow.
- [DTF Format](https://welthulk.github.io/Sparlectra.jl/dtf_format/):
  the cards, the outage records, the validation workflow against FOR002.
- [CGMES Import](https://welthulk.github.io/Sparlectra.jl/cgmes_import/):
  the mapping reference, config keys, placeholder guards, and the
  topology processor.
- [CGMES Export](https://welthulk.github.io/Sparlectra.jl/cgmes_export/):
  profiles written, identity stability, provenance.
- [Workshop tour, basic](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_tour.ipynb)
  and [advanced](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_tour_advanced.ipynb):
  everything you can do with the imported network.
- [Feature Matrix](https://welthulk.github.io/Sparlectra.jl/feature_matrix/)

