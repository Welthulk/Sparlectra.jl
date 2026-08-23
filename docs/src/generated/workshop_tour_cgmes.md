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
purpose: Literate.jl source of the CGMES tour: what an ENTSO-E delivery
         is, reading and analyzing it, bus-branch and node-breaker
         imports (with and without a TP profile, #314), SV validation,
         and the export round trip. Runs on the official conformity
         test sets, fetched on demand.

# The Sparlectra workshop tour: CGMES

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_tour_cgmes.ipynb)

> **Note:** This workshop was created with AI assistance and is reviewed
> and curated by the maintainer; it is not a fully machine-generated text.

> **Level: Expert.** You should be comfortable importing and solving
> networks (basic tour); no prior CGMES knowledge is assumed, building
> that knowledge is the point of this tour.

CGMES (Common Grid Model Exchange Standard) is how European TSOs
exchange grid models: a DELIVERY of several RDF/XML files, each carrying
one PROFILE of the same network. This tour works on the official
ENTSO-E conformity test sets, downloaded on demand (about 22 MB, once):

1. Anatomy of a delivery: profiles, and what a summary shows
2. When an import cannot work: the analysis report
3. Bus-branch import, solve, and validation against the shipped state
4. Node-breaker: with a TP profile, and without one (the topology processor)
5. The export round trip

## Warm-up and the test data

The ENTSO-E conformity package bundles reference networks in several
variants; `ensureCGMESTestConfigurations` downloads and extracts it once
into a local cache and returns the extraction root. Everything below
works on those files.

````@example workshop_tour_cgmes
using Sparlectra

root = Sparlectra.CGMESImporter.ensureCGMESTestConfigurations()
println("conformity sets under: ", root)

# the two study deliveries of this tour
microgrid_be = joinpath(root, "MicroGrid", "BaseCase_BC", "CGMES_v2.4.15_MicroGridTestConfiguration_BC_BE_v2")
microgrid_bd = joinpath(root, "MicroGrid", "BaseCase_BC", "CGMES_v2.4.15_MicroGridTestConfiguration_BD_v2")
minigrid_nb = joinpath(root, "MiniGrid", "NodeBreaker", "CGMES_v2.4.15_MiniGridTestConfiguration_BaseCase_Complete_v3")
minigrid_bd = joinpath(root, "MiniGrid", "NodeBreaker", "CGMES_v2.4.15_MiniGridTestConfiguration_Boundary_v3")
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
ite, erg = runpf!(net, 30, 1e-8, 0; islands_enabled = true)
println("solved: status ", erg, " in ", ite, " iterations")
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

## Where to go next

- [CGMES Import](https://welthulk.github.io/Sparlectra.jl/cgmes_import/):
  the mapping reference, config keys, placeholder guards, and the
  topology processor.
- [CGMES Export](https://welthulk.github.io/Sparlectra.jl/cgmes_export/):
  profiles written, identity stability, provenance.
- [Workshop tour, basic](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_tour.ipynb)
  and [advanced](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_tour_advanced.ipynb):
  everything you can do with the imported network.
- [Feature Matrix](https://welthulk.github.io/Sparlectra.jl/feature_matrix/)

