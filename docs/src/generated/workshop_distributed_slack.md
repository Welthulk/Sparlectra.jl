```@meta
EditURL = "../../lit/workshop_distributed_slack.jl"
```

# Distributed slack: who covers the imbalance?

> **Level: Advanced**, deepens the slack chapters of the basic tour.

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_distributed_slack.ipynb)

> **Note:** This workshop was created with AI assistance and is reviewed
> and curated by the maintainer; it is not a fully machine-generated text.

A classical power flow gives the entire active-power imbalance, every
megawatt of load not covered by scheduled generation, plus all network
losses, to ONE reference bus. Real interconnections do not work that way:
primary and secondary control spread the imbalance over many governors,
each contributing according to its participation factor. The distributed
slack in [Sparlectra.jl](https://github.com/Welthulk/Sparlectra.jl) models
exactly that: the imbalance $\lambda$ becomes an extra solver state, and
every participating generator picks up its share
$\Delta P_i = \alpha_i \, \lambda$ with $\sum_i \alpha_i = 1$.

This notebook walks through the one question everything else hangs on:
WHERE do the weights $\alpha_i$ come from? You will run the same shortfall
case under every weight mode, watch the shares change, and see the honest
failure when no valid participant exists.

> **Note:** On Google Colab the install cell takes a few minutes on a
> fresh session (package download and precompilation). Colab's Julia
> version may change over time; this notebook targets Julia ≥ 1.12.

## Warm-up and shared helpers

Julia compiles each function on first use. This cell loads the package,
collects the helpers of the whole notebook, and warms both solver paths
(classical single slack and distributed slack) on the study case itself,
so the chapters below run at full speed.

````@example workshop_distributed_slack
using Sparlectra
````

The study case as a MATPOWER file: a four-bus ring, drawn as a diagram
in the study-case section below, with a DELIBERATE
20 MW shortfall (70 MW load, only 50 MW scheduled PV generation), so
somebody must visibly cover the gap. The 21st generator column is the
MATPOWER participation factor APF: 0.75 for G2 and 0.25 for G3,
deliberately DIFFERENT from their 30:20 schedule ratio, so the weight
modes below give visibly different shares.

````@example workshop_distributed_slack
function write_demo_case(dir)
  path = joinpath(dir, "case4_distributed_slack.m")
  write(path, """
  function mpc = case4_distributed_slack
  mpc.version = '2';
  mpc.baseMVA = 100;
  mpc.bus = [
  1 3 0 0 0 0 1 1.00 0 110 1 1.1 0.9;
  2 2 0 0 0 0 1 1.01 0 110 1 1.1 0.9;
  3 2 0 0 0 0 1 1.00 0 110 1 1.1 0.9;
  4 1 70 20 0 0 1 1.00 0 110 1 1.1 0.9;
  ];
  mpc.gen = [
  1 0 0 300 -300 1.00 100 1 250 -250 0 0 0 0 0 0 0 0 0 0 0;
  2 30 0 50 -50 1.01 100 1 100 0 0 0 0 0 0 0 0 0 0 0 0.75;
  3 20 0 50 -50 1.00 100 1 100 0 0 0 0 0 0 0 0 0 0 0 0.25;
  ];
  mpc.branch = [
  1 2 0.01 0.08 0.02 0 0 0 0 0 1 -360 360;
  2 3 0.02 0.10 0.02 0 0 0 0 0 1 -360 360;
  3 4 0.01 0.06 0.02 0 0 0 0 0 1 -360 360;
  1 4 0.03 0.12 0.02 0 0 0 0 0 1 -360 360;
  ];
  """)
  return path
end
case_path = write_demo_case(mktempdir())

# fresh import for every run (each run mutates the net in place)
load_case() = createNetFromMatPowerFile(filename = case_path)

# the interesting quantity per run: which bus delivers how much active
# power BEYOND its schedule (uses two solver internals, fine for a demo)
function print_beyond_schedule(net)
  Y = Sparlectra.createYBUS(net = net, sparse = true, printYBUS = false)
  V = [n._vm_pu * cis(deg2rad(n._va_deg)) for n in net.nodeVec]
  extra = (real.(Sparlectra.calc_injections(Y, V)) .- real.(Sparlectra.buildComplexSVec(net))) .* net.baseMVA
  for (i, n) in enumerate(net.nodeVec)
    abs(extra[i]) < 0.005 && continue
    println("    bus ", i, " (", Sparlectra.getCompName(n.comp), ")  +", round(extra[i]; digits = 2), " MW beyond schedule")
  end
end

# participation table from the solver status (bus, alpha, delta P)
function print_participation(net)
  st = Sparlectra.rectangular_pf_status(net)
  println("    lambda = ", round(st.distributed_slack_lambda_p_mw; digits = 2), " MW over ", st.distributed_slack_participants, " participant bus(es), mode ", st.distributed_slack_mode)
  for row in st.distributed_slack_participation
    println("    ", row.bus, ": alpha = ", round(row.alpha; digits = 3), "  ->  dP = ", round(row.dp_mw; digits = 2), " MW (scheduled ", round(row.pg_mw; digits = 1), " MW)")
  end
end

t_first = @elapsed runpf!(load_case(), 30, 1e-8, 0)
t_dslack = @elapsed runpf!(load_case(), 30, 1e-8, 0; distributed_slack_enabled = true)
println("warm: classical ", round(t_first; digits = 2), " s, distributed ", round(t_dslack; digits = 2), " s (first calls compile)")
````

## The study case

```text
  G1 (slack)              G2 (PV, Pg 30 MW, APF 0.75)
   |                       |
   B1 ------- line ------- B2
   |                       |
  line                    line
   |                       |
   B4 ------- line ------- B3
   |                       |
 load 70 MW               G3 (PV, Pg 20 MW, APF 0.25)
```

Scheduled PV generation covers 50 of the 70 MW; the remaining 20 MW plus
the network losses are the imbalance $\lambda$ someone has to deliver.

## Chapter 1: the classical single slack

**Example 1: the classical single slack.** Without distributed slack,
the answer is brutal: the REF bus alone
absorbs everything. The PV buses deliver exactly their schedule, to the
kilowatt, no matter how large the gap is.

````@example workshop_distributed_slack
net = load_case()
runpf!(net, 30, 1e-8, 0)
println("classical single slack:")
print_beyond_schedule(net)
````

Reading aid (Example 1): G1 delivers the 20 MW shortfall plus all
losses. On a small
demo that is harmless; on a real interconnection it concentrates hundreds
of MW on one arbitrary bus, distorts the flows around it, and tells you
nothing about how the fleet would actually respond.

## Chapter 2: pg_weighted, share by schedule (the default)

**Example 2: share by schedule.**
`distributed_slack_enabled = true` turns the imbalance into a solver
state. The weights come from `p_mode`. The default `:pg_weighted` uses
the SCHEDULED Pg of every generator-type prosumer at a REF or PV bus:
G2 with 30 MW and G3 with 20 MW get raw weights 30 and 20, normalized to
$\alpha = 0.6 / 0.4$. Big units take proportionally more, the classical
"participation by size".

````@example workshop_distributed_slack
net = load_case()
runpf!(net, 30, 1e-8, 0; distributed_slack_enabled = true, distributed_slack_p_mode = :pg_weighted)
println("pg_weighted (raw weights 30 : 20):")
print_participation(net)
print_beyond_schedule(net)
````

## Chapter 3: imported, share by declared participation factors

**Example 3: share by declared participation factors.**
Grid data often DECLARES the participation: MATPOWER carries it in the
optional `APF` generator column, CGMES in `normalPF`. Sparlectra imports
both into `ProSumer.participationFactor`, and `p_mode = :imported` uses
them as raw weights. Our case declares 0.75 / 0.25, deliberately NOT the
30:20 schedule ratio, and the shares follow the declaration, not the
schedule:

````@example workshop_distributed_slack
net = load_case()
runpf!(net, 30, 1e-8, 0; distributed_slack_enabled = true, distributed_slack_p_mode = :imported)
println("imported (MATPOWER APF 0.75 : 0.25):")
print_participation(net)
print_beyond_schedule(net)
````

## Chapter 4: headroom_weighted, share by free capacity

**Example 4: share by free capacity.**
`:headroom_weighted` uses `max(maxP - Pg, 0)`: whoever has the most
UNUSED capacity takes the most. Two things change against the earlier
modes (Examples 2 and 3). First, the REFERENCE generator joins the
pool: candidates are the generators at REF or PV buses, and under
pg_weighted (Example 2, Pg = 0) and imported (Example 3, no APF) the
slack unit was dropped as invalid, but its
headroom of 250 MW is a perfectly valid weight. The raw weights are
therefore 250 : 70 : 80. Second, note the reversal among the PV units:
G3, the SMALLER unit, now takes more than G2, because 100 - 20 leaves
more reserve than 100 - 30. There is also `:pmax_weighted` (share by
installed size) for markets that contract by capacity.

````@example workshop_distributed_slack
net = load_case()
runpf!(net, 30, 1e-8, 0; distributed_slack_enabled = true, distributed_slack_p_mode = :headroom_weighted)
println("headroom_weighted (raw weights 250 : 70 : 80):")
print_participation(net)
print_beyond_schedule(net)
````

## Chapter 5: explicit, you decide

**Example 5: explicit weights.**
`p_mode = :explicit` takes the weights from a table you pass in (config
key `power_flow.distributed_slack.weights`, or the keyword below), keyed
by bus name or bus index. Weights need not sum to one; normalization is
automatic. Here bus 2 gets everything:

````@example workshop_distributed_slack
net = load_case()
runpf!(net, 30, 1e-8, 0; distributed_slack_enabled = true, distributed_slack_p_mode = :explicit, distributed_slack_weights = Dict("2" => 1.0))
println("explicit (all on bus 2):")
print_participation(net)
print_beyond_schedule(net)
````

## The rules behind the weights, in one place

Whatever the mode, the same pipeline runs per synchronous island:

- Candidates are the generator-type prosumers at the island's REF or PV
  buses. Fixed PQ injections (HVDC converter stand-ins, boundary
  equivalents) never participate; the bus-type gate excludes them.
- Invalid candidates are dropped, not guessed: a weight that is missing,
  zero, negative, NaN, or Inf removes the unit from the pool (with a
  debug log and a drop counter in the result metadata).
- The surviving raw weights are normalized to $\sum \alpha = 1$ and
  aggregated per bus; the solver then adds one state $\lambda$ and every
  participant moves by $\alpha_i \lambda$. With
  `distributed_slack_respect_p_limits = true` (default) a participant is
  clamped at its P limit and the remainder redistributes over the rest.

## The honest failure: no valid participant

**Example 6: the honest failure.** What if the declared data is missing?
`:imported` (Example 3) on a case WITHOUT any
APF column drops every candidate. The default
`distributed_slack_fallback = :error` refuses to guess and throws an
actionable error; `fallback = :ref_only` warns and solves classically
instead, so batch runs survive a case with missing declarations:

````@example workshop_distributed_slack
no_apf = load_case()
for ps in no_apf.prosumpsVec        ## strip the imported factors
  ps.participationFactor = nothing
end
runpf!(no_apf, 30, 1e-8, 0; distributed_slack_enabled = true, distributed_slack_p_mode = :imported, distributed_slack_fallback = :ref_only)
println("fallback ref_only: solved classically, the slack absorbs everything again")
print_beyond_schedule(no_apf)
````

Reading aid (Example 6): with `:error` (the default) the message names the island,
the mode, and the number of dropped candidates, and points at
`fallback=ref_only`; nothing silently degrades to a single slack.

One configuration note for `run_sparlectra`/YAML users: the same switches
live under `power_flow.distributed_slack.*` (`enabled`, `p_mode`,
`fallback`, `weights`, `respect_p_limits`). Distributed slack and the
non-ideal external-grid source are mutually exclusive by validation:
both answer "who covers the imbalance", so you enable at most one.

## Where to go next

- [Slack types and short circuit](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_slack_short_circuit.ipynb):
  ideal slack vs. external-grid source vs. distributed slack on one
  8-bus network.
- [Solver documentation](https://welthulk.github.io/Sparlectra.jl/solver/):
  the rectangular Newton-Raphson formulation the extra state plugs into.
- [Workshop tour](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_tour.ipynb):
  all workshop examples in one Colab session.

