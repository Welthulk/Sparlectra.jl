```@meta
EditURL = "../../lit/workshop_tour_control.jl"
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

file: docs/lit/workshop_tour_control.jl
purpose: Literate.jl source of the workshop tour part 3, coordinated
         transformer control: why parallel transformers must not be
         regulated independently (circulating reactive power), the
         master/slave group (#322), and the CGMES mapping via shared
         RegulatingControl objects.

# The Sparlectra workshop tour, part 3: coordinated control

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_tour_control.ipynb)

> **Level: Advanced and up.** You should have the
> [basic tour](https://welthulk.github.io/Sparlectra.jl/generated/workshop_tour/)
> behind you (building nets, reading results, chapter 4's tap control);
> the [advanced tour](https://welthulk.github.io/Sparlectra.jl/generated/workshop_tour_advanced/)
> is helpful but not required.

Every real substation has them: two or three transformers in parallel
between the same busbars, regulated as ONE unit. This part of the
workshop is about why "as one unit" is not optional, what goes wrong
when each transformer brings its own controller, and how the
master/slave group models the coordinated regulation, down to the CGMES
data that describes it.

1. The substation, and a first innocent solve
2. The circulating-current trap: misaligned taps
3. Why two independent controllers make it worse
4. The master/slave group
5. Where the data comes from: CGMES RegulatingControl

## Warm-up and the substation

The study network is the classical two-transformer substation: a strong
110-kV grid connection, two identical 110/20-kV units in parallel onto
one 20-kV busbar, and a feeder load behind it.

```text
           HV (110 kV, slack 1.02 pu)
          /  \
        T1    T2      two identical units, ratio taps 0.9..1.1,
          \  /        step 0.0125
           LV (20 kV busbar)
           │
         Feeder (40 MW / 12 MVAr)
```

````@example workshop_tour_control
using Sparlectra

function build_substation(name::String)
  net = Net(name = name, baseMVA = 100.0)
  addBus!(net = net, busName = "HV", vn_kV = 110.0)
  addBus!(net = net, busName = "LV", vn_kV = 20.0)
  addBus!(net = net, busName = "Feeder", vn_kV = 20.0)
  addProsumer!(net = net, busName = "HV", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.02, va_deg = 0.0, referencePri = "HV")
  addProsumer!(net = net, busName = "Feeder", type = "ENERGYCONSUMER", p = 40.0, q = 12.0)
  addPIModelTrafo!(net = net, fromBus = "HV", toBus = "LV", r_pu = 0.006, x_pu = 0.12, b_pu = 0.0, ratio = 1.0, shift_deg = 0.0, status = 1)
  addPIModelTrafo!(net = net, fromBus = "HV", toBus = "LV", r_pu = 0.006, x_pu = 0.12, b_pu = 0.0, ratio = 1.0, shift_deg = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "LV", toBus = "Feeder", r_pu = 0.01, x_pu = 0.05, b_pu = 0.0, status = 1)
  # give both units their ratio-tap machinery
  for br in net.branchVec
    if br.ratio != 0.0
      br.has_ratio_tap = true
      br.tap_min = 0.9
      br.tap_max = 1.1
      br.tap_step = 0.0125
    end
  end
  ok, msg = validate!(net = net)
  ok || error("Network validation failed: $msg")
  return net
end

solve!(net) = begin
  ite, erg = runpf!(net, 30, 1e-8, 0)
  erg == 0 || error("power flow did not converge (status $erg)")
  calcNetLosses!(net)
  ite
end

net0 = build_substation("substation")
trafo_ids = [string(br.branchIdx) for br in net0.branchVec if br.ratio != 0.0]
solve!(net0)
t0 = [br for br in net0.branchVec if br.ratio != 0.0]
println("aligned taps: T1 carries ", round(t0[1].fBranchFlow.pFlow; digits = 2), " MW / ", round(t0[1].fBranchFlow.qFlow; digits = 2), " MVAr, ",
        "T2 carries ", round(t0[2].fBranchFlow.pFlow; digits = 2), " MW / ", round(t0[2].fBranchFlow.qFlow; digits = 2), " MVAr")
````

Reading aid: identical units, identical taps, so the load splits evenly.
This is the healthy state every coordination scheme wants to preserve.

## The circulating-current trap

Now the deliberate mistake: T1 steps four taps UP, T2 four taps DOWN.
The busbar voltage barely moves (the two errors cancel at the LV node),
but the RATIO MISMATCH drives a reactive current around the
HV-T1-LV-T2-HV loop that does nothing except heat both transformers:

````@example workshop_tour_control
net1 = build_substation("misaligned")
t1 = [br for br in net1.branchVec if br.ratio != 0.0]
t1[1].tap_ratio = 1.05; t1[1].ratio = 1.05
t1[2].tap_ratio = 0.95; t1[2].ratio = 0.95
solve!(net1)
println("misaligned taps: T1 carries ", round(t1[1].fBranchFlow.qFlow; digits = 2), " MVAr, T2 carries ", round(t1[2].fBranchFlow.qFlow; digits = 2), " MVAr")
println("LV busbar voltage: ", round(get_bus_vm_pu(net1, "LV"); digits = 4), " pu (vs ", round(get_bus_vm_pu(net0, "LV"); digits = 4), " pu aligned)")
````

Reading aid: the load still needs only 12 MVAr, but the loop now pushes
tens of MVAr through one unit and pulls them back through the other,
with opposite signs: that difference is pure circulation. Nothing in the
busbar voltage betrays it; you have to look at the per-unit flows. THIS
is why parallel transformers are regulated together: every scheme that
can leave the taps apart will eventually produce this picture.

## Why two independent controllers make it worse

The tempting setup, one voltage controller per transformer on the same
busbar target, is exactly such a scheme. Each controller measures the
same voltage but steps its own tap; nothing forces the two taps to move
together. Sparlectra warns at registration when a second independent
tap controller targets an already-regulated bus:

````@example workshop_tour_control
net2 = build_substation("fighting")
addPowerTransformerControl!(net2; trafo = trafo_ids[1], mode = :voltage, target_bus = "Feeder", target_vm_pu = 1.0)
# the second registration triggers the warning naming the trap
addPowerTransformerControl!(net2; trafo = trafo_ids[2], mode = :voltage, target_bus = "Feeder", target_vm_pu = 1.0)
````

Reading aid: the warning (look above) says it plainly: two controllers
on one voltage fight each other via circulating reactive power. On
well-behaved cases they may happen to move in lockstep; on any
asymmetry (unequal starting taps, one unit at a range end, measurement
ordering) they drift apart, and each tap step apart is one more notch
of circulation. The warning is a warning, not an error, because odd but
legitimate topologies exist; for parallel units the supported form is
the group.

## The master/slave group

One controller owns the group: the MASTER runs the ordinary discrete
voltage loop, and every accepted master step is mirrored onto the
FOLLOWERS step-synchronously (whole steps of each follower's own
`tap_step`, clamped to its range). Followers cannot carry their own
ratio controller and cannot follow two groups; both are enforced at
registration.

One subtlety deserves attention: the DEADBAND. Two synchronized units
move the busbar voltage twice as far per step as one unit, so the
group's deadband must cover at least half of that aggregated step
effect, or the loop cannot settle between two group steps:

````@example workshop_tour_control
net3 = build_substation("group")
addPowerTransformerControl!(net3; trafo = trafo_ids[1], followers = [trafo_ids[2]], mode = :voltage, target_bus = "Feeder", target_vm_pu = 1.03, deadband_vm_pu = 5e-3)
res = run_control!(net3)
b1 = Sparlectra._find_trafo_branch(net3, trafo_ids[1])
b2 = Sparlectra._find_trafo_branch(net3, trafo_ids[2])
calcNetLosses!(net3)
println("group: status = ", res.status, ", Vm(Feeder) = ", round(get_bus_vm_pu(net3, "Feeder"); digits = 4), " pu (target 1.03)")
println("  taps: master ", round(b1.tap_ratio; digits = 4), ", follower ", round(b2.tap_ratio; digits = 4), " (aligned)")
println("  reactive split: T1 ", round(b1.fBranchFlow.qFlow; digits = 2), " MVAr, T2 ", round(b2.fBranchFlow.qFlow; digits = 2), " MVAr (no circulation)")
for e in controllableElements(net3)
  println("  element view: ", e.element, " -> ", e.quantity, " @ ", e.target, ", status = ", e.status)
end
````

Reading aid: the master stepped, the follower mirrored, the taps stayed
aligned, and the reactive power splits evenly again: the healthy picture
from the first solve, now WITH the voltage target met. Declaratively the
same group is one `power_transformer` entry with a `followers` list
under `control.controllers`.

## Where the data comes from: CGMES RegulatingControl

CGMES models the group directly. Every tap changer may reference a
`RegulatingControl` object (`TapChanger.TapChangerControl`) carrying the
target (`mode`, `targetValue`, `targetDeadband`, the regulated
`Terminal`); several tap changers referencing ONE such object ARE the
regulated parallel group. There is no explicit master in the standard;
the per-changer `controlEnabled` flag expresses master/follower by tool
convention.

Sparlectra's importer (`importCGMES(...; tap_control = true)`) groups by
the shared object: the first enabled tap changer of a control becomes
the master, every further one joins as a follower, reported in the
import messages as `follows the group of ...`. Without that grouping, a
delivery with a regulated parallel pair would import as two independent,
fighting controllers, the exact trap from the previous chapter, produced
straight from the data with no user error involved. Tap changers with
`controlEnabled = false` stay at their fixed position.

## Where to go next

- [Control Framework](https://welthulk.github.io/Sparlectra.jl/control_framework/):
  the group's theory section, the hook interface, and the declarative
  controller configuration.
- [Transformer taps notebook](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_transformers.ipynb):
  the single-unit tap mechanics this part builds on.
- [Workshop tour, advanced](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_tour_advanced.ipynb):
  FACTS limits, N-1, state estimation, threads.
- [CGMES Import](https://welthulk.github.io/Sparlectra.jl/cgmes_import/):
  how deliveries map onto the network model.

