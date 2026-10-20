# [Dynamic Unit Data](@id dynamic_data)

A generating unit can carry dynamic parameters: inertia, damping, primary
and secondary control, unit type. They sit in an optional
[`DynamicData`](@ref) record on the unit (`ProSumer.dynamics`, default
`nothing`).

**No calculation reads these fields.** Power flow, state estimation, short
circuit, N-1 and QSTS give the same results with and without them. They are
data for a later frequency analysis and for the exchange through SCF. The
only reader is the island slack choice of the CGMES import, which uses the
unit type (see [below](#Unit-type-and-island-slack)).

## Fields

| Field | Unit | Meaning | Default |
|---|---|---|---|
| `inertia_s` | s | inertia constant H on the rating of the unit | `nothing` (no inertia) |
| `damping_pu` | pu/pu | damping coefficient D | 0.0 |
| `rating_mva` | MVA | rating that H and D refer to, when it differs from the unit's `ratedS` | `nothing` (`ratedS` applies) |
| `droop_pu` | pu | primary control droop (statics), for example 0.05 | `nothing` (no primary control) |
| `deadband_hz` | Hz | primary control deadband | 0.0 |
| `primary_reserve_mw` | MW | primary reserve limit, symmetric | 0.0 |
| `governor_time_s` | s | governor/turbine time constant, first order | `nothing` |
| `agc_participation` | pu | share of the secondary control signal, 0 to 1 | 0.0 |
| `secondary_reserve_mw` | MW | secondary reserve limit, symmetric | 0.0 |
| `unit_type` | enum | `SynchronousUnit`, `GridFollowingUnit`, `GridFormingUnit` | `SynchronousUnit` when `inertia_s` is set, else `GridFollowingUnit` |
| `virtual_inertia_s` | s | synthetic inertia of a grid-forming unit | `nothing` |
| `black_start` | | the unit can energise a dead network | `false` |
| `min_stable_mw` | MW | minimum stable generation | `nothing` |

`inertia_s` and `virtual_inertia_s` are separate, so a converter can carry
a rating and a synthetic inertia without counting as a machine.

The constructor checks the values: every number finite, inertia, virtual
inertia, time constant, deadband and reserves not negative, droop and
rating greater than zero when set, `agc_participation` in [0, 1]. A
violation throws [`DynamicDataError`](@ref), which names the field.

```julia
using Sparlectra
net = importSCF("data/scf/sp_case5.scf.json")
gen = first(ps for ps in net.prosumpsVec if isGenerator(ps))
gen.dynamics = DynamicData(inertia_s = 5.0, droop_pu = 0.05, primary_reserve_mw = 20.0)
```

## Sources

| Source | What it fills |
|---|---|
| SCF | every field, block `extra.<machine>.dynamics` ([Sparlectra Case Format](@ref scf_dynamics)) |
| CGMES | inertia, damping, rating and unit type from the DY profile, droop and governor time constant of some turbine governors, with `cgmes_import.read_dynamics = true` ([CGMES Import](@ref cgmes_dynamics)) |
| MATPOWER, PowSyBl IIDM, PGM, DTF | nothing; the field stays `nothing` |

## Unit type and island slack

The CGMES import gives every electrical island that has infeed but no
reference a slack. A unit whose `unit_type` is `GridFollowingUnit` is no
candidate for it: it cannot hold voltage and frequency of an island. An
island whose units are all grid following (static var compensators and
loads aside) cannot stand in operation; the import switches it off with a
warning, like a part without source. Synchronous and grid-forming units
remain candidates. A unit without `DynamicData` is a candidate as before;
for a `PowerElectronicsConnection` without it the warning says that its
type is unknown and assumed grid following.
