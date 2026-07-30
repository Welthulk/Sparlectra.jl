# Short-Circuit Analysis

!!! note "Status"
    Sparlectra currently **reads and reports** short-circuit source data from
    CGMES deliveries; the calculation itself is not yet available. This page
    explains what short-circuit analysis is, how it is commonly done, and how
    Sparlectra approaches it.

## What short-circuit analysis is

A power-flow solution describes the *normal* operating state of a network.
Short-circuit analysis answers a different question: **what current flows
when the network faults** — and therefore what the switchgear must be able
to interrupt, what busbars and conductors must withstand mechanically and
thermally, and whether protection relays will detect the fault and trip
selectively. Both extremes matter: the *maximum* fault current sizes the
equipment; the *minimum* fault current is the one a protection scheme must
still be able to see.

## How it is commonly done: IEC 60909

The industry-standard method is IEC 60909. It is deliberately *not* a
time-domain simulation: it is a quasi-static calculation with defined
safety margins, so results are reproducible and conservative without
detailed dynamic machine models.

**Equivalent voltage source at the fault location.** The entire network is
reduced to the short-circuit impedance seen from the faulted bus, driven by
a single equivalent source. All other sources are removed; loads and line
charging are neglected. The initial symmetrical short-circuit current at
nominal voltage `Un` follows from that impedance `Zk`:

```math
I_k'' = \frac{c \, U_n}{\sqrt{3}\,\lvert Z_k \rvert}
```

**The voltage factor** `c` absorbs the difference between nominal and
actual pre-fault voltage plus a safety margin. It is tabulated by voltage
level and by objective: `c_max` (for example 1.10 in HV networks) for
maximum currents, `c_min` (0.95 or 1.00) for minimum currents.

**Derived quantities.** Everything else follows from the initial current:
the **peak current** (first crest, mechanical stress)

```math
i_p = \kappa \,\sqrt{2}\; I_k'',
\qquad
\kappa \approx 1.02 + 0.98\, e^{-3R/X}
```

— which is why the resistance-to-reactance ratio at the fault matters, not
just the impedance magnitude — plus the **breaking current** at contact
separation and the **thermal equivalent current** for the fault duration.

**Where the fault current comes from.** A synchronous machine feeds the
fault through its subtransient reactance (`x''_d`); a network feeder is
represented by an equivalent impedance derived from its declared maximum
short-circuit power and its R/X ratio; asynchronous motors contribute
briefly and matter mostly for the peak current and for faults close to the
motor — the standard defines when they may be neglected. Neglecting motor
contributions makes the *maximum* current an underestimate, which is the
non-conservative direction; a careful implementation must say so when it
does it.

**Unbalanced faults.** Single line-to-earth — the most common real fault —
line-to-line and two-phase-to-earth faults are computed with **symmetrical
components**: the fault current follows from combinations of the positive-,
negative- and zero-sequence networks. The zero-sequence network is where
transformers dominate: their vector group and star-point earthing decide
whether a zero-sequence current path exists at all.

## How the network is solved: Z-bus, and Takahashi

For one fault location the needed quantity is a single diagonal element of
the bus impedance matrix: with the short-circuit admittance matrix `Y_sc`
(series impedances only, source admittances on the diagonal, no loads or
shunts), the fault impedance at bus `i` is

```math
Z_{ii} = \left( Y_{sc}^{-1} \right)_{ii}
```

Nobody inverts the full matrix for this. The practical method is the
**Z-bus column solve**: factorize `Y_sc` once (sparse LU), then obtain
column `i` of the inverse with one forward/backward substitution against
the `i`-th unit vector. That column contains `Z_ii` for the fault current
*and* the off-diagonal entries needed for the voltage sag and the branch
contributions around the fault.

**What the Takahashi method adds.** For a *sweep over all buses* the column
approach costs one triangular solve per bus. Takahashi, Fagan and Chen
(1973) showed that all inverse elements matching the sparsity pattern of
the LU factors — in particular the **entire diagonal** — can be computed in
a single backward pass over the factors, at a cost comparable to the
factorization itself. This *sparse inverse* (also known as the Takahashi
equations) is the classic power-system technique for exactly this job:
every `Z_ii` of a large network in one pass, without ever forming the
dense inverse. A full-network fault-level study therefore scales like one
factorization, not like `n` solves.

Sparlectra's approach follows this pattern: per-bus faults via the column
solve; the Takahashi sparse inverse is the designated optimization for
all-bus fault-level sweeps on large networks.

## How Sparlectra does it

**Available today: the data, and its completeness.** Every CGMES import
harvests the short-circuit source data into
`CGMESImportResult.shortcircuit` — read, never silently altered, values in
CGMES units:

| Class | Harvested attributes | Role in the calculation |
|---|---|---|
| `ExternalNetworkInjection` | `maxInitialSymShCCurrent`, `minInitialSymShCCurrent`, `maxR1ToX1Ratio`, `minR1ToX1Ratio`, `maxR0ToX0Ratio`, `maxZ0ToZ1Ratio`, `ikSecond`, `governorSCD` | network-feeder equivalent impedance in the positive and zero sequence |
| `SynchronousMachine` | `satDirectSubtransX` (`x''_d`), `satDirectTransX`, `r0`, `x0`, `r2`, `x2`, `earthing`, `ratedS`, `ratedU` | machine contribution on machine base, sequence impedances, earthing |
| `ACLineSegment` | `r0`, `x0`, `b0ch`, `g0ch`, `shortCircuitEndTemperature` | zero-sequence line model; conductor temperature for minimum-current cases |
| `PowerTransformerEnd` | `r0`, `x0`, `grounded`, `rground`, `xground` | zero-sequence transformer model and star-point treatment |
| `EquivalentInjection` | `r`, `x`, `r0`, `x0`, `r2`, `x2` | boundary equivalents in all three sequence networks |

Whether a given delivery could support a calculation is measurable:
`shortCircuitCoverage(result.shortcircuit)` reports, per class, the record
count and per attribute the fill rate; `printShortCircuitCoverage` renders
it, and every CGMES run writes the same view into `cgmes.log`.

**The calculation, as designed.** The first stage is the balanced
three-phase fault: maximum and minimum initial symmetrical current per
selected bus (an all-bus sweep is a loop over the same path), IEC voltage
factors by level, peak current from the R/X ratio at the fault. Machine
data that is incomplete is substituted with documented defaults or skipped
with a warning — and because a short-circuit result is safety-relevant, any
such substitution is **flagged on the affected result rows themselves**,
not only logged. Unbalanced faults and the derived breaking/thermal
quantities follow once the zero-sequence model (transformer vector groups
and earthing, harvested above) is assembled.

## References

- IEC 60909-0, *Short-circuit currents in three-phase a.c. systems — Part 0:
  Calculation of currents*; IEC 60909-4 for worked verification examples.
- K. Takahashi, J. Fagan, M.-S. Chen: *Formation of a sparse bus impedance
  matrix and its application to short circuit study*, PICA 1973.
- ENTSO-E CGMES profile documentation (EquipmentShortCircuit profile).
- [CGMES Import](cgmes_import.md) — how the source data arrives and is
  validated.
