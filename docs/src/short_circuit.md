# Short-Circuit Analysis

!!! note "Status"
    Sparlectra computes **balanced three-phase short-circuit currents**
    (`Ik''` maximum and minimum, peak current `i_p`) for CGMES deliveries —
    programmatically via [`runShortCircuit!`](@ref) and through the Web UI's
    **Short circuit** button. Unbalanced faults (single line-to-earth,
    line-to-line) are not supported yet. This page explains what
    short-circuit analysis is, how it is commonly done, and how Sparlectra
    does it.

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

**The equations.** Write the (scaled, permuted) factorization as
$A = L\,D\,\tilde U$ with $L$ unit lower triangular, $D$ the pivot
diagonal, and $\tilde U$ unit upper triangular, and let $Z = A^{-1}$.
Multiplying $Z A = I$ and $A Z = I$ out and rearranging gives the two
Erisman-Tinney identities

```math
Z = D^{-1}L^{-1} + (I - \tilde U)\,Z,
\qquad
Z = \tilde U^{-1}D^{-1} + Z\,(I - L).
```

Because $D^{-1}L^{-1}$ is lower triangular and $I - \tilde U$ is strictly
upper triangular, the first identity expresses every diagonal and upper
entry of $Z$ through entries of $Z$ with LARGER row index; the second does
the same for the lower entries through larger column indices:

```math
Z_{jj} = \frac{1}{d_j} - \sum_{k > j} \tilde U_{jk} Z_{kj},
\qquad
Z_{ij} = -\sum_{k > i} \tilde U_{ik} Z_{kj} \;\; (i < j),
\qquad
Z_{ij} = -\sum_{k > j} Z_{ik} L_{kj} \;\; (i > j).
```

Processing columns from $n$ down to $1$ therefore needs no entry that has
not been computed yet. The key structural result is that, restricted to
the sparsity pattern of $(L + U)^{\mathsf T}$ (the FILLED factors), these
sums only ever reference entries inside that same pattern - the recurrence
is self-contained, and the whole selected inverse costs one backward pass
over $\mathrm{nnz}(L) + \mathrm{nnz}(U)$ entries, comparable to the
factorization itself. The diagonal of the original matrix stays inside
that pattern as long as the row and column pivot orders coincide, which is
what UMFPACK's symmetric strategy produces on the structurally symmetric
$Y_{sc}$; the implementation checks exactly this and counts every
out-of-pattern reference defensively.

Sparlectra implements both approaches: per-bus faults use the column
solve, and all-bus sweeps can opt into the Takahashi sparse inverse with
`runShortCircuit!(net; sweep_method = :takahashi)`. One selected-inverse
pass per island then replaces the per-bus solves (measured 34x to 264x
over the serial sweep between 2000 and 16000 buses, growing with size).
The results agree with the default `sweep_method = :solves` to machine
precision (about `1e-15` relative) but not bitwise, which is why the
solve-based sweep remains the default; islands where the method does not
apply (an unsymmetric UMFPACK pivot ordering, a pattern-closure
violation) fall back to column solves automatically. Threaded sweeps
(`runtime.parallel.*`) and the Takahashi pass compose: the selected
inverse removes the per-bus solve cost, the threads cover whatever solves
remain.

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
| `AsynchronousMachine` | `iaIrRatio`, `rxLockedRotorRatio`, `efficiency`, `ratedMechanicalPower`, `polePairNumber`, `ratedS`, `ratedU`, `ratedPowerFactor` | motor contribution to the maximum current (locked-rotor impedance) |

Whether a given delivery could support a calculation is measurable:
`shortCircuitCoverage(result.shortcircuit)` reports, per class, the record
count and per attribute the fill rate; `printShortCircuitCoverage` renders
it, and every CGMES run writes the same view into `cgmes.log`.

**The calculation.** [`runShortCircuit!`](@ref) computes the balanced
three-phase initial symmetrical current:

```julia
result = importCGMES(path = ["grid.zip", "boundary.zip"])
sc = runShortCircuit!(result; case = :max)          # or :min; buses = :all default
sc = runShortCircuit!(result.net, result.shortcircuit; buses = ["Bus_1_220"], case = :min)
printShortCircuitResult(sc)
```

Per fault bus you get `Ik''` (kA), `Sk''` (MVA), and `κ`/`i_p` from the R/X
ratio at the fault location (IEC 60909-0 method b, capped at 1.8 below 1 kV
/ 2.0 above). The positive-sequence short-circuit matrix uses series branch
impedances only — loads, line charging and shunt compensators are dropped
per the standard. Four source types feed the fault:

- **synchronous machines**: `x''_d` on machine base, converted to network
  base, with the §6.6.3 fictitious resistance so the R/X ratio at the fault
  stays meaningful;
- **network feeders** (`ExternalNetworkInjection`): equivalent impedance
  from the declared initial short-circuit current and R/X ratio — at the
  connection point the declared current is reproduced exactly;
- **boundary equivalents** with a declared positive-sequence impedance;
- **asynchronous motors**: locked-rotor impedance
  `Z_M = (1/(I_LR/I_rM)) · U_rM²/S_rM` per §6.7, with `S_rM` taken from the
  rated apparent power or formed from mechanical power, efficiency and
  power factor. Motors raise only the **maximum** current; they never enter
  the minimum case.

Voltage factors follow IEC 60909-0 Table 1 by voltage level (`c_max`
1.05/1.10, `c_min` 0.95/1.00); `short_circuit.c_factor` (or the `c_factor`
keyword) is a scalar expert override for verification runs. The per-bus
path is a Z-bus column solve — one sparse LU factorization per island and
case, one triangular solve per fault bus; the Takahashi sparse inverse
remains the designated optimization for all-bus sweeps on large networks.

**Safety flags — read them.** A short-circuit result is safety-relevant, so
every substituted default and every skipped contribution is **flagged on
the affected result rows themselves** (`contains_defaulted_data` plus a
reason list), not only logged. The documented substitutions: a machine
without `x''_d` gets 0.2 pu on machine base; a feeder without an R/X ratio
gets R = 0.1·X; a motor without a locked-rotor R/X ratio gets the §6.7.2
guidance value (0.10/0.15 for MV motors, 0.42 for LV). A motor or feeder
whose impedance cannot be formed at all is skipped and its whole island is
flagged: the maximum current is then a **lower bound** — the
non-conservative direction for equipment rating. Buses in islands without
any source report `status = :no_source` with `NaN` currents instead of
fabricated numbers.

**Current limitations.** The transformer impedance correction `K_T`
(IEC 60909-0 §6.3.3) and the generator correction `K_G` (§6.6.3) are not
applied — the harvested data does not carry the needed rated power factors;
their absence biases `Ik''` slightly high on transformer-/generator-near
buses. The LV `c_max` variant 1.10 (for +10 % voltage-tolerance bands) and
a configurable per-voltage-level `c` table are not available. Unbalanced
faults (single line-to-earth, line-to-line) and the derived
breaking/thermal quantities are not supported yet; they require the
zero-sequence model whose source data (transformer vector groups and
earthing) is already harvested. The reference tests derive their expected
values analytically from the IEC formulas (full derivations in
`test/test_short_circuit.jl`).

**Web UI.** The PowerFlow form's **Short circuit** button runs both cases
(no power-flow solve involved) and writes `short_circuit_max.csv` /
`short_circuit_min.csv` plus a result-page summary; it is only selectable
for CGMES deliveries that carry short-circuit source data — see
[Web UI](webui.md).

## References

- IEC 60909-0, *Short-circuit currents in three-phase a.c. systems — Part 0:
  Calculation of currents*; IEC 60909-4 for worked verification examples.
- K. Takahashi, J. Fagan, M.-S. Chen: *Formation of a sparse bus impedance
  matrix and its application to short circuit study*, PICA 1973.
- ENTSO-E CGMES profile documentation (EquipmentShortCircuit profile).
- [CGMES Import](cgmes_import.md) — how the source data arrives and is
  validated.
