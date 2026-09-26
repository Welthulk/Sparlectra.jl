# Short-Circuit Analysis

Sparlectra computes balanced three-phase short-circuit currents (`Ik''`
maximum and minimum, peak current `i_p`) for CGMES deliveries, via
[`runShortCircuit!`](@ref) and the Web UI's **Short circuit** button.
Unbalanced faults (single line-to-earth, line-to-line) are not supported.

Short-circuit analysis asks what current flows when the network faults:
what the switchgear must interrupt, what busbars and conductors must
withstand, and whether protection trips selectively. The maximum fault
current sizes the equipment; the minimum is the one a protection scheme
must still see.

## How it is commonly done: IEC 60909

IEC 60909 is a quasi-static calculation with defined safety margins:
reproducible and conservative without dynamic machine models.

**Equivalent voltage source at the fault location.** The network is
reduced to the short-circuit impedance `Zk` seen from the faulted bus,
driven by a single equivalent source; other sources are removed, loads and
line charging neglected. The initial symmetrical short-circuit current at
nominal voltage `Un` is

```math
I_k'' = \frac{c \, U_n}{\sqrt{3}\,\lvert Z_k \rvert}
```

**The voltage factor** `c` absorbs the difference between nominal and
pre-fault voltage plus a safety margin, tabulated by voltage level and
objective: `c_max` (for example 1.10 in HV networks) for maximum currents,
`c_min` (0.95 or 1.00) for minimum currents.

**Derived quantities.** The peak current (first crest, mechanical stress)
is

```math
i_p = \kappa \,\sqrt{2}\; I_k'',
\qquad
\kappa \approx 1.02 + 0.98\, e^{-3R/X}
```

so the R/X ratio at the fault matters, not just the impedance magnitude.
Breaking and thermal equivalent currents follow as well.

**Where the fault current comes from.** A synchronous machine feeds the
fault through its subtransient reactance (`x''_d`), a network feeder
through an equivalent impedance from its declared short-circuit power and
R/X ratio; asynchronous motors contribute briefly, mostly to the peak
current. Neglecting motors underestimates the maximum current, the
non-conservative direction, so an implementation must say when it does it.

**Unbalanced faults.** Single line-to-earth (the most common real fault),
line-to-line and two-phase-to-earth faults use symmetrical components;
transformer vector group and star-point earthing decide whether a
zero-sequence path exists at all.

## How the network is solved: Z-bus, and Takahashi

For one fault location the needed quantity is one diagonal element of the
bus impedance matrix: with `Y_sc` the short-circuit admittance matrix
(series impedances only, source admittances on the diagonal, no loads or
shunts), the fault impedance at bus `i` is

```math
Z_{ii} = \left( Y_{sc}^{-1} \right)_{ii}
```

The **Z-bus column solve** factorizes `Y_sc` once (sparse LU) and obtains
column `i` of the inverse with one triangular solve against the `i`-th
unit vector; the column holds `Z_ii` and the off-diagonal entries for
voltage sag and branch contributions.

**What the Takahashi method adds.** Takahashi, Fagan and Chen (1973)
showed that all inverse elements inside the sparsity pattern of the LU
factors, in particular the entire diagonal, follow from one backward pass
over the factors at a cost comparable to the factorization: an all-bus
sweep scales like one factorization, not like `n` solves.

**The equations.** Write the (scaled, permuted) factorization as
$A = L\,D\,\tilde U$ with $L$ unit lower triangular, $D$ the pivot
diagonal, $\tilde U$ unit upper triangular, and let $Z = A^{-1}$.
Multiplying $Z A = I$ and $A Z = I$ out gives the Erisman-Tinney
identities

```math
Z = D^{-1}L^{-1} + (I - \tilde U)\,Z,
\qquad
Z = \tilde U^{-1}D^{-1} + Z\,(I - L).
```

$D^{-1}L^{-1}$ is lower triangular and $I - \tilde U$ strictly upper
triangular, so the first identity expresses every diagonal and upper entry
of $Z$ through entries with larger row index; the second does the same for
the lower entries through larger column indices:

```math
Z_{jj} = \frac{1}{d_j} - \sum_{k > j} \tilde U_{jk} Z_{kj},
\qquad
Z_{ij} = -\sum_{k > i} \tilde U_{ik} Z_{kj} \;\; (i < j),
\qquad
Z_{ij} = -\sum_{k > j} Z_{ik} L_{kj} \;\; (i > j).
```

Processing columns from $n$ down to $1$ needs no entry that is not yet
computed. Restricted to the sparsity pattern of $(L + U)^{\mathsf T}$ (the
filled factors), the sums only reference entries inside that pattern, so
the selected inverse costs one backward pass over
$\mathrm{nnz}(L) + \mathrm{nnz}(U)$ entries. The diagonal stays inside the
pattern as long as row and column pivot orders coincide, which UMFPACK's
symmetric strategy produces on the structurally symmetric $Y_{sc}$; the
implementation checks this and counts every out-of-pattern reference.

Per-bus faults use the column solve, all-bus sweeps the Takahashi sparse
inverse, one pass per island. `sweep_method = :auto` (default) applies the pass
to islands with at least `short_circuit.takahashi_min_buses` buses
(default 50); `:takahashi` and `:solves` force one method for every
island. Results agree with `:solves` to about `1e-15` relative, not
bitwise; islands where the method does not apply (unsymmetric UMFPACK
pivot ordering, pattern-closure violation) fall back to column solves.
The service and Web UI paths honor `short_circuit.sweep_method`. Threaded
sweeps (`runtime.parallel.*`) compose with the pass.

## How Sparlectra does it

| Use | Where |
|---|---|
| Call | [`runShortCircuit!`](@ref) (`case = :max` or `:min`, `buses = :all` or a list), `printShortCircuitResult`; coverage: `shortCircuitCoverage`, `printShortCircuitCoverage` |
| Config key | `short_circuit.c_factor` (or the `c_factor` keyword): scalar expert override of the voltage factor for verification runs; `short_circuit.sweep_method` (`:auto`, `:takahashi`, `:solves`), `short_circuit.takahashi_min_buses` (default 50) |
| Result field | per fault bus `Ik''` (kA), `Sk''` (MVA), `κ`, `i_p`, `status` (`:no_source` with `NaN` currents), `contains_defaulted_data` plus a reason list |
| Artifact | `short_circuit_max.csv`, `short_circuit_min.csv`; coverage view in `cgmes.log` |
| Web UI | PowerFlow form, **Short circuit** button: both cases, no power-flow solve; offered for CGMES deliveries with short-circuit source data ([Web UI](webui.md)) |

**The data.** Every CGMES import harvests the short-circuit source data
into `CGMESImportResult.shortcircuit` (read, never altered, CGMES units):

| Class | Harvested attributes | Role in the calculation |
|---|---|---|
| `ExternalNetworkInjection` | `maxInitialSymShCCurrent`, `minInitialSymShCCurrent`, `maxR1ToX1Ratio`, `minR1ToX1Ratio`, `maxR0ToX0Ratio`, `maxZ0ToZ1Ratio`, `ikSecond`, `governorSCD` | network-feeder equivalent impedance in the positive and zero sequence |
| `SynchronousMachine` | `satDirectSubtransX` (`x''_d`), `satDirectTransX`, `r0`, `x0`, `r2`, `x2`, `earthing`, `ratedS`, `ratedU` | machine contribution on machine base, sequence impedances, earthing |
| `ACLineSegment` | `r0`, `x0`, `b0ch`, `g0ch`, `shortCircuitEndTemperature` | zero-sequence line model; conductor temperature for minimum-current cases |
| `PowerTransformerEnd` | `r0`, `x0`, `grounded`, `rground`, `xground` | zero-sequence transformer model and star-point treatment |
| `EquivalentInjection` | `r`, `x`, `r0`, `x0`, `r2`, `x2` | boundary equivalents in all three sequence networks |
| `AsynchronousMachine` | `iaIrRatio`, `rxLockedRotorRatio`, `efficiency`, `ratedMechanicalPower`, `polePairNumber`, `ratedS`, `ratedU`, `ratedPowerFactor` | motor contribution to the maximum current (locked-rotor impedance) |

`shortCircuitCoverage(result.shortcircuit)` reports per class the record
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

Per fault bus: `Ik''` (kA), `Sk''` (MVA), and `κ`/`i_p` from the R/X ratio
at the fault location (IEC 60909-0 method b, capped at 1.8 below 1 kV /
2.0 above). Four source types feed the fault:

- **synchronous machines**: `x''_d` on machine base, converted to network
  base, with the §6.6.3 fictitious resistance;
- **network feeders** (`ExternalNetworkInjection`): equivalent impedance
  from the declared initial short-circuit current and R/X ratio,
  reproduced exactly at the connection point;
- **boundary equivalents** with a declared positive-sequence impedance;
- **asynchronous motors**: locked-rotor impedance
  `Z_M = (1/(I_LR/I_rM)) · U_rM²/S_rM` per §6.7, with `S_rM` from the rated
  apparent power or from mechanical power, efficiency and power factor;
  motors enter only the maximum case.

Voltage factors follow IEC 60909-0 Table 1 by voltage level (`c_max`
1.05/1.10, `c_min` 0.95/1.00).

**Safety flags.** Every substituted default and every skipped contribution
is flagged on the affected result rows (`contains_defaulted_data` plus a
reason list). Substitutions: a machine without `x''_d` gets 0.2 pu on
machine base; a feeder without an R/X ratio gets R = 0.1·X; a motor
without a locked-rotor R/X ratio gets the §6.7.2 guidance value
(0.10/0.15 for MV motors, 0.42 for LV). A motor or feeder whose impedance
cannot be formed is skipped and its island flagged: the maximum current is
then a lower bound, the non-conservative direction. Buses in islands
without any source report `status = :no_source` with `NaN` currents.

**Limitations.** The transformer impedance correction `K_T` (IEC 60909-0
§6.3.3) and the generator correction `K_G` (§6.6.3) are not applied (the
harvested data lacks the rated power factors), so `Ik''` is biased
slightly high near transformers and generators. The LV `c_max` variant
1.10 (+10 % voltage-tolerance bands) and a per-voltage-level `c` table are
not available. Unbalanced faults and the breaking/thermal quantities need
the zero-sequence model, whose source data (transformer vector groups and
earthing) is already harvested. The reference tests derive their expected
values from the IEC formulas (`test/test_short_circuit.jl`).

## References

- IEC 60909-0, *Short-circuit currents in three-phase a.c. systems, Part 0:
  Calculation of currents*; IEC 60909-4 for worked verification examples.
- K. Takahashi, J. Fagan, M.-S. Chen: *Formation of a sparse bus impedance
  matrix and its application to short circuit study*, PICA 1973.
- ENTSO-E CGMES profile documentation (EquipmentShortCircuit profile).
- [CGMES Import](cgmes_import.md): how the source data arrives and is
  validated.
