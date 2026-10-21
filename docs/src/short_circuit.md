# Short-Circuit Analysis

Sparlectra computes balanced three-phase short-circuit currents (`Ik''`
maximum and minimum, peak current `i_p`) via [`runShortCircuit!`](@ref)
and the Web UI's **Short circuit** button. Unbalanced faults (single
line-to-earth, line-to-line) are not supported.

The maximum fault current sizes switchgear, busbars and conductors; the
minimum is the one a protection scheme must still see.

## How it is commonly done: IEC 60909

IEC 60909 is a quasi-static calculation with defined safety margins,
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

so the R/X ratio at the fault matters, not only the impedance magnitude
(the generic factor; Sparlectra reports the method-b form of the result
section). Sparlectra computes $I_k''$, $i_p$ and $\kappa$ (method b);
breaking and thermal equivalent currents are not computed. None of the
three depends on the system frequency.

**Where the fault current comes from.** A synchronous machine feeds the
fault through its subtransient reactance (`x''_d`), a network feeder
through an equivalent impedance from its declared short-circuit power and
R/X ratio; asynchronous motors contribute briefly, mostly to the peak
current. A unit behind a full-size converter (wind, solar, battery) feeds a
current its control limits; IEC 60909-0:2016 models it as a current source.
Neglecting motors or converters underestimates the maximum current, the
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
factors, the entire diagonal among them, follow from one backward pass
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
triangular, so the first identity gives every diagonal and upper entry of
$Z$ from entries with larger row index, the second the lower entries from
larger column indices:

```math
Z_{jj} = \frac{1}{d_j} - \sum_{k > j} \tilde U_{jk} Z_{kj},
\qquad
Z_{ij} = -\sum_{k > i} \tilde U_{ik} Z_{kj} \;\; (i < j),
\qquad
Z_{ij} = -\sum_{k > j} Z_{ik} L_{kj} \;\; (i > j).
```

Columns processed from $n$ down to $1$ need no entry not yet computed.
Restricted to the sparsity pattern of $(L + U)^{\mathsf T}$ (the filled
factors), the sums reference only entries inside that pattern, so the
selected inverse costs one backward pass over
$\mathrm{nnz}(L) + \mathrm{nnz}(U)$ entries. The diagonal stays inside the
pattern while row and column pivot orders coincide, which UMFPACK's
symmetric strategy produces on the structurally symmetric $Y_{sc}$; the
implementation checks this and counts every out-of-pattern reference.

Per-bus faults use the column solve, all-bus sweeps the Takahashi sparse
inverse, one pass per island. `sweep_method = :auto` (default) applies the
pass to islands with at least `short_circuit.takahashi_min_buses` buses
(default 50); `:takahashi` and `:solves` force one method for every
island. Results agree with `:solves` to about `1e-15` relative, not
bitwise; islands where the method does not apply (unsymmetric UMFPACK
pivot ordering, pattern-closure violation) fall back to column solves. The
service and Web UI paths honor `short_circuit.sweep_method`. Threaded
sweeps (`runtime.parallel.*`) compose with the pass.

## How Sparlectra does it

| Use | Where |
|---|---|
| Call | [`runShortCircuit!`](@ref) (`case = :max` or `:min`, `buses = :all` or a list), `printShortCircuitResult`; coverage: `shortCircuitCoverage`, `printShortCircuitCoverage` |
| Config key | `short_circuit.c_factor` (or the `c_factor` keyword): scalar expert override of the voltage factor for verification runs; `short_circuit.sweep_method` (`:auto`, `:takahashi`, `:solves`), `short_circuit.takahashi_min_buses` (default 50); `short_circuit.converter_current_factor` (default 0, off): converter current without `maxIFault` as a multiple of the rated current |
| Result field | per fault bus `Ik''` (kA), `Sk''` (MVA), `κ`, `i_p`, `status` (`:no_source` with `NaN` currents), `contains_defaulted_data` plus a reason list, `lower_bound`; per run `status` and `status_reason` |
| Artifact | `short_circuit_max.csv`, `short_circuit_min.csv` (one row per bus; `flagged`, `reasons` and `lower_bound` carry the flags); coverage view in `run.log` |
| Web UI | PowerFlow form, **Short circuit** button: both cases, no power-flow solve; which cases offer it: [Web UI](webui.md) |

**The data.** Every CGMES import harvests the short-circuit source data
into `CGMESImportResult.shortcircuit` and onto the network
(`net.sc_sources`), so an export or a run of the network alone has it
(read, never altered, CGMES units). The source classes keep the records
of the objects the import mapped onto the network; a skipped one (terminal
disconnected, out of service) feeds no fault current:

| Class | Harvested attributes | Role in the calculation |
|---|---|---|
| `ExternalNetworkInjection` | `maxInitialSymShCCurrent`, `minInitialSymShCCurrent`, `maxR1ToX1Ratio`, `minR1ToX1Ratio`, `maxR0ToX0Ratio`, `maxZ0ToZ1Ratio`, `ikSecond`, `governorSCD` | network-feeder equivalent impedance in the positive and zero sequence |
| `SynchronousMachine` | `satDirectSubtransX` (`x''_d`), `satDirectTransX`, `r0`, `x0`, `r2`, `x2`, `earthing`, `ratedS`, `ratedU` | machine contribution on machine base, sequence impedances, earthing |
| `ACLineSegment` | `r0`, `x0`, `b0ch`, `g0ch`, `shortCircuitEndTemperature` | zero-sequence line model; conductor temperature for minimum-current cases |
| `PowerTransformerEnd` | `r0`, `x0`, `grounded`, `rground`, `xground` | zero-sequence transformer model and star-point treatment |
| `EquivalentInjection` | `r`, `x`, `r0`, `x0`, `r2`, `x2` | boundary equivalents in all three sequence networks |
| `AsynchronousMachine` | `iaIrRatio`, `rxLockedRotorRatio`, `efficiency`, `ratedMechanicalPower`, `polePairNumber`, `ratedS`, `ratedU`, `ratedPowerFactor` | motor contribution to the maximum current (locked-rotor impedance) |
| `PowerElectronicsConnection` (CGMES 3.0) | `maxIFault` (CIM17, per unit of the rated current; in no CGMES 3.0 profile, read where a delivery carries it), `ratedS`, `ratedU` | converter current source in the maximum current |

`shortCircuitCoverage(result.shortcircuit)` reports the record count per
class and the fill rate per attribute; `printShortCircuitCoverage` renders
it, and every CGMES power-flow run writes the same view as a table into
`cgmes_import_report.md`.

**The calculation.** [`runShortCircuit!`](@ref) computes the balanced
three-phase initial symmetrical current:

```julia
result = importCGMES(path = ["grid.zip", "boundary.zip"])
sc = runShortCircuit!(result; case = :max)          # or :min; buses = :all default
sc = runShortCircuit!(result.net, result.shortcircuit; buses = ["Bus_1_220"], case = :min)
printShortCircuitResult(sc)
```

Per fault bus: `Ik''` (kA), `Sk''` (MVA), and `κ`/`i_p` from the R/X ratio
of $Z_k$ at the fault location (IEC 60909-0 method b: the reported value is
$\kappa = \min\bigl(1.15\,(1.02 + 0.98\, e^{-3R/X}),\ \kappa_{\max}\bigr)$
with $\kappa_{\max} = 1.8$ up to and including 1 kV and $2.0$ above). Five
source types feed the fault:

- synchronous machines: `x''_d` on machine base, converted to network
  base, with the §6.6.3 fictitious resistance;
- network feeders (`ExternalNetworkInjection`): equivalent impedance from
  the declared initial short-circuit current and R/X ratio, reproduced
  exactly at the connection point;
- boundary equivalents with a declared positive-sequence impedance;
- asynchronous motors: locked-rotor impedance
  `Z_M = (1/(I_LR/I_rM)) · U_rM²/S_rM` per §6.7, with `S_rM` from the rated
  apparent power or from mechanical power, efficiency and power factor;
  motors enter only the maximum case;
- converter-connected units (CGMES 3.0 `PowerElectronicsConnection`: wind,
  solar and battery parks, mapped with `cgmes_import.map_power_electronics`):
  current sources (IEC 60909-0:2016 clause 6.9), maximum case only;
- synchronous machines on a wind or solar unit (CGMES `SynchronousMachine`
  on a `WindGeneratingUnit` or `SolarGeneratingUnit`, IIDM generator with
  `energySource` `WIND` or `SOLAR`): the same current sources, see below.

A branch built from a CGMES `EquivalentBranch` enters with its series
impedance as delivered, also with a negative `r`: a network equivalent is no
equipment. A leg of a three-winding transformer in its star equivalent does
the same, since only the sum of two legs is physical. The report names the
number of such branches with negative `r` in the islands of the fault buses,
one line for the equivalent branches and one for the three-winding legs.

A converter feeds $I_{sk} = \mathrm{maxIFault} \cdot S_r / (\sqrt{3}\, U_r)$
from its `ratedS` and `ratedU`. Its share is superposed on the equivalent
voltage source. With $Z$ the impedance matrix of the island in which the
voltage sources are short-circuited (the matrix of the first term) and $j$
over the converter buses, the current at fault bus $i$ is

```math
I_k'' = \frac{c \, U_n}{\sqrt{3}\,\lvert Z_{ii} \rvert}
      + \frac{1}{\lvert Z_{ii} \rvert} \sum_j \lvert Z_{ij} \rvert \, I_{sk,j}
```

The magnitudes are added, which bounds the phasor sum from above for any
angle of the converter currents. `κ` applies to the first term only,
$i_p = \kappa \sqrt{2}\, I''_{vs} + \sqrt{2}\, I''_{conv}$; the reported
`zk_ohm`, R/X and `κ` are those of $Z_{ii}$.

An island fed by converters alone has no voltage source and no $Z_{ii}$.
At a bolted fault IEC 60909 neglects loads and shunt branches, so the whole
current of every converter of the island flows into the fault. Every bus
of such an island reports `status = :ok` with the arithmetic sum
$I_k'' = \sum_j I_{sk,j}$ over the converters of the island, referred to
the voltage of the fault bus, and $i_p = \sqrt{2}\, I_k''$ without `κ`;
`zk_ohm`, R/X and `κ` are `NaN`. The sum assumes all currents in phase and
is therefore an upper bound: the row is flagged with the reason
`converter-only island, arithmetic sum, no voltage source`, but not marked
as a lower bound (a skipped unit in the island still marks it). With no
converter feeding (factor off, no `maxIFault`) the island stays
`:no_source`.

`maxIFault` comes from CIM17. No CGMES 3.0 profile (EQ, SC, SSH, DY, OP)
defines it, so a conforming delivery has none and its converters feed
nothing and flag their island (next paragraph). The current source applies
where a delivery carries the attribute beyond the profile, or where the
source data is built in code.

`short_circuit.converter_current_factor` (default `0.0`, off) fills the
gap: a converter without `maxIFault` then feeds that multiple of its rated
current, $I_{sk} = k \cdot S_r / (\sqrt{3}\, U_r)$, and the rows of its
island carry the substitute as a reason instead of the lower-bound flag. A
stated `maxIFault` always wins over the factor; a stated value that is not
positive feeds nothing and flags the island with that value. The factor
is a case key: the configuration file, the case configuration file and an
API override can set it. The keyword of [`runShortCircuit!`](@ref) is
`converter_current_factor`.

The fault current of a converter is a manufacturer figure, not a standard
value, so the factor has no numeric default and has to be set for the
case. As orientation, converters typically feed about 1.0 to 1.5 times
their rated current. Left off, every converter unit feeds nothing and its
island is reported as a lower bound.

A synchronous machine on a wind or solar unit is a converter unit as well:
a CGMES `SynchronousMachine` whose `GeneratingUnit` is a
`WindGeneratingUnit` or `SolarGeneratingUnit`, an IIDM generator with
`energySource` `WIND` or `SOLAR`. It takes the converter path, not the
machine path: $I_{sk} = k \cdot S_r / (\sqrt{3}\, U_r)$ with $k$ =
`short_circuit.converter_current_factor`, $S_r$ from `ratedS`, else from
the maximum active power of the unit (CGMES
`GeneratingUnit.maxOperatingP`, IIDM `maxP`, taken at power factor 1).
A rated power of 9999 MVA and above and a maximum of 9999 MW and above are
placeholders, not sizes. With the factor off or without a size the machine
feeds nothing and flags its island as a lower bound, like a converter
without `maxIFault` or `ratedS`. Its record carries the unit class in the
field `converter_unit`; an SCF file keeps the field in its `sc_source` row,
so the rule holds after an SCF round trip. The CGMES export writes the
unit with its wind or solar class, so the rule also holds after a CGMES
round trip ([CGMES Export](cgmes_export.md)).

Voltage factors follow IEC 60909-0 Table 1 by voltage level (`c_max`
1.05/1.10, `c_min` 0.95/1.00).

**Closed switches and links.** Buses joined by closed links (busbar
couplers, retained CGMES switches) are contracted into one node, the same
groups the power flow contracts. Every bus of a fed group reports the
`Ik''` and the island of that node, also a bus that hangs on the network
through a switch alone. Machines, feeders and motors on any bus of the
group feed the node, and their reasons appear on every row of the group. A
branch between two buses of one group is shorted by the link and left
out. A group without any branch stays `:isolated`.

**Safety flags.** Every substituted default and every skipped contribution
is flagged on the affected result rows (`contains_defaulted_data` plus a
reason list). A machine without `x''_d` gets 0.2 pu on its size (rated
power, else maximum active power, else the network base with a reason of
its own; a maximum of 9999 MW and above and a rated power of 9999 MVA and
above are placeholders, not sizes, and the reason names a placeholder
rated power; a machine with `x''_d` and a placeholder rated power takes the
network base); a
feeder without an R/X ratio gets R = 0.1·X; a motor without a locked-rotor
R/X ratio gets the §6.7.2 guidance value (0.10/0.15 for MV motors, 0.42
for LV). A motor or feeder whose impedance cannot be formed is skipped and
its island flagged: the maximum current is then a lower bound, the
non-conservative direction. A converter without `maxIFault` (with
`short_circuit.converter_current_factor` off) or without `ratedS` (and a converter-connected unit of a network whose source data has no
record for it, and a machine on a wind or solar unit with the factor off
or without a size) feeds nothing and flags its island the same way in the
maximum case; a missing `ratedU` is replaced by the bus nominal voltage
and flagged. Every row of such an island carries `lower_bound = true`
(CSV column `lower_bound`, count `sc_lower_bound_rows` in the run
metadata); a row on a default or a substitute alone is flagged but has
`lower_bound = false`. The minimum case needs no flag, because a missing source
lowers the current. Buses in islands without any source report
`status = :no_source` with `NaN` currents; an island with feeding
converters alone gets their sum (above). A machine without data counts
as a source on the default in every format, so a run is `succeeded` only
when every source carries its data (next section).

### [Source data and run status](@id short_circuit_source_data)

One rule holds for the source data of every format:

| Source data of the case | Run |
|---|---|
| every source carries its quantity | status `succeeded` |
| some sources carry it | status `warning` (`short_circuit_partial_defaults`): the message counts the sources with data, the others enter with default reactances and the rows they feed are flagged with the reason |
| no source carries it, the machines have a rated power | status `warning` (`short_circuit_defaults_only`): a complete table on default reactances, the statement leads the message and `run.log` |
| no source carries anything an impedance can come from | refused with `short_circuit_data_missing` |

A `succeeded` run with flagged rows gives the reason
`short_circuit_flagged_lower_bound` when at least one row is a lower bound
(a skipped source), else `short_circuit_flagged_defaults` (defaults or a
substituted converter current only).

The rule lives in the library: [`runShortCircuit!`](@ref) stores the
status and its reason on the result (`status`, `status_reason`; `failed`
with `short_circuit_data_missing` also when no bus could be evaluated), and
the service reports them unchanged. The three levels and where they show:

| Level | Printed table, `run.log` | CSV | Run metadata | Web UI |
|---|---|---|---|---|
| row: defaulted or substituted data | column `flagged`, reason list below the table | `flagged`, `reasons` | `sc_flagged_rows` | count on the result page, columns of the table |
| island: `Ik''max` is a lower bound | head line (count), reason list | `lower_bound` | `sc_lower_bound_rows` | count on the result page, column of the table |
| run: status and reason | head line | (none, see `result.json`) | `sc_status_line_max`, `sc_status_line_min`; `status`, `reason` in `result.json` | run status card, head line on the result page and above both tables |

`printShortCircuitResult` starts every table with the head line, for
example `run status: warning (short_circuit_partial_defaults), 4 of 4 rows
flagged, 2 rows with a lower bound`.

The quantity is `x''_d` of a synchronous machine, the short-circuit current
of a feeder, the reactance of an equivalent injection, the locked-rotor
ratio of a motor. Converters do not enter this count; one without
`maxIFault` flags the rows of its island instead (with the factor set, as a
substitute).

**Limitations.** The transformer impedance correction `K_T` (IEC 60909-0
§6.3.3) and the generator correction `K_G` (§6.6.3) are not applied (the
harvested data lacks the rated power factors), so `Ik''` is biased
slightly high near transformers and generators. The LV `c_max` variant
1.10 (+10 % voltage-tolerance bands) and a per-voltage-level `c` table are
not available. Unbalanced faults need the zero-sequence model, whose
source data (transformer vector groups and earthing) is already harvested.
The breaking current $I_b$, the thermal equivalent current and $\kappa$ by
method c are not computed.

## References

- IEC 60909-0, *Short-circuit currents in three-phase a.c. systems, Part 0:
  Calculation of currents*; IEC 60909-4 for worked verification examples.
- K. Takahashi, J. Fagan, M.-S. Chen: *Formation of a sparse bus impedance
  matrix and its application to short circuit study*, PICA 1973.
- ENTSO-E CGMES profile documentation (EquipmentShortCircuit profile).
- [CGMES Import](cgmes_import.md): how the source data arrives and is
  validated.
