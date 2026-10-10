# Slack Bus and External Grid Sources

A power flow needs one bus with a fixed voltage; a real grid connection is a
source behind a finite impedance. IEC 60909-0 gives that impedance from the
declared short-circuit data, and an auxiliary-bus formulation carries it into
the power flow and the short-circuit calculation without changing the solver
(`addExternalGrid!`, see
[Implementation in Sparlectra](#Implementation-in-Sparlectra)).

## Why the load flow needs a slack

For $n$ buses with nodal admittance matrix $Y_\mathrm{bus}$ the injections
follow from the voltages:

```math
I = Y_\mathrm{bus} V, \qquad
S_i(V) = V_i \, \overline{(Y_\mathrm{bus} V)_i}, \qquad i = 1, \dots, n.
```

The load flow asks for the voltages with $S_i(V) = S_{\mathrm{spec},i}$. In
rectangular coordinates ($V_k = V_{r,k} + j V_{i,k}$) that is $2n$ real
equations in $2n$ real unknowns: square, yet not solvable as posed, for two
reasons.

**1. The losses are unknown before the solution.** Summing the bus balances
gives

```math
\sum_{i=1}^{n} S_i(V) = S_\mathrm{loss}(V),
```

where $S_\mathrm{loss}$ collects the series losses and the shunt and
charging consumption of every branch and shunt. With the injection
specified at every bus a solution exists only if
$\sum_i S_{\mathrm{spec},i} = S_\mathrm{loss}(V^\ast)$ at the yet unknown
solution. At least one bus must leave its injection free.

**2. The equations fix no angle reference.** The residuals of PQ and PV
buses are invariant under a uniform rotation of every phasor:

```math
S_i\!\left(e^{j\delta} V\right) = S_i(V), \qquad
\bigl\lvert e^{j\delta} V_i \bigr\rvert = \lvert V_i \rvert .
```

Differentiating at $\delta = 0$ puts the rotation generator $jV$ in the
kernel of the full Jacobian at every point:

```math
J(V)\,(jV) = 0 ,
```

with $jV$ read as the real state vector $[-V_i;\, V_r]$. A PQ/PV-only
system is structurally singular: a solution comes as a one-parameter family
of rotated copies, and Newton has no isolated root. One angle must be
fixed.

**The slack bus repairs both defects.** Its complex voltage is fixed
(magnitude and angle) and its two balance equations are dropped: the fixed
angle removes the rotation kernel, the dropped equations free the slack
injection for losses and imbalance.

| Bus type | Specified | Solved |
|----------|-----------|--------|
| PQ | $P_i$, $Q_i$ | voltage magnitude and angle |
| PV | $P_i$, voltage magnitude | angle, $Q_i$ |
| Slack (REF) | voltage magnitude **and** angle | $P$, $Q$ |

The reduced system has $2(n-1)$ equations in $2(n-1)$ unknowns and is
generically regular. The slack injection is an output:

```math
S_\mathrm{ref} = V_\mathrm{ref}\, \overline{(Y_\mathrm{bus} V)_\mathrm{ref}} .
```

The rectangular solver implements this reduction: the state is
$x = [\,V_r(\text{non-slack});\; V_i(\text{non-slack})\,] \in
\mathbb{R}^{2(n-1)}$, and the slack voltage enters only as data through the
$Y_\mathrm{bus}$ coupling of its neighbors ([Solver Guide](solver.md)).

## Reference priority

A case file names its slack. Which unit takes over when there is none, or
when the named one is gone, comes up in three places: an AC island without a
reference, the outage of the slack unit in N-1, and `auto_slack` on a case
without a slack. Every unit carries a reference priority for that
(`ProSumer.referencePriority`, the `referencePriority` keyword of
`addProsumer!` and `addExternalGrid!`) with the semantics of the CGMES
attribute: 0 no preference, 1 the strongest candidate, a larger number
weaker. It orders candidates; it does not make a unit the slack by itself.
The order is the same in all three places:

1. the smallest positive priority wins;
2. between equal priorities, and among units without one,
   [`reference_candidate_rank`](@ref) decides: a network injection before a
   machine, a unit that regulates the voltage of its own bus before one that
   does not, then the larger unit;
3. equal keys go to the smallest bus index.

An island without a reference considers its voltage-controlled units first
and its other generating units only when it has none; a static var
compensator is never the first choice, it carries no active power.

| Source | Priority |
|---|---|
| MATPOWER | type 3 is priority 1 for every unit of the bus, also on a second type-3 bus that stays PV in an island that already has its reference |
| CGMES import | the `referencePriority` of every `SynchronousMachine` and `ExternalNetworkInjection`, as delivered; the slack choice of the import is unchanged ([CGMES Import](cgmes_import.md)) |
| CGMES export | the stored priority of every machine and network injection; a slack unit without one is written with 1 ([CGMES Export](cgmes_export.md)) |
| SCF | `extra.<appliance>.reference_priority`, optional, absent means 0 ([Case Format](scf.md)) |

The island report `ac_islands.csv` names the chosen bus and the reason
(reference priority or strongest unit), the N-1 result names the bus that
took over. `examples/others/exp_island_reference_priority.jl` shows it on
the shipped two-island case `data/scf/two_islands_prio.scf.json`.

## What the ideal slack idealizes away

A bus whose voltage is constant whatever current it supplies is an **ideal
voltage source with zero internal impedance**:

```math
V(I) = U_\mathrm{ref} \quad \text{for every } I
\qquad\Longleftrightarrow\qquad
Z_\mathrm{th} = 0 .
```

No real grid connection behaves like this. Three artifacts follow:

1. **Infinite voltage stiffness.** The bus voltage shows no reaction to
   loading; the profile around the reference bus comes out too good.
2. **Concentrated balance.** The whole loss-plus-imbalance term lands on one
   machine, although primary control spreads it over many. The distributed
   active-power slack in the [Solver Guide](solver.md) addresses the
   balance and keeps the voltage stiffness ideal.
3. **No short-circuit contribution.** The initial symmetrical
   short-circuit current at a fault with impedance $Z_k$ to the source is
   $I_k'' = c\,U_n / (\sqrt{3}\,\lvert Z_k \rvert)$; for $Z_k \to 0$ it
   diverges. An ideal slack carries no finite short-circuit datum
   ([Effect on the short-circuit calculation](#Effect-on-the-short-circuit-calculation)).

## The external grid as an IEC 60909-0 network feeder

IEC 60909-0 models the connection to a superordinate network, the
**network feeder** (German: *Netzeinspeisung*), as an ideal source behind an
impedance derived from two declared quantities: the initial symmetrical
short-circuit apparent power $S_{kQ}''$ (or current $I_{kQ}''$) at the
connection point Q, and the ratio $R_Q/X_Q$:

```math
I_{kQ}'' = \frac{S_{kQ}''}{\sqrt{3}\; U_{nQ}},
\qquad
Z_Q = \frac{c\, U_{nQ}}{\sqrt{3}\; I_{kQ}''} = \frac{c\, U_{nQ}^2}{S_{kQ}''},
```

with $U_{nQ}$ the nominal voltage at the connection point and $c$ the
voltage factor of the case. The impedance splits by the ratio:

```math
X_Q = \frac{Z_Q}{\sqrt{1 + (R_Q/X_Q)^2}},
\qquad
R_Q = \left(\frac{R_Q}{X_Q}\right) X_Q .
```

Without a known ratio, IEC 60909-0 permits $R_Q/X_Q = 0.1$
($X_Q = 0.995\, Z_Q$) for high-voltage feeders; the short-circuit engine
substitutes the same value and flags the affected result rows.

The finite $Z_Q$ gives the source a finite stiffness. Its terminal
characteristic in normal operation is

```math
V_t = U_\mathrm{ref} - Z_Q\, I_t ,
```

and for a load $P + jQ$ drawn through the feeder the longitudinal voltage
drop is approximately

```math
\Delta V \approx \frac{R_Q\, P + X_Q\, Q}{\lvert V_t \rvert} .
```

A weak grid (small $S_{kQ}''$, large $Z_Q$) shows large voltage swings under
load, which the ideal slack suppresses. The usual grid-strength measure is
$\mathrm{SCR} = S_{kQ}'' / S_\mathrm{load}$.

**Voltage factor in the load flow.** The feeder impedance in the load flow
uses $c = 1$; $c_\mathrm{max}$ and $c_\mathrm{min}$ are short-circuit safety
margins (nominal versus pre-fault voltage) and belong to the short-circuit
cases.

## Taking the source into the load-flow equations

### Variant 0: ideal representation (the slack as implemented)

The connection bus is the slack: $V_t \equiv U_\mathrm{ref}$, internal
impedance neglected. This is the default and correct whenever the feeder is
strong relative to the studied network. The short-circuit data ($S_k''$,
$R/X$) is carried alongside and changes no load-flow result.

### Variant A: augmented equations (explicit source current)

Keep the connection bus $t$ as an ordinary bus and add the source current
$I_s \in \mathbb{C}$ injected at $t$ as a new unknown, coupled by the
terminal equation

```math
U_\mathrm{ref} - Z_Q\, I_s - V_t = 0 .
```

The bus-$t$ residual gains a source term, bilinear in the new unknown,

```math
r_t \;=\; V_t\,\overline{(Y_\mathrm{bus} V)_t}
\;-\; S_{\mathrm{spec},t}
\;-\; V_t\, \overline{I_s} \;=\; 0 ,
```

and the terminal equation, split into real and imaginary parts, adds two
real equations for the two new real unknowns:

```math
\begin{bmatrix} U_{\mathrm{ref},r} \\ U_{\mathrm{ref},i} \end{bmatrix}
- \begin{bmatrix} R_Q & -X_Q \\ X_Q & R_Q \end{bmatrix}
  \begin{bmatrix} I_{s,r} \\ I_{s,i} \end{bmatrix}
- \begin{bmatrix} V_{t,r} \\ V_{t,i} \end{bmatrix}
= \begin{bmatrix} 0 \\ 0 \end{bmatrix} .
```

No bus is eliminated: $2n + 2$ unknowns and $2n + 2$ equations. Unlike the
all-PQ system it is regular: the constraint contains the fixed phasor
$U_\mathrm{ref}$ and anchors the angle, and the free source current absorbs
the loss balance. The price: two new Jacobian columns (derivatives of the
bus-$t$ rows with respect to $I_{s,r}, I_{s,i}$) and two new rows ($-I_2$
with respect to $V_{t,r}, V_{t,i}$, the negative $2 \times 2$ impedance
block with respect to the current), new residual and block types for state
indexing, Jacobian assembly, step control and active-set bookkeeping.

### Variant B: auxiliary internal bus (equivalent, reuses the machinery)

The same physics fits the unmodified formulation when the internal source
node gets its own bus. Add one auxiliary bus $q$, fix
$V_q = U_\mathrm{ref}$ (bus $q$ **is** the slack) and connect it to the
terminal bus $t$ with a series branch of impedance $Z_Q$ and zero shunt
admittance. The branch stamps into the admittance matrix as usual,

```math
Y_{qq} \leftarrow Y_{qq} + \frac{1}{Z_Q}, \qquad
Y_{tt} \leftarrow Y_{tt} + \frac{1}{Z_Q}, \qquad
Y_{qt} = Y_{tq} \leftarrow -\frac{1}{Z_Q},
```

and the terminal bus becomes an ordinary solved bus (PQ, or PV if something
else regulates it).

**Equivalence to Variant A.** The branch current

```math
I_s = \frac{V_q - V_t}{Z_Q} = \frac{U_\mathrm{ref} - V_t}{Z_Q} ,
```

is Variant A's terminal equation solved for $I_s$. Substituted into the
bus-$t$ balance it reproduces Variant A's residual; the constraint rows
become the nodal equations of bus $q$, which the slack reduction removes.
Variant B is Variant A with the source current eliminated analytically:
same solution, two fewer unknowns, no new equation types.

| Formulation | Buses | Real unknowns | New residual types | Solver changes |
|-------------|-------|---------------|--------------------|----------------|
| Variant 0, ideal slack | `n` | `2(n-1)` | none | none |
| Variant A, augmented source | `n` | `2n+2` | source-constraint rows, bilinear injection term | new Jacobian blocks, state bookkeeping |
| Variant B, auxiliary bus | `n+1` | `2n` | none | none |

(Variant B has $2((n{+}1)-1) = 2n$ unknowns: the new bus is the slack and
drops out of the state.) With unchanged equation types, Y-bus assembly,
Jacobian, Q-limit switching, island handling and step control apply
unchanged.

**Per-unit value of the branch.** With system base $S_\mathrm{base}$ and
the voltage base $U_{nQ}$ at the connection point,

```math
z_\mathrm{pu}
= \frac{Z_Q}{Z_\mathrm{base}}
= \frac{c\, U_{nQ}^2 / S_{kQ}''}{U_{nQ}^2 / S_\mathrm{base}}
= \frac{c\, S_\mathrm{base}}{S_{kQ}''} ,
```

with $c = 1$ for the load-flow branch. Example: $S_{kQ}'' = 3000\,\mathrm{MVA}$
on a $100\,\mathrm{MVA}$ base gives $z_\mathrm{pu} = 1/30 \approx 0.0333$,
split by $R_Q/X_Q$. The impedance base is the nominal voltage of the
connection bus itself ($Z_\mathrm{base} = U_{nQ}^2 / S_\mathrm{base}$), so
no $(U_{nQ}/U_\mathrm{base})^2$ factor applies.

### The stiff limit

For $S_{kQ}'' \to \infty$ the branch impedance vanishes,
$z_\mathrm{pu} = c\,S_\mathrm{base}/S_{kQ}'' \to 0$, and the terminal
voltage converges to the internal one:

```math
\lvert V_t - U_\mathrm{ref} \rvert
= \lvert Z_Q \rvert \, \lvert I_s \rvert \;\longrightarrow\; 0 ,
```

since $I_s$ stays bounded (it converges to the injection current of the
ideal slack). Variant B degenerates continuously to Variant 0 (a regression
test: a huge $S_k''$ must reproduce the ideal-slack solution). Numerically
the limit is hostile: the stamped
admittance $1/z_\mathrm{pu} \propto S_{kQ}''$ dominates the rows and columns
of buses $q$ and $t$, and the Jacobian's conditioning degrades linearly in
$S_k''$. Do not emulate an ideal feeder with an enormous $S_k''$; use
Variant 0.

### What changes in the results

* The terminal voltage and angle become load-dependent: the angle reference
  sits on the internal bus $q$, so $V_t$ has a nonzero angle and a magnitude
  below (or above, for reverse flow) $U_\mathrm{ref}$.
* The source's power depends on the measuring point: injection at the
  internal node and power arriving at the terminal differ by the branch
  loss $\Delta S = (R_Q + jX_Q)\,\lvert I_s \rvert^2$.
* The neighborhood sees a realistic voltage profile; the differences decay
  with electrical distance from the connection point.

Unchanged: one angle reference per island, the loss balance absorbed by one
free injection, every other bus keeps its PQ/PV role.

## Implementation in Sparlectra

| Use | Where |
|---|---|
| Call | `addExternalGrid!` (ideal by default, Variant B with `internal_impedance = true`; validates the short-circuit data at add time); `convertSlackToExternalGrid!` for imported networks; `runShortCircuit!(net; buses, case, c_factor)` on a native network with the same engine as a CGMES delivery |
| Config key | `power_flow.external_grid` ([Power-Flow Configuration](powerflow_configuration.md)) |
| Web UI | fieldset **External grid source** |
| Example | `examples/powerflow/exp_external_grid_comparison.jl` |

## Effect on the short-circuit calculation

The IEC 60909-0 method ([`runShortCircuit!`](@ref)) replaces the operating
state by the equivalent voltage source at the fault location: sources enter
only through their internal impedances to ground, loads and line charging
are dropped, and the network reduces to the driving-point impedance
$Z_{ff}$ seen from the fault bus $f$:

```math
I_k''(f) = \frac{c\, U_n(f)}{\sqrt{3}\; \lvert Z_{ff} \rvert} .
```

The per-island short-circuit matrix is assembled from the series branch
impedances; every source contributes a shunt admittance on the diagonal of
its connection bus. The feeder's stamp is

```math
Y_{tt} \leftarrow Y_{tt} + \frac{1}{R_Q + jX_Q} ,
```

with $R_Q, X_Q$ from the equations above at the $c$ of the case. Several
feeders on one bus stack as parallel admittances.

**Consistency.** For a network of a single feeder, a fault at the connection
point sees $Z_{ff} = Z_Q$ and recovers the declared current exactly:

```math
I_k'' = \frac{c\, U_n}{\sqrt{3}} \cdot \frac{S_{kQ}''}{c\, U_n^2}
      = \frac{S_{kQ}''}{\sqrt{3}\, U_n} = I_{kQ}'' .
```

The voltage factor cancels by construction.

**Maximum and minimum case.** The maximum case combines $S_{k,\max}''$ with
$c_\mathrm{max}$ (equipment rating), the minimum case $S_{k,\min}''$ with
$c_\mathrm{min}$ (protection sensitivity). Minimum data is optional: a
feeder without it is skipped in the minimum case and the affected rows are
flagged.

Two consequences:

1. **An ideal slack contributes nothing to a short circuit.** $Z = 0$ has no
   finite admittance stamp; taken literally it would short the equivalent
   source. "This bus is the slack" is load-flow bookkeeping, not
   short-circuit data: an island whose only source is the load-flow slack
   reports `status = :no_source`, which is why an external-grid element
   must carry $S_k''$ and $R/X$.
2. **The load-flow representation does not influence the short-circuit
   result.** IEC 60909-0 ignores the load-flow state (the pre-fault voltage
   is replaced by the $c$-factor convention), so Variant 0 and Variant B
   give the same short-circuit result. In Variant B the auxiliary branch is
   inert: the feeder admittance is stamped at the physical connection bus,
   the internal bus has no path to ground of its own, and a dead-end branch
   carries no fault current. Two rules hold: the feeder record stays
   anchored at the physical connection bus, and the internal bus is
   excluded from fault sweeps. The differing voltage factors ($c = 1$ in
   the load-flow branch, $c_\mathrm{max}/c_\mathrm{min}$ in the stamp) never
   meet in one calculation.

Peak current $i_p$, the $\kappa$ factor, the Z-bus column solve and the
flag semantics: [Short-Circuit Analysis](short_circuit.md).

## Summary

| Aspect | Ideal slack (Variant 0) | External grid source (Variant B) |
|--------|--------------------------|----------------------------------|
| Nature | boundary condition of the equations | physical model of the grid connection |
| Internal impedance | zero (infinitely stiff) | `Un²/Sk''`, finite (`c = 1` in the load flow; `c_max`/`c_min` only in the short-circuit stamp) |
| Terminal voltage | fixed, load-independent | load-dependent drop and angle shift |
| Angle reference | at the connection bus | at the internal (auxiliary) bus |
| Equation system | bus eliminated, `2(n-1)` unknowns | one bus added, `2n` unknowns, no new equation types |
| Solver changes | none | none (that is the point of the auxiliary-bus form) |
| Short-circuit contribution | none possible (no finite datum) | feeder admittance from `Sk''` and `R/X` per IEC 60909-0 |
| Limit relation | none | degenerates to the ideal slack for `Sk'' → ∞` |

## References

* IEC 60909-0:2016, *Short-circuit currents in three-phase a.c. systems,
  Part 0: Calculation of currents* (network feeder model, voltage factor
  $c$, minimum/maximum cases).
* Bergen & Vittal, *Power Systems Analysis*; Oeding & Oswald,
  *Elektrische Kraftwerke und Netze* (slack/PV/PQ split, Newton power flow).
