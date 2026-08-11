# Slack Bus and External Grid Sources

!!! note "Scope"
    This is a theory page. It explains why the load-flow formulation needs a
    slack bus, what the slack idealizes away, how a network feeder
    (external grid) is modeled per IEC 60909-0 as a voltage source behind an
    internal impedance, **how the equation system changes** when that source
    is taken into the load flow, which additional functions an
    implementation needs, and how each representation enters the
    **short-circuit calculation**. For solver internals see the
    [Solver Guide](solver.md); for the short-circuit method see
    [Short-Circuit Analysis](short_circuit.md). Sparlectra implements this
    page's design as the external-grid element (issue #299):
    `addExternalGrid!` (ideal by default, non-ideal via
    `internal_impedance = true`), `convertSlackToExternalGrid!` for imported
    nets, [`runShortCircuit!`](@ref) directly on a `Net`, the
    `power_flow.external_grid` configuration (see
    [Power-Flow Configuration](powerflow_configuration.md)) and the Web UI's
    **External grid source** fieldset. The runnable comparison lives in
    `examples/powerflow/exp_external_grid_comparison.jl`.

## Why the load flow needs a slack

Consider a network of $n$ buses with nodal admittance matrix
$Y_\mathrm{bus}$. The bus injections follow from the voltages:

```math
I = Y_\mathrm{bus} V, \qquad
S_i(V) = V_i \, \overline{(Y_\mathrm{bus} V)_i}, \qquad i = 1, \dots, n.
```

The load-flow problem asks for the complex voltages $V$ such that the
computed injections match the specified ones, $S_i(V) = S_{\mathrm{spec},i}$.
In rectangular coordinates ($V_i = V_{r,i} + j V_{i,i}$) that is $2n$ real
equations in $2n$ real unknowns — formally square, yet not solvable as
posed, for two structural reasons.

**1. The losses are unknown before the solution.** Summing all bus balances
gives an exact conservation identity: the injections add up to what the
network elements absorb,

```math
\sum_{i=1}^{n} S_i(V) = S_\mathrm{loss}(V),
```

where $S_\mathrm{loss}$ collects the series losses and the shunt/charging
consumption of every branch and shunt. If the injection were specified at
*every* bus, a solution could only exist if the specified values happened to
satisfy $\sum_i S_{\mathrm{spec},i} = S_\mathrm{loss}(V^\ast)$ at the (yet
unknown) solution — a measure-zero coincidence. At least one bus must
therefore leave its injection free, to absorb the a-priori unknown
loss-plus-imbalance term.

**2. The equations fix no angle reference.** All residual quantities used at
PQ and PV buses are invariant under a uniform rotation of every voltage
phasor:

```math
S_i\!\left(e^{j\delta} V\right) = S_i(V), \qquad
\bigl\lvert e^{j\delta} V_i \bigr\rvert = \lvert V_i \rvert .
```

Differentiating at $\delta = 0$ shows that the rotation generator $jV$ lies
in the kernel of the full (unreduced) Jacobian at every point:

```math
J(V)\,(jV) = 0 .
```

A PQ/PV-only system is thus *structurally* singular: if a solution exists at
all, it comes as a one-parameter family of rotated copies, and Newton's
method has no isolated root to converge to. One voltage angle must be fixed
to select a representative.

**The slack bus repairs both defects at once.** One bus — the slack, also
called the reference (REF) bus — has its complex voltage fixed entirely
(magnitude *and* angle), and its two balance equations are dropped from the
system. Fixing the angle removes the rotation kernel (an admissible
variation must leave the slack entries untouched, which $jV$ does not);
dropping the balance equations frees the slack injection to absorb losses
and imbalance. The classical bus-type split is:

| Bus type | Specified | Solved |
|----------|-----------|--------|
| PQ | $P_i$, $Q_i$ | voltage magnitude and angle |
| PV | $P_i$, voltage magnitude | angle, $Q_i$ |
| Slack (REF) | voltage magnitude **and** angle | $P$, $Q$ |

The reduced system has $2(n-1)$ equations in $2(n-1)$ unknowns and is
generically regular. The slack injection is an *output*, evaluated after
convergence:

```math
S_\mathrm{ref} = V_\mathrm{ref}\, \overline{(Y_\mathrm{bus} V)_\mathrm{ref}} .
```

Sparlectra's rectangular solver implements exactly this reduction: the state
vector is $x = [\,V_r(\text{non-slack});\; V_i(\text{non-slack})\,] \in
\mathbb{R}^{2(n-1)}$, and the slack voltage is set once and never solved —
it enters the equations only as data through the $Y_\mathrm{bus}$ coupling
terms of its neighbors (see the [Solver Guide](solver.md)).

## What the ideal slack idealizes away

Because its voltage is held constant no matter what current it supplies, the
slack bus is — read as a physical device — an **ideal voltage source with
zero internal impedance**:

```math
V(I) = U_\mathrm{ref} \quad \text{for every } I
\qquad\Longleftrightarrow\qquad
Z_\mathrm{th} = 0 .
```

No real grid connection behaves like this. The idealization produces three
distinct artifacts:

1. **Infinite voltage stiffness.** The bus voltage shows no reaction to
   loading; the voltage profile in the electrical neighborhood of the
   reference bus comes out systematically too good.
2. **Concentrated balance.** The entire loss-plus-imbalance term lands on
   one machine, although in a real interconnection primary control spreads
   it over many units. This artifact is addressed — orthogonally to this
   page — by the distributed active-power slack described in the
   [Solver Guide](solver.md), which distributes the *balance* role but keeps
   the voltage stiffness ideal.
3. **No usable short-circuit contribution.** The initial symmetrical
   short-circuit current at a fault with impedance $Z_k$ to the source is
   $I_k'' = c\,U_n / (\sqrt{3}\,\lvert Z_k \rvert)$; for $Z_k \to 0$ it
   diverges. An ideal slack therefore carries no finite short-circuit
   datum at all — see the short-circuit section below.

## The external grid as an IEC 60909-0 network feeder

IEC 60909-0 models the connection to a superordinate network — the
**network feeder** (German: *Netzeinspeisung*) — as an ideal source behind
an internal impedance derived from two declared quantities: the initial
symmetrical short-circuit apparent power $S_{kQ}''$ (or equivalently the
current $I_{kQ}''$) at the connection point Q, and the resistance-to-
reactance ratio $R_Q/X_Q$:

```math
I_{kQ}'' = \frac{S_{kQ}''}{\sqrt{3}\; U_{nQ}},
\qquad
Z_Q = \frac{c\, U_{nQ}}{\sqrt{3}\; I_{kQ}''} = \frac{c\, U_{nQ}^2}{S_{kQ}''},
```

with $U_{nQ}$ the nominal voltage at the connection point and $c$ the
voltage factor of the considered case. The impedance splits by the declared
ratio:

```math
X_Q = \frac{Z_Q}{\sqrt{1 + (R_Q/X_Q)^2}},
\qquad
R_Q = \left(\frac{R_Q}{X_Q}\right) X_Q .
```

When no exact ratio is known, IEC 60909-0 permits assuming
$R_Q/X_Q = 0.1$ (equivalently $X_Q = 0.995\, Z_Q$) for high-voltage
feeders; Sparlectra's short-circuit engine applies the same substitution and
**flags** the affected result rows rather than substituting silently.

The finite $Z_Q$ gives the source a finite stiffness. Its terminal
characteristic in normal operation is

```math
V_t = U_\mathrm{ref} - Z_Q\, I_t ,
```

and for a load $P + jQ$ drawn through the feeder the longitudinal component
of the voltage drop is approximately

```math
\Delta V \approx \frac{R_Q\, P + X_Q\, Q}{\lvert V_t \rvert} .
```

A weak grid (small $S_{kQ}''$, large $Z_Q$) shows large voltage swings under
load — precisely the behavior the ideal slack suppresses. The ratio of the
feeder's short-circuit power to the connected load,
$\mathrm{SCR} = S_{kQ}'' / S_\mathrm{load}$, is the usual grid-strength
measure.

**Voltage factor in the load flow.** The factor $c$ is a short-circuit
safety concept (it absorbs the difference between nominal and actual
pre-fault voltage plus a margin). When the feeder impedance is used as an
*operating-point* model inside the load flow, $c = 1$ is the appropriate
choice; $c_\mathrm{max}$/$c_\mathrm{min}$ belong exclusively to the
short-circuit cases.

## Taking the source into the load-flow equations

### Variant 0 — ideal representation (the slack as implemented)

The external grid is represented by making its connection bus the slack:
$V_t \equiv U_\mathrm{ref}$, internal impedance neglected. This is today's
behavior, and it remains the correct default whenever the feeder is strong
relative to the studied network. The short-circuit data ($S_k''$, $R/X$) is
then *carried* alongside the network but changes no load-flow result.

### Variant A — augmented equations (explicit source current)

The direct way to put the source *into* the equations keeps the connection
bus $t$ as an ordinary bus and introduces the source current
$I_s \in \mathbb{C}$ injected at $t$ as a new unknown, coupled by the
source's terminal equation:

```math
U_\mathrm{ref} - Z_Q\, I_s - V_t = 0 .
```

Two things change in the system:

* **The bus-$t$ residual gains a source term.** With the extra injection,
  the power balance at $t$ reads

  ```math
  r_t \;=\; V_t\,\overline{(Y_\mathrm{bus} V)_t}
  \;-\; S_{\mathrm{spec},t}
  \;-\; V_t\, \overline{I_s} \;=\; 0 ,
  ```

  so the residual now depends on the new unknown through the bilinear term
  $V_t \overline{I_s}$.

* **Two new real equations and two new real unknowns.** Splitting the
  terminal equation into real and imaginary parts:

  ```math
  \begin{bmatrix} U_{\mathrm{ref},r} \\ U_{\mathrm{ref},i} \end{bmatrix}
  - \begin{bmatrix} R_Q & -X_Q \\ X_Q & R_Q \end{bmatrix}
    \begin{bmatrix} I_{s,r} \\ I_{s,i} \end{bmatrix}
  - \begin{bmatrix} V_{t,r} \\ V_{t,i} \end{bmatrix}
  = \begin{bmatrix} 0 \\ 0 \end{bmatrix} .
  ```

No bus is eliminated in this formulation: the system has $2n + 2$ real
unknowns ($2n$ voltage components plus $I_{s,r}, I_{s,i}$) and $2n + 2$
equations ($2n$ power balances plus the two constraint rows). It is square
and — unlike the all-PQ system — regular, because the constraint contains
the **fixed phasor** $U_\mathrm{ref}$ and is therefore *not* rotation
invariant: it anchors the angle. Likewise the free source current absorbs
the loss balance. Both classical slack roles are taken over by the source
equations.

The price is structural: the Jacobian gains two new columns (the
derivatives of the bus-$t$ power rows with respect to $I_{s,r}, I_{s,i}$,
from $\partial (V_t \overline{I_s})/\partial I_s$) and two new rows (the
constraint derivatives, $-I_2$ with respect to $V_{t,r}, V_{t,i}$ and the
negative of the $2 \times 2$ impedance block with respect to the current).
These are **new residual and block types**: state indexing, Jacobian
assembly, damping/line-search/trust-region trial states, and the active-set
bookkeeping would all need to learn about the extra state. That is a core
solver modification.

### Variant B — auxiliary internal bus (equivalent, reuses the machinery)

The same physics fits into the *unmodified* formulation by giving the
internal source node an explicit bus. Add one auxiliary bus $q$, fix its
voltage $V_q = U_\mathrm{ref}$ — bus $q$ **is** the slack — and connect it
to the terminal bus $t$ with a series branch of impedance $Z_Q$
(zero shunt admittance). The branch stamps into the admittance matrix in
the standard way,

```math
Y_{qq} \leftarrow Y_{qq} + \frac{1}{Z_Q}, \qquad
Y_{tt} \leftarrow Y_{tt} + \frac{1}{Z_Q}, \qquad
Y_{qt} = Y_{tq} \leftarrow -\frac{1}{Z_Q},
```

and the terminal bus becomes an ordinary solved bus (PQ, or PV if something
else regulates it).

**Equivalence to Variant A.** The branch current is

```math
I_s = \frac{V_q - V_t}{Z_Q} = \frac{U_\mathrm{ref} - V_t}{Z_Q} ,
```

which is exactly the terminal equation of Variant A solved for $I_s$.
Substituting it into the bus-$t$ balance reproduces Variant A's modified
residual; the constraint rows themselves become the ordinary nodal
equations of bus $q$, which the slack reduction then removes. Variant B is
Variant A with the source current eliminated analytically — same solution,
two fewer unknowns, and **no new equation types**.

The bookkeeping comparison:

| Formulation | Buses | Real unknowns | New residual types | Solver changes |
|-------------|-------|---------------|--------------------|----------------|
| Variant 0 — ideal slack | `n` | `2(n-1)` | none | none |
| Variant A — augmented source | `n` | `2n+2` | source-constraint rows, bilinear injection term | new Jacobian blocks, state bookkeeping |
| Variant B — auxiliary bus | `n+1` | `2n` | none | none |

(Variant B has $2((n{+}1)-1) = 2n$ unknowns: one bus more, but the new bus
is the slack and drops out of the state.) Because Variant B changes nothing
in the equation *types*, the entire existing machinery — sparse Y-bus
assembly, analytic rectangular Jacobian, Q-limit active-set switching,
island handling, damping/merit/trust-region step control — applies
unchanged. This is the decisive argument for Variant B in an existing
Newton implementation.

**Per-unit value of the branch.** With system base $S_\mathrm{base}$ and
the voltage base equal to $U_{nQ}$ at the connection point,

```math
z_\mathrm{pu}
= \frac{Z_Q}{Z_\mathrm{base}}
= \frac{c\, U_{nQ}^2 / S_{kQ}''}{U_{nQ}^2 / S_\mathrm{base}}
= \frac{c\, S_\mathrm{base}}{S_{kQ}''} ,
```

with $c = 1$ for the load-flow branch as argued above. Example:
$S_{kQ}'' = 3000\,\mathrm{MVA}$ on a $100\,\mathrm{MVA}$ base gives
$z_\mathrm{pu} = 1/30 \approx 0.0333$, split by $R_Q/X_Q$. If the bus
nominal voltage differs from the voltage base, the general form
$z_\mathrm{pu} = c\, (U_{nQ}/U_\mathrm{base})^2\, S_\mathrm{base}/S_{kQ}''$
applies.

### The stiff limit

For $S_{kQ}'' \to \infty$ the branch impedance vanishes,
$z_\mathrm{pu} = c\,S_\mathrm{base}/S_{kQ}'' \to 0$, and the terminal
voltage converges to the internal one:

```math
\lvert V_t - U_\mathrm{ref} \rvert
= \lvert Z_Q \rvert \, \lvert I_s \rvert \;\longrightarrow\; 0 ,
```

since $I_s$ stays bounded (it converges to the injection current of the
ideal slack). Variant B therefore degenerates *continuously* to Variant 0 —
the ideal slack is the exact limit of the non-ideal source, which makes a
sharp regression test: a huge $S_k''$ must reproduce the ideal-slack
solution to tight tolerance.

Numerically the limit is hostile: the stamped admittance
$1/z_\mathrm{pu} \propto S_{kQ}''$ grows without bound, the rows and
columns of buses $q$ and $t$ become dominated by the $\pm 1/z$ entries, and
the Jacobian's conditioning degrades linearly in $S_k''$. An "almost
ideal" feeder via an artificially enormous $S_k''$ is therefore the wrong
tool in production — use the ideal representation (Variant 0) directly when
ideal behavior is wanted.

### What changes in the results

Compared with the ideal slack, the non-ideal source changes:

* **Terminal voltage and angle become load-dependent.** The angle
  reference now sits on the internal bus $q$, so $V_t$ has a nonzero angle
  and a magnitude below (or above, for reverse flow) $U_\mathrm{ref}$.
* **The source's power depends on the measuring point.** The injection at
  the internal node and the power arriving at the terminal differ by the
  branch loss, $\Delta S = (R_Q + jX_Q)\,\lvert I_s \rvert^2$.
* **The neighborhood sees a realistic voltage profile.** The stiffness
  artifact of the ideal slack disappears; differences decay with electrical
  distance from the connection point.

What does **not** change: the network still has exactly one angle reference
per island, the loss balance is still absorbed by exactly one free
injection, and every other bus keeps its PQ/PV role.

## Required additional functions

Putting the above into a package needs a small, well-defined API surface
(implemented in Sparlectra as issue #299):

1. **A constructor, `addExternalGrid!`.** One call that (a) creates the
   load-flow side — in the ideal variant exactly what a manually added
   slack-type injection produces today, via the existing prosumer path, no
   new solver concept; and (b) converts the declared short-circuit data at
   add time into the feeder record the short-circuit engine already
   consumes: $I_{kQ}''$ from $S_{kQ}''$ and the bus nominal voltage, plus
   the $R/X$ ratios, for the maximum and (optionally) minimum case.
   Validation (positive $S_k''$, $S_{k,\min}'' \le S_{k,\max}''$,
   non-negative ratios) happens here, not in the engine.
2. **A native short-circuit data container on the network.** The engine is
   deliberately duck-typed over the *record contract* of the CGMES
   short-circuit harvest; a field-identical native container
   (`NativeShortCircuitData`) lets programmatically built networks — and
   MATPOWER imports, which carry no short-circuit attributes at all — feed
   the same engine without touching it.
3. **A convenience overload `runShortCircuit!(net; buses, case, c_factor)`**
   that reads the native container, so a natively built network runs the
   identical calculation path as a CGMES delivery.
4. **Copy-path safety.** Every code path that reconstructs or copies a
   network must carry the container along; a regression test copies a net
   with an external grid and runs the short circuit on the copy.
5. **Optionally, the non-ideal load-flow variant** as a flag on the
   constructor (`internal_impedance = true`): create the auxiliary bus and
   the series branch through the *existing* bus/branch APIs, mark the
   auxiliary bus as the slack, and tag both as internal so reports can
   distinguish them.

Equally important is what is **not** needed: with Variant B there are no
new residual types, no new Jacobian blocks, no changes to step control or
active-set logic — no solver function is added or modified. The entire
"how does the equation system change" question is answered inside the
network model.

## Effect on the short-circuit calculation

The short-circuit method (IEC 60909-0, as implemented by
[`runShortCircuit!`](@ref)) replaces the operating state by the **equivalent
voltage source at the fault location**: all sources are removed and
represented only by their internal impedances to ground, loads and line
charging are dropped, and the network reduces to the driving-point
impedance $Z_{ff}$ seen from the fault bus $f$:

```math
I_k''(f) = \frac{c\, U_n(f)}{\sqrt{3}\; \lvert Z_{ff} \rvert} .
```

In Sparlectra the per-island short-circuit matrix is assembled from the
series branch impedances only, and every source contributes a shunt
admittance on the diagonal of its connection bus; the feeder's stamp is

```math
Y_{tt} \leftarrow Y_{tt} + \frac{1}{R_Q + jX_Q} ,
```

with $R_Q, X_Q$ from the equations above, evaluated with the $c$ of the
considered case. Several feeders on one bus stack as parallel admittances —
physically parallel infeeds.

**Consistency of the feeder model.** For a network consisting of a single
feeder, a fault at the connection point sees $Z_{ff} = Z_Q$, and the
declared current is recovered *exactly*:

```math
I_k'' = \frac{c\, U_n}{\sqrt{3}} \cdot \frac{S_{kQ}''}{c\, U_n^2}
      = \frac{S_{kQ}''}{\sqrt{3}\, U_n} = I_{kQ}'' .
```

The voltage factor cancels by construction of the feeder equation — the
model returns precisely the datum the grid operator declared.

**Maximum and minimum case.** The maximum case combines the declared
$S_{k,\max}''$ with $c_\mathrm{max}$ (equipment rating); the minimum case
combines $S_{k,\min}''$ with $c_\mathrm{min}$ (protection sensitivity).
Minimum data is optional: a feeder without declared minimum values is
skipped in the minimum case and the affected rows are flagged — the
existing safety-flag contract of the engine, not a new mechanism.

The two central consequences for the slack-versus-source question:

1. **An ideal slack contributes nothing to a short circuit.** Its
   idealization $Z = 0$ has no finite admittance stamp; taken literally it
   would short the equivalent source and produce an unbounded current. The
   property "this bus is the slack" is load-flow bookkeeping, not
   short-circuit data. An island whose only "source" is the load-flow
   slack therefore reports `status = :no_source` — which is exactly why an
   external-grid element must carry $S_k''$ and $R/X$: without them, no
   fault current can be attributed to the grid connection at all.
2. **The load-flow representation does not influence the short-circuit
   result.** IEC 60909-0 deliberately ignores the load-flow state — the
   pre-fault voltage is replaced by the $c$-factor convention — so it is
   irrelevant whether the load flow modeled the external grid ideally
   (Variant 0) or non-ideally (Variant B). In Variant B the auxiliary
   branch is even *inert* in the short-circuit network: the feeder's
   admittance is stamped at the physical connection bus, the internal bus
   has no path to ground of its own, and a dead-end branch carries no
   fault current. Only two rules must hold: the feeder record stays
   anchored at the physical connection bus, and the fictitious internal
   bus is excluded from fault sweeps (a "fault" there has no physical
   meaning). The differing voltage factors — $c = 1$ in the load-flow
   branch, $c_\mathrm{max}/c_\mathrm{min}$ in the short-circuit stamp —
   never conflict, because the two representations are never combined in
   one calculation.

Peak current $i_p$, the $\kappa$ factor, the Z-bus column solve, and the
flag semantics are described in [Short-Circuit Analysis](short_circuit.md).

## Summary

| Aspect | Ideal slack (Variant 0) | External grid source (Variant B) |
|--------|--------------------------|----------------------------------|
| Nature | boundary condition of the equations | physical model of the grid connection |
| Internal impedance | zero (infinitely stiff) | `c·Un²/Sk''`, finite |
| Terminal voltage | fixed, load-independent | load-dependent drop and angle shift |
| Angle reference | at the connection bus | at the internal (auxiliary) bus |
| Equation system | bus eliminated, `2(n-1)` unknowns | one bus added, `2n` unknowns, no new equation types |
| Solver changes | none | none (that is the point of the auxiliary-bus form) |
| Short-circuit contribution | none possible (no finite datum) | feeder admittance from `Sk''` and `R/X` per IEC 60909-0 |
| Limit relation | — | degenerates to the ideal slack for `Sk'' → ∞` |

## References

* IEC 60909-0:2016, *Short-circuit currents in three-phase a.c. systems —
  Part 0: Calculation of currents* (network feeder model, voltage factor
  $c$, minimum/maximum cases).
* Standard load-flow formulations of the slack/PV/PQ split and the Newton
  power flow: Bergen & Vittal, *Power Systems Analysis*; Oeding & Oswald,
  *Elektrische Kraftwerke und Netze*.
