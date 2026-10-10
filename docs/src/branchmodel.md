# Network Branch and Transformer Model

Branches, transformers and phase-shifting transformers (PSTs) share one
four-terminal equivalent circuit. Typed tap-changer models feed it, an
outer-loop control regulates taps, and the models map onto the ENTSO-E
CGMES data model.

## 1. Common branch model

Each branch uses the same four-terminal model with a complex transformation
ratio $N$ on the from side: $N = 1$ for power lines, a complex value for
transformers. The branch admittance matrix $Y_{br}$ is:

```math
Y_{br} = \begin{bmatrix}
    \frac{1}{\tau^2} \left( y_{ser} + y_{0,from} \right) &
    -y_{ser} \frac{1}{\tau e^{-j\phi}} \\
    -y_{ser} \frac{1}{\tau e^{j\phi}} &
    y_{ser} + y_{0,to}
\end{bmatrix}
```

with `y_ser` the series admittance from the resistance `R` and the
reactance `X`, and $y_{0,from}$, $y_{0,to}$ the shunt arms at the two
terminals (conductance `G`, susceptance `B`):

```math
N = \tau e^{j\phi}
```

```math
y_{ser} = \frac{1}{R + jX}
```

```math
y_{0,from} = G_{from} + jB_{from}, \qquad y_{0,to} = G_{to} + jB_{to}, \qquad
y_{shunt} = y_{0,from} + y_{0,to} = G + jB
```

The branch stores the two arms (`g_from_pu`, `b_from_pu`, `g_to_pu`,
`b_to_pu`) and the totals `g_pu`, `b_pu` as their sums. MATPOWER and the
symmetric pi are the case $y_{0,from} = y_{0,to} = (G + jB)/2$, which every
builder produces when given totals only; the from arm sits behind the
ideal transformer, so a transformer whose magnetizing admittance belongs to
one side keeps it there (see
[Where the magnetizing admittance sits](@ref magnetizing_placement)).

The magnitude $\tau$ is the off-nominal tap ratio, the angle $\phi$ the
phase shift. A pure ratio tap changer moves $\tau$, a pure phase shifter
$\phi$, a combined regulator (German *Schrägregler*) both.

### Circuit diagram

```
                    y_ser
      x--┓┏---------###----------x
         ||   |             |
         ||   # y_0,from    # y_0,to
         ||   |             |
      x--┛┗----------------------x
         N = tau * exp(j * phi)
```

### Sign and conjugation convention

With series admittance $y$, the shunt arms $y_{0,from}$ and $y_{0,to}$, and
from-side tap $t = \tau e^{j\phi}$, the stamped entries (`calcAdmittance`)
are:

```math
\begin{aligned}
Y_{ff} &= \frac{y + y_{0,from}}{\lvert t \rvert^2} \\
Y_{ft} &= -\frac{y}{\overline{t}} \\
Y_{tf} &= -\frac{y}{t} \\
Y_{tt} &= y + y_{0,to}
\end{aligned}
```

With the symmetric split both arms are $0.5\,y_{sh}$ and the entries are
the classic MATPOWER stamp.

The tap is applied on the from side. Reversing a PST with $\tau = 1$ is
electrically a sign flip of $\phi$; with an off-nominal ratio, reversing the
orientation also moves the ratio tap and is not a pure sign change.

## 2. Y-bus assembly

### Diagonal entries (π-model + shunts)

The diagonal Y-bus entry of node $i$ is the nodal self-admittance: the sum
of what every branch at the bus stamps on its terminal, plus the explicit
bus shunt:

```math
Y_{ii} = \sum_{k \in \mathcal{N}(i)} \left( y_{ik} + y_{0,i}^{(ik)} \right) + y_i^{sh}
```

with $y_{ik}$ the series admittance of branch $i-k$, $y_{0,i}^{(ik)}$ its
shunt arm at the terminal $i$ ($y_{ik}^{sh}/2$ for the symmetric split) and
$y_i^{sh}$ the explicit shunt admittance at bus $i$. The line is the case
$t = 1$ of the stamp in section 1; a transformer divides the stamp of its
from terminal by $\lvert t \rvert^2$ ($Y_{ff}$ above). The off-diagonal
entry of a line is

```math
Y_{ik} = -y_{ik}
```

(a transformer: $-y_{ik}/\overline{t}$ and $-y_{ik}/t$, see $Y_{ft}$ and
$Y_{tf}$). For a network of lines the diagonal is therefore the negative
off-diagonal row sum plus the shunt arms and the bus shunt:

```math
Y_{ii} = -\sum_{k \neq i} Y_{ik} + \sum_{k \in \mathcal{N}(i)} y_{0,i}^{(ik)} + y_i^{sh}
```

Only without any shunt (no line charging, no bus shunt) does the row sum
vanish, $Y_{ii} = -\sum_{k \neq i} Y_{ik}$.

The branch builders (`addACLine!`, `addPIModelACLine!`, `addPIModelTrafo!`)
stamp series admittance plus half shunt on each side unless the four
terminal values `g_from_pu`, `b_from_pu`, `g_to_pu`, `b_to_pu` are given
(all four or none); explicit shunts are added as nodal shunt terms when
`bus_shunt_model = "admittance"`.

### [Where the magnetizing admittance sits](@id magnetizing_placement)

A transformer's magnetizing branch ($G$, $B$ of the no-load test) is a
one-sided quantity. PowSyBl and OpenLoadFlow put it wholly on side 1,
behind the ideal transformer (`twtSplitShuntAdmittance = false`), and the
PowSyBl importer keeps it there as the from arm with the to arm zero. A
CGMES `PowerTransformerEnd` carries `g`, `b` per end; the importer refers
end 1 to the to-side base and keeps each end's admittance on its own
terminal. MATPOWER has no place for the split: its reader produces the
symmetric half, its writer puts the symmetric part on the branch and the
excess on a bus shunt named after the branch (see
[MATPOWER](matpower.md)).

### Bus-shunt modeling modes

Bus shunts from sources such as the MATPOWER `Gs`/`Bs` columns enter in
one of two ways:

- `"admittance"` (default): the bus shunt admittance $y_i^{sh} = G_i + jB_i$
  is stamped into the Y-bus diagonal as part of $Y_{ii}$.
- `"voltage_dependent_injection"`: the shunt is not stamped into Y-bus; its
  power is evaluated in the nonlinear injection/mismatch path as a local
  voltage-dependent term.

For a bus shunt admittance $y_i^{sh}$ and local voltage magnitude $|V_i|$:

```math
S_i^{sh} = |V_i|^2 \overline{y_i^{sh}}
```

A positive conductance contributes positive active power; the reactive sign
follows the conjugate of the shunt admittance. The rectangular mismatch uses
`S_calc - S_spec`, so the injection mode subtracts $S_i^{sh}$ from the
specified net injection. Both modes give the same solution and never
double-count: each shunt is either stamped into Y-bus or an injection term.

## 3. One-sided open branches

A branch carries two terminal flags `from_status`/`to_status` next to the
aggregate `status`. The aggregate stays the user-facing switch
(`setBranchStatus!` sets all three; `status = 1` iff both terminals are
closed); `setBranchTerminalStatus!(br; from =, to =)` opens or closes
single terminals. The terminal state is one of `:closed`, `:open_from`,
`:open_to`, `:open`.

### The pi reduction

Take the pi model with series admittance $Y_s = 1/(r + jx)$ and the shunt
arms $Y_{0,from}$ at the closed and $Y_{0,to}$ at the open terminal (both
$(g + jb)/2$ for the symmetric split), terminal `to` open, `from` closed,
no load at the open end. Seen from the closed bus the branch collapses
exactly to the input admittance

```math
Y_{in} = Y_{0,from} + \frac{Y_s Y_{0,to}}{Y_s + Y_{0,to}}
```

Because $|Y_s| \gg |Y_{0,to}|$ for any realistic line, $Y_{in} \approx
Y_{0,from} + Y_{0,to} = g + jb$: the one-sided open line draws its **full**
charging, not half of it, and a transformer whose magnetizing arm sits on
the closed side draws exactly that arm.

The implementation uses the equivalent Schur-complement form on the two-port
from `calcAdmittance`, which already carries the complex ratio: with the `to`
end open $Y_{in} = Y_{11} - Y_{12} Y_{21} / Y_{22}$, with the `from` end open
$Y_{in} = Y_{22} - Y_{21} Y_{12} / Y_{11}$. Lines and transformers
(off-nominal ratio, phase shift) are covered uniformly.

### Open-end voltage and power

The voltage at the open terminal follows from the divider (zero current at
the open end): with `to` open $U_{open} = -Y_{21}/Y_{22} \cdot U_{from}$. It
reproduces the Ferranti rise ($|U_{open}| > |U_{from}|$ for $b > 0$) and is
reported as a branch result (`open_end_vm_pu` / `open_end_va_deg`) without
adding a node to the solved system; the open bus itself is isolated unless
other branches feed it.

The closed terminal carries $S = |U|^2 \cdot \overline{Y_{in}}$: the
imaginary part is the charging reactive power, the real part the power
absorbed by $g$ plus the ohmic loss of the charging current in $r$, small
but not zero. The open terminal carries $S = 0$, and the branch loss equals
the closed-terminal power.

### The equivalent dangling-node formulation

Keeping the full branch and attaching an auxiliary zero-injection PQ bus at
the open end is exactly equivalent for a pure pi branch (the test suite uses
it as the correctness anchor), but it adds one bus per open terminal and
would change bus counts, island reports, result tables and the CGMES
roundtrip identity. Equipment at the open end (a shunt, a load) would need
it and is out of scope.

In the solvers the Y-bus stamps only $Y_{in}$ on the diagonal of the closed
bus, the DC power flow ignores the branch ($B'$ carries no shunts), and the
short-circuit matrix drops it like the charging arms of closed branches.
Results mark partial rows `open@to`/`open@from`, count them under
`Open terminals` in the header, and carry `terminal_state` plus the open-end
voltage in `ACPFlowReport.branches` and the detailed CSV. In the bus table
an isolated open-end bus shows the Ferranti voltage in its V/phi columns,
flagged `open-end` in the Control column (skipped when several
one-sided-open branches end at the bus; an energized bus keeps its solved
voltage). Example: `exp_open_terminal_line.jl`; the basic workshop tour
shows it in chapter 2.

## 4. Tap-changer modelling layers

Source-format parsing, tap-changer semantics, equivalent-circuit
calculation and the solver representation are separate layers:

```text
Importer (MatpowerIO, DTFImporter, ...)
    | maps source fields -> tap-changer model structs, no physics formulas
    v
transformer.jl        data types only (no behaviour)
    v
equicircuit.jl        pure functions: model + tap position -> (ratio, shift_deg, x_pu?)
    v
Branch.ratio / Branch.shift_deg / Branch.x_pu
    v
Y-bus stamping / rectangular NR / outer-loop control
```

Every tap/PST formula lives once, in `equicircuit.jl`. Importers only
construct model structs and call the helpers.

### Data types (`transformer.jl`)

All tap-changer models share the supertype `AbstractTapChangerModel`.

**Ratio tap changer**: `PowerTransformerTaps` carries the tap range and step
definition (`step`, `lowStep`, `highStep`, `neutralStep`,
`voltageIncrement_kV` / `tapStepPercent`, `tapSign`) plus nameplate metadata
(`neutralU`, `neutralU_ratio`). Its `convention` field makes the ratio
convention explicit; `:neutral_relative` applies the correction
`corr = 1 + (step - neutralStep)·tapStepPercent/100` as a divisor on the
winding ratio. Side information lives on `PowerTransformer.tapSideNumber`.

**Phase tap changer**: `PhaseTapChangerModel` classifies the PST technology
via a `kind` field:

```text
kind::Symbol      :symmetrical | :asymmetrical | :tabular
                  (quadrature booster = :asymmetrical with ψ = 90°, no own kind)
step, lowStep, highStep, neutralStep
voltage_step_increment            # per-step voltage increment (linear/nonlinear)
step_phase_shift_increment        # per-step phase increment (linear models)
winding_connection_angle_deg      # ψ, only for :asymmetrical
x_min, x_max                      # X(0), X(αmax) for tap-dependent reactance
table::Union{Nothing,Vector{TapTablePoint}}   # used when kind == :tabular
convention::Symbol
```

**Tap table**: a `PhaseTapChangerModel(kind = :tabular)` is backed by a
non-empty vector of `TapTablePoint` rows (`step`, `ratio`, `angle_deg`,
optional `x_pu`) with strictly ascending, unique steps (validated, never
silently sorted). `lowStep`/`highStep` default to the table range,
`neutralStep` must be a step of the table, and the formula parameters
(`voltage_step_increment`, `winding_connection_angle_deg`, `x_min`, `x_max`)
must be `nothing`. A table overrides the formula path whenever present.

Both kinds attach to a winding: `PowerTransformerWinding` has a
`taps::Union{Nothing,PowerTransformerTaps}` slot and a parallel
`phase_taps::Union{Nothing,PhaseTapChangerModel}` slot, as in CGMES, where a
tap changer hangs on a transformer end.

!!! note "What the winding connection angle ψ means"
    `winding_connection_angle_deg` (ψ) is not a symmetrical-components or
    sequence angle; Sparlectra works in the positive sequence, and ψ lives
    there. ψ is the angle at which a regulator's additional voltage is
    injected relative to the base voltage, the geometry of the regulating
    vector in the complex voltage plane:

    ```math
    \text{regulating vector} = 1 + f \cdot e^{j\psi}, \qquad
    f = (\text{step} - \text{neutralStep}) \cdot u
    ```

    ψ decides how a tap move splits between magnitude and phase: at
    ψ = 0° the added voltage is in phase (a pure longitudinal regulator,
    the regulating vector stays real, `shift_deg` is exactly `0`), at
    ψ = 90° in quadrature (a quadrature booster, mainly a phase shift),
    and for 0° < ψ < 90° ratio and phase change in the proportion set by
    ψ (a combined regulator, Schrägregler).

    From ψ and the tap fraction `f`, `calcPhaseTapAngleRatio` derives the
    effective `ratio` and `shift_deg` stamped into the branch (from-side
    convention: `ratio = 1/|v|`, `shift = -arg(v)`).

### Constructing transformers

Ratio (OLTC) tap changer on a winding:

```julia
taps = PowerTransformerTaps(
  Vn_kV = 110.0,
  step = 0, lowStep = -9, highStep = 9, neutralStep = 0,
  voltageIncrement_kV = 1.1,          # per-step voltage increment
)
# convention defaults to :neutral_relative
```

Symmetrical phase shifter (pure angle regulation):

```julia
pst_sym = PhaseTapChangerModel(
  kind = :symmetrical,
  step = 3, lowStep = -10, highStep = 10, neutralStep = 0,
  voltage_step_increment = 0.012,     # per-step, pu of rated voltage
)
```

Asymmetrical phase shifter / combined regulator (ψ ≠ 0):

```julia
pst_skew = PhaseTapChangerModel(
  kind = :asymmetrical,
  step = 7, lowStep = -13, highStep = 13, neutralStep = 0,
  voltage_step_increment = 0.18 / 13, # e.g. 18 % over 13 steps
  winding_connection_angle_deg = 60.0,
)
# quadrature booster is the same with winding_connection_angle_deg = 90.0
```

Tabular phase shifter (table overrides formulas; carries no formula
parameters):

```julia
table = [
  TapTablePoint(step = -1, ratio = 1.00, angle_deg = -3.0, x_pu = 0.045),
  TapTablePoint(step =  0, ratio = 1.00, angle_deg =  0.0, x_pu = 0.040),
  TapTablePoint(step =  1, ratio = 1.00, angle_deg =  3.0, x_pu = 0.045),
]
pst_tab = PhaseTapChangerModel(
  kind = :tabular,
  neutralStep = 0,                    # must exist in the table
  table = table,                      # lowStep/highStep derived from the table
)
```

Attaching a phase-tap model when building a transformer branch:

```julia
addPIModelTrafo!(
  net = net,
  fromBus = "B1", toBus = "B2",
  r_pu = 0.01, x_pu = 0.08, b_pu = 0.0,
  ratio = 1.0, shift_deg = 0.0, status = 1,
  phase_taps = pst_asym,           # the model drives the branch from here on
)
```

### The resolver

A winding that carries a typed model (`taps`, `phase_taps`) drives its
branch. `resolve_branch_taps!` (`equicircuit.jl`) runs when the transformer
branch is built and whenever a tap controller moves a step: it takes
`ratio` and `shift_deg` of the winding as the neutral point, multiplies the
ratio-tap correction and the phase model's effective ratio onto it, adds
the phase model's shift, derives the branch's tap grid (`tap_min`,
`tap_max`, `tap_step`, `phase_min_deg`, `phase_max_deg`, mean degree per
step) from the model ranges, applies `model.tap_changer_model` to
`r_pu`/`x_pu` from the equipment base, and takes the model's own `X(alpha)`
when it carries `x_min`/`x_max` or a table. The branch is then
`taps_derived`.

With a model present the branch always shows the model's values; when the
values the branch was built with differ by more than 1e-9, one line per
transformer in `net.tapModelNotices` says what was replaced and at which
step, and a service run (Web UI, `run_powerflow_api`) writes those lines to
the artifact `tap_models.log`. Without a model the explicit
`ratio`/`shift_deg` and the legacy degree grid stay as they were.
Controllers move the step of a modelled transformer, never `ratio`/`shift`
directly (probe and update pick the step closest to the proposal); the
state estimator's tap fixation rounds a derived phase regulator to the
model's nearest step.

### Behaviour (`equicircuit.jl`)

The equivalent-circuit helpers turn a model plus a tap position into the
branch quantities, one pure function per formula family:

| Function | Purpose |
|---|---|
| `calcRatioTapCorrection(taps; step)` | ratio-tap multiplicative correction `1 + (step - neutralStep)·tapStepPercent/100` |
| `calcRatioTapRange(taps)` | `(tap_min, tap_max, tap_step)` in ratio terms |
| `calcPhaseTapFraction(m; step)` | shared tap fraction `f = (step - neutralStep)·voltage_step_increment` |
| `calcPhaseTapAngleRatio(m; step)` | `(effective_ratio, effective_shift_deg, regulating_vector)` |
| `calcPhaseTapReactance(m, α)` | tap-angle-dependent reactance `X(α)` interpolated between `x_min`/`x_max` |
| `calcPhaseTapTable(m; step)` | exact lookup of a tabular row |

For `calcPhaseTapAngleRatio`, a `:symmetrical` changer computes
`α = 2·atand(f/2)` with magnitude always `1.0`; an `:asymmetrical` changer
maps the regulating vector `1 + f·e^{jψ}` through the low-level primitive
`calcSkewAngleTap`, of which the quadrature booster (ψ = 90°) is a special
case. A `:tabular` model resolves ratio and angle by lookup and reconstructs
the regulating vector from the stored degrees.

### Importer mapping

- DTF builds a `PhaseTapChangerModel(kind = :asymmetrical,
  winding_connection_angle_deg = skew, ...)` and calls
  `calcPhaseTapAngleRatio` for the branch `ratio`/`shift`. The pure
  longitudinal case (ψ = 0) keeps the shift at exactly `0.0`.
- MATPOWER keeps the direct `TAP`/`SHIFT` path, the CGMES "General Case"
  (raw values), without a model struct; shift unit, shift sign and ratio
  convention are import options ([MATPOWER](matpower.md)).

### Tap-impedance correction

Independent of the typed PST models, `model.tap_changer_model` selects
whether the tap changer of an imported case (MATPOWER and DTF) is ideal
(default) or re-refers the transformer R/X through the tapped winding
(`impedance_correction`, implemented centrally in `calcTapCorrectedRX` /
`calcTapImpedanceCorrectionFactor`); see
[Transformer tap-changer model](configuration.md#Transformer-tap-changer-model).

### Three-winding transformers

A three-winding transformer is a star (T) equivalent with an auxiliary
star-point bus: each winding becomes its own `PowerTransformerWinding`,
stamped as a separate branch to the AUX bus (`create3WTWindings!`, MVA
method). Every winding carries its own `taps` and `phase_taps` slots; a
phase-shifting winding places the regulating vector on its branch to the
star point, with the same stamping and ψ interpretation as a two-winding
device.

`create3WTWindings!` accepts an optional `phase_tap_side` (winding index
`1..3`, `0` = none) and `phase_taps::PhaseTapChangerModel` pair, the same
1-based convention as `tap_side`; `phase_tap_side` may equal `tap_side` when
a winding carries both a ratio tap and a phase tap:

```julia
psc = PhaseTapChangerModel(kind = :asymmetrical, step = 0, lowStep = -8, highStep = 8, neutralStep = 0, winding_connection_angle_deg = 60.0)
w1, w2, w3 = create3WTWindings!(u_kV = [110.0, 20.0, 10.0], sn_MVA = [100.0, 80.0, 20.0], addEx_Side = [tmp1, tmp2, tmp3], sh_deg = [0.0, 0.0, 0.0], tap_side = 1, tap = tapSettings, phase_tap_side = 2, phase_taps = psc)
```

Not implemented: resolving `w2.phase_taps` into an effective ratio/shift on
the AUX-bus branch, and addressing a single 3WT winding from the outer-loop
`PowerTransformerControl` framework.

## 5. Transformer control (outer loop)

Transformers are regulated within the branch PI model using the complex tap
`t = τ·e^{jφ}` and without auxiliary nodes: `τ` for voltage control, `φ`
for active-power-flow (PST) control, both together for combined regulation.
Loop mechanics, deadbands, limits and the hook interface:
[Control Framework](control_framework.md). Branch-model specific are the
complex tap in the PI equivalent and the tap-dependent reactance.

### Tap-dependent reactance X(α)

For a PST whose winding carries a typed `PhaseTapChangerModel` with
reactance data, the outer loop couples the series reactance to the tap
angle: every accepted phase-tap move also updates the branch `x_pu`, and the
next outer-loop solve re-stamps the Y-bus from it. Mapping from the
controller's continuous angle to a reactance:

- formula models (`:symmetrical`/`:asymmetrical` with `x_min = X(0)`,
  `x_max = X(αmax)`): `calcPhaseTapReactance` evaluated at the continuous
  angle;
- tabular models: the nearest table row by angle supplies its per-step
  `x_pu` (no interpolation between rows).

The coupling is opt-in per device: a winding without a typed model, a
formula model without `x_min`/`x_max`, or a tabular model without `x_pu`
values keeps its static reactance (MATPOWER general-case PSTs,
CGMES-flattened PSTs). The probe that estimates the tap direction perturbs
the reactance consistently with the apply step, restores both, and refreshes
the branch flows around each probe solve.

### Discrete tap behaviour

```text
tap_ratio_new       = clamp(tap_ratio ± tap_step,             tap_min,       tap_max)
phase_shift_deg_new = clamp(phase_shift_deg ± phase_step_deg, phase_min_deg, phase_max_deg)
```

### Phase-shift control direction (practical probe)

The sign of a phase-shifter control action is probed on the active model,
not hard-coded:

1. Read `P_ab` at the current phase angle.
2. Move the phase tap one step up (`phase_step_deg`, clamped to the
   range; on a modelled transformer one model step), re-solve, read
   `P_ab` again.
3. The sign of `Delta_P_ab` gives the direction in which `P_ab` grows.
4. If `P_ab < target`, move `phi` in that direction; otherwise the
   opposite way.

See `examples/others/exp_pst_reactance_coupling.jl`.

### Controllers

Registration (`addPowerTransformerControl!` with its full keyword set,
master/slave groups via `followers`, declarative YAML entries) and the
controller result surfaces: [Control Framework](control_framework.md).
Branch-model specific is only the inline form, a controller attached while
the transformer branch is created:

```julia
ctrl = PowerTransformerControl(
  trafo = "",
  mode = :voltage,
  target_bus = "B5",
  target_vm_pu = 1.01,
  control_ratio = true,
  control_phase = false,
)

addPIModelTrafo!(
  net = net,
  fromBus = "B1",
  toBus = "B2",
  r_pu = 0.01,
  x_pu = 0.08,
  b_pu = 0.0,
  ratio = 1.0,
  shift_deg = 0.0,
  status = 1,
  controls = [ctrl],
)
```

### Scope and limits

One `target_bus` measurement, one transformer tap as actuator, a
`target_vm_pu ± deadband` objective; parallel transformers on one bus form
a master/slave group ("Master/slave groups for parallel transformers" in
[Control Framework](control_framework.md)). Not covered: auxiliary
transformer nodes, tap variables inside the Newton iteration,
participation-factor allocation and tap-limit redistribution within
transformer groups.

## 6. CGMES / ENTSO-E mapping

The typed tap-changer models follow the CGMES data model:

| Sparlectra | CGMES / CIM |
|---|---|
| `PowerTransformerTaps` | `RatioTapChanger` on a `TransformerEnd` |
| `PowerTransformerTaps.voltageIncrement_kV` / `tapStepPercent` | `RatioTapChanger.stepVoltageIncrement` |
| `PhaseTapChangerModel(kind = :symmetrical)` | `PhaseTapChangerSymmetrical` |
| `PhaseTapChangerModel(kind = :asymmetrical)` | `PhaseTapChangerAsymmetrical` |
| `winding_connection_angle_deg` (ψ) | `PhaseTapChangerAsymmetrical.windingConnectionAngle` |
| quadrature booster (ψ = 90°) | asymmetrical special case (no separate CIM class) |
| `voltage_step_increment` | `PhaseTapChangerNonLinear.voltageStepIncrement` |
| `step_phase_shift_increment` | `PhaseTapChangerLinear.stepPhaseShiftIncrement` |
| `PhaseTapChangerModel(kind = :tabular)` | `PhaseTapChangerTabular` |
| `TapTablePoint` | `PhaseTapChangerTablePoint` / `TapChangerTablePoint` |
| `TapTablePoint.step / ratio / angle_deg / x_pu` | `TapChangerTablePoint.step / ratio / angle / x` |

CGMES recommends exchanging tabular tap data where available instead of
recalculating parameters from technology formulas: a `:tabular` model
overrides the formula path, and the formula kinds (`:symmetrical`,
`:asymmetrical`) are used only without a table. MATPOWER's raw `TAP`/`SHIFT`
corresponds to the CGMES "General Case".

## 7. Related documents

- ENTSO-E, *Phase Shift Transformers Modelling*, CGMES v2.4, 28 May 2014:
  PST technology classification, tap-angle formulas, reactance-versus-angle
  characteristics, tabular exchange.
- IEC 61970-301 (CIM base) and the CGMES profiles: the classes in the
  mapping table above.
- MATPOWER case format documentation: the `TAP` / `SHIFT` branch columns.
- Sparlectra examples (`examples/others/`): `tap_control_demo_grid.jl`,
  `tap_control_schraeg_two_controllers.jl`,
  `exp_pst_reactance_coupling.jl`; see the
  [examples overview](examples_overview.md).
