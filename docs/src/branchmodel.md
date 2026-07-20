# Network Branch and Transformer Model

This page describes how Sparlectra represents branches, transformers and
phase-shifting transformers (PSTs): the common four-terminal equivalent
circuit, the typed tap-changer models that feed it, the outer-loop control that
regulates taps, and how these map onto the ENTSO-E CGMES data model.

## 1. Common branch model

Each branch is treated with the same four-terminal network model, with a
complex transformation ratio $N$ applied on the from side. For power lines the
transmission ratio $N$ is set to 1. For transformers the transformation ratio
$N$ is a complex value.

The branch admittance matrix $Y_{br}$ is:

```math
Y_{br} = \begin{bmatrix}
    \frac{1}{\tau^2} \left( y_{ser} + \frac{y_{shunt}}{2} \right) &
    -y_{ser} \frac{1}{\tau e^{-j\phi}} \\
    -y_{ser} \frac{1}{\tau e^{j\phi}} &
    y_{ser} + \frac{y_{shunt}}{2}
\end{bmatrix}
```

where:
- `y_ser` is the series admittance,
- `y_shunt` is the total branch shunt admittance,
- `R` is the resistance component,
- `X` is the reactance component,
- `G` is the conductance component,
- `B` is the susceptance component,
- $N$ is the complex transformation factor (e.g. 1 for power lines).

```math
N = \tau e^{j\phi}
```

```math
y_{ser} = \frac{1}{R + jX}
```

```math
y_{shunt} = G + jB
```

The magnitude $\tau$ is the off-nominal tap ratio; the angle $\phi$ is the
phase shift. A pure ratio tap changer moves $\tau$, a pure phase shifter moves
$\phi$, and a combined regulator (German *Schrägregler*) moves both. This is the
standard Y-bus-stamped transformer model and covers every transformer of
practical interest — a real transformer or PST always has a finite series
impedance and is therefore a normal branch stamp.

### Circuit diagram

```
                    y_ser
      x--┓┏---------###----------x
         ||   |             |
         ||   # y_shunt     # y_shunt
         ||   |             |
      x--┛┗----------------------x
         N = tau * exp(j * phi)
```

### Sign and conjugation convention

With series admittance $y$, total shunt $y_{sh}$, and from-side tap $t = \tau
e^{j\phi}$, the stamped entries are:

```math
\begin{aligned}
Y_{ff} &= \frac{y + 0.5\,y_{sh}}{\lvert t \rvert^2} \\
Y_{ft} &= -\frac{y}{\overline{t}} \\
Y_{tf} &= -\frac{y}{t} \\
Y_{tt} &= y + 0.5\,y_{sh}
\end{aligned}
```

This is the Sparlectra sign and conjugation convention. The tap is applied on
the from side; reversing the PST orientation with $\tau = 1$ is electrically
equivalent to flipping the sign of $\phi$, but with an off-nominal ratio,
reversing orientation also moves the ratio tap and is not a pure sign change.

## 2. Y-bus assembly

### Diagonal entries (π-model + shunts)

For a node $i$, the diagonal Y-bus entry is the nodal self-admittance:

```math
Y_{ii} = \sum_{k \in \mathcal{N}(i)} y_{ik} + y_i^{sh}
```

with:
- $y_{ik}$: series admittance contribution of branch $i-k$
- $y_i^{sh}$: explicit shunt admittance at bus $i$

For a π-model branch $i-k$, the local diagonal stamp is:

```math
Y_{ii} \mathrel{+}= y_{ik} + \frac{y_{ik}^{sh}}{2}
```

and the off-diagonal relation is:

```math
Y_{ik} = -y_{ik}
```

Hence (without explicit shunts):

```math
Y_{ii} = -\sum_{k \neq i} Y_{ik}
```

Interpretation: $Y_{ii}$ is the full self-admittance seen at bus $i$; it
combines network coupling (series paths) and local shunts, and represents the
current leaving bus $i$ for $V_i = 1\,\mathrm{pu}$ in nodal form. In typical
grids this supports diagonal dominance and helps the numerical robustness of
NR/SE workflows; the real part is usually non-negative, and the imaginary part
reflects the balance of inductive series effects and capacitive or inductive
shunt terms.

Implementation: branch builders (`addACLine!`, `addPIModelACLine!`,
`addPIModelTrafo!`) stamp series admittance plus half shunt on each side
according to the branch model; explicit shunts are added as nodal shunt terms
when `bus_shunt_model = "admittance"`.

### Bus-shunt modeling modes

Sparlectra supports two representations for real bus shunts imported from
sources such as MATPOWER `Gs`/`Bs` columns:

- `"admittance"` (default): the classical treatment. The bus shunt admittance
  $y_i^{sh} = G_i + jB_i$ is stamped into the Y-bus diagonal as part of
  $Y_{ii}$. Preserves existing numerical behavior.
- `"voltage_dependent_injection"`: the bus shunt is not stamped into Y-bus.
  Its power is evaluated in the nonlinear injection/mismatch path as a local
  voltage-dependent term.

For a bus shunt admittance $y_i^{sh}$ and local voltage magnitude $|V_i|$, the
sign convention is:

```math
S_i^{sh} = |V_i|^2 \overline{y_i^{sh}}
```

A positive conductance contributes positive active power, while the reactive
sign follows the complex conjugate of the shunt admittance. The rectangular
mismatch uses `S_calc - S_spec`, so voltage-dependent injection mode subtracts
$S_i^{sh}$ from the specified net injection. This keeps the equations
equivalent to the admittance model while avoiding double-counting: each bus
shunt is either stamped into Y-bus or represented as a voltage-dependent
injection term, never both. The injection mode is useful when a solver
formulation wants the admittance matrix to contain only branch/network coupling
while keeping shunt effects in the nonlinear injection equations.

## 3. Tap-changer modelling layers

Transformer and PST semantics are richer than a single `ratio + shift` branch.
Sparlectra separates the concerns into explicit layers so that source-format
parsing, tap-changer semantics, equivalent-circuit calculation and the solver
representation stay independent:

```text
Importer (MatpowerIO, DTFImporter, ...)
    | maps source fields -> tap-changer model structs, NO physics formulas
    v
transformer.jl        data types only (no behaviour)
    v
equicircuit.jl        pure functions: model + tap position -> (ratio, shift_deg, x_pu?)
    v
Branch.ratio / Branch.shift_deg / Branch.x_pu
    v
Y-bus stamping / rectangular NR / outer-loop control
```

Guiding rule: every tap/PST formula lives exactly once, in `equicircuit.jl`.
Importers only construct model structs and call the helpers; they contain no
transformer physics formulas of their own.

### Data types (`transformer.jl`)

All tap-changer models share the supertype `AbstractTapChangerModel`.

**Ratio tap changer** — `PowerTransformerTaps` is the ratio-tap variant. It
carries the tap range and step definition (`step`, `lowStep`, `highStep`,
`neutralStep`, `voltageIncrement_kV` / `tapStepPercent`, `tapSign`) plus
nameplate metadata (`neutralU`, `neutralU_ratio`). Its `convention` field makes
the ratio convention explicit; `:neutral_relative` applies the correction
`corr = 1 + (step − neutralStep)·tapStepPercent/100` as a divisor on the
winding ratio. Side information lives on `PowerTransformer.tapSideNumber`.

**Phase tap changer** — `PhaseTapChangerModel` classifies the PST technology
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

**Tap table** — `TapTablePoint` holds one discrete tap row: `step`, `ratio`,
`angle_deg`, and optional `x_pu`. A `PhaseTapChangerModel(kind = :tabular)` is
backed by a non-empty vector of these with strictly ascending, unique steps
(validated, never silently sorted). `lowStep`/`highStep` are derived from the
table when omitted, `neutralStep` must be a step present in the table, and a
tabular model carries no formula parameters (`voltage_step_increment`,
`winding_connection_angle_deg`, `x_min`, `x_max` must be `nothing`) — the table
is the single source of truth. A table overrides formula-based reconstruction
whenever present.

Both tap-changer kinds attach to a winding: `PowerTransformerWinding` has a
`taps::Union{Nothing,PowerTransformerTaps}` slot and a parallel
`phase_taps::Union{Nothing,PhaseTapChangerModel}` slot. This mirrors the CGMES
model, where a tap changer hangs on a transformer end.

!!! note "What the winding connection angle ψ means"
    The winding connection angle `winding_connection_angle_deg` (ψ) is **not** a
    symmetrical-components / sequence angle — Sparlectra works in the positive
    sequence, and ψ lives entirely there. ψ is the angle at which a regulator's
    *additional voltage* is injected relative to the base voltage, i.e. the
    geometry of the regulating vector in the complex voltage plane:

    ```math
    \text{regulating vector} = 1 + f \cdot e^{j\psi}, \qquad
    f = (\text{step} - \text{neutralStep}) \cdot u
    ```

    ψ decides **how a tap move splits between magnitude and phase**:

    - **ψ = 0°** — the added voltage is in phase → pure *longitudinal* (ratio)
      regulator: the regulating vector stays real, `shift_deg` is exactly `0`,
      only the ratio changes.
    - **ψ = 90°** — added voltage in quadrature → *quadrature booster*:
      produces mainly a phase shift.
    - **0° < ψ < 90°** — *Schrägregler*: produces both ratio and phase change
      in the proportion set by ψ.

    From ψ and the tap fraction `f`, `calcPhaseTapAngleRatio` derives the
    effective `ratio` and `shift_deg` that are stamped into the branch (from the
    default from-side convention: `ratio = 1/|v|`, `shift = -arg(v)`). So ψ is
    fully used today — it is the parameter that shapes the effective complex
    tap; it is simply not a per-phase or per-sequence quantity.

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

Symmetrical phase-shifter (pure quadrature-type angle regulation):

```julia
pst_sym = PhaseTapChangerModel(
  kind = :symmetrical,
  step = 3, lowStep = -10, highStep = 10, neutralStep = 0,
  voltage_step_increment = 0.012,     # per-step, pu of rated voltage
)
```

Asymmetrical phase-shifter / Schrägregler (ψ ≠ 0):

```julia
pst_skew = PhaseTapChangerModel(
  kind = :asymmetrical,
  step = 7, lowStep = -13, highStep = 13, neutralStep = 0,
  voltage_step_increment = 0.18 / 13, # e.g. 18 % over 13 steps
  winding_connection_angle_deg = 60.0,
)
# quadrature booster is the same with winding_connection_angle_deg = 90.0
```

Tabular phase-shifter (table overrides formulas; carries no formula params):

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
)
# the equivalent-circuit helpers resolve model + step -> ratio/shift for the branch
```

### Behaviour (`equicircuit.jl`)

The equivalent-circuit helpers turn a model plus a tap position into the
branch quantities. They are pure functions, one per formula family:

| Function | Purpose |
|---|---|
| `calcRatioTapCorrection(taps; step)` | ratio-tap multiplicative correction `1 + (step − neutralStep)·tapStepPercent/100` |
| `calcRatioTapRange(taps)` | `(tap_min, tap_max, tap_step)` in ratio terms |
| `calcPhaseTapFraction(m; step)` | shared tap fraction `f = (step − neutralStep)·voltage_step_increment` |
| `calcPhaseTapAngleRatio(m; step)` | `(effective_ratio, effective_shift_deg, regulating_vector)` |
| `calcPhaseTapReactance(m, α)` | tap-angle-dependent reactance `X(α)` interpolated between `x_min`/`x_max` |
| `calcPhaseTapTable(m; step)` | exact lookup of a tabular row |

For `calcPhaseTapAngleRatio`, a `:symmetrical` changer computes
`α = 2·atand(f/2)` with magnitude always `1.0`; an `:asymmetrical` changer maps
the regulating vector `1 + f·e^{jψ}` (with winding connection angle ψ) through
the low-level primitive `calcSkewAngleTap`, of which the quadrature booster
(ψ = 90°) is a special case. A `:tabular` model resolves ratio and angle by
lookup and reconstructs the regulating vector from the stored degrees. The
`calcVKDependence` spline over tap tables is the precedent for tabular
interpolation should a smooth characteristic ever be required.

### Reactance dependence X(α)

For PSTs whose series reactance varies across the tap range,
`calcPhaseTapReactance` evaluates `X(α)` by interpolating between the endpoint
reactances `x_min = X(0)` and `x_max = X(αmax)` per technology, or returns the
tabular `x_pu` of the row for a `:tabular` model. The reactance helper is
available and independently usable; whether a solved operating point tracks
`X(α)` as taps move depends on how the branch reactance is fed into the Y-bus
between control iterations.

### Importer mapping

- **DTF**: builds a `PhaseTapChangerModel(kind = :asymmetrical,
  winding_connection_angle_deg = skew, ...)` and calls `calcPhaseTapAngleRatio`
  to derive the branch `ratio`/`shift`. The skew-angle physics lives in the
  equivalent-circuit layer, not in the parser. The pure-longitudinal case
  (ψ = 0) keeps the shift at exactly `0.0`.
- **MATPOWER**: keeps the direct `TAP`/`SHIFT` path — this is the CGMES
  "General Case" (raw values), and no model struct is required. Branch `SHIFT`
  is interpreted as the phase angle $\phi$ on the from side by default
  (`matpower_shift_unit = "deg"`, `matpower_shift_sign = 1`). Some PEGASE-style
  cases carry small radian-like phase-shifter values; `matpower_import.jl` can
  set `matpower_shift_unit = "rad"` and `matpower_shift_sign = -1` to test or
  apply that convention. MATPOWER branch `TAP` is used as stored by default
  (`matpower_ratio = "normal"`); set `matpower_ratio = "reciprocal"` for input
  files whose off-nominal transformer ratios must be inverted before import.

### Existing tap-impedance correction

Independent of the typed PST models, Sparlectra also offers an imported-case
tap-changer reactance treatment selected by `transformer.tap_changer_model`:
`ideal` (default) keeps the tap changer free of series-impedance feedback,
while `impedance_correction` re-refers transformer R/X through the tapped
winding via `|1 + f·e^{jφ}|²`. It applies to all transformers of an imported
case (both MATPOWER and DTF importers) and is implemented centrally in
`calcTapCorrectedRX` / `calcTapImpedanceCorrectionFactor`.

### Three-winding transformers

A three-winding transformer is modelled as a star (T) equivalent with an
auxiliary star-point bus: each of the three windings becomes its own
`PowerTransformerWinding` and is stamped as a separate branch to the AUX bus
(`create3WTWindings!`, MVA method). Because every winding is a full
`PowerTransformerWinding`, each already carries its own `taps` and `phase_taps`
slots — so a phase-shifting winding (e.g. a three-winding *Schrägregler*) is
represented by placing the regulating vector on the branch from that one
winding to the star point, leaving the other two windings unaffected. The
positive-sequence stamping and the ψ interpretation are exactly the same as for
a two-winding device; the star point simply gives each winding its own branch
to regulate.

`create3WTWindings!` accepts an optional `phase_tap_side` (winding index
`1..3`, `0` = none) and `phase_taps::PhaseTapChangerModel` pair to attach a
PST model to one winding — the same 1-based convention as `tap_side`, and
`phase_tap_side` may equal `tap_side` when a winding carries both a ratio tap
and a phase tap (Schrägregler):

```julia
psc = PhaseTapChangerModel(kind = :asymmetrical, step = 0, lowStep = -8, highStep = 8, neutralStep = 0, winding_connection_angle_deg = 60.0)
w1, w2, w3 = create3WTWindings!(u_kV = [110.0, 20.0, 10.0], sn_MVA = [100.0, 80.0, 20.0], addEx_Side = [tmp1, tmp2, tmp3], sh_deg = [0.0, 0.0, 0.0], tap_side = 1, tap = tapSettings, phase_tap_side = 2, phase_taps = psc)
```

Resolving `w2.phase_taps` into an effective ratio/shift on the AUX-bus branch,
and addressing a single 3WT winding from the outer-loop
`PowerTransformerControl` framework, are not implemented yet — see the
analysis report at `docs/dev/3wt_phase_tap_controller_addressing.md` for the
current gaps.

## 4. Transformer control (outer loop)

Sparlectra regulates transformers within the branch PI model using the complex
tap `t = τ·e^{jφ}` and **without auxiliary nodes**: `τ` for voltage control,
`φ` for active-power-flow (PST) control, both together for a *Schrägregler*.

### Numerical method

Tap control is an outer loop around the power flow:

1. Solve PF with current taps
2. Evaluate the control error
3. Update the tap(s) (continuous or discrete)
4. Re-run PF
5. Stop on convergence, limits, or iteration cap

No augmentation of the Newton system is performed. Taps are supervisory control
updates, not algebraic unknowns. This keeps the Jacobian/state vector
untouched, allows simpler solver backends, centralizes control logic
(deadbands, limits, discrete steps), and gives deterministic post-processing
between iterations.

### Discrete tap behaviour

```text
tap_ratio_new       = clamp(tap_ratio ± tap_step,             tap_min,       tap_max)
phase_shift_deg_new = clamp(phase_shift_deg ± phase_step_deg, phase_min_deg, phase_max_deg)
```

### Phase-shift control direction (practical probe)

The sign of a phase-shifter control action should not be hard-coded; probe it
on the active model:

1. Compute `P_ab(phi = 0 deg)`
2. Compute `P_ab(phi = +5 deg)`
3. Evaluate `Delta_P_ab = P_ab(5 deg) − P_ab(0 deg)`
4. In control: if `P_ab < target`, move `phi` in the direction that increases
   `P_ab`; otherwise move it the opposite way.

See `examples/example_transformer_phase_shift_control.jl`.

### Controller framework

Transformer tap/phase control is implemented as an `AbstractOuterController`.
Controllers are collected by `collect_outer_controllers(net)` (and, for
tap-specific resolution, `_tap_controllers(net)` from
`PowerTransformerWinding.controls`), deduplicated by identity, and used
consistently for execution and reporting. When at least one controller is
present, `run_sparlectra` calls `run_control!` for outer-loop orchestration.
The controllers adjust `ratio`/`shift` directly within their limits; they do
not depend on the CGMES tap-changer classification.

### API

```julia
addPowerTransformerControl!(net;
  trafo = "1",
  mode = :voltage,
  target_bus = "B5",
  target_vm_pu = 1.01,
  control_ratio = true,
  control_phase = false,
  is_discrete = true)
```

Active-power control on a branch:

```julia
addPowerTransformerControl!(net;
  trafo = "1",
  mode = :branch_active_power,
  target_branch = ("B1", "B2"),
  p_target_mw = 250.0,
  control_ratio = false,
  control_phase = true,
  is_discrete = true)
```

Combined voltage + active-power (*Schrägregler*):

```julia
addPowerTransformerControl!(net;
  trafo = "1",
  mode = :voltage_and_branch_active_power,
  target_bus = "B5",
  target_vm_pu = 1.01,
  target_branch = ("B1", "B2"),
  p_target_mw = 250.0,
  control_ratio = true,
  control_phase = true,
  is_discrete = true)
```

Advanced direct use and result inspection:

```julia
result = run_control!(
  net;
  pf_config = powerflow_config(),
  control_config = control_config(),
)

result = latest_control_result(net)
println(result.status)            # outer control-loop terminal state
println(result.outer_iterations)
println(result.powerflow_solves)
println(result.controllers)
println(result.trace)             # machine-readable, no console parsing
```

`result.status` is the outer control-loop terminal state, separate from
numerical PF success/failure.

### Inline controller definition

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

### Remote voltage control — scope

Sparlectra supports basic remote voltage control: a `target_bus` measurement,
one transformer tap as actuator, and a `target_vm_pu ± deadband` objective —
i.e. single-controller remote regulation. A complete implementation as used in
real grid control systems or CGMES-based coordination additionally requires
multiple transformers controlling the same remote bus, coordination between
controllers (participation factors / priority rules), limit handling with
redistribution when a transformer hits its tap limit, anti-hunting mechanisms,
and deterministic group convergence. Currently controllers operate
independently with no grouping or shared objective; coordinated multi-actuator
remote voltage control is not yet implemented.

### Limits / scope

- no auxiliary transformer nodes
- no coupling of tap variables into the Newton iteration
- no coordinated multi-transformer control

## 5. CGMES / ENTSO-E mapping

The typed tap-changer models are aligned with the CGMES data model, so that
CIM-based exchange maps onto Sparlectra with minimal reinterpretation:

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

The CGMES guidance recommends exchanging tabular tap data where available
instead of recalculating parameters from technology formulas. Sparlectra
honours this: a `:tabular` model overrides the formula path, and the formula
kinds (`:symmetrical`, `:asymmetrical`) are used only when no table is
provided. MATPOWER's raw `TAP`/`SHIFT` corresponds to the CGMES "General Case".

## 6. Related documents

- ENTSO-E, *Phase Shift Transformers Modelling*, CGMES v2.4, 28 May 2014 — the
  reference for PST technology classification (symmetrical / asymmetrical),
  the tap-angle formulas, the reactance-versus-angle characteristics, and the
  recommendation to exchange tabular data.
- IEC 61970-301 (CIM base) and the CGMES profiles — the class definitions
  (`RatioTapChanger`, `PhaseTapChanger*`, `TapChangerTablePoint`) referenced in
  the mapping table above.
- MATPOWER case format documentation — the `TAP` / `SHIFT` branch columns and
  the ratio/shift conventions handled at import.
- Sparlectra examples: `examples/tap_control_demo_grid.jl` (OLTC voltage, PST
  active-power, and Schrägregler combined control in one network) and
  `examples/example_transformer_phase_shift_control.jl` (phase-shift direction
  probe).
