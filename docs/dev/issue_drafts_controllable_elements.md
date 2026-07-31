# Issue drafts: controllable-element increments after #227 points 1+2

Status: drafts, 2026-07-31. Points 1 (generic controllable-element
reporting, `controllableElements`) and 2 (SVC-style `ShuntVoltageControl`)
are implemented; the sections below are ready-to-post issue texts for the
remaining increments from
`docs/dev/issue227_controllable_elements_analysis.md`. Each is independent.

---

## Draft 1: TCSC-like series-reactance controller

**Title:** Add a TCSC-like series-reactance controller for line branches

### Background

The outer control framework now supports a series-reactance actuator: the
PST X(α) coupling updates `br.x_pu` on accepted tap moves and the next
outer-loop solve re-stamps the Y-bus (see the tap-dependent reactance
section in `branchmodel.md`). A TCSC (thyristor-controlled series
capacitor) is exactly this actuator on a *line* branch: a variable series
reactance within `[x_min, x_max]`, steered onto a branch-flow target.

### Proposal

- `SeriesReactanceControl <: AbstractOuterController` with
  `addSeriesReactanceControl!(net; branch, p_target_mw, x_min_pu, x_max_pu,
  deadband_p_mw, …)`.
- Actuator: `br.x_pu` of the controlled line branch (continuous), clamped to
  the range; the excluded resonance band of a physical TCSC can be a later
  refinement (two disjoint ranges).
- Evaluate/propose via secant on (x, P) like the machine/SVC controllers;
  the phase-probe pattern (flow refresh around probe solves) applies.
- `controllableElements` descriptor: actuator `:series_x_pu`, quantity
  `:branch_active_power`.
- At-limit semantics: clamped x, honest `at_limit` — the branch then behaves
  as a fixed compensated line.

### Acceptance

- A loop network where the TCSC steers the parallel-path flow onto a target;
  honest `at_limit` when the target is out of range; a branch without the
  controller is bit-identical; element row present; example + tests + docs.

---

## Draft 2: STATCOM variant with current-based reactive limit

**Title:** STATCOM mode for machine voltage control: current-based Q limit

### Background

`MachineVoltageControl` holds a remote (or local) voltage with the machine's
reactive power inside constant `[qmin, qmax]` — the standard machine limit.
A STATCOM is current-limited: at its limit it injects maximum reactive
*current*, so the available Q scales with the terminal voltage,
`Q_lim = V · I_max`, instead of staying constant.

### Proposal

- Add an optional current-limit mode to the machine controller (or a thin
  `StatcomVoltageControl` wrapper): kwargs `i_max_ka` (or `s_max_mva` at
  1.0 p.u.) replacing fixed `qmin/qmax`.
- Per outer iteration, re-evaluate the bounds from the current bus voltage
  before clamping the secant step; `at_limit` then reflects the
  voltage-dependent bound.
- Element descriptor: device "STATCOM (VSC)"; actuator stays
  `:machine_q_mvar` with live bounds.

### Acceptance

- In-range behavior identical to today; at the limit the delivered Q tracks
  `V · I_max` across outer iterations (test with a depressed-voltage case);
  constant-limit mode unchanged bit-for-bit; example + tests + docs.

---

## Draft 3: Back-to-back HVDC pairing controller

**Title:** Pairing controller for back-to-back HVDC converter injections

### Background

CGMES multi-area imports map DC border crossings as two per-side converter
injections (Stage-0: fixed PCC operating points, no coupling). A
back-to-back link is a *controlled* coupling: transfer `P` with
`P1 + P2 + P_loss = 0` and per-terminal Q or V control.

### Proposal

- `HvdcPairControl <: AbstractOuterController` referencing the two
  converter prosumers: kwargs `p_transfer_mw` (signed), `loss_mw` (fixed or
  fraction), per-terminal `q_mvar` or `vset_pu`.
- Apply step keeps the pairing invariant exactly (`P2 = −P1 − P_loss`);
  terminal voltage control reuses the machine-controller secant per side.
- Importer wiring later: recognize the split-border converter pairs
  (`_detectNonCancellingBoundarySides!` already identifies them) and offer
  opt-in controller attachment.

### Acceptance

- A two-area net coupled only by the pair: transfer target reached, area
  balances consistent, island handling unchanged (the areas stay separate
  AC islands); honest `at_limit` on converter ratings; element rows for the
  pair; example + tests + docs.

---

## Draft 4: Control sensitivities from the solved Jacobian

**Title:** Control sensitivities: dx/du for registered controllers

### Background

For a solved power flow, `dx/du = -J⁻¹ · ∂r/∂u` estimates how a control
action changes voltages, angles, and flows. The generic element records
(`controllableElements`) now name every registered actuator `u`; the solver
factorizes `J` anyway. This replaces probe-style finite differences (extra
solves) with one linear solve per actuator.

### Proposal

- Solver-side hook exposing the last factorized Jacobian (or a re-factorize
  API) after convergence.
- Per actuator kind, the residual derivative `∂r/∂u`: shunt susceptance
  (diagonal Q term), machine Q (single bus injection), tap ratio / phase
  shift / series x (branch stamp derivatives).
- Public API per the issue sketch: `compute_control_sensitivity(net, ctrl)`
  → per-bus dVm/du (+ optional branch-flow rows);
  `compute_influence_footprint(net, ctrl; threshold)` → the buses where the
  sensitivity exceeds a threshold; `rank_controllers_for_target(net,
  target_bus)`.
- Validation: compare against finite-difference probes on small nets.

### Acceptance

- Sensitivities match finite differences within tolerance for all four
  actuator kinds; footprint/ranking demonstrated on a mid-size case; no
  change to solve results; docs incl. the influence-area use case
  (tap-change influence areas from the original issue motivation).
