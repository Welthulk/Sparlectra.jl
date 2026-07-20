# Task: Issue #261 Stage 4 — Tabular PST model (TapTablePoint, kind = :tabular)

Repository: Sparlectra.jl (branch off `main`, requires Stage 3 merged:
`PhaseTapChangerModel`, `calcPhaseTapFraction`, `calcPhaseTapAngleRatio`,
`calcPhaseTapReactance`).
Reference: `docs/src/transformer_pst_architecture.md`.
Note: this task was written against the Stage-3 task spec; verify actual
signatures/field order on-site before editing.

Principle (Issue #261): tabular data OVERRIDES formula-based reconstruction
whenever a table is present. Real PSTs are best exchanged as tables
(cim:PhaseTapChangerTabular / cim:PhaseTapChangerTablePoint).

## Step 1 — `src/transformer.jl`

1. New immutable struct (data only), placed before `PhaseTapChangerModel`:
   ```julia
   struct TapTablePoint
     step::Int
     ratio::Float64
     angle_deg::Float64
     x_pu::Union{Nothing,Float64}
   end
   ```
   Keyword constructor with `x_pu = nothing`.
2. Extend `PhaseTapChangerModel`: new field
   `table::Union{Nothing,Vector{TapTablePoint}}` as LAST field, keyword
   default `nothing`.
3. Update constructor validation:
   - REMOVE the Stage-3 rejection of `:tabular` (and adapt the Stage-3 test
     that asserted this rejection — this is the ONLY permitted edit to
     existing tests).
   - `kind == :tabular` requires: `table !== nothing`, non-empty, strictly
     ascending unique `step` values (sort a copy or throw — decide: throw
     `ArgumentError` on unsorted/duplicate steps, do not silently sort).
   - For `kind == :tabular`: `lowStep`/`highStep`/`neutralStep` must be
     consistent with the table: `lowStep == first step`,
     `highStep == last step`, `neutralStep` must exist as a table step.
     Convenience: if the caller passes the table, allow
     `lowStep`/`highStep` to be derived automatically when omitted
     (keyword defaults computed from the table); `neutralStep` stays required.
   - `kind != :tabular` with `table !== nothing` → ArgumentError (no silent
     mixing of formula and table models).
   - `voltage_step_increment` / `winding_connection_angle_deg` /
     `x_min` / `x_max` must be `nothing` for `:tabular` (ArgumentError
     otherwise) — the table is the single source of truth.
4. `Base.show`: for `:tabular` print number of table points and step range
   instead of increment fields.

## Step 2 — `src/equicircuit.jl`

1. New pure lookup:
   ```julia
   """
       calcPhaseTapTable(m::PhaseTapChangerModel; step::Int = m.step)
         -> (effective_ratio, effective_shift_deg, x_pu::Union{Nothing,Float64})

   Exact lookup of the table row for `step`. Throws ArgumentError if `step`
   has no table row or if `m.table === nothing`.
   ```
   Implementation detail: build the lookup via `findfirst` over the (already
   validated, sorted) vector — no Dict caching in this stage.
2. Integrate into the existing resolvers (table overrides formulas):
   - `calcPhaseTapAngleRatio`: `kind == :tabular` → delegate to
     `calcPhaseTapTable`, derive the regulating vector from
     `(ratio, angle_deg)` consistently with the model's `convention`
     (document the derivation in the docstring; must be the exact inverse of
     how `calcSkewAngleTap` maps regulating vector → (ratio, shift) for
     `:reciprocal_from_side`).
   - `calcPhaseTapFraction`: `kind == :tabular` → ArgumentError with a clear
     message (there is no linear fraction for tabular models).
   - `calcPhaseTapReactance`: `kind == :tabular` → return the row's `x_pu`
     (may be `nothing`); the `alpha_deg` argument is ignored for tabular —
     document this.

No interpolation between steps in this stage: taps are discrete; interpolation
(for continuous outer-loop control) is deferred until the X(α)/outer-loop
stage. State this in the docstring.

## Step 3 — Exports

`src/Sparlectra.jl`: export `TapTablePoint`, `calcPhaseTapTable`.

## Step 4 — Tests

New file `test/test_phase_tap_table.jl`, registered in `test/runtests.jl`
(same `controls` group) and `docs/src/tests.md`:

1. Constructor validation: empty table, duplicate steps, unsorted steps,
   `neutralStep` not in table, `table` on non-tabular kind, increment/ψ/x
   fields set on `:tabular` — each throws ArgumentError.
2. Auto-derivation of `lowStep`/`highStep` from the table.
3. `calcPhaseTapTable`: exact lookup incl. `x_pu === nothing` passthrough;
   missing step throws.
4. Override semantics: a `:tabular` model routed through
   `calcPhaseTapAngleRatio` returns exactly the table values (ratio,
   angle) — build a table from values produced by the Stage-3 asymmetrical
   formula and assert the tabular path reproduces the formula path
   bit-identically (round-trip test; also validates the regulating-vector
   derivation).
5. `calcPhaseTapReactance` on `:tabular`: returns row `x_pu`, ignores
   `alpha_deg`.
6. `calcPhaseTapFraction` on `:tabular` throws.
7. Adapt the single Stage-3 test that asserted `:tabular` rejection
   (now asserts successful construction instead).
8. Full suite green otherwise without edits.

## Step 5 — Documentation

1. Append "Stage 4 implemented" note to
   `docs/src/transformer_pst_architecture.md`: TapTablePoint fields with CIM
   mapping (cim:PhaseTapChangerTablePoint.{step,ratio,angle,x}), override
   rule, no-interpolation decision and its rationale.
2. Changelog entry under the current unreleased version, referencing
   Issue #261 (Stage 4).
3. Watch `@ref` targets in docstrings (Documenter fails on links to
   undocumented symbols — known from Stage 2).

## Acceptance criteria

- [ ] Full test suite green; only the one Stage-3 rejection test modified.
- [ ] Round-trip test proves tabular path == formula path for
      formula-generated tables (bit-identical).
- [ ] `:tabular` models carry no formula parameters (constructor enforces).
- [ ] No importer changes in this stage (grep: `TapTablePoint` appears in
      `src/transformer.jl`, `src/equicircuit.jl`, `src/Sparlectra.jl` only).
- [ ] `docs/make.jl` builds.

## Out of scope (do NOT implement)

- Importer support for tables (no source format delivers them yet; DTF and
  MATPOWER untouched).
- Interpolation between table steps.
- Wiring table `x_pu` into `Branch.x_pu` or the outer control loop.
- Ideal transformers / supernode preprocessing.
- Any change to Y-bus stamping, Newton solver, Q-limit logic, control
  framework.
