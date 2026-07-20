# Task: Three-winding transformer phase-tap support (Schrägregler on a 3WT winding)

Repository: Sparlectra.jl (branch off `main`, typed tap-changer models merged:
`AbstractTapChangerModel`, `PhaseTapChangerModel`,
`PowerTransformerWinding.phase_taps`).
Constraint: **zero behaviour change for every existing call site and test.**
This task was written against the pre-merge source for `create3WTWindings!`
(the 3WT code was not touched by the tap-model work) — verify actual
signatures/field order on-site before editing.

## Background

A 3WT is modelled as a star (T) equivalent with an AUX star-point bus; each
winding is a full `PowerTransformerWinding` stamped as its own branch to the
AUX bus (`create3WTWindings!` in `src/transformer.jl`, MVA method). Since
every winding already has a `phase_taps::Union{Nothing,PhaseTapChangerModel}`
slot, the data model needs no change. What is missing:

1. `create3WTWindings!` accepts only `tap::PowerTransformerTaps` (ratio, one
   side via `tap_side`) plus a static `sh_deg` vector — there is no way to
   attach a `PhaseTapChangerModel` to one of the three windings.
2. It is unverified whether the outer-loop `PowerTransformerControl` can
   address a *single winding branch* of a 3WT as its actuator.

## Part 1 — Extend `create3WTWindings!` (implementation)

1. Verify the current keyword signature on-site
   (`src/transformer.jl`, function `create3WTWindings!`; pre-merge it was:
   `u_kV, sn_MVA, addEx_Side, sh_deg, tap_side, tap::PowerTransformerTaps`).
2. Add two optional keywords, defaults preserving today's behaviour exactly:
   - `phase_tap_side::Int = 0` — winding index `1..3` carrying the phase-tap
     model; `0` = none.
   - `phase_taps::Union{Nothing,PhaseTapChangerModel} = nothing`.
   Validation (ArgumentError with clear message):
   - `phase_tap_side == 0` XOR `phase_taps === nothing` mismatches
     (`phase_tap_side != 0` requires a model and vice versa);
   - `phase_tap_side ∈ 0:3`;
   - `phase_tap_side` may equal `tap_side` (ratio + phase on the same winding
     is a valid Schrägregler configuration) — do NOT forbid it.
3. In the winding-construction loop, attach the model to the selected
   winding's `phase_taps` field. NOTE the pre-existing off-by-one style in
   that loop: the ratio tap uses `(side - 1 == tap_side)`. Do NOT copy this
   blindly — determine on-site whether `tap_side` is 0-based or 1-based in
   the existing convention (check call sites and docstring example
   `tap_side = 1`), then make `phase_tap_side` follow the SAME convention as
   `tap_side` and document both in the docstring. If the existing convention
   turns out to be inconsistent between docstring and code, report it —
   do not silently "fix" the ratio-tap behaviour.
4. Do NOT derive or modify `sh_deg` from the model in this task. The static
   `sh_deg` and the phase-tap model coexist; resolving the model to an
   effective `ratio`/`shift` for the AUX-bus branch is intentionally out of
   scope here (see Part 2 findings first). This task only makes the model
   available on the winding.
5. Update the docstring (fields, conventions, one example with a
   `PhaseTapChangerModel` on side 2).

## Part 2 — Controller addressing of a 3WT winding (analysis only, no code)

Investigate and report — do not implement:

1. How does `addPowerTransformerControl!` / `PowerTransformerControl` resolve
   `trafo = "..."` to a branch? (`src/tap_control.jl`,
   `addPowerTransformerControl!` and the branch lookup it performs.)
2. For a 3WT (three branches to the AUX bus): can the current lookup select
   exactly one winding branch? If it matches by transformer name only, which
   branch does it pick, and is that deterministic?
3. Where do `tap_min/tap_max/tap_step` and `phase_min_deg/phase_max_deg/
   phase_step_deg` on the branch come from for 3WT branches today
   (`src/branch.jl` derivation from `w.taps` — does the 3WT construction path
   populate them at all)?
4. Deliver a short written report (markdown, `docs/dev/` or as PR comment)
   answering: what exactly is missing for "outer-loop phase control of one 3WT
   winding", with file/function references — but change no control code.

## Part 3 — Tests

New test file `test/test_3wt_phase_taps.jl`, registered following the existing
pattern in `test/runtests.jl` and `docs/src/tests.md`:

1. Existing behaviour: `create3WTWindings!` called exactly as in the current
   docstring example (no new keywords) produces windings with
   `phase_taps === nothing` on all three sides and identical values to before
   (construct once, compare fields against hand-derived expectations or a
   pre-change snapshot within the test).
2. Attachment: with `phase_tap_side` selecting side 2 and an `:asymmetrical`
   model (ψ = 60°), exactly winding 2 has `phase_taps !== nothing`, sides 1
   and 3 have `nothing`, and the stored model round-trips (kind, ψ, steps).
3. Ratio + phase on the same side (`tap_side == phase_tap_side`): both `taps`
   and `phase_taps` set on that winding.
4. Validation: `phase_tap_side = 2` without model throws; model without side
   throws; `phase_tap_side = 4` throws.
5. Full suite green; no existing test modified.

## Part 4 — Documentation

1. Extend the "Three-winding transformers" subsection of the merged branch/
   transformer model page (`docs/src/branchmodel.md`) with one sentence and
   the new constructor example.
2. Changelog entry under the current unreleased version (new keywords on
   `create3WTWindings!`; controller addressing analysed, reference the report).
3. Watch `@ref` targets in docstrings (Documenter fails on links to
   undocumented symbols).

## Acceptance criteria

- [ ] All existing tests pass unmodified; default-argument calls of
      `create3WTWindings!` are behaviour-identical.
- [ ] New keywords validated as specified; conventions documented and aligned
      with the existing `tap_side` convention.
- [ ] Part-2 report exists and answers the four questions with code
      references.
- [ ] New tests green; `docs/make.jl` builds.

## Out of scope (do NOT implement)

- Resolving a 3WT winding's `PhaseTapChangerModel` into effective
  `ratio`/`shift` on the AUX-bus branch (depends on Part-2 findings).
- Any change to `PowerTransformerControl`, the outer-loop framework, or
  branch tap-limit derivation.
- DTF/MATPOWER importer changes (no source format delivers 3WT phase-tap
  data today).
- X(α) coupling, ideal couplings, solver/Y-bus changes.
