# Task: Validate DTF importer migration (Issue #261 Stages 2–4)

Read-only validation. **Do NOT change any source or test file.** The goal is to
confirm that the DTF importer actually uses the new typed tap-changer models
and that no "programmed around it" duplication was reintroduced. Report
findings; only fix if a check fails (and then only the minimal fix, see final
section).

Repository: Sparlectra.jl, current feature branch (Stages 2–4 merged).

## Part A — Static checks (grep-level acceptance criteria)

Run each and report the actual output. Expected result stated per check.

1. Importer constructs the new model:
   ```bash
   grep -n "PhaseTapChangerModel" src/DTFImporter.jl
   ```
   Expected: at least one hit at the tap-model construction site inside
   `_dtf_effective_transformer_tap`.

2. No direct low-level primitive call left in the importer:
   ```bash
   grep -n "calcSkewAngleTap" src/DTFImporter.jl
   ```
   Expected: **no hits** (the importer must go through
   `calcPhaseTapAngleRatio`, not call `calcSkewAngleTap` directly).

3. Tap fraction formula exists exactly once, in the importer field mapping:
   ```bash
   grep -rn "longitudinal_range_percent / 100" src/
   ```
   Expected: exactly one hit, in `src/DTFImporter.jl`.

4. Ratio-tap correction formula not duplicated (Stage 2):
   ```bash
   grep -rn "neutralStep) \* " src/
   grep -rn "tapStepPercent / 100" src/
   ```
   Expected: the correction `1 + (step - neutralStep) * tapStepPercent/100`
   appears only inside `calcRatioTapCorrection` in `src/equicircuit.jl`
   (docstring mentions are fine; count real code occurrences and report).

5. Importer imports the new helpers:
   ```bash
   grep -n "using ..Sparlectra:" src/DTFImporter.jl
   ```
   Expected: import list includes `PhaseTapChangerModel` and
   `calcPhaseTapAngleRatio` (and `calcPhaseTapFraction` if used); report whether
   `calcSkewAngleTap` is still imported and whether it is still referenced
   anywhere in the file (it should not be, per check 2).

6. Exports present:
   ```bash
   grep -n "AbstractTapChangerModel\|PhaseTapChangerModel\|TapTablePoint\|calcRatioTapCorrection\|calcRatioTapRange\|calcPhaseTapFraction\|calcPhaseTapAngleRatio\|calcPhaseTapReactance\|calcPhaseTapTable" src/Sparlectra.jl
   ```
   Expected: all of these are exported.

Produce a short table: check | command | expected | actual | pass/fail.

## Part B — Runtime check: importer really builds a PhaseTapChangerModel

Write a **throwaway** REPL snippet (do not commit; run it inline via
`julia --project=. -e '...'` or a temp file you delete afterwards) that:

1. Loads a DTF case that contains a phase-shifting / skew-angle transformer.
   Locate a suitable fixture used by the DTF tests:
   ```bash
   grep -rln "added_voltage_angle_deg\|longitudinal_range_percent\|skew" test/ | head
   grep -rln "FOR002\|\.dtf" test/ | head
   ```
   Reuse whatever loader the DTF tests use (find it in the DTF test file's
   includes / helper functions). Do not invent a new import path.

2. After import, confirm the model is actually attached to the transformer
   winding (the Stage-3 field `PowerTransformerWinding.phase_taps`):
   - iterate the imported net's transformers,
   - for the PST/skew transformer, check that the relevant winding's
     `phase_taps` is a `PhaseTapChangerModel` (not `nothing`),
   - print `kind`, `winding_connection_angle_deg`, `step`, `highStep`,
     `convention`.

   NOTE: if the current DTF migration only uses `PhaseTapChangerModel`
   transiently inside `_dtf_effective_transformer_tap` and does NOT persist it
   on the winding, `phase_taps` will be `nothing`. Report this clearly — it is
   an important finding (see Part D), not necessarily a bug: Stage 3 required
   the formula to move into the model, but persisting the model on the winding
   may or may not have been implemented. State exactly which is the case.

3. Independently verify the resolver path is exercised: call
   `calcPhaseTapAngleRatio` on a constructed `:asymmetrical`
   `PhaseTapChangerModel` mirroring the imported case's parameters and confirm
   it returns the same `(ratio, shift_deg)` the imported branch ended up with
   (`br.ratio`, `br.phase_shift_deg` / `br.angle`). Use `isapprox` with a tight
   tolerance; report the max deviation.

## Part C — DTF suite still green

Run the DTF-relevant tests (the profile that includes `dtf_extended` / FOR002).
Determine the correct invocation from `test/runtests.jl` (how profiles/groups
are selected — do not guess flags):

```bash
grep -n "dtf\|FOR002\|profile\|ARGS\|group" test/runtests.jl | head -30
```

Then run that profile and report: total tests, pass/fail count, and any
`@test_skip`/`@test_broken`. Expected per prior summary: all green except one
pre-existing `@test_skip`.

## Part D — Report

Produce a concise findings report answering exactly these questions:

1. Does the DTF importer construct a `PhaseTapChangerModel`? (yes/no + line)
2. Is `calcSkewAngleTap` still called directly by the importer? (must be no)
3. Is the tap-fraction / ratio-correction math free of duplication? (yes/no,
   with the grep counts)
4. Is the `PhaseTapChangerModel` **persisted** on the winding
   (`phase_taps !== nothing`) after import, or only used transiently during
   resolution? State which.
5. Does the DTF test suite pass? (numbers)
6. Do the imported branch results match a direct `calcPhaseTapAngleRatio`
   evaluation? (max deviation)

## If a check fails

- Static duplication reintroduced (Part A 3/4 > 1 real hit) or importer still
  calls `calcSkewAngleTap` directly (A2): report the offending lines; propose
  the minimal fix but do NOT apply it without confirmation.
- If Part B.2 shows `phase_taps === nothing` (model not persisted): do NOT
  "fix" it silently. This is a scoping question — persisting the model on the
  winding is needed for the future X(α) coupling issue but was not necessarily
  in Stage 3's scope. Report it as an open item and stop.

Do not commit anything. Delete any temporary scripts you created.
