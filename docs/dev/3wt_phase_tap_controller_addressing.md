# 3WT winding addressing for outer-loop transformer control (analysis)

Issue #261 follow-up: can `PowerTransformerControl` (the outer control loop,
`src/tap_control.jl`) address a *single winding branch* of a genuine
three-winding transformer built with `create3WTWindings!` (the MVA-method,
single-`PowerTransformer`-with-`side1/side2/side3` model, `src/transformer.jl`)?
Analysis only — no control code was changed for this report.

## 1. How `addPowerTransformerControl!` resolves `trafo = "..."` to a branch

`addPowerTransformerControl!` (`src/tap_control.jl:289-353`) resolves the
`trafo::String` argument in two steps:

1. `br = _find_trafo_branch(net, trafo)` (`src/tap_control.jl:307`) calls the
   helper at `src/tap_control.jl:66-73`:
   ```julia
   function _find_trafo_branch(net::Net, name::String)::Branch
     for br in net.branchVec
       if br.ratio != 0.0 && (br.comp.cName == name || br.comp.cID == name || string(br.branchIdx) == name)
         return br
       end
     end
     error("PowerTransformerControl: transformer $(name) not found. ...")
   end
   ```
   This returns the **first** `Branch` in `net.branchVec` whose own component
   name/id, or branch index, matches `name`. It operates purely on `Branch`
   identity — there is no "side"/"winding" concept here at all.
2. `trafo_obj = findfirst(t -> t.comp.cFrom_bus == br.comp.cFrom_bus && t.comp.cTo_bus == br.comp.cTo_bus, net.trafos)`
   (`src/tap_control.jl:308`) then looks up the *parent* `PowerTransformer`
   object by matching **the transformer's own single `comp.cFrom_bus`/`cTo_bus`**
   against the *resolved branch's* `cFrom_bus`/`cTo_bus`.
3. The new controller is always pushed to `net.trafos[trafo_obj].side1.controls`
   (`src/tap_control.jl:351`) — hard-coded to `side1`, regardless of which
   `Branch` was actually resolved in step 1.

## 2. Can the lookup select exactly one winding branch of a 3WT?

**Branch identity (step 1) is not the blocking problem.** For a 3WT built via
`create3WTWindings!`, each winding is stamped as its own `Branch` from the AUX
bus to a different real bus (`docs/src/branchmodel.md`, "Three-winding
transformers"). `Branch`'s own `comp` is built per-call from `(vn_kV, from,
to, id)` (`getBranchComp`, `src/branch.jl:191-194`), and `from`/`to` differ
per winding (AUX→HV, AUX→MV, AUX→LV), so the three branches get distinct
`cName`/`cID` values. `_find_trafo_branch` can therefore pick a specific one
of the three deterministically **by branch identity** (name, id, or branch
index string) — it does not collide across windings by name alone, contrary
to the concern in the background note. `Branch` construction itself is
side-generic: `src/branch.jl:196`
(`w = (side in [1,2,3]) ? (side==1 ? branch.side1 : ...) : error(...)`) is
exercised identically for `side=1,2,3`, so a 3WT's `side3` branch goes
through the exact same PI-model derivation as any 2WT branch — no separate
3WT code path exists here.

**Step 2 is the actual blocker.** `PowerTransformer.comp` (`src/transformer.jl:656`)
is a single `AbstractComponent` for the whole transformer object (one
`cFrom_bus`/`cTo_bus` pair), even when the object has three windings
(`side1/side2/side3`, `src/transformer.jl:664-666`). The `findfirst` at
`src/tap_control.jl:308` can only ever succeed for whichever one of the three
legs happens to match that single recorded bus pair; for the other two legs,
`trafo_obj` comes back `nothing` and `addPowerTransformerControl!` throws
(`src/tap_control.jl:348-350`, `"no transformer object found for branch ..."`).
So even a `Branch` correctly resolved by name in step 1 cannot reliably be
matched back to its parent `PowerTransformer` in step 2 unless it happens to
be the one leg whose bus pair matches `PowerTransformer.comp`.

**And even where step 2 does resolve, step 3 is wrong for windings 2/3.**
`push!(net.trafos[trafo_obj].side1.controls, ctrl)` (`src/tap_control.jl:351`)
is unconditional — it always attaches the new controller to `side1.controls`,
never to `side2.controls`/`side3.controls`, regardless of which winding/branch
was actually resolved. `_tap_controllers` (`src/tap_control.jl:27-42`) and
`clearTapControllers!` (`src/tap_control.jl:49-58`) do iterate all three
sides when *reading* controllers, so a controller attached to `side2`/`side3`
directly (e.g. constructed by hand and appended to `winding.controls`) would
be picked up correctly by the outer loop — the gap is specifically in
`addPowerTransformerControl!`'s attachment step, not in the outer-loop
collection/execution path.

**Conclusion:** today, `addPowerTransformerControl!` cannot correctly target
"winding 2" or "winding 3" of a true 3WT — not because branch-name lookup is
ambiguous, but because (a) the parent-`PowerTransformer` lookup assumes a
1:1 branch↔transformer relationship that doesn't hold for a 3-winding object,
and (b) the attachment step is hard-coded to `side1` regardless. Note also
that this three-winding, single-`PowerTransformer` model (`create3WTWindings!`)
currently has **no call site wiring it into `net.trafos`/`net.branchVec`** at
all (see §3) — `add3WTPiModelTrafo!` (`src/network.jl:973-1024`), the only
3WT path actually wired into a `Net`, builds three fully independent 2WT
`PowerTransformer` objects instead (via three `add2WTPIModelTrafo!` calls,
each `side3 = nothing`), so this addressing gap is currently latent rather
than user-visible.

## 3. Where do 3WT branch tap ranges come from today?

`Branch`'s `PowerTransformer` constructor branch (`src/branch.jl:183-218`)
is the single derivation path for **any** side (1, 2, or 3):

- Ratio-tap range: defaults `tap_min = 0.9; tap_max = 1.1; tap_step = 0.00625`
  (`src/branch.jl:211-213`), overridden by `calcRatioTapRange(w.taps)`
  (`src/branch.jl:214-216`, delegating to `src/equicircuit.jl:255-259`) **iff**
  `w.taps !== nothing`.
- Phase-tap range: `phase_min_deg = -30.0, phase_max_deg = 30.0,
  phase_step_deg = 1.25` are **hard-coded literals** passed straight into
  `new(...)` at `src/branch.jl:218` — there is no equivalent
  `calcPhaseTapRange`-style branch that inspects `w.phase_taps`. A `grep` for
  `phase_taps` in `branch.jl` returns nothing.

So for a 3WT winding branch (`side=3`) exactly as for any 2WT branch: ratio
range is populated from `w.taps` when present, but phase range is always the
same three hard-coded constants — `w.phase_taps` (whether from a 2WT or a 3WT
winding) has **zero effect** on the branch's `phase_min_deg/phase_max_deg/
phase_step_deg`, or on anything else read by the solver. This confirms the
Part-1 background: attaching a `PhaseTapChangerModel` to a winding via
`create3WTWindings!`/`phase_tap_side` makes the model available on the data
model, but resolving it into effective branch parameters is a distinct,
unimplemented step — out of scope for this task by design.

## Summary — what's missing for outer-loop phase control of one 3WT winding

1. `create3WTWindings!`'s true 3-winding `PowerTransformer` (`side1/side2/side3`
   sharing one object) is not wired into any `Net` builder; `add3WTPiModelTrafo!`
   uses three independent 2WT objects instead (`src/network.jl:973-1024`).
2. `addPowerTransformerControl!`'s second lookup step
   (`src/tap_control.jl:308`) assumes one `Branch` maps to one
   `PowerTransformer.comp` bus pair; a genuine 3-winding `PowerTransformer`
   has only one `comp`, so at most one of its three winding-branches can ever
   resolve.
3. `addPowerTransformerControl!`'s attachment step (`src/tap_control.jl:351`)
   is hard-coded to `side1.controls`; there is no parameter to target
   `side2`/`side3`.
4. `Branch`'s phase-tap range fields (`phase_min_deg/phase_max_deg/
   phase_step_deg`, `src/branch.jl:218`) are hard-coded constants independent
   of `w.phase_taps` for every branch (2WT and 3WT alike) — a
   `PhaseTapChangerModel` attached via `phase_tap_side` has no numerical
   effect until this derivation is implemented.

None of the above was changed by this task; items 1-4 are prerequisites for
any future "outer-loop phase control of one 3WT winding" work.

## Appendix: pre-existing `tap_side` convention/implementation mismatch

While implementing `phase_tap_side` (kept implemented against its own,
directly-1-based `side == phase_tap_side` comparison, deliberately not
reusing the pattern below), the following pre-existing issue in
`create3WTWindings!`'s ratio-tap side selection was found and left
unmodified, per this task's explicit instruction not to silently fix
ratio-tap behaviour:

```julia
for side = 1:3
  tap = (side - 1 == tap_side) ? tap : nothing
  ...
end
```

- The docstring documents `tap_side` as 1-based (`"[1,2,3]. It is 0 if there
  is no tap"`), but the comparison `side - 1 == tap_side` is written as if
  `tap_side` were 0-based.
- Independently of that offset, the loop reassigns the **same** `tap` local
  binding (the function's own keyword argument) on every iteration. On any
  iteration where the condition is false, `tap` is overwritten with
  `nothing` — irreversibly, since it is the same binding read on the next
  iteration. Consequently, only `tap_side == 0` ever preserves the original
  `tap` object (because side `1` is evaluated first, before any overwrite
  can occur), and it attaches to **side 1** — not "no tap" as the docstring
  states for `tap_side == 0`. Every `tap_side` value that documentation
  suggests should attach a tap to a *specific* winding (`1`, `2`, or `3`)
  in fact discards the `tap` object entirely, leaving all three windings
  with `taps === nothing`.
- This was confirmed empirically (see `test/test_3wt_phase_taps.jl`, testset
  `"create3WTWindings! existing docstring example: no phase taps"`, which
  snapshots the current docstring example's `tap_side = 1` call producing
  `taps === nothing` on all three windings).

This is orthogonal to the phase-tap-changer work in this task and was left
untouched; fixing it is a separate, follow-up bugfix.
