# Rectangular refactor audit

## Rectangular entry points and call layers

### 1) Public user-facing entry points

- `runpf!(net::Net; ...)` in `src/jacobian_complex.jl`.
  - **Role:** main public power-flow entry that dispatches by configured method and forwards rectangular-related options through config/runtime arguments.
  - **Called by:** direct user calls and higher-level orchestration (`run_acpflow`).
  - **Placement recommendation:** remain in `src/jacobian_complex.jl` for now as public API glue; later consider a dedicated public solver entry module only after rectangular/other-method dispatch boundaries are stabilized.
  - **Visibility:** public.

- `runpf_rectangular!(...)` in `src/jacobian_complex.jl`.
  - **Role:** rectangular convenience/public-facing wrapper that prepares network-level arguments and calls the network-integrated rectangular solver.
  - **Called by:** `runpf!` rectangular branch and potentially direct user calls.
  - **Placement recommendation:** remain for now; later candidate to move to a dedicated `rectangular_network_entry.jl` (or similar) once the net-loop internals are extracted.
  - **Visibility:** public/semi-public (exported/externally callable compatibility-oriented entry).

### 2) Network-integrated solver

- `runpf_rectangular!(...)` in `src/jacobian_complex.jl`.
  - **Role:** network-integrated rectangular solve loop that still mixes orchestration, active-set/Q-limit handling, status table/workspace writes, result updates, and final status reporting.
  - **Called by:** `runpf_rectangular!`.
  - **Placement recommendation:** should eventually be reduced to orchestration-only in `src/jacobian_complex.jl`, with internal responsibilities extracted into focused helpers under `src/powerflow_rectangular/`.
  - **Visibility:** semi-public/internal interface (named and discoverable, but mainly internal pipeline entry).

### 3) Standalone array-level solver

- `run_complex_nr_rectangular(...)` in `src/powerflow_rectangular/rectangular_standalone_solver.jl`.
  - **Role:** array-level Newton driver on `(Ybus, V, S, bus_types, Vset, ...)` independent of net-level active-set orchestration.
  - **Called by:** `runpf_rectangular!`.
  - **Placement recommendation:** correct current location; keep in split helper module.
  - **Visibility:** internal/semi-public utility for rectangular stack.

### 4) Newton-step helper

- `complex_newton_step_rectangular(...)` in `src/powerflow_rectangular/rectangular_newton_step.jl`.
  - **Role:** computes one rectangular Newton update and applies fixed damping or autodamped backtracking selection.
  - **Called by:** `run_complex_nr_rectangular`.
  - **Placement recommendation:** keep in helper module; no move needed now.
  - **Visibility:** internal.

### 5) Interface/wrapper outside the split

- `run_acpflow(...)` and `_run_acpflow_net!(...)` in `src/run_acpflow.jl`.
  - **Role:** high-level orchestration/config ingestion; forwards method/autodamp/autodamp_min and related settings via `PowerFlowConfig` into `runpf!`.
  - **Called by:** user-facing ACP-flow workflows.
  - **Placement recommendation:** keep here; this is intentionally above solver internals.
  - **Visibility:** public (`run_acpflow`) + internal (`_run_acpflow_net!`).

- `runpf_external!(...)` in `src/solver_interface.jl`.
  - **Role:** external-solver interface path; included as a neighboring wrapper/interface layer but not part of rectangular core call chain.
  - **Called by:** solver interface consumers.
  - **Placement recommendation:** unchanged.
  - **Visibility:** public/semi-public interface.

## Remaining responsibilities in `src/jacobian_complex.jl`

1. **Network-integrated rectangular orchestration**
   - **Functions (approx.):** `runpf_rectangular!`, `runpf_rectangular!`.
   - **Current responsibility:** assembles network arrays, controls iteration policy, calls standalone solver/helper stack, orchestrates post-processing.
   - **Safe to extract now:** yes, incrementally.
   - **Next extraction target:** move non-entry helper blocks used only by `runpf_rectangular!` to dedicated helper file(s) (while keeping function signature stable).
   - **Risk:** medium (many keyword/threaded dependencies).

2. **Active Q-limit switching loop and decision plumbing**
   - **Functions/blocks (approx.):** Q-limit loop logic inside `runpf_rectangular!`, plus related adjustment hooks and counters.
   - **Current responsibility:** dynamic PV/PQ switching decisions, iteration-phase guard handling, adjust-vset interactions, trace/logging invocation.
   - **Safe to extract now:** conditionally yes, but requires careful keyword-forwarding and state object boundaries.
   - **Next extraction target:** dedicated internal helper that receives explicit workspace/state bundle.
   - **Risk:** high (behavior-sensitive convergence path).

3. **Final status construction/reporting and workspace table handling**
   - **Functions/blocks (approx.):** status writes/summaries and calls into `_RectangularPFStatusTable` / workspace helpers.
   - **Current responsibility:** result status normalization and summary output coordination.
   - **Safe to extract now:** yes.
   - **Next extraction target:** isolate status finalization and reporting glue into one helper file.
   - **Risk:** low-to-medium (mostly bookkeeping, but user-visible status strings/tags).

### Finalization extraction progress (current)

- **Completed in this step:** extracted low-risk post-iteration finalization helpers to `rectangular_finalization.jl` for:
  - bus-type sync to network node types,
  - solved-voltage write-back,
  - final injection vector computation,
  - bus P/Q result write-back,
  - total bus power reduction/write-back.
- **Still in `runpf_rectangular!`:**
  - active Q-limit switching workflow,
  - final Q-limit summary acceptance/rejection logic,
  - wrong-branch final diagnostics and status integration,
  - overall rectangular orchestration and status construction.

4. **Public wrapper and config-adapter surface**
   - **Functions (approx.):** `_runpf_with_config!`, `runpf!(...)` overloads, rectangular dispatch branches.
   - **Current responsibility:** compatibility and public API dispatch, runtime-key gating, config-based forwarding.
   - **Safe to extract now:** partially; wrappers should stay stable while internals move.
   - **Next extraction target:** thin helper for rectangular-only option marshalling (without changing API).
   - **Risk:** medium (public-call stability constraints).

5. **Legacy header/comments mismatch with split reality**
   - **Blocks:** top-of-file comments/doc text still describing larger monolithic ownership.
   - **Safe to extract/change now:** yes (doc-only).
   - **Next extraction target:** keep comments synchronized to current split and explicitly mark remaining non-extracted responsibilities.
   - **Risk:** low.

## Damping and autodamping policy check

1. **Keyword-definition locations**
   - `run_complex_nr_rectangular` defines `damp`, `autodamp`, `autodamp_min` in `src/powerflow_rectangular/rectangular_standalone_solver.jl`.
   - `runpf_rectangular!` defines the same trio in `src/jacobian_complex.jl`.
   - `runpf_rectangular!` defines and forwards the same trio in `src/jacobian_complex.jl`.
   - Newton-step/autodamp logic uses these in `complex_newton_step_rectangular` / `choose_rectangular_autodamp` in `src/powerflow_rectangular/rectangular_newton_step.jl`.

2. **Where config/YAML enters rectangular path**
   - YAML/config defaults feed `PowerFlowConfig` in `src/configuration.jl` (`autodamp`, `autodamp_min`).
   - `run_acpflow` and `_run_acpflow_net!` load/forward these into solver config in `src/run_acpflow.jl`.
   - `_runpf_with_config!` in `src/jacobian_complex.jl` forwards config values into `runpf!`/`runpf_rectangular!` path.
   - `damp` is runtime keyword-enabled in `runpf!(net::Net; ...)` (from raw runtime config map via `get(raw, "damp", 1.0)`), not a `PowerFlowConfig` field in current code.

3. **Operational semantics of `damp`**
   - With `autodamp=false`, `complex_newton_step_rectangular` applies fixed `damp` as direct step multiplier.
   - With `autodamp=true`, autodamp search starts from `alpha = damp` and backtracks toward `autodamp_min`; therefore `damp` acts as initial upper trial factor.
   - This behavior is consistent with helper doc comments in `rectangular_newton_step.jl`.

4. **Validation of `autodamp_min <= damp`**
   - `_validate_rectangular_damping` enforces `0 < autodamp_min <= damp <= 1` and is called in both autodamp-selection path and fixed-damp path (with clamped secondary argument in fixed mode).
   - Validation appears consistent in active solver paths.

5. **Documentation clarity**
   - Core helper-level semantics are documented.
   - User-facing docs are still somewhat fragmented across comments/signatures/config; a single concise policy note in user-facing solver docs would improve clarity.

6. **Default consistency check (resolved)**
   - Previously, direct rectangular solver signatures used `autodamp_min=1e-3` while `PowerFlowConfig`/YAML used `0.05`.
   - Canonical user-visible default is now aligned to `autodamp_min=0.05` across config and direct rectangular entry signatures.
   - The operational policy remains unchanged: with `autodamp=false`, `damp` is a fixed step factor; with `autodamp=true`, autodamping backtracks from `damp` down to `autodamp_min`.
   - Validation rule remains `0 < autodamp_min <= damp <= 1`.

7. **Naming/documentation improvement recommendation (no behavior change in this task)**
   - Keep parameter name `damp` (no API rename).
   - Document explicitly in one place: “`damp` is fixed factor when `autodamp=false`; with `autodamp=true`, it is the initial trial factor/upper bound for backtracking.”
   - Add a follow-up documentation note/warning about `autodamp=true` with `damp < 1` (smaller initial search ceiling).
   - Follow-up task should decide whether to align default `autodamp_min` values across config and direct-function signatures.

## Recommended next small tasks

1. **Goal:** reduce `runpf_rectangular!` to orchestration-only by extracting status/finalization glue into one helper.
   - **Files touched:** `src/jacobian_complex.jl`, `src/powerflow_rectangular/rectangular_status_workspace.jl`, `src/powerflow_rectangular/README.md`.
   - **Risk:** medium.
   - **Validation command:** `julia --project=. test/test_solver_interface.jl`.
   - **Behavior changes allowed:** no.

2. **Goal:** extract active Q-limit switching sub-loop into a focused internal helper API (no algorithm change).
   - **Files touched:** `src/jacobian_complex.jl`, new helper under `src/powerflow_rectangular/`.
   - **Risk:** high.
   - **Validation command:** `julia --project=. test/runtests.jl`.
   - **Behavior changes allowed:** no.

3. **Goal:** unify/document damping policy in public-facing docs/comments and clarify `damp` dual role.
   - **Files touched:** `src/powerflow_rectangular/README.md`, relevant doc comments in `src/jacobian_complex.jl` and/or `src/run_acpflow.jl` (doc-only).
   - **Risk:** low.
   - **Validation command:** `julia --project=. -e "using Sparlectra"`.
   - **Behavior changes allowed:** no.

4. **Goal:** keep rectangular damping-policy docs/tests synchronized with the now-aligned `autodamp_min=0.05` default.
   - **Files touched:** `src/configuration.jl`, `src/configuration.yaml.example`, rectangular entry signatures/docs, tests.
   - **Risk:** medium.
   - **Validation command:** `julia --project=. test/runtests.jl`.
   - **Behavior changes allowed:** yes (explicitly controlled follow-up decision).
