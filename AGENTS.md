# AGENTS.md

This repository is a Julia project.

## Rules
- Use Julia for code changes, scripts, diagnostics, and tests.
- Always run Julia commands with `julia --project=.`
- Do not use Python to inspect or modify Julia source unless explicitly requested.
- Prefer minimal, targeted changes.
- Keep changes local to the affected implementation path whenever possible.
- Do not silently change public APIs or example-script behavior in bugfix or minor-improvement tasks.
- For explicit refactoring or breaking-change tasks, prefer the clean target design over backward-compatibility layers.
- Do not add compatibility shims, deprecated wrappers, legacy keyword paths, or transitional duplicate APIs unless explicitly requested.
- When an API, option, file path, or example behavior is intentionally changed, update all in-repository callers, examples, documentation, and tests instead of preserving the old path.

## Agent execution discipline

These rules are mandatory for automated coding agents.

- Treat every task list as a checklist. If a prompt contains multiple tasks, complete all of them unless the user explicitly asks for only one.
- Do not choose an easier subset of the task. If a requested item is unsafe, impossible, ambiguous, or outside the repository state, state that explicitly and leave it unclaimed.
- Do not report "done", "ready", "complete", or "merge-ready" after partial work.
- If work is partial, say `PARTIAL` clearly in the summary and list every unfinished item.
- Do not commit partial cleanup as if it were successful completion. A commit may be made only when the requested acceptance criteria are satisfied, or when the user explicitly asked for a partial checkpoint.
- Before committing, re-read the task acceptance criteria and verify each item.
- Do not use broad refactoring to hide an incomplete targeted fix.
- If a task asks for proof, provide the exact command output or a precise line-by-line classification of remaining matches.
- If a command fails, do not continue as if the task passed. Fix the smallest cause, rerun the failed command, and report both the failure and the rerun result.
- Never claim a test or verification command was run unless it was actually run in the current working tree.


## Repository freshness and current-HEAD gate

These rules are mandatory before any automated coding agent edits files.

- Always work from the current repository state, not from a stale patch, stale mental model, previous Codex run, or copied block from an old commit.
- Before editing, run and report:
  ```bash
  git status --short
  git rev-parse --short HEAD
  git log --oneline -n 5
  ```
- If a remote/upstream branch exists, also run and report:
  ```bash
  git fetch --all --prune
  git status -sb
  git rev-list --left-right --count HEAD...@{u}
  ```
- If the working tree is dirty, stop immediately. Do not edit files. Report the dirty files.
- If the local branch is behind its upstream, stop immediately. Do not auto-merge, auto-rebase, or continue editing unless the user explicitly requested that action.
- If the local branch is ahead of upstream, report that fact before editing.
- If no upstream is configured, report that fact explicitly and continue only from the current local HEAD.
- For bugfix tasks, reproduce the current failure from the current HEAD before editing. Do not rely only on pasted logs from an earlier commit.
- Never claim that a task was performed on the latest commit unless the current HEAD and repository status were checked in the same run.

## Stale-patch and merge-conflict prevention

- Do not apply a patch generated against an older commit without first verifying that it applies cleanly to the current HEAD.
- Do not paste or reconstruct large code blocks from an old commit unless the task explicitly asks for a restore and the complete dependency block has been identified.
- Do not restore symbols one by one from history when they form a call chain. Restore or repair the complete coherent block.
- Before changing a function body, inspect the current function signature, all callers, and all forwarded keywords in the current HEAD.
- If a function body from history references a keyword or helper that is not present in the current signature or current file, do not paste it as-is.
- If a task requires Git conflict resolution, stop and report the conflict. Do not guess the intended resolution.
- If a merge conflict marker is present (`<<<<<<<`, `=======`, `>>>>>>>`), stop and report it unless the task explicitly asks to resolve merge conflicts.
- Do not commit code that only compiles by accident while leaving broken call chains for the next test step.

## Protected core solver files

The following files are protected because they contain central solver, configuration, or execution-path logic:

```text
src/jacobian_complex.jl
src/run_acpflow.jl
src/control_framework.jl
src/configuration.jl
src/Sparlectra.jl
```

For these files:

- Do not perform broad autonomous refactoring unless the user explicitly requested that exact refactoring.
- Do not mix refactoring with feature work or bugfix work.
- Do not make large formatting-only changes.
- Do not restore partial snippets from old commits.
- Do not add compatibility shims, wrapper layers, or duplicate APIs unless explicitly requested.
- Do not change public entry points such as `runpf!`, `runpf_rectangular!`, or `run_acpflow` without explicitly checking all in-repository callers.
- Do not change keyword signatures without checking the complete forwarding chain.
- Do not change solver status metadata without checking all readers, tests, and reporting code.
- Do not edit tests to hide a missing public or internal API binding.

When touching these files, first run:

```bash
julia --project=. -e "using Sparlectra"
```

Then check any relevant binding chain explicitly. For the rectangular power-flow path, use:

```bash
julia --project=. -e "using Sparlectra; for s in Symbol.(String[\"runpf!\", \"runpf_rectangular!\", \"run_complex_nr_rectangular\"]); @assert isdefined(Sparlectra, s) s; @assert !isempty(methods(getfield(Sparlectra, s))) s; println(s, \" => \", length(methods(getfield(Sparlectra, s))), \" method(s)\"); end; @assert isdefined(Main, Symbol(\"runpf!\"))"
```

If any required binding is missing, stop and report the missing binding before editing further.

## Core solver recovery policy

If a protected solver file is syntactically or semantically broken:

- First restore package loadability.
- Do not start feature work while `using Sparlectra` fails.
- Do not split, refactor, or redesign the solver while the baseline is broken.
- If a previous partial restore caused missing symbol chains, identify the complete call chain before editing.
- For the rectangular power-flow path, keep this chain coherent:
  ```text
  runpf! -> runpf_rectangular!
  ```
- If a restored wrapper body references keywords such as `validate_limits_after_pf`, `q_limit_violation_headroom`, `pv_table_rows`, Q-limit options, or wrong-branch options, verify that the active method signature binds them and that callers forward them consistently.
- Do not add global variables or dummy functions for missing keyword names. Missing keywords must be fixed in the relevant method signature and forwarding chain.
- If the correct repair would require restoring a large historical block, stop and ask for explicit permission unless the user already requested a restore.
- Prefer a complete file restore from the last known passing commit over incremental patching when the current file is corrupted.

## Forbidden actions unless explicitly requested

Automated coding agents must not perform the following actions unless the user explicitly requests them:

- active wrong-branch rescue/retry loop,
- Monte-Carlo or multi-start solver scans,
- continuation or homotopy start logic,
- new DC-start or QV-start fallback modes,
- broad solver restructuring mixed with a bugfix,
- large-scale rewrite of `src/jacobian_complex.jl`,
- changing legacy public APIs to hide test failures,
- converting many tests from unqualified public calls to qualified internal calls,
- adding global variables or dummy functions to satisfy missing local keyword errors,
- creating new changelog files,
- creating broad compatibility layers after an intentional API cleanup.

## Required stop conditions

Stop immediately and report the first blocker if any of the following occurs:

- working tree is dirty before editing,
- local branch is behind upstream,
- `using Sparlectra` fails before a task that requires code changes,
- the first failing test cannot be reproduced from the current HEAD,
- a requested source file contains merge conflict markers,
- an intended edit would require guessing between conflicting implementations,
- a protected solver file has missing call-chain dependencies,
- validation reveals a new missing-symbol chain caused by the current edit.

Do not continue to the next task item after a stop condition.

## Multi-task and checklist handling

For prompts containing several tasks, issues, checklist items, or numbered sections:

1. Convert the request into an internal checklist before editing.
2. Work through every checklist item.
3. After editing, revisit the checklist and mark each item as:
   - done,
   - not applicable, with reason,
   - blocked, with reason,
   - not done.
4. The final response must include unresolved or blocked items. Do not omit them.
5. If only some tasks were completed, the final response must not recommend merging or release.

## Acceptance criteria and evidence

Acceptance criteria are binding. Search/evidence commands in a task are also binding.

- If the task provides `rg`, `find`, `julia`, or documentation-build commands, run them exactly unless the user explicitly says not to.
- If a search command still returns matches, classify every remaining match.
- Use this classification vocabulary when applicable:
  - `expected obsolete-key rejection`
  - `legacy/non-production path`
  - `solver-interface compatibility boundary`
  - `historical changelog`
  - `test for removed behavior`
  - `still to fix`
- No unclassified leftover search match is acceptable.
- A branch is not merge-ready if any required command failed or any remaining match is classified as `still to fix`.
- For API removals, update implementation, tests, examples, and documentation in the same change. Do not leave stale callers to fail later.
- For option removals, search the full repository for the option name and update or classify every match.

## Commit discipline

- Make a commit only after the requested verification passes, unless the user explicitly requests an intermediate commit.
- Do not create PR metadata that claims completion when the work is partial.
- Commit messages must describe the actual completed change, not the intended full task if only part was finished.
- If verification is not run, say so in the final response and do not mark the task as verified.
- If verification fails after changes, either fix the failure or leave the branch uncommitted unless explicitly asked otherwise.

## Standard commands
- Verify toolchain:
  `which julia && julia --version`
- Instantiate environment:
  `julia --project=. -e 'using Pkg; Pkg.instantiate()'`
- Run tests:
  `julia --project=. test/runtests.jl`
- Alternative full test:
  `julia --project=. -e 'using Pkg; Pkg.test()'`
- Run with deprecation warnings as errors when investigating Julia 1.12 world-age or Revise issues:
  `julia --project=. --depwarn=error path/to/script.jl`


## Verification policy

Use the smallest meaningful verification first, then broader verification when the task affects shared behavior.

Default verification levels:

- Syntax/package load check:
  `julia --project=. -e 'using Sparlectra'`
- Fast test profile:
  `julia --project=. test/runtests.jl`
- Extended test profile:
  `SPARLECTRA_TEST_PROFILE=extended julia --project=. test/runtests.jl`
- Documentation build:
  `julia --project=docs docs/make.jl`

Required rules:

- If a task changes package loading, module exports, function signatures, or central config, run the package load check.
- If a task changes solver behavior, configuration, tests, examples, MATPOWER runner behavior, or documentation, run the relevant profile or explicitly state why it was not run.
- If a task is intended to make a branch merge-ready, run package load, fast profile, extended profile, and docs build unless the user explicitly excludes one.
- If any verification command fails, the final response must say `NOT READY`.
- Do not hide warnings that indicate future errors, especially Julia world-age warnings under Julia 1.12. Fix them if they are in touched code or test runner code.

PowerShell extended-profile pattern:

```powershell
$env:SPARLECTRA_TEST_PROFILE="extended"
julia --project=. test/runtests.jl
Remove-Item Env:\SPARLECTRA_TEST_PROFILE
```

## Repo conventions
- Keep code comments and docstrings in English.
- Prefer existing Sparlectra naming and style.
- Keep `src/Sparlectra.jl` readable as the public module entry point.
- Important public exports in `src/Sparlectra.jl` must have short inline comments or concise group comments that explain their role at the module-interface level.
- Do not comment every trivial helper in `src/Sparlectra.jl`; focus on important user-facing functions, solver entry points, configuration APIs, import/export APIs, diagnostics, and external interfaces.
- When adding, removing, or renaming important exports in `src/Sparlectra.jl`, update the surrounding comments in the same change.
- Before editing, inspect current tests that cover the affected code.
- After code changes, run the smallest relevant test first, then broader tests if needed.
- Avoid broad refactoring while fixing a targeted bug.
- During planned refactoring, remove obsolete code paths instead of keeping them as compatibility ballast.

## Refactoring and breaking changes

For tasks that explicitly request refactoring, cleanup, simplification, or breaking changes, backward compatibility is **not** required by default.

Use these rules:
- Prefer one clean implementation path over parallel old/new code paths.
- Remove obsolete wrappers, deprecated keyword forwarding, legacy YAML files, compatibility aliases, and unused fallback branches.
- Update internal callers immediately instead of preserving old signatures.
- Update examples and documentation to the new behavior.
- Remove tests that only protect removed compatibility behavior.
- Do not keep deprecation warnings or migration layers unless the task explicitly asks for a transition period.
- If public behavior changes intentionally, state this clearly in the changelog or task summary.
- If a breaking change affects documented usage, update the documentation in the same change.


## Test documentation maintenance

When adding, removing, renaming, moving, or materially changing tests, update the test-suite documentation in `docs/src/tests.md` in the same change.

The documentation must stay aligned with:

- available test profiles,
- test runner commands,
- files loaded by each profile,
- logical test groups,
- what each group verifies,
- offline/extended/all expectations.

Do not add new test files or new profile groups without documenting them.

## Test hygiene / avoiding test explosion

Tests are mandatory for bugfixes, changed behavior, and new features, but the test suite must remain maintainable. Avoid test explosion.

Default policy:
- Extend existing tests first.
- Add a new test file only when no suitable existing test location exists.
- Remove tests for deleted behavior, obsolete compatibility paths, and temporary diagnostics.
- Consolidate overlapping tests into one coherent testset or a table-driven test.
- Prefer fewer, stronger tests over many narrowly duplicated tests.

Before adding a new test file or a new large test block:
- Inspect existing tests for the affected code path.
- Check whether an existing test can be extended instead of adding a new near-duplicate test.
- Prefer merging related regression checks into one coherent testset when they exercise the same function, option path, or failure mode.
- Avoid preserving temporary diagnostic tests that were only useful during bug investigation.
- Remove or consolidate obsolete tests when the underlying failure mode is already covered by a clearer regression test.
- Remove tests that only verify backward compatibility for APIs, files, options, or wrappers intentionally removed by a refactoring task.
- Keep test names explicit about the behavior being protected, not about the historical debugging session.
- Keep test fixtures small; do not add large MATPOWER cases or long-running examples unless they are explicitly required.
- Avoid broad end-to-end tests when a focused unit or integration test can protect the same behavior.
- When adding a test for a bugfix, document the protected failure mode in a short inline comment if it is not obvious from the test name.
- If multiple tests cover the same behavior with only small parameter differences, use a table-driven test or a loop inside one testset instead of duplicating test bodies.

During cleanup or refactoring tasks:
- Revisit recently added tests and decide whether they are still needed.
- Merge redundant tests where possible.
- Delete tests that only assert implementation details and no longer protect public or intended internal behavior.
- Keep at least one regression test for each fixed bug, but avoid multiple overlapping regressions for the same failure mode.
- For breaking changes, keep tests for the new intended behavior and delete tests that only protect the removed behavior.

When reporting work:
- State which smallest relevant test was run.
- If tests were added, state whether an existing test was extended or why a new test was necessary.
- If tests were consolidated or removed, state which behavior remains covered.

## Code hygiene / avoiding implementation bloat

Before completing a task, inspect the touched implementation path for unnecessary growth.

- Look for duplicated helper logic introduced during debugging.
- Remove dead code, stale debug branches, and obsolete temporary diagnostics unless they are intentionally kept behind a clear diagnostic option.
- Prefer one well-named helper over repeated local code blocks.
- Keep logging and diagnostics compact by default; verbose output must be opt-in.
- Avoid adding options that merely paper over one case unless the behavior is clearly documented and reusable.
- Do not retain experimental code paths without tests and documentation.
- If a cleanup would become too broad for the current bugfix, leave a short follow-up issue instead of mixing unrelated refactoring into the fix.

### Julia indexing style
To avoid linter warnings (`IndexFromLength`):

- Do **not** write loops like:
  `for i = 1:size(A,1)`

- Prefer array axes:
  `for i in axes(A,1)`

- Prefer `eachindex` when iterating over arrays:
  `for i in eachindex(A)`

This keeps the code compatible with non-standard array indices and avoids StaticLint warnings.

## For this repo
- Main language: Julia
- Project root must be used as Julia project
- When changing power flow or network logic, prefer adding or updating tests in `test/`
- Prefer extending or consolidating existing tests over creating new narrowly overlapping test files.
- Remove obsolete tests when refactoring deletes the behavior they were written for.

## VS Code / developer workflow notes
- Example programs in `examples/` are primarily written for VS Code developer usage.
- When practical, prefer calling script entry functions via `Base.invokelatest(...)` to reduce world-age issues in interactive development sessions.
- Do not rely only on `if abspath(PROGRAM_FILE) == @__FILE__` for developer workflows, since this may not match all VS Code execution modes.
- Recommend `Revise.jl` for contributors working interactively; keep it out of project dependencies so Sparlectra installation remains fast and lightweight.

## Julia 1.12 / Revise / world-age handling

Julia 1.12 reports stricter warnings for accessing newly defined global bindings from an older world age. These warnings can appear in VS Code or Revise-based workflows, especially in example scripts that define `main()` or benchmark entry functions and call them immediately.

Typical warning:

```text
WARNING: Detected access to binding `Main.main` in a world prior to its definition world.
Julia 1.12 has introduced more strict world age semantics for global bindings.
Hint: Add an appropriate `invokelatest` around the access to this binding.
```

Julia 1.12 world-age warnings are not harmless noise in this repository. Treat them as future errors.

- If a task touches a script entry point, test runner, or dynamically included test file, avoid `getfield(@__MODULE__, name)` followed by a stale call.
- Prefer:
  ```julia
  runner = Base.invokelatest(getfield, @__MODULE__, name)
  Base.invokelatest(runner)
  ```
- If a warning appears in verification output from touched code, fix it before claiming the task is complete.

### Required handling
- For script entry points, call the entry function through `Base.invokelatest`, for example:
  `Base.invokelatest(main)`
- If the entry function requires arguments, use:
  `Base.invokelatest(main, args...)`
- If a benchmark or example script dynamically defines and then calls another local entry function, call that function through `Base.invokelatest` at the call site.
- Do not wrap internal hot-loop functions with `invokelatest`; use it only at script/developer entry boundaries where world-age problems occur.
- In example scripts included into `Main` or ad-hoc modules, qualify Sparlectra enum constants and types with `Sparlectra.` (for example `Sparlectra.Slack`, `Sparlectra.PV`, `Sparlectra.PQ`) unless they are explicitly imported in that file.
- When calling a helper function defined in the same example file from a `main`, `show_once`, benchmark, or redirected-output closure, prefer `Base.invokelatest(getfield(@__MODULE__, :helper_name), args...; kwargs...)` to avoid Julia 1.12/Revise world-age binding warnings.
- When investigating this class of issue, reproduce with:
  `julia --project=. --depwarn=error path/to/script.jl`
  so that Julia produces a stack trace.

### Keyword-forwarding rule for example scripts

Many example scripts forward large keyword sets through wrapper functions. When adding a new keyword option, every wrapper in the call chain must be updated consistently.

If a call fails with an error such as:

```text
MethodError: no method matching bench_run_acpflow(...; autodamp=..., start_projection=...)
got unsupported keyword arguments "autodamp", "start_projection", ...
```

then the task is usually **not** only a world-age issue. It means that a caller forwards keyword arguments which the callee does not accept.

When fixing this:
- Find the complete call chain, for example:
  `main(...) -> run_case(...) -> bench_run_acpflow(...) -> run_net_acpflow(...)`
- Add the new keyword arguments to all relevant wrapper signatures if they are intended to be supported.
- Forward the keywords explicitly to the next layer.
- If a keyword is not meaningful for a specific wrapper, filter it intentionally and document why.
- Do not pass arbitrary keyword dictionaries blindly into functions unless the existing style already does this.
- Prefer explicit keyword signatures for public or developer-facing example entry points.
- Add or update a smoke test if the script is part of tested behavior.
- Before adding a new keyword-forwarding test, check whether an existing wrapper smoke test can be extended.

### Example pattern

Prefer:

```julia
function main(; start_projection::Bool=false, autodamp::Bool=false, kwargs...)
    return Base.invokelatest(run_example; start_projection=start_projection, autodamp=autodamp, kwargs...)
end

Base.invokelatest(main)
```

Only use `kwargs...` when this matches the local style and the downstream function is designed to accept them. Otherwise add explicit keywords.

## Required after code changes
- After modifying Julia code, run the smallest relevant test.
- If no narrower test is obvious, run:
  `julia --project=. test/runtests.jl`
- Report the exact command and the relevant failure output if tests fail.
- Before adding new tests, check whether existing tests can be extended or merged.
- Before finishing, check whether new tests created during debugging can be consolidated.
- Remove tests that only cover removed compatibility paths or deleted legacy behavior.
- Avoid leaving behind redundant tests that only differ in incidental parameters.


## Final response contract for coding tasks

Every coding-agent final response must use this structure when code was changed:

```text
Summary
- ...

Verification
- command:
  result: PASS/FAIL/NOT RUN
  notes: ...

Search evidence
- command:
  result/classification: ...

Unfinished items
- none
```

Rules:

- `Unfinished items` must not be omitted.
- If nothing is unfinished, write `none`.
- If tests were not run, write `NOT RUN`; do not imply they passed.
- If a requested verification command failed, the merge recommendation must be `NOT READY`.
- If only a search/evidence pass was requested, include the exact search command and summarize all remaining matches.

## Power-flow execution preference
- Prefer `run_net_acpflow(...)` over calling `runpf!(...)` directly when practical.
- Reason: `run_net_acpflow(...)` also executes post-processing (e.g. `calcNetLosses!`,
  branch-flow calculation, and link-flow calculation) so reported network losses
  and branch/link flows are populated consistently.

## Definition: Feature vs Improvement vs Bugfix

To ensure consistent review standards, changes are classified as follows:

### Feature
A change is considered a **feature** if it introduces new user-visible capabilities or extends the modeling or solver functionality.

Typical indicators:
- New public API functions or significant extensions of existing APIs
- New modeling elements (e.g. controllers, network components, constraints)
- New solver capabilities or modes (e.g. new algorithms, control strategies)
- New data processing pipelines (e.g. measurement models, converters)
- New outputs or result types not previously available

Requirements for features:
- Inline docstrings (mandatory)
- Tests (mandatory)
- Example in `examples/exp_<feature>.jl` (mandatory)
- Documentation update (mandatory)
- Feature matrix review/update (mandatory)
- Test-hygiene review (mandatory): prefer extending existing tests and avoid duplicate test coverage.

---

### Improvement
A change is an **improvement** if it enhances existing functionality without changing its fundamental behavior or adding new capabilities.

Typical indicators:
- Performance optimizations
- Numerical stability improvements
- Refactoring without API change
- Improved error handling or logging
- Internal restructuring

Requirements for improvements:
- Tests updated if behavior is affected
- Docstrings updated if behavior or usage clarity changes
- No example required (unless it improves clarity significantly)
- No mandatory feature matrix update
- Review whether existing tests can be consolidated if the improvement replaces older behavior.
- If the improvement is an explicit refactoring or cleanup, remove obsolete compatibility code and its tests instead of preserving both old and new paths.

---

### Bugfix
A change is a **bugfix** if it corrects incorrect behavior relative to documented or expected functionality.

Typical indicators:
- Fixing wrong numerical results
- Correcting sign errors, indexing errors, or logic flaws
- Resolving crashes or invalid states

Requirements for bugfixes:
- Add or extend a test that reproduces the bug
- Prefer extending an existing regression test over adding a new test file
- Keep changes minimal and targeted
- No example required
- Documentation update only if the previous behavior was documented incorrectly
- Remove temporary diagnostic tests once a clear regression test exists.

---

### Edge cases
- If a change both fixes a bug and introduces new capabilities → classify as **feature**
- If unsure → default to **feature** to enforce full documentation and examples

## Sparlectra-specific classification examples

### Feature examples
Classify as **feature** when adding or exposing new capabilities, for example:

- New network elements:
  - FACTS devices
  - tap-changing transformers
  - phase-shifting transformers
  - new controller models
  - new load or generator models

- New power-flow behavior:
  - new Q-limit handling modes
  - voltage-setpoint adjustment strategies
  - remote voltage control
  - branch-flow target control
  - new slack or distributed-slack strategies

- New solver functionality:
  - new NR formulation
  - new linear solver option
  - APSLF / HELM-related solver integration
  - NR-polish workflow
  - new convergence or damping strategy

- New state-estimation functionality:
  - new measurement type
  - measurement generator
  - bad-data detection
  - observability analysis extensions
  - deactivate-and-rerun diagnostics

- New import/export capability:
  - MATPOWER conversion extensions
  - CGMES mapping improvements that expose new model semantics
  - new result export format

Feature requirements:
- Inline docstrings
- Targeted tests
- Example program in `examples/exp_<feature>.jl`
- Documentation update
- Feature matrix review/update
- Test-hygiene review: avoid duplicating existing coverage and consolidate related test cases where possible.

---

### Improvement examples
Classify as **improvement** when existing behavior is made better without adding a new capability, for example:

- Faster Y-bus assembly with unchanged results
- Cleaner internal Jacobian construction
- Better logging for existing solver steps
- More robust convergence checks using existing criteria
- Refactoring branch-flow calculation without API change
- Improving existing examples without adding new functionality
- Making existing error messages clearer

Improvement requirements:
- Update tests if behavior or output changes
- Update docstrings if user-facing behavior becomes clearer or changes
- Documentation update only if existing documentation becomes incomplete or misleading
- Consolidate or remove obsolete tests if the improvement replaces previous diagnostic-only behavior.

---

### Bugfix examples
Classify as **bugfix** when correcting wrong behavior, for example:

- Wrong transformer tap sign or phase-shift convention
- Incorrect branch-flow direction
- Incorrect link-flow post-processing
- Incorrect PV/PQ switching behavior
- Indexing errors in bus or branch mapping
- Incorrect treatment of shunt admittances
- Incorrect result tables or loss calculations
- Crashes for valid networks or valid measurements
- Example-script keyword mismatch after adding solver options
- Julia 1.12 world-age warnings/errors in VS Code or Revise workflows

Bugfix requirements:
- Add or extend a regression test
- Prefer an existing testset if it already covers the affected code path
- Keep the change minimal and targeted
- Update documentation only if the previous documented behavior was wrong
- No new example required unless it helps demonstrate the corrected behavior
- Do not keep multiple overlapping regression tests for the same bug.

---

### Sparlectra rule of thumb
If the change creates a new way for users to model, solve, estimate, diagnose, import, export, or interpret a network, classify it as a **feature**.

If it only makes an existing way faster, cleaner, safer, or more robust, classify it as an **improvement**.

If it corrects behavior that was wrong, classify it as a **bugfix**.

## Documentation and examples for new features
- New functions must include inline docstrings explaining purpose, inputs, outputs, and behavior.
- Docstrings must be written in English and follow existing Sparlectra style.

- For **new features** (i.e. functionality that adds a new capability, not bugfixes or minor improvements):
  - Add at least one example program in `examples/`
  - Naming convention: `exp_<feature>.jl`
  - The example should demonstrate typical usage and expected behavior

- For **new features**, also ensure:
  - The existing documentation is reviewed and updated where appropriate
  - The feature is added to a suitable documentation page (theory, usage, or API section)
  - The feature matrix is checked and updated if applicable
  - The tests are reviewed for overlap before adding new files or large new test blocks

- Tests alone are not sufficient for new features; documentation and examples are required.

## Changelog and release notes
- For each new feature, improvement, or bugfix, add a corresponding entry in `CHANGELOG.md` under the appropriate version heading.
- Follow the format:
  - **Feature**: Brief description of the new capability, its purpose, and any important usage notes.
  - **Improvement**: Brief description of the enhancement, what it improves, and any impact on existing behavior.
  - **Bugfix**: Brief description of the issue fixed, its impact, and any relevant details for users.

## Documentation (Documenter.jl / Markdown + Math)

To ensure stable and predictable rendering in Documenter.jl:

### Math usage
- Use `$...$` for **inline math** (variables inside text).
- Use ```math blocks only for **actual equations**, not for single symbols.
- Do **not** start list items with `$...$` → this may be rendered as display math (centered).
  - Instead write:
    - ✅ `- variable $y_{ik}$: ...`
    - ❌ `- $y_{ik}$: ...`

### Lists and definitions
- Avoid placing math expressions alone at the beginning of list items.
- Prefer descriptive text before math symbols for stable layout.
- Use `—` or text instead of relying on punctuation directly after math.

### Code vs. Math
- Use backticks `` `...` `` for:
  - variable names in code context (`y_ser`, `addACLine!`)
- Use math `$...$` only for mathematical notation (`y_{ik}`, `Y_{ii}`)

### LaTeX correctness
- Use valid LaTeX only (e.g. `$Y_{ii}$`, not `$\Y_{ii}$`)
- Units should be written as:
  - `$V_i = 1\,\\mathrm{pu}$`

### General rule
- Inline math = inline context
- Math blocks = equations only
- Never mix both unintentionally

### Goal
Keep documentation:
- left-aligned
- stable across Documenter versions
- consistent with KaTeX rendering

---


## Private project naming / confidentiality

- Do not use or reference the term `APSLF` in public documentation, examples, comments, changelogs, issue discussions, generated texts, or commit messages unless explicitly requested by the repository owner.
- `APSLF` is considered a private/internal project designation.
- Prefer neutral wording such as:
  - `analytical power-series solver`
  - `power-series load-flow method`
  - `analytical embedding solver`
  - `series-based load-flow approach`
- Do not introduce the term `APSLF` into newly generated public-facing repository artifacts.
- Existing historic references may remain unchanged unless the task explicitly includes renaming or cleanup work.
- If unsure whether a text is public-facing, treat it as public-facing by default.

## Do not add new changelog!
- Changelog remains here docs/src/changelog.md
- Add new entries to that file, do not create a new changelog or add entries to the README or other documentation files.
