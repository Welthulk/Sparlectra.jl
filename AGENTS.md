# AGENTS.md

This repository is a Julia project.

## Rules
- Use Julia for code changes, scripts, diagnostics, and tests.
- Always run Julia commands with `julia --project=.`
- Do not use Python to inspect or modify Julia source unless explicitly requested.
- Prefer minimal, targeted changes.
- Keep changes local to the affected implementation path whenever possible.
- Do not silently change public APIs or example-script behavior.

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

## Repo conventions
- Keep code comments and docstrings in English.
- Prefer existing Sparlectra naming and style.
- Before editing, inspect current tests that cover the affected code.
- After code changes, run the smallest relevant test first, then broader tests if needed.
- Avoid broad refactoring while fixing a targeted bug.

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

## VS Code / developer workflow notes
- Example programs in `src/examples/` are primarily written for VS Code developer usage.
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

### Required handling
- For script entry points, call the entry function through `Base.invokelatest`, for example:
  `Base.invokelatest(main)`
- If the entry function requires arguments, use:
  `Base.invokelatest(main, args...)`
- If a benchmark or example script dynamically defines and then calls another local entry function, call that function through `Base.invokelatest` at the call site.
- Do not wrap internal hot-loop functions with `invokelatest`; use it only at script/developer entry boundaries where world-age problems occur.
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
- Example in `src/examples/exp_<feature>.jl` (mandatory)
- Documentation update (mandatory)
- Feature matrix review/update (mandatory)

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

---

### Bugfix
A change is a **bugfix** if it corrects incorrect behavior relative to documented or expected functionality.

Typical indicators:
- Fixing wrong numerical results
- Correcting sign errors, indexing errors, or logic flaws
- Resolving crashes or invalid states

Requirements for bugfixes:
- Add or extend a test that reproduces the bug
- Keep changes minimal and targeted
- No example required
- Documentation update only if the previous behavior was documented incorrectly

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
- Example program in `src/examples/exp_<feature>.jl`
- Documentation update
- Feature matrix review/update

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
- Keep the change minimal and targeted
- Update documentation only if the previous documented behavior was wrong
- No new example required unless it helps demonstrate the corrected behavior

---

### Sparlectra rule of thumb
If the change creates a new way for users to model, solve, estimate, diagnose, import, export, or interpret a network, classify it as a **feature**.

If it only makes an existing way faster, cleaner, safer, or more robust, classify it as an **improvement**.

If it corrects behavior that was wrong, classify it as a **bugfix**.

## Documentation and examples for new features
- New functions must include inline docstrings explaining purpose, inputs, outputs, and behavior.
- Docstrings must be written in English and follow existing Sparlectra style.

- For **new features** (i.e. functionality that adds a new capability, not bugfixes or minor improvements):
  - Add at least one example program in `src/examples/`
  - Naming convention: `exp_<feature>.jl`
  - The example should demonstrate typical usage and expected behavior

- For **new features**, also ensure:
  - The existing documentation is reviewed and updated where appropriate
  - The feature is added to a suitable documentation page (theory, usage, or API section)
  - The feature matrix is checked and updated if applicable

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

# Codex task: Fix Julia 1.12 world-age and keyword-forwarding errors in example scripts

## Context

The repository currently shows Julia 1.12 / Revise-related warnings and a keyword mismatch in an example script, for example:

```text
WARNING: Detected access to binding `Main.bench_run_acpflow` in a world prior to its definition world.
Julia 1.12 has introduced more strict world age semantics for global bindings.
Hint: Add an appropriate `invokelatest` around the access to this binding.
```

and:

```text
MethodError: no method matching bench_run_acpflow(...; autodamp=..., autodamp_min=..., start_projection=..., ...)
got unsupported keyword arguments "autodamp", "autodamp_min", "start_projection",
"start_projection_try_dc_start", "start_projection_try_blend_scan",
"start_projection_blend_lambdas", "start_projection_dc_angle_limit_deg"
```

The affected example is at least:

```text
src/examples/matpower_import.jl
```

## Goal

Make the example script robust under Julia 1.12 and VS Code / Revise workflows, and fix the keyword mismatch so newly introduced solver/start-projection options are either correctly supported or intentionally filtered.

## Required work

1. Inspect the call chain in `src/examples/matpower_import.jl`
   - Identify where `main()` is defined and called.
   - Identify where `bench_run_acpflow` is defined and called.
   - Identify all wrappers between the parsed configuration/options and the final power-flow call.

2. Fix Julia 1.12 world-age warnings
   - Use `Base.invokelatest(...)` at script/developer entry boundaries.
   - Apply it to `main()` invocation if needed.
   - Apply it to dynamically accessed local entry functions such as `bench_run_acpflow` only at the problematic call site.
   - Do not wrap internal numerical kernels or hot loops.

3. Fix keyword forwarding
   - The caller currently passes keywords not accepted by `bench_run_acpflow`.
   - Decide whether these keywords are intended to be supported by this example:
     - `autodamp`
     - `autodamp_min`
     - `start_projection`
     - `start_projection_try_dc_start`
     - `start_projection_try_blend_scan`
     - `start_projection_blend_lambdas`
     - `start_projection_dc_angle_limit_deg`
   - If they are intended to be supported, add them explicitly to the `bench_run_acpflow` keyword signature and forward them to the downstream function.
   - If any are not meaningful in this path, filter them explicitly before calling `bench_run_acpflow` and add a short comment explaining why.
   - Do not leave accidental unsupported keyword forwarding in place.

4. Preserve current behavior
   - Existing behavior must remain unchanged when these options are not provided.
   - Defaults must match the existing solver defaults or current configuration defaults.
   - Avoid broad refactoring.

5. Add a small smoke/regression test if practical
   - Prefer a lightweight test that validates the example call path accepts the new keyword set.
   - The test does not need to run a large MATPOWER case.
   - If this example is not part of the formal test suite, add a minimal direct test for the wrapper function or document the manual test command.

6. Update changelog
   - Add a `Bugfix` entry in `CHANGELOG.md`.
   - Mention Julia 1.12 / Revise world-age handling and keyword-forwarding correction for the MATPOWER example path.

## Acceptance criteria

- Running the affected example no longer emits Julia 1.12 world-age warnings for `main` or `bench_run_acpflow`.
- The keyword mismatch is fixed.
- The following keywords are either accepted and forwarded correctly or intentionally filtered with a short comment:
  - `autodamp`
  - `autodamp_min`
  - `start_projection`
  - `start_projection_try_dc_start`
  - `start_projection_try_blend_scan`
  - `start_projection_blend_lambdas`
  - `start_projection_dc_angle_limit_deg`
- Existing tests still pass.
- The smallest relevant test or manual command is reported.
- If a test fails, report the exact command and relevant failure output.

## Suggested diagnostic command

```bash
julia --project=. --depwarn=error src/examples/matpower_import.jl
```

If the example requires arguments or a config file, use the smallest available MATPOWER case and the minimal options needed to exercise the failing path.

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