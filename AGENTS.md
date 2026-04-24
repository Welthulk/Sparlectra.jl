# AGENTS.md

This repository is a Julia project.

## Rules
- Use Julia for code changes, scripts, diagnostics, and tests.
- Always run Julia commands with `julia --project=.`
- Do not use Python to inspect or modify Julia source unless explicitly requested.
- Prefer minimal, targeted changes.

## Standard commands
- Verify toolchain:
  `which julia && julia --version`
- Instantiate environment:
  `julia --project=. -e 'using Pkg; Pkg.instantiate()'`
- Run tests:
  `julia --project=. test/runtests.jl`
- Alternative full test:
  `julia --project=. -e 'using Pkg; Pkg.test()'`

## Repo conventions
- Keep code comments and docstrings in English.
- Prefer existing Sparlectra naming and style.
- Before editing, inspect current tests that cover the affected code.
- After code changes, run the smallest relevant test first, then broader tests if needed.

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

- Tests alone are not sufficient for new features; documentation and examples are required