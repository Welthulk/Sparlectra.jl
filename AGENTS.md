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
