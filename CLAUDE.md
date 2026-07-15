# CLAUDE.md

Guidance for Claude Code when working in this repository. Sparlectra is a Julia
package for power-flow calculation via Newton-Raphson (see the top-level
`README.md` for the functional overview). The project root is the Julia
project (`julia --project=.`).

## Before editing

- Use Julia for all code changes, scripts, diagnostics, and tests in this repo. Don't reach for Python to inspect or modify Julia source unless explicitly asked.
- Check repository freshness first:
  ```bash
  git status --short
  git rev-parse --short HEAD
  git log --oneline -n 5
  git fetch --all --prune && git status -sb
  ```
  If the working tree is dirty or the local branch is behind upstream, surface that before editing instead of silently continuing.
- Reproduce bugs from the current `HEAD` — don't rely on pasted logs or patches from an older commit. Verify a patch still applies cleanly before using it.
- If a file has merge-conflict markers (`<<<<<<<`, `=======`, `>>>>>>>`), stop and report instead of guessing a resolution.
- When restoring or repairing code from history, restore the whole coherent call chain (signature, all callers, all forwarded keywords) rather than pasting isolated snippets.

## Protected core solver files

These contain central solver, configuration, or execution-path logic — treat changes here with extra care:

```text
src/powerflow_rectangular/rectangular_network_solver.jl
src/run_acpflow.jl
src/control_framework.jl
src/configuration.jl
src/Sparlectra.jl
```

Rules for these files:
- No broad autonomous refactoring unless explicitly requested; don't mix refactoring with feature/bugfix work.
- Don't change public entry points (`runpf!`, `runpf_rectangular!`, `run_acpflow`) or keyword signatures without checking every in-repository caller and the full forwarding chain.
- Don't change solver status metadata without checking all readers, tests, and reporting code.
- Don't add compatibility shims, wrapper layers, or duplicate APIs unless explicitly requested.
- Before touching these files, confirm the package still loads and the key bindings still resolve:
  ```bash
  julia --project=. -e "using Sparlectra"
  julia --project=. -e "using Sparlectra; for s in Symbol.([\"runpf!\", \"runpf_rectangular!\", \"run_complex_nr_rectangular\"]); @assert isdefined(Sparlectra, s) s; println(s, \" => \", length(methods(getfield(Sparlectra, s))), \" method(s)\"); end"
  ```
- If a protected file is broken, restore package loadability first — don't start feature work or refactor while `using Sparlectra` fails. Prefer a full restore from the last known-good commit over incremental patching of corrupted code.

## Actions requiring explicit request

Don't do these unless the user explicitly asks for them:
- active wrong-branch rescue/retry loops, Monte-Carlo or multi-start solver scans, continuation/homotopy start logic, new DC-start or QV-start fallback modes;
- large-scale rewrites of `rectangular_network_solver.jl`, or broad solver restructuring mixed into a bugfix;
- changing public APIs to hide test failures, or converting public calls to qualified internal calls just to dodge a failure;
- adding global variables or dummy functions to paper over a missing keyword — fix the actual method signature and forwarding chain instead;
- a second, independent MATPOWER/Web UI benchmark download/cache workflow (the Web UI already owns case registry/resolution/caching — reuse it via `start_powerflow_run` / `get_powerflow_result` / the existing cache path);
- creating a new changelog file (the changelog is `docs/src/changelog.md` — see below).

## Julia 1.12 world-age / Revise handling

Julia 1.12 warns when a script entry point calls a just-defined global binding from an older world age:
```text
WARNING: Detected access to binding `Main.main` in a world prior to its definition world.
```
Treat this as a future error, not noise, in touched code.
- Call script/example/test entry points through `Base.invokelatest`, e.g. `Base.invokelatest(runner)` or `Base.invokelatest(main, args...)`.
- Only wrap entry-point boundaries this way — don't wrap internal hot-loop functions.
- Reproduce with `julia --project=. --depwarn=error path/to/script.jl` to get a stack trace when investigating.
- In example/test files included into `Main`, qualify Sparlectra enum constants/types (`Sparlectra.Slack`, `Sparlectra.PV`, `Sparlectra.PQ`) unless explicitly imported.
- When a wrapper forwards a large keyword set (common in example/benchmark scripts) and a new keyword is added, update every wrapper in the call chain, not just the innermost function.

## Standard commands

```bash
which julia && julia --version
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. -e 'using Sparlectra'                    # package load / syntax check
julia --project=. test/runtests.jl                          # fast profile (default)
SPARLECTRA_TEST_PROFILE=extended julia --project=. test/runtests.jl   # extended profile
julia --project=. test/runtests.jl all                      # fast + extended
julia --project=docs docs/make.jl                            # docs build
```

Verification depth should match the change: a package-load check for anything touching module structure/exports/signatures/config, the fast profile for regular changes, and fast+extended+docs before calling a branch merge-ready.

Docs build notes:
- `docs/make.jl` is the acceptance gate for anything touching documentation, Documenter config, public docstrings, or autodocs-referenced files. If it reaches Documenter and fails (unresolved `@ref`, missing docstring, doctest failure, invalid Markdown/math, duplicate objects), the task isn't done.
- If it can't even reach Documenter due to a network/proxy/registry failure, say so explicitly rather than claiming the docs build passed.
- A misleading `Base.-` docstring error after adding an `@ref` link usually means a hyphenated anchor label (e.g. `{#reference-api}`) — use underscores (`{#reference_api}`) instead.

## Testing conventions

- Extend or consolidate existing tests before adding a new test file; avoid near-duplicate tests that only differ in incidental parameters (prefer a table-driven test or loop inside one testset).
- Keep test fixtures small — no large MATPOWER cases or long-running examples unless explicitly required.
- Remove tests that only protect deleted/obsolete compatibility behavior when doing a refactor.
- When adding/removing/renaming/moving tests or profile groups, update `docs/src/tests.md` in the same change (profiles, commands, files per profile, what each group verifies).
- Keep at least one regression test per fixed bug, but don't stack multiple overlapping regressions for the same failure mode.

## Code hygiene

- Prefer minimal, targeted changes; keep edits local to the affected implementation path.
- Don't silently change public APIs or example-script behavior in bugfix/improvement tasks. For explicit refactoring/breaking-change tasks, prefer the clean target design over compatibility layers — update all in-repo callers, examples, docs, and tests instead of preserving the old path.
- Remove dead code, stale debug branches, and duplicated helper logic introduced during debugging; prefer one well-named helper over repeated blocks.
- Keep logging/diagnostics compact by default; verbose output is opt-in.
- Julia indexing style — avoid `IndexFromLength` lint warnings:
  - Prefer `for i in axes(A, 1)` over `for i = 1:size(A,1)`.
  - Prefer `for i in eachindex(A)` when iterating over an array.

### CSS / stylesheet formatting

Keep `.css` files human-readable — one selector block per block, one property per line, consistent indentation (e.g. `src/webui/static/sparlectra.css`). Don't collapse or minify existing readable stylesheets, and don't hide large CSS changes in generated/minified output.

## File headers, docstrings, comments

- New source/script/test/doc-support files use the repository's existing license-header style (don't invent new license text).
- Add docstrings (English) for public/semi-public APIs, exported bindings, and anything referenced by Documenter `@ref` or `@autodocs` — purpose, key args, return values, side effects, failure behavior.
- Inline comments (English) are for non-obvious choices: invariants, safety checks, path-containment rules, cache-vs-source distinctions, subtle numerical or Web UI/service decisions. Don't restate the next line of code.
- Important public exports in `src/Sparlectra.jl` should carry short inline/group comments explaining their role at the module-interface level; update these comments when exports are added/removed/renamed.

## Feature / Improvement / Bugfix classification

Used to decide what's required (docstrings, tests, examples, docs, changelog entry):

| Class | Typical indicators | Requirements |
|---|---|---|
| **Feature** | new public API, new modeling elements (FACTS, tap/phase-shifting transformers, controllers), new solver modes, new measurement/estimation capability, new import/export capability | mandatory: docstrings, tests, `examples/exp_<feature>.jl`, docs update, feature-matrix update, test-hygiene review |
| **Improvement** | perf/numerical-stability work, refactor without API change, better logging/error messages, internal restructuring | tests/docstrings updated only if behavior/usage changes; no example required; remove obsolete tests if it replaces prior diagnostic-only behavior |
| **Bugfix** | wrong numerical results, sign/indexing errors, crashes on valid input, example-script keyword mismatches, world-age issues | add/extend a regression test (prefer extending an existing testset), keep the change minimal, docs update only if prior docs were wrong |

If a change both fixes a bug and adds a capability, classify as feature. If unsure, default to feature (fuller documentation/example bar).

Rule of thumb: new way to model/solve/estimate/diagnose/import/export/interpret a network → feature. Existing way made faster/cleaner/safer → improvement. Corrects something that was wrong → bugfix.

## Changelog

- The changelog lives at **`docs/src/changelog.md`** — there is no root `CHANGELOG.md`. Add entries there under the right version heading; don't create a new changelog file or add release notes to the README or other doc pages.
- Format per entry: **Feature** / **Improvement** / **Bugfix** — brief description, purpose/impact, and anything users need to know.

## Documenter / Markdown / math conventions

- Every `@ref` target must resolve when running `docs/make.jl`. If a binding lives in a submodule, include that submodule in the relevant `@autodocs` block (e.g. `Sparlectra.FetchMatpowerCase`, not just `Sparlectra`). Intentionally-internal bindings should use plain code markup, not public `@ref` links.
- Avoid hyphenated custom anchor IDs (`{#reference-api}` / `(@ref reference-api)`) — Documenter can parse the hyphen as subtraction and report a misleading `Base.-` docstring error. Use underscores (`{#reference_api}`, `[API](@ref reference_api)`).
- Math: use `$...$` for inline math, ```` ```math ```` blocks only for real equations (not lone symbols). Don't start a list item with a bare `$...$` (renders as centered display math) — put descriptive text before the symbol.
- Use backticks for code/identifiers (`y_ser`, `addACLine!`), `$...$` only for actual math notation, and valid LaTeX (`$Y_{ii}$`, not `$\Y_{ii}$`).

## Confidentiality

Don't use or introduce the term `APSLF` in anything public-facing (docs, examples, comments, changelog, commit messages, generated text) — it's a private/internal project designation. Prefer neutral wording like "analytical power-series solver" or "series-based load-flow approach". Existing historic references may stay unless a task explicitly includes renaming/cleanup of them. When unsure whether text is public-facing, treat it as public-facing.

## Web UI / MATPOWER case cache

The Web UI already owns MATPOWER case registry, resolution, and local caching. Don't build a second, independent download/cache workflow — benchmark through the existing service path (`start_powerflow_run` → `get_powerflow_result` → `list_powerflow_artifacts` → `resolve_powerflow_artifact`), reuse the same case names/cache directory, and report clearly when a case isn't available through that mechanism rather than fetching it another way. Don't commit downloaded MATPOWER benchmark cases or generated large `.jl` case files.

For Web UI phase instrumentation: keep the operation log high-level (not per-iteration/per-linear-solve); put detailed repeated timings in `performance.log`/`run.log`/structured result metadata instead.

## Developer workflow notes (examples/)

- Example programs in `examples/` are primarily written for interactive (VS Code) use; prefer calling their entry functions via `Base.invokelatest(...)` over relying solely on `if abspath(PROGRAM_FILE) == @__FILE__`.
- `Revise.jl` is recommended for interactive contributor workflows but is intentionally kept out of project dependencies.
