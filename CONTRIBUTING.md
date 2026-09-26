# Contributing to Sparlectra.jl

Sparlectra.jl is primarily a research and educational project with limited
maintainer capacity. To keep the codebase consistent and maintainable,
contributions follow the process below.

## Discuss first

Non-trivial changes require prior agreement. Open an Issue or Discussion,
describe the idea and wait for maintainer feedback before writing code.
Pull requests without prior agreement are usually not reviewed.

Trivial fixes (typos, broken links, obvious one-line bugs) can go directly
to a pull request.

## Workflow

1. Fork the repository.
2. Create a feature branch.
3. Implement the agreed change.
4. Open a pull request that references the Issue or Discussion.

For interactive development, `Revise.jl` from your global Julia environment
avoids restarting Julia after every edit. Do not add `Revise` to the project
dependencies. Where available, run example entry points via
`Base.invokelatest(...)` to avoid world-age issues.

## Code

Sparlectra prioritizes numerical robustness, deterministic behavior,
conceptual clarity and minimal complexity. In practice:

- clear structure and descriptive names
- minimal dependencies
- deterministic algorithms
- mathematical changes documented and justified
- no performance regressions
- comments and docstrings in English

Run artifacts are never committed from the repository root. Files such as
`ac_islands.csv`, `ac_island_solver_summary.csv`, `ac_island_<id>_solver.log`,
`matpower_dcline.csv`, `q_limit.log`, `performance.log`, `run.log` and
`effective_config.yaml` belong in the run output directory or in a
test-owned temporary directory.

Documentation headings with a Documenter label (`## [Text](@id page-slug)`)
are referenced by the Web UI help. Never change or drop an `@id`. After
editing a referenced section, regenerate the help excerpts and commit them.
Details: [DEVELOPER.md](DEVELOPER.md#documentation-anchors-are-web-ui-contracts).

## Tests

Every functional change includes dedicated unit tests in the repository.
Tests must be deterministic and cover relevant edge cases. Critical
numerical paths must be tested explicitly.

New code should not significantly reduce overall coverage. Coverage reports
are welcome but optional.

Document how to run the new tests, the expected results and any assumptions
or limitations in the test file or the documentation, not only in the pull
request description.

## Benchmarks

Changes to solver behavior, numerical routines or performance-critical code
require a reproducible before/after benchmark that the maintainer can run:

- test case (network size, type)
- runtime
- iteration count
- convergence behavior

## Acceptance

A pull request may be rejected for missing prior discussion, missing tests,
missing benchmarks where required, insufficient documentation, or added
complexity without clear benefit.

Reviews are best-effort. Response times vary, and not every pull request
will be accepted.

## Developer reference

Internal architecture notes (documentation anchor contract, Web UI job
lifecycle, outer-loop controller interface) are in
[DEVELOPER.md](DEVELOPER.md).
