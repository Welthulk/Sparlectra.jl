# MATPOWER adapter (developer notes)

## Purpose

Everything MATPOWER: reading `.m` and generated `.jl` case files, building a
`Net` from parsed case data, exporting a `Net` back to a case file, the
config-driven case runner, and on-demand download/caching of published
benchmark cases.

## File responsibilities

| File | Responsibility | Key functions |
|---|---|---|
| `MatpowerIO.jl` | Case parsing into `MatpowerCase`, convention normalization, dcline diagnostics | `read_case`, `legacy_sort_bus` |
| `createnet_powermat.jl` | Build a `Net` from case data | `createNetFromMatPowerCase`, `createNetFromMatPowerFile` |
| `exportMatPower.jl` | Export a `Net` as a MATPOWER case file | `writeMatpowerCasefile` |
| `matpower_runner.jl` | Config-driven case runner with performance profiling | `run_matpower_case` |
| `FetchMatpowerCase.jl` | Download and locally cache published benchmark cases | `fetch_matpower_case`, `ensure_matpower_case` |

## Conventions

- Parsing and construction are separate on purpose: `MatpowerIO.read_case`
  yields the case object, `createNetFromMatPowerCase` builds the network.
  `createNetFromMatPowerFile` is the exported one-call form.
- Convention switches (`matpower_shift_sign`, `matpower_shift_unit`,
  `matpower_ratio`, `tap_changer_model`) exist because published cases do not
  agree; a per-case sidecar can pin them. Fewer Newton iterations after a
  convention change is the practical sign the convention is right.
- The roundtrip marker `mpc.sparlectra.tap_changer_model` prevents reapplying
  the tap-impedance correction on reimport of our own exports.
- The Web UI owns benchmark case resolution and caching; benchmark work goes
  through the service path instead of a second download workflow.

## Reference

Rendered API documentation: `docs/src/reference_adapters.md`. User-facing
documentation: `docs/src/matpower_format.md` and `docs/src/matpower_import.md`.
