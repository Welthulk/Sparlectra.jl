# Sparlectra Case Format (developer notes)

## Purpose

The Sparlectra Case Format (SCF): a self-describing JSON case format
(`<case>.scf.json`) that carries the network, its operating point, case-scoped
configuration, controllers, measurements, and study definitions (contingency,
short circuit). Reader, writer, and the JSON layer live here.

## Include order

`scf.jl` is the include hub:

```julia
include("scf_json.jl")
include("scf_controllers.jl")
include("scf_export.jl")
include("scf_import.jl")
```

## File responsibilities

| File | Responsibility | Key functions |
|---|---|---|
| `scf.jl` | Include hub | (includes only) |
| `scf_json.jl` | Deterministic pretty JSON writer and minimal parser | `scf_json_string`, `scf_number_or_sentinel` |
| `scf_export.jl` | Serialize a `Net` into the format | `exportSCF`, `net_to_scf`, `scf_is_case_config_key` |
| `scf_import.jl` | Read, validate, and build a `Net`; run-feeding accessors | `importSCF`, `scf_to_net`, `scf_case_config`, `scf_case_studies` |
| `scf_controllers.jl` | Controllers (FACTS and regulation) in the format | `_scf_apply_controllers!` |

## Conventions

- The reader enforces the current format revision exactly (deliberate
  pre-release rule, stated in the specification). ALL run-feeding readers
  check the revision themselves; only display-only readers stay lenient.
- The case-scope rule for configuration keys has one source in code,
  `scf_is_case_config_key` (derived from the GUI allow-list); keep the prose
  in the specification in step when it changes.
- Roundtrip fidelity is guarded by a reflection comparison over all fields
  (`test/test_scf.jl`, `scf_roundtrip_field_diffs`) with an exception list
  that names reasons. When a roundtrip question comes up, run that comparison
  first.
- The writer is deterministic (stable key order, stable number formatting);
  diffs of exported files are meaningful and tests rely on that.

## Reference

Rendered API documentation: `docs/src/reference_adapters.md`. Format
specification: `docs/src/scf.md`.
