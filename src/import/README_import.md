# Case import entry point (developer notes)

## Purpose

The one import entry point of the run path: format detection, per-format
construction, and the `ImportedCase` record every service consumes. Format
dispatch used to exist five times across the services; each copy applied a
different subset of the configuration. This directory replaces that with a
single path.

## File responsibilities

| File | Responsibility | Key functions |
|---|---|---|
| `case_import.jl` | `ImportedCase`, format detection, per-format construction, the configured CGMES import | `import_case`, `_detect_case_format`, `importCGMES(config; ...)` |

## Conventions

- Every service imports through `import_case(path, config; run_kind, ...)`.
  `run_kind` matters: the CGMES `hvdc_mode` applies only to `:powerflow`
  runs, every other kind gets plain injections.
- Format policy stays with the caller: a service that does not accept a
  format checks `ImportedCase.format` and words its own refusal.
- The CGMES `start_values` decision resolves here because `auto` can only be
  answered against the delivery actually read.
- Provenance keys (import records, path lists, decisions) feed the service
  artifacts; extend the dictionary rather than adding side channels.

## Reference

Rendered API documentation: `docs/src/reference_adapters.md`.
