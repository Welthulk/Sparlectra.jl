# DTF importer (developer notes)

## Purpose

Native importer for legacy DTF (FOR001) network files: parsing the section
cards, building a `Net`, and outage-list support. Experimental/internal
format; ambiguous `.DAT` files are only read when the FOR001 markers are
present and `case_format = :dtf_for001` is requested.

## File responsibilities

| File | Responsibility | Key functions |
|---|---|---|
| `DTFImporter.jl` | The whole module: parser, net construction, outage handling | `read_dtf`, `build_net`, `apply_single_branch_outage!` |

## Conventions

- The run path constructs through `import_case` (`src/import/case_import.jl`),
  which rejects DC-line-like content before parsing and stamps the configured
  network parameters afterwards.
- Tap conversion and shunt modeling follow the configured
  `transformer.tap_changer_model` and `matpower.bus_shunt_model`; the format
  itself carries no such settings.

## Reference

Rendered API documentation: `docs/src/reference_adapters.md`. Format notes:
`docs/src/dtf_format.md`.
