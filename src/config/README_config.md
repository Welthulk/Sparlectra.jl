# Configuration (developer notes)

## Purpose

The typed configuration surface (`SparlectraConfig` and all `*Config` structs)
plus the dependency-free YAML subset parser. The shipped default configuration
lives here as `configuration.yaml.example`.

## File responsibilities

| File | Responsibility | Key functions |
|---|---|---|
| `configuration.jl` | All config structs, loading, validation, accessors for the active configuration | `load_sparlectra_config`, `set_sparlectra_config!`, `powerflow_config` |
| `yamlparams.jl` | Minimal YAML parser used for configuration and benchmark params | `load_yaml_dict`, `merge_yaml_dict!`, `parse_yaml_scalar` |
| `configuration.yaml.example` | The annotated default configuration; `DEFAULT_SPARLECTRA_CONFIG_PATH` points here | (data) |

## Conventions

- The YAML parser is a deliberate subset: block style only. `key: {}` parses
  as the string `"{}"`, so dict-valued keys need block style plus placeholder
  normalization.
- New configuration keys need: struct field, parser wiring, validation,
  documentation in `docs/src/configuration.md`, and (if GUI-editable) an entry
  in `GUI_EDITABLE_CONFIG_KEYS` (`src/api/config_overrides.jl`).
- Loaded templates are never mutated by runs; run paths copy with
  `_copy_sparlectra_with_powerflow` and friends (`src/acpflow/import_context.jl`).

## Reference

Rendered API documentation: `docs/src/reference_config.md`. User-facing
documentation: `docs/src/configuration.md`.
