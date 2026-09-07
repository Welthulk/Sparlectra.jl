# CGMES importer and exporter (developer notes)

## Purpose

Lean CGMES (ENTSO-E CIM) import for Sparlectra plus the CGMES 2.4.15 RDF/XML
exporter. The importer is built as numbered layers; each layer consumes only
the layers below it.

## Layered include order

`CGMESImporter.jl` is the module root and includes in this order:

```julia
include("cgmes_schema.jl")             # shared namespace/profile tables
include("cgmes_container.jl")          # layer 1: folder/ZIP -> file list
include("cgmes_reader.jl")             # layer 2: RDF/XML property bags
include("cgmes_store.jl")              # layer 3: merged profile store
include("cgmes_import_analysis.jl")    # delivery analysis (no solve)
include("cgmes_voltage_inference.jl")  # nominal-voltage reconstruction
include("cgmes_topology.jl")           # layer 4: bus-branch topology from TP
include("cgmes_keys.jl")               # structural identity keys
include("cgmes_mapping.jl")            # layer 5: CIM objects -> Net
include("cgmes_report.jl")             # diagnostics and SV comparison
include("cgmes_export.jl")             # RDF/XML writer
include("cgmes_testsets.jl")           # conformity test-set access
```

## File responsibilities

| File | Responsibility | Key functions |
|---|---|---|
| `CGMESImporter.jl` | Module root, public import entry | `importCGMES` |
| `cgmes_schema.jl` | Namespaces, profile URIs, version detection | `cgmesVersionFromNamespace` |
| `cgmes_container.jl` | Turn folder/ZIP(s) into a profile file list | `collectCGMESFiles` |
| `cgmes_reader.jl` | Generic RDF/XML property-bag reader | `readCGMESFile!` |
| `cgmes_store.jl` | Merged profile store with by-class access | `loadCGMES`, `ref`, `num` |
| `cgmes_topology.jl` | Bus-branch topology from the TP profile | `buildTopology`, `busOfEquipment` |
| `cgmes_mapping.jl` | Map CIM objects onto the Sparlectra Net | `_attachTapControl!` and the per-class mappers |
| `cgmes_keys.jl` | Structural identity keys for `net.cgmes_ids` | `cgmesKeyPowerTransformer` |
| `cgmes_voltage_inference.jl` | Reconstruct missing nominal voltages | `_inferNominalVoltages` |
| `cgmes_report.jl` | Import summary, SV comparison, SC coverage | `summarizeCGMES`, `compareWithSV` |
| `cgmes_import_analysis.jl` | Why-does-this-not-import analysis | `analyzeCGMES` |
| `cgmes_export.jl` | CGMES 2.4.15 writer (EQ, TP, SSH, SV) | `writeCGMESFiles` |
| `cgmes_testsets.jl` | On-demand conformity test-set download/cache | `fetchCGMESTestSet` |

## Conventions

- The run path imports through `importCGMES(config; path, name)` (defined in
  `src/import/case_import.jl`), which applies the whole `cgmes_import`
  configuration in one place. Do not unpack configuration values at call
  sites.
- Placeholder guards (shunt capacity, tap band) warn and skip by default;
  `cgmes_import.placeholder_guards = :strict` turns them into errors.
- The `start_values` decision (`:auto`/`:sv`/`:flat`) is resolved at import
  time against the delivery actually read, not in configuration code.

## Reference

Rendered API documentation: `docs/src/reference_adapters.md`. User-facing
documentation: `docs/src/cgmes_import.md` and `docs/src/cgmes_export.md`.
