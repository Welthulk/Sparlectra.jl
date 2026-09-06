# Build tools (developer notes)

## Purpose

Fast-start build helpers: the PrecompileTools workload that covers the
MATPOWER import to power-flow path, and the one-call sysimage/app builders.

## File responsibilities

| File | Responsibility | Key functions |
|---|---|---|
| `precompile.jl` | Precompile workload (import -> rectangular/DC solve) so first interactive runs skip most JIT compilation | (workload block) |
| `sysimage_builder.jl` | One-call API around the fast-start sysimage build | `buildSysimage` |

## Conventions

- Sysimage builds are long; they are run on request, not as a side effect.
- The Web UI reads sysimage state through `src/webui/sysimage.jl`; keep the
  metadata contract (manifest path, status fields) in sync with it.
- These files have no reference page; `buildSysimage` renders on
  `docs/src/reference_api.md` (Application entry points).
