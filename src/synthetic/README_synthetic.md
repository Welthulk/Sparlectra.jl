# Synthetic grids (developer notes)

## Purpose

Deterministic synthetic tiled-grid network builder for tests, benchmarks, and
performance work: reproducible networks of arbitrary size without external
case files.

## File responsibilities

| File | Responsibility | Key functions |
|---|---|---|
| `synthetic_grids.jl` | Tiled-grid construction and bus indexing | `build_synthetic_tiled_grid_net`, `synthetic_tiled_grid_bus_index` |

## Conventions

- Determinism is the point: same size in, same network out. Do not add
  randomness without a seed parameter.
- Keep generated fixtures out of version control; tests build them on the
  fly.

## Reference

Rendered API documentation: `docs/src/reference_synthetic.md`.
