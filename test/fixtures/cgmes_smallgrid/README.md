# ENTSO-E SmallGrid (CGMES 2.4.15, NodeBreaker BaseCase)

Unchanged copies of two files of the ENTSO-E CGMES conformity test
configurations (package CASv2.0, `SmallGrid/NodeBreaker`):

- `CGMES_v2.4.15_SmallGridTestConfiguration_BaseCase_Complete_v3.0.0.zip`
- `CGMES_v2.4.15_SmallGridTestConfiguration_Boundary_v3.0.0.zip`

The test suite uses them as the network of the synthetic measurement
fixture `test/fixtures/measurements/smallgrid_synthetic.table.csv`
(built by `tools/gen_smallgrid_measurements.jl`), in the CGMES importer
tests (`element_by_mrid`, import from a loaded store) and in the snapshot
series tests (`test/test_cgmes_series.jl`, and `test/test_cgmes_ssh_series.jl`,
which builds a 96-step series with `tools/cgmes_ssh_series.jl`). This is the
one ENTSO-E package checked into the repository.
