# Network Reports (ACPFlowReport)

`buildACPFlowReport` returns a machine-readable report object next to the
formatted terminal output of `printACPFlowResults`.

## Build a report

```julia
using Sparlectra

net = Net(name = "report_demo", baseMVA = 100.0)

addBus!(net = net, busName = "B1", vn_kV = 110.0)
addBus!(net = net, busName = "B2", vn_kV = 110.0)
addBus!(net = net, busName = "B3", vn_kV = 110.0)
addBus!(net = net, busName = "B4", vn_kV = 110.0)

# line totals in ohm and siemens, or values per km and the length
addACLine!(net = net, fromBus = "B1", toBus = "B2", r_ohm = 0.6, x_ohm = 1.95, b_s = 1.5e-5, length_km = 5.0)
addACLinePerKm!(net = net, fromBus = "B3", toBus = "B4", length_km = 20.0, r_ohm_per_km = 0.12, x_ohm_per_km = 0.39, c_nf_per_km = 9.5)

# Optional bus-link connection (reported in report.links)
addLink!(net = net, fromBus = "B2", toBus = "B3", status = 1)

addProsumer!(net = net, busName = "B1", type = "EXTERNALNETWORKINJECTION", referencePri = "B1")
addProsumer!(net = net, busName = "B2", type = "GENERATOR", p = 10.0, q = 1.0)
addProsumer!(net = net, busName = "B3", type = "ENERGYCONSUMER", p = 15.0, q = 5.0)
addProsumer!(net = net, busName = "B4", type = "ENERGYCONSUMER", p = 25.0, q = 10.0)

cfg = SparlectraConfig(powerflow = PowerFlowConfig(max_iter = 40, tol = 1e-10))
result = run_sparlectra(net = net, config = cfg)
ite, etime = result.iterations, result.elapsed_s

report = buildACPFlowReport(
  net;
  ct = etime,
  ite = ite,
  tol = 1e-10,
  converged = result.final_converged,
  solver = :rectangular,
)
```

The snippet is a complete script for a `julia --project=.` session.

## Report content

`ACPFlowReport` has the fields `metadata`, `nodes`, `branches`, `links`,
`transformer_controls`, `q_limit_events` and `hvdc_links`:

```julia
report.metadata.total_p_loss_MW
report.nodes[1]
report.branches[1]
report.links
```

## DataFrame conversion

Each vector field converts directly:

```julia
using DataFrames

nodes_df = DataFrame(report.nodes)
branches_df = DataFrame(report.branches)
links_df = DataFrame(report.links)
```
