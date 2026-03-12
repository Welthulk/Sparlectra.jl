# Network Reports (ACPFlowReport)

Besides the formatted terminal output (`printACPFlowResults`), Sparlectra provides a
machine-readable report object for downstream processing.

## Build a report

```julia
using Sparlectra

net = Net(name = "report_demo", baseMVA = 100.0)

addBus!(net = net, busName = "B1", busType = "SLACK", vn_kV = 110.0)
addBus!(net = net, busName = "B2", busType = "PQ", vn_kV = 110.0)
addBus!(net = net, busName = "B3", busType = "PQ", vn_kV = 110.0)
addBus!(net = net, busName = "B4", busType = "PQ", vn_kV = 110.0)

addACLine!(net = net, fromBus = "B1", toBus = "B2", length = 5.0, r = 0.05, x = 0.5)
addACLine!(net = net, fromBus = "B3", toBus = "B4", length = 20.0, r = 0.05, x = 0.5)

# Optional bus-link connection (reported in report.links)
addLink!(net = net, fromBus = "B2", toBus = "B3", status = 1)

addProsumer!(net = net, busName = "B1", type = "EXTERNALNETWORKINJECTION", referencePri = "B1")
addProsumer!(net = net, busName = "B2", type = "GENERATOR", p = 10.0, q = 1.0)
addProsumer!(net = net, busName = "B3", type = "ENERGYCONSUMER", p = 15.0, q = 5.0)
addProsumer!(net = net, busName = "B4", type = "ENERGYCONSUMER", p = 25.0, q = 10.0)

ite, erg, etime = run_net_acpflow(
  net = net,
  max_ite = 40,
  tol = 1e-10,
  method = :polar_full,
  opt_sparse = true,
  show_results = false,
)

report = buildACPFlowReport(
  net;
  ct = etime,
  ite = ite,
  tol = 1e-10,
  converged = (erg == 0),
  solver = :polar_full,
)
```

## Report content

`ACPFlowReport` contains:

- `metadata`
- `nodes`
- `branches`
- `links`
- `q_limit_events`

Example access:

```julia
report.metadata.total_p_loss_MW
report.nodes[1]
report.branches[1]
report.links
```

## DataFrame conversion

Each table-like vector can be converted directly:

```julia
using DataFrames

nodes_df = DataFrame(report.nodes)
branches_df = DataFrame(report.branches)
links_df = DataFrame(report.links)
```

## Full runnable example

See `src/examples/using_netreports.jl` for a complete script.
