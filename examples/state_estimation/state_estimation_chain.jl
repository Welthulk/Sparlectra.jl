# Copyright 2023-2026 Udo Schmitz
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# Date: 2026-08-26
# file: examples/state_estimation/state_estimation_chain.jl
# purpose: SE phase 5 walkthrough of the measurement-to-N-1 chain, program
#          side: write/read a measurement CSV v1 (lossless roundtrip), run
#          the SE, then start the power flow from the estimate with both
#          modes (se_state: model authoritative; se_snapshot: balance
#          takeover, immediate convergence, slack pickup at tolerance). The
#          same chain drives the Web UI /stateestimation page via
#          se_state.csv artifacts.

using Sparlectra
using Printf

include(joinpath(@__DIR__, "..", "others", "example_header.jl"))

function _create_demo_net()
  net = Net(name = "se_chain_demo", baseMVA = 100.0)
  for b in ("Slack", "LoadA", "LoadB")
    addBus!(net = net, busName = b, vn_kV = 110.0)
  end
  addACLine!(net = net, fromBus = "Slack", toBus = "LoadA", length = 12.0, r = 0.02, x = 0.2, c_nf_per_km = 0.0, tanδ = 0.0)
  addACLine!(net = net, fromBus = "LoadA", toBus = "LoadB", length = 9.0, r = 0.02, x = 0.2, c_nf_per_km = 0.0, tanδ = 0.0)
  addACLine!(net = net, fromBus = "Slack", toBus = "LoadB", length = 11.0, r = 0.02, x = 0.2, c_nf_per_km = 0.0, tanδ = 0.0)
  addProsumer!(net = net, busName = "Slack", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.02, va_deg = 0.0, referencePri = "Slack")
  addProsumer!(net = net, busName = "LoadA", type = "ENERGYCONSUMER", p = 35.0, q = 12.0)
  addProsumer!(net = net, busName = "LoadB", type = "ENERGYCONSUMER", p = 28.0, q = 9.0)
  ok, msg = validate!(net = net)
  ok || error("demo net invalid: $msg")
  return net
end

function run_state_estimation_chain()
  print_example_banner("examples/state_estimation/state_estimation_chain.jl", "measurement CSV roundtrip, SE, and the PF started from the estimate (se_state / se_snapshot)")

  net = _create_demo_net()
  iteFlat, erg = runpf!(net, 40, 1e-10, 0; method = :rectangular, opt_flatstart = true)
  erg == 0 || error("Power flow did not converge")
  @printf("flat-start PF: %d iterations\n\n", iteFlat)

  # 1) measurement CSV v1: write the synthetic set, wipe, read it back
  dir = mktempdir()
  mfile = joinpath(dir, "demo.measurements.csv")
  append!(net.measurements, generateMeasurementsFromPF(net; noise = false))
  w = writeMeasurementsCSV(net; file = mfile)
  empty!(net.measurements)
  r = readMeasurementsCSV!(net; file = mfile)
  @printf("measurement CSV v1: wrote %d rows, read %d back (atomic, line-precise errors on bad files)\n\n", w.count, r.total)

  # 2) the estimation registers the chain start state
  res = runse!(net; maxIte = 30, tol = 1e-10, updateNet = true)
  @printf("SE: converged in %d iteration(s), J = %.3e (dof %d)\n\n", res.iterations, res.objectiveJ, res.dof)

  # 3) PF from the estimate, both modes
  r1 = runpf_from_se!(net, 40, 1e-10, 0; mode = :se_state, method = :rectangular)
  @printf("se_state    : %d iteration(s) (flat start needed %d); model injections stay authoritative\n", r1.iterations, iteFlat)
  runse!(net; maxIte = 30, tol = 1e-10, updateNet = true)
  r2 = runpf_from_se!(net, 40, 1e-10, 0; mode = :se_snapshot, method = :rectangular)
  @printf("se_snapshot : %d iteration(s), slack pickup %.2e MW; persistent loads untouched\n\n", r2.iterations, abs(r2.slack_pickup_mw))

  # 4) the persistence half the Web UI chain uses (se_state.csv artifact)
  sfile = joinpath(dir, "se_state.csv")
  runse!(net; maxIte = 30, tol = 1e-10, updateNet = true)
  writeSEStateCSV(net; file = sfile)
  net2 = _create_demo_net()
  readSEStateCSV!(net2; file = sfile)
  r3 = runpf_from_se!(net2, 40, 1e-10, 0; mode = :se_snapshot, method = :rectangular)
  @printf("chained via se_state.csv on a fresh net: %d iteration(s)\n", r3.iterations)
  println("(the Web UI runs this chain from /stateestimation: SE run -> \"run power flow from this estimate\" -> N-1)")
  return nothing
end

run_example(run_state_estimation_chain)
