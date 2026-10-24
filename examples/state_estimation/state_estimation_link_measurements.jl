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

# Date: 2026-10-23
# file: examples/state_estimation/state_estimation_link_measurements.jl
# purpose: link (bus coupler) flow measurements as WLS rows: a single
#          coupler whose P meter is 100 MW off the balance flow gets a
#          normalized residual, the bad-data loop eliminates it and the flow
#          stays on the balance value; on a coupler ring the measured links
#          give a circulation state that the flows after the estimate carry.

using Sparlectra
using Printf

include(joinpath(@__DIR__, "..", "others", "example_header.jl"))

const _LINK_DEMO_STD = measurementStdDevs(vm = 2e-3, pinj = 1.0, qinj = 1.0, pflow = 1.0, qflow = 1.0)

# S feeds B1 over a line, the coupler B1-B1a joins B1a, a line runs on to B2
function _single_coupler_net()
  net = Net(name = "se_link_single", baseMVA = 100.0)
  for b in ("S", "B1", "B1a", "B2")
    addBus!(net = net, busName = b, vn_kV = 110.0)
  end
  addPIModelACLine!(net = net, fromBus = "S", toBus = "B1", r_pu = 0.02, x_pu = 0.2, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "B1a", toBus = "B2", r_pu = 0.02, x_pu = 0.2, b_pu = 0.0, status = 1)
  addProsumer!(net = net, busName = "S", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.02, va_deg = 0.0, referencePri = "S")
  addProsumer!(net = net, busName = "B1", type = "ENERGYCONSUMER", p = 10.0, q = 3.0)
  addProsumer!(net = net, busName = "B1a", type = "ENERGYCONSUMER", p = 5.0, q = 2.0)
  addProsumer!(net = net, busName = "B2", type = "ENERGYCONSUMER", p = 20.0, q = 6.0)
  addLink!(net = net, fromBus = "B1", toBus = "B1a", status = 1)
  ok, msg = validate!(net = net)
  ok || error("demo net invalid: $msg")
  return net
end

# three couplers R1 -> R2 -> R3 -> R1 form a ring; lines feed R1 and load R2, R3
function _coupler_ring_net()
  net = Net(name = "se_link_ring", baseMVA = 100.0)
  for b in ("S", "R1", "R2", "R3", "L2", "L3")
    addBus!(net = net, busName = b, vn_kV = 110.0)
  end
  addPIModelACLine!(net = net, fromBus = "S", toBus = "R1", r_pu = 0.02, x_pu = 0.2, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "R2", toBus = "L2", r_pu = 0.02, x_pu = 0.2, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "R3", toBus = "L3", r_pu = 0.02, x_pu = 0.2, b_pu = 0.0, status = 1)
  addProsumer!(net = net, busName = "S", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.02, va_deg = 0.0, referencePri = "S")
  addProsumer!(net = net, busName = "L2", type = "ENERGYCONSUMER", p = 12.0, q = 4.0)
  addProsumer!(net = net, busName = "L3", type = "ENERGYCONSUMER", p = 18.0, q = 5.0)
  addLink!(net = net, fromBus = "R1", toBus = "R2", status = 1)
  addLink!(net = net, fromBus = "R2", toBus = "R3", status = 1)
  addLink!(net = net, fromBus = "R3", toBus = "R1", status = 1)
  ok, msg = validate!(net = net)
  ok || error("demo net invalid: $msg")
  return net
end

function _estimate!(net, meas)
  empty!(net.measurements)
  append!(net.measurements, meas)
  return with_state_estimation_config(max_iter = 30, tol = 1e-6, update_net = true) do
    runse!(net, meas)
  end
end

function run_state_estimation_link_measurements()
  print_example_banner("examples/state_estimation/state_estimation_link_measurements.jl", "link flow measurements in the WLS: gross error on a single coupler, elimination, ring circulation")

  # 1) single coupler: the meter reads 100 MW below the balance flow
  net = _single_coupler_net()
  runpf!(net, 40, 1e-10, 0; method = :rectangular)
  calcNetLosses!(net)
  calcLinkFlowsKCL!(net)
  fBal = net.linkVec[1].pFlow_MW
  meas = generateMeasurementsFromPF(net; noise = false, stddev = _LINK_DEMO_STD)
  addPflowMeasurement!(meas; net = net, value = fBal - 100.0, sigma = 2.0, linkNr = 1, id = "P_coupler")
  @printf("balance flow of the coupler: %.3f MW, meter: %.3f MW (sigma 2 MW)\n", fBal, fBal - 100.0)

  diag = with_state_estimation_config(max_iter = 30, tol = 1e-6) do
    runse_diagnostics(net, meas)
  end
  top = first(diag.diagnostics.measurement_ranking)
  @printf("largest normalized residual: %s, rn = %.1f, wii = %.2f\n", top.id, top.normalized_residual, top.wii)
  for e in diag.eliminations
    @printf("eliminated: %s (rn %.1f before)\n", e.id, e.normalized_residual_before)
  end

  # the final estimate on the cleaned set, as the SE service runs it
  eliminated = Set(e.measurement_index for e in diag.eliminations)
  clean = [k in eliminated ? Measurement(typ = m.typ, value = m.value, sigma = m.sigma, active = false, busIdx = m.busIdx, branchIdx = m.branchIdx, direction = m.direction, id = m.id, linkIdx = m.linkIdx) : m for (k, m) in enumerate(meas)]
  _estimate!(net, clean)
  rows = calcLinkFlowsSE!(net)
  @printf("coupler flow after the elimination: %.3f MW (source %s)\n\n", rows[1].pFlow_MW, string(rows[1].source))

  # 2) coupler ring: the balances leave one circulation free; P meters on
  # the three couplers read the balance split plus 5 MW around the ring
  ring = _coupler_ring_net()
  runpf!(ring, 40, 1e-10, 0; method = :rectangular)
  calcNetLosses!(ring)
  calcLinkFlowsKCL!(ring)
  fKcl = [l.pFlow_MW for l in ring.linkVec]
  measR = generateMeasurementsFromPF(ring; noise = false, stddev = _LINK_DEMO_STD)
  for j in 1:3
    addPflowMeasurement!(measR; net = ring, value = fKcl[j] + 5.0, sigma = 1.0, linkNr = j, id = "P_ring_$(j)")
  end
  res = _estimate!(ring, measR)
  rowsR = calcLinkFlowsSE!(ring)
  @printf("ring: %d circulation state(s), J = %.2e at dof %d\n", length(something(res.linkCirculations, Float64[])), res.objectiveJ, res.dof)
  for (j, r) in enumerate(rowsR)
    @printf("  %s: minimum-norm split %8.3f MW, estimated %8.3f MW (source %s)\n", r.name, fKcl[j], r.pFlow_MW, string(r.source))
  end
  return nothing
end

run_example(run_state_estimation_link_measurements)
