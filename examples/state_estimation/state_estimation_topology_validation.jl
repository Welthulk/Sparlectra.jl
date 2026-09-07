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

# Date: 2026-08-27
# file: examples/state_estimation/state_estimation_topology_validation.jl
# purpose: topology validation (0.10.0): a wrong service state (a
#          transformer recorded closed that is actually open) poisons the
#          estimation differently from bad telemetry. Shows all three
#          ADVISORY stages: the linear pre-checks, the suspected-station
#          classification when elimination cannot cure the data, and the
#          hypothesis test that re-runs the estimation per status toggle
#          on working copies and only RECOMMENDS.

using Sparlectra
using Printf

include(joinpath(@__DIR__, "..", "others", "example_header.jl"))

function _create_topo_net()
  net = Net(name = "topo_demo", baseMVA = 100.0)
  for (b, vn) in (("H1", 110.0), ("H2", 110.0), ("L1", 20.0), ("L2", 20.0))
    addBus!(net = net, busName = b, vn_kV = vn)
  end
  addProsumer!(net = net, busName = "H1", type = "EXTERNALNETWORKINJECTION", referencePri = "H1", vm_pu = 1.02, va_deg = 0.0)
  addProsumer!(net = net, busName = "L1", type = "ENERGYCONSUMER", p = 25.0, q = 8.0)
  addProsumer!(net = net, busName = "L2", type = "ENERGYCONSUMER", p = 15.0, q = 5.0)
  addPIModelACLine!(net = net, fromBus = "H1", toBus = "H2", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "L1", toBus = "L2", r_pu = 0.02, x_pu = 0.1, b_pu = 0.0, status = 1)
  addPIModelTrafo!(net = net, fromBus = "H1", toBus = "L1", r_pu = 0.002, x_pu = 0.06, b_pu = 0.0, ratio = 1.0, shift_deg = 0.0, status = 1)
  addPIModelTrafo!(net = net, fromBus = "H2", toBus = "L2", r_pu = 0.002, x_pu = 0.06, b_pu = 0.0, ratio = 1.0, shift_deg = 0.0, status = 1)
  return net
end

function run_state_estimation_topology_validation()
  # the TRUTH: transformer H2-L2 is open; the measurements come from that
  # state. The MODEL believes it closed.
  tnet = _create_topo_net()
  setBranchStatus!(tnet.branchVec[4], false)
  ite, erg = runpf!(tnet, 40, 1e-12, 0; method = :rectangular)
  erg == 0 || error("true-state power flow did not converge")
  meas = generateMeasurementsFromPF(tnet; includeImag = true, noise = false)
  net = _create_topo_net()

  # stage 1: the linear pre-checks. The open transformer produced NO rows,
  # so the status contradiction cannot fire here; the checks still say
  # what they looked at (the result is always reported, never silent).
  pre = validate_topology(net, meas)
  println("stage 1: ", pre.summary)

  # stage 2: the fingerprint. Elimination exhausts, the band stays high,
  # and the surviving suspects cluster at the transformer's stations:
  # a wrong service state, not bad telemetry.
  diag = runse_diagnostics(net, meas; max_eliminations = 2, maxIte = 100, tol = 1e-8)
  println("stage 2: stop = ", diag.stop_reason, ", band = ", diag.final_diagnostics.objective.reason)
  for f in something(diag.topology_findings, NamedTuple[])
    println("         suspected at station ", f.location, " (", f.evidence, ")")
  end

  # stage 3: the hypothesis test. Each candidate's status is toggled on a
  # WORKING COPY and the estimation re-run; the true hypothesis lands in
  # the band with J near zero and ranks first. Nothing is switched.
  rep = test_topology_hypotheses(net, meas; maxIte = 100, tol = 1e-8)
  println("stage 3: ", rep.n_candidates, " candidate(s), ambiguous = ", rep.ambiguous)
  for r in rep.recommendations
    @printf("         %-24s %-16s J %10.3f -> %10.3f  %s\n", r.element, r.hypothesis, r.j_before, r.j_after, r.verdict)
  end
  println("the model net is untouched: branch 4 status = ", net.branchVec[4].status, " (still closed; recommendations only)")
  return nothing
end

run_example(run_state_estimation_topology_validation)
