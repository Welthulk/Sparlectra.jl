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

# file: docs/lit/warmup/workshop_tour.jl
# purpose: compile warm-up of the workshop tour (docs/lit/workshop_tour.jl):
#          runs every path the chapters use once on throwaway nets, so the
#          chapters do not stall on first-call compilation. Included by the
#          first code cell of the notebook, which defines the helpers used
#          here (solve!, build_ring7, build_grid, build_oltc, build_qu).
#          Not part of the library.

"""
    warmup()

Run every compile-heavy path of the workshop tour once and print how long
the first solve and the rest took. The output of the paths is discarded.
"""
function warmup()
  println("warm-up: compiles the paths of this workshop once; how long it takes depends on the machine (a Colab session is several times slower than a desktop)")
  ## tiny warm-up net: a grid connection WITH declared short-circuit data
  wnet = Net(name = "warmup", baseMVA = 100.0)
  addBus!(net = wnet, busName = "A", vn_kV = 110.0)
  addBus!(net = wnet, busName = "B", vn_kV = 110.0)
  addExternalGrid!(net = wnet, busName = "A", vm_pu = 1.0, sk_max_MVA = 2000.0, sk_min_MVA = 1500.0, rx_max = 0.1, internal_impedance = false)
  addProsumer!(net = wnet, busName = "B", type = "ENERGYCONSUMER", p = 10.0, q = 3.0)
  addPIModelACLine!(net = wnet, fromBus = "A", toBus = "B", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)

  t_first = @elapsed runpf!(wnet, 10, 1e-8, 0)
  t_second = @elapsed runpf!(wnet, 10, 1e-8, 0)
  println("power flow     : first solve ", round(t_first; digits = 2), " s (compiles), second ", round(t_second * 1000; digits = 2), " ms")

  ## every further path of the tour once, output discarded: the result
  ## tables, the model edits with the MATPOWER export, bus links, the short
  ## circuit, the tap control loop and Q(U)
  t_paths = @elapsed redirect_stdout(devnull) do
    wr = build_ring7("warmup_ring")
    etime, ite = solve!(wr)
    printACPFlowResults(wr, etime, ite, 1e-8)
    br = getNetBranchNumberVec(net = wr, fromBus = "B3", toBus = "B6")
    setNetBranchStatus!(net = wr, branchNr = br[1], status = 0)
    removeACLine!(net = wr, fromBus = "B2", toBus = "B5")
    markIsolatedBuses!(net = wr, log = false)
    validate!(net = wr)
    solve!(wr)
    writeMatpowerCasefile(wr, joinpath(mktempdir(), "warmup.m"))
    wl = Net(name = "warmup_links", baseMVA = 100.0)
    for b in ("S", "L1", "L2")
      addBus!(net = wl, busName = b, vn_kV = 110.0)
    end
    addProsumer!(net = wl, busName = "S", type = "EXTERNALNETWORKINJECTION", referencePri = "S", vm_pu = 1.0, va_deg = 0.0)
    addProsumer!(net = wl, busName = "L2", type = "ENERGYCONSUMER", p = 10.0, q = 2.0)
    addPIModelACLine!(net = wl, fromBus = "S", toBus = "L1", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
    addLink!(net = wl, fromBus = "L1", toBus = "L2", status = 1)
    validate!(net = wl)
    solve!(wl)
    calcLinkFlowsKCL!(wl)
    wg = build_grid(:slack)
    etime, ite = solve!(wg)
    printACPFlowResults(wg, etime, ite, 1e-8)
    printShortCircuitResult(runShortCircuit!(wg; case = :max))
    wo = build_oltc()
    addTapController!(wo; trafo = "T1", mode = :voltage, target_bus = "Load", target_vm_pu = 1.0, control_ratio = true, control_phase = false, is_discrete = true, deadband_vm_pu = 5e-3)
    run_sparlectra(net = wo)
    printTapControllerSummary(stdout, wo)
    solve!(build_qu(5.0, 1.0))
    # a machine that hits its Q limit (PV to PQ switching, the Q-limit log)
    wq = Net(name = "warmup_qlimits", baseMVA = 100.0)
    for b in ("Q1", "Q2", "Q3")
      addBus!(net = wq, busName = b, vn_kV = 110.0)
    end
    addProsumer!(net = wq, busName = "Q1", type = "EXTERNALNETWORKINJECTION", referencePri = "Q1", vm_pu = 1.0, va_deg = 0.0)
    addProsumer!(net = wq, busName = "Q2", type = "GENERATOR", p = 20.0, vm_pu = 1.05, qMin = -5.0, qMax = 5.0)
    addProsumer!(net = wq, busName = "Q3", type = "ENERGYCONSUMER", p = 45.0, q = 20.0)
    addPIModelACLine!(net = wq, fromBus = "Q1", toBus = "Q2", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
    addPIModelACLine!(net = wq, fromBus = "Q2", toBus = "Q3", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
    validate!(net = wq)
    solve!(wq)
    printQLimitLog(wq)
    distributeBusResults!(wq)
    # the distributed slack of chapter 3
    solve!(build_grid(:slack); distributed_slack_enabled = true, distributed_slack_p_mode = :pg_weighted)
  end
  println("further paths  : ", round(t_paths; digits = 2), " s (tables, edits, export, links, Q limits, distributed slack, short circuit, taps, Q(U)); everything warm")
  return nothing
end
