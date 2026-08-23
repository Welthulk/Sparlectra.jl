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

# file: examples/others/exp_facts_limit_modes.jl
# purpose: FACTS limit characteristics side by side (issue #297 Drafts A/E/F):
#          constant-Q machine box vs STATCOM (Q = V*S_max, linear) vs SVC
#          (Q = V^2*B, quadratic) on one sagging corridor, plus the SSSC
#          injected-voltage window against the TCSC fixed window on a loop
#          network. Theory: docs/src/facts.md.

using Sparlectra

include(joinpath(@__DIR__, "..", "internal", "example_header.jl"))

"""
    main()

Part 1, shunt side: a weak corridor sags well below 1.0 pu. Three
compensators with the SAME 10-MVAr rating at 1.0 pu try to hold the far
bus, each hits its limit, and each limit behaves differently: the machine
box keeps a constant 10 MVAr, the STATCOM delivers V*10 (linear collapse),
the SVC delivers V^2*10 (quadratic collapse). The point: the ranking shows
exactly under the depressed voltage the devices were installed for.

Part 2, series side: on a two-corridor loop the same flow target is given
to a TCSC (fixed reactance window, reaches the target) and to an SSSC with
a tight injectable voltage (window x_base +- V_inj,max/|I| shrinks with the
branch current, pins at the limit with the injected voltage exhausted).

Part 3, the combined device: the same loop with a machine at M2 gets a
UPFC in ONE call (`addUpfcControl!`, issue #325): the SSSC series converter
steers the corridor flow, the STATCOM shunt converter holds the load-bus
voltage, registered as one named device pair. The result rows show the two
actuators of one device with `at_limit` per converter side. The stationary
quadrature model carries no series active-power injection (see
docs/src/facts.md for the limitation).
"""
function main()
  print_example_banner("examples/others/exp_facts_limit_modes.jl", "FACTS limit characteristics: constant-Q vs STATCOM (V*S_max) vs SVC (V^2*B), and the SSSC injected-voltage window vs the TCSC fixed window")

  rating = 10.0   # MVAr at 1.0 pu, shared by all three shunt-side devices

  # weak corridor: Slack -> Mid -> Load, heavy load, visible sag at Mid/Load
  build_corridor = function (; with_machine::Bool)
    cnet = Net(name = "facts_corridor", baseMVA = 100.0)
    for bus in ("Slack", "Mid", "Load")
      addBus!(net = cnet, busName = bus, vn_kV = 110.0)
    end
    addProsumer!(net = cnet, busName = "Slack", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "Slack")
    addProsumer!(net = cnet, busName = "Load", type = "LOAD", p = 60.0, q = 25.0)
    with_machine && addProsumer!(net = cnet, busName = "Mid", type = "GENERATOR", p = 0.0, q = 0.0)
    addPIModelACLine!(net = cnet, fromBus = "Slack", toBus = "Mid", r_pu = 0.02, x_pu = 0.20, b_pu = 0.0, status = 1)
    addPIModelACLine!(net = cnet, fromBus = "Mid", toBus = "Load", r_pu = 0.02, x_pu = 0.20, b_pu = 0.0, status = 1)
    ok, msg = validate!(net = cnet)
    ok || error("corridor net invalid: $msg")
    return cnet
  end

  # constant-Q machine box: at the limit the injection is a fixed 10 MVAr
  box_net = build_corridor(with_machine = true)
  addMachineVoltageControl!(box_net; bus = "Mid", target_bus = "Load", target_vm_pu = 1.0, qmin_mvar = -rating, qmax_mvar = rating)
  run_control!(box_net)
  box = only([c for c in box_net.machineControls if c isa MachineVoltageControl])

  # STATCOM: current-based limit, delivered Q = V * S_max (live bound)
  st_net = build_corridor(with_machine = true)
  addMachineVoltageControl!(st_net; bus = "Mid", target_bus = "Load", target_vm_pu = 1.0, s_max_mva = rating)
  run_control!(st_net)
  st = only([c for c in st_net.machineControls if c isa MachineVoltageControl])
  v_st = get_bus_vm_pu(st_net, "Mid")

  # SVC: clamped susceptance, delivered Q = V^2 * B through the Y-bus
  svc_net = build_corridor(with_machine = false)
  addShuntVoltageControl!(svc_net; bus = "Mid", target_vm_pu = 1.0, bs_min_mvar = -rating, bs_max_mvar = rating)
  run_sparlectra(net = svc_net)
  svc = only(Sparlectra._shunt_controllers(svc_net))
  v_svc = get_bus_vm_pu(svc_net, "Mid")

  # Part 2: series side on a two-corridor loop, target above the natural split
  build_loop = function ()
    lnet = Net(name = "facts_loop", baseMVA = 100.0)
    for bus in ("A", "M1", "M2", "B")
      addBus!(net = lnet, busName = bus, vn_kV = 110.0)
    end
    addProsumer!(net = lnet, busName = "A", type = "EXTERNALNETWORKINJECTION", referencePri = "A", vm_pu = 1.0, va_deg = 0.0)
    # ENERGYCONSUMER with POSITIVE p/q consumes 80 MW at B
    addProsumer!(net = lnet, busName = "B", type = "ENERGYCONSUMER", p = 80.0, q = 20.0)
    addPIModelACLine!(net = lnet, fromBus = "A", toBus = "M1", r_pu = 0.01, x_pu = 0.10, b_pu = 0.0, status = 1)
    addPIModelACLine!(net = lnet, fromBus = "M1", toBus = "B", r_pu = 0.01, x_pu = 0.10, b_pu = 0.0, status = 1)
    addPIModelACLine!(net = lnet, fromBus = "A", toBus = "M2", r_pu = 0.02, x_pu = 0.20, b_pu = 0.0, status = 1)
    addPIModelACLine!(net = lnet, fromBus = "M2", toBus = "B", r_pu = 0.02, x_pu = 0.20, b_pu = 0.0, status = 1)
    ok, msg = validate!(net = lnet)
    ok || error("loop net invalid: $msg")
    return lnet
  end

  tcsc_net = build_loop()
  tcsc = addSeriesReactanceControl!(tcsc_net; fromBus = "A", toBus = "M2", p_target_mw = 35.0, x_min_pu = 0.02, x_max_pu = 0.30)
  run_control!(tcsc_net)

  sssc_net = build_loop()
  sssc = addSeriesReactanceControl!(sssc_net; fromBus = "A", toBus = "M2", p_target_mw = 35.0, v_inj_max_pu = 0.01)
  run_control!(sssc_net)
  printSeriesReactanceControllerSummary(stdout, sssc_net)

  # Part 3: the UPFC composite, one call for the SSSC+STATCOM pair
  upfc_net = build_loop()
  addProsumer!(net = upfc_net, busName = "M2", type = "GENERATOR", p = 0.0, q = 0.0)
  upfc = addUpfcControl!(upfc_net; fromBus = "A", toBus = "M2", shunt_bus = "M2", target_bus = "B", target_vm_pu = 0.99, p_target_mw = 35.0, v_inj_max_pu = 0.08, s_max_mva = 40.0)
  run_control!(upfc_net)
  println()
  println("UPFC composite ", upfc.name, ": one device, two converter rows:")
  for row in controllableElements(upfc_net)
    println("  ", rpad(row.name, 26), rpad(row.device, 48), "actuator ", rpad(string(row.actuator), 18), "at_limit = ", row.at_limit)
  end

  # Part 4: the FULL UPFC (#326), independent P and Q on one line. Needs the
  # shunt at the sending bus, so a small mesh: S(slack)-I=[UPFC]=J-L(load)
  # plus a parallel S-L path.
  fnet = Net(name = "facts_upfc_full", baseMVA = 100.0)
  for b in ("S", "I", "J", "L")
    addBus!(net = fnet, busName = b, vn_kV = 110.0)
  end
  addProsumer!(net = fnet, busName = "S", type = "EXTERNALNETWORKINJECTION", referencePri = "S", vm_pu = 1.0, va_deg = 0.0)
  addProsumer!(net = fnet, busName = "I", type = "GENERATOR", p = 0.0, q = 0.0)
  addProsumer!(net = fnet, busName = "L", type = "ENERGYCONSUMER", p = 90.0, q = 30.0)
  addPIModelACLine!(net = fnet, fromBus = "S", toBus = "I", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = fnet, fromBus = "I", toBus = "J", r_pu = 0.02, x_pu = 0.18, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = fnet, fromBus = "J", toBus = "L", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = fnet, fromBus = "S", toBus = "L", r_pu = 0.02, x_pu = 0.16, b_pu = 0.0, status = 1)
  validate!(net = fnet)
  full = addUpfcControl!(fnet; model = :full, fromBus = "I", toBus = "J", shunt_bus = "I",
                         p_target_mw = 40.0, q_target_mvar = 10.0, q_shunt_mvar = 0.0,
                         v_inj_max_pu = 0.30, s_max_mva = 120.0,
                         deadband_p_mw = 1e-2, deadband_q_mvar = 1e-2, max_outer_iters = 80)
  run_control!(fnet; control_config = ControlConfig(max_outer_iterations = 80))
  fu = full.upfc

  return (
    upfc_full = (p = fu.achieved_p_mw, q = fu.achieved_q_mvar, p_se = fu.p_se_mw, p_sh = fu.p_sh_mw, dc = fu.dc_residual_mw),

    rating = rating,
    box = (q = box.q_mvar, v = get_bus_vm_pu(box_net, "Mid"), at_limit = box.at_limit),
    statcom = (q = st.q_mvar, v = v_st, at_limit = st.at_limit),
    svc = (q = v_svc^2 * svc.bs_mvar, v = v_svc, bs = svc.bs_mvar, at_limit = svc.at_limit),
    tcsc = (p = tcsc.achieved_p_mw, x = tcsc.x_pu, converged = tcsc.converged),
    sssc = (p = sssc.achieved_p_mw, x = sssc.x_pu, x_base = sssc.x_base_pu, i_pu = sssc.i_pu, v_inj = abs(sssc.x_pu - sssc.x_base_pu) * sssc.i_pu, at_limit = sssc.at_limit),
    upfc = (name = upfc.name, p = upfc.series.achieved_p_mw, series_converged = upfc.series.converged, vm = upfc.shunt.achieved_vm_pu, q = upfc.shunt.q_mvar, shunt_at_limit = upfc.shunt.at_limit),
  )
end

result = run_example(main)
r = result.rating
println()
println("Shunt side, all rated ", r, " MVAr at 1.0 pu, all at their capacitive limit:")
println("  machine box : Q = ", round(result.box.q; digits = 2), " MVAr at V = ", round(result.box.v; digits = 4), " pu (constant)")
println("  STATCOM     : Q = ", round(result.statcom.q; digits = 2), " MVAr at V = ", round(result.statcom.v; digits = 4), " pu (= V * S_max, ", round(100 * result.statcom.q / r; digits = 1), " % of rating)")
println("  SVC         : Q = ", round(result.svc.q; digits = 2), " MVAr at V = ", round(result.svc.v; digits = 4), " pu (= V^2 * B,   ", round(100 * result.svc.q / r; digits = 1), " % of rating)")
println("  ranking under sag: constant box > STATCOM (linear) > SVC (quadratic)")
println()
println("Series side, both steering the A->M2 corridor to 35 MW:")
println("  TCSC (window 0.02..0.30 pu): P = ", round(result.tcsc.p; digits = 2), " MW at x = ", round(result.tcsc.x; digits = 4), " pu, converged = ", result.tcsc.converged)
println("  SSSC (V_inj,max 0.01 pu)   : P = ", round(result.sssc.p; digits = 2), " MW at x = ", round(result.sssc.x; digits = 4), " pu (x_base ", result.sssc.x_base, ", |I| = ", round(result.sssc.i_pu; digits = 3), " pu)")
println("    injected voltage ", round(result.sssc.v_inj; digits = 4), " pu of ", 0.01, " pu available, at_limit = ", result.sssc.at_limit)
println()
println("UPFC composite (", result.upfc.name, "), one call, two converters:")
println("  series: P = ", round(result.upfc.p; digits = 2), " MW toward 35 MW, converged = ", result.upfc.series_converged)
println("  shunt : V = ", round(result.upfc.vm; digits = 4), " pu at B with Q = ", round(result.upfc.q; digits = 2), " MVAr, at_limit = ", result.upfc.shunt_at_limit)
println()
println("FULL UPFC (model = :full), independent P and Q on line I->J:")
println("  line   : P = ", round(result.upfc_full.p; digits = 2), " MW (target 40) and Q = ", round(result.upfc_full.q; digits = 2), " MVAr (target 10), both at once")
println("  DC link: P_se = ", round(result.upfc_full.p_se; digits = 3), " MW, P_sh = ", round(result.upfc_full.p_sh; digits = 3), " MW, residual = ", round(result.upfc_full.dc; digits = 4), " MW")
