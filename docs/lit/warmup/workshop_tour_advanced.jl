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

# file: docs/lit/warmup/workshop_tour_advanced.jl
# purpose: compile warm-up of the advanced workshop tour
#          (docs/lit/workshop_tour_advanced.jl): runs every path the chapters
#          use once on throwaway nets, so the chapters do not stall on
#          first-call compilation. Included by the first code cell of the
#          notebook, which defines the helpers used here (solve!, bus_va_deg
#          and the network builders). Not part of the library.

"""
    warmup()

Run every compile-heavy path of the advanced workshop tour once and print
how long the first solve and the single paths took. The output of the
further paths is discarded.
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

  t_sc = @elapsed runShortCircuit!(wnet; case = :max)
  println("short circuit  : ", round(t_sc; digits = 2), " s")

  ## HVDC pair path: two 2-bus islands coupled by a controller
  whv = Net(name = "warmup_hvdc", baseMVA = 100.0)
  for b in ("W1", "W2", "W3", "W4")
    addBus!(net = whv, busName = b, vn_kV = 110.0)
  end
  addProsumer!(net = whv, busName = "W1", type = "EXTERNALNETWORKINJECTION", referencePri = "W1", vm_pu = 1.0, va_deg = 0.0)
  addProsumer!(net = whv, busName = "W3", type = "EXTERNALNETWORKINJECTION", referencePri = "W3", vm_pu = 1.0, va_deg = 0.0)
  addProsumer!(net = whv, busName = "W2", type = "ENERGYCONSUMER", p = 5.0, q = 1.0)
  addProsumer!(net = whv, busName = "W4", type = "ENERGYCONSUMER", p = 5.0, q = 1.0)
  addPIModelACLine!(net = whv, fromBus = "W1", toBus = "W2", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = whv, fromBus = "W3", toBus = "W4", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
  addProsumer!(net = whv, busName = "W2", type = "GENERATOR", p = -5.0, q = 0.0)
  addProsumer!(net = whv, busName = "W4", type = "GENERATOR", p = 5.0, q = 0.0)
  addHvdcPairControl!(whv; from_bus = "W2", to_bus = "W4", p_transfer_mw = 5.0)
  t_hvdc = @elapsed run_control!(whv; controllers = collect_outer_controllers(whv), pf_config = PowerFlowConfig(method = :rectangular, max_iter = 15, tol = 1e-8), control_config = ControlConfig(max_outer_iterations = 4, trace = false))
  println("HVDC control   : ", round(t_hvdc; digits = 2), " s")

  ## state-estimation path: synthetic measurements plus one WLS run
  setMeasurementsFromPF!(wnet; includeVm = true, includePinj = true, includeQinj = true, includePflow = true, includeQflow = true, noise = false)
  t_se = @elapsed with_state_estimation_config(() -> runse!(wnet); max_iter = 8, tol = 1e-6, flatstart = true, jac_eps = 1e-6, update_net = false)
  println("state estimator: ", round(t_se; digits = 2), " s")

  ## chapter-5 path: one single-case contingency batch on the warm-up net
  t_n1 = @elapsed runContingencies!(wnet, generateN1Branches(wnet))
  println("contingency    : ", round(t_n1; digits = 2), " s")

  ## every further path of the tour once, output discarded: remote voltage
  ## control, the HVDC variants with their result tables, the noisy SE with
  ## observability, the FACTS devices, the N-1 report and export, and the
  ## serial/parallel short-circuit sweep
  t_paths = @elapsed redirect_stdout(devnull) do
    wpf = PowerFlowConfig(max_iter = 30, tol = 1e-9)
    wrc = PowerFlowConfig(method = :rectangular, max_iter = 25, tol = 1e-8)
    ## chapter 1: remote voltage control
    wr = build_rvc(-50.0, 50.0)
    runpf!(wr; config = wpf, verbose = 0)
    addMachineVoltageControl!(wr; bus = "GenBus", target_bus = "Load", target_vm_pu = 1.05, deadband_vm_pu = 5e-4)
    run_control!(wr; controllers = collect_outer_controllers(wr), pf_config = wpf, control_config = ControlConfig(max_outer_iterations = 15), verbose = 0)
    printMachineControllerSummary(stdout, wr)
    ## chapter 2: HVDC pair, grid-forming, droop, meshed
    wb = build_b2b("warmup_b2b")
    addHvdcLink!(wb; from_bus = "A2", to_bus = "C2")
    etime, ite = solve!(wb; islands_enabled = true)
    printACPFlowResults(wb, etime, ite, 1e-8)
    addHvdcPairControl!(wb; from_bus = "A2", to_bus = "C2", p_transfer_mw = 120.0, loss_mw = 4.0, p_rating_mw = 150.0)
    wres = run_control!(wb; controllers = collect_outer_controllers(wb), pf_config = wrc, control_config = ControlConfig(max_outer_iterations = 8, trace = false))
    calcNetLosses!(wb)
    printACPFlowResults(wb, etime, wres.last_pf_iterations, 1e-8)
    bus_va_deg(wb, "A2")
    ws = build_b2b_source("warmup_b2b_source")
    solve!(ws; islands_enabled = true)
    p_c = get_branch_p_from_to_mw(ws, "C2", "C1")
    wm = build_b2b_source("warmup_b2b_mirrored"; sending_mw = -(p_c + 4.0))
    addHvdcLink!(wm; from_bus = "A2", to_bus = "C2")
    etime, ite = solve!(wm; islands_enabled = true)
    printACPFlowResults(wm, etime, ite, 1e-8)
    wf = build_b2b_source("warmup_b2b_grid_forming")
    addHvdcPairControl!(wf; from_bus = "A2", to_bus = "C2", mode = :island_feed, loss_mw = 4.0, p_rating_mw = 150.0)
    wres = run_control!(wf; controllers = collect_outer_controllers(wf), pf_config = wrc, control_config = ControlConfig(max_outer_iterations = 8, trace = false))
    calcNetLosses!(wf)
    printACPFlowResults(wf, etime, wres.last_pf_iterations, 1e-8)
    wd = build_b2b_droop("warmup_b2b_droop"; sending_mw = -(p_c + 4.0))
    etime, ite = solve!(wd; islands_enabled = true)
    printACPFlowResults(wd, etime, ite, 1e-8)
    wt = build_meshed("warmup_meshed"; c1_model = :pv)
    addHvdcPairControl!(wt; from_bus = "A2", to_bus = "C2", p_transfer_mw = 120.0, loss_mw = 4.0, p_rating_mw = 150.0)
    wres = run_control!(wt; controllers = collect_outer_controllers(wt), pf_config = wrc, control_config = ControlConfig(max_outer_iterations = 8, trace = false))
    calcNetLosses!(wt)
    printACPFlowResults(wt, etime, wres.last_pf_iterations, 1e-8)
    ## chapter 3: noisy measurements, observability, SE with net update
    wse = build_ring7("warmup_se")
    runpf!(wse, 40, 1e-10, 0)
    setMeasurementsFromPF!(wse; includeVm = true, includePinj = true, includeQinj = true, includePflow = true, includeQflow = true, noise = true, stddev = measurementStdDevs(vm = 1e-3, pinj = 1.0, qinj = 1.0, pflow = 0.7, qflow = 0.7), rng = MersenneTwister(1))
    ## the chapter calls these at top level in do-block form: the dynamic
    ## entry of that call compiles on its own, so call it through a barrier
    wsecfg = Base.inferencebarrier(with_state_estimation_config)
    wsecfg(flatstart = true, jac_eps = 1e-6) do
      evaluate_global_observability(wse)
    end
    wsest = wsecfg(max_iter = 12, tol = 1e-6, flatstart = true, jac_eps = 1e-6, update_net = true) do
      runse!(wse)
    end
    ## the per-bus table of the estimate, sorted by bus index
    for (name, idx) in sort(collect(wse.busDict); by = last)
      println(name, round(abs(wsest.voltages[idx]); digits = 4), round(rad2deg(angle(wsest.voltages[idx])); digits = 3))
    end
    ## chapter 4: machine box, STATCOM, SVC, switched bank, TCSC, SSSC, UPFC
    wc = build_sag_corridor("warmup_box"; with_machine = true)
    addMachineVoltageControl!(wc; bus = "Mid", target_bus = "Load", target_vm_pu = 1.0, qmin_mvar = -10.0, qmax_mvar = 10.0)
    run_control!(wc)
    wc = build_sag_corridor("warmup_statcom"; with_machine = true)
    addMachineVoltageControl!(wc; bus = "Mid", target_bus = "Load", target_vm_pu = 1.0, s_max_mva = 10.0)
    run_control!(wc)
    controllableElements(wc)
    wc = build_sag_corridor("warmup_msc"; with_machine = false)
    addShuntVoltageControl!(wc; bus = "Mid", target_vm_pu = 0.95, bs_min_mvar = -40.0, bs_max_mvar = 40.0, step_mvar = 10.0)
    run_control!(wc)
    wc = build_sag_corridor("warmup_svc"; with_machine = false)
    addShuntVoltageControl!(wc; bus = "Mid", target_vm_pu = 1.0, bs_min_mvar = -10.0, bs_max_mvar = 10.0)
    run_control!(wc)
    wl = build_facts_loop("warmup_tcsc")
    addSeriesReactanceControl!(wl; fromBus = "A", toBus = "M2", p_target_mw = 35.0, x_min_pu = 0.02, x_max_pu = 0.30)
    run_control!(wl)
    wl = build_facts_loop("warmup_sssc")
    addSeriesReactanceControl!(wl; fromBus = "A", toBus = "M2", p_target_mw = 35.0, v_inj_max_pu = 0.01)
    run_control!(wl)
    wl = build_facts_loop("warmup_upfc")
    addProsumer!(net = wl, busName = "M2", type = "GENERATOR", p = 0.0, q = 0.0)
    addUpfcControl!(wl; fromBus = "A", toBus = "M2", shunt_bus = "M2", target_bus = "B", target_vm_pu = 0.99, p_target_mw = 35.0, v_inj_max_pu = 0.08, s_max_mva = 40.0)
    run_control!(wl)
    controllableElements(wl)
    wu = build_facts_mesh("warmup_upfc_full")
    addUpfcControl!(wu; model = :full, fromBus = "I", toBus = "J", shunt_bus = "I", p_target_mw = 40.0, q_target_mvar = 10.0, q_shunt_mvar = 0.0, v_inj_max_pu = 0.30, s_max_mva = 120.0, deadband_p_mw = 1e-2, deadband_q_mvar = 1e-2, max_outer_iters = 80)
    run_control!(wu; control_config = ControlConfig(max_outer_iterations = 80))
    calcNetLosses!(wu)
    printACPFlowResults(wu, 0.0, 1, 1e-8)
    ## chapter 5: N-1 with a voltage band and a spur pair that an outage
    ## strands without a reference, report, ranking and CSV export
    wn = build_ring7("warmup_n1")
    addBus!(net = wn, busName = "B8", vn_kV = 110.0)
    addBus!(net = wn, busName = "B9", vn_kV = 110.0)
    addPIModelACLine!(net = wn, fromBus = "B4", toBus = "B8", r_pu = 0.02, x_pu = 0.10, b_pu = 0.0, status = 1)
    addPIModelACLine!(net = wn, fromBus = "B8", toBus = "B9", r_pu = 0.02, x_pu = 0.10, b_pu = 0.0, status = 1)
    addProsumer!(net = wn, busName = "B8", type = "LOAD", p = 6.0, q = 2.0)
    addProsumer!(net = wn, busName = "B9", type = "LOAD", p = 4.0, q = 1.0)
    validate!(net = wn)
    wres = runContingencies!(wn, generateN1Branches(wn); vm_min_pu = 0.95, vm_max_pu = 1.05)
    printContingencyResults(wres)
    wranked = sort([r for r in wres if r.converged]; by = r -> r.min_vm_pu)
    println(wranked[1].name, length(wranked[1].voltage_violations))
    writeContingencyResultsCSV(joinpath(mktempdir(), "warmup_n1.csv"), wres)
    ## chapter 6: the short-circuit sweep, serial and parallel, on two
    ## rings large enough (60 buses each) for the selected-inverse pass
    wp = Net(name = "warmup_parallel", baseMVA = 100.0)
    for k in 1:2
      for i in 1:60
        addBus!(net = wp, busName = "P$(k)_B$(i)", vn_kV = 110.0)
      end
      addExternalGrid!(net = wp, busName = "P$(k)_B1", vm_pu = 1.0, sk_max_MVA = 2000.0 + 100.0 * k, sk_min_MVA = 1500.0, rx_max = 0.1, internal_impedance = false)
      for i in 1:60
        addPIModelACLine!(net = wp, fromBus = "P$(k)_B$(i)", toBus = "P$(k)_B$(i == 60 ? 1 : i + 1)", r_pu = 0.001, x_pu = 0.004, b_pu = 0.0, status = 1)
      end
    end
    ## two rings, two slacks: validate! logs an @info about it, not wanted here
    Base.CoreLogging.with_logger(() -> validate!(net = wp), Base.CoreLogging.NullLogger())
    wsc_ser = runShortCircuit!(wp; case = :max, parallel_enabled = false)
    wsc_par = runShortCircuit!(wp; case = :max, parallel_min_work_items = 2)
    println(isequal(wsc_ser.rows, wsc_par.rows))
  end
  println("further paths  : ", round(t_paths; digits = 2), " s (control loops, HVDC, SE, FACTS, N-1 report, result tables); everything warm")
  return nothing
end
