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

# file: docs/lit/warmup/workshop_se_diagnostics.jl
# purpose: compile warm-up of the state estimation diagnostics workshop
#          (docs/lit/workshop_se_diagnostics.jl): runs every path the
#          chapters use once on throwaway nets, so the chapters do not stall
#          on first-call compilation. Included by the first code cell of the
#          notebook, which loads Sparlectra and Random.
#          Not part of the library.

"""
    warmup()

Run every compile-heavy path of the state estimation diagnostics workshop
once and print how long the first solve, the first diagnostics report and
the rest took. The output of the paths is discarded.
"""
function warmup()
  println("warm-up: compiles the paths of this workshop once; how long it takes depends on the machine (a Colab session is several times slower than a desktop)")
  wnet = Net(name = "warmup", baseMVA = 100.0)
  addBus!(net = wnet, busName = "A", vn_kV = 110.0)
  addBus!(net = wnet, busName = "B", vn_kV = 110.0)
  addProsumer!(net = wnet, busName = "A", type = "EXTERNALNETWORKINJECTION", referencePri = "A", vm_pu = 1.0, va_deg = 0.0)
  addProsumer!(net = wnet, busName = "B", type = "ENERGYCONSUMER", p = 10.0, q = 3.0)
  addPIModelACLine!(net = wnet, fromBus = "A", toBus = "B", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
  t_pf = @elapsed runpf!(wnet, 10, 1e-8, 0)
  setMeasurementsFromPF!(wnet; includeVm = true, includePinj = true, includeQinj = true, includePflow = true, includeQflow = true, noise = false)
  t_dg = @elapsed validate_measurements(wnet)
  println("warm: power flow ", round(t_pf; digits = 2), " s, diagnostics ", round(t_dg; digits = 2), " s (first calls compile)")

  ## every further path of the chapters once, output discarded
  t_paths = @elapsed redirect_stdout(devnull) do
    calcNetLosses!(wnet)
    wstd = measurementStdDevs(vm = 1e-3, pinj = 1.0, qinj = 1.0, pflow = 0.7, qflow = 0.7)
    setMeasurementsFromPF!(wnet; includeVm = true, includePinj = true, includeQinj = true, includePflow = true, includeQflow = true, noise = true, stddev = wstd, rng = MersenneTwister(1))
    with_state_estimation_config(max_eliminations = 3) do
      runse_diagnostics(wnet)
    end
    wmeas = Measurement[m for m in wnet.measurements]
    ## tolerance at the finite-difference floor, so the warm-up logs no notice
    with_state_estimation_config(max_iter = 20, tol = 1e-6, update_net = false) do
      runse!(wnet, wmeas)
    end
    with_state_estimation_config(max_iter = 20, tol = 1e-6, update_net = false, robust = true) do
      runse!(wnet, wmeas)
    end
    wp = get_branch_p_from_to_mw(wnet, "A", "B")
    wq = get_branch_q_from_to_mvar(wnet, "A", "B")
    addImagMeasurement!(wnet; fromBus = "A", toBus = "B", direction = :from, value = 1000.0 * sqrt(wp^2 + wq^2) / (sqrt(3.0) * 110.0), sigma = 5.0)
    with_state_estimation_config(flatstart = true, jac_eps = 1e-6) do
      evaluate_global_observability(wnet)
    end
    validate_measurements(wnet)
    ## a meshed throwaway net with a gross flow error, so the report, the
    ## elimination loop and the robust weights take their bad-data branches
    ## (the healthy two-bus set above never leaves the band), plus the result
    ## lines the chapters print from them; the precheck and tolerance notices
    ## of this throwaway net are formatted into devnull (a console logger,
    ## so the notice display compiles too)
    Base.CoreLogging.with_logger(Base.CoreLogging.ConsoleLogger(devnull)) do
      wr = Net(name = "warmup_ring", baseMVA = 100.0)
      addBus!(net = wr, busName = "R1", vn_kV = 110.0, vm_pu = 1.02, va_deg = 0.0)
      for i in 2:4
        addBus!(net = wr, busName = "R$(i)", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
      end
      addPIModelACLine!(net = wr, fromBus = "R1", toBus = "R2", r_pu = 0.010, x_pu = 0.080, b_pu = 0.0, status = 1)
      addPIModelACLine!(net = wr, fromBus = "R2", toBus = "R3", r_pu = 0.011, x_pu = 0.085, b_pu = 0.0, status = 1)
      addPIModelACLine!(net = wr, fromBus = "R3", toBus = "R4", r_pu = 0.012, x_pu = 0.090, b_pu = 0.0, status = 1)
      addPIModelACLine!(net = wr, fromBus = "R4", toBus = "R1", r_pu = 0.010, x_pu = 0.080, b_pu = 0.0, status = 1)
      addPIModelACLine!(net = wr, fromBus = "R1", toBus = "R3", r_pu = 0.009, x_pu = 0.070, b_pu = 0.0, status = 1)
      addProsumer!(net = wr, busName = "R1", type = "EXTERNALNETWORKINJECTION", referencePri = "R1", vm_pu = 1.02, va_deg = 0.0)
      addProsumer!(net = wr, busName = "R2", type = "GENERATOR", p = 60.0, q = 10.0)
      addProsumer!(net = wr, busName = "R3", type = "LOAD", p = 35.0, q = 10.0)
      addProsumer!(net = wr, busName = "R4", type = "LOAD", p = 45.0, q = 15.0)
      validate!(net = wr)
      runpf!(wr, 40, 1e-10, 0)
      calcNetLosses!(wr)
      wvm = [wr.nodeVec[i]._vm_pu for i in eachindex(wr.nodeVec)]
      setMeasurementsFromPF!(wr; includeVm = true, includePinj = true, includeQinj = true, includePflow = true, includeQflow = true, noise = true, stddev = wstd, rng = MersenneTwister(42))
      wbad = findfirst(m -> m.typ == Sparlectra.PflowMeas, wr.measurements)
      m0 = wr.measurements[wbad]
      wr.measurements[wbad] = Measurement(typ = m0.typ, value = m0.value + 25.0, sigma = m0.sigma, busIdx = m0.busIdx, branchIdx = m0.branchIdx, direction = m0.direction, id = m0.id, linkIdx = m0.linkIdx)
      println("corrupted ", m0.id, ": ", round(m0.value; digits = 2), " -> ", round(m0.value + 25.0; digits = 2), " MW (sigma ", m0.sigma, ")")
      wrep = validate_measurements(wr)
      println("J/dof now: ", round(wrep.objective.value / wrep.objective.dof; digits = 1), ", z_wh = ", round(wrep.objective.z_wh; digits = 1), " (reason = ", wrep.objective.reason, ")")
      top = wrep.largest_normalized_residual
      println("prime suspect: ", top.id, "  |r_N| = ", round(top.abs_normalized_residual; digits = 1), ", w_ii = ", round(top.wii; digits = 2), ", localizable = ", top.localizable)
      println("suspicious rows (|r_N| >= 3): ", length(wrep.suspicious_measurements))
      for row in wrep.measurement_ranking[1:2]
        println("  ", rpad(row.id, 22), " |r_N| = ", lpad(round(row.abs_normalized_residual; digits = 1), 6), "   w_ii = ", round(row.wii; digits = 2))
      end
      wdiag = with_state_estimation_config(max_eliminations = 3) do
        runse_diagnostics(wr)
      end
      println("eliminations: ", length(wdiag.eliminations), ", stop reason: ", wdiag.stop_reason)
      for e in wdiag.eliminations
        println("  step ", e.elimination, ": removed ", e.id, " (|r_N| before = ", round(abs(e.normalized_residual_before); digits = 1), ")")
      end
      fin = wdiag.final_diagnostics
      println("after elimination: J/dof = ", round(fin.objective.value / fin.objective.dof; digits = 3), ", z_wh = ", round(fin.objective.z_wh; digits = 3), ", reason = ", fin.objective.reason)
      wbadm = Measurement[m for m in wr.measurements]
      wplain = with_state_estimation_config(max_iter = 20, tol = 1e-8, update_net = false) do
        runse!(wr, wbadm)
      end
      wrob = with_state_estimation_config(max_iter = 20, tol = 1e-8, update_net = false, robust = true) do
        runse!(wr, wbadm)
      end
      werr = maximum(abs.(abs.(wplain.voltages) .- wvm))
      println("max |Vm_est - Vm_true|: plain = ", round(werr; sigdigits = 3), " pu, robust = ", round(maximum(abs.(abs.(wrob.voltages) .- wvm)); sigdigits = 3), " pu")
      for rw in wrob.robustRows
        println("  robust row: ", rw.id, "  stage ", rw.stage, ", t = ", round(rw.t; digits = 1), ", sigma widened by ", round(rw.sigma_factor; digits = 1), "x")
      end
    end
    wsh = Net(name = "warmup_shunt", baseMVA = 100.0)
    addBus!(net = wsh, busName = "A", vn_kV = 110.0)
    addBus!(net = wsh, busName = "B", vn_kV = 110.0)
    addProsumer!(net = wsh, busName = "A", type = "EXTERNALNETWORKINJECTION", referencePri = "A", vm_pu = 1.0, va_deg = 0.0)
    addProsumer!(net = wsh, busName = "B", type = "ENERGYCONSUMER", p = 10.0, q = 3.0)
    addPIModelACLine!(net = wsh, fromBus = "A", toBus = "B", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
    addShunt!(net = wsh, busName = "B", pShunt = 0.0, qShunt = 10.0)
    validate!(net = wsh)
    runpf!(wsh, 40, 1e-10, 0)
    wmsh = generateMeasurementsFromPF(wsh; includeShuntQ = true, noise = false, stddev = measurementStdDevs(vm = 1e-3, pinj = 0.1, qinj = 0.1, pflow = 0.1, qflow = 0.1, shuntq = 0.1))
    setShuntEstimation!(wsh; busName = "B")
    with_state_estimation_config(max_iter = 30, tol = 1e-6, update_net = false, update_shunts = true) do
      runse!(wsh, wmsh)
    end
  end
  println("further paths  : ", round(t_paths; digits = 2), " s (elimination, robust SE, observability, current rows, shunt estimation); everything warm")
  return nothing
end
