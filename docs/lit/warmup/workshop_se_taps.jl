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

# file: docs/lit/warmup/workshop_se_taps.jl
# purpose: compile warm-up of the tap estimation workshop
#          (docs/lit/workshop_se_taps.jl): runs every path the examples use
#          once on throwaway nets, so the examples do not stall on
#          first-call compilation. Included by the first code cell of the
#          notebook, which loads Sparlectra, Random and Printf and defines
#          the helpers used here (build_net, truth_measurements).
#          Not part of the library.

"""
    warmup()

Run every compile-heavy path of the tap estimation workshop once and print
how long the first solve, the first estimation and the rest took. The
output and log notices of the paths are discarded.
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
  t_se = @elapsed runse!(wnet)
  println("warm-up time, first power flow: ", round(t_pf; digits = 2), " s, first estimation: ", round(t_se; digits = 2), " s (both compile)")

  ## every further path of the examples once on the study network, output
  ## and log notices discarded
  t_paths = @elapsed redirect_stdout(devnull) do
    Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do
      wtruth, wmeas = truth_measurements(ratio_step = 1.0, phase_step = 1.0, includeImag = true, includeVa = true, vaBusIdxs = [2, 3], vaRefOffsetDeg = 1.0)
      wv = Measurement[]
      addVmMeasurement!(wv; net = wtruth, busName = "L1", value = 1.0, sigma = 1e-3)
      addPinjMeasurement!(wv; net = wtruth, busName = "L1", value = -25.0, sigma = 1.0)
      addQinjMeasurement!(wv; net = wtruth, busName = "L1", value = -8.0, sigma = 1.0)
      addPflowMeasurement!(wv; net = wtruth, fromBus = "L1", toBus = "L2", value = 1.0, sigma = 0.7, direction = :from)
      addQflowMeasurement!(wv; net = wtruth, fromBus = "L1", toBus = "L2", value = 1.0, sigma = 0.7, direction = :from)
      for mode in (:ratio, :pst, :both)
        wm = build_net()
        setTapEstimation!(wm; trafo = 3, mode = mode)
        wr = with_state_estimation_config(max_iter = 40, tol = 1e-6, update_taps = false) do
          runse!(wm, deepcopy(wmeas))
        end
        filter(r -> r.role === :ratio, wr.tapEstimates)
      end
      ## the radial transformer of Example 6 (frozen, with its notice)
      wb = build_net()
      setTapEstimation!(wb; trafo = 5, mode = :ratio)
      with_state_estimation_config(max_iter = 40, tol = 1e-6) do
        runse!(wb, filter(m -> m.id != "Vm_bus_5", deepcopy(wmeas)))
      end
      wo = build_net()
      append!(wo.measurements, deepcopy(wmeas))
      with_state_estimation_config(flatstart = true, jac_eps = 1e-6) do
        evaluate_global_observability(wo)
      end
      validate_topology(wo, wmeas)
      with_state_estimation_config(robust_mode = :replacement, k_suppress = 4.0, suppression_sigma = 2000.0) do
        runse!(build_net(), deepcopy(wmeas))
      end
      wi = build_net()
      addProsumer!(net = wi, busName = "H2", type = "EXTERNALNETWORKINJECTION", referencePri = "H2", vm_pu = 1.015, va_deg = 0.0)
      for b in (1, 2)
        wi.branchVec[b].status = 0
        wi.branchVec[b].from_status = 0
        wi.branchVec[b].to_status = 0
      end
      validate!(net = wi)
      runpf!(wi, 40, 1e-10, 0; islands_enabled = true)
      calcNetLosses!(wi)
      wmi = generateMeasurementsFromPF(wi; includeVm = true, includePinj = true, includeQinj = true, includePflow = true, includeQflow = true, noise = true, stddev = measurementStdDevs(vm = 1e-3, pinj = 1.0, qinj = 1.0, pflow = 0.7, qflow = 0.7), rng = MersenneTwister(1))
      with_state_estimation_config(max_iter = 40, tol = 1e-6) do
        runse!(wi, wmi)
      end
      @printf("%d %.2f %s\n", 1, 1.0, "warm")
      ## the error against the truth and every result line of the examples:
      ## each @printf format compiles its own method, so they are replayed
      ## here with the argument types of the examples
      wvm = [getNodeVm(n) for n in wtruth.nodeVec]
      werr = maximum(abs.([getNodeVm(n) for n in wtruth.nodeVec] .- wvm))
      wload = wtruth.nodeVec[3]._pƩLoad === nothing ? 0.0 : wtruth.nodeVec[3]._pƩLoad
      j, d, k, b, s, q, f = 1.0, 2, 3, true, "warm", :good, 1.0
      @printf("exactly determined set: %d measurements, %d states\n", k, k)
      @printf("SE quality: J = %.6f, dof = %d\n", j, d)
      @printf("redundant set (%d Qinj rows added): %d measurements, %d states\n", k, k, k)
      @printf("SE quality: J = %.3f, dof = %d, J/dof = %.3f\n", j, d, j)
      @printf("%-20s  %2d  %3d  %7.3f  %6.3f   %.2e\n", s, k, d, j, j, werr)
      @printf("T1 ratio tap not released: model at step %+d, truth at step %+d\n", 0, 4)
      @printf("SE quality: J = %.1f, dof = %d, J inside its 3-sigma band: %s\n", j, d, b)
      @printf("  true value (used to generate measurements):  %+7.0f\n", f)
      @printf("  model before estimation:                     %+7.0f\n", f)
      @printf("  estimated by SE (continuous):                %+7.3f   sigma_s = %.3f steps\n", f, f)
      @printf("  fixed to the nearest step:                   %+7.0f   %.2f sigma_s from the estimate, out of range: %s\n", f, f, b)
      @printf("  model after run (update_taps = false):       %+7.0f   %s\n", f, s)
      @printf("SE quality: J = %.2f, dof = %d\n", j, d)
      @printf("tap estimation switched off, T1 ratio tap at its model step %+.0f: status %s, reason %s\n", f, q, q)
      @printf("SE quality: J = %.1f, dof = %d\n", j, d)
      @printf("T1 phase tap (%s), position [steps, 0 = neutral]:\n", q)
      @printf("  estimated by SE (continuous):                %+7.3f\n", f)
      @printf("  fixed to the nearest step:                   %+7.0f\n", f)
      @printf("  %-7s  %+5.0f  %+13.0f  %+16.3f  %+9.0f\n", q, f, f, f, f)
      @printf("  SE result:                                   %+7.0f   status %s, reason %s\n", f, q, q)
      @printf("  without currents   %12d   %5.3f\n", k, j)
      @printf("  with currents      %12d   %5.3f\n", k, j)
      @printf("largest change of the estimated |V| between the two runs: %.1e pu\n", wload)
      @printf("  without currents       %12d   %-13s   %3d   %.2e\n", k, q, d, werr)
      @printf("  with all currents      %12d   %-13s   %3d   %.2e\n", k, q, d, werr)
      @printf("  true value (added to the PMU angles):  %+7.3f\n", f)
      @printf("  estimated by SE:                       %+7.3f\n", f)
      @printf("SE quality: J = %.2f, dof = %d (%d measurements)\n", j, d, k)
      @printf("tie closed in the model (wrong):   SE quality J = %.1f, dof = %d, J inside its 3-sigma band: %s\n", j, d, b)
      @printf("tie open in the model (corrected): SE quality J = %.1f, dof = %d, J inside its 3-sigma band: %s\n", j, d, b)
      @printf("  %6d  %5d  %-9s  %7.3f  %3d  %s\n", k, k, b, j, d, q)
      @printf("  %6s  %5d  %-9s  %7.3f  %3d\n", s, k, s, j, d)
      @printf("plain WLS, P flow H1-H2 at H1 off by %+.0f MW: SE quality J = %.1f, dof = %d, J inside its 3-sigma band: %s\n", f, j, d, b)
      @printf("  J over all rows, original sigma (honest):   %8.2f, dof = %d\n", j, d)
      @printf("  J_active over the rows the solve trusted:   %8.2f, dof = %d\n", j, d)
      @printf("  rows suppressed:                            %5d\n", k)
    end
  end
  println("warm-up time, further paths: ", round(t_paths; digits = 2), " s (tap estimation, currents, PMU, observability, topology, islands, robust); everything warm")
  return nothing
end
