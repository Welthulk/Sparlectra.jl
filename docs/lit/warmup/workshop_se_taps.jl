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
  println("warm: power flow ", round(t_pf; digits = 2), " s, estimator ", round(t_se; digits = 2), " s (first calls compile)")

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
        mode === :ratio ? setTapEstimation!(wm; trafo = 3, mode = mode) : setTapEstimation!(wm; trafo = 3, mode = mode, alpha_deg = 90.0)
        with_state_estimation_config(max_iter = 40, tol = 1e-10) do
          runse!(wm, deepcopy(wmeas))
        end
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
      with_state_estimation_config(max_iter = 40, tol = 1e-10) do
        runse!(wi, wmi)
      end
      @printf("%d %.2f %s\n", 1, 1.0, "warm")
      ## the error against the truth and every result line of the examples:
      ## each @printf format compiles its own method, so they are replayed
      ## here with the argument types of the examples
      wvm = [getNodeVm(n) for n in wtruth.nodeVec]
      werr = maximum(abs.([getNodeVm(n) for n in wtruth.nodeVec] .- wvm))
      wload = wtruth.nodeVec[3]._pƩLoad === nothing ? 0.0 : wtruth.nodeVec[3]._pƩLoad
      j, d, k, b, s, q = 1.0, 2, 3, true, "warm", :good
      @printf("exactly determined: m = %d, states = %d, dof = %d, J = %.6f\n", k, k, d, j)
      @printf("redundant:          m = %d, dof = %d, J = %.3f, J/dof = %.3f\n", k, d, j, j)
      @printf("m = %2d  dof = %2d  J = %7.3f  J/dof = %6.3f  max_err = %.2e  added: %s\n", k, d, j, j, werr, s)
      @printf("tap not released: J = %.1f, dof = %d, inside band: %s\n", j, d, b)
      @printf("released: J = %.2f, dof = %d | electrical step = %+.3f -> fixed to %+d (truth %+d), out of range: %s\n", j, d, j, k, 4, b)
      @printf("PST: J = %.2f | electrical step = %+.3f -> fixed to %+d (truth %+d)\n", j, j, k, 3)
      @printf("cascade: J = %.2f | ratio %+.3f -> %+d (truth %+d), phase %+.3f -> %+d (truth %+d)\n", j, j, k, -2, j, k, 2)
      @printf("radial T3: J = %.2f | electrical step = %+.3f -> fixed to %+d (truth %+d)\n", j, j, k, 3)
      @printf("complete set:  m = %d -> %d, J/dof = %.3f -> %.3f, state moves by %.1e pu\n", k, k, j, j, wload)
      @printf("no flows on the tie:      m = %d, observability = %s, dof = %d, max_err = %.2e\n", k, q, d, werr)
      @printf("PMU set: m = %d, J = %.2f, dof = %d, estimated reference offset = %+.3f deg (injected %+.1f)\n", k, j, d, j, 4.0)
      @printf("wrong status:  J = %.1f, dof = %d, inside band: %s\n", j, d, b)
      @printf("total: J = %.2f, dof = %d, islands = %d\n", j, d, k)
      @printf("  island %d: %d buses, converged = %s, J = %.2f, dof = %d, band = %s\n", k, k, b, j, d, q)
      @printf("replacement:  J = %.1f (honest), J_active = %.2f over dof = %d, suppressed rows = %d\n", j, j, d, k)
    end
  end
  println("further paths  : ", round(t_paths; digits = 2), " s (tap estimation, currents, PMU, observability, topology, islands, robust); everything warm")
  return nothing
end
