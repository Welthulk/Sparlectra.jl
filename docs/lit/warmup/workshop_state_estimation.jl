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

# file: docs/lit/warmup/workshop_state_estimation.jl
# purpose: compile warm-up of the state estimation workshop
#          (docs/lit/workshop_state_estimation.jl): runs every path the
#          chapters use once on throwaway nets, so the chapters do not stall
#          on first-call compilation. Included by the first code cell of the
#          notebook, which loads Sparlectra and LinearAlgebra.
#          Not part of the library.

"""
    warmup()

Run every compile-heavy path of the state estimation workshop once and
print how long the first solve, the first estimation and the rest took.
The output of the paths is discarded.
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
  t_se = @elapsed with_state_estimation_config(() -> runse!(wnet); max_iter = 8, tol = 1e-6, flatstart = true, jac_eps = 1e-6, update_net = false)
  println("warm: power flow ", round(t_pf; digits = 2), " s, estimator ", round(t_se; digits = 2), " s (first calls compile)")
  @assert t_pf > 0.0 && t_se > 0.0

  ## every further path of the chapters once, output discarded
  t_paths = @elapsed redirect_stdout(devnull) do
    calcNetLosses!(wnet)
    with_state_estimation_config(flatstart = true, jac_eps = 1e-6) do
      evaluate_global_observability(wnet)
      evaluate_local_observability(wnet, [1, 2])
    end
    measurement_jacobian(wnet)
    wh = [1.0 0.0; 0.0 1.0; 1.0 1.0]
    evaluate_observability_matrix(wh)
    evaluate_local_observability_matrix(wh, [1])
    numerical_row_redundant(wh, 1)
    structural_row_redundant(wh, 1)
    nullspace([1.0 -1.0])
    empty!(wnet.measurements)
    addVmMeasurement!(wnet; busName = "A", value = 1.0, sigma = 0.002)
    addPinjMeasurement!(wnet; busName = "B", value = -10.0, sigma = 1.0)
    addQinjMeasurement!(wnet; busName = "B", value = -3.0, sigma = 1.0)
    addPflowMeasurement!(wnet; fromBus = "A", toBus = "B", value = get_branch_p_from_to_mw(wnet, "A", "B"), sigma = 0.8, direction = :from)
    addQflowMeasurement!(wnet; fromBus = "A", toBus = "B", value = get_branch_q_from_to_mvar(wnet, "A", "B"), sigma = 0.8, direction = :from)
    addQflowMeasurement!(wnet; branchNr = 1, value = get_branch_q_from_to_mvar(wnet, "A", "B"), sigma = 0.8, direction = :from)
    addPmuPhasorMeasurement!(wnet; busName = "B", vm_pu = wnet.nodeVec[2]._vm_pu, va_deg = wnet.nodeVec[2]._va_deg + 2.0, sigmaVm = 0.002, sigmaVa = 0.02)
    with_state_estimation_config(max_iter = 15, tol = 1e-6, flatstart = true, jac_eps = 1e-6, update_net = false) do
      runse!(wnet)
    end
    wp = deepcopy(wnet)
    removeProsumer!(net = wp, busName = "B", type = "ENERGYCONSUMER")
    empty!(wp.measurements)
    runpf!(wp, 40, 1e-10, 0)
    addZeroInjectionMeasurements!(wp; sigma = 1e-6)
    ## the study network of Example 1: buses with a start voltage,
    ## generator and load prosumers, the seeded noisy measurement set, the
    ## estimate written back and read per bus
    ws = Net(name = "warmup_study", baseMVA = 100.0)
    addBus!(net = ws, busName = "S1", vn_kV = 110.0, vm_pu = 1.02, va_deg = 0.0)
    addBus!(net = ws, busName = "S2", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
    addBus!(net = ws, busName = "S3", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
    addPIModelACLine!(net = ws, fromBus = "S1", toBus = "S2", r_pu = 0.010, x_pu = 0.080, b_pu = 0.0, status = 1)
    addPIModelACLine!(net = ws, fromBus = "S2", toBus = "S3", r_pu = 0.011, x_pu = 0.085, b_pu = 0.0, status = 1)
    addPIModelACLine!(net = ws, fromBus = "S3", toBus = "S1", r_pu = 0.012, x_pu = 0.090, b_pu = 0.0, status = 1)
    addProsumer!(net = ws, busName = "S1", type = "EXTERNALNETWORKINJECTION", referencePri = "S1", vm_pu = 1.02, va_deg = 0.0)
    addProsumer!(net = ws, busName = "S2", type = "GENERATOR", p = 60.0, q = 10.0)
    addProsumer!(net = ws, busName = "S3", type = "LOAD", p = 35.0, q = 10.0)
    validate!(net = ws)
    runpf!(ws, 40, 1e-10, 0)
    wstd = measurementStdDevs(vm = 1e-3, pinj = 1.0, qinj = 1.0, pflow = 0.7, qflow = 0.7)
    setMeasurementsFromPF!(ws; includeVm = true, includePinj = true, includeQinj = true, includePflow = true, includeQflow = true, noise = true, stddev = wstd, rng = MersenneTwister(42))
    wse = with_state_estimation_config(max_iter = 12, tol = 1e-6, flatstart = true, jac_eps = 1e-6, update_net = true) do
      runse!(ws)
    end
    println("J / dof = ", round(wse.objectiveJ / wse.dof; digits = 3), ", dof ", wse.dof, ", ", wse.jWithin3Sigma)
    for (name, idx) in sort(collect(ws.busDict); by = last)
      v = wse.voltages[idx]
      println(rpad(name, 4), "  Vm = ", round(abs(v); digits = 4), " pu   Va = ", round(rad2deg(angle(v)); digits = 3), "°")
    end
    ## the matrix examples: 3x3 and 5x3 literals, a single flow row with its
    ## null space
    w3 = [
      1.0 0.0 0.0
      0.0 1.0 0.0
      0.0 0.0 1.0
    ]
    wo3 = evaluate_observability_matrix(w3)
    println(wo3.numerical_observable, wo3.dof, wo3.numerical_critical_measurement_indices, wo3.unobservable_state_columns)
    wflow = [1.0 -1.0]
    println(vec(round.(nullspace(wflow); digits = 4)))
    isapprox(1.0, 1.0; atol = 1e-9)
  end
  println("further paths  : ", round(t_paths; digits = 2), " s (observability, manual and zero-injection rows, Jacobian, PMU); everything warm")
  return nothing
end
