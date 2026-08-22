# Copyright 2023–2026 Udo Schmitz
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

# Date: 2026-07-26
# file: examples/run_state_estimation_suite.jl
# purpose: suite runner that executes the state-estimation example programs (WLS, observability, diagnostics, Monte-Carlo study) in fresh subprocesses and reports a summary

include(joinpath(@__DIR__, "internal", "example_suite_runner.jl"))
using Sparlectra
using LinearAlgebra: svdvals

# state_estimation_observability.jl internally includes
# state_estimation_wls.jl; that duplication is contained in its own
# subprocess and log.
const SUITE_SPECS = ExampleSpec[
  ExampleSpec(name = "wls", file = "state_estimation/state_estimation_wls.jl", purpose = "classical WLS state-estimation workflow with synthetic measurements"),
  ExampleSpec(name = "manual_measurements", file = "state_estimation/state_estimation_manual_measurements.jl", purpose = "measurement set built via the manual measurement helpers (addVmMeasurement!, ...)"),
  ExampleSpec(name = "observability", file = "state_estimation/state_estimation_observability.jl", purpose = "progressively deactivates branch-flow measurements and logs observability-redundancy metrics"),
  ExampleSpec(name = "passive_bus_zib_comparison", file = "state_estimation/state_estimation_passive_bus_zib_comparison.jl", purpose = "WLS state estimation with and without zero-injection (ZIB) measurements"),
  ExampleSpec(name = "pmu_angles", file = "state_estimation/state_estimation_pmu_angles.jl", purpose = "PMU voltage-angle measurements with the reference-angle offset state alpha (aligned vs. shifted PMU time base)"),
  ExampleSpec(name = "diagnostics", file = "state_estimation/usage_state_estimation_diagnostics.jl", purpose = "bad-data diagnostics: inject a bad measurement, validate, deactivate and rerun"),
  ExampleSpec(name = "h_matrix_observability", file = "state_estimation/h_matrix_observability_demo.jl", purpose = "small measurement Jacobians H and the public observability helpers"),
  ExampleSpec(name = "mc_study", file = "state_estimation/mc_state_estimation_study.jl", timeout_s = 1200, purpose = "Monte-Carlo WLS state-estimation error study on the 7-bus workshop net"),
]

# Every suite run also documents the STUDY SETUP itself: the measurement
# matrix of the shared 7-bus workshop network with the standard synthetic
# set, plus a stability assessment (rank, redundancy, conditioning).
# Appended to the summary markdown so the suite report answers "what was
# measured where, and how solid is the estimation problem" without
# re-running anything.
function write_measurement_matrix_report(output_dir::AbstractString)
  net = Net(name = "se_suite_matrix", baseMVA = 100.0)
  addBus!(net = net, busName = "B1", vn_kV = 110.0, vm_pu = 1.02, va_deg = 0.0)
  for i in 2:7
    addBus!(net = net, busName = "B$(i)", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
  end
  addPIModelACLine!(net = net, fromBus = "B1", toBus = "B2", r_pu = 0.010, x_pu = 0.080, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "B2", toBus = "B3", r_pu = 0.011, x_pu = 0.085, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "B3", toBus = "B4", r_pu = 0.012, x_pu = 0.090, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "B4", toBus = "B5", r_pu = 0.010, x_pu = 0.080, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "B5", toBus = "B6", r_pu = 0.011, x_pu = 0.085, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "B6", toBus = "B7", r_pu = 0.012, x_pu = 0.090, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "B7", toBus = "B1", r_pu = 0.010, x_pu = 0.080, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "B2", toBus = "B5", r_pu = 0.009, x_pu = 0.070, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "B3", toBus = "B6", r_pu = 0.009, x_pu = 0.070, b_pu = 0.0, status = 1)
  addProsumer!(net = net, busName = "B1", type = "EXTERNALNETWORKINJECTION", referencePri = "B1", vm_pu = 1.02, va_deg = 0.0)
  addProsumer!(net = net, busName = "B3", type = "GENERATOR", p = 60.0, q = 10.0)
  addProsumer!(net = net, busName = "B2", type = "LOAD", p = 35.0, q = 10.0)
  addProsumer!(net = net, busName = "B4", type = "LOAD", p = 45.0, q = 15.0)
  addProsumer!(net = net, busName = "B5", type = "LOAD", p = 25.0, q = 8.0)
  addProsumer!(net = net, busName = "B6", type = "LOAD", p = 30.0, q = 10.0)
  addProsumer!(net = net, busName = "B7", type = "LOAD", p = 20.0, q = 6.0)
  ok, msg = validate!(net = net)
  ok || error("suite matrix net invalid: $msg")
  _, erg = runpf!(net, 40, 1e-10, 0)
  erg == 0 || error("suite matrix reference power flow did not converge")
  setMeasurementsFromPF!(net; includeVm = true, includePinj = true, includeQinj = true, includePflow = true, includeQflow = true, noise = false)

  mj = measurement_jacobian(net)
  gobs = evaluate_global_observability(net; flatstart = true, jacEps = 1e-6)
  sv = svdvals(mj.H)
  cond_H = sv[end] > 0.0 ? sv[1] / sv[end] : Inf
  # heuristic verdict on the WLS problem's conditioning: the normal
  # equations square the condition of H, so cond(H) above ~1e7 means the
  # gain matrix loses more than half the double-precision digits
  verdict = cond_H < 1e4 ? "well conditioned" : cond_H < 1e7 ? "usable, watch the weights" : "ill conditioned (gain matrix loses > 8 digits)"

  path = joinpath(output_dir, "state_estimation_suite_measurement_matrix.md")
  open(path, "w") do io
    println(io, "# Measurement matrix and stability")
    println(io)
    println(io, "Shared 7-bus workshop network, standard synthetic set (noise-free).")
    println(io)
    println(io, "## Stability assessment")
    println(io)
    println(io, "| metric | value |")
    println(io, "|---|---|")
    println(io, "| measurements m | ", gobs.n_measurements, " |")
    println(io, "| states n | ", gobs.n_states, " |")
    println(io, "| numerical rank | ", gobs.numerical_rank, " |")
    println(io, "| structurally / numerically observable | ", gobs.structural_observable, " / ", gobs.numerical_observable, " |")
    println(io, "| redundancy dof (m - n) | ", gobs.dof, " |")
    println(io, "| redundancy ratio m/n | ", round(gobs.redundancy_ratio; digits = 2), " |")
    println(io, "| critical measurements | ", length(gobs.numerical_critical_measurement_indices), " |")
    println(io, "| sigma_max / sigma_min of H | ", round(sv[1]; digits = 4), " / ", round(sv[end]; digits = 6), " |")
    println(io, "| cond(H) | ", round(cond_H; digits = 1), " |")
    println(io, "| quality | ", gobs.quality, " |")
    println(io, "| verdict | ", verdict, " |")
    println(io)
    println(io, "## Occupancy matrix (X = row touches state)")
    println(io)
    println(io, "One row per active measurement; columns are the state vector")
    println(io, "(angles of the non-slack buses, then all magnitudes).")
    println(io)
    println(io, "| # | type | location | ", join(mj.cols, " | "), " |")
    println(io, "|---:|---|---|", repeat("---|", length(mj.cols)))
    tol = 1e-9
    for (i, r) in enumerate(mj.rows)
      marks = [abs(mj.H[i, j]) > tol ? "X" : "" for j in eachindex(mj.cols)]
      println(io, "| ", r.index, " | ", replace(string(r.type), "Meas" => ""), " | ", r.location, " | ", join(marks, " | "), " |")
    end
  end
  return path
end

# SPARLECTRA_EXAMPLE_SUITE_NO_MAIN=1 is a test-only escape hatch: it lets the
# test suite include this file to inspect the registry without running it.
if get(ENV, "SPARLECTRA_EXAMPLE_SUITE_NO_MAIN", "0") != "1"
  run_example(run_example_suite, "state_estimation_suite", "State-estimation example suite", SUITE_SPECS)
  out_dir = normpath(joinpath(@__DIR__, "_out", "state_estimation_suite"))
  if isdir(out_dir)
    matrix_path = Base.invokelatest(write_measurement_matrix_report, out_dir)
    println("  - ", matrix_path)
  end
end
nothing
