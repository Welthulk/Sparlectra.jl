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
#
# file: examples/others/apslf_pv_diagnostic.jl
# purpose: platform diagnostic for the APSLF solver on networks with PV
#          buses. Four checks run in order and print one line each, so a
#          report from another machine names the first failing layer:
#          1. complex sparse LU (SuiteSparse) on a random matrix,
#          2. ldiv! into a column view (the PV recursion writes that way),
#          3. AnalyticLoadFlow's own 9-bus case with PV buses, raw solver,
#          4. Sparlectra: 3-bus ring and sp_case5 (PV buses), NR against
#             APSLF without polish, plus the convergence radius level.
#
# Call: julia --startup-file=no --project=. examples/others/apslf_pv_diagnostic.jl

using Sparlectra
using AnalyticLoadFlow
using LinearAlgebra
using SparseArrays
using Random
using Printf

println("Julia ", VERSION, " on ", Sys.KERNEL, " ", Sys.MACHINE, ", threads ", Threads.nthreads())
println("BLAS: ", BLAS.get_config())
println("AnalyticLoadFlow ", pkgversion(AnalyticLoadFlow), ", Sparlectra ", pkgversion(Sparlectra))

# 1. complex sparse LU
rng = MersenneTwister(1)
n = 400
A = sprandn(rng, ComplexF64, n, n, 0.02) + 10I
b = randn(rng, ComplexF64, n)
x = lu(A) \ b
@printf("1. complex sparse LU residual: %.2e (expected ~1e-14)\n", norm(A * x - b) / norm(b))

# 2. ldiv! into a column view of a dense complex matrix
F = lu(A)
C = zeros(ComplexF64, n, 3)
ldiv!(@view(C[:, 2]), F, b)
@printf("2. ldiv! into a view, residual: %.2e (expected ~1e-14), column 1/3 untouched: %s\n", norm(A * C[:, 2] - b) / norm(b), all(iszero, C[:, 1]) && all(iszero, C[:, 3]))

# 3. AnalyticLoadFlow's own PV case, raw solver
res9 = solve_demo_case(demo_case_9bus(); order = 24, use_pade = true, nr_polish = false, verbose = 0)
res9_mismatch = hasproperty(res9, :residual_inf) ? res9.residual_inf : NaN
@printf("3. AnalyticLoadFlow 9-bus (PV buses), raw solver without polish: converged = %s, mismatch %.2e\n", res9.converged, res9_mismatch)

# 4. Sparlectra: ring3 and sp_case5, NR against APSLF without polish
quiet = OutputConfig(logfile_results = :off, console_summary = false, startup_latency_hint = false)
cfg_nr = SparlectraConfig(powerflow = PowerFlowConfig(solver = :rectangular, rescue = false), output = quiet)
cfg_ap = SparlectraConfig(powerflow = PowerFlowConfig(solver = :apslf), output = quiet)
function ring3()
  net = Net(name = "diag_ring3", baseMVA = 100.0)
  addBus!(net = net, busName = "B1", vn_kV = 110.0, vm_pu = 1.02, va_deg = 0.0)
  addBus!(net = net, busName = "B2", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
  addBus!(net = net, busName = "B3", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
  addPIModelACLine!(net = net, fromBus = "B1", toBus = "B2", r_pu = 0.010, x_pu = 0.080, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "B2", toBus = "B3", r_pu = 0.011, x_pu = 0.085, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "B3", toBus = "B1", r_pu = 0.012, x_pu = 0.090, b_pu = 0.0, status = 1)
  addProsumer!(net = net, busName = "B1", type = "EXTERNALNETWORKINJECTION", referencePri = "B1", vm_pu = 1.02, va_deg = 0.0)
  addProsumer!(net = net, busName = "B2", type = "GENERATOR", p = 20.0, q = 5.0)
  addProsumer!(net = net, busName = "B3", type = "LOAD", p = 30.0, q = 10.0)
  ok, msg = validate!(net = net)
  ok || error("ring3 invalid: $msg")
  return net
end
scf5 = joinpath(pkgdir(Sparlectra), "data", "scf", "sp_case5.scf.json")
for (label, build) in (("ring3", ring3), ("sp_case5", () -> importSCF(scf5)))
  r_nr = run_sparlectra(net = build(), config = cfg_nr)
  r_ap = run_sparlectra(net = build(), config = cfg_ap)
  dvm = maximum(abs.(getfield.(r_nr.net.nodeVec, :_vm_pu) .- getfield.(r_ap.net.nodeVec, :_vm_pu)))
  st = Sparlectra.rectangular_pf_status(r_ap.net)
  @printf("4. %-9s NR %s (%d it, %.1e), APSLF %s (%d pass, %.1e), max |dVm| %.2e pu, %s\n", label, r_nr.outcome, r_nr.iterations, r_nr.final_mismatch, r_ap.outcome, r_ap.iterations, r_ap.final_mismatch, dvm, st.apslf_convergence_line)
end
