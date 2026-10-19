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

# Date: 2026-10-19
# file: examples/state_estimation/state_estimation_zero_injection_constraint.jl
# purpose: zero injection at passive buses as weighted rows (state_estimation.zero_injection = on)
#          against equality constraints (constraint) on the shipped sp_case14: condition number of
#          the matrix the solver factorizes, J, dof and the largest zero-injection residual

using Sparlectra
using LinearAlgebra
using SparseArrays
using Printf
using Random

include(joinpath(@__DIR__, "..", "others", "example_header.jl"))

const _ROOT = normpath(joinpath(@__DIR__, "..", ".."))

# Noisy measurement set; the passive buses carry only the ZI pair, no
# injection telemetry (that is where the two modes differ).
function _zi_example_set(net::Net)
  meas = generateMeasurementsFromPF(net; noise = true, relativeSigma = true, rng = MersenneTwister(7))
  passive = Set(findPassiveBuses(net))
  filter!(m -> !(m.typ in (Sparlectra.PinjMeas, Sparlectra.QinjMeas) && m.busIdx in passive), meas)
  addZeroInjectionMeasurements!(meas; net = net, busIdxs = sort(collect(passive)))
  return meas, length(passive)
end

# Condition number of the matrix the Gauss-Newton step factorizes at the
# estimated state: G = H'WH with the ZI rows weighted (on), or the KKT
# matrix [G aC'; aC 0] with G from the measurement rows and C the ZI rows
# (constraint), scaled with a = max diag(G) / max |C| as the solver does.
function _step_condition(net::Net, mode::Symbol)
  jac = with_state_estimation_config(flatstart = false, zero_injection = mode) do
    measurement_jacobian(net)
  end
  H = sparse(jac.H)
  w = Float64[inv(r.sigma^2) for r in jac.rows]
  isZI = [net.measurements[r.index].zero_injection for r in jac.rows]
  if mode === :on
    return cond(Matrix(H' * Diagonal(w) * H))
  end
  M = findall(.!isZI)
  C = H[isZI, :]
  G = H[M, :]' * Diagonal(w[M]) * H[M, :]
  a = maximum(abs, diag(G)) / maximum(abs, C)
  K = [Matrix(G) a*Matrix(C)'; a*Matrix(C) zeros(size(C, 1), size(C, 1))]
  return cond(K)
end

function run_zero_injection_constraint()
  print_example_banner("examples/state_estimation/state_estimation_zero_injection_constraint.jl", "zero injection as weighted rows (on) against equality constraints (constraint) on sp_case14")
  net = importSCF(joinpath(_ROOT, "data", "scf", "sp_case14.scf.json"))
  _, erg = runpf!(net, 60, 1e-8, 0; method = :rectangular)
  erg == 0 || error("power flow did not converge")
  meas, npassive = _zi_example_set(net)
  @printf("sp_case14: %d measurement rows, %d passive buses with a ZI pair\n\n", length(meas), npassive)
  @printf("%-11s %12s %6s %10s %5s %22s\n", "mode", "cond", "it", "J", "dof", "max |ZI residual| (MW)")
  vm = Dict{Symbol,Vector{Float64}}()
  for mode in (:on, :constraint)
    snet = deepcopy(net)
    empty!(snet.measurements)
    append!(snet.measurements, meas)
    res = with_state_estimation_config(zero_injection = mode) do
      runse!(snet)
    end
    res.converged || error("state estimation with zero_injection = $(mode) did not converge")
    diag = with_state_estimation_config(zero_injection = mode) do
      validate_measurements(snet, meas)
    end
    zimax = maximum(abs(row.residual) for row in diag.measurement_ranking if startswith(String(row.id), "ZI"))
    @printf("%-11s %12.3e %6d %10.4f %5d %22.3e\n", mode, _step_condition(snet, mode), res.iterations, res.objectiveJ, res.dof, zimax)
    vm[mode] = Float64[nd._vm_pu for nd in snet.nodeVec]
  end
  @printf("\nmax |dVm| between the two modes: %.2e pu\n", maximum(abs.(vm[:on] .- vm[:constraint])))
  println("constraint: the ZI rows hold to the solver tolerance, leave J and keep their place in the dof.")
  return nothing
end

run_example(run_zero_injection_constraint)
