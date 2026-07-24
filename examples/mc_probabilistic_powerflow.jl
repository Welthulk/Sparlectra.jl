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

# Date: 2026-07-24
# file: examples/mc_probabilistic_powerflow.jl
# purpose: Monte-Carlo probabilistic power flow on case14 — random load scaling,
#          per-bus voltage statistics from repeated runpf! solves

using Sparlectra
using Random
using Statistics
using Printf

include(joinpath(@__DIR__, "internal", "example_header.jl"))

# --- Study parameters --------------------------------------------------------
const N_SAMPLES = 1000          # Monte-Carlo samples
const LOAD_SIGMA = 0.1          # std of the normal load-scaling factor
const LOAD_FACTOR_MIN = 0.5     # truncation bounds for the scaling factor
const LOAD_FACTOR_MAX = 1.5
const VM_BAND = (0.95, 1.05)    # acceptable voltage band [pu]
const RNG_SEED = 20260724       # reproducibility: same seed -> same numbers
const CASEFILE = "case14.m"

"""Draw a truncated-normal load factor Normal(1, LOAD_SIGMA) in [MIN, MAX] by rejection."""
function drawLoadFactor(rng::AbstractRNG)::Float64
  while true
    f = 1.0 + LOAD_SIGMA * randn(rng)
    LOAD_FACTOR_MIN <= f <= LOAD_FACTOR_MAX && return f
  end
end

function main()
  print_example_banner("examples/mc_probabilistic_powerflow.jl", "Monte-Carlo probabilistic power flow on case14 with random load scaling")

  case_path = ensure_casefile(CASEFILE)
  net = createNetFromMatPowerFile(filename = case_path)
  nb = length(net.nodeVec)

  # Base-case solve.
  ite0, erg0 = runpf!(net, 40, 1e-8, 0)
  erg0 == 0 || error("Base-case power flow did not converge (status=$(erg0)).")
  vm_base = [node._vm_pu for node in net.nodeVec]
  println("Base case solved in $(ite0) iterations.")

  # Load prosumers (non-generators with a specified P or Q). runpf! rebuilds
  # its injection vector from prosumer pVal/qVal on every call
  # (buildComplexSVec), so scaling these fields in place and restoring them
  # afterwards is a full reset — no deepcopy of the network is needed.
  # TODO: Kandidat für scale_load-Port
  load_indices = [k for (k, ps) in enumerate(net.prosumpsVec) if !isGenerator(ps) && (!isnothing(ps.pVal) || !isnothing(ps.qVal))]
  base_p = Dict(k => net.prosumpsVec[k].pVal for k in load_indices)
  base_q = Dict(k => net.prosumpsVec[k].qVal for k in load_indices)
  println("Scaling $(length(load_indices)) load prosumers per sample (P and Q with the same factor).")

  rng = MersenneTwister(RNG_SEED)
  vm_samples = Matrix{Float64}(undef, nb, N_SAMPLES)   # per-bus Vm of converged samples (columns beyond n_converged unused)
  factors = Vector{Float64}(undef, N_SAMPLES)
  n_converged = 0
  n_failed = 0

  for _ in 1:N_SAMPLES
    f = drawLoadFactor(rng)
    for k in load_indices
      ps = net.prosumpsVec[k]
      isnothing(base_p[k]) || (ps.pVal = base_p[k] * f)
      isnothing(base_q[k]) || (ps.qVal = base_q[k] * f)
    end

    _, erg = runpf!(net, 40, 1e-8, 0)
    if erg == 0
      n_converged += 1
      factors[n_converged] = f
      @inbounds for b in 1:nb
        vm_samples[b, n_converged] = net.nodeVec[b]._vm_pu
      end
    else
      n_failed += 1
    end

    # Restore the exact base loads so every sample scales the same base case.
    for k in load_indices
      ps = net.prosumpsVec[k]
      ps.pVal = base_p[k]
      ps.qVal = base_q[k]
    end
  end

  n_converged > 0 || error("No Monte-Carlo sample converged.")
  vm = view(vm_samples, :, 1:n_converged)
  fac = view(factors, 1:n_converged)

  println()
  println("=== Monte-Carlo probabilistic power flow ($(CASEFILE)) ===")
  @printf("samples          : %d\n", N_SAMPLES)
  @printf("converged        : %d (%.1f %%)\n", n_converged, 100.0 * n_converged / N_SAMPLES)
  @printf("failed           : %d\n", n_failed)
  @printf("load factor      : mean=%.4f  min=%.4f  max=%.4f\n", mean(fac), minimum(fac), maximum(fac))
  @printf("seed             : %d\n", RNG_SEED)
  println()
  @printf("%-5s %10s %10s %10s %10s %10s %12s\n", "bus", "Vm_base", "Vm_mean", "Vm_std", "q05", "q95", "out_of_band")
  println("-"^73)
  total_violations = 0
  for b in 1:nb
    row = view(vm, b, :)
    q05, q95 = quantile(row, (0.05, 0.95))
    n_out = count(v -> !(VM_BAND[1] <= v <= VM_BAND[2]), row)
    total_violations += n_out
    @printf("%-5d %10.4f %10.4f %10.5f %10.4f %10.4f %12d\n", b, vm_base[b], mean(row), std(row), q05, q95, n_out)
  end
  println("-"^73)
  n_samples_violating = count(s -> any(b -> !(VM_BAND[1] <= vm[b, s] <= VM_BAND[2]), 1:nb), 1:n_converged)
  @printf("samples with any Vm outside [%.2f, %.2f]: %d of %d (%.1f %%)\n", VM_BAND[1], VM_BAND[2], n_samples_violating, n_converged, 100.0 * n_samples_violating / n_converged)
  @printf("bus-sample violations total               : %d\n", total_violations)

  # Note on parallelization: samples are independent, so this loop could run
  # on Threads.@threads with one deepcopy(net) per thread (never share a Net
  # across threads). Kept serial here for clarity and determinism.
  return nothing
end

main()
