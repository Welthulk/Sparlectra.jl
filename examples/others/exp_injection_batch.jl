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
# file: examples/others/exp_injection_batch.jl
# purpose: the injection batch (runInjectionBatch) on the shipped sp_case188:
#          32 random load scenarios (factor per bus in 0.8 to 1.2, seed 1)
#          and a ramp of 96 steps from 0.8 to 1.2 of the load, Q limits
#          off, with the share of scenarios in which a PV bus would leave
#          its reactive limits; the ramp once more from the solved base
#          case (start = :base). Started in an environment that loads CUDA
#          and CUDSS, the same batches also run on the GPU and are compared,
#          and the random batch runs in the two-phase variants g5 and g6.
#          julia -t 4 --project=. examples/others/exp_injection_batch.jl

using Sparlectra
using Printf

include(joinpath(@__DIR__, "..", "others", "example_header.jl"))

const SCF_CASE = joinpath(dirname(dirname(pathof(Sparlectra))), "data", "scf", "sp_case188.scf.json")

function summary(label, r, seconds)
  conv = count(r.converged)
  viol = count(s -> r.converged[s] && r.qlimit_violations[s] > 0, eachindex(r.converged))
  @printf("%-22s %-4s %3d of %3d converged, iterations %d to %d, Q limit violated in %3d, Vmin %.4f pu, %.3f s\n",
    label, String(r.backend), conv, length(r.converged), minimum(r.iterations), maximum(r.iterations), viol, minimum(r.vm[:, r.converged]), seconds)
  return nothing
end

function main()
  print_example_banner("examples/others/exp_injection_batch.jl", "many load scenarios of one network: random factors and a 96-step ramp, on CPU threads (and the GPU if loaded)")
  net = importSCF(SCF_CASE)
  println("network: ", length(net.nodeVec), " buses, Julia threads: ", Threads.nthreads())

  random = random_load_scenarios(net, 32; band = 0.2, seed = 1)
  ramp = load_ramp_scenarios(net, 96; low = 0.8, high = 1.2)
  batches = [("random, 32 scenarios", random), ("ramp, 96 steps", ramp)]

  cpu = Dict{String,Any}()
  for (label, sc) in batches
    t = @elapsed r = runInjectionBatch(net, sc)
    cpu[label] = r
    summary(label, r, t)
  end

  # the same ramp started from the solved base case instead of the net's
  # start values, as the QSTS steps do
  t = @elapsed warm = runInjectionBatch(net, ramp; start = :base)
  summary("ramp, from the base", warm, t)

  # the ramp step by step: where the first Q-limit violation appears
  r = cpu["ramp, 96 steps"]
  first_viol = findfirst(s -> r.converged[s] && r.qlimit_violations[s] > 0, eachindex(r.converged))
  println(first_viol === nothing ? "no PV bus leaves its Q limits on the ramp" :
          @sprintf("first Q-limit violation on the ramp at step %d (load factor %.3f), %d PV bus(es)", first_viol, ramp.load_scale[1, first_viol], r.qlimit_violations[first_viol]))

  if Sparlectra.gpu_status().status === :available
    println()
    for (label, sc) in batches
      t = @elapsed g = runInjectionBatch(net, sc; backend = :gpu)
      summary(label, g, t)
      ref = cpu[label]
      both = g.converged .& ref.converged
      @printf("  against the CPU batch: max |dVm| %.1e pu, max |dVa| %.1e deg, same convergence %s, same iterations %s%s\n",
        maximum(abs.(g.vm[:, both] .- ref.vm[:, both])), maximum(abs.(g.va_deg[:, both] .- ref.va_deg[:, both])), g.converged == ref.converged, g.iterations == ref.iterations,
        isempty(g.note) ? "" : string(" (", g.note, ")"))
    end
    # two phases: Float32 on the device first, then Float64 steps with the
    # last Float32 LU (g5) or a Float64 Newton on the CPU (g6)
    for v in (:g5, :g6)
      t = @elapsed g = runInjectionBatch(net, random; backend = :gpu, gpu_two_phase = v)
      summary(string("random, two-phase ", v), g, t)
      # a network the GPU batch does not cover ran on the CPU (see `note`)
      g.backend === :gpu && @printf("  phase-1 steps %d to %d, phase-2 steps %d to %d, switch: %s\n", extrema(g.stats[:tp_phase1_steps])..., extrema(g.stats[:tp_phase2_steps])...,
        join(sort(unique(string.(g.stats[:tp_switch_reason]))), ", "))
    end
  else
    println("GPU: ", Sparlectra.gpu_status().reason, ", the GPU batch is not shown")
  end
  return nothing
end

run_example(main)
