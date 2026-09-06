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

# Date: 2026-08-26
# file: examples/state_estimation/bench_takahashi_diagnostics.jl
# purpose: benchmark of the SE diagnostics Omega paths (dense pinv versus the
#          Takahashi selected inverse) over a synthetic-grid size sweep. Run
#          manually (not part of the example suite; minutes at the largest
#          size). Calibration record 2026-08-26 on the reference machine:
#          crossover at about 130 to 150 states -> takahashi_min_states
#          default 200 (conservative, small systems never regress); speedups
#          1.2x at 199 states, 2.0x at 449, 3.2x at 799, 4.6x at 1799
#          (3.68 s -> 0.80 s), agreement <= 5e-13 relative on wii.

using Sparlectra
using Random
using Printf

include(joinpath(@__DIR__, "..", "others", "example_header.jl"))

function _bench_one(nbuses::Int)
  net, _ = build_synthetic_tiled_grid_net(nbuses)
  runpf!(net, 60, 1e-8, 0; method = :rectangular)
  meas = generateMeasurementsFromPF(net; noise = true, rng = MersenneTwister(1))
  res = runse!(net, meas; maxIte = 40, tol = 1e-8, updateNet = false)
  res.converged || return nothing
  V = res.voltages
  nb = length(net.nodeVec)
  slackIdx = Sparlectra._find_slack_idx(net)
  x = Sparlectra._initial_state_vector(net, slackIdx; flatstart = false)
  p = 0
  for i = 1:nb
    if i != slackIdx
      p += 1
      x[p] = angle(V[i])
    end
    x[nb-1+i] = abs(V[i])
  end
  Ybus = createYBUS(net = net)
  am = Sparlectra._active_measurements(meas)
  H, h = Sparlectra._measurement_jacobian_fd(am, net, x, slackIdx, nb, Ybus)
  z = Sparlectra._measurement_vector(am)
  w = Sparlectra._weight_vector(am)
  r = z - h
  # warm both paths, then take the min of two timed runs each
  Sparlectra._residual_diagnostics(H, r, w; minStates = typemax(Int))
  Sparlectra._residual_diagnostics(H, r, w; minStates = 1)
  td = minimum((@elapsed Sparlectra._residual_diagnostics(H, r, w; minStates = typemax(Int))) for _ = 1:2)
  tt = minimum((@elapsed Sparlectra._residual_diagnostics(H, r, w; minStates = 1)) for _ = 1:2)
  d = Sparlectra._residual_diagnostics(H, r, w; minStates = typemax(Int))
  t = Sparlectra._residual_diagnostics(H, r, w; minStates = 1)
  agree = maximum(abs.(d.wii .- t.wii) ./ max.(abs.(d.wii), 1e-12))
  return (states = length(x), m = length(am), dense_s = td, tak_s = tt, agree = agree)
end

function run_bench_takahashi_diagnostics(; sizes = (25, 49, 100, 225, 400, 900))
  print_example_banner("examples/state_estimation/bench_takahashi_diagnostics.jl", "SE diagnostics Omega paths: dense pinv versus Takahashi selected inverse")
  @printf("%8s %8s %8s %12s %14s %9s %10s\n", "buses", "states", "m", "dense [s]", "takahashi [s]", "speedup", "agree")
  for nb in sizes
    b = _bench_one(nb)
    b === nothing && (println(nb, ": SE not converged, skipped"); continue)
    @printf("%8d %8d %8d %12.4f %14.4f %9.2f %10.2e\n", nb, b.states, b.m, b.dense_s, b.tak_s, b.dense_s / b.tak_s, b.agree)
  end
  println("\ncrossover sets state_estimation.takahashi_min_states (default 200).")
  return nothing
end

run_example(run_bench_takahashi_diagnostics)
