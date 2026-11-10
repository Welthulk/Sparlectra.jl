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

# Date: 2026-11-10
# file: examples/state_estimation/state_estimation_multi_area.jl
# purpose: multi-area state estimation on the shipped sp_case118 (three areas
#          from the MATPOWER area column): the partition, the two-level
#          estimate against the integrated one, and the bad-data elimination
#          with a gross error on a tie branch

using Sparlectra
using Printf
using Random

include(joinpath(@__DIR__, "..", "others", "example_header.jl"))

const MPOWER_DIR_MA = joinpath(@__DIR__, "..", "..", "data", "mpower")

function run_multi_area()
  print_example_banner("examples/state_estimation/state_estimation_multi_area.jl", "multi-area state estimation on sp_case118: partition, two-level estimate against the integrated one, bad data on a tie branch")
  net = Sparlectra.createNetFromMatPowerFile(filename = joinpath(MPOWER_DIR_MA, "sp_case118.m"), flatstart = false, bus_shunt_model = :admittance, matpower_shift_sign = 1.0, matpower_shift_unit = :deg, matpower_ratio = :normal, tap_changer_model = :ideal)
  _, erg = runpf!(net, 60, 1e-8, 0)
  erg == 0 || error("sp_case118: reference power flow did not converge")
  meas = generateMeasurementsFromPF(net; noise = true, rng = MersenneTwister(1))

  # the areas come from the MATPOWER bus table
  part = area_partition(net)
  println("areas ", part.areas, ", ", length(part.tie_branches), " tie branches, border buses ", join([getCompName(net.nodeVec[b].comp) for b in part.border_buses], ", "))

  # two-level estimate; the result carries the integrated estimate on the
  # same rows for the comparison
  r = runse_multi_area!(deepcopy(net), meas)
  @printf("converged %s after %d passes, J %.4f, dof %d\n", r.converged, r.passes, r.objectiveJ, r.dof)
  for a in r.areas
    @printf("  area %d: %d buses, %d border buses, %d rows, J %.2f\n", a.area, length(a.buses), length(a.border_buses), a.rows, a.j)
  end
  @printf("largest difference to the integrated estimate: %.2e pu, %.2e deg\n", r.comparison.max_dvm_pu, r.comparison.max_dva_deg)

  # a gross error on the P flow of a tie branch: found and removed as in
  # the integrated estimate
  k = findfirst(m -> m.typ == Sparlectra.PflowMeas && m.branchIdx == first(part.tie_branches) && m.direction === :from, meas)
  m = meas[k]
  bad = copy(meas)
  bad[k] = Measurement(typ = m.typ, value = m.value + 40 * m.sigma, sigma = m.sigma, branchIdx = m.branchIdx, direction = m.direction, id = m.id)
  d = runse_multi_area_diagnostics(deepcopy(net), bad)
  for t in d.trace
    @printf("eliminated %s (area %d), normalized residual %.2f\n", t.id, t.area, t.normalized_residual_before)
  end
  println("stop reason: ", d.stop_reason)
  return r
end

run_example(run_multi_area)
