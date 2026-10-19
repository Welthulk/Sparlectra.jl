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

# file: examples/state_estimation/state_estimation_convert_measurements.jl
# purpose: migrate a measurement CSV v1 to the measurement table with
#          convertMeasurementsV1 on the shipped sp_case14: write a v1 set
#          (voltages, bus balances, flows, currents, angles), convert it,
#          show the first table rows and check that both give the same
#          measurements and the same estimate.

using Sparlectra
using Printf

include(joinpath(@__DIR__, "..", "others", "example_header.jl"))

const _ROOT = normpath(joinpath(@__DIR__, "..", ".."))

function run_state_estimation_convert_measurements()
  print_example_banner("examples/state_estimation/state_estimation_convert_measurements.jl", "measurement CSV v1 to measurement table with convertMeasurementsV1 on sp_case14")

  net = importSCF(joinpath(_ROOT, "data", "scf", "sp_case14.scf.json"))
  _, erg = runpf!(net, 40, 1e-10, 0)
  erg == 0 || error("power flow did not converge")

  # a v1 set as older Sparlectra versions wrote it
  dir = mktempdir()
  v1 = joinpath(dir, "sp_case14.v1.csv")
  empty!(net.measurements)
  append!(net.measurements, generateMeasurementsFromPF(net; noise = true, includeVa = true, includeImag = true, rng = Sparlectra.Random.Xoshiro(1)))
  writeMeasurementsCSV(net; file = v1, headerComments = ["case: sp_case14.scf.json"])

  # the conversion: the case gives the keys (element names here; for a
  # CGMES case pass index = cgmesMeasurementKeys(result) for mRIDs)
  table = joinpath(dir, "sp_case14.table.csv")
  res = convertMeasurementsV1(v1, net; out = table)
  print(res.report)
  println("first rows of the table:")
  for line in filter(l -> !startswith(l, "#"), readlines(table))[1:6]
    println("  ", line)
  end

  # the same measurements both ways: the set written to the v1 file (the
  # writer is lossless; the v1 reader was removed in 0.50.0) against the
  # table import
  a = deepcopy(net)
  b = deepcopy(net)
  applyMeasurementSet!(b, only(importMeasurementTable(readMeasurementTable(table), nameMeasurementKeys(net)).sets); zero_injection = false)
  same = all(x.value === y.value && x.sigma === y.sigma && x.typ == y.typ for (x, y) in zip(a.measurements, b.measurements))
  @printf("measurements: %d v1, %d table, values and sigmas bit for bit equal: %s\n", length(a.measurements), length(b.measurements), same)

  ra = with_state_estimation_config(update_net = false) do
    runse!(a)
  end
  rb = with_state_estimation_config(update_net = false) do
    runse!(b)
  end
  @printf("estimate: J = %.6f (v1) and %.6f (table), max |dV| = %.1e pu\n", ra.objectiveJ, rb.objectiveJ, maximum(abs.(ra.voltages .- rb.voltages)))
  return nothing
end

run_example(run_state_estimation_convert_measurements)
