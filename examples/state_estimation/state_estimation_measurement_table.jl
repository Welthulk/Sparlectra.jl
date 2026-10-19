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

# file: examples/state_estimation/state_estimation_measurement_table.jl
# purpose: the measurement table end to end. On the shipped sp_case14 (name
#          keys): generate a table from the solved power flow, once without
#          and once with noise, read it back and estimate; then cut the load
#          injection rows and fill them with pseudo rows; sigmas per
#          voltage level (state_estimation.sigma_by_voltage_level) for the
#          generator and for rows without a sigma cell. On the ENTSO-E
#          SmallGrid (mRID keys): read the synthetic fixture, print the
#          import report and let the bad-data elimination find the rows with
#          gross errors.

using Sparlectra
using Printf

include(joinpath(@__DIR__, "..", "others", "example_header.jl"))

const _ROOT = normpath(joinpath(@__DIR__, "..", ".."))

# write rows as a table, read the file back and estimate its first set
function _estimate_table(net::Net, rows, index::MeasurementKeyIndex)
  file = joinpath(mktempdir(), "case.table.csv")
  writeMeasurementTable(file, rows; comments = ["example table"])
  imp = importMeasurementTable(readMeasurementTable(file), index)
  set = first(imp.sets)
  work = deepcopy(net)
  applyMeasurementSet!(work, set)
  se = with_state_estimation_config(update_net = false) do
    runse!(work)
  end
  return set, se
end

function run_state_estimation_measurement_table()
  print_example_banner("examples/state_estimation/state_estimation_measurement_table.jl", "measurement table: generator and estimation on sp_case14, bad data on the SmallGrid fixture")

  # 1) sp_case14: a case without CGMES, keys are element names
  net = importSCF(joinpath(_ROOT, "data", "scf", "sp_case14.scf.json"))
  _, erg = runpf!(net, 40, 1e-10, 0)
  erg == 0 || error("power flow did not converge")
  Vpf = buildVoltageVector(net)
  index = nameMeasurementKeys(net)

  exact = generateMeasurementTableFromPF(net; noise = false, include_injections = true)
  println("first rows of the generated table (name keys):")
  for r in exact[1:4]
    @printf("  %-24s %-12s %s %10.4f  sigma %.4f\n", r.terminal_mrid, r.node_mrid, r.type, r.value, r.sigma)
  end
  set, se = _estimate_table(net, exact, index)
  @printf("noise 0: %d measurements, J = %.2e, max |V_se - V_pf| = %.2e pu\n", length(set.measurements), se.objectiveJ, maximum(abs.(se.voltages .- Vpf)))

  noisy = generateMeasurementTableFromPF(net; coverage = 0.8, seed = 11)
  set, se = _estimate_table(net, noisy, index)
  @printf("noise, coverage 0.8: %d measurements, J = %.2f, dof %d, max |V_se - V_pf| = %.2e pu\n", length(set.measurements), se.objectiveJ, se.dof, maximum(abs.(se.voltages .- Vpf)))

  # pseudo rows: drop the injection rows of the loads, then let the
  # generator fill them from the setpoints of the case with a large sigma
  loads = Set(k for (k, t) in index.terminals if t.kind === :injection && t.detail == "load")
  cut = filter(r -> !(r.terminal_mrid in loads), exact)
  pseudo = generatePseudoMeasurementTable(net, cut; sigma = 10.0)
  set, se = _estimate_table(net, vcat(cut, pseudo), index)
  @printf("%d load injection rows cut, %d pseudo rows (source %s): J = %.2e, max |V_se - V_pf| = %.2e pu\n", count(r -> r.terminal_mrid in loads, exact), length(pseudo), first(pseudo).source, se.objectiveJ, maximum(abs.(se.voltages .- Vpf)))

  # sigmas per voltage level in physical units (placeholders; in a case
  # configuration file this is state_estimation.sigma_by_voltage_level).
  # The generator writes them and draws its noise with them; a row without
  # a sigma cell gets the sigma of its level at the import.
  levels = [Sparlectra.SESigmaLevel(vn_kv = 110, p_mw = 1.0, q_mvar = 1.0, u_kv = 0.3, rel_pct = 1.0), Sparlectra.SESigmaLevel(vn_kv = 20, p_mw = 0.2, q_mvar = 0.2, u_kv = 0.06)]
  noisy = generateMeasurementTableFromPF(net; seed = 11, sigma_by_voltage_level = levels)
  rows = [MeasurementRow(i, r.timestamp, r.measurement_mrid, r.terminal_mrid, r.node_mrid, r.type, r.value, nothing, r.source, r.valid, :none, "") for (i, r) in enumerate(noisy)]
  set = first(importMeasurementTable(rows, index; sigma_by_voltage_level = levels).sets)
  work = deepcopy(net)
  applyMeasurementSet!(work, set)
  se = with_state_estimation_config(update_net = false) do
    runse!(work)
  end
  origin_count(o) = count(==(o), values(set.sigma_origin))
  @printf("sigma by voltage level: %d measurements with the level sigma, %d with the default; J = %.2f, dof %d\n\n", origin_count(:level), origin_count(:default), se.objectiveJ, se.dof)

  # 2) SmallGrid: a CGMES delivery, keys are mRIDs; the fixture carries
  #    four rows with gross errors in its first set
  dir = joinpath(_ROOT, "test", "fixtures", "cgmes_smallgrid")
  table = joinpath(_ROOT, "test", "fixtures", "measurements", "smallgrid_synthetic.table.csv")
  result = redirect_stdout(() -> importCGMES(path = sort!(filter(f -> endswith(f, ".zip"), readdir(dir; join = true)))), devnull)
  imp = importMeasurementTable(readMeasurementTable(table), cgmesMeasurementKeys(result))
  # the summary lines of the import report (the problem rows follow there)
  for line in filter(l -> !startswith(l, "  "), split(measurementImportReport(imp), "\n"))
    isempty(line) || println(line)
  end
  applyMeasurementSet!(result.net, first(imp.sets))
  d = with_state_estimation_config(update_net = false, max_eliminations = 10) do
    runse_diagnostics(result.net, Measurement[m for m in result.net.measurements])
  end
  for e in d.eliminations
    @printf("  eliminated %s, normalized residual %.1f\n", e.id, e.normalized_residual_before)
  end
  o = d.final_diagnostics.objective
  @printf("stop: %s, J = %.1f, dof %d, J/dof = %.3f\n", d.stop_reason, o.value, o.dof, o.value / o.dof)
  return nothing
end

run_example(run_state_estimation_measurement_table)
