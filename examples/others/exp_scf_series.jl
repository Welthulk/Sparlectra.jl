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
# Date: 2026-10-14
# file: examples/others/exp_scf_series.jl
# purpose: the shipped SCF time series data/qsts/sp_case14_day.qsts.zip (96
#          steps of 15 minutes on sp_case14): read it, solve one step by
#          hand with apply! and runpf!, then every step with the snapshot
#          run (QSTS mode A), once with the taps of the case and once with
#          its tap controller, print the day at full hours and write the
#          result files; then a few steps with the controllers of
#          data/qsts/sp_case14_day.control.yaml; last, the whole day in order
#          with those controllers (QSTS mode B): taps and bank carried from
#          step to step, the switching list and the metrics.
#          julia --project=. examples/others/exp_scf_series.jl

using Sparlectra
using Printf

include(joinpath(@__DIR__, "..", "others", "example_header.jl"))

const SERIES_ZIP = joinpath(dirname(dirname(pathof(Sparlectra))), "data", "qsts", "sp_case14_day.qsts.zip")
const SERIES_CONTROL = joinpath(dirname(SERIES_ZIP), "sp_case14_day.control.yaml")

# the largest value of a quantity column of one step, ignoring the gaps
_step_max(values::AbstractVector{Float64}) = (v = filter(!isnan, values); isempty(v) ? NaN : maximum(v))

function main()
  print_example_banner("examples/others/exp_scf_series.jl", "an SCF time series: one base case, one delta per step, one zip, every step solved on its own")

  # 1. the series: base case, steps, ops per step
  series = read_scf_series(SERIES_ZIP)
  println(series.name, ": ", length(series.steps), " steps of ", series.step_minutes, " min, ", sum(length, series.ops), " patch operations, base case ",
    series.base.sparlectra.meta["case_name"])
  println("step 57 (", series.steps[57].label, "): ", length(series.ops[57]), " ops, among them ", only(filter(op -> op.op === :status, series.ops[57])))
  println()

  # 2. one step by hand: apply the delta to a copy of the base net and solve
  base = build_net(series)
  net = deepcopy(base)
  apply!(net, series.ops[57], series.index)
  iterations, status = redirect_stdout(() -> runpf!(net, 25, 1e-8, 0), devnull)
  @printf("step 57 by hand: converged %s in %d iterations, lowest voltage %.4f pu\n", status == 0, iterations, minimum(nd._vm_pu for nd in net.nodeVec))
  println()

  # 3. every step on its own (QSTS mode A): the base case first, then each
  #    step from the base solution, taps as in the case
  config = SparlectraConfig(Dict{String,Any}("qsts" => Dict{String,Any}("mode" => "snapshots")))
  result = redirect_stdout(() -> runQsts(series; config = config), devnull)
  println(length(result.steps), " steps solved, ", count(s -> s.status === :solved, result.steps), " converged, base case ", result.base.status)

  # 4. the same with the tap controller of the case in every step
  controlled = SparlectraConfig(Dict{String,Any}("qsts" => Dict{String,Any}("mode" => "snapshots", "snapshot_controls" => true)))
  with_taps = redirect_stdout(() -> runQsts(series; config = controlled), devnull)
  tap_branch = findfirst(k -> !isnan(with_taps.branch_values[k, findfirst(==("tap_ratio"), QSTS_BRANCH_QUANTITIES), 1]), eachindex(with_taps.elements.branch_ids))
  tap_id = with_taps.elements.branch_ids[tap_branch]

  vm = findfirst(==("vm_pu"), QSTS_BUS_QUANTITIES)
  loading = findfirst(==("loading_pct"), QSTS_BRANCH_QUANTITIES)
  losses = qstsValues(result, "", "total_p_loss_mw")
  taps = qstsValues(with_taps, tap_id, "tap_ratio")
  @printf("%-6s %-6s %9s %11s %11s %9s\n", "step", "time", "Vmin[pu]", "loading[%]", "losses[MW]", "tap " * with_taps.elements.branch_names[tap_branch])
  for k in 1:4:length(result.steps)
    vmin = minimum(filter(!isnan, result.bus_values[:, vm, k]))
    @printf("%-6d %-6s %9.4f %11.1f %11.3f %9.4f\n", k, result.steps[k].label, vmin, _step_max(result.branch_values[:, loading, k]), losses[k], taps[k])
  end
  println()

  # 5. the result files: long format, steps, elements, manifest
  out = mktempdir()
  files = writeQstsResults(out, result)
  println("written to ", out, ": ", join(basename.(files), ", "), " (", length(qstsLongTable(result).step), " rows in the long format)")
  println()

  # 6. the example controls: the case file carries the tap controller of
  #    Moorau_110_20, sp_case14_day.control.yaml gives it the control range
  #    tap_min 1.0 (not below the neutral ratio), adds a tap controller for
  #    Wehrden_110_20 and makes the shunt at Moorau_20 a switched bank
  #    (blocks of 2.5 MVAr, 0 is the case value) at priority 1 behind the
  #    Moorau tap, target 1.01 pu, deadband 0.02 pu. Applied once to a fresh
  #    base net, they travel with every copy. Each step solves with the
  #    outer control loop: from 18:00 to 20:30 the Moorau tap stands at 1.0
  #    and the bank switches in one block after it; off peak the bank
  #    stays off.
  cfg = load_sparlectra_config(SERIES_CONTROL; reload = true)
  controlled = build_net(series)
  applyConfiguredControllers!(controlled, cfg.control)
  println("with the example controls (case taps 1.0 and 1.0, Moorau range from 1.0):")
  @printf("%-6s %-6s %11s %12s %11s %10s %11s %10s\n", "step", "time", "tap Moorau", "tap Wehrden", "bank[MVAr]", "V Moorau", "V Wehrden", "control")
  for k in (0, 1, 29, 57, 72, 73, 77, 83, 84, 88)
    step_net = deepcopy(controlled)
    k == 0 || apply!(step_net, series.ops[k], series.index)
    redirect_stdout(() -> run_sparlectra(net = step_net, config = cfg), devnull)
    taps = sort(buildTapControllerReportRows(step_net); by = r -> r.transformer_branch_index)
    bank = only(r for r in latest_control_result(step_net).controllers if hasproperty(r, :bs_mvar))
    @printf("%-6s %-6s %11.5f %12.5f %11.1f %10.4f %11.4f %10s\n", k == 0 ? "base" : string(k), k == 0 ? "-" : series.steps[k].label, taps[1].tap_ratio, taps[2].tap_ratio, bank.bs_mvar,
      get_bus_vm_pu(step_net, "Moorau_20"), get_bus_vm_pu(step_net, "Wehrden_20"), string(latest_control_result(step_net).status))
  end
  println()

  # 7. the whole day in order with the same controllers (QSTS mode B): each
  #    step starts from the taps and the bank of the step before; step 1
  #    takes the values of the series. The switching list holds every net
  #    move, the metrics count them and find hunting.
  sequential = load_sparlectra_config(SERIES_CONTROL; reload = true, overrides = Dict("qsts" => Dict("mode" => "sequential")))
  day = redirect_stdout(() -> runQsts(series; config = sequential), devnull)
  println("the day in order (mode B): ", count(s -> s.status === :solved, day.steps), " of ", length(day.steps), " steps converged")
  for (id, name) in zip(day.elements.controller_ids, day.elements.controller_names)
    position = qstsValues(day, id, "position")
    blocks = String[]
    first_step = 1
    for k in 2:(length(position) + 1)
      (k <= length(position) && position[k] == position[k - 1]) && continue
      push!(blocks, @sprintf("%s-%s: %+d", day.steps[first_step].label, day.steps[k - 1].label, Int(position[k - 1])))
      first_step = k
    end
    println(rpad(name, 20), join(blocks, ", "))
  end
  sw = qstsSwitching(day)
  println("switching list:")
  for i in eachindex(sw.step)
    @printf("  %-6s %-10s %+d -> %+d  %s\n", day.steps[sw.step[i]].label, sw.element_id[i], sw.old[i], sw.new[i], sw.trigger[i])
  end
  m = qstsMetrics(day)
  value(id, metric) = m.value[findfirst(i -> m.element_id[i] == id && m.metric[i] == metric, eachindex(m.metric))]
  @printf("switch count %d, hunting patterns %d, loss energy %.2f MWh\n", value("", "switch_count_total"), value("", "n_hunting_total"), value("", "loss_energy_mwh"))
  files = writeQstsResults(mktempdir(), day)
  println("mode B writes ", join(basename.(files), ", "))
  # the voltage of Moorau_20 (bus id 8) over the day as SVG, with the band
  chart = qstsChart(day, "8", "vm_pu")
  svg = joinpath(mktempdir(), "moorau_20_vm.svg")
  write(svg, chart.svg)
  println("chart of Moorau_20 vm_pu (", count(!isnan, chart.values), " points, band ", chart.band, "): ", svg)
  return nothing
end

run_example(main)
