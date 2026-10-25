# Copyright 2023-2026 Udo Schmitz                                             #src
#                                                                             #src
# Licensed under the Apache License, Version 2.0 (the "License");             #src
# you may not use this file except in compliance with the License.            #src
# You may obtain a copy of the License at                                     #src
#                                                                             #src
#     http://www.apache.org/licenses/LICENSE-2.0                              #src
#                                                                             #src
# Unless required by applicable law or agreed to in writing, software         #src
# distributed under the License is distributed on an "AS IS" BASIS,           #src
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.    #src
# See the License for the specific language governing permissions and         #src
# limitations under the License.                                              #src
#                                                                             #src
# file: docs/lit/workshop_qsts.jl                                             #src
# purpose: Literate.jl source of the QSTS workshop notebook and its           #src
#          Documenter page: the shipped example series sp_case14_day         #src
#          (96 steps of 15 minutes) in mode A (snapshots) and mode B         #src
#          (sequential with the tap and bank controllers of                  #src
#          sp_case14_day.control.yaml), the switch qsts.mode, the accessors  #src
#          qstsValues, qstsSwitching and qstsMetrics, one chart, and what    #src
#          differs between the modes. Every printed number sits next to an   #src
#          @assert, and the test suite RUNS this file, so the notebook       #src
#          cannot drift silently. Regenerate the committed outputs with      #src
#          `julia --project=docs docs/generate_notebooks.jl`.                 #src

# # QSTS: a day on sp_case14 in two modes
#
# > **Level: Advanced**, builds on the OLTC chapter of the basic tour.
#
# [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_qsts.ipynb)
#
# > **Note:** This workshop was created with AI assistance and is reviewed
# > and curated by the maintainer; it is not a fully machine-generated text.
#
# A quasi-static time series (QSTS) solves one power flow per step of a
# series of operating states. Sparlectra has two modes. Mode A
# (`snapshots`) solves every step on its own. Mode B (`sequential`) solves
# the steps in order and carries the tap and bank positions from one step
# to the next, so the controllers act as they would over a day. This
# workshop runs the shipped example day on `sp_case14` in both modes and
# shows where the results differ.
#
# > **Note:** On Google Colab the install cell takes a minute or two on a
# > fresh session (package download and precompilation). Colab's Julia
# > version may change over time; this notebook targets Julia >= 1.12.

#nb # ## Setup (Colab)
#nb # This cell installs Sparlectra from GitHub (branch `main`) into a fresh
#nb # temporary environment. The isolation matters on Colab: the shared
#nb # default environment ships many preinstalled packages, and installing
#nb # anything there triggers precompilation of that whole stack. Run this
#nb # cell first, once per session; it takes a minute or two.
#nb # > **Colab is slow here.** The install cell below fetches Sparlectra and
#nb # > its dependencies and precompiles them on Colab's two cores; that takes a
#nb # > minute or two and prints one line per package as it goes. Sparlectra
#nb # > itself precompiles in well under a minute; the first cell that
#nb # > imports a case and solves it then compiles for about a minute more,
#nb # > every later cell runs at speed. Wait for the install cell to finish
#nb # > before running the next one; a second click on it starts the whole
#nb # > install again.
#nb using Pkg
#nb Pkg.activate(temp = true)
#nb ## no APSLF in this workshop: skip the precompile workload of AnalyticLoadFlow (set before Pkg.add, the precompile processes inherit ENV)
#nb ENV["ANALYTICLOADFLOW_PRECOMPILE_WORKLOAD"] = "off"
#nb Pkg.add(url = "https://github.com/Welthulk/Sparlectra.jl", rev = "main")
#nb ## To test another branch, set rev to its name, e.g. rev = "dev/r0.10.0".
#nb ## For the latest registered release use: Pkg.add("Sparlectra")

# ## Warm-up and the example series
#
# Julia compiles each function on first use. This cell loads the package
# and reads the series, so any chapter can run on its own after it; a code
# cell that uses a name from another cell says so in a comment. Its last
# line, `warmup()`, runs every path the chapters use once on the first two
# steps of the day; its code is in `docs/lit/warmup/workshop_qsts.jl` of
# the package. How long the warm-up takes depends on the machine.
#
# The series is `data/qsts/sp_case14_day.qsts.zip`: the case file
# `sp_case14.scf.json` as base and 96 steps of 15 minutes. Each step is a
# list of patch operations against the base case: the eight loads follow a
# daily shape, the two machines other than the reference follow the load.
# From 14:00 to 15:45 the line `Ostheim_Neuwiese_1` (component id 26) is
# out, one more operation in those eight steps. The taps and the shunt are
# not in the series.

using Sparlectra
using Printf

## the shipped series and the controller file next to it
data = joinpath(dirname(dirname(pathof(Sparlectra))), "data", "qsts")
series = read_scf_series(joinpath(data, "sp_case14_day.qsts.zip"))
control_file = joinpath(data, "sp_case14_day.control.yaml")
println("series ", series.name, ": ", length(series.steps), " steps of ", Int(series.step_minutes), " min (", series.steps[1].label, " to ", series.steps[end].label, "), ", sum(length, series.ops), " patch operations")
outage = only(filter(op -> op.op === :status, series.ops[57]))
println("step 57 (", series.steps[57].label, "), status operation: ", outage.target, " ", outage.id, " status set to ", outage.value, " (0 = out of service)")
@assert length(series.steps) == 96 && series.step_minutes == 15.0
@assert sum(length, series.ops) == 968
@assert [k for k in eachindex(series.ops) if any(op -> op.op === :status, series.ops[k])] == collect(57:64)

## compile every path the chapters use once; the warm-up code is in
## docs/lit/warmup/workshop_qsts.jl of the installed package
include(pkgdir(Sparlectra, "docs", "lit", "warmup", "workshop_qsts.jl"))
warmup()

# ## 1. The switch `qsts.mode`
#
# `runQsts` takes the run mode from the configuration key `qsts.mode`. The
# default is `off`, and with it `runQsts` refuses to run. The two other
# values are `snapshots` (mode A) and `sequential` (mode B).

## uses: series (warm-up cell)
println("default qsts.mode: ", SparlectraConfig().qsts.mode)
refusal = try
    runQsts(series; config = SparlectraConfig())
    "ran"
catch err
    err isa ArgumentError || rethrow()
    sprint(showerror, err)
end
println(refusal)
@assert SparlectraConfig().qsts.mode === :off
@assert occursin("qsts.mode is off", refusal)

# ## 2. Mode A: every step on its own
#
# Mode A solves the base case, then every step on a fresh copy of the base
# network with the delta of the step applied, started from the base
# solution. Nothing is carried from one step to the next, so the steps can
# run in parallel. Without `qsts.snapshot_controls = true` no controller
# runs: the taps stay where the case file puts them.
#
# `qstsValues(result, element_id, quantity)` returns one quantity of one
# element over all steps, `NaN` where a step has no value. On an SCF series
# the element id is the component id of the case file: bus 8 is
# `Moorau_20`, bus 14 `Wehrden_20`, branch 26 the line that goes out at
# 14:00. An empty id selects a step quantity.

## uses: series (warm-up cell)
config_a = SparlectraConfig(Dict{String,Any}("qsts" => Dict{String,Any}("mode" => "snapshots")))
plain = redirect_stdout(() -> runQsts(series; config = config_a), devnull)
labels = [s.label for s in plain.steps]
vm_moorau = qstsValues(plain, "8", "vm_pu")
vm_wehrden = qstsValues(plain, "14", "vm_pu")
line_out = findall(isnan, qstsValues(plain, "26", "p_from_mw"))
println("mode A: ", count(s -> s.status === :solved, plain.steps), " of ", length(plain.steps), " steps solved, controller status ", join(unique(s.control_status for s in plain.steps), ", "))
println("branch 26 Ostheim_Neuwiese_1, P at the from side: no value from ", labels[first(line_out)], " to ", labels[last(line_out)], " (out of service)")
println("bus voltage Vm [pu] over the day, mode A without controllers (power flow result):")
println("  bus                 lowest   highest   steps below 0.99 pu")
for (name, vm) in (("Moorau_20 (bus 8)", vm_moorau), ("Wehrden_20 (bus 14)", vm_wehrden))
    @printf("  %-19s %6.4f   %6.4f    %2d\n", name, minimum(vm), maximum(vm), count(<(0.99), vm))
end
println("  Moorau_20 below 0.99 pu from ", labels[findfirst(<(0.99), vm_moorau)], " to ", labels[findlast(<(0.99), vm_moorau)])
@assert all(s -> s.status === :solved && s.control_status === :none, plain.steps)
@assert line_out == collect(57:64)
@assert findall(<(0.99), vm_moorau) == collect(73:83)
@assert count(<(0.99), vm_wehrden) == 67 && maximum(vm_wehrden) < 1.0

# Without control the evening peak (18:00 to 20:30) pulls `Moorau_20`
# below 0.99 pu, and `Wehrden_20` stays below 1.0 pu all day.

# ## 3. Mode B: the day in order with controllers
#
# The file `sp_case14_day.control.yaml` holds the controllers for this
# day. The case file already has a tap controller on `Moorau_110_20`
# (target 1.0 pu at `Moorau_20`, deadband 0.01 pu); the file gives it the
# range `tap_min: 1.0`, adds a tap controller on `Wehrden_110_20` (1.0 pu
# at `Wehrden_20`) and turns the shunt at `Moorau_20` into a switched bank
# of 2.5 MVAr blocks behind the Moorau tap (`priority: 1`, target 1.01 pu,
# deadband 0.02 pu). `load_sparlectra_config` reads it, the override sets
# the mode.
#
# Each step starts from the solution of the step before and from the taps
# and bank of its end. A controller has the address `<actuator>#<kind>`:
# `30#OLTC` is the Moorau transformer, `31#OLTC` the Wehrden transformer,
# `8#MSC` the bank at bus 8. Its `position` is the tap step (0 at ratio
# 1.0) or, for the bank, the number of blocks.

## uses: control_file, series (warm-up cell)
config_b = load_sparlectra_config(control_file; reload = true, overrides = Dict("qsts" => Dict("mode" => "sequential")))
day = redirect_stdout(() -> runQsts(series; config = config_b), devnull)
println("mode B: ", count(s -> s.status === :solved, day.steps), " of ", length(day.steps), " steps solved")
println("controller position, mode B (OLTC: tap step, 0 at ratio 1.0; MSC: blocks in service):")
println("  controller  name               at 00:00   changed to")
for (id, name) in zip(day.elements.controller_ids, day.elements.controller_names)
    position = qstsValues(day, id, "position")
    changes = [k for k in 2:length(position) if position[k] != position[k - 1]]
    println("  ", rpad(id, 12), rpad(name, 19), lpad(Int(position[1]), 8), "   ", join([string(Int(position[k]), " at ", day.steps[k].label) for k in changes], ", "))
end
@assert day.mode === :sequential && all(s -> s.status === :solved, day.steps)
@assert day.elements.controller_ids == ["30#OLTC", "31#OLTC", "8#MSC"]
@assert qstsValues(day, "30#OLTC", "position") == vcat(fill(0.0, 3), fill(1.0, 23), fill(0.0, 62), fill(1.0, 8))
@assert qstsValues(day, "31#OLTC", "position") == vcat(fill(0.0, 25), fill(-1.0, 7), fill(-2.0, 64))
@assert qstsValues(day, "8#MSC", "bs_mvar") == vcat(fill(0.0, 72), fill(2.5, 24))

# The Moorau tap goes one step up at 00:45, back to its limit 1.0 at
# 06:30 and up again at 22:00. The Wehrden tap goes two steps down in the morning and stays
# there. At 18:00 the Moorau tap is at its limit and the voltage falls
# below 0.99 pu, so the bank switches in one block. It stays in for the
# rest of the day: with the block in, the voltage stays inside the band of
# the bank.
#
# `qstsSwitching` lists every step in which a position changed against the
# end of the step before, with old and new position and the trigger.
# `qstsMetrics` counts the moves, finds hunting (up, down, up within
# `qsts.hunting_window` steps) and sums the loss energy.

## uses: day (3. Mode B: the day in order with controllers)
@isdefined(day) || error("Run the section \"3. Mode B: the day in order with controllers\" first: it sets up day.")
sw = qstsSwitching(day)
println("switching list, mode B:")
println("  time   controller  position old -> new   trigger")
for i in eachindex(sw.step)
    println("  ", rpad(day.steps[sw.step[i]].label, 7), rpad(sw.element_id[i], 12), lpad(sw.old[i], 12), " -> ", rpad(sw.new[i], 6), sw.trigger[i])
end
metrics = qstsMetrics(day)
metric(id, name) = metrics.value[only(i for i in eachindex(metrics.metric) if metrics.element_id[i] == id && metrics.metric[i] == name)]
println("metrics, mode B:")
println("  switching operations, all controllers:     ", Int(metric("", "switch_count_total")))
println("  hunting patterns, all controllers:         ", Int(metric("", "n_hunting_total")))
@printf("  loss energy of the day [MWh]:              %.2f\n", metric("", "loss_energy_mwh"))
println("  30#OLTC Moorau, steps at the tap limit:    ", Int(metric("30#OLTC", "steps_at_limit")), " of ", length(day.steps))
@assert sw.step == [4, 26, 27, 33, 73, 89] && all(==("control"), sw.trigger)
@assert metric("", "switch_count_total") == 6 && metric("", "n_hunting_total") == 0
@assert metric("", "loss_energy_mwh") ≈ sum(qstsValues(day, "", "total_p_loss_mw")) / 4
@assert round(metric("", "loss_energy_mwh"); digits = 2) == 21.91
@assert metric("30#OLTC", "steps_at_limit") == 65

# Six moves in 96 steps and no hunting. The loss energy is the sum of the
# step losses times 15 minutes.

# ## 4. One chart
#
# `qstsChart` draws one quantity of one element over the steps as SVG, one
# point per step with a tooltip, the band dashed. Here `Wehrden_20` in
# mode B with the deadband of its tap controller, 0.99 to 1.01 pu. A
# notebook and the documentation show the returned chart directly; its
# field `svg` holds the text for a file.

## uses: day (3. Mode B: the day in order with controllers)
@isdefined(day) || error("Run the section \"3. Mode B: the day in order with controllers\" first: it sets up day.")
chart = qstsChart(day, "14", "vm_pu"; band = (0.99, 1.01))
println("Wehrden_20 (bus 14), Vm [pu], mode B (power flow result):")
@printf("  lowest of the day:                 %.4f\n", minimum(chart.values))
@printf("  highest of the day:                %.4f\n", maximum(chart.values))
@printf("  band of 31#OLTC (dashed):          %.2f to %.2f\n", chart.band...)
@assert count(!isnan, chart.values) == 96 && chart.band == (0.99, 1.01)
@assert all(v -> 0.99 <= v <= 1.01, chart.values)
chart

# In mode A without control the same bus stays below 0.99 pu in 67 steps
# (section 2). In mode B it stays inside the band all day.

# ## 5. What differs between the modes
#
# Mode A can run the controllers too: `qsts.snapshot_controls = true` runs
# the controllers of the configuration in every step, but every step
# starts again from the taps and the bank of the case file. We run it with
# the same controller file and compare with mode B.

## uses: control_file, series (warm-up cell); day (3. Mode B: the day in order with controllers)
@isdefined(day) || error("Run the section \"3. Mode B: the day in order with controllers\" first: it sets up day.")
config_ac = load_sparlectra_config(control_file; reload = true, overrides = Dict("qsts" => Dict("mode" => "snapshots", "snapshot_controls" => true)))
each = redirect_stdout(() -> runQsts(series; config = config_ac), devnull)
step_changes(v) = count(k -> v[k] != v[k - 1], 2:length(v))
println("step 1 (", day.steps[1].label, "), all bus values equal in mode A and mode B: ", isequal(each.bus_values[:, :, 1], day.bus_values[:, :, 1]))
println("                                                 mode A    mode B")
println("tap ratio changes between consecutive steps:")
for (id, name) in (("30", "Moorau_110_20 (30)"), ("31", "Wehrden_110_20 (31)"))
    @printf("  %-45s %6d    %6d\n", name, step_changes(qstsValues(each, id, "tap_ratio")), step_changes(qstsValues(day, id, "tap_ratio")))
end
k = 22
println(day.steps[k].label, ", power flow result:")
@printf("  %-45s %6.4f    %6.4f\n", "tap ratio Moorau_110_20 (30) [pu]", qstsValues(each, "30", "tap_ratio")[k], qstsValues(day, "30", "tap_ratio")[k])
@printf("  %-45s %6.4f    %6.4f\n", "Vm Moorau_20 (bus 8) [pu]", qstsValues(each, "8", "vm_pu")[k], qstsValues(day, "8", "vm_pu")[k])
k = 84
println(day.steps[k].label, ", power flow result:")
@printf("  %-45s %6.4f    %6.4f\n", "Vm Moorau_20 (bus 8) [pu]", qstsValues(each, "8", "vm_pu")[k], qstsValues(day, "8", "vm_pu")[k])
@printf("  %-45s %6s    %6.1f\n", "bank 8#MSC in service [MVAr]", "n/a", qstsValues(day, "8#MSC", "bs_mvar")[k])
println("whole day:")
@printf("  %-45s %6d    %6d\n", "Newton iterations, sum over 96 steps", sum(s.iterations for s in each.steps), sum(s.iterations for s in day.steps))
@printf("  %-45s %6d    %6d\n", "switching rows", length(each.switching), length(day.switching))
@printf("  %-45s %6d    %6d\n", "metric rows", length(each.metrics), length(day.metrics))
@assert isequal(each.bus_values[:, :, 1], day.bus_values[:, :, 1]) && isequal(each.branch_values[:, :, 1], day.branch_values[:, :, 1])
@assert step_changes(qstsValues(each, "30", "tap_ratio")) == 2 && step_changes(qstsValues(day, "30", "tap_ratio")) == 3
@assert step_changes(qstsValues(each, "31", "tap_ratio")) == 10 && step_changes(qstsValues(day, "31", "tap_ratio")) == 2
@assert (qstsValues(each, "30", "tap_ratio")[22], qstsValues(day, "30", "tap_ratio")[22]) == (1.0, 1.0125)
@assert all(v -> 0.99 <= v <= 1.01, (qstsValues(each, "8", "vm_pu")[22], qstsValues(day, "8", "vm_pu")[22]))
@assert 0.99 <= qstsValues(each, "8", "vm_pu")[84] < 1.0 && qstsValues(day, "8", "vm_pu")[84] > 1.0
@assert sum(s.iterations for s in each.steps) > 2 * sum(s.iterations for s in day.steps)
@assert isempty(each.switching) && isempty(each.metrics)

# Step 1 is the same in both modes: mode B starts the day from the base
# solution, as mode A starts every step. After that the history decides.
#
# - **A deadband allows more than one position.** At 05:15 both Moorau
#   positions hold `Moorau_20` inside 0.99 to 1.01 pu. Mode A starts from
#   the case tap 1.0 and finds no reason to move; mode B still stands at
#   1.0125 from 00:45 and finds no reason either.
# - **In mode A the position follows the load.** The Wehrden tap changes
#   between consecutive steps ten times in mode A and twice in mode B. A
#   change in mode A is no switching operation: each step searches its
#   position from the tap of the case file.
# - **Mode B keeps the bank in.** At 20:45 mode A starts with the bank off;
#   `Moorau_20` at 0.991 pu lies inside the band of the bank, so it stays
#   off. Mode A reports no controller values, so its bank column reads
#   n/a. Mode B carries the block from 18:00 into the night, and the
#   voltage is above 1.0 pu.
# - **Mode B is cheaper per step.** It starts each step from the step
#   before; mode A needs more than twice the Newton iterations.
# - **Only mode B counts.** Switching list and metrics exist in mode B
#   only; `writeQstsResults` writes `qsts_switching.csv` and
#   `qsts_metrics.csv` for it.
#
# | | Mode A, `snapshots` | Mode B, `sequential` |
# |---|---|---|
# | steps | independent, parallel | in order; several series in parallel |
# | controllers | off, or each step from the base case (`qsts.snapshot_controls`) | the `control` section, taps and banks carried |
# | start of a step | base solution | the step before |
# | results | values per step | values, switching list, metrics |
# | use | the set points of the source system per step | controller behaviour over a day |

# ## 6. The result files
#
# `writeQstsResults` writes the long table, the steps, the elements and a
# manifest. Mode B adds the switching list and the metrics.

## uses: day (3. Mode B: the day in order with controllers); plain (2. Mode A: every step on its own)
@isdefined(day) || error("Run the section \"3. Mode B: the day in order with controllers\" first: it sets up day.")
@isdefined(plain) || error("Run the section \"2. Mode A: every step on its own\" first: it sets up plain.")
## format="technical" fixes comma and decimal point; the default follows the regional settings
files_a = writeQstsResults(mktempdir(), plain; format = "technical")
files_b = writeQstsResults(mktempdir(), day; format = "technical")
println("mode A: ", join(basename.(files_a), ", "))
println("mode B: ", join(basename.(files_b), ", "))
@assert basename.(files_a) == ["qsts_values.csv", "qsts_steps.csv", "qsts_elements.csv", "qsts_manifest.json"]
@assert setdiff(basename.(files_b), basename.(files_a)) == ["qsts_switching.csv", "qsts_metrics.csv"]

# ## Where to go from here
#
# - [QSTS: Snapshot Series](../qsts.md): the modes, the series format, the
#   result quantities and the configuration keys.
# - The example `examples/others/exp_scf_series.jl` runs the same day as a
#   script and solves one step by hand with `apply!` and `runpf!`.
# - In the Web UI the button **QSTS** runs a series file and shows the
#   chart of any bus or branch.
