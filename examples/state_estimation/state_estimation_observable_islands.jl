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

# Date: 2026-11-04
# file: examples/state_estimation/state_estimation_observable_islands.jl
# purpose: a measurement set with a region without meters on the shipped sp_case118 and sp_case1354:
#          observable islands, estimation per observable island, closing with scheduled values,
#          measurement placement of each kind, and the suggested meters built in until the set is
#          observable; the run time of each step goes to the console and to examples/_out

using Sparlectra
using Printf
using Random
using Dates

include(joinpath(@__DIR__, "..", "others", "example_header.jl"))

const OUTDIR_ISL = joinpath(@__DIR__, "..", "_out")
const MPOWER_DIR_ISL = joinpath(@__DIR__, "..", "..", "data", "mpower")

# Fixed regions without any meter: every bus row of these buses and every
# flow row of a branch that touches one of them is removed. The noise of the
# remaining rows is seeded, so every run prints the same numbers.
const REGIONS_ISL = Dict("sp_case118.m" => 98:107, "sp_case1354.m" => 50:60)

function _load_case_isl(file::AbstractString)
  net = Sparlectra.createNetFromMatPowerFile(filename = joinpath(MPOWER_DIR_ISL, file), flatstart = false, bus_shunt_model = :admittance, matpower_shift_sign = 1.0, matpower_shift_unit = :deg, matpower_ratio = :normal, tap_changer_model = :ideal)
  _, erg = runpf!(net, 60, 1e-8, 0)
  erg == 0 || error("$(file): reference power flow did not converge")
  return net
end

# The measurement set: Vm, P and Q at every bus and P and Q flows at both
# branch ends from the power flow, without the rows of the region. The
# removed rows are returned as the pool of meters that can be put back.
function _gap_set_isl(net::Net, region)
  ms = generateMeasurementsFromPF(net; includeVm = true, includePinj = true, includeQinj = true, includePflow = true, includeQflow = true, noise = true, rng = MersenneTwister(1))
  touches(m) = (m.busIdx !== nothing && m.busIdx in region) || (m.branchIdx !== nothing && (Int(net.branchVec[m.branchIdx].fromBus) in region || Int(net.branchVec[m.branchIdx].toBus) in region))
  return Measurement[m for m in ms if !touches(m)], Measurement[m for m in ms if touches(m)]
end

# runse! on a copy of the net with one value of on_unobservable
function _estimate_isl(net::Net, meas::Vector{Measurement}, mode::Symbol)
  work = deepcopy(net)
  res = with_state_estimation_config(update_net = true, observability = Sparlectra.ObservabilityConfig(on_unobservable = mode)) do
    runse!(work, meas)
  end
  return res, work
end

# One timed step: runs `f`, prints and records its run time.
function _timed_isl(f, times::Vector{Pair{String,Float64}}, label::AbstractString)
  t = @elapsed r = f()
  push!(times, label => t)
  @printf("  [%s: %.2f s]\n", label, t)
  return r
end

# The chain on one case. `verbose` prints the full island and placement
# reports (sp_case118); otherwise one line per step (sp_case1354).
function _gap_chain_isl(file::AbstractString; verbose::Bool)
  times = Pair{String,Float64}[]
  net = _load_case_isl(file)
  region = REGIONS_ISL[file]
  meas, pool = _gap_set_isl(net, region)
  @printf("%s: %d buses, %d of %d measurement rows, no meter on buses %d to %d and their branches\n", file, length(net.nodeVec), length(meas), length(meas) + length(pool), first(region), last(region))

  # 1. observable islands, deficits, unobservable branches
  obs = _timed_isl(times, "observable islands") do
    evaluate_global_observability(net, meas)
  end
  isl = obs.observable_islands
  obs.quality === :not_observable || error("$(file): the region must leave the set not observable")
  if verbose
    print_observable_islands(isl)
  else
    @printf("  %d observable islands (largest %d buses), rank deficit %d (angles %d, magnitudes %d), %d unobservable branches\n", isl.n_islands, maximum(length, isl.island_buses), isl.rank_deficit, isl.angle_rank_deficit, isl.magnitude_rank_deficit, length(isl.unobservable_branches))
  end

  # 2. every observable island estimated on its own
  se_i, _ = _timed_isl(times, "estimate per observable island") do
    _estimate_isl(net, meas, :observable_islands)
  end
  se_i.converged || error("$(file): the estimate per observable island did not converge")
  est = [r for r in se_i.islands if r.estimated]
  println("  estimated per observable island: ", length(est), " of ", length(se_i.islands), " islands, ", count(.!se_i.observableIslands.estimated), " bus(es) not estimated")
  for r in se_i.islands
    if r.estimated
      @printf("    island %-5s %4d bus(es), reference %s%s, J %.2f, dof %d\n", r.label, length(r.buses), isl.bus_names[r.reference_bus], r.relative_angles ? " (relative angles)" : "", r.objectiveJ, r.dof)
    else
      @printf("    island %-5s %4d bus(es), not estimated (%s)\n", r.label, length(r.buses), replace(string(r.reason), "_" => " "))
    end
  end

  # 3. closing with scheduled values first, the rest per observable island
  se_s, _ = _timed_isl(times, "estimate with scheduled values") do
    _estimate_isl(net, meas, :scheduled)
  end
  se_s.converged || error("$(file): the estimate with scheduled values did not converge")
  println("  restored rows: ", length(se_s.restoredRows))
  for r in se_s.restoredRows
    @printf("    %s at %s: %.1f, sigma %.1f, from %s\n", r.quantity, r.bus_name, r.value, r.sigma, r.source)
  end
  rest = restore_observability(net, meas)
  println("  verdict after restoring: ", replace(string(rest.verdict), "_" => " "), rest.verdict === :weakly_observable ? string(" (at the stored state: ", length(unique(rest.weak.branches)), " branch(es) above the flow bound, ", length(rest.weak.buses), " bus(es) affected)") : "")
  println("  the estimate: ", se_s.observableIslands === nothing ? "one piece" : string(count(r -> r.estimated, se_s.islands), " observable islands"))

  # 4. placement of each kind, 5. the suggested meters built in
  final = Dict{Symbol,Symbol}()
  for kind in (:injection, :flow, :pmu)
    pl = _timed_isl(times, "placement $(kind)") do
      suggest_measurement_placement(net, meas; kind = kind)
    end
    if verbose
      print_measurement_placement(pl)
    else
      @printf("  placement %s: %d steps, target %s\n", kind, length(pl.steps), pl.reached ? "reached" : "not reached ($(pl.stop_reason))")
    end
    built = vcat(meas, pl.rows)
    if !pl.reached
      # a PMU is an angle only: a magnitude gap needs another kind
      more = suggest_measurement_placement(net, built; kind = :injection)
      println("    plus ", length(more.steps), " injection step(s) for what ", kind, " does not close")
      built = vcat(built, more.rows)
    end
    o = evaluate_global_observability(net, built)
    println("    with the suggested meters: ", o.quality)
    round = 0
    while o.quality === :weakly_observable && round < 3
      # the placement closes the deficit, the weak class remains: put back
      # the removed meters there, the flow pair of every branch above the
      # flow bound and the Vm row of every affected bus
      round += 1
      wk = o.weak_observability
      back = Measurement[m for m in pool if (m.typ in (Sparlectra.PflowMeas, Sparlectra.QflowMeas) && m.branchIdx in wk.branches) || (m.typ == Sparlectra.VmMeas && m.busIdx in wk.buses)]
      built = vcat(built, back)
      o = evaluate_global_observability(net, built)
      println("    plus ", length(back), " meter rows at the weak places (", length(unique(wk.branches)), " branch(es), ", length(wk.buses), " bus(es)): ", o.quality)
    end
    o.quality in (:good, :critical) || error("$(file): the suggested meters of kind $(kind) leave the set $(o.quality)")
    final[kind] = o.quality
  end
  return (file = file, region = region, n_meas = length(meas), islands = isl.n_islands, deficit = isl.rank_deficit, restored = length(se_s.restoredRows), restore_verdict = rest.verdict, final = final, times = times)
end

function _write_report_isl(results)
  mkpath(OUTDIR_ISL)
  path = joinpath(OUTDIR_ISL, "state_estimation_observable_islands_$(Dates.format(now(), "yyyymmdd_HHMMSS")).md")
  open(path, "w") do io
    println(io, "# Observable islands, restoration and placement")
    println(io)
    for r in results
      println(io, "## ", r.file)
      println(io)
      println(io, "No meter on buses ", first(r.region), " to ", last(r.region), "; ", r.n_meas, " measurement rows; ", r.islands, " observable islands; rank deficit ", r.deficit, "; ", r.restored, " restored rows (", r.restore_verdict, "); with the suggested meters: ", join(["$(k) $(v)" for (k, v) in sort(collect(r.final))], ", "), ".")
      println(io)
      println(io, "| step | seconds |")
      println(io, "|---|---:|")
      for (label, t) in r.times
        @printf(io, "| %s | %.3f |\n", label, t)
      end
      println(io)
    end
  end
  return path
end

function run_observable_islands()
  print_example_banner("examples/state_estimation/state_estimation_observable_islands.jl", "observable islands, estimation per island, scheduled values and measurement placement on sp_case118 and sp_case1354")
  println("The run time of the first call of a step includes its compilation.")
  println()
  results = Any[]
  push!(results, _gap_chain_isl("sp_case118.m"; verbose = true))
  println()
  push!(results, _gap_chain_isl("sp_case1354.m"; verbose = false))
  path = _write_report_isl(results)
  println()
  println("report: ", normpath(path))
  return results
end

run_example(run_observable_islands)
