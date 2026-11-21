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
# file: benchmark/gpu/batch_study_summary.jl
# purpose: GPU study, stage 2: the Markdown tables of the report from the
#          CSVs of batch_study_gpu.jl (cold and warm series). Throughput
#          per case and workload for C1, C2 and the GPU variants; speedup of
#          the GPU default against the faster CPU path, end to end and
#          against the CPU steady state (without the serial creation of the
#          working copies and the warm-up of the power-mode workers); the
#          break-even network size per workload and the break-even N per
#          case; the CPU task-count scan; device memory and its limits; the
#          mask on and off; first call against the median; correctness and
#          repeat determinism against the C1 batch; the largest final
#          mismatch; the share of scenarios with a Q-limit violation. No
#          package beyond Base.
#
# Run: julia benchmark/gpu/batch_study_summary.jl <batch_cold.csv> [batch_warm.csv ...]

const GPU = ["sequential", "sequential_fast", "sequential_f32", "nonuniform", "uniform", "uniform_nomask"]
const GPU_LABEL = Dict("sequential" => "GPU default", "sequential_fast" => "seq. fast", "sequential_f32" => "seq. f32", "nonuniform" => "nonuniform", "uniform" => "uniform", "uniform_nomask" => "uniform, no mask")

# The status field (column 8) of the memory-limit rows written at e57d363b
# holds a comma; extra fields are joined back into it.
function read_rows(path)
  lines = filter(l -> !isempty(l) && !startswith(l, "#"), readlines(path))
  head = split(lines[1], ",")
  rows = Dict{String,String}[]
  for l in lines[2:end]
    l == lines[1] && continue
    f = split(l, ",")
    extra = length(f) - length(head)
    extra < 0 && continue
    extra > 0 && (f = vcat(f[1:7], join(f[8:8+extra], ","), f[9+extra:end]))
    push!(rows, Dict{String,String}(zip(head, f)))
  end
  return rows
end

num(r, k) = (v = get(r, k, ""); isempty(v) ? NaN : parse(Float64, v))
fmt(x; d = 1) = isnan(x) ? "-" : string(round(x; digits = d))
sci(x) = isnan(x) ? "-" : (x == 0 ? "0" : (e = floor(Int, log10(abs(x))); string(round(x / 10.0^e; digits = 1), "e", e)))
ok(r) = r !== nothing && r["status"] == "ok"

function main()
  rows = reduce(vcat, read_rows.(ARGS))
  get_row(c, s, w, v) = (i = findfirst(r -> r["case"] == c && r["start"] == s && r["workload"] == w && r["variant"] == v, rows); i === nothing ? nothing : rows[i])
  cases = unique(r["case"] for r in rows)
  starts = unique(r["start"] for r in rows)
  steady(r) = num(r, "steady_median_s")

  println("## CPU task-count scan, N = 64 (scenarios per second, end to end / steady)\n")
  tasks = sort(unique(parse(Int, r["tasks"]) for r in rows if r["workload"] == "scan64"))
  println("| Case | variant | ", join(string.(tasks, " tasks"), " | "), " | fastest |")
  println("|---|---|", repeat("---|", length(tasks) + 1))
  for c in cases, v in ("C1", "C2")
    sr = [r for r in rows if r["case"] == c && r["workload"] == "scan64" && r["variant"] == v]
    isempty(sr) && continue
    cell(t) = (i = findfirst(r -> r["tasks"] == string(t), sr); i === nothing ? "-" : string(fmt(num(sr[i], "scen_per_s")), " / ", fmt(64 / steady(sr[i]))))
    best = sr[argmin([num(r, "median_s") for r in sr])]["tasks"]
    println("| ", c, " | ", v, " | ", join(cell.(tasks), " | "), " | ", best, " |")
  end

  for s in starts
    workloads = unique(r["workload"] for r in rows if r["start"] == s && r["workload"] != "scan64")
    println("\n## Start ", s, ": throughput in scenarios per second (end to end, median)\n")
    for w in workloads
      gv = filter(v -> any(r -> r["start"] == s && r["workload"] == w && r["variant"] == v, rows), GPU)
      println("### ", w, "\n")
      println("| Case | C1 (tasks) | C2 (tasks) | ", join([GPU_LABEL[v] for v in gv], " | "), " | speedup default / faster CPU | same, CPU steady |")
      println("|---|---|---|", repeat("---|", length(gv) + 2))
      for c in cases
        c1 = get_row(c, s, w, "C1")
        c2 = get_row(c, s, w, "C2")
        if !ok(c1)
          sk = get_row(c, s, w, "all")
          sk === nothing || println("| ", c, " | ", sk["status"], " |", repeat(" |", length(gv) + 3))
          continue
        end
        cell(r) = r === nothing ? "-" : ok(r) ? fmt(num(r, "scen_per_s")) : startswith(r["status"], "memory limit") ? "out of memory" : first(r["status"], 30)
        g = get_row(c, s, w, "sequential")
        gs = ok(g) ? num(g, "median_s") : NaN
        cpu_best = min(num(c1, "median_s"), num(c2, "median_s"))
        cpu_steady = min(steady(c1), steady(c2))
        gcells = [begin
          rs = [r for r in rows if r["case"] == c && r["start"] == s && r["workload"] == w && r["variant"] == v]
          i = findfirst(ok, rs)
          i !== nothing ? cell(rs[i]) : isempty(rs) ? "-" : any(r -> startswith(r["status"], "memory limit"), rs) ? "out of memory" : first(rs[1]["status"], 30)
        end for v in gv]
        println("| ", c, " | ", fmt(num(c1, "scen_per_s")), " (", c1["tasks"], ") | ", fmt(num(c2, "scen_per_s")), " (", c2["tasks"], ") | ", join(gcells, " | "), " | ", fmt(cpu_best / gs; d = 2), " | ", fmt(cpu_steady / gs; d = 2), " |")
      end
      for (label, st) in (("end to end", false), ("CPU steady", true))
        sp = Tuple{String,Float64}[]
        for c in cases
          c1 = get_row(c, s, w, "C1")
          g = get_row(c, s, w, "sequential")
          (ok(c1) && ok(g)) || continue
          c2 = get_row(c, s, w, "C2")
          cpu = st ? min(steady(c1), steady(c2)) : min(num(c1, "median_s"), num(c2, "median_s"))
          push!(sp, (c, cpu / num(g, "median_s")))
        end
        k = findfirst(i -> all(x -> x[2] >= 1.0, sp[i:end]), eachindex(sp))
        println("\nBreak-even ", label, " (GPU default at least as fast as the faster CPU path from this case on): ", k === nothing ? "none" : sp[k][1])
      end
      println()
    end

    println("### Start ", s, ": break-even N per case (random workloads, GPU default against the faster CPU path)\n")
    println("| Case | speedup per N, end to end | smallest N from which the GPU is at least as fast | same, CPU steady |")
    println("|---|---|---|---|")
    for c in cases
      ns = sort(unique(parse(Int, r["N"]) for r in rows if r["case"] == c && r["start"] == s && startswith(r["workload"], "random") && ok(r) && r["variant"] == "C1"))
      sp = Tuple{Int,Float64,Float64}[]
      for N in ns
        w = string("random", N)
        c1 = get_row(c, s, w, "C1")
        c2 = get_row(c, s, w, "C2")
        g = get_row(c, s, w, "sequential")
        ok(g) || continue
        push!(sp, (N, min(num(c1, "median_s"), num(c2, "median_s")) / num(g, "median_s"), min(steady(c1), steady(c2)) / num(g, "median_s")))
      end
      be(j) = (k = findfirst(i -> all(x -> x[j] >= 1.0, sp[i:end]), eachindex(sp)); k === nothing ? "none" : string(sp[k][1]))
      println("| ", c, " | ", join([string(x[1], ": ", fmt(x[2]; d = 2)) for x in sp], ", "), " | ", be(2), " | ", be(3), " |")
    end
  end

  println("\n## GPU memory (GiB in use on the device, peak of the loop) and limits\n")
  println("| Case | variant | largest peak (workload) | memory limit |")
  println("|---|---|---|---|")
  for c in cases, v in GPU
    rs = [r for r in rows if r["case"] == c && r["variant"] == v && ok(r)]
    isempty(rs) && continue
    i = argmax([num(r, "peak_gpu_gb") for r in rs])
    lim = [replace(r["status"], "memory limit: " => "") for r in rows if r["case"] == c && r["variant"] == v && startswith(r["status"], "memory limit")]
    nmax = maximum(parse(Int, r["N"]) for r in rs)
    println("| ", c, " | ", v, " | ", rs[i]["peak_gpu_gb"], " (", rs[i]["start"], " ", rs[i]["workload"], ") | ", isempty(lim) ? string("not reached up to N = ", nmax) : join(lim, "; "), " |")
  end

  println("\n## GPU default: setup, device loop, first call (median seconds)\n")
  println("| Case | start | workload | total | setup on the CPU | device loop | first call | refactorizations | reordering | matching |")
  println("|---|---|---|---|---|---|---|---|---|---|")
  for c in cases, r in rows
    (r["case"] == c && r["variant"] == "sequential" && ok(r)) || continue
    println("| ", c, " | ", r["start"], " | ", r["workload"], " | ", fmt(num(r, "median_s"); d = 4), " | ", fmt(num(r, "setup_median_s"); d = 4), " | ", fmt(num(r, "gpu_loop_median_s"); d = 4), " | ", fmt(num(r, "first_s"); d = 4), " | ", r["refactorizations"], " | ", r["reordering_alg"], " | ", r["matching_alg"], " |")
  end

  println("\n## Uniform batch with and without the mask (scenarios per second, refactorizations)\n")
  println("| Case | workload | with mask | without mask |")
  println("|---|---|---|---|")
  for c in cases, r in rows
    (r["case"] == c && r["variant"] == "uniform" && ok(r)) || continue
    nm = get_row(c, r["start"], r["workload"], "uniform_nomask")
    println("| ", c, " | ", r["workload"], " | ", fmt(num(r, "scen_per_s")), " (", r["refactorizations"], ") | ", ok(nm) ? string(fmt(num(nm, "scen_per_s")), " (", nm["refactorizations"], ")") : "-", " |")
  end

  println("\n## Correctness against the C1 batch (largest over all workloads of a start)\n")
  println("| Case | start | variant | max abs dVm pu | max abs dVa deg | convergence equal | iterations equal | Q-limit counts equal | repeats bitwise (distinct results min to max) | largest final mismatch pu | C1 = single runpf! |")
  println("|---|---|---|---|---|---|---|---|---|---|---|")
  for c in cases, s in starts, v in vcat(["C1", "C2"], GPU)
    rs = [r for r in rows if r["case"] == c && r["start"] == s && r["variant"] == v && ok(r) && r["workload"] != "scan64"]
    isempty(rs) && continue
    single = v == "C1" ? string(all(r -> r["single_solves_bitwise"] == "true", rs)) : ""
    d = [parse(Int, r["distinct_results"]) for r in rs]
    println("| ", c, " | ", s, " | ", v, " | ", sci(maximum(num(r, "max_dvm_vs_c1") for r in rs)), " | ", sci(maximum(num(r, "max_dva_deg_vs_c1") for r in rs)), " | ",
      all(r -> r["conv_equal_c1"] == "true", rs), " | ", all(r -> r["iter_equal_c1"] == "true", rs), " | ", all(r -> r["qviol_equal_c1"] == "true", rs), " | ",
      all(r -> r["repeat_bitwise"] == "true", rs), " (", minimum(d), " to ", maximum(d), ") | ", sci(maximum(num(r, "final_mismatch_max") for r in rs)), " | ", single, " |")
  end

  println("\n## Convergence and Q limits (C1 batch)\n")
  println("| Case | start | workload | converged | iterations | share with a Q-limit violation | violating PV buses per scenario, mean (max) |")
  println("|---|---|---|---|---|---|---|")
  for c in cases, r in rows
    (r["case"] == c && r["variant"] == "C1" && ok(r) && r["workload"] != "scan64") || continue
    println("| ", c, " | ", r["start"], " | ", r["workload"], " | ", r["converged"], " of ", r["N"], " | ", r["iter_min"], " to ", r["iter_max"], " | ", fmt(100 * num(r, "qviol_share"); d = 0), " % | ", fmt(num(r, "qviol_mean")), " (", r["qviol_max"], ") |")
  end

  tp = [r for r in rows if r["variant"] in ("g5", "g5r", "g6") && ok(r)]
  if !isempty(tp)
    println("\n## Two-phase variants (Float32 phase 1, Float64 phase 2)\n")
    println("Speedup = median time of C1 or of the GPU float32_mixed default (seq. f32) / median time of the variant, same case, start and workload.\n")
    println("| Case | workload | variant | switch | scen/s | speedup vs C1 | vs seq. f32 | phase-1 steps mean (max) | phase-2 steps mean (max) | switch by tol / stagnation / failed | phase-2 refactorizations | to Float64 | phase 1 s | phase 2 s | converged | max abs dVm vs C1 | largest final mismatch | max abs Y |")
    println("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for r in tp
      c1 = get_row(r["case"], r["start"], r["workload"], "C1")
      f32 = get_row(r["case"], r["start"], r["workload"], "sequential_f32")
      p2 = r["variant"] == "g6" ? num(r, "cpu_phase2_median_s") : num(r, "phase2_median_s")
      println("| ", r["case"], " | ", r["workload"], " | ", r["variant"], " | ", r["switch_tol"], " | ", fmt(num(r, "scen_per_s")), " | ",
        ok(c1) ? fmt(num(c1, "median_s") / num(r, "median_s"); d = 2) : "-", " | ", ok(f32) ? fmt(num(f32, "median_s") / num(r, "median_s"); d = 2) : "-", " | ",
        fmt(num(r, "p1_steps_mean")), " (", r["p1_steps_max"], ") | ", fmt(num(r, "p2_steps_mean")), " (", r["p2_steps_max"], ") | ",
        r["switch_by_tol"], " / ", r["switch_by_stagnation"], " / ", r["phase1_failed"], " | ", r["p2_refactorizations"], " | ", r["to_float64"], " | ",
        fmt(num(r, "phase1_median_s"); d = 4), " | ", fmt(p2; d = 4), " | ", r["converged"], " of ", r["N"], " | ", sci(num(r, "max_dvm_vs_c1")), " | ", sci(num(r, "final_mismatch_max")), " | ", fmt(num(r, "max_abs_y"); d = 0), " |")
    end
  end

  other = [r for r in rows if !ok(r) && !startswith(r["status"], "memory limit") && !startswith(r["status"], "skipped")]
  if !isempty(other)
    println("\n## Refused or fallen back\n")
    for r in other
      println("- ", r["case"], " ", r["start"], " ", r["workload"], " ", r["variant"], ": ", r["status"])
    end
  end
end

Base.invokelatest(main)
