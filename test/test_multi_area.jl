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

# file: test/test_multi_area.jl
# purpose: multi-area state estimation, two-level method (1.2.0): the area
#          partition against the area column of the MATPOWER file, the
#          link-group rule, the border information against the Schur
#          complement of the integrated gain matrix, and the estimate
#          against the power flow (noise-free rows) and against runse!
#          on the same noisy rows.

using Test
using Sparlectra
using Random
using LinearAlgebra
using Logging

# The MATPOWER area column (bus table, column 7) read straight from the
# file: the independent truth of the partition.
function _ma_file_areas(path::AbstractString)
  areas = Dict{String,Int}()
  inbus = false
  for line in eachline(path)
    s = strip(line)
    if startswith(s, "mpc.bus = [")
      inbus = true
      continue
    end
    inbus || continue
    startswith(s, "];") && break
    f = split(replace(s, ";" => ""))
    length(f) >= 7 && (areas[f[1]] = parse(Int, f[7]))
  end
  return areas
end

function _ma_case118()
  path = joinpath(pkgdir(Sparlectra), "data", "mpower", "sp_case118.m")
  net = Sparlectra.import_case(path, Sparlectra.load_sparlectra_config()).net
  runpf!(net, 60, 1e-10, 0)
  return net, path
end

# The area runs of one pass on worker threads against the serial loop. Run
# in the se profile and, repeated with four threads, in the threads profile
# (tools/run_gates.sh threads).
function run_multi_area_parallel_tests()
  @testset "multi-area estimation, parallel equals serial" begin (function ()
    net, _ = _ma_case118()
    meas = generateMeasurementsFromPF(net; noise = true, rng = MersenneTwister(5))
    function _ma_with(parallel::Bool)
      active = Sparlectra.active_sparlectra_config()
      pf = Sparlectra._copy_powerflow_with(active.powerflow; islands_parallel_min_buses = 1)
      rt = Sparlectra.RuntimeConfig(parallel = Sparlectra.ParallelRuntimeConfig(enabled = parallel))
      cfg = Sparlectra._sparlectra_config_with(active; powerflow = pf, runtime = rt)
      return Sparlectra.with_sparlectra_config(() -> runse_multi_area!(deepcopy(net), meas; compare_integrated = false, update_net = false), cfg)
    end
    serial = _ma_with(false)
    for _ = 1:5
      par = _ma_with(true)
      @test par.voltages == serial.voltages
      @test par.objectiveJ == serial.objectiveJ
      @test par.border_change == serial.border_change
    end
    parallel_comparison_report("multi-area SE, three areas (5 repeats)")
  end)() end
end

function run_multi_area_tests()
  @testset "multi-area state estimation" begin
    @testset "partition of sp_case118 against the file" begin (function ()
      net, path = _ma_case118()
      part = area_partition(net)
      fileAreas = _ma_file_areas(path)
      @test part.areas == [1, 2, 3]
      # the bus names of the import are the MATPOWER bus numbers
      @test all(part.area_of_bus[b] == fileAreas[string(net.busOrigIdxDict[b])] for b in eachindex(net.nodeVec))
      for (k, br) in enumerate(net.branchVec)
        Sparlectra._branch_terminal_state(br) === :closed || continue
        crosses = part.area_of_bus[br.fromBus] != part.area_of_bus[br.toBus]
        @test crosses == (k in part.tie_branches)
      end
      @test part.border_buses == sort!(unique!(vcat([Int(net.branchVec[k].fromBus) for k in part.tie_branches], [Int(net.branchVec[k].toBus) for k in part.tie_branches])))
      # every active row has exactly one owner; a flow row measured at the
      # to end of a tie branch belongs to the area of that end
      meas = generateMeasurementsFromPF(net)
      @test all(m -> Sparlectra._se_area_owner(net, part, m) in part.areas, meas)
      k = first(part.tie_branches)
      mto = first(m for m in meas if m.branchIdx == k && m.direction === :to)
      @test Sparlectra._se_area_owner(net, part, mto) == part.area_of_bus[net.branchVec[k].toBus]
      @test_throws ArgumentError area_partition(net; areas = Dict("1" => 1))
      # the area list file: header, comment, three separators; bus 28 moved
      # to area 2. The same list through the config key partition.
      list = tempname() * ".csv"
      open(list, "w") do io
        println(io, "bus;area")
        println(io, "# bus 28 moved")
        for (k, nd) in enumerate(net.nodeVec)
          name = getCompName(nd.comp)
          sep = (";", ",", "\t")[mod1(k, 3)]
          println(io, name, sep, name == "28" ? 2 : fileAreas[name])
        end
      end
      fromList = area_partition(net; areas = read_area_list(list))
      b28 = net.busDict["28"]
      @test fromList.area_of_bus[b28] == 2 && count(fromList.area_of_bus .!= part.area_of_bus) == 1
      viaKey = with_state_estimation_config(multi_area = Sparlectra.MultiAreaConfig(partition = list)) do
        Sparlectra._se_configured_partition(net, Sparlectra.state_estimation_config(), nothing)
      end
      @test viaKey.area_of_bus == fromList.area_of_bus
      open(io -> println(io, "28;1"), list, "a")
      @test_throws "named twice" read_area_list(list)
    end)() end

    @testset "a closed link group is never cut" begin (function ()
      # A1-A2-A3-A4 ring with the coupler A3-A3a. The list puts A3a in area
      # 1 and A3 in area 2: the group goes to the area of its representative
      # A3 (lowest bus index), with a note.
      net = Net(name = "ma_link", baseMVA = 100.0)
      for b in ("A1", "A2", "A3", "A3a", "A4")
        addBus!(net = net, busName = b, vn_kV = 110.0)
      end
      addPIModelACLine!(net = net, fromBus = "A1", toBus = "A2", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
      addPIModelACLine!(net = net, fromBus = "A2", toBus = "A3", r_pu = 0.012, x_pu = 0.09, b_pu = 0.0, status = 1)
      addPIModelACLine!(net = net, fromBus = "A3", toBus = "A4", r_pu = 0.011, x_pu = 0.085, b_pu = 0.0, status = 1)
      addPIModelACLine!(net = net, fromBus = "A4", toBus = "A1", r_pu = 0.01, x_pu = 0.07, b_pu = 0.0, status = 1)
      addProsumer!(net = net, busName = "A1", type = "EXTERNALNETWORKINJECTION", referencePri = "A1", vm_pu = 1.02, va_deg = 0.0)
      addProsumer!(net = net, busName = "A2", type = "ENERGYCONSUMER", p = 30.0, q = 10.0)
      addProsumer!(net = net, busName = "A3", type = "ENERGYCONSUMER", p = 12.0, q = 4.0)
      addProsumer!(net = net, busName = "A3a", type = "ENERGYCONSUMER", p = 3.0, q = 1.0)
      addProsumer!(net = net, busName = "A4", type = "ENERGYCONSUMER", p = 20.0, q = 6.0)
      addLink!(net = net, fromBus = "A3", toBus = "A3a", status = 1)
      _, erg = runpf!(net, 40, 1e-10, 0; method = :rectangular)
      @test erg == 0
      part = area_partition(net; areas = Dict("A1" => 1, "A2" => 1, "A3" => 2, "A3a" => 1, "A4" => 2))
      @test part.area_of_bus == [1, 1, 2, 2, 2]
      @test length(part.notes) == 1 && occursin("A3a", only(part.notes))
      Vpf = [nd._vm_pu * cis(deg2rad(nd._va_deg)) for nd in net.nodeVec]
      meas = generateMeasurementsFromPF(net)
      r = runse_multi_area!(deepcopy(net), meas; partition = part)
      @test r.converged
      @test maximum(abs.(r.voltages .- Vpf)) < 1e-5
    end)() end

    @testset "sp_case118: border information, power flow, runse!" begin (function ()
      net, _ = _ma_case118()
      Vpf = [nd._vm_pu * cis(deg2rad(nd._va_deg)) for nd in net.nodeVec]
      part = area_partition(net)
      cfg = Sparlectra.state_estimation_config()
      tol = cfg.tol
      # noise-free rows: the optimum of both estimates is the power flow
      # state (J = 0), an independent truth. Both runs stop once their step
      # is below tol, approaching from different sides: 10 tol per state.
      clean = generateMeasurementsFromPF(net)
      ma = Sparlectra._runse_two_level(net, clean, cfg, part)
      @test ma.converged
      @test maximum(abs.(abs.(ma.V) .- abs.(Vpf))) < 10 * tol
      @test maximum(abs.(angle.(ma.V ./ Vpf))) < 10 * tol
      # the sum of the area models is the Schur complement of the gain
      # matrix of the integrated estimate onto the border variables (rows
      # and interiors of the areas are disjoint), both at the power flow
      # state. The integrated gain comes from the estimator core, not from
      # the area code; the reference angle column of an area without the
      # slack is derived, so this checks it too.
      lin = Ref{Any}(nothing)
      Sparlectra._runse_with_config!(deepcopy(net), clean, cfg; linearization = lin)
      l = lin[]
      nb, s = l.nbus, l.slackIdx
      col(b, q) = q === :vm ? nb - 1 + b : (b < s ? b : b - 1)
      H = l.cap.H
      G = Matrix(H' * Diagonal(l.cap.w) * H)
      bcols = [col(b, q) for (b, q) in ma.vars]
      icols = setdiff(1:size(G, 1), bcols)
      Sint = G[bcols, bcols] - G[bcols, icols] * (G[icols, icols] \ G[icols, bcols])
      Asum, _ = Sparlectra._se_coordinator_system([ma.runs[a].model for a in part.areas], ma.vix, length(ma.vars))
      @test norm(Asum - Sint) / norm(Sint) < 1e-6
      # noisy rows: the fixed point is the estimate of runse! on the same
      # rows, J and dof included
      noisy = generateMeasurementsFromPF(net; noise = true, rng = MersenneTwister(1))
      ri = runse!(deepcopy(net), noisy)
      Vi = ri.voltages
      r = runse_multi_area!(deepcopy(net), noisy)
      @test r.converged && r.passes <= 6
      @test maximum(abs.(abs.(r.voltages) .- abs.(Vi))) < 10 * tol
      @test maximum(abs.(angle.(r.voltages ./ Vi))) < 10 * tol
      @test isapprox(r.objectiveJ, ri.objectiveJ; rtol = 1e-6) && r.dof == ri.dof
      @test r.comparison.max_dvm_pu < 10 * tol && r.comparison.integrated_converged
      @test all(x -> !isfinite(x.difference_sigma) || abs(x.difference_sigma) < 1e-2, r.comparison.rows)
      @test length(r.borders) == length(part.border_buses) && sum(a.rows for a in r.areas) == length(noisy)
      # PMU angles in one area without the slack (its offset is an interior
      # state of that area)
      sb = findfirst(nd -> getNodeType(nd) == Sparlectra.Slack, net.nodeVec)
      pa = first(a for a in part.areas if a != part.area_of_bus[sb])
      pmu = [b for b in eachindex(net.nodeVec) if part.area_of_bus[b] == pa][1:5:end]
      withPmu = generateMeasurementsFromPF(net; noise = true, rng = MersenneTwister(3), includeVa = true, vaBusIdxs = pmu, vaRefOffsetDeg = 7.0)
      rp = runse_multi_area!(deepcopy(net), withPmu; compare_integrated = false)
      rpi = runse!(deepcopy(net), withPmu)
      @test rp.converged && maximum(abs.(rp.voltages .- rpi.voltages)) < 10 * tol
      # refused: PMU rows in two areas, an area without rows
      both = generateMeasurementsFromPF(net; includeVa = true, vaBusIdxs = [b for b in eachindex(net.nodeVec) if !(b in part.border_buses)][1:20:end])
      @test_throws "more than one area" runse_multi_area!(deepcopy(net), both; compare_integrated = false)
      none1 = [m for m in clean if Sparlectra._se_area_owner(net, part, m) != 1]
      @test_throws "area 1 has no active measurement" runse_multi_area!(deepcopy(net), none1; compare_integrated = false)
      # bad data: three gross errors, one on the P flow of a tie branch.
      # The normalized residuals use the exact covariance of the whole
      # estimate, so the elimination sequence, the residuals and the stop
      # reason are those of runse_diagnostics on the same rows.
      bad = copy(noisy)
      kt = findfirst(m -> m.typ == Sparlectra.PflowMeas && m.branchIdx == first(part.tie_branches) && m.direction === :from, bad)
      kp = findfirst(m -> m.typ == Sparlectra.PinjMeas && m.busIdx == 30, bad)
      kq = findfirst(m -> m.typ == Sparlectra.QflowMeas && m.branchIdx == 50, bad)
      for (k, d) in ((kt, 40.0), (kp, -35.0), (kq, 25.0))
        m = bad[k]
        bad[k] = Sparlectra.Measurement(typ = m.typ, value = m.value + d * m.sigma, sigma = m.sigma, busIdx = m.busIdx, branchIdx = m.branchIdx, direction = m.direction, id = m.id)
      end
      di = runse_diagnostics(deepcopy(net), bad)
      dm = runse_multi_area_diagnostics(deepcopy(net), bad)
      @test [t.id for t in dm.trace] == [t.id for t in di.eliminations]
      @test bad[kt].id in [t.id for t in dm.trace]
      @test all(isapprox(a.normalized_residual_before, b.normalized_residual_before; rtol = 1e-4) for (a, b) in zip(dm.trace, di.eliminations))
      @test dm.stop_reason === di.stop_reason === :consistent
    end)() end

    @testset "tap regulators: inside an area estimated, on a tie branch not" begin (function ()
      # sp_case118 with stored ratio changers (1 % per step) on all nine
      # transformers, true steps between -3 and +3 (trafo 162, 8-28, at 0).
      # Bus 28 moved to area 2 makes 162 a tie branch (rule 1 of the
      # concept: its step would couple two areas, reason :tie_branch). The
      # reference is runse! with 162 excluded, the same model. The passes
      # keep the regulators continuous and fix them once at the converged
      # border: fixed every pass, nine regulators settled 1 to 5 steps off.
      stored(ratio, step) = Sparlectra.SETapChangerSet(ratio_n = ratio, changers = [Sparlectra.SETapChanger(name = "RTC", mrid = "_rtc", kind = :ratio, role = :ratio, side = :from, du = 0.01, neutral = 0, low = -8, high = 8, step = step)], origin = :cgmes)
      net, _ = _ma_case118()
      truth = Dict(162 => 0, 163 => 2, 164 => -3, 165 => 1, 166 => -2, 167 => 3, 168 => -1, 169 => 2, 170 => -2)
      for (k, st) in truth
        b = net.branchVec[k]
        b.se_tap_changers = stored(b.ratio, st)
        b.tap_ratio = b.ratio * (1 + st * 0.01)
      end
      runpf!(net, 60, 1e-10, 0)
      meas = generateMeasurementsFromPF(net; noise = true, rng = MersenneTwister(4))
      for k in keys(truth)
        b = net.branchVec[k]
        b.se_tap_changers = stored(b.ratio, 0)
        b.tap_ratio = b.ratio
        setTapEstimation!(net; trafo = k, mode = :ratio)
      end
      areas = Dict(getCompName(nd.comp) => (getCompName(nd.comp) == "28" ? 2 : Int(nd._area)) for nd in net.nodeVec)
      part = area_partition(net; areas = areas)
      @test 162 in part.tie_branches && all(k -> !(k in part.tie_branches), 163:170)
      r = with_logger(NullLogger()) do
        runse_multi_area!(deepcopy(net), meas; partition = part)
      end
      tie = getCompName(net.branchVec[162].comp)
      ri = with_logger(NullLogger()) do
        with_state_estimation_config(tap_estimation = Sparlectra.TapEstimationConfig(exclude = [tie])) do
          runse!(deepcopy(net), meas)
        end
      end
      tol = Sparlectra.state_estimation_config().tol
      @test r.converged && ri.converged
      @test maximum(abs.(r.voltages .- ri.voltages)) < 10 * tol
      @test r.comparison.integrated_converged && r.comparison.max_dvm_pu < 10 * tol
      step(rows, k) = only(t.step_fixed for t in rows if t.name == getCompName(net.branchVec[k].comp))
      @test all(step(r.tapEstimates, k) == step(ri.tapEstimates, k) == truth[k] for k in 163:170)
      @test all(t.reason === :tie_branch && t.status === :not_estimated for t in r.tapEstimates if t.name == tie)
      # the tie transformer sits in two subnets but has one row; the release
      # per area matches the integrated decision for all nine
      @test count(t -> t.name == tie, r.tapEstimates) == 1
      @test length(r.comparison.tap_release) == 9 && r.comparison.tap_release_differences == 0
      # a release that differs is a row with same = false: here the
      # integrated side kept 163 out
      maRows = [t for t in r.tapEstimates]
      intRows = [t.name == getCompName(net.branchVec[163].comp) ? merge(t, (status = :not_estimated, reason = :unobservable)) : t for t in maRows]
      rel = Sparlectra._se_tap_release_compare(maRows, intRows)
      @test count(t -> !t.same, rel) == 1 && only(t.branch for t in rel if !t.same) == 163
    end)() end

    @testset "sp_case10_areas, the shipped small case" begin (function ()
      # data/scf/sp_case10_areas (tools/build_demo_cases.jl): three areas
      # from the bus areas of the file, three tie lines, the coupler
      # Marbach_110 - Marbach_110b inside area 2, the released ratio changer
      # of Oberried_110_20 in area 3. The file holds step 0, the table comes
      # from the power flow at step +2: the step is known by construction.
      dir = joinpath(pkgdir(Sparlectra), "data", "scf")
      ref = Sparlectra.SparlectraConfig(Dict())
      net = Sparlectra.import_case(joinpath(dir, "sp_case10_areas.scf.json"), ref).net
      name(b) = getCompName(net.nodeVec[b].comp)
      part = area_partition(net)
      @test part.areas == [1, 2, 3]
      @test sort([getCompName(net.branchVec[k].comp) for k in part.tie_branches]) == ["Kirchau_Marbach", "Lindwerk_Oberried", "Nesselau_Pfaffing"]
      @test sort(name.(part.border_buses)) == sort(["Kirchau_110", "Marbach_110", "Lindwerk_110", "Oberried_110", "Nesselau_110", "Pfaffing_110"])
      @test part.area_of_bus[net.busDict["Marbach_110b"]] == 2
      trafo = only(k for (k, br) in enumerate(net.branchVec) if getCompName(br.comp) == "Oberried_110_20")
      @test net.branchVec[trafo].tap_est_mode === :ratio
      res = Sparlectra.with_sparlectra_config(ref) do
        readMeasurementTable!(net; file = joinpath(dir, "sp_case10_areas.measurements.csv"))
        meas = Measurement[m for m in net.measurements]
        with_logger(NullLogger()) do
          (ri = runse!(deepcopy(net), meas), rm = runse_multi_area!(deepcopy(net), meas))
        end
      end
      ri, rm = res.ri, res.rm
      tol = Sparlectra.state_estimation_config().tol
      @test ri.converged && rm.converged
      step(rows) = only(t.step_fixed for t in rows if t.name == "Oberried_110_20")
      @test step(ri.tapEstimates) == step(rm.tapEstimates) == 2
      # runse! reports the coupled busbars as one voltage, the comparison
      # inside the result holds the integrated estimate per bus
      @test rm.comparison.max_dvm_pu < 10 * tol && rm.comparison.max_dva_deg < rad2deg(10 * tol)
      @test rm.dof == ri.dof && isapprox(rm.objectiveJ, ri.objectiveJ; rtol = 1e-6)
      @test rm.comparison.tap_release_differences == 0
    end)() end

    @testset "areas of a CGMES import with two operator models" begin (function ()
      # data/cgmes_demo/sp_case6_igm: the import puts each bus in the area of
      # its operator model (cgmes_import.area_source = model_authority), the
      # two X nodes in Nord; the tie halves of Sued are the tie branches. A
      # noisy set from the power flow (seed 1) estimated in two areas equals
      # the integrated estimate on the same rows.
      net = importCGMES(path = joinpath(pkgdir(Sparlectra), "data", "cgmes_demo", "sp_case6_igm")).net
      @test runpf!(net, 30, 1e-10, 0; method = :rectangular)[2] == 0
      part = area_partition(net)
      @test part.areas == [1, 2]
      @test part.area_names == ["http://nord.example/sparlectra", "http://sued.example/sparlectra"]
      @test sort([getCompName(net.nodeVec[b].comp) for b in part.border_buses]) == ["Marbach", "Nesselau", "XN_Holm_Nes", "XN_Rabe_Mar"]
      meas = generateMeasurementsFromPF(net; noise = true, rng = MersenneTwister(1))
      rm = with_logger(NullLogger()) do
        runse_multi_area!(deepcopy(net), meas)
      end
      tol = Sparlectra.state_estimation_config().tol
      @test rm.converged && rm.comparison.integrated_converged
      @test [a.name for a in rm.areas] == part.area_names
      @test rm.comparison.max_dvm_pu < 10 * tol && rm.comparison.max_dva_deg < rad2deg(10 * tol)
    end)() end

    @testset "SE service with multi_area.mode = two_level" begin (function ()
      # the case file carries the bus areas (SCF extra.area) and the rows;
      # the service writes the area tables after the normal estimation, a
      # refusal of the method (zero_injection = constraint) is named and
      # the run stands
      net, _ = _ma_case118()
      net.measurements = generateMeasurementsFromPF(net; noise = true, rng = MersenneTwister(1))
      mktempdir() do d
        case = exportSCF(net; file = joinpath(d, "case.scf.json"), intended_calculations = String["power_flow", "state_estimation"])
        run(extra) = begin
          cfg = joinpath(d, "se_$(length(extra)).yaml")
          write(cfg, "config_version: 1\nstate_estimation:\n  multi_area:\n    mode: two_level\n" * extra)
          out = joinpath(d, "run_$(length(extra))")
          local res = with_logger(NullLogger()) do
            redirect_stdout(() -> SparlectraApp._run_state_estimation_service(case, cfg, out, "ma_$(length(extra))", joinpath(d, "none.csv")), devnull)
          end
          return SparlectraApp.to_dict(res), out
        end
        r, out = run("")
        md = r["metadata"]
        @test r["status"] == "succeeded" && md["se_multi_area"] === true && md["se_multi_area_converged"] === true
        @test md["se_multi_area_areas"] == 3 && md["se_multi_area_max_dvm_pu"] < 10 * Sparlectra.state_estimation_config().tol
        @test all(isfile(joinpath(out, f)) for f in ("se_areas.csv", "se_area_borders.csv", "se_area_compare.csv"))
        @test countlines(joinpath(out, "se_area_compare.csv")) == 1 + 2 * length(net.nodeVec) - 1
        @test occursin("multi-area (two_level): 3 area(s), 4 tie branch(es), converged", read(joinpath(out, "run.log"), String))
        r2, out2 = run("  zero_injection: constraint\n")
        @test r2["status"] == "succeeded" && r2["metadata"]["se_multi_area"] === false
        @test occursin("constraint", r2["metadata"]["se_multi_area_error"]) && !isfile(joinpath(out2, "se_areas.csv"))
        # the result page gets the Areas section for both runs, with the
        # area table when the method ran; none for a run without the mode
        sec = SparlectraApp._webui_se_multi_area_section(r)
        @test count("<tr><td>", sec) == 3 + 8   # three areas, eight border buses
        @test !isempty(SparlectraApp._webui_se_multi_area_section(r2))
        @test isempty(SparlectraApp._webui_se_multi_area_section(Dict{String,Any}("metadata" => Dict{String,Any}("run_mode" => "se"))))
      end
    end)() end
    run_multi_area_parallel_tests()
  end
end
