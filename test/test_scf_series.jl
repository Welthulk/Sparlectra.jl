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

# file: test/test_scf_series.jl
# purpose: the SCF time series (src/adapters/scf/scf_series.jl) on the shipped
#          data/qsts/sp_case14_day.qsts.zip: the check mode of
#          tools/build_scf_series.jl, the entries byte stable through write
#          and read, the load-time checks with the CSV line in the message,
#          parking of a status-0 unit, apply! and restore! bitwise on every
#          step, a power flow per step, the example controls of
#          sp_case14_day.control.yaml at the peak, off peak and in a
#          stressed step, format detection and the power flow entries
#          on the zip; QSTS mode A and mode B (sequential: the example day,
#          carry, restart, hold, resync, dead time, hunting, series in
#          parallel).

isdefined(@__MODULE__, :scf_roundtrip_field_diffs) || include("test_scf_support.jl")

const _SCF_CASE_DATA = normpath(joinpath(@__DIR__, "..", "data", "scf"))
const _SCF_SERIES_DATA = normpath(joinpath(@__DIR__, "..", "data", "qsts"))
const _SCF_SERIES_SHIPPED = joinpath(_SCF_SERIES_DATA, "sp_case14_day.qsts.zip")

# a two-step series zip on sp_case14 with `delta` as delta.csv, for the
# checks of the reader
function _scf_series_test_zip(dir::AbstractString, delta::AbstractString)::String
  base = read_scf_json(joinpath(_SCF_CASE_DATA, "sp_case14.scf.json"))
  steps = [SCFSeriesStep(1, "2026-06-15T00:00Z", "00:00"), SCFSeriesStep(2, "2026-06-15T00:15Z", "00:15")]
  entries = Sparlectra._scf_series_texts(SCFSeries(name = "t", base = base, steps = steps, ops = [PatchOp[], PatchOp[]], step_minutes = 15))
  entries[3] = "delta.csv" => String(delta)
  return Sparlectra._scf_series_write_zip(joinpath(dir, string("t", length(readdir(dir)), ".zip")), entries)
end

# two snapshot runs equal in every value, step and element (NaN equal to
# NaN): the identity a parallel run has to keep with the serial one
function _qsts_results_equal(a::QstsResult, b::QstsResult)::Bool
  same(x, y) = all(f -> isequal(getfield(x, f), getfield(y, f)), fieldnames(typeof(x)))
  length(a.steps) == length(b.steps) || return false
  return all(f -> f in (:base, :steps, :elements) || isequal(getfield(a, f), getfield(b, f)), fieldnames(QstsResult)) &&
         same(a.base, b.base) && same(a.elements, b.elements) && all(same(x, y) for (x, y) in zip(a.steps, b.steps))
end

# two values equal field by field through nested structs and arrays (NaN
# equal to NaN): the identity of a restarted or parallel sequential run
function _qsts_struct_equal(a, b)::Bool
  typeof(a) == typeof(b) || return false
  if a isa AbstractArray
    size(a) == size(b) || return false
    return all(_qsts_struct_equal(x, y) for (x, y) in zip(a, b))
  end
  (a isa Number || a isa AbstractString || a isa Symbol || a === nothing || a isa Tuple || a isa NamedTuple || isempty(fieldnames(typeof(a)))) && return isequal(a, b)
  return all(f -> _qsts_struct_equal(getfield(a, f), getfield(b, f)), fieldnames(typeof(a)))
end

# the example controls with qsts.mode set (and further qsts keys)
function _qsts_example_config(mode::Symbol; kw...)
  cfg = load_sparlectra_config(joinpath(_SCF_SERIES_DATA, "sp_case14_day.control.yaml"); reload = true)
  return Sparlectra._sparlectra_config_with(cfg; qsts = Sparlectra._struct_with(cfg.qsts; mode = mode, kw...))
end

# sp_case14 with every load scaled by factors[k] in step k
function _qsts_load_series(name::String, factors::Vector{Float64})::SCFSeries
  base = read_scf_json(joinpath(_SCF_CASE_DATA, "sp_case14.scf.json"))
  ops = [[PatchOp(op = :scale, target = :load, id = r.id, factor = f) for r in base.data.sym_load] for f in factors]
  return SCFSeries(name = name, base = base, steps = [SCFSeriesStep(k, "", string(k)) for k in eachindex(factors)], ops = ops)
end

# the value of one metric row, NaN without the row
function _qsts_metric(r::QstsResult, id::String, metric::String)::Float64
  m = qstsMetrics(r)
  i = findfirst(j -> m.element_id[j] == id && m.metric[j] == metric, eachindex(m.metric))
  return i === nothing ? NaN : m.value[i]
end

# the message of the ArgumentError `f` throws, "" when it throws none
function _scf_series_error_text(f)::String
  try
    f()
    return ""
  catch e
    e isa ArgumentError || rethrow()
    return sprint(showerror, e)
  end
end

# mode A on worker threads against the serial run `serial` of `series`
# (written to `serial_files`), five times per task count: identical values
# and identical bytes in all four files
function _qsts_mode_a_parallel_check(series::SCFSeries, cfg::SparlectraConfig, serial::QstsResult, serial_files::AbstractVector{<:AbstractString})
  dir = mktempdir()
  for tasks in unique([2, Threads.nthreads()]), rep in 1:5
    par = redirect_stdout(() -> runQsts(series; config = cfg, parallel_enabled = true, parallel_max_tasks = tasks, parallel_min_work_items = 2), devnull)
    @test _qsts_results_equal(serial, par)
    par_files = writeQstsResults(joinpath(dir, string("parallel_", tasks, "_", rep)), par)
    @test all(read(a) == read(b) for (a, b) in zip(serial_files, par_files))
  end
  parallel_comparison_report("QSTS mode A, 96 steps (5 repeats per task count)")
  return nothing
end

function _qsts_mode_b_parallel_testset(series::SCFSeries)
  @testset "mode B series in parallel" begin (function ()
    # three series of 24 steps from the example day, the series in
    # parallel and the steps in order; five runs per task count, identical
    # values and file bytes against the serial run
    local cfg = _qsts_example_config(:sequential)
    local part(r) = SCFSeries(name = string("day ", first(r)), base = series.base, steps = [SCFSeriesStep(i, series.steps[k].timestamp, series.steps[k].label) for (i, k) in enumerate(r)],
                              ops = series.ops[r], step_minutes = 15)
    local parts = [part(1:24), part(25:48), part(49:72)]
    local ticks = Int[]
    local serial = redirect_stdout(() -> runQsts(parts; config = cfg, parallel_enabled = false, progress = (done, total) -> push!(ticks, done)), devnull)
    @test length(serial) == 3 && all(r -> all(s -> s.status === :solved, r.steps), serial)
    @test ticks == collect(1:72)
    local dir = mktempdir()
    local serial_files = [writeQstsResults(joinpath(dir, string("serial_", s)), serial[s]) for s in 1:3]
    for tasks in unique([2, Threads.nthreads()]), rep in 1:5
      local par = redirect_stdout(() -> runQsts(parts; config = cfg, parallel_enabled = true, parallel_max_tasks = tasks, parallel_min_work_items = 2), devnull)
      @test all(_qsts_results_equal(serial[s], par[s]) for s in 1:3)
      local same_bytes = true
      for s in 1:3
        local files = writeQstsResults(joinpath(dir, string("parallel_", tasks, "_", rep, "_", s)), par[s])
        same_bytes = same_bytes && all(read(a) == read(b) for (a, b) in zip(serial_files[s], files))
      end
      @test same_bytes
    end
    parallel_comparison_report("QSTS mode B, three series of 24 steps (5 repeats per task count)")
  end)() end
end

function run_scf_series_tests()
  @testset "SCF time series" begin (function ()
    series = read_scf_series(_SCF_SERIES_SHIPPED)

    @testset "shipped series and the check mode of the build tool" begin (function ()
      @test (length(series.steps), sum(length, series.ops)) == (96, 968)
      @test series.index.origin === :scf_series && isempty(series.parked) && series.step_minutes == 15.0
      @test (series.steps[57].timestamp, series.steps[57].label) == ("2026-06-15T14:00Z", "14:00")
      # the base case is sp_case14 as shipped, byte for byte
      @test Sparlectra._scf_series_read_entries(_SCF_SERIES_SHIPPED)["base.scf.json"] == read(joinpath(_SCF_CASE_DATA, "sp_case14.scf.json"), String)
      @test [k for k in eachindex(series.ops) if any(op -> op.op === :status, series.ops[k])] == collect(57:64)
      local tool = Module(:BuildScfSeriesTool)
      Base.include(tool, normpath(joinpath(@__DIR__, "..", "tools", "build_scf_series.jl")))
      # the binding is newer than this function's world, so it is looked up
      # and called in the latest world
      local check = Base.invokelatest(getfield, tool, :check_shipped)
      @test Base.invokelatest(check) == String[]
    end)() end

    @testset "write and read keep the entries byte for byte" begin (function ()
      local dir = mktempdir()
      local shipped = Sparlectra._scf_series_read_entries(_SCF_SERIES_SHIPPED)
      local once = write_scf_series(joinpath(dir, "a.qsts.zip"), series)
      local twice = write_scf_series(joinpath(dir, "b.qsts.zip"), read_scf_series(once))
      @test Sparlectra._scf_series_read_entries(once) == shipped
      @test Sparlectra._scf_series_read_entries(twice) == shipped
      local r = Sparlectra.ZipArchives.ZipReader(read(once))
      @test Sparlectra.ZipArchives.zip_names(r) == ["series.json", "base.scf.json", "delta.csv"]
      @test all(i -> Sparlectra.ZipArchives.zip_iscompressed(r, i), 1:3)
      # a series without time reference, with an empty value, bus and
      # set ops: the reader returns the same ops, labels default to the step
      local base = read_scf_json(joinpath(_SCF_CASE_DATA, "sp_case14.scf.json"))
      local ops = [
        [PatchOp(op = :set, target = :load, id = 32, field = :p, value = 5.25), PatchOp(op = :set, target = :generator, id = 41, field = :qmin, value = nothing),
          PatchOp(op = :scale, target = :generator, id = 42, field = :p, factor = 0.5), PatchOp(op = :set, target = :bus, id = 1, field = :vm_start, value = 1.01)],
        [PatchOp(op = :status, target = :branch, id = 26, value = 0.0)],
      ]
      local plain = SCFSeries(name = "plain", base = base, steps = [SCFSeriesStep(1, "", "1"), SCFSeriesStep(2, "", "2")], ops = ops)
      local p1 = write_scf_series(joinpath(dir, "c.qsts.zip"), plain)
      local back = read_scf_series(p1)
      @test back.ops == ops && back.step_minutes === nothing && [s.label for s in back.steps] == ["1", "2"]
      @test Sparlectra._scf_series_read_entries(write_scf_series(joinpath(dir, "d.qsts.zip"), back)) == Sparlectra._scf_series_read_entries(p1)
      @test occursin("\n1,set,generator,41,qmin,\n", Sparlectra._scf_series_read_entries(p1)["delta.csv"])
    end)() end

    @testset "load-time checks name the CSV line" begin (function ()
      local dir = mktempdir()
      local head = "step,op,target,id,field,value\n"
      # (rows after the header, CSV line of the error, words of the message)
      local rows = [
        ("1,scale,load,32,0.9\n", 2, "expected 6 fields"),
        ("3,scale,load,32,,0.9\n", 2, "is not a step of series.json"),
        ("2,scale,load,32,,0.9\n1,scale,load,33,,0.9\n", 3, "step 1 after step 2"),
        ("1,drop,load,32,,1\n", 2, "op must be set, scale or status"),
        ("1,scale,tap_controller,46,,1\n", 2, "target \"tap_controller\" is not one of"),
        ("1,scale,load,32,,NaN\n", 2, "is not a finite number"),
        ("1,scale,load,99,,0.9\n", 2, "unknown component id 99"),
        ("1,scale,load,41,,0.9\n", 2, "id 41 is a generator, not a load"),
        ("1,scale,load,32,,0.9\n1,set,load,32,p,5\n", 3, "second write of p in step 1; the first is at delta.csv line 2"),
        ("1,set,load,32,p,\n", 2, "set p needs a value"),
        ("1,set,load,32,power,1\n", 2, "set field must be one of"),
        ("1,status,branch,26,p,0\n", 2, "status takes no field"),
        ("1,status,branch,26,,2\n", 2, "status value must be 0 or 1"),
        ("1,set,load,32,in_service,1\n", 2, "set in_service 1 needs the parked values of id 32"),
        ("1,set,generator,41,in_service,0\n1,set,generator,41,p,10\n", 3, "id 41 is parked in this step by set in_service 0 at delta.csv line 2"),
        ("1,scale,load,32,,0.9\n# late\n", 3, "a comment is allowed only above the header"),
      ]
      for (delta, line, words) in rows
        local msg = _scf_series_error_text(() -> read_scf_series(_scf_series_test_zip(dir, head * delta)))
        @test (delta, occursin(string("delta.csv line ", line, ": "), msg) && occursin(words, msg)) == (delta, true)
      end
      @test occursin("delta.csv line 1: the header must be", _scf_series_error_text(() -> read_scf_series(_scf_series_test_zip(dir, "step,op,target,id,value\n"))))
      # a byte order mark and a comment above the header are fine
      @test read_scf_series(_scf_series_test_zip(dir, "﻿# note\n" * head * "1,scale,load,32,,0.9\n")).ops[1][1].factor == 0.9
      # the entries and the manifest
      local texts = Sparlectra._scf_series_texts(series)
      local bad = entries -> _scf_series_error_text(() -> read_scf_series(Sparlectra._scf_series_write_zip(joinpath(dir, string("m", length(readdir(dir)), ".zip")), entries)))
      local manifest = texts[1][2]
      @test occursin("unexpected entry \"notes.txt\"", bad(vcat(texts, ["notes.txt" => "x"])))
      @test occursin("missing entry \"delta.csv\"", bad(texts[1:2]))
      @test occursin("format_version \"1.1\" is not supported", bad(["series.json" => replace(manifest, "\"1.0\"" => "\"1.1\""), texts[2], texts[3]]))
      @test occursin("series.json step 2: has the number 3", bad(["series.json" => replace(manifest, "\"step\": 2," => "\"step\": 3,"), texts[2], texts[3]]))
      @test occursin("series.json step 2: timestamp 2026-06-15T00:00Z is not after", bad(["series.json" => replace(manifest, "2026-06-15T00:15Z" => "2026-06-15T00:00Z"), texts[2], texts[3]]))
      @test occursin("base.scf.json: SCF JSON parse error", bad([texts[1], "base.scf.json" => "{", texts[3]]))
    end)() end

    @testset "a unit with status 0 in the base case is parked" begin (function ()
      local base = read_scf_json(joinpath(_SCF_CASE_DATA, "sp_case14.scf.json"))
      only(filter(r -> r.id == 33, base.data.sym_load)).status = 0
      local others = [PatchOp(op = :scale, target = :load, id = id, factor = 0.9) for id in (32, 34, 35, 36, 37, 38, 39)]
      local back_on = [PatchOp(op = :set, target = :load, id = 33, field = :in_service, value = 1.0), PatchOp(op = :scale, target = :load, id = 33, factor = 0.5)]
      local s = SCFSeries(name = "parked", base = base, steps = [SCFSeriesStep(k, "", string(k)) for k in 1:3], ops = [others, back_on, PatchOp[]])
      @test s.parked == [33] && s.index.parked_values[33].pVal == 14.0
      local plain = Sparlectra.build_net(base)
      local net = Sparlectra.build_net(s)
      local ref = Sparlectra.build_net(s)
      local k33 = s.index.internal_by_id[33]
      local bus = Sparlectra._scf_prosumer_bus(net.prosumpsVec[k33])
      # the plain reader keeps the unit in service; the series parks it in place
      @test (plain.prosumpsVec[k33].pVal, plain.nodeVec[bus]._pƩLoad) == (14.0, 14.0)
      @test (net.prosumpsVec[k33].pVal, net.prosumpsVec[k33].proSumptionType, net.nodeVec[bus]._pƩLoad) == (0.0, Sparlectra.Consumption, 0.0)
      @test length(net.prosumpsVec) == length(plain.prosumpsVec)
      # step 1 scales the others and leaves the parked unit parked
      local undo = apply!(net, s.ops[1], s.index)
      @test (net.prosumpsVec[k33].pVal, net.nodeVec[bus]._pƩLoad) == (0.0, 0.0)
      @test net.prosumpsVec[s.index.internal_by_id[32]].pVal == 18.0 * 0.9
      restore!(net, undo)
      @test _qsts_bitwise_diffs(ref, net) == String[]
      # step 2 brings it back with its base values, then scales it
      undo = apply!(net, s.ops[2], s.index)
      @test (net.prosumpsVec[k33].pVal, net.prosumpsVec[k33].qVal, net.nodeVec[bus]._pƩLoad) == (7.0, 2.0, 7.0)
      restore!(net, undo)
      @test _qsts_bitwise_diffs(ref, net) == String[]
      # status 0 stays the deletion of the N-1 vocabulary
      undo = apply!(net, [PatchOp(op = :status, target = :load, id = 34, value = 0.0)], s.index)
      @test length(net.prosumpsVec) == length(ref.prosumpsVec) - 1
      restore!(net, undo)
      @test _qsts_bitwise_diffs(ref, net) == String[]
      # the zip keeps the status, the reader parks again
      local back = read_scf_series(write_scf_series(joinpath(mktempdir(), "parked.qsts.zip"), s))
      @test back.parked == [33] && back.ops == s.ops
      # a scale or status 1 on the parked unit needs set in_service 1 in the
      # same step; status 1 alone would leave it parked
      local one_step = ops -> _scf_series_error_text(() -> SCFSeries(name = "x", base = base, steps = [SCFSeriesStep(1, "", "1")], ops = [ops]))
      @test occursin("step 1 op 1: id 33 is parked in the base case", one_step([PatchOp(op = :scale, target = :load, id = 33, factor = 0.5)]))
      @test occursin("step 1 op 1: id 33 is parked in the base case (status 0); status 1 does not bring it back", one_step([PatchOp(op = :status, target = :load, id = 33, value = 1.0)]))
      # a unit the step parks takes no other op: apply! parks first, so a
      # set p after it would make a load of the generator
      local park41 = PatchOp(op = :set, target = :generator, id = 41, field = :in_service, value = 0.0)
      @test one_step([park41]) == ""
      @test occursin("step 1 op 2: id 41 is parked in this step by set in_service 0 at step 1 op 1", one_step([park41, PatchOp(op = :scale, target = :generator, id = 41, factor = 0.5)]))
    end)() end

    @testset "every step: apply! and restore! return the base net bitwise" begin (function ()
      local net = Sparlectra.build_net(series)
      local ref = Sparlectra.build_net(series)
      local bad = Int[]
      for k in eachindex(series.ops)
        local undo = apply!(net, series.ops[k], series.index)
        restore!(net, undo)
        (isempty(_qsts_bitwise_diffs(ref, net)) && isempty(scf_roundtrip_field_diffs(ref, net))) || push!(bad, k)
      end
      @test bad == Int[]
    end)() end

    @testset "mode A: every step on its own, parallel equal to serial" begin (function ()
      local cfg = Sparlectra.SparlectraConfig(Dict{String,Any}("qsts" => Dict{String,Any}("mode" => "snapshots")))
      local serial = redirect_stdout(() -> runQsts(series; config = cfg, parallel_enabled = false), devnull)
      @test serial.base.status === :solved && all(s -> s.status === :solved, serial.steps)
      @test all(s -> s.control_status === :none && s.outer_iterations == 0, serial.steps)
      @test serial.origin === :scf_series && serial.base_step == 0 && serial.step_minutes == 15.0 && (serial.steps[57].label, serial.steps[57].timestamp) == ("14:00", "2026-06-15T14:00Z")
      # the line out from 14:00 to 15:45 has no rows in those steps
      local k26 = series.index.internal_by_id[26]
      @test [k for k in eachindex(serial.steps) if all(isnan, serial.branch_values[k26, :, k])] == collect(57:64)
      local serial_files = writeQstsResults(joinpath(mktempdir(), "serial"), serial)
      _qsts_mode_a_parallel_check(series, cfg, serial, serial_files)
      # step 57 equals a solve by hand bit for bit: a copy of the base net,
      # the voltages of the base solution, the delta, the step configuration
      local base_net = Sparlectra.build_net(series; config = cfg)
      local solved_base = deepcopy(base_net)
      @test redirect_stdout(() -> Sparlectra._execute_sparlectra_powerflow!(solved_base, Sparlectra._qsts_base_config(cfg)), devnull).erg == 0
      local work = deepcopy(base_net)
      Sparlectra._restore_start_voltages!(work, Sparlectra._snapshot_start_voltages(solved_base))
      apply!(work, series.ops[57], series.index)
      Sparlectra.markIsolatedBuses!(net = work, log = false, links = true)
      local ex = redirect_stdout(() -> Sparlectra._execute_sparlectra_powerflow!(work, first(Sparlectra._qsts_step_config(cfg, true))), devnull)
      calcNetLosses!(work)
      @test ex.erg == 0 && ex.iterations == serial.steps[57].iterations
      @test isequal(serial.bus_values[:, 1, 57], [nd._vm_pu for nd in work.nodeVec]) && isequal(serial.bus_values[:, 2, 57], [nd._va_deg for nd in work.nodeVec])
      local closed = [j for (j, br) in enumerate(work.branchVec) if br.status == 1]
      @test isequal(serial.branch_values[closed, 1, 57], [work.branchVec[j].fBranchFlow.pFlow for j in closed])
      @test isequal(serial.branch_values[closed, 4, 57], [work.branchVec[j].tBranchFlow.qFlow for j in closed])
      # the reported transformer flows are those of the solved model (live
      # tap), at both terminals
      local V = ComplexF64[nd._vm_pu * cis(deg2rad(nd._va_deg)) for nd in work.nodeVec]
      local trafo_gap = maximum(max(abs(complex(br.fBranchFlow.pFlow, br.fBranchFlow.qFlow) - Sparlectra._sv_branch_terminal_flow(br, V, :from, work.baseMVA)),
                                    abs(complex(br.tBranchFlow.pFlow, br.tBranchFlow.qFlow) - Sparlectra._sv_branch_terminal_flow(br, V, :to, work.baseMVA)))
                                for br in work.branchVec if br.ratio != 0.0)
      @test trafo_gap <= 1e-9
      # the long table, the CSV and the manifest
      local lt = qstsLongTable(serial)
      @test allequal(length.((lt.step, lt.element_id, lt.quantity, lt.value))) && length(lt.step) == countlines(serial_files[1]) - 1
      @test qstsValues(serial, "", "iterations") == Float64[s.iterations for s in serial.steps]
      @test isequal(qstsValues(serial, "26", "p_from_mw"), serial.branch_values[k26, 1, :])
      if Base.find_package("DataFrames") !== nothing
        local frames = Base.require(Base.identify_package("DataFrames"))
        local df = Base.invokelatest(getfield(frames, :DataFrame), lt)
        @test Base.invokelatest(size, df) == (length(lt.step), 4)
        println("      QSTS long table as DataFrame: RAN")
      else
        println("      QSTS long table as DataFrame: SKIPPED (DataFrames is not in the load path)")
      end
      local manifest = Sparlectra.scf_json_parse(read(serial_files[4], String))
      @test (manifest["format"], manifest["format_version"], manifest["n_steps"]) == ("sparlectra-qsts-result", "1.0", 96)
      @test run_csv(serial_files[1]).header == ["step", "element_id", "quantity", "value"]
      # a snapshot run writes four files and none of the mode B step
      # quantities, controller elements or files
      @test basename.(serial_files) == ["qsts_values.csv", "qsts_steps.csv", "qsts_elements.csv", "qsts_manifest.json"]
      @test isempty(serial.switching) && isempty(serial.metrics) && size(serial.controller_values) == (0, length(QSTS_CONTROLLER_QUANTITIES), 96)
      @test isempty(serial.elements.controller_ids) && serial.mode === :snapshots
      local mode_b_only = Set(vcat(QSTS_STEP_QUANTITIES[(findfirst(==("resync"), QSTS_STEP_QUANTITIES)):end], QSTS_CONTROLLER_QUANTITIES))
      local written = Set(row[3] for row in run_csv(serial_files[1]).rows)
      @test isempty(intersect(written, mode_b_only))
      @test !any(row -> "controller" in row, run_csv(serial_files[3]).rows)
      @test !haskey(manifest, "n_failed") && collect(keys(manifest["qsts"])) == ["snapshot_controls"]
    end)() end

    @testset "mode A: the tap controller only with qsts.snapshot_controls, a failed step" begin (function ()
      local base_net = Sparlectra.build_net(series)
      local moorau = Sparlectra._find_trafo_branch(base_net, only(c.trafo for c in Sparlectra._tap_controllers(base_net)))
      local j = findfirst(br -> br === moorau, base_net.branchVec)
      local tap = findfirst(==("tap_ratio"), QSTS_BRANCH_QUANTITIES)
      local plain = redirect_stdout(() -> runQsts(series; config = Sparlectra.SparlectraConfig(Dict{String,Any}("qsts" => Dict{String,Any}("mode" => "snapshots")))), devnull)
      local ctrl = redirect_stdout(() -> runQsts(series; config = Sparlectra.SparlectraConfig(Dict{String,Any}("qsts" => Dict{String,Any}("mode" => "snapshots", "snapshot_controls" => true)))), devnull)
      # without the key the tap stays at the value of the delta (the base
      # case, the series has no tap op); with it the controller moves it
      @test all(==(moorau.tap_ratio), plain.branch_values[j, tap, :])
      @test any(!=(moorau.tap_ratio), ctrl.branch_values[j, tap, :])
      @test all(s -> s.status === :solved && s.control_status === :converged, ctrl.steps) && any(s -> s.outer_iterations >= 1, ctrl.steps)
      # a tap_pos op on the regulated transformer is the tap of its step
      # without the controllers; with them it is refused, also on a
      # transformer whose controller comes from control.controllers
      local base = read_scf_json(joinpath(_SCF_CASE_DATA, "sp_case14.scf.json"))
      local tap_series(id) = SCFSeries(name = "taps", base = base, steps = [SCFSeriesStep(k, "", string(k)) for k in 1:2],
                                       ops = [PatchOp[], [PatchOp(op = :set, target = :transformer, id = id, field = :tap_pos, value = 5.0)]])
      local trafo_id(br) = only(id for (id, i) in series.index.internal_by_id if series.index.kind_by_id[id] === :transformer && base_net.branchVec[i] === br)
      local recorded = redirect_stdout(() -> runQsts(tap_series(trafo_id(moorau)); config = Sparlectra.SparlectraConfig(Dict{String,Any}("qsts" => Dict{String,Any}("mode" => "snapshots")))), devnull)
      @test all(s -> s.status === :solved, recorded.steps)
      @test recorded.branch_values[j, tap, :] == [moorau.tap_ratio, moorau.ratio / (1.0 + 5.0 * moorau.tap_step)]
      @test_throws "would fight it" runQsts(tap_series(trafo_id(moorau)); config = Sparlectra.SparlectraConfig(Dict{String,Any}("qsts" => Dict{String,Any}("mode" => "snapshots", "snapshot_controls" => true))))
      # a tap that starts outside its control range (case ratio 1.0, range
      # from 1.0125): the warning comes once per controller and run, not
      # once for the base case and once per step
      local range_cfg = Sparlectra.SparlectraConfig(Dict{String,Any}("qsts" => Dict{String,Any}("mode" => "snapshots", "snapshot_controls" => true),
        "control" => Dict{String,Any}("controllers" => Dict{String,Any}("moorau" => Dict{String,Any}("type" => "power_transformer", "trafo" => string(j), "mode" => "voltage", "tap_min" => 1.0125)))))
      local three = SCFSeries(name = "range", base = base, steps = [SCFSeriesStep(k, "", string(k)) for k in 1:3], ops = [PatchOp[], PatchOp[], PatchOp[]])
      local range_logs, _ = Test.collect_test_logs(() -> redirect_stdout(() -> runQsts(three; config = range_cfg), devnull))
      local starts = filter(l -> l.level == Logging.Warn && occursin("starts at ratio", string(l.message)), range_logs)
      @test length(starts) == 1
      @test occursin(r"transformer \S+ starts at ratio 1\.00000 \(tap position \+0\) below its control range \(tap_min = 1\.0125\)", string(only(starts).message))
      local wehrden = only(br for br in base_net.branchVec if br.ratio != 0.0 && br !== moorau)
      local wehrden_ctrl = Dict{String,Any}("type" => "power_transformer", "trafo" => wehrden.comp.cName, "mode" => "voltage", "target_bus" => base_net.nodeVec[wehrden.toBus].comp.cName, "target_vm_pu" => 1.0)
      @test_throws "would fight it" runQsts(tap_series(trafo_id(wehrden)); config = Sparlectra.SparlectraConfig(Dict{String,Any}("qsts" => Dict{String,Any}("mode" => "snapshots", "snapshot_controls" => true),
                                                                                                                  "control" => Dict{String,Any}("controllers" => Dict{String,Any}("wehrden" => wehrden_ctrl)))))
      # a step that does not converge keeps only its step counters, and the
      # run goes on
      local heavy_ops = [PatchOp(op = :scale, target = :load, id = r.id, factor = 50.0) for r in base.data.sym_load]
      local heavy = SCFSeries(name = "heavy", base = base, steps = [SCFSeriesStep(k, "", string(k)) for k in 1:3], ops = [PatchOp[], heavy_ops, PatchOp[]])
      local h = redirect_stdout(() -> runQsts(heavy; config = Sparlectra.SparlectraConfig(Dict{String,Any}("qsts" => Dict{String,Any}("mode" => "snapshots")))), devnull)
      @test [s.status for s in h.steps] == [:solved, :not_converged, :solved]
      local lt = qstsLongTable(h)
      local rows2 = findall(==(2), lt.step)
      @test sort(lt.quantity[rows2]) == sort(["converged", "iterations", "outer_iterations", "n_frozen"]) && all(isempty, lt.element_id[rows2])
      @test qstsValues(h, "", "converged") == [1.0, 0.0, 1.0]
      @test_throws "qsts.mode" runQsts(heavy; config = Sparlectra.SparlectraConfig(Dict()))
      # a step that converges with buses below power_flow.plausibility.vmin_warn
      # (every load x2.2: Torfloh_20 at 0.571 pu, three buses below 0.7 pu)
      # stays solved and gets n_low_vm, one note and one warning naming the
      # bus; the other steps get n_low_vm = 0, the threshold 0 no value
      local low_ops = [PatchOp(op = :scale, target = :load, id = r.id, factor = 2.2) for r in base.data.sym_load]
      local low = SCFSeries(name = "low", base = base, steps = [SCFSeriesStep(k, "", string(k)) for k in 1:3], ops = [PatchOp[], low_ops, PatchOp[]])
      local low_text = r"^QSTS step 2 \(2\): power flow plausibility: lowest voltage 0\.5712 pu at bus Torfloh_20 \(20\.0 kV\), 3 energized bus\(es\) below power_flow\.plausibility\.vmin_warn = 0\.7 pu"
      local lv = @test_logs (:warn, low_text) redirect_stdout(() -> runQsts(low; config = Sparlectra.SparlectraConfig(Dict{String,Any}("qsts" => Dict{String,Any}("mode" => "snapshots")))), devnull)
      @test all(s -> s.status === :solved, lv.steps)
      @test qstsValues(lv, "", "n_low_vm") == [0.0, 3.0, 0.0]
      @test count(t -> occursin(low_text, t), lv.notes) == 1 && !any(t -> startswith(t, "QSTS "), filter(t -> !occursin(low_text, t), lv.notes))
      local off = @test_logs redirect_stdout(() -> runQsts(low; config = Sparlectra.SparlectraConfig(Dict{String,Any}("qsts" => Dict{String,Any}("mode" => "snapshots"),
                                                                                                             "power_flow" => Dict{String,Any}("plausibility" => Dict{String,Any}("vmin_warn" => 0.0))))), devnull)
      @test all(isnan, qstsValues(off, "", "n_low_vm")) && !any(t -> startswith(t, "QSTS "), off.notes)
      # a bus behind a closed link (C at B) keeps its values while its link
      # group is fed, and loses them with it: the import types C Isolated
      # (branches only), and the step marking with the link groups must
      # win over that type
      local lcase = read_scf_json(exportSCF(link_group_net(; ring = false); file = joinpath(mktempdir(), "link_group.scf.json")))
      local lidx = ScenarioIndex(lcase)
      local ab = minimum(id for (id, k) in lidx.kind_by_id if k === :branch)
      local ls = SCFSeries(name = "link group", base = lcase, steps = [SCFSeriesStep(k, "", string(k)) for k in 1:2], ops = [PatchOp[], [PatchOp(op = :status, target = :branch, id = ab, value = 0.0)]])
      local lq = redirect_stdout(() -> runQsts(ls; config = Sparlectra.SparlectraConfig(Dict{String,Any}("qsts" => Dict{String,Any}("mode" => "snapshots")))), devnull)
      local lnet = Sparlectra.build_net(ls)
      local bB, bC = lnet.busDict["B"], lnet.busDict["C"]
      @test all(s -> s.status === :solved, lq.steps)
      @test isfinite(lq.bus_values[bC, 1, 1]) && isequal(lq.bus_values[bC, :, 1], lq.bus_values[bB, :, 1])
      @test all(isnan, lq.bus_values[[bB, bC], :, 2])
      # F2 in a snapshot step: lines 23, 28 and 29 open cut off Eichloh_110,
      # Wehrden_110 and Wehrden_20, fed by the Eichloh unit. With that unit
      # out in step 3 the island has no source and is isolated (no values,
      # its 8 MW load unsupplied); step 4 brings the unit back. Line 22 open
      # in step 5 cuts off three 20 kV buses without a source. Without the
      # rule steps 3 and 5 ended :solver_error and lost every value.
      local cut = [PatchOp(op = :status, target = :branch, id = id, value = 0.0) for id in (23, 28, 29)]
      local iops = [PatchOp[], cut, vcat(cut, [PatchOp(op = :set, target = :generator, id = 42, field = :in_service, value = 0.0)]), cut,
                    [PatchOp(op = :status, target = :branch, id = 22, value = 0.0)], PatchOp[]]
      local iseries = SCFSeries(name = "islands", base = base, steps = [SCFSeriesStep(k, "", string(k)) for k in eachindex(iops)], ops = iops)
      local iq = redirect_stdout(() -> runQsts(iseries; config = Sparlectra.SparlectraConfig(Dict{String,Any}("qsts" => Dict{String,Any}("mode" => "snapshots")))), devnull)
      @test all(s -> s.status === :solved && s.n_frozen == 0, iq.steps)
      local ivm(name) = iq.bus_values[findfirst(==(name), iq.elements.bus_names), 1, :]
      @test all(b -> isnan.(ivm(b)) == [false, false, true, false, false, false], ["Eichloh_110", "Wehrden_110", "Wehrden_20"])
      @test all(b -> isnan.(ivm(b)) == [false, false, false, false, true, false], ["Hassel_20", "Kleinau_20", "Torfloh_20"])
      @test isfinite(ivm("Moorau_20")[3]) && isfinite(ivm("Moorau_20")[5])
      @test isequal(iq.bus_values[:, :, 4], iq.bus_values[:, :, 2])
      @test filter(t -> startswith(t, "QSTS step"), iq.notes) == ["QSTS step 3 (3): 3 bus(es) without a source isolated, 8.0 MW load unsupplied",
                                                                  "QSTS step 5 (5): 3 bus(es) without a source isolated, 12.0 MW load unsupplied"]
      # the cut plus every load times 50: both islands miss convergence in
      # the island-wise solve, which throws for it. The step is
      # :not_converged with the island's iterations and a message naming
      # each island, not :solver_error with 0 iterations and the first line
      # (loads times 50 do not come near a solution: island 1 stops at a
      # mismatch of 76.5, the reason is :numerical_divergence)
      local hops = [PatchOp[], vcat(cut, heavy_ops), PatchOp[]]
      local hq = redirect_stdout(() -> runQsts(SCFSeries(name = "islands heavy", base = base, steps = [SCFSeriesStep(k, "", string(k)) for k in 1:3], ops = hops);
                                               config = Sparlectra.SparlectraConfig(Dict{String,Any}("qsts" => Dict{String,Any}("mode" => "snapshots"),
                                                                                                    "power_flow" => Dict{String,Any}("max_iter" => 8)))), devnull)
      @test [s.status for s in hq.steps] == [:solved, :not_converged, :solved]
      @test hq.steps[2].iterations == 8
      @test occursin(r"^AC island 1 power-flow solve failed: buses=\d+ .*; island 1 not_converged after 8 iterations, reason numerical_divergence, final mismatch [0-9.e+]+; island 2 not_converged after 8 iterations", hq.steps[2].message)
    end)() end

    @testset "the example controls attach and act in a step" begin (function ()
      # sp_case14_day.control.yaml adds the Wehrden tap controller (branch
      # 17), gives the case's tap controller of Moorau (branch 16) the
      # control range from 1.0 and puts the bank at Moorau_20 behind it at
      # rank 1 with target 1.01 and band 0.02, so the shared bus gives no
      # warning (the target lies on the edge of the tap band)
      local cfg = load_sparlectra_config(joinpath(_SCF_SERIES_DATA, "sp_case14_day.control.yaml"); reload = true)
      local net = Sparlectra.build_net(series)
      @test (@test_logs applyConfiguredControllers!(net, cfg.control)) == 2
      @test sort!([(c.trafo, c.target_bus, c.tap_min, c.priority) for c in Sparlectra._tap_controllers(net)]) == [("16", "Moorau_20", 1.0, 0), ("17", "Wehrden_20", nothing, 0)]
      local bank = only(Sparlectra._shunt_controllers(net))
      @test (bank.bus, bank.step_mvar, bank.bs_mvar, bank.bs_min_mvar, bank.bs_max_mvar, bank.target_vm_pu, bank.deadband_vm_pu, bank.priority) == ("Moorau_20", 2.5, 0.0, -5.0, 5.0, 1.01, 0.02, 1)
      local solve_step = function (ops)
        local n = deepcopy(net)
        apply!(n, ops, series.index)
        redirect_stdout(() -> run_sparlectra(net = n, config = cfg), devnull)
        local moorau = only(c for c in Sparlectra._tap_controllers(n) if c.trafo == "16")
        return n, moorau, only(Sparlectra._shunt_controllers(n))
      end
      # step 77 (19:00, evening peak) with the outer control loop: the
      # Moorau tap stands at its range limit 1.0 (without the range it goes
      # to 0.9875), the bank switched in one block after it, which brings
      # Moorau_20 back into the tap band, and the Wehrden tap left the case
      # value
      local n, moorau, nbank = solve_step(series.ops[77])
      @test (n.branchVec[16].tap_ratio, moorau.converged, nbank.bs_mvar, latest_control_result(n).status) == (1.0, true, 2.5, :converged)
      @test 0.99 <= get_bus_vm_pu(n, "Moorau_20") <= 1.01
      @test n.branchVec[17].tap_ratio != net.branchVec[17].tap_ratio
      # off peak (step 29, 07:00) the Moorau tap is in its band at 1.0 and
      # the bank stays off
      n, moorau, nbank = solve_step(series.ops[29])
      @test (n.branchVec[16].tap_ratio, moorau.converged, nbank.bs_mvar, latest_control_result(n).status) == (1.0, true, 0.0, :converged)
      # the peak step under 1.4 times its loads: the bank runs out of blocks
      # behind the tap at its limit, the voltage stays below the band of
      # both, and the closing solve leaves the voltages of a single solve
      local stressed = [op.op === :scale && op.target === :load ? Sparlectra.PatchOp(op.op, op.target, op.id, op.field, op.value, 1.4 * op.factor) : op for op in series.ops[77]]
      n, moorau, nbank = solve_step(stressed)
      @test (n.branchVec[16].tap_ratio, moorau.at_limit, nbank.bs_mvar, nbank.at_limit, latest_control_result(n).status) == (1.0, true, 5.0, true, :blocked)
      @test get_bus_vm_pu(n, "Moorau_20") < 0.99
      local single = deepcopy(n)
      redirect_stdout(() -> runpf!(single, 25, 1e-10, 0), devnull)
      @test maximum(abs(single.nodeVec[i]._vm_pu - n.nodeVec[i]._vm_pu) for i in eachindex(n.nodeVec)) < 1e-6
    end)() end

    @testset "format detection and the power flow entries run the base case" begin (function ()
      local plain_path = joinpath(_SCF_CASE_DATA, "sp_case14.scf.json")
      local smallgrid = joinpath(@__DIR__, "fixtures", "cgmes_smallgrid", "CGMES_v2.4.15_SmallGridTestConfiguration_BaseCase_Complete_v3.0.0.zip")
      @test is_scf_series(_SCF_SERIES_SHIPPED) && !is_scf_series(plain_path) && !is_scf_series(smallgrid)
      @test Sparlectra._detect_case_format(_SCF_SERIES_SHIPPED) === :scf_series
      @test Sparlectra._detect_case_format(smallgrid) === :cgmes
      # the content decides, not the name
      local dir = mktempdir()
      local renamed = joinpath(dir, "day.zip")
      cp(_SCF_SERIES_SHIPPED, renamed)
      @test Sparlectra._detect_case_format(renamed) === :scf_series
      local cfg = Sparlectra.load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH; reload = true)
      local imported = Sparlectra.import_case(renamed, cfg)
      @test imported.format === :scf_series
      @test scf_roundtrip_field_diffs(Sparlectra.build_net(series; config = cfg), imported.net) == String[]
      # the library entry gives the result of the plain base case
      local on_zip = run_sparlectra(casefile = renamed, config = cfg)
      local on_case = run_sparlectra(casefile = plain_path, config = cfg)
      @test on_zip.final_converged && on_zip.iterations == on_case.iterations
      @test [nd._vm_pu for nd in on_zip.net.nodeVec] == [nd._vm_pu for nd in on_case.net.nodeVec]
      # the service run succeeds and run.log names the series
      local rs = start_powerflow_run(Dict{String,Any}("casefile" => renamed, "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => joinpath(dir, "runs")))
      @test rs["status"] == "succeeded"
      @test occursin("SCF time series sp_case14_day: 96 steps", read(joinpath(rs["output_dir"], "run.log"), String))
      # a flat start is a flat start on the QSTS base case too (the base
      # solve skipped the run-level flat start and ran the DC start of the
      # projection before), in both modes: no projection line
      # from the solver, the flat-start line instead, the keys in the notes
      local pf = cfg.powerflow
      local flat = Sparlectra._copy_sparlectra_with_powerflow(cfg, Sparlectra._copy_powerflow_with(pf; start_mode = Sparlectra._copy_start_mode_with(pf.start_mode; flatstart = true, start_projection = true)))
      flat = Sparlectra._sparlectra_config_with(flat; output = Sparlectra._struct_with(flat.output; console_diagnostics = :full))
      local one = SCFSeries(name = "first", base = series.base, steps = series.steps[1:1], ops = series.ops[1:1], step_minutes = 15)
      for mode in (:snapshots, :sequential)
        local mcfg = Sparlectra._sparlectra_config_with(flat; qsts = Sparlectra._struct_with(flat.qsts; mode = mode))
        local logs, r = Test.collect_test_logs(() -> redirect_stdout(() -> mode === :sequential ? only(runQsts([one]; config = mcfg)) : runQsts(one; config = mcfg), devnull); min_level = Logging.Info)
        local msgs = [string(l.message) for l in logs]
        @test r.base.status === :solved
        @test !any(m -> startswith(m, "start projection selected"), msgs)
        @test any(n -> startswith(n, "flat start:") && occursin("power_flow.start_mode.start_projection=false", n), r.notes)
      end
    end)() end

    @testset "mode B on the example" begin (function ()
      # the example day with the controls of sp_case14_day.control.yaml in
      # order, taps and the bank carried from step to step
      local cfg = _qsts_example_config(:sequential)
      local r = redirect_stdout(() -> runQsts(series; config = cfg), devnull)
      @test r.mode === :sequential && r.base.status === :solved && all(s -> s.status === :solved, r.steps)
      @test r.elements.controller_ids == ["30#OLTC", "31#OLTC", "8#MSC"]
      @test r.elements.controller_names == ["16 OLTC", "17 OLTC", "msc_Moorau_20 MSC"]
      @test r.elements.controller_target == ["8", "14", "8"]
      # step 1 is a resync step from the base solution: bit for bit the step
      # of a snapshot run with the same controllers
      local one = SCFSeries(name = "first", base = series.base, steps = series.steps[1:1], ops = series.ops[1:1], step_minutes = 15)
      local a = redirect_stdout(() -> runQsts(one; config = _qsts_example_config(:snapshots; snapshot_controls = true)), devnull)
      @test isequal(a.bus_values[:, :, 1], r.bus_values[:, :, 1]) && isequal(a.branch_values[:, :, 1], r.branch_values[:, :, 1])
      local nA = findfirst(==("resync"), QSTS_STEP_QUANTITIES) - 1
      @test isequal(a.step_values[1:nA, 1], r.step_values[1:nA, 1]) && _qsts_struct_equal(a.steps[1], r.steps[1])
      local position(id) = qstsValues(r, id, "position")
      # the switching list is a chain: every row starts where the position
      # of the step before stood and ends at the position of its step, every
      # change of position has its row, and the switch count is the sum of
      # the counted steps
      local sw = qstsSwitching(r)
      @test !isempty(sw.step) && all(==("control"), sw.trigger)
      for id in r.elements.controller_ids
        local p = position(id)
        local rows = findall(==(id), sw.element_id)
        @test all(sw.old[i] == p[sw.step[i] - 1] && sw.new[i] == p[sw.step[i]] && sw.steps[i] == sw.new[i] - sw.old[i] for i in rows)
        @test [k for k in 2:96 if p[k] != p[k - 1]] == sw.step[rows]
        @test _qsts_metric(r, id, "switch_count") == sum(abs, sw.steps[rows]; init = 0)
      end
      @test (_qsts_metric(r, "", "switch_count_total"), _qsts_metric(r, "", "n_hunting_total"), _qsts_metric(r, "", "n_resync"), _qsts_metric(r, "", "n_failed")) == (sum(abs, sw.steps), 0.0, 1.0, 0.0)
      @test qstsValues(r, "", "n_moves") == Float64[count(==(k), sw.step) for k in 1:96]
      # the loss energy is the sum of the step losses over 15 minutes
      @test _qsts_metric(r, "", "loss_energy_mwh") ≈ sum(qstsValues(r, "", "total_p_loss_mw")) / 4
      # the files: two more than mode A, the manifest with the mode B keys
      local files = writeQstsResults(mktempdir(), r)
      @test basename.(files) == ["qsts_values.csv", "qsts_steps.csv", "qsts_elements.csv", "qsts_switching.csv", "qsts_metrics.csv", "qsts_manifest.json"]
      local switching = run_csv(files[4])
      @test switching.header == ["step", "element_id", "old", "new", "steps", "trigger"] && length(switching.rows) == length(sw.step)
      @test count(row -> "controller" in row, run_csv(files[3]).rows) == 3
      local manifest = Sparlectra.scf_json_parse(read(files[6], String))
      @test (manifest["mode"], manifest["n_failed"], manifest["qsts"]["sequential_start"], manifest["qsts"]["hunting_window"]) == ("sequential", 0, "previous", 3)
      # the long table holds the controller rows, the mode A checks apply
      @test count(==("30#OLTC"), qstsLongTable(r).element_id) == 96 * 6
      # the start from the snapshot solution of each step (run ST) gives the
      # same moves; the voltages differ within the solver tolerance
      local snap = redirect_stdout(() -> runQsts(series; config = _qsts_example_config(:sequential; sequential_start = :snapshot)), devnull)
      @test all(s -> s.status === :solved, snap.steps) && snap.switching == r.switching && isequal(snap.controller_values[:, 1, :], r.controller_values[:, 1, :])
      @test maximum(filter(!isnan, abs.(snap.bus_values[:, 1, :] .- r.bus_values[:, 1, :]))) < 1e-4
      local auto = Sparlectra._sparlectra_config_with(cfg; powerflow = Sparlectra._copy_powerflow_with(cfg.powerflow; mode = :auto))
      @test_throws "a sequential run needs power_flow.mode = manual" runQsts(series; config = auto)
    end)() end

    @testset "mode B carry" begin (function ()
      local cfg = _qsts_example_config(:sequential)
      local input = Sparlectra._qsts_scf_input(series, cfg)
      local opts = Sparlectra._QstsSeqOptions(false, nothing, nothing, nothing, true, false)
      local run_to(stop, initial = nothing) = redirect_stdout(() -> Sparlectra._run_qsts_series(input, cfg, opts, initial, stop), devnull)
      # a restart from step 40 with the carry of step 39 gives the whole run:
      # the carry holds everything one step hands to the next
      local whole = run_to(96)
      local to39 = run_to(39)
      local rest = run_to(96, to39.carry)
      @test _qsts_struct_equal(to39.outs, whole.outs[1:39])
      @test _qsts_struct_equal(rest.outs, whole.outs[40:96]) && _qsts_struct_equal(rest.carry, whole.carry)
      # a step that does not converge holds the carry: every load x50 in
      # step 2 between two steps at x1.2; the Moorau tap moved in step 1
      local heavy = _qsts_load_series("heavy", [1.2, 50.0, 1.2])
      local hcfg = Sparlectra.SparlectraConfig(Dict{String,Any}("qsts" => Dict{String,Any}("mode" => "sequential")))
      local hinput = Sparlectra._qsts_scf_input(heavy, hcfg)
      local hrun(stop) = redirect_stdout(() -> Sparlectra._run_qsts_series(hinput, hcfg, opts, nothing, stop), devnull)
      local h1 = hrun(1)
      local h3 = hrun(3)
      @test [o.out.info.status for o in h3.outs] == [:solved, :not_converged, :solved]
      local h2 = hrun(2)
      @test all(f -> _qsts_struct_equal(getfield(h1.carry, f), getfield(h2.carry, f)), (:taps, :banks, :positions, :vm, :va, :switch_count, :n_hunting, :last_converged))
      @test all(isnan, h3.outs[2].ctrl) && isempty(h3.outs[2].switching)
      # step 3 equals step 2 of the series without the failed step, apart from
      # the jumps against the step before, which a failed step leaves out
      local plain = _qsts_load_series("plain", [1.2, 1.2])
      local pinput = Sparlectra._qsts_scf_input(plain, hcfg)
      local p2 = redirect_stdout(() -> Sparlectra._run_qsts_series(pinput, hcfg, opts, nothing, 2), devnull)
      local o3 = h3.outs[3]
      local o2 = p2.outs[2]
      @test isequal(o3.out.bus, o2.out.bus) && isequal(o3.out.branch, o2.out.branch) && isequal(o3.ctrl, o2.ctrl)
      local jumps = [findfirst(==("max_dvm_step_pu"), QSTS_STEP_QUANTITIES), findfirst(==("max_dp_step_mw"), QSTS_STEP_QUANTITIES)]
      local others = setdiff(eachindex(QSTS_STEP_QUANTITIES), jumps)
      @test isequal(o3.out.step[others], o2.out.step[others]) && all(isnan, o3.out.step[jumps]) && all(isfinite, o2.out.step[jumps])
      local hr = redirect_stdout(() -> runQsts(heavy; config = hcfg), devnull)
      @test _qsts_metric(hr, "", "n_failed") == 1.0 && any(t -> occursin("1 step(s) did not converge", t), hr.notes)
      # resync every 8 steps: steps 1, 9, ..., 89 start from the positions of
      # the series, list their moves as resync and count none of them
      local rs = redirect_stdout(() -> runQsts(series; config = _qsts_example_config(:sequential; resync_every = 8)), devnull)
      local resync_steps = collect(1:8:96)
      @test findall(==(1.0), qstsValues(rs, "", "resync")) == resync_steps && _qsts_metric(rs, "", "n_resync") == 12.0
      local rsw = qstsSwitching(rs)
      @test all(rsw.trigger[i] == "resync" for i in eachindex(rsw.step) if rsw.step[i] in resync_steps)
      @test all(rsw.trigger[i] == "control" for i in eachindex(rsw.step) if !(rsw.step[i] in resync_steps))
      @test all(id -> all(==(0.0), qstsValues(rs, id, "moves")[resync_steps]), rs.elements.controller_ids)
      @test any(==("resync"), rsw.trigger)
      # dead time 2 on a load swing that moves the tap in every step: two
      # counted moves of one controller lie at least three steps apart, and
      # the controller is locked in the two steps after each
      local swing = _qsts_load_series("swing", [isodd(k) ? 1.0 : 1.2 for k in 1:12])
      local dcfg = Sparlectra.SparlectraConfig(Dict{String,Any}("qsts" => Dict{String,Any}("mode" => "sequential", "dead_time_steps" => 2)))
      local d = redirect_stdout(() -> runQsts(swing; config = dcfg), devnull)
      local dsw = qstsSwitching(d)
      local counted = dsw.step[dsw.trigger .== "control"]
      @test length(counted) >= 2 && all(>=(3), diff(counted))
      @test findall(==(1.0), qstsValues(d, "30#OLTC", "locked")) == filter(<=(12), sort!(vcat([[k + 1, k + 2] for k in counted]...)))
    end)() end

    @testset "mode B hunting" begin (function ()
      local h = Sparlectra._qsts_hunting
      @test h([(1, 1), (2, -1), (3, 1)], 3) == 1
      @test h([(1, 1), (2, -1), (3, 1), (4, -1)], 3) == 2
      @test (h([(1, 1), (2, -1), (4, 1)], 3), h([(1, 1), (2, -1), (4, 1)], 4)) == (0, 1)
      @test h([(1, 1), (2, 1), (3, -1)], 3) == 0
      @test h([(1, 1), (2, -1), (3, 0), (3, 1)], 3) == 0
      @test_throws "qsts.hunting_window must be >= 3" h([(1, 1)], 2)
      # the provoked pattern (run H0): loads at x1.0 and x1.2 in turn drive
      # the Moorau tap down and up in every step from step 2 on, seven moves
      # and five overlapping patterns; with one step of dead time none
      local swing = _qsts_load_series("swing", [isodd(k) ? 1.0 : 1.2 for k in 1:8])
      local run(dead) = redirect_stdout(() -> runQsts(swing; config = Sparlectra.SparlectraConfig(Dict{String,Any}("qsts" => Dict{String,Any}("mode" => "sequential", "dead_time_steps" => dead)))), devnull)
      local r = run(0)
      @test qstsValues(r, "30#OLTC", "moves") == [0.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0]
      @test _qsts_metric(r, "30#OLTC", "n_hunting") == 5.0 && _qsts_metric(r, "", "n_hunting_total") == 5.0
      @test _qsts_metric(run(1), "30#OLTC", "n_hunting") == 0.0
      # the counter of the run equals the pure function on its switching
      # list, with a reset at each resync step
      local sw = qstsSwitching(r)
      local moves = Tuple{Int,Int}[]
      for k in eachindex(r.steps)
        qstsValues(r, "", "resync")[k] == 1.0 && push!(moves, (k, 0))
        for i in eachindex(sw.step)
          (sw.step[i] == k && sw.trigger[i] == "control") && push!(moves, (k, sign(sw.steps[i])))
        end
      end
      @test h(moves, 3) == 5 && h(moves, 4) == _qsts_metric(redirect_stdout(() -> runQsts(swing; config = Sparlectra.SparlectraConfig(Dict{String,Any}("qsts" => Dict{String,Any}("mode" => "sequential", "hunting_window" => 4)))), devnull), "30#OLTC", "n_hunting")
    end)() end

    _qsts_mode_b_parallel_testset(series)

    @testset "event list" begin (function ()
      local dir = mktempdir()
      local head = "step,element_id,action,value\n"
      local events(rows; header = head) = (f = joinpath(dir, string("e", length(readdir(dir)), ".csv")); write(f, header * rows); f)
      local msg(f) = _scf_series_error_text(f)
      # the reader: every rule names the file line
      local bad = [
        ("2,26,open\n", 2, "expected 4 fields"),
        ("0,26,open,\n", 2, "step must be a whole number from 1"),
        ("x,26,open,\n", 2, "step must be a whole number from 1"),
        ("2,,open,\n", 2, "element_id is empty"),
        ("2,26,toggle,\n", 2, "action must be open, close, in or out"),
        ("2,26,open,1\n", 2, "value must be empty for action open"),
        ("2,26,open,\n3,22,open,\n2,26,close,\n", 4, "second event of 26 in step 2; the first is at line 2"),
      ]
      for (rows, line, words) in bad
        local m = msg(() -> readQstsEvents(events(rows)))
        @test (rows, occursin(string("line ", line, ": "), m) && occursin(words, m)) == (rows, true)
      end
      @test occursin("line 1: the header must be step,element_id,action,value", msg(() -> readQstsEvents(events("2,26,open,\n"; header = "step,element_mrid,action,value\n"))))
      local base = read_scf_json(joinpath(_SCF_CASE_DATA, "sp_case14.scf.json"))
      # the action must fit the class, and the element must exist
      local series_of(rows) = qstsEventSeries(base, readQstsEvents(events(rows)))
      @test occursin("line 2: action open does not apply to 32, a load", msg(() -> series_of("2,32,open,\n")))
      @test occursin("line 2: action out does not apply to 26, a branch", msg(() -> series_of("2,26,out,\n")))
      @test occursin("line 2: unknown element_id \"99\"", msg(() -> series_of("2,99,open,\n")))
      @test occursin("line 2: step 5 lies beyond n_steps = 3", msg(() -> qstsEventSeries(base, readQstsEvents(events("5,26,open,\n")); n_steps = 3)))
      # the deltas: a state holds until the next event of its element, a
      # unit goes out by parking, an event back to the base state drops the op
      local s = series_of("2,26,open,\n3,22,open,\n3,32,out,\n5,22,close,\n6,32,in,\n")
      @test length(s.steps) == 6 && isempty(s.ops[1]) && s.ops[2] == [PatchOp(op = :status, target = :branch, id = 26, value = 0.0)]
      @test s.ops[3] == [PatchOp(op = :status, target = :branch, id = 22, value = 0.0), PatchOp(op = :status, target = :branch, id = 26, value = 0.0),
                         PatchOp(op = :set, target = :load, id = 32, field = :in_service, value = 0.0)]
      @test s.ops[5] == [PatchOp(op = :status, target = :branch, id = 26, value = 0.0), PatchOp(op = :set, target = :load, id = 32, field = :in_service, value = 0.0)]
      @test s.ops[6] == s.ops[2] && s.provenance["event_list"]["rows"] == 5 && endswith(s.provenance["event_list"]["source"], ".csv")
      # the zip keeps the ops and the event record
      local back = read_scf_series(write_scf_series(joinpath(dir, "events.qsts.zip"), s))
      @test back.ops == s.ops && back.provenance["event_list"] == s.provenance["event_list"]
      # the run (E0): opening Moorau_Hassel in step 3 cuts off Hassel_20,
      # Kleinau_20 and Torfloh_20, which hold load and no source; without the
      # isolation that step ends :solver_error (no reference in the island).
      # The run has the example controllers on: in step 3 the Moorau tap
      # moves after the isolation, and the control pass solves again with
      # the cut-off buses still Isolated (before the bus type fix of 0.41.0
      # that solve saw them as islands without reference).
      local r = redirect_stdout(() -> runQsts(back; config = _qsts_example_config(:sequential)), devnull)
      @test all(st -> st.status === :solved, r.steps)
      @test qstsValues(r, "", "n_moves") == [0.0, 0.0, 1.0, 0.0, 1.0, 0.0]
      @test qstsValues(r, "30#OLTC", "position") == [0.0, 0.0, 1.0, 1.0, 0.0, 0.0]
      @test qstsValues(r, "", "n_islands") == [1.0, 1.0, 2.0, 2.0, 1.0, 1.0]
      @test qstsValues(r, "", "n_isolated") == [0.0, 0.0, 3.0, 3.0, 0.0, 0.0]
      # the load of the cut-off buses
      @test qstsValues(r, "", "shed_load_mw") == [0.0, 0.0, 12.0, 12.0, 0.0, 0.0]
      @test r.elements.bus_names[11:13] == ["Hassel_20", "Kleinau_20", "Torfloh_20"]
      @test all(b -> all(isnan, qstsValues(r, r.elements.bus_ids[b], "vm_pu")[3:4]) && all(isfinite, qstsValues(r, r.elements.bus_ids[b], "vm_pu")[[1, 2, 5, 6]]), 11:13)
      local manifest = Sparlectra.scf_json_parse(read(last(writeQstsResults(joinpath(dir, "run"), r)), String))
      @test manifest["event_list"]["rows"] == 5
      # a snapshot run and the snapshot start refuse an event series
      @test_throws "qsts.mode = sequential only" runQsts(back; config = Sparlectra.SparlectraConfig(Dict{String,Any}("qsts" => Dict{String,Any}("mode" => "snapshots"))))
      @test_throws "qsts.sequential_start = snapshot" runQsts(back; config = Sparlectra.SparlectraConfig(Dict{String,Any}("qsts" => Dict{String,Any}("mode" => "sequential", "sequential_start" => "snapshot"))))
    end)() end

    @testset "mode B chart" begin (function ()
      local r = redirect_stdout(() -> runQsts(series; config = _qsts_example_config(:sequential)), devnull)
      local solved = count(st -> st.status === :solved, r.steps)
      # Moorau_20 (bus 8): one point per solved step, the voltage band of
      # the run
      local c = qstsChart(r, "8", "vm_pu")
      @test count("<circle", c.svg) == solved == 96 && count("<polyline", c.svg) == 1
      @test c.band == (r.qsts.vm_min_pu, r.qsts.vm_max_pu)
      @test sprint(show, MIME"image/svg+xml"(), c) == c.svg && startswith(c.svg, "<svg ")
      # a branch flow has no band; the loading has its upper limit; a
      # controller and a step quantity chart as well
      @test qstsChart(r, "26", "p_from_mw").band === nothing
      @test isequal(qstsChart(r, "26", "loading_pct").band, (NaN, r.qsts.loading_max_pct))
      @test count("<circle", qstsChart(r, "8#MSC", "position").svg) == 96 && count("<circle", qstsChart(r, "", "total_p_loss_mw").svg) == 96
      # the line out from 14:00 to 15:45 breaks the flow line into two parts
      @test count("<polyline", qstsChart(r, "26", "p_from_mw").svg) == 2 && count("<circle", qstsChart(r, "26", "p_from_mw").svg) == 88
      @test qstsChart(r, "8", "vm_pu"; band = (0.95, NaN)).band === (0.95, NaN)
      @test qstsChart(r, "8", "vm_pu"; band = nothing).svg != c.svg
      @test qstsChart(r, "8", "vm_pu").svg == c.svg
      @test_throws "unknown element id" qstsChart(r, "no-such", "vm_pu")
      @test_throws "band must be :auto" qstsChart(r, "8", "vm_pu"; band = :wide)
    end)() end
  end)() end
end

# QSTS mode A and mode B on worker threads against the serial runs. The
# extd profile runs both inside run_scf_series_tests; the threads profile
# (tools/run_gates.sh threads) repeats them with four threads.
function run_scf_series_parallel_tests()
  @testset "QSTS on worker threads" begin (function ()
    local series = read_scf_series(_SCF_SERIES_SHIPPED)
    @testset "mode A: parallel steps equal serial" begin (function ()
      local cfg = Sparlectra.SparlectraConfig(Dict{String,Any}("qsts" => Dict{String,Any}("mode" => "snapshots")))
      local serial = redirect_stdout(() -> runQsts(series; config = cfg, parallel_enabled = false), devnull)
      @test all(s -> s.status === :solved, serial.steps)
      _qsts_mode_a_parallel_check(series, cfg, serial, writeQstsResults(joinpath(mktempdir(), "serial"), serial))
    end)() end
    _qsts_mode_b_parallel_testset(series)
  end)() end
end
