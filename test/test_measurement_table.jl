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

# file: test/test_measurement_table.jl
# purpose: tests the mRID-keyed measurement table on the ENTSO-E SmallGrid
#          (test/fixtures/cgmes_smallgrid): the format rules of the reader,
#          the import report (valid, inactive and problem rows per
#          snapshot), the group sums of injections, the sign conventions
#          (a noise-free table reproduces the SV state), and the bad-data
#          detection on the synthetic fixture
#          (test/fixtures/measurements/smallgrid_synthetic.table.csv,
#          built by tools/gen_smallgrid_measurements.jl), and the name
#          addressing of a case without CGMES, the generator from a
#          solved power flow on the shipped sp_case14, and the pseudo rows
#          for unobserved injections on both

using Test
using Sparlectra

const _SMALLGRID_DIR = joinpath(@__DIR__, "fixtures", "cgmes_smallgrid")
const _SMALLGRID_TABLE = joinpath(@__DIR__, "fixtures", "measurements", "smallgrid_synthetic.table.csv")

_smallgrid_import() = redirect_stdout(() -> importCGMES(path = sort!(filter(f -> endswith(f, ".zip"), readdir(_SMALLGRID_DIR; join = true)))), devnull)

# the id importMeasurementTable gives a terminal or node row
_table_row_id(r) = string(r.type, "@", isempty(r.terminal_mrid) ? r.node_mrid : r.terminal_mrid)

# a solved 110 kV triangle whose bus B2 carries two generators and two loads:
# their generated names (Gen_110, Ld_110.0) repeat at the bus. With
# `shunt_q` bus B3 also carries one shunt per entry (MVAr at 1 pu; with two
# entries their generated names repeat too).
function _shared_name_net(; shunt_q::Vector{Float64} = Float64[])
  net = Net(name = "shared_names", baseMVA = 100.0)
  for b in ("B1", "B2", "B3")
    addBus!(net = net, busName = b, vn_kV = 110.0)
  end
  addPIModelACLine!(net = net, fromBus = "B1", toBus = "B2", r_pu = 0.1, x_pu = 0.4, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "B2", toBus = "B3", r_pu = 0.1, x_pu = 0.4, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "B1", toBus = "B3", r_pu = 0.1, x_pu = 0.4, b_pu = 0.0, status = 1)
  addProsumer!(net = net, busName = "B1", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "B1")
  addProsumer!(net = net, busName = "B2", type = "SYNCHRONOUSMACHINE", p = 20.0, q = 5.0)
  addProsumer!(net = net, busName = "B2", type = "SYNCHRONOUSMACHINE", p = 15.0, q = 2.0)
  addProsumer!(net = net, busName = "B2", type = "ENERGYCONSUMER", p = 30.0, q = 10.0)
  addProsumer!(net = net, busName = "B2", type = "ENERGYCONSUMER", p = 25.0, q = 8.0)
  addProsumer!(net = net, busName = "B3", type = "ENERGYCONSUMER", p = 40.0, q = 12.0)
  for q in shunt_q
    addShunt!(net = net, busName = "B3", pShunt = 0.0, qShunt = q)
  end
  _, erg = runpf!(net, 40, 1e-10, 0)
  erg == 0 || error("_shared_name_net: power flow did not converge")
  return net
end

# write, read, import with the name index and estimate one table; returns
# the set and the SE result
function _estimate_table(net, rows)
  file = test_scratch_path(".table.csv")
  writeMeasurementTable(file, rows)
  set = only(importMeasurementTable(readMeasurementTable(file), nameMeasurementKeys(net)).sets)
  work = deepcopy(net)
  applyMeasurementSet!(work, set)
  se = with_state_estimation_config(update_net = false) do
    runse!(work)
  end
  return set, se
end

# The estimate of a v1 measurement vector `a` and of its converted table set
# `b` on `net`: the same J, state and dof, and the same eliminations of the
# bad-data pass (by the id map of the converter). Returns the number of
# eliminations.
function _converted_same_estimate(net, a, b, ids)
  wa = deepcopy(net)
  empty!(wa.measurements)
  append!(wa.measurements, a)
  wb = deepcopy(net)
  empty!(wb.measurements)
  append!(wb.measurements, b)
  ra, rb = (with_state_estimation_config(update_net = false) do
    runse!(w)
  end for w in (wa, wb))
  @test ra.objectiveJ === rb.objectiveJ && ra.voltages == rb.voltages && ra.dof == rb.dof
  da, db = (with_state_estimation_config(update_net = false, max_eliminations = 10) do
    runse_diagnostics(w, Measurement[m for m in w.measurements])
  end for w in (wa, wb))
  @test [get(ids, e.id, e.id) for e in da.eliminations] == [e.id for e in db.eliminations]
  @test da.final_diagnostics.objective.value === db.final_diagnostics.objective.value
  return length(da.eliminations)
end

function run_measurement_table_tests()
  @testset "measurement table: format rules" begin
    (function ()
      file = tempname() * ".csv"
      write(file, """
      # sparlectra-measurement-table v1
      # a comment above the header
      timestamp,measurement_mrid,terminal_mrid,node_mrid,type,value,sigma,source,valid
      2026-01-01T12:00Z,,t1,,P,1.5,0.1,scada,true
      2026-01-01T12:00Z,,t1,n1,P,1.5,0.1,scada,true
      2026-01-01T12:00Z,,,n1,I,1.5,0.1,scada,true
      2026-01-01T12:00Z,,t1,,P,abc,0.1,scada,true
      2026-01-01T12:00Z,,t1,,P,1.5,-1,scada,true
      2026-01-01T12:00Z,,t1,,P,1.5,0.1,scada,yes
      2020-01-01 12:00,,t1,,P,1.5,0.1,scada,true
      2026-01-01T12:00Z,,t1,,P,1.5
      2026-01-01T12:00Z,m1,,,P,1.5,,scada,true
      ,,,n1,U,110.0,,"op, manual",false
      """)
      @test isMeasurementTable(file)
      rows = readMeasurementTable(file)
      @test [r.problem for r in rows] == [:none, :key_count, :unknown_type, :bad_value, :bad_sigma, :bad_valid, :bad_timestamp, :field_count, :unknown_type, :none]
      @test rows[end].source == "op, manual" && rows[end].valid == false && rows[end].sigma === nothing
      # a file with another header is not a measurement table at all
      other = tempname() * ".csv"
      write(other, "type;bus;value\n")
      @test !isMeasurementTable(other)
      @test_throws ArgumentError readMeasurementTable(other)
      # the version line: missing or another version is refused by name;
      # the CSV v1 marker is not a table version
      body = join(filter(l -> !startswith(l, "# sparlectra-measurement-table"), readlines(file)), "\n")
      unversioned = tempname() * ".csv"
      write(unversioned, body)
      @test !isMeasurementTable(unversioned)
      @test_throws "expected '# sparlectra-measurement-table v1' before the header" readMeasurementTable(unversioned)
      v2 = tempname() * ".csv"
      write(v2, "# sparlectra-measurement-table v2\n" * body)
      @test isMeasurementTable(v2)
      @test_throws "version 'v2', this reader reads v1" readMeasurementTable(v2)
      csv1 = tempname() * ".csv"
      write(csv1, "# sparlectra-measurements v1\n" * body)
      @test !isMeasurementTable(csv1)
      @test_throws ArgumentError readMeasurementTable(csv1)
      # trailing blanks after the version are fine, other text is named
      tail = tempname() * ".csv"
      write(tail, "# sparlectra-measurement-table v1 generated by hand\n" * body)
      @test isMeasurementTable(tail)
      @test_throws "text after the version" readMeasurementTable(tail)
      write(tail, "# sparlectra-measurement-table v1 \t\n" * body)
      @test length(readMeasurementTable(tail)) == length(rows)
      # timestamps: one instant in several ISO 8601 forms is one key, an
      # impossible date or time is refused
      tsf = tempname() * ".csv"
      write(tsf, "# sparlectra-measurement-table v1\ntimestamp,measurement_mrid,terminal_mrid,node_mrid,type,value,sigma,source,valid\n" *
        join(("$(t),,,n1,U,110.0,,,true" for t in ("2026-01-01T12:00Z", "2026-01-01T12:00:00Z", "2026-01-01T12:00:00.000Z", "2026-01-01T12:00:30.500Z", "2024-02-29T00:00Z", "2026-02-30T10:00Z", "2026-01-01T25:61Z")), "\n") * "\n")
      tsr = readMeasurementTable(tsf)
      @test [r.problem for r in tsr] == [:none, :none, :none, :none, :none, :bad_timestamp, :bad_timestamp]
      @test [r.timestamp for r in tsr[1:5]] == ["2026-01-01T12:00Z", "2026-01-01T12:00Z", "2026-01-01T12:00Z", "2026-01-01T12:00:30.5Z", "2024-02-29T00:00Z"]
      # U and I are magnitudes: a negative value is refused, a negative P is
      # a direction
      negf = tempname() * ".csv"
      write(negf, "# sparlectra-measurement-table v1\ntimestamp,measurement_mrid,terminal_mrid,node_mrid,type,value,sigma,source,valid\n,,,n1,U,-110.0,,,true\n,,t1,,I,-5.0,,,true\n,,t1,,P,-5.0,,,true\n")
      @test [r.problem for r in readMeasurementTable(negf)] == [:bad_value, :bad_value, :none]
      # round trip through the writer, quoting included
      back = tempname() * ".csv"
      writeMeasurementTable(back, rows[[1, 10]])
      again = readMeasurementTable(back)
      @test [(r.terminal_mrid, r.node_mrid, r.value, r.source, r.valid) for r in again] == [(r.terminal_mrid, r.node_mrid, r.value, r.source, r.valid) for r in rows[[1, 10]]]
      @test readline(back) == "# sparlectra-measurement-table v1"
      # a comma and a quote survive the round trip; a line break cannot be
      # read back line by line, so the writer refuses it and writes nothing
      odd = MeasurementRow(1, "", "", "t1", "", "P", 1.0, 0.1, "meter \"A\", bay 2", true, :none, "")
      writeMeasurementTable(back, [odd])
      @test only(readMeasurementTable(back)).source == odd.source
      # blanks around a cell: the reader strips them, so the writer does too
      # and a second write gives the same bytes
      padded = MeasurementRow(1, "", "", " t1 ", "", "P", 1.0, 0.1, " meter A ", true, :none, "")
      f1, f2 = tempname() * ".csv", tempname() * ".csv"
      writeMeasurementTable(f1, [padded])
      writeMeasurementTable(f2, readMeasurementTable(f1))
      @test read(f1) == read(f2)
      broken = tempname() * ".csv"
      for r in (MeasurementRow(1, "", "", "t1", "", "P", 1.0, 0.1, "meter\nA", true, :none, ""), odd)
        @test_throws ArgumentError writeMeasurementTable(broken, [r]; comments = r === odd ? ["two\nlines"] : String[])
      end
      @test !isfile(broken)
      # the fixture reads and writes back byte for byte (version line,
      # comments, every number in its shortest round-trip form)
      fixture_comments = [l[3:end] for l in readlines(_SMALLGRID_TABLE) if startswith(l, "# ") && !startswith(l, "# sparlectra-measurement-table")]
      copy_file = tempname() * ".csv"
      writeMeasurementTable(copy_file, readMeasurementTable(_SMALLGRID_TABLE); comments = fixture_comments)
      @test read(copy_file) == read(_SMALLGRID_TABLE)
    end)()
  end

  @testset "measurement table: SmallGrid import report and group sums" begin
    (function ()
      result = _smallgrid_import()
      index = cgmesMeasurementKeys(result)
      rows = readMeasurementTable(_SMALLGRID_TABLE)
      imp = importMeasurementTable(rows, index)
      @test isempty(imp.rejected)
      @test [s.timestamp for s in imp.sets] == ["2026-01-01T12:00Z", "2026-01-01T13:00Z"]
      s1, s2 = imp.sets
      # snapshot 1: one row switched off, three handmade rows with one fault each
      @test s1.counts[:inactive] == 1
      @test (s1.counts[:unknown_key], s1.counts[:not_modelled], s1.counts[:unknown_type]) == (1, 1, 1)
      @test s1.counts[:rows] == s1.counts[:active] + s1.counts[:inactive] + 3
      @test sort!([p.problem for p in s1.problems]) == [:not_modelled, :unknown_key, :unknown_type]
      # snapshot 2: every injection row lands in a bus sum, nothing is lost
      @test length(s2.problems) == 0
      ninj = count(r -> r.timestamp == s2.timestamp && !isempty(r.terminal_mrid) && index.terminals[r.terminal_mrid].kind in (:injection, :shunt), rows)
      sums = filter(m -> m.typ in (Sparlectra.PinjMeas, Sparlectra.QinjMeas, Sparlectra.ShuntQMeas), s2.measurements)
      @test ninj > 0
      @test sum(length(s2.sources[m.id]) for m in sums) == ninj
      # a member whose only value is marked invalid makes the bus sum
      # inactive; a member without any row leaves it partial. Neither drops
      # a row silently.
      bus, members = first(sort!([(k, v) for (k, v) in index.groups if k[1] === :injection && length(v) >= 2]; by = x -> x[1][2]))
      victim = first(sort(members))
      flip(r) = MeasurementRow(r.line, r.timestamp, r.measurement_mrid, r.terminal_mrid, r.node_mrid, r.type, r.value, r.sigma, r.source, false, r.problem, r.detail)
      is_victim(r) = r.timestamp == s2.timestamp && r.terminal_mrid == victim && r.type == "P"
      s2inv = importMeasurementTable([is_victim(r) ? flip(r) : r for r in rows], index).sets[2]
      sum_p = only(m for m in s2inv.measurements if m.typ == Sparlectra.PinjMeas && m.busIdx == bus[2])
      @test !sum_p.active
      @test s2inv.counts[:invalid_rows] == 1 && s2inv.counts[:inactive] == 1
      s2cut = importMeasurementTable(filter(!is_victim, rows), index).sets[2]
      @test s2cut.counts[:partial_injection] == length(members) - 1
      @test !any(m -> m.typ == Sparlectra.PinjMeas && m.busIdx == bus[2], s2cut.measurements)
      # redundancy: U at a busbar section and at its TopologicalNode are two
      # measurements of one bus, both kept; the same key twice is a duplicate
      bbs, bbus = first(sort!([(k, b) for (k, b) in index.nodes if b > 0 && !any(r -> r.node_mrid == k, rows)]; by = first))
      tn_row = first(r for r in rows if r.timestamp == s2.timestamp && r.type == "U" && get(index.nodes, r.node_mrid, 0) == bbus)
      extra = [MeasurementRow(9001, s2.timestamp, "", "", bbs, "U", tn_row.value, tn_row.sigma, "test", true, :none, ""),
        MeasurementRow(9002, tn_row.timestamp, "", "", tn_row.node_mrid, "U", tn_row.value, tn_row.sigma, "test", true, :none, "")]
      s2red = importMeasurementTable(vcat(rows, extra), index).sets[2]
      @test count(m -> m.typ == Sparlectra.VmMeas && m.busIdx == bbus, s2red.measurements) == 2
      @test s2red.counts[:duplicate] == 1 && only(s2red.problems).line == 9002
      # the report names every problem row with its line
      report = measurementImportReport(imp)
      @test all(p -> occursin("line $(p.line): $(p.problem)", report), s1.problems)
    end)()
  end

  @testset "measurement table: U through a terminal follows its node" begin
    (function ()
      # an Analog U on a terminal reads the voltage of the terminal's node:
      # on a switch terminal at an energized node it is a VmMeas there, on a
      # terminal at a de-energized node it is not_modelled (the terminal
      # target keeps its real bus, as cgmesMeasurementKeys builds it)
      T = Sparlectra.MeasurementTarget
      index = Sparlectra.MeasurementKeyIndex(
        Dict("sw_live" => T(:not_modelled, 1, 0, :none, "Breaker is not a modelled equipment"), "ld_dead" => T(:not_modelled, 2, 0, :none, "EnergyConsumer at a node that is not energized")),
        Dict("tn_live" => 1, "tn_dead" => 0),
        Dict("a_live" => ("sw_live", "U", "LineToLineVoltage"), "a_dead" => ("ld_dead", "U", "LineToLineVoltage")),
        Dict{Tuple{Symbol,Int},Vector{String}}(), [110.0, 110.0])
      rows = [MeasurementRow(1, "", "a_live", "", "", "", 111.0, 0.5, "scada", true, :none, ""), MeasurementRow(2, "", "a_dead", "", "", "", 108.0, 0.5, "scada", true, :none, "")]
      set = only(importMeasurementTable(rows, index).sets)
      @test only(set.measurements).typ == Sparlectra.VmMeas && only(set.measurements).busIdx == 1
      @test only(set.problems).line == 2 && only(set.problems).problem === :not_modelled
      # the report names the terminal and the bus of its node, not the
      # terminal as a node
      @test occursin("terminal ld_dead", only(set.problems).detail) && occursin("bus 2", only(set.problems).detail)
      # a negative U through an Analog is refused like a negative U row
      neg_index = Sparlectra.MeasurementKeyIndex(index.terminals, index.nodes, merge(index.analogs, Dict("a_neg" => ("sw_live", "U", "LineToLineVoltage"))), index.groups, index.bus_kV)
      neg = only(importMeasurementTable([MeasurementRow(3, "", "a_neg", "", "", "", -111.0, 0.5, "scada", true, :none, "")], neg_index).sets)
      @test isempty(neg.measurements) && only(neg.problems).problem === :bad_value
    end)()
  end

  @testset "measurement table: invalid row next to a valid one of a member" begin
    (function ()
      # two Analogs on one load terminal, one marked invalid: the bus sum uses
      # the valid one and the invalid row is reported with its line
      T = Sparlectra.MeasurementTarget
      index = Sparlectra.MeasurementKeyIndex(
        Dict("ld1" => T(:injection, 1, 0, :none, "EnergyConsumer")), Dict("tn1" => 1),
        Dict("a1" => ("ld1", "P", "ThreePhaseActivePower"), "a2" => ("ld1", "P", "ThreePhaseActivePower")),
        Dict((:injection, 1) => ["ld1"]), [110.0])
      rows = [MeasurementRow(10, "", "a1", "", "", "", 5.0, 0.1, "a", true, :none, ""), MeasurementRow(11, "", "a2", "", "", "", 50.0, 0.1, "b", false, :none, "")]
      set = only(importMeasurementTable(rows, index).sets)
      m = only(set.measurements)
      @test m.typ == Sparlectra.PinjMeas && m.active && m.value == -5.0 && set.sources[m.id] == [10]
      @test only(set.problems).line == 11 && only(set.problems).problem === :invalid_unused
      @test occursin("line 11: invalid_unused", measurementImportReport(importMeasurementTable(rows, index)))
      # in an inactive sum (the second member has only an invalid row) the
      # first member still uses its valid row; its invalid row is reported
      index2 = Sparlectra.MeasurementKeyIndex(
        Dict("ld1" => T(:injection, 1, 0, :none, "EnergyConsumer"), "ld2" => T(:injection, 1, 0, :none, "EnergyConsumer")), Dict("tn1" => 1),
        Dict("a1" => ("ld1", "P", "ThreePhaseActivePower"), "a2" => ("ld1", "P", "ThreePhaseActivePower"), "b1" => ("ld2", "P", "ThreePhaseActivePower")),
        Dict((:injection, 1) => ["ld1", "ld2"]), [110.0])
      rows2 = [rows; MeasurementRow(12, "", "b1", "", "", "", 1.0, 0.1, "c", false, :none, "")]
      set2 = only(importMeasurementTable(rows2, index2).sets)
      m2 = only(set2.measurements)
      @test !m2.active && m2.value == -6.0 && set2.sources[m2.id] == [10, 12]
      @test only(set2.problems).line == 11 && only(set2.problems).problem === :invalid_unused
    end)()
  end

  @testset "measurement table: noise-free table reproduces the SV state" begin
    (function ()
      # sign and unit conventions: branch flows into the branch, U in kV;
      # an estimate on exact values has to land on the SV voltages
      result = _smallgrid_import()
      rows = generateMeasurementTableFromSV(result; noise = false)
      @test rows == generateMeasurementTableFromSV(result; noise = false)
      @test_throws ArgumentError generateMeasurementTableFromSV(result; sigma_u_rel = 0.0)
      @test_throws ArgumentError generateMeasurementTableFromSV(result; bad_data_sigmas = NaN)
      set = importMeasurementTable([MeasurementRow(i, r.timestamp, r.measurement_mrid, r.terminal_mrid, r.node_mrid, r.type, r.value, r.sigma, r.source, r.valid, :none, "") for (i, r) in enumerate(rows)], cgmesMeasurementKeys(result)).sets[1]
      @test length(set.measurements) == length(rows)
      # without zero-injection rows: the SV state of the delivery itself
      # leaves up to 0.14 MVAr at its passive buses, so ZI rows at 1 kW
      # would move the estimate off the SV voltages this test compares with
      applyMeasurementSet!(result.net, set; zero_injection = false)
      se = with_state_estimation_config(update_net = false) do
        runse!(result.net)
      end
      @test se.converged
      @test se.objectiveJ < 1e-8
      dvm = [abs(abs(se.voltages[m.busIdx]) - m.value) for m in set.measurements if m.typ == Sparlectra.VmMeas]
      @test maximum(dvm) < 1e-8
    end)()
  end

  @testset "measurement table: branch open at one end in a CGMES case" begin
    (function ()
      # runse! refuses flow measurements on a branch open at one terminal,
      # so the CGMES index and the SV generator leave both of its ends out
      # and name the reason; the estimate on the rest runs
      result = _smallgrid_import()
      net = result.net
      degree = zeros(Int, length(net.nodeVec))
      for br in net.branchVec
        br.status == 1 || continue
        degree[br.fromBus] += 1
        degree[br.toBus] += 1
      end
      k = first(k for (k, br) in enumerate(net.branchVec) if br.status == 1 && degree[br.fromBus] >= 3 && degree[br.toBus] >= 3)
      # the importer sets status 0 for such a branch; a net changed in code
      # can carry status 1 with one terminal flag 0 (_branch_terminal_state
      # reads that as open_to), and only that state reached the guard
      net.branchVec[k].to_status = 0
      terms = [t for (t, (bidx, _)) in result.branch_side_of_terminal if bidx == net.branchVec[k].branchIdx]
      @test length(terms) == 2
      index = cgmesMeasurementKeys(result)
      @test all(t -> index.terminals[Sparlectra._measurement_key(t)].kind === :not_modelled && occursin("partially open", index.terminals[Sparlectra._measurement_key(t)].detail), terms)
      rows = generateMeasurementTableFromSV(result; noise = false)
      @test !any(r -> r.terminal_mrid in Set(Sparlectra._measurement_key(t) for t in terms), rows)
      # a row on such a terminal is reported, not handed to the estimator
      table = [MeasurementRow(i, r.timestamp, r.measurement_mrid, r.terminal_mrid, r.node_mrid, r.type, r.value, r.sigma, r.source, r.valid, :none, "") for (i, r) in enumerate(rows)]
      push!(table, MeasurementRow(length(rows) + 1, "", "", Sparlectra._measurement_key(terms[1]), "", "P", 1.0, 0.1, "test", true, :none, ""))
      set = only(importMeasurementTable(table, index).sets)
      @test only(set.problems).problem === :not_modelled && only(set.problems).line == length(rows) + 1
      # no ZI rows, as in the SV test above (SV mismatch at passive buses)
      applyMeasurementSet!(net, set; zero_injection = false)
      se = with_state_estimation_config(update_net = false) do
        runse!(net)
      end
      @test se.converged
      @test se.objectiveJ < 1e-8
    end)()
  end

  @testset "measurement table: name addressing on sp_case14" begin
    (function ()
      net = load_fixture_net("sp_case14")
      _, erg = runpf!(net, 40, 1e-10, 0)
      @test erg == 0
      index = nameMeasurementKeys(net)
      busname = Dict(b => n for (n, b) in net.busDict)
      # the keys: bus names, <branch name>:from/:to, prosumer names; a name
      # several prosumers share is qualified with the bus
      @test all(haskey(index.nodes, n) for n in keys(net.busDict))
      # the file names the two parallel lines Ostheim_Neuwiese_1 and _2, so
      # every branch key is plain. A name two branches share (the generated
      # names of parallel MATPOWER lines repeat) is qualified with the
      # branch number, a unique name stays plain
      bkey(k) = (n = net.branchVec[k].comp.cName; count(b -> b.comp.cName == n, net.branchVec) > 1 ? string(n, "#", k) : n)
      @test all(k -> bkey(k) == net.branchVec[k].comp.cName, eachindex(net.branchVec))
      dup = deepcopy(net)
      dup.branchVec[2].comp.cName = dup.branchVec[1].comp.cName
      dindex = nameMeasurementKeys(dup)
      shared_name = dup.branchVec[1].comp.cName
      for (k, key) in ((1, shared_name * "#1"), (2, shared_name * "#2"), (3, dup.branchVec[3].comp.cName))
        t = dindex.terminals[key * ":to"]
        @test t.kind === :branch && t.branchIdx == k && t.direction === :to
      end
      @test !haskey(dindex.terminals, shared_name * ":to")
      shared = [ps for ps in net.prosumpsVec if count(q -> q.comp.cName == ps.comp.cName, net.prosumpsVec) > 1]
      @test !isempty(shared)
      ps = first(shared)
      t = index.terminals[ps.comp.cName * "@" * busname[Sparlectra.getPosumerBusIndex(ps)]]
      @test t.kind === :injection && t.busIdx == Sparlectra.getPosumerBusIndex(ps)
      @test isempty(index.analogs)
      # the bus-indexed set of the CSV generator, rewritten with names, gives
      # the same measurements through the table
      ref = generateMeasurementsFromPF(net; includePinj = false, includeQinj = false, noise = false)
      Row = NamedTuple{(:timestamp, :measurement_mrid, :terminal_mrid, :node_mrid, :type, :value, :sigma, :source, :valid),Tuple{String,String,String,String,String,Float64,Float64,String,Bool}}
      rows = Row[]
      for m in ref
        if m.typ == Sparlectra.VmMeas
          vn = net.nodeVec[m.busIdx].comp.cVN
          push!(rows, Row(("", "", "", busname[m.busIdx], "U", m.value * vn, m.sigma * vn, "test", true)))
        else
          push!(rows, Row(("", "", bkey(m.branchIdx) * ":" * string(m.direction), "", m.typ == Sparlectra.PflowMeas ? "P" : "Q", m.value, m.sigma, "test", true)))
        end
      end
      # rows a non-CGMES case cannot resolve: an unknown name and an mRID
      push!(rows, Row(("", "", "no_such_branch:from", "", "P", 1.0, 0.1, "test", true)))
      push!(rows, Row(("", "0447ee09-c766-11e1-8775-005056c00008", "", "", "", 1.0, 0.1, "test", true)))
      file = test_scratch_path(".table.csv")
      writeMeasurementTable(file, rows)
      imp = importMeasurementTable(readMeasurementTable(file), index)
      set = only(imp.sets)
      @test length(set.measurements) == length(ref)
      @test set.counts[:unknown_key] == 2
      report = measurementImportReport(imp)
      @test occursin("line $(length(rows) + 1): unknown_key", report) && occursin("line $(length(rows) + 2): unknown_key", report)
      # zero injection (state_estimation.zero_injection, default on): the
      # set has no injection rows, so every passive bus gets a ZI pair on
      # apply, counted in the set and the report; off adds none
      npass = length(findPassiveBuses(net))
      @test npass >= 1
      work = deepcopy(net)
      applyMeasurementSet!(work, set)
      @test set.counts[:zero_injection] == 2 * npass
      @test count(m -> startswith(m.id, "ZI_"), work.measurements) == 2 * npass
      @test length(work.measurements) == length(set.measurements) + 2 * npass
      @test occursin("zero-injection rows $(2 * npass)", measurementImportReport(imp))
      applyMeasurementSet!(work, set; zero_injection = false)
      @test set.counts[:zero_injection] == 0 && !any(m -> startswith(m.id, "ZI_"), work.measurements)
      with_state_estimation_config(zero_injection = :off) do
        applyMeasurementSet!(work, set)
      end
      @test set.counts[:zero_injection] == 0 && length(work.measurements) == length(set.measurements)
      key(m) = (m.typ, m.busIdx, m.branchIdx, m.direction)
      got = Dict(key(m) => m for m in set.measurements)
      @test all(isapprox(got[key(m)].value, m.value; atol = 1e-9) && isapprox(got[key(m)].sigma, m.sigma; atol = 1e-12) for m in ref)
      # the SE service picks the name index for a case without CGMES
      mktempdir() do d
        cfg = joinpath(d, "table_on.yaml")
        write(cfg, "config_version: 1\nstate_estimation:\n  measurement_table: true\n")
        case = abspath(joinpath(dirname(@__DIR__), "data", "scf", "sp_case14.scf.json"))
        res = redirect_stdout(devnull) do
          SparlectraApp._run_state_estimation_service(case, cfg, joinpath(d, "on"), "names", file)
        end
        dres = SparlectraApp.to_dict(res)
        @test dres["status"] == "succeeded"
        @test occursin("unknown_key 2", read(joinpath(d, "on", "measurement_import.txt"), String))
        # a coupler ring: two closed links from a PQ bus to a new busbar
        # with a 4 MW load. The balance fixes only the sum, so the service
        # splits it evenly without link rows (2/2); with the link P rows the
        # estimated ring circulation gives 3/1 (consistent rows, sigma 0.01)
        rnet = deepcopy(net)
        rnet._locked = false
        pq = net.nodeVec[findfirst(n -> getNodeType(n) == Sparlectra.PQ, net.nodeVec)].comp.cName
        addBus!(net = rnet, busName = "RB", vn_kV = rnet.nodeVec[rnet.busDict[pq]].comp.cVN)
        addProsumer!(net = rnet, busName = "RB", type = "ENERGYCONSUMER", p = 4.0, q = 1.0)
        addLink!(net = rnet, fromBus = pq, toBus = "RB")
        addLink!(net = rnet, fromBus = pq, toBus = "RB")
        _, erg_r = runpf!(rnet, 40, 1e-10, 0)
        @test erg_r == 0
        rcase = joinpath(d, "ring.scf.json")
        exportSCF(rnet; file = rcase)
        l1, l2 = (l.cName for l in rnet.linkVec)
        # U at every bus and every branch flow of the solved ring case
        rname = Dict(b => n for (n, b) in rnet.busDict)
        ring_rows = Row[]
        for m in generateMeasurementsFromPF(rnet; includePinj = false, includeQinj = false, noise = false)
          if m.typ == Sparlectra.VmMeas
            vn = rnet.nodeVec[m.busIdx].comp.cVN
            push!(ring_rows, Row(("", "", "", rname[m.busIdx], "U", m.value * vn, m.sigma * vn, "test", true)))
          elseif m.branchIdx !== nothing
            push!(ring_rows, Row(("", "", bkey(m.branchIdx) * ":" * string(m.direction), "", m.typ == Sparlectra.PflowMeas ? "P" : "Q", m.value, m.sigma, "test", true)))
          end
        end
        rfile_kcl = joinpath(d, "ring_kcl.table.csv")
        rfile_meas = joinpath(d, "ring_meas.table.csv")
        writeMeasurementTable(rfile_kcl, ring_rows)
        writeMeasurementTable(rfile_meas, vcat(ring_rows, [Row(("", "", l1 * ":from", "", "P", 3.0, 0.01, "test", true)), Row(("", "", l2 * ":from", "", "P", 1.0, 0.01, "test", true))]))
        link_p(tag, f) = begin
          r = redirect_stdout(devnull) do
            SparlectraApp._run_state_estimation_service(rcase, cfg, joinpath(d, tag), tag, f)
          end
          @test SparlectraApp.to_dict(r)["status"] == "succeeded"
          @test occursin("link flows: se_link_flows.csv (2 closed link(s)", read(joinpath(d, tag, "run.log"), String))
          L = readlines(joinpath(d, tag, "se_link_flows.csv"))
          dl = occursin(';', L[1]) ? ';' : ','
          col = findfirst(==("p_mw"), split(L[1], dl))
          [parse(Float64, replace(split(x, dl)[col], ',' => '.')) for x in L[2:end]]
        end
        pk = link_p("ring_kcl", rfile_kcl)
        @test all(isapprox.(pk, 2.0; atol = 0.01))
        pm = link_p("ring_meas", rfile_meas)
        @test isapprox(pm[1], 3.0; atol = 0.01) && isapprox(pm[2], 1.0; atol = 0.01)
      end
    end)()
  end

  @testset "measurement table: generator from a solved power flow" begin
    (function ()
      net = load_fixture_net("sp_case14")
      _, erg = runpf!(net, 40, 1e-10, 0)
      @test erg == 0
      Vpf = Sparlectra.buildVoltageVector(net)
      index = nameMeasurementKeys(net)
      estimate(rows) = begin
        file = test_scratch_path(".table.csv")
        writeMeasurementTable(file, rows)
        set = only(importMeasurementTable(readMeasurementTable(file), index).sets)
        work = deepcopy(net)
        applyMeasurementSet!(work, set)
        se = with_state_estimation_config(update_net = false) do
          runse!(work)
        end
        (set, se)
      end
      # noise 0 reproduces the solved state (acceptance 1e-6 pu), with and
      # without the injection rows (their sign: power into the equipment)
      for inj in (false, true)
        set, se = estimate(generateMeasurementTableFromPF(net; noise = false, include_injections = inj))
        @test isempty(set.problems)
        @test se.converged
        @test se.objectiveJ < 1e-8
        @test maximum(abs.(se.voltages .- Vpf)) < 1e-6
        @test inj == any(m -> m.typ == Sparlectra.PinjMeas, set.measurements)
      end
      # a sigma whose weight 1/sigma^2 overflows is lifted to the floor of the
      # CSV v1 reader and counted; the estimate keeps a finite J
      exact = generateMeasurementTableFromPF(net; noise = false)
      tiny_u = findfirst(r -> r.type == "U", exact)
      tiny_p = findfirst(r -> r.type == "P", exact)
      set, se = estimate([i in (tiny_u, tiny_p) ? merge(r, (sigma = 1e-200,)) : r for (i, r) in enumerate(exact)])
      @test set.counts[:sigma_lifted] == 2
      @test count(m -> m.sigma == Sparlectra.ZERO_INJECTION_SIGMA, set.measurements) == 2
      @test se.converged && isfinite(se.objectiveJ)
      # the floor follows state_estimation.zero_injection_sigma_mva below
      # ZERO_INJECTION_SIGMA (a zero-injection row written at 1e-7 keeps
      # it) and never rises above it (30 would lift every real row)
      tiny_file = test_scratch_path(".table.csv")
      writeMeasurementTable(tiny_file, [i in (tiny_u, tiny_p) ? merge(r, (sigma = 1e-200,)) : r for (i, r) in enumerate(exact)])
      for (zs, floor) in ((1e-7, 1e-7), (30.0, Sparlectra.ZERO_INJECTION_SIGMA))
        tiny_set = with_state_estimation_config(zero_injection_sigma_mva = zs) do
          only(importMeasurementTable(readMeasurementTable(tiny_file), index).sets)
        end
        @test tiny_set.counts[:sigma_lifted] == 2
        @test count(m -> m.sigma == floor, tiny_set.measurements) == 2
      end
      # the keyword is the sigma of the caller's configuration (the SE
      # service passes its case value) and wins over the active default
      kw_set = only(importMeasurementTable(readMeasurementTable(tiny_file), index; zero_injection_sigma = 1e-7).sets)
      @test count(m -> m.sigma == 1e-7, kw_set.measurements) == 2
      # coverage: the share of buses and branch ends, drawn by the seed
      full = generateMeasurementTableFromPF(net; noise = false)
      half = generateMeasurementTableFromPF(net; noise = false, coverage = 0.5, seed = 3)
      nU(rows) = count(r -> r.type == "U", rows)
      nP(rows) = count(r -> r.type == "P", rows)
      @test nU(full) == length(net.nodeVec) && nP(full) == 2 * length(net.branchVec)
      @test nU(half) == round(Int, 0.5 * nU(full)) && nP(half) == round(Int, 0.5 * nP(full))
      @test half == generateMeasurementTableFromPF(net; noise = false, coverage = 0.5, seed = 3)
      @test half != generateMeasurementTableFromPF(net; noise = false, coverage = 0.5, seed = 4)
      # noise and gross errors follow the seed as well
      noisy = generateMeasurementTableFromPF(net; bad_data_fraction = 0.05, seed = 7)
      @test noisy == generateMeasurementTableFromPF(net; bad_data_fraction = 0.05, seed = 7)
      @test count(r -> r.source == "synthetic;bad_data", noisy) == round(Int, 0.05 * length(noisy))
      @test_throws ArgumentError generateMeasurementTableFromPF(net; coverage = 1.5)
      # a sigma the reader would reject (0 floor, 0 relative U sigma) is refused
      for kw in ((sigma_u_rel = 0.0,), (sigma_pq_floor = 0.0,), (sigma_pq_rel = -0.01,), (bad_data_sigmas = NaN, bad_data_fraction = 0.5), (bad_data_sigmas = -1.0,))
        @test_throws ArgumentError generateMeasurementTableFromPF(net; kw...)
      end
    end)()
  end

  @testset "measurement table: sigma by voltage level" begin
    (function ()
      L(; kw...) = Sparlectra.SESigmaLevel(; kw...)
      # hand calculation, p_mw = 2 and rel_pct = 1: max(2, 1 % of 300 MW) =
      # 3 MW at a reading of 300 MW, max(2, 1 % of 100 MW) = 2 MW at 100 MW
      lv = [L(vn_kv = 110, p_mw = 2.0, rel_pct = 1.0)]
      @test Sparlectra._sigma_from_level(lv, 110.0, "P", 300.0) == 3.0
      @test Sparlectra._sigma_from_level(lv, 110.0, "P", -100.0) == 2.0
      # a column the row does not set and a level without a row: default
      @test Sparlectra._sigma_from_level(lv, 110.0, "Q", 300.0) === nothing
      @test Sparlectra._sigma_from_level(lv, 20.0, "P", 300.0) === nothing

      # sp_case14 has two levels, 110 kV and 20 kV; transformer
      # Moorau_110_20 runs from Moorau_110 to Moorau_20
      net = load_fixture_net("sp_case14")
      _, erg = runpf!(net, 40, 1e-10, 0)
      @test erg == 0
      Vpf = Sparlectra.buildVoltageVector(net)
      index = nameMeasurementKeys(net)
      two = [L(vn_kv = 110, p_mw = 1.5, q_mvar = 2.5, u_kv = 0.4, rel_pct = 0.5), L(vn_kv = 20, p_mw = 0.3, q_mvar = 0.4, u_kv = 0.08)]
      # hand rows: (key, node?, type, value, own sigma)
      hand = [("Moorau_110", true, "U", 115.5, nothing), ("Moorau_20", true, "U", 20.5, nothing),
        ("Moorau_110_20:from", false, "P", 400.0, nothing), ("Moorau_110_20:to", false, "P", -40.0, nothing),
        ("Neuwiese_Moorau:from", false, "Q", 100.0, nothing), ("Moorau_Hassel:from", false, "Q", 5.0, nothing),
        ("Moorau_Dallwitz:from", false, "P", 50.0, 0.7), ("Ld_20.0@Hassel_20", false, "P", -2.0, nothing)]
      hrows = [MeasurementRow(i, "", "", node ? "" : key, node ? key : "", typ, v, s, "test", true, :none, "") for (i, (key, node, typ, v, s)) in enumerate(hand)]
      hset(levels) = (s = only(importMeasurementTable(hrows, index; sigma_by_voltage_level = levels).sets); @test isempty(s.problems); (Dict(m.id => m.sigma for m in s.measurements), s.sigma_origin))
      h, ho = hset(two)
      # U in pu of its bus: max(0.4, 0.5 % of 115.5 kV) = 0.5775 kV over
      # 110 kV, and 0.08 kV over 20 kV
      @test h["U@Moorau_110"] ≈ 0.00525 rtol = 1e-14
      @test h["U@Moorau_20"] ≈ 0.004 rtol = 1e-14
      # P and Q at the bus of the branch end: the 110 kV end of the
      # transformer max(1.5, 0.5 % of 400 MW) = 2 MW, its 20 kV end 0.3 MW
      @test h["P@Moorau_110_20:from"] == 2.0 && h["P@Moorau_110_20:to"] == 0.3
      @test h["Q@Neuwiese_Moorau:from"] == 2.5 && h["Q@Moorau_Hassel:from"] == 0.4
      # the row's own sigma wins
      @test h["P@Moorau_Dallwitz:from"] == 0.7 && ho["P@Moorau_Dallwitz:from"] === :row
      # a bus sum of one member carries the member sigma
      @test h["PinjMeas@bus11"] == 0.3 && ho["PinjMeas@bus11"] === :level && ho["U@Moorau_20"] === :level
      # the bus sum is the root sum of its members: four injections at B2
      # (110 kV), 600 MW gives 3 MW, 100 MW three times 1.5 MW, so
      # sqrt(3^2 + 3 * 1.5^2) = sqrt(15.75)
      shared = _shared_name_net()
      sidx = nameMeasurementKeys(shared)
      members = sort(sidx.groups[(:injection, 2)])
      @test length(members) == 4
      srows = [MeasurementRow(i, "", "", key, "", "P", i == 1 ? 600.0 : 100.0, nothing, "test", true, :none, "") for (i, key) in enumerate(members)]
      sset = only(importMeasurementTable(srows, sidx; sigma_by_voltage_level = two).sets)
      @test only(sset.measurements).sigma ≈ sqrt(15.75) rtol = 1e-15
      # a level missing from the case table: with only the 110 kV row the
      # 20 kV measurements take the default of their type
      one = [L(vn_kv = 110, p_mw = 1.5, q_mvar = 2.5, u_kv = 0.4), L(vn_kv = 380, p_mw = 2.0)]
      h1, ho1 = hset(one)
      defaults = measurementStdDevs()
      @test h1["U@Moorau_20"] == defaults[Sparlectra.VmMeas] && ho1["U@Moorau_20"] === :default
      @test h1["P@Moorau_110_20:to"] == defaults[Sparlectra.PflowMeas] && ho1["P@Moorau_110_20:to"] === :default
      @test h1["P@Moorau_110_20:from"] == 1.5 && ho1["P@Moorau_110_20:from"] === :level

      # the summary of the generated table, every sigma cell emptied except
      # one P row, under the 110 kV row alone: the 20 kV measurements count
      # as default, the configured 380 kV level is not in the case
      exact = generateMeasurementTableFromPF(net; noise = false, include_injections = true)
      keep = findfirst(r -> r.type == "P" && occursin(":from", r.terminal_mrid), exact)
      blank(rows; own = nothing) = [MeasurementRow(i, r.timestamp, r.measurement_mrid, r.terminal_mrid, r.node_mrid, r.type, r.value, i == own ? 0.7 : nothing, r.source, r.valid, :none, "") for (i, r) in enumerate(rows)]
      rows = blank(exact; own = keep)
      set = only(importMeasurementTable(rows, index; sigma_by_voltage_level = two).sets)
      set1 = only(importMeasurementTable(rows, index; sigma_by_voltage_level = one).sets)
      summ = Sparlectra._sigma_level_summary(net, set1.measurements; origins = set1.sigma_origin, levels = one)
      @test summ.unused == [380.0] && summ.unlocated == 0
      s110 = only(filter(l -> l.vn_kv == 110.0, summ.levels))
      s20 = only(filter(l -> l.vn_kv == 20.0, summ.levels))
      @test s110.configured && !s20.configured
      @test s110.row == 1 && s110.default == 0 && s110.level == s110.measurements - 1
      @test s20.default == s20.measurements && s20.level == 0
      @test s110.measurements + s20.measurements == length(set1.measurements)
      # the pu sigmas of the report: P over baseMVA, U over Vn, I unset
      @test s110.sigma_pu.p == 1.5 / net.baseMVA && s110.sigma_pu.u == 0.4 / 110.0 && s110.sigma_pu.i === nothing
      @test s20.sigma_pu.p === nothing

      # without a matching level (an empty table or one that names no
      # voltage of the case) the import and the estimate are bit for bit
      # those without the key
      base_set, base_se = begin
        s0 = only(importMeasurementTable(rows, index).sets)
        w = deepcopy(net)
        applyMeasurementSet!(w, s0)
        (s0, with_state_estimation_config(() -> runse!(w); update_net = false))
      end
      for levels in (Sparlectra.SESigmaLevel[], [L(vn_kv = 999, p_mw = 1.0, u_kv = 1.0)])
        s0 = only(importMeasurementTable(rows, index; sigma_by_voltage_level = levels).sets)
        @test [m.sigma for m in s0.measurements] == [m.sigma for m in base_set.measurements]
        w = deepcopy(net)
        applyMeasurementSet!(w, s0)
        se0 = with_state_estimation_config(() -> runse!(w); update_net = false, sigma_by_voltage_level = levels)
        @test se0.objectiveJ === base_se.objectiveJ && se0.voltages == base_se.voltages && se0.residuals == base_se.residuals
      end
      # the same on the SmallGrid table, every second sigma cell emptied so
      # the default applies to those rows: the run without the key (the
      # configuration default) against runs that set it in the active
      # configuration, import and estimate alike
      sg = _smallgrid_import()
      sg_index = cgmesMeasurementKeys(sg)
      sg_rows = [MeasurementRow((f === :sigma && iseven(i) ? nothing : getfield(r, f) for f in fieldnames(MeasurementRow))...) for (i, r) in enumerate(readMeasurementTable(_SMALLGRID_TABLE))]
      sg_run(; kw...) = with_state_estimation_config(; update_net = false, kw...) do
        s = first(importMeasurementTable(sg_rows, sg_index).sets)
        w = deepcopy(sg.net)
        applyMeasurementSet!(w, s)
        (s, runse!(w))
      end
      sg_set, sg_se = sg_run()
      @test any(m -> sg_set.sigma_origin[m.id] === :default, sg_set.measurements)
      for levels in (Sparlectra.SESigmaLevel[], [L(vn_kv = 999, p_mw = 1.0, u_kv = 1.0)])
        s0, se0 = sg_run(sigma_by_voltage_level = levels)
        @test [m.sigma for m in s0.measurements] == [m.sigma for m in sg_set.measurements]
        @test se0.objectiveJ === sg_se.objectiveJ && se0.voltages == sg_se.voltages && se0.residuals == sg_se.residuals
      end
      # the configured table is the default of the import
      set_cfg = with_state_estimation_config(() -> only(importMeasurementTable(rows, index).sets); sigma_by_voltage_level = two)
      @test [m.sigma for m in set_cfg.measurements] == [m.sigma for m in set.measurements]

      # the generator writes the level sigma (U in kV: max(0.4, 0.5 % of the
      # reading) at 110 kV, the 20 kV row without rel_pct its absolute
      # sigma) and draws its noise with it; noise 0 reproduces the state to
      # 1e-6 pu with or without the table
      gen = generateMeasurementTableFromPF(net; noise = false, include_injections = true, sigma_by_voltage_level = two)
      u110 = only(filter(r -> r.type == "U" && r.node_mrid == "Moorau_110", gen))
      @test u110.sigma == max(0.4, 0.005 * u110.value)
      @test only(filter(r -> r.type == "P" && r.terminal_mrid == "Moorau_110_20:to", gen)).sigma == 0.3
      @test only(filter(r -> r.type == "Q" && r.terminal_mrid == "Moorau_Hassel:from", gen)).sigma == 0.4
      for table in (gen, exact)
        _, se_g = _estimate_table(net, table)
        @test se_g.converged && maximum(abs.(se_g.voltages .- Vpf)) < 1e-6
      end
      noisy_two = generateMeasurementTableFromPF(net; seed = 3, sigma_by_voltage_level = two)
      noisy_def = generateMeasurementTableFromPF(net; seed = 3, sigma_by_voltage_level = Sparlectra.SESigmaLevel[])
      plain = generateMeasurementTableFromPF(net; noise = false, sigma_by_voltage_level = Sparlectra.SESigmaLevel[])
      # the same random draw, scaled by each row's written sigma
      @test all(isapprox((a.value - p.value) / a.sigma, (b.value - p.value) / b.sigma; atol = 1e-8) for (a, b, p) in zip(noisy_two, noisy_def, plain))
      @test any(a.sigma != b.sigma for (a, b) in zip(noisy_two, noisy_def))
      # pseudo rows carry the level sigma where the level has a column: the
      # 20 kV load its absolute sigmas, the 110 kV load max(1.5, 0.5 %)
      loads = Set(k for (k, t) in index.terminals if t.kind === :injection && t.detail == "load")
      cut = filter(r -> !(r.terminal_mrid in loads), exact)
      pseudo = generatePseudoMeasurementTable(net, cut; sigma = 10.0, sigma_by_voltage_level = two)
      prow(key, typ) = only(filter(r -> r.terminal_mrid == key && r.type == typ, pseudo))
      @test prow("Ld_20.0@Hassel_20", "P").sigma == 0.3 && prow("Ld_20.0@Hassel_20", "Q").sigma == 0.4
      p110 = prow("Ld_110.0@Roggau_110", "P")
      @test p110.sigma == max(1.5, 0.005 * abs(p110.value))
      @test all(r -> r.sigma == 10.0, generatePseudoMeasurementTable(net, cut; sigma = 10.0, sigma_by_voltage_level = Sparlectra.SESigmaLevel[]))

      # a 380/110 kV transformer of sp_case188: the P row at its 110 kV end
      # takes the 110 kV row, the one at its 380 kV end the 380 kV row
      net188 = load_fixture_net("sp_case188")
      vn188(b) = net188.nodeVec[b].comp.cVN
      idx188 = nameMeasurementKeys(net188)
      is_trafo(t) = t.kind === :branch && (br = net188.branchVec[t.branchIdx]; br.ratio != 0.0 && Set((vn188(br.fromBus), vn188(br.toBus))) == Set((380.0, 110.0)))
      k = minimum(t.branchIdx for (_, t) in idx188.terminals if is_trafo(t))
      ends = [key for (key, t) in idx188.terminals if t.kind === :branch && t.branchIdx == k]
      trows = [MeasurementRow(i, "", "", key, "", "P", 50.0, nothing, "test", true, :none, "") for (i, key) in enumerate(ends)]
      tset = only(importMeasurementTable(trows, idx188; sigma_by_voltage_level = [L(vn_kv = 380, p_mw = 2.0), L(vn_kv = 110, p_mw = 1.0)]).sets)
      side(m) = m.direction === :from ? net188.branchVec[k].fromBus : net188.branchVec[k].toBus
      @test sort!([(vn188(side(m)), m.sigma) for m in tset.measurements]) == [(110.0, 1.0), (380.0, 2.0)]

      # the SV generator on the SmallGrid takes the table the same way: the
      # level sigma at the configured level, its own sigma elsewhere
      result = _smallgrid_import()
      cindex = cgmesMeasurementKeys(result)
      top = maximum(nd.comp.cVN for nd in result.net.nodeVec)
      sv = generateMeasurementTableFromSV(result; noise = false, sigma_by_voltage_level = [L(vn_kv = top, p_mw = 3.0, q_mvar = 4.0, u_kv = 1.0)])
      sv0 = generateMeasurementTableFromSV(result; noise = false, sigma_by_voltage_level = Sparlectra.SESigmaLevel[])
      @test length(sv) == length(sv0)
      at_top(r) = result.net.nodeVec[isempty(r.node_mrid) ? cindex.terminals[r.terminal_mrid].busIdx : cindex.nodes[r.node_mrid]].comp.cVN == top
      @test any(at_top, sv) && all(r.sigma == (r.type == "U" ? 1.0 : r.type == "P" ? 3.0 : 4.0) for r in sv if at_top(r))
      @test all(r.sigma == r0.sigma for (r, r0) in zip(sv, sv0) if !at_top(r))

      # the Web UI measurement generator (generateMeasurementsFromPF, CSV
      # v1) takes the table of the case configuration: a P, I or U row at
      # 110 kV carries the level sigma of its noise-free value and the same
      # noise draw scaled to it; Q at 110 kV (no column), the 20 kV rows and
      # the zero-injection constraints stay as without the key
      mktempdir() do d
        case = joinpath(d, "sp_case14.scf.json")
        cp(joinpath(dirname(@__DIR__), "data", "scf", "sp_case14.scf.json"), case)
        opts = SparlectraApp.MeasurementGeneratorOptions(; noise = true, gross_k = 0.0, tap_steps = 0.0, include_i = true, sigma_u_pct = 0.5, sigma_i_pct = 1.0, sigma_p_pct = 1.0, sigma_q_pct = 1.0, sigma_ia_deg = 0.0, seed = 42)
        gen_set(name) = begin
          out = joinpath(d, name)
          g = redirect_stdout(() -> SparlectraApp._se_generate_measurement_set(case, out, opts), devnull)
          w = load_fixture_net("sp_case14")
          readMeasurementTable!(w; file = out, zero_injection = false)
          truth = Dict(String(m.captures[1]) => parse(Float64, m.captures[2]) for m in eachmatch(r"truth_value,([^,\s]+),(\S+)", read(out, String)))
          (g, Dict(m.id => m for m in w.measurements), truth, w)
        end
        g0, m0, t0, w = gen_set("plain.csv")
        write(Sparlectra.case_config_path(case), "config_version: 1\nscope: case\ncase: sp_case14.scf.json\nstate_estimation:\n  sigma_by_voltage_level:\n    110:\n      p_mw: 1.5\n      i_a: 20\n      u_kv: 0.8\n      rel_pct: 0.5\n")
        g1, m1, t1, _ = gen_set("level.csv")
        @test g0.sigma_level_rows == 0 && g1.sigma_level_rows > 0
        # the same rows, except the current magnitudes the noise drew below
        # zero, which the table does not take (the draw differs with the sigma)
        @test all(k -> startswith(k, "I@"), symdiff(Set(keys(m0)), Set(keys(m1))))
        @test all(t0[k] == t1[k] for k in intersect(keys(t0), keys(t1)))
        lr = L(vn_kv = 110, p_mw = 1.5, i_a = 20, u_kv = 0.8, rel_pct = 0.5)
        vn_of(m) = w.nodeVec[Sparlectra._sigma_level_bus(w, m)].comp.cVN
        takes_level(id, m, t) = vn_of(m) == 110.0 && Sparlectra._sigma_table_type(m.typ) in ("P", "I", "U") && haskey(t, id) && !startswith(id, "ZI")
        pairs = [(id, a, m1[id]) for (id, a) in m0 if haskey(m1, id)]
        lvl = [(id, a, b) for (id, a, b) in pairs if takes_level(id, a, t0)]
        # one row of each level type by hand: sigma from the noise-free value
        # and the noise draw scaled to it
        for typ in (Sparlectra.PflowMeas, Sparlectra.ImagMeas, Sparlectra.VmMeas)
          id, a, b = first(sort!(filter(x -> x[2].typ == typ, lvl); by = first))
          t = t0[id]
          want = typ == Sparlectra.VmMeas ? max(0.8, 0.005 * t * 110.0) / 110.0 : max(typ == Sparlectra.PflowMeas ? 1.5 : 20.0, 0.005 * abs(t))
          @test b.sigma ≈ want rtol = 1e-12
          @test (b.value - t) / b.sigma ≈ (a.value - t) / a.sigma atol = 1e-6
        end
        # the other rows stay: a Q flow at 110 kV, a 20 kV row, a
        # zero-injection constraint
        rest = [(id, a, b) for (id, a, b) in pairs if !takes_level(id, a, t0)]
        for pick in ((id, a) -> a.typ == Sparlectra.QflowMeas && vn_of(a) == 110.0, (id, a) -> vn_of(a) == 20.0, (id, a) -> startswith(id, "ZI"))
          id, a, b = first(sort!(filter(x -> pick(x[1], x[2]), rest); by = first))
          @test b.sigma == a.sigma && b.value == a.value
        end
        # the rows of the written table that take the level sigma
        @test count(takes_level(id, b, t1) for (id, b) in m1) == g1.sigma_level_rows
      end
    end)()
  end

  @testset "measurement table: pseudo rows for unobserved injections" begin
    (function ()
      # sp_case14 (name keys): the injection rows of the pure load buses are
      # cut; the generator writes exactly those members from the setpoints,
      # which equal the solved values at a load bus, so the noise-free state
      # is reproduced
      net = load_fixture_net("sp_case14")
      _, erg = runpf!(net, 40, 1e-10, 0)
      @test erg == 0
      Vpf = Sparlectra.buildVoltageVector(net)
      index = nameMeasurementKeys(net)
      exact = generateMeasurementTableFromPF(net; noise = false, include_injections = true)
      @test isempty(generatePseudoMeasurementTable(net, exact; sigma = 10.0))
      loads = Set(k for (k, t) in index.terminals if t.kind === :injection && t.detail == "load")
      @test length(loads) >= 3
      cut = filter(r -> !(r.terminal_mrid in loads), exact)
      pseudo = generatePseudoMeasurementTable(net, cut; sigma = 10.0, sigma_q = 20.0)
      @test Set((r.terminal_mrid, r.type) for r in pseudo) == Set((k, t) for k in loads for t in ("P", "Q"))
      @test all(r -> r.source == "pseudo;ssh" && r.valid && r.sigma == (r.type == "P" ? 10.0 : 20.0), pseudo)
      set, se = _estimate_table(net, vcat(cut, pseudo))
      @test isempty(set.problems)
      @test se.converged
      @test se.objectiveJ < 1e-8
      @test maximum(abs.(se.voltages .- Vpf)) < 1e-6
      for kw in ((sigma = 0.0,), (sigma = NaN,), (sigma = 1.0, sigma_q = -1.0))
        @test_throws ArgumentError generatePseudoMeasurementTable(net, cut; kw...)
      end

      # SmallGrid (mRID keys): set 2 observes every injection, so nothing is
      # written for it; set 1 has no injection rows, so every member gets P
      # and Q. Cutting the injection rows of some load buses of set 2 gives
      # pseudo rows for exactly their members, and the estimate on the
      # filled set stays consistent (real rows win at sigma 5 MW).
      result = _smallgrid_import()
      cindex = cgmesMeasurementKeys(result)
      rows = readMeasurementTable(_SMALLGRID_TABLE)
      ts1, ts2 = "2026-01-01T12:00Z", "2026-01-01T13:00Z"
      groups = sort!([(k[2], v) for (k, v) in cindex.groups if k[1] === :injection]; by = first)
      nmembers = sum(length(v) for (_, v) in groups)
      full = generatePseudoMeasurementTable(result, rows; sigma = 5.0)
      @test count(r -> r.timestamp == ts1, full) == 2 * nmembers
      @test count(r -> r.timestamp == ts2, full) == 0
      loadbuses = filter(g -> all(m -> cindex.terminals[m].detail == "EnergyConsumer", g[2]), groups)
      @test any(g -> length(g[2]) == 2, loadbuses)
      victims = Set(m for (_, v) in [loadbuses[1:2]; first(filter(g -> length(g[2]) == 2, loadbuses))] for m in v)
      rows2 = filter(r -> r.timestamp == ts2 && !(r.terminal_mrid in victims), rows)
      pseudo2 = generatePseudoMeasurementTable(result, rows2; sigma = 5.0)
      @test Set((r.terminal_mrid, r.type) for r in pseudo2) == Set((m, t) for m in victims for t in ("P", "Q"))
      @test all(r -> r.timestamp == ts2, pseudo2)
      # through the writer and the reader, like any table
      file = test_scratch_path(".table.csv")
      writeMeasurementTable(file, vcat(rows2, pseudo2))
      set2 = only(importMeasurementTable(readMeasurementTable(file), cindex).sets)
      @test isempty(set2.problems)
      @test count(m -> m.typ == Sparlectra.PinjMeas, set2.measurements) == length(groups)
      work = deepcopy(result.net)
      applyMeasurementSet!(work, set2)
      se2 = with_state_estimation_config(update_net = false) do
        runse!(work)
      end
      @test se2.converged
      @test se2.objectiveJ / se2.dof < 1.5
    end)()
  end

  @testset "measurement table: repeated names at one bus" begin
    (function ()
      # two generators and two loads at B2 share their generated names; each
      # gets an addressable key with its order at the bus and stays a member
      # of the bus sum
      net = _shared_name_net()
      index = nameMeasurementKeys(net)
      gkeys = ["Gen_110@B2#1", "Gen_110@B2#2", "Ld_110.0@B2#1", "Ld_110.0@B2#2"]
      @test all(k -> index.terminals[k].kind === :injection && index.terminals[k].busIdx == 2, gkeys)
      @test sort(index.groups[(:injection, 2)]) == sort(gkeys)
      # noise 0 with the injection rows reproduces the solved state
      Vpf = Sparlectra.buildVoltageVector(net)
      rows = generateMeasurementTableFromPF(net; noise = false, include_injections = true)
      set, se = _estimate_table(net, rows)
      @test isempty(set.problems)
      @test se.converged
      @test maximum(abs.(se.voltages .- Vpf)) < 1e-6
      pinj = only(m for m in set.measurements if m.typ == Sparlectra.PinjMeas && m.busIdx == 2)
      @test isapprox(pinj.value, 20.0 + 15.0 - 30.0 - 25.0; atol = 1e-9)
      # a member without a row gives no sum (partial_injection), never a sum
      # over the remaining members
      cut, _ = _estimate_table(net, filter(r -> r.terminal_mrid != "Gen_110@B2#2", rows))
      @test !any(m -> m.typ in (Sparlectra.PinjMeas, Sparlectra.QinjMeas) && m.busIdx == 2, cut.measurements)
      @test cut.counts[:partial_injection] == 6
      # a prosumer without a name has no key, but it stays a member: the sum
      # at its bus is partial instead of being formed without it
      unnamed = deepcopy(net)
      unnamed.prosumpsVec[3].comp.cName = ""
      uindex = nameMeasurementKeys(unnamed)
      @test length(uindex.groups[(:injection, 2)]) == 4
      uset, _ = _estimate_table(unnamed, generateMeasurementTableFromPF(unnamed; noise = false, include_injections = true))
      @test !any(m -> m.typ in (Sparlectra.PinjMeas, Sparlectra.QinjMeas) && m.busIdx == 2, uset.measurements)
      @test uset.counts[:partial_injection] == 6
    end)()
  end

  @testset "measurement table: name keys keep their exact form" begin
    (function ()
      # bus names "_B9" and "B9" are two buses: the canonical mRID form (no
      # leading "_") must not merge them, or a U row lands on the other bus
      net = Net(name = "names", baseMVA = 100.0)
      for b in ("B1", "_B9", "B9", "B7", " B7")
        addBus!(net = net, busName = b, vn_kV = 110.0)
      end
      for b in ("_B9", "B9", "B7", " B7")
        addACLine!(net = net, fromBus = "B1", toBus = b, r_ohm = 1.0, x_ohm = 4.0)
      end
      addProsumer!(net = net, busName = "B1", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "B1")
      addProsumer!(net = net, busName = "_B9", type = "ENERGYCONSUMER", p = 10.0, q = 2.0)
      addProsumer!(net = net, busName = "B9", type = "ENERGYCONSUMER", p = 30.0, q = 6.0)
      _, erg = runpf!(net, 40, 1e-10, 0)
      @test erg == 0
      index = nameMeasurementKeys(net)
      @test index.nodes["_B9"] == net.busDict["_B9"] && index.nodes["B9"] == net.busDict["B9"]
      R(l, node, v) = MeasurementRow(l, "", "", "", node, "U", v, 0.2, "t", true, :none, "")
      set = only(importMeasurementTable([R(1, "_B9", 109.0), R(2, "B9", 108.0), R(3, "B7", 110.0)], index).sets)
      bus_of(v) = only(m.busIdx for m in set.measurements if isapprox(m.value * 110.0, v))
      @test bus_of(109.0) == net.busDict["_B9"] && bus_of(108.0) == net.busDict["B9"]
      # "B7" and " B7" read as one key (the reader strips cells): ambiguous,
      # a named problem row, never one of the two buses
      p = only(set.problems)
      @test p.line == 3 && p.problem === :not_modelled && occursin("2 buses", p.detail)
      # the mRID index keeps the canonical form: "_tn1" finds "tn1"
      T = Sparlectra.MeasurementTarget
      cidx = Sparlectra.MeasurementKeyIndex(Dict{String,T}(), Dict("tn1" => 1), Dict{String,Tuple{String,String,String}}(), Dict{Tuple{Symbol,Int},Vector{String}}(), [110.0])
      @test only(only(importMeasurementTable([R(1, "_tn1", 110.0)], cidx).sets).measurements).busIdx == 1
    end)()
  end

  @testset "measurement table: two shunts at one bus" begin
    (function ()
      # the power flow stamps both shunts: two shunts of 6 and 3 MVAr solve
      # like one of 9 MVAr
      net = _shared_name_net(shunt_q = [6.0, 3.0])
      Vpf = Sparlectra.buildVoltageVector(net)
      @test maximum(abs.(Vpf .- Sparlectra.buildVoltageVector(_shared_name_net(shunt_q = [9.0])))) < 1e-10
      # the table sums both into one ShuntQMeas, and the estimator predicts
      # both: noise 0 reproduces the state
      rows = generateMeasurementTableFromPF(net; noise = false, include_injections = true)
      set, se = _estimate_table(net, rows)
      @test isempty(set.problems)
      @test se.converged
      @test se.objectiveJ < 1e-8
      @test maximum(abs.(se.voltages .- Vpf)) < 1e-6
      q9 = imag(abs2(Vpf[3]) * conj(9.0im / net.baseMVA) * net.baseMVA)
      @test isapprox(only(m for m in set.measurements if m.typ == Sparlectra.ShuntQMeas).value, q9; atol = 1e-9)
      # the measurement generator writes one ShuntQMeas per bus with the same sum
      gen = filter(m -> m.typ == Sparlectra.ShuntQMeas, generateMeasurementsFromPF(net; includeShuntQ = true, noise = false))
      @test length(gen) == 1 && isapprox(gen[1].value, q9; atol = 1e-9)
    end)()
  end

  @testset "measurement table: bus balances, angles and the shunt bay current" begin
    (function ()
      # one shunt at B3: its bay current is a modelled quantity
      net = _shared_name_net(shunt_q = [6.0])
      V = Sparlectra.buildVoltageVector(net)
      busname = Dict(b => n for (n, b) in net.busDict)
      ref = generateMeasurementsFromPF(net; includeVa = true, includeImag = true, includeIa = true, noise = false)
      bay = Measurement(typ = Sparlectra.ImagMeas, value = 1.0, sigma = 1.0, busIdx = 3)
      ibay = abs(Sparlectra._end_current_A(net, bay, V, nothing))
      Row = NamedTuple{(:timestamp, :measurement_mrid, :terminal_mrid, :node_mrid, :type, :value, :sigma, :source, :valid),Tuple{String,String,String,String,String,Float64,Float64,String,Bool}}
      tname = Dict(Sparlectra.PflowMeas => "P", Sparlectra.QflowMeas => "Q", Sparlectra.ImagMeas => "I", Sparlectra.IaMeas => "Ia", Sparlectra.PinjMeas => "P", Sparlectra.QinjMeas => "Q", Sparlectra.VaMeas => "Va")
      rows = Row[]
      for m in ref
        if m.typ == Sparlectra.VmMeas
          vn = net.nodeVec[m.busIdx].comp.cVN
          push!(rows, Row(("", "", "", busname[m.busIdx], "U", m.value * vn, m.sigma * vn, "test", true)))
        elseif m.busIdx !== nothing
          push!(rows, Row(("", "", "", busname[m.busIdx], tname[m.typ], m.value, m.sigma, "test", true)))
        else
          push!(rows, Row(("", "", net.branchVec[m.branchIdx].comp.cName * ":" * string(m.direction), "", tname[m.typ], m.value, m.sigma, "test", true)))
        end
      end
      shunt = only(net.shuntVec).comp.cName
      push!(rows, Row(("", "", shunt, "", "I", ibay, 0.5, "test", true)))
      # an angle at an injection is no modelled quantity
      gen = only(k for (k, t) in nameMeasurementKeys(net).terminals if t.kind === :injection && t.busIdx == 1)
      push!(rows, Row(("", "", gen, "", "Ia", 10.0, 1.0, "test", true)))
      file = test_scratch_path(".table.csv")
      writeMeasurementTable(file, rows)
      imp = importMeasurementTable(readMeasurementTable(file), nameMeasurementKeys(net))
      set = only(imp.sets)
      @test get(set.counts, :unsupported, 0) == 1 && length(set.problems) == 1
      key(m) = (m.typ, m.busIdx, m.branchIdx, m.direction)
      got = Dict(key(m) => m for m in set.measurements)
      # every type of the generator arrives with the value it was written with
      for m in ref
        g = got[key(m)]
        @test g.value == (m.typ == Sparlectra.VmMeas ? (m.value * net.nodeVec[m.busIdx].comp.cVN) / net.nodeVec[m.busIdx].comp.cVN : m.value)
      end
      @test Set(m.typ for m in set.measurements) == Set([Sparlectra.VmMeas, Sparlectra.VaMeas, Sparlectra.PinjMeas, Sparlectra.QinjMeas, Sparlectra.PflowMeas, Sparlectra.QflowMeas, Sparlectra.ImagMeas, Sparlectra.IaMeas])
      sb = only(m for m in set.measurements if m.typ == Sparlectra.ImagMeas && m.busIdx == 3)
      @test sb.value == ibay && sb.sigma == 0.5
      # the bus balance P of B2 and its prosumer sum would both be redundant
      # meters; here only the balance is written, and the noise-free set
      # reproduces the solved state
      work = deepcopy(net)
      applyMeasurementSet!(work, set; zero_injection = false)
      se = with_state_estimation_config(update_net = false) do
        runse!(work)
      end
      @test se.converged && se.objectiveJ < 1e-8
      @test maximum(abs.(se.voltages .- V)) < 1e-7
      # two shunts at one bus: their bay current is not the sum of the
      # magnitudes, so I rows there are unsupported
      net2 = _shared_name_net(shunt_q = [6.0, 3.0])
      skeys = sort!([k for (k, t) in nameMeasurementKeys(net2).terminals if t.kind === :shunt])
      f2 = test_scratch_path(".table.csv")
      # a zero-injection balance keeps its sigma and gets the ZI id
      writeMeasurementTable(f2, vcat([Row(("", "", k, "", "I", 1.0, 0.5, "test", true)) for k in skeys], [Row(("", "", "", "B3", "P", 0.0, 0.05, "converted;v1;zero_injection", false))]))
      s2 = only(importMeasurementTable(readMeasurementTable(f2), nameMeasurementKeys(net2)).sets)
      @test s2.counts[:unsupported] == 2
      zi = only(s2.measurements)
      @test zi.id == "ZI_PINJ_bus_3" && zi.typ == Sparlectra.PinjMeas && zi.sigma == 0.05 && !zi.active
      # a bus link end is a link flow target: a P row at the to end becomes
      # a link measurement with the opposite sign (positive from the link's
      # from bus), Q at the from end keeps its sign, I is unsupported; the
      # set written back as a table carries both rows at the from end, an
      # open link resolves to not_modelled
      net3 = deepcopy(net)
      net3._locked = false
      addBus!(net = net3, busName = "B4", vn_kV = 110.0)
      addLink!(net = net3, fromBus = "B3", toBus = "B4")
      lname = net3.linkVec[1].cName
      idx3 = nameMeasurementKeys(net3)
      t = idx3.terminals[lname * ":to"]
      @test t.kind === :link && t.branchIdx == 1 && t.direction === :to
      f3 = test_scratch_path(".table.csv")
      writeMeasurementTable(f3, [Row(("", "", lname * ":to", "", "P", 5.0, 0.5, "test", true)), Row(("", "", lname * ":from", "", "Q", -2.0, 0.5, "test", true)), Row(("", "", lname * ":from", "", "I", 10.0, 1.0, "test", true))])
      s3 = only(importMeasurementTable(readMeasurementTable(f3), idx3).sets)
      @test s3.counts[:unsupported] == 1
      lp = only(m for m in s3.measurements if m.typ == Sparlectra.PflowMeas)
      lq = only(m for m in s3.measurements if m.typ == Sparlectra.QflowMeas)
      @test lp.linkIdx == 1 && lp.value == -5.0 && lp.direction === :from && lp.active && lp.branchIdx === nothing
      @test lq.linkIdx == 1 && lq.value == -2.0
      back3 = Sparlectra._measurement_table_from_measurements(net3, s3.measurements, idx3)
      @test isempty(back3.problems)
      @test sort!([(r.terminal_mrid, r.type, r.value) for r in back3.rows]) == [(lname * ":from", "P", "-5.0"), (lname * ":from", "Q", "-2.0")]
      net3.linkVec[1].status = 0
      @test nameMeasurementKeys(net3).terminals[lname * ":from"].kind === :not_modelled
      # CGMES: the closed breaker of sp_casePST is a bus link; its terminals
      # resolve to the link, one at each end
      pst = redirect_stdout(() -> importCGMES(path = cgmes_fixture_dir("sp_casePST")), devnull)
      lt = [t for t in values(cgmesMeasurementKeys(pst).terminals) if t.kind === :link]
      @test length(lt) == 2 && Set(t.direction for t in lt) == Set([:from, :to]) && all(t -> t.branchIdx == 1, lt)
    end)()
  end

  @testset "measurement table: U in kV gives back the pu value bit for bit" begin
    (function ()
      # a Float64 kV value cannot give back every pu value; the kV text can
      net = load_fixture_net("sp_case14")
      busname = Dict(b => n for (n, b) in net.busDict)
      rng = Sparlectra.Random.Xoshiro(7)
      Row = NamedTuple{(:timestamp, :measurement_mrid, :terminal_mrid, :node_mrid, :type, :value, :sigma, :source, :valid),Tuple{String,String,String,String,String,String,String,String,Bool}}
      pus = Tuple{Int,Float64,Float64}[]
      rows = Row[]
      for k = 1:400
        b = mod1(k, length(net.nodeVec))
        vn = net.nodeVec[b].comp.cVN
        v = 0.9 + 0.2 * rand(rng)
        sg = 0.001 + 0.01 * rand(rng)
        push!(pus, (b, v, sg))
        # a fresh timestamp per row keeps every bus measured once per set
        push!(rows, Row((string("2026-01-01T", lpad(k ÷ 60, 2, '0'), ":", lpad(k % 60, 2, '0'), "Z"), "", "", busname[b], "U", Sparlectra._measurement_table_kv_text(v, vn), Sparlectra._measurement_table_kv_text(sg, vn), "test", true)))
      end
      # some values need more digits than a Float64 prints
      @test any(r -> Sparlectra._measurement_table_digits(r.value) > 17, rows)
      file = test_scratch_path(".table.csv")
      writeMeasurementTable(file, rows)
      read_rows = readMeasurementTable(file)
      # the long texts survive a read and a write byte for byte
      file2 = test_scratch_path(".table.csv")
      writeMeasurementTable(file2, read_rows)
      @test read(file) == read(file2)
      imp = importMeasurementTable(read_rows, nameMeasurementKeys(net))
      @test length(imp.sets) == length(pus)
      @test all(only(s.measurements).value == v && only(s.measurements).sigma == sg && only(s.measurements).busIdx == b for (s, (b, v, sg)) in zip(imp.sets, pus))
      # readMeasurementTable! applies one set, the first or the one named
      w = deepcopy(net)
      s5 = readMeasurementTable!(w; file = file, timestamp = "2026-01-01T00:05:00Z", zero_injection = false)
      @test s5.timestamp == "2026-01-01T00:05Z" && only(w.measurements).value == pus[5][2]
      @test readMeasurementTable!(w; file = file, zero_injection = false).timestamp == "2026-01-01T00:01Z"
      @test_throws ArgumentError readMeasurementTable!(w; file = file, timestamp = "2027-01-01T00:00Z")
      # a short cell is divided as the Float64 it parses to, as before
      @test Sparlectra._measurement_table_pu("112.5", 112.5, 110.0) == 112.5 / 110.0
      @test_throws ArgumentError writeMeasurementTable(test_scratch_path(".table.csv"), [Row(("", "", "", "B", "U", "abc", "", "test", true))])
    end)()
  end

  @testset "measurement table: converter from the measurement CSV v1" begin
    (function ()
      net = load_fixture_net("sp_case14")
      _, erg = runpf!(net, 40, 1e-10, 0)
      @test erg == 0
      # a v1 set with every type the table carries: voltages and angles,
      # bus balances, flows, currents and current angles at branch ends,
      # zero-injection rows with a sigma of their own, an inactive row, and
      # one row the table cannot take twice (a duplicate quantity)
      ms = generateMeasurementsFromPF(net; noise = true, includeVa = true, includeImag = true, includeIa = true, rng = Sparlectra.Random.Xoshiro(5))
      passive = findPassiveBuses(net)
      filter!(m -> !(m.typ in (Sparlectra.PinjMeas, Sparlectra.QinjMeas) && m.busIdx in passive), ms)
      addZeroInjectionMeasurements!(ms; net = net, busIdxs = passive, sigma = 0.05)
      ms[3] = Measurement(typ = ms[3].typ, value = ms[3].value, sigma = ms[3].sigma, active = false, busIdx = ms[3].busIdx, id = ms[3].id)
      push!(ms, ms[1])
      empty!(net.measurements)
      append!(net.measurements, ms)
      d = mktempdir()
      v1 = joinpath(d, "sp_case14.v1.csv")
      truth = ["truth_value,$(m.id),$(repr(m.value))" for m in ms[1:5]]
      writeMeasurementsCSV(net; file = v1, headerComments = ["case: sp_case14.scf.json", "critical_rows: $(ms[2].id) $(ms[4].id)"], footerComments = truth)
      tbl = joinpath(d, "sp_case14.table.csv")
      res = convertMeasurementsV1(v1, net; out = tbl)
      # the duplicate and the negative currents the noise drew are listed
      # with their lines, in the result and the table
      neg = [k for (k, m) in enumerate(ms) if m.typ == Sparlectra.ImagMeas && m.value < 0.0]
      @test !isempty(neg)
      @test [(p.line, p.problem) for p in res.problems] == vcat([(k + 4, :bad_value) for k in neg], [(length(ms) + 4, :duplicate)])
      @test res.rows == length(ms) - 1 - length(neg)
      lines = readlines(tbl)
      @test lines[1] == "# sparlectra-measurement-table v1" && lines[2] == "# case: sp_case14.scf.json"
      @test any(l -> startswith(l, "# not converted: line $(length(ms) + 4) ") && occursin("duplicate", l), lines)
      # comment ids follow the rows
      @test "# critical_rows: $(res.ids[ms[2].id]) $(res.ids[ms[4].id])" in lines
      @test all("# truth_value,$(res.ids[m.id]),$(repr(m.value))" in lines for m in ms[1:5])
      @test count(r -> occursin("zero_injection", r.source), readMeasurementTable(tbl)) == 2 * length(passive)
      # back in pu: the import gives the v1 measurements bit for bit, in
      # their order, with the ids the converter announced
      # (the v1 side as the converter's private parse reads it; the public
      # v1 reader was removed in 0.50.0)
      ref_ms = Sparlectra._read_measurements_csv_v1(net, v1).measurements
      @test ref_ms == ms
      a = [m for (k, m) in enumerate(ref_ms[1:(end-1)]) if !(k in neg)]
      set = only(importMeasurementTable(readMeasurementTable(tbl), nameMeasurementKeys(net)).sets)
      b = set.measurements
      @test isempty(set.problems) && length(a) == length(b)
      @test all(x.typ == y.typ && x.value === y.value && x.sigma === y.sigma && x.active == y.active && x.busIdx == y.busIdx && x.branchIdx == y.branchIdx && x.direction == y.direction for (x, y) in zip(a, b))
      @test all(res.ids[x.id] == y.id for (x, y) in zip(a, b))
      @test count(m -> startswith(m.id, "ZI_") && m.sigma == 0.05, b) == 2 * length(passive)
      # the same estimate and the same eliminations
      _converted_same_estimate(net, a, b, res.ids)
      # a CGMES case keyed by mRID: the first SmallGrid set with its gross
      # errors, written as v1 with mRID bus references and converted back
      result = _smallgrid_import()
      index = cgmesMeasurementKeys(result)
      sg = first(importMeasurementTable(readMeasurementTable(_SMALLGRID_TABLE), index).sets)
      sgnet = result.net
      empty!(sgnet.measurements)
      append!(sgnet.measurements, sg.measurements)
      v1sg = joinpath(d, "smallgrid.v1.csv")
      writeMeasurementsCSV(sgnet; file = v1sg, busReference = :mrid)
      ressg = convertMeasurementsV1(v1sg, sgnet; out = joinpath(d, "smallgrid.table.csv"), index = index)
      # the group sums of the fixture are bus balances in v1 and stay so
      @test isempty(ressg.problems)
      refsg_ms = Sparlectra._read_measurements_csv_v1(sgnet, v1sg).measurements
      setsg = only(importMeasurementTable(readMeasurementTable(ressg.file), index).sets)
      @test length(refsg_ms) == length(setsg.measurements)
      @test all(x.typ == y.typ && x.value === y.value && x.sigma === y.sigma && x.active == y.active for (x, y) in zip(refsg_ms, setsg.measurements))
      elim = _converted_same_estimate(sgnet, refsg_ms, setsg.measurements, ressg.ids)
      @test elim >= 1
      # a link flow goes to a terminal row of its link; rows the table has
      # no place for (a current on a link) are listed, not dropped
      link = Measurement(typ = Sparlectra.PflowMeas, value = 1.0, sigma = 0.1, direction = :from, id = "link_row", linkIdx = 1)
      ilink = Measurement(typ = Sparlectra.ImagMeas, value = 1.0, sigma = 0.1, direction = :from, id = "link_i", linkIdx = 1)
      lnet = deepcopy(net)
      lnet._locked = false
      addBus!(net = lnet, busName = "LB", vn_kV = lnet.nodeVec[1].comp.cVN)
      addLink!(net = lnet, fromBus = lnet.nodeVec[2].comp.cName, toBus = "LB")
      conv = Sparlectra._measurement_table_from_measurements(lnet, [link, ilink], nameMeasurementKeys(lnet))
      @test length(conv.rows) == 1 && only(conv.rows).type == "P" && only(conv.problems).problem === :unsupported
    end)()
  end

  @testset "measurement table: the Web UI generator writes the table" begin
    (function ()
      # the generator keeps its set construction (currents, current angles,
      # gross errors, zero-injection constraints) and writes the set as a
      # table, and the SE service runs that table
      mktempdir() do d
        for (case, label) in ((joinpath(dirname(@__DIR__), "data", "scf", "sp_case14.scf.json"), "sp14"), (_SMALLGRID_DIR, "smallgrid"))
          opts = SparlectraApp.MeasurementGeneratorOptions(; noise = true, gross_k = 8.0, gross_count = 2, tap_steps = 0.0, include_i = true, sigma_u_pct = 0.5, sigma_i_pct = 1.0, sigma_p_pct = 1.0, sigma_q_pct = 1.0, sigma_ia_deg = 0.1, passive_as_zi = true, seed = 7)
          tbl = joinpath(d, "$(label).table.csv")
          g = redirect_stdout(() -> SparlectraApp._se_generate_measurement_set(case, tbl, opts), devnull)
          @test readline(tbl) == "# sparlectra-measurement-table v1" && isMeasurementTable(tbl)
          @test g.rows == length(g.measurements) - length(g.not_written)
          @test all(p -> p.problem === :bad_value, g.not_written)
          @test occursin("# truth_value,", read(tbl, String)) && occursin("zero_injection", read(tbl, String))
          # the table carries the generated set: value, sigma and flag of
          # every written row bit for bit, in the order of the set
          imported = SparlectraApp._se_import_case(case, Sparlectra.load_sparlectra_config())
          set = first(importMeasurementTable(readMeasurementTable(tbl), Sparlectra._measurement_key_index(imported)).sets)
          skipped = Set(p.line for p in g.not_written)
          written = [m for (k, m) in enumerate(g.measurements) if !(k in skipped)]
          @test length(written) == length(set.measurements)
          @test all(x.typ == y.typ && x.value === y.value && x.sigma === y.sigma && x.active == y.active for (x, y) in zip(written, set.measurements))
          # the SE service runs the table
          rb = run_with_expected_warnings(() -> SparlectraApp.to_dict(redirect_stdout(devnull) do
            SparlectraApp._run_state_estimation_service(case, Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, joinpath(d, "$(label)_table"), "$(label)_table", tbl)
          end), [r"topology precheck"])
          @test rb["status"] == "succeeded"
        end
      end
    end)()
  end

  @testset "measurement table: SE service switch" begin
    (function ()
      mktempdir() do d
        # off: the table is refused and the key is named
        cfg_off = joinpath(d, "table_off.yaml")
        write(cfg_off, "config_version: 1\nstate_estimation:\n  measurement_table: false\n")
        off = redirect_stdout(devnull) do
          SparlectraApp._run_state_estimation_service(_SMALLGRID_DIR, cfg_off, joinpath(d, "off"), "table_off", _SMALLGRID_TABLE)
        end
        d_off = SparlectraApp.to_dict(off)
        @test d_off["status"] == "failed"
        @test occursin("state_estimation.measurement_table", d_off["message"])
        # on: the first set is estimated, the import report is an artifact.
        # The sigma table names only a level the SmallGrid does not have, so
        # the numbers stay those without it and the run reports the level
        cfg = joinpath(d, "table_on.yaml")
        write(cfg, "config_version: 1\nstate_estimation:\n  measurement_table: true\n  sigma_by_voltage_level:\n    999:\n      p_mw: 1.0\n")
        on = redirect_stdout(devnull) do
          SparlectraApp._run_state_estimation_service(_SMALLGRID_DIR, cfg, joinpath(d, "on"), "table_on", _SMALLGRID_TABLE)
        end
        d_on = SparlectraApp.to_dict(on)
        @test d_on["status"] == "succeeded"
        @test d_on["metadata"]["se_measurement_rows"] == 862
        report = read(joinpath(d, "on", "measurement_import.txt"), String)
        @test occursin("2 set(s)", report) && occursin("unknown_key 1", report)
        # per nominal voltage: the measurements and where their sigma came
        # from, in the metadata, run.log and se_diagnostics.md; the 999 kV
        # row is a note, not an error
        lvl = d_on["metadata"]["se_sigma_by_voltage_level"]
        sg_levels = sort!(unique([nd.comp.cVN for nd in _smallgrid_import().net.nodeVec]); rev = true)
        @test [r["vn_kv"] for r in lvl] == sg_levels
        @test all(r -> r["measurements"] == r["sigma_row"] + r["sigma_level"] + r["sigma_default"] && r["sigma_level"] == 0 && !r["configured"], lvl)
        @test sum(r["measurements"] for r in lvl) == 862
        @test d_on["metadata"]["se_sigma_levels_not_in_case"] == [999.0]
        on_log = read(joinpath(d, "on", "run.log"), String)
        @test occursin("999.0 kV is configured in `state_estimation.sigma_by_voltage_level`, but no bus of this case has that nominal voltage", on_log)
        @test occursin(string(sg_levels[1], " kV: "), on_log)
        @test occursin("## Measurement sigma by voltage level", read(joinpath(d, "on", "se_diagnostics.md"), String))
        # state_estimation.zero_injection on (the default): the passive buses
        # without an injection row get ZI pairs, counted in the report, the
        # log and the metadata, and written to the run's measurements.csv
        # as node rows with the source token zero_injection; off adds none
        # and estimates the 862 rows of the table alone
        zi_on = d_on["metadata"]["se_zero_injection_rows"]
        @test zi_on > 0 && iseven(zi_on)
        @test occursin("zero-injection rows $(zi_on)", report)
        @test occursin("zero injection: $(zi_on) row(s) added", read(joinpath(d, "on", "run.log"), String))
        art = joinpath(d, "on", "measurements.csv")
        art_rows = readMeasurementTable(art)
        @test count(r -> r.source == "se_run;zero_injection", art_rows) == zi_on
        # the artifact is the estimated set as a table: read back with the
        # run's mRID index it gives the first set of the input table plus the
        # zero-injection rows, bitwise in value and sigma (redundant U at
        # busbar sections and bus sums of several prosumers included), and
        # an estimation on the artifact gives the same J
        @test length(unique(r.timestamp for r in art_rows)) == 1
        back = _smallgrid_import()
        back_set = readMeasurementTable!(back.net; file = art, index = cgmesMeasurementKeys(back), zero_injection = false)
        @test isempty(back_set.problems)
        ref = _smallgrid_import()
        applyMeasurementSet!(ref.net, first(importMeasurementTable(readMeasurementTable(_SMALLGRID_TABLE), cgmesMeasurementKeys(ref)).sets); zero_injection = true)
        mkey(m) = (m.id, m.typ, m.busIdx, m.branchIdx, m.direction, m.active, m.zero_injection, reinterpret(UInt64, m.value), reinterpret(UInt64, m.sigma))
        @test sort!(mkey.(ref.net.measurements)) == sort!(mkey.(back.net.measurements))
        on2 = redirect_stdout(devnull) do
          SparlectraApp._run_state_estimation_service(_SMALLGRID_DIR, cfg, joinpath(d, "on2"), "table_on2", art)
        end
        d_on2 = SparlectraApp.to_dict(on2)
        @test d_on2["status"] == "succeeded"
        @test d_on2["metadata"]["se_objective"] === d_on["metadata"]["se_objective"]
        @test d_on2["metadata"]["se_zero_injection_rows"] == 0
        # the verdict: one line per table row, joined by its line number
        verdict = run_csv(joinpath(d, "on", "measurements_verdict.csv"))
        @test length(verdict.rows) == length(art_rows)
        @test [Int(verdict.number(first(row))) for row in verdict.rows] == [r.line for r in art_rows]
        cfg_zoff = joinpath(d, "table_zi_off.yaml")
        write(cfg_zoff, "config_version: 1\nstate_estimation:\n  measurement_table: true\n  zero_injection: off\n")
        zoff = redirect_stdout(devnull) do
          SparlectraApp._run_state_estimation_service(_SMALLGRID_DIR, cfg_zoff, joinpath(d, "zoff"), "table_zoff", _SMALLGRID_TABLE)
        end
        d_zoff = SparlectraApp.to_dict(zoff)
        @test d_zoff["status"] == "succeeded"
        @test d_zoff["metadata"]["se_measurement_rows"] == 862
        @test d_zoff["metadata"]["se_zero_injection_rows"] == 0
        @test occursin("zero-injection rows 0", read(joinpath(d, "zoff", "measurement_import.txt"), String))
        @test !any(r -> occursin("zero_injection", r.source), readMeasurementTable(joinpath(d, "zoff", "measurements.csv")))
        @test !haskey(d_on["metadata"], "se_duplicate_rows")
        # a second U row of the first set on a busbar section of a measured
        # node is redundancy the table allows: no duplicate finding for the run
        index = cgmesMeasurementKeys(_smallgrid_import())
        rows = readMeasurementTable(_SMALLGRID_TABLE)
        t1 = rows[1].timestamp
        tn_row = first(r for r in rows if r.timestamp == t1 && r.type == "U" && r.valid)
        bbs = first(sort!([k for (k, b) in index.nodes if b == index.nodes[tn_row.node_mrid] && k != tn_row.node_mrid && !any(r -> r.node_mrid == k, rows)]))
        red = joinpath(d, "redundant.table.csv")
        writeMeasurementTable(red, vcat(rows, [MeasurementRow(0, t1, "", "", bbs, "U", tn_row.value, tn_row.sigma, "test", true, :none, "")]))
        r2 = SparlectraApp.to_dict(redirect_stdout(devnull) do
          SparlectraApp._run_state_estimation_service(_SMALLGRID_DIR, cfg, joinpath(d, "red"), "table_red", red)
        end)
        @test r2["status"] == "succeeded"
        @test r2["metadata"]["se_measurement_rows"] == 863
        @test !haskey(r2["metadata"], "se_duplicate_rows")
      end
    end)()
  end
end
