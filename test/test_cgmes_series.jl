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

# file: test/test_cgmes_series.jl
# purpose: series check of CGMES snapshots (checkCGMESSeries, qsts.series.*)
#          on snapshots derived in memory from the ENTSO-E SmallGrid fixture:
#          renamed topological nodes and lines are matched again, a removed
#          machine is churn and not structure, removed lines lower the
#          structural overlap and split the series below the limit, a time
#          step above max_gap splits it, the churn limit only flags, report
#          and CSV name every pair, changed line parameters count as
#          parameter drift while switch states, set points and dispatch do not,
#          a power electronics connection and its unit count as injection
#          churn, and a store and a path take model key and time from the
#          headers.

using Dates
import ZipArchives
isdefined(@__MODULE__, :_cgmes_synth3w_dir) || include(joinpath(@__DIR__, "cgmes_synth3w_fixture.jl"))
isdefined(@__MODULE__, :_cgmes_synth_eb_dir) || include(joinpath(@__DIR__, "cgmes_synth_eb_fixture.jl"))

# give one object a new mRID and point every reference at it, as an
# exporter that numbers its objects anew would write it
function _series_rename!(store, old::String, new::String)
  obj = store.objects[old]
  delete!(store.objects, old)
  store.objects[new] = Sparlectra.CGMESImporter.CIMObject(new, obj.class, obj.attrs, obj.refs, obj.source)
  ids = store.byclass[obj.class]
  ids[findfirst(==(old), ids)] = new
  for o in values(store.objects), (k, v) in o.refs
    v == old && (o.refs[k] = new)
  end
  return store
end

# remove one piece of equipment together with its terminals
function _series_drop!(store, mrid::String)
  obj = store.objects[mrid]
  delete!(store.objects, mrid)
  filter!(!=(mrid), store.byclass[obj.class])
  for t in collect(values(store.objects))
    if t.class === :Terminal && get(t.refs, :ConductingEquipment, "") == mrid
      delete!(store.objects, t.mrid)
      filter!(!=(t.mrid), store.byclass[:Terminal])
    end
  end
  return store
end

# add one PowerElectronicsConnection (CGMES 3.0) with one unit (default a
# PhotoVoltaicUnit) on the node of the first load by mRID, SSH p and q in
# load sign; a copy of the store changes them with _series_set_pec!
function _series_add_pec!(store, p::Float64, q::Float64; unit::Symbol = :PhotoVoltaicUnit, unit_attrs::Dict{Symbol,String} = Dict{Symbol,String}(:minP => "0", :maxP => "20"))
  C = Sparlectra.CGMESImporter
  host = first(sort!(copy(store.byclass[:EnergyConsumer])))
  ht = only(t for t in C.objectsOf(store, :Terminal) if get(t.refs, :ConductingEquipment, "") == host)
  add!(m, cls, attrs, refs) = (store.objects[m] = C.CIMObject(m, cls, attrs, refs, "pec_series"); push!(get!(Vector{String}, store.byclass, cls), m))
  add!("_pec_series", :PowerElectronicsConnection, Dict{Symbol,String}(:name => "PEC_SERIES", :p => string(p), :q => string(q), :minQ => "-5", :maxQ => "5"), Dict{Symbol,String}())
  trefs = Dict{Symbol,String}(:ConductingEquipment => "_pec_series")
  for r in (:TopologicalNode, :ConnectivityNode)
    haskey(ht.refs, r) && (trefs[r] = ht.refs[r])
  end
  add!("_pec_series_t", :Terminal, Dict{Symbol,String}(:sequenceNumber => "1", :connected => "true"), trefs)
  add!("_pec_series_unit", unit, merge(Dict{Symbol,String}(:name => "PEC_SERIES_UNIT"), unit_attrs), Dict{Symbol,String}(:PowerElectronicsConnection => "_pec_series"))
  return store
end

function _series_set_pec!(store, p::Float64, q::Float64)
  o = store.objects["_pec_series"]
  o.attrs[:p] = string(p)
  o.attrs[:q] = string(q)
  return store
end

# the classes of the structural overlap as documented for
# qsts.series.min_structural_overlap; TopologicalNode is not among them
const _SERIES_TEST_STRUCTURAL = [:ACLineSegment, :PowerTransformer, :PowerTransformerEnd, :BusbarSection, :Breaker, :Disconnector, :ConnectivityNode]

_series_count(store, classes) = sum(count(m -> !(m in store.boundary), get(store.byclass, c, String[])) for c in classes)

function run_cgmes_series_tests()
  C = Sparlectra.CGMESImporter
  dir = joinpath(@__DIR__, "fixtures", "cgmes_smallgrid")
  files = sort!(filter(f -> endswith(f, ".zip"), readdir(dir; join = true)))
  s1 = redirect_stdout(() -> C.loadCGMES(files), devnull)
  # s2: every topological node of the model and the first ten lines get new
  # mRIDs, one machine is gone; boundary objects keep their mRIDs (they come
  # from the boundary set, and the check does not count them)
  s2 = deepcopy(s1)
  tns = sort!(filter(m -> !(m in s2.boundary), s2.byclass[:TopologicalNode]))
  for m in tns
    _series_rename!(s2, m, "_new" * m)
  end
  for m in sort!(copy(s2.byclass[:ACLineSegment]))[1:10]
    _series_rename!(s2, m, "_new" * m)
  end
  _series_drop!(s2, first(sort!(copy(s2.byclass[:SynchronousMachine]))))
  # and a power electronics connection with its unit is new
  _series_add_pec!(s2, -12.0, -2.0)
  # s3: a fifth of the lines removed, the connection at another p and q
  s3 = deepcopy(s2)
  _series_set_pec!(s3, -18.0, -4.0)
  dropped = sort!(copy(s3.byclass[:ACLineSegment]))[1:div(length(s3.byclass[:ACLineSegment]), 5)]
  for m in dropped
    _series_drop!(s3, m)
  end
  t0 = DateTime(2026, 1, 1, 12, 0, 0)
  inputs = [CGMESSeriesInput(s1; label = "s1", scenario_time = t0, model_key = "m"),
            CGMESSeriesInput(s2; label = "s2", scenario_time = t0 + Minute(60), model_key = "m"),
            CGMESSeriesInput(s3; label = "s3", scenario_time = t0 + Minute(120), model_key = "m")]

  @testset "CGMES series check: matching and splits" begin (function ()
    r = checkCGMESSeries(inputs; config = QstsSeriesConfig())
    p12, p23 = r.pairs
    d12 = Dict(d.class => d for d in p12.classes)
    # renumbered lines and nodes are found again, the overlap stays 1
    @test d12[:ACLineSegment].renamed == 10 && d12[:ACLineSegment].only_a == 0 && d12[:ACLineSegment].only_b == 0
    @test d12[:TopologicalNode].kind === :info
    @test d12[:TopologicalNode].renamed == length(tns) && d12[:TopologicalNode].ambiguous == 0
    @test d12[:TopologicalNode].only_a == 0 && d12[:TopologicalNode].only_b == 0
    @test p12.structural_overlap == 1.0 && p12.split_reason === nothing
    @test p12.gap_minutes == 60.0
    # new mRIDs and removed objects are no parameter drift
    @test all(d -> d.drifted == 0, vcat(p12.classes, p23.classes))
    # the removed machine is churn, not structure
    @test d12[:SynchronousMachine].only_a == 1 && d12[:SynchronousMachine].kind === :injection
    @test p12.injection_churn > 0.0 && !p12.churn_flag
    # the new power electronics connection and its unit are churn as well:
    # of the n1 injection objects of s1 one is gone and two are new
    @test (d12[:PowerElectronicsConnection].kind, d12[:PowerElectronicsConnection].only_b) == (:injection, 1)
    @test (d12[:PhotoVoltaicUnit].kind, d12[:PhotoVoltaicUnit].only_b) == (:injection, 1)
    n1 = sum(d.common + d.renamed + d.only_a for d in p12.classes if d.kind === :injection)
    @test p12.injection_churn ≈ 3 / (n1 + 2)
    # a new SSH p and q is neither churn nor drift
    d23 = Dict(d.class => d for d in p23.classes)
    @test (d23[:PowerElectronicsConnection].common, d23[:PowerElectronicsConnection].drifted) == (1, 0)
    # the removed lines: common over union of the structural classes; in a
    # node-breaker model switches and connectivity nodes dominate that pool,
    # so a fifth of the lines stays above the default limit of 0.95
    n2 = _series_count(s2, _SERIES_TEST_STRUCTURAL)
    @test p23.structural_overlap ≈ (n2 - length(dropped)) / n2
    @test p23.structural_overlap > 0.95 && p23.split_reason === nothing
    @test [s.snapshots for s in r.series] == [[1, 2, 3]]
    # a limit above that overlap splits the series at the removed lines
    rs = checkCGMESSeries(inputs; config = QstsSeriesConfig(min_structural_overlap = 0.99))
    @test [s.snapshots for s in rs.series] == [[1, 2], [3]]
    @test rs.pairs[1].split_reason === nothing && rs.pairs[2].split_reason === :structural_overlap
    # max_gap 30 minutes splits every pair
    rg = checkCGMESSeries(inputs; config = QstsSeriesConfig(max_gap = 30.0))
    @test [s.snapshots for s in rg.series] == [[1], [2], [3]] && all(p -> p.split_reason === :max_gap, rg.pairs)
    # the churn limit flags the pair with the removed machine, the series stays whole
    rc = checkCGMESSeries(inputs; config = QstsSeriesConfig(max_injection_churn = 0.0))
    @test [s.snapshots for s in rc.series] == [[1, 2, 3]]
    @test rc.pairs[1].churn_flag && !rc.pairs[2].churn_flag
  end)() end

  @testset "CGMES series check: report and CSV" begin (function ()
    r = checkCGMESSeries(inputs; config = QstsSeriesConfig())
    io = IOBuffer()
    printCGMESSeriesReport(io, r)
    text = String(take!(io))
    @test occursin("pair s1 -> s2", text) && occursin("pair s2 -> s3", text)
    mktempdir() do d
      csv = writeCGMESSeriesReport(joinpath(d, "series.csv"), r)
      rows = split(chomp(read(csv, String)), '\n')
      @test length(rows) == 1 + sum(length(p.classes) for p in r.pairs)
      @test Set((split(row, ',')[2], split(row, ',')[3]) for row in rows[2:end]) == Set([("s1", "s2"), ("s2", "s3")])
      # a label with a comma is quoted and keeps the column count
      x = inputs[1]
      rq = checkCGMESSeries([CGMESSeriesInput(x.source; label = "s1, base", scenario_time = x.scenario_time, model_key = x.model_key); inputs[2:end]]; config = QstsSeriesConfig())
      qrows = split(chomp(read(writeCGMESSeriesReport(joinpath(d, "quoted.csv"), rq), String)), '\n')
      first_pair = filter(row -> occursin("\"s1, base\"", row), qrows)
      @test length(first_pair) == length(rq.pairs[1].classes)
      @test all(row -> count(==(','), replace(row, "\"s1, base\"" => "s1")) == count(==(','), qrows[1]), first_pair)
    end
  end)() end

  @testset "CGMES series check: parameter drift" begin (function ()
    # d1: three lines get a new r, a fourth a new x within the tolerance,
    # five breakers open, loads, machines, equivalent injections and
    # generating units get new set points and dispatch; d2: a fifth line gets
    # a new mRID and a new x
    lines = sort!(copy(s1.byclass[:ACLineSegment]))
    d1 = deepcopy(s1)
    for ln in lines[1:3]
      o = d1.objects[ln]
      o.attrs[:r] = string(C.num(o, :r) * 1.05 + 0.01)
    end
    o4 = d1.objects[lines[4]]
    @test C.num(o4, :x) != 0.0
    o4.attrs[:x] = string(C.num(o4, :x) * (1 + 1e-9))
    for sw in sort!(copy(d1.byclass[:Breaker]))[1:5]
      d1.objects[sw].attrs[:open] = "true"
    end
    # every numeric state value of the SmallGrid injections; the set check
    # below fails if the fixture stops carrying one of them
    state = [(:EnergyConsumer, (:p, :q)), (:SynchronousMachine, (:p, :q, :referencePriority)),
             (:EquivalentInjection, (:p, :q, :regulationTarget)), (:ThermalGeneratingUnit, (:normalPF, :initialP))]
    scaled = Set{Tuple{Symbol,Symbol}}()
    for (cls, attrs) in state, o in C.objectsOf(d1, cls), a in attrs
      v = C.num(o, a)
      v === nothing && continue
      o.attrs[a] = string(v * 1.1 + 1.0)
      push!(scaled, (cls, a))
    end
    @test scaled == Set((cls, a) for (cls, attrs) in state for a in attrs)
    d2 = deepcopy(d1)
    _series_rename!(d2, lines[5], "_new" * lines[5])
    o5 = d2.objects["_new"*lines[5]]
    o5.attrs[:x] = string(C.num(o5, :x) * 1.2 + 0.01)
    tq = DateTime(2026, 1, 1, 12, 0, 0)
    dinputs = [CGMESSeriesInput(s1; label = "d0", scenario_time = tq, model_key = "m"),
               CGMESSeriesInput(d1; label = "d1", scenario_time = tq + Minute(15), model_key = "m"),
               CGMESSeriesInput(d2; label = "d2", scenario_time = tq + Minute(30), model_key = "m")]
    r = checkCGMESSeries(dinputs; config = QstsSeriesConfig(max_injection_churn = 0.0))
    q1, q2 = r.pairs
    drift1 = Dict(d.class => d.drifted for d in q1.classes)
    drift2 = Dict(d.class => d.drifted for d in q2.classes)
    @test drift1[:ACLineSegment] == 3 && all(v == 0 for (c, v) in drift1 if c !== :ACLineSegment)
    @test only(d for d in q2.classes if d.class === :ACLineSegment).renamed == 1
    @test drift2[:ACLineSegment] == 1 && all(v == 0 for (c, v) in drift2 if c !== :ACLineSegment)
    # drift is reported only: overlap, churn and the series stay as they were
    @test all(q -> q.structural_overlap == 1.0 && q.injection_churn == 0.0 && !q.churn_flag && q.split_reason === nothing, r.pairs)
    @test [s.snapshots for s in r.series] == [[1, 2, 3]]
    # the text report lists the class whose only difference is the drift
    io = IOBuffer()
    printCGMESSeriesReport(io, r)
    text = split(String(take!(io)), '\n')
    i1 = findfirst(startswith("pair d0 -> d1"), text)
    i2 = findfirst(startswith("pair d1 -> d2"), text)
    block1 = filter(startswith("  "), text[(i1 + 1):(i2 - 1)])
    @test length(block1) == 1 && occursin("ACLineSegment", block1[1]) && occursin(r"\s3$", block1[1])
    mktempdir() do tmp
      rows = [split(row, ',') for row in split(chomp(read(writeCGMESSeriesReport(joinpath(tmp, "drift.csv"), r), String)), '\n')]
      head = rows[1]
      col(name) = findfirst(==(name), head)
      @test head[end] == "drifted" && all(row -> length(row) == length(head), rows)
      got = Dict((row[col("a")], row[col("class")]) => parse(Int, row[col("drifted")]) for row in rows[2:end])
      @test got[("d0", "ACLineSegment")] == 3 && got[("d1", "ACLineSegment")] == 1
      @test count(!=(0), values(got)) == 2
    end
    # a battery behind a power electronics connection: stored energy and
    # battery state are CGMES 3.0 SSH state and no drift, a new rated energy
    # (EQ) in the next step is
    b0 = _series_add_pec!(deepcopy(s1), -12.0, -2.0; unit = :BatteryUnit,
                          unit_attrs = Dict{Symbol,String}(:minP => "-10", :maxP => "10", :ratedE => "40", :storedE => "20", :batteryState => "BatteryStateKind.Discharging"))
    b1 = _series_set_pec!(deepcopy(b0), -14.0, -3.0)
    b1.objects["_pec_series_unit"].attrs[:storedE] = "17"
    b1.objects["_pec_series_unit"].attrs[:batteryState] = "BatteryStateKind.Charging"
    b2 = deepcopy(b1)
    b2.objects["_pec_series_unit"].attrs[:ratedE] = "50"
    rb = checkCGMESSeries([CGMESSeriesInput(st; label = l, scenario_time = tq + Minute(15 * (k - 1)), model_key = "m") for (k, (st, l)) in enumerate(zip((b0, b1, b2), ("b0", "b1", "b2")))]; config = QstsSeriesConfig())
    db = [Dict(d.class => d for d in q.classes) for q in rb.pairs]
    @test [(d[:BatteryUnit].common, d[:BatteryUnit].drifted) for d in db] == [(1, 0), (1, 1)]
    @test all(d -> (d[:PowerElectronicsConnection].common, d[:PowerElectronicsConnection].drifted) == (1, 0), db)
  end)() end

  @testset "CGMES series check: header defaults of a store and a path" begin (function ()
    # without keywords both take the authority set of the EQ header and the
    # scenario time of the SSH header
    r = checkCGMESSeries([CGMESSeriesInput(s1), CGMESSeriesInput(files; label = "path")]; config = QstsSeriesConfig())
    eq = only(f for f in s1.files if :EQ in f.profiles)
    a, b = r.snapshots
    @test !isempty(eq.modeling_authority) && a.model_key == eq.modeling_authority && b.model_key == eq.modeling_authority
    @test a.scenario_time == DateTime(2030, 1, 2, 9, 0, 0) && b.scenario_time == a.scenario_time
    @test a.label == "snapshot" && b.label == "path"
    @test a.counts == b.counts && a.has_tp && b.has_tp
    @test length(r.series) == 1 && only(r.pairs).structural_overlap == 1.0 && only(r.pairs).injection_churn == 0.0
  end)() end

  # a switch closed in the SSH while the TP keeps its two sides on two TNs
  # (case 3): counted per snapshot and named in the report
  @testset "CGMES series check: SSH against TP" begin (function ()
    local t2 = deepcopy(s1)
    local opened = _union_open_line!(t2)[2]
    local t3 = deepcopy(t2)
    t3.objects[opened[1]].attrs[:open] = "false"
    local tz = DateTime(2026, 1, 1, 12)
    local ins = [CGMESSeriesInput(st; label = l, scenario_time = tz + Minute(60 * (k - 1)), model_key = "m") for (k, (st, l)) in enumerate(zip((s1, t2, t3), ("s1", "split", "closed")))]
    local rc = checkCGMESSeries(ins; config = QstsSeriesConfig())
    @test [sn.ssh_tp_conflicts for sn in rc.snapshots] == [0, 0, 1]
    local buf = IOBuffer()
    printCGMESSeriesReport(buf, rc)
    @test occursin("snapshot closed: 1 switch(es) closed in the SSH, not retained, between two topological nodes (SSH against TP, case 3)", String(take!(buf)))
  end)() end
end

# union network of a CGMES snapshot series (task stage 3, items 17/18) on
# snapshots derived in memory from the ENTSO-E SmallGrid fixture with one
# power electronics connection added: s2 opens the switches at both ends of
# one line, so its two line-side CNs move to new TNs, scales the loads and
# gives the connection another p and q; s3 drops one machine and scales the
# loads the other way; every step must be reproduced by template plus delta

# open the first switch at each end of the first line (by mRID) that has a
# switch at both of its CNs (`ends = 1:1` opens the first end alone); inside
# each old TN the CNs reachable from the line end over the remaining closed
# switches move to a new TN
function _union_open_line!(s; ends = 1:2)
  C = Sparlectra.CGMESImporter
  by_eq = Dict{String,Vector{Any}}(); by_cn = Dict{String,Vector{Any}}()
  for t in C.objectsOf(s, :Terminal)
    eq = get(t.refs, :ConductingEquipment, nothing); cn = get(t.refs, :ConnectivityNode, nothing)
    eq === nothing || push!(get!(Vector{Any}, by_eq, eq), t)
    cn === nothing || push!(get!(Vector{Any}, by_cn, cn), t)
  end
  sw_classes = Set(C._SWITCH_CLASSES)
  tn_of(t) = something(get(t.refs, :TopologicalNode, nothing), get(s.objects[t.refs[:ConnectivityNode]].refs, :TopologicalNode, nothing), Some(nothing))
  switches_at(cn) = sort!([t.refs[:ConductingEquipment] for t in get(by_cn, cn, []) if s.objects[t.refs[:ConductingEquipment]].class in sw_classes])
  for line in sort!(copy(s.byclass[:ACLineSegment]))
    ts = get(by_eq, line, [])
    length(ts) == 2 || continue
    sws = [switches_at(t.refs[:ConnectivityNode]) for t in ts]
    all(!isempty, sws) || continue
    opened = [first(sws[i]) for i in ends]
    for m in opened
      s.objects[m].attrs[:open] = "true"
    end
    for (i, t) in enumerate(ts)
      i in ends || continue
      old = tn_of(t)
      seen = Set([t.refs[:ConnectivityNode]]); stack = [t.refs[:ConnectivityNode]]
      while !isempty(stack)
        cn = pop!(stack)
        for tt in get(by_cn, cn, [])
          eq = s.objects[tt.refs[:ConductingEquipment]]
          (eq.class in sw_classes && !(eq.mrid in opened) && get(eq.attrs, :open, "false") != "true") || continue
          for o in by_eq[eq.mrid]
            c2 = o.refs[:ConnectivityNode]
            (c2 in seen || tn_of(o) != old) && continue
            push!(seen, c2); push!(stack, c2)
          end
        end
      end
      new = string("_s2TN", i)
      tnobj = s.objects[old]
      s.objects[new] = C.CIMObject(new, :TopologicalNode, copy(tnobj.attrs), copy(tnobj.refs), tnobj.source)
      push!(s.byclass[:TopologicalNode], new)
      for cn in seen
        haskey(s.objects[cn].refs, :TopologicalNode) && (s.objects[cn].refs[:TopologicalNode] = new)
        for tt in get(by_cn, cn, [])
          haskey(tt.refs, :TopologicalNode) && (tt.refs[:TopologicalNode] = new)
        end
      end
    end
    return line, opened
  end
  error("SmallGrid: no line with a switch at both ends")
end

function _union_scale_loads!(s, f)
  C = Sparlectra.CGMESImporter
  n = 0
  for cls in (:EnergyConsumer, :ConformLoad, :NonConformLoad), o in C.objectsOf(s, cls)
    for a in (:p, :q)
      v = C.num(o, a)
      v === nothing && continue
      o.attrs[a] = string(v * f)
      # a zero stays zero and gives no op
      v == 0.0 || (n += 1)
    end
  end
  return n
end

# drop a machine with its terminals and the SV flows on them
function _union_drop_machine!(s, m)
  C = Sparlectra.CGMESImporter
  terms = Set(t.mrid for t in C.objectsOf(s, :Terminal) if get(t.refs, :ConductingEquipment, "") == m)
  for (cls, ids) in s.byclass
    filter!(ids) do x
      o = s.objects[x]
      drop = x == m || x in terms || (cls === :SvPowerFlow && get(o.refs, :Terminal, "") in terms)
      drop && delete!(s.objects, x)
      !drop
    end
  end
  return s
end

# put the terminal of the first load (by mRID) on the connectivity node of
# the next load on another TN, as an export that hangs an equivalent on
# another node: against an unchanged snapshot the terminal becomes a union
# node of its own, with a _TPL coupling to its old node and one to its new
function _union_move_load!(s)
  C = Sparlectra.CGMESImporter
  tn_of(t) = something(get(t.refs, :TopologicalNode, nothing), get(s.objects[t.refs[:ConnectivityNode]].refs, :TopologicalNode, nothing), Some(nothing))
  term = Dict(t.refs[:ConductingEquipment] => t for t in C.objectsOf(s, :Terminal) if haskey(t.refs, :ConductingEquipment) && haskey(t.refs, :ConnectivityNode))
  loads = sort!([m for m in s.byclass[:EnergyConsumer] if haskey(term, m)])
  a = term[loads[1]]
  b = term[first(m for m in loads[2:end] if tn_of(term[m]) != tn_of(a))]
  a.refs[:ConnectivityNode] = b.refs[:ConnectivityNode]
  haskey(a.refs, :TopologicalNode) && (a.refs[:TopologicalNode] = tn_of(b))
  return s
end

# J1 - J2 = JX = J3 - J4: JX has no branch and joins J2 and J3 through two
# closed links. The reference sits at J1, J4 carries a machine without
# regulation, loads at J2 and J3: one island of the power flow, split at JX
# by the branch-only isolation marking of a re-import.
function _union_switch_only_island_net()::Net
  net = Net(name = "switch_only_island", baseMVA = 100.0)
  for b in ("J1", "J2", "JX", "J3", "J4")
    addBus!(net = net, busName = b, vn_kV = 110.0)
  end
  addPIModelACLine!(net = net, fromBus = "J1", toBus = "J2", r_pu = 0.02, x_pu = 0.2, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "J3", toBus = "J4", r_pu = 0.02, x_pu = 0.2, b_pu = 0.0, status = 1)
  addLink!(net = net, fromBus = "J2", toBus = "JX", status = 1)
  addLink!(net = net, fromBus = "JX", toBus = "J3", status = 1)
  addProsumer!(net = net, busName = "J1", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.02, va_deg = 0.0, referencePri = "J1")
  addProsumer!(net = net, busName = "J4", type = "SYNCHRONOUSMACHINE", p = 20.0, q = 0.0)
  addProsumer!(net = net, busName = "J2", type = "ENERGYCONSUMER", p = 12.0, q = 3.0)
  addProsumer!(net = net, busName = "J3", type = "ENERGYCONSUMER", p = 15.0, q = 4.0)
  ok, msg = validate!(net = net)
  ok || error("switch-only island net invalid: $msg")
  return net
end

function run_cgmes_union_tests()
  @testset "CGMES series union" begin (function ()
    C = Sparlectra.CGMESImporter
    dir = joinpath(@__DIR__, "fixtures", "cgmes_smallgrid")
    files = sort!(filter(f -> endswith(f, ".zip"), readdir(dir; join = true)))
    s1 = redirect_stdout(() -> C.loadCGMES(files), devnull)
    _series_add_pec!(s1, -12.0, -2.0)
    cfg = Sparlectra.SparlectraConfig(Dict{String,Any}("cgmes_import" => Dict{String,Any}("map_power_electronics" => true)))
    r1 = redirect_stdout(() -> importCGMES(cfg; path = s1, name = "s1"), devnull)
    s2 = deepcopy(s1)
    line, opened = _union_open_line!(s2)
    nload = _union_scale_loads!(s2, 1.1)
    _series_set_pec!(s2, -18.0, -4.0)
    s3 = deepcopy(s1)
    _union_scale_loads!(s3, 0.9)
    # the first machine by mRID that is not the slack of the s1 import
    machine = first(m for m in sort!(copy(s3.byclass[:SynchronousMachine])) if r1.element_by_mrid[m][1] === :generator &&
                    r1.net.prosumpsVec[r1.element_by_mrid[m][2]].referencePri === nothing)
    _union_drop_machine!(s3, machine)
    t0 = DateTime(2026, 1, 1, 12)
    inputs = [CGMESSeriesInput(s1; label = "s1", scenario_time = t0, model_key = "m"),
              CGMESSeriesInput(s2; label = "s2", scenario_time = t0 + Minute(60), model_key = "m"),
              CGMESSeriesInput(s3; label = "s3", scenario_time = t0 + Minute(120), model_key = "m")]
    u = redirect_stdout(() -> buildCGMESSeriesUnion(inputs; config = cfg), devnull)
    ntn = length(C._tpTNs(s1)); nbd = count(m -> m in s1.boundary, s1.byclass[:TopologicalNode])
    # two line-side CNs split off in s2: two union nodes more than the base TNs
    @test length(u.node_mrid) == ntn + nbd + 2
    @test sort(u.link_switches) == sort(opened)
    @test all(u.link_status[k][m] == (k == 2 ? 0 : 1) for m in opened, k in 1:3)
    @test u.n_tp_links == 0
    @test isempty(u.deltas[1].ops) && isempty(u.deltas[1].residual)
    d2, d3 = u.deltas[2], u.deltas[3]
    @test count(op -> op.target === :link, d2.ops) == 2
    @test !any(op -> op.target === :link, d3.ops)
    # link changes by kind: the two opened switches, no _TPL coupling
    @test u.n_switch_link_changes == [0, 2, 0] && u.n_tpl_link_changes == [0, 0, 0]
    # s2 with one load terminal on the node of another load: two _TPL
    # couplings, both change in step 2 next to the two switch links
    s2m = _union_move_load!(deepcopy(s2))
    um = redirect_stdout(() -> buildCGMESSeriesUnion([inputs[1], CGMESSeriesInput(s2m; label = "s2m", scenario_time = t0 + Minute(60), model_key = "m")]; config = cfg), devnull)
    @test um.n_tp_links == 2 && isempty(um.messages)
    @test um.n_switch_link_changes == [0, 2] && um.n_tpl_link_changes == [0, 2]
    @test count(op -> op.target === :link, um.deltas[2].ops) == 4
    @test isempty(compareSnapshotReproduction(um, 2, redirect_stdout(() -> importCGMES(cfg; path = s2m, name = "s2m"), devnull)))
    @test count(op -> op.field === :in_service && op.value == 0.0, d3.ops) == 1
    @test count(op -> op.target === :load && op.field in (:p, :q), d2.ops) == nload
    # the mapped power electronics connection: a p and a q op in s2 (the
    # negated SSH values), none in s3, which has the values of s1
    pec_id = u.index.id_by_mrid["_pec_series"]
    @test u.index.kind_by_id[pec_id] === :generator
    @test [(op.op, op.field, op.value) for op in d2.ops if op.id == pec_id] == [(:set, :p, 18.0), (:set, :q, 4.0)]
    @test !any(op -> op.id == pec_id, d3.ops)
    # every step is reproduced by template plus delta, the model and the
    # solved state: every union bus that R_k builds carries the voltage of
    # its R_k bus
    rks = [redirect_stdout(() -> importCGMES(cfg; path = s, name = string("s", k)), devnull) for (k, s) in enumerate((s1, s2, s3))]
    for k in 1:3
      @test isempty(compareSnapshotReproduction(u, k, rks[k]))
      net = deepcopy(u.template)
      Sparlectra.apply!(net, vcat(u.deltas[k].ops, u.deltas[k].control_ops), u.index)
      @test last(runpf!(net, 30, 1e-10, 0)) == 0
      @test last(runpf!(rks[k].net, 30, 1e-10, 0)) == 0
      rk_bus = Dict(tn => rks[k].net.busDict[name] for (tn, name) in rks[k].topo.bus_name if haskey(rks[k].net.busDict, name))
      dv = 0.0; da = 0.0; nb = 0
      for (b, nid) in u.node_of_bus
        rb = get(rk_bus, u.node_sig[nid][k], 0)
        (rb == 0 || Sparlectra.isIsolated(net.nodeVec[b])) && continue
        dv = max(dv, abs(net.nodeVec[b]._vm_pu - rks[k].net.nodeVec[rb]._vm_pu))
        da = max(da, abs(net.nodeVec[b]._va_deg - rks[k].net.nodeVec[rb]._va_deg))
        nb += 1
      end
      @test nb > 0 && dv <= 1e-8 && da <= 1e-6
    end
    # the comparison sees a wrong snapshot
    @test !isempty(compareSnapshotReproduction(u, 2, rks[3]))
    # mode A on the union: every step converges, every union bus carries the
    # voltage of the solved direct import, the SV metrics equal compareWithSV
    # of the direct import, and the line switched off in s2 has no rows there
    # the run takes the cgmes_import settings the union was built with
    qcfg = Sparlectra.SparlectraConfig(Dict{String,Any}("qsts" => Dict{String,Any}("mode" => "snapshots"), "power_flow" => Dict{String,Any}("tol" => 1e-10), "cgmes_import" => Dict{String,Any}("map_power_electronics" => true)))
    r = redirect_stdout(() -> runQsts(u; config = qcfg), devnull)
    @test r.base.status === :solved && all(s -> s.status === :solved, r.steps)
    @test [s.is_base for s in r.steps] == [true, false, false] && r.origin === :cgmes_series && r.step_minutes == 60.0
    step_value(q, k) = r.step_values[findfirst(==(q), QSTS_STEP_QUANTITIES), k]
    for k in 1:3
      rk_bus = Dict(tn => rks[k].net.busDict[name] for (tn, name) in rks[k].topo.bus_name if haskey(rks[k].net.busDict, name))
      dv = 0.0; da = 0.0; nb = 0
      for (b, nid) in u.node_of_bus
        rb = get(rk_bus, u.node_sig[nid][k], 0)
        (rb == 0 || isnan(r.bus_values[b, 1, k])) && continue
        dv = max(dv, abs(r.bus_values[b, 1, k] - rks[k].net.nodeVec[rb]._vm_pu))
        da = max(da, abs(r.bus_values[b, 2, k] - rks[k].net.nodeVec[rb]._va_deg))
        nb += 1
      end
      @test nb > 0 && dv <= 1e-8 && da <= 1e-6
      c = compareWithSV(rks[k])
      @test step_value("sv_n_nodes", k) == c.n
      @test all(isapprox(step_value(q, k), getfield(c, f); atol = 1e-8) for (q, f) in (("sv_max_dvm_pu", :max_dvm), ("sv_rms_dvm_pu", :rms_dvm), ("sv_max_dva_deg", :max_dva),
                                                                                      ("sv_rms_dva_deg", :rms_dva), ("sv_va_offset_deg", :va_ref_offset_deg)))
      @test step_value("sv_tn_spread_pu", k) <= 1e-12
    end
    jl = findfirst(==(line), r.elements.branch_ids)
    @test !isnan(r.branch_values[jl, 1, 1]) && all(isnan, r.branch_values[jl, :, 2])
    lt = qstsLongTable(r)
    @test any(i -> lt.step[i] == 1 && lt.element_id[i] == line, eachindex(lt.step)) && !any(i -> lt.step[i] == 2 && lt.element_id[i] == line, eachindex(lt.step))
    # the switch and the import configuration of the union, each named
    @test_throws "qsts.mode" runQsts(u; config = cfg)
    other = Sparlectra.SparlectraConfig(Dict{String,Any}("qsts" => Dict{String,Any}("mode" => "snapshots"), "cgmes_import" => Dict{String,Any}("map_energy_source" => true)))
    @test_throws "cgmes_import.map_energy_source" runQsts(u; config = other)
    # the vocabulary is refused for a series with machine controllers
    mc = Sparlectra.SparlectraConfig(Dict{String,Any}("cgmes_import" => Dict{String,Any}("machine_control" => true)))
    @test_throws C.CGMESSeriesError buildCGMESSeriesUnion(inputs; config = mc)
    # the base step defaults to qsts.base_step of the configuration; a step
    # outside the series is refused before any snapshot is read
    beyond = Sparlectra.SparlectraConfig(Dict{String,Any}("qsts" => Dict{String,Any}("base_step" => 4)))
    @test_throws C.CGMESSeriesError buildCGMESSeriesUnion(inputs; config = beyond)
    @test_throws "base_step 4 outside 1:3" buildCGMESSeriesUnion(inputs; config = beyond)
  end)() end

  # the SmallGrid SV carries no flow at a line terminal, so the snapshots get
  # one SvPowerFlow per line terminal with the model flow of their solved
  # direct import; mode A on the union of these snapshots must give the same
  # flows, which checks terminal, branch position and side of the reference
  @testset "CGMES series union: SV flows" begin (function ()
    C = Sparlectra.CGMESImporter
    dir = joinpath(@__DIR__, "fixtures", "cgmes_smallgrid")
    s1 = redirect_stdout(() -> C.loadCGMES(sort!(filter(f -> endswith(f, ".zip"), readdir(dir; join = true)))), devnull)
    s2 = deepcopy(s1)
    line, _ = _union_open_line!(s2)
    _union_scale_loads!(s2, 1.1)
    cfg = Sparlectra.SparlectraConfig(Dict())
    terminals = 0
    for (k, s) in enumerate((s1, s2))
      rk = redirect_stdout(() -> importCGMES(cfg; path = s, name = string("s", k)), devnull)
      @test last(redirect_stdout(() -> runpf!(rk.net, 30, 1e-10, 0), devnull)) == 0
      V = ComplexF64[nd._vm_pu * cis(deg2rad(nd._va_deg)) for nd in rk.net.nodeVec]
      branch_of = Dict(br.branchIdx => br for br in rk.net.branchVec)
      source = s.objects[first(s.byclass[:SvVoltage])].source
      for t in C.objectsOf(s, :Terminal)
        get(t.refs, :ConductingEquipment, "") in s.byclass[:ACLineSegment] || continue
        haskey(rk.branch_side_of_terminal, t.mrid) || continue
        idx, side = rk.branch_side_of_terminal[t.mrid]
        flow = Sparlectra._sv_branch_terminal_flow(branch_of[idx], V, side, rk.net.baseMVA)
        m = string("_svpf_", lstrip(t.mrid, '_'))
        s.objects[m] = C.CIMObject(m, :SvPowerFlow, Dict{Symbol,String}(:p => string(real(flow)), :q => string(imag(flow))), Dict{Symbol,String}(:Terminal => t.mrid), source)
        push!(get!(Vector{String}, s.byclass, :SvPowerFlow), m)
        terminals += 1
      end
    end
    t0 = DateTime(2026, 1, 1, 12)
    inputs = [CGMESSeriesInput(s1; label = "s1", scenario_time = t0, model_key = "m"), CGMESSeriesInput(s2; label = "s2", scenario_time = t0 + Minute(60), model_key = "m")]
    u = redirect_stdout(() -> buildCGMESSeriesUnion(inputs; config = cfg), devnull)
    # the terminals of the line switched off in s2 sit on dead buses there
    @test [length(ref.branch) for ref in u.sv] == [terminals ÷ 2, terminals ÷ 2 - 2]
    qcfg = Sparlectra.SparlectraConfig(Dict{String,Any}("qsts" => Dict{String,Any}("mode" => "snapshots"), "power_flow" => Dict{String,Any}("tol" => 1e-10)))
    r = redirect_stdout(() -> runQsts(u; config = qcfg), devnull)
    @test all(s -> s.status === :solved, r.steps)
    step_value(q, k) = r.step_values[findfirst(==(q), QSTS_STEP_QUANTITIES), k]
    @test [step_value("sv_n_terminals", k) for k in 1:2] == [terminals ÷ 2, terminals ÷ 2 - 2]
    @test all(step_value(q, k) <= 1e-6 for q in ("sv_max_dp_mw", "sv_max_dq_mvar"), k in 1:2)
    # both sides of a line carry their own difference
    cols = [findfirst(==(q), QSTS_BRANCH_QUANTITIES) for q in ("sv_dp_from_mw", "sv_dq_from_mvar", "sv_dp_to_mw", "sv_dq_to_mvar")]
    jl = findfirst(==(line), r.elements.branch_ids)
    @test all(!isnan, r.branch_values[jl, cols, 1]) && all(isnan, r.branch_values[jl, cols, 2])
  end)() end

  # the star point of a three-winding transformer has no TN and no mRID; the
  # union addresses it as `<PowerTransformer mRID>#star` and isolates it in a
  # step without the transformer (found on a real series: without the
  # address every leg showed as misplaced)
  @testset "CGMES series union: star point of a three-winding transformer" begin (function ()
    C = Sparlectra.CGMESImporter
    s1 = redirect_stdout(() -> C.loadCGMES(_cgmes_synth3w_dir()), devnull)
    s2 = deepcopy(s1)
    ends = Set(m for m in s2.byclass[:PowerTransformerEnd] if get(s2.objects[m].refs, :PowerTransformer, "") == "_t3w")
    drop = union(Set(["_t3w", "_rtc2", "_tcc2"]), ends, Set(m for m in s2.byclass[:Terminal] if get(s2.objects[m].refs, :ConductingEquipment, "") == "_t3w"))
    for (cls, ids) in s2.byclass
      filter!(x -> !(x in drop), ids)
    end
    foreach(m -> delete!(s2.objects, m), drop)
    cfg = Sparlectra.SparlectraConfig(Dict())
    t0 = DateTime(2026, 1, 1, 12)
    inputs = [CGMESSeriesInput(s1; label = "with", scenario_time = t0, model_key = "m"),
              CGMESSeriesInput(s2; label = "without", scenario_time = t0 + Minute(60), model_key = "m")]
    u = redirect_stdout(() -> buildCGMESSeriesUnion(inputs; config = cfg), devnull)
    star = get(u.index.id_by_mrid, "_t3w#star", 0)
    @test star != 0 && u.index.kind_by_id[star] === :bus
    @test any(op -> op.id == star && op.field === :isolated && op.value == 1.0, u.deltas[2].ops)
    @test count(op -> op.target === :transformer && op.field === :status && op.value == 0.0, u.deltas[2].ops) == 3
    for (k, s) in enumerate((s1, s2))
      rk = redirect_stdout(() -> importCGMES(cfg; path = s, name = string("s", k)), devnull)
      @test isempty(compareSnapshotReproduction(u, k, rk))
    end
  end)() end

  # an exporter may reuse a mRID for an object of another class (found on a
  # real series: a GeneratingUnit came back as an EquivalentInjection); the
  # union keeps both, the second under `<mRID>#<class>`
  @testset "CGMES series union: a mRID reused for another class" begin (function ()
    C = Sparlectra.CGMESImporter
    dir = joinpath(@__DIR__, "fixtures", "cgmes_smallgrid")
    s1 = redirect_stdout(() -> C.loadCGMES(sort!(filter(f -> endswith(f, ".zip"), readdir(dir; join = true)))), devnull)
    cfg = Sparlectra.SparlectraConfig(Dict())
    r1 = redirect_stdout(() -> importCGMES(cfg; path = s1, name = "s1"), devnull)
    # a machine that is not the reference, with its generating unit
    sm = first(m for m in sort!(copy(s1.byclass[:SynchronousMachine])) if r1.element_by_mrid[m][1] === :generator &&
               r1.net.prosumpsVec[r1.element_by_mrid[m][2]].referencePri === nothing && haskey(s1.objects[m].refs, :GeneratingUnit))
    g = s1.objects[sm].refs[:GeneratingUnit]
    s2 = deepcopy(s1)
    t = first(o for o in C.objectsOf(s2, :Terminal) if get(o.refs, :ConductingEquipment, "") == sm)
    # s2: the machine and its unit are gone, an equivalent injection with the
    # unit's mRID feeds the same power at the same node (load sign in both)
    _union_drop_machine!(s2, sm)
    filter!(!=(g), s2.byclass[s2.objects[g].class]); delete!(s2.objects, g)
    o = s1.objects[sm]
    s2.objects[g] = C.CIMObject(g, :EquivalentInjection, Dict{Symbol,String}(:name => "EI", :p => get(o.attrs, :p, "0"), :q => get(o.attrs, :q, "0")), Dict{Symbol,String}(), o.source)
    push!(get!(Vector{String}, s2.byclass, :EquivalentInjection), g)
    tm = string(g, "_T1")
    refs = Dict{Symbol,String}(:ConductingEquipment => g)
    for key in (:ConnectivityNode, :TopologicalNode)
      haskey(t.refs, key) && (refs[key] = t.refs[key])
    end
    s2.objects[tm] = C.CIMObject(tm, :Terminal, Dict{Symbol,String}(:sequenceNumber => "1", :connected => "true"), refs, t.source)
    push!(s2.byclass[:Terminal], tm)
    t0 = DateTime(2026, 1, 1, 12)
    inputs = [CGMESSeriesInput(s1; label = "machine", scenario_time = t0, model_key = "m"),
              CGMESSeriesInput(s2; label = "equivalent", scenario_time = t0 + Minute(60), model_key = "m")]
    u = redirect_stdout(() -> buildCGMESSeriesUnion(inputs; config = cfg), devnull)
    ei = string(g, "#EquivalentInjection")
    @test get(u.alias, ei, "") == g
    @test u.index.kind_by_id[u.index.id_by_mrid[ei]] === :external_grid
    d2 = u.deltas[2]
    @test any(op -> op.id == u.index.id_by_mrid[sm] && op.field === :in_service && op.value == 0.0, d2.ops)
    @test any(op -> op.id == u.index.id_by_mrid[ei] && op.field === :in_service && op.value == 1.0, d2.ops)
    @test isempty(d2.residual)
    rks = [redirect_stdout(() -> importCGMES(cfg; path = s, name = string("s", k)), devnull) for (k, s) in enumerate((s1, s2))]
    for k in 1:2
      @test isempty(compareSnapshotReproduction(u, k, rks[k]))
    end
    # the same reuse between two steps that are not the base: the base step
    # has neither machine nor unit, both steps add objects, and in time order
    # the unit of step 1 keeps its mRID and the injection of step 3 takes the
    # alias
    s0 = deepcopy(s1)
    _union_drop_machine!(s0, sm)
    filter!(!=(g), s0.byclass[s0.objects[g].class]); delete!(s0.objects, g)
    inputs3 = [CGMESSeriesInput(s; label = l, scenario_time = t0 + Minute(60 * (k - 1)), model_key = "m")
               for (k, (s, l)) in enumerate(((s1, "machine"), (s0, "neither"), (s2, "equivalent")))]
    u3 = redirect_stdout(() -> buildCGMESSeriesUnion(inputs3; config = cfg, base_step = 2), devnull)
    @test u3.alias == Dict(ei => g)
    sm3 = get(u3.index.id_by_mrid, sm, 0); ei3 = get(u3.index.id_by_mrid, ei, 0)
    @test ei3 != 0 && u3.index.kind_by_id[ei3] === :external_grid
    # both are missing from the base step, so the template parks both
    @test sm3 != 0 && sm3 in u3.parked && ei3 in u3.parked
    r0 = redirect_stdout(() -> importCGMES(cfg; path = s0, name = "s0"), devnull)
    for (k, rk) in enumerate((rks[1], r0, rks[2]))
      @test isempty(compareSnapshotReproduction(u3, k, rk))
    end
  end)() end

  # in an island without a slack the solver promotes a unit to reference and
  # keeps its start voltage, so there the start is part of the state and the
  # delta carries it (found on a real series: two-bus islands kept the start
  # of the template instead of the snapshot's). An island the import gives a
  # slack holds the start angle of that slack bus, which the delta carries too
  @testset "CGMES series union: start values of an island without slack" begin (function ()
    C = Sparlectra.CGMESImporter
    dir = joinpath(@__DIR__, "fixtures", "cgmes_smallgrid")
    s1 = redirect_stdout(() -> C.loadCGMES(sort!(filter(f -> endswith(f, ".zip"), readdir(dir; join = true)))), devnull)
    cfg = Sparlectra.SparlectraConfig(Dict())
    r1 = redirect_stdout(() -> importCGMES(cfg; path = s1, name = "s1"), devnull)
    # the first line (by mRID) whose outage leaves a generator in an island
    # without a slack; a regulating one fixes the magnitude, the start the angle
    local line = ""
    local island = Int[]
    for m in sort!([m for (m, (k, _)) in r1.element_by_mrid if k === :branch])
      n = deepcopy(r1.net)
      br = n.branchVec[r1.element_by_mrid[m][2]]
      br.status = 0; br.from_status = 0; br.to_status = 0
      for c in Sparlectra.electricalIslandComponents(n)
        # a lone bus without a branch is isolated by the import, the measured
        # islands had two buses
        (length(c) >= 2 && !any(b -> Sparlectra.getNodeType(n.nodeVec[b]) == Sparlectra.Slack, c)) || continue
        any(ps -> Sparlectra.isGenerator(ps) && Sparlectra._scf_prosumer_bus(ps) in c, n.prosumpsVec) || continue
        # some outages leave a main island that does not solve; runpf!
        # reports that as an ErrorException, and such a line is skipped
        solves = try
          redirect_stdout(() -> runpf!(deepcopy(n), 50, 1e-8, 0; islands_enabled = true), devnull)
          true
        catch err
          err isa ErrorException || rethrow()
          false
        end
        solves || continue
        line = m; island = c
        break
      end
      isempty(line) || break
    end
    @test !isempty(line)
    # s2: the line is gone and the SV start of the island is shifted
    s2 = deepcopy(s1)
    terms = Set(t.mrid for t in C.objectsOf(s2, :Terminal) if get(t.refs, :ConductingEquipment, "") == line)
    for (cls, ids) in s2.byclass
      filter!(ids) do x
        drop = x == line || x in terms || (cls === :SvPowerFlow && get(s2.objects[x].refs, :Terminal, "") in terms)
        drop && delete!(s2.objects, x)
        !drop
      end
    end
    tn_of_bus = Dict(r1.net.busDict[name] => tn for (tn, name) in r1.topo.bus_name if haskey(r1.net.busDict, name))
    island_tns = Set(tn_of_bus[b] for b in island)
    for sv in C.objectsOf(s2, :SvVoltage)
      get(sv.refs, :TopologicalNode, "") in island_tns || continue
      sv.attrs[:v] = string(C.num(sv, :v) * 1.02)
      sv.attrs[:angle] = string(C.num(sv, :angle, 0.0) + 3.0)
    end
    # s3: as s2, and the machines of the island lose their voltage control.
    # The import then gives the island a slack (one warning); that slack
    # holds its start angle, so the delta carries the start of its bus
    s3 = deepcopy(s2)
    nreg = 0
    for t in C.objectsOf(s3, :Terminal)
      get(t.refs, :TopologicalNode, "") in island_tns || continue
      eq = get(s3.objects, get(t.refs, :ConductingEquipment, ""), nothing)
      (eq === nothing || !haskey(eq.refs, :RegulatingControl)) && continue
      delete!(eq.refs, :RegulatingControl)
      nreg += 1
    end
    @test nreg > 0
    t0 = DateTime(2026, 1, 1, 12)
    for (sk, rule_slack) in ((s2, false), (s3, true))
      inputs = [CGMESSeriesInput(s1; label = "meshed", scenario_time = t0, model_key = "m"),
                CGMESSeriesInput(sk; label = "island", scenario_time = t0 + Minute(60), model_key = "m")]
      # the import of the rule-slack snapshot logs its warning: once in the
      # union, once in the direct import
      quiet(f) = rule_slack ? (@test_logs (:warn, r"^CGMES import: 1 data defect") redirect_stdout(f, devnull)) : redirect_stdout(f, devnull)
      u = quiet(() -> buildCGMESSeriesUnion(inputs; config = cfg))
      island_ids = Set(u.index.id_by_mrid[tn] for tn in island_tns if haskey(u.index.id_by_mrid, tn))
      @test !isempty(island_ids)
      @test !any(op -> op.id in island_ids && op.field === :isolated, u.deltas[2].ops)
      r2 = quiet(() -> importCGMES(cfg; path = sk, name = "s2"))
      rk_bus = Dict(tn => r2.net.busDict[name] for (tn, name) in r2.topo.bus_name if haskey(r2.net.busDict, name))
      if rule_slack
        # the island's slack bus carries its start; the PQ bus of the island
        # needs none
        slack_tns = [tn for tn in island_tns if Sparlectra.getNodeType(r2.net.nodeVec[rk_bus[tn]]) == Sparlectra.Slack]
        @test length(slack_tns) == 1
        @test count(m -> startswith(m, "warning: electrical island of 2 bus(es) "), r2.messages) == 1
        @test all(any(op -> op.id == u.index.id_by_mrid[tn] && op.field === :va_start, u.deltas[2].ops) for tn in slack_tns)
      else
        @test all(any(op -> op.id == i && op.field === :va_start, u.deltas[2].ops) for i in island_ids)
      end
      @test isempty(compareSnapshotReproduction(u, 2, r2))
      # the solved state, island included
      net = deepcopy(u.template)
      Sparlectra.apply!(net, vcat(u.deltas[2].ops, u.deltas[2].control_ops), u.index)
      @test last(redirect_stdout(() -> runpf!(net, 50, 1e-8, 0; islands_enabled = true), devnull)) == 0
      @test last(redirect_stdout(() -> runpf!(r2.net, 50, 1e-8, 0; islands_enabled = true), devnull)) == 0
      dv = 0.0; da = 0.0
      for (b, nid) in u.node_of_bus
        rb = get(rk_bus, u.node_sig[nid][2], 0)
        (rb == 0 || Sparlectra.isIsolated(net.nodeVec[b])) && continue
        dv = max(dv, abs(net.nodeVec[b]._vm_pu - r2.net.nodeVec[rb]._vm_pu))
        da = max(da, abs(net.nodeVec[b]._va_deg - r2.net.nodeVec[rb]._va_deg))
      end
      @test dv <= 1e-8 && da <= 1e-6
    end
  end)() end

  # QSTS mode A starts each CGMES step from the state of its own snapshot,
  # not from the base solution. s2 has the loads 10 %
  # up and as SV its own state solved under the run settings, so its step
  # starts at its solution. Q limits are off: a PV bus at its limit would
  # cost the switching iterations from any start
  @testset "CGMES series union: mode A starts each step from its own SV" begin (function ()
    C = Sparlectra.CGMESImporter
    dir = joinpath(@__DIR__, "fixtures", "cgmes_smallgrid")
    s1 = redirect_stdout(() -> C.loadCGMES(sort!(filter(f -> endswith(f, ".zip"), readdir(dir; join = true)))), devnull)
    cfg = Sparlectra.SparlectraConfig(Dict())
    qcfg = Sparlectra.SparlectraConfig(Dict{String,Any}("qsts" => Dict{String,Any}("mode" => "snapshots"), "power_flow" => Dict{String,Any}("qlimits" => Dict{String,Any}("enabled" => false))))
    s2 = deepcopy(s1)
    _union_scale_loads!(s2, 1.1)
    solved = redirect_stdout(() -> importCGMES(cfg; path = s2, name = "s2"), devnull)
    ex, _ = redirect_stdout(() -> Sparlectra._qsts_solve!(solved.net, qcfg), devnull)
    @test ex !== nothing && ex.erg == 0
    for sv in C.objectsOf(s2, :SvVoltage)
      b = get(solved.net.busDict, get(solved.topo.bus_name, get(sv.refs, :TopologicalNode, ""), ""), 0)
      b == 0 && continue
      nd = solved.net.nodeVec[b]
      sv.attrs[:v] = string(nd._vm_pu * Sparlectra.getNodeVn(nd))
      sv.attrs[:angle] = string(nd._va_deg)
    end
    t0 = DateTime(2026, 1, 1, 12)
    inputs = [CGMESSeriesInput(s1; label = "s1", scenario_time = t0, model_key = "m"), CGMESSeriesInput(s2; label = "s2", scenario_time = t0 + Minute(60), model_key = "m")]
    u = redirect_stdout(() -> buildCGMESSeriesUnion(inputs; config = cfg), devnull)
    # template plus delta starts every live bus at its start in the direct
    # import of s2, which is the SV; the base solution lies away from it
    r2 = redirect_stdout(() -> importCGMES(cfg; path = s2, name = "s2"), devnull)
    rk_bus = Dict(tn => r2.net.busDict[name] for (tn, name) in r2.topo.bus_name if haskey(r2.net.busDict, name))
    d2 = u.deltas[2]
    net = deepcopy(u.template)
    Sparlectra.apply!(net, vcat(d2.ops, d2.control_ops, d2.start_ops), u.index)
    base = deepcopy(u.template)
    @test last(redirect_stdout(() -> runpf!(base, 50, 1e-10, 0), devnull)) == 0
    dv = 0.0; da = 0.0; gap = 0.0; nb = 0
    for (b, nid) in u.node_of_bus
      rb = get(rk_bus, u.node_sig[nid][2], 0)
      (rb == 0 || Sparlectra.isIsolated(r2.net.nodeVec[rb])) && continue
      dv = max(dv, abs(net.nodeVec[b]._vm_pu - r2.net.nodeVec[rb]._vm_pu))
      da = max(da, abs(net.nodeVec[b]._va_deg - r2.net.nodeVec[rb]._va_deg))
      gap = max(gap, abs(base.nodeVec[b]._vm_pu - r2.net.nodeVec[rb]._vm_pu))
      nb += 1
    end
    @test nb > 0 && dv == 0.0 && da == 0.0 && gap > 1e-3
    # mode A: step 2 starts at its solution (one iteration; from the base
    # solution it took four)
    r = redirect_stdout(() -> runQsts(u; config = qcfg), devnull)
    @test all(s -> s.status === :solved, r.steps)
    @test r.steps[2].iterations <= 1
  end)() end

  # the union reads its snapshots one at a time, a path source in both
  # passes; from paths it must equal the union of the same snapshots given as
  # stores. Step 2 is the base case written as an xml folder with the loads
  # of its SSH scaled
  @testset "CGMES series union: path sources" begin (function ()
    C = Sparlectra.CGMESImporter
    dir = joinpath(@__DIR__, "fixtures", "cgmes_smallgrid")
    files = sort!(filter(f -> endswith(f, ".zip"), readdir(dir; join = true)))
    cfg = Sparlectra.SparlectraConfig(Dict())
    t0 = DateTime(2026, 1, 1, 12)
    mktempdir() do d
      zr = ZipArchives.ZipReader(read(only(filter(f -> occursin("BaseCase", f), files))))
      for e in ZipArchives.zip_names(zr)
        xml = String(ZipArchives.zip_readentry(zr, e))
        if occursin("_SSH_", e)
          xml = replace(xml, r"<cim:EnergyConsumer\.[pq]>[^<]*<" => function (tag)
            local head, val = split(chop(tag), '>')
            return string(head, ">", parse(Float64, val) * 1.1, "<")
          end)
        end
        write(joinpath(d, e), xml)
      end
      paths = [files, [d, only(filter(f -> occursin("Boundary", f), files))]]
      stores = [redirect_stdout(() -> C.loadCGMES(p), devnull) for p in paths]
      series(srcs) = [CGMESSeriesInput(srcs[k]; label = string("s", k), scenario_time = t0 + Minute(60 * (k - 1)), model_key = "m") for k in 1:2]
      up = redirect_stdout(() -> buildCGMESSeriesUnion(series(paths); config = cfg), devnull)
      us = redirect_stdout(() -> buildCGMESSeriesUnion(series(stores); config = cfg), devnull)
      @test up.index.id_by_mrid == us.index.id_by_mrid && up.index.kind_by_id == us.index.kind_by_id && up.index.internal_by_id == us.index.internal_by_id
      @test all(isequal(up.deltas[k].ops, us.deltas[k].ops) && isequal(up.deltas[k].residual, us.deltas[k].residual) for k in 1:2)
      # the written loads arrived: one op per nonzero load value
      @test count(op -> op.target === :load && op.field in (:p, :q), up.deltas[2].ops) == _union_scale_loads!(deepcopy(stores[1]), 1.1)
      for k in 1:2
        rk = redirect_stdout(() -> importCGMES(cfg; path = paths[k], name = string("s", k)), devnull)
        @test isempty(compareSnapshotReproduction(up, k, rk))
      end
    end
  end)() end

  # equipment parameters that an export changes between snapshots (measured:
  # lines with other r and x on the next day) become set ops of the delta,
  # not residual
  @testset "CGMES series union: changed branch parameters" begin (function ()
    C = Sparlectra.CGMESImporter
    dir = joinpath(@__DIR__, "fixtures", "cgmes_smallgrid")
    s1 = redirect_stdout(() -> C.loadCGMES(sort!(filter(f -> endswith(f, ".zip"), readdir(dir; join = true)))), devnull)
    s2 = deepcopy(s1)
    lines = sort!(copy(s2.byclass[:ACLineSegment]))[1:3]
    for m in lines, a in (:r, :x, :bch)
      v = C.num(s2.objects[m], a)
      v === nothing || (s2.objects[m].attrs[a] = string(v * 1.1))
    end
    cfg = Sparlectra.SparlectraConfig(Dict())
    t0 = DateTime(2026, 1, 1, 12)
    inputs = [CGMESSeriesInput(s1; label = "before", scenario_time = t0, model_key = "m"),
              CGMESSeriesInput(s2; label = "after", scenario_time = t0 + Minute(60), model_key = "m")]
    u = redirect_stdout(() -> buildCGMESSeriesUnion(inputs; config = cfg), devnull)
    d2 = u.deltas[2]
    ids = Set(u.index.id_by_mrid[m] for m in lines)
    @test isempty(d2.residual)
    @test Set(op.id for op in d2.ops) == ids
    @test all(op -> op.field in (:r_pu, :x_pu, :b_from_pu, :b_to_pu, :g_from_pu, :g_to_pu), d2.ops)
    @test all(any(op -> op.id == i && op.field === :x_pu, d2.ops) for i in ids)
    for (k, s) in enumerate((s1, s2))
      rk = redirect_stdout(() -> importCGMES(cfg; path = s, name = string("s", k)), devnull)
      @test isempty(compareSnapshotReproduction(u, k, rk))
    end
    # a phase-stamped EquivalentBranch (a :branch with ratio 1) whose r21, x21
    # change: the shift of the stamp follows as an op of its own
    e1 = redirect_stdout(() -> C.loadCGMES(_cgmes_synth_eb_dir()), devnull)
    e2 = deepcopy(e1)
    z21 = complex(1.5, 7.0) * cis(deg2rad(-16.0))
    e2.objects["_eb_asym"].attrs[:r21] = string(real(z21))
    e2.objects["_eb_asym"].attrs[:x21] = string(imag(z21))
    einputs = [CGMESSeriesInput(e1; label = "before", scenario_time = t0, model_key = "m"),
               CGMESSeriesInput(e2; label = "after", scenario_time = t0 + Minute(60), model_key = "m")]
    eu = redirect_stdout(() -> buildCGMESSeriesUnion(einputs; config = cfg), devnull)
    ed2 = eu.deltas[2]
    @test isempty(ed2.residual)
    # the snapshot run's element list names them equivalent branches; the
    # union builds every element whatever its SSH state, EB_OPEN included
    @test count(==(:equivalent_branch), Sparlectra._qsts_elements(eu.template, eu.index).branch_kind) == 4
    @test any(op -> op.id == eu.index.id_by_mrid["_eb_asym"] && op.field === :shift_deg && isapprox(op.value, -8.0; atol = 1e-9), ed2.ops)
    for (k, s) in enumerate((e1, e2))
      rk = redirect_stdout(() -> importCGMES(cfg; path = s, name = string("e", k)), devnull)
      @test isempty(compareSnapshotReproduction(eu, k, rk))
    end
    # the stamp comes and goes: EB_SAME gets an r21, x21 of its own (a line
    # in step 1, stamped in step 2), EB_ASYM loses it (stamped, then a line).
    # The delta had no op for the change, so step 2 kept the line without
    # shift and the old shift; each step must solve to its direct import
    e3 = deepcopy(e1)
    z21 = complex(1.21, 6.05) * cis(deg2rad(10.0))
    e3.objects["_eb_same"].attrs[:r21] = string(real(z21))
    e3.objects["_eb_same"].attrs[:x21] = string(imag(z21))
    e3.objects["_eb_asym"].attrs[:r21] = e3.objects["_eb_asym"].attrs[:r]
    e3.objects["_eb_asym"].attrs[:x21] = e3.objects["_eb_asym"].attrs[:x]
    sinputs = [CGMESSeriesInput(e1; label = "before", scenario_time = t0, model_key = "m"),
               CGMESSeriesInput(e3; label = "after", scenario_time = t0 + Minute(60), model_key = "m")]
    su = redirect_stdout(() -> buildCGMESSeriesUnion(sinputs; config = cfg), devnull)
    @test all(d -> isempty(d.residual), su.deltas)
    for (k, s) in enumerate((e1, e3))
      rk = redirect_stdout(() -> importCGMES(cfg; path = s, name = string("s", k)), devnull)
      @test isempty(compareSnapshotReproduction(su, k, rk))
      net = deepcopy(su.template)
      Sparlectra.apply!(net, vcat(su.deltas[k].ops, su.deltas[k].control_ops), su.index)
      @test last(redirect_stdout(() -> runpf!(net, 50, 1e-10, 0), devnull)) == 0
      @test last(redirect_stdout(() -> runpf!(rk.net, 50, 1e-10, 0), devnull)) == 0
      rk_bus = Dict(tn => rk.net.busDict[name] for (tn, name) in rk.topo.bus_name if haskey(rk.net.busDict, name))
      dv = 0.0; da = 0.0; nb = 0
      for (b, nid) in su.node_of_bus
        rb = get(rk_bus, su.node_sig[nid][k], 0)
        (rb == 0 || Sparlectra.isIsolated(rk.net.nodeVec[rb])) && continue
        dv = max(dv, abs(net.nodeVec[b]._vm_pu - rk.net.nodeVec[rb]._vm_pu))
        da = max(da, abs(net.nodeVec[b]._va_deg - rk.net.nodeVec[rb]._va_deg))
        nb += 1
      end
      @test nb == 4 && dv < 1e-9 && da < 1e-7
    end
  end)() end

  # a tap step that changes between snapshots: the SSH step is kept as
  # reference, and the ratio goes to control_ops when the regulator of that
  # transformer is active (a sequential run leaves it to the controller),
  # else to ops
  @testset "CGMES series union: controlled ratio, SSH steps and tap data" begin (function ()
    C = Sparlectra.CGMESImporter
    s1 = redirect_stdout(() -> C.loadCGMES(_cgmes_synth3w_dir()), devnull)
    s2 = deepcopy(s1)
    s2.objects["_rtc2"].attrs[:step] = "13"
    t0 = DateTime(2026, 1, 1, 12)
    inputs = [CGMESSeriesInput(s1; label = "neutral", scenario_time = t0, model_key = "m"),
              CGMESSeriesInput(s2; label = "raised", scenario_time = t0 + Minute(60), model_key = "m")]
    for tap in (false, true)
      cfg = Sparlectra.SparlectraConfig(Dict{String,Any}("cgmes_import" => Dict{String,Any}("tap_control" => tap)))
      u = redirect_stdout(() -> buildCGMESSeriesUnion(inputs; config = cfg), devnull)
      d2 = u.deltas[2]
      leg = u.index.id_by_mrid["_t3w_e2"]
      @test d2.ssh_steps["_rtc2"] == 13.0 && u.deltas[1].ssh_steps["_rtc2"] == 11.0
      @test any(op -> op.id == leg && op.field === :ratio, tap ? d2.control_ops : d2.ops)
      @test !any(op -> op.id == leg && op.field === :ratio, tap ? d2.ops : d2.control_ops)
      rk = redirect_stdout(() -> importCGMES(cfg; path = s2, name = "raised"), devnull)
      @test isempty(compareSnapshotReproduction(u, 2, rk))
    end
    # tap changer data that differs from the template: the importer sets
    # range, step, table ratios and impedance per step for every changer,
    # and a step whose EQ changes the range or the table has other data.
    # Each step must carry the tap fields of its direct import exactly
    # (compared field by field, a record's table by its fields)
    tap_fields = (:has_ratio_tap, :tap_min, :tap_max, :tap_step, :tap_ratio_points, :has_phase_tap, :phase_min_deg, :phase_max_deg, :phase_step_deg, :tap_step_impedance)
    tap_value(br, f) = (v = getfield(br, f); v isa Sparlectra.TapStepImpedance ? [getfield(v, g) for g in fieldnames(typeof(v))] : v)
    function tap_steps(ins, mrid_of)
      tcfg = Sparlectra.SparlectraConfig(Dict{String,Any}("cgmes_import" => Dict{String,Any}("tap_control" => true)))
      u = redirect_stdout(() -> buildCGMESSeriesUnion(ins; config = tcfg), devnull)
      for k in eachindex(ins)
        rk = redirect_stdout(() -> importCGMES(tcfg; path = ins[k].source, name = ins[k].label), devnull)
        m = mrid_of(rk)
        net = deepcopy(u.template)
        Sparlectra.apply!(net, vcat(u.deltas[k].ops, u.deltas[k].control_ops), u.index)
        tb = net.branchVec[u.index.internal_by_id[u.index.id_by_mrid[m]]]
        rb = rk.net.branchVec[rk.element_by_mrid[m][2]]
        @test (ins[k].label, [f for f in tap_fields if !isequal(tap_value(tb, f), tap_value(rb, f))]) == (ins[k].label, Symbol[])
        @test (ins[k].label, u.deltas[k].residual, compareSnapshotReproduction(u, k, rk)) == (ins[k].label, String[], String[])
      end
      return u
    end
    # the regulated leg of the three-winding transformer: a shorter range,
    # then its regulator off with another range (an unregulated changer
    # carries its range as well)
    s3 = deepcopy(s1)
    s3.objects["_rtc2"].attrs[:highStep] = "17"
    s4 = deepcopy(s1)
    s4.objects["_rtc2"].attrs[:controlEnabled] = "false"
    s4.objects["_rtc2"].attrs[:highStep] = "19"
    u = tap_steps([CGMESSeriesInput(s; label = l, scenario_time = t0 + Minute(60 * (k - 1)), model_key = "m") for (k, (s, l)) in enumerate(((s1, "neutral"), (s3, "short range"), (s4, "unregulated")))],
                  rk -> "_t3w_e2")
    leg = u.index.id_by_mrid["_t3w_e2"]
    @test [count(op -> op.id == leg && op.field === :tap_data, d.ops) for d in u.deltas] == [0, 1, 1] && length(u.index.tap_data[leg]) == 2
    # a two-winding transformer whose regulated changer has a table with r,
    # x, g, b deviations: the table of the snapshot has other ratios and
    # twice the deviations, so the table ratios and the impedance per step
    # differ from the template
    w1 = redirect_stdout(() -> C.loadCGMES(_cgmes_synth_tapside_dir(:two_winding; table = (r = 2.0, x = 4.0, g = 1.0, b = 1.5), control_kv = 20.6)), devnull)
    w2 = deepcopy(w1)
    for m in w2.byclass[:RatioTapChangerTablePoint]
      o = w2.objects[m]
      o.attrs[:ratio] = string(1.0 + (C.num(o, :ratio) - 1.0) * 1.2)
      for a in (:r, :x, :g, :b)
        o.attrs[a] = string(2.0 * C.num(o, a))
      end
    end
    wu = tap_steps([CGMESSeriesInput(w1; label = "table", scenario_time = t0, model_key = "m"), CGMESSeriesInput(w2; label = "new table", scenario_time = t0 + Minute(60), model_key = "m")],
                   rk -> only(m for (m, (kind, _)) in rk.element_by_mrid if kind === :transformer))
    @test only(only(values(wu.index.tap_data))).tap_step_impedance !== nothing
  end)() end

  @testset "CGMES series union: a bus behind a closed switch in mode A" begin (function ()
    # C hangs on B by a closed switch alone; the import types it Isolated
    # (branches only) and the delta carries that, but the power flow solves
    # C with B: mode A gives C the values of B, and none once the switch
    # opens in the second snapshot
    C = Sparlectra.CGMESImporter
    net = link_group_net(; ring = true)
    redirect_stdout(() -> runpf!(net, 30, 1e-8, 0), devnull)
    dir = mktempdir()
    redirect_stdout(() -> writeCGMESFiles(net; path = dir), devnull)
    s1 = redirect_stdout(() -> C.loadCGMES(sort!(readdir(dir; join = true))), devnull)
    s2 = deepcopy(s1)
    sw = only(m for cls in C._SWITCH_CLASSES for m in get(s2.byclass, cls, String[]))
    s2.objects[sw].attrs[:open] = "true"
    t0 = DateTime(2026, 1, 1, 12)
    inputs = [CGMESSeriesInput(s1; label = "closed", scenario_time = t0, model_key = "m"),
              CGMESSeriesInput(s2; label = "open", scenario_time = t0 + Minute(60), model_key = "m")]
    cfg = Sparlectra.SparlectraConfig(Dict{String,Any}())
    u = redirect_stdout(() -> buildCGMESSeriesUnion(inputs; config = cfg), devnull)
    bB, bC = u.template.busDict["B"], u.template.busDict["C"]
    @test Sparlectra.isIsolated(u.template.nodeVec[bC])
    r = redirect_stdout(() -> runQsts(u; config = Sparlectra.SparlectraConfig(Dict{String,Any}("qsts" => Dict{String,Any}("mode" => "snapshots")))), devnull)
    @test all(s -> s.status === :solved, r.steps)
    @test isfinite(r.bus_values[bC, 1, 1]) && r.bus_values[bC, 1:2, 1] == r.bus_values[bB, 1:2, 1]
    @test all(isnan, r.bus_values[bC, :, 2]) && isfinite(r.bus_values[bB, 1, 2])
  end)() end

  @testset "CGMES series union: an isolated bus is the smallest of a live switch group" begin (function ()
    # the line C - D opens in the second snapshot, C keeps only its switch
    # to B: the delta isolates C (0 V) and the solver starts the merged
    # node {C, B} from C, its smallest bus. The step starts C from B, so it
    # converges to the high-voltage solution, not to the low-voltage branch
    # a 0 V start ends on.
    net, s1, s2 = isolated_first_cgmes_snapshots()
    t0 = DateTime(2026, 1, 1, 12)
    inputs = [CGMESSeriesInput(s1; label = "closed", scenario_time = t0, model_key = "m"),
              CGMESSeriesInput(s2; label = "open", scenario_time = t0 + Minute(60), model_key = "m")]
    cfg = Sparlectra.SparlectraConfig(Dict{String,Any}())
    u = redirect_stdout(() -> buildCGMESSeriesUnion(inputs; config = cfg), devnull)
    bB, bC = u.template.busDict["B"], u.template.busDict["C"]
    @test bC == 1 && Sparlectra._active_link_representative_map(u.template)[bB] == bC
    @test any(op -> u.index.internal_by_id[op.id] == bC && op.field === :isolated && op.value == 1.0, u.deltas[2].ops)
    r = redirect_stdout(() -> runQsts(u; config = Sparlectra.SparlectraConfig(Dict{String,Any}("qsts" => Dict{String,Any}("mode" => "snapshots")))), devnull)
    @test all(s -> s.status === :solved, r.steps)
    ref = isolated_first_reference(net)
    @test r.bus_values[bB, 1, 2] ≈ ref._vm_pu atol = 1e-6
    @test r.bus_values[bB, 2, 2] ≈ ref._va_deg atol = 1e-4
    @test r.bus_values[bC, 1:2, 2] == r.bus_values[bB, 1:2, 2]
  end)() end

  @testset "CGMES series union: an island split only by a switch-only bus" begin (function ()
    # the re-import isolates JX (no branch of its own), the solver joins
    # both sides at JX into one island with the reference at J1. QSTS and
    # the SV comparison judge the islands of the power flow: the J3 - J4
    # side is not frozen, gets no start values of its own (only the slack
    # bus does) and is aligned with the offset of the whole island
    C = Sparlectra.CGMESImporter
    net = _union_switch_only_island_net()
    redirect_stdout(() -> runpf!(net, 30, 1e-8, 0), devnull)
    dir = mktempdir()
    redirect_stdout(() -> writeCGMESFiles(net; path = dir), devnull)
    s1 = redirect_stdout(() -> C.loadCGMES(sort!(readdir(dir; join = true))), devnull)
    s2 = deepcopy(s1)
    _union_scale_loads!(s2, 1.2)
    # a different SV state in step 2: start values become ops
    for o in C.objectsOf(s2, :SvVoltage)
      o.attrs[:angle] = string(C.num(o, :angle, 0.0) + 0.5)
    end
    t0 = DateTime(2026, 1, 1, 12)
    inputs = [CGMESSeriesInput(s1; label = "base", scenario_time = t0, model_key = "m"),
              CGMESSeriesInput(s2; label = "loads", scenario_time = t0 + Minute(60), model_key = "m")]
    cfg = Sparlectra.SparlectraConfig(Dict{String,Any}())
    rk = redirect_stdout(() -> importCGMES(cfg; path = s1, name = "base"), devnull)
    jx = only(rk.net.isoNodes)
    @test haskey(Sparlectra._linkFedIsolatedBuses(rk.net), jx)
    # the two sides as the branch-only electrical view splits them
    sides = Sparlectra.electricalIslandComponents(rk.net)
    far = only(c for c in sides if !any(b -> Sparlectra.getNodeType(rk.net.nodeVec[b]) == Sparlectra.Slack, c))
    near = only(c for c in sides if c !== far)
    name_of = Dict(i => n for (n, i) in rk.net.busDict)
    u = redirect_stdout(() -> buildCGMESSeriesUnion(inputs; config = cfg), devnull)
    far_u = [u.template.busDict[name_of[b]] for b in far]
    started = Set(u.index.internal_by_id[op.id] for op in u.deltas[2].ops if op.target === :bus && op.field in (:vm_start, :va_start))
    @test !isempty(started) && isempty(intersect(started, far_u))
    r = redirect_stdout(() -> runQsts(u; config = Sparlectra.SparlectraConfig(Dict{String,Any}("qsts" => Dict{String,Any}("mode" => "snapshots")))), devnull)
    qs(q) = findfirst(==(q), Sparlectra.QSTS_STEP_QUANTITIES)
    @test all(s -> s.status === :solved, r.steps)
    @test r.step_values[qs("n_frozen"), :] == [0.0, 0.0]
    @test all(isfinite, r.bus_values[far_u, 1, :])
    # direct import: J3 is one node with J2 and JX and gets their aligned
    # angle difference
    redirect_stdout(() -> runpf!(rk.net, 30, 1e-8, 0), devnull)
    c = compareWithSV(rk)
    row(b) = only(x for x in c.rows if x.bus == name_of[b])
    linked = only(b for b in far if any(l -> l.status == 1 && (Int(l.fromBus) == b || Int(l.toBus) == b), rk.net.linkVec))
    @test row(linked).dva_aligned ≈ row(jx).dva_aligned atol = 1e-9
    # the SV metrics of the base step equal compareWithSV of the direct
    # import, JX included at the voltage of its switch group
    @test r.step_values[qs("sv_n_nodes"), 1] == c.n
    @test r.step_values[qs("sv_max_dvm_pu"), 1] ≈ c.max_dvm atol = 1e-6
    @test r.step_values[qs("sv_max_dva_deg"), 1] ≈ c.max_dva atol = 1e-6
  end)() end

  # a switch inside one topological node is no link of the union, so an event
  # could not open it; extra_link_switches splits the node at the switch. On
  # the SmallGrid the two line-end switches of _union_open_line! become links,
  # both steps are reproduced, and opening one of them on the template gives
  # the solved direct import of the snapshot with that switch open
  @testset "CGMES series union: extra link switches" begin (function ()
    C = Sparlectra.CGMESImporter
    dir = joinpath(@__DIR__, "fixtures", "cgmes_smallgrid")
    s1 = redirect_stdout(() -> C.loadCGMES(sort!(filter(f -> endswith(f, ".zip"), readdir(dir; join = true)))), devnull)
    s2 = deepcopy(s1)
    _union_scale_loads!(s2, 1.1)
    _, both = _union_open_line!(deepcopy(s1))
    so = deepcopy(s1)
    _, one = _union_open_line!(so; ends = 1:1)
    @test length(both) == 2 && one == both[1:1]
    cfg = Sparlectra.SparlectraConfig(Dict{String,Any}())
    t0 = DateTime(2026, 1, 1, 12)
    inputs = [CGMESSeriesInput(s1; label = "s1", scenario_time = t0, model_key = "m"),
              CGMESSeriesInput(s2; label = "s2", scenario_time = t0 + Minute(60), model_key = "m")]
    nnodes = length(C._tpTNs(s1)) + count(m -> m in s1.boundary, s1.byclass[:TopologicalNode])
    u0 = redirect_stdout(() -> buildCGMESSeriesUnion(inputs; config = cfg), devnull)
    @test length(u0.node_mrid) == nnodes && isempty(u0.link_switches)
    u = redirect_stdout(() -> buildCGMESSeriesUnion(inputs; config = cfg, extra_link_switches = both), devnull)
    @test length(u.node_mrid) == nnodes + 2
    @test sort(u.link_switches) == sort(both) && u.n_tp_links == 0
    @test all(u.link_status[k][m] == 1 for m in both, k in 1:2)
    @test isempty(u.messages) && isempty(u.deltas[1].ops)
    for (k, s) in enumerate((s1, s2))
      rk = redirect_stdout(() -> importCGMES(cfg; path = s, name = string("s", k)), devnull)
      @test isempty(compareSnapshotReproduction(u, k, rk))
    end
    # the event: one line end opened by a link op on the template
    closed = deepcopy(u.template)
    net = deepcopy(u.template)
    Sparlectra.apply!(net, [Sparlectra.PatchOp(op = :status, target = :link, id = u.index.id_by_mrid[only(one)], value = 0.0)], u.index)
    @test last(redirect_stdout(() -> runpf!(closed, 30, 1e-10, 0), devnull)) == 0
    @test last(redirect_stdout(() -> runpf!(net, 30, 1e-10, 0), devnull)) == 0
    ro = redirect_stdout(() -> importCGMES(cfg; path = so, name = "open"), devnull)
    @test last(redirect_stdout(() -> runpf!(ro.net, 30, 1e-10, 0), devnull)) == 0
    # union bus -> bus of the open snapshot over a terminal of its node
    term_tn = C._seriesTopo(so, "open").term_tn
    ro_bus = Dict(tn => ro.net.busDict[name] for (tn, name) in ro.topo.bus_name if haskey(ro.net.busDict, name))
    node_terms = Dict{Int,String}()
    for (t, nid) in u.node_of_terminal
      node_terms[nid] = min(get(node_terms, nid, t), t)
    end
    dv = 0.0; da = 0.0; moved = 0.0; nb = 0
    for (b, nid) in u.node_of_bus
      rb = get(ro_bus, get(term_tn, node_terms[nid], ""), 0)
      rb == 0 && continue
      dv = max(dv, abs(net.nodeVec[b]._vm_pu - ro.net.nodeVec[rb]._vm_pu))
      da = max(da, abs(net.nodeVec[b]._va_deg - ro.net.nodeVec[rb]._va_deg))
      moved = max(moved, abs(net.nodeVec[b]._vm_pu - closed.nodeVec[b]._vm_pu))
      nb += 1
    end
    @test nb == length(u.node_of_bus) && dv <= 1e-8 && da <= 1e-6
    @test moved > 1e-3
    # a name that is no switch of the series is refused with the name
    line = first(sort!(copy(s1.byclass[:ACLineSegment])))
    @test_throws "extra_link_switches: $(line) is not a switch of the series" buildCGMESSeriesUnion(inputs; config = cfg, extra_link_switches = [line])
  end)() end
end
