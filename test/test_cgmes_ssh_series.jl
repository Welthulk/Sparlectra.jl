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

# file: test/test_cgmes_ssh_series.jl
# purpose: synthetic SSH series of tools/cgmes_ssh_series.jl on the SmallGrid
#          fixture: 96 bundles from one snapshot with every other file byte
#          for byte, scenario times in 15 minute steps, equal seeds give equal
#          SSH bytes, the more specific role key wins, the thermal balance
#          holds, profile ranges, the range of one role key and thermal
#          bounds from the [range] table are kept and a binding bound is
#          reported per step, PowerElectronicsConnection objects take the
#          role of their unit class, the series check finds one series without
#          churn, and the union of the series reproduces every step with p/q
#          ops plus the link ops of a switch event. The SmallGrid has no open
#          switch and no
#          switch between two topological nodes (measured: all 1107 switches
#          closed, each inside one TN), so the base is derived first: one
#          line-end disconnector opened in the SSH and its line-side node
#          split off in the TP. Events close it at step 40 and open it again
#          at step 60. The derived base also writes the SSH model id as
#          #_<id> and lets the SV depend on it, so the manifest ids have to
#          follow the header reader. A short run from the two fixture zips
#          covers zip input, the zip!entry sources and an identical file
#          given twice.

using Dates
using TOML
using ZipArchives
isdefined(@__MODULE__, :generateCGMESSSHSeries) || include(joinpath(@__DIR__, "..", "tools", "cgmes_ssh_series.jl"))

# replace `from` by `to` inside the block of object `mrid` (rdf:about);
# false when the file has no such block or the block no such text
function _sshs_test_block_replace(text::String, mrid::String, from::String, to::String)
  i = findfirst("rdf:about=\"#$(mrid)\"", text)
  i === nothing && return text, false
  stop = findnext("</cim:", text, last(i))
  j = findnext(from, text, last(i))
  (j === nothing || first(j) > first(stop)) && return text, false
  return string(text[1:(first(j) - 1)], to, text[(last(j) + 1):end]), true
end

# The SmallGrid with one open switch between two topological nodes: at the
# first line end (lines and switches by mRID) whose node splits when one of
# its switches opens, that switch is opened in the SSH, and the connectivity
# nodes that still reach the line end over closed switches move to a new TN
# in the TP. The SSH model id is rewritten from urn:uuid:<id> to #_<id> and
# the SV declares a DependentOn on it. Written as plain XML files into `dir`.
# Returns the switch, the old and the new TN.
function _sshs_test_split_base(src::Vector{String}, dir::String)
  C = Sparlectra.CGMESImporter
  s = redirect_stdout(() -> C.loadCGMES(src), devnull)
  by_eq = Dict{String,Vector{C.CIMObject}}()
  by_cn = Dict{String,Vector{C.CIMObject}}()
  for t in C.objectsOf(s, :Terminal)
    push!(get!(Vector{C.CIMObject}, by_eq, t.refs[:ConductingEquipment]), t)
    haskey(t.refs, :ConnectivityNode) && push!(get!(Vector{C.CIMObject}, by_cn, t.refs[:ConnectivityNode]), t)
  end
  cn_tn = Dict(cn.mrid => cn.refs[:TopologicalNode] for cn in C.objectsOf(s, :ConnectivityNode) if haskey(cn.refs, :TopologicalNode))
  isswitch(m) = s.objects[m].class in C._SWITCH_CLASSES
  found = nothing
  for line in sort!(copy(s.byclass[:ACLineSegment])), t in sort(by_eq[line]; by = x -> C.num(x, :sequenceNumber, 1.0))
    cn = t.refs[:ConnectivityNode]
    old = cn_tn[cn]
    for sw in sort!([x.refs[:ConductingEquipment] for x in by_cn[cn] if isswitch(x.refs[:ConductingEquipment])])
      other = only(x.refs[:ConnectivityNode] for x in by_eq[sw] if x.refs[:ConnectivityNode] != cn)
      reach = Set([cn])
      stack = [cn]
      while !isempty(stack)
        c = pop!(stack)
        for x in by_cn[c], y in by_eq[x.refs[:ConductingEquipment]]
          e = x.refs[:ConductingEquipment]
          c2 = y.refs[:ConnectivityNode]
          (e != sw && isswitch(e) && !(c2 in reach) && get(cn_tn, c2, "") == old) || continue
          push!(reach, c2)
          push!(stack, c2)
        end
      end
      other in reach && continue
      found = (sw, old, reach)
      break
    end
    # a break in a loop over two iterators leaves both
    found === nothing || break
  end
  found === nothing && error("SmallGrid: no line-end switch splits its topological node")
  sw, old, reach = found
  new = "_sshs_test_split_TN"
  moved = sort!(collect(union(reach, Set(x.mrid for c in reach for x in by_cn[c]))))
  mkpath(dir)
  src_files = C.collectCGMESFiles(src)
  heads = [C.readCGMESHeader(f) for f in src_files]
  ssh_id = only(h.model_id for h in heads if :SSH in h.profiles)
  for (f, h) in zip(src_files, heads)
    text = f.content
    profiles = h.profiles
    if :TP in profiles
      n = 0
      for m in moved
        text, ok = _sshs_test_block_replace(text, m, "rdf:resource=\"#$(old)\"", "rdf:resource=\"#$(new)\"")
        n += ok
      end
      n >= length(reach) || error("SmallGrid TP: only $(n) references moved")
      i = findfirst("<cim:TopologicalNode rdf:ID=\"$(old)\">", text)
      j = findnext("</cim:TopologicalNode>", text, last(i))
      block = replace(text[first(i):last(j)], "rdf:ID=\"$(old)\"" => "rdf:ID=\"$(new)\"", r"(<cim:IdentifiedObject.name>[^<]*)" => s"\1 split")
      k = findlast("</rdf:RDF>", text)
      text = string(text[1:(first(k) - 1)], "  ", block, "\n", text[first(k):end])
    elseif :SSH in profiles
      text, ok = _sshs_test_block_replace(text, sw, "<cim:Switch.open>false</cim:Switch.open>", "<cim:Switch.open>true</cim:Switch.open>")
      ok || error("SmallGrid SSH: no Switch.open for $(sw)")
      occursin("rdf:about=\"urn:uuid:$(ssh_id)\"", text) || error("SmallGrid SSH: model id $(ssh_id) not written as urn:uuid")
      text = replace(text, "rdf:about=\"urn:uuid:$(ssh_id)\"" => "rdf:about=\"#_$(ssh_id)\""; count = 1)
    elseif :SV in profiles
      text = replace(text, "</md:FullModel>" => "  <md:Model.DependentOn rdf:resource=\"#_$(ssh_id)\" />\n  </md:FullModel>"; count = 1)
    end
    write(joinpath(dir, basename(f.name)), text)
  end
  return sw, old, new
end

# A copy of `base` with three PowerElectronicsConnection objects (the
# SmallGrid has none): one with a PowerElectronicsWindUnit, one with a
# PhotoVoltaicUnit, one with a BatteryUnit and a PhotoVoltaicUnit. EQ gets
# the connections and their units (each unit points at its connection), SSH
# their p and q in load sign. Written as plain XML files into `dir`.
function _sshs_test_pec_base(base::String, dir::String)
  C = Sparlectra.CGMESImporter
  pecs = [("_sshs_pec_wind", ["PowerElectronicsWindUnit"], "-20", "-2"), ("_sshs_pec_pv", ["PhotoVoltaicUnit"], "-10", "0"), ("_sshs_pec_mix", ["BatteryUnit", "PhotoVoltaicUnit"], "-5", "-1")]
  eq = join(("""
      <cim:PowerElectronicsConnection rdf:ID="$(m)">
        <cim:IdentifiedObject.name>$(m)</cim:IdentifiedObject.name>
      </cim:PowerElectronicsConnection>
    """ * join(("""
      <cim:$(u) rdf:ID="$(m)_u$(i)">
        <cim:PowerElectronicsUnit.PowerElectronicsConnection rdf:resource="#$(m)" />
      </cim:$(u)>
    """ for (i, u) in enumerate(us))) for (m, us, _, _) in pecs))
  ssh = join(("""
      <cim:PowerElectronicsConnection rdf:about="#$(m)">
        <cim:PowerElectronicsConnection.p>$(p)</cim:PowerElectronicsConnection.p>
        <cim:PowerElectronicsConnection.q>$(q)</cim:PowerElectronicsConnection.q>
      </cim:PowerElectronicsConnection>
    """ for (m, _, p, q) in pecs))
  mkpath(dir)
  for f in C.collectCGMESFiles([base])
    profiles = C.readCGMESHeader(f).profiles
    add = (:EQ in profiles && !(:EQ_BD in profiles)) ? eq : (:SSH in profiles ? ssh : "")
    k = findlast("</rdf:RDF>", f.content)
    write(joinpath(dir, basename(f.name)), isempty(add) ? f.content : string(f.content[1:(first(k) - 1)], add, f.content[first(k):end]))
  end
  return [m for (m, _, _, _) in pecs]
end

# the value of one metric row of a QSTS result, NaN without the row
function _qsts_metric_value(r, id::String, metric::String)::Float64
  m = qstsMetrics(r)
  i = findfirst(j -> m.element_id[j] == id && m.metric[j] == metric, eachindex(m.metric))
  return i === nothing ? NaN : m.value[i]
end

# The SmallGrid split into a base in `tmp` and the 96-step SSH series of
# these tests (role file, switch events), base and steps loaded; `src` are
# the fixture zips.
function _sshs_test_series(tmp::String)
  C = Sparlectra.CGMESImporter
  src = sort!(filter(f -> endswith(f, ".zip"), readdir(joinpath(@__DIR__, "fixtures", "cgmes_smallgrid"); join = true)))
  base = joinpath(tmp, "base")
  sw, old_tn, new_tn = _sshs_test_split_base(src, base)
  # the SmallGrid has no wind unit: the plain machine key says wind, the
  # unit-class key says thermal and wins for every machine; the equivalent
  # injections carry the wind ramp
  roles = joinpath(tmp, "roles.toml")
  write(roles, """
    EnergyConsumer = "load"
    EquivalentInjection = "wind"
    SynchronousMachine = "wind"
    "SynchronousMachine:ThermalGeneratingUnit" = "thermal"
    """)
  events = joinpath(tmp, "events.csv")
  write(events, "step,element_id,action,value\n40,$(sw),close,\n60,$(sw),open,\n")
  r = generateCGMESSSHSeries(base, joinpath(tmp, "series"); roles = roles, events = events)
  base_store = redirect_stdout(() -> C.loadCGMES(base), devnull)
  stores = [redirect_stdout(() -> C.loadCGMES(b), devnull) for b in r.bundles]
  return (; src, base, sw, old_tn, new_tn, roles, events, r, base_store, stores)
end

# One tap controller for mode B on the union `u`: transformer
# B_2WT_220_73_39 (branch 180 of the union template) holds its 132 kV side
# (Olive2) at 1.0 pu within 0.005 pu.
function _sshs_olive_controllers(u)
  bus = only(k for (k, v) in u.template.busDict if v == u.template.branchVec[180].toBus)
  return Dict{String,Any}("olive" => Dict{String,Any}("type" => "power_transformer", "trafo" => "180", "mode" => "voltage", "target_bus" => bus,
                                                      "target_vm_pu" => 1.0, "deadband_vm_pu" => 0.005))
end

# mode B starts each step from the step before (sequential_start = previous),
# the own SV start of mode A does not apply
_sshs_sequential_config(ctl) = Sparlectra.SparlectraConfig(Dict{String,Any}("qsts" => Dict{String,Any}("mode" => "sequential", "sequential_start" => "previous"), "control" => Dict{String,Any}("controllers" => ctl)))

# two unions on worker threads against the serial mode B run `r` of `u`
function _sshs_two_unions_check(u, seq, r)
  two = redirect_stdout(() -> runQsts([u, u]; config = seq, parallel_enabled = true, parallel_max_tasks = 2, parallel_min_work_items = 2), devnull)
  @test all(isequal(t.controller_values, r.controller_values) && isequal(t.bus_values, r.bus_values) && t.switching == r.switching for t in two)
  parallel_comparison_report("QSTS mode B on two CGMES unions")
  return two
end

function run_cgmes_ssh_series_tests()
  C = Sparlectra.CGMESImporter
  mktempdir() do tmp
    (; src, base, sw, old_tn, new_tn, roles, events, r, base_store, stores) = _sshs_test_series(tmp)
    ssh_name = only(f.name for f in base_store.files if :SSH in f.profiles)

    @testset "CGMES SSH series: bundles, bytes and manifest" begin (function ()
      @test length(r.bundles) == 96 && all(isfile, r.bundles) && issorted(r.bundles)
      # every entry but the SSH byte for byte as in the base, the SSH changed
      local base_bytes = Dict(f => read(joinpath(base, f)) for f in readdir(base))
      local same = true
      local ssh_changed = true
      for b in r.bundles
        zr = ZipReader(read(b))
        same &= sort(zip_names(zr)) == sort!(collect(keys(base_bytes)))
        for n in zip_names(zr)
          n == ssh_name ? (ssh_changed &= zip_readentry(zr, n) != base_bytes[n]) : (same &= zip_readentry(zr, n) == base_bytes[n])
        end
      end
      @test same && ssh_changed
      # scenario times in 15 minute steps from the base time, a new model id per step
      local head = only(f for f in base_store.files if :SSH in f.profiles)
      local ssh_heads = [only(f for f in st.files if :SSH in f.profiles) for st in stores]
      @test [DateTime(h.scenario_time) for h in ssh_heads] == [DateTime(head.scenario_time) + Minute(15 * (k - 1)) for k in 1:96]
      @test r.scenario_times == [h.scenario_time for h in ssh_heads]
      @test length(unique(h.model_id for h in ssh_heads)) == 96 && !any(h -> h.model_id == head.model_id, ssh_heads)
      local man = TOML.parsefile(r.manifest)
      @test [joinpath(dirname(r.manifest), b["target"]) for b in man["bundle"]] == r.bundles
      @test [b["scenario_time"] for b in man["bundle"]] == r.scenario_times
      @test all(length(b["files"]) == length(base_bytes) for b in man["bundle"])
      # ids as the header reader gives them: the SSH written as #_<id> keeps
      # its underscore in the manifest, and the SV's dependency on the base
      # SSH, which no bundle carries, is listed and noted
      @test startswith(head.model_id, "_") && all(h -> startswith(h.model_id, "_"), ssh_heads)
      @test [only(f["uuid"] for f in b["files"] if f["kind"] == "SSH") for b in man["bundle"]] == [h.model_id for h in ssh_heads]
      @test all(b["missing_dependencies"] == [head.model_id] && any(n -> occursin("depends on the base SSH $(head.model_id)", n), b["notes"]) for b in man["bundle"])
      local gen = man["generator"]
      @test gen["seed"] == 1 && gen["steps"] == 96 && gen["step_minutes"] == 15 && gen["roles"]["SynchronousMachine:ThermalGeneratingUnit"] == "thermal"
      # equal seeds give equal SSH bytes; another seed gives other values
      local again = generateCGMESSSHSeries(base, joinpath(tmp, "again"); roles = roles, events = events)
      local seed2 = generateCGMESSSHSeries(base, joinpath(tmp, "seed2"); roles = roles, events = events, seed = 2)
      ssh_of(p) = zip_readentry(ZipReader(read(p)), ssh_name, String)
      body(s) = s[last(findfirst("</md:FullModel>", s)):end]
      # the base description (a licence notice in the SmallGrid) stays, the
      # generator text follows it
      function descr_of(s)
        local hit = match(r"<md:Model\.description>(.*?)</md:Model\.description>"s, s)
        return hit === nothing ? nothing : hit.captures[1]
      end
      local base_descr = descr_of(String(copy(base_bytes[ssh_name])))
      function descr_ok(k)
        local dk = descr_of(ssh_of(r.bundles[k]))
        return dk !== nothing && occursin("Synthetic SSH step $(k) of 96", dk) && (base_descr === nothing || startswith(dk, base_descr))
      end
      @test all(descr_ok, 1:96)
      @test all(ssh_of(a) == ssh_of(b) for (a, b) in zip(r.bundles, again.bundles))
      @test all(body(ssh_of(a)) != body(ssh_of(b)) for (a, b) in zip(r.bundles, seed2.bundles))
    end)() end

    @testset "CGMES SSH series: roles, profiles and balance" begin (function ()
      local machines = collect(C.objectsOf(base_store, :SynchronousMachine))
      @test all(r.role[o.mrid] == "thermal" for o in machines)
      # the plain machine key matches nothing next to the unit-class key, and the notes say so
      @test any(n -> occursin("role key SynchronousMachine matches no object", n), r.notes)
      @test sort!([m for (m, v) in r.role if v == "load"]) == sort(base_store.byclass[:EnergyConsumer])
      @test sort!([m for (m, v) in r.role if v == "wind"]) == sort(base_store.byclass[:EquivalentInjection])
      # the machine with referencePriority 1 keeps its values in every step
      local refm = only(o.mrid for o in machines if C.num(o, :referencePriority, 0.0) == 1.0)
      pq(st, m) = (C.num(st.objects[m], :p), C.num(st.objects[m], :q))
      @test r.reference == [refm] && all(pq(st, refm) == pq(base_store, refm) for st in stores)
      # balance: the sum of p over every load, injection and machine stays
      # at the base; six decimals per value keep it within 1e-4 MW
      psum(st) = sum(C.num(st.objects[m], :p, 0.0) for m in keys(r.role))
      local sums = [psum(st) for st in stores]
      @test maximum(abs.(sums .- psum(base_store))) <= 1e-4
      @test maximum(abs.(sums .- r.p_sum)) <= 1e-9
      # the profiles move: the load sum spans night and evening peak, the wind
      # ramps up and down, the thermal machines follow one factor
      csum(st, cls) = sum(C.num(o, :p, 0.0) for o in C.objectsOf(st, cls))
      local lr = [csum(st, :EnergyConsumer) / csum(base_store, :EnergyConsumer) for st in stores]
      @test minimum(lr) < 0.7 && maximum(lr) > 1.1
      local wr = [csum(st, :EquivalentInjection) / csum(base_store, :EquivalentInjection) for st in stores]
      @test wr[1] < 0.4 && maximum(wr) > 0.9 && wr[end] < maximum(wr) - 0.2
      @test all(>(0.0), r.thermal_factor) && maximum(r.thermal_factor) - minimum(r.thermal_factor) > 0.2
      local th = first(o.mrid for o in machines if o.mrid != refm && C.num(o, :p, 0.0) != 0.0)
      @test all(C.num(st.objects[th], :p) == round(C.num(base_store.objects[th], :p) * f; digits = 6) for (st, f) in zip(stores, r.thermal_factor))
      @test all(iszero, r.imbalance)

      # range keys: the default ranges written out give the values of the run
      # without them (only the description in the header differs)
      local spec = TOML.parsefile(roles)
      local dflt = generateCGMESSSHSeries(base, joinpath(tmp, "range_default"); roles = merge(spec, Dict("range" => Dict("load" => [0.6, 1.2], "wind" => [0, 1.3]))), events = events)
      ssh_body(p) = (s = zip_readentry(ZipReader(read(p)), ssh_name, String); s[last(findfirst("</md:FullModel>", s)):end])
      @test all(ssh_body(a) == ssh_body(b) for (a, b) in zip(r.bundles, dflt.bundles))
      # load and wind on narrower ranges over one day in two-hour steps, the
      # thermal factor bounded so that it binds on both sides: no step fails,
      # the rest stays with the slack and every step's bundle says how much
      local narrow = merge(spec, Dict("range" => Dict("load" => [0.8, 1.2], "wind" => [0.5, 1.0], "thermal" => [0.85, 1.15])))
      local n = generateCGMESSSHSeries(base, joinpath(tmp, "range_narrow"); roles = narrow, steps = 12, minutes = 120)
      local nst = [redirect_stdout(() -> C.loadCGMES(b), devnull) for b in n.bundles]
      # every value between base times low and base times high, both rounded
      # as the generator rounds; the load sum spans the range
      function within(sts, cls, lo, hi)
        local ok = true
        for st in sts, m in base_store.byclass[cls]
          local p0 = C.num(base_store.objects[m], :p, 0.0)
          local a, b = round(p0 * lo; digits = 6), round(p0 * hi; digits = 6)
          ok &= min(a, b) <= C.num(st.objects[m], :p, 0.0) <= max(a, b)
        end
        return ok
      end
      @test within(nst, :EnergyConsumer, 0.8, 1.2) && within(nst, :EquivalentInjection, 0.5, 1.0)
      local nlr = [csum(st, :EnergyConsumer) / csum(base_store, :EnergyConsumer) for st in nst]
      @test minimum(nlr) < 0.85 && maximum(nlr) > 1.15
      local hit = findall(!=(0.0), n.imbalance)
      @test all(f -> 0.85 <= f <= 1.15, n.thermal_factor) && any(==(0.85), n.thermal_factor) && any(==(1.15), n.thermal_factor)
      @test hit == findall(f -> f in (0.85, 1.15), n.thermal_factor)
      # the reported rest is what the slack has to carry: the change of the sum of p
      @test maximum(abs.(n.p_sum .- n.base_p .- n.imbalance)) <= 1e-4
      local nman = TOML.parsefile(n.manifest)
      @test nman["generator"]["range"] == Dict("load" => [0.8, 1.2], "wind" => [0.5, 1.0], "thermal" => [0.85, 1.15])
      @test [b["imbalance_mw"] for b in nman["bundle"]] == n.imbalance && [b["thermal_factor"] for b in nman["bundle"]] == n.thermal_factor
      @test findall(b -> any(t -> occursin("of the change stay with the slack", t), b["notes"]), nman["bundle"]) == hit
      @test occursin("ranges load=0.8..1.2, thermal=0.85..1.15, wind=0.5..1,", zip_readentry(ZipReader(read(n.bundles[1])), ssh_name, String))
      # a range on one role key: the equivalent injections, here loads, on
      # their own range, the consumers on the load range with the very values
      # of the run without the key range (one draw per object either way)
      local kr_roles = Dict{String,Any}("EnergyConsumer" => "load", "EquivalentInjection" => "load", "SynchronousMachine:ThermalGeneratingUnit" => "thermal")
      local kr0 = generateCGMESSSHSeries(base, joinpath(tmp, "key_range_off"); roles = kr_roles, steps = 6, minutes = 240)
      local kr = generateCGMESSSHSeries(base, joinpath(tmp, "key_range"); roles = merge(kr_roles, Dict("range" => Dict("EquivalentInjection" => [0.9, 1.1]))), steps = 6, minutes = 240)
      local k0st = [redirect_stdout(() -> C.loadCGMES(b), devnull) for b in kr0.bundles]
      local kst = [redirect_stdout(() -> C.loadCGMES(b), devnull) for b in kr.bundles]
      pvals(sts, cls) = [C.num(st.objects[m], :p, 0.0) for st in sts for m in base_store.byclass[cls]]
      @test pvals(kst, :EnergyConsumer) == pvals(k0st, :EnergyConsumer)
      @test minimum(csum(st, :EnergyConsumer) / csum(base_store, :EnergyConsumer) for st in kst) < 0.7
      @test within(kst, :EquivalentInjection, 0.9, 1.1) && !within(k0st, :EquivalentInjection, 0.9, 1.1)
      @test maximum(abs.(kr.p_sum .- kr.base_p)) <= 1e-4
      @test TOML.parsefile(kr.manifest)["generator"]["range"] == Dict("EquivalentInjection" => [0.9, 1.1])
      @test occursin("ranges EquivalentInjection=0.9..1.1,", zip_readentry(ZipReader(read(kr.bundles[1])), ssh_name, String))
      # a range that reaches no object has no effect, and the notes of every
      # bundle say so: range.thermal without thermal objects (and then no
      # thermal_factor in the manifest), range.load when the one load key has
      # its own range
      local idle = generateCGMESSSHSeries(base, joinpath(tmp, "range_idle"); roles = Dict{String,Any}("EnergyConsumer" => "load", "range" => Dict("thermal" => [0.5, 1.5], "load" => [0.8, 1.2], "EnergyConsumer" => [0.9, 1.1])), steps = 2)
      idle_note(k) = any(t -> occursin("range.$(k) reaches no object", t), idle.notes)
      @test idle_note("thermal") && idle_note("load") && !idle_note("EnergyConsumer")
      @test all(b -> !haskey(b, "thermal_factor") && issubset(idle.notes, b["notes"]), TOML.parsefile(idle.manifest)["bundle"])
    end)() end

    @testset "CGMES SSH series: PowerElectronicsConnection roles" begin (function ()
      # keyed by the class of their units like a machine by its
      # GeneratingUnit: the wind unit on the wind ramp, the PV unit on the
      # load profile with a range of its own, the connection with two unit
      # classes on the plain key (fixed) with a note; the thermal machines
      # balance within their bounds, the connections counted in the sum
      local pbase = joinpath(tmp, "base_pec")
      local wind, pv, mix = _sshs_test_pec_base(base, pbase)
      local proles = Dict{String,Any}(
        "EnergyConsumer" => "load",
        "SynchronousMachine:ThermalGeneratingUnit" => "thermal",
        "PowerElectronicsConnection" => "fixed",
        "PowerElectronicsConnection:PowerElectronicsWindUnit" => "wind",
        "PowerElectronicsConnection:PhotoVoltaicUnit" => "load",
        "range" => Dict("PowerElectronicsConnection:PhotoVoltaicUnit" => [0.9, 1.1], "thermal" => [0.5, 1.5]),
      )
      local pr = generateCGMESSSHSeries(pbase, joinpath(tmp, "series_pec"); roles = proles, steps = 12, minutes = 120)
      local pbs = redirect_stdout(() -> C.loadCGMES(pbase), devnull)
      local psts = [redirect_stdout(() -> C.loadCGMES(b), devnull) for b in pr.bundles]
      @test (pr.role[wind], pr.role[pv], pr.role[mix]) == ("wind", "load", "fixed")
      @test any(n -> occursin("PowerElectronicsConnection $(mix) has units of several classes (BatteryUnit, PhotoVoltaicUnit), the plain key PowerElectronicsConnection applies", n), pr.notes)
      @test !any(n -> occursin("matches no object", n) || occursin("reaches no object", n), pr.notes)
      pp(st, m) = C.num(st.objects[m], :p)
      local wf = [pp(st, wind) / pp(pbs, wind) for st in psts]
      @test minimum(wf) < 0.4 && maximum(wf) > 0.9
      @test all(st -> 0.9 <= pp(st, pv) / pp(pbs, pv) <= 1.1, psts) && length(unique(pp(st, pv) for st in psts)) > 1
      @test all(st -> (pp(st, mix), C.num(st.objects[mix], :q)) == (-5.0, -1.0), psts)
      @test maximum(abs.(pr.p_sum .- pr.base_p .- pr.imbalance)) <= 1e-4
      @test maximum(abs.([sum(pp(st, m) for m in keys(pr.role)) for st in psts] .- pr.p_sum)) <= 1e-9
      @test_throws "unknown unit class WindGeneratingUnit of PowerElectronicsConnection" generateCGMESSSHSeries(pbase, joinpath(tmp, "bad"); roles = Dict("PowerElectronicsConnection:WindGeneratingUnit" => "wind"))
      @test_throws "and PowerElectronicsConnection can be keyed by a unit class" generateCGMESSSHSeries(pbase, joinpath(tmp, "bad"); roles = Dict("EnergyConsumer:BatteryUnit" => "load"))
    end)() end

    @testset "CGMES SSH series: series check and union" begin (function ()
      local inputs = [CGMESSeriesInput(st; label = "step $(k)") for (k, st) in enumerate(stores)]
      local rep = checkCGMESSeries(inputs; config = QstsSeriesConfig())
      @test length(rep.series) == 1 && rep.series[1].snapshots == collect(1:96)
      @test all(p -> p.structural_overlap == 1.0 && p.injection_churn == 0.0 && p.gap_minutes == 15.0, rep.pairs)
      # the switch the events close between two TNs (TP not recomputed) is
      # SSH against TP exactly in the event steps, and the report names it
      @test [s.ssh_tp_conflicts for s in rep.snapshots] == [40 <= k <= 59 ? 1 : 0 for k in 1:96]
      local io = IOBuffer()
      printCGMESSeriesReport(io, rep)
      local txt = String(take!(io))
      @test occursin("snapshot step 40: 1 switch(es) closed in the SSH", txt) && !occursin("snapshot step 39:", txt)
      local cfg = Sparlectra.SparlectraConfig(Dict())
      local u = redirect_stdout(() -> buildCGMESSeriesUnion(inputs; config = cfg), devnull)
      # one union node per TN of the base, boundary TNs included, no _TPL coupling
      @test length(u.node_mrid) == length(C._tpTNs(base_store)) + count(m -> m in base_store.boundary, base_store.byclass[:TopologicalNode])
      @test old_tn in u.node_mrid && new_tn in u.node_mrid
      @test u.n_tp_links == 0 && isempty(u.messages)
      # the switch of the events is the one link, closed in steps 40 to 59
      @test u.link_switches == [sw]
      @test [C.boolval(st.objects[sw], :open) for st in stores] == [!(40 <= k <= 59) for k in 1:96]
      pq_op(op) = op.op === :set && op.field in (:p, :q) && op.target in (:load, :generator, :external_grid)
      @test all(all(op -> pq_op(op) || (op.op === :status && op.target === :link), d.ops) for d in u.deltas)
      @test [count(op -> op.target === :link, d.ops) for d in u.deltas] == [40 <= k <= 59 ? 1 : 0 for k in 1:96]
      @test all(op.id == u.index.id_by_mrid[sw] && op.value == 1.0 for d in u.deltas for op in d.ops if op.target === :link)
      @test all(d -> isempty(d.control_ops) && isempty(d.residual), u.deltas)
      @test isempty(u.deltas[1].ops) && all(k -> any(pq_op, u.deltas[k].ops), 2:96)
      # template plus delta reproduces the direct import of every step
      local failed = Int[]
      for (k, st) in enumerate(stores)
        rk = redirect_stdout(() -> importCGMES(cfg; path = st, name = "step $(k)"), devnull)
        isempty(compareSnapshotReproduction(u, k, rk)) || push!(failed, k)
      end
      @test isempty(failed)

      @testset "mode B on the union" begin (function ()
        # QSTS mode B (sequential_start = previous) over the 96 steps with one
        # tap controller from the configuration (the SmallGrid has no
        # TapChangerControl): transformer B_2WT_220_73_39 (branch 180) holds
        # its 132 kV side (Olive2) at 1.0 pu within 0.005 pu. Mode B starts
        # each step from the step before; the own SV start of mode A does
        # not apply. The series keeps every SSH tap, so
        # the reference of every step is the base step of the transformer
        # and dpos_ref is the simulated position minus that base step. The
        # controller steps on the delivered tap changer (range 0.985 to
        # 1.105 in steps of 0.015, the importer gives every changer its
        # range): ratio 0.985 from 1.0 is step -1.
        local br = u.template.branchVec[180]
        local ctl = _sshs_olive_controllers(u)
        local seq = _sshs_sequential_config(ctl)
        local r = redirect_stdout(() -> runQsts(u; config = seq, parallel_enabled = false), devnull)
        local id = only(r.elements.controller_ids)
        @test r.mode === :sequential && r.origin === :cgmes_series && all(s -> s.status === :solved, r.steps)
        @test endswith(id, "#OLTC")
        local base_pos = round(Int, (br.tap_ratio - 1.0) / br.tap_step)
        @test base_pos == -1
        local pos = qstsValues(r, id, "position")
        local dpos = qstsValues(r, id, "dpos_ref")
        @test all(==(base_pos), qstsValues(r, id, "ref_position")) && dpos == pos .- base_pos
        @test (_qsts_metric_value(r, id, "max_abs_dpos_ref"), _qsts_metric_value(r, id, "rms_dpos_ref")) == (maximum(abs, dpos), sqrt(sum(abs2, dpos) / 96))
        # the switching rows chain the positions: each row starts where the
        # step before ended (step 1 at the base step); the resync row comes
        # first, every control row is counted, no hunting
        local sw = qstsSwitching(r)
        @test all(sw.old[i] == (sw.step[i] == 1 ? base_pos : pos[sw.step[i] - 1]) && sw.new[i] == pos[sw.step[i]] && sw.steps[i] == sw.new[i] - sw.old[i] for i in eachindex(sw.step))
        @test sum(sw.steps) == pos[end] - base_pos
        @test sw.trigger[1] == "resync" && count(==("control"), sw.trigger) == length(sw.step) - 1
        @test (_qsts_metric_value(r, id, "switch_count"), _qsts_metric_value(r, id, "n_hunting")) == (Float64(count(==("control"), sw.trigger)), 0.0)
        local vm = filter(!isnan, qstsValues(r, r.elements.bus_ids[br.toBus], "vm_pu"))
        @test length(vm) == 96 && all(v -> 0.995 - 1e-6 <= v <= 1.005 + 1e-6, vm)
        # each step is compared with its SV as in a snapshot run
        @test all(isfinite, qstsValues(r, "", "sv_max_dvm_pu"))
        # the start is the step before: a step that repeats the step before
        # (snapshots 1, 2, 2) starts at its solution with the carried tap and
        # needs one iteration (from the base solution it needs four)
        local u3 = redirect_stdout(() -> buildCGMESSeriesUnion(inputs[[1, 2, 2]]; config = cfg), devnull)
        local r3 = redirect_stdout(() -> runQsts(u3; config = seq, parallel_enabled = false), devnull)
        @test all(s -> s.status === :solved, r3.steps) && r3.steps[3].iterations == 1 && qstsSwitching(r3).step == [1]
        @test maximum(abs(x - y) for (x, y) in zip(r3.bus_values[:, :, 2], r3.bus_values[:, :, 3]) if !isnan(x); init = 0.0) <= 1e-12
        # two unions in parallel equal the serial run; the vector writer puts
        # each into series_<k> with its position in the manifest
        local two = _sshs_two_unions_check(u, seq, r)
        local dir = mktempdir()
        local files = writeQstsResults(dir, two)
        @test length(files) == 12 && all(f -> startswith(f, joinpath(dir, "series_1")) || startswith(f, joinpath(dir, "series_2")), files)
        local m2 = Sparlectra.scf_json_parse(read(joinpath(dir, "series_2", "qsts_manifest.json"), String))
        @test (m2["series"], m2["n_series"], m2["mode"]) == (2, 2, "sequential")
      end)() end

      @testset "event list on the union" begin (function ()
        # 12 steps on the union template: L1 (B_ACL_132_6_7) opens in step
        # 3, which makes L2 (B_ACL_132_69_41) a bridge; L2 opens in step 6 and
        # cuts off three buses with load and no source, and closes in step 9.
        # The sourceless island is isolated (F2): no values, n_isolated and
        # shed_load_mw, and n_islands counts it from step 6 to 8.
        local l1 = "_045d71d4-c766-11e1-8775-005056c00008"
        local l2 = "_04607f15-c766-11e1-8775-005056c00008"
        local j1 = u.index.internal_by_id[u.index.id_by_mrid[l1]]
        local j2 = u.index.internal_by_id[u.index.id_by_mrid[l2]]
        local cut = deepcopy(u.template)
        cut.branchVec[j1].status = 0
        @test !(j2 in Sparlectra._scenario_bridge_branches(u.template)) && j2 in Sparlectra._scenario_bridge_branches(cut)
        local f = joinpath(mktempdir(), "events.csv")
        write(f, "step,element_id,action,value\n3,$(l1),open,\n6,$(l2),open,\n9,$(l2),close,\n")
        local ev = readQstsEvents(f)
        local es = qstsEventSeries(u, ev; n_steps = 12)
        @test length(es.ops) == 12 && isempty(es.ops[2]) && length(es.ops[3]) == 3 && length(es.ops[6]) == 6 && es.ops[9] == es.ops[3]
        local seq = Sparlectra.SparlectraConfig(Dict{String,Any}("qsts" => Dict{String,Any}("mode" => "sequential")))
        local r = redirect_stdout(() -> runQsts(es; config = seq), devnull)
        @test all(s -> s.status === :solved, r.steps) && r.base_step == 0
        local islands = qstsValues(r, "", "n_islands")
        @test islands[6] == islands[5] + 1 && islands[9] == islands[8] - 1 && islands[6:8] == fill(islands[6], 3)
        local isolated = qstsValues(r, "", "n_isolated")
        @test all(>(0), isolated[6:8]) && all(isolated[k] == isolated[5] for k in (9, 10)) && all(>(0), qstsValues(r, "", "shed_load_mw")[6:8])
        @test isolated[6] - isolated[5] == 3.0
        # the cut-off buses have no values from step 6 to 8 and again from 9
        local far = u.template.branchVec[j2]
        local ends = [far.fromBus, far.toBus]
        local vm(b) = qstsValues(r, r.elements.bus_ids[b], "vm_pu")
        @test any(b -> all(isnan, vm(b)[6:8]) && all(isfinite, vm(b)[9:12]), ends)
        local manifest = Sparlectra.scf_json_parse(read(last(writeQstsResults(mktempdir(), r)), String))
        @test manifest["event_list"] == Dict{String,Any}("source" => "events.csv", "rows" => 3)
        # an event series has no snapshots: a snapshot run and the snapshot
        # start are refused with the key named
        @test_throws "qsts.mode = sequential only" runQsts(es; config = Sparlectra.SparlectraConfig(Dict{String,Any}("qsts" => Dict{String,Any}("mode" => "snapshots"))))
        @test_throws "qsts.sequential_start = snapshot" runQsts(es; config = Sparlectra.SparlectraConfig(Dict{String,Any}("qsts" => Dict{String,Any}("mode" => "sequential", "sequential_start" => "snapshot"))))
        # a switch inside one node is no element of the union
        write(f, "step,element_id,action,value\n2,no-such-switch,open,\n")
        @test_throws "line 2: unknown element_id \"no-such-switch\" (a switch must be a link of the union" qstsEventSeries(u, readQstsEvents(f))
      end)() end
    end)() end

    @testset "CGMES SSH series: zip input" begin (function ()
      # the two fixture zips as they are, the boundary zip given twice
      local z = generateCGMESSSHSeries([src; src[2]], joinpath(tmp, "zips"); roles = roles, steps = 2)
      local stem = splitext(basename(src[1]))[1]
      @test basename.(z.bundles) == ["$(stem)_20300102T0900Z.zip", "$(stem)_20300102T0915Z.zip"]
      # entry name => (zip, bytes) of the fixture zips
      local orig = Dict{String,Tuple{String,Vector{UInt8}}}()
      for p in src
        local zr = ZipReader(read(p))
        for n in zip_names(zr)
          orig[n] = (p, zip_readentry(zr, n))
        end
      end
      local same_bytes = true
      for b in z.bundles
        local zb = ZipReader(read(b))
        same_bytes &= sort(zip_names(zb)) == sort!(collect(keys(orig)))
        for n in zip_names(zb)
          n == ssh_name || (same_bytes &= zip_readentry(zb, n) == orig[n][2])
        end
      end
      @test same_bytes
      local man = TOML.parsefile(z.manifest)
      @test all(length(b["files"]) == length(orig) && isempty(b["missing_dependencies"]) for b in man["bundle"])
      @test all(f["source"] == string(orig[f["file"]][1], "!", f["file"]) for b in man["bundle"] for f in b["files"])
      # urn:uuid:<id> in the SSH, the bare id in the manifest
      local ids = [C.readCGMESHeader(only(f for f in C.collectCGMESFiles(b) if f.name == ssh_name)).model_id for b in z.bundles]
      @test [only(f["uuid"] for f in b["files"] if f["kind"] == "SSH") for b in man["bundle"]] == ids
      # the second boundary copy is left out and noted in every bundle
      local bd = [n for (n, (p, _)) in orig if p == src[2]]
      local copies = filter(n -> occursin("left out", n), z.notes)
      @test length(copies) == length(bd) == 2 && all(n -> any(e -> occursin("$(src[2])!$(e)", n), bd), copies)
      @test all(b -> issubset(z.notes, b["notes"]), man["bundle"])
    end)() end

    @testset "CGMES SSH series: input errors" begin (function ()
      local bad = joinpath(tmp, "bad_events.csv")
      write(bad, "step,element_id,action,value\n3,$(sw),toggle,\n")
      @test_throws "line 2: action \"toggle\"" generateCGMESSSHSeries(base, joinpath(tmp, "bad"); roles = roles, events = bad)
      local consumer = first(base_store.byclass[:EnergyConsumer])
      write(bad, "5,$(consumer),open,\n")
      @test_throws "line 1: $(consumer) is not a switch" generateCGMESSSHSeries(base, joinpath(tmp, "bad"); roles = roles, events = bad)
      @test_throws "unknown class Generator" generateCGMESSSHSeries(base, joinpath(tmp, "bad"); roles = Dict("Generator" => "load"))
      # the three equivalent injections cannot take up the load change
      @test_throws "the balance needs a thermal factor of" generateCGMESSSHSeries(base, joinpath(tmp, "bad"); roles = Dict("EnergyConsumer" => "load", "EquivalentInjection" => "thermal"))
      @test_throws "range.load = [1.2, 0.8], expected [low, high]" generateCGMESSSHSeries(base, joinpath(tmp, "bad"); roles = Dict("EnergyConsumer" => "load", "range" => Dict("load" => [1.2, 0.8])))
      # a key range needs a key of the mapping with role load or wind
      @test_throws "range.EnergySource is neither a profile" generateCGMESSSHSeries(base, joinpath(tmp, "bad"); roles = Dict("EnergyConsumer" => "load", "range" => Dict("EnergySource" => [0.8, 1.2])))
      @test_throws "range.EquivalentInjection: the key has role thermal" generateCGMESSSHSeries(base, joinpath(tmp, "bad"); roles = Dict("EnergyConsumer" => "load", "EquivalentInjection" => "thermal", "range" => Dict("EquivalentInjection" => [0.8, 1.2])))
      # one file name, two contents
      local other = joinpath(tmp, "other")
      local tpbd = only(n for n in readdir(base) if occursin("_TP_BD_", n))
      mkpath(other)
      write(joinpath(other, tpbd), read(joinpath(base, tpbd), String) * "\n")
      @test_throws "two different files named $(tpbd)" generateCGMESSSHSeries([base, other], joinpath(tmp, "bad"); roles = roles)
      # every check runs before the output directory exists
      @test !isdir(joinpath(tmp, "bad"))
    end)() end
  end
end

# Mode B on two unions on worker threads against the serial run, alone: the
# threads profile (tools/run_gates.sh threads) repeats it with four threads.
function run_cgmes_ssh_series_parallel_tests()
  @testset "CGMES SSH series: two unions in parallel" begin (function ()
    mktempdir() do tmp
      local fx = _sshs_test_series(tmp)
      local inputs = [CGMESSeriesInput(st; label = "step $(k)") for (k, st) in enumerate(fx.stores)]
      local u = redirect_stdout(() -> buildCGMESSeriesUnion(inputs; config = Sparlectra.SparlectraConfig(Dict())), devnull)
      local seq = _sshs_sequential_config(_sshs_olive_controllers(u))
      local r = redirect_stdout(() -> runQsts(u; config = seq, parallel_enabled = false), devnull)
      @test all(s -> s.status === :solved, r.steps)
      _sshs_two_unions_check(u, seq, r)
    end
  end)() end
end
