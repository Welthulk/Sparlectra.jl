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

# file: test/test_legacy_powsybl_importer.jl
# purpose: legacy profile: the testsets of test/test_powsybl_importer.jl that pin an earlier
#          change (an issue, a bug, a bit-identity or off-path guarantee).
#          They moved here unchanged from test_powsybl_importer.jl; the helpers they use
#          are defined in that file, which the legacy group includes first.
using Sparlectra
using Test

function run_legacy_powsybl_importer_tests()
  # cut from testset: PowSyBl importer
    @testset "bundle round trip (every fixture)" begin (function ()
      for case in _POWSYBL_CASES
        tables = read_powsybl_bundle(_powsybl_fixture(case))
        mktempdir() do dir
          out = joinpath(dir, case * ".powsybl")
          write_powsybl_bundle(tables, out)
          back = read_powsybl_bundle(out)
          for name in Sparlectra.POWSYBL_TABLE_NAMES
            a = Sparlectra.powsybl_table(tables, name)
            b = Sparlectra.powsybl_table(back, name)
            @test keys(a) == keys(b)
            for c in keys(a)
              @test isequal(a[c], b[c])
              @test eltype(a[c]) == eltype(b[c])
            end
          end
          @test back.manifest["case"] == tables.manifest["case"]
          @test Sparlectra.detect(PowsyblAdapter, out)
        end
      end
    end)() end

    @testset "builder structure (every fixture)" begin (function ()
      for case in _POWSYBL_CASES
        tables = read_powsybl_bundle(_powsybl_fixture(case))
        net, report = build_net_from_powsybl(tables, PowsyblAdapterOptions())
        bb = tables.bus_breaker_view_buses
        sw = tables.switches
        skipped_buses = count(s -> s.kind == "bus", report.skipped)
        # an auxiliary bus (star point, boundary, X-node) carries the AuxBus component type
        real_buses = count(nd -> !Sparlectra._is_aux_node(nd), net.nodeVec)
        @test real_buses == length(bb.id) - skipped_buses
        retained = count(identity, sw.retained)
        @test length(net.linkVec) == retained - count(s -> s.kind == "switch", report.skipped)
        # no branch-derived bus shunt (0.20.0): the shunt list holds the
        # shunt compensators only (the SVC is a prosumer)
        @test length(net.shuntVec) == length(tables.shunt_compensators.id)
        @test report.counts["switch"].read == length(sw.id)
        # no retained switch is a branch: every branch joins two buses through
        # an impedance element, links are the only switch representation
        @test length(net.branchVec) == report.counts["line"].built + report.counts["2wt"].built + 3 * report.counts["3wt"].built + report.counts["dangling_line"].built + 2 * report.counts["tie_line"].built
        # islands: one per distinct non-negative synchronous component; the
        # net is raw (retained switches are links, not yet merged), so the
        # electrical view over branches AND closed links is the right count
        comps = Set(c for c in tables.buses.synchronous_component if c >= 0)
        islands = Sparlectra.electricalIslandComponents(net)
        @test length(islands) == length(comps)
        # every paired dangling line is skipped into its tie line
        dl = tables.dangling_lines
        for i in eachindex(dl.id)
          dl.paired[i] || continue
          entry = findfirst(s -> s.id == dl.id[i] && s.kind == "dangling_line", report.skipped)
          @test entry !== nothing
          entry === nothing || @test occursin("tie line", report.skipped[entry].reason)
        end
        # every skipped element carries id, kind and a reason
        @test all(s -> !isempty(s.id) && !isempty(s.kind) && !isempty(s.reason), report.skipped)
        # bus-view id in the external-id channel
        @test net.nodeVec[net.busDict[bb.id[1]]].comp.cID == bb.bus_id[1]
      end
    end)() end

    @testset "patched fixtures (four_substations)" begin (function ()
      base = read_powsybl_bundle(_powsybl_fixture("four_substations"))
      # a line open at side 1 is a one-sided open branch
      t = read_powsybl_bundle(_powsybl_fixture("four_substations"))
      t.lines.connected1[1] = false
      net, report = build_net_from_powsybl(t, PowsyblAdapterOptions())
      f = net.busDict[t.lines.bus_breaker_bus1_id[1]]
      to = net.busDict[t.lines.bus_breaker_bus2_id[1]]
      br = only(b for b in net.branchVec if b.fromBus == f && b.toBus == to)
      @test Sparlectra._branch_terminal_state(br) == :open_from
      # a retained switch set open is an open link
      t2 = read_powsybl_bundle(_powsybl_fixture("four_substations"))
      k = findfirst(identity, t2.switches.retained)
      t2.switches.open[k] = true
      net2, _ = build_net_from_powsybl(t2, PowsyblAdapterOptions())
      @test count(l -> l.status == 0, net2.linkVec) == 1
      @test length(net2.linkVec) == length(net.linkVec)
      # a generator regulating another bus: held local with a report entry, or PQ
      # GH1 is not the slack of its component (GH2 is), so its bus type moves
      t3 = read_powsybl_bundle(_powsybl_fixture("four_substations"))
      g = t3.generators
      i = findfirst(==("GH1"), g.id)
      other = "S2VL1_0"
      g.regulated_bus_id[i] = other
      g.regulated_bus_breaker_bus_id[i] = other
      net3, report3 = build_net_from_powsybl(t3, PowsyblAdapterOptions(remote_regulation = :hold_local))
      @test any(m -> occursin("GH1", m) && occursin(other, m) && occursin(g.bus_id[i], m), report3.messages)
      @test net3.nodeVec[net3.busDict[g.bus_breaker_bus_id[i]]]._nodeType == Sparlectra.PV
      net4, report4 = build_net_from_powsybl(t3, PowsyblAdapterOptions(remote_regulation = :pq))
      @test any(m -> occursin("GH1", m) && occursin("PQ", m), report4.messages)
      @test net4.nodeVec[net4.busDict[g.bus_breaker_bus_id[i]]]._nodeType == Sparlectra.PQ
      # the unpatched build carries no remote notice
      _, report0 = build_net_from_powsybl(base, PowsyblAdapterOptions())
      @test !any(m -> occursin("regulates bus", m), report0.messages)
    end)() end

    # block E.2 (0.20.0): the MATPOWER export of micro_grid_be moves the
    # terminal excess of every asymmetric branch shunt into the bus GS/BS,
    # named per branch; the reimport keeps them as parts of the bus shunt,
    # reproduces the Y-bus, and an outage of the branch removes its part
    @testset "MATPOWER export rule and read-back (micro_grid_be)" begin (function ()
      tables = read_powsybl_bundle(_powsybl_fixture("micro_grid_be"))
      net, _ = build_net_from_powsybl(tables, PowsyblAdapterOptions(remote_regulation = :remote))
      mktempdir() do dir
        mfile = joinpath(dir, "micro_grid_be.m")
        writeMatpowerCasefile(net, mfile; write_solution = false)
        txt = read(mfile, String)
        @test occursin("mpc.sparlectra.branch_shunts", txt)
        @test occursin("% branch-derived shunt of branch", txt)
        @test_throws ArgumentError writeMatpowerCasefile(net, joinpath(dir, "x.m"); write_solution = false, asymmetric_shunts = :drop)
        cfg = load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH; reload = true)
        back = Sparlectra.import_case(mfile, cfg).net
        Y1 = createYBUS(net = net, sparse = true)
        Y2 = createYBUS(net = back, sparse = true)
        @test size(Y1) == size(Y2)
        @test maximum(abs.(Y1 - Y2)) <= 1e-12
        # lines stay lines: the two 225 kV lines 4-6 carry a symmetric g,
        # which a MATPOWER line row cannot hold; it travels as a branch part
        # of the bus shunt. They used to be written as "T" and came back as
        # transformers at ratio 1 (5 lines and 8 transformers instead of 7/6)
        @test (length(back.linesAC), length(back.trafos)) == (length(net.linesAC), length(net.trafos))
        parts = sum(length(sh.branch_parts) for sh in back.shuntVec; init = 0)
        @test parts > 0
        # the outage of a branch with a part takes the part out of the bus shunt
        sh = first(sh for sh in back.shuntVec if !isempty(sh.branch_parts))
        bidx, part = first(sh.branch_parts)
        work = deepcopy(back)
        wsh = work.shuntVec[findfirst(s -> s.busIdx == sh.busIdx, work.shuntVec)]
        y0 = wsh.y_pu_shunt
        Sparlectra._remove_branch_shunt_parts!(work, bidx)
        @test isapprox(wsh.y_pu_shunt, y0 - part; atol = 1e-15)
        name = getCompName(back.branchVec[bidx].comp)
        res = runContingencies!(back, [ContingencyCase(name, :branch, name)]; parallel_enabled = false)
        @test length(res) == 1
        # D2 of #355: a scenario on this reimport (the service converts a
        # MATPOWER case to SCF for its scenario sources and runs the imported
        # net) removes the part on the patch path too, and restores it; the
        # patch path (forced by a neutral op) equals the outage path
        mcase = Sparlectra.convert_case(MatpowerAdapter(), Sparlectra.MatpowerIO.read_case(mfile; legacy_compat = true), Sparlectra.matpower_adapter_options(cfg))
        midx = ScenarioIndex(mcase)
        bid = only(id for (id, k) in midx.internal_by_id if k == bidx && midx.kind_by_id[id] in (:branch, :transformer))
        op_out = PatchOp(op = :status, target = midx.kind_by_id[bid], id = bid, value = 0.0)
        work2 = deepcopy(back)
        wsh2 = work2.shuntVec[findfirst(s -> s.busIdx == sh.busIdx, work2.shuntVec)]
        undo = apply!(work2, [op_out], midx)
        @test isapprox(wsh2.y_pu_shunt, y0 - part; atol = 1e-15)
        restore!(work2, undo)
        @test wsh2.y_pu_shunt == y0
        neutral = first(id for (id, k) in midx.kind_by_id if k === :load)
        rs = redirect_stdout(devnull) do
          runScenarios!(back, [Scenario(name = "outage", ops = [op_out]), Scenario(name = "patched", ops = [op_out, PatchOp(op = :scale, target = :load, id = neutral, factor = 1.0)])]; index = midx, parallel_enabled = false)
        end
        @test all(f -> f in (:name, :kind, :elements) || isequal(getproperty(rs[1], f), getproperty(rs[2], f)), fieldnames(ContingencyResult))
        # the CGMES export follows the same rule: a segment has one charging
        # value, so the five dangling lines write their one-sided admittance
        # as a compensator at the network bus, named per line, and the
        # delivery solves to the voltages of the source (the total on the
        # segment, read back in halves, was 3e-3 pu off). The CGMES import
        # puts half of a transformer's end admittance on each terminal
        # (0.41.0), so the two transformers with their magnetizing arm on
        # side 1 write it as a compensator too
        _, erg = runpf!(net, 60, 1e-10, 0; method = :rectangular, qlimits_enabled = false)
        @test erg == 0
        notes = String[]
        delivery = writeCGMESFiles(net; path = mkpath(joinpath(dir, "cgmes")), created = Dates.DateTime(2026, 1, 1), notices = notes, zip = true)
        @test count(n -> startswith(n, "line ") && occursin("the excess as a LinearShuntCompensator", n), notes) == 5
        @test count(n -> startswith(n, "transformer ") && occursin("the excess as a LinearShuntCompensator", n), notes) == 2
        cim = Sparlectra.importCGMES(path = delivery[end], require_boundary = false).net
        @test length(cim.shuntVec) == length(net.shuntVec) + 7
        _, erg_cim = runpf!(cim, 60, 1e-10, 0; method = :rectangular, qlimits_enabled = false)
        @test erg_cim == 0
        vm_of(n) = (names = Dict(idx => name for (name, idx) in n.busDict); Dict(names[nd.busIdx] => nd._vm_pu for nd in n.nodeVec))
        source_vm, delivery_vm = vm_of(net), vm_of(cim)
        @test keys(source_vm) == keys(delivery_vm)
        @test maximum(abs(source_vm[k] - delivery_vm[k]) for k in keys(source_vm)) <= 1e-10
      end
    end)() end

    @testset "native IIDM reader against the bundles (every fixture)" begin (function ()
      # the state columns hold the file's state in the reader and
      # OpenLoadFlow's solution in a bundle; every other schema column must
      # agree, rows joined by the index columns
      state = Set([:p, :q, :i, :p1, :q1, :i1, :p2, :q2, :i2, :p3, :q3, :i3, :v_mag, :v_angle, :v, :angle, :boundary_p, :boundary_q, :boundary_i, :boundary_v_mag, :boundary_v_angle, :min_q_at_p, :max_q_at_p, :solved_tap_position, :solved_section_count])
      same(a::AbstractString, b::AbstractString) = a == b
      same(a::Float64, b::Float64) = (isnan(a) && isnan(b)) || a == b || (isfinite(a) && isfinite(b) && abs(a - b) <= 1e-9 * max(1.0, abs(a), abs(b)))
      same(a, b) = a == b
      for case in _POWSYBL_CASES
        dir = _powsybl_fixture(case)
        bundle = read_powsybl_bundle(dir)
        native = Sparlectra.read_iidm_tables(joinpath(dir, case * ".xiidm"))
        differing = String[]
        for name in Sparlectra.POWSYBL_TABLE_NAMES
          b = Sparlectra.powsybl_table(bundle, name)
          n = Sparlectra.powsybl_table(native, name)
          schema = Sparlectra.POWSYBL_SCHEMA[name]
          idx = [c.name for c in schema if c.index]
          key(t, i) = Tuple(t[c][i] for c in idx)
          nb = Sparlectra.powsybl_table_rows(b)
          nn = Sparlectra.powsybl_table_rows(n)
          nb == nn || push!(differing, "$(name): $(nb) rows in the bundle, $(nn) native")
          rows = Dict(key(n, j) => j for j in 1:nn)
          for i in 1:nb
            j = get(rows, key(b, i), nothing)
            if j === nothing
              push!(differing, "$(name): row $(key(b, i)) missing")
              continue
            end
            for c in schema
              (haskey(b, c.name) && !(c.name in state)) || continue
              same(b[c.name][i], n[c.name][j]) || push!(differing, "$(name).$(c.name) row $(key(b, i)): bundle $(repr(b[c.name][i])), native $(repr(n[c.name][j]))")
            end
          end
        end
        @test isempty(differing)
        isempty(differing) || println("      powsybl native reader $(case): ", join(first(differing, 5), "; "))
      end
      # micro_grid_be's file state was solved at other tap positions: the
      # reader drops it and the report says so; ieee14's state seeds the start
      stale = Sparlectra.read_iidm_tables(joinpath(_powsybl_fixture("micro_grid_be"), "micro_grid_be.xiidm"))
      @test all(isnan, stale.bus_breaker_view_buses.v_mag)
      @test occursin("starts flat", stale.manifest["state_note"])
      _, stale_report = build_net_from_powsybl(stale, PowsyblAdapterOptions())
      @test any(n -> occursin("starts flat", n), stale_report.messages)
      fresh = Sparlectra.read_iidm_tables(joinpath(_powsybl_fixture("ieee14"), "ieee14.xiidm"))
      @test all(isfinite, fresh.bus_breaker_view_buses.v_mag)
      @test isempty(fresh.manifest["state_note"])
      # refused constructs name themselves
      err = @test_throws ArgumentError Sparlectra.read_iidm_tables(joinpath(_powsybl_fixture("ieee14"), "ieee14.xiidm.bz2"))
      @test occursin("unpack", err.value.msg)
      mktempdir() do dir
        bad = joinpath(dir, "notiidm.xiidm")
        write(bad, "<?xml version=\"1.0\"?><root/>")
        bad_err = @test_throws ArgumentError Sparlectra.read_iidm_tables(bad)
        @test occursin("not <network>", bad_err.value.msg)
      end
    end)() end

  return true
end
