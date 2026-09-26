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

# file: test/test_powsybl_importer.jl
# purpose: the PowSyBl adapter on the five shipped bundles under
#          test/fixtures/powsybl: bundle reader failures and round trip,
#          builder structure, the power flow against the OpenLoadFlow
#          reference voltages, the Y-bus against MATPOWER case14, branch
#          flows against the OLF columns, patched fixtures (open terminal,
#          open switch, remote regulation), the slack override, the Web UI
#          wiring (format hint, option section, help topics, the shipped
#          bundles in the selector), the N-1 outage of the PST, the MATPOWER
#          export rule with its read-back, the SE self-test, the native
#          IIDM reader against every bundle with the power flow from the
#          .xiidm file, and the import_case dispatch on both sources.
#          Profile `adapters`.

using Sparlectra
using SparlectraApp
using Test
using DelimitedFiles
using LinearAlgebra
using Random

const _POWSYBL_FIXTURES = joinpath(@__DIR__, "fixtures", "powsybl")
const _POWSYBL_CASES = ("ieee14", "ieee57", "four_substations", "micro_grid_be", "eurostag_tie_lines")

_powsybl_fixture(case::AbstractString) = joinpath(_POWSYBL_FIXTURES, case * ".powsybl")

# reference_buses.csv of a bundle: bus-breaker id => (v_pu, v_angle_deg, component, bus-view id)
function _powsybl_reference(dir::AbstractString)
  data, header = readdlm(joinpath(dir, "reference_buses.csv"), ',', String; quotes = true, header = true)
  cols = Dict(Symbol(strip(h)) => k for (k, h) in enumerate(vec(header)))
  rows = Dict{String,NamedTuple{(:v_pu, :va, :comp, :bus_view),Tuple{Float64,Float64,Int,String}}}()
  for i in axes(data, 1)
    rows[String(data[i, cols[:bus_breaker_id]])] = (v_pu = parse(Float64, data[i, cols[:v_pu]]), va = parse(Float64, data[i, cols[:v_angle_deg]]), comp = parse(Int, data[i, cols[:synchronous_component]]), bus_view = String(data[i, cols[:id]]))
  end
  return rows
end

# The OLF reference bus per synchronous component, from the manifest.
function _powsybl_reference_buses(tables)
  refs = Dict{Int,String}()
  for c in get(get(tables.manifest, "reference", Dict()), "components", [])
    refs[Int(c["synchronous_component_num"])] = String(c["reference_bus_id"])
  end
  return refs
end

# Largest voltage and angle deviation of a solved net against the OLF
# reference; angles relative to the OLF reference bus of each component.
function _powsybl_voltage_deviation(net, tables, dir)
  ref = _powsybl_reference(dir)
  refbuses = _powsybl_reference_buses(tables)
  offsets = Dict{Int,Float64}()
  for (bus, r) in ref
    haskey(net.busDict, bus) || continue
    r.bus_view == get(refbuses, r.comp, "") && (offsets[r.comp] = net.nodeVec[net.busDict[bus]]._va_deg - r.va)
  end
  max_dv = 0.0
  max_da = 0.0
  worst = ""
  for (bus, r) in ref
    haskey(net.busDict, bus) || continue
    node = net.nodeVec[net.busDict[bus]]
    dv = abs(node._vm_pu - r.v_pu)
    da = abs(node._va_deg - r.va - get(offsets, r.comp, 0.0))
    dv > max_dv && (max_dv = dv; worst = bus)
    max_da = max(max_da, da)
  end
  return (max_dv = max_dv, max_da = max_da, worst = worst)
end

# The run configuration that reproduces OpenLoadFlow's defaults: every
# synchronous component solved, the P imbalance shared the way OLF shares
# it (participation factors of the import, proportional to max_p over the
# units with a nonzero target), and a limit switch without hysteresis.
function _powsybl_olf_config(; hysteresis_pu::Float64 = 1e-6)
  pf = PowerFlowConfig(max_iter = 60, tol = 1e-9, islands_enabled = true, qlimits = Sparlectra.QLimitConfig(hysteresis_pu = hysteresis_pu), distributed_slack = Sparlectra.DistributedSlackConfig(enabled = true, p_mode = :imported))
  return SparlectraConfig(powerflow = pf, output = OutputConfig(logfile_results = :off, startup_latency_hint = false), control = ControlConfig())
end

function run_powsybl_importer_tests()
  @testset "PowSyBl importer" begin (function ()
    @testset "bundle reader failures (ieee14)" begin (function ()
      mktempdir() do dir
        bad = joinpath(dir, "ieee14.powsybl")
        cp(_powsybl_fixture("ieee14"), bad)
        vl = joinpath(bad, "voltage_levels.csv")
        lines = readlines(vl)
        lines[1] = replace(lines[1], "\"nominal_v\"" => "\"nominal_kv\"")
        write(vl, join(lines, "\n") * "\n")
        err = try
          read_powsybl_bundle(bad)
          nothing
        catch e
          e
        end
        @test err isa ArgumentError
        @test occursin("voltage_levels", err.msg) && occursin("nominal_v", err.msg)
        # the version check comes before any table
        m = joinpath(bad, "manifest.json")
        write(m, replace(read(m, String), "\"format_version\": 1" => "\"format_version\": 2"))
        err2 = try
          read_powsybl_bundle(bad)
          nothing
        catch e
          e
        end
        @test err2 isa ArgumentError
        @test occursin("format_version 2", err2.msg)
        # detection reads the manifest marker, not the directory name alone
        @test !Sparlectra.detect(PowsyblAdapter, dir)
      end
    end)() end

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
        # an auxiliary bus (star point, boundary, X-node) carries the Aux_ component name
        real_buses = count(nd -> !startswith(getCompName(nd.comp), "Aux_"), net.nodeVec)
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
        @test occursin("PowSyBl import", format_powsybl_report(report))
        # bus-view id in the external-id channel
        @test net.nodeVec[net.busDict[bb.id[1]]].comp.cID == bb.bus_id[1]
      end
    end)() end

    @testset "power flow against the OpenLoadFlow reference (every fixture)" begin (function ()
      bands = Dict("ieee14" => (1e-6, 1e-4), "ieee57" => (1e-6, 1e-4), "four_substations" => (1e-5, 1e-3), "micro_grid_be" => (1e-5, 1e-3), "eurostag_tie_lines" => (1e-5, 1e-3))
      for case in _POWSYBL_CASES
        dir = _powsybl_fixture(case)
        tables = read_powsybl_bundle(dir)
        # remote regulation as OLF does it (voltageRemoteControl on by default)
        opts = case == "micro_grid_be" ? PowsyblAdapterOptions(remote_regulation = :remote) : PowsyblAdapterOptions()
        net, report = build_net_from_powsybl(tables, opts)
        res = run_sparlectra(net = net, config = _powsybl_olf_config())
        @test res.numerical_converged
        dev = _powsybl_voltage_deviation(net, tables, dir)
        dv_band, da_band = bands[case]
        @test dev.max_dv <= dv_band
        @test dev.max_da <= da_band
        println("      powsybl $(case): max |dV| $(round(dev.max_dv; sigdigits = 3)) pu, max |dtheta| $(round(dev.max_da; sigdigits = 3)) deg (bands $(dv_band) / $(da_band), worst bus $(dev.worst))")
        # the same network converges from a flat start with the default Q-limit hysteresis
        flat, _ = build_net_from_powsybl(tables, opts)
        for nd in flat.nodeVec
          Sparlectra.setVmVa!(node = nd, vm_pu = 1.0, va_deg = 0.0)
        end
        flat_cfg = SparlectraConfig(powerflow = PowerFlowConfig(max_iter = 60, tol = 1e-8, islands_enabled = true), output = OutputConfig(logfile_results = :off, startup_latency_hint = false), control = ControlConfig())
        flat_res = run_sparlectra(net = flat, config = flat_cfg)
        @test flat_res.numerical_converged
      end
    end)() end

    @testset "Y-bus identity ieee14 against MATPOWER case14" begin (function ()
      case14 = joinpath(Sparlectra.large_cases_dir(), "case14.m")
      if !isfile(case14)
        println("      powsybl Y-bus identity: SKIPPED (case14.m not in the large-case directory $(Sparlectra.large_cases_dir()))")
      else
        println("      powsybl Y-bus identity: RAN on $(case14)")
        tables = read_powsybl_bundle(_powsybl_fixture("ieee14"))
        net, _ = build_net_from_powsybl(tables, PowsyblAdapterOptions())
        mp = createNetFromMatPowerFile(filename = case14)
        # bus-breaker ids B1..B14 in bus order, MATPOWER buses 1..14 in bus order
        order = sortperm([parse(Int, match(r"(\d+)$", b).match) for b in tables.bus_breaker_view_buses.id])
        Y = Matrix(createYBUS(net = net))[order, order]
        Ym = Matrix(createYBUS(net = mp))
        # branches 7-8 (14/20 kV) and 7-9 (14/12 kV) carry tap 1 in MATPOWER
        # and are lines between voltage levels in PowSyBl's IEEE-CDF import,
        # with large opposite b1/b2 that turn the physical conductor into the
        # nominal-ratio transformer; the per-terminal arms of the branch
        # reproduce them without any bus shunt, so the identity holds
        # everywhere
        for i in 1:14, j in 1:14
          @test abs(Y[i, j] - Ym[i, j]) <= 1e-9
        end
      end
    end)() end

    @testset "branch flows against the OLF columns (ieee14)" begin (function ()
      dir = _powsybl_fixture("ieee14")
      tables = read_powsybl_bundle(dir)
      net, _ = build_net_from_powsybl(tables, PowsyblAdapterOptions())
      res = run_sparlectra(net = net, config = _powsybl_olf_config())
      @test res.numerical_converged
      ln = tables.lines
      checked = 0
      for i in eachindex(ln.id)
        f = net.busDict[ln.bus_breaker_bus1_id[i]]
        t = net.busDict[ln.bus_breaker_bus2_id[i]]
        br = findfirst(b -> (b.fromBus == f && b.toBus == t) || (b.fromBus == t && b.toBus == f), net.branchVec)
        br === nothing && continue
        flow = getBranchFlow(net.branchVec[br], net.nodeVec[f], net.nodeVec[t])
        (flow === nothing || flow.pFlow === nothing) && continue
        # p1 is positive from the side-1 bus into the branch, so is
        # Sparlectra's from-side flow
        @test abs(flow.pFlow - ln.p1[i]) <= 1e-3
        # the from arm of the branch carries b1 entirely (0.20.0), so q1
        # compares without a bus-shunt correction
        @test abs(flow.qFlow - ln.q1[i]) <= 1e-3
        checked += 1
      end
      @test checked == length(ln.id)
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

    @testset "slack override (ieee14)" begin (function ()
      tables = read_powsybl_bundle(_powsybl_fixture("ieee14"))
      net, report = build_net_from_powsybl(tables, PowsyblAdapterOptions(slack_ids = ["B2-G"]))
      @test length(report.slack) == 1
      @test report.slack[1].generator == "B2-G"
      @test report.slack[1].reason == "override"
      @test net.slackVec == [net.busDict["B2"]]
      _, plain = build_net_from_powsybl(tables, PowsyblAdapterOptions())
      @test plain.slack[1].generator == "B1-G"
      @test plain.slack[1].reason != "override"
      # an unknown id is a notice, the default choice stays
      _, unknown = build_net_from_powsybl(tables, PowsyblAdapterOptions(slack_ids = ["NOPE"]))
      @test unknown.slack[1].generator == "B1-G"
      @test any(m -> occursin("NOPE", m), unknown.messages)
    end)() end

    # the Web UI knows the format: the format hint of a bundle, the option
    # section derived from PowsyblAdapterOptions, and the help topics that
    # point at the configuration section of the docs page
    @testset "Web UI wiring (ieee14)" begin (function ()
      bundle = _powsybl_fixture("ieee14")
      @test SparlectraApp._webui_case_format_hint(bundle) == :powsybl
      @test SparlectraApp._webui_case_format_hint(joinpath(bundle, "ieee14.xiidm")) == :powsybl
      html = SparlectraApp._webui_adapter_options_html(:powsybl, Dict{String,Any}())
      for field in ("powsybl_import_hvdc_mode", "powsybl_import_remote_regulation", "powsybl_import_multi_slack", "powsybl_import_slack_ids", "powsybl_import_base_mva")
        @test occursin("name=\"$(field)\"", html)
        @test SparlectraApp.WEBUI_FORM_HELP_TOPICS[field] == "powsybl_import." * field[length("powsybl_import_")+1:end]
      end
      # block F (0.20.0): the shipped bundles under data/powsybl are offered
      # in the case selector and staged into the case directory as a whole
      root = dirname(@__DIR__)
      offered = SparlectraApp._webui_bundled_scf_options(root)
      @test "ieee14.powsybl" in offered && "four_substations.powsybl" in offered && "micro_grid_be.powsybl" in offered
      mktempdir() do cache
        staged = SparlectraApp._webui_stage_bundled_case!(root, cache, "ieee14.powsybl")
        @test staged == joinpath(cache, "ieee14.powsybl")
        @test isfile(joinpath(staged, "manifest.json"))
        @test Sparlectra.detect(PowsyblAdapter, staged)
      end
      topics = [t for t in keys(SparlectraApp.WEBUI_HELP_TOPICS) if startswith(t, "powsybl_import.")]
      @test length(topics) == 5
      @test all(startswith(String(SparlectraApp.WEBUI_HELP_TOPICS[t].doc), "powsybl_import/#") for t in topics)
      @test occursin("(@id powsybl-import-config)", read(joinpath(dirname(@__DIR__), "docs", "src", "powsybl_import.md"), String))
    end)() end

    # block E.1 (0.20.0): the PST of four_substations carries its magnetizing
    # admittance on the branch, so its outage takes the admittance away and
    # leaves the shunt list alone; the property is the diagonal difference
    # of the Y-bus at both PST buses, exactly the PST's own two-port entries
    @testset "N-1 outage of the PST leaves no branch-derived shunt (four_substations)" begin (function ()
      tables = read_powsybl_bundle(_powsybl_fixture("four_substations"))
      net, _ = build_net_from_powsybl(tables, PowsyblAdapterOptions())
      k = findfirst(br -> occursin("TWT", getCompName(br.comp)) || (br.ratio != 0.0 && br.angle != 0.0), net.branchVec)
      @test k !== nothing
      pst = net.branchVec[k]
      @test length(net.shuntVec) == length(tables.shunt_compensators.id)
      Ybase = createYBUS(net = net, sparse = false)
      out = deepcopy(net)
      setBranchStatus!(out.branchVec[k], false)
      Yout = createYBUS(net = out, sparse = false)
      @test length(out.shuntVec) == length(net.shuntVec)
      y11, _, _, y22 = calcAdmittance(pst, pst.comp.cVN, net.baseMVA)
      f = Int(pst.fromBus)
      t = Int(pst.toBus)
      @test isapprox(Ybase[f, f] - Yout[f, f], y11; atol = 1e-12)
      @test isapprox(Ybase[t, t] - Yout[t, t], y22; atol = 1e-12)
      @test Yout[f, t] == 0.0
      # the engine solves the outage
      cases = [ContingencyCase(getCompName(pst.comp), :branch, getCompName(pst.comp))]
      res = runContingencies!(net, cases; parallel_enabled = false)
      @test length(res) == 1
      @test length(net.shuntVec) == length(tables.shunt_compensators.id)
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
      end
    end)() end

    # block E.3 (0.20.0): with the magnetizing admittances on the branches a
    # measurement set drawn from the solved power flow (Gaussian noise at
    # the stated sigmas, so the objective sits inside the chi-square band)
    # passes the estimator's global test with no suspicious row
    @testset "SE self-test on ieee14 and micro_grid_be" begin (function ()
      std = measurementStdDevs(vm = 1e-4, pinj = 1e-3, qinj = 1e-3, pflow = 1e-3, qflow = 1e-3)
      for (case, opts) in (("ieee14", PowsyblAdapterOptions()), ("micro_grid_be", PowsyblAdapterOptions(remote_regulation = :remote)))
        tables = read_powsybl_bundle(_powsybl_fixture(case))
        net, _ = build_net_from_powsybl(tables, opts)
        res = run_sparlectra(net = net, config = _powsybl_olf_config())
        @test res.numerical_converged
        meas = generateMeasurementsFromPF(net; includeVm = true, includePinj = true, includeQinj = true, includePflow = true, includeQflow = true, noise = true, stddev = std, rng = MersenneTwister(7))
        diag = runse_diagnostics(net, meas)
        summary = summarize_se_diagnostics(diag)
        @test summary.global_consistency
        @test summary.suspicious_count == 0
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
      bands = Dict("ieee14" => (1e-6, 1e-4), "ieee57" => (1e-6, 1e-4), "four_substations" => (1e-5, 1e-3), "micro_grid_be" => (1e-5, 1e-3), "eurostag_tie_lines" => (1e-5, 1e-3))
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
        # the power flow from the .xiidm lands in the bands of the bundle
        opts = case == "micro_grid_be" ? PowsyblAdapterOptions(remote_regulation = :remote) : PowsyblAdapterOptions()
        net, report = build_net_from_powsybl(native, opts)
        res = run_sparlectra(net = net, config = _powsybl_olf_config())
        @test res.numerical_converged
        dev = _powsybl_voltage_deviation(net, bundle, dir)
        dv_band, da_band = bands[case]
        @test dev.max_dv <= dv_band
        @test dev.max_da <= da_band
        println("      powsybl $(case) from .xiidm: max |dV| $(round(dev.max_dv; sigdigits = 3)) pu, max |dtheta| $(round(dev.max_da; sigdigits = 3)) deg")
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
      err = try
        Sparlectra.read_iidm_tables(joinpath(_powsybl_fixture("ieee14"), "ieee14.xiidm.bz2"))
        nothing
      catch e
        e
      end
      @test err isa ArgumentError && occursin("unpack", err.msg)
      mktempdir() do dir
        bad = joinpath(dir, "notiidm.xiidm")
        write(bad, "<?xml version=\"1.0\"?><root/>")
        err = try
          Sparlectra.read_iidm_tables(bad)
          nothing
        catch e
          e
        end
        @test err isa ArgumentError && occursin("not <network>", err.msg)
      end
    end)() end

    @testset "import_case dispatch on the bundle and the .xiidm file (ieee14)" begin (function ()
      cfg = load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH; reload = true)
      bundle = _powsybl_fixture("ieee14")
      ic = Sparlectra.import_case(bundle, cfg)
      @test ic.format == :powsybl
      @test ic.provenance["powsybl_report"] isa PowsyblImportReport
      @test ic.provenance["powsybl_manifest"]["case"] == "ieee14"
      ite, erg = runpf!(ic.net, 40, 1e-8, 0; method = :rectangular)
      @test erg == 0
      # the .xiidm file imports natively, whether or not the extension is loaded
      xiidm = joinpath(bundle, "ieee14.xiidm")
      @test Sparlectra.detect(PowsyblAdapter, xiidm)
      ix = Sparlectra.import_case(xiidm, cfg)
      @test ix.format == :powsybl
      @test length(ix.net.nodeVec) == 14
      @test ix.provenance["powsybl_manifest"]["source"] == "Sparlectra IIDM reader"
      # the service entry accepts the bundle directory and the file
      run_cfg = SparlectraConfig(powerflow = PowerFlowConfig(max_iter = 40, tol = 1e-8), output = OutputConfig(logfile_results = :off, startup_latency_hint = false))
      res = run_sparlectra(casefile = bundle, config = run_cfg)
      @test res.numerical_converged
      res_file = run_sparlectra(casefile = xiidm, config = run_cfg)
      @test res_file.numerical_converged
      @test length(res_file.net.nodeVec) == 14
      # convert_case captures the network with the report in the meta
      case = convert_case(PowsyblAdapter(), bundle)
      @test case.sparlectra !== nothing
      @test haskey(case.sparlectra.meta, "powsybl_import")
    end)() end
  end)() end
end
