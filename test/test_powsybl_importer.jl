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
using Dates

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
    # the shipped IIDM files under data/powsybl are offered in the case
    # selector and staged into the case directory; a table bundle is a
    # test reference and not offered
    root = dirname(@__DIR__)
    offered = SparlectraApp._webui_bundled_scf_options(root)
    @test "ieee14.xiidm" in offered && "four_substations.xiidm" in offered && "micro_grid_be.xiidm" in offered && "ieee14_sc.xiidm" in offered
    @test !any(endswith(".powsybl"), offered)
    @test !SparlectraApp._webui_is_user_selectable_case("ieee14.powsybl")
    mktempdir() do cache
      staged = SparlectraApp._webui_stage_bundled_case!(root, cache, "ieee14.xiidm")
      @test staged == joinpath(cache, "ieee14.xiidm")
      @test isfile(staged)
      @test Sparlectra.detect(PowsyblAdapter, staged)
      @test SparlectraApp._webui_stage_bundled_case!(root, cache, "ieee14.powsybl") === nothing
    end
    # every shipped file is the file of its test fixture, whose frozen
    # OpenLoadFlow reference the band tests compare against
    for case in ("ieee14", "four_substations", "micro_grid_be")
      @test read(joinpath(root, "data", "powsybl", case * ".xiidm")) == read(joinpath(_powsybl_fixture(case), case * ".xiidm"))
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

  # the generators of a PowSyBl case are its short-circuit sources: the
  # reactances of the IIDM extension generatorShortCircuit (ohm at the
  # generator's voltage level) reach the engine in pu on the machine base,
  # a generator without the extension is evaluated with the engine's
  # default and flags the rows, and the columns survive a bundle
  @testset "short-circuit sources from generatorShortCircuit (ieee14)" begin (function ()
    source = read(joinpath(_powsybl_fixture("ieee14"), "ieee14.xiidm"), String)
    ids = [String(m.captures[1]) for m in eachmatch(r"<iidm:generator id=\"([^\"]+)\"", source)]
    @test length(ids) == 5
    extension(id, attrs) = "    <iidm:extension id=\"$(id)\">\n        <gsc:generatorShortCircuit xmlns:gsc=\"http://www.powsybl.org/schema/iidm/ext/generator_short_circuit/1_0\" $(attrs)/>\n    </iidm:extension>\n"
    mktempdir() do dir
      one = joinpath(dir, "one.xiidm")
      write(one, replace(source, "</iidm:network>" => extension(ids[1], "directSubtransX=\"20.0\" directTransX=\"30.0\" stepUpTransformerX=\"5.0\"") * "</iidm:network>"))
      tables = Sparlectra.read_iidm_tables(one)
      row = findfirst(==(ids[1]), tables.generators.id)
      @test tables.generators.direct_subtrans_x[row] == 20.0
      @test tables.generators.direct_trans_x[row] == 30.0
      @test tables.generators.step_up_transformer_x[row] == 5.0
      @test count(isnan, tables.generators.direct_subtrans_x) == 4
      net, _ = build_net_from_powsybl(tables, PowsyblAdapterOptions())
      machines = net.sc_sources.synchronous_machines
      @test length(machines) == 5
      machine = only(m for m in machines if m.mrid == ids[1])
      vn = Sparlectra.getNodeVn(net.nodeVec[net.busDict[machine.bus]])
      base = something(machine.ratedS_MVA, net.baseMVA)
      # the step-up reactance lies in series with the subtransient one
      @test machine.satDirectSubtransX_pu ≈ 25.0 * base / vn^2 rtol = 1e-12
      @test machine.satDirectTransX_pu ≈ 30.0 * base / vn^2 rtol = 1e-12
      @test count(m -> m.satDirectSubtransX_pu === nothing, machines) == 4
      # the four generators without the extension carry no rating either
      # (the file states none), which is a reason of its own
      flagged = run_with_expected_warnings(["has no usable x''_d", "has no usable ratedS", "has neither x''_d nor ratedS"]) do
        runShortCircuit!(net; case = :max)
      end
      @test all(r -> r.status === :ok, flagged.rows)
      @test any(r -> r.contains_defaulted_data, flagged.rows)
      @test any(r -> any(why -> occursin("has neither x''_d nor ratedS", why), r.reasons), flagged.rows)
      # the columns are part of the bundle form
      write_powsybl_bundle(tables, joinpath(dir, "one.powsybl"))
      again = read_powsybl_bundle(joinpath(dir, "one.powsybl"))
      @test again.generators.direct_subtrans_x[findfirst(==(ids[1]), again.generators.id)] == 20.0
      # every generator with data and a rating: nothing is substituted
      all_file = joinpath(dir, "all.xiidm")
      rated = replace(source, r"<iidm:generator id=" => "<iidm:generator ratedS=\"100.0\" id=")
      write(all_file, replace(rated, "</iidm:network>" => join(extension(id, "directSubtransX=\"20.0\"") for id in ids) * "</iidm:network>"))
      full, _ = build_net_from_powsybl(Sparlectra.read_iidm_tables(all_file), PowsyblAdapterOptions())
      clean = run_with_expected_warnings(String[]) do
        runShortCircuit!(full; case = :max)
      end
      @test all(r -> r.status === :ok && !r.contains_defaulted_data, clean.rows)
      # a wind generator is a converter unit: its record carries the mark
      # and the engine takes it as a current source, k times the rated
      # current from ratedS (as a CGMES machine on a WindGeneratingUnit)
      wind_file = joinpath(dir, "wind.xiidm")
      write(wind_file, replace(source, "<iidm:generator id=\"$(ids[2])\" energySource=\"OTHER\"" => "<iidm:generator id=\"$(ids[2])\" energySource=\"WIND\" ratedS=\"60.0\""))
      wnet, _ = build_net_from_powsybl(Sparlectra.read_iidm_tables(wind_file), PowsyblAdapterOptions())
      wrec = only(m for m in wnet.sc_sources.synchronous_machines if m.mrid == ids[2])
      @test wrec.converter_unit == "WindGeneratingUnit"
      @test count(m -> m.converter_unit !== nothing, wnet.sc_sources.synchronous_machines) == 1
      wrow = k -> only(run_with_expected_warnings(() -> runShortCircuit!(wnet; buses = [wrec.bus], case = :max, converter_current_factor = k), ["short-circuit: "]).rows)
      @test isapprox(wrow(1.5).ik_kA - wrow(0.0).ik_kA, 1.5 * 60.0 / (sqrt(3.0) * wrec.ratedU_kV); rtol = 1e-9)
      # the short-circuit data and the ratings change nothing in the power
      # flow: the same voltages as the file without them, bit for bit
      plain, _ = build_net_from_powsybl(Sparlectra.read_iidm_tables(joinpath(_powsybl_fixture("ieee14"), "ieee14.xiidm")), PowsyblAdapterOptions())
      for net in (plain, full)
        _, erg = runpf!(net, 40, 1e-10, 0; method = :rectangular)
        @test erg == 0
      end
      @test [n._vm_pu for n in plain.nodeVec] == [n._vm_pu for n in full.nodeVec]
      @test [n._va_deg for n in plain.nodeVec] == [n._va_deg for n in full.nodeVec]
      # the shipped demo file is this network; it has no reference
      # solution of its own, its power flow is the one of ieee14
      shipped = joinpath(dirname(@__DIR__), "data", "powsybl", "ieee14_sc.xiidm")
      demo_bundle, _ = build_net_from_powsybl(Sparlectra.read_iidm_tables(shipped), PowsyblAdapterOptions())
      _, demo_erg = runpf!(demo_bundle, 40, 1e-10, 0; method = :rectangular)
      @test demo_erg == 0
      @test maximum(abs.([n._vm_pu for n in demo_bundle.nodeVec] .- [n._vm_pu for n in plain.nodeVec])) < 1e-12
      @test count(m -> m.satDirectSubtransX_pu !== nothing && m.ratedS_MVA !== nothing, demo_bundle.sc_sources.synchronous_machines) == 5
      # the CGMES export keeps what the network carries: reactance and
      # rating of every machine come back from the delivery
      mkpath(joinpath(dir, "cgmes"))
      delivery = writeCGMESFiles(full; path = joinpath(dir, "cgmes"), created = Dates.DateTime(2026, 1, 1, 12, 0, 0), notices = String[], zip = true)
      back = Sparlectra.importCGMES(path = delivery[end], require_boundary = false).net
      returned = back.sc_sources.synchronous_machines
      @test length(returned) == 5
      @test all(m -> m.satDirectSubtransX_pu !== nothing && m.ratedS_MVA !== nothing, returned)
      @test sort([m.satDirectSubtransX_pu for m in returned]) ≈ sort([m.satDirectSubtransX_pu for m in full.sc_sources.synchronous_machines]) rtol = 1e-12
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
    # IIDM states no frequency: the conversion writes the frequency of
    # its options (from model.nominal_frequency_hz), not the global one
    cfg60 = SparlectraConfig(model = ModelConfig(nominal_frequency_hz = 60.0))
    @test Sparlectra.powsybl_adapter_options(cfg60).f_nom == 60.0
    case60 = convert_case(PowsyblAdapter(), bundle, Sparlectra.powsybl_adapter_options(cfg60))
    @test case60.sparlectra.meta["f_nom"] == 60.0
    # the c1 written at 60 Hz reads back to the same susceptance
    @test Sparlectra.createYBUS(net = Sparlectra.build_net(case60; config = cfg60), sparse = true) == Sparlectra.createYBUS(net = Sparlectra.build_net(case; config = cfg60), sparse = true)
  end)() end
end
