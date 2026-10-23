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

# file: test/fast/ref_openloadflow.jl
# purpose: profile `fast`: the five shipped PowSyBl
#          bundles (test/fixtures/powsybl) against the OpenLoadFlow solution
#          they carry (bus voltages, line flows, generator outputs). Voltage
#          and flow legs are copies from test_powsybl_importer.jl, the
#          generator leg is new. Since 0.50.0 also the test grids, the CIGRE
#          HV net and sp_case14 against OpenLoadFlow solutions of their
#          MATPOWER export (test/fixtures/openloadflow), since 0.60.0 also
#          MATPOWER case14 and case118.
using Sparlectra
using Test
using DelimitedFiles
using SHA

const _REF_POWSYBL = normpath(joinpath(@__DIR__, "..", "fixtures", "powsybl"))
const _REF_POWSYBL_CASES = ["ieee14", "ieee57", "four_substations", "micro_grid_be", "eurostag_tie_lines"]
_ref_powsybl_fixture(case::AbstractString) = joinpath(_REF_POWSYBL, case * ".powsybl")

# reference_buses.csv of a bundle: bus-breaker id => (v_pu, v_angle_deg, component, bus-view id)
function _ref_powsybl_reference(dir::AbstractString)
  data, header = readdlm(joinpath(dir, "reference_buses.csv"), ',', String; quotes = true, header = true)
  cols = Dict(Symbol(strip(h)) => k for (k, h) in enumerate(vec(header)))
  rows = Dict{String,NamedTuple{(:v_pu, :va, :comp, :bus_view),Tuple{Float64,Float64,Int,String}}}()
  for i in axes(data, 1)
    rows[String(data[i, cols[:bus_breaker_id]])] = (v_pu = parse(Float64, data[i, cols[:v_pu]]), va = parse(Float64, data[i, cols[:v_angle_deg]]), comp = parse(Int, data[i, cols[:synchronous_component]]), bus_view = String(data[i, cols[:id]]))
  end
  return rows
end

# the OLF reference bus per synchronous component, from the manifest
function _ref_powsybl_reference_buses(tables)
  refs = Dict{Int,String}()
  for c in get(get(tables.manifest, "reference", Dict()), "components", [])
    refs[Int(c["synchronous_component_num"])] = String(c["reference_bus_id"])
  end
  return refs
end

# largest voltage and angle deviation against the OLF reference; angles
# relative to the OLF reference bus of each component
function _ref_powsybl_voltage_deviation(net, tables, dir)
  ref = _ref_powsybl_reference(dir)
  refbuses = _ref_powsybl_reference_buses(tables)
  offsets = Dict{Int,Float64}()
  for (bus, r) in ref
    haskey(net.busDict, bus) || continue
    r.bus_view == get(refbuses, r.comp, "") && (offsets[r.comp] = net.nodeVec[net.busDict[bus]]._va_deg - r.va)
  end
  max_dv = 0.0
  max_da = 0.0
  for (bus, r) in ref
    haskey(net.busDict, bus) || continue
    node = net.nodeVec[net.busDict[bus]]
    max_dv = max(max_dv, abs(node._vm_pu - r.v_pu))
    max_da = max(max_da, abs(node._va_deg - r.va - get(offsets, r.comp, 0.0)))
  end
  return (max_dv = max_dv, max_da = max_da)
end

# the run configuration that reproduces OpenLoadFlow's defaults (copy of
# _powsybl_olf_config): every synchronous component solved, the P imbalance
# shared by the imported participation factors, no Q-limit hysteresis
function _ref_olf_config()
  pf = PowerFlowConfig(max_iter = 60, tol = 1e-9, islands_enabled = true, qlimits = Sparlectra.QLimitConfig(hysteresis_pu = 1e-6), distributed_slack = Sparlectra.DistributedSlackConfig(enabled = true, p_mode = :imported))
  return SparlectraConfig(powerflow = pf, output = OutputConfig(logfile_results = :off, startup_latency_hint = false), control = ControlConfig())
end

_ref_olf_options(case) = case == "micro_grid_be" ? PowsyblAdapterOptions(remote_regulation = :remote) : PowsyblAdapterOptions()

# Sparlectra's own test nets against OpenLoadFlow (0.50.0): the shared test
# grids (createTest3BusNet, also with its reactive limit of 15 Mvar binding,
# and createTest5BusNet), the CIGRE HV benchmark (createCIGRE) and the
# shipped sp_case14; since 0.60.0 MATPOWER case14 and case118, whose file
# carries only the rounded IEEE CDF state and so no usable reference. The
# two are read from their fixture .m (the MATPOWER export of the cached
# case), so the test downloads nothing. The limit: the fixture is a file
# Sparlectra wrote (case14 with baseKV 1.0 instead of 0), so the test checks
# the solver against OpenLoadFlow on that model, not the import of the
# original MATPOWER file.
# test/fixtures/openloadflow holds per case the MATPOWER
# model file the net was exported to and the OLF solution of that file; the
# header of each CSV names the tool, its parameters and the sha256 of the
# input. The nets are built here as the other tests build them, so a change
# of a builder or of sp_case14 shows up as a deviation.
const _REF_OLF_NETS = normpath(joinpath(@__DIR__, "..", "fixtures", "openloadflow"))
const _REF_OLF_NET_CASES = [
  (case = "test3bus", build = () -> createTest3BusNet()),
  (case = "test3bus_qlim15", build = () -> createTest3BusNet(qlim_min = -15.0, qlim_max = 15.0)),
  (case = "test5bus", build = () -> createTest5BusNet(qlim_min = -200.0, qlim_max = 200.0)),
  (case = "cigre", build = () -> createCIGRE()),
  (case = "sp_case14", build = () -> _ref_olf_flat!(Sparlectra.import_case(normpath(joinpath(@__DIR__, "..", "..", "data", "scf", "sp_case14.scf.json")), Sparlectra.SparlectraConfig(Dict())).net)),
  (case = "case14", build = () -> _ref_olf_flat!(Sparlectra.import_case(joinpath(_REF_OLF_NETS, "case14.m"), Sparlectra.SparlectraConfig(Dict())).net)),
  (case = "case118", build = () -> _ref_olf_flat!(Sparlectra.import_case(joinpath(_REF_OLF_NETS, "case118.m"), Sparlectra.SparlectraConfig(Dict())).net)),
]

# sp_case14 carries its solved state as start state, so Newton would stop
# after one step; a flat start (setpoints kept on PV and slack buses) makes
# the solver do the work
function _ref_olf_flat!(net)
  for n in net.nodeVec
    n._nodeType == Sparlectra.PQ && (n._vm_pu = 1.0)
    n._va_deg = 0.0
  end
  return net
end

# the run configuration of these nets: one slack, no distributed slack, the
# reactive limits enforced (they bind only in test3bus_qlim15 and case118,
# where OLF runs with its reactive limits on as well; both switch the same
# six buses 19, 32, 34, 92, 103 and 105 of case118 to PQ; OLF with limits
# would also hold the Q limit of the slack, which in case14 binds and which
# Sparlectra and MATPOWER do not apply), controllers off (sp_case14
# declares a tap controller; the fixture is the file state, as the demo-case
# test solves it with runpf!)
function _ref_olf_net_config()
  pf = PowerFlowConfig(max_iter = 50, tol = 1e-10)
  return SparlectraConfig(powerflow = pf, output = OutputConfig(logfile_results = :off, startup_latency_hint = false), control = ControlConfig(enabled = false))
end

# rows of a fixture CSV as Dicts keyed by column name, header lines skipped
function _ref_olf_rows(path::AbstractString)
  data, header = readdlm(path, ',', String; header = true, comments = true, comment_char = '#')
  cols = strip.(vec(header))
  return [Dict(String(c) => String(data[i, k]) for (k, c) in enumerate(cols)) for i in axes(data, 1)]
end

# the input hash the fixture header states, against the shipped input file
# (line endings normalized, so a CRLF checkout compares equal)
function _ref_olf_input_matches(case::AbstractString)
  m = match(r"# input: (\S+) sha256 ([0-9a-f]{64})", read(joinpath(_REF_OLF_NETS, case * "_buses.csv"), String))
  m === nothing && return false
  text = replace(read(joinpath(_REF_OLF_NETS, m.captures[1]), String), "\r\n" => "\n")
  return bytes2hex(SHA.sha256(text)) == m.captures[2]
end

function run_ref_openloadflow_tests()
  @testset "PowSyBl bundles: power flow against OpenLoadFlow ($(case))" for case in _REF_POWSYBL_CASES
    (function ()
      bands = Dict("ieee14" => (1e-6, 1e-4), "ieee57" => (1e-6, 1e-4), "four_substations" => (1e-5, 1e-3), "micro_grid_be" => (1e-5, 1e-3), "eurostag_tie_lines" => (1e-5, 1e-3))
      dir = _ref_powsybl_fixture(case)
      tables = read_powsybl_bundle(dir)
      net, _ = build_net_from_powsybl(tables, _ref_olf_options(case))
      res = run_sparlectra(net = net, config = _ref_olf_config())
      @test res.numerical_converged
      dev = _ref_powsybl_voltage_deviation(net, tables, dir)
      dv_band, da_band = bands[case]
      @test dev.max_dv <= dv_band
      @test dev.max_da <= da_band
      # generator outputs: OLF writes p and q of every generator in load
      # sign; summed per bus they are the bus generation of the result
      # report. Not on four_substations: its generators sit on node-breaker
      # buses joined by closed switches next to VSC and LCC stations, and
      # the split of one electrical bus among its bus-breaker buses is not
      # comparable per bus (measured deviations up to 17 Mvar there with the
      # voltages inside 1e-5 pu)
      case == "four_substations" && return nothing
      g = tables.generators
      gen_p = Dict{String,Float64}()
      gen_q = Dict{String,Float64}()
      for i in eachindex(g.id)
        b = g.bus_breaker_bus_id[i]
        (b isa AbstractString && haskey(net.busDict, b) && isfinite(g.p[i]) && isfinite(g.q[i])) || continue
        gen_p[b] = get(gen_p, b, 0.0) - g.p[i]
        gen_q[b] = get(gen_q, b, 0.0) - g.q[i]
      end
      @test !isempty(gen_p)
      rows = Dict(x.bus => x for x in buildACPFlowReport(net).nodes)
      for (b, p) in gen_p
        row = rows[net.busDict[b]]
        @test isapprox(row.p_gen_MW, p; atol = 1e-2)
        @test isapprox(row.q_gen_MVar, gen_q[b]; atol = 1e-2)
      end
    end)()
  end

  @testset "PowSyBl ieee14: line flows against the OLF columns" begin (function ()
    tables = read_powsybl_bundle(_ref_powsybl_fixture("ieee14"))
    net, _ = build_net_from_powsybl(tables, PowsyblAdapterOptions())
    res = run_sparlectra(net = net, config = _ref_olf_config())
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
      @test abs(flow.pFlow - ln.p1[i]) <= 1e-3
      @test abs(flow.qFlow - ln.q1[i]) <= 1e-3
      checked += 1
    end
    @test checked == length(ln.id)
  end)() end

  # Y-bus of the PowSyBl ieee14 against MATPOWER case14, read from the
  # shipped OLF input case14.m (branch r, x, b, ratio and bus Gs, Bs are the
  # MATPOWER values; the file adds names and flat voltages). Branches 7-8
  # (14/20 kV) and 7-9 (14/12 kV) carry tap 1 in MATPOWER and are lines
  # between voltage levels in PowSyBl's IEEE-CDF import, with large opposite
  # b1/b2 that turn the physical conductor into the nominal-ratio
  # transformer; the per-terminal arms of the branch reproduce them without
  # any bus shunt, so the identity holds everywhere.
  @testset "PowSyBl ieee14: Y-bus against MATPOWER case14" begin (function ()
    tables = read_powsybl_bundle(_ref_powsybl_fixture("ieee14"))
    net, _ = build_net_from_powsybl(tables, PowsyblAdapterOptions())
    mp = createNetFromMatPowerFile(filename = joinpath(_REF_OLF_NETS, "case14.m"))
    # bus-breaker ids B1..B14 in bus order, MATPOWER buses 1..14 in bus order
    order = sortperm([parse(Int, match(r"(\d+)$", b).match) for b in tables.bus_breaker_view_buses.id])
    Y = Matrix(createYBUS(net = net))[order, order]
    Ym = Matrix(createYBUS(net = mp))
    @test maximum(abs.(Y .- Ym)) <= 1e-9
  end)() end

  # Bands: the models are identical (the export script checks the Y-bus of
  # the written rows against createYBUS to 1e-12), so the bands are the
  # existing exact-model bands of ieee14: 1e-6 pu, 1e-4 deg, 1e-3 MW/Mvar.
  @testset "Sparlectra test nets: power flow against OpenLoadFlow ($(c.case))" for c in _REF_OLF_NET_CASES
    (function ()
      @test _ref_olf_input_matches(c.case)
      net = c.build()
      res = run_sparlectra(net = net, config = _ref_olf_net_config())
      @test res.numerical_converged
      node_of = Dict(n.busIdx => n for n in net.nodeVec)
      slack = only(n for n in net.nodeVec if n._nodeType == Sparlectra.Slack)
      rows = Dict(x.bus => x for x in buildACPFlowReport(net).nodes)
      buses = _ref_olf_rows(joinpath(_REF_OLF_NETS, c.case * "_buses.csv"))
      @test length(buses) == length(net.nodeVec)
      for r in buses
        b = parse(Int, r["bus"])
        node = node_of[b]
        @test abs(node._vm_pu - parse(Float64, r["v_pu"])) <= 1e-6
        @test abs(node._va_deg - slack._va_deg - parse(Float64, r["v_angle_deg"])) <= 1e-4
        @test abs(rows[b].p_gen_MW - parse(Float64, r["p_gen_mw"])) <= 1e-3
        @test abs(rows[b].q_gen_MVar - parse(Float64, r["q_gen_mvar"])) <= 1e-3
      end
      branches = _ref_olf_rows(joinpath(_REF_OLF_NETS, c.case * "_branches.csv"))
      @test length(branches) == length(net.branchVec)
      for r in branches
        br = net.branchVec[parse(Int, r["row"])]
        @test (br.fromBus, br.toBus) == (parse(Int, r["from_bus"]), parse(Int, r["to_bus"]))
        @test abs(br.fBranchFlow.pFlow - parse(Float64, r["p_from_mw"])) <= 1e-3
        @test abs(br.fBranchFlow.qFlow - parse(Float64, r["q_from_mvar"])) <= 1e-3
        @test abs(br.tBranchFlow.pFlow - parse(Float64, r["p_to_mw"])) <= 1e-3
        @test abs(br.tBranchFlow.qFlow - parse(Float64, r["q_to_mvar"])) <= 1e-3
      end
    end)()
  end
  return true
end
