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
#          generator leg is new.
using Sparlectra
using Test
using DelimitedFiles

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
  return true
end
