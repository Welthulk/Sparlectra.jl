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

# file: test/fast/ref_entsoe_sv.jl
# purpose: profile `fast`: the ENTSO-E SmallGrid
#          conformity set (test/fixtures/cgmes_smallgrid) against its own SV
#          profile, and MicroGrid BE and NL and ReliCapGrid Galia from the
#          CGMES cache when present;
#          the SmallGrid state estimation with a noise-free table.
using Sparlectra
using Test

const _REF_SMALLGRID_DIR = normpath(joinpath(@__DIR__, "..", "fixtures", "cgmes_smallgrid"))
_ref_smallgrid_import() = redirect_stdout(() -> importCGMES(path = sort!(filter(f -> endswith(f, ".zip"), readdir(_REF_SMALLGRID_DIR; join = true))), name = "smallgrid"), devnull)

# Solves a CGMES import and compares it with the SV profile of its own
# delivery: voltages, angles, branch flows, and the generation of every
# machine bus (SvPowerFlow of the unit terminals, load sign) against the bus
# generation of the result report. Units rows are one bus each, named by
# comp.cName, in these sets (no closed links on a machine bus).
function _ref_check_against_sv(res, bands)
  _, erg = runpf!(res.net, 40, 1e-8, 0)
  @test erg == 0
  cmp = compareWithSV(res)
  @test cmp.n == bands.n
  @test cmp.max_dvm < bands.dvm
  @test cmp.max_dva < bands.dva
  @test cmp.flows.max_dp < bands.dp
  @test cmp.flows.max_dq < bands.dq
  units = [r for r in cmp.flows.rows if r.kind == :units]
  @test length(units) == bands.units
  rows = Dict(x.bus => x for x in buildACPFlowReport(res.net).nodes)
  by_name = Dict(n.comp.cName => k for (k, n) in enumerate(res.net.nodeVec))
  for r in units
    row = rows[by_name[r.bus]]
    @test isapprox(row.p_gen_MW, -r.sv_p; atol = bands.dp)
    @test isapprox(row.q_gen_MVar, -r.sv_q; atol = bands.dq)
  end
  return nothing
end

function run_ref_entsoe_sv_tests()
  @testset "ENTSO-E SmallGrid power flow against the SV profile" begin (function ()
    # bands: the measured state of the SmallGrid import against its SV
    # (max dvm 0.00014 pu, 0.01 deg, flows max 0.1 MW) with headroom
    _ref_check_against_sv(_ref_smallgrid_import(), (n = 118, dvm = 5e-4, dva = 0.05, dp = 0.5, dq = 0.5, units = 22))
  end)() end

  # MicroGrid (BE and NL) is not part of the repository: it is read from the
  # ENTSO-E package in the CGMES cache (data/CGMES or SPARLECTRA_CGMES_CACHE)
  micro = [
    (alias = "microgrid_be", bands = (n = 11, dvm = 5e-4, dva = 0.05, dp = 0.1, dq = 0.7, units = 7)),
    (alias = "microgrid_nl", bands = (n = 10, dvm = 5e-4, dva = 0.05, dp = 2.5, dq = 2.5, units = 7)),
  ]
  for m in micro
    zip = joinpath(Sparlectra.CGMESImporter.cgmesTestSetCacheDir(), "cases", "cgmes_" * m.alias * ".zip")
    if !isfile(zip)
      println("      entsoe_sv: ", m.alias, " against the SV profile SKIPPED (", zip, " not in the CGMES cache)")
      continue
    end
    println("      entsoe_sv: ", m.alias, " against the SV profile RAN")
    @testset "ENTSO-E $(m.alias) power flow against the SV profile" begin (function ()
      _ref_check_against_sv(redirect_stdout(() -> importCGMES(path = [zip], name = m.alias), devnull), m.bands)
    end)() end
  end

  # ReliCapGrid Galia (CGMES 3.0) from the CGMES cache, alone with its border
  # to Belgovia open. Its phase shifter GA-TR2_1 (PhaseTapChangerLinear,
  # xMin 27.9674, xMax 29.7154 ohm) sits at step 3, 6 of at most 40 degrees,
  # and the SV state carries the reactance of the xMin/xMax curve:
  # max dvm 6.9e-6 pu with the curve, 2.6e-5 pu with the end reactance.
  galia_files = Sparlectra.CGMESImporter._rcgMemberFiles("galia")
  if !all(f -> isfile(joinpath(f[3], f[1])), galia_files)
    println("      entsoe_sv: relicap galia against the SV profile SKIPPED (member files not in the CGMES cache)")
  else
    println("      entsoe_sv: relicap galia against the SV profile RAN")
    @testset "ReliCapGrid Galia power flow against the SV profile" begin (function ()
      zip = Sparlectra.CGMESImporter.fetchReliCapGridSet("galia"; outdir = mktempdir())
      res = redirect_stdout(() -> importCGMES(path = [zip], name = "galia", require_boundary = false), devnull)
      @test any(m -> occursin("GA-TR2_1", m) && occursin("xMin/xMax curve", m), res.messages)
      _ref_check_against_sv(res, (n = 12, dvm = 1.5e-5, dva = 0.01, dp = 0.05, dq = 0.05, units = 9))
    end)() end
  end

  @testset "SmallGrid: a noise-free measurement table reproduces the SV state" begin (function ()
    result = _ref_smallgrid_import()
    rows = generateMeasurementTableFromSV(result; noise = false)
    set = importMeasurementTable([MeasurementRow(i, r.timestamp, r.measurement_mrid, r.terminal_mrid, r.node_mrid, r.type, r.value, r.sigma, r.source, r.valid, :none, "") for (i, r) in enumerate(rows)], cgmesMeasurementKeys(result)).sets[1]
    @test length(set.measurements) == length(rows)
    # without zero-injection rows: the SV state leaves up to 0.14 MVAr at
    # the passive buses, ZI rows at 1 kW would pull the estimate off it
    applyMeasurementSet!(result.net, set; zero_injection = false)
    se = with_state_estimation_config(update_net = false) do
      runse!(result.net)
    end
    @test se.converged
    @test se.objectiveJ < 1e-8
    dvm = [abs(abs(se.voltages[m.busIdx]) - m.value) for m in set.measurements if m.typ == Sparlectra.VmMeas]
    @test maximum(dvm) < 1e-8
  end)() end
  return true
end
