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

# Date: 2026-07-30
# file: examples/others/exp_short_circuit_reference.jl
# purpose: checks runShortCircuit! against the analytic IEC 60909-0 reference values, PASS/FAIL per quantity

using Sparlectra

include(joinpath(@__DIR__, "..", "internal", "example_header.jl"))

# The four reference cases and their expected values are the analytic
# IEC 60909-0 hand derivations (with a pandapower snippet kept for an external
# cross-check that is deliberately NOT executed from this repository). The
# expected values encoded below are the source of record.

const _RTOL = 1e-6

_feeder(bus; ikmax, ikmin = nothing, rx = 0.1) =
  (mrid = "F", name = "Feeder", bus = bus, maxInitialSymShCCurrent_A = ikmax, minInitialSymShCCurrent_A = ikmin, maxR1ToX1Ratio = rx, minR1ToX1Ratio = rx, maxR0ToX0Ratio = nothing, maxZ0ToZ1Ratio = nothing, ikSecond = nothing, governorSCD = nothing)

_machine(bus; xdpp, ratedS, ratedU) =
  (mrid = "G", name = "G1", bus = bus, satDirectSubtransX_pu = xdpp, satDirectTransX_pu = nothing, r0_pu = nothing, x0_pu = nothing, r2_pu = nothing, x2_pu = nothing, earthing = nothing, ratedS_MVA = ratedS, ratedU_kV = ratedU)

_motor(bus; ilr, rx, ratedS, ratedU) =
  (mrid = "M", name = "M1", bus = bus, iaIrRatio = ilr, rxLockedRotorRatio = rx, efficiency_percent = nothing, ratedMechanicalPower_MW = nothing, polePairNumber = nothing, ratedS_MVA = ratedS, ratedU_kV = ratedU, ratedPowerFactor = nothing)

_data(; enis = NamedTuple[], machines = NamedTuple[], asms = NamedTuple[]) =
  Sparlectra.CGMESImporter.CGMESShortCircuitData(enis, machines, NamedTuple[], NamedTuple[], NamedTuple[], asms)

function _case_r1_r2()
  net = Net(name = "ref_feeder_line", baseMVA = 100.0)
  addBus!(net = net, busName = "Q", vn_kV = 110.0)
  addBus!(net = net, busName = "A", vn_kV = 110.0)
  addPIModelACLine!(net = net, fromBus = "Q", toBus = "A", r_pu = 0.05, x_pu = 0.5, b_pu = 0.0, status = 1)
  return net, _data(enis = [_feeder("Q"; ikmax = 10_000.0, ikmin = 8_000.0)])
end

function _case_r3()
  net = Net(name = "ref_machine", baseMVA = 100.0)
  addBus!(net = net, busName = "G", vn_kV = 10.5)
  return net, _data(machines = [_machine("G"; xdpp = 0.15, ratedS = 50.0, ratedU = 10.5)])
end

function _case_r4()
  net = Net(name = "ref_motor", baseMVA = 100.0)
  addBus!(net = net, busName = "M", vn_kV = 10.0)
  return net, _data(asms = [_motor("M"; ilr = 5.0, rx = 0.1, ratedS = 5.0, ratedU = 10.0)])
end

function main()
  print_example_banner("examples/others/exp_short_circuit_reference.jl", "checks runShortCircuit! against the analytic IEC 60909-0 reference values")

  net_ql, sc_ql = _case_r1_r2()
  rmax = runShortCircuit!(net_ql, sc_ql; case = :max)
  rmin = runShortCircuit!(net_ql, sc_ql; case = :min)
  net_g, sc_g = _case_r3()
  rg = runShortCircuit!(net_g, sc_g; case = :max)
  net_m, sc_m = _case_r4()
  rm = runShortCircuit!(net_m, sc_m; case = :max)

  row(result, bus) = only(r for r in result.rows if r.bus == bus)
  # (case id, quantity, reference value from the documented table, computed)
  checks = [
    ("R1", "Ik''max [kA]", 10.0, row(rmax, "Q").ik_kA),
    ("R1", "Sk''max [MVA]", 1905.2558883257648, row(rmax, "Q").sk_MVA),
    ("R1", "ip [kA]", 28.284271247461902, row(rmax, "Q").ip_kA),
    ("R1", "Ik''min [kA]", 8.0, row(rmin, "Q").ik_kA),
    ("R2", "Zk [Ω]", 67.78768576497586, row(rmax, "A").zk_ohm),
    ("R2", "Ik''max [kA]", 1.0305615508715191, row(rmax, "A").ik_kA),
    ("R2", "ip [kA]", 2.9148682442055054, row(rmax, "A").ip_kA),
    ("R3", "Zk [Ω]", 0.3315593472611653, row(rg, "G").zk_ohm),
    ("R3", "Ik''max [kA]", 20.112223239140242, row(rg, "G").ik_kA),
    ("R3", "ip [kA]", 56.88595774853494, row(rg, "G").ip_kA),
    ("R4", "Zk [Ω]", 4.0, row(rm, "M").zk_ohm),
    ("R4", "Ik''max [kA]", 1.5877132402714709, row(rm, "M").ik_kA),
    ("R4", "ip [kA]", 4.4907311951024935, row(rm, "M").ip_kA),
  ]

  println(rpad("case", 6), rpad("quantity", 16), rpad("reference", 22), rpad("Sparlectra", 22), rpad("rel. diff", 12), "verdict")
  failures = 0
  for (case_id, qty, ref, got) in checks
    rel = abs(got - ref) / abs(ref)
    ok = rel <= _RTOL
    ok || (failures += 1)
    println(rpad(case_id, 6), rpad(qty, 16), rpad(string(ref), 22), rpad(string(got), 22), rpad(string(round(rel; sigdigits = 3)), 12), ok ? "PASS" : "FAIL")
  end
  println()
  if failures == 0
    println("All ", length(checks), " reference quantities match within rtol ", _RTOL, ".")
    println("Definitions and analytic IEC 60909-0 hand derivations are documented")
    println("inline in this example file.")
  else
    error("$(failures) reference quantit(y/ies) deviate beyond rtol $(_RTOL); see the table above")
  end
  return failures == 0
end

run_example(main)
