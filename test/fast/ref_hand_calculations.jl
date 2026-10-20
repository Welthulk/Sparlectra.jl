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

# file: test/fast/ref_hand_calculations.jl
# purpose: profile `fast`: results against hand
#          calculations whose derivation stands in the comment. Copies of
#          the hand-derived checks of test_dc_powerflow.jl, testgrid.jl
#          (generation of a slack or PV bus with a load),
#          test_external_grid.jl and test_short_circuit.jl, and three closed
#          form cases: a two-bus feeder, a phase shifter, a three-winding
#          transformer.
using Sparlectra
using Test

# 3-bus DC fixture of test_dc_powerflow.jl: slack B1, generator 25 MW and
# load 20 MW at B2, load 45 MW at B3; x12 = 0.08, x13 = 0.12, x23 = 0.10.
function _ref_dc3_net()::Net
  net = Net(name = "dc3", baseMVA = 100.0)
  addBus!(net = net, busName = "B1", vn_kV = 110.0)
  addBus!(net = net, busName = "B2", vn_kV = 110.0)
  addBus!(net = net, busName = "B3", vn_kV = 110.0)
  addPIModelACLine!(net = net, fromBus = "B1", toBus = "B2", r_pu = 0.02, x_pu = 0.08, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "B1", toBus = "B3", r_pu = 0.03, x_pu = 0.12, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "B2", toBus = "B3", r_pu = 0.02, x_pu = 0.10, b_pu = 0.0, status = 1)
  addProsumer!(net = net, busName = "B1", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.02, va_deg = 0.0, referencePri = "B1")
  addProsumer!(net = net, busName = "B2", type = "GENERATOR", p = 25.0, q = 0.0, vm_pu = 1.01)
  addProsumer!(net = net, busName = "B2", type = "ENERGYCONSUMER", p = 20.0, q = 10.0)
  addProsumer!(net = net, busName = "B3", type = "ENERGYCONSUMER", p = 45.0, q = 15.0)
  refreshBusTypesFromProsumers!(net)
  return net
end

# SC record helpers of test_short_circuit.jl
const _REF_SC = Sparlectra.CGMESImporter.CGMESShortCircuitData
_ref_sc_feeder(; bus, ikmax_A = nothing, ikmin_A = nothing, rx_max = nothing, rx_min = nothing, name = "Feeder") =
  (mrid = name, name = name, bus = bus, maxInitialSymShCCurrent_A = ikmax_A, minInitialSymShCCurrent_A = ikmin_A, maxR1ToX1Ratio = rx_max, minR1ToX1Ratio = rx_min, maxR0ToX0Ratio = nothing, maxZ0ToZ1Ratio = nothing, ikSecond = nothing, governorSCD = nothing)
_ref_sc_machine(; bus, xdpp = nothing, ratedS = nothing, ratedU = nothing, name = "G1") =
  (mrid = name, name = name, bus = bus, satDirectSubtransX_pu = xdpp, satDirectTransX_pu = nothing, r0_pu = nothing, x0_pu = nothing, r2_pu = nothing, x2_pu = nothing, earthing = nothing, ratedS_MVA = ratedS, ratedU_kV = ratedU)
_ref_sc_asm(; bus, ilr = nothing, rx = nothing, ratedS = nothing, ratedU = nothing, p_mech = nothing, efficiency = nothing, cosphi = nothing, pole_pairs = nothing, name = "M1") =
  (mrid = name, name = name, bus = bus, iaIrRatio = ilr, rxLockedRotorRatio = rx, efficiency_percent = efficiency, ratedMechanicalPower_MW = p_mech, polePairNumber = pole_pairs, ratedS_MVA = ratedS, ratedU_kV = ratedU, ratedPowerFactor = cosphi)
_ref_sc_pec(; bus, maxIFault = nothing, ratedS = nothing, ratedU = nothing, name = "PEC1") =
  (mrid = name, name = name, bus = bus, maxIFault_pu = maxIFault, ratedS_MVA = ratedS, ratedU_kV = ratedU)
_ref_sc_data(; enis = NamedTuple[], machines = NamedTuple[], asms = NamedTuple[], pecs = NamedTuple[]) = _REF_SC(enis, machines, NamedTuple[], NamedTuple[], NamedTuple[], asms, pecs)

# Two-bus feeder in per unit, load P + jQ behind R + jX from a stiff source
# V1 = 1 at angle 0. With I = conj(S / V2) and V1 = V2 + Z I the magnitude
# follows from a quadratic in |V2|^2,
#   |V2|^2 = (a + sqrt(a^2 - 4 |Z|^2 |S|^2)) / 2,  a = |V1|^2 - 2 (R P + X Q),
# and the angle from tan(-d2) = (X P - R Q) / (|V2|^2 + R P + X Q). The source
# supplies S + Z |I|^2 with |I|^2 = |S|^2 / |V2|^2.
function _ref_two_bus_solution(R, X, P, Q)
  a = 1 - 2 * (R * P + X * Q)
  v2sq = (a + sqrt(a^2 - 4 * (R^2 + X^2) * (P^2 + Q^2))) / 2
  d2 = -atan((X * P - R * Q) / (v2sq + R * P + X * Q))
  i2sq = (P^2 + Q^2) / v2sq
  return (v2 = sqrt(v2sq), d2_deg = rad2deg(d2), p_src = P + R * i2sq, q_src = Q + X * i2sq)
end

_ref_node_rows(net) = Dict(r.bus => r for r in buildACPFlowReport(net).nodes)

function run_ref_hand_calculation_tests()
  @testset "DC power flow, 3-bus hand calculation" begin (function ()
    # Hand derivation (DC, lossless, theta1 = 0, 100 MVA base):
    #   P2 = 0.25 - 0.20 = 0.05 pu, P3 = -0.45 pu
    #   B22 = 1/0.08 + 1/0.10 = 22.5, B33 = 1/0.12 + 1/0.10 = 18.333, B23 = -10
    #   det = 22.5 * 18.333 - 100 = 312.5
    #   theta2 = (0.05 * 18.333 - 10 * 0.45) / 312.5 = -0.011467 rad = -0.656992 deg
    #   theta3 = (-22.5 * 0.45 + 10 * 0.05) / 312.5 = -0.030800 rad = -1.764710 deg
    #   P12 = -theta2 / 0.08 = 14.333 MW, P13 = -theta3 / 0.12 = 25.667 MW,
    #   P23 = (theta2 - theta3) / 0.10 = 19.333 MW
    #   slack = total load - other generation = 65 - 25 = 40 MW
    net = _ref_dc3_net()
    report = rundcpf!(net)
    @test report.metadata.converged
    @test report.metadata.solver === :dc
    va = Dict(row.bus_name => row.va_deg for row in report.nodes)
    @test isapprox(va["B1"], 0.0; atol = 1e-9)
    @test isapprox(va["B2"], -0.656992; atol = 1e-4)
    @test isapprox(va["B3"], -1.764710; atol = 1e-4)
    pf = Dict((row.from_bus, row.to_bus) => row.p_from_MW for row in report.branches)
    @test isapprox(pf[(1, 2)], 14.3333; atol = 1e-3)
    @test isapprox(pf[(1, 3)], 25.6667; atol = 1e-3)
    @test isapprox(pf[(2, 3)], 19.3333; atol = 1e-3)
    slack_row = only(filter(row -> row.bus_name == "B1", report.nodes))
    @test isapprox(slack_row.p_gen_MW, 40.0; atol = 1e-6)
    @test all(isapprox(row.p_to_MW, -row.p_from_MW; atol = 1e-12) for row in report.branches)
  end)() end

  @testset "generation of a slack and a PV bus that carry a load" begin (function ()
    # Hand derivation (from testgrid.jl): A slack with 60 MW load and a
    # 5 MW unit, B PV with a 50 MW machine and a 10 MW load, lossless line
    # x = 0.1 pu between two 1 pu buses. B exports 40 MW to A, so
    # sin(d) = 0.4 * 0.1 and each end draws q_line = (1 - cos(d)) / x.
    # Slack generation 60 - 40 = 20 MW, slack Q 5 + q_line, PV Q 4 + q_line.
    # Of the slack bus, the scheduled 5 MW unit keeps its 5 MW and the
    # reference unit takes the 15 MW rest.
    build() = begin
      net = Net(name = "gen_with_load", baseMVA = 100.0)
      addBus!(net = net, busName = "A", vn_kV = 110.0)
      addBus!(net = net, busName = "B", vn_kV = 110.0)
      addPIModelACLine!(net = net, fromBus = "A", toBus = "B", r_pu = 0.0, x_pu = 0.1, b_pu = 0.0, status = 1)
      addProsumer!(net = net, busName = "A", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "A")
      addProsumer!(net = net, busName = "A", type = "GENERATOR", p = 5.0, q = 1.0)
      addProsumer!(net = net, busName = "A", type = "ENERGYCONSUMER", p = 60.0, q = 5.0)
      addProsumer!(net = net, busName = "B", type = "SYNCHRONOUSMACHINE", p = 50.0, vm_pu = 1.0)
      addProsumer!(net = net, busName = "B", type = "ENERGYCONSUMER", p = 10.0, q = 4.0)
      net
    end
    q_line = (1.0 - sqrt(1.0 - 0.04^2)) / 0.1 * 100.0
    net = build()
    @test last(runpf!(net, 30, 1e-12, 0; method = :rectangular)) == 0
    a, b = net.nodeVec[1], net.nodeVec[2]
    @test isapprox(a._pƩGen, 20.0; atol = 1e-8)
    @test isapprox(a._qƩGen, 5.0 + q_line; atol = 1e-8)
    @test isapprox(b._qƩGen, 4.0 + q_line; atol = 1e-8)
    rows = buildACPFlowReport(net).nodes
    ra = only(r for r in rows if r.bus == 1)
    rb = only(r for r in rows if r.bus == 2)
    @test isapprox(ra.p_gen_MW, 20.0; atol = 1e-8)
    @test isapprox(rb.q_gen_MVar, 4.0 + q_line; atol = 1e-8)
    distributeBusResults!(net)
    @test isapprox(net.prosumpsVec[1].pRes, 15.0; atol = 1e-8)
    @test isapprox(net.prosumpsVec[2].pRes, 5.0; atol = 1e-12)
    @test isapprox(net.prosumpsVec[4].pRes, 50.0; atol = 1e-12)
    # the external-solver write-back on the same solved voltages
    model = buildPfModel(net; include_limits = false)
    V = ComplexF64[n._vm_pu * cis(deg2rad(n._va_deg)) for n in net.nodeVec]
    a._pƩGen = b._qƩGen = 0.0
    applyPfSolution!(net, model, PFSolution(V = V, converged = true, iters = 0, residual_inf = 0.0))
    @test isapprox(a._pƩGen, 20.0; atol = 1e-8)
    @test isapprox(b._qƩGen, 4.0 + q_line; atol = 1e-8)
    # DC: lossless, the same 20 MW
    dnet = build()
    rundcpf!(dnet)
    @test isapprox(dnet.nodeVec[1]._pƩGen, 20.0; atol = 1e-8)
  end)() end

  @testset "SC hand calculation, external grid on a single 110 kV bus" begin (function ()
    # Hand derivation (test_external_grid.jl): Zq = c Un^2 / Sk,
    # Ik'' = c Un / (sqrt(3) |Zq|) = Sk / (sqrt(3) Un), Sk'' = sqrt(3) Un Ik'',
    # kappa = min(1.15 (1.02 + 0.98 exp(-3 R/X)), 2.0), ip = kappa sqrt(2) Ik''.
    for rx in (0.1, 0.3)
      net = Net(name = "sc_hand", baseMVA = 100.0)
      addBus!(net = net, busName = "B1", vn_kV = 110.0)
      addExternalGrid!(net = net, busName = "B1", sk_max_MVA = 3000.0, rx_max = rx)
      res = runShortCircuit!(net; case = :max, c_factor = 1.1)
      row = only(res.rows)
      ik_hand = 3000.0 / (sqrt(3.0) * 110.0)
      kappa_hand = min(1.15 * (1.02 + 0.98 * exp(-3.0 * rx)), 2.0)
      @test row.status === :ok
      @test row.zk_ohm ≈ 1.1 * 110.0^2 / 3000.0 rtol = 1e-6
      @test row.rx_ratio ≈ rx rtol = 1e-6
      @test row.ik_kA ≈ ik_hand rtol = 1e-6
      @test row.sk_MVA ≈ 3000.0 rtol = 1e-6
      @test row.kappa ≈ kappa_hand rtol = 1e-6
      @test row.ip_kA ≈ kappa_hand * sqrt(2.0) * ik_hand rtol = 1e-6
    end
  end)() end

  @testset "SC hand calculations, IEC 60909-0 (test_short_circuit.jl)" begin (function ()
    # feeder: Z_Q = c Un / (sqrt(3) Ik_declared); at the feeder bus Zk = Z_Q,
    # so Ik'' = Ik_declared; Sk'' = sqrt(3) 110 10 MVA; R/X 0.1 gives
    # kappa_b = 1.746, 1.15 kappa_b > 2 so the HV cap 2.0 applies.
    net = Net(name = "sc_single", baseMVA = 100.0)
    addBus!(net = net, busName = "B1", vn_kV = 110.0)
    sc = _ref_sc_data(enis = [_ref_sc_feeder(bus = "B1", ikmax_A = 10_000.0, ikmin_A = 8_000.0, rx_max = 0.1, rx_min = 0.1)])
    rmax = runShortCircuit!(net, sc; case = :max)
    row = only(rmax.rows)
    @test row.status === :ok
    @test isapprox(row.ik_kA, 10.0; rtol = 1e-9)
    @test isapprox(row.sk_MVA, sqrt(3.0) * 110.0 * 10.0; rtol = 1e-9)
    @test isapprox(row.kappa, 2.0; rtol = 1e-12)
    @test isapprox(row.ip_kA, 2.0 * sqrt(2.0) * 10.0; rtol = 1e-9)
    rmin = runShortCircuit!(net, sc; case = :min)
    @test isapprox(rmin.rows[1].ik_kA, 8.0; rtol = 1e-9)
    @test rmin.rows[1].c == 1.00 && rmax.rows[1].c == 1.10

    # line: Z_base 121 ohm, Zk(B2) = Z_Q + (0.05 + j0.5) 121 ohm with
    # X_Q = Z_Q / sqrt(1 + 0.1^2), R_Q = 0.1 X_Q; line charging does not enter
    net2 = Net(name = "sc_two_bus", baseMVA = 100.0)
    addBus!(net = net2, busName = "B1", vn_kV = 110.0)
    addBus!(net = net2, busName = "B2", vn_kV = 110.0)
    addPIModelACLine!(net = net2, fromBus = "B1", toBus = "B2", r_pu = 0.05, x_pu = 0.5, b_pu = 0.2, status = 1)
    r2 = runShortCircuit!(net2, _ref_sc_data(enis = [_ref_sc_feeder(bus = "B1", ikmax_A = 10_000.0, rx_max = 0.1)]); case = :max)
    zq = 1.10 * 110_000.0 / (sqrt(3.0) * 10_000.0)
    xq = zq / sqrt(1.0 + 0.1^2)
    zk = (0.1 * xq + im * xq) + (0.05 + 0.5im) * 121.0
    row2 = only(row for row in r2.rows if row.bus == "B2")
    @test isapprox(row2.zk_ohm, abs(zk); rtol = 1e-9)
    @test isapprox(row2.ik_kA, 1.10 * 110.0 / (sqrt(3.0) * abs(zk)); rtol = 1e-9)
    @test isapprox(row2.rx_ratio, real(zk) / imag(zk); rtol = 1e-9)

    # synchronous machine: x''d 0.15 pu on 50 MVA -> 0.3 pu on 100 MVA,
    # R = 0.07 X (IEC 60909-0 6.6.3), Ik'' = 1.10 10.5 / (sqrt(3) |Z| Z_base)
    netm = Net(name = "sc_machine", baseMVA = 100.0)
    addBus!(net = netm, busName = "M1", vn_kV = 10.5)
    rm = runShortCircuit!(netm, _ref_sc_data(machines = [_ref_sc_machine(bus = "M1", xdpp = 0.15, ratedS = 50.0, ratedU = 10.5)]); case = :max)
    x_pu = 0.15 * (100.0 / 50.0)
    z_ohm = abs(0.07 * x_pu + im * x_pu) * (10.5^2 / 100.0)
    @test rm.rows[1].status === :ok
    @test isapprox(rm.rows[1].ik_kA, 1.10 * 10.5 / (sqrt(3.0) * z_ohm); rtol = 1e-9)

    # asynchronous machine: Z_M = (1/ilr) U^2 / S = (1/5) 10^2 / 5 = 4 ohm,
    # Ik'' = 1.10 10 / (sqrt(3) 4); S from P/(eta cos phi) = 3.6/(0.9 0.8) = 5 MVA
    neta = Net(name = "sc_motor", baseMVA = 100.0)
    addBus!(net = neta, busName = "M1", vn_kV = 10.0)
    ra = runShortCircuit!(neta, _ref_sc_data(asms = [_ref_sc_asm(bus = "M1", ilr = 5.0, rx = 0.1, ratedS = 5.0, ratedU = 10.0)]); case = :max)
    @test isapprox(ra.rows[1].zk_ohm, 4.0; rtol = 1e-9)
    @test isapprox(ra.rows[1].ik_kA, 1.10 * 10.0 / (sqrt(3.0) * 4.0); rtol = 1e-9)
    @test isapprox(ra.rows[1].rx_ratio, 0.1; rtol = 1e-9)
    rp = runShortCircuit!(neta, _ref_sc_data(asms = [_ref_sc_asm(bus = "M1", ilr = 5.0, rx = 0.1, ratedU = 10.0, p_mech = 3.6, efficiency = 90.0, cosphi = 0.8)]); case = :max)
    @test isapprox(rp.rows[1].ik_kA, ra.rows[1].ik_kA; rtol = 1e-12)

    # converter at B2 of the chain B1 - B2 - B3 (feeder at B1, 10 kA, every
    # impedance R/X 0.1): I_sk = 1.2 20 / (sqrt(3) 110) kA; fault at B1 or B2
    # takes all of it, at B3 the share |Z_Q + Z_L| / |Z_Q + 2 Z_L|
    netc = Net(name = "sc_pec", baseMVA = 100.0)
    for b in ("B1", "B2", "B3")
      addBus!(net = netc, busName = b, vn_kV = 110.0)
    end
    addPIModelACLine!(net = netc, fromBus = "B1", toBus = "B2", r_pu = 0.05, x_pu = 0.5, b_pu = 0.0, status = 1)
    addPIModelACLine!(net = netc, fromBus = "B2", toBus = "B3", r_pu = 0.05, x_pu = 0.5, b_pu = 0.0, status = 1)
    feeder = _ref_sc_feeder(bus = "B1", ikmax_A = 10_000.0, ikmin_A = 8_000.0, rx_max = 0.1, rx_min = 0.1)
    base = runShortCircuit!(netc, _ref_sc_data(enis = [feeder]); case = :max)
    with = runShortCircuit!(netc, _ref_sc_data(enis = [feeder], pecs = [_ref_sc_pec(bus = "B2", maxIFault = 1.2, ratedS = 20.0, ratedU = 110.0)]); case = :max)
    isk = 1.2 * 20.0 / (sqrt(3.0) * 110.0)
    zqc = 0.1 * xq + im * xq
    zl = (0.05 + 0.5im) * 121.0
    share = Dict("B1" => 1.0, "B2" => 1.0, "B3" => abs(zqc + zl) / abs(zqc + 2 * zl))
    for (b0, b1) in zip(base.rows, with.rows)
      @test isapprox(b1.ik_kA - b0.ik_kA, share[b1.bus] * isk; rtol = 1e-9)
      @test isapprox(b1.ip_kA, b0.ip_kA + sqrt(2.0) * share[b1.bus] * isk; rtol = 1e-9)
    end

    # source on a switch-only bus (ring A-B-D-A, C on B by a closed link):
    # Z_B = Z_QA + (2/3) Z_L, with the feeder at C Z_B' = Z_QC || Z_B,
    # Ik'' = c Un / (sqrt(3) |Z|); all R/X 0.1, so the feeder at C adds
    # exactly its declared current
    netl = Net(name = "sc_link_group", baseMVA = 100.0)
    for b in ("A", "B", "C", "D")
      addBus!(net = netl, busName = b, vn_kV = 110.0)
    end
    addPIModelACLine!(net = netl, fromBus = "A", toBus = "B", r_pu = 0.01, x_pu = 0.1, b_pu = 0.0, status = 1)
    addPIModelACLine!(net = netl, fromBus = "A", toBus = "D", r_pu = 0.01, x_pu = 0.1, b_pu = 0.0, status = 1)
    addPIModelACLine!(net = netl, fromBus = "B", toBus = "D", r_pu = 0.01, x_pu = 0.1, b_pu = 0.0, status = 1)
    addLink!(net = netl, fromBus = "B", toBus = "C", status = 1)
    markIsolatedBuses!(net = netl)
    fa = _ref_sc_feeder(bus = "A", ikmax_A = 10_000.0, ikmin_A = 8_000.0, rx_max = 0.1, rx_min = 0.1, name = "QA")
    fc = _ref_sc_feeder(bus = "C", ikmax_A = 5_000.0, ikmin_A = 4_000.0, rx_max = 0.1, rx_min = 0.1, name = "QC")
    zqf(c, ik_A) = (z = c * 110_000.0 / (sqrt(3.0) * ik_A); x = z / sqrt(1.0 + 0.1^2); 0.1 * x + im * x)
    zll = (0.01 + 0.1im) * 121.0
    for (case, c, ia, ic) in ((:max, 1.10, 10_000.0, 5_000.0), (:min, 1.00, 8_000.0, 4_000.0))
      zb = zqf(c, ia) + (2.0 / 3.0) * zll
      zb_c = 1.0 / (1.0 / zqf(c, ic) + 1.0 / zb)
      r0 = Dict(row.bus => row for row in runShortCircuit!(netl, _ref_sc_data(enis = [fa]); case = case).rows)
      r1 = Dict(row.bus => row for row in runShortCircuit!(netl, _ref_sc_data(enis = [fa, fc]); case = case).rows)
      @test isapprox(r0["B"].ik_kA, c * 110.0 / (sqrt(3.0) * abs(zb)); rtol = 1e-9)
      @test isapprox(r1["B"].ik_kA, c * 110.0 / (sqrt(3.0) * abs(zb_c)); rtol = 1e-9)
      @test isapprox(r1["B"].ik_kA - r0["B"].ik_kA, ic / 1000.0; rtol = 1e-9)
      @test r1["C"].ik_kA == r1["B"].ik_kA
    end
  end)() end
  @testset "hand calculation: two-bus feeder" begin (function ()
    net = Net(name = "ref_two_bus", baseMVA = 100.0)
    addBus!(net = net, busName = "A", vn_kV = 110.0)
    addBus!(net = net, busName = "B", vn_kV = 110.0)
    addProsumer!(net = net, busName = "A", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "A")
    addProsumer!(net = net, busName = "B", type = "ENERGYCONSUMER", p = 40.0, q = 10.0)
    addPIModelACLine!(net = net, fromBus = "A", toBus = "B", r_pu = 0.02, x_pu = 0.1, b_pu = 0.0, status = 1)
    ok, _ = validate!(net = net)
    @test ok
    _, erg = runpf!(net, 30, 1e-12, 0)
    @test erg == 0
    ref = _ref_two_bus_solution(0.02, 0.1, 0.4, 0.1)
    @test isapprox(net.nodeVec[2]._vm_pu, ref.v2; atol = 1e-10)
    @test isapprox(net.nodeVec[2]._va_deg, ref.d2_deg; atol = 1e-8)
    rows = _ref_node_rows(net)
    @test isapprox(rows[1].p_gen_MW, 100 * ref.p_src; atol = 1e-6)
    @test isapprox(rows[1].q_gen_MVar, 100 * ref.q_src; atol = 1e-6)
  end)() end

  @testset "hand calculation: phase-shifting transformer, MATPOWER shift convention" begin (function ()
    # MATPOWER manual: a branch with ratio 1 and shift phi has
    #   P_f = (V_f V_t / x) sin(theta_f - theta_t - phi),
    #   Q_f = (V_f^2 - V_f V_t cos(theta_f - theta_t - phi)) / x
    # and the same Q_t by symmetry. Bus 1 is the slack (V = 1), bus 2 a PV bus
    # with V = 1 that generates 30 MW, x = 0.1 pu, phi = 10 degrees, lossless.
    # The line carries -0.3 pu out of bus 1, so sin(d) = -0.03 with
    # d = theta_1 - theta_2 - phi, hence theta_2 = -phi + asin(0.03) and
    # Q_f = Q_t = 10 (1 - cos(asin 0.03)) pu.
    mktempdir() do dir
      file = joinpath(dir, "ref_pst.m")
      write(
        file,
        """
function mpc = ref_pst
mpc.version = '2';
mpc.baseMVA = 100;
mpc.bus = [
 1 3 0 0 0 0 1 1.0 0 110 1 1.1 0.9;
 2 2 0 0 0 0 1 1.0 0 110 1 1.1 0.9;
];
mpc.gen = [
 1 0 0 100 -100 1.0 100 1 200 -200;
 2 30 0 100 -100 1.0 100 1 100 0;
];
mpc.branch = [
 1 2 0 0.1 0 100 100 100 0 10 1 -360 360;
];
""",
      )
      net = createNetFromMatPowerFile(filename = file, flatstart = false, enable_pq_gen_controllers = true, bus_shunt_model = :admittance, matpower_shift_sign = 1.0, matpower_shift_unit = :deg, matpower_ratio = :normal, tap_changer_model = :ideal)
      _, erg = runpf!(net, 30, 1e-12, 0)
      @test erg == 0
      d = asin(0.03)
      @test isapprox(net.nodeVec[2]._va_deg, -10.0 + rad2deg(d); atol = 1e-8)
      rows = _ref_node_rows(net)
      @test isapprox(rows[1].p_gen_MW, -30.0; atol = 1e-8)
      @test isapprox(rows[2].p_gen_MW, 30.0; atol = 1e-8)
      q = 100 * 10 * (1 - cos(d))
      @test isapprox(rows[1].q_gen_MVar, q; atol = 1e-8)
      @test isapprox(rows[2].q_gen_MVar, q; atol = 1e-8)
    end
  end)() end

  @testset "hand calculation: three-winding transformer in the star equivalent" begin (function ()
    # add3WTPiModelTrafo! takes r and x of the three windings in ohm referred
    # to the highest rated voltage (the star point sits at that base), here
    # 380 kV. The load sits at the 110 kV winding, the 20 kV winding carries
    # nothing and so stands at the star point voltage. The chain from the
    # 380 kV bus to the load is one series impedance Z = Z_h + Z_m with the
    # two-bus solution; the star point is V_star = 1 - Z_h I.
    net = Net(name = "ref_3wt", baseMVA = 100.0)
    addBus!(net = net, busName = "H", vn_kV = 380.0)
    addBus!(net = net, busName = "M", vn_kV = 110.0)
    addBus!(net = net, busName = "L", vn_kV = 20.0)
    addProsumer!(net = net, busName = "H", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "H")
    add3WTPiModelTrafo!(net = net, HBBus = "H", MBBus = "M", LVBus = "L", r = (0.2, 0.3, 0.4), x = (4.0, 6.0, 10.0), b = (0.0, 0.0, 0.0), ratedU_kV = (380.0, 110.0, 20.0), ratedS_MVA = (300.0, 200.0, 50.0), status = 1)
    addProsumer!(net = net, busName = "M", type = "ENERGYCONSUMER", p = 40.0, q = 10.0)
    ok, _ = validate!(net = net)
    @test ok
    _, erg = runpf!(net, 30, 1e-12, 0)
    @test erg == 0
    zbase = 380.0^2 / 100.0
    zh = (0.2 + 4.0im) / zbase
    zm = (0.3 + 6.0im) / zbase
    ref = _ref_two_bus_solution(real(zh + zm), imag(zh + zm), 0.4, 0.1)
    v_m = ref.v2 * cis(deg2rad(ref.d2_deg))
    v_star = 1.0 - zh * conj((0.4 + 0.1im) / v_m)
    node = Dict(name => net.nodeVec[net.busDict[name]] for name in ("H", "M", "L"))
    @test isapprox(node["M"]._vm_pu, ref.v2; atol = 1e-10)
    @test isapprox(node["M"]._va_deg, ref.d2_deg; atol = 1e-8)
    @test isapprox(node["L"]._vm_pu, abs(v_star); atol = 1e-10)
    @test isapprox(node["L"]._va_deg, rad2deg(angle(v_star)); atol = 1e-8)
    rows = _ref_node_rows(net)
    @test isapprox(rows[net.busDict["H"]].p_gen_MW, 100 * ref.p_src; atol = 1e-6)
    @test isapprox(rows[net.busDict["H"]].q_gen_MVar, 100 * ref.q_src; atol = 1e-6)
  end)() end
  return true
end
