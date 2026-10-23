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
#
# file: test/test_short_circuit.jl
# purpose: short-circuit verification beyond the hand values of
#          fast/ref_hand_calculations.jl: bus selection and argument checks,
#          the safety-flag contract (default, warning, result flag),
#          asynchronous-machine lower-bound flagging, converter sources,
#          short_circuit.c_factor and converter_current_factor config
#          coverage, and the contraction of closed-link groups. This file
#          needs no CGMES fixtures.

using Test
using Sparlectra

const _SC = Sparlectra.CGMESImporter.CGMESShortCircuitData

_sc_feeder(; bus, ikmax_A = nothing, ikmin_A = nothing, rx_max = nothing, rx_min = nothing, name = "Feeder") =
  (mrid = name, name = name, bus = bus, maxInitialSymShCCurrent_A = ikmax_A, minInitialSymShCCurrent_A = ikmin_A, maxR1ToX1Ratio = rx_max, minR1ToX1Ratio = rx_min, maxR0ToX0Ratio = nothing, maxZ0ToZ1Ratio = nothing, ikSecond = nothing, governorSCD = nothing)

_sc_machine(; bus, xdpp = nothing, ratedS = nothing, ratedU = nothing, name = "G1") =
  (mrid = name, name = name, bus = bus, satDirectSubtransX_pu = xdpp, satDirectTransX_pu = nothing, r0_pu = nothing, x0_pu = nothing, r2_pu = nothing, x2_pu = nothing, earthing = nothing, ratedS_MVA = ratedS, ratedU_kV = ratedU)

_sc_asm(; bus, ilr = nothing, rx = nothing, ratedS = nothing, ratedU = nothing, p_mech = nothing, efficiency = nothing, cosphi = nothing, pole_pairs = nothing, name = "M1") =
  (mrid = name, name = name, bus = bus, iaIrRatio = ilr, rxLockedRotorRatio = rx, efficiency_percent = efficiency, ratedMechanicalPower_MW = p_mech, polePairNumber = pole_pairs, ratedS_MVA = ratedS, ratedU_kV = ratedU, ratedPowerFactor = cosphi)

_sc_pec(; bus, maxIFault = nothing, ratedS = nothing, ratedU = nothing, name = "PEC1") =
  (mrid = name, name = name, bus = bus, maxIFault_pu = maxIFault, ratedS_MVA = ratedS, ratedU_kV = ratedU)

_sc_data(; enis = NamedTuple[], machines = NamedTuple[], asms = NamedTuple[], pecs = NamedTuple[]) = _SC(enis, machines, NamedTuple[], NamedTuple[], NamedTuple[], asms, pecs)

# 110 kV feeder bus with one line to a second 110 kV bus.
function _sc_two_bus_net(; r_pu = 0.05, x_pu = 0.5)
  net = Net(name = "sc_two_bus", baseMVA = 100.0)
  addBus!(net = net, busName = "B1", vn_kV = 110.0)
  addBus!(net = net, busName = "B2", vn_kV = 110.0)
  # b_pu deliberately non-zero: IEC drops line charging, so it must not
  # change the result (asserted in fast/ref_hand_calculations.jl).
  addPIModelACLine!(net = net, fromBus = "B1", toBus = "B2", r_pu = r_pu, x_pu = x_pu, b_pu = 0.2, status = 1)
  return net
end

# 110 kV chain B1 - B2 - B3, both lines (0.05 + j0.5) pu, and the feeder at
# B1 (10 kA, R/X 0.1) of the converter tests.
function _sc_pec_chain()
  net = Net(name = "sc_pec", baseMVA = 100.0)
  for b in ("B1", "B2", "B3")
    addBus!(net = net, busName = b, vn_kV = 110.0)
  end
  addPIModelACLine!(net = net, fromBus = "B1", toBus = "B2", r_pu = 0.05, x_pu = 0.5, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "B2", toBus = "B3", r_pu = 0.05, x_pu = 0.5, b_pu = 0.0, status = 1)
  feeder = _sc_feeder(bus = "B1", ikmax_A = 10_000.0, ikmin_A = 8_000.0, rx_max = 0.1, rx_min = 0.1)
  return net, feeder
end

# multi-island SC fixture: n feeder-fed rings of m buses with declared Sk''
function _sc_sweep_net(n, m)
  net = Net(name = "sc_sweep_$(n)x$(m)", baseMVA = 100.0)
  for k in 1:n
    name = i -> "S$(k)_B$(i)"
    for i in 1:m
      addBus!(net = net, busName = name(i), vn_kV = 110.0)
    end
    addExternalGrid!(net = net, busName = name(1), vm_pu = 1.0, sk_max_MVA = 2000.0 + 100.0 * k, sk_min_MVA = 1500.0, rx_max = 0.1, internal_impedance = false)
    for i in 1:m
      addPIModelACLine!(net = net, fromBus = name(i), toBus = name(i == m ? 1 : i + 1), r_pu = 0.001, x_pu = 0.004, b_pu = 0.0, status = 1)
    end
  end
  ok, _ = validate!(net = net)
  @test ok
  return net
end

function run_short_circuit_tests()
  @testset "Short circuit (IEC 60909-0)" begin (function ()
    @testset "bus selection, argument checks and system frequency" begin (function ()
      # the hand values of this two-bus net are in fast/ref_hand_calculations.jl;
      # here only what the reference does not cover
      net = _sc_two_bus_net()
      sc = _sc_data(enis = [_sc_feeder(bus = "B1", ikmax_A = 10_000.0, rx_max = 0.1)])
      r = runShortCircuit!(net, sc; case = :max)
      row2 = only(row for row in r.rows if row.bus == "B2")
      @test row2.status === :ok
      # buses list variant selects a single bus and gives the all-bus value
      rb = runShortCircuit!(net, sc; buses = ["B2"], case = :max)
      @test length(rb.rows) == 1 && rb.rows[1].bus == "B2"
      @test isapprox(rb.rows[1].ik_kA, row2.ik_kA; rtol = 1e-12)
      @test_throws ArgumentError runShortCircuit!(net, sc; buses = ["NOPE"], case = :max)
      @test_throws ArgumentError runShortCircuit!(net, sc; case = :typical)
      # System frequency: Ik'' uses R and X as given (X already at the network
      # frequency), drops the line capacitance, and kappa by method b,
      # 1.02 + 0.98·e^(-3R/X), contains no f. A line from per-km capacitance
      # gets 6/5 times b at 60 Hz; Ik'' and kappa stay bit identical.
      sc_at(f) = begin
        n = Net(name = "sc_freq", baseMVA = 100.0, nominal_frequency_hz = f)
        addBus!(net = n, busName = "B1", vn_kV = 110.0)
        addBus!(net = n, busName = "B2", vn_kV = 110.0)
        addACLinePerKm!(net = n, fromBus = "B1", toBus = "B2", length_km = 30.0, r_ohm_per_km = 0.12, x_ohm_per_km = 0.39, c_nf_per_km = 9.5)
        (n, runShortCircuit!(n, sc; case = :max))
      end
      n50, r50 = sc_at(50.0)
      n60, r60 = sc_at(60.0)
      @test isapprox(n60.branchVec[1].b_pu / n50.branchVec[1].b_pu, 6 / 5; rtol = 1e-14)
      @test [row.ik_kA for row in r60.rows] == [row.ik_kA for row in r50.rows]
      @test [row.kappa for row in r60.rows] == [row.kappa for row in r50.rows]
      @test [row.ip_kA for row in r60.rows] == [row.ip_kA for row in r50.rows]
    end)() end

    @testset "safety-flag contract: default + warning + result flag" begin (function ()
      # A machine without x''d triggers all three: the documented default
      # (0.2 pu machine base), a log warning, and the per-row flag.
      net = Net(name = "sc_flag", baseMVA = 100.0)
      addBus!(net = net, busName = "M1", vn_kV = 10.5)
      sc = _sc_data(machines = [_sc_machine(bus = "M1", xdpp = nothing, ratedS = 50.0, ratedU = 10.5)])
      r = @test_logs (:warn, r"no usable x''_d, default 0.2 pu") match_mode = :any runShortCircuit!(net, sc; case = :max)
      row = r.rows[1]
      @test row.status === :ok
      @test row.contains_defaulted_data == true
      # a default reactance is flagged, but the machine feeds: no lower bound
      @test row.lower_bound == false
      @test any(occursin("default 0.2 pu", reason) for reason in row.reasons)
      @test any(occursin("default 0.2 pu", m) for m in r.messages)
      # documented default actually applied: x = 0.2·(100/50) = 0.4 pu
      x_pu = 0.2 * (100.0 / 50.0)
      z_ohm = abs(0.07 * x_pu + im * x_pu) * (10.5^2 / 100.0)
      @test isapprox(row.ik_kA, 1.10 * 10.5 / (sqrt(3.0) * z_ohm); rtol = 1e-9)
      # printer surfaces the flag and the reason
      out = sprint(io -> printShortCircuitResult(io, r))
      @test occursin("yes", out)
      @test occursin("default 0.2 pu", out)

      # The size the default refers to: rated power, else the maximum
      # active power (from the record, else from the unit of the network at
      # the bus), else the network base with a reason of its own. A stated
      # reactance without a rating stays on the network base: the maximum
      # is no base for a value the file gives.
      ik_on = (x_machine, size) -> begin
        x = x_machine * (100.0 / size)
        1.10 * 10.5 / (sqrt(3.0) * abs((size >= 100.0 ? 0.05 : 0.07) * x + im * x) * (10.5^2 / 100.0))
      end
      with_unit = Net(name = "sc_size", baseMVA = 100.0)
      addBus!(net = with_unit, busName = "M1", vn_kV = 10.5)
      addProsumer!(net = with_unit, busName = "M1", type = "SYNCHRONOUSMACHINE", p = 20.0, q = 0.0, pMax = 40.0)
      sized = record -> merge(_sc_machine(bus = "M1", ratedU = 10.5), record)
      for (label, network, record, size, reason) in [
        ("maximum in the record", net, sized((pmax_MW = 25.0,)), 25.0, "default 0.2 pu on its maximum active power 25.0 MW"),
        ("maximum of the unit at the bus", with_unit, sized(NamedTuple()), 40.0, "default 0.2 pu on its maximum active power 40.0 MW"),
        ("rated power before the maximum", with_unit, sized((ratedS_MVA = 50.0,)), 50.0, "default 0.2 pu (machine base)"),
        ("no size at all", net, sized(NamedTuple()), 100.0, "rests on defaults alone (0.2 pu on the network base 100.0 MVA)"),
        # 9999 MW is how a case file says that no limit is stated
        ("placeholder maximum", net, sized((pmax_MW = 9999.0,)), 100.0, "nor a usable maximum active power"),
        ("placeholder maximum with one digit more", net, sized((pmax_MW = 99999.0,)), 100.0, "nor a usable maximum active power"),
        # the same for the rated power: 9999 MVA and above is no size, the
        # search goes on with the maximum and the reason names the value;
        # sqrt(2)*99999.9 MVA is what a delivery with +-99999.9 limits carries
        ("placeholder rated power, maximum of the unit at the bus", with_unit, sized((ratedS_MVA = 141421.2,)), 40.0, "has no usable x''_d and a placeholder ratedS of 141421.2 MVA (9999.0 MVA and above is no size): default 0.2 pu on its maximum active power 40.0 MW"),
        ("placeholder rated power, no other size", net, sized((ratedS_MVA = 9999.0,)), 100.0, "has no x''_d, a placeholder ratedS of 9999.0 MVA (9999.0 MVA and above is no size) and no usable maximum active power"),
        ("rated power just below the placeholder", with_unit, sized((ratedS_MVA = 9998.0,)), 9998.0, "default 0.2 pu (machine base)"),
      ]
        result = run_with_expected_warnings(() -> runShortCircuit!(network, _sc_data(machines = [record]); case = :max), ["SynchronousMachine G1"])
        @testset "$(label)" begin
          @test isapprox(result.rows[1].ik_kA, ik_on(0.2, size); rtol = 1e-9)
          @test any(why -> occursin(reason, why), result.rows[1].reasons)
        end
      end
      stated = run_with_expected_warnings(() -> runShortCircuit!(with_unit, _sc_data(machines = [sized((satDirectSubtransX_pu = 0.3,))]); case = :max), ["has no usable ratedS"])
      @test isapprox(stated.rows[1].ik_kA, ik_on(0.3, 100.0); rtol = 1e-9)
      # a stated reactance next to a placeholder rated power: network base,
      # and the reason names the placeholder
      stated_ph = run_with_expected_warnings(() -> runShortCircuit!(with_unit, _sc_data(machines = [sized((satDirectSubtransX_pu = 0.3, ratedS_MVA = 141421.2))]); case = :max), ["has a placeholder ratedS"])
      @test isapprox(stated_ph.rows[1].ik_kA, ik_on(0.3, 100.0); rtol = 1e-9)
      @test "short-circuit: SynchronousMachine G1 has a placeholder ratedS of 141421.2 MVA (9999.0 MVA and above is no size): network base 100.0 MVA substituted" in stated_ph.rows[1].reasons
      # two islands: a placeholder machine next to a feeder (default
      # reactance, flagged, no lower bound) and a feeder without a current
      # next to one with it (skipped source, lower bound in the maximum case
      # only, where the engine sets the typed flag on every row of the island)
      two = Net(name = "sc_two_islands", baseMVA = 100.0)
      foreach(b -> addBus!(net = two, busName = b, vn_kV = 110.0), ("A", "B", "C", "D"))
      addPIModelACLine!(net = two, fromBus = "A", toBus = "B", r_pu = 0.01, x_pu = 0.05, b_pu = 0.0, status = 1)
      addPIModelACLine!(net = two, fromBus = "C", toBus = "D", r_pu = 0.01, x_pu = 0.05, b_pu = 0.0, status = 1)
      two_sc = _sc_data(enis = [_sc_feeder(bus = "A", ikmax_A = 10_000.0, ikmin_A = 10_000.0, rx_max = 0.1, rx_min = 0.1, name = "F_A"), _sc_feeder(bus = "C", ikmax_A = 8_000.0, ikmin_A = 8_000.0, rx_max = 0.1, rx_min = 0.1, name = "F_C"), _sc_feeder(bus = "D", name = "F_D")], machines = [_sc_machine(bus = "B", ratedS = 141421.2, ratedU = 110.0, name = "G_PH")])
      for case in (:max, :min)
        two_res = run_with_expected_warnings(() -> runShortCircuit!(two, two_sc; case = case), ["G_PH", "F_D"])
        rows = Dict(row.bus => row for row in two_res.rows)
        @test all(rows[b].contains_defaulted_data for b in ("A", "B", "C", "D"))
        @test [rows[b].lower_bound for b in ("A", "B", "C", "D")] == (case === :max ? [false, false, true, true] : [false, false, false, false])
        # two of four sources without data: the run status is a warning, and
        # the printed table starts with it and the two counts
        @test (two_res.status, two_res.status_reason) == (:warning, "short_circuit_partial_defaults")
        @test first(split(sprint(io -> printShortCircuitResult(io, two_res)), '\n')) == "run status: warning (short_circuit_partial_defaults), 4 of 4 rows flagged, $(case === :max ? 2 : 0) rows with a lower bound"
      end
      # every source with its data: succeeded without a reason; no source at
      # all: no bus can be evaluated, failed as in the service
      clean = runShortCircuit!(two, _sc_data(enis = [_sc_feeder(bus = "A", ikmax_A = 10_000.0, rx_max = 0.1, name = "F_A"), _sc_feeder(bus = "C", ikmax_A = 8_000.0, rx_max = 0.1, name = "F_C")]); case = :max)
      @test (clean.status, clean.status_reason) == (:succeeded, nothing)
      @test first(split(sprint(io -> printShortCircuitResult(io, clean)), '\n')) == "run status: succeeded, 0 of 4 rows flagged, 0 rows with a lower bound"
      empty_res = runShortCircuit!(two, _sc_data(); case = :max)
      @test (empty_res.status, empty_res.status_reason) == (:failed, "short_circuit_data_missing")
      # the service outcome reads the typed flag: rows on a default alone
      # (a feeder without its R/X ratio) are no lower bound
      rx_default = run_with_expected_warnings(() -> runShortCircuit!(two, _sc_data(enis = [_sc_feeder(bus = "A", ikmax_A = 10_000.0, name = "F_A"), _sc_feeder(bus = "C", ikmax_A = 8_000.0, rx_max = 0.1, name = "F_C")]); case = :max), ["no usable R/X ratio"])
      ok_state = (total = 2, with_data = 2, without_any = 0, refused = false, defaults_only = false, partial = false, warning = false)
      @test SparlectraApp._sc_outcome(ok_state, rx_default, "t.csv").reason == "short_circuit_flagged_defaults"
      # the service counts a machine with only a placeholder rated power as
      # a source without any data, one with a real rating as defaultable
      @test SparlectraApp._sc_source_data_state(_sc_data(machines = [sized((ratedS_MVA = 141421.2,))])).refused
      @test !SparlectraApp._sc_source_data_state(_sc_data(machines = [sized((ratedS_MVA = 50.0,))])).refused
      # a machine on a wind or solar unit is a current source, not counted
      # as a machine, as a PowerElectronicsConnection is not
      @test SparlectraApp._sc_source_data_state(_sc_data(machines = [sized((ratedS_MVA = 141421.2, converter_unit = "WindGeneratingUnit"))])).total == 0
    end)() end

    @testset "asynchronous machine: default R/X, missing data, minimum case" begin (function ()
      # the hand motor impedance Z_M and S from P/(eta cos phi) are in
      # fast/ref_hand_calculations.jl
      net = Net(name = "sc_motor", baseMVA = 100.0)
      addBus!(net = net, busName = "M1", vn_kV = 10.0)

      # Missing locked-rotor R/X: §6.7.2 default (0.10 — MV, ≥ 1 MW per pole
      # pair) is substituted AND flagged; the motor still contributes.
      sc_rx = _sc_data(asms = [_sc_asm(bus = "M1", ilr = 5.0, ratedS = 5.0, ratedU = 10.0, p_mech = 4.0, pole_pairs = 2.0)])
      rrx = @test_logs (:warn, r"no usable rxLockedRotorRatio") match_mode = :any runShortCircuit!(net, sc_rx; case = :max)
      @test rrx.rows[1].status === :ok
      @test rrx.rows[1].contains_defaulted_data == true
      @test isapprox(rrx.rows[1].rx_ratio, 0.10; rtol = 1e-9)

      # Incomplete |Z_M| data (no iaIrRatio): skipped with the island-wide
      # lower-bound flag — never a silent guess.
      sc_bad = _sc_data(asms = [_sc_asm(bus = "M1", ratedS = 5.0, ratedU = 10.0)])
      rbad = @test_logs (:warn, r"no usable iaIrRatio") match_mode = :any runShortCircuit!(net, sc_bad; case = :max)
      @test rbad.rows[1].status === :no_source   # motor was the only source
      @test rbad.rows[1].contains_defaulted_data == true
      @test any(occursin("lower bound", reason) for reason in rbad.rows[1].reasons)

      # Motors never contribute to the minimum case: with a feeder present
      # the max/min gap at the motor bus exceeds the pure c-factor ratio.
      net2 = _sc_two_bus_net()
      sc2 = _sc_data(enis = [_sc_feeder(bus = "B1", ikmax_A = 10_000.0, ikmin_A = 10_000.0, rx_max = 0.1, rx_min = 0.1)], asms = [_sc_asm(bus = "B2", ilr = 5.0, rx = 0.1, ratedS = 20.0, ratedU = 110.0)])
      r2max = runShortCircuit!(net2, sc2; case = :max)
      r2min = runShortCircuit!(net2, sc2; case = :min)
      b2max = only(row for row in r2max.rows if row.bus == "B2")
      b2min = only(row for row in r2min.rows if row.bus == "B2")
      # identical declared feeder currents → without the motor the ratio would
      # be exactly c_max/c_min = 1.10; the motor lifts only the max case
      @test b2max.ik_kA / b2min.ik_kA > 1.2
      @test b2max.contains_defaulted_data == false && b2min.contains_defaulted_data == false
    end)() end

    @testset "asynchronous machine: skipped with island-wide lower-bound flag (max only)" begin (function ()
      net = _sc_two_bus_net()
      addProsumer!(net = net, busName = "B2", type = "ASYNCHRONOUSMACHINE", p = 5.0, q = 2.0)
      # same component-type marker the CGMES importer sets for motors (the
      # prosumer constructor collapses the type string to a generic Load)
      net.prosumpsVec[end].comp.cTyp = Sparlectra.AsynchronousMachine
      sc = _sc_data(enis = [_sc_feeder(bus = "B1", ikmax_A = 10_000.0, ikmin_A = 8_000.0, rx_max = 0.1, rx_min = 0.1)])
      rmax = @test_logs (:warn, r"AsynchronousMachine .* lower bound") match_mode = :any runShortCircuit!(net, sc; case = :max)
      # the whole island is a lower bound, not just the motor bus
      @test all(row.contains_defaulted_data for row in rmax.rows)
      @test all(row.lower_bound for row in rmax.rows)
      ok_state = (total = 1, with_data = 1, without_any = 0, refused = false, defaults_only = false, partial = false, warning = false)
      @test SparlectraApp._sc_outcome(ok_state, rmax, "t.csv").reason == "short_circuit_flagged_lower_bound"
      @test all(any(occursin("lower bound", reason) for reason in row.reasons) for row in rmax.rows)
      # motors never contribute to (and never flag) the minimum case
      rmin = runShortCircuit!(net, sc; case = :min)
      @test all(!row.contains_defaulted_data for row in rmin.rows)
      @test !any(row.lower_bound for row in rmin.rows)
    end)() end

    @testset "converter current source per hand-derived superposition" begin (function ()
      # Hand derivation on a 110 kV chain B1 - B2 - B3, feeder at B1 (10 kA,
      # R/X 0.1), both lines (0.05 + j0.5)·121 Ω (R/X 0.1 as well, so every
      # impedance has one angle and the magnitudes add like the phasors).
      # Converter at B2: I_sk = maxIFault·ratedS/(√3·ratedU)
      #   = 1.2·20/(√3·110) kA. With the feeder short-circuited it sees
      # Z_Q + Z_L towards B1 and Z_L + fault towards B3:
      #   fault at B1 or B2: all of I_sk flows into the fault;
      #   fault at B3: |Z_32|/|Z_33| = |Z_Q + Z_L|/|Z_Q + 2·Z_L| of it.
      net, feeder = _sc_pec_chain()
      base = runShortCircuit!(net, _sc_data(enis = [feeder]); case = :max)
      with = runShortCircuit!(net, _sc_data(enis = [feeder], pecs = [_sc_pec(bus = "B2", maxIFault = 1.2, ratedS = 20.0, ratedU = 110.0)]); case = :max)
      isk = 1.2 * 20.0 / (sqrt(3.0) * 110.0)
      zq = 1.10 * 110_000.0 / (sqrt(3.0) * 10_000.0)
      xq = zq / sqrt(1.0 + 0.1^2)
      zqc = 0.1 * xq + im * xq
      zl = (0.05 + 0.5im) * 121.0
      share = Dict("B1" => 1.0, "B2" => 1.0, "B3" => abs(zqc + zl) / abs(zqc + 2 * zl))
      for (b0, b1) in zip(base.rows, with.rows)
        @test b1.status === :ok && !b1.contains_defaulted_data
        @test isapprox(b1.ik_kA - b0.ik_kA, share[b1.bus] * isk; rtol = 1e-9)
        @test isapprox(b1.sk_MVA, sqrt(3.0) * 110.0 * b1.ik_kA; rtol = 1e-12)
        # Z_k, R/X and kappa stay those of the voltage source; the converter
        # part enters the peak without kappa
        @test (b1.zk_ohm, b1.rx_ratio, b1.kappa) == (b0.zk_ohm, b0.rx_ratio, b0.kappa)
        @test isapprox(b1.ip_kA, b0.ip_kA + sqrt(2.0) * share[b1.bus] * isk; rtol = 1e-9)
      end
      # the all-bus Takahashi sweep agrees to machine precision (the threaded
      # sweep against the serial one: run_short_circuit_parallel_tests)
      pec_data = _sc_data(enis = [feeder], pecs = [_sc_pec(bus = "B2", maxIFault = 1.2, ratedS = 20.0, ratedU = 110.0)])
      serial = runShortCircuit!(net, pec_data; case = :max, parallel_enabled = false)
      tak = runShortCircuit!(net, pec_data; case = :max, sweep_method = :takahashi, parallel_enabled = false)
      @test all(isapprox(a.ik_kA, b.ik_kA; rtol = 1e-10) for (a, b) in zip(serial.rows, tak.rows))
      # the rated current comes from ratedU of the unit; without it the bus
      # nominal voltage stands in, and the row says so
      for (ru, flagged) in ((115.0, false), (nothing, true))
        r = run_with_expected_warnings(() -> runShortCircuit!(net, _sc_data(enis = [feeder], pecs = [_sc_pec(bus = "B2", maxIFault = 1.2, ratedS = 20.0, ratedU = ru)]); buses = ["B2"], case = :max), ["PowerElectronicsConnection PEC1 has no usable ratedU"])
        @test isapprox(r.rows[1].ik_kA - base.rows[2].ik_kA, 1.2 * 20.0 / (sqrt(3.0) * something(ru, 110.0)); rtol = 1e-9)
        @test r.rows[1].contains_defaulted_data == flagged
      end
      # the minimum case: no converter share and no flag
      bmin = runShortCircuit!(net, _sc_data(enis = [feeder]); case = :min)
      wmin = runShortCircuit!(net, _sc_data(enis = [feeder], pecs = [_sc_pec(bus = "B2", maxIFault = 1.2, ratedS = 20.0, ratedU = 110.0)]); case = :min)
      @test isequal(bmin.rows, wmin.rows)
      # without maxIFault the configured factor stands in: k times the rated
      # current with the same shares, each row naming the substitute and none
      # the lower bound; a stated maxIFault wins over the factor
      k = 1.5
      fac = run_with_expected_warnings(() -> runShortCircuit!(net, _sc_data(enis = [feeder], pecs = [_sc_pec(bus = "B2", ratedS = 20.0, ratedU = 110.0)]); case = :max, converter_current_factor = k), ["PowerElectronicsConnection PEC1 has no maxIFault: short_circuit.converter_current_factor 1.5 times the rated current substituted"])
      for (b0, b1) in zip(base.rows, fac.rows)
        @test isapprox(b1.ik_kA - b0.ik_kA, share[b1.bus] * k * 20.0 / (sqrt(3.0) * 110.0); rtol = 1e-9)
        @test b1.contains_defaulted_data && any(why -> occursin("converter_current_factor 1.5 times the rated current substituted", why), b1.reasons)
        @test !any(why -> occursin("lower bound", why), b1.reasons)
        @test !b1.lower_bound
      end
      @test isequal(runShortCircuit!(net, pec_data; case = :max, converter_current_factor = k).rows, with.rows)
      # only an absent maxIFault takes the factor; a stated one that is no
      # positive current feeds nothing and flags the island with its value
      zero = run_with_expected_warnings(() -> runShortCircuit!(net, _sc_data(enis = [feeder], pecs = [_sc_pec(bus = "B2", maxIFault = 0.0, ratedS = 20.0, ratedU = 110.0)]); case = :max, converter_current_factor = k), ["PowerElectronicsConnection PEC1 has maxIFault 0.0, not a usable current"])
      @test [r.ik_kA for r in zero.rows] == [r.ik_kA for r in base.rows]
      @test all(any(why -> occursin("PEC1 has maxIFault 0.0, not a usable current: contribution skipped; Ik''max in its island is a lower bound", why), r.reasons) for r in zero.rows)
      @test all(r -> r.lower_bound, zero.rows)
      # the service outcome: rows on the substitute are no lower bound
      ok_state = (total = 1, with_data = 1, without_any = 0, refused = false, defaults_only = false, partial = false, warning = false)
      @test SparlectraApp._sc_outcome(ok_state, fac, "t.csv").reason == "short_circuit_flagged_defaults"
      @test SparlectraApp._sc_outcome(ok_state, zero, "t.csv").reason == "short_circuit_flagged_lower_bound"
      @test_throws ArgumentError runShortCircuit!(net, pec_data; converter_current_factor = -1.0)
      # without maxIFault, without ratedS, and a converter prosumer without
      # a record: no share, the island flagged as a lower bound
      addProsumer!(net = net, busName = "B3", type = "GENERATOR", p = 5.0, q = 0.0)
      net.prosumpsVec[end].comp.cTyp = Sparlectra.PowerElectronicsConnection
      for (label, pecs, reason) in [
        ("no maxIFault", [_sc_pec(bus = "B2", ratedS = 20.0, ratedU = 110.0), _sc_pec(bus = "B3", maxIFault = 1.0, ratedS = 5.0, ratedU = 110.0, name = "PEC3")], "PowerElectronicsConnection at bus B2: no converter short-circuit model, contribution skipped; Ik''max in its island is a lower bound"),
        ("no ratedS", [_sc_pec(bus = "B2", maxIFault = 1.2, ratedU = 110.0), _sc_pec(bus = "B3", maxIFault = 1.0, ratedS = 5.0, ratedU = 110.0, name = "PEC3")], "PowerElectronicsConnection PEC1 has maxIFault but no usable ratedS"),
        ("prosumer without a record", NamedTuple[], "PowerElectronicsConnection at bus B3: no converter short-circuit model"),
      ]
        r = run_with_expected_warnings(() -> runShortCircuit!(net, _sc_data(enis = [feeder], pecs = pecs); case = :max), ["short-circuit: PowerElectronicsConnection"])
        @testset "$(label)" begin
          @test all(row.contains_defaulted_data && any(why -> occursin(reason, why), row.reasons) for row in r.rows)
          # only the recorded, usable converter at B3 adds its current
          b3 = isempty(pecs) ? 0.0 : 5.0 / (sqrt(3.0) * 110.0)
          @test isapprox(r.rows[3].ik_kA - base.rows[3].ik_kA, b3; rtol = 1e-9)
        end
      end
    end)() end

    @testset "island without a source reports :no_source, isolated buses :isolated" begin (function ()
      net = Net(name = "sc_islands", baseMVA = 100.0)
      addBus!(net = net, busName = "A1", vn_kV = 110.0)
      addBus!(net = net, busName = "A2", vn_kV = 110.0)
      addBus!(net = net, busName = "C1", vn_kV = 110.0)
      addPIModelACLine!(net = net, fromBus = "A1", toBus = "A2", r_pu = 0.01, x_pu = 0.1, b_pu = 0.0, status = 1)
      sc = _sc_data(enis = [_sc_feeder(bus = "A1", ikmax_A = 5_000.0, rx_max = 0.1)])
      r = runShortCircuit!(net, sc; case = :max)
      by_bus = Dict(row.bus => row for row in r.rows)
      @test by_bus["A1"].status === :ok
      @test by_bus["A2"].status === :ok
      @test by_bus["C1"].status === :no_source
      @test isnan(by_bus["C1"].ik_kA)
      @test by_bus["C1"].contains_defaulted_data == true
      @test any(occursin("no short-circuit source", reason) for reason in by_bus["C1"].reasons)
      # an island of converter units alone:
      # a PowerElectronicsConnection without maxIFault at C1 (20 MVA) and a
      # machine on a wind unit at C2 (30 MVA), both 110 kV, k = 1.2. No
      # voltage source, so every bus gets the arithmetic sum of the
      # converter currents, Ik'' = k (S1 + S2) / (sqrt(3) 110 kV)
      # = 1.2 * 50 / 190.526 = 0.31492 kA, ip = sqrt(2) Ik'' without kappa,
      # flagged, no lower bound; with the factor off nothing feeds and the
      # island stays :no_source
      addBus!(net = net, busName = "C2", vn_kV = 110.0)
      addPIModelACLine!(net = net, fromBus = "C1", toBus = "C2", r_pu = 0.01, x_pu = 0.1, b_pu = 0.0, status = 1)
      conv = _sc_data(enis = [_sc_feeder(bus = "A1", ikmax_A = 5_000.0, rx_max = 0.1)], pecs = [_sc_pec(bus = "C1", ratedS = 20.0, ratedU = 110.0)], machines = [merge(_sc_machine(bus = "C2", ratedS = 30.0, ratedU = 110.0, name = "W1"), (converter_unit = "WindGeneratingUnit",))])
      ik_hand = 1.2 * (20.0 + 30.0) / (sqrt(3.0) * 110.0)
      on = Dict(row.bus => row for row in run_with_expected_warnings(() -> runShortCircuit!(net, conv; case = :max, converter_current_factor = 1.2), ["short-circuit: "]).rows)
      for b in ("C1", "C2")
        @test on[b].status === :ok
        @test isapprox(on[b].ik_kA, ik_hand; rtol = 1e-12)
        @test isapprox(on[b].ip_kA, sqrt(2.0) * ik_hand; rtol = 1e-12)
        @test on[b].contains_defaulted_data && !on[b].lower_bound
        @test "short-circuit: converter-only island, arithmetic sum, no voltage source" in on[b].reasons
      end
      @test on["A1"].ik_kA == by_bus["A1"].ik_kA
      off = Dict(row.bus => row for row in run_with_expected_warnings(() -> runShortCircuit!(net, conv; case = :max), ["short-circuit: "]).rows)
      @test all(off[b].status === :no_source && off[b].lower_bound for b in ("C1", "C2"))
    end)() end

    @testset "closed-link groups are one node, as in the power flow" begin (function ()
      # Ring A-B-D-A (ring = false drops B-D) with bus C on B through a
      # closed link only, plus E-F: two buses joined by a closed link and no
      # branch. Feeder at A (10 kA max, 8 kA min), optional feeder at C
      # (5 kA max, 4 kA min), R/X 0.1 for both.
      function link_net(ring::Bool)
        net = Net(name = "sc_link_group", baseMVA = 100.0)
        for b in ("A", "B", "C", "D", "E", "F")
          addBus!(net = net, busName = b, vn_kV = 110.0)
        end
        addPIModelACLine!(net = net, fromBus = "A", toBus = "B", r_pu = 0.01, x_pu = 0.1, b_pu = 0.0, status = 1)
        addPIModelACLine!(net = net, fromBus = "A", toBus = "D", r_pu = 0.01, x_pu = 0.1, b_pu = 0.0, status = 1)
        ring && addPIModelACLine!(net = net, fromBus = "B", toBus = "D", r_pu = 0.01, x_pu = 0.1, b_pu = 0.0, status = 1)
        addLink!(net = net, fromBus = "B", toBus = "C", status = 1)
        addLink!(net = net, fromBus = "E", toBus = "F", status = 1)
        # the marking the power flow sees: C, E and F have no branch
        markIsolatedBuses!(net = net)
        return net
      end
      feeder_a = _sc_feeder(bus = "A", ikmax_A = 10_000.0, ikmin_A = 8_000.0, rx_max = 0.1, rx_min = 0.1, name = "QA")
      feeder_c = _sc_feeder(bus = "C", ikmax_A = 5_000.0, ikmin_A = 4_000.0, rx_max = 0.1, rx_min = 0.1, name = "QC")
      rows(net, enis, case; sweep = :solves) = Dict(row.bus => row for row in runShortCircuit!(net, _sc_data(enis = enis); case = case, sweep_method = sweep).rows)

      @testset "switch-only bus reports the Ik'' of its group ($(ring ? "ring" : "radial"))" for ring in (true, false)
        net = link_net(ring)
        @test net.isoNodes == sort([net.busDict[b] for b in ("C", "E", "F")])
        for case in (:max, :min)
          r = rows(net, [feeder_a], case)
          @test r["C"].status === :ok
          @test r["C"].ik_kA == r["B"].ik_kA
          @test r["C"].island == r["B"].island
          @test isfinite(r["C"].ik_kA)
          ring || continue
          # a source on the switch-only bus C: the hand Ik'' with and without
          # it is in fast/ref_hand_calculations.jl. Here the impedance of the
          # contracted node, Z_B' = Z_QC parallel (Z_QA + (2/3) Z_L) with
          # Z_Q = c Un/(sqrt(3) Ik_declared) at R/X 0.1 and Z_L = (0.01 + j0.1) 121 ohm,
          # and the selected-inverse sweep on the same contracted matrix
          c, ia, ic = case === :max ? (1.10, 10_000.0, 5_000.0) : (1.00, 8_000.0, 4_000.0)
          zq = ik_A -> (z = c * 110_000.0 / (sqrt(3.0) * ik_A); x = z / sqrt(1.0 + 0.1^2); 0.1 * x + im * x)
          zb_c = 1.0 / (1.0 / zq(ic) + 1.0 / (zq(ia) + (2.0 / 3.0) * (0.01 + 0.1im) * 121.0))
          r1 = rows(net, [feeder_a, feeder_c], case)
          @test isapprox(r1["B"].zk_ohm, abs(zb_c); rtol = 1e-9)
          @test r1["C"].ik_kA == r1["B"].ik_kA
          rt = rows(net, [feeder_a, feeder_c], case; sweep = :takahashi)
          @test all(isapprox(rt[b].ik_kA, r1[b].ik_kA; rtol = 1e-10) for b in ("A", "B", "C", "D"))
        end
      end

      @testset "group without a branch stays isolated" begin
        net = link_net(true)
        r = rows(net, [feeder_a, _sc_feeder(bus = "E", ikmax_A = 5_000.0, rx_max = 0.1, name = "QE")], :max)
        for b in ("E", "F")
          @test r[b].status === :isolated
          @test isnan(r[b].ik_kA)
        end
        # one group, one island
        @test r["E"].island == r["F"].island
        @test r["E"].island != r["A"].island
      end
    end)() end

    @testset "WebUI data-check gating helper" begin (function ()
      # Byte scan over the delivery contents (fast fixture: plain XML folder;
      # the real ZIP path goes through the same collectCGMESFiles reader and
      # is exercised by the extended CGMES service test).
      # The gate follows the rule of the run: a machine counts with its
      # reactance or at least its rated power, not by its mere presence.
      with_data = mktempdir()
      write(joinpath(with_data, "grid_EQ.xml"), "<rdf:RDF><cim:SynchronousMachine rdf:ID=\"g\"><cim:RotatingMachine.ratedS>100</cim:RotatingMachine.ratedS></cim:SynchronousMachine></rdf:RDF>")
      without_data = mktempdir()
      write(joinpath(without_data, "grid_EQ.xml"), "<rdf:RDF><cim:EnergyConsumer rdf:ID=\"l\"/></rdf:RDF>")
      bare_machine = mktempdir()
      write(joinpath(bare_machine, "grid_EQ.xml"), "<rdf:RDF><cim:SynchronousMachine rdf:ID=\"g\"/></rdf:RDF>")
      @test SparlectraApp._webui_case_has_short_circuit_data(with_data) == true
      @test SparlectraApp._webui_case_has_short_circuit_data(without_data) == false
      @test SparlectraApp._webui_case_has_short_circuit_data(bare_machine) == false
      # cached second call returns the same verdict
      @test SparlectraApp._webui_case_has_short_circuit_data(without_data) == false
      # unresolvable paths must NOT lock the button (service explains instead)
      @test SparlectraApp._webui_case_has_short_circuit_data(joinpath(without_data, "missing.zip")) == true
    end)() end

    @testset "short_circuit.c_factor config coverage" begin (function ()
      @test Sparlectra.ShortCircuitConfig().c_factor == 0.0
      ok_yaml = test_scratch_path(".yaml")
      write(ok_yaml, "config_version: 1\nshort_circuit:\n  c_factor: 1.05\n")
      @test Sparlectra.load_sparlectra_config(ok_yaml; reload = true).shortcircuit.c_factor == 1.05
      bad_yaml = test_scratch_path(".yaml")
      write(bad_yaml, "config_version: 1\nshort_circuit:\n  c_factor: 1.4\n")
      err = @test_throws ArgumentError Sparlectra.load_sparlectra_config(bad_yaml; reload = true)
      @test occursin("short_circuit.c_factor", sprint(showerror, err.value))
      # the run keyword: a scalar c factor replaces the IEC table (1.10 for
      # the maximum case at 10.5 kV) linearly, 1.4 is out of range
      netm = Net(name = "sc_machine", baseMVA = 100.0)
      addBus!(net = netm, busName = "M1", vn_kV = 10.5)
      scm = _sc_data(machines = [_sc_machine(bus = "M1", xdpp = 0.15, ratedS = 50.0, ratedU = 10.5)])
      ik_table = runShortCircuit!(netm, scm; case = :max).rows[1].ik_kA
      @test isapprox(runShortCircuit!(netm, scm; case = :max, c_factor = 1.0).rows[1].ik_kA, ik_table / 1.10; rtol = 1e-9)
      @test_throws ArgumentError runShortCircuit!(netm, scm; c_factor = 1.4)
      # converter current factor: off by default, a positive multiple
      # accepted, a negative one rejected with the key named
      @test Sparlectra.ShortCircuitConfig().converter_current_factor == 0.0
      write(ok_yaml, "config_version: 1\nshort_circuit:\n  converter_current_factor: 1.2\n")
      @test Sparlectra.load_sparlectra_config(ok_yaml; reload = true).shortcircuit.converter_current_factor == 1.2
      write(bad_yaml, "config_version: 1\nshort_circuit:\n  converter_current_factor: -0.5\n")
      bad = @test_throws ArgumentError Sparlectra.load_sparlectra_config(bad_yaml; reload = true)
      @test occursin("short_circuit.converter_current_factor", sprint(showerror, bad.value))
      # a case key: the case configuration file and an API override carry it
      case_m = test_scratch_path(".m")
      write(Sparlectra.case_config_path(case_m), "config_version: 1\nscope: case\ncase: $(basename(case_m))\nshort_circuit:\n  converter_current_factor: 1.5\n")
      @test Sparlectra.resolve_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, case_m).config.shortcircuit.converter_current_factor == 1.5
      @test Sparlectra.validate_gui_config_overrides(Dict{String,Any}("short_circuit.converter_current_factor" => 1.5)) == Dict{String,Any}("short_circuit" => Dict{String,Any}("converter_current_factor" => 1.5))
      @test_throws ArgumentError Sparlectra.validate_gui_config_overrides(Dict{String,Any}("short_circuit.converter_current_factor" => -1.0))
    end)() end

    run_short_circuit_parallel_tests()
  end)() end
end

# The parallel sweeps against the serial sweep. run_short_circuit_tests runs
# them in the pf profile; the threads profile (tools/run_gates.sh threads)
# repeats them with four threads.
function run_short_circuit_parallel_tests()
  @testset "converter sweep: parallel equals serial" begin (function ()
    # the threaded sweep only reads the converter sums formed before the
    # spawn and gives the serial rows bitwise, over repeated runs
    net, feeder = _sc_pec_chain()
    pec_data = _sc_data(enis = [feeder], pecs = [_sc_pec(bus = "B2", maxIFault = 1.2, ratedS = 20.0, ratedU = 110.0)])
    serial = runShortCircuit!(net, pec_data; case = :max, parallel_enabled = false)
    for _ in 1:5
      par = runShortCircuit!(net, pec_data; case = :max, parallel_enabled = true, parallel_max_tasks = Threads.nthreads(), parallel_min_work_items = 2)
      @test isequal(par.rows, serial.rows)
    end
    parallel_comparison_report("sc converter sweep (5 repeats)")
  end)() end

  @testset "parallel all-bus sweep identity (Phase 3)" begin (function ()
    net = _sc_sweep_net(3, 40)
    for case in (:max, :min)
      serial = runShortCircuit!(net; case = case, parallel_enabled = false)
      # row-by-row equality, == on all fields, NaN-aware via isequal,
      # for the capped-serial fallback AND the genuinely parallel sweep
      for max_tasks in (1, Threads.nthreads())
        par = runShortCircuit!(net; case = case, parallel_enabled = true, parallel_max_tasks = max_tasks, parallel_min_work_items = 2)
        @test length(par.rows) == length(serial.rows)
        @test isequal(serial.rows, par.rows)
      end
    end
    parallel_comparison_report("sc all-bus sweep, 3 rings of 40 buses")

    # sweep_method = :takahashi: one selected-inverse pass per island;
    # must agree with :solves to machine precision (documented: not
    # bitwise) with identical statuses/flags, serial and parallel alike
    @test_throws ArgumentError runShortCircuit!(net; sweep_method = :bogus)
    @test_throws ArgumentError runShortCircuit!(net; takahashi_min_buses = 0)
    solves_ref = runShortCircuit!(net; case = :max, sweep_method = :solves, parallel_enabled = false)
    # sweep_method = :auto (the default): below the island threshold it is
    # bitwise the solves sweep, with threshold 1 bitwise the takahashi
    # sweep; the strategy choice itself never changes a row
    auto_small = runShortCircuit!(net; case = :max, sweep_method = :auto, takahashi_min_buses = 10_000, parallel_enabled = false)
    @test isequal(auto_small.rows, solves_ref.rows)
    auto_tak = runShortCircuit!(net; case = :max, sweep_method = :auto, takahashi_min_buses = 1, parallel_enabled = false)
    tak_ref = runShortCircuit!(net; case = :max, sweep_method = :takahashi, parallel_enabled = false)
    @test isequal(auto_tak.rows, tak_ref.rows)
    for parallel in (false, true)
      tak = runShortCircuit!(net; case = :max, sweep_method = :takahashi, parallel_enabled = parallel, parallel_min_work_items = 2)
      @test length(tak.rows) == length(solves_ref.rows)
      for i in eachindex(solves_ref.rows)
        a, b = solves_ref.rows[i], tak.rows[i]
        @test a.status === b.status
        @test a.contains_defaulted_data == b.contains_defaulted_data
        if a.status === :ok
          @test isapprox(a.ik_kA, b.ik_kA; rtol = 1e-10)
          @test isapprox(a.zk_ohm, b.zk_ohm; rtol = 1e-10)
          @test isapprox(a.sk_MVA, b.sk_MVA; rtol = 1e-10)
        else
          @test isequal(a.ik_kA, b.ik_kA)
        end
      end
    end
    # the takahashi sweep is deterministic: repeated runs are bitwise equal
    tak1 = runShortCircuit!(net; case = :max, sweep_method = :takahashi, parallel_enabled = false)
    tak2 = runShortCircuit!(net; case = :max, sweep_method = :takahashi, parallel_enabled = false)
    @test isequal(tak1.rows, tak2.rows)

    # parallel sweep equality with buses = :all on the shipped operated
    # grid (load_fixture_net: runs on every install; the case declares
    # its own IEC 60909 source, so no hand-attached external grid)
    net60 = load_fixture_net("sp_case60")
    s60 = runShortCircuit!(net60; buses = :all, case = :max, parallel_enabled = false)
    p60 = runShortCircuit!(net60; buses = :all, case = :max, parallel_enabled = true, parallel_min_work_items = 2)
    @test isequal(s60.rows, p60.rows)
    parallel_comparison_report("sc all-bus sweep, sp_case60")
  end)() end
end
