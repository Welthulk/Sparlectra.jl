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
# purpose: short-circuit verification — hand-derived IEC 60909-0 reference cases
#          (analytic derivations in the comments, computed independently of
#          the implementation), the safety-flag contract (default + warning
#          + result flag), asynchronous-machine lower-bound flagging, and
#          short_circuit.c_factor config coverage. The MicroGrid
#          plausibility sweep lives in test_cgmes_importer.jl (extended,
#          cache-gated); this file needs no CGMES fixtures.

using Test
using Sparlectra

const _SC = Sparlectra.CGMESImporter.CGMESShortCircuitData

_sc_feeder(; bus, ikmax_A = nothing, ikmin_A = nothing, rx_max = nothing, rx_min = nothing, name = "Feeder") =
  (mrid = name, name = name, bus = bus, maxInitialSymShCCurrent_A = ikmax_A, minInitialSymShCCurrent_A = ikmin_A, maxR1ToX1Ratio = rx_max, minR1ToX1Ratio = rx_min, maxR0ToX0Ratio = nothing, maxZ0ToZ1Ratio = nothing, ikSecond = nothing, governorSCD = nothing)

_sc_machine(; bus, xdpp = nothing, ratedS = nothing, ratedU = nothing, name = "G1") =
  (mrid = name, name = name, bus = bus, satDirectSubtransX_pu = xdpp, satDirectTransX_pu = nothing, r0_pu = nothing, x0_pu = nothing, r2_pu = nothing, x2_pu = nothing, earthing = nothing, ratedS_MVA = ratedS, ratedU_kV = ratedU)

_sc_asm(; bus, ilr = nothing, rx = nothing, ratedS = nothing, ratedU = nothing, p_mech = nothing, efficiency = nothing, cosphi = nothing, pole_pairs = nothing, name = "M1") =
  (mrid = name, name = name, bus = bus, iaIrRatio = ilr, rxLockedRotorRatio = rx, efficiency_percent = efficiency, ratedMechanicalPower_MW = p_mech, polePairNumber = pole_pairs, ratedS_MVA = ratedS, ratedU_kV = ratedU, ratedPowerFactor = cosphi)

_sc_data(; enis = NamedTuple[], machines = NamedTuple[], asms = NamedTuple[]) = _SC(enis, machines, NamedTuple[], NamedTuple[], NamedTuple[], asms)

# Single 110 kV bus fed by one network feeder.
function _sc_single_bus_net()
  net = Net(name = "sc_single", baseMVA = 100.0)
  addBus!(net = net, busName = "B1", vn_kV = 110.0)
  return net
end

# 110 kV feeder bus with one line to a second 110 kV bus.
function _sc_two_bus_net(; r_pu = 0.05, x_pu = 0.5)
  net = Net(name = "sc_two_bus", baseMVA = 100.0)
  addBus!(net = net, busName = "B1", vn_kV = 110.0)
  addBus!(net = net, busName = "B2", vn_kV = 110.0)
  # b_pu deliberately non-zero: IEC drops line charging, so it must not
  # change the result (asserted below).
  addPIModelACLine!(net = net, fromBus = "B1", toBus = "B2", r_pu = r_pu, x_pu = x_pu, b_pu = 0.2, status = 1)
  return net
end

function run_short_circuit_tests()
  @testset "Short circuit (IEC 60909-0)" begin
    @testset "feeder reproduces its declared current at the connection point" begin
      # Hand derivation: Z_Q = c·Un/(√3·Ik_declared); at the feeder bus the
      # whole short-circuit impedance IS Z_Q, so Ik'' = c·Un/(√3·Z_Q)
      # = Ik_declared exactly — the c factor cancels. This is the IEC
      # semantics of the declared feeder current and independent of every
      # per-unit conversion in the implementation.
      net = _sc_single_bus_net()
      sc = _sc_data(enis = [_sc_feeder(bus = "B1", ikmax_A = 10_000.0, ikmin_A = 8_000.0, rx_max = 0.1, rx_min = 0.1)])
      rmax = runShortCircuit!(net, sc; case = :max)
      @test length(rmax.rows) == 1
      row = rmax.rows[1]
      @test row.status === :ok
      @test row.contains_defaulted_data == false
      @test isapprox(row.ik_kA, 10.0; rtol = 1e-9)
      # Sk'' = √3·Un·Ik'' = √3·110·10 = 1905.2558... MVA (hand value)
      @test isapprox(row.sk_MVA, sqrt(3.0) * 110.0 * 10.0; rtol = 1e-9)
      # R/X = 0.1 → κ_b = 1.02 + 0.98·e^(−0.3) = 1.746002…;
      # 1.15·κ_b = 2.0079… → HV cap 2.0 applies; i_p = 2.0·√2·10 = 28.284… kA
      @test isapprox(row.rx_ratio, 0.1; rtol = 1e-9)
      @test isapprox(row.kappa, 2.0; rtol = 1e-12)
      @test isapprox(row.ip_kA, 2.0 * sqrt(2.0) * 10.0; rtol = 1e-9)
      rmin = runShortCircuit!(net, sc; case = :min)
      @test isapprox(rmin.rows[1].ik_kA, 8.0; rtol = 1e-9)
      @test rmin.rows[1].c == 1.00   # IEC Table 1, HV c_min
      @test rmax.rows[1].c == 1.10   # IEC Table 1, HV c_max
    end

    @testset "line impedance adds per hand-derived series model" begin
      # Hand derivation on 100 MVA / 110 kV (Z_base = 121 Ω):
      #   Z_Q  = 1.10·110000/(√3·10000) Ω, split X_Q = Z_Q/√(1+0.1²),
      #          R_Q = 0.1·X_Q  (independent arithmetic below)
      #   Z_L  = (0.05 + j0.5)·121 Ω
      #   Zk(B2) = Z_Q + Z_L, Ik'' = 1.10·110/(√3·|Zk|) kA.
      # Line charging (b_pu = 0.2) must NOT enter — IEC drops shunts.
      net = _sc_two_bus_net()
      sc = _sc_data(enis = [_sc_feeder(bus = "B1", ikmax_A = 10_000.0, rx_max = 0.1)])
      r = runShortCircuit!(net, sc; case = :max)
      zq = 1.10 * 110_000.0 / (sqrt(3.0) * 10_000.0)
      xq = zq / sqrt(1.0 + 0.1^2)
      rq = 0.1 * xq
      zl = (0.05 + 0.5im) * 121.0
      zk = (rq + im * xq) + zl
      ik_expected = 1.10 * 110.0 / (sqrt(3.0) * abs(zk))
      row2 = only(row for row in r.rows if row.bus == "B2")
      @test row2.status === :ok
      @test isapprox(row2.zk_ohm, abs(zk); rtol = 1e-9)
      @test isapprox(row2.ik_kA, ik_expected; rtol = 1e-9)
      @test isapprox(row2.rx_ratio, real(zk) / imag(zk); rtol = 1e-9)
      # buses list variant selects a single bus
      rb = runShortCircuit!(net, sc; buses = ["B2"], case = :max)
      @test length(rb.rows) == 1 && rb.rows[1].bus == "B2"
      @test isapprox(rb.rows[1].ik_kA, ik_expected; rtol = 1e-12)
      @test_throws ArgumentError runShortCircuit!(net, sc; buses = ["NOPE"], case = :max)
      @test_throws ArgumentError runShortCircuit!(net, sc; case = :typical)
    end

    @testset "synchronous machine per hand-derived x''d conversion" begin
      # Hand derivation (10.5 kV machine bus, Sbase 100 MVA):
      #   x''d = 0.15 pu on 50 MVA / 10.5 kV → x_net = 0.15·(100/50)·(10.5/10.5)² = 0.3 pu
      #   R_Gf = 0.07·x  (HV machine < 100 MVA, IEC 60909-0 §6.6.3)
      #   Z_base = 10.5²/100 Ω; Ik'' = 1.10·10.5/(√3·|Z|·Z_base) — independent below.
      net = Net(name = "sc_machine", baseMVA = 100.0)
      addBus!(net = net, busName = "M1", vn_kV = 10.5)
      sc = _sc_data(machines = [_sc_machine(bus = "M1", xdpp = 0.15, ratedS = 50.0, ratedU = 10.5)])
      r = runShortCircuit!(net, sc; case = :max)
      x_pu = 0.15 * (100.0 / 50.0)
      z_pu = 0.07 * x_pu + im * x_pu
      z_ohm = abs(z_pu) * (10.5^2 / 100.0)
      ik_expected = 1.10 * 10.5 / (sqrt(3.0) * z_ohm)
      row = r.rows[1]
      @test row.status === :ok
      @test row.contains_defaulted_data == false
      @test isapprox(row.ik_kA, ik_expected; rtol = 1e-9)
      # scalar c-factor override replaces the table linearly
      r2 = runShortCircuit!(net, sc; case = :max, c_factor = 1.0)
      @test isapprox(r2.rows[1].ik_kA, ik_expected / 1.10; rtol = 1e-9)
      @test_throws ArgumentError runShortCircuit!(net, sc; c_factor = 1.4)
    end

    @testset "safety-flag contract: default + warning + result flag" begin
      # A machine without x''d triggers all three: the documented default
      # (0.2 pu machine base), a log warning, and the per-row flag.
      net = Net(name = "sc_flag", baseMVA = 100.0)
      addBus!(net = net, busName = "M1", vn_kV = 10.5)
      sc = _sc_data(machines = [_sc_machine(bus = "M1", xdpp = nothing, ratedS = 50.0, ratedU = 10.5)])
      r = @test_logs (:warn, r"no usable x''_d — default 0.2 pu") match_mode = :any runShortCircuit!(net, sc; case = :max)
      row = r.rows[1]
      @test row.status === :ok
      @test row.contains_defaulted_data == true
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
    end

    @testset "asynchronous machine per hand-derived IEC §6.7 motor impedance" begin
      # Hand derivation (10 kV motor bus, Sbase 100 MVA):
      #   Z_M = (1/ilr)·U_rM²/S_rM = (1/5)·10²/5 = 4 Ω, R/X = 0.1
      #   → X = 4/√1.01, R = 0.1·X; Ik'' = 1.10·10/(√3·|Z_M|) at the bus.
      net = Net(name = "sc_motor", baseMVA = 100.0)
      addBus!(net = net, busName = "M1", vn_kV = 10.0)
      sc = _sc_data(asms = [_sc_asm(bus = "M1", ilr = 5.0, rx = 0.1, ratedS = 5.0, ratedU = 10.0)])
      r = runShortCircuit!(net, sc; case = :max)
      row = r.rows[1]
      @test row.status === :ok
      @test row.contains_defaulted_data == false
      @test isapprox(row.zk_ohm, 4.0; rtol = 1e-9)
      @test isapprox(row.ik_kA, 1.10 * 10.0 / (sqrt(3.0) * 4.0); rtol = 1e-9)
      @test isapprox(row.rx_ratio, 0.1; rtol = 1e-9)

      # S_rM via mechanical power: P_rM/(η·cosφ) = 3.6/(0.9·0.8) = 5 MVA —
      # must reproduce the ratedS path exactly.
      sc_p = _sc_data(asms = [_sc_asm(bus = "M1", ilr = 5.0, rx = 0.1, ratedU = 10.0, p_mech = 3.6, efficiency = 90.0, cosphi = 0.8)])
      rp = runShortCircuit!(net, sc_p; case = :max)
      @test isapprox(rp.rows[1].ik_kA, row.ik_kA; rtol = 1e-12)
      @test rp.rows[1].contains_defaulted_data == false

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
    end

    @testset "asynchronous machine: skipped with island-wide lower-bound flag (max only)" begin
      net = _sc_two_bus_net()
      addProsumer!(net = net, busName = "B2", type = "ASYNCHRONOUSMACHINE", p = 5.0, q = 2.0)
      # same component-type marker the CGMES importer sets for motors (the
      # prosumer constructor collapses the type string to a generic Load)
      net.prosumpsVec[end].comp.cTyp = Sparlectra.AsynchronousMachine
      sc = _sc_data(enis = [_sc_feeder(bus = "B1", ikmax_A = 10_000.0, ikmin_A = 8_000.0, rx_max = 0.1, rx_min = 0.1)])
      rmax = @test_logs (:warn, r"AsynchronousMachine .* lower bound") match_mode = :any runShortCircuit!(net, sc; case = :max)
      # the whole island is a lower bound, not just the motor bus
      @test all(row.contains_defaulted_data for row in rmax.rows)
      @test all(any(occursin("lower bound", reason) for reason in row.reasons) for row in rmax.rows)
      # motors never contribute to (and never flag) the minimum case
      rmin = runShortCircuit!(net, sc; case = :min)
      @test all(!row.contains_defaulted_data for row in rmin.rows)
    end

    @testset "island without a source reports :no_source, isolated buses :isolated" begin
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
    end

    @testset "WebUI data-check gating helper" begin
      # Byte scan over the delivery contents (fast fixture: plain XML folder;
      # the real ZIP path goes through the same collectCGMESFiles reader and
      # is exercised by the extended CGMES service test).
      with_data = mktempdir()
      write(joinpath(with_data, "grid_EQ.xml"), "<rdf:RDF><cim:SynchronousMachine rdf:ID=\"g\"/></rdf:RDF>")
      without_data = mktempdir()
      write(joinpath(without_data, "grid_EQ.xml"), "<rdf:RDF><cim:EnergyConsumer rdf:ID=\"l\"/></rdf:RDF>")
      @test Sparlectra._webui_case_has_short_circuit_data(with_data) == true
      @test Sparlectra._webui_case_has_short_circuit_data(without_data) == false
      # cached second call returns the same verdict
      @test Sparlectra._webui_case_has_short_circuit_data(without_data) == false
      # unresolvable paths must NOT lock the button (service explains instead)
      @test Sparlectra._webui_case_has_short_circuit_data(joinpath(without_data, "missing.zip")) == true
    end

    @testset "short_circuit.c_factor config coverage" begin
      @test Sparlectra.ShortCircuitConfig().c_factor == 0.0
      ok_yaml = tempname() * ".yaml"
      write(ok_yaml, "short_circuit:\n  c_factor: 1.05\n")
      @test Sparlectra.load_sparlectra_config(ok_yaml; reload = true).shortcircuit.c_factor == 1.05
      bad_yaml = tempname() * ".yaml"
      write(bad_yaml, "short_circuit:\n  c_factor: 1.4\n")
      err = try
        Sparlectra.load_sparlectra_config(bad_yaml; reload = true)
        nothing
      catch e
        e
      end
      @test err isa ArgumentError
      @test occursin("short_circuit.c_factor", sprint(showerror, err))
    end

    @testset "parallel all-bus sweep identity (Phase 3)" begin
      # multi-island SC fixture: n feeder-fed rings with declared Sk''
      function build_sweep_net(n, m)
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

      net = build_sweep_net(3, 40)
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
      if Threads.nthreads() > 1
        println("sc parallel sweep: RAN with ", Threads.nthreads(), " threads")
      else
        println("sc parallel sweep: fallback-only run (single-threaded test process); the threaded path runs in the --threads=4 battery")
      end

      # sweep_method = :takahashi: one selected-inverse pass per island;
      # must agree with :solves to machine precision (documented: not
      # bitwise) with identical statuses/flags, serial and parallel alike
      @test_throws ArgumentError runShortCircuit!(net; sweep_method = :bogus)
      solves_ref = runShortCircuit!(net; case = :max, parallel_enabled = false)
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

      # case14 with buses = :all through the MATPOWER import (bundled case)
      case14 = joinpath(Sparlectra.MPOWER_DIR, "case14.m")
      if isfile(case14)
        net14 = createNetFromMatPowerFile(filename = case14)
        first_bus = argmin(name -> net14.busDict[name], collect(keys(net14.busDict)))
        addExternalGrid!(net = net14, busName = first_bus, vm_pu = 1.06, sk_max_MVA = 3000.0, sk_min_MVA = 2000.0, rx_max = 0.1, internal_impedance = false)
        s14 = runShortCircuit!(net14; buses = :all, case = :max, parallel_enabled = false)
        p14 = runShortCircuit!(net14; buses = :all, case = :max, parallel_enabled = true, parallel_min_work_items = 2)
        @test isequal(s14.rows, p14.rows)
        println("sc parallel sweep case14: RAN")
      else
        println("sc parallel sweep case14: SKIPPED (bundled case14.m not found)")
      end
    end
  end
end
