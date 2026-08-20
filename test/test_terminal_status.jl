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

# file: test/test_terminal_status.jl
# purpose: per-terminal branch status (r0.9.10, one-sided open branches):
#          the dangling-node reference anchor for lines and transformers,
#          bitwise equivalence of both-open with status 0, reuse-backend pattern
#          invalidation on terminal toggles, isolation of the open-side
#          bus, MATPOWER export shunt substitution, and the result surface
#          (flows, open-end voltage, table markers).

# Two-bus fixture: slack at A, the branch under test A -> B, optional load
# at B. The reference formulation keeps the FULL branch and adds a
# zero-injection auxiliary PQ bus at the open end (dangling-node
# formulation); for a pure pi branch both are exactly equivalent.
function _ts_line_net(name::String; to_status::Int = 1, from_status::Int = 1, load_at_b::Bool = false)
  net = Net(name = name, baseMVA = 100.0)
  addBus!(net = net, busName = "A", vn_kV = 380.0)
  addBus!(net = net, busName = "B", vn_kV = 380.0)
  addProsumer!(net = net, busName = "A", type = "EXTERNALNETWORKINJECTION", referencePri = "A", vm_pu = 1.02, va_deg = 0.0)
  load_at_b && addProsumer!(net = net, busName = "B", type = "ENERGYCONSUMER", p = 20.0, q = 5.0)
  addPIModelACLine!(net = net, fromBus = "A", toBus = "B", r_pu = 0.02, x_pu = 0.16, b_pu = 0.9, g_pu = 0.004, status = 1, from_status = from_status, to_status = to_status)
  ok, msg = validate!(net = net)
  ok || error("fixture validation failed: $msg")
  return net
end

function run_terminal_status_tests()
  @testset "per-terminal branch status (r0.9.10)" begin
    @testset "state helper and setters" begin
      net = _ts_line_net("ts_helper")
      br = net.branchVec[1]
      @test Sparlectra._branch_terminal_state(br) == :closed
      setBranchTerminalStatus!(br; to = false)
      @test Sparlectra._branch_terminal_state(br) == :open_to
      @test br.status == 0
      setBranchTerminalStatus!(br; to = true)
      @test Sparlectra._branch_terminal_state(br) == :closed
      @test br.status == 1
      setBranchTerminalStatus!(br; from = false, to = false)
      @test Sparlectra._branch_terminal_state(br) == :open
      setBranchStatus!(br, true)
      @test Sparlectra._branch_terminal_state(br) == :closed
      # legacy direct aggregate write without touching the flags reads :open
      br.status = 0
      @test Sparlectra._branch_terminal_state(br) == :open
      setBranchStatus!(br, true)
    end

    @testset "dangling-node reference anchor: line, open at to" begin
      # reduction under test
      net = _ts_line_net("ts_anchor"; to_status = 0)
      ite, erg = runpf!(net, 30, 1e-12, 0)
      @test erg == 0
      calcNetLosses!(net)
      br = net.branchVec[1]
      # reference: full branch plus zero-injection auxiliary PQ bus at the
      # open end (the solver computes U_open itself)
      ref = Net(name = "ts_anchor_ref", baseMVA = 100.0)
      addBus!(net = ref, busName = "A", vn_kV = 380.0)
      addBus!(net = ref, busName = "AUX", vn_kV = 380.0)
      addProsumer!(net = ref, busName = "A", type = "EXTERNALNETWORKINJECTION", referencePri = "A", vm_pu = 1.02, va_deg = 0.0)
      addProsumer!(net = ref, busName = "AUX", type = "ENERGYCONSUMER", p = 0.0, q = 0.0)
      addPIModelACLine!(net = ref, fromBus = "A", toBus = "AUX", r_pu = 0.02, x_pu = 0.16, b_pu = 0.9, g_pu = 0.004, status = 1)
      ok, _ = validate!(net = ref)
      @test ok
      ite_r, erg_r = runpf!(ref, 30, 1e-12, 0)
      @test erg_r == 0
      calcNetLosses!(ref)
      rbr = ref.branchVec[1]
      # correctness anchor: U_from, S_from, and U_open agree to 1e-10
      @test isapprox(get_bus_vm_pu(net, "A"), get_bus_vm_pu(ref, "A"); atol = 1e-10)
      @test isapprox(br.fBranchFlow.pFlow, rbr.fBranchFlow.pFlow; atol = 1e-10)
      @test isapprox(br.fBranchFlow.qFlow, rbr.fBranchFlow.qFlow; atol = 1e-10)
      @test isapprox(br.open_end_vm_pu, get_bus_vm_pu(ref, "AUX"); atol = 1e-10)
      @test isapprox(br.open_end_va_deg, ref.nodeVec[ref.busDict["AUX"]]._va_deg; atol = 1e-10)
      # Ferranti rise and FULL charging (not half): |Q_from| close to b*V^2
      @test br.open_end_vm_pu > get_bus_vm_pu(net, "A")
      @test abs(br.fBranchFlow.qFlow) > 0.75 * 0.9 * 100.0
      # the open terminal carries zero by definition
      @test br.tBranchFlow.pFlow == 0.0 && br.tBranchFlow.qFlow == 0.0
      # loss bookkeeping: branch loss equals the closed-terminal power
      @test isapprox(br.pLosses, br.fBranchFlow.pFlow; atol = 1e-12)
    end

    @testset "dangling-node reference anchor: transformer, both sides" begin
      for (fs, ts) in ((1, 0), (0, 1))
        net = Net(name = "ts_trafo_$(fs)$(ts)", baseMVA = 100.0)
        addBus!(net = net, busName = "H", vn_kV = 380.0)
        addBus!(net = net, busName = "L", vn_kV = 110.0)
        # keep the reference on whichever side stays closed
        refbus = fs == 1 ? "H" : "L"
        addProsumer!(net = net, busName = refbus, type = "EXTERNALNETWORKINJECTION", referencePri = refbus, vm_pu = 1.0, va_deg = 0.0)
        addPIModelTrafo!(net = net, fromBus = "H", toBus = "L", r_pu = 0.004, x_pu = 0.1, b_pu = 0.06, status = 1, ratio = 1.05, shift_deg = 5.0, from_status = fs, to_status = ts)
        ok, _ = validate!(net = net)
        @test ok
        ite, erg = runpf!(net, 30, 1e-12, 0)
        @test erg == 0
        calcNetLosses!(net)
        br = net.branchVec[1]
        st = Sparlectra._branch_terminal_state(br)
        @test st == (ts == 0 ? :open_to : :open_from)
        # reference with an auxiliary zero-injection PQ bus on the open side
        ref = Net(name = "ts_trafo_ref_$(fs)$(ts)", baseMVA = 100.0)
        addBus!(net = ref, busName = "H", vn_kV = 380.0)
        addBus!(net = ref, busName = "L", vn_kV = 110.0)
        addProsumer!(net = ref, busName = refbus, type = "EXTERNALNETWORKINJECTION", referencePri = refbus, vm_pu = 1.0, va_deg = 0.0)
        openbus = fs == 1 ? "L" : "H"
        addProsumer!(net = ref, busName = openbus, type = "ENERGYCONSUMER", p = 0.0, q = 0.0)
        addPIModelTrafo!(net = ref, fromBus = "H", toBus = "L", r_pu = 0.004, x_pu = 0.1, b_pu = 0.06, status = 1, ratio = 1.05, shift_deg = 5.0)
        ok2, _ = validate!(net = ref)
        @test ok2
        _, erg_r = runpf!(ref, 30, 1e-12, 0)
        @test erg_r == 0
        calcNetLosses!(ref)
        u_open_ref = get_bus_vm_pu(ref, openbus)
        @test isapprox(br.open_end_vm_pu, u_open_ref; atol = 1e-10)
        closed_flow = ts == 0 ? br.fBranchFlow : br.tBranchFlow
        ref_flow = ts == 0 ? ref.branchVec[1].fBranchFlow : ref.branchVec[1].tBranchFlow
        @test isapprox(closed_flow.pFlow, ref_flow.pFlow; atol = 1e-10)
        @test isapprox(closed_flow.qFlow, ref_flow.qFlow; atol = 1e-10)
      end
    end

    @testset "both flags open equals status 0 bit for bit" begin
      net_a = _ts_line_net("ts_bits_a"; load_at_b = true)
      net_b = _ts_line_net("ts_bits_b"; load_at_b = true)
      setBranchStatus!(net_a.branchVec[1], false)
      setBranchTerminalStatus!(net_b.branchVec[1]; from = false, to = false)
      markIsolatedBuses!(net = net_a, log = false)
      markIsolatedBuses!(net = net_b, log = false)
      Ya = createYBUS(net = net_a)
      Yb = createYBUS(net = net_b)
      @test Ya == Yb
    end

    @testset "terminal toggle between two solves (reuse-backend pattern invalidation)" begin
      # three buses so the net stays solvable with the branch open: A - B
      # (under test) and A - C carrying a load
      net = Net(name = "ts_toggle", baseMVA = 100.0)
      for (b, v) in (("A", 380.0), ("B", 380.0), ("C", 380.0))
        addBus!(net = net, busName = b, vn_kV = v)
      end
      addProsumer!(net = net, busName = "A", type = "EXTERNALNETWORKINJECTION", referencePri = "A", vm_pu = 1.0, va_deg = 0.0)
      addProsumer!(net = net, busName = "B", type = "ENERGYCONSUMER", p = 10.0, q = 2.0)
      addProsumer!(net = net, busName = "C", type = "ENERGYCONSUMER", p = 15.0, q = 3.0)
      addPIModelACLine!(net = net, fromBus = "A", toBus = "B", r_pu = 0.01, x_pu = 0.08, b_pu = 0.2, status = 1)
      addPIModelACLine!(net = net, fromBus = "A", toBus = "C", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
      ok, _ = validate!(net = net)
      @test ok
      _, erg1 = runpf!(net, 30, 1e-10, 0)
      @test erg1 == 0
      vm_b_closed = get_bus_vm_pu(net, "B")
      # open the to terminal of A-B: B loses its supply, becomes isolated
      setBranchTerminalStatus!(net.branchVec[1]; to = false)
      markIsolatedBuses!(net = net, log = false)
      _, erg2 = runpf!(net, 30, 1e-10, 0)
      @test erg2 == 0
      calcNetLosses!(net)
      br = net.branchVec[1]
      @test br.open_end_vm_pu !== nothing
      # charging draw at A present, no through flow to B
      @test br.tBranchFlow.pFlow == 0.0
      # reclose and solve again: pattern changes back, result matches run 1
      setBranchTerminalStatus!(net.branchVec[1]; to = true)
      markIsolatedBuses!(net = net, log = false)
      refreshBusTypesFromProsumers!(net)
      _, erg3 = runpf!(net, 30, 1e-10, 0)
      @test erg3 == 0
      @test isapprox(get_bus_vm_pu(net, "B"), vm_b_closed; atol = 1e-9)
    end

    @testset "island report: open-side bus is isolated" begin
      net = _ts_line_net("ts_iso"; to_status = 0, load_at_b = false)
      @test net.busDict["B"] in net.isoNodes
      report = Sparlectra.detect_ac_islands(net)
      @test length(report.rows) == 1
      @test report.rows[1].n_bus == 1
    end

    @testset "MATPOWER export: partial branch as status 0 plus exact Yin shunt" begin
      # three buses: A stays energized through A-C after the reimport turns
      # the partial branch fully off (its Yin travels as a bus shunt at A)
      net = Net(name = "ts_mpexport", baseMVA = 100.0)
      for b in ("A", "B", "C")
        addBus!(net = net, busName = b, vn_kV = 380.0)
      end
      addProsumer!(net = net, busName = "A", type = "EXTERNALNETWORKINJECTION", referencePri = "A", vm_pu = 1.02, va_deg = 0.0)
      addProsumer!(net = net, busName = "C", type = "ENERGYCONSUMER", p = 25.0, q = 6.0)
      addPIModelACLine!(net = net, fromBus = "A", toBus = "B", r_pu = 0.02, x_pu = 0.16, b_pu = 0.9, g_pu = 0.004, status = 1, to_status = 0)
      addPIModelACLine!(net = net, fromBus = "A", toBus = "C", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
      ok, _ = validate!(net = net)
      @test ok
      _, erg = runpf!(net, 30, 1e-12, 0)
      @test erg == 0
      calcNetLosses!(net)
      vm_a = get_bus_vm_pu(net, "A")
      vm_c = get_bus_vm_pu(net, "C")
      dir = mktempdir()
      path = joinpath(dir, "ts_mpexport.m")
      writeMatpowerCasefile(net, path)
      txt = read(path, String)
      @test occursin("open_terminal=to", txt)
      net2 = createNetFromMatPowerFile(filename = path)
      _, erg2 = runpf!(net2, 30, 1e-12, 0)
      @test erg2 == 0
      # same voltages as the native net (shunt substitution correct); the
      # partial state itself is lost by design (branch back as out of
      # service plus a shunt)
      vms2 = sort([n._vm_pu for n in net2.nodeVec if !(n.busIdx in net2.isoNodes)])
      vms1 = sort([vm_a, vm_c])
      @test isapprox(vms2, vms1; atol = 1e-8)
      @test all(br -> Sparlectra._branch_terminal_state(br) in (:closed, :open), net2.branchVec)
    end

    @testset "CGMES roundtrip preserves the terminal flags" begin
      # a net with one partially open line and one partially open 2WT plus
      # enough closed topology to stay importable
      net = Net(name = "ts_cgmes", baseMVA = 100.0)
      for (b, v) in (("A", 380.0), ("B", 380.0), ("C", 380.0), ("L", 110.0))
        addBus!(net = net, busName = b, vn_kV = v)
      end
      addProsumer!(net = net, busName = "A", type = "EXTERNALNETWORKINJECTION", referencePri = "A", vm_pu = 1.0, va_deg = 0.0)
      addProsumer!(net = net, busName = "C", type = "ENERGYCONSUMER", p = 15.0, q = 4.0)
      addPIModelACLine!(net = net, fromBus = "A", toBus = "C", r_pu = 0.01, x_pu = 0.08, b_pu = 0.1, status = 1)
      addPIModelACLine!(net = net, fromBus = "A", toBus = "B", r_pu = 0.02, x_pu = 0.16, b_pu = 0.6, status = 1, to_status = 0)
      addPIModelTrafo!(net = net, fromBus = "C", toBus = "L", r_pu = 0.004, x_pu = 0.1, b_pu = 0.02, status = 1, ratio = 1.0, to_status = 0)
      ok, _ = validate!(net = net)
      @test ok
      _, erg = runpf!(net, 30, 1e-10, 0)
      @test erg == 0
      calcNetLosses!(net)
      dir = mktempdir()
      files = writeCGMESFiles(net; path = dir)
      res = importCGMES(path = dir, name = "ts_cgmes_rt")
      states = sort([Sparlectra._branch_terminal_state(br) for br in res.net.branchVec])
      @test count(==(:closed), states) == 1
      @test count(s -> s in (:open_to, :open_from), states) == 2
      @test any(m -> occursin("one open terminal", m), res.messages)
    end

    @testset "result surface: table marker, header count, report columns" begin
      net = _ts_line_net("ts_results"; to_status = 0)
      _, erg = runpf!(net, 30, 1e-10, 0)
      @test erg == 0
      calcNetLosses!(net)
      mktempdir() do d
        cd(d) do
          printACPFlowResults(net, 0.0, 1, 1e-10, true)
          txt = read("result_$(net.name).txt", String)
          @test occursin("open@to", txt)
          @test occursin("Open terminals :         1", txt)
        end
      end
      rep = buildACPFlowReport(net; ite = 1, converged = true)
      row = only(rep.branches)
      @test row.terminal_state == :open_to
      @test !ismissing(row.open_end_vm_pu) && row.open_end_vm_pu > 1.0
      # SE guard: a flow measurement on the partial branch is rejected
      meas = [Sparlectra.Measurement(typ = Sparlectra.PflowMeas, value = 0.0, sigma = 0.01, branchIdx = 1, direction = :from, id = "bad")]
      cfg = StateEstimationConfig(max_iter = 10, tol = 1e-6)
      @test_throws ErrorException runse!(net, meas, cfg)
    end
  end
  return true
end
