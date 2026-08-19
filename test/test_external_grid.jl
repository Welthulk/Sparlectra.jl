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
# file: test/test_external_grid.jl
# purpose: external-grid element (issue #299) — addExternalGrid! as ideal
#          slack (PF invariance of carried SC data), the native feeder record
#          against hand-derived IEC 60909-0 values, min-case and R/X flag
#          semantics, input validation, net-copy regression, and the
#          internal-impedance (non-ideal source) variant incl. the stiff
#          sk->inf limit against the ideal representation.

using Test
using Sparlectra

# Four-bus meshed 110 kV fixture with one generator and two loads: small by
# test convention, but structurally a "real" net (mesh, PV bus, mixed loads)
# so the PF-invariance assertion is not trivially empty.
function _eg_fixture(; slack_bus::String = "B1", vm_slack::Float64 = 1.02)
  net = Net(name = "extgrid_fixture", baseMVA = 100.0)
  for b in ("B1", "B2", "B3", "B4")
    addBus!(net = net, busName = b, vn_kV = 110.0)
  end
  addPIModelACLine!(net = net, fromBus = "B1", toBus = "B2", r_pu = 0.01, x_pu = 0.06, b_pu = 0.02, status = 1)
  addPIModelACLine!(net = net, fromBus = "B2", toBus = "B3", r_pu = 0.02, x_pu = 0.09, b_pu = 0.02, status = 1)
  addPIModelACLine!(net = net, fromBus = "B3", toBus = "B4", r_pu = 0.015, x_pu = 0.08, b_pu = 0.02, status = 1)
  addPIModelACLine!(net = net, fromBus = "B4", toBus = "B1", r_pu = 0.012, x_pu = 0.07, b_pu = 0.02, status = 1)
  addProsumer!(net = net, busName = slack_bus, type = "EXTERNALNETWORKINJECTION", referencePri = slack_bus, vm_pu = vm_slack, va_deg = 0.0)
  addProsumer!(net = net, busName = "B3", type = "GENERATOR", p = 30.0, vm_pu = 1.01, isRegulated = true)
  addProsumer!(net = net, busName = "B2", type = "ENERGYCONSUMER", p = 40.0, q = 12.0)
  addProsumer!(net = net, busName = "B4", type = "ENERGYCONSUMER", p = 35.0, q = 8.0)
  return net
end

_eg_state(net) = [(getNodeVm(n), n._va_deg) for n in net.nodeVec]

function run_external_grid_tests()
  @testset "External grid (issue #299)" begin
    @testset "duck-typing contract: field-identical to the CGMES container" begin
      # The SC engine is duck-typed over the record contract — the two
      # containers must stay field-identical by name, order, and type. This
      # guard enforces machine-checkably what the native_sc_data.jl comment
      # can only ask for.
      @test fieldnames(NativeShortCircuitData) == fieldnames(Sparlectra.CGMESImporter.CGMESShortCircuitData)
      @test fieldtypes(NativeShortCircuitData) == fieldtypes(Sparlectra.CGMESImporter.CGMESShortCircuitData)
    end

    @testset "PF invariance: carried SC data changes no power-flow result" begin
      ref = _eg_fixture()
      ite_ref, erg_ref = runpf!(ref, 30, 1e-8, 0)
      @test erg_ref == 0

      # Same fixture, but the slack is created through addExternalGrid! with
      # finite short-circuit data. The PF side must be exactly the manual
      # slack-prosumer path.
      net = Net(name = "extgrid_fixture", baseMVA = 100.0)
      for b in ("B1", "B2", "B3", "B4")
        addBus!(net = net, busName = b, vn_kV = 110.0)
      end
      addPIModelACLine!(net = net, fromBus = "B1", toBus = "B2", r_pu = 0.01, x_pu = 0.06, b_pu = 0.02, status = 1)
      addPIModelACLine!(net = net, fromBus = "B2", toBus = "B3", r_pu = 0.02, x_pu = 0.09, b_pu = 0.02, status = 1)
      addPIModelACLine!(net = net, fromBus = "B3", toBus = "B4", r_pu = 0.015, x_pu = 0.08, b_pu = 0.02, status = 1)
      addPIModelACLine!(net = net, fromBus = "B4", toBus = "B1", r_pu = 0.012, x_pu = 0.07, b_pu = 0.02, status = 1)
      addExternalGrid!(net = net, busName = "B1", vm_pu = 1.02, sk_max_MVA = 3000.0, rx_max = 0.1)
      addProsumer!(net = net, busName = "B3", type = "GENERATOR", p = 30.0, vm_pu = 1.01, isRegulated = true)
      addProsumer!(net = net, busName = "B2", type = "ENERGYCONSUMER", p = 40.0, q = 12.0)
      addProsumer!(net = net, busName = "B4", type = "ENERGYCONSUMER", p = 35.0, q = 8.0)

      @test net.slackVec == ref.slackVec
      ite, erg = runpf!(net, 30, 1e-8, 0)
      @test erg == 0
      @test ite == ite_ref
      for (a, b) in zip(_eg_state(net), _eg_state(ref))
        @test a[1] ≈ b[1] atol = 1e-12
        @test a[2] ≈ b[2] atol = 1e-12
      end

      # Conversion into the feeder record: full ENI contract, at add time.
      @test length(net.sc_sources.external_network_injections) == 1
      eni = net.sc_sources.external_network_injections[1]
      @test eni.bus == "B1"
      @test eni.mrid == "native-eni-B1"
      @test eni.maxInitialSymShCCurrent_A ≈ 3000.0e6 / (sqrt(3.0) * 110.0e3) rtol = 1e-12
      @test eni.minInitialSymShCCurrent_A === nothing
      @test eni.maxR1ToX1Ratio == 0.1
      @test eni.minR1ToX1Ratio === nothing
      @test eni.maxR0ToX0Ratio === nothing && eni.governorSCD === nothing
    end

    @testset "SC hand calculation (single 110 kV bus, c = 1.1)" begin
      # Hand-derived independently of the engine: Zq = c·Un²/Sk,
      # Ik'' = c·Un/(√3·|Zq|) = Sk/(√3·Un) (the c cancels by construction),
      # Sk'' = √3·Un·Ik'', κ = min(1.15·(1.02 + 0.98·e^(−3·R/X)), 2.0),
      # ip = κ·√2·Ik''.
      for rx in (0.1, 0.3)
        net = Net(name = "sc_hand", baseMVA = 100.0)
        addBus!(net = net, busName = "B1", vn_kV = 110.0)
        addExternalGrid!(net = net, busName = "B1", sk_max_MVA = 3000.0, rx_max = rx)
        res = runShortCircuit!(net; case = :max, c_factor = 1.1)
        row = only(res.rows)
        ik_hand = 3000.0 / (sqrt(3.0) * 110.0)
        kappa_hand = min(1.15 * (1.02 + 0.98 * exp(-3.0 * rx)), 2.0)
        @test row.status === :ok
        @test row.contains_defaulted_data == false
        @test row.zk_ohm ≈ 1.1 * 110.0^2 / 3000.0 rtol = 1e-6
        @test row.rx_ratio ≈ rx rtol = 1e-6
        @test row.ik_kA ≈ ik_hand rtol = 1e-6
        @test row.sk_MVA ≈ 3000.0 rtol = 1e-6
        @test row.kappa ≈ kappa_hand rtol = 1e-6
        @test row.ip_kA ≈ kappa_hand * sqrt(2.0) * ik_hand rtol = 1e-6
      end
    end

    @testset "min case and R/X flag semantics" begin
      # With sk_min: the :min case uses it; rx_min defaults to rx_max
      # (task decision) so the deliberately declared minimum is NOT flagged.
      net = Net(name = "sc_min", baseMVA = 100.0)
      addBus!(net = net, busName = "B1", vn_kV = 110.0)
      addExternalGrid!(net = net, busName = "B1", sk_max_MVA = 3000.0, sk_min_MVA = 1500.0, rx_max = 0.2)
      eni = only(net.sc_sources.external_network_injections)
      @test eni.minR1ToX1Ratio == 0.2
      resmin = runShortCircuit!(net; case = :min)
      rowmin = only(resmin.rows)
      @test rowmin.status === :ok
      @test rowmin.contains_defaulted_data == false
      @test rowmin.ik_kA ≈ 1500.0 / (sqrt(3.0) * 110.0) rtol = 1e-6

      # An explicit rx_min wins over the rx_max default.
      net2 = Net(name = "sc_min2", baseMVA = 100.0)
      addBus!(net = net2, busName = "B1", vn_kV = 110.0)
      addExternalGrid!(net = net2, busName = "B1", sk_max_MVA = 3000.0, sk_min_MVA = 1500.0, rx_max = 0.2, rx_min = 0.35)
      @test only(net2.sc_sources.external_network_injections).minR1ToX1Ratio == 0.35
      resmin2 = runShortCircuit!(net2; case = :min)
      @test only(resmin2.rows).rx_ratio ≈ 0.35 rtol = 1e-6
    end

    @testset "min case without sk_min skips the feeder with the engine flag" begin
      net = Net(name = "sc_nomin", baseMVA = 100.0)
      addBus!(net = net, busName = "B1", vn_kV = 110.0)
      addExternalGrid!(net = net, busName = "B1", sk_max_MVA = 3000.0)
      resmax = runShortCircuit!(net; case = :max)
      @test only(resmax.rows).status === :ok
      resmin = runShortCircuit!(net; case = :min)
      rowmin = only(resmin.rows)
      @test rowmin.status === :no_source
      @test rowmin.contains_defaulted_data == true
      @test any(occursin("minInitialSymShCCurrent", r) for r in rowmin.reasons)
    end

    @testset "parallel feeders stack; mrids stay unique" begin
      net = Net(name = "sc_par", baseMVA = 100.0)
      addBus!(net = net, busName = "B1", vn_kV = 110.0)
      addExternalGrid!(net = net, busName = "B1", sk_max_MVA = 1000.0, rx_max = 0.1)
      addExternalGrid!(net = net, busName = "B1", sk_max_MVA = 1000.0, rx_max = 0.1)
      mrids = [e.mrid for e in net.sc_sources.external_network_injections]
      @test length(unique(mrids)) == 2
      res = runShortCircuit!(net; case = :max)
      # Two identical feeders in parallel halve Zk: Ik'' doubles exactly.
      @test only(res.rows).ik_kA ≈ 2000.0 / (sqrt(3.0) * 110.0) rtol = 1e-6

      # Removal safety: the suffix continues after the highest surviving id,
      # not after the record count — deleting the first record and adding a
      # third feeder must not recycle an id that is still in use.
      deleteat!(net.sc_sources.external_network_injections, 1)
      addExternalGrid!(net = net, busName = "B1", sk_max_MVA = 1000.0, rx_max = 0.1)
      mrids_after = [e.mrid for e in net.sc_sources.external_network_injections]
      @test length(unique(mrids_after)) == length(mrids_after)

      # A same-prefix neighbor bus (B10 vs B1) must not leak into B1's
      # suffix scan.
      addBus!(net = net, busName = "B10", vn_kV = 110.0)
      addExternalGrid!(net = net, busName = "B10", sk_max_MVA = 500.0)
      @test any(e -> e.mrid == "native-eni-B10", net.sc_sources.external_network_injections)
    end

    @testset "validation errors" begin
      base = () -> begin
        net = Net(name = "sc_val", baseMVA = 100.0)
        addBus!(net = net, busName = "B1", vn_kV = 110.0)
        net
      end
      @test_throws ArgumentError addExternalGrid!(net = base(), busName = "B1", sk_max_MVA = 0.0)
      @test_throws ArgumentError addExternalGrid!(net = base(), busName = "B1", sk_max_MVA = -10.0)
      @test_throws ArgumentError addExternalGrid!(net = base(), busName = "B1", sk_max_MVA = NaN)
      @test_throws ArgumentError addExternalGrid!(net = base(), busName = "B1", sk_max_MVA = 100.0, sk_min_MVA = 0.0)
      @test_throws ArgumentError addExternalGrid!(net = base(), busName = "B1", sk_max_MVA = 100.0, sk_min_MVA = 200.0)
      @test_throws ArgumentError addExternalGrid!(net = base(), busName = "B1", sk_max_MVA = 100.0, rx_max = -0.1)
      @test_throws ArgumentError addExternalGrid!(net = base(), busName = "B1", sk_max_MVA = 100.0, rx_min = -0.1)
    end

    @testset "copy regression: sc_sources survives net copying" begin
      net = _eg_fixture()
      addExternalGrid!(net = net, busName = "B2", sk_max_MVA = 2000.0, rx_max = 0.15)
      copy_net = deepcopy(net)
      @test length(copy_net.sc_sources.external_network_injections) == 1
      res = runShortCircuit!(copy_net; case = :max)
      row = first(r for r in res.rows if r.bus == "B2")
      @test row.status === :ok
      @test row.ik_kA ≈ first(r for r in runShortCircuit!(net; case = :max).rows if r.bus == "B2").ik_kA
    end

    @testset "external grid and distributed slack are mutually exclusive" begin
      # Both decide who covers the imbalance; combined, the source's import
      # would be forced to its zero participation share. The configuration
      # must reject the pair (same pattern as autodamp vs. trust region).
      mktempdir() do dir
        bad = joinpath(dir, "both.yaml")
        write(bad, "power_flow:\n  external_grid:\n    enabled: true\n  distributed_slack:\n    enabled: true\n")
        @test_throws ArgumentError load_sparlectra_config(bad; reload = true)
        ok = joinpath(dir, "one.yaml")
        write(ok, "power_flow:\n  external_grid:\n    enabled: true\n")
        cfg = load_sparlectra_config(ok; reload = true)
        @test cfg.powerflow.external_grid.enabled
        @test !cfg.powerflow.distributed_slack.enabled
      end
    end

    @testset "convertSlackToExternalGrid! demotes and stays consistent" begin
      net = _eg_fixture()
      note = convertSlackToExternalGrid!(net = net, sk_max_MVA = 500.0, rx_max = 0.1)
      @test occursin("slack bus B1", note)
      @test occursin("Sk'' = 500.0 MVA", note)
      @test occursin("B1__extgrid_int", note)
      # The demotion re-derives bus types itself (not only as a side effect
      # of the addExternalGrid! call): the terminal bus is an ordinary PQ
      # bus, the internal bus the only slack.
      aux = geNetBusIdx(net = net, busName = "B1__extgrid_int")
      @test net.slackVec == [aux]
      @test net.nodeVec[geNetBusIdx(net = net, busName = "B1")]._nodeType == Sparlectra.PQ
      _, erg = runpf!(net, 30, 1e-8, 0)
      @test erg == 0
      @test getNodeVm(net.nodeVec[geNetBusIdx(net = net, busName = "B1")]) < 1.02

      @test_throws ArgumentError convertSlackToExternalGrid!(net = Net(name = "noslack", baseMVA = 100.0), sk_max_MVA = 100.0)
      net2 = _eg_fixture()
      @test_throws ArgumentError convertSlackToExternalGrid!(net = net2, busName = "B2", sk_max_MVA = 100.0)
    end

    @testset "result print states the connection: slack vs. source" begin
      # Plain slack: the header names the slack bus, no source wording.
      plain = _eg_fixture()
      _, erg_p = runpf!(plain, 30, 1e-8, 0)
      @test erg_p == 0
      mktempdir() do dir
        printACPFlowResults(plain, 0.0, 3, 1e-8, true, dir)
        txt = read(joinpath(dir, "result_$(plain.name).txt"), String)
        @test occursin("Grid connection: slack bus B1", txt)
        @test !occursin("external-grid source", txt)
      end

      # Converted source: prose line with feeder data + internal slack, and
      # the type column reports SOURCE instead of SLACK on the internal bus.
      net = _eg_fixture()
      convertSlackToExternalGrid!(net = net, sk_max_MVA = 500.0, rx_max = 0.1)
      _, erg = runpf!(net, 30, 1e-8, 0)
      @test erg == 0
      mktempdir() do dir
        printACPFlowResults(net, 0.0, 3, 1e-8, true, dir)
        txt = read(joinpath(dir, "result_$(net.name).txt"), String)
        @test occursin("Grid connection: external-grid source at B1 (Sk'' = 500.0 MVA, R/X = 0.1; internal slack: B1__extgrid_int)", txt)
        @test occursin("SOURCE", txt)
        @test !occursin("| SLACK", txt)
        # header counts the reference models separately: this net's only
        # reference is the external-grid source, no plain slack remains
        @test occursin("Slack: 0 Source: 1", txt)
      end
      report = buildACPFlowReport(net; ite = 3, converged = true)
      aux_row = only(r for r in report.nodes if r.bus_name == "B1__extgrid_int")
      @test aux_row.type == "SOURCE"
    end

    @testset "internal impedance: non-ideal source (stage 2)" begin
      build = (; internal, sk) -> begin
        net = Net(name = "extgrid_int", baseMVA = 100.0)
        addBus!(net = net, busName = "B1", vn_kV = 110.0)
        addBus!(net = net, busName = "B2", vn_kV = 110.0)
        addPIModelACLine!(net = net, fromBus = "B1", toBus = "B2", r_pu = 0.01, x_pu = 0.05, b_pu = 0.0, status = 1)
        addProsumer!(net = net, busName = "B2", type = "ENERGYCONSUMER", p = 50.0, q = 10.0)
        addExternalGrid!(net = net, busName = "B1", vm_pu = 1.02, sk_max_MVA = sk, rx_max = 0.1, internal_impedance = internal)
        net
      end

      ideal = build(internal = false, sk = 3000.0)
      _, erg_i = runpf!(ideal, 30, 1e-8, 0)
      @test erg_i == 0
      vm_ideal = getNodeVm(ideal.nodeVec[1])

      # Stiff limit: a huge Sk must reproduce the ideal slack. The series
      # impedance scales as baseMVA/Sk, so sk = 1e10 puts the residual
      # voltage drop safely below the 1e-8 pu assertion (measured ~1.6e-8
      # at the task's original 1e9 — impedance, not tolerance, is the knob).
      stiff = build(internal = true, sk = 1.0e10)
      @test length(stiff.nodeVec) == 3
      @test stiff.slackVec == [3]                       # slack moved to the internal bus
      @test startswith(getCompName(stiff.nodeVec[3].comp), "Aux_")
      _, erg_s = runpf!(stiff, 30, 1e-8, 0)
      @test erg_s == 0
      @test abs(getNodeVm(stiff.nodeVec[1]) - vm_ideal) < 1e-8

      # Realistic Sk: terminal voltage droops and the angle shifts — the
      # terminal bus is now an ordinary solved bus.
      weak = build(internal = true, sk = 500.0)
      _, erg_w = runpf!(weak, 30, 1e-8, 0)
      @test erg_w == 0
      @test getNodeVm(weak.nodeVec[1]) < vm_ideal
      @test abs(weak.nodeVec[1]._va_deg) > 0.1

      # The feeder record stays anchored at the physical connection bus and
      # the auxiliary branch is inert in the SC network: the fault at B1
      # reproduces the declared current exactly.
      res = runShortCircuit!(weak; case = :max)
      row = first(r for r in res.rows if r.bus == "B1")
      @test row.ik_kA ≈ 500.0 / (sqrt(3.0) * 110.0) rtol = 1e-6

      # One internal-impedance external grid per bus.
      @test_throws ArgumentError addExternalGrid!(net = weak, busName = "B1", sk_max_MVA = 100.0, internal_impedance = true)
    end
  end
  return nothing
end
