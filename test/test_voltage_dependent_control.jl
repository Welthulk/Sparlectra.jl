# Copyright 2023–2026 Udo Schmitz

using Test
using Sparlectra
using Printf

function run_voltage_dependent_control_tests()
  @testset "Voltage dependent P(U)/Q(U) control" begin
    @testset "Characteristic evaluation" begin
      ch = PiecewiseLinearCharacteristic([(0.95, -0.2), (1.0, 0.0), (1.05, 0.2)])

      v_mid, dv_mid = evaluate_characteristic(ch, 0.975)
      @test isapprox(v_mid, -0.1; atol = 1e-12)
      @test isapprox(dv_mid, 4.0; atol = 1e-12)

      v_lo, dv_lo = evaluate_characteristic(ch, 0.90)
      @test isapprox(v_lo, -0.2; atol = 1e-12)
      @test dv_lo == 0.0

      v_hi, dv_hi = evaluate_characteristic(ch, 1.10)
      @test isapprox(v_hi, 0.2; atol = 1e-12)
      @test dv_hi == 0.0

      qu = QUController(ch, -0.05, 0.05)
      q_val, q_slope = evaluate_controller(qu, 1.05)
      @test isapprox(q_val, 0.05; atol = 1e-12)
      @test q_slope == 0.0
    end

    @testset "Physical-unit controller inputs (kV, MW, MVAr)" begin
      ch_qu = make_characteristic([(104.5, 30.0), (110.0, 0.0), (115.5, -20.0)]; voltage_unit = :kV, value_unit = :MVAr, vn_kV = 110.0, sbase_MVA = 100.0)
      ch_pu = make_characteristic([(104.5, 20.0), (110.0, 10.0), (115.5, 0.0)]; voltage_unit = :kV, value_unit = :MW, vn_kV = 110.0, sbase_MVA = 100.0)

      @test ch_qu.points == [(0.95, 0.3), (1.0, 0.0), (1.05, -0.2)]
      @test ch_pu.points == [(0.95, 0.2), (1.0, 0.1), (1.05, 0.0)]

      qu = QUController(ch_qu; qmin_MVAr = -50.0, qmax_MVAr = 50.0, sbase_MVA = 100.0)
      pu = PUController(ch_pu; pmin_MW = 0.0, pmax_MW = 50.0, sbase_MVA = 100.0)

      @test isapprox(qu.qmin_pu, -0.5; atol = 1e-12)
      @test isapprox(qu.qmax_pu, 0.5; atol = 1e-12)
      @test isapprox(pu.pmin_pu, 0.0; atol = 1e-12)
      @test isapprox(pu.pmax_pu, 0.5; atol = 1e-12)
    end

    @testset "No-control compatibility and solver integration" begin
      net_plain = createTest3BusNet()
      V_plain = buildVoltageVector(net_plain)
      S_plain = buildComplexSVec(net_plain)
      S_eval_plain, dP_plain, dQ_plain = buildControlledSVec(net_plain, V_plain)
      @test S_plain ≈ S_eval_plain
      @test all(iszero, dP_plain)
      @test all(iszero, dQ_plain)

      net = Net(name = "qu_control_case", baseMVA = 100.0)
      addBus!(net = net, busName = "B1", vn_kV = 110.0)
      addBus!(net = net, busName = "B2", vn_kV = 110.0)
      addACLine!(net = net, fromBus = "B1", toBus = "B2", length = 20.0, r = 0.05, x = 0.4)

      addProsumer!(net = net, busName = "B1", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "B1")
      qu_curve = PiecewiseLinearCharacteristic([(0.95, 0.4), (1.0, 0.0), (1.05, -0.2)])
      addProsumer!(net = net, busName = "B2", type = "SYNCHRONOUSMACHINE", p = 20.0, q = 0.0, qu_controller = QUController(qu_curve, -0.5, 0.5))
      addProsumer!(net = net, busName = "B2", type = "ENERGYCONSUMER", p = 60.0, q = 20.0)

      ite, erg = runpf!(net, 30, 1e-9, 0; method = :rectangular, opt_sparse = true)
      @test erg == 0
      @test ite <= 30

      V = buildVoltageVector(net)
      vm = abs(V[2])
      q_ctrl_pu, _ = evaluate_controller(net.prosumpsVec[2].quController, vm)
      S_eval, _, _ = buildControlledSVec(net, V)
      qnet_expected = q_ctrl_pu * net.baseMVA - 20.0
      @test isapprox(imag(S_eval[2]) * net.baseMVA, qnet_expected; atol = 1e-8)
    end

    @testset "Result printout shows Control column" begin
      net = Net(name = "control_print_case", baseMVA = 100.0)
      for b in ("B1", "B2", "B3", "B4")
        addBus!(net = net, busName = b, vn_kV = 110.0)
      end

      addACLine!(net = net, fromBus = "B1", toBus = "B2", length = 20.0, r = 0.02, x = 0.2)
      addACLine!(net = net, fromBus = "B1", toBus = "B3", length = 20.0, r = 0.02, x = 0.2)
      addACLine!(net = net, fromBus = "B1", toBus = "B4", length = 20.0, r = 0.02, x = 0.2)

      addProsumer!(net = net, busName = "B1", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "B1")
      addProsumer!(net = net, busName = "B2", type = "SYNCHRONOUSMACHINE", p = 15.0, q = 0.0, qu_controller = QUController(make_characteristic([(0.95, 0.2), (1.05, -0.2)]), -0.5, 0.5))
      addProsumer!(net = net, busName = "B3", type = "SYNCHRONOUSMACHINE", p = 0.0, q = 5.0, pu_controller = PUController(make_characteristic([(0.95, 0.2), (1.05, 0.0)]), -0.5, 0.5))
      addProsumer!(net = net, busName = "B4", type = "SYNCHRONOUSMACHINE", p = 0.0, q = 0.0, qu_controller = QUController(make_characteristic([(0.95, 0.1), (1.05, -0.1)]), -0.5, 0.5), pu_controller = PUController(make_characteristic([(0.95, 0.1), (1.05, 0.0)]), -0.5, 0.5))
      addProsumer!(net = net, busName = "B2", type = "ENERGYCONSUMER", p = 30.0, q = 12.0)
      addProsumer!(net = net, busName = "B3", type = "ENERGYCONSUMER", p = 20.0, q = 8.0)
      addProsumer!(net = net, busName = "B4", type = "ENERGYCONSUMER", p = 10.0, q = 4.0)

      _, erg = runpf!(net, 30, 1e-8, 0; method = :rectangular)
      @test erg == 0

      tmp = mktempdir()
      printACPFlowResults(net, 0.0, 1, 1e-8, true, tmp; converged = true, solver = :rectangular)
      out = read(joinpath(tmp, "result_$(net.name).txt"), String)
      V = buildVoltageVector(net)
      controlled_idx = findfirst(ps -> !isnothing(ps.puController) && !isnothing(ps.quController), net.prosumpsVec)
      @test !isnothing(controlled_idx)
      ps_ctrl = net.prosumpsVec[controlled_idx]
      vm_ctrl = abs(V[getPosumerBusIndex(ps_ctrl)])
      p_ctrl_pu, _ = evaluate_controller(ps_ctrl.puController, vm_ctrl)
      q_ctrl_pu, _ = evaluate_controller(ps_ctrl.quController, vm_ctrl)
      p_ctrl_mw = p_ctrl_pu * net.baseMVA
      q_ctrl_mvar = q_ctrl_pu * net.baseMVA

      @test occursin("Control", out)
      @test occursin("Q(U)", out)
      @test occursin("P(U)", out)
      @test occursin("Q(U), P(U)", out)
      @test occursin(" - ", out)
      @test occursin(@sprintf("%10.3f", p_ctrl_mw), out)
      @test occursin(@sprintf("%10.3f", q_ctrl_mvar), out)
    end
  end

  return true
end
