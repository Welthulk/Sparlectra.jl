# file: test/test_piline_g.jl
# Tests for the optional shunt conductance (g / g_pu) on PI-model AC lines
# (CGMES import precondition B1, task_plan_cgmes_import.md §7.6).

function run_piline_g_tests()
  @testset "PI-model line shunt conductance" begin
    @testset "ACLineSegment explicit g (PI path)" begin
      seg = ACLineSegment(vn_kv = 110.0, from = 1, to = 2, length = 1.0, r = 0.01, x = 0.1, b = 0.02, g = 0.004, paramsBasedOnLength = false, isPIModel = true)
      r, x, b, g = getLineRXBG(seg)
      @test r == 0.01 && x == 0.1
      @test b == 0.02
      @test g == 0.004
    end

    @testset "edge case g > 0 with b = 0 (network equivalents)" begin
      seg = ACLineSegment(vn_kv = 110.0, from = 1, to = 2, length = 1.0, r = 0.01, x = 0.1, b = 0.0, g = 0.004, paramsBasedOnLength = false, isPIModel = true)
      _, _, b, g = getLineRXBG(seg)
      @test b == 0.0
      @test g == 0.004
    end

    @testset "default stays 0.0 (backward compatible)" begin
      seg = ACLineSegment(vn_kv = 110.0, from = 1, to = 2, length = 1.0, r = 0.01, x = 0.1, b = 0.02, paramsBasedOnLength = false, isPIModel = true)
      _, _, _, g = getLineRXBG(seg)
      @test g == 0.0
    end

    @testset "tanδ derivation unaffected, explicit g wins" begin
      # legacy tanδ route still derives g from c_nf_per_km
      derived = ACLineSegment(vn_kv = 110.0, from = 1, to = 2, length = 10.0, r = 0.1, x = 0.4, c_nf_per_km = 10.0, tanδ = 0.02)
      @test derived.g > 0.0
      # explicit g suppresses the derivation
      explicit = ACLineSegment(vn_kv = 110.0, from = 1, to = 2, length = 10.0, r = 0.1, x = 0.4, b = 3.0e-6, g = 5.0e-7, c_nf_per_km = 10.0, tanδ = 0.02)
      @test explicit.g == 5.0e-7
      @test explicit.b == 3.0e-6
    end

    @testset "addPIModelACLine! end-to-end into branch and YBUS" begin
      net = Net(name = "piline_g_test", baseMVA = 100.0)
      addBus!(net = net, busName = "B1", vn_kV = 110.0)
      addBus!(net = net, busName = "B2", vn_kV = 110.0)
      addPIModelACLine!(net = net, fromBus = "B1", toBus = "B2", r_pu = 0.01, x_pu = 0.1, b_pu = 0.02, g_pu = 0.004, status = 1)
      br = net.branchVec[end]
      @test br.g_pu == 0.004
      @test br.b_pu == 0.02

      Y = createYBUS(net = net, sparse = false)
      yser = 1.0 / (0.01 + 0.1im)
      @test isapprox(Y[1, 1], yser + 0.5 * (0.004 + 0.02im); atol = 1e-12)
      @test isapprox(Y[1, 2], -yser; atol = 1e-12)
    end

    @testset "addPIModelACLine! without g_pu unchanged" begin
      net = Net(name = "piline_g_default", baseMVA = 100.0)
      addBus!(net = net, busName = "B1", vn_kV = 110.0)
      addBus!(net = net, busName = "B2", vn_kV = 110.0)
      addPIModelACLine!(net = net, fromBus = "B1", toBus = "B2", r_pu = 0.01, x_pu = 0.1, b_pu = 0.02, status = 1)
      @test net.branchVec[end].g_pu == 0.0
    end
  end
end
