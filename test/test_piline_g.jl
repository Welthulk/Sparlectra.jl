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

# file: test/test_piline_g.jl
# purpose: tests the optional shunt conductance g on PI-model AC lines:
#          explicit g on ACLineSegment and addPIModelACLine!, defaults, and
#          the end-to-end effect on branch and YBUS
# Tests for the optional shunt conductance (g / g_pu) on PI-model AC lines
# (CGMES import precondition; see docs/src/cgmes_import.md).

function run_piline_g_tests()
  @testset "PI-model line shunt conductance" begin (function ()
    @testset "ACLineSegment explicit g (PI path)" begin (function ()
      seg = ACLineSegment(vn_kv = 110.0, from = 1, to = 2, length = 1.0, r = 0.01, x = 0.1, b = 0.02, g = 0.004, paramsBasedOnLength = false, isPIModel = true)
      r, x, b, g = getLineRXBG(seg)
      @test r == 0.01 && x == 0.1
      @test b == 0.02
      @test g == 0.004
    end)() end

    @testset "edge case g > 0 with b = 0 (network equivalents)" begin (function ()
      seg = ACLineSegment(vn_kv = 110.0, from = 1, to = 2, length = 1.0, r = 0.01, x = 0.1, b = 0.0, g = 0.004, paramsBasedOnLength = false, isPIModel = true)
      _, _, b, g = getLineRXBG(seg)
      @test b == 0.0
      @test g == 0.004
    end)() end

    @testset "default stays 0.0 (backward compatible)" begin (function ()
      seg = ACLineSegment(vn_kv = 110.0, from = 1, to = 2, length = 1.0, r = 0.01, x = 0.1, b = 0.02, paramsBasedOnLength = false, isPIModel = true)
      _, _, _, g = getLineRXBG(seg)
      @test g == 0.0
    end)() end

    @testset "tanδ derivation unaffected, explicit g wins" begin (function ()
      # legacy tanδ route still derives g from c_nf_per_km
      derived = ACLineSegment(vn_kv = 110.0, from = 1, to = 2, length = 10.0, r = 0.1, x = 0.4, c_nf_per_km = 10.0, tanδ = 0.02)
      @test derived.g > 0.0
      # explicit g suppresses the derivation
      explicit = ACLineSegment(vn_kv = 110.0, from = 1, to = 2, length = 10.0, r = 0.1, x = 0.4, b = 3.0e-6, g = 5.0e-7, c_nf_per_km = 10.0, tanδ = 0.02)
      @test explicit.g == 5.0e-7
      @test explicit.b == 3.0e-6
    end)() end

    @testset "addPIModelACLine! end-to-end into branch and YBUS" begin (function ()
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
    end)() end

    @testset "addPIModelACLine! without g_pu unchanged" begin (function ()
      net = Net(name = "piline_g_default", baseMVA = 100.0)
      addBus!(net = net, busName = "B1", vn_kV = 110.0)
      addBus!(net = net, busName = "B2", vn_kV = 110.0)
      addPIModelACLine!(net = net, fromBus = "B1", toBus = "B2", r_pu = 0.01, x_pu = 0.1, b_pu = 0.02, status = 1)
      @test net.branchVec[end].g_pu == 0.0
    end)() end

    # per-terminal charging admittance (0.20.0): the symmetric half is the
    # default, an explicit split stamps each arm on its own end, the from
    # arm through |t|^2, and a partial keyword set is refused
    @testset "asymmetric shunt split on lines and transformers" begin (function ()
      net = Net(name = "piline_split", baseMVA = 100.0)
      addBus!(net = net, busName = "B1", vn_kV = 110.0)
      addBus!(net = net, busName = "B2", vn_kV = 110.0)
      addBus!(net = net, busName = "B3", vn_kV = 110.0)
      addPIModelACLine!(net = net, fromBus = "B1", toBus = "B2", r_pu = 0.01, x_pu = 0.1, b_pu = 0.0, status = 1, g_from_pu = 0.004, b_from_pu = 0.03, g_to_pu = 0.0, b_to_pu = 0.01)
      ln = net.branchVec[end]
      @test (ln.g_from_pu, ln.b_from_pu, ln.g_to_pu, ln.b_to_pu) == (0.004, 0.03, 0.0, 0.01)
      # the totals are the sums, whatever total the caller passed
      @test ln.b_pu == 0.04 && ln.g_pu == 0.004
      @test !Sparlectra.has_symmetric_shunt(ln)
      addPIModelTrafo!(net = net, fromBus = "B2", toBus = "B3", r_pu = 0.005, x_pu = 0.05, b_pu = 0.0, status = 1, ratio = 1.05, shift_deg = 0.0, g_from_pu = 0.002, b_from_pu = -0.01, g_to_pu = 0.0, b_to_pu = 0.0)
      tr = net.branchVec[end]
      @test tr.b_pu == -0.01 && tr.g_pu == 0.002
      Y = createYBUS(net = net, sparse = false)
      yl = 1.0 / (0.01 + 0.1im)
      yt = 1.0 / (0.005 + 0.05im)
      @test isapprox(Y[1, 1], yl + (0.004 + 0.03im); atol = 1e-12)
      @test isapprox(Y[2, 2], yl + 0.01im + (yt + (0.002 - 0.01im)) / 1.05^2; atol = 1e-12)
      @test isapprox(Y[3, 3], yt; atol = 1e-12)
      @test isapprox(Y[1, 2], -yl; atol = 1e-12)
      # the flow helpers use the arm of the end they leave: the reactive
      # charging at a flat 1 pu profile differs by the arm difference
      V = fill(1.0 + 0.0im, 3)
      s12 = Sparlectra._closed_branch_flow_pu(V, 1, 2, ln, 1)
      s21 = Sparlectra._closed_branch_flow_pu(V, 2, 1, ln, 2)
      @test isapprox(imag(s12) - imag(s21), -(0.03 - 0.01); atol = 1e-12)
      @test isapprox(real(s12) - real(s21), 0.004; atol = 1e-12)
      # the one writer keeps the totals in step
      Sparlectra.set_branch_shunt!(ln; g_from_pu = 0.0, b_from_pu = 0.02, g_to_pu = 0.0, b_to_pu = 0.02)
      @test ln.b_pu == 0.04 && Sparlectra.has_symmetric_shunt(ln)
      Sparlectra.set_branch_shunt_total!(ln; g_pu = 0.002, b_pu = 0.05)
      @test (ln.g_from_pu, ln.b_from_pu, ln.g_to_pu, ln.b_to_pu) == (0.001, 0.025, 0.001, 0.025)
      # a partial keyword set names the branch
      @test_throws ArgumentError addPIModelACLine!(net = net, fromBus = "B1", toBus = "B3", r_pu = 0.01, x_pu = 0.1, b_pu = 0.02, status = 1, b_from_pu = 0.01)
      @test_throws ArgumentError BranchModel(r_pu = 0.01, x_pu = 0.1, b_pu = 0.02, g_pu = 0.0, ratio = 0.0, angle = 0.0, g_to_pu = 0.0)
      bm = BranchModel(r_pu = 0.01, x_pu = 0.1, b_pu = 0.02, g_pu = 0.0, ratio = 0.0, angle = 0.0)
      @test (bm.b_from_pu, bm.b_to_pu) == (0.01, 0.01)
      # no code path writes a branch total alone: the setters in branch.jl
      # are the only assignments to b_pu and g_pu of a branch
      src_root = joinpath(dirname(@__DIR__), "src")
      offenders = String[]
      for (root, _, files) in walkdir(src_root), f in files
        endswith(f, ".jl") || continue
        path = joinpath(root, f)
        relpath(path, src_root) == "branch.jl" && continue
        for (i, line) in enumerate(eachline(path))
          occursin(r"\.(b_pu|g_pu)\s*=[^=]", line) && push!(offenders, "$(relpath(path, src_root)):$(i)")
        end
      end
      @test isempty(offenders)
    end)() end
  end)() end
end
