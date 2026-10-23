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
# purpose: tests the AC line builders: the shunt conductance g of
#          addPIModelACLine! on branch and YBUS, the per-terminal shunt split,
#          and the physical input of
#          addACLine! (ohm and siemens totals) and addACLinePerKm! (values per
#          km, b from the capacitance at the system frequency of the net)

function run_piline_g_tests()
  # 0.50.0: addACLine! takes ohm and siemens totals, addACLinePerKm! values
  # per km. Before, addACLine! read r, x, b as per unit, multiplied them by
  # the length up to 1 km, read the siemens value derived from the
  # capacitance as per unit (Z_base too small) and set b = 0 without tanδ.
  @testset "physical line input: ohm, siemens and per km" begin (function ()
    vn, sbase = 110.0, 100.0
    zb = vn^2 / sbase
    newnet() = begin
      n = Net(name = "acline_physical", baseMVA = sbase)
      addBus!(net = n, busName = "B1", vn_kV = vn)
      addBus!(net = n, busName = "B2", vn_kV = vn)
      n
    end
    ybus(n) = createYBUS(net = n, sparse = false)

    # ohm and siemens against the same line in per unit: identical YBUS
    r_ohm, x_ohm, b_s, g_s = 2.4, 7.8, 5.7e-5, 1.2e-7
    phys = newnet()
    addACLine!(net = phys, fromBus = "B1", toBus = "B2", r_ohm = r_ohm, x_ohm = x_ohm, b_s = b_s, g_s = g_s, length_km = 20.0)
    seg = phys.linesAC[end]
    @test !isLinePIModel(seg)
    @test getLineRXBG(seg) == (r_ohm, x_ohm, b_s, g_s)
    @test seg.length_km == 20.0
    r_pu, x_pu, b_pu, g_pu = toPU_RXBG(r = r_ohm, x = x_ohm, g = g_s, b = b_s, v_kv = vn, baseMVA = sbase)
    @test getLineRXBG_pu(seg, vn, sbase) == (r_pu, x_pu, b_pu, g_pu)
    pu = newnet()
    addPIModelACLine!(net = pu, fromBus = "B1", toBus = "B2", r_pu = r_pu, x_pu = x_pu, b_pu = b_pu, g_pu = g_pu, status = 1)
    @test ybus(phys) == ybus(pu)
    # and the conversion is the one by hand: b in pu is b in siemens times Z_base
    @test isapprox(phys.branchVec[end].r_pu, r_ohm / zb; rtol = 1e-14)
    @test isapprox(phys.branchVec[end].b_pu, b_s * zb; rtol = 1e-14)
    # the length is information only
    nolen = newnet()
    addACLine!(net = nolen, fromBus = "B1", toBus = "B2", r_ohm = r_ohm, x_ohm = x_ohm, b_s = b_s, g_s = g_s)
    @test nolen.linesAC[end].length_km === nothing
    @test ybus(nolen) == ybus(phys)

    # per km against the totals: identical YBUS, for tanδ and for g per km
    len, r_km, x_km, c_km = 20.0, 0.12, 0.39, 9.5
    b_tot = 2.0 * pi * 50.0 * c_km * 1.0e-9 * len
    for (gkw, g_tot) in (((tanδ = 0.002,), 0.002 * b_tot), ((g_us_per_km = 0.05,), 0.05 * 1.0e-6 * len))
      perkm = newnet()
      addACLinePerKm!(; net = perkm, fromBus = "B1", toBus = "B2", length_km = len, r_ohm_per_km = r_km, x_ohm_per_km = x_km, c_nf_per_km = c_km, gkw...)
      tot = newnet()
      addACLine!(net = tot, fromBus = "B1", toBus = "B2", r_ohm = r_km * len, x_ohm = x_km * len, b_s = b_tot, g_s = g_tot, length_km = len)
      @test ybus(perkm) == ybus(tot)
      @test getLineRXBG(perkm.linesAC[end]) == getLineRXBG(tot.linesAC[end])
    end

    # b from the capacitance without tanδ is not zero (it was before 0.50.0)
    noloss = newnet()
    addACLinePerKm!(net = noloss, fromBus = "B1", toBus = "B2", length_km = len, r_ohm_per_km = r_km, x_ohm_per_km = x_km, c_nf_per_km = c_km)
    br = noloss.branchVec[end]
    @test br.g_pu == 0.0
    # 2π·50 Hz·9.5 nF/km·20 km = 59.69 µS, times Z_base 121 ohm
    @test isapprox(br.b_pu, 2.0 * pi * 50.0 * 9.5e-9 * 20.0 * 121.0; rtol = 1e-14)
    @test isapprox(br.b_pu, 0.0072225; rtol = 1e-4)

    # the frequency is a property of the network (net.nominal_frequency_hz):
    # 60 Hz gives 6/5 times b, x is untouched (no input in inductance exists)
    @test noloss.nominal_frequency_hz == 50.0
    n60 = Net(name = "acline_60hz", baseMVA = sbase, nominal_frequency_hz = 60.0)
    addBus!(net = n60, busName = "B1", vn_kV = vn)
    addBus!(net = n60, busName = "B2", vn_kV = vn)
    addACLinePerKm!(net = n60, fromBus = "B1", toBus = "B2", length_km = len, r_ohm_per_km = r_km, x_ohm_per_km = x_km, c_nf_per_km = c_km)
    @test isapprox(n60.branchVec[end].b_pu / br.b_pu, 6 / 5; rtol = 1e-14)
    @test n60.branchVec[end].x_pu == br.x_pu
    @test n60.branchVec[end].r_pu == br.r_pu
    # a new Net takes model.nominal_frequency_hz; the line then reads the net
    cfg60 = SparlectraConfig(model = ModelConfig(nominal_frequency_hz = 60.0))
    n60cfg = with_sparlectra_config(newnet, cfg60)
    @test n60cfg.nominal_frequency_hz == 60.0
    addACLinePerKm!(net = n60cfg, fromBus = "B1", toBus = "B2", length_km = len, r_ohm_per_km = r_km, x_ohm_per_km = x_km, c_nf_per_km = c_km)
    @test n60cfg.branchVec[end].b_pu == n60.branchVec[end].b_pu
    # the run header names the frequency only where a capacitance became a
    # susceptance: a line in ohm and siemens and a per-km line with c = 0
    # never read it
    @test Sparlectra.nominal_frequency_line(noloss) == "50 Hz (configuration)"
    @test Sparlectra.nominal_frequency_line(n60cfg) == "60 Hz (configuration)"
    @test Sparlectra.nominal_frequency_source(phys) == "not_used"
    zero_c = newnet()
    addACLinePerKm!(net = zero_c, fromBus = "B1", toBus = "B2", length_km = len, r_ohm_per_km = r_km, x_ohm_per_km = x_km, c_nf_per_km = 0.0)
    @test Sparlectra.nominal_frequency_line(zero_c) == "not used (no capacitance in farad)"
    @test_throws ArgumentError Net(name = "f0", baseMVA = sbase, nominal_frequency_hz = 0.0)
    @test_throws ArgumentError Net(name = "fneg", baseMVA = sbase, nominal_frequency_hz = -50.0)
    @test Sparlectra.ModelConfig(Dict{String,Any}("model" => Dict{String,Any}("nominal_frequency_hz" => 60))).nominal_frequency_hz == 60.0
    @test_throws ArgumentError Sparlectra.ModelConfig(Dict{String,Any}("model" => Dict{String,Any}("nominal_frequency_hz" => 0.0)))
    # the SCF writer states the frequency of the net, whatever the active
    # configuration says, and turns b into the capacitance with it; a file
    # without f_nom is read with the configured one
    root = Sparlectra.net_to_scf(n60)
    @test root["sparlectra"]["meta"]["f_nom"] == 60.0
    @test isapprox(only(root["data"]["line"])["c1"], c_km * 1.0e-9 * len; rtol = 1e-12)
    delete!(root["sparlectra"]["meta"], "f_nom")
    back = with_sparlectra_config(() -> Sparlectra.scf_to_net(root), cfg60)
    @test isapprox(back.branchVec[end].b_pu, n60.branchVec[end].b_pu; rtol = 1e-14)

    # mixed or incomplete input is refused
    @test_throws ArgumentError addACLine!(net = newnet(), fromBus = "B1", toBus = "B2", r_ohm = 1.0, x_ohm = 2.0, c_nf_per_km = 9.5)
    @test_throws ArgumentError addACLinePerKm!(net = newnet(), fromBus = "B1", toBus = "B2", length_km = len, r_ohm_per_km = r_km, x_ohm_per_km = x_km, c_nf_per_km = c_km, tanδ = 0.002, g_us_per_km = 0.05)
    @test_throws ArgumentError addACLinePerKm!(net = newnet(), fromBus = "B1", toBus = "B2", length_km = 0.0, r_ohm_per_km = r_km, x_ohm_per_km = x_km, c_nf_per_km = c_km)
  end)() end

  # Y11 = y_ser + (g + jb)/2 and Y12 = -y_ser by hand; the second row is
  # g > 0 with b = 0 (network equivalents), which must reach the matrix too
  @testset "addPIModelACLine! end-to-end into branch and YBUS" begin (function ()
    yser = 1.0 / (0.01 + 0.1im)
    for (b_pu, g_pu) in ((0.02, 0.004), (0.0, 0.004))
      net = Net(name = "piline_g_test", baseMVA = 100.0)
      addBus!(net = net, busName = "B1", vn_kV = 110.0)
      addBus!(net = net, busName = "B2", vn_kV = 110.0)
      addPIModelACLine!(net = net, fromBus = "B1", toBus = "B2", r_pu = 0.01, x_pu = 0.1, b_pu = b_pu, g_pu = g_pu, status = 1)
      br = net.branchVec[end]
      @test br.g_pu == g_pu
      @test br.b_pu == b_pu

      Y = createYBUS(net = net, sparse = false)
      @test isapprox(Y[1, 1], yser + 0.5 * (g_pu + b_pu * im); atol = 1e-12)
      @test isapprox(Y[1, 2], -yser; atol = 1e-12)
    end
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
end
