# Copyright 2023–2026 Udo Schmitz
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

# file: test/test_transformer_phase_shift.jl

function run_transformer_phase_shift_tests()
  function _make_pi_trafo_branch(; shift_deg::Float64, ratio::Float64 = 1.0)::Branch
    net = Net(name = "phase_shift_probe", baseMVA = 100.0)
    addBus!(net = net, busName = "B1", vn_kV = 110.0)
    addBus!(net = net, busName = "B2", vn_kV = 110.0)
    addPIModelTrafo!(net = net, fromBus = "B1", toBus = "B2", r_pu = 0.01, x_pu = 0.10, b_pu = 0.02, status = 1, ratio = ratio, shift_deg = shift_deg)
    return only(net.branchVec)
  end

  @testset "Transformer tap/phase shift sign convention" begin
    y_ser = inv(0.01 + 0.10im)
    y_sh = 0.0 + 0.02im
    t = 1.05 * cis(deg2rad(7.5))

    expected_ff = (y_ser + 0.5 * y_sh) / abs2(t)
    expected_ft = -y_ser / conj(t)
    expected_tf = -y_ser / t
    expected_tt = y_ser + 0.5 * y_sh

    branch0 = _make_pi_trafo_branch(shift_deg = 7.5, ratio = 1.05)

    Y_ff, Y_ft, Y_tf, Y_tt = calcAdmittance(branch0, 110.0, 100.0)

    @test isapprox(Y_ff, expected_ff; atol = 1e-12, rtol = 0.0)
    @test isapprox(Y_ft, expected_ft; atol = 1e-12, rtol = 0.0)
    @test isapprox(Y_tf, expected_tf; atol = 1e-12, rtol = 0.0)
    @test isapprox(Y_tt, expected_tt; atol = 1e-12, rtol = 0.0)
  end

  @testset "Empirical phase-shift direction check (2-bus, Δφ=+5°)" begin
    # Minimal 2-bus setup for P_ab direction validation.
    # Fixed terminal voltages isolate transformer phase-shift effect.
    V = ComplexF64[1.0 + 0.0im, 1.0 + 0.0im]

    branch_phi0 = _make_pi_trafo_branch(shift_deg = 0.0, ratio = 1.0)
    branch_phi5 = _make_pi_trafo_branch(shift_deg = 5.0, ratio = 1.0)

    function p_from_to(branch::Branch, V::Vector{ComplexF64})
      Y_ff, Y_ft, _, _ = calcAdmittance(branch, 110.0, 100.0)
      I_f = Y_ff * V[1] + Y_ft * V[2]
      return real(V[1] * conj(I_f))
    end

    p0 = p_from_to(branch_phi0, V)
    p5 = p_from_to(branch_phi5, V)
    ΔP_from_to = p5 - p0

    # Observed with Sparlectra sign convention: +5° decreases P_ab.
    @test ΔP_from_to < 0.0

    # Corrected controller rule: do not assume sign convention.
    # Use empirical direction from the Δφ probe.
    phase_step = 0.5
    increase_direction = sign(ΔP_from_to)
    target_above = p0 + 0.05
    target_below = p0 - 0.05
    ctrl_step_when_below_target = (p0 < target_above) ? increase_direction * phase_step : -increase_direction * phase_step
    ctrl_step_when_above_target = (p0 < target_below) ? increase_direction * phase_step : -increase_direction * phase_step

    @test ctrl_step_when_below_target < 0.0
    @test ctrl_step_when_above_target > 0.0
  end

  return nothing
end
