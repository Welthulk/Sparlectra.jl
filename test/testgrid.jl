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

# file: test/testgrid.jl
# purpose: grid and power-flow tests: transformer and network validation,
#          MATPOWER import, PF scenarios, links and reporting, plus the
#          extended parser group

# Author: Udo Schmitz (https://github.com/Welthulk)

# Nameplate to equivalent circuit (0.32.0): the short-circuit loss path carried
# a unit factor of 1e-9 instead of 1e-3 (r_k of 1.5e-6 ohm instead of 1.5
# ohm), the iron-loss conductance was taken in kW/kV^2 (a thousand times too
# large) so the magnetizing susceptance fell into a silent catch and came out
# as 0, and add2WTrafo! passed pfe_kw as the short-circuit loss with p0_kW = 0,
# so the magnetizing branch it promised was always empty. Expected values from
# the textbook formulas for a 110 kV, 40 MVA unit: v_k 12 %, v_kr 0.5 %
# (P_k 200 kW gives the same r_k), P_fe 30 kW, i_0 0.5 %.
function test_transformer_nameplate_conversion()::Bool
  zb = 110.0^2 / 40.0                      # 302.5 ohm
  r_expected = 0.005 * zb                  # 1.5125 ohm
  x_expected = sqrt((0.12 * zb)^2 - r_expected^2)
  g_expected = 30e3 / 110e3^2              # 2.479e-6 S
  y_expected = 0.005 * 40e6 / 110e3^2      # 1.653e-5 S
  b_expected = -sqrt(y_expected^2 - g_expected^2)
  for md in (
    Sparlectra.TransformerModelParameters(sn_MVA = 40.0, vk_percent = 12.0, vkr_percent = 0.5, i0_percent = 0.5, p0_kW = 30.0),
    Sparlectra.TransformerModelParameters(sn_MVA = 40.0, vk_percent = 12.0, pk_kW = 200.0, i0_percent = 0.5, p0_kW = 30.0),
  )
    r, x, g, b = Sparlectra.calcTransformerRXGB(110.0, md)
    isapprox(r, r_expected; rtol = 1e-12) || return false
    isapprox(x, x_expected; rtol = 1e-12) || return false
    isapprox(g, g_expected; rtol = 1e-12) || return false
    isapprox(b, b_expected; rtol = 1e-12) || return false
  end
  # an inconsistent nameplate (no-load current below the iron-loss current)
  # says so instead of swallowing a DomainError
  md_bad = Sparlectra.TransformerModelParameters(sn_MVA = 40.0, vk_percent = 12.0, vkr_percent = 0.5, i0_percent = 0.001, p0_kW = 30.0)
  _, _, g_bad, b_bad = @test_logs (:warn, r"no-load current i0 = 0.001 %") Sparlectra.calcTransformerRXGB(110.0, md_bad)
  (isapprox(g_bad, g_expected; rtol = 1e-12) && b_bad == 0.0) || return false
  # the same nameplate through add2WTrafo!: the branch carries the magnetizing
  # branch on the 110 kV / 100 MVA base (g_pu = g * Vn^2 / S_base)
  net = Net(name = "nameplate", baseMVA = 100.0)
  addBus!(net = net, busName = "HV", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
  addBus!(net = net, busName = "LV", vn_kV = 20.0, vm_pu = 1.0, va_deg = 0.0)
  add2WTrafo!(net = net, fromBus = "HV", toBus = "LV", sn_mva = 40.0, vk_percent = 12.0, vkr_percent = 0.5, pfe_kw = 30.0, i0_percent = 0.5)
  br = net.branchVec[1]
  zb100 = 110.0^2 / 100.0
  isapprox(br.r_pu, r_expected / zb100; rtol = 1e-12) || return false
  isapprox(br.x_pu, x_expected / zb100; rtol = 1e-12) || return false
  isapprox(br.g_pu, g_expected * zb100; rtol = 1e-12) || return false
  isapprox(br.b_pu, b_expected * zb100; rtol = 1e-12) || return false
  # a unit without loss data stays lossless and without a magnetizing branch
  add2WTrafo!(net = net, fromBus = "HV", toBus = "LV", sn_mva = 40.0, vk_percent = 12.0, vkr_percent = 0.0, pfe_kw = 0.0, i0_percent = 0.0)
  br0 = net.branchVec[2]
  return br0.r_pu == 0.0 && br0.g_pu == 0.0 && br0.b_pu == 0.0 && isapprox(br0.x_pu, 0.12 * zb / zb100; rtol = 1e-12)
end

# updatePQShunt! must set the same per-unit admittance as the constructor path
# (addShunt!): it used to store G and B in Siemens (p / vn_kV^2) next to a
# per-unit y, so the two writers of these fields disagreed. Both shunts must
# give the same power flow, and `show` must write its whole text to its own
# stream (a `print(IO, ...)` sent part of it to stdout).
function test_shunt_update_matches_constructor()::Bool
  function shunt_net(p, q)
    net = Net(name = "shunt_update", baseMVA = 100.0)
    addBus!(net = net, busName = "S", vn_kV = 110.0)
    addBus!(net = net, busName = "L", vn_kV = 110.0)
    addPIModelACLine!(net = net, fromBus = "S", toBus = "L", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
    addProsumer!(net = net, busName = "S", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "S")
    addProsumer!(net = net, busName = "L", type = "ENERGYCONSUMER", p = 40.0, q = 10.0)
    addShunt!(net = net, busName = "L", pShunt = p, qShunt = q)
    return net
  end
  a = shunt_net(2.0, 8.0)
  b = shunt_net(0.0, 0.0)
  updatePQShunt!(b.shuntVec[1], 2.0, 8.0)
  sa, sb = a.shuntVec[1], b.shuntVec[1]
  same_y = sa.y_pu_shunt == sb.y_pu_shunt && sa.G_shunt == sb.G_shunt && sa.B_shunt == sb.B_shunt
  # 8 MVar at 1 pu on a 100 MVA base is 0.08 pu (the Siemens value would be 8 / 110^2)
  per_unit = isapprox(sb.G_shunt, 0.02; atol = 1e-15) && isapprox(sb.B_shunt, 0.08; atol = 1e-15)
  _, ea = runpf!(a, 30, 1e-10, 0)
  _, eb = runpf!(b, 30, 1e-10, 0)
  dvm = maximum(abs(x._vm_pu - y._vm_pu) for (x, y) in zip(a.nodeVec, b.nodeVec))
  dva = maximum(abs(x._va_deg - y._va_deg) for (x, y) in zip(a.nodeVec, b.nodeVec))
  same_pf = ea == 0 && eb == 0 && dvm < 1e-12 && dva < 1e-10
  shown = sprint(show, sb)
  whole_show = occursin(string(sb.y_pu_shunt), shown)
  ok = same_y && per_unit && same_pf && whole_show
  ok || @warn "updatePQShunt! and the constructor disagree" same_y per_unit same_pf dvm dva whole_show
  return ok
end

function test_matpower_flatstart_uses_generator_voltage_setpoints()::Bool
  mpc = (
    name = "case_flatstart_vg_setpoints",
    baseMVA = 100.0,
    bus = [
      1 3 0.0 0.0 0.0 0.0 1 0.97 12.0 110.0 1 1.1 0.9
      2 1 20.0 8.0 0.0 0.0 1 0.93 -7.0 110.0 1 1.1 0.9
      3 2 0.0 0.0 0.0 0.0 1 0.96 5.0 110.0 1 1.1 0.9
    ],
    gen = [
      1 30.0 0.0 999.0 -999.0 1.04 100.0 1 999.0 0.0 0 0 0 0 0 0 0 0 0 0 0
      3 25.0 0.0 999.0 -999.0 1.02 100.0 1 999.0 0.0 0 0 0 0 0 0 0 0 0 0 0
    ],
    branch = [
      1 2 0.01 0.05 0.0 9999.0 0.0 0.0 0.0 0.0 1 -60.0 60.0
      2 3 0.01 0.05 0.0 9999.0 0.0 0.0 0.0 0.0 1 -60.0 60.0
    ],
    gencost = nothing,
    bus_name = nothing,
  )

  net = Sparlectra.createNetFromMatPowerCase(mpc = mpc, log = false, flatstart = true)
  V0, slack_idx = Sparlectra.initialVrect(net; flatstart = true)
  busVec, _ = Sparlectra.getBusData(net.nodeVec, net.baseMVA, true; net = net)

  # reference_vm_pu / reference_va_deg replace the slack values of the file
  net_ref = Sparlectra.createNetFromMatPowerCase(mpc = mpc, log = false, flatstart = true, reference_vm_pu = 1.01, reference_va_deg = -3.0)
  V0_ref, slack_ref = Sparlectra.initialVrect(net_ref; flatstart = true)
  busVec_ref, _ = Sparlectra.getBusData(net_ref.nodeVec, net_ref.baseMVA, true; net = net_ref)
  reference_override_ok = slack_ref == 1 && isapprox(abs(V0_ref[1]), 1.01; atol = 1e-12) && isapprox(angle(V0_ref[1]), deg2rad(-3.0); atol = 1e-12) && isapprox(busVec_ref[1].vm_pu, 1.01; atol = 1e-12) && isapprox(busVec_ref[1].va_rad, deg2rad(-3.0); atol = 1e-12)

  return slack_idx == 1 &&
         isapprox(abs(V0[1]), 1.04; atol = 1e-12) &&
         isapprox(angle(V0[1]), deg2rad(12.0); atol = 1e-12) &&
         isapprox(abs(V0[2]), 1.0; atol = 1e-12) &&
         isapprox(angle(V0[2]), deg2rad(12.0); atol = 1e-12) &&
         isapprox(abs(V0[3]), 1.02; atol = 1e-12) &&
         isapprox(angle(V0[3]), deg2rad(12.0); atol = 1e-12) &&
         isapprox(busVec[1].vm_pu, 1.04; atol = 1e-12) &&
         isapprox(busVec[1].va_rad, deg2rad(12.0); atol = 1e-12) &&
         isapprox(busVec[2].vm_pu, 1.0; atol = 1e-12) &&
         isapprox(busVec[2].va_rad, 0.0; atol = 1e-12) &&
         isapprox(busVec[3].vm_pu, 1.02; atol = 1e-12) &&
         isapprox(busVec[3].va_rad, 0.0; atol = 1e-12) &&
         reference_override_ok
end

function test_matpower_import_uses_bus_type_for_regulation()::Bool
  mpc = (
    name = "case_bus_type_controls_regulation",
    baseMVA = 100.0,
    bus = [
      1 3 0.0 0.0 0.0 0.0 1 1.0 0.0 110.0 1 1.1 0.9
      2 1 0.0 0.0 0.0 0.0 1 1.0 0.0 110.0 1 1.1 0.9
      3 2 0.0 0.0 0.0 0.0 1 1.02 0.0 110.0 1 1.1 0.9
    ],
    gen = [
      1 10.0 0.0 999.0 -999.0 1.0 100.0 1 999.0 0.0 0 0 0 0 0 0 0 0 0 0 0
      2 30.0 0.0 999.0 -999.0 1.01 100.0 1 999.0 0.0 0 0 0 0 0 0 0 0 0 0 0
      3 40.0 0.0 999.0 -999.0 1.02 100.0 1 999.0 0.0 0 0 0 0 0 0 0 0 0 0 0
    ],
    branch = [
      1 2 0.01 0.05 0.0 9999.0 0.0 0.0 0.0 0.0 1 -60.0 60.0
      2 3 0.01 0.05 0.0 9999.0 0.0 0.0 0.0 0.0 1 -60.0 60.0
    ],
    gencost = nothing,
    bus_name = nothing,
  )

  net = nothing
  @test_logs (:info, r"MATPOWER import: PQ generator limits mapped to constant P\(U\)/Q\(U\) controllers") match_mode = :any begin
    net = Sparlectra.createNetFromMatPowerCase(mpc = mpc, log = false, flatstart = false)
  end

  b2 = geNetBusIdx(net = net, busName = "2")
  b3 = geNetBusIdx(net = net, busName = "3")
  b2_prosumers = getBusProsumers(net, b2)
  b3_prosumers = getBusProsumers(net, b3)
  b2_gen = only(filter(isGenerator, b2_prosumers))
  b3_gen = only(filter(isGenerator, b3_prosumers))

  return getNodeType(net.nodeVec[b2]) == Sparlectra.PQ && getNodeType(net.nodeVec[b3]) == Sparlectra.PV && !isnothing(b2_gen.puController) && !isnothing(b2_gen.quController) && isnothing(b3_gen.puController) && isnothing(b3_gen.quController)
end

function test_prosumer_aggregation_preserves_bus_types_and_injections()::Bool
  net = Net(name = "aggregation_regression", baseMVA = 100.0)
  addBus!(net = net, busName = "Slack", vn_kV = 110.0)
  addBus!(net = net, busName = "PV", vn_kV = 110.0)
  addBus!(net = net, busName = "PQ", vn_kV = 110.0)

  addPIModelACLine!(net = net, fromBus = "Slack", toBus = "PV", r_pu = 0.01, x_pu = 0.10, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "PV", toBus = "PQ", r_pu = 0.01, x_pu = 0.10, b_pu = 0.0, status = 1)

  addProsumer!(net = net, busName = "Slack", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "Slack", defer_bus_type_refresh = true)
  addProsumer!(net = net, busName = "PV", type = "SYNCHRONOUSMACHINE", p = 50.0, q = 5.0, vm_pu = 1.01, defer_bus_type_refresh = true)
  addProsumer!(net = net, busName = "PV", type = "ENERGYCONSUMER", p = 10.0, q = 3.0, defer_bus_type_refresh = true)
  addProsumer!(net = net, busName = "PQ", type = "ENERGYCONSUMER", p = 20.0, q = 7.0, defer_bus_type_refresh = true)
  addProsumer!(net = net, busName = "PQ", type = "GENERATOR", p = 2.0, q = 1.0, defer_bus_type_refresh = true)

  refreshBusTypesFromProsumers!(net)
  S = buildComplexSVec(net)

  slack = geNetBusIdx(net = net, busName = "Slack")
  pv = geNetBusIdx(net = net, busName = "PV")
  pq = geNetBusIdx(net = net, busName = "PQ")

  return getNodeType(net.nodeVec[slack]) == Sparlectra.Slack && getNodeType(net.nodeVec[pv]) == Sparlectra.PV && getNodeType(net.nodeVec[pq]) == Sparlectra.PQ && isapprox(S[pv], ComplexF64(0.40, 0.02); atol = 1e-12) && isapprox(S[pq], ComplexF64(-0.18, -0.06); atol = 1e-12)
end

function test_rectangular_autodamp_backtracks_oversized_step()::Bool
  Y = ComplexF64[0.0 - 10.0im 0.0 + 10.0im; 0.0 + 10.0im 0.0 - 10.0im]
  V = ComplexF64[1.0+0.0im, 1.0+0.0im]
  S = ComplexF64[0.0+0.0im, -1.0-0.2im]
  bus_types = [:Slack, :PQ]
  Vset = [1.0, 1.0]
  F0 = Sparlectra.mismatch_rectangular(Y, V, S, bus_types, Vset, 1)

  oversized_delta = [0.0, -1.0]
  alpha, Vtrial, trial_mismatch = Sparlectra.choose_rectangular_autodamp(Y, V, S, oversized_delta, F0; slack_idx = 1, damp = 1.0, autodamp_min = 1e-3, bus_types = bus_types, Vset = Vset)
  full_step_mismatch = maximum(abs.(Sparlectra.mismatch_rectangular(Y, ComplexF64[1.0+0.0im, 1.0-1.0im], S, bus_types, Vset, 1)))

  return alpha < 1.0 && trial_mismatch < maximum(abs.(F0)) && trial_mismatch < full_step_mismatch && Vtrial[1] == V[1]
end

# The polar update applies the tangent step dV as magnitude times
# 1 + Re(dV/V) and angle plus Im(dV/V): a rotation keeps the magnitude,
# where the linear update V + dV inflates it by sqrt(1 + b^2) (the flat-start
# divergence on case9241pegase). A zero-voltage bus has no
# polar form and takes the linear step.
function test_rectangular_polar_update_rotates_without_inflation()::Bool
  V = ComplexF64[1.0+0.0im, 1.0+0.0im, 0.0+0.0im]
  # bus 2: pure rotation by 0.5 rad; bus 3 (zero voltage): linear step
  delta = [0.0, 0.1, 0.5, 0.2]
  non_slack = [2, 3]
  Vp = Sparlectra._apply_rectangular_delta(V, delta, 1, non_slack, 1.0; update = :polar)
  Vr = Sparlectra._apply_rectangular_delta(V, delta, 1, non_slack, 1.0)
  ok = isapprox(abs(Vp[2]), 1.0; atol = 1e-14) && isapprox(angle(Vp[2]), 0.5; atol = 1e-14)
  ok &= isapprox(abs(Vr[2]), sqrt(1.25); atol = 1e-14)
  ok &= Vp[3] == ComplexF64(0.1, 0.2) && Vp[1] == V[1]
  # a scaling step (dV parallel to V) is the same in both updates
  Vs = Sparlectra._apply_rectangular_delta(V, [-0.2, 0.0, 0.0, 0.0], 1, non_slack, 1.0; update = :polar)
  ok &= isapprox(Vs[2], 0.8 + 0.0im; atol = 1e-14)
  # half a step scales the polar update like the linear one (autodamp alpha)
  Vh = Sparlectra._apply_rectangular_delta(V, delta, 1, non_slack, 0.5; update = :polar)
  ok &= isapprox(angle(Vh[2]), 0.25; atol = 1e-14) && isapprox(abs(Vh[2]), 1.0; atol = 1e-14)
  return ok
end

# Both updates are Newton on the same equations and reach the same state; the
# polar update from the flat start of the shipped 60-bus case converges, and
# the unsupported value is refused.
# Power mode (0.30.1): the Ybus, the symbolic analysis and the work arrays
# stay on the network between solves. The solution with the switch equals
# the solution without it to 1e-12 pu on two shipped cases (warm solve on
# the same net, Ybus and context reused), a branch taken out of service
# after a power-mode solve gives the solution a fresh net gives (the Ybus
# fingerprint moves, the context re-analyses), and the cache does not
# survive a deepcopy (native factorization memory is never shared).
function test_power_mode_same_solution_and_topology_change()::Bool
  _state(net) = [n._vm_pu * cis(deg2rad(n._va_deg)) for n in net.nodeVec]
  for f in ("sp_case118.m", "sp_case300.m")
    path = joinpath(dirname(@__DIR__), "data", "mpower", f)
    ref_net = Sparlectra.createNetFromMatPowerFile(filename = path)
    it_ref, erg_ref = runpf!(ref_net, 30, 1e-10, 0; qlimits_enabled = false)
    erg_ref == 0 || return false
    net = Sparlectra.createNetFromMatPowerFile(filename = path)
    runpf!(net, 30, 1e-10, 0; qlimits_enabled = false, power_mode = true)
    # the warm solve: restart from the file columns on the same net
    fresh = Sparlectra.createNetFromMatPowerFile(filename = path)
    for (a, b) in zip(net.nodeVec, fresh.nodeVec)
      a._vm_pu = b._vm_pu
      a._va_deg = b._va_deg
    end
    it_pm, erg_pm = runpf!(net, 30, 1e-10, 0; qlimits_enabled = false, power_mode = true)
    erg_pm == 0 && it_pm == it_ref || return false
    maximum(abs.(_state(net) .- _state(ref_net))) <= 1e-12 || return false
    cache = net._power_cache
    cache isa Sparlectra.PowerModeCache && cache.ybus_reuse_count >= 1 && cache.context_reuse_count >= 1 || return false
    # a copy starts without the cache
    copied = deepcopy(net)
    copied._power_cache isa Sparlectra.PowerModeCache && copied._power_cache.ybus === nothing || return false
    # topology change: the first in-service non-radial branch out, solve again
    br = net.branchVec[1]
    Sparlectra.setBranchStatus!(br, false)
    Sparlectra.setBranchStatus!(ref_net.branchVec[1], false)
    it_t, erg_t = runpf!(net, 30, 1e-10, 0; qlimits_enabled = false, power_mode = true)
    it_r, erg_r = runpf!(ref_net, 30, 1e-10, 0; qlimits_enabled = false)
    erg_t == erg_r || return false
    if erg_r == 0
      maximum(abs.(_state(net) .- _state(ref_net))) <= 1e-10 || return false
      cache.ybus_reuse_count == 1 || return false   # the fingerprint moved, the Ybus was rebuilt
    end
  end
  # 0.30.2: N-1 with power mode equals N-1 without it, outage by outage.
  # The kept Jacobian assembly replays recorded positions; it was reused
  # when the bus types matched the previous solve's START types although
  # (a) that solve had ended with other types and (b) an outage had rebuilt
  # the Ybus with another pattern. With the warm active set the types match
  # from outage to outage, which made 55 of 186 sp_case118 outages differ
  # (one diverged) before the fix.
  path118 = joinpath(dirname(@__DIR__), "data", "mpower", "sp_case118.m")
  cases = Sparlectra.generateN1Branches(Sparlectra.createNetFromMatPowerFile(filename = path118))
  plain = Sparlectra.runContingencies!(Sparlectra.createNetFromMatPowerFile(filename = path118), cases; parallel_enabled = false, warm_active_set = true, warm_note = Ref(""))
  pm = Sparlectra.runContingencies!(Sparlectra.createNetFromMatPowerFile(filename = path118), cases; parallel_enabled = false, warm_active_set = true, power_mode = true, warm_note = Ref(""))
  all(i -> plain[i].converged == pm[i].converged && plain[i].iterations == pm[i].iterations, eachindex(plain)) || return false
  # #457: the worker reset must not wipe the worker's power-mode cache when
  # the template carries none (a rescued or flat base case)
  template = Sparlectra.createNetFromMatPowerFile(filename = path118)
  template._power_cache = nothing
  work = deepcopy(template)
  runpf!(work, 30, 1e-8, 0; power_mode = true)
  kept = work._power_cache
  kept isa Sparlectra.PowerModeCache || return false
  Sparlectra._reset_scenario_worker!(work, template)
  work._power_cache === kept || return false
  runpf!(work, 30, 1e-8, 0; power_mode = true)
  kept.context_reuse_count >= 1 || return false
  return true
end

# Dishonest Newton (power_flow.jacobian_reuse, 0.30.2): the next step solves
# with the previous factorization while the mismatch falls by the minimum
# reduction. The converged state equals the honest solution within the
# solver tolerance on a shipped case, with and without power mode, and the
# reuse actually happens (reused steps > 0, refactorisations counted). With
# the Q-limit active set on, the PQ set is the same as with honest Newton.
# On the `umfpack` backend there is no factorization to reuse: no error, no
# reused step, and the run log says so.
function test_jacobian_reuse_same_solution()::Bool
  _state(net) = [n._vm_pu * cis(deg2rad(n._va_deg)) for n in net.nodeVec]
  path(f) = joinpath(dirname(@__DIR__), "data", "mpower", f)
  for pm in (false, true)
    honest = Sparlectra.createNetFromMatPowerFile(filename = path("sp_case1354.m"))
    _, erg_h = runpf!(honest, 30, 1e-8, 0; qlimits_enabled = false, opt_flatstart = true, power_mode = pm)
    net = Sparlectra.createNetFromMatPowerFile(filename = path("sp_case1354.m"))
    prof = Dict{Symbol,Any}()
    _, erg = runpf!(net, 30, 1e-8, 0; qlimits_enabled = false, opt_flatstart = true, power_mode = pm, jacobian_reuse = true, performance_profile = prof)
    erg_h == 0 && erg == 0 || return false
    maximum(abs.(_state(net) .- _state(honest))) <= 1e-8 || return false
    prof[:jacobian_reuse_steps] > 0 && prof[:jacobian_reuse_refactorisations] >= 2 || return false
    st = Sparlectra.rectangular_pf_status(net)
    st.jacobian_reuse && st.jacobian_reuse_steps == prof[:jacobian_reuse_steps] || return false
  end
  # active set: same PQ set and the same state as honest Newton
  honest = Sparlectra.createNetFromMatPowerFile(filename = path("sp_case300.m"))
  _, erg_h = runpf!(honest, 30, 1e-8, 0; opt_flatstart = true)
  net = Sparlectra.createNetFromMatPowerFile(filename = path("sp_case300.m"))
  _, erg = runpf!(net, 30, 1e-8, 0; opt_flatstart = true, jacobian_reuse = true)
  erg_h == 0 && erg == 0 || return false
  [Sparlectra.getNodeType(n) for n in net.nodeVec] == [Sparlectra.getNodeType(n) for n in honest.nodeVec] || return false
  maximum(abs.(_state(net) .- _state(honest))) <= 1e-8 || return false
  # umfpack keeps no factorization: no reuse, not silent
  net = Sparlectra.createNetFromMatPowerFile(filename = path("sp_case118.m"))
  prof = Dict{Symbol,Any}()
  log = mktemp() do tmp, io
    redirect_stdout(io) do
      runpf_rectangular!(net; maxiter = 30, tol = 1e-8, verbose = 1, qlimits_enabled = false, opt_flatstart = true, linear_solver = :umfpack, jacobian_reuse = true, performance_profile = prof)
    end
    close(io)
    read(tmp, String)
  end
  occursin("Jacobian reuse ignored: linear_solver umfpack keeps no factorisation", log) || return false
  prof[:jacobian_reuse_steps] == 0 || return false
  return true
end

function test_rectangular_polar_update_reaches_the_same_solution()::Bool
  net_r = load_fixture_net("sp_case60")
  net_p = load_fixture_net("sp_case60")
  it_r, erg_r = runpf!(net_r, 30, 1e-10, 0; method = :rectangular, opt_flatstart = true, qlimits_enabled = false, newton_update = :rectangular)
  it_p, erg_p = runpf!(net_p, 30, 1e-10, 0; method = :rectangular, opt_flatstart = true, qlimits_enabled = false, newton_update = :polar)
  erg_r == 0 && erg_p == 0 || return false
  it_p <= 30 || return false
  for (a, b) in zip(net_r.nodeVec, net_p.nodeVec)
    isapprox(a._vm_pu, b._vm_pu; atol = 1e-7) || return false
    isapprox(a._va_deg, b._va_deg; atol = 1e-6) || return false
  end
  refused = try
    runpf!(load_fixture_net("sp_case60"), 30, 1e-10, 0; method = :rectangular, newton_update = :cartesian)
    false
  catch err
    err isa ArgumentError && occursin("newton_update", sprint(showerror, err))
  end
  return refused
end

function test_rectangular_nonfinite_mismatch_diagnostics_use_finite_history()::Bool
  Y = ComplexF64[0.0 - 10.0im 0.0 + 10.0im; 0.0 + 10.0im 0.0 - 10.0im]
  S = ComplexF64[0.0+0.0im, -1.0-0.2im]
  bus_types = [:Slack, :PQ]
  Vset = [1.0, 1.0]
  V1 = ComplexF64[1.0+0.0im, 1.0+0.0im]
  V2 = ComplexF64[1.0+0.0im, 0.9-0.1im]
  diagnostics = Sparlectra._rectangular_mismatch_diagnostics(
    Y,
    ComplexF64[1.0+0.0im, NaN+NaN*im],
    S,
    bus_types,
    Vset,
    1,
    0.0;
    history = [239.0, 1609.0, 5.55e290, NaN],
    step_diagnostics = [(alpha = 0.05, trial_mismatch = 1609.0, accepted_improvement = false) for _ = 1:8],
    best_finite_iteration = 1,
    best_finite_voltage = V1,
    last_finite_iteration = 3,
    last_finite_voltage = V2,
    first_nonfinite_iteration = 4,
    first_nonfinite_voltage = ComplexF64[1.0+0.0im, NaN+NaN*im],
  )

  # the final status reads the same history: best mismatch without the NaN,
  # reason nr_nonfinite
  status =
    Sparlectra._build_rectangular_final_status(Net(name = "diagnostic_status", baseMVA = 100.0), false, false, false, :nr_mismatch_not_converged, nothing, nothing, 0.0, [239.0, 1609.0, 5.55e290, NaN], 0, 0, 0, Int[], Sparlectra._wrong_branch_not_checked_result(), :warn, false, :disabled, NamedTuple()).status

  return diagnostics.first_nonfinite_iteration == 4 &&
         diagnostics.last_finite_mismatch == 5.55e290 &&
         diagnostics.best_mismatch == 239.0 &&
         diagnostics.mismatch_history_trend == :diverging_to_nonfinite &&
         length(diagnostics.top_mismatch_rows_by_iteration) == 3 &&
         diagnostics.top_mismatch_rows_by_iteration[1].label == :best_finite_iteration &&
         diagnostics.top_mismatch_rows_by_iteration[2].label == :last_finite_iteration &&
         diagnostics.autodamp_failure === true &&
         status.reason == :nr_nonfinite &&
         status.status == :nr_nonfinite &&
         status.best_mismatch == 239.0
end

function test_rectangular_start_projection_dc_candidate_matches_rundcpf()::Bool
  # Regression (grid-bench study 7.1, case6495rte): the projection's DC
  # candidate solved imag(Ybus) with the bus injections only and missed the
  # phase-shifter injections. With the net's DC model (assemble_dc_bbus, the
  # model rundcpf! solves) the candidate angles equal rundcpf!'s; the path
  # without the injection vector is off by exactly Bred \ Pbusinj.
  shift_deg = 10.0
  net = Net(name = "dc3_shift", baseMVA = 100.0)
  addBus!(net = net, busName = "B1", vn_kV = 110.0)
  addBus!(net = net, busName = "B2", vn_kV = 110.0)
  addBus!(net = net, busName = "B3", vn_kV = 110.0)
  addPIModelACLine!(net = net, fromBus = "B1", toBus = "B2", r_pu = 0.02, x_pu = 0.08, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "B1", toBus = "B3", r_pu = 0.03, x_pu = 0.12, b_pu = 0.0, status = 1)
  addPIModelTrafo!(net = net, fromBus = "B2", toBus = "B3", r_pu = 0.02, x_pu = 0.10, b_pu = 0.0, status = 1, ratio = 1.0, shift_deg = shift_deg)
  addProsumer!(net = net, busName = "B1", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.02, va_deg = 0.0, referencePri = "B1")
  addProsumer!(net = net, busName = "B2", type = "GENERATOR", p = 25.0, q = 0.0, vm_pu = 1.01)
  addProsumer!(net = net, busName = "B2", type = "ENERGYCONSUMER", p = 20.0, q = 10.0)
  addProsumer!(net = net, busName = "B3", type = "ENERGYCONSUMER", p = 45.0, q = 15.0)
  refreshBusTypesFromProsumers!(net)

  bus_types, Vset, slack_idx = Sparlectra.extract_bus_types_and_vset(net)
  n = length(bus_types)
  Y = createYBUS(net = net, sparse = true, printYBUS = false)
  S = buildComplexSVec(net)
  Vraw = ComplexF64[Vset[k] + 0.0im for k in 1:n]
  model = Sparlectra._projection_dc_model(net)
  model isa Sparlectra._ProjectionDcModel || return false
  any(!iszero, model.Pbusinj) || return false

  Vdc_model = Sparlectra._dc_angle_start_rectangular(Y, Vraw, S, bus_types, Vset, slack_idx; dc_angle_limit_deg = 60.0, dc_model = model)
  Vdc_ybus = Sparlectra._dc_angle_start_rectangular(Y, Vraw, S, bus_types, Vset, slack_idx; dc_angle_limit_deg = 60.0)
  no_injection = Sparlectra._ProjectionDcModel(model.Bbus, zeros(Float64, n))
  Vdc_noinj = Sparlectra._dc_angle_start_rectangular(Y, Vraw, S, bus_types, Vset, slack_idx; dc_angle_limit_deg = 60.0, dc_model = no_injection)

  report = rundcpf!(net)
  report.metadata.converged || return false
  theta_dc = deg2rad.([row.va_deg for row in report.nodes])
  theta_model = angle.(Vdc_model)
  theta_ybus = angle.(Vdc_ybus)
  theta_noinj = angle.(Vdc_noinj)

  # The shift term: the reduced Bbus applied to the shift injections.
  non_slack = Sparlectra.non_slack_indices(n, slack_idx)
  pos = Sparlectra.build_pos_map(non_slack, n)
  Bred = Sparlectra.reduce_susceptance_to_nonslack(model.Bbus, non_slack, pos, slack_idx, length(non_slack); value_of = x -> -x)
  shift_term = Bred \ model.Pbusinj[non_slack]

  # The projection itself with a requested DC angle start hands the model
  # candidate through and names the model in its summary.
  profile = Dict{Symbol,Any}(:enabled => true)
  Vproj = Sparlectra.project_rectangular_start(Y, Vraw, S, bus_types, Vset, slack_idx; enabled = true, try_dc_start = true, try_blend_scan = false, dc_angle_limit_deg = 60.0, requested_angle_mode = :dc, performance_profile = profile, dc_model = model)
  summary = profile[:start_projection_summary]

  return maximum(abs.(theta_model .- theta_dc)) <= 1e-12 &&
         maximum(abs.(theta_ybus .- theta_dc)) > 1e-4 &&
         maximum(abs.((theta_noinj[non_slack] .- theta_model[non_slack]) .- shift_term)) <= 1e-12 &&
         maximum(abs.(shift_term)) > 1e-4 &&
         Vproj == Vdc_model &&
         summary.selected === :dc_start &&
         summary.dc_model === :dc_bbus
end

function test_rectangular_start_projection_keeps_raw_without_finite_improvement()::Bool
  Ydense = ComplexF64[0.0 - 10.0im 0.0 + 10.0im; 0.0 + 10.0im 0.0 - 10.0im]
  Y = sparse(Ydense)
  Vraw = ComplexF64[1.0+0.0im, 1.0+0.0im]
  bus_types = [:Slack, :PQ]
  Vset = [1.0, 1.0]

  # an optional DC candidate (no start mode requested) replaces the seed
  # when it lowers the mismatch
  Sload = ComplexF64[0.0+0.0im, -1.0-0.2im]
  raw_mismatch = maximum(abs.(Sparlectra.mismatch_rectangular(Y, Vraw, Sload, bus_types, Vset, 1)))
  Voptional = Sparlectra.project_rectangular_start(Y, Vraw, Sload, bus_types, Vset, 1; enabled = true, try_dc_start = true, try_blend_scan = true, blend_lambdas = [0.25, 0.5, 0.75], dc_angle_limit_deg = 60.0, performance_profile = Dict{Symbol,Any}(:enabled => true))
  optional_mismatch = maximum(abs.(Sparlectra.mismatch_rectangular(Y, Voptional, Sload, bus_types, Vset, 1)))
  Vdc = Sparlectra._dc_angle_start_rectangular(Y, Vraw, Sload, bus_types, Vset, 1; dc_angle_limit_deg = 60.0)
  optional_dc_ok = Voptional[1] == Vraw[1] && optional_mismatch < raw_mismatch && Vdc == Voptional

  unmeasured_profile = Dict{Symbol,Any}(:enabled => true)
  Vunmeasured = Sparlectra.project_rectangular_start(Y, Vraw, ComplexF64[0.0+0.0im, -1.0-0.2im], bus_types, Vset, 1; enabled = true, try_dc_start = true, try_blend_scan = false, measure_candidates = false, performance_profile = unmeasured_profile)

  nonfinite_profile = Dict{Symbol,Any}(:enabled => true)
  Vnonfinite = Sparlectra.project_rectangular_start(Y, Vraw, ComplexF64[0.0+0.0im, NaN+0.0im], bus_types, Vset, 1; enabled = true, try_dc_start = true, try_blend_scan = false, measure_candidates = true, performance_profile = nonfinite_profile)

  # a requested DC start replaces the seed when its residual 2-norm is at
  # least 10 percent below the seed's (#465): with load at bus 2 the flat
  # seed has 1.02, the DC start 0.25
  Srequested = ComplexF64[0.0+0.0im, -1.0-0.2im]
  requested_dc_profile = Dict{Symbol,Any}(:enabled => true)
  Vrequested_dc = Sparlectra.project_rectangular_start(
    Y,
    Vraw,
    Srequested,
    bus_types,
    Vset,
    1;
    enabled = true,
    try_dc_start = true,
    try_blend_scan = false,
    measure_candidates = true,
    requested_angle_mode = :dc,
    requested_voltage_mode = :profile_blend,
    performance_profile = requested_dc_profile,
  )

  # without injections the flat seed is the exact solution and the DC start
  # equals it: no candidate is 10 percent better, the seed stays (#465)
  tie_profile = Dict{Symbol,Any}(:enabled => true)
  Sparlectra.project_rectangular_start(Y, Vraw, ComplexF64[0.0+0.0im, 0.0+0.0im], bus_types, Vset, 1; enabled = true, try_dc_start = true, try_blend_scan = false, measure_candidates = true, requested_angle_mode = :dc, requested_voltage_mode = :profile_blend, performance_profile = tie_profile)
  tie_summary = tie_profile[:start_projection_summary]

  invalid_dc_profile = Dict{Symbol,Any}(:enabled => true)
  Vinvalid_dc = Sparlectra.project_rectangular_start(
    Y,
    Vraw,
    ComplexF64[0.0+0.0im, NaN+0.0im],
    bus_types,
    Vset,
    1;
    enabled = true,
    try_dc_start = false,
    try_blend_scan = false,
    measure_candidates = true,
    requested_angle_mode = :dc,
    requested_voltage_mode = :profile_blend,
    performance_profile = invalid_dc_profile,
  )

  # Regression: a seed that is already a solved state (mismatch below the
  # DC start's) is kept although DC angles were requested; the projection
  # used to replace the file start of case6495rte (0.03 pu) with the DC
  # start (575 pu), which did not converge.
  Ssolved = ComplexF64[0.0+0.0im, -1.0-0.2im]
  Vsolved = copy(Vraw)
  for _ in 1:20
    Vsolved[2] = 1.0 + 0.0im + conj(Ssolved[2] / Vsolved[2]) / (0.0 - 10.0im)
  end
  solved_seed_profile = Dict{Symbol,Any}(:enabled => true)
  Vsolved_seed = Sparlectra.project_rectangular_start(
    Y,
    Vsolved,
    Ssolved,
    bus_types,
    Vset,
    1;
    enabled = true,
    try_dc_start = true,
    try_blend_scan = true,
    measure_candidates = true,
    requested_angle_mode = :dc,
    requested_voltage_mode = :profile_blend,
    performance_profile = solved_seed_profile,
  )
  solved_seed_summary = solved_seed_profile[:start_projection_summary]

  requested_summary = requested_dc_profile[:start_projection_summary]
  invalid_summary = invalid_dc_profile[:start_projection_summary]

  # Regression: optional DC candidates still require measured improvement, but an explicitly
  # requested DC angle start is applied whenever it passes the hard finite-voltage guard.
  return optional_dc_ok &&
         Vunmeasured == Vraw &&
         unmeasured_profile[:start_projection_summary].selected === :raw &&
         unmeasured_profile[:start_projection_summary].reason === :candidate_mismatch_not_measured &&
         ismissing(unmeasured_profile[:start_projection_summary].best_mismatch) &&
         Vnonfinite == Vraw &&
         nonfinite_profile[:start_projection_summary].selected === :raw &&
         nonfinite_profile[:start_projection_summary].reason === :no_finite_improvement &&
         ismissing(nonfinite_profile[:start_projection_summary].best_mismatch) &&
         Vrequested_dc == Sparlectra._dc_angle_start_rectangular(Y, Vraw, Srequested, bus_types, Vset, 1; dc_angle_limit_deg = 60.0) &&
         tie_summary.selected === :raw &&
         tie_summary.reason === :raw_seed_closer_than_dc_start &&
         requested_summary.selected === :dc_start &&
         requested_summary.reason === :requested_dc_angle_start &&
         requested_summary.dc_angle_start_built === true &&
         requested_summary.dc_angle_start_valid === true &&
         requested_summary.dc_angle_start_applied === true &&
         requested_summary.fallback_to_raw === false &&
         requested_summary.requested_dc_worse_than_raw === (requested_summary.dc_mismatch > requested_summary.raw_mismatch) &&
         Vsolved_seed == Vsolved &&
         solved_seed_summary.selected === :raw &&
         solved_seed_summary.reason === :raw_seed_closer_than_dc_start &&
         solved_seed_summary.fallback_to_raw === true &&
         solved_seed_summary.dc_angle_start_applied === false &&
         solved_seed_summary.raw_mismatch < solved_seed_summary.dc_mismatch &&
         Vinvalid_dc == Vraw &&
         invalid_summary.selected === :explicit_fallback_raw &&
         invalid_summary.reason === :invalid_dc_angle_start &&
         invalid_summary.dc_angle_start_built === true &&
         invalid_summary.dc_angle_start_valid === false &&
         invalid_summary.dc_angle_start_applied === false &&
         invalid_summary.fallback_to_raw === true &&
         invalid_summary.fallback_reason === :invalid_dc_angle_start
end

# The ratio-profile start candidate (power_flow.start_mode.ratio_profile):
# on a flat start the projection also measures the flat profile with every
# PQ magnitude scaled by the off-nominal ratios on its path from the
# reference. Pinned: the exact two-bus case (the candidate IS the solution),
# no candidate without factors, a requested DC start keeping precedence,
# the breadth-first factors on a shipped case, and the forwarding chain
# from the keyword and from the configuration down to the projection.
function test_rectangular_start_projection_ratio_profile()::Bool
  # one branch with the tap t on its from side and no load: the solution is
  # V2 = V1 / t, and the flat start puts the whole ratio offset across the
  # series impedance
  t = 0.95
  y = 1.0 / (0.01 + 0.1im)
  Y = sparse(ComplexF64[y/t^2 -y/t; -y/t y])
  Vraw = ComplexF64[1.0, 1.0]
  S = zeros(ComplexF64, 2)
  bus_types = [:Slack, :PQ]
  Vset = [1.0, 1.0]
  factors = [1.0, 1.0 / t]
  run_projection(; kw...) = (prof = Dict{Symbol,Any}(:enabled => true); V = Sparlectra.project_rectangular_start(Y, Vraw, S, bus_types, Vset, 1; enabled = true, try_dc_start = true, try_blend_scan = true, performance_profile = prof, kw...); (V, prof[:start_projection_summary]))
  V1, s1 = run_projection(; ratio_profile_factors = factors)
  exact = s1.selected === :ratio_profile && s1.reason === :finite_improvement && V1[2] == ComplexF64(1.0 / t, 0.0) && V1[1] == Vraw[1] && s1.ratio_profile_built === true && s1.ratio_profile_mismatch < 1e-12 && s1.candidates == 6
  # no factors: nothing is built and the selection is the one without the
  # candidate (the flat seed against DC and blends, all equal here)
  V0, s0 = run_projection()
  no_candidate = s0.ratio_profile_built === false && ismissing(s0.ratio_profile_mismatch) && s0.selected === :raw && V0 == Vraw && s0.candidates == 5
  # a requested DC angle start keeps precedence over the candidate, as it
  # does over the blends; the candidate is still measured and reported
  Vdc, sdc = run_projection(; ratio_profile_factors = factors, requested_angle_mode = :dc)
  dc_first = sdc.selected !== :ratio_profile && sdc.ratio_profile_built === true && Vdc != V1

  # breadth-first factors on a shipped case: the shipped transformers sit at
  # their nominal ratio, so one of them is moved off nominal (0.95) here.
  # The reference has 1.0, every other bus is its search parent's factor
  # times one branch step, and the level behind that transformer differs
  # from 1.0
  net = load_fixture_net("sp_case14")
  itr = findfirst(br -> Sparlectra.getNodeVn(net.nodeVec[Int(br.fromBus)]) != Sparlectra.getNodeVn(net.nodeVec[Int(br.toBus)]), net.branchVec)
  itr === nothing && return false
  net.branchVec[itr].ratio = 0.95
  net.branchVec[itr].tap_ratio = 0.95
  slack = findfirst(nd -> Sparlectra.getNodeType(nd) == Sparlectra.Slack, net.nodeVec)
  f = Sparlectra._ratio_profile_factors(net, slack)
  steps(k) = [f[Int(br.fromBus)] * (1.0 / abs(Sparlectra.calcBranchRatio(br))) for br in net.branchVec if br.status == 1 && Int(br.toBus) == k] ∪ [f[Int(br.toBus)] * abs(Sparlectra.calcBranchRatio(br)) for br in net.branchVec if br.status == 1 && Int(br.fromBus) == k]
  bfs = f[slack] == 1.0 && all(k -> k == slack || f[k] in steps(k), eachindex(f)) && any(x -> x != 1.0, f) && all(x -> isfinite(x) && x > 0.0, f)

  # forwarding: the solver builds the candidate only when the projection
  # runs on a flat seed with the keyword on; without the projection nothing
  # is measured at all (a run-level flat start); the configuration default
  # reaches it
  built(net0; kw...) = (prof = Dict{Symbol,Any}(:enabled => true); Sparlectra.runpf_rectangular!(net0; maxiter = 30, tol = 1e-8, damp = 1.0, qlimits_enabled = false, performance_profile = prof, kw...); s = get(prof, :start_projection_summary, nothing); s === nothing ? nothing : s.ratio_profile_built)
  chain = built(load_fixture_net("sp_case14"); opt_flatstart = true, start_projection = true, start_projection_ratio_profile = true) === true &&
          built(load_fixture_net("sp_case14"); opt_flatstart = true, start_projection = false, start_projection_ratio_profile = true) === nothing &&
          built(load_fixture_net("sp_case14"); opt_flatstart = true, start_projection = true) === false &&
          built(load_fixture_net("sp_case14"); opt_flatstart = false, start_projection = true, start_projection_ratio_profile = true) === false
  cfg_prof = Dict{Symbol,Any}(:enabled => true)
  runpf!(load_fixture_net("sp_case14"); config = Sparlectra.PowerFlowConfig(start_mode = Sparlectra.StartModeConfig(flatstart = true, start_projection = true)), performance_profile = cfg_prof)
  config_default = Sparlectra.StartModeConfig().ratio_profile === true && cfg_prof[:start_projection_summary].ratio_profile_built === true

  # auxiliary-bus flat start (the same branch steps, limited to auxiliary
  # buses, applied without the projection): a net without auxiliary buses
  # keeps its flat start vector bit-identical (sp_case14 with the moved
  # transformer from above)
  Vflat, _ = Sparlectra.initialVrect(net; flatstart = true)
  Vkept = copy(Vflat)
  no_aux = Sparlectra._aux_bus_flat_start!(Vkept, net, slack) == 0 && isequal(Vkept, Vflat)
  # a star point between the reference (stiff winding, ratio 0.95 on the
  # star side) and a load bus (weaker winding, ratio 1): the star point starts
  # at the current-free value of the stiff winding, the visible buses stay
  # flat, and the solver applies it without the projection: the initial
  # mismatch is the 0.05 pu offset across the weak winding (about 0.5 pu),
  # where the plain flat start would put it across the stiff one (over 40 pu)
  star_net() = begin
    sn = Net(name = "aux_star", baseMVA = 100.0)
    addBus!(net = sn, busName = "R", vn_kV = 110.0)
    addBus!(net = sn, busName = "AUX3WT_T", vn_kV = 110.0, isAux = true)
    addBus!(net = sn, busName = "L", vn_kV = 110.0)
    addProsumer!(net = sn, busName = "R", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "R")
    addPIModelTrafo!(net = sn, fromBus = "AUX3WT_T", toBus = "R", r_pu = 0.0005, x_pu = 0.001, b_pu = 0.0, status = 1, ratio = 0.95, ratedS = 100.0)
    addPIModelTrafo!(net = sn, fromBus = "AUX3WT_T", toBus = "L", r_pu = 0.01, x_pu = 0.1, b_pu = 0.0, status = 1, ratio = 1.0, ratedS = 100.0)
    addProsumer!(net = sn, busName = "L", type = "ENERGYCONSUMER", p = 5.0, q = 1.0)
    sn
  end
  sn = star_net()
  ka, kr, kl = sn.busDict["AUX3WT_T"], sn.busDict["R"], sn.busDict["L"]
  Vs, sslack = Sparlectra.initialVrect(sn; flatstart = true)
  rule = Sparlectra._aux_bus_flat_start!(Vs, sn, sslack) == 1 && isapprox(Vs[ka], ComplexF64(abs(Sparlectra.calcBranchRatio(sn.branchVec[1])), 0.0); atol = 1e-14) && Vs[kr] == 1.0 + 0.0im && Vs[kl] == 1.0 + 0.0im
  sn2 = star_net()
  _, erg = Sparlectra.runpf_rectangular!(sn2; maxiter = 30, tol = 1e-8, damp = 1.0, newton_update = :polar, opt_flatstart = true, start_projection = false, qlimits_enabled = false, verbose = 0)
  solver_chain = erg == 0 && Sparlectra.rectangular_pf_status(sn2).initial_mismatch < 1.0

  return exact && no_candidate && dc_first && bfs && chain && config_default && no_aux && rule && solver_chain
end

# A run-level flat start is a flat start (before, a flat start selected in the Web UI on sp_case118 still ran from the DC
# start of the projection the configuration template switches on). With
# power_flow.flatstart the run ignores the stored VM/VA and runs no
# projection, no DC seed and no ratio-profile candidate; only auxiliary buses
# (three-winding star points) start at their ratio, every angle at the slack
# angle; with power_flow.start_mode.rotate_behind_pst buses behind a phase
# shifter at the slack angle plus the shift. A start from stored voltages
# keeps the projection.
function test_run_level_flat_start_is_flat()::Bool
  mp_dir = joinpath(dirname(@__DIR__), "data", "mpower")
  # the template's start settings: projection on, DC angle and blended
  # voltage modes, ratio profile on (its default)
  template_start(flat; kw...) = Sparlectra.StartModeConfig(; flatstart = flat, start_projection = true, try_dc_start = true, try_blend_scan = true, angle_mode = :dc, voltage_mode = :profile_blend, ratio_profile = true, kw...)
  function start_of(start_mode; casefile = "sp_case118.m", net = nothing)
    cfg = Sparlectra.SparlectraConfig(powerflow = Sparlectra.PowerFlowConfig(start_mode = start_mode), output = Sparlectra.OutputConfig(logfile_results = :off, startup_latency_hint = false))
    prof = Dict{Symbol,Any}(:enabled => true, :capture_initial_residual_rows => true)
    r = net === nothing ? run_sparlectra(casefile = casefile, path = mp_dir, config = cfg, performance_profile = prof) : run_sparlectra(net = net, config = cfg, performance_profile = prof)
    # the first solve's rows only (a Q-limit round would append more)
    rows = prof[:initial_residual_rows][1:length(r.net.nodeVec)]
    return (r = r, rows = rows, summary = get(prof, :start_projection_summary, nothing), timings = get(prof, :timings, Dict{Symbol,Any}()))
  end
  # every non-auxiliary bus flat: PQ at 1.0 pu, PV and slack at their
  # setpoints, every angle at the slack angle (the reference angle from the
  # case, 0.99 deg on micro_grid_be); with `turned` (rotate_behind_pst) plus
  # the phase shifts of the flat-start angles
  shift_deg(net) = (th = Sparlectra._flat_start_shift_angles(net, findfirst(nd -> Sparlectra.getNodeType(nd) == Sparlectra.Slack, net.nodeVec)); th === nothing ? zeros(length(net.nodeVec)) : rad2deg.(th))
  flat_except_aux(x, net; turned = false) = (vset = Sparlectra._bus_voltage_setpoints_from_prosumers(net); sh = turned ? shift_deg(net) : zeros(length(net.nodeVec));
                                             va_s = x.rows[findfirst(r -> r.bus_type === :Slack, x.rows)].va_deg_start;
                                             all(k -> Sparlectra._is_aux_node(net.nodeVec[k]) || ((x.rows[k].bus_type === :Slack || isapprox(x.rows[k].va_deg_start, va_s + sh[k]; atol = 1e-12)) && x.rows[k].vm_pu_start == (x.rows[k].bus_type === :PQ ? 1.0 : vset[k])), eachindex(x.rows)))

  # sp_case118 from a flat start although the file carries a solved state:
  # no projection, no candidate, no projection phase in the timings
  flat = start_of(template_start(true))
  no_projection = flat.r.final_converged && flat.summary === nothing && !haskey(flat.timings, :start_projection) && flat_except_aux(flat, flat.r.net)
  # the same start as a bare flat start with the ratio profile off
  bare = start_of(Sparlectra.StartModeConfig(flatstart = true, start_projection = false, ratio_profile = false))
  same_start = bare.summary === nothing && [(x.vm_pu_start, x.va_deg_start) for x in bare.rows] == [(x.vm_pu_start, x.va_deg_start) for x in flat.rows]
  # the DC seed is off under a flat start as well
  seeded = start_of(template_start(true; dc_seed_unconditional = true))
  no_seed = Sparlectra.dc_pf_status(seeded.r.net) === nothing && [(x.vm_pu_start, x.va_deg_start) for x in seeded.rows] == [(x.vm_pu_start, x.va_deg_start) for x in flat.rows]
  # a start from the stored voltages keeps the projection
  stored = start_of(template_start(false))
  projected = stored.summary !== nothing && stored.summary.dc_angle_start_built === true && haskey(stored.timings, :start_projection)

  # case14 (shipped sp_case14, its two transformers moved off nominal so the
  # ratio profile would scale the 20 kV level): every bus starts flat. The
  # control: the factors the ratio profile would apply differ from 1.0 on
  # PQ buses, so the check would fail if the candidate were applied.
  c14 = load_fixture_net("sp_case14")
  for (br, t) in zip([br for br in c14.branchVec if Sparlectra.getNodeVn(c14.nodeVec[Int(br.fromBus)]) != Sparlectra.getNodeVn(c14.nodeVec[Int(br.toBus)])], (0.95, 0.93))
    br.ratio = t
    br.tap_ratio = t
  end
  c14_unsolved = deepcopy(c14)  # the solve below changes the net
  slack14 = findfirst(nd -> Sparlectra.getNodeType(nd) == Sparlectra.Slack, c14.nodeVec)
  would_scale = any(k -> Sparlectra.getNodeType(c14.nodeVec[k]) == Sparlectra.PQ && Sparlectra._ratio_profile_factors(c14, slack14)[k] != 1.0, eachindex(c14.nodeVec))
  c14_flat = start_of(template_start(true); net = c14)
  case14_flat = would_scale && c14_flat.r.final_converged && c14_flat.summary === nothing && !any(nd -> Sparlectra._is_aux_node(nd), c14_flat.r.net.nodeVec) && flat_except_aux(c14_flat, c14_flat.r.net)
  # the rescue strategies of a flat-start run that project from the flat
  # profile (settled_qlimits, dc_seed) see the ratio profile switched off as
  # well. Control: settled_qlimits with the candidate forced on keeps it and
  # raises PQ buses, so the check would fail if the run config kept the key.
  function rescue_start(sm, ratio)
    prof = Dict{Symbol,Any}(:enabled => true, :capture_initial_residual_rows => true)
    net = deepcopy(c14_unsolved)
    Sparlectra.runpf_rectangular!(net; maxiter = 30, tol = 1e-8, damp = 1.0, qlimits_enabled = false, performance_profile = prof,
      opt_flatstart = sm.flatstart, start_projection = sm.start_projection, start_projection_try_dc_start = sm.try_dc_start,
      start_projection_try_blend_scan = sm.try_blend_scan, start_projection_requested_angle_mode = sm.angle_mode,
      start_projection_requested_voltage_mode = sm.voltage_mode, start_projection_ratio_profile = ratio,
      start_projection_accept_unmeasured_dc_start = sm.accept_unmeasured_dc_start)
    rows = prof[:initial_residual_rows][1:length(net.nodeVec)]
    return (built = prof[:start_projection_summary].ratio_profile_built, pq_flat = all(r -> r.bus_type !== :PQ || r.vm_pu_start == 1.0, rows))
  end
  run_cfg, _ = Sparlectra._flatstart_profile_config(Sparlectra.SparlectraConfig(powerflow = Sparlectra.PowerFlowConfig(start_mode = template_start(true))))
  flat_rescue = [pf.start_mode for (name, pf) in Sparlectra._rescue_config_variants(run_cfg.powerflow) if name in (:settled_qlimits, :dc_seed)]
  rescue_flat = length(flat_rescue) == 2 && all(sm -> sm.flatstart && sm.start_projection && !sm.ratio_profile, flat_rescue) &&
                all(sm -> (x = rescue_start(sm, sm.ratio_profile); !x.built && x.pq_flat), flat_rescue) &&
                !rescue_start(flat_rescue[1], true).pq_flat

  # auxiliary buses keep the ratio start (micro_grid_be, a three-winding
  # transformer whose star point has a 0.95 winding): the star point starts
  # at 0.95 pu, every other bus flat
  dir = joinpath(@__DIR__, "fixtures", "powsybl", "micro_grid_be.powsybl")
  be, _ = Sparlectra.build_net_from_powsybl(Sparlectra.read_powsybl_bundle(dir), Sparlectra.PowsyblAdapterOptions(remote_regulation = :remote))
  aux = [k for k in eachindex(be.nodeVec) if Sparlectra._is_aux_node(be.nodeVec[k])]
  be_flat = start_of(template_start(true); net = be)
  # the slack angle of micro_grid_be (0.99 deg) is not 0, so the angle check
  # is not only a check of zeros. Its one phase shifter (-2.15 deg) lies in a
  # mesh, so rotate_behind_pst turns no bus: the start with the switch on is
  # the start without it
  be_slack_va = be_flat.rows[findfirst(r -> r.bus_type === :Slack, be_flat.rows)].va_deg_start
  be2, _ = Sparlectra.build_net_from_powsybl(Sparlectra.read_powsybl_bundle(dir), Sparlectra.PowsyblAdapterOptions(remote_regulation = :remote))
  be_turn = start_of(template_start(true; rotate_behind_pst = true); net = be2)
  be_shifted = abs(be_slack_va) > 0.5 && any(br -> br.status == 1 && angle(Sparlectra.calcBranchRatio(br)) != 0.0, be2.branchVec) &&
               flat_except_aux(be_turn, be_turn.r.net; turned = true) && [(x.vm_pu_start, x.va_deg_start) for x in be_turn.rows] == [(x.vm_pu_start, x.va_deg_start) for x in be_flat.rows]
  aux_ratio = any(k -> isapprox(be_flat.rows[k].vm_pu_start, 0.95; atol = 1e-12), aux) && flat_except_aux(be_flat, be_flat.r.net) && be_shifted && be_flat.summary === nothing

  return no_projection && same_start && no_seed && projected && case14_flat && rescue_flat && aux_ratio
end

# The flat start behind phase shifters (power_flow.start_mode.rotate_behind_pst,
# keyword start_rotate_behind_pst, default false; 0.60.0). Net: slack S
# (110 kV, 1.02 pu) - line 0.01 + j0.08 - F (110 kV) - transformer F-T
# (ratio 1.05, z 0.002 + j0.06, shift 30 or 90 deg) - T (20 kV, load 30 MW,
# 10 Mvar). With every angle at 0 (switch off, the default) the polar Newton
# converged at 30 and 90 deg to the low-voltage root (0.044 pu at T) and
# reported erg 0.
# Hand calculation of both roots: refer the line to the T side of the ideal
# transformer, E = 1.02 / 1.05, Z = zl / 1.05^2 + zt = R + jX; a shifter
# only turns the angle at T, so with P = 0.3, Q = 0.1 the magnitude solves
# |U|^4 + (2(PR + QX) - E^2) |U|^2 + (P^2 + Q^2) |Z|^2 = 0,
# roots 0.95319002 pu and 0.04413162 pu. With the switch on the flat start
# must give T the angle -shift, and the solve the high root with the
# iterations and the angle at T of the solve at shift 0 (the shift only turns
# T). The same for the stored start of a hand-built net (addBus! defaults
# 1 pu, 0 deg), which holds no angle; a stored start with an angle is left as
# it is. With the switch off T starts at 0.
# Mesh: S - line - F - line - T and a 30 deg shifter S-T parallel to that
# path, plus a 20 deg shifter T-R. T is reached without crossing a shifter
# and keeps 0 (in a mesh the parallel path fixes its angle, the turn would
# put the whole shift across the lines); R is reached across one shifter
# (T-R, -20) and not across two (S-T-R, -50), which a search for the fewest
# branches would pick.
# Slack angle 10 deg: every bus starts at 10 deg, T with the switch on at
# 10 - shift, and the solution is the one at slack angle 0 turned by 10 deg.
# A star point (three-winding) with a shift on its stiffest winding starts
# at the turned value V_aux = t V_R with the switch on and passes the shift on
# to the bus behind. The model-only MATPOWER export writes the default flat
# start (star point at its ratio, no turn); read back without flat start,
# that start holds no angle, so a reader with the switch on turns it
# (PV bus G behind a second 90 deg shifter at F included) and reaches the
# solution of the flat start. With slack angle 10 deg the file holds every
# bus at 10 deg; the reader takes that as a start without angles as well and
# starts T at 10 - 90 = -80 deg.
# Wrong-branch hint (0.60.0): a voltage finding after a start without angles
# on a net with a shifter and the switch off names
# power_flow.start_mode.rotate_behind_pst = true in a warning and in the
# result summary line; no hint with the switch on or without a shifter.
function test_flat_start_carries_phase_shift()::Bool
  function shifted_net(shift; va_slack = 0.0)
    net = Net(name = "flat_shift", baseMVA = 100.0)
    addBus!(net = net, busName = "S", vn_kV = 110.0)
    addBus!(net = net, busName = "F", vn_kV = 110.0)
    addBus!(net = net, busName = "T", vn_kV = 20.0)
    addProsumer!(net = net, busName = "S", type = "EXTERNALNETWORKINJECTION", referencePri = "S", vm_pu = 1.02, va_deg = va_slack)
    addProsumer!(net = net, busName = "T", type = "ENERGYCONSUMER", p = 30.0, q = 10.0)
    addPIModelACLine!(net = net, fromBus = "S", toBus = "F", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
    addPIModelTrafo!(net = net, fromBus = "F", toBus = "T", r_pu = 0.002, x_pu = 0.06, b_pu = 0.0, ratio = 1.05, shift_deg = shift, status = 1)
    # the reference angle as an imported state sets it (addProsumer! keeps
    # the bus angle it finds)
    va_slack == 0.0 || setVmVa!(node = net.nodeVec[net.busDict["S"]], vm_pu = 1.02, va_deg = va_slack)
    return net
  end
  E = 1.02 / 1.05
  Z = (0.01 + 0.08im) / 1.05^2 + (0.002 + 0.06im)
  b = 2 * (0.3 * real(Z) + 0.1 * imag(Z)) - E^2
  c = (0.3^2 + 0.1^2) * abs2(Z)
  u_high = sqrt((-b + sqrt(b^2 - 4c)) / 2)
  ref = shifted_net(0.0)
  ite0, erg0 = runpf!(ref, 40, 1e-12, 0; method = :rectangular, opt_flatstart = true)
  va0 = ref.nodeVec[ref.busDict["T"]]._va_deg
  ok = erg0 == 0 && Sparlectra.StartModeConfig().rotate_behind_pst === false
  for shift in (30.0, 90.0)
    net = shifted_net(shift)
    t = net.busDict["T"]
    # switch off: T starts at 0 like every other bus
    Voff, _ = Sparlectra.initialVrect(net; flatstart = true)
    ok &= angle(Voff[t]) == 0.0 && angle(Voff[net.busDict["F"]]) == 0.0
    V0, _ = Sparlectra.initialVrect(net; flatstart = true, rotate_behind_pst = true)
    ok &= isapprox(rad2deg(angle(V0[t])), -shift; atol = 1e-12) && angle(V0[net.busDict["F"]]) == 0.0 && abs(V0[t]) == 1.0
    ite, erg = runpf!(net, 40, 1e-12, 0; method = :rectangular, opt_flatstart = true, start_rotate_behind_pst = true)
    nd = net.nodeVec[t]
    ok &= erg == 0 && ite == ite0 && isapprox(nd._vm_pu, u_high; atol = 1e-9) && isapprox(nd._va_deg, va0 - shift; atol = 1e-9) && Sparlectra.rectangular_pf_status(net).branch_quality_status === :ok
    # stored start without angles: turned like the flat start, high root
    sn3 = shifted_net(shift)
    Vst, _ = Sparlectra.initialVrect(sn3; flatstart = false, rotate_behind_pst = true)
    ok &= isapprox(rad2deg(angle(Vst[t])), -shift; atol = 1e-12)
    _, erg_s = runpf!(sn3, 40, 1e-12, 0; method = :rectangular, start_rotate_behind_pst = true)
    nds = sn3.nodeVec[t]
    ok &= erg_s == 0 && isapprox(nds._vm_pu, u_high; atol = 1e-9) && isapprox(nds._va_deg, va0 - shift; atol = 1e-9) && Sparlectra.rectangular_pf_status(sn3).branch_quality_status === :ok
    # model-only export: the default flat start, T not listed
    fs = Sparlectra._matpower_flat_start(shifted_net(shift))
    ok &= !haskey(fs, t) && !haskey(fs, net.busDict["F"])
  end
  # wrong-branch hint: with the switch off the flat start at 30 deg ends on
  # the low root (voltage_collapse, erg 0), and the warning and the result
  # summary line name the switch; with the switch on (high root) no hint.
  # A net without shifter that ends on the low root from a stored start
  # without angles (T at 0.05 pu, 0 deg) gets no hint either.
  hn = shifted_net(30.0)
  _, erg_h = @test_logs (:warn, r"power_flow\.start_mode\.rotate_behind_pst = true") match_mode = :any runpf!(hn, 40, 1e-12, 0; method = :rectangular, opt_flatstart = true)
  st_h = Sparlectra.rectangular_pf_status(hn)
  ok &= erg_h == 0 && st_h.wrong_branch_reason === :voltage_collapse && occursin("rotate_behind_pst = true", sprint(io -> Sparlectra._print_wrong_branch_summary_line(io, hn)))
  hint_logs, _ = Test.collect_test_logs() do
    runpf!(shifted_net(30.0), 40, 1e-12, 0; method = :rectangular, opt_flatstart = true, start_rotate_behind_pst = true)
  end
  ok &= !any(l -> occursin("rotate_behind_pst", string(l.message)), hint_logs)
  # the no-angle test reads the vector Newton starts from: a flat start that
  # the DC start projection replaces (T at -32.4 deg after the solve, high
  # root) does not count as a start without angles
  hd = shifted_net(30.0)
  runpf!(hd, 40, 1e-12, 0; method = :rectangular, opt_flatstart = true, start_projection = true, start_projection_requested_angle_mode = :dc, start_projection_accept_unmeasured_dc_start = true)
  ok &= st_h.start_without_angles_behind_pst === true && Sparlectra.rectangular_pf_status(hd).start_without_angles_behind_pst === false
  ns = shifted_net(0.0)
  setVmVa!(node = ns.nodeVec[ns.busDict["T"]], vm_pu = 0.05, va_deg = 0.0)
  _, erg_ns = runpf!(ns, 40, 1e-12, 0; method = :rectangular, opt_flatstart = false)
  st_ns = Sparlectra.rectangular_pf_status(ns)
  ok &= erg_ns == 0 && st_ns.wrong_branch_reason === :voltage_collapse && isempty(Sparlectra._wrong_branch_pst_start_hint(st_ns))
  # the configuration key reaches the solver (YAML-style dict, run_sparlectra)
  cfg_rot = SparlectraConfig(Dict("power_flow" => Dict("max_iter" => 40, "tol" => 1e-12, "flatstart" => true, "start_mode" => Dict("rotate_behind_pst" => true)),
                                  "output" => Dict("logfile_results" => "off", "console_summary" => false)))
  r_rot = run_sparlectra(net = shifted_net(30.0), config = cfg_rot)
  ok &= cfg_rot.powerflow.start_mode.rotate_behind_pst === true && r_rot.final_converged && isapprox(r_rot.net.nodeVec[r_rot.net.busDict["T"]]._vm_pu, u_high; atol = 1e-9)
  # slack angle 10 deg: the profile counts from it
  sa = shifted_net(30.0; va_slack = 10.0)
  ta, fa = sa.busDict["T"], sa.busDict["F"]
  Va_off, _ = Sparlectra.initialVrect(sa; flatstart = true)
  Va_on, _ = Sparlectra.initialVrect(sa; flatstart = true, rotate_behind_pst = true)
  ok &= isapprox(rad2deg(angle(Va_off[ta])), 10.0; atol = 1e-12) && isapprox(rad2deg(angle(Va_off[fa])), 10.0; atol = 1e-12) &&
        isapprox(rad2deg(angle(Va_on[ta])), -20.0; atol = 1e-12) && isapprox(rad2deg(angle(Va_on[fa])), 10.0; atol = 1e-12)
  _, erg_a = runpf!(sa, 40, 1e-12, 0; method = :rectangular, opt_flatstart = true, start_rotate_behind_pst = true)
  ok &= erg_a == 0 && isapprox(sa.nodeVec[ta]._vm_pu, u_high; atol = 1e-9) && isapprox(sa.nodeVec[ta]._va_deg, va0 - 30.0 + 10.0; atol = 1e-9)
  # a stored start with an angle is not touched
  kept = shifted_net(30.0)
  kept.nodeVec[kept.busDict["T"]]._va_deg = -5.0
  Vk, _ = Sparlectra.initialVrect(kept; flatstart = false, rotate_behind_pst = true)
  ok &= isapprox(rad2deg(angle(Vk[kept.busDict["T"]])), -5.0; atol = 1e-12) && angle(Vk[kept.busDict["F"]]) == 0.0
  # mesh: the shifter S-T parallel to the lines S-F-T keeps T at 0, R behind
  # the shifter T-R gets -20 (fewest shifters), not -50 (fewest branches)
  mesh = Net(name = "flat_shift_mesh", baseMVA = 100.0)
  foreach(bn -> addBus!(net = mesh, busName = bn, vn_kV = 110.0), ("S", "F", "T", "R"))
  addProsumer!(net = mesh, busName = "S", type = "EXTERNALNETWORKINJECTION", referencePri = "S", vm_pu = 1.0, va_deg = 0.0)
  addProsumer!(net = mesh, busName = "R", type = "ENERGYCONSUMER", p = 10.0, q = 2.0)
  addPIModelACLine!(net = mesh, fromBus = "S", toBus = "F", r_pu = 0.01, x_pu = 0.1, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = mesh, fromBus = "F", toBus = "T", r_pu = 0.01, x_pu = 0.1, b_pu = 0.0, status = 1)
  addPIModelTrafo!(net = mesh, fromBus = "S", toBus = "T", r_pu = 0.002, x_pu = 0.06, b_pu = 0.0, ratio = 1.0, shift_deg = 30.0, status = 1)
  addPIModelTrafo!(net = mesh, fromBus = "T", toBus = "R", r_pu = 0.002, x_pu = 0.06, b_pu = 0.0, ratio = 1.0, shift_deg = 20.0, status = 1)
  Vm_on, _ = Sparlectra.initialVrect(mesh; flatstart = true, rotate_behind_pst = true)
  ok &= angle(Vm_on[mesh.busDict["T"]]) == 0.0 && angle(Vm_on[mesh.busDict["F"]]) == 0.0 && isapprox(rad2deg(angle(Vm_on[mesh.busDict["R"]])), -20.0; atol = 1e-12)
  _, erg_m = runpf!(mesh, 40, 1e-12, 0; method = :rectangular, opt_flatstart = true, start_rotate_behind_pst = true)
  ok &= erg_m == 0 && minimum(nd -> nd._vm_pu, mesh.nodeVec) > 0.9
  # star point: R - (AUX-R, ratio 0.95, shift 30, stiff) - AUX - (AUX-L, ratio 1) - L
  sn = Net(name = "aux_star_shift", baseMVA = 100.0)
  addBus!(net = sn, busName = "R", vn_kV = 110.0)
  addBus!(net = sn, busName = "AUX3WT_T", vn_kV = 110.0, isAux = true)
  addBus!(net = sn, busName = "L", vn_kV = 110.0)
  addProsumer!(net = sn, busName = "R", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "R")
  addPIModelTrafo!(net = sn, fromBus = "AUX3WT_T", toBus = "R", r_pu = 0.0005, x_pu = 0.001, b_pu = 0.0, status = 1, ratio = 0.95, shift_deg = 30.0, ratedS = 100.0)
  addPIModelTrafo!(net = sn, fromBus = "AUX3WT_T", toBus = "L", r_pu = 0.01, x_pu = 0.1, b_pu = 0.0, status = 1, ratio = 1.0, ratedS = 100.0)
  addProsumer!(net = sn, busName = "L", type = "ENERGYCONSUMER", p = 5.0, q = 1.0)
  Vs, slack = Sparlectra.initialVrect(sn; flatstart = true, rotate_behind_pst = true)
  Sparlectra._aux_bus_flat_start!(Vs, sn, slack; rotate_behind_pst = true)
  ka, kl = sn.busDict["AUX3WT_T"], sn.busDict["L"]
  ok &= isapprox(Vs[ka], 0.95 * cis(deg2rad(30.0)); atol = 1e-14) && isapprox(Vs[kl], cis(deg2rad(30.0)); atol = 1e-14)
  Vsd, _ = Sparlectra.initialVrect(sn; flatstart = true)
  Sparlectra._aux_bus_flat_start!(Vsd, sn, slack)
  ok &= isapprox(Vsd[ka], ComplexF64(0.95, 0.0); atol = 1e-14) && angle(Vsd[ka]) == 0.0 && Vsd[kl] == 1.0 + 0.0im
  fsn = Sparlectra._matpower_flat_start(sn)
  ok &= isapprox(fsn[ka][1], 0.95; atol = 1e-14) && fsn[ka][2] == 0.0 && !haskey(fsn, kl)
  # PV bus G (1.0 pu, 10 MW) behind a second 90 degree shifter at F: the
  # model-only export of the unsolved net, read back without flat start and
  # solved with the switch on, gives the solution of the flat start
  pv_net() = (n = shifted_net(90.0); addBus!(net = n, busName = "G", vn_kV = 20.0);
              addProsumer!(net = n, busName = "G", type = "SYNCHRONOUSMACHINE", p = 10.0, vm_pu = 1.0);
              addPIModelTrafo!(net = n, fromBus = "F", toBus = "G", r_pu = 0.002, x_pu = 0.06, b_pu = 0.0, ratio = 1.0, shift_deg = 90.0, status = 1); n)
  pv_ref = pv_net()
  _, erg_pv = runpf!(pv_ref, 40, 1e-12, 0; method = :rectangular, opt_flatstart = true, start_rotate_behind_pst = true)
  ok &= erg_pv == 0
  mktempdir() do dir
    f = joinpath(dir, "pv_shift.m")
    Sparlectra.writeMatpowerCasefile(pv_net(), f; write_solution = false)
    back = Sparlectra.createNetFromMatPowerFile(filename = f, flatstart = false)
    kt, kg = (findfirst(nd -> nd.busIdx == pv_ref.busDict[b], back.nodeVec) for b in ("T", "G"))
    Vb, _ = Sparlectra.initialVrect(back; rotate_behind_pst = true)
    ok &= isapprox(rad2deg(angle(Vb[kt])), -90.0; atol = 1e-9) && isapprox(rad2deg(angle(Vb[kg])), -90.0; atol = 1e-9)
    _, erg_b = runpf!(back, 40, 1e-12, 0; method = :rectangular, start_rotate_behind_pst = true)
    for (k, name) in ((kt, "T"), (kg, "G"))
      nr = pv_ref.nodeVec[pv_ref.busDict[name]]
      ok &= erg_b == 0 && isapprox(back.nodeVec[k]._vm_pu, nr._vm_pu; atol = 1e-9) && isapprox(back.nodeVec[k]._va_deg, nr._va_deg; atol = 1e-7)
    end
    # slack angle 10 deg: the export writes every bus at 10 deg, which the
    # reader with the switch on still takes as a start without angles
    f10 = joinpath(dir, "shift_va10.m")
    net10 = shifted_net(90.0; va_slack = 10.0)
    Sparlectra.writeMatpowerCasefile(net10, f10; write_solution = false)
    back10 = Sparlectra.createNetFromMatPowerFile(filename = f10, flatstart = false)
    kt10, kf10 = (findfirst(nd -> nd.busIdx == net10.busDict[b], back10.nodeVec) for b in ("T", "F"))
    V10, _ = Sparlectra.initialVrect(back10; rotate_behind_pst = true)
    ok &= isapprox(rad2deg(angle(V10[kt10])), -80.0; atol = 1e-9) && isapprox(rad2deg(angle(V10[kf10])), 10.0; atol = 1e-9)
    _, erg10 = runpf!(back10, 40, 1e-12, 0; method = :rectangular, start_rotate_behind_pst = true)
    ok &= erg10 == 0 && isapprox(back10.nodeVec[kt10]._vm_pu, u_high; atol = 1e-9)
  end
  return ok
end

function test_matpower_vmva_selfcheck_ignores_slack_pq_spec()::Bool
  mpc = Sparlectra.MatpowerIO.MatpowerCase(
    "case_slack_spec_ignored",
    100.0,
    [
      1 3 0.0 0.0 0.0 0.0 1 1.0 0.0 110.0 1 1.1 0.9
      2 1 0.0 0.0 0.0 0.0 1 1.0 0.0 110.0 1 1.1 0.9
    ],
    [
    # very large slack Pg/Qg that should NOT be enforced by PF equations
    1 500.0 300.0 999.0 -999.0 1.0 100.0 1 999.0 0.0 0 0 0 0 0 0 0 0 0 0 0],
    [
    # no branch -> Scalc = 0 at all buses
    # keep out-of-service branch to satisfy matrix width expectations
    1 2 0.01 0.05 0.0 9999.0 0.0 0.0 0.0 0.0 0 -60.0 60.0],
    nothing,
    nothing,
  )

  stats = Sparlectra.MatpowerIO.vmva_power_mismatch_stats(mpc)
  slack_ok = get(stats, :ok, false) && isapprox(get(stats, :max_p_mis_pu, NaN), 0.0; atol = 1e-12) && isapprox(get(stats, :max_q_mis_pu, NaN), 0.0; atol = 1e-12)

  # a non-contiguous bus number (9001) is mapped, one P and one Q row counted
  mpc_noncontig = Sparlectra.MatpowerIO.MatpowerCase(
    "case_noncontig",
    100.0,
    [
      1 3 0.0 0.0 0.0 0.0 1 1.00 0.0 110.0 1 1.1 0.9
      9001 1 0.0 0.0 0.0 0.0 1 0.99 -3.0 110.0 1 1.1 0.9
    ],
    [1 0.0 0.0 999.0 -999.0 1.0 100.0 1 999.0 0.0 0 0 0 0 0 0 0 0 0 0 0],
    [1 9001 0.01 0.05 0.0 9999.0 0.0 0.0 0.0 0.0 1 -60.0 60.0],
    nothing,
    nothing,
  )
  stats_nc = Sparlectra.MatpowerIO.vmva_power_mismatch_stats(mpc_noncontig)
  noncontig_ok = get(stats_nc, :ok, false) && get(stats_nc, :n_p, 0) == 1 && get(stats_nc, :n_q, 0) == 1 && isfinite(get(stats_nc, :max_p_mis_pu, NaN)) && isfinite(get(stats_nc, :max_q_mis_pu, NaN))
  return slack_ok && noncontig_ok
end

function test_matpower_build_ybus_returns_sparse_expected_values()::Bool
  base = 100.0
  bus = [
    1 3 0.0 0.0 10.0 -5.0 1 1.0 0.0 110.0 1 1.1 0.9
    2 1 0.0 0.0 0.0 0.0 1 1.0 0.0 110.0 1 1.1 0.9
  ]
  branch = [1 2 0.01 0.05 0.02 9999.0 0.0 0.0 2.0 30.0 1 -60.0 60.0]
  ybus = Sparlectra.MatpowerIO.build_ybus_matpower(bus, branch, base; matpower_shift_unit = :deg, matpower_shift_sign = 1.0, matpower_ratio = :normal)
  y = inv(0.01 + 0.05im)
  ysh = 0.01im
  tap = 2.0 * cis(pi / 6)
  expected11 = (10.0 - 5.0im) / base + (y + ysh) / (tap * conj(tap))
  expected22 = y + ysh
  expected12 = -y / conj(tap)
  expected21 = -y / tap

  return SparseArrays.issparse(ybus) && ybus isa SparseArrays.SparseMatrixCSC{ComplexF64,Int} && isapprox(ybus[1, 1], expected11; atol = 1e-12) && isapprox(ybus[2, 2], expected22; atol = 1e-12) && isapprox(ybus[1, 2], expected12; atol = 1e-12) && isapprox(ybus[2, 1], expected21; atol = 1e-12)
end

function test_matpower_compare_vmva_wraps_angle_differences()::Bool
  mpc = Sparlectra.MatpowerIO.MatpowerCase(
    "case_angle_wrap_compare",
    100.0,
    [
      1 3 0.0 0.0 0.0 0.0 1 1.00 0.0 110.0 1 1.1 0.9
      2 1 0.0 0.0 0.0 0.0 1 0.99 179.0 110.0 1 1.1 0.9
    ],
    [1 0.0 0.0 999.0 -999.0 1.0 100.0 1 999.0 0.0 0 0 0 0 0 0 0 0 0 0 0],
    [1 2 0.01 0.05 0.0 9999.0 0.0 0.0 0.0 0.0 0 -60.0 60.0],
    nothing,
    nothing,
  )
  net = Sparlectra.createNetFromMatPowerCase(mpc = mpc, log = false, flatstart = false)
  net.nodeVec[1]._vm_pu = 1.0
  net.nodeVec[1]._va_deg = 0.0
  net.nodeVec[2]._vm_pu = 0.99
  net.nodeVec[2]._va_deg = -179.0

  ok, stats = redirect_stdout(devnull) do
    Sparlectra.MatpowerIO.compare_vm_va(net, mpc; show_diff = false, tol_vm = 1e-12, tol_va = 2.1)
  end
  return ok && isapprox(get(stats, :max_dva, NaN), 2.0; atol = 1e-12)
end

function _synthetic_pv_vg_mismatch_case(; gen_rows = nothing)
  gens = gen_rows === nothing ? [
    1 0.0 0.0 999.0 -999.0 1.0 100.0 1 999.0 0.0 zeros(1, 11)...
    2 50.0 0.0 999.0 -999.0 1.04 100.0 1 999.0 0.0 zeros(1, 11)...
  ] : gen_rows
  return Sparlectra.MatpowerIO.MatpowerCase(
    "case_pv_vg_mismatch",
    100.0,
    [
      1 3 0.0 0.0 0.0 0.0 1 1.00 0.0 110.0 1 1.1 0.9
      2 2 0.0 0.0 0.0 0.0 1 1.02 0.0 110.0 1 1.1 0.9
      3 1 20.0 5.0 0.0 0.0 1 0.98 -2.0 110.0 1 1.1 0.9
    ],
    gens,
    [
      1 2 0.01 0.05 0.0 9999.0 0.0 0.0 0.0 0.0 1 -60.0 60.0
      2 3 0.01 0.05 0.0 9999.0 0.0 0.0 0.0 0.0 1 -60.0 60.0
    ],
    nothing,
    nothing,
  )
end

function test_matpower_pv_voltage_source_and_compare_modes()::Bool
  mpc = _synthetic_pv_vg_mismatch_case()
  net_gen = Sparlectra.createNetFromMatPowerCase(mpc = mpc, matpower_pv_voltage_source = :gen_vg)
  net_bus = Sparlectra.createNetFromMatPowerCase(mpc = mpc, matpower_pv_voltage_source = :bus_vm)

  pv_idx = geNetBusIdx(net = net_gen, busName = "2")
  source_ok = isapprox(net_gen.nodeVec[pv_idx]._vm_pu, 1.04; atol = 1e-12) && isapprox(net_bus.nodeVec[pv_idx]._vm_pu, 1.02; atol = 1e-12)

  net_gen.nodeVec[pv_idx]._vm_pu = 1.04
  bus_ok, bus_stats = redirect_stdout(devnull) do
    Sparlectra.MatpowerIO.compare_vm_va(net_gen, mpc; compare_voltage_reference = :bus_vm, tol_vm = 0.01, tol_va = 180.0)
  end
  imported_ok, imported_stats = redirect_stdout(devnull) do
    Sparlectra.MatpowerIO.compare_vm_va(net_gen, mpc; compare_voltage_reference = :imported_setpoint, matpower_pv_voltage_source = :gen_vg, tol_vm = 0.01, tol_va = 180.0)
  end

  rows = Sparlectra.MatpowerIO.pv_voltage_reference_rows(mpc; net = net_gen, matpower_pv_voltage_source = :gen_vg)
  diagnostic_ok = any(row -> row.busI == 2 && isapprox(row.dvg_bus, 0.02; atol = 1e-12) && isapprox(row.dvm_vset, 0.0; atol = 1e-12), rows)

  # hybrid compare: a PV bus switched to PQ by its Q limit is compared
  # against the bus Vm of the file (0.01 pu off at 1.03 against 1.02)
  net_q = Sparlectra.createNetFromMatPowerCase(mpc = mpc, matpower_pv_voltage_source = :gen_vg)
  Sparlectra.setNodeType!(net_q.nodeVec[pv_idx], "PQ")
  net_q.qLimitEvents[pv_idx] = :max
  net_q.nodeVec[pv_idx]._vm_pu = 1.03
  net_q.nodeVec[pv_idx]._va_deg = 0.0
  hybrid_ok, hybrid_stats = redirect_stdout(devnull) do
    Sparlectra.MatpowerIO.compare_vm_va(net_q, mpc; show_diff = true, compare_voltage_reference = :hybrid, matpower_pv_voltage_source = :gen_vg, tol_vm = 1e-12, tol_va = 180.0)
  end
  counts = get(hybrid_stats, :vm_ref_kind_counts, Dict{Symbol,Int}())
  final_pq_summary = get(get(hybrid_stats, :vm_ref_kind_summary, Dict{Symbol,Any}()), :final_pq_after_qlimit, (count = 0, max_dvm = NaN, max_dva = NaN))
  hybrid_rows_ok =
    hybrid_ok &&
    get(hybrid_stats, :compare_status, :fail) == :warn &&
    isapprox(get(hybrid_stats, :max_dvm, NaN), 0.01; atol = 1e-12) &&
    isapprox(final_pq_summary.max_dvm, 0.01; atol = 1e-12) &&
    get(counts, :final_pq_after_qlimit, 0) == 1 &&
    get(counts, :pq_bus_vm, 0) == 1 &&
    get(counts, :ref_slack_imported_setpoint, 0) == 1 &&
    get(counts, :active_pv_imported_setpoint, 0) == 0

  # several generators on one bus: the first online Vg wins; an offline
  # generator gives no setpoint and the bus Vm of the file is used
  mpc_multi = _synthetic_pv_vg_mismatch_case(gen_rows = [
    1 0.0 0.0 999.0 -999.0 1.0 100.0 1 999.0 0.0 zeros(1, 11)...
    2 25.0 0.0 999.0 -999.0 1.04 100.0 1 999.0 0.0 zeros(1, 11)...
    2 25.0 0.0 999.0 -999.0 1.05 100.0 1 999.0 0.0 zeros(1, 11)...
  ])
  net_multi = redirect_stderr(devnull) do
    Sparlectra.createNetFromMatPowerCase(mpc = mpc_multi, matpower_pv_voltage_source = :gen_vg, matpower_pv_voltage_mismatch_tol_pu = 1e-4)
  end
  first_online_ok = isapprox(net_multi.nodeVec[geNetBusIdx(net = net_multi, busName = "2")]._vm_pu, 1.04; atol = 1e-12)
  mpc_fallback = _synthetic_pv_vg_mismatch_case(gen_rows = [
    1 0.0 0.0 999.0 -999.0 1.0 100.0 1 999.0 0.0 zeros(1, 11)...
    2 50.0 0.0 999.0 -999.0 1.04 100.0 0 999.0 0.0 zeros(1, 11)...
  ])
  net_fallback = redirect_stderr(devnull) do
    Sparlectra.createNetFromMatPowerCase(mpc = mpc_fallback, matpower_pv_voltage_source = :gen_vg)
  end
  fallback_ok = isapprox(net_fallback.nodeVec[geNetBusIdx(net = net_fallback, busName = "2")]._vm_pu, 1.02; atol = 1e-12)

  return source_ok && !bus_ok && get(bus_stats, :max_dvm, 0.0) > 0.019 && imported_ok && get(imported_stats, :max_dvm, 1.0) < 0.011 && diagnostic_ok && hybrid_rows_ok && first_online_ok && fallback_ok
end

function testISOBusses()
  Sbase_MVA = 1000.0
  netName = "isobus"
  net = Net(name = netName, baseMVA = Sbase_MVA)
  @debug "Creating $netName test network"
  addBus!(net = net, busName = "B1", vn_kV = 220.0)
  addBus!(net = net, busName = "B2", vn_kV = 220.0)
  addBus!(net = net, busName = "B3", vn_kV = 220.0)
  addBus!(net = net, busName = "B4", vn_kV = 220.0)
  addBus!(net = net, busName = "B5", vn_kV = 220.0)

  add_test_line_pu_c!(net = net, fromBus = "B1", toBus = "B2", length_km = 100.0, r_pu = 0.0653, x_pu = 0.398, c_nf_per_km = 9.08, tanδ = 0.0, ratedS = 250.0, status = 1)
  add_test_line_pu_c!(net = net, fromBus = "B1", toBus = "B3", length_km = 100.0, r_pu = 0.0653, x_pu = 0.398, c_nf_per_km = 9.08, tanδ = 0.0, ratedS = 250.0, status = 1)
  add_test_line_pu_c!(net = net, fromBus = "B2", toBus = "B4", length_km = 100.0, r_pu = 0.0653, x_pu = 0.398, c_nf_per_km = 9.08, tanδ = 0.0, ratedS = 250.0, status = 1)
  add_test_line_pu_c!(net = net, fromBus = "B3", toBus = "B4", length_km = 100.0, r_pu = 0.0653, x_pu = 0.398, c_nf_per_km = 9.08, tanδ = 0.0, ratedS = 250.0, status = 1)
  add_test_line_pu_c!(net = net, fromBus = "B5", toBus = "B4", length_km = 300.0, r_pu = 0.0653, x_pu = 0.398, c_nf_per_km = 9.08, tanδ = 0.0, ratedS = 250.0, status = 1)
  add_test_line_pu_c!(net = net, fromBus = "B5", toBus = "B4", length_km = 300.0, r_pu = 0.0653, x_pu = 0.398, c_nf_per_km = 9.08, tanδ = 0.0, ratedS = 250.0, status = 1)

  addShunt!(net = net, busName = "B3", pShunt = 0.0, qShunt = 180.0)
  addShunt!(net = net, busName = "B2", pShunt = 0.0, qShunt = 180.0)
  addShunt!(net = net, busName = "B4", pShunt = 0.0, qShunt = 180.0)

  addProsumer!(net = net, busName = "B3", type = "ENERGYCONSUMER", p = 285.0, q = 200.0)
  addProsumer!(net = net, busName = "B4", type = "ENERGYCONSUMER", p = 103.0, q = 62.0)

  addProsumer!(net = net, busName = "B5", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.03, va_deg = 0.0, referencePri = "B5")
  addProsumer!(net = net, busName = "B2", type = "SYNCHRONOUSMACHINE", p = 600.0, vm_pu = 1.03, va_deg = 0.0)

  result, msg = validate!(net = net, log = true)
  if !result
    @warn msg
    return false
  end
  tol = 1e-6
  maxIte = 10
  verbose = 0       # 0: no output, 1: iteration norm, 2: + Y-Bus, 3: + Jacobian, 4: + Power Flow

  etime = @elapsed begin
    ite, erg = runpf!(net, maxIte, tol, verbose)
  end
  if erg != 0
    @warn "Power flow did not converge"
    printACPFlowResults(net, etime, ite, tol; converged = false)
    return false
  end

  brVec = getNetBranchNumberVec(net = net, fromBus = "B1", toBus = "B2")
  setNetBranchStatus!(net = net, branchNr = brVec[1], status = 0)
  brVec = getNetBranchNumberVec(net = net, fromBus = "B1", toBus = "B3")
  setNetBranchStatus!(net = net, branchNr = brVec[1], status = 0)

  if length(net.isoNodes) != 1
    @warn "Expected 1 isolated node, found: $(length(net.isoNodes))"
    return false
  end
  if net.isoNodes[1] != 1
    @warn "Expected isolated node 1, found: $(net.isoNodes[1])"
    return false
  end

  return true
end

function test_bus_shunt_model_modes()::Bool
  y_expected = ComplexF64(10.0, 20.0) / 100.0

  net_adm = Net(name = "bus_shunt_admittance", baseMVA = 100.0, bus_shunt_model = "admittance")
  addBus!(net = net_adm, busName = "B1", vn_kV = 110.0)
  addShunt!(net = net_adm, busName = "B1", pShunt = 10.0, qShunt = 20.0)
  Y_adm = createYBUS(net = net_adm, sparse = false, printYBUS = false)

  net_vdi = Net(name = "bus_shunt_vdi", baseMVA = 100.0, bus_shunt_model = "voltage_dependent_injection")
  addBus!(net = net_vdi, busName = "B1", vn_kV = 110.0)
  addShunt!(net = net_vdi, busName = "B1", pShunt = 10.0, qShunt = 20.0)
  Y_vdi = createYBUS(net = net_vdi, sparse = false, printYBUS = false)

  S1, dP1, dQ1 = buildControlledSVec(net_vdi, ComplexF64[1.0+0.0im])
  S2, dP2, dQ2 = buildControlledSVec(net_vdi, ComplexF64[2.0+0.0im])

  mpc = (
    name = "case2_shunt_modes",
    baseMVA = 100.0,
    bus = [
      1 3 0.0 0.0 10.0 20.0 1 1.0 0.0 110.0 1 1.10 0.90
      2 1 0.0 0.0 0.0 0.0 1 1.0 0.0 110.0 1 1.10 0.90
    ],
    gen = [1 0.0 0.0 100.0 -100.0 1.0 100.0 1 100.0 0.0 0 0 0 0 0 0 0 0 0 0 0],
    branch = [1 2 0.02 0.06 0.0 0 0 0 0.0 0.0 1 -360 360],
    gencost = nothing,
    bus_name = ["BUS1", "BUS2"],
  )
  mp_adm = Sparlectra.createNetFromMatPowerCase(mpc = mpc, bus_shunt_model = "admittance")
  mp_vdi = Sparlectra.createNetFromMatPowerCase(mpc = mpc, bus_shunt_model = "voltage_dependent_injection")

  invalid_ok = false
  try
    Sparlectra.normalize_bus_shunt_model("not_a_mode")
  catch err
    invalid_ok = err isa ErrorException && occursin("Invalid bus_shunt_model", sprint(showerror, err))
  end

  return isapprox(Y_adm[1, 1], y_expected; atol = 1e-12, rtol = 0.0) &&
         isapprox(Y_vdi[1, 1], 0.0 + 0.0im; atol = 1e-12, rtol = 0.0) &&
         isapprox(S1[1], -conj(y_expected); atol = 1e-12, rtol = 0.0) &&
         isapprox(S2[1], -4.0 * conj(y_expected); atol = 1e-12, rtol = 0.0) &&
         isapprox(dP1[1], -2.0 * real(conj(y_expected)); atol = 1e-12, rtol = 0.0) &&
         isapprox(dQ2[1], -4.0 * imag(conj(y_expected)); atol = 1e-12, rtol = 0.0) &&
         mp_adm.shuntVec[1].model == :Y &&
         mp_vdi.shuntVec[1].model == :VoltageDependentInjection &&
         isapprox(createYBUS(net = mp_adm, sparse = false, printYBUS = false)[1, 1] - createYBUS(net = mp_vdi, sparse = false, printYBUS = false)[1, 1], y_expected; atol = 1e-12, rtol = 0.0) &&
         invalid_ok
end

function test_matpower_file_import_honors_explicit_overrides()::Bool
  tmpdir = mktempdir()
  case_path = joinpath(tmpdir, "case_inline.m")
  write(
    case_path,
    """
function mpc = case_inline
mpc.version = '2';
mpc.baseMVA = 100;
mpc.bus = [
1 3 0 0 0 0 1 1.0 0 110 1 1.1 0.9;
2 1 90 30 5 10 1 1.0 0 110 1 1.1 0.9;
];
mpc.gen = [
1 100 0 300 -300 1.02 100 1 300 0;
];
mpc.branch = [
1 2 0.01 0.05 0.0 999 999 999 1.0 20 1 -360 360;
];
""",
  )
  prev_cfg = Sparlectra.active_sparlectra_config()
  global_cfg = Sparlectra.SparlectraConfig(Dict("power_flow" => Dict("start_mode" => Dict("flatstart" => true)), "model" => Dict("bus_shunt_model" => "voltage_dependent_injection"), "matpower_import" => Dict("shift_sign" => -1.0, "shift_unit" => "rad")))
  run_cfg = Sparlectra.SparlectraConfig(Dict("power_flow" => Dict("start_mode" => Dict("flatstart" => false)), "model" => Dict("bus_shunt_model" => "admittance"), "matpower_import" => Dict("shift_sign" => 1.0, "shift_unit" => "deg")))
  try
    Sparlectra.set_sparlectra_config!(global_cfg)

    net_direct = Sparlectra.createNetFromMatPowerFile(filename = case_path, flatstart = false, bus_shunt_model = :admittance, matpower_shift_unit = :deg)
    direct_ok = (net_direct.flatstart == false) && (net_direct.bus_shunt_model == :admittance)

    result = Sparlectra.run_sparlectra(casefile = basename(case_path), path = dirname(case_path), config = run_cfg)
    cfg_ok = (result isa Sparlectra.SparlectraRunResult) && (result.net.flatstart == false) && (result.net.bus_shunt_model == :admittance)
    return direct_ok && cfg_ok
  finally
    Sparlectra.set_sparlectra_config!(prev_cfg)
  end
end

function test_run_sparlectra_forwards_wrong_branch_config()::Bool
  tmpdir = mktempdir()
  case_path = joinpath(tmpdir, "case_wrong_branch_forwarding.m")
  write(
    case_path,
    """
function mpc = case_wrong_branch_forwarding
mpc.version = '2';
mpc.baseMVA = 100;
mpc.bus = [
1 3 0 0 0 0 1 1.0 0 110 1 1.1 0.9;
2 1 90 30 5 10 1 1.0 0 110 1 1.1 0.9;
];
mpc.gen = [
1 100 0 300 -300 1.02 100 1 300 0;
];
mpc.branch = [
1 2 0.01 0.05 0.0 999 999 999 1.0 20 1 -360 360;
];
""",
  )

  function run_with_detection(mode::Symbol)
    cfg = Sparlectra.SparlectraConfig(Dict("power_flow" => Dict("wrong_branch_detection" => String(mode), "wrong_branch_max_branch_angle_deg" => 0.001), "output" => Dict("logfile_results" => "full")))
    profile = Dict{Symbol,Any}(:enabled => true)
    result, output = mktemp() do output_path, io
      result = redirect_stdout(io) do
        Sparlectra.run_sparlectra(casefile = basename(case_path), path = dirname(case_path), config = cfg, performance_profile = profile)
      end
      flush(io)
      return result, read(output_path, String)
    end
    return result, Sparlectra.rectangular_pf_status(result.net), profile, output
  end

  result_warn, st_warn, profile_warn, output_warn = run_with_detection(:warn)
  row_warn = Sparlectra._sparlectra_status_row(result_warn)
  warn_ok =
    st_warn.wrong_branch_detection === :warn &&
    st_warn.branch_quality_metrics.max_branch_angle_deg == 0.001 &&
    st_warn.wrong_branch_status === :warn &&
    st_warn.wrong_branch_reason === :branch_angle_exceeded &&
    st_warn.wrong_branch_branch_angle_violation_count > 0 &&
    st_warn.final_converged === true &&
    result_warn.numerical_converged === true &&
    result_warn.solution_available === true &&
    result_warn.final_converged === true &&
    row_warn.numerical_converged === true &&
    row_warn.converged === true &&
    row_warn.erg == 0 &&
    row_warn.solution_available === true &&
    row_warn.q_limit_active_set_ok === st_warn.q_limit_active_set_ok &&
    row_warn.active_set_converged === st_warn.active_set_converged &&
    row_warn.pv_q_limit_violations == st_warn.pv_q_limit_violations &&
    row_warn.ref_q_limit_violations == st_warn.ref_q_limit_violations &&
    isequal(row_warn.final_pv_voltage_residual, st_warn.final_pv_voltage_residual) &&
    haskey(profile_warn[:timings], :postprocess_losses_and_flows) &&
    occursin("AC Power Flow Results", output_warn)

  result_fail, st_fail, profile_fail, output_fail = run_with_detection(:fail)
  row_fail = Sparlectra._sparlectra_status_row(result_fail)
  fail_ok =
    st_fail.wrong_branch_detection === :fail &&
    st_fail.branch_quality_metrics.max_branch_angle_deg == 0.001 &&
    st_fail.wrong_branch_status === :fail &&
    st_fail.wrong_branch_reason === :wrong_branch_detected &&
    st_fail.wrong_branch_branch_angle_violation_count > 0 &&
    st_fail.numerical_converged === true &&
    st_fail.final_converged === false &&
    st_fail.status === :wrong_branch_detected &&
    st_fail.reason === :wrong_branch_detected &&
    result_fail.numerical_converged === true &&
    result_fail.solution_available === false &&
    result_fail.final_converged === false &&
    row_fail.numerical_converged === true &&
    row_fail.final_converged === false &&
    row_fail.converged === false &&
    row_fail.erg == 1 &&
    row_fail.solution_available === false &&
    row_fail.converged === row_fail.final_converged &&
    row_fail.erg == (row_fail.final_converged ? 0 : 1) &&
    row_fail.numerical_converged === result_fail.numerical_converged &&
    row_fail.solution_available === result_fail.solution_available &&
    row_fail.outcome === result_fail.outcome &&
    row_fail.reason === result_fail.reason &&
    row_fail.wrong_branch_detection === st_fail.wrong_branch_detection &&
    row_fail.wrong_branch_status === st_fail.wrong_branch_status &&
    row_fail.wrong_branch_reason === st_fail.wrong_branch_reason &&
    row_fail.wrong_branch_angle_spread_deg == st_fail.wrong_branch_angle_spread_deg &&
    row_fail.wrong_branch_max_branch_angle_deg == st_fail.wrong_branch_max_branch_angle_deg &&
    row_fail.wrong_branch_branch_angle_violation_count == st_fail.wrong_branch_branch_angle_violation_count &&
    row_fail.wrong_branch_worst_branch_angle_deg == st_fail.wrong_branch_worst_branch_angle_deg &&
    row_fail.wrong_branch_rescue_attempted === st_fail.wrong_branch_rescue_attempted &&
    row_fail.wrong_branch_rescue_reason === st_fail.wrong_branch_rescue_reason &&
    !haskey(profile_fail[:timings], :postprocess_losses_and_flows) &&
    !haskey(profile_fail[:timings], :result_output) &&
    !occursin("AC Power Flow Results", output_fail)

  # the check runs on every numerically converged AC result, Q limits off
  # included (it sat inside the Q-limit block before 0.20.5, so a flat
  # start that converged onto a low-voltage branch with Q limits off was
  # reported as not_checked)
  cfg_noq = Sparlectra.SparlectraConfig(Dict("power_flow" => Dict("wrong_branch_detection" => "warn", "wrong_branch_max_branch_angle_deg" => 0.001, "qlimits" => Dict("enabled" => false)), "output" => Dict("logfile_results" => "off")))
  result_noq = redirect_stdout(devnull) do
    Sparlectra.run_sparlectra(casefile = basename(case_path), path = dirname(case_path), config = cfg_noq)
  end
  st_noq = Sparlectra.rectangular_pf_status(result_noq.net)
  noq_ok = st_noq.wrong_branch_status === :warn && st_noq.wrong_branch_reason === :branch_angle_exceeded && result_noq.final_converged === true

  # issue #462: the reference-branch bound travels the same chain
  # (configuration -> runpf! -> final status). The branch-angle rule stays
  # at its 90 degrees, so the 1-2 branch at the reference is caught by the
  # new rule alone; warn only, the run stays converged
  cfg_ref = Sparlectra.SparlectraConfig(Dict("power_flow" => Dict("wrong_branch_detection" => "warn", "wrong_branch_max_reference_branch_angle_deg" => 0.001, "qlimits" => Dict("enabled" => false)), "output" => Dict("logfile_results" => "off")))
  result_ref = redirect_stdout(devnull) do
    Sparlectra.run_sparlectra(casefile = basename(case_path), path = dirname(case_path), config = cfg_ref)
  end
  st_ref = Sparlectra.rectangular_pf_status(result_ref.net)
  ref_ok = st_ref.wrong_branch_status === :warn && st_ref.wrong_branch_reason === :reference_branch_angle_exceeded && st_ref.wrong_branch_reference_branch == "1-2" && st_ref.wrong_branch_reference_branch_angle_deg > 0.001 && result_ref.final_converged === true

  # the status reaches the report metadata: ok on a clean solve, not_checked
  # with the check off, warn with a band no solved bus meets
  function report_status(detection; kwargs...)
    net = createTest3BusNet()
    _, erg = runpf!(net, 40, 1e-8, 0; method = :rectangular, wrong_branch_detection = detection, kwargs...)
    erg == 0 || return nothing
    return buildACPFlowReport(net; converged = true, solver = :rectangular).metadata
  end
  md_clean = report_status(:warn)
  md_off = report_status(:off)
  md_suspect = report_status(:warn; wrong_branch_min_vm_pu = 1.20)
  report_ok =
    md_clean !== nothing && md_clean.wrong_branch_status == :ok &&
    md_off !== nothing && md_off.wrong_branch_status == :not_checked &&
    md_suspect !== nothing && md_suspect.wrong_branch_status == :warn && md_suspect.wrong_branch_reason == :low_voltage_magnitude

  # and the AC island diagnostics CSV of an island-wise run carries it per island
  island_net = Net(name = "wrong_branch_islands", baseMVA = 100.0)
  for busName in ("A1", "A2", "B1", "B2")
    addBus!(net = island_net, busName = busName, vn_kV = 110.0)
  end
  addPIModelACLine!(net = island_net, fromBus = "A1", toBus = "A2", r_pu = 0.01, x_pu = 0.10, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = island_net, fromBus = "B1", toBus = "B2", r_pu = 0.01, x_pu = 0.10, b_pu = 0.0, status = 1)
  addProsumer!(net = island_net, busName = "A1", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "A1")
  addProsumer!(net = island_net, busName = "A2", type = "ENERGYCONSUMER", p = 10.0, q = 3.0)
  addProsumer!(net = island_net, busName = "B1", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "B1")
  addProsumer!(net = island_net, busName = "B2", type = "ENERGYCONSUMER", p = 8.0, q = 2.0)
  refreshBusTypesFromProsumers!(island_net)
  island_header, island_rows = mktempdir() do tmpdir
    profile = Dict{Symbol,Any}(:output_dir => tmpdir)
    island_cfg = SparlectraConfig(powerflow = PowerFlowConfig(max_iter = 40, islands_enabled = true, wrong_branch_detection = :warn, wrong_branch_min_vm_pu = 1.20))
    island_result = run_sparlectra(net = island_net, config = island_cfg, performance_profile = profile)
    island_result.final_converged || return (String[], Vector{String}[])
    csv = run_csv(joinpath(tmpdir, "ac_island_solver_summary.csv"))
    return (csv.header, csv.rows)
  end
  status_col = findfirst(==("wrong_branch_status"), island_header)
  csv_ok = status_col !== nothing && any(row -> row[status_col] == "warn", island_rows)

  return warn_ok && fail_ok && noq_ok && ref_ok && report_ok && csv_ok
end

# The voltage band of the wrong-branch check judges every energised level,
# the angle heuristics the highest one: a 110 kV bus at 0.30 pu behind a
# clean 380 kV level is a wrong branch (case13659pegase's 150 kV collapse
# under three clean 750 kV buses read as ok before 0.20.5).
function test_wrong_branch_voltage_band_every_level()::Bool
  net = Net(name = "wb_levels", baseMVA = 100.0)
  addBus!(net = net, busName = "HV1", vn_kV = 380.0)
  addBus!(net = net, busName = "HV2", vn_kV = 380.0)
  addBus!(net = net, busName = "MV", vn_kV = 110.0)
  addPIModelACLine!(net = net, fromBus = "HV1", toBus = "HV2", r_pu = 0.01, x_pu = 0.1, b_pu = 0.0, status = 1)
  addPIModelTrafo!(net = net, fromBus = "HV2", toBus = "MV", r_pu = 0.005, x_pu = 0.05, b_pu = 0.0, status = 1)
  addProsumer!(net = net, busName = "HV1", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "HV1")
  addProsumer!(net = net, busName = "MV", type = "ENERGYCONSUMER", p = 10.0, q = 3.0)
  refreshBusTypesFromProsumers!(net)
  V = ComplexF64[1.0 + 0.0im, 0.99 - 0.01im, 0.30 - 0.05im]
  bus_types = [:Slack, :PQ, :PQ]
  Vset = [1.0, 1.0, 1.0]
  res = Sparlectra._check_wrong_branch_solution(V, bus_types, Vset, 1; net = net, min_vm_pu = 0.70, max_vm_pu = 1.30, max_angle_spread_deg = 180.0, max_branch_angle_deg = 90.0, min_low_vm_count = 1)
  # the 110 kV level (one bus, below the band) trips the level-share rule
  # first; the level and its counts are reported
  share_ok = res.status === :warn && res.reason === :low_voltage_magnitude && res.level_kV == 110.0 && res.level_low_vm_count == 1 && res.level_bus_count == 1 && res.low_vm_count == 1 && res.min_vm_pu < 0.31 && 3 in res.lowest_buses
  # with the single-bus count above the low buses present the level-share
  # rule carries the finding (the 110 kV level is fully below the band)
  res_share = Sparlectra._check_wrong_branch_solution(V, bus_types, Vset, 1; net = net, min_vm_pu = 0.70, max_vm_pu = 1.30, max_angle_spread_deg = 180.0, max_branch_angle_deg = 90.0, min_low_vm_count = 5, low_vm_share = 0.05)
  share_ok = share_ok && res_share.status === :warn && res_share.reason === :low_voltage_level_share
  # with the floor above 110 kV the low bus leaves the judged set: clean
  res_floor = Sparlectra._check_wrong_branch_solution(V, bus_types, Vset, 1; net = net, min_vm_pu = 0.70, max_vm_pu = 1.30, max_angle_spread_deg = 180.0, max_branch_angle_deg = 90.0, min_low_vm_count = 1, min_vn_kV = 220.0)
  floor_ok = res_floor.status === :ok && res_floor.low_vm_count == 0
  # below the floor a collapsed bus is still a finding through the collapse
  # threshold (case2848rte's branch sits on its 63 kV level at 0.02 pu)
  res_collapse = Sparlectra._check_wrong_branch_solution(V, bus_types, Vset, 1; net = net, min_vm_pu = 0.70, max_vm_pu = 1.30, max_angle_spread_deg = 180.0, max_branch_angle_deg = 90.0, min_low_vm_count = 1, min_vn_kV = 220.0, collapse_vm_pu = 0.5)
  floor_ok = floor_ok && res_collapse.status === :warn && res_collapse.reason === :voltage_collapse && res_collapse.low_vm_count == 1
  # a single bus angle beyond the bound is a finding of its own
  Vang = ComplexF64[1.0 + 0.0im, 0.99 - 0.01im, -0.5 + 0.86im]
  res_ang = Sparlectra._check_wrong_branch_solution(Vang, bus_types, Vset, 1; net = net, min_vm_pu = 0.70, max_vm_pu = 1.30, max_angle_spread_deg = 180.0, max_branch_angle_deg = 1e9, min_low_vm_count = 1, max_bus_angle_deg = 120.0)
  ang_ok = res_ang.status === :warn && res_ang.reason === :bus_angle_exceeded && res_ang.max_bus_angle_deg > 119.0
  # plain steps above the bound are reported without changing the status
  res_steps = Sparlectra._check_wrong_branch_solution(ComplexF64[1.0 + 0.0im, 0.99 - 0.01im, 0.98 - 0.02im], bus_types, Vset, 1; net = net, min_vm_pu = 0.70, max_vm_pu = 1.30, max_angle_spread_deg = 180.0, max_branch_angle_deg = 90.0, min_low_vm_count = 1, max_plain_steps = 20, plain_steps = 25)
  steps_ok = res_steps.status === :ok && res_steps.plain_steps == 25 && res_steps.plain_steps_exceeded === true
  return share_ok && floor_ok && ang_ok && steps_ok
end

# Issue #462: the angle across the reference bus's own branches. The shipped
# sp_casePST has the topology of case13659pegase's reference: the machine at
# bus 1 (15.5 kV) reaches the 230 kV network through its own transformer,
# a branch the highest-level branch-angle rule never sees. Its operating
# root stays quiet; the same state rotated as a whole by 164 degrees
# against the reference (case13659pegase's rotated root: magnitudes of the
# operating state, every other bus 164.5 degrees away) is flagged with the
# new reason, ahead of bus_angle_exceeded, and naming the transformer.
function test_wrong_branch_reference_branch_angle()::Bool
  case = joinpath(@__DIR__, "..", "data", "scf", "sp_casePST.scf.json")
  net = Sparlectra.importSCF(case)
  cfg = Sparlectra.SparlectraConfig(Dict("output" => Dict("logfile_results" => "off", "console_summary" => false)))
  r = redirect_stdout(devnull) do
    Sparlectra.run_sparlectra(net = net, config = cfg)
  end
  st = Sparlectra.rectangular_pf_status(net)
  # template default 90 degrees; the operating root sits at 1.4 degrees
  operating_ok = r.final_converged && st.wrong_branch_status === :ok && st.wrong_branch_reference_branch == "1-2" && st.wrong_branch_reference_branch_angle_deg < 10.0
  slack = findfirst(n -> Sparlectra.getNodeType(n) == Sparlectra.Slack, net.nodeVec)
  V = ComplexF64[n._vm_pu * cis(deg2rad(n._va_deg)) for n in net.nodeVec]
  Vrot = ComplexF64[i == slack ? V[i] : V[i] * cis(deg2rad(164.0)) for i in eachindex(V)]
  bus_types = [i == slack ? :Slack : :PQ for i in eachindex(V)]
  check(Vs, bound) = Sparlectra._check_wrong_branch_solution(Vs, bus_types, ones(length(Vs)), slack; net = net, min_vm_pu = 0.70, max_vm_pu = 1.30, max_angle_spread_deg = 180.0, max_branch_angle_deg = 90.0, min_low_vm_count = 1, max_bus_angle_deg = 120.0, collapse_vm_pu = 0.5, max_reference_branch_angle_deg = bound)
  rot = check(Vrot, 90.0)
  rotated_ok = rot.status === :warn && rot.reason === :reference_branch_angle_exceeded && rot.reference_branch == "1-2" && 160.0 < rot.reference_branch_angle_deg < 170.0
  # the rule off (0): the same state falls back to the bus-angle rule, the
  # angle is still measured
  off = check(Vrot, 0.0)
  off_ok = off.reason === :bus_angle_exceeded && off.reference_branch_angle_deg == rot.reference_branch_angle_deg
  return operating_ok && rotated_ok && off_ok
end

function test_link_kcl_simple()
  net = Net(name = "link_kcl", baseMVA = 100.0)
  addBus!(net = net, busName = "B1", vn_kV = 110.0)
  addBus!(net = net, busName = "B2", vn_kV = 110.0)

  addBusGenPower!(net = net, busName = "B1", p = 10.0, q = 2.0)
  addBusLoadPower!(net = net, busName = "B2", p = 10.0, q = 2.0)

  addLink!(net = net, fromBus = "B1", toBus = "B2", status = 1)
  calcLinkFlowsKCL!(net)

  l = net.linkVec[1]
  iexp = hypot(10.0, 2.0) / (Sparlectra.Wurzel3 * 110.0)
  type_ok = (l isa AbstractComponent) && l.cTyp == Sparlectra.LinkC && toComponentTyp("BUSLINK") == Sparlectra.LinkC && toComponentTyp("LINK") == Sparlectra.LinkC
  return type_ok && isapprox(l.pFlow_MW, 10.0; atol = 1e-8) && isapprox(l.qFlow_MVar, 2.0; atol = 1e-8) && isapprox(l.iFrom_kA, iexp; atol = 1e-10) && isapprox(l.iTo_kA, iexp; atol = 1e-10)
end

# addLink! refuses a link at a slack bus and between buses of different
# types (PV and PQ) with an AssertionError naming the reason.
function test_link_rejects_slack_connection()
  function refusal(type1, kw1, needle)
    net = Net(name = "link_reject", baseMVA = 100.0)
    addBus!(net = net, busName = "B1", vn_kV = 110.0)
    addBus!(net = net, busName = "B2", vn_kV = 110.0)
    addProsumer!(; net = net, busName = "B1", type = type1, kw1...)
    addProsumer!(net = net, busName = "B2", type = "ENERGYCONSUMER", p = 5.0, q = 1.0)
    try
      addLink!(net = net, fromBus = "B1", toBus = "B2", status = 1)
      return false
    catch err
      err isa AssertionError || rethrow()
      return occursin(needle, lowercase(string(err)))
    end
  end
  slack_ok = refusal("EXTERNALNETWORKINJECTION", (vm_pu = 1.0, va_deg = 0.0, referencePri = "B1"), "slack")
  mixed_ok = refusal("SYNCHRONOUSMACHINE", (p = 10.0, vm_pu = 1.02), "same bus type")
  return slack_ok && mixed_ok
end

function test_link_bus_merge_pf()
  tol = 1e-8
  maxIte = 30

  net_link = Net(name = "link_merge_pf", baseMVA = 100.0)
  addBus!(net = net_link, busName = "B0", vn_kV = 110.0)
  addBus!(net = net_link, busName = "B1", vn_kV = 110.0)
  addBus!(net = net_link, busName = "B2", vn_kV = 110.0)
  addBus!(net = net_link, busName = "B3", vn_kV = 110.0)

  add_test_line_pu_c!(net = net_link, fromBus = "B0", toBus = "B1", length_km = 5.0, r_pu = 0.05, x_pu = 0.5, c_nf_per_km = 10.0, tanδ = 0.0)
  add_test_line_pu_c!(net = net_link, fromBus = "B2", toBus = "B3", length_km = 20.0, r_pu = 0.05, x_pu = 0.5, c_nf_per_km = 10.0, tanδ = 0.0)
  addLink!(net = net_link, fromBus = "B1", toBus = "B2", status = 1)

  addProsumer!(net = net_link, busName = "B0", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "B0")
  addProsumer!(net = net_link, busName = "B1", type = "GENERATOR", p = 10.0, q = 1.0)
  addProsumer!(net = net_link, busName = "B2", type = "ENERGYCONSUMER", p = 15.0, q = 5.0)
  addProsumer!(net = net_link, busName = "B3", type = "ENERGYCONSUMER", p = 25.0, q = 10.0)

  ite_link, erg_link = runpf!(net_link, maxIte, tol, 0; method = :rectangular)
  if erg_link != 0
    @warn "Link network did not converge" ite = ite_link
    return false
  end

  b1_idx = geNetBusIdx(net = net_link, busName = "B1")
  b2_idx = geNetBusIdx(net = net_link, busName = "B2")
  b3_idx = geNetBusIdx(net = net_link, busName = "B3")

  if !isapprox(net_link.nodeVec[b1_idx]._vm_pu, net_link.nodeVec[b2_idx]._vm_pu; atol = 1e-10)
    @warn "Merged link buses do not share voltage magnitude" vm1 = net_link.nodeVec[b1_idx]._vm_pu vm2 = net_link.nodeVec[b2_idx]._vm_pu
    return false
  end
  if !isapprox(net_link.nodeVec[b1_idx]._va_deg, net_link.nodeVec[b2_idx]._va_deg; atol = 1e-8)
    @warn "Merged link buses do not share voltage angle" va1 = net_link.nodeVec[b1_idx]._va_deg va2 = net_link.nodeVec[b2_idx]._va_deg
    return false
  end

  net_eq = Net(name = "merge_reference", baseMVA = 100.0)
  addBus!(net = net_eq, busName = "B0", vn_kV = 110.0)
  addBus!(net = net_eq, busName = "B12", vn_kV = 110.0)
  addBus!(net = net_eq, busName = "B3", vn_kV = 110.0)

  add_test_line_pu_c!(net = net_eq, fromBus = "B0", toBus = "B12", length_km = 5.0, r_pu = 0.05, x_pu = 0.5, c_nf_per_km = 10.0, tanδ = 0.0)
  add_test_line_pu_c!(net = net_eq, fromBus = "B12", toBus = "B3", length_km = 20.0, r_pu = 0.05, x_pu = 0.5, c_nf_per_km = 10.0, tanδ = 0.0)
  addProsumer!(net = net_eq, busName = "B0", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "B0")
  addProsumer!(net = net_eq, busName = "B12", type = "GENERATOR", p = 10.0, q = 1.0)
  addProsumer!(net = net_eq, busName = "B12", type = "ENERGYCONSUMER", p = 15.0, q = 5.0)
  addProsumer!(net = net_eq, busName = "B3", type = "ENERGYCONSUMER", p = 25.0, q = 10.0)

  ite_eq, erg_eq = runpf!(net_eq, maxIte, tol, 0; method = :rectangular)
  if erg_eq != 0
    @warn "Reference network did not converge" ite = ite_eq
    return false
  end

  b12_idx = geNetBusIdx(net = net_eq, busName = "B12")
  b3_eq_idx = geNetBusIdx(net = net_eq, busName = "B3")

  vm_link_b3 = net_link.nodeVec[b3_idx]._vm_pu
  va_link_b3 = net_link.nodeVec[b3_idx]._va_deg
  vm_eq_b3 = net_eq.nodeVec[b3_eq_idx]._vm_pu
  va_eq_b3 = net_eq.nodeVec[b3_eq_idx]._va_deg

  return isapprox(net_link.nodeVec[b1_idx]._vm_pu, net_eq.nodeVec[b12_idx]._vm_pu; atol = 1e-8) && isapprox(net_link.nodeVec[b1_idx]._va_deg, net_eq.nodeVec[b12_idx]._va_deg; atol = 1e-7) && isapprox(vm_link_b3, vm_eq_b3; atol = 1e-8) && isapprox(va_link_b3, va_eq_b3; atol = 1e-7)
end

# Builds a net whose slack sits on B2, which the link merge contracts onto B1, so
# the slack is deliberately *not* the cluster representative. The link is added
# before the slack prosumer, otherwise addLink!'s slack assertion rejects it.
function _build_non_representative_slack_net(name::String)
  net = Net(name = name, baseMVA = 100.0)
  addBus!(net = net, busName = "B0", vn_kV = 110.0)
  addBus!(net = net, busName = "B1", vn_kV = 110.0)
  addBus!(net = net, busName = "B2", vn_kV = 110.0)
  addBus!(net = net, busName = "B3", vn_kV = 110.0)

  add_test_line_pu_c!(net = net, fromBus = "B0", toBus = "B1", length_km = 5.0, r_pu = 0.05, x_pu = 0.5, c_nf_per_km = 10.0, tanδ = 0.0)
  add_test_line_pu_c!(net = net, fromBus = "B2", toBus = "B3", length_km = 20.0, r_pu = 0.05, x_pu = 0.5, c_nf_per_km = 10.0, tanδ = 0.0)
  addLink!(net = net, fromBus = "B1", toBus = "B2", status = 1)

  addProsumer!(net = net, busName = "B2", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "B2")
  addProsumer!(net = net, busName = "B0", type = "ENERGYCONSUMER", p = 15.0, q = 5.0)
  addProsumer!(net = net, busName = "B3", type = "ENERGYCONSUMER", p = 25.0, q = 10.0)
  return net
end

# Reference for the net above with the link cluster already collapsed by hand.
function _build_non_representative_slack_reference()
  net = Net(name = "non_rep_slack_reference", baseMVA = 100.0)
  addBus!(net = net, busName = "B0", vn_kV = 110.0)
  addBus!(net = net, busName = "B12", vn_kV = 110.0)
  addBus!(net = net, busName = "B3", vn_kV = 110.0)

  add_test_line_pu_c!(net = net, fromBus = "B0", toBus = "B12", length_km = 5.0, r_pu = 0.05, x_pu = 0.5, c_nf_per_km = 10.0, tanδ = 0.0)
  add_test_line_pu_c!(net = net, fromBus = "B12", toBus = "B3", length_km = 20.0, r_pu = 0.05, x_pu = 0.5, c_nf_per_km = 10.0, tanδ = 0.0)

  addProsumer!(net = net, busName = "B12", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "B12")
  addProsumer!(net = net, busName = "B0", type = "ENERGYCONSUMER", p = 15.0, q = 5.0)
  addProsumer!(net = net, busName = "B3", type = "ENERGYCONSUMER", p = 25.0, q = 10.0)
  return net
end

# Regression for the slackVec remap in _merged_pf_net. Contracting a closed link
# picks the lower bus index as the cluster representative; a slack on the higher
# index is therefore neutralized and marked Isolated, and the Slack type moves to
# the representative. slackVec has to follow, otherwise it points at a bus the
# solver no longer sees.
function test_link_merge_remaps_non_representative_slack()::Bool
  net = _build_non_representative_slack_net("non_rep_slack_merge")
  # Bus typing is derived from prosumers, so slackVec is only populated here.
  Sparlectra.refreshBusTypesFromProsumers!(net)
  if net.slackVec != [3]
    @warn "Precondition failed: slack is not on bus 3" slackVec = net.slackVec
    return false
  end

  wnet, reps, has_merges = Sparlectra._merged_pf_net(net)
  # reps[3] == 2 is the whole point of the fixture: the slack is merged away.
  if !has_merges || reps[3] != 2
    @warn "Precondition failed: bus 3 is not merged onto bus 2" has_merges reps
    return false
  end

  # _merged_pf_net must return a self-consistent net on its own -- deliberately
  # checked before refreshBusTypesFromProsumers! runs, since that call would
  # rebuild slackVec from the prosumers and hide a stale entry.
  if wnet.slackVec != [2]
    @warn "Merged net keeps a stale slack index" slackVec = wnet.slackVec
    return false
  end
  if Sparlectra.getNodeType(wnet.nodeVec[3]) != Sparlectra.Isolated || Sparlectra.getNodeType(wnet.nodeVec[2]) != Sparlectra.Slack
    @warn "Merged net has inconsistent bus types" types = [n._nodeType for n in wnet.nodeVec]
    return false
  end

  tol = 1e-8
  ite, erg = runpf!(net, 30, tol, 0; method = :rectangular)
  if erg != 0
    @warn "Non-representative slack network did not converge" ite
    return false
  end

  b1 = geNetBusIdx(net = net, busName = "B1")
  b2 = geNetBusIdx(net = net, busName = "B2")
  if !isapprox(net.nodeVec[b1]._vm_pu, net.nodeVec[b2]._vm_pu; atol = 1e-10) || !isapprox(net.nodeVec[b1]._va_deg, net.nodeVec[b2]._va_deg; atol = 1e-8)
    @warn "Merged link buses do not share the solved voltage"
    return false
  end
  # The slack bus itself must still hold the reference voltage after back-sync.
  if !isapprox(net.nodeVec[b2]._vm_pu, 1.0; atol = 1e-10) || !isapprox(net.nodeVec[b2]._va_deg, 0.0; atol = 1e-8)
    @warn "Slack bus lost its reference voltage" vm = net.nodeVec[b2]._vm_pu va = net.nodeVec[b2]._va_deg
    return false
  end

  # Beyond the bookkeeping: the solved state must equal the network in which the
  # link cluster was collapsed by hand. This catches a remap that keeps slackVec
  # tidy but shifts the reference onto the wrong bus.
  ref = _build_non_representative_slack_reference()
  ite_ref, erg_ref = runpf!(ref, 30, tol, 0; method = :rectangular)
  if erg_ref != 0
    @warn "Reference network did not converge" ite_ref
    return false
  end

  ok = true
  # B1/B2 of the link net correspond to the single bus B12 of the reference.
  for (name, ref_name) in (("B0", "B0"), ("B1", "B12"), ("B3", "B3"))
    i = geNetBusIdx(net = net, busName = name)
    j = geNetBusIdx(net = ref, busName = ref_name)
    if !isapprox(net.nodeVec[i]._vm_pu, ref.nodeVec[j]._vm_pu; atol = 1e-8) || !isapprox(net.nodeVec[i]._va_deg, ref.nodeVec[j]._va_deg; atol = 1e-7)
      @warn "Solved state deviates from the hand-merged reference" bus = name vm = net.nodeVec[i]._vm_pu vm_ref = ref.nodeVec[j]._vm_pu va = net.nodeVec[i]._va_deg va_ref = ref.nodeVec[j]._va_deg
      ok = false
    end
  end
  return ok
end

# A bus B1 joined by a closed link to B2, which carries a load and a
# generator of `p2` MW, and with `linked = false` the same net with B1 and B2
# merged by hand into B1. Without `pv` B1 holds the slack and B2's unit is a
# PQ generator; with `pv` both buses carry a voltage-controlled machine with
# Q limits and the slack sits on B0.
function _link_power_net(linked::Bool; p2::Float64, pv::Bool = false, vm_pu::Float64 = 1.02)
  net = Net(name = linked ? "link_powers" : "link_powers_ref", baseMVA = 100.0)
  a = "B1"
  b = linked ? "B2" : "B1"
  for name in (linked ? ("B0", "B1", "B2", "B3") : ("B0", "B1", "B3"))
    addBus!(net = net, busName = name, vn_kV = 110.0)
  end
  # charging of 10 nF/km in per unit on 110 kV and 100 MVA (Z_base 121 ohm)
  b_pu(length_km) = 2π * 50.0 * 10.0e-9 * length_km * 110.0^2 / 100.0
  addPIModelACLine!(net = net, fromBus = "B0", toBus = a, r_pu = 0.05, x_pu = 0.5, b_pu = b_pu(5.0), status = 1)
  addPIModelACLine!(net = net, fromBus = b, toBus = "B3", r_pu = 0.05, x_pu = 0.5, b_pu = b_pu(20.0), status = 1)
  # the link goes in before the units: addLink! refuses a slack bus and
  # buses of different types
  linked && addLink!(net = net, fromBus = "B1", toBus = "B2", status = 1)
  if pv
    addProsumer!(net = net, busName = "B0", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "B0")
    addProsumer!(net = net, busName = a, type = "SYNCHRONOUSMACHINE", p = 20.0, vm_pu = vm_pu, qMin = -10.0, qMax = 30.0)
    addProsumer!(net = net, busName = b, type = "SYNCHRONOUSMACHINE", p = p2, vm_pu = vm_pu, qMin = -20.0, qMax = 20.0)
  else
    addProsumer!(net = net, busName = a, type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = a)
    addProsumer!(net = net, busName = b, type = "GENERATOR", p = p2, q = 1.0)
  end
  addProsumer!(net = net, busName = a, type = "ENERGYCONSUMER", p = 5.0, q = 1.0)
  addProsumer!(net = net, busName = b, type = "ENERGYCONSUMER", p = 15.0, q = 5.0)
  addProsumer!(net = net, busName = "B3", type = "ENERGYCONSUMER", p = 25.0, q = 10.0)
  addProsumer!(net = net, busName = "B0", type = "ENERGYCONSUMER", p = 8.0, q = 2.0)
  return net
end

# runpf! solves a net with a closed link on a merged copy. The solved bus
# powers have to come back to the members of the link cluster: the slack
# balance on the slack bus, every other member keeps its own generation and
# load, and a voltage-controlled pair shares its Q inside its limits. Only
# the voltages came back before, the slack bus kept its input value and the
# total bus power stayed empty.
function test_link_merge_returns_bus_powers()::Bool
  ok = true
  idx(net, name) = geNetBusIdx(net = net, busName = name)
  # B2 balanced in P (generation = load) and unbalanced
  for p2 in (15.0, 10.0)
    net = _link_power_net(true; p2 = p2)
    ref = _link_power_net(false; p2 = p2)
    ok &= last(runpf!(net, 30, 1e-10, 0; method = :rectangular)) == 0
    ok &= last(runpf!(ref, 30, 1e-10, 0; method = :rectangular)) == 0
    n1, n2 = net.nodeVec[idx(net, "B1")], net.nodeVec[idx(net, "B2")]
    r1 = ref.nodeVec[idx(ref, "B1")]
    # the slack bus kept its input value (here `nothing`) before the fix
    (n1._pƩGen isa Float64 && n1._qƩGen isa Float64) || return false
    # the slack bus value is the generation of its bus, so B2's own
    # scheduled generation (p2 MW, 1 MVAr) is the only difference to the
    # merged bus
    ok &= isapprox(n1._pƩGen, r1._pƩGen - p2; atol = 1e-8)
    ok &= isapprox(n1._qƩGen, r1._qƩGen - 1.0; atol = 1e-8)
    # each member keeps its own load and generation
    ok &= (n1._pƩLoad, n1._qƩLoad) == (5.0, 1.0)
    ok &= (n2._pƩGen, n2._qƩGen, n2._pƩLoad, n2._qƩLoad) == (p2, 1.0, 15.0, 5.0)
    ok &= all(isapprox.(getTotalBusPower(net = net), getTotalBusPower(net = ref); atol = 1e-8))
    # the units take the bus results: the slack unit that of its bus
    distributeBusResults!(net)
    ok &= isapprox(net.prosumpsVec[1].pRes, n1._pƩGen; atol = 1e-12) && net.prosumpsVec[2].pRes == p2
    # the merged-away member B2 keeps its own row in the report
    ok &= any(row -> row.bus_name == "B2", buildACPFlowReport(net; converged = true, solver = :rectangular).nodes)
  end
  # two linked PV machines: their Q sum is that of the merged bus, both at
  # the same relative position in their own Q bands
  net = _link_power_net(true; p2 = 10.0, pv = true)
  ref = _link_power_net(false; p2 = 10.0, pv = true)
  ok &= last(runpf!(net, 30, 1e-10, 0; method = :rectangular)) == 0
  ok &= last(runpf!(ref, 30, 1e-10, 0; method = :rectangular)) == 0
  n1, n2 = net.nodeVec[idx(net, "B1")], net.nodeVec[idx(net, "B2")]
  r1 = ref.nodeVec[idx(ref, "B1")]
  ok &= isapprox(n1._qƩGen + n2._qƩGen, r1._qƩGen; atol = 1e-8)
  ok &= (n1._pƩGen, n2._pƩGen) == (20.0, 10.0)
  # the bus values are the generator Q
  g1, g2 = n1._qƩGen, n2._qƩGen
  ok &= isapprox((g1 + 10.0) / 40.0, (g2 + 20.0) / 40.0; atol = 1e-10)
  # a setpoint the pair cannot hold: the cluster clamps at the sum of both
  # bands (50 MVAr) like the merged bus, at the same voltage, each machine at
  # its own Qmax. The merged copy kept the band of the representative alone
  # and clamped at its 30 MVAr, at another voltage.
  net = _link_power_net(true; p2 = 10.0, pv = true, vm_pu = 1.15)
  ref = _link_power_net(false; p2 = 10.0, pv = true, vm_pu = 1.15)
  ok &= last(runpf!(net, 30, 1e-10, 0; method = :rectangular)) == 0
  ok &= last(runpf!(ref, 30, 1e-10, 0; method = :rectangular)) == 0
  n1, n2 = net.nodeVec[idx(net, "B1")], net.nodeVec[idx(net, "B2")]
  r1 = ref.nodeVec[idx(ref, "B1")]
  ok &= get(ref.qLimitEvents, idx(ref, "B1"), nothing) === :max
  ok &= isapprox(r1._qƩGen, 50.0; atol = 1e-8) && isapprox(n1._vm_pu, r1._vm_pu; atol = 1e-8)
  ok &= isapprox(n1._qƩGen, 30.0; atol = 1e-8) && isapprox(n2._qƩGen, 20.0; atol = 1e-8)
  return ok
end

# Island-wise solving combined with link merges. Both features were tested
# separately before, never together, and the combination crashed in
# _sync_merged_results_to_original! *after* every island had converged. It is the
# realistic shape for CGMES deliveries: retained switches become links, and a
# partial delivery decomposes into several islands.
#
# The second island is a separate two-bus network with its own reference, so the
# island decomposition has to produce exactly two islands.
function test_island_solve_with_link_merge()::Bool
  net = _build_non_representative_slack_net("island_link_merge")
  addBus!(net = net, busName = "C0", vn_kV = 110.0)
  addBus!(net = net, busName = "C1", vn_kV = 110.0)
  add_test_line_pu_c!(net = net, fromBus = "C0", toBus = "C1", length_km = 10.0, r_pu = 0.05, x_pu = 0.5, c_nf_per_km = 10.0, tanδ = 0.0)
  addProsumer!(net = net, busName = "C0", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "C0")
  addProsumer!(net = net, busName = "C1", type = "ENERGYCONSUMER", p = 8.0, q = 2.0)

  Sparlectra.refreshBusTypesFromProsumers!(net)
  if length(Sparlectra.electricalIslandComponents(net)) != 2
    @warn "Precondition failed: expected two electrical islands" islands = Sparlectra.electricalIslandComponents(net)
    return false
  end

  # slackVec remapping itself is covered by the test above; here the merge only
  # has to be active so that the island path exercises the back-sync.
  if !Sparlectra._merged_pf_net(net)[3]
    @warn "Precondition failed: no link merge in the island network"
    return false
  end

  tol = 1e-8
  ite, erg = runpf!(net, 30, tol, 0; method = :rectangular, islands_enabled = true)
  if erg != 0
    @warn "Island-wise solve with link merges did not converge" ite
    return false
  end

  # B2 only receives a value through the representative map; a back-sync that
  # silently did not run would leave it at its flat-start value.
  b1 = geNetBusIdx(net = net, busName = "B1")
  b2 = geNetBusIdx(net = net, busName = "B2")
  if !isapprox(net.nodeVec[b1]._vm_pu, net.nodeVec[b2]._vm_pu; atol = 1e-10) || !isapprox(net.nodeVec[b1]._va_deg, net.nodeVec[b2]._va_deg; atol = 1e-8)
    @warn "Merged link buses were not synced back after island-wise solving"
    return false
  end

  # Islands are solved independently, so the second island must not shift the
  # first island's solution by a single digit. This also pins that the island
  # split did not accidentally move buses between islands.
  single = _build_non_representative_slack_net("island_link_merge_single")
  _, erg_single = runpf!(single, 30, tol, 0; method = :rectangular)
  if erg_single != 0
    @warn "Single-island comparison network did not converge"
    return false
  end

  ok = true
  for name in ("B0", "B1", "B2", "B3")
    i = geNetBusIdx(net = net, busName = name)
    j = geNetBusIdx(net = single, busName = name)
    if !isapprox(net.nodeVec[i]._vm_pu, single.nodeVec[j]._vm_pu; atol = 1e-8) || !isapprox(net.nodeVec[i]._va_deg, single.nodeVec[j]._va_deg; atol = 1e-7)
      @warn "Island-wise result deviates from the single-island solve" bus = name
      ok = false
    end
  end
  return ok
end

function test_link_closed_keeps_shunt_reporting_on_original_bus()
  net = Net(name = "link_shunt_reporting", baseMVA = 100.0)
  addBus!(net = net, busName = "B0", vn_kV = 110.0)
  addBus!(net = net, busName = "B1", vn_kV = 110.0)
  addBus!(net = net, busName = "B2", vn_kV = 110.0)
  addBus!(net = net, busName = "B3", vn_kV = 110.0)

  add_test_line_pu_c!(net = net, fromBus = "B0", toBus = "B1", length_km = 5.0, r_pu = 0.05, x_pu = 0.5, c_nf_per_km = 10.0, tanδ = 0.0)
  add_test_line_pu_c!(net = net, fromBus = "B2", toBus = "B3", length_km = 20.0, r_pu = 0.05, x_pu = 0.5, c_nf_per_km = 10.0, tanδ = 0.0)
  addLink!(net = net, fromBus = "B1", toBus = "B2", status = 1)
  addShunt!(net = net, busName = "B2", pShunt = 0.0, qShunt = 5.0)

  addProsumer!(net = net, busName = "B0", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "B0")
  addProsumer!(net = net, busName = "B1", type = "GENERATOR", p = 10.0, q = 1.0)
  addProsumer!(net = net, busName = "B2", type = "ENERGYCONSUMER", p = 15.0, q = 5.0)
  addProsumer!(net = net, busName = "B3", type = "ENERGYCONSUMER", p = 25.0, q = 10.0)

  _, erg = runpf!(net, 25, 1e-8, 0; method = :rectangular)
  erg == 0 || return false

  b2_idx = geNetBusIdx(net = net, busName = "B2")
  sh = net.shuntVec[1]
  vm = net.nodeVec[b2_idx]._vm_pu
  expected_q = -vm^2 * 5.0

  return isapprox(sh.q_shunt, expected_q; atol = 1e-8) && isapprox(net.nodeVec[b2_idx]._qShunt, expected_q; atol = 1e-8)
end

# A ring of three links: A f = b with the singular incidence matrix A of the
# ring, the flows are the minimum-norm solution pinv(A) b. Second row: a
# reactor at B2 (qShunt -6, MATPOWER sign: a positive value injects) draws
# 6 MVar at vm 1.0 pu, so b_B2 = -4 - 6 = -10 and the ring balances
# (12 - 10 - 2 = 0). Until 0.60.0 the allocation added the draw (b_B2 = +2
# with this reactor, and a 6 MVar capacitor gave -10).
function test_link_kcl_ring_allocation()
  A = [
    1.0 0.0 -1.0
    -1.0 1.0 0.0
    0.0 -1.0 1.0
  ]
  ok = true
  for (q1, shunt_q, bQ) in ((6.0, 0.0, [6.0, -4.0, -2.0]), (12.0, -6.0, [12.0, -10.0, -2.0]))
    net = Net(name = "link_kcl_ring", baseMVA = 100.0)
    addBus!(net = net, busName = "B1", vn_kV = 110.0)
    addBus!(net = net, busName = "B2", vn_kV = 110.0)
    addBus!(net = net, busName = "B3", vn_kV = 110.0)
    # P balance + 15 - 10 - 5 = 0
    addBusGenPower!(net = net, busName = "B1", p = 15.0, q = q1)
    addBusLoadPower!(net = net, busName = "B2", p = 10.0, q = 4.0)
    addBusLoadPower!(net = net, busName = "B3", p = 5.0, q = 2.0)
    if shunt_q != 0.0
      addShunt!(net = net, busName = "B2", pShunt = 0.0, qShunt = shunt_q)
      for node in net.nodeVec
        node._vm_pu = 1.0
      end
      updateShuntPowers!(net = net)
    end
    # B1 -> B2 -> B3 -> B1 forms a ring
    addLink!(net = net, fromBus = "B1", toBus = "B2", status = 1)
    addLink!(net = net, fromBus = "B2", toBus = "B3", status = 1)
    addLink!(net = net, fromBus = "B3", toBus = "B1", status = 1)
    # no power flow needed for the KCL allocation
    calcLinkFlowsKCL!(net)

    bP = [15.0, -10.0, -5.0]
    flowsP = [l.pFlow_MW for l in net.linkVec]
    flowsQ = [l.qFlow_MVar for l in net.linkVec]
    ok &= all(isapprox.(flowsP, pinv(A) * bP; atol = 1e-10, rtol = 0.0)) && all(isapprox.(flowsQ, pinv(A) * bQ; atol = 1e-10, rtol = 0.0))
    ok &= all(isapprox.(A * flowsP, bP; atol = 1e-10, rtol = 0.0)) && all(isapprox.(A * flowsQ, bQ; atol = 1e-10, rtol = 0.0))
  end
  return ok
end

function test_link_ring_pf_stability()
  # Build a mixed grid: radial AC lines from slack plus a meshed 3-link ring
  # among PQ buses where link flows are recovered via KCL.
  net = Net(name = "link_ring_pf_stability", baseMVA = 100.0)
  addBus!(net = net, busName = "B0", vn_kV = 110.0)
  addBus!(net = net, busName = "B1", vn_kV = 110.0)
  addBus!(net = net, busName = "B2", vn_kV = 110.0)
  addBus!(net = net, busName = "B3", vn_kV = 110.0)

  # Two physical AC branches feed the ring from the slack side.
  add_test_line_pu_c!(net = net, fromBus = "B0", toBus = "B1", length_km = 10.0, r_pu = 0.05, x_pu = 0.4, c_nf_per_km = 10.0, tanδ = 0.0)
  add_test_line_pu_c!(net = net, fromBus = "B0", toBus = "B2", length_km = 12.0, r_pu = 0.05, x_pu = 0.4, c_nf_per_km = 10.0, tanδ = 0.0)

  # Three links form a closed loop (meshed/singular incidence structure).
  addLink!(net = net, fromBus = "B1", toBus = "B2", status = 1)
  addLink!(net = net, fromBus = "B2", toBus = "B3", status = 1)
  addLink!(net = net, fromBus = "B3", toBus = "B1", status = 1)

  addProsumer!(net = net, busName = "B0", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "B0")
  addProsumer!(net = net, busName = "B1", type = "GENERATOR", p = 6.0, q = 1.0)
  addProsumer!(net = net, busName = "B2", type = "ENERGYCONSUMER", p = 18.0, q = 5.0)
  addProsumer!(net = net, busName = "B3", type = "ENERGYCONSUMER", p = 8.0, q = 3.0)

  # Solve PF first, then derive link flows from nodal balances.
  _, erg = runpf!(net, 40, 1e-9, 0; method = :rectangular)
  erg == 0 || return false

  calcLinkFlowsKCL!(net)

  b1 = geNetBusIdx(net = net, busName = "B1")
  b2 = geNetBusIdx(net = net, busName = "B2")
  b3 = geNetBusIdx(net = net, busName = "B3")

  # In this symmetric setup, ring buses should settle to (nearly) equal V magnitude/angle.
  vm_equal = isapprox(net.nodeVec[b1]._vm_pu, net.nodeVec[b2]._vm_pu; atol = 1e-10) && isapprox(net.nodeVec[b1]._vm_pu, net.nodeVec[b3]._vm_pu; atol = 1e-10)
  va_equal = isapprox(net.nodeVec[b1]._va_deg, net.nodeVec[b2]._va_deg; atol = 1e-8) && isapprox(net.nodeVec[b1]._va_deg, net.nodeVec[b3]._va_deg; atol = 1e-8)

  # Build link incidence matrix A for (B1,B2,B3) in link order.
  A = zeros(Float64, 3, 3)
  for (j, l) in enumerate(net.linkVec)
    A[l.fromBus-1, j] += 1.0
    A[l.toBus-1, j] -= 1.0
  end

  # Rebuild the KCL right-hand side used by link allocation:
  #   b = (gen - load - shunt draw) - branch_out
  P_inj = zeros(Float64, 3)
  Q_inj = zeros(Float64, 3)
  P_shunt = zeros(Float64, 3)
  Q_shunt = zeros(Float64, 3)
  P_branch_out = zeros(Float64, 3)
  Q_branch_out = zeros(Float64, 3)

  for (k, busidx) in enumerate((b1, b2, b3))
    node = net.nodeVec[busidx]
    P_inj[k] = (isnothing(node._pƩGen) ? 0.0 : node._pƩGen) - (isnothing(node._pƩLoad) ? 0.0 : node._pƩLoad)
    Q_inj[k] = (isnothing(node._qƩGen) ? 0.0 : node._qƩGen) - (isnothing(node._qƩLoad) ? 0.0 : node._qƩLoad)
    P_shunt[k] = isnothing(node._pShunt) ? 0.0 : node._pShunt
    Q_shunt[k] = isnothing(node._qShunt) ? 0.0 : node._qShunt
  end

  for br in net.branchVec
    br.status == 1 || continue
    if br.fromBus in (b1, b2, b3)
      k = br.fromBus == b1 ? 1 : (br.fromBus == b2 ? 2 : 3)
      if !isnothing(br.fBranchFlow)
        P_branch_out[k] += br.fBranchFlow.pFlow
        Q_branch_out[k] += br.fBranchFlow.qFlow
      end
    end
    if br.toBus in (b1, b2, b3)
      k = br.toBus == b1 ? 1 : (br.toBus == b2 ? 2 : 3)
      if !isnothing(br.tBranchFlow)
        P_branch_out[k] += br.tBranchFlow.pFlow
        Q_branch_out[k] += br.tBranchFlow.qFlow
      end
    end
  end

  bP = P_inj .- P_shunt .- P_branch_out
  bQ = Q_inj .- Q_shunt .- Q_branch_out
  # Remove any tiny component residual to match solver behavior.
  bP .-= sum(bP) / 3
  bQ .-= sum(bQ) / 3

  flowsP = [l.pFlow_MW for l in net.linkVec]
  flowsQ = [l.qFlow_MVar for l in net.linkVec]

  return vm_equal && va_equal && isapprox.(A * flowsP, bP; atol = 1e-7, rtol = 0.0) |> all && isapprox.(A * flowsQ, bQ; atol = 1e-7, rtol = 0.0) |> all
end

function test_run_sparlectra_normalizes_projected_matpower_starts()::Bool
  tmpdir = mktempdir()
  case_path = joinpath(tmpdir, "case_projected_start.m")
  write(
    case_path,
    """
function mpc = case_projected_start
mpc.version = '2';
mpc.baseMVA = 100;
mpc.bus = [
101 3 0 0 0 0 1 1.04 3.0 110 1 1.1 0.9;
205 2 10 5 0 0 1 1.03 2.0 110 1 1.1 0.9;
309 1 90 30 0 0 1 0.97 -4.0 110 1 1.1 0.9;
];
mpc.gen = [
101 100 0 300 -300 1.04 100 1 300 0;
205 20 0 100 -100 1.03 100 1 100 0;
];
mpc.branch = [
101 205 0.01 0.05 0.0 999 999 999 0 0 1 -360 360;
205 309 0.01 0.05 0.0 999 999 999 0 0 1 -360 360;
];
""",
  )

  cases = (
    (voltage_mode = :profile_blend, angle_mode = :dc, expected_flatstart = false),
    (voltage_mode = :classic, angle_mode = :classic, expected_flatstart = true),
    (voltage_mode = :classic, angle_mode = :matpower_va, expected_flatstart = false),
    (voltage_mode = :classic, angle_mode = :bus_va_blend, expected_flatstart = false),
    (voltage_mode = :all_bus_vm, angle_mode = :classic, expected_flatstart = false),
  )
  for case in cases
    cfg = Sparlectra.SparlectraConfig(
      powerflow = Sparlectra.PowerFlowConfig(start_mode = Sparlectra.StartModeConfig(flatstart = true, voltage_mode = case.voltage_mode, angle_mode = case.angle_mode), qlimits = Sparlectra.QLimitConfig(lock_pv_to_pq_buses = [205])),
      output = Sparlectra.OutputConfig(console_summary = false, logfile_results = :off),
    )
    ctx = Sparlectra._import_sparlectra_context(basename(case_path), dirname(case_path), cfg)
    effective_cfg = Sparlectra._resolve_matpower_powerflow_ids_after_import(ctx.net, ctx.config)
    effective_cfg.powerflow.start_mode.flatstart == case.expected_flatstart || return false
    ctx.projected_start_applied == !case.expected_flatstart || return false
    effective_cfg.powerflow.qlimits.lock_pv_to_pq_buses == [2] || return false
    cfg.powerflow.start_mode.flatstart == true || return false
    cfg.powerflow.qlimits.lock_pv_to_pq_buses == [205] || return false
  end
  # lock ids are MATPOWER bus numbers: an internal position passes, an
  # unknown id is dropped, duplicates collapse
  lock_cfg = Sparlectra.SparlectraConfig(powerflow = Sparlectra.PowerFlowConfig(qlimits = Sparlectra.QLimitConfig(lock_pv_to_pq_buses = [205])), output = Sparlectra.OutputConfig(console_summary = false, logfile_results = :off))
  imported_net = Sparlectra._import_sparlectra_net(basename(case_path), dirname(case_path), lock_cfg)
  resolve(ids) = Sparlectra._resolve_matpower_lock_pv_to_pq_buses(imported_net, ids)
  return resolve([205]) == [2] && resolve([2]) == [2] && resolve([999999]) == Int[] && resolve([999999, 205, 2, 205]) == [2]
end

function test_matpower_read_case_m_postprocessing()::Bool
  mfile = test_scratch_path(".m")
  open(mfile, "w") do io
    write(
      io,
      """
function mpc = case_postproc
mpc.version = '2';
mpc.baseMVA = 10;
mpc.bus = [
  1 3 1000 0 0 0 1 1.0 0.0 12.47 1 1.1 0.9;
  2 1 500 0 0 0 1 1.0 0.0 12.47 1 1.1 0.9;
];
mpc.gen = [
  1 0 0 100 -100 1.0 100 1 100 0 0 0 0 0 0 0 0 0 0 0 0;
];
mpc.branch = [
  1 2 0.0577 0.0409 0 0 0 0 0 0 1 -360 360;
];
mpc.gencost = [
  2 0 0 3 0 20 0;
];
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
Vbase = mpc.bus(1, BASE_KV) * 1e3;
Sbase = mpc.baseMVA * 1e6;
mpc.branch(:, [BR_R BR_X]) = mpc.branch(:, [BR_R BR_X]) / (Vbase^2 / Sbase);
mpc.bus(:, [PD, QD]) = mpc.bus(:, [PD, QD]) / 1e3;
pf = 0.85;
mpc.bus(:, QD) = mpc.bus(:, PD) * sin(acos(pf));
mpc.bus(:, PD) = mpc.bus(:, PD) * pf;
end
""",
    )
  end

  mpc = Sparlectra.MatpowerIO.read_case_m(mfile; legacy_compat = false)
  rm(mfile; force = true)

  pd_expected = 0.85
  qd_expected = 1.0 * sin(acos(0.85))

  zbase = (12.47e3)^2 / (10e6)
  r_expected = 0.0577 / zbase
  x_expected = 0.0409 / zbase

  ok_pd = isapprox(mpc.bus[1, 3], pd_expected; atol = 1e-12, rtol = 0.0)
  ok_qd = isapprox(mpc.bus[1, 4], qd_expected; atol = 1e-12, rtol = 0.0)
  ok_r = isapprox(mpc.branch[1, 3], r_expected; atol = 1e-12, rtol = 0.0)
  ok_x = isapprox(mpc.branch[1, 4], x_expected; atol = 1e-12, rtol = 0.0)

  if !(ok_pd && ok_qd && ok_r && ok_x)
    @warn "MATPOWER post-processing parse mismatch" pd = mpc.bus[1, 3] qd = mpc.bus[1, 4] r = mpc.branch[1, 3] x = mpc.branch[1, 4]
    return false
  end
  return true
end

function test_matpower_matrix_block_scanner_whitespace()::Bool
  variants = (
    """
mpc.bus = [
1 2 3;
4 5 6;
];
""",
    """
mpc.bus=[
1 2 3;
4 5 6;
];
""",
    """
mpc.bus    =    [
1 2 3;
4 5 6;
];
""",
  )

  expected = [1.0 2.0 3.0; 4.0 5.0 6.0]
  for txt in variants
    parsed = Sparlectra.MatpowerIO.parse_matrix_block(txt, "mpc.bus")
    parsed == expected || return false
  end
  return true
end

function test_matpower_matrix_block_scanner_large_body()::Bool
  row_count = 25000
  rows = join(("$(i) $(i + 1) $(i + 2);" for i = 1:row_count), '\n')
  txt = string("function mpc = synthetic_large\n", "mpc.version = '2';\n", "mpc.baseMVA = 100;\n", "mpc.bus = [\n", rows, "\n];\n", "mpc.gen = [\n1 0 0;\n];\n")

  parsed = Sparlectra.MatpowerIO.parse_matrix_block(txt, "mpc.bus"; ncols = 3)
  return size(parsed) == (row_count, 3) && parsed[1, :] == [1.0, 2.0, 3.0] && parsed[end, :] == [row_count, row_count + 1, row_count + 2]
end

function _assign_vset_controller!(net::Net, bus_name::String; vstep_pu::Union{Nothing,Float64}, tap_steps_down::Union{Nothing,Int}, tap_steps_up::Union{Nothing,Int}, all_generators::Bool = false)
  bus = geNetBusIdx(net = net, busName = bus_name)
  assigned = 0
  for ps in net.prosumpsVec
    isGenerator(ps) || continue
    getPosumerBusIndex(ps) == bus || continue
    ps.vstep_pu = vstep_pu
    ps.tap_steps_down = tap_steps_down
    ps.tap_steps_up = tap_steps_up
    assigned += 1
    all_generators || break
  end
  return assigned
end

# qlimit_mode on the 3-bus and 5-bus test nets: :adjust_vset keeps a bus
# with a voltage-step controller PV at a lowered setpoint, falls back to PQ
# without a controller or with no steps left, refuses two controllers on one
# bus and an incomplete step configuration; the default mode is
# :switch_to_pq.
function test_q_limit_adjust_vset_modes()::Bool
  station(net) = net.nodeVec[geNetBusIdx(net = net, busName = "STATION1")]
  # the controller set-up is checked before the first Newton step and
  # refused with error(); any other exception is a defect and propagates
  function refused(net)
    try
      runpf!(net, 20, 1e-8, 0; method = :rectangular, qlimit_mode = :adjust_vset)
    catch err
      err isa ErrorException || rethrow()
      return true
    end
    return false
  end

  net = createTest3BusNet(cooldown = 0, hyst_pu = 0.0, qlim_min = -30.0, qlim_max = 30.0)
  _assign_vset_controller!(net, "STATION1"; vstep_pu = 0.005, tap_steps_down = 2, tap_steps_up = 2)
  _, erg = runpf!(net, 40, 1e-8, 0; method = :rectangular, qlimit_mode = :adjust_vset, qlimit_max_outer = 20)
  adjusted_ok = erg == 0 && getNodeType(station(net)) == Sparlectra.PV && station(net)._vm_pu < 1.025

  net = createTest3BusNet(cooldown = 0, hyst_pu = 0.0, qlim_min = -15.0, qlim_max = 15.0)
  _, erg = runpf!(net, 30, 1e-8, 0; method = :rectangular, qlimit_mode = :adjust_vset)
  no_controller_ok = erg == 0 && getNodeType(station(net)) == Sparlectra.PQ

  net = createTest3BusNet(cooldown = 0, hyst_pu = 0.0, qlim_min = -15.0, qlim_max = 15.0)
  _assign_vset_controller!(net, "STATION1"; vstep_pu = 0.005, tap_steps_down = 0, tap_steps_up = 0)
  _, erg = runpf!(net, 30, 1e-8, 0; method = :rectangular, qlimit_mode = :adjust_vset, qlimit_max_outer = 1)
  exhausted_ok = erg == 0 && getNodeType(station(net)) == Sparlectra.PQ

  net = createTest5BusNet(cooldown = 0, hyst_pu = 0.0, qlim_min = -15.0, qlim_max = 15.0, mul_gens = true)
  _assign_vset_controller!(net, "B3"; vstep_pu = 0.005, tap_steps_down = 2, tap_steps_up = 2, all_generators = true)
  two_controllers_refused = refused(net)

  net = createTest3BusNet(cooldown = 0, hyst_pu = 0.0, qlim_min = -20.0, qlim_max = 20.0)
  _assign_vset_controller!(net, "STATION1"; vstep_pu = 0.005, tap_steps_down = nothing, tap_steps_up = nothing)
  incomplete_refused = refused(net)

  net = createTest5BusNet(cooldown = 0, hyst_pu = 0.0, qlim_min = -15.0, qlim_max = 15.0, mul_gens = true)
  _assign_vset_controller!(net, "B3"; vstep_pu = 0.005, tap_steps_down = 1, tap_steps_up = 1, all_generators = false)
  _, erg = runpf!(net, 40, 1e-8, 0; method = :rectangular, qlimit_mode = :adjust_vset)
  one_of_several_ok = erg == 0

  net_default = createTest3BusNet(cooldown = 0, hyst_pu = 0.0, qlim_min = -15.0, qlim_max = 15.0)
  net_explicit = createTest3BusNet(cooldown = 0, hyst_pu = 0.0, qlim_min = -15.0, qlim_max = 15.0)
  _, erg_default = runpf!(net_default, 30, 1e-8, 0; method = :rectangular)
  _, erg_explicit = runpf!(net_explicit, 30, 1e-8, 0; method = :rectangular, qlimit_mode = :switch_to_pq)
  default_ok = erg_default == 0 && erg_explicit == 0 && getNodeType(station(net_default)) == getNodeType(station(net_explicit)) && isapprox(station(net_default)._vm_pu, station(net_explicit)._vm_pu; atol = 1e-8, rtol = 0.0)

  return adjusted_ok && no_controller_ok && exhausted_ok && two_controllers_refused && incomplete_refused && one_of_several_ok && default_ok
end

# The classic Q-limit modes: with a base power flow that does not converge
# classic_simultaneous switches nothing and says why; classic_one_at_a_time
# keeps the log entry of its switching round (the inner solves used to reset
# net.qLimitLog every outer round, so a real PV->PQ switch left an empty log).
function test_q_limit_matpower_mode_base_failure_does_not_switch()::Bool
  net = createTest3BusNet(cooldown = 0, hyst_pu = 0.0, qlim_min = -15.0, qlim_max = 15.0)
  bus = geNetBusIdx(net = net, busName = "STATION1")
  _, erg = runpf!(net, 1, 1e-14, 0; method = :rectangular, qlimit_enforcement_mode = :classic_simultaneous, qlimit_max_outer = 5)
  st = Sparlectra.rectangular_pf_status(net)
  base_failure_ok =
    erg == 1 &&
    getNodeType(net.nodeVec[bus]) == Sparlectra.PV &&
    isempty(net.qLimitLog) &&
    !isnothing(st) &&
    st.qlimit_enforcement_mode == :classic_simultaneous &&
    st.base_pf_converged === false &&
    st.qlimit_enforcement_started === false &&
    st.final_outcome == :base_pf_not_converged

  net = createTest3BusNet(cooldown = 2, hyst_pu = 0.01, qlim_min = -15.0, qlim_max = 15.0)
  runpf!(net, 30, 1e-8, 0; method = :rectangular, qlimit_enforcement_mode = :classic_one_at_a_time, qlimit_max_outer = 5)
  st = Sparlectra.rectangular_pf_status(net)
  one_at_a_time_ok = !isnothing(st) && st.qlimit_enforcement_mode == :classic_one_at_a_time && st.qlimit_reenable_events == 0 && length(net.qLimitLog) == 1 && net.qLimitLog[1].side == :max

  return base_failure_ok && one_at_a_time_ok
end

# #375: the PQ->PV release decides on the voltage side once the solver hands
# the voltages and setpoints in (rectangular_qlimit_iteration.jl does): a bus
# clamped at Qmax goes back to PV when Vm > Vset + v_hyst_pu, one clamped at
# Qmin when Vm < Vset - v_hyst_pu. The reactive injection sits ON the clamp
# in every call below, so the old Q band test would never release; a release
# proves the voltage rule. Every case starts from a fresh clamp because the
# one-retry guard keeps a bus clamped after its second hit.
function test_q_limit_reenable_voltage_rule()::Bool
  nb = 1
  qmin_pu = [0.0]
  qmax_pu = [1.0]
  pv_orig_mask = trues(nb)
  # a fresh net with bus 1 clamped at `side` (:max by qreq 1.11 over Qmax 1.0,
  # :min by qreq -0.05 under Qmin 0.0); returns the net and its bus types
  function clamped(side::Symbol)
    net = Net(name = "q_limit_voltage_rule", baseMVA = 100.0)
    bus_types = [:PV]
    qreq = side == :max ? 1.11 : -0.05
    changed, _ = Sparlectra.active_set_q_limits!(
      net, 1, nb;
      get_qreq_pu = _ -> qreq,
      is_pv = bus -> (bus_types[bus] == :PV),
      make_pq! = (bus, _qclamp, _side) -> (bus_types[bus] = :PQ),
      make_pv! = bus -> (bus_types[bus] = :PV),
      qmin_pu = qmin_pu, qmax_pu = qmax_pu, pv_orig_mask = pv_orig_mask,
      allow_reenable = false, q_hyst_pu = 0.0, cooldown_iters = 0,
    )
    ok = changed && bus_types[1] == :PQ && get(net.qLimitEvents, 1, nothing) === side && length(net.qLimitLog) == 1
    return net, bus_types, ok
  end
  # the release call of the next iteration: the injection sits on the clamp,
  # only the voltages can free the bus
  function release_step(net, bus_types; vm::Float64, vset::Float64 = 1.0, v_hyst_pu::Float64 = 1e-4, voltages::Bool = true, qreq::Float64 = bus_types[1] == :PQ && get(net.qLimitEvents, 1, :max) === :max ? 1.0 : 0.0)
    changed, reenabled = Sparlectra.active_set_q_limits!(
      net, 2, nb;
      get_qreq_pu = _ -> qreq,
      is_pv = bus -> (bus_types[bus] == :PV),
      make_pq! = (bus, _qclamp, _side) -> (bus_types[bus] = :PQ),
      make_pv! = bus -> (bus_types[bus] = :PV),
      qmin_pu = qmin_pu, qmax_pu = qmax_pu, pv_orig_mask = pv_orig_mask,
      allow_reenable = true, q_hyst_pu = 0.0, cooldown_iters = 0,
      get_vm_pu = voltages ? (_ -> vm) : nothing,
      get_vset_pu = voltages ? (_ -> vset) : nothing,
      v_hyst_pu = v_hyst_pu,
    )
    return (changed = changed, reenabled = reenabled, pv = bus_types[1] == :PV, events = length(net.qLimitEvents))
  end
  released(r) = !r.changed && r.reenabled && r.pv && r.events == 0
  held(r) = !r.changed && !r.reenabled && !r.pv && r.events == 1

  # 1 and 2: clamped at Qmax, voltage above the setpoint -> released
  net, bt, ok1 = clamped(:max)
  ok2 = released(release_step(net, bt; vm = 1.02))
  # 3: at Qmax with the voltage below the setpoint the clamp is right, no release
  net, bt, ok3a = clamped(:max)
  ok3 = ok3a && held(release_step(net, bt; vm = 0.99))
  # 4: mirror image at Qmin: released below the setpoint, held above it
  net, bt, ok4a = clamped(:min)
  ok4 = ok4a && released(release_step(net, bt; vm = 0.98))
  net, bt, ok4b = clamped(:min)
  ok4 = ok4 && ok4b && held(release_step(net, bt; vm = 1.02))
  # 5: the margin: 5e-5 above the setpoint is inside 1e-4, outside 1e-5
  net, bt, ok5a = clamped(:max)
  ok5 = ok5a && held(release_step(net, bt; vm = 1.00005, v_hyst_pu = 1e-4))
  net, bt, ok5b = clamped(:max)
  ok5 = ok5 && ok5b && released(release_step(net, bt; vm = 1.00005, v_hyst_pu = 1e-5))
  # 6: without voltages the Q band test decides as before: an injection on
  # the clamp holds, one inside the band releases
  net, bt, ok6a = clamped(:max)
  ok6 = ok6a && held(release_step(net, bt; vm = 1.02, voltages = false, qreq = 1.0))
  net, bt, ok6b = clamped(:max)
  ok6 = ok6 && ok6b && released(release_step(net, bt; vm = 0.99, voltages = false, qreq = 0.5))

  return ok1 && ok2 && ok3 && ok4 && ok5 && ok6
end

# Q-limit start modes on the 3-bus net: qlimit_start_iter 10 with three
# iterations switches nothing; qlimit_start_mode auto with an unbounded
# Q delta switches STATION1 to PQ.
function test_q_limit_start_iter_delays_switching()::Bool
  net = createTest3BusNet(cooldown = 0, hyst_pu = 0.0, qlim_min = -15.0, qlim_max = 15.0)
  _, erg = runpf!(net, 3, 1e-12, 0; method = :rectangular, qlimit_start_iter = 10)
  bus = geNetBusIdx(net = net, busName = "STATION1")
  delayed_ok = erg == 1 && getNodeType(net.nodeVec[bus]) == Sparlectra.PV && isempty(net.qLimitLog)

  net = createTest3BusNet(cooldown = 0, hyst_pu = 0.0, qlim_min = -15.0, qlim_max = 15.0)
  _, erg = runpf!(net, 30, 1e-8, 0; method = :rectangular, qlimit_start_mode = :auto, qlimit_auto_q_delta_pu = Inf)
  auto_ok = erg == 0 && getNodeType(net.nodeVec[bus]) == Sparlectra.PQ && !isempty(net.qLimitLog)

  return delayed_ok && auto_ok
end

function test_q_limit_hysteresis_delays_small_pv_to_pq_overshoot()::Bool
  net = Net(name = "q_limit_hysteresis_deadband", baseMVA = 100.0)
  nb = 1
  qmin_pu = [-1.0]
  qmax_pu = [1.0]
  bus_types = [:PV]
  pv_orig_mask = trues(nb)

  changed, _ = Sparlectra.active_set_q_limits!(
    net,
    1,
    nb;
    get_qreq_pu = _ -> 1.05,
    is_pv = bus -> (bus_types[bus] == :PV),
    make_pq! = (bus, _qclamp, _side) -> (bus_types[bus] = :PQ),
    make_pv! = bus -> (bus_types[bus] = :PV),
    qmin_pu = qmin_pu,
    qmax_pu = qmax_pu,
    pv_orig_mask = pv_orig_mask,
    allow_reenable = false,
    q_hyst_pu = 0.10,
    cooldown_iters = 0,
  )
  delayed = !changed && bus_types[1] == :PV && isempty(net.qLimitLog)

  changed, _ = Sparlectra.active_set_q_limits!(
    net,
    2,
    nb;
    get_qreq_pu = _ -> 1.11,
    is_pv = bus -> (bus_types[bus] == :PV),
    make_pq! = (bus, _qclamp, _side) -> (bus_types[bus] = :PQ),
    make_pv! = bus -> (bus_types[bus] = :PV),
    qmin_pu = qmin_pu,
    qmax_pu = qmax_pu,
    pv_orig_mask = pv_orig_mask,
    allow_reenable = false,
    q_hyst_pu = 0.10,
    cooldown_iters = 0,
  )

  first_switch = changed && bus_types[1] == :PQ && length(net.qLimitLog) == 1

  changed, reenabled = Sparlectra.active_set_q_limits!(
    net,
    3,
    nb;
    get_qreq_pu = _ -> 0.0,
    is_pv = bus -> (bus_types[bus] == :PV),
    make_pq! = (bus, _qclamp, _side) -> (bus_types[bus] = :PQ),
    make_pv! = bus -> (bus_types[bus] = :PV),
    qmin_pu = qmin_pu,
    qmax_pu = qmax_pu,
    pv_orig_mask = pv_orig_mask,
    allow_reenable = true,
    q_hyst_pu = 0.10,
    cooldown_iters = 0,
  )
  first_reenable = !changed && reenabled && bus_types[1] == :PV && isempty(net.qLimitEvents)
  # the release is in the log, with the side the bus was clamped at
  first_reenable &= [(e.iter, e.side, e.kind) for e in net.qLimitLog] == [(2, :max, :hit), (3, :max, :release)]

  changed, _ = Sparlectra.active_set_q_limits!(
    net,
    4,
    nb;
    get_qreq_pu = _ -> 1.11,
    is_pv = bus -> (bus_types[bus] == :PV),
    make_pq! = (bus, _qclamp, _side) -> (bus_types[bus] = :PQ),
    make_pv! = bus -> (bus_types[bus] = :PV),
    qmin_pu = qmin_pu,
    qmax_pu = qmax_pu,
    pv_orig_mask = pv_orig_mask,
    allow_reenable = false,
    q_hyst_pu = 0.10,
    cooldown_iters = 0,
  )
  second_switch = changed && bus_types[1] == :PQ && Sparlectra.qlimit_pv_to_pq_event_count(net) == 2

  changed, reenabled = Sparlectra.active_set_q_limits!(
    net,
    5,
    nb;
    get_qreq_pu = _ -> 0.0,
    is_pv = bus -> (bus_types[bus] == :PV),
    make_pq! = (bus, _qclamp, _side) -> (bus_types[bus] = :PQ),
    make_pv! = bus -> (bus_types[bus] = :PV),
    qmin_pu = qmin_pu,
    qmax_pu = qmax_pu,
    pv_orig_mask = pv_orig_mask,
    allow_reenable = true,
    q_hyst_pu = 0.10,
    cooldown_iters = 0,
  )
  second_reenable_blocked = !changed && !reenabled && bus_types[1] == :PQ && length(net.qLimitLog) == 3

  # #475: on a converged iterate the twice-clamped bus gets the joint extra
  # release; when it violates again it goes back as a group revert (logged,
  # but not counted as a hit), and the solve makes no further extra release
  reenable_call(it; converged) = Sparlectra.active_set_q_limits!(
    net,
    it,
    nb;
    get_qreq_pu = _ -> 0.0,
    is_pv = bus -> (bus_types[bus] == :PV),
    make_pq! = (bus, _qclamp, _side) -> (bus_types[bus] = :PQ),
    make_pv! = bus -> (bus_types[bus] = :PV),
    qmin_pu = qmin_pu,
    qmax_pu = qmax_pu,
    pv_orig_mask = pv_orig_mask,
    allow_reenable = true,
    q_hyst_pu = 0.10,
    cooldown_iters = 0,
    converged_iterate = converged,
  )
  changed, reenabled = reenable_call(6; converged = true)
  converged_reenable = !changed && reenabled && bus_types[1] == :PV
  changed, _ = Sparlectra.active_set_q_limits!(
    net,
    7,
    nb;
    get_qreq_pu = _ -> 1.11,
    is_pv = bus -> (bus_types[bus] == :PV),
    make_pq! = (bus, _qclamp, _side) -> (bus_types[bus] = :PQ),
    make_pv! = bus -> (bus_types[bus] = :PV),
    qmin_pu = qmin_pu,
    qmax_pu = qmax_pu,
    pv_orig_mask = pv_orig_mask,
    allow_reenable = false,
    q_hyst_pu = 0.10,
    cooldown_iters = 0,
  )
  group_revert = changed && bus_types[1] == :PQ && net.qLimitLog[end].kind === :group_revert && Sparlectra.qlimit_switch_count(net, 1) == 2
  changed, reenabled = reenable_call(8; converged = true)
  no_second_extra = !changed && !reenabled && bus_types[1] == :PQ

  # the guard in lock_pq mode switches an overshoot inside the hysteresis
  # band (1.05 against Qmax 1.0 and q_hyst_pu 0.10) at once
  guard_net = Net(name = "q_limit_violation_guard", baseMVA = 100.0)
  guard_types = [:PV]
  changed, reenabled = Sparlectra.active_set_q_limits!(
    guard_net,
    2,
    nb;
    get_qreq_pu = _ -> 1.05,
    is_pv = bus -> (guard_types[bus] == :PV),
    make_pq! = (bus, _qclamp, _side) -> (guard_types[bus] = :PQ),
    make_pv! = bus -> (guard_types[bus] = :PV),
    qmin_pu = [-1.0],
    qmax_pu = [1.0],
    pv_orig_mask = trues(nb),
    allow_reenable = false,
    q_hyst_pu = 0.10,
    cooldown_iters = 0,
    qlimit_guard_violation_mode = :lock_pq,
    qlimit_guard_violation_threshold_pu = 1e-4,
  )
  guard_switch = changed && !reenabled && guard_types[1] == :PQ && length(guard_net.qLimitLog) == 1 && guard_net.qLimitLog[1].side == :max

  return delayed && first_switch && first_reenable && second_switch && second_reenable_blocked && converged_reenable && group_revert && no_second_extra && guard_switch
end


function test_solve_sparse_system_small_and_large_backends()::Bool
  smallA = sparse([1.0 1.0; 2.0 2.0])
  smallb = [1.0, 2.0]
  small_diag = Dict{Symbol,Any}(:enabled => true)
  smallx = Sparlectra.solve_sparse_system(smallA, smallb; context = :test_small, diagnostics = small_diag)
  small_ok = all(isfinite, smallx) && isapprox(smallA * smallx, smallb; atol = 1e-10, rtol = 0.0) && get(small_diag, :backend, :missing) in (:sparse_backslash, :sparse_qr_fallback, :dense_svd_fallback_small_system)

  n = Sparlectra.DEFAULT_DENSE_FALLBACK_MAX_N + 4
  largeA = spdiagm(0 => fill(4.0, n), 1 => fill(-1.0, n - 1), -1 => fill(-1.0, n - 1))
  largeb = ones(Float64, n)
  large_diag = Dict{Symbol,Any}(:enabled => true)
  largex = Sparlectra.solve_sparse_system(largeA, largeb; context = :test_large, diagnostics = large_diag)
  large_ok = all(isfinite, largex) && isapprox(largeA * largex, largeb; atol = 1e-10, rtol = 1e-10) && get(large_diag, :backend, :missing) != :dense_svd_fallback_small_system

  return small_ok && large_ok
end

function test_rectangular_jacobian_and_dc_start_matrices_are_sparse()::Bool
  Ybus = sparse([1, 1, 2, 2, 2, 3, 3], [1, 2, 1, 2, 3, 2, 3], ComplexF64[8im, -8im, -8im, 14im, -6im, -6im, 6im], 3, 3)
  V = ComplexF64[1.0+0im, 1.0+0im, 1.0+0im]
  S = ComplexF64[0.0+0im, -0.3-0.1im, -0.2-0.08im]
  bus_types = [:Slack, :PQ, :PQ]
  Vset = [1.0, 1.0, 1.0]
  slack_idx = 1
  J = Sparlectra.build_rectangular_jacobian_pq_pv(Ybus, V, bus_types, Vset, slack_idx)
  jac_ok = J isa SparseMatrixCSC

  perf = Dict{Symbol,Any}(:enabled => true)
  _ = Sparlectra._dc_angle_start_rectangular_profiled(Ybus, V, S, bus_types, Vset, slack_idx; performance_profile = perf)
  dc_ok = get(perf, :dc_matrix_is_sparse, false) == true

  return dc_ok && jac_ok
end

function test_rectangular_singular_step_returns_nonconvergence()::Bool
  net = Net(name = "singular_rectangular_step", baseMVA = 100.0)
  addBus!(net = net, busName = "Slack", vn_kV = 110.0)
  addBus!(net = net, busName = "Load", vn_kV = 110.0)
  addProsumer!(net = net, busName = "Slack", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "Slack")
  addProsumer!(net = net, busName = "Load", type = "ENERGYCONSUMER", p = 10.0, q = 5.0)
  _, erg = runpf!(net, 5, 1e-8, 0; method = :rectangular)
  return erg == 1
end

function test_bus_type_resolution_from_prosumers()::Bool
  net = Net(name = "bus_type_resolution", baseMVA = 100.0)
  addBus!(net = net, busName = "LoadBus", vn_kV = 110.0)
  addBus!(net = net, busName = "PVBus", vn_kV = 110.0)
  addBus!(net = net, busName = "GenPQBus", vn_kV = 110.0)
  addBus!(net = net, busName = "MixBus", vn_kV = 110.0)
  addBus!(net = net, busName = "SlackBus", vn_kV = 110.0)
  addBus!(net = net, busName = "NoProsumerBus", vn_kV = 110.0)

  addProsumer!(net = net, busName = "LoadBus", type = "ENERGYCONSUMER", p = 10.0, q = 3.0)
  addProsumer!(net = net, busName = "PVBus", type = "SYNCHRONOUSMACHINE", p = 20.0, vm_pu = 1.02, isRegulated = true)
  addProsumer!(net = net, busName = "GenPQBus", type = "SYNCHRONOUSMACHINE", p = 15.0, q = 2.0)
  addProsumer!(net = net, busName = "MixBus", type = "SYNCHRONOUSMACHINE", p = 12.0, q = 1.0)
  addProsumer!(net = net, busName = "MixBus", type = "SYNCHRONOUSMACHINE", p = 8.0, vm_pu = 1.01, isRegulated = true)
  addProsumer!(net = net, busName = "SlackBus", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "SlackBus")
  addProsumer!(net = net, busName = "NoProsumerBus", type = "ENERGYCONSUMER", p = 5.0, q = 1.0)

  pv_bus = geNetBusIdx(net = net, busName = "PVBus")
  genpq_bus = geNetBusIdx(net = net, busName = "GenPQBus")
  pv_regulation_ok = any(isGenerator(ps) && isRegulating(ps) for ps in getBusProsumers(net, pv_bus))
  genpq_regulation_ok = all(!(isGenerator(ps) && isRegulating(ps)) for ps in getBusProsumers(net, genpq_bus))

  ok_types =
    getEffectiveBusType(net = net, busName = "LoadBus") == Sparlectra.PQ &&
    getEffectiveBusType(net = net, busName = "PVBus") == Sparlectra.PV &&
    getEffectiveBusType(net = net, busName = "GenPQBus") == Sparlectra.PQ &&
    getEffectiveBusType(net = net, busName = "MixBus") == Sparlectra.PV &&
    getEffectiveBusType(net = net, busName = "SlackBus") == Sparlectra.Slack &&
    getEffectiveBusType(net = net, busName = "NoProsumerBus") == Sparlectra.PQ

  # two slack prosumers on one bus (a MATPOWER bus with several generators)
  # give one Slack bus
  net2 = Net(name = "multi_slack_same_bus", baseMVA = 100.0)
  addBus!(net = net2, busName = "SlackBus", vn_kV = 110.0)
  addProsumer!(net = net2, busName = "SlackBus", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "SlackBus")
  addProsumer!(net = net2, busName = "SlackBus", type = "GENERATOR", p = 15.0, q = 2.0, referencePri = "SlackBus")
  two_slack_bus = geNetBusIdx(net = net2, busName = "SlackBus")
  two_slack_ok = getEffectiveBusType(net2, two_slack_bus) == Sparlectra.Slack && count(isSlack, getBusProsumers(net2, two_slack_bus)) == 2

  # addProsumer! marks a machine regulated when it gets vm_pu or tap steps,
  # and its bus becomes PV
  addBus!(net = net, busName = "PVVm", vn_kV = 110.0)
  addBus!(net = net, busName = "PVTap", vn_kV = 110.0)
  addBus!(net = net, busName = "MachPQ", vn_kV = 110.0)
  addProsumer!(net = net, busName = "PVVm", type = "SYNCHRONOUSMACHINE", p = 10.0, vm_pu = 1.01)
  addProsumer!(net = net, busName = "PVTap", type = "SYNCHRONOUSMACHINE", p = 8.0, vstep_pu = 0.005, tap_steps_down = 2, tap_steps_up = 2)
  addProsumer!(net = net, busName = "MachPQ", type = "SYNCHRONOUSMACHINE", p = 6.0, q = 1.5)
  ps_vm = only(getBusProsumers(net, geNetBusIdx(net = net, busName = "PVVm")))
  ps_tap = only(getBusProsumers(net, geNetBusIdx(net = net, busName = "PVTap")))
  ps_pq = only(getBusProsumers(net, geNetBusIdx(net = net, busName = "MachPQ")))
  auto_flags_ok =
    ps_vm.isRegulated && ps_tap.isRegulated && !ps_pq.isRegulated &&
    getEffectiveBusType(net = net, busName = "PVVm") == Sparlectra.PV &&
    getEffectiveBusType(net = net, busName = "PVTap") == Sparlectra.PV &&
    getEffectiveBusType(net = net, busName = "MachPQ") == Sparlectra.PQ

  return ok_types && pv_regulation_ok && genpq_regulation_ok && two_slack_ok && auto_flags_ok
end

function test_regulated_generator_bus_targets_include_unregulated_generators()::Bool
  net = Net(name = "regulated_generator_bus_targets", baseMVA = 100.0)
  addBus!(net = net, busName = "Slack", vn_kV = 110.0)
  addBus!(net = net, busName = "PVBus", vn_kV = 110.0)

  addProsumer!(net = net, busName = "Slack", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "Slack")

  # Regulated generator defines PV behavior at the bus.
  addProsumer!(net = net, busName = "PVBus", type = "SYNCHRONOUSMACHINE", p = 20.0, q = 3.0, vm_pu = 1.02, isRegulated = true, qMin = -30.0, qMax = 40.0)
  # Additional unregulated generator at the same bus.
  addProsumer!(net = net, busName = "PVBus", type = "SYNCHRONOUSMACHINE", p = 7.5, q = 1.0, qMin = -10.0, qMax = 15.0)

  # The Q band of a PV bus bounds its total generator Q. An unregulated unit
  # runs at its fixed Q: it enters both ends with that Q, its own limits do
  # not widen the band, and a missing limit does not open it (PVFix). A
  # regulating unit without limits opens it, whichever unit comes first
  # (PVOrder). A load at the bus does not count. A unit with a Q(U)
  # characteristic is fixed at its Q at the held voltage (PVQU): 1.02 pu on
  # the line from (0.95, 0.10) to (1.05, -0.10) gives -0.04 pu, neither its
  # scheduled q nor its own limits count.
  addBus!(net = net, busName = "PVFix", vn_kV = 110.0)
  addProsumer!(net = net, busName = "PVFix", type = "SYNCHRONOUSMACHINE", p = 20.0, q = 3.0, vm_pu = 1.02, isRegulated = true, qMin = -30.0, qMax = 40.0)
  addProsumer!(net = net, busName = "PVFix", type = "GENERATOR", p = 5.0, q = -2.0)
  addProsumer!(net = net, busName = "PVFix", type = "ENERGYCONSUMER", p = 4.0, q = 1.0)
  addBus!(net = net, busName = "PVOrder", vn_kV = 110.0)
  addProsumer!(net = net, busName = "PVOrder", type = "SYNCHRONOUSMACHINE", p = 10.0, vm_pu = 1.01, isRegulated = true)
  addProsumer!(net = net, busName = "PVOrder", type = "SYNCHRONOUSMACHINE", p = 10.0, vm_pu = 1.01, isRegulated = true, qMin = -20.0, qMax = 20.0)
  addBus!(net = net, busName = "PVQU", vn_kV = 110.0)
  addProsumer!(net = net, busName = "PVQU", type = "SYNCHRONOUSMACHINE", p = 20.0, q = 3.0, vm_pu = 1.02, isRegulated = true, qMin = -30.0, qMax = 40.0)
  qu_curve = make_characteristic([(0.95, 0.10), (1.05, -0.10)])
  addProsumer!(net = net, busName = "PVQU", type = "GENERATOR", p = 5.0, q = 2.0, qMin = -50.0, qMax = 50.0, qu_controller = QUController(qu_curve; qmin_MVAr = -50.0, qmax_MVAr = 50.0, sbase_MVA = 100.0))

  pv_bus = geNetBusIdx(net = net, busName = "PVBus")
  fix_bus = geNetBusIdx(net = net, busName = "PVFix")
  order_bus = geNetBusIdx(net = net, busName = "PVOrder")
  qu_bus = geNetBusIdx(net = net, busName = "PVQU")
  node = net.nodeVec[pv_bus]
  qmin_pu, qmax_pu = getQLimits_pu(net)

  p_target_ok = isapprox(node._pƩGen, 27.5; atol = 1e-9, rtol = 0.0)
  q_target_ok = isapprox(node._qƩGen, 4.0; atol = 1e-9, rtol = 0.0)
  # (-30 + 1) / 100 and (40 + 1) / 100
  q_limit_ok = isapprox(qmin_pu[pv_bus], -0.29; atol = 1e-9, rtol = 0.0) && isapprox(qmax_pu[pv_bus], 0.41; atol = 1e-9, rtol = 0.0)
  # (-30 - 2) / 100 and (40 - 2) / 100
  fix_ok = isapprox(qmin_pu[fix_bus], -0.32; atol = 1e-9, rtol = 0.0) && isapprox(qmax_pu[fix_bus], 0.38; atol = 1e-9, rtol = 0.0)
  order_ok = qmin_pu[order_bus] == -Inf && qmax_pu[order_bus] == Inf
  # -0.30 - 0.04 and 0.40 - 0.04
  qu_ok = isapprox(qmin_pu[qu_bus], -0.34; atol = 1e-9, rtol = 0.0) && isapprox(qmax_pu[qu_bus], 0.36; atol = 1e-9, rtol = 0.0)
  bus_type_ok = getEffectiveBusType(net = net, busName = "PVBus") == Sparlectra.PV

  return p_target_ok && q_target_ok && q_limit_ok && fix_ok && order_ok && qu_ok && bus_type_ok
end

# The summary mode of printACPFlowResults returns early and closes the
# result file itself; without that close the file stays unflushed and empty.
# The Q-limit lines carry the counts of qLimitEvents and qLimitLog.
function test_summary_result_output_closes_file_handle()::Bool
  net = Net(name = "pv2pq_label_smoke", baseMVA = 100.0)
  addBus!(net = net, busName = "Slack", vn_kV = 110.0)
  addBus!(net = net, busName = "Load", vn_kV = 110.0)
  addPIModelACLine!(net = net, fromBus = "Slack", toBus = "Load", r_pu = 0.01, x_pu = 0.05, b_pu = 0.0, status = 1)
  addProsumer!(net = net, busName = "Slack", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "Slack")
  addProsumer!(net = net, busName = "Load", type = "ENERGYCONSUMER", p = 1.0, q = 0.3)
  net.qLimitEvents[2] = :max
  push!(net.qLimitLog, (it = 1, bus = 2, q = 0.4, qmin = -0.2, qmax = 0.2, action = "switch_to_pq"))

  outdir = mktempdir()
  Sparlectra.printACPFlowResults(net, 0.01, 2, 1e-6, true, outdir; converged = true, solver = :rectangular, solver_time_s = 0.008, result_mode = :summary, max_rows = 5)
  result_file = joinpath(outdir, "result_$(net.name).txt")
  if !isfile(result_file)
    return false
  end
  out = read(result_file, String)
  return !isempty(out) && occursin(r"PV→PQ locks\s*:\s*1\n", out) && occursin(r"PV→PQ events\s*:\s*1\n", out)
end

# Condition-number diagnostics (condition_number.jl): the Hager estimate
# must match exactly computable cases, stay within its documented factor of
# the true 1-norm condition number, flag near-singular systems, and reject
# invalid input.
function test_condition_number_estimator()::Bool
  # diagonal case: kappa_1 is exactly max|d| * max|1/d|
  D = Matrix(Diagonal([1.0, 2.0, 4.0]))
  isapprox(condestJacobian(D), 4.0; rtol = 1e-12) || return false
  isapprox(condestJacobian(D; exact = true), 4.0; rtol = 1e-12) || return false

  # sparse diagonally dominant system: the estimate is a lower bound on the
  # true 1-norm condition number and stays within its documented factor
  A = spdiagm(0 => fill(4.0, 40), 1 => fill(-1.0, 39), -1 => fill(-1.0, 39))
  est = condestJacobian(A)
  exact1 = opnorm(Matrix(A), 1) * opnorm(inv(Matrix(A)), 1)
  (est <= exact1 * (1 + 1e-9) && est >= exact1 / 3) || return false

  # near-singular system must be flagged as such
  S = Matrix(Diagonal([1.0, 1.0, 1e-13]))
  condestJacobian(S) > 1e12 || return false

  # reportCondition returns the same estimate it prints
  isapprox(reportCondition(D), condestJacobian(D); rtol = 1e-12) || return false

  # invalid input is rejected
  try
    condestJacobian(zeros(2, 3))
    return false
  catch err
    err isa ArgumentError || return false
  end
  try
    condestJacobian(zeros(0, 0))
    return false
  catch err
    err isa ArgumentError || return false
  end

  # net method: solver Jacobian at the current voltage state of a solved
  # small net gives a finite, plausible estimate
  net = Net(name = "cond_check", baseMVA = 100.0)
  addBus!(net = net, busName = "Slack", vn_kV = 110.0)
  addBus!(net = net, busName = "Load", vn_kV = 110.0)
  addPIModelACLine!(net = net, fromBus = "Slack", toBus = "Load", r_pu = 0.01, x_pu = 0.05, b_pu = 0.0, status = 1)
  addProsumer!(net = net, busName = "Slack", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "Slack")
  addProsumer!(net = net, busName = "Load", type = "ENERGYCONSUMER", p = 10.0, q = 3.0)
  ite, erg = runpf!(net, 20, 1e-8, 0; method = :rectangular)
  erg == 0 || return false
  knet = condestJacobian(net)
  (isfinite(knet) && knet >= 1.0) || return false

  # solver-status path: the rectangular solve stores a lazy, cached
  # condition-estimate closure; repeated calls return the same value
  # without refactorizing
  st = Sparlectra.rectangular_pf_status(net)
  st !== nothing || return false
  hasproperty(st, :jacobian_condest) || return false
  kthunk = st.jacobian_condest()
  (isfinite(kthunk) && kthunk >= 1.0) || return false
  st.jacobian_condest() == kthunk || return false

  # the full result log always carries the condition line
  outdir = mktempdir()
  Sparlectra.printACPFlowResults(net, 0.01, ite, 1e-8, true, outdir; converged = true, solver = :rectangular)
  out_on = read(joinpath(outdir, "result_$(net.name).txt"), String)
  occursin("Jacobian cond. : kappa1(J)", out_on) || return false

  # overview availability rule: without a stored rectangular thunk (unsolved
  # net, APSLF, DC) the metadata path reports nothing instead of
  # reconstructing a Jacobian the run never factored
  net_unsolved = Net(name = "cond_no_thunk", baseMVA = 100.0)
  addBus!(net = net_unsolved, busName = "Slack", vn_kV = 110.0)
  addBus!(net = net_unsolved, busName = "Load", vn_kV = 110.0)
  addPIModelACLine!(net = net_unsolved, fromBus = "Slack", toBus = "Load", r_pu = 0.01, x_pu = 0.05, b_pu = 0.0, status = 1)
  addProsumer!(net = net_unsolved, busName = "Slack", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "Slack")
  addProsumer!(net = net_unsolved, busName = "Load", type = "ENERGYCONSUMER", p = 10.0, q = 3.0)
  Sparlectra._jacobian_condest(net_unsolved; warn_on_failure = false, require_thunk = true) === nothing || return false

  # regression: de-energized buses reduce the Ybus dimension below the bus
  # count. The net-level estimate used to fail with a DimensionMismatch
  # (reduced Ybus against full-length state vectors); it must now describe
  # the active subsystem, and the result log must carry the line instead
  # of the skip warning.
  net_iso = Net(name = "cond_check_iso", baseMVA = 100.0)
  addBus!(net = net_iso, busName = "Slack", vn_kV = 110.0)
  addBus!(net = net_iso, busName = "Load", vn_kV = 110.0)
  addBus!(net = net_iso, busName = "Dead", vn_kV = 110.0)
  addPIModelACLine!(net = net_iso, fromBus = "Slack", toBus = "Load", r_pu = 0.01, x_pu = 0.05, b_pu = 0.0, status = 1)
  addProsumer!(net = net_iso, busName = "Slack", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "Slack")
  addProsumer!(net = net_iso, busName = "Load", type = "ENERGYCONSUMER", p = 10.0, q = 3.0)
  markIsolatedBuses!(net = net_iso, log = false)
  length(net_iso.isoNodes) == 1 || return false
  ite_iso, erg_iso = runpf!(net_iso, 20, 1e-8, 0; method = :rectangular)
  erg_iso == 0 || return false
  k_iso = condestJacobian(net_iso)
  (isfinite(k_iso) && k_iso >= 1.0) || return false
  st_iso = Sparlectra.rectangular_pf_status(net_iso)
  (st_iso !== nothing && hasproperty(st_iso, :jacobian_condest)) || return false
  k_iso_solver = st_iso.jacobian_condest()
  (isfinite(k_iso_solver) && k_iso_solver >= 1.0) || return false
  Sparlectra.printACPFlowResults(net_iso, 0.01, ite_iso, 1e-8, true, outdir; converged = true, solver = :rectangular)
  out_iso = read(joinpath(outdir, "result_$(net_iso.name).txt"), String)
  occursin("Jacobian cond. : kappa1(J)", out_iso) || return false
  return true
end

# Compact remove-functions coverage (replaces the former legacy/remove
# group): one small net, every public remover once,
# acceptance and rejection paths with network integrity afterwards. The
# rejection paths log @error by design; NullLogger keeps the run output
# clean exactly like the retired testremove.jl did.
function test_remove_functions_compact()
  net = Net(name = "remove_compact", baseMVA = 100.0)
  addBus!(net = net, busName = "B1", vn_kV = 110.0)
  addBus!(net = net, busName = "B2", vn_kV = 110.0)
  addBus!(net = net, busName = "B3", vn_kV = 110.0)
  addBus!(net = net, busName = "B5", vn_kV = 110.0)
  addBus!(net = net, busName = "B6", vn_kV = 20.0)
  add_test_line_pu_c!(net = net, fromBus = "B1", toBus = "B2", length_km = 10.0, r_pu = 0.01, x_pu = 0.1, c_nf_per_km = 10.0, tanδ = 0.0, ratedS = 100.0)
  add_test_line_pu_c!(net = net, fromBus = "B2", toBus = "B3", length_km = 15.0, r_pu = 0.01, x_pu = 0.1, c_nf_per_km = 10.0, tanδ = 0.0, ratedS = 100.0)
  add_test_line_pu_c!(net = net, fromBus = "B2", toBus = "B5", length_km = 15.0, r_pu = 0.01, x_pu = 0.1, c_nf_per_km = 10.0, tanδ = 0.0, ratedS = 100.0)
  add2WTrafo!(net = net, fromBus = "B1", toBus = "B6", sn_mva = 100.0, vk_percent = 10.0, vkr_percent = 0.5, pfe_kw = 20.0, i0_percent = 0.1)
  addShunt!(net = net, busName = "B2", pShunt = 0.0, qShunt = 10.0)
  addProsumer!(net = net, busName = "B3", type = "ENERGYCONSUMER", p = 50.0, q = 20.0)
  addProsumer!(net = net, busName = "B5", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.03, va_deg = 0.0, referencePri = "B5")
  ok, msg = validate!(net = net)
  @test ok
  n_branch = length(net.branchVec)
  Logging.with_logger(Logging.NullLogger()) do
    @test removeACLine!(net = net, fromBus = "B2", toBus = "B3")
    @test length(net.branchVec) == n_branch - 1
    @test !removeACLine!(net = net, fromBus = "B2", toBus = "B3")   # already gone
    @test removeTrafo!(net = net, fromBus = "B1", toBus = "B6")
    @test !removeTrafo!(net = net, fromBus = "B1", toBus = "B6")
    @test removeShunt!(net = net, busName = "B2")
    @test isempty(net.shuntVec)
    @test !removeShunt!(net = net, busName = "B2")
    n_pros = length(net.prosumpsVec)
    @test removeProsumer!(net = net, busName = "B3", type = "ENERGYCONSUMER")
    @test length(net.prosumpsVec) == n_pros - 1
    # removeBus! validates only: slack and connected buses are refused,
    # a fully detached bus reports removable
    @test !removeBus!(net = net, busName = "B5")   # slack
    @test !removeBus!(net = net, busName = "B2")   # still carries branches
    @test removeBus!(net = net, busName = "B3")    # detached now
  end
  return true
end

# #375: the PQ->PV release decides on the voltage side. The shipped
# Zeng/Chiang 14-bus case (loads +25 percent, tight Q bands) ended on a
# non-physical solution under active_set: three machines at Qmax with the
# voltage above the setpoint, zero releases, the Q-V check reporting all
# three. With the voltage rule the run reaches the physical solution the
# one-at-a-time mode finds (bus 3 clamped at Qmax, bus 1 at Qmin, bus 2
# back on PV) and the Q-V check is clean. The feeder
# with a machine clamped at Qmax and its voltage BELOW the setpoint must
# stay clamped.
function test_active_set_voltage_side_release()::Bool
  @testset "Active set re-enable on the voltage side (#375)" begin (function ()
    scf_dir = joinpath(dirname(@__DIR__), "data", "scf")
    zeng = joinpath(scf_dir, "case14_zeng_p306_activeSet_A.scf.json")
    cfg = Sparlectra.resolve_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, zeng).config
    @test cfg.powerflow.qlimits.enforcement_mode === :active_set
    r = run_sparlectra(casefile = zeng, config = cfg)
    @test r.final_converged
    st = Sparlectra.rectangular_pf_status(r.net)
    @test st.qlimit_reenable_events >= 1
    # bus 2 (Qmax, voltage above its setpoint) is released; bus 3 stays at
    # Qmax with the voltage below the setpoint, bus 1 at Qmin
    @test !haskey(r.net.qLimitEvents, 2)
    @test get(r.net.qLimitEvents, 3, nothing) === :max && get(r.net.qLimitEvents, 1, nothing) === :min
    @test isempty(Sparlectra.qvCharacteristicViolations(r.net))
    # Q-limit stability check (Zhao et al. 2008, (6)/(7)) after the run: no
    # violation under the active set; buses 3 and 6 sit at Qmax below U_set
    # under the one-retry lock (two hits)
    @test r.diagnostics.qlimit_stability_violations == 0
    @test [(row.bus, row.end_type, row.reached_limit, row.locked_by) for row in r.diagnostics.qlimit_stability_report] ==
          [(1, :PQ, :min, :none), (2, :PV, :none, :none), (3, :PQ, :max, :oscillation), (6, :PQ, :max, :oscillation)]
    # clamped at Qmax with the voltage below the setpoint: no release
    feeder = joinpath(scf_dir, "feeder3_hardpv_A.scf.json")
    fcfg = Sparlectra.resolve_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, feeder).config
    rf = run_sparlectra(casefile = feeder, config = fcfg)
    @test rf.final_converged
    b2 = rf.net.busDict["B2"]
    @test get(rf.net.qLimitEvents, b2, nothing) === :max
    @test Sparlectra.rectangular_pf_status(rf.net).qlimit_reenable_events == 0
    # the margin is a configuration key that reaches the net
    @test rf.net.reenable_v_hyst_pu == fcfg.powerflow.qlimits.reenable_v_hyst_pu

    # with the voltage-side release switched off (margin 1.0) the run ends
    # non-physical; fast/ref_api_smoke.jl checks that case. The Q-V rows of
    # that run are the event buses of the stability report on their limit on
    # the wrong side, so every one of them fails the stability check.
    cfg_a = Sparlectra.resolve_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, zeng, Dict{String,Any}("power_flow.qlimits.reenable_v_hyst_pu" => 1.0)).config
    ra = @test_logs (:warn, r"Q-limit stability check") match_mode = :any run_sparlectra(casefile = zeng, config = cfg_a)
    qv_a = Set(row.bus for row in Sparlectra.qvCharacteristicViolations(ra.net) if row.significant)
    @test !isempty(qv_a)
    @test issubset(qv_a, Set(row.bus for row in ra.diagnostics.qlimit_stability_report if !row.ok))
    # the classic one-at-a-time mode ends on a physical point as well
    cfg_c = Sparlectra.resolve_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, zeng, Dict{String,Any}("power_flow.qlimits.enforcement_mode" => "classic_one_at_a_time")).config
    # the stability check flags bus 6, PV 0.68 MVAr over its Qmax (the final
    # check accepts it within the hysteresis, see below)
    rc = @test_logs (:warn, r"Q-limit stability check: 1 voltage-regulated bus\(es\).*6 \(PV outside its band, U 1\.07000 pu, U_set 1\.07000 pu, Q 24\.678 MVAr outside -6\.000 \.\. 24\.000\)") match_mode = :any run_sparlectra(casefile = zeng, config = cfg_c)
    @test rc.final_converged
    @test rc.diagnostics.qlimit_stability_violations == 1
    @test [row.bus for row in rc.diagnostics.qlimit_stability_report if !row.ok] == [6]
    # classic_simultaneous ends on the abnormal solution of the paper: buses
    # 2, 3 and 6 at Qmax with the voltage above U_set violate (7); the
    # classic modes do not switch back.
    cfg_s = Sparlectra.resolve_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, zeng, Dict{String,Any}("power_flow.qlimits.enforcement_mode" => "classic_simultaneous")).config
    rs = @test_logs (:warn, r"Q-limit stability check: 3 voltage-regulated bus\(es\)") match_mode = :any run_sparlectra(casefile = zeng, config = cfg_s)
    @test rs.final_converged
    @test rs.diagnostics.qlimit_stability_violations == 3
    bad_s = [row for row in rs.diagnostics.qlimit_stability_report if !row.ok]
    @test [(row.bus, row.reached_limit, row.check_6, row.check_7, row.locked_by) for row in bad_s] == [(2, :max, true, false, :none), (3, :max, true, false, :none), (6, :max, true, false, :none)]
    @test all(row -> row.u_pu > row.uset_pu, bad_s)
    @test qlimit_stability_line(rs.diagnostics.qlimit_stability_status, rs.diagnostics.qlimit_stability_report) == "Q-limit stability violations: 3 (bus 2, 3, 6)"
    @test isempty(Sparlectra.qvCharacteristicViolations(rc.net))
    # The final Q-limit check judges both modes alike: the classic run keeps
    # bus 6 on PV over Qmax, inside the 1 MVAr hysteresis, and the check says
    # so instead of failing the run; the active-set run is
    # accepted with an ok or within_hysteresis verdict (up to 0.17.2 a
    # released machine inside the band failed it).
    st_c = Sparlectra.rectangular_pf_status(rc.net)
    @test st_c.final_q_check_status === :within_hysteresis
    row6 = only(filter(row -> row.busI == 6, st_c.final_q_check_rows))
    @test row6.side === :high && row6.class === :within_hysteresis
    @test st.final_q_check_status in (:ok, :within_hysteresis)
    @test st.final_q_check_status !== :remaining_pv_q_limit_violations
    # with Q-limit handling off a converged run has nothing to check: the
    # status says so instead of claiming there was no converged solution
    cfg_d = Sparlectra.resolve_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, zeng, Dict{String,Any}("power_flow.qlimits.enabled" => false)).config
    rd = run_sparlectra(casefile = zeng, config = cfg_d)
    @test rd.final_converged
    @test Sparlectra.rectangular_pf_status(rd.net).final_q_check_status === :qlimits_disabled
    @test (rd.diagnostics.qlimit_stability_status, rd.diagnostics.qlimit_stability_violations) == (:qlimits_disabled, 0)
    @test Sparlectra.final_q_check_line(:qlimits_disabled, NamedTuple[], 100.0) == "Final Q-limit check: not evaluated (Q-limit handling disabled)."
    @test Sparlectra.final_q_check_line(:not_evaluated, NamedTuple[], 100.0) == "Final Q-limit check: not evaluated (no converged solution)."
    # No bound on the voltage difference between the two modes on purpose:
    # switching variants need not end on the same point (measured:
    # 1.3e-3 pu apart at bus 6, the machine on its Qmax edge, clamped again by
    # the active set and kept on PV inside the hysteresis band by the classic
    # mode). What the test guards is that each mode ends on a PHYSICAL point:
    # the Q-V check above is clean for both.
    # on the shipped case118 the margin does not change the clamped set, and
    # nothing oscillates
    c118 = joinpath(dirname(@__DIR__), "data", "mpower", "sp_case118.m")
    clamped118 = Vector{Int}[]
    for margin in (1.0, 1e-4)
      cfg118 = Sparlectra.resolve_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, c118, Dict{String,Any}("power_flow.qlimits.reenable_v_hyst_pu" => margin)).config
      r118 = run_sparlectra(casefile = c118, config = cfg118)
      @test r118.final_converged
      push!(clamped118, sort!(collect(keys(r118.net.qLimitEvents))))
      @test Sparlectra.rectangular_pf_status(r118.net).oscillating_buses == 0
      @test Sparlectra.rectangular_pf_status(r118.net).final_q_check_status === :ok
    end
    @test clamped118[1] == clamped118[2]
    # the feeder machine sits ON its clamp, so the final check has nothing to say
    @test Sparlectra.rectangular_pf_status(rf.net).final_q_check_status === :ok
  end)() end
  return true
end

# The final Q-limit check every enforcement mode ends with (task
# final-q-check): judged by the SIZE of the overshoot with two
# thresholds, the switching hysteresis and final_q_accept_pu. Synthetic
# injections on the 3-bus net, so every class is hit exactly, plus the
# configuration default and validation and the warning of the bounded class.
function test_final_q_limit_classification()::Bool
  @testset "Final Q-limit check: two thresholds" begin (function ()
    net = createTest3BusNet(cooldown = 1, hyst_pu = 0.01, qlim_min = -15.0, qlim_max = 15.0)
    pv = geNetBusIdx(net = net, busName = "STATION1")
    nb = length(net.nodeVec)
    bus_types = [Sparlectra.getNodeType(node) == Sparlectra.Slack ? :Slack : Sparlectra.getNodeType(node) == Sparlectra.PV ? :PV : :PQ for node in net.nodeVec]
    @test bus_types[pv] == :PV
    qmin_pu = fill(-Inf, nb)
    qmax_pu = fill(Inf, nb)
    qmin_pu[pv] = -0.15
    qmax_pu[pv] = 0.15
    qload = Sparlectra.build_qload_pu(net)
    # the injection that puts the machine at Qmax + dev (or Qmin - dev)
    injections(dev; side = :high) = [ComplexF64(0.0, i == pv ? ((side === :high ? qmax_pu[pv] + dev : qmin_pu[pv] - dev) - qload[pv]) : 0.0) for i in 1:nb]
    classify(dev; hyst = 0.01, accept = 0.02, side = :high) = classify_final_q_limits(net, injections(dev; side), bus_types, qmin_pu, qmax_pu; q_hyst_pu = hyst, final_q_accept_pu = accept, tol = 1e-8)
    # at the limit within the tolerance: no row, ok
    @test classify(0.0).status === :ok
    @test isempty(classify(0.0).rows)
    # inside the hysteresis: a row, accepted, class within_hysteresis
    within = classify(0.005)
    @test within.status === :within_hysteresis
    @test length(within.rows) == 1 && within.rows[1].class === :within_hysteresis && within.rows[1].side === :high
    @test isapprox(within.rows[1].dev_pu, 0.005; atol = 1e-12)
    @test within.max_dev_pu == within.rows[1].dev_pu
    @test occursin("bus $(within.rows[1].busI) over Qmax by 0.0050 pu / 0.50 MVAr (within_hysteresis)", final_q_check_line(within.status, within.rows, net.baseMVA))
    # beyond the hysteresis, inside the bound: bounded
    bounded = classify(0.015)
    @test bounded.status === :bounded_q_limit_violation
    @test bounded.bounded == 1 && bounded.violations == 0
    @test bounded.buses == string(bounded.rows[1].busI, ":bounded")
    # beyond the bound: a remaining violation, the run is not accepted
    @test classify(0.03).status === :remaining_pv_q_limit_violations
    @test classify(0.03).violations == 1
    # the mirror image under Qmin
    low = classify(0.015; side = :low)
    @test low.rows[1].side === :low && low.rows[1].class === :bounded
    # regression: both thresholds at zero is the strict check of 0.17.2
    @test classify(0.005; hyst = 0.0, accept = 0.0).status === :remaining_pv_q_limit_violations
    # a PQ bus is not judged, whatever its injection
    bus_types_pq = copy(bus_types)
    bus_types_pq[pv] = :PQ
    @test classify_final_q_limits(net, injections(0.03), bus_types_pq, qmin_pu, qmax_pu; q_hyst_pu = 0.01, final_q_accept_pu = 0.02, tol = 1e-8).status === :ok

    # the acceptance decision: within and bounded are accepted, bounded with
    # a warning naming bus and overshoot, a violation rejects the run
    finalize(dev) = Sparlectra._finalize_rectangular_qlimit_summary(net, ones(ComplexF64, nb), injections(dev), bus_types, qmin_pu, qmax_pu, true, :none;
      verbose = 0, qlimit_trace_enabled = false, q_hyst_pu = 0.01, final_q_accept_pu = 0.02, tol = 1e-8, pv_table_rows = 30,
      qlimit_guard_accept_bounded_violations = false, qlimit_guard_max_remaining_violations = 0)
    @test finalize(0.005).converged
    @test finalize(0.005).final_q_check.status === :within_hysteresis
    accepted = @test_logs (:warn, r"bounded violation accepted") finalize(0.015)
    @test accepted.converged
    @test accepted.rejection_reason === :bounded_q_limit_violation
    rejected = finalize(0.03)
    @test !rejected.converged
    @test rejected.rejection_reason === :remaining_pv_q_limit_violations

    # the key: default twice the hysteresis, never below it, stamped on the net
    cfg = Sparlectra.QLimitConfig(Dict{String,Any}("hysteresis_pu" => 0.01))
    @test cfg.final_q_accept_pu == 0.02
    @test Sparlectra.QLimitConfig(Dict{String,Any}("hysteresis_pu" => 0.0)).final_q_accept_pu == 0.0
    @test Sparlectra.QLimitConfig(Dict{String,Any}("hysteresis_pu" => 0.01, "final_q_accept_pu" => 0.05)).final_q_accept_pu == 0.05
    @test_throws ArgumentError Sparlectra.QLimitConfig(Dict{String,Any}("hysteresis_pu" => 0.01, "final_q_accept_pu" => 0.005))
    # auto (the template value) follows the hysteresis, so a user file that
    # raises hysteresis_pu alone is not rejected by an inherited literal
    @test Sparlectra.QLimitConfig(Dict{String,Any}("hysteresis_pu" => 0.03, "final_q_accept_pu" => "auto")).final_q_accept_pu == 0.06
    @test "power_flow.qlimits.final_q_accept_pu" in Sparlectra.GUI_EDITABLE_CONFIG_KEYS
    @test Net(name = "n", baseMVA = 100.0, q_hyst_pu = 0.03).final_q_accept_pu == 0.06

    # the same status for every mode at the same operating point: sp_case14
    # under the three enforcement modes
    case = joinpath(dirname(@__DIR__), "data", "scf", "sp_case14.scf.json")
    statuses = Symbol[]
    for mode in ("active_set", "classic_one_at_a_time", "classic_simultaneous")
      cfg14 = Sparlectra.resolve_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, case, Dict{String,Any}("power_flow.qlimits.enforcement_mode" => mode)).config
      r14 = run_sparlectra(casefile = case, config = cfg14)
      @test r14.final_converged
      push!(statuses, Sparlectra.rectangular_pf_status(r14.net).final_q_check_status)
    end
    @test length(unique(statuses)) == 1
    @test statuses[1] in (:ok, :within_hysteresis)
  end)() end
  return true
end

function run_grid_fast_tests()
  @testset "Remove functions" begin (function ()
    @test test_remove_functions_compact() == true
  end)() end
  @testset "Transformer and network validation" begin (function ()
    @test test_transformer_nameplate_conversion() == true
    @test test_shunt_update_matches_constructor() == true
    @test testISOBusses() == true
  end)() end

  @testset "MATPOWER import/export" begin (function ()
    @test test_matpower_import_uses_bus_type_for_regulation() == true
    @test test_matpower_flatstart_uses_generator_voltage_setpoints() == true
    @test test_prosumer_aggregation_preserves_bus_types_and_injections() == true
    @test test_matpower_build_ybus_returns_sparse_expected_values() == true
    @test test_matpower_vmva_selfcheck_ignores_slack_pq_spec() == true
    @test test_matpower_compare_vmva_wraps_angle_differences() == true
    @test test_matpower_pv_voltage_source_and_compare_modes() == true
    @test test_bus_shunt_model_modes() == true
    @test test_matpower_read_case_m_postprocessing() == true
    @test test_matpower_matrix_block_scanner_whitespace() == true
    @test test_matpower_file_import_honors_explicit_overrides() == true
    @test test_run_sparlectra_forwards_wrong_branch_config() == true
    @test test_wrong_branch_voltage_band_every_level() == true
    @test test_wrong_branch_reference_branch_angle() == true
    @test test_run_sparlectra_normalizes_projected_matpower_starts() == true
  end)() end

  @testset "Power flow scenarios" begin (function ()
    @test test_rectangular_autodamp_backtracks_oversized_step() == true
    @test test_rectangular_polar_update_rotates_without_inflation() == true
    @test test_rectangular_polar_update_reaches_the_same_solution() == true
    @test test_power_mode_same_solution_and_topology_change() == true
    @test test_jacobian_reuse_same_solution() == true
    @test test_rectangular_nonfinite_mismatch_diagnostics_use_finite_history() == true
    @test test_rectangular_start_projection_dc_candidate_matches_rundcpf() == true
    @test test_rectangular_start_projection_keeps_raw_without_finite_improvement() == true
    @test test_rectangular_start_projection_ratio_profile() == true
    @test test_run_level_flat_start_is_flat() == true
    @test test_flat_start_carries_phase_shift() == true
    @test test_q_limit_adjust_vset_modes() == true
    @test test_q_limit_matpower_mode_base_failure_does_not_switch() == true
    @test test_q_limit_reenable_voltage_rule() == true
    @test test_q_limit_start_iter_delays_switching() == true
    @test test_q_limit_hysteresis_delays_small_pv_to_pq_overshoot() == true
    @test test_solve_sparse_system_small_and_large_backends() == true
    @test test_rectangular_jacobian_and_dc_start_matrices_are_sparse() == true
    @test test_rectangular_singular_step_returns_nonconvergence() == true
    @test test_condition_number_estimator() == true
    @test test_bus_type_resolution_from_prosumers() == true
    @test test_regulated_generator_bus_targets_include_unregulated_generators() == true
    @test test_final_q_limit_classification() == true
    @test test_active_set_voltage_side_release() == true
  end)() end

  @testset "Link behaviour and reporting" begin (function ()
    @test test_link_kcl_simple() == true
    @test test_link_rejects_slack_connection() == true
    @test test_link_bus_merge_pf() == true
    @test test_link_merge_remaps_non_representative_slack() == true
    @test test_link_merge_returns_bus_powers() == true
    @test test_island_solve_with_link_merge() == true
    @test test_link_closed_keeps_shunt_reporting_on_original_bus() == true
    @test test_link_kcl_ring_allocation() == true
    @test test_link_ring_pf_stability() == true
    @test test_summary_result_output_closes_file_handle() == true
  end)() end
end

function run_grid_extended_tests()
  @testset "large/stress MATPOWER helpers" begin (function ()
    @test test_matpower_matrix_block_scanner_large_body() == true
  end)() end
end

run_grid_tests() = run_grid_fast_tests()
