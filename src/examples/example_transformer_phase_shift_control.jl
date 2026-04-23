using Sparlectra

"""
    build_phase_probe_branch(; phase_shift_deg, tap_ratio=1.0) -> Branch

Builds a minimal 2-bus network with one PI-model transformer branch and returns
that branch. The setup is intentionally tiny so phase-shift sign behavior can be
inspected in isolation.
"""
function build_phase_probe_branch(; phase_shift_deg::Float64, tap_ratio::Float64 = 1.0)::Branch
  net = Net(name = "phase_shift_control_demo", baseMVA = 100.0)
  addBus!(net = net, busName = "A", vn_kV = 110.0)
  addBus!(net = net, busName = "B", vn_kV = 110.0)
  addPIModelTrafo!(
    net = net,
    fromBus = "A",
    toBus = "B",
    r_pu = 0.01,
    x_pu = 0.10,
    b_pu = 0.02,
    status = 1,
    ratio = tap_ratio,
    shift_deg = phase_shift_deg,
  )
  return only(net.branchVec)
end

"""
    p_from_to(branch, V_from, V_to)

Computes active power flow P_ab from bus A -> B using branch admittance terms:
`I_from = Y_ff*V_from + Y_ft*V_to`, `S_from = V_from*conj(I_from)`, `P_ab = real(S_from)`.
"""
function p_from_to(branch::Branch, V_from::ComplexF64, V_to::ComplexF64)::Float64
  Y_ff, Y_ft, _, _ = calcAdmittance(branch, 110.0, 100.0)
  I_from = Y_ff * V_from + Y_ft * V_to
  return real(V_from * conj(I_from))
end

"""
    recommend_phase_step(current_p, target_p, delta_p_probe; step_deg=0.5)

Controller rule without fixed sign assumption:
1) Probe with `Δφ > 0` and measure `delta_p_probe = P(φ+Δφ) - P(φ)`.
2) If `current_p < target_p`, move φ in the direction that increases P.
3) Otherwise move φ in the opposite direction.
"""
function recommend_phase_step(current_p::Float64, target_p::Float64, delta_p_probe::Float64; step_deg::Float64 = 0.5)::Float64
  increase_direction = sign(delta_p_probe)
  return (current_p < target_p) ? increase_direction * step_deg : -increase_direction * step_deg
end

function main()
  println("=== Transformer phase-shift control demo (Sparlectra convention) ===")
  println("t = tap_ratio * cis(deg2rad(phase_shift_deg))")

  # fixed terminal voltages for clear sensitivity probe
  V_a = 1.0 + 0.0im
  V_b = 1.0 + 0.0im

  br0 = build_phase_probe_branch(phase_shift_deg = 0.0)
  br5 = build_phase_probe_branch(phase_shift_deg = 5.0)

  p0 = p_from_to(br0, V_a, V_b)
  p5 = p_from_to(br5, V_a, V_b)
  delta_p = p5 - p0

  println("P_ab(phi=0°) = $(round(p0, digits = 8)) pu")
  println("P_ab(phi=5°) = $(round(p5, digits = 8)) pu")
  println("ΔP_ab = P(5°)-P(0°) = $(round(delta_p, digits = 8)) pu")
  println(delta_p > 0 ? "Result: +5° increases P_ab." : "Result: +5° decreases P_ab.")

  target = p0 + 0.05
  step = recommend_phase_step(p0, target, delta_p; step_deg = 0.5)
  println("Target example: target = P_ab(phi=0°)+0.05 pu")
  println("Recommended phase update step = $(step)°")
end

Base.invokelatest(main)
