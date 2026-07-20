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

# Date: 2026-07-20
# file: examples/exp_3wt_phase_taps.jl
# purpose: builds a 3WT in three tap configurations (OLTC only, PST/Schraegregler only, both combined) and solves each with runpf!, demonstrating create3WTWindings!'s phase_tap_side/phase_taps keywords

using Sparlectra

include(joinpath(@__DIR__, "internal", "example_header.jl"))

# Physical parameters of the three AUX-bus legs (Ohm/Siemens/MVA), shared by
# build_case_net and print_case_topology so both stay in sync.
const LEG_R_OHM = (0.20, 0.30, 0.40)
const LEG_X_OHM = (4.00, 6.00, 10.00)
const LEG_B_SIEMENS = (0.0, 0.0, 0.0)
const LEG_S_MVA = (1000.0, 500.0, 200.0)
const LEG_NAMES = ("HV (AUX→B2)", "MV (AUX→B3)", "LV (AUX→B4)")

# --------------------------------------------------------------------------
# add3WTPiModelTrafo! (the only 3WT path actually wired into a Net) always
# builds all three AUX-bus legs with ratio=1.0 / shift_deg=0.0 — it has no
# hook for a tap position, and create3WTWindings!'s PowerTransformerTaps /
# PhaseTapChangerModel objects are not wired into any Net at all (see
# docs/dev/3wt_phase_tap_controller_addressing.md, §3).
#
# This helper mirrors add3WTPiModelTrafo!'s own internal logic 1:1 (same
# AUX-bus creation, same three add2WTPIModelTrafo! calls), but adds the one
# thing missing: it lets the caller override ratio/shift_deg on exactly one
# leg. That override is computed *outside*, with the real tap-changer
# formulas (calcRatioTapCorrection for an OLTC, calcPhaseTapAngleRatio for a
# PST) — no new solver logic, only the missing link between an existing tap
# model and an existing branch parameter.
function add_3wt_with_leg_override!(;
  net::Net,
  HBBus::String,
  MBBus::String,
  LVBus::String,
  r::NTuple{3,Float64},
  x::NTuple{3,Float64},
  b::NTuple{3,Float64},
  ratedU_kV::NTuple{3,Float64},
  ratedS_MVA::NTuple{3,Float64},
  status::Int = 1,
  override_leg::Int = 0,          # 1=HV, 2=MV, 3=LV; 0 = no override
  override_ratio::Float64 = 1.0,
  override_shift_deg::Float64 = 0.0,
)
  @assert override_leg in (0, 1, 2, 3) "override_leg must be 0, 1, 2 or 3"

  U_aux_kV = maximum(ratedU_kV)
  aux_bus = "Aux3WT_" * HBBus * "_" * MBBus * "_" * LVBus
  if !hasBusInNet(net = net, busName = aux_bus)
    addBus!(net = net, busName = aux_bus, vn_kV = U_aux_kV, isAux = true)
  end

  ratios = [1.0, 1.0, 1.0]
  shifts = [0.0, 0.0, 0.0]
  if override_leg != 0
    ratios[override_leg] = override_ratio
    shifts[override_leg] = override_shift_deg
  end

  buses = (HBBus, MBBus, LVBus)
  for i = 1:3
    add2WTPIModelTrafo!(net = net, fromBus = aux_bus, toBus = buses[i], side = 1, r = r[i], x = x[i], b = b[i], status = status, ratedU = U_aux_kV, ratedS = ratedS_MVA[i], ratio = ratios[i], shift_deg = shifts[i])
  end

  return aux_bus
end

# Same small test grid for all three cases — only the leg override changes,
# and (for Case 3) the winding-voltage layout.
#
# Normal layout (Case 1 OLTC, Case 2 PST asymmetrical): HV/MV/LV genuinely
# different (380/110/20 kV) — a classic autotransformer-style 3WT.
#
# PST-symmetrical layout (Case 3): HV=MV (here 380/380 kV) — the phase
# shifter sits in-line between two busbars of EQUAL voltage, only the third
# (excitation) winding sits at a different voltage level. This is not a
# special case, it's the real-world design of a phase-shifting transformer
# / quadrature booster: it never changes voltage magnitude by definition
# (effective_ratio for :symmetrical is always 1.0, see
# calcPhaseTapAngleRatio), so a voltage jump between the two series windings
# would not make sense anyway.
#
#   HV/MV kV                   AUX*        LV kV
#   Slack o--AC-Line--o-------------o----------o--AC-Line--o  Load: 80 MW, 30 MVar
#         B1          B2            | LV kV
#                                    o B4  Shunt: 5 MVar
function build_case_net(; hv_kV::Float64 = 380.0, mv_kV::Float64 = 110.0, lv_kV::Float64 = 20.0, override_leg::Int = 0, override_ratio::Float64 = 1.0, override_shift_deg::Float64 = 0.0)
  net = Net(name = "3wt_tap_case_demo", baseMVA = 1000.0)

  addBus!(net = net, busName = "B1", vn_kV = hv_kV)
  addBus!(net = net, busName = "B2", vn_kV = hv_kV)  # HV side of the 3WT
  addBus!(net = net, busName = "B3", vn_kV = mv_kV)  # MV side
  addBus!(net = net, busName = "B4", vn_kV = lv_kV)  # LV side (excitation / shunt bus)
  addBus!(net = net, busName = "B5", vn_kV = mv_kV)  # load bus (at MV level)

  addProsumer!(net = net, busName = "B1", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "B1")
  addACLine!(net = net, fromBus = "B1", toBus = "B2", length = 1.0, r = 0.01, x = 0.10)

  add_3wt_with_leg_override!(net = net, HBBus = "B2", MBBus = "B3", LVBus = "B4", r = LEG_R_OHM, x = LEG_X_OHM, b = LEG_B_SIEMENS, ratedU_kV = (hv_kV, mv_kV, lv_kV), ratedS_MVA = LEG_S_MVA, override_leg = override_leg, override_ratio = override_ratio, override_shift_deg = override_shift_deg)

  addACLine!(net = net, fromBus = "B3", toBus = "B5", length = 1.0, r = 0.01, x = 0.10)
  addProsumer!(net = net, busName = "B5", type = "ENERGYCONSUMER", p = 80.0, q = 30.0)
  addShunt!(net = net, busName = "B4", pShunt = 0.0, qShunt = 5.0)

  return net
end

# Prints the test-network topology for one case — buses, voltage levels and
# the three AUX legs with their r/x and (possibly overridden) ratio/shift_deg
# — so the network being solved is visible, not just the result table.
function print_case_topology(; hv_kV::Float64, mv_kV::Float64, lv_kV::Float64, override_leg::Int, override_ratio::Float64, override_shift_deg::Float64)
  aux_kV = maximum((hv_kV, mv_kV, lv_kV))
  ratios = [1.0, 1.0, 1.0]
  shifts = [0.0, 0.0, 0.0]
  if override_leg != 0
    ratios[override_leg] = override_ratio
    shifts[override_leg] = override_shift_deg
  end

  println("  Test network:")
  println("    B1 (", hv_kV, " kV, Slack)  --AC line--  B2 (", hv_kV, " kV, 3WT HV terminal)")
  println("    AUX (", aux_kV, " kV, star point)  --  3 PI-model legs:")
  for i = 1:3
    tag = (i == override_leg) ? "   <- tap applied here" : ""
    println("      ", LEG_NAMES[i], "   r=", LEG_R_OHM[i], " Ω  x=", LEG_X_OHM[i], " Ω  ratio=", round(ratios[i], digits = 5), "  shift=", round(shifts[i], digits = 3), "°", tag)
  end
  println("    B3 (", mv_kV, " kV, 3WT MV terminal)  --AC line--  B5 (", mv_kV, " kV, Load 80 MW / 30 MVar)")
  println("    B4 (", lv_kV, " kV, 3WT LV terminal)  --  Shunt 5 MVar")
  println()
end

# Builds the case net, runs runpf!, and prints the full result table.
function run_case(label::String; hv_kV::Float64 = 380.0, mv_kV::Float64 = 110.0, lv_kV::Float64 = 20.0, override_leg::Int, override_ratio::Float64, override_shift_deg::Float64, note::Union{Nothing,String} = nothing)
  leg_name = override_leg == 1 ? "HV leg (AUX→B2)" : override_leg == 2 ? "MV leg (AUX→B3)" : "LV leg (AUX→B4)"
  println(label)
  println("  ", leg_name, "  ratio=", round(override_ratio, digits = 5), "  shift_deg=", round(override_shift_deg, digits = 3), "°")
  println()

  if !isnothing(note)
    println("  >>> ", note)
    println()
  end

  print_case_topology(hv_kV = hv_kV, mv_kV = mv_kV, lv_kV = lv_kV, override_leg = override_leg, override_ratio = override_ratio, override_shift_deg = override_shift_deg)

  net = build_case_net(hv_kV = hv_kV, mv_kV = mv_kV, lv_kV = lv_kV, override_leg = override_leg, override_ratio = override_ratio, override_shift_deg = override_shift_deg)

  result, msg = validate!(net = net)
  if !result
    @warn msg
    return nothing
  end

  etime = @elapsed begin
    ite, erg = runpf!(net, 50, 1e-9, 0)
  end
  if erg != 0
    @warn "Power flow did not converge"
    return nothing
  end

  calcNetLosses!(net)
  distributeBusResults!(net)
  printACPFlowResults(net, etime, ite, 1e-9)
  println()

  return net
end

function print_task_description()
  println("="^78)
  println("Task: the same small 3WT test grid is solved three times (one runpf! each).")
  println("Topology, load, shunt and slack stay identical — only the tap setting on one")
  println("leg of the 3WT star-equivalent changes between cases.")
  println()
  println("  Case 1 – pure OLTC on the 3WT:")
  println("           ratio tap on the HV leg (AUX→B2), computed with")
  println("           calcRatioTapCorrection.")
  println("  Case 2 – PST asymmetrical (Schrägregler):")
  println("           ratio+phase tap on the MV leg (AUX→B3), computed with")
  println("           calcPhaseTapAngleRatio, ψ=60°.")
  println("  Case 3 – PST symmetrical:")
  println("           phase-only tap on the MV leg (AUX→B3), computed with")
  println("           calcPhaseTapAngleRatio, ratio stays 1.0. Here HV=MV=380kV")
  println("           (an in-line phase shifter between two 380kV busbars), only the")
  println("           LV winding (excitation) stays at 20kV — Case 1/2 still use")
  println("           HV/MV/LV = 380/110/20 kV.")
  println()
  println("Note: create3WTWindings!'s PowerTransformerTaps/PhaseTapChangerModel objects")
  println("are not wired into add3WTPiModelTrafo! (the only Net-capable 3WT path) — see")
  println("docs/dev/3wt_phase_tap_controller_addressing.md. The ratio/shift_deg below are")
  println("computed with the real tap-changer formulas and applied directly to the")
  println("affected AUX leg, as a manual stand-in for that missing wiring.")
  println("="^78)
  println()
end

"""
    main()

Runs the three 3WT tap cases (pure OLTC, PST asymmetrical/Schrägregler, PST
symmetrical) described by [`print_task_description`](@ref) and returns the
solved `Net` for each (`nothing` for a case whose power flow did not
converge).
"""
function main()
  print_example_banner("examples/exp_3wt_phase_taps.jl", "builds a 3WT in three tap configurations (OLTC only, PST/Schraegregler only, both combined) and solves each with runpf!, demonstrating create3WTWindings!'s phase_tap_side/phase_taps keywords")
  print_task_description()

  # Case 1: pure OLTC on the HV leg.
  oltc_tap = PowerTransformerTaps(Vn_kV = 380.0, step = 5, lowStep = -9, highStep = 9, neutralStep = 0, voltageIncrement_kV = 3.8)
  ratio_oltc = calcRatioTapCorrection(oltc_tap)
  net_oltc = run_case("Case 1: 3WT with pure OLTC (step=$(oltc_tap.step))"; override_leg = 1, override_ratio = ratio_oltc, override_shift_deg = 0.0)

  # Case 2: PST asymmetrical on the MV leg.
  pst_asym = PhaseTapChangerModel(kind = :asymmetrical, step = 5, lowStep = -8, highStep = 8, neutralStep = 0, voltage_step_increment = 0.0125, winding_connection_angle_deg = 60.0)
  tap_asym = calcPhaseTapAngleRatio(pst_asym)
  net_pst_asym = run_case("Case 2: 3WT with PST asymmetrical (Schrägregler, step=$(pst_asym.step))"; override_leg = 2, override_ratio = tap_asym.effective_ratio, override_shift_deg = tap_asym.effective_shift_deg)

  # Case 3: PST symmetrical on the MV leg — HV=MV=380kV (in-line phase shifter
  # between two 380kV busbars), only the LV winding (excitation) stays at 20kV.
  pst_sym = PhaseTapChangerModel(kind = :symmetrical, step = 5, lowStep = -8, highStep = 8, neutralStep = 0, voltage_step_increment = 0.0125)
  tap_sym = calcPhaseTapAngleRatio(pst_sym)
  net_pst_sym = run_case(
    "Case 3: 3WT with PST symmetrical (step=$(pst_sym.step))";
    hv_kV = 380.0,
    mv_kV = 380.0,
    lv_kV = 20.0,
    override_leg = 2,
    override_ratio = tap_sym.effective_ratio,
    override_shift_deg = tap_sym.effective_shift_deg,
    note = "Special case: this 3-winding model explicitly represents the excitation branch as a 3rd winding on its own voltage level (20kV). A textbook symmetric PST is normally built as a plain 2-winding transformer on ONE voltage level (both terminals equal), affecting only the active power flow P — its voltage-magnitude ratio never changes, so there is no reactive/voltage effect to speak of. The 3-winding view here is the more detailed, physical picture, not the standard simplified one.",
  )

  return (oltc = net_oltc, pst_asymmetrical = net_pst_asym, pst_symmetrical = net_pst_sym)
end

run_example(main)