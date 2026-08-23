# Copyright 2023-2026 Udo Schmitz                                             #src
#                                                                             #src
# Licensed under the Apache License, Version 2.0 (the "License");             #src
# you may not use this file except in compliance with the License.            #src
# You may obtain a copy of the License at                                     #src
#                                                                             #src
#     http://www.apache.org/licenses/LICENSE-2.0                              #src
#                                                                             #src
# Unless required by applicable law or agreed to in writing, software         #src
# distributed under the License is distributed on an "AS IS" BASIS,           #src
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.    #src
# See the License for the specific language governing permissions and         #src
# limitations under the License.                                              #src
#                                                                             #src
# file: docs/lit/workshop_transformers.jl                                     #src
# purpose: Literate.jl source of the transformer-taps workshop notebook and   #src
#          its Documenter page: ratio taps (OLTC) vs phase taps (PST /        #src
#          Schraegregler), the tap-changer device math                        #src
#          (calcRatioTapCorrection, calcPhaseTapAngleRatio), the             #src
#          phase-tap outer control loop, and the 3WT star equivalent with a   #src
#          tap on one leg. Consolidates exp_3wt_phase_taps.jl and             #src
#          exp_pst_reactance_coupling.jl into one notebook. Regenerate the    #src
#          committed outputs with                                             #src
#          `julia --project=docs docs/generate_notebooks.jl`.                 #src

# # Transformer taps: moving voltage, moving power
#
# > **Level: Advanced**, companion of basic-tour chapter 4.
#
# [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_transformers.ipynb)
#
# > **Note:** This workshop was created with AI assistance and is reviewed
# > and curated by the maintainer; it is not a fully machine-generated text.
#
# A transformer enters the power flow with exactly two adjustable numbers:
# the winding ratio $\tau$ and the phase shift $\varphi$. The ratio tap
# (OLTC) moves $\tau$ and thereby VOLTAGE; the phase tap (PST, in German a
# Schrägregler) moves $\varphi$ and thereby ACTIVE POWER. In this notebook
# you compute real tap positions with the device formulas of
# [Sparlectra.jl](https://github.com/Welthulk/Sparlectra.jl), watch each
# tap type move its quantity, close the loop with the outer phase-tap
# controller, and finish with a three-winding transformer carrying a tap
# on one leg.
#
# > **Note:** On Google Colab the install cell takes a few minutes on a
# > fresh session (package download and precompilation). Colab's Julia
# > version may change over time; this notebook targets Julia ≥ 1.12.

#nb # ## Setup (Colab)
#nb # This cell installs Sparlectra from GitHub (branch `main`) into a fresh
#nb # temporary environment. The isolation matters on Colab: the shared
#nb # default environment ships many preinstalled packages, and installing
#nb # anything there triggers precompilation of that whole stack. Run this
#nb # cell first, once per session; it takes a few minutes.
#nb using Pkg
#nb Pkg.activate(temp = true)
#nb Pkg.add(url = "https://github.com/Welthulk/Sparlectra.jl", rev = "main")
#nb ## To test another branch, set rev to its name, e.g. rev = "dev/r0.9.8".
#nb ## For the latest registered release use: Pkg.add("Sparlectra")
#nb ## Switching versions in a running session? A "[loaded: ...]" note means
#nb ## the old version is still active — restart the runtime, then rerun
#nb ## this cell.

# ## Warm-up and shared helpers
#
# Julia compiles each function on first use. This cell loads the package,
# defines the helpers all chapters share, and warms the two paths the
# notebook exercises (plain power flow and the outer control loop), so
# the chapters below run at full speed.

using Sparlectra

## solve helper: converged power flow or an error
function solve!(net)
  ite, erg = runpf!(net, 50, 1e-9, 0)
  erg == 0 || error("Power flow did not converge (status = $erg)")
  calcNetLosses!(net)
  return ite
end

bus_vm(net, bus) = round(get_bus_vm_pu(net, bus); digits = 4)

## warm both paths on a throwaway 2-bus net with a transformer
wnet = Net(name = "warmup", baseMVA = 100.0)
addBus!(net = wnet, busName = "A", vn_kV = 110.0)
addBus!(net = wnet, busName = "B", vn_kV = 20.0)
addProsumer!(net = wnet, busName = "A", type = "EXTERNALNETWORKINJECTION", referencePri = "A", vm_pu = 1.0, va_deg = 0.0)
addProsumer!(net = wnet, busName = "B", type = "ENERGYCONSUMER", p = 10.0, q = 3.0)
addPIModelTrafo!(net = wnet, fromBus = "A", toBus = "B", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, ratio = 1.0, shift_deg = 0.0, status = 1)
t_pf = @elapsed solve!(wnet)
addPowerTransformerControl!(wnet; trafo = "1", mode = :branch_active_power, target_branch = ("A", "B"), p_target_mw = 10.0, control_ratio = false, control_phase = true, deadband_p_mw = 5.0)
wnet.branchVec[1].has_phase_tap = true
wnet.branchVec[1].phase_min_deg = -10.0
wnet.branchVec[1].phase_max_deg = 10.0
wnet.branchVec[1].phase_step_deg = 0.5
t_ctrl = @elapsed run_sparlectra(net = wnet)
println("warm: power flow ", round(t_pf; digits = 2), " s, control loop ", round(t_ctrl; digits = 2), " s (first calls compile)")

# ## Two taps, two jobs
#
# Every transformer branch is stamped into the Y-bus as (see the
# [Branch Model](https://welthulk.github.io/Sparlectra.jl/branchmodel/)
# page):
#
# ```math
# Y_{br} = \begin{bmatrix}
#   \frac{1}{\tau^2}\left(y_{ser} + \frac{y_{shunt}}{2}\right) &
#   -y_{ser}\,\frac{1}{\tau e^{-j\varphi}} \\
#   -y_{ser}\,\frac{1}{\tau e^{j\varphi}} &
#   y_{ser} + \frac{y_{shunt}}{2}
# \end{bmatrix}
# ```
#
# $\tau$ scales voltage magnitude across the transformer, so the RATIO tap
# is the voltage lever. $\varphi$ shifts the voltage ANGLE, and since
# $P \sim \sin(\delta_1 - \delta_2)$, the PHASE tap is the active-power
# lever. Sparlectra ships the device math for both:
# `PowerTransformerTaps` with `calcRatioTapCorrection` for the OLTC, and
# `PhaseTapChangerModel` with `calcPhaseTapAngleRatio` for the PST.
#
# ## Chapter 1: the ratio tap (OLTC) moves voltage
#
# **Example 1: the ratio tap moves voltage.** A 380/110 kV transformer
# feeds a 60 MW load. The OLTC has 19 positions
# (-9 to +9) of 3.8 kV each, exactly 1 percent of 380 kV per step. We
# compute $\tau$ for a few positions with the device formula and watch the
# load-side voltage follow:
#
# ```text
#   S (slack, 380 kV) --- line --- A ==OLTC== B (110 kV, load 60 MW)
# ```

function feeder(ratio)
  net = Net(name = "oltc_feeder", baseMVA = 100.0)
  addBus!(net = net, busName = "S", vn_kV = 380.0)
  addBus!(net = net, busName = "A", vn_kV = 380.0)
  addBus!(net = net, busName = "B", vn_kV = 110.0)
  addProsumer!(net = net, busName = "S", type = "EXTERNALNETWORKINJECTION", referencePri = "S", vm_pu = 1.0, va_deg = 0.0)
  addProsumer!(net = net, busName = "B", type = "ENERGYCONSUMER", p = 60.0, q = 20.0)
  addPIModelACLine!(net = net, fromBus = "S", toBus = "A", r_pu = 0.01, x_pu = 0.06, b_pu = 0.0, status = 1)
  addPIModelTrafo!(net = net, fromBus = "A", toBus = "B", r_pu = 0.01, x_pu = 0.10, b_pu = 0.0, ratio = ratio, shift_deg = 0.0, status = 1)
  validate!(net = net)
  solve!(net)
  return net
end

for step in (-5, 0, 5)
  taps = PowerTransformerTaps(Vn_kV = 380.0, step = step, lowStep = -9, highStep = 9, neutralStep = 0, voltageIncrement_kV = 3.8)
  tau = calcRatioTapCorrection(taps)
  net = feeder(tau)
  println("step ", lpad(step, 3), ": tau = ", round(tau; digits = 4), "  ->  Vm(B) = ", bus_vm(net, "B"), " pu")
end

# Reading aid (Example 1): five steps move tau by 5 percent and the load voltage by
# roughly the same amount, in the opposite direction: raising tau means
# more turns on the primary, so the secondary voltage drops. This manual
# sweep is exactly what the OLTC CONTROLLER automates against a voltage
# target; the closed loop is chapter 4 of the
# [workshop tour](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_tour.ipynb).
#
# ## Chapter 2: the phase tap (PST) moves power
#
# **Example 2: the phase tap moves power.** Now the loop network every
# PST lives for: a transformer in parallel with
# a line. Without a phase shift the flow splits by impedance and stays
# put. Every phase-tap step tilts the angle across the transformer and
# REROUTES active power between the two parallel paths:
#
# ```text
#        +=== PST (tau, phi) ===+
#   S ---|                      |--- M --- line --- L (load 70 MW)
#        +------- line ---------+
# ```
#
# The asymmetrical PST (the classical Schrägregler) injects its boost
# voltage at a winding angle, here 60 degrees; each step therefore changes
# BOTH the angle and, slightly, the ratio. `calcPhaseTapAngleRatio`
# returns the effective pair for a given step:

function pst_loop(ratio, shift_deg)
  net = Net(name = "pst_loop", baseMVA = 100.0)
  for b in ("S", "M", "L")
    addBus!(net = net, busName = b, vn_kV = 110.0)
  end
  addProsumer!(net = net, busName = "S", type = "EXTERNALNETWORKINJECTION", referencePri = "S", vm_pu = 1.0, va_deg = 0.0)
  addProsumer!(net = net, busName = "L", type = "ENERGYCONSUMER", p = 70.0, q = 20.0)
  addPIModelTrafo!(net = net, fromBus = "S", toBus = "M", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, ratio = ratio, shift_deg = shift_deg, status = 1)
  addPIModelACLine!(net = net, fromBus = "S", toBus = "M", r_pu = 0.03, x_pu = 0.20, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "M", toBus = "L", r_pu = 0.02, x_pu = 0.12, b_pu = 0.0, status = 1)
  validate!(net = net)
  solve!(net)
  return net
end

pst = step -> PhaseTapChangerModel(kind = :asymmetrical, step = step, lowStep = -8, highStep = 8, neutralStep = 0, voltage_step_increment = 0.0125, winding_connection_angle_deg = 60.0)
for step in (-6, 0, 6)
  tap = calcPhaseTapAngleRatio(pst(step))
  net = pst_loop(tap.effective_ratio, tap.effective_shift_deg)
  println("step ", lpad(step, 3), ": phi = ", lpad(round(tap.effective_shift_deg; digits = 2), 6), " deg, tau = ", round(tap.effective_ratio; digits = 4), "  ->  P(trafo) = ", round(get_branch_p_from_to_mw(net, "S", "M"); digits = 1), " MW of 70")
end

# Reading aid (Example 2): at neutral the transformer takes the larger share (it has
# the lower reactance). Six steps of phase shift swing tens of MW from
# one parallel path to the other while the total delivery stays 70 MW
# plus losses: the PST does not produce power, it REROUTES it. The
# `:symmetrical` kind does the same with tau pinned to exactly 1.0 (pure
# angle, no voltage side effect), which is the in-line quadrature booster.
#
# ## Chapter 3: closing the loop, and the X(alpha) subtlety
#
# **Example 3: closing the loop.** In operation nobody dials raw steps: a
# target flow is given and the
# outer control loop moves the tap, here on the same S/M/L parallel-path
# loop as Example 2 (diagram there), with the source raised to 1.02 pu.
# `addPowerTransformerControl!` with
# `control_phase = true` regulates the active power through the
# transformer. One physical subtlety makes the difference to Example 2:
# a real PST's series reactance is not constant, it follows the device
# characteristic $X(\alpha)$ as the boost winding moves. Attach a
# `PhaseTapChangerModel` with `x_min`/`x_max` to the winding and every
# accepted tap move re-stamps the branch with the reactance at the new
# angle:

net = Net(name = "pst_ctrl", baseMVA = 100.0)
for b in ("S", "M", "L")
  addBus!(net = net, busName = b, vn_kV = 110.0)
end
addProsumer!(net = net, busName = "S", type = "EXTERNALNETWORKINJECTION", referencePri = "S", vm_pu = 1.02, va_deg = 0.0)
addProsumer!(net = net, busName = "L", type = "ENERGYCONSUMER", p = 70.0, q = 20.0)
addPIModelTrafo!(net = net, fromBus = "S", toBus = "M", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, ratio = 1.0, shift_deg = 0.0, status = 1)
addPIModelACLine!(net = net, fromBus = "S", toBus = "M", r_pu = 0.03, x_pu = 0.20, b_pu = 0.0, status = 1)
addPIModelACLine!(net = net, fromBus = "M", toBus = "L", r_pu = 0.02, x_pu = 0.12, b_pu = 0.01, status = 1)
tbr = getNetBranch(net = net, fromBus = "S", toBus = "M")
tbr.has_phase_tap = true
tbr.phase_min_deg = -10.0
tbr.phase_max_deg = 10.0
tbr.phase_step_deg = 0.5
## the device characteristic: X grows from x_min at neutral to x_max at
## the range end
net.trafos[1].side1.phase_taps = Sparlectra.PhaseTapChangerModel(kind = :symmetrical, step = 0, lowStep = -10, highStep = 10, neutralStep = 0, voltage_step_increment = 0.01, x_min = 0.08, x_max = 0.16)

run_sparlectra(net = net)
p0 = get_branch_p_from_to_mw(net, "S", "M")
addPowerTransformerControl!(net; trafo = string(tbr.branchIdx), mode = :branch_active_power, target_branch = ("S", "M"), p_target_mw = p0 - 8.0, control_ratio = false, control_phase = true, deadband_p_mw = 4.0)
run_sparlectra(net = net)
println("controller: P(trafo) ", round(p0; digits = 1), " -> ", round(get_branch_p_from_to_mw(net, "S", "M"); digits = 1), " MW (target ", round(p0 - 8.0; digits = 1), ")")
println("            settled at phi = ", tbr.phase_shift_deg, " deg with x_pu = ", round(tbr.x_pu; digits = 4), " (device characteristic, not the 0.08 import value)")

# Reading aid (Example 3): the converged reactance sits ABOVE the import-time 0.08 pu
# because the controller left neutral and the branch follows
# $X(\alpha)$. Without the typed model the loop would steer the same
# target with a slightly wrong (static) reactance; the runnable
# comparison of both variants is `exp_pst_reactance_coupling.jl` in the
# repository.
#
# ## Chapter 4: the three-winding transformer
#
# **Example 4: the three-winding transformer.** A 3WT is solved as its
# STAR EQUIVALENT: one auxiliary bus in the middle
# and three two-winding legs, one per winding. A tap changer always sits
# on ONE winding, so in the equivalent it lands on one leg. The test
# network hangs the star between a slack feeder and an MV load:
#
# ```text
#   B1 (slack) ---- B2 (380 kV, HV) ==leg 1==+
#                                            |
#   B5 (load) ----- B3 (110 kV, MV) ==leg 2==+ AUX (star point)
#                                            |
#                   B4 ( 20 kV, LV) ==leg 3==+
#
#   ----  AC line    ==leg==  2WT leg of the star equivalent
#   B5 carries the 80 MW load, B4 a 5 MVAr shunt
# ```

function build_3wt(; oltc_step::Int)
  net = Net(name = "3wt_demo", baseMVA = 1000.0)
  addBus!(net = net, busName = "B1", vn_kV = 380.0)
  addBus!(net = net, busName = "B2", vn_kV = 380.0)
  addBus!(net = net, busName = "B3", vn_kV = 110.0)
  addBus!(net = net, busName = "B4", vn_kV = 20.0)
  addBus!(net = net, busName = "B5", vn_kV = 110.0)
  addProsumer!(net = net, busName = "B1", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "B1")
  addACLine!(net = net, fromBus = "B1", toBus = "B2", length = 1.0, r = 0.01, x = 0.10)
  ## the star equivalent by hand: AUX bus plus three 2WT legs; the OLTC
  ## ratio (device formula) is applied to the HV leg only
  taps = PowerTransformerTaps(Vn_kV = 380.0, step = oltc_step, lowStep = -9, highStep = 9, neutralStep = 0, voltageIncrement_kV = 3.8)
  addBus!(net = net, busName = "AUX", vn_kV = 380.0, isAux = true)
  add2WTPIModelTrafo!(net = net, fromBus = "AUX", toBus = "B2", side = 1, r = 0.20, x = 4.00, b = 0.0, status = 1, ratedU = 380.0, ratedS = 1000.0, ratio = calcRatioTapCorrection(taps), shift_deg = 0.0)
  add2WTPIModelTrafo!(net = net, fromBus = "AUX", toBus = "B3", side = 1, r = 0.30, x = 6.00, b = 0.0, status = 1, ratedU = 380.0, ratedS = 500.0, ratio = 1.0, shift_deg = 0.0)
  add2WTPIModelTrafo!(net = net, fromBus = "AUX", toBus = "B4", side = 1, r = 0.40, x = 10.00, b = 0.0, status = 1, ratedU = 380.0, ratedS = 200.0, ratio = 1.0, shift_deg = 0.0)
  addACLine!(net = net, fromBus = "B3", toBus = "B5", length = 1.0, r = 0.01, x = 0.10)
  addProsumer!(net = net, busName = "B5", type = "ENERGYCONSUMER", p = 80.0, q = 30.0)
  addShunt!(net = net, busName = "B4", pShunt = 0.0, qShunt = 5.0)
  validate!(net = net)
  solve!(net)
  return net
end

for step in (0, 5)
  net3wt = build_3wt(oltc_step = step)
  println("OLTC step ", step, " on the HV leg: Vm(B3) = ", bus_vm(net3wt, "B3"), " pu, Vm(B4) = ", bus_vm(net3wt, "B4"), " pu")
end

# Reading aid (Example 4): the tap sits on the HV leg, so BOTH output windings move
# together; a tap on the MV leg would move only B3. One honest note: the
# tap objects that `create3WTWindings!` can carry are not yet wired into
# the Net-building 3WT path, so this notebook computes the correction
# with the device formula and applies it to the leg explicitly, which is
# the same arithmetic the wiring would do. The full three-case study
# (OLTC, PST asymmetrical, PST symmetrical on a 3WT) is
# `exp_3wt_phase_taps.jl` in the repository.
#
# ## Where to go next
#
# - [Transformer Control](https://welthulk.github.io/Sparlectra.jl/transformer_control/):
#   the OLTC and phase-tap controllers, targets, deadbands, and limits.
# - [Branch Model](https://welthulk.github.io/Sparlectra.jl/branchmodel/):
#   how tau and phi enter the admittance stamp.
# - [Workshop tour](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_tour.ipynb),
#   chapter 3: the closed-loop OLTC holding a bus voltage.
