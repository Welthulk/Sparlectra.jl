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

# file: docs/lit/workshop_tour_advanced.jl
# purpose: Literate.jl source of the ADVANCED workshop tour (Expert and
#          Beyond tiers): remote voltage control, HVDC pairing, state
#          estimation, FACTS limits, N-1 contingency batches, and threaded
#          sweeps. Companion of workshop_tour.jl (Newcomer to Advanced).

# # The Sparlectra workshop tour, advanced
#
# [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_tour_advanced.ipynb)
#
# This is the second half of the Sparlectra workshop: it picks up where the
# [basic tour](https://welthulk.github.io/Sparlectra.jl/generated/workshop_tour/)
# (Newcomer to Advanced: first network, model work, slack types and short
# circuit, tap control, Q(U)) leaves off. You should be comfortable
# building a net and reading a power-flow result; everything else is
# introduced here.
#
# After the warm-up the chapters climb two tiers:
#
# **Expert**
# 1. Remote voltage control by a machine
# 2. A steerable HVDC link (back-to-back pairing, incl. meshed operation)
# 3. State estimation
# 4. FACTS devices and their limits (STATCOM vs SVC, switched banks, TCSC vs SSSC)
# 5. N-1 contingency analysis
#
# **Beyond**
# 6. Using your cores: parallel sweeps on Julia threads
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
#nb ## To test another branch, set rev to its name, e.g. rev = "dev/r0.9.12".
#nb ## For the latest registered release use: Pkg.add("Sparlectra")
#nb ## Switching versions in a running session? A "[loaded: ...]" note means
#nb ## the old version is still active: restart the runtime, then rerun
#nb ## this cell.

# ## Warm-up and shared helpers
#
# Julia compiles each function on first use. This one cell warms EVERY
# path the chapters exercise, so nothing stalls mid-tour: the
# Newton-Raphson solver, the HVDC pair controller (chapter 2), the WLS
# state estimator (chapter 3), the contingency batch (chapter 5), and the
# IEC 60909 short circuit (chapter 6 sweeps it under threads). The
# `using` clauses and the small helpers of the whole tour live here too,
# collected up top so they cannot be missed.

using Sparlectra
using Random

## solve helper used by all chapters (25 iterations, tolerance 1e-8)
function solve!(net; kwargs...)
  etime = @elapsed begin
    ite, erg = runpf!(net, 25, 1e-8, 0; kwargs...)
  end
  erg == 0 || error("Power flow did not converge (status = $erg)")
  calcNetLosses!(net)
  return etime, ite
end

## peek into the solved state (chapter 2 reads bus angles with it)
bus_va_deg(net, bus) = net.nodeVec[net.busDict[bus]]._va_deg

## the 7-bus double ring from the basic tour's chapter 1, packed as a
## helper: two chapters here (state estimation, N-1) reuse it
function build_ring7(name::String)
  net = Net(name = name, baseMVA = 100.0)
  addBus!(net = net, busName = "B1", vn_kV = 110.0, vm_pu = 1.02, va_deg = 0.0)
  for i in 2:7
    addBus!(net = net, busName = "B$(i)", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
  end
  addPIModelACLine!(net = net, fromBus = "B1", toBus = "B2", r_pu = 0.010, x_pu = 0.080, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "B2", toBus = "B3", r_pu = 0.011, x_pu = 0.085, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "B3", toBus = "B4", r_pu = 0.012, x_pu = 0.090, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "B4", toBus = "B5", r_pu = 0.010, x_pu = 0.080, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "B5", toBus = "B6", r_pu = 0.011, x_pu = 0.085, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "B6", toBus = "B7", r_pu = 0.012, x_pu = 0.090, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "B7", toBus = "B1", r_pu = 0.010, x_pu = 0.080, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "B2", toBus = "B5", r_pu = 0.009, x_pu = 0.070, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "B3", toBus = "B6", r_pu = 0.009, x_pu = 0.070, b_pu = 0.0, status = 1)
  addProsumer!(net = net, busName = "B1", type = "EXTERNALNETWORKINJECTION", referencePri = "B1", vm_pu = 1.02, va_deg = 0.0)
  addProsumer!(net = net, busName = "B3", type = "GENERATOR", p = 60.0, q = 10.0)
  addProsumer!(net = net, busName = "B2", type = "LOAD", p = 35.0, q = 10.0)
  addProsumer!(net = net, busName = "B4", type = "LOAD", p = 45.0, q = 15.0)
  addProsumer!(net = net, busName = "B5", type = "LOAD", p = 25.0, q = 8.0)
  addProsumer!(net = net, busName = "B6", type = "LOAD", p = 30.0, q = 10.0)
  addProsumer!(net = net, busName = "B7", type = "LOAD", p = 20.0, q = 6.0)
  ok, msg = validate!(net = net)
  ok || error("Network validation failed: $msg")
  return net
end

## tiny warm-up net: a grid connection WITH declared short-circuit data
wnet = Net(name = "warmup", baseMVA = 100.0)
addBus!(net = wnet, busName = "A", vn_kV = 110.0)
addBus!(net = wnet, busName = "B", vn_kV = 110.0)
addExternalGrid!(net = wnet, busName = "A", vm_pu = 1.0, sk_max_MVA = 2000.0, sk_min_MVA = 1500.0, rx_max = 0.1, internal_impedance = false)
addProsumer!(net = wnet, busName = "B", type = "ENERGYCONSUMER", p = 10.0, q = 3.0)
addPIModelACLine!(net = wnet, fromBus = "A", toBus = "B", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)

t_first = @elapsed runpf!(wnet, 10, 1e-8, 0)
t_second = @elapsed runpf!(wnet, 10, 1e-8, 0)
println("power flow     : first solve ", round(t_first; digits = 2), " s (compiles), second ", round(t_second * 1000; digits = 2), " ms")

t_sc = @elapsed runShortCircuit!(wnet; case = :max)
println("short circuit  : ", round(t_sc; digits = 2), " s")

## HVDC pair path: two 2-bus islands coupled by a controller
whv = Net(name = "warmup_hvdc", baseMVA = 100.0)
for b in ("W1", "W2", "W3", "W4")
  addBus!(net = whv, busName = b, vn_kV = 110.0)
end
addProsumer!(net = whv, busName = "W1", type = "EXTERNALNETWORKINJECTION", referencePri = "W1", vm_pu = 1.0, va_deg = 0.0)
addProsumer!(net = whv, busName = "W3", type = "EXTERNALNETWORKINJECTION", referencePri = "W3", vm_pu = 1.0, va_deg = 0.0)
addProsumer!(net = whv, busName = "W2", type = "ENERGYCONSUMER", p = 5.0, q = 1.0)
addProsumer!(net = whv, busName = "W4", type = "ENERGYCONSUMER", p = 5.0, q = 1.0)
addPIModelACLine!(net = whv, fromBus = "W1", toBus = "W2", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
addPIModelACLine!(net = whv, fromBus = "W3", toBus = "W4", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
addProsumer!(net = whv, busName = "W2", type = "GENERATOR", p = -5.0, q = 0.0)
addProsumer!(net = whv, busName = "W4", type = "GENERATOR", p = 5.0, q = 0.0)
addHvdcPairControl!(whv; from_bus = "W2", to_bus = "W4", p_transfer_mw = 5.0)
t_hvdc = @elapsed run_control!(whv; controllers = collect_outer_controllers(whv), pf_config = PowerFlowConfig(method = :rectangular, max_iter = 15, tol = 1e-8), control_config = ControlConfig(max_outer_iterations = 4, trace = false))
println("HVDC control   : ", round(t_hvdc; digits = 2), " s")

## state-estimation path: synthetic measurements plus one WLS run
setMeasurementsFromPF!(wnet; includeVm = true, includePinj = true, includeQinj = true, includePflow = true, includeQflow = true, noise = false)
t_se = @elapsed runse!(wnet; maxIte = 8, tol = 1e-6, flatstart = true, jacEps = 1e-6, updateNet = false)
println("state estimator: ", round(t_se; digits = 2), " s")

## chapter-5 path: one single-case contingency batch on the warm-up net
t_n1 = @elapsed runContingencies!(wnet, generateN1Branches(wnet))
println("contingency    : ", round(t_n1; digits = 2), " s, everything warm")

# ## Part I: Expert
#
# ## Chapter 1: remote voltage control by a machine
#
# A machine regulates the voltage at a **different** bus via its reactive
# output, the counterpart of a CGMES `RegulatingControl` at a foreign
# terminal. The outer loop drives the machine's Q until the remote target
# is met, and parks honestly `at_limit` when the reactive range is too
# small. Details:
# [Remote Voltage Control](https://welthulk.github.io/Sparlectra.jl/remote_voltage_control/).

function build_rvc(qmin_mvar::Float64, qmax_mvar::Float64)
  net = Net(name = "tour_rvc", baseMVA = 100.0)
  addBus!(net = net, busName = "Slack", vn_kV = 110.0)
  addBus!(net = net, busName = "GenBus", vn_kV = 110.0)
  addBus!(net = net, busName = "Load", vn_kV = 110.0)
  addProsumer!(net = net, busName = "Slack", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.02, va_deg = 0.0, referencePri = "Slack")
  addProsumer!(net = net, busName = "GenBus", type = "SYNCHRONOUSMACHINE", p = 30.0, q = 0.0, qMin = qmin_mvar, qMax = qmax_mvar)
  addProsumer!(net = net, busName = "Load", type = "ENERGYCONSUMER", p = -70.0, q = -20.0)
  addPIModelACLine!(net = net, fromBus = "Slack", toBus = "GenBus", r_pu = 0.02, x_pu = 0.12, b_pu = 0.01, status = 1)
  addPIModelACLine!(net = net, fromBus = "GenBus", toBus = "Load", r_pu = 0.02, x_pu = 0.12, b_pu = 0.01, status = 1)
  return net
end

pf = PowerFlowConfig(max_iter = 30, tol = 1e-9)

net_rvc = build_rvc(-50.0, 50.0)
runpf!(net_rvc; config = pf, verbose = 0)
println("uncontrolled: Vm(Load) = ", round(get_bus_vm_pu(net_rvc, "Load"); digits = 4), " pu")

addMachineVoltageControl!(net_rvc; bus = "GenBus", target_bus = "Load", target_vm_pu = 1.05, deadband_vm_pu = 5e-4)
result = run_control!(net_rvc; controllers = collect_outer_controllers(net_rvc), pf_config = pf, control_config = ControlConfig(max_outer_iterations = 15), verbose = 0)
println("controlled:   Vm(Load) = ", round(get_bus_vm_pu(net_rvc, "Load"); digits = 4), " pu (target 1.05)")
println("loop: status = ", result.status, ", outer iterations = ", result.outer_iterations, ", pf solves = ", result.powerflow_solves)
printMachineControllerSummary(stdout, net_rvc)

# The honest failure mode: cut the reactive range to ±2 MVAr and the same
# target is out of reach. The controller parks at its limit and says so
# instead of pretending convergence.

net_rvc2 = build_rvc(-2.0, 2.0)
runpf!(net_rvc2; config = pf, verbose = 0)
addMachineVoltageControl!(net_rvc2; bus = "GenBus", target_bus = "Load", target_vm_pu = 1.05, deadband_vm_pu = 5e-4)
result2 = run_control!(net_rvc2; controllers = collect_outer_controllers(net_rvc2), pf_config = pf, control_config = ControlConfig(max_outer_iterations = 15), verbose = 0)
println("limited:      Vm(Load) = ", round(get_bus_vm_pu(net_rvc2, "Load"); digits = 4), " pu (target 1.05)")
printMachineControllerSummary(stdout, net_rvc2)

# ## Chapter 2: a steerable HVDC link (back-to-back pairing)
#
# Two AC areas joined ONLY by an HVDC converter pair: no AC tie, no angle
# coupling, so the areas stay two separate electrical islands with their
# own references. The transfer through the link is a control setpoint, not
# the result of an angle difference. First the Stage-0 view: two fixed
# injections reproduce a snapshot of 80 MW transfer with 4 MW converter
# loss.
#
# Note that EACH island keeps a classical reference of its own (`A1` and
# `C1`); the converters at `A2`/`C2` are plain injections the controller
# will steer. The result header therefore reports `Slack: 2`, and both
# reference buses appear as `SLACK` rows in the bus table:
#
# ```text
#      island A                          island C
#   A1 -------- A2  ===== DC link =====  C2 -------- C1
# (slack)    load 40 MW              load 50 MW    (slack)
#            + converter             + converter
#            (from side)             (to side)
# ```

function build_b2b(name::String)
  net = Net(name = name, baseMVA = 100.0)
  for b in ("A1", "A2", "C1", "C2")
    addBus!(net = net, busName = b, vn_kV = 380.0)
  end
  addPIModelACLine!(net = net, fromBus = "A1", toBus = "A2", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "C1", toBus = "C2", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
  addProsumer!(net = net, busName = "A1", type = "EXTERNALNETWORKINJECTION", referencePri = "A1", vm_pu = 1.0, va_deg = 0.0)
  addProsumer!(net = net, busName = "C1", type = "EXTERNALNETWORKINJECTION", referencePri = "C1", vm_pu = 1.0, va_deg = 0.0)
  addProsumer!(net = net, busName = "A2", type = "ENERGYCONSUMER", p = 40.0, q = 10.0)
  addProsumer!(net = net, busName = "C2", type = "ENERGYCONSUMER", p = 50.0, q = 12.0)
  addProsumer!(net = net, busName = "A2", type = "GENERATOR", p = -80.0, q = 0.0)  ## converter, exports
  addProsumer!(net = net, busName = "C2", type = "GENERATOR", p = 76.0, q = 0.0)   ## converter, receives 80 - 4
  return net
end

net6 = build_b2b("tour_b2b")
## register the hand-built link so the result tables report it (importers
## and addHvdcPairControl! do this automatically)
addHvdcLink!(net6; from_bus = "A2", to_bus = "C2")
etime, ite = solve!(net6; islands_enabled = true)
println("two islands solved in ", ite, " iteration(s)")
printACPFlowResults(net6, etime, ite, 1e-8)

# Reading aid: both areas balance on their own reference; the link carries
# whatever the snapshot says. To make it steerable, pair the two converter
# injections: `addHvdcPairControl!` enforces the invariant
# $P_\text{to} = P_\text{transfer} - P_\text{loss}$ exactly and lets you
# retarget the transfer.

addHvdcPairControl!(net6; from_bus = "A2", to_bus = "C2", p_transfer_mw = 120.0, loss_mw = 4.0, p_rating_mw = 150.0)
result6 = run_control!(net6; controllers = collect_outer_controllers(net6), pf_config = PowerFlowConfig(method = :rectangular, max_iter = 25, tol = 1e-8), control_config = ControlConfig(max_outer_iterations = 8, trace = false))
calcNetLosses!(net6)
printACPFlowResults(net6, etime, result6.last_pf_iterations, 1e-8)

# Reading aid: the HVDC pair reports inside the `Control` section of the
# classical result output, in the same aligned label/value layout as the
# transformer, machine, and TCSC summaries (`printHvdcPairControllerSummary`
# prints the same block standalone). The link itself has its own `HVDC
# Link Flows` table right after the link table: ordered transfer,
# delivered power, loss, and the controller status per link. And look at
# the LAST row of the branch table: the link appears there too, typed
# `Link` and marked "HVDC, not a branch", so the topology reads in one
# table; it looks like a branch (two buses, power moves), but no Y-bus
# element exists behind it, Q is zero, and its Pv column is the converter
# loss. The link now carries 120 MW instead of 80: the line A1 -> A2
# supplies 40 MW more (area A's reference generates the export), while
# C1 -> C2 turns around and carries the received power away from the
# converter bus.
#
# **Why does the solver still report two islands?** An AC voltage angle is
# only defined *within* one synchronous island, relative to that island's
# own reference. The link transfers power but no angle information (there
# is no branch, no admittance, no angle coupling between the areas), so
# each island keeps its own reference pinned at 0 degrees. The two-island
# report is the model telling you the areas are asynchronous; it would be
# wrong for it to disappear. The peek below makes that visible: both
# reference buses sit at exactly 0.0 deg, and comparing an A-side angle
# with a C-side angle carries no information, because each is measured
# against a different zero.

## bus_va_deg comes from the warm-up cell (shared helpers up top)
for (bus, role) in (("A1", "reference of island A"), ("A2", "converter, exports 120 MW"), ("C1", "reference of island C"), ("C2", "converter, receives 116 MW"))
  println(rpad(bus, 4), rpad(role, 27), ": Vm = ", round(get_bus_vm_pu(net6, bus); digits = 4), " pu, Va = ", round(bus_va_deg(net6, bus); digits = 3), " deg")
end

# Reading aid: within island A the angle falls toward A2 (the converter
# bus imports 120 MW plus the local load from the reference), within
# island C it rises toward C2 (the converter bus feeds the island). Each
# gradient is meaningful only against its own 0-degree reference. The same
# controller attaches automatically on import when a MATPOWER case sets
# `matpower_dcline_mode = paired_control` or a CGMES delivery is loaded
# with `hvdc_mode = paired_control`. Theory:
# [HVDC Back-to-Back](https://welthulk.github.io/Sparlectra.jl/hvdc_back_to_back/).

# ### HVDC as the island's source (grid-forming)
#
# Can the converter itself BE the reference (slack) of island C? Not as
# one side of the paired controller: the pairing treats the transfer as a
# *setpoint* on both sides, while a slack's power is the *outcome* of its
# island's balance (load plus losses). One injection cannot be both at
# once, so `addHvdcPairControl!` refuses a reference bus by design.
#
# The converse is a perfectly valid model though, called a grid-forming
# (Vf) converter: the receiving converter IS the island's source. It
# holds voltage and angle at its PCC, and its power output follows from
# whatever the island draws. Think of an offshore platform or an
# asynchronously supplied island grid.
#
# A word on terms: "slack" is the solver's name for an island's reference
# node, and every island has exactly one, however it is modeled. The ideal
# slack of the basic tour's chapter 3, the external-grid SOURCE behind an impedance, and
# the grid-forming converter here are three MODELS of that one reference.
# The result output keeps them apart: `SOURCE` in the bus table (and
# `Source: m` in the header) is reserved for the external-grid feeder
# element, because there the reference voltage sits BEHIND an impedance.
# The grid-forming converter is an ideal reference directly at its PCC, so
# its row honestly reads `SLACK`; its role is marked in the `Control`
# column instead: `B2B src` for the grid-forming reference, `B2B` for a
# steered converter injection. Island C below has NO source of its own;
# its reference moves onto the converter bus C2 (compare the sketch with
# the setpoint variant above):
#
# ```text
#      island A                          island C
#   A1 -------- A2  ===== DC link =====  C2 -------- C1
# (slack)    load 40 MW            grid-forming    load 50 MW
#            + sending             converter
#            converter             (= island C reference)
# ```

function build_b2b_source(name::String; sending_mw::Float64 = 0.0)
  net = Net(name = name, baseMVA = 100.0)
  for b in ("A1", "A2", "C1", "C2")
    addBus!(net = net, busName = b, vn_kV = 380.0)
  end
  addPIModelACLine!(net = net, fromBus = "A1", toBus = "A2", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "C2", toBus = "C1", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
  addProsumer!(net = net, busName = "A1", type = "EXTERNALNETWORKINJECTION", referencePri = "A1", vm_pu = 1.0, va_deg = 0.0)
  addProsumer!(net = net, busName = "A2", type = "ENERGYCONSUMER", p = 40.0, q = 10.0)
  ## island C: no classical slack. The receiving converter is the
  ## grid-forming source, holding 1.0 pu / 0 deg at its PCC bus C2.
  addProsumer!(net = net, busName = "C2", type = "EXTERNALNETWORKINJECTION", referencePri = "C2", vm_pu = 1.0, va_deg = 0.0)
  addProsumer!(net = net, busName = "C1", type = "ENERGYCONSUMER", p = 50.0, q = 12.0)
  ## sending-side converter; mirrored below once the island balance is known
  addProsumer!(net = net, busName = "A2", type = "GENERATOR", p = sending_mw, q = 0.0)
  return net
end

net7 = build_b2b_source("tour_b2b_source")
solve!(net7; islands_enabled = true)
p_island_c = get_branch_p_from_to_mw(net7, "C2", "C1")
println("grid-forming converter at C2 delivers ", round(p_island_c; digits = 3), " MW (island C load + line loss)")
for bus in ("C2", "C1")
  println("  ", bus, ": Vm = ", round(get_bus_vm_pu(net7, bus); digits = 4), " pu, Va = ", round(bus_va_deg(net7, bus); digits = 3), " deg")
end

# Reading aid: the transfer is no longer a setpoint. Island C decides how
# much it draws (load plus line loss), the converter delivers exactly
# that, and the island's reference sits at the converter PCC: 1.0 pu and
# 0 degrees at C2, while the load bus C1 hangs below it. The sending side
# must now *mirror* the island draw plus the converter loss:

net8 = build_b2b_source("tour_b2b_mirrored"; sending_mw = -(p_island_c + 4.0))
addHvdcLink!(net8; from_bus = "A2", to_bus = "C2")  ## Stage-0 record for the result tables
etime8, ite8 = solve!(net8; islands_enabled = true)
println("sending side A1 -> A2 carries ", round(get_branch_p_from_to_mw(net8, "A1", "A2"); digits = 3), " MW (40 MW local load + ", round(p_island_c + 4.0; digits = 3), " MW export + line loss)")
printACPFlowResults(net8, etime8, ite8, 1e-8)

# Reading aid: compare the SLACK rows with the setpoint variant. The
# references are now `A1` and `C2`: the receiving converter itself is
# island C's reference, its `Pg` is the island draw, and `C1` carries only
# load.

# And the refusal from above, demonstrated: pairing the grid-forming
# converter is rejected, because its injection is the island balance, not
# a setpoint:

try
  addHvdcPairControl!(net8; from_bus = "A2", to_bus = "C2", p_transfer_mw = 50.0)
catch err
  println(sprint(showerror, err))
end

# The pairing controller automates exactly this mirror: `mode =
# :island_feed` reads the island balance after each solve and keeps the
# sending side matched, with an honest `at_limit` once the island draw
# exceeds `p_rating_mw`. No transfer setpoint is given, the island
# decides:

net9 = build_b2b_source("tour_b2b_grid_forming")
addHvdcPairControl!(net9; from_bus = "A2", to_bus = "C2", mode = :island_feed, loss_mw = 4.0, p_rating_mw = 150.0)
result9 = run_control!(net9; controllers = collect_outer_controllers(net9), pf_config = PowerFlowConfig(method = :rectangular, max_iter = 25, tol = 1e-8), control_config = ControlConfig(max_outer_iterations = 8, trace = false))
calcNetLosses!(net9)
printACPFlowResults(net9, etime, result9.last_pf_iterations, 1e-8)

# Reading aid, and the direct comparison with the setpoint variant above.
# Both versions report `Slack: 2` in the header (one reference per island,
# always), but WHERE the references sit is the whole difference:
#
# - **Setpoint variant:** references at `A1` and `C1`, each island brings
#   its own source. Both converters are plain PQ injections steered by the
#   controller (`-120 / +116 MW` in the bus table), and the transfer is an
#   order the link follows. Look at C1's `Pg`: it is negative, island C's
#   own slack absorbs the surplus the link pumps in.
# - **Grid-forming variant (this table):** references at `A1` and `C2`.
#   The receiving converter itself is island C's `SLACK` row, marked
#   `B2B src` in the `Control` column; its `Pg` column is the island
#   balance outcome (load plus line loss, no setpoint anywhere), `C1`
#   carries only load, and the HVDC block in the `Control` section shows
#   the transfer as "mirrored from island draw". The `HVDC Link Flows`
#   table lists the same link with `mode = island_feed` and the mirrored
#   transfer in its `P_from` column.
#
# Try `p_rating_mw = 40.0`: the sending side pins at the rating with
# `at_limit = true` and `converged = false`. The island's reference still
# balances in the model (a power flow cannot show the collapse), so the
# honest flag is what marks the undeliverable draw.

# ### Variant: grid-forming with droop (SOURCE model)
#
# The ideal reference above holds exactly 1.0 pu at the PCC no matter what
# the island draws. A real VSC has finite control stiffness. The
# external-grid element from the basic tour's chapter 3 models exactly that: the reference
# voltage sits BEHIND the impedance $Z_Q = U_n^2 / S_k''$, so the PCC
# voltage droops under load. Declaring the converter that way finally
# makes the bus table say `SOURCE`, and the header counts
# `Slack: 1 Source: 1`:

function build_b2b_droop(name::String; sending_mw::Float64 = 0.0, sk_mva::Float64 = 800.0)
  net = Net(name = name, baseMVA = 100.0)
  for b in ("A1", "A2", "C1", "C2")
    addBus!(net = net, busName = b, vn_kV = 380.0)
  end
  addPIModelACLine!(net = net, fromBus = "A1", toBus = "A2", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "C2", toBus = "C1", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
  addProsumer!(net = net, busName = "A1", type = "EXTERNALNETWORKINJECTION", referencePri = "A1", vm_pu = 1.0, va_deg = 0.0)
  addProsumer!(net = net, busName = "A2", type = "ENERGYCONSUMER", p = 40.0, q = 10.0)
  addProsumer!(net = net, busName = "C1", type = "ENERGYCONSUMER", p = 50.0, q = 12.0)
  ## the grid-forming converter as a non-ideal source: reference voltage
  ## behind Z_Q, declared by its short-circuit power like a real feeder
  addExternalGrid!(net = net, busName = "C2", vm_pu = 1.0, sk_max_MVA = sk_mva, sk_min_MVA = sk_mva, rx_max = 0.1, internal_impedance = true)
  addProsumer!(net = net, busName = "A2", type = "GENERATOR", p = sending_mw, q = 0.0)
  return net
end

net10 = build_b2b_droop("tour_b2b_droop"; sending_mw = -(p_island_c + 4.0))
etime10, ite10 = solve!(net10; islands_enabled = true)
printACPFlowResults(net10, etime10, ite10, 1e-8)

# Reading aid: the header reads `Slack: 1 Source: 1`, the grid connection
# line names the source with its feeder data, and the hidden anchor bus
# behind the impedance is the `SOURCE` row. The PCC bus `C2` is now a
# plain PQ bus whose voltage sags below 1.0 pu under the island load: that
# sag is the droop, and its size follows from the declared `sk_max_MVA`
# (stiffer converter = higher $S_k''$ = less droop). Pairing this
# source-model reference with the island_feed controller is a possible
# follow-up of the pairing controller.

# ### Meshed operation: the two areas get an AC tie
#
# So far the link was the ONLY connection. Now close an AC branch between
# `A1` and `C1`: the two areas become one synchronous island, and one
# synchronous island carries exactly ONE angle reference. The link
# transfers power, the tie transfers the angle.
#
# ```text
#      one synchronous island
#   A1 -------- A2  ===== DC link =====  C2 -------- C1
#   |                                                 |
#   +------------------- AC tie ---------------------+
# ```

function build_meshed(name::String; c1_model::Symbol)
  net = Net(name = name, baseMVA = 100.0)
  for b in ("A1", "A2", "C1", "C2")
    addBus!(net = net, busName = b, vn_kV = 380.0)
  end
  addPIModelACLine!(net = net, fromBus = "A1", toBus = "A2", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "C1", toBus = "C2", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "A1", toBus = "C1", r_pu = 0.02, x_pu = 0.16, b_pu = 0.0, status = 1)
  addProsumer!(net = net, busName = "A1", type = "EXTERNALNETWORKINJECTION", referencePri = "A1", vm_pu = 1.0, va_deg = 0.0)
  if c1_model === :reference
    ## scene 1: island C keeps its old reference although the tie closed
    addProsumer!(net = net, busName = "C1", type = "EXTERNALNETWORKINJECTION", referencePri = "C1", vm_pu = 1.0, va_deg = 0.0)
  else
    ## scenes 2+: demoted to a voltage-regulated generator (PV)
    addProsumer!(net = net, busName = "C1", type = "GENERATOR", p = 20.0, q = 0.0, vm_pu = 1.0, isRegulated = true)
  end
  addProsumer!(net = net, busName = "A2", type = "ENERGYCONSUMER", p = 40.0, q = 10.0)
  addProsumer!(net = net, busName = "C2", type = "ENERGYCONSUMER", p = 50.0, q = 12.0)
  addProsumer!(net = net, busName = "A2", type = "GENERATOR", p = -80.0, q = 0.0)  ## converter, exports
  addProsumer!(net = net, busName = "C2", type = "GENERATOR", p = 76.0, q = 0.0)   ## converter, receives
  return net
end

# Keeping BOTH old references fails fast, with the buses named. The solver
# never demotes a reference on its own; you decide which one survives:

netm = build_meshed("tour_meshed_two_refs"; c1_model = :reference)
try
  solve!(netm; islands_enabled = true)
catch err
  println(sprint(showerror, err))
end

# Reading aid: this is the same one-reference-per-island rule from the
# beginning of the chapter, now seen from the other side. Demote `C1` to a
# voltage-regulated generator and the meshed net solves; the pair keeps
# its setpoint and the tie carries the balance:

netm2 = build_meshed("tour_meshed"; c1_model = :pv)
addHvdcPairControl!(netm2; from_bus = "A2", to_bus = "C2", p_transfer_mw = 120.0, loss_mw = 4.0, p_rating_mw = 150.0)
resultm = run_control!(netm2; controllers = collect_outer_controllers(netm2), pf_config = PowerFlowConfig(method = :rectangular, max_iter = 25, tol = 1e-8), control_config = ControlConfig(max_outer_iterations = 8, trace = false))
calcNetLosses!(netm2)
printACPFlowResults(netm2, etime, resultm.last_pf_iterations, 1e-8)

# Reading aid: ONE island, ONE `SLACK` row (`A1`), `C1` is an ordinary PV
# generator now. The branch table shows the parallel paths side by side:
# three real AC branches plus the `Link` row marked "HVDC, not a branch".
# The `HVDC Link Flows` table still shows the ordered 120 MW with 4 MW
# loss, and the tie `A1 -> C1` carries whatever the link over- or
# under-delivers relative to what area C draws. Retarget the pair and
# watch the exchange move between link and tie:

ctrlm = only(collect_outer_controllers(netm2))
for target in (120.0, 40.0)
  ctrlm.p_transfer_mw = target
  ctrlm.p_applied = false
  run_control!(netm2; controllers = collect_outer_controllers(netm2), pf_config = PowerFlowConfig(method = :rectangular, max_iter = 25, tol = 1e-8), control_config = ControlConfig(max_outer_iterations = 8, trace = false))
  calcNetLosses!(netm2)
  println("transfer ", target, " MW ordered: AC tie A1 -> C1 carries ", round(get_branch_p_from_to_mw(netm2, "A1", "C1"); digits = 3), " MW")
end

# Reading aid: the pair keeps its order exactly (that is what a setpoint
# means), so every retarget shows up one-to-one in the tie flow. And
# `mode = :island_feed`? A grid-forming converter inside a synchronous
# grid is a different device model and out of scope: with the tie closed
# the registration is rejected by the same one-reference rule, and a
# demoted reference afterwards makes the controller report
# `invalid_topology` instead of silently changing modes.

# ## Chapter 3: state estimation
#
# Close the loop: solve a reference power flow on the basic tour's 7-bus ring (build_ring7, see the warm-up),
# derive a noisy synthetic measurement set from it, check observability,
# and let the WLS estimator reconstruct the state. The full narrative is
# the
# [state-estimation notebook](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_state_estimation.ipynb).

net_se = build_ring7("tour_se")
ite_pf, status_pf = runpf!(net_se, 40, 1e-10, 0)
status_pf == 0 || error("Power flow did not converge")

std = measurementStdDevs(vm = 1e-3, pinj = 1.0, qinj = 1.0, pflow = 0.7, qflow = 0.7)
setMeasurementsFromPF!(
  net_se;
  includeVm = true,
  includePinj = true,
  includeQinj = true,
  includePflow = true,
  includeQflow = true,
  noise = true,
  stddev = std,
  rng = MersenneTwister(42),
)

gobs = evaluate_global_observability(net_se; flatstart = true, jacEps = 1e-6)
println("observability: ", gobs.quality, " (", gobs.n_measurements, " measurements, ", gobs.n_states, " states)")

se = runse!(net_se; maxIte = 12, tol = 1e-6, flatstart = true, jacEps = 1e-6, updateNet = true)
println("SE converged: ", se.converged, " in ", se.iterations, " iterations")
println("objective J:  ", round(se.objectiveJ; digits = 2), " (dof ", se.dof, ", within 3σ: ", se.jWithin3Sigma, ")")
for (name, idx) in sort(collect(net_se.busDict); by = last)
  v = se.voltages[idx]
  println(rpad(name, 4), "  Vm = ", round(abs(v); digits = 4), " pu   Va = ", round(rad2deg(angle(v)); digits = 3), "°")
end

# Reading aid: with mild noise the estimate reproduces the reference state
# to a few 1e-3 pu, and $J$ lands near the degrees of freedom, the
# textbook health check for a WLS estimator.
#
# ## Chapter 4: FACTS devices and their limits
#
# FACTS devices use power electronics to control voltage and flow. The
# workshop has met two members already: the phase-shifting tap (basic
# tour, chapter 4's family) and the HVDC pair (chapter 2). This chapter is about the
# LIMIT characteristics, because that is where the devices actually
# differ. In range, every shunt compensator holds its voltage target the
# same way; at the limit:
#
# - a classical machine keeps its constant reactive box `[Qmin, Qmax]`,
# - a STATCOM is current-limited, it delivers $Q = V \cdot S_{max}$
#   (LINEAR collapse under a sag),
# - an SVC is susceptance-limited, it delivers $Q = V^2 \cdot B$
#   (QUADRATIC collapse).
#
# The ranking shows exactly under the depressed voltage the compensator
# was installed for. We build one weak corridor that sags to about
# 0.92 pu and give all three devices the SAME 10-MVAr rating at 1.0 pu:

facts_rating = 10.0
function build_sag_corridor(name::String; with_machine::Bool)
  cnet = Net(name = name, baseMVA = 100.0)
  for bus in ("Slack", "Mid", "Load")
    addBus!(net = cnet, busName = bus, vn_kV = 110.0)
  end
  addProsumer!(net = cnet, busName = "Slack", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "Slack")
  addProsumer!(net = cnet, busName = "Load", type = "LOAD", p = 60.0, q = 25.0)
  with_machine && addProsumer!(net = cnet, busName = "Mid", type = "GENERATOR", p = 0.0, q = 0.0)
  addPIModelACLine!(net = cnet, fromBus = "Slack", toBus = "Mid", r_pu = 0.02, x_pu = 0.20, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = cnet, fromBus = "Mid", toBus = "Load", r_pu = 0.02, x_pu = 0.20, b_pu = 0.0, status = 1)
  ok, msg = validate!(net = cnet)
  ok || error("corridor net invalid: $msg")
  return cnet
end

## classical machine: constant reactive box, the outer loop parks at_limit
box_net = build_sag_corridor("tour_facts_box"; with_machine = true)
addMachineVoltageControl!(box_net; bus = "Mid", target_bus = "Load", target_vm_pu = 1.0, qmin_mvar = -facts_rating, qmax_mvar = facts_rating)
run_control!(box_net)
box_ctrl = only([c for c in box_net.machineControls if c isa MachineVoltageControl])

## STATCOM: the SAME controller with a converter rating instead of the box;
## the bound Q_lim = V * S_max is refreshed from the solved terminal
## voltage every outer iteration, so the delivered Q tracks the sag
st_net = build_sag_corridor("tour_facts_statcom"; with_machine = true)
addMachineVoltageControl!(st_net; bus = "Mid", target_bus = "Load", target_vm_pu = 1.0, s_max_mva = facts_rating)
run_control!(st_net)
st_ctrl = only([c for c in st_net.machineControls if c isa MachineVoltageControl])
v_st = get_bus_vm_pu(st_net, "Mid")

## SVC: continuous susceptance; at the clamp the Y-bus stamp makes the
## delivered Q follow V^2 all by itself
svc_net = build_sag_corridor("tour_facts_svc"; with_machine = false)
addShuntVoltageControl!(svc_net; bus = "Mid", target_vm_pu = 1.0, bs_min_mvar = -facts_rating, bs_max_mvar = facts_rating)
run_control!(svc_net)
svc_ctrl = only([c for c in svc_net.machineControls if c isa ShuntVoltageControl])
v_svc = get_bus_vm_pu(svc_net, "Mid")

println("all three at their capacitive limit, rated ", facts_rating, " MVAr at 1.0 pu:")
println("  machine box : Q = ", round(box_ctrl.q_mvar; digits = 2), " MVAr (constant)")
println("  STATCOM     : Q = ", round(st_ctrl.q_mvar; digits = 2), " MVAr = V*S_max at V = ", round(v_st; digits = 4), " pu (", round(100 * st_ctrl.q_mvar / facts_rating; digits = 1), " % of rating)")
println("  SVC         : Q = ", round(v_svc^2 * svc_ctrl.bs_mvar; digits = 2), " MVAr = V^2*B at V = ", round(v_svc; digits = 4), " pu (", round(100 * v_svc^2 * svc_ctrl.bs_mvar / facts_rating; digits = 1), " % of rating)")

# Most voltage control in real grids is not an SVC but a mechanically
# SWITCHED bank (MSC/MSR): whole capacitor blocks go in or out, nothing in
# between. The same controller models that with `step_mvar` (issue #324).
# Two design rules make a switched bank usable in a solver loop: it moves
# in WHOLE blocks truncated toward the target (approaching from one side,
# never overshooting, so it cannot hunt between two adjacent steps), and
# when no whole block improves the voltage further it PARKS on the reached
# step, deliberately the last one BEFORE crossing the target
# (under-compensation instead of a possible overvoltage):

## the same corridor and target once with a continuous SVC, once as a bank
cont_net = build_sag_corridor("tour_facts_cont"; with_machine = false)
addShuntVoltageControl!(cont_net; bus = "Mid", target_vm_pu = 0.95, bs_min_mvar = -40.0, bs_max_mvar = 40.0)
run_control!(cont_net)
cont_ctrl = only([c for c in cont_net.machineControls if c isa ShuntVoltageControl])
println("continuous SVC: Bs = ", round(cont_ctrl.bs_mvar; digits = 2), " MVAr, Vm(Mid) = ", round(get_bus_vm_pu(cont_net, "Mid"); digits = 4), " pu, status = ", cont_ctrl.status)

msc_net = build_sag_corridor("tour_facts_msc"; with_machine = false)
addShuntVoltageControl!(msc_net; bus = "Mid", target_vm_pu = 0.95, bs_min_mvar = -40.0, bs_max_mvar = 40.0, step_mvar = 10.0)
msc_res = run_control!(msc_net)
msc_ctrl = only([c for c in msc_net.machineControls if c isa ShuntVoltageControl])
println("MSC bank (4 x 10 MVAr blocks): Bs = ", msc_ctrl.bs_mvar, " MVAr (block ", round(Int, msc_ctrl.bs_mvar / 10.0), " of 4), Vm(Mid) = ", round(get_bus_vm_pu(msc_net, "Mid"); digits = 4), " pu")
println("  status = ", msc_ctrl.status, " after ", msc_res.outer_iterations, " outer iterations (parked on the step grid, no hunting)")

# Reading aid: the continuous SVC settles wherever the secant sends it
# (about 24.7 MVAr here) and converges; the bank can only offer 20 or 30
# and stops at 20, the last block before crossing the 0.95-pu target.
# `status = :parked` is the bank's third state next to `converged` and
# `at_limit`: the step resolution, not the rating, ended the movement.

# The series side has the same split. A TCSC owns a FIXED reactance window
# and keeps it at any loading; an SSSC injects a series voltage, so its
# usable window $|x - x_{base}| \le V_{inj,max}/|I|$ SHRINKS with the
# branch current: it saturates exactly at high transfer. Same two-corridor
# loop, same 35-MW target, once per device:

function build_facts_loop(name::String)
  lnet = Net(name = name, baseMVA = 100.0)
  for bus in ("A", "M1", "M2", "B")
    addBus!(net = lnet, busName = bus, vn_kV = 110.0)
  end
  addProsumer!(net = lnet, busName = "A", type = "EXTERNALNETWORKINJECTION", referencePri = "A", vm_pu = 1.0, va_deg = 0.0)
  addProsumer!(net = lnet, busName = "B", type = "ENERGYCONSUMER", p = 80.0, q = 20.0)
  addPIModelACLine!(net = lnet, fromBus = "A", toBus = "M1", r_pu = 0.01, x_pu = 0.10, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = lnet, fromBus = "M1", toBus = "B", r_pu = 0.01, x_pu = 0.10, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = lnet, fromBus = "A", toBus = "M2", r_pu = 0.02, x_pu = 0.20, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = lnet, fromBus = "M2", toBus = "B", r_pu = 0.02, x_pu = 0.20, b_pu = 0.0, status = 1)
  ok, msg = validate!(net = lnet)
  ok || error("loop net invalid: $msg")
  return lnet
end

tcsc_net = build_facts_loop("tour_facts_tcsc")
tcsc_ctrl = addSeriesReactanceControl!(tcsc_net; fromBus = "A", toBus = "M2", p_target_mw = 35.0, x_min_pu = 0.02, x_max_pu = 0.30)
run_control!(tcsc_net)
println("TCSC window 0.02..0.30 pu : P = ", round(tcsc_ctrl.achieved_p_mw; digits = 2), " MW at x = ", round(tcsc_ctrl.x_pu; digits = 4), " pu, converged = ", tcsc_ctrl.converged)

sssc_net = build_facts_loop("tour_facts_sssc")
sssc_ctrl = addSeriesReactanceControl!(sssc_net; fromBus = "A", toBus = "M2", p_target_mw = 35.0, v_inj_max_pu = 0.01)
run_control!(sssc_net)
println("SSSC V_inj,max 0.01 pu    : P = ", round(sssc_ctrl.achieved_p_mw; digits = 2), " MW at x = ", round(sssc_ctrl.x_pu; digits = 4), " pu, at_limit = ", sssc_ctrl.at_limit)
println("  live window [", round(sssc_ctrl.x_min_pu; digits = 4), ", ", round(sssc_ctrl.x_max_pu; digits = 4), "] pu around x_base ", sssc_ctrl.x_base_pu, ", injected voltage ", round(abs(sssc_ctrl.x_pu - sssc_ctrl.x_base_pu) * sssc_ctrl.i_pu; digits = 4), " pu of 0.01 available")

# Reading aid: the TCSC reaches the 35-MW target (x moves from 0.20 down
# to about 0.07 pu), the SSSC pins at its injectable voltage with the flow
# stuck near 28.6 MW: its whole window is $\pm 0.033$ pu at this loading.
# Both devices report the miss honestly as `at_limit` instead of
# pretending convergence.
#
# IN range, the converter devices are ordinary voltage/flow controllers:
# with a rating that never binds, the STATCOM settles into the same
# deadband as the classical constant-Q controller. All registered devices
# describe themselves through one generic view, `controllableElements`:
# element, device, actuator with its LIVE range (for the STATCOM that is
# $\pm V \cdot S_{max}$ of the last solved operating point, not the
# nameplate), target, and status:

ok_net = build_sag_corridor("tour_facts_inrange"; with_machine = true)
addMachineVoltageControl!(ok_net; bus = "Mid", target_bus = "Load", target_vm_pu = 0.95, s_max_mva = 200.0)
run_control!(ok_net)
for e in controllableElements(ok_net)
  println(e.device, " on ", e.element, ": ", e.actuator, " in [", round(e.actuator_min; digits = 2), ", ", round(e.actuator_max; digits = 2), "], status = ", e.status)
end
println("in range: Vm(Load) = ", round(get_bus_vm_pu(ok_net, "Load"); digits = 4), " pu (target 0.95)")

# The same devices are available declaratively, as named entries under
# `control.controllers` in a configuration file (types `machine_voltage`
# with `s_max_mva`, `series_reactance` with `v_inj_max_pu`); see the
# [Control Framework](https://welthulk.github.io/Sparlectra.jl/control_framework/)
# page. The full taxonomy (incl. why a UPFC is
# deliberately NOT implemented and what to use instead) is on the
# [FACTS Devices](https://welthulk.github.io/Sparlectra.jl/facts/) page;
# `examples/others/exp_facts_limit_modes.jl` runs these contrasts as a
# script, and the TCSC has its own
# [notebook](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_series_compensation.ipynb).
#
# ## Chapter 5: N-1 contingency analysis
#
# Would the network survive the loss of any single branch? An N-1 batch
# answers that by outaging every branch in turn and re-solving. The API
# does the bookkeeping that is easy to get wrong: every case works on its
# own copy of the SOLVED base case (warm start, the base net is never
# mutated), and failures are REPORTED per case, never thrown.
#
# The 7-bus ring is too forgiving for a demo (every single outage
# leaves a connected path), so we hang a two-bus spur B8-B9 on a single
# line from B4. Losing that line strands the PAIR as an island with a live
# branch but no voltage reference, and the report has to say so instead of
# crashing. (A single stranded bus would not do: `markIsolatedBuses!`
# inside the batch marks buses without any in-service branch as isolated
# and solves the rest, so a one-bus stub converges cleanly.)

net_n1 = build_ring7("tour_n1")
addBus!(net = net_n1, busName = "B8", vn_kV = 110.0)
addBus!(net = net_n1, busName = "B9", vn_kV = 110.0)
addPIModelACLine!(net = net_n1, fromBus = "B4", toBus = "B8", r_pu = 0.02, x_pu = 0.10, b_pu = 0.0, status = 1)
addPIModelACLine!(net = net_n1, fromBus = "B8", toBus = "B9", r_pu = 0.02, x_pu = 0.10, b_pu = 0.0, status = 1)
addProsumer!(net = net_n1, busName = "B8", type = "LOAD", p = 6.0, q = 2.0)
addProsumer!(net = net_n1, busName = "B9", type = "LOAD", p = 4.0, q = 1.0)
ok_n1, msg_n1 = validate!(net = net_n1)
ok_n1 || error("N-1 net invalid: $msg_n1")

cases = generateN1Branches(net_n1)
println(length(cases), " single-branch outage cases")
results = runContingencies!(net_n1, cases; vm_min_pu = 0.95, vm_max_pu = 1.05)
printContingencyResults(results)

# Reading aid, three things to check in the table:
#
# - the `B_ACL_110_4_8` outage (the B4-B8 spur line) reports
#   `islanded without reference`: the B8-B9 pair loses its only feed and
#   the case is a reported failure, not a crash;
# - the ring outages converge, and the tight 0.95-pu floor flags the
#   voltage sags an outage causes (the `V-viol` column counts the buses
#   outside the band per case; the B1-B2 outage sags the ring to 0.86 pu);
# - iterations stay small: every case starts from the solved base-case
#   voltages (warm start), not from flat.
#
# The result rows are plain data, so post-processing is ordinary Julia:
# rank the converged cases by their voltage sag and export everything as
# CSV for a report or a spreadsheet:

ranked = sort([r for r in results if r.converged]; by = r -> r.min_vm_pu)
println("worst three outages by post-outage Vmin:")
for r in ranked[1:3]
  println("  ", rpad(r.name, 16), " Vmin = ", round(r.min_vm_pu; digits = 4), " pu, ", length(r.voltage_violations), " bus(es) below the band")
end
csv_n1 = joinpath(mktempdir(), "n1_results.csv")
writeContingencyResultsCSV(csv_n1, results)
println("exported: ", csv_n1, " (", filesize(csv_n1), " bytes)")

# Two knobs worth knowing: `retry_flat_start = true` grants one bounded
# flat-start retry per failed case (the rescue ladder is deliberately NOT
# in this loop), and `generateN1Branches(net; include_transformers = false)`
# restricts the case list to AC lines. Parallel circuits are disambiguated
# as `name#branchIdx`, and imported MATPOWER FOR001 outage lists work as
# case sources too; see
# [N-1 Contingency Analysis](https://welthulk.github.io/Sparlectra.jl/contingency/).
# On a real case this batch is exactly what the threads in the next
# chapter are for.
#
# ## Part II: Beyond
#
# ## Chapter 6: using your cores
#
# Independent work items (island solves, short-circuit fault sweeps, N-1
# contingency batches) fan out over Julia THREADS since 0.9.10, gated by
# `runtime.parallel.*` (default on). Threads are fixed at process start:
# `julia --threads=auto` uses all cores, and without the flag everything
# runs serially through the identical code path. First: how many do we
# have right now?

println("this session runs on ", Threads.nthreads(), " Julia thread(s)")

# The mechanics on a fault sweep: eight feeder-fed rings, every bus a
# fault location, once serial and once parallel. The results must be
# IDENTICAL: parallelism only changes the wall clock, never a number.

net_par = Net(name = "tour_parallel", baseMVA = 100.0)
for k in 1:8
  bus = i -> "P$(k)_B$(i)"
  for i in 1:500
    addBus!(net = net_par, busName = bus(i), vn_kV = 110.0)
  end
  addExternalGrid!(net = net_par, busName = bus(1), vm_pu = 1.0, sk_max_MVA = 2000.0 + 100.0 * k, sk_min_MVA = 1500.0, rx_max = 0.1, internal_impedance = false)
  for i in 1:500
    addPIModelACLine!(net = net_par, fromBus = bus(i), toBus = bus(i == 500 ? 1 : i + 1), r_pu = 0.001, x_pu = 0.004, b_pu = 0.0, status = 1)
  end
end
validate!(net = net_par)

runShortCircuit!(net_par; case = :max, parallel_enabled = false)      ## warm both paths
runShortCircuit!(net_par; case = :max, parallel_min_work_items = 2)
t_ser = @elapsed sc_ser = runShortCircuit!(net_par; case = :max, parallel_enabled = false)
t_par = @elapsed sc_par = runShortCircuit!(net_par; case = :max, parallel_min_work_items = 2)
println("fault sweep over ", length(sc_ser.rows), " buses: serial ", round(t_ser * 1000; digits = 1), " ms, parallel ", round(t_par * 1000; digits = 1), " ms (", round(t_ser / t_par; digits = 2), "x)")
println("rows identical: ", isequal(sc_ser.rows, sc_par.rows))

# Reading aid: on Colab's free tier you usually get 1-2 vCPUs, so the
# factor here stays modest (with one thread the parallel call falls back
# to the very same serial function). The real effect wants your local
# machine: `julia --threads=auto --project=.
# examples/run_parallel_suite.jl` runs the three dedicated demos, island
# solving (`power_flow.islands.mode: solve_parallel`), this fault sweep at
# 8000 buses, and a full N-1 contingency batch on case1354pegase (measured
# 71.7 s serial vs 17.6 s on 16 threads), each asserting serial/parallel
# identity. The N-1 batch itself is chapter 5; on a large case it is the
# prime customer of these threads.


# ## Where to go next
#
# The basic tour (Newcomer to Advanced) and the focused notebooks:
#
# - [Workshop tour, basic](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_tour.ipynb)
# - [Slack types and short circuit](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_slack_short_circuit.ipynb)
# - [Distributed slack](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_distributed_slack.ipynb)
# - [Transformer taps](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_transformers.ipynb)
# - [TCSC flow steering](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_series_compensation.ipynb)
# - [State estimation](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_state_estimation.ipynb)
#
# And the documentation for going further:
#
# - [FACTS Devices](https://welthulk.github.io/Sparlectra.jl/facts/)
# - [Control Framework](https://welthulk.github.io/Sparlectra.jl/control_framework/)
# - [Remote Voltage Control](https://welthulk.github.io/Sparlectra.jl/remote_voltage_control/)
# - [HVDC Back-to-Back](https://welthulk.github.io/Sparlectra.jl/hvdc_back_to_back/)
# - [State Estimation](https://welthulk.github.io/Sparlectra.jl/state_estimation/)
# - [N-1 Contingency Analysis](https://welthulk.github.io/Sparlectra.jl/contingency/)
# - [Feature Matrix](https://welthulk.github.io/Sparlectra.jl/feature_matrix/)
