# Converted from FOR001.DAT fixed-column dataset.
# Conversion basis: 100 MVA.
# R/X are converted from Ohm using the from-side bus base kV.
# Branch shunt B is converted from Siemens; branch shunt G is preserved as split bus shunt conductance.
using Sparlectra

const FOR001_BASE_MVA = 100.0

function build_for001_testnetz13(; bus_shunt_model = "admittance")
  net = Net(name = "FOR001_Testnetz13", baseMVA = FOR001_BASE_MVA, bus_shunt_model = bus_shunt_model)

  # Buses. PF bus types are intentionally not set here; Sparlectra derives them from attached prosumers.
  addBus!(net = net, busName = "DELTA2S1", vn_kV = 230.0, vm_pu = 1.0, va_deg = 0.0)
  addBus!(net = net, busName = "BETA2 S1", vn_kV = 230.0, vm_pu = 1.0, va_deg = 0.0)
  addBus!(net = net, busName = "ASTADTS1", vn_kV = 230.0, vm_pu = 1.0, va_deg = 0.0)
  addBus!(net = net, busName = "NORD  S1", vn_kV = 230.0, vm_pu = 1.0, va_deg = 0.0)
  addBus!(net = net, busName = "OST   S1", vn_kV = 230.0, vm_pu = 1.0, va_deg = 0.0)
  addBus!(net = net, busName = "BSTADTS1", vn_kV = 230.0, vm_pu = 1.0, va_deg = 0.0)
  addBus!(net = net, busName = "WEILERS1", vn_kV = 230.0, vm_pu = 1.0, va_deg = 0.0)
  addBus!(net = net, busName = "BURG  S1", vn_kV = 230.0, vm_pu = 1.0, va_deg = 0.0)
  addBus!(net = net, busName = "GAMMA S1", vn_kV = 230.0, vm_pu = 1.0, va_deg = 0.0)
  addBus!(net = net, busName = "SUED  S1", vn_kV = 230.0, vm_pu = 1.0, va_deg = 0.0)
  addBus!(net = net, busName = "ALPHA S1", vn_kV = 400.0, vm_pu = 1.0, va_deg = 0.0)
  addBus!(net = net, busName = "BETA1 S1", vn_kV = 400.0, vm_pu = 1.0, va_deg = 0.0)
  addBus!(net = net, busName = "DELTA1S1", vn_kV = 400.0, vm_pu = 1.0, va_deg = 0.0)

  # Branches
  # L1- ALPHA S1 -> DELTA1S1
  addPIModelACLine!(net = net, fromBus = "ALPHA S1", toBus = "DELTA1S1", r_pu = 0.00063125, x_pu = 0.00606875, b_pu = 2.368, status = 1, ratedS = 2639.645431)
  # L1- ALPHA S1 -> BETA1 S1
  addPIModelACLine!(net = net, fromBus = "ALPHA S1", toBus = "BETA1 S1", r_pu = 0.000376875, x_pu = 0.0033125, b_pu = 1.288, status = 1, ratedS = 2639.645431)
  # L1B BETA1 S1 -> DELTA1S1
  addPIModelACLine!(net = net, fromBus = "BETA1 S1", toBus = "DELTA1S1", r_pu = 0.000725, x_pu = 0.00625, b_pu = 2.432, status = 1, ratedS = 2639.645431)
  # L1A BETA1 S1 -> DELTA1S1
  addPIModelACLine!(net = net, fromBus = "BETA1 S1", toBus = "DELTA1S1", r_pu = 0.0007125, x_pu = 0.0061375, b_pu = 2.4, status = 1, ratedS = 2639.645431)
  # L1C BETA2 S1 -> BURG  S1
  addPIModelACLine!(net = net, fromBus = "BETA2 S1", toBus = "BURG  S1", r_pu = 0.00491493383743, x_pu = 0.0206049149338, b_pu = 0.0391989, status = 1, ratedS = 288.8194722)
  # L1A BETA2 S1 -> BURG  S1
  addPIModelACLine!(net = net, fromBus = "BETA2 S1", toBus = "BURG  S1", r_pu = 0.0047258979206, x_pu = 0.0206049149338, b_pu = 0.0394105, status = 1, ratedS = 288.8194722)
  # L1B BETA2 S1 -> BURG  S1
  addPIModelACLine!(net = net, fromBus = "BETA2 S1", toBus = "BURG  S1", r_pu = 0.00491493383743, x_pu = 0.0206049149338, b_pu = 0.0391989, status = 1, ratedS = 288.8194722)
  # L1B BURG  S1 -> GAMMA S1
  addPIModelACLine!(net = net, fromBus = "BURG  S1", toBus = "GAMMA S1", r_pu = 0.00232514177694, x_pu = 0.00984877126654, b_pu = 0.020102, status = 1, ratedS = 557.72036)
  # L1A BURG  S1 -> GAMMA S1
  addPIModelACLine!(net = net, fromBus = "BURG  S1", toBus = "GAMMA S1", r_pu = 0.00232514177694, x_pu = 0.00984877126654, b_pu = 0.020102, status = 1, ratedS = 557.72036)
  # L1- BSTADTS1 -> OST   S1
  addPIModelACLine!(net = net, fromBus = "BSTADTS1", toBus = "OST   S1", r_pu = 0.00232514177694, x_pu = 0.00984877126654, b_pu = 0.020102, status = 1, ratedS = 557.72036)
  # L1- ASTADTS1 -> OST   S1
  addPIModelACLine!(net = net, fromBus = "ASTADTS1", toBus = "OST   S1", r_pu = 0.00232514177694, x_pu = 0.00984877126654, b_pu = 0.020102, status = 1, ratedS = 557.72036)
  # L1- DELTA2S1 -> ASTADTS1
  addPIModelACLine!(net = net, fromBus = "DELTA2S1", toBus = "ASTADTS1", r_pu = 0.0153119092628, x_pu = 0.065595463138, b_pu = 0.141243, status = 1, ratedS = 557.72036)
  # L1- WEILERS1 -> ASTADTS1
  addPIModelACLine!(net = net, fromBus = "WEILERS1", toBus = "ASTADTS1", r_pu = 0.00113421550095, x_pu = 0.00555765595463, b_pu = 0.0184621, status = 1, ratedS = 557.72036)
  # L1- DELTA2S1 -> WEILERS1
  addPIModelACLine!(net = net, fromBus = "DELTA2S1", toBus = "WEILERS1", r_pu = 0.00241965973535, x_pu = 0.01202268431, b_pu = 0.0398866, status = 1, ratedS = 557.72036)
  # L1- ASTADTS1 -> NORD  S1
  addPIModelACLine!(net = net, fromBus = "ASTADTS1", toBus = "NORD  S1", r_pu = 0.00241965973535, x_pu = 0.01202268431, b_pu = 0.0398866, status = 1, ratedS = 557.72036)
  # L1B ASTADTS1 -> BSTADTS1
  addPIModelACLine!(net = net, fromBus = "ASTADTS1", toBus = "BSTADTS1", r_pu = 0.00253308128544, x_pu = 0.0127032136106, b_pu = 0.0062422, status = 1, ratedS = 278.86018)
  # L1A ASTADTS1 -> BSTADTS1
  addPIModelACLine!(net = net, fromBus = "ASTADTS1", toBus = "BSTADTS1", r_pu = 0.00253308128544, x_pu = 0.0127032136106, b_pu = 0.0062422, status = 1, ratedS = 278.86018)
  # L1- BSTADTS1 -> WEILERS1
  addPIModelACLine!(net = net, fromBus = "BSTADTS1", toBus = "WEILERS1", r_pu = 0.0103969754253, x_pu = 0.0438563327032, b_pu = 0.08464, status = 1, ratedS = 557.72036)
  # L1- BSTADTS1 -> BURG  S1
  addPIModelACLine!(net = net, fromBus = "BSTADTS1", toBus = "BURG  S1", r_pu = 0.00132325141777, x_pu = 0.00879017013233, b_pu = 0.029095, status = 1, ratedS = 557.72036)
  # L1- BETA2 S1 -> WEILERS1
  addPIModelACLine!(net = net, fromBus = "BETA2 S1", toBus = "WEILERS1", r_pu = 0.00132325141777, x_pu = 0.00879017013233, b_pu = 0.029095, status = 1, ratedS = 557.72036)
  # L1B BSTADTS1 -> SUED  S1
  addPIModelACLine!(net = net, fromBus = "BSTADTS1", toBus = "SUED  S1", r_pu = 0.00232514177694, x_pu = 0.00984877126654, b_pu = 0.020102, status = 1, ratedS = 557.72036)
  # L1A BSTADTS1 -> SUED  S1
  addPIModelACLine!(net = net, fromBus = "BSTADTS1", toBus = "SUED  S1", r_pu = 0.00232514177694, x_pu = 0.00984877126654, b_pu = 0.020102, status = 1, ratedS = 557.72036)
  # T1C BETA1 S1 -> BETA2 S1
  addPIModelTrafo!(net = net, fromBus = "BETA1 S1", toBus = "BETA2 S1", r_pu = 0.0001125, x_pu = 0.0051375, b_pu = -0.02944, status = 1, ratedU = 400.0, ratedS = 1143.153533, ratio = 0.995670995671, shift_deg = 0.0)
  # T1A BETA1 S1 -> BETA2 S1
  addPIModelTrafo!(net = net, fromBus = "BETA1 S1", toBus = "BETA2 S1", r_pu = 0.0001125, x_pu = 0.0051375, b_pu = -0.02944, status = 1, ratedU = 400.0, ratedS = 1143.153533, ratio = 0.995670995671, shift_deg = 0.0)
  # T1B BETA1 S1 -> BETA2 S1
  addPIModelTrafo!(net = net, fromBus = "BETA1 S1", toBus = "BETA2 S1", r_pu = 0.0001125, x_pu = 0.0051375, b_pu = -0.02944, status = 1, ratedU = 400.0, ratedS = 1143.153533, ratio = 0.995670995671, shift_deg = 0.0)
  # T1B DELTA1S1 -> DELTA2S1
  addPIModelTrafo!(net = net, fromBus = "DELTA1S1", toBus = "DELTA2S1", r_pu = 6.5e-05, x_pu = 0.0051, b_pu = -0.01616, status = 1, ratedU = 400.0, ratedS = 1731.357987, ratio = 0.995670995671, shift_deg = 0.0)
  # T1A DELTA1S1 -> DELTA2S1
  addPIModelTrafo!(net = net, fromBus = "DELTA1S1", toBus = "DELTA2S1", r_pu = 6.5e-05, x_pu = 0.0051, b_pu = -0.01616, status = 1, ratedU = 400.0, ratedS = 1731.357987, ratio = 0.995670995671, shift_deg = 0.0)

  # Split branch shunt conductance that standard MATPOWER/Sparlectra PI calls cannot attach as branch G.
  addShuntMatpower!(net = net, busName = "DELTA2S1", Gs = 0.8096, Bs = 0.0)
  addShuntMatpower!(net = net, busName = "BETA2 S1", Gs = 1.4856, Bs = 0.0)
  addShuntMatpower!(net = net, busName = "BETA1 S1", Gs = 1.4856, Bs = 0.0)
  addShuntMatpower!(net = net, busName = "DELTA1S1", Gs = 0.8096, Bs = 0.0)

  # Loads
  addProsumer!(net = net, busName = "ASTADTS1", type = "ENERGYCONSUMER", p = 300.0, q = 100.0)
  addProsumer!(net = net, busName = "OST   S1", type = "ENERGYCONSUMER", p = 200.0, q = 50.0)
  addProsumer!(net = net, busName = "WEILERS1", type = "ENERGYCONSUMER", p = 150.0, q = 30.0)
  addProsumer!(net = net, busName = "BURG  S1", type = "ENERGYCONSUMER", p = 200.0, q = 50.0)
  addProsumer!(net = net, busName = "GAMMA S1", type = "ENERGYCONSUMER", p = 300.0, q = 100.0)
  addProsumer!(net = net, busName = "ALPHA S1", type = "ENERGYCONSUMER", p = 400.0, q = 200.0)

  # Generators. Only the REF bus is regulating here; the other generators are fixed PQ injections because no PV/Q-limit data is present.
  addProsumer!(net = net, busName = "BETA2 S1", type = "GENERATOR", p = 600.0, q = 200.0)
  addProsumer!(net = net, busName = "NORD  S1", type = "GENERATOR", p = 100.0, q = 40.0)
  addProsumer!(net = net, busName = "BSTADTS1", type = "SYNCHRONOUSMACHINE", p = 300.0, q = 100.0, referencePri = "BSTADTS1", vm_pu = 1.0, va_deg = 0.0)
  addProsumer!(net = net, busName = "WEILERS1", type = "GENERATOR", p = 400.0, q = 200.0)
  addProsumer!(net = net, busName = "SUED  S1", type = "GENERATOR", p = 150.0, q = 50.0)

  valid, msg = validate!(net = net)
  valid || error(msg)
  return net
end

# Example:
# net = build_for001_testnetz13()
# ite, erg = runpf!(net, 10, 5e-3, 0)
