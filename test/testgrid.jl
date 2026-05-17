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

# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 20.6.2023
# CIGRE HV network
# Source of this network can be found here (Task Force C6.04.02 ): https://www.researchgate.net/publication/271963972_TF_C60402_TB_575_--_Benchmark_Systems_for_Network_Integration_of_Renewable_and_Distributed_Energy_Resources

function test_2WTPITrafo()
  Sbase_MVA = 1000.0
  netName = "trafo_2W_PIT"
  net = Net(name = netName, baseMVA = Sbase_MVA)
  @debug "Creating $netName test network"
  addBus!(net = net, busName = "B1", vn_kV = 220.0)
  addBus!(net = net, busName = "B2", vn_kV = 20.0)

  add2WTPIModelTrafo!(net = net, fromBus = "B1", toBus = "B2", r = 0.0, x = 0.4, ratedU = 220.0, ratedS = 1000.0)
  addProsumer!(net = net, busName = "B2", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "B1")
  result, msg = validate!(net = net)
  if !result
    @warn msg
    return false
  end
  return result
end

function test_3WTPITrafo(verbose::Int = 0; method::Symbol = :rectangular, opt_fd::Bool = false, opt_sparse::Bool = true)::Bool
  #
  #                                                  380kV      AUX*        110kV
  # 380kV              -----------------------AC-LINE--o----------o----------o----------------AC-LINE---------o  Load: p=80 MW, q=30 MVar
  #                 Slack                                    |
  #                                                          | 20kV
  #                                                          o  Shunt: p=0 MW, q=5 MVar
  #

  Sbase_MVA = 1000.0
  netName = "trafo_3W_PIT"
  net = Net(name = netName, baseMVA = Sbase_MVA)
  @debug "Creating $netName test network"

  # --- buses ---
  addBus!(net = net, busName = "B1", vn_kV = 380.0)
  addBus!(net = net, busName = "B2", vn_kV = 380.0)  # HV side of 3W trafo
  addBus!(net = net, busName = "B3", vn_kV = 110.0)  # MV side
  addBus!(net = net, busName = "B4", vn_kV = 20.0)   # LV side (shunt bus)
  addBus!(net = net, busName = "B5", vn_kV = 110.0)  # load bus at 110kV

  # --- slack injection ---
  addProsumer!(net = net, busName = "B1", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "B1")

  # --- line 380 kV: B1 -- B2 ---
  addACLine!(net = net, fromBus = "B1", toBus = "B2", length = 1.0, r = 0.01, x = 0.10)

  # --- 3-winding transformer: (B2=380kV, B3=110kV, B4=20kV) ---
  # Use simple physical values; add3WTPiModelTrafo! is expected to:
  #   - create AUX bus (if missing)
  #   - convert to PU internally via add2WTPIModelTrafo!
  aux_bus = add3WTPiModelTrafo!(
    net = net,
    HBBus = "B2",
    MBBus = "B3",
    LVBus = "B4",
    r = (0.20, 0.30, 0.40),                 # Ohm (example)
    x = (4.00, 6.00, 10.00),                # Ohm (example)
    b = (0.0, 0.0, 0.0),                    # Siemens (or use 0.0 if not modeled)
    ratedU_kV = (380.0, 110.0, 20.0),       # kV
    ratedS_MVA = (1000.0, 500.0, 200.0),    # MVA (example)
    status = 1,
  )

  # --- line 110 kV: B3 -- B5 (to remote load bus) ---
  addACLine!(net = net, fromBus = "B3", toBus = "B5", length = 1.0, r = 0.01, x = 0.10)

  # --- load at B5: 80 MW / 30 MVar ---
  addProsumer!(net = net, busName = "B5", type = "ENERGYCONSUMER", p = 80.0, q = 30.0)

  # --- shunt at B4: q = 5 MVar (capacitive/inductive sign depends on your convention) ---
  addShunt!(net = net, busName = "B4", pShunt = 0.0, qShunt = 5.0, in_service = 1)

  # --- structural checks ---
  if !hasBusInNet(net = net, busName = aux_bus)
    @warn "AUX bus was not created/found by add3WTPiModelTrafo!" aux_bus = aux_bus
    return false
  end

  # Expect at least:
  #  - 2 AC lines  (B1-B2, B3-B5)
  #  - 3 trafo branches (AUX-B2, AUX-B3, AUX-B4)
  expected_min_branches = 5
  if length(net.branchVec) < expected_min_branches
    @warn "Unexpectedly low branch count" n = length(net.branchVec) expected_min = expected_min_branches
    return false
  end

  # --- final validate ---
  result, msg = validate!(net = net)
  if !result
    @warn msg
    return false
  end

  tol = 1e-9
  maxIte = 50
  print_results = (verbose > 0)
  result = true

  etim = 0.0
  etim = @elapsed begin
    ite, erg = runpf!(net, maxIte, tol, verbose, method = method, opt_fd = opt_fd, opt_sparse = opt_sparse)
    if erg != 0
      @info "Full-system power flow did not converge"
      result = false
    end
  end

  V = buildVoltageVector(net)
  calcNetLosses!(net, V)
  distributeBusResults!(net)
  if print_results
    printACPFlowResults(net, etim, ite, tol; converged = result, solver = method)
    printProsumerResults(net)
    printQLimitLog(net; sort_by = :bus)
  end

  return result
end

function test_acpflow(verbose::Int = 0; lLine_6a6b::Float64 = 0.01, damp::Float64 = 1.0, method::Symbol = :rectangular, opt_sparse = true)::Bool
  net = createCIGRE(lLine_6a6b)
  tol = 1e-6
  maxIte = 25
  print_results = (verbose > 0)
  result = true
  etime = @elapsed begin
    ite, erg = runpf!(net, maxIte, tol, verbose; method = method, damp = damp, opt_sparse = opt_sparse)
  end
  if erg != 0
    @warn "Power flow did not converge"
    result = false
  end

  if print_results
    V = buildVoltageVector(net)
    calcNetLosses!(net, V)
    distributeBusResults!(net)
    printACPFlowResults(net, etime, ite, tol; converged = (erg == 0), solver = method)
    printProsumerResults(net)
    printQLimitLog(net; sort_by = :bus)
  end

  return result
end

function test_NBI_MDO()
  myNet = createCIGRE()
  result = true

  nodeNumberVec = Vector{Int}()
  branchTupleSet = Set{Tuple}()

  for n in myNet.nodeVec
    if n.busIdx in nodeNumberVec
      @warn "node", n, "is already in set!"
    else
      push!(nodeNumberVec, n.busIdx)
    end
  end

  for b in myNet.branchVec
    tuple = (b.fromBus, b.toBus)
    if tuple ∉ branchTupleSet
      push!(branchTupleSet, tuple)
    end
  end

  A = getNBI(nodeNumberVec, branchTupleSet)
  order = mdoRCM(length(nodeNumberVec), branchTupleSet)
  if order == [2, 5, 8, 7, 13, 3, 12, 11, 1, 4, 6, 9, 10]
    result = result & true
  else
    @warn "expected MDO-Order: [2, 5, 8, 7, 13, 3, 12, 11, 1, 4, 6, 9, 10], result order: $(order)"
    result = false
  end

  return result
end

function testNetwork()::Bool
  myNet = createCIGRE()
  result, msg = validate!(net = myNet)
  return result
end

function getTestFilePathName()
  filename = "cigre.m"
  jpath = joinpath(pwd(), "data", "mpower", filename)
  return jpath
end

function testExportMatpower()
  myNet = createCIGRE()
  case = myNet.name
  writeMatpowerCasefile(myNet, getTestFilePathName())
  return true
end

function testImportMatpower()::Bool
  mpc = (
    name = "case2_inline",
    baseMVA = 100.0,

    # bus: [bus_i type Pd Qd Gs Bs area Vm Va baseKV zone Vmax Vmin]
    bus = [
      1 3 0.0 0.0 0.0 0.0 1 1.0 0.0 110.0 1 1.1 0.9
      2 1 50.0 30.0 0.0 0.0 1 1.0 0.0 110.0 1 1.1 0.9
    ],

    # gen: 21 cols (MATPOWER v2)
    gen = [
      1 0.0 0.0 999.0 -999.0 1.0 100.0 1 999.0 0.0 0 0 0 0 0 0 0 0 0 0 0;
    ],

    # branch: [fbus tbus r x b rateA rateB rateC ratio angle status angmin angmax]
    branch = [
      1 2 0.01 0.05 0.0 9999.0 0.0 0.0 0.0 0.0 1 -60.0 60.0;
    ],
    gencost = nothing,
    bus_name = nothing,
  )

  net = Sparlectra.createNetFromMatPowerCase(mpc = mpc, log = false, flatstart = false)

  if strip(net.name) != "case2_inline"
    @warn "Unexpected net.name" got = net.name
    return false
  end
  if length(net.nodeVec) != 2
    @warn "Expected 2 nodes" got = length(net.nodeVec)
    return false
  end
  if length(net.branchVec) != 1
    @warn "Expected 1 branch" got = length(net.branchVec)
    return false
  end

  ok, msg = validate!(net = net)
  ok || (@warn msg; return false)

  return true
end

function test_matpower_import_defaults_no_reenable()::Bool
  mpc = (
    name = "case2_defaults",
    baseMVA = 100.0,
    bus = [
      1 3 0.0 0.0 0.0 0.0 1 1.0 0.0 110.0 1 1.1 0.9
      2 1 20.0 10.0 0.0 0.0 1 1.0 0.0 110.0 1 1.1 0.9
    ],
    gen = [
      1 30.0 0.0 999.0 -999.0 1.0 100.0 1 999.0 0.0 0 0 0 0 0 0 0 0 0 0 0
    ],
    branch = [
      1 2 0.01 0.05 0.0 9999.0 0.0 0.0 0.0 0.0 1 -60.0 60.0
    ],
    gencost = nothing,
    bus_name = nothing,
  )

  net = Sparlectra.createNetFromMatPowerCase(mpc = mpc, log = false, flatstart = false)

  if net.cooldown_iters != 0
    @warn "Expected default cooldown_iters = 0 for MATPOWER import" got = net.cooldown_iters
    return false
  end

  if !isapprox(net.q_hyst_pu, 0.0; atol = 0.0, rtol = 0.0)
    @warn "Expected default q_hyst_pu = 0.0 for MATPOWER import" got = net.q_hyst_pu
    return false
  end

  return true
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

  return slack_idx == 1 &&
         isapprox(abs(V0[1]), 1.04; atol = 1e-12) && isapprox(angle(V0[1]), deg2rad(12.0); atol = 1e-12) &&
         isapprox(abs(V0[2]), 1.0; atol = 1e-12) && isapprox(angle(V0[2]), 0.0; atol = 1e-12) &&
         isapprox(abs(V0[3]), 1.02; atol = 1e-12) && isapprox(angle(V0[3]), 0.0; atol = 1e-12) &&
         isapprox(busVec[1].vm_pu, 1.04; atol = 1e-12) && isapprox(busVec[1].va_rad, deg2rad(12.0); atol = 1e-12) &&
         isapprox(busVec[2].vm_pu, 1.0; atol = 1e-12) && isapprox(busVec[2].va_rad, 0.0; atol = 1e-12) &&
         isapprox(busVec[3].vm_pu, 1.02; atol = 1e-12) && isapprox(busVec[3].va_rad, 0.0; atol = 1e-12)
end

function test_matpower_reference_override_controls_flatstart()::Bool
  mpc = (
    name = "case_reference_override",
    baseMVA = 100.0,
    bus = [
      1 3 0.0 0.0 0.0 0.0 1 0.97 12.0 110.0 1 1.1 0.9
      2 1 20.0 8.0 0.0 0.0 1 0.93 -7.0 110.0 1 1.1 0.9
    ],
    gen = [
      1 30.0 0.0 999.0 -999.0 1.04 100.0 1 999.0 0.0 0 0 0 0 0 0 0 0 0 0 0
    ],
    branch = [
      1 2 0.01 0.05 0.0 9999.0 0.0 0.0 0.0 0.0 1 -60.0 60.0
    ],
    gencost = nothing,
    bus_name = nothing,
  )

  net = Sparlectra.createNetFromMatPowerCase(mpc = mpc, log = false, flatstart = true, reference_vm_pu = 1.01, reference_va_deg = -3.0)
  V0, slack_idx = Sparlectra.initialVrect(net; flatstart = true)
  busVec, _ = Sparlectra.getBusData(net.nodeVec, net.baseMVA, true; net = net)

  return slack_idx == 1 &&
         isapprox(abs(V0[1]), 1.01; atol = 1e-12) && isapprox(angle(V0[1]), deg2rad(-3.0); atol = 1e-12) &&
         isapprox(busVec[1].vm_pu, 1.01; atol = 1e-12) && isapprox(busVec[1].va_rad, deg2rad(-3.0); atol = 1e-12)
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

  return getNodeType(net.nodeVec[b2]) == Sparlectra.PQ &&
         getNodeType(net.nodeVec[b3]) == Sparlectra.PV &&
         !isnothing(b2_gen.puController) &&
         !isnothing(b2_gen.quController) &&
         isnothing(b3_gen.puController) &&
         isnothing(b3_gen.quController)
end

function test_prosumer_aggregation_preserves_bus_types_and_injections()::Bool
  net = Net(name = "aggregation_regression", baseMVA = 100.0)
  addBus!(net = net, busName = "Slack", vn_kV = 110.0)
  addBus!(net = net, busName = "PV", vn_kV = 110.0)
  addBus!(net = net, busName = "PQ", vn_kV = 110.0)

  addACLine!(net = net, fromBus = "Slack", toBus = "PV", length = 1.0, r = 0.01, x = 0.10)
  addACLine!(net = net, fromBus = "PV", toBus = "PQ", length = 1.0, r = 0.01, x = 0.10)

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

  return getNodeType(net.nodeVec[slack]) == Sparlectra.Slack &&
         getNodeType(net.nodeVec[pv]) == Sparlectra.PV &&
         getNodeType(net.nodeVec[pq]) == Sparlectra.PQ &&
         isapprox(S[pv], ComplexF64(0.40, 0.02); atol = 1e-12) &&
         isapprox(S[pq], ComplexF64(-0.18, -0.06); atol = 1e-12)
end

function test_rectangular_autodamp_backtracks_oversized_step()::Bool
  Y = ComplexF64[0.0 - 10.0im 0.0 + 10.0im; 0.0 + 10.0im 0.0 - 10.0im]
  V = ComplexF64[1.0 + 0.0im, 1.0 + 0.0im]
  S = ComplexF64[0.0 + 0.0im, -1.0 - 0.2im]
  bus_types = [:Slack, :PQ]
  Vset = [1.0, 1.0]
  F0 = Sparlectra.mismatch_rectangular(Y, V, S, bus_types, Vset, 1)

  oversized_delta = [0.0, -1.0]
  alpha, Vtrial, trial_mismatch = Sparlectra.choose_rectangular_autodamp(
    Y,
    V,
    S,
    oversized_delta,
    F0;
    slack_idx = 1,
    damp = 1.0,
    autodamp_min = 1e-3,
    bus_types = bus_types,
    Vset = Vset,
  )
  full_step_mismatch = maximum(abs.(Sparlectra.mismatch_rectangular(Y, ComplexF64[1.0 + 0.0im, 1.0 - 1.0im], S, bus_types, Vset, 1)))

  return alpha < 1.0 &&
         trial_mismatch < maximum(abs.(F0)) &&
         trial_mismatch < full_step_mismatch &&
         Vtrial[1] == V[1]
end

function test_rectangular_start_projection_improves_dc_seed()::Bool
  Ydense = ComplexF64[0.0 - 10.0im 0.0 + 10.0im; 0.0 + 10.0im 0.0 - 10.0im]
  Y = sparse(Ydense)
  Vraw = ComplexF64[1.0 + 0.0im, 1.0 + 0.0im]
  S = ComplexF64[0.0 + 0.0im, -1.0 - 0.2im]
  bus_types = [:Slack, :PQ]
  Vset = [1.0, 1.0]
  profile = Dict{Symbol,Any}(:enabled => true)

  raw_mismatch = maximum(abs.(Sparlectra.mismatch_rectangular(Y, Vraw, S, bus_types, Vset, 1)))
  Vproj = Sparlectra.project_rectangular_start(
    Y,
    Vraw,
    S,
    bus_types,
    Vset,
    1;
    enabled = true,
    try_dc_start = true,
    try_blend_scan = true,
    blend_lambdas = [0.25, 0.5, 0.75],
    dc_angle_limit_deg = 60.0,
    performance_profile = profile,
  )
  projected_mismatch = maximum(abs.(Sparlectra.mismatch_rectangular(Y, Vproj, S, bus_types, Vset, 1)))

  return Vproj[1] == Vraw[1] &&
         projected_mismatch < raw_mismatch &&
         profile[:dc_matrix_size] == (1, 1) &&
         profile[:dc_matrix_nnz] == 1 &&
         profile[:dc_solve_backend] === :sparse_lu_umfpack &&
         profile[:dc_solve_reduced_dimension] == 1
end

function test_matpower_vmva_selfcheck_noncontiguous_bus_numbers()::Bool
  mpc = Sparlectra.MatpowerIO.MatpowerCase(
    "case_noncontig",
    100.0,
    # bus_i includes non-contiguous high id (9001)
    [
      1 3 0.0 0.0 0.0 0.0 1 1.00 0.0 110.0 1 1.1 0.9
      9001 1 0.0 0.0 0.0 0.0 1 0.99 -3.0 110.0 1 1.1 0.9
    ],
    [
      1 0.0 0.0 999.0 -999.0 1.0 100.0 1 999.0 0.0 0 0 0 0 0 0 0 0 0 0 0
    ],
    [
      1 9001 0.01 0.05 0.0 9999.0 0.0 0.0 0.0 0.0 1 -60.0 60.0
    ],
    nothing,
    nothing,
  )

  stats = Sparlectra.MatpowerIO.vmva_power_mismatch_stats(mpc)
  return get(stats, :ok, false) && get(stats, :n_p, 0) == 1 && get(stats, :n_q, 0) == 1 && isfinite(get(stats, :max_p_mis_pu, NaN)) && isfinite(get(stats, :max_q_mis_pu, NaN))
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
      1 500.0 300.0 999.0 -999.0 1.0 100.0 1 999.0 0.0 0 0 0 0 0 0 0 0 0 0 0
    ],
    [
      # no branch -> Scalc = 0 at all buses
      # keep out-of-service branch to satisfy matrix width expectations
      1 2 0.01 0.05 0.0 9999.0 0.0 0.0 0.0 0.0 0 -60.0 60.0
    ],
    nothing,
    nothing,
  )

  stats = Sparlectra.MatpowerIO.vmva_power_mismatch_stats(mpc)
  return get(stats, :ok, false) && isapprox(get(stats, :max_p_mis_pu, NaN), 0.0; atol = 1e-12) && isapprox(get(stats, :max_q_mis_pu, NaN), 0.0; atol = 1e-12)
end

function test_matpower_compare_vmva_wraps_angle_differences()::Bool
  mpc = Sparlectra.MatpowerIO.MatpowerCase(
    "case_angle_wrap_compare",
    100.0,
    [
      1 3 0.0 0.0 0.0 0.0 1 1.00 0.0 110.0 1 1.1 0.9
      2 1 0.0 0.0 0.0 0.0 1 0.99 179.0 110.0 1 1.1 0.9
    ],
    [
      1 0.0 0.0 999.0 -999.0 1.0 100.0 1 999.0 0.0 0 0 0 0 0 0 0 0 0 0 0
    ],
    [
      1 2 0.01 0.05 0.0 9999.0 0.0 0.0 0.0 0.0 0 -60.0 60.0
    ],
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
  return source_ok && !bus_ok && get(bus_stats, :max_dvm, 0.0) > 0.019 && imported_ok && get(imported_stats, :max_dvm, 1.0) < 0.011 && diagnostic_ok
end

function test_matpower_compare_vmva_final_pq_after_qlimit_uses_bus_vm()::Bool
  mpc = _synthetic_pv_vg_mismatch_case()
  net = Sparlectra.createNetFromMatPowerCase(mpc = mpc, matpower_pv_voltage_source = :gen_vg)
  pv_idx = geNetBusIdx(net = net, busName = "2")

  Sparlectra.setNodeType!(net.nodeVec[pv_idx], "PQ")
  net.qLimitEvents[pv_idx] = :max
  net.nodeVec[pv_idx]._vm_pu = 1.03
  net.nodeVec[pv_idx]._va_deg = 0.0

  text_path, text_io = mktemp()
  close(text_io)
  ok, stats = open(text_path, "w") do io
    redirect_stdout(io) do
      Sparlectra.MatpowerIO.compare_vm_va(net, mpc; show_diff = true, compare_voltage_reference = :hybrid, matpower_pv_voltage_source = :gen_vg, tol_vm = 1e-12, tol_va = 180.0)
    end
  end
  text = read(text_path, String)
  rm(text_path; force = true)
  counts = get(stats, :vm_ref_kind_counts, Dict{Symbol,Int}())
  kind_summary = get(stats, :vm_ref_kind_summary, Dict{Symbol,Any}())
  final_pq_summary = get(kind_summary, :final_pq_after_qlimit, (count = 0, max_dvm = NaN, max_dva = NaN))

  return ok &&
         get(stats, :compare_status, :fail) == :warn &&
         isapprox(get(stats, :max_dvm, NaN), 0.01; atol = 1e-12) &&
         isapprox(final_pq_summary.max_dvm, 0.01; atol = 1e-12) &&
         get(counts, :final_pq_after_qlimit, 0) == 1 &&
         get(counts, :pq_bus_vm, 0) == 1 &&
         get(counts, :ref_slack_imported_setpoint, 0) == 1 &&
         get(counts, :active_pv_imported_setpoint, 0) == 0 &&
         occursin("final_pq_after_qlimit", text) &&
         occursin("pq_bus_vm", text) &&
         occursin("ref_slack_imported_setpoint", text) &&
         occursin("status           : WARN", text) &&
         occursin("Voltage deviations on final_pq_after_qlimit buses are expected", text)
end

function test_matpower_multi_generator_vg_and_fallback()::Bool
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
  return first_online_ok && fallback_ok
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

  addACLine!(net = net, fromBus = "B1", toBus = "B2", length = 100.0, r = 0.0653, x = 0.398, c_nf_per_km = 9.08, tanδ = 0.0, ratedS = 250.0, status = 1)
  addACLine!(net = net, fromBus = "B1", toBus = "B3", length = 100.0, r = 0.0653, x = 0.398, c_nf_per_km = 9.08, tanδ = 0.0, ratedS = 250.0, status = 1)
  addACLine!(net = net, fromBus = "B2", toBus = "B4", length = 100.0, r = 0.0653, x = 0.398, c_nf_per_km = 9.08, tanδ = 0.0, ratedS = 250.0, status = 1)
  addACLine!(net = net, fromBus = "B3", toBus = "B4", length = 100.0, r = 0.0653, x = 0.398, c_nf_per_km = 9.08, tanδ = 0.0, ratedS = 250.0, status = 1)
  addACLine!(net = net, fromBus = "B5", toBus = "B4", length = 300.0, r = 0.0653, x = 0.398, c_nf_per_km = 9.08, tanδ = 0.0, ratedS = 250.0, status = 1)
  addACLine!(net = net, fromBus = "B5", toBus = "B4", length = 300.0, r = 0.0653, x = 0.398, c_nf_per_km = 9.08, tanδ = 0.0, ratedS = 250.0, status = 1)

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

function createCIGRE(lLine_6a6b = 0.01)::Net
  Sbase_MVA = 1000.0
  netName = "cigre"
  net = Net(name = netName, baseMVA = Sbase_MVA)
  @debug "Creating $netName test network"
  # A R E A 1
  # Bus 1, 220kV  
  addBus!(net = net, busName = "B1", vn_kV = 220.0)
  # Bus 9, 22kV  
  addBus!(net = net, busName = "B9", vn_kV = 22.0)
  # Bus 7, 22kV  
  addBus!(net = net, busName = "B7", vn_kV = 380.0)
  # Bus 2, 220kV  
  addBus!(net = net, busName = "B2", vn_kV = 220.0)
  # Bus 10, 22kV  
  addBus!(net = net, busName = "B10", vn_kV = 22.0)
  # A R E A 2
  # B6a, 220kV  
  addBus!(net = net, busName = "B6a", vn_kV = 220.0)
  # B6b, 220kV   
  addBus!(net = net, busName = "B6b", vn_kV = 220.0)
  # Bus 12, 22kV  
  addBus!(net = net, busName = "B12", vn_kV = 22.0)
  # A R E A 3
  # Bus 5, 220kV  
  addBus!(net = net, busName = "B5", vn_kV = 220.0)
  # Bus 4, 220kV  
  addBus!(net = net, busName = "B4", vn_kV = 220.0)
  # Bus 3, 220kV
  addBus!(net = net, busName = "B3", vn_kV = 220.0)
  # Bus 8, 380kV  
  addBus!(net = net, busName = "B8", vn_kV = 380.0)
  # Bus 11, 22kV  
  addBus!(net = net, busName = "B11", vn_kV = 22.0)

  ##### Shunt    
  addShunt!(net = net, busName = "B6a", pShunt = 0.0, qShunt = 180.0)
  addShunt!(net = net, busName = "B4", pShunt = 0.0, qShunt = 160.0)
  addShunt!(net = net, busName = "B5", pShunt = 0.0, qShunt = 80.0)

  ###### ACLines   
  # 220kV-Lines
  # (bus1, bus2, length_km=100)
  addACLine!(net = net, fromBus = "B1", toBus = "B2", length = 100.0, r = 0.0653, x = 0.398, c_nf_per_km = 9.08, tanδ = 0.0)
  # (bus1, bus6a, length_km=300)
  addACLine!(net = net, fromBus = "B1", toBus = "B6a", length = 300.0, r = 0.0653, x = 0.398, c_nf_per_km = 9.08, tanδ = 0.0)
  # (bus2, bus5, length_km=300)
  addACLine!(net = net, fromBus = "B2", toBus = "B5", length = 300.0, r = 0.0653, x = 0.398, c_nf_per_km = 9.08, tanδ = 0.0)
  # (bus3, bus4, length_km=100)
  addACLine!(net = net, fromBus = "B3", toBus = "B4", length = 100.0, r = 0.0653, x = 0.398, c_nf_per_km = 9.08, tanδ = 0.0)
  # (bus3, bus4, length_km=100)
  addACLine!(net = net, fromBus = "B3", toBus = "B4", length = 300.0, r = 0.0653, x = 0.398, c_nf_per_km = 9.08, tanδ = 0.0)
  # (bus4, bus5, length_km=300)
  addACLine!(net = net, fromBus = "B4", toBus = "B5", length = 300.0, r = 0.0653, x = 0.398, c_nf_per_km = 9.08, tanδ = 0.0)
  # (bus4, bus6a, length_km=300)
  addACLine!(net = net, fromBus = "B4", toBus = "B6a", length = 300.0, r = 0.0653, x = 0.398, c_nf_per_km = 9.08, tanδ = 0.0)
  #(bus6a, bus6b, length_km=length_km_6a_6b)
  addACLine!(net = net, fromBus = "B6a", toBus = "B6b", length = lLine_6a6b, r = 0.0653, x = 0.398, c_nf_per_km = 9.08, tanδ = 0.0)

  # 380kV-Lines
  # (bus7, bus8, length_km=600)
  addACLine!(net = net, fromBus = "B7", toBus = "B8", length = 600.0, r = 0.0328, x = 0.312, c_nf_per_km = 11.5, tanδ = 0.0)

  #### Trafos  
  # 'Trafo 1-7' (bus7, bus1, sn_mva=1000,  vn_hv_kv=380, vn_lv_kv=220, vkr_percent=0.0, vk_percent=13.0, pfe_kw=0, i0_percent=0,   shift_degree=0.0)
  add2WTrafo!(net = net, fromBus = "B7", toBus = "B1", sn_mva = 1000.0, vk_percent = 13.0, vkr_percent = 0.28, pfe_kw = 20.0, i0_percent = 0.06)
  #Trafo 3-8, (bus8, bus3, sn_mva=1000,  vn_hv_kv=380, vn_lv_kv=220, vkr_percent=0.0, vk_percent=13.0, pfe_kw=0, i0_percent=0, shift_degree=0.0)
  add2WTrafo!(net = net, fromBus = "B8", toBus = "B3", sn_mva = 1000.0, vk_percent = 13.0, vkr_percent = 0.0, pfe_kw = 0.0, i0_percent = 0.0)
  # Trafo 9-1 (bus1, bus9, sn_mva=1000,vn_hv_kv=220, vn_lv_kv=22, vkr_percent=0.0,vk_percent=13.0, pfe_kw=0, i0_percent=0,shift_degree=330.0)
  add2WTrafo!(net = net, fromBus = "B1", toBus = "B9", sn_mva = 1000.0, vk_percent = 13.0, vkr_percent = 0.0, pfe_kw = 0.0, i0_percent = 0.0)
  # Trafo 10-2 (bus2, bus10, sn_mva=1000, vn_hv_kv=220, vn_lv_kv=22, vkr_percent=0.0, vk_percent=13.0, pfe_kw=0, i0_percent=0, shift_degree=330.0)
  add2WTrafo!(net = net, fromBus = "B2", toBus = "B10", sn_mva = 1000.0, vk_percent = 13.0, vkr_percent = 0.0, pfe_kw = 0.0, i0_percent = 0.0)
  # Trafo 11-3 (bus3, bus11, sn_mva=1000, vn_hv_kv=220, vn_lv_kv=22, vkr_percent=0.0, vk_percent=13.0, pfe_kw=0, i0_percent=0, shift_degree=330.0)
  add2WTrafo!(net = net, fromBus = "B3", toBus = "B11", sn_mva = 1000.0, vk_percent = 13.0, vkr_percent = 0.0, pfe_kw = 0.0, i0_percent = 0.0)
  # 'Trafo 12-6b' (bus6b, bus12, sn_mva=500, vn_hv_kv=220, vn_lv_kv=22, vkr_percent=0.0, vk_percent=13.0, pfe_kw=0, i0_percent=0, shift_degree=330.0)
  add2WTrafo!(net = net, fromBus = "B6b", toBus = "B12", sn_mva = 500.0, vk_percent = 13.0, vkr_percent = 0.0, pfe_kw = 0.0, i0_percent = 0.0)

  # Loads
  # 'Load 2' (bus2, p_mw=285, q_mvar=200)  
  # addProsumer!(net = net, busName = "B3", type = "ENERGYCONSUMER", p = 285.0, q = 200.0)
  addProsumer!(net = net, busName = "B2", type = "ENERGYCONSUMER", p = 285.0, q = 200.0)
  # 'Load 3' (bus3, p_mw=103.0, q_mvar=62.0)
  addProsumer!(net = net, busName = "B3", type = "ENERGYCONSUMER", p = 80.0, q = 30.0)
  # 'Load 4' (net_cigre_hv, bus4, p_mw=326.0, q_mvar=244.0, name='Load 4')
  addProsumer!(net = net, busName = "B4", type = "ENERGYCONSUMER", p = 320.0, q = 200.0)
  # 'Load 5' (net_cigre_hv, bus5, p_mw=103, q_mvar=62, name='Load 5')	
  addProsumer!(net = net, busName = "B5", type = "ENERGYCONSUMER", p = 80.0, q = 30.0)
  # 'Load 6a' (net_cigre_hv, bus6a, p_mw=435, q_mvar=296, name='Load 6a')
  addProsumer!(net = net, busName = "B6a", type = "ENERGYCONSUMER", p = 435.0, q = 296.0)

  # Generators
  # 'Generator 9' (bus9, vm_pu=1.03, va_degree=0)
  addProsumer!(net = net, busName = "B9", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.03, va_deg = 0.0, referencePri = "B9")
  # 'Generator 10' (bus10, vm_pu=1.03, p_mw=500)
  addProsumer!(net = net, busName = "B10", type = "SYNCHRONOUSMACHINE", p = 500.0, vm_pu = 1.03, va_deg = 0.0)
  # 'Generator 11' (bus11, vm_pu=1.03, p_mw=200)
  addProsumer!(net = net, busName = "B11", type = "SYNCHRONOUSMACHINE", p = 200.0, vm_pu = 1.03, va_deg = 0.0)
  # 'Generator 12' (bus12, vm_pu=1.03, p_mw=300)
  addProsumer!(net = net, busName = "B12", type = "SYNCHRONOUSMACHINE", p = 300.0, vm_pu = 1.03, va_deg = 0.0)

  return net
end

function createTest5BusNet(; cooldown = 0, hyst_pu = 0.0, qlim_min = nothing, qlim_max = nothing, pq_only = false, mul_gens::Bool = true)::Net
  Sbase_MVA = 100.0
  netName = "test5bus"
  r = 0.05
  x = 0.5
  c_nf_per_km = 10.0
  tanδ = 0.0

  @debug "Creating $netName test network with qlim_min=$qlim_min, qlim_max=$qlim_max, pq_only=$pq_only"

  Bus5Net = Net(name = netName, baseMVA = Sbase_MVA, cooldown_iters = cooldown, q_hyst_pu = hyst_pu)
  addBus!(net = Bus5Net, busName = "B1", vn_kV = 110.0)
  addBus!(net = Bus5Net, busName = "B2", vn_kV = 110.0)
  addBus!(net = Bus5Net, busName = "B3", vn_kV = 110.0)
  addBus!(net = Bus5Net, busName = "B4", vn_kV = 110.0)
  addBus!(net = Bus5Net, busName = "B5", vn_kV = 110.0)

  addACLine!(net = Bus5Net, fromBus = "B1", toBus = "B3", length = 20.0, r = r, x = x, c_nf_per_km = c_nf_per_km, tanδ = tanδ)
  addACLine!(net = Bus5Net, fromBus = "B3", toBus = "B5", length = 50.0, r = r, x = x, c_nf_per_km = c_nf_per_km, tanδ = tanδ)
  addACLine!(net = Bus5Net, fromBus = "B5", toBus = "B4", length = 100.0, r = r, x = x, c_nf_per_km = c_nf_per_km, tanδ = tanδ)
  addACLine!(net = Bus5Net, fromBus = "B4", toBus = "B1", length = 100.0, r = r, x = x, c_nf_per_km = c_nf_per_km, tanδ = tanδ)
  addACLine!(net = Bus5Net, fromBus = "B4", toBus = "B2", length = 100.0, r = r, x = x, c_nf_per_km = c_nf_per_km, tanδ = tanδ)
  addACLine!(net = Bus5Net, fromBus = "B2", toBus = "B1", length = 20.0, r = r, x = x, c_nf_per_km = c_nf_per_km, tanδ = tanδ)

  addProsumer!(net = Bus5Net, busName = "B1", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "B1")
  if mul_gens
    addProsumer!(net = Bus5Net, busName = "B2", type = "ENERGYCONSUMER", p = 25.0, q = 8.0)
    addProsumer!(net = Bus5Net, busName = "B2", type = "ENERGYCONSUMER", p = 25.0, q = 7.0)
    vm_b3 = pq_only ? nothing : 1.03
    addProsumer!(net = Bus5Net, busName = "B3", type = "SYNCHRONOUSMACHINE", p = 20.0, q = 15.0, vm_pu = vm_b3, va_deg = 0.0, qMax = 0.6 * qlim_max, qMin = 0.6 * qlim_min)
    addProsumer!(net = Bus5Net, busName = "B3", type = "SYNCHRONOUSMACHINE", p = 10.0, q = 10.0, vm_pu = vm_b3, va_deg = 0.0, qMax = 0.4 * qlim_max, qMin = 0.4 * qlim_min)
    addProsumer!(net = Bus5Net, busName = "B4", type = "ENERGYCONSUMER", p = 25.0, q = 7.0)
    addProsumer!(net = Bus5Net, busName = "B4", type = "ENERGYCONSUMER", p = 25.0, q = 8.0)
    addProsumer!(net = Bus5Net, busName = "B5", type = "ENERGYCONSUMER", p = 10.0, q = 5.0)
    addProsumer!(net = Bus5Net, busName = "B5", type = "ENERGYCONSUMER", p = 15.0, q = 1.0)
    if !pq_only
      addProsumer!(net = Bus5Net, busName = "B5", type = "SYNCHRONOUSMACHINE", p = 0.0, vm_pu = 1.0, qMax = qlim_max, qMin = qlim_min)
    end
  else
    addProsumer!(net = Bus5Net, busName = "B2", type = "ENERGYCONSUMER", p = 50.0, q = 15.0)
    vm_b3 = pq_only ? nothing : 1.0
    addProsumer!(net = Bus5Net, busName = "B3", type = "SYNCHRONOUSMACHINE", p = 30.0, q = 25.0, vm_pu = vm_b3, va_deg = 0.0, qMax = qlim_max, qMin = qlim_min)
    addProsumer!(net = Bus5Net, busName = "B4", type = "ENERGYCONSUMER", p = 50.0, q = 15.0)
    addProsumer!(net = Bus5Net, busName = "B5", type = "ENERGYCONSUMER", p = 25.0, q = 10.0)
  end

  return Bus5Net
end

function createTest3BusNet(; cooldown = 0, hyst_pu = 0.0, qlim_min = nothing, qlim_max = nothing)::Net
  # Simple 3-bus network
  #
  #  ASTADT        STATION1
  # <--|---------------|<--- Generator 
  #    |-------       |
  #            |      |
  #            --------|<---- EXTERNALNETWORKINJECTION
  #                 VERBUND  
  Sbase_MVA = 100.0
  netName = "test3bus"

  r = 0.0
  x = 0.4
  s = 25.0
  c_nf_per_km = 9.55
  tanδ = 0.0

  vm_pu_STATION1 = 1.027273
  vm_pu_VERBUND = 1.018182

  @debug "Creating $netName test network with qlim_min=$qlim_min, qlim_max=$qlim_max"

  Bus3Net = Net(name = netName, baseMVA = Sbase_MVA, cooldown_iters = cooldown, q_hyst_pu = hyst_pu)

  addBus!(net = Bus3Net, busName = "ASTADT", vn_kV = 110.0)
  addBus!(net = Bus3Net, busName = "STATION1", vn_kV = 110.0)
  addBus!(net = Bus3Net, busName = "VERBUND", vn_kV = 110.0)

  addACLine!(net = Bus3Net, fromBus = "ASTADT", toBus = "STATION1", length = s, r = r, x = x, c_nf_per_km = c_nf_per_km, tanδ = tanδ)
  addACLine!(net = Bus3Net, fromBus = "ASTADT", toBus = "VERBUND", length = s, r = r, x = x, c_nf_per_km = c_nf_per_km, tanδ = tanδ)
  addACLine!(net = Bus3Net, fromBus = "VERBUND", toBus = "STATION1", length = s, r = r, x = x, c_nf_per_km = c_nf_per_km, tanδ = tanδ)

  addProsumer!(net = Bus3Net, busName = "VERBUND", type = "EXTERNALNETWORKINJECTION", vm_pu = vm_pu_VERBUND, va_deg = 0.0, referencePri = "VERBUND")
  addProsumer!(net = Bus3Net, busName = "STATION1", type = "SYNCHRONOUSMACHINE", p = 70.0, q = 33.2, vm_pu = vm_pu_STATION1, qMax = qlim_max, qMin = qlim_min)
  addProsumer!(net = Bus3Net, busName = "ASTADT", type = "ENERGYCONSUMER", p = 100.0, q = 30.0)

  return Bus3Net
end

function createTest2BusNet(; cooldown = 0, hyst_pu = 0.0, qlim_min = nothing, qlim_max = nothing)::Net
  # Simple 2-bus network
  #   Bus1            Bus2
  #  ->|---------------|<--- Generator
  #   Slack           PV
  Sbase_MVA = 100.0
  netName = "test2bus"
  r = 0.05
  x = 0.5
  c_nf_per_km = 10.0
  tanδ = 0.0

  @debug "Creating $netName test network with qlim_min=$qlim_min, qlim_max=$qlim_max"

  Bus2Net = Net(name = netName, baseMVA = Sbase_MVA, cooldown_iters = cooldown, q_hyst_pu = hyst_pu)
  addBus!(net = Bus2Net, busName = "B1", vn_kV = 110.0)
  addBus!(net = Bus2Net, busName = "B2", vn_kV = 110.0)

  addACLine!(net = Bus2Net, fromBus = "B1", toBus = "B2", length = 20.0, r = r, x = x, c_nf_per_km = c_nf_per_km, tanδ = tanδ)

  addProsumer!(net = Bus2Net, busName = "B1", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "B1")
  addProsumer!(net = Bus2Net, busName = "B2", type = "SYNCHRONOUSMACHINE", p = 20.0, vm_pu = 1.05, va_deg = 0.0, qMax = qlim_max, qMin = qlim_min)

  return Bus2Net
end

function test_3BusNet(verbose::Int = 0, qlim::Float64 = 15.0, method::Symbol = :rectangular, opt_fd::Bool = false, opt_sparse::Bool = false)
  net = createTest3BusNet(cooldown = 2, hyst_pu = 0.01, qlim_min = -qlim, qlim_max = qlim)
  tol = 1e-12
  maxIte = 50
  print_results = (verbose > 0)
  result = true

  pv_names = ["STATION1"]
  etim = 0.0
  etim = @elapsed begin
    ite, erg = runpf!(net, maxIte, tol, verbose, method = method, opt_fd = opt_fd, opt_sparse = opt_sparse)
    if erg != 0
      @info "Full-system power flow did not converge"
      result = false
    end
  end

  #

  hit = pv_hit_q_limit(net, pv_names)
  calcNetLosses!(net)
  distributeBusResults!(net)
  if print_results
    printACPFlowResults(net, etim, ite, tol; converged = result, solver = method)
    printProsumerResults(net)
    printQLimitLog(net; sort_by = :bus)
  end

  tp = getTotalBusPower(net = net)
  tl = getTotalLosses(net = net)
  if verbose > 1
    @info "Total Power: P=$(tp[1]), Q=$(tp[2])"
    @info "Total Losses: P=$(tl[1]), Q=$(tl[2])"
  end

  @test all(isapprox.(tp, tl; atol = 1e-6))

  if qlim < 33.2
    return hit == true
  else
    return hit == false
  end
end

function test_5BusNet(verbose::Int = 0, qlim::Float64 = 20.0, method::Symbol = :rectangular, opt_fd::Bool = false, opt_sparse::Bool = false)
  net = createTest5BusNet(cooldown = 2, hyst_pu = 0.01, qlim_min = -qlim, qlim_max = qlim)
  tol = 1e-9
  maxIte = 50
  print_results = (verbose > 0)
  result = true

  pv_names = ["B3"]
  etim = 0.0
  etim = @elapsed begin
    ite, erg = runpf!(net, maxIte, tol, verbose, method = method, opt_fd = opt_fd, opt_sparse = opt_sparse)
    if erg != 0
      @info "Full-system power flow did not converge"
      result = false
    end
  end

  #

  hit = pv_hit_q_limit(net, pv_names)

  if print_results
    V = buildVoltageVector(net)
    calcNetLosses!(net, V)
    distributeBusResults!(net)
    printACPFlowResults(net, etim, ite, tol; converged = result, solver = method)
    printProsumerResults(net)
    printQLimitLog(net; sort_by = :bus)
  end

  return hit == true
end

# -----------------------------------------------------------------------------
# Test intent: MATPOWER inline case vs manually assembled Net (with bus shunts)
#
# Why this exists:
# - We want to verify that `createNetFromMatPowerCase` builds the same electrical
#   model as the explicit Sparlectra API calls (`addBus!`, `addACLine!`,
#   `addProsumer!`, `addShunt!`).
# - This specifically guards the bus-shunt mapping from MATPOWER bus columns
#   `Gs/Bs` into Sparlectra's internal shunt representation.
#
# What is compared:
# 1) Converged PF state (complex bus voltages V)
# 2) Slack-bus complex injection (P/Q)
#
# If these are equal within tolerance, both model construction paths are
# considered equivalent for this scenario.
# -----------------------------------------------------------------------------
function test_mp_inline_vs_manual_shunt(verbose::Int = 0; method::Symbol = :rectangular, opt_sparse::Bool = true)::Bool
  tol_pf = 1e-10
  maxIte = 40
  tol_cmp = 5e-4
  print_results = (verbose > 0)

  # -------------------------
  # 1) MATPOWER inline case
  # -------------------------
  mpc = (
    name = "case3_mp_vs_manual_shunt",
    baseMVA = 100.0,

    # bus: [bus_i type Pd Qd Gs Bs area Vm Va baseKV zone Vmax Vmin]
    # We place non-zero Bs on buses 2 and 3 and keep line charging at 0.
    # This isolates the test to *bus shunts* (no branch-charging side effects).
    bus = [
      1 3 0.0 0.0 0.0 0.0 1 1.06 0.0 110.0 1 1.10 0.90
      2 2 20.0 5.0 0.0 6.0 1 1.03 0.0 110.0 1 1.10 0.90
      3 1 90.0 30.0 0.0 12.0 1 1.00 0.0 110.0 1 1.10 0.90
    ],

    # gen: 21 cols (MATPOWER v2)
    # bus Pg Qg Qmax Qmin Vg mBase status Pmax Pmin ...
    gen = [
      1 120.0 0.0 300.0 -300.0 1.06 100.0 1 300.0 0.0 0 0 0 0 0 0 0 0 0 0 0
      2 60.0 0.0 200.0 -200.0 1.03 100.0 1 200.0 0.0 0 0 0 0 0 0 0 0 0 0 0
    ],

    # branch: [fbus tbus r x b rateA rateB rateC ratio angle status angmin angmax]
    # NOTE: line charging b=0 here (we test shunts via bus Gs/Bs).
    branch = [
      1 2 0.02 0.06 0.0 0 0 0 0.0 0.0 1 -360 360
      1 3 0.08 0.24 0.0 0 0 0 0.0 0.0 1 -360 360
      2 3 0.06 0.18 0.0 0 0 0 0.0 0.0 1 -360 360
    ],
    gencost = nothing,
    bus_name = ["BUS1 (Slack)", "BUS2 (PV)", "BUS3 (PQ)"],
  )

  # Build network via MATPOWER import path.
  net_mp = Sparlectra.createNetFromMatPowerCase(mpc = mpc, log = false, flatstart = false)

  ok, msg = validate!(net = net_mp)
  ok || (@warn msg; return false)

  # Run PF for MATPOWER-derived net.
  erg_mp = 0
  ite_mp = 0
  et_mp = @elapsed begin
    ite_mp, erg_mp = runpf!(net_mp, maxIte, tol_pf, verbose; method = method, opt_sparse = opt_sparse)
  end
  if erg_mp != 0
    @warn "MATPOWER-derived net: power flow did not converge"
    return false
  end

  updateShuntPowers!(net = net_mp)

  Y_mp = createYBUS(net = net_mp, sparse = opt_sparse, printYBUS = false)
  V_mp = buildVoltageVector(net_mp)
  Sinj_mp = calc_injections(Y_mp, V_mp)
  # -------------------------
  # 2) Manual Net build
  # -------------------------
  netName = "case3_manual_vs_mp_shunt"
  net_man = Net(name = netName, baseMVA = 100.0)
  addBus!(net = net_man, busName = "B1", vn_kV = 110.0, vm_pu = 1.06, va_deg = 0.0)
  addBus!(net = net_man, busName = "B2", vn_kV = 110.0, vm_pu = 1.03, va_deg = 0.0)
  addBus!(net = net_man, busName = "B3", vn_kV = 110.0, vm_pu = 1.00, va_deg = 0.0)
  # Lines without line charging, so only bus shunts differ/drive this check.
  addACLine!(net = net_man, fromBus = "B1", toBus = "B2", length = 1.0, r = 0.02, x = 0.06)
  addACLine!(net = net_man, fromBus = "B1", toBus = "B3", length = 1.0, r = 0.08, x = 0.24)
  addACLine!(net = net_man, fromBus = "B2", toBus = "B3", length = 1.0, r = 0.06, x = 0.18)

  # Loads from MATPOWER Pd/Qd mapped to ENERGYCONSUMER.
  addProsumer!(net = net_man, busName = "B3", type = "ENERGYCONSUMER", p = 90.0, q = 30.0)

  # MATPOWER bus 2 also has load: Pd=20 MW, Qd=5 MVAr
  addProsumer!(net = net_man, busName = "B2", type = "ENERGYCONSUMER", p = 20.0, q = 5.0)

  addProsumer!(net = net_man, busName = "B1", type = "EXTERNALNETWORKINJECTION", referencePri = "B1")

  # PV generator at B2 (P + Vm setpoint)
  addProsumer!(net = net_man, busName = "B2", type = "SYNCHRONOUSMACHINE", p = 60.0, vm_pu = 1.03, va_deg = 0.0, qMin = -200.0, qMax = 200.0)

  # Manual shunts chosen to match the MATPOWER case after import.
  # NOTE: this test validates *equivalence of assembled network models*,
  # not a standalone sign-convention proof. The numerical checks below
  # (V and slack P/Q) are the actual acceptance criteria.
  addShunt!(net = net_man, busName = "B2", pShunt = 0.0, qShunt = 6.0, in_service = 1)
  addShunt!(net = net_man, busName = "B3", pShunt = 0.0, qShunt = 12.0, in_service = 1)

  ok2, msg2 = validate!(net = net_man)
  ok2 || (@warn msg2; return false)

  erg_man = 0
  ite_man = 0

  et_man = @elapsed begin
    ite_man, erg_man = runpf!(net_man, maxIte, tol_pf, verbose; method = method, opt_sparse = opt_sparse)
  end
  if erg_man != 0
    @warn "Manual net: power flow did not converge"
    return false
  end
  updateShuntPowers!(net = net_man)

  Y_man = createYBUS(net = net_man, sparse = opt_sparse, printYBUS = false)
  V_man = buildVoltageVector(net_man)
  Sinj_man = calc_injections(Y_man, V_man)
  # -------------------------
  # 3) Compare results
  # -------------------------
  # Compare complex voltages directly (magnitude and angle are implicit in V).
  if length(V_mp) != length(V_man)
    @warn "Voltage vector length mismatch" n_mp = length(V_mp) n_man = length(V_man)
    return false
  end

  # Bus order assumption:
  # - MATPOWER buses are indexed 1..3
  # - Manual net is created in order B1, B2, B3
  # So vectors are expected to be aligned index-wise.
  if !all(isapprox.(V_mp, V_man; atol = tol_cmp, rtol = 0.0))
    @warn "Voltage mismatch between MATPOWER-derived and manual net" V_mp = V_mp V_man = V_man
    return false
  end

  # Compare slack injections (P,Q) in MW/MVAr as system-level consistency check.
  slack_idx = 1
  Psl_mp = real(Sinj_mp[slack_idx]) * mpc.baseMVA
  Qsl_mp = imag(Sinj_mp[slack_idx]) * mpc.baseMVA
  Psl_man = real(Sinj_man[slack_idx]) * mpc.baseMVA
  Qsl_man = imag(Sinj_man[slack_idx]) * mpc.baseMVA

  if !(isapprox(Psl_mp, Psl_man; atol = 1e-4, rtol = 0.0) && isapprox(Qsl_mp, Qsl_man; atol = 1e-4, rtol = 0.0))
    @warn "Slack P/Q mismatch" Psl_mp = Psl_mp Qsl_mp = Qsl_mp Psl_man = Psl_man Qsl_man = Qsl_man
    return false
  end

  if print_results
    @info "Qsh MAN" B2 = net_man.nodeVec[2]._qShunt B3 = net_man.nodeVec[3]._qShunt
    @info "Qsh MP" B2 = net_mp.nodeVec[2]._qShunt B3 = net_mp.nodeVec[3]._qShunt
    @info "manual setpoints" B1_vm = net_man.nodeVec[1]._vm_pu B1_va = net_man.nodeVec[1]._va_deg B2_vm = net_man.nodeVec[2]._vm_pu
    @info "mp setpoints" B1_vm = net_mp.nodeVec[1]._vm_pu B1_va = net_mp.nodeVec[1]._va_deg B2_vm = net_mp.nodeVec[2]._vm_pu

    @info "MP net PF ok" et = et_mp ite = ite_mp
    @info "MAN net PF ok" et = et_man ite = ite_man
    printACPFlowResults(net_mp, et_mp, ite_mp, tol_pf; converged = true, solver = method)
    printACPFlowResults(net_man, et_man, ite_man, tol_pf; converged = true, solver = method)
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

  S1, dP1, dQ1 = buildControlledSVec(net_vdi, ComplexF64[1.0 + 0.0im])
  S2, dP2, dQ2 = buildControlledSVec(net_vdi, ComplexF64[2.0 + 0.0im])

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
  return isapprox(l.pFlow_MW, 10.0; atol = 1e-8) && isapprox(l.qFlow_MVar, 2.0; atol = 1e-8) && isapprox(l.iFrom_kA, iexp; atol = 1e-10) && isapprox(l.iTo_kA, iexp; atol = 1e-10)
end

function test_link_component_type()
  net = Net(name = "link_component", baseMVA = 100.0)
  addBus!(net = net, busName = "B1", vn_kV = 110.0)
  addBus!(net = net, busName = "B2", vn_kV = 110.0)
  addLink!(net = net, fromBus = "B1", toBus = "B2", status = 1)

  l = net.linkVec[1]
  return (l isa AbstractComponent) && l.cTyp == Sparlectra.LinkC && toComponentTyp("BUSLINK") == Sparlectra.LinkC && toComponentTyp("LINK") == Sparlectra.LinkC
end

function test_link_rejects_slack_connection()
  net = Net(name = "link_slack_reject", baseMVA = 100.0)
  addBus!(net = net, busName = "B1", vn_kV = 110.0)
  addBus!(net = net, busName = "B2", vn_kV = 110.0)
  addProsumer!(net = net, busName = "B1", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "B1")
  addProsumer!(net = net, busName = "B2", type = "ENERGYCONSUMER", p = 5.0, q = 1.0)

  try
    addLink!(net = net, fromBus = "B1", toBus = "B2", status = 1)
    return false
  catch err
    return err isa AssertionError && occursin("slack", lowercase(string(err)))
  end
end

function test_link_rejects_mixed_bus_types()
  net = Net(name = "link_type_reject", baseMVA = 100.0)
  addBus!(net = net, busName = "B1", vn_kV = 110.0)
  addBus!(net = net, busName = "B2", vn_kV = 110.0)
  addProsumer!(net = net, busName = "B1", type = "SYNCHRONOUSMACHINE", p = 10.0, vm_pu = 1.02)
  addProsumer!(net = net, busName = "B2", type = "ENERGYCONSUMER", p = 5.0, q = 1.0)

  try
    addLink!(net = net, fromBus = "B1", toBus = "B2", status = 1)
    return false
  catch err
    return err isa AssertionError && occursin("same bus type", lowercase(string(err)))
  end
end

function test_link_bus_merge_pf()
  tol = 1e-8
  maxIte = 30

  net_link = Net(name = "link_merge_pf", baseMVA = 100.0)
  addBus!(net = net_link, busName = "B0", vn_kV = 110.0)
  addBus!(net = net_link, busName = "B1", vn_kV = 110.0)
  addBus!(net = net_link, busName = "B2", vn_kV = 110.0)
  addBus!(net = net_link, busName = "B3", vn_kV = 110.0)

  addACLine!(net = net_link, fromBus = "B0", toBus = "B1", length = 5.0, r = 0.05, x = 0.5, c_nf_per_km = 10.0, tanδ = 0.0)
  addACLine!(net = net_link, fromBus = "B2", toBus = "B3", length = 20.0, r = 0.05, x = 0.5, c_nf_per_km = 10.0, tanδ = 0.0)
  addLink!(net = net_link, fromBus = "B1", toBus = "B2", status = 1)

  addProsumer!(net = net_link, busName = "B0", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "B0")
  addProsumer!(net = net_link, busName = "B1", type = "GENERATOR", p = 10.0, q = 1.0)
  addProsumer!(net = net_link, busName = "B2", type = "ENERGYCONSUMER", p = 15.0, q = 5.0)
  addProsumer!(net = net_link, busName = "B3", type = "ENERGYCONSUMER", p = 25.0, q = 10.0)

  ite_link, erg_link = runpf!(net_link, maxIte, tol, 0; method = :polar_full, opt_sparse = true)
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

  addACLine!(net = net_eq, fromBus = "B0", toBus = "B12", length = 5.0, r = 0.05, x = 0.5, c_nf_per_km = 10.0, tanδ = 0.0)
  addACLine!(net = net_eq, fromBus = "B12", toBus = "B3", length = 20.0, r = 0.05, x = 0.5, c_nf_per_km = 10.0, tanδ = 0.0)
  addProsumer!(net = net_eq, busName = "B0", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "B0")
  addProsumer!(net = net_eq, busName = "B12", type = "GENERATOR", p = 10.0, q = 1.0)
  addProsumer!(net = net_eq, busName = "B12", type = "ENERGYCONSUMER", p = 15.0, q = 5.0)
  addProsumer!(net = net_eq, busName = "B3", type = "ENERGYCONSUMER", p = 25.0, q = 10.0)

  ite_eq, erg_eq = runpf!(net_eq, maxIte, tol, 0; method = :polar_full, opt_sparse = true)
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

function test_link_bus_merge_pf_default_method()
  return test_link_bus_merge_pf()
end

function test_rectangular_merge_fallback_suppresses_polar_deprecation_warning()::Bool
  net = Net(name = "link_merge_warning_gate", baseMVA = 100.0)
  addBus!(net = net, busName = "B0", vn_kV = 110.0)
  addBus!(net = net, busName = "B1", vn_kV = 110.0)
  addBus!(net = net, busName = "B2", vn_kV = 110.0)
  addBus!(net = net, busName = "B3", vn_kV = 110.0)

  addACLine!(net = net, fromBus = "B0", toBus = "B1", length = 5.0, r = 0.05, x = 0.5, c_nf_per_km = 10.0, tanδ = 0.0)
  addACLine!(net = net, fromBus = "B2", toBus = "B3", length = 20.0, r = 0.05, x = 0.5, c_nf_per_km = 10.0, tanδ = 0.0)
  addLink!(net = net, fromBus = "B1", toBus = "B2", status = 1)
  addProsumer!(net = net, busName = "B0", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "B0")
  addProsumer!(net = net, busName = "B1", type = "GENERATOR", p = 10.0, q = 1.0)
  addProsumer!(net = net, busName = "B2", type = "ENERGYCONSUMER", p = 15.0, q = 5.0)
  addProsumer!(net = net, busName = "B3", type = "ENERGYCONSUMER", p = 25.0, q = 10.0)

  old_warned = Sparlectra._warned_full_solver_deprecated[]
  Sparlectra._warned_full_solver_deprecated[] = false
  try
    @test_logs min_level = Logging.Warn (:warn, r"runpf!: rectangular solver detected internal Isolated buses from active-link merges; using rectangular FD fallback instead of :polar_full") match_mode = :all begin
      _, erg = redirect_stdout(devnull) do
        runpf!(net, 30, 1e-8, 1; method = :rectangular, opt_sparse = true)
      end
      if erg != 0
        return false
      end
    end
  finally
    Sparlectra._warned_full_solver_deprecated[] = old_warned
  end

  return true
end

function test_link_closed_keeps_shunt_reporting_on_original_bus()
  net = Net(name = "link_shunt_reporting", baseMVA = 100.0)
  addBus!(net = net, busName = "B0", vn_kV = 110.0)
  addBus!(net = net, busName = "B1", vn_kV = 110.0)
  addBus!(net = net, busName = "B2", vn_kV = 110.0)
  addBus!(net = net, busName = "B3", vn_kV = 110.0)

  addACLine!(net = net, fromBus = "B0", toBus = "B1", length = 5.0, r = 0.05, x = 0.5, c_nf_per_km = 10.0, tanδ = 0.0)
  addACLine!(net = net, fromBus = "B2", toBus = "B3", length = 20.0, r = 0.05, x = 0.5, c_nf_per_km = 10.0, tanδ = 0.0)
  addLink!(net = net, fromBus = "B1", toBus = "B2", status = 1)
  addShunt!(net = net, busName = "B2", pShunt = 0.0, qShunt = 5.0)

  addProsumer!(net = net, busName = "B0", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "B0")
  addProsumer!(net = net, busName = "B1", type = "GENERATOR", p = 10.0, q = 1.0)
  addProsumer!(net = net, busName = "B2", type = "ENERGYCONSUMER", p = 15.0, q = 5.0)
  addProsumer!(net = net, busName = "B3", type = "ENERGYCONSUMER", p = 25.0, q = 10.0)

  _, erg = runpf!(net, 25, 1e-8, 0; method = :polar_full, opt_sparse = true)
  erg == 0 || return false

  b2_idx = geNetBusIdx(net = net, busName = "B2")
  sh = net.shuntVec[1]
  vm = net.nodeVec[b2_idx]._vm_pu
  expected_q = -vm^2 * 5.0

  return isapprox(sh.q_shunt, expected_q; atol = 1e-8) && isapprox(net.nodeVec[b2_idx]._qShunt, expected_q; atol = 1e-8)
end

function test_link_kcl_ring_allocation()
  net = Net(name = "link_kcl_ring", baseMVA = 100.0)
  # 3 busses
  addBus!(net = net, busName = "B1", vn_kV = 110.0)
  addBus!(net = net, busName = "B2", vn_kV = 110.0)
  addBus!(net = net, busName = "B3", vn_kV = 110.0)

  # Zero-Bilance: + 15 -10 -5 = 0;   +6 -4 -2 = 0
  addBusGenPower!(net = net, busName = "B1", p = 15.0, q = 6.0)
  addBusLoadPower!(net = net, busName = "B2", p = 10.0, q = 4.0)
  addBusLoadPower!(net = net, busName = "B3", p = 5.0, q = 2.0)

  # B1 -> B2 -> B3 -> B1 forms a ring
  addLink!(net = net, fromBus = "B1", toBus = "B2", status = 1)
  addLink!(net = net, fromBus = "B2", toBus = "B3", status = 1)
  addLink!(net = net, fromBus = "B3", toBus = "B1", status = 1)

  # no power flow needed for KCL test
  calcLinkFlowsKCL!(net)

  # manual KCL solution for ring with 3 busses and 3 links:
  # formula: A * f = b, where A is the incidence matrix of the ring and b is the net injection vector (P or Q)
  # in a ring A is singular, so we use the pseudo-inverse to get the minimum-norm solution for the flows
  A = [
    1.0 0.0 -1.0
    -1.0 1.0 0.0
    0.0 -1.0 1.0
  ]
  bP = [15.0, -10.0, -5.0]
  bQ = [6.0, -4.0, -2.0]

  # calculate expected flows using pseudo-inverse (since A is singular for a ring, but we want the minimum-norm solution)
  fP_expected = pinv(A) * bP
  fQ_expected = pinv(A) * bQ

  flowsP = [net.linkVec[1].pFlow_MW, net.linkVec[2].pFlow_MW, net.linkVec[3].pFlow_MW]
  flowsQ = [net.linkVec[1].qFlow_MVar, net.linkVec[2].qFlow_MVar, net.linkVec[3].qFlow_MVar]

  return isapprox.(flowsP, fP_expected; atol = 1e-10, rtol = 0.0) |> all && isapprox.(flowsQ, fQ_expected; atol = 1e-10, rtol = 0.0) |> all && isapprox.(A * flowsP, bP; atol = 1e-10, rtol = 0.0) |> all && isapprox.(A * flowsQ, bQ; atol = 1e-10, rtol = 0.0) |> all
end

function test_link_kcl_ring_allocation_with_shunt()
  net = Net(name = "link_kcl_ring_shunt", baseMVA = 100.0)
  addBus!(net = net, busName = "B1", vn_kV = 110.0)
  addBus!(net = net, busName = "B2", vn_kV = 110.0)
  addBus!(net = net, busName = "B3", vn_kV = 110.0)

  addBusGenPower!(net = net, busName = "B1", p = 15.0, q = 12.0)
  addBusLoadPower!(net = net, busName = "B2", p = 10.0, q = 4.0)
  addBusLoadPower!(net = net, busName = "B3", p = 5.0, q = 2.0)

  # Y-model shunt at B2; at vm=1.0 pu this contributes Q_shunt = -6 MVar.
  addShunt!(net = net, busName = "B2", pShunt = 0.0, qShunt = 6.0)
  net.nodeVec[1]._vm_pu = 1.0
  net.nodeVec[2]._vm_pu = 1.0
  net.nodeVec[3]._vm_pu = 1.0
  updateShuntPowers!(net = net)

  addLink!(net = net, fromBus = "B1", toBus = "B2", status = 1)
  addLink!(net = net, fromBus = "B2", toBus = "B3", status = 1)
  addLink!(net = net, fromBus = "B3", toBus = "B1", status = 1)

  calcLinkFlowsKCL!(net)

  A = [
    1.0 0.0 -1.0
    -1.0 1.0 0.0
    0.0 -1.0 1.0
  ]
  bP = [15.0, -10.0, -5.0]
  bQ = [12.0, -10.0, -2.0]
  fP_expected = pinv(A) * bP
  fQ_expected = pinv(A) * bQ

  flowsP = [net.linkVec[1].pFlow_MW, net.linkVec[2].pFlow_MW, net.linkVec[3].pFlow_MW]
  flowsQ = [net.linkVec[1].qFlow_MVar, net.linkVec[2].qFlow_MVar, net.linkVec[3].qFlow_MVar]

  return isapprox.(flowsP, fP_expected; atol = 1e-10, rtol = 0.0) |> all && isapprox.(flowsQ, fQ_expected; atol = 1e-10, rtol = 0.0) |> all && isapprox.(A * flowsP, bP; atol = 1e-10, rtol = 0.0) |> all && isapprox.(A * flowsQ, bQ; atol = 1e-10, rtol = 0.0) |> all
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
  addACLine!(net = net, fromBus = "B0", toBus = "B1", length = 10.0, r = 0.05, x = 0.4, c_nf_per_km = 10.0, tanδ = 0.0)
  addACLine!(net = net, fromBus = "B0", toBus = "B2", length = 12.0, r = 0.05, x = 0.4, c_nf_per_km = 10.0, tanδ = 0.0)

  # Three links form a closed loop (meshed/singular incidence structure).
  addLink!(net = net, fromBus = "B1", toBus = "B2", status = 1)
  addLink!(net = net, fromBus = "B2", toBus = "B3", status = 1)
  addLink!(net = net, fromBus = "B3", toBus = "B1", status = 1)

  addProsumer!(net = net, busName = "B0", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "B0")
  addProsumer!(net = net, busName = "B1", type = "GENERATOR", p = 6.0, q = 1.0)
  addProsumer!(net = net, busName = "B2", type = "ENERGYCONSUMER", p = 18.0, q = 5.0)
  addProsumer!(net = net, busName = "B3", type = "ENERGYCONSUMER", p = 8.0, q = 3.0)

  # Solve PF first, then derive link flows from nodal balances.
  _, erg = runpf!(net, 40, 1e-9, 0; method = :polar_full, opt_sparse = false)
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
  #   b = (gen - load + shunt) - branch_out
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

  bP = P_inj .+ P_shunt .- P_branch_out
  bQ = Q_inj .+ Q_shunt .- Q_branch_out
  # Remove any tiny component residual to match solver behavior.
  bP .-= sum(bP) / 3
  bQ .-= sum(bQ) / 3

  flowsP = [l.pFlow_MW for l in net.linkVec]
  flowsQ = [l.qFlow_MVar for l in net.linkVec]

  return vm_equal && va_equal && isapprox.(A * flowsP, bP; atol = 1e-7, rtol = 0.0) |> all && isapprox.(A * flowsQ, bQ; atol = 1e-7, rtol = 0.0) |> all
end

function test_acpflow_report_object()
  net = createTest3BusNet(cooldown = 2, hyst_pu = 0.01, qlim_min = -15.0, qlim_max = 15.0)
  ite, erg = runpf!(net, 50, 1e-10, 0; method = :polar_full, opt_sparse = true)
  erg == 0 || return false

  calcNetLosses!(net)
  calcLinkFlowsKCL!(net)

  report = buildACPFlowReport(net; ct = 0.123, ite = ite, tol = 1e-10, converged = true, solver = :polar_full)

  node_count_ok = length(report.nodes) == length(net.nodeVec)
  branch_count_ok = length(report.branches) == length(net.branchVec)
  link_count_ok = length(report.links) == length(net.linkVec)
  ctrl_count_ok = length(report.transformer_controls) == 0

  total_losses_ok = isapprox(report.metadata.total_p_loss_MW, getTotalLosses(net = net)[1]; atol = 1e-8) && isapprox(report.metadata.total_q_loss_MVar, getTotalLosses(net = net)[2]; atol = 1e-8)

  first_node = report.nodes[1]
  first_branch = report.branches[1]

  has_node_keys = haskey(first_node, :bus) && haskey(first_node, :vm_pu) && haskey(first_node, :q_limit_hit)
  has_branch_keys = haskey(first_branch, :from_bus) && haskey(first_branch, :p_from_MW) && haskey(first_branch, :overloaded)

  return node_count_ok && branch_count_ok && link_count_ok && ctrl_count_ok && total_losses_ok && has_node_keys && has_branch_keys
end

function test_report_uses_user_bus_names_and_pf_node_count()
  net = Net(name = "report_link_names", baseMVA = 100.0)
  addBus!(net = net, busName = "Slack", vn_kV = 110.0)
  addBus!(net = net, busName = "Bus1", vn_kV = 110.0)
  addBus!(net = net, busName = "Bus1a", vn_kV = 110.0)

  addACLine!(net = net, fromBus = "Slack", toBus = "Bus1", length = 5.0, r = 0.05, x = 0.5, c_nf_per_km = 10.0, tanδ = 0.0)
  addLink!(net = net, fromBus = "Bus1", toBus = "Bus1a", status = 1)

  addProsumer!(net = net, busName = "Slack", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "Slack")
  addProsumer!(net = net, busName = "Bus1", type = "ENERGYCONSUMER", p = 10.0, q = 2.0)
  addProsumer!(net = net, busName = "Bus1a", type = "ENERGYCONSUMER", p = 5.0, q = 1.0)

  ite, erg = runpf!(net, 30, 1e-10, 0; method = :polar_full, opt_sparse = true)
  erg == 0 || return false

  calcNetLosses!(net)
  calcLinkFlowsKCL!(net)

  report = buildACPFlowReport(net; ct = 0.0, ite = ite, tol = 1e-10, converged = true, solver = :polar_full)
  report_text = mktempdir() do tmpdir
    printACPFlowResults(net, 0.0, ite, 1e-10, true, tmpdir; converged = true, solver = :polar_full)
    read(joinpath(tmpdir, "result_$(net.name).txt"), String)
  end

  has_bus1a_row = any(row -> row.bus_name == "Bus1a", report.nodes)
  merged_pf_count_ok = occursin("PF Nodes", report_text) && occursin("2 (after active-link merge)", report_text)
  controllers_line_ok = occursin("Controllers", report_text)
  link_names_ok = occursin("Bus1", report_text) && occursin("Bus1a", report_text)
  branch_connection_ok = occursin("Slack -> Bus1", report_text) && occursin("Line", report_text)

  return has_bus1a_row && merged_pf_count_ok && controllers_line_ok && link_names_ok && branch_connection_ok
end

function test_matpower_read_case_m_postprocessing()::Bool
  mfile = tempname() * ".m"
  open(mfile, "w") do io
    write(io, """
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
""")
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

function test_q_limit_adjust_vset_success()::Bool
  net = createTest3BusNet(cooldown = 0, hyst_pu = 0.0, qlim_min = -30.0, qlim_max = 30.0)
  _assign_vset_controller!(net, "STATION1"; vstep_pu = 0.005, tap_steps_down = 2, tap_steps_up = 2)

  _, erg = runpf!(net, 40, 1e-8, 0; method = :rectangular, qlimit_mode = :adjust_vset, qlimit_max_outer = 20)
  bus = geNetBusIdx(net = net, busName = "STATION1")
  return erg == 0 && getNodeType(net.nodeVec[bus]) == Sparlectra.PV && net.nodeVec[bus]._vm_pu < 1.025
end

function test_q_limit_adjust_vset_no_controller_switches()::Bool
  net = createTest3BusNet(cooldown = 0, hyst_pu = 0.0, qlim_min = -15.0, qlim_max = 15.0)
  _, erg = runpf!(net, 30, 1e-8, 0; method = :rectangular, qlimit_mode = :adjust_vset)
  bus = geNetBusIdx(net = net, busName = "STATION1")
  return erg == 0 && getNodeType(net.nodeVec[bus]) == Sparlectra.PQ
end

function test_q_limit_adjust_vset_multiple_controllers_error()::Bool
  net = createTest5BusNet(cooldown = 0, hyst_pu = 0.0, qlim_min = -15.0, qlim_max = 15.0, mul_gens = true)
  _assign_vset_controller!(net, "B3"; vstep_pu = 0.005, tap_steps_down = 2, tap_steps_up = 2, all_generators = true)
  try
    runpf!(net, 20, 1e-8, 0; method = :rectangular, qlimit_mode = :adjust_vset)
  catch
    return true
  end
  return false
end

function test_q_limit_adjust_vset_invalid_config_error()::Bool
  net = createTest3BusNet(cooldown = 0, hyst_pu = 0.0, qlim_min = -20.0, qlim_max = 20.0)
  _assign_vset_controller!(net, "STATION1"; vstep_pu = 0.005, tap_steps_down = nothing, tap_steps_up = nothing)
  try
    runpf!(net, 20, 1e-8, 0; method = :rectangular, qlimit_mode = :adjust_vset)
  catch
    return true
  end
  return false
end

function test_q_limit_adjust_vset_step_exhaustion_fallback()::Bool
  net = createTest3BusNet(cooldown = 0, hyst_pu = 0.0, qlim_min = -15.0, qlim_max = 15.0)
  _assign_vset_controller!(net, "STATION1"; vstep_pu = 0.005, tap_steps_down = 0, tap_steps_up = 0)
  _, erg = runpf!(net, 30, 1e-8, 0; method = :rectangular, qlimit_mode = :adjust_vset, qlimit_max_outer = 1)
  bus = geNetBusIdx(net = net, busName = "STATION1")
  return erg == 0 && getNodeType(net.nodeVec[bus]) == Sparlectra.PQ
end

function test_q_limit_adjust_vset_multigen_single_controller()::Bool
  net = createTest5BusNet(cooldown = 0, hyst_pu = 0.0, qlim_min = -15.0, qlim_max = 15.0, mul_gens = true)
  _assign_vset_controller!(net, "B3"; vstep_pu = 0.005, tap_steps_down = 1, tap_steps_up = 1, all_generators = false)
  _, erg = runpf!(net, 40, 1e-8, 0; method = :rectangular, qlimit_mode = :adjust_vset)
  return erg == 0
end

function test_q_limit_default_behavior_unchanged()::Bool
  net_default = createTest3BusNet(cooldown = 0, hyst_pu = 0.0, qlim_min = -15.0, qlim_max = 15.0)
  net_explicit = createTest3BusNet(cooldown = 0, hyst_pu = 0.0, qlim_min = -15.0, qlim_max = 15.0)
  _, erg_default = runpf!(net_default, 30, 1e-8, 0; method = :rectangular)
  _, erg_explicit = runpf!(net_explicit, 30, 1e-8, 0; method = :rectangular, qlimit_mode = :switch_to_pq)
  b = geNetBusIdx(net = net_default, busName = "STATION1")
  same_type = getNodeType(net_default.nodeVec[b]) == getNodeType(net_explicit.nodeVec[b])
  same_vm = isapprox(net_default.nodeVec[b]._vm_pu, net_explicit.nodeVec[b]._vm_pu; atol = 1e-8, rtol = 0.0)
  return erg_default == 0 && erg_explicit == 0 && same_type && same_vm
end


function test_q_limit_start_iter_delays_switching()::Bool
  net = createTest3BusNet(cooldown = 0, hyst_pu = 0.0, qlim_min = -15.0, qlim_max = 15.0)
  _, erg = runpf!(net, 3, 1e-12, 0; method = :rectangular, qlimit_start_iter = 10)
  bus = geNetBusIdx(net = net, busName = "STATION1")
  return erg == 1 && getNodeType(net.nodeVec[bus]) == Sparlectra.PV && isempty(net.qLimitLog)
end

function test_q_limit_auto_q_delta_accepts_switching()::Bool
  net = createTest3BusNet(cooldown = 0, hyst_pu = 0.0, qlim_min = -15.0, qlim_max = 15.0)
  _, erg = runpf!(net, 30, 1e-8, 0; method = :rectangular, qlimit_start_mode = :auto_q_delta, qlimit_auto_q_delta_pu = Inf)
  bus = geNetBusIdx(net = net, busName = "STATION1")
  return erg == 0 && getNodeType(net.nodeVec[bus]) == Sparlectra.PQ && !isempty(net.qLimitLog)
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
  second_switch = changed && bus_types[1] == :PQ && length(net.qLimitLog) == 2

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
  second_reenable_blocked = !changed && !reenabled && bus_types[1] == :PQ && length(net.qLimitLog) == 2

  return delayed && first_switch && first_reenable && second_switch && second_reenable_blocked
end

function test_q_limit_violation_guard_bypasses_hysteresis_for_strong_overshoot()::Bool
  net = Net(name = "q_limit_violation_guard", baseMVA = 100.0)
  nb = 1
  qmin_pu = [-1.0]
  qmax_pu = [1.0]
  bus_types = [:PV]
  pv_orig_mask = trues(nb)

  changed, reenabled = Sparlectra.active_set_q_limits!(
    net,
    2,
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
    qlimit_guard_violation_mode = :lock_pq,
    qlimit_guard_violation_threshold_pu = 1e-4,
  )

  return changed && !reenabled && bus_types[1] == :PQ && length(net.qLimitLog) == 1 && net.qLimitLog[1].side == :max
end

function test_solve_linear_singular_sparse_qr_fallback()::Bool
  A = sparse([1.0 1.0; 2.0 2.0])
  b = [1.0, 2.0]
  x = Sparlectra.solve_linear(A, b)
  return all(isfinite, x) && isapprox(A * x, b; atol = 1e-10, rtol = 0.0)
end
function test_solve_linear_regularized_sparse_fallback()::Bool
  A = spzeros(Float64, 3, 3)
  b = [1.0, -2.0, 3.0]
  x = Sparlectra._regularized_normal_solve(A, b)
  return length(x) == 3 && all(isfinite, x) && norm(x) <= 1e-12
end

function test_rectangular_singular_step_returns_nonconvergence()::Bool
  net = Net(name = "singular_rectangular_step", baseMVA = 100.0)
  addBus!(net = net, busName = "Slack", vn_kV = 110.0)
  addBus!(net = net, busName = "Load", vn_kV = 110.0)
  addProsumer!(net = net, busName = "Slack", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "Slack")
  addProsumer!(net = net, busName = "Load", type = "ENERGYCONSUMER", p = 10.0, q = 5.0)
  _, erg = runpf!(net, 5, 1e-8, 0; method = :rectangular, opt_sparse = true)
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

  return ok_types && pv_regulation_ok && genpq_regulation_ok
end

function test_multiple_slack_prosumers_same_bus_supported()::Bool
  net = Net(name = "multi_slack_same_bus", baseMVA = 100.0)
  addBus!(net = net, busName = "SlackBus", vn_kV = 110.0)
  addBus!(net = net, busName = "LoadBus", vn_kV = 110.0)

  addProsumer!(net = net, busName = "SlackBus", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "SlackBus")
  # Second in-service generator at the same slack bus (MATPOWER-like multi-gen slack bus case).
  addProsumer!(net = net, busName = "SlackBus", type = "GENERATOR", p = 15.0, q = 2.0, referencePri = "SlackBus")
  addProsumer!(net = net, busName = "LoadBus", type = "ENERGYCONSUMER", p = 10.0, q = 3.0)

  slack_bus = geNetBusIdx(net = net, busName = "SlackBus")
  load_bus = geNetBusIdx(net = net, busName = "LoadBus")
  return getEffectiveBusType(net, slack_bus) == Sparlectra.Slack &&
         getEffectiveBusType(net, load_bus) == Sparlectra.PQ &&
         count(isSlack, getBusProsumers(net, slack_bus)) == 2
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

  pv_bus = geNetBusIdx(net = net, busName = "PVBus")
  node = net.nodeVec[pv_bus]
  qmin_pu, qmax_pu = getQLimits_pu(net)

  p_target_ok = isapprox(node._pƩGen, 27.5; atol = 1e-9, rtol = 0.0)
  q_target_ok = isapprox(node._qƩGen, 4.0; atol = 1e-9, rtol = 0.0)
  q_limit_ok = isapprox(qmin_pu[pv_bus], -0.4; atol = 1e-9, rtol = 0.0) &&
               isapprox(qmax_pu[pv_bus], 0.55; atol = 1e-9, rtol = 0.0)
  bus_type_ok = getEffectiveBusType(net = net, busName = "PVBus") == Sparlectra.PV

  return p_target_ok && q_target_ok && q_limit_ok && bus_type_ok
end

function test_addprosumer_auto_regulated_flags()::Bool
  net = Net(name = "addprosumer_auto_regulated_flags", baseMVA = 100.0)
  addBus!(net = net, busName = "Slack", vn_kV = 110.0)
  addBus!(net = net, busName = "PVVm", vn_kV = 110.0)
  addBus!(net = net, busName = "PVTap", vn_kV = 110.0)
  addBus!(net = net, busName = "PQ", vn_kV = 110.0)

  addProsumer!(net = net, busName = "Slack", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "Slack")
  addProsumer!(net = net, busName = "PVVm", type = "SYNCHRONOUSMACHINE", p = 10.0, vm_pu = 1.01)
  addProsumer!(net = net, busName = "PVTap", type = "SYNCHRONOUSMACHINE", p = 8.0, vstep_pu = 0.005, tap_steps_down = 2, tap_steps_up = 2)
  addProsumer!(net = net, busName = "PQ", type = "SYNCHRONOUSMACHINE", p = 6.0, q = 1.5)

  ps_vm = only(getBusProsumers(net, geNetBusIdx(net = net, busName = "PVVm")))
  ps_tap = only(getBusProsumers(net, geNetBusIdx(net = net, busName = "PVTap")))
  ps_pq = only(getBusProsumers(net, geNetBusIdx(net = net, busName = "PQ")))

  flags_ok = ps_vm.isRegulated && ps_tap.isRegulated && !ps_pq.isRegulated
  bus_types_ok =
    getEffectiveBusType(net = net, busName = "PVVm") == Sparlectra.PV &&
    getEffectiveBusType(net = net, busName = "PVTap") == Sparlectra.PV &&
    getEffectiveBusType(net = net, busName = "PQ") == Sparlectra.PQ

  return flags_ok && bus_types_ok
end

function run_grid_tests()
  @testset "Grid and power-flow regression tests" begin
    @testset "Transformer and network validation" begin
      @test test_2WTPITrafo() == true
      @test test_3WTPITrafo() == true
      @test testNetwork() == true
      @test test_NBI_MDO() == true
      @test testISOBusses() == true
    end

    @testset "MATPOWER import/export" begin
      @test testImportMatpower() == true
      @test test_matpower_import_defaults_no_reenable() == true
      @test test_matpower_import_uses_bus_type_for_regulation() == true
      @test test_matpower_flatstart_uses_generator_voltage_setpoints() == true
      @test test_matpower_reference_override_controls_flatstart() == true
      @test test_prosumer_aggregation_preserves_bus_types_and_injections() == true
      @test test_matpower_vmva_selfcheck_noncontiguous_bus_numbers() == true
      @test test_matpower_vmva_selfcheck_ignores_slack_pq_spec() == true
      @test test_matpower_compare_vmva_wraps_angle_differences() == true
      @test test_matpower_pv_voltage_source_and_compare_modes() == true
      @test test_matpower_compare_vmva_final_pq_after_qlimit_uses_bus_vm() == true
      @test test_matpower_multi_generator_vg_and_fallback() == true
      @test test_mp_inline_vs_manual_shunt() == true
      @test test_bus_shunt_model_modes() == true
      @test test_matpower_read_case_m_postprocessing() == true
    end

    @testset "Power flow scenarios" begin
      @test test_5BusNet(0, 10.0) == true
      @test test_3BusNet(0, 150.0, :rectangular, false, false) == true
      @test test_3BusNet(0, 150.0, :polar_full, false, false) == true
      @test test_acpflow(0; lLine_6a6b = 0.01, damp = 1.0, method = :classic, opt_sparse = true) == true
      @test test_acpflow(0; lLine_6a6b = 0.01, damp = 1.0, method = :rectangular, opt_sparse = true) == true
      @test test_acpflow(0; lLine_6a6b = 0.01, damp = 1.0, method = :rectangular, opt_sparse = false) == true
      @test test_acpflow(0; lLine_6a6b = 0.01, damp = 1.0, method = :polar_full, opt_sparse = true) == true
      @test test_rectangular_autodamp_backtracks_oversized_step() == true
      @test test_rectangular_start_projection_improves_dc_seed() == true
      @test test_q_limit_adjust_vset_success() == true
      @test test_q_limit_adjust_vset_no_controller_switches() == true
      @test test_q_limit_adjust_vset_multiple_controllers_error() == true
      @test test_q_limit_adjust_vset_invalid_config_error() == true
      @test test_q_limit_adjust_vset_step_exhaustion_fallback() == true
      @test test_q_limit_adjust_vset_multigen_single_controller() == true
      @test test_q_limit_default_behavior_unchanged() == true
      @test test_q_limit_start_iter_delays_switching() == true
      @test test_q_limit_auto_q_delta_accepts_switching() == true
      @test test_q_limit_hysteresis_delays_small_pv_to_pq_overshoot() == true
      @test test_q_limit_violation_guard_bypasses_hysteresis_for_strong_overshoot() == true
      @test test_solve_linear_singular_sparse_qr_fallback() == true
      @test test_rectangular_singular_step_returns_nonconvergence() == true
      @test test_bus_type_resolution_from_prosumers() == true
      @test test_multiple_slack_prosumers_same_bus_supported() == true
      @test test_regulated_generator_bus_targets_include_unregulated_generators() == true
      @test test_addprosumer_auto_regulated_flags() == true
    end

    @testset "Link behaviour and reporting" begin
      @test test_link_kcl_simple() == true
      @test test_link_component_type() == true
      @test test_link_rejects_slack_connection() == true
      @test test_link_rejects_mixed_bus_types() == true
      @test test_link_bus_merge_pf() == true
      @test test_link_bus_merge_pf_default_method() == true
      @test test_rectangular_merge_fallback_suppresses_polar_deprecation_warning() == true
      @test test_link_closed_keeps_shunt_reporting_on_original_bus() == true
      @test test_link_kcl_ring_allocation() == true
      @test test_link_kcl_ring_allocation_with_shunt() == true
      @test test_link_ring_pf_stability() == true
      @test test_acpflow_report_object() == true
      @test test_report_uses_user_bus_names_and_pf_node_count() == true
    end
  end
end
