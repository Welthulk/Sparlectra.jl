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

# file: test/test_experments.jl
# for expermimental work and testing of new features
using Test
using Sparlectra
using Printf
using Logging

include("testgrid.jl")
global_logger(ConsoleLogger(stderr, Logging.Debug))

function importCaseFile()
  case = "case5.m"
  run_acpflow(casefile = case, opt_sparse = true, method = :polar_full)
  run_acpflow(casefile = case, opt_fd = true, method = :polar_full)
  run_acpflow(casefile = case, opt_fd = true, opt_sparse = true, method = :polar_full)
  run_acpflow(casefile = case, opt_sparse = true, opt_fd = false, method = :rectangular)
  run_acpflow(casefile = case, opt_sparse = true, opt_fd = false, method = :classic)
end

function cross_check()
  net = createCIGRE()
  open("network_dump.txt", "w") do io
    println(io, "=== TRANSFORMERS ===")
    show(io, "text/plain", net.trafos)

    println(io, "\n\n=== AC LINES ===")
    show(io, "text/plain", net.linesAC)
  end

  # 1) Referenz: polar_full (should convergen)
  iters_polar, status_polar = runpf!(net, 10, 1e-8, 2; method = :polar_full, opt_fd = false)
  println("polar_full: iters = $iters_polar, status = $status_polar")
  println("---------------------------------------------------")

  ## 2) FD-Jacobian
  #net = createCIGRE()
  iters_fd, status_fd = runpf!(net, 10, 1e-8, 1; method = :rectangular, opt_fd = true)
  println("rectangular FD: iters = $iters_fd, status = $status_fd")
  println("---------------------------------------------------")
  # 3) analytic Jacobian
  net = createCIGRE()
  iters_ana, status_ana = runpf!(net, 10, 1e-8, 1; method = :rectangular, opt_fd = false)
  println("rectangular analytic: iters = $iters_ana, status = $status_ana")
  println("---------------------------------------------------")
  # 4) classic Jacobian
  net = createCIGRE()
  iters_cl, status_cl = runpf!(net, 10, 1e-8, 1; method = :classic, opt_fd = false)
  println("jacobian classic: iters = $iters_cl, status = $status_cl")
  println("---------------------------------------------------")
end

function PGMTestBusNet()::Net
  # Simple 2-bus network
  #   Bus1            Bus2
  #  S->|---------------|
  #     # 
  #  Slack           

  Sbase_MVA = 1.0
  netName = "PGMBus"
  r = 0.01
  x = 0.1
  c_nf_per_km = 0.0
  tanδ = 0.0
  length_km = 20.0
  method = :polar_full
  @debug "Creating $netName test network for PGM experiments"

  Bus2Net = Net(name = netName, baseMVA = Sbase_MVA, cooldown_iters = 0, q_hyst_pu = 0.1)

  addBus!(net = Bus2Net, busName = "B1", busType = "Slack", vn_kV = 110.0)
  addBus!(net = Bus2Net, busName = "B2", busType = "PQ", vn_kV = 110.0)

  addACLine!(net = Bus2Net, fromBus = "B1", toBus = "B2", length = length_km, r = r, x = x, c_nf_per_km = c_nf_per_km, tanδ = tanδ)

  addProsumer!(net = Bus2Net, busName = "B1", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "B1")
  addShunt!(net = Bus2Net, busName = "B1", pShunt = 0.0, qShunt = 5.0)

  tol = 1e-9
  maxIte = 50
  print_results = true
  result = true
  verbose = 2
  pv_names = ["B3"]

  etim = 0.0
  etim = @elapsed begin
    ite, erg = runpf!(Bus2Net, maxIte, tol, verbose, method = :polar_full, opt_fd = false, opt_sparse = false)
    if erg != 0
      @info "Full-system power flow did not converge"
      result = false
    end
  end

  #

  #hit = pv_hit_q_limit(net, pv_names)

  calcNetLosses!(Bus2Net)
  distributeBusResults!(Bus2Net)
  if print_results
    printACPFlowResults(Bus2Net, etim, ite, tol; converged = result, solver = method)
    println("---------------------------------------------------")
    println(Bus2Net)
    printProsumerResults(Bus2Net)
    #printQLimitLog(Bus2Net; sort_by = :bus)
  end

  return Bus2Net
end

#test_acpflow(1;lLine_6a6b = 0.01, damp = 1.0, method = :rectangular, opt_sparse = true)
#test_5BusNet(1, 500.0, :rectangular, false, false)
#test_5BusNet(1, 500.0, :polar_full, false, false)
#test_5BusNet(1, 500.0, :classic, false, false)
#test_3BusNet(2, 50.0, :rectangular, false, false)
#test_3BusNet(2, 30.0, :rectangular, true, true)
#test_3BusNet(2, 30.0, :rectangular, true, false)
#test_3BusNet(2, 30.0, :rectangular, false, false)
#test_3BusNet(2, 30.0, :classic, false, false)
#test_3BusNet(2, 30.0, :polar_full, false, false)
#PGMTestBusNet()
test_2WTPITrafo()
test_3WTPITrafo(2; method = :rectangular, opt_fd = false, opt_sparse = true)
