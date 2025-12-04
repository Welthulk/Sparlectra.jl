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

function test_5BusNet(verbose::Int = 0, qlim::Float64 = 20.0, opt_fd::Bool = false, opt_sparse::Bool = true)
  net = createTest5BusNet(cooldown = 2, hyst_pu = 0.000, qlim_min = -qlim, qlim_max = qlim)
  tol = 1e-9
  maxIte = 50
  print_results = (verbose > 0)
  result = true

  pv_names = ["B3"]
  etim = 0.0
  etim = @elapsed begin
    ite, erg = runpf!(net, maxIte, tol, verbose; method = :rectangular, opt_fd = opt_fd, opt_sparse = opt_sparse)
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
    distribute_all_bus_results!(net)
    printACPFlowResults(net, etim, ite, tol)
    printProsumerResults(net)
    printQLimitLog(net; sort_by = :bus)
  end

  return hit==true
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

#test_5BusNet(1, 20.0 )
#test_5BusNet(1, 5.0, :polar_full, opt_fd, opt_sparse)
#test_5BusNet(1, 5.0, :classic, opt_fd, opt_sparse)
#test_acpflow(1;lLine_6a6b = 0.01, damp = 1.0, method = :rectangular, opt_sparse = true)
#test_5BusNet(1, 25.0, false, false)
