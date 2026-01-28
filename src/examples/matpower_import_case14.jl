using Sparlectra
#case = "case14.m"
case = "case14.m"
flatstart = true
filename = joinpath(@__DIR__, "..", "..", "data", "mpower", case)
log = true
#net = createNetFromMatPowerFile(filename = filename, log = false, flatstart = true)

#showNet(net, verbose = true) 

#run_acpflow(casefile = case, opt_sparse = true, opt_fd = false, method = :classic, opt_flatstart = flatstart)

#status, iters = runpf!(net, 25)

#println("status = ", status, "   iters = ", iters)
#println("converged = ", status == 1)

#mynet = run_acpflow(casefile = case, opt_sparse = true, method = :polar_full, opt_flatstart = flatstart)
#printQLimitLog(mynet; sort_by = :bus)
#run_acpflow(casefile = case, opt_fd = true, method = :polar_full, opt_flatstart = flatstart)
run_acpflow(casefile = case, opt_fd = true, opt_sparse = true, method = :polar_full, opt_flatstart = flatstart)
#run_acpflow(casefile = case, opt_sparse = true, opt_fd = false, method = :rectangular, opt_flatstart = flatstart)
#run_acpflow(casefile = case, opt_sparse = true, opt_fd = false, method = :classic, opt_flatstart = flatstart)



