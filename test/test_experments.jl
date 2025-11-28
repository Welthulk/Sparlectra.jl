# file: test/test_expermints.jl
# for expermimental work and testing of new features
using Test
using Sparlectra
using Printf
using Logging

include("testgrid.jl")
global_logger(ConsoleLogger(stderr, Logging.Info))
#=
case = "case5.m"   # oder irgendein anderes .m im Pfad

# Test 1 – Standard Newton (polar), aber mit Sparse:
run_acpflow(casefile = case, opt_sparse = true, method = :polar_full)

# Test 2 – FD-Jacobian aktivieren:
run_acpflow(casefile = case, opt_fd = true, method = :polar_full)

# Test 3 – Polar-NR: FD + Sparse zusammen:
run_acpflow(casefile = case, opt_fd = true, opt_sparse = true, method = :polar_full)

# Test 4 – Vergleich: Rectangular-NR (als Referenz)
run_acpflow(casefile = case, opt_sparse = true, opt_fd = false, method = :rectangular)

# Test 4 – Vergleich: Classic-NR (als Referenz)
run_acpflow(casefile = case, opt_sparse = true, opt_fd = false, method = :classic)
=#


function test_5BusNet(verbose::Int = 0, qlim::Float64 = 20.0, methode::Symbol = :rectangular, opt_fd::Bool = false, opt_sparse::Bool = false)
    net = createTest5BusNet(cooldown = 2, hyst_pu = 0.000, qlim_min = -qlim, qlim_max = qlim)
    tol = 1e-9
    maxIte = 50
    print_results = (verbose > 0)
    result = true
    
    pv_names = ["B2"]        
    etim = 0.0
    etim = @elapsed begin        
        ite, erg = runpf!(net, maxIte, tol, verbose, method = methode, opt_fd = opt_fd, opt_sparse = opt_sparse)
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
        printACPFlowResults(net, etim, ite, tol)
        printQLimitLog(net; sort_by=:bus)
    end    

    return hit==true
end
opt_fd = true
opt_sparse = true


test_5BusNet(1,5.0, :polar_full, opt_fd, opt_sparse)
test_5BusNet(1,5.0, :rectangular, opt_fd, opt_sparse)
test_5BusNet(1,5.0, :classic, opt_fd, opt_sparse)