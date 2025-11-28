# file: test/test_expermints.jl
# for expermimental work and testing of new features
using Test
using Sparlectra
using Printf
using Logging

include("testgrid.jl")
global_logger(ConsoleLogger(stderr, Logging.Info))


function test_5BusNet(verbose::Int = 0, qlim::Float64 = 20.0, methode::Symbol = :rectangular    )
    net = createTest5BusNet(cooldown = 2, hyst_pu = 0.000, qlim_min = -qlim, qlim_max = qlim)
    tol = 1e-9
    maxIte = 50
    print_results = (verbose > 0)
    result = true
    
    pv_names = ["B2"]        
    etim = 0.0
    etim = @elapsed begin        
        ite, erg = runpf!(net, maxIte, tol, verbose, method = methode)
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

test_5BusNet(1,5.0, :polar_full)
test_5BusNet(1,5.0, :rectangular)
test_5BusNet(1,5.0, :classic)