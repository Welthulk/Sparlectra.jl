using Sparlectra
using Printf
using Logging

include("testgrid.jl")
global_logger(ConsoleLogger(stderr, Logging.Info))

#--- Test helper:  ---
function test_2BusNet(verbose::Int = 0, qlim::Float64 = 20.0)
    net = createTest2BusNet(cooldown = 2, hyst_pu = 0.000, qlim_min = -qlim, qlim_max = qlim)
    tol = 1e-9
    maxIte = 50
    print_results = (verbose > 0)
    result = true
    
    pv_names = ["B2"]    
    
    etim = 0.0
    etim = @elapsed begin
        ite, erg = runpf_full!(net, maxIte, tol, verbose)
        if erg != 0
            @info "Full-system power flow did not converge"
            result = false
        end
    end  
          
    hit = pv_hit_q_limit(net, pv_names)  

    if print_results
        calcNetLosses!(net)
        printACPFlowResults(net, etim, ite, tol)
        printQLimitLog(net; sort_by=:bus)
    end
    

    return hit==true
end

function test_5BusNet(verbose::Int = 0, qlim::Float64 = 20.0)
    net = createTest5BusNet(cooldown = 2, hyst_pu = 0.000, qlim_min = -qlim, qlim_max = qlim)
    tol = 1e-9
    maxIte = 50
    print_results = (verbose > 0)
    result = true
    
    pv_names = ["B2"]    
    
    etim = 0.0
    etim = @elapsed begin
        ite, erg = runpf_full!(net, maxIte, tol, verbose)
        if erg != 0
            @info "Full-system power flow did not converge"
            result = false
        end
    end  
          
    hit = pv_hit_q_limit(net, pv_names)  

    if print_results
        calcNetLosses!(net)
        printACPFlowResults(net, etim, ite, tol)
        printQLimitLog(net; sort_by=:bus)
    end
    

    return hit==true
end

function testSmallGridPF()
    file = "SmallGrid.m"
    path="C:/Users/scud/.julia/dev/Sparlectra/data/mpower"

    printResultToFile = false
    tol = 1e-6
    ite = 10
    verbose = 0       # 0: no output, 1: iteration norm, 2: + Y-Bus, 3: + Jacobian, 4: + Power Flow

    net = run_acpflow(max_ite= ite,tol = tol, path=path, casefile=file, verbose=verbose, printResultToFile = printResultToFile)

end

#testSmallGridPF()
#test_2BusNet(3,15.0)
test_5BusNet(1,15.0)