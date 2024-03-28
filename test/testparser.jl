using Test                   # Import Test module for unit testing
using Sparlectra             # Import Sparlectra module for power flow analysis
using Sparlectra.ResDataTypes  # Import data types from Sparlectra module
using Sparlectra.SparlectraNet # Import network functions from Sparlectra module
using BenchmarkTools         # Import BenchmarkTools module for benchmarking
using Logging                # Import Logging module for logging

global_logger(ConsoleLogger(stderr, Logging.Warn))  # Set global logger to log to stderr with INFO level

function acpflow(casefile::String, iterations::Int=10, verbose::Int=0, mdo::Bool=false, tol::Float64=1e-6, base_MVA::Float64=0.0, printResultToFile::Bool=false, printResultAnyCase::Bool=false)
    """
    Function to perform AC power flow analysis.

    Parameters:
    - casefile: String, the file path of the network data file.
    - iterations: Int, the maximum number of iterations for the power flow algorithm (default: 10).
    - verbose: Int, verbosity level for output (default: 0).
    - mdo: Bool, flag to enable Modo solution (default: false).
    - tol: Float64, tolerance for convergence criterion (default: 1e-6).
    - base_MVA: Float64, base MVA for the network (default: 0.0).
    - printResultToFile: Bool, flag to print results to a file (default: false).
    - printResultAnyCase: Bool, flag to print results even if the power flow fails (default: false).
    """

    ext = splitext(casefile)[2]  # Get the file extension
    myNet = nothing              # Initialize myNet variable

    if ext == ".m"
        # Read network data from Matpower .m file
        path = "mpower"
        jpath = joinpath(pwd(), "data", path, strip(casefile))
        if !isfile(jpath)
            @error "File $(jpath) not found"
            return
        end
        myNet = createNetFromMatPowerFile(jpath, base_MVA, (verbose >= 1), mdo)
    else
        @error "File extension $(ext) not supported"
        return
    end

    # Run power flow
    ite = 0
    etime = @elapsed begin
        ite, erg = runpf!(myNet, iterations, tol, verbose) 
    end

    if erg == 0 || printResultAnyCase
        # Calculate network losses and print results
        calcNetLosses!(myNet)
        jpath = printResultToFile ? joinpath(pwd(), "data", path) : ""
        printACPFlowResults(myNet, etime, ite, tol, printResultToFile, jpath)
    elseif erg == 1
        @warn "Newton-Raphson did not converge"
    else
        @error "Error during calculation of Newton-Raphson"
    end
end

# Define input parameters
file = "SmallGrid.m"
printResultToFile = false
printResultAnyCase = false
mdo = false
tol = 1e-8
ite = 25
verbose = 0       # 0: no output, 1: iteration norm, 2: + Y-Bus, 3: + Jacobian, 4: + Power Flow
base_MVA = 0.0    # Set a value to a non-zero value to override the corresponding value in the case file.

# Call acpflow function with input parameters and measure execution time
@time acpflow(file, ite, verbose, mdo, tol, base_MVA, printResultToFile, printResultAnyCase)
