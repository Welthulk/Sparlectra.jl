"""
Function to perform AC power flow analysis.

Parameters:
- max_ite: Int, the maximum number of iterations for the power flow algorithm (default: 10).
- tol: Float64, tolerance for convergence criterion (default: 1e-6).
- casefile: String, the name of the case file to load.
- path: Union{Nothing,String}, the path to the case file (default: nothing).
- verbose: Int, verbosity level for output (default: 0).
- printResultToFile: Bool, flag to print results to a file (default: false).
- printResultAnyCase: Bool, flag to print results even if the power flow fails (default: false).

Returns:
- Net, the network object.
"""
function run_acpflow(; max_ite::Int = 10, tol::Float64 = 1e-6, casefile::String, path::Union{Nothing,String} = nothing, verbose::Int = 0, printResultToFile::Bool = false, printResultAnyCase::Bool = false, opt_fd::Bool = false, opt_sparse::Bool = false, method::Symbol = :polar_full)::Net
  ext = splitext(casefile)[2]  # Get the file extension
  myNet = nothing              # Initialize myNet variable
  in_path = nothing
  out_path = nothing
  if ext == ".m"
    # Read network data from Matpower .m file    
    if path === nothing
      in_path = joinpath(pwd(), "data", "mpower", strip(casefile))
      out_path = joinpath(pwd(), "data", "mpower")
    else
      in_path = joinpath(path, strip(casefile))
      out_path = joinpath(path)
    end
    if !isfile(in_path)
      println("Error: File $(in_path) not found")
      return
    end
    myNet = createNetFromMatPowerFile(in_path, (verbose > 0))
  else
    println("Error: File extension $(ext) not supported!")
    return
  end

  # Run power flow
  ite = 0
  etime = @elapsed begin
    ite, erg = runpf!(myNet, max_ite, tol, verbose; opt_fd = opt_fd, opt_sparse = opt_sparse, method = method)
  end

  if erg == 0 || printResultAnyCase
    # Calculate network losses and print results
    calcNetLosses!(myNet)
    jpath = printResultToFile ? out_path : ""
    printACPFlowResults(myNet, etime, ite, tol, printResultToFile, jpath; converged = (erg == 0), solver = method)
  elseif erg == 1
    println("Newton-Raphson did not converge")
  else
    @error "Errors during calculation of Newton-Raphson"
  end

  return myNet
end

"""
Function to perform AC power flow analysis.

Parameters:
- net: Net, the network object.
- max_ite: Int, the maximum number of iterations for the power flow algorithm (default: 10).
- tol: Float64, tolerance for convergence criterion (default: 1e-6).
- verbose: Int, verbosity level for output (default: 0).
- printResultToFile: Bool, flag to print results to a file (default: false).
- printResultAnyCase: Bool, flag to print results even if the power flow fails (default: false).
"""
function run_net_acpflow(; net::Net, max_ite::Int = 10, tol::Float64 = 1e-6, verbose::Int = 0, printResultToFile::Bool = false, printResultAnyCase::Bool = false, opt_fd::Bool = false, opt_sparse::Bool = false, method::Symbol = :polar_full)

  # Run power flow
  ite = 0
  etime = @elapsed begin
    ite, erg = runpf!(net, max_ite, tol, verbose; opt_fd = opt_fd, opt_sparse = opt_sparse, method = method)
  end

  if erg == 0 || printResultAnyCase
    # Calculate network losses and print results
    calcNetLosses!(net)
    jpath = printResultToFile ? out_path : ""
    printACPFlowResults(net, etime, ite, tol, printResultToFile, jpath; converged = (erg == 0), solver = method)
  elseif erg == 1
    println("Newton-Raphson did not converge")
  else
    @error "Errors during calculation of Newton-Raphson"
  end
end
