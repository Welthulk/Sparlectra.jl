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

"""
Function to perform AC power flow analysis.

Parameters:
- max_ite: Int, the maximum number of iterations for the power flow algorithm (default: 30).
- tol: Float64, tolerance for convergence criterion (default: 1e-6).
- casefile: String, the name of the case file to load.
- path: Union{Nothing,String}, the path to the case file (default: nothing).
- verbose: Int, verbosity level for output (default: 0).
- printResultToFile: Bool, flag to print results to a file (default: false).
- printResultAnyCase: Bool, flag to print results even if the power flow fails (default: false).

Returns:
- Net, the network object.
"""
function run_acpflow(;
  max_ite::Int = 30,
  tol::Float64 = 1e-6,
  casefile::String,
  path::Union{Nothing,String} = nothing,
  verbose::Int = 0,  
  printResultToFile::Bool = false,
  printResultAnyCase::Bool = false,
  opt_fd::Bool = false,
  opt_sparse::Bool = false,
  method::Symbol = :polar_full,
  opt_flatstart::Bool = false,
  show_results::Bool = true,
  cooldown_iters::Int = 0,
  q_hyst_pu::Float64 = 0.0,
  pv_table_rows::Int = 30,
  check_q_limit_signs::Bool = false,
  autocorrect_q_limit_signs::Bool = false,
  validate_limits_after_pf::Bool = false,
  q_limit_violation_headroom::Float64 = 0.20,
  lock_pv_to_pq_buses::AbstractVector{Int} = Int[],
)::Net  
  ext = lowercase(splitext(casefile)[2])
  myNet = nothing              # Initialize myNet variable
  in_path = nothing
  out_path = nothing
  if ext in (".m", ".jl")
    if path === nothing
      in_path  = joinpath(pwd(), "data", "mpower", strip(casefile))
      out_path = joinpath(pwd(), "data", "mpower")
    else
      in_path  = joinpath(path, strip(casefile))
      out_path = joinpath(path)
    end

    if !isfile(in_path)
      error("File $(in_path) not found")
    end

    myNet = createNetFromMatPowerFile(filename = in_path, log = (verbose > 0), flatstart = opt_flatstart, cooldown = cooldown_iters, q_hyst_pu = q_hyst_pu)
    
    if verbose > 1
      # --- DEBUG START ---
      @info "DEBUG Full net after import:"
      @printf "DEBUG Full net after import:\n"
      showNet(myNet, verbose = true)

      Y = createYBUS(net=myNet, sparse=false, printYBUS=false)  # dense for inspection      
      V0, slack = initialVrect(myNet; flatstart=myNet.flatstart)
      S = buildComplexSVec(myNet)

      @printf "DEBUG Yabs_max: %.6e\n" maximum(abs.(Y))
      @printf "DEBUG Ydiag_max: %.6e\n" maximum(abs.(diag(Y)))
      @printf "DEBUG Ydiag_imag_max: %.6e\n" maximum(abs.(imag.(diag(Y))))

      @printf "DEBUG slack bus index: %d\n" slack
      @printf "DEBUG initial V0: %s\n" string(V0)
      @printf "DEBUG case: name=%s baseMVA=%.1f flatstart=%s slack=%d\n" myNet.name myNet.baseMVA myNet.flatstart slack
      @printf "DEBUG sums (pu): sumP=%.6f sumQ=%.6f\n" sum(real.(S)) sum(imag.(S))
      @printf "DEBUG V0: Vslack=%s Vmin=%.6f Vmax=%.6f\n" string(V0[slack]) minimum(abs.(V0)) maximum(abs.(V0))
      @printf "DEBUG Y: Ydiag_min=%.6e Ydiag_max=%.6e Yabs_max=%.6e\n" minimum(abs.(diag(Y))) maximum(abs.(diag(Y))) maximum(abs.(Y))
      # --- DEBUG END ---
    end

  else
    error("File extension $(ext) not supported! Use .m or .jl")
  end

  # Resolve PV->PQ lock list (prefer original bus IDs; fallback to internal indices).
  lock_pv_to_pq_buses_resolved = Int[]
  if !isempty(lock_pv_to_pq_buses)
    orig_to_net = Dict{Int,Int}()
    for (net_idx, orig_idx) in myNet.busOrigIdxDict
      orig_to_net[orig_idx] = net_idx
    end
    for bus in lock_pv_to_pq_buses
      if haskey(orig_to_net, bus)
        push!(lock_pv_to_pq_buses_resolved, orig_to_net[bus])
      elseif 1 <= bus <= length(myNet.nodeVec)
        push!(lock_pv_to_pq_buses_resolved, bus)
      elseif verbose > 0
        @warn "lock_pv_to_pq_buses entry not found (ignored): $bus"
      end
    end
    unique!(lock_pv_to_pq_buses_resolved)
    sort!(lock_pv_to_pq_buses_resolved)
  end

  # Run power flow
  ite = 0
  etime = @elapsed begin
    ite, erg = runpf!(myNet, max_ite, tol, verbose; opt_fd = opt_fd, opt_sparse = opt_sparse, method = method, opt_flatstart = opt_flatstart, pv_table_rows = pv_table_rows, validate_limits_after_pf = validate_limits_after_pf, q_limit_violation_headroom = q_limit_violation_headroom, lock_pv_to_pq_buses = lock_pv_to_pq_buses_resolved)
  end
  
  if erg == 0 || printResultAnyCase
    # Calculate network losses and print results
    calcNetLosses!(myNet)
    calcLinkFlowsKCL!(myNet)
    jpath = printResultToFile ? out_path : ""
    if show_results || printResultAnyCase
      printACPFlowResults(myNet, etime, ite, tol, printResultToFile, jpath; converged = (erg == 0), solver = method)
    end
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
- max_ite: Int, the maximum number of iterations for the power flow algorithm (default: 30).
- tol: Float64, tolerance for convergence criterion (default: 1e-6).
- verbose: Int, verbosity level for output (default: 0).
- printResultToFile: Bool, flag to print results to a file (default: false).
- printResultAnyCase: Bool, flag to print results even if the power flow fails (default: false).
"""
function run_net_acpflow(; net::Net, max_ite::Int = 30, tol::Float64 = 1e-6, verbose::Int = 0, printResultToFile::Bool = false, printResultAnyCase::Bool = false, opt_fd::Bool = false, opt_sparse::Bool = false, method::Symbol = :polar_full, show_results::Bool = true, lock_pv_to_pq_buses::AbstractVector{Int} = Int[])

  # Run power flow
  ite = 0
  etime = @elapsed begin
    ite, erg = runpf!(net, max_ite, tol, verbose; opt_fd = opt_fd, opt_sparse = opt_sparse, method = method, lock_pv_to_pq_buses = lock_pv_to_pq_buses)
  end

  if erg == 0 || printResultAnyCase
    # Calculate network losses and print results
    calcNetLosses!(net)
    calcLinkFlowsKCL!(net)
    jpath = ""
    if show_results
      printACPFlowResults(net, etime, ite, tol, printResultToFile, jpath; converged = (erg == 0), solver = method)
    end
  elseif erg == 1
    println("Newton-Raphson did not converge")
  else
    @error "Errors during calculation of Newton-Raphson"
  end

  return (ite, erg, etime)
end
