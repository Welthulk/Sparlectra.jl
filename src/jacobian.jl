# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 09.8.2023
# include-file jacobian.jl

# classical Newton-Raphson method

debug = false
angle_limit = false

function setJacobianAngleLimit(value::Bool)
  global angle_limit = value
end

function setJacobianDebug(value::Bool)
  global debug = value
end

function printVector(power::Vector{Float64}, busVec::Vector{BusData}, first::String = "p", second::String = "q", eng::Bool = false, delta_min::Float64 = 0.0)
  print("[")
  i = 0
  n = length(busVec)

  for (index, bus) in enumerate(busVec)
    if bus.type == PQ
      i += 1
      index_p = i
      i += 1
      index_q = i
      p_val = round(power[index_p], digits = 3)
      q_val = round(power[index_q], digits = 3)
      if eng
        p = @sprintf("%.2e", power[index_p])
        q = @sprintf("%.2e", power[index_q])
      else
        p = round(power[index_p], digits = 3)
        q = round(power[index_q], digits = 3)
      end
      if index != n && abs(p_val) >= delta_min || abs(q_val) >= delta_min
        print("(bus:$(index) [$(index_p),$(index_q)] $(first)=$(p), $(second)=$(q)), ")
      elseif abs(p_val) >= delta_min || abs(q_val) >= delta_min
        print("(bus:$(index) [$(index_p),$(index_q)] $(first)=$(p), $(second)=$(q))")
      end
    elseif bus.type == PV
      i += 1
      index_p = i
      p = round(power[index_p], digits = 2)
      if index != n && abs(p) >= delta_min
        print("(bus:$(index) [$(index_p)] $(first)=$(p)), ")
      elseif abs(p) >= delta_min
        print("(bus:$(index) [$(index_p)] $(first)=$(p))")
      end
    end
  end
  print("]'\n")
end

function getPowerFeeds(busVec::Vector{BusData}, n_pq::Int, n_pv::Int)::Vector{Float64}
  size = n_pq * 2 + n_pv

  power_feeds = zeros(Float64, size)

  if debug
    n = n_pq + n_pv
    println("\nPower Feeds: Elements = $(n), Vector-Size = $(size)")
  end

  vindex = 0 # index of vector of power flows    
  for bus in busVec
    if bus.type == PQ
      vindex += 1
      index_p = vindex
      power_feeds[index_p] = bus.pƩ
      vindex += 1
      index_q = vindex
      power_feeds[index_q] = bus.qƩ
    elseif bus.type == PV
      vindex += 1
      index_p = vindex
      power_feeds[index_p] = bus.pƩ
    end
  end

  if debug
    printVector(power_feeds, busVec)
  end
  return power_feeds
end

#=
Matrix formulation of calculation of node powers.
S = Vdiag * conj(Y*V)

- `Y::AbstractMatrix{ComplexF64}`: Y-Bus matrix.
- `busVec::Vector{BusData}`: Vector of bus data.
- `feeders::Vector{Float64}`: Vector of feeders.
- `n_pq::Int`: Number of PQ buses.
- `n_pv::Int`: Number of PV buses.
- `log::Bool`: Whether to log intermediate results.

# Returns
A vector representing the residuum.
=#
function residuum(Y::AbstractMatrix{ComplexF64}, busVec::Vector{BusData}, feeders::Vector{Float64}, n_pq::Int, n_pv::Int, log::Bool)::Vector{Float64}
  # create complex vector of voltages
  V = [bus.vm_pu * exp(im * bus.va_rad) for bus in busVec]
  # create diagonal matrix of voltages
  Vdiag = Diagonal(V)

  # Power Calculation (Knotenleistung)
  S = Vdiag * conj(Y * V)
  size = n_pq * 2 + n_pv
  Δpq = zeros(Float64, size)

  if log
    println("\nresiduum: Vector-Size = $(size)")
  end

  vindex = 0
  pfIdx = 0
  @inbounds for bus in busVec
    pfIdx += 1
    if bus.type == PQ
      vindex += 1
      index_p = vindex
      Δpq[index_p] = feeders[index_p] - real(S[pfIdx])
      vindex += 1
      index_q = vindex
      Δpq[index_q] = feeders[index_q] - imag(S[pfIdx])
      bus._pRes = real(S[pfIdx])
      bus._qRes = imag(S[pfIdx])
    elseif bus.type == PV
      vindex += 1
      index_p = vindex
      Δpq[index_p] = feeders[index_p] - real(S[pfIdx])
      bus._pRes = real(S[pfIdx])
      bus._qRes = imag(S[pfIdx])
    elseif bus.type == Slack
      bus._pRes = real(S[pfIdx])
      bus._qRes = imag(S[pfIdx])
    end
  end

  if log
    printVector(Δpq, busVec, "Δp", "Δq", false, 0.0)
    norm_delta_p = norm(Δpq)
    if debug
      println("residuum: $(norm_delta_p)")
      println("S: $S")
      println("feeders: $feeders")
    end
  end

  return Δpq
end

#=
purpose: Calculation of Jacobian-Matrix
Date: 16.8.23
 Ns : Number of State Variables, a slack node needs to be removed
 Ns = 2*(Nk-1) - n_pv
 Ns = 2*(n_pq+n_pv +1 -1) -n_pv
 Ns = 2*n_pq + n_pv
 Ns-> size of jacobian matrix

    (2*i-1,2*j-1) | (2*i-1,2*j)  
    ∂p/∂𝜑=H       | V*∂p/∂v=N
    ------------------------------
    (2*i,2*j-1)   | (2*i,2*j)
    ∂q/∂v=J       | V*∂q/∂v=L

Realpower:
Hij = ∂pi/∂𝜑j    = +vi * [ gij * vj * sin(𝜑i-𝜑j-αij)] 
Hii = ∂pi/∂𝜑i    = -vi * [∑ gik*uk*sin(𝜑i-𝜑j-αik)] + vi*[gii * ui * sin(𝜑i-𝜑i-αi)] : k = 1 .. n -> neighbors

Nij = Vj*∂pi/∂vj = +vi * [ gij * vj * cos(𝜑i-𝜑j-αij) ]
Nii = Vi*∂pi/∂vi = +vi * [∑ gik* vk * cos(𝜑i-𝜑j-αik)] + vi * ( gii * vi * cos(𝜑i-𝜑i-αii) )

Reaktive Power:
Jij = ∂qi/∂𝜑j    = -vi * [gij * vj * cos(𝜑i-𝜑j-αij)]
Jii = ∂qi/∂𝜑i    = +vi * [∑ gik * vk * cos(𝜑i-𝜑j-αik)] - vi*[gii * ui * cos(𝜑i-𝜑i-αi)] : k = 1 .. n -> neighbors

Lij = Vj*∂qi/∂vj = +vi * [ gij * vj * sin(𝜑i-𝜑j-αij) ] }
Lii = Vi*∂qi/∂vi = +vi * [∑ gik* vk * sin(𝜑i-𝜑j-αik)] + vi*( gii * vi * sin(𝜑i-𝜑i-αii) )}
=#
function calcJacobian(Y::AbstractMatrix{ComplexF64}, busVec::Vector{BusData}, adjBranch::Vector{Vector{Int}}, busTypeVec::Vector{NodeType}, slackIdx::Int, n_pq::Int, n_pv::Int, log::Bool = false, sparse::Bool = true)
  global debug
  function printdebug(case::String, grad::String, i::Int, j::Int, val::Float64 = 0.0, bus_i::Int = 0, bus_j::Int = 0)
    if debug
      value = round(val, digits = 3)
      println("$(case): $(grad) [$(i),$(j),$(value)] bus[$(bus_i),$(bus_j)]")
    end
  end

  function shiftIJ(idx::Int)
    i = idx
    if idx >= slackIdx
      i = idx - 1
    end
    return i
  end

  n = n_pq + n_pv
  num_pv_nodes = n_pv
  num_q_nodes = n_pq

  size_j = n_pq * 2 + n_pv
  size_i = size_j

  if debug
    println("\nJacobian: $(size_j) x $(size_i), number of PQ nodes:$num_q_nodes,  number of PV nodes: $num_pv_nodes")
  end

  # Initialization of the Jacobian matrix
  if sparse
    jacobian = spzeros(Float64, size_j, size_i) # Jacobi is not complex but real
  else
    jacobian = zeros(Float64, size_j, size_i) # Jacobi is not complex but real
  end

  # i: Index of Bus
  # j: Index of Bus
  # @inbounds: no bounds checking
  @inbounds for (i, b) in enumerate(busVec)
    vm_i = busVec[i].vm_pu
    va_i = busVec[i].va_rad
    bus_type_i = busVec[i].type

    cnt_pv_i = countNodes(busTypeVec, i, PV)
    for j in adjBranch[i]
      bus_type_j = busVec[j].type
      cnt_pv_j = countNodes(busTypeVec, j, PV)

      i1 = (2 * shiftIJ(i) - 1) - cnt_pv_i
      i2 = 2 * shiftIJ(i) - cnt_pv_i
      j1 = (2 * shiftIJ(j) - 1) - cnt_pv_j
      j2 = 2 * shiftIJ(j) - cnt_pv_j

      if bus_type_i == PQ
        if bus_type_j == PQ
          case = "PQ + PQ"
          if i == j
            Yii = abs(Y[i, i])
            arg = -angle(Y[i, i])

            # Hii = ∂pi/∂𝜑i = vi*[gii * ui * sin(𝜑i-𝜑i-αi)] -vi * [∑ gik*uk*sin(𝜑i-𝜑k-αik)] +  : k = 1 .. n -> neighbors                    
            jacobian[i1, j1] = vm_i * (Yii * vm_i * sin(arg)) - vm_i * sum(abs(Y[i, k]) * busVec[k].vm_pu * sin(va_i - busVec[k].va_rad - angle(Y[i, k])) for k in adjBranch[i])
            # Jii = ∂qi/∂𝜑i    = - vi*[gii * vi * cos(𝜑i-𝜑i-αi)] + vi * [∑ gik * vk * cos(𝜑i-𝜑j-αik)]  : k = 1 .. n -> neighbors                    
            jacobian[i2, j1] = -vm_i * (Yii * vm_i * cos(arg)) + vm_i * sum(abs(Y[i, k]) * busVec[k].vm_pu * cos(va_i - busVec[k].va_rad - angle(Y[i, k])) for k in adjBranch[i])
            # Nii = Vi*∂pi/∂vi = vi * (gii * vi * cos(𝜑i-𝜑i-αii)) +vi * [∑ gik* vk * cos(𝜑i-𝜑j-αik)]
            jacobian[i1, j2] = vm_i * (Yii * vm_i * cos(arg)) + vm_i * sum(abs(Y[i, k]) * busVec[k].vm_pu * cos(va_i - busVec[k].va_rad - angle(Y[i, k])) for k in adjBranch[i])
            # Lii = Vi*∂qi/∂vi = + vi*( gii * vi * sin(𝜑i-𝜑i-αii))  + vi * [∑ gik* vk * sin(𝜑i-𝜑j-αik)]          
            jacobian[i2, j2] = vm_i * (Yii * vm_i * sin(arg)) + vm_i * sum(abs(Y[i, k]) * busVec[k].vm_pu * sin(va_i - busVec[k].va_rad - angle(Y[i, k])) for k in adjBranch[i])

            printdebug(case, "Hii", i1, j1, jacobian[i1, j1], i, j)
            printdebug(case, "Jii", i2, j1, jacobian[i2, j1], i, j)
            printdebug(case, "Nii", i1, j2, jacobian[i1, j2], i, j)
            printdebug(case, "Lii", i2, j2, jacobian[i2, j2], i, j)

          else
            vm_j = busVec[j].vm_pu
            va_j = busVec[j].va_rad
            Yij = abs(Y[i, j])
            arg = va_i - va_j - angle(Y[i, j])

            # Hij = ∂pi/∂𝜑j    = +vi * [ gij * vj * sin(𝜑i-𝜑j-αij)] 
            jacobian[i1, j1] = vm_i * Yij * vm_j * sin(arg)
            # Jij = ∂qi/∂𝜑j    = -vi * [gij * vj * cos(𝜑i-𝜑j-αij)]                
            jacobian[i2, j1] = -vm_i * Yij * vm_j * cos(arg)
            # Nij = Vj*∂pi/∂vj = +vi * gij * vj * cos(𝜑i-𝜑j-αij)                
            jacobian[i1, j2] = vm_i * Yij * vm_j * cos(arg)
            # Lij = Vj*∂qi/∂vj = +vi * [ gij * vj * sin(𝜑i-𝜑j-αij) ] }                
            jacobian[i2, j2] = vm_i * Yij * vm_j * sin(arg)
            printdebug(case, "Hij", i1, j1, jacobian[i1, j1], i, j)
            printdebug(case, "Jij", i2, j1, jacobian[i2, j1], i, j)
            printdebug(case, "Nij", i1, j2, jacobian[i1, j2], i, j)
            printdebug(case, "Lij", i2, j2, jacobian[i2, j2], i, j)
          end
        elseif bus_type_j == PV
          case = "PQ + PV"

          if i == j
            Yii = abs(Y[i, i])
            arg = -angle(Y[i, i])

            # Hii = ∂pi/∂𝜑i    vi*[gii * ui * sin(𝜑i-𝜑i-αi)] -vi * [∑ gik*uk*sin(𝜑i-𝜑k-αik)] +  : k = 1 .. n -> neighbors
            jacobian[i1, j1] = vm_i * (Yii * vm_i * sin(arg)) - vm_i * sum(abs(Y[i, k]) * busVec[k].vm_pu * sin(va_i - busVec[k].va_rad - angle(Y[i, k])) for k in adjBranch[i])
            # Jii = ∂qi/∂𝜑i    = - vi*[gii * vi * cos(𝜑i-𝜑i-αi)] + vi * [∑ gik * vk * cos(𝜑i-𝜑j-αik)]  : k = 1 .. n -> neighbors
            jacobian[i2, j1] = -vm_i * (Yii * vm_i * cos(arg)) + vm_i * sum(abs(Y[i, k]) * busVec[k].vm_pu * cos(va_i - busVec[k].va_rad - angle(Y[i, k])) for k in adjBranch[i])

            printdebug(case, "Hii", i1, j1, jacobian[i1, j1], i, j)
            printdebug(case, "Jii", i2, j1, jacobian[i2, j1], i, j)
          else
            vm_j = busVec[j].vm_pu
            va_j = busVec[j].va_rad
            Yij = abs(Y[i, j])
            arg = va_i - va_j - angle(Y[i, j])

            # Hij = ∂pi/∂𝜑j    = +vi * [ gij * vj * sin(𝜑i-𝜑j-αij)] 
            jacobian[i1, j1] = vm_i * Yij * vm_j * sin(arg)
            # Jij = ∂qi/∂𝜑j    = -vi * [gij * vj * cos(𝜑i-𝜑j-αij)]                
            jacobian[i2, j1] = -vm_i * Yij * vm_j * cos(arg)

            printdebug(case, "Hij", i1, j1, jacobian[i1, j1], i, j)
            printdebug(case, "Jij", i2, j1, jacobian[i2, j1], i, j)
          end
        end # if bus_type_j == PQ           
      elseif bus_type_i == PV
        if bus_type_j == PQ
          case = "PV + PQ"
          # H, N
          if i == j
            Yii = abs(Y[i, i])
            arg = -angle(Y[i, i])

            # Hii = ∂pi/∂𝜑i    vi*[gii * ui * sin(𝜑i-𝜑i-αi)] -vi * [∑ gik*uk*sin(𝜑i-𝜑k-αik)] +  : k = 1 .. n -> neighbors
            jacobian[i1, j1] = vm_i * (Yii * vm_i * sin(arg)) - vm_i * sum(abs(Y[i, k]) * busVec[k].vm_pu * sin(va_i - busVec[k].va_rad - angle(Y[i, k])) for k in adjBranch[i])
            # Nii = Vi*∂pi/∂vi = vi * (gii * vi * cos(𝜑i-𝜑i-αii)) +vi * [∑ gik* vk * cos(𝜑i-𝜑j-αik)]
            jacobian[i2, j1] = vm_i * (Yii * vm_i * cos(arg)) + vm_i * sum(abs(Y[i, k]) * busVec[k].vm_pu * cos(va_i - busVec[k].va_rad - angle(Y[i, k])) for k in adjBranch[i])

            printdebug(case, "Hii", i1, j1, jacobian[i1, j1], i, j)
            printdebug(case, "Nii", i1, j2, jacobian[i1, j2], i, j)
          else
            vm_j = busVec[j].vm_pu
            va_j = busVec[j].va_rad
            Yij = abs(Y[i, j])
            arg = va_i - va_j - angle(Y[i, j])

            # Hij = ∂pi/∂𝜑j    = +vi * [ gij * vj * sin(𝜑i-𝜑j-αij)] 
            jacobian[i1, j1] = vm_i * Yij * vm_j * sin(arg)
            # Nij = Vj*∂pi/∂vj = +vi * gij * vj * cos(𝜑i-𝜑j-αij)                
            jacobian[i1, j2] = vm_i * Yij * vm_j * cos(arg)

            printdebug(case, "Hij", i1, j1, jacobian[i1, j1], i, j)
            printdebug(case, "Nij", i1, j2, jacobian[i1, j2], i, j)
          end
        elseif bus_type_j == PV
          case = "PV + PV"

          if i == j
            Yii = abs(Y[i, i])
            arg = -angle(Y[i, i])

            # Hii = ∂pi/∂𝜑i    vi*[gii * ui * sin(𝜑i-𝜑i-αi)] -vi * [∑ gik*uk*sin(𝜑i-𝜑k-αik)] +  : k = 1 .. n -> neighbors
            jacobian[i1, j1] = vm_i * (Yii * vm_i * sin(arg)) - vm_i * sum(abs(Y[i, k]) * busVec[k].vm_pu * sin(va_i - busVec[k].va_rad - angle(Y[i, k])) for k in adjBranch[i])
            printdebug(case, "Hii", i1, j1, jacobian[i1, j1], i, j)
          else
            vm_j = busVec[j].vm_pu
            va_j = busVec[j].va_rad
            Yij = abs(Y[i, j])
            arg = va_i - va_j - angle(Y[i, j])

            # Hij = ∂pi/∂𝜑j   = +vi * [ gij * vj * sin(𝜑i-𝜑j-αij)] 
            jacobian[i1, j1] = vm_i * Yij * vm_j * sin(arg)

            printdebug(case, "Hij", i1, j1, jacobian[i1, j1], i, j)
          end # if i == j
        end # if bus_type_j == PV           
      end # if bus_type_i == PQ                      
    end # for j in adjBranch[i] 
  end # for i in 1:n

  if log
    println("\nJacobian:\n")
    for j = 1:size_i
      print("\t$(j):")
    end
    println()
    for i = 1:size_i
      print("$(i):\t")
      for j = 1:size_j
        mag = jacobian[i, j]
        mag = round(mag, digits = 1)
        print("$(mag)\t")
      end
      println()
    end
  end

  return jacobian
end

#=
Performs the Newton-Raphson power flow calculation.

# Arguments
- `Y::AbstractMatrix{ComplexF64}`: Y-Bus matrix.
- `nodes::Vector{Node}`: Vector of nodes.
- `Sbase_MVA::Float64`: Base MVA of the system.
- `maxIte::Int`: Maximum number of iterations.
- `tolerance::Float64 = 1e-6`: Tolerance for convergence (default: 1e-6).
- `verbose::Int = 0`: Verbosity level (default: 0).
- `sparse::Bool = false`: Whether to use sparse matrices (default: false).

# Returns
A tuple containing the number of iterations and the result of the calculation:
- `0`: Convergence reached.
- `1`: No convergence.
- `2`: Unsolvable system of equations.
- `3`: Error.
=#
function calcNewtonRaphson!(net::Net, Y::AbstractMatrix{ComplexF64}, maxIte::Int, tolerance::Float64 = 1e-6, verbose::Int = 0, sparse::Bool = false, flatStart::Bool = false, angle_limit::Bool = false, debug::Bool = false)
  nodes = net.nodeVec
  isoNodes = net.isoNodes
  Sbase_MVA = net.baseMVA
  setJacobianAngleLimit(angle_limit)
  setJacobianDebug(debug)

  busVec, slackNum = getBusData(nodes, Sbase_MVA, flatStart)

  adjBranch = adjacentBranches(Y, debug)
  num_pv_nodes = count(bus -> bus.type == PV, busVec)
  num_pq_nodes = count(bus -> bus.type == PQ, busVec)

  size = num_pq_nodes * 2 + num_pv_nodes

  busTypeVec, slackIdx = getBusTypeVec(busVec)
  @assert slackIdx == slackNum "Error: Number of slack nodes is not equal to the number of slack nodes in the bus vector ($(slackIdx) /= $(slackNum))."

  erg = 1 # result of calculation
  # 0 convergence reached
  # 1 no convergence
  # 2 unsolvable system of equations
  # 3 error

  if debug
    tn = num_pq_nodes + num_pv_nodes
    println("\ncalcNewtonRaphson: Matrix-Size: $(size) x $(size), Sbase_MVA: $(Sbase_MVA), total number of nodes: $(tn), number of PQ nodes: $num_pq_nodes, number of PV nodes: $num_pv_nodes\n")
  end

  iteration_count = 0
  power_feeds = getPowerFeeds(busVec, num_pq_nodes, num_pv_nodes)

  while iteration_count <= maxIte
    # Calculation of the residual
    delta_P = residuum(Y, busVec, power_feeds, num_pq_nodes, num_pv_nodes, (verbose > 2))
    if verbose > 2
      println("\ndelta_s: Iteration $(iteration_count)")
      printVector(delta_P, busVec, "p", "q", false, 0.0)
    end

    norm_p = norm(delta_P)
    if verbose > 0
      @printf " norm %e, tol %e, ite %d\n" norm_p tolerance iteration_count
    end

    if norm_p < tolerance
      if (verbose > 0)
        println("\nConvergence is reached after $(iteration_count) iterations")
      end
      erg = 0
      break
    end

    # Calculation of the Jacobian matrix      
    jacobian = calcJacobian(Y, busVec, adjBranch, busTypeVec, slackIdx, num_pq_nodes, num_pv_nodes, (verbose >= 3), sparse)

    # Solution of the linear system of equations for delta_x        
    delta_x = nothing
    try
      # Solver => \
      delta_x = jacobian \ delta_P

      if verbose > 3
        println("\ndelta_v: Iteration $(iteration_count)")
        printVector(delta_x, busVec, "Δφ", "ΔU/U", true, 0.0)
      end
    catch err
      D = det(jacobian)
      println("\nUnsolvable system of equations. Determinate J = $(D), Error = $(err))")
      erg = 2
      break
    end

    try
      # Update of the voltage values
      kidx = 0 # counter for PV nodes
      sidx = 0 # counter for slack node

      for (i, bus) in enumerate(busVec)
        if bus.type == Slack
          sidx += 1 # for slack node
          continue
        end

        index_1 = 2 * (i - sidx) - 1 - kidx
        index_2 = 2 * (i - sidx) - kidx
        if bus.type == PQ
          if angle_limit
            busVec[i].va_rad += min(delta_x[index_1], 0.7)
          else
            busVec[i].va_rad += delta_x[index_1]
          end
          busVec[i].vm_pu += busVec[i].vm_pu * delta_x[index_2]
        elseif bus.type == PV
          if angle_limit
            busVec[i].va_rad += min(delta_x[index_1], 0.7)
          else
            busVec[i].va_rad += delta_x[index_1]
          end
          kidx += 1
        end

        V = busVec[i].vm_pu * exp(im * busVec[i].va_rad)
        busVec[i].vm_pu = abs(V)
        busVec[i].va_rad = angle(V)
      end
    catch err
      println("\nError: $(err)")
      erg = 3
      break
    end

    iteration_count += 1
  end

  if verbose > 2
    println("\nSolution: Iteration $(iteration_count)\n")
  end

  lastNode = length(nodes)
  for bus in busVec
    idx = bus.idx
    if idx in isoNodes
      continue
    end
    idx += count(i -> i < idx, isoNodes)
    vm_pu = bus.vm_pu
    va_deg = rad2deg(bus.va_rad)

    setVmVa!(node = nodes[idx], vm_pu = vm_pu, va_deg = va_deg)
    if bus.type == PV
      nodes[idx]._qƩGen = bus._qRes * Sbase_MVA
    elseif bus.type == Slack
      nodes[idx]._pƩGen = bus._pRes * Sbase_MVA
      nodes[idx]._qƩGen = bus._qRes * Sbase_MVA
    end

    if verbose > 3
      vm_pu = round(nodes[idx]._vm_pu, digits = 3)
      va_deg = round(nodes[idx]._va_deg, digits = 3)
      if va_deg > 180
        va_deg -= 360
        va_deg = round(va_deg, digits = 3)
      end
      if idx == lastNode
        println("(Bus=$(idx), vm=$(vm_pu), va=$(va_deg))")
      else
        print("(Bus=$(idx), vm=$(vm_pu), va=$(va_deg)), ")
      end
    end
  end

  totalBusP = sum(bus -> bus._pRes, busVec)
  totalBusQ = sum(bus -> bus._qRes, busVec)
  setTotalBusPower!(net = net, p = totalBusP, q = totalBusQ)

  return iteration_count, erg
end

"""
Runs the power flow calculation using the Newton-Raphson method.

# Arguments
- `net::Net`: Network data structure.
- `maxIte::Int`: Maximum number of iterations.
- `tolerance::Float64 = 1e-6`: Tolerance for convergence (default: 1e-6).
- `verbose::Int = 0`: Verbosity level (default: 0).

# Returns
A tuple containing the number of iterations and the result of the calculation:
- `0`: Convergence reached.
- `1`: No convergence.
- `2`: Unsolvable system of equations.
- `3`: Error.
"""
function runpf_classic!(net::Net, maxIte::Int, tolerance::Float64 = 1e-6, verbose::Int = 0)
  sparse = false
  printYBus = false
  if length(net.nodeVec) == 0
    @error "No nodes found in $(jpath)"
    return
  elseif length(net.nodeVec) > 60
    sparse = true
  end

  if length(net.nodeVec) < 20
    printYBus = (verbose > 1)
  end

  Y = createYBUS(net = net, sparse = sparse, printYBUS = printYBus)

  return calcNewtonRaphson!(net, Y, maxIte, tolerance, verbose, sparse)
end



"""
    runpf!(net, maxIte, tolerance=1e-6, verbose=0; method=:polar_full)

Unified AC power flow interface.

Arguments:
- `net::Net`: network
- `maxIte::Int`: maximum iterations
- `tolerance::Float64`: mismatch tolerance
- `verbose::Int`: verbosity level
- `method::Symbol`: `:polar_full` (default) or `:rectangular`

Returns:
    (iterations::Int, status::Int)

where `status == 0` indicates convergence.
"""
function runpf!(net::Net,
                maxIte::Int,
                tolerance::Float64=1e-6,
                verbose::Int=0;
                method::Symbol = :polar_full)

    if method === :polar_full        
        return runpf_full!(net, maxIte, tolerance, verbose)
    elseif method === :rectangular        
        return runpf_rectangular!(net, maxIte, tolerance, verbose)
    elseif method === :classic
       return runpf_classic!(net, maxIte, tolerance, verbose)
    else
        error("runpf!: unknown method $(method). Use :polar_full or :rectangular.")
    end
end
