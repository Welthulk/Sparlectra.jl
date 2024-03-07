# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 09.8.2023
# include-file jacobian.jl

# classical Newton-Raphson method

debug = false

mutable struct BusData
  idx::Int        # index of the bus, necessary for sorting  
  vm_pu::Float64  # voltage in pu
  va_rad::Float64 # angle in rad
  pƩ::Float64     # sum of real power
  qƩ::Float64     # sum of reactive power  
  type::ResDataTypes.NodeType

  function BusData(idx::Int, vm_pu::Float64, va_deg::Float64, sumP::Float64, sumQ::Float64, type::ResDataTypes.NodeType)
    new(idx, vm_pu, va_deg, sumP, sumQ, type)
  end

  function Base.show(io::IO, bus::BusData)
    va_deg = round(rad2deg(bus.va_rad), digits = 3)
    print(io, "BusData($(bus.idx), $(bus.vm_pu), $(va_deg)°), $(bus.pƩ), $(bus.qƩ), $(bus.type))")
  end
end

# return vector of BusData
function getBusData(nodes::Vector{ResDataTypes.Node}, Sbase_MVA::Float64, verbose::Int, flatStart::Bool = false)
  busVec = Vector{BusData}()

  slackIdx = 0
  if verbose > 2
    println("\ngetBusData:")
  end

  sumLoad_p = 0.0
  sumLoad_q = 0.0
  sumGen_p = 0.0
  sumGen_q = 0.0

  for (i, n) in enumerate(nodes)
    if n._nodeType == ResDataTypes.Slack
      slackIdx = i
    end
    type = n._nodeType

    p = 0
    q = 0
    #=
    @show n.busIdx
    @show n._pƩLoad
    @show n._qƩLoad
    @show n._pƩGen
    @show n._qƩGen
    =#
    p += n._pƩLoad === nothing ? 0.0 : n._pƩLoad * -1.0
    q += n._qƩLoad === nothing ? 0.0 : n._qƩLoad * -1.0

    sumLoad_p += n._pƩLoad === nothing ? 0.0 : n._pƩLoad * -1.0
    sumLoad_q += n._qƩLoad === nothing ? 0.0 : n._qƩLoad * -1.0
    # Hint: Shunts are considered in Y-Bus Matrix

    p += n._pƩGen === nothing ? 0.0 : n._pƩGen
    q += n._qƩGen === nothing ? 0.0 : n._qƩGen

    sumGen_p += n._pƩGen === nothing ? 0.0 : n._pƩGen
    sumGen_q += n._qƩGen === nothing ? 0.0 : n._qƩGen

    p = p / Sbase_MVA
    q = q / Sbase_MVA

    if flatStart
      if type == ResDataTypes.PQ
        vm_pu = 1.0
        va_deg = 0.0
      elseif type == ResDataTypes.PV
        vm_pu = n._vm_pu === nothing ? 1.0 : n._vm_pu
        va_deg = 0.0
      elseif type == ResDataTypes.Slack
        vm_pu = n._vm_pu === nothing ? 1.0 : n._vm_pu
        va_deg = n._va_deg === nothing ? 0.0 : deg2rad(n._va_deg)
      end
    else
      vm_pu = n._vm_pu === nothing ? 1.0 : n._vm_pu
      va_deg = n._va_deg === nothing ? 0.0 : deg2rad(n._va_deg)
    end
    b = BusData(n.busIdx, vm_pu, va_deg, p, q, type)
    push!(busVec, b)
  end

  if slackIdx == 0
    throw("No slack node found")
  end

  sort!(busVec, by = x -> x.idx)

  sumLoad_p = sumLoad_p * -1.0
  sumLoad_q = sumLoad_q * -1.0
  delta_p = round((sumGen_p - sumLoad_p), digits = 3)
  delta_q = round((sumGen_q - sumLoad_q), digits = 3)

  if verbose > 0
    println("\n∑Load: [$(sumLoad_p), $(sumLoad_q)], ∑Gen [$(sumGen_p), $(sumGen_q)]  Δp, Δq: [$(delta_p), $(delta_q)]")
  end

  if verbose > 2
    println("\nslack bus: $(slackIdx)")
    for b in busVec
      println("$(b)")
    end
  end

  return busVec, slackIdx
end # getBusData

# helper function to count number of nodes of type = value [PQ, PV]
function getBusTypeVec(busVec::Vector{BusData}, log::Bool = false)
  busTypeVec = Vector{ResDataTypes.NodeType}()
  slackIdx = 0
  for (i, bus) in enumerate(busVec) #1:n    
    if bus.type == ResDataTypes.Slack
      slackIdx = i
    end
    push!(busTypeVec, bus.type)
  end
  if log
    for (i, bus) in enumerate(busTypeVec)
      print("busVec[$(i)]: $(bus), ")
    end
  end
  return busTypeVec, slackIdx
end

# count number of nodes of type = value [PQ, PV]
function countNodes(busTypeVec::Vector{ResDataTypes.NodeType}, pos, value::ResDataTypes.NodeType)
  sum = 0

  for (index, bus_type) in enumerate(busTypeVec)
    if index >= pos
      break
    end

    if bus_type == value
      sum += 1
    end
  end

  return sum
end

function printVector(power::Vector{Float64}, busVec::Vector{BusData}, first::String = "p", second::String = "q", eng::Bool = false, delta_min::Float64 = 0.0)
  print("[")
  i = 0
  n = length(busVec)

  for (index, bus) in enumerate(busVec)
    if bus.type == ResDataTypes.PQ
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
    elseif bus.type == ResDataTypes.PV
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

function getPowerFeeds(busVec::Vector{BusData}, n_pq::Int, n_pv::Int, log::Bool = false)::Vector{Float64}
  size = n_pq * 2 + n_pv

  power_feeds = zeros(Float64, size)

  if log
    n = n_pq + n_pv
    println("\nPower Feeds: Elements = $(n), Vector-Size = $(size)")
  end

  vindex = 0 # index of vector of power flows    
  for bus in busVec
    if bus.type == ResDataTypes.PQ
      vindex += 1
      index_p = vindex
      power_feeds[index_p] = bus.pƩ
      vindex += 1
      index_q = vindex
      power_feeds[index_q] = bus.qƩ
    elseif bus.type == ResDataTypes.PV
      vindex += 1
      index_p = vindex
      power_feeds[index_p] = bus.pƩ
    end
  end

  if log
    printVector(power_feeds, busVec)
  end
  return power_feeds
end

"""
 Matrix formulation of calculation of node powers
 S = Vdiag * conj(Y*V) 
 V: vector of complex voltages, 
 Y: Y-Bus Matrix
"""
function residuum(Y::AbstractMatrix{ComplexF64}, busVec::Vector{BusData}, feeders::Vector{Float64}, n_pq::Int, n_pv::Int, log::Bool)::Tuple{Vector{Float64}, ComplexF64}
  # create complex vector of voltages
  V = [bus.vm_pu * exp(im * bus.va_rad) for bus in busVec]
  # create diagonal matrix of voltages
  Vdiag = Diagonal(V)

  # Power Calculation (Knotenleistung)
  S = Vdiag * conj(Y * V)
  S_slack = 0.0 + 0.0im
  size = n_pq * 2 + n_pv
  Δpq = zeros(Float64, size)

  if log
    println("\nresiduum: Vector-Size = $(size)")
  end

  vindex = 0
  pfIdx = 0
  for bus in busVec
    pfIdx += 1
    if bus.type == ResDataTypes.PQ
      vindex += 1
      index_p = vindex
      Δpq[index_p] = feeders[index_p] - real(S[pfIdx])
      vindex += 1
      index_q = vindex
      Δpq[index_q] = feeders[index_q] - imag(S[pfIdx])
    elseif bus.type == ResDataTypes.PV
      vindex += 1
      index_p = vindex
      Δpq[index_p] = feeders[index_p] - real(S[pfIdx])
    elseif bus.type == ResDataTypes.Slack
     S_slack = S[pfIdx]
    end
  end

  if log
    printVector(Δpq, busVec, "Δp", "Δq", false, 0.0)
    norm_delta_p = norm(Δpq)
    println("residuum: $(norm_delta_p)")
    println("S: $S")
    println("feeders: $feeders")
  end

  return Δpq, S_slack
end

"""
classical formulation of calculation of node powers
si = ui ∑ gik**uik* (k = 1...n) # only neighbors!
ui = vi exp(j𝜑i)
uk = vk exp(j𝜑k)
gik* = gik exp(-jαk) = Gik - jBik

si = vi*[∑ vk*(Gik*cos(𝜑i-𝜑k)+Bik*sin(𝜑i-𝜑k)) +jvk*(Gik*sin(𝜑i-𝜑k) - Bik*cos(𝜑i-𝜑k))]

p = vi* (∑ vk * (Gik*cos(𝜑i-𝜑k) + Bik*sin(𝜑i-𝜑k)))
q = vi* (∑ vk * (Gik*sin(𝜑i-𝜑k) - Bik*cos(𝜑i-𝜑k)))
"""
function calcPowerFlow(Y::AbstractMatrix{ComplexF64}, busVec::Vector{BusData}, adjBranch::Vector{Vector{Int}}, n_pq::Int, n_pv::Int, log::Bool = false)::Vector{Float64}
  size = n_pq * 2 + n_pv

  if log
    println("\nPower Flow: Vector-Size = $(size)")
  end

  power_flows = zeros(Float64, size)

  vindex = 0 # index of vector of power flows    
  for (i, bus) in enumerate(busVec)
    if bus.type == ResDataTypes.PQ
      vindex += 1
      index_p = vindex
      vindex += 1
      index_q = vindex

      power_flows[index_p] = bus.vm_pu * sum(abs(Y[i, k]) * busVec[k].vm_pu * cos(busVec[i].va_rad - busVec[k].va_rad - angle(Y[i, k])) for k in adjBranch[i]) # k: index of adjanced bus
      power_flows[index_q] = bus.vm_pu * sum(abs(Y[i, k]) * busVec[k].vm_pu * sin(busVec[i].va_rad - busVec[k].va_rad - angle(Y[i, k])) for k in adjBranch[i]) # k: index of adjanced bus
    elseif bus.type == ResDataTypes.PV
      vindex += 1
      index_p = vindex
      power_flows[index_p] = bus.vm_pu * sum(abs(Y[i, k]) * busVec[k].vm_pu * cos(busVec[i].va_rad - busVec[k].va_rad - angle(Y[i, k])) for k in adjBranch[i]) # k: index of adjanced bus         
    end
  end

  if log
    printVector(power_flows, busVec)
  end
  return power_flows
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
function calcJacobian(Y::AbstractMatrix{ComplexF64}, busVec::Vector{BusData}, adjBranch::Vector{Vector{Int}}, busTypeVec::Vector{ResDataTypes.NodeType}, slackIdx::Int, n_pq::Int, n_pv::Int, log::Bool = false, sparse::Bool = true)
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

  if log
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
  for (i, b) in enumerate(busVec)
    vm_i = busVec[i].vm_pu
    va_i = busVec[i].va_rad
    bus_type_i = busVec[i].type

    cnt_pv_i = countNodes(busTypeVec, i, ResDataTypes.PV)
    for j in adjBranch[i]
      bus_type_j = busVec[j].type
      cnt_pv_j = countNodes(busTypeVec, j, ResDataTypes.PV)

      i1 = (2 * shiftIJ(i) - 1) - cnt_pv_i
      i2 = 2 * shiftIJ(i) - cnt_pv_i
      j1 = (2 * shiftIJ(j) - 1) - cnt_pv_j
      j2 = 2 * shiftIJ(j) - cnt_pv_j

      if bus_type_i == ResDataTypes.PQ
        if bus_type_j == ResDataTypes.PQ
          case = "PQ + PQ"
          @debug "i:$(i), j:$(j), i1:$(i1), i2:$(i2), j1:$(j1), j2:$(j2)"
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
        elseif bus_type_j == ResDataTypes.PV
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
        end # if bus_type_j == ResDataTypes.PQ           
      elseif bus_type_i == ResDataTypes.PV
        if bus_type_j == ResDataTypes.PQ
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
        elseif bus_type_j == ResDataTypes.PV
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
        end # if bus_type_j == ResDataTypes.PV           
      end # if bus_type_i == ResDataTypes.PQ                      
    end # for j in adjBranch[i] 
  end # for i in 1:n

  if log
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

# main function for calculation Newton-Raphson
function calcNewtonRaphson!(Y::AbstractMatrix{ComplexF64}, nodes::Vector{ResDataTypes.Node}, Sbase_MVA::Float64, maxIte::Int, tolerance::Float64 = 1e-6, verbose::Int = 0, sparse::Bool = false)
  @debug global debug = true

  busVec, slackNum = getBusData(nodes, Sbase_MVA, verbose)
  
  adjBranch = adjacentBranches(Y, (verbose > 1))
  num_pv_nodes = count(bus -> bus.type == ResDataTypes.PV, busVec)
  num_pq_nodes = count(bus -> bus.type == ResDataTypes.PQ, busVec)

  size = num_pq_nodes * 2 + num_pv_nodes

  busTypeVec, slackIdx = getBusTypeVec(busVec, (verbose >= 2))
  @assert slackIdx == slackNum "Error: Number of slack nodes is not equal to the number of slack nodes in the bus vector ($(slackIdx) /= $(slackNum))."

  erg = 1 # result of calculation
  # 0 convergence reached
  # 1 no convergence
  # 2 unsolvable system of equations
  # 3 error

  if verbose > 0
    tn = num_pq_nodes + num_pv_nodes
    println("\ncalcNewtonRaphson: Matrix-Size: $(size) x $(size), Sbase_MVA: $(Sbase_MVA), total number of nodes: $(tn), number of PQ nodes: $num_pq_nodes, number of PV nodes: $num_pv_nodes\n")
  end

  iteration_count = 0
  power_feeds = getPowerFeeds(busVec, num_pq_nodes, num_pv_nodes, (verbose > 1))
  s_slack = 0.0 + 0.0im
  while iteration_count <= maxIte
    # Calculation of the power flows
    #=
    power_flows = calcPowerFlow(Y, busVec, adjBranch, num_pq_nodes, num_pv_nodes, (verbose >= 2))      
    if verbose > 0
      printVector(power_flows, busVec, "pFlow","qFlow", false, 0.0)  
    end 
    delta_P = power_feeds - power_flows
    =#
    
    # Calculation of the residual
    delta_P, s_slack = residuum(Y, busVec, power_feeds, num_pq_nodes, num_pv_nodes, (verbose > 1))
    if verbose > 1
      println("\ndelta_P: Iteration $(iteration_count)")
      printVector(delta_P, busVec, "p", "q", false, 0.0)
    end

    norm_p = norm(delta_P)
    if verbose > 0
      @printf " norm %e, tol %e, ite %d\n" norm_p tolerance iteration_count
    end

    if norm_p < tolerance
      println("\nConvergence is reached after $(iteration_count) iterations")
      erg = 0
      break
    end

    # Calculation of the Jacobian matrix      
    jacobian = calcJacobian(Y, busVec, adjBranch, busTypeVec, slackIdx, num_pq_nodes, num_pv_nodes, (verbose > 2), sparse)

    # Solution of the linear system of equations for delta_x        
    delta_x = nothing
    try
      delta_x = jacobian \ delta_P

      #LU-Decomposition:
      # lu_decomposition = lu(jacobian)       
      #delta_x = lu_decomposition \ delta_P
      #iterative Solvers:      
      #delta_x = cg(jacobian, delta_P)

      if verbose > 1
        println("\ndelta_x: Iteration $(iteration_count)")
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
        if bus.type == ResDataTypes.Slack
          sidx += 1 # for slack node
          continue
        end

        index_1 = 2 * (i - sidx) - 1 - kidx
        index_2 = 2 * (i - sidx) - kidx
        if bus.type == ResDataTypes.PQ
          busVec[i].va_rad += delta_x[index_1]
          busVec[i].vm_pu += busVec[i].vm_pu * delta_x[index_2]
        elseif bus.type == ResDataTypes.PV
          busVec[i].va_rad += delta_x[index_1]
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

  if verbose > 1
    println("\nSolution: Iteration $(iteration_count)\n")
  end

  lastNode = length(nodes)
  for bus in busVec
    i = bus.idx
    vm_pu = bus.vm_pu
    va_deg = rad2deg(bus.va_rad)

    setVmVa!(nodes[i], vm_pu, va_deg)
    if i == slackIdx
      s_slack = s_slack*Sbase_MVA
      setGenPower!(nodes[i], real(s_slack), imag(s_slack))    
    end
    if verbose > 1
      vm_pu = round(nodes[i]._vm_pu, digits = 3)
      va_deg = round(nodes[i]._va_deg, digits = 3)
      if va_deg > 180
        va_deg -= 360
        va_deg = round(va_deg, digits = 3)
      end
      if i == lastNode
        println("(Bus=$(i), vm=$(vm_pu), va=$(va_deg))")
      else
        print("(Bus=$(i), vm=$(vm_pu), va=$(va_deg)), ")
      end
    end
  end

  return iteration_count, erg
end
