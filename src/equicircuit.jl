# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 22.05.2023
# include-file equicircuit.jl
function cubicSplineCoefs(x, y)
  n = length(x)
  h = diff(x)

  # Create a system of equations to calculate the second derivatives
  A = zeros(n, n)
  A[1, 1] = 1.0
  A[n, n] = 1.0

  for i = 2:n-1
    A[i, i-1] = h[i-1]
    A[i, i] = 2 * (h[i-1] + h[i])
    A[i, i+1] = h[i]
  end

  b = zeros(n)

  for i = 2:n-1
    b[i] = 3 * ((y[i+1] - y[i]) / h[i] - (y[i] - y[i-1]) / h[i-1])
  end

  # Solve the system of equations
  c = A \ b

  # Calculate the remaining coefficients
  a = y[1:n-1]
  b = (y[2:n] - y[1:n-1]) ./ h - h .* (c[2:n] + 2 * c[1:n-1]) / 3
  d = (c[2:n] - c[1:n-1]) ./ (3 * h)

  return a, b, c, d
end

function evaluateCubicSpline(x, a, b, c, d, x_new)
  n = length(x)
  y_new = similar(x_new)
  range = 1:length(x_new)
  for i in range
    if x_new[i] < x[1] || x_new[i] > x[n]
      error("x_new value is out of range")
    end

    index = searchsortedlast(x, x_new[i])

    if index == 0
      index = 1
    elseif index == n
      index = n - 1
    end

    t = (x_new[i] - x[index]) / (x[index+1] - x[index])

    y_new[i] = a[index] + b[index] * t + c[index] * t^2 + d[index] * t^3
  end

  return y_new
end

function calcVKDependence(xTaps::Vector{Int}, yVKs::Vector{Float64}, tapPos::Float64)::Float64
  a, b, c, d = cubicSplineCoefs(xTaps, yVKs)
  vk = evaluateCubicSpline(xTaps, a, b, c, d, [tapPos])[1]
  return vk
end

function calcComplexRatio(tapRatio::Float64, ShiftAngle_grad::Float64)::ComplexF64
  tr = tapRatio * cosd(ShiftAngle_grad)
  ti = tapRatio * sind(ShiftAngle_grad)
  return tr + ti * im
end

function calcNeutralU(neutralU_ratio::Float64, vn_hv::Float64, tap_min::Integer, tap_max::Integer, tap_step_percent::Float64)::Float64
  return round(neutralU_ratio * vn_hv + (tap_max - tap_min) * tap_step_percent / 100.0, digits = 0)
end

"""
    toPU_RXGB(; r::Float64, x::Float64, g::Union{Nothing, Float64}=nothing, b::Union{Nothing, Float64}=nothing, v_kv::Float64, baseMVA::Float64)::Tuple{Float64, Float64, Float64, Float64}

Converts the resistance, reactance, conductance, and susceptance from physical units to per unit.

# Arguments
- `r::Float64`: The resistance in Ohm.
- `x::Float64`: The reactance in Ohm.
- `g::Union{Nothing, Float64}`: The conductance in S. It can be `Nothing` or a `Float64` value.
- `b::Union{Nothing, Float64}`: The susceptance in S. It can be `Nothing` or a `Float64` value.
- `v_kv::Float64`: The voltage in kV.
- `baseMVA::Float64`: The base power in MVA.

# Returns
- `r_pu::Float64`: The per unit resistance.
- `x_pu::Float64`: The per unit reactance.
- `g_pu::Float64`: The per unit conductance.
- `b_pu::Float64`: The per unit susceptance.

# Example
```julia
toPU_RXGB(r = 0.01, x = 0.1, g = 0.02, b = 0.02, v_kv = 110.0, baseMVA = 100.0)
```
"""
function toPU_RXGB(; r::Float64, x::Float64, g::Union{Nothing,Float64} = nothing, b::Union{Nothing,Float64} = nothing, v_kv::Float64, baseMVA::Float64)::Tuple{Float64,Float64,Float64,Float64}
  z_base = (v_kv * v_kv) / baseMVA
  y_base = 1.0 / z_base
  r_pu = r / z_base
  x_pu = x / z_base

  g_pu = b_pu = 0.0
  !isnothing(g) ? g_pu = g * y_base : 0.0
  !isnothing(b) ? b_pu = b * y_base : 0.0

  return r_pu, x_pu, g_pu, b_pu
end

"""
    to_RXGB(r_pu::Float64, x_pu::Float64, g_pu::Union{Nothing,Float64} = nothing, b_pu::Union{Nothing,Float64} = nothing, v_kv::Float64, baseMVA::Float64)::Tuple{Float64,Float64,Float64,Float64}

Converts the resistance, reactance, conductance, and susceptance from per unit to physical units.

# Arguments
- `r_pu::Float64`: The per unit resistance.
- `x_pu::Float64`: The per unit reactance.
- `g_pu::Union{Nothing, Float64}`: The per unit conductance. It can be `Nothing` or a `Float64` value.
- `b_pu::Union{Nothing, Float64}`: The per unit susceptance. It can be `Nothing` or a `Float64` value.
- `v_kv::Float64`: The voltage in kV.
- `baseMVA::Float64`: The base power in MVA.

# Returns
- `r::Float64`: The resistance in Ohm.
- `x::Float64`: The reactance in Ohm.
- `g::Float64`: The conductance in S.
- `b::Float64`: The susceptance in S.

# Example
```julia
to_RXGB(r_pu = 0.01, x_pu = 0.1, g_pu = 0.02, b_pu = 0.02, v_kv = 110.0, baseMVA = 100.0)
```
"""
function to_RXGB(; r_pu::Float64, x_pu::Float64, g_pu::Union{Nothing,Float64} = nothing, b_pu::Union{Nothing,Float64} = nothing, v_kv::Float64, baseMVA::Float64)::Tuple{Float64,Float64,Float64,Float64}
  z_base = (v_kv^2) / baseMVA
  r = r_pu * z_base
  x = x_pu * z_base
  g = b = 0.0
  !isnothing(g_pu) ? g = g_pu / z_base : 0.0
  !isnothing(b_pu) ? b = b_pu / z_base : 0.0
  return r, x, g, b
end
#=
function removeIsolatedNodesFromYBUS(Y::AbstractMatrix{ComplexF64}, isolated_nodes::Vector{Int})
    if issparse(Y)
        # Entferne die Zeilen und Spalten der isolierten Knoten aus der Sparse-Matrix
        n = size(Y, 1)
        rows_to_keep = setdiff(1:n, isolated_nodes)
        Y = Y[rows_to_keep, rows_to_keep]
    else
        # Entferne die Zeilen und Spalten der isolierten Knoten aus der Dense-Matrix
        Y = Y[setdiff(1:end, isolated_nodes), setdiff(1:end, isolated_nodes)]
    end
    return Y
end
=#


"""
    createYBUS(branchVec::Vector{Branch}, shuntVec::Vector{Shunt}, isoNodes::Vector{Int}, sparse::Bool = true, printYBUS::Bool = false)

Creates the bus admittance matrix (YBUS) of the network.

# Arguments
- `branchVec::Vector{Branch}`: The vector of branches in the network.
- `shuntVec::Vector{Shunt}`: The vector of shunts in the network.
- `isoNodes::Vector{Int}`: The vector of isolated nodes in the network.
- `sparse::Bool`: A flag to indicate if the YBUS matrix should be sparse. Default is `true`.
- `printYBUS::Bool`: A flag to indicate if the YBUS matrix should be printed. Default is `false`.

# Returns
- `Y::Matrix{ComplexF64}`: The bus admittance matrix (YBUS).

"""

function createYBUS(;net::Net, sparse::Bool = true, printYBUS::Bool = false)
  
  # Bestimme die maximale Busnummer im Netzwerk, unter Berücksichtigung isolierter Busse
  max_bus = maximum(max(branch.fromBus, branch.toBus) for branch in net.branchVec)
  n = max_bus - length(net.isoNodes)

  @debug "Dimension YBus:", n

  if sparse
    Y = spzeros(ComplexF64, n, n)
  else
    Y = zeros(ComplexF64, n, n)
  end

  if debug
    if sparse
      t = "sparse"
    else
      t = "normal"
    end
    println("\nYBUS: Size = $(n)x$(n) ($(t))\n")
  end

  for branch in net.branchVec
    fromNode = branch.fromBus
    toNode = branch.toBus

    # Überspringe Zweige, die außer Betrieb sind
    if branch.status == 0
      @debug "createYBUS: Branch $(branch) out of service, skipping "
      continue
    end

    # Überspringe Zweige, die isolierte Busse verbinden
    if fromNode in net.isoNodes || toNode in net.isoNodes
      continue
    end

    r = branch.r_pu
    x = branch.x_pu
    b = branch.b_pu
    g = branch.g_pu

    yik = inv((r + x * im))
    susceptance = 0.5 * (g + b * im) # pi-model

    t = 1.0 + 0.0 * im
    ratio = branch.ratio
    shift_degree = branch.angle
    if ratio != 0.0 || shift_degree != 0.0
      t = calcComplexRatio(ratio, shift_degree)
    end

    # Korrigiere die Indizes basierend auf den isolierten Knoten
    fromNode -= count(i -> i < fromNode, net.isoNodes)
    toNode -= count(i -> i < toNode, net.isoNodes)
    
    if branch.fromBusSwitch == 1
      Y[fromNode, fromNode] += (yik + susceptance) / abs2(t)
    else
      Y[fromNode, fromNode] += susceptance + (yik*susceptance/(yik + susceptance))      
    end  
    
    if branch.toBusSwitch == 1
      Y[toNode, toNode] += (yik + susceptance)
    else
      Y[toNode, toNode] += susceptance + (yik*susceptance/(yik + susceptance))      
    end
    
    if branch.fromBusSwitch == 1 
      Y[fromNode, toNode] += (-1.0 * yik / conj(t))
    end
    
    if branch.toBusSwitch == 1  
      Y[toNode, fromNode] += (-1.0 * yik / t)
    end



  end

  for sh in net.shuntVec
    node = sh.busIdx
    if node in net.isoNodes      
      continue
    end
    node -= count(i -> i < node, net.isoNodes)    
    y = sh.y_pu_shunt
    Y[node, node] += y
  end

  if printYBUS
    println("\nYBUS:\n")
    red_text = "\x1b[31m"  # ANSI escape code for red text
    reset_text = "\x1b[0m"  # ANSI escape code to reset text color
    green_text = "\x1b[32m"  # ANSI escape code for green text

    for j = 1:n
      farbe = reset_text
      if j > 1
        print("\t\t$farbe$j:$reset_text")
      else
        print("\t$farbe$j:$reset_text")
      end
    end
    println()

    for i = 1:n
      farbe = reset_text
      print("$farbe$i$reset_text:\t")
      for j = 1:n
        mag = abs(Y[i, j])
        real_part = real(Y[i, j])
        if real_part <= 0
          mag = -mag
          phase = angle(-1 * Y[i, j])
        else
          phase = angle(Y[i, j])
        end
        mag = round(mag, digits = 2)
        phase_deg = round(rad2deg(phase), digits = 3)

        farbe = reset_text

        if isnan(mag) || isnan(phase_deg)
          print("$red_text NaN \t\t$reset_text")
        elseif mag == 0
          print("$(farbe)0.0$(reset_text) \t\t")  # Print extra space for pure zeros
        else
          print("$(farbe)$(mag)∠$(phase_deg)$(reset_text)\t")
        end
      end
      println()
    end
  end

  return Y
end

"""
adjacentBranches: Find adjacent branches for each node in the network.

Parameters:
- `Y::AbstractMatrix{ComplexF64}`: Admittance matrix of the network.
- `log::Bool = false`: Optional parameter indicating whether to print the adjacent branches (default is false).

Returns:
- `adjList::Vector{Vector{Int}}`: Vector of vectors containing the indices of adjacent branches for each node.

"""
function adjacentBranches(Y::AbstractMatrix{ComplexF64}, log::Bool = false)::Vector{Vector{Int}}
  n = size(Y, 1)
  adjList = Vector{Vector{Int}}(undef, n)

  for i = 1:n
    adjList[i] = Vector{Int}()

    for j = 1:n
      if Y[i, j] != 0.0
        push!(adjList[i], j)
      end
    end
  end
  if log
    println("\nadjacentBranches:\n")
    for (i, adj) in enumerate(adjList)
      println("node $i -> adjacent: $adj")
    end
  end
  return adjList
end
