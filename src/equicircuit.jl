# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 22.05.2023
# include-file equicircuit.jl
"""
    cubicSplineCoefs(x::Vector{Float64}, y::Vector{Float64})::Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}}

Calculates the coefficients of the cubic spline interpolation for the given data points.

# Arguments
- `x::Vector{Float64}`: A vector of x-coordinates of the data points.
- `y::Vector{Float64}`: A vector of y-coordinates of the data points.

# Returns
- `a::Vector{Float64}`: The coefficients for the cubic term.
- `b::Vector{Float64}`: The coefficients for the quadratic term.
- `c::Vector{Float64}`: The coefficients for the linear term.
- `d::Vector{Float64}`: The coefficients for the constant term.

# Example
```julia
x = [1.0, 2.0, 3.0, 4.0]
y = [1.0, 4.0, 9.0, 16.0]
a, b, c, d = cubicSplineCoefs(x, y)
"""
function cubicSplineCoefs(x, y)
  n = length(x)
  h = diff(x)

  # Create a system of equations to calculate the second derivatives
  A = zeros(n, n)
  A[1, 1] = 1.0
  A[n, n] = 1.0

  for i = 2:(n-1)
    A[i, i-1] = h[i-1]
    A[i, i] = 2 * (h[i-1] + h[i])
    A[i, i+1] = h[i]
  end

  b = zeros(n)

  for i = 2:(n-1)
    b[i] = 3 * ((y[i+1] - y[i]) / h[i] - (y[i] - y[i-1]) / h[i-1])
  end

  # Solve the system of equations
  c = A \ b

  # Calculate the remaining coefficients
  a = y[1:(n-1)]
  b = (y[2:n] - y[1:(n-1)]) ./ h - h .* (c[2:n] + 2 * c[1:(n-1)]) / 3
  d = (c[2:n] - c[1:(n-1)]) ./ (3 * h)

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

"""
    calcVKDependence(xTaps::Vector{Int}, yVKs::Vector{Float64}, tapPos::Float64)::Float64

Calculates the voltage dependence on the tap position using cubic spline interpolation.

# Arguments
- `xTaps::Vector{Int}`: A vector of tap positions.
- `yVKs::Vector{Float64}`: A vector of corresponding voltage values.
- `tapPos::Float64`: The current tap position for which the voltage is to be calculated.

# Returns
- `Float64`: The interpolated voltage value at the given tap position.

# Example
```julia
xTaps = [1, 2, 3, 4, 5]
yVKs = [1.0, 1.1, 1.2, 1.3, 1.4]
tapPos = 2.5
voltage = calcVKDependence(xTaps, yVKs, tapPos)
"""
function calcVKDependence(xTaps::Vector{Int}, yVKs::Vector{Float64}, tapPos::Float64)::Float64
  a, b, c, d = cubicSplineCoefs(xTaps, yVKs)
  vk = evaluateCubicSpline(xTaps, a, b, c, d, [tapPos])[1]
  return vk
end

function calcComplexRatio(; tapRatio::Float64, angleInDegrees::Float64)::ComplexF64
  tr = tapRatio * cosd(angleInDegrees)
  ti = tapRatio * sind(angleInDegrees)
  return tr + ti * im
end

"""
    calcNeutralU(neutralU_ratio::Float64, vn_hv::Float64, tap_min::Integer, tap_max::Integer, tap_step_percent::Float64)::Float64

Calculates the neutral voltage of a transformer based on the given parameters.

# Arguments
- `neutralU_ratio::Float64`: The ratio of the neutral voltage to the rated high voltage.
- `vn_hv::Float64`: The rated high voltage of the transformer.
- `tap_min::Integer`: The minimum tap position.
- `tap_max::Integer`: The maximum tap position.
- `tap_step_percent::Float64`: The percentage change in voltage per tap step.

# Returns
- `Float64`: The calculated neutral voltage.

# Example
```julia
neutral_voltage = calcNeutralU(1.0, 110.0, -10, 10, 1.25)
"""
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
function toPU_RXBG(; r::T, x::T, g::Union{Nothing,T} = nothing, b::Union{Nothing,T} = nothing, v_kv::T, baseMVA::T)::NTuple{4,T} where {T<:Real}
  z_base = (v_kv * v_kv) / baseMVA
  y_base = inv(z_base)

  r_pu = r * y_base          # r / z_base
  x_pu = x * y_base          # x / z_base
  g_pu = isnothing(g) ? zero(T) : g * z_base   # g / y_base
  b_pu = isnothing(b) ? zero(T) : b * z_base   # b / y_base

  return r_pu, x_pu, b_pu, g_pu
end

"""
    fromPU_RXBG(r_pu::Float64, x_pu::Float64, g_pu::Union{Nothing,Float64} = nothing, b_pu::Union{Nothing,Float64} = nothing, v_kv::Float64, baseMVA::Float64)::Tuple{Float64,Float64,Float64,Float64}

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
fromPU_RXBG(r_pu = 0.01, x_pu = 0.1, g_pu = 0.02, b_pu = 0.02, v_kv = 110.0, baseMVA = 100.0)
```
"""

function fromPU_RXBG(; r_pu::T, x_pu::T, g_pu::Union{Nothing,T} = nothing, b_pu::Union{Nothing,T} = nothing, v_kv::T, baseMVA::T)::NTuple{4,T} where {T<:Real}
  z_base = (v_kv * v_kv) / baseMVA
  y_base = inv(z_base)

  r = r_pu * z_base
  x = x_pu * z_base
  g = isnothing(g_pu) ? zero(T) : g_pu * y_base
  b = isnothing(b_pu) ? zero(T) : b_pu * y_base

  return r, x, b, g
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
function createYBUS(; net::Net, sparse::Bool = true, printYBUS::Bool = false)

  # Determine the maximum bus number in the network, taking isolated buses into account
  max_bus = maximum(max(branch.fromBus, branch.toBus) for branch in net.branchVec)
  n = max_bus - length(net.isoNodes)

  Y = sparse ? spzeros(ComplexF64, n, n) : zeros(ComplexF64, n, n)

  for branch in net.branchVec
    @assert branch isa AbstractBranch "branch $(branch) is not of type AbstractBranch"
    fromNode = branch.fromBus
    toNode = branch.toBus
    # Skip branches that are out of service
    if branch.status == 0
      @debug "createYBUS: Branch $(branch) out of service, skipping "
      continue
    end
    # skip isolated nodes
    if fromNode in net.isoNodes || toNode in net.isoNodes
      continue
    end

    # correct the indices based on the isolated nodes
    fromNode -= count(i -> i < fromNode, net.isoNodes)
    toNode -= count(i -> i < toNode, net.isoNodes)
    # Y-Matrix
    y_11, y_12, y_21, y_22 = calcAdmittance(branch, branch.comp.cVN, net.baseMVA)
    Y[fromNode, fromNode] += y_11
    Y[toNode, toNode] += y_22
    Y[fromNode, toNode] += y_12
    Y[toNode, fromNode] += y_21
  end

  for sh in net.shuntVec
    sh.status == 0 && continue
    bus = sh.busIdx
    bus in net.isoNodes && continue
    bus -= count(i -> i < bus, net.isoNodes)

    # Shunt-Admittanz diagonal addieren
    Y[bus, bus] += sh.y_pu_shunt
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
