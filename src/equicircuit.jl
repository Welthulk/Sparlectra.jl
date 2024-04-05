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
createYBUS: Create the admittance matrix YBUS from the branch vector and the slack index.

Parameters:
- `branchVec::Vector{Branch}`: Vector of branch objects representing the network branches.
- `shuntVec::Vector{Shunt}`: Vector of shunt objects representing shunt elements in the network.
- `sparse::Bool = true`: Optional parameter indicating whether to create a sparse YBUS matrix (default is true).
- `printYBUS::Bool = false`: Optional parameter indicating whether to print the YBUS matrix (default is false).

Returns:
- `Y::SparseMatrixCSC{ComplexF64}` or `Array{ComplexF64}`: The YBUS admittance matrix.
"""
function createYBUS(branchVec::Vector{Branch}, shuntVec::Vector{Shunt}, sparse::Bool = true, printYBUS::Bool = false)
  @assert length(branchVec) > 0 "branchVec must not be empty"

  n = maximum(max(branch.fromBus, branch.toBus) for branch in branchVec)

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

  for branch in branchVec
    fromNode = branch.fromBus
    toNode = branch.toBus
    r = x = b = g = 0.0
    if Int(branch.status) == 0
      @debug "createYBUS: Branch $(branch) out of service, skipping "
      continue
    else
      r = branch.r_pu
      x = branch.x_pu
      b = branch.b_pu
      g = branch.g_pu
    end

    yik = inv((r + x * im))
    susceptance = 0.5 * (g + b * im) # pi-model

    t = 1.0 + 0.0 * im
    ratio = branch.ratio
    shift_degree = branch.angle
    if ratio != 0.0 || shift_degree != 0.0
      t = calcComplexRatio(ratio, shift_degree)
    end

    Y[fromNode, fromNode] = Y[fromNode, fromNode] + ((yik + susceptance)) / abs2(t)
    Y[toNode, toNode] = Y[toNode, toNode] + (yik + susceptance)

    Y[fromNode, toNode] = Y[fromNode, toNode] + (-1.0 * yik / conj(t))
    Y[toNode, fromNode] = Y[toNode, fromNode] + (-1.0 * yik / t)
  end

  for sh in shuntVec
    node = sh.busIdx
    y = sh.y_pu_shunt
    Y[node, node] = Y[node, node] + y
  end

  if printYBUS
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
