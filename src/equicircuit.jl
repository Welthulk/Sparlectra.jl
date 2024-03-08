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

"""
purpose: calculate the voltage dependence of a tap position
+ xTaps: vector of tap positions
+ yVKs: vector of voltage ratios
+ tapPos: tap position
"""
function calcVKDependence(xTaps::Vector{Int}, yVKs::Vector{Float64}, tapPos::Float64)::Float64
  a, b, c, d = cubicSplineCoefs(xTaps, yVKs)
  vk = evaluateCubicSpline(xTaps, a, b, c, d, [tapPos])[1]
  return vk
end

"""
purpose: calculate the complex tap ratio
+ tapRatio: tap ratio
+ ShiftAngle_grad: shift angle in degree
"""
function calcComplexRatio(tapRatio::Float64, ShiftAngle_grad::Float64)::ComplexF64
  tr = tapRatio * cosd(ShiftAngle_grad)
  ti = tapRatio * sind(ShiftAngle_grad)
  return tr + ti * im
end

"""
purpose: calculate of neutral voltage. Tap position in neutral position.
+ neutralU_ratio: neutral voltage ratio
+ vn_hv: rated voltage high voltage side in kV
+ tap_min: minimum tap position
+ tap_max: maximum tap position
+ tap_step_percent: tap step in percent
"""
function calcNeutralU(neutralU_ratio::Float64, vn_hv::Float64, tap_min::Integer, tap_max::Integer, tap_step_percent::Float64)::Float64
  return round(neutralU_ratio * vn_hv + (tap_max - tap_min) * tap_step_percent / 100.0, digits = 0)
end


"""
purpose: recalculation model data of transformer
input:
baseMVA: 
Sn_MVA: rated power of transformer, MVA
ratedU_kV: rated voltage, kV
r_pu: short circuit resistance, p.u.
x_pu: short circuit reaktance, p.u.
bm: open loop magnitacation, p.u.
hint: gm is set to zero, g_pu = 0.0!
"""
function recalc_trafo_model_data(; baseMVA::Float64, Sn_MVA::Float64, ratedU_kV::Float64, r_pu::Float64, x_pu::Float64, b_pu::Float64, isPerUnit::Bool)::Tuple{Float64,Float64,Float64,Float64}
  wurz3 = sqrt(3.0)
  Pfe_kW = 0.0 # -> g_pu = 0.0
  if isPerUnit
    base_power_3p = baseMVA
    base_power_3p = baseMVA
    base_power_1p = base_power_3p / 3.0
    base_i_to = base_power_3p / (ratedU_kV * wurz3)
    base_y_to = base_i_to * base_i_to / base_power_1p
    base_z_to = 1.0 / base_y_to

    # g_pu = 0.0!
    y_shunt = b_pu
    Y_shunt = y_shunt * base_y_to
    i0 = 100.0 * Y_shunt * ratedU_kV^2 / Sn_MVA

    z_ser_pu = sqrt(r_pu^2 + x_pu^2)
    Z_Series_abs = z_ser_pu * base_z_to
    uk = Z_Series_abs * Sn_MVA / ratedU_kV^2
    Rk = r_pu * base_z_to
    P_kW = 100.0 * Rk * Sn_MVA^2 / ratedU_kV^2
  else
    u_quad = ratedU_kV^2
    z_base = u_quad / Sn_MVA
    y_base = Sn_MVA / u_quad

    Y_shunt = b_pu
    i0 = Y_shunt * z_base * 100.0

    Z_Series_abs = sqrt(r_pu^2 + x_pu^2)
    uk = Z_Series_abs * y_base

    Rk = r_pu
    P_kW = 100.0 * Rk * Sn_MVA^2 / u_quad
  end

  return uk, P_kW, i0, Pfe_kW
end


"""
createYBUS: Create the admittance matrix YBUS from the branch vector and the slack index
#cs: counter system (or meetering system?) (GER: "Zählpfeilsystem") [1= Consumer, 2= Producer]
#tc: tap Changer side in equicircuit: [fromBus=1, toBus=2]
"""
function createYBUS(branchVec::Vector{Branch}, shuntVec::Vector{ResDataTypes.Shunt}, sparse::Bool = true, printYBUS::Bool = false)
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
    susceptance = 0.5*(g + b * im) # pi-model

    t = 1.0 + 0.0 * im
    ratio = branch.ratio
    shift_degree = branch.angle
    if ratio != 0.0 || shift_degree != 0.0
      t = calcComplexRatio(ratio, shift_degree)
    end

    Y[toNode, toNode] = Y[toNode, toNode] +  (yik + susceptance)
    Y[fromNode, fromNode] = Y[fromNode, fromNode] + ((yik + susceptance)) / abs2(t)
    
    Y[fromNode, toNode] = Y[fromNode, toNode] + (-1.0 * yik / conj(t))
    Y[toNode, fromNode] = Y[toNode, fromNode] + (-1.0* yik / t)
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

