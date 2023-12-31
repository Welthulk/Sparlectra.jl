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
purpos: calculation of a two port in p.u.
+ baseV: base voltage in kV
+ baseS: base power in MVA
+ r: resistance in Ohm
+ x: reactance in Ohm
+ b: susceptance in Ohm
+ g: conductance in Ohm
"""
function calcTwoPortPU(baseV::Float64, baseS::Float64, r::Float64, x::Float64, b::Float64, g::Float64)
  baseZ = (baseV)^2 / baseS
  r_pu = r / baseZ
  x_pu = x / baseZ
  b_pu = b * baseZ
  g_pu = g * baseZ
  return r_pu, x_pu, b_pu, g_pu
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
purpose: calculate of ratio with correction of tap position
+ tapPos: tap position
+ tapNeutral: neutral tap position
+ tapStepPercent: tap step in percent
"""
function calcTapCorr(tapPos::Int, tapNeutral::Int, tapStepPercent::Float64)::Float64
  return 1 + (tapPos - tapNeutral) * tapStepPercent / 100
end

"""
purpose: calculate of ratio with correction of tap position
+ bus_hv_kv: high voltage bus voltage in kV
+ vn_hv_kv: rated voltage high voltage side in kV
+ bus_lv_kv: low voltage bus voltage in kV
+ vn_lv_kv: rated voltage low voltage side in kV
+ tapPos: tap position
+ tapNeutral: neutral tap position
+ tapStepPercent: tap step in percent
"""
function calcRatio(bus_hv_kv::Float64, vn_hv_kv::Float64, bus_lv_kv::Float64, vn_lv_kv::Float64, tapPos::Integer, tapNeutral::Integer, tapStepPercent::Float64, side::Integer)::Float64
  @assert side in [1, 2] "side must be 1 or 2"

  ratio = bus_hv_kv * vn_lv_kv / (bus_lv_kv * vn_hv_kv)

  if tapPos == tapNeutral
    return ratio
  end

  corr = calcTapCorr(tapPos, tapNeutral, tapStepPercent)

  return ratio / corr
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
 purpose: calculates R, X, rfe and xmu of a transformer
 + sn_mva: rated power in MVA
 + vn_hv_kv: rated voltage high voltage side in kV
 + vk_percent: short circuit voltage in percent
 + vkr_percent: real part short circuit voltage in percent  
 + pfe_kw: iron losses in kW
 + i0_percent: no load current in percent (open circuit)
 
"""
function calcTrafoParams(sn_mva::Float64, vn_hv_kv::Float64, vk_percent::Float64, vkr_percent::Float64, pfe_kw::Float64 = 0.0, i0_percent::Float64 = 0.0)
  @assert sn_mva > 0.0 "sn_mva must be > 0.0"
  @assert vn_hv_kv > 0.0 "vn_hv_kv must be > 0.0"
  @assert vk_percent > 0.0 "vk_percent must be > 0.0"
  @assert vkr_percent > 0.0 "vkr_percent must be > 0.0"

  fak = vn_hv_kv^2 / sn_mva

  # Impedanz
  zk = vk_percent / 100.0 * fak

  # Resistanz
  rk = vkr_percent / 100.0 * fak

  # Reaktanz
  xk = sqrt(zk^2 - rk^2) # Reaktanz

  # Suszeptanz
  if pfe_kw > 0.0 && i0_percent > 0.0
    pfe = pfe_kw / 1000.0
    v_quad = vn_hv_kv^2
    Yfe = pfe / v_quad
    Y0 = i0_percent / 100.0 * sn_mva / v_quad
    gm = Yfe
    try
      bm = -1.0 * sqrt(Y0^2 - Yfe^2)
    catch
      @warn "bm is set to 0.0 (Y0^2 - Yfe^2 < 0.0)"
      bm = 0.0
    end
  else
    gm = 0.0
    bm = 0.0
  end

  return rk, xk, bm, gm
end

#=
purpose: calculate transformer parameters for 3WT-Transformer
+ side: side of transformer (1, 2 or 3)
+ sn_hv_mva: rated power in MVA
+ sn_mv_mva: rated power in MVA
+ sn_lv_mva: rated power in MVA
+ vn_hv_kv: rated voltage high voltage side in kV
+ vn_mv_kv: rated voltage medium voltage side in kV
+ vn_lv_kv: rated voltage low voltage side in kV
+ vk_hv_percent: short circuit voltage high voltage side in percent
+ vk_mv_percent: short circuit voltage medium voltage side in percent
+ vk_lv_percent: short circuit voltage low voltage side in percent
+ vkr_hv_percent: real part short circuit voltage high voltage side in percent
+ vkr_mv_percent: real part short circuit voltage medium voltage side in percent
+ vkr_lv_percent: real part short circuit voltage low voltage side in percent
+ pfe_kw: iron losses in kW
+ i0_percent: no load current in percent (open circuit)
further information see: "MVA Method for Three-Winding Transformer" by Ver Pangonilo
Model: 3WT-Transformer

#                              hv_bus
#                                |
#                                T1                               
#                                |    
#                                *auxilary_bus -> V = hv_bus
#                               / \
#                              /   \ 
#                             T2    T3
#                            /       \  
#                         mv_bus    lv_bus   
=#
function calc3WTParams(
  side::Integer,
  sn_hv_mva::Float64,
  sn_mv_mva::Float64,
  sn_lv_mva::Float64,
  vn_hv_kv::Float64,
  vn_mv_kv::Float64,
  vn_lv_kv::Float64,
  vk_hv_percent::Float64,
  vk_mv_percent::Float64,
  vk_lv_percent::Float64,
  vkr_hv_percent::Float64,
  vkr_mv_percent::Float64,
  vkr_lv_percent::Float64,
  pfe_kw::Float64 = 0.0,
  i0_percent::Float64 = 0.0,
)
  @assert side in [1, 2, 3] "side must be 1, 2 or 3"
  @assert sn_hv_mva > 0.0 "sn_hv_mva must be > 0.0"
  @assert sn_mv_mva > 0.0 "sn_mv_mva must be > 0.0"
  @assert sn_lv_mva > 0.0 "sn_lv_mva must be > 0.0"
  @assert vn_hv_kv > 0.0 "vn_hv_kv must be > 0.0"
  @assert vn_mv_kv > 0.0 "vn_mv_kv must be > 0.0"
  @assert vn_lv_kv > 0.0 "vn_lv_kv must be > 0.0"
  @assert vk_hv_percent > 0.0 "vk_hv_percent must be > 0.0"
  @assert vk_mv_percent > 0.0 "vk_mv_percent must be > 0.0"
  @assert vk_lv_percent > 0.0 "vk_lv_percent must be > 0.0"
  @assert vkr_hv_percent > 0.0 "vkr_hv_percent must be > 0.0"
  @assert vkr_mv_percent > 0.0 "vkr_mv_percent must be > 0.0"
  @assert vkr_lv_percent > 0.0 "vkr_lv_percent must be > 0.0"

  # magnitude         
  vk1_hm = vk_hv_percent * sn_hv_mva / (min(sn_hv_mva, sn_mv_mva))
  vk1_ml = vk_mv_percent * sn_hv_mva / (min(sn_mv_mva, sn_lv_mva))
  vk1_lh = vk_lv_percent * sn_hv_mva / (min(sn_hv_mva, sn_lv_mva))

  vk1_T1 = 0.5 * (vk1_hm + vk1_lh - vk1_ml)
  vk1_T2 = 0.5 * (vk1_ml + vk1_hm - vk1_lh)
  vk1_T3 = 0.5 * (vk1_ml + vk1_lh - vk1_hm)

  vk_T1 = vk1_T1
  vk_T2 = vk1_T2 * sn_mv_mva / sn_hv_mva
  vk_T3 = vk1_T3 * sn_lv_mva / sn_hv_mva

  # real part
  vkr1_hm = vkr_hv_percent * sn_hv_mva / (min(sn_hv_mva, sn_mv_mva))
  vkr1_ml = vkr_mv_percent * sn_hv_mva / (min(sn_mv_mva, sn_lv_mva))
  vkr1_lh = vkr_lv_percent * sn_hv_mva / (min(sn_hv_mva, sn_lv_mva))

  vkr1_T1 = 0.5 * (vkr1_hm + vkr1_lh - vkr1_ml)
  vkr1_T2 = 0.5 * (vkr1_ml + vkr1_hm - vkr1_lh)
  vkr1_T3 = 0.5 * (vkr1_ml + vkr1_lh - vkr1_hm)

  vkr_T1 = vkr1_T1
  vkr_T2 = vkr1_T2 * sn_mv_mva / sn_hv_mva
  vkr_T3 = vkr1_T3 * sn_lv_mva / sn_hv_mva

  # calculate transformer parameters
  if side == 1
    rk_T1, xk_T1, bm_T1, gm_T1 = calcTrafoParams(sn_hv_mva, vn_hv_kv, vk_T1, vkr_T1, pfe_kw, i0_percent)
    return rk_T1, xk_T1, bm_T1, gm_T1
  elseif side == 2
    rk_T2, xk_T2, bm_T2, gm_T2 = calcTrafoParams(sn_mv_mva, vn_mv_kv, vk_T2, vkr_T2, 0.0, 0.0)
    return rk_T2, xk_T2, bm_T2, gm_T2
  elseif side == 3
    rk_T3, xk_T3, bm_T3, gm_T3 = calcTrafoParams(sn_lv_mva, vn_lv_kv, vk_T3, vkr_T3, 0.0, 0.0)
    return rk_T3, xk_T3, bm_T3, gm_T3
  end
end

function calcYShunt(pShunt_MW::Float64, qShunt_MVA::Float64, ratio::Float64, SBase_MVA::Float64)
  p_shunt = pShunt_MW * ratio
  q_shunt = qShunt_MVA * ratio
  Sh_ref = Complex(p_shunt, q_shunt)
  y_pu = Sh_ref / SBase_MVA
  return y_pu
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

  sign = 1.0

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
    if (branch.status == 0)
      @debug "createYBUS: Branch $(branch.comp.cName) not in use, skipping "
      continue
    end

    if branch.isParallel
      @debug "createYBUS: Branch $(branch.comp.cName) is parallel, skipping "
      continue
    end

    r = branch.r_pu
    x = branch.x_pu
    b = branch.b_pu
    g = branch.g_pu

    yik = inv((r + x * im))
    susceptance = g / 2 + b / 2 * im # pi-model

    t = 1.0 + 0.0 * im
    ratio = branch.ratio
    shift_degree = branch.angle
    if ratio != 0.0 || shift_degree != 0.0
      t = calcComplexRatio(ratio, shift_degree)
    end

    Y[fromNode, fromNode] = Y[fromNode, fromNode] + sign * ((yik + susceptance)) / abs2(t)
    Y[toNode, toNode] = Y[toNode, toNode] + sign * (yik + susceptance)
    Y[fromNode, toNode] = -sign * yik / conj(t)
    Y[toNode, fromNode] = -sign * yik / t
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

function adjacentBranches(Y::AbstractMatrix{ComplexF64}, log::Bool = false)
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

function setParallelBranches!(branches::Vector{Branch})
  branchTupleSet = Set{Tuple}()
  branchDict = Dict{Tuple{Integer,Integer},Vector{Branch}}()

  for b in branches
    tupple = (b.fromBus, b.toBus)

    if tupple in branchTupleSet
      b.isParallel = true
      @debug "Branch $(b) is parallel!"

      existing_branches = branchDict[tupple]
      b_total = 0.0
      g_total = 0.0
      r_inv = 0.0
      x_inv = 0.0
      for existing_b in existing_branches
        b_total += existing_b.b_pu
        g_total += existing_b.g_pu
        r_inv += 1 / existing_b.r_pu
        x_inv += 1 / existing_b.x_pu
      end
      for existing_b in existing_branches
        existing_b.b_pu = b_total
        existing_b.g_pu = g_total
        existing_b.r_pu = 1 / r_inv
        existing_b.x_pu = 1 / x_inv
      end
    else
      branchDict[tupple] = [b]
      push!(branchTupleSet, tupple)
    end
  end
end
