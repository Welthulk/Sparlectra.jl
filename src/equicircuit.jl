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

function calc_y_pu(Y::ComplexF64, baseS::Float64, baseV::Float64)
  return Y * baseS / baseV^2
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
purpose: calculation of change of volatage per tap position in percent
+vs: change per tap
+vn: nominal voltage
"""
function calcTapStepPercent(vs::Float64, vn::Float64)
  return (vs / vn) * 100.0
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
function calcRatio(bus_hv_kv::Float64, vn_hv_kv::Float64, bus_lv_kv::Float64, vn_lv_kv::Float64, tapPos::Integer, tapNeutral::Integer, tapStepPercent::Float64)::Float64
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
purpose: calculate of neutral voltage. Tap position in neutral position.
+ sbase_VA: base power in VA
+ u2_V: rated voltage in V (to_side, low voltage side)
+ uk: short circuit voltage in percent
+ sn_VA: rated power in VA
+ pk_W: short circuit power in W (Cu-losses)
+ i0_A: relativ no load current
"""
function calcTrafoParamsSI(sbase_VA::Float64, u2_V::Float64, uk::Float64, sn_VA::Float64, pk_W::Float64, i0::Float64, p0_W::Float64)::Tuple{Float64,Float64,Float64,Float64}
  return calcTrafoParams(sn_mva = sn_VA * 1e-6, vn_hv_kv = u2_V * 1e-3, vk_percent = uk * 100.0, pk_kw = pk_W * 1e-3, pfe_kw = p0_W * 1e-3, i0_percent = i0 * 100.0)
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
 purpose: calculates R, X, rfe and xmu of a transformer
 + sn_mva: rated power in MVA
 + vn_hv_kv: rated voltage high voltage side in kV
 + vk_percent: short circuit voltage in percent
 + vkr_percent: real part short circuit voltage in percent  
 + pfe_kw: iron losses in kW
 + i0_percent: no load current in percent (open circuit)
 
"""
function calcTrafoParams(; sn_mva::Float64, vn_hv_kv::Float64, vk_percent::Float64, pk_kw::Union{Nothing,Float64} = nothing, vkr_percent::Union{Nothing,Float64} = nothing, pfe_kw::Union{Nothing,Float64} = nothing, i0_percent::Union{Nothing,Float64} = nothing)
  @assert sn_mva > 0.0 "sn_mva must be > 0.0"
  @assert vn_hv_kv > 0.0 "vn_hv_kv must be > 0.0"
  @assert vk_percent > 0.0 "vk_percent must be > 0.0"
  @assert !(pk_kw === nothing && vkr_percent === nothing) "At least one of the parameters pk_kw=$(pk_kw) or vkr_percent=$(vkr_percent) must be set"

  z_base = vn_hv_kv^2 / sn_mva
  # Impedanz
  zk = vk_percent / 100.0 * z_base
  # Resistanz
  if !isnothing(vkr_percent)
    rk = vkr_percent / 100.0 * z_base
  else
    rk = (pk_kw / sn_mva) / z_base
  end
  # Reaktanz
  xk = sqrt(zk^2 - rk^2) # Reaktanz
  # Suszeptanz
  if !isnothing(pfe_kw) && pfe_kw > 0.0 && !isnothing(i0_percent) && i0_percent > 0.0
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

function calcGB_Shunt(p_shunt::Float64, q_shunt::Float64, vn_kv::Float64)
  g = p_shunt / vn_kv^2
  b = q_shunt / vn_kv^2
  return g, b
end

function calcPQ_Shunt(g1::Float64, b1::Float64, vn_kv::Float64)
  p_shunt = g1 * vn_kv^2
  q_shunt = b1 * vn_kv^2
  return p_shunt, q_shunt
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
    elseif branch.skipYBus
      @debug "createYBUS: Branch $(branch) skipping "
      continue
    elseif branch.isParallel
      @debug "createYBUS: Branch $(branch) is parallel, using adjustes parameters "
      r = branch.adjRXGB.r_pu
      x = branch.adjRXGB.x_pu
      b = branch.adjRXGB.b_pu
      g = branch.adjRXGB.g_pu
    else
      r = branch.r_pu
      x = branch.x_pu
      b = branch.b_pu
      g = branch.g_pu
    end

    yik = inv((r + x * im))
    susceptance = g / 2 + b / 2 * im # pi-model

    t = 1.0 + 0.0 * im
    ratio = branch.ratio
    shift_degree = branch.angle
    if ratio != 0.0 || shift_degree != 0.0
      t = calcComplexRatio(ratio, shift_degree)
    end

    Y[toNode, toNode] = Y[toNode, toNode] +  (yik + susceptance)
    Y[fromNode, fromNode] = Y[fromNode, fromNode] + ((yik + susceptance)) / abs2(t)
    
    Y[fromNode, toNode] = -1.0 * yik / conj(t)
    Y[toNode, fromNode] = -1.0* yik / t
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
      existing_branches = branchDict[tupple]
      push!(existing_branches, b)
    else
      branchDict[tupple] = [b]
      push!(branchTupleSet, tupple)
    end
  end

  for (k, b_vec) in branchDict
    if length(b_vec) > 1
      sum_b_pu = 0.0
      sum_g_pu = 0.0
      sum_z = 0.0

      for b in b_vec
        b.isParallel = true
        b.skipYBus = true
        if b.status == 1          
          sum_b_pu += b.b_pu
          sum_g_pu += b.g_pu
          sum_z += (b.r_pu - b.x_pu * im) / (b.r_pu^2 + b.x_pu^2)
        end
      end

      z_total = 1.0 / sum_z
      r_pu = real(z_total)
      x_pu = imag(z_total)

      last_b = b_vec[end]
      for (i,b) in enumerate(b_vec)
        if b.status == 1
          adjP = AdjElecParams(r_pu = r_pu, x_pu = x_pu, b_pu = sum_b_pu, g_pu = sum_g_pu)
          setAdjElecParam!(adjP, b)
          @debug "branch (2): $(b)"
          if i==1
            @debug "last parallel element (3): $(b)"
            b.skipYBus = false
          end
        end
      end
    end
  end
end
