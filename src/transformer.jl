# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 10.05.2023
# include-file transformer.jl

mutable struct TransformerModelParameters
  sn_MVA::Float64 # PGM-Parameter sn in VA
  vk_percent::Float64 # PGM-Parameter sn in percent  
  pk_kW::Union{Nothing,Float64} # PGM-Parameter sn in W
  i0_percent::Float64 # PGM-Parameter sn in A
  p0_kW::Float64 # PGM-Parameter sn in W
  vkr_percent::Union{Nothing,Float64} # real Part of the short-circuit voltage in percent of the rated voltage, nothing = not given

  function TransformerModelParameters(; sn_MVA::Float64, vk_percent::Float64, vkr_percent::Union{Nothing,Float64} = nothing, pk_kW::Union{Nothing,Float64} = nothing, i0_percent::Float64, p0_kW::Float64)
    @assert sn_MVA > 0.0 "sn_mva must be > 0.0"
    @assert vk_percent > 0.0 "uk_percent must be > 0.0"
    @assert !(pk_kW === nothing && ukr_percent === nothing) "At least one of the parameters pk_kw=$(pk_kw) or ukr_percent=$(ukr_percent) must be set"

    new(sn_MVA, vk_percent, pk_kW, i0_percent, p0_kW, vkr_percent)
  end

  function Base.show(io::IO, x::TransformerModelParameters)
    print(io, "TransformerModelParameters(")
    print(io, "sn_MVA=$(x.sn_MVA), ")
    print(io, "vk_percent=$(x.vk_percent), ")
    print(io, "pk_kW=$(x.pk_kW), ")
    print(io, "i0_percent=$(x.i0_percent), ")
    print(io, "p0_kW=$(x.p0_kW)")
    if !isnothing(x.vkr_percent)
      print(io, ", vkr_percent=$(x.vkr_percent)")
    end
    print(io, ")")
  end
end

#helper
#=
function tap_adjust_impedance(tap_pos, tap_min, tap_max, tap_nom, xk, xk_min, xk_max)
    if tap_pos <= max(tap_nom, tap_max) && tap_pos >= min(tap_nom, tap_max)
        if tap_max == tap_nom
            return xk
        end

        xk_increment_per_tap = (xk_max - xk) / (tap_max - tap_nom)
        return xk + (tap_pos - tap_nom) * xk_increment_per_tap
    end

    if tap_min == tap_nom
        return xk
    end

    xk_increment_per_tap = (xk_min - xk) / (tap_min - tap_nom)
    return xk + (tap_pos - tap_nom) * xk_increment_per_tap
end
=#

function calcTransformerRXGB(Vn_kV::Float64, modelData::TransformerModelParameters)::Tuple{Float64,Float64,Float64,Float64}
  z_base = Vn_kV^2 / modelData.sn_MVA
  # Impedanz  
  zk = modelData.vk_percent * 1e-2 * z_base

  # Resistanz
  if !isnothing(modelData.vkr_percent) && modelData.vkr_percent > 0.0
    rk = modelData.vkr_percent * 1e-2 * z_base
  else
    rk = (Vn_kV^2 / modelData.sn_MVA^2) * modelData.pk_kW * 1e-9
  end

  # Reaktanz
  xk = 0.0
  if rk < zk
    xk = sqrt(zk^2 - rk^2)
  else
    @warn "rk >= zk, xk is set to 0.0"
  end

  # Suszeptanz
  if !isnothing(modelData.p0_kW) && modelData.p0_kW > 0.0 && !isnothing(modelData.i0_percent) && modelData.i0_percent > 0.0
    pfe = modelData.p0_kW
    v_quad = Vn_kV^2
    Yfe = pfe / v_quad # realpart
    Yabs = modelData.i0_percent * 1e-2 * modelData.sn_MVA / v_quad # = sqrt(Yfe^2 + Yx^2)
    gm = Yfe
    try
      bm = -1.0 * sqrt(Yabs^2 - Yfe^2)
    catch
      @debug "bm is set to 0.0 (Yabs^2 - Yfe^2 < 0.0)"
      bm = 0.0
    end
  else
    gm = 0.0
    bm = 0.0
  end

  return (rk, xk, gm, bm)
end

mutable struct PowerTransformerTaps
  step::Int                          # actual step/position
  lowStep::Int                       # cim:TapChanger.lowStep
  highStep::Int                      # cim:TapChanger.highStep
  neutralStep::Int                   # cim:TapChanger.neutralStep # 
  voltageIncrement_kV::Float64       # cim:RatioTapChanger.stepVoltageIncrement in kV

  neutralU::Float64                  # cim:TapChanger.neutralU, usually = ratedU of the transformer end, but can deviate
  neutralU_ratio::Float64
  tapStepPercent::Float64
  tapSign::Integer

  function PowerTransformerTaps(; Vn_kV::Float64, step::Int, lowStep::Int, highStep::Int, neutralStep::Int, voltageIncrement_kV::Float64, neutralU::Union{Nothing,Float64} = nothing, neutralU_ratio::Union{Nothing,Float64} = nothing)
    @assert Vn_kV > 0.0 "Vn_kV must be > 0.0"
    @assert highStep != lowStep "tap high position equals low position $(highStep) = $(lowStep)"

    tapStepPercent = 0
    if voltageIncrement_kV != 0.0
      tapStepPercent = (Vn_kV / voltageIncrement_kV) * 1e-2
    end

    #FIXME: calculation of neutralU is not correct
    if isnothing(neutralU)
      if isnothing(neutralU_ratio)
        neutralU_ratio = 1.0
        neutralU = Vn_kV
      end
      neutralU = Vn_kV
      # neutralU = round(neutralU_ratio * Vn_kV + (highStep - lowStep) * tapStepPercent * 1e2, digits = 0)
    end
    tapSign = (highStep > lowStep) ? 1 : -1
    new(step, lowStep, highStep, neutralStep, voltageIncrement_kV, neutralU, neutralU_ratio, tapStepPercent, tapSign)
  end

  function Base.show(io::IO, x::PowerTransformerTaps)
    print(io, "PowerTransformerTaps(")
    print(io, "step=$(x.step), ")
    print(io, "lowStep=$(x.lowStep), ")
    print(io, "hightStep=$(x.highStep), ")
    print(io, "neutralStep=$(x.neutralStep), ")
    print(io, "voltageIncrement_kV=$(x.voltageIncrement_kV), ")
    print(io, "neutralU=$(x.neutralU), ")
    print(io, "neutralU_ratio=$(x.neutralU_ratio), ")
    print(io, "tapStepPercent=$(x.tapStepPercent) ")
    print(io, ")")
  end
end

#= TODO not implemented
function adjustVkDep(Vk::Float64, Vkmax::Float64, Vkmin::Float64, tap::PowerTransformerTaps)
  if isnothing(tap)
    return Vk
  end

   #if tap.neutralStep <= max(tap.neutralStep, tap.highStep) >= min(tap.neutralStep, tap.highStep)

    if tap.step <= max(tap.neutralStep, tap.highStep) && tap.step >= min(tap.neutralStep, tap.highStep)
        if tap_max == tap_nom
            return xk
        end
        cor = (xk_max - xk) / (tap_max - tap_nom)
        return xk + (tap_pos - tap_nom) * xk_increment_per_tap
    end
    if tap_min == tap_nom
        return xk
    end
    xk_increment_per_tap = (xk_min - xk) / (tap_min - tap_nom)
    return xk + (tap_pos - tap_nom) * xk_increment_per_tap
end
=#

mutable struct PowerTransformerWinding
  Vn::Float64                                # cim:PowerTransformerEnd.ratedU in kV
  r::Float64                                 # cim:PowerTransformerEnd.r in Ohm
  x::Float64                                 # cim:PowerTransformerEnd.x in Ohm
  b::Union{Nothing,Float64}                  # cim:PowerTransformerEnd.b in S
  g::Union{Nothing,Float64}                  # cim:PowerTransformerEnd.g in S
  shift_degree::Union{Nothing,Float64}       # cim:PowerTransformerEnd.phaseAngleClock 
  ratedU::Union{Nothing,Float64}             # cim:PowerTransformerEnd.ratedU
  ratedS::Union{Nothing,Float64}             # cim:PowerTransformerEnd.ratedS  
  taps::Union{Nothing,PowerTransformerTaps}
  isPu_RXGB::Union{Nothing,Bool}          # r,x,b in p.u. given? nothing = false
  modelData::Union{Nothing,TransformerModelParameters}
  _isEmpty::Bool

  function PowerTransformerWinding(
    Vn::Float64,
    r::Float64,
    x::Float64,
    b::Union{Nothing,Float64} = nothing,
    g::Union{Nothing,Float64} = nothing,
    shift_degree::Union{Nothing,Float64} = nothing,
    ratedU::Union{Nothing,Float64} = nothing,
    ratedS::Union{Nothing,Float64} = nothing,
    taps::Union{Nothing,PowerTransformerTaps} = nothing,
    isPu_RXGB::Union{Nothing,Bool} = nothing,
    modelData::Union{Nothing,TransformerModelParameters} = nothing,
  )
    new(Vn, r, x, b, g, shift_degree, ratedU, ratedS, taps, isPu_RXGB, modelData, isnothing(modelData))
  end

  function PowerTransformerWinding(;
    Vn_kV::Float64,
    modelData::Union{Nothing,TransformerModelParameters} = nothing,
    shift_degree::Union{Nothing,Float64} = nothing,
    ratedU::Union{Nothing,Float64} = nothing,
    ratedS::Union{Nothing,Float64} = nothing,
    taps::Union{Nothing,PowerTransformerTaps} = nothing,
  )
    @assert Vn_kV > 0.0 "Vn_kv must be > 0.0"
    if isnothing(ratedU)
      ratedU = Vn_kV
    end

    if isnothing(modelData)
      new(Vn_kV, 0.0, 0.0, 0.0, 0.0, shift_degree, ratedU, ratedS, taps, false, nothing, true)  
    else
      r, x, b, g = calcTransformerRXGB(ratedU, modelData)
      new(Vn_kV, r, x, b, g, shift_degree, ratedU, ratedS, taps, false, modelData, false)
    end
  end

  function Base.show(io::IO, x::PowerTransformerWinding)
    print(io, "PowerTransformerWinding(")
    print(io, "Vn=$(x.Vn), ")
    print(io, "r=$(x.r), ")
    print(io, "x=$(x.x), ")
    if !(isnothing(x.b))
      print(io, "b=$(x.b), ")
    end
    if !(isnothing(x.g))
      print(io, "g=$(x.g), ")
    end
    if !(isnothing(x.shift_degree))
      print(io, "shift_degree=$(x.shift_degree), ")
    end
    if !(isnothing(x.ratedU))
      print(io, "ratedU=$(x.ratedU), ")
    end
    if !(isnothing(x.ratedS))
      print(io, "ratedS=$(x.ratedS), ")
    end
    if !(isnothing(x.taps))
      print(io, "taps=$(x.taps), ")
    end
    if !(isnothing(x.isPu_RXGB))
      print(io, "isPu_RXGB=$(x.isPu_RXGB), ")
    end
    if !(isnothing(x.modelData))
      print(io, "$(x.modelData), ")
    end
    print(io, "$(x._isEmpty) ")
    print(io, ")")
  end
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
function recalc_trafo_model_data(; baseMVA::Float64, Sn_MVA::Float64, ratedU_kV::Float64, r_pu::Float64, x_pu::Float64, b_pu::Float64, isPerUnit::Bool)::TransformerModelParameters
  Pfe_kW = 0.0 # -> g_pu = 0.0
  if isPerUnit
    base_power_3p = baseMVA    
    base_power_1p = base_power_3p / 3.0
    base_i_to = base_power_3p / (ratedU_kV * Wurzel3)
    base_y_to = base_i_to * base_i_to / base_power_1p
    base_z_to = 1.0 / base_y_to

    # g_pu = 0.0!
    y_shunt = b_pu
    Y_shunt = y_shunt * base_y_to
    i0 = 100.0 * Y_shunt * ratedU_kV^2 / Sn_MVA

    z_ser_pu = sqrt(r_pu^2 + x_pu^2)
    Z_Series_abs = z_ser_pu * base_z_to
    uk = 100.0*Z_Series_abs * Sn_MVA / ratedU_kV^2
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
  
  return TransformerModelParameters(sn_MVA = Sn_MVA, vk_percent = uk, pk_kW = P_kW, i0_percent = i0, p0_kW = Pfe_kW)
end

function getRXBG(o::PowerTransformerWinding)::Tuple{Float64,Float64,Union{Nothing,Float64},Union{Nothing,Float64}}
  return (o.r, o.x, o.b, o.g)
end
# helper
function hasTaps(x::Union{Nothing,PowerTransformerWinding})::Bool
  if isnothing(x)
    return false
  else
    return !(isnothing(x.taps))
  end
end

function isPerUnit_RXGB(o::PowerTransformerWinding)
  if isnothing(o.isPu_RXGB)
    return false
  else
    return o.isPu_RXGB
  end
end

mutable struct PowerTransformer <: AbstractBranch
  comp::AbstractComponent
  trafoTyp::TrafoTyp
  isControlled::Bool              # cim:TapChanger.controlEnabled 
  nS::Integer                     # Number of (aktive) Sides >=2   
  nController::Integer            # Number of Sides controlled
  isBiWinder::Bool
  HVSideNumber::Integer           # Number of HV Side [1,2,3], for 3WT HV =1, MV=2, LV=3 
  tapSideNumber::Integer          # Number of Tap Side [1,2,3], = 0 if no Tap  
  side1::PowerTransformerWinding
  side2::PowerTransformerWinding
  side3::Union{Nothing,PowerTransformerWinding}
  _equiParms::Integer             # Winding-Side with Short-Circuit-Parameters, 1 for side1,2 for side2, 3 for all sides (3WT)

  function PowerTransformer(comp::AbstractComponent, tapEnable::Bool, s1::PowerTransformerWinding, s2::PowerTransformerWinding, s3::Union{Nothing,PowerTransformerWinding} = nothing, trafoTyp::TrafoTyp = Sparlectra.Ratio)
    equiParms = 0
    if isnothing(s3)
      n = 2
      isBi = true
      if s1.r != 0.0 || s1.x != 0.0
        equiParms = 1
      elseif s2.r != 0.0 || s2.x != 0.0
        equiParms = 2
      end
    else
      isBi = false
      n = 3
      equiParms = 3
    end

    controller = 0
    TapSideNumber = 0
    for x in [s1, s2, s3]
      if hasTaps(x)
        controller += 1
        if controller == 1
          TapSideNumber = x == s1 ? 1 : (x == s2 ? 2 : 3)
        end
      end
    end

    v1 = s1.Vn
    v2 = s2.Vn
    if !isnothing(s3)
      v3 = s3.Vn
    else
      v3 = 0.0
    end
    # find HV Side

    if v1 >= v2 && v1 >= v3
      HVSideNumber = 1
    elseif v2 >= v1 && v2 >= v3
      HVSideNumber = 2
    else
      HVSideNumber = 3
    end

    new(comp, trafoTyp, tapEnable, n, controller, isBi, HVSideNumber, TapSideNumber, s1, s2, s3, equiParms)
  end

  function Base.show(io::IO, x::PowerTransformer)
    print(io, "PowerTransformer(")
    print(io, x.comp)
    print(io, " trafoTyp=", toString(x.trafoTyp), ", ")
    print(io, "isControlled=$(x.isControlled), ")
    print(io, "ns=$(x.nS), ")
    print(io, "nController=$(x.nController), ")
    print(io, "isBiWinder=$(x.isBiWinder), ")
    print(io, "HVSideNumber=$(x.HVSideNumber), ")
    print(io, "tapSideNumber=$(x.tapSideNumber), ")
    print(io, "equiParms=$(x._equiParms), ")
    println(io, "") # next line
    println(io, "side1=$(x.side1), ")
    if (isnothing(x.side3))
      print(io, "side2=$(x.side2) ")
    else
      println(io, "side2=$(x.side2), ")
      print(io, "side3=$(x.side3) ")
    end

    print(io, ")")
  end
end

function getWinding2WT(x::PowerTransformer)
  @assert x.isBiWinder "Transformer is not a 2WT"
  return !x.side1._isEmpty ? x.side1 : x.side2
end

function getSideNumber2WT(x::PowerTransformer)
  @assert x.isBiWinder "Transformer is not a 2WT"
  return !x.side1._isEmpty ? 1 : 2
end

function getVn2WT(x::PowerTransformer)
  side = getSideNumber2WT(x)
  vn = (side == 1) ? x.side1.Vn : x.side2.Vn
  return vn
end

function getTapStepPercent(x::PowerTransformer)
  @assert x.isBiWinder "Transformer is not a 2WT"
  tap = x.tapSideNumber == 1 ? x.side1.taps : (x.tapSideNumber == 2 ? x.side2.taps : nothing)
  if isnothing(tap)
    return (false, 0.0)
  end
  vn = getVn2WT(x)
  return (true, (tap.voltageIncrement_kV / vn) * 100.0)
end

function isTapInNeutralPosition(x::PowerTransformer)
  @assert x.isBiWinder "Transformer is not a 2WT"
  tap = x.tapSideNumber == 1 ? x.side1.taps : (x.tapSideNumber == 2 ? x.side2.taps : nothing)
  if isnothing(tap)
    return true
  else
    return (tap.step == tap.neutralStep) ? true : false
  end
end

function calcTransformerRatio(x::PowerTransformer)
  @assert x.isBiWinder "Transformer is not a 2WT"
  @assert !isnothing(x.side1.Vn) "Vn must be set for winding 1"
  @assert !isnothing(x.side2.Vn) "Vn must be set for winding 2"
  @assert !isnothing(x.side1.ratedU) "ratedU must be set for winding 1"
  @assert !isnothing(x.side2.ratedU) "ratedU must be set for winding 2"

  v1 = x.side1.Vn
  v2 = x.side2.Vn
  r1 = x.side1.ratedU
  r2 = x.side2.ratedU
  ratio = (v1 * r2) / (v2 * r1)
  taps = x.tapSideNumber == 1 ? x.side1.taps : (x.tapSideNumber == 2 ? x.side2.taps : nothing)
  if isTapInNeutralPosition(x)
    return ratio
  else
    stat, tapStepPercent = getTapStepPercent(x)
    if stat
      tapStepPercent = (taps.voltageIncrement_kV / getVn2WT(x)) * 100.0
      corr = 1 + (taps.step - taps.neutralStep) * tapStepPercent / 100
      return ratio / corr
    else
      @warn "no taps, but tapSideNumber != 0!"
      return ratio
    end
  end
end

function toString(o::TrafoTyp)::String
  if o == Ratio
    return "Ratio"
  elseif o == PhaseShifter
    return "PhaseShifter"
  elseif o == PhaseTapChanger
    return "PhaseTapChanger"
  else
    return "UnknownT"
  end
end

function create2WTRatioTransformerNoTaps(;from::Int, to::Int, vn_hv_kv::Float64, vn_lv_kv::Float64, sn_mva::Float64,vk_percent::Float64, vkr_percent::Float64, pfe_kw::Float64, i0_percent::Float64)::PowerTransformer
  
  c = getTrafoImpPGMComp(false, vn_hv_kv, from, to)

  modelData = TransformerModelParameters(sn_MVA = sn_mva, vk_percent = vk_percent, vkr_percent = vkr_percent, pk_kW = pfe_kw, i0_percent = i0_percent, p0_kW = 0.0)
  w1 = PowerTransformerWinding(Vn_kV = vn_hv_kv, modelData = modelData)
  w2 = PowerTransformerWinding(Vn_kV = vn_lv_kv, modelData = modelData)
  
  trafo = PowerTransformer(c, false, w1, w2, nothing, Sparlectra.Ratio)  
  return trafo
  
end

# purpose: creates Windings for 3WT using MVA Method, see: "MVA Method for Three-Winding Transformer" by Ver Pangonilo
#
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
#

#ajustVkTap = 1.0 # FIXME: implement adjustVkTap!
function create3WTWindings!(; u_kV::Array{Float64,1}, sn_MVA::Array{Float64,1}, addEx_Side::Array{TransformerModelParameters,1}, sh_deg::Array{Float64,1}, tap_side::Int, tap::PowerTransformerTaps)::Tuple{PowerTransformerWinding,PowerTransformerWinding,PowerTransformerWinding}
  for side = 1:3
    if isnothing(addEx_Side[side].vkr_percent)
      @assert !isnothing(addEx_Side[side].pk_kW) "Either pk_kW or vkr_percent must be set!"
      # Rk = (U^2/S^2) * Pk; Rk = Ukr * U^2/S * 1e-2 => Ukr = Pk * 100 / S
      # pk -> kW, S -> MVA; kW/MVA = 1e-3 => faktor = 1e-3*100 = 0.1
      addEx_Side[side].vkr_percent = addEx_Side[side].pk_kW * 0.1 / addEx_Side[side].sn_MVA
    end
  end

  vkVec = []
  vkrVec = []
  for side = 1:3
    fak = addEx_Side[1].sn_MVA / (min(addEx_Side[side].sn_MVA, addEx_Side[(side%3)+1].sn_MVA))
    push!(vkVec, (addEx_Side[side].vk_percent * fak))
    push!(vkrVec, (addEx_Side[side].vkr_percent * fak))
  end

  vkVec[1] = 0.5 * (vkVec[1] + vkVec[3] - vkVec[2]) * 1.0
  vkVec[2] = 0.5 * (vkVec[2] + vkVec[1] - vkVec[3]) * addEx_Side[2].sn_MVA / addEx_Side[1].sn_MVA
  vkVec[3] = 0.5 * (vkVec[2] + vkVec[3] - vkVec[1]) * addEx_Side[3].sn_MVA / addEx_Side[1].sn_MVA

  vkrVec[1] = 0.5 * (vkrVec[1] + vkrVec[3] - vkVec[2]) * 1.0
  vkrVec[2] = 0.5 * (vkrVec[2] + vkrVec[1] - vkVec[3]) * addEx_Side[2].sn_MVA / addEx_Side[1].sn_MVA
  vkrVec[3] = 0.5 * (vkrVec[2] + vkrVec[3] - vkVec[1]) * addEx_Side[3].sn_MVA / addEx_Side[1].sn_MVA

  pVec = []
  for side = 1:3
    push!(pVec, calcTransformerRXGB(u_kV[side], addEx_Side[side]))
  end

  wVec = []
  for side = 1:3
    tap = (side - 1 == tap_side) ? tap : nothing
    w = PowerTransformerWinding(u_kV[side], pVec[side][1], pVec[side][2], pVec[side][3], pVec[side][4], sh_deg[side], u_kV[side], sn_MVA[side], tap, false, addEx_Side[side])
    push!(wVec, w)
  end

  return (wVec[1], wVec[2], wVec[3])
end

function getWT2BusID(Vn::Float64, from::Int, to::Int)
  name = "2WT_$(string(convert(Int,trunc(Vn))))"
  id = "#$name\\_$from\\_$to"
  return name, id
end

function getWT3BusID(Vn::Float64, from::Int, to::Int, to3::Int)
  name = "3WT_$(string(convert(Int,trunc(Vn))))"
  id = "#$name\\_$from\\_$to\\_$to3"
  return name, id
end

# in the case of parallel transformers, the `to3` connections should always be distinct
function getWT3AuxBusID(Vn::Float64, from::Int, to::Int, to3::Int)
  name = "3WT_Aux_$(string(convert(Int, trunc(Vn))))"
  id = "#$name\\_$from\\_$to\\_$to3"
  return name, id
end

function getTrafoImpPGMComp(aux::Bool, Vn::Float64, from::Int, to::Int, to3::Union{Nothing,Int} = nothing)
  cName, cID = aux ? getWT3AuxBusID(Vn, from, to, to3) : isnothing(to3) ? getWT2BusID(Vn, from, to) : getWT3BusID(Vn, from, to, to3)
  return aux ? ImpPGMComp3WT(cID, cName, toComponentTyp("POWERTRANSFORMER"), Vn, from, to, to3) : ImpPGMComp(cID, cName, toComponentTyp("POWERTRANSFORMER"), Vn, from, to)
end
