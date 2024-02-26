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

  function TransformerModelParameters(; sn_MVA::Float64, vk_percent::Float64, vkr_percent::Union{Nothing,Float64}, pk_kW::Union{Nothing,Float64}, i0_percent::Float64, p0_kW::Float64)
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
function calcTransformerRXGB(Vn_kV::Float64, modelData::TransformerModelParameters)::Tuple{Float64,Float64,Float64,Float64}
  z_base = Vn_kV^2 / modelData.sn_MVA
  # Impedanz
  zk = modelData.vk_percent * 1e-2 * z_base
  # Resistanz
  if !isnothing(modelData.vkr_percent) && modelData.vkr_percent > 0.0
    rk = modelData.vkr_percent * 1e-2 * z_base
  else    
    rk = (Vn_kV^2/modelData.sn_MVA^2)*modelData.pk_kW    
  end
  # Reaktanz
  xk = sqrt(zk^2 - rk^2) # Reaktanz
  # Suszeptanz
  if !isnothing(modelData.p0_kW) && modelData.p0_kW > 0.0 && !isnothing(modelData.i0_percent) && modelData.i0_percent > 0.0
    pfe = modelData.p0_kW
    v_quad = Vn_kV^2
    Yfe = pfe / v_quad
    Y0 = modelData.i0_percent / 100.0 * modelData.sn_MVA / v_quad
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

  function PowerTransformerTaps(; Vn_kV::Float64, step::Int, lowStep::Int, highStep::Int, neutralStep::Int, voltageIncrement_kV::Float64, neutralU::Union{Nothing,Float64}=nothing, neutralU_ratio::Union{Nothing,Float64}=nothing)
    @assert Vn_kV > 0.0 "Vn_kV must be > 0.0"
    @assert voltageIncrement_kV > 0.0 "voltageIncrement must be > 0.0"

    tapStepPercent = (Vn_kV / voltageIncrement_kV) * 1e-2

    #FIXME: calculation of neutralU is not correct
    if isnothing(neutralU)
      if isnothing(neutralU_ratio)
        neutralU_ratio = 1.0
        neutralU = Vn_kV
      end
      neutralU = Vn_kV
    # neutralU = round(neutralU_ratio * Vn_kV + (highStep - lowStep) * tapStepPercent * 1e2, digits = 0)
    end

    new(step, lowStep, highStep, neutralStep, voltageIncrement_kV, neutralU, neutralU_ratio, tapStepPercent)
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

  #=
  function PowerTransformerWinding(;
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
    new(Vn, r, x, b, g, shift_degree, ratedU, ratedS, taps, isPu_RXGB, modelData)
  end
  =#
  function PowerTransformerWinding(;
    Vn_kV::Float64,
    modelData::Union{Nothing,TransformerModelParameters} = nothing,
    shift_degree::Union{Nothing,Float64} = nothing,
    ratedU::Union{Nothing,Float64} = nothing,
    ratedS::Union{Nothing,Float64} = nothing,
    taps::Union{Nothing,PowerTransformerTaps} = nothing,
  )
    @assert Vn_kV > 0.0 "Vn_kv must be > 0.0"
    if isnothing(modelData)
      new(Vn_kV, 0.0, 0.0, nothing, nothing, nothing, nothing, nothing, nothing, false, nothing, true)
    else
      if isnothing(ratedU)
        ratedU = Vn_kV
      end
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

function getRXBG(o::PowerTransformerWinding)::Tuple{Float64,Float64,Float64,Float64}
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

mutable struct PowerTransformer
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

  function PowerTransformer(comp::AbstractComponent, tapEnable::Bool, s1::PowerTransformerWinding, s2::PowerTransformerWinding, s3::Union{Nothing,PowerTransformerWinding} = nothing, trafoTyp::TrafoTyp = ResDataTypes.Ratio)
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

  if !x.side1._isEmpty
    return x.side1
  else
    return x.side2
  end
end

#TODO: implement
#function getRatio2WT(x::PowerTransformer, bus_HV_Vn::Float64, bus_LV_Vn::Float64) 
#   @assert x.isBiWinder "Transformer is not a 2WT"
  
#end
# helper

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
