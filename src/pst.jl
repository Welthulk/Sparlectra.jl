"""
  mutable struct SymmetricalPhaseShifter <: AbstractBranch

The `SymmetricalPhaseShifter` struct represents a symmetrical phase shifter component.

## Fields
- `comp::AbstractComponent`: The component associated with the phase shifter.
- `neutralStep::Int`: The neutral step or position of the phase shifter.
- `step::Int`: The actual step or position of the phase shifter.
- `δu::Float64`: The voltage step increment of the phase shifter.
- `x_0::Float64`: The minimum value of the phase shifter's position.
- `x_α_max::Float64`: The maximum value of the phase shifter's position.
- `x_α::Float64`: The position of the phase shifter.
- `α_max_rad::Float64`: The maximum angle of the phase shifter in radians.
- `α_rad::Float64`: The angle of the phase shifter in radians.

## Constructors
- `SymmetricalPhaseShifter(; comp::AbstractComponent, neutralStep::Int, step::Int, δu::Float64, x_0::Float64 = 0.0, x_α_max::Float64, α_max_deg::Float64)`: Constructs a new `SymmetricalPhaseShifter` object.

"""

mutable struct SymmetricalPhaseShifter <: AbstractBranch
  comp::AbstractComponent            # component
  sn_mva::Float64                    # cim:PowerTransformer.ratedS, rated power in MVA
  neutralStep::Int                   # cim:TapChanger.neutralStep, neutral step/position
  maxStep::Int                       # maximum step/position 
  step::Int                          # cim:PhaseTapChangerTablePoint.step, actual step/position
  δu::Float64                        # cim:PhaseTapChangerNonLinear.voltageStepIncrement    
  x_0::Float64                       # cim:PhaseTapChangerLinear.xMin or cim:PhaseTapChangerNonLinear.xMin
  x_α_max::Float64                   # cim:PhaseTapChangerLinear.xMax or cim:PhaseTapChangerNonLinear.xMax  
  x_α::Float64                       # cim:TapChangerTablePoint.x  
  α_max_rad::Float64
  α_rad::Float64                     # cim:PhaseTapChangerTablePoint.angle in radians  

  
  function SymmetricalPhaseShifter(; sn_mva::Float64, vn_kV::Float64, from::Int, to::Int, neutralStep::Int, maxStep::Int, step::Int, δu::Float64, x_0::Float64 = 0.0, x_α_max::Float64)
    comp = getPSTImpPGMComp(vn_kV::Float64, from, to)
    obj = new(comp, sn_mva, neutralStep, maxStep, step, δu, x_0, x_α_max, 0.0, 0.0, 0.0)
    calcAlphaMax!(obj)
    setCurrentStep!(obj, step)
    calcCurrentAngle!(obj)
    calcCurrentX!(obj)
    
    return obj
  end

  function Base.show(io::IO, x::SymmetricalPhaseShifter)
    print(io, "SymmetricalPhaseShifter(")
    print(io, "comp=$(x.comp), ")
    print(io, "sn_mva=$(x.sn_mva), ")
    print(io, "neutralStep=$(x.neutralStep), ")
    print(io, "maxStep=$(x.maxStep), ")
    print(io, "step=$(x.step), ")
    print(io, "δu=$(x.δu), ")
    print(io, "x_0=$(x.x_0), ")
    print(io, "x_α_max=$(x.x_α_max), ")
    print(io, "x_α=$(x.x_α), ")
    print(io, "α_max_rad=$(x.α_max_rad), ")
    print(io, "α_rad=$(x.α_rad)")
    print(io, ")")
  end

end

function setCurrentStep!(o::SymmetricalPhaseShifter, step::Int)
  o.step = step
end

function calcAlphaMax!(o::SymmetricalPhaseShifter)
  o.α_max_rad = 2.0 * atan((o.maxStep - o.neutralStep) * o.δu / 2.0)  
end

function calcCurrentAngle!(o::SymmetricalPhaseShifter)
  o.α_rad = 2.0 * atan((o.step - o.neutralStep) * o.δu / 2.0)
end

function calcCurrentX!(o::SymmetricalPhaseShifter)
  o.x_α = o.x_0 + (o.x_α_max - o.x_0) * (sin(o.α_rad / 2.0) / sin(o.α_max_rad / 2.0))^2
end

"""
  setCurrentStep(o::SymmetricalPhaseShifter, n::Int)

Set the current step of a symmetrical phase shifter.

# Arguments
- `o::SymmetricalPhaseShifter`: The symmetrical phase shifter object.
- `n::Int`: The step to set.

"""
function setCurrentStep(o::SymmetricalPhaseShifter, n::Int)
  @assert n >= 0 && n <= o.maxStep
  setCurrentStep!(o,n)  
  calcCurrentAngle!(o)
  calcCurrentX!(o)
end

function getCurrentAngle(o::SymmetricalPhaseShifter)::Float64
  return o.α_rad
end

function getCurrentX(o::SymmetricalPhaseShifter)::Float64
  return o.x_α
end

function calcAdmittance(branch::SymmetricalPhaseShifter, v_kV::Float64, baseMVA::Float64)::Tuple{ComplexF64,ComplexF64,ComplexF64,ComplexF64}
  z_base = (v_kV*v_kV) / baseMVA
  z_ser_pu = getCurrentX(branch)/z_base
  angle = getCurrentAngle(branch)
  y_11 = 1.0 / getCurrentX(z_ser_pu)  
  y_12= -z_ser_pu*exp(-im*angle)
  y_21 = -z_ser_pu*exp(im*angle)  
  y_22 = y_11
  return (y_11, y_12, y_21, y_22)
end

function getPSTImpPGMComp(Vn::Float64, from::Int, to::Int)
  cName = "PST_$(string(round(Vn,digits=1)))"
  cID = "#$cName\\_$from\\_$to"  
  return ImpPGMComp(cID, cName, toComponentTyp("SYMMETRICALPHASESHIFTER"), Vn, from, to)
end

#=


=#