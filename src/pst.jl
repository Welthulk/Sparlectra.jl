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
  neutralStep::Int                   # cim:TapChanger.neutralStep, neutral step/position
  step::Int                          # cim:PhaseTapChangerTablePoint.step, actual step/position
  δu::Float64                        # cim:PhaseTapChangerNonLinear.voltageStepIncrement    
  x_0::Float64                       # cim:PhaseTapChangerLinear.xMin or cim:PhaseTapChangerNonLinear.xMin
  x_α_max::Float64                   # cim:PhaseTapChangerLinear.xMax or cim:PhaseTapChangerNonLinear.xMax  
  x_α::Float64                       # cim:TapChangerTablePoint.x  
  α_max_rad::Float64
  α_rad::Float64                     # cim:PhaseTapChangerTablePoint.angle in radians  

  function SymmetricalPhaseShifter(; comp::AbstractComponent, neutralStep::Int, step::Int, δu::Float64, x_0::Float64 = 0.0, x_α_max::Float64, α_max_deg::Float64)
    α_max_rad = α_max_deg * π / 180.0
    obj = new(comp, neutralStep, step, δu, x_0, x_α_max, 0.0, α_max_rad, 0.0)
    setCurrentStep(obj, step)
    return obj
  end
end

function calcCurrentAngle(o::SymmetricalPhaseShifter)
  o.α_rad = 2.0 * atan((o.step - o.neutralStep) * o.δu / 2.0)
end

function calcCurrentX(o::SymmetricalPhaseShifter)
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
  o.step = n
  calcCurrentAngle(o)
  calcCurrentX(o)
end

function getCurrentAngle(o::SymmetricalPhaseShifter)::Float64
  return o.α_rad
end

function getCurrentX(o::SymmetricalPhaseShifter)::Float64
  return o.x_α
end
