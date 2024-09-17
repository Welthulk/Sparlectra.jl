# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 10.05.2023
# include-file component.jl

@enum ComponentTyp begin
  UnknownC
  NodeC
  LineC
  Trafo
  Busbarsection
  Generator
  Load
  ExternalNetworkInjection
  AsynchronousMachine
  SynchronousMachine
  EnergyConsumer
  LinearShuntCompensator
  StaticVarCompensator
  AuxBus
  BranchC
  SymmetricalPhaseShifterC
end

@enum TrafoTyp UnknownT = 0 Ratio = 1 PhaseShifter = 2 PhaseTapChanger = 3 PIModel = 4

@enum NodeType UnknownN = 0 PQ = 1 PV = 2 Slack = 3 Isolated = 4

@enum ProSumptionType UnknownP = 0 Injection = 1 Consumption = 2

abstract type AbstractComponent end
"""
    Component

A structure representing a component in a power system.

# Fields
- `cID::String`: The ID of the component.
- `cName::String`: The name of the component.
- `cTyp::ComponentTyp`: The type of the component.
- `cVN::Float64`: The nominal voltage of the component in kV.

# Constructors
- `Component(id::String, name::String, typ::ComponentTyp, vn::Float64)`: Creates a new `Component` instance with a specified nominal voltage.
- `Component(id::String, name::String, typ::ComponentTyp)`: Creates a new `Component` instance with a nominal voltage of 1 kV.
- `Component(id, name, cmp::String, Unenn::Float64)`: Creates a new `Component` instance with a specified nominal voltage and component type as a string.

# Methods
- `Base.show(io::IO, x::Component)`: Prints the `Component` instance.

# Example
```julia
Component("1", "Generator", "GENERATOR", 110.0)
```
"""
struct Component <: AbstractComponent
  cID::String
  cName::String
  cTyp::ComponentTyp
  cVN::Float64           # nominal voltage in kV

  function Component(id::String, name::String, typ::ComponentTyp, vn::Float64)
    new(id, name, typ, vn)
  end

  function Component(id::String, name::String, typ::ComponentTyp)
    new(id, name, typ, 1)
  end

  function Component(id, name, cmp::String, Unenn::Float64)
    typ = toComponentTyp(cmp)
    new(id, name, typ, Unenn)
  end

  function Base.show(io::IO, x::Component)
    print(io, "Component(")
    print(io, "ID=$(x.cID), ")
    print(io, "Name=$(x.cName), ")
    print(io, "Typ=$(x.cTyp), ")
    print(io, "Vn=$(x.cVN)")
    print(io, ")")
  end
end

struct ImpPGMComp <: AbstractComponent
  cID::String
  cName::String
  cTyp::ComponentTyp
  cVN::Float64           # nominal voltage in kV  
  cFrom_bus::Union{Nothing,Int64}
  cTo_bus::Union{Nothing,Int64}

  function ImpPGMComp(id::String, name::String, typ::ComponentTyp, vn::Float64)
    new(id, name, typ, vn, nothing, nothing)
  end

  function ImpPGMComp(id::String, name::String, typ::ComponentTyp, vn::Float64, fBus::Int64, tBus::Int64)
    new(id, name, typ, vn, fBus, tBus)
  end

  function ImpPGMComp(cmp::Component, fBus::Int64, tBus::Int64, newID::Union{Nothing,String} = nothing)
    if isnothing(newID)
      new(cmp.cID, cmp.cName, cmp.cTyp, cmp.cVN, fBus, tBus)
    else
      new(newID, cmp.cName, cmp.cTyp, cmp.cVN, fBus, tBus)
    end
  end

  function Base.show(io::IO, x::ImpPGMComp)
    print(io, "Component(")
    print(io, "ID=$(x.cID), ")
    print(io, "Name=$(x.cName), ")
    print(io, "Typ=$(x.cTyp), ")
    print(io, "Vn=$(x.cVN), ")
    if !isnothing(x.cFrom_bus)
      print(io, "From_bus=$(x.cFrom_bus), ")
    end
    if !isnothing(x.cTo_bus)
      print(io, "To_bus=$(x.cTo_bus), ")
    end
    print(io, ")")
  end
end

struct ImpPGMComp3WT <: AbstractComponent
  cID::String
  cName::String
  cTyp::ComponentTyp
  cVN::Float64           # nominal voltage in kV  
  cHV_bus::Int64
  cMV_bus::Int64
  cLV_bus::Int64

  function ImpPGMComp3WT(id::String, name::String, typ::ComponentTyp, vn::Float64, hvBus::Int64, mvBus::Int64, lvBus::Int64)
    new(id, name, typ, vn, hvBus, mvBus, lvBus)
  end

  function Base.show(io::IO, x::ImpPGMComp3WT)
    print(io, "Component(")
    print(io, "ID=$(x.cID), ")
    print(io, "Name=$(x.cName), ")
    print(io, "Typ=$(x.cTyp), ")
    print(io, "Vn=$(x.cVN), ")
    print(io, "HV_bus=$(x.cHV_bus), ")
    print(io, "MV_bus=$(x.cMV_bus), ")
    print(io, "LV_bus=$(x.cLV_bus) ")
    print(io, ")")
  end
end

# helper
function toComponentTyp(o::String)::ComponentTyp
  val = uppercase(o)
  if val == "ACLINESEGMENT"
    return LineC
  elseif val == "BUS" || val == "BUSBARSECTION"
    return Busbarsection
  elseif val == "POWERTRANSFORMER" || val == "TRANSFORMER"
    return Trafo
  elseif val == "ASYNCHRONOUSMACHINE"
    return AsynchronousMachine
  elseif val == "SYNCHRONOUSMACHINE"
    return SynchronousMachine
  elseif val == "GENERATOR"
    return Generator
  elseif val == "LOAD"
    return Load
  elseif val == "EXTERNALNETWORKINJECTION"
    return ExternalNetworkInjection
  elseif val == "ENERGYCONSUMER"
    return EnergyConsumer
  elseif val == "SHUNT"
    return LinearShuntCompensator
  elseif val == "LINEARSHUNTCOMPENSATOR"
    return LinearShuntCompensator
  elseif val == "STATICVARCOMPENSATOR"
    return StaticVarCompensator
  elseif val == "NODE"
    return NodeC
  elseif val == "AUXBUS"
    return AuxBus
  elseif val == "BRANCH"
    return BranchC
  elseif val == "SYMMETRICALPHASESHIFTER"  
    return SymmetricalPhaseShifterC
  else
    return UnknownC
  end
end

function getCompName(c::AbstractComponent)::String
  return c.cName
end

function getCompID(c::AbstractComponent)::String
  return c.cID
end