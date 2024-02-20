# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 10.05.2023
# include-file component.jl

# Data type to describe a component
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
end

@enum SeitenTyp UnknownS = 0 Seite1 = 1 Seite2 = 2 Seite3 = 3

@enum TrafoTyp UnknownT = 0 Ratio = 1 PhaseShifter = 2 PhaseTapChanger = 3

@enum NodeType UnknownN = 0 PQ = 1 PV = 2 Slack = 3 Isolated = 4
@enum ProSumptionType UnknownP = 0 Injection = 1 Consumption = 2

abstract type AbstractComponent end
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
  else
    return UnknownC
  end
end
