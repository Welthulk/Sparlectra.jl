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
