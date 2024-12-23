# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 10.05.2023
# include-file prosumer.jl

# Data type to describe producers and consumers
"""
    ProSumer

A mutable structure representing a prosumer in a power system. A prosumer is an entity that either produces or consumes power.

# Fields
- `comp::AbstractComponent`: The component of the prosumer.
- `busIdx::Int`: The index of the bus where the prosumer is connected.
- `pGen::Float64`: The active power generation of the prosumer.
- `qGen::Float64`: The reactive power generation of the prosumer.
- `pLoad::Float64`: The active power consumption of the prosumer.
- `qLoad::Float64`: The reactive power consumption of the prosumer.
- `status::Int`: The status of the prosumer. 1 = in service, 0 = out of service.

# Constructors
- `ProSumer(comp::AbstractComponent, busIdx::Int, pGen::Float64, qGen::Float64, pLoad::Float64, qLoad::Float64, status::Int)`: Creates a new `ProSumer` instance.

# Methods
- `Base.show(io::IO, prosumer::ProSumer)`: Prints the `ProSumer` instance.

# Example
```julia
prosumer = ProSumer(comp, 1, 100.0, 50.0, 80.0, 40.0, 1)
"""
mutable struct ProSumer
  comp::AbstractComponent
  ratedS::Union{Nothing,Float64}
  ratedU::Union{Nothing,Float64}
  qPercent::Union{Nothing,Float64}
  pVal::Union{Nothing,Float64}      # aktive power
  qVal::Union{Nothing,Float64}      # reactive power
  maxP::Union{Nothing,Float64}
  minP::Union{Nothing,Float64}
  maxQ::Union{Nothing,Float64}
  minQ::Union{Nothing,Float64}
  ratedPowerFactor::Union{Nothing,Float64}
  referencePri::Union{Nothing,Integer}
  vm_pu::Union{Nothing,Float64}
  va_deg::Union{Nothing,Float64}
  proSumptionType::ProSumptionType
  isAPUNode::Bool
  qGenRepl::Union{Nothing,Float64}

  # @enum NodeType UnknownN=0 PQ=1 PV=2 Slack=3
  # @enum ProSumptionType UnknownP=0 Injection=1 Consumption=2 

  function ProSumer(;
    vn_kv::Float64,
    busIdx::Int,
    oID::Int,
    type::ProSumptionType,
    ratedS::Union{Nothing,Float64} = nothing,
    ratedU::Union{Nothing,Float64} = nothing,
    qPercent::Union{Nothing,Float64} = nothing,
    p::Union{Nothing,Float64} = nothing,
    q::Union{Nothing,Float64} = nothing,
    maxP::Union{Nothing,Float64} = nothing,
    minP::Union{Nothing,Float64} = nothing,
    maxQ::Union{Nothing,Float64} = nothing,
    minQ::Union{Nothing,Float64} = nothing,
    ratedPowerFactor::Union{Nothing,Float64} = nothing,
    referencePri::Union{Nothing,Integer} = nothing,
    vm_pu::Union{Nothing,Float64} = nothing,
    va_deg::Union{Nothing,Float64} = nothing,
    isAPUNode::Bool = false,
  )
    comp = getProSumPGMComp(vn_kv, busIdx, isGenerator(type), oID)

    if isnothing(vm_pu)
      vm_pu = 1.0
    end

    if isnothing(va_deg)
      va_deg = 0.0
    end

    new(comp, ratedS, ratedU, qPercent, p, q, maxP, minP, maxQ, minQ, ratedPowerFactor, referencePri, vm_pu, va_deg, type, isAPUNode, nothing)
  end

  function Base.show(io::IO, prosumption::ProSumer)
    print(io, "ProSumption( ")
    print(io, prosumption.comp, ", ")

    if (!isnothing(prosumption.ratedS))
      print(io, "ratedS: ", prosumption.ratedS, " MVA, ")
    end

    if (!isnothing(prosumption.ratedU))
      print(io, "ratedU: ", prosumption.ratedU, " kV, ")
    end

    if (!isnothing(prosumption.qPercent))
      print(io, "qPercent: ", prosumption.qPercent, " %, ")
    end

    if (!isnothing(prosumption.pVal))
      print(io, "pVal: ", prosumption.pVal, " MW, ")
    end

    if (!isnothing(prosumption.qVal))
      print(io, "qVal: ", prosumption.qVal, " MVar, ")
    end

    if (!isnothing(prosumption.maxP))
      print(io, "maxP: ", prosumption.maxP, " MW, ")
    end

    if (!isnothing(prosumption.minP))
      print(io, "minP: ", prosumption.minP, " MW, ")
    end

    if (!isnothing(prosumption.maxQ))
      print(io, "maxQ: ", prosumption.maxQ, " MVar, ")
    end

    if (!isnothing(prosumption.minQ))
      print(io, "minQ: ", prosumption.minQ, " MVar, ")
    end

    if (!isnothing(prosumption.ratedPowerFactor))
      print(io, "ratedPowerFactor: ", prosumption.ratedPowerFactor, ", ")
    end

    if (!isnothing(prosumption.referencePri))
      print(io, "referencePri: ", prosumption.referencePri, ", ")
    end

    if (!isnothing(prosumption.vm_pu))
      print(io, "vm_pu: ", prosumption.vm_pu, ", ")
    end

    if (!isnothing(prosumption.va_deg))
      print(io, "va_deg: ", prosumption.va_deg, ", ")
    end

    if (!isnothing(prosumption.qGenRepl))
      print(io, "qGenRepl: ", prosumption.qGenRepl, ", ")
    end

    print(io, "proSumptionType: ", prosumption.proSumptionType, ")")
  end
end

function isAPUNode(o::ProSumer)
  return o.isAPUNode
end

function isGenerator(c::Sparlectra.ProSumptionType)::Bool
  if c == Injection
    return true
  else
    return false
  end
end # isGenerator

function isGenerator(o::ProSumer)::Bool
  return isGenerator(o.proSumptionType)
end # isGenerator

function getProSumPGMComp(Vn::Float64, from::Int, isGen::Bool, id::Int)
  cTyp = isGen ? toComponentTyp("GENERATOR") : toComponentTyp("LOAD")
  cName = isGen ? "Gen_$(string(convert(Int,trunc(Vn))))" : "Ld_$(string(round(Vn,digits=1)))"
  cID = "#$cName\\_$from\\_#$(string(id))"
  return ImpPGMComp(cID, cName, cTyp, Vn, from, from)
end

function setQGenReplacement!(o::ProSumer, q::Float64)
  o.qGenRepl = q
end

function getQGenReplacement(o::ProSumer)::Union{Nothing,Float64}
  return o.qGenRepl
end

function updatePQ!(o::ProSumer, p::Union{Nothing,Float64}, q::Union{Nothing,Float64})
  if !isnothing(p)
    o.pVal = p
  end
  if !isnothing(q)
    o.qVal = q
  end
end

function isSlack(o::ProSumer)
  if !isnothing(o.referencePri) && o.referencePri > 0
    return true
  else
    return false
  end
end

function toProSumptionType(o::ComponentTyp)::ProSumptionType
  if o == Generator || o == ExternalNetworkInjection || o == SynchronousMachine
    return Injection
  elseif o == AsynchronousMachine || o == Load || o == EnergyConsumer
    return Consumption
  else
    return UnknownP
  end
end

function toProSumptionType(o::String)::ProSumptionType
  val = uppercase(o)
  if val == "GENERATOR" || val == "EXTERNALNETWORKINJECTION" || val == "SYNCHRONOUSMACHINE"
    return Injection
  elseif val == "ASYNCHRONOUSMACHINE" || val == "LOAD" || val == "ENERGYCONSUMER"
    return Consumption
  else
    return UnknownP
  end
end

# helper
function toString(o::ProSumptionType)::String
  if o == Injection
    return "Injection"
  elseif o == Consumption
    return "Consumption"
  else
    return "UnknownP"
  end
end
