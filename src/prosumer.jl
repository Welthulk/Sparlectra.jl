# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 10.05.2023
# include-file prosumer.jl

# Data type to describe producers and consumers
mutable struct ProSumer
  comp::AbstractComponent
  nodeID::Union{Nothing,String}  # ID for mapping to node, could set later
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
  vm_degree::Union{Nothing,Float64}
  proSumptionType::ProSumptionType
  isAPUNode::Bool

  # @enum NodeType UnknownN=0 PQ=1 PV=2 Slack=3
  # @enum ProSumptionType UnknownP=0 Injection=1 Consumption=2 

  function ProSumer(
    equipment::Component,
    nodeID::String,
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
    vm_degree::Union{Nothing,Float64} = nothing,
  )

    proSumptionType = toProSumptionType(equipment.cTyp)
    isAPUNode = isPUNode(equipment.cTyp)

    if isnothing(vm_pu)
      vm_pu = 1.0
    end

    if isnothing(vm_degree)
      vm_degree = 0.0
    end

    new(equipment, nodeID, ratedS, ratedU, qPercent, p, q, maxP, minP, maxQ, minQ, ratedPowerFactor, referencePri, vm_pu, vm_degree, proSumptionType, isAPUNode)
  end

  function Base.show(io::IO, prosumption::ProSumer)
    print(io, "ProSumption( ")
    print(io, prosumption.comp, ", ")
    print(io, "NodeID: ", prosumption.nodeID, ", ")


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

    if (!isnothing(prosumption.vm_degree))
      print(io, "vm_degree: ", prosumption.vm_degree, ", ")
    end

    print(io, "proSumptionType: ", prosumption.proSumptionType, ")")
  end

end

