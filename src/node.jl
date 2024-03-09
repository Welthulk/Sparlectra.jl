# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 10.05.2023
# include-file nodes.jl

# Data type to describe the topology
mutable struct Node
  comp::AbstractComponent  
  busIdx::Integer  
  _nodeType::NodeType
  _ratedS::Union{Nothing,Float64}
  _lZone::Union{Nothing,Integer} # loss Zone
  _area::Union{Nothing,Integer}

  _vm_pu::Union{Nothing,Float64} # Voltage magnitude in p.u
  _va_deg::Union{Nothing,Float64} # Voltage angle in degree

  _pƩLoad::Union{Nothing,Float64}  # Ʃ active power
  _qƩLoad::Union{Nothing,Float64}  # Ʃ reactive power
  _pShunt::Union{Nothing,Float64} # Ʃ active power
  _qShunt::Union{Nothing,Float64} # Ʃ reactive power
  _pƩGen::Union{Nothing,Float64}  # Ʃ active power injected
  _qƩGen::Union{Nothing,Float64}  # Ʃ reactive power injected

  function Node(;  
    busIdx::Integer,    
    Vn_kV::Float64,
    nodeType::NodeType,        
    ratedS::Union{Nothing,Float64} = nothing,
    zone::Union{Nothing,Integer} = nothing,
    area::Union{Nothing,Integer} = nothing,
    vm_pu::Union{Nothing,Float64} = nothing,
    va_deg::Union{Nothing,Float64} = nothing,
    pƩLoad::Union{Nothing,Float64} = nothing,
    qƩLoad::Union{Nothing,Float64} = nothing,
    pShunt::Union{Nothing,Float64} = nothing,
    qShunt::Union{Nothing,Float64} = nothing,
    pƩGen::Union{Nothing,Float64} = nothing,
    qƩGen::Union{Nothing,Float64} = nothing,
    isAux::Bool = false,
  )
    c = getNocdeComp(Vn_kV, busIdx, nodeType, isAux)
    new(c, busIdx, nodeType, ratedS, zone, area, vm_pu, va_deg, pƩLoad, qƩLoad, pShunt, qShunt, pƩGen, qƩGen)
  end


  function Base.show(io::IO, node::Node)
    print(io, "Node( ")
    print(io, node.comp, ", ")
    print(io, "BusIndex: ", node.busIdx, ", ")
    
    print(io, "nodeType: ", node._nodeType, ", ")
    if (!isnothing(node._lZone))
      print(io, "Loss Zone: ", node._lZone, ", ")
    end
    if (!isnothing(node._area))
      print(io, "Area: ", node._area, ", ")
    end
    if (!isnothing(node._vm_pu))
      print(io, "Vm: ", node._vm_pu, ", ")
    end
    if (!isnothing(node._va_deg))
      print(io, "Va: ", node._va_deg, ", ")
    end
    if (!isnothing(node._ratedS))
      print(io, "RatedS: ", node._ratedS, ", ")
    end
    if (!isnothing(node._pƩLoad))
      print(io, "ƩPLoad: ", node._pƩLoad, ", ")
    end
    if (!isnothing(node._qƩLoad))
      print(io, "ƩQLoad: ", node._qƩLoad, ", ")
    end
    if (!isnothing(node._pShunt))
      print(io, "pShunt: ", node._pShunt, ", ")
    end
    if (!isnothing(node._qShunt))
      print(io, "qShunt: ", node._qShunt, ", ")
    end
    if (!isnothing(node._pƩGen))
      print(io, "ƩPGen: ", node._pƩGen, ", ")
    end
    if (!isnothing(node._qƩGen))
      print(io, "ƩQGen: ", node._qƩGen, ", ")
    end
    println(io, ")")
  end
end

function getNocdeComp(Vn_kV::Float64, node_idx::Int, nodeType, isAux::Bool = false)::ImpPGMComp  
  cTyp = toComponentTyp("Busbarsection")
  if (nodeType == Slack)
    name = "Bus_$(Int(node_idx))_$(string(convert(Int,trunc(Vn_kV))))*"
  elseif isAux
    name = "Aux_$(Int(node_idx))_$(string(convert(Int,trunc(Vn_kV))))"
  else
    name = "Bus_$(Int(node_idx))_$(string(convert(Int,trunc(Vn_kV))))"
  end
  cID = "#"*name*"#"
  return ImpPGMComp(cID, name, cTyp, Vn_kV, node_idx, node_idx)
end

mutable struct NodeParameters
  kIdx::Integer # Knoten Index
  ratedS::Union{Nothing,Float64}
  vm_pu::Union{Nothing,Float64}   # Voltage magnitude in p.u
  va_deg::Union{Nothing,Float64}  # Voltage angle in degree  
  pƩLoad::Union{Nothing,Float64}  # Ʃ active power
  qƩLoad::Union{Nothing,Float64}  # Ʃ reactive power
  pShunt::Union{Nothing,Float64}  # Ʃ active power
  qShunt::Union{Nothing,Float64}  # Ʃ reactive power
  pƩGen::Union{Nothing,Float64}   # Ʃ active power injected
  qƩGen::Union{Nothing,Float64}   # Ʃ reactive power injected
  function NodeParameters(
    kIdx::Integer,
    ratedS::Union{Nothing,Float64} = nothing,
    vm_pu::Union{Nothing,Float64} = nothing,
    va_deg::Union{Nothing,Float64} = nothing,
    pƩLoad::Union{Nothing,Float64} = nothing,
    qƩLoad::Union{Nothing,Float64} = nothing,
    pShunt::Union{Nothing,Float64} = nothing,
    qShunt::Union{Nothing,Float64} = nothing,
    pƩGen::Union{Nothing,Float64} = nothing,
    qƩGen::Union{Nothing,Float64} = nothing,
  )
    new(kIdx, ratedS, vm_pu, va_deg, pƩLoad, qƩLoad, pShunt, qShunt, pƩGen, qƩGen)
  end

  function Base.show(io::IO, nodeParam::NodeParameters)
    print(io, "NodeParameters( ")
    print(io, "kIdx: ", nodeParam.kIdx, ", ")
    if (!isnothing(nodeParam.ratedS))
      print(io, "RatedS: ", nodeParam.ratedS, ", ")
    end
    if (!isnothing(nodeParam.vm_pu))
      print(io, "Vm: ", nodeParam.vm_pu, ", ")
    end
    if (!isnothing(nodeParam.va_deg))
      print(io, "Va: ", nodeParam.va_deg, ", ")
    end
    if (!isnothing(nodeParam.pƩLoad))
      print(io, "ƩP: ", nodeParam.pƩLoad, ", ")
    end
    if (!isnothing(nodeParam.qƩLoad))
      print(io, "ƩQ: ", nodeParam.qƩLoad, ", ")
    end
    if (!isnothing(nodeParam.pShunt))
      print(io, "pShunt: ", nodeParam.pShunt, ", ")
    end
    if (!isnothing(nodeParam.qShunt))
      print(io, "qShunt: ", nodeParam.qShunt, ", ")
    end
    if (!isnothing(nodeParam.pƩGen))
      print(io, "ƩPGen: ", nodeParam.pƩGen, ", ")
    end
    if (!isnothing(nodeParam.qƩGen))
      print(io, "ƩQGen: ", nodeParam.qƩGen, ", ")
    end
    println(io, ")")
  end
end

"""
Purpose: Set node parameters
"""
function setNodeParameters!(node::Node, nodeParam::NodeParameters)
  if (!isnothing(nodeParam.ratedS))
    node._ratedS = nodeParam.ratedS
  end

  if (!isnothing(nodeParam.vm_pu))
    node._vm_pu = nodeParam.vm_pu
  end

  if (!isnothing(nodeParam.va_deg))
    node._va_deg = nodeParam.va_deg
  end

  if (!isnothing(nodeParam.pƩLoad))
    node._pƩLoad = nodeParam.pƩLoad
  end

  if (!isnothing(nodeParam.qƩLoad))
    node._qƩLoad = nodeParam.qƩLoad
  end

  if (!isnothing(nodeParam.pShunt))
    node._pShunt = nodeParam.pShunt
  end

  if (!isnothing(nodeParam.qShunt))
    node._qShunt = nodeParam.qShunt
  end

  if (!isnothing(nodeParam.pƩGen))
    node._pƩGen = nodeParam.pƩGen
  end

  if (!isnothing(nodeParam.qƩGen))
    node._qƩGen = nodeParam.qƩGen
  end
end

function addShuntPower!(;node::Node, p::Float64, q::Float64)
  if (isnothing(node._pShunt))
    node._pShunt = 0.0
  end
  node._pShunt = node._pShunt + p

  if (isnothing(node._qShunt))
    node._qShunt = 0.0
  end
  node._qShunt = node._qShunt + q  
end

function addLoadPower!(;node::Node, p::Float64, q::Float64)
  if (isnothing(node._pƩLoad))
    node._pƩLoad = 0.0
  end
  node._pƩLoad = node._pƩLoad + p

  if (isnothing(node._qƩLoad))
    node._qƩLoad = 0.0
  end
  node._qƩLoad = node._qƩLoad + q
end

function addGenPower!(;node::Node, p::Float64, q::Float64)
  if (isnothing(node._pƩGen))
    node._pƩGen = 0.0
  end
  node._pƩGen = node._pƩGen + p

  if (isnothing(node._qƩGen))
    node._qƩGen = 0.0
  end
  node._qƩGen = node._qƩGen + q
end

function busComparison(node1::Node, node2::Node)
  node1.busIdx < node2.busIdx
end

function setVmVa!(;node::Node, vm_pu::Float64, va_deg::Float64)
  node._vm_pu = vm_pu
  node._va_deg = va_deg
end

function isSlack(o::ResDataTypes.Node)
  if o._nodeType == ResDataTypes.Slack
    return true
  else
    return false
  end  
end

function toNodeType(o::Int)::NodeType
  if o == 1
    return PQ
  elseif o == 2
    return PV
  elseif o == 3
    return Slack
  elseif o == 4
    return Isolated
  else
    return UnknownN
  end
end

function toNodeType(o::String)::NodeType
  val = uppercase(o)
  if val == "PQ"
    return PQ
  elseif val == "PV"
    return PV
  elseif val == "SLACK"
    return Slack
  elseif val == "ISOLATED"
    return Isolated
  else
    return UnknownN
  end
end

function toString(o::NodeType)::String
  if o == PQ
    return "PQ"
  elseif o == PV
    return "PV"
  elseif o == Slack
    return "SLACK"
  elseif o == Isolated
    return "ISOLATED"
  else
    return "UNKNOWN"
  end
end