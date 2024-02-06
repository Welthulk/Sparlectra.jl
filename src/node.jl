# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 10.05.2023
# include-file nodes.jl

# Data type to describe the topology
mutable struct Node
  comp::AbstractComponent
  terminals::Vector{Terminal}
  busIdx::Integer
  _kidx::Integer # original busnumber
  _nodeType::NodeType
  _auxNodeID::Union{Nothing,String} # auxiliary node ID for mapping to node

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

  function Node(
    c::Component,
    t::Vector{Terminal},
    busIdx::Integer,
    nodeType::NodeType,
    auxNodeID::Union{Nothing,String} = nothing,
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
  )
    new(c, t, busIdx, busIdx, nodeType, auxNodeID, ratedS, zone, area, vm_pu, va_deg, pƩLoad, qƩLoad, pShunt, qShunt, pƩGen, qƩGen)
  end

  function Node(
    c::Component,
    t::Vector{Terminal},
    busIdx::Integer,
    kIdx::Integer,
    nodeType::NodeType,
    auxNodeID::Union{Nothing,String} = nothing,
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
  )
    new(c, t, busIdx, kIdx, nodeType, auxNodeID, ratedS, zone, area, vm_pu, va_deg, pƩLoad, qƩLoad, pShunt, qShunt, pƩGen, qƩGen)
  end

  function Node(
    id::String,
    name::String,
    vn::Float64,
    terminals::Vector{Terminal},
    auxNodeID::Union{Nothing,String} = nothing,
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
  )
    nodeType = ResDataTypes.Busbarsection
    c = Component(id, name, nodeType, vn)

    for term in terminals
      if c.cVN != term.comp.cVN
        @warn "Voltage levels of terminals are not equal Vn = $(c.cVN) != $(term.comp.cVN), Terminal: $(term.comp.cName), Node-Name: $(name)"
      end
    end

    t = copy(terminals)
    sort!(t, by = x -> x.seite)
    new(c, t, 0, 0, ResDataTypes.UnknownN, auxNodeID, ratedS, zone, area, vm_pu, va_deg, pƩLoad, qƩLoad, pShunt, qShunt, pƩGen, qƩGen)
  end

  function Base.show(io::IO, node::Node)
    print(io, "Node( ")
    print(io, "ID: ", node.comp.cID, ", ")
    print(io, "Name: ", node.comp.cName, ", ")
    print(io, "Typ: ", node.comp.cTyp, ", ")
    print(io, "Vn: ", node.comp.cVN, ", ")
    print(io, "BusIndex: ", node.busIdx, ", ")
    print(io, "NodeIndex: ", node._kidx, ", ")
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
    for t in node.terminals
      println(io, t)
    end
  end
end

#helper 
function getNodeName(node::Node)::String
  return node.comp.cName
end

function getNodeID(node::Node)::String
  return node.comp.cID
end

function getNodeVn(node::Node)::Float64
  return node.comp.cVN
end

function getNodeType(node::Node)::NodeType
  return node._nodeType
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

function setNodeIdx!(node::Node, kIdx::Integer)
  node._kidx = kIdx
end

function setBusIdx!(node::Node, bIdx::Integer)
  node.busIdx = bIdx
end

function setNodeType!(node::Node, nodeType::NodeType)
  node._nodeType = nodeType
end

function setNodeType!(node::Node, type::String)
  nodeType = toNodeType(type)
  node._nodeType = nodeType
end

function setRatedS!(node::Node, ratedS::Float64)
  node._ratedS = ratedS
end

function setNodePQ!(node::Node, p::Union{Nothing,Float64}, q::Union{Nothing,Float64})
  if !isnothing(p)
    node._pƩLoad = p
  end
  if !isnothing(q)
    node._qƩLoad = q
  end
end

function setShuntPower!(node::Node, p::Union{Nothing,Float64}, q::Union{Nothing,Float64})
  if !isnothing(p)
    node._pShunt = p
  end

  if !isnothing(q)
    node._qShunt = q
  end
end

function setGenPower!(node::Node, p::Union{Nothing,Float64}, q::Union{Nothing,Float64})
  if !isnothing(p)
    node._pƩGen = p
  end

  if !isnothing(q)
    node._qƩGen = q
  end
end

function addAktivePower!(node::Node, p::Float64)
  if (isnothing(node._pƩLoad))
    node._pƩLoad = 0.0
  end
  node._pƩLoad = node._pƩLoad + p
end

function addReaktivePower!(node::Node, q::Float64)
  if (isnothing(node._qƩLoad))
    node._qƩLoad = 0.0
  end
  node._qƩLoad = node._qƩLoad + q
end

function addGenAktivePower!(node::Node, p::Float64)
  if (isnothing(node._pƩGen))
    node._pƩGen = 0.0
  end
  node._pƩGen = node._pƩGen + p
end

function addGenReaktivePower!(node::Node, q::Float64)
  if (isnothing(node._qƩGen))
    node._qƩGen = 0.0
  end
  node._qƩGen = node._qƩGen + q
end

function hasPowerInjection(node::Node)
  if (isnothing(node._pƩLoad) && isnothing(node._qƩLoad))
    return false
  end
  return true
end

function hasShuntInjection(node::Node)
  if (isnothing(node._pShunt) && isnothing(node._qShunt))
    return false
  end
  return true
end

# helper for sorting
function nodeComparison(node1::Node, node2::Node)
  node1._kidx < node2._kidx
end

function busComparison(node1::Node, node2::Node)
  node1.busIdx < node2.busIdx
end

function setLossZone!(node::Node, lZone::Integer)
  node._lZone = lZone
end

function setArea!(node::Node, area::Integer)
  node._area = area
end

function setVm!(node::Node, Vm::Float64)
  node._Vm_pu = Vm
end

function setVa!(node::Node, Va::Float64)
  node._Va_deg = Va
end

function setVmVa!(node::Node, vm_pu::Float64, va_deg::Float64)
  node._vm_pu = vm_pu
  node._va_deg = va_deg
end
