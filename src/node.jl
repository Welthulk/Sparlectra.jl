# Copyright 2023–2025 Udo Schmitz
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 10.05.2023
# include-file nodes.jl

# Data type to describe the topology
"""
    Node

A mutable structure representing a node in a power system.

# Fields
- `comp::AbstractComponent`: The component of the node.
- `busIdx::Integer`: The index of the bus.
- `_nodeType::NodeType`: The type of the node.
- `_ratedS::Union{Nothing,Float64}`: The rated power of the node.
- `_lZone::Union{Nothing,Integer}`: The loss zone of the node.
- `_area::Union{Nothing,Integer}`: The area of the node.
- `_vm_pu::Union{Nothing,Float64}`: The voltage magnitude of the node in per unit.
- `_va_deg::Union{Nothing,Float64}`: The voltage angle of the node in degrees.
- `_pƩLoad::Union{Nothing,Float64}`: The total active power load at the node.
- `_qƩLoad::Union{Nothing,Float64}`: The total reactive power load at the node.
- `_pShunt::Union{Nothing,Float64}`: The total active power shunt at the node.
- `_qShunt::Union{Nothing,Float64}`: The total reactive power shunt at the node.
- `_pƩGen::Union{Nothing,Float64}`: The total active power generation at the node.
- `_qƩGen::Union{Nothing,Float64}`: The total reactive power generation at the node.
- `_vmin_pu::Union{Nothing,Float64}`: The minimum voltage magnitude at the node in per unit.
- `_vmax_pu::Union{Nothing,Float64}`: The maximum voltage magnitude at the node in per unit.

# Constructors
- `Node(;
    busIdx::Integer,
    vn_kV::Float64,
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
    vmin_pu::Union{Nothing,Float64} = nothing,
    vmax_pu::Union{Nothing,Float64} = nothing,
    isAux::Bool = false,
    oBusIdx::Union{Nothing,Int} = nothing,
  )`: Creates a new `Node` instance.

# Methods
- `Base.show(io::IO, node::Node)`: Prints the `Node` instance.
"""
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
  _vmin_pu::Union{Nothing,Float64} # Minimum voltage magnitude in p.u
  _vmax_pu::Union{Nothing,Float64} # Maximum voltage magnitude in p.u

  function Node(;
    busIdx::Integer,
    vn_kV::Float64,
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
    vmin_pu::Union{Nothing,Float64} = nothing,
    vmax_pu::Union{Nothing,Float64} = nothing,
    isAux::Bool = false,
    oBusIdx::Union{Nothing,Int} = nothing,
  )
    bIdx = busIdx
    if !isnothing(oBusIdx)
      bIdx = oBusIdx
    end
    c = getNodeComp(vn_kV, bIdx, nodeType, isAux)

    new(c, busIdx, nodeType, ratedS, zone, area, vm_pu, va_deg, pƩLoad, qƩLoad, pShunt, qShunt, pƩGen, qƩGen, vmin_pu, vmax_pu)
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
    if (!isnothing(node._vmin_pu))
      print(io, "Vmin: ", node._vmin_pu, ", ")
    end
    if (!isnothing(node._vmax_pu))
      print(io, "Vmax: ", node._vmax_pu, ", ")
    end
    println(io, ")")
  end
end

function getNodeComp(Vn_kV::Float64, node_idx::Int, nodeType, isAux::Bool = false)::ImpPGMComp
  cTyp = toComponentTyp("Busbarsection")
  if (nodeType == Slack)
    name = "Bus_$(Int(node_idx))_$(string(convert(Int,trunc(Vn_kV))))*"
  elseif isAux
    name = "Aux_$(Int(node_idx))_$(string(convert(Int,trunc(Vn_kV))))"
  else
    name = "Bus_$(Int(node_idx))_$(string(convert(Int,trunc(Vn_kV))))"
  end
  cID = "#" * name * "#"
  return ImpPGMComp(cID, name, cTyp, Vn_kV, node_idx, node_idx)
end

function addShuntPower!(; node::Node, p::Float64, q::Float64)
  if !isnothing(p)
    if (isnothing(node._pShunt))
      node._pShunt = 0.0
    end
    node._pShunt = node._pShunt + p
  end

  if !isnothing(q)
    if (isnothing(node._qShunt))
      node._qShunt = 0.0
    end
    node._qShunt = node._qShunt + q
  end
end

function addLoadPower!(; node::Node, p::Union{Nothing,Float64}, q::Union{Nothing,Float64})
  if !isnothing(p)
    if (isnothing(node._pƩLoad))
      node._pƩLoad = 0.0
    end
    node._pƩLoad = node._pƩLoad + p
  end
  if !isnothing(q)
    if (isnothing(node._qƩLoad))
      node._qƩLoad = 0.0
    end
    node._qƩLoad = node._qƩLoad + q
  end
end

function addGenPower!(; node::Node, p::Union{Nothing,Float64}, q::Union{Nothing,Float64})
  if !isnothing(p)
    if (isnothing(node._pƩGen))
      node._pƩGen = 0.0
    end
    node._pƩGen = node._pƩGen + p
  end

  if !isnothing(q)
    if (isnothing(node._qƩGen))
      node._qƩGen = 0.0
    end
    node._qƩGen = node._qƩGen + q
  end
end

function busComparison(node1::Node, node2::Node)
  node1.busIdx < node2.busIdx
end

function setVmVa!(; node::Node, vm_pu::Float64, va_deg::Union{Nothing,Float64} = nothing)
  node._vm_pu = vm_pu
  if !isnothing(va_deg)
    node._va_deg = va_deg
  end
end

function isSlack(o::Node)
  if o._nodeType == Slack
    return true
  else
    return false
  end
end

function getNodeVm(o::Node)::Float64
  return o._vm_pu
end

function getNodeVn(o::Node)::Float64
  return o.comp.cVN
end

function isPQNode(o::Node)
  if o._nodeType == PQ
    return true
  else
    return false
  end
end

function isPVNode(o::Node)
  if o._nodeType == PV || o._nodeType == Slack
    return true
  else
    return false
  end
end

function getNodeType(o::Node)::NodeType
  return o._nodeType
end

function isIsolated(o::Node)
  if o._nodeType == Isolated
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

function setNodeType!(o::Node, typ::String)
  o._nodeType = toNodeType(typ)
end
