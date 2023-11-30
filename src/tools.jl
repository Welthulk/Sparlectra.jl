# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 23.6.2023
# include-file SparlectraTools.jl


function isGenerator(c::ResDataTypes.Component)
  if c.cTyp == ResDataTypes.Generator || c.cTyp == ResDataTypes.SynchronousMachine
    return true
  else
    return false
  end
end # isGenerator

function isMotor(c::ResDataTypes.Component)
  if c.cTyp == ResDataTypes.AsynchronousMachine
    return true
  else
    return false
  end
end # isMotor

function isLoad(c::ResDataTypes.Component)
  if c.cTyp == ResDataTypes.Load || c.cTyp == ResDataTypes.ExternalNetworkInjection || c.cTyp == ResDataTypes.EnergyConsumer
    return true
  else
    return false
  end
end # isLoad

function isExternalNetworkInjection(c::ResDataTypes.Component)
  if c.cTyp == ResDataTypes.ExternalNetworkInjection
    return true
  else
    return false
  end
end # isExternalNetworkInjection

function isShunt(c::ResDataTypes.Component)
  if c.cTyp == ResDataTypes.LinearShuntCompensator || c.cTyp == ResDataTypes.StaticVarCompensator
    return true
  else
    return false
  end
end

roundUpToNearest100(number) = ceil(number / 100) * 100
"""
  purpose: reports the crusial data for summarizing the network
"""
function reportCrusialData(acseg::ResDataTypes.ACLineSegment)
  @info "Line-Name: $(acseg.comp.cName), ID: $(acseg.comp.cID) , Voltage: $(acseg.comp.cVN), Length: $(acseg.length)"
end #   

function reportCrusialData(trafo::ResDataTypes.PowerTransformer)
  typ = ResDataTypes.toString(trafo.trafoTyp)
  @info "Trafo-Name: $(trafo.comp.cName), ID: $(trafo.comp.cID), TR: $(trafo.side1.Vn)/$(trafo.side2.Vn), Typ: $(typ), Controlled: $(trafo.isControlled), BiWinder: $(trafo.isBiWinder)"
end # reportCrusialData

function reportCrusialData(node::ResDataTypes.Node)

  anzZW = 0
  for _ in node.terminals
    anzZW += 1
  end

  if !isnothing(node._ratedS)
    s = roundUpToNearest100(node._ratedS)
  else
    s = "***"
  end
  @info "Node-Name: $(node.comp.cName), ID: $(node.comp.cID), Voltage: $(node.comp.cVN), S: $(s) MVA, Type: $(node._nodeType), Terminals: $(anzZW)"
end # reportCrusialData

function reportCrusialData(o::ResDataTypes.ProSumer)
  faktor = 1.0
  sber = 1.0
  ber = false

  pqavail = false

  s = "***"

  pf = "***"

  p = "***"

  q = "***"


  if !isnothing(o.pVal)
    p = o.pVal
  end

  if !isnothing(o.qVal)
    q = o.qVal
  end

  if p != "***" && q != "***"
    pqavail = true
  end

  if !isnothing(o.ratedS)
    s = roundUpToNearest100(o.ratedS)
  elseif pqavail
    sber = sqrt(p^2 + q^2)
    ber = true
    s = round(sber, digits = 2)
  end

  if !isnothing(o.ratedPowerFactor)
    pf = round(o.ratedPowerFactor, digits = 2)
  elseif ber
    pf = round(p / sber, digits = 2)
  end

  typ = ResDataTypes.toString(o.proSumptionType)
  name = o.comp.cName
  id = o.comp.cID
  V = o.comp.cVN
  art = o.comp.cTyp
  if !ber
    @info "ProSumption-Name: $(name), ID: $(id), Voltage: $(V), Typ: $(art), Node-Typ: $(typ),  S: $(s) MVA, fak: $(pf), p: $(p), q: $(q)"
  else
    @info "ProSumption-Name: $(name), ID: $(id), Voltage: $(V), Typ: $(art), Node-Typ: $(typ),  S: $(s) MVA (calc), fak: $(pf) (calc), p: $(p), q: $(q)"
  end
end # reportCrusialData

function getTrafoById(id::String, trafos::Vector{ResDataTypes.PowerTransformer})
  for trafo in trafos
    if trafo.comp.cID == id
      return trafo
    end
  end
end # GetTrafoById

function getTrafoByName(name::String, trafos::Vector{ResDataTypes.PowerTransformer})
  for trafo in trafos
    if trafo.comp.cName == name
      return trafo
    end
  end
end # GetTrafoByName   

function getNodeByName(name::String, nodes::Vector{ResDataTypes.Node})
  for node in nodes
    if node.comp.cName == name
      return node
    end
  end
end # GetNodeByName

function getTerminalsByName(name::String, terminals::Vector{ResDataTypes.Terminal})
  for terminal in terminals
    if terminal.comp.cName == name
      return terminal
    end
  end
end # getTerminalsByName

#depricated
function getNodesAndBranches(nodes::Vector{Node})

  KnotenVec = Vector{Int32}()
  _ZweigVec = Vector{String}()
  i = 0 # Knoten
  for node in nodes
    i += 1
    push!(KnotenVec, i)
    terminals = node.Terminals

    for term in terminals
      if term.comp.cTyp == ResDataTypes.LineC || term.comp.cTyp == ResDataTypes.Trafo
        push!(_ZweigVec, term.comp.Name)
      end
    end
  end
  ZweigVec = unique(_ZweigVec)
  return KnotenVec, ZweigVec
end

"""
Purpose: Checks uniqueness of line ids
"""
function checkLineIdIsUnique(asec::Vector{ResDataTypes.ACLineSegment})
  for l in asec
    id = l.comp.cID
    for l2 in asec
      if l2.comp.cID == id && l2 != l
        @warn "line id $id is not unique"
        return false
      end
    end
  end
  return true
end # checkLineIdIsUnique

"""
  Purpose: Checks uniqueness of transformer ids
"""
function checkTransformerIdIsUnique(trafo::Vector{ResDataTypes.PowerTransformer})
  for t in trafo
    id = t.comp.cID
    for t2 in trafo
      if t2.comp.cID == id && t2 != t
        @warn "transformer id $id is not unique"
        return false
      end
    end
  end
  return true
end # checkTransformerIdIsUnique

"""
  Purpose: Checks uniqueness of prosumer ids
"""
function checkProSumerIdIsUnique(prosumer::Vector{ResDataTypes.ProSumer})
  result = true
  for p in prosumer
    id = p.comp.cID
    for p2 in prosumer
      if p2.comp.cID == id && p2 != p
        @warn "prosumer id $id is not unique"
        result = false
      end
    end
  end
  return result
end # checkProSumerIdIsUnique

"""
  Purpose: Checks uniqueness of node ids
"""
function checkNodeIdsUnique(nodes::Vector{ResDataTypes.Node})::Bool
  for n in nodes
    id = n.comp.cID
    for n2 in nodes
      if n2.comp.cID == id && n2 != n
        @warn "node id $id is not unique"
        return false
      end
    end
  end
  return true
end

function checkNodeNamesUnique(nodes::Vector{ResDataTypes.Node})::Bool
  for n in nodes
    name = n.comp.cName
    for n2 in nodes
      if n2.comp.cName == name && n2 != n
        @warn "node name $name is not unique"
        return false
      end
    end
  end
  return true
end

function searchTerminalByName(name::String, terminals::Vector{ResDataTypes.Terminal})
  for term in terminals
    if term.comp.cName == name
      return term
    end
  end
  return nothing
end

function searchTerminalById(id::String, terminals::Vector{ResDataTypes.Terminal})
  for term in terminals
    if term.comp.cID == id
      return term
    end
  end
  return nothing
end

function searchNodeByName(name::String, nodes::Vector{ResDataTypes.Node})
  for node in nodes
    if node.comp.cName == name
      return node
    end
  end
  return nothing
end

function searchTerminalByNameAndSide(name::String, nodes::Vector{ResDataTypes.Node}, side::ResDataTypes.SeitenTyp)
  for node in nodes
    for term in node.terminals
      if term.comp.cName == name
        if term.seite == side
          return term
        end
      end
    end
  end  
  return nothing
end

function searchNodeById(id::String, nodes::Vector{ResDataTypes.Node})
  for node in nodes
    if node.comp.cID == id
      return node
    end
  end
  return nothing
end

"""
  Purpose: Search for a node by terminal id and side
"""
function searchNodeByTerminalId(id::String, nodes::Vector{ResDataTypes.Node}, side::ResDataTypes.SeitenTyp)
  for node in nodes
    for term in node.terminals
      if term.comp.cID == id
        if term.seite == side
          return node, term
        end
      end
    end
  end

  return nothing, nothing
end

"""
  Purpose: Checks if the given terminal has two sides
"""
function checkTerminalHasTwoSides(cID::String, nodes::Vector{ResDataTypes.Node})::Bool

  if isnothing(nodes)
    @warn "nodes is nothing"
    return false
  end

  if isnothing(cID)
    @warn "cID is nothing"
    return false
  end

  n1, t1 = searchNodeByTerminalId(cID, nodes, ResDataTypes.Seite1)
  if isnothing(n1)
    @warn "first terminal $(cID) not found!"
    return false
  end

  n2, t2 = searchNodeByTerminalId(cID, nodes, ResDataTypes.Seite2)
  if isnothing(n2)
    @warn "second terminal $(cID) not found!"
    return false
  end

  return true
end # checkTerminalHasTwoSides

"""
  Purpose: Checks if the given terminal has two sides
"""
function checkTerminalOnlyOnce(cID::String, nodes::Vector{ResDataTypes.Node})::Bool
  if isnothing(nodes)
    @warn "nodes is nothing"
    return false
  end

  if isnothing(cID)
    @warn "cID is nothing"
    return false
  end

  side1 = true
  side2 = true

  n1, t1 = searchNodeByTerminalId(cID, nodes, ResDataTypes.Seite1)
  if isnothing(n1) && isnothing(t1)
    side1 = false
  end
  n2, t2 = searchNodeByTerminalId(cID, nodes, ResDataTypes.Seite2)
  if isnothing(n2) && isnothing(t2)
    side2 = false
  end

  if !side1 && !side2
    @warn "terminal $(cID) is not connected to any node!"
    return false
  end

  if side1 & side2
    @warn "terminal $(cID) is connected to two nodes!"
    return false
  end

  return true
end # checkTerminalOnlyOnce

function checkNodeConnections(nodeVec::Vector{ResDataTypes.Node})
  result = true
  for node in nodeVec
    for term in node.terminals
      typ = term.comp.cTyp

      if typ == ResDataTypes.LineC || typ == ResDataTypes.Trafo

        id = term.comp.cID
        # Check Net Trail
        ok1 = checkTerminalHasTwoSides(id, nodeVec)
        if !ok1        
          @warn "Terminal $id has not two sides"
          result = false
        end
      elseif typ == ResDataTypes.Generator ||
             typ == ResDataTypes.Load ||
             typ == ResDataTypes.ExternalNetworkInjection ||
             typ == ResDataTypes.SynchronousMachine ||
             typ == ResDataTypes.AsynchronousMachine ||
             typ == ResDataTypes.EnergyConsumer ||
             typ == ResDataTypes.LinearShuntCompensator ||
             typ == ResDataTypes.StaticVarCompensator

        id = term.comp.cID
        ok2 = checkTerminalOnlyOnce(id, nodeVec)
        if !ok2
          @warn "Terminal $id has more than one sides"
          result = false
        end
      end # if typ
    end # for term   
  end # for node
  return result
end # function

"""
purpose  : Function to identify net tails, nodes with only two terminals (BusbarSection+other) and no power injection or consumption
date     : 8.8.23
"""
function findNetTrails(branches::Vector{ResDataTypes.Branch})
  connected_buses = Set{Integer}()

  for branch in branches
    push!(connected_buses, branch.fromBus)
    push!(connected_buses, branch.toBus)
  end

  all_buses = Set{Integer}([branch.fromBus for branch in branches] |> unique |> collect)

  network_outliers = setdiff(all_buses, connected_buses)

  return network_outliers
end
"""
purpose  : Function to identify net tails, branches with only one terminal
date     : 8.8.23
"""
function findSingleConnectedBuses(branches::Vector{ResDataTypes.Branch}, nodes::Vector{ResDataTypes.Node})
  for n in nodes
    if n._nodeType == ResDataTypes.PQ
      if length(n.Terminals) == 2
        # set node isolated if no power injection
        if !(hasPowerInjection(n) || hasShuntInjection(n))
          bus = n.busIdx
          @info "node ", n.comp.cID, ", Bus ", bus, " is isolated"
          setNodeType!(n, ResDataTypes.Isolated)
          # search connected branch and set out of service
          for b in branches
            if b.toBus == bus
              @info "branch (1)", b.branchC.ID, ", from ", b.fromBus, ", toBus", b.toBus, " is out of service"
              b.status = 0
              break
            end
          end
        end
      end
    end
  end
end

"""
purpose: removes isolated buses, change fromNode and toNodes indexes
"""
function renumberNodeIdx!(branchVec::Vector{Branch}, nodes::Vector{ResDataTypes.Node}, log::Bool = false)
  kIdxBusIdxMap = Dict{Int,Int}()
  idx = 0

  for n in nodes
    if n._nodeType != ResDataTypes.Isolated
      idx += 1
      n._kidx = idx
      n.busIdx = idx
      kIdxBusIdxMap[n.busIdx] = idx
    else
      n.busIdx = 0
      n._kidx = 0
    end
  end
  #@show kIdxBusIdxMap  

  for b in branchVec
    if haskey(kIdxBusIdxMap, b.fromBus) # fromBus .eq. _kidx-old        
      #if log && (b.fromBus !=kIdxBusIdxMap[b.fromBus]);  @info "fromBus-old: ", b.fromBus, " to: ", kIdxBusIdxMap[b.fromBus]; end        
      b.fromBus = kIdxBusIdxMap[b.fromBus]

    end

    if haskey(kIdxBusIdxMap, b.toBus)
      #if log && (b.toBus !=kIdxBusIdxMap[b.toBus]); @info "toBus: old toNode: ", b.toBus, " new toNode: ", kIdxBusIdxMap[b.toBus]; end
      b.toBus = kIdxBusIdxMap[b.toBus]
    end
  end
end







