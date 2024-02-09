module Sparlectra
# Author: Udo Schmitz (https://github.com/Welthulk)
# Purpose: network calculation

## R E S D A T A T Y P E S ##############################################################################################################################################################################################################
module ResDataTypes

# resource data types for working with Sparlectra
const Wurzel3 = 1.7320508075688772
const SparlectraVersion = VersionNumber("0.4.100")
export
  # constants
  Wurzel3,
  # classes
  ComponentTyp,
  Component,
  ImpPGMComp,
  Terminal,
  Node,
  NodeParameters,
  ACLineSegment,
  TrafoTyp,
  PowerTransformerTaps,
  PowerTransformerWinding,
  PowerTransformer,
  TransformesAdditionalParameters,
  ProSumer,
  BranchFlow,
  Branch,
  Shunt,
  Net,
  # functions
  addGenAktivePower!,
  addGenReaktivePower!,
  toComponentTyp,
  toNodeType,
  nodeComparison,
  busComparison,
  getNodeName,
  getNodeID,
  getNodeVn,
  getNodeType,
  setGenPower!,
  setNodeIdx!,
  setBusIdx!,
  setNodeType!,
  setRatedS!,
  setNodePQ!,
  setShuntPower!,
  setVmVa!,
  setNodeParameters!,
  addAktivePower!,
  addReaktivePower!,
  hasPowerInjection,
  hasShuntInjection,
  setBranchFlow!,
  setBranchStatus!,
  get_line_parameters

include("component.jl")

# helper
function toString(o::NodeType)::String
  if o == PQ
    return "PQ"
  elseif o == PV
    return "PV"
  elseif o == Slack
    return "Slack"
  elseif o == Isolated
    return "Isolated"
  else
    return "UnknownN"
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

function isPUNode(type::ComponentTyp)::Bool
  if type == SynchronousMachine
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

function toNodeType(o::ComponentTyp, r::Integer)::NodeType
  if r == 1
    return Slack
  elseif o == Generator || o == ExternalNetworkInjection || o == SynchronousMachine
    return PV
  elseif o == AsynchronousMachine || o == Load || o == EnergyConsumer
    return PQ
  else
    return UnknownP
  end
end

# helper
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
function toProSumptionType(o::ComponentTyp)::ProSumptionType
  if o == Generator || o == ExternalNetworkInjection || o == SynchronousMachine
    return Injection
  elseif o == AsynchronousMachine || o == Load || o == EnergyConsumer
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

include("lines.jl")
include("transformer.jl")

# helper
function toString(o::AbstractComponent)::String
  return "Name: " * o.cName * ", Typ: " * string(o.cTyp)
end

# helper
function toSeitenTyp(o::String)::SeitenTyp
  if occursin("1", o)
    return Seite1
  elseif occursin("2", o)
    return Seite2
  elseif occursin("3", o)
    return Seite3
  else
    return UnknownS
  end
end

include("terminal.jl")
include("prosumer.jl")

struct invalidException <: Exception
  message::String
end

include("node.jl")
include("branch.jl")
include("shunt.jl")

struct Net
  name::String
  baseMVA::Float64
  slack::Integer
  nodeVec::Vector{ResDataTypes.Node}
  linesAC::Vector{ResDataTypes.ACLineSegment}
  trafos::Vector{ResDataTypes.PowerTransformer}
  branchVec::Vector{ResDataTypes.Branch}
  prosumpsVec::Vector{ResDataTypes.ProSumer}
  shuntVec::Vector{ResDataTypes.Shunt}

  function Net(
    name,
    baseMVA::Float64,
    slack::Integer,
    nodeVec::Vector{ResDataTypes.Node},
    linesAC::Vector{ResDataTypes.ACLineSegment},
    trafos::Vector{ResDataTypes.PowerTransformer},
    branchVec::Vector{ResDataTypes.Branch},
    prosumpsVec::Vector{ResDataTypes.ProSumer},
    shuntVec::Vector{ResDataTypes.Shunt},
  )
    new(name, baseMVA, slack, nodeVec, linesAC, trafos, branchVec, prosumpsVec, shuntVec)
  end
end

end # module ResDataTypes
## S P A R Q L Q U E R Y C G M E S ###############################################################################################################################################################################################
module SparqlQueryCGMES

using HTTP
using JSON
using Sparlectra.ResDataTypes

export
  # constants
  # classes
  # functions
  getNodes,
  getLines,
  getTrafos,
  getProSumption

const NodeVector = Vector{ResDataTypes.Node}
const ACLineVector = Vector{ResDataTypes.ACLineSegment}
const TrafoVector = Vector{ResDataTypes.PowerTransformer}
const ProSumptionVector = Vector{ResDataTypes.ProSumer}
const ShuntVector = Vector{ResDataTypes.Shunt}

# helper
function isStringEmpty(str::AbstractString)
  if str == "" || length(str) == 0
    return true
  else
    return false
  end
end

# helper
function getFloatValue(b::AbstractDict, k::String, v::String = "value")
  try
    str = b[k][v]
    result = parse(Float64, str)
  catch
    result = nothing
  end
end
# helper
function getIntegerValue(b::AbstractDict, k::String, v::String = "value")
  try
    str = b[k][v]
    result = parse(Int64, str)
  catch
    result = nothing
  end
end
# helper
function getBoolValue(b::AbstractDict, k::String, v::String = "value")
  try
    str = b[k][v]
    result = parse(Bool, str)
  catch
    result = nothing
  end
end
# helper
function getStringValue(b::AbstractDict, k::String, v::String = "value")
  try
    str = b[k][v]
    result = str
  catch
    result = ""
  end
end

include("cgmesprosumptionquery.jl")
include("cgmestrafoquery.jl")
include("cgmesnodesquery.jl")
include("cgmeslinesquery.jl")

header2 = Dict("Accept" => "application/sparql-results+json", "Content-Type" => "application/sparql-query")

function getNodes(url::String, debug::Bool)::NodeVector
  myurl = url #endpoint
  mynodes = NodeVector()

  success = SparqlQueryCGMES.QueryNodes!(myurl, mynodes, debug)
  if !success
    @warn "No result found"
  end
  return mynodes
end

function getLines(url::String, debug::Bool)
  myurl = url #endpoint
  ACLines = ACLineVector()
  success = QueryLines!(myurl, ACLines, debug)
  if success
    if debug
      for acseg in ACLines
        @show acseg
      end
    end
    return ACLines
  else
    #println("No result found")
    return nothing
  end
end

function getTrafos(url::String, debug::Bool)
  myurl = url #endpoint
  trafosVec = TrafoVector()
  success = QueryTrafos!(myurl, trafosVec, debug)
  if success
    if debug
      for trafo in trafosVec
        @show trafo
      end
    end
    return trafosVec
  else
    return nothing
  end
end

function getProSumption(url::String, debug::Bool)
  myurl = url #endpoint
  prosumptionVec = ProSumptionVector()
  success = QueryProSumption!(myurl, prosumptionVec, debug)
  if success
    if debug
      for prosumption in prosumptionVec
        @show prosumption
      end
    end
    return prosumptionVec
  else
    #println("kein Ergebnis")
    return nothing
  end
end

end # module SparqlQueryCGMES
module SparlectraTools
using Sparlectra.ResDataTypes
export
  # constants
  # classes
  # functions
  getNodeByName,
  getTerminalsByName,
  reportCrusialData,
  checkNodeIdsUnique,
  checkNodeNamesUnique,
  checkTerminalHasTwoSides,
  checkTerminalOnlyOnce,
  roundUpToNearest100,
  isGenerator,
  isExternalNetworkInjection,
  isShunt,
  isLoad,
  isMotor,
  isSlack,
  searchNodeByTerminalId,
  findNetTrails,
  findSingleConnectedBuses,
  checkLineIdIsUnique,
  searchNodeByName,
  searchTerminalByNameAndSide,
  checkTransformerIdIsUnique,
  checkProSumerIdIsUnique,
  checkNodeConnections,
  renumberNodeIdx!

include("tools.jl")

end # module SparlectraTools
module SparlectraImport
using Sparlectra
using Sparlectra.ResDataTypes
using Sparlectra.SparlectraTools
using RegularExpressions
using UUIDs
using JSON

export 
  # constants
  # classes
  # functions
  jsonparser, 
  casefileparser,
  pgmparser

include("import.jl")


end # module SparlectraImport


module SparlectraExport
using Sparlectra
using Sparlectra.ResDataTypes
using Sparlectra.SparlectraTools
using JSON
using DataStructures

export
# constants
# classes
# functions
  writeMatpowerCasefile,
  exportPGM

include("exportMatPower.jl")
include("exportPGM.jl")
include("equicircuit.jl")
end # module SparlectraExport


module SparlectraNet
using Sparlectra
using Sparlectra.ResDataTypes
using Sparlectra.SparlectraImport
using Sparlectra.SparqlQueryCGMES
using Sparlectra.SparlectraTools
using UUIDs
using LinearAlgebra
using SparseArrays
using Printf
using Logging

export
  # constants
  # classes
  # functions
  createNetFromFile,
  createNetFromMatPowerFile,
  createNetFromTripleStore,
  createNetFromPGM,
  casefileparser,
  calc_y_pu,
  calcPQ_Shunt,
  calcGB_Shunt,
  calcTwoPortPU,
  calcVKDependence,
  calcComplexRatio,
  calcTapStepPercent,
  calcTapCorr,
  calcRatio,
  calcNeutralU,
  calcTrafoParamsSI,
  calcTrafoParams,
  calc3WTParams,
  calcYShunt,
  createYBUS,
  adjacentBranches,
  getNBI,
  mdoRCM,
  createBranchVectorFromNodeVector,
  fixSequenceNumberInNodeVec!,
  setParallelBranches!,
  calcJacobian,
  getBusData,
  calcPowerFlow,
  getPowerFeeds,
  adjacentBranches,
  calcNewtonRaphson!,
  calcNetLosses!

include("equicircuit.jl")
include("jacobian.jl")
include("losses.jl")
include("nbi.jl")
include("createnet.jl")
include("readnetfromfile.jl")
include("readpowermat.jl")


end


module SparlectraResult
using Sparlectra
using Sparlectra.ResDataTypes
using Printf

export printACPFlowResults

include("results.jl")

end
end # module Sparlectra
