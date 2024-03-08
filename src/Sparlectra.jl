module Sparlectra
# Author: Udo Schmitz (https://github.com/Welthulk)
# Purpose: network calculation

## R E S D A T A T Y P E S ##############################################################################################################################################################################################################
module ResDataTypes

# resource data types for working with Sparlectra
const Wurzel3 = 1.7320508075688772
const SparlectraVersion = VersionNumber("0.5.003")
export
  # constants
  Wurzel3,
  # classes
  ComponentTyp,
  Component,
  ImpPGMComp,
  ImpPGMComp3WT,
  Terminal,
  Node,
  NodeParameters,
  ACLineSegment,
  TrafoTyp,
  PowerTransformerTaps,
  PowerTransformerWinding,
  PowerTransformer,
  TransformerModelParameters,  
  ProSumer,
  BranchModel,
  BranchFlow,
  Branch,
  Shunt,
  Net,
  # functions  
  getRXBG,
  # Transformers
  getSideNumber2WT,
  getWinding2WT,
  calcTransformerRatio,
  create3WTWindings!,
  getTrafoImpPGMComp,  
  getWT3AuxBusID,
  isPerUnit_RXGB,
  # nodes
  addGenAktivePower!,
  addGenReaktivePower!,
  toComponentTyp,
  toNodeType,
  nodeComparison,
  busComparison,
  setGenPower!,
  
  
  
  setRatedS!,
  setNodePQ!,
  getShuntPGMComp,
  setShuntPower!,
  setVmVa!,
  setNodeParameters!,
  addAktivePower!,
  addReaktivePower!,
  getProSumPGMComp,
  hasPowerInjection,
  hasShuntInjection,
  # shunts
  setBranchFlow!,
  setBranchStatus!,
  setAdjElecParam!,
  # Shunt
  getGBShunt,
  getPQShunt,
  # ACLineSegment
  get_line_parameters,
  # ProSumer
  isSlack,
  isGenerator,
  toProSumptionType

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



abstract type AbstractBranch end

include("lines.jl")
include("transformer.jl")
include("prosumer.jl")
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

module SparlectraImport
using Sparlectra
using Sparlectra.ResDataTypes
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
using LinearAlgebra
using SparseArrays
using Printf
using Logging

export
  # constants
  # classes
  # functions
  createNetFromPGM,
  createNetFromMatPowerFile,
  casefileparser,
  calcNeutralU,
  recalc_trafo_model_data,
  createYBUS,

  getNBI,
  mdoRCM,
  createBranchVectorFromNodeVector!,
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
include("createnet_pgm.jl")
include("createnet_powermat.jl")


end


module SparlectraResult
using Sparlectra
using Sparlectra.ResDataTypes
using Printf

export printACPFlowResults

include("results.jl")

end
end # module Sparlectra
