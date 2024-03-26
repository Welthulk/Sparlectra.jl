module Sparlectra
# Author: Udo Schmitz (https://github.com/Welthulk)
# Purpose: network calculation

## resource data types for working with Sparlectra
module ResDataTypes

# resource data types for working with Sparlectra
const Wurzel3 = 1.7320508075688772
const SparlectraVersion = VersionNumber("0.4.10")
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
  # Compomnent
  toComponentTyp,
  getRXBG,
  # Transformers
  getSideNumber2WT,
  getWinding2WT,
  calcTransformerRatio,
  recalc_trafo_model_data,
  create3WTWindings!,
  getTrafoImpPGMComp,
  getWT3AuxBusID,
  isPerUnit_RXGB,
  # Nodes  
  setRatedS!,
  setVmVa!,
  setNodeParameters!,
  addShuntPower!,
  addLoadPower!,
  addGenPower!,    
  toNodeType,
  busComparison,  
  toString,
  # Branch
  setBranchFlow!,
  setBranchStatus!,
  # Shunt
  getGBShunt,
  getPQShunt,
  # ACLineSegment
  get_line_parameters,
  # ProSumer
  isSlack,
  isGenerator,
  isAPUNode,
  setQGenReplacement!,
  getQGenReplacement,
  toProSumptionType

include("component.jl")

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

module SparlectraNet
using Sparlectra
using Sparlectra.ResDataTypes
using LinearAlgebra
using SparseArrays
using Printf
using Logging

export
  # constants
  # classes
  # functions  
  casefileparser,
  createNetFromMatPowerFile,
  calcComplexRatio,
  calcNeutralU,  
  createYBUS,
  getNBI,
  mdoRCM,
  setJacobianDebug,
  calcJacobian,
  calcPowerFlow,
  adjacentBranches,
  calcNewtonRaphson!,
  runpf!,
  calcNetLosses!,
  writeMatpowerCasefile,
  printACPFlowResults, 
  convertPVtoPQ!

include("import.jl")
include("equicircuit.jl")
include("jacobian.jl")
include("losses.jl")
include("nbi.jl")
include("createnet_powermat.jl")
include("exportMatPower.jl")
include("results.jl")
end

end # module Sparlectra
