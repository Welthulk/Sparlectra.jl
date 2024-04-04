module Sparlectra
# Author: Udo Schmitz (https://github.com/Welthulk)
# Purpose: network calculation

# Naming Conventions:
# The project follows the Julia Naming Conventions for the most part, 
# but it's important to note that the naming convention for functions might deviate. 
# In this module, functions are written in CamelCase with a lowercase initial letter. 
using LinearAlgebra
using SparseArrays
using Printf
using Logging


# resource data types for working with Sparlectra
const Wurzel3 = 1.7320508075688772
const SparlectraVersion = VersionNumber("0.4.10")
abstract type AbstractBranch end

export
  # constants
  Wurzel3, ComponentTyp,
  # classes  
  # Component
  AbstractComponent, Component, ImpPGMComp, ImpPGMComp3WT, 
  # Node
  Node, 
  # Line
  ACLineSegment,
  # Trafo
  TrafoTyp, PowerTransformerTaps,  PowerTransformerWinding,  PowerTransformer, TransformerModelParameters,
  # ProSumer
  ProSumer,
  # Branch
  AbstractBranch, Branch,BranchModel,  BranchFlow, getBranchFlow, setBranchFlow!, setBranchStatus!,getBranchLosses, setBranchLosses!,
  # Shunt
  Shunt,
  # Net
  Net,
  
  # functions  
  # Compomnent
  toComponentTyp,  getRXBG,
  # Transformers
  getSideNumber2WT,  getWinding2WT,  calcTransformerRatio, recalc_trafo_model_data, create2WTRatioTransformerNoTaps, create3WTWindings!,
  getTrafoImpPGMComp,  getWT3AuxBusID,  isPerUnit_RXGB,
  # Nodes  
  setRatedS!,  setVmVa!,  addShuntPower!,  addLoadPower!,  addGenPower!,  getNodeVn,  isSlack,  isPVNode,  isPQNode, toNodeType,
  busComparison,   toString,
  # Branch
  setBranchFlow!,  setBranchStatus!,
  # Shunt
  getGBShunt,  getPQShunt,
  # ACLineSegment
  get_line_parameters,
  # ProSumer
  isSlack, isGenerator, isAPUNode, setQGenReplacement!, getQGenReplacement, toProSumptionType,
  # Net
  addBus!, addShunt!, addACLine!, add2WTrafo!, addPIModellTrafo!, addProsumer!, lockNet!, validate, hasBusInNet, getNetOrigBusIdx, geNetBusIdx, setTotalLosses!, getTotalLosses, getBusType, get_bus_vn_kV, get_vn_kV,
  # create_powermat.jl
  casefileparser, createNetFromMatPowerFile,
  # exportMatPower.jl
  writeMatpowerCasefile,
  # equicircuit.jl
  calcComplexRatio, calcNeutralU,  createYBUS, adjacentBranches,
  # nbi.jl
  getNBI, mdoRCM,
  # jacobian.jl
  setJacobianDebug, runpf!,
  # losses.jl
  calcNetLosses!,
  
  # results.jl
  printACPFlowResults, convertPVtoPQ!,

  # run_acpflow.jl
  run_acpflow


include("component.jl")
include("lines.jl")
include("transformer.jl")
include("prosumer.jl")
include("node.jl")
include("branch.jl")
include("shunt.jl")
include("network.jl")
include("import.jl")
include("equicircuit.jl")
include("jacobian.jl")
include("losses.jl")
include("nbi.jl")
include("createnet_powermat.jl")
include("exportMatPower.jl")
include("results.jl")
include("run_acpflow.jl")


end # module Sparlectra
