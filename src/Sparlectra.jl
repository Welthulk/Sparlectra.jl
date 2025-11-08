# Author: Udo Schmitz (https://github.com/Welthulk)
# Purpose: network calculation

# Naming Conventions:
# The project follows the Julia Naming Conventions for the most part, 
# but it's important to note that the naming convention for functions might deviate. 
# In this module, functions are written in CamelCase with a lowercase initial letter. 

#! format: off

module Sparlectra
using LinearAlgebra, SparseArrays, Printf, Logging

# resource data types for working with Sparlectra
const Wurzel3 = 1.7320508075688772
const SparlectraVersion = v"0.4.23"
version() = SparlectraVersion
abstract type AbstractBranch end

const _MODULE_DOC = """
    Sparlectra $(SparlectraVersion)

Sparlectra is a Julia package for the calculation of electrical networks. 
It is designed to be used in the context of power system analysis and optimization. 

- GitHub Repository: https://github.com/welthulk/Sparlectra.jl
- Website: https://welthulk.github.io/Sparlectra.jl

**Naming Conventions:**
The project follows the Julia Naming Conventions for the most part, 
but it's important to note that the naming convention for functions might deviate. 
In this module, functions are written in CamelCase with a lowercase initial letter.
"""
@doc _MODULE_DOC Sparlectra
export
  version,
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
  AbstractBranch, Branch, BranchModel, BranchFlow, getBranchFlow, setBranchFlow!, getBranchNumber, getBranchLosses, setBranchLosses!, setBranchStatus!,
  # Shunt
  Shunt,
  # Net
  Net,
  
  # functions  
  # Compomnent
  toComponentTyp, getCompName, getCompID, 
  # Transformers
  getSideNumber2WT,  getWinding2WT,  calcTransformerRatio, recalc_trafo_model_data, create2WTRatioTransformerNoTaps, create3WTWindings!,
  getTrafoImpPGMComp,  getWT3AuxBusID,  isPerUnit_RXGB, getWindingRatedS, getTrafoRXBG, getTrafoRXBG_pu, 
  # Nodes  
  setRatedS!,  setVmVa!,  addShuntPower!,  addLoadPower!,  addGenPower!,  getNodeVn,  isSlack,  isPVNode,  isPQNode, isIsolated, toNodeType, setNodeType!, setBusType!,
  busComparison, toString,
  # Branch
  getBranchIdx, calcBranchYser, calcBranchYshunt, calcBranchRatio, calcAdmittance,
  # Shunt
  getGBShunt,  getPQShunt, updatePQShunt!,
  # ACLineSegment
  get_line_parameters, isLinePIModel, getLineRXBG, getLineRXBG_pu,
  # ProSumer
  isSlack, isGenerator, isAPUNode, setQGenReplacement!, getQGenReplacement, toProSumptionType, updatePQ!,
  # Net
  addBus!, addShunt!, addACLine!, addPIModelACLine!, add2WTrafo!, addPIModelTrafo!, addProsumer!, lockNet!, validate!, hasBusInNet, addBusGenPower!, addBusLoadPower!, addBusShuntPower!, setNodeVoltage!, setNodeAngle!,
  getNetOrigBusIdx, geNetBusIdx, setNetBranchStatus!, getNetBranch, getNetBranchNumberVec, setTotalLosses!, getTotalLosses, getBusType, get_bus_vn_kV, get_vn_kV, updateBranchParameters!, hasShunt!, 
  getShunt!, markIsolatedBuses!,setTotalBusPower!, setPVGeneratorQLimitsAll!, setPVBusVset!,setPVBusQLimits!,setAllPVVset!,rebuildQLimits!,
  # remove_functions.jl
  removeBus!, removeBranch!, removeACLine!, removeTrafo!, removeShunt!, removeProsumer!, clearIsolatedBuses!,
  # import.jl
  casefileparser, createNetFromMatPowerFile,
  # exportMatPower.jl
  writeMatpowerCasefile,
  # equicircuit.jl
  calcComplexRatio, calcNeutralU,  createYBUS, adjacentBranches, toPU_RXBG, to_RXBG, 
  # nbi.jl
  getNBI, mdoRCM,
  # jacobian.jl
  setJacobianDebug, setJacobianAngleLimit, runpf!,
  # jacobian_full.jl (neu)
  getPowerFeeds_full, residuum_full_withPV, calcJacobian_withPVIdentity, calcNewtonRaphson_withPVIdentity!, runpf_full!,
  # limits.jl
  printQLimitLog,logQLimitHit!, lastQLimitIter, get_Q_limits_pu, buildQLimits!,logQLimitHit!,lastQLimitIter, resetQLimitLog!,
  # losses.jl
  calcNetLosses!,
  
  # results.jl
  printACPFlowResults, convertPVtoPQ!,

  # run_acpflow.jl
  run_acpflow, run_net_acpflow

include("component.jl")
include("lines.jl")
include("transformer.jl")
include("prosumer.jl")
include("node.jl")
include("branch.jl")
include("shunt.jl")
include("busdata.jl")
include("network.jl")
include("import.jl")
include("equicircuit.jl")
include("jacobian.jl")
include("limits.jl")
include("jacobian_full.jl")
include("losses.jl")
include("nbi.jl")
include("createnet_powermat.jl")
include("exportMatPower.jl")
include("results.jl")
include("run_acpflow.jl")
include("remove_functions.jl") 

#! format: on
end # module Sparlectra
