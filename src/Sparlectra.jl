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
const SparlectraVersion = v"0.4.25"
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
  # utilities.jl
  zero_row!, print_jacobian,
  # BusData
  BusData, getBusData, getBusTypeVec, countNodes, map_NR_voltage_to_net!,buildVoltageVector_from_busVec,
  # Compomnent
  toComponentTyp, getCompName, getCompID, 
  # Transformers
  getSideNumber2WT,  getWinding2WT,  calcTransformerRatio, recalc_trafo_model_data, create2WTRatioTransformerNoTaps, create3WTWindings!,
  getTrafoImpPGMComp,  getWT3AuxBusID,  isPerUnit_RXGB, getWindingRatedS, getTrafoRXBG, getTrafoRXBG_pu, 
  # Nodes  
  setRatedS!,  setVmVa!,  addShuntPower!,  addLoadPower!,  addGenPower!,  getNodeVn,  isSlack,  isPVNode,  isPQNode, isIsolated, toNodeType, setNodeType!, 
  busComparison, toString,
  # Branch
  getBranchIdx, calcBranchYser, calcBranchYshunt, calcBranchRatio, calcAdmittance,
  # Shunt
  getGBShunt,  getPQShunt, updatePQShunt!,
  # ACLineSegment
  get_line_parameters, isLinePIModel, getLineRXBG, getLineRXBG_pu,
  # ProSumer
  isSlack, isGenerator, isAPUNode, setQGenReplacement!, getQGenReplacement, toProSumptionType, updatePQ!, getPosumerBusIndex,
  # Network
  addBus!, addShunt!, addACLine!, addPIModelACLine!, add2WTrafo!, addPIModelTrafo!, addProsumer!, lockNet!, validate!, hasBusInNet, addBusGenPower!, addBusLoadPower!, addBusShuntPower!, setNodeVoltage!, setNodeAngle!,
  getNetOrigBusIdx, geNetBusIdx, setNetBranchStatus!, getNetBranch, getNetBranchNumberVec, setTotalLosses!, getTotalLosses, getBusType, get_bus_vn_kV, get_vn_kV, updateBranchParameters!, hasShunt!, 
  getShunt!, markIsolatedBuses!,setTotalBusPower!, setPVBusVset!, setQLimits!, 
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
  setJacobianDebug, setJacobianAngleLimit, runpf!, runpf_full!,
  # jacobian_full.jl
  getPowerFeeds_full, residuum_full_withPV, calcJacobian_withPVIdentity, calcNewtonRaphson_withPVIdentity!, runpf_full!, residuum_state_full_withPV,
  # limits.jl
  printQLimitLog,logQLimitHit!, lastQLimitIter, getQLimits_pu, logQLimitHit!,lastQLimitIter, resetQLimitLog!, pv_hit_q_limit,
  # losses.jl
  calcNetLosses!, buildVoltageVector,
  # results.jl
  printACPFlowResults, convertPVtoPQ!,
  # run_acpflow.jl
  run_acpflow, run_net_acpflow,
  # complex_state_nr.jl
  complex_newton_step, build_complex_state_vector, run_complex_nr

include("utilities.jl")
include("component.jl")
include("lines.jl")
include("transformer.jl")
include("prosumer.jl")
include("node.jl")
include("branch.jl")
include("shunt.jl")
include("network.jl")
include("busdata.jl")
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
include("complex_state_nr.jl")
include("jacobian_complex.jl")


#! format: on
end # module Sparlectra
