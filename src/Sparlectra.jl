# Copyright 2023–2026 Udo Schmitz
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
# Purpose: network calculation

# Naming Conventions:
# The project follows the Julia Naming Conventions for the most part, 
# but it's important to note that the naming convention for functions might deviate. 
# In this module, functions are written in CamelCase with a lowercase initial letter. 

# file: src/Sparlectra.jl
#! format: off

module Sparlectra
using LinearAlgebra, Dates, SparseArrays, Printf, Logging

const MPOWER_DIR = normpath(joinpath(pkgdir(@__MODULE__), "data", "mpower"))


# resource data types for working with Sparlectra
const Wurzel3 = 1.7320508075688772
const SparlectraVersion = v"0.4.33"
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
  setRatedS!,  setVmVa!,  addShuntPower!,  addLoadPower!,  addGenPower!,  getNodeVn,  isSlack,  isPVNode,  isPQNode, isIsolated, toNodeType, setNodeType!,getNodeType, 
  busComparison, toString,
  # Branch
  getBranchIdx, calcBranchYser, calcBranchYshunt, calcBranchRatio, calcAdmittance,
  # Shunt
  getGBShunt,  getPQShunt, updatePQShunt!,
  # ACLineSegment
  get_line_parameters, isLinePIModel, getLineRXBG, getLineRXBG_pu,
  # ProSumer
  isSlack, isGenerator, isAPUNode, setQGenReplacement!, getQGenReplacement, toProSumptionType, updatePQ!, getPosumerBusIndex, setPQResult!,
  # Network
  addBus!, addShunt!, addACLine!, addPIModelACLine!, add2WTrafo!, addPIModelTrafo!, addProsumer!, lockNet!, validate!, hasBusInNet, addBusGenPower!, addBusLoadPower!, addBusShuntPower!, setNodeVoltage!, setNodeAngle!,
  getNetOrigBusIdx, geNetBusIdx, setNetBranchStatus!, getNetBranch, getNetBranchNumberVec, setTotalLosses!, getTotalLosses, getBusType, get_bus_vn_kV, get_vn_kV, updateBranchParameters!, hasShunt!, 
  getShunt!, markIsolatedBuses!,setTotalBusPower!, setPVBusVset!, setQLimits!, getNodeVm,distributeBusResults!, getTotalBusPower, getTotalLosses, buildVoltageVector,initialVrect, buildComplexSVec,
  add2WTPIModelTrafo!, add3WTPiModelTrafo!,showNet,
  # remove_functions.jl
  removeBus!, removeBranch!, removeACLine!, removeTrafo!, removeShunt!, removeProsumer!, clearIsolatedBuses!,
  # import.jl
  createNetFromMatPowerFile, _createDict,
  # exportMatPower.jl
  writeMatpowerCasefile,
  # equicircuit.jl
  calcComplexRatio, calcNeutralU,  createYBUS, adjacentBranches, toPU_RXBG, fromPU_RXBG, 
  # nbi.jl
  getNBI, mdoRCM,
  # jacobian.jl
  runpf!, setJacobianDebug, setJacobianAngleLimit,
  # jacobian_full.jl
  runpf_full!, 
  # limits.jl
  printQLimitLog,logQLimitHit!, lastQLimitIter, getQLimits_pu, logQLimitHit!,lastQLimitIter, resetQLimitLog!, pv_hit_q_limit,
  # losses.jl
  calcNetLosses!, 
  # results.jl
  printACPFlowResults, printProsumerResults,
  # run_acpflow.jl
  run_acpflow, run_net_acpflow,
  # jacobian_complex.jl
  runpf_rectangular!, mismatch_rectangular, complex_newton_step_rectangular_fd,
  # External solver interface
  PFModel, PFSolution, AbstractExternalSolver,
  buildPfModel, mismatchInf, applyPfSolution!, solvePf, runpf_external!,
  ensure_casefile

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
include("MatpowerIO.jl")   
include("createnet_powermat.jl")
include("equicircuit.jl")
include("jacobian.jl")
include("limits.jl")
include("jacobian_full.jl")
include("losses.jl")
include("nbi.jl")
include("exportMatPower.jl")
include("results.jl")
include("run_acpflow.jl")
include("remove_functions.jl") 
include("jacobian_complex.jl")
include("jacobian_fd.jl")
include("solver_interface.jl")
include("FetchMatpowerCase.jl")
#! format: on
end # module Sparlectra
