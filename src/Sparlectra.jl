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

using BenchmarkTools
using Dates
using LinearAlgebra
using Logging
using Printf
using SparseArrays

const MPOWER_DIR = normpath(joinpath(pkgdir(@__MODULE__), "data", "mpower"))


# resource data types for working with Sparlectra
const Wurzel3 = 1.7320508075688772
const SparlectraVersion = v"0.8.2"
version() = SparlectraVersion
abstract type AbstractBranch end

const _MODULE_DOC = """
    Sparlectra $(SparlectraVersion)

Sparlectra is a Julia package for the calculation of electrical networks.
It is designed to be used in the context of power system analysis and optimization.

- GitHub Repository: https://github.com/welthulk/Sparlectra.jl
- Website: https://welthulk.github.io/Sparlectra.jl

"""
@doc _MODULE_DOC Sparlectra
export
  version,                                # Return the package version.

  # constants
  Wurzel3,
  ComponentTyp,

  # classes

  # Component
  AbstractComponent,
  Component,
  ImpPGMComp,
  ImpPGMComp3WT,

  # Node
  Node,

  # Line
  ACLineSegment,

  # Trafo
  TrafoTyp,
  PowerTransformerTaps,
  PowerTransformerWinding,
  PowerTransformer,
  TransformerModelParameters,

  # ProSumer
  ProSumer,
  AbstractVoltageDependentController,
  PiecewiseLinearCharacteristic,
  QUController,
  PUController,
  VoltageAdjustConfig,

  # Tap controller
  AbstractOuterController,
  AbstractControlState,
  AbstractControlUpdate,
  PowerTransformerControl,

  # Branch
  AbstractBranch,
  Branch,
  BranchModel,
  BranchFlow,
  getBranchFlow,
  setBranchFlow!,
  getBranchNumber,
  getBranchLosses,
  setBranchLosses!,
  setBranchStatus!,

  # Link
  BusLink,
  setLinkStatus!,
  setLinkFlow!,
  setLinkCurrent!,

  # Shunt
  Shunt,

  # Net
  Net,

  # functions

  # utilities.jl
  zero_row!,
  print_jacobian,

  # yamlparams.jl
  parse_yaml_scalar,
  load_yaml_dict,
  merge_yaml_dict!,
  as_bool,
  as_int_vector,

  # configuration.jl
  StartModeConfig,
  QLimitConfig,
  PowerFlowConfig,
  ControlConfig,
  ObservabilityConfig,
  StateEstimationConfig,
  MatpowerImportConfig,
  PerformanceConfig,
  BenchmarkConfig,
  RuntimeConfig,
  DiagnosticsConfig,
  OutputConfig,
  SparlectraConfig,
  load_sparlectra_config,                 # Load a configuration file without activating it.
  load_sparlectra_config!,                # Load and activate a configuration file.
  set_sparlectra_config!,                 # Replace the active global configuration.
  active_sparlectra_config,               # Return the currently active configuration.
  powerflow_config,                       # Access the active power-flow settings.
  matpower_import_config,                 # Access MATPOWER import settings.
  state_estimation_config,
  diagnostics_config,
  output_config,
  performance_config,
  benchmark_config,
  runtime_config,
  print_effective_config,                 # Print the resolved effective configuration.
  configuration_path_from_inputs,         # Resolve the configuration file path.
  run_matpower_case,                      # Run a MATPOWER case through the high-level workflow.
  run_synthetic_tiled_grid_pf_perf,       # Run synthetic-grid PF performance tests.
  run_voltage_dependent_control_demo,

  # synthetic_grids.jl
  synthetic_tiled_grid_bus_index,
  build_synthetic_tiled_grid_net,         # Build a synthetic benchmark network.
  build_tiled_grid_net,

  # BusData
  BusData,
  getBusData,
  getBusTypeVec,
  countNodes,
  map_NR_voltage_to_net!,                 # Copy NR voltage results back into the network.
  buildVoltageVector_from_busVec,         # Build a voltage vector from bus data.

  # Compomnent
  toComponentTyp,
  getCompName,
  getCompID,

  # Transformers
  getSideNumber2WT,
  getWinding2WT,
  calcTransformerRatio,                   # Compute transformer off-nominal ratio.
  recalc_trafo_model_data,                # Recompute transformer equivalent model data.
  create2WTRatioTransformerNoTaps,
  create3WTWindings!,
  getTrafoImpPGMComp,
  getWT3AuxBusID,
  isPerUnit_RXGB,
  getWindingRatedS,
  getTrafoRXBG,
  getTrafoRXBG_pu,

  # Nodes
  setRatedS!,
  setVmVa!,                               # Set voltage magnitude and angle on a node.
  addShuntPower!,
  addLoadPower!,
  addGenPower!,
  getNodeVn,
  isSlack,
  isPVNode,
  isPQNode,
  isIsolated,
  toNodeType,
  setNodeType!,
  getNodeType,
  busComparison,
  toString,

  # Branch
  getBranchIdx,
  calcBranchYser,
  calcBranchYshunt,
  calcBranchRatio,
  calcAdmittance,

  # Link
  BusLink,
  setLinkStatus!,
  setLinkFlow!,
  setLinkCurrent!,

  # Shunt
  getGBShunt,
  getPQShunt,
  updatePQShunt!,

  # ACLineSegment
  get_line_parameters,
  isLinePIModel,
  getLineRXBG,
  getLineRXBG_pu,

  # ProSumer
  isSlack,
  isGenerator,
  isAPUNode,
  isRegulating,
  setQGenReplacement!,
  getQGenReplacement,
  toProSumptionType,
  updatePQ!,
  getPosumerBusIndex,
  setPQResult!,
  evaluate_characteristic,                # Evaluate a voltage-dependent characteristic.
  evaluate_controller,                    # Evaluate a voltage-dependent controller.
  make_characteristic,                    # Create a controller characteristic.
  has_qu_controller,
  has_pu_controller,

  # Network
  addBus!,                                # Add a bus to the network.
  addShunt!,                              # Add a shunt element.
  addACLine!,                             # Add an AC line.
  addPIModelACLine!,                      # Add an AC line as a PI model.
  add2WTrafo!,                            # Add a two-winding transformer.
  addPIModelTrafo!,                       # Add a transformer PI model.
  addProsumer!,                           # Add generator/load prosumer data.
  lockNet!,                               # Finalize network topology before solving.
  validate!,                              # Validate network consistency.
  hasBusInNet,
  addBusGenPower!,
  addBusLoadPower!,
  addBusShuntPower!,
  setNodeVoltage!,
  setNodeAngle!,
  getNetOrigBusIdx,
  geNetBusIdx,
  setNetBranchStatus!,
  getNetBranch,
  getNetBranchNumberVec,
  setTotalLosses!,
  getTotalLosses,
  getBusType,
  getEffectiveBusType,
  getBusProsumers,
  refreshBusTypesFromProsumers!,          # Rebuild effective bus types from prosumers.
  get_bus_vn_kV,
  get_vn_kV,
  updateBranchParameters!,
  hasShunt!,
  getShunt!,
  markIsolatedBuses!,
  setTotalBusPower!,
  setPVBusVset!,                          # Set PV voltage target.
  setQLimits!,                            # Set generator reactive-power limits.
  getNodeVm,
  distributeBusResults!,
  getTotalBusPower,
  getTotalLosses,
  buildVoltageVector,                     # Build complex bus-voltage vector.
  initialVrect,                           # Build the rectangular initial voltage vector.
  buildComplexSVec,                       # Build specified complex power vector.
  buildControlledSVec,                    # Build controlled specified power vector.
  has_voltage_dependent_control,
  addShuntMatpower!,
  normalize_bus_shunt_model,              # Normalize bus-shunt modeling mode.
  bus_shunt_totals_pu,
  log_bus_shunt_model,
  add2WTPIModelTrafo!,
  add3WTPiModelTrafo!,
  showNet,                                # Print a network summary.
  buildQLimits!,                          # Build Q-limit arrays for the solver.
  updateShuntPowers!,
  addLink!,
  setNetLinkStatus!,
  getNetLinks,
  calcLinkFlowsKCL!,                      # Compute bus-link flows from KCL.
  collect_outer_controllers,
  run_control!,                           # Run outer-loop control workflow.
  latest_control_result,
  ControlRunResult,
  addPowerTransformerControl!,            # Add transformer control data.
  addTapController!,                      # Add a tap-controller wrapper.
  clearTapControllers!,
  get_bus_vm_pu,
  get_branch_p_from_to_mw,
  get_branch_q_from_to_mvar,
  buildTapControllerReportRows,
  printTapControllerSummary,

  # remove_functions.jl
  removeBus!,                             # Remove a bus and dependent data.
  removeBranch!,                          # Remove a branch.
  removeACLine!,                          # Remove an AC line.
  removeTrafo!,                           # Remove a transformer.
  removeShunt!,                           # Remove a shunt.
  removeProsumer!,                        # Remove prosumer data.
  clearIsolatedBuses!,                    # Remove isolated imported buses.
  apply_mp_isolated_buses!,

  # import.jl
  createNetFromMatPowerFile,              # Import a MATPOWER case file as Net.
  _createDict,
  apply_matpower_bus_voltage!,            # Apply MATPOWER bus voltage data.
  apply_mp_bus_vmva_init!,                # Initialize Vm/Va from MATPOWER data.

  # exportMatPower.jl
  writeMatpowerCasefile,                  # Export a network as MATPOWER case file.

  # equicircuit.jl
  calcComplexRatio,
  calcNeutralU,
  createYBUS,                             # Build the network admittance matrix.
  adjacentBranches,
  toPU_RXBG,
  fromPU_RXBG,
  branchFlow_pu,                          # Compute branch flow in per unit.

  # nbi.jl
  getNBI,                                 # Compute node-branch incidence data.
  mdoRCM,                                 # Run RCM ordering helper.

  # jacobian.jl
  runpf!,                                 # Run the classic Newton-Raphson PF solver.
  setJacobianDebug,
  setJacobianAngleLimit,

  # jacobian_full.jl
  runpf_full!,                            # Run full/polar PF variant.

  # limits.jl
  printQLimitLog,
  printPVQLimitsTable,                    # Print PV Q-limit diagnostics.
  printFinalLimitValidation,              # Print final Q-limit validation.
  validate_q_limit_signs!,                # Validate Q-limit sign conventions.
  logQLimitHit!,
  lastQLimitIter,
  getQLimits_pu,
  logQLimitHit!,
  lastQLimitIter,
  resetQLimitLog!,                        # Clear Q-limit event log.
  pv_hit_q_limit,
  has_q_limits,
  active_set_q_limits!,                   # Apply active-set PV/PQ Q-limit logic.

  # losses.jl
  calcNetLosses!,                         # Compute total network losses.

  # results.jl
  ACPFlowReport,
  buildACPFlowReport,                     # Build structured AC PF report.
  printACPFlowResults,                    # Print AC PF result report.
  printProsumerResults,

  # run_acpflow.jl
  run_acpflow,                            # Public high-level AC PF runner.
  run_matpower_case,                      # Run a MATPOWER case through the high-level workflow.

  # solver_core.jl
  calc_injections,                        # Compute complex bus injections.
  calc_currents,                          # Compute bus currents.
  solve_linear,
  solve_sparse_system,                    # Solve a sparse linear system.
  build_pos_map,
  slack_elimination_indices,
  extract_bus_types_and_vset,             # Extract solver bus types and V targets.
  build_qload_pu,
  build_voltage_vector,                   # Build solver voltage vector.
  compute_sbus_and_totals,                # Build Sbus and aggregate P/Q totals.

  # measurements.jl / state_estimation.jl
  MeasurementType,
  Measurement,
  measurementStdDevs,
  generateMeasurementsFromPF,             # Generate synthetic SE measurements from PF.
  setMeasurementsFromPF!,                 # Replace measurements from PF results.
  addMeasurement!,
  addVmMeasurement!,                      # Add voltage-magnitude measurement.
  addPinjMeasurement!,                    # Add active-power injection measurement.
  addQinjMeasurement!,                    # Add reactive-power injection measurement.
  addPflowMeasurement!,                   # Add active branch-flow measurement.
  addQflowMeasurement!,                   # Add reactive branch-flow measurement.
  findPassiveBuses,
  addZeroInjectionMeasurements!,          # Add zero-injection constraints.
  SEResult,
  runse!,                                 # Run weighted least-squares state estimation.
  numeric_rank,
  runse_diagnostics,                      # Run SE diagnostics.
  validate_measurements,
  summarize_se_diagnostics,
  print_se_diagnostics,
  numerical_observable,
  structural_observable,
  numerical_row_redundant,
  structural_row_redundant,
  evaluate_observability_matrix,
  evaluate_local_observability_matrix,
  evaluate_global_observability,          # Evaluate global SE observability.
  evaluate_local_observability,           # Evaluate local SE observability.

  # jacobian_complex.jl
  runpf_rectangular!,                     # Run rectangular complex-state NR PF.
  mismatch_rectangular,                   # Compute rectangular PF mismatch vector.
  project_rectangular_start,              # Build/select rectangular start candidate.

  # External solver interface
  PFModel,
  PFSolution,
  AbstractExternalSolver,
  buildPfModel,                           # Build external PF solver model.
  mismatchInf,                            # Compute infinity-norm mismatch.
  applyPfSolution!,                       # Apply an external PF solution to Net.
  solvePf,                                # Solve through an external solver interface.
  runpf_external!,                        # Run PF via external solver interface.
  ensure_casefile                         # Resolve or fetch a case file.


include("utilities.jl")
include("yamlparams.jl")
include("control_framework.jl")
include("configuration.jl")
include("component.jl")
include("lines.jl")
include("transformer.jl")
include("prosumer.jl")
include("node.jl")
include("branch.jl")
include("link.jl")
include("shunt.jl")
include("network.jl")
include("synthetic_grids.jl")
include("tap_control.jl")
include("busdata.jl")
include("MatpowerIO.jl")
include("createnet_powermat.jl")
include("equicircuit.jl")
include("limits.jl")
include("losses.jl")
include("exportMatPower.jl")
include("results.jl")
include("run_acpflow.jl")
include("matpower_runner.jl")
include("remove_functions.jl")
include("solver_core.jl")
include("powerflow_rectangular/rectangular_core_equations.jl")
include("powerflow_rectangular/rectangular_wrong_branch.jl")
include("jacobian_complex.jl")
include("solver_interface.jl")
include("FetchMatpowerCase.jl")
include("measurements.jl")
include("state_estimation.jl")
#! format: on
end # module Sparlectra
