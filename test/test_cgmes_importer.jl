# Copyright 2023-2026 Udo Schmitz
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

# file: test/test_cgmes_importer.jl
# purpose: tests the CGMES importer end to end: RDF reader semantics and
#          profile classification on synthetic fixtures, mapping of taps,
#          controllers, machines and shunts, and SV comparison runs
# Tests for the CGMES importer (see docs/src/cgmes_import.md):
# generic reader semantics on a synthetic in-memory fixture, plus
# summarizeCGMES assertions on the ENTSO-E MicroGrid when the local test-set
# cache exists (fetched by examples/experimental/cgmes_fetch_testsets.jl;
# skipped otherwise — no network access from tests).

using Sparlectra.CGMESImporter:
  CGMESFile, CGMESStore, collectCGMESFiles, loadCGMES, summarizeCGMES, objectsOf, countOf, num, str, boolval, enumval, ref, unresolvedReferences, readCGMESFile!, CIMObject, importFailureAnalysis
import ZipArchives

const _CGMES_SYNTH_EQ = """<?xml version="1.0" encoding="UTF-8"?>
<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:cim="http://iec.ch/TC57/2013/CIM-schema-cim16#" xmlns:md="http://iec.ch/TC57/61970-552/ModelDescription/1#">
<md:FullModel rdf:about="urn:uuid:eq-model">
<md:Model.profile>http://entsoe.eu/CIM/EquipmentCore/3/1</md:Model.profile>
<md:Model.profile>http://entsoe.eu/CIM/EquipmentShortCircuit/3/1</md:Model.profile>
</md:FullModel>
<cim:BaseVoltage rdf:ID="_bv110">
<cim:BaseVoltage.nominalVoltage>110</cim:BaseVoltage.nominalVoltage>
</cim:BaseVoltage>
<cim:ACLineSegment rdf:ID="_line1">
<cim:IdentifiedObject.name>L1</cim:IdentifiedObject.name>
<cim:ACLineSegment.r>2.5</cim:ACLineSegment.r>
<cim:ACLineSegment.x>10.0</cim:ACLineSegment.x>
<cim:ACLineSegment.gch>0.001</cim:ACLineSegment.gch>
<cim:ACLineSegment.r0>7.5</cim:ACLineSegment.r0>
<cim:Conductor.length>12.0</cim:Conductor.length>
<cim:ConductingEquipment.BaseVoltage rdf:resource="#_bv110"/>
</cim:ACLineSegment>
<cim:SynchronousMachine rdf:ID="_sm1">
<cim:IdentifiedObject.name>G1</cim:IdentifiedObject.name>
<cim:SynchronousMachine.operatingMode rdf:resource="http://iec.ch/TC57/2013/CIM-schema-cim16#SynchronousMachineOperatingMode.generator"/>
<cim:RotatingMachine.ratedS>150</cim:RotatingMachine.ratedS>
<cim:SynchronousMachine.earthing>false</cim:SynchronousMachine.earthing>
<cim:Equipment.EquipmentContainer rdf:resource="#_missing_container"/>
</cim:SynchronousMachine>
<cim:TopologicalIsland rdf:ID="_island1">
<cim:IdentifiedObject.name>ISL1</cim:IdentifiedObject.name>
<cim:TopologicalIsland.TopologicalNodes rdf:resource="#_tn_a"/>
<cim:TopologicalIsland.TopologicalNodes rdf:resource="#_tn_b"/>
<cim:TopologicalIsland.TopologicalNodes rdf:resource="#_tn_c"/>
</cim:TopologicalIsland>
</rdf:RDF>
"""

const _CGMES_SYNTH_SSH = """<?xml version="1.0" encoding="UTF-8"?>
<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:cim="http://iec.ch/TC57/2013/CIM-schema-cim16#" xmlns:md="http://iec.ch/TC57/61970-552/ModelDescription/1#">
<md:FullModel rdf:about="urn:uuid:ssh-model">
<md:Model.profile>http://entsoe.eu/CIM/SteadyStateHypothesis/1/1</md:Model.profile>
</md:FullModel>
<cim:SynchronousMachine rdf:about="#_sm1">
<cim:RotatingMachine.p>-120</cim:RotatingMachine.p>
<cim:RotatingMachine.q>-30</cim:RotatingMachine.q>
</cim:SynchronousMachine>
</rdf:RDF>
"""

const _CGMES_SYNTH_DIFF = """<?xml version="1.0" encoding="UTF-8"?>
<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:cim="http://iec.ch/TC57/2013/CIM-schema-cim16#" xmlns:dm="http://iec.ch/TC57/61970-552/DifferenceModel/1#" xmlns:md="http://iec.ch/TC57/61970-552/ModelDescription/1#">
<dm:DifferenceModel rdf:about="urn:uuid:diff-model">
<md:Model.profile>http://entsoe.eu/CIM/EquipmentCore/3/1</md:Model.profile>
</dm:DifferenceModel>
<cim:ACLineSegment rdf:ID="_should_not_appear">
<cim:ACLineSegment.r>1.0</cim:ACLineSegment.r>
</cim:ACLineSegment>
</rdf:RDF>
"""

# Synthetic 3-bus star delivery (EQ+SSH+TP+SV): slack ENI at HV, loads at
# MV/LV, one 3W transformer whose end-2 RatioTapChanger carries an enabled
# TapChangerControl regulating the MV bus. Purpose-built for #294 point 4 —
# controllers on star-equivalent 3W legs — because no cached ENTSO-E set
# ships a controlled 3W changer. The target 114.4 kV (1.04 pu) is chosen so
# the controller MUST move the tap from the SSH neutral position (the
# uncontrolled solve lands near 1.0135 pu), and the 2.6 kV deadband spans a
# bit more than one 1.25 % step so discrete stepping can settle.
const _CGMES_SYNTH3W_EQ = """<?xml version="1.0" encoding="UTF-8"?>
<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:cim="http://iec.ch/TC57/2013/CIM-schema-cim16#" xmlns:md="http://iec.ch/TC57/61970-552/ModelDescription/1#">
<md:FullModel rdf:about="urn:uuid:eq3w"><md:Model.profile>http://entsoe.eu/CIM/EquipmentCore/3/1</md:Model.profile></md:FullModel>
<cim:BaseVoltage rdf:ID="_bv220"><cim:BaseVoltage.nominalVoltage>220</cim:BaseVoltage.nominalVoltage></cim:BaseVoltage>
<cim:BaseVoltage rdf:ID="_bv110"><cim:BaseVoltage.nominalVoltage>110</cim:BaseVoltage.nominalVoltage></cim:BaseVoltage>
<cim:BaseVoltage rdf:ID="_bv20"><cim:BaseVoltage.nominalVoltage>20</cim:BaseVoltage.nominalVoltage></cim:BaseVoltage>
<cim:ExternalNetworkInjection rdf:ID="_eni">
<cim:IdentifiedObject.name>GRID</cim:IdentifiedObject.name>
</cim:ExternalNetworkInjection>
<cim:Terminal rdf:ID="_t_eni"><cim:Terminal.ConductingEquipment rdf:resource="#_eni"/><cim:ACDCTerminal.sequenceNumber>1</cim:ACDCTerminal.sequenceNumber></cim:Terminal>
<cim:EnergyConsumer rdf:ID="_load_mv"><cim:IdentifiedObject.name>LOAD_MV</cim:IdentifiedObject.name></cim:EnergyConsumer>
<cim:Terminal rdf:ID="_t_load_mv"><cim:Terminal.ConductingEquipment rdf:resource="#_load_mv"/><cim:ACDCTerminal.sequenceNumber>1</cim:ACDCTerminal.sequenceNumber></cim:Terminal>
<cim:EnergyConsumer rdf:ID="_load_lv"><cim:IdentifiedObject.name>LOAD_LV</cim:IdentifiedObject.name></cim:EnergyConsumer>
<cim:Terminal rdf:ID="_t_load_lv"><cim:Terminal.ConductingEquipment rdf:resource="#_load_lv"/><cim:ACDCTerminal.sequenceNumber>1</cim:ACDCTerminal.sequenceNumber></cim:Terminal>
<cim:PowerTransformer rdf:ID="_t3w"><cim:IdentifiedObject.name>T3W</cim:IdentifiedObject.name></cim:PowerTransformer>
<cim:PowerTransformerEnd rdf:ID="_t3w_e1">
<cim:PowerTransformerEnd.PowerTransformer rdf:resource="#_t3w"/>
<cim:TransformerEnd.endNumber>1</cim:TransformerEnd.endNumber>
<cim:PowerTransformerEnd.ratedU>220</cim:PowerTransformerEnd.ratedU>
<cim:PowerTransformerEnd.ratedS>200</cim:PowerTransformerEnd.ratedS>
<cim:PowerTransformerEnd.r>0.4</cim:PowerTransformerEnd.r>
<cim:PowerTransformerEnd.x>24.0</cim:PowerTransformerEnd.x>
<cim:TransformerEnd.Terminal rdf:resource="#_t_t3w_1"/>
</cim:PowerTransformerEnd>
<cim:PowerTransformerEnd rdf:ID="_t3w_e2">
<cim:PowerTransformerEnd.PowerTransformer rdf:resource="#_t3w"/>
<cim:TransformerEnd.endNumber>2</cim:TransformerEnd.endNumber>
<cim:PowerTransformerEnd.ratedU>110</cim:PowerTransformerEnd.ratedU>
<cim:PowerTransformerEnd.ratedS>200</cim:PowerTransformerEnd.ratedS>
<cim:PowerTransformerEnd.r>0.1</cim:PowerTransformerEnd.r>
<cim:PowerTransformerEnd.x>6.0</cim:PowerTransformerEnd.x>
<cim:TransformerEnd.Terminal rdf:resource="#_t_t3w_2"/>
</cim:PowerTransformerEnd>
<cim:PowerTransformerEnd rdf:ID="_t3w_e3">
<cim:PowerTransformerEnd.PowerTransformer rdf:resource="#_t3w"/>
<cim:TransformerEnd.endNumber>3</cim:TransformerEnd.endNumber>
<cim:PowerTransformerEnd.ratedU>20</cim:PowerTransformerEnd.ratedU>
<cim:PowerTransformerEnd.ratedS>50</cim:PowerTransformerEnd.ratedS>
<cim:PowerTransformerEnd.r>0.02</cim:PowerTransformerEnd.r>
<cim:PowerTransformerEnd.x>0.8</cim:PowerTransformerEnd.x>
<cim:TransformerEnd.Terminal rdf:resource="#_t_t3w_3"/>
</cim:PowerTransformerEnd>
<cim:Terminal rdf:ID="_t_t3w_1"><cim:Terminal.ConductingEquipment rdf:resource="#_t3w"/><cim:ACDCTerminal.sequenceNumber>1</cim:ACDCTerminal.sequenceNumber></cim:Terminal>
<cim:Terminal rdf:ID="_t_t3w_2"><cim:Terminal.ConductingEquipment rdf:resource="#_t3w"/><cim:ACDCTerminal.sequenceNumber>2</cim:ACDCTerminal.sequenceNumber></cim:Terminal>
<cim:Terminal rdf:ID="_t_t3w_3"><cim:Terminal.ConductingEquipment rdf:resource="#_t3w"/><cim:ACDCTerminal.sequenceNumber>3</cim:ACDCTerminal.sequenceNumber></cim:Terminal>
<cim:RatioTapChanger rdf:ID="_rtc2">
<cim:IdentifiedObject.name>T3W_OLTC</cim:IdentifiedObject.name>
<cim:RatioTapChanger.TransformerEnd rdf:resource="#_t3w_e2"/>
<cim:RatioTapChanger.stepVoltageIncrement>1.25</cim:RatioTapChanger.stepVoltageIncrement>
<cim:TapChanger.lowStep>1</cim:TapChanger.lowStep>
<cim:TapChanger.highStep>21</cim:TapChanger.highStep>
<cim:TapChanger.neutralStep>11</cim:TapChanger.neutralStep>
<cim:TapChanger.TapChangerControl rdf:resource="#_tcc2"/>
</cim:RatioTapChanger>
<cim:TapChangerControl rdf:ID="_tcc2">
<cim:IdentifiedObject.name>T3W_OLTC_CTRL</cim:IdentifiedObject.name>
<cim:RegulatingControl.mode rdf:resource="http://iec.ch/TC57/2013/CIM-schema-cim16#RegulatingControlModeKind.voltage"/>
<cim:RegulatingControl.Terminal rdf:resource="#_t_load_mv"/>
</cim:TapChangerControl>
<cim:NonlinearShuntCompensator rdf:ID="_nlsh">
<cim:IdentifiedObject.name>NLSH</cim:IdentifiedObject.name>
</cim:NonlinearShuntCompensator>
<cim:Terminal rdf:ID="_t_nlsh"><cim:Terminal.ConductingEquipment rdf:resource="#_nlsh"/><cim:ACDCTerminal.sequenceNumber>1</cim:ACDCTerminal.sequenceNumber></cim:Terminal>
<cim:NonlinearShuntCompensatorPoint rdf:ID="_nlsh_p1">
<cim:NonlinearShuntCompensatorPoint.NonlinearShuntCompensator rdf:resource="#_nlsh"/>
<cim:NonlinearShuntCompensatorPoint.sectionNumber>1</cim:NonlinearShuntCompensatorPoint.sectionNumber>
<cim:NonlinearShuntCompensatorPoint.b>1.0e-4</cim:NonlinearShuntCompensatorPoint.b>
<cim:NonlinearShuntCompensatorPoint.g>0</cim:NonlinearShuntCompensatorPoint.g>
</cim:NonlinearShuntCompensatorPoint>
<cim:NonlinearShuntCompensatorPoint rdf:ID="_nlsh_p2">
<cim:NonlinearShuntCompensatorPoint.NonlinearShuntCompensator rdf:resource="#_nlsh"/>
<cim:NonlinearShuntCompensatorPoint.sectionNumber>2</cim:NonlinearShuntCompensatorPoint.sectionNumber>
<cim:NonlinearShuntCompensatorPoint.b>0.8e-4</cim:NonlinearShuntCompensatorPoint.b>
<cim:NonlinearShuntCompensatorPoint.g>0</cim:NonlinearShuntCompensatorPoint.g>
</cim:NonlinearShuntCompensatorPoint>
<cim:NonlinearShuntCompensatorPoint rdf:ID="_nlsh_p3">
<cim:NonlinearShuntCompensatorPoint.NonlinearShuntCompensator rdf:resource="#_nlsh"/>
<cim:NonlinearShuntCompensatorPoint.sectionNumber>3</cim:NonlinearShuntCompensatorPoint.sectionNumber>
<cim:NonlinearShuntCompensatorPoint.b>0.6e-4</cim:NonlinearShuntCompensatorPoint.b>
<cim:NonlinearShuntCompensatorPoint.g>0</cim:NonlinearShuntCompensatorPoint.g>
</cim:NonlinearShuntCompensatorPoint>
</rdf:RDF>
"""

const _CGMES_SYNTH3W_SSH = """<?xml version="1.0" encoding="UTF-8"?>
<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:cim="http://iec.ch/TC57/2013/CIM-schema-cim16#" xmlns:md="http://iec.ch/TC57/61970-552/ModelDescription/1#">
<md:FullModel rdf:about="urn:uuid:ssh3w"><md:Model.profile>http://entsoe.eu/CIM/SteadyStateHypothesis/1/1</md:Model.profile></md:FullModel>
<cim:ExternalNetworkInjection rdf:about="#_eni"><cim:ExternalNetworkInjection.p>-45</cim:ExternalNetworkInjection.p><cim:ExternalNetworkInjection.q>-12</cim:ExternalNetworkInjection.q></cim:ExternalNetworkInjection>
<cim:EnergyConsumer rdf:about="#_load_mv"><cim:EnergyConsumer.p>40</cim:EnergyConsumer.p><cim:EnergyConsumer.q>10</cim:EnergyConsumer.q></cim:EnergyConsumer>
<cim:EnergyConsumer rdf:about="#_load_lv"><cim:EnergyConsumer.p>5</cim:EnergyConsumer.p><cim:EnergyConsumer.q>1</cim:EnergyConsumer.q></cim:EnergyConsumer>
<cim:RatioTapChanger rdf:about="#_rtc2"><cim:TapChanger.step>11</cim:TapChanger.step><cim:TapChanger.controlEnabled>true</cim:TapChanger.controlEnabled></cim:RatioTapChanger>
<cim:TapChangerControl rdf:about="#_tcc2"><cim:RegulatingControl.enabled>true</cim:RegulatingControl.enabled><cim:RegulatingControl.targetValue>114.4</cim:RegulatingControl.targetValue><cim:RegulatingControl.targetDeadband>2.6</cim:RegulatingControl.targetDeadband></cim:TapChangerControl>
<cim:Terminal rdf:about="#_t_eni"><cim:ACDCTerminal.connected>true</cim:ACDCTerminal.connected></cim:Terminal>
<cim:Terminal rdf:about="#_t_load_mv"><cim:ACDCTerminal.connected>true</cim:ACDCTerminal.connected></cim:Terminal>
<cim:Terminal rdf:about="#_t_load_lv"><cim:ACDCTerminal.connected>true</cim:ACDCTerminal.connected></cim:Terminal>
<cim:Terminal rdf:about="#_t_t3w_1"><cim:ACDCTerminal.connected>true</cim:ACDCTerminal.connected></cim:Terminal>
<cim:Terminal rdf:about="#_t_t3w_2"><cim:ACDCTerminal.connected>true</cim:ACDCTerminal.connected></cim:Terminal>
<cim:Terminal rdf:about="#_t_t3w_3"><cim:ACDCTerminal.connected>true</cim:ACDCTerminal.connected></cim:Terminal>
<cim:NonlinearShuntCompensator rdf:about="#_nlsh"><cim:ShuntCompensator.sections>2</cim:ShuntCompensator.sections></cim:NonlinearShuntCompensator>
<cim:Terminal rdf:about="#_t_nlsh"><cim:ACDCTerminal.connected>true</cim:ACDCTerminal.connected></cim:Terminal>
</rdf:RDF>
"""

const _CGMES_SYNTH3W_TP = """<?xml version="1.0" encoding="UTF-8"?>
<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:cim="http://iec.ch/TC57/2013/CIM-schema-cim16#" xmlns:md="http://iec.ch/TC57/61970-552/ModelDescription/1#">
<md:FullModel rdf:about="urn:uuid:tp3w"><md:Model.profile>http://entsoe.eu/CIM/Topology/4/1</md:Model.profile></md:FullModel>
<cim:TopologicalNode rdf:ID="_tn_hv"><cim:IdentifiedObject.name>HV_BUS</cim:IdentifiedObject.name><cim:TopologicalNode.BaseVoltage rdf:resource="#_bv220"/></cim:TopologicalNode>
<cim:TopologicalNode rdf:ID="_tn_mv"><cim:IdentifiedObject.name>MV_BUS</cim:IdentifiedObject.name><cim:TopologicalNode.BaseVoltage rdf:resource="#_bv110"/></cim:TopologicalNode>
<cim:TopologicalNode rdf:ID="_tn_lv"><cim:IdentifiedObject.name>LV_BUS</cim:IdentifiedObject.name><cim:TopologicalNode.BaseVoltage rdf:resource="#_bv20"/></cim:TopologicalNode>
<cim:Terminal rdf:about="#_t_eni"><cim:Terminal.TopologicalNode rdf:resource="#_tn_hv"/></cim:Terminal>
<cim:Terminal rdf:about="#_t_load_mv"><cim:Terminal.TopologicalNode rdf:resource="#_tn_mv"/></cim:Terminal>
<cim:Terminal rdf:about="#_t_load_lv"><cim:Terminal.TopologicalNode rdf:resource="#_tn_lv"/></cim:Terminal>
<cim:Terminal rdf:about="#_t_t3w_1"><cim:Terminal.TopologicalNode rdf:resource="#_tn_hv"/></cim:Terminal>
<cim:Terminal rdf:about="#_t_t3w_2"><cim:Terminal.TopologicalNode rdf:resource="#_tn_mv"/></cim:Terminal>
<cim:Terminal rdf:about="#_t_t3w_3"><cim:Terminal.TopologicalNode rdf:resource="#_tn_lv"/></cim:Terminal>
<cim:Terminal rdf:about="#_t_nlsh"><cim:Terminal.TopologicalNode rdf:resource="#_tn_mv"/></cim:Terminal>
</rdf:RDF>
"""

const _CGMES_SYNTH3W_SV = """<?xml version="1.0" encoding="UTF-8"?>
<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:cim="http://iec.ch/TC57/2013/CIM-schema-cim16#" xmlns:md="http://iec.ch/TC57/61970-552/ModelDescription/1#">
<md:FullModel rdf:about="urn:uuid:sv3w"><md:Model.profile>http://entsoe.eu/CIM/StateVariables/4/1</md:Model.profile></md:FullModel>
<cim:SvVoltage rdf:ID="_sv_hv"><cim:SvVoltage.TopologicalNode rdf:resource="#_tn_hv"/><cim:SvVoltage.v>220</cim:SvVoltage.v><cim:SvVoltage.angle>0</cim:SvVoltage.angle></cim:SvVoltage>
<cim:SvVoltage rdf:ID="_sv_mv"><cim:SvVoltage.TopologicalNode rdf:resource="#_tn_mv"/><cim:SvVoltage.v>111.5</cim:SvVoltage.v><cim:SvVoltage.angle>-1.5</cim:SvVoltage.angle></cim:SvVoltage>
<cim:SvVoltage rdf:ID="_sv_lv"><cim:SvVoltage.TopologicalNode rdf:resource="#_tn_lv"/><cim:SvVoltage.v>19.9</cim:SvVoltage.v><cim:SvVoltage.angle>-2.0</cim:SvVoltage.angle></cim:SvVoltage>
<cim:TopologicalIsland rdf:ID="_sv_isl"><cim:IdentifiedObject.name>ISL</cim:IdentifiedObject.name><cim:TopologicalIsland.TopologicalNodes rdf:resource="#_tn_hv"/><cim:TopologicalIsland.TopologicalNodes rdf:resource="#_tn_mv"/><cim:TopologicalIsland.TopologicalNodes rdf:resource="#_tn_lv"/></cim:TopologicalIsland>
</rdf:RDF>
"""

function _cgmes_synth3w_dir()::String
  dir = mktempdir()
  write(joinpath(dir, "synth3w_EQ.xml"), _CGMES_SYNTH3W_EQ)
  write(joinpath(dir, "synth3w_SSH.xml"), _CGMES_SYNTH3W_SSH)
  write(joinpath(dir, "synth3w_TP.xml"), _CGMES_SYNTH3W_TP)
  write(joinpath(dir, "synth3w_SV.xml"), _CGMES_SYNTH3W_SV)
  return dir
end

# Synthetic 3-bus chain (EQ+SSH+TP+SV) for #294 point 3 — remote-regulating
# machines: slack ENI at bus A, machine G_B at bus B whose voltage
# RegulatingControl terminal sits at the LOAD bus C (112.2 kV = 1.02 pu),
# load at C. A second machine G_C at bus C regulates locally but is parked
# behind a disconnected terminal; the `sm2` toggle reconnects it so the same
# delivery also exercises the "target bus already voltage-held" fallback.
const _CGMES_SYNTH_RVC_EQ = """<?xml version="1.0" encoding="UTF-8"?>
<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:cim="http://iec.ch/TC57/2013/CIM-schema-cim16#" xmlns:md="http://iec.ch/TC57/61970-552/ModelDescription/1#">
<md:FullModel rdf:about="urn:uuid:eqrvc"><md:Model.profile>http://entsoe.eu/CIM/EquipmentCore/3/1</md:Model.profile></md:FullModel>
<cim:BaseVoltage rdf:ID="_bv110"><cim:BaseVoltage.nominalVoltage>110</cim:BaseVoltage.nominalVoltage></cim:BaseVoltage>
<cim:ExternalNetworkInjection rdf:ID="_eni">
<cim:IdentifiedObject.name>GRID</cim:IdentifiedObject.name>
<cim:ExternalNetworkInjection.referencePriority>1</cim:ExternalNetworkInjection.referencePriority>
<cim:ExternalNetworkInjection.maxP>200</cim:ExternalNetworkInjection.maxP>
</cim:ExternalNetworkInjection>
<cim:Terminal rdf:ID="_t_eni"><cim:Terminal.ConductingEquipment rdf:resource="#_eni"/><cim:ACDCTerminal.sequenceNumber>1</cim:ACDCTerminal.sequenceNumber></cim:Terminal>
<cim:SynchronousMachine rdf:ID="_smb">
<cim:IdentifiedObject.name>G_B</cim:IdentifiedObject.name>
<cim:RotatingMachine.ratedS>100</cim:RotatingMachine.ratedS>
<cim:SynchronousMachine.minQ>-50</cim:SynchronousMachine.minQ>
<cim:SynchronousMachine.maxQ>50</cim:SynchronousMachine.maxQ>
<cim:RegulatingCondEq.RegulatingControl rdf:resource="#_rc_smb"/>
</cim:SynchronousMachine>
<cim:Terminal rdf:ID="_t_smb"><cim:Terminal.ConductingEquipment rdf:resource="#_smb"/><cim:ACDCTerminal.sequenceNumber>1</cim:ACDCTerminal.sequenceNumber></cim:Terminal>
<cim:RegulatingControl rdf:ID="_rc_smb">
<cim:IdentifiedObject.name>G_B_RC</cim:IdentifiedObject.name>
<cim:RegulatingControl.mode rdf:resource="http://iec.ch/TC57/2013/CIM-schema-cim16#RegulatingControlModeKind.voltage"/>
<cim:RegulatingControl.Terminal rdf:resource="#_t_load"/>
</cim:RegulatingControl>
<cim:SynchronousMachine rdf:ID="_smc">
<cim:IdentifiedObject.name>G_C</cim:IdentifiedObject.name>
<cim:RotatingMachine.ratedS>60</cim:RotatingMachine.ratedS>
<cim:SynchronousMachine.minQ>-30</cim:SynchronousMachine.minQ>
<cim:SynchronousMachine.maxQ>30</cim:SynchronousMachine.maxQ>
<cim:RegulatingCondEq.RegulatingControl rdf:resource="#_rc_smc"/>
</cim:SynchronousMachine>
<cim:Terminal rdf:ID="_t_smc"><cim:Terminal.ConductingEquipment rdf:resource="#_smc"/><cim:ACDCTerminal.sequenceNumber>1</cim:ACDCTerminal.sequenceNumber></cim:Terminal>
<cim:RegulatingControl rdf:ID="_rc_smc">
<cim:IdentifiedObject.name>G_C_RC</cim:IdentifiedObject.name>
<cim:RegulatingControl.mode rdf:resource="http://iec.ch/TC57/2013/CIM-schema-cim16#RegulatingControlModeKind.voltage"/>
<cim:RegulatingControl.Terminal rdf:resource="#_t_smc"/>
</cim:RegulatingControl>
<cim:EnergyConsumer rdf:ID="_load"><cim:IdentifiedObject.name>LOAD_C</cim:IdentifiedObject.name></cim:EnergyConsumer>
<cim:Terminal rdf:ID="_t_load"><cim:Terminal.ConductingEquipment rdf:resource="#_load"/><cim:ACDCTerminal.sequenceNumber>1</cim:ACDCTerminal.sequenceNumber></cim:Terminal>
<cim:ACLineSegment rdf:ID="_l_ab">
<cim:IdentifiedObject.name>L_AB</cim:IdentifiedObject.name>
<cim:ACLineSegment.r>2.42</cim:ACLineSegment.r>
<cim:ACLineSegment.x>12.1</cim:ACLineSegment.x>
<cim:ACLineSegment.bch>0</cim:ACLineSegment.bch>
<cim:ConductingEquipment.BaseVoltage rdf:resource="#_bv110"/>
</cim:ACLineSegment>
<cim:Terminal rdf:ID="_t_lab_1"><cim:Terminal.ConductingEquipment rdf:resource="#_l_ab"/><cim:ACDCTerminal.sequenceNumber>1</cim:ACDCTerminal.sequenceNumber></cim:Terminal>
<cim:Terminal rdf:ID="_t_lab_2"><cim:Terminal.ConductingEquipment rdf:resource="#_l_ab"/><cim:ACDCTerminal.sequenceNumber>2</cim:ACDCTerminal.sequenceNumber></cim:Terminal>
<cim:ACLineSegment rdf:ID="_l_bc">
<cim:IdentifiedObject.name>L_BC</cim:IdentifiedObject.name>
<cim:ACLineSegment.r>2.42</cim:ACLineSegment.r>
<cim:ACLineSegment.x>12.1</cim:ACLineSegment.x>
<cim:ACLineSegment.bch>0</cim:ACLineSegment.bch>
<cim:ConductingEquipment.BaseVoltage rdf:resource="#_bv110"/>
</cim:ACLineSegment>
<cim:Terminal rdf:ID="_t_lbc_1"><cim:Terminal.ConductingEquipment rdf:resource="#_l_bc"/><cim:ACDCTerminal.sequenceNumber>1</cim:ACDCTerminal.sequenceNumber></cim:Terminal>
<cim:Terminal rdf:ID="_t_lbc_2"><cim:Terminal.ConductingEquipment rdf:resource="#_l_bc"/><cim:ACDCTerminal.sequenceNumber>2</cim:ACDCTerminal.sequenceNumber></cim:Terminal>
</rdf:RDF>
"""

const _CGMES_SYNTH_RVC_SSH_TEMPLATE = """<?xml version="1.0" encoding="UTF-8"?>
<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:cim="http://iec.ch/TC57/2013/CIM-schema-cim16#" xmlns:md="http://iec.ch/TC57/61970-552/ModelDescription/1#">
<md:FullModel rdf:about="urn:uuid:sshrvc"><md:Model.profile>http://entsoe.eu/CIM/SteadyStateHypothesis/1/1</md:Model.profile></md:FullModel>
<cim:ExternalNetworkInjection rdf:about="#_eni"><cim:ExternalNetworkInjection.p>-42</cim:ExternalNetworkInjection.p><cim:ExternalNetworkInjection.q>-10</cim:ExternalNetworkInjection.q></cim:ExternalNetworkInjection>
<cim:SynchronousMachine rdf:about="#_smb"><cim:RotatingMachine.p>-20</cim:RotatingMachine.p><cim:RotatingMachine.q>0</cim:RotatingMachine.q></cim:SynchronousMachine>
<cim:SynchronousMachine rdf:about="#_smc"><cim:RotatingMachine.p>-5</cim:RotatingMachine.p><cim:RotatingMachine.q>0</cim:RotatingMachine.q></cim:SynchronousMachine>
<cim:RegulatingControl rdf:about="#_rc_smb"><cim:RegulatingControl.enabled>true</cim:RegulatingControl.enabled><cim:RegulatingControl.targetValue>112.2</cim:RegulatingControl.targetValue></cim:RegulatingControl>
<cim:RegulatingControl rdf:about="#_rc_smc"><cim:RegulatingControl.enabled>true</cim:RegulatingControl.enabled><cim:RegulatingControl.targetValue>112.2</cim:RegulatingControl.targetValue></cim:RegulatingControl>
<cim:EnergyConsumer rdf:about="#_load"><cim:EnergyConsumer.p>60</cim:EnergyConsumer.p><cim:EnergyConsumer.q>15</cim:EnergyConsumer.q></cim:EnergyConsumer>
<cim:Terminal rdf:about="#_t_eni"><cim:ACDCTerminal.connected>true</cim:ACDCTerminal.connected></cim:Terminal>
<cim:Terminal rdf:about="#_t_smb"><cim:ACDCTerminal.connected>true</cim:ACDCTerminal.connected></cim:Terminal>
<cim:Terminal rdf:about="#_t_smc"><cim:ACDCTerminal.connected>__SM2CONN__</cim:ACDCTerminal.connected></cim:Terminal>
<cim:Terminal rdf:about="#_t_load"><cim:ACDCTerminal.connected>true</cim:ACDCTerminal.connected></cim:Terminal>
<cim:Terminal rdf:about="#_t_lab_1"><cim:ACDCTerminal.connected>true</cim:ACDCTerminal.connected></cim:Terminal>
<cim:Terminal rdf:about="#_t_lab_2"><cim:ACDCTerminal.connected>true</cim:ACDCTerminal.connected></cim:Terminal>
<cim:Terminal rdf:about="#_t_lbc_1"><cim:ACDCTerminal.connected>true</cim:ACDCTerminal.connected></cim:Terminal>
<cim:Terminal rdf:about="#_t_lbc_2"><cim:ACDCTerminal.connected>true</cim:ACDCTerminal.connected></cim:Terminal>
</rdf:RDF>
"""

const _CGMES_SYNTH_RVC_TP = """<?xml version="1.0" encoding="UTF-8"?>
<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:cim="http://iec.ch/TC57/2013/CIM-schema-cim16#" xmlns:md="http://iec.ch/TC57/61970-552/ModelDescription/1#">
<md:FullModel rdf:about="urn:uuid:tprvc"><md:Model.profile>http://entsoe.eu/CIM/Topology/4/1</md:Model.profile></md:FullModel>
<cim:TopologicalNode rdf:ID="_tn_a"><cim:IdentifiedObject.name>BUS_A</cim:IdentifiedObject.name><cim:TopologicalNode.BaseVoltage rdf:resource="#_bv110"/></cim:TopologicalNode>
<cim:TopologicalNode rdf:ID="_tn_b"><cim:IdentifiedObject.name>BUS_B</cim:IdentifiedObject.name><cim:TopologicalNode.BaseVoltage rdf:resource="#_bv110"/></cim:TopologicalNode>
<cim:TopologicalNode rdf:ID="_tn_c"><cim:IdentifiedObject.name>BUS_C</cim:IdentifiedObject.name><cim:TopologicalNode.BaseVoltage rdf:resource="#_bv110"/></cim:TopologicalNode>
<cim:Terminal rdf:about="#_t_eni"><cim:Terminal.TopologicalNode rdf:resource="#_tn_a"/></cim:Terminal>
<cim:Terminal rdf:about="#_t_smb"><cim:Terminal.TopologicalNode rdf:resource="#_tn_b"/></cim:Terminal>
<cim:Terminal rdf:about="#_t_smc"><cim:Terminal.TopologicalNode rdf:resource="#_tn_c"/></cim:Terminal>
<cim:Terminal rdf:about="#_t_load"><cim:Terminal.TopologicalNode rdf:resource="#_tn_c"/></cim:Terminal>
<cim:Terminal rdf:about="#_t_lab_1"><cim:Terminal.TopologicalNode rdf:resource="#_tn_a"/></cim:Terminal>
<cim:Terminal rdf:about="#_t_lab_2"><cim:Terminal.TopologicalNode rdf:resource="#_tn_b"/></cim:Terminal>
<cim:Terminal rdf:about="#_t_lbc_1"><cim:Terminal.TopologicalNode rdf:resource="#_tn_b"/></cim:Terminal>
<cim:Terminal rdf:about="#_t_lbc_2"><cim:Terminal.TopologicalNode rdf:resource="#_tn_c"/></cim:Terminal>
</rdf:RDF>
"""

const _CGMES_SYNTH_RVC_SV = """<?xml version="1.0" encoding="UTF-8"?>
<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:cim="http://iec.ch/TC57/2013/CIM-schema-cim16#" xmlns:md="http://iec.ch/TC57/61970-552/ModelDescription/1#">
<md:FullModel rdf:about="urn:uuid:svrvc"><md:Model.profile>http://entsoe.eu/CIM/StateVariables/4/1</md:Model.profile></md:FullModel>
<cim:SvVoltage rdf:ID="_sv_a"><cim:SvVoltage.TopologicalNode rdf:resource="#_tn_a"/><cim:SvVoltage.v>113.3</cim:SvVoltage.v><cim:SvVoltage.angle>0</cim:SvVoltage.angle></cim:SvVoltage>
<cim:SvVoltage rdf:ID="_sv_b"><cim:SvVoltage.TopologicalNode rdf:resource="#_tn_b"/><cim:SvVoltage.v>112.5</cim:SvVoltage.v><cim:SvVoltage.angle>-0.8</cim:SvVoltage.angle></cim:SvVoltage>
<cim:SvVoltage rdf:ID="_sv_c"><cim:SvVoltage.TopologicalNode rdf:resource="#_tn_c"/><cim:SvVoltage.v>111.8</cim:SvVoltage.v><cim:SvVoltage.angle>-1.4</cim:SvVoltage.angle></cim:SvVoltage>
</rdf:RDF>
"""

function _cgmes_synth_rvc_dir(; sm2_connected::Bool = false)::String
  dir = mktempdir()
  write(joinpath(dir, "synthrvc_EQ.xml"), _CGMES_SYNTH_RVC_EQ)
  write(joinpath(dir, "synthrvc_SSH.xml"), replace(_CGMES_SYNTH_RVC_SSH_TEMPLATE, "__SM2CONN__" => sm2_connected ? "true" : "false"))
  write(joinpath(dir, "synthrvc_TP.xml"), _CGMES_SYNTH_RVC_TP)
  write(joinpath(dir, "synthrvc_SV.xml"), _CGMES_SYNTH_RVC_SV)
  return dir
end

function _cgmes_synth_store()::CGMESStore
  dir = mktempdir()
  write(joinpath(dir, "synth_EQ.xml"), _CGMES_SYNTH_EQ)
  write(joinpath(dir, "synth_SSH.xml"), _CGMES_SYNTH_SSH)
  write(joinpath(dir, "synth_DIFF.xml"), _CGMES_SYNTH_DIFF)
  return loadCGMES(dir)
end

function run_cgmes_importer_tests()
  @testset "CGMES importer" begin
    store = _cgmes_synth_store()

    @testset "profile classification and version" begin
      @test store.version == "2.4.15"
      eqinfo = only(filter(f -> occursin("EQ", f.name), store.files))
      @test eqinfo.header == :FullModel
      @test eqinfo.profiles == Set([:EQ, :EQ_SC])   # per-file profile *set* (A2/D-7)
      @test !eqinfo.skipped
    end

    @testset "DifferenceModel files are skipped with reason (D-6)" begin
      diffinfo = only(filter(f -> occursin("DIFF", f.name), store.files))
      @test diffinfo.skipped
      @test diffinfo.header == :DifferenceModel
      @test !isempty(diffinfo.skip_reason)
      @test !haskey(store.objects, "_should_not_appear")
    end

    @testset "import failure analysis names supplied models and gaps" begin
      # The synthetic store carries a model id per file and known unresolved
      # references (_missing_container, the island's _tn_* members).
      report = importFailureAnalysis(store)
      @test occursin("Supplied models:", report)
      @test occursin("eq-model", report)
      @test occursin("Unresolved references:", report)
      @test occursin("Verdict:", report)

      # A file declaring a prerequisite that is not part of the input must be
      # called out with the exact missing model id — the boundary-set case.
      depdir = mktempdir()
      dep_eq = replace(
        _CGMES_SYNTH_EQ,
        "<md:Model.profile>http://entsoe.eu/CIM/EquipmentCore/3/1</md:Model.profile>" =>
          "<md:Model.profile>http://entsoe.eu/CIM/EquipmentCore/3/1</md:Model.profile>\n<md:Model.DependentOn rdf:resource=\"urn:uuid:boundary-model-123\"/>",
      )
      write(joinpath(depdir, "dep_EQ.xml"), dep_eq)
      dep_report = importFailureAnalysis(loadCGMES(depdir))
      @test occursin("MISSING", dep_report)
      @test occursin("boundary-model-123", dep_report)
      @test occursin("prerequisite model(s)", dep_report)

      # A dependency satisfied by a supplied file is not reported as missing.
      satdir = mktempdir()
      sat_eq = replace(
        _CGMES_SYNTH_EQ,
        "<md:Model.profile>http://entsoe.eu/CIM/EquipmentCore/3/1</md:Model.profile>" =>
          "<md:Model.profile>http://entsoe.eu/CIM/EquipmentCore/3/1</md:Model.profile>\n<md:Model.DependentOn rdf:resource=\"urn:uuid:ssh-model\"/>",
      )
      write(joinpath(satdir, "sat_EQ.xml"), sat_eq)
      write(joinpath(satdir, "sat_SSH.xml"), _CGMES_SYNTH_SSH)
      sat_report = importFailureAnalysis(loadCGMES(satdir))
      @test occursin("Declared dependencies: all satisfied", sat_report)
    end

    @testset "infer_base_voltages reconstructs missing nominal voltages" begin
      # The RVC delivery with its BaseVoltage catalog stripped — the exact
      # shape of a real delivery whose catalog lives in a missing boundary EQ.
      nobv_dir = mktempdir()
      write(joinpath(nobv_dir, "nobv_EQ.xml"), replace(_CGMES_SYNTH_RVC_EQ, "<cim:BaseVoltage rdf:ID=\"_bv110\"><cim:BaseVoltage.nominalVoltage>110</cim:BaseVoltage.nominalVoltage></cim:BaseVoltage>" => ""))
      write(joinpath(nobv_dir, "nobv_TP.xml"), replace(_CGMES_SYNTH_RVC_TP, "<cim:TopologicalNode.BaseVoltage rdf:resource=\"#_bv110\"/>" => ""))
      write(joinpath(nobv_dir, "nobv_SSH.xml"), replace(_CGMES_SYNTH_RVC_SSH_TEMPLATE, "__SM2CONN__" => "false"))
      write(joinpath(nobv_dir, "nobv_SV.xml"), _CGMES_SYNTH_RVC_SV)

      # Without the option the import aborts with the typed analysis error.
      @test_throws Sparlectra.CGMESImporter.CGMESImportError importCGMES(path = nobv_dir, name = "nobv_off")

      # With the option every bus level is reconstructed from the SV state
      # (111.8–115.5 kV snap to 110), the substitution is summarized as one
      # warning, and the network still solves.
      res = importCGMES(path = nobv_dir, name = "nobv_on", infer_base_voltages = true)
      @test length(res.net.nodeVec) == 3
      @test all(Sparlectra.getNodeVn(b) == 110.0 for b in res.net.nodeVec)
      @test count(m -> occursin("inferred base voltages", m), res.messages) == 1
      @test last(runpf!(res.net, 30, 1e-8, 0)) == 0
    end

    @testset "self-loop line maps to a shunt notice, not a branch" begin
      # Both terminals of a line on ONE topological node (busbar link modeled
      # as a line) used to trip the branch constructor's from!=to assertion.
      loop_dir = mktempdir()
      loop_eq = replace(
        _CGMES_SYNTH_RVC_EQ,
        "<cim:ACLineSegment rdf:ID=\"_l_ab\">" =>
          "<cim:ACLineSegment rdf:ID=\"_l_loop\"><cim:IdentifiedObject.name>L_LOOP</cim:IdentifiedObject.name><cim:ACLineSegment.r>0.01</cim:ACLineSegment.r><cim:ACLineSegment.x>0.05</cim:ACLineSegment.x></cim:ACLineSegment>\n<cim:Terminal rdf:ID=\"_t_loop_1\"><cim:Terminal.ConductingEquipment rdf:resource=\"#_l_loop\"/><cim:ACDCTerminal.sequenceNumber>1</cim:ACDCTerminal.sequenceNumber></cim:Terminal>\n<cim:Terminal rdf:ID=\"_t_loop_2\"><cim:Terminal.ConductingEquipment rdf:resource=\"#_l_loop\"/><cim:ACDCTerminal.sequenceNumber>2</cim:ACDCTerminal.sequenceNumber></cim:Terminal>\n<cim:ACLineSegment rdf:ID=\"_l_ab\">",
      )
      loop_tp = replace(
        _CGMES_SYNTH_RVC_TP,
        "<cim:Terminal rdf:about=\"#_t_lab_1\">" =>
          "<cim:Terminal rdf:about=\"#_t_loop_1\"><cim:Terminal.TopologicalNode rdf:resource=\"#_tn_a\"/></cim:Terminal>\n<cim:Terminal rdf:about=\"#_t_loop_2\"><cim:Terminal.TopologicalNode rdf:resource=\"#_tn_a\"/></cim:Terminal>\n<cim:Terminal rdf:about=\"#_t_lab_1\">",
      )
      write(joinpath(loop_dir, "loop_EQ.xml"), loop_eq)
      write(joinpath(loop_dir, "loop_TP.xml"), loop_tp)
      write(joinpath(loop_dir, "loop_SSH.xml"), replace(_CGMES_SYNTH_RVC_SSH_TEMPLATE, "__SM2CONN__" => "false"))
      write(joinpath(loop_dir, "loop_SV.xml"), _CGMES_SYNTH_RVC_SV)
      loop_res = importCGMES(path = loop_dir, name = "self_loop")
      @test count(m -> occursin("connects bus", m) && occursin("to itself", m), loop_res.messages) == 1
      # only the two real lines become branches
      @test length(loop_res.net.linesAC) == 2
      @test last(runpf!(loop_res.net, 30, 1e-8, 0)) == 0
    end

    @testset "rdf:ID creates, literals and inherited attributes" begin
      @test countOf(store, :ACLineSegment) == 1
      line = only(objectsOf(store, :ACLineSegment))
      @test str(line, :name) == "L1"                 # IdentifiedObject.name → :name
      @test num(line, :r) == 2.5
      @test num(line, :gch) == 0.001                 # CGMES conductance present
      @test num(line, :r0) == 7.5                    # short-circuit attribute read (§7.7)
      @test num(line, :length) == 12.0               # Conductor.length → :length
      @test num(line, :missing_attr, 99.0) == 99.0
    end

    @testset "references and enums" begin
      line = only(objectsOf(store, :ACLineSegment))
      bv = ref(store, line, :BaseVoltage)
      @test bv !== nothing && bv.class == :BaseVoltage
      @test num(bv, :nominalVoltage) == 110.0
      sm = only(objectsOf(store, :SynchronousMachine))
      @test enumval(sm, :operatingMode) == "SynchronousMachineOperatingMode.generator"
      @test boolval(sm, :earthing) == false
    end

    @testset "multi-valued references overflow into refsAll (#294 point 9)" begin
      island = only(objectsOf(store, :TopologicalIsland))
      # legacy shape untouched: refs carries the FIRST value under the short
      # key (plus the dotted-name fallback), never a silent last-wins
      @test island.refs[:TopologicalNodes] == "_tn_a"
      @test Sparlectra.CGMESImporter.refsAll(store, island, :TopologicalNodes) == ["_tn_a", "_tn_b", "_tn_c"]
      # single-valued references stay single, absent ones come back empty
      line = only(objectsOf(store, :ACLineSegment))
      @test Sparlectra.CGMESImporter.refsAll(store, line, :BaseVoltage) == ["_bv110"]
      @test Sparlectra.CGMESImporter.refsAll(store, line, :NoSuchRef) == String[]
    end

    @testset "rdf:about overlays EQ object (SSH)" begin
      sm = only(objectsOf(store, :SynchronousMachine))
      @test num(sm, :p) == -120.0                    # from SSH overlay
      @test num(sm, :ratedS) == 150.0                # EQ value kept
      @test countOf(store, :SynchronousMachine) == 1 # no duplicate object
    end

    @testset "unresolved references reported" begin
      unresolved = unresolvedReferences(store)
      @test any(u -> u.target == "_missing_container", unresolved)
    end

    @testset "summarizeCGMES on synthetic set" begin
      dir = mktempdir()
      write(joinpath(dir, "synth_EQ.xml"), _CGMES_SYNTH_EQ)
      s = summarizeCGMES(path = dir)
      @test s.object_count == 4
      @test s.unresolved_count >= 1
      @test !s.boundary_missing_hint                 # dangling ref is not a topology class
      @test (:ACLineSegment => 1) in s.class_histogram
    end

    @testset "ReactiveCapabilityCurve interpolation (#294 point 1)" begin
      # MicroGrid BE-G1 shape: P −100 → ±200, P 0 → ±300, P 100 → ±200 MVAr.
      pts = [(-100.0, -200.0, 200.0), (0.0, -300.0, 300.0), (100.0, -200.0, 200.0)]
      # interior interpolation (the machine's actual operating point)
      @test Sparlectra.CGMESImporter._curveQHull(pts, -90.0) == (-210.0, 210.0)
      @test Sparlectra.CGMESImporter._curveQHull(pts, 0.0) == (-300.0, 300.0)
      @test Sparlectra.CGMESImporter._curveQHull(pts, 50.0) == (-250.0, 250.0)
      # clamping outside the curve's P domain
      @test Sparlectra.CGMESImporter._curveQHull(pts, -500.0) == (-200.0, 200.0)
      @test Sparlectra.CGMESImporter._curveQHull(pts, 500.0) == (-200.0, 200.0)
      # the sign-convention hull: a curve written in the load convention
      # (y1/y2 swapped in sign) yields the same band
      swapped = [(-100.0, 200.0, -200.0), (0.0, 300.0, -300.0), (100.0, 200.0, -200.0)]
      @test Sparlectra.CGMESImporter._curveQHull(swapped, -90.0) == (-210.0, 210.0)
      # degenerate inputs are refused, the caller falls back to scalars
      @test Sparlectra.CGMESImporter._curveQHull(NTuple{3,Float64}[], 0.0) === nothing
      @test Sparlectra.CGMESImporter._curveQHull([(0.0, 0.0, 0.0)], 0.0) === nothing
      # single-point curve with a real range works
      @test Sparlectra.CGMESImporter._curveQHull([(0.0, -50.0, 50.0)], 25.0) == (-50.0, 50.0)
    end

    @testset "Result tables stay aligned with long names" begin
      # CGMES bus and branch identifiers routinely exceed the column widths;
      # @sprintf pads but never truncates, so without fitting every following
      # column would shift (reported from a MicroGrid Assembled run).
      @test Sparlectra._fitColumn("short", 25) == "short"
      @test length(Sparlectra._fitColumn("TN_Border_ST23 -> BE-Busbar_2", 25)) == 25
      @test endswith(Sparlectra._fitColumn("TN_Border_ST23 -> BE-Busbar_2", 25), "…")
      @test Sparlectra._fitColumn("exactly_twenty_five_chars", 25) == "exactly_twenty_five_chars"
    end

    # The Web UI docs reader serves an allowlist — the CGMES page and the
    # contextual help for its options must be reachable from the interface.
    @testset "Web UI documentation wiring" begin
      page = Sparlectra.resolve_webui_doc_page("cgmes_import")
      @test page !== nothing && page.file == "cgmes_import.md"
      @test isfile(joinpath(dirname(@__DIR__), "docs", "src", page.file))
      # The page must carry the option reference the docs reader links to.
      text = read(joinpath(dirname(@__DIR__), "docs", "src", page.file), String)
      for key in ("cgmes_import.path", "cgmes_import.base_mva", "cgmes_import.require_boundary", "cgmes_import.tap_control", "cgmes_import.ignore_connected")
        @test occursin(key, text)
      end
    end

    # A single ReliCapGrid model is one area of a multi-area system, so its
    # border nodes hang free when imported alone. The combined aliases fetch
    # several areas plus their shared boundary files as ONE delivery. Table
    # consistency is checked here without touching the network; the fetch itself
    # needs GitHub and is exercised by the fixture block below when cached.
    @testset "ReliCapGrid combined aliases" begin
      combined = Sparlectra.CGMESImporter.RELICAPGRID_COMBINED
      singles = Sparlectra.CGMESImporter.RELICAPGRID_ALIASES
      @test haskey(combined, "relicapgrid_cgm")
      @test haskey(combined, "svedala_neighbours")
      # Every member must be a known single model, else the fetch would fail
      # only at download time.
      for (name, members) in combined
        @test !isempty(members)
        for m in members
          @test haskey(singles, m)
        end
      end
      # relicapgrid_cgm spans every single model we offer.
      @test sort(combined["relicapgrid_cgm"]) == sort(collect(keys(singles)))

      # A border file is named Boundary_Border-<AreaA>-<AreaB>.xml. The border is
      # closed by a combination when BOTH its areas are members; otherwise the
      # nodes on that border still hang free.
      border_areas(file) = split(replace(replace(file, "Boundary_Border-" => ""), ".xml" => ""), "-")
      function open_borders(members)
        model_names = Set(lowercase(singles[m].model) for m in members)
        borders = Set{String}()
        for m in members
          union!(borders, singles[m].boundary)
        end
        return sort([b for b in borders if !all(a -> lowercase(a) in model_names, border_areas(b))])
      end

      # The whole point of the CGM alias: not a single border left open.
      @test isempty(open_borders(combined["relicapgrid_cgm"]))
      # The cheap variant closes Svedala's own two borders but inherits its
      # neighbours' outer borders — documented, not accidental.
      @test open_borders(combined["svedala_neighbours"]) == ["Boundary_Border-Espheim-Portheim.xml", "Boundary_Border-Galia-Belgovia.xml"]
      # Svedala alone: both of its borders hang free. That is the situation the
      # combined aliases exist to fix.
      @test length(open_borders(["svedala"])) == 2
      # Aliases must be discoverable through the public alias list.
      all_aliases = Sparlectra.CGMESImporter.allCGMESTestSetAliases()
      @test "relicapgrid_cgm" in all_aliases
      @test "svedala_neighbours" in all_aliases
      @test "portheim" in all_aliases
      # convenience shorthand + the MiniGrid pair (regression: cgmes:minigrid
      # was not resolvable although the set ships with the conformity package)
      @test "microgrid" in all_aliases
      @test "minigrid" in all_aliases
      @test "minigrid_nb" in all_aliases
      @test Sparlectra.CGMESImporter.CGMES_TESTSET_ALIASES["microgrid"] == Sparlectra.CGMESImporter.CGMES_TESTSET_ALIASES["microgrid_be"]
    end

    # --- real fixtures (optional, offline-safe) -----------------------------
    cache = get(ENV, "SPARLECTRA_CGMES_CACHE", joinpath(dirname(@__DIR__), "data", "CGMES"))
    bc = joinpath(cache, "extracted", "MicroGrid", "BaseCase_BC")
    be = joinpath(bc, "CGMES_v2.4.15_MicroGridTestConfiguration_BC_BE_v2")
    asm = joinpath(bc, "CGMES_v2.4.15_MicroGridTestConfiguration_BC_Assembled_v2")
    bd = joinpath(bc, "CGMES_v2.4.15_MicroGridTestConfiguration_BD_v2")
    if isdir(be) && isdir(bd)
      @testset "MicroGrid BC BE (ENTSO-E fixture)" begin
        s = summarizeCGMES(path = [be, bd])
        @test s.version == "2.4.15"
        @test s.unresolved_count == 0
        @test !s.boundary_missing_hint
        hist = Dict(s.class_histogram)
        @test hist[:TopologicalNode] == 12
        @test hist[:ACLineSegment] == 7
        @test hist[:PowerTransformerEnd] == 9
        # boundary detection: without BD the hint must fire
        s2 = summarizeCGMES(path = be)
        @test s2.unresolved_count > 0
        @test s2.boundary_missing_hint
      end

      # Stage-1 exit criterion: converge and reproduce the SV profile.
      # vm tolerance 2e-4 (the SvVoltage values carry ~6 significant digits,
      # so ~1.5e-4 is the data noise floor), va tolerance 0.05°.
      @testset "Stage 1: import + runpf! + compareWithSV" begin
        res = importCGMES(path = [be, bd], name = "MicroGrid_BE")
        @test length(res.net.nodeVec) == 12
        @test res.slack_bus != ""
        # ReactiveCapabilityCurve (#294 point 1): BE-G1 writes the degenerate
        # 0/0 scalar pair precisely because its real limits live in the curve.
        # Evaluated at its scheduled P = −90 MW the curve gives ±210 MVAr
        # (linear between the −100 → ±200 and 0 → ±300 points), replacing the
        # former ±wide fallback.
        @test any(m -> occursin("BE-G1", m) && occursin("ReactiveCapabilityCurve", m) && occursin("[-210.0, 210.0]", m), res.messages)
        @test any(ps -> Sparlectra.isGenerator(ps) && something(ps.minQ, 0.0) == -210.0 && something(ps.maxQ, 0.0) == 210.0, res.net.prosumpsVec)
        # Distributed slack (#192): MicroGrid SSH carries normalPF only with
        # value 0, which maps to "unknown" — only positive factors count as
        # imported participation.
        @test all(ps.participationFactor === nothing for ps in res.net.prosumpsVec if Sparlectra.isGenerator(ps))
        # short-circuit harvest (§7.7): read, not evaluated
        @test length(res.shortcircuit.ac_line_segments) == 7
        @test length(res.shortcircuit.synchronous_machines) == 2
        @test !isempty(res.shortcircuit.equivalent_injections)
        _, erg = runpf!(res.net, 30, 1e-8, 0)
        @test erg == 0
        cmp = compareWithSV(res)
        @test cmp.n == 11
        @test cmp.max_dvm < 2e-4
        @test cmp.max_dva < 0.05
        # SvPowerFlow comparison (D5b): the conformity sets carry injection
        # terminals only; tolerances reflect the SV data noise floor
        @test cmp.flows.n >= 12
        @test cmp.flows.max_dp < 3.0
        @test cmp.flows.max_dq < 2.0
        @test any(r -> r.kind == :shunt, cmp.flows.rows)
        @test any(r -> r.kind == :units, cmp.flows.rows)

        # boundary error path — the typed abort carries the import analysis
        @test_throws Sparlectra.CGMESImporter.CGMESImportError importCGMES(path = be)
      end

      if isdir(asm)
        @testset "Stage 1: assembled model (both sides + X-node EI skip)" begin
          res = importCGMES(path = [asm, bd], name = "MicroGrid_Assembled")
          @test count(m -> occursin("assembled boundary node", m), res.messages) == 10
          @test count(m -> occursin("bus link", m), res.messages) >= 1   # NL busbar coupler
          _, erg = runpf!(res.net, 30, 1e-8, 0)
          @test erg == 0
          cmp = compareWithSV(res)
          @test cmp.max_dvm < 2e-4
          @test cmp.max_dva < 0.05
        end
      end

      # WebUI text-entry path: alias → combined single ZIP (base case plus
      # boundary) in the case directory; importing that ZIP must reproduce
      # the folder import exactly. Offline-safe: the extracted cache exists
      # here, so no download happens.
      @testset "Test-set fetch: alias to combined ZIP" begin
        out = mktempdir()
        z = Sparlectra.CGMESImporter.fetchCGMESTestSet("microgrid_be"; outdir = out)
        @test isfile(z) && basename(z) == "cgmes_microgrid_be.zip"
        rz = importCGMES(path = z, name = "from_zip")
        @test length(rz.net.nodeVec) == 12
        _, erg = runpf!(rz.net, 30, 1e-8, 0)
        @test erg == 0
        @test compareWithSV(rz).max_dvm < 2e-4
        @test_throws ErrorException Sparlectra.CGMESImporter.fetchCGMESTestSet("no_such_set"; outdir = out)
        # The WebUI runs through the case selector and the service resolver,
        # not through run_sparlectra_api directly — both must accept the ZIP
        # (regression: the delivery imported fine but was invisible in the
        # selector and rejected by the resolver).
        @test Sparlectra._webui_is_user_selectable_case(z)
        @test basename(Sparlectra._resolve_powerflow_casefile(basename(z), out)) == basename(z)
        @test basename(Sparlectra._resolve_powerflow_casefile(z, out)) == basename(z)
        @test_throws ArgumentError Sparlectra._resolve_powerflow_casefile("does_not_exist.zip", out)

        # Upload-time completeness feedback: the user must learn at import
        # time whether a delivery can be computed, not first at run time.
        sgdir = joinpath(cache, "extracted", "SmallGrid", "BusBranch")
        basezip = joinpath(sgdir, "CGMES_v2.4.15_SmallGridTestConfiguration_BaseCase_Complete_v3.0.0.zip")
        bndzip = joinpath(sgdir, "CGMES_v2.4.15_SmallGridTestConfiguration_Boundary_v3.0.0.zip")
        if isfile(basezip) && isfile(bndzip)
          lonely = mktempdir()
          cp(basezip, joinpath(lonely, "base.zip"))
          # With the local test-set cache present the boundary is supplied
          # automatically; the message must say so and the file must appear.
          lonely_role = Sparlectra._webui_cgmes_upload_role(joinpath(lonely, "base.zip"), lonely)
          @test occursin("boundary supplied automatically", lonely_role) || occursin("BOUNDARY SET MISSING", lonely_role)
          if occursin("boundary supplied automatically", lonely_role)
            @test length(filter(n -> occursin(r"(?i)boundary", n), readdir(lonely))) == 1
            rl = run_sparlectra_api(casefile = joinpath(lonely, "base.zip"), output_dir = mktempdir())
            @test rl.status == :succeeded
          end
          # the pre-analysis must state version, profiles and element counts
          @test occursin("CGMES 2.4.15", lonely_role)
          @test occursin("nodes", lonely_role) && occursin("lines", lonely_role) && occursin("transformers", lonely_role)
          together = mktempdir()
          cp(basezip, joinpath(together, "base.zip"))
          cp(bndzip, joinpath(together, "boundary.zip"))
          @test occursin("ready to compute", Sparlectra._webui_cgmes_upload_role(joinpath(together, "base.zip"), together))
          @test occursin("boundary set", Sparlectra._webui_cgmes_upload_role(joinpath(together, "boundary.zip"), together))
          # and the run picks the neighbouring boundary up automatically
          rr = run_sparlectra_api(casefile = joinpath(together, "base.zip"), output_dir = mktempdir())
          @test rr.status == :succeeded
          @test rr.metadata["cgmes_buses"] == 118
        end

        # Boundary-only deliveries are companions, not runnable cases:
        # they must stay out of the case selector.
        bonly = filter(n -> occursin(r"(?i)boundary", n), readdir(out))
        for b in bonly
          @test !Sparlectra._webui_is_user_selectable_case(joinpath(out, b))
        end
        @test Sparlectra._webui_is_user_selectable_case(z)   # the combined set stays selectable

        # WebUI resolve handler on the cgmes: prefix
        resp = Sparlectra.handle_powerflow_case_resolve(Dict("casefile" => "cgmes:microgrid_be"); output_root = mktempdir(), case_directory = out)
        @test resp.status == 303
        headers = resp.headers isa AbstractDict ? resp.headers : Dict(resp.headers)
        @test occursin("cgmes_microgrid_be.zip", get(headers, "Location", ""))
      end

      # D4b: the framework/API path — CGMES delivery as casefile, boundary
      # set from the cgmes_import config block, dispatch by auto-detection
      @testset "D4b: run_sparlectra_api dispatch on a CGMES delivery" begin
        out = mktempdir()
        cfgfile = joinpath(out, "cfg.yaml")
        write(cfgfile, "cgmes_import:\n  path: \"" * bd * "\"\n")
        r = run_sparlectra_api(casefile = be, config_file = cfgfile, output_dir = out)
        @test r.status == :succeeded
        @test r.metadata["input_format_detected"] == "cgmes"
        @test r.metadata["cgmes_version"] == "2.4.15"
        @test r.metadata["cgmes_buses"] == 12
        @test r.metadata["cgmes_slack_bus"] == "BE-Busbar_4"
        @test r.metadata["cgmes_messages"] > 0
        # a dedicated cgmes.log documents the import next to run.log
        logfile = joinpath(out, "cgmes.log")
        @test isfile(logfile)
        text = read(logfile, String)
        @test occursin("# CGMES import report", text)
        @test occursin("## Class histogram", text)
        @test occursin("## Network built", text)
        @test occursin("## Short-circuit source data", text)
        @test occursin("## Importer messages", text)
      end

      # The report must survive a failed solve — it is written right after the
      # import, not after the power flow (regression: a diagnose run whose
      # power flow failed produced no cgmes.log at all).
      @testset "cgmes.log is written even when the run fails" begin
        out2 = mktempdir()
        cfg2 = joinpath(out2, "c.yaml")
        write(cfg2, "power_flow:\n  max_iter: 1\n  tol: 1.0e-14\n")
        rf = run_sparlectra_api(casefile = be, config_file = cfg2, output_dir = out2)
        @test isfile(joinpath(out2, "cgmes.log"))
        @test occursin("# CGMES import report", read(joinpath(out2, "cgmes.log"), String))
      end

      # Short-circuit plausibility on the complete MicroGrid records (machines and
      # feeders fully attributed per the compendium): finite positive Ik'' on
      # every energized bus, max >= min per bus, and no motor-skip flag —
      # consistent with shortCircuitCoverage (BE carries no
      # AsynchronousMachine). The analytic reference and safety-flag tests
      # live in test/test_short_circuit.jl (fast profile).
      @testset "runShortCircuit! MicroGrid plausibility" begin
        res = importCGMES(path = [be, bd], name = "MicroGrid_SC")
        smax = runShortCircuit!(res; case = :max)
        smin = runShortCircuit!(res; case = :min)
        @test length(smax.rows) == 12
        @test all(row.status === :ok for row in smax.rows)
        @test all(isfinite(row.ik_kA) && row.ik_kA > 0.0 for row in smax.rows)
        @test all(isfinite(row.ip_kA) && row.ip_kA >= sqrt(2.0) * row.ik_kA for row in smax.rows)
        mins = Dict(row.bus_idx => row.ik_kA for row in smin.rows)
        @test all(row.ik_kA >= mins[row.bus_idx] for row in smax.rows)
        # no AsynchronousMachine in the delivery -> no motor skip, no flags
        @test !any(ps.comp.cTyp == Sparlectra.AsynchronousMachine for ps in res.net.prosumpsVec)
        @test all(!row.contains_defaulted_data for row in smax.rows)
        @test isempty(smax.messages)
      end

      # MiniGrid carries the full motor attribute set (iaIrRatio,
      # rxLockedRotorRatio, rated data on all six machines) — the motors are
      # therefore computed per IEC 60909-0 §6.7 and the maximum case carries
      # NO lower-bound flags (regression: before the motor formula every
      # MiniGrid island was flagged).
      @testset "runShortCircuit! MiniGrid: motors computed, no flags" begin
        mini = Sparlectra.CGMESImporter.fetchCGMESTestSet("minigrid"; outdir = mktempdir())
        res = importCGMES(path = mini, name = "MiniGrid_SC")
        @test !isempty(res.shortcircuit.asynchronous_machines)
        @test all(m.iaIrRatio !== nothing && m.rxLockedRotorRatio !== nothing for m in res.shortcircuit.asynchronous_machines)
        smax = runShortCircuit!(res; case = :max)
        ok = [row for row in smax.rows if row.status === :ok]
        @test !isempty(ok)
        @test all(isfinite(row.ik_kA) && row.ik_kA > 0.0 for row in ok)
        @test all(!row.contains_defaulted_data for row in ok)
        @test !any(occursin("lower bound", m) for m in smax.messages)
      end

      # The Web UI "Short circuit" button's service path — CGMES
      # import + runShortCircuit! max/min without any power-flow solve,
      # CSV artifacts, registry-compatible result. Negative cases must fail
      # with explicit reasons, never with empty tables.
      @testset "short_circuit_mode service run" begin
        root = mktempdir()
        cfgsc = joinpath(root, "c.yaml")
        write(cfgsc, "power_flow:\n  max_iter: 40\n")
        # the service path accepts case FILES; use the combined alias ZIP
        # (base + boundary in one delivery, packed from the local cache)
        zipcase = Sparlectra.CGMESImporter.fetchCGMESTestSet("microgrid_be"; outdir = mktempdir())
        resp = start_powerflow_run(Dict("casefile" => zipcase, "config_file" => cfgsc, "output_root" => root, "short_circuit_mode" => true))
        @test resp["success"] === true
        @test resp["metadata"]["run_mode"] == "short_circuit"
        rid = resp["run_id"]
        @test isfile(joinpath(root, rid, "short_circuit_max.csv"))
        @test isfile(joinpath(root, rid, "short_circuit_min.csv"))
        max_rows = readlines(joinpath(root, rid, "short_circuit_max.csv"))
        @test max_rows[1] == "bus,vn_kV,island,status,c,zk_ohm,rx_ratio,ik_kA,sk_MVA,kappa,ip_kA,flagged,reasons"
        @test length(max_rows) - 1 == 12
        @test resp["metadata"]["sc_flagged_rows"] == 0
        @test resp["metadata"]["sc_max_ik_kA"] > 0.0
        @test occursin("Short-circuit run", read(joinpath(root, rid, "run.log"), String))
        # run.log carries the coverage report so a flagged run can be audited
        @test occursin("SynchronousMachine", read(joinpath(root, rid, "run.log"), String))

        # MATPOWER case: explicit rejection, no artifacts
        mp = start_powerflow_run(Dict("casefile" => ensure_casefile("case14.m"), "config_file" => cfgsc, "output_root" => root, "short_circuit_mode" => true))
        @test mp["success"] === false
        @test mp["reason"] == "short_circuit_requires_cgmes"

        # diagnose_mode and short_circuit_mode are mutually exclusive
        both = start_powerflow_run(Dict("casefile" => zipcase, "config_file" => cfgsc, "output_root" => root, "short_circuit_mode" => true, "diagnose_mode" => true))
        @test both["success"] === false
        @test occursin("mutually exclusive", both["message"])

        # data-check gating helper recognizes the real delivery
        @test Sparlectra._webui_case_has_short_circuit_data(be) == true
      end

      # The Web UI "Analyze import" button's service path — parse-only
      # pre-check, import_analysis.txt artifact, Failed when the delivery
      # would not import (missing declared dependency).
      @testset "import_analysis_mode service run" begin
        root = mktempdir()
        cfgia = joinpath(root, "c.yaml")
        write(cfgia, "power_flow:\n  max_iter: 40\n")

        # Healthy delivery: succeeded, importable, artifact written.
        okcase = Sparlectra.CGMESImporter.fetchCGMESTestSet("microgrid_be"; outdir = mktempdir())
        okresp = start_powerflow_run(Dict("casefile" => okcase, "config_file" => cfgia, "output_root" => root, "import_analysis_mode" => true))
        @test okresp["success"] === true
        @test okresp["metadata"]["run_mode"] == "import_analysis"
        @test okresp["metadata"]["import_analysis_missing_dependencies"] == 0
        okreport = read(joinpath(root, okresp["run_id"], "import_analysis.txt"), String)
        @test occursin("Supplied models:", okreport)
        @test occursin("Verdict:", okreport)

        # A delivery declaring an absent prerequisite: FAILED with the
        # explicit reason, analysis artifact still written. The service
        # accepts case FILES, so the synthetic delivery is packed as a zip.
        badzip = joinpath(mktempdir(), "bad_delivery.zip")
        ZipArchives.ZipWriter(badzip) do w
          ZipArchives.zip_newfile(w, "bad_EQ.xml")
          write(w, replace(
            _CGMES_SYNTH_EQ,
            "<md:Model.profile>http://entsoe.eu/CIM/EquipmentCore/3/1</md:Model.profile>" =>
              "<md:Model.profile>http://entsoe.eu/CIM/EquipmentCore/3/1</md:Model.profile>\n<md:Model.DependentOn rdf:resource=\"urn:uuid:absent-boundary\"/>",
          ))
        end
        badresp = start_powerflow_run(Dict("casefile" => badzip, "config_file" => cfgia, "output_root" => root, "import_analysis_mode" => true))
        @test badresp["success"] === false
        @test badresp["reason"] == "import_analysis_not_importable"
        @test badresp["metadata"]["import_analysis_missing_dependencies"] == 1
        @test occursin("absent-boundary", read(joinpath(root, badresp["run_id"], "import_analysis.txt"), String))

        # Non-CGMES case: explicit rejection.
        mpia = start_powerflow_run(Dict("casefile" => ensure_casefile("case14.m"), "config_file" => cfgia, "output_root" => root, "import_analysis_mode" => true))
        @test mpia["success"] === false
        @test mpia["reason"] == "import_analysis_requires_cgmes"

        # Mode exclusivity.
        exia = start_powerflow_run(Dict("casefile" => okcase, "config_file" => cfgia, "output_root" => root, "import_analysis_mode" => true, "short_circuit_mode" => true))
        @test exia["success"] === false
        @test occursin("excludes", exia["message"])
      end

      # cgmes_import.start_values (WebUI-selectable start state) + the
      # mandatory SV comparison artifact. Effectiveness is asserted via the
      # logged decision line (the contract), not by re-deriving solver
      # internals; precedence over a hostile power_flow.flatstart is part of
      # the same assertions.
      @testset "start_values selection + mandatory SV comparison" begin
        for (mode, hostile_flatstart) in (("flat", "false"), ("sv", "true"))
          out = mktempdir()
          cfg = joinpath(out, "c.yaml")
          # hostile flatstart value opposes what the mode needs — start_values must win
          write(cfg, "power_flow:\n  flatstart: " * hostile_flatstart * "\ncgmes_import:\n  start_values: " * mode * "\n  path: \"" * bd * "\"\n")
          r = run_sparlectra_api(casefile = be, config_file = cfg, output_dir = out, case_format = :cgmes)
          @test r.status == :succeeded
          run_log = read(joinpath(out, "run.log"), String)
          cgmes_log = read(joinpath(out, "cgmes.log"), String)
          @test occursin("CGMES start values: " * mode, run_log)
          @test occursin("CGMES start values: " * mode, cgmes_log)
          @test occursin("overrides: power_flow.flatstart=" * hostile_flatstart, run_log)
          mode == "sv" && @test occursin("start-value machines forced off", run_log)
          # SV comparison runs for both modes: artifacts, log block, metadata
          @test isfile(joinpath(out, "sv_compare.csv"))
          @test isfile(joinpath(out, "sv_compare_flows.csv"))
          @test occursin("## SV comparison (run status: converged, start values: " * mode * ")", cgmes_log)
          sv_rows = readlines(joinpath(out, "sv_compare.csv"))
          @test sv_rows[1] == "bus,vm_pu,sv_vm_pu,dvm,va_deg,sv_va_deg,dva,dva_aligned"
          @test length(sv_rows) - 1 == 11
          @test r.metadata["cgmes_start_values"] == mode
          @test r.metadata["cgmes_sv_compare_status"] == "converged"
          @test r.metadata["cgmes_sv_compare_n"] == 11
          for key in ("cgmes_sv_compare_max_dvm", "cgmes_sv_compare_rms_dvm", "cgmes_sv_compare_max_dva", "cgmes_sv_compare_rms_dva", "cgmes_sv_compare_va_ref_offset_deg")
            @test isfinite(r.metadata[key])
          end
          # solved-from-SV and solved-from-flat both land on the SV solution
          @test r.metadata["cgmes_sv_compare_max_dvm"] < 2e-4
        end

        # auto (the shipped default) resolves per delivery: this one carries an
        # SvVoltage state, so it must solve from it and say so.
        auto_out = mktempdir()
        auto_cfg = joinpath(auto_out, "c.yaml")
        write(auto_cfg, "cgmes_import:\n  start_values: auto\n  path: \"" * bd * "\"\n")
        auto_r = run_sparlectra_api(casefile = be, config_file = auto_cfg, output_dir = auto_out, case_format = :cgmes)
        @test auto_r.status == :succeeded
        auto_log = read(joinpath(auto_out, "run.log"), String)
        @test occursin("CGMES start values: sv (auto: delivery carries SvVoltage for", auto_log)
        @test occursin("start-value machines forced off", auto_log)
        @test auto_r.metadata["cgmes_start_values"] == "sv"

        # Negative: a MATPOWER run ignores the key completely — no decision
        # line, no artifacts, no metadata keys.
        mp_out = mktempdir()
        mp_cfg = joinpath(mp_out, "c.yaml")
        write(mp_cfg, "cgmes_import:\n  start_values: sv\n")
        mp = run_sparlectra_api(casefile = ensure_casefile("case14.m"), config_file = mp_cfg, output_dir = mp_out)
        @test mp.status == :succeeded
        @test !isfile(joinpath(mp_out, "sv_compare.csv"))
        @test !haskey(mp.metadata, "cgmes_start_values")
        @test !haskey(mp.metadata, "cgmes_sv_compare_status")
        @test !occursin("CGMES start values", read(joinpath(mp_out, "run.log"), String))
      end

      # Fixed-reference self-check on a CGMES delivery: the SV voltages must
      # reach the solver verbatim (a base config with flatstart: true —
      # observed in a real WebUI configuration — must not wipe them), and the
      # residual/attribution artifacts must be written. MicroGrid's SV state
      # carries a known residual of ~2.6 pu, i.e. clearly not the flat-start
      # value, which pins the start point actually used.
      @testset "CGMES fixed-reference self-check (SV start, artifacts)" begin
        out3 = mktempdir()
        cfg3 = joinpath(out3, "c.yaml")
        write(cfg3, "power_flow:\n  flatstart: true\ncgmes_import:\n  path: \"" * bd * "\"\n")
        rsc = run_fixed_reference_self_check(casefile = be, config_file = cfg3, output_dir = out3, case_format = :cgmes)
        @test rsc.raw_result !== nothing
        @test rsc.raw_result.iterations == 1
        @test rsc.metadata["cgmes_no_sv_buses"] == 0
        summary = read(joinpath(out3, "self_check.log"), String)
        @test occursin("start values taken verbatim from import", summary)
        m = match(r"start_state_residual_inf: ([0-9.eE+-]+)", summary)
        @test m !== nothing
        # SV-state residual of this fixture; a flat start would report ~4.3.
        @test isapprox(parse(Float64, m.captures[1]), 2.6205572054022745; rtol = 1e-6)
        residuals = readlines(joinpath(out3, "self_check_residuals.csv"))
        @test residuals[1] == "bus_id,bus_name,vn_kV,bus_type,vm_pu_start,va_deg_start,p_residual,q_residual,has_sv,n_transformer_terminals,n_shunts"
        @test length(residuals) - 1 == 12
        @test all(split(l, ',')[9] == "true" for l in residuals[2:end])
        # transformer adjacency must be populated (B_2WT_ markers)
        @test any(parse(Int, split(l, ',')[10]) > 0 for l in residuals[2:end])
        @test occursin("no-SV buses: 0", read(joinpath(out3, "cgmes.log"), String))
      end

      # Stage 2 (Phase E): start from the SSH tap positions and let the
      # outer-loop controllers find the operating point. Note that the
      # reference SvTapStep positions are NOT reproduced exactly — the CGMES
      # target deadbands of this data set are wide (PST: 35 MW, OLTC: 0.5 kV
      # on 10.5 kV), so a controller legitimately stops one step short. The
      # testable properties are: controllers are created, the loop converges,
      # targets are met within their deadband, and the tap moves toward the
      # SV reference rather than away from it.
      @testset "Stage 2: tap controllers from TapChangerControl" begin
        stage1 = importCGMES(path = [be, bd], name = "s1")
        @test isempty(collect(Sparlectra._tap_controllers(stage1.net)))

        res = importCGMES(path = [be, bd], tap_control = true, name = "s2")
        ctrls = collect(Sparlectra._tap_controllers(res.net))
        @test length(ctrls) == 2
        @test any(c -> c.mode == :branch_active_power && c.control_phase, ctrls)
        @test any(c -> c.mode == :voltage && c.control_ratio, ctrls)

        # the OLTC of this set regulates the slack busbar → disabled by the guard
        @test count(m -> occursin("disabled — target bus is", m), res.messages) == 1
        @test count(c -> !c.enabled, ctrls) == 1

        cfg = SparlectraConfig(powerflow = PowerFlowConfig(max_iter = 30, tol = 1e-8), output = OutputConfig(logfile_results = :off), control = ControlConfig())
        runres = run_sparlectra(net = res.net, config = cfg)
        @test runres.numerical_converged
        cres = latest_control_result(res.net)
        @test cres.status == :converged

        pst = only(filter(c -> c.mode == :branch_active_power, ctrls))
        @test pst.converged
        @test abs(pst.achieved_p_mw - pst.p_target_mw) <= pst.deadband_p_mw

        # Taps are discrete, so the distance to the reference is measured in
        # STEPS, not in degrees: project the controlled angle back onto the
        # changer's step grid (same evaluation the importer uses, so the
        # asymmetrical PST's non-uniform degree spacing cannot distort the
        # result) and compare integers. The solver moves in a constant
        # `phase_step_deg`, i.e. a linear approximation of that grid — if the
        # approximation ever became too coarse, this assertion is where it
        # would surface.
        pstbr = Sparlectra._find_trafo_branch(res.net, pst.trafo)
        ptc = only(objectsOf(res.store, :PhaseTapChangerAsymmetrical))
        low = Int(round(num(ptc, :lowStep, 0.0)))
        high = Int(round(num(ptc, :highStep, 0.0)))
        nearest_step(angle) = argmin(
          [abs(Sparlectra.CGMESImporter._phaseTapRatioShift_atstep(res.store, ptc, s)[2] - angle) for s = low:high],
        ) + low - 1
        step_controlled = nearest_step(pstbr.phase_shift_deg)
        step_sv = 16         # SvTapStep position of BE-TR2_1
        @test abs(step_controlled - step_sv) <= 1
      end
    else
      @info "CGMES MicroGrid fixture not cached — skipping ENTSO-E fixture tests (run examples/experimental/cgmes_fetch_testsets.jl to enable)"
    end

    @testset "3W-leg tap controller from TapChangerControl (#294 point 4)" begin
      dir = _cgmes_synth3w_dir()
      # without tap_control: no controllers, taps stay at their SSH position
      plain = importCGMES(path = dir, name = "synth3w_plain")
      @test isempty(collect(Sparlectra._tap_controllers(plain.net)))
      # with tap_control: the star-equivalent leg gets the controller
      res = importCGMES(path = dir, name = "synth3w", tap_control = true)
      @test !any(m -> occursin("not wired in Stage 2", m), res.messages)
      @test any(m -> occursin("tap control: T3W (leg 2)", m), res.messages)
      ctrls = collect(Sparlectra._tap_controllers(res.net))
      @test length(ctrls) == 1
      oltc = only(ctrls)
      @test oltc.mode == :voltage && oltc.control_ratio && oltc.enabled
      cfg = SparlectraConfig(powerflow = PowerFlowConfig(max_iter = 30, tol = 1e-8), output = OutputConfig(logfile_results = :off), control = ControlConfig())
      runres = run_sparlectra(net = res.net, config = cfg)
      @test runres.numerical_converged
      @test latest_control_result(res.net).status == :converged
      @test oltc.converged
      mv = res.net.nodeVec[res.net.busDict["MV_BUS"]]
      # target 1.04 pu ± 0.0118 pu — and the tap really moved: the
      # uncontrolled operating point sits at ≈1.0135 pu, outside the band
      @test abs(mv._vm_pu - 1.04) <= 2.6 / 110.0 / 2.0 + 1e-9
      @test mv._vm_pu > 1.025
    end

    @testset "Remote-regulating machine (#294 point 3)" begin
      # default: the remote control stays the Stage-1 held-PV fallback
      plain = importCGMES(path = _cgmes_synth_rvc_dir(), name = "synthrvc_plain")
      @test any(m -> occursin("G_B has a remote voltage RegulatingControl", m) && occursin("held PV at its own bus", m), plain.messages)
      @test isempty(Sparlectra._machine_controllers(plain.net))
      b_plain = plain.net.nodeVec[plain.net.busDict["BUS_B"]]
      @test b_plain._nodeType == Sparlectra.PV

      # machine_control = true: G_B becomes a PQ injection with an outer-loop
      # controller regulating BUS_C to 112.2 kV = 1.02 pu
      res = importCGMES(path = _cgmes_synth_rvc_dir(), name = "synthrvc", machine_control = true)
      @test any(m -> occursin("machine control: G_B", m) && occursin("1.02 pu at BUS_C", m), res.messages)
      ctrls = Sparlectra._machine_controllers(res.net)
      @test length(ctrls) == 1
      ctrl = only(ctrls)
      @test ctrl.target_bus == "BUS_C"
      @test isapprox(ctrl.target_vm_pu, 112.2 / 110.0; atol = 1e-9)
      @test ctrl.qmin_mvar == -50.0 && ctrl.qmax_mvar == 50.0
      @test res.net.nodeVec[res.net.busDict["BUS_B"]]._nodeType == Sparlectra.PQ
      # the control loop drives the remote bus onto the target
      cfg = SparlectraConfig(powerflow = PowerFlowConfig(max_iter = 30, tol = 1e-8), output = OutputConfig(logfile_results = :off), control = ControlConfig())
      runres = run_sparlectra(net = res.net, config = cfg)
      @test runres.numerical_converged
      @test latest_control_result(res.net).status == :converged
      @test ctrl.converged
      @test abs(res.net.nodeVec[res.net.busDict["BUS_C"]]._vm_pu - 112.2 / 110.0) <= ctrl.deadband_vm_pu + 1e-9

      # fallback: with G_C reconnected the target bus is already voltage-held
      # locally — the plan reverts to held-PV with a notice
      held = importCGMES(path = _cgmes_synth_rvc_dir(sm2_connected = true), name = "synthrvc_held", machine_control = true)
      @test any(m -> occursin("G_B", m) && occursin("remote voltage control not attachable", m) && occursin("already voltage-held", m), held.messages)
      @test isempty(Sparlectra._machine_controllers(held.net))
      @test held.net.nodeVec[held.net.busDict["BUS_B"]]._nodeType == Sparlectra.PV
    end

    @testset "NonlinearShuntCompensator mapping (#294 point 7)" begin
      dir = _cgmes_synth3w_dir()
      res = importCGMES(path = dir, name = "synth3w_nlsh")
      # sections = 2 of 3 points: B = (1.0 + 0.8)e-4 S * 110² = 2.178 MVAr
      @test any(m -> occursin("nonlinear shunt NLSH", m) && occursin("sections 2/3", m) && occursin("2.178 MVAr", m), res.messages)
      @test !any(m -> occursin("skip: NonlinearShuntCompensator", m), res.messages)
      # multi-valued reference guard (#294 point 9): the SV TopologicalIsland
      # membership list triggers exactly one class/property notice
      @test count(m -> occursin("multi-valued reference TopologicalIsland.TopologicalNodes on 1 object(s)", m), res.messages) == 1
      @test length(res.net.shuntVec) == 1
      _, erg = runpf!(res.net, 30, 1e-8, 0)
      @test erg == 0
    end

    # PST tabular wiring + angle-sign convention (#294 point 2). Sparlectra
    # folds an end-2 (to-side) tap angle NEGATED into the branch shift
    # (θ_eff = θ1 − θ2, standard end-referral; deciding fixture: RealGrid,
    # asserted below). KNOWN AMBIGUITY: the ENTSO-E PSEI PTE2 sets encode the
    # same parameters on end 2 but their SV expects the UNflipped angle —
    # those two sets therefore deviate by dva ≈ 0.29° / dp ≈ 0.79 MW by
    # design of the chosen convention. PTE1 (end 1) is convention-independent
    # and must match SV exactly.
    pst_sets = ("PST_PTChLin_PTE1_PSEI", "PST_PTChLin_PTE2_PSEI", "PST_PTChTab_PTE2_PSEI")
    if all(isdir(joinpath(cache, "extracted", s)) for s in pst_sets)
      @testset "PST tabular wiring and tap-side sign (PSEI sets)" begin
        res1 = importCGMES(path = joinpath(cache, "extracted", "PST_PTChLin_PTE1_PSEI"), name = "pte1")
        _, erg1 = runpf!(res1.net, 30, 1e-10, 0)
        @test erg1 == 0
        @test [b.phase_shift_deg for b in res1.net.branchVec if b.phase_shift_deg != 0.0] == [-5.0]
        cmp1 = compareWithSV(res1)
        @test cmp1.max_dva < 1e-4
        @test cmp1.flows.max_dp < 1e-3
        for set in ("PST_PTChLin_PTE2_PSEI", "PST_PTChTab_PTE2_PSEI")
          res2 = importCGMES(path = joinpath(cache, "extracted", set), name = set)
          _, erg2 = runpf!(res2.net, 30, 1e-10, 0)
          @test erg2 == 0
          # end-2 flip: table/formula angle −5 arrives as +5 on the branch
          @test [b.phase_shift_deg for b in res2.net.branchVec if b.phase_shift_deg != 0.0] == [5.0]
          # documented deviation bound of the convention choice — NOT noise;
          # shrinking this bound means the convention handling changed.
          # Measured on the RAW deltas: the reference alignment (median
          # offset removal) would absorb part of this systematic deviation.
          cmp2 = compareWithSV(res2)
          @test 0.2 < maximum(abs(r.dva) for r in cmp2.rows) < 0.4
        end
        # tabular-specific: the mapping message reports the resolved table row
        res_tab = importCGMES(path = joinpath(cache, "extracted", "PST_PTChTab_PTE2_PSEI"), name = "ptchtab")
        @test any(m -> occursin("tabular phase tap", m) && occursin("step 2", m) && occursin("-5.0", m), res_tab.messages)
        @test !any(m -> occursin("skip: PhaseTapChangerTabular", m), res_tab.messages)
      end
    else
      @info "CGMES PSEI PST fixtures not cached — PST sign/tabular tests skipped (fetch via examples/experimental/cgmes_fetch_testsets.jl)"
    end

    # RealGrid is the deciding data-faithfulness fixture for both the tabular
    # wiring and the end-2 angle flip: its nine PhaseTapChangerTabular
    # changers (four with real flow) reproduce the delivered SV state only
    # with table+flip applied — SV-start max|F| drops from ≈197 pu (taps
    # skipped) / ≈395 pu (taps unflipped) to ≈5 pu. The bound 10 pu holds
    # that ordering with margin.
    rgdir = joinpath(cache, "extracted", "RealGrid", "CGMES_v2.4.15_RealGridTestConfiguration_v2")
    if isdir(rgdir)
      @testset "RealGrid tabular PSTs reproduce SV (tap-side sign)" begin
        res = importCGMES(path = rgdir, name = "RealGrid")
        @test count(m -> occursin("tabular phase tap", m), res.messages) == 9
        net = res.net
        n = length(net.nodeVec)
        Yred = Sparlectra.createYBUS(net = net, sparse = true, printYBUS = false)
        Yb = size(Yred, 1) == n ? Yred : Sparlectra._expand_ybus_for_isolated_nodes(Yred, n, net.isoNodes)
        V0, slack = Sparlectra.initialVrect(net; flatstart = false)
        S = Sparlectra.buildComplexSVec(net)
        bt = [Sparlectra.getNodeType(nd) == Sparlectra.Slack ? :Slack : Sparlectra.getNodeType(nd) == Sparlectra.PV ? :PV : :PQ for nd in net.nodeVec]
        Vset = [something(nd._vm_pu, 1.0) for nd in net.nodeVec]
        F = Sparlectra.mismatch_rectangular(Yb, V0, S, bt, Vset, slack)
        @test maximum(abs.(F)) < 10.0

        # #294 point 10: the SvPowerFlow comparison must stay available in the
        # presence of isolated buses (RealGrid ships 158 de-energized ones) —
        # the reduced Ybus is re-embedded and iso rows are filtered instead of
        # switching the whole flow comparison off. Bounds are generous
        # envelopes over the measured values (n=6861, rms dp 0.026 MW,
        # max dp 2.1 MW, max dvm 0.031, max dva 0.96°).
        @test !isempty(net.isoNodes)
        _, erg = runpf!(net, 60, 1e-8, 0)
        @test erg == 0
        cmp = compareWithSV(res)
        @test cmp.n > 5000
        @test cmp.max_dvm < 0.05
        @test cmp.max_dva < 2.0
        @test cmp.flows.n > 5000
        @test cmp.flows.rms_dp < 1.0
        @test cmp.flows.max_dp < 20.0
      end
    else
      @info "CGMES RealGrid fixture not cached — tabular-PST SV regression skipped (fetch via examples/experimental/cgmes_fetch_testsets.jl)"
    end

    # Distributed slack (#192): positive GeneratingUnit.normalPF must arrive
    # as ProSumer.participationFactor through the full import path. MiniGrid
    # BusBranch is the smallest cached set with a positive factor (one unit
    # carries normalPF = 1; the zero-valued rest maps to unknown).
    mini = joinpath(cache, "extracted", "MiniGrid", "BusBranch")
    mini_bc = joinpath(mini, "CGMES_v2.4.15_MiniGridTestConfiguration_BaseCase_v3")
    mini_bd = joinpath(mini, "CGMES_v2.4.15_MiniGridTestConfiguration_Boundary_v3")
    if isdir(mini_bc) && isdir(mini_bd)
      @testset "normalPF arrives as participationFactor (MiniGrid)" begin
        res = importCGMES(path = [mini_bc, mini_bd], name = "MiniGrid_BC")
        pfs = [ps.participationFactor for ps in res.net.prosumpsVec if Sparlectra.isGenerator(ps) && ps.participationFactor !== nothing]
        @test !isempty(pfs)
        @test all(==(1.0), pfs)
      end
    else
      @info "CGMES MiniGrid fixture not cached — normalPF arrival test skipped (fetch via examples/experimental/cgmes_fetch_testsets.jl)"
    end

    # FullGrid is ENTSO-E's import/export completeness set: it extends
    # MicroGrid T3 so that every class defined in the CGMES profiles appears
    # at least once (HVDC converters, SVC, AsynchronousMachine, all four PST
    # types). It is explicitly NOT an analytical fixture — no SV accuracy is
    # asserted here. Its job is to prove that the reader walks past unknown
    # classes instead of stumbling, and that they are counted, not dropped.
    fg = joinpath(cache, "extracted", "FullGrid")
    fgbb = joinpath(fg, "CGMES_v2.4.15_FullGridTestConfiguration_BB_BE_v2")
    fgnb = joinpath(fg, "CGMES_v2.4.15_FullGridTestConfiguration_NB_BE_v4")
    fgbd = joinpath(fg, "CGMES_v2.4.15_FullGridTestConfiguration_BD_v1")
    if isdir(fgbb) && isdir(fgnb) && isdir(fgbd)
      @testset "FullGrid: reader completeness across all CGMES classes" begin
        s = summarizeCGMES(path = [fgbb, fgbd])
        @test s.version == "2.4.15"
        @test length(s.class_histogram) > 100      # every profile class present
        @test !s.boundary_missing_hint
        hist = Dict(s.class_histogram)
        # classes far outside the Stage-1 mapping set must still be counted
        for cls in (:VsConverter, :DCLineSegment, :StaticVarCompensator, :AsynchronousMachine, :PhaseTapChangerTabular)
          @test haskey(hist, cls)
        end
        # boundary detection still fires on the richer set
        s2 = summarizeCGMES(path = fgnb)
        @test s2.boundary_missing_hint

        # #297 Draft B: the DC topology (ACDCConverterDCTerminal, DCNode,
        # DCLineSegment) groups the four converters into two links; both are
        # named in the default-mode messages and attached as controllers
        # under paired_control (FullGrid SSH has a zero side, so the derived
        # loss equals the transfer; the notice states the numbers honestly).
        store_fg = Sparlectra.CGMESImporter.loadCGMES([fgbb, fgbd])
        fg_pairs = Sparlectra.CGMESImporter._detectHvdcPairs(store_fg)
        @test count(p -> p.extra == 0, fg_pairs) == 2
        rp = importCGMES(path = [fgbb, fgbd], name = "FullGrid_paired", hvdc_mode = :paired_control)
        @test length(Sparlectra._hvdc_pair_controllers(rp.net)) == 2
        @test count(m -> occursin("HVDC pair attached as controller", m), rp.messages) == 2

        # Node-breaker sets ship the TP profile too, so they read as
        # bus-branch today: the NB variant must yield the same network as
        # the BB variant of the same grid (no node-breaker stage needed).
        rbb = importCGMES(path = [fgbb, fgbd], name = "FullGrid_BB")
        rnb = importCGMES(path = [fgnb, fgbd], name = "FullGrid_NB")
        @test length(rbb.net.nodeVec) == length(rnb.net.nodeVec)
        @test length(rbb.net.branchVec) == length(rnb.net.branchVec)
        @test length(rbb.net.linkVec) == length(rnb.net.linkVec)

        # StaticVarCompensator: mapped as a P = 0 reactive injection
        # (MATPOWER-style); FullGrid's placeholder Ω ratings (±0.99 Ω at
        # 225 kV ↔ ±51 GVAr) trigger the plausibility clamp, and the droop
        # slope is reported as ignored
        @test any(m -> occursin("StaticVarCompensator", m) && occursin("clamped", m), rbb.messages)
        @test any(m -> occursin("StaticVarCompensator", m) && occursin("slope", m), rbb.messages)
        @test any(ps -> something(ps.pVal, -1.0) == 0.0 && something(ps.minQ, 0.0) == -1000.0 && something(ps.maxQ, 0.0) == 1000.0, rbb.net.prosumpsVec)
        @test string(Sparlectra.toProSumptionType("STATICVARCOMPENSATOR")) == "Injection"

        # D6: FullGrid populates the ExternalNetworkInjection short-circuit
        # attributes that MicroGrid lacks — harvested, not evaluated
        eni = only(rbb.shortcircuit.external_network_injections)
        @test eni.maxInitialSymShCCurrent_A !== nothing
        @test eni.minInitialSymShCCurrent_A !== nothing
        @test eni.maxR1ToX1Ratio !== nothing
        @test eni.maxZ0ToZ1Ratio !== nothing
        @test eni.bus !== nothing

        # completeness view of the harvested data (#277 input): per class the
        # record count and per attribute the fill rate
        cov = shortCircuitCoverage(rbb.shortcircuit)
        # six classes since the AsynchronousMachine motor attributes joined
        # the harvest (locked-rotor data for the short-circuit evaluation)
        @test length(cov) == 6
        @test any(r -> r.class == "AsynchronousMachine", cov)
        eni_cov = only(filter(r -> r.class == "ExternalNetworkInjection", cov))
        @test eni_cov.records == 1
        @test any(a -> a.attribute == :maxInitialSymShCCurrent_A && a.filled == a.total == 1, eni_cov.attributes)
        buf = IOBuffer()
        printShortCircuitCoverage(buf, rbb.shortcircuit)
        rendered = String(take!(buf))
        @test occursin("SynchronousMachine", rendered)
        @test occursin("✓", rendered)

        # Placeholder guards (2026-07-30): FullGrid's systematic X.99 family —
        # the tabular PST table row ratio 9.99 and a nonlinear shunt point
        # with b = g = 0.99 S (a 50-GW shunt at 225 kV) — is detected, warned
        # about, and kept out of the model; the AsynchronousMachine maps as a
        # Stage-0 PQ operating point (#294 point 6). With that the network
        # itself SOLVES from a flat start. The shipped SV stays internally
        # inconsistent (14.5° across a 0.3 Ω line), so the SV-based start and
        # the SV comparison remain meaningless for this set — deliberately
        # not asserted here.
        @test any(m -> startswith(m, "warning:") && occursin("implausible tap correction", m), rbb.messages)
        @test any(m -> startswith(m, "warning:") && occursin("implausible admittance", m) && occursin("shunt skipped", m), rbb.messages)
        @test any(m -> occursin("AsynchronousMachine ASM_1", m) && occursin("fixed PQ operating point", m), rbb.messages)
        _, ergflat = runpf!(rbb.net, 60, 1e-8, 0; opt_flatstart = true)
        @test ergflat == 0

        # cgmes_import.placeholder_guards = strict: the same filler values
        # abort the import instead of being skipped
        strict_err = try
          importCGMES(path = [fgbb, fgbd], name = "FullGrid_strict", strict_placeholder_guards = true)
          nothing
        catch e
          e
        end
        @test strict_err isa ErrorException
        @test occursin("placeholder_guards = strict", sprint(showerror, strict_err))
      end
    end

    mgbase = joinpath(cache, "extracted", "MiniGrid", "BusBranch")
    mgbc2 = joinpath(mgbase, "CGMES_v2.4.15_MiniGridTestConfiguration_BaseCase_v3")
    mgbd2 = joinpath(mgbase, "CGMES_v2.4.15_MiniGridTestConfiguration_Boundary_v3")
    if isdir(mgbc2) && isdir(mgbd2)
      @testset "MiniGrid: AsynchronousMachine mapping closes the SV gap (#294 point 6)" begin
        res = importCGMES(path = [mgbc2, mgbd2], name = "MiniGrid_BB")
        asm = filter(m -> occursin("AsynchronousMachine", m) && occursin("fixed PQ operating point", m), res.messages)
        @test length(asm) == 3          # M3, M2a, M2b — 9 MW / 5 MVAr of motor load at bus 7
        _, erg = runpf!(res.net, 30, 1e-8, 0)
        @test erg == 0
        cmp = compareWithSV(res)
        # without the motors the solve missed SV by max dvm 0.031 / flow rms
        # 5.24 MW; with them MiniGrid is SV-tight
        @test cmp.max_dvm < 1e-3
        @test cmp.flows.rms_dp < 0.05
      end
    else
      @info "CGMES MiniGrid fixture not cached — AsynchronousMachine SV test skipped (fetch via examples/experimental/cgmes_fetch_testsets.jl)"
    end

    sgbase = joinpath(cache, "extracted", "SmallGrid", "BusBranch")
    sgbc = joinpath(sgbase, "CGMES_v2.4.15_SmallGridTestConfiguration_BaseCase_Complete_v3.0.0")
    sgbd = joinpath(sgbase, "CGMES_v2.4.15_SmallGridTestConfiguration_Boundary_v3.0.0")
    if isdir(sgbc) && isdir(sgbd)
      @testset "Stage 1: SmallGrid bus-branch" begin
        res = importCGMES(path = [sgbc, sgbd], name = "SmallGrid_BB")
        @test length(res.net.nodeVec) == 118
        _, erg = runpf!(res.net, 30, 1e-8, 0)
        @test erg == 0
        cmp = compareWithSV(res)
        @test cmp.max_dvm < 2e-4
        @test cmp.max_dva < 0.05
        @test cmp.flows.n >= 100
        @test cmp.flows.max_dp < 0.5
        @test cmp.flows.max_dq < 0.5
      end
    end

    # --- ReliCapGrid / Svedala: second source, first CGMES 3.0 delivery -----
    #
    # Svedala comes from ENTSO-E's ReliCapGrid repository, not the conformity
    # package, and its layout differs fundamentally: individual CIM/XML files
    # per model, one boundary file per border, plus a shared commonData file
    # (which alone takes the unresolved references from 571 down to 20).
    #
    # Its job in this suite is threefold:
    #   1. it is the only CGMES 3.0 delivery we have — the 3.0 path was
    #      prepared in the schema layer but never exercised by real data;
    #   2. it pins the SSH `Equipment.inService` semantics (a CGMES 3.0
    #      addition): ReliCapGrid parks out-of-service units and switches with
    #      connected=true, inService=false, no SvVoltage — ignoring the flag
    #      imported 413 MW of phantom generation and falsely closed switches
    #      between separately-solved network parts;
    #   3. its SV declares two angle references in two TopologicalIslands, and
    #      with inService honored our electrical-island view AGREES with that
    #      partition — one reference per island under multi_slack.
    #
    # Historical note: this fixture long asserted NO converged power flow
    # (border nodes hanging free). The accumulated fixes changed that — with
    # Equipment.inService honored, per-island references and the single-sided
    # border equivalents standing in for the absent neighbours, Svedala alone
    # now solves AND reproduces its SV state; the testset asserts it.
    # The per-model cache holds only the grid files; the border files and the
    # commonData catalog are cached once in _shared (see _rcgMemberFiles). A
    # raw model folder is therefore not importable on its own — resolve the
    # delivery through the alias mechanism, which packs grid + borders +
    # commonData into one ZIP and works without network when the cache is warm.
    rcg_root = joinpath(cache, "relicapgrid", "cgmes-3.0_ncp-2.5_tc-2.0")
    if isdir(joinpath(rcg_root, "Svedala")) && isdir(joinpath(rcg_root, "_shared"))
      svedala = Sparlectra.CGMESImporter.fetchCGMESTestSet("svedala"; outdir = mktempdir())
      @testset "ReliCapGrid Svedala (CGMES 3.0)" begin
        s = summarizeCGMES(path = svedala)
        @test s.version == "3.0"
        @test !s.boundary_missing_hint            # the alias ZIP carries the border files
        @test sum(v for (_, v) in s.class_histogram) == 14506
        @test length(s.class_histogram) == 52
        # 20 references stay unresolved even with commonData: they point into
        # the neighbouring areas that this single-area delivery does not carry.
        @test s.unresolved_count == 20

        res = importCGMES(path = svedala, name = "Svedala")
        net = res.net
        # buses are created lazily, so equipment skipped as out of service
        # (32 objects here) never materializes its bus: 132, not the 200 of the
        # raw TP profile. Only 1 of the 47 retained switches is actually in
        # service AND closed.
        @test length(net.nodeVec) == 132
        @test length(net.branchVec) == 150
        @test length(net.linkVec) == 1
        @test count(m -> occursin("out of service (SSH inService=false)", m), res.messages) == 32
        @test any(m -> occursin("VATTENDRAGET_G1", m) && occursin("out of service", m), res.messages)

        # The SV profile names the reference the exporting tool actually used;
        # that outranks referencePriority in the slack chain.
        @test any(m -> occursin("SV declares", m) && occursin("BP_SD-BO3", m) && occursin("HÄLLAN_G1_CNODE", m), res.messages)
        @test any(m -> occursin("slack: taken from the SV angle reference", m), res.messages)

        # The two angle references sit in two different electrical islands —
        # matching the delivery's own two TopologicalIslands. (With inService
        # ignored, falsely-closed switches merged the islands and made the
        # second reference look physically wrong; that was an artifact.)
        bp = net.busDict["BP_SD-BO3"]
        hg = net.busDict["HÄLLAN_G1_CNODE"]
        island_of = Sparlectra.electricalIslandOfBus(net)
        @test island_of[bp] != island_of[hg]
        @test length(Sparlectra.electricalIslandComponents(net)) == 2
        # branch-only view still differs from the electrical view: the one
        # closed link joins two branch-components into one electrical island
        @test length(Sparlectra._active_ac_island_components(net)) == 3

        # multi_slack is the default: each electrical island gets its own
        # SV-declared reference — matching the delivery's TopologicalIslands,
        # and required for the island-wise power flow to run at all (an island
        # without any reference fails _validate_island_references!).
        @test net.slackVec == [bp, hg]
        @test Sparlectra.getNodeType(net.nodeVec[bp]) == Sparlectra.Slack
        @test any(m -> occursin("additional reference", m) && occursin("BP_SD-BO3", m), res.messages)
        # the legacy single-reference behavior stays available
        res_ss = importCGMES(path = svedala, name = "Svedala_single_slack", multi_slack = false)
        @test res_ss.net.slackVec == [hg]
        @test Sparlectra.getNodeType(res_ss.net.nodeVec[res_ss.net.busDict["BP_SD-BO3"]]) != Sparlectra.Slack
        @test any(m -> occursin("BP_SD-BO3", m) && occursin("not promoted to slack", m), res_ss.messages)

        # The single area solves and reproduces its own SV state: the kept
        # single-sided border equivalents stand in for the absent neighbours.
        # (Controlled solver settings — the Q-limit path has its own open
        # solver issues and is not what this fixture is accountable for.)
        ite_sv, erg_sv = runpf!(net, 100, 1e-6, 0; method = :rectangular, islands_enabled = true,
                                islands_diagnostic_continue_after_failure = true,
                                qlimits_enabled = false, opt_flatstart = true)
        @test erg_sv == 0
        cmp_sv = compareWithSV(res)
        @test cmp_sv.n > 100
        @test cmp_sv.rms_dva < 5.0     # measured 2.05°
        @test cmp_sv.rms_dvm < 0.08    # measured 0.037
      end

      # Svedala's five parked generators also carry placeholder voltage targets
      # (RegulatingControl.targetValue = 0.001 kV ≈ 5e-5 pu). With inService
      # honored they are skipped outright — but the plausibility band must
      # still catch such values when the units are forced alive, otherwise a
      # delivery that parks units some OTHER way produces a reference bus at
      # ~0 V and a power flow that "converges" onto an all-zero solution.
      @testset "implausible voltage setpoints are rejected, not obeyed" begin
        # default import: parked units are skipped before their setpoint is read
        res = importCGMES(path = svedala, name = "Svedala_vset")
        @test isempty([m for m in res.messages if startswith(m, "warning:")])

        # ignore_connected revives them — now the band has to fire
        res_ic = importCGMES(path = svedala, name = "Svedala_vset_ic", ignore_connected = true)
        warnings = [m for m in res_ic.messages if startswith(m, "warning:")]
        @test length(warnings) == 5
        @test all(m -> occursin("implausible voltage target", m), warnings)
        # The substituted value has to be derived from the nominal data, and the
        # message must name the unit so the defect is traceable to the source.
        @test any(m -> occursin("VATTENDRAGET_G1", m) && occursin("0.001 kV", m) && occursin("17.0 kV", m), warnings)

        # The affected buses must end up at a usable voltage rather than at zero.
        for name in ("VATTENDRAGET_G1_CNODE", "ATOMSBORG_G1_CNODE", "HÄSTSJÖ_G2_CNODE", "STUPET_G3_CNODE", "RUTHUVUD_G3_CNODE")
          idx = res_ic.net.busDict[name]
          @test res_ic.net.nodeVec[idx]._vm_pu > 0.5
        end

        # The band is a default, not law: widening it accepts the delivery's own
        # values again, which is what a data set that legitimately regulates
        # outside the default band needs.
        res_wide = importCGMES(path = svedala, name = "Svedala_vset_wide", ignore_connected = true, vset_min_pu = 0.0, vset_max_pu = 1.0e9)
        @test isempty([m for m in res_wide.messages if startswith(m, "warning:")])
        @test_throws ArgumentError importCGMES(path = svedala, name = "Svedala_vset_bad", vset_min_pu = 1.2, vset_max_pu = 0.8)
      end

      # The combined seven-area delivery exercises the whole multi-area chain:
      # AC X-node equivalents discarded, the DC border crossing split per side
      # (Stage-0 HVDC), inline VSC converters mapped as fixed PCC injections,
      # and — unlike every single-area model — a power flow that actually
      # converges and reproduces the delivery's own SV state. Offline-safe:
      # fetchCGMESTestSet builds the ZIP from the per-model cache without
      # network when all member folders are present.
      rcg_members = ("Svedala", "Espheim", "Portheim", "Belgovia", "Galia", "Britheim", "Nordheim")
      if all(isdir(joinpath(rcg_root, m)) for m in rcg_members) && isdir(joinpath(rcg_root, "_shared"))
        @testset "ReliCapGrid combined model (multi-area assembly + Stage-0 HVDC)" begin
          out = mktempdir()
          z = Sparlectra.CGMESImporter.fetchCGMESTestSet("relicapgrid_cgm"; outdir = out)
          res = importCGMES(path = z, name = "ReliCapGrid_CGM", multi_slack = true)
          msgs = res.messages

          # cancelling AC X-node pairs are discarded ...
          @test count(m -> occursin("has both sides present", m), msgs) == 25
          # ... the non-cancelling DC crossing is split per side instead
          @test count(m -> occursin("node split per side", m), msgs) == 1
          @test any(m -> occursin("BP_SD-EH_DC1", m) && occursin("node split per side", m), msgs)
          @test haskey(res.net.busDict, "BP_SD-EH_DC1") && haskey(res.net.busDict, "BP_SD-EH_DC1@2")
          # inline VSC stations become fixed PCC injections
          @test count(m -> occursin("HVDC converter: fixed PCC injection", m), msgs) == 2
          # Stage 0 names the detected pair without attaching a controller
          @test count(m -> occursin("HVDC pair detected", m), msgs) == 1
          @test isempty(Sparlectra._hvdc_pair_controllers(res.net))

          # #297 Draft B: paired_control derives the DC-crossing transfer and
          # loss from the two SSH operating points (-109.118 / +100.02 MW)
          res_paired = importCGMES(path = z, name = "ReliCapGrid_CGM_paired", multi_slack = true, hvdc_mode = :paired_control)
          paired = Sparlectra._hvdc_pair_controllers(res_paired.net)
          @test length(paired) == 1
          @test isapprox(only(paired).p_transfer_mw, 109.118; atol = 1e-9)
          @test isapprox(only(paired).loss_mw, 9.098; atol = 1e-9)
          @test any(m -> occursin("HVDC pair attached as controller", m), res_paired.messages)

          # the assembled model must solve and reproduce the delivery's SV
          # state. Q-limits stay off here: the active-set path currently cannot
          # handle this case (documented solver follow-up), and the quality
          # gates below are what the import is accountable for.
          net = res.net
          ite, erg = runpf!(net, 100, 1e-6, 0; method = :rectangular, islands_enabled = true,
                            islands_diagnostic_continue_after_failure = true,
                            qlimits_enabled = false, opt_flatstart = true)
          @test erg == 0
          cmp = compareWithSV(res)
          @test cmp.n > 250
          @test cmp.rms_dva < 10.0          # measured 6.3° — angle profile reproduced
          @test cmp.flows.n > 500           # flow comparison works despite isolated buses
          @test cmp.flows.rms_dp < 25.0     # measured 16.3 MW
        end
      end
    else
      @info "ReliCapGrid Svedala not cached — skipping CGMES 3.0 fixture (fetchCGMESTestSet(\"svedala\") enables it)"
    end
  end
end
