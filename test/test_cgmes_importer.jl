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
# purpose: tests the CGMES importer on synthetic in-memory deliveries: RDF
#          reader semantics and profile classification, import-failure
#          analysis, base-voltage inference, mapping of taps, controllers,
#          machines and shunts, and the service runs that take a delivery
#          ZIP. Every input is built in memory (see docs/src/cgmes_import.md);
#          the tests on downloaded ENTSO-E and ReliCapGrid deliveries were
#          removed and return on self-built deliveries from the exporter.

using Sparlectra.CGMESImporter:
  CGMESStore, loadCGMES, summarizeCGMES, objectsOf, countOf, num, str, boolval, enumval, ref, unresolvedReferences, importFailureAnalysis
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

    @testset "SV angle alignment reaches the parent island detection" begin
      # Regression 2026-09-06. compareWithSV removes ONE angle offset per
      # island, because an angle is only defined up to a constant per island.
      # That call sits in the submodule CGMESImporter while detect_ac_islands
      # lives in the parent, so it must be qualified: an unqualified name
      # compiles fine and throws UndefVarError at RUNTIME, where the
      # surrounding catch turned it into a single global offset without any
      # visible sign. The per-island alignment was dead from the day it was
      # written, and it only surfaced when the reported improvement was
      # reproduced on the run that motivated it instead of recalculated.
      @test Core.eval(Sparlectra.CGMESImporter, :(Sparlectra.detect_ac_islands)) isa Function
      # the unqualified name is genuinely absent, which is why it has to be
      # written out; if this ever becomes true the qualification stops being
      # load-bearing and the test above still holds
      @test !isdefined(Sparlectra.CGMESImporter, :detect_ac_islands)
      src = read(joinpath(pkgdir(Sparlectra), "src", "adapters", "cgmes", "cgmes_report.jl"), String)
      @test occursin("Sparlectra.detect_ac_islands(net)", src)
      # and the fallback must not be silent any more: a comparison that lost
      # the island map reports secondary islands as tens of degrees off, so
      # the degraded mode has to announce itself
      @test occursin("island detection failed", src)
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
      res = run_with_expected_warnings(() -> importCGMES(path = nobv_dir, name = "nobv_on", infer_base_voltages = true), ["data defect"])
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
    # consistency is checked here without touching the network or the cache;
    # the fetch itself needs GitHub and is not exercised by the tests.
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

    # The service legs below used to run on fetched ENTSO-E ZIPs. The tests
    # on downloaded deliveries are gone (self-built deliveries from the
    # exporter take their place in a later round); what stays is what an
    # in-memory delivery can carry: the synthetic 3W delivery packed as a
    # ZIP is the healthy case, the synthetic EQ with an absent prerequisite
    # the broken one, and a MATPOWER case covers the non-CGMES rejections.
    @testset "import_analysis_mode service run (synthetic deliveries)" begin
      root = mktempdir()
      cfgia = joinpath(root, "c.yaml")
      write(cfgia, "config_version: 1\npower_flow:\n  max_iter: 40\n")

      # Healthy delivery: succeeded, importable, artifact written. The
      # service accepts case FILES, so the synthetic delivery is packed as
      # a zip.
      okcase = joinpath(mktempdir(), "synth3w.zip")
      ZipArchives.ZipWriter(okcase) do w
        for f in readdir(_cgmes_synth3w_dir(); join = true)
          ZipArchives.zip_newfile(w, basename(f))
          write(w, read(f, String))
        end
      end
      okresp = start_powerflow_run(Dict("casefile" => okcase, "config_file" => cfgia, "output_root" => root, "import_analysis_mode" => true))
      @test okresp["success"] === true
      @test okresp["metadata"]["run_mode"] == "import_analysis"
      @test okresp["metadata"]["import_analysis_missing_dependencies"] == 0
      okreport = read(joinpath(root, okresp["run_id"], "import_analysis.txt"), String)
      @test occursin("Supplied models:", okreport)
      @test occursin("Verdict:", okreport)

      # A delivery declaring an absent prerequisite: FAILED with the
      # explicit reason, analysis artifact still written.
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
      mpcase = abspath(joinpath(dirname(@__DIR__), "data", "mpower", "warmup_casePST.m"))
      mpia = start_powerflow_run(Dict("casefile" => mpcase, "config_file" => cfgia, "output_root" => root, "import_analysis_mode" => true))
      @test mpia["success"] === false
      @test mpia["reason"] == "import_analysis_requires_cgmes"

      # Mode exclusivity (checked before the case is touched).
      exia = start_powerflow_run(Dict("casefile" => okcase, "config_file" => cfgia, "output_root" => root, "import_analysis_mode" => true, "short_circuit_mode" => true))
      @test exia["success"] === false
      @test occursin("excludes", exia["message"])
    end

    # The Web UI "Short circuit" button's service path: negative cases must
    # fail with explicit reasons, never with empty tables. The plausibility
    # sweep on a real delivery left with the downloaded deliveries.
    @testset "short_circuit_mode service run rejects non-CGMES cases and mode mixing" begin
      root = mktempdir()
      cfgsc = joinpath(root, "c.yaml")
      write(cfgsc, "config_version: 1\npower_flow:\n  max_iter: 40\n")
      mpcase = abspath(joinpath(dirname(@__DIR__), "data", "mpower", "warmup_casePST.m"))
      # MATPOWER case: explicit rejection, no artifacts
      mp = start_powerflow_run(Dict("casefile" => mpcase, "config_file" => cfgsc, "output_root" => root, "short_circuit_mode" => true))
      @test mp["success"] === false
      @test mp["reason"] == "short_circuit_requires_cgmes"
      # diagnose_mode and short_circuit_mode are mutually exclusive
      both = start_powerflow_run(Dict("casefile" => mpcase, "config_file" => cfgsc, "output_root" => root, "short_circuit_mode" => true, "diagnose_mode" => true))
      @test both["success"] === false
      @test occursin("mutually exclusive", both["message"])
    end

    # cgmes_import.start_values is a CGMES-only key: a MATPOWER run ignores
    # it completely (no decision line, no SV artifacts, no metadata keys).
    @testset "cgmes_import.start_values has no effect on a MATPOWER run" begin
      mp_out = mktempdir()
      mp_cfg = joinpath(mp_out, "c.yaml")
      write(mp_cfg, "config_version: 1\ncgmes_import:\n  start_values: sv\n")
      mp = run_sparlectra_api(casefile = abspath(joinpath(dirname(@__DIR__), "data", "mpower", "warmup_casePST.m")), config_file = mp_cfg, output_dir = mp_out)
      @test mp.status == :succeeded
      @test !isfile(joinpath(mp_out, "sv_compare.csv"))
      @test !haskey(mp.metadata, "cgmes_start_values")
      @test !haskey(mp.metadata, "cgmes_sv_compare_status")
      @test !occursin("CGMES start values", read(joinpath(mp_out, "run.log"), String))
    end

    # The voltage-setpoint plausibility band is validated before any file is
    # read: an inverted band is an argument error, not a silent no-op.
    @testset "inverted vset band is rejected" begin
      @test_throws ArgumentError importCGMES(path = _cgmes_synth3w_dir(), name = "vset_bad", vset_min_pu = 1.2, vset_max_pu = 0.8)
    end
  end
end
