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
# purpose: tests the CGMES importer on synthetic in-memory deliveries (RDF
#          reader semantics and profile classification, import-failure
#          analysis, base-voltage inference, mapping of taps, controllers,
#          machines and shunts) and on the checked-in deliveries under
#          test/fixtures/cgmes that tools/gen_cgmes_fixtures.jl exports
#          from the shipped cases: folder and zip import, SV and source
#          reproduction, the Stage-2 OLTC controller, the PST sign, the
#          Q-limit hull, SCF export, the service runs (API dispatch,
#          import analysis, short circuit, fixed-reference self-check).
#          Nothing is downloaded (see docs/src/cgmes_import.md).

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

# The service layer takes case FILES: pack the four profile files of a
# checked-in delivery (cgmes_fixture_dir) into one zip at test time; no zip
# is checked in.
function _pack_cgmes_fixture_zip(case::AbstractString)::String
  z = joinpath(mktempdir(), string(case, "_CGMES.zip"))
  ZipArchives.ZipWriter(z) do w
    for f in sort(readdir(cgmes_fixture_dir(case); join = true))
      ZipArchives.zip_newfile(w, basename(f))
      write(w, read(f))
    end
  end
  return z
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

    # The service legs on in-memory deliveries: the synthetic 3W delivery
    # packed as a ZIP is the healthy case, the synthetic EQ with an absent
    # prerequisite the broken one, and a MATPOWER case covers the non-CGMES
    # rejections (the healthy runs on a checked-in delivery are above).
    # --- self-built deliveries (test/fixtures/cgmes) --------------------------
    #
    # The three checked-in deliveries are exports of shipped cases, written
    # by tools/gen_cgmes_fixtures.jl: sp_case14 (OLTC voltage controller),
    # sp_case118 (54 machines with asymmetric Q limits, synchronous
    # condensers), sp_casePST (phase-shifting transformer, bus link). Their
    # SV profile is the source's solution with every regulated unit holding
    # its target (Q-limits off in the fixture solve), so the profile is the
    # solution of the delivery itself and an import that rebuilds the model
    # reproduces it to floating-point noise. Every solve below keeps
    # Q-limits off for the same reason; what the Q-limit switching makes of
    # the delivery is a solver question, not an importer one.
    fixture_cases = [
      # (case, buses, branches, links, slack, source loader)
      ("sp_case14", 14, 17, 0, "Ostheim_110", () -> importSCF(joinpath(dirname(@__DIR__), "data", "scf", "sp_case14.scf.json"))),
      ("sp_case118", 118, 186, 0, "69", () -> createNetFromMatPowerFile(filename = joinpath(dirname(@__DIR__), "data", "mpower", "sp_case118.m"), flatstart = false, enable_pq_gen_controllers = true, bus_shunt_model = :admittance, matpower_shift_sign = 1.0, matpower_shift_unit = :deg, matpower_ratio = :normal, tap_changer_model = :ideal)),
      ("sp_casePST", 9, 8, 1, "1", () -> importSCF(joinpath(dirname(@__DIR__), "data", "scf", "sp_casePST.scf.json"))),
    ]
    solve_fixture!(net; kwargs...) = runpf!(net, 60, 1e-8, 0; method = :rectangular, qlimits_enabled = false, kwargs...)

    @testset "fixture deliveries import from the folder and from a zip" begin
      for (case, nbus, nbr, nlink, slack, _) in fixture_cases
        dir = cgmes_fixture_dir(case)
        s = summarizeCGMES(path = dir)
        @test s.version == "2.4.15"
        @test s.unresolved_count == 0
        @test !s.boundary_missing_hint
        hist = Dict(s.class_histogram)
        @test hist[:TopologicalNode] == nbus
        @test hist[:SvVoltage] == nbus
        res = importCGMES(path = dir, name = case)
        @test length(res.net.nodeVec) == nbus
        @test length(res.net.branchVec) == nbr
        @test length(res.net.linkVec) == nlink
        @test res.slack_bus == slack
        @test isempty(res.no_sv_buses)
        # a self-built delivery is complete: no importer warning, no skip
        @test !any(m -> startswith(m, "warning:") || occursin("skip", m), res.messages)
        # the service layer takes case FILES, so the zip packed from the
        # same four profiles must rebuild the same network
        rz = importCGMES(path = _pack_cgmes_fixture_zip(case), name = string(case, "_zip"))
        @test length(rz.net.nodeVec) == nbus
        @test length(rz.net.branchVec) == nbr
        @test Set(keys(rz.net.busDict)) == Set(keys(res.net.busDict))
        @test rz.slack_bus == slack
      end
    end

    # The importer's exit criterion on its own exporter's output: the SV
    # profile and the source net's solution are reproduced to noise. Both
    # start points are checked, the SV state (the solve must not move) and
    # a flat start (the model itself must reach the same state).
    @testset "fixture SV profile and source solution are reproduced" begin
      for (case, nbus, _, _, _, loader) in fixture_cases
        res = importCGMES(path = cgmes_fixture_dir(case), name = case)
        ite, erg = solve_fixture!(res.net)
        @test erg == 0
        @test ite <= 2
        cmp = compareWithSV(res)
        println("      ", case, ": SV comparison max |dvm| = ", cmp.max_dvm, " pu, max |dva| = ", cmp.max_dva, " deg, flows max |dp| = ", cmp.flows.max_dp, " MW")
        @test cmp.n == nbus
        @test cmp.max_dvm < 1e-6
        @test cmp.max_dva < 1e-6
        @test cmp.flows.n > nbus
        @test cmp.flows.max_dp < 1e-6
        @test cmp.flows.max_dq < 1e-6
        flat = importCGMES(path = cgmes_fixture_dir(case), name = string(case, "_flat"))
        _, erg_flat = solve_fixture!(flat.net; opt_flatstart = true)
        @test erg_flat == 0
        @test compareWithSV(flat).max_dvm < 1e-6
        # and the source case, solved the same way, lands on the same state
        # bus by bus (the TopologicalNode names are the source bus names)
        src = loader()
        @test solve_fixture!(src)[2] == 0
        @test Set(keys(src.busDict)) == Set(keys(res.net.busDict))
        for (bus, i) in src.busDict
          j = res.net.busDict[bus]
          @test isapprox(src.nodeVec[i]._vm_pu, res.net.nodeVec[j]._vm_pu; atol = 1e-6)
          @test isapprox(src.nodeVec[i]._va_deg, res.net.nodeVec[j]._va_deg; atol = 1e-6)
        end
      end
    end

    # Stage 2 semantics on the exporter's TapChangerControl: without
    # tap_control the delivery imports as a fixed-tap model, with it the
    # OLTC of sp_case14 comes back as the voltage controller the source case
    # carries (same target bus, target and deadband), and the control loop
    # meets the target within its deadband.
    @testset "OLTC controller of sp_case14 arrives as a tap controller (Stage 2)" begin
      dir = cgmes_fixture_dir("sp_case14")
      stage1 = importCGMES(path = dir, name = "s1")
      @test isempty(collect(Sparlectra._tap_controllers(stage1.net)))
      res = importCGMES(path = dir, tap_control = true, name = "s2")
      ctrls = collect(Sparlectra._tap_controllers(res.net))
      @test length(ctrls) == 1
      c = only(ctrls)
      @test c.mode == :voltage
      @test c.control_ratio
      @test !c.control_phase
      @test c.enabled
      @test isempty(c.followers)
      @test any(m -> startswith(m, "tap control:") && occursin("voltage", m), res.messages)
      src = importSCF(joinpath(dirname(@__DIR__), "data", "scf", "sp_case14.scf.json"))
      sc = only(collect(Sparlectra._tap_controllers(src)))
      @test c.target_bus == sc.target_bus
      @test isapprox(something(c.target_vm_pu, NaN), sc.target_vm_pu; atol = 1e-9)
      @test isapprox(c.deadband_vm_pu, sc.deadband_vm_pu; atol = 1e-9)
      # the controlled transformer connects the same buses as in the source
      # and keeps its tap machinery (range and step from the RatioTapChanger)
      names = Sparlectra._bus_name_by_idx(res.net)
      src_names = Sparlectra._bus_name_by_idx(src)
      br = Sparlectra._find_trafo_branch(res.net, c.trafo)
      sbr = Sparlectra._find_trafo_branch(src, sc.trafo)
      @test Set([names[br.fromBus], names[br.toBus]]) == Set([src_names[sbr.fromBus], src_names[sbr.toBus]])
      @test br.has_ratio_tap
      @test br.tap_step > 0.0
      cfg = SparlectraConfig(powerflow = PowerFlowConfig(max_iter = 30, tol = 1e-8), output = OutputConfig(logfile_results = :off), control = ControlConfig())
      runres = run_sparlectra(net = res.net, config = cfg)
      @test runres.numerical_converged
      @test latest_control_result(res.net).status == :converged
      @test c.converged
      @test abs(c.achieved_vm_pu - c.target_vm_pu) <= c.deadband_vm_pu
    end

    # The PST of sp_casePST travels as a single-step PhaseTapChangerLinear on
    # end 1; its shift comes back with the source sign, and the closed
    # Breaker comes back as the bus link.
    @testset "PST of sp_casePST keeps its tap-side sign, bus link survives" begin
      res = importCGMES(path = cgmes_fixture_dir("sp_casePST"), name = "pst")
      src = importSCF(joinpath(dirname(@__DIR__), "data", "scf", "sp_casePST.scf.json"))
      shifted(net) = [br for br in net.branchVec if br.phase_shift_deg != 0.0]
      @test length(shifted(src)) == 1
      @test length(shifted(res.net)) == 1
      sbr, ibr = only(shifted(src)), only(shifted(res.net))
      @test sbr.phase_shift_deg < 0.0
      @test isapprox(ibr.phase_shift_deg, sbr.phase_shift_deg; atol = 1e-9)
      @test isapprox(ibr.angle, sbr.angle; atol = 1e-9)
      names = Sparlectra._bus_name_by_idx(res.net)
      src_names = Sparlectra._bus_name_by_idx(src)
      @test (names[ibr.fromBus], names[ibr.toBus]) == (src_names[sbr.fromBus], src_names[sbr.toBus])
      @test any(m -> occursin("bus link", m) && occursin("Breaker", m), res.messages)
      @test length(res.net.linkVec) == 1
      @test only(res.net.linkVec).status == 1
    end

    # sp_case118 carries 54 synchronous machines, among them condensers
    # (P = 0, voltage regulated) and asymmetric minQ/maxQ pairs. The
    # condensers arrive as regulated zero-P units. The Q limits arrive as
    # the Stage-1 hull of both sign readings (the ENTSO-E sets are
    # inconsistent about the convention), so an asymmetric pair widens to
    # its symmetric envelope; this pins that contract rather than limit
    # fidelity.
    @testset "sp_case118 machines: synchronous condensers and the Q-limit hull" begin
      res = importCGMES(path = cgmes_fixture_dir("sp_case118"), name = "c118")
      src = fixture_cases[2][6]()
      by_bus(net) = Dict(Sparlectra._bus_name_by_idx(net)[Sparlectra.getPosumerBusIndex(ps)] => ps for ps in net.prosumpsVec if Sparlectra.isGenerator(ps))
      sg, ig = by_bus(src), by_bus(res.net)
      @test length(sg) == 54
      @test Set(keys(ig)) == Set(keys(sg))
      condensers = [b for (b, ps) in sg if something(ps.pVal, 0.0) == 0.0 && ps.isRegulated == true]
      @test !isempty(condensers)
      for b in condensers
        @test something(ig[b].pVal, NaN) == 0.0
        @test ig[b].isRegulated == true
        @test isapprox(something(ig[b].vm_pu, NaN), sg[b].vm_pu; atol = 1e-9)
      end
      asymmetric = 0
      for (b, ps) in sg
        (ps.minQ === nothing || ps.maxQ === nothing) && continue
        abs(ps.minQ + ps.maxQ) > 1e-9 && (asymmetric += 1)
        @test ig[b].minQ == min(ps.minQ, -ps.maxQ)
        @test ig[b].maxQ == max(ps.maxQ, -ps.minQ)
      end
      @test asymmetric > 0
    end

    # SCF export of a delivery (issue #342): a delivery-sourced net reaches
    # the case format with its source identity intact, the file is
    # deterministic, and the case re-imports to the same solution.
    @testset "SCF export of a fixture delivery round-trips" begin
      res = importCGMES(path = cgmes_fixture_dir("sp_case14"), name = "scf14")
      @test solve_fixture!(res.net)[2] == 0
      root = Sparlectra.net_to_scf(res.net; source_format = "cgmes", source_reference = "sp_case14 fixture")
      extra = root["sparlectra"]["extra"]
      tagged = [e for e in values(extra) if get(e, "source_id_kind", "") == "cgmes_mrid"]
      @test !isempty(tagged)
      @test all(e -> !startswith(String(e["external_id"]), "#"), tagged)
      names = Set(String(e["name"]) for e in values(extra))
      @test all(bus -> bus in names, keys(res.net.busDict))
      d = mktempdir()
      a = exportSCF(res.net; file = joinpath(d, "a.scf.json"), source_format = "cgmes")
      b = exportSCF(res.net; file = joinpath(d, "b.scf.json"), source_format = "cgmes")
      @test read(a, String) == read(b, String)
      back = importSCF(a)
      @test length(back.nodeVec) == 14
      @test solve_fixture!(back)[2] == 0
      for (bus, i) in res.net.busDict
        j = back.busDict[bus]
        @test isapprox(res.net.nodeVec[i]._vm_pu, back.nodeVec[j]._vm_pu; atol = 1e-6)
        @test isapprox(res.net.nodeVec[i]._va_deg, back.nodeVec[j]._va_deg; atol = 1e-6)
      end
    end

    # The framework/API path: a delivery zip as casefile, dispatch by
    # auto-detection, cgmes.log next to run.log, and the SV comparison
    # artifacts of the default start (auto resolves to sv on a delivery
    # that carries SvVoltage).
    @testset "run_sparlectra_api dispatch on a fixture delivery zip" begin
      z = _pack_cgmes_fixture_zip("sp_case14")
      out = mktempdir()
      r = run_sparlectra_api(casefile = z, output_dir = out)
      @test r.status == :succeeded
      @test r.metadata["input_format_detected"] == "cgmes"
      @test r.metadata["cgmes_version"] == "2.4.15"
      @test r.metadata["cgmes_buses"] == 14
      @test r.metadata["cgmes_slack_bus"] == "Ostheim_110"
      @test r.metadata["cgmes_messages"] > 0
      logfile = joinpath(out, "cgmes.log")
      @test isfile(logfile)
      text = read(logfile, String)
      for section in ("# CGMES import report", "## Class histogram", "## Network built", "## Short-circuit source data", "## Importer messages", "## SV comparison")
        @test occursin(section, text)
      end
      @test r.metadata["cgmes_start_values"] == "sv"
      @test r.metadata["cgmes_sv_compare_status"] == "converged"
      @test r.metadata["cgmes_sv_compare_n"] == 14
      @test r.metadata["cgmes_sv_compare_max_dvm"] < 1e-6
      @test isfile(joinpath(out, "sv_compare.csv"))
      @test isfile(joinpath(out, "sv_compare_flows.csv"))
      @test occursin("CGMES start values: sv (auto: delivery carries SvVoltage for", read(joinpath(out, "run.log"), String))

      # cgmes_import.start_values selects the start state and wins over a
      # hostile power_flow.flatstart; both starts land on the SV solution
      for (mode, hostile_flatstart) in (("flat", "false"), ("sv", "true"))
        out_m = mktempdir()
        cfg_m = joinpath(out_m, "c.yaml")
        write(cfg_m, "config_version: 1\npower_flow:\n  flatstart: " * hostile_flatstart * "\n  tol: 1.0e-10\ncgmes_import:\n  start_values: " * mode * "\n")
        rm_ = run_sparlectra_api(casefile = z, config_file = cfg_m, output_dir = out_m, case_format = :cgmes)
        @test rm_.status == :succeeded
        run_log = read(joinpath(out_m, "run.log"), String)
        @test occursin("CGMES start values: " * mode, run_log)
        @test occursin("overrides: power_flow.flatstart=" * hostile_flatstart, run_log)
        @test rm_.metadata["cgmes_start_values"] == mode
        @test rm_.metadata["cgmes_sv_compare_max_dvm"] < 1e-6
        sv_rows = readlines(joinpath(out_m, "sv_compare.csv"))
        @test sv_rows[1] == "bus,vm_pu,sv_vm_pu,dvm,va_deg,sv_va_deg,dva,dva_aligned"
        @test length(sv_rows) - 1 == 14
      end

      # the report survives a failed solve: it is written right after the
      # import, not after the power flow (one flat-start iteration cannot
      # reach 1e-14)
      out_f = mktempdir()
      cfg_f = joinpath(out_f, "c.yaml")
      write(cfg_f, "config_version: 1\npower_flow:\n  max_iter: 1\n  tol: 1.0e-14\ncgmes_import:\n  start_values: flat\n")
      rf = run_sparlectra_api(casefile = z, config_file = cfg_f, output_dir = out_f)
      @test rf.status != :succeeded
      @test isfile(joinpath(out_f, "cgmes.log"))
      @test occursin("# CGMES import report", read(joinpath(out_f, "cgmes.log"), String))
    end

    # The Web UI "Analyze import" and "Short circuit" buttons on a fixture.
    # The exporter writes SynchronousMachine objects without x''_d, ratedS
    # or ratedU (its sc_source keyword takes a CGMES harvest only), so the
    # short-circuit run finds machines but no usable source data on them:
    # it completes with every row flagged as a lower bound and says so by
    # name, both in the result reason and in the substitution warnings.
    @testset "import_analysis_mode and short_circuit_mode service runs on a fixture" begin
      root = mktempdir()
      cfg = joinpath(root, "c.yaml")
      write(cfg, "config_version: 1\npower_flow:\n  max_iter: 40\n")
      z = _pack_cgmes_fixture_zip("sp_case14")

      ia = start_powerflow_run(Dict("casefile" => z, "config_file" => cfg, "output_root" => root, "import_analysis_mode" => true))
      @test ia["success"] === true
      @test ia["metadata"]["run_mode"] == "import_analysis"
      @test ia["metadata"]["import_analysis_missing_dependencies"] == 0
      report = read(joinpath(root, ia["run_id"], "import_analysis.txt"), String)
      @test occursin("Supplied models:", report)
      @test occursin("Verdict:", report)

      # the button gate sees machines, the run reports what they lack
      @test Sparlectra._webui_case_has_short_circuit_data(cgmes_fixture_dir("sp_case14"))
      sc = run_with_expected_warnings(() -> start_powerflow_run(Dict("casefile" => z, "config_file" => cfg, "output_root" => root, "short_circuit_mode" => true)), ["has no usable x''_d", "has no usable ratedS", "has no usable ratedU"])
      @test sc["success"] === true
      @test sc["reason"] == "short_circuit_flagged_lower_bound"
      @test occursin("lower bound", sc["message"])
      @test sc["metadata"]["run_mode"] == "short_circuit"
      @test sc["metadata"]["sc_case_rows"] == 14
      @test sc["metadata"]["sc_flagged_rows"] == 14
      @test sc["metadata"]["sc_max_ik_kA"] > 0.0
      rid = sc["run_id"]
      max_rows = readlines(joinpath(root, rid, "short_circuit_max.csv"))
      @test max_rows[1] == "bus,vn_kV,island,status,c,zk_ohm,rx_ratio,ik_kA,sk_MVA,kappa,ip_kA,flagged,reasons"
      @test length(max_rows) - 1 == 14
      @test all(occursin("true", split(l, ',')[12]) for l in max_rows[2:end])
      @test isfile(joinpath(root, rid, "short_circuit_min.csv"))
      run_log = read(joinpath(root, rid, "run.log"), String)
      @test occursin("Short-circuit run", run_log)
      @test occursin("SynchronousMachine", run_log)
    end

    # Fixed-reference self-check on a delivery: the SV voltages reach the
    # solver verbatim (a base config with flatstart: true must not wipe
    # them) and the residual/attribution artifacts are written. The fixture
    # SV is a solution, so the start residual is noise; a flat start would
    # leave MW-sized residuals on the loaded buses.
    @testset "CGMES fixed-reference self-check (SV start, artifacts)" begin
      out = mktempdir()
      cfg = joinpath(out, "c.yaml")
      write(cfg, "config_version: 1\npower_flow:\n  flatstart: true\n")
      rsc = run_fixed_reference_self_check(casefile = cgmes_fixture_dir("sp_casePST"), config_file = cfg, output_dir = out, case_format = :cgmes)
      @test rsc.status == :succeeded
      @test rsc.raw_result !== nothing
      @test rsc.raw_result.iterations == 1
      @test rsc.metadata["cgmes_no_sv_buses"] == 0
      summary = read(joinpath(out, "self_check.log"), String)
      @test occursin("start values taken verbatim from import", summary)
      m = match(r"start_state_residual_inf: ([0-9.eE+-]+)", summary)
      @test m !== nothing
      @test parse(Float64, m.captures[1]) < 1e-8
      residuals = readlines(joinpath(out, "self_check_residuals.csv"))
      @test residuals[1] == "bus_id,bus_name,vn_kV,bus_type,vm_pu_start,va_deg_start,p_residual,q_residual,has_sv,n_transformer_terminals,n_shunts"
      @test length(residuals) - 1 == 9
      @test all(split(l, ',')[9] == "true" for l in residuals[2:end])
      @test any(parse(Int, split(l, ',')[10]) > 0 for l in residuals[2:end])
      @test occursin("no-SV buses: 0", read(joinpath(out, "cgmes.log"), String))
    end

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
    # fail with explicit reasons, never with empty tables (the run on a
    # checked-in delivery is above).
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
