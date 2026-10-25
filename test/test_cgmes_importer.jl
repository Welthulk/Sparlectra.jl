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
#          data/cgmes_demo that tools/gen_cgmes_fixtures.jl exports
#          from the shipped cases: folder and zip import, SV and source
#          reproduction, the Stage-2 OLTC controller, the PST sign, the
#          Q-limit hull, SCF export, the service runs (API dispatch,
#          import analysis, short circuit, fixed-reference self-check).
#          On the ENTSO-E SmallGrid under test/fixtures/cgmes_smallgrid:
#          the mRID-to-element map and the import from a loaded store.
#          Nothing is downloaded (see docs/src/cgmes_import.md).

using Sparlectra.CGMESImporter:
  CGMESStore, loadCGMES, summarizeCGMES, objectsOf, countOf, num, str, boolval, enumval, ref, unresolvedReferences, importFailureAnalysis
import ZipArchives
isdefined(@__MODULE__, :scf_roundtrip_field_diffs) || include("test_scf_support.jl")

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

isdefined(@__MODULE__, :_cgmes_synth3w_dir) || include(joinpath(@__DIR__, "cgmes_synth3w_fixture.jl"))
isdefined(@__MODULE__, :_cgmes_synth_eb_dir) || include(joinpath(@__DIR__, "cgmes_synth_eb_fixture.jl"))

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

# Synthetic 2W delivery (EQ+SSH+TP+SV) with a RatioTapChangerTable: slack ENI
# at a 110 kV bus, a 110/20 kV transformer whose end-1 RatioTapChanger stands
# at step 5 of 1..5, load at the 20 kV bus. The table ratio at step 5 (1.045)
# differs from the linear step formula (1 + 2 * 1.25 % = 1.025), and the
# table is not uniform. The SV profile is the closed-form two-bus solution
# with the table ratio: ideal ratio at the from side, series impedance on the
# end-2 base (the importer's 2W branch model), slack at 1.0 pu. `target_kV`
# is the TapChangerControl target, `point_x` a reactance correction in
# percent written on every table point (the importer applies it to the series
# reactance; the SV profile is computed without it, so a delivery with
# `point_x` is not checked against its SV). `control_enabled = false` writes
# the changer with controlEnabled false, so the import regulates it not.
const _RTC_TABLE_RATIOS = Dict(1 => 0.97, 2 => 0.985, 3 => 1.0, 4 => 1.02, 5 => 1.045)

function _cgmes_synth_rtc_table_dir(; target_kV::Float64 = 20.4, point_x::Float64 = 0.0, control_enabled::Bool = true)::String
  t = _RTC_TABLE_RATIOS[5]
  z = complex(0.6, 12.0) / (110.0^2 / 100.0)
  y = 1.0 / z
  s_load = complex(20.0, 5.0) / 100.0
  v1 = complex(1.0, 0.0)
  v2 = v1 / t
  for _ in 1:200
    v2 = v1 / t - z * conj(s_load / v2)
  end
  i = (v1 / t - v2) * y
  # power into the equipment at each terminal, MW/Mvar
  s_from = (v1 / t) * conj(i) * 100.0
  s_to = -v2 * conj(i) * 100.0
  fmt(x) = string(round(x; digits = 12))
  points = join(["<cim:RatioTapChangerTablePoint rdf:ID=\"_rtct_p$(k)\"><cim:TapChangerTablePoint.step>$(k)</cim:TapChangerTablePoint.step><cim:TapChangerTablePoint.ratio>$(r)</cim:TapChangerTablePoint.ratio><cim:TapChangerTablePoint.x>$(point_x)</cim:TapChangerTablePoint.x><cim:RatioTapChangerTablePoint.RatioTapChangerTable rdf:resource=\"#_rtct\"/></cim:RatioTapChangerTablePoint>" for (k, r) in sort(collect(_RTC_TABLE_RATIOS))], "\n")
  hdr = """<?xml version="1.0" encoding="UTF-8"?>
<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:cim="http://iec.ch/TC57/2013/CIM-schema-cim16#" xmlns:md="http://iec.ch/TC57/61970-552/ModelDescription/1#">"""
  eq = """$(hdr)
<md:FullModel rdf:about="urn:uuid:eqrtct"><md:Model.profile>http://entsoe.eu/CIM/EquipmentCore/3/1</md:Model.profile></md:FullModel>
<cim:BaseVoltage rdf:ID="_bv110"><cim:BaseVoltage.nominalVoltage>110</cim:BaseVoltage.nominalVoltage></cim:BaseVoltage>
<cim:BaseVoltage rdf:ID="_bv20"><cim:BaseVoltage.nominalVoltage>20</cim:BaseVoltage.nominalVoltage></cim:BaseVoltage>
<cim:ExternalNetworkInjection rdf:ID="_eni"><cim:IdentifiedObject.name>GRID</cim:IdentifiedObject.name><cim:ExternalNetworkInjection.referencePriority>1</cim:ExternalNetworkInjection.referencePriority></cim:ExternalNetworkInjection>
<cim:Terminal rdf:ID="_t_eni"><cim:Terminal.ConductingEquipment rdf:resource="#_eni"/><cim:ACDCTerminal.sequenceNumber>1</cim:ACDCTerminal.sequenceNumber></cim:Terminal>
<cim:EnergyConsumer rdf:ID="_load"><cim:IdentifiedObject.name>LOAD_LV</cim:IdentifiedObject.name></cim:EnergyConsumer>
<cim:Terminal rdf:ID="_t_load"><cim:Terminal.ConductingEquipment rdf:resource="#_load"/><cim:ACDCTerminal.sequenceNumber>1</cim:ACDCTerminal.sequenceNumber></cim:Terminal>
<cim:PowerTransformer rdf:ID="_t2w"><cim:IdentifiedObject.name>T2W</cim:IdentifiedObject.name></cim:PowerTransformer>
<cim:PowerTransformerEnd rdf:ID="_t2w_e1"><cim:PowerTransformerEnd.PowerTransformer rdf:resource="#_t2w"/><cim:TransformerEnd.endNumber>1</cim:TransformerEnd.endNumber><cim:PowerTransformerEnd.ratedU>110</cim:PowerTransformerEnd.ratedU><cim:PowerTransformerEnd.ratedS>40</cim:PowerTransformerEnd.ratedS><cim:PowerTransformerEnd.r>0.6</cim:PowerTransformerEnd.r><cim:PowerTransformerEnd.x>12.0</cim:PowerTransformerEnd.x><cim:TransformerEnd.Terminal rdf:resource="#_t_t2w_1"/></cim:PowerTransformerEnd>
<cim:PowerTransformerEnd rdf:ID="_t2w_e2"><cim:PowerTransformerEnd.PowerTransformer rdf:resource="#_t2w"/><cim:TransformerEnd.endNumber>2</cim:TransformerEnd.endNumber><cim:PowerTransformerEnd.ratedU>20</cim:PowerTransformerEnd.ratedU><cim:PowerTransformerEnd.ratedS>40</cim:PowerTransformerEnd.ratedS><cim:PowerTransformerEnd.r>0</cim:PowerTransformerEnd.r><cim:PowerTransformerEnd.x>0</cim:PowerTransformerEnd.x><cim:TransformerEnd.Terminal rdf:resource="#_t_t2w_2"/></cim:PowerTransformerEnd>
<cim:Terminal rdf:ID="_t_t2w_1"><cim:Terminal.ConductingEquipment rdf:resource="#_t2w"/><cim:ACDCTerminal.sequenceNumber>1</cim:ACDCTerminal.sequenceNumber></cim:Terminal>
<cim:Terminal rdf:ID="_t_t2w_2"><cim:Terminal.ConductingEquipment rdf:resource="#_t2w"/><cim:ACDCTerminal.sequenceNumber>2</cim:ACDCTerminal.sequenceNumber></cim:Terminal>
<cim:RatioTapChanger rdf:ID="_rtc"><cim:IdentifiedObject.name>T2W_OLTC</cim:IdentifiedObject.name><cim:RatioTapChanger.TransformerEnd rdf:resource="#_t2w_e1"/><cim:RatioTapChanger.stepVoltageIncrement>1.25</cim:RatioTapChanger.stepVoltageIncrement><cim:TapChanger.lowStep>1</cim:TapChanger.lowStep><cim:TapChanger.highStep>5</cim:TapChanger.highStep><cim:TapChanger.neutralStep>3</cim:TapChanger.neutralStep><cim:RatioTapChanger.RatioTapChangerTable rdf:resource="#_rtct"/><cim:TapChanger.TapChangerControl rdf:resource="#_tcc"/></cim:RatioTapChanger>
<cim:RatioTapChangerTable rdf:ID="_rtct"><cim:IdentifiedObject.name>T2W_TABLE</cim:IdentifiedObject.name></cim:RatioTapChangerTable>
$(points)
<cim:TapChangerControl rdf:ID="_tcc"><cim:IdentifiedObject.name>T2W_CTRL</cim:IdentifiedObject.name><cim:RegulatingControl.mode rdf:resource="http://iec.ch/TC57/2013/CIM-schema-cim16#RegulatingControlModeKind.voltage"/><cim:RegulatingControl.Terminal rdf:resource="#_t_load"/></cim:TapChangerControl>
</rdf:RDF>
"""
  ssh = """$(hdr)
<md:FullModel rdf:about="urn:uuid:sshrtct"><md:Model.profile>http://entsoe.eu/CIM/SteadyStateHypothesis/1/1</md:Model.profile></md:FullModel>
<cim:ExternalNetworkInjection rdf:about="#_eni"><cim:ExternalNetworkInjection.p>$(fmt(-real(s_from)))</cim:ExternalNetworkInjection.p><cim:ExternalNetworkInjection.q>$(fmt(-imag(s_from)))</cim:ExternalNetworkInjection.q></cim:ExternalNetworkInjection>
<cim:EnergyConsumer rdf:about="#_load"><cim:EnergyConsumer.p>20</cim:EnergyConsumer.p><cim:EnergyConsumer.q>5</cim:EnergyConsumer.q></cim:EnergyConsumer>
<cim:RatioTapChanger rdf:about="#_rtc"><cim:TapChanger.step>5</cim:TapChanger.step><cim:TapChanger.controlEnabled>$(control_enabled)</cim:TapChanger.controlEnabled></cim:RatioTapChanger>
<cim:TapChangerControl rdf:about="#_tcc"><cim:RegulatingControl.enabled>true</cim:RegulatingControl.enabled><cim:RegulatingControl.targetValue>$(target_kV)</cim:RegulatingControl.targetValue><cim:RegulatingControl.targetDeadband>0.5</cim:RegulatingControl.targetDeadband></cim:TapChangerControl>
<cim:Terminal rdf:about="#_t_eni"><cim:ACDCTerminal.connected>true</cim:ACDCTerminal.connected></cim:Terminal>
<cim:Terminal rdf:about="#_t_load"><cim:ACDCTerminal.connected>true</cim:ACDCTerminal.connected></cim:Terminal>
<cim:Terminal rdf:about="#_t_t2w_1"><cim:ACDCTerminal.connected>true</cim:ACDCTerminal.connected></cim:Terminal>
<cim:Terminal rdf:about="#_t_t2w_2"><cim:ACDCTerminal.connected>true</cim:ACDCTerminal.connected></cim:Terminal>
</rdf:RDF>
"""
  tp = """$(hdr)
<md:FullModel rdf:about="urn:uuid:tprtct"><md:Model.profile>http://entsoe.eu/CIM/Topology/4/1</md:Model.profile></md:FullModel>
<cim:TopologicalNode rdf:ID="_tn_hv"><cim:IdentifiedObject.name>HV_BUS</cim:IdentifiedObject.name><cim:TopologicalNode.BaseVoltage rdf:resource="#_bv110"/></cim:TopologicalNode>
<cim:TopologicalNode rdf:ID="_tn_lv"><cim:IdentifiedObject.name>LV_BUS</cim:IdentifiedObject.name><cim:TopologicalNode.BaseVoltage rdf:resource="#_bv20"/></cim:TopologicalNode>
<cim:Terminal rdf:about="#_t_eni"><cim:Terminal.TopologicalNode rdf:resource="#_tn_hv"/></cim:Terminal>
<cim:Terminal rdf:about="#_t_load"><cim:Terminal.TopologicalNode rdf:resource="#_tn_lv"/></cim:Terminal>
<cim:Terminal rdf:about="#_t_t2w_1"><cim:Terminal.TopologicalNode rdf:resource="#_tn_hv"/></cim:Terminal>
<cim:Terminal rdf:about="#_t_t2w_2"><cim:Terminal.TopologicalNode rdf:resource="#_tn_lv"/></cim:Terminal>
</rdf:RDF>
"""
  sv = """$(hdr)
<md:FullModel rdf:about="urn:uuid:svrtct"><md:Model.profile>http://entsoe.eu/CIM/StateVariables/4/1</md:Model.profile></md:FullModel>
<cim:SvVoltage rdf:ID="_sv_hv"><cim:SvVoltage.TopologicalNode rdf:resource="#_tn_hv"/><cim:SvVoltage.v>110</cim:SvVoltage.v><cim:SvVoltage.angle>0</cim:SvVoltage.angle></cim:SvVoltage>
<cim:SvVoltage rdf:ID="_sv_lv"><cim:SvVoltage.TopologicalNode rdf:resource="#_tn_lv"/><cim:SvVoltage.v>$(fmt(20.0 * abs(v2)))</cim:SvVoltage.v><cim:SvVoltage.angle>$(fmt(rad2deg(angle(v2))))</cim:SvVoltage.angle></cim:SvVoltage>
<cim:SvTapStep rdf:ID="_sv_tap"><cim:SvTapStep.TapChanger rdf:resource="#_rtc"/><cim:SvTapStep.position>5</cim:SvTapStep.position></cim:SvTapStep>
<cim:SvPowerFlow rdf:ID="_svpf_1"><cim:SvPowerFlow.Terminal rdf:resource="#_t_t2w_1"/><cim:SvPowerFlow.p>$(fmt(real(s_from)))</cim:SvPowerFlow.p><cim:SvPowerFlow.q>$(fmt(imag(s_from)))</cim:SvPowerFlow.q></cim:SvPowerFlow>
<cim:SvPowerFlow rdf:ID="_svpf_2"><cim:SvPowerFlow.Terminal rdf:resource="#_t_t2w_2"/><cim:SvPowerFlow.p>$(fmt(real(s_to)))</cim:SvPowerFlow.p><cim:SvPowerFlow.q>$(fmt(imag(s_to)))</cim:SvPowerFlow.q></cim:SvPowerFlow>
</rdf:RDF>
"""
  dir = mktempdir()
  write(joinpath(dir, "synthrtct_EQ.xml"), eq)
  write(joinpath(dir, "synthrtct_SSH.xml"), ssh)
  write(joinpath(dir, "synthrtct_TP.xml"), tp)
  write(joinpath(dir, "synthrtct_SV.xml"), sv)
  return dir
end

# The service layer takes case FILES: pack the four profile files of a
# checked-in delivery (cgmes_fixture_dir) into one zip at test time; no zip
# is checked in.
_pack_cgmes_fixture_zip(case::AbstractString)::String = _pack_cgmes_dir_zip(cgmes_fixture_dir(case), case)

# the same for a delivery written to a directory at test time
function _pack_cgmes_dir_zip(dir::AbstractString, stem::AbstractString)::String
  z = joinpath(mktempdir(), string(stem, "_CGMES.zip"))
  ZipArchives.ZipWriter(z) do w
    for f in sort(readdir(dir; join = true))
      ZipArchives.zip_newfile(w, basename(f))
      write(w, read(f))
    end
  end
  return z
end

# A copy of the sp_case14 delivery whose SV states an SvTapStep four steps
# off the SSH step of the uncontrolled ratio tap changer; the SV flows were
# solved at the SSH step.
function _sp_case14_off_sv_tap_step()
  dir = mktempdir()
  for f in readdir(cgmes_fixture_dir("sp_case14"))
    cp(joinpath(cgmes_fixture_dir("sp_case14"), f), joinpath(dir, f))
  end
  rtc = only(tc for tc in objectsOf(loadCGMES(dir), :RatioTapChanger) if !haskey(tc.refs, :TapChangerControl))
  ssh_step = Int(num(rtc, :step))
  svfile = joinpath(dir, "sp_case14_SV.xml")
  write(svfile, replace(read(svfile, String), "</rdf:RDF>" => """<cim:SvTapStep rdf:ID="_sv_tap_step"><cim:SvTapStep.position>$(ssh_step - 4)</cim:SvTapStep.position><cim:SvTapStep.TapChanger rdf:resource="#$(rtc.mrid)"/></cim:SvTapStep></rdf:RDF>"""))
  return (dir = dir, ssh_step = ssh_step, sv_step = ssh_step - 4)
end

function _cgmes_synth_store()::CGMESStore
  dir = mktempdir()
  write(joinpath(dir, "synth_EQ.xml"), _CGMES_SYNTH_EQ)
  write(joinpath(dir, "synth_SSH.xml"), _CGMES_SYNTH_SSH)
  write(joinpath(dir, "synth_DIFF.xml"), _CGMES_SYNTH_DIFF)
  return loadCGMES(dir)
end

# A copy of a loaded SmallGrid store with DY objects added in memory under
# a file of profile DY (SmallGrid ships no DY):
#   TwinBrch SM (200 MVA): SynchronousMachineTimeConstantReactance with
#     inertia 4.25 and damping 0.5, GovSteam1 with k 25, t1 0.4, mwbase 100,
#     an ExcIEEEAC4A on that machine model
#   Breed SM (750 MVA): SynchronousMachineSimplified with inertia 6, no
#     damping, GovSteam0 with r 0.05 and mwbase 0
#   DeerCrk SM: a disabled SynchronousMachineTimeConstantReactance
#   a SynchronousMachineEquivalentCircuit and a GovHydro1 on a machine the
#   delivery does not have
#   a LoadAggregate with a LoadStatic, referred to by an EQ EnergyConsumer,
#   and a LoadAggregate nothing refers to
function _cgmes_dy_store(loaded::CGMESStore)::CGMESStore
  CI = Sparlectra.CGMESImporter
  s = deepcopy(loaded)
  add!(mrid, class, attrs, refs = Dict{Symbol,String}()) = (s.objects[mrid] = CI.CIMObject(mrid, class, Dict{Symbol,String}(attrs), Dict{Symbol,String}(refs), "dy_test.xml"); push!(get!(Vector{String}, s.byclass, class), mrid))
  sm(name) = only(o for o in objectsOf(s, :SynchronousMachine) if str(o, :name) == name).mrid
  add!("_dy_m1", :SynchronousMachineTimeConstantReactance, Dict(:enabled => "true", :inertia => "4.25", :damping => "0.5"), Dict(:SynchronousMachine => sm("TwinBrch SM")))
  add!("_dy_g1", :GovSteam1, Dict(:enabled => "true", :k => "25", :t1 => "0.4", :mwbase => "100"), Dict(:SynchronousMachineDynamics => "_dy_m1"))
  add!("_dy_x1", :ExcIEEEAC4A, Dict(:enabled => "true", :ka => "200"), Dict(:SynchronousMachineDynamics => "_dy_m1"))
  add!("_dy_m2", :SynchronousMachineSimplified, Dict(:enabled => "true", :inertia => "6"), Dict(:SynchronousMachine => sm("Breed SM")))
  add!("_dy_g2", :GovSteam0, Dict(:enabled => "true", :r => "0.05", :t1 => "0.3", :mwbase => "0"), Dict(:SynchronousMachineDynamics => "_dy_m2"))
  add!("_dy_m3", :SynchronousMachineTimeConstantReactance, Dict(:enabled => "false", :inertia => "3"), Dict(:SynchronousMachine => sm("DeerCrk SM")))
  add!("_dy_m4", :SynchronousMachineEquivalentCircuit, Dict(:enabled => "true", :inertia => "5"), Dict(:SynchronousMachine => "_dy_no_such_machine"))
  add!("_dy_g4", :GovHydro1, Dict(:enabled => "true", :rperm => "0.04", :mwbase => "100"), Dict(:SynchronousMachineDynamics => "_dy_m4"))
  add!("_dy_la", :LoadAggregate, Dict(:name => "LA"))
  add!("_dy_ls", :LoadStatic, Dict(:name => "LS"), Dict(:LoadAggregate => "_dy_la"))
  add!("_dy_la_orphan", :LoadAggregate, Dict(:name => "LA2"))
  first(objectsOf(s, :EnergyConsumer)).refs[:LoadDynamics] = "_dy_la"
  push!(s.files, CI.CGMESFileInfo("dy_test.xml", :FullModel, ["http://entsoe.eu/CIM/Dynamics/3/1"], Set([:DY]), "2.4.15", 11, false, ""))
  return s
end

# A copy of a loaded SmallGrid store with small islands on new nodes, each
# node with an SvVoltage of 1.03 pu, for the islands without a slack:
#   A: a2 - a1, a load of 11.3 MW and an EquivalentInjection with P = 0 and
#      Q = 0.7 MVAr: no source
#   B: b1 - b2, a machine with referencePriority 7 and 5 MW at b1, an
#      EquivalentInjection with 20 MW, a load and a machine with
#      referencePriority 1 drawing 3 MW at b2
#   C: c1 - c2, EquivalentInjections with +10 MW and -70.3 MW, no maxP:
#      every bus has degree 1 (the fragment of a network equivalent)
#   D: d1 - d2 - d3, network injections at d1 (10 of 100 MW) and d2 (20 of
#      50 MW), an EquivalentInjection with 40 MW at d3, a load and a network
#      injection drawing 5 MW with maxP 200 MW at d2
#   E: e1 - e2 - e3, no maxP: an EquivalentInjection with +10 MW at e1
#      (degree 1), one drawing 70.3 MW at e2 (degree 2), a load at e3
#   F: f1 - f2, an EquivalentInjection drawing 5 MW and a load: no source
#   G: g1 - g2, a machine drawing 3 MW at g1 whose voltage RegulatingControl
#      (1.02 pu) sits on the terminal of a load at g2
#   H: h1 - h2, an EquivalentInjection drawing 8 MW and a load that feeds
#      10 MW in: no source
# None of these units but the machine on G regulates its voltage.
# With `linked = true` the store gets only island J instead:
#   J: j1 - j2 = jx = j3 - j4, closed Breakers j2 - jx and jx - j3, jx has
#      no branch: EquivalentInjections with 10 MW at j1 and 20 MW at j4,
#      loads at j2 and j3, no regulation
function _cgmes_f2_island_store(loaded::CGMESStore; linked::Bool = false)::CGMESStore
  CI = Sparlectra.CGMESImporter
  s = deepcopy(loaded)
  host = first(objectsOf(s, :EnergyConsumer))
  host_t = only(t for t in objectsOf(s, :Terminal) if get(t.refs, :ConductingEquipment, "") == host.mrid)
  bv = s.objects[host_t.refs[:TopologicalNode]].refs[:BaseVoltage]
  vn = num(s.objects[bv], :nominalVoltage)
  add!(mrid, class, attrs, refs = Dict{Symbol,String}()) = (s.objects[mrid] = CI.CIMObject(mrid, class, Dict{Symbol,String}(attrs), Dict{Symbol,String}(refs), "f2_test"); push!(get!(Vector{String}, s.byclass, class), mrid))
  term!(eq, tn, seq = 1) = add!(string(eq, "_t", seq), :Terminal, Dict(:sequenceNumber => string(seq), :connected => "true"), Dict(:ConductingEquipment => eq, :TopologicalNode => tn))
  function line!(mrid, a, b)
    add!(mrid, :ACLineSegment, Dict(:name => uppercase(mrid), :r => "1.0", :x => "5.0", :bch => "0"), Dict(:BaseVoltage => bv))
    term!(mrid, a, 1)
    term!(mrid, b, 2)
  end
  unit!(mrid, class, tn, attrs) = (add!(mrid, class, merge(Dict(:name => uppercase(mrid)), attrs)); term!(mrid, tn))
  function node!(tn)
    add!(tn, :TopologicalNode, Dict(:name => uppercase(tn)), Dict(:BaseVoltage => bv))
    add!(string(tn, "_sv"), :SvVoltage, Dict(:v => string(1.03 * vn), :angle => "2.0"), Dict(:TopologicalNode => tn))
  end
  if linked
    foreach(node!, ("_f2_j1", "_f2_j2", "_f2_jx", "_f2_j3", "_f2_j4"))
    line!("_f2_lj1", "_f2_j1", "_f2_j2")
    line!("_f2_lj2", "_f2_j3", "_f2_j4")
    for (sw, a, b) in (("_f2_sw_j2", "_f2_j2", "_f2_jx"), ("_f2_sw_j3", "_f2_jx", "_f2_j3"))
      add!(sw, :Breaker, Dict(:name => uppercase(sw), :open => "false"))
      term!(sw, a, 1)
      term!(sw, b, 2)
    end
    unit!("_f2_ei_j1", :EquivalentInjection, "_f2_j1", Dict(:p => "-10", :q => "0"))
    unit!("_f2_ei_j4", :EquivalentInjection, "_f2_j4", Dict(:p => "-20", :q => "0"))
    unit!("_f2_load_j2", :EnergyConsumer, "_f2_j2", Dict(:p => "12", :q => "0"))
    unit!("_f2_load_j3", :EnergyConsumer, "_f2_j3", Dict(:p => "15", :q => "0"))
    return s
  end
  foreach(node!, ("_f2_a1", "_f2_a2", "_f2_b1", "_f2_b2", "_f2_c1", "_f2_c2", "_f2_d1", "_f2_d2", "_f2_d3", "_f2_e1", "_f2_e2", "_f2_e3", "_f2_f1", "_f2_f2", "_f2_g1", "_f2_g2", "_f2_h1", "_f2_h2"))
  line!("_f2_la", "_f2_a1", "_f2_a2")
  unit!("_f2_ei_a", :EquivalentInjection, "_f2_a1", Dict(:p => "0", :q => "-0.7"))
  unit!("_f2_load_a", :EnergyConsumer, "_f2_a2", Dict(:p => "11.3", :q => "2"))
  # CGMES load sign convention: a unit with p < 0 injects
  line!("_f2_lb", "_f2_b1", "_f2_b2")
  unit!("_f2_sm_b", :SynchronousMachine, "_f2_b1", Dict(:p => "-5", :q => "0", :ratedS => "50", :minQ => "-20", :maxQ => "20", :referencePriority => "7"))
  unit!("_f2_ei_b", :EquivalentInjection, "_f2_b2", Dict(:p => "-20", :q => "0"))
  unit!("_f2_load_b", :EnergyConsumer, "_f2_b2", Dict(:p => "24", :q => "1"))
  unit!("_f2_sm_b2", :SynchronousMachine, "_f2_b2", Dict(:p => "3", :q => "0", :ratedS => "50", :minQ => "-20", :maxQ => "20", :referencePriority => "1"))
  line!("_f2_lc", "_f2_c1", "_f2_c2")
  unit!("_f2_ei_c1", :EquivalentInjection, "_f2_c1", Dict(:p => "-10", :q => "0"))
  unit!("_f2_ei_c2", :EquivalentInjection, "_f2_c2", Dict(:p => "70.3", :q => "0"))
  line!("_f2_ld1", "_f2_d1", "_f2_d2")
  line!("_f2_ld2", "_f2_d2", "_f2_d3")
  unit!("_f2_eni_d1", :ExternalNetworkInjection, "_f2_d1", Dict(:p => "-10", :q => "0", :minP => "0", :maxP => "100"))
  unit!("_f2_eni_d2", :ExternalNetworkInjection, "_f2_d2", Dict(:p => "-20", :q => "0", :minP => "0", :maxP => "50"))
  unit!("_f2_ei_d3", :EquivalentInjection, "_f2_d3", Dict(:p => "-40", :q => "0"))
  unit!("_f2_load_d", :EnergyConsumer, "_f2_d2", Dict(:p => "75", :q => "3"))
  unit!("_f2_eni_d2b", :ExternalNetworkInjection, "_f2_d2", Dict(:p => "5", :q => "0", :minP => "-200", :maxP => "200"))
  line!("_f2_le1", "_f2_e1", "_f2_e2")
  line!("_f2_le2", "_f2_e2", "_f2_e3")
  unit!("_f2_ei_e1", :EquivalentInjection, "_f2_e1", Dict(:p => "-10", :q => "0"))
  unit!("_f2_ei_e2", :EquivalentInjection, "_f2_e2", Dict(:p => "70.3", :q => "0"))
  unit!("_f2_load_e", :EnergyConsumer, "_f2_e3", Dict(:p => "1", :q => "0"))
  line!("_f2_lf", "_f2_f1", "_f2_f2")
  unit!("_f2_ei_f", :EquivalentInjection, "_f2_f1", Dict(:p => "5", :q => "0"))
  unit!("_f2_load_f", :EnergyConsumer, "_f2_f2", Dict(:p => "2", :q => "0"))
  line!("_f2_lg", "_f2_g1", "_f2_g2")
  unit!("_f2_load_g", :EnergyConsumer, "_f2_g2", Dict(:p => "2", :q => "0"))
  add!("_f2_rc_g", :RegulatingControl, Dict(:name => "_F2_RC_G", :mode => "http://iec.ch/TC57/2013/CIM-schema-cim16#RegulatingControlModeKind.voltage", :enabled => "true", :targetValue => string(1.02 * vn)), Dict(:Terminal => "_f2_load_g_t1"))
  add!("_f2_sm_g", :SynchronousMachine, Dict(:name => "_F2_SM_G", :p => "3", :q => "0", :ratedS => "50", :minQ => "-20", :maxQ => "20"), Dict(:RegulatingControl => "_f2_rc_g"))
  term!("_f2_sm_g", "_f2_g1")
  line!("_f2_lh", "_f2_h1", "_f2_h2")
  unit!("_f2_ei_h", :EquivalentInjection, "_f2_h1", Dict(:p => "8", :q => "0"))
  unit!("_f2_load_h", :EnergyConsumer, "_f2_h2", Dict(:p => "-10", :q => "0"))
  return s
end

function run_cgmes_importer_tests()
  @testset "CGMES importer" begin (function ()
    store = _cgmes_synth_store()

    @testset "profile classification and version" begin (function ()
      @test store.version == "2.4.15"
      eqinfo = only(filter(f -> occursin("EQ", f.name), store.files))
      @test eqinfo.header == :FullModel
      @test eqinfo.profiles == Set([:EQ, :EQ_SC])   # per-file profile *set* (A2/D-7)
      @test !eqinfo.skipped
    end)() end

    @testset "DifferenceModel files are skipped with reason (D-6)" begin (function ()
      diffinfo = only(filter(f -> occursin("DIFF", f.name), store.files))
      @test diffinfo.skipped
      @test diffinfo.header == :DifferenceModel
      @test !isempty(diffinfo.skip_reason)
      @test !haskey(store.objects, "_should_not_appear")
    end)() end

    @testset "import failure analysis names supplied models and gaps" begin (function ()
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
    end)() end

    @testset "infer_base_voltages reconstructs missing nominal voltages" begin (function ()
      # The RVC delivery with its BaseVoltage catalog stripped — the exact
      # shape of a real delivery whose catalog lives in a missing boundary EQ.
      nobv_dir = mktempdir()
      write(joinpath(nobv_dir, "nobv_EQ.xml"), replace(_CGMES_SYNTH_RVC_EQ, "<cim:BaseVoltage rdf:ID=\"_bv110\"><cim:BaseVoltage.nominalVoltage>110</cim:BaseVoltage.nominalVoltage></cim:BaseVoltage>" => ""))
      write(joinpath(nobv_dir, "nobv_TP.xml"), replace(_CGMES_SYNTH_RVC_TP, "<cim:TopologicalNode.BaseVoltage rdf:resource=\"#_bv110\"/>" => ""))
      write(joinpath(nobv_dir, "nobv_SSH.xml"), replace(_CGMES_SYNTH_RVC_SSH_TEMPLATE, "__SM2CONN__" => "false"))
      write(joinpath(nobv_dir, "nobv_SV.xml"), _CGMES_SYNTH_RVC_SV)

      # Without the option the import aborts with the typed analysis error.
      @test_throws Sparlectra.CGMESImporter.CGMESImportError importCGMES(path = nobv_dir, name = "nobv_off")
      # A service run writes the analysis of that error into the import
      # report, the file its message points at: the power-flow service adds
      # what summarizeCGMES still reads, the short-circuit service the
      # analysis alone.
      local pf_out = mktempdir()
      local pf_fail = run_sparlectra_api(casefile = nobv_dir, output_dir = pf_out, case_format = :cgmes)
      @test pf_fail.status == :failed
      @test occursin("cgmes_import_report.md", something(pf_fail.message, ""))
      local pf_text = cgmes_report_plain(pf_out)
      @test startswith(pf_text, "CGMES import report: IMPORT FAILED\n")
      @test occursin("\nImport analysis\nCGMES import analysis\n", pf_text)
      @test occursin("\nnobv_EQ.xml | EQ | ", pf_text)
      local sc_out = mktempdir()
      local sc_fail = SparlectraApp._run_short_circuit_service(nobv_dir, Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, sc_out, "sc_nobv")
      @test sc_fail.status == :failed
      @test occursin("\nImport analysis\nCGMES import analysis\n", cgmes_report_plain(sc_out))

      # With the option every bus level is reconstructed from the SV state
      # (111.8–115.5 kV snap to 110), the substitution is summarized as one
      # warning, and the network still solves.
      res = run_with_expected_warnings(() -> importCGMES(path = nobv_dir, name = "nobv_on", infer_base_voltages = true), ["data defect"])
      @test length(res.net.nodeVec) == 3
      @test all(Sparlectra.getNodeVn(b) == 110.0 for b in res.net.nodeVec)
      @test count(m -> occursin("inferred base voltages", m), res.messages) == 1
      @test last(runpf!(res.net, 30, 1e-8, 0)) == 0
    end)() end

    @testset "self-loop line maps to a shunt notice, not a branch" begin (function ()
      # Both terminals of a line on ONE topological node (busbar link modeled
      # as a line) used to trip the branch constructor's from!=to assertion.
      # Its charging stays as a shunt at the bus.
      loop_dir = mktempdir()
      loop_eq = replace(
        _CGMES_SYNTH_RVC_EQ,
        "<cim:ACLineSegment rdf:ID=\"_l_ab\">" =>
          "<cim:ACLineSegment rdf:ID=\"_l_loop\"><cim:IdentifiedObject.name>L_LOOP</cim:IdentifiedObject.name><cim:ACLineSegment.r>0.01</cim:ACLineSegment.r><cim:ACLineSegment.x>0.05</cim:ACLineSegment.x><cim:ACLineSegment.bch>1e-5</cim:ACLineSegment.bch></cim:ACLineSegment>\n<cim:Terminal rdf:ID=\"_t_loop_1\"><cim:Terminal.ConductingEquipment rdf:resource=\"#_l_loop\"/><cim:ACDCTerminal.sequenceNumber>1</cim:ACDCTerminal.sequenceNumber></cim:Terminal>\n<cim:Terminal rdf:ID=\"_t_loop_2\"><cim:Terminal.ConductingEquipment rdf:resource=\"#_l_loop\"/><cim:ACDCTerminal.sequenceNumber>2</cim:ACDCTerminal.sequenceNumber></cim:Terminal>\n<cim:ACLineSegment rdf:ID=\"_l_ab\">",
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
      # the charging shunt is addressed under the derived key, the line
      # itself has no element: B = 1e-5 S * 110 kV^2 = 0.121 MVAr at BUS_A
      @test !haskey(loop_res.element_by_mrid, "_l_loop")
      kind, k = loop_res.element_by_mrid["_l_loop#charging"]
      @test kind == :shunt
      sh = loop_res.net.shuntVec[k]
      @test sh.busIdx == loop_res.net.busDict["BUS_A"]
      @test isapprox(sh.B_shunt, 0.121 / loop_res.net.baseMVA; rtol = 1e-12)
      @test last(runpf!(loop_res.net, 30, 1e-8, 0)) == 0
    end)() end

    @testset "EquivalentBranch on one level, across levels, phase-stamped and as self-loop" begin (function ()
      # Without these branches BUS_C has no element and BUS_D hangs on its
      # load alone, cut off from the slack.
      res = importCGMES(path = _cgmes_synth_eb_dir(), name = "synth_eb")
      net = res.net
      @test count(m -> occursin("EquivalentBranch EB_LOOP", m) && occursin("to itself", m), res.messages) == 1
      summary = filter(m -> startswith(m, "notice: EquivalentBranch:"), res.messages)
      @test length(summary) == 1
      @test occursin("3 mapped (2 same voltage level, 1 across voltage levels), 2 skipped (1 on one topological node, 1 open or out of service, 0 without topological node), 2 asymmetric (r21, x21 differ from r, x) phase-stamped (ratio 1, shift arg(z21/z12)/2), 0 of them with |z21| other than |z12|", only(summary))
      @test !occursin("max ||z21|/|z12| - 1|", only(summary))
      @test "_eb_loop" in res.skipped_equipment && "_eb_open" in res.skipped_equipment
      @test !haskey(res.element_by_mrid, "_eb_loop") && !haskey(res.element_by_mrid, "_eb_open")
      bus(name) = net.busDict[name]
      # same level: pi branch without charging, r and x on the 110 kV base
      kind, ks = res.element_by_mrid["_eb_same"]
      @test kind == :branch
      bs = net.branchVec[ks]
      @test (bs.fromBus, bs.toBus) == (bus("BUS_B"), bus("BUS_C"))
      @test isapprox(bs.r_pu, 1.21 * net.baseMVA / 110.0^2; rtol = 1e-12) && isapprox(bs.x_pu, 6.05 * net.baseMVA / 110.0^2; rtol = 1e-12)
      @test bs.b_pu == 0.0 && bs.g_pu == 0.0
      # r21, x21 a pure rotation of r, x, on one level (EB_ASYM, -12 degrees)
      # and across levels (EB_CROSS, +8 degrees): ratio 1, shift
      # arg(z21/z12)/2, all on the base of terminal 1. For |z21| = |z12| the
      # geometric-mean stamp equals the stamp with Y11, Y12, Y22 exact
      # (series z12 e^{j phi}, arms (1 - e^{-j phi})/z12 and
      # 1/z21 - e^{-j phi}/z12) to 1e-15
      z21s = _cgmes_synth_eb_z21(1.0)
      pu(z) = (p = Sparlectra.toPU_RXBG(r = real(z), x = imag(z), v_kv = 110.0, baseMVA = net.baseMVA); complex(p[1], p[2]))
      for (m, b1, b2, z12, z21, phi_deg) in (("_eb_asym", "BUS_B", "BUS_C", pu(complex(1.5, 7.0)), pu(z21s[1]), -6.0), ("_eb_cross", "BUS_C", "BUS_D", pu(complex(0.6, 4.84)), pu(z21s[2]), 4.0))
        kind, k = res.element_by_mrid[m]
        @test kind == :branch
        br = net.branchVec[k]
        @test (br.fromBus, br.toBus) == (bus(b1), bus(b2))
        @test br.tap_ratio == 1.0 && isapprox(br.phase_shift_deg, phi_deg; atol = 1e-9)
        # a transformer-type branch without tap changer: an empty set, so
        # the tap estimation never releases it
        @test isempty(br.se_tap_changers.changers) && Sparlectra._tap_release_mode(br) === nothing
        phi = deg2rad(br.phase_shift_deg)
        y12 = abs(1 / z12)
        @test abs(complex(br.r_pu, br.x_pu) - z12 * cis(phi)) <= 1e-15 * abs(z12)
        @test abs(complex(br.g_from_pu, br.b_from_pu) - (1 - cis(-phi)) / z12) <= 1e-15 * y12
        @test abs(complex(br.g_to_pu, br.b_to_pu) - (1 / z21 - cis(-phi) / z12)) <= 1e-15 * y12
      end
      bc = net.branchVec[res.element_by_mrid["_eb_cross"][2]]
      @test res.branch_side_of_terminal["_t_ebs_1"] == (bs.branchIdx, :from) && res.branch_side_of_terminal["_t_ebc_2"] == (bc.branchIdx, :to)
      # a stamped branch is no tap-controlled transformer
      @test isempty(Sparlectra.buildTapControllerReportRows(net))
      # the delivered values stay on the branch index as export metadata,
      # also in a copy of the net
      md = net.equivalent_branch_metadata
      @test length(md) == 3 && deepcopy(net).equivalent_branch_metadata == md
      @test md[ks] == (r = 1.21, x = 6.05, r21 = 1.21, x21 = 6.05, vn_kV = 110.0, stamped = false)
      mc = md[res.element_by_mrid["_eb_cross"][2]]
      @test mc.stamped && (mc.r, mc.x, mc.vn_kV) == (0.6, 4.84, 110.0) && complex(mc.r21, mc.x21) == z21s[2]
      # the line EB_SAME and the stamped EB_ASYM share the bus pair and k = 1
      # of their two counters; the keys stay apart
      ids = net.cgmes_ids
      @test ids[Sparlectra.CGMESImporter.cgmesKeyEquivalentBranch("BUS_B", "BUS_C", 1)] == "eb_same"
      @test ids[Sparlectra.CGMESImporter.cgmesKeyEquivalentBranchTransformer("BUS_B", "BUS_C", 1)] == "eb_asym"
      @test ids[Sparlectra.CGMESImporter.cgmesKeyEquivalentBranchTransformer("BUS_C", "BUS_D", 1)] == "eb_cross"
      # a class the importer maps: with every EquivalentBranch out of service
      # the importer skips each one and says so, and the table of classes
      # not mapped has no row for the class
      local eb_store = deepcopy(res.store)
      foreach(m -> (eb_store.objects[m].attrs[:inService] = "false"), eb_store.byclass[:EquivalentBranch])
      local eb_off = importCGMES(path = eb_store, name = "synth_eb_off")
      @test count(m -> m in eb_off.skipped_equipment, eb_store.byclass[:EquivalentBranch]) == 5
      @test !any(r -> r.class === :EquivalentBranch, Sparlectra.CGMESImporter._unmappedClasses(eb_off.store, eb_off.element_by_mrid))
      @test !any(m -> startswith(m, "notice: class not mapped: EquivalentBranch"), eb_off.messages)
      # The SV profile is the solved exact two-port model, so the branch
      # flows match at the SV state and the solved import stays there.
      for stage in (:sv_start, :solved)
        stage == :solved && @test last(runpf!(net, 30, 1e-8, 0)) == 0
        cmp = compareWithSV(res)
        rows = filter(r -> r.name in ("EB_SAME", "EB_ASYM", "EB_CROSS"), cmp.flows.rows)
        @test length(rows) == 6
        @test all(r -> abs(r.dp) < 1e-6 && abs(r.dq) < 1e-6, rows)
        @test cmp.max_dvm < 1e-6 && cmp.max_dva < 1e-6
      end
      @test isempty(net.isoNodes)

      # |z21| = 1.01 |z12|: built with the geometric mean, every entry is off
      # by sqrt(1.01) or its inverse and both row sums stay zero. At the SV
      # state terminal 1 misses by |V1| |V1 - V2| (1 - 1/sqrt(1.01)) / |z12|
      # and terminal 2 by |V2| |V1 - V2| (sqrt(1.01) - 1) / |z21|, the error
      # scales with the voltage difference. The first stamp (Y11, Y12, Y22
      # exact, the whole error on Y21, commit bc30df69) missed terminal 2 by
      # |V1| |V2| 0.01 / |z21| and ended 0.069 pu and 0.054 degrees from the
      # SV, because that error acts on V1. Bound of the solved import below:
      # below 2e-4 pu and 0.02 degrees (measured 1.76e-4 pu, 0.018 degrees).
      res1 = importCGMES(path = _cgmes_synth_eb_dir(mag = 1.01), name = "synth_eb_mag")
      @test count(m -> occursin("2 asymmetric (r21, x21 differ from r, x) phase-stamped (ratio 1, shift arg(z21/z12)/2), 2 of them with |z21| other than |z12| (max ||z21|/|z12| - 1| = 0.01, built with the geometric mean of |z12| and |z21|)", m), res1.messages) == 1
      V = _cgmes_synth_eb_state(1.01)
      z21m = _cgmes_synth_eb_z21(1.01)
      cmp1 = compareWithSV(res1)
      for (name, b1, b2, i, j, z12, z21) in (("EB_ASYM", "BUS_B", "BUS_C", 2, 3, complex(1.5, 7.0), z21m[1]), ("EB_CROSS", "BUS_C", "BUS_D", 3, 4, complex(0.6, 4.84), z21m[2]))
        rows = filter(r -> r.name == name, cmp1.flows.rows)
        @test length(rows) == 2
        d1 = abs(complex(only(filter(r -> r.bus == b1, rows)).dp, only(filter(r -> r.bus == b1, rows)).dq))
        d2 = abs(complex(only(filter(r -> r.bus == b2, rows)).dp, only(filter(r -> r.bus == b2, rows)).dq))
        @test isapprox(d1, abs(V[i]) * abs(V[i] - V[j]) * (1 - 1 / sqrt(1.01)) / abs(z12); rtol = 1e-6)
        @test isapprox(d2, abs(V[j]) * abs(V[i] - V[j]) * (sqrt(1.01) - 1) / abs(z21); rtol = 1e-6)
      end
      @test last(runpf!(res1.net, 30, 1e-8, 0)) == 0
      cmp1 = compareWithSV(res1)
      @test cmp1.max_dvm < 2e-4 && cmp1.max_dva < 0.02
    end)() end

    # A network equivalent is no equipment: a negative r or x is normal in
    # Ward and Kron reductions. The import takes it over unchanged, the CGMES
    # export and the short circuit accept it (one count line in the short
    # circuit report), and the phase stamp holds for z21 near -z12, where
    # arg(z21/z12) is near 180 degrees and the stamp angle near 90 degrees.
    @testset "EquivalentBranch with negative parameters" begin (function ()
      local C = Sparlectra.CGMESImporter
      local store = C.loadCGMES(_cgmes_synth_eb_dir())
      store.objects["_eb_same"].attrs[:r] = "-0.5"
      store.objects["_eb_same"].attrs[:r21] = "-0.5"
      local res = importCGMES(path = store, name = "eb_negative")
      local net = res.net
      local br = net.branchVec[res.element_by_mrid["_eb_same"][2]]
      @test br.r_pu == Sparlectra.toPU_RXBG(r = -0.5, x = 6.05, v_kv = 110.0, baseMVA = net.baseMVA)[1] && br.r_pu < 0.0
      @test last(runpf!(net, 30, 1e-8, 0)) == 0
      # the results name every equivalent branch as such and report the loss
      # of the negative one as computed, below zero
      Sparlectra.calcNetLosses!(net)
      local rows = Sparlectra.buildACPFlowReport(net).branches
      @test count(r -> r.branch_kind == "equivalent_branch", rows) == 3
      @test only(filter(r -> r.branch_index == br.branchIdx, rows)).p_loss_MW < 0.0
      @test Sparlectra._branch_kind_label(net, br) == "EqBr"
      @test length(writeCGMESFiles(net; path = mktempdir(), notices = String[])) == 4
      local sc = runShortCircuit!(net, res.shortcircuit)
      @test count(m -> startswith(m, Sparlectra._SC_EB_NEGATIVE_R * "1 "), sc.messages) == 1
      @test occursin(Sparlectra._SC_EB_NEGATIVE_R * "1 ", sprint(Sparlectra.printShortCircuitResult, sc))

      # three buses: slack B1, lines B1-B2 and B2-B3, the stamped equivalent
      # B1-B3 and loads at B2 and B3; the reference state solves the exact
      # two-port I1 = (V1 - V3)/z12, I3 = (V3 - V1)/z21 by Gauss-Seidel
      local z12 = 0.05 + 0.5im
      local zl = 0.01 + 0.1im
      for f in (-1.0, -0.99)
        local z21 = f * z12
        local st = C._ebPhaseStamp(z12, z21)
        @test isapprox(abs(st.phi), pi / 2; atol = 1e-12)
        local t = Sparlectra.Net(name = "eb3", baseMVA = 100.0)
        foreach(b -> addBus!(net = t, busName = b, vn_kV = 110.0), ("B1", "B2", "B3"))
        addPIModelACLine!(net = t, fromBus = "B1", toBus = "B2", r_pu = real(zl), x_pu = imag(zl), b_pu = 0.0, status = 1)
        addPIModelACLine!(net = t, fromBus = "B2", toBus = "B3", r_pu = real(zl), x_pu = imag(zl), b_pu = 0.0, status = 1)
        Sparlectra._addPIModelTrafo_by_idx!(net = t, from = 1, to = 3, r_pu = real(st.zs), x_pu = imag(st.zs), b_pu = 0.0, status = 1, ratio = 1.0, shift_deg = rad2deg(st.phi),
          g_from_pu = real(st.y0f), b_from_pu = imag(st.y0f), g_to_pu = real(st.y0t), b_to_pu = imag(st.y0t))
        addProsumer!(net = t, busName = "B1", type = "EXTERNALNETWORKINJECTION", referencePri = "B1", vm_pu = 1.0, va_deg = 0.0)
        addProsumer!(net = t, busName = "B2", type = "ENERGYCONSUMER", p = 30.0, q = 10.0)
        addProsumer!(net = t, busName = "B3", type = "ENERGYCONSUMER", p = 20.0, q = 5.0)
        local Y11, Y12, Y21, Y22 = Sparlectra.calcAdmittance(t.branchVec[end], 1.0, 1.0)
        # Y12 Y21 exact, Y11 and Y22 off by sqrt(|z12|/|z21|) and its inverse
        @test isapprox(Y12 * Y21, 1 / (z12 * z21); rtol = 1e-12)
        @test isapprox(Y11 * z12, sqrt(abs(z12) / abs(z21)); rtol = 1e-12) && isapprox(Y22 * z21, sqrt(abs(z21) / abs(z12)); rtol = 1e-12)
        local Y = zeros(ComplexF64, 3, 3)
        for (i, j) in ((1, 2), (2, 3))
          Y[i, i] += 1 / zl; Y[j, j] += 1 / zl; Y[i, j] -= 1 / zl; Y[j, i] -= 1 / zl
        end
        Y[1, 1] += 1 / z12; Y[1, 3] -= 1 / z12; Y[3, 3] += 1 / z21; Y[3, 1] -= 1 / z21
        local sload = ComplexF64[0, -(0.3 + 0.1im), -(0.2 + 0.05im)]
        local V = ComplexF64[1, 1, 1]
        for _ in 1:5000
          for i in 2:3
            V[i] = (conj(sload[i] / V[i]) - sum(Y[i, j] * V[j] for j in 1:3 if j != i)) / Y[i, i]
          end
        end
        @test maximum(abs.(Y * V .- conj.(sload ./ V))[2:3]) < 1e-12
        @test last(runpf!(t, 30, 1e-10, 0)) == 0
        @test maximum(k -> abs(t.nodeVec[k]._vm_pu - abs(V[k])), 1:3) < 1e-3
        f == -1.0 && @test maximum(k -> abs(t.nodeVec[k]._vm_pu * cis(deg2rad(t.nodeVec[k]._va_deg)) - V[k]), 1:3) < 1e-9
      end
    end)() end

    # An EquivalentBranch with |z| = 1e-8 pu, below the default threshold of
    # cgmes_import.equivalent_branch_min_z_pu (1e-6), builds no branch: its
    # buses are joined as a link, like a closed switch, with a notice. With
    # the threshold 0 (keyword) or 1e-9 (configuration) it stays a branch.
    @testset "EquivalentBranch below the impedance threshold" begin (function ()
      local C = Sparlectra.CGMESImporter
      local zstore = C.loadCGMES(_cgmes_synth_eb_dir())
      local o = zstore.objects["_eb_same"]
      o.attrs[:r] = "0"
      o.attrs[:r21] = "0"
      o.attrs[:x] = string(1e-8 * 110.0^2 / 100.0)
      local res = importCGMES(path = zstore, name = "eb_tiny")
      local kind, k = res.element_by_mrid["_eb_same"]
      @test kind == :link
      local lk = res.net.linkVec[k]
      @test Set((lk.fromBus, lk.toBus)) == Set((res.net.busDict["BUS_B"], res.net.busDict["BUS_C"])) && lk.status == 1
      @test count(m -> occursin("EquivalentBranch EB_SAME has an impedance below 1.0e-6 pu, no branch created, bus link", m), res.messages) == 1
      local summary = only(filter(m -> startswith(m, "notice: EquivalentBranch:"), res.messages))
      @test startswith(summary, "notice: EquivalentBranch: 2 mapped (1 same voltage level, 1 across voltage levels)") && occursin("1 below 1.0e-6 pu as bus links", summary)
      @test !haskey(res.net.equivalent_branch_metadata, length(res.net.branchVec) + 1) && length(res.net.equivalent_branch_metadata) == 2
      @test last(runpf!(res.net, 30, 1e-8, 0)) == 0
      @test importCGMES(path = zstore, name = "eb_tiny_off", equivalent_branch_min_z_pu = 0.0).element_by_mrid["_eb_same"][1] == :branch
      local cfg = Sparlectra.SparlectraConfig(Dict{String,Any}("cgmes_import" => Dict{String,Any}("equivalent_branch_min_z_pu" => 1e-9)))
      @test cfg.cgmes.equivalent_branch_min_z_pu == 1e-9 && Sparlectra.SparlectraConfig(Dict{String,Any}()).cgmes.equivalent_branch_min_z_pu == 1e-6
      @test importCGMES(cfg; path = zstore, name = "eb_tiny_cfg").element_by_mrid["_eb_same"][1] == :branch
      @test_throws ArgumentError Sparlectra.SparlectraConfig(Dict{String,Any}("cgmes_import" => Dict{String,Any}("equivalent_branch_min_z_pu" => -1.0)))
    end)() end

    # The star equivalent of a three-winding transformer can have a negative
    # leg r for a real unit; only the pair sums are physical. The CGMES and
    # MATPOWER exports and the short circuit take the leg over unchanged,
    # while a line or a two-winding transformer with r < 0 still aborts.
    @testset "three-winding leg with negative r" begin (function ()
      local C = Sparlectra.CGMESImporter
      local store = C.loadCGMES(_cgmes_synth3w_dir())
      # end 3: -0.002 + j0.8 ohm at 20 kV; the pair sums with end 1
      # (0.4 ohm at 220 kV, 0.0033 ohm at 20 kV) and end 2 stay positive
      store.objects["_t3w_e3"].attrs[:r] = "-0.002"
      merge!(store.objects["_eni"].attrs, Dict(:maxInitialSymShCCurrent => "10000", :maxR1ToX1Ratio => "0.1"))
      local res = importCGMES(path = store, name = "t3w_negative")
      local net = res.net
      local leg3 = net.branchVec[res.element_by_mrid["_t3w_e3"][2]]
      local leg1 = net.branchVec[res.element_by_mrid["_t3w_e1"][2]]
      local z1 = complex(Sparlectra.toPU_RXBG(r = 0.4, x = 24.0, v_kv = 220.0, baseMVA = net.baseMVA)[1:2]...)
      local z3 = complex(Sparlectra.toPU_RXBG(r = -0.002, x = 0.8, v_kv = 20.0, baseMVA = net.baseMVA)[1:2]...)
      @test leg3.r_base_pu == real(z3) < 0.0
      @test Int(leg3.branchIdx) in Sparlectra._negative_r_allowed_branches(net)
      @test last(runpf!(net, 30, 1e-8, 0)) == 0
      local dir = mktempdir()
      @test length(writeCGMESFiles(net; path = dir, notices = String[])) == 4
      local back = importCGMES(path = dir, name = "t3w_back")
      @test isapprox(back.net.branchVec[back.element_by_mrid["_t3w_e3"][2]].r_pu, leg3.r_pu; rtol = 1e-9)
      @test last(runpf!(back.net, 30, 1e-8, 0)) == 0
      local common = [n for n in keys(back.net.busDict) if haskey(net.busDict, n)]
      @test length(common) == length(net.nodeVec)
      @test maximum(n -> abs(back.net.nodeVec[back.net.busDict[n]]._vm_pu - net.nodeVec[net.busDict[n]]._vm_pu), common) < 1e-8
      writeMatpowerCasefile(net, joinpath(dir, "t3w.m"); write_solution = false)
      @test isapprox(minimum(Sparlectra.MatpowerIO.read_case(joinpath(dir, "t3w.m")).branch[:, 3]), real(z3); rtol = 1e-6)
      # short circuit: one count line, and at the LV bus the grid impedance
      # seen at the HV bus plus legs 1 and 3 (both with ratio 1)
      local sc = runShortCircuit!(net, res.shortcircuit)
      @test count(m -> startswith(m, Sparlectra._SC_3W_NEGATIVE_R * "1 "), sc.messages) == 1
      @test occursin(Sparlectra._SC_3W_NEGATIVE_R * "1 ", sprint(Sparlectra.printShortCircuitResult, sc))
      local hv = only(filter(r -> r.bus_idx == Int(leg1.toBus), sc.rows))
      local lv = only(filter(r -> r.bus_idx == Int(leg3.toBus), sc.rows))
      @test hv.status == :ok && lv.status == :ok && lv.vn_kV == 20.0
      local zg = hv.zk_ohm / (220.0^2 / net.baseMVA) * (hv.rx_ratio + im) / sqrt(1 + hv.rx_ratio^2)
      local zk = zg + z1 + z3
      @test isapprox(lv.zk_ohm / (20.0^2 / net.baseMVA), abs(zk); rtol = 1e-9)
      @test isapprox(lv.rx_ratio, real(zk) / imag(zk); rtol = 1e-9) && lv.rx_ratio > 0.0

      # a line and a two-winding transformer with r < 0 are equipment
      for kind in (:line, :transformer)
        local t = Sparlectra.Net(name = "neg_" * string(kind), baseMVA = 100.0)
        foreach(b -> addBus!(net = t, busName = b, vn_kV = 110.0), ("A", "B"))
        if kind === :line
          addPIModelACLine!(net = t, fromBus = "A", toBus = "B", r_pu = -0.001, x_pu = 0.1, b_pu = 0.0, status = 1)
        else
          Sparlectra._addPIModelTrafo_by_idx!(net = t, from = 1, to = 2, r_pu = -0.001, x_pu = 0.1, b_pu = 0.0, status = 1, ratio = 1.0, shift_deg = 0.0)
        end
        addProsumer!(net = t, busName = "A", type = "EXTERNALNETWORKINJECTION", referencePri = "A", vm_pu = 1.0, va_deg = 0.0)
        @test isempty(Sparlectra._negative_r_allowed_branches(t))
        @test_throws ErrorException writeCGMESFiles(t; path = mktempdir(), notices = String[])
        @test_throws ErrorException writeMatpowerCasefile(t, joinpath(mktempdir(), "neg.m"); write_solution = false)
        @test_throws ErrorException runShortCircuit!(t, Sparlectra.NativeShortCircuitData())
      end
    end)() end

    @testset "rdf:ID creates, literals and inherited attributes" begin (function ()
      @test countOf(store, :ACLineSegment) == 1
      line = only(objectsOf(store, :ACLineSegment))
      @test str(line, :name) == "L1"                 # IdentifiedObject.name → :name
      @test num(line, :r) == 2.5
      @test num(line, :gch) == 0.001                 # CGMES conductance present
      @test num(line, :r0) == 7.5                    # short-circuit attribute read (§7.7)
      @test num(line, :length) == 12.0               # Conductor.length → :length
      @test num(line, :missing_attr, 99.0) == 99.0
    end)() end

    @testset "references and enums" begin (function ()
      line = only(objectsOf(store, :ACLineSegment))
      bv = ref(store, line, :BaseVoltage)
      @test bv !== nothing && bv.class == :BaseVoltage
      @test num(bv, :nominalVoltage) == 110.0
      sm = only(objectsOf(store, :SynchronousMachine))
      @test enumval(sm, :operatingMode) == "SynchronousMachineOperatingMode.generator"
      @test boolval(sm, :earthing) == false
      @test any(u -> u.target == "_missing_container", unresolvedReferences(store))
    end)() end

    @testset "rdf:about overlays EQ object (SSH)" begin (function ()
      sm = only(objectsOf(store, :SynchronousMachine))
      @test num(sm, :p) == -120.0                    # from SSH overlay
      @test num(sm, :ratedS) == 150.0                # EQ value kept
      @test countOf(store, :SynchronousMachine) == 1 # no duplicate object
    end)() end

    @testset "summarizeCGMES on synthetic set" begin (function ()
      dir = mktempdir()
      write(joinpath(dir, "synth_EQ.xml"), _CGMES_SYNTH_EQ)
      s = summarizeCGMES(path = dir)
      @test s.object_count == 4
      @test s.unresolved_count >= 1
      @test !s.boundary_missing_hint                 # dangling ref is not a topology class
      @test (:ACLineSegment => 1) in s.class_histogram
    end)() end

    @testset "Result tables stay aligned with long names" begin (function ()
      # CGMES bus and branch identifiers routinely exceed the column widths;
      # @sprintf pads but never truncates, so without fitting every following
      # column would shift (reported from a MicroGrid Assembled run).
      @test Sparlectra._fitColumn("short", 25) == "short"
      @test length(Sparlectra._fitColumn("TN_Border_ST23 -> BE-Busbar_2", 25)) == 25
      @test endswith(Sparlectra._fitColumn("TN_Border_ST23 -> BE-Busbar_2", 25), "…")
      @test Sparlectra._fitColumn("exactly_twenty_five_chars", 25) == "exactly_twenty_five_chars"
    end)() end

    @testset "RatioTapChangerTable ratio at the tap position" begin (function ()
      dir = _cgmes_synth_rtc_table_dir()
      res = importCGMES(path = dir, name = "synth_rtct")
      br = res.net.branchVec[last(res.element_by_mrid["_t2w"])]
      # the table point of step 5, not the linear formula (1.025)
      @test br.ratio ≈ _RTC_TABLE_RATIOS[5] atol = 1e-12
      @test !any(m -> occursin("T2W_OLTC", m) && occursin("table", m), res.messages)
      # SV start: both transformer terminals match their SvPowerFlow
      rows0 = filter(r -> r.kind in (:branch_from, :branch_to), compareWithSV(res).flows.rows)
      @test length(rows0) == 2
      @test maximum(r -> max(abs(r.dp), abs(r.dq)), rows0) < 1e-6
      # the solve lands on the SV state
      @test last(runpf!(res.net, 30, 1e-10, 0)) == 0
      cmp = compareWithSV(res)
      @test cmp.n == 2
      @test cmp.max_dvm < 1e-8
      @test cmp.max_dva < 1e-6
      # no SvTapStep and no SSH step: the changer sits at
      # its middle position, the neutral step (moved to 2 here), or without
      # a neutral step the band middle 3 of 1..5; the ratio is the table
      # point of that step, the set of the tap estimation marks it, and one
      # notice names it
      function nostep_copy(src, neutral)
        nd = mktempdir()
        for f in readdir(src)
          txt = read(joinpath(src, f), String)
          nstep = neutral === nothing ? "" : "<cim:TapChanger.neutralStep>$(neutral)</cim:TapChanger.neutralStep>"
          write(joinpath(nd, f), replace(txt, r"<cim:SvTapStep .*?</cim:SvTapStep>" => "", "<cim:TapChanger.step>5</cim:TapChanger.step>" => "", "<cim:TapChanger.neutralStep>3</cim:TapChanger.neutralStep>" => nstep))
        end
        return nd
      end
      for (neutral, step) in ((2, 2), (nothing, 3))
        nd = nostep_copy(dir, neutral)
        nr = importCGMES(path = nd, name = "synth_rtct_nostep")
        nbr = nr.net.branchVec[last(nr.element_by_mrid["_t2w"])]
        @test nbr.ratio ≈ _RTC_TABLE_RATIOS[step] atol = 1e-12
        nch = only(nbr.se_tap_changers.changers)
        @test (nch.step, nch.step_assumed) == (Float64(step), true)
        @test count(==("notice: tap changer T2W/T2W_OLTC: step not transferred, middle position $(step) assumed"), nr.messages) == 1
        # the export has no mark for it: the ratio of the middle position
        # goes out and one notice names the transformer
        enotes = String[]
        writeCGMESFiles(nr.net; path = mktempdir(), notices = enotes)
        @test count(==("transformer 2WT_110: ratio step not transferred, the export writes the assumed middle position"), enotes) == 1
      end
      @test !only(br.se_tap_changers.changers).step_assumed && !any(m -> occursin("step not transferred", m), res.messages)
      rnotes = String[]
      writeCGMESFiles(res.net; path = mktempdir(), notices = rnotes)
      @test !any(n -> occursin("step not transferred", n), rnotes)
      # tap_control: start ratio, range and step come from the table ratios
      ctl = importCGMES(path = dir, name = "synth_rtct_ctl", tap_control = true)
      cbr = ctl.net.branchVec[last(ctl.element_by_mrid["_t2w"])]
      @test length(collect(Sparlectra._tap_controllers(ctl.net))) == 1
      @test cbr.has_ratio_tap
      @test cbr.tap_ratio ≈ _RTC_TABLE_RATIOS[5] atol = 1e-12
      @test cbr.tap_min ≈ minimum(values(_RTC_TABLE_RATIOS)) atol = 1e-12
      @test cbr.tap_max ≈ maximum(values(_RTC_TABLE_RATIOS)) atol = 1e-12
      @test cbr.tap_step ≈ (maximum(values(_RTC_TABLE_RATIOS)) - minimum(values(_RTC_TABLE_RATIOS))) / 4 atol = 1e-12
      @test cbr.tap_ratio_points ≈ sort(collect(values(_RTC_TABLE_RATIOS))) atol = 1e-12
      # a target between the voltages of two table points: the controller
      # steps from point to point and ends on a table ratio (the mean step
      # would stop at 1.0075 or 1.02625)
      cfg = SparlectraConfig(powerflow = PowerFlowConfig(max_iter = 30, tol = 1e-8), output = OutputConfig(logfile_results = :off), control = ControlConfig())
      for target_kV in (20.4, 19.73, 19.4)
        tc = importCGMES(path = _cgmes_synth_rtc_table_dir(target_kV = target_kV), name = "synth_rtct_tgt", tap_control = true)
        tbr = tc.net.branchVec[last(tc.element_by_mrid["_t2w"])]
        @test run_sparlectra(net = tc.net, config = cfg).numerical_converged
        @test latest_control_result(tc.net).status == :converged
        @test any(r -> isapprox(tbr.tap_ratio, r; atol = 1e-12), values(_RTC_TABLE_RATIOS))
      end
      # a changer without a transferred step that the controller moves off
      # its middle position 3 loses the mark of the stored set, and the
      # export names no assumed step; one it leaves there keeps both
      for (target_kV, moves) in ((19.4, true), (19.73, false))
        mc = importCGMES(path = nostep_copy(_cgmes_synth_rtc_table_dir(target_kV = target_kV), nothing), name = "synth_rtct_nostep_ctl", tap_control = true)
        mbr = mc.net.branchVec[last(mc.element_by_mrid["_t2w"])]
        @test only(mbr.se_tap_changers.changers).step_assumed
        @test run_sparlectra(net = mc.net, config = cfg).numerical_converged
        @test (mbr.tap_ratio != _RTC_TABLE_RATIOS[3], only(mbr.se_tap_changers.changers).step_assumed) == (moves, !moves)
        mnotes = String[]
        writeCGMESFiles(mc.net; path = mktempdir(), notices = mnotes)
        @test count(n -> occursin("step not transferred", n), mnotes) == (moves ? 0 : 1)
      end
      # per-step r/x/g/b corrections of the table are named once and the point
      # of the tap position (x + 5 %) applies to the series reactance
      xr = importCGMES(path = _cgmes_synth_rtc_table_dir(point_x = 5.0), name = "synth_rtct_x")
      @test count(m -> occursin("T2W_OLTC", m) && occursin("r/x/g/b", m), xr.messages) == 1
      xbr = xr.net.branchVec[last(xr.element_by_mrid["_t2w"])]
      @test xbr.x_pu ≈ br.x_pu * 1.05 rtol = 1e-12
      @test xbr.r_pu ≈ br.r_pu atol = 1e-15
      # a changer the import does not regulate (no tap_control, or
      # controlEnabled false) carries range, table ratios and impedance per
      # step as well, the start step unchanged; a controller configured
      # afterwards steps on the table
      pts = sort(collect(values(_RTC_TABLE_RATIOS)))
      for (tap_control, enabled) in ((false, true), (true, false))
        fr = importCGMES(path = _cgmes_synth_rtc_table_dir(point_x = 5.0, control_enabled = enabled), name = "synth_rtct_free", tap_control = tap_control)
        fbr = fr.net.branchVec[last(fr.element_by_mrid["_t2w"])]
        @test isempty(Sparlectra._tap_controllers(fr.net))
        @test fbr.tap_ratio ≈ _RTC_TABLE_RATIOS[5] atol = 1e-12
        @test fbr.x_pu ≈ xbr.x_pu rtol = 1e-12
        @test fbr.tap_ratio_points !== nothing && fbr.tap_ratio_points ≈ pts
        @test (fbr.tap_min, fbr.tap_max) == (pts[1], pts[end])
        t = fbr.tap_step_impedance
        @test t !== nothing && t.ratio_values ≈ pts && all(==((1.0, 1.05, 1.0, 1.0)), t.ratio_factors) && t.x_pu ≈ br.x_pu
        tap_control && continue
        lv = Sparlectra._bus_name_by_idx(fr.net)[fbr.toBus]
        addPowerTransformerControl!(fr.net; trafo = string(fbr.branchIdx), mode = :voltage, target_bus = lv, target_vm_pu = 19.73 / 20.0, deadband_vm_pu = 0.0125)
        @test run_sparlectra(net = fr.net, config = cfg).numerical_converged
        @test fbr.tap_ratio != _RTC_TABLE_RATIOS[5] && any(r -> isapprox(fbr.tap_ratio, r; atol = 1e-12), pts)
      end
    end)() end

    @testset "tap at its own terminal, impedance on the inner side" begin (function ()
      # the SV of both deliveries is computed with the tap at the terminal
      # of its end and the impedance inside; a model that puts the tap on
      # the other side of the impedance misses the terminal flows by MW
      flow_rows(res) = [r for r in compareWithSV(res).flows.rows if r.kind in (:branch_from, :branch_to)]
      # 2W tapped at end 2 only: built from end 2 to end 1, the rows match
      # at the SV state itself
      res2 = importCGMES(path = _cgmes_synth_tapside_dir(:two_winding), name = "tap2w")
      br = res2.net.branchVec[last(res2.element_by_mrid["_pt"])]
      names = Sparlectra._bus_name_by_idx(res2.net)
      @test (names[br.fromBus], names[br.toBus]) == ("BUS2", "BUS1")
      @test res2.branch_side_of_terminal["_t_pt2"][2] == :from
      rows = flow_rows(res2)
      @test length(rows) == 2
      @test all(r -> abs(r.dp) < 1e-6 && abs(r.dq) < 1e-6, rows)
      # 3W tapped at end 1: the star point has no SV, so the rows are
      # compared at the solved state, which is the SV state
      res3 = importCGMES(path = _cgmes_synth_tapside_dir(:three_winding), name = "tap3w")
      leg1 = res3.net.branchVec[last(res3.element_by_mrid["_pt_e1"])]
      @test Sparlectra._bus_name_by_idx(res3.net)[leg1.fromBus] == "BUS1"
      @test runpf!(res3.net, 30, 1e-10, 0)[2] == 0
      @test compareWithSV(res3).max_dvm < 1e-9
      rows = flow_rows(res3)
      @test length(rows) == 3
      @test all(r -> abs(r.dp) < 1e-6 && abs(r.dq) < 1e-6, rows)
      # a RatioTapChangerTable point with r, x, g, b deviations: the SV is
      # computed with them on the series impedance and the magnetizing
      # admittance, the branch carries the factors of the step
      table = (r = 4.0, x = 10.0, g = -20.0, b = 15.0)
      rest = importCGMES(path = _cgmes_synth_tapside_dir(:two_winding; table = table), name = "tap2w_table")
      brt = rest.net.branchVec[last(rest.element_by_mrid["_pt"])]
      @test brt.x_pu ≈ br.x_pu * 1.1 && brt.r_pu ≈ br.r_pu * 1.04
      @test brt.b_pu ≈ br.b_pu * 1.15 && brt.g_pu ≈ br.g_pu * 0.8
      rows = flow_rows(rest)
      @test length(rows) == 2
      @test all(r -> abs(r.dp) < 1e-6 && abs(r.dq) < 1e-6, rows)
      # under a voltage controller the deviations belong to the step: the
      # import starts at the SSH step, the controller moves the tap and the
      # branch carries the table point of the step it reached (deviations
      # of the fixture: those of `table` times step / -4)
      ctl = importCGMES(path = _cgmes_synth_tapside_dir(:two_winding; table = table, control_kv = 21.2), name = "tap2w_table_ctl", tap_control = true)
      brc = ctl.net.branchVec[last(ctl.element_by_mrid["_pt"])]
      @test brc.x_pu ≈ brt.x_pu && brc.g_pu ≈ brt.g_pu
      cfg = SparlectraConfig(powerflow = PowerFlowConfig(max_iter = 30, tol = 1e-8), output = OutputConfig(logfile_results = :off), control = ControlConfig())
      run_sparlectra(net = ctl.net, config = cfg)
      # the end-2 tap sits at the from side: ratio = rated ratio * (1 + 1.5 % * step)
      s = round(Int, (brc.tap_ratio * (115.0 / 21.0) * (20.0 / 110.0) - 1.0) / 0.015)
      @test s != -4
      w = s / -4
      @test brc.x_pu ≈ br.x_pu * (1 + 0.10 * w) && brc.r_pu ≈ br.r_pu * (1 + 0.04 * w)
      @test brc.b_pu ≈ br.b_pu * (1 + 0.15 * w) && brc.g_pu ≈ br.g_pu * (1 - 0.20 * w)
    end)() end

    # The import stores per transformer branch the tap changers of the state
    # estimation (se_tap_changers, origin :cgmes). The truth is a
    # re-import of the same delivery with one changer at another step (SSH
    # step and SvTapStep set in a copy of the store): at the import step and
    # at two other steps the stored set gives the ratio and shift of the
    # re-imported branch (1e-12) and, where the changer carries a per-step
    # impedance, its r, x, g and b relative to the import step. Deliveries:
    # the own export of a net with typed asymmetrical (step 16 of 1..25,
    # neutral 13, 1.25 %, psi 90), symmetrical (with an xMin/xMax curve) and
    # tabular phase models, a RatioTapChanger, a transformer without changer,
    # a fixed shift (the export's single-step PhaseTapChangerLinear) and a
    # transformer tapped on both ends; the synthetic 2W delivery tapped on
    # end 2 with a RatioTapChangerTable carrying r, x, g, b on every point
    # (the writer has no RTC table); the synthetic 3W delivery with two
    # tapped legs.
    @testset "tap changer sets of the tap estimation" begin (function ()
      S = Sparlectra
      # a re-import of `store` with the tap changer `mrid` at step `k`
      function at_step(store, mrid, k)
        st = deepcopy(store)
        st.objects[mrid].attrs[:step] = string(k)
        for o in objectsOf(st, :SvTapStep)
          get(o.refs, :TapChanger, "") == mrid && (o.attrs[:position] = string(k))
        end
        return importCGMES(path = st, name = "at_step", require_boundary = false)
      end
      # set of branch `k` of `res` against the branch `k` of a re-import with
      # changer `c` at `steps` (absolute): ratio, shift, impedance factors
      function agrees(res, k, c, steps)
        set = res.net.branchVec[k].se_tap_changers
        ch = set.changers[c]
        s0 = Float64[x.step - x.neutral for x in set.changers]
        F0 = S._tap_impedance_factors(set, s0)
        b0 = res.net.branchVec[k]
        ok = true
        for st in steps
          bk = at_step(res.store, ch.mrid, st).net.branchVec[k]
          s = copy(s0)
          s[c] = st - ch.neutral
          R, P = S._tap_set_ratio_shift(set, s)
          F = S._tap_impedance_factors(set, s)
          ok &= isapprox(R, bk.tap_ratio; rtol = 1e-12) && abs(P - bk.phase_shift_deg) <= 1e-12
          ok &= isapprox(b0.r_pu * F[1] / F0[1], bk.r_pu; rtol = 1e-12) && isapprox(b0.x_pu * F[2] / F0[2], bk.x_pu; rtol = 1e-12)
          ok &= isapprox(b0.g_from_pu * F[3] / F0[3], bk.g_from_pu; rtol = 1e-12, atol = 1e-15) && isapprox(b0.b_from_pu * F[4] / F0[4], bk.b_from_pu; rtol = 1e-12, atol = 1e-15)
        end
        ok || println("      tap changer set of branch ", k, ", changer ", ch.name, ": no match with the re-import at steps ", steps)
        return ok
      end
      fields(ch) = (ch.kind, ch.role, ch.side, ch.du, ch.alpha_deg, ch.neutral, ch.low, ch.high, ch.step, ch.impedance)

      # --- the own export of typed models
      net = Net(name = "sets", baseMVA = 100.0)
      for b in ("A", "B", "C", "D", "E", "F", "G", "H")
        addBus!(net = net, busName = b, vn_kV = 110.0)
      end
      addProsumer!(net = net, busName = "A", type = "EXTERNALNETWORKINJECTION", referencePri = "A", vm_pu = 1.0, va_deg = 0.0)
      addProsumer!(net = net, busName = "H", type = "ENERGYCONSUMER", p = 10.0, q = 2.0)
      asym() = PhaseTapChangerModel(kind = :asymmetrical, step = 16, lowStep = 1, highStep = 25, neutralStep = 13, voltage_step_increment = 0.0125, winding_connection_angle_deg = 90.0, convention = :direct_regulating_vector)
      sym = PhaseTapChangerModel(kind = :symmetrical, step = 3, lowStep = -8, highStep = 8, neutralStep = 0, voltage_step_increment = 0.02, x_min = 0.05, x_max = 0.06, convention = :direct_regulating_vector)
      tab = PhaseTapChangerModel(kind = :tabular, step = 2, neutralStep = 0, convention = :direct_regulating_vector, table = [TapTablePoint(step = k, ratio = 1.0 + 0.001 * k^2, angle_deg = 2.0 * k) for k in -3:3])
      rtc = PowerTransformerTaps(Vn_kV = 110.0, step = 2, lowStep = -5, highStep = 5, neutralStep = 0, voltageIncrement_kV = 1.1)
      tr(f, t; shift = 0.0, phase_taps = nothing, taps = nothing) = addPIModelTrafo!(net = net, fromBus = f, toBus = t, r_pu = 0.002, x_pu = 0.05, b_pu = 0.0, status = 1, ratio = 1.0, shift_deg = shift, phase_taps = phase_taps, taps = taps)
      tr("A", "B"; phase_taps = asym())
      tr("B", "C"; phase_taps = sym)
      tr("C", "D"; phase_taps = tab)
      tr("D", "E"; taps = rtc)
      tr("E", "F")
      tr("F", "G"; shift = -3.0)
      tr("G", "H"; phase_taps = asym())
      # an RTC on end 2 next to the phase changer on end 1 (both ends tapped)
      net.trafos[end].side2.taps = PowerTransformerTaps(Vn_kV = 110.0, step = -1, lowStep = -4, highStep = 4, neutralStep = 0, voltageIncrement_kV = 1.65)
      dir = mktempdir()
      writeCGMESFiles(net; path = dir)
      res = importCGMES(path = dir, name = "sets", require_boundary = false)
      sets = [br.se_tap_changers for br in res.net.branchVec]
      @test all(st -> st !== nothing && st.origin === :cgmes, sets) && all(st -> st.ratio_n == 1.0, sets[1:6])
      # the export writes ratedU1 with the end-2 correction of step -1
      # (1 - 1.5 %) so that the re-import keeps the branch ratio 1
      @test isapprox(sets[7].ratio_n, 0.985; rtol = 1e-12)
      @test fields(only(sets[1].changers)) == (:asymmetrical, :phase, :from, 0.0125, 90.0, 13, 1, 25, 16.0, :none)
      @test fields(only(sets[2].changers))[1:9] == (:symmetrical, :phase, :from, 0.02, nothing, 0, -8, 8, 3.0) && only(sets[2].changers).impedance === :x_curve
      t3 = only(sets[3].changers)
      @test (t3.kind, t3.role, t3.side, t3.neutral, t3.step) == (:tabular, :phase, :from, 0, 2.0) && [(p.step, p.ratio, p.angle_deg) for p in t3.table] == [(k, 1.0 + 0.001 * k^2, 2.0 * k) for k in -3:3]
      @test fields(only(sets[4].changers)) == (:ratio, :ratio, :from, 0.01, nothing, 0, -5, 5, 2.0, :none)
      # no changer: an empty set, never released
      @test isempty(sets[5].changers) && S._tap_release_mode(res.net.branchVec[5]) === nothing
      # the export's single-step PhaseTapChangerLinear is the fixed shift
      @test isempty(sets[6].changers) && sets[6].shift_n_deg == -3.0 && res.net.branchVec[6].phase_shift_deg == -3.0
      @test count(m -> occursin("is the single-step phase changer of a CGMES export: its shift of -3.0 deg stays fixed", m), res.messages) == 1
      # tapped on both ends: the end-1 changer multiplies, the end-2 one divides
      @test [fields(c) for c in sets[7].changers] == [(:asymmetrical, :phase, :from, 0.0125, 90.0, 13, 1, 25, 16.0, :none), (:ratio, :ratio, :to, 0.015, nothing, 0, -4, 4, -1.0, :none)]
      @test S._tap_release_mode(res.net.branchVec[7]) === :both
      @test agrees(res, 1, 1, (16, 3, 25))
      @test agrees(res, 2, 1, (3, -8, 5))
      @test agrees(res, 3, 1, (2, -3, 1))
      @test agrees(res, 4, 1, (2, -5, 4))
      @test agrees(res, 7, 1, (16, 10, 20)) && agrees(res, 7, 2, (-1, -3, 2))
      # the pattern needs all five values and the suffix: renamed, or with
      # another normalStep, the changer is a linear phase changer again
      ptc = only(tc for tc in objectsOf(res.store, :PhaseTapChangerLinear))
      for (key, value) in ((:name, "SHIFTER"), (:normalStep, "0"))
        st = deepcopy(res.store)
        st.objects[ptc.mrid].attrs[key] = value
        net6 = importCGMES(path = st, name = "not_the_pattern", require_boundary = false).net
        set6 = net6.branchVec[6].se_tap_changers
        @test set6.shift_n_deg == 0.0 && [(c.kind, c.dphi_deg, c.low, c.high, c.step) for c in set6.changers] == [(:linear, -3.0, 0, 1, 1.0)]
        # a degree-form changer (alpha 0) is released by its kind: the
        # mass release takes it as a phase shifter
        @test S._tap_release_mode(net6.branchVec[6]) === :pst
        @test setTapEstimation!(deepcopy(net6); trafo = 6, mode = :pst).s0 == [1.0]
      end

      # --- a changer whose table the estimation cannot use (a point at a
      # step that is not whole), next to a usable one: the import applies
      # both, so the set keeps the left-out one as a fixed part and gives
      # the branch at the import steps without moving the usable changer
      net2 = Net(name = "sets_fixed", baseMVA = 100.0)
      addBus!(net = net2, busName = "A", vn_kV = 110.0)
      addBus!(net = net2, busName = "B", vn_kV = 110.0)
      addProsumer!(net = net2, busName = "A", type = "EXTERNALNETWORKINJECTION", referencePri = "A", vm_pu = 1.0, va_deg = 0.0)
      addProsumer!(net = net2, busName = "B", type = "ENERGYCONSUMER", p = 10.0, q = 2.0)
      addPIModelTrafo!(net = net2, fromBus = "A", toBus = "B", r_pu = 0.002, x_pu = 0.05, b_pu = 0.0, status = 1, ratio = 1.0, shift_deg = 0.0, phase_taps = tab)
      net2.trafos[end].side2.taps = PowerTransformerTaps(Vn_kV = 110.0, step = -1, lowStep = -4, highStep = 4, neutralStep = 0, voltageIncrement_kV = 1.65)
      dir2 = mktempdir()
      writeCGMESFiles(net2; path = dir2)
      st2 = importCGMES(path = dir2, name = "sets_fixed", require_boundary = false).store
      pt3 = only(pt for pt in objectsOf(st2, :PhaseTapChangerTablePoint) if num(pt, :step) == -3.0)
      pt3.attrs[:step] = "-3.4"
      rf = importCGMES(path = st2, name = "sets_fixed", require_boundary = false)
      bf = rf.net.branchVec[1]
      sf = bf.se_tap_changers
      @test [(c.kind, c.side, c.du, c.step) for c in sf.changers] == [(:ratio, :to, 0.015, -1.0)]
      # ratio_n: the structural ratio 0.985 of the export (see set 7) times
      # the table ratio 1.004 at step 2
      @test isapprox(sf.ratio_n, 0.985 * 1.004; rtol = 1e-12) && sf.shift_n_deg == 4.0
      R, P = S._tap_set_ratio_shift(sf, [-1.0])
      @test isapprox(R, bf.tap_ratio; rtol = 1e-12) && abs(P - bf.phase_shift_deg) <= 1e-12
      @test S._tap_s_from_branch(sf, bf) == [-1.0]
      @test count(m -> occursin("is not estimable, its PhaseTapChangerTable has steps that are not whole and strictly ascending: its ratio 1.004 and shift 4.0 deg at step 2.0 stay fixed", m), rf.messages) == 1

      # --- RatioTapChangerTable with r, x, g, b on every point, end 2 only
      # (the branch is turned, the changer sits at its from side)
      rt = importCGMES(path = _cgmes_synth_tapside_dir(:two_winding; table = (r = 4.0, x = 10.0, g = -20.0, b = 15.0), control_kv = 21.2), name = "tap2w_table")
      kt = last(rt.element_by_mrid["_pt"])
      ct = only(rt.net.branchVec[kt].se_tap_changers.changers)
      @test (ct.kind, ct.role, ct.side, ct.neutral, ct.low, ct.high, ct.step, ct.outside_du, ct.impedance, length(ct.table)) == (:tabular, :ratio, :from, 0, -10, 10, -4.0, 0.015, :table, 21)
      @test agrees(rt, kt, 1, (-4, 0, 3))
      # a point at a step that is not whole, away from the import step: the
      # import still applies the table, the set keeps it as the fixed part
      st_rt = deepcopy(rt.store)
      pt7 = only(pt for pt in objectsOf(st_rt, :RatioTapChangerTablePoint) if num(pt, :step) == 7.0)
      pt7.attrs[:step] = "7.5"
      rt2 = importCGMES(path = st_rt, name = "tap2w_table_fixed")
      b2 = rt2.net.branchVec[kt]
      @test isempty(b2.se_tap_changers.changers) && S._tap_release_mode(b2) === nothing
      @test (b2.tap_ratio, b2.phase_shift_deg) == (rt.net.branchVec[kt].tap_ratio, rt.net.branchVec[kt].phase_shift_deg)
      @test isapprox(b2.se_tap_changers.ratio_n, b2.tap_ratio; rtol = 1e-12) && b2.se_tap_changers.shift_n_deg == 0.0

      # --- 3W: two tapped legs (from the side bus to the star), one leg
      # without changer
      r3 = importCGMES(path = _cgmes_synth3w_dir(end1_tap_step = 13), name = "synth3w_sets")
      legs = [last(r3.element_by_mrid[e]) for e in ("_t3w_e1", "_t3w_e2", "_t3w_e3")]
      leg_sets = [r3.net.branchVec[k].se_tap_changers for k in legs]
      @test [(c.kind, c.side, c.du, c.neutral, c.step) for c in leg_sets[1].changers] == [(:ratio, :from, 0.01, 11, 13.0)]
      @test [(c.kind, c.side, c.du, c.neutral, c.step) for c in leg_sets[2].changers] == [(:ratio, :from, 0.0125, 11, 11.0)]
      @test isempty(leg_sets[3].changers)
      @test agrees(r3, legs[1], 1, (13, 9, 15)) && agrees(r3, legs[2], 1, (11, 8, 14))
    end)() end

    # sp_casePST in CGMES and SCF: the SCF file carries the nameplate
    # Delta-u phase shifter of branch 17 (alpha 90, 1 % per step, -10..10,
    # angle base -5.729578 deg). The CGMES export has no
    # class for it and writes the fixed angle with the single-step
    # PhaseTapChangerLinear 2WT_230_PTC, which the import recognises: in
    # CGMES the transformer has no regulator, in SCF it has one. The
    # RatioTapChanger of the machine transformer is a changer in both and
    # stays out of the mass release (machine transformer). The lossy export
    # is a follow-up outside the tap estimation.
    @testset "sp_casePST: the CGMES and the SCF tap changer sets differ" begin (function ()
      S = Sparlectra
      scf = importSCF(joinpath(dirname(@__DIR__), "data", "scf", "sp_casePST.scf.json"))
      cg = importCGMES(path = cgmes_fixture_dir("sp_casePST"), name = "sp_casePST")
      pst(net) = only(k for (k, br) in enumerate(net.branchVec) if br.ratio != 0.0 && br.phase_shift_deg != 0.0)
      ks, kc = pst(scf), pst(cg.net)
      ss = S._effective_tap_changers(scf.branchVec[ks])
      @test [(c.kind, c.side, c.du, c.alpha_deg, c.low, c.high) for c in ss.changers] == [(:asymmetrical, :to, 0.01, 90.0, -10, 10)]
      @test S._tap_release_mode(scf.branchVec[ks]) === :pst
      cs = cg.net.branchVec[kc].se_tap_changers
      @test cs.origin === :cgmes && isempty(cs.changers) && cs.shift_n_deg == -5.729578
      @test S._tap_release_mode(cg.net.branchVec[kc]) === nothing
      @test scf.branchVec[ks].phase_shift_deg == cg.net.branchVec[kc].phase_shift_deg == -5.729578
      for net in (scf, cg.net)
        km = only(k for (k, br) in enumerate(net.branchVec) if br.ratio != 0.0 && k != pst(net))
        ch = only(S._effective_tap_changers(net.branchVec[km]).changers)
        @test (ch.kind, ch.du, ch.low, ch.high) == (:ratio, 0.025, -2, 2) && S._is_machine_transformer(net, km)
      end
    end)() end

    # A converged solution on the low-voltage branch stays converged and
    # carries the plausibility warning (power_flow.plausibility.vmin_warn,
    # default 0.7 pu) in the result diagnostics, the log, the service
    # metadata and run.log. The heavy 3W delivery has two solutions, MV bus
    # 0.954 pu and 0.456 pu. Before the star-point start of the importer and
    # the tap-side fix, the HV-side star start led Newton from the operating
    # SV to the low-voltage solution; on the corrected model it no longer
    # does, so the low-voltage solution is reached from its own SV
    # (`sv = :low`), a synthetic trigger.
    @testset "converged low-voltage solution gets the plausibility warning" begin (function ()
      cfg = SparlectraConfig(powerflow = PowerFlowConfig(max_iter = 30, tol = 1e-8), output = OutputConfig(logfile_results = :off, console_summary = false))
      res = importCGMES(path = _cgmes_synth3w_heavy_dir(sv = :low), name = "synth3w_heavy_low")
      r = @test_logs (:warn, r"lowest voltage 0\.4560 pu at bus MV_BUS \(110\.0 kV\), 2 energized bus\(es\) below power_flow\.plausibility\.vmin_warn = 0\.7 pu") match_mode = :any run_sparlectra(net = res.net, config = cfg)
      # a warning, no abort
      @test r.numerical_converged && r.final_converged && r.solution_available
      d = r.diagnostics
      @test d.plausibility_status === :low_voltage
      @test d.plausibility_min_vm_bus == "MV_BUS"
      @test isapprox(d.plausibility_min_vm_pu, 0.4560; atol = 1e-4)
      @test d.plausibility_low_vm_count == 2
      @test occursin("MV_BUS", d.plausibility_warning)
      # the threshold comes from the configuration: at 0.4 the same
      # solution passes, 0 switches the check off
      cfg04 = SparlectraConfig(Dict("power_flow" => Dict("max_iter" => 30, "tol" => 1e-8, "plausibility" => Dict("vmin_warn" => 0.4)), "output" => Dict("logfile_results" => "off", "console_summary" => false)))
      @test cfg04.powerflow.plausibility.vmin_warn == 0.4
      r04 = @test_logs min_level = Logging.Warn run_sparlectra(net = res.net, config = cfg04)
      @test r04.diagnostics.plausibility_status === :ok
      @test isapprox(r04.diagnostics.plausibility_min_vm_pu, 0.4560; atol = 1e-4)
      @test isempty(r04.diagnostics.plausibility_warning)
      cfg0 = SparlectraConfig(powerflow = PowerFlowConfig(max_iter = 30, tol = 1e-8, plausibility = Sparlectra.PowerFlowPlausibilityConfig(vmin_warn = 0.0)), output = OutputConfig(logfile_results = :off, console_summary = false))
      @test (@test_logs min_level = Logging.Warn run_sparlectra(net = res.net, config = cfg0)).diagnostics.plausibility_status === :not_checked
      # the same injections from a flat start reach the operating solution:
      # checked, no warning
      flat = importCGMES(path = _cgmes_synth3w_heavy_dir(), name = "synth3w_heavy_flat")
      cfg_flat = SparlectraConfig(powerflow = PowerFlowConfig(max_iter = 30, tol = 1e-8, start_mode = StartModeConfig(flatstart = true)), output = OutputConfig(logfile_results = :off, console_summary = false))
      rf = @test_logs min_level = Logging.Warn run_sparlectra(net = flat.net, config = cfg_flat)
      @test rf.diagnostics.plausibility_status === :ok
      @test isapprox(rf.diagnostics.plausibility_min_vm_pu, 0.9536; atol = 1e-4)
      # the service run on the delivery with the low-voltage SV: run.log and
      # the result metadata
      out = mktempdir()
      rs = run_sparlectra_api(casefile = _pack_cgmes_dir_zip(_cgmes_synth3w_heavy_dir(sv = :low), "synth3w_heavy_low"), output_dir = out)
      @test rs.status == :succeeded
      @test rs.metadata["plausibility_status"] == "low_voltage"
      @test rs.metadata["plausibility_min_vm_bus"] == "MV_BUS"
      @test isapprox(rs.metadata["plausibility_min_vm_pu"], 0.4560; atol = 1e-4)
      @test occursin("power flow plausibility: lowest voltage 0.4560 pu at bus MV_BUS", read(joinpath(out, "run.log"), String))
    end)() end

    @testset "isolated bus as the smallest bus of a fed switch group" begin (function ()
      # C keeps only its switch to the loaded bus B and is the smallest bus
      # of the group: the import types it Isolated, and the solver starts
      # the merged node {C, B} from C. C starts at the voltage of B, so the
      # solve ends on the high-voltage solution, not on the low-voltage
      # branch a 0 V start reaches
      net, _, s2 = isolated_first_cgmes_snapshots()
      rk = redirect_stdout(() -> importCGMES(path = s2, name = "isolated_first"), devnull)
      n = rk.net
      bB, bC = n.busDict["B"], n.busDict["C"]
      @test bC == 1 && Sparlectra._active_link_representative_map(n)[bB] == bC
      @test Sparlectra.isIsolated(n.nodeVec[bC]) && bC in n.isoNodes
      @test (n.nodeVec[bC]._vm_pu, n.nodeVec[bC]._va_deg) == (n.nodeVec[bB]._vm_pu, n.nodeVec[bB]._va_deg)
      redirect_stdout(() -> runpf!(n, 30, 1e-8, 0), devnull)
      ref = isolated_first_reference(net)
      @test n.nodeVec[bB]._vm_pu ≈ ref._vm_pu atol = 1e-6
      @test n.nodeVec[bB]._va_deg ≈ ref._va_deg atol = 1e-4
    end)() end

    # The service legs on in-memory deliveries: the synthetic 3W delivery
    # packed as a ZIP is the healthy case, the synthetic EQ with an absent
    # prerequisite the broken one, and a MATPOWER case covers the non-CGMES
    # rejections (the healthy runs on a checked-in delivery are above).
    # --- self-built deliveries (data/cgmes_demo) --------------------------
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

    # The importer's exit criterion on its own exporter's output: the SV
    # profile and the source net's solution are reproduced to noise. Both
    # start points are checked, the SV state (the solve must not move) and
    # a flat start (the model itself must reach the same state).
    @testset "fixture SV profile and source solution are reproduced" begin (function ()
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
        # the service layer takes case FILES, so the zip packed from the
        # same four profiles must rebuild the same network
        rz = importCGMES(path = _pack_cgmes_fixture_zip(case), name = string(case, "_zip"))
        @test Set(keys(rz.net.busDict)) == Set(keys(res.net.busDict))
        @test length(rz.net.branchVec) == length(res.net.branchVec)
        @test rz.slack_bus == res.slack_bus
      end
      # an SvVoltage with v = 0 counts as no SV, as in compareWithSV: the
      # bus does not start at zero voltage but at the start of the neighbour
      # over the branch or closed link of the smallest |z|, and is listed in
      # no_sv_buses
      plain = importCGMES(path = cgmes_fixture_dir("sp_case14"), name = "v0_plain")
      v0_bus, v0_tn = first(sort([(bus, tn) for (tn, bus) in plain.topo.bus_name if plain.net.nodeVec[plain.net.busDict[bus]]._nodeType == Sparlectra.PQ]))
      v0_store = loadCGMES(cgmes_fixture_dir("sp_case14"))
      only(sv for sv in objectsOf(v0_store, :SvVoltage) if sv.refs[:TopologicalNode] == v0_tn).attrs[:v] = "0"
      v0 = importCGMES(path = v0_store, name = "v0")
      v0_i = v0.net.busDict[v0_bus]
      v0_node = v0.net.nodeVec[v0_i]
      v0_nbrs = vcat([(abs(complex(br.r_pu, br.x_pu)), br.fromBus == v0_i ? br.toBus : br.fromBus) for br in v0.net.branchVec if br.status == 1 && v0_i in (br.fromBus, br.toBus)], [(0.0, l.fromBus == v0_i ? l.toBus : l.fromBus) for l in v0.net.linkVec if l.status == 1 && v0_i in (l.fromBus, l.toBus)])
      v0_src = v0.net.nodeVec[last(minimum(v0_nbrs))]
      @test (v0_node._vm_pu, v0_node._va_deg) == (v0_src._vm_pu, v0_src._va_deg)
      @test v0.no_sv_buses == [v0_bus]
      @test count(m -> startswith(m, "notice: 1 bus(es) without a usable SvVoltage start at the voltage of a neighbour"), v0.messages) == 1
      @test !any(r -> r.bus == v0_bus, compareWithSV(v0).rows)
      # a leaf bus with SvVoltage v = 0 behind an EquivalentBranch of 0.02 pu
      # with a P = Q = 0 EquivalentInjection, on a copy whose SV angles are
      # all shifted by -40 deg: started at 1.0 pu / 0 deg, Newton ends with
      # the leaf at 0 V, a short circuit that pulls its neighbour down to
      # 0.3 pu; from the neighbour's start it solves at the neighbour voltage
      leaf_store = loadCGMES(cgmes_fixture_dir("sp_case14"))
      for sv in objectsOf(leaf_store, :SvVoltage)
        sv.attrs[:angle] = string(num(sv, :angle, 0.0) - 40.0)
      end
      leaf_host = first(sort(collect(objectsOf(leaf_store, :SvVoltage)); by = sv -> num(sv, :angle, 0.0))).refs[:TopologicalNode]
      leaf_add!(mrid, class, attrs, refs = Dict{Symbol,String}()) = (leaf_store.objects[mrid] = Sparlectra.CGMESImporter.CIMObject(mrid, class, Dict{Symbol,String}(attrs), Dict{Symbol,String}(refs), "leaf_test"); push!(get!(Vector{String}, leaf_store.byclass, class), mrid))
      leaf_term!(eq, tn, seq) = leaf_add!(string(eq, "_t", seq), :Terminal, Dict(:sequenceNumber => string(seq), :connected => "true"), Dict(:ConductingEquipment => eq, :TopologicalNode => tn))
      leaf_add!("_leaf_tn", :TopologicalNode, Dict(:name => "LEAF"), Dict(:BaseVoltage => leaf_store.objects[leaf_host].refs[:BaseVoltage]))
      leaf_add!("_leaf_sv", :SvVoltage, Dict(:v => "0", :angle => "0"), Dict(:TopologicalNode => "_leaf_tn"))
      leaf_add!("_leaf_eb", :EquivalentBranch, Dict(:name => "LEAF_EB", :r => "0.3", :x => "2.4"))
      leaf_term!("_leaf_eb", leaf_host, 1)
      leaf_term!("_leaf_eb", "_leaf_tn", 2)
      leaf_add!("_leaf_ei", :EquivalentInjection, Dict(:name => "LEAF_EI", :p => "0", :q => "0"))
      leaf_term!("_leaf_ei", "_leaf_tn", 1)
      leaf = importCGMES(path = leaf_store, name = "leaf")
      leaf_bus = leaf.topo.bus_name["_leaf_tn"]
      @test leaf.no_sv_buses == [leaf_bus]
      leaf_node = leaf.net.nodeVec[leaf.net.busDict[leaf_bus]]
      host_node = leaf.net.nodeVec[leaf.net.busDict[leaf.topo.bus_name[leaf_host]]]
      @test host_node._va_deg < -40.0
      @test (leaf_node._vm_pu, leaf_node._va_deg) == (host_node._vm_pu, host_node._va_deg)
      solve_fixture!(leaf.net)
      @test leaf_node._vm_pu > 0.9
      @test isapprox(leaf_node._vm_pu, host_node._vm_pu; atol = 1e-6)
      @test isapprox(leaf_node._va_deg, host_node._va_deg; atol = 1e-4)
      # a listed star point is left alone, and a listed bus reachable only
      # through it keeps its start (it used to look up bus index 0)
      star_net = Net(name = "star_no_sv", baseMVA = 100.0)
      for b in ("H", "S", "M")
        addBus!(net = star_net, busName = b, vn_kV = 110.0)
      end
      addPIModelACLine!(net = star_net, fromBus = "H", toBus = "S", r_pu = 0.02, x_pu = 0.2, b_pu = 0.0, status = 1)
      addPIModelACLine!(net = star_net, fromBus = "S", toBus = "M", r_pu = 0.02, x_pu = 0.2, b_pu = 0.0, status = 1)
      Sparlectra.setVmVa!(node = star_net.nodeVec[star_net.busDict["H"]], vm_pu = 1.01, va_deg = -20.0)
      star_ctx = Sparlectra.CGMESImporter._MapCtx()
      push!(star_ctx.star_buses, star_net.busDict["S"])
      star_start = [(nd._vm_pu, nd._va_deg) for nd in star_net.nodeVec]
      @test Sparlectra.CGMESImporter._startNoSvBusesFromNeighbours!(star_net, star_ctx, ["M", "S"]) == 0
      @test [(nd._vm_pu, nd._va_deg) for nd in star_net.nodeVec] == star_start
    end)() end

    # sp_case118 carries 54 synchronous machines, among them condensers
    # (P = 0, voltage regulated) and asymmetric minQ/maxQ pairs. The
    # condensers arrive as regulated zero-P units. The Q limits arrive as
    # the Stage-1 hull of both sign readings (the ENTSO-E sets are
    # inconsistent about the convention), so an asymmetric pair widens to
    # its symmetric envelope; this pins that contract rather than limit
    # fidelity.
    @testset "sp_case118 machines: synchronous condensers and the Q-limit hull" begin (function ()
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
    end)() end

    # Fixed-reference self-check on a delivery: the SV voltages reach the
    # solver verbatim (a base config with flatstart: true must not wipe
    # them) and the residual/attribution artifacts are written. The fixture
    # SV is a solution, so the start residual is noise; a flat start would
    # leave MW-sized residuals on the loaded buses.
    @testset "CGMES fixed-reference self-check (SV start, artifacts)" begin (function ()
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
      residuals = run_csv(joinpath(out, "self_check_residuals.csv"))
      @test residuals.header == ["bus_id", "bus_name", "vn_kV", "bus_type", "vm_pu_start", "va_deg_start", "p_residual", "q_residual", "has_sv", "n_transformer_terminals", "n_shunts"]
      @test length(residuals.rows) == 9
      @test all(row[9] == "true" for row in residuals.rows)
      @test any(residuals.number(row[10]) > 0 for row in residuals.rows)
      @test occursin("No-SV buses: 0", cgmes_report_plain(out))
      # a delivery without DY profile gets no DY section
      @test !occursin("Dynamic data (DY profile)", cgmes_report_plain(out))
      # the tap changer of the import notice (SvTapStep four steps off the
      # SSH step) is listed in self_check.log with both steps
      off = _sp_case14_off_sv_tap_step()
      tout = mktempdir()
      rtap = run_fixed_reference_self_check(casefile = off.dir, config_file = cfg, output_dir = tout, case_format = :cgmes)
      fit = only(rtap.metadata["cgmes_sv_tap_step_fits"])
      @test (fit.sv_step, fit.best_step) == (off.sv_step, off.ssh_step)
      tap_summary = read(joinpath(tout, "self_check.log"), String)
      @test occursin("cgmes_sv_tap_step_fits: 1 ", tap_summary)
      @test occursin(": SvTapStep $(off.sv_step) mismatch ", tap_summary) && occursin(", best step $(off.ssh_step) mismatch ", tap_summary)
      @test occursin("notice: 1 tap changer(s) whose SvTapStep fits", cgmes_report_plain(tout))
    end)() end

    @testset "import_analysis_mode service run (synthetic deliveries)" begin (function ()
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
    end)() end

    # The Web UI "Short circuit" button's service path: negative cases must
    # fail with explicit reasons, never with empty tables (the run on a
    # checked-in delivery is above).
    @testset "short_circuit_mode service run rejects non-CGMES cases and mode mixing" begin (function ()
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
    end)() end

    # cgmes_import.start_values is a CGMES-only key: a MATPOWER run ignores
    # it completely (no decision line, no SV artifacts, no metadata keys).
    @testset "cgmes_import.start_values has no effect on a MATPOWER run" begin (function ()
      mp_out = mktempdir()
      mp_cfg = joinpath(mp_out, "c.yaml")
      write(mp_cfg, "config_version: 1\ncgmes_import:\n  start_values: sv\n")
      mp = run_sparlectra_api(casefile = abspath(joinpath(dirname(@__DIR__), "data", "mpower", "warmup_casePST.m")), config_file = mp_cfg, output_dir = mp_out)
      @test mp.status == :succeeded
      @test !isfile(joinpath(mp_out, "sv_compare.csv"))
      @test !haskey(mp.metadata, "cgmes_start_values")
      @test !haskey(mp.metadata, "cgmes_sv_compare_status")
      @test !occursin("CGMES start values", read(joinpath(mp_out, "run.log"), String))
    end)() end

    # The voltage-setpoint plausibility band is validated before any file is
    # read: an inverted band is an argument error, not a silent no-op.
    @testset "inverted vset band is rejected" begin (function ()
      @test_throws ArgumentError importCGMES(path = _cgmes_synth3w_dir(), name = "vset_bad", vset_min_pu = 1.2, vset_max_pu = 0.8)
    end)() end

    # ENTSO-E SmallGrid (118 buses, boundary set, node-breaker with TP).
    # The importer names elements after its own bus indices, not after
    # IdentifiedObject.name, so the identity check is topological: the
    # addressed element sits on the TopologicalNode buses of the object's
    # terminals and carries its SSH power. Links, 3W legs, tap controllers
    # and the self-loop charging shunt are checked on the synthetic and
    # self-built deliveries above; SmallGrid has none of them.
    @testset "SmallGrid: element_by_mrid and import from a loaded store" begin (function ()
      CI = Sparlectra.CGMESImporter
      files = sort!(filter(f -> endswith(f, ".zip"), readdir(joinpath(@__DIR__, "fixtures", "cgmes_smallgrid"); join = true)))
      res = importCGMES(path = files, name = "smallgrid")
      net, store, topo, addr = res.net, res.store, res.topo, res.element_by_mrid
      @test isempty(res.skipped_equipment)
      names = Sparlectra._bus_name_by_idx(net)
      tbus(obj, seq) = CI.busOfEquipment(topo, obj, seq).bus
      ebus(e) = topo.bus_name[CI.tnOfTerminal(store.objects[e.refs[:Terminal]])]
      # every created object is addressed, at the element on its buses
      for line in objectsOf(store, :ACLineSegment)
        kind, k = addr[line.mrid]
        br = net.branchVec[k]
        @test kind == :branch
        @test (names[br.fromBus], names[br.toBus]) == (tbus(line, 1), tbus(line, 2))
      end
      # a transformer with tap changers only at end 2 runs end 2 -> end 1
      # (SmallGrid has one)
      n_turned = 0
      for pt in objectsOf(store, :PowerTransformer)
        kind, k = addr[pt.mrid]
        ends = CI._endsOfTransformer(store, pt)
        br = net.branchVec[k]
        @test kind == :transformer
        turned = isempty(CI._tapChangersOfEnd(store, ends[1].mrid)) && !isempty(CI._tapChangersOfEnd(store, ends[2].mrid))
        n_turned += turned
        @test (names[br.fromBus], names[br.toBus]) == (turned ? (ebus(ends[2]), ebus(ends[1])) : (ebus(ends[1]), ebus(ends[2])))
      end
      @test n_turned == 1
      # the magnetizing admittance of a 2W transformer sits half on each
      # terminal (0.41.0), as in the SV of SmallGrid: the unit at the end of
      # its one transformer with an admittance matches at the SV state; the
      # per-end placement of 0.20.0 missed it by 2.29 MVAr
      pinevlle = only(r for r in compareWithSV(res).flows.rows if r.kind == :units && occursin("Pinevlle", r.name))
      @test abs(pinevlle.dq) < 0.01
      for (cls, kind, sign) in ((:EnergyConsumer, :load, 1.0), (:SynchronousMachine, :generator, -1.0), (:EquivalentInjection, :external_grid, -1.0))
        for obj in objectsOf(store, cls)
          obj_kind, k = addr[obj.mrid]
          ps = net.prosumpsVec[k]
          @test obj_kind == kind
          @test names[Sparlectra.getPosumerBusIndex(ps)] == tbus(obj, 1)
          @test ps.pVal == sign * num(obj, :p, 0.0)
        end
      end
      # an EquivalentInjection takes its Q limits like a machine: the
      # ReactiveCapabilityCurve at its SSH p, else the minQ/maxQ hull, and
      # none from a missing or degenerate pair. SmallGrid states none; a copy
      # of the store gets a pair on the first, a curve and a pair it outranks
      # on the second and the pair 0/0 on the third. Its minP/maxP arrive as
      # delivered, in generator sign (none from 0/0).
      ei_q(r, o) = (u -> (u.minQ, u.maxQ))(r.net.prosumpsVec[last(r.element_by_mrid[o.mrid])])
      ei_p(r, o) = (u -> (u.minP, u.maxP))(r.net.prosumpsVec[last(r.element_by_mrid[o.mrid])])
      @test all(o -> ei_q(res, o) == (nothing, nothing), objectsOf(store, :EquivalentInjection))
      @test all(o -> ei_p(res, o) == (nothing, nothing), objectsOf(store, :EquivalentInjection))
      local eq_store = deepcopy(store)
      local eis = collect(objectsOf(eq_store, :EquivalentInjection))
      merge!(eis[1].attrs, Dict(:minQ => "-40", :maxQ => "30", :minP => "0", :maxP => "300"))
      merge!(eis[2].attrs, Dict(:minQ => "-90", :maxQ => "90", :minP => "-50", :maxP => "80"))
      eis[2].refs[:ReactiveCapabilityCurve] = "_ei_rcc"
      merge!(eis[3].attrs, Dict(:minQ => "0", :maxQ => "0", :minP => "0", :maxP => "0"))
      eq_store.objects["_ei_rcc"] = CI.CIMObject("_ei_rcc", :ReactiveCapabilityCurve, Dict{Symbol,String}(), Dict{Symbol,String}(), "ei_q_test")
      push!(get!(Vector{String}, eq_store.byclass, :ReactiveCapabilityCurve), "_ei_rcc")
      for (k, x) in enumerate(("-1000", "1000"))
        eq_store.objects["_ei_rcc_$(k)"] = CI.CIMObject("_ei_rcc_$(k)", :CurveData, Dict{Symbol,String}(:xvalue => x, :y1value => "-25", :y2value => "25"), Dict{Symbol,String}(:Curve => "_ei_rcc"), "ei_q_test")
        push!(get!(Vector{String}, eq_store.byclass, :CurveData), "_ei_rcc_$(k)")
      end
      local eq_res = importCGMES(path = eq_store, name = "smallgrid_ei_q")
      @test [ei_q(eq_res, o) for o in eis] == [(-40.0, 40.0), (-25.0, 25.0), (nothing, nothing)]
      @test [ei_p(eq_res, o) for o in eis] == [(0.0, 300.0), (-50.0, 80.0), (nothing, nothing)]
      @test count(m -> occursin("EquivalentInjection $(something(str(eis[2], :name), eis[2].mrid)): Q limits [-25.0, 25.0] MVAr from its ReactiveCapabilityCurve", m), eq_res.messages) == 1
      for sh in objectsOf(store, :LinearShuntCompensator)
        kind, k = addr[sh.mrid]
        @test kind == :shunt
        @test names[net.shuntVec[k].busIdx] == tbus(sh, 1)
      end
      # the closed switches all sit inside one TopologicalNode: no link
      @test isempty(net.linkVec)
      @test !any(v -> first(v) == :link, values(addr))
      # indices are unique per kind, and every branch, prosumer and shunt of
      # this delivery is addressed
      by_kind = Dict{Symbol,Vector{Int}}()
      for (kind, k) in values(addr)
        push!(get!(Vector{Int}, by_kind, kind), k)
      end
      @test all(allunique, values(by_kind))
      @test sort(vcat(by_kind[:branch], by_kind[:transformer])) == collect(eachindex(net.branchVec))
      @test sort(vcat(by_kind[:load], by_kind[:generator], by_kind[:external_grid])) == collect(eachindex(net.prosumpsVec))
      @test sort(by_kind[:shunt]) == collect(eachindex(net.shuntVec))
      # a loaded store builds the same network as the path and stays as it was
      fingerprint(s) = (Dict(m => (o.class, copy(o.attrs), copy(o.refs)) for (m, o) in s.objects), deepcopy(s.byclass))
      loaded = loadCGMES(files)
      before = fingerprint(loaded)
      from_store = importCGMES(path = loaded, name = "smallgrid")
      @test isempty(scf_roundtrip_field_diffs(net, from_store.net))
      @test from_store.element_by_mrid == addr
      @test from_store.store === loaded
      @test fingerprint(loaded) == before
      # without TP and SV the topology processor derives the buses and
      # writes them into the store, so the import has to work on a copy
      nb_store = loadCGMES(files; profile_filter = setdiff(CI.IMPORT_PROFILE_TAGS, Set([:TP, :SV])))
      nb_before = fingerprint(nb_store)
      nb = importCGMES(path = nb_store, name = "smallgrid_nb")
      @test any(m -> startswith(m, "topology processor:"), nb.messages)
      @test length(nb.net.nodeVec) == length(net.nodeVec)
      @test fingerprint(nb_store) == nb_before
      # the configured import takes ignore_connected from the configuration
      # unless the caller passes it
      cfg = SparlectraConfig()
      @test !any(m -> occursin("ignore_connected = true", m), importCGMES(cfg; path = loaded, name = "sg_cfg").messages)
      @test any(m -> occursin("ignore_connected = true", m), importCGMES(cfg; path = loaded, name = "sg_cfg_ic", ignore_connected = true).messages)

      # P limits of 9999 MW and above are placeholders (for example
      # +-99999.9 MW on units without a limit) and arrive as no limit. In a copy
      # DeerCrk SM (maxOperatingP 40 MW) gets +-99999.9 on its unit and its
      # Q pair, and the first EquivalentInjection a P pair with one
      # placeholder side. Taken as a limit, the placeholder made DeerCrk the
      # largest unit of the net: the replacement slack of ensureSlack! and
      # about 0.95 of the pmax weights of the distributed slack.
      @testset "placeholder P limits" begin (function ()
        local ph_store = deepcopy(loaded)
        local deer = only(sm for sm in objectsOf(ph_store, :SynchronousMachine) if str(sm, :name) == "DeerCrk SM")
        merge!(ph_store.objects[deer.refs[:GeneratingUnit]].attrs, Dict(:minOperatingP => "-99999.9", :maxOperatingP => "99999.9"))
        merge!(deer.attrs, Dict(:minQ => "-99999.9", :maxQ => "99999.9"))
        local ei1 = first(objectsOf(ph_store, :EquivalentInjection))
        merge!(ei1.attrs, Dict(:minP => "-99999.9", :maxP => "300"))
        local ph = importCGMES(path = ph_store, name = "sg_pmax_placeholder")
        local k_deer = last(ph.element_by_mrid[deer.mrid])
        local ps_deer = ph.net.prosumpsVec[k_deer]
        @test (ps_deer.minP, ps_deer.maxP) == (nothing, nothing)
        # the Q pair stays as delivered: an open band, where a missing pair
        # would give +-max(ratedS, baseMVA)
        @test (ps_deer.minQ, ps_deer.maxQ) == (-99999.9, 99999.9)
        local ps_ei1 = ph.net.prosumpsVec[last(ph.element_by_mrid[ei1.mrid])]
        @test (ps_ei1.minP, ps_ei1.maxP) == (nothing, nothing)
        @test "notice: P limits of 9999 MW and above read as no limit (placeholder): 1 EquivalentInjection, 1 GeneratingUnit" in ph.messages
        @test !any(m -> occursin("read as no limit", m), res.messages)
        @test ph.slack_bus == res.slack_bus
        # the slack unit drops out (no reference, no priority): ensureSlack!,
        # the replacement of an outaged slack in N-1, takes a real 800 MW
        # machine. The three EquivalentInjections are external units and
        # outrank every machine, so the copy leaves them out.
        local n1 = deepcopy(ph.net)
        filter!(p -> p.comp.cTyp != Sparlectra.ExternalNetworkInjection, n1.prosumpsVec)
        local k_slack = only(k for (k, p) in enumerate(n1.prosumpsVec) if Sparlectra.isSlack(p))
        n1.prosumpsVec[k_slack].referencePri = nothing
        n1.prosumpsVec[k_slack].referencePriority = 0
        local picked = ensureSlack!(n1; log = false)
        @test picked != Sparlectra.getPosumerBusIndex(ps_deer)
        @test only(p for p in n1.prosumpsVec if Sparlectra.isSlack(p)).maxP == 800.0
        # distributed slack: no pmax or headroom weight for the placeholder unit
        local bt = [Sparlectra.getNodeType(nd) == Sparlectra.Slack ? :Slack : Sparlectra.getNodeType(nd) == Sparlectra.PV ? :PV : :PQ for nd in ph.net.nodeVec]
        for mode in (:pmax_weighted, :headroom_weighted)
          local ds = Sparlectra.build_distributed_slack_state(ph.net, bt; p_mode = mode, fallback = :error)
          @test !any(g -> g.gen_index == k_deer, ds.gen_table)
          @test sum(ds.alpha) ≈ 1.0
        end
      end)() end

      # RegulatingCondEq.controlEnabled of a SynchronousMachine. All 19
      # SmallGrid machines carry true, each alone on its bus with a local
      # control whose target equals the SV voltage. In a copy of the store one
      # machine is switched off (with an SSH q of its own, so the sign shows),
      # another loses the attribute, and the slack machine is switched off
      # with its control target moved away from the SV voltage. The same flag
      # on a StaticVarCompensator and an ExternalNetworkInjection follows below.
      @testset "RegulatingCondEq controlEnabled" begin (function ()
        local ce_store = deepcopy(loaded)
        local sm_named(nm) = only(sm for sm in objectsOf(ce_store, :SynchronousMachine) if str(sm, :name) == nm)
        local sm_off, sm_kept, sm_slack = sm_named("Breed SM"), sm_named("DeerCrk SM"), sm_named("Sporn SM")
        sm_off.attrs[:controlEnabled] = "false"
        sm_off.attrs[:q] = "-37.5"
        delete!(sm_kept.attrs, :controlEnabled)
        sm_slack.attrs[:controlEnabled] = "false"
        local slack_tn = CI.busOfEquipment(topo, sm_slack, 1).tn
        ce_store.objects[sm_slack.refs[:RegulatingControl]].attrs[:targetValue] = string(1.02 * topo.vn_kV[slack_tn])
        local ce = importCGMES(path = ce_store, name = "sg_ctrl")
        local ps_of(r, sm) = r.net.prosumpsVec[last(r.element_by_mrid[sm.mrid])]
        local node_of(r, sm) = r.net.nodeVec[Sparlectra.getPosumerBusIndex(ps_of(r, sm))]
        # off: a PQ unit at the negated SSH p and q, Q limits as in the plain import
        local ps_off = ps_of(ce, sm_off)
        @test !Sparlectra.isRegulating(ps_off) && !Sparlectra.isSlack(ps_off)
        @test node_of(ce, sm_off)._nodeType == Sparlectra.PQ
        @test (ps_off.pVal, ps_off.qVal) == (-num(sm_off, :p), 37.5) == (450.0, 37.5)
        @test (ps_off.minQ, ps_off.maxQ) == (ps_of(res, sm_off).minQ, ps_of(res, sm_off).maxQ) == (-640.0, 640.0)
        # attribute removed: regulating as in the plain import
        @test Sparlectra.isRegulating(ps_of(ce, sm_kept))
        @test node_of(ce, sm_kept)._nodeType == Sparlectra.PV
        @test ps_of(ce, sm_kept).vm_pu == ps_of(res, sm_kept).vm_pu
        # the slack choice ignores the flag; the slack machine holds the SV
        # voltage of its bus, not the moved target
        @test ce.slack_bus == res.slack_bus
        @test Sparlectra.isSlack(ps_of(ce, sm_slack))
        @test ps_of(ce, sm_slack).vm_pu == CI._svVoltageByTN(ce_store)[slack_tn][1] / topo.vn_kV[slack_tn] == ps_of(res, sm_slack).vm_pu
        @test "notice: 2 SynchronousMachine(s) with controlEnabled = false, no voltage regulation: Breed SM, Sporn SM" in ce.messages
        @test count(m -> startswith(m, "notice: slack SynchronousMachine Sporn SM has controlEnabled = false"), ce.messages) == 1
        @test !any(m -> occursin("controlEnabled", m), res.messages)

        # StaticVarCompensator (SmallGrid has none): one added on the node of
        # a load whose bus is PQ, with an enabled voltage RegulatingControl at
        # its own terminal and the SV voltage as target, so the SV guard of the
        # SVC mapping keeps it PV. Flag absent: PV at the target. Flag false:
        # PQ at its SSH q, one notice.
        local ld = first(ec for ec in objectsOf(loaded, :EnergyConsumer) if net.nodeVec[net.busDict[tbus(ec, 1)]]._nodeType == Sparlectra.PQ)
        local ld_t = only(t for t in objectsOf(loaded, :Terminal) if get(t.refs, :ConductingEquipment, "") == ld.mrid)
        local svc_tn = ld_t.refs[:TopologicalNode]
        local svc_vm = CI._svVoltageByTN(loaded)[svc_tn][1] / topo.vn_kV[svc_tn]
        local svc_import(flag) = begin
          local st = deepcopy(loaded)
          local add!(mrid, class, attrs, refs) = (st.objects[mrid] = CI.CIMObject(mrid, class, attrs, refs, "svc_test"); push!(get!(Vector{String}, st.byclass, class), mrid))
          add!("_svc_rc", :RegulatingControl, Dict{Symbol,String}(:mode => "RegulatingControlModeKind.voltage", :enabled => "true", :targetValue => string(svc_vm * topo.vn_kV[svc_tn])), Dict{Symbol,String}(:Terminal => "_svc_t"))
          local attrs = Dict{Symbol,String}(:name => "SVC_TEST", :q => "-12.5", :capacitiveRating => "500", :inductiveRating => "-500")
          flag === nothing || (attrs[:controlEnabled] = flag)
          add!("_svc", :StaticVarCompensator, attrs, Dict{Symbol,String}(:RegulatingControl => "_svc_rc"))
          add!("_svc_t", :Terminal, Dict{Symbol,String}(:sequenceNumber => "1", :connected => "true"), Dict{Symbol,String}(:ConductingEquipment => "_svc", :TopologicalNode => svc_tn, :ConnectivityNode => ld_t.refs[:ConnectivityNode]))
          local r = importCGMES(path = st, name = "sg_svc")
          (r, r.net.prosumpsVec[last(r.element_by_mrid["_svc"])])
        end
        local svc_on, ps_svc_on = svc_import(nothing)
        @test Sparlectra.isRegulating(ps_svc_on) && isapprox(ps_svc_on.vm_pu, svc_vm; atol = 1e-12)
        @test svc_on.net.nodeVec[Sparlectra.getPosumerBusIndex(ps_svc_on)]._nodeType == Sparlectra.PV
        local svc_off, ps_svc_off = svc_import("false")
        @test !Sparlectra.isRegulating(ps_svc_off)
        @test svc_off.net.nodeVec[Sparlectra.getPosumerBusIndex(ps_svc_off)]._nodeType == Sparlectra.PQ
        @test (ps_svc_off.pVal, ps_svc_off.qVal) == (0.0, 12.5)
        @test "notice: 1 StaticVarCompensator(s) with controlEnabled = false, no voltage regulation: SVC_TEST" in svc_off.messages

        # ExternalNetworkInjection: sp_case14 has one, Gen_110 on Ostheim_110,
        # the slack by referencePriority 1, local control target 1.02 pu equal
        # to the SV voltage, no Q limits. Flag false with the target moved to
        # 1.05 pu: still the slack, held at the SV voltage. Flag false with
        # the priority moved to a machine: a PQ bus at the SSH p/q.
        local s14 = loadCGMES(cgmes_fixture_dir("sp_case14"))
        local eni = only(objectsOf(s14, :ExternalNetworkInjection))
        local r14 = importCGMES(path = s14, name = "s14")
        local eni_tn = CI.busOfEquipment(r14.topo, eni, 1).tn
        local eni_sv = CI._svVoltageByTN(s14)[eni_tn][1] / r14.topo.vn_kV[eni_tn]
        local ps_eni(r) = r.net.prosumpsVec[last(r.element_by_mrid[eni.mrid])]
        @test Sparlectra.isSlack(ps_eni(r14)) && ps_eni(r14).vm_pu == eni_sv == 1.02
        local s14_slack = deepcopy(s14)
        s14_slack.objects[eni.mrid].attrs[:controlEnabled] = "false"
        s14_slack.objects[eni.refs[:RegulatingControl]].attrs[:targetValue] = string(1.05 * r14.topo.vn_kV[eni_tn])
        local e_slack = importCGMES(path = s14_slack, name = "s14_eni_slack")
        @test e_slack.slack_bus == r14.slack_bus == "Ostheim_110"
        @test Sparlectra.isSlack(ps_eni(e_slack)) && ps_eni(e_slack).vm_pu == eni_sv
        @test "notice: 1 ExternalNetworkInjection(s) with controlEnabled = false, no voltage regulation: Gen_110" in e_slack.messages
        @test count(m -> startswith(m, "notice: slack ExternalNetworkInjection Gen_110 has controlEnabled = false"), e_slack.messages) == 1
        local s14_pq = deepcopy(s14)
        s14_pq.objects[eni.mrid].attrs[:controlEnabled] = "false"
        delete!(s14_pq.objects[eni.mrid].attrs, :referencePriority)
        local sm14 = first(objectsOf(s14_pq, :SynchronousMachine))
        sm14.attrs[:referencePriority] = "1"
        local e_pq = importCGMES(path = s14_pq, name = "s14_eni_pq")
        @test e_pq.slack_bus == CI.busOfEquipment(r14.topo, sm14, 1).bus != "Ostheim_110"
        @test !Sparlectra.isRegulating(ps_eni(e_pq)) && !Sparlectra.isSlack(ps_eni(e_pq))
        @test e_pq.net.nodeVec[e_pq.net.busDict["Ostheim_110"]]._nodeType == Sparlectra.PQ
        @test (ps_eni(e_pq).pVal, ps_eni(e_pq).qVal) == (-num(eni, :p), -num(eni, :q))

        # EquivalentInjection: no controlEnabled, it regulates only with EQ
        # regulationCapability and SSH regulationStatus both true. The three
        # SmallGrid units carry both false. In a copy each gets the SV voltage
        # of its node as target: the first with both flags true (PV), the
        # second with status true and capability false, the third with
        # capability true and no status (both PQ, the third counted).
        local ei_store = deepcopy(loaded)
        local eis = sort!(collect(objectsOf(ei_store, :EquivalentInjection)); by = o -> o.mrid)
        @test length(eis) == 3 && all(o -> boolval(o, :regulationCapability) === boolval(o, :regulationStatus) === false, eis)
        for o in eis
          local o_tn = CI.busOfEquipment(topo, o, 1).tn
          o.attrs[:regulationTarget] = string(CI._svVoltageByTN(loaded)[o_tn][1])
        end
        eis[1].attrs[:regulationCapability] = "true"
        eis[1].attrs[:regulationStatus] = "true"
        eis[2].attrs[:regulationStatus] = "true"
        eis[3].attrs[:regulationCapability] = "true"
        delete!(eis[3].attrs, :regulationStatus)
        local ei_r = importCGMES(path = ei_store, name = "sg_ei_reg")
        local ps_ei(o) = ei_r.net.prosumpsVec[last(ei_r.element_by_mrid[o.mrid])]
        local ei_node(o) = ei_r.net.nodeVec[Sparlectra.getPosumerBusIndex(ps_ei(o))]
        @test Sparlectra.isRegulating(ps_ei(eis[1])) && ei_node(eis[1])._nodeType == Sparlectra.PV
        @test !Sparlectra.isRegulating(ps_ei(eis[2])) && ei_node(eis[2])._nodeType == Sparlectra.PQ
        @test !Sparlectra.isRegulating(ps_ei(eis[3])) && ei_node(eis[3])._nodeType == Sparlectra.PQ
        @test count(m -> occursin("regulationStatus", m), ei_r.messages) == 1
        @test "notice: 1 EquivalentInjection(s) with regulationCapability = true and no SSH regulationStatus, no voltage regulation" in ei_r.messages
        @test !any(m -> occursin("regulationStatus", m), res.messages)
      end)() end

      # cgmes_import.map_energy_source. SmallGrid has no EnergySource, so two
      # are added to a copy of the store on the nodes of two loads, each with
      # a connected terminal and an SvPowerFlow equal to its SSH values (both
      # in load sign convention). The second one is out of service.
      @testset "EnergySource behind cgmes_import.map_energy_source" begin (function ()
        local es_store = deepcopy(loaded)
        local es_hosts = collect(objectsOf(es_store, :EnergyConsumer))[1:2]
        local es_tn = Dict{String,String}()
        for (i, (mrid, nm, p, q, in_service)) in enumerate((("_es_test_1", "ES_TEST_1", -25.5, -4.25, true), ("_es_test_2", "ES_TEST_2", -7.0, 1.5, false)))
          local host_t = only(t for t in objectsOf(es_store, :Terminal) if get(t.refs, :ConductingEquipment, "") == es_hosts[i].mrid)
          local es_attrs = Dict{Symbol,String}(:name => nm, :activePower => string(p), :reactivePower => string(q))
          in_service || (es_attrs[:inService] = "false")
          local t_mrid = string(mrid, "_t")
          local sv_mrid = string(mrid, "_sv")
          es_tn[mrid] = host_t.refs[:TopologicalNode]
          es_store.objects[mrid] = CI.CIMObject(mrid, :EnergySource, es_attrs, Dict{Symbol,String}(), "es_test")
          es_store.objects[t_mrid] = CI.CIMObject(t_mrid, :Terminal, Dict{Symbol,String}(:sequenceNumber => "1", :connected => "true"), Dict{Symbol,String}(:ConductingEquipment => mrid, :TopologicalNode => host_t.refs[:TopologicalNode], :ConnectivityNode => host_t.refs[:ConnectivityNode]), "es_test")
          es_store.objects[sv_mrid] = CI.CIMObject(sv_mrid, :SvPowerFlow, Dict{Symbol,String}(:p => string(p), :q => string(q)), Dict{Symbol,String}(:Terminal => t_mrid), "es_test")
          push!(get!(Vector{String}, es_store.byclass, :EnergySource), mrid)
          push!(es_store.byclass[:Terminal], t_mrid)
          push!(es_store.byclass[:SvPowerFlow], sv_mrid)
        end
        local es_key = Sparlectra._measurement_key("_es_test_1_t")

        # off (the default): no prosumer, both listed as not mapped, one notice
        # that names the switch and carries no number of its own
        local es_off = importCGMES(path = es_store, name = "sg_es_off")
        @test length(es_off.net.prosumpsVec) == length(net.prosumpsVec)
        @test all(m -> !haskey(es_off.element_by_mrid, m) && m in es_off.skipped_equipment, ("_es_test_1", "_es_test_2"))
        @test count(m -> startswith(m, "notice: EnergySource not mapped"), es_off.messages) == 1
        @test "notice: EnergySource not mapped, switch cgmes_import.map_energy_source (count and SSH P sum in the table of classes not mapped)" in es_off.messages
        @test !any(r -> r.kind == :units && occursin("ES_TEST", r.name), compareWithSV(es_off).flows.rows)
        @test CI.cgmesMeasurementKeys(es_off).terminals[es_key].kind == :not_modelled
        # the report has one row per class with SSH active power and no
        # object in the network: here EnergySource, count and SSH P sum, below
        # one header that states the load sign. SmallGrid itself maps every
        # such class, so it has neither header nor row.
        @test CI._unmappedClasses(es_off.store, es_off.element_by_mrid) == [(class = :EnergySource, count = 2, p_mw = -32.5, no_p = 0)]
        @test count(m -> startswith(m, "notice: class not mapped:"), es_off.messages) == 1
        @test "notice: class not mapped: EnergySource, 2 object(s), P sum -32.5 MW" in es_off.messages
        @test count(==("notice: classes not mapped, P sum from SSH in load sign (positive means flow out of the node):"), es_off.messages) == 1
        @test !any(m -> startswith(m, "notice: class not mapped:") || startswith(m, "notice: classes not mapped"), res.messages)

        # on: a fixed PQ injection on the load's bus whose P and Q are the
        # negated SvPowerFlow at its terminal; no regulation, no Q limits, no slack
        local es_on = importCGMES(path = es_store, name = "sg_es_on", map_energy_source = true)
        @test length(es_on.net.prosumpsVec) == length(net.prosumpsVec) + 1
        local es_addr = es_on.element_by_mrid["_es_test_1"]
        local es_ps = es_on.net.prosumpsVec[last(es_addr)]
        local es_sv = es_store.objects["_es_test_1_sv"]
        @test first(es_addr) == :generator
        @test (es_ps.pVal, es_ps.qVal) == (-num(es_sv, :p), -num(es_sv, :q)) == (25.5, 4.25)
        @test Sparlectra._bus_name_by_idx(es_on.net)[Sparlectra.getPosumerBusIndex(es_ps)] == es_on.topo.bus_name[es_tn["_es_test_1"]]
        @test Sparlectra.isGenerator(es_ps)
        @test !Sparlectra.isRegulating(es_ps)
        @test !Sparlectra.isSlack(es_ps)
        @test es_ps.minQ === nothing && es_ps.maxQ === nothing
        @test es_on.slack_bus == res.slack_bus
        @test any(r -> r.kind == :units && occursin("ES_TEST_1", r.name), compareWithSV(es_on).flows.rows)
        @test CI.cgmesMeasurementKeys(es_on).terminals[es_key].kind == :injection
        # inService = false: skipped with its own message and counted
        @test !haskey(es_on.element_by_mrid, "_es_test_2") && "_es_test_2" in es_on.skipped_equipment
        @test any(m -> startswith(m, "skip: EnergySource ES_TEST_2 (out of service"), es_on.messages)
        @test any(m -> startswith(m, "notice: map_energy_source = true: 1 EnergySource object(s) mapped as PQ injections (P = 25.5 MW, Q = 4.25 MVAr injected), 1 skipped"), es_on.messages)
        @test !any(m -> startswith(m, "notice: electrical island"), es_on.messages)
        # one EnergySource is in the network, so the class has no line
        @test isempty(CI._unmappedClasses(es_on.store, es_on.element_by_mrid))
        @test !any(m -> startswith(m, "notice: class not mapped:"), es_on.messages)
        # the same EnergySource moved onto the node of a regulating machine
        # (Breed SM, alone on its PV bus, Q limits -640/640 MVAr): the bus
        # band is the machine's limits shifted by the source's fixed Q, the
        # source without Q limits does not open it to -Inf/Inf
        local pv_store = deepcopy(es_store)
        local breed = only(sm for sm in objectsOf(pv_store, :SynchronousMachine) if str(sm, :name) == "Breed SM")
        local breed_t = only(t for t in objectsOf(pv_store, :Terminal) if get(t.refs, :ConductingEquipment, "") == breed.mrid)
        for key in (:TopologicalNode, :ConnectivityNode)
          pv_store.objects["_es_test_1_t"].refs[key] = breed_t.refs[key]
        end
        local es_pv = importCGMES(path = pv_store, name = "sg_es_pv", map_energy_source = true)
        local breed_ps = es_pv.net.prosumpsVec[last(es_pv.element_by_mrid[breed.mrid])]
        local breed_bus = Sparlectra.getPosumerBusIndex(breed_ps)
        @test Sparlectra.getPosumerBusIndex(es_pv.net.prosumpsVec[last(es_pv.element_by_mrid["_es_test_1"])]) == breed_bus
        @test es_pv.net.nodeVec[breed_bus]._nodeType == Sparlectra.PV
        @test (breed_ps.minQ, breed_ps.maxQ) == (-640.0, 640.0)
        local qmin_pv, qmax_pv = Sparlectra.getQLimits_pu(es_pv.net)
        @test isapprox(qmin_pv[breed_bus] * es_pv.net.baseMVA, -640.0 + 4.25; atol = 1e-9)
        @test isapprox(qmax_pv[breed_bus] * es_pv.net.baseMVA, 640.0 + 4.25; atol = 1e-9)
        # boundary objects do not count; an object without an SSH value is
        # counted and adds nothing to the sum
        local pec(mrid, attrs) = CI.CIMObject(mrid, :PowerElectronicsConnection, attrs, Dict{Symbol,String}(), "pec_test")
        local pec_store = CI.CGMESStore(
          Dict("_pec_1" => pec("_pec_1", Dict(:p => "-10")), "_pec_2" => pec("_pec_2", Dict{Symbol,String}()), "_pec_bd" => pec("_pec_bd", Dict(:p => "-99"))),
          Dict(:PowerElectronicsConnection => ["_pec_1", "_pec_2", "_pec_bd"]),
          CI.CGMESFileInfo[],
          "3.0",
          Set(["_pec_bd"]),
          CI.CGMESMultiRefs(),
        )
        @test CI._unmappedClasses(pec_store, Dict{String,Tuple{Symbol,Int}}()) == [(class = :PowerElectronicsConnection, count = 2, p_mw = -10.0, no_p = 1)]
        # with map_power_electronics on the class is one the importer maps:
        # objects without a connected terminal are skipped on purpose, no row
        @test isempty(CI._unmappedClasses(pec_store, Dict{String,Tuple{Symbol,Int}}(); map_power_electronics = true))
        # every class with terminals counts, also without SSH active power: a
        # Ground (no SSH value, the importer does not map the class) on the
        # node of a load gets a row without P sum. A StationSupply out of
        # service is a class the importer maps and skips, so it gets none.
        # SmallGrid's own switches and busbar sections are consumed by the
        # topology (above: no row on the plain import).
        local gr_store = deepcopy(loaded)
        local gr_host = first(objectsOf(gr_store, :EnergyConsumer))
        local gr_ht = only(t for t in objectsOf(gr_store, :Terminal) if get(t.refs, :ConductingEquipment, "") == gr_host.mrid)
        for (mrid, cls, attrs) in (("_gr_test", :Ground, Dict{Symbol,String}(:name => "GR_TEST")), ("_ss_test", :StationSupply, Dict{Symbol,String}(:name => "SS_TEST", :p => "1.5", :q => "0.5", :inService => "false")))
          local t_mrid = string(mrid, "_t")
          gr_store.objects[mrid] = CI.CIMObject(mrid, cls, attrs, Dict{Symbol,String}(), "gr_test")
          gr_store.objects[t_mrid] = CI.CIMObject(t_mrid, :Terminal, Dict{Symbol,String}(:sequenceNumber => "1", :connected => "true"), Dict{Symbol,String}(:ConductingEquipment => mrid, :TopologicalNode => gr_ht.refs[:TopologicalNode], :ConnectivityNode => gr_ht.refs[:ConnectivityNode]), "gr_test")
          push!(get!(Vector{String}, gr_store.byclass, cls), mrid)
          push!(gr_store.byclass[:Terminal], t_mrid)
        end
        local gr = importCGMES(path = gr_store, name = "sg_ground")
        @test "_ss_test" in gr.skipped_equipment
        @test CI._unmappedClasses(gr.store, gr.element_by_mrid) == [(class = :Ground, count = 1, p_mw = nothing, no_p = 1)]
        @test filter(m -> occursin("classes not mapped", m) || startswith(m, "notice: class not mapped:"), gr.messages) ==
              ["notice: classes not mapped, P sum from SSH in load sign (positive means flow out of the node):", "notice: class not mapped: Ground, 1 object(s)"]

        # an island whose only source is an EnergySource: two new nodes, one
        # line, the source on one node and a load on the other. The import
        # makes the source the island's slack (largest infeed, both buses of
        # degree 1) with a warning. With multi_slack = false it gives the
        # island no slack and names it; the island-wise solve then promotes
        # the source to the island reference (detect_ac_islands)
        local isl_store = deepcopy(es_store)
        local isl_bv = isl_store.objects[es_tn["_es_test_1"]].refs[:BaseVoltage]
        local isl_add!(mrid, class, attrs, refs) = (isl_store.objects[mrid] = CI.CIMObject(mrid, class, attrs, refs, "es_test"); push!(get!(Vector{String}, isl_store.byclass, class), mrid))
        local isl_term!(mrid, eq, tn, seq) = isl_add!(mrid, :Terminal, Dict{Symbol,String}(:sequenceNumber => string(seq), :connected => "true"), Dict{Symbol,String}(:ConductingEquipment => eq, :TopologicalNode => tn))
        for tn in ("_isl_tn_1", "_isl_tn_2")
          isl_add!(tn, :TopologicalNode, Dict{Symbol,String}(:name => uppercase(tn)), Dict{Symbol,String}(:BaseVoltage => isl_bv))
        end
        isl_add!("_isl_line", :ACLineSegment, Dict{Symbol,String}(:name => "ISL_LINE", :r => "1.0", :x => "5.0", :bch => "0"), Dict{Symbol,String}(:BaseVoltage => isl_bv))
        isl_term!("_isl_line_t1", "_isl_line", "_isl_tn_1", 1)
        isl_term!("_isl_line_t2", "_isl_line", "_isl_tn_2", 2)
        isl_add!("_es_test_3", :EnergySource, Dict{Symbol,String}(:name => "ES_TEST_3", :activePower => "-5", :reactivePower => "-1"), Dict{Symbol,String}())
        isl_term!("_es_test_3_t", "_es_test_3", "_isl_tn_1", 1)
        isl_add!("_isl_load", :EnergyConsumer, Dict{Symbol,String}(:name => "ISL_LOAD", :p => "4", :q => "1"), Dict{Symbol,String}())
        isl_term!("_isl_load_t", "_isl_load", "_isl_tn_2", 1)
        local es_isl_slack = @test_logs (:warn, r"^CGMES import: 1 data defect") importCGMES(path = isl_store, name = "sg_es_island_slack", map_energy_source = true)
        @test Sparlectra.isSlack(es_isl_slack.net.prosumpsVec[last(es_isl_slack.element_by_mrid["_es_test_3"])])
        @test count(m -> startswith(m, "warning: electrical island of 2 bus(es) around ") && occursin("ES_TEST_3 at ", m) && occursin("(largest infeed, P = 5.0 MW, only degree-1 buses carry infeed)", m), es_isl_slack.messages) == 1
        @test !any(m -> startswith(m, "notice: electrical island"), es_isl_slack.messages)
        local es_isl = importCGMES(path = isl_store, name = "sg_es_island", map_energy_source = true, multi_slack = false)
        local isl_bus = Sparlectra.getPosumerBusIndex(es_isl.net.prosumpsVec[last(es_isl.element_by_mrid["_es_test_3"])])
        @test !(isl_bus in es_isl.net.isoNodes)
        @test count(m -> startswith(m, "notice: electrical island of 2 bus(es) around ") && occursin("no source but EnergySource (ES_TEST_3); the island-wise solve makes it the island reference", m), es_isl.messages) == 1
        @test count(m -> startswith(m, "notice: electrical island"), es_isl.messages) == 1
        local isl_row = only(r for r in Sparlectra.detect_ac_islands(es_isl.net; promote_generators = true).rows if isl_bus in r.buses)
        @test (isl_row.status, isl_row.chosen_ref_bus) == ("promote_generator_ref", isl_bus)

        # the configured import forwards the switch, and the loader reads it
        @test !SparlectraConfig().cgmes.map_energy_source
        local es_cfg = SparlectraConfig(cgmes = Sparlectra.CGMESImportConfig(map_energy_source = true))
        @test first(importCGMES(es_cfg; path = es_store, name = "sg_es_cfg").element_by_mrid["_es_test_1"]) == :generator
        @test Sparlectra.CGMESImportConfig(Dict{String,Any}("cgmes_import" => Dict{String,Any}("map_energy_source" => true))).map_energy_source
      end)() end

      # PowerElectronicsConnection (CGMES 3.0) behind
      # cgmes_import.map_power_electronics. No shipped delivery carries one
      # (SmallGrid is 2.4.15, the cgmes_demo exports have none), so
      # connections are added to a copy of the store on load buses the plain
      # import leaves PQ, each with a connected terminal, SSH p and q in load
      # sign and an SvPowerFlow equal to them; some with a voltage
      # RegulatingControl at a load terminal.
      @testset "PowerElectronicsConnection behind cgmes_import.map_power_electronics" begin (function ()
        local pe_store = deepcopy(loaded)
        local add!(s, mrid, class, attrs, refs) = (s.objects[mrid] = CI.CIMObject(mrid, class, attrs, refs, "pec_test"); push!(get!(Vector{String}, s.byclass, class), mrid))
        local ce_terms = Dict{String,Vector{CI.CIMObject}}()
        for t in objectsOf(pe_store, :Terminal)
          push!(get!(Vector{CI.CIMObject}, ce_terms, get(t.refs, :ConductingEquipment, "")), t)
        end
        # one load terminal per PQ bus of the plain import
        local pq_terms = CI.CIMObject[]
        for ec in objectsOf(pe_store, :EnergyConsumer)
          local t = only(ce_terms[ec.mrid])
          local tn = t.refs[:TopologicalNode]
          (any(u -> u.refs[:TopologicalNode] == tn, pq_terms) || net.nodeVec[net.busDict[topo.bus_name[tn]]]._nodeType != Sparlectra.PQ) && continue
          push!(pq_terms, t)
        end
        local tn_of(t) = t.refs[:TopologicalNode]
        local kv(t, vm) = string(vm * topo.vn_kV[tn_of(t)])
        local function pec!(s, mrid, nm, host, p, q; extra = Dict{Symbol,String}(), rc = nothing)
          local refs = Dict{Symbol,String}()
          if rc !== nothing
            add!(s, string(mrid, "_rc"), :RegulatingControl, merge(Dict{Symbol,String}(:mode => "RegulatingControlModeKind.voltage"), rc.attrs), Dict{Symbol,String}(:Terminal => rc.terminal.mrid))
            refs[:RegulatingControl] = string(mrid, "_rc")
          end
          add!(s, mrid, :PowerElectronicsConnection, merge(Dict{Symbol,String}(:name => nm, :p => string(p), :q => string(q)), extra), refs)
          local trefs = Dict{Symbol,String}(:ConductingEquipment => mrid, :TopologicalNode => tn_of(host))
          haskey(host.refs, :ConnectivityNode) && (trefs[:ConnectivityNode] = host.refs[:ConnectivityNode])
          add!(s, string(mrid, "_t"), :Terminal, Dict{Symbol,String}(:sequenceNumber => "1", :connected => "true"), trefs)
          add!(s, string(mrid, "_sv"), :SvPowerFlow, Dict{Symbol,String}(:p => string(p), :q => string(q)), Dict{Symbol,String}(:Terminal => string(mrid, "_t")))
          return nothing
        end
        # PQ, with Q limits and two units (P limits -5 + 0, 5 + 10)
        pec!(pe_store, "_pec_pq", "PEC_PQ", pq_terms[1], -12.5, -3.0; extra = Dict{Symbol,String}(:minQ => "-8", :maxQ => "6", :ratedS => "15"))
        add!(pe_store, "_pec_pq_pv", :PhotoVoltaicUnit, Dict{Symbol,String}(:minP => "0", :maxP => "10"), Dict{Symbol,String}(:PowerElectronicsConnection => "_pec_pq"))
        add!(pe_store, "_pec_pq_bat", :BatteryUnit, Dict{Symbol,String}(:minP => "-5", :maxP => "5"), Dict{Symbol,String}(:PowerElectronicsConnection => "_pec_pq"))
        # local voltage control at 1.02 pu, no minQ/maxQ
        pec!(pe_store, "_pec_pv", "PEC_PV", pq_terms[2], -20.0, -1.0; extra = Dict{Symbol,String}(:ratedS => "30", :controlEnabled => "true"), rc = (terminal = pq_terms[2], attrs = Dict{Symbol,String}(:enabled => "true", :targetValue => kv(pq_terms[2], 1.02))))
        # the same control, switched off at the connection or at the control
        pec!(pe_store, "_pec_ctrl_off", "PEC_CTRL_OFF", pq_terms[3], -4.0, 0.0; extra = Dict{Symbol,String}(:controlEnabled => "false"), rc = (terminal = pq_terms[3], attrs = Dict{Symbol,String}(:enabled => "true", :targetValue => kv(pq_terms[3], 1.02))))
        pec!(pe_store, "_pec_rc_off", "PEC_RC_OFF", pq_terms[4], -3.0, 0.0; extra = Dict{Symbol,String}(:controlEnabled => "true"), rc = (terminal = pq_terms[4], attrs = Dict{Symbol,String}(:enabled => "false", :targetValue => kv(pq_terms[4], 1.02))))
        # voltage control of another load bus, at its SV voltage
        local remote_vm = net.nodeVec[net.busDict[topo.bus_name[tn_of(pq_terms[6])]]]._vm_pu
        pec!(pe_store, "_pec_remote", "PEC_REMOTE", pq_terms[5], -6.0, 0.0; extra = Dict{Symbol,String}(:controlEnabled => "true", :minQ => "-10", :maxQ => "10"), rc = (terminal = pq_terms[6], attrs = Dict{Symbol,String}(:enabled => "true", :targetValue => kv(pq_terms[6], remote_vm))))
        pec!(pe_store, "_pec_oos", "PEC_OOS", pq_terms[1], -1.0, 0.0; extra = Dict{Symbol,String}(:inService => "false"))
        local pec_mrids = ("_pec_pq", "_pec_pv", "_pec_ctrl_off", "_pec_rc_off", "_pec_remote", "_pec_oos")

        # off (the default): no prosumer, every connection listed as not
        # mapped, one notice that names the switch, and the row in the table
        # of classes not mapped with all six objects and their SSH P sum
        local pe_off = importCGMES(path = pe_store, name = "sg_pec_off")
        @test length(pe_off.net.prosumpsVec) == length(net.prosumpsVec)
        @test all(m -> !haskey(pe_off.element_by_mrid, m) && m in pe_off.skipped_equipment, pec_mrids)
        @test count(m -> occursin("PowerElectronicsConnection", m), pe_off.messages) == 2
        @test "notice: PowerElectronicsConnection not mapped, switch cgmes_import.map_power_electronics (count and SSH P sum in the table of classes not mapped)" in pe_off.messages
        @test "notice: class not mapped: PowerElectronicsConnection, 6 object(s), P sum -46.5 MW" in pe_off.messages
        @test !any(r -> r.kind == :units && occursin("PEC_", r.name), compareWithSV(pe_off).flows.rows)
        @test CI.cgmesMeasurementKeys(pe_off).terminals[Sparlectra._measurement_key("_pec_pq_t")].kind == :not_modelled
        @test pe_off.slack_bus == res.slack_bus

        # on: the mapping
        local pe = importCGMES(path = pe_store, name = "sg_pec", map_power_electronics = true)
        local ps_of(m) = pe.net.prosumpsVec[last(pe.element_by_mrid[m])]
        local node_of(m) = pe.net.nodeVec[Sparlectra.getPosumerBusIndex(ps_of(m))]
        @test length(pe.net.prosumpsVec) == length(net.prosumpsVec) + 5
        @test all(m -> first(pe.element_by_mrid[m]) == :generator, ("_pec_pq", "_pec_pv", "_pec_ctrl_off", "_pec_rc_off", "_pec_remote"))
        # PQ: the negated SSH values on the host's bus, the minQ/maxQ hull and
        # the summed unit limits
        local pq = ps_of("_pec_pq")
        @test (pq.pVal, pq.qVal) == (-num(pe_store.objects["_pec_pq_sv"], :p), -num(pe_store.objects["_pec_pq_sv"], :q)) == (12.5, 3.0)
        @test Sparlectra._bus_name_by_idx(pe.net)[Sparlectra.getPosumerBusIndex(pq)] == pe.topo.bus_name[tn_of(pq_terms[1])]
        @test Sparlectra.isGenerator(pq) && !Sparlectra.isRegulating(pq) && !Sparlectra.isSlack(pq)
        @test (pq.minQ, pq.maxQ, pq.minP, pq.maxP) == (-8.0, 8.0, -5.0, 15.0)
        @test node_of("_pec_pq")._nodeType == Sparlectra.PQ
        # an enabled voltage control at its own bus: PV at the target, the wide
        # Q limits max(ratedS, baseMVA), and a notice
        local pv = ps_of("_pec_pv")
        @test Sparlectra.isRegulating(pv) && !Sparlectra.isSlack(pv)
        @test node_of("_pec_pv")._nodeType == Sparlectra.PV
        @test pv.vm_pu ≈ 1.02
        @test (pv.minQ, pv.maxQ) == (-100.0, 100.0)
        @test any(m -> startswith(m, "notice: PowerElectronicsConnection PEC_PV regulates voltage and has no usable minQ/maxQ"), pe.messages)
        # controlEnabled = false or RegulatingControl.enabled = false: PQ
        for m in ("_pec_ctrl_off", "_pec_rc_off")
          @test !Sparlectra.isRegulating(ps_of(m))
          @test node_of(m)._nodeType == Sparlectra.PQ
        end
        @test (ps_of("_pec_ctrl_off").pVal, ps_of("_pec_rc_off").pVal) == (4.0, 3.0)
        # a PQ connection without minQ/maxQ keeps no Q limits (only a
        # regulating one gets the wide fallback)
        @test all(m -> (ps_of(m).minQ, ps_of(m).maxQ) == (nothing, nothing), ("_pec_ctrl_off", "_pec_rc_off"))
        # a remote target: held PV at its own bus without machine_control
        @test node_of("_pec_remote")._nodeType == Sparlectra.PV
        @test any(m -> occursin("PEC_REMOTE has a remote voltage RegulatingControl", m), pe.messages)
        # inService = false: skipped and counted
        @test !haskey(pe.element_by_mrid, "_pec_oos") && "_pec_oos" in pe.skipped_equipment
        # the short-circuit harvest keeps the mapped connections only, on the
        # result and on the network; off, it keeps none
        @test sort!([r.mrid for r in pe.shortcircuit.power_electronics_connections]) == sort!(["_pec_pq", "_pec_pv", "_pec_ctrl_off", "_pec_rc_off", "_pec_remote"])
        @test pe.net.sc_sources.power_electronics_connections == pe.shortcircuit.power_electronics_connections
        @test isempty(pe_off.shortcircuit.power_electronics_connections)
        @test "notice: 5 PowerElectronicsConnection object(s) mapped (2 voltage-regulating, 3 PQ; P = 45.5 MW, Q = 4.0 MVAr injected), 1 skipped" in pe.messages
        @test pe.slack_bus == res.slack_bus
        # the class reached the network, so the report has no line for it
        @test isempty(CI._unmappedClasses(pe.store, pe.element_by_mrid))
        @test !any(m -> startswith(m, "notice: class not mapped:"), pe.messages)
        @test any(r -> r.kind == :units && occursin("PEC_PQ", r.name), compareWithSV(pe).flows.rows)
        @test CI.cgmesMeasurementKeys(pe).terminals[Sparlectra._measurement_key("_pec_pq_t")].kind == :injection
        @test last(runpf!(pe.net, 30, 1e-8, 0)) == 0
        @test abs(node_of("_pec_pv")._vm_pu - 1.02) < 1e-9

        # with machine_control the remote target gets an outer-loop controller
        local pe_mc = importCGMES(path = pe_store, name = "sg_pec_mc", machine_control = true, map_power_electronics = true)
        @test any(m -> startswith(m, "machine control: PEC_REMOTE → voltage"), pe_mc.messages)
        @test pe_mc.net.nodeVec[Sparlectra.getPosumerBusIndex(pe_mc.net.prosumpsVec[last(pe_mc.element_by_mrid["_pec_remote"])])]._nodeType == Sparlectra.PQ

        # an island whose only source is a PQ connection is treated like an
        # EnergySource island: the import makes the connection the island's
        # slack (largest infeed, both buses of degree 1) with a warning, and
        # with multi_slack = false gives the island no slack and the notice.
        # The warning also shows that the SV angle reference below did not
        # make it the slack: that island would have had a reference already
        local isl_store = deepcopy(loaded)
        local isl_bv = isl_store.objects[tn_of(pq_terms[1])].refs[:BaseVoltage]
        local isl_term!(mrid, eq, tn, seq) = add!(isl_store, mrid, :Terminal, Dict{Symbol,String}(:sequenceNumber => string(seq), :connected => "true"), Dict{Symbol,String}(:ConductingEquipment => eq, :TopologicalNode => tn))
        for tn in ("_pisl_tn_1", "_pisl_tn_2")
          add!(isl_store, tn, :TopologicalNode, Dict{Symbol,String}(:name => uppercase(tn)), Dict{Symbol,String}(:BaseVoltage => isl_bv))
        end
        add!(isl_store, "_pisl_line", :ACLineSegment, Dict{Symbol,String}(:name => "PISL_LINE", :r => "1.0", :x => "5.0", :bch => "0"), Dict{Symbol,String}(:BaseVoltage => isl_bv))
        isl_term!("_pisl_line_t1", "_pisl_line", "_pisl_tn_1", 1)
        isl_term!("_pisl_line_t2", "_pisl_line", "_pisl_tn_2", 2)
        add!(isl_store, "_pisl_pec", :PowerElectronicsConnection, Dict{Symbol,String}(:name => "PEC_ISL", :p => "-5", :q => "-1"), Dict{Symbol,String}())
        isl_term!("_pisl_pec_t", "_pisl_pec", "_pisl_tn_1", 1)
        add!(isl_store, "_pisl_load", :EnergyConsumer, Dict{Symbol,String}(:name => "PISL_LOAD", :p => "4", :q => "1"), Dict{Symbol,String}())
        isl_term!("_pisl_load_t", "_pisl_load", "_pisl_tn_2", 1)
        # the SV angle reference of the island is the connection's node, the
        # one place where the importer would give a unit the slack role
        add!(isl_store, "_pisl_island", :TopologicalIsland, Dict{Symbol,String}(:name => "PISL_ISLAND"), Dict{Symbol,String}(:AngleRefTopologicalNode => "_pisl_tn_1"))
        local pe_isl_slack = @test_logs (:warn, r"^CGMES import: 1 data defect") importCGMES(path = isl_store, name = "sg_pec_island_slack", map_power_electronics = true)
        @test Sparlectra.isSlack(pe_isl_slack.net.prosumpsVec[last(pe_isl_slack.element_by_mrid["_pisl_pec"])])
        @test count(m -> startswith(m, "warning: electrical island of 2 bus(es) around ") && occursin("PEC_ISL at ", m) && occursin("(largest infeed, P = 5.0 MW, only degree-1 buses carry infeed)", m) && endswith(m, "its unit type is unknown and assumed grid-following, an island of grid-following units only cannot stand in operation"), pe_isl_slack.messages) == 1
        @test !any(m -> startswith(m, "notice: electrical island"), pe_isl_slack.messages)
        local pe_isl = importCGMES(path = isl_store, name = "sg_pec_island", map_power_electronics = true, multi_slack = false)
        @test count(m -> startswith(m, "notice: electrical island"), pe_isl.messages) == 1
        @test any(m -> startswith(m, "notice: electrical island of 2 bus(es) around ") && occursin("no source but PowerElectronicsConnection (PEC_ISL); the island-wise solve makes it the island reference", m), pe_isl.messages)
        local pisl = pe_isl.net.prosumpsVec[last(pe_isl.element_by_mrid["_pisl_pec"])]
        @test !Sparlectra.isSlack(pisl)
        @test pe_isl.net.nodeVec[Sparlectra.getPosumerBusIndex(pisl)]._nodeType == Sparlectra.PQ

        # short circuit: a connection is no synchronous source. On the bus of
        # a machine without ratedS and maxOperatingP the maxP of the
        # connection's unit is not taken as the machine's size (Ik'' and the
        # machine's reason as without the connection), and the connection
        # flags its island as a lower bound; with the switch off it is no
        # network element and flags nothing
        local sc_base = deepcopy(loaded)
        local sm = first(objectsOf(sc_base, :SynchronousMachine))
        delete!(sm.attrs, :ratedS)
        local gu = get(sm.refs, :GeneratingUnit, nothing)
        gu === nothing || delete!(sc_base.objects[gu].attrs, :maxOperatingP)
        local sm_t = only(ce_terms[sm.mrid])
        local sc_pec = deepcopy(sc_base)
        pec!(sc_pec, "_pec_sc", "PEC_SC", sm_t, -10.0, 0.0)
        add!(sc_pec, "_pec_sc_pv", :PhotoVoltaicUnit, Dict{Symbol,String}(:minP => "0", :maxP => "35"), Dict{Symbol,String}(:PowerElectronicsConnection => "_pec_sc"))
        local sc_bus = topo.bus_name[tn_of(sm_t)]
        local sc_row(store, on) = only(run_with_expected_warnings(() -> runShortCircuit!(importCGMES(path = store, name = "sg_pec_sc", map_power_electronics = on); buses = [sc_bus], case = :max), ["short-circuit: "]).rows)
        local row_without, row_with, row_off = sc_row(sc_base, true), sc_row(sc_pec, true), sc_row(sc_pec, false)
        @test row_with.ik_kA == row_without.ik_kA
        @test any(r -> occursin("SynchronousMachine $(str(sm, :name)) has neither x''_d nor ratedS nor a usable maximum active power", r), row_with.reasons)
        @test "short-circuit: PowerElectronicsConnection at bus $(sc_bus): no converter short-circuit model, contribution skipped; Ik''max in its island is a lower bound" in row_with.reasons
        @test !any(r -> occursin("PowerElectronicsConnection", r), row_without.reasons)
        @test (row_off.ik_kA, row_off.reasons) == (row_without.ik_kA, row_without.reasons)
        # with maxIFault, ratedS and ratedU the connection is a current
        # source: at its own bus all of I_sk = maxIFault·ratedS/(√3·ratedU)
        # flows into the fault, its flag is gone, the machine's stays
        local sc_vn = topo.vn_kV[tn_of(sm_t)]
        local sc_ifault = deepcopy(sc_base)
        pec!(sc_ifault, "_pec_sc", "PEC_SC", sm_t, -10.0, 0.0; extra = Dict{Symbol,String}(:maxIFault => "1.1", :ratedS => "30", :ratedU => string(sc_vn)))
        local row_ifault = sc_row(sc_ifault, true)
        @test isapprox(row_ifault.ik_kA - row_without.ik_kA, 1.1 * 30.0 / (sqrt(3.0) * sc_vn); rtol = 1e-9)
        @test !any(r -> occursin("PowerElectronicsConnection", r), row_ifault.reasons)
        @test row_ifault.reasons == row_without.reasons
        @test (row_ifault.zk_ohm, row_ifault.kappa) == (row_without.zk_ohm, row_without.kappa)
        @test (sc_row(sc_ifault, false).ik_kA) == row_without.ik_kA
        # the configured import forwards the switch, and the YAML loader reads it
        @test !SparlectraConfig().cgmes.map_power_electronics
        local pe_cfg = SparlectraConfig(cgmes = Sparlectra.CGMESImportConfig(map_power_electronics = true))
        @test first(importCGMES(pe_cfg; path = pe_store, name = "sg_pec_cfg").element_by_mrid["_pec_pq"]) == :generator
        local yaml_path = joinpath(mktempdir(), "pec.yaml")
        write(yaml_path, "cgmes_import:\n  map_power_electronics: true\n")
        @test Sparlectra.CGMESImportConfig(Sparlectra.load_yaml_dict(yaml_path)).map_power_electronics
      end)() end

      # A SynchronousMachine on a WindGeneratingUnit or SolarGeneratingUnit
      # is a converter unit: in the short
      # circuit a current source like a PowerElectronicsConnection,
      # short_circuit.converter_current_factor times the rated current from
      # ratedS, else from GeneratingUnit.maxOperatingP, else no contribution
      # and a lower bound. SmallGrid's units are ThermalGeneratingUnits; in a
      # copy of the store the unit of one machine changes its class.
      @testset "SynchronousMachine on a wind or solar unit is a converter current source" begin (function ()
        local sm = first(objectsOf(loaded, :SynchronousMachine))
        local gu_mrid = sm.refs[:GeneratingUnit]
        local sm_bus = CI.busOfEquipment(topo, sm, 1).bus
        local function variant(class; ratedS = :keep, maxP = :keep)
          local s = deepcopy(loaded)
          local g = s.objects[gu_mrid]
          local attrs = copy(g.attrs)
          maxP === :keep || (maxP === nothing ? delete!(attrs, :maxOperatingP) : (attrs[:maxOperatingP] = string(maxP)))
          s.objects[gu_mrid] = CI.CIMObject(gu_mrid, class, attrs, copy(g.refs), g.source)
          filter!(!=(gu_mrid), s.byclass[g.class])
          push!(get!(Vector{String}, s.byclass, class), gu_mrid)
          local m = s.objects[sm.mrid]
          ratedS === :keep || (ratedS === nothing ? delete!(m.attrs, :ratedS) : (m.attrs[:ratedS] = string(ratedS)))
          return s
        end
        local imported(s) = importCGMES(path = s, name = "sg_wind")
        local row(res, k) = only(run_with_expected_warnings(() -> runShortCircuit!(res; buses = [sm_bus], case = :max, converter_current_factor = k), ["short-circuit: "]).rows)
        local wind = imported(variant(:WindGeneratingUnit; ratedS = 40.0))
        local rec = only(r for r in wind.shortcircuit.synchronous_machines if r.mrid == sm.mrid)
        @test rec.converter_unit == "WindGeneratingUnit"
        # the rated current refers to ratedU, else to the bus voltage
        local urated = something(rec.ratedU_kV, topo.vn_kV[CI.tnOfTerminal(only(t for t in objectsOf(loaded, :Terminal) if get(t.refs, :ConductingEquipment, "") == sm.mrid))])
        local label = "SynchronousMachine $(str(sm, :name)) on a WindGeneratingUnit"
        # factor off: no converter model, no contribution, a lower bound
        local off = row(wind, 0.0)
        @test off.lower_bound
        @test "short-circuit: $(label) is a converter unit without a converter short-circuit model (short_circuit.converter_current_factor off): contribution skipped; Ik''max in its island is a lower bound" in off.reasons
        # ratedS: k times the rated current, all of it into a fault at its bus
        local with_rated = row(wind, 1.2)
        @test isapprox(with_rated.ik_kA - off.ik_kA, 1.2 * 40.0 / (sqrt(3.0) * urated); rtol = 1e-9)
        @test !with_rated.lower_bound && with_rated.contains_defaulted_data
        @test "short-circuit: $(label) is a converter unit: short_circuit.converter_current_factor 1.2 times the rated current from its ratedS 40.0 MVA substituted" in with_rated.reasons
        @test (with_rated.zk_ohm, with_rated.kappa) == (off.zk_ohm, off.kappa)
        # a converter unit never enters the minimum case
        @test !any(r -> occursin(label, r), only(run_with_expected_warnings(() -> runShortCircuit!(wind; buses = [sm_bus], case = :min, converter_current_factor = 1.2), ["short-circuit: "]).rows).reasons)
        # no ratedS: the maximum active power of the unit (solar the same)
        local solar = imported(variant(:SolarGeneratingUnit; ratedS = nothing, maxP = 25.0))
        @test only(r for r in solar.shortcircuit.synchronous_machines if r.mrid == sm.mrid).pmax_MW == 25.0
        local from_p = row(solar, 1.2)
        @test isapprox(from_p.ik_kA - off.ik_kA, 1.2 * 25.0 / (sqrt(3.0) * urated); rtol = 1e-9)
        @test "short-circuit: SynchronousMachine $(str(sm, :name)) on a SolarGeneratingUnit is a converter unit: short_circuit.converter_current_factor 1.2 times the rated current from its maximum active power 25.0 MW substituted" in from_p.reasons
        # placeholders in both (the values of a delivery with +-99999.9 MW
        # limits): no size, no contribution, a lower bound
        local ph = row(imported(variant(:WindGeneratingUnit; ratedS = sqrt(2.0) * 99999.9, maxP = 99999.9)), 1.2)
        @test ph.ik_kA == off.ik_kA && ph.lower_bound
        @test any(r -> startswith(r, "short-circuit: $(label) is a converter unit without a usable ratedS or maximum active power (ratedS ") && endswith(r, " MVA is a placeholder), so its rated current is unknown: contribution skipped; Ik''max in its island is a lower bound"), ph.reasons)
        # a thermal machine is unchanged: no converter mark, no maximum in
        # the record, and the factor changes no row of the plain import
        local plain = imported(deepcopy(loaded))
        @test all(r -> r.converter_unit === nothing && r.pmax_MW === nothing, plain.shortcircuit.synchronous_machines)
        local all_rows(k) = run_with_expected_warnings(() -> runShortCircuit!(plain; case = :max, converter_current_factor = k), ["short-circuit: "]).rows
        local thermal = all_rows(0.0)
        @test isequal(thermal, all_rows(1.2))
        @test only(r for r in thermal if r.bus == sm_bus).ik_kA > with_rated.ik_kA
        # the SCF file carries the mark in its source row: the same Ik'' after
        # export and import
        local scf_file = joinpath(mktempdir(), "sg_wind.scf.json")
        exportSCF(wind.net; file = scf_file)
        local back = importSCF(scf_file)
        @test only(r for r in back.sc_sources.synchronous_machines if r.mrid == sm.mrid).converter_unit == "WindGeneratingUnit"
        local back_row = only(run_with_expected_warnings(() -> runShortCircuit!(back; buses = [sm_bus], case = :max, converter_current_factor = 1.2), ["short-circuit: "]).rows)
        @test back_row.ik_kA ≈ with_rated.ik_kA rtol = 1e-9
        # the CGMES export keeps the class as delivered: one machine on an
        # offshore WindGeneratingUnit, another
        # on a SolarGeneratingUnit come back with the same class, the same
        # unit mRID and windGenUnitType, and as current sources again
        local sm2 = collect(objectsOf(loaded, :SynchronousMachine))[2]
        local gu2_mrid = sm2.refs[:GeneratingUnit]
        local both = variant(:WindGeneratingUnit; ratedS = 40.0)
        both.objects[gu_mrid].attrs[:windGenUnitType] = "WindGenUnitKind.offshore"
        local g2 = both.objects[gu2_mrid]
        both.objects[gu2_mrid] = CI.CIMObject(gu2_mrid, :SolarGeneratingUnit, copy(g2.attrs), copy(g2.refs), g2.source)
        filter!(!=(gu2_mrid), both.byclass[g2.class])
        push!(get!(Vector{String}, both.byclass, :SolarGeneratingUnit), gu2_mrid)
        local first_import = imported(both)
        local cgmes_dir = mktempdir()
        writeCGMESFiles(first_import.net; path = cgmes_dir, notices = String[])
        local reread = loadCGMES(cgmes_dir)
        local unit_of(m) = reread.objects[reread.objects[m].refs[:GeneratingUnit]]
        @test (unit_of(sm.mrid).mrid, unit_of(sm.mrid).class) == (gu_mrid, :WindGeneratingUnit)
        @test str(unit_of(sm.mrid), :windGenUnitType) == "WindGenUnitKind.offshore"
        @test (unit_of(sm2.mrid).mrid, unit_of(sm2.mrid).class) == (gu2_mrid, :SolarGeneratingUnit)
        local cgmes_back = importCGMES(path = reread, name = "sg_wind_back")
        @test Dict(r.mrid => r.converter_unit for r in cgmes_back.shortcircuit.synchronous_machines if r.converter_unit !== nothing) == Dict(sm.mrid => "WindGeneratingUnit", sm2.mrid => "SolarGeneratingUnit")
        local cgmes_row = only(run_with_expected_warnings(() -> runShortCircuit!(cgmes_back; buses = [sm_bus], case = :max, converter_current_factor = 1.2), ["short-circuit: "]).rows)
        # the export names the machine after its network element, so the
        # reason is matched without the name
        @test any(r -> endswith(r, " on a WindGeneratingUnit is a converter unit: short_circuit.converter_current_factor 1.2 times the rated current from its ratedS 40.0 MVA substituted"), cgmes_row.reasons)
        @test cgmes_row.ik_kA ≈ row(first_import, 1.2).ik_kA rtol = 1e-6
      end)() end

      # bank regulators behind cgmes_import.shunt_control:
      # the SmallGrid carries 14 LinearShuntCompensators, each with SSH
      # controlEnabled and an enabled voltage RegulatingControl at its own
      # terminal; W.Lancst SC sits on a bus a machine holds PV
      @testset "bank regulators behind cgmes_import.shunt_control" begin (function ()
        local plain = importCGMES(path = loaded, name = "sg_banks_off")
        @test isempty(Sparlectra._shunt_controllers(plain.net)) && !any(m -> occursin("shunt control", m), plain.messages)
        local on = importCGMES(path = loaded, name = "sg_banks", shunt_control = true)
        local banks = Sparlectra._shunt_controllers(on.net)
        # the regulators take over the imported shunts, none is added
        @test length(banks) == 13 && length(on.net.shuntVec) == length(plain.net.shuntVec) == 14
        @test sort!([b.shunt_idx for b in banks]) == setdiff(1:14, [5])
        @test "notice: shunt control: W.Lancst SC stays a fixed shunt, bus W.Lancst is PV, its voltage is held otherwise" in on.messages
        @test "notice: cgmes_import.shunt_control = true: 13 bank regulator(s) attached, 1 controlled bank(s) kept as fixed shunts" in on.messages
        # Sprigg SC: bPerSection 1.148e-4 S at 132 kV, 10 sections, 5 in, target
        # 129.935 kV, deadband 1 kV
        local sprigg = only(b for b in banks if b.name == "Sprigg SC")
        local blk = 1.148e-4 * 132.0^2
        @test (sprigg.bus, sprigg.step_mvar, sprigg.bs_min_mvar, sprigg.bs_max_mvar) == ("Sprigg", blk, 0.0, 10.0 * blk)
        @test sprigg.bs_mvar ≈ 5.0 * blk && sprigg.target_vm_pu == 129.935 / 132.0 && sprigg.deadband_vm_pu == 1.0 / 132.0 / 2.0
        # the taken-over shunt keeps its admittance
        @test on.net.shuntVec[sprigg.shunt_idx].y_pu_shunt == plain.net.shuntVec[sprigg.shunt_idx].y_pu_shunt
        # the run with the controllers converges, every bank inside its band
        local res = redirect_stdout(() -> run_sparlectra(net = on.net, config = SparlectraConfig()), devnull)
        @test res.final_converged && latest_control_result(on.net).status === :converged
        @test all(b -> abs(get_bus_vm_pu(on.net, b.bus) - b.target_vm_pu) <= b.deadband_vm_pu, banks)
        # the configured import forwards the switch, the YAML loader reads it
        @test !SparlectraConfig().cgmes.shunt_control
        @test length(Sparlectra._shunt_controllers(importCGMES(SparlectraConfig(cgmes = Sparlectra.CGMESImportConfig(shunt_control = true)); path = loaded, name = "sg_banks_cfg").net)) == 13
        @test Sparlectra.CGMESImportConfig(Dict{String,Any}("cgmes_import" => Dict{String,Any}("shunt_control" => true))).shunt_control
        @test Sparlectra.cgmes_adapter_options(SparlectraConfig(cgmes = Sparlectra.CGMESImportConfig(shunt_control = true))).shunt_control
        # addShuntVoltageControl! with shunt_idx refuses a shunt it cannot take over
        local net = deepcopy(plain.net)
        local olive = net.shuntVec[3]
        local olive_bus = only(k for (k, v) in net.busDict if v == olive.busIdx)
        local other = only(k for (k, v) in net.busDict if v == net.shuntVec[1].busIdx)
        local add(; kw...) = Sparlectra.addShuntVoltageControl!(net; target_vm_pu = 1.0, bs_min_mvar = 0.0, bs_max_mvar = 80.0, step_mvar = 8.0, kw...)
        @test_throws "shunt_idx 99 out of range" add(bus = olive_bus, shunt_idx = 99)
        @test_throws "is not connected to bus" add(bus = other, shunt_idx = 3)
        net.shuntVec[3].y_pu_shunt = ComplexF64(0.01, imag(olive.y_pu_shunt))
        @test_throws "has a conductance" add(bus = olive_bus, shunt_idx = 3)
        net.shuntVec[3].y_pu_shunt = ComplexF64(0.0, imag(olive.y_pu_shunt))
        # a start between two blocks is snapped and written to the shunt
        local c = add(bus = olive_bus, shunt_idx = 3)
        @test c.shunt_idx == 3 && c.bs_mvar == 40.0 && imag(net.shuntVec[3].y_pu_shunt) == 40.0 / net.baseMVA && length(net.shuntVec) == 14
      end)() end

      # the DY profile behind cgmes_import.read_dynamics (#506), on the
      # store copy of _cgmes_dy_store and on a DY file next to the SmallGrid
      # zips: off, nothing is read and the network is the plain one; on,
      # inertia, damping, rating and droop land on the right machines, the
      # report counts every DY class, and the data survive SCF
      @testset "dynamic data from the DY profile" begin (function ()
        local plain = importCGMES(path = loaded, name = "sg_dy_plain")
        local dys = _cgmes_dy_store(loaded)
        local ps_of(r, name) = r.net.prosumpsVec[last(r.element_by_mrid[only(o for o in objectsOf(loaded, :SynchronousMachine) if str(o, :name) == name).mrid])]
        # no DY profile: no DY summary, no report line
        @test plain.dynamics === nothing && isempty(dynamicsReportLines(plain))
        # off: no DynamicData, one DY line naming the switch (also the whole
        # report section), the network of the plain import
        local off = importCGMES(path = dys, name = "sg_dy_off")
        @test all(ps -> ps.dynamics === nothing, off.net.prosumpsVec)
        local off_line = "DY profile present (1 file(s)), not read: cgmes_import.read_dynamics = false"
        @test filter(m -> startswith(m, "notice: DY"), off.messages) == ["notice: " * off_line]
        @test dynamicsReportLines(off) == [off_line]
        @test isempty(scf_roundtrip_field_diffs(plain.net, off.net))
        # on: TwinBrch carries inertia, damping, its EQ rating and the droop
        # of GovSteam1 (1/k on mwbase 100, converted to 200 MVA) with its t1
        local on = importCGMES(path = dys, name = "sg_dy_on", read_dynamics = true)
        local tw = ps_of(on, "TwinBrch SM").dynamics
        @test tw === DynamicData(inertia_s = 4.25, damping_pu = 0.5, rating_mva = 200.0, droop_pu = (1 / 25) * 200.0 / 100.0, governor_time_s = 0.4, unit_type = SynchronousUnit)
        # Breed: inertia without damping; the droop of GovSteam0 is left out
        # (mwbase 0), and its t1 is a steam bowl, not the governor lag
        @test ps_of(on, "Breed SM").dynamics === DynamicData(inertia_s = 6.0, rating_mva = 750.0)
        @test ps_of(on, "DeerCrk SM").dynamics === nothing
        @test count(ps -> ps.dynamics !== nothing, on.net.prosumpsVec) == 2
        # nothing else of the network changed
        @test all(m -> startswith(m, "prosumer.dynamics"), scf_roundtrip_field_diffs(plain.net, on.net))
        for line in ("notice: DY profile read (cgmes_import.read_dynamics): 2 unit(s) with dynamic data",
                     "notice: DY class SynchronousMachineTimeConstantReactance: 1 mapped, 1 not mapped, 0 without counterpart in EQ",
                     "notice: DY class SynchronousMachineSimplified: 1 mapped, 0 not mapped, 0 without counterpart in EQ",
                     "notice: DY class SynchronousMachineEquivalentCircuit: 0 mapped, 0 not mapped, 1 without counterpart in EQ",
                     "notice: DY class GovSteam1: 1 mapped, 0 not mapped, 0 without counterpart in EQ",
                     "notice: DY class GovSteam0: 0 mapped, 1 not mapped, 0 without counterpart in EQ",
                     "notice: DY class GovHydro1: 0 mapped, 0 not mapped, 1 without counterpart in EQ",
                     "notice: DY class ExcIEEEAC4A: 0 mapped, 1 not mapped, 0 without counterpart in EQ",
                     "notice: DY class LoadAggregate: 0 mapped, 1 not mapped, 1 without counterpart in EQ",
                     "notice: DY class LoadStatic: 0 mapped, 1 not mapped, 0 without counterpart in EQ",
                     "notice: DY SynchronousMachineTimeConstantReactance: 1 object(s), function block disabled",
                     "notice: DY GovSteam0: 1 object(s), droop left out: mwbase 0 or absent, the base of the droop is not stated",
                     "notice: DY GovSteam0: 1 object(s), neither a usable droop nor a governor time constant")
          @test line in on.messages
        end
        @test count(m -> startswith(m, "notice: DY class "), on.messages) == 9
        # the report section: the notices of the import, the DY file and the
        # two units under their SynchronousMachine names
        local rep = dynamicsReportLines(on)
        @test filter(l -> !startswith(l, "DY files: ") && !startswith(l, "units with dynamic data: "), rep) == [m[9:end] for m in on.messages if startswith(m, "notice: DY")]
        @test "DY files: dy_test.xml" in rep
        @test Set(split(only(l for l in rep if startswith(l, "units with dynamic data: "))[26:end], ", ")) == Set(["TwinBrch SM", "Breed SM"])
        # a record against the DynamicData rules is dropped with a warning
        local neg = deepcopy(dys)
        neg.objects["_dy_m2"].attrs[:inertia] = "-6"
        local on_neg = @test_logs (:warn, r"^CGMES import: 1 data defect") importCGMES(path = neg, name = "sg_dy_neg", read_dynamics = true)
        @test ps_of(on_neg, "Breed SM").dynamics === nothing
        @test "warning: DY data of $(Sparlectra.getCompName(ps_of(on_neg, "Breed SM").comp)) not taken: DynamicDataError: inertia_s = -6.0 must be >= 0" in on_neg.messages
        @test "notice: DY class SynchronousMachineSimplified: 0 mapped, 1 not mapped, 0 without counterpart in EQ" in on_neg.messages
        # a machine model without inertia still makes the unit synchronous
        # (the keyword default would say grid-following), and GovHydro1 on
        # it gives rperm 0.04 on mwbase 100 converted to 750 MVA: 0.3 pu
        local hyd = deepcopy(dys)
        delete!(hyd.objects["_dy_m2"].attrs, :inertia)
        hyd.objects["_dy_g4"].refs[:SynchronousMachineDynamics] = "_dy_m2"
        local on_hyd = importCGMES(path = hyd, name = "sg_dy_hyd", read_dynamics = true)
        @test ps_of(on_hyd, "Breed SM").dynamics === DynamicData(rating_mva = 750.0, droop_pu = 0.04 * 750.0 / 100.0, unit_type = SynchronousUnit)
        @test "notice: DY class GovHydro1: 1 mapped, 0 not mapped, 0 without counterpart in EQ" in on_hyd.messages
        # CGMES with DY to SCF and back: every value bit for bit
        local f = exportSCF(on.net; file = joinpath(mktempdir(), "sg_dy.scf.json"))
        local back = importSCF(f)
        @test [ps.dynamics for ps in back.prosumpsVec] == [ps.dynamics for ps in on.net.prosumpsVec]
        @test all(k -> back.prosumpsVec[k].dynamics === on.net.prosumpsVec[k].dynamics, eachindex(on.net.prosumpsVec))
        # from files: off, the DY file is skipped by name; on, it is read
        local dir = mktempdir()
        local tw_mrid = only(o for o in objectsOf(loaded, :SynchronousMachine) if str(o, :name) == "TwinBrch SM").mrid
        local ec_mrid = first(objectsOf(loaded, :EnergyConsumer)).mrid
        write(joinpath(dir, "sg_DY.xml"), """<?xml version="1.0" encoding="UTF-8"?>
        <rdf:RDF xmlns:cim="http://iec.ch/TC57/2013/CIM-schema-cim16#" xmlns:md="http://iec.ch/TC57/61970-552/ModelDescription/1#" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
          <md:FullModel rdf:about="urn:uuid:5f0e4c1a-0000-4000-8000-000000000506">
            <md:Model.profile>http://entsoe.eu/CIM/Dynamics/3/1</md:Model.profile>
          </md:FullModel>
          <cim:SynchronousMachineTimeConstantReactance rdf:ID="_dyf_m1">
            <cim:DynamicsFunctionBlock.enabled>true</cim:DynamicsFunctionBlock.enabled>
            <cim:RotatingMachineDynamics.inertia>4.25</cim:RotatingMachineDynamics.inertia>
            <cim:RotatingMachineDynamics.damping>0.5</cim:RotatingMachineDynamics.damping>
            <cim:SynchronousMachineDynamics.SynchronousMachine rdf:resource="#$(tw_mrid)"/>
          </cim:SynchronousMachineTimeConstantReactance>
          <cim:EnergyConsumer rdf:about="#$(ec_mrid)">
            <cim:EnergyConsumer.LoadDynamics rdf:resource="#_dyf_la"/>
          </cim:EnergyConsumer>
          <cim:LoadAggregate rdf:ID="_dyf_la">
            <cim:IdentifiedObject.name>LA</cim:IdentifiedObject.name>
          </cim:LoadAggregate>
        </rdf:RDF>
        """)
        local paths = vcat(sort!(filter(f -> endswith(f, ".zip"), readdir(joinpath(@__DIR__, "fixtures", "cgmes_smallgrid"); join = true))), [joinpath(dir, "sg_DY.xml")])
        local f_off = importCGMES(path = paths, name = "sg_dyf_off")
        @test count(m -> startswith(m, "file skipped: sg_DY.xml") && endswith(m, "profile not in import set"), f_off.messages) == 1
        @test !any(o -> o.source == "sg_DY.xml", values(f_off.store.objects))
        @test dynamicsReportLines(f_off) == ["DY profile present (1 file(s)), not read: cgmes_import.read_dynamics = false"]
        @test isempty(scf_roundtrip_field_diffs(plain.net, f_off.net))
        local f_on = importCGMES(path = paths, name = "sg_dyf_on", read_dynamics = true)
        @test ps_of(f_on, "TwinBrch SM").dynamics === DynamicData(inertia_s = 4.25, damping_pu = 0.5, rating_mva = 200.0)
        @test "notice: DY class LoadAggregate: 0 mapped, 1 not mapped, 0 without counterpart in EQ" in f_on.messages
        @test all(m -> startswith(m, "prosumer.dynamics"), scf_roundtrip_field_diffs(plain.net, f_on.net))
        # the service run writes the same lines into the DY section of
        # cgmes_import_report.md: with the switch the counts, without it the
        # one line
        for (rd, expected) in ((true, dynamicsReportLines(f_on)), (false, dynamicsReportLines(f_off)))
          local out = mktempdir()
          local cfg = joinpath(out, "c.yaml")
          write(cfg, "config_version: 1\ncgmes_import:\n  path: \"" * paths[2] * ";" * paths[3] * "\"\n  read_dynamics: " * string(rd) * "\n")
          local r = run_sparlectra_api(casefile = paths[1], config_file = cfg, output_dir = out, case_format = :cgmes)
          @test r.status == :succeeded
          local text = cgmes_report_plain(out)
          @test occursin("\nDynamic data (DY profile)\n" * join(expected, "\n") * "\nShort-circuit source data", text)
        end
        # the configured import, the YAML loader and the adapter options
        # forward the switch
        @test !SparlectraConfig().cgmes.read_dynamics
        @test ps_of(importCGMES(SparlectraConfig(cgmes = Sparlectra.CGMESImportConfig(read_dynamics = true)); path = dys, name = "sg_dy_cfg"), "TwinBrch SM").dynamics === tw
        @test Sparlectra.CGMESImportConfig(Dict{String,Any}("cgmes_import" => Dict{String,Any}("read_dynamics" => true))).read_dynamics
        @test Sparlectra.cgmes_adapter_options(SparlectraConfig(cgmes = Sparlectra.CGMESImportConfig(read_dynamics = true))).read_dynamics
      end)() end

      # The DY profile of the ENTSO-E MicroGrid (CGMES cache, not in the
      # repository): BE_DY_2.xml holds two and NL_DY_2.xml three
      # SynchronousMachineTimeConstantReactance records, each on a machine of
      # the network, so two and three units get DynamicData; without the
      # switch the report has the one line naming it
      for (alias, n_units) in (("microgrid_be", 2), ("microgrid_nl", 3))
        local zip = joinpath(Sparlectra.CGMESImporter.cgmesTestSetCacheDir(), "cases", "cgmes_" * alias * ".zip")
        if !isfile(zip)
          println("      cgmes_importer: ", alias, " DY profile SKIPPED (", zip, " not in the CGMES cache)")
          continue
        end
        println("      cgmes_importer: ", alias, " DY profile RAN")
        @testset "DY profile of $(alias)" begin (function ()
          local on = importCGMES(path = [zip], name = alias, read_dynamics = true)
          @test count(ps -> ps.dynamics !== nothing, on.net.prosumpsVec) == n_units
          @test length(on.dynamics.units) == n_units && length(on.dynamics.files) == 1
          @test first(dynamicsReportLines(on)) == "DY profile read (cgmes_import.read_dynamics): $(n_units) unit(s) with dynamic data"
          @test "DY class SynchronousMachineTimeConstantReactance: $(n_units) mapped, 0 not mapped, 0 without counterpart in EQ" in dynamicsReportLines(on)
          @test dynamicsReportLines(importCGMES(path = [zip], name = alias)) == ["DY profile present (1 file(s)), not read: cgmes_import.read_dynamics = false"]
        end)() end
      end

      # Two boundary sets next to a base case (the Web UI case directory with
      # the MicroGrid and the MiniGrid test configurations): the run pairs the
      # one that resolves the case. By name alone the MicroGrid BD sorts
      # first and the MiniGrid import failed with unresolved topology
      # references.
      local ext = joinpath(Sparlectra.CGMESImporter.cgmesTestSetCacheDir(), "extracted")
      local mini_case = joinpath(ext, "MiniGrid", "NodeBreaker", "CGMES_v2.4.15_MiniGridTestConfiguration_BaseCase_Complete_v3")
      local mini_bd = joinpath(ext, "MiniGrid", "NodeBreaker", "CGMES_v2.4.15_MiniGridTestConfiguration_Boundary_v3")
      local micro_bd = joinpath(ext, "MicroGrid", "BaseCase_BC", "CGMES_v2.4.15_MicroGridTestConfiguration_BD_v2")
      if all(isdir, (mini_case, mini_bd, micro_bd))
        println("      cgmes_importer: boundary pairing with two boundary sets RAN")
        @testset "boundary pairing picks the set that resolves the case" begin (function ()
          local dir = mktempdir()
          for p in (mini_case, mini_bd, micro_bd)
            symlink(p, joinpath(dir, basename(p)))
          end
          local paths, auto = Sparlectra._cgmes_delivery_paths(joinpath(dir, basename(mini_case)), Sparlectra.CGMESImportConfig())
          @test auto
          @test basename.(paths) == [basename(mini_case), basename(mini_bd)]
        end)() end
      else
        println("      cgmes_importer: boundary pairing with two boundary sets SKIPPED (MiniGrid or MicroGrid boundary not in the CGMES cache)")
      end

      # cgmes_import.equivalent_injection_as_source (default off): the
      # configured import, the YAML loader and the adapter options forward it
      @test !SparlectraConfig().cgmes.equivalent_injection_as_source
      @test Sparlectra.CGMESImportConfig(Dict{String,Any}("cgmes_import" => Dict{String,Any}("equivalent_injection_as_source" => true))).equivalent_injection_as_source
      @test Sparlectra.cgmes_adapter_options(SparlectraConfig(cgmes = Sparlectra.CGMESImportConfig(equivalent_injection_as_source = true))).equivalent_injection_as_source

      # MicroGrid BE alone (CGMES cache, not in the repository): no
      # referencePriority, no ExternalNetworkInjection, and the SV names the
      # X-node TN_Border_ST23 as angle reference of its TopologicalIsland
      # "Slack". The node carries only the EquivalentInjection BE-Inj-XZE_ST23
      # (r = x = 0 in the EQ), so the slack sits there, at its SV voltage and
      # angle (-5.58021 deg), and not at the largest machine BE-G1.
      local be_zip = joinpath(Sparlectra.CGMESImporter.cgmesTestSetCacheDir(), "cases", "cgmes_microgrid_be.zip")
      if !isfile(be_zip)
        println("      cgmes_importer: microgrid_be angle reference from SV SKIPPED (", be_zip, " not in the CGMES cache)")
      else
        println("      cgmes_importer: microgrid_be angle reference from SV RAN")
        @testset "MicroGrid BE: angle reference from SV, EquivalentInjection as source" begin (function ()
          local x_ei = "_87ea56f3-962a-427a-85d6-13b1f9295174"
          local x_va = -5.580210
          local ps_x(r) = r.net.prosumpsVec[last(r.element_by_mrid[x_ei])]
          local node_type(r, bus) = r.net.nodeVec[r.net.busDict[bus]]._nodeType
          # a copy of the delivery: without the files `drop` names, with the
          # EQ text changed by `edit`
          local repack(dst; drop = n -> false, edit = (n, s) -> s) = begin
            local rd = ZipArchives.ZipReader(read(be_zip))
            ZipArchives.ZipWriter(dst) do w
              for n in ZipArchives.zip_names(rd)
                drop(n) && continue
                ZipArchives.zip_newfile(w, n)
                write(w, edit(n, String(ZipArchives.zip_readentry(rd, n))))
              end
            end
            dst
          end
          # bands of the SV comparison: max |dVm| was 8.08e-5 pu with the
          # slack at BE-G1 (0.60.0 before this change) and is 7.56e-5 pu with
          # the slack on the X-node. The reference bus sits exactly at its SV
          # angle. The raw angle offset (compareWithSV, median of the angle
          # differences) is -1.42e-3 deg with the slack on the X-node and was
          # -5.9e-4 deg with the slack at BE-G1: it is not zero, so the band
          # below bounds it and the exact check is on the reference bus
          local check_sv(r) = begin
            @test last(runpf!(r.net, 40, 1e-8, 0)) == 0
            local cmp = compareWithSV(r)
            @test cmp.n == 11
            @test cmp.max_dvm < 8.1e-5
            @test abs(cmp.va_ref_offset_deg) < 5e-3
            @test abs(only(x for x in cmp.rows if x.bus == "TN_Border_ST23").dva) < 1e-9
            cmp
          end

          # 4a: switch off
          local off = importCGMES(path = [be_zip], name = "be_off")
          @test off.slack_bus == "TN_Border_ST23"
          @test Sparlectra.isSlack(ps_x(off)) && ps_x(off).va_deg == x_va
          @test "notice: island Slack: angle reference TN_Border_ST23 taken from SV" in off.messages
          @test !any(m -> startswith(m, "warning:"), off.messages)
          # the other injections keep their own bus type: BE-G1 regulates
          # its bus and stays PV, the other X-nodes stay PQ
          @test node_type(off, "BE-Busbar_4") == Sparlectra.PV
          @test count(nd -> nd._nodeType == Sparlectra.Slack, off.net.nodeVec) == 1
          check_sv(off)

          # 2b on a second SV TopologicalIsland "Second" whose angle
          # reference is the node of BE-G1, in the same electrical island: the
          # machine outranks the EquivalentInjection, and the line names the
          # holder's island and the island of the reference not promoted
          local two = loadCGMES(be_zip)
          local g1_tn = only(o.mrid for o in objectsOf(two, :TopologicalNode) if str(o, :name) == "BE-Busbar_4")
          two.objects["_second_isl"] = Sparlectra.CGMESImporter.CIMObject("_second_isl", :TopologicalIsland, Dict{Symbol,String}(:name => "Second"), Dict{Symbol,String}(:AngleRefTopologicalNode => g1_tn), "sv_test")
          push!(get!(Vector{String}, two.byclass, :TopologicalIsland), "_second_isl")
          local r2 = importCGMES(path = two, name = "be_two_refs")
          @test r2.slack_bus == "BE-Busbar_4"
          @test count(m -> m == "notice: island Second: angle reference BE-Busbar_4, other injections keep their own bus type; the SV angle reference TN_Border_ST23 of island Slack lies in the same electrical island and is not promoted to slack", r2.messages) == 1
          @test count(nd -> nd._nodeType == Sparlectra.Slack, r2.net.nodeVec) == 1

          # 4b: switch on, every X-node of the delivery carries a source at
          # its SV voltage and angle; the reference stays where it was
          local on = importCGMES(path = [be_zip], name = "be_on", equivalent_injection_as_source = true)
          @test on.slack_bus == "TN_Border_ST23"
          local srcs = [on.net.prosumpsVec[k] for (m, (kind, k)) in on.element_by_mrid if kind == :external_grid]
          @test length(srcs) == 5 && all(p -> p.isRegulated && p.comp.cTyp == Sparlectra.ExternalNetworkInjection, srcs)
          local by_bus = Dict(nd.busIdx => nd for nd in on.net.nodeVec)
          @test all(p -> p.vm_pu == by_bus[Sparlectra.getPosumerBusIndex(p)]._vm_pu && p.va_deg == by_bus[Sparlectra.getPosumerBusIndex(p)]._va_deg, srcs)
          @test ps_x(on).va_deg == x_va
          @test any(m -> startswith(m, "notice: 5 EquivalentInjection(s) imported as source at the SV voltage of their node"), on.messages)
          @test any(m -> startswith(m, "notice: 5 EquivalentInjection(s) imported as source without internal impedance") && occursin("BE-Inj-XZE_ST23", m), on.messages)
          check_sv(on)

          # A6: the X-node injection with r = 4.84 ohm, x = 48.4 ohm at
          # Un = 220 kV: |Z| = 48.4 * sqrt(1.01) ohm, so the SCF source states
          # sk = Un^2/|Z| = 220e3^2 / (48.4 * sqrt(1.01)) = 1e9/sqrt(1.01) VA
          # (995.04 MVA) and rx_ratio = r/x = 0.1
          local d = mktempdir()
          local rx_zip = repack(joinpath(d, "be_rx.zip"); edit = (n, s) -> (occursin("_EQ_", n) && !occursin("_BD_", n)) ? replace(s, r"(<cim:EquivalentInjection rdf:ID=\"_87ea56f3-962a-427a-85d6-13b1f9295174\">.*?<cim:EquivalentInjection.r>)0e\+000(</cim:EquivalentInjection.r>\s*<cim:EquivalentInjection.x>)0e\+000"s => s"\g<1>4.84\g<2>48.4") : s)
          local rx = importCGMES(path = [rx_zip], name = "be_rx", equivalent_injection_as_source = true)
          @test any(m -> startswith(m, "notice: 4 EquivalentInjection(s) imported as source without internal impedance") && !occursin("BE-Inj-XZE_ST23", m), rx.messages)
          local a = exportSCF(rx.net; file = joinpath(d, "a.scf.json"), source_format = "cgmes")
          local srow = only(r for r in Sparlectra.scf_json_parse(read(a, String))["data"]["source"] if haskey(r, "sk"))
          @test isapprox(srow["sk"], 1e9 / sqrt(1.01); rtol = 1e-12)
          @test isapprox(srow["rx_ratio"], 0.1; rtol = 1e-12)
          # power_flow.external_grid (source auto) takes the same feeder for
          # the slack from the EquivalentInjection record: 1000/sqrt(1.01) MVA
          local fd = Sparlectra._declared_slack_feeder(rx.net, rx.net.sc_sources)
          @test fd !== nothing && isapprox(fd.sk_MVA, 1e3 / sqrt(1.01); rtol = 1e-12) && isapprox(fd.rx, 0.1; rtol = 1e-12)
          @test srow["u_ref"] == ps_x(rx).vm_pu && rad2deg(srow["u_ref_angle"]) == x_va
          # 4d: the SCF round trip is a fixed point and keeps the one record
          # of the injection: the reader takes it from the sc_source block
          # and adds no second feeder record from the source's sk
          local back = importSCF(a)
          @test isempty(back.sc_sources.external_network_injections)
          @test count(e -> e.bus == "TN_Border_ST23", back.sc_sources.equivalent_injections) == 1
          local fb = Sparlectra._declared_slack_feeder(back, back.sc_sources)
          @test fb !== nothing && isapprox(fb.sk_MVA, 1e3 / sqrt(1.01); rtol = 1e-12)
          local b = exportSCF(back; file = joinpath(d, "b.scf.json"), source_format = "cgmes")
          @test read(b, String) == read(exportSCF(importSCF(b); file = joinpath(d, "c.scf.json"), source_format = "cgmes"), String)

          # 4c: the same delivery without its SV file has no angle reference
          # in SV or SSH: the largest machine takes the slack, and the
          # notice names the electrical island
          local nosv = importCGMES(path = [repack(joinpath(d, "be_nosv.zip"); drop = n -> occursin("_SV_", n))], name = "be_nosv")
          @test nosv.slack_bus == "BE-Busbar_4"
          @test any(m -> occursin(r"^notice: electrical island \d+: no angle reference in SV or SSH, largest machine BE-G1 at BE-Busbar_4 as slack$", m), nosv.messages)
        end)() end
      end

      # The short-circuit harvest keeps the sources the import mapped. A
      # machine of SmallGrid without ratedS and maxOperatingP is the live
      # source at its bus (its row rests on defaults). A machine with a
      # disconnected terminal and a network injection out of service, both
      # with full short-circuit data on that bus, are skipped and leave no
      # record, so Ik'' and the reasons there stay as without them; their
      # records, put back, would feed that bus.
      @testset "short-circuit harvest keeps only mapped sources" begin (function ()
        local add!(s, mrid, class, attrs, refs) = (s.objects[mrid] = CI.CIMObject(mrid, class, attrs, refs, "sc_test"); push!(get!(Vector{String}, s.byclass, class), mrid))
        local tn_of(t) = t.refs[:TopologicalNode]
        local sc_base = deepcopy(loaded)
        local sm = first(objectsOf(sc_base, :SynchronousMachine))
        delete!(sm.attrs, :ratedS)
        local gu = get(sm.refs, :GeneratingUnit, nothing)
        gu === nothing || delete!(sc_base.objects[gu].attrs, :maxOperatingP)
        local sm_t = only(t for t in objectsOf(sc_base, :Terminal) if get(t.refs, :ConductingEquipment, "") == sm.mrid)
        local sc_bus = topo.bus_name[tn_of(sm_t)]
        local row_without = only(run_with_expected_warnings(() -> runShortCircuit!(importCGMES(path = sc_base, name = "sg_sc_base"); buses = [sc_bus], case = :max), ["short-circuit: "]).rows)
        local sc_skip = deepcopy(sc_base)
        add!(sc_skip, "_sm_off", :SynchronousMachine, Dict{Symbol,String}(:name => "SM_OFF", :p => "-10", :q => "0", :satDirectSubtransX => "0.2", :ratedS => "100", :ratedU => string(topo.vn_kV[tn_of(sm_t)])), Dict{Symbol,String}())
        add!(sc_skip, "_sm_off_t", :Terminal, Dict{Symbol,String}(:sequenceNumber => "1", :connected => "false"), Dict{Symbol,String}(:ConductingEquipment => "_sm_off", :TopologicalNode => tn_of(sm_t)))
        add!(sc_skip, "_eni_off", :ExternalNetworkInjection, Dict{Symbol,String}(:name => "ENI_OFF", :p => "0", :q => "0", :inService => "false", :maxInitialSymShCCurrent => "20000", :maxR1ToX1Ratio => "0.1"), Dict{Symbol,String}())
        add!(sc_skip, "_eni_off_t", :Terminal, Dict{Symbol,String}(:sequenceNumber => "1", :connected => "true"), Dict{Symbol,String}(:ConductingEquipment => "_eni_off", :TopologicalNode => tn_of(sm_t)))
        local skip_res = importCGMES(path = sc_skip, name = "sg_sc_skip")
        @test "_sm_off" in skip_res.skipped_equipment && "_eni_off" in skip_res.skipped_equipment
        @test !any(r -> r.mrid in ("_sm_off", "_eni_off"), vcat(skip_res.shortcircuit.synchronous_machines, skip_res.shortcircuit.external_network_injections))
        @test skip_res.net.sc_sources.synchronous_machines == skip_res.shortcircuit.synchronous_machines
        @test skip_res.net.sc_sources.external_network_injections == skip_res.shortcircuit.external_network_injections
        local row_skip = only(run_with_expected_warnings(() -> runShortCircuit!(skip_res; buses = [sc_bus], case = :max), ["short-circuit: "]).rows)
        @test (row_skip.ik_kA, row_skip.reasons) == (row_without.ik_kA, row_without.reasons)
        local all_recs = CI.collectShortCircuitData(sc_skip, skip_res.topo)
        local sc_back = CI.CGMESShortCircuitData(all_recs.external_network_injections, all_recs.synchronous_machines, skip_res.shortcircuit.ac_line_segments, skip_res.shortcircuit.transformer_ends, skip_res.shortcircuit.equivalent_injections, skip_res.shortcircuit.asynchronous_machines, skip_res.shortcircuit.power_electronics_connections)
        @test any(r -> r.mrid == "_sm_off" && r.bus == sc_bus, sc_back.synchronous_machines)
        @test only(run_with_expected_warnings(() -> runShortCircuit!(skip_res.net, sc_back; buses = [sc_bus], case = :max), ["short-circuit: "]).rows).ik_kA > row_without.ik_kA + 1.0
      end)() end

      # islands without a slack, on new nodes next to SmallGrid
      # (_cgmes_f2_island_store); SmallGrid itself has none
      @testset "islands without a slack" begin (function ()
        # one logged warning collects the warning: messages of the import
        local f2 = @test_logs (:warn, r"^CGMES import: 4 data defect") importCGMES(path = _cgmes_f2_island_store(loaded), name = "sg_f2")
        local f2_bus(tn) = f2.net.busDict[f2.topo.bus_name[tn]]
        local f2_ps(mrid) = f2.net.prosumpsVec[last(f2.element_by_mrid[mrid])]
        # A, F and H have no source: A a load and an injection with P = 0 and
        # no regulation, F a load and an injection that draws power, H an
        # injection that draws power and a load that feeds in (a load is no
        # source). All are de-energized, voltage 0, load unsupplied; on F and
        # H the unsupplied P includes the power the injection drew, and the
        # 10 MW the load on H fed in are reported apart, not netted
        for (tns, units, line, tail) in ((("_f2_a1", "_f2_a2"), ("_f2_ei_a", "_f2_load_a"), "_f2_la", "load P = 11.3 MW unsupplied"),
                                         (("_f2_f1", "_f2_f2"), ("_f2_ei_f", "_f2_load_f"), "_f2_lf", "load P = 7.0 MW unsupplied"),
                                         (("_f2_h1", "_f2_h2"), ("_f2_ei_h", "_f2_load_h"), "_f2_lh", "load P = 8.0 MW unsupplied, load infeed P = 10.0 MW dropped"))
          for tn in tns
            @test f2_bus(tn) in f2.net.isoNodes
            @test f2.net.nodeVec[f2_bus(tn)]._vm_pu == 0.0
          end
          @test all(m -> (f2_ps(m).pVal, f2_ps(m).qVal) == (0.0, 0.0), units)
          @test !Sparlectra.isSlack(f2_ps(first(units)))
          @test f2.net.branchVec[last(f2.element_by_mrid[line])].status == 0
          @test count(m -> startswith(m, "notice: de-energized a component of 2 bus(es) around ") && endswith(m, "2 injection(s) zeroed, $(tail)"), f2.messages) == 1
        end
        @test count(m -> startswith(m, "notice: de-energized a component"), f2.messages) == 3
        @test "notice: 3 component(s) without source de-energized: 6 bus(es), load P = 26.3 MW unsupplied, load infeed P = 10.0 MW dropped" in f2.messages

        # B to E: infeed but no slack and no voltage-regulating unit; the
        # import makes one unit with infeed the slack, at the start (SV)
        # voltage of its bus, with one warning per island and one notice with
        # the count. A unit that draws power is no candidate: not the machine
        # with the stronger referencePriority on B, not the network injection
        # with the larger reserve on D, not the injection on the inner bus of
        # E, where the infeed on the degree-1 bus wins
        local f2_warn(mrid, reason) = count(m -> startswith(m, "warning: electrical island of ") && occursin("has infeed but no reference; $(uppercase(mrid)) at ", m) && occursin("made its slack ($(reason)), voltage setpoint 1.03 pu", m), f2.messages)
        for (slack, others, reason) in (("_f2_sm_b", ("_f2_ei_b", "_f2_sm_b2"), "referencePriority 7"),
                                        ("_f2_ei_c1", ("_f2_ei_c2",), "largest infeed, P = 10.0 MW, only degree-1 buses carry infeed"),
                                        ("_f2_eni_d2", ("_f2_eni_d1", "_f2_ei_d3", "_f2_eni_d2b"), "largest active power reserve off the degree-1 buses, 30.0 MW"),
                                        ("_f2_ei_e1", ("_f2_ei_e2",), "largest infeed, P = 10.0 MW, only degree-1 buses carry infeed"))
          local ps = f2_ps(slack)
          @test Sparlectra.isSlack(ps) && ps.isRegulated && ps.vm_pu ≈ 1.03
          @test Sparlectra.getNodeType(f2.net.nodeVec[Sparlectra.getPosumerBusIndex(ps)]) == Sparlectra.Slack
          @test !any(m -> Sparlectra.isSlack(f2_ps(m)), others)
          @test f2_warn(slack, reason) == 1
        end
        @test count(m -> startswith(m, "warning: electrical island of "), f2.messages) == 4
        @test "notice: 4 electrical island(s) without reference got a slack from the importer (one warning each)" in f2.messages
        @test f2.slack_bus == res.slack_bus
        # the islands solve with their own slack: on C the slack covers the
        # 70.3 MW the other injection draws plus the line loss, and keeps
        # its setpoint
        local f2_net = deepcopy(f2.net)
        @test last(redirect_stdout(() -> runpf!(f2_net, 50, 1e-8, 0; islands_enabled = true), devnull)) == 0
        local c1 = f2_net.nodeVec[f2_bus("_f2_c1")]
        @test c1._vm_pu ≈ 1.03
        @test 70.3 < c1._pƩGen < 71.0
        @test f2_net.prosumpsVec[last(f2.element_by_mrid["_f2_ei_c2"])].pVal == -70.3

        # multi_slack = false keeps the single-reference import: no island
        # slack, no warning
        local f2_single = importCGMES(path = _cgmes_f2_island_store(loaded), name = "sg_f2_single", multi_slack = false)
        @test !any(m -> Sparlectra.isSlack(f2_single.net.prosumpsVec[last(f2_single.element_by_mrid[m])]), ("_f2_sm_b", "_f2_ei_c1", "_f2_eni_d2", "_f2_ei_e1"))
        @test !any(m -> startswith(m, "warning:"), f2_single.messages)

        # G: the drawing machine regulates the voltage of the load bus, so
        # G has a source and stays energized, also with machine_control =
        # true, where the regulation is remote. Being no unit with infeed,
        # the machine is no slack candidate; without another reference it
        # gets no controller and is held PV at its own bus, as by default
        for mc in (false, true)
          local g = mc ? (@test_logs (:warn, r"^CGMES import: 4 data defect") importCGMES(path = _cgmes_f2_island_store(loaded), name = "sg_f2_mc", machine_control = true)) : f2
          local g_ps = g.net.prosumpsVec[last(g.element_by_mrid["_f2_sm_g"])]
          local g1 = g.net.busDict[g.topo.bus_name["_f2_g1"]]
          @test !any(tn -> g.net.busDict[g.topo.bus_name[tn]] in g.net.isoNodes, ("_f2_g1", "_f2_g2"))
          @test g_ps.pVal == -3.0 && g_ps.isRegulated && !Sparlectra.isSlack(g_ps)
          @test Sparlectra.getNodeType(g.net.nodeVec[g1]) == Sparlectra.PV
          @test isempty(Sparlectra._machine_controllers(g.net))
          @test count(m -> occursin("_F2_SM_G", m) && occursin("held PV at its own bus", m), g.messages) == 1
          @test mc == ("notice: machine control _F2_SM_G not attached (its island has no slack, no PV bus and no infeed); held PV at its own bus" in g.messages)
        end
      end)() end

      @testset "island joined through a bus without branches" begin (function ()
        # J: jx has no branch of its own and joins j2 and j3 through two
        # closed switches, one node in the power flow. The marking counts
        # branches only and isolates jx, but the island slack rule works on
        # the islands of the power flow, where J is one island: one slack,
        # the larger infeed on j4. Split at jx, each side got a slack: two
        # references in one solver island, and the solve aborted before
        # Newton
        local fj = @test_logs (:warn, r"^CGMES import: 1 data defect") importCGMES(path = _cgmes_f2_island_store(loaded; linked = true), name = "sg_f2_j")
        local fj_bus(tn) = fj.net.busDict[fj.topo.bus_name[tn]]
        local fj_ps(mrid) = fj.net.prosumpsVec[last(fj.element_by_mrid[mrid])]
        @test fj_bus("_f2_jx") in fj.net.isoNodes
        @test haskey(Sparlectra._linkFedIsolatedBuses(fj.net), fj_bus("_f2_jx"))
        # the node of jx keeps its bus for a voltage measurement
        @test Sparlectra.CGMESImporter.cgmesMeasurementKeys(fj).nodes[Sparlectra._measurement_key("_f2_jx")] == fj_bus("_f2_jx")
        @test Sparlectra.isSlack(fj_ps("_f2_ei_j4")) && !Sparlectra.isSlack(fj_ps("_f2_ei_j1"))
        @test count(m -> startswith(m, "warning: electrical island of "), fj.messages) == 1
        local wnet = first(Sparlectra._merged_pf_net(fj.net))
        Sparlectra.refreshBusTypesFromProsumers!(wnet)
        @test maximum(row -> row.n_ref, Sparlectra.detect_ac_islands(wnet).rows) == 1
        local fj_net = deepcopy(fj.net)
        @test last(redirect_stdout(() -> runpf!(fj_net, 50, 1e-8, 0; islands_enabled = true), devnull)) == 0
        @test fj_net.nodeVec[fj_bus("_f2_jx")]._vm_pu == fj_net.nodeVec[fj_bus("_f2_j2")]._vm_pu == fj_net.nodeVec[fj_bus("_f2_j3")]._vm_pu
      end)() end

      # The island slack choice reads DynamicData.unit_type (#506): a unit
      # stated grid-following is no candidate, an island of grid-following
      # units only is switched off, synchronous and grid-forming units stay
      # candidates. The import fills only the synchronous type (DY profile),
      # so the types are set on an import without island slacks and the
      # island step runs on that net.
      @testset "island slack choice reads the unit type" begin (function ()
        local CI = Sparlectra.CGMESImporter
        local base = importCGMES(path = _cgmes_f2_island_store(loaded), name = "sg_f2_types", multi_slack = false)
        local net = deepcopy(base.net)
        local ps(m) = net.prosumpsVec[last(base.element_by_mrid[m])]
        local bus(tn) = net.busDict[base.topo.bus_name[tn]]
        # B: the machine with referencePriority 7 is grid-following, the
        # EquivalentInjection synchronous; C: both units grid-following; D:
        # the unit with the largest reserve grid-forming; E: the unit with
        # infeed grid-following, the one drawing power of unknown type
        ps("_f2_sm_b").dynamics = DynamicData(unit_type = GridFollowingUnit)
        ps("_f2_ei_b").dynamics = DynamicData(inertia_s = 3.0)
        ps("_f2_ei_c1").dynamics = DynamicData(unit_type = GridFollowingUnit)
        ps("_f2_ei_c2").dynamics = DynamicData(unit_type = GridFollowingUnit)
        ps("_f2_eni_d2").dynamics = DynamicData(virtual_inertia_s = 2.0, unit_type = GridFormingUnit)
        ps("_f2_ei_e1").dynamics = DynamicData(unit_type = GridFollowingUnit)
        local ctx = CI._MapCtx()
        merge!(ctx.element_by_mrid, base.element_by_mrid)
        CI._slackForIslandsWithoutReference!(net, base.store, ctx)
        local warns = filter(m -> startswith(m, "warning: electrical island of "), ctx.messages)
        @test length(warns) == 4
        # B: the synchronous unit takes the slack, not the stronger priority
        @test Sparlectra.isSlack(ps("_f2_ei_b")) && !Sparlectra.isSlack(ps("_f2_sm_b"))
        @test count(m -> occursin("_F2_EI_B at ", m) && occursin("(largest infeed, P = 20.0 MW, only degree-1 buses carry infeed)", m), warns) == 1
        # D: the grid-forming unit stays the candidate it was
        @test Sparlectra.isSlack(ps("_f2_eni_d2"))
        @test count(m -> occursin("_F2_ENI_D2 at ", m) && occursin("(largest active power reserve off the degree-1 buses, 30.0 MW)", m), warns) == 1
        # C: switched off, buses isolated at 0 pu, the line out of service,
        # the 10 MW infeed dropped and the 70.3 MW drawn unsupplied
        for tn in ("_f2_c1", "_f2_c2")
          @test bus(tn) in net.isoNodes
          @test net.nodeVec[bus(tn)]._vm_pu == 0.0
        end
        @test all(m -> (ps(m).pVal, ps(m).qVal) == (0.0, 0.0) && !Sparlectra.isSlack(ps(m)), ("_f2_ei_c1", "_f2_ei_c2"))
        @test net.branchVec[last(base.element_by_mrid["_f2_lc"])].status == 0
        @test count(m -> occursin("has infeed from grid-following units only and cannot stand in operation; switched off, 2 injection(s) zeroed, generation P = 10.0 MW dropped, load P = 70.3 MW unsupplied", m), warns) == 1
        # E: no candidate left, no slack, the island stays energized
        @test !Sparlectra.isSlack(ps("_f2_ei_e1")) && !Sparlectra.isSlack(ps("_f2_ei_e2"))
        @test !(bus("_f2_e1") in net.isoNodes)
        @test count(m -> endswith(m, "has infeed but no reference, and every unit with infeed is grid-following; no slack chosen"), warns) == 1
        @test "notice: 2 electrical island(s) without reference got a slack from the importer (one warning each)" in ctx.messages
        @test "notice: 1 electrical island(s) with grid-following units only switched off (one warning each)" in ctx.messages
      end)() end

      @testset "SV angle alignment on the islands of the power flow" begin (function ()
        # J again: the power flow solves j1 .. j4 as one island through the
        # switch-only bus jx, so compareWithSV removes one angle offset for
        # all of J. The SV of the j3/j4 side is written 2 degrees off the
        # solved state. With one offset for J the 2 degrees stay visible on
        # one side; split at jx, each side took its own offset and the
        # deviation was reported as 0.
        local js = _cgmes_f2_island_store(loaded; linked = true)
        local j_tns = ("_f2_j1", "_f2_j2", "_f2_jx", "_f2_j3", "_f2_j4")
        local solved = @test_logs (:warn, r"^CGMES import: 1 data defect") importCGMES(path = js, name = "sg_f2_j_solve")
        @test last(redirect_stdout(() -> runpf!(solved.net, 50, 1e-8, 0; islands_enabled = true), devnull)) == 0
        for tn in j_tns
          local nd = solved.net.nodeVec[solved.net.busDict[solved.topo.bus_name[tn]]]
          local shift = tn in ("_f2_j3", "_f2_j4") ? 2.0 : 0.0
          js.objects[string(tn, "_sv")].attrs[:v] = string(nd._vm_pu * solved.topo.vn_kV[tn])
          js.objects[string(tn, "_sv")].attrs[:angle] = string(nd._va_deg + shift)
        end
        local res_j = @test_logs (:warn, r"^CGMES import: 1 data defect") importCGMES(path = js, name = "sg_f2_j_sv")
        @test last(redirect_stdout(() -> runpf!(res_j.net, 50, 1e-8, 0; islands_enabled = true), devnull)) == 0
        local isl = Sparlectra._power_flow_island_of_bus(res_j.net)
        local jb = [res_j.net.busDict[res_j.topo.bus_name[tn]] for tn in j_tns]
        @test length(unique(isl[b] for b in jb)) == 1
        local rows = Dict(r.bus => r for r in compareWithSV(res_j).rows)
        local dva(tn) = rows[res_j.topo.bus_name[tn]].dva_aligned
        @test all(tn -> abs(dva(tn)) < 1e-6, ("_f2_j1", "_f2_j2"))
        @test all(tn -> isapprox(dva(tn), -2.0; atol = 1e-6), ("_f2_j3", "_f2_j4"))
      end)() end
    end)() end

    # The deliveries tested state no network frequency (no BaseFrequency),
    # so the key decides: the lines arrive in siemens (b already at the
    # network frequency), the power flow is bit identical at 50 and 60 Hz,
    # and the conversion to SCF writes the key value as meta.f_nom.
    @testset "SmallGrid: system frequency from the key" begin (function ()
      dir = joinpath(@__DIR__, "fixtures", "cgmes_smallgrid")
      cfg_at(f) = Sparlectra.SparlectraConfig(model = Sparlectra.ModelConfig(nominal_frequency_hz = f))
      imp(f) = redirect_stdout(() -> Sparlectra.import_case(dir, cfg_at(f)), devnull)
      n50 = imp(50.0).net
      n60 = imp(60.0).net
      @test (n50.nominal_frequency_hz, n60.nominal_frequency_hz) == (50.0, 60.0)
      @test Sparlectra.createYBUS(net = n60, sparse = true) == Sparlectra.createYBUS(net = n50, sparse = true)
      @test runpf!(n60, 40, 1e-8, 0) == runpf!(n50, 40, 1e-8, 0)
      @test [nd._vm_pu for nd in n60.nodeVec] == [nd._vm_pu for nd in n50.nodeVec]
      @test [nd._va_deg for nd in n60.nodeVec] == [nd._va_deg for nd in n50.nodeVec]
      @test Sparlectra.cgmes_adapter_options(cfg_at(60.0)).f_nom == 60.0
      case60 = redirect_stdout(() -> convert_case(CGMESAdapter(), dir, Sparlectra.CGMESAdapterOptions(f_nom = 60.0)), devnull)
      @test case60.sparlectra.meta["f_nom"] == 60.0
    end)() end
  end)() end
end
