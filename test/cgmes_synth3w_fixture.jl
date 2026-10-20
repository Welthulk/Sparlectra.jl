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

# file: test/cgmes_synth3w_fixture.jl
# purpose: the synthetic three-winding delivery shared by the importer and the
#          series union tests (written to a temporary directory on demand)

# Synthetic 3-bus star delivery (EQ+SSH+TP+SV): slack ENI at HV, loads at
# MV/LV, one 3W transformer whose end-2 RatioTapChanger carries an enabled
# TapChangerControl regulating the MV bus. Purpose-built for #294 point 4,
# controllers on star-equivalent 3W legs, because no cached ENTSO-E set
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

# An unregulated RatioTapChanger on end 1 (1 % per step, neutral 11), added
# by `_cgmes_synth3w_dir(end1_tap_step = k)` at SSH step k. Off neutral it
# moves the star voltage away from the HV voltage.
const _CGMES_SYNTH3W_RTC1_EQ = """<cim:RatioTapChanger rdf:ID="_rtc1">
<cim:IdentifiedObject.name>T3W_HV_TAP</cim:IdentifiedObject.name>
<cim:RatioTapChanger.TransformerEnd rdf:resource="#_t3w_e1"/>
<cim:RatioTapChanger.stepVoltageIncrement>1.0</cim:RatioTapChanger.stepVoltageIncrement>
<cim:TapChanger.lowStep>1</cim:TapChanger.lowStep>
<cim:TapChanger.highStep>21</cim:TapChanger.highStep>
<cim:TapChanger.neutralStep>11</cim:TapChanger.neutralStep>
</cim:RatioTapChanger>
"""

# Magnetizing admittance g, b (siemens on the end's ratedU) per end, written
# into the EQ by `_cgmes_synth3w_dir(magnetizing = true)`.
const _CGMES_SYNTH3W_MAG = [(1.0e-6, -4.0e-6), (4.0e-6, -1.6e-5), (1.0e-4, -5.0e-4)]

# `open_end = k` writes the terminal of end k with connected = false in the
# SSH, the transformer stays in service.
function _cgmes_synth3w_dir(; end1_tap_step::Union{Nothing,Int} = nothing, magnetizing::Bool = false, open_end::Union{Nothing,Int} = nothing)::String
  dir = mktempdir()
  eq, ssh = _CGMES_SYNTH3W_EQ, _CGMES_SYNTH3W_SSH
  if open_end !== nothing
    tline = "<cim:Terminal rdf:about=\"#_t_t3w_$(open_end)\"><cim:ACDCTerminal.connected>true<"
    occursin(tline, ssh) || error("no SSH terminal for end $(open_end)")
    ssh = replace(ssh, tline => replace(tline, "true" => "false"))
  end
  if magnetizing
    for (k, (g, b)) in enumerate(_CGMES_SYNTH3W_MAG)
      tline = "<cim:TransformerEnd.Terminal rdf:resource=\"#_t_t3w_$(k)\"/>"
      eq = replace(eq, tline => "<cim:PowerTransformerEnd.g>$(g)</cim:PowerTransformerEnd.g>\n<cim:PowerTransformerEnd.b>$(b)</cim:PowerTransformerEnd.b>\n" * tline)
    end
  end
  if end1_tap_step !== nothing
    eq = replace(eq, "</rdf:RDF>" => _CGMES_SYNTH3W_RTC1_EQ * "</rdf:RDF>")
    ssh = replace(ssh, "</rdf:RDF>" => """<cim:RatioTapChanger rdf:about="#_rtc1"><cim:TapChanger.step>$(end1_tap_step)</cim:TapChanger.step><cim:TapChanger.controlEnabled>false</cim:TapChanger.controlEnabled></cim:RatioTapChanger>\n</rdf:RDF>""")
  end
  write(joinpath(dir, "synth3w_EQ.xml"), eq)
  write(joinpath(dir, "synth3w_SSH.xml"), ssh)
  write(joinpath(dir, "synth3w_TP.xml"), _CGMES_SYNTH3W_TP)
  write(joinpath(dir, "synth3w_SV.xml"), _CGMES_SYNTH3W_SV)
  return dir
end

# An SV profile for the synthetic delivery from given values: one SvVoltage
# per TopologicalNode (`"_tn_hv" => (kV, degrees)`) and one SvPowerFlow per
# terminal (`"_t_t3w_1" => (MW, Mvar)`, flow into the equipment), written at
# full precision so a re-import reads back the same numbers.
function _cgmes_synth3w_sv(voltages, flows)::String
  io = IOBuffer()
  println(io, """<?xml version="1.0" encoding="UTF-8"?>""")
  println(io, """<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:cim="http://iec.ch/TC57/2013/CIM-schema-cim16#" xmlns:md="http://iec.ch/TC57/61970-552/ModelDescription/1#">""")
  println(io, """<md:FullModel rdf:about="urn:uuid:sv3w"><md:Model.profile>http://entsoe.eu/CIM/StateVariables/4/1</md:Model.profile></md:FullModel>""")
  for (tn, (v, a)) in voltages
    println(io, """<cim:SvVoltage rdf:ID="_sv$(tn)"><cim:SvVoltage.TopologicalNode rdf:resource="#$(tn)"/><cim:SvVoltage.v>$(repr(v))</cim:SvVoltage.v><cim:SvVoltage.angle>$(repr(a))</cim:SvVoltage.angle></cim:SvVoltage>""")
  end
  for (t, (p, q)) in flows
    println(io, """<cim:SvPowerFlow rdf:ID="_svpf$(t)"><cim:SvPowerFlow.Terminal rdf:resource="#$(t)"/><cim:SvPowerFlow.p>$(repr(p))</cim:SvPowerFlow.p><cim:SvPowerFlow.q>$(repr(q))</cim:SvPowerFlow.q></cim:SvPowerFlow>""")
  end
  println(io, "</rdf:RDF>")
  return String(take!(io))
end

# Synthetic deliveries with one tap changer off its neutral step: a 3W
# transformer tapped at end 1, or a 2W transformer tapped at end 2 (end 1
# without changer). The SV state is computed here in physical units with the
# arrangement CGMES uses for a tapped end: the ideal transformer of the tap
# sits at the terminal of its own end, the series impedance (referred
# through the rated voltages only) lies on the inner side. A 3W leg has the
# magnetizing admittance of its end at the star node (end 1 tapped, end 3
# untapped carry one); a 2W transformer has the total of both ends, half on
# each side of the impedance (the SV conventions the importer follows since
# 0.41.0 for two-winding and 0.42.0 for three-winding transformers).
# Loads take the terminal flows, so the SV is a solved state of that
# arrangement. `kind` is :three_winding or :two_winding.
# `table` (2W only) gives the changer a RatioTapChangerTable whose point at
# the step carries r, x, g, b deviations in percent, `(r = .., x = .., g =
# .., b = ..)`; the point ratio is the linear correction, and the SV applies
# the deviations to the series impedance and the magnetizing admittance.
# `control_kv` (with `table`) enables a voltage TapChangerControl at BUS2
# with that target and writes a table point for every step, its deviations
# those of `table` scaled by step / start step (zero at the neutral step).
function _cgmes_synth_tapside_dir(kind::Symbol; table::Union{Nothing,NamedTuple} = nothing, control_kv::Union{Nothing,Float64} = nothing)::String
  kind in (:three_winding, :two_winding) || throw(ArgumentError("kind $(kind)"))
  three = kind === :three_winding
  (three && table !== nothing) && throw(ArgumentError("a tap table is built for the two-winding delivery only"))
  (control_kv !== nothing && table === nothing) && throw(ArgumentError("control_kv needs a table"))
  fmt(x) = string(round(x; digits = 12))
  pol(v, deg) = v * cis(deg2rad(deg))
  if three
    vn = [220.0, 110.0, 20.0]
    rated = [225.0, 115.0, 21.0]
    z = [complex(0.4, 24.0), complex(0.1, 6.0), complex(0.02, 0.8)]   # ohm on the end's ratedU
    y = [complex(1.0e-6, -4.0e-6), 0.0im, complex(1.0e-4, -5.0e-4)]    # siemens on the end's ratedU
    tap_end, neutral, step, incr = 1, 11, 15, 1.25
    V = [pol(220.0, 0.0), pol(112.3, -1.7), pol(20.4, -2.3)]
  else
    vn = [110.0, 20.0]
    rated = [115.0, 21.0]
    z = [complex(0.2, 8.0), complex(0.01, 0.4)]
    y = [complex(2.0e-6, -1.0e-5), complex(1.0e-4, -5.0e-4)]
    tap_end, neutral, step, incr = 2, 0, -4, 1.5
    V = [pol(113.0, 0.0), pol(20.6, -3.1)]
  end
  corr = 1.0 + (step - neutral) * incr / 100.0
  n = length(vn)
  S = zeros(ComplexF64, n)
  if three
    # leg k: node k, ideal ratio (rated_k * c_k / rated_1), inner point m_k,
    # impedance to the star; the magnetizing admittance of every end at the
    # star node, referred through the rated voltages
    c = [k == tap_end ? corr : 1.0 for k in 1:n]
    Vm = [V[k] / (rated[k] * c[k] / rated[1]) for k in 1:n]
    zs = [z[k] * (rated[1] / rated[k])^2 for k in 1:n]
    ys = [y[k] * (rated[k] / rated[1])^2 for k in 1:n]
    Vs = sum(Vm[k] / zs[k] for k in 1:n) / (sum(1.0 / zs[k] for k in 1:n) + sum(ys))
    for k in 1:n
      S[k] = Vm[k] * conj((Vm[k] - Vs) / zs[k])
    end
  else
    # node 2, ideal ratio (rated_2 * c / rated_1), inner point m, impedance
    # (both ends on the end-1 rated voltage) to node 1; the magnetizing
    # admittance of both ends on the end-1 rated voltage, half at m and
    # half at node 1
    Vm = V[2] / (rated[2] * corr / rated[1])
    zt = z[1] + z[2] * (rated[1] / rated[2])^2
    yh = (y[1] + y[2] * (rated[2] / rated[1])^2) / 2
    if table !== nothing
      zt = complex(real(zt) * (1 + table.r / 100), imag(zt) * (1 + table.x / 100))
      yh = complex(real(yh) * (1 + table.g / 100), imag(yh) * (1 + table.b / 100))
    end
    S[2] = Vm * conj(yh * Vm + (Vm - V[1]) / zt)
    S[1] = V[1] * conj(yh * V[1] + (V[1] - Vm) / zt)
  end
  tag = three ? "tap3w" : "tap2w"
  hdr(p, prof) = """<?xml version="1.0" encoding="UTF-8"?>
<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:cim="http://iec.ch/TC57/2013/CIM-schema-cim16#" xmlns:md="http://iec.ch/TC57/61970-552/ModelDescription/1#">
<md:FullModel rdf:about="urn:uuid:$(tag)_$(p)"><md:Model.profile>$(prof)</md:Model.profile></md:FullModel>
"""
  eq = IOBuffer(); ssh = IOBuffer(); tp = IOBuffer(); sv = IOBuffer()
  print(eq, hdr("eq", "http://entsoe.eu/CIM/EquipmentCore/3/1"))
  print(ssh, hdr("ssh", "http://entsoe.eu/CIM/SteadyStateHypothesis/1/1"))
  print(tp, hdr("tp", "http://entsoe.eu/CIM/Topology/4/1"))
  print(sv, hdr("sv", "http://entsoe.eu/CIM/StateVariables/4/1"))
  terminal(io, id, eqid, seq) = println(io, """<cim:Terminal rdf:ID="$(id)"><cim:Terminal.ConductingEquipment rdf:resource="#$(eqid)"/><cim:ACDCTerminal.sequenceNumber>$(seq)</cim:ACDCTerminal.sequenceNumber></cim:Terminal>""")
  wire(id, k) = (println(ssh, """<cim:Terminal rdf:about="#$(id)"><cim:ACDCTerminal.connected>true</cim:ACDCTerminal.connected></cim:Terminal>"""); println(tp, """<cim:Terminal rdf:about="#$(id)"><cim:Terminal.TopologicalNode rdf:resource="#_tn$(k)"/></cim:Terminal>"""))
  flow(id, s) = println(sv, """<cim:SvPowerFlow rdf:ID="$(id)_pf"><cim:SvPowerFlow.Terminal rdf:resource="#$(id)"/><cim:SvPowerFlow.p>$(fmt(real(s)))</cim:SvPowerFlow.p><cim:SvPowerFlow.q>$(fmt(imag(s)))</cim:SvPowerFlow.q></cim:SvPowerFlow>""")
  for k in 1:n
    println(eq, """<cim:BaseVoltage rdf:ID="_bv$(k)"><cim:BaseVoltage.nominalVoltage>$(vn[k])</cim:BaseVoltage.nominalVoltage></cim:BaseVoltage>""")
    println(tp, """<cim:TopologicalNode rdf:ID="_tn$(k)"><cim:IdentifiedObject.name>BUS$(k)</cim:IdentifiedObject.name><cim:TopologicalNode.BaseVoltage rdf:resource="#_bv$(k)"/></cim:TopologicalNode>""")
    println(sv, """<cim:SvVoltage rdf:ID="_sv$(k)"><cim:SvVoltage.TopologicalNode rdf:resource="#_tn$(k)"/><cim:SvVoltage.v>$(fmt(abs(V[k])))</cim:SvVoltage.v><cim:SvVoltage.angle>$(fmt(rad2deg(angle(V[k]))))</cim:SvVoltage.angle></cim:SvVoltage>""")
  end
  # slack infeed at bus 1, a load at every other bus taking the terminal flow
  println(eq, """<cim:ExternalNetworkInjection rdf:ID="_eni"><cim:IdentifiedObject.name>GRID</cim:IdentifiedObject.name></cim:ExternalNetworkInjection>""")
  terminal(eq, "_t_eni", "_eni", 1); wire("_t_eni", 1); flow("_t_eni", -S[1])
  println(ssh, """<cim:ExternalNetworkInjection rdf:about="#_eni"><cim:ExternalNetworkInjection.p>$(fmt(-real(S[1])))</cim:ExternalNetworkInjection.p><cim:ExternalNetworkInjection.q>$(fmt(-imag(S[1])))</cim:ExternalNetworkInjection.q></cim:ExternalNetworkInjection>""")
  for k in 2:n
    println(eq, """<cim:EnergyConsumer rdf:ID="_load$(k)"><cim:IdentifiedObject.name>LOAD$(k)</cim:IdentifiedObject.name></cim:EnergyConsumer>""")
    terminal(eq, "_t_load$(k)", "_load$(k)", 1); wire("_t_load$(k)", k); flow("_t_load$(k)", -S[k])
    println(ssh, """<cim:EnergyConsumer rdf:about="#_load$(k)"><cim:EnergyConsumer.p>$(fmt(-real(S[k])))</cim:EnergyConsumer.p><cim:EnergyConsumer.q>$(fmt(-imag(S[k])))</cim:EnergyConsumer.q></cim:EnergyConsumer>""")
  end
  println(eq, """<cim:PowerTransformer rdf:ID="_pt"><cim:IdentifiedObject.name>$(three ? "T3W" : "T2W")</cim:IdentifiedObject.name></cim:PowerTransformer>""")
  for k in 1:n
    println(eq, """<cim:PowerTransformerEnd rdf:ID="_pt_e$(k)"><cim:PowerTransformerEnd.PowerTransformer rdf:resource="#_pt"/><cim:TransformerEnd.endNumber>$(k)</cim:TransformerEnd.endNumber><cim:PowerTransformerEnd.ratedU>$(rated[k])</cim:PowerTransformerEnd.ratedU><cim:PowerTransformerEnd.ratedS>100</cim:PowerTransformerEnd.ratedS><cim:PowerTransformerEnd.r>$(real(z[k]))</cim:PowerTransformerEnd.r><cim:PowerTransformerEnd.x>$(imag(z[k]))</cim:PowerTransformerEnd.x><cim:PowerTransformerEnd.g>$(real(y[k]))</cim:PowerTransformerEnd.g><cim:PowerTransformerEnd.b>$(imag(y[k]))</cim:PowerTransformerEnd.b><cim:TransformerEnd.Terminal rdf:resource="#_t_pt$(k)"/></cim:PowerTransformerEnd>""")
    terminal(eq, "_t_pt$(k)", "_pt", k); wire("_t_pt$(k)", k); flow("_t_pt$(k)", S[k])
  end
  tcc = control_kv === nothing ? "" : "<cim:TapChanger.TapChangerControl rdf:resource=\"#_tcc\"/>"
  println(eq, """<cim:RatioTapChanger rdf:ID="_rtc"><cim:IdentifiedObject.name>RTC</cim:IdentifiedObject.name>$(tcc)<cim:RatioTapChanger.TransformerEnd rdf:resource="#_pt_e$(tap_end)"/><cim:RatioTapChanger.stepVoltageIncrement>$(incr)</cim:RatioTapChanger.stepVoltageIncrement><cim:TapChanger.lowStep>$(neutral - 10)</cim:TapChanger.lowStep><cim:TapChanger.highStep>$(neutral + 10)</cim:TapChanger.highStep><cim:TapChanger.neutralStep>$(neutral)</cim:TapChanger.neutralStep>$(table === nothing ? "" : "<cim:RatioTapChanger.RatioTapChangerTable rdf:resource=\"#_rtct\"/>")</cim:RatioTapChanger>""")
  if table !== nothing
    println(eq, """<cim:RatioTapChangerTable rdf:ID="_rtct"><cim:IdentifiedObject.name>RTC_TABLE</cim:IdentifiedObject.name></cim:RatioTapChangerTable>""")
    for s in (control_kv === nothing ? (step:step) : ((neutral - 10):(neutral + 10)))
      w = (s - neutral) / (step - neutral)
      println(eq, """<cim:RatioTapChangerTablePoint rdf:ID="_rtct_p$(s - neutral + 10)"><cim:RatioTapChangerTablePoint.RatioTapChangerTable rdf:resource="#_rtct"/><cim:TapChangerTablePoint.step>$(s)</cim:TapChangerTablePoint.step><cim:TapChangerTablePoint.ratio>$(1.0 + (s - neutral) * incr / 100.0)</cim:TapChangerTablePoint.ratio><cim:TapChangerTablePoint.r>$(table.r * w)</cim:TapChangerTablePoint.r><cim:TapChangerTablePoint.x>$(table.x * w)</cim:TapChangerTablePoint.x><cim:TapChangerTablePoint.g>$(table.g * w)</cim:TapChangerTablePoint.g><cim:TapChangerTablePoint.b>$(table.b * w)</cim:TapChangerTablePoint.b></cim:RatioTapChangerTablePoint>""")
    end
  end
  if control_kv !== nothing
    println(eq, """<cim:TapChangerControl rdf:ID="_tcc"><cim:IdentifiedObject.name>TCC</cim:IdentifiedObject.name><cim:RegulatingControl.mode rdf:resource="http://iec.ch/TC57/2013/CIM-schema-cim16#RegulatingControlModeKind.voltage"/><cim:RegulatingControl.Terminal rdf:resource="#_t_load2"/></cim:TapChangerControl>""")
    println(ssh, """<cim:TapChangerControl rdf:about="#_tcc"><cim:RegulatingControl.enabled>true</cim:RegulatingControl.enabled><cim:RegulatingControl.targetValue>$(control_kv)</cim:RegulatingControl.targetValue><cim:RegulatingControl.targetDeadband>0.05</cim:RegulatingControl.targetDeadband><cim:RegulatingControl.discrete>true</cim:RegulatingControl.discrete></cim:TapChangerControl>""")
  end
  println(ssh, """<cim:RatioTapChanger rdf:about="#_rtc"><cim:TapChanger.step>$(step)</cim:TapChanger.step><cim:TapChanger.controlEnabled>$(control_kv !== nothing)</cim:TapChanger.controlEnabled></cim:RatioTapChanger>""")
  println(sv, """<cim:SvTapStep rdf:ID="_svtap"><cim:SvTapStep.TapChanger rdf:resource="#_rtc"/><cim:SvTapStep.position>$(step)</cim:SvTapStep.position></cim:SvTapStep>""")
  for io in (eq, ssh, tp, sv)
    println(io, "</rdf:RDF>")
  end
  dir = mktempdir()
  for (p, io) in (("EQ", eq), ("SSH", ssh), ("TP", tp), ("SV", sv))
    write(joinpath(dir, string(tag, "_", p, ".xml")), take!(io))
  end
  return dir
end

# A heavily loaded variant of the same delivery with two power-flow
# solutions: the tap changer sits on end 1 at step 1 (ten steps off
# neutral), the MV load is 425 MW. The operating solution has the MV bus at
# 0.954 pu, the low-voltage solution at 0.456 pu (solved on this model,
# rounded to the written digits). `sv = :high` writes the operating solution
# as the SV profile, `sv = :low` the low-voltage one; the SSH injections are
# the same in both. Started from its SV, Newton stays on that solution. Used
# by the plausibility check of a converged run.
function _cgmes_synth3w_heavy_dir(; sv::Symbol = :high)::String
  sv in (:high, :low) || throw(ArgumentError("sv must be :high or :low"))
  eq = replace(_CGMES_SYNTH3W_EQ, "<cim:RatioTapChanger.TransformerEnd rdf:resource=\"#_t3w_e2\"/>" => "<cim:RatioTapChanger.TransformerEnd rdf:resource=\"#_t3w_e1\"/>")
  ssh = replace(_CGMES_SYNTH3W_SSH,
    "<cim:ExternalNetworkInjection.p>-45</cim:ExternalNetworkInjection.p><cim:ExternalNetworkInjection.q>-12</cim:ExternalNetworkInjection.q>" => "<cim:ExternalNetworkInjection.p>-460.2</cim:ExternalNetworkInjection.p><cim:ExternalNetworkInjection.q>-401.5</cim:ExternalNetworkInjection.q>",
    "<cim:EnergyConsumer.p>40</cim:EnergyConsumer.p><cim:EnergyConsumer.q>10</cim:EnergyConsumer.q>" => "<cim:EnergyConsumer.p>425</cim:EnergyConsumer.p><cim:EnergyConsumer.q>85</cim:EnergyConsumer.q>",
    "<cim:EnergyConsumer.p>5</cim:EnergyConsumer.p><cim:EnergyConsumer.q>1</cim:EnergyConsumer.q>" => "<cim:EnergyConsumer.p>30</cim:EnergyConsumer.p><cim:EnergyConsumer.q>6</cim:EnergyConsumer.q>",
    "<cim:TapChanger.step>11</cim:TapChanger.step>" => "<cim:TapChanger.step>1</cim:TapChanger.step>")
  mv, lv = sv === :high ? (("104.896", "-23.492"), ("20.187", "-14.385")) : (("50.158", "-55.481"), ("13.946", "-22.634"))
  svx = replace(_CGMES_SYNTH3W_SV,
    "<cim:SvVoltage.v>111.5</cim:SvVoltage.v><cim:SvVoltage.angle>-1.5</cim:SvVoltage.angle>" => "<cim:SvVoltage.v>$(mv[1])</cim:SvVoltage.v><cim:SvVoltage.angle>$(mv[2])</cim:SvVoltage.angle>",
    "<cim:SvVoltage.v>19.9</cim:SvVoltage.v><cim:SvVoltage.angle>-2.0</cim:SvVoltage.angle>" => "<cim:SvVoltage.v>$(lv[1])</cim:SvVoltage.v><cim:SvVoltage.angle>$(lv[2])</cim:SvVoltage.angle>")
  dir = mktempdir()
  write(joinpath(dir, "synth3w_EQ.xml"), eq)
  write(joinpath(dir, "synth3w_SSH.xml"), ssh)
  write(joinpath(dir, "synth3w_TP.xml"), _CGMES_SYNTH3W_TP)
  write(joinpath(dir, "synth3w_SV.xml"), svx)
  return dir
end
