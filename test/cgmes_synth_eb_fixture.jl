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

# file: test/cgmes_synth_eb_fixture.jl
# purpose: the synthetic EquivalentBranch delivery shared by the importer
#          and the export tests (written to a temporary directory on demand)

# Four buses (EQ+SSH+TP+SV): slack ENI at BUS_A, a line A-B, loads at BUS_B
# and at the 220 kV BUS_D. BUS_C and BUS_D hang on equivalent branches only:
# - EB_SAME, B-C on 110 kV, with r21 equal to r (symmetric),
# - EB_ASYM, B-C on 110 kV parallel to EB_SAME, z21 = mag * z12 * e^{-j12deg}
#   (asymmetric, on the same bus pair as a line EquivalentBranch),
# - EB_CROSS, C (terminal 1, 110 kV) to D (terminal 2, 220 kV),
#   z21 = mag * z12 * e^{j8deg} (asymmetric),
# - EB_LOOP, both terminals on BUS_A,
# - EB_OPEN, A-C with terminal 2 disconnected.
# `mag` is |z21|/|z12| of the two asymmetric branches: 1 is a pure rotation,
# which the phase stamp of the importer reproduces exactly.
# The SV profile is the solved target model (see _cgmes_synth_eb_state),
# with the exact two-port I1 = (V1 - V2)/z12, I2 = (V2 - V1)/z21 for every
# branch: every SvVoltage and every SvPowerFlow of a connected terminal
# comes from that solution.
_cgmes_synth_eb_z21(mag::Float64)::NTuple{2,ComplexF64} = (mag * complex(1.5, 7.0) * cis(deg2rad(-12.0)), mag * complex(0.6, 4.84) * cis(deg2rad(8.0)))

_cgmes_synth_eb_eq(mag::Float64)::String = (z21 = _cgmes_synth_eb_z21(mag); """<?xml version="1.0" encoding="UTF-8"?>
<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:cim="http://iec.ch/TC57/2013/CIM-schema-cim16#" xmlns:md="http://iec.ch/TC57/61970-552/ModelDescription/1#">
<md:FullModel rdf:about="urn:uuid:eqeb"><md:Model.profile>http://entsoe.eu/CIM/EquipmentCore/3/1</md:Model.profile></md:FullModel>
<cim:BaseVoltage rdf:ID="_bv110"><cim:BaseVoltage.nominalVoltage>110</cim:BaseVoltage.nominalVoltage></cim:BaseVoltage>
<cim:BaseVoltage rdf:ID="_bv220"><cim:BaseVoltage.nominalVoltage>220</cim:BaseVoltage.nominalVoltage></cim:BaseVoltage>
<cim:ExternalNetworkInjection rdf:ID="_eni">
<cim:IdentifiedObject.name>GRID</cim:IdentifiedObject.name>
<cim:ExternalNetworkInjection.referencePriority>1</cim:ExternalNetworkInjection.referencePriority>
</cim:ExternalNetworkInjection>
<cim:Terminal rdf:ID="_t_eni"><cim:Terminal.ConductingEquipment rdf:resource="#_eni"/><cim:ACDCTerminal.sequenceNumber>1</cim:ACDCTerminal.sequenceNumber></cim:Terminal>
<cim:EnergyConsumer rdf:ID="_load_b"><cim:IdentifiedObject.name>LOAD_B</cim:IdentifiedObject.name></cim:EnergyConsumer>
<cim:Terminal rdf:ID="_t_load_b"><cim:Terminal.ConductingEquipment rdf:resource="#_load_b"/><cim:ACDCTerminal.sequenceNumber>1</cim:ACDCTerminal.sequenceNumber></cim:Terminal>
<cim:EnergyConsumer rdf:ID="_load_d"><cim:IdentifiedObject.name>LOAD_D</cim:IdentifiedObject.name></cim:EnergyConsumer>
<cim:Terminal rdf:ID="_t_load_d"><cim:Terminal.ConductingEquipment rdf:resource="#_load_d"/><cim:ACDCTerminal.sequenceNumber>1</cim:ACDCTerminal.sequenceNumber></cim:Terminal>
<cim:ACLineSegment rdf:ID="_l_ab">
<cim:IdentifiedObject.name>L_AB</cim:IdentifiedObject.name>
<cim:ACLineSegment.r>2.42</cim:ACLineSegment.r>
<cim:ACLineSegment.x>12.1</cim:ACLineSegment.x>
<cim:ACLineSegment.bch>0</cim:ACLineSegment.bch>
<cim:ConductingEquipment.BaseVoltage rdf:resource="#_bv110"/>
</cim:ACLineSegment>
<cim:Terminal rdf:ID="_t_lab_1"><cim:Terminal.ConductingEquipment rdf:resource="#_l_ab"/><cim:ACDCTerminal.sequenceNumber>1</cim:ACDCTerminal.sequenceNumber></cim:Terminal>
<cim:Terminal rdf:ID="_t_lab_2"><cim:Terminal.ConductingEquipment rdf:resource="#_l_ab"/><cim:ACDCTerminal.sequenceNumber>2</cim:ACDCTerminal.sequenceNumber></cim:Terminal>
<cim:EquivalentNetwork rdf:ID="_eqn"><cim:IdentifiedObject.name>EQN</cim:IdentifiedObject.name></cim:EquivalentNetwork>
<cim:EquivalentBranch rdf:ID="_eb_same">
<cim:IdentifiedObject.name>EB_SAME</cim:IdentifiedObject.name>
<cim:EquivalentEquipment.EquivalentNetwork rdf:resource="#_eqn"/>
<cim:EquivalentBranch.r>1.21</cim:EquivalentBranch.r>
<cim:EquivalentBranch.x>6.05</cim:EquivalentBranch.x>
<cim:EquivalentBranch.r21>1.21</cim:EquivalentBranch.r21>
</cim:EquivalentBranch>
<cim:Terminal rdf:ID="_t_ebs_1"><cim:Terminal.ConductingEquipment rdf:resource="#_eb_same"/><cim:ACDCTerminal.sequenceNumber>1</cim:ACDCTerminal.sequenceNumber></cim:Terminal>
<cim:Terminal rdf:ID="_t_ebs_2"><cim:Terminal.ConductingEquipment rdf:resource="#_eb_same"/><cim:ACDCTerminal.sequenceNumber>2</cim:ACDCTerminal.sequenceNumber></cim:Terminal>
<cim:EquivalentBranch rdf:ID="_eb_asym">
<cim:IdentifiedObject.name>EB_ASYM</cim:IdentifiedObject.name>
<cim:EquivalentEquipment.EquivalentNetwork rdf:resource="#_eqn"/>
<cim:EquivalentBranch.r>1.5</cim:EquivalentBranch.r>
<cim:EquivalentBranch.x>7.0</cim:EquivalentBranch.x>
<cim:EquivalentBranch.r21>$(real(z21[1]))</cim:EquivalentBranch.r21>
<cim:EquivalentBranch.x21>$(imag(z21[1]))</cim:EquivalentBranch.x21>
</cim:EquivalentBranch>
<cim:Terminal rdf:ID="_t_eba_1"><cim:Terminal.ConductingEquipment rdf:resource="#_eb_asym"/><cim:ACDCTerminal.sequenceNumber>1</cim:ACDCTerminal.sequenceNumber></cim:Terminal>
<cim:Terminal rdf:ID="_t_eba_2"><cim:Terminal.ConductingEquipment rdf:resource="#_eb_asym"/><cim:ACDCTerminal.sequenceNumber>2</cim:ACDCTerminal.sequenceNumber></cim:Terminal>
<cim:EquivalentBranch rdf:ID="_eb_cross">
<cim:IdentifiedObject.name>EB_CROSS</cim:IdentifiedObject.name>
<cim:EquivalentEquipment.EquivalentNetwork rdf:resource="#_eqn"/>
<cim:EquivalentBranch.r>0.6</cim:EquivalentBranch.r>
<cim:EquivalentBranch.x>4.84</cim:EquivalentBranch.x>
<cim:EquivalentBranch.r21>$(real(z21[2]))</cim:EquivalentBranch.r21>
<cim:EquivalentBranch.x21>$(imag(z21[2]))</cim:EquivalentBranch.x21>
</cim:EquivalentBranch>
<cim:Terminal rdf:ID="_t_ebc_1"><cim:Terminal.ConductingEquipment rdf:resource="#_eb_cross"/><cim:ACDCTerminal.sequenceNumber>1</cim:ACDCTerminal.sequenceNumber></cim:Terminal>
<cim:Terminal rdf:ID="_t_ebc_2"><cim:Terminal.ConductingEquipment rdf:resource="#_eb_cross"/><cim:ACDCTerminal.sequenceNumber>2</cim:ACDCTerminal.sequenceNumber></cim:Terminal>
<cim:EquivalentBranch rdf:ID="_eb_loop">
<cim:IdentifiedObject.name>EB_LOOP</cim:IdentifiedObject.name>
<cim:EquivalentEquipment.EquivalentNetwork rdf:resource="#_eqn"/>
<cim:EquivalentBranch.r>0.1</cim:EquivalentBranch.r>
<cim:EquivalentBranch.x>0.5</cim:EquivalentBranch.x>
</cim:EquivalentBranch>
<cim:Terminal rdf:ID="_t_ebl_1"><cim:Terminal.ConductingEquipment rdf:resource="#_eb_loop"/><cim:ACDCTerminal.sequenceNumber>1</cim:ACDCTerminal.sequenceNumber></cim:Terminal>
<cim:Terminal rdf:ID="_t_ebl_2"><cim:Terminal.ConductingEquipment rdf:resource="#_eb_loop"/><cim:ACDCTerminal.sequenceNumber>2</cim:ACDCTerminal.sequenceNumber></cim:Terminal>
<cim:EquivalentBranch rdf:ID="_eb_open">
<cim:IdentifiedObject.name>EB_OPEN</cim:IdentifiedObject.name>
<cim:EquivalentEquipment.EquivalentNetwork rdf:resource="#_eqn"/>
<cim:EquivalentBranch.r>2.0</cim:EquivalentBranch.r>
<cim:EquivalentBranch.x>9.0</cim:EquivalentBranch.x>
</cim:EquivalentBranch>
<cim:Terminal rdf:ID="_t_ebo_1"><cim:Terminal.ConductingEquipment rdf:resource="#_eb_open"/><cim:ACDCTerminal.sequenceNumber>1</cim:ACDCTerminal.sequenceNumber></cim:Terminal>
<cim:Terminal rdf:ID="_t_ebo_2"><cim:Terminal.ConductingEquipment rdf:resource="#_eb_open"/><cim:ACDCTerminal.sequenceNumber>2</cim:ACDCTerminal.sequenceNumber></cim:Terminal>
</rdf:RDF>
""")

const _CGMES_SYNTH_EB_SSH = """<?xml version="1.0" encoding="UTF-8"?>
<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:cim="http://iec.ch/TC57/2013/CIM-schema-cim16#" xmlns:md="http://iec.ch/TC57/61970-552/ModelDescription/1#">
<md:FullModel rdf:about="urn:uuid:ssheb"><md:Model.profile>http://entsoe.eu/CIM/SteadyStateHypothesis/1/1</md:Model.profile></md:FullModel>
<cim:ExternalNetworkInjection rdf:about="#_eni"><cim:ExternalNetworkInjection.p>-50</cim:ExternalNetworkInjection.p><cim:ExternalNetworkInjection.q>-15</cim:ExternalNetworkInjection.q></cim:ExternalNetworkInjection>
<cim:EnergyConsumer rdf:about="#_load_b"><cim:EnergyConsumer.p>30</cim:EnergyConsumer.p><cim:EnergyConsumer.q>10</cim:EnergyConsumer.q></cim:EnergyConsumer>
<cim:EnergyConsumer rdf:about="#_load_d"><cim:EnergyConsumer.p>20</cim:EnergyConsumer.p><cim:EnergyConsumer.q>5</cim:EnergyConsumer.q></cim:EnergyConsumer>
<cim:Terminal rdf:about="#_t_ebo_2"><cim:ACDCTerminal.connected>false</cim:ACDCTerminal.connected></cim:Terminal>
</rdf:RDF>
"""

const _CGMES_SYNTH_EB_TP = """<?xml version="1.0" encoding="UTF-8"?>
<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:cim="http://iec.ch/TC57/2013/CIM-schema-cim16#" xmlns:md="http://iec.ch/TC57/61970-552/ModelDescription/1#">
<md:FullModel rdf:about="urn:uuid:tpeb"><md:Model.profile>http://entsoe.eu/CIM/Topology/4/1</md:Model.profile></md:FullModel>
<cim:TopologicalNode rdf:ID="_tn_a"><cim:IdentifiedObject.name>BUS_A</cim:IdentifiedObject.name><cim:TopologicalNode.BaseVoltage rdf:resource="#_bv110"/></cim:TopologicalNode>
<cim:TopologicalNode rdf:ID="_tn_b"><cim:IdentifiedObject.name>BUS_B</cim:IdentifiedObject.name><cim:TopologicalNode.BaseVoltage rdf:resource="#_bv110"/></cim:TopologicalNode>
<cim:TopologicalNode rdf:ID="_tn_c"><cim:IdentifiedObject.name>BUS_C</cim:IdentifiedObject.name><cim:TopologicalNode.BaseVoltage rdf:resource="#_bv110"/></cim:TopologicalNode>
<cim:TopologicalNode rdf:ID="_tn_d"><cim:IdentifiedObject.name>BUS_D</cim:IdentifiedObject.name><cim:TopologicalNode.BaseVoltage rdf:resource="#_bv220"/></cim:TopologicalNode>
<cim:Terminal rdf:about="#_t_eni"><cim:Terminal.TopologicalNode rdf:resource="#_tn_a"/></cim:Terminal>
<cim:Terminal rdf:about="#_t_load_b"><cim:Terminal.TopologicalNode rdf:resource="#_tn_b"/></cim:Terminal>
<cim:Terminal rdf:about="#_t_load_d"><cim:Terminal.TopologicalNode rdf:resource="#_tn_d"/></cim:Terminal>
<cim:Terminal rdf:about="#_t_lab_1"><cim:Terminal.TopologicalNode rdf:resource="#_tn_a"/></cim:Terminal>
<cim:Terminal rdf:about="#_t_lab_2"><cim:Terminal.TopologicalNode rdf:resource="#_tn_b"/></cim:Terminal>
<cim:Terminal rdf:about="#_t_ebs_1"><cim:Terminal.TopologicalNode rdf:resource="#_tn_b"/></cim:Terminal>
<cim:Terminal rdf:about="#_t_ebs_2"><cim:Terminal.TopologicalNode rdf:resource="#_tn_c"/></cim:Terminal>
<cim:Terminal rdf:about="#_t_eba_1"><cim:Terminal.TopologicalNode rdf:resource="#_tn_b"/></cim:Terminal>
<cim:Terminal rdf:about="#_t_eba_2"><cim:Terminal.TopologicalNode rdf:resource="#_tn_c"/></cim:Terminal>
<cim:Terminal rdf:about="#_t_ebc_1"><cim:Terminal.TopologicalNode rdf:resource="#_tn_c"/></cim:Terminal>
<cim:Terminal rdf:about="#_t_ebc_2"><cim:Terminal.TopologicalNode rdf:resource="#_tn_d"/></cim:Terminal>
<cim:Terminal rdf:about="#_t_ebl_1"><cim:Terminal.TopologicalNode rdf:resource="#_tn_a"/></cim:Terminal>
<cim:Terminal rdf:about="#_t_ebl_2"><cim:Terminal.TopologicalNode rdf:resource="#_tn_a"/></cim:Terminal>
<cim:Terminal rdf:about="#_t_ebo_1"><cim:Terminal.TopologicalNode rdf:resource="#_tn_a"/></cim:Terminal>
<cim:Terminal rdf:about="#_t_ebo_2"><cim:Terminal.TopologicalNode rdf:resource="#_tn_c"/></cim:Terminal>
</rdf:RDF>
"""

# Topological nodes in solver order, with their nominal voltage (kV).
const _CGMES_SYNTH_EB_TN = [("_tn_a", 110.0), ("_tn_b", 110.0), ("_tn_c", 110.0), ("_tn_d", 220.0)]

# Branches of the target model: terminals, node indices, z12 and z21 in ohm.
# EB_CROSS carries both on the base of terminal 1 and connects BUS_D through
# a per-unit ratio of 1.
function _cgmes_synth_eb_branches(mag::Float64)
  z21 = _cgmes_synth_eb_z21(mag)
  return [
    ("_t_lab_1", "_t_lab_2", 1, 2, complex(2.42, 12.1), complex(2.42, 12.1)),
    ("_t_ebs_1", "_t_ebs_2", 2, 3, complex(1.21, 6.05), complex(1.21, 6.05)),
    ("_t_eba_1", "_t_eba_2", 2, 3, complex(1.5, 7.0), z21[1]),
    ("_t_ebc_1", "_t_ebc_2", 3, 4, complex(0.6, 4.84), z21[2]),
  ]
end

# Solves the target model with its own Gauss-Seidel iteration, independent of
# the importer: BUS_A held at 113.3 kV and 0 degrees (the slack keeps the SV
# voltage of its bus), the loads of the SSH profile, voltages in kV referred
# to 110 kV (BUS_D through its nominal voltage), powers in MVA. Every branch
# is the two-port I1 = (V1 - V2)/z12, I2 = (V2 - V1)/z21. Returns the complex
# node voltages in that reference.
function _cgmes_synth_eb_state(mag::Float64)::Vector{ComplexF64}
  Y = zeros(ComplexF64, 4, 4)
  for (_, _, i, j, z12, z21) in _cgmes_synth_eb_branches(mag)
    Y[i, i] += 1 / z12
    Y[i, j] -= 1 / z12
    Y[j, j] += 1 / z21
    Y[j, i] -= 1 / z21
  end
  s = ComplexF64[0, -(30 + 10im), 0, -(20 + 5im)]
  V = ComplexF64[113.3, 110, 110, 110]
  for _ in 1:10_000
    dmax = 0.0
    for i in 2:4
      vi = (conj(s[i] / V[i]) - sum(Y[i, j] * V[j] for j in 1:4 if j != i)) / Y[i, i]
      dmax = max(dmax, abs(vi - V[i]))
      V[i] = vi
    end
    dmax < 1e-13 && return V
  end
  error("synthetic EquivalentBranch fixture: Gauss-Seidel did not converge")
end

function _cgmes_synth_eb_dir(; mag::Float64 = 1.0)::String
  dir = mktempdir()
  V = _cgmes_synth_eb_state(mag)
  io = IOBuffer()
  println(io, """<?xml version="1.0" encoding="UTF-8"?>
  <rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:cim="http://iec.ch/TC57/2013/CIM-schema-cim16#" xmlns:md="http://iec.ch/TC57/61970-552/ModelDescription/1#">
  <md:FullModel rdf:about="urn:uuid:sveb"><md:Model.profile>http://entsoe.eu/CIM/StateVariables/4/1</md:Model.profile></md:FullModel>""")
  for (k, (tn, vn)) in enumerate(_CGMES_SYNTH_EB_TN)
    v = abs(V[k]) * vn / 110.0
    a = rad2deg(angle(V[k]))
    println(io, "<cim:SvVoltage rdf:ID=\"_sv$(tn)\"><cim:SvVoltage.TopologicalNode rdf:resource=\"#$(tn)\"/><cim:SvVoltage.v>$(v)</cim:SvVoltage.v><cim:SvVoltage.angle>$(a)</cim:SvVoltage.angle></cim:SvVoltage>")
  end
  # flows in load sign: into the branch at each terminal, the loads as given,
  # the slack unit with the negative of its injection
  flows = Tuple{String,ComplexF64}[("_t_load_b", 30 + 10im), ("_t_load_d", 20 + 5im)]
  slack = 0.0im
  for (t1, t2, i, j, z12, z21) in _cgmes_synth_eb_branches(mag)
    i12 = (V[i] - V[j]) / z12
    push!(flows, (t1, V[i] * conj(i12)), (t2, V[j] * conj((V[j] - V[i]) / z21)))
    i == 1 && (slack += V[i] * conj(i12))
  end
  push!(flows, ("_t_eni", -slack))
  for (t, s) in flows
    println(io, "<cim:SvPowerFlow rdf:ID=\"_pf$(t)\"><cim:SvPowerFlow.Terminal rdf:resource=\"#$(t)\"/><cim:SvPowerFlow.p>$(real(s))</cim:SvPowerFlow.p><cim:SvPowerFlow.q>$(imag(s))</cim:SvPowerFlow.q></cim:SvPowerFlow>")
  end
  println(io, "</rdf:RDF>")
  write(joinpath(dir, "synth_eb_EQ.xml"), _cgmes_synth_eb_eq(mag))
  write(joinpath(dir, "synth_eb_SSH.xml"), _CGMES_SYNTH_EB_SSH)
  write(joinpath(dir, "synth_eb_TP.xml"), _CGMES_SYNTH_EB_TP)
  write(joinpath(dir, "synth_eb_SV.xml"), String(take!(io)))
  return dir
end
