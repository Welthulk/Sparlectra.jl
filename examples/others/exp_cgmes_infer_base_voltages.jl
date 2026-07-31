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

# file: examples/others/exp_cgmes_infer_base_voltages.jl
# purpose: cgmes_import.infer_base_voltages — import a delivery whose
#          BaseVoltage catalog is missing (it lives in the absent boundary
#          EQ) by reconstructing nominal voltages from the SV state.

using Sparlectra

include(joinpath(@__DIR__, "..", "internal", "example_header.jl"))

# A two-bus delivery WITHOUT any BaseVoltage object: grid injection at BUS_1,
# load at BUS_2, one line. The SV voltages (112.2 / 108.9 kV) are the only
# voltage-level evidence — exactly the situation of a real TSO delivery
# imported without its boundary files.
const _INFER_DEMO_EQ = """<?xml version="1.0" encoding="UTF-8"?>
<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:cim="http://iec.ch/TC57/2013/CIM-schema-cim16#" xmlns:md="http://iec.ch/TC57/61970-552/ModelDescription/1#">
<md:FullModel rdf:about="urn:uuid:inferdemo-eq"><md:Model.profile>http://entsoe.eu/CIM/EquipmentCore/3/1</md:Model.profile></md:FullModel>
<cim:ExternalNetworkInjection rdf:ID="_eni"><cim:IdentifiedObject.name>GRID</cim:IdentifiedObject.name><cim:ExternalNetworkInjection.referencePriority>1</cim:ExternalNetworkInjection.referencePriority><cim:ExternalNetworkInjection.maxP>200</cim:ExternalNetworkInjection.maxP></cim:ExternalNetworkInjection>
<cim:Terminal rdf:ID="_t_eni"><cim:Terminal.ConductingEquipment rdf:resource="#_eni"/><cim:ACDCTerminal.sequenceNumber>1</cim:ACDCTerminal.sequenceNumber></cim:Terminal>
<cim:EnergyConsumer rdf:ID="_load"><cim:IdentifiedObject.name>LOAD</cim:IdentifiedObject.name></cim:EnergyConsumer>
<cim:Terminal rdf:ID="_t_load"><cim:Terminal.ConductingEquipment rdf:resource="#_load"/><cim:ACDCTerminal.sequenceNumber>1</cim:ACDCTerminal.sequenceNumber></cim:Terminal>
<cim:ACLineSegment rdf:ID="_line"><cim:IdentifiedObject.name>L_12</cim:IdentifiedObject.name><cim:ACLineSegment.r>1.2</cim:ACLineSegment.r><cim:ACLineSegment.x>9.7</cim:ACLineSegment.x></cim:ACLineSegment>
<cim:Terminal rdf:ID="_t_l1"><cim:Terminal.ConductingEquipment rdf:resource="#_line"/><cim:ACDCTerminal.sequenceNumber>1</cim:ACDCTerminal.sequenceNumber></cim:Terminal>
<cim:Terminal rdf:ID="_t_l2"><cim:Terminal.ConductingEquipment rdf:resource="#_line"/><cim:ACDCTerminal.sequenceNumber>2</cim:ACDCTerminal.sequenceNumber></cim:Terminal>
</rdf:RDF>
"""

const _INFER_DEMO_TP = """<?xml version="1.0" encoding="UTF-8"?>
<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:cim="http://iec.ch/TC57/2013/CIM-schema-cim16#" xmlns:md="http://iec.ch/TC57/61970-552/ModelDescription/1#">
<md:FullModel rdf:about="urn:uuid:inferdemo-tp"><md:Model.profile>http://entsoe.eu/CIM/Topology/4/1</md:Model.profile></md:FullModel>
<cim:TopologicalNode rdf:ID="_tn_1"><cim:IdentifiedObject.name>BUS_1</cim:IdentifiedObject.name></cim:TopologicalNode>
<cim:TopologicalNode rdf:ID="_tn_2"><cim:IdentifiedObject.name>BUS_2</cim:IdentifiedObject.name></cim:TopologicalNode>
<cim:Terminal rdf:about="#_t_eni"><cim:Terminal.TopologicalNode rdf:resource="#_tn_1"/></cim:Terminal>
<cim:Terminal rdf:about="#_t_load"><cim:Terminal.TopologicalNode rdf:resource="#_tn_2"/></cim:Terminal>
<cim:Terminal rdf:about="#_t_l1"><cim:Terminal.TopologicalNode rdf:resource="#_tn_1"/></cim:Terminal>
<cim:Terminal rdf:about="#_t_l2"><cim:Terminal.TopologicalNode rdf:resource="#_tn_2"/></cim:Terminal>
</rdf:RDF>
"""

const _INFER_DEMO_SSH = """<?xml version="1.0" encoding="UTF-8"?>
<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:cim="http://iec.ch/TC57/2013/CIM-schema-cim16#" xmlns:md="http://iec.ch/TC57/61970-552/ModelDescription/1#">
<md:FullModel rdf:about="urn:uuid:inferdemo-ssh"><md:Model.profile>http://entsoe.eu/CIM/SteadyStateHypothesis/1/1</md:Model.profile></md:FullModel>
<cim:Terminal rdf:about="#_t_eni"><cim:ACDCTerminal.connected>true</cim:ACDCTerminal.connected></cim:Terminal>
<cim:Terminal rdf:about="#_t_load"><cim:ACDCTerminal.connected>true</cim:ACDCTerminal.connected></cim:Terminal>
<cim:Terminal rdf:about="#_t_l1"><cim:ACDCTerminal.connected>true</cim:ACDCTerminal.connected></cim:Terminal>
<cim:Terminal rdf:about="#_t_l2"><cim:ACDCTerminal.connected>true</cim:ACDCTerminal.connected></cim:Terminal>
<cim:EnergyConsumer rdf:about="#_load"><cim:EnergyConsumer.p>40</cim:EnergyConsumer.p><cim:EnergyConsumer.q>12</cim:EnergyConsumer.q></cim:EnergyConsumer>
</rdf:RDF>
"""

const _INFER_DEMO_SV = """<?xml version="1.0" encoding="UTF-8"?>
<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:cim="http://iec.ch/TC57/2013/CIM-schema-cim16#" xmlns:md="http://iec.ch/TC57/61970-552/ModelDescription/1#">
<md:FullModel rdf:about="urn:uuid:inferdemo-sv"><md:Model.profile>http://entsoe.eu/CIM/StateVariables/4/1</md:Model.profile></md:FullModel>
<cim:SvVoltage rdf:ID="_sv_1"><cim:SvVoltage.TopologicalNode rdf:resource="#_tn_1"/><cim:SvVoltage.v>112.2</cim:SvVoltage.v><cim:SvVoltage.angle>0.0</cim:SvVoltage.angle></cim:SvVoltage>
<cim:SvVoltage rdf:ID="_sv_2"><cim:SvVoltage.TopologicalNode rdf:resource="#_tn_2"/><cim:SvVoltage.v>108.9</cim:SvVoltage.v><cim:SvVoltage.angle>-2.1</cim:SvVoltage.angle></cim:SvVoltage>
</rdf:RDF>
"""

"""
    main()

Write the boundary-less delivery to a temp folder. Without the option the
import aborts (no resolvable `BaseVoltage`, the analysis names the gap).
With `infer_base_voltages = true` both bus levels are reconstructed from
the SV voltages (112.2 / 108.9 kV snap to the 110 kV standard level), the
substitution is summarized as a warning, and the network solves.
"""
function main()
  print_example_banner("examples/others/exp_cgmes_infer_base_voltages.jl", "reconstruct missing nominal voltages from the SV state (cgmes_import.infer_base_voltages)")
  dir = mktempdir()
  write(joinpath(dir, "demo_EQ.xml"), _INFER_DEMO_EQ)
  write(joinpath(dir, "demo_TP.xml"), _INFER_DEMO_TP)
  write(joinpath(dir, "demo_SSH.xml"), _INFER_DEMO_SSH)
  write(joinpath(dir, "demo_SV.xml"), _INFER_DEMO_SV)

  aborted = try
    importCGMES(path = dir, name = "infer_demo_off")
    false
  catch
    true
  end

  res = importCGMES(path = dir, name = "infer_demo_on", infer_base_voltages = true)
  vns = sort!(unique(Sparlectra.getNodeVn(b) for b in res.net.nodeVec))
  warning = only(m for m in res.messages if occursin("inferred base voltages", m))
  ite, erg = runpf!(res.net, 30, 1e-8, 0)
  return (aborted_without_option = aborted, bus_levels_kV = vns, warning = warning, converged = erg == 0, iterations = ite)
end

result = run_example(main)
println()
println("without option: import aborted = ", result.aborted_without_option)
println("with option:    bus levels ", result.bus_levels_kV, " kV, converged = ", result.converged, " (", result.iterations, " iterations)")
println(result.warning)
