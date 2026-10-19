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

# file: examples/others/exp_cgmes_energy_source.jl
# purpose: cgmes_import.map_energy_source: a delivery that carries its
#          generation as EnergySource, imported with the key off and on,
#          and the same park as a CGMES 3.0 PowerElectronicsConnection
#          with cgmes_import.map_power_electronics, written to CGMES
#          2.4.15 as EquivalentInjection and imported again.

using Sparlectra

include(joinpath(@__DIR__, "..", "others", "example_header.jl"))

# Two buses at 110 kV: the grid infeed GRID at BUS_A (the slack), a 30 MW
# park PV_PARK modelled as EnergySource and a 40 MW load at BUS_B, one line.
# CGMES states EnergySource.activePower in load sign convention, so the
# generating park carries -30 MW in the SSH and in its SvPowerFlow.
const _ES_DEMO_EQ = """<?xml version="1.0" encoding="UTF-8"?>
<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:cim="http://iec.ch/TC57/2013/CIM-schema-cim16#" xmlns:md="http://iec.ch/TC57/61970-552/ModelDescription/1#">
<md:FullModel rdf:about="urn:uuid:esdemo-eq"><md:Model.profile>http://entsoe.eu/CIM/EquipmentCore/3/1</md:Model.profile></md:FullModel>
<cim:BaseVoltage rdf:ID="_bv110"><cim:BaseVoltage.nominalVoltage>110</cim:BaseVoltage.nominalVoltage></cim:BaseVoltage>
<cim:ExternalNetworkInjection rdf:ID="_eni"><cim:IdentifiedObject.name>GRID</cim:IdentifiedObject.name><cim:ExternalNetworkInjection.referencePriority>1</cim:ExternalNetworkInjection.referencePriority><cim:ExternalNetworkInjection.maxP>200</cim:ExternalNetworkInjection.maxP></cim:ExternalNetworkInjection>
<cim:Terminal rdf:ID="_t_eni"><cim:Terminal.ConductingEquipment rdf:resource="#_eni"/><cim:ACDCTerminal.sequenceNumber>1</cim:ACDCTerminal.sequenceNumber></cim:Terminal>
<cim:EnergySource rdf:ID="_es"><cim:IdentifiedObject.name>PV_PARK</cim:IdentifiedObject.name></cim:EnergySource>
<cim:Terminal rdf:ID="_t_es"><cim:Terminal.ConductingEquipment rdf:resource="#_es"/><cim:ACDCTerminal.sequenceNumber>1</cim:ACDCTerminal.sequenceNumber></cim:Terminal>
<cim:EnergyConsumer rdf:ID="_load"><cim:IdentifiedObject.name>LOAD_B</cim:IdentifiedObject.name></cim:EnergyConsumer>
<cim:Terminal rdf:ID="_t_load"><cim:Terminal.ConductingEquipment rdf:resource="#_load"/><cim:ACDCTerminal.sequenceNumber>1</cim:ACDCTerminal.sequenceNumber></cim:Terminal>
<cim:ACLineSegment rdf:ID="_line"><cim:IdentifiedObject.name>L_AB</cim:IdentifiedObject.name><cim:ACLineSegment.r>2.42</cim:ACLineSegment.r><cim:ACLineSegment.x>12.1</cim:ACLineSegment.x><cim:ACLineSegment.bch>0</cim:ACLineSegment.bch><cim:ConductingEquipment.BaseVoltage rdf:resource="#_bv110"/></cim:ACLineSegment>
<cim:Terminal rdf:ID="_t_l1"><cim:Terminal.ConductingEquipment rdf:resource="#_line"/><cim:ACDCTerminal.sequenceNumber>1</cim:ACDCTerminal.sequenceNumber></cim:Terminal>
<cim:Terminal rdf:ID="_t_l2"><cim:Terminal.ConductingEquipment rdf:resource="#_line"/><cim:ACDCTerminal.sequenceNumber>2</cim:ACDCTerminal.sequenceNumber></cim:Terminal>
</rdf:RDF>
"""

const _ES_DEMO_TP = """<?xml version="1.0" encoding="UTF-8"?>
<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:cim="http://iec.ch/TC57/2013/CIM-schema-cim16#" xmlns:md="http://iec.ch/TC57/61970-552/ModelDescription/1#">
<md:FullModel rdf:about="urn:uuid:esdemo-tp"><md:Model.profile>http://entsoe.eu/CIM/Topology/4/1</md:Model.profile></md:FullModel>
<cim:TopologicalNode rdf:ID="_tn_a"><cim:IdentifiedObject.name>BUS_A</cim:IdentifiedObject.name><cim:TopologicalNode.BaseVoltage rdf:resource="#_bv110"/></cim:TopologicalNode>
<cim:TopologicalNode rdf:ID="_tn_b"><cim:IdentifiedObject.name>BUS_B</cim:IdentifiedObject.name><cim:TopologicalNode.BaseVoltage rdf:resource="#_bv110"/></cim:TopologicalNode>
<cim:Terminal rdf:about="#_t_eni"><cim:Terminal.TopologicalNode rdf:resource="#_tn_a"/></cim:Terminal>
<cim:Terminal rdf:about="#_t_es"><cim:Terminal.TopologicalNode rdf:resource="#_tn_b"/></cim:Terminal>
<cim:Terminal rdf:about="#_t_load"><cim:Terminal.TopologicalNode rdf:resource="#_tn_b"/></cim:Terminal>
<cim:Terminal rdf:about="#_t_l1"><cim:Terminal.TopologicalNode rdf:resource="#_tn_a"/></cim:Terminal>
<cim:Terminal rdf:about="#_t_l2"><cim:Terminal.TopologicalNode rdf:resource="#_tn_b"/></cim:Terminal>
</rdf:RDF>
"""

const _ES_DEMO_SSH = """<?xml version="1.0" encoding="UTF-8"?>
<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:cim="http://iec.ch/TC57/2013/CIM-schema-cim16#" xmlns:md="http://iec.ch/TC57/61970-552/ModelDescription/1#">
<md:FullModel rdf:about="urn:uuid:esdemo-ssh"><md:Model.profile>http://entsoe.eu/CIM/SteadyStateHypothesis/1/1</md:Model.profile></md:FullModel>
<cim:ExternalNetworkInjection rdf:about="#_eni"><cim:ExternalNetworkInjection.p>-10</cim:ExternalNetworkInjection.p><cim:ExternalNetworkInjection.q>-7</cim:ExternalNetworkInjection.q></cim:ExternalNetworkInjection>
<cim:EnergySource rdf:about="#_es"><cim:EnergySource.activePower>-30</cim:EnergySource.activePower><cim:EnergySource.reactivePower>-5</cim:EnergySource.reactivePower></cim:EnergySource>
<cim:EnergyConsumer rdf:about="#_load"><cim:EnergyConsumer.p>40</cim:EnergyConsumer.p><cim:EnergyConsumer.q>12</cim:EnergyConsumer.q></cim:EnergyConsumer>
<cim:Terminal rdf:about="#_t_eni"><cim:ACDCTerminal.connected>true</cim:ACDCTerminal.connected></cim:Terminal>
<cim:Terminal rdf:about="#_t_es"><cim:ACDCTerminal.connected>true</cim:ACDCTerminal.connected></cim:Terminal>
<cim:Terminal rdf:about="#_t_load"><cim:ACDCTerminal.connected>true</cim:ACDCTerminal.connected></cim:Terminal>
<cim:Terminal rdf:about="#_t_l1"><cim:ACDCTerminal.connected>true</cim:ACDCTerminal.connected></cim:Terminal>
<cim:Terminal rdf:about="#_t_l2"><cim:ACDCTerminal.connected>true</cim:ACDCTerminal.connected></cim:Terminal>
</rdf:RDF>
"""

const _ES_DEMO_SV = """<?xml version="1.0" encoding="UTF-8"?>
<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:cim="http://iec.ch/TC57/2013/CIM-schema-cim16#" xmlns:md="http://iec.ch/TC57/61970-552/ModelDescription/1#">
<md:FullModel rdf:about="urn:uuid:esdemo-sv"><md:Model.profile>http://entsoe.eu/CIM/StateVariables/4/1</md:Model.profile></md:FullModel>
<cim:SvVoltage rdf:ID="_sv_a"><cim:SvVoltage.TopologicalNode rdf:resource="#_tn_a"/><cim:SvVoltage.v>113.3</cim:SvVoltage.v><cim:SvVoltage.angle>0</cim:SvVoltage.angle></cim:SvVoltage>
<cim:SvVoltage rdf:ID="_sv_b"><cim:SvVoltage.TopologicalNode rdf:resource="#_tn_b"/><cim:SvVoltage.v>112.0</cim:SvVoltage.v><cim:SvVoltage.angle>-0.6</cim:SvVoltage.angle></cim:SvVoltage>
<cim:SvPowerFlow rdf:ID="_svpf_es"><cim:SvPowerFlow.Terminal rdf:resource="#_t_es"/><cim:SvPowerFlow.p>-30</cim:SvPowerFlow.p><cim:SvPowerFlow.q>-5</cim:SvPowerFlow.q></cim:SvPowerFlow>
</rdf:RDF>
"""

# The same delivery in CGMES 3.0 with the park as PowerElectronicsConnection:
# namespace and profile URIs of 3.0, p and q in place of activePower and
# reactivePower (same load sign). Applied in this order, so the class name
# is replaced only after the attribute names.
const _PEC_DEMO_RENAMES = [
  "http://iec.ch/TC57/2013/CIM-schema-cim16#" => "http://iec.ch/TC57/CIM100#",
  "http://entsoe.eu/CIM/EquipmentCore/3/1" => "http://iec.ch/TC57/ns/CIM/CoreEquipment-EU/3.0",
  "http://entsoe.eu/CIM/Topology/4/1" => "http://iec.ch/TC57/ns/CIM/Topology-EU/3.0",
  "http://entsoe.eu/CIM/SteadyStateHypothesis/1/1" => "http://iec.ch/TC57/ns/CIM/SteadyStateHypothesis-EU/3.0",
  "http://entsoe.eu/CIM/StateVariables/4/1" => "http://iec.ch/TC57/ns/CIM/StateVariables-EU/3.0",
  "EnergySource.activePower" => "PowerElectronicsConnection.p",
  "EnergySource.reactivePower" => "PowerElectronicsConnection.q",
  "cim:EnergySource" => "cim:PowerElectronicsConnection",
]
_pec_demo_xml(xml::AbstractString) = foldl((s, r) -> replace(s, r), _PEC_DEMO_RENAMES; init = String(xml))

# import with the keys given, solve, and read what the slack has to cover;
# `class` names the class whose import notice is returned: the one that
# names the key when the class is off, the one that counts the mapped
# objects when it is on
function _es_demo_run(dir::AbstractString, name::AbstractString; map_energy_source::Bool = false, map_power_electronics::Bool = false, class::AbstractString = "EnergySource")
  res = importCGMES(path = dir, name = name, map_energy_source = map_energy_source, map_power_electronics = map_power_electronics)
  notice = only(m for m in res.messages if startswith(m, "notice: $(class) not mapped, switch") || occursin("$(class) object(s) mapped", m))
  park = get(res.element_by_mrid, "_es", nothing)
  park_p = park === nothing ? nothing : res.net.prosumpsVec[last(park)].pVal
  iterations, erg = redirect_stdout(devnull) do
    runpf!(res.net, 30, 1e-8, 0)
  end
  slack_p = res.net.nodeVec[res.net.busDict[res.slack_bus]]._pƩGen
  return (notice = notice, park_p = park_p, converged = erg == 0, iterations = iterations, slack_p = slack_p)
end

# the CGMES 3.0 import with the park as PowerElectronicsConnection, written
# as CGMES 2.4.15 (no class for it: EquivalentInjection with the
# description PowerElectronicsConnection) and imported again
function _pec_export_demo(dir30::AbstractString)
  res = importCGMES(path = dir30, name = "pec_demo_export", map_power_electronics = true)
  redirect_stdout(devnull) do
    runpf!(res.net, 30, 1e-8, 0)
  end
  out = mktempdir()
  notices = String[]
  files = writeCGMESFiles(res.net; path = out, notices = notices)
  back = importCGMES(path = out, name = "pec_demo_reimport")
  eis = collect(Sparlectra.CGMESImporter.objectsOf(back.store, :EquivalentInjection))
  park = back.net.prosumpsVec[last(back.element_by_mrid[only(eis).mrid])]
  return (notice = only(notices), classes_in_eq = count("<cim:EquivalentInjection rdf:ID", read(files[1], String)), park_type = park.comp.cTyp, park_p = park.pVal)
end

"""
    main()

Write the two-bus delivery to a temp folder and import it twice. With the
default `map_energy_source = false` the park is left out, one notice names
the key, and the slack covers the whole load. With `map_energy_source = true` the
park is a 30 MW PQ injection (`p = -activePower`, equal to `-SvPowerFlow.p`
at its terminal) and the slack covers only the rest. A third import reads
the same delivery in CGMES 3.0 with the park as `PowerElectronicsConnection`
and `map_power_electronics = true`, which maps it the same way
(`p = -PowerElectronicsConnection.p`). That import is then written as
CGMES 2.4.15, where the park becomes an `EquivalentInjection`, and read
again: the same P, and the converter type comes back from its description.
"""
function main()
  print_example_banner("examples/others/exp_cgmes_energy_source.jl", "EnergySource generation in a CGMES delivery (cgmes_import.map_energy_source)")
  dir = mktempdir()
  write(joinpath(dir, "demo_EQ.xml"), _ES_DEMO_EQ)
  write(joinpath(dir, "demo_TP.xml"), _ES_DEMO_TP)
  write(joinpath(dir, "demo_SSH.xml"), _ES_DEMO_SSH)
  write(joinpath(dir, "demo_SV.xml"), _ES_DEMO_SV)
  dir30 = mktempdir()
  for (file, xml) in (("demo_EQ.xml", _ES_DEMO_EQ), ("demo_TP.xml", _ES_DEMO_TP), ("demo_SSH.xml", _ES_DEMO_SSH), ("demo_SV.xml", _ES_DEMO_SV))
    write(joinpath(dir30, file), _pec_demo_xml(xml))
  end
  return (
    off = _es_demo_run(dir, "es_demo_off"),
    on = _es_demo_run(dir, "es_demo_on"; map_energy_source = true),
    pec = _es_demo_run(dir30, "pec_demo_on"; map_power_electronics = true, class = "PowerElectronicsConnection"),
    pec_export = _pec_export_demo(dir30),
  )
end

result = run_example(main)
println()
for (label, r) in (("map_energy_source = false", result.off), ("map_energy_source = true ", result.on), ("CGMES 3.0, map_power_electronics = true", result.pec))
  println(label, ": park P = ", r.park_p === nothing ? "not mapped" : string(r.park_p, " MW"), ", slack P = ", round(r.slack_p; digits = 2), " MW, converged = ", r.converged, " (", r.iterations, " iterations)")
  println("  ", r.notice)
end
println()
println("CGMES 2.4.15 export of the CGMES 3.0 import: ", result.pec_export.classes_in_eq, " EquivalentInjection in the EQ")
println("  ", result.pec_export.notice)
println("  re-imported park: P = ", result.pec_export.park_p, " MW, component type ", result.pec_export.park_type)
