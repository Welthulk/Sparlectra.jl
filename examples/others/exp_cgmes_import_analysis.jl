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

# file: examples/others/exp_cgmes_import_analysis.jl
# purpose: analyzeCGMES on a deliberately incomplete delivery — the report
#          names the missing declared dependency (the boundary-set case)
#          instead of a bare unresolved-reference count.

using Sparlectra

include(joinpath(@__DIR__, "..", "internal", "example_header.jl"))

# A minimal EQ file whose header declares a prerequisite model that is not
# part of the input — exactly how a real delivery references its boundary set.
const _ANALYSIS_DEMO_EQ = """<?xml version="1.0" encoding="UTF-8"?>
<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:cim="http://iec.ch/TC57/2013/CIM-schema-cim16#" xmlns:md="http://iec.ch/TC57/61970-552/ModelDescription/1#">
<md:FullModel rdf:about="urn:uuid:demo-eq-model">
<md:Model.profile>http://entsoe.eu/CIM/EquipmentCore/3/1</md:Model.profile>
<md:Model.DependentOn rdf:resource="urn:uuid:demo-boundary-model"/>
</md:FullModel>
<cim:ACLineSegment rdf:ID="_demo_line">
<cim:IdentifiedObject.name>DemoLine</cim:IdentifiedObject.name>
<cim:ConductingEquipment.BaseVoltage rdf:resource="#_bv_in_missing_boundary"/>
</cim:ACLineSegment>
</rdf:RDF>
"""

"""
    main()

Write the incomplete delivery to a temp folder and run `analyzeCGMES` on it.
The printed report lists the supplied models, calls out the missing
`md:Model.DependentOn` prerequisite by id, aggregates the unresolved
references, and closes with a plain-language verdict. The same report is
printed automatically when `importCGMES` aborts on such a delivery.
"""
function main()
  print_example_banner("examples/others/exp_cgmes_import_analysis.jl", "analyzeCGMES: name the missing boundary dependency of an incomplete CGMES delivery")
  dir = mktempdir()
  write(joinpath(dir, "demo_EQ.xml"), _ANALYSIS_DEMO_EQ)
  report = Sparlectra.CGMESImporter.analyzeCGMES(path = dir)
  return (missing_named = occursin("demo-boundary-model", report), has_verdict = occursin("Verdict:", report))
end

result = run_example(main)
println()
println("missing dependency named: ", result.missing_named, ", verdict present: ", result.has_verdict)
