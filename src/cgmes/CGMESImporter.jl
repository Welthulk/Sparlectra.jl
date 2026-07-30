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

# file: src/cgmes/CGMESImporter.jl
# purpose: lean CGMES (ENTSO-E CIM) import for Sparlectra — module root.
# Stage 0: container, generic RDF/XML reader, merged store, summarizeCGMES.
# Concept: task_concept_cgmes_import.md; plan: task_plan_cgmes_import.md.

"""
    CGMESImporter

Lean CGMES 2.4.15 / 3.0 import layer. Stage 0 provides reading and
diagnostics (`loadCGMES`, `summarizeCGMES`); Net mapping follows in later
stages.
"""
module CGMESImporter

import EzXML
import ZipArchives
import Downloads
import ..Sparlectra

export CGMESFile, CIMObject, CGMESFileInfo, CGMESStore, CGMESSummary
export collectCGMESFiles, loadCGMES, summarizeCGMES
export objectsOf, countOf, num, str, boolval, enumval, ref, refsAll, unresolvedReferences
export CGMESTopology, buildTopology
export fetchCGMESTestSet, ensureCGMESTestConfigurations, CGMES_TESTSET_ALIASES
export fetchReliCapGridSet, RELICAPGRID_ALIASES, RELICAPGRID_COMBINED, allCGMESTestSetAliases
export CGMESImportResult, CGMESShortCircuitData, importCGMES, createNetFromCGMES, compareWithSV, shortCircuitCoverage, printShortCircuitCoverage

include("cgmes_schema.jl")
include("cgmes_container.jl")
include("cgmes_reader.jl")
include("cgmes_store.jl")
include("cgmes_topology.jl")
include("cgmes_mapping.jl")
include("cgmes_report.jl")
include("cgmes_testsets.jl")

end # module CGMESImporter
