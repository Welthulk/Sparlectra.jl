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

# file: src/adapters/cgmes/glue.jl
# purpose: adapter-side convenience overloads that dispatch core entry
#          points on CGMES adapter types (stage 5: the core block must not
#          reference adapter types at definition time, so this glue lives
#          in the adapters block)

"""
    runShortCircuit!(result::CGMESImporter.CGMESImportResult; kwargs...) -> ShortCircuitResult

Convenience overload: run the balanced IEC 60909-0 short-circuit sweep
directly on a CGMES import result (net plus harvested short-circuit
data). See the core `runShortCircuit!` methods for the keyword surface.
"""
runShortCircuit!(result::CGMESImporter.CGMESImportResult; kwargs...) = runShortCircuit!(result.net, result.shortcircuit; kwargs...)
