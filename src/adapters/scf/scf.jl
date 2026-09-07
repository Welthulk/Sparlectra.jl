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

# file: src/adapters/scf/scf.jl
# purpose: include hub for the Sparlectra Case Format (SCF, issue #342):
#          the PGM-compatible single-file case format: writer (exportSCF)
#          and reader (importSCF), sharing the JSON encoding module.

include("scf_json.jl")
include("scf_case.jl")
include("scf_controllers.jl")
include("scf_export.jl")
include("scf_import.jl")
