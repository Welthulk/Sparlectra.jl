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

using Sparlectra
import Sparlectra: MatpowerIO

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

case = "case14.jl"      # or "case14.m"
flatstart = false

# -----------------------------------------------------------------------------
# Ensure case is available locally (download on demand)
# -----------------------------------------------------------------------------

# This will:
# - download case14.m into data/mpower/ if missing
# - generate case14.jl if requested and missing
local_case = Sparlectra.FetchMatpowerCase.ensure_casefile(
    case;
    outdir = Sparlectra.MPOWER_DIR,
    to_jl = true,
    overwrite = false,
)

# -----------------------------------------------------------------------------
# Read case (Julia or MATPOWER)
# -----------------------------------------------------------------------------

mpc = MatpowerIO.read_case(local_case)

# -----------------------------------------------------------------------------
# Run power flow
# -----------------------------------------------------------------------------

run_acpflow(
    casefile = basename(local_case),   # run_acpflow expects name under data/mpower
    opt_fd = true,
    opt_sparse = true,
    method = :polar_full,
    opt_flatstart = flatstart,
)
