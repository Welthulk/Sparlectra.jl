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

case = "case30.m"      # or "case14.m"
flatstart = true
verbose = 1
println("Sparlectra version: ", Sparlectra.version(), "\n")
println("Importing MATPOWER case file: $case\n")
println("DEMO: Matpower_import_case14.jl\n  ------------------------------\n")


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

net = createNetFromMatPowerFile(filename = local_case, log = (verbose > 0), flatstart = flatstart)

# -----------------------------------------------------------------------------
# Run power flow
# -----------------------------------------------------------------------------
iters, status = runpf!(net, 50, 1e-8, 2; method=:rectangular, damp=0.2, opt_fd=true, opt_sparse=true, opt_flatstart=true)
print(showNet(net, verbose = true))
print("Power flow  $iters iterations with status $status\n")
run_acpflow(
    casefile = basename(local_case),   # run_acpflow expects name under data/mpower
    opt_fd = true,
    opt_sparse = true,
    method = :polar_full,
    opt_flatstart = flatstart,
    show_results = true,
    verbose = verbose,
)
