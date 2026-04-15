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

"""
Example: import a (parameterizable) MATPOWER case and run AC load flow.

Usage examples:
- julia --project=. src/examples/import_case_xxxyyyy_loadflow.jl
- julia --project=. src/examples/import_case_xxxyyyy_loadflow.jl case14.m
- julia --project=. src/examples/import_case_xxxyyyy_loadflow.jl xxxyyyy.m
"""

using Sparlectra

const DEFAULT_CASEFILE = "case14.m"
const show_results = false

function run_import_case_xxxyyyy_loadflow_example(; casefile::AbstractString = DEFAULT_CASEFILE, method::Symbol = :rectangular, max_ite::Int = 30, tol::Float64 = 1e-6, verbose::Int = 1)
  local_case = Sparlectra.FetchMatpowerCase.ensure_casefile(String(casefile); outdir = Sparlectra.MPOWER_DIR, to_jl = false, overwrite = false)

  println("Loaded case file: ", local_case)
  println("Running AC power flow with method=", method, ", max_ite=", max_ite, ", tol=", tol)
  # Time network creation
  net_creation_time = @elapsed begin
    net = createNetFromMatPowerFile(filename = local_case, log = (verbose > 0), flatstart = false)
  end

  # Time power flow execution  
  pf_time = @elapsed begin
    ite, erg, etime = run_net_acpflow(net = net, max_ite = max_ite, tol = tol, verbose = verbose, opt_sparse = true, opt_fd = false, method = method, show_results = true)
  end

  converged = (erg == 0)
  println("Summary: converged=", converged, ", iterations=", ite, ", solver_time=", round(etime; digits = 6), ", total_pf_time=", round(pf_time; digits = 6), " net_creation_time=", round(net_creation_time; digits = 6), " seconds  ")

  return nothing
  #return (net = net, converged = converged, iterations = ite, elapsed_s = etime)
end

function main()
  casefile = isempty(ARGS) ? DEFAULT_CASEFILE : ARGS[1]
  Base.invokelatest(() -> run_import_case_xxxyyyy_loadflow_example(casefile = casefile))
end

main()
