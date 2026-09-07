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

# file: tools/app_workload_runtime.jl
# purpose: precompile execution workload for the RUNTIME app flavor
#          (buildApp(flavor = :runtime)): the library paths only, no Web UI
#          server, no sockets, no network. Power flow on the embedded
#          warm-up cases small to large, state estimation with bad-data
#          diagnostics, an N-1 sweep, losses, and the result printers.
#          Never run this file directly; PackageCompiler executes it in a
#          child process during the build.

using Sparlectra

# measurement guard (RP1 worktree incident 2026-09-03): this script must run
# against THE repository checkout it lives in, never a stale worktree or
# another depot copy picked up through a wrong --project
let expected = normpath(joinpath(@__DIR__, "..", "src")), actual = normpath(String(pathof(Sparlectra)))
  startswith(actual, expected) || error("tools guard: Sparlectra loaded from " * actual * ", expected under " * expected * "; start julia with --project=" * normpath(joinpath(@__DIR__, "..")))
end

if Base.find_package("AnalyticLoadFlow") !== nothing
  @eval using AnalyticLoadFlow
end

function run_runtime_workload()
  pkgroot = pkgdir(Sparlectra)
  nets = Sparlectra.Net[]
  for name in ("warmup_case3.jl", "warmup_case14.jl", "warmup_case118.jl")
    path = joinpath(pkgroot, "data", "webui", name)
    isfile(path) || (@warn "runtime workload: warm-up case missing" path; continue)
    nt = include(path)
    net = Sparlectra.createNetFromMatPowerCase(mpc = nt)
    ite, erg = runpf!(net, 60, 1e-8, 0; method = :rectangular)
    erg == 0 || @warn "runtime workload: warm-up case did not converge" name
    push!(nets, net)
  end
  isempty(nets) && error("runtime workload: no warm-up case solved")
  # losses + human-facing printers on the mid-size case
  net = length(nets) >= 2 ? nets[2] : nets[end]
  calcNetLosses!(net)
  buildACPFlowReport(net; ct = 0.0, ite = 1, converged = true)
  redirect_stdout(devnull) do
    printACPFlowResults(net, 0.0, 1, 1e-8)
  end
  # state estimation with diagnostics and the SE printer
  try
    append!(net.measurements, generateMeasurementsFromPF(net; noise = false))
    runse!(net)
    diag = runse_diagnostics(net; max_eliminations = 0)
    print_se_diagnostics(devnull, diag.diagnostics)
  catch err
    @warn "runtime workload: SE trace skipped" exception = err
  end
  # N-1 branch sweep on the mid-size case
  try
    cases = Sparlectra.generateN1Branches(net)
    runContingencies!(net, cases)
  catch err
    @warn "runtime workload: N-1 trace skipped" exception = err
  end
  return nothing
end

Base.invokelatest(run_runtime_workload)

# Test-suite trace: run the fast-profile test cases so everything they touch
# is compiled into the app as well. The webui group is skipped explicitly,
# this flavor ships no GUI. A test failure must not abort the build: warn
# and continue, the suite is a trace here, not a gate (the gates run in CI
# and the developer workflow).
try
  ENV["SPARLECTRA_TEST_PROFILE"] = "fast"
  # webui left the fast profile in the risk-based resort, so the runtime
  # flavor needs no group skip any more
  Base.invokelatest(include, joinpath(pkgdir(Sparlectra), "test", "runtests.jl"))
catch err
  @warn "runtime workload: test-suite trace failed (build continues)" exception = err
end
