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

# file: tools/sysimage_workload.jl
# purpose: precompile execution workload for the fast-start sysimage build
#          (tools/build_sysimage.jl). Exercises the interactive paths beyond
#          the shipped PrecompileTools workload: Web UI start with the
#          reserved warm-up machinery (no run-history entries), one page
#          request through the real socket handler, one full power-flow
#          service run, one N-1 contingency service run (the Contingency button
#          path), one full short-circuit service run (CGMES MiniGrid),
#          one run_sparlectra call on case14, then a clean shutdown. Never
#          run this file directly; PackageCompiler executes it in a child
#          process during the build.

using Sparlectra
using Sockets

# Load AnalyticLoadFlow exactly like start_webui.jl does when it is
# installed. This matters beyond the APSLF paths themselves: a package
# loaded AFTER the sysimage invalidates precompiled methods inside the
# image, and the first run then recompiles them (measured: 36 s instead of
# 1 s for the first case118 service run). Loading it during the build means
# the image is compiled in the post-AnalyticLoadFlow method world, so the
# runtime load becomes a no-op for compilation.
if Base.find_package("AnalyticLoadFlow") !== nothing
  @eval using AnalyticLoadFlow
end

# nonstandard ports so a Web UI the maintainer has open on 8080 does not
# make the build silently skip the server paths; the first free candidate
# wins (a REPL session may hold an earlier workload port)
const _WORKLOAD_PORTS = 8091:8097

# One service run into a throwaway output root; config_file and output_root
# are filled in here so call sites only state what differs. Failures are
# logged, not fatal: a workload gap costs first-click latency, not the build.
function _workload_service_run(request::Dict{String,Any}; label::String)
  mktempdir() do outdir
    request["config_file"] = Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH
    request["output_root"] = outdir
    result = start_powerflow_run(request)
    status = String(get(result, "status", ""))
    status == "failed" && @warn "sysimage workload: $(label) service run failed" message = get(result, "message", "?")
    run_id = get(result, "run_id", nothing)
    run_id === nothing || get_powerflow_result(String(run_id))
  end
  return nothing
end

function _workload_request(port::Int, target::String)
  sock = Sockets.connect("127.0.0.1", port)
  try
    write(sock, string("GET ", target, " HTTP/1.1\r\nHost: 127.0.0.1:", port, "\r\nConnection: close\r\n\r\n"))
    return String(read(sock))
  finally
    close(sock)
  end
end

function run_workload()
  server = nothing
  port = 0
  for candidate in _WORKLOAD_PORTS
    server = try
      start_sparlectra_webui(open_browser = false, port = candidate, warmup = true, warmup_store_result = false)
    catch err
      err isa ArgumentError && occursin("already in use", sprint(showerror, err)) ? nothing : rethrow()
    end
    if server !== nothing
      port = candidate
      break
    end
  end
  server === nothing && error("sysimage workload: no free port in $(_WORKLOAD_PORTS)")
  try
    # the warm-up waits for the first served page; request it through the
    # real socket path so the handler chain is part of the compile trace
    _workload_request(port, "/")
    deadline = time() + 180.0
    while server.runtime.warmup_state !== :done && time() < deadline
      sleep(0.2)
    end
    server.runtime.warmup_state === :done || @warn "sysimage workload: warm-up did not finish within 180 s; continuing with the paths compiled so far"
    # after the warm-up the same route serves the real PowerFlow form; that
    # render path must be part of the trace, otherwise the first form view
    # after a fast start pays it as JIT time (measured ~11 s)
    _workload_request(port, "/powerflow")
    _workload_request(port, "/webui/fast-start")
    # the N-1 weights editor seeds a table from the case's element names (builds
    # the net); trace it so the first "edit N-1 weights" open is not JIT. case14
    # must live in the server's case directory for the seeded-name path.
    try
      cp(Sparlectra.ensure_casefile("case14.m"), joinpath(server.runtime.case_directory, "case14.m"); force = true)
    catch err
      @warn "sysimage workload: could not stage case14 for the weights editor" exception = err
    end
    _workload_request(port, "/powerflow/contingency-weights?case=case14.m")
    # real service runs through the full pipeline (case resolution,
    # artifact writers for the CSVs, result.json, run index, result
    # lookup): the first real Web UI run otherwise pays exactly this
    # compilation. Temp output roots keep the user's run history clean.
    # 1. MATPOWER power flow (case14).
    _workload_service_run(Dict{String,Any}(
      "casefile" => Sparlectra.ensure_casefile("case14.m"),
      "config_overrides" => Dict{String,Any}("output.logfile_results" => "off", "benchmark.enabled" => false),
    ); label = "matpower power-flow")
    # N-1 contingency (#331 Phase 5): the "Contingency (N-1)" button path through
    # runContingencies!, the CSV/report writers, and the result registry, so the
    # first click after a fast start is not paid as JIT. Branch kind on case14
    # covers the shared code; the generator kind reuses the same service.
    _workload_service_run(Dict{String,Any}(
      "casefile" => Sparlectra.ensure_casefile("case14.m"),
      "contingency_mode" => true,
      "contingency_kind" => "branch",
    ); label = "matpower contingency n-1")
    # The MiniGrid CGMES delivery covers the CGMES compile paths (ZIP and
    # XML reading, profile harvesting, net construction, control mapping).
    # It comes from the same case cache the Web UI uses; when absent it is
    # fetched once through the existing registry (network). A failed fetch
    # only skips these traces, never the build, but the gap is logged so
    # the coverage loss is visible.
    sc_zip = joinpath(server.runtime.case_directory, "cgmes_minigrid.zip")
    if !isfile(sc_zip)
      sc_zip = try
        Sparlectra.CGMESImporter.fetchCGMESTestSet("minigrid"; outdir = server.runtime.case_directory)
      catch err
        @warn "sysimage workload: CGMES service paths NOT traced (MiniGrid fetch failed; the first CGMES run will pay JIT)" exception = err
        nothing
      end
    end
    if sc_zip !== nothing
      # 2. CGMES power flow: import plus the PF branch specific to CGMES
      # nets (tap/machine control mapping, SV handling).
      _workload_service_run(Dict{String,Any}(
        "casefile" => sc_zip,
        "config_overrides" => Dict{String,Any}("output.logfile_results" => "off", "benchmark.enabled" => false),
      ); label = "cgmes power-flow")
      # 3. CGMES short circuit: Z-bus solve for both c-factor cases, CSV
      # artifact writers, coverage report (the Short-circuit button path).
      _workload_service_run(Dict{String,Any}(
        "casefile" => sc_zip,
        "short_circuit_mode" => true,
      ); label = "cgmes short-circuit")
    end
    # one full interactive solve (the case is pre-fetched by the build
    # script); logfile output stays off so the workload never leaves
    # run_case*.log files in examples/_out, regardless of any user
    # configuration
    workload_cfg = SparlectraConfig(output = OutputConfig(logfile_results = :off))
    result = run_sparlectra(casefile = "case14.m", config = workload_cfg)
    result.final_converged || @warn "sysimage workload: case14 run did not converge" outcome = result.outcome
  finally
    close(server)
  end
  return nothing
end

Base.invokelatest(run_workload)
