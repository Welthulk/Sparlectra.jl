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
# purpose: precompile execution workload for the sysimage build
#          (tools/build_sysimage.jl). Exercises the interactive paths beyond
#          the shipped PrecompileTools workload: Web UI start with the
#          reserved warm-up machinery (no run-history entries), one page
#          request through the real socket handler, one full power-flow
#          service run, one N-1 contingency service run (the Contingency button
#          path), one full short-circuit service run (CGMES MiniGrid),
#          one run_sparlectra call on the tracked MATPOWER case, and a power flow plus state
#          estimation for EVERY import format (MATPOWER, SCF, CGMES, DTF),
#          then a clean shutdown. Never
#          run this file directly; PackageCompiler executes it in a child
#          process during the build.

using Sparlectra

# measurement guard (RP1 worktree incident 2026-09-03): this script must run
# against THE repository checkout it lives in, never a stale worktree or
# another depot copy picked up through a wrong --project
let expected = normpath(joinpath(@__DIR__, "..", "src")), actual = normpath(String(pathof(Sparlectra)))
  startswith(actual, expected) || error("tools guard: Sparlectra loaded from " * actual * ", expected under " * expected * "; start julia with --project=" * normpath(joinpath(@__DIR__, "..")))
end
using Sockets
using Logging

# --- console versus log ------------------------------------------------------
# The workload is a compile TRACE, and it used to stream every line of it to
# the console: the whole fast test profile, the Web UI group, every service
# run and every warning any of them produced. On Windows that buried the one
# line that mattered (maintainer, 2026-09-07).
#
# From here the console gets a progress line per phase and nothing else; the
# detail goes to sysimage_workload.log next to the image, where a failure can
# be read afterwards. `_CONSOLE` is captured BEFORE any redirection, because
# `stdout` inside a redirected block is the log file.
const _CONSOLE = stdout
const _WORKLOAD_LOG_PATH = try
  joinpath(dirname(String(Base.invokelatest(Sparlectra.webui_sysimage_path))), "sysimage_workload.log")
catch
  joinpath(tempdir(), "sparlectra_sysimage_workload.log")
end
const _WORKLOAD_LOG = try
  mkpath(dirname(_WORKLOAD_LOG_PATH))
  open(_WORKLOAD_LOG_PATH, "w")
catch
  devnull
end

_console(msg) = (println(_CONSOLE, msg); flush(_CONSOLE))

## Run one phase of the trace with its output captured. NOTHING here may
## abort the build: the workload is a trace, and a gap costs first-click
## latency, not correctness. That was the intent all along, but the guard
## was missing at two levels and a Windows-only failure (a temp directory
## that cannot be removed while a file in it is still open) took the whole
## build down with a stacktrace.
function _phase(label::AbstractString, f::Function)
  print(_CONSOLE, "  ", rpad(label, 40), " ")
  flush(_CONSOLE)
  t0 = time()
  outcome = "ok"
  try
    redirect_stdio(stdout = _WORKLOAD_LOG, stderr = _WORKLOAD_LOG) do
      with_logger(SimpleLogger(_WORKLOAD_LOG)) do
        f()
      end
    end
  catch err
    outcome = "FAILED"
    println(_WORKLOAD_LOG, "\n=== phase ", label, " failed ===")
    showerror(_WORKLOAD_LOG, err, catch_backtrace())
    println(_WORKLOAD_LOG)
  end
  flush(_WORKLOAD_LOG)
  println(_CONSOLE, rpad(outcome, 8), round(time() - t0; digits = 1), " s")
  flush(_CONSOLE)
  return outcome == "ok"
end

# AnalyticLoadFlow.jl is a required dependency, so `using Sparlectra` already
# brought it in; the explicit import keeps the APSLF paths in the traced world
# age of the image.
using AnalyticLoadFlow

# nonstandard ports so a Web UI the maintainer has open on 8080 does not
# make the build silently skip the server paths; the first free candidate
# wins (a REPL session may hold an earlier workload port)
const _WORKLOAD_PORTS = 8091:8097

# One service run into a throwaway output root; config_file and output_root
# are filled in here so call sites only state what differs. Failures are
# logged, not fatal: a workload gap costs first-click latency, not the build.
#
# Every run prints one line, successful ones included. That is deliberate:
# a silent success is indistinguishable from a trace that never ran, and
# exactly that cost an hour of hunting a "missing" DTF trace that had in
# fact executed (2026-09-05). The build log now names what was traced.
function _workload_service_run(request::Dict{String,Any}; label::String)
  t0 = time()
  # The temp directory is created and removed BY HAND, not through
  # `mktempdir() do ... end`: on Windows a directory cannot be removed while
  # a file inside it is still open, and a service run leaves log and artifact
  # handles behind. The do-block form raises on cleanup, and that exception
  # used to escape the whole workload and abort the build with a stacktrace
  # (Windows only, maintainer 2026-09-07). A leftover temp directory costs
  # nothing; a failed build costs the user the sysimage.
  outdir = mktempdir(; cleanup = false)
  try
    request["config_file"] = Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH
    request["output_root"] = outdir
    result = start_powerflow_run(request)
    status = String(get(result, "status", ""))
    status == "failed" && @warn "sysimage workload: $(label) service run failed" message = get(result, "message", "?")
    run_id = get(result, "run_id", nothing)
    run_id === nothing || get_powerflow_result(String(run_id))
    println("sysimage workload: ", rpad(label, 28), rpad(isempty(status) ? "?" : status, 12), round(time() - t0; digits = 1), " s")
  catch err
    @warn "sysimage workload: $(label) trace failed (build continues)" exception = (err, catch_backtrace())
  finally
    try
      rm(outdir; recursive = true, force = true)
    catch err
      @warn "sysimage workload: temporary output root could not be removed (harmless)" dir = outdir exception = err
    end
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

# The MATPOWER traces run on the TRACKED case, not on a downloaded one: a
# sysimage build must not depend on the network, and case14.m is not
# shipped with the package. warmup_casePST.m brings its own measurement set
# along, which the state-estimation trace uses.
_workload_matpower_case() = joinpath(pkgdir(Sparlectra), "data", "mpower", "warmup_casePST.m")

function run_workload()
  server = nothing
  port = 0
  for candidate in _WORKLOAD_PORTS
    # ANY failure means "this port is not usable", never "abort the build".
    # The old form rethrew unless the message said "already in use", and the
    # wording differs per platform: on Windows the same condition arrives as
    # an IOError with a different text, so the build died on a busy port.
    server = try
      start_sparlectra_webui(open_browser = false, port = candidate)
    catch err
      @warn "sysimage workload: port $(candidate) not usable" exception = err
      nothing
    end
    if server !== nothing
      port = candidate
      break
    end
  end
  if server === nothing
    @warn "sysimage workload: no free port in $(_WORKLOAD_PORTS); the Web UI paths are not traced (the build continues, the first page view pays the compilation)"
    return nothing
  end
  try
    # request the pages through the real socket path so the whole handler
    # chain is part of the compile trace, not just the render functions.
    # Without this the first form view paid it as JIT time (measured ~11 s).
    _workload_request(port, "/")
    _workload_request(port, "/powerflow")
    # the N-1 weights editor seeds a table from the case's element names (builds
    # the net); trace it so the first "edit N-1 weights" open is not JIT. The case
    # must live in the server's case directory for the seeded-name path.
    try
      cp(_workload_matpower_case(), joinpath(server.runtime.case_directory, "warmup_casePST.m"); force = true)
    catch err
      @warn "sysimage workload: could not stage the MATPOWER case for the weights editor" exception = err
    end
    _workload_request(port, "/powerflow/contingency-weights?case=warmup_casePST.m")
    # real service runs through the full pipeline (case resolution,
    # artifact writers for the CSVs, result.json, run index, result
    # lookup): the first real Web UI run otherwise pays exactly this
    # compilation. Temp output roots keep the user's run history clean.
    # 1. MATPOWER power flow (the tracked warm-up PST case).
    _workload_service_run(Dict{String,Any}(
      "casefile" => _workload_matpower_case(),
      "config_overrides" => Dict{String,Any}("output.logfile_results" => "off", "benchmark.enabled" => false),
    ); label = "matpower power-flow")
    # N-1 contingency (#331 Phase 5): the "Contingency (N-1)" button path through
    # runContingencies!, the CSV/report writers, and the result registry, so the
    # first click after a start on the image is not paid as JIT. Branch kind on the case
    # covers the shared code; the generator kind reuses the same service.
    _workload_service_run(Dict{String,Any}(
      "casefile" => _workload_matpower_case(),
      "contingency_mode" => true,
      "contingency_kind" => "branch",
    ); label = "matpower contingency n-1")
    # State estimation: since stage 4A block 4 the SE section renders on
    # the Runs page (/stateestimation only redirects there, which this
    # request still exercises), plus one SE service run (measurement CSV
    # read, observability, WLS solve, bad-data diagnostics, se_state.csv
    # chain anchor), so the first SE click on the image does not pay
    # the estimator JIT. A failure only costs the trace, never the build.
    _workload_request(port, "/stateestimation")
    _workload_request(port, "/powerflow")
    try
      se_case = _workload_matpower_case()
      if isfile(se_case)
        cfgw = Sparlectra.load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH; reload = true)
        wnet = Sparlectra._import_sparlectra_net(se_case, nothing, cfgw)
        runpf!(wnet, 40, 1e-8, 0; method = :rectangular)
        append!(wnet.measurements, generateMeasurementsFromPF(wnet; noise = false))
        se_meas = joinpath(server.runtime.case_directory, "warmup_casePST.measurements.csv")
        writeMeasurementsCSV(wnet; file = se_meas)
        _workload_service_run(Dict{String,Any}("casefile" => se_case, "se_mode" => true, "measurement_file" => se_meas); label = "state estimation")
      end
    catch err
      @warn "sysimage workload: SE trace skipped" exception = err
    end
    # losses and the human-facing print/report paths: calcNetLosses!, the
    # report builder, the result printer, and the SE diagnostics printer
    # (into devnull) — otherwise the first "show results" click after a
    # start on the image pays their JIT
    try
      se_case2 = _workload_matpower_case()
      if isfile(se_case2)
        cfg2 = Sparlectra.load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH; reload = true)
        pnet = Sparlectra._import_sparlectra_net(se_case2, nothing, cfg2)
        ite2, erg2 = runpf!(pnet, 40, 1e-8, 0; method = :rectangular)
        if erg2 == 0
          calcNetLosses!(pnet)
          buildACPFlowReport(pnet; ct = 0.0, ite = ite2, converged = true)
          redirect_stdout(devnull) do
            printACPFlowResults(pnet, 0.0, ite2, 1e-8)
          end
          append!(pnet.measurements, generateMeasurementsFromPF(pnet; noise = false))
          diag2 = runse_diagnostics(pnet; max_eliminations = 0)
          print_se_diagnostics(devnull, diag2.diagnostics)
        end
      end
    catch err
      @warn "sysimage workload: losses/print trace skipped" exception = err
    end
    # Every IMPORT FORMAT gets a power flow and a state estimation here, not
    # only MATPOWER: the importers and the estimator specialize per format, so
    # a format missing from this trace pays its JIT on the user's first click
    # (maintainer 2026-09-05: "vor dem Sysimage muessen matpower, cgmes und
    # dtf PF + SE + sp_cases gemacht werden"). Each block fails softly - a
    # missing fixture costs the trace, never the build.
    try
      scf_case = joinpath(pkgdir(Sparlectra), "data", "scf", "sp_case14.scf.json")
      if isfile(scf_case)
        snet = importSCF(scf_case)
        runpf!(snet, 40, 1e-8, 0; method = :rectangular)
        # the SCF case carries its own measurement set, which is the path a
        # user takes when running SE straight from a case file
        if !isempty(snet.measurements)
          runse!(snet, Vector{Sparlectra.Measurement}(snet.measurements), Sparlectra.StateEstimationConfig(max_iter = 8))
        end
        exportSCF(snet; file = joinpath(server.runtime.case_directory, "workload_roundtrip.scf.json"))
        _workload_service_run(Dict{String,Any}("casefile" => scf_case); label = "SCF power flow")
        _workload_service_run(Dict{String,Any}("casefile" => scf_case, "se_mode" => true, "measurement_file" => "case"); label = "SCF state estimation")
      end
    catch err
      @warn "sysimage workload: SCF trace skipped" exception = err
    end
    try
      dtf_case = joinpath(pkgdir(Sparlectra), "data", "DTF", "FOR001.DAT")
      if isfile(dtf_case)
        # The DTF parser is a large body of code and one service run does not
        # reach all of it: measured after a first attempt, the first DTF run
        # in a fresh image still cost 2.25 s against 0.007 s for the second.
        # So the native path is exercised directly as well - reader, net
        # builder and case summary - and on a second file, because the
        # branch/outage shapes differ per file.
        for f in ("FOR001.DAT", "FOR001B.DAT")
          fp = joinpath(pkgdir(Sparlectra), "data", "DTF", f)
          isfile(fp) || continue
          try
            dcase = Sparlectra.DTFImporter.read_dtf(fp)
            Sparlectra.DTFImporter.case_summary(dcase)
            dnet = Sparlectra.DTFImporter.build_net(dcase)
            runpf!(dnet, 40, 1e-8, 0; method = :rectangular)
            calcNetLosses!(dnet)
          catch err
            @warn "sysimage workload: DTF direct trace skipped for $(f)" exception = err
          end
        end
        # .DAT alone is ambiguous (FOR002 reference files share it), so the
        # native path has to be named explicitly - the same way the Web UI
        # form does it
        _workload_service_run(Dict{String,Any}("casefile" => dtf_case, "case_format" => "dtf_for001"); label = "DTF power flow")
        # PF and SE per format is the requirement, so DTF gets its estimator
        # run too: the measurements come from the solved net, the service
        # then reads them back through the CSV path like a real request.
        dse = Sparlectra.DTFImporter.build_net(Sparlectra.DTFImporter.read_dtf(dtf_case))
        runpf!(dse, 40, 1e-8, 0; method = :rectangular)
        append!(dse.measurements, generateMeasurementsFromPF(dse; noise = false))
        dmeas = joinpath(server.runtime.case_directory, "dtf_for001.measurements.csv")
        writeMeasurementsCSV(dse; file = dmeas)
        _workload_service_run(Dict{String,Any}("casefile" => dtf_case, "case_format" => "dtf_for001", "se_mode" => true, "measurement_file" => dmeas); label = "DTF state estimation")
      end
    catch err
      @warn "sysimage workload: DTF trace skipped" exception = err
    end
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
      # 4. CGMES state estimation: the estimator on a CGMES-built net, whose
      # tap/machine controls and SV start state differ from the MATPOWER
      # path, so its specializations are separate. Measurements are generated
      # from the solved state, the same route the Web UI generator takes.
      try
        cnet = Sparlectra._se_import_case_net(sc_zip, Sparlectra.load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH; reload = true))
        runpf!(cnet, 40, 1e-8, 0; method = :rectangular)
        append!(cnet.measurements, generateMeasurementsFromPF(cnet; noise = false))
        cmeas = joinpath(server.runtime.case_directory, "cgmes_minigrid.measurements.csv")
        writeMeasurementsCSV(cnet; file = cmeas)
        _workload_service_run(Dict{String,Any}("casefile" => sc_zip, "se_mode" => true, "measurement_file" => cmeas); label = "cgmes state estimation")
      catch err
        @warn "sysimage workload: CGMES SE trace skipped" exception = err
      end
    end
    # one full interactive solve (the case is pre-fetched by the build
    # script); logfile output stays off so the workload never leaves
    # run_case*.log files in examples/_out, regardless of any user
    # configuration
    workload_cfg = SparlectraConfig(output = OutputConfig(logfile_results = :off))
    result = run_sparlectra(casefile = _workload_matpower_case(), config = workload_cfg)
    result.final_converged || @warn "sysimage workload: MATPOWER run did not converge" outcome = result.outcome
  finally
    close(server)
  end
  return nothing
end

_console("[workload] tracing the interactive paths; detail in " * _WORKLOAD_LOG_PATH)

_phase("service and Web UI paths", () -> Base.invokelatest(run_workload))

# Test-suite trace: run the fast-profile test cases so everything they touch
# is compiled into the image as well. Since the risk-based profile resort the
# fast profile carries numerics and model groups only; the full sysimage is
# built FOR the Web UI, so the webui group is traced explicitly afterwards.
# A test failure must not abort the build: the suite is a TRACE here, not a
# gate (the gates run in CI and the developer workflow), and a failing case
# has still compiled everything it touched on the way to failing.
_phase("test profile (fast)", function ()
  ENV["SPARLECTRA_TEST_PROFILE"] = "fast"
  Base.invokelatest(include, joinpath(pkgdir(Sparlectra), "test", "runtests.jl"))
end)

_phase("Web UI test group", function ()
  Base.invokelatest(include, joinpath(pkgdir(Sparlectra), "test", "test_webui.jl"))
  runner = Base.invokelatest(getfield, Main, :run_webui_fast_tests)
  Base.invokelatest(runner)
end)

_console("[workload] done; a FAILED phase above is a trace gap, not a broken image")
try
  close(_WORKLOAD_LOG)
catch
end
