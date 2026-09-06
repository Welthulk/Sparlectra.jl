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

# file: tools/build_app.jl
# purpose: build a relocatable Sparlectra application (a real executable
#          with an embedded Julia runtime, no Julia installation needed on
#          the target machine) with PackageCompiler.create_app in the
#          shared @sparlectra-sysimage-build environment. Flavors:
#          --flavor=full     bin/sparlectra with webui/run/se/n1 commands
#          --flavor=runtime  library runtime only: run/se/n1, no GUI
#          --script=<path>   embed a user script as the app entry point
#          Library users call this through Sparlectra.buildApp(); the app
#          is built locally and is NOT a distribution artifact.
#          No cross-compilation: build on the OS you target.

const _BUILD_ENV = "sparlectra-sysimage-build"
const _DRY_RUN = "--dry-run" in ARGS
const _REPO_ROOT = normpath(joinpath(@__DIR__, ".."))

function _arg(name::String, default::String)::String
  prefix = string("--", name, "=")
  for a in ARGS
    startswith(a, prefix) && return String(a[(length(prefix)+1):end])
  end
  return default
end

function _log(msg::AbstractString)
  println("[build_app] ", msg)
  flush(stdout)
end

const _FLAVOR = Symbol(_arg("flavor", "full"))
const _SCRIPT = _arg("script", "")
const _TARGET = _arg("target", "")

# The generated app package: a tiny wrapper whose julia_main dispatches the
# CLI. The user script variant embeds the script SOURCE as a string, so the
# built app needs no loose files next to the binary.
function _write_app_package(appdir::String; flavor::Symbol, script::String)
  mkpath(joinpath(appdir, "src"))
  script_block = ""
  if !isempty(script)
    content = read(script, String)
    script_block = string("const USER_SCRIPT = ", repr(content), "\nconst USER_SCRIPT_NAME = ", repr(basename(script)), "\n")
  end
  # build timestamp baked into the executable: the Web UI header shows it
  # ("standalone app, built ..."), stamped into ENV at every start
  built_block = string("const APP_BUILT_AT = ", repr(Libc.strftime("%Y-%m-%dT%H:%M:%S", time())), "\n")

  # user-facing CLI help; the webui line exists only in the full flavor
  usage = string(
    "usage: sparlectra <command> [options]\n",
    "\n",
    "commands:\n",
    flavor == :full ? "  webui [--port=N] [--no-browser]  start the local Web UI (also the default without arguments)\n" : "",
    "  run <casefile>              AC power flow (Newton-Raphson)\n",
    "  se <casefile>               state estimation with bad-data diagnostics\n",
    "  n1 <casefile>               N-1 branch/transformer contingency screening\n",
    "  export <casefile>           write the case as a Sparlectra Case Format file\n",
    "  version                     print the Sparlectra version\n",
    "  help                        print this help\n",
    "\n",
    "options for run / se / n1:\n",
    "  --config=<file.yaml>        Sparlectra YAML configuration (library format;\n",
    "                              keys the file does not set keep their defaults)\n",
    "  --set key=value             override one configuration key (repeatable),\n",
    "                              e.g. --set power_flow.max_iter=50\n",
    "  --results[=compact|classic|full]  print the result tables to the console\n",
    "\n",
    "options for export:\n",
    "  --pgm                       write the PLAIN power-grid-model dataset\n",
    "                              (no sparlectra block: no names, slack roles,\n",
    "                              tap nameplates, measurements or configuration)\n",
    "  --out=<file.json>           output path (default: next to the case file)\n",
    "\n",
    "options for se:\n",
    "  --measurements=<file.csv>   measurement set to estimate from; default is\n",
    "                              <case>.measurements.csv next to the case file,\n",
    "                              otherwise synthetic noise-free measurements are\n",
    "                              generated from the power flow\n",
    "\n",
    "options for n1:\n",
    "  --scenario-source=<src>     file_block (the case file's own scenarios\n",
    "                              block), external_file, n1_all, n1_branches\n",
    "                              (default), n1_generators\n",
    "  --scenario-file=<file.json> scenario JSON for --scenario-source=external_file\n",
    "  --screening-mode=<mode>     off, flag, or only (default from the\n",
    "                              contingency.screening configuration)\n",
  )

  # the webui command: a real server in the full flavor, a clear refusal in
  # the runtime flavor (library users who never open the GUI)
  webui_branch = flavor == :full ? """
    if cmd == "webui"
      port = 8080
      open_browser = true
      for a in rest
        if startswith(a, "--port=")
          port = parse(Int, a[(length("--port=")+1):end])
        elseif a == "--no-browser"
          open_browser = false
        else
          println(stderr, "unknown webui option: " * a)
          return 1
        end
      end
      server = Sparlectra.start_sparlectra_webui(open_browser = open_browser, port = port)
      wait(server.task)
      return 0
    end
  """ : """
    if cmd == "webui"
      println(stderr, "this runtime build has no Web UI entry; use the full flavor or the library")
      return 1
    end
  """

  # without arguments the full flavor opens the Web UI (double-click
  # friendly); the runtime flavor prints the help instead
  noargs_branch = flavor == :full ? """
    if isempty(args)
      cmd = "webui"
      rest = String[]
    else
      cmd = args[1]
      rest = args[2:end]
    end
  """ : """
    if isempty(args)
      println(stderr, USAGE)
      return 1
    end
    cmd = args[1]
    rest = args[2:end]
  """

  # a deleted working directory (the shell sat in a folder that was removed,
  # e.g. by a rebuild of this very app) would turn every relative path and
  # abspath call into an IOError, in the webui as a 500 on every page;
  # continue from the home directory instead of failing
  cwd_guard = """
      try
        pwd()
      catch
        println(stderr, "working directory no longer exists; continuing from " * homedir())
        cd(homedir())
      end
      # the Web UI header reads this to show "standalone app, built ..."
      ENV["SPARLECTRA_APP_BUILT"] = APP_BUILT_AT
  """

  main_body = isempty(script) ? """
  function julia_main()::Cint
    args = ARGS
    try
  $(cwd_guard)
  $(noargs_branch)
      if cmd in ("help", "--help", "-h")
        println(USAGE)
        return 0
      end
      if cmd == "version"
        v = Base.pkgversion(Sparlectra)
        println("Sparlectra ", v === nothing ? "unknown version" : string(v))
        return 0
      end
  $(webui_branch)
      if cmd == "run" || cmd == "se" || cmd == "n1" || cmd == "export"
        opts = _prepare(rest)
        cmd == "run" && return _cmd_run(opts)
        cmd == "se" && return _cmd_se(opts)
        cmd == "export" && return _cmd_export(opts)
        return _cmd_n1(opts)
      end
      println(stderr, "unknown command: " * cmd)
      println(stderr, USAGE)
      return 1
    catch err
      println(stderr, sprint(showerror, err))
      return 1
    end
  end
  """ : """
  function julia_main()::Cint
    try
  $(cwd_guard)
      # embedded user script: the whole point of this build
      Base.invokelatest(include_string, Main, USER_SCRIPT, USER_SCRIPT_NAME)
      return 0
    catch err
      println(stderr, sprint(showerror, err))
      return 1
    end
  end
  """

  # command implementations and option parsing; only generated for the CLI
  # variant (a user script brings its own logic)
  helpers = isempty(script) ? """
  const USAGE = $(repr(usage))

  # "true"/"42"/"1.5" become typed scalars, everything else stays a string
  # (the config validation reports wrong types with the offending key)
  function _scalar(v::AbstractString)
    lv = lowercase(v)
    lv == "true" && return true
    lv == "false" && return false
    i = tryparse(Int, v)
    i === nothing || return i
    f = tryparse(Float64, v)
    f === nothing || return f
    return String(v)
  end

  # nest a dotted key ("power_flow.tol") into the override dictionary
  function _nest!(d::Dict{String,Any}, key::AbstractString, value)
    parts = split(String(key), '.')
    cur = d
    for p in parts[1:(end-1)]
      nxt = get!(cur, String(p), Dict{String,Any}())
      nxt isa Dict{String,Any} || error("--set path conflicts with another override: " * key)
      cur = nxt
    end
    cur[String(parts[end])] = value
    return d
  end

  function _set_kv!(overrides::Dict{String,Any}, kv::AbstractString)
    eq = findfirst('=', kv)
    eq === nothing && error("--set expects key=value, got: " * kv)
    _nest!(overrides, kv[1:(eq-1)], _scalar(kv[(eq+1):end]))
    return overrides
  end

  # split positionals from options and ACTIVATE the configuration: --config
  # loads a user YAML over the packaged defaults, --set overrides single
  # keys, --results is shorthand for --set output.logfile_results=<mode>.
  # The active config then steers import, power flow, SE, and output alike.
  function _prepare(rest::Vector{String})
    positionals = String[]
    overrides = Dict{String,Any}()
    config_path = ""
    measurements = ""
    results_mode = ""
    strict_pgm = false
    out_path = ""
    scenario_source = ""
    scenario_file = ""
    screening_mode = ""
    i = 1
    while i <= length(rest)
      a = rest[i]
      if a == "--results"
        results_mode = "classic"
      elseif startswith(a, "--results=")
        results_mode = String(a[(length("--results=")+1):end])
      elseif startswith(a, "--config=")
        config_path = String(a[(length("--config=")+1):end])
      elseif startswith(a, "--measurements=")
        measurements = String(a[(length("--measurements=")+1):end])
      elseif a == "--pgm"
        strict_pgm = true
      elseif startswith(a, "--out=")
        out_path = String(a[(length("--out=")+1):end])
      elseif startswith(a, "--scenario-source=")
        scenario_source = String(a[(length("--scenario-source=")+1):end])
      elseif startswith(a, "--scenario-file=")
        scenario_file = String(a[(length("--scenario-file=")+1):end])
      elseif startswith(a, "--screening-mode=")
        screening_mode = String(a[(length("--screening-mode=")+1):end])
      elseif a == "--set"
        i < length(rest) || error("--set expects key=value")
        i += 1
        _set_kv!(overrides, rest[i])
      elseif startswith(a, "--set=")
        _set_kv!(overrides, String(a[(length("--set=")+1):end]))
      elseif startswith(a, "--")
        error("unknown option: " * a)
      else
        push!(positionals, String(a))
      end
      i += 1
    end
    isempty(results_mode) || _nest!(overrides, "output.logfile_results", results_mode)
    if !isempty(config_path)
      isfile(config_path) || error("config file not found: " * config_path)
      Sparlectra.load_sparlectra_config!(abspath(config_path); reload = true, overrides = overrides)
    elseif !isempty(overrides)
      Sparlectra.load_sparlectra_config!(; reload = true, overrides = overrides)
    end
    return (positionals = positionals, measurements = measurements, results_mode = results_mode, strict_pgm = strict_pgm, out_path = out_path, scenario_source = scenario_source, scenario_file = scenario_file, screening_mode = screening_mode)
  end

  function _casefile(positionals::Vector{String})
    isempty(positionals) && error("missing case file argument (see sparlectra help)")
    casefile = positionals[1]
    isfile(casefile) || error("case file not found: " * casefile)
    return casefile
  end

  function _cmd_run(opts)
    casefile = _casefile(opts.positionals)
    result = Sparlectra.run_sparlectra(casefile = casefile)
    println("converged: ", result.final_converged, ", iterations: ", result.iterations)
    return result.final_converged ? 0 : 1
  end

  function _cmd_se(opts)
    casefile = _casefile(opts.positionals)
    # the power flow supplies the imported network and the SE start state;
    # measurement source: --measurements, else the sidecar CSV next to the
    # case, else synthetic noise-free measurements from this power flow
    result = Sparlectra.run_sparlectra(casefile = casefile)
    result.final_converged || println(stderr, "warning: the initial power flow did not converge")
    net = result.net
    mfile = opts.measurements
    if isempty(mfile)
      sidecar = string(first(splitext(casefile)), ".measurements.csv")
      isfile(sidecar) && (mfile = sidecar)
    end
    if isempty(mfile)
      println("no measurement file found; generating synthetic noise-free measurements from the power flow")
      append!(net.measurements, Sparlectra.generateMeasurementsFromPF(net; noise = false))
    else
      println("measurements: " * mfile)
      Sparlectra.readMeasurementsCSV!(net; file = mfile)
    end
    se = Sparlectra.runse!(net)
    println("SE converged: ", se.converged, ", iterations: ", se.iterations, ", J = ", round(se.objectiveJ, digits = 3), ", dof = ", se.dof, ", J in 3-sigma band: ", se.jWithin3Sigma)
    diag = Sparlectra.runse_diagnostics(net)
    Sparlectra.print_se_diagnostics(stdout, diag.diagnostics)
    if !isempty(opts.results_mode)
      # with the default update_net=true the net carries the estimated
      # state, so the standard result tables show the ESTIMATE
      Sparlectra.calcNetLosses!(net)
      Sparlectra.printACPFlowResults(net, 0.0, se.iterations, Sparlectra.state_estimation_config().tol)
    end
    return se.converged ? 0 : 1
  end

  # write the case as a Sparlectra Case Format file; --pgm drops the
  # namespaced block and writes the plain power-grid-model dataset
  function _cmd_export(opts)
    casefile = _casefile(opts.positionals)
    # import only: an export describes the equipment, so there is no reason
    # to solve the case first (and a non-converging case still exports)
    net = Sparlectra._import_sparlectra_net(String(abspath(casefile)), nothing, Sparlectra.load_sparlectra_config())
    sidecar = string(first(splitext(abspath(casefile))), ".measurements.csv")
    isfile(sidecar) && Sparlectra.readMeasurementsCSV!(net; file = sidecar)
    stem = first(splitext(basename(casefile)))
    out = isempty(opts.out_path) ? joinpath(dirname(abspath(casefile)), stem * (opts.strict_pgm ? ".pgm.json" : ".scf.json")) : opts.out_path
    path = Sparlectra.exportSCF(net; file = out, case_name = stem, source_reference = basename(casefile), strict_pgm = opts.strict_pgm)
    println("wrote ", path, " (", round(filesize(path) / 1024; digits = 1), " kB)")
    opts.strict_pgm && println("plain power-grid-model dataset: names, slack roles, tap nameplates, measurements and configuration are NOT in this file")
    return 0
  end

  # N-1 / scenario batch (scenario task step 5): the scenario source picks
  # the case list (default: generated branch outages, the historical CLI
  # behavior); file_block and external_file run the scenario model addressed
  # by SCF component ids, so they need an SCF case file. The screening mode
  # defaults to the contingency.screening configuration.
  function _cmd_n1(opts)
    casefile = _casefile(opts.positionals)
    source = isempty(opts.scenario_source) ? "n1_branches" : opts.scenario_source
    source in ("file_block", "external_file", "n1_all", "n1_branches", "n1_generators") || error("--scenario-source must be file_block, external_file, n1_all, n1_branches, or n1_generators")
    screening = isempty(opts.screening_mode) ? Sparlectra.contingency_config().screening_mode : Symbol(lowercase(opts.screening_mode))
    margin = Sparlectra.contingency_config().screening_margin_pct
    ladder = Sparlectra.contingency_config().rescue_ladder
    result = Sparlectra.run_sparlectra(casefile = casefile)
    if !result.final_converged
      println(stderr, "the base-case power flow did not converge")
      return 1
    end
    net = result.net
    results = if source in ("file_block", "external_file")
      endswith(lowercase(casefile), ".json") || error("--scenario-source=" * source * " resolves component ids against the typed case and needs a Sparlectra Case Format case (.scf.json)")
      scfcase = Sparlectra.read_scf_json(abspath(casefile))
      set = if source == "file_block"
        s = Sparlectra.scf_case_scenarios(scfcase)
        s === nothing && error("the case file carries neither a scenarios nor a contingencies block")
        s
      else
        isfile(opts.scenario_file) || error("--scenario-source=external_file needs --scenario-file=<file.json>")
        Sparlectra.scenario_set_from_dict(Sparlectra.scf_json_parse(read(opts.scenario_file, String)))
      end
      println("running ", length(set.scenarios), " scenario(s) from ", source, " (screening ", screening, ")...")
      Sparlectra.runScenarios!(net, set; index = Sparlectra.ScenarioIndex(scfcase), rescue_ladder = ladder, screening_mode = screening, screening_margin_pct = margin)
    else
      cases = source == "n1_generators" ? Sparlectra.generateN1Generators(net) : source == "n1_all" ? vcat(Sparlectra.generateN1Branches(net), Sparlectra.generateN1Generators(net)) : Sparlectra.generateN1Branches(net)
      println("running ", length(cases), " N-1 case(s) (screening ", screening, ")...")
      Sparlectra.runContingencies!(net, cases; rescue_ladder = ladder, screening_mode = screening, screening_margin_pct = margin)
    end
    screening === :off || println(count(r -> r.screened, results), " of ", length(results), " row(s) screened (first-order estimates, no full solve)")
    report = Sparlectra.buildContingencyReport(results)
    Sparlectra.printContingencyReport(stdout, report)
    return 0
  end
  """ : ""

  open(joinpath(appdir, "src", "SparlectraApp.jl"), "w") do io
    print(
      io,
      """
      module SparlectraApp

      using Sparlectra
      $(script_block)
      $(built_block)
      $(helpers)
      $(main_body)
      end # module
      """,
    )
  end
  open(joinpath(appdir, "Project.toml"), "w") do io
    println(io, "name = \"SparlectraApp\"")
    println(io, "uuid = \"7c1d43f2-64d1-4f38-9f6a-2a52c0ffee01\"")
    println(io, "version = \"0.1.0\"")
    println(io, "")
    println(io, "[deps]")
  end
  return nothing
end

function main()
  flavor = _FLAVOR
  flavor in (:full, :runtime) || error("build_app: --flavor must be full or runtime (got $(flavor))")
  script = _SCRIPT
  !isempty(script) && !isfile(script) && error("build_app: --script file not found: $(script)")
  target = isempty(_TARGET) ? joinpath(homedir(), "sparlectra-app") : _TARGET
  if _DRY_RUN
    _log("dry run: would generate the app package (flavor $(flavor)$(isempty(script) ? "" : ", user script $(basename(script))"))")
    _log("would build into: $(target)")
    _log("executable would be: $(joinpath(target, "bin", Sys.iswindows() ? "sparlectra.exe" : "sparlectra"))")
    _log("dry run finished")
    return nothing
  end
  started = time()
  @eval using Pkg
  pkgm = Base.invokelatest(getfield, @__MODULE__, :Pkg)
  _log("activating shared build environment @$(_BUILD_ENV)")
  Base.invokelatest(pkgm.activate, _BUILD_ENV; shared = true)
  if Base.find_package("PackageCompiler") === nothing
    _log("installing PackageCompiler into the build environment (one-time)")
    Base.invokelatest(pkgm.add, "PackageCompiler")
  end
  @eval using PackageCompiler
  pc = Base.invokelatest(getfield, @__MODULE__, :PackageCompiler)

  appdir = mktempdir(; prefix = "sparlectra_app_")
  _log("generating the app package in $(appdir) (flavor $(flavor))")
  _write_app_package(appdir; flavor = flavor, script = script)
  # resolve the app's dependencies against THIS checkout/installation; the
  # APSLF solver arrives with Sparlectra (required dependency), so the app
  # itself declares nothing extra
  Base.invokelatest(pkgm.activate, appdir)
  Base.invokelatest(pkgm.develop; path = _REPO_ROOT)
  Base.invokelatest(pkgm.instantiate)

  workload = flavor == :runtime ? joinpath(@__DIR__, "app_workload_runtime.jl") : joinpath(@__DIR__, "sysimage_workload.jl")
  _log("building the app (this typically takes 20 to 40 minutes)...")
  _log("target: $(target)")
  Base.invokelatest(getglobal(pc, :create_app), appdir, target; executables = ["sparlectra" => "julia_main"], precompile_execution_file = workload, force = true, include_lazy_artifacts = true)
  exe = joinpath(target, "bin", Sys.iswindows() ? "sparlectra.exe" : "sparlectra")
  elapsed = round((time() - started) / 60; digits = 1)
  _log("done in $(elapsed) min: $(exe)")
  _log("the app folder is relocatable; no Julia installation is needed to run it")
  _log("no cross-compilation: this binary runs on $(Sys.iswindows() ? "Windows" : Sys.isapple() ? "macOS" : "Linux") only")
  return nothing
end

Base.invokelatest(main)
