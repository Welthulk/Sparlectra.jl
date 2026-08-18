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

# file: src/webui/webui.jl
# purpose: local Web UI server (start_sparlectra_webui): TCP request loop,
#          output-root defaults, browser launch, warmup, and lifecycle and
#          shutdown handling
using Markdown
using Sockets

mutable struct _SparlectraWebUIRuntime
  listener::Union{Sockets.TCPServer,Nothing}
  case_directory::String
  config_file::String
  operation_log::String
  startup_config_error::Union{String,Nothing}
  runner
  shutdown_on_browser_close::Bool
  heartbeat_received::Bool
  last_heartbeat::Float64
  active_requests::Int
  shutdown_reason::Union{Symbol,Nothing}
  lifecycle_io::IO
  lock::ReentrantLock
  warmup_state::Symbol
end

"""
    default_webui_output_root() -> String

Return the operating-system-specific, user-writable directory used for Web UI
runs, downloaded MATPOWER cases, and the operation log.
"""
function default_webui_output_root()::String
  if Sys.iswindows()
    root = get(ENV, "LOCALAPPDATA", "")
    return isempty(root) ? joinpath(pwd(), "sparlectra_webui_runs") : joinpath(root, "Sparlectra", "WebUI", "runs")
  elseif Sys.isapple()
    home = homedir()
    return isempty(home) ? joinpath(pwd(), "sparlectra_webui_runs") : joinpath(home, "Library", "Application Support", "Sparlectra", "WebUI", "runs")
  elseif Sys.islinux()
    root = get(ENV, "XDG_STATE_HOME", "")
    isempty(root) && (root = joinpath(homedir(), ".local", "state"))
    return joinpath(root, "sparlectra", "webui", "runs")
  end
  return joinpath(pwd(), "sparlectra_webui_runs")
end

"""Return the user-writable default Web UI configuration path."""
default_webui_config_path(output_root::AbstractString = default_webui_output_root())::String =
  joinpath(dirname(abspath(output_root)), "config", "configuration.yaml")

"""Return the user-writable MATPOWER case cache used by the Web UI."""
default_webui_case_cache_dir(output_root::AbstractString = default_webui_output_root())::String =
  joinpath(dirname(abspath(output_root)), "data", "mpower")

"""Return the Web UI operation-log path beneath the user application root."""
default_webui_operation_log_path(output_root::AbstractString = default_webui_output_root())::String =
  joinpath(dirname(abspath(output_root)), "logs", WEBUI_OPERATION_LOG_FILENAME)

function _sparlectra_package_path()::String
  pkg_path = pathof(@__MODULE__)
  return pkg_path === nothing ? "" : String(pkg_path)
end

function _sparlectra_git_commit_sha()::Union{String,Nothing}
  package_path = _sparlectra_package_path()
  isempty(package_path) && return nothing
  root = abspath(joinpath(dirname(package_path), ".."))
  return _git_head_commit_sha(root)
end

"""
    _git_head_commit_sha(root) -> Union{String,Nothing}

Resolve the HEAD commit SHA of the git checkout at `root` by reading git's
on-disk files (no `git` executable required). Handles plain checkouts
(`.git` directory), `git worktree` checkouts (`.git` is a `gitdir: <path>`
pointer file whose shared refs live in the `commondir`), detached HEADs, and
refs packed into `packed-refs`. Returns `nothing` when `root` is not a git
checkout or the SHA cannot be determined.
"""
function _git_head_commit_sha(root::AbstractString)::Union{String,Nothing}
  dotgit = joinpath(root, ".git")
  gitdir = ""
  if isdir(dotgit)
    gitdir = dotgit
  elseif isfile(dotgit)
    pointer = strip(read(dotgit, String))
    startswith(pointer, "gitdir:") || return nothing
    target = strip(pointer[(sizeof("gitdir:") + 1):end])
    isempty(target) && return nothing
    gitdir = isabspath(target) ? target : abspath(joinpath(root, target))
  else
    return nothing
  end

  head_path = joinpath(gitdir, "HEAD")
  isfile(head_path) || return nothing
  head = strip(read(head_path, String))
  isempty(head) && return nothing
  startswith(head, "ref:") || return head

  ref = strip(chop(head; head = 4, tail = 0))
  isempty(ref) && return nothing

  common = gitdir
  commondir_file = joinpath(gitdir, "commondir")
  if isfile(commondir_file)
    target = strip(read(commondir_file, String))
    isempty(target) || (common = isabspath(target) ? target : abspath(joinpath(gitdir, target)))
  end

  ref_path = joinpath(common, split(ref, '/')...)
  if isfile(ref_path)
    sha = strip(read(ref_path, String))
    return isempty(sha) ? nothing : sha
  end

  packed = joinpath(common, "packed-refs")
  isfile(packed) || return nothing
  # open(...) do closes the handle on every exit path. eachline(path) would
  # keep the file open until GC after the early return below — on Windows the
  # still-open packed-refs then blocks directory removal (EBUSY, issue #290
  # test cleanup) and every banner-SHA resolution leaks one descriptor.
  return open(packed) do io
    for line in eachline(io)
      entry = strip(line)
      (isempty(entry) || startswith(entry, '#') || startswith(entry, '^')) && continue
      parts = split(entry; limit = 2)
      length(parts) == 2 || continue
      if strip(parts[2]) == ref
        sha = strip(parts[1])
        return isempty(sha) ? nothing : sha
      end
    end
    return nothing
  end
end

include("sysimage.jl")
include("options.jl")
include("forms.jl")
include("docs.jl")
include("operations.jl")
include("views.jl")
include("handlers.jl")
include("routes.jl")

mutable struct SparlectraWebUIServer
  listener::Sockets.TCPServer
  task::Task
  url::String
  runtime::_SparlectraWebUIRuntime
  browser_monitor_task::Union{Task,Nothing}
end

function _webui_status_text(status::Int)
  return get(Dict(200 => "OK", 204 => "No Content", 303 => "See Other", 400 => "Bad Request", 404 => "Not Found", 500 => "Internal Server Error"), status, "OK")
end

function _webui_read_request(socket::Sockets.TCPSocket)
  request_line = chomp(readline(socket))
  parts = split(request_line)
  length(parts) >= 2 || throw(ArgumentError("Invalid HTTP request line."))
  method, target = parts[1], parts[2]
  headers = Dict{String,String}()
  while !eof(socket)
    line = strip(readline(socket), ['\r', '\n'])
    isempty(line) && break
    header = split(line, ':'; limit = 2)
    length(header) == 2 && (headers[lowercase(strip(header[1]))] = strip(header[2]))
  end
  length_bytes = parse(Int, get(headers, "content-length", "0"))
  content_type_raw = get(headers, "content-type", "")
  if length_bytes > WEBUI_CASE_IMPORT_MAX_REQUEST_BYTES && occursin("multipart/form-data", lowercase(content_type_raw))
    return method, target, Dict{String,Any}("casefiles" => WebUICaseUpload[], "case_import_request_oversized" => "true")
  end
  body_bytes = length_bytes > 0 ? read(socket, length_bytes) : UInt8[]
  content_type = lowercase(content_type_raw)
  form = if startswith(content_type, "application/x-www-form-urlencoded")
    _webui_parse_pairs(String(body_bytes))
  elseif startswith(content_type, "multipart/form-data")
    _webui_parse_multipart_form(body_bytes, content_type_raw)
  else
    Dict{String,String}()
  end
  return method, target, form
end

function _webui_parse_multipart_form(body::Vector{UInt8}, content_type::AbstractString)
  boundary_match = match(r"boundary=([^;]+)", content_type)
  boundary_match === nothing && return Dict{String,Any}()
  boundary = "--" * strip(boundary_match.captures[1], ['"'])
  text = String(body)
  form = Dict{String,Any}()
  uploads = WebUICaseUpload[]
  for part in split(text, boundary)
    occursin("\r\n\r\n", part) || continue
    header_text, content = split(part, "\r\n\r\n"; limit = 2)
    content = replace(content, r"\r\n--$" => "")
    content = replace(content, r"\r\n$" => "")
    name_match = match(r"name=\"([^\"]+)\"", header_text)
    name_match === nothing && continue
    filename_match = match(r"filename=\"([^\"]*)\"", header_text)
    field = name_match.captures[1]
    if filename_match !== nothing
      push!(uploads, WebUICaseUpload(filename_match.captures[1], Vector{UInt8}(codeunits(content))))
    else
      form[field] = content
    end
  end
  form["casefiles"] = uploads
  return form
end

function _webui_write_response(socket::Sockets.TCPSocket, response::SparlectraWebUIResponse)
  headers = copy(response.headers)
  push!(headers, "Content-Length" => string(length(response.body)))
  push!(headers, "Connection" => "close")
  write(socket, "HTTP/1.1 $(response.status) $(_webui_status_text(response.status))\r\n")
  for (name, value) in headers
    write(socket, name, ": ", value, "\r\n")
  end
  write(socket, "\r\n")
  write(socket, response.body)
  flush(socket)
  return nothing
end

function _webui_serve_client(socket::Sockets.TCPSocket, output_root::String, runtime::_SparlectraWebUIRuntime)
  _webui_begin_request!(runtime)
  try
    method, target, form = _webui_read_request(socket)
    response = route_sparlectra_webui(method, target, form; output_root, runtime)
    _webui_write_response(socket, response)
  catch err
    record_webui_operation!(output_root, "internal_error"; method = "UNKNOWN", route = "request", status = 500, message = sprint(showerror, err), user_action = false)
    try
      _webui_write_response(socket, _webui_html(render_webui_error(500, sprint(showerror, err)); status = 500))
    catch
    end
  finally
    _webui_finish_request!(runtime)
    close(socket)
  end
  return nothing
end

function _webui_platform()::Symbol
  return Sys.iswindows() ? :windows : Sys.isapple() ? :macos : :linux
end

function _webui_first_executable(candidates; executable_lookup = Sys.which, path_exists = isfile)
  for candidate in candidates
    if isabspath(candidate)
      path_exists(candidate) && return candidate
    else
      executable = executable_lookup(candidate)
      executable === nothing || return String(executable)
    end
  end
  return nothing
end

function _webui_windows_browser_candidates(environment)::Vector{String}
  candidates = String[]
  locations = (
    ("PROGRAMFILES", joinpath("Microsoft", "Edge", "Application", "msedge.exe")),
    ("PROGRAMFILES(X86)", joinpath("Microsoft", "Edge", "Application", "msedge.exe")),
    ("LOCALAPPDATA", joinpath("Microsoft", "Edge", "Application", "msedge.exe")),
    ("PROGRAMFILES", joinpath("Google", "Chrome", "Application", "chrome.exe")),
    ("PROGRAMFILES(X86)", joinpath("Google", "Chrome", "Application", "chrome.exe")),
    ("LOCALAPPDATA", joinpath("Google", "Chrome", "Application", "chrome.exe")),
    ("PROGRAMFILES", joinpath("BraveSoftware", "Brave-Browser", "Application", "brave.exe")),
    ("LOCALAPPDATA", joinpath("BraveSoftware", "Brave-Browser", "Application", "brave.exe")),
  )
  for (variable, suffix) in locations
    root = get(environment, variable, "")
    isempty(root) || push!(candidates, joinpath(root, suffix))
  end
  append!(candidates, ("msedge.exe", "chrome.exe", "brave.exe"))
  return candidates
end

const _WEBUI_APP_WINDOW_SIZE = "1500,950"

function _webui_app_window_command(url::String; platform::Symbol = _webui_platform(), executable_lookup = Sys.which, path_exists = ispath, environment = ENV)::Union{Cmd,Nothing}
  if platform == :macos
    applications = (
      "/Applications/Microsoft Edge.app",
      "/Applications/Google Chrome.app",
      "/Applications/Chromium.app",
      "/Applications/Brave Browser.app",
    )
    application = findfirst(path_exists, applications)
    application === nothing && return nothing
    return `open -na $(applications[application]) --args --app=$url --window-size=$(_WEBUI_APP_WINDOW_SIZE)`
  end

  candidates = platform == :windows ? _webui_windows_browser_candidates(environment) : String[
    "microsoft-edge",
    "microsoft-edge-stable",
    "google-chrome",
    "google-chrome-stable",
    "chromium",
    "chromium-browser",
    "brave-browser",
  ]
  executable = _webui_first_executable(candidates; executable_lookup, path_exists)
  executable === nothing && return nothing
  return `$executable --app=$url --window-size=$(_WEBUI_APP_WINDOW_SIZE)`
end

function _webui_generic_open_command(url::String; platform::Symbol = _webui_platform(), executable_lookup = Sys.which)::Union{Tuple{Cmd,Symbol},Nothing}
  platform == :linux || return nothing

  xdg_open = executable_lookup("xdg-open")
  xdg_open === nothing || return (`$(String(xdg_open)) $url`, :xdg_open)

  gio = executable_lookup("gio")
  gio === nothing || return (`$(String(gio)) open $url`, :gio_open)

  sensible_browser = executable_lookup("sensible-browser")
  sensible_browser === nothing || return (`$(String(sensible_browser)) $url`, :sensible_browser)

  return nothing
end

function _webui_browser_open_command(url::String; platform::Symbol = _webui_platform(), executable_lookup = Sys.which, path_exists = ispath, environment = ENV)::Union{Tuple{Cmd,Symbol},Nothing}
  app_command = _webui_app_window_command(url; platform, executable_lookup, path_exists, environment)
  app_command === nothing || return (app_command, :app_window)
  return _webui_generic_open_command(url; platform, executable_lookup)
end

function _webui_app_command(url::String; platform::Symbol = _webui_platform(), executable_lookup = Sys.which, path_exists = ispath, environment = ENV)::Union{Cmd,Nothing}
  selected = _webui_browser_open_command(url; platform, executable_lookup, path_exists, environment)
  selected === nothing && return nothing
  return selected[1]
end

# Probe whether the process holding host:port is a Sparlectra Web UI: fetch
# /powerflow and look for the page marker. Bounded by a timeout so a silent
# non-HTTP listener cannot hang the startup path; any failure means "not ours".
function _webui_port_holds_sparlectra(host_string::AbstractString, port::Integer; timeout_seconds::Real = 3.0)::Bool
  address = host_string == "localhost" ? "127.0.0.1" : String(host_string)
  sock = try
    Sockets.connect(address, UInt16(port))
  catch
    return false
  end
  try
    write(sock, "GET /powerflow HTTP/1.1\r\nHost: $(address):$(port)\r\nConnection: close\r\n\r\n")
    reader = @async try
      read(sock)
    catch
      UInt8[]
    end
    if timedwait(() -> istaskdone(reader), Float64(timeout_seconds)) === :ok
      return occursin("Sparlectra", String(fetch(reader)))
    end
    return false
  catch
    return false
  finally
    try
      close(sock)
    catch
    end
  end
end

function _webui_open_browser(url::String)
  selected = _webui_browser_open_command(url)
  if selected === nothing
    @info "Selected Web UI browser-opening strategy" strategy = "manual_only" url
    @warn "Could not find an app-capable browser; open the Web UI URL manually" url
    return nothing
  end
  command, strategy = selected
  @info "Selected Web UI browser-opening strategy" strategy = String(strategy) url
  try
    return run(command; wait = false)
  catch err
    @warn "Could not open the Web UI app window" url exception = (err, catch_backtrace())
  end
  return nothing
end

function _webui_set_warmup_state!(runtime::_SparlectraWebUIRuntime, state::Symbol)
  lock(runtime.lock) do
    runtime.warmup_state = state
  end
  return nothing
end

function _webui_warmup_in_progress(runtime::_SparlectraWebUIRuntime)::Bool
  return lock(runtime.lock) do
    # :waiting_first_page → server up, warm-up solve deliberately deferred
    # until the warm-up page has been served once; :page_served → page went
    # out, solve about to start; :warming → solve running.
    runtime.warmup_state in (:waiting_first_page, :page_served, :warming)
  end
end

"""Flip the warm-up gate once the warm-up page has been rendered for the
browser, so the CPU-bound warm-up solve only starts after the user can see
the "warming up" message instead of a blank tab."""
function _webui_mark_warmup_page_served!(runtime::_SparlectraWebUIRuntime)
  lock(runtime.lock) do
    runtime.warmup_state === :waiting_first_page && (runtime.warmup_state = :page_served)
  end
  return nothing
end

_webui_mark_warmup_page_served!(::Any) = nothing

function _webui_warmup_page_served(runtime::_SparlectraWebUIRuntime)::Bool
  return lock(runtime.lock) do
    runtime.warmup_state !== :waiting_first_page
  end
end

# Some tests pass a lightweight NamedTuple stand-in for `runtime` (only the
# fields a particular handler needs) instead of a full _SparlectraWebUIRuntime;
# such stand-ins never warm up.
_webui_warmup_in_progress(::Any)::Bool = false

function _webui_record_heartbeat!(runtime::_SparlectraWebUIRuntime)
  lock(runtime.lock) do
    runtime.heartbeat_received = true
    runtime.last_heartbeat = time()
  end
  return nothing
end

function _webui_begin_request!(runtime::_SparlectraWebUIRuntime)
  lock(runtime.lock) do
    runtime.active_requests += 1
  end
  return nothing
end

function _webui_finish_request!(runtime::_SparlectraWebUIRuntime)
  lock(runtime.lock) do
    runtime.active_requests -= 1
    runtime.heartbeat_received && (runtime.last_heartbeat = time())
  end
  return nothing
end

function _webui_shutdown_reason_text(reason::Union{Symbol,Nothing})::String
  reason === :explicit_shutdown && return " by explicit shutdown"
  reason === :browser_window_close && return " by browser window close"
  reason === :ctrl_c && return " by Ctrl+C"
  reason === :server_closed && return " by server close"
  reason === :process_exit && return " by process exit"
  reason === :hard_reset && return " after hard reset"
  return ""
end

function _webui_lifecycle_println(runtime::_SparlectraWebUIRuntime, parts...)
  println(runtime.lifecycle_io, parts...)
  flush(runtime.lifecycle_io)
  return nothing
end

function _webui_startup_log(io::IO, event::AbstractString; operation_log::Union{Nothing,AbstractString} = nothing, fields...)
  _ = io
  operation_log === nothing || record_webui_operation!(operation_log, event; route = "/powerflow", method = "START", user_action = false, fields...)
  return nothing
end

function _webui_startup_failure!(io::IO, err, bt; operation_log::Union{Nothing,AbstractString} = nothing, phase::AbstractString = "startup")
  message = sprint(showerror, err)
  _webui_startup_log(io, "webui_start_failed"; operation_log, status = "failed", phase, message)
  showerror(io, err, bt)
  println(io)
  flush(io)
  return nothing
end

function _webui_request_shutdown!(runtime::_SparlectraWebUIRuntime; reason::Union{Symbol,Nothing} = nothing)
  listener = lock(runtime.lock) do
    if runtime.shutdown_reason === nothing && reason !== nothing
      runtime.shutdown_reason = reason
    end
    runtime.listener
  end
  listener === nothing || (isopen(listener) && close(listener))
  return nothing
end

function _webui_monitor_browser_process(runtime::_SparlectraWebUIRuntime, browser_process)
  try
    wait(browser_process)
  catch err
    err isa InvalidStateException || @debug "Web UI browser process monitor ended without shutdown" exception = (err, catch_backtrace())
    return nothing
  end
  should_shutdown = lock(runtime.lock) do
    runtime.shutdown_on_browser_close && runtime.listener !== nothing && isopen(runtime.listener)
  end
  should_shutdown && _webui_request_shutdown!(runtime; reason = :browser_window_close)
  return nothing
end

"""
    start_sparlectra_webui(; host="127.0.0.1", port=8080,
                            output_root=nothing,
                            config_file=DEFAULT_SPARLECTRA_CONFIG_PATH,
                            open_browser=false,
                            shutdown_on_browser_close=false,
                            warmup=false, warmup_casefile=nothing,
                            warmup_store_result=false,
                            warmup_delay_seconds=2.0) -> Union{Nothing,SparlectraWebUIServer}

Start the loopback-only PowerFlow interface and load its persistent run registry
before accepting requests. The returned handle can be stopped with
`close(server)` or the browser's **Stop Web UI** button. When the requested
port is already held by another **Sparlectra** Web UI (probed via its
`/powerflow` page), no error is raised: the running instance is opened in
the browser (with `open_browser = true`) and `nothing` is returned — stop
that instance first to actually restart. A foreign process on the port
still raises the explicit `ArgumentError`. Browser-process
lifetime is not used for automatic shutdown by default because common browsers
may return a short-lived launcher process instead of a reliably owned window.
`warmup` defaults to `nothing`: the configuration decides (`webui.warmup`,
default `true`, editable under the form's Advanced options); an explicit
`true`/`false` always wins. When warm-up is active, hidden asynchronous runs
compile the common import/API/solver path plus the short-circuit path.
When `warmup=true`, hidden asynchronous runs compile the common import/API/
solver path. By default the bundled synthetic 3-bus case runs first (plumbing),
followed by the bundled synthetic 118-bus case (realistic sparse-solve and
reporting paths), each in a temporary output directory, so no run-history
entries are added and no artifacts retained. An explicit `warmup_casefile`
replaces the whole sequence. Set `warmup_store_result=true` to retain warm-up
artifacts beneath the configured output root (one subdirectory per case).

The warm-up run is scheduled on a background task, but Julia's single-threaded
cooperative scheduler still means a CPU-bound solve of that task can delay the
HTTP server from responding to the *first* browser request until the solve
reaches a yield point. The warm-up solve therefore only starts after the
"warming up" page has actually been served once (with a 30 s fallback when no
browser connects), followed by a `warmup_delay_seconds` grace period (default
`2.0`) so the browser can paint it. After the solve, the PowerFlow form render
path is pre-compiled as well, so the automatic refresh that replaces the
warm-up page is fast.
"""
function _run_sparlectra_webui_warmup(output_root::AbstractString; warmup_casefile::Union{Nothing,AbstractString} = nothing, warmup_store_result::Bool = false, runner = run_sparlectra_api)
  # Default warm-up sequence: the 3-bus case warms the plumbing (parser, config,
  # result pipeline), the synthetic 118-bus case afterwards warms the realistic
  # sparse-solve and reporting paths, so the user's first real run does not pay
  # that compilation. An explicit warmup_casefile replaces the whole sequence.
  casefiles = if warmup_casefile === nothing
    [joinpath(_WEBUI_PACKAGE_ROOT, "data", "webui", "warmup_case3.jl"), joinpath(_WEBUI_PACKAGE_ROOT, "data", "webui", "warmup_case118.jl")]
  else
    [abspath(warmup_casefile)]
  end
  for f in casefiles
    isfile(f) || throw(ArgumentError("Web UI warm-up case file not found: $(f)"))
  end
  execute(casefile, output_dir) = runner(
    casefile = casefile,
    config_file = DEFAULT_SPARLECTRA_CONFIG_PATH,
    output_dir = output_dir,
    config_overrides = Dict("output.logfile_results" => "off", "benchmark.enabled" => false),
    performance_timing = :off,
    run_diagnostics = false,
  )
  # Each case gets its own output directory (the runner writes fixed-name
  # artifacts). The last result is returned; a failed run aborts the sequence
  # and returns its result so the caller sees the failure.
  result = nothing
  for casefile in casefiles
    result = if warmup_store_result
      output_dir = joinpath(abspath(output_root), "webui-warmup", first(splitext(basename(casefile))))
      mkpath(output_dir)
      execute(casefile, output_dir)
    else
      mktempdir() do output_dir
        execute(casefile, output_dir)
      end
    end
    (result !== nothing && hasproperty(result, :success) && !result.success) && return result
  end
  _warm_up_short_circuit_path()
  return result
end

# The Short-circuit button hits a code path the power-flow warm-up never
# touches (Z-bus column solve, IEC c-factor table, result/coverage types),
# so its first click paid the full compilation. Warm it on a throwaway
# two-bus net with one feeder — no CGMES delivery needed, `runShortCircuit!`
# accepts any `Net` plus harvested data. Failures stay silent: warm-up must
# never affect startup.
function _warm_up_short_circuit_path()
  try
    net = Net(name = "webui_warmup_sc", baseMVA = 100.0)
    addBus!(net = net, busName = "WU1", vn_kV = 110.0)
    addBus!(net = net, busName = "WU2", vn_kV = 110.0)
    addPIModelACLine!(net = net, fromBus = "WU1", toBus = "WU2", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
    addProsumer!(net = net, busName = "WU1", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "WU1")
    addProsumer!(net = net, busName = "WU2", type = "ENERGYCONSUMER", p = 10.0, q = 2.0)
    feeder = (; mrid = "warmup_feeder", name = "WU_FEEDER", bus = "WU1",
      maxInitialSymShCCurrent_A = 20000.0, minInitialSymShCCurrent_A = 15000.0,
      maxR1ToX1Ratio = 0.1, minR1ToX1Ratio = 0.1)
    sc_data = CGMESImporter.CGMESShortCircuitData([feeder], NamedTuple[], NamedTuple[], NamedTuple[], NamedTuple[], NamedTuple[])
    for case in (:max, :min)
      runShortCircuit!(net, sc_data; case = case)
    end
  catch err
    @debug "Web UI short-circuit warm-up skipped" exception = (err, catch_backtrace())
  end
  return nothing
end

function _provision_webui_runtime!(root::AbstractString, config_file::Union{Nothing,AbstractString})
  configuration = config_file === nothing ? default_webui_config_path(root) : abspath(config_file)
  case_directory = default_webui_case_cache_dir(root)
  operation_log = default_webui_operation_log_path(root)
  mkpath.(unique((abspath(root), dirname(configuration), case_directory, dirname(operation_log))))
  config_file === nothing && !isfile(configuration) && cp(DEFAULT_SPARLECTRA_CONFIG_PATH, configuration)
  isfile(configuration) || throw(ArgumentError("Web UI configuration file not found: $(configuration)"))
  for warmup_name in ("warmup_case3.jl", "warmup_case118.jl")
    source = joinpath(_WEBUI_PACKAGE_ROOT, "data", "webui", warmup_name)
    isfile(source) || continue
    destination = joinpath(case_directory, basename(source))
    isfile(destination) || cp(source, destination)
  end
  return (config_file = configuration, case_directory, operation_log)
end

function _webui_validate_startup_config(configuration::AbstractString)
  try
    load_sparlectra_config(configuration; reload = true)
    return nothing
  catch err
    return "Configuration error in $(configuration): $(sprint(showerror, err))"
  end
end

function start_sparlectra_webui(; host::AbstractString = "127.0.0.1", port::Integer = 8080, output_root::Union{Nothing,AbstractString} = nothing, config_file::Union{Nothing,AbstractString} = nothing, open_browser::Bool = false, shutdown_on_browser_close::Bool = false, auto_shutdown_on_browser_close::Union{Nothing,Bool} = nothing, browser_heartbeat_timeout_seconds::Real = 15.0, warmup::Union{Nothing,Bool} = nothing, warmup_casefile::Union{Nothing,AbstractString} = nothing, warmup_store_result::Bool = false, warmup_delay_seconds::Real = 2.0, _test_runner = start_powerflow_run, _lifecycle_io::IO = stdout, _browser_opener = _webui_open_browser)::Union{Nothing,SparlectraWebUIServer}
  host_string = String(host)
  host_string in ("127.0.0.1", "localhost", "::1") || (err = ArgumentError("Sparlectra Web UI only accepts loopback hosts: 127.0.0.1, localhost, or ::1."); _webui_startup_failure!(_lifecycle_io, err, catch_backtrace(); phase = "validate_arguments"); throw(err))
  1 <= port <= 65535 || (err = ArgumentError("Web UI port must be between 1 and 65535."); _webui_startup_failure!(_lifecycle_io, err, catch_backtrace(); phase = "validate_arguments"); throw(err))
  timeout = Float64(browser_heartbeat_timeout_seconds)
  isfinite(timeout) && timeout > 0 || (err = ArgumentError("Browser heartbeat timeout must be a positive finite number."); _webui_startup_failure!(_lifecycle_io, err, catch_backtrace(); phase = "validate_arguments"); throw(err))
  root = abspath(output_root === nothing ? default_webui_output_root() : String(output_root))
  paths = nothing
  try
    paths = _provision_webui_runtime!(root, config_file)
  catch err
    wrapped = ArgumentError("Could not create Web UI output directory $(root): $(sprint(showerror, err))")
    _webui_startup_failure!(_lifecycle_io, wrapped, catch_backtrace(); phase = "provision_runtime")
    throw(wrapped)
  end
  # Prune once at startup so normal requests keep append-only operation logging
  # unless a size safety cap is reached later.
  try
    _prune_webui_operation_log!(paths.operation_log; _webui_operation_log_options()...)
  catch err
    @warn "Could not prune Web UI operation log during startup" exception = (err, catch_backtrace())
  end
  _webui_startup_log(_lifecycle_io, "webui_start_requested"; operation_log = paths.operation_log, host = host_string, port = Int(port))
  config_error = _webui_validate_startup_config(paths.config_file)
  _webui_startup_log(_lifecycle_io, "webui_config_loaded"; operation_log = paths.operation_log, status = config_error === nothing ? "loaded" : "error_visible", config_file = paths.config_file, message = config_error)
  recovery = refresh_powerflow_run_registry!(root)
  _webui_startup_log(_lifecycle_io, "webui_routes_registered"; operation_log = paths.operation_log, status = "registered")
  record_webui_operation!(paths.operation_log, "webui_start"; route = "/powerflow", method = "START", status = "started", user_action = false, output_root = root, config_file = paths.config_file, case_cache_dir = paths.case_directory, operation_log = paths.operation_log)
  for result in get(recovery, "stale_recovered_runs", SparlectraApiResult[])
    record_webui_operation!(paths.operation_log, "webui_stale_active_run_recovered"; route = "/powerflow", method = "START", status = "interrupted_unknown", user_action = false, run_id = result.run_id, last_known_phase = get(result.metadata, "last_phase", nothing))
  end
  address = host_string == "localhost" ? ip"127.0.0.1" : parse(Sockets.IPAddr, host_string)
  listener = try
    Sockets.listen(address, UInt16(port))
  catch err
    if err isa Base.IOError && err.code == Base.UV_EADDRINUSE
      # The raw libuv IOError ("listen: address already in use (EADDRINUSE)")
      # gives no hint about *why*. The most common cause is a Sparlectra Web
      # UI already running (an earlier start in this session, or another
      # Julia process) — probe for it and, instead of aborting, hand the
      # user the running instance.
      if _webui_port_holds_sparlectra(host_string, port)
        url = "http://$(host_string == "localhost" ? "127.0.0.1" : host_string):$(Int(port))/powerflow"
        @info "Sparlectra Web UI is already running at $(url) — opening it instead of starting a second instance. To restart it, stop the running server first (close(server), the Stop Web UI button, or end its Julia process)."
        _webui_startup_log(_lifecycle_io, "webui_existing_instance_opened"; operation_log = paths.operation_log, status = "reused", host = host_string, port = Int(port))
        open_browser && _browser_opener(url)
        return nothing
      end
      wrapped = ArgumentError("Cannot start Sparlectra Web UI: $(host_string):$(port) is already in use by a non-Sparlectra process. Stop that process, or pass a different port=... to start_sparlectra_webui.")
      _webui_startup_failure!(_lifecycle_io, wrapped, catch_backtrace(); operation_log = paths.operation_log, phase = "bind")
      throw(wrapped)
    end
    _webui_startup_failure!(_lifecycle_io, err, catch_backtrace(); operation_log = paths.operation_log, phase = "bind")
    rethrow()
  end
  _webui_startup_log(_lifecycle_io, "webui_server_bound"; operation_log = paths.operation_log, status = "bound", host = host_string, port = Int(port))
  effective_shutdown_on_browser_close = auto_shutdown_on_browser_close === nothing ? shutdown_on_browser_close : Bool(auto_shutdown_on_browser_close)
  # The library default stays "no warm-up" (a programmatic caller should not
  # silently pay compile runs); `webui.warmup` governs the user-facing
  # launchers, which read it and pass the value explicitly.
  effective_warmup = warmup === nothing ? false : Bool(warmup)
  runtime = _SparlectraWebUIRuntime(listener, paths.case_directory, paths.config_file, paths.operation_log, config_error, _test_runner, effective_shutdown_on_browser_close, false, 0.0, 0, nothing, _lifecycle_io, ReentrantLock(), effective_warmup ? :waiting_first_page : :disabled)
  task = @async begin
    try
      while isopen(listener)
        try
          socket = accept(listener)
          @async _webui_serve_client(socket, root, runtime)
        catch err
          isopen(listener) && @error "Sparlectra Web UI accept loop failed" exception = (err, catch_backtrace())
        end
      end
    finally
      reason = lock(runtime.lock) do
        runtime.shutdown_reason
      end
      record_webui_operation!(runtime.operation_log, "webui_stopped"; route = "/powerflow", method = "STOP", status = "stopped", user_action = false, reason = string(something(reason, :process_exit)))
      _webui_lifecycle_println(runtime, "Sparlectra Web UI stopped$(_webui_shutdown_reason_text(reason)).")
    end
  end
  url_host = host_string == "::1" ? "[::1]" : host_string
  url = "http://$(url_host):$(port)/powerflow"
  # Pre-warm the case-selector scan in the background: the first form render
  # otherwise pays the content classification of the whole case directory
  # (large CGMES ZIPs make that seconds). The memo cache in forms.jl keeps
  # every later render fast; failures stay silent, startup must not depend
  # on the cache directory being readable.
  @async try
    # short delay so the cooperative scheduler returns control to the
    # startup path first: the server banner and the browser open within a
    # couple of seconds, and the one-time scan then runs in the background
    # (or inside the first form request, whichever comes first; the memo
    # lock shares the per-file results between the two)
    sleep(1.0)
    _webui_casefile_options_in_directory(paths.case_directory)
    _webui_for002_reference_options_in_directory(paths.case_directory)
  catch
  end
  browser_monitor_task = nothing
  if effective_warmup
    @async try
      # Serve the warm-up page before the CPU-bound warm-up solve starts: on
      # Julia's single-threaded cooperative scheduler a busy warm-up task
      # delays the HTTP response to the first browser request, leaving the
      # user staring at a blank tab. Wait until the warm-up page has actually
      # been rendered once (with a fallback timeout in case no browser
      # connects), then give the browser a short grace period to paint it.
      waited = 0.0
      while waited < 30.0 && !_webui_warmup_page_served(runtime)
        sleep(0.2)
        waited += 0.2
      end
      delay = Float64(warmup_delay_seconds)
      isfinite(delay) && delay > 0 && sleep(delay)
      _webui_set_warmup_state!(runtime, :warming)
      warmup_result = _run_sparlectra_webui_warmup(root; warmup_casefile, warmup_store_result)
      warmup_result.success || @warn "Sparlectra Web UI warm-up run did not converge" reason = warmup_result.reason message = warmup_result.message
      # Pre-compile the PowerFlow form render path too, so the automatic
      # refresh that replaces the warm-up page does not pay the first-render
      # JIT cost while the user is waiting.
      render_powerflow_form(output_root = root, case_directory = paths.case_directory, operation_log = paths.operation_log, selected_config_file = paths.config_file)
    catch err
      @warn "Sparlectra Web UI warm-up failed; normal runs remain available" exception = (err, catch_backtrace())
    finally
      _webui_set_warmup_state!(runtime, :done)
    end
  end
  if open_browser
    browser_process = _browser_opener(url)
    if browser_process !== nothing && effective_shutdown_on_browser_close
      record_webui_operation!(paths.operation_log, "browser_close_monitor_skipped"; route = "/powerflow", method = "START", status = "skipped", user_action = false, reason = "not_reliable_on_this_platform")
    end
  end
  server = SparlectraWebUIServer(listener, task, url, runtime, browser_monitor_task)
  record_webui_operation!(paths.operation_log, "webui_started"; route = "/powerflow", method = "START", status = "started", user_action = false, output_root = root, config_file = paths.config_file, case_cache_dir = paths.case_directory, operation_log = paths.operation_log)
  _webui_lifecycle_println(runtime, "Sparlectra Web UI is available at ", url)
  _webui_lifecycle_println(runtime, "Stop: use Stop Web UI in the browser, close(server), or Ctrl+C here.")
  _webui_lifecycle_println(runtime, "Operation log: ", paths.operation_log)
  @info "Sparlectra Web UI started" url output_root = abspath(root)
  return server
end

function _wait_sparlectra_webui(server::SparlectraWebUIServer; wait_for_task = wait)::Bool
  interrupted = false
  try
    wait_for_task(server.task)
  catch err
    err isa InterruptException || rethrow()
    interrupted = true
  finally
    _webui_request_shutdown!(server.runtime; reason = interrupted ? :ctrl_c : :server_closed)
  end
  return interrupted
end

function Base.close(server::SparlectraWebUIServer)
  _webui_request_shutdown!(server.runtime; reason = :server_closed)
  return nothing
end
