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
  git_head = joinpath(root, ".git", "HEAD")
  isfile(git_head) || return nothing
  head = strip(read(git_head, String))
  if startswith(head, "ref:")
    ref = strip(chop(head; head = 4, tail = 0))
    ref_path = joinpath(root, ".git", split(ref, '/')...)
    isfile(ref_path) || return nothing
    sha = strip(read(ref_path, String))
    return isempty(sha) ? nothing : sha
  end
  return isempty(head) ? nothing : head
end

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
                            warmup_store_result=false) -> SparlectraWebUIServer

Start the loopback-only PowerFlow interface and load its persistent run registry
before accepting requests. The returned handle can be stopped with
`close(server)` or the browser's **Stop Web UI** button. Browser-process
lifetime is not used for automatic shutdown by default because common browsers
may return a short-lived launcher process instead of a reliably owned window.
When `warmup=true`, a hidden asynchronous run compiles the common import/API/
solver path. By default it uses the bundled original synthetic case and a
temporary output directory, so it does not add a run-history entry or retain
artifacts. Set `warmup_store_result=true` to retain warm-up artifacts beneath
the configured output root.
"""
function _run_sparlectra_webui_warmup(output_root::AbstractString; warmup_casefile::Union{Nothing,AbstractString} = nothing, warmup_store_result::Bool = false, runner = run_sparlectra_api)
  casefile = warmup_casefile === nothing ? joinpath(_WEBUI_PACKAGE_ROOT, "data", "webui", "warmup_case3.jl") : abspath(warmup_casefile)
  isfile(casefile) || throw(ArgumentError("Web UI warm-up case file not found: $(casefile)"))
  execute(output_dir) = runner(
    casefile = casefile,
    config_file = DEFAULT_SPARLECTRA_CONFIG_PATH,
    output_dir = output_dir,
    config_overrides = Dict("output.logfile_results" => "off", "benchmark.enabled" => false),
    performance_timing = :off,
    run_diagnostics = false,
  )
  if warmup_store_result
    output_dir = joinpath(abspath(output_root), "webui-warmup")
    mkpath(output_dir)
    return execute(output_dir)
  end
  return mktempdir() do output_dir
    execute(output_dir)
  end
end

function _provision_webui_runtime!(root::AbstractString, config_file::Union{Nothing,AbstractString})
  configuration = config_file === nothing ? default_webui_config_path(root) : abspath(config_file)
  case_directory = default_webui_case_cache_dir(root)
  operation_log = default_webui_operation_log_path(root)
  mkpath.(unique((abspath(root), dirname(configuration), case_directory, dirname(operation_log))))
  config_file === nothing && !isfile(configuration) && cp(DEFAULT_SPARLECTRA_CONFIG_PATH, configuration)
  isfile(configuration) || throw(ArgumentError("Web UI configuration file not found: $(configuration)"))
  source = joinpath(_WEBUI_PACKAGE_ROOT, "data", "webui", "warmup_case3.jl")
  if isfile(source)
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

function start_sparlectra_webui(; host::AbstractString = "127.0.0.1", port::Integer = 8080, output_root::Union{Nothing,AbstractString} = nothing, config_file::Union{Nothing,AbstractString} = nothing, open_browser::Bool = false, shutdown_on_browser_close::Bool = false, auto_shutdown_on_browser_close::Union{Nothing,Bool} = nothing, browser_heartbeat_timeout_seconds::Real = 15.0, warmup::Bool = false, warmup_casefile::Union{Nothing,AbstractString} = nothing, warmup_store_result::Bool = false, _test_runner = start_powerflow_run, _lifecycle_io::IO = stdout, _browser_opener = _webui_open_browser)::SparlectraWebUIServer
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
    _webui_startup_failure!(_lifecycle_io, err, catch_backtrace(); operation_log = paths.operation_log, phase = "bind")
    rethrow()
  end
  _webui_startup_log(_lifecycle_io, "webui_server_bound"; operation_log = paths.operation_log, status = "bound", host = host_string, port = Int(port))
  effective_shutdown_on_browser_close = auto_shutdown_on_browser_close === nothing ? shutdown_on_browser_close : Bool(auto_shutdown_on_browser_close)
  runtime = _SparlectraWebUIRuntime(listener, paths.case_directory, paths.config_file, paths.operation_log, config_error, _test_runner, effective_shutdown_on_browser_close, false, 0.0, 0, nothing, _lifecycle_io, ReentrantLock())
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
  browser_monitor_task = nothing
  if warmup
    @async try
      warmup_result = _run_sparlectra_webui_warmup(root; warmup_casefile, warmup_store_result)
      warmup_result.success || @warn "Sparlectra Web UI warm-up run did not converge" reason = warmup_result.reason message = warmup_result.message
    catch err
      @warn "Sparlectra Web UI warm-up failed; normal runs remain available" exception = (err, catch_backtrace())
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
