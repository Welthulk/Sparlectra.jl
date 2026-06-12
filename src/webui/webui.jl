using Markdown
using Sockets

mutable struct _SparlectraWebUIRuntime
  listener::Union{Sockets.TCPServer,Nothing}
  auto_shutdown_on_browser_close::Bool
  heartbeat_timeout_seconds::Float64
  heartbeat_received::Bool
  last_heartbeat::Float64
  active_requests::Int
  lock::ReentrantLock
end

include("forms.jl")
include("docs.jl")
include("views.jl")
include("handlers.jl")
include("routes.jl")

mutable struct SparlectraWebUIServer
  listener::Sockets.TCPServer
  task::Task
  url::String
  runtime::_SparlectraWebUIRuntime
  heartbeat_task::Task
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
  body = length_bytes > 0 ? String(read(socket, length_bytes)) : ""
  form = get(headers, "content-type", "") |> lowercase |> content_type -> startswith(content_type, "application/x-www-form-urlencoded") ? _webui_parse_pairs(body) : Dict{String,String}()
  return method, target, form
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

function _webui_app_command(url::String; platform::Symbol = _webui_platform(), executable_lookup = Sys.which, path_exists = ispath, environment = ENV)::Union{Cmd,Nothing}
  if platform == :macos
    applications = (
      "/Applications/Microsoft Edge.app",
      "/Applications/Google Chrome.app",
      "/Applications/Chromium.app",
      "/Applications/Brave Browser.app",
    )
    application = findfirst(path_exists, applications)
    application === nothing && return nothing
    return `open -na $(applications[application]) --args --app=$url`
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
  return `$executable --app=$url`
end

function _webui_open_browser(url::String)
  command = _webui_app_command(url)
  if command === nothing
    @warn "Could not find an app-capable browser; open the Web UI URL manually" url
    return nothing
  end
  try
    run(command; wait = false)
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

function _webui_request_shutdown!(runtime::_SparlectraWebUIRuntime)
  listener = lock(runtime.lock) do
    runtime.listener
  end
  listener === nothing || (isopen(listener) && close(listener))
  return nothing
end

function _webui_monitor_heartbeat(runtime::_SparlectraWebUIRuntime)
  while true
    listener, should_stop = lock(runtime.lock) do
      current_listener = runtime.listener
      expired = runtime.auto_shutdown_on_browser_close && runtime.heartbeat_received && runtime.active_requests == 0 && time() - runtime.last_heartbeat > runtime.heartbeat_timeout_seconds
      (current_listener, expired)
    end
    (listener === nothing || !isopen(listener)) && return nothing
    if should_stop
      _webui_request_shutdown!(runtime)
      return nothing
    end
    sleep(min(0.25, runtime.heartbeat_timeout_seconds / 4))
  end
end

"""
    start_sparlectra_webui(; host="127.0.0.1", port=8080,
                            output_root="results/powerflow_service",
                            open_browser=false,
                            auto_shutdown_on_browser_close=true,
                            browser_heartbeat_timeout_seconds=15.0,
                            warmup=false, warmup_casefile=nothing,
                            warmup_store_result=false) -> SparlectraWebUIServer

Start the loopback-only PowerFlow interface and load its persistent run registry
before accepting requests. The returned handle can be stopped with
`close(server)` or the browser's **Stop Web UI** button. When
`auto_shutdown_on_browser_close=true`, the server stops after the first browser
heartbeat has been seen and then expires; a server started without a browser
continues running because the timeout is not armed before that first heartbeat.
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

function start_sparlectra_webui(; host::AbstractString = "127.0.0.1", port::Integer = 8080, output_root::AbstractString = "results/powerflow_service", open_browser::Bool = false, auto_shutdown_on_browser_close::Bool = true, browser_heartbeat_timeout_seconds::Real = 15.0, warmup::Bool = false, warmup_casefile::Union{Nothing,AbstractString} = nothing, warmup_store_result::Bool = false)::SparlectraWebUIServer
  host_string = String(host)
  host_string in ("127.0.0.1", "localhost", "::1") || throw(ArgumentError("Sparlectra Web UI only accepts loopback hosts: 127.0.0.1, localhost, or ::1."))
  1 <= port <= 65535 || throw(ArgumentError("Web UI port must be between 1 and 65535."))
  timeout = Float64(browser_heartbeat_timeout_seconds)
  isfinite(timeout) && timeout > 0 || throw(ArgumentError("Browser heartbeat timeout must be a positive finite number."))
  root = String(output_root)
  refresh_powerflow_run_registry!(root)
  address = host_string == "localhost" ? ip"127.0.0.1" : parse(Sockets.IPAddr, host_string)
  listener = Sockets.listen(address, UInt16(port))
  runtime = _SparlectraWebUIRuntime(listener, auto_shutdown_on_browser_close, timeout, false, 0.0, 0, ReentrantLock())
  task = @async begin
    while isopen(listener)
      try
        socket = accept(listener)
        @async _webui_serve_client(socket, root, runtime)
      catch err
        isopen(listener) && @error "Sparlectra Web UI accept loop failed" exception = (err, catch_backtrace())
      end
    end
  end
  heartbeat_task = @async _webui_monitor_heartbeat(runtime)
  url_host = host_string == "::1" ? "[::1]" : host_string
  url = "http://$(url_host):$(port)/powerflow"
  server = SparlectraWebUIServer(listener, task, url, runtime, heartbeat_task)
  if warmup
    @async try
      warmup_result = _run_sparlectra_webui_warmup(root; warmup_casefile, warmup_store_result)
      warmup_result.success || @warn "Sparlectra Web UI warm-up run did not converge" reason = warmup_result.reason message = warmup_result.message
    catch err
      @warn "Sparlectra Web UI warm-up failed; normal runs remain available" exception = (err, catch_backtrace())
    end
  end
  open_browser && _webui_open_browser(url)
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
    close(server)
  end
  return interrupted
end

function Base.close(server::SparlectraWebUIServer)
  _webui_request_shutdown!(server.runtime)
  return nothing
end
