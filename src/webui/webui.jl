using Markdown
using Sockets

include("forms.jl")
include("docs.jl")
include("views.jl")
include("handlers.jl")
include("routes.jl")

mutable struct SparlectraWebUIServer
  listener::Sockets.TCPServer
  task::Task
  url::String
end

function _webui_status_text(status::Int)
  return get(Dict(200 => "OK", 303 => "See Other", 400 => "Bad Request", 404 => "Not Found", 500 => "Internal Server Error"), status, "OK")
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

function _webui_serve_client(socket::Sockets.TCPSocket, output_root::String)
  try
    method, target, form = _webui_read_request(socket)
    response = route_sparlectra_webui(method, target, form; output_root = output_root)
    _webui_write_response(socket, response)
  catch err
    try
      _webui_write_response(socket, _webui_html(render_webui_error(500, sprint(showerror, err)); status = 500))
    catch
    end
  finally
    close(socket)
  end
  return nothing
end

function _webui_open_browser(url::String)
  command = Sys.iswindows() ? `cmd /c start $url` : Sys.isapple() ? `open $url` : `xdg-open $url`
  @async try
    run(command)
  catch err
    @warn "Could not open the default browser" url exception = (err, catch_backtrace())
  end
  return nothing
end

"""
    start_sparlectra_webui(; host="127.0.0.1", port=8080,
                            output_root="results/powerflow_service",
                            open_browser=false) -> SparlectraWebUIServer

Start the local browser-based PowerFlow interface. The server accepts loopback
hosts only (`127.0.0.1`, `localhost`, or `::1`), defaults to port `8080`, and
uses the existing PowerFlow service for execution, history, and artifact
resolution. The call returns immediately with a server handle; call
`close(server)` to stop it. Set `open_browser=true` to attempt opening the local
URL with the operating system's default browser.
"""
function start_sparlectra_webui(; host::AbstractString = "127.0.0.1", port::Integer = 8080, output_root::AbstractString = "results/powerflow_service", open_browser::Bool = false)::SparlectraWebUIServer
  host_string = String(host)
  host_string in ("127.0.0.1", "localhost", "::1") || throw(ArgumentError("Sparlectra Web UI only accepts loopback hosts: 127.0.0.1, localhost, or ::1."))
  1 <= port <= 65535 || throw(ArgumentError("Web UI port must be between 1 and 65535."))
  address = host_string == "localhost" ? ip"127.0.0.1" : parse(Sockets.IPAddr, host_string)
  listener = Sockets.listen(address, UInt16(port))
  root = String(output_root)
  task = @async begin
    while isopen(listener)
      try
        socket = accept(listener)
        @async _webui_serve_client(socket, root)
      catch err
        isopen(listener) && @error "Sparlectra Web UI accept loop failed" exception = (err, catch_backtrace())
      end
    end
  end
  url_host = host_string == "::1" ? "[::1]" : host_string
  url = "http://$(url_host):$(port)/powerflow"
  server = SparlectraWebUIServer(listener, task, url)
  open_browser && _webui_open_browser(url)
  @info "Sparlectra Web UI started" url output_root = abspath(root)
  return server
end

function Base.close(server::SparlectraWebUIServer)
  isopen(server.listener) && close(server.listener)
  return nothing
end
