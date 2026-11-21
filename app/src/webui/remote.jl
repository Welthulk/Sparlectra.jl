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

# file: app/src/webui/remote.jl
# purpose: remote mode of the Web UI: how a server is reached (local, SSH
#          tunnel, direct bind), the connection lines for the console and the
#          Info menu, and the access token (direct mode, optional in SSH mode)

# The connection of the server whose request is being rendered. The page
# renderers get no runtime, so route_sparlectra_webui binds it for one
# request and the header (Info menu, Stop confirmation) and the footer read
# it. A scoped value, not a global: each request (and every task it spawns)
# sees its own server, also with two servers in one process. nothing = a
# render outside a request of a running server.
const _WEBUI_REQUEST_CONNECTION = Base.ScopedValues.ScopedValue{Union{Nothing,_WebUIConnection}}(nothing)

# Name of the cookie that keeps the access token after the first request.
# Browsers do not separate cookies by port (there is no port-scoped cookie),
# and SameSite=Strict counts every port of one host name as the same site:
# the cookie, and with it the token, also reaches every other HTTP service
# the browser opens on this host name. The port in the name only keeps two
# Web UI servers on one host from overwriting each other's token.
_webui_token_cookie_name(connection::_WebUIConnection)::String = string("sparlectra_webui_token_", connection.port)

# 128 random bits from the operating system, as 32 hex characters.
_webui_generate_access_token()::String = bytes2hex(rand(Random.RandomDevice(), UInt8, 16))

# Comparison whose time does not depend on where the first difference lies.
function _webui_tokens_equal(a::AbstractString, b::AbstractString)::Bool
  x, y = codeunits(a), codeunits(b)
  length(x) == length(y) || return false
  difference = 0x00
  for i in eachindex(x, y)
    difference |= x[i] ⊻ y[i]
  end
  return difference == 0x00
end

function _webui_cookie_value(headers::AbstractDict, name::AbstractString)::String
  for part in split(get(headers, "cookie", ""), ';')
    pair = split(strip(part), '='; limit = 2)
    length(pair) == 2 && pair[1] == name && return String(pair[2])
  end
  return ""
end

# The request target without its token parameter: the address the browser is
# sent to after the token was accepted, so the token does not stay in the
# address bar, the history or a Referer header.
function _webui_target_without_token(target::AbstractString)::String
  parts = split(String(target), '?'; limit = 2)
  length(parts) == 1 && return String(target)
  kept = [p for p in split(parts[2], '&') if !isempty(p) && !startswith(p, "token=")]
  return isempty(kept) ? String(parts[1]) : string(parts[1], "?", join(kept, "&"))
end

# The page a request without a valid token gets. It names nothing of the
# server (no case, no path, no version); "Sparlectra" stays in the text so a
# second start recognizes the running instance (_webui_port_holds_sparlectra).
const _WEBUI_ACCESS_DENIED_PAGE = """<!doctype html>
<html lang="en"><head><meta charset="utf-8"><title>Sparlectra Web UI</title></head>
<body><h1>Sparlectra Web UI: access token required</h1>
<p>This server accepts requests only with its access token.
Open the address with the token that the host console printed at start
(<code>http://&lt;address&gt;:&lt;port&gt;/powerflow?token=&lt;token&gt;</code>). The browser keeps the token afterwards.</p>
</body></html>"""

"""
    _webui_access_check(connection, method, target, headers) -> Union{Nothing,SparlectraWebUIResponse}

The access check of a server with an access token (direct mode, and SSH
mode with `webui.remote_ssh_require_token`), run before a request is routed.
Returns `nothing` when the request may pass: the connection has no token
(local mode, SSH mode without token), or the request carries the token
cookie. A `GET` with the right `token` query parameter gets a `303` redirect
to the same address without the token and a cookie that holds it
(`HttpOnly`, `SameSite=Strict`), so the token is typed once. The cookie is
not limited to the port: the browser also sends it to every other HTTP
service on the same host name. Every other request gets `401` with a page
that names nothing of the server.
"""
function _webui_access_check(connection::_WebUIConnection, method::AbstractString, target::AbstractString, headers::AbstractDict)
  token = connection.access_token
  isempty(token) && return nothing
  cookie_name = _webui_token_cookie_name(connection)
  _webui_tokens_equal(_webui_cookie_value(headers, cookie_name), token) && return nothing
  _, query = _webui_split_target(target)
  if uppercase(String(method)) == "GET" && _webui_tokens_equal(String(get(query, "token", "")), token)
    cookie = string(cookie_name, "=", token, "; Path=/; HttpOnly; SameSite=Strict")
    return SparlectraWebUIResponse(303, ""; headers = ["Location" => _webui_target_without_token(target), "Set-Cookie" => cookie])
  end
  return SparlectraWebUIResponse(401, _WEBUI_ACCESS_DENIED_PAGE)
end

# The client address for the operation log line of a refused request.
function _webui_peer_text(socket::Sockets.TCPSocket)::String
  try
    address, port = Sockets.getpeername(socket)
    return string(address, ":", Int(port))
  catch err
    # expected failure: the client closed the connection already
    err isa Base.IOError || rethrow()
    return "unknown"
  end
end

_webui_is_loopback(host::AbstractString)::Bool = host in ("127.0.0.1", "localhost", "::1")

# An address as it appears in a URL or an ssh -L forward: IPv6 in brackets.
_webui_url_host(address::AbstractString)::String = occursin(':', address) ? string("[", address, "]") : String(address)

"""
    _webui_resolve_connection(webui; remote, host, port, hostname) -> _WebUIConnection

The connection a start uses. `webui` is the `WebUIConfig` of the server's
configuration file. Without `remote` the server is local: `host` (default
`127.0.0.1`, loopback only) and `port` (default 8080). With `remote` the
configuration decides: `webui.remote_mode = ssh` keeps a loopback address,
`direct` binds `webui.remote_bind_address` and carries the access token
(`webui.remote_access_token`, or a new one when that is empty). The SSH mode
carries the token too when `webui.remote_ssh_require_token` is set. The port is
`webui.remote_port` in both remote modes. Explicit `host` and `port`
keywords win over the configuration. Throws an `ArgumentError` that names
the key for an address that is not loopback in local or SSH mode, or not an
IP address in direct mode.
"""
function _webui_resolve_connection(webui; remote::Bool, host::Union{Nothing,AbstractString}, port::Union{Nothing,Integer}, hostname::AbstractString = gethostname())::_WebUIConnection
  if !remote
    bind = host === nothing ? "127.0.0.1" : String(host)
    _webui_is_loopback(bind) || throw(ArgumentError("Sparlectra Web UI only accepts loopback hosts: 127.0.0.1, localhost, or ::1. For access from another machine start the remote mode (start_webui.jl --remote, see the Web UI documentation)."))
    return _WebUIConnection(mode = :local, bind_address = bind, port = port === nothing ? 8080 : Int(port), hostname = String(hostname))
  end
  remote_port = port === nothing ? webui.remote_port : Int(port)
  if webui.remote_mode === :direct
    bind = host === nothing ? webui.remote_bind_address : String(host)
    try
      bind == "localhost" || parse(Sockets.IPAddr, bind)
    catch err
      # expected failure: a host name or a typo instead of an IP address
      err isa ArgumentError || rethrow()
      throw(ArgumentError("webui.remote_bind_address must be an IP address (0.0.0.0 for every interface); got $(bind)."))
    end
    return _WebUIConnection(mode = :direct, bind_address = bind, port = remote_port, hostname = String(hostname), access_token = _webui_configured_token(webui))
  end
  bind = host === nothing ? "127.0.0.1" : String(host)
  _webui_is_loopback(bind) || throw(ArgumentError("webui.remote_mode = ssh keeps the server on a loopback host (127.0.0.1, localhost, ::1); got $(bind). Set webui.remote_mode: direct to bind a network address."))
  # a loopback address is open to every user and process on the host; the
  # token keeps them out where the host is shared
  token = webui.remote_ssh_require_token ? _webui_configured_token(webui) : ""
  return _WebUIConnection(mode = :ssh, bind_address = bind, port = remote_port, hostname = String(hostname), access_token = token)
end

# The configured token, or a new one per start when the key is empty.
_webui_configured_token(webui)::String = isempty(webui.remote_access_token) ? _webui_generate_access_token() : webui.remote_access_token

"""
    _webui_connection_info(connection) -> NamedTuple

The connection data of a server: `hostname`, `bind_address`, `port`, `mode`
(`:local`, `:ssh`, `:direct`), `url` (the address to open; in SSH mode the
end of the tunnel on the user's machine, in direct mode the server address;
never with the token), `ssh_command` (the tunnel command with a `<user>`
placeholder, `""` in direct mode) and `token_required` (the server asks for
the access token).
"""
function _webui_connection_info(connection::_WebUIConnection)
  host = isempty(connection.hostname) ? gethostname() : connection.hostname
  port = connection.port
  if connection.mode === :direct
    # a wildcard bind is reached through the host name, a single interface
    # through its own address
    url_host = connection.bind_address in ("0.0.0.0", "::") ? host : _webui_url_host(connection.bind_address)
    return (hostname = host, bind_address = connection.bind_address, port = port, mode = :direct, url = "http://$(url_host):$(port)/powerflow", ssh_command = "", token_required = !isempty(connection.access_token))
  end
  # the tunnel forwards to the address the server listens on; the browser
  # opens the near end of the tunnel on the user's machine
  target = connection.bind_address == "localhost" ? "localhost" : _webui_url_host(connection.bind_address)
  ssh_command = "ssh -N -L $(port):$(target):$(port) <user>@$(host)"
  url = connection.mode === :ssh ? "http://127.0.0.1:$(port)/powerflow" : "http://$(target):$(port)/powerflow"
  return (hostname = host, bind_address = connection.bind_address, port = port, mode = connection.mode, url = url, ssh_command = ssh_command, token_required = !isempty(connection.access_token))
end

"""
    webui_connection_info(server) -> NamedTuple

How a running Web UI is reached, the same data the Info menu shows under
**Server**: `hostname`, `bind_address`, `port`, `mode` (`:local` for a
normal start, `:ssh` and `:direct` for the two remote modes), `url`,
`ssh_command` (tunnel command with a `<user>` placeholder, empty in direct
mode) and `token_required` (direct mode, or SSH mode with
`webui.remote_ssh_require_token`). `url` never carries the token; the token
is printed on the host console at start.
"""
webui_connection_info(server::SparlectraWebUIServer) = _webui_connection_info(server.runtime.connection)

function _webui_mode_text(connection::_WebUIConnection)::String
  connection.mode === :direct && return "remote, direct (access token, not encrypted)"
  connection.mode === :ssh && return isempty(connection.access_token) ? "remote, SSH tunnel (no token: every user and process on the host can connect)" : "remote, SSH tunnel with access token"
  return "local, loopback access only"
end

# The footer names how the page's server is reached; a render outside a
# running server claims nothing.
function _webui_footer_mode_text(connection::Union{Nothing,_WebUIConnection} = _WEBUI_REQUEST_CONNECTION[])::String
  connection === nothing && return "PowerFlow Web UI"
  connection.mode === :local && return "Local PowerFlow Web UI · loopback access only"
  host = isempty(connection.hostname) ? gethostname() : connection.hostname
  return string("PowerFlow Web UI on ", host, " · ", _webui_mode_text(connection))
end

# The lines the host console prints at the start of a remote mode. Only the
# console gets the token: it is the one place a person on the host reads.
function _webui_connection_lines(connection::_WebUIConnection)::Vector{String}
  info = _webui_connection_info(connection)
  if connection.mode === :ssh
    info.token_required || return ["Remote mode (SSH tunnel) on $(info.hostname), port $(info.port). On your own machine run",
      "    $(info.ssh_command)",
      "and open $(info.url)",
      "Every user with a login on this host can open the port too (webui.remote_ssh_require_token: true asks for a token)."]
    return ["Remote mode (SSH tunnel, access token) on $(info.hostname), port $(info.port). On your own machine run",
      "    $(info.ssh_command)",
      "and open",
      "    $(info.url)?token=$(connection.access_token)",
      "The browser keeps the token after the first request."]
  elseif connection.mode === :direct
    return ["Remote mode (direct) on $(info.hostname), listening on $(_webui_url_host(info.bind_address)):$(info.port). Open",
      "    $(info.url)?token=$(connection.access_token)",
      "The browser keeps the token after the first request. The traffic is not encrypted (plain HTTP): use this mode only in a network you trust."]
  end
  return String[]
end

# Placeholder in the Info menu for the server section: the menu of most pages
# is memoized (_webui_default_info_menu), the server section must follow the
# server that renders the page, so _webui_layout fills it in on every render.
const _WEBUI_SERVER_INFO_SLOT = "<!-- sparlectra-server-info -->"

function _webui_server_info_html(connection::Union{Nothing,_WebUIConnection} = _WEBUI_REQUEST_CONNECTION[])::String
  connection === nothing && return ""
  info = _webui_connection_info(connection)
  rows = String[
    "<dt>Host</dt><dd><code>$(_webui_escape(info.hostname))</code></dd>",
    "<dt>Bind address</dt><dd><code>$(_webui_escape(info.bind_address))</code></dd>",
    "<dt>Port</dt><dd><code>$(info.port)</code></dd>",
    "<dt>Mode</dt><dd>$(_webui_escape(_webui_mode_text(connection)))</dd>",
  ]
  # the tunnel row belongs to the SSH mode only: a local start is used on
  # its own host, a direct bind needs no tunnel
  info.mode === :ssh && push!(rows, "<dt>SSH tunnel</dt><dd><code>$(_webui_escape(info.ssh_command))</code></dd>")
  where_text = info.mode === :ssh ? "on your own machine: " : ""
  token_text = info.token_required ? "?token=&lt;token&gt;</code>, the token is on the host console" : "</code>"
  push!(rows, "<dt>URL</dt><dd>$(where_text)<code>$(_webui_escape(info.url))$(token_text)</dd>")
  docs = string(_webui_docs_base_url(), "webui/#webui_remote")
  return """<h2>Server</h2>
<dl>
$(join(rows, "\n"))
</dl>
<p class=\"field-hint\">Remote use and its settings: <a href=\"$(_webui_escape(docs))\" target=\"_blank\" rel=\"noopener noreferrer\">documentation</a>.</p>
"""
end

# The onsubmit attribute of the Stop Web UI form: in a remote mode other
# people may use the server, so the click asks first. Empty in local mode.
function _webui_stop_confirm_attribute(connection::Union{Nothing,_WebUIConnection} = _WEBUI_REQUEST_CONNECTION[])::String
  (connection === nothing || connection.mode === :local) && return ""
  # a host name holds letters, digits, dots and hyphens; anything else is
  # dropped so the text cannot break out of the JavaScript string
  host = filter(c -> isletter(c) || isdigit(c) || c in ".-_", isempty(connection.hostname) ? gethostname() : connection.hostname)
  return " onsubmit=\"return confirm('Stop the Sparlectra Web UI on $(host)? Every user of this server loses the connection.');\""
end
