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

# file: app/src/webui/update.jl
# purpose: "Check for updates" of the Info panel and of start_webui.jl: one
#          request to the latest-release page of GitHub on the click or at an
#          interactive start, the comparison with the running version, how
#          this copy was installed, and for a copy made by tools/install_webui.*
#          the unattended update in a detached process that waits for this
#          Web UI to stop

# The page GitHub redirects to the tag of the latest release. Not the REST
# API: that one allows 60 anonymous requests per hour per public IP, shared
# by every machine behind one company NAT, and the check then failed with
# "rate limit reached" for users who never clicked it before.
const WEBUI_UPDATE_LATEST_URL = "https://github.com/Welthulk/Sparlectra.jl/releases/latest"
const WEBUI_UPDATE_RELEASES_PAGE = "https://github.com/Welthulk/Sparlectra.jl/releases"
const WEBUI_UPDATE_TIMEOUT_SECONDS = 10.0
# the check at start waits less: an offline start should not stall
const WEBUI_CONSOLE_UPDATE_TIMEOUT_SECONDS = 5.0
# how long the question at start waits for an answer before it counts as no
const WEBUI_CONSOLE_UPDATE_ANSWER_SECONDS = 60.0
# start_webui.jl flag that skips the check at start (also SPARLECTRA_UPDATE_CHECK=0)
const WEBUI_CONSOLE_UPDATE_FLAG = "--no-update-check"
# the one-line install of docs/src/webui_reference.md, per platform
const WEBUI_UPDATE_ONE_LINE_SH = "curl -fsSL https://raw.githubusercontent.com/Welthulk/Sparlectra.jl/main/tools/install_webui.sh | sh"
const WEBUI_UPDATE_ONE_LINE_WINDOWS = "iwr -useb https://raw.githubusercontent.com/Welthulk/Sparlectra.jl/main/tools/install_webui.bat -OutFile install_webui.bat; .\\install_webui.bat"
# how long the updater waits for this Web UI process to end before it gives up
const WEBUI_UPDATE_WAIT_SECONDS = 600

"""
    _webui_update_http_get(url; timeout) -> NamedTuple

One HEAD request through the Downloads standard library, redirects
followed: `(status, headers, url, error, timed_out)`. `url` is the address
the redirects ended at. `status` is 0 when no HTTP answer came back (no
network, name resolution, TLS, timeout); `error` then carries the curl
message and `timed_out` says whether the time ran out. Header names are
lowercased.
"""
function _webui_update_http_get(url::AbstractString; timeout::Real = WEBUI_UPDATE_TIMEOUT_SECONDS)
  # throw = false: a transport failure comes back as a RequestError value
  # instead of an exception, so no catch is needed here
  answer = Downloads.request(String(url); method = "HEAD", output = devnull, timeout = Float64(timeout), throw = false)
  if answer isa Downloads.RequestError
    # curl code 28 is CURLE_OPERATION_TIMEDOUT
    return (status = 0, headers = Dict{String,String}(), url = "", error = String(answer.message), timed_out = answer.code == 28)
  end
  headers = Dict{String,String}(lowercase(String(first(h))) => String(last(h)) for h in answer.headers)
  return (status = Int(answer.status), headers = headers, url = String(answer.url), error = "", timed_out = false)
end

# The fetch the check uses. Tests put a function of the same shape here, so
# no test ever reaches the network.
const _WEBUI_UPDATE_FETCH = Ref{Any}(_webui_update_http_get)

"""
    _webui_update_release_from_url(url) -> Union{Nothing,NamedTuple}

`(version, tag, page)` from the address the latest-release page redirected
to, or `nothing` when it is not a release tag page of the repository of the
form `.../releases/tag/vX.Y.Z` (a repository without a release ends on
`.../releases`, a proxy or login page somewhere else).
"""
function _webui_update_release_from_url(url::AbstractString)
  m = match(r"^https://github\.com/Welthulk/Sparlectra\.jl/releases/tag/(v(\d+\.\d+\.\d+))$", url)
  m === nothing && return nothing
  return (version = VersionNumber(m.captures[2]), tag = String(m.captures[1]), page = String(url))
end

"""
    webui_check_for_update(; fetch, current, url, timeout) -> NamedTuple

Ask GitHub for the latest release (the redirect of its latest-release page,
not the REST API) and compare it with the
running version (`current`, default the Sparlectra version). Called by
the "Check for updates" button and by `webui_console_update_offer`
at an interactive start. Returns `state` (`:up_to_date`, `:newer` or
`:failed`), `current`, `latest` (`nothing` on failure), `tag`,
`release_url`, and for a failure `failure` (`:network`, `:timeout`,
`:rate_limit`, `:http`, `:malformed`) with a `message` naming it. Never
throws for a network or answer problem.
"""
function webui_check_for_update(; fetch = _WEBUI_UPDATE_FETCH[], current::VersionNumber = version(), url::AbstractString = WEBUI_UPDATE_LATEST_URL, timeout::Real = WEBUI_UPDATE_TIMEOUT_SECONDS)
  failed(failure, message) = (state = :failed, current = current, latest = nothing, tag = "", release_url = WEBUI_UPDATE_RELEASES_PAGE, failure = failure, message = String(message))
  answer = fetch(url; timeout = timeout)
  if answer.status == 0
    answer.timed_out && return failed(:timeout, "no answer from GitHub within $(round(Int, timeout)) s (timeout)")
    return failed(:network, "no network connection to github.com ($(answer.error))")
  end
  # the release pages carry no API quota; a 429 is GitHub (or a proxy)
  # asking to slow down, and it says so
  answer.status == 429 && return failed(:rate_limit, "GitHub answered HTTP 429 (too many requests), try again later")
  answer.status == 404 && return failed(:http, "GitHub lists no release for Welthulk/Sparlectra.jl (HTTP 404)")
  answer.status == 200 || return failed(:http, "GitHub answered HTTP $(answer.status)")
  release = _webui_update_release_from_url(answer.url)
  release === nothing && return failed(:malformed, "GitHub did not lead to a release tag of the form vX.Y.Z (ended at $(isempty(answer.url) ? "no address" : answer.url))")
  state = release.version > current ? :newer : :up_to_date
  return (state = state, current = current, latest = release.version, tag = release.tag, release_url = release.page, failure = :none, message = "")
end

"""
    webui_install_layout(root = SPARLECTRA_ROOT; flavor_kind) -> Symbol

How the running copy was installed, which decides what "Install update"
may do:

  * `:git_checkout`: `root` holds `.git` (a directory, or the pointer file
    of a git worktree). A development checkout is updated by hand.
  * `:installed_copy`: no `.git`, the folder is named `Sparlectra` and holds
    `Project.toml`, `start_webui.jl` and the install script of this
    platform under `tools/`. That is the layout `tools/install_webui.*`
    creates from a release download and updates in place (old copy to
    `Sparlectra.old`).
  * `:unknown`: anything else (the standalone app, a registry installation,
    a renamed download).
"""
function webui_install_layout(root::AbstractString = SPARLECTRA_ROOT; flavor_kind::Symbol = webui_runtime_flavor().kind)::Symbol
  # the standalone executable is rebuilt, not updated by the script
  flavor_kind === :app && return :unknown
  dir = rstrip(String(root), ['/', '\\'])
  ispath(joinpath(dir, ".git")) && return :git_checkout
  basename(dir) == "Sparlectra" || return :unknown
  all(isfile, (joinpath(dir, "Project.toml"), joinpath(dir, "start_webui.jl"), joinpath(dir, "tools", _webui_update_installer_name()))) || return :unknown
  return :installed_copy
end

_webui_update_installer_name(; windows::Bool = Sys.iswindows())::String = windows ? "install_webui.bat" : "install_webui.sh"

# The install action needs this process to END after the shutdown, because
# the updater waits for its process id. A Web UI started by start_webui.jl
# (directly, through start_webui.sh/.bat or after a sysimage relaunch) does;
# one started from a REPL keeps the Julia session alive.
_webui_started_by_start_script()::Bool = basename(PROGRAM_FILE) == "start_webui.jl"

# One random value per Web UI process. The install form carries it, and the
# install action refuses a POST without it: the server binds to loopback but
# has no other cross-site protection, and any page in the same browser can
# post a form to 127.0.0.1. A new process (after the update) has a new value.
# The value is kept with the process id: a value made while a package image
# or a sysimage was built would otherwise be shared by every later process.
const _WEBUI_UPDATE_TOKEN = Ref((0, ""))
function _webui_update_token()::String
  pid, token = _WEBUI_UPDATE_TOKEN[]
  pid == getpid() && return token
  token = bytes2hex(rand(Random.RandomDevice(), UInt8, 16))
  _WEBUI_UPDATE_TOKEN[] = (Int(getpid()), token)
  return token
end

"""
    webui_update_dir(operation_log) -> String

The folder of the update files next to the operation log: the updater
script, the copy of the install script it runs, and `update.log`.
"""
webui_update_dir(operation_log::AbstractString)::String = joinpath(dirname(webui_operation_log_path(operation_log)), "update")

webui_update_log_path(operation_log::AbstractString)::String = joinpath(webui_update_dir(operation_log), "update.log")

# the result of the last update, written by the updater after the install
# script ended: "ok" or "failed", the tag, the reason, the time
webui_update_status_path(operation_log::AbstractString)::String = joinpath(webui_update_dir(operation_log), "last_update.txt")

"""
    webui_last_update(operation_log) -> Union{Nothing,NamedTuple}

The result of the last update the Web UI started, read from the file the
updater writes when the install script has ended: `(ok, tag, reason,
finished, log)`, or `nothing` when no update ran yet. A failed update stays
visible (Info panel, start in a terminal) until a later update succeeds.
"""
function webui_last_update(operation_log::AbstractString)
  path = webui_update_status_path(operation_log)
  isfile(path) || return nothing
  # strip: the Windows updater writes CRLF and the reason after a blank
  lines = [String(strip(l)) for l in readlines(path)]
  isempty(lines) && return nothing
  line(i) = i <= length(lines) ? lines[i] : ""
  return (ok = line(1) == "ok", tag = line(2), reason = line(3), finished = line(4), log = webui_update_log_path(operation_log))
end

# the Info panel line of a failed last update (nothing after a good one)
function _webui_info_last_update(operation_log::AbstractString)::String
  text = _webui_last_update_failure_text(operation_log)
  isempty(text) && return ""
  return "<div class=\"alert alert-error update-last-failed\" role=\"alert\">$(_webui_escape(text))</div>"
end

# one line on a failed last update, "" otherwise
function _webui_last_update_failure_text(operation_log::AbstractString)::String
  last = webui_last_update(operation_log)
  (last === nothing || last.ok) && return ""
  to = isempty(last.tag) ? "" : " to $(last.tag)"
  at = isempty(last.finished) ? "" : " ($(last.finished))"
  reason = isempty(last.reason) ? "no reason recorded" : last.reason
  return "The last update$(to) failed$(at): $(reason). This copy was not changed by it. Log: $(last.log)"
end

"""
    _webui_update_answers(; flavor_kind, sysimage_path, app_root, home, windows) -> NamedTuple

The answers of the unattended install, taken from the current state so that
the update forces nothing new: `build_sysimage` when this session runs on
the sysimage or an image of this checkout is on disk (otherwise the stale
image would be removed at the next start without a terminal), and
`create_shortcut` when the launcher the install script writes exists.
"""
function _webui_update_answers(; flavor_kind::Symbol = webui_runtime_flavor().kind, sysimage_path::AbstractString = webui_sysimage_path(), app_root::AbstractString = SPARLECTRA_APP_ROOT, home::AbstractString = homedir(), windows::Bool = Sys.iswindows(), apple::Bool = Sys.isapple())
  build = flavor_kind === :sysimage || _webui_update_own_sysimage(sysimage_path, app_root)
  launchers = if windows
    desktops = [joinpath(home, "Desktop")]
    onedrive = get(ENV, "OneDrive", "")
    isempty(onedrive) || push!(desktops, joinpath(onedrive, "Desktop"))
    [joinpath(d, "Sparlectra Web UI.lnk") for d in desktops]
  elseif apple
    [joinpath(home, "Desktop", "Sparlectra Web UI")]
  else
    [joinpath(home, ".local", "share", "applications", "sparlectra-webui.desktop"), joinpath(home, "Desktop", "sparlectra-webui.desktop"), joinpath(home, "sparlectra-webui")]
  end
  return (build_sysimage = build, create_shortcut = any(p -> ispath(p) || islink(p), launchers))
end

# an image on disk whose metadata names this checkout as its owner
function _webui_update_own_sysimage(image::AbstractString, app_root::AbstractString)::Bool
  meta_path = joinpath(dirname(image), "sysimage_meta.toml")
  (isfile(image) && isfile(meta_path)) || return false
  meta = try
    TOML.parsefile(meta_path)
  catch err
    # expected failure: a truncated or hand-edited metadata file; such an
    # image is not used at the next start anyway
    err isa TOML.ParserError || rethrow()
    return false
  end
  owner = get(meta, "project_dir", "")
  owner isa AbstractString && !isempty(owner) || return false
  same(dir) = rstrip(isdir(dir) ? realpath(dir) : normpath(abspath(dir)), ['/', '\\'])
  return same(owner) == same(app_root)
end

_webui_sh_quote(s::AbstractString)::String = string("'", replace(String(s), "'" => "'\\''"), "'")
# inside "..." of a batch file only % needs escaping
_webui_bat_text(s::AbstractString)::String = replace(String(s), "%" => "%%")

"""
    webui_update_plan(; root, current, tag, operation_log, pid, answers, windows) -> NamedTuple

Write the updater script and a copy of the install script into
`webui_update_dir(operation_log)`, start a fresh `update.log` there, and
return `(command, updater, installer, log, install_dir, root, answers)`.
Nothing is started: `command` is the detached command
`start_webui_update!` runs. The updater file carries the time of the
request in its name, one file per update.

The updater waits until the process `pid` (this Web UI) has ended, runs the
copy of the install script with `SPARLECTRA_UPDATE=1`, the release `tag`
of the check as `SPARLECTRA_UPDATE_TAG` (the install script then downloads
that tag and does not ask GitHub again), the sysimage and
shortcut answers of `answers`, `SPARLECTRA_INSTALL_DIR` (the folder that
holds `Sparlectra`) and `SPARLECTRA_START_WEBUI=0`, and then starts the
Web UI through `start_webui.sh` / `start_webui.bat` of the copy. The install
script itself keeps the previous copy as `Sparlectra.old`. A copy of the
install script is run because the original moves with the old copy; on
Windows cmd.exe reads a batch file by its path line by line and would
continue in the new file.

On Linux and macOS the updater runs without a terminal, all output (also of
the restarted Web UI) goes to `update.log`. On Windows it opens a console
window of its own; it first waits until Windows lets the folder be moved
(the console of the old Web UI holds it as working directory until it is
closed) and writes the output of the install script to `update.log`.
"""
function webui_update_plan(; root::AbstractString = SPARLECTRA_ROOT, current::VersionNumber = version(), tag::AbstractString, operation_log::AbstractString, pid::Integer = getpid(), answers = _webui_update_answers(), windows::Bool = Sys.iswindows(), requested_in::AbstractString = "in the Web UI")
  dir = rstrip(String(root), ['/', '\\'])
  install_dir = dirname(dir)
  update_dir = webui_update_dir(operation_log)
  mkpath(update_dir)
  log = joinpath(update_dir, "update.log")
  status = webui_update_status_path(operation_log)
  # a result from an earlier update must not stand for this one
  rm(status; force = true)
  installer = joinpath(update_dir, _webui_update_installer_name(; windows))
  cp(joinpath(dir, "tools", _webui_update_installer_name(; windows)), installer; force = true)
  requested = Dates.now()
  stamp = Dates.format(requested, "yyyy-mm-dd HH:MM:SS")
  # one updater file per update: on Windows the console of an earlier
  # update runs the Web UI from its updater and reads that batch file by its
  # path again when the Web UI ends; rewriting it would run foreign lines
  file_stamp = Dates.format(requested, "yyyymmdd_HHMMSS")
  write(log, "Sparlectra update requested $(requested_in) at $(stamp): v$(current) to $(tag), copy $(dir)\n")
  b(x) = x ? "1" : "0"
  if windows
    updater = joinpath(update_dir, "sparlectra_update_$(file_stamp).bat")
    q(s) = _webui_bat_text(s)
    lines = String[
      "@echo off",
      "rem Sparlectra update, written by the Web UI (Info, Check for updates,",
      "rem Install update) at $(stamp). It runs once and stays as a record.",
      "setlocal",
      "title Sparlectra update",
      "set \"LOG=$(q(log))\"",
      # paths in quotes: an & or a parenthesis in a folder name would end the echo
      "echo Sparlectra update from v$(current) to $(tag). Log: \"%LOG%\"",
      "echo Waiting for the Web UI (process $(pid)) to stop ...",
      "set /a n=0",
      ":wait_webui",
      "tasklist /FI \"PID eq $(pid)\" /NH 2>nul | find \" $(pid) \" >nul",
      "if errorlevel 1 goto webui_stopped",
      "set /a n+=1",
      "if %n% gtr $(WEBUI_UPDATE_WAIT_SECONDS) goto webui_still_running",
      "ping -n 2 127.0.0.1 >nul",
      "goto wait_webui",
      ":webui_still_running",
      "echo The Web UI did not stop within $(WEBUI_UPDATE_WAIT_SECONDS ÷ 60) minutes, nothing was changed.>>\"%LOG%\"",
      "echo The Web UI did not stop within $(WEBUI_UPDATE_WAIT_SECONDS ÷ 60) minutes, nothing was changed.",
      "pause",
      "exit /b 1",
      ":webui_stopped",
      "rem Windows does not move a folder a program still uses, for example the",
      "rem console of the old Web UI (its working directory) until it is closed",
      ":probe_folder",
      "move \"$(q(dir))\" \"$(q(dir)).update-probe\" >nul 2>nul",
      "if not errorlevel 1 goto folder_free",
      "echo The folder \"$(q(dir))\" is still in use. Close the window of the old Web UI, the update then goes on. Closing this window cancels the update.",
      "ping -n 4 127.0.0.1 >nul",
      "goto probe_folder",
      ":folder_free",
      "move \"$(q(dir)).update-probe\" \"$(q(dir))\" >nul",
      "rem set by the sysimage relaunch of the old Web UI; the new start has to check the image again",
      "set \"SPARLECTRA_SYSIMAGE_CHECKED=\"",
      "set \"SPARLECTRA_INSTALL_DIR=$(q(install_dir))\"",
      "set \"SPARLECTRA_UPDATE=1\"",
      "set \"SPARLECTRA_UPDATE_TAG=$(q(tag))\"",
      "set \"SPARLECTRA_BUILD_SYSIMAGE=$(b(answers.build_sysimage))\"",
      "set \"SPARLECTRA_CREATE_SHORTCUT=$(b(answers.create_shortcut))\"",
      "set \"SPARLECTRA_START_WEBUI=0\"",
      "echo Running the install script, this takes a few minutes. Its output goes to \"%LOG%\" ...",
      "call \"$(q(installer))\" >>\"%LOG%\" 2>&1 <nul",
      "set \"RC=%ERRORLEVEL%\"",
      "type \"%LOG%\"",
      "if not \"%RC%\"==\"0\" goto update_failed",
      ">\"$(q(status))\" echo ok",
      ">>\"$(q(status))\" echo $(q(tag))",
      ">>\"$(q(status))\" echo.",
      ">>\"$(q(status))\" echo %DATE% %TIME%",
      "echo Update to $(q(tag)) finished.",
      "goto update_done",
      ":update_failed",
      "rem the reason the install script printed, else its exit code; the reason keeps its leading blank",
      "set \"REASON= the install script ended with exit code %RC%\"",
      "for /f \"tokens=1,* delims=:\" %%A in ('findstr /b /c:\"Update FAILED:\" \"%LOG%\"') do set \"REASON=%%B\"",
      ">\"$(q(status))\" echo failed",
      ">>\"$(q(status))\" echo $(q(tag))",
      ">>\"$(q(status))\" echo%REASON%",
      ">>\"$(q(status))\" echo %DATE% %TIME%",
      "echo.",
      "echo Update FAILED:%REASON%",
      "echo The previous version starts again. Log: \"%LOG%\"",
      "pause",
      ":update_done",
      "set \"SPARLECTRA_INSTALL_DIR=\"",
      "set \"SPARLECTRA_UPDATE=\"",
      "set \"SPARLECTRA_UPDATE_TAG=\"",
      "set \"SPARLECTRA_BUILD_SYSIMAGE=\"",
      "set \"SPARLECTRA_CREATE_SHORTCUT=\"",
      "set \"SPARLECTRA_START_WEBUI=\"",
      "echo Starting the Web UI ...",
      "rem this console has a keyboard; the restart after an update does not ask again",
      "set \"SPARLECTRA_UPDATE_CHECK=0\"",
      "call \"$(q(joinpath(dir, "start_webui.bat")))\"",
    ]
    # CRLF: cmd.exe finds labels in an LF-only batch file unreliably
    write(updater, join(lines, "\r\n") * "\r\n")
    comspec = get(ENV, "COMSPEC", "cmd.exe")
    # `start` opens the console window of the update and returns at once;
    # verbatim, because cmd.exe parses this command line itself and the
    # quoting Julia applies for other programs would break the title argument
    command = Cmd(Cmd([comspec, "/c", "start", "\"Sparlectra update\"", "\"$(updater)\""]); windows_verbatim = true, dir = update_dir)
  else
    updater = joinpath(update_dir, "sparlectra_update_$(file_stamp).sh")
    qs = _webui_sh_quote
    script = """
    #!/bin/sh
    # Sparlectra update, written by the Web UI (Info, Check for updates,
    # Install update) at $(stamp). It runs once and stays as a record.
    exec >>$(qs(log)) 2>&1 </dev/null
    echo "Waiting for the Web UI (process $(pid)) to stop ..."
    n=0
    while kill -0 $(pid) 2>/dev/null
    do
      n=\$((n + 1))
      if [ \$n -gt $(WEBUI_UPDATE_WAIT_SECONDS) ]; then
        echo "The Web UI did not stop within $(WEBUI_UPDATE_WAIT_SECONDS ÷ 60) minutes, nothing was changed."
        exit 1
      fi
      sleep 1
    done
    # set by the sysimage relaunch of the old Web UI; the new start has to
    # check the image again
    unset SPARLECTRA_SYSIMAGE_CHECKED
    SPARLECTRA_INSTALL_DIR=$(qs(install_dir))
    SPARLECTRA_UPDATE=1
    SPARLECTRA_UPDATE_TAG=$(qs(tag))
    SPARLECTRA_BUILD_SYSIMAGE=$(b(answers.build_sysimage))
    SPARLECTRA_CREATE_SHORTCUT=$(b(answers.create_shortcut))
    SPARLECTRA_START_WEBUI=0
    export SPARLECTRA_INSTALL_DIR SPARLECTRA_UPDATE SPARLECTRA_UPDATE_TAG SPARLECTRA_BUILD_SYSIMAGE SPARLECTRA_CREATE_SHORTCUT SPARLECTRA_START_WEBUI
    sh $(qs(installer))
    rc=\$?
    if [ \$rc -eq 0 ]; then
      printf 'ok\\n%s\\n\\n%s\\n' $(qs(tag)) "\$(date '+%Y-%m-%d %H:%M')" > $(qs(status))
      echo "Update to $(tag) finished."
    else
      # the reason the install script printed, else its exit code
      reason=\$(grep '^Update FAILED: ' $(qs(log)) | tail -n 1 | sed 's/^Update FAILED: //')
      [ -n "\$reason" ] || reason="the install script ended with exit code \$rc"
      printf 'failed\\n%s\\n%s\\n%s\\n' $(qs(tag)) "\$reason" "\$(date '+%Y-%m-%d %H:%M')" > $(qs(status))
      echo "Update FAILED: \$reason. The previous version starts again."
    fi
    unset SPARLECTRA_INSTALL_DIR SPARLECTRA_UPDATE SPARLECTRA_UPDATE_TAG SPARLECTRA_BUILD_SYSIMAGE SPARLECTRA_CREATE_SHORTCUT SPARLECTRA_START_WEBUI
    echo "Starting the Web UI: $(joinpath(dir, "start_webui.sh"))"
    exec sh $(qs(joinpath(dir, "start_webui.sh")))
    """
    write(updater, script)
    command = Cmd(`sh $(updater)`; dir = update_dir)
  end
  return (command = command, updater = updater, installer = installer, log = log, install_dir = install_dir, root = dir, answers = answers)
end

"""
    start_webui_update!(plan; spawn) -> NamedTuple

Start `plan.command` detached, with no streams of its own (the updater
writes its log itself), and return `(started, message)`. `spawn` is the
function that runs the command; tests replace it.
"""
function start_webui_update!(plan; spawn = cmd -> run(pipeline(detach(cmd); stdout = devnull, stderr = devnull); wait = false))
  try
    spawn(plan.command)
  catch err
    # expected failure: the shell (sh, cmd.exe) cannot be started
    err isa Base.IOError || rethrow()
    return (started = false, message = "Could not start the update: $(sprint(showerror, err))")
  end
  return (started = true, message = "")
end

# The question of webui_console_update_offer: one line from `input`, or ""
# when nobody answers within `seconds`, at the end of the input, or when the
# input is closed while the question waits. Only an explicit y or yes counts.
function _webui_console_read_answer(prompt::AbstractString; io::IO = stdout, input::IO = stdin, seconds::Real = WEBUI_CONSOLE_UPDATE_ANSWER_SECONDS)::String
  print(io, prompt)
  flush(io)
  reader = @async try
    readline(input)
  catch err
    # expected failure: the input is closed while the question waits
    err isa EOFError || err isa Base.IOError || rethrow()
    ""
  end
  if timedwait(() -> istaskdone(reader), Float64(seconds); pollint = 0.1) !== :ok
    println(io)
    return ""
  end
  return fetch(reader)
end

"""
    webui_console_update_offer(; args, env, interactive, fetch, current, timeout, layout, root, ask, port_busy, sysimage_busy, operation_log, answers, spawn, io, windows) -> Symbol

The update check of `start_webui.jl`, before the server starts. Runs
`webui_check_for_update` with a short `timeout` (5 s) and, for a
newer release on an installed copy (`webui_install_layout`), asks on the
console whether to install it. Only an answer `y` or `yes` (any case)
installs, through the same updater as the Info panel (`webui_update_plan`,
`start_webui_update!`); the updater waits for this process to end, runs the
install script and starts the new version. Any other answer, an empty line,
the end of the input or no answer within 60 s means no.

Returns what happened, so the caller knows whether to start the server:

  * `:started`: the updater runs; the caller must exit without starting
    the server.
  * `:disabled`: `--no-update-check` in `args`, `SPARLECTRA_UPDATE_CHECK=0`
    or `JULIA_PKG_OFFLINE=true` in `env`; no request.
  * `:not_interactive`: stdin is not a terminal (`interactive`), no request
    and never an install.
  * `:up_to_date`, `:failed` (offline, timeout, rate limit, one line printed),
    `:manual` (git checkout or unknown layout: one line on how to update),
    `:busy` (a Web UI already runs on port 8080 or a sysimage build runs),
    `:declined`, `:start_failed` (the updater could not be written or
    started; nothing was changed).

Never throws for a network problem. `fetch`, `ask`, `port_busy`,
`sysimage_busy` and `spawn` are replaced in tests.
"""
function webui_console_update_offer(; args::AbstractVector{<:AbstractString} = ARGS, env::AbstractDict = ENV, interactive::Bool = isa(stdin, Base.TTY), fetch = _WEBUI_UPDATE_FETCH[], current::VersionNumber = version(), timeout::Real = WEBUI_CONSOLE_UPDATE_TIMEOUT_SECONDS, layout::Union{Nothing,Symbol} = nothing, root::AbstractString = SPARLECTRA_ROOT, ask = _webui_console_read_answer, port_busy = () -> _webui_port_holds_sparlectra("127.0.0.1", 8080), sysimage_busy = () -> sysimage_build_active(), operation_log::AbstractString = default_webui_operation_log_path(), answers = nothing, spawn = nothing, io::IO = stdout, windows::Bool = Sys.iswindows())::Symbol
  # a failed last update is named at every start, also when the check is off
  # (the restart after an update switches it off)
  failed_update = _webui_last_update_failure_text(operation_log)
  isempty(failed_update) || println(io, failed_update)
  opted_out = WEBUI_CONSOLE_UPDATE_FLAG in args || get(env, "SPARLECTRA_UPDATE_CHECK", "1") == "0" || lowercase(get(env, "JULIA_PKG_OFFLINE", "")) == "true"
  opted_out && return :disabled
  # piped input, CI or a service: no request and no question, so nothing is
  # ever installed without a person. The restart after an update does not
  # ask either: the unix updater starts it with stdin on /dev/null, the
  # Windows updater sets SPARLECTRA_UPDATE_CHECK=0 for its console.
  interactive || return :not_interactive
  check = webui_check_for_update(; fetch, current, timeout)
  if check.state === :failed
    println(io, "Update check: $(check.message).")
    return :failed
  end
  check.state === :up_to_date && return :up_to_date
  println(io, "Sparlectra $(check.latest) is available (running $(current)): $(check.release_url)")
  kind = layout === nothing ? webui_install_layout(root) : layout
  dir = rstrip(String(root), ['/', '\\'])
  if kind === :git_checkout
    println(io, "This copy is a git checkout and is updated by hand: git fetch --tags, then git checkout $(check.tag) in $(dir).")
    return :manual
  elseif kind !== :installed_copy
    println(io, "This copy was not made by the install script; update it with: $(windows ? WEBUI_UPDATE_ONE_LINE_WINDOWS : WEBUI_UPDATE_ONE_LINE_SH)")
    return :manual
  end
  # the updater waits for this process only; a second Web UI on the same
  # folder or a running sysimage build would have the folder replaced
  # under it
  if port_busy()
    println(io, "A Sparlectra Web UI already runs on port 8080; stop it to install the update.")
    return :busy
  end
  if sysimage_busy()
    println(io, "A sysimage build is running; install the update when it has finished.")
    return :busy
  end
  reply = something(ask("Install $(check.latest) now? [y/N] "), "")
  if !(lowercase(strip(reply)) in ("y", "yes"))
    println(io, "Not installed; starting $(current).")
    return :declined
  end
  plan = try
    webui_update_plan(; root = dir, current, tag = check.tag, operation_log, answers = answers === nothing ? _webui_update_answers() : answers, windows, requested_in = "from start_webui.jl")
  catch err
    # expected failure: the update folder next to the operation log cannot
    # be created or written
    err isa Base.IOError || err isa SystemError || rethrow()
    println(io, "Could not prepare the update ($(sprint(showerror, err))); starting $(current).")
    return :start_failed
  end
  started = spawn === nothing ? start_webui_update!(plan) : start_webui_update!(plan; spawn)
  record_webui_operation!(operation_log, "update_install_requested"; route = "start_webui.jl", method = "CONSOLE", user_action = true, status = started.started ? "started" : "failed", layout = String(kind), tag = check.tag, log = plan.log, message = started.started ? nothing : started.message)
  if !started.started
    println(io, "$(started.message); starting $(current).")
    return :start_failed
  end
  println(io, "Installing $(check.tag). The previous copy is kept as $(joinpath(plan.install_dir, "Sparlectra.old")); log: $(plan.log)")
  if windows
    println(io, "The update runs in a console window of its own and starts the new version. Close this window so that it can replace the folder.")
  else
    println(io, "The update runs in the background and starts the new version when it has finished.")
  end
  return :started
end
"""
    render_webui_update_check(check; layout, root, operation_log, active_job, sysimage_busy, started_by_script, windows) -> String

The page after "Check for updates": the result line, and for a newer
release what "Install update" does on this copy (installed copy), the git
commands (git checkout), or the one-line install (unknown layout).
"""
function render_webui_update_check(check; layout::Symbol, root::AbstractString = SPARLECTRA_ROOT, operation_log::AbstractString, active_job::Bool = false, sysimage_busy::Bool = false, started_by_script::Bool = _webui_started_by_start_script(), answers = nothing, windows::Bool = Sys.iswindows())::String
  e = _webui_escape
  dir = rstrip(String(root), ['/', '\\'])
  parts = String[]
  failed_update = _webui_last_update_failure_text(operation_log)
  isempty(failed_update) || push!(parts, "<div class=\"alert alert-error update-last-failed\" role=\"alert\">$(e(failed_update))</div>")
  if check.state === :failed
    push!(parts, "<div class=\"alert alert-error update-check-result\" role=\"alert\"><span>Update check failed: $(e(check.message)). Running $(e(string(check.current))). <a href=\"$(e(WEBUI_UPDATE_RELEASES_PAGE))\" target=\"_blank\" rel=\"noopener noreferrer\">Releases on GitHub</a></span></div>")
  elseif check.state === :up_to_date
    push!(parts, "<div class=\"alert alert-info update-check-result\" role=\"status\">Up to date ($(e(string(check.current)))).</div>")
  else
    push!(parts, "<div class=\"alert alert-info update-check-result\" role=\"status\"><strong>$(e(string(check.latest))) available</strong> (running $(e(string(check.current)))). <a href=\"$(e(check.release_url))\" target=\"_blank\" rel=\"noopener noreferrer\">Release notes</a></div>")
    push!(parts, _webui_update_install_section(check; layout, dir, operation_log, active_job, sysimage_busy, started_by_script, answers, windows))
  end
  push!(parts, "<p class=\"field-help\">The check asks GitHub for its latest-release page once per click, and start_webui.jl asks once per start from a terminal (off with --no-update-check or SPARLECTRA_UPDATE_CHECK=0); nothing about this installation is sent.</p>")
  return _webui_layout("Check for updates", "<section class=\"panel update-check\">$(join(parts, ""))</section>"; show_back = true)
end

function _webui_update_install_section(check; layout::Symbol, dir::AbstractString, operation_log::AbstractString, active_job::Bool, sysimage_busy::Bool, started_by_script::Bool, answers, windows::Bool)::String
  e = _webui_escape
  one_line = windows ? WEBUI_UPDATE_ONE_LINE_WINDOWS : WEBUI_UPDATE_ONE_LINE_SH
  if layout === :git_checkout
    cd_line = windows ? "cd /d \"$(dir)\"" : "cd \"$(dir)\""
    commands = join([cd_line, "git fetch --tags", "git checkout $(check.tag)"], "\n")
    return "<h2>Update by hand</h2><p>This copy is a git checkout (<code>$(e(dir))</code>). A development checkout is updated by hand; the Web UI does not change it. Stop the Web UI, then:</p><pre class=\"update-commands\">$(e(commands))</pre><p>Start the Web UI again afterwards. Local changes stay as they are; <code>git checkout</code> refuses a tag that would overwrite them.</p>"
  end
  if layout !== :installed_copy || !started_by_script
    where = layout === :installed_copy ? dirname(dir) : ""
    reason = layout === :installed_copy ?
             "This Web UI runs inside a Julia session that does not end when the Web UI stops, so it cannot update itself. Stop it and run the install command in <code>$(e(where))</code>, the folder that holds <code>Sparlectra</code>:" :
             "This copy was not installed by the install script, so the Web UI does not update it. The one-line install puts the latest release into <code>Sparlectra/</code> in the current folder, or updates the copy there:"
    return "<h2>Update by hand</h2><p>$(reason)</p><pre class=\"update-commands\">$(e(one_line))</pre>"
  end
  install_dir = dirname(dir)
  ans = answers === nothing ? _webui_update_answers() : answers
  log = webui_update_log_path(operation_log)
  script = joinpath("tools", _webui_update_installer_name(; windows))
  facts = string(
    "<dl class=\"sysimage-facts\">",
    "<dt>Copy</dt><dd><code>$(e(dir))</code></dd>",
    "<dt>Kept as</dt><dd><code>$(e(joinpath(install_dir, "Sparlectra.old")))</code></dd>",
    "<dt>Sysimage</dt><dd>", ans.build_sysimage ? "rebuilt after the update (this copy uses one)" : "not built (this copy has none)", "</dd>",
    "<dt>Desktop launcher</dt><dd>", ans.create_shortcut ? "written again (it exists)" : "not created (there is none)", "</dd>",
    "<dt>Log</dt><dd><code>$(e(log))</code></dd>",
    "</dl>",
  )
  how = windows ?
        "The update runs in a console window of its own. If it reports the folder as in use, close the window of the old Web UI." :
        "The update runs in the background without a terminal; the log named below shows its progress, and the new Web UI writes its output there too."
  head = "<h2>Install update</h2><p><strong>Install update</strong> stops this Web UI, runs <code>$(e(script))</code> of this copy unattended and starts the new version; the browser window opens again when it is ready. $(how)</p>$(facts)"
  active_job && return head * "<div class=\"alert warning\" role=\"status\">A run is active. The update can be installed when it has finished or was aborted.</div>"
  sysimage_busy && return head * "<div class=\"alert warning\" role=\"status\">A sysimage build is running. The update can be installed when it has finished.</div>"
  form = string(
    "<form method=\"post\" action=\"/webui/update/install\" class=\"update-install-form\">",
    "<input type=\"hidden\" name=\"token\" value=\"$(e(_webui_update_token()))\">",
    "<input type=\"hidden\" name=\"tag\" value=\"$(e(check.tag))\">",
    "<label class=\"check\"><input type=\"checkbox\" name=\"confirm_update\" value=\"yes\" required> Stop the Web UI and install $(e(string(check.latest)))</label> ",
    "<button type=\"submit\">Install update</button>",
    "</form>",
  )
  return head * form
end

"""
    handle_webui_update_check(; operation_log, fetch, layout, ...) -> SparlectraWebUIResponse

POST /webui/update/check: run the check once and render its page. The
result goes into the operation log as `update_check`.
"""
function handle_webui_update_check(; operation_log::AbstractString, fetch = _WEBUI_UPDATE_FETCH[], layout::Symbol = webui_install_layout(), root::AbstractString = SPARLECTRA_ROOT, current::VersionNumber = version(), output_root::AbstractString = default_webui_output_root(), started_by_script::Bool = _webui_started_by_start_script(), answers = nothing)::SparlectraWebUIResponse
  check = webui_check_for_update(; fetch, current)
  record_webui_operation!(operation_log, "update_check"; route = "/webui/update/check", method = "POST", user_action = true, status = String(check.state), latest = check.latest === nothing ? nothing : string(check.latest), failure = check.state === :failed ? String(check.failure) : nothing, message = isempty(check.message) ? nothing : check.message, layout = String(layout))
  page = render_webui_update_check(check; layout, root, operation_log, active_job = _webui_active_job() !== nothing, sysimage_busy = sysimage_build_active(; output_root), started_by_script, answers)
  return _webui_html(page)
end

"""
    handle_webui_update_install(form; operation_log, shutdown, spawn, ...) -> SparlectraWebUIResponse

POST /webui/update/install. Refuses (and starts nothing) unless the copy is
an installed copy started by start_webui.jl, the confirm box is ticked, the
form carries this process's token and a release tag, no run is active and
no sysimage build runs. Then it writes the plan (`webui_update_plan`),
starts the updater detached and stops this Web UI through `shutdown`
shortly after the answer is sent. A git checkout or an unknown layout never
starts a process.
"""
function handle_webui_update_install(form::AbstractDict; operation_log::AbstractString, shutdown, spawn = nothing, layout::Symbol = webui_install_layout(), root::AbstractString = SPARLECTRA_ROOT, output_root::AbstractString = default_webui_output_root(), started_by_script::Bool = _webui_started_by_start_script(), active_job = _webui_active_job, answers = nothing, windows::Bool = Sys.iswindows())::SparlectraWebUIResponse
  field(name) = strip(String(something(_webui_form_value(form, name, ""), "")))
  function refuse(status::Integer, message::AbstractString)
    record_webui_operation!(operation_log, "update_install_requested"; route = "/webui/update/install", method = "POST", user_action = true, status = "rejected", layout = String(layout), message = String(message))
    return _webui_html(render_webui_error(status, message); status)
  end
  layout === :installed_copy || return refuse(409, "This copy is not an installation made by the install script; the Web UI does not update it. Check for updates shows how to update it.")
  started_by_script || return refuse(409, "This Web UI runs inside a Julia session that does not end when it stops; stop it and run the install script by hand.")
  field("confirm_update") == "yes" || return refuse(400, "Tick the confirmation box to install the update.")
  field("token") == _webui_update_token() || return refuse(403, "The install form is out of date. Run Check for updates again.")
  tag = field("tag")
  occursin(r"^v\d+\.\d+\.\d+$", tag) || return refuse(400, "The install form names no release.")
  active_job() === nothing || return refuse(409, "A run is active. Install the update when it has finished or was aborted.")
  sysimage_build_active(; output_root) && return refuse(409, "A sysimage build is running. Install the update when it has finished.")
  shutdown === nothing && return refuse(503, "The update can only be installed from a running Web UI server.")
  plan = webui_update_plan(; root, tag, operation_log, answers = answers === nothing ? _webui_update_answers() : answers, windows)
  started = spawn === nothing ? start_webui_update!(plan) : start_webui_update!(plan; spawn)
  if !started.started
    record_webui_operation!(operation_log, "update_install_requested"; route = "/webui/update/install", method = "POST", user_action = true, status = "failed", layout = String(layout), message = started.message, log = plan.log)
    return _webui_html(render_webui_error(500, started.message); status = 500)
  end
  record_webui_operation!(operation_log, "update_install_requested"; route = "/webui/update/install", method = "POST", user_action = true, status = "started", layout = String(layout), tag, log = plan.log, updater = plan.updater)
  @async begin
    # let the answer reach the browser before the listener closes
    sleep(0.2)
    shutdown()
  end
  return _webui_html(render_webui_update_started(plan; tag, windows))
end

# A page of its own, without the header: the server stops right after it.
function render_webui_update_started(plan; tag::AbstractString, windows::Bool = Sys.iswindows())::String
  e = _webui_escape
  how = windows ?
        "The update runs in its own console window and starts the new version when it is done. If that window reports the folder as in use, close the window of the old Web UI." :
        "The update runs in the background and starts the new version when it is done; the browser window opens again by itself."
  return "<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width,initial-scale=1\"><title>Installing update · Sparlectra</title><link rel=\"stylesheet\" href=\"/static/sparlectra.css\"></head><body><main><section class=\"panel shutdown-message\"><h1>Installing $(e(tag))</h1><p>This Web UI stops now. $(e(how))</p><p>Log: <code>$(e(plan.log))</code></p><p>The previous copy is kept at <code>$(e(joinpath(plan.install_dir, "Sparlectra.old")))</code>; delete it once the new version runs.</p><p>If no new window opens, read the log and start the Web UI with the desktop launcher or <code>$(e(joinpath(plan.root, windows ? "start_webui.bat" : "start_webui.sh")))</code>.</p></section></main></body></html>"
end
