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
#
# file: start_webui.jl
# purpose: start the Sparlectra Web UI, and own the sysimage decision for
#          every platform. The start scripts only call this file, so Linux,
#          macOS and Windows go through exactly the same code: check for a
#          usable sysimage (tools/sysimage_launcher.jl), offer to build one
#          when it is missing or outdated, then relaunch through the image.
#          Started from a terminal, it then asks GitHub for a newer release
#          and offers to install it (an installed copy only); only an answer
#          y or yes installs. `--no-update-check` or SPARLECTRA_UPDATE_CHECK=0
#          skips the check for unattended starts, and so does
#          JULIA_PKG_OFFLINE=true.

# --- sysimage handling ------------------------------------------------------
# This runs BEFORE `using Sparlectra`, deliberately: the launcher module sits
# on plain Base, so deciding whether to relaunch through the sysimage costs
# no package load. Julia can only take an image at process start (-J), so a
# process that wants one has to start itself again.
include(joinpath(@__DIR__, "tools", "sysimage_launcher.jl"))
using .SysimageLauncher: handle_sysimage, console_answer, prepare_environments, repair_environment, running_on_sysimage, warm_webui_page, ENV_ONLY_FLAG

const _SPARLECTRA_LIBRARY_DIR = abspath(@__DIR__)
const _SPARLECTRA_APP_DIR = abspath(joinpath(@__DIR__, "app"))

# The environments FIRST, then the sysimage question. A fresh checkout (git
# clone or release download) has no Manifest.toml in either directory; both
# are resolved, installed and compiled here, and the user is told what is
# happening and how long it took. A dependency below its compat bound (a
# manifest that predates a version bump) is repaired the same way: it is the
# one state that loads without complaint and can compute wrong numbers.
# `--env-only` stops after this block and the compile below: the first-start
# work without a server, for a test or for an installer.
Base.invokelatest(prepare_environments, _SPARLECTRA_LIBRARY_DIR, _SPARLECTRA_APP_DIR)
const _ENV_ONLY = ENV_ONLY_FLAG in ARGS

_ENV_ONLY || handle_sysimage(filter(!=(ENV_ONLY_FLAG), copy(ARGS)), @__FILE__, _SPARLECTRA_APP_DIR)

# --- the application package ------------------------------------------------
# It lives in app/ with its own environment; putting that directory on the
# load path makes `using SparlectraApp` resolve from any starting project
# without loading Pkg first. Loading compiles what the environment setup
# above did not (a machine where the compile cache was wiped), so the first
# start says so and reports the time; a second start finds everything
# compiled. The repair below is the second line of defence: a depot with a
# missing artifact or a half-written package directory gets here instead of
# failing deep inside Base.
pushfirst!(LOAD_PATH, _SPARLECTRA_APP_DIR)
let app_id = Base.PkgId(Base.UUID("2e39659a-4f4c-4326-b34e-e3df865cea74"), "SparlectraApp"),
    lib_id = Base.PkgId(Base.UUID("31ce9bba-fd9d-44a1-b005-f5f509afda38"), "Sparlectra")
  compiled = Base.isprecompiled(lib_id) && Base.isprecompiled(app_id)
  compiled || println("First start: compiling the library and the application; this takes a while.")
  started = time()
  try
    @eval using SparlectraApp
  catch err
    Base.invokelatest(repair_environment, _SPARLECTRA_APP_DIR,
      "Loading SparlectraApp failed (" * first(sprint(showerror, err), 160) * ")."; label = "application")
    @eval using SparlectraApp
  end
  compiled ? println("Packages already compiled.") : println("Compiled in ", round(Int, time() - started), " s.")
end

if _ENV_ONLY
  println("Environment check finished.")
  exit(0)
end

# --- update check -----------------------------------------------------------
# After the sysimage decision, because the check lives in SparlectraApp: on a
# sysimage it runs in the relaunched process, whose id the updater waits for.
# Only an installed copy is offered the install, and only an explicit yes
# installs. When the updater runs, this process ends without a server; the
# updater starts the new version. The question reads the console through
# console_answer like the sysimage question, so a reader left over from that
# question after its timeout takes the answer to this one instead of racing
# with a second reader for the same line.
function offer_update()
  ask = prompt -> console_answer(prompt; seconds = SparlectraApp.WEBUI_CONSOLE_UPDATE_ANSWER_SECONDS)
  SparlectraApp.webui_console_update_offer(; ask) === :started && exit(0)
  return nothing
end
Base.invokelatest(Base.invokelatest(getfield, @__MODULE__, :offer_update))

function main()
  # On a sysimage the pages are compiled already and the window opens at
  # once. Without one, the package image of the application carries the
  # compiled start page (app/src/precompile.jl); the launcher starts the
  # server without a window, renders the start page once through the
  # server's own routing call (same runtime, same argument types as a GET
  # request, so the window's request finds everything ready) and opens the
  # window when the page is ready, a few seconds later.
  on_sysimage = running_on_sysimage()
  server = SparlectraApp.start_sparlectra_webui(open_browser = on_sysimage)
  # nothing = a Sparlectra Web UI already runs on the port; on a sysimage it
  # was opened in the browser already, without one it is opened here (that
  # instance compiled its pages long ago); no new server task to wait on.
  if server === nothing
    on_sysimage || SparlectraApp._webui_open_browser("http://127.0.0.1:8080/powerflow")
    return nothing
  end
  # A script run ends the process on SIGINT by default: Julia prints the
  # signal and a native backtrace, and no handler below ever runs. With this
  # switch Ctrl+C arrives as an InterruptException in this task, and the
  # server stops through the same path as Stop Web UI (reason :ctrl_c).
  Base.exit_on_sigint(false)
  try
    if !on_sysimage
      # the request line split into SubStrings like the server reads it, so
      # the call meets the code the package image compiled for a real GET
      warm_webui_page(() -> SparlectraApp.route_sparlectra_webui(split("GET /powerflow")..., Dict{String,String}(); output_root = abspath(SparlectraApp.default_webui_output_root()), runtime = server.runtime).status == 200)
      SparlectraApp._webui_open_browser(server.url)
    end
    SparlectraApp._wait_sparlectra_webui(server)
  catch err
    # Ctrl+C during the warmup or the browser start, before the wait
    err isa InterruptException || rethrow()
    SparlectraApp._webui_request_shutdown!(server.runtime; reason = :ctrl_c)
  end
  # The accept loop aborts running jobs and prints "Sparlectra Web UI
  # stopped by Ctrl+C." when it ends; the process waits for that line. A
  # second Ctrl+C while running jobs are aborted ends the wait at once.
  try
    wait(server.task)
  catch err
    err isa InterruptException || rethrow()
  end
  return nothing
end

# VS Code's inline eval runs the whole file as one world age. Since Julia
# 1.12.7 even a direct getfield counts as a binding access from the old
# world, so the lookup itself must also go through invokelatest.
Base.invokelatest(Base.invokelatest(getfield, @__MODULE__, :main))
