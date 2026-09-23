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

# --- sysimage handling ------------------------------------------------------
# This runs BEFORE `using Sparlectra`, deliberately: the launcher module sits
# on plain Base, so deciding whether to relaunch through the sysimage costs
# no package load. Julia can only take an image at process start (-J), so a
# process that wants one has to start itself again.
include(joinpath(@__DIR__, "tools", "sysimage_launcher.jl"))
using .SysimageLauncher: handle_sysimage, prepare_environments, repair_environment, ENV_ONLY_FLAG

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

function main()
  # No warm-up: either this process runs on the sysimage, where the code is
  # already compiled, or the user chose to compile on first use and was told
  # so. A hidden warm-up run on top of that only delayed the first page.
  server = SparlectraApp.start_sparlectra_webui(open_browser = true)
  # nothing = a Sparlectra Web UI already runs on the port; it was opened in
  # the browser instead, so there is no new server task to wait on.
  server === nothing && return nothing

  try
    wait(server.task)
  catch err
    if err isa InterruptException
      @info "Ctrl+C received; closing Sparlectra Web UI."
      close(server)
      return nothing
    end
    rethrow()
  end
end

# VS Code's inline eval runs the whole file as one world age. Since Julia
# 1.12.7 even a direct getfield counts as a binding access from the old
# world, so the lookup itself must also go through invokelatest.
Base.invokelatest(Base.invokelatest(getfield, @__MODULE__, :main))
