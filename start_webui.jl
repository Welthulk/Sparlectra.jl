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
using .SysimageLauncher: handle_sysimage

handle_sysimage(copy(ARGS), @__FILE__, abspath(@__DIR__))

# --- package environment -----------------------------------------------------
# `using Sparlectra` needs a RESOLVED environment, and Manifest.toml is not
# tracked. Two situations lead here, and they need different remedies:
#
#   * a fresh checkout has no manifest at all;
#   * an OLD manifest is present and predates a dependency change, so it never
#     learned about a package the Project.toml now requires.
#
# Julia names neither in a way a user can act on. The dependency scan fails
# deep inside Base with
#
#   KeyError: key Base.PkgId(UUID("19ecf91d-..."), "AnalyticLoadFlow") not found
#
# and a stacktrace through Base.Precompilation.
#
# RESOLVE, then instantiate, and in that order. `instantiate` installs what the
# manifest lists and cannot add a package the manifest never mentioned: on a
# stale manifest it fails with "AnalyticLoadFlow is a direct dependency, but
# does not appear in the manifest ... run Pkg.resolve()". This block used to
# call instantiate alone, so it turned one unhelpful error into another
# (reported from a Windows 11 checkout, 2026-09-07). `resolve` rewrites the
# manifest from Project.toml and covers the missing-manifest case as well.
try
  @eval using Sparlectra
catch err
  println("Preparing the package environment; this happens once after a checkout or a dependency change.")
  println("(", first(sprint(showerror, err), 160), ")")
  @eval using Pkg
  pkgm = Base.invokelatest(getfield, @__MODULE__, :Pkg)
  try
    Base.invokelatest(pkgm.resolve)
    Base.invokelatest(pkgm.instantiate)
  catch resolve_err
    println()
    println("Could not prepare the dependencies of this checkout.")
    println("Run this once in the checkout directory and start again:")
    println("    julia --project=. -e \"using Pkg; Pkg.resolve(); Pkg.instantiate()\"")
    println()
    println("If that fails too, delete Manifest.toml and repeat. It is not tracked,")
    println("and a manifest left over from an older Sparlectra is the usual reason.")
    println()
    rethrow(resolve_err)
  end
  @eval using Sparlectra
end

function main()
  # No warm-up: either this process runs on the sysimage, where the code is
  # already compiled, or the user chose to compile on first use and was told
  # so. A hidden warm-up run on top of that only delayed the first page.
  server = Sparlectra.start_sparlectra_webui(open_browser = true)
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
