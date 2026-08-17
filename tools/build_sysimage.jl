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

# file: tools/build_sysimage.jl
# purpose: build the optional fast-start sysimage for the Web UI with
#          PackageCompiler in the shared build environment
#          @sparlectra-sysimage-build (the package Project.toml stays free
#          of PackageCompiler). The image lands below the Web UI user root
#          next to its sysimage_meta.toml validity contract; the launchers
#          pick it up automatically. Run from the checkout:
#          `julia tools/build_sysimage.jl` (a `--dry-run` flag prints the
#          plan without building; the test suite uses it as a parse smoke).

const _REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const _BUILD_ENV = "sparlectra-sysimage-build"
const _DRY_RUN = ("--dry-run" in ARGS) || get(ENV, "SPARLECTRA_SYSIMAGE_DRY_RUN", "0") == "1"

# Pkg is deliberately NOT loaded at the top level: the dry-run path (used by
# the test-suite smoke check) must work in environments whose load path has
# no Pkg, such as the Pkg.test sandbox. The real build path loads Pkg inside
# main() and resolves the module dynamically (Julia 1.12 world age).

function _log(msg::AbstractString)
  println("[build_sysimage] ", msg)
  flush(stdout)
end

function _prepare_build_env(pkgm::Module)::Bool
  _log("activating shared build environment @$(_BUILD_ENV)")
  Base.invokelatest(pkgm.activate, _BUILD_ENV; shared = true)
  if Base.find_package("PackageCompiler") === nothing
    _log("installing PackageCompiler into the build environment (one-time)")
    Base.invokelatest(pkgm.add, "PackageCompiler")
  end
  # match how the script was invoked: from a checkout the repo itself is
  # developed into the build environment, otherwise the released package
  project_file = joinpath(_REPO_ROOT, "Project.toml")
  in_checkout = isfile(project_file) && occursin("name = \"Sparlectra\"", read(project_file, String))
  if in_checkout
    _log("developing Sparlectra from $(_REPO_ROOT)")
    Base.invokelatest(pkgm.develop; path = _REPO_ROOT)
  else
    _log("adding the released Sparlectra package")
    Base.invokelatest(pkgm.add, "Sparlectra")
  end
  # When the user has AnalyticLoadFlow installed (APSLF solver), it must be
  # part of the image: start_webui.jl loads it at startup, and loading a
  # package AFTER the sysimage invalidates precompiled methods inside the
  # image — the first run then silently recompiles them (measured 36 s
  # instead of 1 s on case118). find_package sees the user's default
  # environment through the load path, so this detects a plain
  # `Pkg.add("AnalyticLoadFlow")` installation.
  with_alf = Base.find_package("AnalyticLoadFlow") !== nothing
  if with_alf
    _log("AnalyticLoadFlow detected; baking it into the image as well")
    try
      Base.invokelatest(pkgm.add, "AnalyticLoadFlow")
    catch err
      _log("could not add AnalyticLoadFlow to the build environment ($(sprint(showerror, err))); building without it — the first run after start will recompile invalidated methods")
      with_alf = false
    end
  end
  Base.invokelatest(pkgm.instantiate)
  return with_alf
end

function main()
  started = time()
  if _DRY_RUN
    # parse-and-plan mode for the test suite: resolve everything that does
    # not touch the package environment or PackageCompiler
    _log("dry run: no build environment changes, no PackageCompiler call")
    if Base.find_package("Sparlectra") !== nothing
      @eval using Sparlectra
      # the module binding is newer than this function invocation (Julia
      # 1.12 world age); resolve it dynamically via invokelatest
      spar = Base.invokelatest(getfield, @__MODULE__, :Sparlectra)
      img = Base.invokelatest(getglobal(spar, :webui_sysimage_path))
      meta = Base.invokelatest(getglobal(spar, :webui_sysimage_meta_path))
      _log("would build: $(img)")
      _log("would write: $(meta)")
    else
      _log("Sparlectra not loadable in this environment; dry run ends after parsing")
    end
    _log("dry run finished")
    return nothing
  end

  @eval using Pkg
  pkgm = Base.invokelatest(getfield, @__MODULE__, :Pkg)
  with_alf = _prepare_build_env(pkgm)
  @eval using PackageCompiler
  @eval using Sparlectra
  # same world-age pattern as above: both modules were loaded inside this
  # function invocation
  spar = Base.invokelatest(getfield, @__MODULE__, :Sparlectra)
  pc = Base.invokelatest(getfield, @__MODULE__, :PackageCompiler)

  img = Base.invokelatest(getglobal(spar, :webui_sysimage_path))
  meta_path = Base.invokelatest(getglobal(spar, :webui_sysimage_meta_path))
  manifest = Base.invokelatest(getglobal(spar, :webui_sysimage_manifest_path))
  mkpath(dirname(img))

  # pre-fetch the workload cases so the child process never needs the
  # network: case14.m for the power-flow service run, the MiniGrid CGMES
  # delivery for the short-circuit service run. A failed MiniGrid fetch
  # only costs that trace (the workload logs the gap), never the build.
  _log("ensuring case14.m is cached for the workload")
  Base.invokelatest(getglobal(spar, :ensure_casefile), "case14.m")
  case_cache = Base.invokelatest(getglobal(spar, :default_webui_case_cache_dir))
  if !isfile(joinpath(case_cache, "cgmes_minigrid.zip"))
    _log("fetching the MiniGrid CGMES delivery for the short-circuit workload")
    try
      cgmes = getglobal(spar, :CGMESImporter)
      Base.invokelatest(getglobal(cgmes, :fetchCGMESTestSet), "minigrid"; outdir = case_cache)
    catch err
      _log("MiniGrid fetch failed ($(sprint(showerror, err))); the short-circuit service path will not be traced")
    end
  end

  packages = with_alf ? ["Sparlectra", "AnalyticLoadFlow"] : ["Sparlectra"]
  workload = joinpath(@__DIR__, "sysimage_workload.jl")
  _log("building sysimage (this typically takes 10 to 20 minutes)...")
  _log("target: $(img)")
  _log("packages: $(join(packages, ", "))")
  Base.invokelatest(getglobal(pc, :create_sysimage), packages; sysimage_path = img, precompile_execution_file = workload)

  Base.invokelatest(getglobal(spar, :write_sysimage_meta), meta_path; manifest_path = manifest)
  elapsed = round((time() - started) / 60; digits = 1)
  size_mb = round(filesize(img) / 1024^2; digits = 1)
  _log("done in $(elapsed) min: $(img) ($(size_mb) MB)")
  _log("start_webui.sh / start_webui.bat use the image automatically on the next start")
  _log("bypass with SPARLECTRA_NO_SYSIMAGE=1; rebuild after package updates")
  return nothing
end

Base.invokelatest(main)
