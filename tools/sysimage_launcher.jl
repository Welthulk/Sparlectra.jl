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
# file: tools/sysimage_launcher.jl
# purpose: the decisions every Sparlectra launcher has to make BEFORE the
#          package can be loaded: is this checkout's environment resolvable
#          at all, and is there a usable sysimage to relaunch through. Deliberately a module on plain Base (TOML and SHA are
#          stdlib): start_webui.jl includes it BEFORE `using Sparlectra`,
#          because the whole point is to decide whether this process should
#          be replaced by one that starts from the image. Being a module also
#          makes the decision testable without starting a Web UI.

module SysimageLauncher

using TOML
using SHA

export handle_sysimage, unresolved_dependencies, repair_environment

const REBUILD_FLAG = "--rebuild-sysimage"
const NO_IMAGE_FLAG = "--no-sysimage"

"Path of the Web UI sysimage for this platform (mirrors default_webui_output_root)."
function sysimage_path()::String
  if Sys.isapple()
    return joinpath(homedir(), "Library", "Application Support", "Sparlectra", "WebUI", "sysimage", "sparlectra.dylib")
  elseif Sys.iswindows()
    return joinpath(get(ENV, "LOCALAPPDATA", joinpath(homedir(), "AppData", "Local")), "Sparlectra", "WebUI", "sysimage", "sparlectra.dll")
  end
  state = get(ENV, "XDG_STATE_HOME", joinpath(homedir(), ".local", "state"))
  return joinpath(state, "sparlectra", "webui", "sysimage", "sparlectra.so")
end

"""
Return `nothing` when the image at `path` can be used, otherwise a short
reason naming what is wrong. Three things are checked, and the third one is
the expensive lesson: the meta file pins the Julia version and the Manifest
hash, i.e. the DEPENDENCIES, and says nothing about this package's own code.
Editing `src/` leaves the Manifest untouched, so an image built before the
edit still looked fresh and the Web UI silently served the old code. A
source file newer than the image therefore counts as outdated.

`Sparlectra.webui_sysimage_problem` makes the same decision for the Web UI's
sysimage page and cannot be called from here: this function decides whether
the current process should be REPLACED by one booting from the image, which
is settled before the package is loaded. `test/test_webui.jl` asserts that
both implementations agree, so the duplication cannot drift unnoticed.
"""
function sysimage_problem(path::AbstractString, project_dir::AbstractString)::Union{Nothing,String}
  meta_path = joinpath(dirname(path), "sysimage_meta.toml")
  (isfile(path) && isfile(meta_path)) || return "no sysimage found"
  meta = try
    TOML.parsefile(meta_path)
  catch
    return "the sysimage metadata is unreadable"
  end
  built_for = get(meta, "julia_version", "")
  built_for == string(VERSION) || return "the sysimage was built for Julia $(built_for), this is $(VERSION)"
  manifest = joinpath(project_dir, "Manifest.toml")
  if isfile(manifest)
    current = bytes2hex(open(sha256, manifest))
    get(meta, "manifest_sha256", "") == current || return "the sysimage does not match the current Manifest.toml"
  end
  image_time = mtime(path)
  src = joinpath(project_dir, "src")
  if isdir(src)
    for (root, _, files) in walkdir(src), f in files
      endswith(f, ".jl") || continue
      mtime(joinpath(root, f)) > image_time && return "the sysimage is older than $(src)"
    end
  end
  return nothing
end

"""
Ask `prompt` and return true unless the user answers no. No answer means
yes: an unattended start (nobody at the keyboard, or no terminal at all)
should end up with the image rather than quietly without it.
"""
function ask_build(prompt::AbstractString; seconds::Real = 30.0)::Bool
  isa(stdin, Base.TTY) || return true
  print(prompt)
  flush(stdout)
  answer = Ref{String}("")
  reader = @async begin
    answer[] = try
      readline(stdin)
    catch
      ""
    end
  end
  waited = 0.0
  while !istaskdone(reader) && waited < seconds
    sleep(0.1)
    waited += 0.1
  end
  istaskdone(reader) || println()
  return !(lowercase(strip(answer[])) in ("n", "no", "nein"))
end

function build_sysimage(project_dir::AbstractString)::Bool
  exe = joinpath(Sys.BINDIR, Base.julia_exename())
  # No banner here: the build script prints its own header and then keeps ONE
  # progress line up to date. Anything printed around it would scroll that
  # line away, which is exactly what the quiet console is for.
  # The build script directly, not through `using Sparlectra; buildSysimage()`:
  # that intermediate process would load (and possibly precompile) the whole
  # package just to spawn the same script. A checkout without tools/ (an
  # installation from the registry) falls back to the package entry point.
  script = joinpath(project_dir, "tools", "build_sysimage.jl")
  cmd = isfile(script) ? `$(exe) --startup-file=no --project=$(project_dir) $(script)` :
        `$(exe) --startup-file=no --project=$(project_dir) -e "using Sparlectra; buildSysimage()"`
  try
    run(cmd)
    return true
  catch err
    @warn "The sysimage build failed; starting without it" exception = err
    return false
  end
end

"""
    relaunch_with_sysimage(image, script, project_dir, args)

Start `script` again in a new process that boots from `image`, then exit
with that process's result.

Threads need care. They are fixed at process start, so the child has to be
told. JULIA_NUM_THREADS is inherited and may say `auto` or carry an
interactive pool (`auto,1`), which `Threads.nthreads()` cannot express: it
only counts the default pool. So the variable wins where it exists, and
only a caller who passed `-t` on the command line gets an explicit
`--threads`.
"""
function relaunch_with_sysimage(image::AbstractString, script::AbstractString, project_dir::AbstractString, args::Vector{String})
  exe = joinpath(Sys.BINDIR, Base.julia_exename())
  ENV["SPARLECTRA_SYSIMAGE_CHECKED"] = "1"
  thread_flag = String[]
  haskey(ENV, "JULIA_NUM_THREADS") || push!(thread_flag, "--threads=$(Threads.nthreads())")
  # NO startup file in the child, unless the caller insists. The child is the
  # process that RUNS on the image, and a personal startup.jl usually loads
  # Revise: measured, loading Revise invalidates 1530 precompiled method
  # instances of this image, which then have to be inferred again on first
  # use. Building an image and then throwing a third of it away on startup is
  # the "still slow with the sysimage" report of 2026-09-07.
  startup = get(ENV, "SPARLECTRA_STARTUP_FILE", "no")
  cmd = `$(exe) -J$(image) --startup-file=$(startup) --project=$(project_dir) $(thread_flag) $(script) $(args)`
  println("Sysimage: ", image)
  # Ctrl+C reaches the whole process group, so it hits this process too. The
  # child prints its own "closing Sparlectra Web UI"; a stack trace from the
  # parent on top of that would look like a crash. ignorestatus keeps a
  # non-zero child exit from raising here as well - it is forwarded instead.
  code = try
    run(ignorestatus(cmd)).exitcode
  catch err
    err isa InterruptException || rethrow()
    0
  end
  exit(code)
end

"""
    unresolved_dependencies(project_dir) -> Vector{String}

Direct dependencies of `Project.toml` that `Manifest.toml` does not know, or
`["<no Manifest.toml>"]` when there is no manifest at all. Empty means the
environment can be loaded.

Two TOML reads, no package load, a few milliseconds. That cheapness is the
point: it runs BEFORE the sysimage question, and the order matters. Answering
"no" to that question and then walking into an unloadable environment is what
happened on a Windows 11 checkout (2026-09-07); answering "yes" would have been
worse, because the build works in the SEPARATE @sparlectra-sysimage-build
environment and would have spent eleven minutes before the checkout's own
manifest turned out to be the problem.
"""
function unresolved_dependencies(project_dir::AbstractString)::Vector{String}
  project = joinpath(project_dir, "Project.toml")
  manifest = joinpath(project_dir, "Manifest.toml")
  isfile(project) || return String[]
  isfile(manifest) || return ["<no Manifest.toml>"]
  deps, known = try
    d = get(TOML.parsefile(project), "deps", Dict{String,Any}())
    m = TOML.parsefile(manifest)
    # manifest_format 2.0 nests everything under [deps]; 1.0 puts the packages
    # at the top level next to the metadata keys
    n = haskey(m, "deps") ? keys(m["deps"]) : setdiff(keys(m), ("julia_version", "manifest_format", "project_hash", "manifest_version"))
    (keys(d), n)
  catch err
    # expected failure: a truncated or hand-edited TOML. Saying "unresolved"
    # here is right - the environment cannot be trusted either way - and the
    # repair below reports what Pkg makes of it.
    println("Could not read the package environment (", first(sprint(showerror, err), 120), "); trying to repair it.")
    return ["<unreadable Project.toml or Manifest.toml>"]
  end
  return sort!([String(d) for d in deps if !(d in known)])
end

"""
Bring the environment into a loadable state, printing only when something is
actually wrong. `resolve` before `instantiate`, in that order: instantiate
installs what the manifest lists and cannot add a package the manifest never
mentioned.
"""
function repair_environment(project_dir::AbstractString, reason::AbstractString)
  println(reason)
  println("Resolving the package environment; this happens once after a checkout or a dependency change.")
  @eval using Pkg
  pkgm = Base.invokelatest(getfield, @__MODULE__, :Pkg)
  try
    Base.invokelatest(pkgm.resolve)
    Base.invokelatest(pkgm.instantiate)
    println("Package environment resolved.")
  catch err
    println()
    println("Could not prepare the dependencies of this checkout.")
    println("Run this once in the checkout directory and start again:")
    println("    julia --project=. -e \"using Pkg; Pkg.resolve(); Pkg.instantiate()\"")
    println()
    println("If that fails too, delete Manifest.toml and repeat. It is not tracked,")
    println("and a manifest left over from an older Sparlectra is the usual reason.")
    println()
    rethrow(err)
  end
  return nothing
end

"""
    handle_sysimage(args, script, project_dir)

The whole decision, shared by every platform: skip when asked to
(`SPARLECTRA_NO_SYSIMAGE=1`, `--no-sysimage`) or when this process already
booted from the image; relaunch through a usable image; otherwise offer the
build and relaunch through the fresh one. `--rebuild-sysimage` builds even
when the current image is fine. Returns `nothing` when the caller should
just continue in this process.
"""
function handle_sysimage(args::Vector{String}, script::AbstractString, project_dir::AbstractString)
  get(ENV, "SPARLECTRA_NO_SYSIMAGE", "0") == "1" && return nothing
  NO_IMAGE_FLAG in args && return nothing
  # set by relaunch_with_sysimage: this process IS the relaunched one
  get(ENV, "SPARLECTRA_SYSIMAGE_CHECKED", "0") == "1" && return nothing

  image = sysimage_path()
  rebuild = REBUILD_FLAG in args
  passthrough = filter(a -> a != REBUILD_FLAG && a != NO_IMAGE_FLAG, args)
  problem = sysimage_problem(image, project_dir)

  if !rebuild && problem === nothing
    relaunch_with_sysimage(image, script, project_dir, passthrough)
    return nothing
  end

  if rebuild
    println("Rebuilding the sysimage on request ($(REBUILD_FLAG)).")
  else
    println("Sysimage: ", problem, ".")
    println("Building one takes a few minutes. Afterwards the Web UI starts in seconds")
    println("and no page has to compile on first use.")
    if !ask_build("Build the sysimage now? [Y/n] ")
      println("Starting without a sysimage: Julia compiles every code path on first")
      println("use, so the first click on each page takes noticeably longer.")
      println("You can build one later with:  $(basename(script)) $(REBUILD_FLAG)")
      return nothing
    end
  end

  if build_sysimage(project_dir)
    problem = sysimage_problem(image, project_dir)
    if problem === nothing
      relaunch_with_sysimage(image, script, project_dir, passthrough)
      return nothing
    end
    println("The build produced no usable image ($(problem)); starting without it.")
  end
  return nothing
end

end # module
