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

export handle_sysimage, unresolved_dependencies, compat_lower_bound, outdated_dependencies, repair_environment, prepare_environments, sysimage_source_roots, sysimage_build_script, external_sysimage, remove_stale_sysimage, ENV_ONLY_FLAG

const REBUILD_FLAG = "--rebuild-sysimage"
const NO_IMAGE_FLAG = "--no-sysimage"
# prepare the environments and compile, then leave: what a fresh checkout
# needs before the first real start, and what a test can run headless
const ENV_ONLY_FLAG = "--env-only"

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
    sysimage_source_roots(project_dir) -> Vector{String}

The source trees an image built for `project_dir` depends on: the
project's own `src`, and, when the project is the application package
inside a library checkout (`app/` next to the library's `Project.toml`),
the library's `src` as well. An edit in either tree outdates the image.
"""
function sysimage_source_roots(project_dir::AbstractString)::Vector{String}
  roots = [joinpath(project_dir, "src")]
  parent = dirname(abspath(project_dir))
  isfile(joinpath(parent, "Project.toml")) && push!(roots, joinpath(parent, "src"))
  return roots
end

# the build script lives in the repository's tools directory: next to the
# library project, one level above the application project
function sysimage_build_script(project_dir::AbstractString)::String
  own = joinpath(project_dir, "tools", "build_sysimage.jl")
  isfile(own) && return own
  return joinpath(dirname(abspath(project_dir)), "tools", "build_sysimage.jl")
end

# the metadata file that belongs to the managed image, next to it
sysimage_meta_path(image::AbstractString)::String = joinpath(dirname(image), "sysimage_meta.toml")

"""
    external_sysimage(managed; current = <image of this process>) -> Union{Nothing,String}

The image this process was started on when it is neither the stock Julia
image nor the managed one: an image given from outside with `-J`. Such an
image is not checked, not rebuilt and never removed; `nothing` otherwise.
"""
function external_sysimage(managed::AbstractString; current::AbstractString = unsafe_string(Base.JLOptions().image_file))::Union{Nothing,String}
  isempty(current) && return nothing
  startswith(basename(current), "sys.") && return nothing
  same = try
    isfile(current) && isfile(managed) && realpath(current) == realpath(managed)
  catch
    false
  end
  return same ? nothing : String(current)
end

"""
    remove_stale_sysimage(image) -> Symbol

Delete the managed image and its metadata file, nothing else: the two known
paths, no wildcard. Returns `:removed`, `:absent` when neither exists, or
`:locked` when the file cannot be deleted (Windows keeps an image open
while another Julia runs on it); then a message says so, the start goes on
without an image and the next start tries again.
"""
function remove_stale_sysimage(image::AbstractString)::Symbol
  meta = sysimage_meta_path(image)
  (isfile(image) || isfile(meta)) || return :absent
  try
    isfile(image) && rm(image)
    isfile(meta) && rm(meta)
    return :removed
  catch err
    println("The outdated sysimage could not be removed (", first(sprint(showerror, err), 120), "); it stays in place and is not used. The next start tries again.")
    return :locked
  end
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
  for src in sysimage_source_roots(project_dir)
    isdir(src) || continue
    for (root, _, files) in walkdir(src), f in files
      endswith(f, ".jl") || continue
      mtime(joinpath(root, f)) > image_time && return "the sysimage is older than $(src)"
    end
  end
  return nothing
end

"""
Ask `prompt` and return true only when the user answers yes. The build is
opt-in: Enter, a timeout and a start without a terminal all mean no. A build
takes minutes, on a machine with a virus scanner on the compile cache far
longer, and a user who did not ask for it should not wait for it;
`--rebuild-sysimage` and a plain `y` start one.
"""
function ask_build(prompt::AbstractString; seconds::Real = 30.0)::Bool
  isa(stdin, Base.TTY) || return false
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
  return lowercase(strip(answer[])) in ("y", "yes", "j", "ja")
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
  script = sysimage_build_script(project_dir)
  cmd = isfile(script) ? `$(exe) --startup-file=no --project=$(project_dir) $(script)` :
        `$(exe) --startup-file=no --project=$(project_dir) -e "using SparlectraApp; buildSysimage()"`
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
whose manifest version lies below the lower bound of the `[compat]` entry
(reported as `"<name> <manifest version> < <bound>"`), or
`["<no Manifest.toml>"]` when there is no manifest at all. Empty means the
environment can be loaded.

The compat case matters because Julia does not check a Manifest against
the compat: a checkout whose manifest predates a dependency bump loads the
old version without a message, and an old solver version can compute wrong
numbers. Only the plain caret form of a compat entry (`"0.9.15"`,
`"^0.9.15"`, `"~0.9.15"`, `"=0.9.15"`, first entry of a comma list) is read;
ranges and inequalities need Pkg and are left to the resolve.

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
  deps, known, compat, entries = try
    pt = TOML.parsefile(project)
    d = get(pt, "deps", Dict{String,Any}())
    m = TOML.parsefile(manifest)
    # manifest_format 2.0 nests everything under [deps]; 1.0 puts the packages
    # at the top level next to the metadata keys
    e = haskey(m, "deps") ? m["deps"] : Dict{String,Any}(k => v for (k, v) in m if !(k in ("julia_version", "manifest_format", "project_hash", "manifest_version")))
    (keys(d), keys(e), get(pt, "compat", Dict{String,Any}()), e)
  catch err
    # expected failure: a truncated or hand-edited TOML. Saying "unresolved"
    # here is right - the environment cannot be trusted either way - and the
    # repair below reports what Pkg makes of it.
    println("Could not read the package environment (", first(sprint(showerror, err), 120), "); trying to repair it.")
    return ["<unreadable Project.toml or Manifest.toml>"]
  end
  missing = [String(d) for d in deps if !(d in known)]
  # a known dependency whose manifest version is below the compat bound
  for d in deps
    d in known || continue
    bound = compat_lower_bound(get(compat, d, ""))
    bound === nothing && continue
    entry = entries[d]
    # manifest 2.0 lists one dict per package name; stdlibs carry no version
    versions = [get(x, "version", "") for x in (entry isa AbstractVector ? entry : [entry]) if x isa AbstractDict]
    for v in versions
      isempty(v) && continue
      parsed = tryparse(VersionNumber, String(v))
      parsed === nothing && continue
      parsed < bound && push!(missing, string(d, " ", parsed, " < ", bound))
    end
  end
  return sort!(missing)
end

"""
    compat_lower_bound(spec) -> Union{Nothing,VersionNumber}

The lower bound of a `[compat]` entry in its plain caret form (`"0.9.15"`,
`"^0.9.15"`, `"~0.9.15"`, `"=0.9.15"`; a comma list contributes its first
entry). Missing components count as zero. Anything else (ranges, `>=`,
hyphens) returns `nothing`: Pkg's parser is the authority for those and it
is not loaded here on purpose.
"""
function compat_lower_bound(spec::AbstractString)::Union{Nothing,VersionNumber}
  first_entry = strip(first(split(spec, ','; limit = 2)))
  isempty(first_entry) && return nothing
  first_entry = lstrip(first_entry, ('^', '~', '='))
  occursin(r"^\d+(\.\d+){0,2}$", first_entry) || return nothing
  return VersionNumber(first_entry)
end

"""
Bring the environment into a loadable state, printing only when something is
actually wrong. `resolve` before `instantiate`, in that order: instantiate
installs what the manifest lists and cannot add a package the manifest never
mentioned.

`update` names packages whose manifest version lies below the compat bound
(the `"<name> <version> < <bound>"` entries of `unresolved_dependencies`).
They are updated FIRST: `resolve` treats every manifest version as an
explicit requirement and fails with "Unsatisfiable requirements" on such a
manifest instead of lifting the version; `Pkg.update` of the named packages
lifts it.
"""
function repair_environment(project_dir::AbstractString, reason::AbstractString; update::Vector{String} = String[], label::AbstractString = "package")
  println(reason)
  println("Resolving and compiling the ", label, " environment; this happens once after a checkout or a dependency change and takes a while.")
  @eval using Pkg
  pkgm = Base.invokelatest(getfield, @__MODULE__, :Pkg)
  started = time()
  try
    # the environment to repair is the one named, not whichever project the
    # process was started with (the launcher runs from the library checkout
    # and repairs the application environment app/)
    Base.invokelatest(pkgm.activate, project_dir; io = devnull)
    isempty(update) || Base.invokelatest(pkgm.update, update)
    Base.invokelatest(pkgm.resolve)
    # instantiate installs and precompiles what the manifest lists
    Base.invokelatest(pkgm.instantiate)
    println("The ", label, " environment is ready after ", round(Int, time() - started), " s.")
  catch err
    println()
    println("Could not prepare the dependencies of this checkout.")
    println("Run this once in the checkout directory and start again:")
    if isempty(update)
      println("    julia --project=$(project_dir) -e \"using Pkg; Pkg.resolve(); Pkg.instantiate()\"")
    else
      println("    julia --project=$(project_dir) -e \"using Pkg; Pkg.update([" * join(("\\\"" * u * "\\\"" for u in update), ", ") * "]); Pkg.instantiate()\"")
    end
    println()
    println("If that fails too, delete Manifest.toml and repeat. It is not tracked,")
    println("and a manifest left over from an older Sparlectra is the usual reason.")
    println()
    rethrow(err)
  end
  return nothing
end

"""
    prepare_environments(library_dir, app_dir) -> (library = Bool, application = Bool)

Bring the library environment and the application environment into a
loadable state, in that order (the application carries the library by
path), and say what was done. Each flag is true when that environment was
resolved and compiled here. A second start finds both up to date and does
nothing beyond the two TOML reads.
"""
function prepare_environments(library_dir::AbstractString, app_dir::AbstractString)
  repaired = Bool[]
  for (label, dir) in (("library", library_dir), ("application", app_dir))
    missing_deps = unresolved_dependencies(dir)
    if isempty(missing_deps)
      push!(repaired, false)
      continue
    end
    println("First start: setting up the ", label, " environment (", dir, "): ", join(missing_deps, ", "), ".")
    repair_environment(dir, "This checkout has no ready " * label * " environment yet."; update = outdated_dependencies(missing_deps), label = label)
    push!(repaired, true)
  end
  any(repaired) || println("Environments up to date: library and application.")
  return (library = repaired[1], application = repaired[2])
end

"""
    outdated_dependencies(unresolved) -> Vector{String}

The package names among the entries of `unresolved_dependencies` that report
a manifest version below the compat bound; the argument for `update` of
`repair_environment`.
"""
function outdated_dependencies(unresolved::AbstractVector{<:AbstractString})::Vector{String}
  return [String(first(split(d, ' '))) for d in unresolved if occursin(" < ", d)]
end

# `image`, `ask`, `relaunch` and `current_image` are the seams a test uses:
# a fixture image, a fixed answer, a recorded relaunch instead of a new
# process, and the image this process runs on.
"""
    handle_sysimage(args, script, project_dir)

The whole decision, shared by every platform: skip when asked to
(`SPARLECTRA_NO_SYSIMAGE=1`, `--no-sysimage`) or when this process already
booted from the image; relaunch through a usable image; otherwise offer the
build and relaunch through the fresh one. `--rebuild-sysimage` builds even
when the current image is fine. Returns `nothing` when the caller should
just continue in this process.
"""
function handle_sysimage(args::Vector{String}, script::AbstractString, project_dir::AbstractString;
                         image::AbstractString = sysimage_path(), ask = ask_build, relaunch = relaunch_with_sysimage,
                         current_image::AbstractString = unsafe_string(Base.JLOptions().image_file))
  get(ENV, "SPARLECTRA_NO_SYSIMAGE", "0") == "1" && return nothing
  NO_IMAGE_FLAG in args && return nothing
  # set by relaunch_with_sysimage: this process IS the relaunched one
  get(ENV, "SPARLECTRA_SYSIMAGE_CHECKED", "0") == "1" && return nothing

  # an image given from outside (-J) is neither judged nor touched
  external = external_sysimage(image; current = current_image)
  if external !== nothing
    println("Running on a sysimage given from outside (", external, "): it is not checked, not rebuilt and not removed.")
    return nothing
  end

  rebuild = REBUILD_FLAG in args
  passthrough = filter(a -> a != REBUILD_FLAG && a != NO_IMAGE_FLAG, args)
  problem = sysimage_problem(image, project_dir)

  if !rebuild && problem === nothing
    relaunch(image, script, project_dir, passthrough)
    return nothing
  end

  # an image that exists but cannot be used (Manifest or a source tree
  # changed, other Julia, unreadable metadata, metadata missing) is stale;
  # the Web UI never starts on it, and unless a new one replaces it now, it
  # goes
  stale = problem !== nothing && (isfile(image) || isfile(sysimage_meta_path(image)))

  if rebuild
    println("Rebuilding the sysimage on request ($(REBUILD_FLAG)).")
  else
    println("Sysimage: ", problem, ".")
    println("Building one takes a few minutes. Afterwards the Web UI starts in seconds")
    println("and no page has to compile on first use.")
    if !ask("Build the sysimage now? [y/N] ")
      if stale && remove_stale_sysimage(image) === :removed
        println("Sysimage is out of date and was removed. Start with $(basename(script)) $(REBUILD_FLAG) to build a new one.")
      end
      println("Starting without a sysimage: Julia compiles every code path on first")
      println("use, so the first click on each page takes noticeably longer.")
      println("You can build one later with:  $(basename(script)) $(REBUILD_FLAG)")
      return nothing
    end
  end

  # the build writes to a staging file and replaces the image only when it
  # succeeded; a failed build leaves the old image, which is stale and goes
  if build_sysimage(project_dir)
    problem = sysimage_problem(image, project_dir)
    if problem === nothing
      relaunch(image, script, project_dir, passthrough)
      return nothing
    end
    println("The build produced no usable image ($(problem)); starting without it.")
  end
  if stale && remove_stale_sysimage(image) === :removed
    println("The outdated sysimage was removed; start with $(basename(script)) $(REBUILD_FLAG) to build a new one.")
  end
  return nothing
end

end # module
