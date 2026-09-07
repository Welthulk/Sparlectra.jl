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
# purpose: the sysimage decision every Sparlectra launcher shares: is there a
#          usable image, should one be built, and relaunching the process
#          through it. Deliberately a module on plain Base (TOML and SHA are
#          stdlib): start_webui.jl includes it BEFORE `using Sparlectra`,
#          because the whole point is to decide whether this process should
#          be replaced by one that starts from the image. Being a module also
#          makes the decision testable without starting a Web UI.

module SysimageLauncher

using TOML
using SHA

export handle_sysimage

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
  println("Building the sysimage (this window stays busy until it is done)...")
  # The build script directly, not through `using Sparlectra; buildSysimage()`:
  # that intermediate process would load (and possibly precompile) the whole
  # package just to spawn the same script. A checkout without tools/ (an
  # installation from the registry) falls back to the package entry point.
  script = joinpath(project_dir, "tools", "build_sysimage.jl")
  cmd = isfile(script) ? `$(exe) --project=$(project_dir) $(script)` :
        `$(exe) --project=$(project_dir) -e "using Sparlectra; buildSysimage()"`
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
  cmd = `$(exe) -J$(image) --project=$(project_dir) $(thread_flag) $(script) $(args)`
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
    println("Building one takes roughly 10 to 15 minutes. Afterwards the Web UI starts")
    println("in seconds and no page has to compile on first use.")
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
