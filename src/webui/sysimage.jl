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

# file: src/webui/sysimage.jl
# purpose: shared sysimage location and metadata contract: path helpers
#          below the Web UI user root, the sysimage_meta.toml writer, and
#          how a session reports the way it was started (app, sysimage,
#          native) for the Web UI header. The validity check itself lives in
#          start_webui.jl, which must decide BEFORE loading this package
#          whether to relaunch through the image. PackageCompiler never
#          appears here; building happens in tools/build_sysimage.jl in a
#          separate environment and process.

"""Return the sysimage directory beneath the Web UI user application root."""
webui_sysimage_dir(output_root::AbstractString = default_webui_output_root())::String = joinpath(dirname(abspath(output_root)), "sysimage")

"""Return the platform sysimage file extension: so (Linux), dylib (macOS), dll (Windows)."""
webui_sysimage_ext()::String = Sys.iswindows() ? "dll" : (Sys.isapple() ? "dylib" : "so")

"""Return the full sysimage path `<user root>/sysimage/sparlectra.<ext>`."""
webui_sysimage_path(output_root::AbstractString = default_webui_output_root())::String = joinpath(webui_sysimage_dir(output_root), string("sparlectra.", webui_sysimage_ext()))

"""Return the metadata path next to the sysimage (`sysimage_meta.toml`)."""
webui_sysimage_meta_path(output_root::AbstractString = default_webui_output_root())::String = joinpath(webui_sysimage_dir(output_root), "sysimage_meta.toml")

"""
Return the Manifest.toml whose state the sysimage is built against: the
checkout manifest when running from a dev checkout (the manifest the
launchers run with), otherwise the manifest of the active environment
(registered installations ship no Manifest.toml in the package directory).
"""
function webui_sysimage_manifest_path()::String
  checkout = joinpath(SPARLECTRA_ROOT, "Manifest.toml")
  isfile(checkout) && return checkout
  proj = Base.active_project()
  proj === nothing && return checkout
  return joinpath(dirname(proj), "Manifest.toml")
end

"""SHA-256 hex digest of a file, the manifest fingerprint of the metadata contract."""
_sysimage_file_sha256(path::AbstractString)::String = bytes2hex(sha256(read(path)))

# short OS label; informational next to the hard julia-version/manifest checks
_sysimage_os_label()::String = Sys.iswindows() ? "windows" : (Sys.isapple() ? "macos" : (Sys.islinux() ? "linux" : "other"))

"""
    write_sysimage_meta(meta_path; manifest_path) -> Dict

Write `sysimage_meta.toml`, the validity contract of the image:
Sparlectra version, Julia `VERSION` string, OS label and architecture,
SHA-256 and path of `manifest_path` (best effort, omitted when no manifest
exists, e.g. registered installations), and the build timestamp. Returns the
written dictionary. The launchers compare the Julia version string and the
manifest hash with shell tools only, so both are stored as plain strings.
"""
function write_sysimage_meta(meta_path::AbstractString; manifest_path::AbstractString = webui_sysimage_manifest_path())
  meta = Dict{String,Any}("sparlectra_version" => string(SparlectraVersion), "julia_version" => string(VERSION), "os" => _sysimage_os_label(), "arch" => string(Sys.ARCH), "built_at" => Dates.format(Dates.now(), "yyyy-mm-ddTHH:MM:SS"))
  # The manifest fingerprint is best effort: registered installations may
  # legitimately lack a resolvable Manifest.toml. The sparlectra_version
  # check in sysimage_status covers the package-update case there; a missing
  # manifest must never fail the build after a 10-plus-minute compile.
  if isfile(manifest_path)
    meta["manifest_sha256"] = _sysimage_file_sha256(manifest_path)
    meta["manifest_path"] = abspath(manifest_path)
  else
    @warn "write_sysimage_meta: no Manifest.toml found; staleness detection falls back to the Sparlectra/Julia version checks" manifest_path
  end
  mkpath(dirname(abspath(meta_path)))
  open(meta_path, "w") do io
    TOML.print(io, meta)
  end
  return meta
end

"""
    webui_runtime_flavor() -> NamedTuple

How this Julia session was started, for the Web UI header: `kind` is
`:app` (standalone executable; the app entry point stamps
`SPARLECTRA_APP_BUILT` with its build time), `:sysimage` (the running
system image IS the Sparlectra image), or `:native` (plain session).
`built` carries the build timestamp for app and sysimage, `nothing` for
native. Detection never throws; unreadable metadata falls back to the
image file time.
"""
function webui_runtime_flavor()
  app_built = strip(get(ENV, "SPARLECTRA_APP_BUILT", ""))
  isempty(app_built) || return (kind = :app, built = String(app_built))
  img = try
    unsafe_string(Base.JLOptions().image_file)
  catch
    ""
  end
  simg = try
    webui_sysimage_path()
  catch
    ""
  end
  # realpath on both sides: the flatpak state root is a symlink onto
  # ~/.local/state/sparlectra, and JLOptions carries the resolved path
  same = false
  try
    same = !isempty(img) && !isempty(simg) && isfile(img) && isfile(simg) && realpath(img) == realpath(simg)
  catch
    same = false
  end
  if same
    built = nothing
    try
      meta = webui_sysimage_meta_path()
      if isfile(meta)
        m = match(r"built_at\s*=\s*\"([^\"]+)\"", read(meta, String))
        m !== nothing && (built = String(m.captures[1]))
      end
    catch
    end
    built === nothing && (built = Libc.strftime("%Y-%m-%dT%H:%M:%S", mtime(img)))
    return (kind = :sysimage, built = built)
  end
  return (kind = :native, built = nothing)
end

# --- refreshing the image from the Web UI ------------------------------------
# A session can work for hours on an image that was correct when it was built
# and still pay JIT, simply because a path was never traced. "Refresh the
# image" is therefore a normal user action, not a maintenance chore, and it
# has to be reachable without leaving the browser.

"""Return the build log the sysimage build writes next to the image."""
webui_sysimage_build_log_path(output_root::AbstractString = default_webui_output_root())::String = joinpath(webui_sysimage_dir(output_root), "sysimage_build.log")

"""Return the progress file a running sysimage build updates while it works."""
webui_sysimage_progress_path(output_root::AbstractString = default_webui_output_root())::String = joinpath(webui_sysimage_dir(output_root), "sysimage_build_progress.toml")

"""
    webui_sysimage_problem(; image_path, project_dir) -> Union{Nothing,String}

Return `nothing` when the image on disk can be used, otherwise a short reason
naming what is wrong: wrong Julia version, dependencies changed
(Manifest.toml), or a source file newer than the image.

The third check is the expensive lesson: the metadata pins the Julia version
and the Manifest hash, i.e. the DEPENDENCIES, and says nothing about this
package's own code. Editing `src/` leaves the Manifest untouched, so an image
built before the edit still looked fresh and the Web UI silently served the
old code.

!!! note "Kept in sync with the launcher by a test"
    `SysimageLauncher.sysimage_problem` in `tools/sysimage_launcher.jl` makes
    the same decision and cannot call this one: it runs BEFORE the package is
    loaded, because what it decides is whether this process should be
    replaced by one that boots from the image. The two implementations are
    held together by `test/test_webui.jl`, which asserts they agree.
"""
function webui_sysimage_problem(; image_path::AbstractString = webui_sysimage_path(), project_dir::AbstractString = SPARLECTRA_ROOT)::Union{Nothing,String}
  meta_path = joinpath(dirname(image_path), "sysimage_meta.toml")
  (isfile(image_path) && isfile(meta_path)) || return "no sysimage found"
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
  image_time = mtime(image_path)
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
    read_sysimage_build_progress(; output_root) -> Union{Nothing,Dict{String,Any}}

Read the progress file `tools/build_sysimage.jl` keeps while it builds, or
`nothing` when no build has ever run here. Keys: `state` (`running`, `done`,
`failed`), `step`/`steps`, `phase`, `detail`, `message`, `elapsed_seconds`,
`updated_at`, `pid`, `log`.
"""
function read_sysimage_build_progress(; output_root::AbstractString = default_webui_output_root())::Union{Nothing,Dict{String,Any}}
  path = webui_sysimage_progress_path(output_root)
  isfile(path) || return nothing
  return try
    TOML.parsefile(path)
  catch err
    # expected failure: the builder is writing the file right now, so a read
    # can land on a partial line. Reporting "no progress yet" for one poll is
    # correct; the next refresh a second later gets the complete file.
    @debug "read_sysimage_build_progress: progress file not parseable yet" path exception = err
    nothing
  end
end

"""Seconds since a running build last reported progress; `nothing` without a build."""
function _sysimage_build_age(progress)::Union{Nothing,Float64}
  progress === nothing && return nothing
  updated = get(progress, "updated_at", nothing)
  updated isa Number || return nothing
  return max(0.0, time() - Float64(updated))
end

# The builder refreshes the progress file about every five seconds even while
# nothing changes, so a build that stopped reporting for a minute is gone
# (killed, crashed hard, machine rebooted) rather than merely busy.
const _SYSIMAGE_BUILD_STALE_SECONDS = 60.0

# The `starting` state has no heartbeat yet: it is written before the build
# process exists, and the entry stands until that process has booted Julia and
# written its own first entry. On a cold Windows start that takes tens of
# seconds, so this window is deliberately much wider than the running one.
const _SYSIMAGE_BUILD_STARTING_SECONDS = 180.0

"""
    sysimage_build_active(; output_root) -> Bool

Whether a sysimage build is running right now, which includes the seconds
between pressing the button and the build process reporting for the first
time. A build that stopped updating its progress file counts as gone, not as
active: a crashed builder must never lock the refresh button forever.
"""
function sysimage_build_active(; output_root::AbstractString = default_webui_output_root())::Bool
  progress = read_sysimage_build_progress(; output_root)
  progress === nothing && return false
  state = get(progress, "state", "")
  state in ("running", "starting") || return false
  age = _sysimage_build_age(progress)
  age === nothing && return true
  return age < (state == "starting" ? _SYSIMAGE_BUILD_STARTING_SECONDS : _SYSIMAGE_BUILD_STALE_SECONDS)
end

"""
    _write_sysimage_build_progress(output_root; state, phase, message)

Write a minimal progress entry from the Web UI side.

This exists for the gap between pressing the refresh button and the build
process being alive: the page only polls while a build is active, and until
2026-09-07 "active" required the child's first progress entry. Julia needs
seconds to boot before it can write one, tens of seconds on a cold Windows,
and during that window the page sat there looking dead (reported by the
maintainer, who asked whether there is a status display at all).

The state written here is `starting`, not `running`, and that distinction is
load-bearing: the build script refuses to start when it finds a FRESH
`running` entry, so writing `running` here would make the Web UI block the
very build it just launched.
"""
function _write_sysimage_build_progress(output_root::AbstractString; state::AbstractString, phase::AbstractString, message::AbstractString = "")
  esc(v) = replace(String(v), "\\" => "\\\\", "\"" => "\\\"", "\n" => " ")
  path = webui_sysimage_progress_path(output_root)
  try
    mkpath(dirname(path))
    open(path, "w") do io
      println(io, "state = \"", esc(state), "\"")
      println(io, "step = 0")
      println(io, "steps = 4")
      println(io, "phase = \"", esc(phase), "\"")
      println(io, "detail = \"\"")
      println(io, "message = \"", esc(message), "\"")
      println(io, "elapsed_seconds = 0.0")
      println(io, "updated_at = ", round(time(); digits = 1))
      println(io, "pid = 0")
      println(io, "log = \"", esc(webui_sysimage_build_log_path(output_root)), "\"")
    end
  catch err
    # expected failure: an unwritable or full state directory. The build still
    # runs, only its progress cannot be shown, so this is reported and not
    # turned into a failed rebuild.
    @warn "Could not write the sysimage build progress file; the build runs, but the page cannot show its progress." path exception = err
  end
  return nothing
end

"""
    start_sysimage_rebuild!(; output_root) -> NamedTuple

Start `tools/build_sysimage.jl` in a DETACHED child process and return
immediately: `(started, message, log)`. The build takes minutes, so it must
outlive both the request and the Web UI session that triggered it.

Two details are load-bearing:

  * The child is started from the plain Julia executable, never through
    `Base.julia_cmd()`. That command carries the `-J` of the CURRENT session,
    and when the Web UI itself runs on the Sparlectra image, PackageCompiler
    would build the new image incrementally on top of the old one.
  * The build writes to a staging file and moves it into place at the end, so
    this very session can keep running on the image being replaced (`mv`
    swaps the directory entry and leaves the old inode mapped).
"""
function start_sysimage_rebuild!(; output_root::AbstractString = default_webui_output_root())
  log = webui_sysimage_build_log_path(output_root)
  sysimage_build_active(; output_root) && return (started = false, message = "A sysimage build is already running.", log = log)
  pkgroot = pkgdir(@__MODULE__)
  pkgroot === nothing && return (started = false, message = "Cannot locate the Sparlectra package directory.", log = log)
  script = joinpath(pkgroot, "tools", "build_sysimage.jl")
  isfile(script) || return (started = false, message = "The build script is missing at $(script); a registry installation without tools/ cannot rebuild the image.", log = log)
  exe = joinpath(Sys.BINDIR, Base.julia_exename())
  cmd = Cmd(`$(exe) --startup-file=no --project=$(pkgroot) $(script)`; dir = pkgroot)
  # BEFORE spawning: the page polls only while a build is active, and the child
  # needs seconds to boot before it can report anything. Writing the entry
  # first also means a spawn that throws leaves a truthful `failed` entry
  # rather than a phantom build.
  _write_sysimage_build_progress(output_root; state = "starting", phase = "starting the build process")
  try
    # Progress and log go to files, so the child needs no streams of its own.
    run(pipeline(detach(cmd); stdout = devnull, stderr = devnull); wait = false)
  catch err
    reason = "Could not start the build: $(sprint(showerror, err))"
    _write_sysimage_build_progress(output_root; state = "failed", phase = "starting the build process", message = reason)
    return (started = false, message = reason, log = log)
  end
  return (started = true, message = "The sysimage build has started. It runs in the background and takes a few minutes.", log = log)
end

# Which image file this SESSION is running on, remembered the first time it is
# asked. A refresh replaces the image by `mv`, so the running process keeps its
# old mapped inode while the path already points at the new file: after a
# successful rebuild the Web UI still serves the code it started with, and it
# has to say so instead of letting the user wonder why the JIT pause is back.
const _SESSION_IMAGE_MTIME = Ref(0.0)
const _SESSION_IMAGE_SEEN = Ref(false)

function _session_image_mtime()::Float64
  if !_SESSION_IMAGE_SEEN[]
    img = try
      unsafe_string(Base.JLOptions().image_file)
    catch
      ""
    end
    _SESSION_IMAGE_MTIME[] = isempty(img) || !isfile(img) ? 0.0 : mtime(img)
    _SESSION_IMAGE_SEEN[] = true
  end
  return _SESSION_IMAGE_MTIME[]
end

"""
    sysimage_restart_pending(; output_root) -> Bool

True when the image on disk is newer than the one this session booted from,
i.e. a rebuild finished and only a restart of the Web UI will actually use it.
False for a native session, which has no image to be behind.
"""
function sysimage_restart_pending(; output_root::AbstractString = default_webui_output_root())::Bool
  flavor = try
    webui_runtime_flavor()
  catch
    (kind = :native, built = nothing)
  end
  flavor.kind === :sysimage || return false
  stamp = _session_image_mtime()
  stamp == 0.0 && return false
  img = webui_sysimage_path(output_root)
  # one second of slack: mtime resolution differs per filesystem
  return isfile(img) && mtime(img) > stamp + 1.0
end
