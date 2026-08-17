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
# purpose: shared sysimage location and metadata contract for the fast-start
#          feature: path helpers below the Web UI user root, the
#          sysimage_meta.toml writer and reader, and the validity check that
#          the launchers mirror with shell tools (Julia version string and
#          Manifest.toml SHA-256). PackageCompiler itself never appears
#          here; building happens in tools/build_sysimage.jl in a separate
#          environment and process.

"""Return the sysimage directory beneath the Web UI user application root."""
webui_sysimage_dir(output_root::AbstractString = default_webui_output_root())::String = joinpath(dirname(abspath(output_root)), "sysimage")

"""Return the platform sysimage file extension: so (Linux), dylib (macOS), dll (Windows)."""
webui_sysimage_ext()::String = Sys.iswindows() ? "dll" : (Sys.isapple() ? "dylib" : "so")

"""Return the full sysimage path `<user root>/sysimage/sparlectra.<ext>`."""
webui_sysimage_path(output_root::AbstractString = default_webui_output_root())::String = joinpath(webui_sysimage_dir(output_root), string("sparlectra.", webui_sysimage_ext()))

"""Return the metadata path next to the sysimage (`sysimage_meta.toml`)."""
webui_sysimage_meta_path(output_root::AbstractString = default_webui_output_root())::String = joinpath(webui_sysimage_dir(output_root), "sysimage_meta.toml")

"""Return the Manifest.toml whose state the sysimage is built against (the checkout manifest the launchers run with)."""
webui_sysimage_manifest_path()::String = joinpath(SPARLECTRA_ROOT, "Manifest.toml")

"""SHA-256 hex digest of a file, the manifest fingerprint of the metadata contract."""
_sysimage_file_sha256(path::AbstractString)::String = bytes2hex(sha256(read(path)))

# short OS label; informational next to the hard julia-version/manifest checks
_sysimage_os_label()::String = Sys.iswindows() ? "windows" : (Sys.isapple() ? "macos" : (Sys.islinux() ? "linux" : "other"))

"""
    write_sysimage_meta(meta_path; manifest_path) -> Dict

Write `sysimage_meta.toml`, the validity contract of the fast-start image:
Sparlectra version, Julia `VERSION` string, OS label and architecture,
SHA-256 of `manifest_path`, and the build timestamp. Returns the written
dictionary. The launchers compare the Julia version string and the manifest
hash with shell tools only, so both are stored as plain strings.
"""
function write_sysimage_meta(meta_path::AbstractString; manifest_path::AbstractString = webui_sysimage_manifest_path())
  isfile(manifest_path) || error("write_sysimage_meta: manifest not found at $(manifest_path)")
  meta = Dict{String,Any}(
    "sparlectra_version" => string(SparlectraVersion),
    "julia_version" => string(VERSION),
    "os" => _sysimage_os_label(),
    "arch" => string(Sys.ARCH),
    "manifest_sha256" => _sysimage_file_sha256(manifest_path),
    "built_at" => Dates.format(Dates.now(), "yyyy-mm-ddTHH:MM:SS"),
  )
  mkpath(dirname(abspath(meta_path)))
  open(meta_path, "w") do io
    TOML.print(io, meta)
  end
  return meta
end

"""Read `sysimage_meta.toml`; returns the dictionary or `nothing` when missing or unparsable."""
function read_sysimage_meta(meta_path::AbstractString)::Union{Nothing,Dict{String,Any}}
  isfile(meta_path) || return nothing
  return try
    TOML.parsefile(meta_path)
  catch
    nothing
  end
end

"""
    sysimage_status(; output_root, manifest_path) -> NamedTuple

Evaluate the fast-start image state for the Web UI page and the tests:
`present` (image file exists), `meta_present`, `valid` (all checks pass),
`reasons` (English mismatch descriptions, empty when valid), plus size,
build timestamp, and the stored versions. The checks mirror the launcher
logic exactly: Julia version string equality and Manifest.toml SHA-256
equality, plus OS/arch equality that the launchers get for free by running
on the machine that built the image.
"""
function sysimage_status(; output_root::AbstractString = default_webui_output_root(), manifest_path::AbstractString = webui_sysimage_manifest_path())
  image_path = webui_sysimage_path(output_root)
  meta_path = webui_sysimage_meta_path(output_root)
  present = isfile(image_path)
  meta = read_sysimage_meta(meta_path)
  meta_present = meta !== nothing
  reasons = String[]
  present || push!(reasons, "sysimage file missing")
  meta_present || push!(reasons, "sysimage_meta.toml missing or unreadable")
  if meta_present
    stored_julia = get(meta, "julia_version", "")
    stored_julia == string(VERSION) || push!(reasons, "Julia version changed (image: $(stored_julia), running: $(VERSION))")
    stored_os = get(meta, "os", "")
    stored_os == _sysimage_os_label() || push!(reasons, "operating system changed (image: $(stored_os))")
    stored_arch = get(meta, "arch", "")
    stored_arch == string(Sys.ARCH) || push!(reasons, "architecture changed (image: $(stored_arch))")
    stored_sha = get(meta, "manifest_sha256", "")
    if !isfile(manifest_path)
      push!(reasons, "Manifest.toml missing at $(manifest_path)")
    elseif stored_sha != _sysimage_file_sha256(manifest_path)
      push!(reasons, "Manifest.toml changed since the build (package update); rebuild via tools/build_sysimage.jl")
    end
  end
  return (
    image_path = image_path,
    meta_path = meta_path,
    present = present,
    meta_present = meta_present,
    valid = isempty(reasons),
    reasons = reasons,
    size_bytes = present ? filesize(image_path) : 0,
    built_at = meta_present ? get(meta, "built_at", "") : "",
    julia_version = meta_present ? get(meta, "julia_version", "") : "",
    sparlectra_version = meta_present ? get(meta, "sparlectra_version", "") : "",
    manifest_sha256 = meta_present ? get(meta, "manifest_sha256", "") : "",
  )
end

# Background-build registry: one build at a time, mirroring the job-dict
# style of the PowerFlow web jobs reduced to a single slot. The build always
# runs as a separate Julia process (tools/build_sysimage.jl); PackageCompiler
# must never become a runtime dependency and a build inside the server
# process would block it.
const _SYSIMAGE_BUILD_LOCK = ReentrantLock()
const _SYSIMAGE_BUILD_JOB = Ref{Union{Nothing,Dict{String,Any}}}(nothing)

"""Return the build log path `<user root>/sysimage/sysimage_build.log`."""
webui_sysimage_build_log_path(output_root::AbstractString = default_webui_output_root())::String = joinpath(webui_sysimage_dir(output_root), "sysimage_build.log")

"""
    webui_sysimage_active() -> Bool

True when this process runs on the fast-start sysimage (the launcher passed
`-J <user root>/sysimage/sparlectra.<ext>`). The start script uses this to
skip the JIT warm-up, whose work the image already contains.
"""
function webui_sysimage_active()::Bool
  img = try
    unsafe_string(Base.JLOptions().image_file)
  catch
    return false
  end
  isempty(img) && return false
  return abspath(img) == abspath(webui_sysimage_path())
end

"""
    sysimage_build_state() -> NamedTuple

State of the background fast-start build for the Web UI page: `status`
(`:idle`, `:running`, `:completed`, `:failed`), `running`, and the start and
finish timestamps of the most recent build in this server session.
"""
function sysimage_build_state()
  return lock(_SYSIMAGE_BUILD_LOCK) do
    job = _SYSIMAGE_BUILD_JOB[]
    job === nothing && return (status = :idle, running = false, started_at = "", finished_at = "", elapsed_seconds = 0.0)
    status = Symbol(job["status"])
    elapsed = status === :running ? time() - Float64(get(job, "started_epoch", time())) : 0.0
    (status = status, running = status === :running, started_at = String(get(job, "started_at", "")), finished_at = String(get(job, "finished_at", "")), elapsed_seconds = elapsed)
  end
end

"""
    start_sysimage_build!(output_root; _test_command) -> NamedTuple

Spawn `tools/build_sysimage.jl` as a separate Julia process with stdout and
stderr redirected to the `sysimage_build.log` artifact, register the single
build slot, and emit the operation-log events `sysimage_build_started` and
`sysimage_build_finished` or `sysimage_build_failed`. Returns
`(ok, message)`; a second build while one runs is rejected, never queued.
`_test_command` substitutes the spawned command for the tests (the suite
must not run a real PackageCompiler build).
"""
function start_sysimage_build!(output_root::AbstractString = default_webui_output_root(); operation_log::AbstractString = output_root, _test_command::Union{Nothing,Cmd} = nothing)
  return lock(_SYSIMAGE_BUILD_LOCK) do
    job = _SYSIMAGE_BUILD_JOB[]
    if job !== nothing && job["status"] == "running"
      return (ok = false, message = "A fast-start image build is already running.")
    end
    script = joinpath(SPARLECTRA_ROOT, "tools", "build_sysimage.jl")
    if _test_command === nothing && !isfile(script)
      return (ok = false, message = "Build script not found at $(script); building needs a Sparlectra checkout.")
    end
    log_path = webui_sysimage_build_log_path(output_root)
    mkpath(dirname(log_path))
    julia_bin = joinpath(Sys.BINDIR, Base.julia_exename())
    cmd = _test_command === nothing ? Cmd([julia_bin, "--startup-file=no", script]) : _test_command
    process = run(pipeline(cmd; stdout = log_path, stderr = log_path); wait = false)
    new_job = Dict{String,Any}("status" => "running", "started_at" => Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"), "started_epoch" => time(), "finished_at" => "", "log_path" => log_path, "process" => process)
    _SYSIMAGE_BUILD_JOB[] = new_job
    record_webui_operation!(operation_log, "sysimage_build_started"; log = basename(log_path))
    @async begin
      ok = try
        success(process)
      catch
        false
      end
      lock(_SYSIMAGE_BUILD_LOCK) do
        new_job["status"] = ok ? "completed" : "failed"
        new_job["finished_at"] = Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS")
      end
      record_webui_operation!(operation_log, ok ? "sysimage_build_finished" : "sysimage_build_failed"; log = basename(log_path))
    end
    return (ok = true, message = "Build started.")
  end
end
