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
