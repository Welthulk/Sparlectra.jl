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

# file: src/build/sysimage_builder.jl
# purpose: one-call APIs around the sysimage build
#          (tools/build_sysimage.jl) and the relocatable app build
#          (tools/build_app.jl). Each build runs in a CHILD process so
#          the caller's package environment is never activated away and
#          PackageCompiler stays out of the package dependencies (it lives
#          in the shared @sparlectra-sysimage-build environment).

"""
    buildSysimage(; dry_run = false, quiet = false) -> NamedTuple

Build the Sparlectra sysimage with one call.

The build precompiles the interactive paths through a workload
(power flow on the embedded cases small to large: case3, case14,
case118; a state-estimation run with bad-data diagnostics; an N-1
contingency run; a CGMES import with short circuit; network losses and the
result/diagnostics printers), so the Web UI afterwards starts in seconds
instead of paying JIT compilation on every first click.

**This takes a while**: typically 10 to 20 minutes of CPU time, a few GB of
RAM, and roughly half a GB of disk for the image. Run it once after
installing or updating Sparlectra, not before every session. The launch
scripts (`start_webui.sh` / `start_webui.bat`) pick the image up
automatically on the next start; bypass it with `SPARLECTRA_NO_SYSIMAGE=1`.

`PackageCompiler` is installed into the shared `@sparlectra-sysimage-build`
environment on first use (one-time network access); the package's own
dependencies stay untouched. The build runs in a child Julia process, so
the calling session keeps its active project.

Keyword arguments:
- `dry_run = true` only reports the target paths without building.
- `quiet = true` suppresses the build output.

Returns `(sysimage_path, meta_path, built, size_mb)`.
"""
function buildSysimage(; dry_run::Bool = false, quiet::Bool = false)
  pkgroot = pkgdir(@__MODULE__)
  pkgroot === nothing && error("buildSysimage: cannot locate the Sparlectra package directory")
  script = joinpath(pkgroot, "tools", "build_sysimage.jl")
  isfile(script) || error("buildSysimage: build script not found at $(script)")
  img = webui_sysimage_path()
  meta = webui_sysimage_meta_path()
  if !dry_run && !quiet
    @info "Building the Sparlectra sysimage. This typically takes 10 to 20 minutes; afterwards the Web UI starts in seconds instead of paying JIT on every first click." target = img
  end
  # run the script against THIS package's project: without it the child
  # would load whatever Sparlectra version the default environment holds
  # (seen: an older registered version without the sysimage helpers).
  #
  # The plain executable, NOT Base.julia_cmd(): that command repeats the -J of
  # the calling session, and when the caller already runs on the Sparlectra
  # image (the Web UI does, and it offers a refresh button) PackageCompiler
  # would build the new image incrementally on top of the old one instead of
  # on the stock Julia image.
  exe = joinpath(Sys.BINDIR, Base.julia_exename())
  cmd = dry_run ? `$(exe) --startup-file=no --project=$(pkgroot) $(script) --dry-run` : `$(exe) --startup-file=no --project=$(pkgroot) $(script)`
  io_out = quiet ? devnull : stdout
  io_err = quiet ? devnull : stderr
  run(pipeline(Cmd(cmd; dir = pkgroot); stdout = io_out, stderr = io_err))
  built = !dry_run && isfile(img)
  return (sysimage_path = img, meta_path = meta, built = built, size_mb = built ? round(filesize(img) / 1024^2; digits = 1) : 0.0)
end
