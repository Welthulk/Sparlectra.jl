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
# purpose: build the sysimage for the Web UI with PackageCompiler in the
#          shared build environment @sparlectra-sysimage-build (the package
#          Project.toml stays free of PackageCompiler). The image lands below
#          the Web UI user root next to its sysimage_meta.toml validity
#          contract; the launchers pick it up automatically. Run from the
#          checkout: `julia tools/build_sysimage.jl`.
#
#          The console shows ONE progress line and nothing else; the full
#          output of Pkg, the workload and PackageCompiler goes to
#          sysimage_build.log next to the image. The same run writes
#          sysimage_build_progress.toml, which the Web UI reads to show the
#          progress of a build it started itself.
#
#          Flags: `--dry-run` prints the plan without building (the test
#          suite uses it as a parse smoke), `--verbose` streams everything to
#          the console instead of the progress line.
#          Library users call the same build through the exported
#          Sparlectra.buildSysimage() one-liner, which runs this script in
#          a child process.

const _REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const _BUILD_ENV = "sparlectra-sysimage-build"
const _DRY_RUN = ("--dry-run" in ARGS) || get(ENV, "SPARLECTRA_SYSIMAGE_DRY_RUN", "0") == "1"
const _VERBOSE = ("--verbose" in ARGS) || get(ENV, "SPARLECTRA_SYSIMAGE_VERBOSE", "0") == "1"

# Pkg is deliberately NOT loaded at the top level: the dry-run path (used by
# the test-suite smoke check) must work in environments whose load path has
# no Pkg, such as the Pkg.test sandbox. The real build path loads Pkg inside
# main() and resolves the module dynamically (Julia 1.12 world age).

# The launcher module owns the sysimage location on plain Base, which is what
# this script needs BEFORE `using Sparlectra`: the log and the progress file
# have to exist from the first second, and loading the package to learn where
# they go would put the slowest step ahead of the first progress line.
# Base.include with an explicit module, not a bare `include`: the test suite
# runs this script's dry-run path inside a fresh `Module()`, which has no
# `include` of its own.
Base.include(@__MODULE__, joinpath(@__DIR__, "sysimage_launcher.jl"))

# --- console versus log ------------------------------------------------------
# `_CONSOLE` is captured BEFORE any redirection, because `stdout` inside the
# teeing block below IS the pipe. Everything the build prints (Pkg resolution,
# the workload trace, PackageCompiler, the C compiler) is written to the log;
# the console only ever carries the progress line and the final summary.
const _CONSOLE = stdout
const _CONSOLE_IS_TTY = isa(_CONSOLE, Base.TTY)

const _SYSIMAGE_DIR = dirname(SysimageLauncher.sysimage_path())
const _BUILD_LOG_PATH = joinpath(_SYSIMAGE_DIR, "sysimage_build.log")
const _PROGRESS_PATH = joinpath(_SYSIMAGE_DIR, "sysimage_build_progress.toml")

# The workload announces every step it enters with this prefix on stdout. It
# is the ONLY channel between the two processes: PackageCompiler runs the
# workload in a child whose output this script already tees, so a marker line
# costs nothing and needs no second file.
const PROGRESS_MARKER = "@@SPARLECTRA_STEP@@"

# Last line the workload prints; see tools/sysimage_workload.jl.
const WORKLOAD_DONE_MARKER = "workload done"

const _BUILD_STEPS = 4

mutable struct _Progress
  state::String        # running | done | failed
  step::Int
  phase::String
  detail::String
  message::String
  started::Float64
  step_started::Float64
end

const _PROGRESS = _Progress("running", 0, "starting", "", "", time(), time(), )

_fmt_clock(seconds::Real) = string(lpad(Int(fld(seconds, 60)), 2, '0'), ":", lpad(Int(floor(seconds)) % 60, 2, '0'))

_toml_escape(s::AbstractString) = replace(String(s), "\\" => "\\\\", "\"" => "\\\"", "\n" => " ")

"""
Write the progress file the Web UI polls. Failures are ignored on purpose and
named here: the state directory can be read-only or full, and a build that
cannot report its progress is still a build worth finishing. The console line
is unaffected, so nothing becomes invisible.
"""
# The progress file is rewritten every few seconds; one warning is enough.
const _PROGRESS_WRITE_WARNED = Ref(false)

function _write_progress()
  try
    mkpath(_SYSIMAGE_DIR)
    open(_PROGRESS_PATH, "w") do io
      println(io, "state = \"", _toml_escape(_PROGRESS.state), "\"")
      println(io, "step = ", _PROGRESS.step)
      println(io, "steps = ", _BUILD_STEPS)
      println(io, "phase = \"", _toml_escape(_PROGRESS.phase), "\"")
      println(io, "detail = \"", _toml_escape(_PROGRESS.detail), "\"")
      println(io, "message = \"", _toml_escape(_PROGRESS.message), "\"")
      println(io, "elapsed_seconds = ", round(time() - _PROGRESS.started; digits = 1))
      println(io, "updated_at = ", round(time(); digits = 1))
      println(io, "pid = ", getpid())
      println(io, "log = \"", _toml_escape(_BUILD_LOG_PATH), "\"")
    end
  catch err
    # Expected failure: unwritable or full state directory. Said ONCE - the
    # file is rewritten every few seconds, so a message per attempt would
    # bury the build; saying nothing at all would leave the Web UI's progress
    # display mysteriously empty with no reason anywhere.
    if !_PROGRESS_WRITE_WARNED[]
      _PROGRESS_WRITE_WARNED[] = true
      _console("Cannot write $(_PROGRESS_PATH) ($(sprint(showerror, err))); the build continues, but the Web UI cannot show its progress.")
    end
  end
  return nothing
end

# Without a terminal the in-place progress line is invisible (`\r` needs one),
# and a build that prints nothing for ten minutes is indistinguishable from a
# hung one. So a redirected or piped console gets ONE line per phase change
# instead - about fifteen lines for a whole build, which is still quiet.
function _announce(text::AbstractString)
  (_VERBOSE || !_CONSOLE_IS_TTY) || return nothing
  println(_CONSOLE, "  [", max(_PROGRESS.step, 1), "/", _BUILD_STEPS, "] ", text)
  flush(_CONSOLE)
  return nothing
end

function _set_step(step::Int, phase::AbstractString)
  _PROGRESS.step = step
  _PROGRESS.phase = String(phase)
  _PROGRESS.detail = ""
  _PROGRESS.step_started = time()
  _write_progress()
  _announce(phase)
  return nothing
end

function _set_detail(detail::AbstractString)
  detail == _PROGRESS.detail && return nothing
  _PROGRESS.detail = String(detail)
  _write_progress()
  _announce(string(_PROGRESS.phase, " - ", detail))
  return nothing
end

# Width of the last console line, so the next one can overwrite it completely;
# a shorter phase name must not leave the tail of a longer one on screen.
const _LINE_WIDTH = Ref(0)

function _paint_progress()
  _CONSOLE_IS_TTY || return nothing
  elapsed = _fmt_clock(time() - _PROGRESS.started)
  detail = isempty(_PROGRESS.detail) ? "" : string(" - ", _PROGRESS.detail)
  line = string("  [", max(_PROGRESS.step, 1), "/", _BUILD_STEPS, "] ", _PROGRESS.phase, detail, "  ", elapsed)
  length(line) > 100 && (line = string(first(line, 97), "..."))
  print(_CONSOLE, "\r", rpad(line, max(_LINE_WIDTH[], length(line))))
  flush(_CONSOLE)
  _LINE_WIDTH[] = length(line)
  return nothing
end

"Clear the in-place progress line so the next real output starts on a clean row."
function _clear_progress_line()
  _CONSOLE_IS_TTY || return nothing
  _LINE_WIDTH[] == 0 && return nothing
  print(_CONSOLE, "\r", " "^_LINE_WIDTH[], "\r")
  flush(_CONSOLE)
  _LINE_WIDTH[] = 0
  return nothing
end

"Print a line on the console, above the progress line."
function _console(msg::AbstractString)
  _clear_progress_line()
  println(_CONSOLE, msg)
  flush(_CONSOLE)
  return nothing
end

"""
Classify a captured output line. `nothing` means "log only"; a string is the
detail to show next to the phase. Only lines that state PROGRESS qualify:
the workload's own step markers and the few PackageCompiler lines that mark a
transition. Everything else stays in the log, which is the point of the
exercise (maintainer, 2026-09-07: "wenig bis keine Ausgaben in der cmd").
"""
function _progress_detail(line::AbstractString)
  startswith(line, PROGRESS_MARKER) && return strip(line[(length(PROGRESS_MARKER)+1):end])
  # Two independent signs that the trace is over and the compilation has
  # started. PackageCompiler's own line is the precise one; the workload's
  # last line is the fallback, so a reworded PackageCompiler message leaves
  # the progress stuck at "tracing" instead of never moving at all.
  occursin("PackageCompiler: compiling", line) && return "compiling"
  startswith(line, WORKLOAD_DONE_MARKER) && return "compiling"
  occursin("Precompiling", line) && return "precompiling dependencies"
  occursin(r"^\s*(Installed|Resolving|Updating|Added|Downloaded)\b", line) && return "resolving packages"
  return nothing
end

"Tail of the build log, for the console when the build fails."
function _log_tail(n::Int = 25)
  isfile(_BUILD_LOG_PATH) || return String[]
  lines = try
    readlines(_BUILD_LOG_PATH)
  catch err
    # expected failure: the log was removed or is unreadable. This runs on
    # the failure path, where an empty tail would look like "no reason given"
    return ["(build log unreadable: $(sprint(showerror, err)))"]
  end
  filter!(l -> !isempty(strip(l)), lines)
  return lines[max(1, end - n + 1):end]
end

function _prepare_build_env(pkgm::Module)::Bool
  println("activating shared build environment @", _BUILD_ENV)
  Base.invokelatest(pkgm.activate, _BUILD_ENV; shared = true)
  if Base.find_package("PackageCompiler") === nothing
    println("installing PackageCompiler into the build environment (one-time)")
    Base.invokelatest(pkgm.add, "PackageCompiler")
  end
  # match how the script was invoked: from a checkout the repo itself is
  # developed into the build environment, otherwise the released package
  project_file = joinpath(_REPO_ROOT, "Project.toml")
  in_checkout = isfile(project_file) && occursin("name = \"Sparlectra\"", read(project_file, String))
  if in_checkout
    println("developing Sparlectra from ", _REPO_ROOT)
    Base.invokelatest(pkgm.develop; path = _REPO_ROOT)
  else
    println("adding the released Sparlectra package")
    Base.invokelatest(pkgm.add, "Sparlectra")
  end
  # AnalyticLoadFlow (the APSLF solver) is a REQUIRED dependency, so resolving
  # Sparlectra brings it into the build environment and the image contains it.
  # That matters beyond convenience: loading a package AFTER the sysimage
  # invalidates precompiled methods inside the image, and the first run then
  # silently recompiles them (measured 36 s instead of 1 s on case118).
  Base.invokelatest(pkgm.instantiate)
  return true
end

"""
Raised when the user configuration carries legacy keys that cannot be migrated
automatically. The build stops on it: an image built against a configuration
the user still has to edit by hand is an image they cannot use.
"""
struct ConfigNotMigratable <: Exception
  path::String
  reasons::Vector{String}
end

"""
Bring the user's Web UI configuration up to the current key layout, in place.

At most ONE console line, and only when something actually changed (maintainer,
2026-09-07). A version-0 file otherwise produces an alias notice at every
single start, and telling the user to run a function by hand is work the build
can simply do: it already writes to that directory and it already takes
minutes. A timestamped backup is kept.

Returns normally when the configuration is usable, and throws
`ConfigNotMigratable` when it is not: duplicate YAML keys need human judgement,
and guessing which of two values was meant is the one thing a migration must
never do.
"""
function _refresh_user_config(spar::Module)
  path = Base.invokelatest(getglobal(spar, :default_webui_config_path))
  isfile(path) || return nothing          # nothing provisioned yet, nothing to migrate
  result = try
    Base.invokelatest(getglobal(spar, :refresh_sparlectra_config_file), path; write = true, backup = true)
  catch err
    # expected failure: an unreadable or syntactically broken YAML. The user
    # has to fix that themselves, so it is reported, not swallowed.
    throw(ConfigNotMigratable(path, [sprint(showerror, err)]))
  end
  result.changed || return nothing
  result.written || throw(ConfigNotMigratable(path, isempty(result.warnings) ? ["the refreshed file could not be written"] : result.warnings))
  backup = result.backup_path === nothing ? "" : " (previous version kept as $(basename(result.backup_path)))"
  _console("Configuration updated to the current key names: $(path)$(backup)")
  return nothing
end

"""
The build itself, with stdout/stderr already teed into the log by `main`.
Returns the finished image path.

The image is compiled to a STAGING file next to the target and moved into
place at the very end. Two reasons, and the second one is why the Web UI can
offer a refresh at all: a build that dies halfway leaves the previous working
image untouched, and `mv` replaces the directory entry without writing into
the old inode, so a Web UI currently RUNNING on that image keeps its mapped
pages instead of taking a SIGBUS mid-request.
"""
function _run_build()
  @eval using Pkg
  pkgm = Base.invokelatest(getfield, @__MODULE__, :Pkg)
  _set_step(1, "preparing the build environment")
  _prepare_build_env(pkgm)

  @eval using PackageCompiler
  @eval using Sparlectra
  # same world-age pattern as above: both modules were loaded inside this
  # function invocation
  spar = Base.invokelatest(getfield, @__MODULE__, :Sparlectra)
  pc = Base.invokelatest(getfield, @__MODULE__, :PackageCompiler)

  # before anything expensive: the image is built FOR this configuration
  _refresh_user_config(spar)

  img = Base.invokelatest(getglobal(spar, :webui_sysimage_path))
  meta_path = Base.invokelatest(getglobal(spar, :webui_sysimage_meta_path))
  manifest = Base.invokelatest(getglobal(spar, :webui_sysimage_manifest_path))
  mkpath(dirname(img))
  # keep the platform extension: PackageCompiler and the linker both key off it
  staging = string(img, ".building.", Base.invokelatest(getglobal(spar, :webui_sysimage_ext)))
  isfile(staging) && rm(staging; force = true)

  # The MATPOWER, SCF and DTF traces run on TRACKED cases (data/mpower,
  # data/scf, data/DTF), so the build needs no network for them. It used to
  # fetch case14.m here, which is not shipped: on a fresh checkout that made
  # the build depend on a download, and on this machine it only worked
  # because an older run had left the file in the cache.
  # The MiniGrid CGMES delivery is the one exception, because a CGMES
  # delivery cannot be assembled from the raw model directories; a failed
  # fetch only costs that trace (the workload logs the gap), never the build.
  case_cache = Base.invokelatest(getglobal(spar, :default_webui_case_cache_dir))
  if !isfile(joinpath(case_cache, "cgmes_minigrid.zip"))
    println("fetching the MiniGrid CGMES delivery for the short-circuit workload")
    try
      cgmes = getglobal(spar, :CGMESImporter)
      Base.invokelatest(getglobal(cgmes, :fetchCGMESTestSet), "minigrid"; outdir = case_cache)
    catch err
      println("MiniGrid fetch failed (", sprint(showerror, err), "); the short-circuit service path will not be traced")
    end
  end

  # both are baked in: AnalyticLoadFlow is a required dependency, and a
  # package loaded after the image would invalidate its precompiled methods
  packages = ["Sparlectra", "AnalyticLoadFlow"]
  workload = joinpath(@__DIR__, "sysimage_workload.jl")
  _set_step(2, "tracing the interactive paths")
  println("target: ", img)
  println("staging: ", staging)
  println("packages: ", join(packages, ", "))
  Base.invokelatest(getglobal(pc, :create_sysimage), packages; sysimage_path = staging, precompile_execution_file = workload)

  _set_step(4, "installing the image")
  # mv over the live image: atomic, and the old inode survives for whoever
  # still has it mapped (see the docstring)
  mv(staging, img; force = true)
  Base.invokelatest(getglobal(spar, :write_sysimage_meta), meta_path; manifest_path = manifest)
  return img
end

function main()
  started = time()
  _PROGRESS.started = started
  if _DRY_RUN
    # parse-and-plan mode for the test suite: resolve everything that does
    # not touch the package environment or PackageCompiler
    println("[build_sysimage] dry run: no build environment changes, no PackageCompiler call")
    if Base.find_package("Sparlectra") !== nothing
      @eval using Sparlectra
      # the module binding is newer than this function invocation (Julia
      # 1.12 world age); resolve it dynamically via invokelatest
      spar = Base.invokelatest(getfield, @__MODULE__, :Sparlectra)
      img = Base.invokelatest(getglobal(spar, :webui_sysimage_path))
      meta = Base.invokelatest(getglobal(spar, :webui_sysimage_meta_path))
      println("[build_sysimage] would build: ", img)
      println("[build_sysimage] would write: ", meta)
    else
      println("[build_sysimage] Sparlectra not loadable in this environment; dry run ends after parsing")
    end
    println("[build_sysimage] dry run finished")
    return nothing
  end

  mkpath(_SYSIMAGE_DIR)

  # Refuse a second concurrent build. There are now two ways to start one -
  # the launcher question and the Web UI's refresh button - and they would
  # otherwise compile into the same staging file and overwrite each other's
  # log and progress. A builder that stopped reporting for more than a
  # minute counts as gone, so a killed build cannot block the next one
  # (same rule as Sparlectra.sysimage_build_active).
  if isfile(_PROGRESS_PATH)
    running = try
      other = SysimageLauncher.TOML.parsefile(_PROGRESS_PATH)
      updated = get(other, "updated_at", 0.0)
      get(other, "state", "") == "running" && updated isa Number && time() - Float64(updated) < 60.0
    catch err
      # expected failure: the other builder is writing the file right now, so
      # a read can land mid-line. Treat an unreadable file as "no build".
      _console("Ignoring an unreadable build-progress file ($(sprint(showerror, err))).")
      false
    end
    if running
      _console("A sysimage build is already running (see $(_BUILD_LOG_PATH)); this one does nothing.")
      exit(0)
    end
  end

  log = try
    open(_BUILD_LOG_PATH, "w")
  catch err
    # expected failure: unwritable state directory. Without a log there is
    # nothing to tee into, so the build runs loud rather than blind.
    _console("Cannot write $(_BUILD_LOG_PATH) ($(sprint(showerror, err))); building with full output.")
    devnull
  end

  _console("Building the Sparlectra sysimage. Detail: " * _BUILD_LOG_PATH)
  _console("This takes a few minutes; the line below updates while it works.")
  _write_progress()

  # Tee: everything the build writes goes to the log, and the few lines that
  # state progress update the console line. redirect_stdio replaces the
  # process file descriptors, so the PackageCompiler child process and the
  # workload it runs are captured as well (verified: a `run(...)` inside the
  # block lands in the log).
  pipe = Pipe()
  Base.link_pipe!(pipe; reader_supports_async = true, writer_supports_async = true)
  reader = @async begin
    for line in eachline(pipe)
      println(log, line)
      # flushed per line, not at the end: the log is what the Web UI offers
      # to look at WHILE the build runs, and it is the only evidence left
      # when a build is killed. A buffered log is empty exactly when it is
      # needed. The write cost is nothing next to compiling.
      flush(log)
      detail = _progress_detail(line)
      if detail !== nothing
        # the trace is over once PackageCompiler starts compiling
        detail == "compiling" ? _set_step(3, "compiling the system image") : _set_detail(detail)
      end
      _VERBOSE && println(_CONSOLE, "    | ", line)
    end
    flush(log)
  end
  # One task owns the console line AND the heartbeat in the progress file.
  # The heartbeat is what lets the Web UI tell a build that is busy compiling
  # (no output for minutes) from one that died: a progress file older than a
  # minute means gone, not slow.
  painter = @async begin
    ticks = 0
    while _PROGRESS.state == "running"
      # in verbose mode the raw lines scroll past, so the in-place line would
      # only fight with them; the heartbeat below still runs
      _VERBOSE || _paint_progress()
      ticks += 1
      ticks % 5 == 0 && _write_progress()
      sleep(1.0)
    end
  end

  img = nothing
  failure = nothing
  try
    redirect_stdio(stdout = pipe.in, stderr = pipe.in) do
      img = Base.invokelatest(_run_build)
    end
  catch err
    failure = err
    println(log, "\n=== sysimage build failed ===")
    showerror(log, err, catch_backtrace())
    println(log)
  finally
    try
      close(pipe.in)
    catch err
      # expected failure: the write end is already closed because the build
      # process died with it. Nothing to do, but not worth hiding either.
      _VERBOSE && _console("Closing the capture pipe raised: $(sprint(showerror, err))")
    end
    try
      wait(reader)
    catch err
      # the tee task itself failed, so the log below is incomplete: say so
      # rather than letting a truncated log look like the whole story
      _console("The build-log capture failed ($(sprint(showerror, err))); $(_BUILD_LOG_PATH) is incomplete.")
    end
    _PROGRESS.state = failure === nothing ? "done" : "failed"
    try
      wait(painter)
    catch
    end
    _clear_progress_line()
    try
      flush(log)
      log === devnull || close(log)
    catch
    end
  end

  elapsed = round((time() - started) / 60; digits = 1)
  if failure isa ConfigNotMigratable
    # not a build defect: the user has to edit their file, and until then the
    # Web UI has to run without an image
    _PROGRESS.state = "failed"
    _PROGRESS.message = "the configuration needs manual changes"
    _write_progress()
    _console("No sysimage was built: the configuration needs changes that cannot be made automatically.")
    _console("  File: $(failure.path)")
    for reason in failure.reasons
      _console("  " * reason)
    end
    _console("Start without an image until it is fixed:  SPARLECTRA_NO_SYSIMAGE=1  (or --no-sysimage),")
    _console("edit the file, then build again.")
    exit(2)
  end
  if failure !== nothing
    _PROGRESS.message = first(sprint(showerror, failure), 300)
    _write_progress()
    _console("The sysimage build FAILED after $(elapsed) min. Last lines of $(_BUILD_LOG_PATH):")
    for line in _log_tail()
      _console("  " * line)
    end
    _console("The previous image (if any) was left untouched.")
    exit(1)
  end

  size_mb = round(filesize(img) / 1024^2; digits = 1)
  _PROGRESS.message = "$(size_mb) MB in $(elapsed) min"
  _write_progress()
  _console("Sysimage ready in $(elapsed) min: $(img) ($(size_mb) MB)")
  _console("It is used automatically on the next start; bypass with SPARLECTRA_NO_SYSIMAGE=1.")
  return nothing
end

Base.invokelatest(main)
