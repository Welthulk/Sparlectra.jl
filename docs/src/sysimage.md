# Sysimage

Starting the Web UI pays Julia's just-in-time compilation cost twice: once
for loading the package and once for the first solves. Sparlectra ships two
complementary answers:

1. **PrecompileTools workload (always on).** The package precompiles the
   solver hot path at install time (MATPOWER parse, rectangular
   Newton-Raphson with UMFPACK, losses, DC power flow; issue #288). This
   needs no setup and already covers most of the solver latency.
2. **PackageCompiler sysimage (this page).** A system image bakes
   Sparlectra and its dependencies into one ahead-of-time compiled shared
   library that Julia loads via `-J`. With it, a Web UI start reaches a
   served PowerFlow page in a few seconds; the remaining JIT time
   disappears.

## The start offers to build it

`start_webui.sh` / `start_webui.bat` only launch Julia; the decision itself
lives in `start_webui.jl` (through `tools/sysimage_launcher.jl`), so all
platforms behave identically. On every start it checks whether a usable
image exists:

- **usable:** the process starts itself again with `-J <image>` and prints
  one `Sysimage: ...` line.
- **missing or outdated:** it names the reason and asks
  `Build the sysimage now? [Y/n]`. **No answer means yes** (30 seconds, or
  no terminal at all): an unattended start should end up with the image
  rather than quietly without it. Answering `n` starts without one, with a
  note that every code path then compiles on first use.

Two flags and one environment variable steer it:

```bash
./start_webui.sh --rebuild-sysimage   # build a fresh image even if the current one is fine
./start_webui.sh --no-sysimage        # this start ignores the image and the question
SPARLECTRA_NO_SYSIMAGE=1 ./start_webui.sh   # same, for callers that cannot pass arguments
```

## Building the image directly

The simplest way, from any Julia session with Sparlectra installed:

```julia
using Sparlectra
buildSysimage()
```

One call, no arguments. **Plan for it**: the build takes a few minutes,
needs a few GB of RAM, and the image is roughly half a GB on disk; run it
once after installing or updating Sparlectra, not before every session. It
runs the build in a child process (your session keeps its active project) and
installs PackageCompiler into the shared `@sparlectra-sysimage-build`
environment on first use. `buildSysimage(dry_run = true)` only reports the
target paths.

Equivalent alternatives from a checkout:

- **Script:** `julia tools/build_sysimage.jl`. The build runs in the shared
  environment `@sparlectra-sysimage-build` (PackageCompiler is installed
  there, never into the package `Project.toml`), executes the workload in
  `tools/sysimage_workload.jl`, and writes the image below the Web UI user
  root.

### What the workload traces

The workload runs the important paths ONCE each, and one trace per import
format so no format compiles on its first click: a power flow and a state
estimation each for MATPOWER (`warmup_casePST`), SCF (the shipped demo
case, including its export), DTF (`FOR001.DAT`, reader and net builder
directly as well) and CGMES (MiniGrid, fetched once into the regular case
cache if missing), plus a CGMES short circuit, an N-1 contingency run, the
network losses and the result/diagnostics printers, the Web UI start with
its page renders through the real socket handler, one `run_sparlectra`
call, and a clean shutdown. No run-history entries are created.

It deliberately does NOT run the test suite. Tracing the whole fast profile
plus the Web UI test group is what the build used to do, and it dominated
the build time by a wide margin without reaching a Web UI path the steps
above miss. For a maintainer comparison the old behavior is still one
variable away:

```bash
SPARLECTRA_SYSIMAGE_TRACE_TESTS=1 julia --project=. tools/build_sysimage.jl
```

A step that fails is a gap in the trace, never a failed build: the affected
path simply compiles on its first use. Every step reports its outcome in
the build log, successful ones included, because a silent success cannot be
told apart from a step that never ran.

### What the build prints

The console gets one line that rewrites itself:

```text
  [3/4] compiling the system image  04:21
```

Everything else - Pkg resolution, the workload trace, PackageCompiler, the
linker - goes to `sysimage_build.log` next to the image, and is written as
it happens rather than at the end, so it can be watched and survives a
killed build. When the console is not a terminal (output redirected to a
file, a CI log), the progress line is replaced by one line per phase
change, about fifteen for a whole build. `--verbose` streams everything to
the console instead.

A build that fails prints the tail of the log directly, and leaves the
PREVIOUS image untouched: the new one is compiled to a staging file and
moved into place only after it is complete.

### Refreshing the image from the Web UI

The Web UI's **Sysimage** page (`/webui/sysimage`, linked from the Info
panel) shows the same validity verdict and starts the same build in the
background. Because the finished image is moved into place with `mv`, the
directory entry is replaced while the old file stays mapped, so a Web UI
running on the image being replaced keeps working. It also keeps running
the OLD code: the page says so afterwards and asks for a restart. See
[Web UI](webui.md).

When AnalyticLoadFlow (the APSLF solver) is installed, the build detects it
and bakes it into the image as well. This is not only about the APSLF paths:
`start_webui.jl` loads AnalyticLoadFlow at startup, and any package loaded
*after* the sysimage invalidates precompiled methods inside the image: the
first run then silently recompiles them (measured: 36 s instead of 1 s for
the first `case118` service run). With the package inside the image the
startup load is a no-op for compilation.
A build takes a few minutes and produces an image of a few hundred MB. The
time splits roughly evenly between resolving and precompiling the build
environment, the workload trace, and the image compilation itself; a cold
build environment adds the one-time PackageCompiler installation.
With the image the server is up
after about 2 s; the first PowerFlow form view still pays a one-time scan
of the case directory (prewarmed in the background, afterwards a few
milliseconds per view). Without the image the same start took about 45 s.
A freshly built image is used by the relaunch that follows the build, so it
takes effect immediately.

## Standalone executable: a tool, not a product feature

Sparlectra runs on an installed Julia. A relocatable executable with an
embedded Julia runtime (for a machine without Julia) is therefore not
part of the package and not offered in the Web UI; the build script lives
in the checkout and is run directly:

```bash
julia --project=. tools/build_app.jl --flavor=full      # Web UI plus command line
julia --project=. tools/build_app.jl --flavor=runtime   # library runtime, no GUI
julia --project=. tools/build_app.jl --script=my_workflow.jl
```

Plan for 20 to 40 minutes and an app folder of 0.5 to 1 GB. There is no
cross-compilation: build on the operating system you target, and rebuild
after updating Sparlectra. The resulting `bin/sparlectra` exposes the
commands `webui`, `run`, `se`, `n1` and `version`, takes the normal
library YAML through `--config`/`--set`, and `bin/sparlectra help` prints
its own overview.

### The launchers read no startup file

`start_webui.sh` / `start_webui.bat` and the relaunch through the image start
Julia with `--startup-file=no`. The Web UI is a server process, not a REPL, and
a personal `startup.jl` typically loads `Revise`: loading Revise on top of the
Sparlectra image INVALIDATES a large part of it, and the invalidated methods
are then inferred again on first use. Building an image and throwing much of it
away at startup is worse than not building one. It was paid twice, because the
launcher and the relaunched child both read the file.

Set `SPARLECTRA_STARTUP_FILE=yes` to get the old behavior back for a session
where Revise in the server process is actually wanted.

### The configuration is brought up to date

A build first rewrites the Web UI configuration to the current key layout
(`refresh_sparlectra_config_file`, keeping a timestamped backup) and says so in
one line. A configuration that cannot be migrated automatically, a file with
duplicate YAML keys for example, stops the build with the reason named: an
image built against a configuration the user still has to edit by hand is an
image they cannot use. Start with `SPARLECTRA_NO_SYSIMAGE=1` until the file is
fixed, then build again.

## Where the image lives

`<user root>/sysimage/sparlectra.<ext>` next to its metadata
`sysimage_meta.toml`, with the platform extension `so` (Linux), `dylib`
(macOS), `dll` (Windows). The user root is the same application root that
already holds runs, logs, configuration, and the MATPOWER cache:

| Platform | Location |
|---|---|
| Linux | `$XDG_STATE_HOME/sparlectra/webui/sysimage/` (default `~/.local/state/...`) |
| macOS | `~/Library/Application Support/Sparlectra/WebUI/sysimage/` |
| Windows | `%LOCALAPPDATA%\Sparlectra\WebUI\sysimage\` |

Note that the user root follows the environment: a Flatpak-packaged IDE
(for example VSCodium) sets its own `XDG_STATE_HOME`, so a build started
from its integrated terminal lands under
`~/.var/app/<app-id>/.local/state/...` while a plain terminal uses
`~/.local/state/...`. Two separate roots mean two separate configurations,
case caches, run histories, and images, and settings silently drift apart
between them. The recommended permanent fix is one shared root via a
symlink (merge any content you want to keep first):

```bash
mv ~/.var/app/<app-id>/.local/state/sparlectra ~/.var/app/<app-id>/.local/state/sparlectra.backup
ln -s ~/.local/state/sparlectra ~/.var/app/<app-id>/.local/state/sparlectra
```

The run registry resolves symlinks during validation, so runs written from
either environment stay visible in both. Alternatively copy
`sparlectra.<ext>` plus `sysimage_meta.toml` between the roots, or rebuild
from the environment you start from.

## Staleness rules

`sysimage_meta.toml` is the validity contract: Sparlectra version, Julia
version, OS and architecture, the SHA-256 of the checkout `Manifest.toml`,
and the build timestamp. `start_webui.jl` checks these before every start,
and the Web UI's Sysimage page shows the same verdict, naming the one that
fails:

| Check | Message |
|---|---|
| image and metadata present | `no sysimage found` |
| Julia version matches | `the sysimage was built for Julia X, this is Y` |
| `Manifest.toml` hash matches | `the sysimage does not match the current Manifest.toml` |
| no `src/*.jl` newer than the image | `the sysimage is older than <checkout>/src` |

The last row is the one that cost an evening. The metadata pins the
DEPENDENCIES, so editing the package source leaves it untouched: an image
built before the edit still looked fresh, and the Web UI silently served
the old code while the change was "not there". A source file newer than the
image now disables it, which is slower but never wrong. Released
installations never hit this (package updates always change the manifest);
in a development checkout it means a rebuild after every edit, or
`--no-sysimage` while working.

Package and Julia updates disable the image the same way, until the next
rebuild. There is no silent automatic rebuild: the start says what is wrong
and asks.

## Escape hatch

`SPARLECTRA_NO_SYSIMAGE=1` makes the launchers skip the image
unconditionally. Use it to rule the sysimage out when debugging suspected
invalidation or precompilation issues:

```bash
SPARLECTRA_NO_SYSIMAGE=1 ./start_webui.sh
```
