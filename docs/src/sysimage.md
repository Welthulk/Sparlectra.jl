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

One call, no arguments. **Plan for it**: the build typically takes 10 to
20 minutes, needs a few GB of RAM, and the image is roughly half a GB on
disk; run it once after installing or updating Sparlectra, not before
every session. It runs the build in a child process (your session keeps
its active project), installs PackageCompiler into the shared
`@sparlectra-sysimage-build` environment on first use, and executes the
full workload: power flow on the embedded cases from
small to large (case3, case14, case118), a state-estimation run with
bad-data diagnostics, an N-1 contingency run, a CGMES import with short
circuit, network losses, and the result/diagnostics printers; the build
then also traces the fast-profile test suite, so everything the tests
touch is compiled into the image as well.
`buildSysimage(dry_run = true)` only reports the target paths.

Equivalent alternatives from a checkout:

- **Script:** `julia tools/build_sysimage.jl`. The build runs in the shared
  environment `@sparlectra-sysimage-build` (PackageCompiler is installed
  there, never into the package `Project.toml`), executes the workload in
  `tools/sysimage_workload.jl`, and writes the image below the Web UI user
  root. The workload covers the paths a session actually hits, one trace
  per import format so no format compiles on its first click: a power flow
  and a state estimation each for MATPOWER (`case14`), SCF (the shipped
  demo case, including its export), DTF (`FOR001.DAT`, reader and net
  builder directly as well) and CGMES (MiniGrid, fetched once into the
  regular case cache if missing), plus a CGMES short circuit, an N-1 run,
  the Web UI start with its form renders, and a clean shutdown. No
  run-history entries are created.

When AnalyticLoadFlow (the APSLF solver) is installed, the build detects it
and bakes it into the image as well. This is not only about the APSLF paths:
`start_webui.jl` loads AnalyticLoadFlow at startup, and any package loaded
*after* the sysimage invalidates precompiled methods inside the image: the
first run then silently recompiles them (measured: 36 s instead of 1 s for
the first `case118` service run). With the package inside the image the
startup load is a no-op for compilation.
Expect a build to take roughly 6 to 20 minutes and an image of a few
hundred MB (reference measurement on a Linux workstation, Julia 1.12:
6.5 minutes build time, 537 MB image). With the image the server is up
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
and the build timestamp. `start_webui.jl` checks three things before every
start and names the one that fails:

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
