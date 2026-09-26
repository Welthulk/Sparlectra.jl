# Sysimage

A fresh Julia process compiles Sparlectra and every solver path on first
use. A PackageCompiler sysimage (Sparlectra and its dependencies as one
ahead-of-time compiled library, loaded with `-J`) removes that wait:
`using Sparlectra` returns at once and a Web UI start serves its page in
seconds. The Web UI start offers to build the image; the
**Sysimage** page refreshes it from the browser.

## Build

1. Start the Web UI: `start_webui.sh` / `start_webui.bat` run
   `start_webui.jl`, deciding through `tools/sysimage_launcher.jl`.
   A usable image restarts the process with `-J <image>` and prints one
   `Sysimage: ...` line.
2. A missing or outdated image is named and the start asks
   `Build the sysimage now? [y/N]`: `y` or `--rebuild-sysimage` builds, no
   answer (Enter, 30 seconds, no terminal) means no, and every path
   compiles on first use.
3. The build first migrates the Web UI configuration to the current key
   layout (`refresh_sparlectra_config_file`, timestamped backup);
   duplicate YAML keys, for example, stop it with the reason named: start
   with `SPARLECTRA_NO_SYSIMAGE=1` until fixed.
4. A child process installs PackageCompiler into the shared
   `@sparlectra-sysimage-build` environment, not the package
   `Project.toml`, includes AnalyticLoadFlow (the APSLF solver) when
   installed, and needs minutes, a few GB of RAM and half a GB on disk,
   once per install or update.
5. The relaunch uses the fresh image; the first PowerFlow form view still
   scans the case directory once, prewarmed in the background.

| Entry | Build call |
|---|---|
| Command line | `buildSysimage()` with `SparlectraApp` loaded (`dry_run = true` reports the target paths only, `target_dir` chooses another output directory), or `julia --project=app tools/build_sysimage.jl` in a checkout (workload `tools/sysimage_workload.jl`) |
| Web UI | the **Sysimage** page (`/webui/sysimage`, linked from the Info panel) shows the verdict of the start and runs the same build in the background |
| One-line installers | build the image only with `SPARLECTRA_BUILD_SYSIMAGE=1` |

```bash
./start_webui.sh --rebuild-sysimage   # build a fresh image even if the current one is fine
./start_webui.sh --no-sysimage        # this start ignores the image and the question
```

**What the build prints.** One progress line
(`[3/4] compiling the system image  04:21`; a line per phase change
without a terminal; everything with `--verbose`); the full log is
`sysimage_build.log` next to the image, written as it happens. A failed
build prints its tail and leaves the previous image untouched.

**Refreshing from the Web UI.** The new image is moved into place, so a
Web UI running on the old one keeps working, on the old code; the page
says so and asks for a restart. A REPL session (or `julia --project=app`)
has no image and cannot switch; the page shows the line to use: the start
script, or
`julia -J <image> --startup-file=no --project=<app>` ([Web UI](webui.md)).

**Your own scripts.** Any Julia call takes the image with `--sysimage`:

```bash
# buildSysimage() prints the image path. Rebuild after updating Sparlectra;
# an image cannot be moved between operating systems.
julia --sysimage /path/to/sparlectra.so --project=/path/to/project my_script.jl
julia --sysimage /path/to/sparlectra.so --project=/path/to/project -e 'using Sparlectra; runpf!(importSCF("case.scf.json"))'
# the first runpf!, runse! or runShortCircuit! runs at full speed: batch script,
# REPL, notebook kernel or scheduled job
```

## Workload

Each important path is traced once per import format, so no format
compiles on its first click. A failed step is a gap in the trace, not a
failed build; every step logs its outcome; no run-history entries are
created.

| Step | What is traced |
|---|---|
| MATPOWER | a power flow and a state estimation on `warmup_casePST`; one `run_sparlectra` call |
| SCF | power flow and state estimation on the shipped demo case, including its export |
| DTF | `FOR001.DAT`: power flow and state estimation, reader and net builder directly as well |
| CGMES | MiniGrid (fetched once into the regular case cache if missing): power flow, state estimation and a short circuit |
| N-1 | one contingency run |
| results | the network losses and the result and diagnostics printers |
| Web UI | the start with its page renders through the real socket handler, then a clean shutdown |
| APSLF | the solver and the APSLF-seeded rectangular start (AnalyticLoadFlow brings its own PrecompileTools workload) |
| test suite | not traced by default: it dominates the build time without reaching a Web UI path the steps above miss; `SPARLECTRA_SYSIMAGE_TRACE_TESTS=1` adds it for a comparison |

## Staleness

`sysimage_meta.toml` next to the image is the validity contract (Sparlectra
version, Julia version, OS and architecture, SHA-256 of the checkout
`Manifest.toml`, build timestamp); `start_webui.jl` checks it before every
start, the Sysimage page shows the same verdict, nothing is rebuilt
silently.

| Trigger | What happens |
|---|---|
| image or `sysimage_meta.toml` missing | `no sysimage found`; the start asks whether to build |
| another Julia version | `the sysimage was built for Julia X, this is Y`; the image is never started. Unless a new one is built now, the launcher deletes the image and its metadata file (only these two managed paths) and says `Sysimage is out of date and was removed` |
| `Manifest.toml` hash differs (package update) | `the sysimage does not match the current Manifest.toml`; deleted as above |
| a `src/*.jl` file newer than the image (development checkout) | `the sysimage is older than <checkout>/src`; deleted as above. Released installations never hit this, package updates change the manifest; in a checkout it means a rebuild after every edit, or `--no-sysimage` while working |
| metadata unreadable | treated as outdated |
| failed build | the failed image is removed; the previous one was never replaced, the build compiles to a staging file and moves it into place only when complete |
| Windows, image file kept open by another Julia | cannot be deleted; the start goes on without an image and the next start tries again |
| image given from outside with `-J` | neither checked, rebuilt nor removed |

## Environment

| Variable or flag | Effect |
|---|---|
| `SPARLECTRA_PRECOMPILE_WORKLOAD` = `core` or `full` | set before the installation. Unset, the package precompiles only its module. `core` warms MATPOWER and SCF import, the rectangular Newton-Raphson solve with losses and `run_sparlectra` on a small network; `full` also warms PGM import, state estimation, APSLF and the hybrid start, DC power flow, the service layer and the tap control loop, with a much longer precompile. The sysimage build uses `full` |
| `SPARLECTRA_NO_SYSIMAGE=1` | the launchers skip the image and the question unconditionally (same as `--no-sysimage`, for callers that cannot pass arguments); rules the image out when debugging suspected invalidation or precompilation issues |
| `SPARLECTRA_BUILD_SYSIMAGE=1` | the one-line installers build the image ([Web UI](webui.md)) |
| `SPARLECTRA_STARTUP_FILE=yes` | the launchers and the relaunch through the image read `startup.jl` again (default `--startup-file=no`), for a session where Revise in the server process is wanted |
| `SPARLECTRA_SYSIMAGE_TRACE_TESTS=1` | `tools/build_sysimage.jl` also traces the test suite, for a comparison only |
| `--rebuild-sysimage` | build a fresh image even if the current one is fine |
| `--no-sysimage` | this start ignores the image and the question |
| `--verbose` | `tools/build_sysimage.jl` streams the whole build to the console |

## Where the image lives

`<user root>/sysimage/sparlectra.<ext>` (`so` on Linux, `dylib` on macOS,
`dll` on Windows) next to `sysimage_meta.toml`; the user root also holds
runs, logs, configuration and the MATPOWER cache:

| Platform | Location |
|---|---|
| Linux | `$XDG_STATE_HOME/sparlectra/webui/sysimage/` (default `~/.local/state/...`) |
| macOS | `~/Library/Application Support/Sparlectra/WebUI/sysimage/` |
| Windows | `%LOCALAPPDATA%\Sparlectra\WebUI\sysimage\` |
| Flatpak-packaged IDE (VSCodium, for example) | its own `XDG_STATE_HOME`, so its terminal uses `~/.var/app/<app-id>/.local/state/...` and a plain terminal `~/.local/state/...`: two roots whose configurations, case caches, run histories and images drift apart. Fix: one shared root via the symlink below (merge content you want to keep first; the run registry resolves symlinks, so runs from either environment stay visible in both), or copy `sparlectra.<ext>` plus `sysimage_meta.toml` between the roots, or rebuild where you start from |

```bash
mv ~/.var/app/<app-id>/.local/state/sparlectra ~/.var/app/<app-id>/.local/state/sparlectra.backup
ln -s ~/.local/state/sparlectra ~/.var/app/<app-id>/.local/state/sparlectra
```

## Standalone executable

A relocatable executable with an embedded Julia runtime (for a machine
without Julia) is a checkout tool, not part of the package and not offered
in the Web UI. Plan for 20 to 40 minutes and an app folder of 0.5 to 1 GB;
no cross-compilation: build on the target operating system, rebuild after
updating Sparlectra.

```bash
julia --project=app tools/build_app.jl --flavor=full      # Web UI plus command line
julia --project=app tools/build_app.jl --flavor=runtime   # library runtime, no GUI
julia --project=app tools/build_app.jl --script=my_workflow.jl
# bin/sparlectra exposes the commands webui, run, se, n1 and version, takes the
# library YAML through --config/--set; `bin/sparlectra help` prints its overview
```

!!! details "Why it is built this way"
    **PrecompileTools versus sysimage.** A PrecompileTools workload runs at
    every install and update, and a virus scanner on the compile cache
    makes that write slow: the default warms nothing beyond the module,
    `SPARLECTRA_PRECOMPILE_WORKLOAD` is an opt-in, and the sysimage
    (dependencies included, loaded once at process start) is built on
    demand instead. Anything loaded after the image invalidates
    precompiled methods inside it, recompiled silently on the first run:
    hence AnalyticLoadFlow is compiled in and loaded by `start_webui.jl`
    at startup, and the launchers pass `--startup-file=no`, since a
    personal `startup.jl` usually loads `Revise`, which invalidates a
    large part of the image in the launcher and again in the relaunched
    child.

    **Staleness by source time.** The metadata pins the dependencies, not
    the source, so an image older than a source file would silently serve
    old code: disabled means slower but never wrong. The image is replaced
    by `mv` because a running process keeps its mapped file, so a refresh
    cannot break the session that started it.
