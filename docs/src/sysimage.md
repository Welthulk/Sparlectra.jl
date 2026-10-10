# Sysimage

A fresh Julia process compiles Sparlectra and every solver path on first
use; a PackageCompiler sysimage (Sparlectra and its dependencies as one
ahead-of-time compiled library, loaded with `-J`) removes that wait. The
Web UI start offers to build the image; the **Sysimage** page refreshes
it from the browser.

## Build

`start_webui.jl` (run by the start scripts) decides through
`tools/sysimage_launcher.jl`: a usable image restarts the process with
`-J <image>` (one `Sysimage: ...` line); a missing or outdated image is
named and the start asks `Build the sysimage now? [y/N]`. `y` or
`--rebuild-sysimage` builds; no answer (Enter, 30 seconds, no terminal)
means no, and every path then compiles on first use except the start
page, which the application's package image carries compiled
(`Preparing the first page ... ready in N s`). The build migrates the Web
UI configuration first (`refresh_sparlectra_config_file`, timestamped
backup; a duplicate YAML key stops it, start with
`SPARLECTRA_NO_SYSIMAGE=1` until fixed), then a child process installs
PackageCompiler into the shared `@sparlectra-sysimage-build` environment
and compiles Sparlectra, SparlectraApp and AnalyticLoadFlow: minutes, a
few GB of RAM, half a GB on disk, once per install or update.

| Entry | Build call |
|---|---|
| Command line | `buildSysimage()` with `SparlectraApp` loaded (`dry_run = true` reports the target paths only), or `julia --project=app tools/build_sysimage.jl` in a checkout (workload `tools/sysimage_workload.jl`) |
| Web UI | the **Sysimage** page (`/webui/sysimage`, linked from the Info panel) shows the verdict of the start and runs the same build in the background |
| One-line installers | build the image only with `SPARLECTRA_BUILD_SYSIMAGE=1` |

```bash
./start_webui.sh --rebuild-sysimage   # build a fresh image even if the current one is fine
./start_webui.sh --no-sysimage        # this start ignores the image and the question
```

The build prints one progress line (`[3/4] compiling the system image
04:21`); the full log is `sysimage_build.log` next to the image, and a
failed build leaves the previous image untouched. A refresh from the Web
UI moves the new image into place; a Web UI running on the old one keeps
working on the old code until restarted. A REPL session (or
`julia --project=app`) has no image; the page shows the line to use (the
start script, or `julia -J <image> --startup-file=no --project=<app>`).
Your own scripts take the image with `--sysimage`:

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
compiles on its first click; a failed step is a gap in the trace, not a
failed build.

| Step | What is traced |
|---|---|
| MATPOWER | power flow and state estimation on `warmup_casePST`; one `run_sparlectra` call |
| SCF | power flow, state estimation and export on the shipped demo case |
| DTF | `FOR001.DAT`: power flow and state estimation, reader and net builder |
| CGMES | MiniGrid (fetched once if missing): power flow, state estimation, short circuit |
| N-1 | one contingency run |
| results | network losses, result and diagnostics printers |
| Web UI | the start, its page renders and a clean shutdown |
| APSLF | the solver and the APSLF-seeded rectangular start |
| test suite | not traced by default; `SPARLECTRA_SYSIMAGE_TRACE_TESTS=1` adds it for a comparison |

## Staleness

`sysimage_meta.toml` next to the image records Sparlectra and Julia
version, OS and architecture, the `Manifest.toml` SHA-256 and the build
time; the launcher checks it before every start, the Sysimage page shows
the same verdict, nothing is rebuilt silently.

| Trigger | What happens |
|---|---|
| image or `sysimage_meta.toml` missing | `no sysimage found`; the start asks whether to build |
| another Julia version | `the sysimage was built for Julia X, this is Y`; unless rebuilt now, the launcher deletes image and metadata (`Sysimage is out of date and was removed`) |
| `Manifest.toml` hash differs (package update) | `the sysimage does not match the current Manifest.toml`; deleted as above |
| a `src/*.jl` file newer than the image (development checkout) | `the sysimage is older than <checkout>/src`; deleted as above (rebuild after every edit, or `--no-sysimage` while working) |
| metadata unreadable | treated as outdated |
| failed build | the previous image stays: the build compiles to a staging file and moves it into place only when complete |
| Windows, image file kept open by another Julia | the start goes on without an image and the next start tries again |
| image given from outside with `-J` | neither checked, rebuilt nor removed |

## Environment

| Variable or flag | Effect |
|---|---|
| `SPARLECTRA_PRECOMPILE_WORKLOAD` = `core` or `full` | set before the installation (unset: only the module is precompiled). `core` warms MATPOWER and SCF import, the rectangular solve and `run_sparlectra` on a small network; `full` also PGM import, state estimation, APSLF and the hybrid start, DC power flow, the service layer and the tap control loop (much longer precompile; the sysimage build uses it) |
| `SPARLECTRA_NO_SYSIMAGE=1` | the launchers skip the image and the question (same as `--no-sysimage`, for callers that cannot pass arguments) |
| `SPARLECTRA_BUILD_SYSIMAGE=1` | the one-line installers build the image ([Web UI](webui.md)) |
| `SPARLECTRA_STARTUP_FILE=yes` | the launchers and the relaunch read `startup.jl` again (default `--startup-file=no`), for Revise in the server process |
| `SPARLECTRA_SYSIMAGE_TRACE_TESTS=1` | `tools/build_sysimage.jl` also traces the test suite |
| `--verbose` | `tools/build_sysimage.jl` streams the whole build to the console (otherwise one progress line) |

## Where the image lives

`<user root>/sysimage/sparlectra.<ext>` (`so` on Linux, `dylib` on macOS,
`dll` on Windows) next to `sysimage_meta.toml`; the user root also holds
runs, logs, configuration and the MATPOWER cache:

| Platform | Location |
|---|---|
| Linux | `$XDG_STATE_HOME/sparlectra/webui/sysimage/` (default `~/.local/state/...`) |
| macOS | `~/Library/Application Support/Sparlectra/WebUI/sysimage/` |
| Windows | `%LOCALAPPDATA%\Sparlectra\WebUI\sysimage\` |
| Flatpak-packaged IDE (VSCodium, for example) | its own `XDG_STATE_HOME` (`~/.var/app/<app-id>/.local/state/...`): IDE terminal and plain terminal use two roots that drift apart. Fix: one shared root via the symlink below (merge first), or copy image plus `sysimage_meta.toml` between the roots, or rebuild where you start from |

```bash
mv ~/.var/app/<app-id>/.local/state/sparlectra ~/.var/app/<app-id>/.local/state/sparlectra.backup
ln -s ~/.local/state/sparlectra ~/.var/app/<app-id>/.local/state/sparlectra
```

## Standalone executable

A relocatable executable with an embedded Julia runtime (for a machine
without Julia) is a checkout tool, not offered in the Web UI: 20 to 40
minutes, an app folder of 0.5 to 1 GB, built on the target operating
system (no cross-compilation), rebuilt after updating Sparlectra.

```bash
julia --project=app tools/build_app.jl --flavor=full      # Web UI plus command line
julia --project=app tools/build_app.jl --flavor=runtime   # library runtime, no GUI
julia --project=app tools/build_app.jl --script=my_workflow.jl
# bin/sparlectra exposes the commands webui, run, se, n1, export and version, takes
# the library YAML through --config/--set; `bin/sparlectra help` prints its overview
```

!!! details "Why it is built this way"
        A PrecompileTools workload runs at every install and update, and a
    virus scanner on the compile cache makes that slow, so the default
    warms only the module and the sysimage is built on demand. Anything
    loaded after the image invalidates methods inside it, hence
    AnalyticLoadFlow is compiled in and `startup.jl` (usually `Revise`) is
    skipped. An image older than a source file is disabled rather than
    trusted. `mv` replaces the image because a running process keeps its
    mapped file.
