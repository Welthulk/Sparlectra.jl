# Fast Start (Sysimage)

Starting the Web UI pays Julia's just-in-time compilation cost twice: once
for loading the package and once for the first solves. Sparlectra ships two
complementary answers:

1. **PrecompileTools workload (always on).** The package precompiles the
   solver hot path at install time (MATPOWER parse, rectangular
   Newton-Raphson with UMFPACK, losses, DC power flow; issue #288). This
   needs no setup and already covers most of the solver latency.
2. **PackageCompiler sysimage (optional, this page).** A system image bakes
   Sparlectra and its dependencies into one ahead-of-time compiled shared
   library that Julia loads via `-J`. With it, a Web UI start reaches a
   served PowerFlow page in a few seconds; the remaining JIT time
   disappears.

## Building the image

From a checkout, either:

- **Script:** `julia tools/build_sysimage.jl`. The build runs in the shared
  environment `@sparlectra-sysimage-build` (PackageCompiler is installed
  there, never into the package `Project.toml`), executes the workload in
  `tools/sysimage_workload.jl` (Web UI start with the reserved warm-up, one
  `case14` run, clean shutdown; no run-history entries), and writes the
  image below the Web UI user root.
- **Web UI:** the **Fast start** page (navigation bar) shows the current
  image state and offers a **Build fast-start image** button. The build
  runs as a separate Julia process; the Web UI stays usable, the progress
  lands in the `sysimage_build.log` artifact next to the image, and the
  operation log records `sysimage_build_started` and
  `sysimage_build_finished` or `sysimage_build_failed`. One build at a
  time; the trigger is POST-only and loopback-only like every other route.

Expect a build to take roughly 7 to 20 minutes and an image of a few
hundred MB (reference measurement on a Linux workstation, Julia 1.12:
7.1 minutes build time, 537 MB image; `start_webui.sh` reached the served
PowerFlow page after 15 s with the image versus 45 s without). The image
takes effect on the **next start**.

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

## Staleness rules

`sysimage_meta.toml` is the validity contract: Sparlectra version, Julia
version, OS and architecture, the SHA-256 of the checkout `Manifest.toml`,
and the build timestamp. `start_webui.sh` / `start_webui.bat` compare the
Julia version string and the manifest hash with plain shell tools before
every start:

- **Match:** the launcher appends `-J <image>` and prints one
  `Fast start: using sysimage ...` line.
- **Mismatch or missing image:** the launcher prints one warning
  (`sysimage stale, starting without it; rebuild via
  tools/build_sysimage.jl`) and starts normally. A stale image never
  breaks a start.

Package updates (a changed `Manifest.toml`) and Julia updates therefore
disable the image until the next rebuild; the Web UI page shows the exact
mismatch reasons. There is no automatic rebuild by design: detect and warn
only.

## Escape hatch

`SPARLECTRA_NO_SYSIMAGE=1` makes the launchers skip the image
unconditionally. Use it to rule the sysimage out when debugging suspected
invalidation or precompilation issues:

```bash
SPARLECTRA_NO_SYSIMAGE=1 ./start_webui.sh
```
