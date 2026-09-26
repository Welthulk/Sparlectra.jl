# Local PowerFlow Web UI

The Web UI is a local browser interface on top of the
[PowerFlow service](powerflow_service.md). It holds presentation, form
parsing and route handling only; every calculation runs through
`start_powerflow_run` and `run_sparlectra_api`. It binds to loopback
(`127.0.0.1`, `localhost`, `::1`) and has no authentication: a single-user
tool on the machine that computes; a wider binding would expose the run
directory and the configuration editor.

## Start

### One-line install

Installs Julia if missing (via juliaup), downloads the latest tagged release
into `Sparlectra/` in the current directory, and starts the Web UI:

```sh
curl -fsSL https://raw.githubusercontent.com/Welthulk/Sparlectra.jl/main/tools/install_webui.sh | sh
```

```powershell
iwr -useb https://raw.githubusercontent.com/Welthulk/Sparlectra.jl/main/tools/install_webui.bat -OutFile install_webui.bat; .\install_webui.bat
```

The same scripts, `tools/install_webui.sh` / `tools/install_webui.bat`, ask
three questions:

- **Update** when an existing copy is older than the latest release; the old
  copy is kept as `Sparlectra.old`.
- **Sysimage**: the installer leaves the [sysimage](sysimage.md) build to the
  Web UI start and builds it only with `SPARLECTRA_BUILD_SYSIMAGE=1`.
- **Desktop shortcut** for restarting the Web UI: Windows `Sparlectra Web
  UI.lnk`; Linux an application-menu entry plus a `.desktop` file on the
  desktop (GNOME needs a one-time right-click "Allow Launching"); macOS a
  desktop symlink.

Unattended installs answer them with `SPARLECTRA_UPDATE=1/0`,
`SPARLECTRA_BUILD_SYSIMAGE=1/0` and `SPARLECTRA_CREATE_SHORTCUT=1/0`.

### Start script

From a checkout or an installed copy:

```sh
julia --project=. start_webui.jl
```

`start_webui.sh` / `start_webui.bat` run the same script (and point at the
install script when Julia is missing); a Windows desktop shortcut is
right-click `start_webui.bat`, *Send to > Desktop (create shortcut)*.
`start_webui.jl` delegates to `start_sparlectra_webui`. It first compares
the direct dependencies in `Project.toml` with `Manifest.toml`, for the
library and for `app/` (four TOML reads, no package load); a missing
dependency or manifest is resolved, instantiated and compiled once, with the
time reported, and the next start compiles nothing. `Manifest.toml` is not
tracked, so a fresh clone and an old manifest are both covered.
`julia --project=. start_webui.jl --env-only` does only that first-start
work.

Then it looks for a usable sysimage and offers to build one; the build,
its flags, the staleness rules and the **Sysimage** page are on
[Sysimage](sysimage.md).

### From the Julia REPL

```julia
using Pkg
Pkg.activate("path/to/Sparlectra/app")
Pkg.instantiate()   # first time only
using SparlectraApp

server = start_sparlectra_webui(open_browser = true)
wait(server.task)
```

Such a session runs without the sysimage; how to start a REPL session on
the image is on [Sysimage](sysimage.md).

`start_sparlectra_webui` returns a `SparlectraWebUIServer` handle and serves
as soon as the socket is bound, with no hidden warm-up run.
`open_browser = true` opens an app-style window (Microsoft Edge, Google
Chrome, Chromium or Brave); without one of them the URL
`http://127.0.0.1:8080/powerflow` is logged. If the port is occupied, stop
the old Julia process or pass another `port`.

### Directories and files

Results go to `%LOCALAPPDATA%\Sparlectra\WebUI\runs` (Windows),
`$XDG_STATE_HOME/sparlectra/webui/runs`, default
`~/.local/state/sparlectra/webui/runs` (Linux) or
`~/Library/Application Support/Sparlectra/WebUI/runs` (macOS), created
automatically; the operation log is in the sibling `logs`, downloaded and
generated cases in the sibling `data/mpower`. That case directory is shared
with `ensure_casefile` and the test suite; `SPARLECTRA_LARGE_CASES_DIR`
overrides it for all three. The bundled `warmup_case*.jl` precompile
workloads are copied there on first start; the reserved `warmup_` prefix
hides them from the case selector (the sysimage workload runs the shipped
`sp_case5`/`sp_case60`).

The first start copies the package configuration template to the
user-writable `config/configuration.yaml`; `output_root = "..."` and
`config_file = "..."` override the defaults, and an explicit configuration
file is never overwritten. Next to it the Web UI keeps:

- `configuration.template.yaml`, the template the file was last aligned
  with: a key still at an old template default follows a new default at
  start (reported with key, old and new value; previous file kept as
  `configuration.yaml.template-follow.bak`), a changed value stays;
- `configuration.yaml.user-keys.txt`, the keys the Web UI itself saved
  (settings save, notice dismiss, configuration editor), which are never
  followed; deleting it or the template copy only loses that memory;
- `configuration.yaml.settings-save.bak`, the backup of a settings save.

The pages show the effective configuration, output root, case cache and
operation log; the browser cannot change the output root. Header and footer
show the logo (`docs/src/assets/logo.png`) and `Sparlectra.version()`;
header and operation log also name the package path and Git commit. The
**Info** panel names the serving build (`sysimage, built ...` or `native
session`) and links to the Sysimage page.

## Importing cases

The **Case** page carries the case chooser, upload, export and the
per-format import options.

### [Case selector](@id webui-case-selector)

The selector lists the MATPOWER `.m` files and runnable DTF `.DAT`
candidates of the case directory, the shipped demo cases `sp_case5` to
`sp_case188` ([Shipped Demo Cases](demo_cases.md)), staged into the cache
with their sidecars (per-case configuration, measurement CSVs) on first
use, the three CGMES deliveries exported by Sparlectra itself
(`data/cgmes_demo`: `sp_case14`, `sp_case118`, `sp_casePST`, four profile
files each) as `<case>_cgmes.zip`, and PowSyBl table bundles
(`<case>.powsybl` directories) and IIDM files (`.xiidm`). Generated `.jl` cache files and
`warmup_` files are hidden. FOR002-like `.DAT` files are not primary cases;
they go into the optional **FOR002 reference** field (absolute path, path
in the case cache, or an offered candidate), used only for legacy reference
comparison.

**Or type case file path** overrides the selector: a bare name such as
`case14.m`, `case118.m` or `case9241pegase.m` is downloaded into the case
directory, an existing absolute or relative `.m`/`.jl` path is used as
given; path-like missing inputs and URLs are rejected. A generated `.jl` cache file resolves back to its `.m` source
(the bypass is recorded) and is rejected without one: large `.jl` literal
cases can fail while Julia loads them, so `.m` stays the canonical source.

**Case input format** defaults to **Auto**: MATPOWER files, CGMES
deliveries (folders and ZIPs), PowSyBl bundles and IIDM files, and, where
the FOR001 markers are unambiguous, native DTF input. **CGMES (ENTSO-E,
folder or ZIP)** forces the CGMES importer and is preselected for a `.zip`
or a directory; **PowSyBl (IIDM file or .powsybl bundle)** is preselected
for those files and shows the PowSyBl import options (HVDC mode, remote
regulation, slack choice, base MVA; see
[PowSyBl Import](powsybl_import.md)). A `.xiidm` file is read in Julia,
no Python is involved.
**Sparlectra Case Format (.scf.json)** and **power-grid-model JSON
(input.json)** name the same reader (the `sparlectra` block is optional);
Auto already sends every `.json` there, so they matter only when the
extension does not say what the file is. The native DTF path is
experimental, for diagnostics and validation; its **Selected DTF outage
labels/indices** field takes one label or index at a time, and the result
page shows a compact outage summary while the rows stay in artifacts.

The format-bound option sections derive from the adapters' option structs,
so the form cannot drift from what a conversion accepts; the Basic set (keys
used in a workshop, or set in the example configuration without a default)
renders directly, the rest and the DTF diagnostics on the Runs page fold
into **Advanced**. Deleting a case removes all its companion files
(configuration, weights, measurement CSVs).

A collapsible MATPOWER acknowledgement next to the case inputs links to the
project, its citation guidance and the 2011 paper DOI, and notes that
ACTIVSg, PEGASE, RTE and other cases ask for additional citations in their
headers.

### [Import case files](@id webui-import-case-files)

**Import case files** opens the native file picker and accepts several files
at once: MATPOWER `.m`/`.M`, DTF `.dat`/`.DAT`, CGMES `.zip` deliveries,
CGMES profile files (`.xml`) and PowSyBl IIDM files (`.xiidm`); the
server validates the extension again. An uploaded `.xiidm` file runs as
it is. A `.powsybl` bundle is a directory and is copied into the case
directory by hand.
Several `.xml` files (the EQ, SSH, TP and SV profiles of one delivery, such
as a `data/cgmes_demo` folder) are packed into one `<stem>_cgmes.zip` named
after the common stem; the set must contain the EQ profile. A CGMES ZIP may
hold nested ZIPs; the importer opens them in memory. A Sparlectra Case
Format `.json` file is parsed and validated before it is stored.

Importing is copy-only: no run, no run directory, no parsing of `.m` code.
Files land in the case directory the form shows (`data/mpower`, in a
checkout too), the selector is rebuilt from disk, and the first runnable
file may be preselected; the run still needs **Start PowerFlow run**.
Limits: 100 MiB per file, 250 MiB per request. An existing file is never
overwritten (`already exists`, the other files still import), and names
that would resolve outside the case directory are refused. Enter in **Or
type case file path** resolves the value the same copy-only way: a bare
case name via [`ensure_casefile`](@ref), `cgmes:<alias>` as below, or a full
local path copied with the same validation; the file then appears in the
**Existing case file** selector.

Re-importing a saved set brings the `.config.yaml` sidecar along with the
`.scf.json` and its measurement set(s) when all are selected together;
the `.yaml` is validated as a case-scope configuration file
(`scope: case`) before it is stored, and rejected otherwise.

### Export and download

The export row writes the open case as `.scf.json` or plain PGM `.pgm.json`,
replacing a format suffix the case already carries (`case14.scf.json` as
PGM gives `case14.pgm.json`); both products appear in the selector as
runnable cases. **Download selected case** serves what the selector shows,
after an export also the file just written: only plain files inside the
case directory (bare name, or an absolute path pointing into it), compared
by real path because the state directory can sit behind a symlink (Flatpak);
anything outside, and a CGMES delivery directory, is refused.

### Save case as

**Save case as**, next to the export buttons (also linked from the State
Estimation section), saves the current case under a new name instead of
exporting in place: one network with several variants, one file per
variant. It writes three files into the case directory:

- `<name>.scf.json` via [`exportSCF`](@ref), with `meta.case_name = <name>`
  and `meta.source_reference` noting the source case; a MATPOWER or DTF
  source is saved as SCF, the original file is never touched.
- `<name>.config.yaml`: the source case's saved case-scope settings merged
  with any unsaved change made on the page; the source case's sidecar file
  is never modified.
- `<name>.measurements.csv` for each measurement set bound to the source
  case, with its `# case:` header rewritten to the new name; a second bound
  set is copied as `<name>_<original stem>.measurements.csv`.

An existing `<name>.scf.json` is refused unless "overwrite" is ticked. The
optional "start from the solved state" checkbox solves the case once with
the effective settings and writes the solved voltages as the copy's
`sparlectra.start_state`. Installation-scope settings (`output.*`,
`benchmark.*`, `runtime.*`, ...) are never written into the case
configuration file. The file layout is described in
[Sparlectra Case Format](scf.md).

### CGMES test configurations

`cgmes:<alias>` in **Or type case file path** downloads the official ENTSO-E
test-configuration package once (about 22 MB) into the CGMES cache
(`data/CGMES`, overridable with `SPARLECTRA_CGMES_CACHE`) and packs the
requested configuration with its boundary set into `cgmes_<alias>.zip` in
the case directory. Aliases: `microgrid_be`, `microgrid_nl`,
`microgrid_assembled`, `smallgrid`, `smallgrid_nb`, `fullgrid`,
`fullgrid_nb`, `realgrid`. Repeated requests reuse the ZIP. If the download
fails, the message names the path where the package can be placed by hand.
The test data is never committed.

### Import analysis on upload

Every uploaded CGMES delivery is checked at once: a complete one (or one
whose boundary is found or supplied automatically) reports "ready to
compute", an incomplete one gets the full import analysis: the message names
the missing declared `md:Model.DependentOn` dependencies by model id
(typically the boundary set), and the report (supplied models, dependency
matching, unresolved-reference histogram, verdict) is written next to the
case as `<case>.import_analysis.txt`. The same analysis is appended to
`cgmes.log` whenever a run's CGMES import aborts. Scripted use: the run mode
`import_analysis_mode` of `start_powerflow_run` (a non-importable delivery
ends as a failed run with reason `import_analysis_not_importable`) or
`analyzeCGMES`.

**Require boundary set** (CGMES section, saved per case) submits
`cgmes_import.require_boundary`; unchecked, an incomplete delivery imports
where possible, but buses without a resolvable `BaseVoltage` still abort,
with the analysis explaining why.

### Long actions on a case

Generating a measurement set or adding noise runs inside the request and
can take a while on a large case (with **critical measurements** the
generator re-reads the criticality after every removed row). Meanwhile the
case is busy: a second generator action or a run is refused with a message
and the buttons show a spinner; a queued or running run blocks the
generator the same way.

## Starting a run

Run-independent options live on the **Settings** page (solver, start,
output and expert options, plus the configuration block: file display,
check and refresh, saved case settings, the ignore switch) and reach a run
through the configuration precedence (`resolve_config`): saving with target
"this case" writes case-scope keys into the case's configuration file,
target "configuration file" merges into the general YAML (backup
`.settings-save.bak`). The **Runs** page keeps the run parameters (modes,
N-1 kind, scenario source, screening, DTF outage fields, benchmark trigger),
the state-estimation section (estimator options and measurement-set
selection submitted with `se_mode`; measurement generator, set details with
inline editors and measurement upload as fold-out tabs) and the N-1 editors
as tabs: the scenario editor embeds directly (SCF cases), the weights
editor loads on first open (seeding its element names builds the net once).
The SE section uses the shared case selection.

The measurement generator's fields are the generator options of
[State Estimation Measurements](state_estimation_measurements.md#Measurement-generator-v2) (sigmas per
class, noise, gross errors, tap deviations, truth source, flow ends, passive
nodes as zero-injection constraints, **critical measurements**); the
removed and critical rows are named in the set's comment lines and the
estimation run log.

### [Form options](@id webui-form-options)

The settings and run forms offer:

- case name or local path, configuration template file, the read-only
  output root (set only through `start_sparlectra_webui(; output_root)`);
- tolerance (decimals or scientific notation such as `1e-8`; the exponent
  spinner, buttons or Up/Down keys, steps the exponent and keeps the
  mantissa: `1e-5` down gives `1e-6`, up gives `1e-4`), maximum iterations,
  autodamping with its minimum factor;
- **Q-limit handling**: the `power_flow.qlimits.enabled` checkbox and the
  enforcement-mode selector in one block; the selector also offers **off**
  and shows it whenever the handling is off;
- **Solver**: one radio group `power_flow_solver` with the peer values
  `rectangular` (AC Newton-Raphson), `apslf` (AnalyticLoadFlow) and `dc`
  (linear screening model), writing `power_flow.solver`, `power_flow.apslf.*`,
  `power_flow.apslf_start.*` and `power_flow.dc.*`
  ([solver selection](powerflow_configuration.md#solver-selection-rectangular-vs-apslf),
  [DC power flow](powerflow_configuration.md#solver-selection-dc-power-flow));
  radio buttons are never disabled, so the choice is always submitted;
- wrong-branch detection mode, start angle and voltage modes;
- **Advanced start values**: the current-iteration pre-solve fields
  (`power_flow.start_current_iteration.*`,
  [guarded pre-solve](powerflow_configuration.md#guarded-current-iteration-start-pre-solve)),
  a preconditioner between start projection and Newton-Raphson, not a new
  start mode;
- **Merit-function line search**: `power_flow.merit.enabled`,
  `power_flow.merit.armijo_c1`, `power_flow.merit.fallback_max_mismatch`
  ([options](powerflow_configuration.md#merit-function-line-search-options)),
  requires `power_flow.autodamp = true`; `scale_p`/`scale_q`/`scale_v` are
  YAML-only; a diagnostic run writes `merit_linesearch.log`;
- **Non-convergence handling** (Advanced): `power_flow.rescue` (retry
  ladder) and `power_flow.dc.fallback` (standalone DC result when AC has no
  solution), see [Power-Flow Configuration](powerflow_configuration.md);
- MATPOWER import conventions: auto-profile mode (`off`, `recommend`,
  `apply`) plus transformer ratio, phase shift, bus shunt, PV voltage source
  and comparison reference;
- logfile output mode, performance timing, post-run diagnostics, benchmark
  enablement, samples and seconds.

Only keys in `GUI_EDITABLE_CONFIG_KEYS` are submitted, the selected template
is never modified, and each run writes `effective_config.yaml`. Every control
carries a hover text with its operation and effect; a **?** next to a
control with something technical behind it (solver options, estimator
diagnostics, Q-limits, formats) opens its help page in a new tab: the hint,
the matching section of this documentation as shipped with the running
version (the beginning of a long section), and a link to the section on
the published site. `webui.docs_base_url` (configuration file) points that
link at a local docs build or a pinned version. **Back** returns through
the browser history with the entered values intact (fallback `/powerflow`).

Fields that do not apply to the selected solver or step control
(autodamp/trust-region) are grayed out in place, not hidden; a grayed field
is `disabled` and not submitted.

### AC, APSLF and DC

`AC (Newton-Raphson, rectangular)` reveals the start-value block (**Flat
start**, **Use APSLF start values** with its order field, **Use DC start
values**) and every AC-only option: tolerance, autodamping/merit/trust-region
step control, Q-limit handling, maximum iterations, wrong-branch detection,
start angle/voltage mode, the current-iteration pre-solve and the
transformer tap-changer model. By default
(`power_flow.start_mode.angle_mode = dc`) Newton-Raphson seeds its start
angles from an internal DC pre-solve, controlled by **Start angle mode**;
that pre-solve is unrelated to the standalone `DC` solver.

`APSLF (AnalyticLoadFlow)` reveals the **APSLF solver options** (highest
coefficient/order, Padé evaluation, NR polish). **Use APSLF start values**
in the Newton-Raphson block seeds the rectangular solver instead and is
mutually exclusive with the APSLF solver (the configuration rejects both).
An unchecked checkbox (APSLF, Q-limit, current-iteration) submits `false`
rather than omitting the key. The APSLF solver
needs the optional AnalyticLoadFlow.jl package in the server process;
without it the run fails with a pre-solve error instead of falling back to
the rectangular solver.

`DC` replaces Newton-Raphson with the linear model
([DC power flow](powerflow_configuration.md#solver-selection-dc-power-flow)),
the `dc` value of `power_flow.solver`, and grays out every AC-only option
above. A DC run uses the same job, abort, status, history and artifact
machinery; the result page marks it with a **DC solution** badge next to
**Solver** and the history column shows `dc`, so the implicit
`Vm = 1.0 pu` and lossless flows are never taken for an AC result;
`iterations` is `1` (a direct solve).

The **Solver** column names the method that ran: the power-flow solver for
a plain power flow, one started from an estimate and every contingency
case; `wls` for a state estimation; `iec60909` for a short circuit; empty
for an import analysis (derived from the run kind when a run carries no
method).

#### [Flat start](@id webui-flat-start)

`power_flow.flatstart` starts the Newton-Raphson solve at 1.0 pu and 0
degrees on every bus and ignores the imported start voltages (MATPOWER
`VM`/`VA`, CGMES `SvVoltage`, SCF `start_state`). Off, the imported values
seed the solve; a delivery is built around its own operating point, so that
is the better start for real networks and the default.

This checkbox is the one start switch: while it is on, a run switches
**Use APSLF start values**, **Use DC start values**, the current-iteration
pre-solve and the start projection off and treats both start modes as
`classic`, so the start really is flat and no projection or pre-solve moves
it; the run log names what was switched off. The greyed controls keep their
saved values and come back when the flat start is unchecked.

A CGMES run honours the flat start under **CGMES start values** = `auto`;
an explicit `sv` on the Case page still starts from the delivery state
([Power flow configuration](powerflow_configuration.md#pf-solver-core)).

### Case-specific settings

A successful result page offers **Save settings for this case**: only the
form options of that run plus traceability metadata go into a sanitized
YAML profile below the output root (not `effective_config.yaml`, solver
internals, artifact paths or convergence diagnostics). A non-converged run
shows **Save these settings anyway** and records that override; nothing is
saved automatically. Reopening the case prefills the form from the profile
with a notice: defaults first, the global configuration untouched, only
editable fields prefilled, and a manual edit wins for the submitted run.

With a case selected, every page (Case, Runs, Settings) shows the values a
run of that case will use: case settings over the configuration file over
the defaults. Without a case, the configuration file. **Show the
configuration file's values** (`?config_view=1`) switches the Settings page
to the file alone, **Show the case settings** back. A save with target
"this case" writes what the form shows; a save to the configuration file
says how many of the saved keys the selected case overrides. Between configuration file and saved profile the last edit wins: a
newer YAML file takes precedence for the keys it sets on the next page load
(the notice says so), other fields keep their saved values; a page refresh
is enough.

With a Sparlectra Case Format file selected, the form is seeded from the
case file: it posts a value for every option it shows (explicit
overrides, the highest precedence level), so selecting a case file moves
its controls to the file's values, and the page says which settings came
from the case. A saved settings profile for the case still wins over the
file, and any control the user edits wins over both.

### Scenarios and N-1

The action row of the run form carries the N-1 controls: outage kind
(branch/generator), **scenario source** (generated N-1 lists, or the
scenarios block of an SCF case), **screening mode** (configured / off /
flag / only) with a **margin** field (empty means the configured
`contingency.screening.margin_pct`), and links to the weights editor and,
for SCF cases, the scenario editor. Screening is off by default
([N-1 Contingency Analysis](contingency.md)). The run
summary names the screened count and the case-list source; the result page
tabulates the cases for every format (name, weight, convergence,
iterations, start, voltage envelope, worst loading, severity, islands, shed
load, screened marker), ranked by severity with failures first and capped
at 100, the linked CSV keeping the complete list in input order. A screened
row shows its one-step estimate (no full solve), a flagged row the estimate
next to the full-run values, a failed case its error text.

**Weights.** The "edit N-1 weights" link next to the outage-kind selector
opens a per-case weights editor. Weights live next to the case as
`<case-stem>.contingency-weights.csv`, the two-column format
`readContingencyWeightsCSV` parses; the file is hidden from the case list
and deleted with the case. The editor seeds a table with the case's element
names and offers a raw-CSV text area and a file upload; an upload replaces
the existing file after validation, a malformed CSV is rejected with the
line number. Rows left at `1.0` are omitted on save. A run picks the
weights up whenever the file exists; names that match no element are
reported in `run.log`, never fatal. A weight file applies to case-list runs
(the outage-kind selector and the `n1_*` scenario sources); a scenario run
from a case file's block or an external scenario JSON uses the per-scenario
`weight` of the block and ignores the weight file. A weight only reorders
the severity ranking and never skips a case
([N-1 Contingency Analysis](contingency.md)).

Scenarios need an SCF case; for other formats the form offers the generated
N-1 sources and says "Scenarios need an SCF case; export this case as SCF
first" next to the export action. The **scenario editor** lists the open
case's scenarios (name, weight, op count; edit, duplicate, delete, new) and
edits one as op rows: op (status/set/scale), target class, a component
selector filled from the case as `id: name (from-to)` and filtered by
class, then the fields the op needs (value for status, field plus value for
set, factor for scale). Validation is server side with the scenario rules;
errors name the scenario and op index, and a tap patch on a regulated
transformer is rejected with the controller named. **Save into case file**
writes the scenarios block through the SCF writer; nothing stays in the
browser.

### Diagnose

**Diagnose** next to **Start PowerFlow run** evaluates the mismatch at the
case's own stored operating point (MATPOWER `VM`/`VA`, or the `SvVoltage`
state of a CGMES delivery) without a corrective Newton step. Every
start-value machine is forced off (`flatstart = false`,
`start_projection = false`, `dc_seed_unconditional = false`,
`start_current_iteration.enabled = false`, `apslf_start.enabled = false`,
plus `max_iter = 1` and `qlimits.enabled = false`), so the reported residual
reflects the imported model, not the start guess or the step control. The
forced settings outrank the form and a case configuration file next to the
case. The run uses the normal pipeline (history, artifact viewer, the
`diagnose.log` report) and adds `self_check.log` (forced settings, start
residual and, for CGMES, the count of buses without a usable `SvVoltage`,
which start at the flat `1.0 pu / 0°` fallback), `self_check_residuals.csv`
(per-bus P/Q residuals with SV coverage, transformer-terminal and shunt
counts) and `diagnose_self_check_config.yaml`. Programmatically:
[`run_fixed_reference_self_check`](@ref).

One step practically never converges; the residual is the measurement, so
history and result page label the run **diagnosed** on a neutral badge (raw
status in the tooltip) with a one-line explanation of the single iteration
and non-zero mismatch. A diagnose run that could not run at all keeps the
failure vocabulary.

### Short circuit

**Short circuit** next to Diagnose evaluates the balanced three-phase
currents (IEC 60909-0, [`runShortCircuit!`](@ref)) for every bus of a CGMES
delivery, maximum and minimum case in one run, without a power-flow solve.
The button is enabled only when the delivery carries short-circuit source
data (synchronous machines, feeder short-circuit currents or equivalent
impedances); the server checks the contents and otherwise disables it with a
tooltip. The run writes `short_circuit_max.csv` and `short_circuit_min.csv`
(per-bus `Ik''`, `Sk''`, `κ`, `i_p`, safety flag and reasons), a `run.log`
with the harvested-data coverage report and a summary row; rows with
defaulted or skipped data get a warning badge because a flagged `Ik''max`
is a lower bound. A delivery without usable sources fails with
`short_circuit_data_missing` (coverage report in `run.log`), a non-CGMES
case with `short_circuit_requires_cgmes`.

### State estimation

The state-estimation section of the Runs page (`/powerflow#state-estimation`;
`/stateestimation` redirects there) uses the shared case selection, selects a
measurement set and runs observability (traffic light, structural-island
note), the WLS solve and the diagnostics of
[State Estimation](state_estimation.md). Uploaded `.csv` files are offered
when the content sniff finds the v1 version comment. Preselection: only a set
bound to the selected case is preselected (cases with a bound set are
starred, foreign sets are labeled); a case file with its own measurements
offers those first and says how many it carries, otherwise the page asks for
a set. A set bound to the case an SCF file was exported from stays usable,
because the file records its source.

The set tab offers download, re-upload and an inline editor (small files,
saved atomically, first line must stay the version comment) and renders the
`sparlectra-taps v1` tap table of a generated file, offered for download with
the file. The page names the state of a carried set (ideal or measured);
**Add noise to this set** perturbs each value with the sigma its own row
declares. **Reset saved settings** under the generator deletes the per-case
settings sidecar. **Estimate taps** (on by default) releases every
in-service transformer with a ratio tap changer except machine transformers;
the result page shows the per-transformer table (model step, fixed step and
the change between them; frozen taps carry their reason in the status
column) and the fixation J drop. The form warns when
`k_suppress < k_eliminate`, because suppressed rows then rarely reach the
elimination.

Artifacts (`measurements.csv`, `se_diagnostics.md`, `se_view.md`,
`se_state.csv`, `shunt_estimates.csv`, `se_tap_estimates.csv`,
`se_bad_data.csv`, `se_deltas.csv`) land in the run history with kind `se`.
The result page offers **Run power flow from this estimate** and the
topology hypothesis test behind a button (writes `topology_hypotheses.md`).
The summary shows `J`, the expected value (`E[J] = dof` for healthy noise)
and the Wilson-Hilferty 3-sigma band verdict (with the `:high`/`:low`
reason), always for the state after the elimination: eliminated rows are
out of `J` and `dof`, suppressed rows stay in the honest `J` and surface as
`J_active`.

## Running jobs

A submission starts a background worker and redirects to the status page:
elapsed time as `HH:MM:SS` beside the status, and a details table with
`elapsed_seconds`, requested and resolved case paths, start time and a
manual refresh link. Queued, running and aborting
pages refresh every two seconds (these `autorefresh=1` requests are not
logged as user actions); terminal pages stop refreshing. The start page and
the run history show a banner with **Open status** and **Abort** while a
job is queued or running. A queued, running or aborting run cannot be
deleted; a terminal one can.

For large cases, the result page, `run.log`, `result.json` and
`performance.log` show time per phase (reader, network builder, solver,
artifacts).

### Abort

A queued or running job can be aborted from its status page or from the
banner on the start page. Abort is cooperative: an operation that cannot be
interrupted, such as a sparse factorization, finishes first; the status page
shows the current phase. An aborted run keeps its run directory and is never
reported as success.

If a run stays in `aborting` for more than 60 s, the status page offers
**Hard reset Web UI**. It marks the run invalid and stops the server.
Restart with `julia --project=. start_webui.jl`.

### Feedback

Formulas on the documentation pages are rendered by a bundled copy of
KaTeX (served from the package, no network needed).

One-off messages (validation errors, import summary) open as a dismissible
popup. Recent errors stay listed on the **Last errors** page.

## Results and artifacts

Every finished run gets a result page: run ID, schema version, status,
convergence and solution flags, iterations, final mismatch, reason and
message, input paths, output directory and a **Solver** entry. Artifacts come from
`list_powerflow_artifacts` and are resolved by exact name through
`resolve_powerflow_artifact`, never by joining browser input to a path.
A CSV artifact opens as a table (delimiter taken from the header line), a
Markdown report (`se_diagnostics.md`, `se_view.md`, DTF summaries) rendered;
**Raw text** on the page shows the file as written. Other text artifacts
(JSON, YAML, logs, HTML) open as escaped text in a scrollable panel that
keeps long lines; other files download, and every artifact page offers a
download. **Download all artifacts as ZIP** packs
the exposed artifacts into `sparlectra_run_<run_id>_artifacts.zip`, skipping
missing optional artifacts and unsafe names.

### [Output modes](@id webui-output-modes)

- **Logfile output mode**: `classic` writes the result report plus a
  compact timing and status summary (`solver_time`, `representative_time`,
  iterations, final mismatch, outcome; benchmark median and sample count when
  benchmarking is on); `full` adds **Full run details** with the effective
  typed configuration, artifact choices and status diagnostics. MATPOWER
  `.m` runs also print the "Original/Final effective MATPOWER import
  options" and "MATPOWER auto-profile recommendations" tables, because `.m`
  files leave `shift_sign`, `shift_unit` and `ratio` ambiguous; DTF `.DAT`
  runs never print them, so their shorter `run.log` is expected.
- **Performance timing** (`off`, `compact`, `full`) writes
  `performance.log` with the phases of one request (request parsing, case
  resolution, configuration, case loading and network construction, solve,
  postprocessing, artifact writing, solver and total time; `full` adds
  internal profile entries), unlike `benchmark.enabled`, which measures
  repeated solves and their median.
- **Export case as CGMES delivery (EQ+TP+SSH+SV, ZIP)** writes the case as
  one re-importable delivery into the run's artifacts, for every case format
  and also on non-converged runs; see [CGMES Export](cgmes_export.md).
- **Write bus/branch CSV files** (API default off, because large networks
  produce large files) writes `bus_voltages_complex.csv` (per bus `vm_pu`,
  `va_deg`, numeric `v_re` and `v_im`, readable `v_complex`, nominal and
  actual voltage, generation, load, Q-limit and control columns) and
  `branch_flows.csv` (active and reactive power at both ends, losses,
  rating, status, overload) from the `ACPFlowReport` rows.
- **CSV format** is the machine-scope key `output.csv_format` (**Save
  settings** writes it to the configuration file, a per-case save leaves it
  out) and applies to every CSV of every run type, state estimation
  included: `technical` (default: comma delimiter, decimal point, no
  grouping), `excel_de` (semicolon, decimal comma, thousands dot),
  `excel_us` (comma, decimal point, thousands comma, grouped numbers
  quoted). The Excel formats avoid exponent notation where practical;
  `v_complex` follows the decimal notation, `v_re`/`v_im` stay numeric.
  Excel may still auto-convert identifiers like `1E5`; use **Data > From
  Text/CSV** with text types when that matters.

### Diagnostic artifacts

`diagnose.log` is written by **Diagnose** or with `run_diagnostics = true`,
never by a normal run. On a non-converged run it is a report: "Diagnosis"
names the worst-mismatch bus and equation and classifies the mismatch trend
(monotonic, oscillatory, stagnant, diverging to non-finite) and the autodamp
health; "Branch anomalies at worst-mismatch bus" scans the incident branches
for zero impedance, off-nominal taps, large phase shifts or an outlying
reactance; "Recommendations" closes with next steps. It reuses the Q-limit
event, PV-limit and limit-validation printers, and a diagnostic exception
stays in the file without changing a successful result. Old run directories
may still hold `diagnose.txt`, which the viewer lists as well.
A run that attempts the current-iteration pre-solve may add
`current_iteration_start.log`: whether the candidate was attempted, accepted,
rejected by a guard, or rejected with the original start values restored.

### Comparing two runs

Tick **Compare** on two history rows (two runs of the same case under
different settings: a Q-limit mode, a solver, a start strategy) and press
**Compare selected runs**; the page shows side by side:

- case file, status, converged flag, iterations and final mismatch, with a
  link to each result page;
- the configuration keys the runs disagree on, read from their
  `effective_config.yaml` (the `_config_sources` bookkeeping is left out);
- the buses each run clamped, from `q_limit_events.csv`, and which only one
  of them clamped;
- total active and reactive losses from `branch_flows.csv`, with the
  difference;
- the largest voltage deviation and the ten buses that differ most, when both
  runs wrote `bus_voltages_complex.csv` (the detailed CSV export).

Everything comes from what the runs wrote, so earlier sessions compare as
well; exactly two runs are required. The box appears only on finished
power-flow runs (plain, diagnose, or started from an estimate) whose result
is still on disk, not on state estimation, short circuit, N-1 or import
analysis runs. Network size is not a criterion (the same case with an
external-grid source and with a slack is exactly the pair one wants); the
page says so when the runs share no bus.

## Run history and operation log

`start_sparlectra_webui` refreshes `powerflow_runs_index.json` under the
output root before serving, so runs of an earlier process appear at once;
**Refresh registry** reloads by hand, and missing, corrupt or unsafe entries
are skipped without hiding the valid ones. History is newest first: local
date and time, run ID, status text and badge (green successful, yellow
warning/partial, red failed, gray unknown, blue running, always with visible
text), solver summary fields and actions; indexes without timestamps use
the `result.json` modification time. **Case file** and **Config file** show
the file name, the full path in the tooltip. **Delete** removes one run,
**Delete all runs** every registered run under the configured root; both
update registry and index and delete only validated run directories, never
unrelated files.

The **Operation Log** page shows and downloads `logs/webui_operations.jsonl`,
a JSON Lines support log independent of the run directories, for error
reports. It records page opens, submissions and validation failures, run
lifecycle changes, artifact views and downloads, abort requests, history
refreshes, deletions, shutdown requests and enabled diagnostics or timing
modes, not static CSS or image requests. Entries carry route, method,
status, run, case, artifact, message and timing fields where available,
`sparlectra_version` and a UTC timestamp `yyyy-mm-ddTHH:MM:SS.sssZ`, never
file contents or configuration bodies. Logging is best effort and cannot
fail a request; phase events stay high level, per-iteration timings belong
to `performance.log`.

Every start drops entries older than `webui.operation_log_retention_days`
(default 10, `0` keeps only the current session) from every operation log
it knows, including one a service call wrote next to the runs;
`SPARLECTRA_WEBUI_OPERATION_LOG_RETENTION_DAYS` overrides the key for
headless runs. On append, a log above 1 MiB is compacted and more than
5,000 valid entries are cut to the newest 1,000, dropping malformed lines.
The page states its size (entries and kB) and offers **Clear operation
log**, which empties the file and writes one entry recording that. A large
log is a matter of entry count, not age; a retention of two or three days
is the lasting remedy.

## [Configuration](@id webui-configuration)

**Check configuration** and **Refresh configuration** bring a local file up
to a template that gained options. Check is a dry run against `src/config/configuration.yaml.example`: missing keys,
known deprecated aliases, duplicate YAML keys and a preview of the refreshed
YAML, nothing written. Refresh keeps existing values, adds missing keys with
template defaults and may normalize deprecated aliases (start-mode settings,
legacy `matpower_*` Q-limit mode names); it never runs at startup, writes a
timestamped backup first, refuses on duplicate keys, and offers uploaded or
pasted YAML as a download instead of rewriting it.

The **Configuration Editor** opens the active YAML in a textarea, validates
it with the same parser and duplicate-key checks, writes only after
validation with a timestamped backup, and warns when a case configuration
file still overrides the global YAML through the prefilled form. After a
save the server reloads the file, so the form shows the new values on the
next page load.

Form values override the YAML configuration; the order is described in
[Configuration](configuration.md). Each run writes `effective_config.yaml`.
**Ignore Web UI settings and use configuration defaults** runs with the YAML
values only. `model.auto_profile: apply` may still adjust MATPOWER import
conventions after form and API values are assembled;
`effective_config.yaml` records that.

**Help** in the header opens this page inside the Web UI, as shipped with
the running version; **Project Docs** opens the published documentation
(`webui.docs_base_url`). The help pages behind the **?** icons carry a
section of this page in full, or the lead paragraph of a library section
with the link to it.

## Shutdown

The server stops with **Stop Web UI**, `close(server)` or `Ctrl+C`. With
`auto_shutdown_on_browser_close = true` it also stops when the browser has been
gone for `browser_heartbeat_timeout_seconds` (best effort).

## Limits

The Web UI is local by design: loopback only, no authentication, no
public-server mode, no database. There is no push channel (WebSockets or
server-sent events), so a running job reports its phase through the page;
there is no topology view and no plotting. The HTTP layer is a compact Julia
`Sockets` implementation, chosen to avoid a web-framework dependency.
