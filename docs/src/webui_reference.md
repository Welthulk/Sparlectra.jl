# Web UI Reference

Operation and details of the [Web UI](webui.md): installation,
directories, every form option, the run kinds, artifacts, history and
configuration.

The Web UI is presentation, form parsing and route handling on top of the
[PowerFlow service](powerflow_service.md); every calculation runs through
`start_powerflow_run` and `run_sparlectra_api`. It binds to loopback
(`127.0.0.1`, `localhost`, `::1`) and has no authentication; a wider
binding would expose the run directory and the configuration editor.

## Installation and start

### One-line install

Installs Julia if missing (via juliaup), downloads the latest tagged release
into `Sparlectra/` in the current directory and starts the Web UI:

```sh
curl -fsSL https://raw.githubusercontent.com/Welthulk/Sparlectra.jl/main/tools/install_webui.sh | sh
```

```powershell
iwr -useb https://raw.githubusercontent.com/Welthulk/Sparlectra.jl/main/tools/install_webui.bat -OutFile install_webui.bat; .\install_webui.bat
```

The scripts (`tools/install_webui.sh`, `tools/install_webui.bat`) ask
three questions:

- **Update** when an existing copy is older than the latest release; the old
  copy is kept as `Sparlectra.old`.
- **Sysimage**: built only with `SPARLECTRA_BUILD_SYSIMAGE=1`, otherwise
  left to the Web UI start ([Sysimage](sysimage.md)).
- **Desktop shortcut**: Windows `Sparlectra Web UI.lnk`; Linux an
  application-menu entry plus a `.desktop` file on the desktop (GNOME needs
  a one-time right-click "Allow Launching"); macOS a desktop symlink.

Unattended installs answer with `SPARLECTRA_UPDATE=1/0`,
`SPARLECTRA_BUILD_SYSIMAGE=1/0` and `SPARLECTRA_CREATE_SHORTCUT=1/0`.

### Start script

From a checkout or an installed copy:

```sh
julia --project=. start_webui.jl
```

`start_webui.sh` / `start_webui.bat` run the same script and point at the
install script when Julia is missing; a Windows desktop shortcut is
right-click `start_webui.bat`, *Send to > Desktop (create shortcut)*. The
script resolves, instantiates and compiles a missing dependency or
manifest once, for the library and for `app/` (`--env-only` does only
that), then looks for a usable sysimage and offers to build one
([Sysimage](sysimage.md)). The start scripts set `JULIA_NUM_THREADS=auto`
unless it is set; a case with several islands then solves its large
islands in parallel ([Parallel Execution](parallel_execution.md)).

### From the Julia REPL

```julia
using Pkg
Pkg.activate("path/to/Sparlectra/app")
Pkg.instantiate()   # first time only
using SparlectraApp

server = start_sparlectra_webui(open_browser = true)
wait(server.task)
```

A REPL session runs without the sysimage ([Sysimage](sysimage.md)). `start_sparlectra_webui` returns a
`SparlectraWebUIServer` handle and serves as soon as the socket is bound.
`open_browser = true` opens an app-style window (Microsoft Edge, Google
Chrome, Chromium or Brave); without one of them the URL
`http://127.0.0.1:8080/powerflow` is logged. If the port is occupied, stop
the old Julia process or pass another `port`.

### Directories and files

Results go to `%LOCALAPPDATA%\Sparlectra\WebUI\runs` (Windows),
`$XDG_STATE_HOME/sparlectra/webui/runs`, default
`~/.local/state/sparlectra/webui/runs` (Linux) or
`~/Library/Application Support/Sparlectra/WebUI/runs` (macOS). The
operation log is in the sibling `logs`, downloaded and generated cases in
the sibling `data/mpower`, shared with `ensure_casefile` and the test
suite (`SPARLECTRA_LARGE_CASES_DIR` overrides it for all three). The
bundled `warmup_case*.jl` workloads land there too, hidden from the
selector.

The first start copies the package configuration template to the
user-writable `config/configuration.yaml` (`output_root = "..."` and
`config_file = "..."` override the defaults; an explicit file is never
overwritten). Next to it: `configuration.template.yaml` (a key still at
an old template default follows a new default at start, reported, the
previous file kept as `configuration.yaml.template-follow.bak`),
`configuration.yaml.user-keys.txt` (keys the Web UI itself saved, never
followed) and `configuration.yaml.settings-save.bak` (backup of a
settings save). The browser cannot change the output root. The **Info**
panel names the serving build (`sysimage, built ...` or `native session`)
and links to the Sysimage page.

## Importing cases

The **Case** page carries the case chooser, upload, export and the
per-format import options.

### [Case selector](@id webui-case-selector)

**Case file** is one combobox. Its list holds the MATPOWER `.m` files and
runnable DTF `.DAT` candidates of the case directory, the shipped demo
cases ([Shipped Demo Cases](demo_cases.md)), the three CGMES deliveries
under `data/cgmes_demo` as `<case>_cgmes.zip`, and the PowSyBl IIDM files
(`.xiidm`). Generated `.jl` cache files and `warmup_` files
are hidden. FOR002-like `.DAT` files go into the optional **FOR002
reference** field, for legacy reference comparison only.

A typed value overrides the list: a bare name such as `case14.m` or
`case9241pegase.m` is downloaded into the case directory, an existing
`.m`/`.jl` path is used as given; path-like missing inputs and URLs are
rejected. A generated `.jl` cache file resolves back to its `.m` source
and is rejected without one.

**Case input format** defaults to **Auto** (MATPOWER, CGMES folders and
ZIPs, PowSyBl IIDM, and native DTF where the reader takes the file as a
network). **CGMES (ENTSO-E, folder or ZIP)** is preselected for a `.zip`
or a directory, **PowSyBl (IIDM file)** for `.xiidm` (it shows the
PowSyBl import options, [PowSyBl IIDM files](@ref webui_powsybl));
**Sparlectra Case Format (.scf.json)** and **power-grid-model JSON
(input.json)** name the same reader, which Auto already uses for every
`.json`. The native DTF path is experimental; **Selected DTF outage
labels/indices** takes one label or index at a time.

The format-bound options of the Basic set render directly, the rest and
the DTF diagnostics fold into **Advanced**. Deleting a case removes its
companion files (configuration, weights, measurement CSVs). A collapsible
MATPOWER acknowledgement names the citations the cases ask for.

### [Import case files](@id webui-import-case-files)

**Import case files** opens the native file picker and accepts several
files at once: MATPOWER `.m`/`.M`, DTF `.dat`/`.DAT`, CGMES `.zip`
deliveries, CGMES profile files (`.xml`), PowSyBl IIDM files (`.xiidm`)
and Sparlectra Case Format `.json` files (validated before they are
stored). Several `.xml` files (the EQ, SSH, TP and SV profiles of one
delivery) are packed into one `<stem>_cgmes.zip`; the set must contain the
EQ profile. A CGMES ZIP may hold nested ZIPs.

Importing is copy-only: no run, no parsing of `.m` code; files land in
the case directory the form shows. Limits: 100 MiB per file, 250 MiB per
request. An existing file is never overwritten (`already exists`, the
other files still import); names outside the case directory are refused.
Enter in **Case file** resolves a typed value the same way: a bare case
name via [`ensure_casefile`](@ref), `cgmes:<alias>` as below, or a full
local path. Re-importing a saved set brings the `.config.yaml` sidecar (a
case-scope file, `scope: case`) along with the `.scf.json` and its
measurement sets when all are selected together.

### Export and download

The export row writes the open case as `.scf.json` or plain PGM
`.pgm.json`, replacing a format suffix the case already carries
(`case14.scf.json` as PGM gives `case14.pgm.json`); both appear in the
selector as runnable cases. **Download selected case** serves what the
selector shows, after an export also the file just written: only plain
files inside the case directory, compared by real path (the state
directory can sit behind a symlink); anything outside, and a CGMES
delivery directory, is refused.

### Save case as

**Save case as**, next to the export buttons, saves the current case under
a new name: one network with several variants, one file per variant. It
writes into the case directory:

- `<name>.scf.json` via [`exportSCF`](@ref), with `meta.case_name = <name>`
  and `meta.source_reference` naming the source; a MATPOWER or DTF source
  is saved as SCF, the original file is never touched;
- `<name>.config.yaml`: the source case's saved settings merged with
  unsaved changes on the page;
- `<name>.measurements.csv` for each measurement set bound to the source
  case, with its `# case:` header rewritten; a second bound set becomes
  `<name>_<original stem>.measurements.csv`.

An existing `<name>.scf.json` is refused unless "overwrite" is ticked.
"Start from the solved state" solves the case once with the effective
settings and writes the solved voltages as the copy's
`sparlectra.start_state`. Installation-scope settings (`output.*`,
`benchmark.*`, `runtime.*`) are never written into the case configuration
file. File layout: [Sparlectra Case Format](scf.md).

### [PowSyBl IIDM files](@id webui_powsybl)

IIDM is the network format of the PowSyBl framework; its XML form has the
extension `.xiidm` and is read as it is.

| You have | Do this |
|---|---|
| no file yet | pick one of the four shipped networks: `ieee14.xiidm`, `four_substations.xiidm`, `micro_grid_be.xiidm`, `ieee14_sc.xiidm` |
| a `.xiidm` file | **Import case files**; it appears in the selector |
| a compressed file (`.xiidm.bz2`, `.gz`) | unpack it first, then import the `.xiidm` |

Every run kind and export works on a PowSyBl case; scenarios need an SCF
export of it. With a PowSyBl case selected, **Input format** shows five
options; the defaults fit most files.

| Option | Default | Change it when |
|---|---|---|
| System base MVA | 100 | you compare per-unit values with a tool on another base |
| HVDC converters | fixed injections | never for now: the paired controller is not available for PowSyBl files |
| Remote voltage regulation | hold the target at the unit's own bus | a generator regulates another bus that shall be held exactly: choose the outer-loop control |
| One slack per synchronous component | on | only the first component shall be solved |
| Slack generator ids | empty | another unit than the chosen one shall carry the reference; several ids are separated by semicolons |

Every power-flow run on a PowSyBl case writes `powsybl_import.log`:
elements read and built per type, skipped elements with their reason, the
reference generator of each synchronous component, and whether the solve
started from the file's voltages or flat. Short circuit uses the connected
generators with the reactance of the extension `generatorShortCircuit`; a
generator without it enters with a default reactance and the rows it feeds
are marked; without any extension or rated power the **Short circuit**
button is not offered (of the shipped networks only `ieee14_sc.xiidm`
carries complete data, `micro_grid_be.xiidm` runs on defaults). Not read: compressed files, time
series, the regulation data of tap changers. There is no IIDM export.
Conventions: [PowSyBl Import](powsybl_import.md).

### [CGMES test configurations](@id webui_cgmes_test_configurations)

`cgmes:<alias>` typed into **Case file** downloads the official ENTSO-E
test-configuration package once (about 22 MB) into the CGMES cache
(`data/CGMES`, overridable with `SPARLECTRA_CGMES_CACHE`) and packs the
requested configuration with its boundary set into `cgmes_<alias>.zip` in
the case directory. Aliases: `microgrid_be`, `microgrid_nl`,
`microgrid_assembled`, `smallgrid`, `smallgrid_nb`, `fullgrid`,
`fullgrid_nb`, `realgrid`. Repeated requests reuse the ZIP. If the download
fails, the message names the path where the package can be placed by hand.

### Import analysis on upload

Every uploaded CGMES delivery is checked at once. A complete one reports
"ready to compute"; an incomplete one gets the import analysis: the
message names the missing declared `md:Model.DependentOn` dependencies by
model id (typically the boundary set), and the report is written next to
the case as `<case>.import_analysis.txt`. The same analysis is appended to
`cgmes.log` whenever a run's CGMES import aborts. Scripted: the run mode
`import_analysis_mode` of `start_powerflow_run` (a non-importable delivery
ends as a failed run with reason `import_analysis_not_importable`) or
`analyzeCGMES`. **Require boundary set** (CGMES section, saved per case)
submits `cgmes_import.require_boundary`; unchecked, an incomplete delivery
imports where possible, but buses without a resolvable `BaseVoltage` still
abort.

### Long actions on a case

Generating a measurement set or adding noise runs inside the request and
can take a while on a large case (with **critical measurements** the
generator re-reads the criticality after every removed row). Meanwhile the
case is busy: a second generator action or a run is refused with a
message, and a queued or running run blocks the generator the same way.

## Starting a run

Run-independent options live on the **Settings** page and reach a run
through the configuration precedence (`resolve_config`): saving with
target "this case" writes case-scope keys into the case's configuration
file, target "configuration file" merges into the general YAML (backup
`.settings-save.bak`). The **Runs** page keeps the run parameters (modes,
N-1 kind, scenario source, screening, DTF outage fields), the
state-estimation section (submitted with `se_mode`) and the N-1 editors
as tabs. The measurement generator's
fields are the generator options of
[State Estimation Measurements](state_estimation_measurements.md#se-generator-v2).

### [Form options](@id webui-form-options)

The settings and run forms offer:

- case name or local path, configuration template file, the read-only
  output root (set only through `start_sparlectra_webui(; output_root)`),
  tolerance (decimals or scientific notation such as `1e-8`; the spinner
  steps the exponent and keeps the mantissa), maximum iterations,
  autodamping with its minimum factor;
- **Q-limit handling**: the `power_flow.qlimits.enabled` checkbox and the
  enforcement-mode selector, which also offers **off**; **PV/PQ
  switching** holds the hysteresis, the final check bound
  (`final_q_accept_pu`, `auto` or a value in pu), the classic pass limit
  (`classic_max_passes`), the first switching iteration (`start_iter`),
  the switch cap per bus (`guard.max_switches`), freezing after repeated
  switching and the guard switch (`guard.enabled`); **Advanced** holds the
  rest of `power_flow.qlimits`
  ([Q-limit options and guard](powerflow_configuration.md#pf-qlimits)).
  Per-unit values take scientific notation, bus lists comma-separated case
  bus numbers. Controls a mode does not read are grayed: everything while
  the handling is off, the start iteration and the guard rules in the
  classic modes, the range rules without the guard;
- **Solver**: one radio group `power_flow_solver` with `rectangular` (AC
  Newton-Raphson), `apslf` (AnalyticLoadFlow) and `dc` (linear screening
  model)
  ([solver selection](powerflow_configuration.md#pf-solver-selection),
  [DC power flow](powerflow_configuration.md));
- wrong-branch detection mode, start angle and voltage modes;
- **Advanced start values**: the current-iteration pre-solve fields
  (`power_flow.start_current_iteration.*`,
  [guarded pre-solve](powerflow_configuration.md#pf-current-iteration-start));
- **Merit-function line search**: `power_flow.merit.enabled`,
  `power_flow.merit.armijo_c1`, `power_flow.merit.fallback_max_mismatch`
  ([options](powerflow_configuration.md#pf-merit)),
  requires `power_flow.autodamp = true`;
- **Non-convergence handling** (Advanced): `power_flow.rescue`,
  `power_flow.auto_slack` and `power_flow.dc.fallback`
  ([Power-Flow Configuration](powerflow_configuration.md));
- MATPOWER import conventions: auto-profile mode (`off`, `recommend`,
  `apply`) plus transformer ratio, phase shift, bus shunt, PV voltage source
  and comparison reference;
- logfile output mode, performance timing and post-run diagnostics.

Only keys in `GUI_EDITABLE_CONFIG_KEYS` are submitted, and each run writes
`effective_config.yaml`. A **?** next to a control opens its help page:
the hint, the matching section of this documentation as shipped with the
running version, and a link to the published site (`webui.docs_base_url`
points it at a local build or a pinned version). **Back** keeps the
entered values. A field that does not apply to the selected solver or
step control is grayed in place, `disabled` and not submitted.

### AC, APSLF and DC

`AC (Newton-Raphson, rectangular)` reveals the start-value block (**Flat
start**, **Use APSLF start values** with its order field, **Use DC start
values**) and every AC-only option. By default
(`power_flow.start_mode.angle_mode = dc`) Newton-Raphson seeds its start
angles from an internal DC pre-solve (**Start angle mode**), unrelated to
the standalone `DC` solver.

`APSLF (AnalyticLoadFlow)` reveals the **APSLF solver options** (highest
coefficient (order), convergence radius, NR polish). **Use APSLF start
values** in the Newton-Raphson block seeds the rectangular solver instead
and is mutually exclusive with the APSLF solver (the configuration rejects
both). An unchecked checkbox (APSLF, Q-limit, current-iteration) submits
`false` rather than omitting the key.

`DC` replaces Newton-Raphson with the linear model
([DC power flow](powerflow_configuration.md))
and grays every AC-only option. The result page marks a DC run with a **DC
solution** badge next to **Solver** and the history column shows `dc`, so
the implicit `Vm = 1.0 pu` and lossless flows are never taken for an AC
result; `iterations` is `1`.

The **Solver** column names the method that ran: the power-flow solver for
a plain power flow, one started from an estimate and every contingency
case; `wls` for a state estimation; `iec60909` for a short circuit; empty
for an import analysis.

#### [Flat start](@id webui-flat-start)

`power_flow.flatstart` starts the Newton-Raphson solve at 1.0 pu and 0
degrees on every bus and ignores the imported start voltages (MATPOWER
`VM`/`VA`, CGMES `SvVoltage`, SCF `start_state`). Off, the imported values
seed the solve; a delivery is built around its own operating point, so
that is the default.

The checkbox sets only the start profile; the APSLF and DC start values,
the current-iteration pre-solve and the start projection follow their own
settings and run from the flat profile. While the flat start is on, the
two start modes count as `classic`, the page greys them (saved values
stay) and the run log names the override. A bare flat start needs the
start machines off (`power_flow.start_mode.start_projection` is a
configuration-file key).

**Ratio profile as a start candidate** (`power_flow.start_mode.ratio_profile`,
on by default) acts only on a flat start: the start projection also tries
the flat profile with each voltage level scaled by the transformer ratios
on its path from the slack, and keeps it when its mismatch is at least 10
percent lower. It lets a flat start converge on networks with stiff
off-nominal windings ([Start strategies](start_strategies.md)).

A CGMES run honours the flat start under **CGMES start values** = `auto`;
an explicit `sv` on the Case page still starts from the delivery state
([Power flow configuration](powerflow_configuration.md#pf-solver-core)).

### Case-specific settings

A successful result page offers **Save settings for this case**: the form
options of that run plus traceability metadata go into a sanitized YAML
profile below the output root. A non-converged run shows **Save these
settings anyway**; nothing is saved automatically. Reopening the case
prefills the form from the profile with a notice; a manual edit wins for
the submitted run.

With a case selected, every page shows the values a run of that case will
use: case settings over the configuration file over the defaults; without
a case, the configuration file. **Show the configuration file's values**
(`?config_view=1`) switches the Settings page to the file alone, **Show
the case settings** back. A save to the configuration file says how many
of the saved keys the selected case overrides. Between configuration file
and saved profile the last edit wins: a newer YAML file takes precedence
for the keys it sets on the next page load, and the notice says so.

With a Sparlectra Case Format file selected, the form is seeded from the
case file (the page says which settings came from it); a saved profile
for the case wins over the file, and an edited control wins over both.

### Scenarios and N-1

The action row of the run form carries the N-1 controls: outage kind
(branch/generator), **scenario source** (generated N-1 lists, or the
scenarios block of an SCF case), **screening mode** (configured / off /
flag / only) with a **margin** field (empty means the configured
`contingency.screening.margin_pct`), and links to the weights editor and,
for SCF cases, the scenario editor. The run summary names the screened
count and the case-list source; the result page tabulates the cases
ranked by severity with failures first, capped at 100 rows (the linked
CSV keeps the complete list in input order). A screened row shows its
one-step estimate, a flagged row the estimate next to the full-run
values, a failed case its error text. Screening, weights and the warm
start themselves: [N-1 Contingency Analysis](contingency.md).

**Weights.** "edit N-1 weights" next to the outage-kind selector opens a
per-case editor. Weights live next to the case as
`<case-stem>.contingency-weights.csv` (the two-column format of
`readContingencyWeightsCSV`), hidden from the case list and deleted with
the case. The editor seeds a table with the case's element names and
offers a raw-CSV text area and a file upload; a malformed CSV is rejected
with the line number, rows left at `1.0` are omitted on save. A run picks
the weights up whenever the file exists; names that match no element are
reported in `run.log`. The file applies to case-list runs (the outage-kind
selector and the `n1_*` scenario sources); a scenario run from a case
file's block or an external scenario JSON uses the per-scenario `weight`.

Scenarios need an SCF case; for other formats the form says so next to the
export action. The **scenario editor** lists the open case's scenarios
(edit, duplicate, delete, new) and edits one as op rows (status, set or
scale on a component chosen by class).
Validation is server side; errors name the scenario and op index, and a
tap patch on a regulated transformer is rejected with the controller
named. **Save into case file** writes the scenarios block through the SCF
writer.

### Diagnose

**Diagnose** next to **Start PowerFlow run** evaluates the mismatch at the
case's stored operating point (MATPOWER `VM`/`VA`, the `SvVoltage` state
of a CGMES delivery) without a corrective Newton step. Every start-value
machine is forced off (`flatstart = false`, `start_projection = false`,
`dc_seed_unconditional = false`, `start_current_iteration.enabled = false`,
`apslf_start.enabled = false`, plus `max_iter = 1` and
`qlimits.enabled = false`), so the residual reflects the imported model;
the forced settings outrank the form and the case configuration file. The
run uses the normal pipeline and adds `self_check.log` (forced settings,
start residual, for CGMES the count of buses without a usable
`SvVoltage`, started at `1.0 pu / 0°`), `self_check_residuals.csv`
(per-bus P/Q residuals) and `diagnose_self_check_config.yaml`. Programmatically:
[`run_fixed_reference_self_check`](@ref).

One step practically never converges; the residual is the measurement, so
history and result page label the run **diagnosed** on a neutral badge
(raw status in the tooltip). A diagnose run that could not run at all
keeps the failure vocabulary.

### Short circuit

**Short circuit** next to Diagnose evaluates the balanced three-phase
currents (IEC 60909-0, [`runShortCircuit!`](@ref)) for every bus, maximum
and minimum case in one run, without a power-flow solve. The button is
enabled only when the case carries short-circuit source data; otherwise
it is disabled with a tooltip. The run writes `short_circuit_max.csv` and
`short_circuit_min.csv` (per-bus `Ik''`, `Sk''`, `κ`, `i_p`, safety flag
and reasons) and a `run.log` with the data coverage report; rows with
defaulted or skipped data get a warning badge because a flagged `Ik''max`
is a lower bound. A MATPOWER or DTF case carries no source data and fails
with `short_circuit_requires_cgmes`. Which data counts as a source in each
format: [Short-Circuit Analysis](@ref short_circuit_source_data).

### State estimation

The state-estimation section of the Runs page
(`/powerflow#state-estimation`; `/stateestimation` redirects there) uses
the shared case selection, selects a measurement set and runs
observability (traffic light, structural-island note), the WLS solve and
the diagnostics of [State Estimation](state_estimation.md). Uploaded
`.csv` files are offered when they carry the v1 version comment. Only a
set bound to the selected case is preselected (such cases are starred); a
case file with its own measurements offers those first, and a set bound
to the case an SCF file was exported from stays usable.

Every generation and every "add noise" writes a new time-stamped file
(`<case>.measurements.<yyyymmdd-HHMMSS>.csv`, `<case>.noisy.<stamp>.csv`)
and arms it for the next run; later the newest set bound to the case is
armed. **Delete measurement file** removes the armed set (never a case
file). The set tab offers download, re-upload and an inline editor (the
first line must stay the version comment) and renders the
`sparlectra-taps v1` tap table of a generated file. **Add noise to this
set** perturbs each value with the sigma its own row declares. **Reset
saved settings** under the generator deletes the per-case settings
sidecar. **estimate taps** (on by default) releases every in-service
transformer with a ratio tap changer except machine transformers; the
result page shows the per-transformer table and the fixation J drop. The
form warns when `k_suppress < k_eliminate`, because suppressed rows then
rarely reach the elimination.

Artifacts (`measurements.csv`, `se_diagnostics.md`, `se_view.md`,
`se_state.csv`, `shunt_estimates.csv`, `se_tap_estimates.csv`,
`se_bad_data.csv`, `se_deltas.csv`) land in the run history with kind `se`.
The result page offers **Run power flow from this estimate** and the
topology hypothesis test behind a button (`topology_hypotheses.md`). The
summary shows `J`, `E[J] = dof` and the Wilson-Hilferty band verdict for
the state after the elimination; suppressed rows stay in `J` and surface
as `J_active`.

## Running jobs

A submission starts a background worker and redirects to the status page
(elapsed time, case paths, start time, a manual refresh link). Queued, running and aborting pages
refresh every two seconds (`autorefresh=1` requests are not logged as
user actions); terminal pages stop. The start page and the run history
show a banner with **Open status** and **Abort** while a job is queued or
running; such a run cannot be deleted.

An N-1 or scenario run shows how far its batch is: **Outages done** reads
`N-1: 37 / 186 outages` (a scenario set reads **Scenarios done**,
`Scenarios: 3 / 10 scenarios`). The counter appears once the first outage
has finished, after the import and the base-case solve; screened outages
count as done. For large cases, the result page, `run.log`, `result.json`
and `performance.log` show time per phase (reader, network builder,
solver, artifacts).

### Abort

A queued or running job can be aborted from its status page or from the
banner. Abort is cooperative: an operation that cannot be interrupted,
such as a sparse factorization, finishes first; the status page shows the
current phase. An aborted run keeps its run directory and is never
reported as success. If a run stays in `aborting` for more than 60 s, the
status page offers **Hard reset Web UI**: it marks the run invalid and
stops the server. Restart with `julia --project=. start_webui.jl`.

### Feedback

One-off messages (validation errors, import summary) open as a dismissible
popup; recent errors stay listed on the **Last errors** page.

## Results and artifacts

Every finished run gets a result page: run ID, schema version, status,
convergence and solution flags, iterations, final mismatch, reason and
message, input paths, output directory and a **Solver** entry. A CSV artifact opens as a table
(delimiter taken from the header line), a Markdown report rendered;
**Raw text** shows the file as written. Other text artifacts (JSON, YAML,
logs, HTML) open as escaped text, other files download, and every
artifact page offers a download. **Download all artifacts as ZIP** packs
the exposed artifacts into `sparlectra_run_<run_id>_artifacts.zip`.

### [Output modes](@id webui-output-modes)

- **Logfile output mode**: `classic` writes the result report plus a
  compact timing and status summary (`solver_time`, `representative_time`,
  iterations, final mismatch, outcome); `full` adds **Full run details**
  with the effective typed configuration, artifact choices and status
  diagnostics. MATPOWER `.m` runs also print the import-option and
  auto-profile tables; DTF `.DAT` runs never do.
- **Performance timing** (`off`, `compact`, `full`) writes
  `performance.log` with the phases of one request (from request parsing to artifact
  writing, solver and total time; `full` adds internal profile entries). A
  run from the Web UI solves once; repeated timing is the benchmark mode
  of `run_matpower_case` in a script
  ([Performance Profiling](performance_profiling.md#perf-benchmark)).
- **Export case as CGMES delivery (EQ+TP+SSH+SV, ZIP)** writes the case as
  one re-importable delivery into the run's artifacts, for every case format
  and also on non-converged runs; see [CGMES Export](cgmes_export.md).
- **Write bus/branch CSV files** (API default off, because large networks
  produce large files) writes `bus_voltages_complex.csv` and
  `branch_flows.csv` from the `ACPFlowReport` rows; columns on
  [Local PowerFlow Service](powerflow_service.md).
- **CSV format** is the machine-scope key `output.csv_format` (**Save
  settings** writes it to the configuration file, a per-case save leaves it
  out) and applies to every CSV of every run type; `auto` (the default)
  follows the regional settings of the machine. Values and delimiters:
  [Configuration](configuration.md). Excel may still auto-convert
  identifiers like `1E5`; use **Data > From Text/CSV** with text types when
  that matters.

### Diagnostic artifacts

`diagnose.log` is written by **Diagnose** or with `run_diagnostics = true`,
never by a normal run. On a non-converged run it names the worst-mismatch
bus and equation, classifies the mismatch trend (monotonic, oscillatory,
stagnant, diverging to non-finite) and the autodamp health, scans the
branches at the worst bus for zero impedance, off-nominal taps, large
phase shifts or an outlying reactance, and closes with recommendations. A
diagnostic exception stays in the file without changing a successful
result. A run that attempts the current-iteration pre-solve may add
`current_iteration_start.log`: whether the candidate was attempted,
accepted, or rejected (by a guard, or with the original start restored).

### Comparing two runs

Tick **Compare** on two history rows (two runs of the same case under
different settings) and press **Compare selected runs**; the page shows
side by side: case file, status, converged flag, iterations and final
mismatch with a link to each result page; the configuration keys the runs
disagree on (from `effective_config.yaml`); the buses each run clamped
(from `q_limit_events.csv`) and which only one of them clamped; total
losses (from `branch_flows.csv`) with the difference; the largest voltage
deviation and the ten buses that differ most, when both runs wrote
`bus_voltages_complex.csv`. Everything comes from what the runs wrote, so
earlier sessions compare as well; exactly two runs are required. The box
appears only on finished power-flow runs (plain, diagnose, or started
from an estimate) whose result is still on disk; the page says so when
the runs share no bus.

## Run history and operation log

`start_sparlectra_webui` refreshes `powerflow_runs_index.json` under the
output root before serving, so runs of an earlier process appear at once;
**Refresh registry** reloads by hand, and missing, corrupt or unsafe entries
are skipped without hiding the valid ones. History is newest first: local
date and time, run ID, status text and badge, solver summary fields and
actions. **Case file** and **Config file** show the file name,
the full path in the tooltip. **Delete** removes one run, **Delete all
runs** every registered run under the configured root; both delete only
validated run directories.

The **Operation Log** page shows and downloads `logs/webui_operations.jsonl`,
a JSON Lines support log independent of the run directories: page opens,
submissions, validation failures, run lifecycle changes, artifact views
and downloads, aborts, history refreshes, deletions and shutdown requests
(not static CSS or image requests), each with route, method, status, run,
case, artifact, message and timing fields where available,
`sparlectra_version` and a UTC timestamp, never file contents or
configuration bodies. Logging cannot fail a request; per-iteration
timings belong to `performance.log`.

Every start drops entries older than `webui.operation_log_retention_days`
(default 10, `0` keeps only the current session);
`SPARLECTRA_WEBUI_OPERATION_LOG_RETENTION_DAYS` overrides the key for
headless runs. On append, a log above 1 MiB is compacted and more than
5,000 valid entries are cut to the newest 1,000. The page states its size
and offers **Clear operation log**, which empties the file and writes one
entry recording that. A large log is a matter of entry count, not age; a
retention of two or three days is the lasting remedy.

## [Configuration](@id webui-configuration)

**Check configuration** and **Refresh configuration** bring a local file up
to a template that gained options. Check is a dry run against
`src/config/configuration.yaml.example` (missing keys, deprecated
aliases, duplicate YAML keys, a preview of the refreshed YAML), nothing
written. Refresh keeps existing values, adds missing keys with template
defaults and may normalize deprecated aliases; it never runs at startup,
writes a timestamped backup first, refuses on duplicate keys, and offers
uploaded or pasted YAML as a download instead of rewriting it.

The **Configuration Editor** opens the active YAML in a textarea, validates
it with the same parser and duplicate-key checks, writes only after
validation with a timestamped backup, and warns when a case configuration
file still overrides the global YAML. After a save the server reloads the
file. The search box above the text finds plain text (**Next** / **Previous**,
Enter / Shift+Enter) or jumps to a dotted key such as
`power_flow.qlimits.guard.enabled`.

Form values override the YAML configuration; the order is described in
[Configuration](configuration.md). Each run writes `effective_config.yaml`.
**Ignore Web UI settings and use configuration defaults** runs with the YAML
values only. `model.auto_profile: apply` may still adjust MATPOWER import
conventions after form and API values are assembled;
`effective_config.yaml` records that.

**Help** in the header opens the help ([Web UI](webui.md)) inside the Web
UI, as shipped with the running version; **Project Docs** opens the
published documentation (`webui.docs_base_url`). The help pages behind
the **?** icons carry a section of the help or of this reference in full,
or the lead paragraph of a library section with the link to it.

### [Search the help](@id webui_help_search)

The search field stands at the top of the help, of this reference and in
the header of every page and looks through the help, the reference and
the help pages of the controls. Every word entered must occur (a part of
a word is enough, case does not matter); hits appear while you type, hits
in a heading first, each with a marked line of text. `/` focuses the
field, Enter opens the first hit, Esc closes the list. Without a hit the
page links to the published documentation. The index ships with the
application; with scripts switched off the field is a plain form.

## Shutdown

The server stops with **Stop Web UI**, `close(server)` or `Ctrl+C`. With
`auto_shutdown_on_browser_close = true` it also stops when the browser has
been gone for `browser_heartbeat_timeout_seconds` (best effort).

## Limits

Loopback only, no authentication, no public-server mode, no database.
No push channel: a running job reports its phase through the page
refresh. No topology view, no plotting.
