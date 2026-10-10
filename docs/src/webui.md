# Web UI

The Web UI is a local browser interface for Sparlectra: choose a case,
start a calculation, read the result. It runs on your own machine, is
reachable from that machine only and has no login. Every control, route
and output mode: [Web UI Reference](webui_reference.md).

## Start

Install and start in one line; Julia is installed if it is missing:

```sh
curl -fsSL https://raw.githubusercontent.com/Welthulk/Sparlectra.jl/main/tools/install_webui.sh | sh
```

```powershell
iwr -useb https://raw.githubusercontent.com/Welthulk/Sparlectra.jl/main/tools/install_webui.bat -OutFile install_webui.bat; .\install_webui.bat
```

From an installed copy or a checkout:

```sh
julia --project=. start_webui.jl
```

`start_webui.sh` and `start_webui.bat` do the same. The browser opens
`http://127.0.0.1:8080/powerflow`. The first start prepares the
environment once and offers to build a system image for fast later starts
([Sysimage](sysimage.md)). Installer questions, REPL start and
directories: [reference](webui_reference.md#Installation-and-start).

## [Case formats at a glance](@id webui-case-formats)

| Format | How the case gets in | Download | Good to know |
|---|---|---|---|
| MATPOWER (`.m`) | **Import case files**, or a name typed into **Case file** | Yes: a standard name such as `case14.m` is fetched once from the MATPOWER project | The first download needs an internet connection |
| CGMES | **Import case files**: the delivery ZIP, or the profile files (EQ, SSH, TP, SV) selected together | Yes: `cgmes:<alias>` fetches the ENTSO-E test configurations once (about 22 MB); aliases under [CGMES test configurations](@ref webui_cgmes_test_configurations) | Every uploaded delivery is checked for completeness |
| PowSyBl IIDM (`.xiidm`) | **Import case files** | No | Unpack a compressed file first; see [PowSyBl IIDM files](@ref webui_powsybl) |
| Sparlectra Case Format, power-grid-model (`.json`) | **Import case files**, or **Export as SCF case file** from any selected case | No | Scenario files and the scenario editor need this format |
| DTF (`.DAT`) | **Import case files** | No | Experimental, for diagnostics and validation |
| Shipped demo cases | Already in the selector | Not needed | Work offline |

Nothing is downloaded without a name you typed. A selected or uploaded
file is used as it is.

## Loading a case

The **Case** page holds the case selector, the import, the exports and
the import options of the selected format. **Case file** picks a case of
the case directory or a shipped demo case
([Shipped Demo Cases](demo_cases.md)); a typed name downloads, a typed
path loads a local file. **Import case files** copies files into the
case directory and never overwrites one. **Case input format** stays on
**Auto** unless the extension is ambiguous. **Save case as** keeps
variants of one network side by side. Details:
[Case selector](@ref webui-case-selector),
[Import case files](@ref webui-import-case-files),
[Save case as](webui_reference.md#Save-case-as).

## Running

The **Runs** page starts the calculations for the selected case; every
start opens a status page that follows the run and can abort it. The
options of a run are on the **Settings** page.

### Power flow

**Start PowerFlow run** solves the network with Newton-Raphson, APSLF or
the linear DC model. Start values come from the case unless **Flat
start** is set; reactive limits follow **Q-limit handling**
([Power-Flow Configuration](powerflow_configuration.md)).

### Diagnose

**Diagnose** evaluates the mismatch at the operating point the case
brings along, without a corrective step, bus by bus. Run it first when a
case does not converge ([Diagnose](webui_reference.md#Diagnose)).

### N-1

The outage kind (branches or generators) selects the generated outage
list; a Sparlectra Case Format case can carry its own scenarios. The
result ranks the cases by severity, failures first
([N-1 Contingency Analysis](contingency.md)).

### State estimation

The state-estimation section selects a measurement set or generates one
from a solved power flow, checks the observability, estimates the state
and reports suspicious measurements; the result page offers a power flow
from the estimate ([State Estimation](state_estimation.md)).

### Short circuit

**Short circuit** computes the balanced three-phase currents of every
bus, maximum and minimum case, without a power-flow solve. The button is
offered when the case carries source data; a run is `succeeded` only when
every source carries its data
([Short-Circuit Analysis](short_circuit.md)).

## Results

Every finished run has a result page with status, iterations, final
mismatch and its artifacts (CSV as a table, Markdown rendered,
**Download all artifacts as ZIP**). **Run history** lists the runs and
compares two runs of one case. The **Operation Log** records what the
Web UI did; attach it to an error report. See [Results and artifacts](webui_reference.md#Results-and-artifacts)
and [Output modes](@ref webui-output-modes).

## Settings

The **Settings** page holds the options that do not belong to a single
run: solver, start, output. **Save settings** writes for the selected
case or into the configuration file; with a case selected, every page
shows the values a run of that case will use (case settings before the
configuration file before the defaults). The **Configuration Editor**
edits the file as text and checks it before saving. See
[Form options](@ref webui-form-options) and
[Configuration](@ref webui-configuration).

**Help** in the header opens this page inside the Web UI, **Search the
help** searches this page, the reference and the help of the controls,
and the **?** next to a control opens its help.

## Limits

One user, the own machine, no login, no public-server mode, no topology
view, no plotting. Technical limits: [reference](webui_reference.md#Limits).
