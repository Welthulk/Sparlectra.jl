# Web UI

The Web UI is a local browser interface for Sparlectra: choose a case,
start a calculation, read the result. It runs on your own machine, is
reachable from that machine only and has no login. Operation and details
are on [Web UI Reference](webui_reference.md).

## Start

Install and start in one line. Julia is installed if it is missing:

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
environment once and offers to build a system image, which makes every
later start fast ([Sysimage](sysimage.md)). The questions of the
installer, the start from the Julia REPL and the directories the Web UI
uses are described in the
[reference](webui_reference.md#Installation-and-start).

## [Case formats at a glance](@id webui-case-formats)

| Format | How the case gets in | Download possible | Good to know |
|---|---|---|---|
| MATPOWER (`.m`) | **Import case files**, or type a name into **Or type case file path** | Yes: a standard name such as `case14.m` or `case118.m` is fetched once from the MATPOWER project | The download needs an internet connection the first time, afterwards the case is local |
| CGMES | **Import case files** with the delivery ZIP, or with the profile files (EQ, SSH, TP, SV) selected together | Yes: `cgmes:<alias>` fetches the ENTSO-E test configurations once (about 22 MB) | Every uploaded delivery is checked for completeness; the aliases are listed under [CGMES test configurations](@ref webui_cgmes_test_configurations) |
| PowSyBl IIDM (`.xiidm`) | **Import case files** | No | The file is read as it is; unpack a compressed file first. See [PowSyBl IIDM files](@ref webui_powsybl) |
| Sparlectra Case Format, power-grid-model (`.json`) | **Import case files**, or **Export as SCF case file** from any selected case | No | Scenario files and the scenario editor need this format |
| DTF (`.DAT`) | **Import case files** | No | Experimental, for diagnostics and validation |
| Shipped demo cases | Already in the selector | Not needed | They work offline |

Nothing is downloaded without a name you typed. A selected or uploaded
file is used as it is.

## Loading a case

The **Case** page holds the case selector, the import, the exports and
the import options of the selected format.

- **Existing case file** lists the cases of the case directory and the
  shipped demo cases ([Shipped Demo Cases](demo_cases.md)).
- **Import case files** takes several files at once and copies them into
  the case directory. It starts no calculation and never overwrites an
  existing file.
- **Or type case file path** takes a name to download (see the table) or
  the path of a local file.
- **Case input format** stays on **Auto** unless the extension does not
  say what the file is.
- **Export as SCF case file**, **Save case as** and **Download selected
  case** write or deliver the case. **Save case as** keeps variants of one
  network side by side.

Details: [Case selector](@ref webui-case-selector),
[Import case files](@ref webui-import-case-files),
[Save case as](webui_reference.md#Save-case-as).

## Running

The **Runs** page starts the calculations for the selected case. Every
start opens a status page that follows the run, and a running job can be
aborted there. The options of a run are on the **Settings** page.

### Power flow

**Start PowerFlow run** solves the network with the selected solver:
Newton-Raphson, APSLF or the linear DC model. The start values come from
the case unless **Flat start** is set. Reactive limits of the machines are
enforced in the mode chosen under **Q-limit handling**. See
[Power-Flow Configuration](powerflow_configuration.md).

### Diagnose

**Diagnose** evaluates the mismatch at the operating point the case
brings along, without a corrective step. The result shows where the
imported model does not balance, bus by bus. It is the first thing to
run when a case does not converge. See
[Diagnose](webui_reference.md#Diagnose).

### N-1

The outage kind (branches or generators) selects the generated outage
list; a case in the Sparlectra Case Format can carry scenarios of its
own. The result ranks the cases by severity, failures first. Screening
and weights are optional. See
[N-1 Contingency Analysis](contingency.md).

### State estimation

The state-estimation section selects a measurement set, or generates one
from a solved power flow. A run checks the observability, estimates the
state and reports suspicious measurements. The result page offers a
power flow from the estimate. See
[State Estimation](state_estimation.md).

### Short circuit

**Short circuit** computes the balanced three-phase currents of every
bus, maximum and minimum case, without a power-flow solve. The button is
offered when the case carries source data; a run is `succeeded` only when
every source carries its data. See
[Short-Circuit Analysis](short_circuit.md).

## Results

Every finished run has a result page with status, iterations, final
mismatch and the list of its artifacts. A CSV artifact opens as a table,
a Markdown report rendered, and **Download all artifacts as ZIP** packs
the run. **Run history** lists the runs, newest first; two runs of the
same case can be compared side by side. The **Operation Log** records
what the Web UI did and is the file to attach to an error report. See
[Results and artifacts](webui_reference.md#Results-and-artifacts) and
[Output modes](@ref webui-output-modes).

## Settings

The **Settings** page holds the options that do not belong to a single
run: solver, start, output. **Save settings** writes either for the
selected case or into the configuration file. With a case selected, every
page shows the values a run of that case will use: case settings before
the configuration file before the defaults. The **Configuration Editor**
edits the file as text and checks it before it is saved. See
[Form options](@ref webui-form-options) and
[Configuration](@ref webui-configuration).

**Help** in the header opens this page inside the Web UI. **Search the
help** looks through this page, the reference and the help of the
controls. The **?** next to a control opens the help of that control.

## Limits

The Web UI is local by design: one user, the own machine, no login and
no public-server mode. A running job reports its progress through the
status page, which refreshes by itself. There is no topology view and no
plotting.
