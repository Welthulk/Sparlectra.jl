# Web UI

The Web UI is a local browser interface for Sparlectra: choose a case,
start a calculation, read the result. It runs on your own machine and is
reachable from that machine only, unless it is started in remote mode
([Remote use](@ref webui_remote)). Every control, route and output mode:
[Web UI Reference](webui_reference.md).

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

## [Remote use](@id webui_remote)

The Web UI can run on another machine, for example a GPU server, and be
used in the browser of your own machine. Start it on that host with the
remote switch:

```sh
./start_webui.sh --remote
julia --project=. start_webui.jl --remote
```

On a Windows host: `start_webui.bat --remote`. A remote start opens no
browser on the host and runs until it is stopped: **Stop Web UI** (asks
before it stops), `Ctrl+C` in its console, or the end of the process. To
keep it running after you log out of the host, start it in `tmux` or with

```sh
(umask 077; nohup ./start_webui.sh --remote > webui.log 2>&1 &)
```

When the server asks for a token, the console prints it, so it ends up in
`webui.log`; `umask 077` creates the log readable for you only. It applies
to a new file only: remove an old `webui.log` first, otherwise it keeps its
old permissions.

At start the console prints the lines to connect. The **Info** menu shows
them under **Server**: host name, bind address, port, mode and the address
to open, in the SSH mode also the tunnel command. Where the server asks for
a token, the address there shows `?token=<token>` as a placeholder; the
token itself is only on the host console. The page footer names the mode
as well.

### SSH tunnel (default)

The server stays on `127.0.0.1` of the host, and nothing is open to the
network. Without a token the server does not check who connects: every
user with a login on the host and every process on it can open
`127.0.0.1:8080` too. On a host you share with others, set

```yaml
webui:
  remote_ssh_require_token: true
```

and the SSH mode asks for the access token like the direct mode (see
there), or run the server on a host only you use. On your own machine open
a tunnel and keep it open while you work.

Linux, macOS and Windows Git Bash:

```sh
ssh -N -L 8080:127.0.0.1:8080 <user>@<host>
```

Windows PowerShell (the OpenSSH client of Windows 10 and 11):

```powershell
ssh.exe -N -L 8080:127.0.0.1:8080 <user>@<host>
```

Then open `http://127.0.0.1:8080/powerflow`. `-N` opens no shell, and
`Ctrl+C` ends the tunnel. If port 8080 is taken on your machine, use
another local port: `-L 9000:127.0.0.1:8080` and
`http://127.0.0.1:9000/powerflow`. If PowerShell does not find `ssh.exe`,
add *OpenSSH Client* under Settings, System, Optional features.

### Direct access

In a network you trust the server can bind a network address instead. Set
in the configuration file of the Web UI (the **Info** menu names it under
**Config file**; keys in [Configuration](configuration.md)):

```yaml
webui:
  remote_mode: direct
  remote_bind_address: 0.0.0.0
  remote_port: 8080
  remote_access_token: ""
```

and start with `--remote`. The server then accepts a request only with its
access token. With an empty `webui.remote_access_token` a new token is
generated at every start; a fixed one needs at least 16 characters and is
written in quotes (`remote_access_token: "0123456789abcdef"`), because the
configuration reader takes unquoted digits for a number. The console prints
the address with the token:

```text
http://<host>:8080/powerflow?token=<token>
```

Open it once: the browser keeps the token in a cookie, and the address bar
drops it. Browsers do not separate cookies by port, so the browser also
sends this cookie, and with it the token, to every other web service it
opens on the same host name. Open nothing else on that host name in this
browser that you do not trust with the token. A request without the token
gets `401`, and the operation log records it as `webui_access_denied`. The
traffic is plain HTTP, not encrypted: token, cases and results cross the
network readable, which is why the SSH tunnel is the default. The firewall
of the host must let the port through (firewalld:
`sudo firewall-cmd --add-port=8080/tcp`).

From the REPL: `start_sparlectra_webui(remote = true)` starts the remote
mode, and `webui_connection_info(server)` returns host, bind address, port,
mode, address and tunnel command.

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

### QSTS

**QSTS** runs an SCF time series step by step, `sequential` (taps and
banks carried from step to step) or `snapshots`; the button is offered
when the case is a series. The controllers come from
`<series>.control.yaml` next to the series. The result page shows the
switching, the steps with violations and a time series chart
([QSTS](@ref webui_qsts)).

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

No user accounts and no public-server mode: the remote mode is meant for
a host in your own network (SSH tunnel, or a direct bind with one shared
token and no encryption). No topology view, no plotting apart from the
QSTS chart. Technical limits: [reference](webui_reference.md#Limits).
