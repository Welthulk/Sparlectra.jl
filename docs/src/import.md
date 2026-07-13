# Import and Export Network Configuration

Sparlectra provides functionality to import network models from Matpower case files and export networks to Matpower format. This page documents the functions for reading, importing, and exporting network configurations.

## Importing Matpower Files

### Basic Import

The simplest way to import a Matpower case file is using the `run_sparlectra` function, which both imports the network and runs a power flow analysis:

```julia
using Sparlectra

case_path = ensure_casefile("case14.m")

# Import a MATPOWER case file and run the configured framework workflow.
result = run_sparlectra(
    casefile = basename(case_path),
    path = dirname(case_path),
)

net = result.net
```

### Configuration-driven MATPOWER batches

Keep `run_sparlectra` as the single-case API when one `SparlectraRunResult` is
expected. For deterministic sequential execution of several configured cases,
list them in YAML and call `run_sparlectra_cases`:

```yaml
matpower_import:
  cases: [case14.m, case118.m]
```

```julia
using Sparlectra

cfg = load_sparlectra_config("path/to/configuration.yaml")
results = run_sparlectra_cases(config = cfg)

for result in results
    println(result.net.name, ": ", result.outcome)
end
```

Configured order is preserved. A non-empty `matpower_import.cases` list wins
over `matpower_import.case`; the single-case setting is used as a fallback when
the list is empty. Bare standard MATPOWER case names are resolved through
`ensure_casefile` and downloaded on demand. Pass `path` when the listed cases
are local fixtures or site-specific files. Batch-level performance-profile
aggregation is intentionally unsupported; profile individual `run_sparlectra`
calls instead.

### Import Without Power Flow

If you want to import a network without immediately running power flow, you can use the `createNetFromMatPowerFile` function:

```julia
using Sparlectra

# Path to the Matpower file
file_path = "path/to/case5.m"

# Import the network configuration
net = createNetFromMatPowerFile(filename = file_path, log = false)
```

### Ensuring Test Cases Exist Locally

When working with standard MATPOWER test cases (for example `case14.m`), you can
let Sparlectra download the file on demand:

```julia
using Sparlectra

# Downloads to Sparlectra.data/mpower if missing and returns absolute path
case_path = ensure_casefile("case14.m")
net = createNetFromMatPowerFile(filename = case_path, log = false)
```

MATPOWER transformer convention options are configured through `SparlectraConfig.matpower` for framework runs or passed directly to the importer for import-only calls. The default `matpower_ratio = "normal"` uses branch `TAP` values as
stored, while `matpower_ratio = "reciprocal"` imports their reciprocal:

```julia
net = createNetFromMatPowerFile(
    filename = case_path,
    matpower_shift_unit = "rad",
    matpower_shift_sign = -1.0,
    matpower_ratio = "reciprocal",
)
```

You can also request a generated `.jl` companion file by passing a `.jl` name:

```julia
case_jl_path = ensure_casefile("case14.jl")
```

### MATPOWER auto-profile pre-run

The `examples/matpower_import.jl` workflow can run a compact diagnostic
pre-run before the main power-flow solve:

It also supports explicit Julia-thread startup control for script workflows:

```bash
julia --threads=8 --project=. examples/matpower_import.jl
julia --project=. examples/matpower_import.jl --julia-threads=8
```

```powershell
$env:SPARLECTRA_JULIA_THREADS = "8"
julia --project=. examples/matpower_import.jl
```

`runtime.julia_threads` is resolved after CLI/env overrides; if the requested
value differs from active `Threads.nthreads()`, the script can re-exec once
with the requested `--threads` setting.

```yaml
matpower_import:
  auto_profile: recommend  # off, recommend, or apply
  auto_profile_log: true
```

`recommend` logs a profile without changing the active run. `apply` applies
only safe import-convention recommendations; solver-start and Q-limit
recommendations are logged as skipped recommendations unless the corresponding
runtime options are configured directly. The pre-run checks transformer ratio
and phase-shift conventions, bus-shunt residual impact, PV/REF voltage-source
handling, large-case flat-start strategy, and conservative PV→PQ Q-limit
settings. For large-case runtime profiles it keeps start projection disabled,
uses DC-angle and blended-voltage
flat-start seeds, recommends practical validation tolerance, and disables
expensive diagnostics unless explicitly requested. The terminal summary and
logfile list the evidence, selected recommendations, applied options, and
explicit YAML values that were preserved.

### Import Parser

The `casefileparser` function parses Matpower case files and returns the raw data arrays:

### MATPOWER metadata and DC lines

By default Sparlectra preserves the historical numeric bus naming behavior for
MATPOWER imports. Set `matpower_import.apply_bus_names: true` (or pass
`apply_bus_names = true`) to use the standard optional MATPOWER `mpc.bus_name`
field as imported bus names while preserving original `BUS_I` values in
`busOrigIdxDict`. Sparlectra also supports user-defined MATPOWER-compatible
metadata fields `mpc.branch_name`, `mpc.branch_kind`, and
`mpc.for001_contingencies`. Enable them with
`apply_branch_names = true`, `apply_branch_kind = true`, and
`import_for001_contingencies = true`; branch names/kinds are attached as
`net.matpower_branch_metadata`, and FOR001 contingency names are copied to
`net.for001Contingencies`. `mpc.branch_kind` values `L`, `LINE`, and `ACL`
force line import, while `T`, `TRAFO`, `TRANSFORMER`, and `2WT` force
transformer import. Missing or length-mismatched metadata falls back to the
legacy electrical heuristic.

FOR/DTF transformer records can carry active no-load conductance (`G`) and
reactive shunt susceptance (`B`) that standard MATPOWER branch rows cannot
represent with transformer-local semantics. The native DTF importer converts
nonzero transformer `G` to equivalent bus shunt conductance split equally across
the transformer terminals, preserving the original branch label, terminal names,
raw values, converted per-unit values, tap/stage data, and an explicit allocation
marker in branch metadata. MATPOWER export writes this richer data in a
versioned `mpc.sparlectra.transformer_losses` extension block and adds a visible
warning comment near the top of the `.m` file. Plain MATPOWER may ignore this
block, so transformer active losses can differ outside Sparlectra; Sparlectra
reimports the block and avoids adding the equivalent `Gs/Bs` shunts twice when
they are already present in the standard bus table.

`matpower_import.matpower_dcline_mode` controls `mpc.dcline` handling. The
default `:pf_injections` imports active DC-line rows as fixed terminal injections. Explicit `:reject_active` remains available for strict validation and
(`status != 0`) abort before solving with `failure_reason =
unsupported_matpower_dcline`, while empty or inactive-only tables are ignored.
`:ignore_inactive` has the same active-row rejection behavior and documents the
inactive-row ignore policy. `:pf_injections` uses a MATPOWER
`toggle_dcline`-compatible power-flow approximation by adding two
generator-like terminal prosumers for each active row. The from terminal uses
`F_BUS`, `PF`, `QF`, `VF`, and `QMINF`/`QMAXF` with active injection
`PG = -PF`. The to terminal uses `T_BUS`, `QT`, `VT`, and
`QMINT`/`QMAXT`; received active power is recomputed as
`PF - (LOSS0 + LOSS1 * PF)` when `LOSS0`/`LOSS1` are present, otherwise the
input `PT` column is used. Terminal buses with voltage setpoints are treated as
voltage-controlled where MATPOWER would make them PV, without demoting
reference buses or activating isolated buses. API/Web UI runs write a compact
`matpower_dcline.csv` artifact describing the mapping. This is not a full HVDC
converter or DC-grid model; OPF constraints, converter controls,
`dclinecost`, and DC-line optimization remain unsupported. DC terminal injections do not create Ybus branches, so disconnected AC components are detected and solved island-wise by default; this is a power-flow approximation, not a complete HVDC converter or DC-grid model.

```julia
using Sparlectra

file = "case9.m"

# Parse the file and get raw data
case_name, baseMVA, busData, genData, branchData = casefileparser(file)

# Now you can work with the raw data arrays
println("Case name: \$case_name")
println("Base MVA: \$baseMVA")
println("Number of buses: \$(size(busData, 1))")
println("Number of generators: \$(size(genData, 1))")
println("Number of branches: \$(size(branchData, 1))")
```

## Exporting Networks to Matpower Format

You can export a Sparlectra network to a Matpower case file using the `writeMatpowerCasefile` function:

```julia
using Sparlectra

# First create or import a network
net = Net(name = "export_example", baseMVA = 100.0)
# Add components to the network...

# Export the network to a Matpower case file
filepath = "path/to/output/export_example.m"
writeMatpowerCasefile(net, filepath)
```

## Running Power Flow on Imported Networks

After importing a network, you can run a power flow analysis:

```julia
using Sparlectra

# Import a network from a Matpower file
net = createNetFromMatPowerFile(filename = "case5.m", log = false)

# Run power flow
tol = 1e-6
max_ite = 10
verbose = 0

# Run power flow
ite, erg = runpf!(net, max_ite, tol, verbose)

# Check results and calculate losses
if erg == 0
  calcNetLosses!(net)
  printACPFlowResults(net, 0.0, ite, tol)
else
  @warn "Power flow did not converge after \$ite iterations"
end
```

## Working with Custom File Paths

You can specify custom paths for import and export operations:

```julia
using Sparlectra

# Import from a specific path
path = "C:/Users/YourUsername/Documents"
file = "custom_case.m"
cfg = load_sparlectra_config("examples/configuration.yaml")
result = run_sparlectra(casefile = file, path = path, config = cfg)
net = result.net

# Export to a specific path
output_path = joinpath(path, "exported_case.m")
writeMatpowerCasefile(net, output_path)
```

## Matpower File Format Notes

- Sparlectra currently reads bus, generator, and branch data from Matpower files
- Additional Matlab functions within the .m file are not supported
- The Matpower format version supported is version 2

!!! Note 
    Please note that providing support for individualized network data issues is beyond the scope of this project, as it is not maintained by an organization. Users are encouraged to take initiative in resolving such issues independently and sharing their results with the community.

## Web UI case-file import

The local Web UI can copy downloaded case files into its case storage from the PowerFlow page. Use **Import case files** to open the browser's native picker, choose one or more MATPOWER `.m`/`.M` files or DFT `.dat`/`.DAT` files, and submit the import form. The import form is separate from **Start PowerFlow run**: uploading files only stores them and refreshes the selector; it does not parse uploaded MATPOWER code or start a calculation.

Imported files go to the same effective case directory used by the Web UI case selector. In a writable repository checkout this is the Web UI `data/mpower` case directory; in installed-package or restricted launch contexts the Web UI uses its user-writable application data `data/mpower` directory next to the run output root. The manual full-path field remains available for advanced workflows.

The Web UI rejects unsupported extensions such as `.jl`, `.zip`, `.yaml`, executables, scripts, and files without an extension. Uploads are limited to 100 MiB per file and 250 MiB per request. Existing files are not overwritten: the import result reports `already exists` for conflicts and still imports other valid files from the same selection. FOR002 reference `.DAT` files remain hidden from the normal runnable selector even when copied into the case directory.
