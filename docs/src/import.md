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
safe recommendations unless the YAML file already contains an explicit value for
that option. The pre-run checks transformer ratio and phase-shift conventions,
bus-shunt residual impact, PV/REF voltage-source handling, large-case flat-start
strategy, and conservative PV→PQ Q-limit settings. For large-case runtime
profiles it keeps start projection disabled, uses DC-angle and blended-voltage
flat-start seeds, recommends practical validation tolerance, and disables
expensive diagnostics unless explicitly requested. The terminal summary and
logfile list the evidence, selected recommendations, applied options, and
explicit YAML values that were preserved.

### Import Parser

The `casefileparser` function parses Matpower case files and returns the raw data arrays:

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
