[![Documentation](https://github.com/Welthulk/Sparlectra.jl/actions/workflows/jekyll-gh-pages.yml/badge.svg)](https://welthulk.github.io/Sparlectra.jl/)


# SPARLECTRA
<a href="https://github.com/Welthulk/Sparlectra.jl/tree/main/"><img align="left" width="100" src="docs/src/assets/logo.png" style="margin-right: 20px" /></a>

This package contains tools for subsequent network calculations. It primarily features a program for calculating load flow using the Newton-Raphson method. The focus is to provide valuable insights into load flow calculations for both students and ambitious professionals.

## Features

- AC power flow with multiple internal Newton-Raphson formulations (`:polar_full`, `:rectangular`, `:classic`).
- Network Modeling: Comprehensive modeling of electrical network components
    Buses (PQ, PV, Slack)
    Transmission lines
    Transformers (2-winding and 3-winding)
    Generators and loads
    Shunts
- Transmission lines and transformers can be represented using π-equivalent models,
allowing direct use of CIM and manufacturer data without additional model conversions.    
- Topological bus links (`addLink!`) for busbar couplers and sectionalizer modeling.
- Canonical external solver interface (`PFModel`/`PFSolution`) to integrate third-party solvers.
- MATPOWER-compatible import/export utilities and local casefile helper workflow.

---

## Documentation

- User documentation: <https://welthulk.github.io/Sparlectra.jl/>
- API reference: <https://welthulk.github.io/Sparlectra.jl/reference/>

## Installation
```julia
using Pkg
Pkg.add("Sparlectra")
```

## Development Toolchain

This repository tracks a `Manifest.toml` generated with Julia `1.12.5`.
For development and CI parity, use Julia `1.12.5` (see `.julia-version`).

If you intentionally work with another Julia minor version, regenerate the
manifest locally with that version before running `Pkg.instantiate()`.

## Quick Start

```julia
using Sparlectra

# Ensure case file exists locally (downloads on demand into data/mpower)
case_path = ensure_casefile("case14.m")

# Build network and run Newton-Raphson power flow
net = createNetFromMatPowerFile(case_path, false)
ite, erg = runpf!(net, 10, 1e-6, 0)

if erg == 0
    calcNetLosses!(net)
    printACPFlowResults(net, 0.0, ite, 1e-6)
end
```

### Network Creation
This package supports the import and export of Matpower .m files, although currently it only reads bus, generator, and branch data from these files. Please note that additional Matlab functions within the .m file are not supported. Additionally, you can modify the imported Matpower files or you can create your own network using easy-to-use functions provided by the package.

### Test Cases and Benchmarks

Sparlectra does not ship third-party power system test cases by default.

Instead, MATPOWER-compatible case files (e.g. `case14`, `case118`) can be
downloaded **on demand** using helper scripts provided with the package.
Downloaded files are stored locally and are not part of the Sparlectra source
distribution.

This approach keeps the repository lightweight and avoids bundling external
data while still allowing reproducible experiments and benchmarks.


### License
This project is licensed under the Apache License, Version 2.0.
[The license file](LICENSE) contains the complete licensing information.








