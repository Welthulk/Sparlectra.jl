# Local MATPOWER-Compatible Test Cases

This directory is used for **locally downloaded** MATPOWER-compatible test cases.

## Important
- This directory is **not populated by default**.
- Files placed here are **not part of the repository distribution**.
- Contents are expected to be excluded from version control (e.g. via `.gitignore`).

## Usage
Test cases may be downloaded automatically by example scripts or helper
functions (e.g. `FetchMatpowerCase.ensure_matpower_case`) when required.

Typical contents:

    data/mpower/
        ├─ case14.m
        ├─ case14.jl
        ├─ case118.m
        └─ ...

Both original files and locally generated derivatives (e.g. `.jl` conversions)
are considered **local artifacts**.

## License Responsibility
The license and usage terms of any downloaded test case are determined by its
original source. Users are responsible for reviewing and complying with those
terms before use or redistribution.
