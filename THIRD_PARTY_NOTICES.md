# Third-Party Notices (Data & Test Cases)

This repository does **not** ship third-party power system test cases by default.

Instead, external test cases (e.g. MATPOWER-compatible `.m` files) may be
**downloaded on demand** by helper scripts provided in this repository.

## General Principles
- Third-party data is **not redistributed** as part of the source repository.
- Any downloaded data remains subject to its original license terms.
- Users are responsible for ensuring compliance with the applicable licenses
  of downloaded datasets.
- The repository provides tooling for reproducible retrieval, not bundled data.

## Downloaded Data
Downloaded test cases are stored locally (outside version control) under:

    data/mpower/

This directory is intended for **local, user-managed data only** and is not
meant to be committed to the repository.

## Provenance
Download scripts may record basic provenance metadata (e.g. source URL,
retrieval date, hash) to aid traceability. Such metadata does not imply
license verification or endorsement.
