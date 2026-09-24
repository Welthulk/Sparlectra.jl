# Security

## Supported install paths

- `Pkg.add("Sparlectra")` from the Julia General registry. The registry pins
  the git tree hash of every registered version, and Pkg verifies it on
  install. This is the recommended path for the library.
- A GitHub release (tag `vX.Y.Z`) for the application under `app/`: the
  source archives, the SBOM (`Sparlectra.spdx.json`) and the one-line
  installers. Every release carries `SHA256SUMS` over these files.
- A clone of the `main` branch for development. `main` is protected against
  force pushes.

Only the latest release is supported.

## Verifying a download

- Compare a downloaded file with `SHA256SUMS` of the release:
  `sha256sum -c SHA256SUMS` on Linux and macOS, `Get-FileHash` on Windows.
- The one-line installers fetch `tools/install_webui.sh` and
  `tools/install_webui.bat` from `main`; their checksums at the release tag
  are in `SHA256SUMS`, so a copy can be checked before it is run.

## Reporting a vulnerability

Please report security issues privately through GitHub's
"Report a vulnerability" form on the Security tab of this repository, not
in a public issue. You will get an answer within a few days.
