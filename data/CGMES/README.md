# CGMES test deliveries (local cache)

This directory is a **local download cache**. It is empty in a fresh
checkout, ignored by Git (only this file is tracked), and nothing in the
package or the test suite depends on it: the tests run on checked-in cases
only, and CGMES deliveries from third parties are never part of the
repository.

The deliveries are useful for trying the CGMES importer on real data. Two
public sources are known to Sparlectra and can be fetched on demand:

- the ENTSO-E CGMES conformity test configurations (CGMES 2.4.15:
  MicroGrid, SmallGrid, FullGrid, RealGrid, MiniGrid, the PST sets), one
  package that is downloaded once and extracted here;
- the ENTSO-E ReliCapGrid models (CGMES 3.0: `svedala`, `espheim`,
  `belgovia`, `galia`, `britheim`, `nordheim`, `portheim`, the combined
  `relicapgrid_cgm` and `svedala_neighbours`), fetched file by file from
  GitHub and packed into one ZIP per alias.

```julia
using Sparlectra
Sparlectra.CGMESImporter.allCGMESTestSetAliases()          # every alias
zip = Sparlectra.CGMESImporter.fetchCGMESTestSet("microgrid_be"; outdir = mktempdir())
res = importCGMES(path = zip, name = "microgrid_be")
```

`fetchCGMESTestSet` downloads on first use and reuses the cache afterwards;
a ReliCapGrid ZIP is packed only when every member file (grid profiles,
boundary files, commonData) is in the cache, and a ZIP from an older cache
layout is refreshed. The cache location defaults to this directory and can
be moved with the environment variable `SPARLECTRA_CGMES_CACHE`.

The Web UI accepts the same ZIPs and folders as case input (see the CGMES
import page of the documentation). Both sources carry their own licenses
and terms; check them before redistributing anything from here.
