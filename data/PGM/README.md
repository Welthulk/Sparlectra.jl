# Plain power-grid-model input files

Cases in the plain power-grid-model (PGM) JSON input format, without the
`sparlectra` block of the Sparlectra Case Format. They exercise the
import path that reads a PGM dataset as delivered by other tools; the
Case Format counterpart of each file lives under `data/scf`.

- `feeder3_hardpv_pgm.json`: the three-bus 110 kV feeder of the workshop
  tour (chapter 5) as a plain PGM dataset: three nodes, two lines, a
  source at B1 (`u_ref 1.02`), a 10 MW generator at B2 with a voltage
  regulator, a load at B3. `feeder3_hardpv_pgm.config.yaml` is its case
  configuration sidecar (rectangular solver, `tol 1e-6`, external-grid
  feeder enabled). `data/scf/feeder3_hardpv_vsPGM.scf.json` carries the
  same data part with the Sparlectra additions, for a side-by-side
  comparison of what the Case Format adds.
