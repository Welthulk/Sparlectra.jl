# [Import and Export Reference](@id reference_import_export)

```@autodocs
Modules = [Sparlectra]
Pages = [
    "import_export.jl",
    "exportMatPower.jl"
]
```

```@autodocs
Modules = [Sparlectra.FetchMatpowerCase]
Pages = [
    "FetchMatpowerCase.jl"
]
```

## Per-terminal branch status in import and export (r0.10.0)

- **CGMES:** `Terminal.connected` maps per terminal onto
  `from_status`/`to_status` on import. A line or two-winding transformer
  with exactly one open terminal becomes a partially open branch (its
  exact pi reduction acts at the closed bus); the earlier substitute
  shunt with half the charging is gone, and one-sided open transformers
  are no longer dropped. `Equipment.inService = false` still counts as
  both terminals open. The export writes `ACDCTerminal.connected` per
  terminal from the flags, so a roundtrip preserves the partial state.
- **MATPOWER:** `BR_STATUS` has no partial state. A partially open branch
  exports as `BR_STATUS = 0` plus its exact input admittance `Y_in` as a
  bus shunt (`GS`/`BS`) at the closed bus, with an `open_terminal=to` (or
  `from`) marker in the per-branch comment. A MATPOWER solve reproduces
  the same Y-bus; the roundtrip loses the partial state by design (the
  branch comes back out of service plus a shunt).
- **DTF:** outage cards identify whole branches; the format carries no
  per-end switch field (see the [DTF format notes](dtf_format.md)).
