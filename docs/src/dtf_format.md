# DTF legacy input format

Sparlectra reads the fixed-column DTF files of the legacy Testnetz13
validation examples natively, without routing through MATPOWER:
`DTFImporter.read_dtf` parses a file, `DTFImporter.build_net` or
`createNetFromDTFFile` builds the `Net`. The FOR001-specific information
(outage cards, trailing records, transformer control fields, nameplate
voltages, branch identity metadata) stays on typed records, raw lines
included, for audit. FOR002 text reports are validation references, not
model input.

## File structure

| Card | Fields | Notes |
|---|---|---|
| Parameter/text cards | legacy run parameters, free-text descriptors | Kept for diagnostics only. |
| Nominal-voltage card | voltage bases referenced by voltage-level indices | Branch R/X/G/B per-unit conversion uses the branch voltage-level index; a 231 kV level is not collapsed to 230 kV. |
| Size/count card | bus, branch, compensation and transformer-control counts, named slack bus | The counts split the fixed-column sections before optional outage or trailing data is parsed. |
| Branch cards (`L`, `T`) | lines and transformers | `T` marks a transformer, every other kind is an AC line. Transformer transverse conductance and susceptance become the PI branch fields `g_pu` and `b_pu` (branch shunt arms, not node shunts). |
| Compensation cards | fixed compensation records of the base network | |
| Transformer-control cards | winding/nameplate voltages, longitudinal tap range and step, optional Schraegregler skew-angle fields | Matched to the transformer by terminals and parallel identifier. The metadata holds the parsed winding values, ratio convention, neutral base ratio, relative tap fraction, skew angle, effective ratio and effective phase shift. |
| Bus cards | bus type, voltage-level index, name, start voltage, angle, load, generation | |
| Outage records | between `AUSFALL` and `ENDE`, after the bus section | Optional, see below. |
| Trailing branch-echo records | a standalone `A` marker followed by a branch-like `L` or `T` record | Duplicates existing branch information for legacy diagnostics; kept as metadata, neither an additional branch nor an outage. |

## Outage records

Every non-empty record between `AUSFALL` and `ENDE` identifies one branch
outage by branch kind, voltage-level index, parallel identifier, from-bus
name and to-bus name. Applying one requires an unambiguous match to exactly
one base-network branch; a missing or ambiguous match produces a diagnostic
and removes nothing. Branch cards carry no per-end switch field, so a
partially open branch cannot be expressed: the aggregate branch status is
used.

## PV, PQ and slack bus interpretation

The slack bus comes from the size card and becomes a Sparlectra slack
generator; PV-style buses keep their voltage-regulating role; PQ
generator/load buses keep fixed injections and are not promoted to PV
because non-zero generation appears on the card.

## Transformer ratio convention

Winding values such as 400/231 kV are nameplate metadata, not a permanent
off-nominal tap; tap positions and skew angles apply as deviations from
neutral. `transformer_ratio_mode` selects the neutral base ratio:

- `:neutral_one` (default, the observed FOR002 legacy semantics): the base
  ratio at neutral tap is `1.0`.
- `:winding_over_network`: the winding/nameplate ratio divided by the network
  nominal-voltage ratio, for example `(400/231) / (400/230)`.

The nameplate ratio stays visible in the metadata fields
`nominal_unregulated_kv`, `nominal_regulated_kv`, `from_bus_vn_kV`,
`to_bus_vn_kV`, `winding_over_network_base_ratio`, `transformer_ratio_mode`,
`base_ratio_used`, `tap_fraction`, `skew_angle_deg`, `effective_ratio` and
`effective_shift_deg`.

With `model.tap_changer_model = impedance_correction` (see [Transformer
tap-changer model](configuration.md#transformer-tap-changer-model)) the
`tap_fraction`/`skew_angle_deg` values are also folded into the transformer
series impedance (`Branch.r_pu`/`Branch.x_pu`), scaling R and X with
`|1 + f·e^(jφ)|²`; a subsequent `writeMatpowerCasefile` export writes the
corrected values, not the raw card impedance.

**Schraegregler skew angle.** The 60-degree field of the Schraegregler
example is the skew angle of the regulating voltage, not the final phase
shift: the complex tap follows from tap range, tap position and skew angle,
converted to the from-side off-nominal convention (reciprocal magnitude,
negative regulating-vector angle).

## Transformer shunt conductance and losses

DTF transformer `G` is stored as `Branch.g_pu` next to `Branch.b_pu`,
`r_pu`, `x_pu`, tap ratio and phase shift, not as terminal bus shunts.
`calcNetLosses!` sums branch-end powers (`S_from + S_to`), so total losses
include the longitudinal `R` and the voltage-dependent `G` losses; the
separate I²R helper reports the longitudinal component only. The
MATPOWER transformer-loss extension ([MATPOWER cases](matpower.md))
preserves `g_pu` across export/reimport round trips; standard MATPOWER
readers ignore it.

## Validation against FOR002

The validation examples compare native solves against FOR002 reports for
cases A-E in both `:neutral_one` and `:winding_over_network` mode.

## Web UI `.DAT` roles

The Web UI classifies `.DAT` uploads by content (`FOR002.DAT`-style names
are hints only) before offering them in the PowerFlow selector:

| Role | Meaning |
|---|---|
| `dtf_network_case` | Complete DTF network case with mandatory cards and counts; selectable as the primary case. |
| `dtf_network_case_with_outages` | Also runnable; exposes the embedded outage records to the outage controls. |
| `dtf_outage_or_reference` | Kept for the DTF outage/reference controls, not runnable as a primary case. |
| `unknown_dat` | Kept only when it passes the upload safety checks, never offered as a runnable case. |
