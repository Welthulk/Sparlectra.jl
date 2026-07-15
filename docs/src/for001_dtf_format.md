# FOR001 / DTF legacy input format

Sparlectra includes a native reader for the legacy DTF/FOR001 files used by K. F. Schaefer's Testnetz13 examples. The importer is intended to reproduce the legacy network-model semantics used by the matching FOR002 ground-load-flow reports while preserving source fields for audit and diagnostics.

## Overview

A FOR001/DTF file is a fixed-column power-flow input file. Sparlectra parses it directly with `DTFImporter.read_dtf` and builds a `Net` with `DTFImporter.build_net` or `createNetFromDTFFile`. The native path does not route through MATPOWER and keeps FOR001-specific information such as outage cards, trailing records, transformer control fields, nameplate voltages and branch identity metadata.

## File structure

The Schaefer files follow this high-level order:

1. parameter/card text records,
2. a voltage-level card,
3. a size card,
4. branch cards,
5. optional compensation cards,
6. transformer-control cards,
7. bus cards,
8. optional outage records and trailing branch-echo records.

The parser preserves raw lines on typed records so diagnostics can compare the electrical model against the original fixed-column payload.

## Voltage-level card

The voltage-level card lists the nominal voltage levels used by subsequent branch and bus records. Branch impedance and admittance conversion uses the branch voltage-level index as the per-unit reference. Sparlectra does not globally collapse a 231 kV level to 230 kV in the default DTF path.

## Size card

The size card gives the expected counts for buses, branches, compensation records and transformer-control records, plus the legacy slack-bus name. Sparlectra uses these counts to split the fixed-column sections before it parses optional outage or trailing data.

## Branch cards

Branch cards describe lines and transformers. The card kind identifies transformer records with `T`; other branch records are treated as AC lines. The branch voltage-level index controls R/X/G/B per-unit conversion. Transformer transverse conductance and susceptance are represented by the native PI branch fields `g_pu` and `b_pu`; the resulting `G/2` arms are internal branch shunt arms, not separate node shunts.

## Transformer-control cards

Transformer-control cards carry the transformer winding/nameplate voltages, longitudinal tap range, tap step, and optional Schraegregler/skew-angle fields. Sparlectra matches a control card to its transformer by terminals and parallel identifier. The resulting metadata includes the parsed winding values, selected ratio convention, neutral base ratio, relative tap fraction, skew angle, effective ratio and effective phase shift.

## Bus cards

Bus cards define the legacy bus type, voltage-level index, name, start voltage, angle, load and generation quantities. Sparlectra maps the Schaefer slack, PV and PQ conventions into its internal bus and prosumer model without turning every non-zero generator row into a regulating bus.

## Outage cards and trailing A branch-echo records

Outage sections are parsed as structured outage metadata. In the Schaefer B-E files, post-bus standalone `A` markers followed by branch-like `L` or `T` records are branch-echo records, not additional electrical branches and not outages. Sparlectra preserves them for diagnostics and does not add them to the network topology.

## PV, PQ and slack bus interpretation

The slack bus is taken from the DTF size/slack information and mapped to a Sparlectra slack generator. PV-style buses keep their voltage-regulating role. PQ generator/load buses retain fixed injections and are not promoted to voltage-regulating buses only because non-zero generation appears on the card.

## Transformer ratio convention

In DTF/FOR001 compatibility mode, transformer winding values such as 400/231 kV are preserved as nameplate/control metadata. They are not automatically interpreted as a permanent off-nominal transformer tap at neutral position.

At neutral tap position, Sparlectra's DTF importer uses neutral-one behaviour by default, matching the observed FOR002 legacy semantics. Tap positions and Schraegregler/skew-angle controls are applied as deviations from this neutral position.

The executable `transformer_ratio_mode` values are:

- `:neutral_one` — the default DTF/FOR001 compatibility mode. The base ratio at neutral tap is `1.0`; longitudinal and skew-angle controls are applied relative to neutral.
- `:winding_over_network` — the previous Sparlectra diagnostic/compatibility interpretation. The neutral base ratio is the winding/nameplate ratio divided by the network nominal-voltage ratio, for example `(400/231) / (400/230)`.

The nameplate/winding ratio remains visible in metadata such as `nominal_unregulated_kv`, `nominal_regulated_kv`, `from_bus_vn_kV`, `to_bus_vn_kV`, `winding_over_network_base_ratio`, `transformer_ratio_mode`, `base_ratio_used`, `tap_fraction`, `skew_angle_deg`, `effective_ratio` and `effective_shift_deg`.

## Schraegregler / skew-angle transformer controls

The DTF 60-degree field in the Schraegregler example is a skew angle of the regulating voltage, not the final transformer phase-shift angle. Sparlectra computes the resulting complex tap from tap range, actual tap position and skew angle. The complex regulating vector is converted to the from-side off-nominal transformer convention by using its reciprocal magnitude and the negative regulating-vector angle.

## Known compatibility notes

The native DTF path is deliberately format-aware. It preserves winding/nameplate data for diagnostics even when neutral-one mode does not apply that value as a neutral off-nominal tap. Branch R/X/G/B conversion remains tied to the DTF branch voltage-level index. FOR002 text reports are treated as legacy validation references, not as additional model input.

## Validation against FOR002

The Schaefer validation examples compare native DTF solves against FOR002 reports for cases A-E. The report workflow executes both `:neutral_one` and `:winding_over_network` transformer ratio modes so voltage bias, slack injections, losses and branch-flow anchors can be compared as true alternative network builds.

## Transformer shunt conductance and losses

DTF transformer `G` is stored on the Sparlectra transformer branch as `Branch.g_pu`, alongside `Branch.b_pu`, `r_pu`, `x_pu`, tap ratio and phase shift. It is not converted to separate terminal bus shunts. `calcNetLosses!` sums branch-end powers (`S_from + S_to`), so total branch losses include both longitudinal `R` losses and voltage-dependent branch-shunt `G` losses. The separate I²R helper reports only the longitudinal copper/winding component. Sparlectra's proprietary MATPOWER transformer-loss extension preserves `g_pu` for Sparlectra export/reimport round trips; standard MATPOWER readers ignore that extension.
