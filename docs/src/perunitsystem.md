# Per Unit System

The Per Unit System is a method for simplifying power and voltage calculations in electrical networks by referencing all values to a standardized base quantity, rendering them dimensionless in the result.

In the Per Unit System, nominal values of voltage, current, and power are used as units, making the absolute size of components irrelevant. Instead, the emphasis is on ratios and relationships between components.

The Per Unit System facilitates the comparison and analysis of networks of different scales as it is independent of specific nominal values and relates to standardized base quantities. Calculations involving square root of 3 are often avoided.

## Base Quantities

- **Sbase**: The pivotal parameter is the reference or base apparent power.
  - Determined by network design, the most powerful generator, or is predefined (usually 1000 MVA).
  - The base apparent power is constant across all regions of the network.

- **Ubase**: Nominal voltage, reference voltage (Base voltage).

- **Ibase**: Reference current, derived quantity.

- **Zbase**: Nominal impedance, reference impedance, derived quantity.

## Calculations

- \( I_{\text{Base}} = \frac{S_{\text{Base}}}{\sqrt{3} \cdot U_{\text{Base}}} \)
- \( Z_{\text{Base}} = \frac{{U_{\text{Base}}}^2}{S_{\text{Base}}} \)
- \( u_\text{pu} = \frac{U_{\text{akt}}}{U_{\text{Base}}} \)
- \( i_\text{pu} = \frac{I_{\text{akt}}}{I_{\text{Base}}} = \frac{I_{\text{akt}}}{S_{\text{Base}}/\sqrt{3} \cdot U_{\text{Base}}} \)
- \( z_\text{pu} = \frac{Z}{Z_{\text{Base}}} \)
- \( s_\text{pu} = \text{upu} \cdot \text{ipu} \)
- \( p_\text{pu} = (\text{ipu})^2 \cdot \text{rpu} \)

## Example

- A-B: 10,000 kVA, 13.8/138 kV, Leakage Resistance 10%
- B-C: 10,000 kVA, 138/69 kV, Leakage Resistance 8%

### Region C

- \( R = 300 \, \Omega \), \( U_{\text{Base}} = 69 \, \text{kV} \), \( S_{\text{Base}} = 10,000 \, \text{kVA} \)
- \( Z_{\text{Base}} = \frac{{69 \, \text{kV}}^2}{10,000 \, \text{kVA}} = 476 \, \Omega \)
- \( \text{zpu} = \frac{300}{476} = 0.63 \, \text{p.u.} \)

### Region A, B

- \( \text{zpu(A-B)} = j0.1 \, \text{p.u.} \)
- \( \text{zpu(B-C)} = j0.08 \, \text{p.u.} \)

## Equivalent Circuit

- \( \text{zges,pu} = j(0.1 + 0.08) + 0.63 \, \text{p.u.} = 0.63 \, \text{p.u.} + j0.18 \, \text{p.u.} \)
- Through the Per Unit System, all voltage levels have been eliminated (transformers are eliminated)!

## Change of Reference Voltage

If the reference voltage and reference nominal apparent power change, the impedances can be converted as follows:

\[ \text{zpu}_{\text{new}} = \text{zpu}_{\text{old}} \cdot \left(\frac{U_{\text{Base,old}}}{U_{\text{Base,new}}}\right)^2 \cdot \left(\frac{S_{\text{Base,new}}}{S_{\text{Base,old}}}\right) \]

## Graphics

- [Graphic 1: Example Network](tbd1.png)
- [Graphic 2: Equivalent Circuit](tbd2.png)
