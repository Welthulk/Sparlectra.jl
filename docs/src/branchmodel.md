# Network Branch Model

Each branch is treated with the same four-terminal network model. It is a four-terminal network with an ideal transformer connected upstream. For power lines, the transmission ratio $N$ is set to 1. For transformers, the transformation ratio $N$ is given as a complex value.

The branch admittance matrix $Y_{br}$ is:

```math
Y_{br} = \begin{bmatrix}
    \frac{1}{\tau^2} \left( y_{ser} + \frac{y_{shunt}}{2} \right) &
    -y_{ser} \frac{1}{\tau e^{-j\phi}} \\
    -y_{ser} \frac{1}{\tau e^{j\phi}} &
    y_{ser} + \frac{y_{shunt}}{2}
\end{bmatrix}
```

where:
- `y_ser` is the series admittance,
- `y_shunt` is the total branch shunt admittance,
- `R` is the resistance component,
- `X` is the reactance component,
- `G` is the conductance component,
- `B` is the susceptance component,
- $N$ is the complex transformation factor (e.g. 1 for power lines).

```math
N = \tau e^{j\phi}
```

```math
y_{ser} = \frac{1}{R + jX}
```

```math
y_{shunt} = G + jB
```

## Y-bus diagonal entries (π-model + shunts)

For a node $i$, the diagonal Y-bus entry is the nodal self-admittance:

```math
Y_{ii} = \sum_{k \in \mathcal{N}(i)} y_{ik} + y_i^{sh}
```

with:
- variable $y_{ik}$: series admittance contribution of branch $i-k$
- variable $y_i^{sh}$: explicit shunt admittance at bus $i$
  
  For a π-model branch $i-k$, the local diagonal stamp is:

```math
Y_{ii} \mathrel{+}= y_{ik} + \frac{y_{ik}^{sh}}{2}
```

and the off-diagonal relation is:

```math
Y_{ik} = -y_{ik}
```

Hence (without explicit shunts):

```math
Y_{ii} = -\sum_{k \neq i} Y_{ik}
```

### Interpretation

- Variable $Y_{ii}$ is the full self-admittance seen at bus $i$
- it combines network coupling (series paths) and local shunts
- it represents current leaving bus $i$ for $V_i = 1\,\mathrm{pu}$ in nodal form

### Practical notes

- In typical grids this supports diagonal dominance and helps numerical robustness of NR/SE workflows
- The real part is usually non-negative
- The imaginary part reflects the balance of inductive series effects and capacitive or inductive shunt terms

### Implementation in Sparlectra

- branch builders (`addACLine!`, `addPIModelACLine!`, `addPIModelTrafo!`) stamp series admittance plus half shunt on each side according to the branch model
- explicit shunts are added as nodal shunt terms when `bus_shunt_model = "admittance"`

## Bus-shunt modeling modes

Sparlectra supports two representations for real bus shunts imported from sources such as MATPOWER `Gs`/`Bs` columns:

- mode `"admittance"`: the classical treatment. The bus shunt admittance $y_i^{sh} = G_i + jB_i$ is stamped into the Y-bus diagonal as part of $Y_{ii}$. This is the default mode and preserves existing numerical behavior.
- mode `"voltage_dependent_injection"`: the bus shunt is not stamped into Y-bus. Instead, its power is evaluated in the nonlinear injection/mismatch path as a local voltage-dependent term.

For a bus shunt admittance $y_i^{sh}$ and local voltage magnitude $|V_i|$, Sparlectra uses the sign convention:

```math
S_i^{sh} = |V_i|^2 \overline{y_i^{sh}}
```

This follows the same sign convention as admittance stamping: a positive conductance contributes positive active power, while the reactive sign follows the complex conjugate of the shunt admittance. The rectangular mismatch uses `S_calc - S_spec`, so voltage-dependent injection mode subtracts $S_i^{sh}$ from the specified net injection. This keeps the equations equivalent to the admittance model while avoiding double-counting: each bus shunt is either stamped into Y-bus or represented as a voltage-dependent injection term, never both.

The voltage-dependent injection mode is useful when a solver formulation wants the network admittance matrix to contain only branch/network coupling while keeping shunt effects in the nonlinear injection equations.

## Phase-shift control direction (practical probe)

For transformer-flow control with a phase shifter, do not hard-code the sign of the controller action. Probe it on the active model:

1. Compute `P_ab(phi = 0 deg)`
2. Compute `P_ab(phi = +5 deg)`
3. Evaluate `Delta_P_ab = P_ab(5 deg) - P_ab(0 deg)`
4. Use the measured direction in control:
   - if `P_ab < target`, move `phi` in the direction that increases `P_ab`
   - otherwise move `phi` in the opposite direction

See executable example:

`src/examples/example_transformer_phase_shift_control.jl`

## Circuit diagram

```
                    y_ser
      x--┓┏---------###----------x
         ||   |             |
         ||   # y_shunt     # y_shunt
         ||   |             |
      x--┛┗----------------------x
         N = tau * exp(j * phi)
```
