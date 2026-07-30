# Short-circuit reference table (developer validation record)

Reference values for `runShortCircuit!` (balanced 3-phase Ik'' per
IEC 60909-0). The table is the source of record for
`examples/others/exp_short_circuit_reference.jl`, which recomputes every
quantity with Sparlectra and reports PASS/FAIL against the values below
(rtol 1e-6, part of `examples/run_short_circuit_suite.jl`).

Provenance of the reference column: **analytic hand evaluation of the IEC
60909-0 formulas** — derivations below, arithmetic independent of the
Sparlectra implementation. The external-tool column is deliberately left
open: no third-party tool is installed or downloaded from this repository.
The pandapower snippet at the end reproduces the same four cases; run it in
any environment that has pandapower and enter the numbers here.

## Case definitions

All cases use `baseMVA = 100`, IEC 60909-0 Table-1 voltage factors
(`c_max = 1.10`, `c_min = 1.00` above 1 kV), and the method-b peak factor
`κ = min(1.15·(1.02 + 0.98·e^(−3·R/X)), 2.0)`.

- **R1 — network feeder, connection point.** One 110 kV bus `Q`; feeder with
  declared `Ik''max = 10 kA`, `Ik''min = 8 kA`, `R/X = 0.1`.
  Derivation: `Z_Q = c·U_n/(√3·Ik_decl)`; at the connection point the whole
  short-circuit impedance is `Z_Q`, so `Ik'' = Ik_decl` exactly (the c
  factor cancels — the IEC semantics of the declared feeder current).
  `Sk'' = √3·U_n·Ik''`; `R/X = 0.1 → κ_b = 1.02+0.98·e^(−0.3) = 1.746002…`,
  `1.15·κ_b = 2.0079… → cap 2.0`, `i_p = 2.0·√2·Ik''`.
- **R2 — feeder + line.** R1's feeder bus `Q` plus a line to bus `A`
  (both 110 kV): `r = 0.05 pu`, `x = 0.5 pu` on `Z_base = 121 Ω`
  → `Z_L = 6.05 + j60.5 Ω`. Feeder split: `X_Q = Z_Q/√(1+0.1²)`,
  `R_Q = 0.1·X_Q`. `Z_k(A) = Z_Q + Z_L`;
  `Ik'' = 1.10·110/(√3·|Z_k|)`. Note `R/X(A) = 0.1` exactly (feeder and
  line share the ratio).
- **R3 — synchronous machine.** One 10.5 kV bus `G`; machine
  `x''_d = 0.15 pu` on 50 MVA / 10.5 kV. Network base:
  `x = 0.15·(100/50)·(10.5/10.5)² = 0.3 pu`; fictitious resistance
  `R_G = 0.07·x` (HV machine < 100 MVA, §6.6.3).
  `Z_k = |0.021 + j0.3|·(10.5²/100) Ω`; `Ik'' = 1.10·10.5/(√3·Z_k)`.
- **R4 — asynchronous motor.** One 10 kV bus `M`; motor
  `I_LR/I_rM = 5`, `S_rM = 5 MVA`, `U_rM = 10 kV`, locked-rotor
  `R/X = 0.1` (§6.7): `Z_M = (1/5)·10²/5 = 4 Ω` exactly;
  `Ik'' = 1.10·10/(√3·4)`.

## Reference table

| Case | Bus | Quantity | Analytic reference | Sparlectra | External tool |
|---|---|---|---|---|---|
| R1 | Q | Ik''max [kA] | 10.0 (exact) | 10.0 | *(open)* |
| R1 | Q | Sk''max [MVA] | 1905.2558883257648 | matches (rtol 1e-6) | *(open)* |
| R1 | Q | i_p [kA] | 28.284271247461902 | matches | *(open)* |
| R1 | Q | Ik''min [kA] | 8.0 (exact) | 8.0 | *(open)* |
| R2 | A | Z_k [Ω] | 67.78768576497586 | matches | *(open)* |
| R2 | A | Ik''max [kA] | 1.0305615508715191 | matches | *(open)* |
| R2 | A | i_p [kA] | 2.9148682442055054 | matches | *(open)* |
| R3 | G | Z_k [Ω] | 0.3315593472611653 | matches | *(open)* |
| R3 | G | Ik''max [kA] | 20.112223239140242 | matches | *(open)* |
| R3 | G | i_p [kA] | 56.88595774853494 | matches | *(open)* |
| R4 | M | Z_k [Ω] | 4.0 (exact) | 4.0 | *(open)* |
| R4 | M | Ik''max [kA] | 1.5877132402714709 | matches | *(open)* |
| R4 | M | i_p [kA] | 4.4907311951024935 | matches | *(open)* |

"matches" = verified by `exp_short_circuit_reference.jl` (rtol 1e-6); the
suite fails loudly if any quantity drifts.

Known modeling differences to keep in mind when filling the external
column: Sparlectra does not yet apply the transformer correction `K_T` or
the generator correction `K_G` (no transformers/`K_G`-relevant data in
these cases, so R1–R4 are unaffected), and `i_p` uses method b with the
1.15 factor and the 2.0 cap — some tools default to method c (equivalent
frequency), which yields slightly different peak values.

## External cross-check snippet (pandapower)

Run this OUTSIDE this repository in any Python environment with pandapower
(not installed or downloaded here by policy), then enter the resulting
numbers in the table above. Sketch — adjust to your pandapower version:

```python
import pandapower as pp
import pandapower.shortcircuit as sc

# R1/R2: feeder + line (110 kV)
net = pp.create_empty_network(sn_mva=100.0)
q = pp.create_bus(net, vn_kv=110.0, name="Q")
a = pp.create_bus(net, vn_kv=110.0, name="A")
# s_sc from the declared current: S = sqrt(3)*110 kV*10 kA = 1905.2559 MVA
pp.create_ext_grid(net, q, s_sc_max_mva=1905.2558883, rx_max=0.1,
                   s_sc_min_mva=1524.2047107, rx_min=0.1)
# 0.05+j0.5 pu on 121 ohm base -> 6.05 + j60.5 ohm (length 1 km)
pp.create_line_from_parameters(net, q, a, length_km=1.0,
                               r_ohm_per_km=6.05, x_ohm_per_km=60.5,
                               c_nf_per_km=0.0, max_i_ka=10.0)
sc.calc_sc(net, case="max", ip=True)   # compare ikss_ka/ip_ka at Q and A
sc.calc_sc(net, case="min")            # compare ikss_ka at Q

# R3: synchronous generator (10.5 kV), xdss 0.15 pu @ 50 MVA, rdss = 0.07*xdss
net3 = pp.create_empty_network(sn_mva=100.0)
g = pp.create_bus(net3, vn_kv=10.5)
pp.create_gen(net3, g, p_mw=0.0, vn_kv=10.5, sn_mva=50.0,
              xdss_pu=0.15, rdss_ohm=0.07*0.15*(10.5**2/50.0))
sc.calc_sc(net3, case="max", ip=True)

# R4: asynchronous motor (10 kV), I_LR/I_r = 5, 5 MVA, R/X = 0.1
net4 = pp.create_empty_network(sn_mva=100.0)
m = pp.create_bus(net4, vn_kv=10.0)
pp.create_motor(net4, m, pn_mech_mw=4.0, cos_phi=0.9, efficiency_n_percent=88.9,
                lrc_pu=5.0, rx=0.1, vn_kv=10.0, sn_mva=5.0)
sc.calc_sc(net4, case="max", ip=True)
```

Alternative references: the worked examples of IEC 60909-4 (or the German
Beiblatt 1 to DIN EN 60909-0) — a fixture slot for published example
numbers exists in `test/test_short_circuit.jl`.
