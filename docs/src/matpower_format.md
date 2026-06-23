# MATPOWER format

Sparlectra can read power-flow test cases in the MATPOWER case format. A MATPOWER case is usually a `.m` file that defines or returns a MATLAB/Octave structure named `mpc`. In the current case format, this structure contains `baseMVA`, `bus`, `gen`, `branch` and optionally `gencost`.[^matpower-data-format]

For a normal power-flow run in Sparlectra, the most relevant parts are `baseMVA`, `bus`, `gen` and `branch`. The optional `gencost` table is mainly relevant for optimal power flow and is usually not needed for a plain load-flow calculation.[^matpower-caseformat]

## Minimal structure

A typical case file has this shape:

```matlab
function mpc = case_example
mpc.version = '2';
mpc.baseMVA = 100;

mpc.bus = [
    % bus_i type Pd Qd Gs Bs area Vm Va baseKV zone Vmax Vmin
      1     3   0  0  0  0  1    1  0  110    1    1.1  0.9;
      2     1  50 30  0  0  1    1  0  110    1    1.1  0.9;
];

mpc.gen = [
    % bus Pg Qg Qmax Qmin Vg mBase status Pmax Pmin
      1   80 0  100  -100 1.0 100   1      100  0;
];

mpc.branch = [
    % fbus tbus r     x      b      rateA rateB rateC ratio angle status angmin angmax
      1    2    0.02  0.06   0.03   100   100   100   0     0     1      -360   360;
];
```

The exact number and meaning of columns is fixed by MATPOWER. The official MATPOWER documentation is the authoritative reference for the complete column lists.[^matpower-data-format]

## `baseMVA`

`baseMVA` is the common power base of the case. MATPOWER uses it for per-unit quantities. Typical examples use values such as `100`, but the correct value depends on the case data.

## `bus`: network nodes

`mpc.bus` describes the electrical buses. Each row is one bus. The most important columns for users are:

| Column | Meaning |
|---|---|
| `BUS_I` | Bus number. Branches and generators refer to this number. |
| `BUS_TYPE` | Bus type: `1` = PQ, `2` = PV, `3` = reference/slack, `4` = isolated. |
| `PD`, `QD` | Active and reactive demand in MW and MVAr. |
| `VM`, `VA` | Voltage magnitude in p.u. and voltage angle in degrees. Often used as initial or reference values. |
| `BASE_KV` | Nominal voltage level in kV. |
| `VMAX`, `VMIN` | Upper and lower voltage magnitude limits in p.u. |

A PQ bus normally represents a load bus. A PV bus usually has generation with a specified active power and voltage magnitude. The reference bus supplies the angular reference for the island.

## `gen`: generators and voltage setpoints

`mpc.gen` describes generators or equivalent injections. Each row is one generator connected to a bus.

Important columns are:

| Column | Meaning |
|---|---|
| `GEN_BUS` | Bus number where the generator is connected. |
| `PG`, `QG` | Active and reactive generator output in MW and MVAr. |
| `QMAX`, `QMIN` | Reactive power limits. These are relevant for Q-limit handling. |
| `VG` | Voltage magnitude setpoint in p.u. |
| `GEN_STATUS` | `> 0` means in service, `<= 0` means out of service. |
| `PMAX`, `PMIN` | Active power limits in MW. |

For load-flow studies, `PG`, `VG`, `QMAX` and `QMIN` are usually the values users inspect first.

## `branch`: lines and transformers

`mpc.branch` describes network connections. Each row is one line, cable, transformer or equivalent branch.

Important columns are:

| Column | Meaning |
|---|---|
| `F_BUS`, `T_BUS` | From-bus and to-bus numbers. |
| `BR_R`, `BR_X` | Series resistance and reactance in p.u. |
| `BR_B` | Total line charging susceptance in p.u. |
| `RATE_A`, `RATE_B`, `RATE_C` | Thermal ratings. `0` is commonly used for unlimited. |
| `TAP` | Transformer off-nominal tap ratio. `0` indicates no off-nominal transformer tap and is mathematically treated like a ratio of `1`. |
| `SHIFT` | Transformer phase-shift angle in degrees. |
| `BR_STATUS` | `1` means in service, `0` means out of service. |
| `ANGMIN`, `ANGMAX` | Optional voltage angle difference limits. |

In MATPOWER, lines and transformers are represented in the same `branch` matrix. A non-zero tap ratio or phase shift indicates transformer behavior.

## `gencost`: optional cost data

`mpc.gencost` contains generator cost data for optimal power-flow studies. It can be present in many public MATPOWER cases, but it is usually not needed for a plain power-flow run.

## Practical checks before importing a case

Before running a case in Sparlectra, check the following points:

1. The file should define or return `mpc`.
2. `mpc.version` should normally be `'2'`.
3. Bus numbers in `branch` and `gen` must exist in `mpc.bus`.
4. At least one in-service reference bus is required per connected island.
5. Decimal values must use a decimal point, not a decimal comma.
6. Inactive generators or branches should have status `0`.
7. Per-unit branch impedances must be consistent with `baseMVA` and the voltage base of the case.

## What Sparlectra uses

For a standard Web UI power-flow run, Sparlectra primarily uses the topology and electrical data from `bus`, `gen` and `branch`. The result page and log files then report whether the numerical power flow converged, how many iterations were needed and what final mismatch remained.

When Q-limit handling is enabled, the generator reactive power limits from `QMAX` and `QMIN` become important. If a generator cannot maintain its voltage setpoint within those limits, the solver may have to treat the bus differently depending on the selected configuration.

## Full reference

This page is intentionally short. For the complete MATPOWER format, use the official MATPOWER documentation:

- MATPOWER Manual, Appendix “Data File Format”: <https://matpower.app/manual/matpower/DataFileFormat.html>
- MATPOWER `caseformat` reference: <https://matpower.org/documentation/ref-manual/legacy/functions/caseformat.html>
- MATPOWER named column constants: <https://matpower.org/docs/ref/matpower6.0/define_constants.html>

[^matpower-data-format]: MATPOWER Manual, “Data File Format”. It describes version 2 case files as an `mpc` structure with `baseMVA`, `bus`, `branch`, `gen` and optional `gencost`; it also defines the bus, generator, branch and cost columns. <https://matpower.app/manual/matpower/DataFileFormat.html>
[^matpower-caseformat]: MATPOWER `caseformat` reference. It defines a MATPOWER case file as an M-file or MAT-file that returns or defines an `mpc` case structure and documents version 2, optional fields and solver-added result columns. <https://matpower.org/documentation/ref-manual/legacy/functions/caseformat.html>
