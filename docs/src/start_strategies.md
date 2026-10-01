# [Start strategies by case](@id start_strategies_page)

A flat start sets every bus to 1 pu and 0 degrees, PV and slack buses to
their setpoint magnitudes. The file's VM and VA columns are the state the
case author stored; on some files a solved state of the case, on others a
rough profile that satisfies nothing. Whether plain Newton converges from
a flat start depends on the distance between that start and the solution
of the file, and on how the Newton step is applied
([Newton update](newton_update.md)); it is not a property of the solver
alone, and it is not a property of the case alone.

## The ladder

The start strategies, from the cheapest to the most expensive, each with
the key that switches it on (every key under `power_flow`):

1. **Flat start, plain Newton**: `flatstart: true`, `autodamp: false`.
2. **More iterations**: `max_iter: 80` instead of 30; it decides nothing
   on the cases below, a flat start that diverges in 16 steps diverges in
   80 too.
3. **Adaptive damping**: `autodamp: true`; the step is shortened whenever
   the full step does not reduce the mismatch.
4. **DC-seed angles**: `flatstart: false`, `start_mode.angle_mode: dc`,
   `start_mode.dc_seed_unconditional: true`; flat magnitudes with the
   angles of a DC power flow.
5. **Projection with blend scan**: `flatstart: false`,
   `start_mode.start_projection: true`, `try_dc_start: true`,
   `try_blend_scan: true`, `voltage_mode: profile_blend`; the start that
   has the smallest mismatch among the raw state, the DC candidate and the
   blends is taken.
6. **Current-iteration pre-solve**: `start_current_iteration.enabled:
   true`; a few guarded current-injection iterations before Newton.
7. **APSLF start values**: `apslf_start.enabled: true`; the analytic
   power-series solution as the Newton start.
8. **The file's columns**: `flatstart: false` with no start machine; the
   stored VM and VA as the start.
9. **Rescue**: `rescue: true`; after a non-converged solve, the ladder
   alternate start, damping, DC seed, settled Q limits is tried in turn.
10. **Auto mode**: `mode: auto`; the configuration is filled from the
    imported network, the stages above are escalated as needed.

`flatstart: true` forces every start machine off; a flat start with a
machine is written `flatstart: false` plus the machine, which then
overrides the stored columns.

## Recipes

Plain Newton from the flat start with the polar update and no damping,
30 iterations, and what to set where that is not enough. Iterations are
mismatch evaluations (Newton steps plus one). The APSLF column is the
analytic power-series solver (`solver: apslf`, order 24 with the Newton
polish).

| case | buses | plain Newton from the flat start | what to set when not | iterations | APSLF solver | remark |
|---|---:|---|---|---:|---|---|
| case14 | 14 | converges | | 5 | converges | |
| case118 | 118 | converges | | 5 | converges | |
| case300 | 300 | converges | | 6 | converges | |
| case1354pegase | 1354 | converges | | 6 | converges | |
| case2848rte | 2848 | converges to a low-voltage branch (min Vm 0.02 pu, 893.6 MW losses) | `flatstart: false` (the file's columns), or the DC seed, for the operating solution (min Vm 0.89 pu, 607.4 MW) | 10 | converges, operating solution | the wrong-branch check flags the low-voltage branch (`voltage_collapse`) |
| case2869pegase | 2869 | converges | | 6 | converges | |
| case3120sp | 3120 | converges | | 7 | converges | |
| case9241pegase | 9241 | converges | | 7 | converges | diverged with the rectangular update; needed `autodamp: true` (13) up to 0.20.5 |
| case1888rte | 1888 | does not converge | DC-seed angles (recipe 4) | 6 | converges, operating solution | with the rectangular update the flat start converged to a low-voltage branch (min Vm 0.06 pu, 1396.9 MW against 980.7 MW) |
| case6495rte | 6495 | does not converge | DC-seed angles (recipe 4) | 8 | does not converge: phase shifters of 6 to 10 degrees on x = 3e-4 pu; not yet supported in the APSLF path | 2543.8 MW; the file's columns converge in 3 |
| case4_dist | 4 | converges | | 4 | converges | |
| case18 | 18 | converges | | 5 | converges | |
| case33bw | 33 | converges | | 4 | converges | |
| mvlv1004 | 1004 | converges | | 7 | converges | |
| mvlv10616 | 10616 | converges | | 6 | converges | |
| mvlv29840 | 29840 | converges | | 6 | converges | |
| case13659pegase | 13659 | does not converge | projection with blend scan (recipe 5) from the flat start, or the file's columns | 9 | does not converge: the reference bus (a 42 MW machine behind one 0.14 pu transformer) carries the file's dispatch surplus along the series path | 8737.2 MW, min Vm 0.84 pu; the undamped rungs with the DC seed land on a second root with 164 degrees across the slack transformer, which the wrong-branch check rejects (`bus_angle_exceeded`) |

A converged state is not always the operating state. The wrong-branch
check (`power_flow.wrong_branch_detection`, see
[Power-Flow Configuration](powerflow_configuration.md)) judges every
voltage level of the result and names the reason; a low-voltage branch
satisfies the power-flow equations to the tolerance like the operating
solution does, so only the voltages tell them apart.

## Where the cases come from

Sparlectra does not ship any of the cases above; the table names them so
a reader can reproduce the rows. Obtain them from their sources under
their own license, citation and redistribution terms:

- `case14`, `case118`, `case300`, `case1354pegase`, `case2869pegase`,
  `case9241pegase`, `case13659pegase`, `case1888rte`, `case2848rte`,
  `case6495rte`, `case3120sp`, `case4_dist`, `case18`, `case33bw`: MATPOWER
  case files, `data/<name>.m` in the MATPOWER repository
  (<https://github.com/MATPOWER/matpower>). Cite MATPOWER as its guidance
  asks, see [Citation and case-file usage](@ref matpower-citation). The
  PEGASE and RTE files request additional case-specific citations in
  their headers (Fliscounakis et al. 2013 and Josz et al. 2016 for the
  PEGASE cases, Josz et al. 2016 for the RTE snapshots); `case3120sp` is
  the Polish system model distributed with MATPOWER; the small
  distribution feeders carry the references of their original papers in
  their headers.
- `mvlv1004`, `mvlv10616`, `mvlv29840`: synthetic radial MV/LV
  distribution grids, generated, not measured. The generator is the port
  of power-grid-model's `FictionalGridGenerator` published in gridoxide
  0.0.2 (Apache-2.0); the benchmark-grids repository
  (<https://github.com/m-mirz/benchmark-grids>, directory `generated/`,
  script `scripts/generate_distribution.py`, fixed seed, byte-identical
  on every run) generates them and exports them to MATPOWER format. Their
  loading is power-grid-model's benchmark loading (voltages 0.67 to 0.89
  pu), not an operating point.

## Conventions

The standard MATPOWER reading (shift in degrees, sign +1, ratio as
stored, bus shunt as admittance) satisfies the power balance of the file
on every case above; the alternative readings (radians, sign -1,
reciprocal ratio) fail wherever a file has a phase shift or an
off-nominal tap. The stored VM and VA columns of the PEGASE files are not
a solved state under any reading (best fit 0.5 to 3.9 pu at the worst
bus; a solved state scores below 0.1 pu), which is why a scan of those
columns must not be used to change the reading. The convention overrides
and the scan are an experimental feature, off by default; see
[MATPOWER conventions (experimental)](@ref matpower_conventions_experimental).

## Q limits

With reactive limits on, the enforcement method selects the solution. On
case300 from the flat start the in-iteration active set converts 48 PV
buses and ends at 422.7 MW of losses, the classic simultaneous and the
classic one-at-a-time methods convert 28 and end at 395.4 and 457.9 MW;
without limits the case solves at 408.3 MW. All are solutions of the
equations with their limit sets. A comparison with limits on therefore
needs the method named (`power_flow.qlimits.enforcement_mode`, see
[Q-limit Switching Strategy](q_limit_switching_strategy.md)).
