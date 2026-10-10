# [Start strategies by case](@id start_strategies_page)

A flat start sets every bus to 1 pu and 0 degrees, PV and slack buses to
their setpoint magnitudes. The file's VM and VA columns are what the case
author stored: a solved state on some files, a rough profile on others.
Whether plain Newton converges from a flat start depends on the distance
between start and solution and on how the step is applied
([Newton update](newton_update.md)), not on the solver or the case alone.

## The ladder

From the cheapest to the most expensive, each with its key under
`power_flow`:

1. **Flat start, plain Newton**: `flatstart: true`, `autodamp: false`,
   `start_mode.start_projection: false`.
2. **More iterations**: `max_iter: 80` instead of 30; decides nothing
   below, a flat start that diverges in 16 steps diverges in 80 too.
3. **Adaptive damping**: `autodamp: true`; the step is shortened whenever
   the full step does not reduce the mismatch.
4. **DC-seed angles**: `flatstart: false`, `start_mode.angle_mode: dc`,
   `start_mode.dc_seed_unconditional: true`; flat magnitudes with the
   angles of a DC power flow.
5. **Projection with blend scan**: `flatstart: false`,
   `start_mode.start_projection: true`, `try_dc_start: true`,
   `try_blend_scan: true`, `voltage_mode: profile_blend`; the start with
   the smallest mismatch among the raw state, the DC candidate and the
   blends.
6. **Current-iteration pre-solve**: `start_current_iteration.enabled:
   true`; a few guarded current-injection iterations before Newton.
7. **APSLF start values**: `apslf_start.enabled: true`; the analytic
   power-series solution as the Newton start.
8. **The file's columns**: `flatstart: false` with no start machine.
9. **Rescue**: `rescue: true`; after a non-converged solve, alternate
   start, damping, DC seed and settled Q limits are tried in turn.
10. **Auto mode**: `mode: auto`; the configuration is filled from the
    imported network and the stages above are escalated as needed.

`flatstart: true` sets only the start profile: every start machine follows
its own key and runs from the flat profile. The two start modes
(`angle_mode`, `voltage_mode`) choose a profile themselves and count as
`classic` under a flat start; the DC angles stay available through the
projection's DC candidate. Recipe 1 therefore needs the projection off as
well.

With the projection on, a flat start also tries the ratio profile
(`start_mode.ratio_profile`, default `true`): the flat magnitudes scaled per
voltage level by the off-nominal transformer ratios on the path from the
reference. A flat start drives a circulating current through every
off-nominal transformer, and a stiff winding off nominal (the star
equivalent of a three-winding transformer) can dominate the start: on the
CGMES MiniGrid it is 164 pu of active power on a 9 MW case and the polar
update diverges; the ratio profile starts at 0.09 pu and converges in three
steps to the delivery's own state.

Independent of the projection, every flat start sets the auxiliary buses
(the star points of three-winding transformers) to the current-free value
of their lowest-impedance winding, the visible buses keep 1.0 pu and 0
degrees; that alone lets the MiniGrid converge from the bare flat start
(polar update, projection off) in four steps. A network without auxiliary
buses starts as before.

## Recipes

Plain Newton from the flat start with the polar update, no damping, 30
iterations, and what to set where that is not enough. Iterations are
mismatch evaluations (Newton steps plus one). The APSLF column is
`solver: apslf`, order 24 with the Newton polish.

| case | buses | plain Newton from the flat start | what to set when not | iterations | APSLF solver | remark |
|---|---:|---|---|---:|---|---|
| case14 | 14 | converges | | 5 | converges | |
| case118 | 118 | converges | | 5 | converges | |
| case300 | 300 | converges | | 6 | converges | |
| case1354pegase | 1354 | converges | | 6 | converges | |
| case2848rte | 2848 | converges to a low-voltage branch (min Vm 0.02 pu, 893.6 MW losses) | `flatstart: false` (the file's columns), or the DC seed, for the operating solution (min Vm 0.89 pu, 607.4 MW) | 10 | converges, operating solution | the wrong-branch check flags the low-voltage branch (`voltage_collapse`) |
| case2869pegase | 2869 | converges | | 6 | converges | |
| case3120sp | 3120 | converges | | 7 | converges | |
| case9241pegase | 9241 | converges | | 7 | converges | with the rectangular update it diverged and needed `autodamp: true` (13 iterations) |
| case1888rte | 1888 | does not converge | DC-seed angles (recipe 4) | 6 | converges, operating solution | with the rectangular update the flat start converged to a low-voltage branch (min Vm 0.06 pu, 1396.9 MW against 980.7 MW) |
| case6495rte | 6495 | does not converge | DC-seed angles (recipe 4) | 8 | does not converge, with or without its phase shifters (all 17 set to 0 degrees: still not converged) | 2543.8 MW; the file's columns converge in 3 |
| case4_dist | 4 | converges | | 4 | converges | |
| case18 | 18 | converges | | 5 | converges | |
| case33bw | 33 | converges | | 4 | converges | |
| mvlv1004 | 1004 | converges | | 7 | converges | |
| mvlv10616 | 10616 | converges | | 6 | converges | |
| mvlv29840 | 29840 | converges | | 6 | converges | |
| case13659pegase | 13659 | does not converge | projection with blend scan (recipe 5) from the flat start, or the file's columns | 9 | does not converge: the reference bus (a 42 MW machine behind one 0.14 pu transformer) carries the file's dispatch surplus along the series path | 8737.2 MW, min Vm 0.84 pu; the undamped rungs with the DC seed land on a second root rotated by 164 degrees against the reference bus (170 degrees across the reference transformer against 6 on the operating root), which the wrong-branch check flags (`reference_branch_angle_exceeded`) |

A converged state is not always the operating state: a low-voltage branch
satisfies the equations to the tolerance like the operating solution, only
the voltages tell them apart. The wrong-branch check
(`power_flow.wrong_branch_detection`,
[Power-Flow Configuration](powerflow_configuration.md)) judges every
voltage level of the result and names the reason.

## Where the cases come from

Sparlectra ships none of the cases above. Obtain them from their sources
under their own license and citation terms:

- `case14`, `case118`, `case300`, `case1354pegase`, `case2869pegase`,
  `case9241pegase`, `case13659pegase`, `case1888rte`, `case2848rte`,
  `case6495rte`, `case3120sp`, `case4_dist`, `case18`, `case33bw`: MATPOWER
  case files, `data/<name>.m` in <https://github.com/MATPOWER/matpower>;
  cite MATPOWER ([Citation and case-file usage](@ref matpower-citation)).
  The PEGASE and RTE files ask for case-specific citations in their headers
  (Fliscounakis et al. 2013 and Josz et al. 2016 for PEGASE, Josz et al.
  2016 for RTE); `case3120sp` is the Polish system model distributed with
  MATPOWER; the distribution feeders name their papers in their headers.
- `mvlv1004`, `mvlv10616`, `mvlv29840`: synthetic radial MV/LV distribution
  grids, generated with a fixed seed by `scripts/generate_distribution.py`
  of <https://github.com/m-mirz/benchmark-grids> (directory `generated/`;
  the generator is the port of power-grid-model's `FictionalGridGenerator`
  in gridoxide 0.0.2, Apache-2.0) and exported to MATPOWER format. Their
  loading is power-grid-model's benchmark loading (0.67 to 0.89 pu), not an
  operating point.

## Conventions

The standard MATPOWER reading (shift in degrees, sign +1, ratio as stored,
bus shunt as admittance) satisfies the power balance of the file on every
case above; the alternative readings (radians, sign -1, reciprocal ratio)
fail wherever a file has a phase shift or an off-nominal tap. The PEGASE
VM and VA columns are not a solved state under any reading (best fit 0.5 to
3.9 pu at the worst bus; a solved state scores below 0.1 pu), so a scan of
those columns must not change the reading. Overrides and scan are
experimental and off by default:
[MATPOWER conventions (experimental)](@ref matpower_conventions_experimental).

## Q limits

With reactive limits on, the enforcement method selects the solution: on
case300 from the flat start the in-iteration active set converts 48 PV
buses and ends at 422.7 MW of losses, the classic simultaneous and
one-at-a-time methods convert 28 and end at 395.4 and 457.9 MW, without
limits the case solves at 408.3 MW. All satisfy the equations with their
limit sets, so a comparison with limits on names the method
(`power_flow.qlimits.enforcement_mode`,
[Q-limit Switching Strategy](q_limit_switching_strategy.md)).
