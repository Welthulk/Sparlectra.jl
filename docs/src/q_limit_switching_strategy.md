# Q-limits in power flow: why switching strategy matters

Reactive-power limits of generators look like a minor modelling detail.
In practice they decide whether a network case converges robustly or
enters a switching cascade.

A CGMES `ReactiveCapabilityCurve` (Q limits as a function of active
power) is evaluated once at the machine's scheduled P during import and
then enforced by the switching strategy, not through the
voltage-dependent Q(U) path; see [CGMES Import](cgmes_import.md).

A PV bus keeps its voltage magnitude fixed; the reactive power needed to
hold it is determined by the solution. If that reactive power exceeds the
generator's admissible range, the bus must become a PQ bus with reactive
power fixed at the violated limit. Two strategies for this delicate
conversion:

- The **active-set strategy** reacts to violations during the nonlinear
  iteration itself. Flexible, but with many generators of narrow or
  nearly identical Q ranges, many PV-to-PQ and PQ-to-PV events can occur,
  and the algorithm is then solving a discrete switching problem on top
  of the smooth power-flow problem.
- The **classical strategy** first solves the base power flow without
  switching. Only after it has converged are violations checked;
  violating generators are clamped to `Qmin` or `Qmax`, the buses
  converted to PQ, and the power flow solved again, either for all
  violating generators at once or step by step for the largest violation.
  A bus converted to PQ is not converted back within the same enforcement
  loop.

Large synthetic or aggregated networks expose the difference: many PV
buses, narrow Q bands and mixed start values make it hard to distinguish
voltage control from reactive-power limitation, and the switching
strategy, not the model, may be the decisive factor.

| Use | Where |
|---|---|
| Config key | `power_flow.qlimits.enforcement_mode`: `active_set` (default), `classic_simultaneous`, `classic_one_at_a_time` |
| Config key | `power_flow.qlimits.reenable_v_hyst_pu` (default `1e-4` pu): voltage hysteresis for releasing a clamped machine back to PV |
| Result field | solver status `converged_limits_failed` with reason `remaining_pv_q_limit_violations` |

- `active_set`: in-iteration switching with guards such as hysteresis,
  cooldown, narrow-range locking and repeated-switching protection.
- `classic_simultaneous`: base power flow with switching disabled; if it
  converges, all detected violations are clamped and converted in one
  outer-loop pass.
- `classic_one_at_a_time`: same principle, but only the largest violation
  per outer-loop pass, which makes the switching sequence easier to
  inspect.

The active-set mode also releases a clamped machine back to PV, decided on
the voltage side: a machine at Qmax whose voltage sits above its setpoint
by more than the hysteresis needs less than Qmax to hold the setpoint and
goes back to PV; the mirror image applies at Qmin. Cooldown and the
one-retry guard still apply, so a machine flipping between the clamp and
the voltage constraint is held after its first retry.

A run can converge numerically and still fail to hold every reactive
limit: the status `converged_limits_failed` says that no admissible active
set was reached, and the run counts as unsuccessful although the bus
balances are satisfied. Example: on `case300` the default `active_set`
ends this way, while `classic_simultaneous` reaches a limit-respecting
solution with a few more switching events. That status is the signal to
try a classical mode.

Comparing modes is the practical diagnostic: a case that fails in
`active_set` with a repeatedly changing active set but behaves in a
classical mode is dominated by discrete switching; a base power flow that
does not converge in a classical mode has its problem upstream of Q-limit
enforcement. Q-limits change the structure of the power-flow problem, so
large-network analysis needs a controlled strategy for the discrete
switching between PV and PQ.
