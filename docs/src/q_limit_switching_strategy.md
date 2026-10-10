# Q-limits in power flow: why switching strategy matters

A PV bus holds its voltage magnitude; the reactive power that takes is part
of the solution. When it leaves the generator's admissible range, the bus
becomes a PQ bus with reactive power fixed at the violated limit. On
networks with many PV buses and narrow or nearly identical Q ranges, how
and when that conversion happens decides whether the case converges or
enters a switching cascade. Two strategies:

- The **active-set strategy** switches during the Newton iteration. With
  many generators of narrow or nearly identical Q ranges, many PV-to-PQ and
  PQ-to-PV events can occur, and the algorithm solves a discrete switching
  problem on top of the smooth power-flow problem.
- The **classical strategy** first solves the power flow without switching,
  then clamps the violating generators to `Qmin` or `Qmax`, converts their
  buses to PQ and solves again, either all violations at once or the
  largest one per pass. A bus converted to PQ is not converted back within
  the same enforcement loop.

A CGMES `ReactiveCapabilityCurve` (Q limits as a function of active power)
is evaluated once at the machine's scheduled P during import and then
enforced by the switching strategy, not through the Q(U) path
([CGMES Import](cgmes_import.md)).

| Use | Where |
|---|---|
| Config key | `power_flow.qlimits.enforcement_mode`: `active_set` (default), `classic_simultaneous`, `classic_one_at_a_time`; the pass limit of the classic loop is `classic_max_passes` (default `30`); all keys in [Q-limit options and guard](@ref pf-qlimits) |
| Result field | run metadata `q_limit_classic_outer_loop_passes` and `q_limit_classic_outer_loop_stop` (`converged`, `max_outer_iterations`, `pf_not_converged_after_qlimit_update`, `no_reference_bus_remaining`, `base_pf_not_converged`); `run.log` names both on one line |
| Result field | solver status `converged_limits_failed` with reason `remaining_pv_q_limit_violations`, or with reason `max_outer_iterations` when the classic loop stopped at its pass limit (the run is then not converged, on the state of its last solve) |

- `active_set`: in-iteration switching with hysteresis, cooldown,
  narrow-range locking and repeated-switching protection; a clamped machine
  is released back to PV on the voltage side
  ([Powerlimits Guide](powerlimits.md)).
- `classic_simultaneous`: base power flow with switching disabled; if it
  converges, all detected violations are clamped and converted in one pass.
- `classic_one_at_a_time`: the same, but only the largest violation per
  pass, which makes the switching sequence easier to inspect.

A run can converge numerically and still fail to hold every reactive limit:
`converged_limits_failed` says that no admissible active set was reached,
and the run counts as unsuccessful although the bus balances are satisfied.
On `case300` the default `active_set` ends this way, while
`classic_simultaneous` reaches a limit-respecting solution with a few more
switching events. That status is the signal to try a classical mode.

Comparing modes is the diagnostic: a case that fails in `active_set` with a
repeatedly changing active set but behaves in a classical mode is dominated
by discrete switching; a base power flow that does not converge in a
classical mode has its problem upstream of Q-limit enforcement.
