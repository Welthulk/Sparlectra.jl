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

## Stability constraints of the final solution

Bus balances and the limits $Q_{min} \le Q \le Q_{max}$ do not make a
solution with Q limits a stable one. Zhao, Chiang, Li and Ju (PSCC 2008)
give two further conditions for every voltage-controlled generator bus:

```math
(U - U_{set})\,(Q_{max} - Q) \ge 0 \tag{6}
```

```math
(U - U_{set})\,(Q - Q_{min}) \le 0 \tag{7}
```

A bus with its voltage below the setpoint must be at $Q_{max}$, a bus with
its voltage above the setpoint at $Q_{min}$, and a bus inside its band must
hold the setpoint. This is the condition of a positive U-Q sensitivity.

The paper calls the switching rule that meets both conditions the original
logic (OL): a PV bus beyond a limit turns PQ at that limit; a PQ bus at
$Q_{max}$ returns to PV when $U > U_{set}$ and its Q is inside the band, one
at $Q_{min}$ when $U < U_{set}$. The voltage of a returning bus can be
damped, $U_{new} = U_{old} + \alpha\,(U_{set} - U_{old})$ with
$0 < \alpha \le 1$. Three supplementary logics stop oscillating bus types:

- SL1: a bus that switched from PV to PQ more than N times stays PQ; its
  return is no longer checked.
- SL2: as SL1, but per iteration only the bus with the largest violation is
  fixed.
- SL3: solve to convergence with fixed bus types, then switch by OL and
  solve again.

SL1 and SL2 can return a converged solution that violates (6) or (7) on
cases without a stable operating point, which the paper confirms with a
continuation power flow. Such a solution satisfies the equations but does
not describe a possible operating state.

| Mode | Return PQ to PV | Fixed buses | Damping of the return | Accepted with a type change in the last iteration | Class |
|---|---|---|---|---|---|
| `active_set` | yes, on the voltage side: at $Q_{max}$ when $U > U_{set}$ + `reenable_v_hyst_pu`, at $Q_{min}$ when $U < U_{set}$ - `reenable_v_hyst_pu`; a released bus that violates again is clamped again | a bus with two limit hits stays PQ, apart from one joint release of all such buses on a converged iterate, undone when one of them hits again; a narrow-range bus under `guard.narrow_range_mode: lock_pq` stays PQ from its first hit; buses in `lock_pv_to_pq_buses` run PQ from the start | no factor of its own: the release only sets the bus type back to PV, and the next Newton step moves the bus with the step scale of all buses (`damp`, autodamp or trust region). The bus takes the full Newton correction toward its setpoint ($\alpha = 1$) only with `damp = 1` and neither autodamp nor trust region, the defaults of the configuration path and of `runpf_rectangular!(net, maxIte, ...)`; the keyword method `runpf_rectangular!(net; ...)` defaults to `damp = 0.2` | no, an iterate in which a bus type changed is not accepted | OL with an SL1 lock after the second hit (N = 1, N = 2 for the one joint release) |
| `classic_simultaneous` | no | every clamped generator stays at its limit, and its bus turns PQ once no regulating generator is left on it; `lock_pv_to_pq_buses` are converted before the base solve | no return | no, every pass ends with a converged solve at fixed types, and the loop stops only when that solve leaves no violation | one-sided SL3: PV to PQ only, no return check |
| `classic_one_at_a_time` | no | as `classic_simultaneous` | no return | as `classic_simultaneous` | one-sided SL3, one generator per pass: the one with the largest violation in the network |

Notes on the table:

- The return in `active_set` needs `hysteresis_pu > 0` or
  `cooldown_iters > 0` (defaults `0.01` and `1`); with both at zero it is
  one-sided like the classic modes.
- The PV to PQ switch of `active_set` uses one of two margins. Under
  `guard.violation_mode: lock_pq` (the default) a machine switches once its
  Q exceeds a limit by `guard.violation_threshold_pu` (default `1e-4` pu),
  so `hysteresis_pu` (default `0.01` pu) does not delay the switch; it is
  the switching margin only under `delayed_switch`, or when it is set below
  the threshold. Its other uses stay: it enables the return (see above), it
  is the within-hysteresis band of the final Q check, and the classic modes
  take `max(tol, hysteresis_pu)` as their violation tolerance.
- `guard.max_switches` (default `3`) with
  `guard.freeze_after_repeated_switching` stops switching a bus that
  reached that many hits, and the run is then reported as not converged
  (`max_switching_exceeded`). It does not produce a converged solution with
  fixed buses.
- A zero-range machine ($Q_{min} = Q_{max}$) runs as PQ at its only value
  under the guard and meets (6) and (7) by construction. A bus run as PQ
  inside its band (`guard.narrow_range_mode: prefer_pq`, or a bus in
  `lock_pv_to_pq_buses` whose scheduled Q lies inside its limits) meets
  them only when its voltage equals the setpoint.
- A fixed bus stays PQ to the end. The run status does not check it
  against (6) or (7). The Q-V characteristic check
  (`qvCharacteristicViolations`, printed with the results) lists the buses
  at a limit on the wrong side of their setpoint, read from the stability
  check below; it does not change the status. The release margin `reenable_v_hyst_pu` (default `1e-4` pu)
  leaves a clamped bus up to that much on the wrong side without a
  release.
- After every converged run `run_sparlectra` checks (6) and (7) at every
  voltage-regulated bus, fixed or not, and a PV bus against its band
  (`diagnostics.qlimit_stability_report`, one warning per run with a
  violation, see [Local PowerFlow Service](powerflow_service.md)). The
  classic modes and a fixed bus can end on a state that fails it; the run
  status stays as it is.
- The classic modes are one-sided on purpose: they follow MATPOWER, which
  converts a violating generator and does not check a return. A converged
  classic run can therefore end at $Q_{max}$ with the voltage above the
  setpoint, or at $Q_{min}$ below it; the stability check reports such end
  states. A bus held PQ by a lock (the two-hit lock or the narrow-range
  lock of `active_set`, `lock_pv_to_pq_buses` in every mode) that fails (6)
  or (7) is named in the warning with the lock that holds it (`locked_by`
  in the report).

## Divergence at the iteration limit

Zhao, Chiang, Li and Ju (PSCC 2008) separate two ways a Newton power flow
with Q limits fails. In a numerical divergence the mismatch does not come
down. In an identification divergence the mismatch is small, but the type of
some regulated buses keeps changing between PV and PQ; the paper reads it as
a sign that no stable operating point exists, the counterpart of a
limit-induced bifurcation.

A rectangular run that uses up `power_flow.max_iter` names the case in its
`reason`, read from its last three iterations; the status stays
`not_converged`:

- `type_identification_divergence`: a bus type changed, and the smallest
  mismatch measured before the switch is at most 10 times
  `power_flow.tol`;
- `iteration_limit`: no type change, the mismatch fell from each
  iteration to the next, and the last mismatch is below the first mismatch
  of the run; the run was still converging;
- `numerical_divergence`: everything else, also a run that climbed away
  from its start and falls again at the end.

The classic modes switch between passes, at a converged state; a pass that
reaches `max_iter` has no type change of its own and ends as
`iteration_limit` or `numerical_divergence`. A run with a type change in its
last iteration is never accepted as converged (table above), so the
identification divergence appears only as a stop at `max_iter`.

On `case14_zeng_p306_activeSet_A` with Q-limit switching held back to a
converged iterate (`qlimits.start_mode: iteration`, `start_iter: 99`),
`max_iter: 3` stops at a mismatch of 7.5e-3 that is still falling
(`iteration_limit`); with `max_iter: 5` the run converges in iteration 5,
four machines switch to PQ, and the run stops there
(`type_identification_divergence`). The last iterate is checked against
(6) and (7) as well (`qlimit_stability_status = last_iterate`).
