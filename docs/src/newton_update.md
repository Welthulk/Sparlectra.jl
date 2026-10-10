# [Newton update: rectangular versus polar](@id newton_update_page)

Every Newton step solves `J dx = -f` for a correction and applies it.
Sparlectra's rectangular solver and a polar solver such as MATPOWER's
`newtonpf` solve the same linear system (one Jacobian in different
variables), so every bus receives the same correction `dV`; they differ in
how `dV` is applied, and on a large transmission case that decides whether
a flat start converges.

## The two updates

Relative to the current voltage, `dV / V = a + jb` asks for a magnitude
change by the fraction `a` and a turn by the angle `b`.

| update | new voltage | magnitude after the step |
|---|---|---|
| rectangular | `V_new = V (1 + a + jb)` | `\|V\| sqrt((1 + a)^2 + b^2)` |
| polar | `V_new = \|V\| (1 + a) e^{j(theta + b)}` | `\|V\| (1 + a)` |

The polar update applies exactly that. The rectangular update adds the
correction as a vector, and a turn by `b` inflates the magnitude by
`sqrt(1 + b^2)` although none was asked for. For small corrections the two
agree, so both converge to the same state; for a large turn they do not: a
bus at 1.00 pu that has to turn by 60 degrees (1.05 rad) ends at 1.00 pu
and 60 degrees with the polar update, and at `1 + 1.05j`, 1.45 pu at 46
degrees, with the rectangular one.

## Why the flat start is where it matters

A flat start sets every angle to 0, so the first step is close to a DC
power flow for the angles and turns whole regions by the solution angles at
once (case9241pegase: -61 to +70 degrees). After that one step:

| case | rectangular update | polar update |
|---|---|---|
| case9241pegase | 9226 of 9241 buses above 1.1 pu, the worst at 7 pu, PV buses up to 5.9 pu off their setpoints | every bus between 0.91 and 1.21 pu, PV buses on their setpoints |
| case3120sp | 3088 buses above 1.1 pu, the worst at 2.25 pu | 0.90 to 1.13 pu |
| case2869pegase | 2444 buses above 1.1 pu, the worst at 2.37 pu | 0.97 to 1.15 pu |

From thousands of buses at 2 to 7 pu Newton does not come back on
case9241pegase (mismatch 533, 10900, 1.4e6, 9.6e10, non-finite after 16
steps); on case3120sp and case2869pegase it comes back at the price of two
to three extra steps. From the file's stored voltages the corrections are
small and both updates behave alike.

## What changes for the user

`power_flow.newton_update` selects `polar` (default) or `rectangular` (the
update of releases up to 0.20.5, kept for comparisons). In the Web UI the
select sits in the Experimental block at the end of the Advanced options on
the Settings page, greyed until "Enable experimental settings" is ticked;
while greyed the configured default applies. With `polar` the rectangular
Jacobian reproduces MATPOWER's `newtonpf` iterates step for step.

- From the file's stored voltages both updates reach the same state in the
  same number of steps.
- From a flat start `polar` converges where `rectangular` diverged
  (case9241pegase in 7 iterations without damping) and saves one to three
  iterations on transmission cases; on distribution feeders, whose angles
  are a few degrees, it costs at most one.
- Damping (`autodamp`, the merit line search, the trust region) scales the
  correction before it is applied, so every start machine and step control
  works with either update
  ([Start strategies by case](start_strategies.md)).

## The active set and the path

A PV bus is switched to PQ on the iterates along the way, so a different
update can give a different switching sequence and, with binding limits, a
different limited solution. The default `power_flow.qlimits.start_iter: 3`
with `start_mode: iteration_or_auto` lets the iterate settle before the
first switch and protects both updates; switching from the second iterate
is not safe with either
([Q-limit Switching Strategy](q_limit_switching_strategy.md)).
