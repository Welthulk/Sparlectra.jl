# [Newton update: rectangular versus polar](@id newton_update_page)

Newton's method solves, in every step, the linear system `J dx = -f` for a
correction of the state and then applies that correction. The
rectangular solver of Sparlectra and a polar solver such as MATPOWER's
`newtonpf` linearise the same power-balance equations and solve the same
linear system (the two Jacobians are one matrix written in different
variables: angle and magnitude against real and imaginary part), so the
correction `dV` every bus receives is the same. What differs is how `dV`
is applied to the voltage, and that choice decides whether a flat start
converges on a large transmission case.

## The two updates

Written relative to the current voltage, the correction `dV / V = a + jb`
says two things: change the magnitude by the fraction `a` and turn the
phasor by the angle `b`.

| update | new voltage | magnitude after the step |
|---|---|---|
| rectangular | `V_new = V (1 + a + jb)` | `\|V\| sqrt((1 + a)^2 + b^2)` |
| polar | `V_new = \|V\| (1 + a) e^{j(theta + b)}` | `\|V\| (1 + a)` |

The polar update applies exactly what the correction asked for. The
rectangular update adds the correction as a vector, and a turn by `b`
inflates the magnitude by `sqrt(1 + b^2)` although no magnitude change
was asked for. For small corrections the two agree, which is why both
converge to the same state at the end; for a large turn they do not.

The example: a bus at 1.00 pu that has to turn by 60 degrees (1.05 rad)
with no magnitude change. The polar update gives 1.00 pu at 60 degrees.
The rectangular update gives `1 + 1.05j`, a phasor of 1.45 pu at 46
degrees: the magnitude grew by 45 percent and the angle fell short.

## Why the flat start is where it matters

A flat start sets every angle to 0. The first Newton step is then close to
a DC power flow for the angles, and on a transmission case it has to turn
whole regions by the solution angles at once (case9241pegase: -61 to +70
degrees). After that one step:

| case | rectangular update | polar update |
|---|---|---|
| case9241pegase | 9226 of 9241 buses above 1.1 pu, the worst at 7 pu, PV buses up to 5.9 pu off their setpoints | every bus between 0.91 and 1.21 pu, PV buses on their setpoints |
| case3120sp | 3088 buses above 1.1 pu, the worst at 2.25 pu | 0.90 to 1.13 pu |
| case2869pegase | 2444 buses above 1.1 pu, the worst at 2.37 pu | 0.97 to 1.15 pu |

From thousands of buses at 2 to 7 pu Newton does not come back on
case9241pegase (the mismatch grows from 533 to 10900, 1.4e6, 9.6e10 and
is non-finite after 16 steps). On case3120sp and case2869pegase it does
come back, at the price of two to three extra steps. From the file's
stored voltages the corrections are small and both updates behave alike.

## What changes for the user

`power_flow.newton_update` selects the update: `polar` is the default
since 0.30.0, `rectangular` (the update up to 0.20.5) remains available
for comparisons; in the Web UI the select sits on the Settings page in
the Experimental block at the end of the Advanced options, greyed out
until "Enable experimental settings" is ticked (while greyed the
configured default applies).
With the polar update the rectangular Jacobian reproduces MATPOWER's
`newtonpf` iterates step for step.

- From the file's stored voltages both updates reach the same state in
  the same number of steps.
- From a flat start the polar update converges where the rectangular
  update diverged (case9241pegase in 7 iterations without damping) and
  saves one to three iterations on transmission cases; on small
  distribution feeders, whose angles are a few degrees, it costs at most
  one.
- Damping (`autodamp`, the merit line search, the trust region) scales
  the correction before it is applied, so every start machine and step
  control works with either update. See
  [Start strategies by case](start_strategies.md) for the recipes.

## The active set and the path

Switching a PV bus to PQ inside the Newton loop looks at the iterates on
the way, so a different update can give a different switching sequence
and, with limits binding, a different limited solution. The default
`power_flow.qlimits.start_iter: 3` with `start_mode: iteration_or_auto`
lets the iterate settle before the first switch and protects both
updates; switching from the second iterate is not safe with either. See
[Q-limit Switching Strategy](q_limit_switching_strategy.md).
