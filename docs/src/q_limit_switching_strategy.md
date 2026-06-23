# Q-limits in power flow: why switching strategy matters

In power-flow analysis, reactive-power limits of generators may look like a minor modelling detail at first. In practice, however, they can decide whether a network case converges robustly or enters a switching cascade.

A PV bus keeps its voltage magnitude fixed in the classical power-flow formulation. The reactive power required to maintain that voltage is determined by the solution. If the required reactive power exceeds the generator's admissible range, the bus can no longer be treated as an unconstrained voltage-controlled PV bus. It must be converted into a PQ bus with reactive power fixed at the violated limit.

This PV-to-PQ conversion is numerically delicate. Two basic strategies are useful to distinguish.

The **active-set strategy** reacts to Q-limit violations during the nonlinear iteration itself. This can be flexible, but it also carries a risk. If many generators have narrow or nearly identical Q ranges, many PV-to-PQ and PQ-to-PV events can occur. In that situation, the algorithm is no longer only solving a smooth nonlinear power-flow problem. It is also solving a discrete switching problem.

The **classical strategy** is more conservative. First, the base power flow is solved without active Q-limit switching. Only after this base run has converged are reactive-power limit violations checked. Violating generators are clamped to `Qmin` or `Qmax`, the affected buses are converted to PQ, and the power flow is solved again. This can be done simultaneously for all violating generators or step by step for the largest violation. A bus that has been converted from PV to PQ is not automatically converted back to PV within the same enforcement loop.

Large synthetic or aggregated networks can expose the difference clearly. Many PV buses, narrow Q bands, different voltage start values, and generator setpoints may make it difficult for the numerical iteration to distinguish cleanly between voltage control and reactive-power limitation. In such cases, the network model is not necessarily wrong. The switching strategy itself may be the decisive factor.

Sparlectra therefore provides more than one Q-limit enforcement mode.

`power_flow.qlimits.enforcement_mode` selects one of these canonical modes:

- `active_set` is the default dynamic mode. It performs in-iteration Q-limit switching and can use guards such as hysteresis, cooldown, narrow-range locking, and repeated-switching protection.
- `classic_simultaneous` is a classical reference mode. It first solves the base power flow with Q-limit switching disabled. If the base solution converges, all detected Q-limit violations are clamped and converted in one outer-loop pass.
- `classic_one_at_a_time` follows the same classical principle, but handles only the largest violation per outer-loop pass. This can make the switching sequence easier to inspect.

A practical diagnostic workflow is to compare modes. If a case fails in `active_set` because the active set changes repeatedly, but behaves more clearly in a classical mode, the problem may be dominated by discrete switching rather than by the continuous Newton iteration alone. If the base power flow itself does not converge in a classical mode, the issue is upstream of Q-limit enforcement and should be investigated separately.

The key point is that Q-limits are not just post-processing. They change the structure of the power-flow problem. Reliable large-network analysis therefore needs not only a Newton solver, but also a controlled strategy for the discrete switching between PV and PQ.
