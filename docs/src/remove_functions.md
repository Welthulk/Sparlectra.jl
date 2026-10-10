# Component Removal

## Conceptual model

- `removeBus!` only checks whether a bus can be removed and returns the
  answer; it does not remove the bus. A slack bus, or a bus with connected
  branches, prosumers or shunts, cannot be removed.
- `removeBranch!`, `removeACLine!` and `removeTrafo!` change the network;
  afterwards it can fall into several parts or leave single buses
  disconnected.
- `markIsolatedBuses!` detects the current topology, `clearIsolatedBuses!`
  removes the buses that became removable, for example after several
  branch removals or contingency-style edits.
- Run `validate!` after any structural change and before the next solver
  call, so topology problems surface before the power flow or the state
  estimation runs.

## Where to find what

- Editing sequence with code: [workshop tour](generated/workshop_tour.md).
- API docs of `removeBus!`, `removeBranch!`, `removeACLine!`,
  `removeTrafo!`, `removeShunt!`, `removeProsumer!` and
  `clearIsolatedBuses!`: [Function Reference](reference.md).
