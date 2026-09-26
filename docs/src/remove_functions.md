# Component Removal

## Conceptual model

### 1. Bus removal is a validity check first

The `Net` struct is immutable at the top level, so `removeBus!` mainly
checks whether a bus is safe to remove. A bus cannot be removed while it
is a slack bus or a connection point for branches, prosumers or shunts.

### 2. Branch-oriented removal changes topology

`removeBranch!`, `removeACLine!` and `removeTrafo!` modify the active
network; afterwards it may split into several parts or leave single buses
disconnected.

### 3. Isolation handling is part of the workflow

`markIsolatedBuses!` detects the current topology, `clearIsolatedBuses!`
removes buses that have become removable; useful after several branch
removals or contingency-style edits.

### 4. Validation closes the loop

Run `validate!` after any structural change and before the next solver
call, so topology problems surface before Newton-Raphson or
state-estimation routines run.

## Where to find what

- Editing sequence with code: [workshop tour](generated/workshop_tour.md).
- API docs of `removeBus!`, `removeBranch!`, `removeACLine!`,
  `removeTrafo!`, `removeShunt!`, `removeProsumer!` and
  `clearIsolatedBuses!`: [Function Reference](reference.md).
