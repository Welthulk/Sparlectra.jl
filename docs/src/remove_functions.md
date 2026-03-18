# Component Removal

## Conceptual model

### 1. Bus removal is a validity check first

`removeBus!` is intentionally different from the other helpers. Because the
`Net` struct is immutable at the top level, the function is mainly used to
check whether a bus is safe to remove at all.

In practice this means a bus cannot be removed when it is still used as:

- a slack bus,
- a connection point for branches,
- a connection point for prosumers, or
- a connection point for shunts.

### 2. Branch-oriented removal changes topology

`removeBranch!`, `removeACLine!`, and `removeTrafo!` modify the active network
contents. After one of these operations, the network may split into multiple
parts or leave single buses disconnected from the energized graph.

### 3. Isolation handling is part of the workflow

Use:

- `markIsolatedBuses!` to detect and inspect the current topology, and
- `clearIsolatedBuses!` to remove buses that have become removable.

This is especially helpful after several branch removals or contingency-style
editing steps.

### 4. Validation closes the loop

After any structural change, run `validate!` before the next solver call. That
keeps the editing workflow predictable and makes topology problems visible
before Newton-Raphson or state-estimation routines are executed.

## Where to find what

- For a practical editing sequence with code examples, see the
  [Workshop](workshop.md).
- For generated API docs of `removeBus!`, `removeBranch!`, `removeACLine!`,
  `removeTrafo!`, `removeShunt!`, `removeProsumer!`, and
  `clearIsolatedBuses!`, see the [Function Reference](reference.md).
