# Component Types

This page documents the main component types used in Sparlectra for power system modeling.

## Network Components

### Net

The main container for a power system network.

```
Sparlectra.Net
```

### Node

Represents a bus in the power system.

```
Sparlectra.Node
```

### Branch Components

```
Sparlectra.Branch
Sparlectra.BranchFlow
Sparlectra.BranchModel
```

Since r0.10.0 a `Branch` carries per-terminal service flags
`from_status`/`to_status` (1 = closed, 0 = open) next to the aggregate
`status`. Invariant: `status = 1` iff both terminals are closed;
`setBranchStatus!` (the user-facing switch) sets all three consistently,
`setBranchTerminalStatus!(br; from =, to =)` toggles individual terminals
and recomputes the aggregate. A branch open at exactly one terminal stays
in the model as its exact pi reduction at the closed bus and reports the
open-end voltage (`open_end_vm_pu`/`open_end_va_deg`) as a result, see the
"One-sided open branches" section of the [branch model](branchmodel.md).

### Prosumer Components

```
Sparlectra.ProSumer
```

### Shunt Components

```
Sparlectra.Shunt
```

### Transformer Components

```
Sparlectra.PowerTransformer
Sparlectra.PowerTransformerWinding
Sparlectra.PowerTransformerTaps
```

### Line Components

```
Sparlectra.ACLineSegment
```

### Link Components

```
Sparlectra.BusLink
```

## Basic Components

```
Sparlectra.Component
Sparlectra.ImpPGMComp
Sparlectra.ImpPGMComp3WT
```

## Enumerations

These enumerations are used to define types and states for various components:

```
Sparlectra.ComponentTyp
Sparlectra.TrafoTyp
Sparlectra.NodeType
Sparlectra.ProSumptionType
```
