# Component Types

## Network Components

### Net

The container of a power system network.

```
Sparlectra.Net
```

### Node

A bus of the power system.

```
Sparlectra.Node
```

### Branch Components

```
Sparlectra.Branch
Sparlectra.BranchFlow
Sparlectra.BranchModel
```

A `Branch` carries the per-terminal flags `from_status`/`to_status`
(1 = closed, 0 = open) next to the aggregate `status` (`1` iff both
terminals are closed); `setBranchStatus!` sets all three,
`setBranchTerminalStatus!(br; from =, to =)` toggles single terminals.
Branches open at one terminal: "One-sided open branches" in the
[branch model](branchmodel.md).

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

```
Sparlectra.ComponentTyp
Sparlectra.TrafoTyp
Sparlectra.NodeType
Sparlectra.ProSumptionType
```
