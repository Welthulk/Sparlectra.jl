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