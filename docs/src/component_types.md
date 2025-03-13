# Component Types

This page documents the main component types used in Sparlectra for power system modeling.

## Network Components

### Net

The main container for a power system network.

```@docs
Sparlectra.Net
```

### Node

Represents a bus in the power system.

```@docs
Sparlectra.Node
```

### Branch Components

```@docs
Sparlectra.Branch
Sparlectra.BranchFlow
Sparlectra.BranchModel
```

### Prosumer Components

```@docs
Sparlectra.ProSumer
```

### Shunt Components

```@docs
Sparlectra.Shunt
```

### Transformer Components

```@docs
Sparlectra.PowerTransformer
Sparlectra.PowerTransformerWinding
Sparlectra.PowerTransformerTaps
```

### Line Components

```@docs
Sparlectra.ACLineSegment
```

## Basic Components

```@docs
Sparlectra.Component
Sparlectra.ImpPGMComp
Sparlectra.ImpPGMComp3WT
```

## Enumerations

These enumerations are used to define types and states for various components:

```@docs
Sparlectra.ComponentTyp
Sparlectra.TrafoTyp
Sparlectra.NodeType
Sparlectra.ProSumptionType
```