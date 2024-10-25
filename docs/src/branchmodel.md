Network Branch Model
=============
Each branch is treated with the same four-terminal network model. It is a four-terminal network with an ideal transformer connected upstream. For power lines, the transmission ratio \(N\) is set to 1. For transformers, the transformation ratio \(N\) is given as a complex value. The admittance matrix \(Y\) looks like this:

```math
Y_{br} = \begin{bmatrix}
    \frac{1}{{\tau^2}} \cdot (y_{ser} + 0.5 \cdot y_{shunt}) & -y_{ser} \cdot \frac{1}{{\tau e^{-j\phi}}} \\
    -y_{ser} \cdot \frac{1}{{\tau e^{j\phi}}} & (y_{ser} + 0.5 \cdot y_{shunt})
\end{bmatrix}
```

where:
-  ``y_ser``  is the series admittance,
-  ``y_shunt``  is the shunt admittance,
-  ``R``  is the resistance component, and
-  ``X``  is the reactance component,
-  ``G``  is the conductance component, and
-  ``B``  is the susceptance component,
-  ``N `` is the complex transformation factor (eg 1 for power lines)
```math
N = \tau \cdot e^{j\phi}
```

```math
y_{ser} = \frac{1}{R + jX}
```
```math
y_{shunt} = G + j \cdot B
```


#### Circuit diagram
```

                    y_ser
      x--┓┏---------###----------x
         ||   |             |
         ||   # y_shunt     # y_shunt
         ||   |             |
      x--┛┗----------------------x 
         N = tau * e^(j*phi)
```



