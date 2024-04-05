Network Branch Model
=============
Each branch is treated with the same four-terminal network model. It is a four-terminal network with an ideal transformer connected upstream. For power lines, the transmission ratio \(N\) is set to 1. For transformers, the transformation ratio \(N\) is given as a complex value. The admittance matrix \(Y\) looks like this:


\[
Y_{br} = \begin{bmatrix}
    \frac{1}{{\tau^2}} \cdot (y_s + j\frac{b_c}{2}) & -y_s \cdot \frac{1}{{\tau e^{-j\phi}}} \\
    -y_s \cdot \frac{1}{{\tau e^{j\phi}}} & (y_s + j\frac{b_c}{2})
\end{bmatrix}
\]




where:
- \( y_s = \frac{1}{R + jX} \) is the series admittance,
- \( R \) is the resistance component, and
- \( X \) is the reactance component,
- \( b_c \) is the transverse admittance,
- \( N = \tau \cdot e^{j\phi}\) is a complex transformation factor (eg 1 for power lines)



##### Circuit diagram

```plaintext	   
                    ys
      x--┓┏---------###----------x
         ||   |             |
         ||   # jbc/2       # jbc/2
         ||   |             |
      x--┛┗----------------------x 
         N (complex number)             

```
<!-- Dies ist ein auskommentierter Abschnitt -->
<!--┏
<! ┣
<! ┗
<! ┓
<! ┃
<! ┛
\( Y_{i0} \) fgfdgdfgfdg
-->


