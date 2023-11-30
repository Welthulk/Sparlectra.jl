
*Text: Udo Schmitz*
*Date: 22.10.2023*
*Version: 1.0*

# Ybus
The Ybus matrix serves as a comprehensive representation of the electrical network, enabling the calculation of power flow within the system. Also known as the nodal admittance matrix, Ybus is a square matrix whose dimensions correspond to the total number of nodes in the network. Comprising complex values, Ybus features real admittance values along its diagonal and imaginary admittance values off-diagonally. The matrix is sparse due to the significantly lower number of connections between nodes compared to the total number of nodes in the network. Additionally, Ybus is symmetric, reflecting equal admittance values between two nodes in both directions.
## Example
The following example shows a simple network with 3 nodes and 3 branches. The network is represented by the following graph:


    1 ┃-----------------L1_2------------------┃ 2
      ┃-┓                                   ┏-┃
        |                                   |
        ┗---L1_3-------┓  ┏-----L3_2--------┛
                     __|__|__                      
                        3

Each line is described by the following parameters:

- Resistance \( R \)
- Inductance  \( X \)
- Capacity  \( C \)
- leakage conductance \( G \)

The four-terminal network representation of a line segment is as follows:



\[ \begin{bmatrix} U_1 \\ U_2 \end{bmatrix} = \begin{bmatrix} Y_{11} & Y_{12} \\ Y_{21} & Y_{22} \end{bmatrix} \cdot \begin{bmatrix} I_1 \\ I_2 \end{bmatrix} \]

Neglecting \(G\) results in:

\[ \begin{bmatrix} U_1 \\ U_2 \end{bmatrix} = \begin{bmatrix} \frac{1}{R} + j\omega L & -j\omega C \\ -j\omega C & \frac{1}{R} + j\omega L \end{bmatrix} \cdot \begin{bmatrix} I_1 \\ I_2 \end{bmatrix} \]

For the example above the Ybus has 3 rows and 3 columns and is therefore square. The rows represent the buses (in our example 1 to 3), the columns represent the connections between the buses. The matrix is perfectly symmetric if passive equivalent circuits exist for the branches. This is always the case for lines and longitudinal regulator transformers. Then only the upper triangular matrix needs to be calculated.



\[ Y = \begin{bmatrix} Y_{11} & Y_{12} & Y_{13} \\ Y_{21} & Y_{22} & Y_{23} \\ Y_{31} & Y_{32} & Y_{33} \end{bmatrix} \]

If the line between bus 1 and bus 3 were to fail, the Ybus would look as follows:
\[ Y = \begin{bmatrix} Y_{11} & Y_{12} & 0 \\ Y_{21} & Y_{22} & Y_{23} \\ 0 & Y_{32} & Y_{33} \end{bmatrix} \]


The non-diagonal elements represent the connections between the buses. The admittance is formed from the sum of the longitudinal admittances of the branch. The diagonal elements are formed from the sum of the admittances of the connected connections including the transverse admittances.


The rule for forming the diagonal elements (\(Y_{ii}\)) is:

\[ Y_{ii} = Y_{i0} + \sum_{j \neq i} Y_{ij} \]

This means that the diagonal element \(Y_{ii}\) is equal to the sum of all admittances between bus \(i\) and all other buses (\(j \neq i\)) including the conductance \( Y_{i0} \).  

For the non-diagonal follows:
\[ Y_{ij} = Y_{ji} = -\frac{1}{Z_{ij}} \]

>**Note**:
For real transmission networks, the admittance matrix is often sparsely populated, meaning many elements are zero. Therefore, it is advisable to use sparse matrix storage for efficient memory utilization. For instance, Julia provides support for sparse matrices.

### Numerical example
The admittance matrix with specific values for the 3-bus system with a failed line \(L_{13}\) is as follows:


##### Line Parameters:
Length: \(1 \, \text{km}\)
Frequency: \(50 \, \text{Hertz}\)
\(B' = 3.18 \times 10^{-7} \, \Omega/km\)
\[
Y_{line} = \frac{1}{Z} = \frac{1}{j \omega B'} = -j \omega B' = -2\pi \times 50\text{Hz} \times B' \times 1000 \, \text{m} = -0.1 Siemens
\]

\(L_{13}\) = \(L_{23}\) = 0
\(L_{12}\) = \(L_{12}\) = \(L_{23}\) =\(L_{32}\) =-0.1 Siemens

######result:
\[ Y = \begin{bmatrix} Y_{11} & Y_{12} & Y_{13} \\ Y_{21} & Y_{22} & Y_{23} \\ Y_{31} & Y_{32} & Y_{33} \end{bmatrix} = \begin{bmatrix} j0.2 & -j0.1 & 0 \\ -j0.1 & j0.2 & -j0.1 \\ 0 & -j0.1 & j0.2 \end{bmatrix} =
\begin{bmatrix} 0.1 & -0.1 & 0 \\ 0.1 & 0.2 & -0.1 \\ 0 & -0.1 & 0.2 \end{bmatrix} \cdot e^{j\frac{\pi}{2}}
\]

### Function to calculate Ybus in Sparlectra
The Ybus matrix is computed in the `CreateYBUS` function within the [`equicircuit.jl` file](../../src/equicircuit.jl). The function takes the following arguments:

- `branchVec::Vector{Branch}`: A vector containing branches.
- `shuntVec::Vector{ResDataTypes.Shunt}`: A vector containing shunts.
- `sparse::Bool = true`: A boolean flag indicating whether to use sparse matrix storage. Default is `true`.
- `printYBUS::Bool = false`: A boolean flag indicating whether to print the YBUS matrix. Default is `false`."





<!-- Dies ist ein auskommentierter Abschnitt -->
<!--┏
<! ┣
<! ┗
<! ┓
<! ┃
<! ┛
\( Y_{i0} \) fgfdgdfgfdg
-->





     
