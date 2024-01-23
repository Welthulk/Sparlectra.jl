*Text: Udo Schmitz*
*Date: 01.11.2023*
*Version: 1.0*

# Structure of the Jacobian Matrix

In power system analysis, the Jacobian matrix is a fundamental tool used to represent the partial derivatives of a system of equations describing the power flow in an electrical network. The Jacobian matrix is particularly crucial in numerical methods like the Newton-Raphson algorithm for solving nonlinear power flow equations.

## Purpose

The purpose of the Jacobian matrix is to facilitate the calculation of the derivatives of real and reactive power injections with respect to the system variables, such as bus voltages (\(V\)) and phase angles (\(\phi\)). The Jacobian matrix is structured using superelements to enhance its clarity.

## Example

Consider the following graph:

    1 ----- 2
    |       |
    |       |
    3 ----- 4----5 (Slack-Bus)

For the 5-bus system, the Jacobian matrix \(\mathbf{J}\) has a size of \(\mathbf{8 \times 8}\), excluding the slack bus (Bus 5), and follows this structure:




\[ \mathbf{J} = \begin{bmatrix}
\mathbf{H}_{11} & \mathbf{J}_{12} & \mathbf{H}_{13} & \mathbf{J}_{14} & \mathbf{H}_{15} & \mathbf{J}_{16} & 0 & 0\\
\mathbf{N}_{21} & \mathbf{L}_{22} & \mathbf{N}_{23} & \mathbf{L}_{24} & \mathbf{N}_{25} & \mathbf{L}_{26} &  0 & 0\\
\mathbf{H}_{31} & \mathbf{J}_{32} & \mathbf{H}_{33} & \mathbf{J}_{34} & 0 & 0 &\mathbf{H}_{37} & \mathbf{J}_{38} \\
\mathbf{N}_{41} & \mathbf{L}_{42} & \mathbf{N}_{43} & \mathbf{L}_{44} & 0 & 0 & \mathbf{N}_{47} & \mathbf{L}_{48} \\
\mathbf{H}_{51} & \mathbf{J}_{52} & 0 & 0 & \mathbf{H}_{55} & \mathbf{J}_{56} & \mathbf{H}_{37} & \mathbf{J}_{38} \\
\mathbf{N}_{61} & \mathbf{L}_{62} & 0 & 0 & \mathbf{N}_{65} & \mathbf{L}_{66} &  \mathbf{N}_{47} & \mathbf{L}_{48} \\
0 & 0 &  \mathbf{H}_{73} & \mathbf{J}_{74} & \mathbf{H}_{75} & \mathbf{J}_{76} & \mathbf{H}_{77} &\mathbf{J}_{78} \\
0 & 0 &  \mathbf{N}_{83} & \mathbf{L}_{84} & \mathbf{N}_{85} &  \mathbf{L}_{86} & \mathbf{N}_{87} &\mathbf{L}_{88} \\
\end{bmatrix} \]

where:

- \(\mathbf{H}_{ij}\) represents the submatrix of \(\frac{\partial \mathbf{P}_i}{\partial \boldsymbol{\phi}_j}\).
- \(\mathbf{J}_{ij}\) represents the submatrix of \(\frac{\partial \mathbf{Q}_i}{\partial \boldsymbol{\phi}_j}\).
- \(\mathbf{N}_{ji}\) represents the submatrix of \(V_j \cdot \frac{\partial \mathbf{P}_i}{\partial \mathbf{V}_j}\).
- \(\mathbf{L}_{ij}\) represents the submatrix of \(V_j \cdot \frac{\partial \mathbf{Q}_i}{\partial \mathbf{V}_j}\).

Each node of the graph forms a supercell in the Jacobian matrix with elements \(\mathbf{H}_{ij}\), \(\mathbf{J}_{ij}\), \(\mathbf{N}_{ji}\), and \(\mathbf{L}_{ij}\). The nodes are arranged in the order of bus numbers. The node-branch relationships are depicted in the following table:



  | Node  | **1**                  | **2**                  | **3**                  | **4**                  |  
  |-------|------------------------|------------------------|------------------------|------------------------|
  | **1** | \(H_{11}\)  \(J_{12}\) | \(H_{13}\)  \(J_{14}\) | \(H_{15}\)  \(J_{16}\) | \(0\)   \(0\)  |
  |       | \(N_{21}\)  \(L_{22}\) | \(N_{23}\)  \(L_{24}\) | \(N_{25}\)  \(L_{26}\) | \(0\)   \(0\)  |  
  | **2** | \(H_{31}\)  \(J_{32}\) | \(H_{33}\)  \(J_{34}\) | \(0\)   \(0\)  | \(H_{37}\)  \(J_{38}\) |
  |       | \(N_{41}\)  \(L_{42}\) | \(N_{43}\)  \(L_{44}\) | \(0\)   \(0\)  | \(N_{47}\)  \(L_{48}\) |
  | **3** | \(H_{51}\)  \(J_{52}\) | \(0\)   \(0\)  | \(H_{55}\)  \(J_{56}\) | \(H_{57}\)  \(J_{58}\) |
  |       | \(N_{61}\)  \(L_{62}\) | \(0\)   \(0\)  | \(N_{65}\)  \(L_{66}\) | \(N_{67}\)  \(L_{68}\) |
  | **4** | \(0\)   \(0\)  | \(H_{73}\)  \(J_{74}\) | \(H_{75}\)  \(J_{76}\) | \(H_{77}\)  \(J_{78}\) |
  |       | \(0\)   \(0\)  | \(N_{83}\)  \(L_{84}\) | \(N_{85}\)  \(L_{86}\) | \(N_{87}\)  \(L_{88}\) |


An example for this graph is presented in the following table, showcasing the results of the first iteration through the utilization of the `printJacobian` function:
```Julia
Jacobian: 8 x 8, number of PQ nodes:4,  number of PV nodes: 0
        1:      2:      3:      4:      5:      6:      7:      8:
1:      19.6    10.1    -9.8    -5.0    -9.8    -5.0    0.0     0.0
2:      -10.1   19.6    5.0     -9.8    5.0     -9.8    0.0     0.0
3:      -9.8    -5.0    19.6    10.0    0.0     0.0     -9.8    -5.0
4:      5.0     -9.8    -10.1   19.6    0.0     0.0     5.0     -9.8
5:      -9.8    -5.0    0.0     0.0     19.6    10.0    -9.8    -5.0
6:      5.0     -9.8    0.0     0.0     -10.1   19.6    5.0     -9.8
7:      0.0     0.0     -9.8    -5.0    -9.8    -5.0    29.4    15.1
8:      0.0     0.0     5.0     -9.8    5.0     -9.8    -15.1   29.4
```
>**Note**:
The Jacobian matrix is real-valued.
## Bus Types
There are essentially three types of bus types in a power system. The PQ node with the unknown voltage V and the unknown angle phi, voltage-controlled PV nodes with the unknown reactive power Q and angle phi, and the reference node (Slack) where voltage and angle are specified. For voltage-regulated nodes, the derivatives with respect to voltage are omitted. Therefore, the elements N and J are zero. These rows can be removed from the Jacobian matrix. The slack node is not considered in the Jacobian matrix. The following table shows the bus types and the corresponding unknowns and Jacobian elements:

| Node Type  | Unknowns                   | Jacobi-Elements|
|------------|----------------------------|----------------------------
| PQ Node    | V, phi                     | Hij, Jij, Nij, Lij                 |
| PV Node    | Q, phi                     | Hij, Jij, Jii=0, Nij, Lij =0  |
| Slack Node | V (prescribed), phi (prescribed)| -

### Example 2
    1 ----- 2 (PV)
    |       |
    |       |
    3 ----- 4----5 (Slack-Bus)
           (PV)

\[ \mathbf{J} = \begin{bmatrix}
\mathbf{H}_{11} & \mathbf{J}_{12} & \mathbf{H}_{13} & \mathbf{0} & \mathbf{H}_{15} & \mathbf{J}_{16} & 0 & 0\\
\mathbf{N}_{21} & \mathbf{L}_{22} & \mathbf{N}_{23} & \mathbf{0} & \mathbf{N}_{25} & \mathbf{L}_{26} &  0 & 0\\
\mathbf{H}_{31} & \mathbf{J}_{32} & \mathbf{H}_{33} & \mathbf{0} & 0 & 0 &\mathbf{H}_{37} & \mathbf{0} \\
\mathbf{0}  & \mathbf{0} & \mathbf{0} & \mathbf{0} & 0 & 0 & \mathbf{0} & \mathbf{0} \\
\mathbf{H}_{51} & \mathbf{J}_{52} & 0 & 0 & \mathbf{H}_{55} & \mathbf{J}_{56} & \mathbf{H}_{37} & \mathbf{0} \\
\mathbf{N}_{61} & \mathbf{L}_{62} & 0 & 0 & \mathbf{N}_{65} & \mathbf{L}_{66} &  \mathbf{N}_{67} & \mathbf{0} \\
0 & 0 &  \mathbf{H}_{73} & \mathbf{0} & \mathbf{H}_{75} & \mathbf{J}_{76} & \mathbf{H}_{77} &\mathbf{0} \\
0 & 0 &  \mathbf{0} & \mathbf{0} & \mathbf{0} &  \mathbf{0} & \mathbf{0} &\mathbf{0} \\
\end{bmatrix} \]

All rows with \(0\) values can be eliminated, resulting in:

\[ \mathbf{J} = \begin{bmatrix}
\mathbf{H}_{11} & \mathbf{J}_{12} & \mathbf{H}_{13} & \mathbf{H}_{15} & \mathbf{J}_{16} & 0\\
\mathbf{N}_{21} & \mathbf{L}_{22} & \mathbf{N}_{23}  & \mathbf{N}_{25} & \mathbf{L}_{26} & 0\\
\mathbf{H}_{31} & \mathbf{J}_{32} & \mathbf{H}_{33}  & 0 & 0 &\mathbf{H}_{37}  \\
\mathbf{H}_{51} & \mathbf{J}_{52} & 0 &  \mathbf{H}_{55} & \mathbf{J}_{56} & \mathbf{H}_{37}  \\
\mathbf{N}_{61} & \mathbf{L}_{62} & 0 &  \mathbf{N}_{65} & \mathbf{L}_{66} &  \mathbf{N}_{67}  \\
0 & 0 &  \mathbf{H}_{73}  & \mathbf{H}_{75} & \mathbf{J}_{76} & \mathbf{H}_{77}  \\
\end{bmatrix} \]


An example for this graph is presented in the following table, showcasing the results of the first iteration using the `printJacobian` function. The displayed output corresponds to the example with 2 PV nodes:
```Julia
Jacobian: 6 x 6, number of PQ nodes:2,  number of PV nodes: 2
        1:      2:      3:      4:      5:      6:
1:      19.6    10.0    -9.8    -9.8    -5.0    0.0
2:      -10.1   19.6    5.0     5.0     -9.8    0.0
3:      -9.8    -5.0    19.6    0.0     0.0     -9.8
4:      -9.8    -5.0    0.0     19.6    10.0    -9.8
5:      5.0     -9.8    0.0     -10.1   19.6    5.0
6:      0.0     0.0     -9.8    -9.8    -5.0    29.5
```
## Number of Superelements
The size of the Jacobian matrix is determined by the relationship: 
\[ N_s = 2 \cdot (N - 1) - N_{\text{PV}} \] 
where \( N \) is the total number of buses and \( N_{\text{PV}} \) is the number of PV buses.




<!-- 
## Superelements Formulation

For a given bus \(i\) and its neighboring bus \(j\), the superelement matrices are calculated as follows:

\[
\begin{align*}
\mathbf{H}_{ij} &= \frac{\partial P_i}{\partial \phi_j} = +v_i \cdot [g_{ij} \cdot v_j \cdot \sin(\phi_i - \phi_j - \alpha_{ij})] \\
\mathbf{J}_{ij} &= \frac{\partial Q_i}{\partial \phi_j} = -v_i \cdot [g_{ij} \cdot v_j \cdot \cos(\phi_i - \phi_j - \alpha_{ij})] \\
\mathbf{N}_{ji} &= V_j \cdot \frac{\partial P_i}{\partial V_j} = +v_i \cdot [g_{ij} \cdot v_j \cdot \cos(\phi_i - \phi_j - \alpha_{ij})] \\
\mathbf{L}_{ij} &= V_j \cdot \frac{\partial Q_i}{\partial V_j} = +v_i \cdot [g_{ij} \cdot v_j \cdot \sin(\phi_i - \phi_j - \alpha_{ij})]
\end{align*}
\]
-->
## Relationship to Ybus Matrix

It's noteworthy that the Jacobian matrix has a similar structure to the Ybus matrix when considering the superelement as a single element. The Ybus matrix represents the nodal admittances in the power system network, and the Jacobian matrix characterizes how changes in bus voltages and phase angles influence real and reactive power injections.

## `calcJacobian` Function
### Function Signature:

```Julia
function calcJacobian(
    Y::AbstractMatrix{ComplexF64},
    busVec::Vector{BusData},
    adjBranch::Vector{Vector{Int}},
    busTypeVec::Vector{ResDataTypes.NodeType},
    slackIdx::Int,
    n_pq::Int,
    n_pv::Int,
    log::Bool = false,
    sparse::Bool = true
)
```

### Parameters:

1. **Y::AbstractMatrix{ComplexF64}**: The admittance matrix of the power system in complex numbers.

2. **busVec::Vector{BusData}**: A vector containing information about each bus in the power system. Each element is of type `BusData`.

3. **adjBranch::Vector{Vector{Int}}**: A vector of vectors representing the adjacency information of branches in the power system. Each sub-vector contains the indices of buses connected by a branch.

4. **busTypeVec::Vector{ResDataTypes.NodeType}**: A vector specifying the type of each bus, indicating whether it's a PQ bus, a PV bus, or a slack bus. Each element is of type `ResDataTypes.NodeType`.

5. **slackIdx::Int**: The index of the slack bus in the system.

6. **n_pq::Int**: The total number of PQ (constant power) buses.

7. **n_pv::Int**: The total number of PV (voltage-controlled power) buses.

8. **log::Bool = false**: A boolean flag indicating whether to log additional information during the calculation.

9. **sparse::Bool = true**: A boolean flag indicating whether to use sparse matrix representations for computational efficiency. 