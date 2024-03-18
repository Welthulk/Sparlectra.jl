*Text: Udo Schmitz*
*Date: 01.11.2023*
*Version: 1.0*

# Network Matrices in Electrical Systems

## 1. Adjacency Matrix

In the context of electrical networks, consider buses as nodes and branches (lines and transformers) as connections between these nodes. The adjacency matrix represents the relationships between buses and branches. If there is a connection (branch) between two buses, the corresponding matrix entry is set to 1; otherwise, it remains at 0.

### Example

Consider an electrical network with buses (Nodes: A, B, C) and branches (Branches: AB, BC, CA). The adjacency matrix would look like this:

\[
\begin{bmatrix}
   & AB & BC & CA \\
A & 0 & 1 & 1 \\
B & 1 & 0 & 1 \\
C & 1 & 1 & 0 \\
\end{bmatrix}
\]

In this example, there is a connection between A and B, B and C, and C and A.

## 2. Node-Branch Incidence Matrix

The node-branch incidence matrix considers the direction of the branches. If a branch is connected to a node and the current flows out, the matrix entry is set to 1. If the current flows into the node, the entry is set to -1. Connections to the ground are not considered.

### Example 1

Consider an electrical network with nodes (Buses: A, B, C) and branches (Branches: AB, BC, CA). The node-branch incidence matrix would look like this:

\[
\begin{bmatrix}
   & AB & BC & CA \\
A & -1 & 0 & 1 \\
B & 1 & -1 & 0 \\
C & 0 & 1 & -1 \\
\end{bmatrix}
\]

In this example, the first column represents the branch AB, and since the current flows from A to B, the entry for A is -1, and for B, it's 1. Similar considerations apply to other branches.

## 3. YBus Matrix Construction

### Using Incidence Matrix and Diagonal Admittance Matrix

The YBus matrix represents the admittance of an electrical network and is crucial in power system analysis. It can be constructed using the node-branch incidence matrix (\(A\)) and the diagonal matrix of branch admittances (\(Y_d\)).

#### Steps to Construct YBus

1. **Compute Diagonal Admittance Matrix (\(Y_d\)):**

   \(Y_d\) is a diagonal matrix where each diagonal element corresponds to the admittance of the respective branch.

   \[
   Y_d =
   \begin{bmatrix}
      Y_{AB} & 0 & 0 \\
      0 & Y_{BC} & 0 \\
      0 & 0 & Y_{CA} \\
   \end{bmatrix}
   \]

2. **Calculate YBus:**

   The YBus matrix is calculated using the formula \(Y_{\text{Bus}} = A \cdot Y_d \cdot A^T\).

   \[
   Y_{\text{Bus}} = A \cdot Y_d \cdot A^T
   \]

   Where:
   - \(A\) is the node-branch incidence matrix.
   - \(Y_d\) is the diagonal admittance matrix.
   - \(A^T\) is the transpose of matrix \(A\).

#### Example 2

Let's consider an example with three buses (A, B, C) and three branches (AB, BC, CA). The incidence matrix (\(A\)) and diagonal admittance matrix (\(Y_d\)) are as follows:

   \[
   A =
   \begin{bmatrix}
      -1 & 0 & 1 \\
      1 & -1 & 0 \\
      0 & 1 & -1 \\
   \end{bmatrix}
   \]

   \[
   Y_d =
   \begin{bmatrix}
      Y_{AB} & 0 & 0 \\
      0 & Y_{BC} & 0 \\
      0 & 0 & Y_{CA} \\
   \end{bmatrix}
   \]

   The YBus matrix can be obtained by performing the matrix multiplication as described above.

   \[
   Y_{\text{Bus}} = A \cdot Y_d \cdot A^T
   \]

   Substituting the matrices with the example values:

   \[
   Y_{\text{Bus}} =
   \begin{bmatrix}
      -1 & 0 & 1 \\
      1 & -1 & 0 \\
      0 & 1 & -1 \\
   \end{bmatrix}
   \cdot
   \begin{bmatrix}
      Y_{AB} & 0 & 0 \\
      0 & Y_{BC} & 0 \\
      0 & 0 & Y_{CA} \\
   \end{bmatrix}
   \cdot
   \begin{bmatrix}
      -1 & 1 & 0 \\
      0 & -1 & 1 \\
      1 & 0 & -1 \\
   \end{bmatrix}
   \]

   Multiplying these matrices, we get:

   \[
   Y_{\text{Bus}} =
   \begin{bmatrix}
      Y_{AB} + Y_{CA} & -Y_{AB} & 0 \\
      -Y_{AB} & Y_{AB} + Y_{BC} & -Y_{BC} \\
      0 & -Y_{BC} & Y_{BC} + Y_{CA} \\
   \end{bmatrix}
   \]

This YBus matrix represents the network's admittance and is essential for power flow and system stability analysis.

> **Note:** In the Sparlectra project, the YBus matrix is created directly. Details about the creation of the YBus matrix and its implementation can be found in the `ybus.md` file within the Sparlectra project. For the calculation of the Jacobian matrix, it is recommended to use the 'adjacentBranches()' function.

### Functions to create NBI and Adjacency Matrices in Sparlectra

The function `adjacentBranches`, used for creating the adjacency matrix, can be found in the file `equicircuit.jl`. The function expects the following parameters:

- `Y`: An AbstractMatrix of ComplexF64 representing the admittance matrix of the electrical network.
- `log`: An optional Boolean parameter (default is `false`) that, when set to `true`, prints additional information about the adjacent branches.

The function `getNBI`, used for constructing the Node-Branch Incidence (NBI) matrix, is implemented in the file `nbi.jl`. The function expects the following parameters:

- `nodeNumberVec`: A vector containing the numbers of nodes in the electrical network.
- `branchTupleVec`: A vector of tuples representing the branches in the electrical network.

The function returns an integer matrix representing the Node-Branch Incidence (NBI) matrix, where rows correspond to nodes, columns correspond to branches, and the entries indicate the direction of the branches in relation to the nodes. Positive entries denote branches entering the node, and negative entries denote branches leaving the node.
