*Text: Udo Schmitz*
*Date: 04.12.2023*
*Version: 1.0*

# Minimum Degree Ordering in Graphs

Minimum Degree Ordering is a graph ordering algorithm used to reduce the fill-in during the factorization of sparse matrices. This technique is often applied in numerical linear algebra and plays a crucial role in optimizing the efficiency of algorithms like Cholesky factorization and solving systems of linear equations.

### Fill-in Example
\[
\text{A} = \begin{bmatrix}
    1 & 1 & 1 & 1 \\
    1 & 1 & 0 & 0 \\
    1 & 0 & 1 & 0 \\
    1 & 0 & 0 & 1 \\
\end{bmatrix}
\]

The LU factorization of A is:
```Julia
F=lu(A)
LU{Float64, Matrix{Float64}, Vector{Int64}}
L factor:
4×4 Matrix{Float64}:
 1.0   0.0  0.0  0.0
 1.0   1.0  0.0  0.0
 1.0  -0.0  1.0  0.0
 1.0   1.0  1.0  1.0
U factor:
4×4 Matrix{Float64}:
 1.0   1.0   1.0   1.0
 0.0  -1.0   0.0  -1.0
 0.0   0.0  -1.0  -1.0
 0.0   0.0   0.0   2.0
```

upper and lower triangular matrix are almost fully occupied.

An optimized order is:

\[
\text{B} = \begin{bmatrix}
    1 & 0 & 0 & 1 \\
    0 & 1 & 0 & 1 \\
    0 & 0 & 1 & 1 \\
    1 & 1 & 1 & 1 \\
\end{bmatrix}
\]

```Julia
F=lu(B)
LU{Float64, Matrix{Float64}, Vector{Int64}}
L factor:
4×4 Matrix{Float64}:
 1.0  0.0  0.0  0.0
 0.0  1.0  0.0  0.0
 0.0  0.0  1.0  0.0
 1.0  1.0  1.0  1.0
U factor:
4×4 Matrix{Float64}:
 1.0  0.0  0.0   0.0
 0.0  1.0  0.0   1.0
 0.0  0.0  1.0   1.0
 0.0  0.0  0.0  -1.0
```




Minimum Degree Order algorithms are employed to determine an almost optimal elimination sequence. Finding the truly optimal solution is often impractical, as the computational cost of identifying the optimal order can surpass the time it takes to solve the system of equations using a less optimized sequence. Therefore, the Minimum Degree Ordering algorithm is used to find a good approximation of the optimal order in a reasonable amount of time.


## Algorithm Overview

The Minimum Degree Ordering algorithm follows these basic steps:

1. **Initialize Degrees**: Compute the degree of each node in the graph, where the degree is the number of edges connected to a node.

2. **Elimination List**: Create an elimination list containing all nodes.

3. **Ordering Iteration**: In each iteration, select the node with the minimum degree from the elimination list.

4. **Elimination**: Eliminate the selected node by updating the degrees of its neighbors.

5. **Repeat**: Continue the process until all nodes are eliminated.

## Minimum Degree Order (MDO) Algorithm

1. **Input**: Graph represented by `nodeNumberVec` and `branchTupleVec`.
2. **Output**: An array `order` representing the nodes in the order to minimize fill-in.

### Pseudo-Code:

```plaintext
function minimumDegreeOrder(graph) -> Array
    N = Number of nodes in the graph
    degree = List of node degrees (degree[i] is the degree of node i)
    order = Empty list for the order of nodes

    while N > 0:
        // Find the node with the minimum degree
        minDegreeNode = findMinDegreeNode(degree)

        // Add the node to the order
        order.append(minDegreeNode)

        // Reduce the degree of neighboring nodes
        for neighbor in Neighbors of minDegreeNode:
            degree[neighbor] = degree[neighbor] - 1

        // Mark the considered node as removed
        degree[minDegreeNode] = -1

        N = N - 1

    return order
```
## Example

Consider a simple graph:

    1 -- 2
    | \  | 
    3 -- 4

Let's represent this graph using a Julia array:

```julia
nodeNumberVec = [1, 2, 3, 4]
branchTupleVec = [(1, 2), (1, 3), (1, 4), (2, 1), (2,4), (3, 1), (3,4), (4, 1), (4,3), (4,2)]
n = length(nodeNumberVec)
```
The Adjacency Matrix of the graph will be:
\[
\text{A} = \begin{bmatrix}
    1 & 1 & 1 & 1 \\
    1 & 1 & 0 & 1 \\
    1 & 0 & 1 & 1 \\
    1 & 1 & 1 & 1 \\
\end{bmatrix}
\]

Now, apply the Minimum Degree Ordering algorithm:
```julia
order = mdoRCM(n, branchTupleVec)
@show order
```
The output will be:
```plaintext
2
1
3
4
```
The rearanged Matrix will be:
\[
\text{A} = \begin{bmatrix}
    1 & 1 & 0 & 1 \\
    1 & 1 & 1 & 1 \\
    1 & 0 & 1 & 1 \\
    1 & 1 & 1 & 1 \\
\end{bmatrix}
\]

>**Note**:
The true strengths of the RCM algorithm become evident when dealing with large sparse matrices. The example above is merely intended to illustrate the approach.
##Implementation in Sparlectra
A implementation of the Minimum Degree Ordering algorithm is provided in the Sparlectra package. The implementation is based on the Reverse Cuthill-McKee algorithm (see: https://en.wikipedia.org/wiki/Cuthill%E2%80%93McKee_algorithm), which is a variation of the Minimum Degree Ordering algorithm. The Reverse Cuthill-McKee algorithm is a heuristic algorithm that is often used to find a good approximation of the optimal order in a reasonable amount of time. The algorithm is encapsulated in the file `nbi.jl`. The key function, `mdoRCM`, takes the network represented by `nodeNumberVec` and `branchTupleVec` as input and returns an array representing the nodes ordered to minimize fill-in.

### Integration with Matpower-Import

In Sparlectra, the solution of the equation system utilizes the `\` operator, where the impact of the Minimum Degree Order, based on my investigations, is relatively minimal. However, for users who wish to experiment and explore the effects of Minimum Degree Ordering, there is an option available during the Matpower import process.

When using the `createNetFromMatPowerFile` function, users can set the `mdo` flag to `true` to enable Minimum Degree Ordering. This allows users to conduct experiments and observe the influence of Minimum Degree Ordering on the network analysis based on their specific use cases.

```julia
using sparlectra
using sparlectra.ResDataTypes
using sparlectra.sparlectraimport
using sparlectra.sparlectranet

# Example Matpower import with Minimum Degree Ordering
base_MVA = 0.0 # # Set a value to a non-zero value to override the corresponding value in the case file.
casefile = "a_case.m"
jpath = joinpath(pwd(), "sparlectra", "data", "mpower", casefile)

myNet = createNetFromMatPowerFile(jpath, base_MVA, (verbose >= 1), mdo)
```



