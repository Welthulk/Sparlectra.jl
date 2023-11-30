*Text: Udo Schmitz*
*Date: 05.12.2023*
*Version: 1.01*


# Newton-Raphson Method

The Newton-Raphson method, often simply called the Newton method, is an iterative numerical technique used for finding the roots (or zeros) of a real-valued function. It is particularly useful for solving nonlinear equations and systems of equations.

**Basic Idea**

The method is based on linearizing the function at each iteration by using its tangent line. The iterative formula for finding the next approximation \(x_{n+1}\) is given by:

\[x_{n+1} =  x_n - \frac{f(x_n)}{f'(x_n)}\]
\[x_{n+1} -  x_n = \Delta x_n = \frac{f(x_n)}{f'(x_n)}\]

Here, \(f(x)\) is the function whose root is being sought, \(f'(x)\) is its derivative, and \(x_n\) is the current approximation.

**Iterative Process**

1. **Initial Guess**: Start with an initial guess \(x_0\).
2. **Iteration**: For each iteration \(n\), compute the next approximation using the formula \(x_{n+1} = x_n - \frac{f(x_n)}{f'(x_n)}\).
3. **Convergence Check**: Repeat the process until the difference between consecutive approximations is smaller than a predefined threshold  \(\epsilon\) or until a maximum number of iterations is reached: \(|\Delta x| \leq \epsilon\)

**Convergence**

The Newton-Raphson method has quadratic convergence, meaning that with each iteration, the number of correct digits roughly doubles. However, convergence is not guaranteed for all functions or initial guesses. It may fail to converge or converge to a local minimum or maximum instead of a root.The Newton-Raphson method is widely used in various fields, including physics, engineering, and numerical analysis. In power system analysis, it is commonly employed to solve the nonlinear power flow equations to determine the steady-state operating conditions of an electrical network.


**Pros and Cons**

- **Efficiency**: (+) The method often converges rapidly for well-behaved functions.
- **Local Convergence**: (+) It converges quickly to a solution near the initial guess.
- **Sensitivity to Initial Guess**: (-) The method's success can depend on the choice of the initial guess.
- **Convergence Issues**: (-) It may fail to converge or converge to an incorrect solution for certain functions.

**Example**

Consider the function \(f(x) = x^2 - 4\) with the goal of finding its root. The derivative is \(f'(x) = 2x\). Applying the Newton-Raphson formula, the iteration becomes \(x_{n+1} = \frac{1}{2}(x_n + \frac{4}{x_n})\).

This method is a fundamental tool in numerical analysis and plays a crucial role in solving nonlinear problems. Its implementation often involves careful consideration of convergence criteria and handling special cases to ensure robustness.

For more in-depth information, including proofs and detailed mathematical concepts, you can refer to the [Wikipedia page on the Newton-Raphson method](https://en.wikipedia.org/wiki/Newton%27s_method).

<div style="page-break-after: always;"></div>

## Application for Power Flow Analysis

In power flow analysis, complex power \(S_i\) is specified at each node \(i\) as a setpoint (demanded power). The injected power has a positive sign, while the absorbed power (loads) is counted as negative (generator reference system). The power equation (system) can be expressed as follows:

\[ \begin{bmatrix} \underline{S}_1 \\ .\\.\\. \\ \underline{S}_N \end{bmatrix} = \begin{bmatrix} \underline{V}_{1} & & & \\ &.\\ & &. \\ & & \\ & && \underline{V}_{N} \end{bmatrix} \cdot \begin{bmatrix} \underline{I}_1 \\ .\\.\\ .\\\underline{I}_N \end{bmatrix}^*\]

where:

\[\underline{V}_i = V_i \cdot e^{j*\phi_i}\]

or simply:

\[\mathbf{s}= diag(\mathbf{v}) \cdot \mathbf{i}^*\]

The currents \(i^*\) are typically not explicitly specified but are measured solely in terms of their magnitudes. As a result, the flows \(P\) and \(Q\) are provided and can be measured in the field. Consequently, the equation requires further rearrangement. By utilizing the node admittance matrix \(Y\) (refer to: [Ybus](ybus.md)), we obtain:
\[\mathbf{s}= diag(\mathbf{v})\cdot (\mathbf{Y} \cdot  \mathbf{v})^*\]

> **Note**: The factor \(3\) is omitted if the bus voltage (\(V\)) is expressed as the phase voltage.

By solving this equation and subtracting the (given!) demanded power (\(S_{D}\)), we can define the power flow problem as follows:

\[f(x)= 0 = s_D - s(x) =  \Delta s \]

This relationship can be divided into real and imaginary parts:
\[ \begin{bmatrix} \Delta p \\ \Delta q \end{bmatrix} = \begin{bmatrix} \text{Re}\{\Delta s\} \\ \text{Im}\{\Delta s\} \end{bmatrix} \]


The mathematical challenge in power flow analysis is the dependence of node power on voltage, resulting in a nonlinear relationship.

The power balance at a node is calculated as follows:

\[S_i = V_i^2 \cdot \underline{y}_{ii}^* + \underline{V}_i^* \sum_{k \in \mathcal{N}} (\underline{y}_{ik} \cdot \underline{V}_{k})^*\]

where:

\[ \underline{V}_i = V_i \cdot e^{j*\phi}\]

The Newton-Raphson method becomes applicable when the solution approach is expanded to systems of equations or matrices:
\[ X_{n+1} = X_n - \frac{F(X_n)}{F'(X_n)} \]
Here, \(X_n\) is the current estimate (is the vector of unknown powers at the nodes), \(F(X_n)\) is the vector given by the node balance equations, and \(F'(X_n)\) is the derivative of \(F(X_n)\) the so called Jacobian matrix \(J\) (refer: [Jacobi](jacobi.md) ).

**Problem Formulation:**
In compact matrix notation, the problem is formulated as follows:


\[ F \cdot \begin{bmatrix} \Delta \varphi \\ \frac{\Delta V}{\text{diag}(V)} \end{bmatrix} = \begin{bmatrix} \Delta p \\ \Delta q \end{bmatrix} \]



> **Note**:
It can be demonstrated that the introduction of \( \Delta V/\text{diag}(V) \) avoids the computation of cosine and sine functions, thus accelerating the calculation of very large matrices.

**Goal**
The goal is to find a solution for the angle \(\varphi\) and the voltage \(V\) for each bus. Once the values of \(\varphi\) and \(V\) are determined, all branch flows can be subsequently calculated in a post-processing step (refer: [Losses](losses.md)).


**Calculating the Right Side**
The right side of the equation is determined as follows:

For \(i\)-th row of \(\Delta p\):
   \[ \Delta p_i = P_{D_i} - P_i \]
   where \(P_{D_i}\) is the demanded active power at node \(i\).

For \(i\)-th row of \(\Delta q\):
   \[ \Delta q_i = Q_{D_i} - Q_i \]
   where \(Q_{D_i}\) is the demanded reactive power at node \(i\).

These components on the right side represent the difference between the demanded and actual active/reactive powers at each node, forming the vector \(\begin{bmatrix} \Delta p \\ \Delta q \end{bmatrix}\) in the Newton-Raphson power flow equation.


The active power \(P_i\) and reactive power \(Q_i\) for node \(i\) can be calculated as follows:

The complex power \(S_i\) injected at node \(i\) is given by:
\[ S_i = V_i \sum_{k=1}^{n} (G_{ik}\cos(\varphi_i - \varphi_k) + B_{ik}\sin(\varphi_i - \varphi_k)) + jV_i \sum_{k=1}^{n} (G_{ik}\sin(\varphi_i - \varphi_k) - B_{ik}\cos(\varphi_i - \varphi_k)) \]

where:
- \(V_i\) is the voltage magnitude at node \(i\),
- \(\phi_i\) is the voltage angle at node \(i\),
- \(G_{ik}\) is the conductance (real part) of the admittance between nodes \(i\) and \(k\),
- \(B_{ik}\) is the susceptance (imaginary part) of the admittance between nodes \(i\) and \(k\),
- \(n\) is the total number of neighboring nodes.

Now, the active power \(P_i\) and reactive power \(Q_i\) can be obtained by separating the real and imaginary parts of \(S_i\):
\[ P_i = V_i \sum_{k=1}^{n} (G_{ik}\cos(\phi_i - \phi_k) + B_{ik}\sin(\phi_i - \phi_k)) \]
\[ Q_i = V_i \sum_{k=1}^{n} (G_{ik}\sin(\phi_i - \phi_k) - B_{ik}\cos(\phi_i - \phi_k)) \]

These expressions represent the active and reactive power contributions from neighboring nodes to the total power injection at node \(i\).


**Calculating the Left Side**
Let's consider the left side of the problem, \(F \cdot \begin{bmatrix} \Delta \varphi \\ \frac{\Delta V}{\text{diag}(V)} \end{bmatrix}\), where \(F\) is the functional matrix or Jacobian matrix.

The functional matrix \(F\) is derived by taking partial derivatives of the power injection equations with respect to voltage angles (\(\varphi\)) and magnitudes (\(V\)). Let's denote \(F_{ij}\) as the element in the \(i\)-th row and \(j\)-th column of the Jacobian matrix.

The real power injection equation for node \(i\) (\(P_i\)) with respect to voltage angles and magnitudes is given by:
\[ \frac{\partial P_i}{\partial \varphi_i} = -Q_i - V_i^2 B_{ii} - \sum_{k=1}^{n} V_i V_k (G_{ik} \sin(\varphi_i - \varphi_k) - B_{ik} \cos(\varphi_i - \varphi_k)) \]
\[ \frac{\partial P_i}{\partial V_i} = 2V_i G_{ii} + \sum_{k=1}^{n} V_k (G_{ik} \cos(\varphi_i - \varphi_k) + B_{ik} \sin(\varphi_i - \varphi_k)) \]

<div style="page-break-after: always;"></div>

The reactive power injection equation for node \(i\) (\(Q_i\)) with respect to voltage angles and magnitudes is given by:
\[ \frac{\partial Q_i}{\partial \varphi_i} = P_i - V_i^2 G_{ii} - \sum_{k=1}^{n} V_i V_k (G_{ik} \cos(\varphi_i - \varphi_k) + B_{ik} \sin(\varphi_i - \varphi_k)) \]
\[ \frac{\partial Q_i}{\partial V_i} = -2V_i B_{ii} + \sum_{k=1}^{n} V_k (G_{ik} \sin(\varphi_i - \varphi_k) - B_{ik} \cos(\varphi_i - \varphi_k)) \]

This process is repeated for each node \(i\) in the power system. The resulting \(F\) matrix is then used in the Newton-Raphson iterative process to solve the power flow equations.

Refer [Jacobi](jacobi.md) for detailed informations.

**Solving the Equations of the form \(A \mathbf{x} = \mathbf{b}\)**
If we express the problem in the form \(A \mathbf{x} = \mathbf{b}\), where \(\mathbf{x}\) is the vector we seek (in this case, \(\begin{bmatrix} \Delta \phi \\ \frac{\Delta V}{\text{diag}(V)} \end{bmatrix}\)), then the Jacobian matrix \(J\) can be considered as matrix \(A\) and the vector on the right side of the equation as \(\mathbf{b}\).

The Gauss elimination method or LU decomposition can be applied to solve this system of equations efficiently. The LU decomposition is often preferred as it decomposes the matrix into lower and upper triangular matrices, allowing for more efficient solving of systems of equations.

> **Note**:
In Julia, one could use the solver function such as `\` for LU decomposition to solve linear systems of equations.

**Inital Values**
In practice, the initial values for the voltage \(V\) and angle \(\varphi\) were chosen as follows:
\(V\) = 1.0
\(\varphi\) = 0.0

# `calcNewtonRaphson!`

**Description:**
Main function for calculating the Newton-Raphson power flow.

**Syntax:**
```
function calcNewtonRaphson!(Y::AbstractMatrix{ComplexF64}, nodes::Vector{ResDataTypes.Node}, Sbase_MVA::Float64, maxIte::Int, tolerance::Float64 = 1e-6, verbose::Int = 0, sparse::Bool = false)
```

**Parameters:**
- `Y`: Admittance matrix of the power system.
- `nodes`: Vector of nodes representing the power system.
- `Sbase_MVA`: Base power value in MVA.
- `maxIte`: Maximum number of iterations allowed.
- `tolerance`: Convergence tolerance. Defaults to `1e-6`.
- `verbose`: Verbosity level for printing debug information. Defaults to `0`.
- `sparse`: Use sparse matrix representation. Defaults to `false`.

**Returns:**
- Tuple containing the number of iterations and the convergence result.
  - `iteration_count`: Number of iterations performed.
  - `erg`: Convergence result (0: Convergence reached, 1: No convergence, 2: Unsolvable system of equations, 3: Error).

<div style="page-break-after: always;"></div>

**Usage:**
```
iteration_count, erg = calcNewtonRaphson!(Y, nodes, Sbase_MVA, maxIte, tolerance, verbose, sparse)
```
**Algorithm Steps:**

1. **Initialization:**
   - Extract bus data, determine adjacent branches, and initialize counters.
   - Print system information if `verbose` is greater than 0.

2. **Main Iteration Loop:**
   - Perform iterations until convergence or reaching the maximum iteration limit.
     - Calculate power feeds and residuals.
     - Check for convergence by comparing the norm of residuals with the specified tolerance.

3. **Convergence Check:**
   - If convergence is reached, print a convergence message and break out of the loop.

4. **Jacobian Calculation:**
   - Calculate the Jacobian matrix for the system of equations.

5. **Linear System Solution:**
   - Attempt to solve the linear system of equations using the calculated Jacobian and residuals.
   - Handle unsolvable system errors.

6. **Update Voltage Values:**
   - Update the voltage values of nodes based on the solution.
   - Handle errors during the update process.

7. **Print Solution:**
   - Print final node voltage values if `verbose` is greater than 1.

8. **Output:**
   - Return the number of iterations and convergence result.

**Example:**

```
Y = ...
nodes = ...
Sbase_MVA = ...
maxIte = 8
tolerance = 1e-8
verbose = 1

iteration_count, convergence_result = calcNewtonRaphson!(Y, nodes, Sbase_MVA, maxIte, tolerance, verbose)
```
*refer: `testparser.jl`
**Note:**
- The function performs the Newton-Raphson power flow calculation, updating the nodes with the final voltage values.
- Use appropriate values for `maxIte` and `tolerance` based on the specific power system and convergence requirements.
- Set `verbose` to control the amount of debug information printed.


