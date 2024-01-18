# Getting Started Calculation of AC Power Flow
# Installation
## First-time Julia Users
SPARLECTRA requires the Julia programming language. You can download Julia from [here](https://julialang.org/downloads/).


## Installing Program Files
Press ALT GR + ] in the Julia REPL to go to the Package Manager. To add and register the Sparlectra package, execute the following command in the Julia Package Manager (REPL):

```
] add https://github.com/Welthulk/Sparlectra.jl
```

The Julia package manager will automatically download the required packages.

>**Hint**:
For own development you need to clone the package from the Git repository and place it in the user directory or home directory (e.g., C:\Users\USERNAME\.julia or %HOME%/.julia) under the `dev` directory.

## Quick Start Example 1 using VS-Code
Open VS-Code, create an Julia-File (e.g. `quickstart.jl`) and place the following code in it:
```julia
using Sparlectra.SparlectraNet
using Sparlectra.SparlectraResult
using BenchmarkTools

# Specify the power system case file
case = "bsp7.json"

# Create the network from the specified case file
myNet = createNetFromFile(case)

# Create the YBUS matrix
Y = createYBUS(myNet.branchVec, myNet.shuntVec, true, true)

# Measure the execution time of the following block
etim = @elapsed begin
    # Run Newton-Raphson load flow analysis
    ite, erg = calcNewtonRaphson!(Y, myNet.nodeVec, myNet.baseMVA, 8, 1e-6, 1, true)
end

# Check if the Newton-Raphson analysis converged successfully
if erg == 0
    # Calculate losses and print results
    calcNetLosses!(myNet.nodeVec, myNet.branchVec, myNet.baseMVA, false)
    printACPFlowResults(myNet, etim, ite)
end
```
Before execution, copy the file `bsp7.json` from the `test` directory (`<%HOME%>/.julia/packages/Sparlectra/<xymz>/test`) to the same directory as the `quickstart.jl` file. Then, execute the code in VS Code. The output should look like this:
```bash
create network from file: bsp7.json
        1:              2:              3:              4:              5:
1:      22.08∠-62.84    -11.04∠-62.85   -11.04∠-62.85   0.0             0.0
2:      -11.04∠-62.85   22.08∠-62.84    0.0             -11.04∠-62.85   0.0
3:      -11.04∠-62.85   0.0             22.08∠-62.84    -11.04∠-62.85   0.0
4:      0.0             -11.04∠-62.85   -11.04∠-62.85   33.12∠-62.84    -11.04∠-62.85
5:      0.0             0.0             0.0             -11.04∠-62.85   11.04∠-62.84

∑Load: [3.0, 6.0], ∑Gen [4.1, 2.0]  Δp, Δq: [1.1, -4.0]

calcNewtonRaphson: Matrix-Size: 8 x 8, Sbase_MVA: 100.0, total number of nodes: 4, number of PQ nodes: 4, number of PV nodes: 0
norm 2.659656e-02, tol 1.000000e-06, ite 0
norm 4.258203e-05, tol 1.000000e-06, ite 1
norm 3.611743e-10, tol 1.000000e-06, ite 2
Convergence is reached after 2 iterations.
================================================================================
| SPARLECTRA Version 00.04.03 - AC Power Flow Results                          |
================================================================================
Date         :  18-Jan-24 17:54:30
Iterations   :         2
Converged in :  0.008790 seconds
Case         :           bsp7
BaseMVA      :       100
Nodes        :         5 (PV: 0 PQ: 4 Slack: 1)
Branches     :         5
Lines        :         5
Trafos       :         0
Generators   :         2
Loads        :         3
Shunts       :         0

total losses (I^2*Z): P =      0.003 [MW], Q =      0.005 [MVar]

======================================================================================================
| Bus                  | Vn [kV]    | V [pu]     | phi [deg]  | P [MW]     | Q [MVar]   | Type       |
======================================================================================================
| Bus_1                | 110        | 0.999      | -0.127     | 0.100      | -0.000     | PQ         |
| Bus_2                | 110        | 0.998      | -0.118     | -1.000     | -2.000     | PQ         |
| Bus_3                | 110        | 0.998      | -0.118     | -1.000     | -2.000     | PQ         |
| Bus_4                | 110        | 0.999      | -0.090     | -0.000     | 0.000      | PQ         |
| Bus_5                | 110        | 1.000      | 0.000      | 1.903      | -0.522     | Slack      |
------------------------------------------------------------------------------------------------------

===========================================================================================================================
| Branch                    | From  | To    | P [MW]     | Q [MVar]   | P [MW]     | Q [MVar]   | Pv [MW]    | Qv [MVar]  |
===========================================================================================================================
| L1_2                      | 1     | 2     | 0.050      | -0.000     | -0.050     | -0.905     | 0.000      |  0.000     |
| L1_3                      | 1     | 3     | 0.050      | -0.000     | -0.050     | -0.905     | 0.000      |  0.000     |
| L2_4                      | 2     | 4     | -0.950     | -1.095     | 0.951      | 0.191      | 0.001      |  0.001     |
| L3_4                      | 3     | 4     | -0.950     | -1.095     | 0.951      | 0.191      | 0.001      |  0.001     |
| L4_5                      | 4     | 5     | -1.901     | -0.382     | 1.903      | -0.522     | 0.001      |  0.003     |
---------------------------------------------------------------------------------------------------------------------------
```

## Quick Start Example 2
For a practical example, refer to the the file [`testparser.jl`](../../test/testparser.jl). Open this file in Visual Studio Code and execute the code to observe how the `acpflow` function is implemented. To use it from the command line use the Windows Batch-File `runacpflow.bat`.

## Creating a Network
1. **Network Data Preparation**: Prepare your network data in JSON format (see [network_configuration_guide.md](network_configuration_guide.md)). Ensure that the file includes essential information about nodes and branches.

2. **Setting Up the Function Call**: providing the necessary parameters:
   - `sparse`: Boolean value for sparse matrix calculations.
   - `casefile::String`: Filepath to your network dataset in JSON format.
   - `writeCase::Bool`: Decide whether to export your networw to an Matpower case file.
   - `iterations::Int`: Set the number of iterations for the Newton-Raphson algorithm, good values are between 3 and 8.
   - `verbose::Int`: Adjust the verbosity level for detailed output.
   - `tol::Float`: tolerance for Newton-Raphson algorithm convergence, good values are between 1e-6 and 1e-9.
3. **Create Network**: Create the network object using the `createNetFromFile` function.
```js
myNet = createNetFromFile(jpath,base_MVA, (verbose>2))  
```
4. **Create Y-Bus Matrix**: Create the Y-Bus matrix of the network using the `createYbus` function.
```js
Y = createYBUS(myNet.branchVec, myNet.shuntVec, (verbose>2), sparse)
```
5. **Run Load Flow Calculation**: Run the load flow calculation using the `newtonRaphson` function.
```js
ite,erg = calcNewtonRaphson!(Y, myNet.nodeVec, myNet.baseMVA, iterations, tol, verbose, sparse)      
```
6. **Calc network losses**: After running the function, call function `calcNetLosses!` losses in the network 
```js
calcNetLosses!(myNet.nodeVec, myNet.branchVec, myNet.baseMVA, log)    
```
7. **Interpreting the Results**: Examine the results. The function `printACPFlowResults` will output information such as losses in the network and execution time.
```js
printACPFlowResults(myNet, etime, ite)
```

## Example Output
```
Convergence is reached after 3 iterations.
================================================================================
| SPARLECTRA Version 00.04.01 - AC Power Flow Results                          |
================================================================================
Date         :  20-Dec-23 21:50:20
Iterations   :         3
Converged in :  0.03365 seconds
Case         :           bsp3
BaseMVA      :       100
Nodes        :         3 (PV: 1 PQ: 1 Slack: 1)
Branches     :         3
Lines        :         3
Trafos       :         0
Generators   :         2
Loads        :         1
Shunts       :         0

total losses (I^2*Z): P =      0.000 [MW], Q =      4.761 [MVar]

======================================================================================================
| Bus                  | Vn [kV]    | V [pu]     | phi [deg]  | P [MW]     | Q [MVar]   | Type       |
======================================================================================================
| ASTADT               | 110        | 1.012      | -1.983     | -100.000   | -30.000    | PQ         |
| STADION1             | 110        | 1.030      | 0.596      | 70.000     | 34.898     | PV         |
| VERBUND              | 110        | 1.020      | 0.000      | 30.000     | -2.974     | Slack      |
------------------------------------------------------------------------------------------------------

===========================================================================================================================
| Branch                    | From  | To    | P [MW]     | Q [MVar]   | P [MW]     | Q [MVar]   | Pv [MW]    | Qv [MVar]  |
===========================================================================================================================
| L1_2                      | 1     | 2     | -56.771    | -20.860    | 56.771     | 22.848     | 0.000      |  2.935     |
| L2_3                      | 2     | 3     | 13.229     | 12.050     | -13.229    | -12.745    | 0.000      |  0.259     |
| L3_1                      | 3     | 1     | 43.229     | 9.771      | -43.229    | -9.140     | 0.000      |  1.568     |
---------------------------------------------------------------------------------------------------------------------------
```