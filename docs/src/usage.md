# Getting Started Calculation of AC Power Flow
## Quick Start
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