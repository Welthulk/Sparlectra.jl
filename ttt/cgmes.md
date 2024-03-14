# Prerequisites
<p>An installation (or HTTP connection) of the Apache Jena Fuseki Graph database is required. For additional information see here: <br>
https://jena.apache.org/documentation/fuseki2/fuseki-quick-start.html<br>

## Usage:
Create a dataset named "explore2" in Jena Fuseki (the name of the dataset can be changed in the source code) and import the following RDF files from the directory ../sparlectra/data/MiniGrid:
- MiniGridTestConfiguration_BC_EQ_v3.0.0.xml
- MiniGridTestConfiguration_BC_SSH_v3.0.0.xml
- MiniGridTestConfiguration_BC_SV_v3.0.0.xml
- MiniGridTestConfiguration_BC_TP_v3.0.0.xml

<p>Hint:<br>
To check the installation you can enter the following command in the query explorer:<br>

```sql
SELECT ?subject ?predicate ?object
WHERE {
  ?subject ?predicate ?object
}
LIMIT 25
```


# Example
Go to the directory ../test and run file:
```
testsparql.jl
```
The following output should appear:

```Script
[ Info: createNetFromTripleStore: slackBusName - HG2, slackBusIdx: 10
[ Info: fix wrong sequence number for node: 1_Node-1 trafo: T4 from: Seite3 to: Seite1
[ Info: fix wrong sequence number for node: 4_Node-11 trafo: T1 from: Seite2 to: Seite1
[ Info: fix wrong sequence number for node: H_Node-4 trafo: T4 from: Seite1 to: Seite3
[ Info: fix wrong sequence number for node: HG1_Node-10 trafo: T1 from: Seite1 to: Seite2
  7.435554 seconds (11.79 M allocations: 786.351 MiB, 3.90% gc time, 99.15% compilation time)
convertion to Matpower CASE-Files, Testcase: (MiniGrid), Filename: (C:\Users\scud\.julia\dev\Sparlectra\data\MiniGrid\MiniGrid.m)
Convergence is reached after 3 iterations.
================================================================================
| SPARLECTRA Version 00.31.19 - AC Power Flow Results                          |
================================================================================
Date         :  30-Nov-23 18:31:13
Iterations   :         3
Converged in :  3.814679 seconds
Case         :       MiniGrid
BaseMVA      :       350
Nodes        :        13 (PV: 2 PQ: 10 (Aux: 2) Slack: 1
Branches     :        17
Lines        :         8
Trafos       :         6
Generators   :         5
Loads        :         3
Shunts       :         0

total losses (I^2*Z): P =      0.189 [MW], Q =      2.086 [MVar]

======================================================================================================
| Bus                  | Vn [kV]    | V [pu]     | phi [deg]  | P [MW]     | Q [MVar]   | Type       |
======================================================================================================
| 1_Node-1             | 380        | 0.971      | 1.552      | -0.000     | 0.000      | PQ         |
| 3_Node-5             | 110        | 0.968      | 1.438      | 0.000      | 0.000      | PQ         |
| 4_Node-11            | 110        | 0.967      | 1.588      | -0.000     | 0.000      | PQ         |
| 2_Node-3             | 110        | 0.971      | 1.552      | -1.808     | -1.177     | PQ         |
| 5_Node-6             | 110        | 0.971      | 1.566      | -4.631     | -6.219     | PQ         |
| 8_Node-2             | 30         | 0.971      | 1.552      | -4.631     | -6.219     | PQ         |
| H_Node-4             | 30         | 0.971      | 1.552      | 0.000      | 0.000      | PQ         |
| HG1_Node-10          | 21         | 1.000      | 1.911      | 5.000      | -10.841    | PV         |
| 7_Node-8             | 10         | 1.012      | 3.284      | 9.000      | 5.000      | PQ         |
| HG2_Node-9           | 10         | 1.000      | 0.000      | -17.829    | -3.241     | Slack      |
| 6_Node-7             | 10         | 1.000      | 3.078      | 10.458     | 18.564     | PV         |
| aux_#12              | 380        | 0.971      | 1.552      | 0.000      | -0.000     | PQ         |
| aux_#13              | 380        | 0.971      | 1.552      | -0.000     | 0.000      | PQ         |
------------------------------------------------------------------------------------------------------

===========================================================================================================================
| Branch                    | From  | To    | P [MW]     | Q [MVar]   | P [MW]     | Q [MVar]   | Pv [MW]    | Qv [MVar]  |
===========================================================================================================================
| L2                        | 2     | 3     | -6.252     | 4.467      | 6.258      | -4.447     | 0.006      |  0.020     |
| L1                        | 4     | 2     | 3.616      | 2.355      | -3.612     | -2.342     | 0.004      |  0.013     |
| L3_a                      | 4     | 5     | -3.616     | -2.355     | 3.617      | 2.356      | 0.000      |  0.002     |
| L3_b                      | 4     | 5     | -1.808     | -1.177     | 1.808      | 1.178      | 0.000      |  0.001     |
| L4                        | 5     | 2     | 7.992      | 5.868      | -7.984     | -5.835     | 0.008      |  0.033     |
| L5                        | 5     | 3     | 1.270      | 6.569      | -1.263     | -6.546     | 0.007      |  0.023     |
| L6                        | 11    | 9     | -8.915     | -4.911     | 9.000      | 5.000      | 0.085      |  0.089     |
| aux_#12                   | 1     | 12    | -0.000     | 0.000      | 0.000      | -0.000     | 0.000      |  0.000     |
| aux_#12                   | 12    | 4     | 0.000      | 0.000      | 0.000      | 0.000      | 0.000      |  0.000     |
| aux_#12                   | 12    | 6     | 0.000      | 0.000      | 0.000      | 0.000      | 0.000      |  0.000     |
| aux_#13                   | 1     | 13    | 0.000      | 0.000      | -0.000     | -0.000     | 0.000      |  0.000     |
| aux_#13                   | 13    | 4     | -0.000     | 0.000      | 0.000      | -0.000     | 0.000      |  0.000     |
| aux_#13                   | 13    | 7     | -0.000     | -0.000     | 0.000      | 0.000      | 0.000      |  0.000     |
| T2                        | 2     | 10    | 17.849     | 3.710      | -17.829    | -3.241     | 0.020      |  0.469     |
| T1                        | 3     | 8     | -4.995     | 10.993     | 5.000      | -10.841    | 0.005      |  0.152     |
| T5                        | 5     | 11    | -12.879    | -14.794    | 12.915     | 15.650     | 0.036      |  0.856     |
| T6                        | 5     | 11    | -6.440     | -7.397     | 6.458      | 7.825      | 0.018      |  0.428     |
---------------------------------------------------------------------------------------------------------------------------
```