# Introduction
The `calcNetLosses!` function in Julia is designed to calculate network losses using the \(I^2 \cdot Z \)
 method. The function takes as input the network topology represented by a vector of nodes (`nodes`) and a vector of branches (`branchVec`). Additionally, the base MVA (`Sbase_MVA`) is provided for normalization (pu). The function calculates branch power flows, accumulates losses, and updates node power information.

## Loss Calculation Formula
The complex losses (\(S_{\text{ij}}\)) for each branch are calculated in this project using the formula:

\[ S_{ij} = P_{ij} + jQ_{ij} = vi \cdot \exp(j \cdot \phi_i) \cdot \left( \left(vi \cdot \exp(j \cdot \phi_i) - vk \cdot \exp(j \cdot \phi_k) \cdot Y_{ik} + vi \cdot \exp(j \cdot \phi_i) \cdot Y_{0ik}\right) \right) \]

where \( Y_{0ik} \) is the transconductance.


## function calcNetLosses!
The `calcNetLosses!` function is implemented in this file to compute network losses using the I^2*Z method. It iterates through the branches, calculates power flows, and updates the power information for nodes. The detailed implementation and calculations are done within the function. For further details, refer to the function documentation and the comments within the code.

## Function Parameters
- `nodes`: Vector of nodes representing the electrical network.
- `branchVec`: Vector of branches representing the electrical network.
- `Sbase_MVA`: Base MVA for normalization of power values.
- `log`: Optional Boolean parameter (default is `false`) for enabling debug logging.
