# Getting Started 

The `acpflow` function (an example is provided in the file [`testparser.jl`](../../test/testparser.jl).) is designed to perform a load flow calculation. If you're looking to program your own load flow calculations, here's a guide on using this function.

## calculation of acpflow

To use the `acpflow` function for load flow calculations, follow these steps:

1. **Network Data Preparation**: Prepare your network data in JSON format. Ensure that the file includes essential information about nodes, branches, and shunts.

2. **Setting Up the Function Call**: Use the `acpflow` function, providing the necessary parameters:
   - `casefile::String`: Filepath to your network dataset in JSON format.
   - `writeCase::Bool`: Decide whether to write the Matpower case file.
   - `iterations::Int`: Set the number of iterations for the Newton-Raphson algorithm.
   - `verbose::Int`: Adjust the verbosity level for detailed output.

```julia
@time acpflow("your_network_data.json", true, 10, 1)
```

3. **Interpreting the Results**: After running the function, examine the results. The function will output information such as losses in the network and execution time.

## Example Usage

For a practical example, refer to the `testparser.jl` file located in the `test` directory. Open this file in Visual Studio Code and execute the code to observe how the `acpflow` function is implemented.

## Understanding the Code

The `acpflow` function internally performs the following steps:

1. **Load and Create Network**: Reads the network data from the specified JSON file and creates the network object.

2. **Matpower Export (Optional)**: Exports the network to the Matpower file format if `writeCase` is set to `true`.

3. **Y-Bus Matrix Generation**: Creates the Y-Bus matrix of the network.

4. **Newton-Raphson Calculation**: Executes the Newton-Raphson load flow calculation.

5. **Result Output**: Displays the calculated results, including losses in the network.

## Troubleshooting

- If the file is not found, an error message will be displayed.
- A warning will be issued if the Newton-Raphson algorithm fails to converge.
- An error message will be shown in case of any issues during the calculation.

## Customization Options

Adjust the following options to customize your load flow calculation:
- `sparse`: Boolean value for sparse matrix calculations.
- `tol`: Tolerance for Newton-Raphson algorithm convergence.
- `base_MVA`: Base MVA for the network data.
- `log`: Enable or disable logging.

By following these guidelines, you can effectively utilize the `acpflow` function for your load flow calculations.

