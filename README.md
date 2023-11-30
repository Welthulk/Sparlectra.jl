# SPARLECTRA

This package contains tools for subsequent network calculations. It primarily features a program for calculating load flow using the Newton-Raphson method. The program has been developed through self-directed learning and is designed to provide valuable insights into load flow calculations for both students and ambitious professionals. 



## About the Program

The load flow calculation will be expanded in the future to include additional methods. My goal is to create a comprehensive platform that not only caters to students but also provides insights and resources for advanced load flow calculations to experienced professionals.

## Network Data Files

This package supports the use of various network data file formats, allowing users to leverage different sources for load flow calculations. Below are the supported formats:

1. **CGMES Version 2.4.15:**
   The network data files should be loaded into a Jena Fuseki Graph Database. The CGMES data format is partially supported, with testing conducted on the ENTSOE example MiniGrid.

2. **Matpower .m-Files:**
   This package also supports Matpower .m-files. Currently, only bus, generator, and branch data are read from these files. Note that additional Matlab functions within the .m file are not supported.

3. **Proprietary JSON Format:**
   There is a proprietary network data format in JSON. Users can leverage this format for their load flow calculations.

It's important to note that while CGMES Version 2.4.15 is partially supported, testing has been conducted primarily on the ENTSOE example MiniGrid. Additionally, Matpower .m-files are supported for bus, generator, and branch data, with limitations on additional Matlab functions within the files. Users are encouraged to review the documentation for specific details on supported features and limitations.


## Contribution Guidelines

This project encourages and welcomes contributions. If you have improvements or new features in mind, feel free to submit a pull request. Collaborative efforts are key to enhancing the functionality and robustness of the tool.

### Naming Conventions
The project follows the Julia Naming Conventions for the most part, but it's important to note that the naming convention for functions might deviate. In this module, functions are written in CamelCase with a lowercase initial letter.

Contributors are welcome to adhere to the Julia Naming Conventions for functions or choose to follow the CamelCase convention used in this module; the choice is yours. Consistency within the module is encouraged, but flexibility is provided to accommodate different preferences.

### How to Contribute

1. Fork the repository.
2. Create a new branch for your feature or bug fix.
3. Make your changes and ensure they adhere to the chosen naming convention.
4. Test your changes thoroughly.
5. Submit a pull request.

### Project Ownership

Sparlectra.jl is a project developed by a hobbyist, and it is not maintained or supported by any organization. As such, contributions are voluntary, and there is no official organization overseeing the project.

### Code Style

- Use consistent indentation.
- Follow Julia best practices where possible.
- Document your code for better understanding.

### Network Data

While contributions to the project are appreciated, please note that providing support for individualized network data issues is beyond the scope of this project, as it is not maintained by an organization. Users are encouraged to take initiative in resolving such issues independently and sharing their results with the community.

Your understanding and collaboration are crucial in maintaining a positive and open-source development environment.


## License
This project is licensed under the BSD-3-Clause - [The license file](LICENSE) contains the complete licensing information.

## Installation
Installation guide can be found in [development](/docs/src/development.md)

## Additional information
Additional information can be found in the [documentation](sparlectra/docs/src/)








