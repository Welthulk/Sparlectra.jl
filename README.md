# SPARLECTRA
<a href="https://github.com/Welthulk/Sparlectra.jl/tree/main/"><img align="left" width="100" src="docs/src/png/logo.png" style="margin-right: 20px" /></a>

This package contains tools for subsequent network calculations. It primarily features a program for calculating load flow using the Newton-Raphson method. The focus is to provide valuable insights into load flow calculations for both students and ambitious professionals.


---

## Installation
```julia
using Pkg
Pkg.add("Sparlectra")
```


### Network Creation
This package supports the import and export of Matpower .m files, although currently it only reads bus, generator, and branch data from these files. Please note that additional Matlab functions within the .m file are not supported. Additionally, you can create your own network using easy-to-use functions provided by the package.



### License
This project is licensed under the BSD-3-Clause - [The license file](LICENSE) contains the complete licensing information.









