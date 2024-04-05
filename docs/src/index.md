Sparlectra
=============

Sparlectra is a Julia package for the simulation of electrical power systems. It primarily features a program for calculating load flow using the Newton-Raphson method. The focus is to provide valuable insights into load flow calculations for both students and ambitious professionals. The package supports the import and export of Matpower .m files, although currently it only reads bus, generator, and branch data from these files. Please note that additional Matlab functions within the .m file are not supported. Additionally, you can create your own network using easy-to-use functions provided by the package.

---

#### Installation
For installation, run the following command in the Julia REPL:
```julia
import Pkg
Pkg.add("Sparlectra")
```
---

#### Contributors
- [Udo Schmitz](https://www.linkedin.com/in/udo-schmitz-315536250/) - Electrical Engineer