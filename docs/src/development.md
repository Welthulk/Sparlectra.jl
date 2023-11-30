# Installation
## First-time Julia Users
SPARLECTRA requires the Julia programming language. You can download Julia from [here](https://julialang.org/downloads/).

## Placing Program Files
Since SPARLECTRA is not registered yet, executing the following command in the package manager won't work:

```
] add Sparlectra
```

Therefore, you need to clone the package from the Git repository and place it in the user directory or home directory (e.g., C:\Users\USERNAME\.julia or %HOME%/.julia) under the `dev` directory.

### Installing Required Julia Packages for Windows Users
The missing Julia packages will be automatically downloaded by the Julia package manager during program execution. Alternatively, you can manually install the packages by running the Windows batch script `install.bat`.

## First Run
Windows users can execute the Windows batch script `runtest.bat`. Optionally, you can specify the path to the installed Julia version.

Command:
```bash
runtest.bat
```

Output:
```bash
Convergence is reached after 3 iterations.
Test Summary: | Pass  Total  Time
Sparlectra.jl |    3      3  6.3s
```
>**Note**:
Note that the long execution is due to the runtime compilation of the program.
# Advantage Development
## Local Development (without GIT-Repo)
If you wish to further develop or use SPARLECTRA without Git, you need to add the following lines to the `Sparlectra.jl` file:

```julia
using Pkg
_path = joinpath(pwd(), "..", "Sparlectra")
Pkg.develop(PackageSpec(path=_path))
```

The variable `_path` indicates where the program package is located.

## Development with GIT-Repo in a different directory other than `dev`
If a different location than the `dev` directory is desired, the Git repository must be made available manually. To do this, it needs to be added to the Julia environment.

For development with a Git repo, it is necessary to add the package to your environment. This can be done with the following command. To load the packages automatically, add the following lines to the file `<programm>/etc/julia/startup.jl`:

```julia
import Pkg
Pkg.activate("Sparlectra")
```

## Adding GIT-Repo
```bash
pkg> add https://github.com/Welthulk/Sparlectra.jl:Sparlectra
```