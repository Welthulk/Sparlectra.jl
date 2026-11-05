# Copyright 2023-2026 Udo Schmitz
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# file: tools/install_workshops.jl
# purpose: the workshop mode of tools/install_webui.sh and .bat
#          (`--workshops`): sets up a folder that is a Julia environment with
#          Sparlectra and IJulia, installs a Jupyter kernel for it, copies the
#          workshop notebooks of the installed Sparlectra into it and starts
#          Jupyter. The installers only make sure Julia 1.12 or newer is
#          there and then run this file, so the steps are the same on every
#          platform. Runs without a project (Pkg is a standard library):
#            julia install_workshops.jl [--dir <folder>]
#                  [--rev <branch> | --path <checkout>] [--no-start]

using Pkg

const REPO_URL = "https://github.com/Welthulk/Sparlectra.jl"
const KERNEL_NAME = "Sparlectra workshops"

const USAGE = """
usage: julia install_workshops.jl [options]
  --dir <folder>     workshop folder (default: ~/sparlectra-workshops)
  --rev <branch>     install Sparlectra from this branch, tag or commit of
                     $(REPO_URL) instead of the registered release
  --url <git url>    repository for --rev (default: the one above)
  --path <checkout>  install Sparlectra from a local checkout (Pkg.develop)
  --no-start         set up only, do not start Jupyter
"""

function parse_options(args::Vector{String})
  dir = joinpath(homedir(), "sparlectra-workshops")
  rev = nothing
  url = REPO_URL
  path = nothing
  start = true
  i = 1
  value(i) = i < length(args) ? args[i+1] : (print(USAGE); error("option $(args[i]) needs a value"))
  while i <= length(args)
    a = args[i]
    if a == "--dir"
      dir = value(i); i += 2
    elseif a == "--rev"
      rev = value(i); i += 2
    elseif a == "--url"
      url = value(i); i += 2
    elseif a == "--path"
      path = value(i); i += 2
    elseif a == "--no-start"
      start = false; i += 1
    else
      print(USAGE)
      error("unknown option $(a)")
    end
  end
  rev !== nothing && path !== nothing && error("--rev and --path exclude each other")
  return (dir = abspath(expanduser(dir)), rev = rev, url = url, path = path === nothing ? nothing : abspath(expanduser(path)), start = start)
end

# The notebooks are generated with the generic kernel name "julia". The copy
# names the kernel installed for this folder, so Jupyter and VS Code open it
# without asking. The kernelspec object is flat (no nested braces), which the
# pattern relies on.
function with_kernel(text::String, spec::String, display::String)
  kernel = "\"kernelspec\": {\"name\": \"$(spec)\", \"display_name\": \"$(display)\", \"language\": \"julia\"}"
  return replace(text, r"\"kernelspec\"\s*:\s*\{[^{}]*\}" => kernel)
end

function main(args::Vector{String})
  o = parse_options(args)
  mkpath(o.dir)
  println("Workshop folder: ", o.dir)
  Pkg.activate(o.dir)

  # Sparlectra: removed first when the folder already has it, so a second run
  # takes the source asked for now and, for the registered release, the
  # newest version
  haskey(Pkg.project().dependencies, "Sparlectra") && Pkg.rm("Sparlectra")
  if o.path !== nothing
    println("Installing Sparlectra from ", o.path)
    Pkg.develop(path = o.path)
  elseif o.rev !== nothing
    println("Installing Sparlectra from ", o.url, " at ", o.rev)
    Pkg.add(url = o.url, rev = o.rev)
  else
    println("Installing the registered Sparlectra release")
    Pkg.add("Sparlectra")
  end
  Pkg.add("IJulia")
  Pkg.precompile()

  # Jupyter kernel that runs in this folder's environment, whatever folder
  # a notebook is opened from. IJulia is loaded at run time, after the
  # activation above, so every access goes through invokelatest.
  ijulia(name) = Base.invokelatest(getfield, Base.require(Main, :IJulia), name)
  kdir = Base.invokelatest(ijulia(:installkernel), KERNEL_NAME, "--project=$(o.dir)")
  spec = basename(kdir)
  m = match(r"\"display_name\"\s*:\s*\"([^\"]*)\"", read(joinpath(kdir, "kernel.json"), String))
  display = m === nothing ? KERNEL_NAME : String(m.captures[1])

  # the notebooks of exactly the Sparlectra just installed
  sp = only(d for d in values(Pkg.dependencies()) if d.name == "Sparlectra")
  src = joinpath(sp.source, "notebooks")
  files = sort(filter(f -> startswith(f, "workshop_") && endswith(f, ".ipynb"), readdir(src)))
  isempty(files) && error("no workshop notebooks in $(src)")
  backups = 0
  for f in files
    text = with_kernel(read(joinpath(src, f), String), spec, display)
    dst = joinpath(o.dir, f)
    # a notebook the user ran or edited is kept next to the new copy
    if isfile(dst) && read(dst, String) != text
      cp(dst, dst * ".bak"; force = true)
      backups += 1
    end
    write(dst, text)
  end
  println("Sparlectra ", sp.version, ": ", length(files), " notebooks copied to ", o.dir)
  backups > 0 && println(backups, " changed notebook(s) kept as <name>.ipynb.bak")
  println()
  println("Open the notebooks:")
  println("  Jupyter:  cd \"", o.dir, "\" and run  julia --project=. -e \"using IJulia; notebook(dir = pwd())\"")
  println("  VS Code:  open the folder, open a notebook, select the kernel \"", display, "\"")
  if !o.start
    println("Jupyter is not started (--no-start).")
    return 0
  end
  # IJulia offers to install Jupyter (through Conda) when none is found
  println("Starting Jupyter (stop it with Ctrl+C)...")
  Base.invokelatest(ijulia(:notebook); dir = o.dir)
  return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
  exit(Base.invokelatest(main, ARGS))
end
