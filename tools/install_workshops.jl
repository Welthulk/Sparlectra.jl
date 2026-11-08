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

using Pkg, SHA, TOML, Dates

const REPO_URL = "https://github.com/Welthulk/Sparlectra.jl"
# what the installer last wrote into the folder (notebook name => SHA-256),
# so a rerun can tell a notebook the user changed from an old copy
const RECORD_FILE = ".sparlectra_workshops.toml"

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

# One kernel per workshop folder: the spec name carries the folder name and a
# short hash of its full path, so a second --dir gets its own kernel instead
# of redirecting the first folder's notebooks into the second environment.
function kernel_names(dir::String)
  h = bytes2hex(sha256(dir))[1:8]
  base = replace(lowercase(basename(dir)), r"[^a-z0-9._-]" => "-")
  isempty(base) && (base = "folder")
  return (spec = "sparlectra-$(base)-$(h)", display = "Sparlectra workshops ($(basename(dir)) $(h))")
end

# Sparlectra is removed first when the folder already has it, so a second run
# takes the source asked for now and, for the registered release, the newest
# version. Project.toml and Manifest.toml are saved before and put back when
# the removal or the new install fails, so a failed run leaves the folder as
# it was instead of without Sparlectra.
function install_sparlectra(o)
  saved = Dict{String,String}()
  for f in ("Project.toml", "Manifest.toml")
    p = joinpath(o.dir, f)
    isfile(p) && (saved[p] = read(p, String))
  end
  try
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
  catch
    # not a handler: any failure (resolver, network, bad --rev) is rethrown
    # unchanged, only the environment files are restored first
    for (p, text) in saved
      write(p, text)
    end
    isempty(saved) || println("Install failed; the previous Project.toml and Manifest.toml are restored.")
    rethrow()
  end
  return nothing
end

# Copies the workshop notebooks `files` from `src` into `dir`, with the
# kernel of the folder. A notebook the user ran or edited is kept next to the
# new copy as <name>.ipynb.<date>-<time>.bak; a copy that is still what the
# installer wrote last time (RECORD_FILE) is replaced without a backup.
# Without a record (a folder set up by an older installer) every differing
# copy counts as changed. Returns the number of replaced notebooks and the
# names of the backups.
function copy_notebooks(src::String, dir::String, files::Vector{String}, spec::String, display::String; stamp::String = Dates.format(now(), "yyyymmdd-HHMMSS"))
  record_path = joinpath(dir, RECORD_FILE)
  record = isfile(record_path) ? TOML.parsefile(record_path) : Dict{String,Any}()
  written = Dict{String,Any}()
  backups = String[]
  updated = 0
  for f in files
    text = with_kernel(read(joinpath(src, f), String), spec, display)
    dst = joinpath(dir, f)
    if isfile(dst)
      old = read(dst, String)
      if old != text
        updated += 1
        if get(record, f, "") != bytes2hex(sha256(old))
          bak = string(dst, ".", stamp, ".bak")
          cp(dst, bak)
          push!(backups, basename(bak))
        end
      end
    end
    write(dst, text)
    written[f] = bytes2hex(sha256(text))
  end
  open(io -> TOML.print(io, written), record_path, "w")
  return (updated = updated, backups = backups)
end

function main(args::Vector{String})
  o = parse_options(args)
  mkpath(o.dir)
  println("Workshop folder: ", o.dir)
  Pkg.activate(o.dir)

  install_sparlectra(o)
  Pkg.add("IJulia")
  Pkg.precompile()

  # Jupyter kernel that runs in this folder's environment, whatever folder
  # a notebook is opened from. IJulia is loaded at run time, after the
  # activation above, so every access goes through invokelatest.
  ijulia(name) = Base.invokelatest(getfield, Base.require(Main, :IJulia), name)
  kn = kernel_names(o.dir)
  kdir = Base.invokelatest(ijulia(:installkernel), kn.display, "--project=$(o.dir)"; specname = kn.spec, displayname = kn.display)
  spec = basename(kdir)
  display = kn.display
  println("Jupyter kernel: \"", display, "\" (", spec, ")")

  # the notebooks of exactly the Sparlectra just installed
  sp = only(d for d in values(Pkg.dependencies()) if d.name == "Sparlectra")
  src = joinpath(sp.source, "notebooks")
  files = sort(filter(f -> startswith(f, "workshop_") && endswith(f, ".ipynb"), readdir(src)))
  isempty(files) && error("no workshop notebooks in $(src)")
  c = copy_notebooks(src, o.dir, files, spec, display)
  println("Sparlectra ", sp.version, ": ", length(files), " notebooks in ", o.dir, " (", c.updated, " updated)")
  isempty(c.backups) || println("Notebooks you changed are kept as: ", join(c.backups, ", "))
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
