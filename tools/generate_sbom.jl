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

# file: tools/generate_sbom.jl
# purpose: generate the SPDX Software Bill of Materials for Sparlectra with
#          PkgToSoftwareBOM in the shared build environment
#          @sparlectra-sbom-build (the package Project.toml stays free of
#          PkgToSoftwareBOM, same rule as PackageCompiler). CI runs this on
#          every published release and attaches the result as a release
#          asset; the SBOM is intentionally NOT checked into the repository.
#          Run from the checkout: `julia tools/generate_sbom.jl [outfile]`
#          (default outfile Sparlectra.spdx.json in the current directory;
#          a `--dry-run` flag prints the plan without loading
#          PkgToSoftwareBOM; SPARLECTRA_SBOM_NO_LICENSESCAN=1 disables the
#          per-package license scan).

const _REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const _BUILD_ENV = "sparlectra-sbom-build"
const _DRY_RUN = ("--dry-run" in ARGS) || get(ENV, "SPARLECTRA_SBOM_DRY_RUN", "0") == "1"
const _NO_LICENSESCAN = get(ENV, "SPARLECTRA_SBOM_NO_LICENSESCAN", "0") == "1"
const _SPARLECTRA_UUID = Base.UUID("31ce9bba-fd9d-44a1-b005-f5f509afda38")

# Pkg is deliberately NOT loaded at the top level: the dry-run path (used by
# the test-suite smoke check) must work in environments whose load path has
# no Pkg, such as the Pkg.test sandbox. The real path loads Pkg inside
# main() and resolves the module dynamically (Julia 1.12 world age).

function _log(msg::AbstractString)
  println("[generate_sbom] ", msg)
  flush(stdout)
end

# first positional (non-flag) argument names the output file
function _output_path()::String
  args = filter(a -> !startswith(a, "--"), ARGS)
  return abspath(isempty(args) ? "Sparlectra.spdx.json" : args[1])
end

function _prepare_build_env(pkgm::Module)
  _log("activating shared build environment @$(_BUILD_ENV)")
  Base.invokelatest(pkgm.activate, _BUILD_ENV; shared = true)
  if Base.find_package("PkgToSoftwareBOM") === nothing
    _log("installing PkgToSoftwareBOM into the build environment (one-time)")
    Base.invokelatest(pkgm.add, "PkgToSoftwareBOM")
  end
  # match how the script was invoked: from a checkout the repo itself is
  # developed into the build environment, otherwise the released package
  project_file = joinpath(_REPO_ROOT, "Project.toml")
  in_checkout = isfile(project_file) && occursin("name = \"Sparlectra\"", read(project_file, String))
  if in_checkout
    _log("developing Sparlectra from $(_REPO_ROOT)")
    Base.invokelatest(pkgm.develop; path = _REPO_ROOT)
  else
    _log("adding the released Sparlectra package")
    Base.invokelatest(pkgm.add, "Sparlectra")
  end
  Base.invokelatest(pkgm.instantiate)
  return nothing
end

function main()
  started = time()
  outpath = _output_path()
  if _DRY_RUN
    # parse-and-plan mode for the test suite: no environment changes, no
    # PkgToSoftwareBOM load, no network
    _log("dry run: no build environment changes, no PkgToSoftwareBOM load")
    _log("would use shared build environment @$(_BUILD_ENV)")
    _log("would describe root package Sparlectra from $(_REPO_ROOT)")
    _log("would write: $(outpath)")
    _log("license scan " * (_NO_LICENSESCAN ? "disabled (SPARLECTRA_SBOM_NO_LICENSESCAN=1)" : "enabled"))
    _log("dry run finished")
    return nothing
  end

  @eval using Pkg
  pkgm = Base.invokelatest(getfield, @__MODULE__, :Pkg)
  _prepare_build_env(pkgm)
  @eval using PkgToSoftwareBOM
  # world-age pattern: the module was loaded inside this function invocation
  sb = Base.invokelatest(getfield, @__MODULE__, :PkgToSoftwareBOM)

  originator = Base.invokelatest(getglobal(sb, :SpdxCreatorV2), "Person", "Udo Schmitz", "")
  declared_license = Base.invokelatest(getglobal(sb, :SpdxLicenseExpressionV2), "Apache-2.0")
  instructions = Base.invokelatest(
    getglobal(sb, :spdxPackageInstructions);
    spdxfile_toexclude = [basename(outpath)],
    originator = originator,
    declaredLicense = declared_license,
    copyright = "Copyright 2023-2026 Udo Schmitz",
    name = "Sparlectra",
  )
  roots = filter(p -> p.first == "Sparlectra", Base.invokelatest(getglobal(sb, :environment_rootpackages)))
  isempty(roots) && error("Sparlectra is not a root package of @$(_BUILD_ENV); run the script from a checkout or install Sparlectra first.")

  # no Manifest.toml is checked into the repository, so the described
  # dependency graph is the resolution at generation time; the comment makes
  # that explicit for SBOM consumers
  creator_comment = "Sparlectra checks in no Manifest.toml; this SBOM describes the dependency graph as freshly resolved at generation time on Julia $(VERSION)."
  creation_data = Base.invokelatest(
    getglobal(sb, :spdxCreationData);
    Name = "Sparlectra.jl Developer SBOM",
    Creators = [originator],
    CreatorComment = creator_comment,
    NamespaceURL = "https://github.com/Welthulk/Sparlectra.jl/releases",
    rootpackages = roots,
    packageInstructions = Dict(_SPARLECTRA_UUID => instructions),
    licenseScan = !_NO_LICENSESCAN,
  )
  _NO_LICENSESCAN && _log("license scan disabled (SPARLECTRA_SBOM_NO_LICENSESCAN=1)")
  _log("generating the SPDX document (scans the full dependency tree; takes a few minutes)...")
  sbom = Base.invokelatest(getglobal(sb, :generateSPDX), creation_data)
  Base.invokelatest(getglobal(sb, :writespdx), sbom, outpath)
  elapsed = round((time() - started) / 60; digits = 1)
  _log("done in $(elapsed) min: $(outpath)")
  return nothing
end

Base.invokelatest(main)
