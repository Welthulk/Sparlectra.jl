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
# file: tools/measure_precompile.ps1
# purpose: measure the precompile time of the library and, where present,
#          the application package for one or more git refs on this
#          machine. Each ref is cloned into its own directory under a work
#          folder, the environments are instantiated first (downloads and
#          the automatic precompile of instantiate stay out of the timed
#          part), then Pkg.precompile() is timed for . and for app/ with a
#          fresh compile cache each. Output: one table row per ref.
#
# usage:   powershell -ExecutionPolicy Bypass -File tools/measure_precompile.ps1 v0.16.2 dev/r0.17.0
#          optional: -WorkDir <folder>  (default: a folder under %TEMP%)
#          optional: -Repo <url>        (default: the public repository)
#          optional: -Depot <folder>    a Julia depot of its own for the runs
#
# depot:   without -Depot the script empties the compile cache of the two
#          packages in the shared depot (~/.julia/compiled), so the next
#          start of any other Sparlectra checkout on this machine
#          precompiles again (minutes on Windows). With -Depot everything
#          (registry, packages, cache) lands in that folder instead, which
#          costs the downloads once and leaves the shared depot untouched.
#
# note:    Pkg.instantiate() precompiles on its own and prints
#          "N dependencies successfully precompiled in X seconds"; that X is
#          the same work as the timed Pkg.precompile() here, on a cache that
#          was empty. The script empties the compile cache of the two
#          packages before the timed call, so both numbers describe a cold
#          precompile. A ref without app/ (v0.16.2 and earlier) yields the
#          library number only.

param(
    [Parameter(Mandatory = $true, Position = 0, ValueFromRemainingArguments = $true)]
    [string[]] $Refs,
    [string] $WorkDir = (Join-Path $env:TEMP "sparlectra_precompile"),
    [string] $Repo = "https://github.com/Welthulk/Sparlectra.jl.git",
    [string] $Depot = ""
)

$ErrorActionPreference = "Stop"
New-Item -ItemType Directory -Force -Path $WorkDir | Out-Null
$depotRoot = Join-Path $env:USERPROFILE ".julia"
if ($Depot -ne "") {
    New-Item -ItemType Directory -Force -Path $Depot | Out-Null
    $env:JULIA_DEPOT_PATH = $Depot
    $depotRoot = $Depot
    Write-Host "Using the depot $Depot for every Julia call of this script."
} else {
    Write-Host "Note: the compile cache of Sparlectra and SparlectraApp in $depotRoot is emptied before each timed run; other checkouts on this machine precompile again afterwards."
}
$compiled = Join-Path $depotRoot "compiled"

function Clear-CompileCache([string] $package) {
    # every Julia minor version has its own cache folder; drop the package's
    # entries in all of them so the timed precompile starts cold
    Get-ChildItem -Path $compiled -Directory -ErrorAction SilentlyContinue | ForEach-Object {
        $dir = Join-Path $_.FullName $package
        if (Test-Path $dir) { Remove-Item -Recurse -Force $dir }
    }
}

function Measure-Precompile([string] $project) {
    $t = Measure-Command {
        & julia --startup-file=no --project=$project -e "using Pkg; Pkg.precompile()" | Out-Null
    }
    return [math]::Round($t.TotalSeconds, 1)
}

$rows = @()
foreach ($ref in $Refs) {
    $name = ($ref -replace '[\\/:]', '_')
    $dir = Join-Path $WorkDir $name
    if (Test-Path $dir) { Remove-Item -Recurse -Force $dir }
    Write-Host "== $ref -> $dir"
    & git clone --quiet --depth 1 --branch $ref $Repo $dir
    if ($LASTEXITCODE -ne 0) { throw "git clone of $ref failed" }

    # downloads and the automatic precompile of instantiate are not measured
    & julia --startup-file=no --project=$dir -e "using Pkg; Pkg.instantiate()"
    $hasApp = Test-Path (Join-Path $dir "app\Project.toml")
    if ($hasApp) {
        & julia --startup-file=no --project=(Join-Path $dir "app") -e "using Pkg; Pkg.instantiate()"
    }

    Clear-CompileCache "Sparlectra"
    $library = Measure-Precompile $dir

    $app = "-"
    if ($hasApp) {
        Clear-CompileCache "SparlectraApp"
        $app = Measure-Precompile (Join-Path $dir "app")
    }
    $rows += [pscustomobject]@{ Version = $ref; "Library [s]" = $library; "App [s]" = $app }
}

Write-Host ""
Write-Host ("Julia " + (& julia --startup-file=no -e "print(VERSION)") + " on " + $env:COMPUTERNAME)
$rows | Format-Table -AutoSize
