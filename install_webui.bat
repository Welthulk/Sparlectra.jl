@echo off
rem Copyright 2023-2026 Udo Schmitz
rem
rem Licensed under the Apache License, Version 2.0 (the "License");
rem you may not use this file except in compliance with the License.
rem You may obtain a copy of the License at
rem
rem     http://www.apache.org/licenses/LICENSE-2.0
rem
rem Unless required by applicable law or agreed to in writing, software
rem distributed under the License is distributed on an "AS IS" BASIS,
rem WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
rem See the License for the specific language governing permissions and
rem limitations under the License.
rem
rem file: install_webui.bat
rem purpose: install-and-start for the Sparlectra Web UI (Windows): installs
rem          Julia via winget/juliaup when missing, obtains Sparlectra at its
rem          latest tagged release (uses this checkout when the script lies
rem          inside one, otherwise downloads next to the script - no git
rem          required), then starts the Web UI.

setlocal
set "DIR=%~dp0"
set "REPO=Welthulk/Sparlectra.jl"

rem --- 1. Julia --------------------------------------------------------------
where julia >nul 2>nul
if not errorlevel 1 goto have_julia
echo Julia not found - installing via winget ^(juliaup^)...
winget install --id 9NJNWW8PVKMN --source msstore --accept-package-agreements --accept-source-agreements
if errorlevel 1 (
  echo winget installation failed. Install Julia manually: https://julialang.org/downloads/
  pause
  exit /b 1
)
rem juliaup puts julia into %USERPROFILE%\.juliaup\bin - make it visible NOW.
set "PATH=%USERPROFILE%\.juliaup\bin;%PATH%"
where julia >nul 2>nul
if errorlevel 1 (
  echo Julia was installed but is not on PATH yet. Open a NEW terminal and run start_webui.bat.
  pause
  exit /b 1
)
:have_julia
julia --version

rem --- 2. Sparlectra at the latest tagged release ----------------------------
rem The pair Project.toml + start_webui.jl marks a Sparlectra checkout; the
rem script then uses it in place instead of downloading a second copy.
set "SPARLECTRA_DIR=%DIR%."
if exist "%DIR%Project.toml" if exist "%DIR%start_webui.jl" goto have_repo

set "SPARLECTRA_DIR=%DIR%Sparlectra"
if exist "%SPARLECTRA_DIR%\Project.toml" (
  echo Using existing copy at %SPARLECTRA_DIR%
  goto have_repo
)
echo Determining the latest Sparlectra release...
for /f "usebackq delims=" %%T in (`powershell -NoProfile -Command "(Invoke-RestMethod 'https://api.github.com/repos/%REPO%/tags')[0].name"`) do set "LATEST=%%T"
if not defined LATEST (
  echo Could not determine the latest release tag.
  pause
  exit /b 1
)
echo Downloading Sparlectra %LATEST% ^(no git required^)...
powershell -NoProfile -Command "Invoke-WebRequest 'https://github.com/%REPO%/archive/refs/tags/%LATEST%.zip' -OutFile '%TEMP%\sparlectra_release.zip'; Expand-Archive '%TEMP%\sparlectra_release.zip' '%TEMP%\sparlectra_release' -Force; Move-Item (Get-ChildItem '%TEMP%\sparlectra_release' -Directory | Select-Object -First 1).FullName '%SPARLECTRA_DIR%'"
if not exist "%SPARLECTRA_DIR%\Project.toml" (
  echo Download failed - %SPARLECTRA_DIR% is incomplete.
  pause
  exit /b 1
)

:have_repo
rem --- 3. Dependencies + start -----------------------------------------------
echo Resolving Julia dependencies...
julia --project="%SPARLECTRA_DIR%" -e "using Pkg; Pkg.instantiate()"
julia --project="%SPARLECTRA_DIR%" "%SPARLECTRA_DIR%\start_webui.jl"
pause
