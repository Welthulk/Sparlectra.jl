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
rem file: start_webui.bat
rem purpose: start the Sparlectra Web UI from this checkout (Windows).
rem          Requires an installed Julia - if none is found, this script
rem          points at install_webui.bat, which installs Julia and fetches
rem          the latest tagged Sparlectra release before starting the Web UI.

setlocal
set "DIR=%~dp0"

where julia >nul 2>nul
if errorlevel 1 (
  echo Julia is not installed ^(no 'julia' on PATH^).
  echo Run   %DIR%install_webui.bat   instead:
  echo it installs Julia ^(via juliaup^), fetches the latest tagged Sparlectra release, and starts the Web UI.
  pause
  exit /b 1
)

rem A fresh checkout has no Manifest.toml yet - resolve dependencies once.
if not exist "%DIR%Manifest.toml" (
  echo First start: resolving Julia dependencies ^(one-time^)...
  julia --project="%DIR%." -e "using Pkg; Pkg.instantiate()"
)

julia --project="%DIR%." "%DIR%start_webui.jl"
pause
