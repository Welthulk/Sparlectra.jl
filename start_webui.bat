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

rem Fast start: use the optional PackageCompiler sysimage when it matches the
rem running Julia and the checkout Manifest (contract in sysimage_meta.toml,
rem written by tools\build_sysimage.jl). The path mirrors
rem Sparlectra.default_webui_output_root in src\webui\webui.jl; keep both in
rem sync. SPARLECTRA_NO_SYSIMAGE=1 skips the image unconditionally. A stale
rem or missing image never blocks the start.
set "USE_SYSIMAGE="
if "%SPARLECTRA_NO_SYSIMAGE%"=="1" goto run_webui
set "SYSIMAGE=%LOCALAPPDATA%\Sparlectra\WebUI\sysimage\sparlectra.dll"
set "SYSMETA=%LOCALAPPDATA%\Sparlectra\WebUI\sysimage\sysimage_meta.toml"
if not exist "%SYSIMAGE%" goto run_webui
if not exist "%SYSMETA%" goto run_webui
set "META_JULIA="
for /f "tokens=3" %%a in ('findstr /b "julia_version" "%SYSMETA%"') do set "META_JULIA=%%~a"
set "RUN_JULIA="
for /f "tokens=3" %%a in ('julia --version') do set "RUN_JULIA=%%a"
set "META_SHA="
for /f "tokens=3" %%a in ('findstr /b "manifest_sha256" "%SYSMETA%"') do set "META_SHA=%%~a"
set "CUR_SHA="
for /f "skip=1 tokens=1" %%a in ('certutil -hashfile "%DIR%Manifest.toml" SHA256') do (
  if not defined CUR_SHA set "CUR_SHA=%%a"
)
if not defined CUR_SHA goto sysimage_stale
if not "%META_JULIA%"=="%RUN_JULIA%" goto sysimage_stale
if /i not "%META_SHA%"=="%CUR_SHA%" goto sysimage_stale
set "USE_SYSIMAGE=1"
goto run_webui

:sysimage_stale
echo Fast start: sysimage stale, starting without it; rebuild via tools\build_sysimage.jl

:run_webui
REM Multi-core by default: the threaded surfaces need Julia THREADS, fixed at
REM process start. "auto" uses all cores; an explicit user setting wins.
if not defined JULIA_NUM_THREADS set "JULIA_NUM_THREADS=auto"
if defined USE_SYSIMAGE (
  echo Fast start: using sysimage %SYSIMAGE%
  julia "-J%SYSIMAGE%" --project="%DIR%." "%DIR%start_webui.jl"
) else (
  julia --project="%DIR%." "%DIR%start_webui.jl"
)
pause
