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
rem          points at tools\install_webui.bat, which installs Julia and fetches
rem          the latest tagged Sparlectra release before starting the Web UI.

setlocal
set "DIR=%~dp0"

where julia >nul 2>nul
if errorlevel 1 (
  echo Julia is not installed ^(no 'julia' on PATH^).
  echo Run   %DIR%tools\install_webui.bat   instead:
  echo it installs Julia ^(via juliaup^), fetches the latest tagged Sparlectra release, and starts the Web UI.
  pause
  exit /b 1
)

rem A fresh checkout has no Manifest.toml yet - resolve dependencies once.
if not exist "%DIR%Manifest.toml" (
  echo First start: resolving Julia dependencies ^(one-time^)...
  julia --project="%DIR%." -e "using Pkg; Pkg.instantiate()"
)

REM Multi-core by default: the threaded surfaces need Julia THREADS, fixed at
REM process start. "auto" uses all cores; an explicit user setting wins.
REM Note: "auto" follows the process affinity mask, not the raw core count.
REM On machines with E-cores or a preset affinity you may get fewer threads
REM than expected - set an explicit number to compare.
if not defined JULIA_NUM_THREADS set "JULIA_NUM_THREADS=auto"

REM Parallel garbage collection. Default is a quarter of the compute threads,
REM which is often the real bottleneck on allocation-heavy runs.
if not defined JULIA_NUM_GC_THREADS set "JULIA_NUM_GC_THREADS=4"

REM No startup file. The Web UI is a server process, not a REPL, and a personal
REM startup.jl usually loads Revise: measured on the Sparlectra sysimage,
REM loading Revise INVALIDATES 1530 precompiled method instances of that image,
REM which are then inferred again on first use. That is the "still slow with the
REM image" report from 2026-09-07, and it was paid twice, because the launcher
REM and the relaunched child both read the file. Set SPARLECTRA_STARTUP_FILE=yes
REM to get the old behavior back.
if not defined SPARLECTRA_STARTUP_FILE set "SPARLECTRA_STARTUP_FILE=no"

REM Process priority. HIGH is the highest class that is safe for a normal
REM desktop; REALTIME needs admin rights and can starve the OS. The child
REM process that start_webui.jl relaunches through the sysimage inherits this.
REM Set SPARLECTRA_PRIORITY=NORMAL (or ABOVENORMAL, BELOWNORMAL) to override.
if not defined SPARLECTRA_PRIORITY set "SPARLECTRA_PRIORITY=HIGH"

rem start_webui.jl checks for a usable sysimage, offers to build one when it
rem is missing or outdated, and relaunches itself through the image. The same
rem code runs on Linux, macOS and Windows. Arguments are passed on:
rem --rebuild-sysimage forces a fresh build, --no-sysimage skips it.
rem The empty "" is the window title start expects - without it start would
rem swallow the next quoted argument. /B keeps it in this console, /WAIT keeps
rem the pause below meaningful.
start "" /%SPARLECTRA_PRIORITY% /B /WAIT julia --startup-file=%SPARLECTRA_STARTUP_FILE% --project="%DIR%." "%DIR%start_webui.jl" %*
pause