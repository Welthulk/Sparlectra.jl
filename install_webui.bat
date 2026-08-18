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
rem          required), optionally builds the fast-start sysimage after
rem          asking the user, then starts the Web UI.
rem          Also works straight from the docs without a GitHub checkout;
rem          run in PowerShell in the folder where Sparlectra should live:
rem            iwr -useb https://raw.githubusercontent.com/Welthulk/Sparlectra.jl/main/install_webui.bat -OutFile install_webui.bat; .\install_webui.bat

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
rem --- 3. Dependencies ---------------------------------------------------------
echo Resolving Julia dependencies...
julia --project="%SPARLECTRA_DIR%" -e "using Pkg; Pkg.instantiate()"

rem --- 4. Optional fast-start sysimage -----------------------------------------
rem One ahead-of-time compiled image removes Julia's JIT warm-up from every
rem Web UI start (see docs: Fast Start). Optional, one-time build of roughly
rem 6 to 20 minutes. SPARLECTRA_BUILD_SYSIMAGE=1 or =0 decides without asking
rem (unattended installs); otherwise the user is asked, default is no.
set "BUILD_IMG=%SPARLECTRA_BUILD_SYSIMAGE%"
if defined BUILD_IMG goto sysimage_decided
set /p BUILD_IMG=Build the optional fast-start sysimage now? One-time build of 6-20 minutes; later Web UI starts skip the warm-up. [y/N] 
:sysimage_decided
if /i "%BUILD_IMG%"=="y" goto build_sysimage
if /i "%BUILD_IMG%"=="yes" goto build_sysimage
if "%BUILD_IMG%"=="1" goto build_sysimage
echo Skipping the sysimage build. Build it anytime with:
echo   julia "%SPARLECTRA_DIR%\tools\build_sysimage.jl"
echo or from the Web UI: Fast start page, button "Build fast-start image".
goto shortcut_prompt

:build_sysimage
echo Building the fast-start sysimage ^(this can take a while^)...
rem A failed build must never block the installation: the launcher only uses
rem an image whose metadata validates, so starting without one is always safe.
julia "%SPARLECTRA_DIR%\tools\build_sysimage.jl"
if errorlevel 1 (
  echo Sysimage build failed - continuing without it. Retry later with:
  echo   julia "%SPARLECTRA_DIR%\tools\build_sysimage.jl"
)

:shortcut_prompt
rem --- 5. Optional desktop shortcut ----------------------------------------------
rem A desktop shortcut restarts the Web UI without hunting for start_webui.bat
rem in the install folder. Default YES: the shortcut costs nothing and deleting
rem one file reverses it. SPARLECTRA_CREATE_SHORTCUT=1 or =0 decides without
rem asking (unattended installs).
set "MAKE_SHORTCUT=%SPARLECTRA_CREATE_SHORTCUT%"
if defined MAKE_SHORTCUT goto shortcut_decided
set /p MAKE_SHORTCUT=Create a desktop shortcut to start the Web UI? [Y/n]
:shortcut_decided
if /i "%MAKE_SHORTCUT%"=="n" goto skip_shortcut
if /i "%MAKE_SHORTCUT%"=="no" goto skip_shortcut
if "%MAKE_SHORTCUT%"=="0" goto skip_shortcut
rem Desktop path via .NET, never a hardcoded %%USERPROFILE%%\Desktop:
rem OneDrive-redirected desktops resolve differently. An existing shortcut of
rem the same name is overwritten so a re-install points it at the current
rem location. A failure never blocks the installation or the Web UI start.
powershell -NoProfile -Command "$ws = New-Object -ComObject WScript.Shell; $lnk = $ws.CreateShortcut((Join-Path ([Environment]::GetFolderPath('Desktop')) 'Sparlectra Web UI.lnk')); $lnk.TargetPath = '%SPARLECTRA_DIR%\start_webui.bat'; $lnk.WorkingDirectory = '%SPARLECTRA_DIR%'; $lnk.Description = 'Start the Sparlectra Web UI'; $lnk.Save()"
if errorlevel 1 (
  echo Could not create the desktop shortcut - start the Web UI manually with:
  echo   "%SPARLECTRA_DIR%\start_webui.bat"
) else (
  echo Desktop shortcut created: "Sparlectra Web UI.lnk"
)
goto start_webui
:skip_shortcut
echo Skipping the desktop shortcut.

:start_webui
rem --- 6. Start -----------------------------------------------------------------
julia --project="%SPARLECTRA_DIR%" "%SPARLECTRA_DIR%\start_webui.jl"
pause