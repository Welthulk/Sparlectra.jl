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
rem file: tools/install_webui.bat
rem purpose: install-and-start for the Sparlectra Web UI (Windows): installs
rem          Julia via winget/juliaup when missing, obtains Sparlectra at its
rem          latest tagged release (uses this checkout when the script lies
rem          inside one, otherwise downloads next to the script - no git
rem          required), offers the update when an existing copy is older
rem          than the latest release, and leaves the sysimage build to the
rem          Web UI start, then starts the Web UI.
rem          Also works straight from the docs without a GitHub checkout;
rem          run in PowerShell in the folder where Sparlectra should live:
rem            iwr -useb https://raw.githubusercontent.com/Welthulk/Sparlectra.jl/main/tools/install_webui.bat -OutFile install_webui.bat; .\install_webui.bat

setlocal
set "DIR=%~dp0"
rem the script ships in tools\ inside the checkout; anchor at the checkout
if exist "%DIR%..\Project.toml" if exist "%DIR%..\start_webui.jl" set "DIR=%DIR%..\"
rem SPARLECTRA_INSTALL_DIR names the folder that holds (or gets) Sparlectra\.
rem The update from the Web UI runs a copy of this script from its log folder
rem and points it back at the installation with this variable.
if defined SPARLECTRA_INSTALL_DIR set "DIR=%SPARLECTRA_INSTALL_DIR%\"
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
if exist "%SPARLECTRA_DIR%\Project.toml" goto existing_copy
goto fresh_download

:existing_copy
echo Using existing copy at %SPARLECTRA_DIR%
rem Update check: an existing copy stays at the release it was installed at,
rem and re-running this installer is the documented way to update. Detect a
rem newer tagged release and offer the update (default yes;
rem SPARLECTRA_UPDATE=1 or =0 decides without asking). A failed check or a
rem failed download never blocks the start with the existing copy.
set "INSTALLED="
for /f "usebackq delims=" %%V in (`powershell -NoProfile -Command "((Get-Content '%SPARLECTRA_DIR%\Project.toml') -match '^version')[0].Split('\"')[1]"`) do set "INSTALLED=%%V"
rem The latest release comes from the redirect of the latest-release page,
rem not from the GitHub REST API: that allows 60 anonymous requests per hour
rem per public IP, shared by every machine behind one company NAT. An update
rem started from the Web UI passes the tag it found as SPARLECTRA_UPDATE_TAG.
rem An update that was asked for (SPARLECTRA_UPDATE=1, or the tag the Web UI
rem passes) ends visibly: every failure prints "Update FAILED: <reason>" and
rem ends with exit code 1, so the updater reports it instead of starting the
rem old copy as if nothing happened.
set "UPDATE_REQUESTED="
if "%SPARLECTRA_UPDATE%"=="1" set "UPDATE_REQUESTED=1"
if defined SPARLECTRA_UPDATE_TAG set "UPDATE_REQUESTED=1"
set "LATEST=%SPARLECTRA_UPDATE_TAG%"
if not defined LATEST for /f "usebackq delims=" %%T in (`powershell -NoProfile -Command "$r=[System.Net.WebRequest]::Create('https://github.com/%REPO%/releases/latest'); $r.Method='HEAD'; $x=$r.GetResponse(); $u=$x.ResponseUri.AbsoluteUri; $x.Close(); if ($u -match '/releases/tag/(v[0-9]+\.[0-9]+\.[0-9]+)$') { $Matches[1] }" 2^>nul`) do set "LATEST=%%T"
if not defined LATEST (
  if defined UPDATE_REQUESTED (
    echo Update FAILED: could not determine the latest release ^(no answer from github.com^)
    exit /b 1
  )
  echo Update check skipped ^(could not reach github.com^).
  goto have_repo
)
if "v%INSTALLED%"=="%LATEST%" (
  echo Already at the latest release ^(v%INSTALLED%^).
  goto have_repo
)
set "DO_UPDATE=%SPARLECTRA_UPDATE%"
if defined DO_UPDATE goto update_decided
set /p DO_UPDATE=Installed: v%INSTALLED%, latest release: %LATEST%. Update now? [Y/n]
:update_decided
if /i "%DO_UPDATE%"=="n" goto keep_installed
if /i "%DO_UPDATE%"=="no" goto keep_installed
if "%DO_UPDATE%"=="0" goto keep_installed
echo Downloading Sparlectra %LATEST%...
rem a folder left by an earlier failed update would be picked up below
if exist "%TEMP%\sparlectra_update" rd /s /q "%TEMP%\sparlectra_update"
rem Extract into a temp dir first: the existing copy is only replaced after a
rem complete download, and it survives as Sparlectra.old (one backup slot) so
rem the update is reversible.
powershell -NoProfile -Command "Invoke-WebRequest 'https://github.com/%REPO%/archive/refs/tags/%LATEST%.zip' -OutFile '%TEMP%\sparlectra_update.zip'; Expand-Archive '%TEMP%\sparlectra_update.zip' '%TEMP%\sparlectra_update' -Force"
if not exist "%TEMP%\sparlectra_update" (
  echo Update FAILED: download of %LATEST% failed ^(https://github.com/%REPO%/archive/refs/tags/%LATEST%.zip^), the previous copy v%INSTALLED% is unchanged
  exit /b 1
)
if exist "%DIR%Sparlectra.old" rd /s /q "%DIR%Sparlectra.old"
move "%SPARLECTRA_DIR%" "%DIR%Sparlectra.old" >nul
rem Windows refuses to move a folder that a program still uses (an open
rem console whose working directory it is); the copy then stays as it is
if errorlevel 1 (
  echo Update FAILED: could not move the old copy aside, a program still uses its folder; the previous copy v%INSTALLED% is unchanged
  exit /b 1
)
powershell -NoProfile -Command "Move-Item (Get-ChildItem '%TEMP%\sparlectra_update' -Directory | Select-Object -First 1).FullName '%SPARLECTRA_DIR%'"
if not exist "%SPARLECTRA_DIR%\Project.toml" (
  if exist "%SPARLECTRA_DIR%" rd /s /q "%SPARLECTRA_DIR%"
  move "%DIR%Sparlectra.old" "%SPARLECTRA_DIR%" >nul
  echo Update FAILED: the download of %LATEST% holds no Sparlectra copy, the previous copy is restored
  exit /b 1
)
echo Updated to %LATEST%. The previous copy is kept at %DIR%Sparlectra.old ^(delete it once the new version runs^).
echo Note: an existing sysimage is stale after an update; the next Web UI start offers a rebuild.
goto have_repo
:keep_installed
echo Keeping v%INSTALLED%.
goto have_repo

:fresh_download
echo Determining the latest Sparlectra release...
for /f "usebackq delims=" %%T in (`powershell -NoProfile -Command "$r=[System.Net.WebRequest]::Create('https://github.com/%REPO%/releases/latest'); $r.Method='HEAD'; $x=$r.GetResponse(); $u=$x.ResponseUri.AbsoluteUri; $x.Close(); if ($u -match '/releases/tag/(v[0-9]+\.[0-9]+\.[0-9]+)$') { $Matches[1] }"`) do set "LATEST=%%T"
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
if errorlevel 1 if defined UPDATE_REQUESTED (
  echo Update FAILED: resolving the Julia dependencies of %SPARLECTRA_DIR% failed
  exit /b 1
)
rem the Web UI and the service layer are the SparlectraApp package under app\
julia --project="%SPARLECTRA_DIR%\app" -e "using Pkg; Pkg.instantiate()"
if errorlevel 1 if defined UPDATE_REQUESTED (
  echo Update FAILED: resolving the Julia dependencies of %SPARLECTRA_DIR%\app failed
  exit /b 1
)

rem --- 4. Sysimage -------------------------------------------------------------
rem The image is NOT built here. start_webui.jl owns that decision on every
rem platform: it checks whether a usable image exists and offers to build one,
rem so there is one implementation and one wording instead of an installer
rem copy that drifts. SPARLECTRA_BUILD_SYSIMAGE=1 builds it right away, which
rem unattended installs use to get the wait over with before the first start.
if not "%SPARLECTRA_BUILD_SYSIMAGE%"=="1" goto shortcut_prompt
echo Building the sysimage ^(this can take a while^)...
rem A failed build must never block the installation: the launcher only uses
rem an image whose metadata validates, so starting without one is safe.
julia "%SPARLECTRA_DIR%\tools\build_sysimage.jl"
if errorlevel 1 (
  echo Sysimage build failed, continuing without it. The Web UI will offer
  echo the build again at its next start.
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
rem SPARLECTRA_START_WEBUI=0 stops here; the update from the Web UI starts the
rem new copy itself, through start_webui.bat
if "%SPARLECTRA_START_WEBUI%"=="0" (
  echo Installation finished, the Web UI is not started ^(SPARLECTRA_START_WEBUI=0^).
  exit /b 0
)
julia --project="%SPARLECTRA_DIR%" "%SPARLECTRA_DIR%\start_webui.jl"
pause