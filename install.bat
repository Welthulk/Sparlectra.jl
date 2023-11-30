@echo off
setlocal enabledelayedexpansion

REM Get the path to the Julia program directory if not passed as a parameter
if "%1" neq "" (
    set JULIA_EXECUTABLE_PATH=%1
) else (
    for /f "tokens=*" %%i in ('where julia.exe') do set JULIA_EXECUTABLE_PATH=%%~dpi
)

REM Adjust paths if necessary
set JULIA_EXECUTABLE="%JULIA_EXECUTABLE_PATH%julia.exe"
set JULIA_PROJECT_PATH=%USERPROFILE%\.julia\dev\Sparlectra

REM Call the separate Julia script to install packages
@echo on
call "%JULIA_EXECUTABLE%" --project="%JULIA_PROJECT_PATH%" install_packages.jl
@echo off

endlocal
