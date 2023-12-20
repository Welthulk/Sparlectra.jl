@echo off
setlocal enabledelayedexpansion

REM Determine the path to the Julia program directory
for /f "tokens=*" %%i in ('where julia.exe') do set JULIA_EXECUTABLE_PATH=%%~dpi

REM Check whether at least one Julia executable was found
if "%JULIA_EXECUTABLE_PATH%"=="" (
    echo Julia nicht gefunden. Stellen Sie sicher, dass Julia installiert ist und im System-PATH enthalten ist.
    exit /b 1
)

REM Select the first Julia executable found
set JULIA_EXECUTABLE="%JULIA_EXECUTABLE_PATH%julia.exe"

REM Adjust paths if necessary
set JULIA_PROJECT_PATH=%USERPROFILE%\.julia\dev\Sparlectra

REM Run the Julia script with the specified parameters
@echo on
"%JULIA_EXECUTABLE%" --color=yes --project="%JULIA_PROJECT_PATH%" "%JULIA_PROJECT_PATH%\test\testparser.jl"
@echo off

pause

endlocal
