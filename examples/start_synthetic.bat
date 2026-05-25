@REM filepath: c:\Users\scud\.julia\dev\Sparlectra\run_perf_benchmark.bat
@echo off
setlocal enabledelayedexpansion

set SCRIPT_DIR=%~dp0
set CONFIG_FILE=%SCRIPT_DIR%src\examples\exp_synthetic_tiled_grid_pf_perf.yaml

if exist "%CONFIG_FILE%" (
    echo Starte Benchmark mit Konfiguration: %CONFIG_FILE%
    julia --project="%SCRIPT_DIR%" "%SCRIPT_DIR%src\examples\exp_synthetic_tiled_grid_pf_perf.jl" "%CONFIG_FILE%" %*
) else (
    echo Starte Benchmark mit Standard-Konfiguration
    julia --project="%SCRIPT_DIR%" "%SCRIPT_DIR%src\examples\exp_synthetic_tiled_grid_pf_perf.jl" %*
)