@echo off
cd /d "%~dp0"
julia --project=. .\start_webui.jl
pause