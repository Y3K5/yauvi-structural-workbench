@echo off
cd /d "%~dp0"
where py >nul 2>nul
if %errorlevel% equ 0 (
  py -3 start_workbench.py
) else (
  python start_workbench.py
)
pause
