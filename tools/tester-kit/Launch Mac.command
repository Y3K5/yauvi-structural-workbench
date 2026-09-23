#!/bin/sh
cd -- "$(dirname -- "$0")" || exit 1
if command -v python3 >/dev/null 2>&1; then
  python3 start_workbench.py
else
  echo "Python 3.10 or newer is needed for the full workbench. Open START_HERE.html for the installation-free showcase."
fi
printf '\nPress Enter to close this window. '
read -r answer
