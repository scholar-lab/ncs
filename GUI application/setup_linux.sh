#!/usr/bin/env bash
# Creates (or updates) the "ncs-annotator" environment with Miniforge / conda / mamba.
# Usage:  bash setup_linux.sh
set -euo pipefail
cd "$(dirname "$0")"
ENV_NAME="ncs-annotator"

# Use mamba if present (faster), otherwise conda.
if command -v mamba >/dev/null 2>&1; then PM=mamba
elif command -v conda >/dev/null 2>&1; then PM=conda
else
  echo "Could not find conda or mamba on your PATH."
  echo "Open a new terminal (so Miniforge is initialised) or run:  source ~/miniforge3/etc/profile.d/conda.sh"
  exit 1
fi
echo "Using: $PM"

if conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
  echo "Environment '$ENV_NAME' already exists - updating it."
  $PM env update -n "$ENV_NAME" -f environment.yml --prune
else
  echo "Creating environment '$ENV_NAME' ..."
  $PM env create -n "$ENV_NAME" -f environment.yml
fi

echo "Checking the installation ..."
conda run -n "$ENV_NAME" python -c "
import tkinter, numpy, scipy, mne, edfio
print('OK  tkinter', tkinter.TkVersion, '| numpy', numpy.__version__, '| scipy', scipy.__version__,
      '| mne', mne.__version__, '| edfio', edfio.__version__)"

chmod +x run_app.sh
echo
echo "Setup finished. Start the program with:   ./run_app.sh"
