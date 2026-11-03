#!/usr/bin/env bash
# Starts the NCS / NCSE Annotator window inside the conda environment.
cd "$(dirname "$0")"
exec conda run --no-capture-output -n ncs-annotator python ncs_app.py "$@"
