NCS / NCSE ANNOTATOR
====================
Pick one or more EDF+ EEG recordings, press Start. For each one the program writes a COPY
named  annotated_<date>_<original name>.edf  (plus a short _report.txt) that you open in
EDFbrowser. Your original files are never changed.

What is marked in the copy: NCS episodes, NCSE patterns, near-miss runs (> 2.5 ED/s but < 10 s),
stretches that could not be judged, and every detected discharge ("ED").
This is a screening aid. Every detection must be reviewed by a clinical neurophysiologist.

RUNNING IT (recommended: nothing to "install", nothing for antivirus to flag)
-----------------------------------------------------------------------------
1. Install Python 3.10-3.12 from python.org (keep "tcl/tk and IDLE" ticked; tick "Add python.exe to PATH").
2. Open a Command Prompt in this folder and run once:   pip install -r requirements.txt
3. Double-click run_app.bat   (or run:  pyw ncs_app.py)

Because only the normal, signed python.exe runs, Windows Defender has nothing unusual to scan.

OPTIONAL: A STAND-ALONE PROGRAM
-------------------------------
Double-click build_exe.bat (needs step 1 above). It creates dist\NCS_Annotator\NCS_Annotator.exe.
Programs made with PyInstaller are often flagged by antivirus tools even when harmless, and nobody can
promise they will not be. To keep the risk low:
  * use the folder build above (not a single-file .exe) and do not compress it with UPX (already off);
  * build it yourself on your own PC with the latest PyInstaller;
  * sign it with a code-signing certificate (this is what removes "Windows protected your PC" /
    SmartScreen warnings; unsigned programs get them regardless of what they do);
  * if Defender still flags it, submit the file as a false positive to Microsoft's Security Intelligence
    portal; other vendors have similar forms.

PRIVACY
-------
The program makes no network connections and sends nothing anywhere. It only reads the files you pick
and writes the new files next to them (or in the folder you choose). Tick "Remove patient name and ID
from the copy" before sharing annotated files.

FILES
-----
ncs_app.py        the window (built-in tkinter only)
ncs_core.py       the detection engine (same method as the notebook)
requirements.txt  packages: numpy, scipy, mne, edfio
