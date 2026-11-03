@echo off
rem Optional: builds a folder-style Windows program (not a single .exe) with PyInstaller.
rem Run this on your own Windows PC with Python 3.10-3.12 installed.
py -3 -m venv build_env || goto :err
call build_env\Scripts\activate.bat
python -m pip install --upgrade pip
pip install -r requirements.txt pyinstaller || goto :err
pyinstaller --noconfirm --clean --onedir --windowed --noupx --name NCS_Annotator --hidden-import ncs_core --collect-data mne ncs_app.py || goto :err
echo.
echo Done. The program is in dist\NCS_Annotator\  (start NCS_Annotator.exe there).
echo Keep the whole folder together. See README.txt about antivirus warnings.
pause
exit /b 0
:err
echo Something went wrong - see the messages above.
pause
exit /b 1
