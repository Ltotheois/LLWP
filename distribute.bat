set home_dir="."

pyinstaller --onefile --windowed --icon=LLWP.ico --distpath %home_dir% %home_dir%\\LLWP.py
copy %home_dir%\\LLWP.py %home_dir%\\LLWP.pyw
REM tar.exe -a -c -f %home_dir%\\LLWP.zip LLWP.exe LLWP.svg
pip freeze >%home_dir%\\executables_graveyard\\%date:~-10,2%%date:~-7,2%%date:~-4,4%.freeze
copy %home_dir%\\LLWP.exe %home_dir%\\executables_graveyard\\%date:~-10,2%%date:~-7,2%%date:~-4,4%.exe
copy %home_dir%\\LLWP.py %home_dir%\\executables_graveyard\\%date:~-10,2%%date:~-7,2%%date:~-4,4%.py

del %home_dir%\\LLWP.spec
rmdir /S /Q %home_dir%\\build 

ECHO "FINISHED"
PAUSE

