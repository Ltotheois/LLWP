set home_dir="."

pyinstaller --onefile --windowed --icon=LLWP.ico --distpath %home_dir%\\versions %home_dir%\\LLWP\\llwp\\LLWP.py

copy %home_dir%\\LLWP\\llwp\\LLWP.py %home_dir%\\versions\\LLWP.pyw
copy %home_dir%\\LLWP\\llwp\\LLWP.py %home_dir%\\versions\\%date:~-10,2%%date:~-7,2%%date:~-4,4%.py

pip freeze >%home_dir%\\versions\\%date:~-10,2%%date:~-7,2%%date:~-4,4%.freeze

del %home_dir%\\LLWP.spec
rmdir /S /Q %home_dir%\\build 

ECHO "FINISHED"
PAUSE

