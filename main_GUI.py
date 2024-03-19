# To set up a virtual environment
# -------------------------------

# For Mac OSX
# python3 -m venv .venv_ibeat           
# source .venv_ibeat/bin/activate

# Windows
# py -3 -m venv .venv           
# .venv/Scripts/activate 

# Install requirements
# pip install -r requirements.txt

# Install editable versions of requirements under development, such as:
# pip install -e C:\Users\steve\Dropbox\Software\QIB-Sheffield\dbdicom
# pip install -e C:\Users\steve\Dropbox\Software\QIB-Sheffield\wezel
# pip install -e C:\Users\steve\Dropbox\Software\QIB-Sheffield\mdreg
# pip install -e C:\Users\steve\Dropbox\Software\dcmri

# to build an executable:
# -----------------------
# distribution mode: splash screen, single file and no console
# pyinstaller --name iBEAt --clean --onefile --noconsole --additional-hooks-dir=. --splash ibeat-logo.png exec.py

# distribution mode for installer: splash screen and no console
# pyinstaller --name iBEAt --clean --noconsole --additional-hooks-dir=. --splash ibeat-logo.png exec.py

# debugging mode: multiple files & no splash & no console included.
# pyinstaller --name iBEAt --clean --noconsole --additional-hooks-dir=. exec.py

# debugging mode: multiple files & no splash & console included.
# pyinstaller --name iBEAt --clean --additional-hooks-dir=. exec.py

# Hack suggested by pyinstaller to overcome issue with infinite recursion. First manually add this line in iBEAt.spec
# import sys ; sys.setrecursionlimit(sys.getrecursionlimit() * 5)
# Then run pyinstaller as:
# pyinstaller iBEAt.spec


import gui


if __name__ == "__main__":

    gui.iBEAt.launch()
