# To set up a virtual environment
# -------------------------------

# For Mac OSX
# python3 -m venv .venv_ibeat           
# source .venv_ibeat/bin/activate

# Windows
# py -3 -m venv .venv_ibeat           
# .venv_ibeat/Scripts/activate 

# Install requirements
# pip install -r requirements.txt

# Install editable versions of requirements under development
# pip install -e C:\Users\steve\Dropbox\Software\QIB-Sheffield\dbdicom
# pip install -e C:\Users\steve\Dropbox\Software\QIB-Sheffield\wezel
# pip install -e C:\Users\steve\Dropbox\Software\QIB-Sheffield\mdreg
# pip install -e C:\Users\steve\Dropbox\Software\QIB-Sheffield\iBEAt\MONAI

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


import wezel
from wezel.plugins import (
    pyvista,
    scipy,
    measure,
    transform,
    segment,
    align,
)
import gui


def build_iBEAt_menu(wzl):
    wzl.add_menu(scipy.menu_filter)
    wzl.add_menu(segment.menu)
    wzl.add_menu(align.menu)
    wzl.add_menu(transform.menu)
    wzl.add_menu(measure.menu)
    wzl.add_menu(gui.iBEAt.menu)
    wzl.add_menu(wezel.menubar.about.menu)

    wzl.add_separator(menu='View', position=5)
    wzl.add_action(pyvista.action_show_mask_surface, menu='View', position=6)
    wzl.add_action(pyvista.action_show_mask_surfaces, menu='View', position=7)
    wzl.add_action(pyvista.action_show_mask_surfaces_with_reference, menu='View', position=8)
    
    wzl.add_action(gui.about.ibeat, menu='About')
    wzl.add_action(gui.about.beat_dkd, menu='About')
    


if __name__ == "__main__":

    wzl = wezel.app(project='iBEAt')
    build_iBEAt_menu(wzl)
    wzl.show()