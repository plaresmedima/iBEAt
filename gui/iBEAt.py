from wezel.gui import Menu
from . import contrasts
from . import steps
from . import mdr
from . import mapping

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


def launch():

    # Build iBEAt menu
    ibeat_menu = Menu('iBEAt')
    ibeat_menu.add(steps.action_rename)
    ibeat_menu.add(steps.menu_segment)
    ibeat_menu.add_separator()
    ibeat_menu.add(contrasts.asl_perfusion)
    ibeat_menu.add(contrasts.menu_T1)
    ibeat_menu.add_separator()
    ibeat_menu.add(mdr.action_T2star)
    ibeat_menu.add(mdr.action_T1)
    ibeat_menu.add(mdr.action_T2)
    ibeat_menu.add(mdr.action_DTI)
    ibeat_menu.add(mdr.action_MT)
    ibeat_menu.add(mdr.action_DCE)
    ibeat_menu.add_separator()
    ibeat_menu.add(mapping.action_T2star)
    ibeat_menu.add(mapping.action_T1)
    ibeat_menu.add(mapping.action_DTI)
    ibeat_menu.add(mapping.action_MT)
    ibeat_menu.add(mapping.action_DCE)
    ibeat_menu.add_separator()
    ibeat_menu.add(steps.action_mdr)
    ibeat_menu.add(steps.action_mapping)
    ibeat_menu.add(steps.action_upload)

    # Build iBEAt wezel project
    wzl = wezel.app(project='iBEAt')

    # Add non-default menus
    wzl.add_menu(scipy.menu_filter)
    wzl.add_menu(segment.menu)
    wzl.add_menu(align.menu)
    wzl.add_menu(transform.menu)
    wzl.add_menu(measure.menu)
    wzl.add_menu(ibeat_menu)
    wzl.add_menu(wezel.menubar.about.menu)

    # Extend default view menu
    wzl.add_separator(menu='View', position=5)
    wzl.add_action(pyvista.action_show_mask_surface, menu='View', position=6)
    wzl.add_action(pyvista.action_show_mask_surfaces, menu='View', position=7)
    wzl.add_action(pyvista.action_show_mask_surfaces_with_reference, menu='View', position=8)

    # Extend about menu
    wzl.add_action(gui.about.ibeat, menu='About')
    wzl.add_action(gui.about.beat_dkd, menu='About')

    # Launch wezel
    wzl.show()




