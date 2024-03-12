
import wezel
from wezel.gui import Menu
from wezel.plugins import (
    pyvista,
    scipy,
    measure,
    transform,
    segment,
    align,
)
from gui import buttons, about


def launch():

    # Build iBEAt menu
    ibeat_menu = Menu('iBEAt')
    ibeat_menu.add(buttons.action_rename)
    ibeat_menu.add(buttons.action_harmonize_seqs)
    ibeat_menu.add(buttons.menu_segment)
    ibeat_menu.add(buttons.menu_mdreg)
    ibeat_menu.add(buttons.menu_map)
    ibeat_menu.add(buttons.menu_align)
    ibeat_menu.add(buttons.menu_measure)
    ibeat_menu.add_separator()
    ibeat_menu.add(buttons.action_upload)

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
    wzl.add_action(about.ibeat, menu='About')
    wzl.add_action(about.beat_dkd, menu='About')

    # Launch wezel
    wzl.show()




