from wezel.gui import Menu
from . import segment
from . import monai
#from . import process

menu = Menu('iBEAt')
menu.add(segment.whole_kidney_mask)
menu.add(segment.renal_sinus_fat)
menu.add(monai.segment)
#menu.add_separator()
#menu.add(process.action_rename)
# menu.add(process.action_mdr)
# menu.add(process.action_mapping)
# menu.add(process.action_segmentation)
# menu.add(process.action_upload)


