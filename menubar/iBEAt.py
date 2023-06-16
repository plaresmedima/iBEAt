from wezel.gui import Menu
from . import segment
from . import pipelines
#from . import monai
from . import process
from . import mdr

menu = Menu('iBEAt')
menu.add(segment.whole_kidney_mask)
menu.add_separator()
menu.add(segment.renal_sinus_fat)
menu.add(pipelines.asl_perfusion)
#menu.add(monai.segment)
menu.add_separator()
menu.add(mdr.action_T2star)
menu.add(mdr.action_T1)
menu.add(mdr.action_T2)
menu.add(mdr.action_DTI)
menu.add(mdr.action_MT)
menu.add(mdr.action_DCE)
menu.add_separator()
menu.add(process.action_rename)
menu.add(process.action_mdr)
menu.add(process.action_mapping)
menu.add(process.action_segmentation)
menu.add(process.action_upload)



