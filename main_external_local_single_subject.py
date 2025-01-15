""" 
@author: Joao Periquito 
iBEAt Analysis MAIN LOCAL Script
2022
Load a previously downloaded XNAT dataset -> Name Standardization ->    Execute MDR (motion correction)   ->  Apply UNETR to kidney segmentation (whole kidney)  ->  Custom Moddeling (T1, T2...) ->  Align all scans with Dixon  ->   Biomarker extraction   
  select_save_folder.py                   -> pipelines.rename.py  ->             pipelines.mdr.py         ->            piplines.apply_AI_segmentation           ->      pipelines.mapping.py     ->      pipelines.align.py      ->     pipelines.export_ROI_stats -> scripts.upload.py
"""

# Windows
# py -3 -m venv .venv           
# .venv/Scripts/activate 

from scripts.single_subject_external_analysis import single_subject
import utilities.select_folder_to_save as select_save_folder

if __name__ == '__main__':

    path = select_save_folder.main()
    
    single_subject(path)