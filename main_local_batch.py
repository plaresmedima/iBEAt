""" 
@author: Joao Periquito 
iBEAt Analysis MAIN LOCAL Scrpit
2022
Download XNAT dataset -> Name Standardization ->    Execute MDR   ->  Apply UNETR to kidney segmentation  ->  Custom Moddeling (T1, T2...) ->     Biomarker extraction   -> Google Drive Upload
  pipelines.xnat.py   -> pipelines.rename.py  -> pipelines.mdr.py ->    piplines.apply_AI_segmentation    ->    pipelines.mapping.py       -> pipelines.export_ROI_stats -> scripts.upload.py
"""

from scripts.subject_all import single_subject
import scripts.XNAT_credentials as XNAT_cred
import scripts.select_folder_to_save as select_save_folder

if __name__ == '__main__':

    number_of_cases_to_run = 20

    path = select_save_folder.main()
    
    #XNAT Credentials
    username, password = XNAT_cred.main()

    #SELECT YOUR DATASET
    #dataset = [site, study, dataset] see below "EXAMPLE DATASET SELECTION"
    
    for i in range (number_of_cases_to_run):
      dataset = [2,1,i]

    ################################################# EXAMPLE DATASET SELECTION #############################################################
    #DATASET CODE FOR LEEDS
    #  (FIRST NUMBER)                 (SECOND NUMBER)                                 (THIRD NUMBER - INNPUT from --num when you run the main script: python main_cluster.py --num n)
    #  2: BEAt-DKD-WP4-Bordeaux    (selected) BEAt-DKD-WP4-Leeds                 (selected) BEAt-DKD-WP4-Leeds -> (selected) Leeds_Patients
    #  3: BEAt-DKD-WP4-Exeter       ->0: Leeds_Patients                            0: Leeds_Patient_4128001
    #  4: BEAt-DKD-WP4-Turku          1: Leeds_volunteer_repeatability_study       1: Leeds_Patient_4128002
    #  5: BEAt-DKD-WP4-Bari           2: Leeds_Phantom_scans                       2: Leeds_Patient_4128004
    #->6: BEAt-DKD-WP4-Leeds          3: Leeds_RAVE_Reconstructions                       ..........
    #  7: BEAt-DKD-WP4-Sheffield      4: Leeds_setup_scans                      ->14: Leeds_Patient_4128015
    #########################################################################################################################################

    #function responsable for ibeat analysis of a single subject (processed images, logs, and csv results are exported to Google Drive)
      single_subject(username, password, path, dataset)