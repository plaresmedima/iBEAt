""" 
@author: Joao Periquito 
iBEAt CLUSTER MAIN Scrpit
2022
Download XNAT dataset -> Name Standardization -> Execute MDR    -> Custom Moddeling (DCE, T2*, DCE)  -> T1 & T2 modelling with parallelization (done in the main)
    XNAT_cluster.py   ->  RENAME_Cluster.py   -> MDR_Cluster.py -> MODELLING_cluster.py
(T1 & T2 modelling are done in the main due to parallelization requirements)

TO RUN THE SCRIPT YOU USE: python main_cluster.py --num n (WHERE n is an integer with the value of the XNAT dataset)
"""

# To develop the application
# --------------------------
# py -3 -m venv .venv           # create virtual environment
# .venv/Scripts/activate        # activate virtual environment
import os
import datetime
import dbdicom as db

import pipelines.rename as rename
import pipelines.mdr as mdr
import pipelines.mapping as map
import scripts.upload as upload
import scripts.segmentation as segmentation

if __name__ == '__main__':

    #################### INPUT ######################
    username = "**********"
    password = "**********"
    #path = "//mnt//fastdata//" + username #CLUSTER PATH TO SAVE DATA, ADD YOUR LOCAL PATH IF YOU WANT TO RUN IT LOCALLY
    path = "C://Users//md1jdsp//Desktop"
    #################################################

    # parser = argparse.ArgumentParser()
    # parser.add_argument('--num',
    #                     dest='num',
    #                     help='Define the XNAT dataset',
    #                     type=int)

    # args = parser.parse_args()

    #dataset = [2,1,args.num]

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

    #ExperimentName = xnat.main(username, password, path, dataset)
    ExperimentName = "iBE-2128_008_baseline"
    pathScan = path + "//" + ExperimentName
    filename_log = pathScan + datetime.datetime.now().strftime('%Y%m%d_%H%M_') + "MDRauto_LogFile.txt" #TODO FIND ANOTHER WAY TO GET A PATH
    
    try: 
        UsedCores = int(len(os.sched_getaffinity(0)))
    except: 
        UsedCores = int(os.cpu_count())

    folder = db.database(path=pathScan)
    folder.set_log(filename_log)
    folder.log("Analysis of " + pathScan.split('//')[-1] + " has started!")
    folder.log("CPU cores: " + str(UsedCores))
    
    try:
        rename.main(folder)
    except Exception as e:
        folder.log("Renaming was NOT completed; error: " + str(e))

    try:
        mdr.main(folder)
    except Exception as e:
        folder.log("Renaming was NOT completed; error: " + str(e))

    try:
        map.main(folder)
    except Exception as e:
        folder.log("Modelling was NOT completed; error: " + str(e))

    try:
        segmentation.main(folder,pathScan)
    except Exception as e:
        folder.log("Segmentation was NOT completed; error: " + str(e))

    upload.main(pathScan, filename_log)
