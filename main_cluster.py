""" 
@author: Joao Periquito 
iBEAt CLUSTER MAIN Scrpit
2022
Download XNAT dataset -> Name Standardization -> Execute MDR    -> Custom Moddeling (DCE, T2*, DCE)  -> T1 & T2 modelling with parallelization (done in the main)
    XNAT_cluster.py   ->  RENAME_Cluster.py   -> MDR_Cluster.py -> MODELLING_cluster.py
(T1 & T2 modelling are done in the main due to parallelization requirements)

TO RUN THE SCRIPT YOU USE: python main_cluster.py --num n (WHERE n is an integer with the value of the XNAT dataset)
"""
import argparse
from scripts.subject_all import single_subject


if __name__ == '__main__':

    #################### INPUT ######################
    username = "*******"
    password = "*******"
    path = "//mnt//fastdata//" + username #CLUSTER PATH TO SAVE DATA, ADD YOUR LOCAL PATH IF YOU WANT TO RUN IT LOCALLY

    #################################################

    parser = argparse.ArgumentParser()
    parser.add_argument('--num',
                        dest='num',
                        help='Define the XNAT dataset',
                        type=int)

    args = parser.parse_args()

    #SELECT YOUR DATASET
    #dataset = [site, study, dataset]
    dataset = [2,1,args.num]

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

    single_subject(username, password, path, dataset)

