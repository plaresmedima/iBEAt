""" 
@author: Joao Periquito 
iBEAt MODELLING Scrpit
2022
Find iBEAt motion corrected pulse sequence name (MDR output: *_mdr_moco) and execute custom modelling for: DCE, DTI, T2*
(T1 & T2 modelling are done in the main due to parallelization requirements)
"""

import datetime
import time
import pipelines.mapping as map


def main(folder,filename_log):

    list_of_series = folder.series()

    current_study = list_of_series[0].parent()
    study = list_of_series[0].new_pibling(StudyDescription=current_study.StudyDescription + '_ModellingResults')

    for i,series in enumerate(list_of_series):

        if series["SequenceName"] is not None:

            if series['SeriesDescription'] == "T2star_map_kidneys_cor-oblique_mbh_magnitude_mdr_moco":
                try:
                    start_time = time.time()
                    file = open(filename_log, 'a')
                    file.write("\n"+str(datetime.datetime.now())[0:19] + ": T2* mapping has started")
                    file.close()

                    map.T2s_Modelling(series, study=study)

                    file = open(filename_log, 'a')
                    file.write("\n"+str(datetime.datetime.now())[0:19] + ": T2* mapping was completed --- %s seconds ---" % (int(time.time() - start_time))) 
                    file.close()   
                except Exception as e: 

                    file = open(filename_log, 'a')
                    file.write("\n"+str(datetime.datetime.now())[0:19] + ": T2* mapping was NOT completed; error: "+str(e)) 
                    file.close()

            elif series['SeriesDescription'] == "DTI_kidneys_cor-oblique_fb_mdr_moco":
                try:
                    start_time = time.time()
                    file = open(filename_log, 'a')
                    file.write("\n"+str(datetime.datetime.now())[0:19] + ": DTI-FA & ADC mapping has started")
                    file.close()
                    
                    map.DTI_Modelling(series, study=study)

                    file = open(filename_log, 'a')
                    file.write("\n"+str(datetime.datetime.now())[0:19] + ": DTI-FA & ADC mapping was completed --- %s seconds ---" % (int(time.time() - start_time))) 
                    file.close()   

                except Exception as e: 
                    file = open(filename_log, 'a')
                    file.write("\n"+str(datetime.datetime.now())[0:19] + ": DTI-FA mapping was NOT completed; error: "+str(e)) 
                    file.close()

            elif series['SeriesDescription'] == "DCE_kidneys_cor-oblique_fb_mdr_moco":
                try:
                    start_time = time.time()
                    file = open(filename_log, 'a')
                    file.write("\n"+str(datetime.datetime.now())[0:19] + ": DCE-MAX mapping has started")
                    file.close()
                    
                    map.DCE_MAX_Modelling(series, study=study)

                    file = open(filename_log, 'a')
                    file.write("\n"+str(datetime.datetime.now())[0:19] + ": DCE-MAX mapping was completed --- %s seconds ---" % (int(time.time() - start_time))) 
                    file.close()   

                except Exception as e: 
                    file = open(filename_log, 'a')
                    file.write("\n"+str(datetime.datetime.now())[0:19] + ": DCE-MAX mapping was NOT completed; error: "+str(e)) 
                    file.close()

            elif series.SeriesDescription == 'MT_ON_kidneys_cor-oblique_bh_mdr_moco':
                try:
                    start_time = time.time()
                    file = open(filename_log, 'a')
                    file.write("\n"+str(datetime.datetime.now())[0:19] + ": MTR mapping has started")
                    file.close()
                    
                    map.MTR_Modelling(series, study=study)

                    file = open(filename_log, 'a')
                    file.write("\n"+str(datetime.datetime.now())[0:19] + ": MTR mapping was completed --- %s seconds ---" % (int(time.time() - start_time))) 
                    file.close()   

                except Exception as e: 
                    file = open(filename_log, 'a')
                    file.write("\n"+str(datetime.datetime.now())[0:19] + ": MTR mapping was NOT completed; error: "+str(e)) 
                    file.close()
            
            elif series.SeriesDescription == 'T1_map_kidneys_cor-oblique_mbh_magnitude_mdr_moco':
                if series.Manufacturer == 'SIEMENS':
                        
                    try:
                        start_time = time.time()
                        file = open(filename_log, 'a')
                        file.write("\n"+str(datetime.datetime.now())[0:19] + ": T1 mapping (Philips) has started")
                        file.close()
                        
                        map.T1_Modelling_Philips(series, study=study)

                        file = open(filename_log, 'a')
                        file.write("\n"+str(datetime.datetime.now())[0:19] + ": T1 mapping (Philips) was completed --- %s seconds ---" % (int(time.time() - start_time))) 
                        file.close()   

                    except Exception as e: 
                        file = open(filename_log, 'a')
                        file.write("\n"+str(datetime.datetime.now())[0:19] + ": T1 mapping (Philips) was NOT completed; error: "+str(e)) 
                        file.close()

    folder.save()
