""" 
@author: Joao Periquito 
iBEAt MDR Scrpit
2022
Find iBEAt standard pulse sequence name  and execute MDR for: DCE, DTI, T1, T2, T2*, MT
"""

import datetime
import time

import pipelines.mdr as mdr



def main(folder,filename_log):

    current_study = folder.series()[0].parent()
    study = folder.series()[0].new_pibling(StudyDescription=current_study.StudyDescription + '_MDRresults')

    for series in folder.series():

        SeqName = series["SeriesDescription"]

        if SeqName is not None:
            print(SeqName)

            if SeqName == "T2star_map_kidneys_cor-oblique_mbh_magnitude":
                try:
                    start_time = time.time()
                    file = open(filename_log, 'a')
                    file.write("\n"+str(datetime.datetime.now())[0:19] + ": T2* motion correction has started")
                    #file.write("\n"+str(datetime.datetime.now())[0:19] + ": RAM Used (GB): " + str(psutil.virtual_memory()[3]/1000000000))
                    file.close()

                    mdr.MDRegT2star(series, study=study)

                    file = open(filename_log, 'a')
                    file.write("\n"+str(datetime.datetime.now())[0:19] + ": T2* motion correction was completed --- %s seconds ---" % (int(time.time() - start_time))) 
                    #file.write("\n"+str(datetime.datetime.now())[0:19] + ": RAM Used (GB): " + str(psutil.virtual_memory()[3]/1000000000))
                    file.close()   

                except Exception as e: 
                    file = open(filename_log, 'a')
                    file.write("\n"+str(datetime.datetime.now())[0:19] + ": T2* motion correction was NOT completed; error: "+str(e)) 
                    file.close()

            elif SeqName == "T1map_kidneys_cor-oblique_mbh_magnitude":
                try:
                    start_time = time.time()
                    file = open(filename_log, 'a')
                    file.write("\n"+str(datetime.datetime.now())[0:19] + ": T1 motion correction has started")
                    #file.write("\n"+str(datetime.datetime.now())[0:19] + ": RAM Used (GB): " + str(psutil.virtual_memory()[3]/1000000000))
                    file.close()

                    mdr.MDRegT1(series, study=study)

                    file = open(filename_log, 'a')
                    file.write("\n"+str(datetime.datetime.now())[0:19] + ": T1 motion correction was completed --- %s seconds ---" % (int(time.time() - start_time))) 
                    #file.write("\n"+str(datetime.datetime.now())[0:19] + ": RAM Used (GB): " + str(psutil.virtual_memory()[3]/1000000000))
                    file.close()   

                except Exception as e: 
                    file = open(filename_log, 'a')
                    file.write("\n"+str(datetime.datetime.now())[0:19] + ": T1 motion correction was NOT completed; error: "+str(e)) 
                    file.close()

            elif SeqName == "T2map_kidneys_cor-oblique_mbh_magnitude":
                try:
                    start_time = time.time()
                    file = open(filename_log, 'a')
                    file.write("\n"+str(datetime.datetime.now())[0:19] + ": T2 motion correction has started")
                    #file.write("\n"+str(datetime.datetime.now())[0:19] + ": RAM Used (GB): " + str(psutil.virtual_memory()[3]/1000000000))
                    file.close()

                    mdr.MDRegT2(series, study=study)

                    file = open(filename_log, 'a')
                    file.write("\n"+str(datetime.datetime.now())[0:19] + ": T2 motion correction was completed --- %s seconds ---" % (int(time.time() - start_time))) 
                    #file.write("\n"+str(datetime.datetime.now())[0:19] + ": RAM Used (GB): " + str(psutil.virtual_memory()[3]/1000000000))
                    file.close()   

                except Exception as e: 
                    file = open(filename_log, 'a')
                    file.write("\n"+str(datetime.datetime.now())[0:19] + ": T2 motion correction was NOT completed; error: "+str(e)) 
                    file.close()   

            elif SeqName == "DTI_kidneys_cor-oblique_fb":
                try:
                    start_time = time.time()
                    file = open(filename_log, 'a')
                    file.write("\n"+str(datetime.datetime.now())[0:19] + ": DTI motion correction has started")
                    file.close()

                    mdr.MDRegDTI(series, study=study)

                    file = open(filename_log, 'a')
                    file.write("\n"+str(datetime.datetime.now())[0:19] + ": DTI motion correction was completed --- %s seconds ---" % (int(time.time() - start_time))) 
                    file.close()   

                except Exception as e: 
                    file = open(filename_log, 'a')
                    file.write("\n"+str(datetime.datetime.now())[0:19] + ": DTI motion correction was NOT completed; error: "+str(e)) 
                    file.close()
            
            elif SeqName == "MT_OFF_kidneys_cor-oblique_bh":

                try:
                    start_time = time.time()
                    file = open(filename_log, 'a')
                    file.write("\n"+str(datetime.datetime.now())[0:19] + ": MT motion correction has started")
                    file.close()

                    MT_OFF = series
                    for series in folder.series():
                            if series['SeriesDescription'] == "MT_ON_kidneys_cor-oblique_bh":
                                MT_ON = series
                                break
                    mdr.MDRegMT([MT_OFF, MT_ON], study=study)

                    file = open(filename_log, 'a')
                    file.write("\n"+str(datetime.datetime.now())[0:19] + ": MT motion correction was completed --- %s seconds ---" % (int(time.time() - start_time))) 
                    file.close()   

                except Exception as e: 
                    file = open(filename_log, 'a')
                    file.write("\n"+str(datetime.datetime.now())[0:19] + ": MT motion correction was NOT completed; error: "+str(e)) 
                    file.close()

            elif SeqName == "DCE_kidneys_cor-oblique_fb":
                try:
                    start_time = time.time()
                    file = open(filename_log, 'a')
                    file.write("\n"+str(datetime.datetime.now())[0:19] + ": DCE motion correction has started")
                    file.close()

                    mdr.MDRegDCE(series, study=study)

                    file = open(filename_log, 'a')
                    file.write("\n"+str(datetime.datetime.now())[0:19] + ": DCE motion correction was completed --- %s seconds ---" % (int(time.time() - start_time))) 
                    file.close()   

                except Exception as e: 
                    file = open(filename_log, 'a')
                    file.write("\n"+str(datetime.datetime.now())[0:19] + ": DCE motion correction was NOT completed; error: "+str(e)) 
                    file.close()

    folder.save()


