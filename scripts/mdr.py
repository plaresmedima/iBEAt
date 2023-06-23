""" 
@author: Joao Periquito 
iBEAt MDR Scrpit
2022
Find iBEAt standard pulse sequence name  and execute MDR for: DCE, DTI, T1, T2, T2*, MT
"""

import time

import pipelines.mdr as mdr


def main(folder):
    start_time = time.time()
    folder.log("MDR has started!")

    current_study = folder.series()[0].parent()
    study = folder.series()[0].new_pibling(StudyDescription=current_study.StudyDescription + '_MDRresults')

    for series in folder.series():

        SeqName = series["SeriesDescription"]

        if SeqName is not None:
            print(SeqName)

            if SeqName == "T2star_map_kidneys_cor-oblique_mbh_magnitude":
                try:
                    mdr.MDRegT2star(series, study=study)
                except Exception as e: 
                    folder.log("T2* motion correction was NOT completed; error: "+str(e))

            elif SeqName == "T1map_kidneys_cor-oblique_mbh_magnitude":
                try:
                    mdr.MDRegT1(series, study=study)
                except Exception as e: 
                    folder.log("T1 motion correction was NOT completed; error: "+str(e))

            elif SeqName == "T2map_kidneys_cor-oblique_mbh_magnitude":
                try:
                    mdr.MDRegT2(series, study=study)
                except Exception as e: 
                    folder.log("T2 motion correction was NOT completed; error: "+str(e))   

            elif SeqName == "DTI_kidneys_cor-oblique_fb":
                try:
                    mdr.MDRegDTI(series, study=study) 
                except Exception as e: 
                    folder.log("DTI motion correction was NOT completed; error: "+str(e))
            
            elif SeqName == "MT_OFF_kidneys_cor-oblique_bh":

                try:
                    MT_OFF = series
                    for series in folder.series():
                        if series['SeriesDescription'] == "MT_ON_kidneys_cor-oblique_bh":
                            MT_ON = series
                            break
                    mdr.MDRegMT([MT_OFF, MT_ON], study=study) 
                except Exception as e: 
                    folder.log("MT motion correction was NOT completed; error: "+str(e))

            elif SeqName == "DCE_kidneys_cor-oblique_fb":
                try:
                    mdr.MDRegDCE(series, study=study)   
                except Exception as e: 
                    folder.log("DCE motion correction was NOT completed; error: "+str(e))

    folder.save()
    folder.log("MDR was completed --- %s seconds ---" % (int(time.time() - start_time)))


