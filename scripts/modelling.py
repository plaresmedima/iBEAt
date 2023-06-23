""" 
@author: Joao Periquito 
iBEAt MODELLING Scrpit
2022
Find iBEAt motion corrected pulse sequence name (MDR output: *_mdr_moco) and execute custom modelling for: DCE, DTI, T2*
(T1 & T2 modelling are done in the main due to parallelization requirements)
"""
import time
import pipelines.mapping as map


def main(folder):

    start_time = time.time()
    folder.log(": Modelling has started!")

    list_of_series = folder.series()

    current_study = list_of_series[0].parent()
    study = list_of_series[0].new_pibling(StudyDescription=current_study.StudyDescription + '_ModellingResults')

    for i,series in enumerate(list_of_series):

        if series["SequenceName"] is not None:

            if series['SeriesDescription'] == "T2star_map_kidneys_cor-oblique_mbh_magnitude_mdr_moco":
                try:
                    map.T2s(series, study=study)
                except Exception as e: 
                    series.log("T2* mapping was NOT completed; error: "+str(e))

            elif series['SeriesDescription'] == "DTI_kidneys_cor-oblique_fb_mdr_moco":
                try:
                    map.DTI(series, study=study)
                except Exception as e: 
                    series.log("DTI-FA & ADC mapping was NOT completed; error: "+str(e))

            elif series['SeriesDescription'] == "DCE_kidneys_cor-oblique_fb_mdr_moco":
                try:
                    map.DCE_MAX(series, study=study)
                except Exception as e: 
                    series.log("DCE-MAX mapping was NOT completed; error: "+str(e))

            elif series.SeriesDescription == 'MT_ON_kidneys_cor-oblique_bh_mdr_moco':
                try:
                    map.MTR(series, study=study)
                except Exception as e: 
                    series.log("MTR mapping was NOT completed; error: "+str(e))

            elif series.SeriesDescription == 'T1_map_kidneys_cor-oblique_mbh_magnitude_mdr_moco':
                if series.Manufacturer == 'SIEMENS':
                    try:
                        map.T1_Philips(series, study=study)  
                    except Exception as e: 
                        series.log("T1 mapping was NOT completed; error: "+str(e))

    folder.save()
    folder.log("Modelling was completed --- %s seconds ---" % (int(time.time() - start_time)))
