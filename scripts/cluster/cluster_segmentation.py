""" 
@author: Joao Periquito 
iBEAt SEGMENTATION Scrpit
2023
Find iBEAt fat and opposed phase scan and execute kmeans segmentation
"""
from dbdicom.wrappers import sklearn, skimage, scipy, dipy
import datetime
import time
import os

def export_segmentation(folder,path):

    resultsFolder = 'segmentation_results'
    resultsPath = path + '_' + resultsFolder

    os.mkdir(resultsPath)

    fat_desc        = 'T1w_abdomen_dixon_cor_bh_fat_post_contrast' 
    out_desc        = 'T1w_abdomen_dixon_cor_bh_out_phase_post_contrast'
    in_desc         = 'T1w_abdomen_dixon_cor_bh_in_phase_post_contrast'
    water_desc      = 'T1w_abdomen_dixon_cor_bh_water_post_contrast'
    k_means1_desc   = 'KMeans cluster 1'
    k_means2_desc   = 'KMeans cluster 2'

    fat             = folder.series(SeriesDescription=fat_desc)
    out_ph          = folder.series(SeriesDescription=out_desc)
    in_ph           = folder.series(SeriesDescription=in_desc)
    water           = folder.series(SeriesDescription=water_desc)
    k_means1        = folder.series(SeriesDescription=k_means1_desc)
    k_means2        = folder.series(SeriesDescription=k_means2_desc)

    exportToFolder = fat + out_ph + in_ph + water + k_means1 + k_means2
    
    for series in exportToFolder:
        print(series.SeriesDescription)    
        series.export_as_dicom(resultsPath)
    
    print('done')


def main(folder,path,filename_log):

    fat_desc = 'T1w_abdomen_dixon_cor_bh_fat_post_contrast' 
    out_desc = 'T1w_abdomen_dixon_cor_bh_out_phase_post_contrast'

    fat = folder.series(SeriesDescription=fat_desc)
    out = folder.series(SeriesDescription=out_desc)
    msg = 'OK'
    if fat==[] or out==[]:
        msg = 'Not all sequences have been found - please check the list'
    elif len(fat)>1 or len(out)>1:
        msg = 'Some sequences have been repeated - please select the appropriate for analysis'
    if msg == 'OK':
        features = fat + out

    try:
        start_time = time.time()
        file = open(filename_log, 'a')
        file.write("\n"+str(datetime.datetime.now())[0:19] + ": Sequential kmeans has started")
        file.close()
        
        sklearn.sequential_kmeans(features, n_clusters=2, multiple_series=True)
        folder.save()

        file = open(filename_log, 'a')
        file.write("\n"+str(datetime.datetime.now())[0:19] + ": Sequential kmeans was completed --- %s seconds ---" % (int(time.time() - start_time))) 
        file.close() 
    
    except Exception as e: 

        file = open(filename_log, 'a')
        file.write("\n"+str(datetime.datetime.now())[0:19] + ": Sequential kmeans was NOT completed; error: "+str(e)) 
        file.close()

    try:
        start_time = time.time()
        file = open(filename_log, 'a')
        file.write("\n"+str(datetime.datetime.now())[0:19] + ": Export segmentation results has started")
        file.close()
        
        export_segmentation(folder,path)

        file = open(filename_log, 'a')
        file.write("\n"+str(datetime.datetime.now())[0:19] + ": Export segmentation results was completed --- %s seconds ---" % (int(time.time() - start_time))) 
        file.close() 
    
    except Exception as e: 

        file = open(filename_log, 'a')
        file.write("\n"+str(datetime.datetime.now())[0:19] + ": Export segmentation results was NOT completed; error: "+str(e)) 
        file.close()

    
