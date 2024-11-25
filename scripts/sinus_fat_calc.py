"""
Batch script to overlay right- and left kidney masks on an fat of phase DIXON image.

The result is saved in a new directory along with a copy of the fat of phase images.
"""
import os
import time
import datetime
import dbdicom as db
from dbdicom.extensions import skimage, scipy, dipy
import csv
import pandas as pd

# Define import and export directories
data_dir = 'C:\\Users\\md1jdsp\\Desktop\\'
source_dir = os.path.join(data_dir, 'sinus_sample')
results_dir = os.path.join(data_dir, 'output_sinus_fat')

# Series descriptions of all data
out_desc        = 'T1w_abdomen_dixon_cor_bh_fat_post_contrast'

left_mask_desc  = 'LK'
right_mask_desc = 'RK'

folder = db.database(path=source_dir)

filename_log = results_dir + '//' + datetime.datetime.now().strftime('%Y%m%d_%H%M_') + "SinusFat_LogFile.txt"

for patient in folder.patients():

    df = pd.DataFrame([])

    file = open(filename_log, 'a')
    file.write("\n")
    file.write("\n"+str(datetime.datetime.now())[0:19] + ": Doing patient " + patient.uid)
    file.close()

    fat = patient.series(SeriesDescription=out_desc)

    if fat != []:
        file = open(filename_log, 'a')
        file.write("\n"+str(datetime.datetime.now())[0:19] + ": Fat has been Found")
        file.close()

    lk = patient.series(SeriesDescription='LK')
    rk = patient.series(SeriesDescription='RK')

    if (lk !=[] and rk !=[]):
        file = open(filename_log, 'a')
        file.write("\n"+str(datetime.datetime.now())[0:19] + ": LK mask has been Found: " + lk[0].SeriesDescription)
        file.write("\n"+str(datetime.datetime.now())[0:19] + ": RK mask has been Found: " + rk[0].SeriesDescription)
        file.close()
    else:
        lk = patient.series(SeriesDescription='LK [overlay]')
        rk = patient.series(SeriesDescription='LK [overlay]')

    if (lk !=[] and rk !=[]):
        file = open(filename_log, 'a')
        file.write("\n"+str(datetime.datetime.now())[0:19] + ": LK mask has been Found: " + lk[0].SeriesDescription)
        file.write("\n"+str(datetime.datetime.now())[0:19] + ": RK mask has been Found: " + rk[0].SeriesDescription)
        file.close()
    
    if (lk !=[] and rk ==[]):
        file = open(filename_log, 'a')
        file.write("\n"+str(datetime.datetime.now())[0:19] + ": ONLY LK mask has been Found: " + lk[0].SeriesDescription)
        file.close()
    
    if (lk ==[] and rk !=[]):
        file = open(filename_log, 'a')
        file.write("\n"+str(datetime.datetime.now())[0:19] + ": ONLY RK mask has been Found: " + rk[0].SeriesDescription)
        file.close()

    if (lk ==[] or rk ==[]):
        continue

    kidneys = lk+rk

    cleanup = True

    sf_series = []

    try:
        fat_image_masked, fat_mask = dipy.median_otsu(fat[0], median_radius=1, numpass=1)
    except:
        file = open(filename_log, 'a')
        file.write("\n"+str(datetime.datetime.now())[0:19] + ": Sinus Fat has NOT been calculated: ")
        file.close()
        continue

    for kidney in kidneys:
        # Pipeline calculation
        kidney_hull = skimage.convex_hull_image_3d(kidney)
        # sinus = scipy.image_calculator(kidney_hull, kidney, 'series 1 - series 2', integer=True)
        # sinus_fat = scipy.image_calculator(fat_mask, sinus, 'series 1 * series 2', integer=True)
        sinus_fat = scipy.image_calculator(fat_mask, kidney_hull, 'series 1 * series 2', integer=True)
        sinus_fat_largest = scipy.extract_largest_cluster_3d(sinus_fat)
        sinus_fat_largest.SeriesDescription = kidney.instance().SeriesDescription + 'SF'
        # Append and display
        sf_series.append(sinus_fat_largest)
                # Remove intermediate results
        if cleanup:
            kidney_hull.remove()
            #sinus.remove()
            sinus_fat.remove()
    
    fat_image_masked.remove()
    if cleanup:
        fat_mask.remove()

    if sf_series !=[] :
        file = open(filename_log, 'a')
        file.write("\n"+str(datetime.datetime.now())[0:19] + ": Sinus Fat has been calculated")
        file.close()
    else:
        file = open(filename_log, 'a')
        file.write("\n"+str(datetime.datetime.now())[0:19] + ": Sinus Fat has NOT been calculated")
        file.close()
        continue
    try:
        df = pd.DataFrame(skimage.volume_features(sf_series))
    except:
        continue

    if (df.empty == False):
         
        df.to_csv(results_dir + '//' + 'sinus_fat.csv', header=None, mode='a')
        
        file = open(filename_log, 'a')
        file.write("\n"+str(datetime.datetime.now())[0:19] + ": Metrics have been collected and saved in CSV ")
        file.close()
    else:
        file = open(filename_log, 'a')
        file.write("\n"+str(datetime.datetime.now())[0:19] + ": Metrics have NOT been collected ")
        file.close()
        continue

    for series in sf_series:
        series.export_as_dicom(results_dir)

        file = open(filename_log, 'a')
        file.write("\n"+str(datetime.datetime.now())[0:19] + ": " + series.SeriesDescription + "  have been exported")
        file.close()

    
    





