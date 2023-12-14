import pandas as pd
from dbdicom.wrappers import skimage, scipy, dipy


def kidney_volumetrics(master_table, folder):
    
    left  = folder.series(SeriesDescription='LK')
    right = folder.series(SeriesDescription='RK')

    kidneys = [left, right]
    features = skimage.volume_features(kidneys)

    master_table = pd.concat([master_table, features],ignore_index=True)

    return master_table

def sinus_fat_volumetrics(master_table,folder):

    out_desc        = 'T1w_abdomen_dixon_cor_bh_fat_post_contrast'

    fat = folder.series(SeriesDescription=out_desc)
    lk  = folder.series(SeriesDescription='LK')
    rk  = folder.series(SeriesDescription='RK')

    kidneys = lk+rk
    cleanup = True
    sf_series = []

    fat_image_masked, fat_mask = dipy.median_otsu(fat[0], median_radius=1, numpass=1)

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

    features = pd.DataFrame(skimage.volume_features(sf_series))

    master_table = pd.concat([master_table, features],ignore_index=True)

    return master_table

def main(master_table, folder):

    master_table = kidney_volumetrics(master_table, folder)
    
    master_table = sinus_fat_volumetrics(master_table,folder)

    return master_table
