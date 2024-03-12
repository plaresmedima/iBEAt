import os
import os.path

import numpy as np
import cv2
from skimage.segmentation import flood

from dbdicom.extensions import skimage, scipy, dipy, sklearn
from dbdicom.pipelines import input_series
import models.UNETR_kidneys_v1 as unetr



export_study = 'Segmentations'


def kidneys(database, weights):

    # Get weights file and check if valid 
    # if not os.path.isfile(weights):
    #     msg = 'The weights file ' + weights + ' has not been found. \n'
    #     msg += 'Please check that the file with model weights is in the folder, and is named ' + unetr.filename
    #     database.dialog.information(msg)
    #     return

    database.message('Segmenting kidneys. This could take a few minutes. Please be patient..')

    # Get appropriate series and check if valid
    #series = database.series(SeriesDescription=unetr.trained_on)
    sery, study = input_series(database, unetr.trained_on, export_study)
    if sery is None:
        return    

    # Loop over the series and create the mask
    #desc = sery.instance().SeriesDescription
    array, header = sery.array(['SliceLocation'], pixels_first=True, first_volume=True)

    # Calculate predictions 
    masks = unetr.apply(array, weights)
    left_kidney, right_kidney = unetr.kidney_masks(masks)

    # Save UNETR output
    result = study.new_child(SeriesDescription = 'BK')
    result.set_array(masks, header, pixels_first=True)
    # result[['WindowCenter','WindowWidth']] = [1.0, 2.0]

    # Save and display left kidney data
    left = study.new_child(SeriesDescription = 'LK')
    left.set_array(left_kidney, header, pixels_first=True)
    # left[['WindowCenter','WindowWidth']] = [1.0, 2.0]
    
    # Save and display right kidney data
    right = study.new_child(SeriesDescription = 'RK')
    right.set_array(right_kidney, header, pixels_first=True)
    # right[['WindowCenter','WindowWidth']] = [1.0, 2.0]

    return left, right


def renal_sinus_fat(folder):

    fat = folder.series(SeriesDescription='T1w_abdomen_dixon_cor_bh_fat_post_contrast')
    lk  = folder.series(SeriesDescription='LK')
    rk  = folder.series(SeriesDescription='RK')

    kidneys = lk+rk
    sf_series = []

    if len(kidneys)==[]:
        msg = 'Cannot perform renal sinus fat segmentation: no kidney masks are available.'
        raise RuntimeError(msg)

    fat_image_masked, fat_mask = dipy.median_otsu(fat[0], median_radius=1, numpass=1)

    for kidney in kidneys:
        kidney_hull = skimage.convex_hull_image_3d(kidney)
        sinus_fat = scipy.image_calculator(fat_mask, kidney_hull, 'series 1 * series 2', integer=True)
        #sinus_fat_open = skimage.opening_3d(sinus_fat)
        sinus_fat_largest = scipy.extract_largest_cluster_3d(sinus_fat)
        sinus_fat_largest.SeriesDescription = kidney.instance().SeriesDescription + 'SF'
        sf_series.append(sinus_fat_largest)
        # Cleanup
        kidney_hull.remove()
        #sinus.remove()
        sinus_fat.remove()
    
    fat_image_masked.remove()
    fat_mask.remove()   

    return sf_series


def compute_whole_kidney_canvas(database):
    series_desc = [
        'T1w_abdomen_dixon_cor_bh_fat_post_contrast',
        'T1w_abdomen_dixon_cor_bh_out_phase_post_contrast'
    ] 
    features, study = input_series(database, series_desc, export_study)
    if features is None:
        return
    clusters = sklearn.sequential_kmeans(features, n_clusters=2, multiple_series=True)
    for c in clusters:
        c.move_to(study)
    return clusters




#### THIS NEEDS CLEANING UP

class Point(object):
    def __init__(self,x,y):
        self.x = x
        self.y = y

    def getX(self):
        return self.x
    def getY(self):
        return self.y

def getGrayDiff(img,currentPoint,tmpPoint):
    return abs(int(img[currentPoint.x,currentPoint.y]) - int(img[tmpPoint.x,tmpPoint.y]))

def selectConnects(p):
 if p != 0:
  connects = [Point(-1, -1), Point(0, -1), Point(1, -1), Point(1, 0), Point(1, 1), \
     Point(0, 1), Point(-1, 1), Point(-1, 0)]
 else:
  connects = [ Point(0, -1), Point(1, 0),Point(0, 1), Point(-1, 0)]
 return connects

def regionGrow(img,seeds,thresh,p = 1):
 height, weight = img.shape
 seedMark = np.zeros(img.shape)
 seedList = []
 for seed in seeds:
  seedList.append(seed)
 label = 1
 connects = selectConnects(p)
 while(len(seedList)>0):
  currentPoint = seedList.pop(0)

  seedMark[currentPoint.x,currentPoint.y] = label
  for i in range(8):
   tmpX = currentPoint.x + connects[i].x
   tmpY = currentPoint.y + connects[i].y
   if tmpX < 0 or tmpY < 0 or tmpX >= height or tmpY >= weight:
    continue
   grayDiff = getGrayDiff(img,currentPoint,Point(tmpX,tmpY))
   if grayDiff < thresh and seedMark[tmpX,tmpY] == 0:
    seedMark[tmpX,tmpY] = label
    seedList.append(Point(tmpX,tmpY))
 return seedMark


def aorta_on_dce(folder):

    desc = "DCE_aorta_axial_fb"
    series, study = input_series(folder, desc, export_study)
    if series is None:
        raise RuntimeError('Cannot create DCE-AIF mask: series ' + desc + ' does not exist. ')

    axial, header = series.array(['AcquisitionTime'], pixels_first=True, first_volume=True)

    cutRatio=0.25             #create a window around the center of the image where the aorta is
    filter_kernel=(15,15)     #gaussian kernel for smoothing the image to destroy any noisy single high intensity filter
    threshold = 2     #threshold for the region growing algorithm

    # Calculate max signal enhancement over a window around the center
    cy, cx = int(axial.shape[0]/2), int(axial.shape[1]/2)
    y0, y1 = int(cy-cy*cutRatio), int(cy+cy*cutRatio)
    x0, x1 = int(cx-cx*cutRatio), int(cx+cx*cutRatio)

    axwin = np.zeros(axial.shape)
    axwin[y0:y1, x0:x1,:] = axial[y0:y1, x0:x1,:]
    axenh = np.max(axwin,axis=2) - np.min(axwin,axis=2)

    # Get 3 seed points with maximal enhhancement values
    axenh_blur = cv2.GaussianBlur(axenh, filter_kernel,cv2.BORDER_DEFAULT)
    _, _, _, p1 = cv2.minMaxLoc(axenh_blur)
    axenh_blur[p1[1],p1[0]] = 0
    _, _, _, p2 = cv2.minMaxLoc(axenh_blur)
    axenh_blur[p2[1],p2[0]] = 0
    _, _, _, p3 = cv2.minMaxLoc(axenh_blur)

    seeds = [Point(p1[1],p1[0]), Point(p2[1],p2[0]), Point(p3[1],p3[0])]
    max_iteration = 20
    for i in range(max_iteration):
        aif_mask = regionGrow(axenh,seeds,threshold)
        if len(aif_mask[aif_mask==1]) >= 25:
            break
        threshold += 1

    aif_mask_series = study.new_series(SeriesDescription='DCE-AIF')
    aif_mask_series.set_array(aif_mask, header[0], pixels_first=True)

    return aif_mask_series

