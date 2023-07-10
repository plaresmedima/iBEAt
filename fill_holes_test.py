import dbdicom as db
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import distance_transform_edt
from scipy.ndimage import map_coordinates
from utilities import fill_kidney_holes

pathScan = "C://Users//md1jdsp//Desktop//iBE-2128_013_baseline"

folder = db.database(path=pathScan)

for series in folder.series():
    SeqName = series["SeriesDescription"]
    print(SeqName)

    if SeqName == "T1map_kidneys_cor-oblique_mbh_magnitude_mdr_par_T1 [sbs rigid]":
        array_Data, header_Data = series.array(['SliceLocation'], pixels_first=True)
    if SeqName == "LK [overlay]":
        array_Mask, header_Mask = series.array(['SliceLocation'], pixels_first=True)

array_Data = np.squeeze(array_Data)
array_Mask = np.squeeze(array_Mask)

array_Corrected = fill_kidney_holes.main(array_Data, array_Mask)

plt.subplot(2,1,1)
plt.imshow(np.transpose(array_Data[:,:,90]),vmin=0, vmax=2000)
plt.subplot(2,1,2)
plt.imshow(np.transpose(array_Corrected[:,:,90]),vmin=0, vmax=2000)
plt.show()