import dbdicom as db
from pipelines import mapping
from utilities import map_to_dixon_mask
from utilities import fill_kidney_holes_interp_v2
import numpy as np


if __name__ == '__main__':
    pathScan = "C://Users//md1jdsp//Desktop//iBE-4128-055"
    #pathScan = "C://Users//md1jdsp//Desktop//iBE-3128-027"

    folder = db.database(path=pathScan)

    for series in folder.series():
        SeqName = series["SeriesDescription"]
        print(SeqName)
        if SeqName == "T1map_kidneys_cor-oblique_mbh_magnitude_mdr_moco":
        # if SeqName == "T1map_kidneys_cor-oblique_mbh_magnitude_mdr_par_T1":
            T1_map = series
        #if SeqName == "T2map_kidneys_cor-oblique_mbh_magnitude_mdr_par":
        #if SeqName == "T2map_kidneys_cor-oblique_mbh_magnitude_mdr_par_T2":
        if SeqName == "T2map_kidneys_cor-oblique_mbh_magnitude_mdr_moco": 
            T2_map = series
        if SeqName == "LK":
            Kidney_mask = series
            array_Mask, header_Mask = series.array(['SliceLocation'], pixels_first=True)

    mapping.T1_then_T2([T1_map,T2_map], Kidney_mask)
    folder.save()

    # params = map_to_dixon_mask.find_sbs_rigid_transformation(T2_map, Kidney_mask, tolerance=0.1, metric='mutual information', region=Kidney_mask,resolutions=[4,2,1])
    # T2_map_coreg = map_to_dixon_mask.apply_sbs_rigid_transformation(T2_map, params, target =Kidney_mask)

    # array_Data, header_Data = T2_map_coreg.array(['SliceLocation'], pixels_first=True)
    # array_Corrected = fill_kidney_holes_interp_v2.main(array_Data, array_Mask)

    # study = T2_map_coreg.parent()
    # Corrected_map_series = series.SeriesDescription + "_T2map_" + "filled"
    # M0_map_series = study.new_series(SeriesDescription=Corrected_map_series)
    # M0_map_series.set_array(array_Corrected,np.squeeze(header_Data[:,0]),pixels_first=True)



    # T2w_imgs_coreg_series = T1_map_to_dixon_mask.main(T2_map,Kidney_mask)
    # T2w_imgs_coreg_array, T2w_imgs_coreg_header = T2w_imgs_coreg_series.array(['SliceLocation', 'AcquisitionTime'], pixels_first=True)
    # T2w_imgs_coreg_array = np.squeeze(T2w_imgs_coreg_array)


