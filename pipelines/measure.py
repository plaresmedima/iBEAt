import os
import pandas as pd
import numpy as np

from dbdicom.extensions import skimage, vreg


def _master_table_file(folder):
    # Get filename and create folder if needed.
    results_path = folder.path() + '_output'
    if not os.path.exists(results_path):
        os.mkdir(results_path)
    return os.path.join(results_path, 'biomarkers.csv')

def _load_master_table(folder):
    # Read master table
    filename_csv = _master_table_file(folder)
    if os.path.isfile(filename_csv):
        return pd.read_csv(filename_csv)
    # If the master table does not exist, create it
    row_headers = ['PatientID', 'SeriesDescription', 'Region of Interest', 'Parameter', 'Value', 'Unit', 'Biomarker', 'StudyDescription']
    master_table = pd.DataFrame(columns=row_headers)
    master_table.to_csv(filename_csv, index=False)
    return master_table

def read_master_table(folder, biomarker):
    table = _load_master_table(folder)
    table = table[table.Biomarker == biomarker]
    if table.empty:
        raise ValueError('Biomarker ' + biomarker +' has not yet been calculated.')
    return table.Value.values[0]

def _update_master_table(folder, table):
    master_table = _load_master_table(folder)
    master_table = pd.concat([master_table, table], ignore_index=True)
    master_table = master_table.drop_duplicates(subset='Biomarker', keep='last', ignore_index=True)
    filename_csv = _master_table_file(folder)
    master_table.to_csv(filename_csv, index=False)

def add_rows(folder, rows):
    row_headers = ['PatientID', 'SeriesDescription', 'Region of Interest', 'Parameter', 'Value', 'Unit', 'Biomarker', 'StudyDescription']
    table = pd.DataFrame(data=rows, columns=row_headers)
    _update_master_table(folder, table)


def kidney_volumetrics(folder):
    left  = folder.series(SeriesDescription='LK')
    right = folder.series(SeriesDescription='RK')
    kidneys = [left, right]
    features = skimage.volume_features(kidneys)
    features['Biomarker'] = features['SeriesDescription'] + '-' + features['Parameter']
    features['Region of Interest'] = 'Kidney'
    _update_master_table(folder, features)
    return features
    

def sinus_fat_volumetrics(folder):
    left  = folder.series(SeriesDescription='LKSF')
    right = folder.series(SeriesDescription='RKSF')
    sinus_fat = [left, right]
    features = skimage.volume_features(sinus_fat)
    features['Biomarker'] = features['SeriesDescription'] + '-' + features['Parameter']
    features['Region of Interest'] = 'Renal Sinus Fat'
    _update_master_table(folder, features)
    return features


def t1_maps(folder):
    seq = 'T1map_kidneys_cor-oblique_mbh_magnitude_mdr_moco'
    pars = ['T1', 'FA']
    units = ['msec', 'deg']
    return features(folder, seq, pars, units)

def t2_maps(folder):
    seq = 'T2map_kidneys_cor-oblique_mbh_magnitude_mdr_moco'
    pars = ['T2']
    units = ['msec']
    return features(folder, seq, pars, units)

def t2star_maps(folder):
    seq = 'T2star_map_kidneys_cor-oblique_mbh_magnitude_mdr_moco'
    pars = [ 'T2s','fw']
    units = ['','']
    return features(folder, seq, pars, units)

def mt_maps(folder):
    seq = 'MT_kidneys_cor-oblique_bh_mdr_moco'
    pars = ['MTR']
    units = ['']
    return features(folder, seq, pars, units)

def ivim_maps(folder):
    seq = 'IVIM_kidneys_cor-oblique_fb_mdr_moco'
    pars = ['MD','Df','ff']
    units = ['','','']
    return features(folder, seq, pars, units)

def dti_maps(folder):
    seq = 'DTI_kidneys_cor-oblique_fb_mdr_moco'
    pars = ['MD','RD','AD','Planarity','Linearity','Sphericity','FA']
    units = ['','','','','','','']
    return features(folder, seq, pars, units)

def dce_maps(folder):
    seq = 'DCE_kidneys_cor-oblique_fb_mdr_moco'
    pars = ['AVD', 'RPF', 'MTT']
    units = ['mL/100mL', 'mL/min/100mL', 'sec']
    return features(folder, seq, pars, units)

# Needs to be brought in line with the others - same format
def asl_maps(folder):
    for ROI in ['LK','RK']:
        kidney  = folder.series(SeriesDescription=ROI)
        rbf = folder.series(SeriesDescription='RBF - ' + ROI)
        features = vreg.mask_statistics(kidney, rbf)
        features['Biomarker'] = features['SeriesDescription'] + '-' + features['Parameter']
        features['Region of Interest'] = 'Kidney'
        _update_master_table(folder, features)
    return features


def features(folder, seq, pars, units):
    for p, par in enumerate(pars):
        vals = []
        for ROI in ['LK','RK']: 
            folder.message('Exporting metrics for parameter ' + par)
            desc = seq + '_' + par + '_map_' + ROI + '_align'
            kidney = folder.series(SeriesDescription=ROI)[0]
            series = folder.series(SeriesDescription=desc)[0]
            kidney_vals = vreg.mask_values(kidney, series)
            vals.append(kidney_vals)

            # update master table
            features = vreg.mask_data_statistics(kidney_vals, kidney, series)
            features['Biomarker'] = [ROI + '-' + par + '-' + metric for metric in features['Parameter'].values]
            features['Unit'] = units[p]
            features['SeriesDescription'] = ROI
            features['Region of Interest'] = 'Kidney'
            _update_master_table(folder, features)

        ROI = 'BK' 
        vals = np.concatenate(vals)
        features = vreg.mask_data_statistics(vals, kidney, series)
        features['Biomarker'] = [ROI + '-' + par + '-' + metric for metric in features['Parameter'].values]
        features['Unit'] = units[p]
        features['SeriesDescription'] = ROI
        features['Region of Interest'] = 'Kidney'
        _update_master_table(folder, features)

    return features


