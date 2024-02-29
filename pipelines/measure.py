import os
import pandas as pd

from dbdicom.extensions import skimage

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
    row_headers = ['PatientID', 'SeriesDescription', 'Region of Interest', 'Parameter', 'Value', 'Unit', 'Biomarker']
    master_table = pd.DataFrame(columns=row_headers)
    master_table.to_csv(filename_csv, index=False)
    return master_table


def _update_master_table(folder, table):
    master_table = _load_master_table(folder)
    # This needs modifying: if the values are already in the table, then update instead of adding
    master_table = pd.concat([master_table, table], ignore_index=True)
    master_table = master_table.drop_duplicates(subset='Biomarker', keep='last', ignore_index=True)
    filename_csv = _master_table_file(folder)
    master_table.to_csv(filename_csv, index=False)


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