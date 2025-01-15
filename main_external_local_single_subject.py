"""
Author: Joao Periquito
Year: 2022

iBEAt Analysis - MAIN EXTERNAL LOCAL
Description:
    For users outside of the University of Sheffield (no access to the UoS Google Drive).

Workflow:                                                      ->     Pipeline Modules:
    1. Load a previously downloaded XNAT dataset.              ->         utilities.select_save_folder.py : Select the save folder for output files.
    2. Standardize naming conventions.                         ->         pipelines.rename.py             : Standardize file and dataset naming conventions.
    3. Execute MDR (motion correction).                        ->         pipelines.mdr.py                : Perform motion correction.
    4. Apply UNETR for kidney segmentation (whole kidney).     ->         pipelines.segment               : Apply AI-based whole kidney segmentation (UNETR).
    5. Perform custom modeling (e.g., T1, T2, etc.).           ->         pipelines.mapping.py            : Perform custom modeling (T1, T2, etc.).
    6. Align all scans with Dixon.                             ->         pipelines.align.py              : Align scans with Dixon.
    7. Extract biomarkers.                                     ->         pipelines.measure               : Extract ROI statistics and biomarkers.

INPUT:
    - Unzipped, downloaded XNAT dataset.

OUTPUT:
    1. Output folder:
        - CSV file containing extracted biomarkers.
        - QC images (motion correction results, masks, maps, alignments, ROI models).
        - Editable canvas for masks.
    2. DICOM folder:
        - Original DICOM files.
        - Newly created DICOM files.

SETUP INSTRUCTIONS:
    Windows:
        1. Create a virtual environment:
           py -3 -m venv .venv
        2. Activate the virtual environment:
           .venv/Scripts/activate
        3. Install python libraries
           pip install -r requirements.txt
        4. Download and unzip an iBEAt dataset from XNAT
        5. Run this script

    Mac OSX:
        1. Create a virtual environment:
           python3 -m venv .venv_ibeat
        2. Activate the virtual environment:
           source .venv_ibeat/bin/activate
        3. Install python libraries
           pip install -r requirements.txt
        4. Download and unzip an iBEAt dataset from XNAT
        5. Run this script
"""

from scripts.single_subject_external_analysis import single_subject
import utilities.select_folder_to_save as select_save_folder

if __name__ == '__main__':

    path = select_save_folder.external()
    
    single_subject(path)