"""
Author: Joao Periquito
Year: 2022

iBEAt Analysis - MAIN EXTERNAL LOCAL
Description:
    For users outside of the University of Sheffield (no access to the UoS Google Drive)

SETUP INSTRUCTIONS:

Step 1. Ensure you have access to Sheffield XNAT (https://qib.shef.ac.uk/app/template/Login.vm#!)
Step 2. Download and zip a complete study for a single iBEAt subject from XNAT (ONLY SIEMENS SITES for now: Bordeaux, Exeter and Leeds)
Step 3. Install the relevant python libraries:

    For Windows:
        1. Create a virtual environment:
           py -3 -m venv .venv
        2. Activate the virtual environment:
           .venv/Scripts/activate
        3. Install python libraries
           pip install -r requirements.txt

    For Mac OSX:
        1. Create a virtual environment:
           python3 -m venv .venv_ibeat
        2. Activate the virtual environment:
           source .venv_ibeat/bin/activate
        3. Install python libraries
           pip install -r requirements.txt

Step 4. Run this script 
        - warning, the computation is intended to run on a cloud or HPC, on a standard laptop/computer this can take several days to complete 
        - for an intemediate update on the progress of the computation feel free to open the logfile during the computation
        - if you interrup the computation at any point all steps completed until that point will produce valid results, 
        - you can resume an interrupted computation by comenting out the steps that have been completed in single_subject_external_analysis.py and run this script again

                   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
OUTPUT:            !!!!!!!!!!!!!!! WARNING: downloaded folder WILL BE MODIFIED by the analysis !!!!!!!!!!!!!!!!
                   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    1. In the original study folder: 
        - Original DICOM files.
        - Newly created DICOM files with calculated maps, segmentations, motion correction...etc 
        - To view this output use any DICOM viewer

    2. New output folder: <original study folder name> + "_output":
        - Created in the same folder as the selected study folder
        - Includes .csv file containing extracted biomarkers (KEY RESULT).
        - Includes images for quality control in .png or .gif format (motion correction results, masks, maps, alignments, ROI models).
    
    3. Log file (.txt):
        - Highlights any error occured during the analysis (for each step)
        - Shows the date and duration of each step


ANALYSIS STAGES:                                               ->     Pipeline Modules:
    1. Load a previously downloaded XNAT dataset.              ->         utilities.select_save_folder.py : Select the save folder for output files.
    2. Standardize naming conventions.                         ->         pipelines.rename.py             : Standardize file and dataset naming conventions.
    3. Execute MDR (motion correction).                        ->         pipelines.mdr.py                : Perform motion correction.
    4. Apply UNETR for kidney segmentation (whole kidney).     ->         pipelines.segment.py            : Apply AI-based whole kidney segmentation (UNETR).
    5. Perform custom modeling (e.g., T1, T2, etc.).           ->         pipelines.mapping.py            : Perform custom modeling (T1, T2, etc.).
    6. Align all scans with Dixon.                             ->         pipelines.align.py              : Align scans with Dixon.
    7. Extract biomarkers.                                     ->         pipelines.measure.py            : Extract ROI statistics and biomarkers.



"""

from scripts.single_subject_external_analysis import single_subject
import utilities.select_folder_to_save as select_save_folder

if __name__ == '__main__':

    path = select_save_folder.external()
    
    single_subject(path)