import os
import shutil

from pydrive.auth import GoogleAuth
from pydrive.drive import GoogleDrive

import utilities.download_folder_from_drive as download_drive
import utilities.zenodo_link as UNETR_zenodo
import requests



import os

def find_mask_in_local_rep(database):
    # Define the mask folder path
    mask_folder = '//mnt//fastdata//md1jdsp//Leeds_Vol_Masks'
    #mask_folder = 'C://Users//md1jdsp//Desktop//Leeds_Vol_Masks'
    
    # Get the name of the dataset folder from the database path
    dataset_folder_name = os.path.basename(database.path())
    
    # List all folders in the mask directory
    mask_folder_name_list = os.listdir(mask_folder)
    
    # Check if the dataset folder name exists in the mask folder list
    if dataset_folder_name in mask_folder_name_list:
        dataset_folder_path = os.path.join(mask_folder, dataset_folder_name)
        for folder in os.listdir(dataset_folder_path):
            try:
                # Get the path to the 'dbdicom' directory
                dbdicom_path = os.path.join(dataset_folder_path, folder, 'dbdicom')
                
                # Check if 'dbdicom' directory exists
                if os.path.exists(dbdicom_path):
                    # List all files in the 'dbdicom' directory
                    maskfiles = [os.path.join(dbdicom_path, file) for file in os.listdir(dbdicom_path)]
                    
                    # Import DICOM files into the database
                    database.import_dicom(maskfiles)
                    database.log(folder + " was found in the repeatability mask local folder")
                else:
                    database.log(folder + " was NOT found in the repeatability mask local folder")
            except Exception as e:
                database.log(f"Error processing folder {folder}: {str(e)}")
        
        # Save the database after processing all folders
        database.save()
    else:
        database.log(f"{dataset_folder_name} was NOT found in the mask folder list")


def find_mask_in_drive(database):

    gauth = GoogleAuth()
    gauth.LocalWebserverAuth() 
    drive = GoogleDrive (gauth)

    file_list = drive.ListFile({'q': f"'{'1AiY11Sn2pqj4lPGzPTx2B3zop3w0HD1v'}' in parents and trashed=false"}).GetList()
    # for file in file_list:
    #     print(f"Title: {file['title']}, ID: {file['id']}, Type: {file['mimeType']}")

    #folder.FetchContent()

    subject_name = database.PatientName
    LK_name = subject_name+'_LK'
    RK_name = subject_name+'_RK'

    download_path = os.path.dirname(database.path())

    LK_drive = 0
    RK_drive = 0

    for item in file_list:
        #print(f"Title: {item['title']}, ID: {item['id']}, Type: {item['mimeType']}")

        if item['title'] == LK_name:
            download_drive.download_folder_from_drive(item['id'], download_path)
            LK_drive = 1

            try:
                files_LK = [os.path.join(os.path.join(download_path,LK_name,'dbdicom'), file) for file in os.listdir(os.path.join(download_path,LK_name,'dbdicom'))]
                database.import_dicom(files_LK)
                shutil.rmtree(os.path.join(download_path,LK_name))
                database.log(LK_name + " was found in the mask drive folder")
            except:
                database.log(LK_name + " was NOT found in the mask drive folder")
            database.save()

        if item['title'] == RK_name:
            download_drive.download_folder_from_drive(item['id'], download_path)
            RK_drive = 1
            try:
                files_RK = [os.path.join(os.path.join(download_path,RK_name,'dbdicom'), file) for file in os.listdir(os.path.join(download_path,RK_name,'dbdicom'))]
                database.import_dicom(files_RK)
                shutil.rmtree(os.path.join(download_path,RK_name))
                database.log(RK_name + " was found in the mask drive folder")
            except:
                database.log(RK_name + " was NOT found in the mask drive folder")
            database.save()
    
    if LK_drive == 1 and RK_drive == 1:
         database.log("Both masks were found in the mask drive folder")
         
# Dummy placeholder function
def kidney_masks(folder):
    # Check for masks on google drive and import into the database if exists.
    # Alternative: save them on XNAT on a separate project and download from there

    find_mask_in_drive(folder)
    #find_mask_in_local_rep(folder)





