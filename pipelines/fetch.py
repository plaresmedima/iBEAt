import os
import shutil

from pydrive.auth import GoogleAuth
from pydrive.drive import GoogleDrive

import utilities.download_folder_from_drive as download_drive
import utilities.zenodo_link as UNETR_zenodo
import requests


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
                database.log("LK Mask was found in the mask drive folder")
            except:
                database.log("LK Mask was NOT found in the mask drive folder")
            database.save()

        if item['title'] == RK_name:
            download_drive.download_folder_from_drive(item['id'], download_path)
            RK_drive = 1
            try:
                files_RK = [os.path.join(os.path.join(download_path,RK_name,'dbdicom'), file) for file in os.listdir(os.path.join(download_path,RK_name,'dbdicom'))]
                database.import_dicom(files_RK)
                shutil.rmtree(os.path.join(download_path,RK_name))
                database.log("RK Mask was found in the mask drive folder")
            except:
                database.log("RK Mask was NOT found in the mask drive folder")
            database.save()
    
    if LK_drive == 1 and RK_drive == 1:
         database.log("Both masks were found in the mask drive folder")
         
# Dummy placeholder function
def kidney_masks(folder):
    # Check for masks on google drive and import into the database if exists.
    # Alternative: save them on XNAT on a separate project and download from there

    find_mask_in_drive(folder)

def dl_models(database):

    unetr, unetr_link= UNETR_zenodo.main()
    record_id = unetr_link.split('.')[-1]

    unetr_path = os.path.join(os.path.dirname(database.path()),unetr)

    if os.path.exists(unetr_path):
        database.log("UNETR was found in the local folder")
        return
    else:   

        zenodo_url = f"https://zenodo.org/records/{record_id}/files/{unetr}?download=1"

        with requests.get(zenodo_url) as req:
                    with open(os.path.join(os.path.dirname(database.path()),unetr), 'wb') as f:
                        for chunk in req.iter_content(chunk_size=8192):
                            if chunk:
                                f.write(chunk)





