import os

from pydrive.auth import GoogleAuth
from pydrive.drive import GoogleDrive

import utilities.download_folder_from_drive as download_drive


def find_mask_in_drive(database):

    gauth = GoogleAuth()
    gauth.LocalWebserverAuth() 
    drive = GoogleDrive (gauth)

    file_list = drive.ListFile({'q': f"'{'1AiY11Sn2pqj4lPGzPTx2B3zop3w0HD1v'}' in parents and trashed=false"}).GetList()
    # for file in file_list:
    #     print(f"Title: {file['title']}, ID: {file['id']}, Type: {file['mimeType']}")

    #folder.FetchContent()

    series_list = database.series()
    subject_name = series_list[0]['PatientID']
    LK_name = subject_name+'_LK'
    RK_name = subject_name+'_RK'

    download_path = database.path()

    for item in file_list:
        #print(f"Title: {item['title']}, ID: {item['id']}, Type: {item['mimeType']}")

        if item['title'] == LK_name:
            print(item['id'])
            download_drive.download_folder_from_drive(item['id'], download_path)
            
            try:
                database.import_dicom(os.path.join(download_path,LK_name))
                database.log("LK Mask was found in the mask drive folder")
            except:
                database.log("LK Mask was NOT found in the mask drive folder")
            database.save()

        if item['title'] == RK_name:
            download_drive.download_folder_from_drive(item['id'], download_path)
            try:
                database.import_dicom(os.path.join(download_path,RK_name))
                database.log("RK Mask was found in the mask drive folder")
            except:
                database.log("RK Mask was NOT found in the mask drive folder")
            database.save()
    return 1


# Dummy placeholder function
def kidney_masks(folder):
    # Check for masks on google drive and import into the database if exists.
    # Alternative: save them on XNAT on a separate project and download from there
    pass

# Dummy placeholder function
def dl_models(folder):
    # Download DL models from google drive to the local folder
    pass


