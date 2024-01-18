import os
import os.path
import time
import datetime
import zipfile
import shutil

from pydrive.auth import GoogleAuth
from pydrive.drive import GoogleDrive

def zipdir(path, ziph):
    # ziph is zipfile handle
    for root, dirs, files in os.walk(path):
        for file in files:
            ziph.write(os.path.join(root, file), 
                       os.path.relpath(os.path.join(root, file), 
                                       os.path.join(path, '..')))

def GoogleDrive_Upload(pathScan,filename_log,filename_csv):

    gauth = GoogleAuth()
    drive = GoogleDrive (gauth)

    #Export results csv
    gfile = drive.CreateFile({'title':os.path.basename(filename_csv) ,'parents': [{'id': '1yavW2rvTxSbByY0PLUf4a_S1rjl9uQbg'}]})
    gfile.SetContentFile(filename_csv)
    gfile.Upload(param={'supportsTeamDrives': True})

    #Export QC
    gfile = drive.CreateFile({'title':os.path.basename(pathScan) ,'parents': [{'id': '1WSD_FPJSpjLWUakwExm1hjwMHG8EmzjB'}],'mimeType': 'application/vnd.google-apps.folder'})
    gfile.Upload(param={'supportsTeamDrives': True})
    folder_id_parent = gfile['id']

    gfile = drive.CreateFile({'title':'Images_To_Check' ,'parents': [{'id': folder_id_parent}],'mimeType': 'application/vnd.google-apps.folder'})
    gfile.Upload(param={'supportsTeamDrives': True})
    folder_id_parent_imgcheck = gfile['id']

    gfile = drive.CreateFile({'title':'Maps' ,'parents': [{'id': folder_id_parent_imgcheck}],'mimeType': 'application/vnd.google-apps.folder'})
    gfile.Upload(param={'supportsTeamDrives': True})
    folder_id_parent_imgcheck_maps = gfile['id']

    gfile = drive.CreateFile({'title':'Aligments' ,'parents': [{'id': folder_id_parent_imgcheck}],'mimeType': 'application/vnd.google-apps.folder'})
    gfile.Upload(param={'supportsTeamDrives': True})
    folder_id_parent_imgcheck_aligments = gfile['id']

    extensions = ['.gif','.png','.xlsx']

    selected_files = []
    for filename in os.listdir(pathScan):
        if any(filename.endswith(ext) for ext in extensions):
            selected_files.append(os.path.join(pathScan, filename))

    for upload_file in selected_files:
        try:
            if 'alignment' in os.path.basename(upload_file):
                gfile = drive.CreateFile({'title':os.path.basename(upload_file) ,'parents': [{'id': folder_id_parent_imgcheck_aligments}]})
                gfile.SetContentFile(upload_file)
                gfile.Upload(param={'supportsTeamDrives': True})
            elif 'Masks' in os.path.basename(upload_file):
                gfile = drive.CreateFile({'title':os.path.basename(upload_file) ,'parents': [{'id': folder_id_parent}]})
                gfile.SetContentFile(upload_file)
                gfile.Upload(param={'supportsTeamDrives': True})
            elif '.xlsx' in os.path.basename(upload_file):
                gfile = drive.CreateFile({'title':os.path.basename(upload_file) ,'parents': [{'id': folder_id_parent}]})
                gfile.SetContentFile(upload_file)
                gfile.Upload(param={'supportsTeamDrives': True})
            elif '.csv' in os.path.basename(upload_file):
                gfile = drive.CreateFile({'title':os.path.basename(upload_file) ,'parents': [{'id': folder_id_parent}]})
                gfile.SetContentFile(upload_file)
                gfile.Upload(param={'supportsTeamDrives': True})
            else:
                gfile = drive.CreateFile({'title':os.path.basename(upload_file) ,'parents': [{'id': folder_id_parent_imgcheck_maps}]})
                gfile.SetContentFile(upload_file)
                gfile.Upload(param={'supportsTeamDrives': True})
        except:
            continue

    pathSegmentation = os.path.join(pathScan, 'segmentation_canvas')

    with zipfile.ZipFile(pathSegmentation + '.zip', 'w', zipfile.ZIP_DEFLATED) as zipf:
        zipdir(pathSegmentation, zipf)

    pathMask = os.path.join(pathScan, 'masks')
    
    with zipfile.ZipFile(pathMask + '.zip', 'w', zipfile.ZIP_DEFLATED) as zipf:
        zipdir(pathMask, zipf)

    upload_file_list = [pathSegmentation + '.zip', pathMask + '.zip',filename_log]

    for upload_file in upload_file_list:
        try:
            gfile = drive.CreateFile({'title':os.path.basename(upload_file) ,'parents': [{'id': folder_id_parent}]})
            gfile.SetContentFile(upload_file)
            gfile.Upload(param={'supportsTeamDrives': True}) # Upload the file.
        except:
            continue

    os.remove(pathMask + '.zip')
    os.remove(pathSegmentation + '.zip')

    #Export Complete Dataset
    with zipfile.ZipFile(pathScan + '.zip', 'w', zipfile.ZIP_DEFLATED) as zipf:
        zipdir(pathScan, zipf)

    gfile = drive.CreateFile({'title':os.path.basename(pathScan + '.zip') ,'parents': [{'id': '1TuoclmCv2TsEgSGIIBYxtFuC7T32fmnW'}]})
    gfile.SetContentFile(pathScan + '.zip')
    gfile.Upload(param={'supportsTeamDrives': True})

    os.remove(pathScan + '.zip')

def main(pathScan,filename_log,filename_csv):

    try:
        file = open(filename_log, 'a')
        file.write("\n"+str(datetime.datetime.now())[0:19] + ": Uploading to Google Drive has started")
        file.close()

        GoogleDrive_Upload(pathScan,filename_log,filename_csv)

    except Exception as e: 
        file = open(filename_log, 'a')
        file.write("\n"+str(datetime.datetime.now())[0:19] + ": Uploading to Google Drive was NOT completed ; error: "+str(e)) 
        file.close()


