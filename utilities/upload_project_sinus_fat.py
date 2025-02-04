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

    pathSegmentation = os.path.join(pathScan + '_output')

    with zipfile.ZipFile(pathSegmentation + '.zip', 'w', zipfile.ZIP_DEFLATED) as zipf:
        zipdir(pathSegmentation, zipf)

    upload_file_list = [pathSegmentation + '.zip',filename_log]

    for upload_file in upload_file_list:
        try:
            gfile = drive.CreateFile({'title':os.path.basename(upload_file) ,'parents': [{'id': '1boIUH8ndGx4ecvNfN-gSdYwzzuj6jYkI'}]})
            gfile.SetContentFile(upload_file)
            gfile.Upload(param={'supportsTeamDrives': True}) # Upload the file.
        except:
            continue

    os.remove(pathSegmentation + '.zip')

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


