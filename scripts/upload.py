import os
import os.path
import time
import datetime
import zipfile

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

    upload_file_list = [filename_log,pathScan + '.zip', filename_csv]
    #upload_file_list = [filename_log,pathScan + '.zip',pathScan + '_segmentation_results' + '.zip']
    #upload_file_list = [filename_log,pathScan + '_segmentation_results' + '.zip', filename_csv]
    #upload_file_list = [filename_log]

    for upload_file in upload_file_list:
        #gfile = drive.CreateFile({'title':os.path.basename(upload_file) ,'parents': [{'id': '1NvbNw00NaHpritRiYPKGC1l-4mOonIs4'}]})

        try:
            gfile = drive.CreateFile({'title':os.path.basename(upload_file) ,'parents': [{'id': '0AF9WRJGysgcAUk9PVA'}]})
            # Read file and set it as the content of this instance.
            gfile.SetContentFile(upload_file)
            gfile.Upload(param={'supportsTeamDrives': True}) # Upload the file.
        except:
            continue


def main(pathScan,filename_log):

    pathSegmentation = pathScan + '_segmentation_results'

    try:
        start_time = time.time()
        file = open(filename_log, 'a')
        file.write("\n"+str(datetime.datetime.now())[0:19] + ": Compressing into a zip file has started")
        file.close()

        with zipfile.ZipFile(pathScan + '.zip', 'w', zipfile.ZIP_DEFLATED) as zipf:
            zipdir(pathScan, zipf)

        # with zipfile.ZipFile(pathSegmentation + '.zip', 'w', zipfile.ZIP_DEFLATED) as zipf:
        #     zipdir(pathSegmentation, zipf)

        file = open(filename_log, 'a')
        file.write("\n"+str(datetime.datetime.now())[0:19] + ": Compressing into a zip file was completed --- %s seconds ---" % (int(time.time() - start_time))) 
        file.close()   

    except Exception as e: 
        file = open(filename_log, 'a')
        file.write("\n"+str(datetime.datetime.now())[0:19] + ": Compressing into a zip file was NOT completed ; error: "+str(e)) 
        file.close()

    try:
        start_time = time.time()
        file = open(filename_log, 'a')
        file.write("\n"+str(datetime.datetime.now())[0:19] + ": Uploading to Google Drive has started")
        file.close()

        GoogleDrive_Upload(pathScan,filename_log)

    except Exception as e: 
        file = open(filename_log, 'a')
        file.write("\n"+str(datetime.datetime.now())[0:19] + ": Uploading to Google Drive was NOT completed ; error: "+str(e)) 
        file.close()


