import pipelines.segment as seg
import pipelines.volumetrics as volumetrics
import os
import os.path
import time
import datetime
import zipfile
import shutil

from pydrive.auth import GoogleAuth
from pydrive.drive import GoogleDrive


def main(master_table,folder):

    mask = 0

    #mask = find_mask_in_drive(folder)

    if mask==0:
        weights = 'UNETR_kidneys_v1.pth'
        seg.segment_kidneys(folder, weights)

    try:
        print('Starting volumetrics')
        master_table = volumetrics.main(master_table, folder)
    except Exception as e: 
        folder.log("Volumetrics was NOT completed; error: "+str(e))

    try:
        print('Starting exporting masks')
        seg.export_masks(folder)
    except Exception as e: 
        folder.log("Masks were NOT exported; error: "+str(e))

    try:
        print('Starting computing canvas')
        seg.compute_whole_kidney_canvas(folder)
    except Exception as e: 
        folder.log("Kidney canvas was NOT computed; error: "+str(e))

    try:
        print('Starting exporting canvas')
        seg.export_whole_kidney_canvas(folder)
    except Exception as e: 
        folder.log("Canvas was NOT exported; error: "+str(e))

    return master_table


def find_mask_in_drive(database):

    gauth = GoogleAuth()
    drive = GoogleDrive (gauth)

    folder = drive.CreateFile({'id': '1NvbNw00NaHpritRiYPKGC1l-4mOonIs4'})
    folder.FetchContent()

    for item in folder:
        print(f"Title: {item['title']}, ID: {item['id']}, Type: {item['mimeType']}")