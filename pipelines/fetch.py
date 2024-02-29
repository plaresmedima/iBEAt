from pydrive.auth import GoogleAuth
from pydrive.drive import GoogleDrive

# Work in progress
def find_mask_in_drive(database):

    gauth = GoogleAuth()
    drive = GoogleDrive (gauth)

    folder = drive.CreateFile({'id': '1NvbNw00NaHpritRiYPKGC1l-4mOonIs4'})
    folder.FetchContent()

    for item in folder:
        print(f"Title: {item['title']}, ID: {item['id']}, Type: {item['mimeType']}")


# Dummy placeholder function
def kidney_masks(folder):
    # Check for masks on google drive and import into the database if exists.
    # Alternative: save them on XNAT on a separate project and download from there
    pass


