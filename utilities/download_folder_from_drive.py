from pydrive.auth import GoogleAuth
from pydrive.drive import GoogleDrive
import os

def download_folder_from_drive(folder_id, download_path):
    # Authenticate and create GoogleDrive instance
    gauth = GoogleAuth()
    gauth.LocalWebserverAuth()  # Creates local webserver and auto handles authentication
    drive = GoogleDrive(gauth)

    # Function to recursively download folders and files
    def download_folder_recursively(folder_id, download_path):
        # Get the folder metadata
        folder = drive.CreateFile({'id': folder_id})
        folder_title = folder['title']
        folder_path = os.path.join(download_path, folder_title)
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)

        # Query for files in the folder
        file_list = drive.ListFile({'q': f"'{folder_id}' in parents and trashed=false"}).GetList()

        # Download files in the folder
        for file in file_list:
            if file['mimeType'] == 'application/vnd.google-apps.folder':
                # If it's a folder, recursively download its contents
                download_folder_recursively(file['id'], folder_path)
            else:
                # If it's a file, download it
                file.GetContentFile(os.path.join(folder_path, file['title']))

    # Call the recursive function to download the folder and its contents
    download_folder_recursively(folder_id, download_path)

# Example usage:
# Specify the ID of the folder on Google Drive and the local download path
# Replace 'folder_id' and 'download_path' with appropriate values
# folder_id = '1L5KQkPWqULjop1RkcUbXxpqnTjDGBX1m'
# download_path = 'C:\\Users\\md1jdsp\\Downloads'
# download_folder_from_drive(folder_id, download_path)