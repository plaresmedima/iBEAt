from pydrive.auth import GoogleAuth
from pydrive.drive import GoogleDrive
import os


def upload_folder_to_drive(folder_path, parent_folder_id):
    # Authenticate and create GoogleDrive instance
    gauth = GoogleAuth()
    gauth.LocalWebserverAuth()  # Authenticate using local webserver
    drive = GoogleDrive(gauth)

    # Function to recursively upload folders and files
    def upload_folder_recursively(folder_path, parent_folder_id):
        # Create folder in Google Drive
        folder_name = os.path.basename(folder_path)
        current_folder = drive.CreateFile({
            'title': folder_name,
            'parents': [{'kind': 'drive#fileLink', 'id': parent_folder_id}],
            'mimeType': 'application/vnd.google-apps.folder'
        })
        current_folder.Upload()

        # Upload files in current directory
        for filename in os.listdir(folder_path):
            file_path = os.path.join(folder_path, filename)
            if os.path.isfile(file_path):
                drive_file = drive.CreateFile({
                    'title': filename,
                    'parents': [{'kind': 'drive#fileLink', 'id': current_folder['id']}]
                })
                drive_file.SetContentFile(file_path)
                drive_file.Upload()

        # Recursively upload subfolders
        for subdir in os.listdir(folder_path):
            subdir_path = os.path.join(folder_path, subdir)
            if os.path.isdir(subdir_path):
                upload_folder_recursively(subdir_path, current_folder['id'])

    # Call the recursive function to upload the folder and its contents
    upload_folder_recursively(folder_path, parent_folder_id)

# folder_path = 'C:\\Users\\md1jdsp\\Desktop\\Masks'
# parent_folder_id = '19t-am8P3vyWkqqw4JINOSPiFXDxstbo6'
# upload_folder_to_drive(folder_path, parent_folder_id)
