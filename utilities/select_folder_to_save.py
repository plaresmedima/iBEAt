import tkinter as tk
from tkinter import filedialog

def internal():
    # Create a Tkinter root window (it will not be shown)
    root = tk.Tk()
    root.withdraw()  # Hide the main window

    # Ask the user to select a folder
    selected_folder = filedialog.askdirectory(title="Select a folder to SAVE iBEAt results or to LOAD a previously downloaded dataset")

    # Print the selected folder
    print(f"Selected Folder: {selected_folder}")
    
    return selected_folder

def external():
    # Create a Tkinter root window (it will not be shown)
    root = tk.Tk()
    root.withdraw()  # Hide the main window

    # Ask the user to select a folder
    selected_folder = filedialog.askdirectory(title="Select a folder to LOAD a previously downloaded dataset from XNAT")

    # Print the selected folder
    print(f"Selected Folder: {selected_folder}")
    
    return selected_folder