import tkinter as tk
from tkinter import filedialog

def main():
    # Create a Tkinter root window (it will not be shown)
    root = tk.Tk()
    root.withdraw()  # Hide the main window

    # Ask the user to select a folder
    selected_folder = filedialog.askdirectory(title="Select a folder to save iBEAt results")

    # Print the selected folder
    print(f"Selected Folder: {selected_folder}")
    
    return selected_folder