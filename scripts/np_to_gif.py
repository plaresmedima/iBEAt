import numpy as np
import imageio

def main(folder,array_3d,path,duration):

    list_of_series = folder.series()

    for i,series in enumerate(list_of_series):
        print(series['SeriesDescription'])
        if series["SequenceName"] is not None:

            if series['SeriesDescription'] == '*_mdr_moco_*':
                # Create frames from 3D array slices



                frames = []
                for slice in range(array_3d.shape[2]):
                    frame = (array_3d[:,:,slice]/np.max(array_3d[:,:,slice]) * 255).astype(np.uint8)  # Ensure data type is uint8
                    frames.append(frame) 

                # Save the frames as a GIF
                imageio.mimsave(path, frames, duration=duration)


