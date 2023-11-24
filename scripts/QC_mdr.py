import os
import numpy as np
import imageio

def main(database, duration = 100):
    list_of_series = database.series()
    for series in (list_of_series):
        if 'mdr_moco' in series['SeriesDescription']:

            array_temp, header_temp = series.array(['SliceLocation', 'AcquisitionTime'], pixels_first=True)#
            np.squeeze(array_temp)
            if len(array_temp) > 3: 
                array_temp = reshape_to_3d(array_temp)

            frames = []
            for slice in range(array_temp.shape[2]):
                frame = (array_temp[:,:,slice]/np.max(array_temp[:,:,slice]) * 255).astype(np.uint8)  # Ensure data type is uint8
                frame = np.transpose(frame)
                frames.append(frame) 

            # Save the frames as a GIF
            imageio.mimsave(os.path.join(database.path(),series['SeriesDescription'] + '.gif'), frames, duration=duration)
            
def reshape_to_3d(arr):
    # Calculate the size for the third dimension using np.prod
    new_third_dimension_size = np.prod(arr.shape[2:])

    # Reshape the array while preserving the first and second dimensions
    reshaped_array = np.reshape(arr, (arr.shape[0], arr.shape[1], new_third_dimension_size))
    
    return reshaped_array