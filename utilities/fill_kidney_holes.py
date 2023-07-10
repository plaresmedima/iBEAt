import numpy as np
import matplotlib.pyplot as plt

def main(input_array, array_mask):

    array_Data_masked = input_array*array_mask
    array_Data_masked[array_Data_masked<=1] = 0 #to guarantee 0 are zeros and no 10^-13

    #Detect the regions in the kidney to fill with values and replace by nans
    a = np.zeros(array_mask.shape) 
    a [array_mask==1] = 1

    b = np.zeros(array_Data_masked.shape) 
    b [array_Data_masked==0] = 1

    c= np.logical_and(a, b)

    array_Data_masked[c==1] = np.nan

    #find the valid indices to take information (ignoring nans and zeros)
    arr = array_Data_masked

    nan_indices = np.isnan(arr)
    non_nan_indices = ~nan_indices
    nonzero_indices = (arr != 0)

    print("Number of Nans in the kidney: " + str(np.count_nonzero(np.isnan(arr))))

    valid_indices = np.logical_and(non_nan_indices, nonzero_indices)

    # Loop through every pixel and if it finds a nan look for all 26 neighbours. At the end replace the nan pixel with the average of the valid neighbours  
    radius = 1

    for x in range (arr.shape[0]):
        for y in range(arr.shape[1]):
            for z in range(arr.shape[2]):
                if np.isnan(arr[x,y,z]):
                    possible_value_to_mean =[]
                    if (arr[x+radius,y,z] !=0) and (valid_indices[x+radius,y,z]==1): #Right 
                        possible_value_to_mean.append(arr[x+radius,y,z])
                    if (arr[x-radius,y,z] !=0) and (valid_indices[x-radius,y,z]==1): #Left
                        possible_value_to_mean.append(arr[x-radius,y,z])
                    if (arr[x,y+radius,z] !=0) and (valid_indices[x,y+radius,z]==1): #Top
                        possible_value_to_mean.append(arr[x,y+radius,z])
                    if (arr[x,y-radius,z] !=0) and (valid_indices[x,y-radius,z]==1): #Bottom
                        possible_value_to_mean.append(arr[x,y-radius,z])
                    if (arr[x,y,z+radius] !=0) and (valid_indices[x,y,z+radius]==1): #Front
                        possible_value_to_mean.append(arr[x,y,z+radius])
                    if (arr[x,y,z-radius] !=0) and (valid_indices[x,y,z-radius]==1): #Back
                        possible_value_to_mean.append(arr[x,y,z-radius])                   
                    if (arr[x+radius,y+radius,z] !=0) and (valid_indices[x+radius,y+radius,z]==1): #Top Diagonal Right
                        possible_value_to_mean.append(arr[x+radius,y+radius,z])
                    if (arr[x-radius,y-radius,z] !=0) and (valid_indices[x-radius,y-radius,z]==1): #Bottom Diagonal Left
                        possible_value_to_mean.append(arr[x-radius,y-radius,z])
                    if (arr[x+radius,y-radius,z] !=0) and (valid_indices[x+radius,y-radius,z]==1): #Bottom Diagonal Right
                        possible_value_to_mean.append(arr[x+radius,y-radius,z])
                    if (arr[x-radius,y+radius,z] !=0) and (valid_indices[x-radius,y+radius,z]==1): #Top Diagonal Left
                        possible_value_to_mean.append(arr[x-radius,y+radius,z])
                    #Z+radius
                    if (arr[x+radius,y+radius,z+radius] !=0) and (valid_indices[x+radius,y+radius,z+radius]==1): #Top Diagonal Right
                        possible_value_to_mean.append(arr[x+radius,y+radius,z+radius])
                    if (arr[x-radius,y-radius,z+radius] !=0) and (valid_indices[x-radius,y-radius,z+radius]==1): #Bottom Diagonal Left
                        possible_value_to_mean.append(arr[x-radius,y-radius,z+radius])
                    if (arr[x+radius,y-radius,z+radius] !=0) and (valid_indices[x+radius,y-radius,z+radius]==1): #Bottom Diagonal Right
                        possible_value_to_mean.append(arr[x+radius,y-radius,z+radius])
                    if (arr[x-radius,y+radius,z+radius] !=0) and (valid_indices[x-radius,y+radius,z+radius]==1): #Top Diagonal Left
                        possible_value_to_mean.append(arr[x-radius,y+radius,z+radius])
                    if (arr[x+radius,y,z+radius] !=0) and (valid_indices[x+radius,y,z+radius]==1):
                        possible_value_to_mean.append(arr[x+radius,y,z+radius])
                    if (arr[x-radius,y,z+radius] !=0) and (valid_indices[x-radius,y,z+radius]==1):
                        possible_value_to_mean.append(arr[x-radius,y,z+radius])
                    if (arr[x,y+radius,z+radius] !=0) and (valid_indices[x,y+radius,z+radius]==1):
                        possible_value_to_mean.append(arr[x,y+radius,z+radius])
                    if (arr[x-radius,y+radius,z+radius] !=0) and (valid_indices[x-radius,y+radius,z+radius]==1):
                        possible_value_to_mean.append(arr[x-radius,y+radius,z+radius])
                    #Z-radius
                    if (arr[x+radius,y+radius,z-radius] !=0) and (valid_indices[x+radius,y+radius,z-radius]==1): #Top Diagonal Right
                        possible_value_to_mean.append(arr[x+radius,y+radius,z-radius])
                    if (arr[x-radius,y-radius,z-radius] !=0) and (valid_indices[x-radius,y-radius,z-radius]==1): #Bottom Diagonal Left
                        possible_value_to_mean.append(arr[x-radius,y-radius,z-radius])
                    if (arr[x+radius,y-radius,z-radius] !=0) and (valid_indices[x+radius,y-radius,z-radius]==1): #Bottom Diagonal Right
                        possible_value_to_mean.append(arr[x+radius,y-radius,z-radius])
                    if (arr[x-radius,y+radius,z-radius] !=0) and (valid_indices[x-radius,y+radius,z-radius]==1): #Top Diagonal Left
                        possible_value_to_mean.append(arr[x-radius,y+radius,z-radius])
                    if (arr[x+radius,y,z-radius] !=0) and (valid_indices[x+radius,y,z-radius]==1):
                        possible_value_to_mean.append(arr[x+radius,y,z-radius])
                    if (arr[x-radius,y,z-radius] !=0) and (valid_indices[x-radius,y,z-radius]==1):
                        possible_value_to_mean.append(arr[x-radius,y,z-radius])
                    if (arr[x,y+radius,z-radius] !=0) and (valid_indices[x,y+radius,z-radius]==1):
                        possible_value_to_mean.append(arr[x,y+radius,z-radius])
                    if (arr[x-radius,y+radius,z-radius] !=0) and (valid_indices[x-radius,y+radius,z-radius]==1):
                        possible_value_to_mean.append(arr[x-radius,y+radius,z-radius])
        
                    arr[x,y,z] = np.average(possible_value_to_mean)


    # plt.imshow(arr[:,:,90])
    # plt.show()

    print("Number of Nans in the kidney: " + str(np.count_nonzero(np.isnan(arr))))

    nan_indices = np.isnan(arr)

    # Find indices of non-NaN and non-zero values
    non_nan_indices = ~nan_indices
    nonzero_indices = (arr != 0)

    # Combine non-NaN and non-zero indices
    valid_indices = np.logical_and(non_nan_indices, nonzero_indices)

    x_range = range (arr.shape[0])
    y_range = range (arr.shape[1])
    z_range = range (arr.shape[2])

    for x in reversed (x_range):
        for y in reversed (y_range):
            for z in reversed(z_range):
                if np.isnan(arr[x,y,z]):
                    possible_value_to_mean =[]
                    if (arr[x+radius,y,z] !=0) and (valid_indices[x+radius,y,z]==1): #Right 
                        possible_value_to_mean.append(arr[x+radius,y,z])
                    if (arr[x-radius,y,z] !=0) and (valid_indices[x-radius,y,z]==1): #Left
                        possible_value_to_mean.append(arr[x-radius,y,z])
                    if (arr[x,y+radius,z] !=0) and (valid_indices[x,y+radius,z]==1): #Top
                        possible_value_to_mean.append(arr[x,y+radius,z])
                    if (arr[x,y-radius,z] !=0) and (valid_indices[x,y-radius,z]==1): #Bottom
                        possible_value_to_mean.append(arr[x,y-radius,z])
                    if (arr[x,y,z+radius] !=0) and (valid_indices[x,y,z+radius]==1): #Front
                        possible_value_to_mean.append(arr[x,y,z+radius])
                    if (arr[x,y,z-radius] !=0) and (valid_indices[x,y,z-radius]==1): #Back
                        possible_value_to_mean.append(arr[x,y,z-radius])                   
                    if (arr[x+radius,y+radius,z] !=0) and (valid_indices[x+radius,y+radius,z]==1): #Top Diagonal Right
                        possible_value_to_mean.append(arr[x+radius,y+radius,z])
                    if (arr[x-radius,y-radius,z] !=0) and (valid_indices[x-radius,y-radius,z]==1): #Bottom Diagonal Left
                        possible_value_to_mean.append(arr[x-radius,y-radius,z])
                    if (arr[x+radius,y-radius,z] !=0) and (valid_indices[x+radius,y-radius,z]==1): #Bottom Diagonal Right
                        possible_value_to_mean.append(arr[x+radius,y-radius,z])
                    if (arr[x-radius,y+radius,z] !=0) and (valid_indices[x-radius,y+radius,z]==1): #Top Diagonal Left
                        possible_value_to_mean.append(arr[x-radius,y+radius,z])
                    #Z+radius
                    if (arr[x+radius,y+radius,z+radius] !=0) and (valid_indices[x+radius,y+radius,z+radius]==1): #Top Diagonal Right
                        possible_value_to_mean.append(arr[x+radius,y+radius,z+radius])
                    if (arr[x-radius,y-radius,z+radius] !=0) and (valid_indices[x-radius,y-radius,z+radius]==1): #Bottom Diagonal Left
                        possible_value_to_mean.append(arr[x-radius,y-radius,z+radius])
                    if (arr[x+radius,y-radius,z+radius] !=0) and (valid_indices[x+radius,y-radius,z+radius]==1): #Bottom Diagonal Right
                        possible_value_to_mean.append(arr[x+radius,y-radius,z+radius])
                    if (arr[x-radius,y+radius,z+radius] !=0) and (valid_indices[x-radius,y+radius,z+radius]==1): #Top Diagonal Left
                        possible_value_to_mean.append(arr[x-radius,y+radius,z+radius])
                    if (arr[x+radius,y,z+radius] !=0) and (valid_indices[x+radius,y,z+radius]==1):
                        possible_value_to_mean.append(arr[x+radius,y,z+radius])
                    if (arr[x-radius,y,z+radius] !=0) and (valid_indices[x-radius,y,z+radius]==1):
                        possible_value_to_mean.append(arr[x-radius,y,z+radius])
                    if (arr[x,y+radius,z+radius] !=0) and (valid_indices[x,y+radius,z+radius]==1):
                        possible_value_to_mean.append(arr[x,y+radius,z+radius])
                    if (arr[x-radius,y+radius,z+radius] !=0) and (valid_indices[x-radius,y+radius,z+radius]==1):
                        possible_value_to_mean.append(arr[x-radius,y+radius,z+radius])
                    #Z-radius
                    if (arr[x+radius,y+radius,z-radius] !=0) and (valid_indices[x+radius,y+radius,z-radius]==1): #Top Diagonal Right
                        possible_value_to_mean.append(arr[x+radius,y+radius,z-radius])
                    if (arr[x-radius,y-radius,z-radius] !=0) and (valid_indices[x-radius,y-radius,z-radius]==1): #Bottom Diagonal Left
                        possible_value_to_mean.append(arr[x-radius,y-radius,z-radius])
                    if (arr[x+radius,y-radius,z-radius] !=0) and (valid_indices[x+radius,y-radius,z-radius]==1): #Bottom Diagonal Right
                        possible_value_to_mean.append(arr[x+radius,y-radius,z-radius])
                    if (arr[x-radius,y+radius,z-radius] !=0) and (valid_indices[x-radius,y+radius,z-radius]==1): #Top Diagonal Left
                        possible_value_to_mean.append(arr[x-radius,y+radius,z-radius])
                    if (arr[x+radius,y,z-radius] !=0) and (valid_indices[x+radius,y,z-radius]==1):
                        possible_value_to_mean.append(arr[x+radius,y,z-radius])
                    if (arr[x-radius,y,z-radius] !=0) and (valid_indices[x-radius,y,z-radius]==1):
                        possible_value_to_mean.append(arr[x-radius,y,z-radius])
                    if (arr[x,y+radius,z-radius] !=0) and (valid_indices[x,y+radius,z-radius]==1):
                        possible_value_to_mean.append(arr[x,y+radius,z-radius])
                    if (arr[x-radius,y+radius,z-radius] !=0) and (valid_indices[x-radius,y+radius,z-radius]==1):
                        possible_value_to_mean.append(arr[x-radius,y+radius,z-radius])
        
                    arr[x,y,z] = np.average(possible_value_to_mean)

    print("Number of Nans in the kidney: " + str(np.count_nonzero(np.isnan(arr))))

    # plt.imshow(arr[:,:,90])
    # plt.show()

    nan_indices = np.isnan(arr)

    # Find indices of non-NaN and non-zero values
    non_nan_indices = ~nan_indices
    nonzero_indices = (arr != 0)

    # Combine non-NaN and non-zero indices
    valid_indices = np.logical_and(non_nan_indices, nonzero_indices)

    x_range = range (arr.shape[0])
    y_range = range (arr.shape[1])
    z_range = range (arr.shape[2])

    for z in reversed (z_range):
        for y in reversed (y_range):
            for x in reversed(x_range):
                if np.isnan(arr[x,y,z]):
                    possible_value_to_mean =[]
                    if (arr[x+radius,y,z] !=0) and (valid_indices[x+radius,y,z]==1): #Right 
                        possible_value_to_mean.append(arr[x+radius,y,z])
                    if (arr[x-radius,y,z] !=0) and (valid_indices[x-radius,y,z]==1): #Left
                        possible_value_to_mean.append(arr[x-radius,y,z])
                    if (arr[x,y+radius,z] !=0) and (valid_indices[x,y+radius,z]==1): #Top
                        possible_value_to_mean.append(arr[x,y+radius,z])
                    if (arr[x,y-radius,z] !=0) and (valid_indices[x,y-radius,z]==1): #Bottom
                        possible_value_to_mean.append(arr[x,y-radius,z])
                    if (arr[x,y,z+radius] !=0) and (valid_indices[x,y,z+radius]==1): #Front
                        possible_value_to_mean.append(arr[x,y,z+radius])
                    if (arr[x,y,z-radius] !=0) and (valid_indices[x,y,z-radius]==1): #Back
                        possible_value_to_mean.append(arr[x,y,z-radius])                   
                    if (arr[x+radius,y+radius,z] !=0) and (valid_indices[x+radius,y+radius,z]==1): #Top Diagonal Right
                        possible_value_to_mean.append(arr[x+radius,y+radius,z])
                    if (arr[x-radius,y-radius,z] !=0) and (valid_indices[x-radius,y-radius,z]==1): #Bottom Diagonal Left
                        possible_value_to_mean.append(arr[x-radius,y-radius,z])
                    if (arr[x+radius,y-radius,z] !=0) and (valid_indices[x+radius,y-radius,z]==1): #Bottom Diagonal Right
                        possible_value_to_mean.append(arr[x+radius,y-radius,z])
                    if (arr[x-radius,y+radius,z] !=0) and (valid_indices[x-radius,y+radius,z]==1): #Top Diagonal Left
                        possible_value_to_mean.append(arr[x-radius,y+radius,z])
                    #Z+radius
                    if (arr[x+radius,y+radius,z+radius] !=0) and (valid_indices[x+radius,y+radius,z+radius]==1): #Top Diagonal Right
                        possible_value_to_mean.append(arr[x+radius,y+radius,z+radius])
                    if (arr[x-radius,y-radius,z+radius] !=0) and (valid_indices[x-radius,y-radius,z+radius]==1): #Bottom Diagonal Left
                        possible_value_to_mean.append(arr[x-radius,y-radius,z+radius])
                    if (arr[x+radius,y-radius,z+radius] !=0) and (valid_indices[x+radius,y-radius,z+radius]==1): #Bottom Diagonal Right
                        possible_value_to_mean.append(arr[x+radius,y-radius,z+radius])
                    if (arr[x-radius,y+radius,z+radius] !=0) and (valid_indices[x-radius,y+radius,z+radius]==1): #Top Diagonal Left
                        possible_value_to_mean.append(arr[x-radius,y+radius,z+radius])
                    if (arr[x+radius,y,z+radius] !=0) and (valid_indices[x+radius,y,z+radius]==1):
                        possible_value_to_mean.append(arr[x+radius,y,z+radius])
                    if (arr[x-radius,y,z+radius] !=0) and (valid_indices[x-radius,y,z+radius]==1):
                        possible_value_to_mean.append(arr[x-radius,y,z+radius])
                    if (arr[x,y+radius,z+radius] !=0) and (valid_indices[x,y+radius,z+radius]==1):
                        possible_value_to_mean.append(arr[x,y+radius,z+radius])
                    if (arr[x-radius,y+radius,z+radius] !=0) and (valid_indices[x-radius,y+radius,z+radius]==1):
                        possible_value_to_mean.append(arr[x-radius,y+radius,z+radius])
                    #Z-radius
                    if (arr[x+radius,y+radius,z-radius] !=0) and (valid_indices[x+radius,y+radius,z-radius]==1): #Top Diagonal Right
                        possible_value_to_mean.append(arr[x+radius,y+radius,z-radius])
                    if (arr[x-radius,y-radius,z-radius] !=0) and (valid_indices[x-radius,y-radius,z-radius]==1): #Bottom Diagonal Left
                        possible_value_to_mean.append(arr[x-radius,y-radius,z-radius])
                    if (arr[x+radius,y-radius,z-radius] !=0) and (valid_indices[x+radius,y-radius,z-radius]==1): #Bottom Diagonal Right
                        possible_value_to_mean.append(arr[x+radius,y-radius,z-radius])
                    if (arr[x-radius,y+radius,z-radius] !=0) and (valid_indices[x-radius,y+radius,z-radius]==1): #Top Diagonal Left
                        possible_value_to_mean.append(arr[x-radius,y+radius,z-radius])
                    if (arr[x+radius,y,z-radius] !=0) and (valid_indices[x+radius,y,z-radius]==1):
                        possible_value_to_mean.append(arr[x+radius,y,z-radius])
                    if (arr[x-radius,y,z-radius] !=0) and (valid_indices[x-radius,y,z-radius]==1):
                        possible_value_to_mean.append(arr[x-radius,y,z-radius])
                    if (arr[x,y+radius,z-radius] !=0) and (valid_indices[x,y+radius,z-radius]==1):
                        possible_value_to_mean.append(arr[x,y+radius,z-radius])
                    if (arr[x-radius,y+radius,z-radius] !=0) and (valid_indices[x-radius,y+radius,z-radius]==1):
                        possible_value_to_mean.append(arr[x-radius,y+radius,z-radius])
        
                    arr[x,y,z] = np.average(possible_value_to_mean)

    print("Number of Nans in the kidney: " + str(np.count_nonzero(np.isnan(arr))))

    # plt.imshow(arr[:,:,90])
    # plt.show()

    nan_indices = np.isnan(arr)

    # Find indices of non-NaN and non-zero values
    non_nan_indices = ~nan_indices
    nonzero_indices = (arr != 0)

    # Combine non-NaN and non-zero indices
    valid_indices = np.logical_and(non_nan_indices, nonzero_indices)

    x_range = range (arr.shape[0])
    y_range = range (arr.shape[1])
    z_range = range (arr.shape[2])

    for z in reversed (z_range):
        for x in reversed (x_range):
            for y in reversed(y_range):
                if np.isnan(arr[x,y,z]):
                    possible_value_to_mean =[]
                    if (arr[x+radius,y,z] !=0) and (valid_indices[x+radius,y,z]==1): #Right 
                        possible_value_to_mean.append(arr[x+radius,y,z])
                    if (arr[x-radius,y,z] !=0) and (valid_indices[x-radius,y,z]==1): #Left
                        possible_value_to_mean.append(arr[x-radius,y,z])
                    if (arr[x,y+radius,z] !=0) and (valid_indices[x,y+radius,z]==1): #Top
                        possible_value_to_mean.append(arr[x,y+radius,z])
                    if (arr[x,y-radius,z] !=0) and (valid_indices[x,y-radius,z]==1): #Bottom
                        possible_value_to_mean.append(arr[x,y-radius,z])
                    if (arr[x,y,z+radius] !=0) and (valid_indices[x,y,z+radius]==1): #Front
                        possible_value_to_mean.append(arr[x,y,z+radius])
                    if (arr[x,y,z-radius] !=0) and (valid_indices[x,y,z-radius]==1): #Back
                        possible_value_to_mean.append(arr[x,y,z-radius])                   
                    if (arr[x+radius,y+radius,z] !=0) and (valid_indices[x+radius,y+radius,z]==1): #Top Diagonal Right
                        possible_value_to_mean.append(arr[x+radius,y+radius,z])
                    if (arr[x-radius,y-radius,z] !=0) and (valid_indices[x-radius,y-radius,z]==1): #Bottom Diagonal Left
                        possible_value_to_mean.append(arr[x-radius,y-radius,z])
                    if (arr[x+radius,y-radius,z] !=0) and (valid_indices[x+radius,y-radius,z]==1): #Bottom Diagonal Right
                        possible_value_to_mean.append(arr[x+radius,y-radius,z])
                    if (arr[x-radius,y+radius,z] !=0) and (valid_indices[x-radius,y+radius,z]==1): #Top Diagonal Left
                        possible_value_to_mean.append(arr[x-radius,y+radius,z])
                    #Z+radius
                    if (arr[x+radius,y+radius,z+radius] !=0) and (valid_indices[x+radius,y+radius,z+radius]==1): #Top Diagonal Right
                        possible_value_to_mean.append(arr[x+radius,y+radius,z+radius])
                    if (arr[x-radius,y-radius,z+radius] !=0) and (valid_indices[x-radius,y-radius,z+radius]==1): #Bottom Diagonal Left
                        possible_value_to_mean.append(arr[x-radius,y-radius,z+radius])
                    if (arr[x+radius,y-radius,z+radius] !=0) and (valid_indices[x+radius,y-radius,z+radius]==1): #Bottom Diagonal Right
                        possible_value_to_mean.append(arr[x+radius,y-radius,z+radius])
                    if (arr[x-radius,y+radius,z+radius] !=0) and (valid_indices[x-radius,y+radius,z+radius]==1): #Top Diagonal Left
                        possible_value_to_mean.append(arr[x-radius,y+radius,z+radius])
                    if (arr[x+radius,y,z+radius] !=0) and (valid_indices[x+radius,y,z+radius]==1):
                        possible_value_to_mean.append(arr[x+radius,y,z+radius])
                    if (arr[x-radius,y,z+radius] !=0) and (valid_indices[x-radius,y,z+radius]==1):
                        possible_value_to_mean.append(arr[x-radius,y,z+radius])
                    if (arr[x,y+radius,z+radius] !=0) and (valid_indices[x,y+radius,z+radius]==1):
                        possible_value_to_mean.append(arr[x,y+radius,z+radius])
                    if (arr[x-radius,y+radius,z+radius] !=0) and (valid_indices[x-radius,y+radius,z+radius]==1):
                        possible_value_to_mean.append(arr[x-radius,y+radius,z+radius])
                    #Z-radius
                    if (arr[x+radius,y+radius,z-radius] !=0) and (valid_indices[x+radius,y+radius,z-radius]==1): #Top Diagonal Right
                        possible_value_to_mean.append(arr[x+radius,y+radius,z-radius])
                    if (arr[x-radius,y-radius,z-radius] !=0) and (valid_indices[x-radius,y-radius,z-radius]==1): #Bottom Diagonal Left
                        possible_value_to_mean.append(arr[x-radius,y-radius,z-radius])
                    if (arr[x+radius,y-radius,z-radius] !=0) and (valid_indices[x+radius,y-radius,z-radius]==1): #Bottom Diagonal Right
                        possible_value_to_mean.append(arr[x+radius,y-radius,z-radius])
                    if (arr[x-radius,y+radius,z-radius] !=0) and (valid_indices[x-radius,y+radius,z-radius]==1): #Top Diagonal Left
                        possible_value_to_mean.append(arr[x-radius,y+radius,z-radius])
                    if (arr[x+radius,y,z-radius] !=0) and (valid_indices[x+radius,y,z-radius]==1):
                        possible_value_to_mean.append(arr[x+radius,y,z-radius])
                    if (arr[x-radius,y,z-radius] !=0) and (valid_indices[x-radius,y,z-radius]==1):
                        possible_value_to_mean.append(arr[x-radius,y,z-radius])
                    if (arr[x,y+radius,z-radius] !=0) and (valid_indices[x,y+radius,z-radius]==1):
                        possible_value_to_mean.append(arr[x,y+radius,z-radius])
                    if (arr[x-radius,y+radius,z-radius] !=0) and (valid_indices[x-radius,y+radius,z-radius]==1):
                        possible_value_to_mean.append(arr[x-radius,y+radius,z-radius])
        
                    arr[x,y,z] = np.average(possible_value_to_mean)

    print("Number of Nans in the kidney: " + str(np.count_nonzero(np.isnan(arr))))

    # plt.imshow(arr[:,:,90])
    # plt.show()

    nan_indices = np.isnan(arr)

    # Find indices of non-NaN and non-zero values
    non_nan_indices = ~nan_indices
    nonzero_indices = (arr != 0)

    # Combine non-NaN and non-zero indices
    valid_indices = np.logical_and(non_nan_indices, nonzero_indices)

    x_range = range (arr.shape[0])
    y_range = range (arr.shape[1])
    z_range = range (arr.shape[2])

    for y in reversed (y_range):
        for z in reversed (z_range):
            for x in reversed(x_range):
                if np.isnan(arr[x,y,z]):
                    possible_value_to_mean =[]
                    if (arr[x+radius,y,z] !=0) and (valid_indices[x+radius,y,z]==1): #Right 
                        possible_value_to_mean.append(arr[x+radius,y,z])
                    if (arr[x-radius,y,z] !=0) and (valid_indices[x-radius,y,z]==1): #Left
                        possible_value_to_mean.append(arr[x-radius,y,z])
                    if (arr[x,y+radius,z] !=0) and (valid_indices[x,y+radius,z]==1): #Top
                        possible_value_to_mean.append(arr[x,y+radius,z])
                    if (arr[x,y-radius,z] !=0) and (valid_indices[x,y-radius,z]==1): #Bottom
                        possible_value_to_mean.append(arr[x,y-radius,z])
                    if (arr[x,y,z+radius] !=0) and (valid_indices[x,y,z+radius]==1): #Front
                        possible_value_to_mean.append(arr[x,y,z+radius])
                    if (arr[x,y,z-radius] !=0) and (valid_indices[x,y,z-radius]==1): #Back
                        possible_value_to_mean.append(arr[x,y,z-radius])                   
                    if (arr[x+radius,y+radius,z] !=0) and (valid_indices[x+radius,y+radius,z]==1): #Top Diagonal Right
                        possible_value_to_mean.append(arr[x+radius,y+radius,z])
                    if (arr[x-radius,y-radius,z] !=0) and (valid_indices[x-radius,y-radius,z]==1): #Bottom Diagonal Left
                        possible_value_to_mean.append(arr[x-radius,y-radius,z])
                    if (arr[x+radius,y-radius,z] !=0) and (valid_indices[x+radius,y-radius,z]==1): #Bottom Diagonal Right
                        possible_value_to_mean.append(arr[x+radius,y-radius,z])
                    if (arr[x-radius,y+radius,z] !=0) and (valid_indices[x-radius,y+radius,z]==1): #Top Diagonal Left
                        possible_value_to_mean.append(arr[x-radius,y+radius,z])
                    #Z+radius
                    if (arr[x+radius,y+radius,z+radius] !=0) and (valid_indices[x+radius,y+radius,z+radius]==1): #Top Diagonal Right
                        possible_value_to_mean.append(arr[x+radius,y+radius,z+radius])
                    if (arr[x-radius,y-radius,z+radius] !=0) and (valid_indices[x-radius,y-radius,z+radius]==1): #Bottom Diagonal Left
                        possible_value_to_mean.append(arr[x-radius,y-radius,z+radius])
                    if (arr[x+radius,y-radius,z+radius] !=0) and (valid_indices[x+radius,y-radius,z+radius]==1): #Bottom Diagonal Right
                        possible_value_to_mean.append(arr[x+radius,y-radius,z+radius])
                    if (arr[x-radius,y+radius,z+radius] !=0) and (valid_indices[x-radius,y+radius,z+radius]==1): #Top Diagonal Left
                        possible_value_to_mean.append(arr[x-radius,y+radius,z+radius])
                    if (arr[x+radius,y,z+radius] !=0) and (valid_indices[x+radius,y,z+radius]==1):
                        possible_value_to_mean.append(arr[x+radius,y,z+radius])
                    if (arr[x-radius,y,z+radius] !=0) and (valid_indices[x-radius,y,z+radius]==1):
                        possible_value_to_mean.append(arr[x-radius,y,z+radius])
                    if (arr[x,y+radius,z+radius] !=0) and (valid_indices[x,y+radius,z+radius]==1):
                        possible_value_to_mean.append(arr[x,y+radius,z+radius])
                    if (arr[x-radius,y+radius,z+radius] !=0) and (valid_indices[x-radius,y+radius,z+radius]==1):
                        possible_value_to_mean.append(arr[x-radius,y+radius,z+radius])
                    #Z-radius
                    if (arr[x+radius,y+radius,z-radius] !=0) and (valid_indices[x+radius,y+radius,z-radius]==1): #Top Diagonal Right
                        possible_value_to_mean.append(arr[x+radius,y+radius,z-radius])
                    if (arr[x-radius,y-radius,z-radius] !=0) and (valid_indices[x-radius,y-radius,z-radius]==1): #Bottom Diagonal Left
                        possible_value_to_mean.append(arr[x-radius,y-radius,z-radius])
                    if (arr[x+radius,y-radius,z-radius] !=0) and (valid_indices[x+radius,y-radius,z-radius]==1): #Bottom Diagonal Right
                        possible_value_to_mean.append(arr[x+radius,y-radius,z-radius])
                    if (arr[x-radius,y+radius,z-radius] !=0) and (valid_indices[x-radius,y+radius,z-radius]==1): #Top Diagonal Left
                        possible_value_to_mean.append(arr[x-radius,y+radius,z-radius])
                    if (arr[x+radius,y,z-radius] !=0) and (valid_indices[x+radius,y,z-radius]==1):
                        possible_value_to_mean.append(arr[x+radius,y,z-radius])
                    if (arr[x-radius,y,z-radius] !=0) and (valid_indices[x-radius,y,z-radius]==1):
                        possible_value_to_mean.append(arr[x-radius,y,z-radius])
                    if (arr[x,y+radius,z-radius] !=0) and (valid_indices[x,y+radius,z-radius]==1):
                        possible_value_to_mean.append(arr[x,y+radius,z-radius])
                    if (arr[x-radius,y+radius,z-radius] !=0) and (valid_indices[x-radius,y+radius,z-radius]==1):
                        possible_value_to_mean.append(arr[x-radius,y+radius,z-radius])
        
                    arr[x,y,z] = np.average(possible_value_to_mean)

    print("Number of Nans in the kidney: " + str(np.count_nonzero(np.isnan(arr))))

    # plt.imshow(arr[:,:,90])
    # plt.show()

    nan_indices = np.isnan(arr)

    # Find indices of non-NaN and non-zero values
    non_nan_indices = ~nan_indices
    nonzero_indices = (arr != 0)

    # Combine non-NaN and non-zero indices
    valid_indices = np.logical_and(non_nan_indices, nonzero_indices)

    x_range = range (arr.shape[0])
    y_range = range (arr.shape[1])
    z_range = range (arr.shape[2])

    for x in reversed (x_range):
        for z in reversed (z_range):
            for y in reversed(y_range):
                if np.isnan(arr[x,y,z]):
                    possible_value_to_mean =[]
                    if (arr[x+radius,y,z] !=0) and (valid_indices[x+radius,y,z]==1): #Right 
                        possible_value_to_mean.append(arr[x+radius,y,z])
                    if (arr[x-radius,y,z] !=0) and (valid_indices[x-radius,y,z]==1): #Left
                        possible_value_to_mean.append(arr[x-radius,y,z])
                    if (arr[x,y+radius,z] !=0) and (valid_indices[x,y+radius,z]==1): #Top
                        possible_value_to_mean.append(arr[x,y+radius,z])
                    if (arr[x,y-radius,z] !=0) and (valid_indices[x,y-radius,z]==1): #Bottom
                        possible_value_to_mean.append(arr[x,y-radius,z])
                    if (arr[x,y,z+radius] !=0) and (valid_indices[x,y,z+radius]==1): #Front
                        possible_value_to_mean.append(arr[x,y,z+radius])
                    if (arr[x,y,z-radius] !=0) and (valid_indices[x,y,z-radius]==1): #Back
                        possible_value_to_mean.append(arr[x,y,z-radius])                   
                    if (arr[x+radius,y+radius,z] !=0) and (valid_indices[x+radius,y+radius,z]==1): #Top Diagonal Right
                        possible_value_to_mean.append(arr[x+radius,y+radius,z])
                    if (arr[x-radius,y-radius,z] !=0) and (valid_indices[x-radius,y-radius,z]==1): #Bottom Diagonal Left
                        possible_value_to_mean.append(arr[x-radius,y-radius,z])
                    if (arr[x+radius,y-radius,z] !=0) and (valid_indices[x+radius,y-radius,z]==1): #Bottom Diagonal Right
                        possible_value_to_mean.append(arr[x+radius,y-radius,z])
                    if (arr[x-radius,y+radius,z] !=0) and (valid_indices[x-radius,y+radius,z]==1): #Top Diagonal Left
                        possible_value_to_mean.append(arr[x-radius,y+radius,z])
                    #Z+radius
                    if (arr[x+radius,y+radius,z+radius] !=0) and (valid_indices[x+radius,y+radius,z+radius]==1): #Top Diagonal Right
                        possible_value_to_mean.append(arr[x+radius,y+radius,z+radius])
                    if (arr[x-radius,y-radius,z+radius] !=0) and (valid_indices[x-radius,y-radius,z+radius]==1): #Bottom Diagonal Left
                        possible_value_to_mean.append(arr[x-radius,y-radius,z+radius])
                    if (arr[x+radius,y-radius,z+radius] !=0) and (valid_indices[x+radius,y-radius,z+radius]==1): #Bottom Diagonal Right
                        possible_value_to_mean.append(arr[x+radius,y-radius,z+radius])
                    if (arr[x-radius,y+radius,z+radius] !=0) and (valid_indices[x-radius,y+radius,z+radius]==1): #Top Diagonal Left
                        possible_value_to_mean.append(arr[x-radius,y+radius,z+radius])
                    if (arr[x+radius,y,z+radius] !=0) and (valid_indices[x+radius,y,z+radius]==1):
                        possible_value_to_mean.append(arr[x+radius,y,z+radius])
                    if (arr[x-radius,y,z+radius] !=0) and (valid_indices[x-radius,y,z+radius]==1):
                        possible_value_to_mean.append(arr[x-radius,y,z+radius])
                    if (arr[x,y+radius,z+radius] !=0) and (valid_indices[x,y+radius,z+radius]==1):
                        possible_value_to_mean.append(arr[x,y+radius,z+radius])
                    if (arr[x-radius,y+radius,z+radius] !=0) and (valid_indices[x-radius,y+radius,z+radius]==1):
                        possible_value_to_mean.append(arr[x-radius,y+radius,z+radius])
                    #Z-radius
                    if (arr[x+radius,y+radius,z-radius] !=0) and (valid_indices[x+radius,y+radius,z-radius]==1): #Top Diagonal Right
                        possible_value_to_mean.append(arr[x+radius,y+radius,z-radius])
                    if (arr[x-radius,y-radius,z-radius] !=0) and (valid_indices[x-radius,y-radius,z-radius]==1): #Bottom Diagonal Left
                        possible_value_to_mean.append(arr[x-radius,y-radius,z-radius])
                    if (arr[x+radius,y-radius,z-radius] !=0) and (valid_indices[x+radius,y-radius,z-radius]==1): #Bottom Diagonal Right
                        possible_value_to_mean.append(arr[x+radius,y-radius,z-radius])
                    if (arr[x-radius,y+radius,z-radius] !=0) and (valid_indices[x-radius,y+radius,z-radius]==1): #Top Diagonal Left
                        possible_value_to_mean.append(arr[x-radius,y+radius,z-radius])
                    if (arr[x+radius,y,z-radius] !=0) and (valid_indices[x+radius,y,z-radius]==1):
                        possible_value_to_mean.append(arr[x+radius,y,z-radius])
                    if (arr[x-radius,y,z-radius] !=0) and (valid_indices[x-radius,y,z-radius]==1):
                        possible_value_to_mean.append(arr[x-radius,y,z-radius])
                    if (arr[x,y+radius,z-radius] !=0) and (valid_indices[x,y+radius,z-radius]==1):
                        possible_value_to_mean.append(arr[x,y+radius,z-radius])
                    if (arr[x-radius,y+radius,z-radius] !=0) and (valid_indices[x-radius,y+radius,z-radius]==1):
                        possible_value_to_mean.append(arr[x-radius,y+radius,z-radius])
        
                    arr[x,y,z] = np.average(possible_value_to_mean)

    print("Number of Nans in the kidney: " + str(np.count_nonzero(np.isnan(arr))))

    # plt.imshow(arr[:,:,90])
    # plt.show()

    nan_indices = np.isnan(arr)

    # Find indices of non-NaN and non-zero values
    non_nan_indices = ~nan_indices
    nonzero_indices = (arr != 0)

    # Combine non-NaN and non-zero indices
    valid_indices = np.logical_and(non_nan_indices, nonzero_indices)

    for z in range (arr.shape[2]):
        for y in range(arr.shape[1]):
            for x in range(arr.shape[0]):
                if np.isnan(arr[x,y,z]):
                    possible_value_to_mean =[]
                    if (arr[x+radius,y,z] !=0) and (valid_indices[x+radius,y,z]==1): #Right 
                        possible_value_to_mean.append(arr[x+radius,y,z])
                    if (arr[x-radius,y,z] !=0) and (valid_indices[x-radius,y,z]==1): #Left
                        possible_value_to_mean.append(arr[x-radius,y,z])
                    if (arr[x,y+radius,z] !=0) and (valid_indices[x,y+radius,z]==1): #Top
                        possible_value_to_mean.append(arr[x,y+radius,z])
                    if (arr[x,y-radius,z] !=0) and (valid_indices[x,y-radius,z]==1): #Bottom
                        possible_value_to_mean.append(arr[x,y-radius,z])
                    if (arr[x,y,z+radius] !=0) and (valid_indices[x,y,z+radius]==1): #Front
                        possible_value_to_mean.append(arr[x,y,z+radius])
                    if (arr[x,y,z-radius] !=0) and (valid_indices[x,y,z-radius]==1): #Back
                        possible_value_to_mean.append(arr[x,y,z-radius])                   
                    if (arr[x+radius,y+radius,z] !=0) and (valid_indices[x+radius,y+radius,z]==1): #Top Diagonal Right
                        possible_value_to_mean.append(arr[x+radius,y+radius,z])
                    if (arr[x-radius,y-radius,z] !=0) and (valid_indices[x-radius,y-radius,z]==1): #Bottom Diagonal Left
                        possible_value_to_mean.append(arr[x-radius,y-radius,z])
                    if (arr[x+radius,y-radius,z] !=0) and (valid_indices[x+radius,y-radius,z]==1): #Bottom Diagonal Right
                        possible_value_to_mean.append(arr[x+radius,y-radius,z])
                    if (arr[x-radius,y+radius,z] !=0) and (valid_indices[x-radius,y+radius,z]==1): #Top Diagonal Left
                        possible_value_to_mean.append(arr[x-radius,y+radius,z])
                    #Z+radius
                    if (arr[x+radius,y+radius,z+radius] !=0) and (valid_indices[x+radius,y+radius,z+radius]==1): #Top Diagonal Right
                        possible_value_to_mean.append(arr[x+radius,y+radius,z+radius])
                    if (arr[x-radius,y-radius,z+radius] !=0) and (valid_indices[x-radius,y-radius,z+radius]==1): #Bottom Diagonal Left
                        possible_value_to_mean.append(arr[x-radius,y-radius,z+radius])
                    if (arr[x+radius,y-radius,z+radius] !=0) and (valid_indices[x+radius,y-radius,z+radius]==1): #Bottom Diagonal Right
                        possible_value_to_mean.append(arr[x+radius,y-radius,z+radius])
                    if (arr[x-radius,y+radius,z+radius] !=0) and (valid_indices[x-radius,y+radius,z+radius]==1): #Top Diagonal Left
                        possible_value_to_mean.append(arr[x-radius,y+radius,z+radius])
                    if (arr[x+radius,y,z+radius] !=0) and (valid_indices[x+radius,y,z+radius]==1):
                        possible_value_to_mean.append(arr[x+radius,y,z+radius])
                    if (arr[x-radius,y,z+radius] !=0) and (valid_indices[x-radius,y,z+radius]==1):
                        possible_value_to_mean.append(arr[x-radius,y,z+radius])
                    if (arr[x,y+radius,z+radius] !=0) and (valid_indices[x,y+radius,z+radius]==1):
                        possible_value_to_mean.append(arr[x,y+radius,z+radius])
                    if (arr[x-radius,y+radius,z+radius] !=0) and (valid_indices[x-radius,y+radius,z+radius]==1):
                        possible_value_to_mean.append(arr[x-radius,y+radius,z+radius])
                    #Z-radius
                    if (arr[x+radius,y+radius,z-radius] !=0) and (valid_indices[x+radius,y+radius,z-radius]==1): #Top Diagonal Right
                        possible_value_to_mean.append(arr[x+radius,y+radius,z-radius])
                    if (arr[x-radius,y-radius,z-radius] !=0) and (valid_indices[x-radius,y-radius,z-radius]==1): #Bottom Diagonal Left
                        possible_value_to_mean.append(arr[x-radius,y-radius,z-radius])
                    if (arr[x+radius,y-radius,z-radius] !=0) and (valid_indices[x+radius,y-radius,z-radius]==1): #Bottom Diagonal Right
                        possible_value_to_mean.append(arr[x+radius,y-radius,z-radius])
                    if (arr[x-radius,y+radius,z-radius] !=0) and (valid_indices[x-radius,y+radius,z-radius]==1): #Top Diagonal Left
                        possible_value_to_mean.append(arr[x-radius,y+radius,z-radius])
                    if (arr[x+radius,y,z-radius] !=0) and (valid_indices[x+radius,y,z-radius]==1):
                        possible_value_to_mean.append(arr[x+radius,y,z-radius])
                    if (arr[x-radius,y,z-radius] !=0) and (valid_indices[x-radius,y,z-radius]==1):
                        possible_value_to_mean.append(arr[x-radius,y,z-radius])
                    if (arr[x,y+radius,z-radius] !=0) and (valid_indices[x,y+radius,z-radius]==1):
                        possible_value_to_mean.append(arr[x,y+radius,z-radius])
                    if (arr[x-radius,y+radius,z-radius] !=0) and (valid_indices[x-radius,y+radius,z-radius]==1):
                        possible_value_to_mean.append(arr[x-radius,y+radius,z-radius])
        
                    arr[x,y,z] = np.average(possible_value_to_mean)

    print("Number of Nans in the kidney: " + str(np.count_nonzero(np.isnan(arr))))

    # plt.imshow(arr[:,:,90])
    # plt.show()

    output_array = arr

    return output_array