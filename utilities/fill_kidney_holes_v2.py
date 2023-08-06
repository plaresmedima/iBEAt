import numpy as np
import matplotlib.pyplot as plt

def main(input_array, array_mask, slice):

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

    plt.imshow(np.transpose(arr[:,:,slice]))
    plt.show()

    valid_indices = np.logical_and(non_nan_indices, nonzero_indices)

    # Loop through every pixel and if it finds a nan look for all 26 neighbours. At the end replace the nan pixel with the average of the valid neighbours  
    radius = 1
    nan_left = np.count_nonzero(np.isnan(arr))
    nan_end_of_loop = []
    nan_end_of_loop.append(nan_left)
    i =0

    while nan_left != 0:

        for z in range (int(arr.shape[2]/2),arr.shape[2]):
            for x in range(arr.shape[0]):
                for y in range(arr.shape[1]):
                    if np.isnan(arr[x,y,z]):
                        possible_value_to_mean =[]
                        try:
                            if (arr[x+radius,y,z] !=0) and (valid_indices[x+radius,y,z]==1): #Right 
                                possible_value_to_mean.append(arr[x+radius,y,z])
                        except:
                            print("out of bounds")
                        try:
                            if (arr[x-radius,y,z] !=0) and (valid_indices[x-radius,y,z]==1): #Left
                                possible_value_to_mean.append(arr[x-radius,y,z])
                        except:
                            print("out of bounds")
                        try:        
                            if (arr[x,y+radius,z] !=0) and (valid_indices[x,y+radius,z]==1): #Top
                                possible_value_to_mean.append(arr[x,y+radius,z])
                        except:
                            print("out of bounds")
                        try:   
                            if (arr[x,y-radius,z] !=0) and (valid_indices[x,y-radius,z]==1): #Bottom
                                possible_value_to_mean.append(arr[x,y-radius,z])
                        except:
                            print("out of bounds")
                        try:
                            if (arr[x,y,z+radius] !=0) and (valid_indices[x,y,z+radius]==1): #Front
                                possible_value_to_mean.append(arr[x,y,z+radius])
                        except:
                            print("out of bounds")
                        try:
                            if (arr[x,y,z-radius] !=0) and (valid_indices[x,y,z-radius]==1): #Back
                                possible_value_to_mean.append(arr[x,y,z-radius]) 
                        except:
                            print("out of bounds")
                        try:             
                            if (arr[x+radius,y+radius,z] !=0) and (valid_indices[x+radius,y+radius,z]==1): #Top Diagonal Right
                                possible_value_to_mean.append(arr[x+radius,y+radius,z])
                        except:
                            print("out of bounds")
                        try:
                            if (arr[x-radius,y-radius,z] !=0) and (valid_indices[x-radius,y-radius,z]==1): #Bottom Diagonal Left
                                possible_value_to_mean.append(arr[x-radius,y-radius,z])
                        except:
                            print("out of bounds")
                        try:
                            if (arr[x+radius,y-radius,z] !=0) and (valid_indices[x+radius,y-radius,z]==1): #Bottom Diagonal Right
                                possible_value_to_mean.append(arr[x+radius,y-radius,z])
                        except:
                            print("out of bounds")
                        try:
                            if (arr[x-radius,y+radius,z] !=0) and (valid_indices[x-radius,y+radius,z]==1): #Top Diagonal Left
                                possible_value_to_mean.append(arr[x-radius,y+radius,z])
                        except:
                            print("out of bounds")

                        #Z+radius

                        try:
                            if (arr[x+radius,y+radius,z+radius] !=0) and (valid_indices[x+radius,y+radius,z+radius]==1): #Top Diagonal Right
                                possible_value_to_mean.append(arr[x+radius,y+radius,z+radius])
                        except:
                            print("out of bounds")
                        try:
                            if (arr[x-radius,y-radius,z+radius] !=0) and (valid_indices[x-radius,y-radius,z+radius]==1): #Bottom Diagonal Left
                                possible_value_to_mean.append(arr[x-radius,y-radius,z+radius])
                        except:
                            print("out of bounds")
                        try:
                            if (arr[x+radius,y-radius,z+radius] !=0) and (valid_indices[x+radius,y-radius,z+radius]==1): #Bottom Diagonal Right
                                possible_value_to_mean.append(arr[x+radius,y-radius,z+radius])
                        except:
                            print("out of bounds")
                        try:    
                            if (arr[x-radius,y+radius,z+radius] !=0) and (valid_indices[x-radius,y+radius,z+radius]==1): #Top Diagonal Left
                                possible_value_to_mean.append(arr[x-radius,y+radius,z+radius])
                        except:
                            print("out of bounds")
                        try:        
                            if (arr[x+radius,y,z+radius] !=0) and (valid_indices[x+radius,y,z+radius]==1):
                                possible_value_to_mean.append(arr[x+radius,y,z+radius])
                        except:
                            print("out of bounds")                    
                        try:
                            if (arr[x-radius,y,z+radius] !=0) and (valid_indices[x-radius,y,z+radius]==1):
                                possible_value_to_mean.append(arr[x-radius,y,z+radius])
                        except:
                            print("out of bounds")
                        try:
                            if (arr[x,y+radius,z+radius] !=0) and (valid_indices[x,y+radius,z+radius]==1):
                                possible_value_to_mean.append(arr[x,y+radius,z+radius])
                        except:
                            print("out of bounds")
                        try:
                            if (arr[x-radius,y+radius,z+radius] !=0) and (valid_indices[x-radius,y+radius,z+radius]==1):
                                possible_value_to_mean.append(arr[x-radius,y+radius,z+radius])
                        except:
                            print("out of bounds")
                        #Z-radius
                        try:
                            if (arr[x+radius,y+radius,z-radius] !=0) and (valid_indices[x+radius,y+radius,z-radius]==1): #Top Diagonal Right
                                possible_value_to_mean.append(arr[x+radius,y+radius,z-radius])
                        except:
                            print("out of bounds")
                        try:
                            if (arr[x-radius,y-radius,z-radius] !=0) and (valid_indices[x-radius,y-radius,z-radius]==1): #Bottom Diagonal Left
                                possible_value_to_mean.append(arr[x-radius,y-radius,z-radius])
                        except:
                            print("out of bounds")
                        try:        
                            if (arr[x+radius,y-radius,z-radius] !=0) and (valid_indices[x+radius,y-radius,z-radius]==1): #Bottom Diagonal Right
                                possible_value_to_mean.append(arr[x+radius,y-radius,z-radius])
                        except:
                            print("out of bounds")
                        try:
                            if (arr[x-radius,y+radius,z-radius] !=0) and (valid_indices[x-radius,y+radius,z-radius]==1): #Top Diagonal Left
                                possible_value_to_mean.append(arr[x-radius,y+radius,z-radius])
                        except:
                            print("out of bounds")
                        try:
                            if (arr[x+radius,y,z-radius] !=0) and (valid_indices[x+radius,y,z-radius]==1):
                                possible_value_to_mean.append(arr[x+radius,y,z-radius])
                        except:
                            print("out of bounds")
                        try:
                            if (arr[x-radius,y,z-radius] !=0) and (valid_indices[x-radius,y,z-radius]==1):
                                possible_value_to_mean.append(arr[x-radius,y,z-radius])
                        except:
                            print("out of bounds")
                        try:
                            if (arr[x,y+radius,z-radius] !=0) and (valid_indices[x,y+radius,z-radius]==1):
                                possible_value_to_mean.append(arr[x,y+radius,z-radius])
                        except:
                            print("out of bounds")
                        try:
                            if (arr[x-radius,y+radius,z-radius] !=0) and (valid_indices[x-radius,y+radius,z-radius]==1):
                                possible_value_to_mean.append(arr[x-radius,y+radius,z-radius])
                        except:
                            print("out of bounds")

                        arr[x,y,z] = np.average(possible_value_to_mean)


        plt.imshow(np.transpose(arr[:,:,slice]))
        plt.show()

        print("Number of Nans in the kidney: " + str(np.count_nonzero(np.isnan(arr))))

        nan_indices = np.isnan(arr)

        # Find indices of non-NaN and non-zero values
        non_nan_indices = ~nan_indices
        nonzero_indices = (arr != 0)

        # Combine non-NaN and non-zero indices
        valid_indices = np.logical_and(non_nan_indices, nonzero_indices)
        

        for z in range (int(arr.shape[2]/2),-1,-1):
            for x in range(arr.shape[0]):
                for y in range(arr.shape[1]):
                    if np.isnan(arr[x,y,z]):
                        possible_value_to_mean =[]
                        try:
                            if (arr[x+radius,y,z] !=0) and (valid_indices[x+radius,y,z]==1): #Right 
                                possible_value_to_mean.append(arr[x+radius,y,z])
                        except:
                            print("out of bounds")
                        try:
                            if (arr[x-radius,y,z] !=0) and (valid_indices[x-radius,y,z]==1): #Left
                                possible_value_to_mean.append(arr[x-radius,y,z])
                        except:
                            print("out of bounds")
                        try:        
                            if (arr[x,y+radius,z] !=0) and (valid_indices[x,y+radius,z]==1): #Top
                                possible_value_to_mean.append(arr[x,y+radius,z])
                        except:
                            print("out of bounds")
                        try:   
                            if (arr[x,y-radius,z] !=0) and (valid_indices[x,y-radius,z]==1): #Bottom
                                possible_value_to_mean.append(arr[x,y-radius,z])
                        except:
                            print("out of bounds")
                        try:
                            if (arr[x,y,z+radius] !=0) and (valid_indices[x,y,z+radius]==1): #Front
                                possible_value_to_mean.append(arr[x,y,z+radius])
                        except:
                            print("out of bounds")
                        try:
                            if (arr[x,y,z-radius] !=0) and (valid_indices[x,y,z-radius]==1): #Back
                                possible_value_to_mean.append(arr[x,y,z-radius]) 
                        except:
                            print("out of bounds")
                        try:             
                            if (arr[x+radius,y+radius,z] !=0) and (valid_indices[x+radius,y+radius,z]==1): #Top Diagonal Right
                                possible_value_to_mean.append(arr[x+radius,y+radius,z])
                        except:
                            print("out of bounds")
                        try:
                            if (arr[x-radius,y-radius,z] !=0) and (valid_indices[x-radius,y-radius,z]==1): #Bottom Diagonal Left
                                possible_value_to_mean.append(arr[x-radius,y-radius,z])
                        except:
                            print("out of bounds")
                        try:
                            if (arr[x+radius,y-radius,z] !=0) and (valid_indices[x+radius,y-radius,z]==1): #Bottom Diagonal Right
                                possible_value_to_mean.append(arr[x+radius,y-radius,z])
                        except:
                            print("out of bounds")
                        try:
                            if (arr[x-radius,y+radius,z] !=0) and (valid_indices[x-radius,y+radius,z]==1): #Top Diagonal Left
                                possible_value_to_mean.append(arr[x-radius,y+radius,z])
                        except:
                            print("out of bounds")

                        #Z+radius

                        try:
                            if (arr[x+radius,y+radius,z+radius] !=0) and (valid_indices[x+radius,y+radius,z+radius]==1): #Top Diagonal Right
                                possible_value_to_mean.append(arr[x+radius,y+radius,z+radius])
                        except:
                            print("out of bounds")
                        try:
                            if (arr[x-radius,y-radius,z+radius] !=0) and (valid_indices[x-radius,y-radius,z+radius]==1): #Bottom Diagonal Left
                                possible_value_to_mean.append(arr[x-radius,y-radius,z+radius])
                        except:
                            print("out of bounds")
                        try:
                            if (arr[x+radius,y-radius,z+radius] !=0) and (valid_indices[x+radius,y-radius,z+radius]==1): #Bottom Diagonal Right
                                possible_value_to_mean.append(arr[x+radius,y-radius,z+radius])
                        except:
                            print("out of bounds")
                        try:    
                            if (arr[x-radius,y+radius,z+radius] !=0) and (valid_indices[x-radius,y+radius,z+radius]==1): #Top Diagonal Left
                                possible_value_to_mean.append(arr[x-radius,y+radius,z+radius])
                        except:
                            print("out of bounds")
                        try:        
                            if (arr[x+radius,y,z+radius] !=0) and (valid_indices[x+radius,y,z+radius]==1):
                                possible_value_to_mean.append(arr[x+radius,y,z+radius])
                        except:
                            print("out of bounds")                    
                        try:
                            if (arr[x-radius,y,z+radius] !=0) and (valid_indices[x-radius,y,z+radius]==1):
                                possible_value_to_mean.append(arr[x-radius,y,z+radius])
                        except:
                            print("out of bounds")
                        try:
                            if (arr[x,y+radius,z+radius] !=0) and (valid_indices[x,y+radius,z+radius]==1):
                                possible_value_to_mean.append(arr[x,y+radius,z+radius])
                        except:
                            print("out of bounds")
                        try:
                            if (arr[x-radius,y+radius,z+radius] !=0) and (valid_indices[x-radius,y+radius,z+radius]==1):
                                possible_value_to_mean.append(arr[x-radius,y+radius,z+radius])
                        except:
                            print("out of bounds")
                        #Z-radius
                        try:
                            if (arr[x+radius,y+radius,z-radius] !=0) and (valid_indices[x+radius,y+radius,z-radius]==1): #Top Diagonal Right
                                possible_value_to_mean.append(arr[x+radius,y+radius,z-radius])
                        except:
                            print("out of bounds")
                        try:
                            if (arr[x-radius,y-radius,z-radius] !=0) and (valid_indices[x-radius,y-radius,z-radius]==1): #Bottom Diagonal Left
                                possible_value_to_mean.append(arr[x-radius,y-radius,z-radius])
                        except:
                            print("out of bounds")
                        try:        
                            if (arr[x+radius,y-radius,z-radius] !=0) and (valid_indices[x+radius,y-radius,z-radius]==1): #Bottom Diagonal Right
                                possible_value_to_mean.append(arr[x+radius,y-radius,z-radius])
                        except:
                            print("out of bounds")
                        try:
                            if (arr[x-radius,y+radius,z-radius] !=0) and (valid_indices[x-radius,y+radius,z-radius]==1): #Top Diagonal Left
                                possible_value_to_mean.append(arr[x-radius,y+radius,z-radius])
                        except:
                            print("out of bounds")
                        try:
                            if (arr[x+radius,y,z-radius] !=0) and (valid_indices[x+radius,y,z-radius]==1):
                                possible_value_to_mean.append(arr[x+radius,y,z-radius])
                        except:
                            print("out of bounds")
                        try:
                            if (arr[x-radius,y,z-radius] !=0) and (valid_indices[x-radius,y,z-radius]==1):
                                possible_value_to_mean.append(arr[x-radius,y,z-radius])
                        except:
                            print("out of bounds")
                        try:
                            if (arr[x,y+radius,z-radius] !=0) and (valid_indices[x,y+radius,z-radius]==1):
                                possible_value_to_mean.append(arr[x,y+radius,z-radius])
                        except:
                            print("out of bounds")
                        try:
                            if (arr[x-radius,y+radius,z-radius] !=0) and (valid_indices[x-radius,y+radius,z-radius]==1):
                                possible_value_to_mean.append(arr[x-radius,y+radius,z-radius])
                        except:
                            print("out of bounds")

                        arr[x,y,z] = np.average(possible_value_to_mean)

        for z in range (int(arr.shape[2]/2),arr.shape[2]):
            for y in range(arr.shape[1]):
                for x in range(arr.shape[0]):
                    if np.isnan(arr[x,y,z]):
                        possible_value_to_mean =[]
                        try:
                            if (arr[x+radius,y,z] !=0) and (valid_indices[x+radius,y,z]==1): #Right 
                                possible_value_to_mean.append(arr[x+radius,y,z])
                        except:
                            print("out of bounds")
                        try:
                            if (arr[x-radius,y,z] !=0) and (valid_indices[x-radius,y,z]==1): #Left
                                possible_value_to_mean.append(arr[x-radius,y,z])
                        except:
                            print("out of bounds")
                        try:        
                            if (arr[x,y+radius,z] !=0) and (valid_indices[x,y+radius,z]==1): #Top
                                possible_value_to_mean.append(arr[x,y+radius,z])
                        except:
                            print("out of bounds")
                        try:   
                            if (arr[x,y-radius,z] !=0) and (valid_indices[x,y-radius,z]==1): #Bottom
                                possible_value_to_mean.append(arr[x,y-radius,z])
                        except:
                            print("out of bounds")
                        try:
                            if (arr[x,y,z+radius] !=0) and (valid_indices[x,y,z+radius]==1): #Front
                                possible_value_to_mean.append(arr[x,y,z+radius])
                        except:
                            print("out of bounds")
                        try:
                            if (arr[x,y,z-radius] !=0) and (valid_indices[x,y,z-radius]==1): #Back
                                possible_value_to_mean.append(arr[x,y,z-radius]) 
                        except:
                            print("out of bounds")
                        try:             
                            if (arr[x+radius,y+radius,z] !=0) and (valid_indices[x+radius,y+radius,z]==1): #Top Diagonal Right
                                possible_value_to_mean.append(arr[x+radius,y+radius,z])
                        except:
                            print("out of bounds")
                        try:
                            if (arr[x-radius,y-radius,z] !=0) and (valid_indices[x-radius,y-radius,z]==1): #Bottom Diagonal Left
                                possible_value_to_mean.append(arr[x-radius,y-radius,z])
                        except:
                            print("out of bounds")
                        try:
                            if (arr[x+radius,y-radius,z] !=0) and (valid_indices[x+radius,y-radius,z]==1): #Bottom Diagonal Right
                                possible_value_to_mean.append(arr[x+radius,y-radius,z])
                        except:
                            print("out of bounds")
                        try:
                            if (arr[x-radius,y+radius,z] !=0) and (valid_indices[x-radius,y+radius,z]==1): #Top Diagonal Left
                                possible_value_to_mean.append(arr[x-radius,y+radius,z])
                        except:
                            print("out of bounds")

                        #Z+radius

                        try:
                            if (arr[x+radius,y+radius,z+radius] !=0) and (valid_indices[x+radius,y+radius,z+radius]==1): #Top Diagonal Right
                                possible_value_to_mean.append(arr[x+radius,y+radius,z+radius])
                        except:
                            print("out of bounds")
                        try:
                            if (arr[x-radius,y-radius,z+radius] !=0) and (valid_indices[x-radius,y-radius,z+radius]==1): #Bottom Diagonal Left
                                possible_value_to_mean.append(arr[x-radius,y-radius,z+radius])
                        except:
                            print("out of bounds")
                        try:
                            if (arr[x+radius,y-radius,z+radius] !=0) and (valid_indices[x+radius,y-radius,z+radius]==1): #Bottom Diagonal Right
                                possible_value_to_mean.append(arr[x+radius,y-radius,z+radius])
                        except:
                            print("out of bounds")
                        try:    
                            if (arr[x-radius,y+radius,z+radius] !=0) and (valid_indices[x-radius,y+radius,z+radius]==1): #Top Diagonal Left
                                possible_value_to_mean.append(arr[x-radius,y+radius,z+radius])
                        except:
                            print("out of bounds")
                        try:        
                            if (arr[x+radius,y,z+radius] !=0) and (valid_indices[x+radius,y,z+radius]==1):
                                possible_value_to_mean.append(arr[x+radius,y,z+radius])
                        except:
                            print("out of bounds")                    
                        try:
                            if (arr[x-radius,y,z+radius] !=0) and (valid_indices[x-radius,y,z+radius]==1):
                                possible_value_to_mean.append(arr[x-radius,y,z+radius])
                        except:
                            print("out of bounds")
                        try:
                            if (arr[x,y+radius,z+radius] !=0) and (valid_indices[x,y+radius,z+radius]==1):
                                possible_value_to_mean.append(arr[x,y+radius,z+radius])
                        except:
                            print("out of bounds")
                        try:
                            if (arr[x-radius,y+radius,z+radius] !=0) and (valid_indices[x-radius,y+radius,z+radius]==1):
                                possible_value_to_mean.append(arr[x-radius,y+radius,z+radius])
                        except:
                            print("out of bounds")
                        #Z-radius
                        try:
                            if (arr[x+radius,y+radius,z-radius] !=0) and (valid_indices[x+radius,y+radius,z-radius]==1): #Top Diagonal Right
                                possible_value_to_mean.append(arr[x+radius,y+radius,z-radius])
                        except:
                            print("out of bounds")
                        try:
                            if (arr[x-radius,y-radius,z-radius] !=0) and (valid_indices[x-radius,y-radius,z-radius]==1): #Bottom Diagonal Left
                                possible_value_to_mean.append(arr[x-radius,y-radius,z-radius])
                        except:
                            print("out of bounds")
                        try:        
                            if (arr[x+radius,y-radius,z-radius] !=0) and (valid_indices[x+radius,y-radius,z-radius]==1): #Bottom Diagonal Right
                                possible_value_to_mean.append(arr[x+radius,y-radius,z-radius])
                        except:
                            print("out of bounds")
                        try:
                            if (arr[x-radius,y+radius,z-radius] !=0) and (valid_indices[x-radius,y+radius,z-radius]==1): #Top Diagonal Left
                                possible_value_to_mean.append(arr[x-radius,y+radius,z-radius])
                        except:
                            print("out of bounds")
                        try:
                            if (arr[x+radius,y,z-radius] !=0) and (valid_indices[x+radius,y,z-radius]==1):
                                possible_value_to_mean.append(arr[x+radius,y,z-radius])
                        except:
                            print("out of bounds")
                        try:
                            if (arr[x-radius,y,z-radius] !=0) and (valid_indices[x-radius,y,z-radius]==1):
                                possible_value_to_mean.append(arr[x-radius,y,z-radius])
                        except:
                            print("out of bounds")
                        try:
                            if (arr[x,y+radius,z-radius] !=0) and (valid_indices[x,y+radius,z-radius]==1):
                                possible_value_to_mean.append(arr[x,y+radius,z-radius])
                        except:
                            print("out of bounds")
                        try:
                            if (arr[x-radius,y+radius,z-radius] !=0) and (valid_indices[x-radius,y+radius,z-radius]==1):
                                possible_value_to_mean.append(arr[x-radius,y+radius,z-radius])
                        except:
                            print("out of bounds")

                        arr[x,y,z] = np.average(possible_value_to_mean)


        plt.imshow(np.transpose(arr[:,:,slice]))
        plt.show()

        print("Number of Nans in the kidney: " + str(np.count_nonzero(np.isnan(arr))))

        nan_indices = np.isnan(arr)

        # Find indices of non-NaN and non-zero values
        non_nan_indices = ~nan_indices
        nonzero_indices = (arr != 0)

        # Combine non-NaN and non-zero indices
        valid_indices = np.logical_and(non_nan_indices, nonzero_indices)
        

        for z in range (int(arr.shape[2]/2),-1,-1):
            for y in range(arr.shape[1]):
                for x in range(arr.shape[0]):
                    if np.isnan(arr[x,y,z]):
                        possible_value_to_mean =[]
                        try:
                            if (arr[x+radius,y,z] !=0) and (valid_indices[x+radius,y,z]==1): #Right 
                                possible_value_to_mean.append(arr[x+radius,y,z])
                        except:
                            print("out of bounds")
                        try:
                            if (arr[x-radius,y,z] !=0) and (valid_indices[x-radius,y,z]==1): #Left
                                possible_value_to_mean.append(arr[x-radius,y,z])
                        except:
                            print("out of bounds")
                        try:        
                            if (arr[x,y+radius,z] !=0) and (valid_indices[x,y+radius,z]==1): #Top
                                possible_value_to_mean.append(arr[x,y+radius,z])
                        except:
                            print("out of bounds")
                        try:   
                            if (arr[x,y-radius,z] !=0) and (valid_indices[x,y-radius,z]==1): #Bottom
                                possible_value_to_mean.append(arr[x,y-radius,z])
                        except:
                            print("out of bounds")
                        try:
                            if (arr[x,y,z+radius] !=0) and (valid_indices[x,y,z+radius]==1): #Front
                                possible_value_to_mean.append(arr[x,y,z+radius])
                        except:
                            print("out of bounds")
                        try:
                            if (arr[x,y,z-radius] !=0) and (valid_indices[x,y,z-radius]==1): #Back
                                possible_value_to_mean.append(arr[x,y,z-radius]) 
                        except:
                            print("out of bounds")
                        try:             
                            if (arr[x+radius,y+radius,z] !=0) and (valid_indices[x+radius,y+radius,z]==1): #Top Diagonal Right
                                possible_value_to_mean.append(arr[x+radius,y+radius,z])
                        except:
                            print("out of bounds")
                        try:
                            if (arr[x-radius,y-radius,z] !=0) and (valid_indices[x-radius,y-radius,z]==1): #Bottom Diagonal Left
                                possible_value_to_mean.append(arr[x-radius,y-radius,z])
                        except:
                            print("out of bounds")
                        try:
                            if (arr[x+radius,y-radius,z] !=0) and (valid_indices[x+radius,y-radius,z]==1): #Bottom Diagonal Right
                                possible_value_to_mean.append(arr[x+radius,y-radius,z])
                        except:
                            print("out of bounds")
                        try:
                            if (arr[x-radius,y+radius,z] !=0) and (valid_indices[x-radius,y+radius,z]==1): #Top Diagonal Left
                                possible_value_to_mean.append(arr[x-radius,y+radius,z])
                        except:
                            print("out of bounds")

                        #Z+radius

                        try:
                            if (arr[x+radius,y+radius,z+radius] !=0) and (valid_indices[x+radius,y+radius,z+radius]==1): #Top Diagonal Right
                                possible_value_to_mean.append(arr[x+radius,y+radius,z+radius])
                        except:
                            print("out of bounds")
                        try:
                            if (arr[x-radius,y-radius,z+radius] !=0) and (valid_indices[x-radius,y-radius,z+radius]==1): #Bottom Diagonal Left
                                possible_value_to_mean.append(arr[x-radius,y-radius,z+radius])
                        except:
                            print("out of bounds")
                        try:
                            if (arr[x+radius,y-radius,z+radius] !=0) and (valid_indices[x+radius,y-radius,z+radius]==1): #Bottom Diagonal Right
                                possible_value_to_mean.append(arr[x+radius,y-radius,z+radius])
                        except:
                            print("out of bounds")
                        try:    
                            if (arr[x-radius,y+radius,z+radius] !=0) and (valid_indices[x-radius,y+radius,z+radius]==1): #Top Diagonal Left
                                possible_value_to_mean.append(arr[x-radius,y+radius,z+radius])
                        except:
                            print("out of bounds")
                        try:        
                            if (arr[x+radius,y,z+radius] !=0) and (valid_indices[x+radius,y,z+radius]==1):
                                possible_value_to_mean.append(arr[x+radius,y,z+radius])
                        except:
                            print("out of bounds")                    
                        try:
                            if (arr[x-radius,y,z+radius] !=0) and (valid_indices[x-radius,y,z+radius]==1):
                                possible_value_to_mean.append(arr[x-radius,y,z+radius])
                        except:
                            print("out of bounds")
                        try:
                            if (arr[x,y+radius,z+radius] !=0) and (valid_indices[x,y+radius,z+radius]==1):
                                possible_value_to_mean.append(arr[x,y+radius,z+radius])
                        except:
                            print("out of bounds")
                        try:
                            if (arr[x-radius,y+radius,z+radius] !=0) and (valid_indices[x-radius,y+radius,z+radius]==1):
                                possible_value_to_mean.append(arr[x-radius,y+radius,z+radius])
                        except:
                            print("out of bounds")
                        #Z-radius
                        try:
                            if (arr[x+radius,y+radius,z-radius] !=0) and (valid_indices[x+radius,y+radius,z-radius]==1): #Top Diagonal Right
                                possible_value_to_mean.append(arr[x+radius,y+radius,z-radius])
                        except:
                            print("out of bounds")
                        try:
                            if (arr[x-radius,y-radius,z-radius] !=0) and (valid_indices[x-radius,y-radius,z-radius]==1): #Bottom Diagonal Left
                                possible_value_to_mean.append(arr[x-radius,y-radius,z-radius])
                        except:
                            print("out of bounds")
                        try:        
                            if (arr[x+radius,y-radius,z-radius] !=0) and (valid_indices[x+radius,y-radius,z-radius]==1): #Bottom Diagonal Right
                                possible_value_to_mean.append(arr[x+radius,y-radius,z-radius])
                        except:
                            print("out of bounds")
                        try:
                            if (arr[x-radius,y+radius,z-radius] !=0) and (valid_indices[x-radius,y+radius,z-radius]==1): #Top Diagonal Left
                                possible_value_to_mean.append(arr[x-radius,y+radius,z-radius])
                        except:
                            print("out of bounds")
                        try:
                            if (arr[x+radius,y,z-radius] !=0) and (valid_indices[x+radius,y,z-radius]==1):
                                possible_value_to_mean.append(arr[x+radius,y,z-radius])
                        except:
                            print("out of bounds")
                        try:
                            if (arr[x-radius,y,z-radius] !=0) and (valid_indices[x-radius,y,z-radius]==1):
                                possible_value_to_mean.append(arr[x-radius,y,z-radius])
                        except:
                            print("out of bounds")
                        try:
                            if (arr[x,y+radius,z-radius] !=0) and (valid_indices[x,y+radius,z-radius]==1):
                                possible_value_to_mean.append(arr[x,y+radius,z-radius])
                        except:
                            print("out of bounds")
                        try:
                            if (arr[x-radius,y+radius,z-radius] !=0) and (valid_indices[x-radius,y+radius,z-radius]==1):
                                possible_value_to_mean.append(arr[x-radius,y+radius,z-radius])
                        except:
                            print("out of bounds")

                        arr[x,y,z] = np.average(possible_value_to_mean)

        plt.imshow(np.transpose(arr[:,:,slice]))
        plt.show()

        print("Number of Nans in the kidney: " + str(np.count_nonzero(np.isnan(arr))))

        nan_indices = np.isnan(arr)

        # Find indices of non-NaN and non-zero values
        non_nan_indices = ~nan_indices
        nonzero_indices = (arr != 0)

        # Combine non-NaN and non-zero indices
        valid_indices = np.logical_and(non_nan_indices, nonzero_indices)

        nan_left = np.count_nonzero(np.isnan(arr))

        output_array = arr

        i=i+1
        nan_end_of_loop.append(nan_left)

        if nan_end_of_loop[i] == nan_end_of_loop[i-1]:
            break
        
    return output_array