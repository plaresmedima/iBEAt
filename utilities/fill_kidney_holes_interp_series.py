import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

def main(input_series, mask_series):

    input_array, input_header = input_series.array(['SliceLocation','AcquisitionTime'], pixels_first=True)
    mask_array , mask_header  = mask_series.array (['SliceLocation'],                   pixels_first=True)

    input_array = np.squeeze(input_array)
    mask_array = np.squeeze(mask_array)

    array_Data_masked = input_array*mask_array
    array_Data_masked[array_Data_masked<=1] = 0 #to guarantee 0 are zeros and no 10^-13

    #Detect the regions in the kidney to fill with values and replace by nans
    a = np.zeros(mask_array.shape) 
    a [mask_array==1] = 1

    b = np.zeros(array_Data_masked.shape) 
    b [array_Data_masked==0] = 1

    c= np.logical_and(a, b)

    array_Data_masked[c==1] = np.nan

    arr = array_Data_masked

    # Create a masked array to ignore NaNs and zeros during interpolation
    masked_arr_zero = np.ma.masked_equal(arr, 0)
    masked_arr = np.ma.masked_invalid(masked_arr_zero)

    # Create arrays of coordinates for interpolation
    coords = np.indices(arr.shape)
    x, y, z = coords[0], coords[1], coords[2]

    # Find the valid (non-masked) elements and corresponding coordinates
    valid_points = np.column_stack([x[~masked_arr.mask], y[~masked_arr.mask], z[~masked_arr.mask]])
    data_valid = masked_arr[~masked_arr.mask]

    # Interpolate using griddata
    interpolated_values = griddata(valid_points, data_valid, (x, y, z), method='linear')

    # Fill NaN values in the original array with interpolated values
    arr[np.isnan(arr)] = interpolated_values[np.isnan(arr)]

    output_array = interpolated_values

    output_array = output_array*mask_array

    output_series = input_series.SeriesDescription + "Coreg_filled_" + input_series.SeriesDescription
    output_series = input_series.new_series(SeriesDescription="Coreg_filled")
    output_series.set_array(np.squeeze(output_array),np.squeeze(mask_header[:,0]),pixels_first=True)

    output_series = input_series.new_sibling(SeriesDescription = input_series.instance().SeriesDescription + ' _filled')
    output_series.set_array(np.squeeze(output_array),np.squeeze(mask_header[:,0]),pixels_first=True)
    output_series.set_affine(mask_series.affine())

    return output_series