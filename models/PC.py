# phase contrat analysis iBEAt study:right and left PC arteries: Exeter site analysis Siemens
# Author: Ning Wang @2024
# Disclaimer: For internal research purpose only and not to be distributed outside of the iBEAt project analysis

import os
import numpy as np
from skimage.feature import canny
from skimage.filters import sato
from skimage import morphology
from scipy import ndimage as ndi
from skimage.measure import regionprops
from utilities.improc import region_grow_range

def normalize(arr):
    arr_min = arr.min()
    arr_max = arr.max()
    normalized_arr = (arr - arr_min) / (arr_max - arr_min)
    return normalized_arr


def screen(arr):
    height, width = arr.shape[:2]
    center_row = height // 2
    center_col = width // 2
    center = (center_row, center_col)
    radius = 50
    y, x = np.ogrid[:arr.shape[0], :arr.shape[1]]
    current_slice = arr[:, :]
    mask = ((x - center[1]) ** 2 + (y - center[0]) ** 2 <= radius ** 2)
    current_slice_circular = current_slice.copy()
    current_slice_circular[~mask] = 0
    arr[:, :] = current_slice_circular

    return arr


def transpose_data(arr):
    arr_transpose = np.transpose(arr, (2, 1, 0))
    arr_mean = np.mean(arr_transpose, axis=0)
    arr_max = np.max(arr_transpose, axis=0)
    return arr_transpose, arr_mean, arr_max


def run_sato(arr, sigmas):
    black_ridges = False
    mode = 'reflect'
    cval = 0
    filtered_image = sato(arr, sigmas=sigmas, black_ridges=black_ridges, mode=mode, cval=cval)
    return filtered_image


def run_canny(arr, sigma):
    low_threshold = 0.1
    high_threshold = 0.2
    detected_image = canny(arr, sigma=sigma, low_threshold=low_threshold, high_threshold=high_threshold)
    return detected_image


def filter_regions(arr_label):
    regions = regionprops(arr_label)
    for region in regions:
        minr, minc, maxr, maxc = region.bbox
        width = maxc - minc
        height = maxr - minr
        aspect_ratio = width / height
        if aspect_ratio < 0.5:
            arr_label[arr_label == region.label] = 0
    return arr_label


def fill_holes_and_label(arr_mean_canny):
    filled_holes = ndi.binary_fill_holes(arr_mean_canny)
    cleaned_holes = morphology.remove_small_objects(filled_holes, min_size=30)
    labels, num_labels = ndi.label(cleaned_holes)
    return labels, num_labels


def process_left_phase_masks(arr_mean, holes_label_new):
    mask_left = arr_mean * holes_label_new
    left_label, left_number = ndi.label(mask_left)
    mean_values = np.array([
        np.mean(mask_left[left_label == label])
        for label in range(1, left_number + 1)
    ])

    if np.all(np.isnan(mean_values)):
        mask_left = np.zeros_like(mask_left)
    else:
        max_mean_label = np.nanargmax(mean_values) + 1
        mask_left[left_label != max_mean_label] = 0
        mask_left[left_label == max_mean_label] = 1

    return mask_left


def process_right_phase_masks(arr_mean, holes_label_new):
    mask_right = arr_mean * holes_label_new
    label_right, right_number = ndi.label(mask_right)
    mean_values = np.array([
        np.mean(mask_right[label_right == label])
        for label in range(1, right_number + 1)
    ])

    if np.all(np.isnan(mean_values)):
        mask_right = np.zeros_like(mask_right)
    else:
        min_mean_label = np.nanargmin(mean_values) + 1
        mask_right[label_right != min_mean_label] = 0
        mask_right[label_right == min_mean_label] = 1

    return mask_right


def process_masks_by_sato(arr_mean, holes_label_new):
    sato_mask = arr_mean * holes_label_new
    sato_label, sato_number = ndi.label(sato_mask)
    mean_values = np.array([
        np.mean(sato_mask[sato_label == label])
        for label in range(1, sato_number + 1)
    ])

    if np.all(np.isnan(mean_values)):
        sato_mask = np.zeros_like(sato_mask)
    else:
        max_mean_label = np.nanargmax(mean_values) + 1
        sato_mask[sato_label != max_mean_label] = 0
        sato_mask[sato_label == max_mean_label] = 1

    return sato_mask


def find_center_of_mass(mask):
    labeled_mask, num_labels = ndi.label(mask)
    centers = []
    for label in range(1, num_labels + 1):
        labeled_region = labeled_mask == label
        center = ndi.center_of_mass(labeled_region)
        centers.append(center)
    return centers


def find_closest_mask(phase_mask_velocity, sato_mask_velocity):
    # Find center of mass for phase_mask_velocity and sato_mask_velocity
    phase_centers = find_center_of_mass(phase_mask_velocity)
    sato_centers = find_center_of_mass(sato_mask_velocity)

    # Calculate center of whole image
    image_center = np.array(phase_mask_velocity.shape) / 2

    # Calculate distances from image center to centers of masks
    phase_distances = [np.linalg.norm(np.array(center) - image_center) for center in phase_centers]
    sato_distances = [np.linalg.norm(np.array(center) - image_center) for center in sato_centers]

    # Find which center is closer to the center of the whole image
    closest_phase_center = phase_centers[np.argmin(phase_distances)]
    closest_sato_center = sato_centers[np.argmin(sato_distances)]

    # Determine which mask is closer to the image center
    closest_mask = phase_mask_velocity if np.linalg.norm(np.array(closest_phase_center) - image_center) <= np.linalg.norm(np.array(closest_sato_center) - image_center) else sato_mask_velocity

    return closest_mask


def get_curve(arr, mask, pixel_spacing=0.607, venc=120):
    # arr = phase array
    blood_velocity_all = []
    blood_flow_rate_all = []

    for slice_index in range(arr.shape[2]):
        arr_phase = arr
        phase_time_mask = arr_phase[:, :, slice_index] * mask
        blood_velocity = []
        blood_flow_rate = []
        PC_labels, label_counts = np.unique(mask, return_counts=True)
        for label, pixel_number in zip(np.unique(mask), label_counts):
            if label != 0:
                pixel_values = phase_time_mask[mask == label]
                mean_pixel_value = np.mean(pixel_values)
                velocity_mean = abs((mean_pixel_value / 4096) * venc)  # velocity encode for Siemens cm/sec
                blood_velocity.append(velocity_mean)
                lumen_area = pixel_number * pixel_spacing**2 # mm2
                rbf = abs(velocity_mean * (lumen_area / 100) * 60)  # in ml/min
                blood_flow_rate.append(rbf)
        blood_velocity_all.append(blood_velocity)
        blood_flow_rate_all.append(blood_flow_rate)
        vel = np.concatenate(blood_velocity_all)
        flow = np.concatenate(blood_flow_rate_all)
    return vel, flow


def create_mask_left(arr_velocity, arr_phase, sigmas_sato=range(1, 6, 1), sigma_canny=3.0):

    # (396, 576, 20)

    arr_velocity_20_576_400, _, _ = transpose_data(arr_velocity)
    _, arr_phase_20_576_400_mean, _ = transpose_data(arr_phase)

    # (20, 576, 396), (576, 396)

    arr_velocity_576_20_400 = np.transpose(arr_velocity_20_576_400, (1, 2, 0)) 
    
    # (576, 396, 20)

    arr_velocity_sato = run_sato(arr_velocity_576_20_400, sigmas_sato)
    arr_velocity_sato_20_576_400 = np.transpose(arr_velocity_sato, (2, 0, 1))

    # (20, 576, 396)

    arr_sato_velocity_mean = np.mean(arr_velocity_sato_20_576_400, axis=0)

    # (576, 396)

    arr_sato_velocity_mean_screen = screen(arr_sato_velocity_mean)

    arr_sato_velocity_mean_screen_nor = normalize(arr_sato_velocity_mean_screen)
    arr_sato_velocity_mean_canny = run_canny(arr_sato_velocity_mean_screen_nor, sigma_canny)
    arr_sato_velocity_holes_label, arr_sato_velocity_holes_label_number = fill_holes_and_label(arr_sato_velocity_mean_canny)
    arr_sato_velocity_holes_label = filter_regions(arr_sato_velocity_holes_label)
    holes_label_velocity_new = arr_sato_velocity_holes_label
    holes_label_velocity_new, label_velocity_number_new = ndi.label(holes_label_velocity_new)
    phase_mask_velocity = process_left_phase_masks(arr_phase_20_576_400_mean, holes_label_velocity_new)
    sato_mask_velocity = process_masks_by_sato(arr_sato_velocity_mean, holes_label_velocity_new)
    mask = find_closest_mask(phase_mask_velocity, sato_mask_velocity)

    return mask.T


def create_mask_right(arr_velocity, arr_phase, sigmas_sato=range(1, 6, 1), sigma_canny=3.0):

    # (396, 576, 20)

    arr_velocity_20_576_400, _, _ = transpose_data(arr_velocity)
    _, arr_phase_20_576_400_mean, _ = transpose_data(arr_phase)

    # (20, 576, 396), (576, 396)

    arr_velocity_576_20_400 = np.transpose(arr_velocity_20_576_400, (1, 2, 0)) 
    
    # (576, 396, 20)

    arr_velocity_sato = run_sato(arr_velocity_576_20_400, sigmas_sato)
    arr_velocity_sato_20_576_400 = np.transpose(arr_velocity_sato, (2, 0, 1))

    # (20, 576, 396)

    arr_sato_velocity_mean = np.mean(arr_velocity_sato_20_576_400, axis=0)

    # (576, 396)

    arr_sato_velocity_mean_screen = screen(arr_sato_velocity_mean)

    arr_sato_velocity_mean_screen_nor = normalize(arr_sato_velocity_mean_screen)
    arr_sato_velocity_mean_canny = run_canny(arr_sato_velocity_mean_screen_nor, sigma_canny)
    arr_sato_velocity_holes_label, arr_sato_velocity_holes_label_number = fill_holes_and_label(arr_sato_velocity_mean_canny)
    arr_sato_velocity_holes_label = filter_regions(arr_sato_velocity_holes_label)
    holes_label_velocity_new = arr_sato_velocity_holes_label
    holes_label_velocity_new, label_velocity_number_new = ndi.label(holes_label_velocity_new)
    phase_mask_velocity = process_right_phase_masks(arr_phase_20_576_400_mean, holes_label_velocity_new)
    sato_mask_velocity = process_masks_by_sato(arr_sato_velocity_mean, holes_label_velocity_new)
    mask = find_closest_mask(phase_mask_velocity, sato_mask_velocity)

    return mask.T




def params(time, vel, flow):
    v_mean = np.mean(vel)
    psv = np.amax(vel) # peak systolic velocity
    edv = np.amin(vel) # end diastolic velocity
    ri = (psv-edv)/psv # resistive index
    pi = (psv-edv)/v_mean # pulsatility index
    f_mean = np.mean(flow)
    psf = np.amax(flow) # peak systolic flow
    edf = np.amin(flow) # end diastolic flow
    pif = (psf-edf)/f_mean
    return (
        v_mean,
        psv,
        edv,
        ri,
        pi,
        f_mean,
        psf,
        edf,
        pif,
    )

def pars():
    return (
        'Mean velocity',
        'Peak systolic velocity',
        'End diastolic velocity',
        'Resistive index',
        'Pulsatility index',
        'Mean blood flow',
        'Peak systolic flow',
        'End diastolic flow',
        'Flow pulsatility index',
    )

def units():
    return (
        'cm/sec',
        'cm/sec',
        'cm/sec',
        '',
        '',
        'mL/min',
        'mL/min',
        'mL/min',
        '',
    )


def curve(vel, mask, pixel_spacing=0.607):
    nt = vel.shape[-1]
    velocity = np.zeros(nt)
    flow = np.zeros(nt)
    for t in range (nt):
        vt = vel[:,:,t][mask>0.5]
        velocity[t] = np.amax(vt) # cm/sec
        flow[t] = np.sum(vt) * pixel_spacing**2 * 60/100 # mL/min
    return velocity, flow



def renal_artery_mask(magn, vel, pixel_spacing=0.607, vel_min=30):
    
    width = 5 # cm

    # Find maximum over the cardiac cycle
    magn = np.amax(magn, axis=-1)
    vel = np.amax(vel, axis=-1)
    
    # Find the pixel with maximum magnitude 
    # in a square with given width around the center
    xc, yc = int(magn.shape[0]/2), int(magn.shape[1]/2)
    width = int(width*10/pixel_spacing)
    magn_max = -1
    for x in range(xc-width, xc+width):
        for y in range(yc-width, yc+width):
            if magn[x,y] > magn_max:
                magn_max = magn[x,y]
                p_max = [x,y]

    # Grow region from pixel with maximum magnitude
    # select pixels with max velocity > vel_min cm/sec
    # return mask array
    #return region_grow_range(vel, [p_max], np.mean(vel[p_max[0]-5:p_max[0]+5, p_max[1]-5:p_max[1]+5]), 10*vel[p_max[0], p_max[1]])
    return region_grow_range(vel, [p_max], vel_min, 10*vel[p_max[0], p_max[1]])









