"""
Script to apply the pretrained UNETR for 3D kidney segmentation.

You will need:

- dixon: A 3D numpy array with post-contrast out of phase DIXON data
- weights: A file with pretrained model weights

Apply the model:

- segments = UNETR_kidneys_v1.apply(dixon, weights)

Clean the results and extract masks as numpy arrays:

- left_kidney, right_kidney = UNETR_kidneys_v1.kidney_masks(segments)

or do both in one go:

- left_kidney, right_kidney = UNETR_kidneys_v1.kidneys(dixon, weights)
"""

import os
import numpy as np 
import scipy.ndimage as ndi
import torch
from monai.inferers import sliding_window_inference
from monai.networks.nets import UNETR

# The default model name must be the same as the filename of this module
filename = os.path.basename(os.path.abspath(__file__))
filename = filename[:-3] + ".pth"

# Setup device
os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
device_str = "cuda" if torch.cuda.is_available() else "cpu"
device = torch.device(device_str)

# Define model architecture
model = UNETR(
    in_channels=1,
    out_channels=3, # BACKGROUND, RIGHT KIDNEY (left on image), LEFT KIDNEY (right on image)
    img_size=(80, 80, 80),
    feature_size=16,
    hidden_size=768,
    mlp_dim=3072,
    num_heads=12,
    pos_embed="perceptron",
    norm_name="instance",
    res_block=True,
    dropout_rate=0.0,
).to(device)

# Required - DICOM series description of validated data
trained_on = "Dixon_post_contrast_out_phase"


def largest_cluster(array:np.ndarray)->np.ndarray:
    """Given a mask array, return a new mask array containing only the largesr cluster.

    Args:
        array (np.ndarray): mask array with values 1 (inside) or 0 (outside)

    Returns:
        np.ndarray: mask array with only a single connect cluster of pixels.
    """
    # Label all features in the array
    label_img, cnt = ndi.label(array)
    # Find the label of the largest feature
    labels = range(1,cnt+1)
    size = [np.count_nonzero(label_img==l) for l in labels]
    max_label = labels[size.index(np.amax(size))]
    # Return a mask corresponding to the largest feature
    return label_img==max_label


def kidney_masks(output_array:np.ndarray)->tuple:
    """Extract kidney masks from the output array of the UNETR

    Args:
        output_array (np.ndarray): 3D numpy array (x,y,z) with integer labels (0=background, 1=right kidney, 2=left kidney)

    Returns:
        tuple: A tuple of 3D numpy arrays (left_kidney, right_kidney) with masks for the kidneys.
    """
    left_kidney = largest_cluster(output_array == 2)
    right_kidney = largest_cluster(output_array == 1)
    return left_kidney, right_kidney


# Required
def apply(input_array:np.ndarray, file:str, overlap=0.3)->np.ndarray:
    """apply UNETR model to DIXON out of phase volume.

    Args:
        input_array (np.ndarray): 3D numpy array (x,y,z) with DIXON data.
        file (str): filepath to file with model weights.
        overlap (float, optional): optimization parameter. Defaults to 0.2.

    Returns:
        np.ndarray: 3D numpy array (x,y,z) with integer labels (0=background, 1=right kidney, 2=left kidney)
    """

    # Normalize data
    input_array = (input_array-np.average(input_array))/np.std(input_array)
    
    # Convert to NCHW[D] format: (1,1,y,x,z)
    # NCHW[D] stands for: batch N, channels C, height H, width W, depth D
    input_array = input_array.transpose(1,0,2) # from (x,y,z) to (y,x,z)
    input_array = np.expand_dims(input_array, axis=(0, 1))

    # Convert to tensor
    input_tensor = torch.tensor(input_array)

    # Load model weights
    weights = torch.load(file, map_location=device)
    model.load_state_dict(weights)
    model.eval()

    # Calculate model output (decrease overlap parameter for faster but less accurate results)
    with torch.no_grad():
        output_tensor = sliding_window_inference(input_tensor, (80,80,80), 4, model, overlap=overlap, device=device_str, progress=True) 

    # From probabilities for each channel to label image
    output_tensor = torch.argmax(output_tensor, dim=1)

    # Convert to numpy
    output_array = output_tensor.numpy(force=True)[0,:,:,:]
        
    # Transpose to original shape
    output_array = output_array.transpose(1,0,2) #from (y,x,z) to (x,y,z)

    return output_array


def kidneys(input_array:np.ndarray, file:str, **kwargs)->tuple:
    """return kdiney masks directly from data

    Args:
        input_array (np.ndarray): 3D numpy array (x,y,z) with DIXON data.
        file (str): filepath to file with model weights.
        overlap (float, optional): optimization parameter. Defaults to 0.2.

    Returns:
        tuple: A tuple of 3D numpy arrays (left_kidney, right_kidney) with masks for the kidneys.
    """
    masks = apply(input_array, file, **kwargs)
    return kidney_masks(masks)
