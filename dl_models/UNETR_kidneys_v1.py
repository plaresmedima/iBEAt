import os
import time
import numpy as np 
import torch
from monai.inferers import sliding_window_inference
from monai.networks.nets import UNETR

# The model name must be the same as the filename of this module
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
trained_on = "T1w_abdomen_dixon_cor_bh_out_phase_post_contrast"


# Required
def apply(input_array:np.ndarray, file:str, overlap=0.3)->np.ndarray:
    """apply UNETR model to DIXON out of phase volume.

    Args:
        input_array (np.ndarray): 3D numpy array (x,y,z) with DIXON data
        file (str): filepath to file with model weights
        overlap (float, optional): optimization parameter. Defaults to 0.2.

    Returns:
        np.ndarray: 3D numpy array (x,y,z) with integer labels (0=background, 1=right kidney, 2=left kidney)
    """

    tic = time.time()

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

    toc = time.time()
    print('Calculation time for UNETR prediction: ' + str (toc-tic) + " seconds")

    return output_array
