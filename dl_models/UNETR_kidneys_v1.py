import os
import time
import numpy as np 
import torch
from monai.inferers import sliding_window_inference
from monai.networks.nets import UNETR
from monai.data import MetaTensor

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
def weights_file(weights_dir:str=None):
    # If a folder is not provided, use default folder
    if weights_dir is None:
        weights_dir = os.path.dirname(os.path.abspath(__file__))
    # Check if file with model weights is valid
    weights_file = os.path.join(weights_dir, filename)
    if not os.path.isfile(weights_file):
        msg = 'The file ' + weights_file + ' has not been found. \n'
        msg += 'Please check that the file with model weights is in the folder, and is named ' + filename
        raise FileNotFoundError(msg)
    return weights_file

# Required
def apply(input_array:np.ndarray, file:str=None, overlap=0.2):

    tic = time.time()
    
    # Load model weights
    if file is None:
        file = weights_file()
    weights = torch.load(file, map_location=device)
    model.load_state_dict(weights)
    model.eval()

    # Normalize numpy array
    input_array = input_array.transpose(4,3,1,0,2) # from (x,y,z,1,1) to (1,1,y,x,z)
    input_array = (input_array-np.average(input_array))/np.std(input_array)

    # Convert numpy array to monai tensor
    input_tensor = torch.tensor(input_array)
    meta = {"some": "info"}
    input_tensor = MetaTensor(input_tensor, meta=meta)

    with torch.no_grad():
        
        # Calculating model output (decrease overlap parameter for faster but less accurate results)
        output_tensor = sliding_window_inference(input_tensor, (80,80,80), 4, model, overlap=overlap, device=device_str) 
    
        # Converting monai tensor to numpy array.
        output_array = np.zeros(output_tensor[0,0,:,:,:].size(),dtype=int)
        output_array[:,:,:] = torch.argmax(output_tensor, dim=1).detach().cpu()[0,:,:,:]
        
    output_array = output_array.transpose(1,0,2) #from (y,x,z) to (x,y,z)

    toc = time.time()
    print('Calculation time for UNETR prediction: ' + str (toc-tic) + " seconds")

    return output_array
