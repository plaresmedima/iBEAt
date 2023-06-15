import os
import time
import numpy as np 
import torch
from monai.inferers import sliding_window_inference
from monai.networks.nets import UNETR
from monai.data import MetaTensor

root_dir = os.path.dirname(os.path.abspath(__file__))
filename = os.path.basename(os.path.abspath(__file__))
modelname = filename[:-3]
weights = os.path.join(root_dir, modelname + ".pth")

# Required 
trained_on = "T1w_abdomen_dixon_cor_bh_out_phase_post_contrast"

# Required
def apply(input_array):
    tic = time.time()
    os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
    device_str = "cuda" if torch.cuda.is_available() else "cpu"
    device = torch.device(device_str)

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

    # Normalizing numpy array
    input_array = input_array.transpose(4,3,1,0,2) # from (x,y,z,1,1) to (1,1,y,x,z)
    input_array = (input_array-np.average(input_array))/np.std(input_array)

    # converting numpy array to monai tensor...
    input_tensor = torch.tensor(input_array)
    meta = {"some": "info"}
    input_tensor = MetaTensor(input_tensor, meta=meta)

    # Loading DL model.
    model.load_state_dict(torch.load(weights, map_location=device))
    model.eval()

    with torch.no_grad():
        
        # Calculating model output
        # Note: decrease overlap parameter for faster (less accurate) results
        output_tensor = sliding_window_inference(input_tensor, (80,80,80), 4, model, overlap=0.3, device=device_str) 
    
        # Converting monai tensor to numpy array.
        output_array = np.zeros(output_tensor[0,0,:,:,:].size(),dtype=int)
        output_array[:,:,:] = torch.argmax(output_tensor, dim=1).detach().cpu()[0,:,:,:]
        
    toc = time.time()
    print('Calculation time for UNETR prediction: ' + str (toc-tic) + " seconds")

    output_array = output_array.transpose(1,0,2) #from (y,x,z) to (x,y,z)

    return output_array
