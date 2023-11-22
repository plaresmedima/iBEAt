import numpy as np
import imageio

def main(array_3d,path,duration):

    # Create frames from 3D array slices
    frames = []
    for slice_2d in array_3d:
        frame = (slice_2d * 255).astype(np.uint8)  # Ensure data type is uint8
        frames.append(frame)

    # Save the frames as a GIF
    imageio.mimsave(path, frames, duration)