import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main(arr,path):

    arr = arr.transpose((1,0,2))
    montage = arr.reshape((arr.shape[0], -1), order='F')
    plt.imsave(path + '.png', montage, cmap='viridis')
    np.save(path + '.npy', arr)
