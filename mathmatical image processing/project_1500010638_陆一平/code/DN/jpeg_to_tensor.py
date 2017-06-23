"""
turn jpegs in a folder to tensors
"""
__author__="2prime"

import sys
import glob
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image  as mpimg
import time
from skimage import transform,data

import warnings

from add_noise import *

warnings.filterwarnings("ignore")

def rgb2gray(rgb):
    '''
    :param rgb: rgb images
    :return: gray images
    '''
    return np.dot(rgb[...,:3], [0.299, 0.587, 0.114])

def is_img(filename):
    return filename[-4:] in [".jpg",".JPG",".png",".bmp"]


def get_data_tensor(data_path,batch,size=(1296, 1728),resize=True ,sigma=25,rgb_factor=True):
    print("turning dataset:",data_path,"into tensor")
    #data_file = filter(is_img, os.listdir(data_path))
    data_tensor = map(lambda x: mpimg.imread(data_path + x), batch)
    print("loaded pictures and turning into tensor")

    #trun into gray pictures
    if rgb_factor:data_tensor = map(rgb2gray, data_tensor)
    if rgb_factor: print("turned into gray picture")


    #resize the pciture
    if resize:
        data_tensor = map(lambda x: transform.resize(x, size), data_tensor)

    tic = time.time()
    data_tensor = np.array(list(data_tensor),dtype=np.float32)
    print("geting data tensor,using time:", time.time() - tic)

    print("noising the picture")
    #ufunc_noise = np.frompyfunc(float_add_noise,1,1)
    #data_noise = ufunc_noise(data_tensor)
    data_noise = sigma*np.random.standard_normal(data_tensor.shape)+data_tensor

    print("end processing")


    return data_tensor,data_noise



if __name__ == "__main__":
    data_path = "./data/MSRA-TD500/train/"
    data_tensor = get_data_tensor(data_path)
    np.save("data_MSRA_TD500_TRAIN.npy",data_tensor)

    for image in data_tensor:
        img = plt.imshow(image)
        img.set_cmap('gray')
        break
    plt.show()
