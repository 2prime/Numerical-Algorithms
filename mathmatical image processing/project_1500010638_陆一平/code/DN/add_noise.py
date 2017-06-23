import numpy as np
import matplotlib
import random
import matplotlib.pyplot as plt
import matplotlib.image  as mpimg

def rgb2gray(rgb):
    '''
    :param rgb: rgb images
    :return: gray images
    '''
    return np.dot(rgb[...,:3], [0.299, 0.587, 0.114])

def float_add_noise(x,sigma=25):
    return x+random.normalvariate(0,sigma)

def add_noise(picture,sigma=25):
    '''
    :param picture: input picture
    :param sigma: var
    :return: noised image
    '''
    size = picture.shape
    return sigma*np.random.standard_normal(size)+picture

if __name__ == "__main__":
    data_path = "./data/MSRA-TD500/train/IMG_0063.JPG"
    image = mpimg.imread(data_path)
    image = rgb2gray(image)
    plt.subplot(121)
    img1 = plt.imshow(image)
    img1.set_cmap('gray')
    plt.axis('off')
    plt.subplot(122)
    img2 = plt.imshow(add_noise(image))
    img2.set_cmap('gray')
    plt.axis('off')

    plt.show()

