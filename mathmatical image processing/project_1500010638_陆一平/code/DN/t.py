import tensorflow as tf
import numpy as np
import matplotlib.pyplot as plt
import random
from tensorflow.python.platform import gfile

from jpeg_to_tensor import *
from model import *


import os


os.environ["CUDA_VISIBLE_DEVICES"]="1"

data_path = "./data/train/"
data_path_test = "./data/test/"

print(os.listdir(data_path))

dir_train = list(filter(is_img, os.listdir(data_path)))
dir_test = list(filter(is_img, os.listdir(data_path_test)))

print(dir_test)


#size=(1296, 1728)
size=(256,256)


#get batch
position_now = 0
batch_size = 20
def get_batch(batch_size,data_size,data1):
    global position_now
    position_now = (position_now + batch_size) % data_size
    print(len(range(position_now,min(position_now+batch_size,data_size))))
    return data1[position_now:min(position_now+batch_size,data_size)]
# get data


data_num = len(dir_train)

data_tensor, data_noise = get_data_tensor(data_path, get_batch(batch_size, data_num, dir_train), size)


#test

test_tensor,test_noise = get_data_tensor(data_path_test,dir_test,size)
print(len(test_tensor))
for i in range(len(test_tensor)):
    plt.imsave('original'+str(i)+'.jpg',test_tensor[i],cmap="gray")
    plt.imsave('test'+str(i)+'.jpg',test_noise[i],cmap="gray")
print('save')

#construct Network
Diffusion_Net = Denoise_Net(data_tensor.shape[-2:])
Diffusion_Net.save()

#load trained Network
#Diffusion_Net.load()

Diffusion_Net.process(test_noise,test_tensor[0])

for i in range(200001):
    if not i%50 and i>400:
        #save trained Network
        #test_tensor, test_noise = get_data_tensor(data_path_test, [dir_test[index]], size)
        Diffusion_Net.process(test_noise, test_tensor[0])
        Diffusion_Net.save()
        if not i%20000:
            for i in range(Diffusion_Net.layer_number):
                Diffusion_Net.plot_kernel(i)
        #index = random.randint(0,test_noise.shape[0]-1)
        #Diffusion_Net.process([test_noise[index]], test_tensor[index])

    data_tensor, data_noise = get_data_tensor(data_path, get_batch(batch_size, data_num, dir_train), size)
    Diffusion_Net.opt(data_tensor,data_noise)

for i in range(Diffusion_Net.layer_number):
    Diffusion_Net.plot_kernel(i,b_print=True)


