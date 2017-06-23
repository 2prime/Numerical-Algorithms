from PDENet import *
from math import exp,pi
from random import random
import numpy as np
import os
from add_noise import *


os.environ["CUDA_VISIBLE_DEVICES"]="4"

#Linear_PDENET([2,2,3]).print()

data_size = 200
batch = 20
size = 256
position = 0

seed = [[random()-0.5,random()-0.5,random()-0.5] for i in range(data_size)]

#generate data

a= [[[[exp(-((x+seed[i][0])**2+(y+seed[i][1])**2)/((size**2)*4*((0.005*t+seed[i][2])**2)))/(0.005*t++seed[i][2]) for x in range(size+1)] for y in range(size+1)] for t in range(1,10)] for i in range(data_size)]
a= list(map(lambda x:add_noise(np.array(x)),a))

a=np.array(a)

#data = np.array([])
print(a.shape)


Net = Linear_PDENET(a.shape[1:])
Net.plot_kernel()

for i in range(100000):
    if (i) and (not (i%100)):
        Net.plot_kernel()
    Net.opt(a[position%100:min(position%data_size+batch,data_size)])
    position += batch

Net.print_kernel()