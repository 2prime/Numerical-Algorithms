import numpy as np
import tensorflow as tf
import matplotlib.pyplot as plt
import matplotlib.image  as mpimg
from tensorflow.python.platform import gfile
from scipy import misc
from vm import *

import os
__author__ = "2prime"


os.environ["CUDA_VISIBLE_DEVICES"]="2"

def weight_variable(shape):
  initial = tf.truncated_normal(shape, stddev=0.1)
  return tf.Variable(initial)

def bias_variable(shape):
  initial = tf.constant(0.1, shape=shape)
  return tf.Variable(initial)

def conv2d(x, W):
  return tf.nn.conv2d(x, W, strides=[1, 1, 1, 1], padding='VALID')

class Linear_PDENET(object):
    '''
    Learning linearize PDE no more than 3 order
    '''
    def __init__(self,size):
        self.conv_num = 1
        self.size = size
        # data placeholder for scale1
        self.data_scale1 = tf.placeholder(tf.float32, [None, size[0], size[1],size[2]])
        # data placeholder for scale2
        self.data_scale2 = tf.placeholder(tf.float32, [None, size[0], size[1]//2+1, size[2]//2+1])
        self.data_scale3 = tf.placeholder(tf.float32, [None, size[0], size[1]//4+1, size[2]//4+1])


        # geting the boundary value
        self.start_scale1 = self.data_scale1[:, 0, :, :]
        self.start_scale2 = self.data_scale2[:, 0, :, :]
        self.start_scale3 = self.data_scale3[:, 0, :, :]

        self.start_scale1 = tf.reshape(self.start_scale1,[-1,size[1],size[2],1])
        self.start_scale2 = tf.reshape(self.start_scale2, [-1, size[1]//2+1, size[2]//2+1, 1])
        self.start_scale3 = tf.reshape(self.start_scale3, [-1, size[1]//4+1, size[2]//4+1, 1])

        self.conv_1 = []#0order
        self.bias_1 = []
        self.conv_2 = []#1order
        self.bias_2 = []
        self.conv_3 = []#2order
        self.bias_3 = []

        for i in range(size[0]-1):
            dif1 = []
            dif2 = []
            self.conv_1.append(tf.Variable(tf.random_normal([5, 5, 1, self.conv_num], mean=0,stddev= 0.0001)))
            self.conv_2.append(tf.Variable(tf.random_normal([5, 5, 1, self.conv_num], mean=0,stddev= 0.0001)))
            self.conv_3.append(tf.Variable(tf.random_normal([5, 5, 1, self.conv_num], mean=0, stddev= 0.01)))
            self.bias_1.append(bias_variable([size[1],size[2], self.conv_num]))
            self.bias_2.append(bias_variable([size[1],size[2], self.conv_num]))
            self.bias_3.append(bias_variable([size[1],size[2], self.conv_num]))



            self.dif_1_scale1 = conv2d(tf.pad(self.start_scale1, [[0, 0], [2, 2], [2, 2], [0, 0]], "REFLECT"),
                                       self.conv_1[0])
            self.dif_1_scale2 = conv2d(tf.pad(self.start_scale2, [[0, 0], [2, 2], [2, 2], [0, 0]], "REFLECT"),
                                       self.conv_1[0])
            self.dif_1_scale3 = conv2d(tf.pad(self.start_scale3, [[0, 0], [2, 2], [2, 2], [0, 0]], "REFLECT"),
                                       self.conv_1[0])



            self.dif_2_scale1 = conv2d(tf.pad(self.start_scale1, [[0, 0], [2, 2], [2, 2], [0, 0]], "REFLECT"),
                                       self.conv_2[0])
            self.dif_2_scale2 = conv2d(tf.pad(self.start_scale2, [[0, 0], [2, 2], [2, 2], [0, 0]], "REFLECT"),
                                       2 * self.conv_2[0])
            self.dif_2_scale3 = conv2d(tf.pad(self.start_scale3, [[0, 0], [2, 2], [2, 2], [0, 0]], "REFLECT"),
                                       4 * self.conv_2[0])



            self.dif_3_scale1 = conv2d(tf.pad(self.start_scale1, [[0, 0], [2, 2], [2, 2], [0, 0]], "REFLECT"),
                                       self.conv_3[0])
            self.dif_3_scale2 = conv2d(tf.pad(self.start_scale2, [[0, 0], [2, 2], [2, 2], [0, 0]], "REFLECT"),
                                       4 * self.conv_3[0])
            self.dif_3_scale3 = conv2d(tf.pad(self.start_scale3, [[0, 0], [2, 2], [2, 2], [0, 0]], "REFLECT"),
                                       16 * self.conv_3[0])

            self.start_scale1 = self.start_scale1 + self.dif_1_scale1 + self.dif_2_scale1 + self.dif_3_scale1
            self.start_scale2 = self.start_scale2 + self.dif_1_scale2 + self.dif_2_scale2 + self.dif_3_scale2
            self.start_scale3 = self.start_scale3 + self.dif_1_scale3 + self.dif_2_scale3 + self.dif_3_scale3

            if not i:
                self.err = tf.reduce_mean(tf.square(self.start_scale1-(tf.reshape(self.data_scale1[:, i+1, :, :],[-1,size[1],size[2],1]))))
                self.err = self.err + tf.reduce_mean(tf.square(self.start_scale2 - (tf.reshape(self.data_scale2[:, i + 1, :, :], [-1, size[1]//2+1, size[2]//2+1, 1]))))
                self.err = self.err + tf.reduce_mean(tf.square(self.start_scale3 - (tf.reshape(self.data_scale3[:, i + 1, :, :], [-1, size[1]//4+1, size[2]//4+1, 1]))))
            else:
                self.err = self.err + tf.reduce_mean(tf.square(self.start_scale1 - (tf.reshape(self.data_scale1[:, i + 1, :, :], [-1, size[1], size[2], 1]))))
                self.err = self.err + tf.reduce_mean(tf.square(self.start_scale2 - (tf.reshape(self.data_scale2[:, i + 1, :, :], [-1, size[1] // 2 + 1, size[2] // 2 + 1, 1]))))
                self.err = self.err + tf.reduce_mean(tf.square(self.start_scale3 - (tf.reshape(self.data_scale3[:, i + 1, :, :], [-1, size[1] // 4 + 1, size[2] // 4 + 1, 1]))))

            # vanish momet regularization


            m1 = vanish_moment([5, 5], [0, 0])
            m2 = vanish_moment([5, 5], [1, 0])
            m3 = vanish_moment([5, 5], [0, 1])
            m4 = vanish_moment([5, 5], [1, 1])
            m5 = vanish_moment([5, 5], [2, 0])
            m6 = vanish_moment([5, 5], [0, 2])
            self.err = self.err + tf.abs(tf.reduce_sum(self.conv_1[i][:, :, 0, 0] * m1) * 0.001)

            self.err = self.err + tf.abs(tf.reduce_sum(self.conv_2[i][:, :, 0, 0] * m1) * 0.001)
            self.err = self.err + tf.abs(tf.reduce_sum(self.conv_2[i][:, :, 0, 0] * m2) * 0.001)
            self.err = self.err + tf.abs(tf.reduce_sum(self.conv_2[i][:, :, 0, 0] * m3) * 0.001)

            self.err = self.err + tf.abs(tf.reduce_sum(self.conv_3[i][:, :, 0, 0] * m1) * 0.001)
            self.err = self.err + tf.abs(tf.reduce_sum(self.conv_3[i][:, :, 0, 0] * m2) * 0.001)
            self.err = self.err + tf.abs(tf.reduce_sum(self.conv_3[i][:, :, 0, 0] * m3) * 0.001)
            self.err = self.err + tf.abs(tf.reduce_sum(self.conv_3[i][:, :, 0, 0] * m4) * 0.001)
            self.err = self.err + tf.abs(tf.reduce_sum(self.conv_3[i][:, :, 0, 0] * m5) * 0.001)
            self.err = self.err + tf.abs(tf.reduce_sum(self.conv_3[i][:, :, 0, 0] * m6) * 0.001)



        self.opt1 = tf.train.ProximalAdagradOptimizer(1e-3).minimize(self.err)
        self.opt2 = tf.train.AdamOptimizer(5e-6).minimize(self.err)
        self.opt3 = tf.train.AdamOptimizer(1e-6).minimize(self.err)
        self.opt4 = tf.train.GradientDescentOptimizer(1e-7).minimize(self.err)

        self.sess = tf.Session()
        self.sess.run(tf.global_variables_initializer())
        print("model constructed")
        self.iter_time = 0


    def opt(self,data):
        data_2 = data[:, :, ::2, ::2]
        data_3 = data[:, :, ::4, ::4]

        print("begin opt")
        print("error",self.err.eval(session=self.sess,feed_dict={self.data_scale1:data,self.data_scale2:data_2,self.data_scale3:data_3}))

        if self.iter_time < 20:
            self.opt1.run(session=self.sess,feed_dict={self.data_scale1:data,self.data_scale2:data_2,self.data_scale3:data_3})
        elif self.iter_time < 100:
            self.opt2.run(session=self.sess,feed_dict={self.data_scale1:data,self.data_scale2:data_2,self.data_scale3:data_3})
        elif self.iter_time <300:
            self.opt3.run(session=self.sess,feed_dict={self.data_scale1:data,self.data_scale2:data_2,self.data_scale3:data_3})
        else:
            self.opt4.run(session=self.sess,feed_dict={self.data_scale1:data,self.data_scale2:data_2,self.data_scale3:data_3})

        self.iter_time += 1

    def print_kernel(self):
        for i in range(self.size[0]-1):
            print('conv1:')
            print(self.sess.run(self.conv_1[i]))
            print('conv2:')
            print(self.sess.run(self.conv_2[i]))
            print('conv3:')
            print(self.sess.run(self.conv_3[i]))

    def plot_kernel(self, b_print=False):
        conv1 = []
        conv2 = []
        conv3 = []
        for i in range(1):
            conv1.append(self.sess.run(self.conv_1[i]))
            conv2.append(self.sess.run(self.conv_2[i]))
            conv3.append(self.sess.run(self.conv_3[i]))

            m1 = vanish_moment([5,5],[0,0])
            m2 = vanish_moment([5,5],[1,0])
            m3 = vanish_moment([5,5],[0,1])
            m4 = vanish_moment([5, 5], [1, 1])
            m5 = vanish_moment([5, 5], [2, 0])
            m6 = vanish_moment([5, 5], [0, 2])

            fig, axes= plt.subplots(nrows=3,ncols=1)
            index = 0
            for ax in axes.flat:
                if index == 0:
                    im = ax.imshow(conv1[i][:,:,0,0],cmap="gray")

                    a1 = tf.reduce_sum(conv1[i][:,:,0,0] * m1)
                    v1 = self.sess.run(a1)
                    a2 = tf.reduce_sum(conv1[i][:, :, 0, 0] * m2)
                    v2 = self.sess.run(a2)
                    a3 = tf.reduce_sum(conv1[i][:, :, 0, 0] * m3)
                    v3 = self.sess.run(a3)
                    a4 = tf.reduce_sum(conv1[i][:, :, 0, 0] * m4)
                    v4 = self.sess.run(a4)
                    a5 = tf.reduce_sum(conv1[i][:, :, 0, 0] * m5)
                    v5 = self.sess.run(a5)
                    a6 = tf.reduce_sum(conv1[i][:, :, 0, 0] * m6)
                    v6 = self.sess.run(a6)

                    print("kernel order 0 vanishing moment:%f,%f,%f,%f,%f,%f"%(v1,v2,v3,v4,v5,v6))

                    #ax.set_xlabel("zero order filter")
                if index == 1:
                    im = ax.imshow(conv2[i][:,:,0,0],cmap="gray")

                    a1 = tf.reduce_sum(conv2[i][:, :, 0, 0] * m1)
                    v1 = self.sess.run(a1)
                    a2 = tf.reduce_sum(conv2[i][:, :, 0, 0] * m2)
                    v2 = self.sess.run(a2)
                    a3 = tf.reduce_sum(conv2[i][:, :, 0, 0] * m3)
                    v3 = self.sess.run(a3)
                    a4 = tf.reduce_sum(conv2[i][:, :, 0, 0] * m4)
                    v4 = self.sess.run(a4)
                    a5 = tf.reduce_sum(conv2[i][:, :, 0, 0] * m5)
                    v5 = self.sess.run(a5)
                    a6 = tf.reduce_sum(conv2[i][:, :, 0, 0] * m6)
                    v6 = self.sess.run(a6)

                    print("kernel order 1 vanishing moment:%f,%f,%f,%f,%f,%f" % (v1, v2, v3, v4, v5, v6))

                    #ax.set_xlabel("first order filter")
                if index == 2:
                    im = ax.imshow(conv3[i][:, :, 0, 0], cmap="gray")

                    a1 = tf.reduce_sum(conv3[i][:, :, 0, 0] * m1)
                    v1 = self.sess.run(a1)
                    a2 = tf.reduce_sum(conv3[i][:, :, 0, 0] * m2)
                    v2 = self.sess.run(a2)
                    a3 = tf.reduce_sum(conv3[i][:, :, 0, 0] * m3)
                    v3 = self.sess.run(a3)
                    a4 = tf.reduce_sum(conv3[i][:, :, 0, 0] * m4)
                    v4 = self.sess.run(a4)
                    a5 = tf.reduce_sum(conv3[i][:, :, 0, 0] * m5)
                    v5 = self.sess.run(a5)
                    a6 = tf.reduce_sum(conv3[i][:, :, 0, 0] * m6)
                    v6 = self.sess.run(a6)

                    print("kernel order 2 vanishing moment:%f,%f,%f,%f,%f,%f" % (v1, v2, v3, v4, v5, v6))

                    #ax.set_xlabel("second order filter")
                index += 1


            fig.colorbar(im,ax=axes.ravel().tolist())
            fig.savefig("kernel"+str(i)+".png")









