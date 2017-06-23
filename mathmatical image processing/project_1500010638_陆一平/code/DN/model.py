import tensorflow as tf
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image  as mpimg
from tensorflow.python.platform import gfile
from scipy import misc

__author__ = "2prime"

conv_number = 96

def weight_variable(shape):
  initial = tf.truncated_normal(shape, stddev=0.1)
  return tf.Variable(initial)

def bias_variable(shape):
  initial = tf.constant(0.1, shape=shape)
  return tf.Variable(initial)

def conv2d(x, W):
  return tf.nn.conv2d(x, W, strides=[1, 1, 1, 1], padding='VALID')

def _tf_fspecial_gauss(size, sigma):
    """Function to mimic the 'fspecial' gaussian MATLAB function
    """
    x_data, y_data = np.mgrid[-size//2 + 1:size//2 + 1, -size//2 + 1:size//2 + 1]

    x_data = np.expand_dims(x_data, axis=-1)
    x_data = np.expand_dims(x_data, axis=-1)

    y_data = np.expand_dims(y_data, axis=-1)
    y_data = np.expand_dims(y_data, axis=-1)

    x = tf.constant(x_data, dtype=tf.float32)
    y = tf.constant(y_data, dtype=tf.float32)

    g = tf.exp(-((x**2 + y**2)/(2.0*sigma**2)))
    return g / tf.reduce_sum(g)

def tf_ssim(img1, img2, cs_map=False, mean_metric=True, size=11, sigma=1.5):
    window = _tf_fspecial_gauss(size, sigma) # window shape [size, size]
    K1 = 0.01
    K2 = 0.03
    L = 1  # depth of image (255 in case the image has a differnt scale)
    C1 = (K1*L)**2
    C2 = (K2*L)**2
    mu1 = tf.nn.conv2d(img1, window, strides=[1,1,1,1], padding='VALID')
    mu2 = tf.nn.conv2d(img2, window, strides=[1,1,1,1],padding='VALID')
    mu1_sq = mu1*mu1
    mu2_sq = mu2*mu2
    mu1_mu2 = mu1*mu2
    sigma1_sq = tf.nn.conv2d(img1*img1, window, strides=[1,1,1,1],padding='VALID') - mu1_sq
    sigma2_sq = tf.nn.conv2d(img2*img2, window, strides=[1,1,1,1],padding='VALID') - mu2_sq
    sigma12 = tf.nn.conv2d(img1*img2, window, strides=[1,1,1,1],padding='VALID') - mu1_mu2
    if cs_map:
        value = (((2*mu1_mu2 + C1)*(2*sigma12 + C2))/((mu1_sq + mu2_sq + C1)*
                    (sigma1_sq + sigma2_sq + C2)),
                (2.0*sigma12 + C2)/(sigma1_sq + sigma2_sq + C2))
    else:
        value = ((2*mu1_mu2 + C1)*(2*sigma12 + C2))/((mu1_sq + mu2_sq + C1)*
                    (sigma1_sq + sigma2_sq + C2))

    if mean_metric:
        value = tf.reduce_mean(1-value)
    return value



class Denoise_Net(object):
    '''
    ResNet for Denoising images
    '''
    def __init__(self,size):
        global conv_number

        self.layer_number= 5
        self.size =size
        self.image_noise = tf.placeholder(tf.float32, [None, size[0],size[1]])
        self.image_noise_1 = tf.reshape(self.image_noise,[-1, size[0],size[1],1])
        self.image = tf.placeholder(tf.float32, [None, size[0],size[1]])

        self.dif_conv1 = []
        self.dif_conv2 = []
        self.syn_conv1 = []
        self.syn_conv2 = []

        for i in range(self.layer_number):
            # W_high:analysis operator
            self.dif_conv1.append(tf.Variable(tf.random_normal([7, 7, 1, conv_number], mean=0,
                                                               stddev=0.01)))  # conv1 operator for the second differential operator
            #self.dif_conv2.append(tf.Variable(tf.random_normal([7, 7, 1, conv_number], mean=0,
                                                               #stddev=10)))  # conv2 operator for the second differential operator

            # W_high^T:synthesis operator
            self.syn_conv1.append(tf.Variable(tf.random_normal([7, 7, conv_number, 1], mean=0, stddev=2)))
            #self.syn_conv2.append(tf.Variable(tf.random_normal([7, 7, conv_number, 1], mean=0, stddev=1)))




        self.b_conv = []
        #self.b_conv2 = bias_variable([size[0],size[1], 1])



        k = tf.constant([[1,-1],[1,-1]],dtype=tf.float32)
        self.canny1 = tf.reshape(k,[2, 2, 1, 1])

        k = tf.constant([[1, 1], [-1, -1]], dtype=tf.float32)
        self.canny2 = tf.reshape(k, [2, 2, 1, 1])

        self.kappa =[]
        self.k1 = []
        self.k2 = []

        self.iter_time = 0

        self.tmp = tf.reshape(self.image, [-1, size[0], size[1], 1])

        #resdual learning
        self.res = self.image_noise_1
        for i in range(self.layer_number):
            #self.k1.append(tf.Variable([1.0]))
            #self.k2.append(tf.Variable([1.0]))
            self.kappa.append(tf.Variable([-0.8]))
            self.b_conv.append(bias_variable([size[0],size[1], conv_number]))

            #shrinkage
            self.dif1 = tf.nn.relu(conv2d(tf.pad(self.res, [[0, 0], [3, 3], [3, 3], [0, 0]], "REFLECT"), self.dif_conv1[i])+self.b_conv[i]) #D_x u


            #self.dif1 = tf.exp(-1*tf.square(self.dif1))

            #self.dif2 = conv2d(self.res, self.dif_conv2[i]) #D_y u

            #self.nabla = tf.square(self.dif1)+tf.square(self.dif2) #|\nabla u|^2

            # g(|\nabla u|^2),g(x)=1/(1+x)
            #self.g_nabla =tf.div(np.array([1.],dtype=np.float32),tf.add(np.array([1.],dtype=np.float32),self.nabla))

            #self.shrinkage1 = self.g_nabla * self.dif1 + self.k1[i]*self.dif1
            #self.shrinkage2 = self.g_nabla * self.dif2 + self.k2[i]*self.dif2

            #synthesis
            self.res_tmp1 = conv2d(tf.pad(self.dif1, [[0, 0], [3, 3], [3, 3], [0, 0]], "REFLECT"),self.syn_conv1[i])
            #self.res_tmp2 = conv2d(self.shrinkage2,self.syn_conv2[i])


            self.res_tmp3 = self.kappa[i] *(self.res - self.image_noise_1)

            self.res = self.res_tmp1 + self.res_tmp3 + self.res
            #greedy learning
            if not i :
                self.err = tf.reduce_mean(tf.square(self.res-self.tmp)) + tf_ssim(self.res,self.tmp)*500
            else:
                self.err = self.err + tf.reduce_mean(tf.square(self.res-self.tmp)) + tf_ssim(self.res,self.tmp)*500




        # l1 soblev error
        #self.err2 = tf.reduce_mean(tf.abs(conv2d(self.res,self.canny1)-conv2d(self.tmp,self.canny1)))
        #self.err3 = tf.reduce_mean(tf.abs(conv2d(self.res, self.canny2) - conv2d(self.tmp, self.canny2)))

        self.err = self.err + tf_ssim(self.res,self.tmp)*500 #+ self.err2 + self.err3


        self.opt_step1 = tf.train.ProximalAdagradOptimizer(1e-3).minimize(self.err)
        self.opt_step2 = tf.train.ProximalAdagradOptimizer(1e-4).minimize(self.err)
        self.opt_step3 = tf.train.AdamOptimizer(1e-5).minimize(self.err)
        self.opt_step4 = tf.train.AdamOptimizer(1e-6).minimize(self.err)

        self.sess = tf.Session()
        self.sess.run(tf.global_variables_initializer())
        print("model constructed")

    def opt(self,input_image,input_noise_image):
        '''
        :param input_image: train image without noise
        :param input_noise_image: train image with noise
        Training the ResNet
        '''
        print("begin opt:,iter time:",self.iter_time)
        self.train_error = self.err.eval(session=self.sess,feed_dict={self.image:input_image,self.image_noise:input_noise_image})
        print("train error:",self.train_error)
        if self.iter_time < 200:
            self.opt_step1.run(session=self.sess,feed_dict={self.image:input_image,self.image_noise:input_noise_image})
        elif self.iter_time <500:
            self.opt_step2.run(session=self.sess,feed_dict={self.image:input_image,self.image_noise:input_noise_image})
        elif self.iter_time <1000:
            self.opt_step3.run(session=self.sess,feed_dict={self.image:input_image,self.image_noise:input_noise_image})
        else:
            self.opt_step4.run(session=self.sess,feed_dict={self.image:input_image,self.image_noise:input_noise_image})
        self.iter_time += 1




    def process(self,input_noise_image,input_image):
        '''
        :param input_noise_image:
        :return: denoised image
        '''
        output = self.res.eval(session=self.sess,feed_dict={self.image_noise:input_noise_image})


        for i in range(len(output)):
            plt.imsave('result_'+str(i)+'_'+str(self.iter_time)+'.jpg',output[i,:,:,0],cmap="gray")
        #plt.imshow(np.concatenate((input_noise_image[0],output[0,:,:,0],input_image),axis=0),cmap='gray')
        #plt.imsave('result'+str(self.iter_time)+'.jpg',np.concatenate((input_noise_image[0],output[0,:,:,0],input_image),axis=0),cmap='gray')
        #plt.show()
        #np.save(np.concatenate((input_noise_image[0],output[0,:,:,0],input_image),axis=0),'result'+str(self.iter_time)+'.npy')
        return output

    def plot_pic(self,image,image_noise,image3):
        plt.subplot(131)
        img1 = plt.imshow(image)
        img1.set_cmap('gray')
        plt.axis('off')
        plt.subplot(132)
        img2 = plt.imshow(image_noise)
        img2.set_cmap('gray')
        plt.axis('off')
        plt.subplot(133)
        img3 = plt.imshow(image3)
        img3.set_cmap('gray')
        plt.axis('off')
        plt.show()


    def plot_kernel(self,index,b_print=False):
        conv1 = []
        conv2 = []
        syn1 = []
        syn2 = []
        for i in range(self.layer_number):
            conv1.append(self.sess.run(self.dif_conv1[i]))
            conv2.append(self.sess.run(self.dif_conv2[i]))
            syn1.append(self.sess.run(self.syn_conv1[i]))
            syn2.append(self.sess.run(self.syn_conv2[i]))

        conv1 = conv1[index]
        conv2 = conv2[index]
        syn1 = syn1[index]
        syn2 = syn2[index]

        for i in range(conv_number):
            plt.subplot(conv_number, 4, i+1)
            plt.imshow(conv1[:, :, 0, i],cmap="gray")

        for i in range(conv_number):
            plt.subplot(conv_number, 4, i+conv_number+1)
            plt.imshow(conv2[:, :, 0, i],cmap="gray")

        for i in range(conv_number):
            if b_print:
                print("conv1 first vanishing moment:",self.sess.run(tf.reduce_mean(self.dif_conv1[index][:, :, 0, i])))
                print("conv2 first vanishing moment:",self.sess.run(tf.reduce_mean(self.dif_conv2[index][:, :, 0, i])))

        for i in range(conv_number):
            plt.subplot(conv_number, 4, i+ 2 * conv_number + 1)
            plt.imshow(syn1[:,:,i,0],cmap="gray")

        for i in range(conv_number):
            plt.subplot(conv_number, 4, i + 3 * conv_number + 1)
            plt.imshow(syn2[:, :, i, 0],cmap="gray")

        for i in range(conv_number):
            if b_print:
                print("syn1 first vanishing moment:", self.sess.run(tf.reduce_mean(self.syn_conv1[index][:, :, i, 0])))
                print("syn2 first vanishing moment:", self.sess.run(tf.reduce_mean(self.syn_conv2[index][:, :, i, 0])))

        plt.savefig('kernel'+ str(self.iter_time) + '_' + str(index) + '.png')

    def save(self):

        conv1 = []
        syn1 = []

        for i in range(self.layer_number):
            conv1.append(self.sess.run(self.dif_conv1[i]))
            syn1.append(self.sess.run(self.syn_conv1[i]))

        np.save('conv1.npy',np.array(conv1))
        np.save('syn1.npy', np.array(syn1))

        #tf.train.write_graph(self.sess.graph_def, "/tmp/load", False)

        #saver = tf.train.Saver(tf.global_variables())


    def load(self):
        self.iter_time = 100
        conv1 = np.load('conv1.npy')
        conv2 = np.load('conv2.npy')
        syn1 = np.load('syn1.npy')
        syn2 = np.load('syn2.npy')
        for i in range(self.layer_number):
            self.sess.run(tf.assign(self.dif_conv1[0], conv1[0][i]))
            self.sess.run(tf.assign(self.dif_conv2[0], conv2[0][i]))
            self.sess.run(tf.assign(self.syn_conv1[0], syn1[0][i]))
            self.sess.run(tf.assign(self.syn_conv2[0], syn2[0][i]))


        kappa = np.load('kappa.npy')
        k1 = np.load('k1.npy')
        k2 = np.load('k2.npy')
        for i in range(self.layer_number):
            self.sess.run(tf.assign(self.kappa[i],kappa[i]))
            self.sess.run(tf.assign(self.k1[i], k1[i]))
            self.sess.run(tf.assign(self.k2[i], k2[i]))



    def test(self,image):
        self.test_ = conv2d(tf.reshape(self.image,[-1, self.size[0],self.size[1],1]),self.canny1)
        image_ = self.sess.run(self.test_,feed_dict={self.image:image})
        plt.imshow(image_[0,:,:,0],cmap='gray')
        plt.show()
        return image_
