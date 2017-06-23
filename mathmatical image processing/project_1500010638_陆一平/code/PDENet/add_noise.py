import  numpy as np

def add_noise(np_array,db=0.05):
    shape = np_array.shape
    noise = np.random.rand(shape[0],shape[1],shape[2])
    f_max = np.max(np_array)
    return np_array+noise*f_max*db

if __name__ == "__main__":
    a = np.array([[[1,2],[3,4]]])
    print(add_noise(a))
