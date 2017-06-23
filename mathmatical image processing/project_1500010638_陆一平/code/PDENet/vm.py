import numpy as np


def vanish_moment(shape,dircetion):
    dircetion_matrix = np.array([[(i**dircetion[0])*(j**dircetion[1]) for j in range(shape[1])] for i in range(shape[0])])
    return dircetion_matrix

if __name__ == "__main__":
    print(vanish_moment([5,5],[1,1]))