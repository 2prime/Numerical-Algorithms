"""
V循环多重网格方法
=======================================
细网格
---------------------------------------
grid1------0------------0--------------
grid2--------0--------0----------------
grid3----------0----0------------------
grid4-------------0--------------------
---------------------------------------
粗网格
=======================================
"""
import time

import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np

from poisson import *
from double_grid import *
import pol


debug = False
error_cal = True

times1 = 8
times2 = 8



t = time.time()


small_grids = True
grid = poisson(lambda x, y: 2 * pi * pi * sin(pi * x) * sin(pi * y), 64 if small_grids else 128)
'''
#64
for i in range(times1):
    grid.Gauss()
    if debug:
        print("complete")

#128
grid = Lifting_operator(grid)
for i in range(times1):
    grid.Gauss()
    if debug:
        print("complete")
'''
#256
grid = Lifting_operator(grid)
for i in range(times1):
    grid.Gauss()
    if debug:
        print("complete")

for times in range(2):
    #128
    grid_1 = Restriction_operator(grid)
    for i in range(times1):
        grid_1.Gauss()
        if debug:
            print("complete")

    #64
    grid_2 = Restriction_operator(grid_1)
    for i in range(times1):
        grid_2.Gauss()
        if debug:
            print("complete")

    #32
    grid_3 = Restriction_operator(grid_2)
    for i in range(times1):
        grid_3.Gauss()
        if debug:
            print("complete")

    #16
    grid_4 = Restriction_operator(grid_3)
    for i in range(times1):
        grid_4.Gauss()
        if debug:
            print("complete")

    #8
    grid_5 = Restriction_operator(grid_4)
    if small_grids:
        grid_5.getanswer()
    else:
        for i in range(2):
            grid_5.Gauss()
            if debug:
                print("complete")

        #4
        grid_6 = Restriction_operator(grid_5)
        grid_6.getanswer()
        if debug:
            print("Complete")


        grid_5.update_error(Lifting_operator(grid_6))
        for i in range(2):
            grid_5.Gauss()
            if debug:
                print("Complete")


    grid_4.update_error(Lifting_operator(grid_5))
    for i in range(times2):
        grid_4.Gauss()
        if debug:
            print("Complete")

    grid_3.update_error(Lifting_operator(grid_4))
    for i in range(times2):
        grid_3.Gauss()
        if debug:
            print("Complete")

    grid_2.update_error(Lifting_operator(grid_3))
    for i in range(times2):
        grid_2.Gauss()
        if debug:
            print("Complete")


    grid_1.update_error(Lifting_operator(grid_2))
    for i in range(times2):
        grid_1.Gauss()
        if debug:
            print("Complete")

    grid.update_error(Lifting_operator(grid_1))
    for i in range(times2):
        grid.Gauss()
        if debug:
            print("Complete")


error_rlt = lambda x, y:  1-(grid.get_value(x, y) / (sin(pi * x) * sin(pi * y)))
error = lambda x, y:  (grid.get_value(x, y) - (sin(pi * x) * sin(pi * y)))
standard = lambda x,y:(sin(pi * x) * sin(pi * y))
if error_cal:
    n = 20
    def error_y(y):
        pol1 = pol.Chebyshev_interpolation(lambda x:error(x,y)**2, -1, 1, n)
        pol2 = pol1.inter()
        return pol2.eval(1)-pol2.eval(-1)
    pol1 = pol.Chebyshev_interpolation(error_y, -1, 1, n)
    pol2 = pol1.inter()
    print("所求误差为%f"%((pol2.eval(1) - pol2.eval(-1))**(1/2)))

    def std_y(y):
        pol1 = pol.Chebyshev_interpolation(lambda x:(sin(pi * x) * sin(pi * y))**2, -1, 1, n)
        pol2 = pol1.inter()
        return pol2.eval(1)-pol2.eval(-1)
    pol1 = pol.Chebyshev_interpolation(std_y, -1, 1, n)
    pol2 = pol1.inter()
    print("标准值积分%f"%((pol2.eval(1) - pol2.eval(-1))**(1/2)))

else:
    for i in range(300):
        x = random() * 2 - 1
        y = random() * 2 - 1
        print("测试点:(%10f,%10f),|| 相对误差:%10f,|| 准确值:%10f,数值解:%10f" % (x, y, error_rlt(x, y),standard(x,y),grid.get_value(x, y)))

print("用时:%f"%(time.time()-t))


fig = plt.figure()
ax = fig.gca(projection='3d')
X = np.arange(-1, 1, 0.025)
Y = np.arange(-1, 1, 0.025)
X, Y = np.meshgrid(X, Y)

R = np.array([[grid.get_value(X[i][j],Y[i][j])for j in range(80)]for i in range(80)])
#R = np.array([[10*error_rlt(X[i][j],Y[i][j])for j in range(80)]for i in range(80)])
#R = np.sin(R)
surf = ax.plot_surface(X, Y, R, rstride=1, cstride=1, cmap=cm.jet,
        linewidth=0, antialiased=False)
ax.set_zlim(-1.01, 1.01)

ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()
