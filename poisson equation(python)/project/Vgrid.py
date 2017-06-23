import time

if __name__ == "__main__":
    import numpy as np
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm
    from matplotlib.ticker import LinearLocator, FormatStrFormatter
    import matplotlib.pyplot as plt


from poisson import *
from double_grid import *
import pol



debug = True
plot_show = True



times1 = 2
times2 = 2

V_time = 10

largest_grid = 9
restriciton_time = 3
V_times = 3
grid = []
grid.append(poisson(lambda x, y: 2 * pi * pi * sin(pi * x) * sin(pi * y), 2**largest_grid))


def V_cyclic_multigrid(grid,times1,times2,restriciton_time=3,V_times=1):
    """实现一次V循环"""
    for i in range(V_times):
        # 限制过程
        for i_1 in range(largest_grid-restriciton_time):
            for j in range(times1):
                grid[-1].Gauss()
                if grid[-1].N >= 512 and debug:
                    print("%5d阶网格Gauss变换完成%d次" % (grid[-1].N, j + 1))
            grid.append(Restriction_operator(grid[-1]))
            for j in range(times1):
                grid[-1].Gauss()
                if grid[-1].N >=512 and debug:
                    print("%5d阶网格Gauss变换完成%d次"%(grid[-1].N,j+1))
            if debug:
                print("限制过程完成%d次"%(i_1+1))

        # 计算最粗网络的标准值
        last_grid = Restriction_operator(grid[-1])
        if debug:
            print("V字顶点的剖分次数:%d"%(last_grid.N))
            print("到达V字顶点")
        if last_grid.N <= 4:
            last_grid.getanswer()
        else:
            for j in range(5):
                last_grid.Gauss()
        if debug:
            print("V字精确值计算结束")
        grid[-1].update_error(Lifting_operator(last_grid))

        # 提升过程
        for i_1 in range(largest_grid-restriciton_time):
            grid_new = grid.pop()
            for j in range(times2):
                grid_new.Gauss()
                if grid[-1].N >= 512 and debug:
                    print("%5d阶网格Gauss变换完成%d次" % (grid[-1].N, j + 1))
            grid[-1].update_error(Lifting_operator(grid_new))
            if debug:
                print("提升过程完成%d次"%(i_1+1))
        if debug:
            print("V循环完成%d次"%(i+1))

def getError(grid,n=20,func=lambda x,y:(sin(pi * x) * sin(pi * y))):
    # 计算相对误差
    error_rlt = lambda x, y:  1-(grid[0].get_value(x, y) / func(x,y))
    # 计算绝对误差
    error = lambda x, y:  (grid[0].get_value(x, y) - func(x,y))

    # 计算n价Gauss积分
    def error_y(y):
        pol1 = pol.Chebyshev_interpolation(lambda x: error(x, y) ** 2, -1, 1, n)
        pol2 = pol1.inter()
        return pol2.eval(1) - pol2.eval(-1)

    pol1 = pol.Chebyshev_interpolation(error_y, -1, 1, n)
    pol2 = pol1.inter()
    print("所求误差为%f" % ((pol2.eval(1) - pol2.eval(-1)) ** (1 / 2)))

    def std_y(y):
        pol1 = pol.Chebyshev_interpolation(lambda x: (sin(pi * x) * sin(pi * y)) ** 2, -1, 1, n)
        pol2 = pol1.inter()
        return pol2.eval(1) - pol2.eval(-1)

    pol1 = pol.Chebyshev_interpolation(std_y, -1, 1, n)
    pol2 = pol1.inter()
    print("标准值积分%f" % ((pol2.eval(1) - pol2.eval(-1)) ** (1 / 2)))

if __name__ == "__main__":

    V_cyclic_multigrid(grid,times1,times2,3,V_times)

    getError(grid)


    # 计算相对误差
    error_rlt = lambda x, y:  1-(grid[0].get_value(x, y) / (sin(pi * x) * sin(pi * y)))
    # 计算绝对误差
    error = lambda x, y:  (grid[0].get_value(x, y) - (sin(pi * x) * sin(pi * y)))
    # 计算标准函数值
    standard = lambda x,y:(sin(pi * x) * sin(pi * y))



    if plot_show:
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        X = np.arange(-1, 1, 0.025)
        Y = np.arange(-1, 1, 0.025)
        X, Y = np.meshgrid(X, Y)

        #R = np.array([[grid[0].get_value(X[i][j], Y[i][j]) for j in range(80)] for i in range(80)])
        R = np.array([[100*error(X[i][j],Y[i][j])for j in range(80)]for i in range(80)])
        surf = ax.plot_surface(X, Y, R, rstride=1, cstride=1, cmap=cm.jet,
                               linewidth=0, antialiased=False)
        ax.set_zlim(-1.01, 1.01)

        ax.zaxis.set_major_locator(LinearLocator(10))
        ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

        fig.colorbar(surf, shrink=0.5, aspect=5)

        plt.show()
