from math import sin,pi
from random import random
import time

import square


def print_matrix(matrix):
    n = len(matrix)
    for i in range(n):
        for j in range(n):
            if j in matrix[i]:
                print("%2d"%(matrix[i][j]),end="  ")
            else:
                print("%2d"%(0),end="  ")
        print("")

def turnMatrix(matrix):
    n = len(matrix)
    result = []
    for i in range(n):
        result_col = []
        for j in range(n):
            if j in matrix[i]:
                result_col.append(matrix[i][j])
            else:
                result_col.append(0)
        result.append(result_col)
    return square.square_matrix(result)


class poisson():
    """used to calculate the numerical solution of PDE u_xx+u_yy=f(x,y)"""
    def __init__(self, func, N,state =True,Residual=[]):
        self.N = N
        # Ture为对func进行GS迭代,False为对残量进行GS迭代
        self.state = state
        self.Residual = Residual
        self.grid_length = 2 / N
        h = 1/(self.grid_length**2)
        self.func = func
        self.grids_value = [0 for i in range((N - 1)**2)]
        # 记录差分格式的矩阵,第一个列表记录非零点的值,第二个列表记录非零点的位置
        """
            | A_n    cI_n                            |
            | cI_n   A_n   cI_n                      |
        A = |        cI_n  A_n   cI_n                |
            |                                        |
            |                       cI_n   A_n   cI_n|
            |                              cI_n   A_n|
        """
        A_value_1 = [[4*h,-h,-h] if not i else[-h,4*h,-h,-h] if (i+1)%(N-1)else[-h,4*h,-h] for i in range(N - 1)]
        A_value_2 = [[-h,4*h,-h,-h] if not i else[-h,-h,4*h,-h,-h] if (i+1)%(N-1)else[-h,-h,4*h,-h] for i in range(N - 1)]
        A_value_3 = [[-h,4*h,-h] if not i else[-h,-h,4*h,-h] if (i+1)%(N-1)else[-h,-h,4*h] for i in range(N - 1)]
        self.matrix_value = A_value_1+A_value_2*(N-3)+A_value_3
        self.matrix_position=[]
        def position(i,j):
            """:return i*(N-1)+j行非零元素的位置"""

            def A_postion(i,j):
                """:return 第i行的A_n的第j行的非零元素的位置"""
                if j:
                    if (j+1)%(N-1):
                        return [i * (N - 1) + j - 1, i * (N - 1) + j, i * (N - 1) + j + 1]
                    return [i * (N - 1) + j - 1, i * (N - 1) + j]
                return [i * (N - 1) + j, i * (N - 1) + j + 1]

            if i :
                if (i+1)%(N-1):
                    return [(i-1) * (N - 1) + j]+A_postion(i,j)+[(i+1) * (N - 1) + j]
                return [(i-1) * (N - 1) + j]+A_postion(i,j)
            return A_postion(i,j)+[(i+1) * (N - 1) + j]
        for i in range(N-1):
            for j in range(N-1):
                self.matrix_position.append(position(i,j))
        self.matrix=[]
        for i,j in zip(self.matrix_position,self.matrix_value):
            result = {}
            for x,y in zip(i,j):
                result[x]=y
            self.matrix.append(result)
        #print("计算完了矩阵")
        if __name__ == "__main__" and self.N <= 8:
            print_matrix(self.matrix)

    def get_func_value(self,x):
        return self.func(-1+self.grid_length*(x//(self.N-1)+1),-1+self.grid_length*(x%(self.N-1)+1))


    def Gauss(self):
        """一次高斯迭代更新格点上的值"""
        if self.state:
            for i in range((self.N - 1) ** 2):
                result = 0
                for j in self.matrix[i]:
                    if i != j:
                        result += self.matrix[i][j] * self.grids_value[j]
                self.grids_value[i] = (self.get_func_value(i) - result) / self.matrix[i][i]
                #print("计算完了1行")
        else:
            self.Gauss_Residual()

    def Gauss_Residual(self):
        """误差提升以后,对残量计算GS迭代"""
        for i in range((self.N - 1) ** 2):
            result = 0
            for j in self.matrix[i]:
                if i != j:
                    result += self.matrix[i][j] * self.grids_value[j]
            self.grids_value[i] = (self.Residual[i] - result) / self.matrix[i][i]
            #print("计算完了1行")
        return self.grids_value


    def get_value(self,x,y):
        grid_x,x_0=int((x+1)/self.grid_length),(x+1)-(int((x+1)/self.grid_length)*self.grid_length)
        grid_y,y_0=int((y+1)/self.grid_length),(y+1)-(int((y+1)/self.grid_length)*self.grid_length)
        def get_value(x,y):
            if x == -1 or x == self.N-1 or y == -1 or y == self.N-1:
                return 0
            return self.grids_value[x*(self.N-1)+y]
        leftdown = get_value(grid_x-1,grid_y-1)
        rightdown = get_value(grid_x,grid_y-1)
        rightup = get_value(grid_x,grid_y)
        leftup = get_value(grid_x,grid_y-1)
        """
        leftdown = sin(pi*(-1+grid_x*self.grid_length))*sin(pi*(-1+grid_y*self.grid_length))
        rightdown = sin(pi*(-1+(grid_x+1)*self.grid_length))*sin(pi*(-1+grid_y*self.grid_length))
        rightup = sin(pi*(-1+(grid_x+1)*self.grid_length))*sin(pi*(-1+(grid_y+1)*self.grid_length))
        leftup = sin(pi*(-1+grid_x*self.grid_length))*sin(pi*(-1+(grid_y+1)*self.grid_length))
"""
        x1 = self.grid_length - x_0
        y1 = self.grid_length - y_0
        return (leftdown*x1*y1+rightdown*x_0*y1+rightup*x_0*y_0+leftup*x1*y_0)/(self.grid_length**2)

    def getanswer(self):
        answer = turnMatrix(self.matrix).result([self.Residual])
        assert len(answer)==len(self.grids_value)
        self.grids_value = answer

    def getResidual(self):
        """残量F-Au的计算,return一个列表"""
        result_value = []
        for i in range((self.N - 1) ** 2):
            result = 0
            for j in self.matrix[i]:
                    result += self.matrix[i][j]*self.grids_value[j]
            if self.state:
                result_value.append(self.get_func_value(i)-result)
            else:
                result_value.append(self.Residual[i]-result)
        assert len(result_value) == len(self.grids_value)
        #print(result_value)
        #return [0 for i in range(len(result_value))]
        return result_value

    def update_error(self,lifted_error):
        assert self.N == lifted_error.N
        for i in range((self.N-1)**2):
            self.grids_value[i] += lifted_error.grids_value[i]
        return True

    def __str__(self):
        num = 1
        string = ""
        for i in self.grids_value:
            string += "%f  "%(i)
            if not num%(self.N-1):
                string += "\n"
            num+=1
        return string



if __name__ == "__main__":
    grid = poisson(lambda x,y:2*pi*pi*sin(pi*x)*sin(pi*y),4)
    #print(grid)
    #grid = poisson(lambda x,y:sin(pi*x)*sin(pi*y),256)
    for i in range(1):
        t = time.time()
        grid.Gauss()
        print("完成了Gauss迭代%d次,用时间%f"%(i+1,time.time()-t))
        #print(grid)
    error = lambda x,y:1-(grid.get_value(x,y)/(sin(pi*x)*sin(pi*y)))
    for i in range(300):
        x = random()*2-1
        y = random()*2-1
        print("测试点:(%f,%f),相对误差:%f"%(x,y,error(x,y)))




