from poisson import  poisson
from math import sin,pi
from random import randrange


def Restriction_operator(poisson_2):
    assert  not poisson_2.N%2
    result_value = []
    Residual = poisson_2.getResidual()
    get_value = lambda x, y: Residual[x * (poisson_2.N - 1) + y]
    for i in range(poisson_2.N//2-1):
        for j in range(poisson_2.N//2-1):
            value = (get_value(2*i,2*j)+get_value(2*i+2,2*j)+get_value(2*i+2,2*j+2)+get_value(2*i,2*j+2)+\
                     4*get_value(2*i+1,2*j+1)+\
                     2*(get_value(2*i+1,2*j)+get_value(2*i+1,2*j+2)+get_value(2*i,2*j+1)+get_value(2*i+2,2*j+1)))/16
            result_value.append(value)
    result = poisson(poisson_2.func, poisson_2.N // 2, False, result_value)
    assert len(result_value) == len(result.grids_value)
    return result

def Lifting_operator(poisson_1):
    result = poisson(poisson_1.func, poisson_1.N * 2)
    result_value = []
    def isborder(i,j):
        return i==-1 or i == poisson_1.N-1 or j == -1 or j == poisson_1.N-1
    def get_value(i,j):
        if isborder(i,j):
            return  0
        return poisson_1.grids_value[i * (poisson_1.N - 1) + j]
    for i in range(poisson_1.N*2 - 1):
        for j in range(poisson_1.N*2 - 1):
            b_isodd1 = not i%2
            b_isodd2 = not j%2
            if not b_isodd1 and not b_isodd2:
                result_value.append(get_value(int(i/2),int(j/2)))
            if not b_isodd1 and b_isodd2:
                result_value.append((get_value(int(i/2),int(j/2))+get_value(int(i/2),int(j/2)-1))/2)
            if b_isodd1 and not b_isodd2:
                result_value.append((get_value(int(i/2),int(j/2))+get_value(int(i/2)-1,int(j/2)))/2)
            if b_isodd2 and b_isodd1:
                result_value.append((get_value(int(i/2),int(j/2)-1)+get_value(int(i/2)-1,int(j/2)-1)+ \
                                     get_value(int(i/2)-1, int(j/2))+get_value(int(i/2),int(j/2)))/4)
    assert len(result_value)==len(result.grids_value)
    result.grids_value = result_value
    return result




if __name__ == "__main__":
    grid = poisson(lambda x, y:200000*sin(pi*x)*sin(pi*y), 6)
    grid.grids_value=[randrange(10) for i in range(25)]
    print(grid)
    print(Restriction_operator(grid))

    grid = poisson(lambda x, y: 200000 * sin(pi * x) * sin(pi * y), 3)
    grid.grids_value=[randrange(10) for i in range(4)]
    print(grid)
    print(Lifting_operator(grid))








