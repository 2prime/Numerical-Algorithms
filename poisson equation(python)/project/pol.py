'''
多项式:
    多项式加法,乘法,幂次
    多项式求值,求导,积分
    多项式taylor展开
切比雪夫多项式
多项式插值
用切比雪夫多项式的根做的多项式插值
'''
from math import cos,pi

def pop_zero(list):
    '''删掉列表末尾的zero'''
    list.reverse()
    # 一个做迭代器一个做结果
    inter =  list.copy()
    result = list.copy()
    list.reverse()
    for i in inter:
        if i != 0:
            break
        result = result[1:]
    result.reverse()
    return result

class poly(object):
    def __init__(self,list,x = 0):
        self.poly = list
        self.x = x

    def __add__(self, other):
        '''Return self+other.'''
        # 做拷贝防止修改
        poly_1 = self.poly.copy()
        poly_2 = other.poly.copy()
        pol1 = poly(poly_1,self.x).taylor(0)
        pol2 = poly(poly_2,other.x).taylor(0)
        poly_1 = pop_zero(pol1.poly.copy())
        poly_2 = pop_zero(pol2.poly.copy())

        # 把短的变长
        if len(poly_1) >= len(poly_2):
            poly_2 += [0]*(len(poly_1)-len(poly_2))
        else:
            poly_1 += [0]*(len(poly_2)-len(poly_1))

        # result存储结果
        result = []
        for i,j in zip(poly_1,poly_2):
            result.append(i + j)
        return poly(result)

    def __eq__(self, other):
        '''Return self == other.'''
        poly1 = pop_zero(self.poly)
        poly2 = pop_zero(other.poly)
        # 在0处打开
        pol1 = poly(poly1,self.x).taylor(0)
        pol2 = poly(poly2,other.x).taylor(0)
        poly1 = pop_zero(pol1.poly)
        poly2 = pop_zero(pol2.poly)
        if len(poly1) != len(poly2):
            return False
        for i,j in zip(poly1,poly2):
            if i != j:
                return False
        return True

    def __ne__(self, other):
        '''Return self != other.'''
        poly1 = pop_zero(self.poly)
        poly2 = pop_zero(other.poly)
        if len(poly1) != len(poly2):
            return True
        for i,j in zip(poly1,poly2):
            if i != j:
                return True
        return False

    def __mul__(self, other):
        '''Return self * other'''
        # 做拷贝
        poly1 = pop_zero(self.poly)
        poly2 = pop_zero(other.poly)
        # 在0处打开
        pol1 = poly(poly1,self.x).taylor(0)
        pol2 = poly(poly2,other.x).taylor(0)
        poly1 = pop_zero(pol1.poly)
        poly2 = pop_zero(pol2.poly)
        len1 = len(poly1)
        len2 = len(poly2)
        # 存储答案
        result = []
        for i in range(len1 + len2 - 1):
            num = 0
            j = 0
            while j <= i:
                if j < len1 and i - j < len2:
                    num += poly1[j] * poly2[i-j]
                j += 1
            result.append(num)
        return poly(result)

    def __pow__(self, power, modulo=None):
        '''Return self**power'''
        result = poly([1])
        for i in range(power):
            result *= self
        return result

    def __gt__(self, other):
        '''Return self > other'''
        if self.deg() > other.deg():
            return True
        if self.deg() < other.deg():
            return False
        poly1 = self.poly.copy()
        poly2 = self.poly.copy()
        poly1.reverse()
        poly2.reverse()
        return poly1 > poly2

    def __str__(self):
        '''return str(self)'''
        if self.x == 0:
            para = "x"
        else:
            para = "(x-%.2f)"%(self.x) if self.x >= 0 else "(x+%.2f)"%(-self.x)
        self.poly = pop_zero(self.poly)
        if self.poly[0] != 0:
            result = str(self.poly[0]) + "+"
        else:
            result = ""
        for i in range(1,len(self.poly)-1):
            result += str(self.poly[i])+"%s^%d+"%(para,i)

        return result + str(self.poly[len(self.poly)-1])+"%s^%d"%(para,len(self.poly)-1) if len(self.poly)!=1 else result[:-1]

    def deg(self):
        '''计算多项式的次数'''
        return len(pop_zero(self.poly)) - 1

    def eval(self, x):
        '''计算多项式在x时候的值'''
        x = x - self.x
        val = 0
        for i in range(0, len(self.poly)):
            val += self.poly[i] * x ** i
        return val

    def der(self):
        '''多项式求导'''
        poly_Original = pop_zero(self.poly)
        result = []
        for i in range(1,len(poly_Original)):
            result.append(i * poly_Original[i])
        return poly(result,self.x)

    def inter(self):
        '''多项式积分'''
        poly_ori = pop_zero(self.poly)
        result = [0]
        for i in range(len(poly_ori)):
            result.append(poly_ori[i]/(i + 1))
        return poly(result,self.x)

    def taylor(self,x):
        '''taylor(self,x) -> self 在x的taylor打开'''

        # 阶乘函数
        def factor(runingtime,par = 1):
            return par if runingtime <= 1 else factor(runingtime - 1,runingtime * par)

        pol = poly(self.poly.copy(),self.x)
        result = []
        for i in range(self.deg() + 1):
            result.append(pol.eval(x)/factor(i))
            pol = pol.der()
        return poly(result,x)

def add(*args):
    sum = poly([0])
    for item in args:
        sum += item
    return sum

# horner算法计算函数茶之
def interpolation(points,values):
    '''
    :param points: 插值点
    :param values: 对应函数值
    :return: 插值多项式
    '''
    n = len(points)
    assert n == len(values) ,"插值位置列表和插值取值列表应该长度相等"
    newton_coefficient = []
    newton_coefficient.append(values[0])
    for k in range(1,n):
        dis = points[k] - points[k-1]
        u = newton_coefficient[k-1]
        #计算P(k)(x_k)
        for i in range(k-2,-1,-1):
            u = u * (points[k]-points[i]) + newton_coefficient[i]
            dis *= (points[k]-points[i])
        newton_coefficient.append((values[k]-u)/dis)
    #print(newton_coefficient)

    def P(t):
        pol = poly([0,1],points[0])
        for i in range(t):
            pol *= poly([0,1],points[i+1])
        return pol

    result = poly([newton_coefficient[0]])*poly([1])
    for i in range(n-1):
        result += poly([newton_coefficient[i+1]])*P(i)
    return result


# 用remez算法计算最佳一致逼近



def Chebyshev(runingtimes,pol1=poly([1]),pol2=poly([0,1])):
    '''
    :param runingtimes: 多项式次数
    :return: 切比雪夫多项式
    '''
    return pol1 if runingtimes == 0 else Chebyshev(runingtimes-1,pol2,pol1*poly([-1])+pol2*poly([0,2]))

def Chebyshev_root(n,a=-1,b=1):
    '''
    :param n: 点数
    :param a: 左端点
    :param b: 右端点
    :return: 切比雪夫插值点
    '''
    result = []
    for i in range(n+1):
        x = cos((2*i+1)*(pi/(2*n+2)))
        result.append((x+1)*(b-a)/2+a)
    return result

def Chebyshev_interpolation(func,a,b,n):
    '''
    :param func: 被插函数
    :param a: 插值左端点
    :param b: 插值右端点
    :param n: 插值次数
    :return:插值函数
    '''
    points = Chebyshev_root(n,a,b)
    return interpolation(points,list(map(func,points)))



if __name__ ==  "__main__":
    print(Chebyshev(6).eval(Chebyshev_root(5)[3]))
