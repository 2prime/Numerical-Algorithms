from poisson import *
from grid import *
from time import time
import pol


debug = True
error_cal = True

t = time()


grid = poisson(lambda x, y: 2 * pi * pi * sin(pi * x) * sin(pi * y), 32)
for i in range(100):
    grid.Gauss()

grid = Lifting_operator(grid)
for i in range(100):
    grid.Gauss()


grid = Lifting_operator(grid)
for i in range(50):
    grid.Gauss()


grid = Restriction_operator(grid)
for i in range(100):
    grid.Gauss()

grid = Lifting_operator(grid)
for i in range(50):
    grid.Gauss()




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
    print("标准值积分为%f"%((pol2.eval(1) - pol2.eval(-1))**(1/2)))

else:
    for i in range(300):
        x = random() * 2 - 1
        y = random() * 2 - 1
        print("测试点:(%10f,%10f),|| 相对误差:%10f,|| 准确值:%10f,数值解:%10f" % (x, y, error_rlt(x, y),standard(x,y),grid.get_value(x, y)))

print("用时:%f"%(time()-t))
