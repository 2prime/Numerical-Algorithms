from math import pi,sin

from Vgrid import V_cyclic_multigrid,getError
from poisson import poisson

times1 = 1
times2 = 1

largest_grid = 9
restriciton_time = 3
V_times = 3
grid = []
grid.append(poisson(lambda x, y: 2 * pi * pi * sin(pi * x) * sin(pi * y), 2**largest_grid))

V_cyclic_multigrid(grid, times1, times2, 3, V_times)
getError(grid)