from LazyCallback import *
from get_tsplib import get_tsplib

name = 'fl417'
coords, dist = get_tsplib(name)

tour_lazy, obj_lazy, time_lazy = tsp_lazy_cplex(coords, dist)
plot_tour(coords, tour_lazy, title=f'Optimal Tour {name}')