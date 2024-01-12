import numpy as np
from scipy.optimize import dual_annealing, differential_evolution

from cshanty.cases import BENCHMARK_TRANSFER, OGURI_CASE_G
from cshanty.optimization_funcs import tof_wrt_weights_obj


cfg = OGURI_CASE_G

x = dual_annealing(
    tof_wrt_weights_obj,
    bounds=((0.1, 10), (0.1, 10), (0.1, 10), (0.1, 10), (0.1, 10)),
    args=(cfg,),
    maxiter=1000,
    x0=np.array([1, 1, 1, 1, 1]),
)

# x = differential_evolution(
#     tof_wrt_weights_obj,
#     bounds=((0.1, 10), (0.1, 10), (0.1, 10), (0.1, 10), (0.1, 10)),
#     args=(cfg,),
#     maxiter=1000,
#     x0=np.array([1, 1, 1, 1, 1]),
# )

print(x)
