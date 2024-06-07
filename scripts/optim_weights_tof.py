import sys

import numpy as np
from scipy.optimize import dual_annealing

from cshanty import cases
from cshanty.optimization_funcs import tof_wrt_all_obj

cfg = cases.BENCHMARK_TRANSFER
cfg.sail_sigma = 0.002

sys.stdout = open(f"optim_weights_tof_{cfg.sail_sigma}.txt", "w", buffering=1)

x = dual_annealing(
    tof_wrt_all_obj,
    bounds=((55, 90), (0.1, 10), (0.1, 10), (0.1, 10), (0.1, 10), (0.1, 10)),
    args=(cfg,),
    maxiter=200,
    x0=np.array([60, 1, 1, 1, 1, 1]),
)

print(x)
