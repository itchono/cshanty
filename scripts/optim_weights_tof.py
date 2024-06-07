import sys

import numpy as np
from scipy.optimize import dual_annealing

from cshanty import cases
from cshanty.optimization_funcs import tof_wrt_all_obj

# capture stdout and write to file
# sys.stdout = open("optim_weights_tof.txt", "w")

cfg = cases.OGURI_CASE_G

x = dual_annealing(
    tof_wrt_all_obj,
    bounds=((30, 90), (0.1, 10), (0.1, 10), (0.1, 10), (0.1, 10), (0.1, 10)),
    args=(cfg,),
    maxiter=1000,
    x0=np.array([60, 1, 1, 1, 1, 1]),
    # x0=np.array(
    #     [45.16114668, 2.89471723, 0.19654136, 3.99686416, 9.81140716, 2.21168978]
    # ),
)

print(x)
