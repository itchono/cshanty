import sys

import numpy as np
from scipy.optimize import dual_annealing

from cshanty import cases
from cshanty.optimization_funcs import tof_wrt_all_obj

if __name__ == "__main__":
    cfgs = [cases.BENCHMARK_TRANSFER, cases.OGURI_CASE_E]
    cfg_names = ["CASE_A", "CASE_E"]

    for cfg, cfg_name in zip(cfgs, cfg_names):
        for sail_sigma in [0.002, 0.005, 0.01, 0.05, 1 / 540]:
            cfg.sail_sigma = sail_sigma
            cfg.t_span = (0, 1.5e8)

            sys.stdout = open(
                f"optim_weights_tof_{cfg_name}_{sail_sigma:.4f}.txt", "w", buffering=1
            )

            x = dual_annealing(
                tof_wrt_all_obj,
                bounds=(
                    (55, 90),
                    (0.1, 10),
                    (0.1, 10),
                    (0.1, 10),
                    (0.1, 10),
                    (0.1, 10),
                ),
                args=(cfg,),
                maxiter=200,
                x0=np.array([60, 1, 1, 1, 1, 1]),
            )

            print(x)
