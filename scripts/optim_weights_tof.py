import sys

import numpy as np
from scipy.optimize import dual_annealing

from cshanty import cases
from cshanty.optimization_funcs import tof_wrt_all_obj
import argparse

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument("case")
    
    args = parser.parse_args()
    
    cfgdict = {
        "CASE_A": cases.BENCHMARK_TRANSFER,
        "CASE_E": cases.OGURI_CASE_E
    }
    
    cfg = cfgdict[args.case]

    for sail_sigma in [0.05, 1 / 540 / 0.85]:
        cfg.sail_sigma = sail_sigma
        cfg.t_span = (0, 1.5e9)

        sys.stdout = open(
            f"optim_weights_tof_{args.case}_{sail_sigma:.4f}.txt", "w", buffering=1
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
