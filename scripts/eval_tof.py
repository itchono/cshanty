import sys

import numpy as np
from scipy.optimize import dual_annealing

from cshanty import cases
from cshanty.optimization_funcs import tof_wrt_all_obj
from cshanty.wrapper import (
    ConfigStruct,
    run_mission,
)

base_cfg = cases.OGURI_CASE_E
params = np.array(
    [57.85374499, 4.29433902, 2.75374472, 1.59842671, 2.75374472, 6.28984581]
)
base_cfg.kappa_degraded = np.deg2rad(params[0])
base_cfg.guidance_weights = params[1:]
base_cfg.sail_sigma = 0.05
base_cfg.t_span = (0, 1.5e9)
sol = run_mission(base_cfg)

tof_d = sol.t[-1] / 86400
n_rev = sol.y[-1, -1] / (2 * np.pi)
print(f"ToF: {tof_d:.0f} d, nRev: {n_rev:.0f}")
