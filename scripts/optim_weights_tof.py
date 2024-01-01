import numpy as np
from scipy.optimize import dual_annealing, differential_evolution

from cshanty.wrapper import (
    ODESolver,
    SteeringLaw,
    PropulsionModel,
    ConfigStruct,
)
from cshanty.optimization_funcs import tof_wrt_weights_obj


cfg = ConfigStruct(
    y0=np.array([20000e3, 0.5, -0.2, 0.5, 0, 0]),
    y_target=np.array([25000e3, 0.2, 0.5, 0, 0.3, 0]),
    propulsion_model=PropulsionModel.SAIL_THRUST,
    solver=ODESolver.RK67,
    steering_law=SteeringLaw.LYAPUNOV,
    t_span=(0, 1e8),
    ode_rel_tol=1e-6,
    ode_h0=1e2,
    guidance_tol=3e-2,
    guidance_weights=np.array([1, 1, 1, 1, 1]),
    penalty_param=1,
    min_pe=6878e3,
    penalty_weight=0,
    kappa_degraded=np.deg2rad(70),
    kappa_feathered=np.deg2rad(91),
)


# x = dual_annealing(
#     tof_wrt_weights_obj,
#     bounds=((0, 10), (0, 10), (0, 10), (0, 10), (0, 10)),
#     args=(cfg,),
#     maxiter=1000,
#     x0=np.array([1, 1, 1, 1, 1]),
# )

x = differential_evolution(
    tof_wrt_weights_obj,
    bounds=((0, 10), (0, 10), (0, 10), (0, 10), (0, 10)),
    args=(cfg,),
    maxiter=1000,
    x0=np.array([1, 1, 1, 1, 1]),
)

print(x)
