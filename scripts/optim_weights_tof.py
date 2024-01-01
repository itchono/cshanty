import numpy as np
from scipy.optimize import minimize

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
    solver=ODESolver.RK89,
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


x = minimize(
    tof_wrt_weights_obj,
    np.array([1, 1, 1, 1, 1]),
    args=(cfg,),
    method="Nelder-Mead",
    bounds=((0, None), (0, None), (0, None), (0, None), (0, None)),
    options={"maxiter": 1000, "disp": True},
)

print(x)
