import numpy as np
from scipy.optimize import minimize_scalar

from cshanty.wrapper import (
    ODESolver,
    SteeringLaw,
    PropulsionModel,
    ConfigStruct,
)
from cshanty.optimization_funcs import tof_wrt_kappa_obj


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


def tof_wrt_just_kappa_d(kappa_d):
    return tof_wrt_kappa_obj(kappa_d, cfg.kappa_feathered, cfg)


sol = minimize_scalar(
    tof_wrt_just_kappa_d,
    bounds=(np.deg2rad(30), np.deg2rad(85)),
    method="bounded",
    options={"maxiter": 1000, "disp": True},
)

print(sol)

print(f"Optimal Degraded Angle for Min TOF is {np.rad2deg(sol.x):.5f} deg")
