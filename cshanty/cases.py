import numpy as np

from cshanty.wrapper import ConfigStruct, ODESolver, SteeringLaw

BENCHMARK_TRANSFER = ConfigStruct(
    y0=np.array([20000e3, 0.5, -0.2, 0.5, 0, 0]),
    y_target=np.array([25000e3, 0.2, 0.5, 0, 0.3]),
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
    kappa_degraded=np.deg2rad(64),
    kappa_feathered=np.deg2rad(90),
    sail_sigma=0.005,
)

OGURI_CASE_G = ConfigStruct(
    y0=np.array([1.162498631250000e7, 0.725, 0, 0, 0, 0]),
    y_target=np.array([42165000, 0, 0, 0, -1]),
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
    kappa_degraded=np.deg2rad(64),
    kappa_feathered=np.deg2rad(90),
    sail_sigma=0.005,
)
