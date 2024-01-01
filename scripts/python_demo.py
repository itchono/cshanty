from cshanty.wrapper import (
    ConfigStruct,
    ODESolver,
    SteeringLaw,
    PropulsionModel,
    run_mission,
)
import numpy as np
from matplotlib import pyplot as plt

# construct mission case
cfg = ConfigStruct(
    y0=np.array([1.162498631250000e7, 0.725, 0, 0, 0, 0]),
    y_target=np.array([42165000, 0, 0, 0, -1]),
    propulsion_model=PropulsionModel.SAIL_THRUST,
    solver=ODESolver.RK67,
    steering_law=SteeringLaw.LYAPUNOV,
    t_span=(0, 1e8),
    ode_rel_tol=1e-6,
    ode_h0=1e2,
    guidance_tol=1e-2,
    guidance_weights=np.array([1, 1, 1, 1, 1]),
    penalty_param=1,
    min_pe=6878e3,
    penalty_weight=0,
    kappa_degraded=np.deg2rad(64),
    kappa_feathered=np.deg2rad(91),
)

# run the mission
sol = run_mission(cfg)
t = sol.t
y = sol.y

# print stats
print(f"n: {sol.n}")
print(f"n_fev: {sol.n_fev}")
print(f"n_step_fail: {sol.n_step_fail}")
print(f"t_final: {sol.t[-1]}")
print(f"Num Revolutions: {sol.y[-1, 5] / (2 * np.pi):.0f}")

# plots
plt.figure()
plt.subplot(311)
plt.plot(t, y[:, 0] / 1e3, label="p")
plt.axhline(cfg.y_target[0] / 1e3, color="k", linestyle="--", label="p_target")

plt.subplot(312)
plt.plot(t, y[:, 1], label="f")
plt.plot(t, y[:, 2], label="g")
plt.axhline(cfg.y_target[1], color="k", linestyle="--", label="f_target")
plt.axhline(cfg.y_target[2], color="k", linestyle="--", label="g_target")

plt.subplot(313)
plt.plot(t, y[:, 3], label="h")
plt.plot(t, y[:, 4], label="k")
plt.axhline(cfg.y_target[3], color="k", linestyle="--", label="h_target")
plt.axhline(cfg.y_target[4], color="k", linestyle="--", label="k_target")

plt.savefig("python_demo.png")
plt.show()
