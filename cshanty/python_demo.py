from cshanty.backend import ffi, lib
import numpy as np
from matplotlib import pyplot as plt

# construct a ConfigStruct
cfg = ffi.new("ConfigStruct *")
cfg.y0 = ffi.new("double[6]", [20000e3, 0.5, -0.2, 0.5, 0, 0])
cfg.y_target = ffi.new("double[6]", [25000e3, 0.2, 0.5, 0, 0.3, 0])
cfg.propulsion_model = ffi.addressof(lib, "sail_thrust")
cfg.solver = ffi.addressof(lib, "rk89")
cfg.steering_law = ffi.addressof(lib, "lyapunov_steering")
cfg.t_span = ffi.new("double[2]", [0, 1e8])
cfg.ode_rel_tol = 1e-6
cfg.ode_h0 = 1e2
cfg.guidance_tol = 3e-2
cfg.guidance_weights = ffi.new("double[5]", [5, 1, 1, 1, 1])
cfg.penalty_param = 1
cfg.min_pe = 6878e3
cfg.penalty_weight = 0

# run the mission
sol = lib.run_mission(cfg)

# copy the solution to a numpy array
t = np.frombuffer(ffi.buffer(sol.t, sol.n * ffi.sizeof("double")), dtype=np.float64)
y = np.frombuffer(
    ffi.buffer(sol.y, sol.n * 6 * ffi.sizeof("double")), dtype=np.float64
).reshape((sol.n, 6))

# print stats
print(f"n: {sol.n}")
print(f"n_fev: {sol.n_fev}")
print(f"n_step_fail: {sol.n_step_fail}")
print(f"t_final: {t[-1]}")
print(f"Num Revolutions: {y[-1, 5] / (2 * np.pi):.0f}")


# plots
plt.figure()
plt.subplot(131)
plt.plot(t, y[:, 0] / 1e3, label="p")
plt.axhline(cfg.y_target[0] / 1e3, color="k", linestyle="--", label="p_target")

plt.subplot(132)
plt.plot(t, y[:, 1], label="f")
plt.plot(t, y[:, 2], label="g")
plt.axhline(cfg.y_target[1], color="k", linestyle="--", label="f_target")
plt.axhline(cfg.y_target[2], color="k", linestyle="--", label="g_target")

plt.subplot(133)
plt.plot(t, y[:, 3], label="h")
plt.plot(t, y[:, 4], label="k")
plt.axhline(cfg.y_target[3], color="k", linestyle="--", label="h_target")
plt.axhline(cfg.y_target[4], color="k", linestyle="--", label="k_target")

plt.savefig("python_demo.png")
plt.show()

# free the memory allocated by C
lib.free(sol.t)
print("freed sol.t")
lib.free(sol.y)
print("freed sol.y")
lib.free(sol)
print("freed sol")
