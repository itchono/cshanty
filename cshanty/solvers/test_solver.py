from cshanty.backend import ffi, lib
import numpy as np
from matplotlib import pyplot as plt

solution = lib.rk67(ffi.addressof(lib, "test_ode"), 0, 10, [1, 0, 0, 0, 0, 0], 5, 1e-5)
print(solution.n)

t_buf = ffi.buffer(solution.t, solution.n * 8)
t = np.frombuffer(t_buf, dtype=np.float64).copy()
t.shape = (solution.n,)

y_buf = ffi.buffer(solution.y, solution.n * 8 * 6)
y = np.frombuffer(y_buf, dtype=np.float64).copy()
y.shape = (solution.n, 6)

lib.free(solution.y)
lib.free(solution.t)
print("Freed")

print(t)
print(y[:, 0])
