import numpy as np
import pytest

from cshanty.machinery import lyapunov_steering, ndf_heuristic
from cshanty.wrapper import (
    ConfigStruct,
    ODESolver,
    SteeringLaw,
)


@pytest.fixture
def dummy_cfg():
    return ConfigStruct(
        y0=np.array([20000e3, 0.5, -0.2, 0.5, 0, 0]),
        y_target=np.array([25000e3, 0.2, 0.5, 0, 0.3]),
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
        kappa_degraded=np.deg2rad(64),
        kappa_feathered=np.deg2rad(91),
        sail_sigma=0.005,
    )


@pytest.mark.parametrize(
    ("y", "y_tgt", "expected"),
    (
        (
            np.array([25000e3, 0.1, 0.5, 0, 0.3, 0]),
            np.array([25000e3, 0.2, 0.5, 0, 0.3]),
            np.array([0, 0.068076458257999]),
        ),
        (
            np.array([1, 1, 1, 1, 1, 1]),
            np.array([25000e3, 0.2, 0.5, 0, 0.3]),
            np.array([-0.066162809112105, -0.148247720664508]),
        ),
    ),
)
def test_lyapunov_steering(y, y_tgt, expected, dummy_cfg: ConfigStruct):
    dummy_cfg.y_target = y_tgt

    angles = lyapunov_steering(0, y, dummy_cfg)
    assert angles == pytest.approx(expected)


def test_ndf_heuristic(dummy_cfg: ConfigStruct):
    y = np.array([20000e3, 0.5, -0.2, 0.5, 0, 0])
    ideal_angles = np.array([-0.066162809112105, -0.148247720664508])
    actual_angles = ndf_heuristic(0, y, ideal_angles, dummy_cfg)

    assert actual_angles == pytest.approx(
        np.array([-0.459015332945718, -0.133381657639300])
    )
