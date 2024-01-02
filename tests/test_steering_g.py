from cshanty.machinery import lyapunov_steering, ndf_heuristic
from cshanty.wrapper import (
    ConfigStruct,
    ODESolver,
    SteeringLaw,
    PropulsionModel,
)
import pytest
import numpy as np


@pytest.fixture
def dummy_cfg():
    return ConfigStruct(
        y0=np.array([1.162498631250000e7, 0.725, 0, 0, 0, 0]),
        y_target=np.array([42165000, 0, 0, 0, -1]),
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
        penalty_weight=5,
        kappa_degraded=np.deg2rad(64),
        kappa_feathered=np.deg2rad(91),
    )


@pytest.mark.parametrize(
    ("y", "y_tgt", "expected"),
    (
        (
            np.array([1, 1, 1, 1, 1, 1]),
            np.array([42165000, 0, 0, 0, -1]),
            np.array([-5.225468916656679e-07, -3.238746887282259e-06]),
        ),
    ),
)
def test_lyapunov_steering(y, y_tgt, expected, dummy_cfg: ConfigStruct):
    dummy_cfg.y_target = y_tgt

    angles = lyapunov_steering(0, y, dummy_cfg)
    assert angles == pytest.approx(expected)
