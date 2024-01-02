from cshanty.machinery import pe_penalty
import pytest
import numpy as np


def test_pe_penalty():
    y = np.array([1.162498631250000e7, 0.725, 0, 0, 0, 0])
    pen_param = 1
    rpmin = 6878e3
    P, dPdy = pe_penalty(y, pen_param, rpmin)

    assert P == pytest.approx(1.020396781144553)
    assert dPdy == pytest.approx(
        np.array([-0.000000086003833, 0.579590368512825, 0, 0, 0])
    )
