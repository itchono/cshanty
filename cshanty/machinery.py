# Wrapper for lower level functions
from cshanty.backend import ffi, lib
from cshanty.wrapper import ConfigStruct

import numpy as np
import numpy.typing as npt


def lyapunov_steering(
    t: float, y: npt.NDArray[np.floating], cfg: ConfigStruct
) -> npt.NDArray[np.floating]:
    """
    Wrapper for the C function lyapunov_steering.
    """
    angles = ffi.new("double[2]")
    y_c = ffi.new("double[6]", y.tolist())
    lib.lyapunov_steering(t, y_c, cfg._cstruct, angles)

    return np.array([angles[0], angles[1]])


def ndf_heuristic(
    t: float,
    y: npt.NDArray[np.floating],
    ideal_angles: npt.NDArray[np.floating],
    cfg: ConfigStruct,
) -> npt.NDArray[np.floating]:
    """
    Wrapper for the C function ndf_heuristic.
    """
    angles = ffi.new("double[2]")
    y_c = ffi.new("double[6]", y.tolist())
    ideal_angles_c = ffi.new("double[2]", ideal_angles.tolist())
    lib.ndf_heuristic(t, y_c, ideal_angles_c, cfg._cstruct, angles)
    return np.array([angles[0], angles[1]])


def sail_thrust(
    t: float, y: npt.NDArray[np.floating], angles: npt.NDArray[np.floating]
) -> npt.NDArray[np.floating]:
    """
    Wrapper for the C function sail_thrust.
    """
    acceleration = ffi.new("double[3]")
    y_c = ffi.new("double[6]", y.tolist())
    angles_c = ffi.new("double[2]", angles.tolist())
    lib.sail_thrust(t, y_c, angles_c, acceleration)
    return np.array([acceleration[0], acceleration[1], acceleration[2]])


def pe_penalty(
    y: npt.NDArray[np.floating], pen_param: float, rpmin: float
) -> tuple[float, npt.NDArray[np.floating]]:
    """
    Wrapper for the C function pe_penalty.
    """
    P = ffi.new("double*")
    dPdy = ffi.new("double[5]")
    y_c = ffi.new("double[6]", y.tolist())
    lib.pe_penalty(y_c, pen_param, rpmin, P, dPdy)

    return P[0], np.array([dPdy[0], dPdy[1], dPdy[2], dPdy[3], dPdy[4]])
