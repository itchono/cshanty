import numpy as np
import numpy.typing as npt

from cshanty.wrapper import (
    ConfigStruct,
    run_mission,
)


def tof_wrt_weights_obj(
    weights: npt.NDArray[np.floating], base_cfg: ConfigStruct
) -> float:
    """
    Time of flight associated with a certain set of guidance weights.

    Parameters
    ----------
    weights : npt.NDArray[np.floating]
        The guidance weights to use.
    base_cfg : ConfigStruct
        The base configuration for the mission.

    Returns
    -------
    float
        The time of flight associated with the given guidance weights.
    """
    base_cfg.guidance_weights = weights
    # run the mission
    sol = run_mission(base_cfg)
    tof = sol.t[-1]

    # determine if the mission was successful
    if sol.fault or not sol.halt:
        print("ERROR: fault or not halt")
        return np.inf

    print(f"Tof: {tof:.0f} sec using weights {weights}")
    return tof


def tof_wrt_kappa_obj(
    kappa_degraded: float, kappa_feathered: float, base_cfg: ConfigStruct
) -> float:
    """
    Time of flight associated with NDF heuristic parameters

    Parameters
    ----------
    kappa_degraded : float
        The degraded sail angle.
    kappa_feathered : float
        The feathered sail angle.
    base_cfg : ConfigStruct
        The base configuration for the mission.

    Returns
    -------
    float
        The time of flight associated with the given guidance weights.
    """
    base_cfg.kappa_degraded = kappa_degraded
    base_cfg.kappa_feathered = kappa_feathered
    # run the mission
    sol = run_mission(base_cfg)
    tof = sol.t[-1]

    # determine if the mission was successful
    if sol.fault or not sol.halt:
        print("fault or not halt")
        return np.inf

    print(
        f"Tof: {tof:.0f} sec using kappa_d {np.rad2deg(kappa_degraded):.2f} and kappa_f {np.rad2deg(kappa_feathered):.2f}"
    )
    return tof


def tof_wrt_all_obj(params: npt.NDArray[np.floating], base_cfg: ConfigStruct) -> float:
    """
    Time of flight associated with a certain set of guidance weights AND degraded angle.

    Parameters
    ----------
    params : npt.NDArray[np.floating]
        angle [0] and weights [1:5]
    base_cfg : ConfigStruct
        The base configuration for the mission.

    Returns
    -------
    float
        The time of flight associated with the given guidance weights.
    """
    base_cfg.kappa_degraded = np.deg2rad(params[0])
    base_cfg.guidance_weights = params[1:]

    # run the mission
    sol = run_mission(base_cfg)
    tof = sol.t[-1]

    if tof - 1000 > base_cfg.t_span[1]:
        print("ERROR: tof too long")
        return np.inf

    # determine if the mission was successful
    if sol.fault or not sol.halt:
        print("ERROR: fault or not halt")
        return np.inf

    print(f"Tof: {tof:.0f} sec = {(tof / 86400):.0f} days using weights {params}")

    return tof
