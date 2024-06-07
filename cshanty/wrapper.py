from dataclasses import dataclass
from enum import Enum

import numpy as np
import numpy.typing as npt

from cshanty.backend import ffi, lib


class ODESolver(Enum):
    """
    Enum of the available ODE solvers.
    """

    RK89 = ffi.addressof(lib, "rk89")
    RK78 = ffi.addressof(lib, "rk78")
    RK67 = ffi.addressof(lib, "rk67")
    RK56 = ffi.addressof(lib, "rk56")
    RK810 = ffi.addressof(lib, "rk810")
    RK1012 = ffi.addressof(lib, "rk1012")
    RK1214 = ffi.addressof(lib, "rk1214")


class SteeringLaw(Enum):
    """
    Enum of the available steering laws.
    """

    LYAPUNOV = ffi.addressof(lib, "lyapunov_steering")


@dataclass
class ConfigStruct:
    """
    Python analog of the C struct ConfigStruct.
    """

    y0: npt.NDArray[np.floating]
    y_target: npt.NDArray[np.floating]
    solver: ODESolver
    steering_law: SteeringLaw
    t_span: tuple[float, float]
    ode_rel_tol: float
    ode_h0: float
    guidance_tol: float
    guidance_weights: npt.NDArray[np.floating]
    penalty_param: float
    min_pe: float
    penalty_weight: float
    kappa_degraded: float
    kappa_feathered: float
    sail_sigma: float

    @property
    def _cstruct(self):
        """
        Convert this ConfigStruct to a C struct.
        """
        return ffi.new(
            "ConfigStruct *",
            {
                "y0": ffi.new("double[6]", self.y0.tolist()),
                "y_target": ffi.new("double[5]", self.y_target.tolist()),
                "solver": self.solver.value,
                "steering_law": self.steering_law.value,
                "t_span": ffi.new("double[2]", self.t_span),
                "ode_rel_tol": self.ode_rel_tol,
                "ode_h0": self.ode_h0,
                "guidance_tol": self.guidance_tol,
                "guidance_weights": ffi.new(
                    "double[5]", self.guidance_weights.tolist()
                ),
                "penalty_param": self.penalty_param,
                "min_pe": self.min_pe,
                "penalty_weight": self.penalty_weight,
                "kappa_degraded": self.kappa_degraded,
                "kappa_feathered": self.kappa_feathered,
                "sail_sigma": self.sail_sigma,
            },
        )


@dataclass
class MissionResult:
    """
    Python analog of the C struct RKSolution.
    """

    t: npt.NDArray[np.floating]
    y: npt.NDArray[np.floating]
    n: int
    n_fev: int
    n_step_fail: int
    halt: bool
    fault: bool

    @classmethod
    def from_cstruct(cls, cstruct):
        """
        Convert a C struct RKSolution to a MissionResult.
        """
        return cls(
            t=np.frombuffer(
                ffi.buffer(cstruct.t, cstruct.n * ffi.sizeof("double")),
                dtype=np.float64,
            ).copy(),
            y=np.frombuffer(
                ffi.buffer(cstruct.y, cstruct.n * 6 * ffi.sizeof("double")),
                dtype=np.float64,
            )
            .reshape((cstruct.n, 6))
            .copy(),
            n=cstruct.n,
            n_fev=cstruct.n_fev,
            n_step_fail=cstruct.n_step_fail,
            halt=cstruct.halt,
            fault=cstruct.fault,
        )


def run_mission(cfg: ConfigStruct) -> MissionResult:
    """
    Run the mission defined by the given ConfigStruct.
    """
    result_cstruct = lib.run_mission(cfg._cstruct)
    result = MissionResult.from_cstruct(result_cstruct)

    # free memory allocated by C code
    lib.free(result_cstruct.t)
    lib.free(result_cstruct.y)
    lib.free(result_cstruct)

    return result
