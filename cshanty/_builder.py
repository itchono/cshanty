from cffi import FFI

ffibuilder = FFI()

ffibuilder.cdef(
    """
typedef struct RKSolution
{
    double *t;       // timesteps
    double (*y)[6];  // function value
    int n;           // number of steps taken
    int n_fev;       // number of function evaluations
    int n_step_fail; // number of failed steps
    bool halt;       // whether the integration was halted
    bool fault;      // whether an error occurred
} RKSolution;

typedef struct ConfigStruct ConfigStruct;

typedef RKSolution(*(*ODESolver)(void(double, double *, double *, bool *, bool *, ConfigStruct *), double, double, ConfigStruct *));
typedef void (*SteeringLaw)(double, double *, ConfigStruct *, double[2]);

struct ConfigStruct
{
    double y0[6];
    double y_target[5];
    ODESolver solver;
    SteeringLaw steering_law;
    double t_span[2];
    double ode_rel_tol;
    double ode_h0;
    double guidance_tol;
    double guidance_weights[5];
    double penalty_param;
    double min_pe;
    double penalty_weight;
    double kappa_degraded;
    double kappa_feathered;
    double sail_sigma;
};

void free(void* ptr); // free memory allocated by C code

RKSolution *rk89(void (*f)(double, double[6], double[6], bool*, bool*, ConfigStruct*), double t0, double tf, ConfigStruct* cfg);
RKSolution *rk78(void (*f)(double, double[6], double[6], bool*, bool*, ConfigStruct*), double t0, double tf, ConfigStruct* cfg);
RKSolution *rk67(void (*f)(double, double[6], double[6], bool*, bool*, ConfigStruct*), double t0, double tf, ConfigStruct* cfg);
RKSolution *rk56(void (*f)(double, double[6], double[6], bool*, bool*, ConfigStruct*), double t0, double tf, ConfigStruct* cfg);
RKSolution *rk810(void (*f)(double, double[6], double[6], bool*, bool*, ConfigStruct*), double t0, double tf, ConfigStruct* cfg);
RKSolution *rk1012(void (*f)(double, double[6], double[6], bool*, bool*, ConfigStruct*), double t0, double tf, ConfigStruct* cfg);
RKSolution *rk1214(void (*f)(double, double[6], double[6], bool*, bool*, ConfigStruct*), double t0, double tf, ConfigStruct* cfg);

void sail_thrust(double t, double y[6], double angles[2], double acceleration[3], double sail_sigma);
void lyapunov_steering(double t, double y[6], ConfigStruct *cfg, double angles[2]);
void ndf_heuristic(double t, double y[6], double ideal_angles[2], ConfigStruct *cfg, double adapted_angles[2]);
void pe_penalty(double y[6], double pen_param, double rpmin, double *P, double dPdy[5]);

RKSolution *run_mission(ConfigStruct *cfg);
"""
)


ffibuilder.set_source(
    "cshanty.backend",
    """
    #include "driver.h"
    #include "constants.h"
    #include "solvers.h"
    """,
    include_dirs=["./include"],
    # extra_compile_args=["-std=c99", "-O3", "-march=native", "-ffast-math", "-lm"],
    extra_compile_args=["/O2"],
    
)

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
