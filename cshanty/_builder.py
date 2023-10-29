from cffi import FFI

ffibuilder = FFI()

# two public functions: orbit and stateFromKeplerian
ffibuilder.cdef(
    """
typedef struct RKSolution
{
    double *t;             // timesteps
    double (*y)[6]; // function value
    int n;                 // number of steps taken
    int n_fev;             // number of function evaluations
    int n_step_fail;       // number of failed steps
} RKSolution;

typedef RKSolution(*(*ODESolver)(void(double, double *, double *), double, double, double *, double, double));
void free(void* ptr); // free memory allocated by C code

RKSolution *rk89(void (*f)(double, double[6], double[6]), double t0, double tf, double y0[6], double h0, double tol);
RKSolution *rk78(void (*f)(double, double[6], double[6]), double t0, double tf, double y0[6], double h0, double tol);
RKSolution *rk67(void (*f)(double, double[6], double[6]), double t0, double tf, double y0[6], double h0, double tol);
RKSolution *rk56(void (*f)(double, double[6], double[6]), double t0, double tf, double y0[6], double h0, double tol);
RKSolution *rk810(void (*f)(double, double[6], double[6]), double t0, double tf, double y0[6], double h0, double tol);
RKSolution *rk1012(void (*f)(double, double[6], double[6]), double t0, double tf, double y0[6], double h0, double tol);
RKSolution *rk1214(void (*f)(double, double[6], double[6]), double t0, double tf, double y0[6], double h0, double tol);


void test_ode(double t, double *y, double *dydt);
"""
)


ffibuilder.set_source(
    "cshanty.backend",
    """
    #include "solver_weights.h"
    #include "solvers.h"

    void test_ode(double t, double *y, double *dydt)
    {
        dydt[0] = -y[0];
    }
    """,
    include_dirs=["./include"],
    extra_compile_args=["-std=c99", "-O3", "-march=native", "-ffast-math", "-lm"],
)

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
