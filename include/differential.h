/**
 * Numerical differential calculus
 */

#define PY_SSIZE_T_CLEAN
#include <Python.h>


static struct FiniteDifference {
    double *(*first)(PyObject *f, double *x, double h, unsigned int d);
    double *(*second)(PyObject *f, double *x, double h, unsigned int d);
    double *(*nth)(PyObject *f, double *x, double h, unsigned int n, unsigned int d);
};
static enum FinDiffRule { FORWARD, BACKWARD, CENTRAL };

static double *forward_first(PyObject *f, double *x, double h, unsigned int d);
static double *forward_second(PyObject *f, double *x, double h, unsigned int d);
static double *forward_nth(PyObject *f, double *x, double h, unsigned int n, unsigned int d);
static struct FiniteDifference forward = { forward_first, forward_second, forward_nth };

static double *backward_first(PyObject *f, double *x, double h, unsigned int d);
static double *backward_second(PyObject *f, double *x, double h, unsigned int d);
static double *backward_nth(PyObject *f, double *x, double h, unsigned int n, unsigned int d);
static struct FiniteDifference backward = { backward_first, backward_second, backward_nth };

static double *central_first(PyObject *f, double *x, double h, unsigned int d);
static double *central_second(PyObject *f, double *x, double h, unsigned int d);
static double *central_nth(PyObject *f, double *x, double h, unsigned int n, unsigned int d);
static struct FiniteDifference central = { central_first, central_second, central_nth };

static double *dquotient(
    PyObject *f, double *x, double h, unsigned int n, unsigned int d,
    struct FiniteDifference *findiff
);

static PyObject *differential_dquotient(PyObject *self, PyObject *args);

static PyMethodDef DifferentialMethods[] = {
    {"dquotient", differential_dquotient, METH_VARARGS, NULL},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef differential_module = {
    PyModuleDef_HEAD_INIT, "differential", NULL, -1, DifferentialMethods
};

PyMODINIT_FUNC PyInit_differential() { return PyModule_Create(&differential_module); }
