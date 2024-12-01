/**
 * Numerical differential calculus
 */

#define PY_SSIZE_T_CLEAN
#include "Python.h"

#ifndef TYPES_H
#define TYPES_H
#include "../include/types.h"
#endif


struct FiniteDifference {
    double *(*first)(struct RealFunction *f, double *x, double h, unsigned int d);
    double *(*second)(struct RealFunction *f, double *x, double h, unsigned int d);
    double *(*nth)(struct RealFunction *f, double *x, double h, unsigned int n, unsigned int d);
};

double *forward_first(struct RealFunction *f, double *x, double h, unsigned int d);
double *forward_second(struct RealFunction *f, double *x, double h, unsigned int d);
double *forward_nth(struct RealFunction *f, double *x, double h, unsigned int n, unsigned int d);

double *backward_first(struct RealFunction *f, double *x, double h, unsigned int d);
double *backward_second(struct RealFunction *f, double *x, double h, unsigned int d);
double *backward_nth(struct RealFunction *f, double *x, double h, unsigned int n, unsigned int d);

double *central_first(struct RealFunction *f, double *x, double h, unsigned int d);
double *central_second(struct RealFunction *f, double *x, double h, unsigned int d);
double *central_nth(struct RealFunction *f, double *x, double h, unsigned int n, unsigned int d);

struct FiniteDifference forward = { forward_first, forward_second, forward_nth };
struct FiniteDifference backward = { backward_first, backward_second, backward_nth };
struct FiniteDifference central = { central_first, central_second, central_nth };

double *dquotient(
    struct RealFunction *f, double *x, double h, unsigned int n, unsigned int d,
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

PyMODINIT_FUNC PyInit_differential() { return PyModuleCreate_(&differential_module); }
