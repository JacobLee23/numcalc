/**
 * Numerical differential calculus
 */

#define PY_SSIZE_T_CLEAN
#include <Python.h>

#ifndef TYPES_H
#define TYPES_H
#include "../include/types.h"
#endif


static struct FiniteDifference {
    double *(*first)(struct RealFunction *f, double *x, double h, unsigned int d);
    double *(*second)(struct RealFunction *f, double *x, double h, unsigned int d);
    double *(*nth)(struct RealFunction *f, double *x, double h, unsigned int n, unsigned int d);
};

static double *forward_first(struct RealFunction *f, double *x, double h, unsigned int d);
static double *forward_second(struct RealFunction *f, double *x, double h, unsigned int d);
static double *forward_nth(struct RealFunction *f, double *x, double h, unsigned int n, unsigned int d);
static struct FiniteDifference forward = { forward_first, forward_second, forward_nth };

static double *backward_first(struct RealFunction *f, double *x, double h, unsigned int d);
static double *backward_second(struct RealFunction *f, double *x, double h, unsigned int d);
static double *backward_nth(struct RealFunction *f, double *x, double h, unsigned int n, unsigned int d);
static struct FiniteDifference backward = { backward_first, backward_second, backward_nth };

static double *central_first(struct RealFunction *f, double *x, double h, unsigned int d);
static double *central_second(struct RealFunction *f, double *x, double h, unsigned int d);
static double *central_nth(struct RealFunction *f, double *x, double h, unsigned int n, unsigned int d);
static struct FiniteDifference central = { central_first, central_second, central_nth };

static enum FinDiffRule { FORWARD, BACKWARD, CENTRAL };

static double *dquotient(
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

PyMODINIT_FUNC PyInit_differential() { return PyModule_Create(&differential_module); }
