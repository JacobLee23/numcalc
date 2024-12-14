/**
 * Numerical integral calculus
 */

#define PY_SSIZE_T_CLEAN
#include <Python.h>


struct Interval {
    double lower;
    double upper;
    unsigned int n;
};
struct Interval *parse_interval(PyObject *ob_interval);
struct Interval **parse_intervals(PyObject *ob_intervals, unsigned int *d);

double delta(struct Interval **intervals, unsigned int d);
short endpoint(struct Interval *interval, unsigned int i, double *x);

typedef short (*RiemannRule)(struct Interval *interval, unsigned int i, double *x);
short left(struct Interval *interval, unsigned int i, double *x);
short right(struct Interval *interval, unsigned int i, double *x);
short midpoint(struct Interval *interval, unsigned int i, double *x);

enum RiemannRules { LEFT, RIGHT, MIDPOINT };
enum RiemannRules *parse_rrules(PyObject *ob_rrules);

short riemann(
    PyObject *f, struct Interval **intervals, enum RiemannRules *rrules, unsigned int d, double *res
);
short trapezoidal(PyObject *f, struct Interval **intervals, unsigned int d, double *res);

typedef struct {
    PyObject_HEAD
    double lower;
    double upper;
    unsigned int n;
} IntervalObject;

static PyTypeObject IntervalType = {
    .ob_base = PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "integral.Interval",
    .tp_doc = NULL,
    .tp_basicsize = sizeof(IntervalObject),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_new = PyType_GenericNew,
};

static PyObject *integral_delta(PyObject *self, PyObject *args);
static PyObject *integral_endpoint(PyObject *self, PyObject *args);

static PyObject *integral_left(PyObject *self, PyObject *args);
static PyObject *integral_right(PyObject *self, PyObject *args);
static PyObject *integral_midpoint(PyObject *self, PyObject *args);

static PyObject *integral_riemann(PyObject *self, PyObject *args);
static PyObject *integral_trapezoidal(PyObject *self, PyObject *args);

static PyMethodDef IntegralMethods[] = {
    {"delta", integral_delta, METH_VARARGS, NULL},
    {"endpoint", integral_endpoint, METH_VARARGS, NULL},
    {"left", integral_left, METH_VARARGS, NULL},
    {"right", integral_right, METH_VARARGS, NULL},
    {"midpoint", integral_midpoint, METH_VARARGS, NULL},
    {"riemann", integral_riemann, METH_VARARGS, NULL},
    {"trapezoidal", integral_trapezoidal, METH_VARARGS, NULL},
    {NULL, NULL, 0, NULL}
};

static PyModuleDef integral_module = {
    PyModuleDef_HEAD_INIT, "integral", NULL, -1, IntegralMethods
};

PyMODINIT_FUNC PyInit_integral() {

    PyObject *m = PyModule_Create(&integral_module);
    if (!m) { return NULL; }

    if (
        PyType_Ready(&IntervalType) < 0
        || PyModule_AddObjectRef(m, "Interval", (PyObject *) &IntervalType) < 0
    ) {
        Py_DECREF(m);
        return NULL;
    }
    if (
        PyModule_AddIntConstant(m, "LEFT", LEFT) < 0
        || PyModule_AddIntConstant(m, "RIGHT", RIGHT) < 0
        || PyModule_AddIntConstant(m, "MIDPOINT", MIDPOINT) < 0
    ) {
        Py_DECREF(m);
        return NULL;
    }

    return m;

}
