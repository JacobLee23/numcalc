/**
 * Source file for "../include/integral.h"
 */

#include <stdlib.h>

#include "../include/integral.h"


double delta(struct Interval **intervals, unsigned int d) {

    double a = 1., b = 1.;
    for (int i = 0; i < d; ++i) {
        a *= (*(intervals + i))->upper - (*(intervals + i))->lower;
        b *= (*(intervals + i))->n;
    }

    return a / b;

}

short endpoint(struct Interval *interval, unsigned int i, double *x) {

    if (!(0 <= i <= interval->n + 1)) { return -1; }

    double dx = (interval->upper - interval->lower) / interval->n;
    *x = interval->lower + i * dx;

    return 0;

}

short left(struct Interval *interval, unsigned int i, double *x) {

    if (!(0 <= i <= interval->n)) { return -1; }

    double dx = (interval->upper - interval->lower) / interval->n;
    *x = interval->lower + i * dx;

    return 0;

}

short right(struct Interval *interval, unsigned int i, double *x) {

    if (!(0 <= i <= interval->n)) { return -1; }

    double dx = (interval->upper - interval->lower) / interval->n;
    *x = interval->lower + (i + 1) * dx;

    return 0;

}

short midpoint(struct Interval *interval, unsigned int i, double *x) {

    if (!(0 <= i <= interval->n)) { return -1; }

    double dx = (interval->upper - interval->lower) / interval->n;
    *x = interval->lower + (i + 1.0 / 2) * dx;

    return 0;

}

static unsigned short inbounds(struct Interval **intervals, unsigned int radix, unsigned int d) {

    unsigned int bound = 1;
    for (int i = 0; i < d; ++i) { bound *= (*(intervals + i))->n; }

    return 0 <= radix <= bound;

}

static unsigned int increment(struct Interval **intervals, unsigned int radix, unsigned int d) {

    unsigned int base = (*intervals)->n;
    unsigned int msd = radix / base;

    if (d == 1) {
        return (msd >= base - 1 ? 0 : msd + 1);
    }
    return (msd >= base - 1 ? 0 : (msd + 1) * base + increment(intervals + 1, radix % base, d - 1));

}

static unsigned int *unpack(struct Interval **intervals, unsigned int radix, unsigned int d) {

    unsigned int *index;
    if (!(index = (unsigned int *)calloc(d, sizeof(unsigned int)))) { return NULL; }

    for (int i = 0; i < d; ++i) {
        *(index + i) = radix / (*(intervals + i))->n;
        radix %= (*(intervals + i))->n;
    }

    return index;

}

static short xvalue(
    struct Interval **intervals, RiemannRule *rules, unsigned int radix, unsigned int d,
    double *x
) {

    unsigned int *index;
    if (!(index = unpack(intervals, radix, d))) { return NULL; }

    for (int i = 0; i < d; ++i) {
        if ((*(rules + i))(*(intervals + i), *(index + i), x) == -1) { return -1; }
    }

    return 0;

}

short riemann(RealFunction f, struct Interval **intervals, RiemannRule *rules, unsigned int d, double *res) {

    double *x;
    if (!(x = (double *)calloc(d, sizeof(double)))) {
        res = NULL;
        return -1;
    }

    const double dv = delta(intervals, d);
    unsigned int radix = 0;

    *res = 0.;

    while (inbounds(intervals, radix, d)) {

        if (xvalue(intervals, rules, radix, d, x) == -1) {
            free(x);
            res = NULL;
            return -1;
        }

        *res += dv * f(x, d);
        radix = increment(intervals, radix, d);

    }

    free(x);

    return 0;

}

short trapezoidal(RealFunction f, struct Interval **intervals, unsigned int d, double *res) {

    RiemannRule *rules;
    if (!(rules = (RiemannRule *)calloc(d, sizeof(RiemannRule)))) {
        res = NULL;
        return -1;
    }
    for (int i = 0; i < d; ++i) {
        *(rules + i) = endpoint;
    }

    double *x;
    if (!(x = (double *)calloc(d, sizeof(double)))) {
        res = NULL;
        return -1;
    }

    const double dv = delta(intervals, d);
    const unsigned int nlegs = 2 << d;
    unsigned int nborders;
    unsigned int radix = 0;

    *res = 0.;

    while (inbounds(intervals, radix, d)) {

        if (xvalue(intervals, rules, radix, d, x) == -1) {
            free(x);
            res = NULL;
            return -1;
        }

        nborders = 0;
        for (int i = 0; i < d; ++i) {
            if (
                *(x + i) == (*(intervals + i))->lower || *(x + i) == (*(intervals + i))->upper
            ) {
                ++nborders;
            }
        }

        *res += dv * (1 << (d - nborders)) * f(x, d) / nlegs;
        radix = increment(intervals, radix, d);

    }
    
    free(x);

    return 0;

}

static PyObject *integral_delta(PyObject *self, PyObject *args) {

    PyObject *ob_intervals;
    if (!PyArg_ParseTuple(args, "OI", &ob_intervals)) { return NULL; }

    if (!PySequence_Check(ob_intervals)) {
        PyErr_SetString(PyExc_TypeError, "Expected a sequence of 'Interval' objects");
        return NULL;
    }

    unsigned int d = (unsigned int)PySequence_Size(ob_intervals);

    struct Interval **intervals;
    if (!(intervals = (struct Interval **)calloc(d, sizeof(struct Interval *)))) {
        PyErr_SetString(PyExc_MemoryError, "Failed to allocate memory");
        return NULL;
    }

    for (unsigned int i = 0; i < d; ++i) {

        PyObject *item = PySequence_GetItem(ob_intervals, i);
        if (!PyObject_TypeCheck(item, &IntervalType)) {

            PyErr_SetString(PyExc_TypeError, "Expected a sequence of 'Interval' objects");
            Py_XDECREF(item);

            for (unsigned int j = 0; j < i; ++j) { free(*(intervals + j)); }
            free(intervals);

            return NULL;

        }

        if (!(*(intervals + i) = (struct Interval *)calloc(1, sizeof(struct Interval)))) {

            PyErr_SetString(PyExc_MemoryError, "Failed to allocate memory");
            Py_DECREF(item);

            for (unsigned int j = 0; j < i; ++j) { free(*(intervals + j)); }
            free(intervals);

            return NULL;

        }

        IntervalObject *pyinterval = (IntervalObject *)item;
        (*(intervals + i))->lower = pyinterval->lower;
        (*(intervals + i))->upper = pyinterval->upper;
        (*(intervals + i))->n = pyinterval->n;
        Py_DECREF(item);

    }

    double result = delta(intervals, d);

    for (unsigned int i = 0; i < d; ++i) { free(*(intervals + i)); }
    free(intervals);

    return PyFloat_FromDouble(result);

}

static PyObject *riemann_rule(PyObject *self, PyObject *args, RiemannRule rule) {

    PyObject *ob_interval;
    unsigned int i;
    if (!PyArg_ParseTuple(args, "OI", &ob_interval, &i)) { return NULL; }

    if (!PyObject_TypeCheck(ob_interval, &IntervalType)) {
        PyErr_SetString(PyExc_TypeError, "Expected an 'Interval' object");
        return NULL;
    }
    struct Interval *py_interval = (struct Interval *)ob_interval;

    struct Interval *interval;
    if (!(interval = (struct Interval *)malloc(sizeof(struct Interval)))) {
        PyErr_SetString(PyExc_MemoryError, "Failed to allocate memory");
        return NULL;
    }

    interval->lower = py_interval->lower;
    interval->upper = py_interval->upper;
    interval->n = py_interval->n;

    double *x;
    if (!(x = (double *)malloc(sizeof(double)))) {
        PyErr_SetString(PyExc_MemoryError, "Failed to allocate memory");
        free(interval);
        return NULL;
    }

    short err;
    if ((err = rule(interval, i, x)) == -1) {
        PyErr_SetString(PyExc_ValueError, "Index of bounds of interval");
        free(interval);
        free(x);
        return NULL;
    }

    PyObject *res = PyFloat_FromDouble(*x);
    free(interval);
    free(x);

    return res;

}

static PyObject *integral_endpoint(
    PyObject *self, PyObject *args
) { return riemann_rule(self, args, endpoint); }

static PyObject *integral_left(
    PyObject *self, PyObject *args
) { return riemann_rule(self, args, left); }

static PyObject *integral_right(
    PyObject *self, PyObject *args
) { return riemann_rule(self, args, right); }

static PyObject *integral_midpoint(
    PyObject *self, PyObject *args
) { return riemman_rule(self, args, midpoint); }
