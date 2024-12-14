/**
 * Source file for "../include/integral.h"
 */

#include <stdlib.h>

#include "../include/functions.h"
#include "../include/integral.h"


struct Interval *parse_interval(PyObject *ob_interval) {

    if (!PyObject_TypeCheck(ob_interval, &IntervalType)) {
        PyErr_SetString(PyExc_TypeError, "Expected an 'Interval' object");
        return NULL;
    }
    IntervalObject *py_interval = (IntervalObject *)ob_interval;

    struct Interval *interval = (struct Interval *)malloc(sizeof(struct Interval));
    if (!interval) {
        PyErr_SetString(PyExc_MemoryError, "Failed to allocate memory");
        return NULL;
    }

    interval->lower = py_interval->lower;
    interval->upper = py_interval->upper;
    interval->n = py_interval->n;

    return interval;

}

struct Interval **parse_intervals(PyObject *ob_intervals, unsigned int *d) {

    if (!PySequence_Check(ob_intervals)) {
        PyErr_SetString(PyExc_TypeError, "Expected a sequence of 'Interval' objects");
        return NULL;
    }

    Py_ssize_t size_intervals = PySequence_Size(ob_intervals);
    if (size_intervals < 0) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to determine sequence size");
        return NULL;
    }
    *d = (unsigned int)size_intervals;

    struct Interval **intervals = (struct Interval **)calloc(*d, sizeof(struct Interval *));
    if (!intervals) {
        PyErr_SetString(PyExc_MemoryError, "Failed to allocate memory");
        return NULL;
    }
    for (unsigned int i = 0; i < *d; ++i) {
        *(intervals + i) = (struct Interval *)malloc(sizeof(struct Interval));
        if (!(*(intervals + i))) {
            PyErr_SetString(PyExc_MemoryError, "Failed to allocate memory");
            for (unsigned int j = 0; j < i; ++j) { free(*(intervals + i)); }
            free(intervals);
            return NULL;
        }
    }
    
    for (unsigned int i = 0; i < *d; ++i) {

        PyObject *item = PySequence_GetItem(ob_intervals, i);
        if (!item) {
            PyErr_SetString(PyExc_RuntimeError, "Failed to access sequence contents");
            PY_XDECREF(item);
            for (unsigned int j = 0; j < i; ++j) { free(*(intervals + i)); }
            free(intervals);
            return NULL;
        }
        if (!PyObject_TypeCheck(item, &IntervalType)) {
            PyErr_SetString(PyExc_TypeError, "Expected an 'Interval' object");
            Py_DECREF(item);
            for (unsigned int j = 0; j < i; ++j) { free(*(intervals + i)); }
            free(intervals);
            return NULL;
        }
        IntervalObject *interval = (IntervalObject *)item;

        (*(intervals + i))->lower = interval->lower;
        (*(intervals + i))->upper = interval->upper;
        (*(intervals + i))->n = interval->n;

        Py_DECREF(item);

    }

    return intervals;

}


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

enum RiemannRules *parse_rrules(PyObject *ob_rrules) {

    if (!PySequence_Check(ob_rrules)) {
        PyErr_SetString(PyExc_TypeError, "Expected a sequence of 'int' objects");
        return NULL;
    }

    Py_ssize_t size_rules = PySequence_Size(ob_rrules);
    if (size_rules < 0) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to determine sequence size");
        return NULL;
    }
    unsigned int d = (unsigned int)size_rules;

    enum RiemannRules *rrules = (enum RiemannRules *)calloc(d, sizeof(enum RiemannRules));
    if (!rrules) {
        PyErr_SetString(PyExc_MemoryError, "Failed to allocate memory");
        free(rrules);
        return NULL;
    }

    for (unsigned int i = 0; i < d; ++i) {
        PyObject *item = PySequence_GetItem(ob_rrules, i);
        if (!item) {
            PyErr_SetString(PyExc_RuntimeError, "Failed to access sequence contents");
            PY_XDECREF(item);
            free(rrules);
            return NULL;
        }
        if (!PyLong_Check(item)) {
            PyErr_SetString(PyExc_TypeError, "Expected a sequence of 'int' objects");
            Py_DECREF(item);
            free(rrules);
            return NULL;
        }
        *(rrules + i) = (enum RiemannRules)PyLong_AsInt(item);
        Py_DECREF(item);
    }

    return rrules;

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

short riemann(struct RealFunction *f, struct Interval **intervals, enum RiemannRules *rrules, unsigned int d, double *res) {

    RiemannRule *rules = (RiemannRule *)calloc(d, sizeof(RiemannRule));
    double *x = (double *)calloc(d, sizeof(double));
    if (!rules || !x) {
        free(rules); free(x);
        res = NULL;
        return -1;
    }
    
    for (unsigned int i = 0; i < d; ++i) {
        switch (*(rrules + i)) {
            case LEFT:
                *(rules + i) = left;
                break;
            case RIGHT:
                *(rules + i) = right;
                break;
            case MIDPOINT:
                *(rules + i) = midpoint;
                break;
            default:
                free(rules); free(x);
                res = NULL;
                return -1;
        }
    }

    const double dv = delta(intervals, d);
    unsigned int radix = 0;
    *res = 0.;
    while (inbounds(intervals, radix, d)) {

        if (xvalue(intervals, rules, radix, d, x) == -1) {
            free(rules); free(x);
            res = NULL;
            return -1;
        }

        *res += dv * (f->c ? f->f.c(x, d) : eval(f->f.py, x, d));
        radix = increment(intervals, radix, d);

    }

    free(rules); free(x);

    return 0;

}

short trapezoidal(struct RealFunction *f, struct Interval **intervals, unsigned int d, double *res) {

    RiemannRule *rules = (RiemannRule *)calloc(d, sizeof(RiemannRule));
    double *x = (double *)calloc(d, sizeof(double));
    if (!rules || !x) {
        free(rules); free(x);
        res = NULL;
        return -1;
    }
    for (int i = 0; i < d; ++i) { *(rules + i) = endpoint; }

    const double dv = delta(intervals, d);
    const unsigned int nlegs = 2 << d;
    unsigned int nborders;
    unsigned int radix = 0;

    *res = 0.;

    while (inbounds(intervals, radix, d)) {

        if (xvalue(intervals, rules, radix, d, x) == -1) {
            free(rules); free(x);
            res = NULL;
            return -1;
        }

        nborders = 0;
        for (int i = 0; i < d; ++i) {
            if (*(x + i) == (*(intervals + i))->lower || *(x + i) == (*(intervals + i))->upper) {
                ++nborders;
            }
        }

        *res += dv * (1 << (d - nborders)) * (f->c ? f->f.c(x, d) : eval(f->f.py, x, d)) / nlegs;
        radix = increment(intervals, radix, d);

    }
    
    free(rules); free(x);

    return 0;

}

static PyObject *integral_delta(PyObject *self, PyObject *args) {

    PyObject *ob_intervals;
    if (!PyArg_ParseTuple(args, "OI", &ob_intervals)) { return NULL; }

    if (!PySequence_Check(ob_intervals)) {
        PyErr_SetString(PyExc_TypeError, "Expected a sequence of 'Interval' objects");
        return NULL;
    }

    unsigned int d;
    struct Interval **intervals = parse_intervals(ob_intervals, &d);
    if (!intervals) { return NULL; }

    double result = delta(intervals, d);

    for (unsigned int i = 0; i < d; ++i) { free(*(intervals + i)); }
    free(intervals);

    return PyFloat_FromDouble(result);

}

static PyObject *riemann_rule(PyObject *self, PyObject *args, RiemannRule rule) {

    PyObject *ob_interval;
    unsigned int i;
    if (!PyArg_ParseTuple(args, "OI", &ob_interval, &i)) { return NULL; }

    struct Interval *interval = parse_interval(ob_interval);
    double *x = (double *)malloc(sizeof(double));
    if (!interval || !x) {
        PyErr_SetString(PyExc_MemoryError, "Failed to allocate memory");
        free(interval); free(x);
        return NULL;
    }

    short err = rule(interval, i, x);
    if (err < 0) {
        PyErr_SetString(PyExc_ValueError, "Index of bounds of interval");
        free(interval); free(x);
        return NULL;
    }

    PyObject *res = PyFloat_FromDouble(*x);
    free(interval); free(x);

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

static PyObject *integral_riemann(PyObject *self, PyObject *args) {

    PyObject *ob_f;
    PyObject *ob_intervals;
    PyObject *ob_rrules;
    if (!PyArg_ParseTuple(args, "OOO", &ob_f, &ob_intervals, &ob_rrules)) { return NULL; }

    unsigned int d;
    struct RealFunction *f = parse_function(ob_f);
    struct Interval **intervals = parse_intervals(ob_intervals, &d);
    enum RiemannRules *rrules = parse_rrules(ob_rrules);
    double *res = (double *)malloc(sizeof(double));
    if (!f || !intervals || !rrules || !res) {
        for (unsigned int i = 0; i < d; ++i) { free(*(intervals + i)); }
        free(f); free(intervals); free(rrules); free(res);
    }

    PyObject *value = (
        riemann(f, intervals, rrules, d, res) ? NULL : PyFloat_FromDouble(*res)
    );

    for (unsigned int i = 0; i < d; ++i) { free(*(intervals + i)); }
    free(f); free(intervals); free(rrules); free(res);

    return value;

}

static PyObject *integral_trapezoidal(PyObject *self, PyObject *args) {

    PyObject *ob_f;
    PyObject *ob_intervals;
    if (!PyArg_ParseTuple(args, "OO", &ob_f, &ob_intervals)) { return NULL; }

    unsigned int d;
    struct RealFunction *f = parse_function(ob_f);
    struct Interval **intervals = parse_intervals(ob_intervals, &d);
    double *res = (double *)malloc(sizeof(double));
    if(!f || !intervals || !res) {
        for (unsigned int i = 0; i < d; ++i) { free(*(intervals + i)); }
        free(f); free(intervals); free(res);
    }

    PyObject *value = (
        trapezoidal(f, intervals, d, res) ? NULL : PyFloat_FromDouble(*res)
    );

    for (unsigned int i = 0; i < d; ++i) { free(*(intervals + i)); }
    free(f); free(intervals); free(res);

    return value;

}
