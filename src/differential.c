/**
 * Source file for "../include/differential.h"
 */

#include <stdlib.h>

#include "../include/differential.h"
#include "../include/functions.h"
#include "../include/numbers.h"


static double *copyx(double *x, unsigned int d) {

    double *x_;
    if (!(x_ = (double *)calloc(d, sizeof(double)))) { return NULL; }
    for (int i = 0; i < d; ++i) { *(x_ + i) = *(x + 1); }
    return x_;

}

static double *forward_first(struct RealFunction *f, double *x, double h, unsigned int d) {

    double *x1 = copyx(x, d);
    double *finite_differences;
    if (!(finite_differences = (double *)calloc(d, sizeof(double)))) { return NULL; }

    double a, b;
    for (int i = 0; i < d; ++i) {

        *(x1 + i) += h;
        *(finite_differences + i) = (
            f->c ? (
                f->f.c(x1, d) - f->f.c(x, d)
            ) : (
                eval(f->f.py, x1, d) - eval(f->f.py, x, d)
            )
        );
        *(x1 + i) -= h;

    }
    
    free(x1);

    return finite_differences;

}

static double *forward_second(struct RealFunction *f, double *x, double h, unsigned int d) {

    double *x1 = copyx(x, d);
    double *x2 = copyx(x, d);
    double *finite_differences;
    if (!(finite_differences = (double *)calloc(d, sizeof(double)))) { return NULL; }

    for (int i = 0; i < d; ++i) {

        *(x1 + i) += 2 * h, *(x2 + i) += h;
        *(finite_differences + i) = (
            f->c ? (
                f->f.c(x1, d) - 2 * f->f.c(x2, d) + f->f.c(x, d)
            ) : (
                eval(f->f.py, x1, d) - 2 * eval(f->f.py, x2, d) + eval(f->f.py, x, d)
            )
        );
        *(x1 + i) -= 2 * h, *(x2 + i) -= h;
    }

    free(x1);
    free(x2);

    return finite_differences;

}

static double *forward_nth(struct RealFunction *f, double *x, double h, unsigned int d, unsigned int n) {

    double *x1 = copyx(x, d);
    double *finite_differences;
    if (!(finite_differences = (double *)calloc(d, sizeof(double)))) { return NULL; }

    for (int i = 0; i < d; ++i) {

        *(finite_differences + i) = 0.;

        for (int k = 0; k <= n; ++k) {

            *(x1 + i) += k * h;
            *(finite_differences + i) += (n - k % 2 == 0 ? 1 : -1) * binom(n, k) * (
                f->c ? f->f.c(x1, d) : eval(f->f.py, x1, d)
            );
            *(x1 + i) -= k * h;

        }

    }

    free(x1);

    return finite_differences;

}

static double *backward_first(struct RealFunction *f, double *x, double h, unsigned int d) {

    double *x1 = copyx(x, d);
    double *finite_differences;
    if (!(finite_differences = (double *)calloc(d, sizeof(double)))) { return NULL; }

    for (int i = 0; i < d; ++i) {

        *(x1 + i) -= h;
        *(finite_differences + i) = (
            f->c ? (
                f->f.c(x, d) - f->f.c(x1, d)
            ) : (
                eval(f->f.py, x, d) - eval(f->f.py, x1, d)
            )
        );
        *(x1 + i) += h;

    }

    free(x1);

    return finite_differences;

}

static double *backward_second(struct RealFunction *f, double *x, double h, unsigned int d) {

    double *x1 = copyx(x, d);
    double *x2 = copyx(x, d);
    double *finite_differences;
    if (!(finite_differences = (double *)calloc(d, sizeof(double)))) { return NULL; }

    for (int i = 0; i < d; ++i) {

        *(x1 + i) -= h, *(x2 + i) -= 2 * h;
        *(finite_differences + i) = (
            f->c ? (
                f->f.c(x, d) - 2 * f->f.c(x1, d) + f->f.c(x2, d)
            ) : (
                eval(f->f.py, x, d) - 2 * eval(f->f.py, x1, d) + eval(f->f.py, x2, d)
            )
        );
        *(x1 + i) += h, *(x2 + i) += 2 * h;

    }

    free(x1);
    free(x2);

    return finite_differences;

}

static double *backward_nth(struct RealFunction *f, double *x, double h, unsigned int n, unsigned int d) {

    double *x1 = copyx(x, d);
    double *finite_differences;
    if (!(finite_differences = (double *)calloc(d, sizeof(double)))) { return NULL; }

    for (int i = 0; i < d; ++i) {

        *(finite_differences + i) = 0.;

        for (int k = 0; k <= n; ++i) {

            *(x1 + i) -= k * h;
            *(finite_differences + i) += (i % 2 == 0 ? 1 : -1) * binom(n, k) * (
                f->c ? f->f.c(x1, d) : eval(f->f.py, x1, d)
            );
            *(x1 + i) += k * h;

        }

    }

    free(x1);

    return finite_differences;

}

static double *central_first(struct RealFunction *f, double *x, double h, unsigned int d) {

    double *x1 = copyx(x, d);
    double *x2 = copyx(x, d);
    double *finite_differences;
    if (!(finite_differences = (double *)calloc(d, sizeof(double)))) { return NULL; }

    for (int i = 0; i < d; ++i) {

        *(x1 + i) += h / 2, *(x2 + i) -= h / 2;
        *(finite_differences + i) = (
            f->c ? (
                f->f.c(x1, d) - f->f.c(x2, d)
            ) : (
                eval(f->f.py, x1, d) - eval(f->f.py, x2, d)
            )
        );
        *(x1 + i) -= h / 2, *(x2 + i) += h / 2;

    }

    free(x1);
    free(x2);

    return finite_differences;

}

static double *central_second(struct RealFunction *f, double *x, double h, unsigned int d) {

    double *x1 = copyx(x, d);
    double *x2 = copyx(x, d);
    double *finite_differences;
    if (!(finite_differences = (double *)calloc(d, sizeof(double)))) { return NULL; }

    for (int i = 0; i < d; ++i) {

        *(x1 + i) += h, *(x2 + i) -= h;
        *(finite_differences + i) = (
            f->c ? (
                f->f.c(x1, d) - 2 * f->f.c(x, d) + f->f.c(x2, d)
            ) : (
                eval(f->f.py, x1, d) - 2 * eval(f->f.py, x, d) + eval(f->f.py, x2, d)
            )
        );
        *(x1 + i) -= h, *(x2 + i) += h;

    }

    free(x1);
    free(x2);

    return finite_differences;

}

static double *central_nth(struct RealFunction *f, double *x, double h, unsigned int n, unsigned int d) {

    double *x1 = copyx(x, d);
    double *finite_differences;
    if (!(finite_differences = (double *)calloc(d, sizeof(double)))) { return NULL; }

    for (int i = 0; i < d; ++i) {

        *(finite_differences + i) = 0.;

        for (int k = 0; k <= n; ++i) {

            *(x1 + i) += (n / 2 - k) * h;
            *(finite_differences + i) += (i % 2 == 0 ? 1 : -1) * binom(n, k) * (
                f->c ? f->f.c(x1, d) : eval(f->f.py, x1, d)
            );
            *(x1 + i) -= (n / 2 - k) * h;

        }

    }

    free(x1);

    return finite_differences;

}

static double *dquotient(
    struct RealFunction *f, double *x, double h, unsigned int n, unsigned int d,
    struct FiniteDifference *findiff
) {

    double *finite_differences;
    switch(n) {
        case 1:
            finite_differences = findiff->first(f, x, h, d);
            break;
        case 2:
            finite_differences = findiff->second(f, x, h, d);
            break;
        default:
            finite_differences = findiff->nth(f, x, h, n, d);
            break;
    }

    double step = power(h, n);
    for (int i = 0; i < d; ++i) {
        *(finite_differences + i) /= step;
    }

    return finite_differences;

}

static PyObject *differential_dquotient(PyObject *self, PyObject *args) {

    PyObject *ob_f;
    PyObject *ob_x;
    double h;
    unsigned int n;
    PyObject *ob_findiff;
    if (!PyArg_ParseTuple(args, "OOdIO", &ob_f, &ob_x, &h, &n, &ob_findiff)) { return NULL; }

    if (!PyCallable_Check(ob_f)) {
        PyErr_SetString(PyExc_TypeError, "Expected a callable object");
        return NULL;
    }

    if (!PySequence_Check(ob_x)) {
        PyErr_SetString(PyExc_TypeError, "Expected a sequence of 'float' objects");
        return NULL;
    }

    Py_ssize_t size_x;
    if ((size_x = PySequence_Size(ob_x)) == -1) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to determine sequence length");
        return NULL;
    }
    unsigned int d = (unsigned int)size_x;

    double *x;
    if (!(x = (double *)calloc(d, sizeof(double)))) {
        PyErr_SetString(PyExc_MemoryError, "Failed to allocate memory");
        return NULL;
    }

    for (unsigned int i = 0; i < d; ++i) {

        PyObject *item;
        if (!(item = PySequence_GetItem(ob_x, i))) {
            PyErr_SetString(PyExc_TypeError, "Expected a sequence of 'float' objects");
            PY_XDECREF(item);
            free(x);
            return NULL;
        }

        if (!PyFloat_Check(item)) {
            PyErr_SetString(PyExc_TypeError, "Expected a sequence of 'float' objects");
            PY_DECREF(item);
            free(x);
            return NULL;
        }

        *(x + i) = PyFloat_AsDouble(item);

        Py_DECREF(item);

    }

}
