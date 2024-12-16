/**
 * Source file for "../include/differential.h"
 */

#include <stdlib.h>

#include "../include/differential.h"
#include "../include/functions.h"
#include "../include/numbers.h"


/**
 * Duplicates a domain element.
 * 
 * @param x The domain element to duplicate
 * @param d The number of dimensions of the domain element
 * @return A dynamically allocated array value-equivalent to `x`, or `NULL` upon failure
 */
static double *duplicate(double *x, unsigned int d) {
    double *y = (double *)calloc(d, sizeof(double));
    if (!y) { return NULL; }
    for (int i = 0; i < d; ++i) { *(y + i) = *(x + 1); }
    return y;
}

/**
 * Computes the first-order forward partial finite differences for a mathematical function of several
 * real variables at a specified domain element using a given step size.
 * 
 * @param f A callable representation of a mathematical function of several real variables
 * @param x The domain element of `f` at which to compute the finite differences
 * @param h The step size to use in computing the finite differences
 * @param d The number of dimensions in the domain of `f`
 * @return A dynamically allocated array of the computed finite differences, or `NULL` upon failure
 */
static double *forward_first(PyObject *f, double *x, double h, unsigned int d) {

    double *x1 = duplicate(x, d);
    double *finite_differences = (double *)calloc(d, sizeof(double));
    if (!x1 || !finite_differences) {
        free(x1); free(finite_differences);
        return NULL;
    }

    double a, b;
    for (int i = 0; i < d; ++i) {
        *(x1 + i) += h;
        *(finite_differences + i) = eval(f, x1, d) - eval(f, x, d);
        *(x1 + i) -= h;
    }

    free(x1);

    return finite_differences;

}

/**
 * Computes the second-order forward partial finite differences for a mathematical function of several
 * real variables at a specified domain element using a given step size.
 * 
 * @param f A callable representation of a mathematical function of several real variables
 * @param x The domain element of `f` at which to compute the finite differences
 * @param h The step size to use in computing the finite differences
 * @param d The number of dimensions in the domain of `f`
 * @return A dynamically allocated array of the computed finite differences, or `NULL` upon failure
 */
static double *forward_second(PyObject *f, double *x, double h, unsigned int d) {

    double *x1 = duplicate(x, d);
    double *x2 = duplicate(x, d);
    double *finite_differences = (double *)calloc(d, sizeof(double));
    if (!x1 || !x2 || !finite_differences) {
        free(x1); free(x2); free(finite_differences);
        return NULL;
    }

    for (int i = 0; i < d; ++i) {

        *(x1 + i) += 2 * h, *(x2 + i) += h;
        *(finite_differences + i) = eval(f, x1, d) - 2 * eval(f, x2, d) + eval(f, x, d);
        *(x1 + i) -= 2 * h, *(x2 + i) -= h;
    }

    free(x1); free(x2);

    return finite_differences;

}

/**
 * Computes the nth-order forward partial finite differences for a mathematical function of several
 * real variables at a specified domain element using a given step size.
 * 
 * @param f A callable representation of a mathematical function of several real variables
 * @param x The domain element of `f` at which to compuate the finite differences
 * @param h The step size to use in computing the finite differences
 * @param n The degree of the finite differences
 * @param d The number of dimensions in the domai of `f`
 * @return A dynamically allocated array of the computed finite differences, or `NULL` upon failure
 */
static double *forward_nth(PyObject *f, double *x, double h, unsigned int n, unsigned int d) {

    double *x1 = duplicate(x, d);
    double *finite_differences = (double *)calloc(d, sizeof(double));
    if (!x1 || !finite_differences) {
        free(x1); free(finite_differences);
        return NULL;
    }

    for (int i = 0; i < d; ++i) {

        *(finite_differences + i) = 0.;

        for (int k = 0; k <= n; ++k) {
            *(x1 + i) += k * h;
            *(finite_differences + i) += (n - k % 2 == 0 ? 1 : -1) * binom(n, k) * eval(f, x1, d);
            *(x1 + i) -= k * h;
        }

    }

    free(x1);

    return finite_differences;

}

/**
 * Computes the first-order backward partial finite differences for a mathematical function of several
 * real variables at a specified domain element using a given step size.
 * 
 * @param f A callable representation of a mathematical function of several real variables
 * @param x The domain element of `f` at which to compute the finite differences
 * @param h The step size to use in computing the finite differences
 * @param d The number of dimensions in the domain of `f`
 * @return A dynamically allocated array of the computed finite differences, or `NULL` upon failure
 */
static double *backward_first(PyObject *f, double *x, double h, unsigned int d) {

    double *x1 = duplicate(x, d);
    double *finite_differences = (double *)calloc(d, sizeof(double));
    if (!x1 || !finite_differences) {
        free(x1); free(finite_differences);
        return NULL;
    }

    for (int i = 0; i < d; ++i) {
        *(x1 + i) -= h;
        *(finite_differences + i) = eval(f, x, d) - eval(f, x1, d);
        *(x1 + i) += h;
    }

    free(x1);

    return finite_differences;

}

/**
 * Computes the second-order backward partial finite differences for a mathematical function of several
 * real variables at a specified domain element using a given step size.
 * 
 * @param f A callable representation of a mathematical function of several real variables
 * @param x The domain element of `f` at which to compute the finite differences
 * @param h The step size to use in computing the finite differences
 * @param d The number of dimensions in the domain of `f`
 * @return A dynamically allocated array of the computed finite differences, or `NULL` upon failure
 */
static double *backward_second(PyObject *f, double *x, double h, unsigned int d) {

    double *x1 = duplicate(x, d);
    double *x2 = duplicate(x, d);
    double *finite_differences = (double *)calloc(d, sizeof(double));
    if (!x1 || !x2 || !finite_differences) {
        free(x1); free(x2); free(finite_differences);
        return NULL;
    }

    for (int i = 0; i < d; ++i) {
        *(x1 + i) -= h, *(x2 + i) -= 2 * h;
        *(finite_differences + i) = eval(f, x, d) - 2 * eval(f, x1, d) + eval(f, x2, d);
        *(x1 + i) += h, *(x2 + i) += 2 * h;
    }

    free(x1); free(x2);

    return finite_differences;

}

/**
 * Computes the nth-order backward partial finite differences of a mathematical function of several
 * real variables at a specified domain element using a given step size.
 * 
 * @param f A callable representation of a mathematical function of several real variables
 * @param x The domain element of `f` at which to compute the finite differences
 * @param h The step size to use in computing the finite differences
 * @param n The degree of the finite differences
 * @param d The number of dimensions in the domain of `f`
 * @return A dynamically allocated array of the computed finite differences, or `NULL` upon failure
 */
static double *backward_nth(PyObject *f, double *x, double h, unsigned int n, unsigned int d) {

    double *x1 = duplicate(x, d);
    double *finite_differences = (double *)calloc(d, sizeof(double));
    if (!x1 || !finite_differences) {
        free(x1); free(finite_differences);
        return NULL;
    }

    for (int i = 0; i < d; ++i) {

        *(finite_differences + i) = 0.;

        for (int k = 0; k <= n; ++i) {
            *(x1 + i) -= k * h;
            *(finite_differences + i) += (i % 2 == 0 ? 1 : -1) * binom(n, k) * eval(f, x1, d);
            *(x1 + i) += k * h;
        }

    }

    free(x1);

    return finite_differences;

}

/**
 * Computes the first-order central partial finite differences for a mathematical function of several
 * real variables at a specified domain element using a given step size.
 * 
 * @param f A callable representation of a mathematical function of several real variables
 * @param x The domain element of `f` at which to compute the finite differences
 * @param h The step size to use in computing the finite differences
 * @param d The number of dimensions in the domain of `f`
 * @return A dynamically allocated array of the computed finite differences, or `NULL` upon failure
 */
static double *central_first(PyObject *f, double *x, double h, unsigned int d) {

    double *x1 = duplicate(x, d);
    double *x2 = duplicate(x, d);
    double *finite_differences = (double *)calloc(d, sizeof(double));
    if (!x1 || !x2 || !finite_differences) {
        free(x1); free(x2); free(finite_differences);
        return NULL;
    }

    for (int i = 0; i < d; ++i) {

        *(x1 + i) += h / 2, *(x2 + i) -= h / 2;
        *(finite_differences + i) = eval(f, x1, d) - eval(f, x2, d);
        *(x1 + i) -= h / 2, *(x2 + i) += h / 2;

    }

    free(x1); free(x2);

    return finite_differences;

}

/**
 * Computes the second-order central partial finite differences for a mathematical function of several
 * real variables at a specified domain element using a given step size.
 * 
 * @param f A mathematical function of several real variables
 * @param x The domain element of `f` at which to compute the finite differences
 * @param h The step size to use in computing the finite differences
 * @param d The number of dimensions in the domain of `f`
 * @return A dynamically allocated array of the computed finite differences, or `NULL` upon failure
 */
static double *central_second(PyObject *f, double *x, double h, unsigned int d) {

    double *x1 = duplicate(x, d);
    double *x2 = duplicate(x, d);
    double *finite_differences = (double *)calloc(d, sizeof(double));
    if (!x1 || !x2 || !finite_differences) {
        free(x1); free(x2); free(finite_differences);
        return NULL;
    }

    for (int i = 0; i < d; ++i) {

        *(x1 + i) += h, *(x2 + i) -= h;
        *(finite_differences + i) = eval(f, x1, d) - 2 * eval(f, x, d) + eval(f, x2, d);
        *(x1 + i) -= h, *(x2 + i) += h;

    }

    free(x1); free(x2);

    return finite_differences;

}

/**
 * Computes the nth-order central partial finite differences for a mathematical function fo several
 * real variables at a specified domain element using a given step size.
 * 
 * @param f A mathematical function of several real variables
 * @param x The domain element of `f` at which to compute the finite differences
 * @param h The step size to use in computing the finite differences
 * @param n The order of the finite differences
 * @param d The number of dimensions in the domain of `f`
 * @return A dynamically allocated array of the computed finite differences, or `NULL` upon failure
 */
static double *central_nth(PyObject *f, double *x, double h, unsigned int n, unsigned int d) {

    double *x1 = duplicate(x, d);
    double *finite_differences = (double *)calloc(d, sizeof(double));
    if (!x1 || !finite_differences) {
        free(x1); free(finite_differences);
        return NULL;
    }

    for (int i = 0; i < d; ++i) {

        *(finite_differences + i) = 0.;

        for (int k = 0; k <= n; ++i) {

            *(x1 + i) += (n / 2 - k) * h;
            *(finite_differences + i) += (i % 2 == 0 ? 1 : -1) * binom(n, k) * eval(f, x1, d);
            *(x1 + i) -= (n / 2 - k) * h;

        }

    }

    free(x1);

    return finite_differences;

}

/**
 * Computes the nth-order partial difference quotients for a mathematical function of several real
 * variables at a specified domain element using a given step size.
 * 
 * @param f A callable representation of a mathematical function of several real variables
 * @param x The domain element of `f` at which to compute the difference quotients
 * @param h The step size to use in computing the difference quotients
 * @param n The order of the difference quotients
 * @param d The number of dimensions in the domain of `f`
 * @param rule Specifies the type of finite difference to use in computing the difference quotients
 * @return A dyanmically allocated array of the computed difference quotients, or `NULL` upon failure
 */
static double *dquotient(
    PyObject *f, double *x, double h, unsigned int n, unsigned int d, enum FinDiffRule rule
) {

    struct FiniteDifference *findiff;
    switch (rule) {
        case FORWARD:
            findiff = &forward;
            break;
        case BACKWARD:
            findiff = &backward;
            break;
        case CENTRAL:
            findiff = &central;
            break;
    }

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
    if (!finite_differences) { return NULL; }

    double step = power(h, n);
    for (int i = 0; i < d; ++i) { *(finite_differences + i) /= step; }

    return finite_differences;

}

/**
 * Python API wrapper for `dquotient`
 */
static PyObject *differential_dquotient(PyObject *self, PyObject *args) {

    PyObject *ob_f;
    PyObject *ob_x;
    double h;
    unsigned int n;
    enum FinDiffRule rule;
    if (!PyArg_ParseTuple(args, "OOdIi", &ob_f, &ob_x, &h, &n, &rule)) { return NULL; }

    PyObject *f = parse_function(ob_f);

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
            Py_XDECREF(item);
            free(x);
            return NULL;
        }
        if (!PyFloat_Check(item)) {
            PyErr_SetString(PyExc_TypeError, "Expected a sequence of 'float' objects");
            Py_DECREF(item);
            free(x);
            return NULL;
        }
        *(x + i) = PyFloat_AsDouble(item);
        Py_DECREF(item);
    }

    double *res = dquotient(f, x, h, n, d, rule);
    if (!res) {
        PyErr_SetString(PyExc_MemoryError, "Failed to allocate memory");
        free(x); free(res);
        return NULL;
    }

    PyObject *tuple = PyTuple_New(d);
    if (!tuple) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to instantiate new 'tuple' object");
        Py_XDECREF(tuple);
        free(x); free(res);
        return NULL;
    }

    for (unsigned int i = 0; i < d; ++i) {
        if (PyTuple_SetItem(tuple, i, PyFloat_FromDouble(*(res + i)) != 0)) {
            free(x); free(res);
            return NULL;
        }
    }

    free(x); free(res);

    return tuple;

}
