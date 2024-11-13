/**
 * Source file for "../include/numbers.h"
 */

#include <stdlib.h>

#include "../include/numbers.h"


double power(double b, unsigned int p) {

    if (p == 0) { return 1.; }

    if (p % 2 == 0) {
        double temp = power(b, p / 2);
        return temp * temp;
    } else {
        return b * power(b, p - 1);
    }

}

unsigned long ipower(unsigned int b, unsigned int p) {

    unsigned long res = 1L;

    while (p > 0) {

        if (p & 1) { res *= b; }
        res *= res;
        p >>= 1;

    }

    return res;

}

unsigned long factorial_(
    unsigned int n, unsigned long res
) { return (n == 0 || n == 1 ? res : factorial_(n - 1, n * res)); }

unsigned long factorial(unsigned int n) { return factorial_(n, 1); }

long binom(int alpha, unsigned int n) {

    if (n == 0) { return 1; }
    long res = 1L;
    for (int k = 1; k <= n; ++k) { res *= (alpha - k + 1) / k; }
    return res;

}

static PyObject *numbers_power(PyObject *self, PyObject *args) {

    const double b;
    const unsigned int p;
    if (!PyArg_ParseTuple(args, "dI", &b, &p)) { return NULL; }
    return PyFloat_FromDouble(power(b, p));

}

static PyObject *numbers_ipower(PyObject *self, PyObject *args) {

    const unsigned int b;
    const unsigned int p;
    if (!PyArg_ParseTuple(args, "II", &b, &p)) { return NULL; }
    return PyLong_FromUnsignedLong(ipower(b, p));

}

static PyObject *numbers_factorial(PyObject *self, PyObject *args) {

    const unsigned int n;
    if (!PyArg_ParseTuple(args, "I", &n)) { return NULL; }
    return PyLong_FromUnsignedLong(factorial(n));

}

static PyObject *numbers_binom(PyObject *self, PyObject *args) {

    const int alpha;
    const unsigned int k;
    if (!PyArg_ParseTuple(args, "iI")) { return NULL; }
    return PyLong_FromLong(binom(alpha, k));

}
