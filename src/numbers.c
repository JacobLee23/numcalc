/**
 * Source file for "../include/numbers.h"
 */

#include <stdlib.h>

#include "../include/numbers.h"


/**
 * Evaluates a real base `b` raised to the power of a non-negative integer `p`.
 */
static double power(double b, unsigned int p) {

    if (p == 0) { return 1.; }

    if (p % 2 == 0) {
        double temp = power(b, p / 2);
        return temp * temp;
    } else {
        return b * power(b, p - 1);
    }

}

/**
 * Evaluates a non-negative base `b` raised to the power of a non-negative integer `p`.
 */
static unsigned long ipower(unsigned int b, unsigned int p) {

    unsigned long res = 1L;

    while (p > 0) {

        if (p & 1) { res *= b; }
        res *= res;
        p >>= 1;

    }

    return res;

}

/**
 * Recursively evaluates the factorial of a non-negative integer `n`.
 */
static unsigned long factorialr(
    unsigned int n, unsigned long res
) { return (n == 0 || n == 1 ? res : factorialr(n - 1, n * res)); }

/**
 * Evaluates the factorial of a non-negative integer `n`.
 */
unsigned long factorial(unsigned int n) { return factorialr(n, 1); }

/**
 * Computes the `n`th binomial number of integer `alpha`.
 */
long binom(int alpha, unsigned int n) {

    if (n == 0) { return 1; }
    long res = 1L;
    for (int k = 1; k <= n; ++k) { res *= (alpha - k + 1) / k; }
    return res;

}

/**
 * Python wrapper for `power`.
 */
static PyObject *numbers_power(PyObject *self, PyObject *args) {

    const double b;
    const unsigned int p;
    if (!PyArg_ParseTuple(args, "dI", &b, &p)) { return NULL; }
    return PyFloat_FromDouble(power(b, p));

}

/**
 * Python wrapper for `ipower`.
 */
static PyObject *numbers_ipower(PyObject *self, PyObject *args) {

    const unsigned int b;
    const unsigned int p;
    if (!PyArg_ParseTuple(args, "II", &b, &p)) { return NULL; }
    return PyLong_FromUnsignedLong(ipower(b, p));

}

/**
 * Python wrapper for `factorial`.
 */
static PyObject *numbers_factorial(PyObject *self, PyObject *args) {

    const unsigned int n;
    if (!PyArg_ParseTuple(args, "I", &n)) { return NULL; }
    return PyLong_FromUnsignedLong(factorial(n));

}

/**
 * Python wrapper for `binom`.
 */
static PyObject *numbers_binom(PyObject *self, PyObject *args) {

    const int alpha;
    const unsigned int k;
    if (!PyArg_ParseTuple(args, "iI")) { return NULL; }
    return PyLong_FromLong(binom(alpha, k));

}
