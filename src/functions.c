/**
 * Source file for "../include/functions.c"
 */

#include "../include/functions.h"
#include "../include/numbers.h"


/**
 * Evaluates a mathematical function of several real variables at a specified domain element.
 *
 * @param f A callable representation of a mathematical function of several real variables
 * @param x The domain element at which to evaluate `f`
 * @param d The number of dimensions in the domain of `f`
 * @return The value of `f` at `x`, or `DNAN` upon failure
 */
double eval(PyObject *f, double *x, unsigned int d) {

    if (!PyCallable_Check(f)) {
        PyErr_SetString(PyExc_TypeError, "Expected a callable object");
        return DNAN;
    }

    PyObject *ob_x;
    if (!PyTuple_New(d)) {
        PyErr_SetString(PyExc_MemoryError, "Failed to instantiate 'tuple' object");
        return DNAN;
    }

    for (unsigned int i = 0; i < d; ++i) {
        PyObject *item = PyFloat_FromDouble(*(x + i));
        if (!item) {
            PyErr_SetString(PyExc_RuntimeError, "Failed to convert 'float' object to double");
            Py_DECREF(item);
            return DNAN;
        }
        PyTuple_SetItem(ob_x, (Py_ssize_t)i, item);
    }

    PyObject *res = PyObject_CallOneArg(f, ob_x);
    Py_DECREF(ob_x);
    if (!res) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to call callable object");
        return DNAN;
    }
    if (!PyFloat_Check(res)) {
        PyErr_SetString(PyExc_TypeError, "Expected callable object to return a 'float' object");
        Py_DECREF(res);
        return DNAN;
    }
    double value = PyFloat_AsDouble(res);
    Py_DECREF(res);

    return value;

}
