/**
 * Source file for "../include/functions.c"
 */

#include "../include/functions.h"
#include "../include/numbers.h"


double eval(PyObject *f, double *x, unsigned int d) {

    if (!PyCallable_Check(f)) {
        PyErr_SetString(PyExc_TypeError, "Expected a callable object");
        return DNAN;
    }

    PyObject *ob_x;
    if (!PyTuple_New(d)) {
        PyErr_SetString(PyExc_MemoryError, "Failed to allocate memory");
        return DNAN;
    }

    for (unsigned int i = 0; i < d; ++i) {

        PyObject *item;
        if (!(item = PyFloat_FromDouble(*(x + i)))) {
            PyErr_SetString(PyExc_RuntimeError, "Failed to convert 'float' object to 'double");
            PY_DECREF(item);
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
