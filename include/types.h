/**
 * Type definitions for mathematical entities
 */

#define PY_SSIZE_T_CLEAR
#include <Python.h>
#include <stdbool.h>

struct RealFunction {
    bool c;
    union {
        double (*c)(double *x, unsigned int d);
        PyObject *py;
    } f;
};

struct RealFunction *parse_function(PyObject *ob_f);
