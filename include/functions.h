/**
 * Evaluation of mathematical function representations
 */

#define PY_SSIZE_T_CLEAN
#include <Python.h>


double eval(PyObject *f, double *x, unsigned int d);
