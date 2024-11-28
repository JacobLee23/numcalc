/**
 * Auxiliary mathematical computation
 */

#define PY_SSIZE_T_CLEAN
#include <Python.h>

#define DNAN 0. / 0.


double power(double b, unsigned int p);
unsigned long ipower(unsigned int b, unsigned int p);
unsigned long factorial(unsigned int n);
long binom(int alpha, unsigned int k);

static PyObject *numbers_power(PyObject *self, PyObject *args);
static PyObject *numbers_ipower(PyObject *self, PyObject *args);
static PyObject *numbers_factorial(PyObject *self, PyObject *args);
static PyObject *numbers_binom(PyObject *self, PyObject *args);

static PyMethodDef NumbersMethods[] = {
    {"power", numbers_power, METH_VARARGS, NULL},
    {"intpower", numbers_ipower, METH_VARARGS, NULL},
    {"factorial", numbers_factorial, METH_VARARGS, NULL},
    {"binom", numbers_binom, METH_VARARGS, NULL},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef numbers_module = {
    PyModuleDef_HEAD_INIT, "numbers", NULL, -1, NumbersMethods
};

PyMODINIT_FUNC PyInit_numbers() { return PyModule_Create(&numbers_module); }
