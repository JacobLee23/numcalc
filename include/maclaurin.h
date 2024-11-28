/**
 * Maclaurin series computation
 */

#define PY_SSIZE_T_CLEAN
#include <Python.h>

double exponential(double x);
double ln(double x);
double geometric(double x, unsigned int alpha);
double binomial(double x, unsigned int alpha);
double root(double x, unsigned int alpha);
double invroot(double x, unsigned int alpha);

double sine(double x);
double cosine(double x);
double tangent(double x);
double secant(double x);
double cosecant(double x);
double cotangent(double x);

double arcsine(double x);
double arccosine(double x);
double arctangent(double x);
double arcsecant(double x);
double arccosecant(double x);
double arccotangent(double x);

double sineh(double x);
double cosineh(double x);
double tangenth(double x);
double secanth(double x);
double cosecanth(double x);
double cotangenth(double x);

double arcsineh(double x);
double arccosh(double x);
double arctangenth(double x);
double arcsecanth(double x);
double arccosecanth(double x);
double arccotangenth(double x);

static PyObject *maclaurin_exp(PyObject *self, PyObject *args);
static PyObject *maclaurin_ln(PyObject *self, PyObject *args);
static PyObject *maclaurin_geometric(PyObject *self, PyObject *args);
static PyObject *maclaurin_binomial(PyObject *self, PyObject *args);
static PyObject *maclaurin_root(PyObject *self, PyObject *args);
static PyObject *maclaurin_invroot(PyObject *self, PyObject *args);

static PyObject *maclaurin_sin(PyObject *self, PyObject *args);
static PyObject *maclaurin_cos(PyObject *self, PyObject *args);
static PyObject *maclaurin_tan(PyObject *self, PyObject *args);
static PyObject *maclaurin_sec(PyObject *self, PyObject *args);
static PyObject *maclaurin_csc(PyObject *self, PyObject *args);
static PyObject *maclaurin_cot(PyObject *self, PyObject *args);

static PyObject *maclaurin_arcsin(PyObject *self, PyObject *args);
static PyObject *maclaurin_arccos(PyObject *self, PyObject *args);
static PyObject *maclaurin_arctan(PyObject *self, PyObject *args);
static PyObject *maclaurin_arcsec(PyObject *self, PyObject *args);
static PyObject *maclaurin_arccsc(PyObject *self, PyObject *args);
static PyObject *maclaurin_arccot(PyObject *self, PyObject *args);

static PyObject *maclaurin_sinh(PyObject *self, PyObject *args);
static PyObject *maclaurin_cosh(PyObject *self, PyObject *args);
static PyObject *maclaurin_tanh(PyObject *self, PyObject *args);
static PyObject *maclaurin_sech(PyObject *self, PyObject *args);
static PyObject *maclaurin_csch(PyObject *self, PyObject *args);
static PyObject *maclaurin_coth(PyObject *self, PyObject *args);

static PyObject *maclaurin_arcsinh(PyObject *self, PyObject *args);
static PyObject *maclaurin_arccosh(PyObject *self, PyObject *args);
static PyObject *maclaurin_arctanh(PyObject *self, PyObject *args);
static PyObject *maclaurin_arcsech(PyObject *self, PyObject *args);
static PyObject *maclaurin_arccsch(PyObject *self, PyObject *args);
static PyObject *maclaurin_arccoth(PyObject *self, PyObject *args);

static PyMethodDef MaclaurinMethods[] = {
    {"exp", maclaurin_exp, METH_VARARGS, NULL},
    {"ln", maclaurin_ln, METH_VARARGS, NULL},
    {"geometric", maclaurin_geometric, METH_VARARGS, NULL},
    {"binomial", maclaurin_binomial, METH_VARARGS, NULL},
    {"root", maclaurin_root, METH_VARARGS, NULL},
    {"invroot", maclaurin_invroot, METH_VARARGS, NULL},
    {"sin", maclaurin_sin, METH_VARARGS, NULL},
    {"cos", maclaurin_cos, METH_VARARGS, NULL},
    {"tan", maclaurin_tan, METH_VARARGS, NULL},
    {"sec", maclaurin_sec, METH_VARARGS, NULL},
    {"csc", maclaurin_csc, METH_VARARGS, NULL},
    {"cot", maclaurin_cot, METH_VARARGS, NULL},
    {"arcsin", maclaurin_arcsin, METH_VARARGS, NULL},
    {"arccos", maclaurin_arccos, METH_VARARGS, NULL},
    {"arctan", maclaurin_arctan, METH_VARARGS, NULL},
    {"arcsec", maclaurin_arcsec, METH_VARARGS, NULL},
    {"arccsc", maclaurin_arccsc, METH_VARARGS, NULL},
    {"arccot", maclaurin_arccot, METH_VARARGS, NULL},
    {"sinh", maclaurin_sinh, METH_VARARGS, NULL},
    {"cosh", maclaurin_cosh, METH_VARARGS, NULL},
    {"tanh", maclaurin_tanh, METH_VARARGS, NULL},
    {"sech", maclaurin_sech, METH_VARARGS, NULL},
    {"csch", maclaurin_csch, METH_VARARGS, NULL},
    {"coth", maclaurin_coth, METH_VARARGS, NULL},
    {"arcsinh", maclaurin_arcsinh, METH_VARARGS, NULL},
    {"arccosh", maclaurin_arccosh, METH_VARARGS, NULL},
    {"arctanh", maclaurin_arctanh, METH_VARARGS, NULL},
    {"arcsech", maclaurin_arcsech, METH_VARARGS, NULL},
    {"arccsch", maclaurin_arccsch, METH_VARARGS, NULL},
    {"arccoth", maclaurin_arccoth, METH_VARARGS, NULL},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef maclaurin_module = {
    PyModuleDef_HEAD_INIT, "maclaurin", NULL, -1, MaclaurinMethods
};

PyMODINIT_FUNC PyInit_maclaurin() { return PyModule_Create(&maclaurin_module); }
