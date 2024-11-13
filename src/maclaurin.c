/**
 * Source file for "../include/maclaurin.h"
 */

#include "../include/maclaurin.h"
#include "../include/numbers.h"


static double exponential_(
    double x, unsigned int n, double term
) { return (n == 0 ? 1 : x / n * term); }

double exponential(double x) {

    unsigned int n = 0;
    double term, res = 0.;
    while ((term = exponential_(x, n++, term)) != 0) { res += term; }
    return res;

}

static double ln_(double x, unsigned int n, double term) {

    if (!(-1 < x <= 1)) { return DNAN; }
    return (n == 1 ? x : -x * (n - 1) / n * term);

}

double ln(double x) {

    if (x <= 0.) { return DNAN; }
    unsigned int n = 1;
    double term, res = 0.;
    while ((term = ln_((x <= 2 ? x - 1 : (1 - x) / x), n++, term)) != 0) { res += term; }
    return (x <= 2 ? res : -res);

}

static double geometric_(double x, unsigned int alpha, unsigned int n, double term) {

    if (!(-1 < x < 1)) { return DNAN; }
    return (n == (alpha - 1) ? alpha - 1 : x * n / (n - alpha + 1) * term);

}

double geometric(double x, unsigned int alpha) {

    if (x == 0.) { return DNAN; }
    unsigned int n = alpha - 1;
    double term, res = 0.;
    while ((term = geometric_(x + 1, alpha, n++, term)) != 0) { res += term; }
    return res;

}

static double binomial_(double x, unsigned int alpha, unsigned int n, double term) {

    if (!(-1 < x < 1)) { return DNAN; }
    return (n == 0 ? alpha * x : x * (alpha - n + 1) / n * term);

}

double binomial(double x, unsigned int alpha) {

    unsigned int n = 0;
    double term, res = 0.;
    while ((term = binomial_(x - 1, alpha, n++, term)) != 0) { res += term; }
    return res;

}

static double root_(
    double x, unsigned int n, unsigned int alpha, double term
) {

    if (!(-1 < x < 1)) { return DNAN; }
    return (n == 0 ? 1 : x * (1 - (n - 1) * alpha) / (n * alpha) * term);

}

double root(double x, unsigned int alpha) {

    if (x < 0.) { return DNAN; }
    unsigned int n = 0;
    double term, res = 0.;
    while ((term = root_(x - 1, alpha, n++, term)) != 0) { res += term; }
    return res;

}

static double invroot_(
    double x, unsigned int n, unsigned int alpha, double term
) {

    if (!(-1 < x < 1)) { return DNAN; }
    return (n == 0 ? 1 : x * (-1 - (n - 1) * alpha) / (n * alpha) * term);

}

double invroot(double x, unsigned int alpha) {

    unsigned int n = 0;
    double term, res = 0.;
    while ((term = invroot_(x - 1, alpha, n++, term)) != 0) { res += term; }
    return res;

}

static double sine_(
    double x, unsigned int n, double term
) { return (n == 0 ? x : -power(x, 2) / ((2 * n) * (2 * n + 1)) * term); }

double sine(double x) {

    unsigned int n = 0;
    double term, res = 0.;
    while ((term = sine_(x, n++, term)) != 0) { res += term; }
    return res;

}

static double cosine_(
    double x, unsigned int n, double term
) { return (n == 0 ? 1 : -power(x, 2) / ((2 * n) * (2 * n - 1)) * term); }

double cosine(double x) {

    unsigned int n = 0;
    double term, res = 0.;
    while ((term = cosine_(x, n++, term)) != 0) { res += term; }
    return res;

}

double tangent(double x) { return sine(x) / cosine(x); }
double secant(double x) { return 1 / cosine(x); }
double cosecant(double x) { return 1 / sine(x); }
double cotangent(double x) { return sine(x) / cosine(x); }

static double arcsine_(double x, unsigned int n, double term) {

    if (!(-1 <= x <= 1)) { return DNAN; }
    return (
        n == 0 ? x : power(x, 2) * (
            ((2 * n) * ipower(2 * n - 1, 2)) / (4 * ipower(n, 2) * (2 * n + 1))
        ) * term
    );

}

double arcsine(double x) {

    if (!(-1 <= x <= 1)) { return DNAN; }

    unsigned int n = 0;
    double term, res = 0.;
    while ((term = arcsine_(x, n++, term)) != 0) { res += term; }
    return res;

}

double arccosine(double x) { return arcsine(1) - arcsine(x); }

static double arctangent_(double x, unsigned int n, double term) {

    if (!(-1 <= x <= 1)) { return DNAN; }
    return (n == 0 ? x : -power(x, 2) * (2 * n - 1) / (2 * n + 1) * term);

}

double arctangent(double x) {

    unsigned int n = 0;
    double term, res = 0.;
    while ((term = arctangent_(x, n++, term)) != 0) { res += term; }

}

double arcsecant(double x) { return arccosine(1 / x); }
double arccosecant(double x) { return arcsine(1 / x); }
double arccotangent(double x) { return arctangent(1 / x); }

static double sineh_(
    double x, unsigned int n, double term
) { return (n == 0 ? x: power(x, 2) / ((2 * n) * (2 * n + 1)) * term); }

double sineh(double x) {

    unsigned int n = 0;
    double term, res = 0.;
    while ((term = sineh_(x, n++, term)) != 0) { res += term; }
    return res;

}

static double cosineh_(
    double x, unsigned int n, double term
) { return (n == 0 ? 1 : power(x, 2) / ((2 * n) * (2 * n - 1)) * term); }

double cosineh(double x) {

    unsigned int n = 0;
    double term, res = 0.;
    while ((term = cosineh_(x, n++, term)) != 0) { res += term; }
    return res;

}

double tangenth(double x) { return sineh(x) / cosineh(x); }
double secanth(double x) { return 1 / cosineh(x); }
double cosecanth(double x) { return 1 / secanth(x); }
double cotangenth(double x) { return cosineh(x) / sineh(x); }

static double arcsineh_(double x, unsigned int n, double term) {

    if (!(-1 <= x <= 1)) { return DNAN; }
    return (
        n == 0 ? x : -power(x, 2) * (
            (2 * n) * ipower(2 * n - 1, 2)) / (4 * ipower(n, 2) * (2 * n + 1)
        ) * term
    );

}

double arcsineh(double x) {

    if (!(-1 <= x <= 1)) { return DNAN; }
    unsigned int n = 0;
    double term, res = 0.;
    while ((term = arcsineh_(x, n++, term)) != 0) { res += term; }
    return res;

}

double arccosineh(double x) {

    double res = arcsinh(root(power(x, 2), 2) - 1);
    return (res >= 0. ? res : -res);

}

static double arctangenth_(double x, unsigned int n, double term) {

    if (!(-1 < x < 1)) { return DNAN; }
    return (n == 0 ? x : power(x, 2) * (2 * n - 1) / (2 * n + 1) * term);

}

double arctangenth(double x) {

    unsigned int n = 0;
    double term, res = 0.;
    while ((term = arctangenth_(x, n++, term)) != 0) { res += term; }
    return res;

}

double arcsecanth(double x) { return arccosineh(1 / x); }
double arccosecanth(double x) { return arccsecanth(1 / x); }
double arccotangenth(double x) { return arctangenth(1 / x); }

static PyObject *maclaurin_(double (*func)(double), PyObject *args) {

    static double x;
    if (!PyArg_ParseTuple(args, "d", &x)) { return NULL; }
    return PyFloat_FromDouble(func(x));

}

static PyObject *maclaurina_(double (*func)(double, unsigned int), PyObject *args) {

    static double x;
    static unsigned int alpha;
    if (!PyArg_ParseTuple(args, "dI", &x, &alpha)) { return NULL; }
    return PyFloat_FromDouble(func(x, alpha));

}

static PyObject *maclaurin_exp(PyObject *self, PyObject *args) { return maclaurin_(exp, args); }
static PyObject *maclaurin_ln(PyObject *self, PyObject *args) { return maclaurin_(ln, args); }
static PyObject *maclaurin_geometric(PyObject *self, PyObject *args) { return maclaurina_(geometric, args); }
static PyObject *maclaurin_binomial(PyObject *self, PyObject *args) { return maclaurina_(binomial, args); }
static PyObject *maclaurin_root(PyObject *self, PyObject *args) { return maclaurina_(root, args); }
static PyObject *maclaurin_invroot(PyObject *self, PyObject *args) { return maclaurina_(invroot, args); }

static PyObject *maclaurin_sin(PyObject *self, PyObject *args) { return maclaurin_(sine, args); }
static PyObject *maclaurin_cos(PyObject *self, PyObject *args) { return maclaurin_(cosine, args); }
static PyObject *maclaurin_tan(PyObject *self, PyObject *args) { return maclaurin_(tangent, args); }
static PyObject *maclaurin_sec(PyObject *self, PyObject *args) { return maclaurin_(secant, args); }
static PyObject *maclaurin_csc(PyObject *self, PyObject *args) { return maclaurin_(cosecant, args); }
static PyObject *maclaurin_cot(PyObject *self, PyObject *args) { return maclaurin_(cotangent, args); }

static PyObject *maclaurin_arcsin(PyObject *self, PyObject *args) { return maclaurin_(arcsine, args); }
static PyObject *maclaurin_arccos(PyObject *self, PyObject *args) { return maclaurin_(arccosine, args); }
static PyObject *maclaurin_arctan(PyObject *self, PyObject *args) { return maclaurin_(arctangent, args); }
static PyObject *maclaurin_arcsec(PyObject *self, PyObject *args) { return maclaurin_(arcsecant, args); }
static PyObject *maclaurin_arccsc(PyObject *self, PyObject *args) { return maclaurin_(arccosecant, args); }
static PyObject *maclaurin_arccot(PyObject *self, PyObject *args) { return maclaurin_(arccotangent, args); }

static PyObject *maclaurin_sinh(PyObject *self, PyObject *args) { return maclaurin_(sineh, args); }
static PyObject *maclaurin_cosh(PyObject *self, PyObject *args) { return maclaurin_(cosineh, args); }
static PyObject *maclaurin_tanh(PyObject *self, PyObject *args) { return maclaurin_(tangenth, args); }
static PyObject *maclaurin_sech(PyObject *self, PyObject *args) { return maclaurin_(secanth, args); }
static PyObject *maclaurin_csch(PyObject *self, PyObject *args) { return maclaurin_(cosecanth, args); }
static PyObject *maclaurin_coth(PyObject *self, PyObject *args) { return maclaurin_(cotangenth, args); }

static PyObject *maclaurin_arcsinh(PyObject *self, PyObject *args) { return maclaurin_(arcsineh, args); }
static PyObject *maclaurin_arccosh(PyObject *self, PyObject *args) { return maclaurin_(arccosineh, args); }
static PyObject *maclaurin_arctanh(PyObject *self, PyObject *args) { return maclaurin_(arctangenth, args); }
static PyObject *maclaurin_arcsech(PyObject *self, PyObject *args) { return maclaurin_(arcsecanth, args); }
static PyObject *maclaurin_arccsch(PyObject *self, PyObject *args) { return maclaurin_(arccosecanth, args); }
static PyObject *maclaurin_arccoth(PyObject *self, PyObject *args) { return maclaurin_(arccotangenth, args); }
