/**
 * Source file for "../include/maclaurin.h"
 */

#include "../include/maclaurin.h"
#include "../include/numbers.h"


int exp_(double x, unsigned int n, double *term) {

    *term = n == 0 ? 1 : x / n * *term;
    return 0;

}


double exp(double x) {
    
    unsigned int n = 0;
    double term, res = 0.;
    while (!exp_(x, n++, &term) && term != 0) { res += term; }

    return res;

}


int ln_(double x, unsigned int n, double *term) {

    if (x < -1) {
        *term = n == 1 ? (1. / x + 1) : (1. / x + 1) * (n - 1) / n * *term;
    } else if (-2 <= x < 0) {
        *term = n == 1 ? -(x + 1) : (x + 1) * (n - 1) / n * *term;
    } else if (0 < x <= 2) {
        *term = n == 1 ? (x - 1) : -(x - 1) * (n - 1) / n * *term;
    } else if (x > 2) {
        *term = n == 1 ? -(1. / x - 1) : -(1. / x - 1) * (n - 1) / n * *term;
    } else {
        term = 0;
        return -1;
    }

    return 0;

}


double ln(double x) {

    unsigned int n = 1;
    double term, res = 0.;
    while (!ln_(x, n++, &term) && term != 0) { res += term; }

    return res;

}


int geometric_(double x, unsigned int alpha, unsigned int n, double *term) {

    if (alpha <= 0) {
        term = 0;
        return -1;
    }

    *term = n == (alpha - 1) ? alpha - 1 : n / (n - alpha + 1) * x * *term;
    return 0;

}


double geometric(double x, unsigned int alpha) {

    unsigned int n = alpha - 1;
    double term, res = 0.;
    while (!geometric_(x, alpha, n++, &term) && term != 0) { res += term; }

    return res;

}


int binomial_(double x, unsigned int alpha, unsigned int n, double *term) {

    *term = n == 0 ? alpha * x : (alpha - n + 1) / n * *term;
    return 0;

}


double binomial(double x, unsigned int alpha) {

    unsigned int n = 0;
    double term, res = 0.;
    while (!binomial_(x, alpha, n++, &term) && term != 0) { res += term; }

    return res;

}


int sqrt_(double x, unsigned int n, double *term) {

    *term = n == 0 ? 1 : -x * (2 * n - 3) / (2 * n) * *term;
    return 0;

}


double sqrt(double x) {

    unsigned int n = 0;
    double term, res = 0.;
    while (!sqrt_(x, n++, &term) && term != 0) { res += term; }

    return res;

}


int invsqrt_(double x, unsigned int n, double *term) {

    *term = n == 0 ? 1 : -x * (2 * n - 1) / (2 * n) * *term;
    return 0;

}


double invsqrt(double x) {

    unsigned int n = 0;
    double term, res = 0.;
    while (!invsqrt_(x, n++, &term) && term != 0) { res += term; }

    return res;

}


int sin_(double x, unsigned int n, double *term) {

    *term = n == 0 ? x : -(x * x) / ((2 * n) * (2 * n + 1)) * *term;
    return 0;

}


double sin(double x) {

    unsigned int n = 0;
    double term, res = 0.;
    while (!sin_(x, n++, &term) && term != 0) { res += term; }

    return res;

}


int cos_(double x, unsigned int n, double *term) {

    *term = n == 0 ? 1 : -(x * x) / ((2 * n) * (2 * n - 1)) * *term;
    return 0;

}


double cos(double x) {

    unsigned int n = 0;
    double term, res = 0.;
    while (!cos_(x, n++, &term) && term != 0) { res += term; }

    return res;

}


double tan(double x) { return sin(x) / cos(x); }


double sec(double x) { return 1 / cos(x); }


double csc(double x) { return 1 / sin(x); }


double cot(double x) { return cos(x) / sin(x); }


int arcsin_(double x, unsigned int n, double *term) {

    if (!(-1 <= x <= 1)) {
        term = 0;
        return -1;
    }

    *term = n == 0 ? x : ((2 * n - 1) * (2 * n) * (2 * n - 1) * (x * x)) / (4 * (n * n) * (2 * n + 1));
    return 0;

}


double arcsin(double x) {

    unsigned int n = 0;
    double term, res = 0.;
    while (!arcsin_(x, n++, &term) && term != 0) { res += term; }

    return res;

}


double arccos(double x) { return arcsin(1) - arcsin(x); }


int arctan_(double x, unsigned int n, double *term) {

    if (!(-1 <= x <= 1)) {
        term = 0;
        return -1;
    }

    *term = n == 0 ? x : -((2 * n - 1) * (x * x)) / (2 * n + 1);
    return 0;

}


double arctan(double x) {

    unsigned int n = 0;
    double term, res = 0.;
    while (!arctan_(x, n++, &term) && term != 0) { res += term; }

    return res;

}


double arcsec(double x) { return arccos(1 / x); }


double arccsc(double x) { return arcsin(1 / x); }


double arccot(double x) { return arctan(1 / x); }


int sinh_(double x, unsigned int n, double *term) {

    *term = n == 0 ? x : (x * x) / ((2 * n) * (2 * n + 1));
    return 0;

}


double sinh(double x) {

    unsigned int n = 0;
    double term, res = 0.;
    while (!sinh_(x, n++, &term) && term != 0) { res += term; }

    return res;

}


int cosh_(double x, unsigned int n, double *term) {

    *term = n == 0 ? 1 : (x * x) / ((2 * n - 1) * (2 * n));
    return 0;

}


double cosh(double x) {

    unsigned int n = 0;
    double term, res = 0.;
    while (!cosh_(x, n++, &term) && term != 0) { res += term; }

    return res;

}


double tanh(double x) { return sinh(x) / cosh(x); }


double sech(double x) { return 1 / cosh(x); }


double csch(double x) { return 1 / sech(x); }


double coth(double x) { return cosh(x) / sinh(x); }


int arcsinh_(double x, unsigned int n, double *term) {

    if (!(-1 <= x <= 1)) {
        term = 0;
        return -1;
    }

    *term = n == 0 ? x : -((2 * n - 1) * (2 * n) * (2 * n - 1) * (x * x)) / (4 * (n * n) * (2 * n + 1));
    return 0;

}


double arcsinh(double x) {

    unsigned int n = 0;
    double term, res = 0.;
    while (!arcsinh_(x, n++, &term) && term != 0) { res += term; }

    return res;

}


int arccosh_(double x, unsigned int n, double *term) {

    if (!(x >= 1)) {
        term = 0;
        return -1;
    }

    *term = n == 1 ? 1 / (4 * (x * x)) : ((2 * n - 1) * (2 * n - 2)) / (4 * (n * n) * (x * x)) * *term;
    return 0;

}


double arccosh(double x) {

    unsigned int n = 1;
    double term, res = 0.;
    while (!arccosh_(x, n++, &term) && term != 0) { res += term; }

    return ln(2 * x) - res;

}


int arctanh_(double x, unsigned int n, double *term) {

    if (!(-1 < x < 1)) {
        term = 0;
        return -1;
    }

    *term = n == 0 ? x : ((2 * n - 1) * (x * x)) / (2 * n + 1);
    return 0;

}


double arctanh(double x) {

    unsigned int n = 0;
    double term, res = 0.;
    while (!arctanh_(x, n++, &term) && term != 0) { res += term; }

    return res;

}


double arcsech(double x) { return arccosh(1 / x); }


double arccsch(double x) { return arccsch(1 / x); }


double arccoth(double x) { return arctanh(1 / x); }
