/**
 * Source file for "../include/integral.h"
 */

#include <stdlib.h>

#include "../include/integral.h"


double delta(struct Interval **intervals, unsigned int d) {

    double a = 1., b = 1.;
    for (int i = 0; i < d; ++i) {
        a *= (*(intervals + i))->upper - (*(intervals + i))->lower;
        b *= (*(intervals + i))->n;
    }

    return a / b;

}


int endpoint(struct Interval *interval, unsigned int i, double *x) {

    if (!(0 <= i <= interval->n + 1)) { return -1; }

    double dx = (interval->upper - interval->lower) / interval->n;
    *x = interval->lower + i * dx;

    return 0;

}


int left(struct Interval *interval, unsigned int i, double *x) {

    if (!(0 <= i <= interval->n)) { return -1; }

    double dx = (interval->upper - interval->lower) / interval->n;
    *x = interval->lower + i * dx;

    return 0;

}


int right(struct Interval *interval, unsigned int i, double *x) {

    if (!(0 <= i <= interval->n)) { return -1; }

    double dx = (interval->upper - interval->lower) / interval->n;
    *x = interval->lower + (i + 1) * dx;

    return 0;

}


int midpoint(struct Interval *interval, unsigned int i, double *x) {

    if (!(0 <= i <= interval->n)) { return -1; }

    double dx = (interval->upper - interval->lower) / interval->n;
    *x = interval->lower + (i + 1.0 / 2) * dx;

    return 0;

}


static unsigned short int inbounds(struct Interval **intervals, unsigned int radix, unsigned int d) {

    unsigned int bound = 1;
    for (int i = 0; i < d; ++i) { bound *= (*(intervals + i))->n; }

    return 0 <= radix <= bound;

}


static unsigned int increment(struct Interval **intervals, unsigned int radix, unsigned int d) {

    unsigned int base = (*intervals)->n;
    unsigned int msd = radix / base;

    if (d == 1) {
        return (msd >= base - 1 ? 0 : msd + 1);
    }
    return (msd >= base - 1 ? 0 : (msd + 1) * base + increment(intervals + 1, radix % base, d - 1));

}


static unsigned int *unpack(struct Interval **intervals, unsigned int radix, unsigned int d) {

    unsigned int *index;
    if (!(index = (unsigned int *)calloc(d, sizeof(unsigned int)))) { return NULL; }

    for (int i = 0; i < d; ++i) {
        *(index + i) = radix / (*(intervals + i))->n;
        radix %= (*(intervals + i))->n;
    }

    return index;

}


static short int xvalue(
    struct Interval **intervals, RiemannRule *rules, unsigned int radix, unsigned int d,
    double *x
) {

    unsigned int *index;
    if (!(index = unpack(intervals, radix, d))) { return NULL; }

    for (int i = 0; i < d; ++i) {
        if ((*(rules + i))(*(intervals + i), *(index + i), x) == -1) { return -1; }
    }

    return 0;

}


short int riemann_sum(RealFunction f, struct Interval **intervals, RiemannRule *rules, unsigned int d, double *res) {

    double *x;
    if (!(x = (double *)calloc(d, sizeof(double)))) {
        res = NULL;
        return -1;
    }

    const double dv = delta(intervals, d);
    unsigned int radix = 0;

    *res = 0.;

    while (inbounds(intervals, radix, d)) {

        if (xvalue(intervals, rules, radix, d, x) == -1) {
            free(x);
            res = NULL;
            return -1;
        }

        *res += dv * f(x, d);
        radix = increment(intervals, radix, d);

    }

    free(x);

    return 0;

}


short int trapezoidal_rule(RealFunction f, struct Interval **intervals, unsigned int d, double *res) {

    RiemannRule *rules;
    if (!(rules = (RiemannRule *)calloc(d, sizeof(RiemannRule)))) {
        res = NULL;
        return -1;
    }
    for (int i = 0; i < d; ++i) {
        *(rules + i) = endpoint;
    }

    double *x;
    if (!(x = (double *)calloc(d, sizeof(double)))) {
        res = NULL;
        return -1;
    }

    const double dv = delta(intervals, d);
    const unsigned int nlegs = 2 << d;
    unsigned int nborders;
    unsigned int radix = 0;

    *res = 0.;

    while (inbounds(intervals, radix, d)) {

        if (xvalue(intervals, rules, radix, d, x) == -1) {
            free(x);
            res = NULL;
            return -1;
        }

        nborders = 0;
        for (int i = 0; i < d; ++i) {
            if (
                *(x + i) == (*(intervals + i))->lower || *(x + i) == (*(intervals + i))->upper
            ) {
                ++nborders;
            }
        }

        *res += dv * (1 << (d - nborders)) * f(x, d) / nlegs;
        radix = increment(intervals, radix, d);

    }
    
    free(x);

    return 0;

}
