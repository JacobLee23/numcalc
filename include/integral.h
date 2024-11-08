/**
 * Numerical integral calculus
 */

#ifndef TYPES_H
#define TYPES_H
#include "types.h"
#endif


struct Interval {
    double lower;
    double upper;
    unsigned int n;
};

double delta(struct Interval **intervals, unsigned int d);

int endpoint(struct Interval *interval, unsigned int i, double *x);
int left(struct Interval *interval, unsigned int i, double *x);
int right(struct Interval *interval, unsigned int i, double *x);
int midpoint(struct Interval *interval, unsigned int i, double *x);
typedef int (*RiemannRule)(struct Interval *interval, unsigned int i, double *x);

unsigned short int inbounds(struct Interval **intervals, unsigned int radix, unsigned int d);
unsigned int increment(struct Interval **intervals, unsigned int radix, unsigned int d);
unsigned int *decode(struct Interval **intervals, unsigned int radix, unsigned int d);
short int xvalue(
    struct Interval **intervals, RiemannRule *rules, unsigned int radix, unsigned int d, double *x
);

short int riemann_sum(
    RealFunction f, struct Interval **intervals, RiemannRule *rules, unsigned int d, double *res
);
short int trapezoidal_rule(
    RealFunction f, struct Interval **intervals, unsigned int d, double *res
);
