/**
 * Source file for "../include/differential.h"
 */

#include <stdlib.h>

#include "../include/differential.h"
#include "../include/numbers.h"


static double *copyx(double *x, unsigned int d) {

    double *x_;
    if (!(x_ = (double *)calloc(d, sizeof(double)))) { return NULL; }
    for (int i = 0; i < d; ++i) { *(x_ + i) = *(x + 1); }
    return x_;

}

double *forward_first(RealFunction f, double *x, double h, unsigned int d) {

    double *x1 = copyx(x, d);
    double *finite_differences;
    if (!(finite_differences = (double *)calloc(d, sizeof(double)))) { return NULL; }

    double a, b;
    for (int i = 0; i < d; ++i) {

        *(x1 + i) += h;
        *(finite_differences + i) = f(x1, d) - f(x, d);
        *(x1 + i) -= h;

    }
    
    free(x1);

    return finite_differences;

}

double *forward_second(RealFunction f, double *x, double h, unsigned int d) {

    double *x1 = copyx(x, d);
    double *x2 = copyx(x, d);
    double *finite_differences;
    if (!(finite_differences = (double *)calloc(d, sizeof(double)))) { return NULL; }

    for (int i = 0; i < d; ++i) {

        *(x1 + i) += 2 * h, *(x2 + i) += h;
        *(finite_differences + i) = f(x1, d) - 2 * f(x2, d) + f(x, d);
        *(x1 + i) -= 2 * h, *(x2 + i) -= h;
    }

    free(x1);
    free(x2);

    return finite_differences;

}

double *forward_nth(RealFunction f, double *x, double h, unsigned int d, unsigned int n) {

    double *x1 = copyx(x, d);
    double *finite_differences;
    if (!(finite_differences = (double *)calloc(d, sizeof(double)))) { return NULL; }

    for (int i = 0; i < d; ++i) {

        *(finite_differences + i) = 0.;

        for (int k = 0; k <= n; ++k) {

            *(x1 + i) += k * h;
            *(finite_differences + i) += (n - k % 2 == 0 ? 1 : -1) * binom(n, k) * f(x1, d);
            *(x1 + i) -= k * h;

        }

    }

    free(x1);

    return finite_differences;

}

double *backward_first(RealFunction f, double *x, double h, unsigned int d) {

    double *x1 = copyx(x, d);
    double *finite_differences;
    if (!(finite_differences = (double *)calloc(d, sizeof(double)))) { return NULL; }

    for (int i = 0; i < d; ++i) {

        *(x1 + i) -= h;
        *(finite_differences + i) = f(x, d) - f(x1, d);
        *(x1 + i) += h;

    }

    free(x1);

    return finite_differences;

}

double *backward_second(RealFunction f, double *x, double h, unsigned int d) {

    double *x1 = copyx(x, d);
    double *x2 = copyx(x, d);
    double *finite_differences;
    if (!(finite_differences = (double *)calloc(d, sizeof(double)))) { return NULL; }

    for (int i = 0; i < d; ++i) {

        *(x1 + i) -= h, *(x2 + i) -= 2 * h;
        *(finite_differences + i) = f(x, d) - 2 * f(x1, d) + f(x2, d);
        *(x1 + i) += h, *(x2 + i) += 2 * h;

    }

    free(x1);
    free(x2);

    return finite_differences;

}

double *backward_nth(RealFunction f, double *x, double h, unsigned int n, unsigned int d) {

    double *x1 = copyx(x, d);
    double *finite_differences;
    if (!(finite_differences = (double *)calloc(d, sizeof(double)))) { return NULL; }

    for (int i = 0; i < d; ++i) {

        *(finite_differences + i) = 0.;

        for (int k = 0; k <= n; ++i) {

            *(x1 + i) -= k * h;
            *(finite_differences + i) += (i % 2 == 0 ? 1 : -1) * binom(n, k) * f(x1, d);
            *(x1 + i) += k * h;

        }

    }

    free(x1);

    return finite_differences;

}

double *central_first(RealFunction f, double *x, double h, unsigned int d) {

    double *x1 = copyx(x, d);
    double *x2 = copyx(x, d);
    double *finite_differences;
    if (!(finite_differences = (double *)calloc(d, sizeof(double)))) { return NULL; }

    for (int i = 0; i < d; ++i) {

        *(x1 + i) += h / 2, *(x2 + i) -= h / 2;
        *(finite_differences + i) = f(x1, d) - f(x2, d);
        *(x1 + i) -= h / 2, *(x2 + i) += h / 2;

    }

    free(x1);
    free(x2);

    return finite_differences;

}

double *central_second(RealFunction f, double *x, double h, unsigned int d) {

    double *x1 = copyx(x, d);
    double *x2 = copyx(x, d);
    double *finite_differences;
    if (!(finite_differences = (double *)calloc(d, sizeof(double)))) { return NULL; }

    for (int i = 0; i < d; ++i) {

        *(x1 + i) += h, *(x2 + i) -= h;
        *(finite_differences + i) = f(x1, d) - 2 * f(x, d) + f(x2, d);
        *(x1 + i) -= h, *(x2 + i) += h;

    }

    free(x1);
    free(x2);

    return finite_differences;

}

double *central_nth(RealFunction f, double *x, double h, unsigned int n, unsigned int d) {

    double *x1 = copyx(x, d);
    double *finite_differences;
    if (!(finite_differences = (double *)calloc(d, sizeof(double)))) { return NULL; }

    for (int i = 0; i < d; ++i) {

        *(finite_differences + i) = 0.;

        for (int k = 0; k <= n; ++i) {

            *(x1 + i) += (n / 2 - k) * h;
            *(finite_differences + i) += (i % 2 == 0 ? 1 : -1) * binom(n, k) * f(x1, d);
            *(x1 + i) -= (n / 2 - k) * h;

        }

    }

    free(x1);

    return finite_differences;

}


double *difference_quotient(
    RealFunction f, double *x, double h, unsigned int n, unsigned int d,
    struct FiniteDifference *findiff
) {

    double *finite_differences;
    switch(n) {
        case 1:
            finite_differences = findiff->first(f, x, h, d);
            break;
        case 2:
            finite_differences = findiff->second(f, x, h, d);
            break;
        default:
            finite_differences = findiff->nth(f, x, h, n, d);
            break;
    }

    double step = power(h, n);
    for (int i = 0; i < d; ++i) {
        *(finite_differences + i) /= step;
    }

    return finite_differences;

}
