/**
 * Source file for "../include/numbers.h"
 */

#include <stdlib.h>

#include "../include/numbers.h"


double power(double b, unsigned int p) {

    double res = 1.;
    for (int i = 0; i < p; ++i) {
        res *= b;
    }
    return res;

}


unsigned long factorial_(unsigned int n, unsigned long res) {

    if (n == 0 || n == 1) {
        return res;
    }

    return factorial_(n - 1, n * res);

}


unsigned long factorial(unsigned int n) { return factorial_(n, 1); }


long int binom(int alpha, unsigned int n) {

    if (n == 0) { return 1; }

    long int res = 1L;
    for (int k = 1; k <= n; ++k) {
        res *= (alpha - k + 1) / k;
    }

    return res;

}
