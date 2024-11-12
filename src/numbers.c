/**
 * Source file for "../include/numbers.h"
 */

#include <stdlib.h>

#include "../include/numbers.h"


double power(double b, unsigned int p) {

    if (p == 0) { return 1.; }

    if (p % 2 == 0) {
        double temp = power(b, p / 2);
        return temp * temp;
    } else {
        return b * power(b, p - 1);
    }

}

unsigned long int intpower(unsigned int b, unsigned int p) {

    unsigned long int res = 1L;

    while (p > 0) {

        if (p & 1) { res *= b; }
        res *= res;
        p >>= 1;

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
