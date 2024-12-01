/**
 * Type definitions for mathematical entities
 */

#include <stdbool.h>

struct RealFunction {
    bool c;
    union {
        double (*c)(double *x, unsigned int d);
        PyObject *py;
    } f;
};
