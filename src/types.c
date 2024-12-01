/**
 * Source file for "../include/types.h"
 */

#include "../include/types.h"

struct RealFunction *parse_function(PyObject *ob_f) {

    if (!PYCallable_Check(ob_f)) {
        PyErr_SetString(PyExc_TypeError, "Expected a callable object");
        return NULL;
    }

    struct RealFunction *f = (struct RealFunction *)malloc(sizeof(struct RealFunction));
    if (!f) {
        PyErr_MemoryError("Failed to allocated memory");
        return NULL;
    }

    f->c = false;
    f->f.py = ob_f;

    return f;

}
