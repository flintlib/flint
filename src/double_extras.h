/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef DOUBLE_EXTRAS_H
#define DOUBLE_EXTRAS_H

#ifdef DOUBLE_EXTRAS_INLINES_C
#define DOUBLE_EXTRAS_INLINE
#else
#define DOUBLE_EXTRAS_INLINE static inline
#endif

#include <math.h>
#include "flint.h"

#ifdef __cplusplus
 extern "C" {
#endif

#define D_BITS 53
#define D_EPS 2.2204460492503130808e-16
#define D_INF HUGE_VAL
#define D_NAN (HUGE_VAL - HUGE_VAL)

double d_randtest(flint_rand_t state);

double d_randtest_signed(flint_rand_t state, slong minexp, slong maxexp);

double d_randtest_special(flint_rand_t state, slong minexp, slong maxexp);

DOUBLE_EXTRAS_INLINE
double d_polyval(const double * poly, int len, double x)
{
    double t;
    int i;

    for (t = poly[len-1], i = len-2; i >= 0; i--)
        t = poly[i] + x * t;

    return t;
}

double d_lambertw(double x);

int d_is_nan(double x);

double d_log2(double x);

#ifdef __cplusplus
}
#endif

#endif
