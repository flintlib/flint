/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef DOUBLE_EXTRAS_H
#define DOUBLE_EXTRAS_H

#ifdef DOUBLE_EXTRAS_INLINES_C
#define DOUBLE_EXTRAS_INLINE
#else
#define DOUBLE_EXTRAS_INLINE static inline
#endif

#include <stdint.h>
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

DOUBLE_EXTRAS_INLINE
int d_is_nan(double x)
{
    return x != x;
}

double d_log2(double x);

typedef union
{
    double f;
    uint64_t i;
}
double_uint64_u;

#define D_MIN_NORMAL_EXPONENT -1022
#define D_MAX_NORMAL_EXPONENT 1023
#define D_EXPONENT_BIAS 1023
#define D_EXPONENT_SHIFT 52

/* Assumes that 2^i is in the normal exponent range. */
FLINT_FORCE_INLINE double d_mul_2exp_inrange(double x, int i)
{
    FLINT_ASSERT(i >= D_MIN_NORMAL_EXPONENT && i <= D_MAX_NORMAL_EXPONENT);
    double_uint64_u u;
    u.i = ((int64_t) (i + D_EXPONENT_BIAS)) << D_EXPONENT_SHIFT;
    return x * u.f;
}

/* Assumes that 2^i, x and x*2^i are all in the normal exponent range. */
/* In particular, also assumes x != 0. */
FLINT_FORCE_INLINE double d_mul_2exp_inrange2(double x, int i)
{
    FLINT_ASSERT(i >= D_MIN_NORMAL_EXPONENT && i <= D_MAX_NORMAL_EXPONENT);
    FLINT_ASSERT(x != 0);

    double_uint64_u u;
    u.f = x;
    u.i += ((int64_t) i) << D_EXPONENT_SHIFT;
    return u.f;
}

FLINT_FORCE_INLINE double d_mul_2exp(double x, int i)
{
    if (i >= D_MIN_NORMAL_EXPONENT && i <= D_MAX_NORMAL_EXPONENT)
        return d_mul_2exp_inrange(x, i);
    else
        return ldexp(x, i);
}

#ifdef __cplusplus
}
#endif

#endif
