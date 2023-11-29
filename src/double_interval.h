/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef DOUBLE_INTERVAL_H
#define DOUBLE_INTERVAL_H

#ifdef DOUBLE_INTERVAL_INLINES_C
#define DOUBLE_INTERVAL_INLINE
#else
#define DOUBLE_INTERVAL_INLINE static inline
#endif

#include <math.h>
#include "double_extras.h"
#include "arb_types.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
    double a;
    double b;
}
di_t;

#define DI_CHECK(__x) \
    if (!(__x.a <= __x.b)) \
    { \
        flint_throw(FLINT_ERROR, "di_t endpoints %g, %g not ordered\n", __x.a, __x.b); \
    } \

DOUBLE_INTERVAL_INLINE
di_t di_interval(double a, double b)
{
    di_t res;

    if (!(a <= b))
    {
        flint_throw(FLINT_ERROR, "di_interval endpoints %g, %g not ordered\n", a, b);
    }

    res.a = a;
    res.b = b;
    return res;
}

DOUBLE_INTERVAL_INLINE
double _di_below(double x)
{
    double t;

    if (x <= 1e300)
    {
        t = x;
        if (t < 0.0)
            t = -t;

        t += 1e-300;
        return x - t * 4.440892098500626e-16;
    }
    else
    {
        if (x != x)
            return -D_INF;

        return 1e300;
    }
}

DOUBLE_INTERVAL_INLINE
double _di_above(double x)
{
    double t;

    if (x >= -1e300)
    {
        t = x;
        if (t < 0.0)
            t = -t;

        t += 1e-300;
        return x + t * 4.440892098500626e-16;
    }
    else
    {
        if (x != x)
            return D_INF;

        return -1e300;
    }
}

DOUBLE_INTERVAL_INLINE
di_t di_neg(di_t x)
{
    di_t res;
    res.a = -x.b;
    res.b = -x.a;
    return res;
}

DOUBLE_INTERVAL_INLINE
di_t di_fast_add(di_t x, di_t y)
{
    di_t res;
    res.a = _di_below(x.a + y.a);
    res.b = _di_above(x.b + y.b);
    return res;
}

DOUBLE_INTERVAL_INLINE
di_t di_fast_sub(di_t x, di_t y)
{
    di_t res;
    res.a = _di_below(x.a - y.b);
    res.b = _di_above(x.b - y.a);
    return res;
}

di_t di_fast_mul(di_t x, di_t y);
di_t di_fast_sqr(di_t x);
di_t di_fast_div(di_t x, di_t y);

DOUBLE_INTERVAL_INLINE
di_t di_fast_add_d(di_t x, double y)
{
    return di_fast_add(x, di_interval(y, y));
}

DOUBLE_INTERVAL_INLINE
di_t di_fast_sub_d(di_t x, double y)
{
    return di_fast_sub(x, di_interval(y, y));
}

DOUBLE_INTERVAL_INLINE
di_t di_fast_mul_d(di_t x, double y)
{
    return di_fast_mul(x, di_interval(y, y));
}

DOUBLE_INTERVAL_INLINE
di_t di_fast_div_d(di_t x, double y)
{
    return di_fast_div(x, di_interval(y, y));
}

di_t di_fast_log_nonnegative(di_t x);

DOUBLE_INTERVAL_INLINE
di_t di_fast_mid(di_t x)
{
    di_t a, b;
    if (x.a == -D_INF || x.b == D_INF)
        return di_interval(-D_INF, D_INF);
    a = di_interval(x.a, x.a);
    b = di_interval(x.b, x.b);
    return di_fast_mul_d(di_fast_add(a, b), 0.5);
}

DOUBLE_INTERVAL_INLINE
double di_fast_ubound_radius(di_t x)
{
    return _di_above((x.b - x.a) * 0.5);
}

DOUBLE_INTERVAL_INLINE
void di_print(di_t x)
{
    flint_printf("[%.17g, %.17g]", x.a, x.b);
}

di_t arb_get_di(const arb_t x);
void arb_set_di(arb_t res, di_t x, slong prec);

DOUBLE_INTERVAL_INLINE
double d_randtest2(flint_rand_t state)
{
    double x;

    x = d_randtest(state);
    if (n_randint(state, 2))
        x = -x;

    return ldexp(x, n_randint(state, 2400) - 1200);
}

DOUBLE_INTERVAL_INLINE
di_t di_randtest(flint_rand_t state)
{
    di_t res;

    res.a = d_randtest2(state);
    res.b = d_randtest2(state);

    if (res.a > res.b)
    {
        double t = res.a;
        res.a = res.b;
        res.b = t;
    }

    return res;
}

#ifdef __cplusplus
}
#endif

#endif

