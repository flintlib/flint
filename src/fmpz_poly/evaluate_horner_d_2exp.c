/*
    Copyright (C) 2016 William Hart
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include <math.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"

/* Naive double+exponent arithmetic; not designed to deal with
   underflow/overflow. */
typedef struct
{
    double m;
    slong e;
}
dpe_t;

/* Do up to ADJUSTMENT_DELAY steps without adjusting mantissa.
   This may cause the mantissa to drift from the normalized range [0.5, 1)
   by a factor about 2^ADJUSTMENT_DELAY. */
#define ADJUSTMENT_DELAY 16

/* Standardize nonzero mantissa to +/- [0.5, 1). */
#define DPE_ADJUST(x) \
    do { \
        int _e; \
        (x).m = frexp((x).m, &_e); \
        (x).e += _e; \
    } while (0)

static dpe_t dpe_set_d_exp(double x, slong e)
{
    dpe_t res;
    res.m = x;
    res.e = e;
    DPE_ADJUST(res);
    return res;
}

static dpe_t dpe_set_fmpz(const fmpz_t x)
{
    dpe_t res;
    res.m = fmpz_get_d_2exp(&res.e, x);
    return res;
}

static dpe_t
dpe_add(dpe_t x, dpe_t y)
{
    dpe_t res;
    slong d;

    d = x.e - y.e;

    if (x.m == 0.0)
        return y;

    if (y.m == 0.0)
        return x;

    if (d >= 0)
    {
        if (d > 53 + ADJUSTMENT_DELAY)
            return x;

        res.m = x.m + ldexp(y.m, -d);
        res.e = x.e;
    }
    else
    {
        d = -d;

        if (d > 53 + ADJUSTMENT_DELAY)
            return y;

        res.m = y.m + ldexp(x.m, -d);
        res.e = y.e;
    }

    /* We delay adjustments */
    /* DPE_ADJUST(res); */
    return res;
}

static dpe_t
dpe_mul(dpe_t x, dpe_t y)
{
    dpe_t res;
    res.m = x.m * y.m;
    res.e = x.e + y.e;

    /* We delay adjustments */
    /* DPE_ADJUST(res); */
    return res;
}

double _fmpz_poly_evaluate_horner_d_2exp2(slong * exp, const fmpz * poly, slong n, double d, slong dexp)
{
    dpe_t s, t, x;
    slong i;

    if (d == 0.0)
        return fmpz_get_d_2exp(exp, poly + 0);

    x = dpe_set_d_exp(d, dexp);
    s = dpe_set_fmpz(poly + n - 1);

    for (i = n - 2; i >= 0; i--)
    {
        s = dpe_mul(s, x);

        if (!fmpz_is_zero(poly + i))
        {
            t = dpe_set_fmpz(poly + i);
            s = dpe_add(s, t);
        }

        /* Delayed adjustments. */
        if (i % ADJUSTMENT_DELAY == 0)
            DPE_ADJUST(s);
    }

    DPE_ADJUST(s);

    *exp = s.e;
    return s.m;
}

double _fmpz_poly_evaluate_horner_d_2exp(slong * exp, const fmpz * poly, slong n, double d)
{
    return _fmpz_poly_evaluate_horner_d_2exp2(exp, poly, n, d, 0);
}

double fmpz_poly_evaluate_horner_d_2exp2(slong * exp, const fmpz_poly_t poly, double d, slong dexp)
{
   if (poly->length == 0)
   {
      *exp = 0;
      return 0.0;
   }

   return _fmpz_poly_evaluate_horner_d_2exp2(exp, poly->coeffs, poly->length, d, dexp);
}


double fmpz_poly_evaluate_horner_d_2exp(slong * exp, const fmpz_poly_t poly, double d)
{
   if (poly->length == 0)
   {
      *exp = 0;
      return 0.0;
   }

   return _fmpz_poly_evaluate_horner_d_2exp(exp, poly->coeffs, poly->length, d);
}
