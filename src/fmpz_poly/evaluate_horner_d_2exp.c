/*
    Copyright (C) 2016 William Hart
    Copyright (C) 2021, 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include <double_extras.h>
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

FLINT_FORCE_INLINE dpe_t dpe_set_d_exp(double x, slong e)
{
    dpe_t res;
    res.m = x;
    res.e = e;
    DPE_ADJUST(res);
    return res;
}

#if 0
FLINT_FORCE_INLINE dpe_t dpe_set_fmpz(const fmpz_t x)
{
    dpe_t res;
    res.m = fmpz_get_d_2exp(&res.e, x);
    return res;
}
#endif

FLINT_FORCE_INLINE dpe_t
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

        res.m = x.m + d_mul_2exp_inrange(y.m, -d);
        res.e = x.e;
    }
    else
    {
        d = -d;

        if (d > 53 + ADJUSTMENT_DELAY)
            return y;

        res.m = y.m + d_mul_2exp_inrange(x.m, -d);
        res.e = y.e;
    }

    /* We delay adjustments */
    /* DPE_ADJUST(res); */
    return res;
}

FLINT_FORCE_INLINE dpe_t
dpe_mul(dpe_t x, dpe_t y)
{
    dpe_t res;
    res.m = x.m * y.m;
    res.e = x.e + y.e;

    /* We delay adjustments */
    /* DPE_ADJUST(res); */
    return res;
}

#if defined(__GNUC__)
# pragma GCC diagnostic push
# pragma GCC diagnostic ignored "-Wmissing-prototypes"
#endif
/* Used in CLD_bound.c */
double _fmpz_poly_evaluate_horner_d_2exp2_precomp(slong * exp, const double * poly, const slong * poly_exp, slong n, double d, slong dexp)
{
    dpe_t s, t, x;
    slong i;

    if (d == 0.0)
    {
        *exp = poly_exp[0];
        return poly[0];
    }

    x = dpe_set_d_exp(d, dexp);

    s.m = poly[n - 1];
    s.e = poly_exp[n - 1];

    for (i = n - 2; i >= 0; i--)
    {
        s = dpe_mul(s, x);

        if (poly[i] != 0.0)
        {
            t.m = poly[i];
            t.e = poly_exp[i];
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
#if defined(__GNUC__)
# pragma GCC diagnostic pop
#endif

double _fmpz_poly_evaluate_horner_d_2exp2(slong * exp, const fmpz * poly, slong n, double d, slong dexp)
{
    double * p;
    slong * p_exp;
    slong i;
    double v;
    TMP_INIT;

    if (d == 0.0)
        return fmpz_get_d_2exp(exp, poly + 0);

    TMP_START;

    p = TMP_ALLOC(n * sizeof(double));
    p_exp = TMP_ALLOC(n * sizeof(slong));

    for (i = 0; i < n; i++)
    {
        p[i] = fmpz_get_d_2exp(&p_exp[i], poly + i);
    }

    /* GCC really wants to complain about this one */
#if defined(__GNUC__) && !defined(__clang__)
# pragma GCC diagnostic push
# pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
    v = _fmpz_poly_evaluate_horner_d_2exp2_precomp(exp, p, p_exp, n, d, dexp);
#if defined(__GNUC__) && !defined(__clang__)
# pragma GCC diagnostic pop
#endif

    TMP_END;

    return v;
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
