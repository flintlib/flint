/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Conversions to and from arf/arb, strings, and random generation,
   for any n. */

#include "fp_contract.h"
#include <float.h>
#include <stdio.h>
#include <string.h>
#include "dfloat.h"
#include "double_extras.h"
#include "arf.h"
#include "arb.h"
#include "mag.h"
#include "gr.h"

/* exact: the components are doubles */
void
_dfloat_get_arf(arf_t res, const double * x, int n)
{
    arf_t t;
    int i;
    arf_set_d(res, x[0]);
    if (n > 1)
    {
        arf_init(t);
        for (i = 1; i < n; i++)
        {
            arf_set_d(t, x[i]);
            arf_add(res, res, t, ARF_PREC_EXACT, ARF_RND_DOWN);
        }
        arf_clear(t);
    }
}

/* greedy rounding: res[i] = round(x - res[0] - ... - res[i-1]); the
   remainder after n components is added to *err (an upper bound) when
   err is not NULL.  Nonfinite x gives a nonfinite head (the ball
   versions check for that before calling). */
void
_dfloat_set_arf(double * res, int n, double * err, const arf_t x)
{
    arf_t t, u;
    int i;

    if (!arf_is_finite(x))
    {
        res[0] = arf_get_d(x, ARF_RND_NEAR);
        for (i = 1; i < n; i++)
            res[i] = 0.0;
        if (err != NULL)
            *err = D_INF;
        return;
    }

    arf_init(t);
    arf_init(u);
    arf_set(t, x);
    for (i = 0; i < n; i++)
    {
        res[i] = arf_get_d(t, ARF_RND_NEAR);
        /* an infinite head means the value is out of range; the
           callers treat that as nonfinite */
        if (!(fabs(res[i]) <= DBL_MAX))
        {
            for (i = i + 1; i < n; i++)
                res[i] = 0.0;
            if (err != NULL)
                *err = D_INF;
            arf_clear(t);
            arf_clear(u);
            return;
        }
        arf_set_d(u, res[i]);
        arf_sub(t, t, u, ARF_PREC_EXACT, ARF_RND_DOWN);
    }

    if (err != NULL && !arf_is_zero(t))
    {
        mag_t m;
        mag_init(m);
        arf_get_mag(m, t);
        *err += mag_get_d(m);
        mag_clear(m);
    }

    arf_clear(t);
    arf_clear(u);
}

void
_dfloat_get_arb(arb_t res, const double * x, int n, double rad)
{
    _dfloat_get_arf(arb_midref(res), x, n);
    if (rad == 0.0)
        mag_zero(arb_radref(res));
    else if (rad == D_INF)
        mag_inf(arb_radref(res));
    else
        mag_set_d(arb_radref(res), rad);   /* rounds up */
}

void
_dfloat_set_arb(double * res, double * rad, int n, const arb_t x)
{
    double err = 0.0;
    int i;

    if (!arf_is_finite(arb_midref(x)))
    {
        for (i = 0; i < n; i++)
            res[i] = 0.0;
        *rad = D_INF;
        return;
    }

    _dfloat_set_arf(res, n, &err, arb_midref(x));
    err += mag_get_d(arb_radref(x));   /* an upper bound */
    *rad = _dfloat_rad_finish(err);

    if (!(_dfloat_abs_sum(res, n) + *rad <= DBL_MAX))
    {
        for (i = 0; i < n; i++)
            res[i] = 0.0;
        *rad = D_INF;
    }
}

/* strings go through arb: decimal input is not exact in general, so
   the plain types round and the balls carry the conversion error */
int
_dfloat_set_str(double * res, double * rad, int n, const char * s)
{
    arb_t t;
    int status;
    double r;

    arb_init(t);
    status = arb_set_str(t, s, 53 * n + 64);
    if (status == 0)
    {
        _dfloat_set_arb(res, &r, n, t);
        if (rad != NULL)
            *rad = r;
    }
    arb_clear(t);
    return status ? GR_UNABLE : GR_SUCCESS;
}

/* rad < 0 means a plain value (printed without radius) */
char *
_dfloat_get_str(const double * x, int n, double rad, slong digits)
{
    arb_t t;
    char * s;

    if (digits <= 0)
        digits = (slong) (53 * n * 0.30102999566398) + 1;

    arb_init(t);
    if (!(rad < 0.0) && _dfloat_ball_is_whole(x, n, rad))
        arb_zero_pm_inf(t);          /* every form of the whole line */
    else
        _dfloat_get_arb(t, x, n, (rad < 0.0) ? 0.0 : rad);
    if (rad < 0.0 && !arf_is_finite(arb_midref(t)))
    {
        s = flint_malloc(8);
        strcpy(s, arf_is_nan(arb_midref(t)) ? "nan" : (arf_sgn(arb_midref(t)) > 0 ? "inf" : "-inf"));
    }
    else
        s = arb_get_str(t, digits, (rad < 0.0) ? ARB_STR_NO_RADIUS : 0);
    arb_clear(t);
    return s;
}

/* random expansion: a random head with a moderate exponent, and
   tails that fill the precision, or a short (exactly representable)
   value; the results are renormalized */
void
_dfloat_randtest(double * res, int n, flint_rand_t state)
{
    double t[DFLOAT_MAX_N];
    int i, k;

    switch (n_randint(state, 8))
    {
        case 0:
            /* small integer */
            t[0] = (double) (slong) n_randint(state, 20) - 10;
            for (i = 1; i < n; i++)
                t[i] = 0.0;
            break;
        case 1:
            /* zero */
            for (i = 0; i < n; i++)
                t[i] = 0.0;
            break;
        case 2:
            /* short dyadic */
            t[0] = d_randtest_signed(state, -10, 10);
            for (i = 1; i < n; i++)
                t[i] = 0.0;
            break;
        default:
            k = 1 + n_randint(state, n);
            t[0] = d_randtest_signed(state, -30, 30);
            for (i = 1; i < k; i++)
                t[i] = t[i - 1] * d_randtest_signed(state, -56, -52);
            for (; i < n; i++)
                t[i] = 0.0;
            break;
    }

    _dfloat_renorm(res, n, NULL, t, n);
}

/* extreme exponents, subnormals, the ends of the double range */
void
_dfloat_randtest_special(double * res, int n, flint_rand_t state)
{
    double t[DFLOAT_MAX_N];
    int i, k;

    switch (n_randint(state, 6))
    {
        case 0:
            _dfloat_randtest(res, n, state);
            return;
        case 1:
            t[0] = d_randtest_signed(state, -1074, 1023);
            for (i = 1; i < n; i++)
                t[i] = 0.0;
            break;
        case 2:
            k = 1 + n_randint(state, n);
            t[0] = d_randtest_signed(state, -1074, 1023);
            for (i = 1; i < k; i++)
                t[i] = t[i - 1] * d_randtest_signed(state, -60, -52);
            for (; i < n; i++)
                t[i] = 0.0;
            break;
        case 3:
            /* huge head, tiny tails */
            k = 1 + n_randint(state, n);
            t[0] = d_randtest_signed(state, -1074, 1023);
            for (i = 1; i < k; i++)
                t[i] = d_randtest_signed(state, -1074, -1000);
            for (; i < n; i++)
                t[i] = 0.0;
            break;
        case 4:
            t[0] = (n_randint(state, 2) ? DBL_MAX : -DBL_MAX) * (n_randint(state, 2) ? 1.0 : 0.5);
            for (i = 1; i < n; i++)
                t[i] = n_randint(state, 2) ? 0.0 : d_randtest_signed(state, 900, 960);
            break;
        default:
            for (i = 0; i < n; i++)
                t[i] = n_randint(state, 2) ? 0.0 : d_randtest_signed(state, -1074, -1020);
            break;
    }

    _dfloat_renorm(res, n, NULL, t, n);
}

/* a random radius: 0, tiny, moderate, huge or infinite */
double
_dfloat_randtest_rad(flint_rand_t state)
{
    switch (n_randint(state, 8))
    {
        case 0:
        case 1:
        case 2:
            return 0.0;
        case 3:
            return fabs(d_randtest_signed(state, -220, -50));
        case 4:
            return fabs(d_randtest_signed(state, -60, 0));
        case 5:
            return fabs(d_randtest_signed(state, -1074, 1023));
        case 6:
            return D_INF;
        default:
            return fabs(d_randtest_signed(state, 0, 30));
    }
}

/* conversions between formats: rounding to smaller n, exact to larger */
void
dfloat_set_dfloat(double * res, int nres, const double * x, int nx)
{
    double t[DFLOAT_MAX_N];
    int i;
    for (i = 0; i < nx; i++)
        t[i] = x[i];
    for (; i < nres; i++)
        t[i] = 0.0;
    _dfloat_renorm(res, nres, NULL, t, FLINT_MAX(nx, nres));
}

void
dfloat_ball_set_ball(double * res, int nres, const double * x, int nx)
{
    double t[DFLOAT_MAX_N], err = 0.0;
    int i;
    if (_dfloat_ball_is_whole(x, nx, x[nx]))
    {
        for (i = 0; i < nres; i++)
            res[i] = 0.0;
        res[nres] = D_INF;
        return;
    }
    for (i = 0; i < nx; i++)
        t[i] = x[i];
    for (; i < nres; i++)
        t[i] = 0.0;
    _dfloat_renorm(res, nres, &err, t, FLINT_MAX(nx, nres));
    res[nres] = _dfloat_rad_finish(x[nx] + err);
    if (!(_dfloat_abs_sum(res, nres) + res[nres] <= DBL_MAX))
    {
        for (i = 0; i < nres; i++)
            res[i] = 0.0;
        res[nres] = D_INF;
    }
}

/* a complex number from its parts (a radius < 0: plain), in the format
   of the acb ring: "re", "im*I" or "(re + im*I)", with " - " for a
   negative exact imaginary part */
char *
_dfloat_complex_get_str(const double * re, const double * im, int n, double rrad, double irad, slong digits)
{
    char * a, * b, * s;
    int i, re_zero = (rrad <= 0.0), im_zero = (irad <= 0.0), neg;
    double t[DFLOAT_MAX_N] = {0.0, 0.0, 0.0, 0.0};

    for (i = 0; i < n; i++)
    {
        if (re[i] != 0.0) re_zero = 0;
        if (im[i] != 0.0) im_zero = 0;
    }

    if (im_zero)
        return _dfloat_get_str(re, n, rrad, digits);

    neg = (irad <= 0.0) && im[0] < 0.0;
    for (i = 0; i < n; i++)
        t[i] = neg ? -im[i] : im[i];
    b = _dfloat_get_str(t, n, irad, digits);

    if (re_zero)
    {
        s = flint_malloc(strlen(b) + 4);
        sprintf(s, "%s%s*I", neg ? "-" : "", b);
    }
    else
    {
        a = _dfloat_get_str(re, n, rrad, digits);
        s = flint_malloc(strlen(a) + strlen(b) + 8);
        sprintf(s, "(%s %s %s*I)", a, neg ? "-" : "+", b);
        flint_free(a);
    }
    flint_free(b);
    return s;
}

/* a part (n components, and a radius for a ball) of an element of a
   dfloat ring as an arb: 0 for a nonfinite plain part (not a number),
   the whole line for any nonfinite ball */
int
_dfloat_part_get_arb(arb_t res, const double * x, int n, int ball)
{
    if (ball)
    {
        if (_dfloat_ball_is_whole(x, n, x[n]))
            arb_zero_pm_inf(res);
        else
            _dfloat_get_arb(res, x, n, x[n]);
        return 1;
    }
    if (!(_dfloat_abs_sum(x, n) <= DBL_MAX))
        return 0;
    _dfloat_get_arb(res, x, n, 0.0);
    return 1;
}
