/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* The exponential function.  The per-N kernel dN__exp_core lives in
   template.inc; it is a pure floating-point computation whose
   worst-case relative error DFLOAT_EXP_EPS_N is derived once and for
   all in dev/dfloat_exp_bound.py, so that the ball version needs no error
   tracking: it evaluates the kernel on the midpoint and adds the
   static bound and the propagated input radius.

   Kernel, for x given as a canonical N-term expansion (the input is
   canonicalised first):

     k = round(x0 / log 2), h = x0 - k L0 (exact: L0 is log 2 to 41
     bits, so k L0 is exact), i = round(256 h), j = round(65536 (h -
     i/256)), s0 = h - i/256 - j/65536 (exact, |s0| <= 2^-17);
     s = s0 + (x1 + ... + x_{N-1}) - k (L1 + ... + L_N), summed with
     the plain N-term kernels in an order that keeps everything
     significant, |s| <= 2^-17 + 2^-31;
     exp(x) = 2^k exp(i/256) exp(j/65536) (1 + p(s)), p(s) the Taylor
     polynomial of expm1 of degree 3, 6, 8, 11, by Horner's rule with
     tapered precision (the inner levels use fewer components), the
     tables and coefficients being nearest-rounded expansions.

   All exactness claims that the analysis relies on (the three
   two-sums in the reduction and the canonical form of the input) are
   verified at run time, and the kernel defers to arb in the (never
   observed) cases where they fail.

   The input radius r enters at the end: for |x' - m| <= r,
   |exp(x') - exp(m)| <= exp(m) (exp(r) - 1), with exp(m) bounded by
   the computed value and exp(r) - 1 by r + r^2 for r <= 1 and by a
   power of two otherwise. */

#include "fp_contract.h"
#include <float.h>
#include "dfloat.h"
#include "double_extras.h"
#include "arb.h"
#include "gr.h"

void
_dfloat_exp_slow(double * res, int n, const double * x)
{
    arb_t t;
    slong prec;
    double err = 0.0;

    arb_init(t);
    _dfloat_get_arb(t, x, n, 0.0);
    if (!arf_is_finite(arb_midref(t)))
    {
        slong i;
        for (i = 0; i < n; i++)
            res[i] = (i == 0) ? D_NAN : 0.0;
        arb_clear(t);
        return;
    }
    for (prec = 53 * n + 64; ; prec *= 2)
    {
        arb_exp(t, t, prec);
        if (arb_rel_accuracy_bits(t) >= 53 * n + 40)
            break;
        _dfloat_get_arb(t, x, n, 0.0);
    }
    _dfloat_set_arf(res, n, &err, arb_midref(t));
    arb_clear(t);
}

void
_dfloat_ball_exp(double * res, double * rrad, int n, const double * x, double rx)
{
    double b[DFLOAT_MAX_N + 1], r[DFLOAT_MAX_N + 1];
    int i;
    for (i = 0; i < n; i++)
        b[i] = x[i];
    b[n] = rx;
    switch (n)
    {
        case 1: d1b_exp((d1b_ptr) r, (d1b_srcptr) b); break;
        case 2: d2b_exp((d2b_ptr) r, (d2b_srcptr) b); break;
        case 3: d3b_exp((d3b_ptr) r, (d3b_srcptr) b); break;
        default: d4b_exp((d4b_ptr) r, (d4b_srcptr) b); break;
    }
    for (i = 0; i < n; i++)
        res[i] = r[i];
    *rrad = r[n];
}

void
_dfloat_exp(double * res, int n, const double * x)
{
    switch (n)
    {
        case 1: d1_exp((d1_ptr) res, (d1_srcptr) x); break;
        case 2: d2_exp((d2_ptr) res, (d2_srcptr) x); break;
        case 3: d3_exp((d3_ptr) res, (d3_srcptr) x); break;
        default: d4_exp((d4_ptr) res, (d4_srcptr) x); break;
    }
}

void
_dfloat_sin_cos(double * sn, double * cs, int n, const double * x, int want)
{
    double s[DFLOAT_MAX_N], c[DFLOAT_MAX_N];
    int i;
    if (want == 1)
    {
        switch (n)
        {
            case 1: d1_sin((d1_ptr) sn, (d1_srcptr) x); break;
            case 2: d2_sin((d2_ptr) sn, (d2_srcptr) x); break;
            case 3: d3_sin((d3_ptr) sn, (d3_srcptr) x); break;
            default: d4_sin((d4_ptr) sn, (d4_srcptr) x); break;
        }
    }
    else if (want == 2)
    {
        switch (n)
        {
            case 1: d1_cos((d1_ptr) cs, (d1_srcptr) x); break;
            case 2: d2_cos((d2_ptr) cs, (d2_srcptr) x); break;
            case 3: d3_cos((d3_ptr) cs, (d3_srcptr) x); break;
            default: d4_cos((d4_ptr) cs, (d4_srcptr) x); break;
        }
    }
    else
    {
        switch (n)
        {
            case 1: d1_sin_cos((d1_ptr) s, (d1_ptr) c, (d1_srcptr) x); break;
            case 2: d2_sin_cos((d2_ptr) s, (d2_ptr) c, (d2_srcptr) x); break;
            case 3: d3_sin_cos((d3_ptr) s, (d3_ptr) c, (d3_srcptr) x); break;
            default: d4_sin_cos((d4_ptr) s, (d4_ptr) c, (d4_srcptr) x); break;
        }
        for (i = 0; i < n; i++)
            sn[i] = s[i], cs[i] = c[i];
    }
}

void
_dfloat_ball_sin_cos(double * sn, double * srad, double * cs, double * crad,
    int n, const double * x, double rx, int want)
{
    double b[DFLOAT_MAX_N + 1], s[DFLOAT_MAX_N + 1], c[DFLOAT_MAX_N + 1];
    int i;
    for (i = 0; i < n; i++)
        b[i] = x[i];
    b[n] = rx;
    if (want == 1)
    {
        switch (n)
        {
            case 1: d1b_sin((d1b_ptr) s, (d1b_srcptr) b); break;
            case 2: d2b_sin((d2b_ptr) s, (d2b_srcptr) b); break;
            case 3: d3b_sin((d3b_ptr) s, (d3b_srcptr) b); break;
            default: d4b_sin((d4b_ptr) s, (d4b_srcptr) b); break;
        }
    }
    else if (want == 2)
    {
        switch (n)
        {
            case 1: d1b_cos((d1b_ptr) c, (d1b_srcptr) b); break;
            case 2: d2b_cos((d2b_ptr) c, (d2b_srcptr) b); break;
            case 3: d3b_cos((d3b_ptr) c, (d3b_srcptr) b); break;
            default: d4b_cos((d4b_ptr) c, (d4b_srcptr) b); break;
        }
    }
    else
    {
        switch (n)
        {
            case 1: d1b_sin_cos((d1b_ptr) s, (d1b_ptr) c, (d1b_srcptr) b); break;
            case 2: d2b_sin_cos((d2b_ptr) s, (d2b_ptr) c, (d2b_srcptr) b); break;
            case 3: d3b_sin_cos((d3b_ptr) s, (d3b_ptr) c, (d3b_srcptr) b); break;
            default: d4b_sin_cos((d4b_ptr) s, (d4b_ptr) c, (d4b_srcptr) b); break;
        }
    }
    if (want & 1)
    {
        for (i = 0; i < n; i++)
            sn[i] = s[i];
        *srad = s[n];
    }
    if (want & 2)
    {
        for (i = 0; i < n; i++)
            cs[i] = c[i];
        *crad = c[n];
    }
}

/* sine and cosine via arb, to the accuracy of the static bounds (which
   are absolute: a relative accuracy of 53 n + 40 bits on the arguments
   of size at most 2^21, or an absolute one where the result is small,
   suffices) */
void
_dfloat_trig_slow(double * sn, double * cs, int n, const double * x, int want)
{
    arb_t t, s, c;
    slong prec;
    double err;

    arb_init(t);
    arb_init(s);
    arb_init(c);
    _dfloat_get_arb(t, x, n, 0.0);
    if (!arb_is_finite(t))
    {
        if (want & 1)
            _dfloat_set_arf(sn, n, &err, arb_midref(t)), sn[0] = D_NAN;
        if (want & 2)
            _dfloat_set_arf(cs, n, &err, arb_midref(t)), cs[0] = D_NAN;
    }
    else
    {
        for (prec = 53 * n + 64; ; prec *= 2)
        {
            arb_sin_cos(s, c, t, prec);
            if (mag_cmp_2exp_si(arb_radref(s), -53 * n - 40) < 0 && mag_cmp_2exp_si(arb_radref(c), -53 * n - 40) < 0)
                break;
        }
        if (want & 1)
            _dfloat_set_arf(sn, n, &err, arb_midref(s));
        if (want & 2)
            _dfloat_set_arf(cs, n, &err, arb_midref(c));
    }
    arb_clear(t);
    arb_clear(s);
    arb_clear(c);
}

void
_dfloat_ball_trig_slow(double * sn, double * srad, double * cs, double * crad,
    int n, const double * x, double rx, int want)
{
    arb_t t, s, c;
    slong prec = 53 * n + 40;

    if (_dfloat_ball_is_whole(x, n, rx))
    {
        /* sin and cos of the whole line: [0 +/- 1] */
        int i;
        for (i = 0; i < n; i++)
        {
            if (want & 1) sn[i] = 0.0;
            if (want & 2) cs[i] = 0.0;
        }
        if (want & 1) *srad = 1.0;
        if (want & 2) *crad = 1.0;
        return;
    }

    arb_init(t);
    arb_init(s);
    arb_init(c);
    _dfloat_get_arb(t, x, n, rx);
    arb_sin_cos(s, c, t, prec);
    if (want & 1)
        _dfloat_set_arb(sn, srad, n, s);
    if (want & 2)
        _dfloat_set_arb(cs, crad, n, c);
    arb_clear(t);
    arb_clear(s);
    arb_clear(c);
}

/* ---- logarithm and arctangent: dispatchers and arb fallbacks ---- */

void
_dfloat_log(double * res, int n, const double * x)
{
    switch (n)
    {
        case 1: d1_log((d1_ptr) res, (d1_srcptr) x); break;
        case 2: d2_log((d2_ptr) res, (d2_srcptr) x); break;
        case 3: d3_log((d3_ptr) res, (d3_srcptr) x); break;
        default: d4_log((d4_ptr) res, (d4_srcptr) x); break;
    }
}

void
_dfloat_atan(double * res, int n, const double * x)
{
    switch (n)
    {
        case 1: d1_atan((d1_ptr) res, (d1_srcptr) x); break;
        case 2: d2_atan((d2_ptr) res, (d2_srcptr) x); break;
        case 3: d3_atan((d3_ptr) res, (d3_srcptr) x); break;
        default: d4_atan((d4_ptr) res, (d4_srcptr) x); break;
    }
}

int
_dfloat_ball_log(double * res, double * rrad, int n, const double * x, double rx)
{
    double b[DFLOAT_MAX_N + 1], r[DFLOAT_MAX_N + 1];
    int i, status;
    for (i = 0; i < n; i++)
        b[i] = x[i];
    b[n] = rx;
    switch (n)
    {
        case 1: status = d1b_log((d1b_ptr) r, (d1b_srcptr) b); break;
        case 2: status = d2b_log((d2b_ptr) r, (d2b_srcptr) b); break;
        case 3: status = d3b_log((d3b_ptr) r, (d3b_srcptr) b); break;
        default: status = d4b_log((d4b_ptr) r, (d4b_srcptr) b); break;
    }
    if (status == GR_SUCCESS)
    {
        for (i = 0; i < n; i++)
            res[i] = r[i];
        *rrad = r[n];
    }
    return status;
}

void
_dfloat_ball_atan(double * res, double * rrad, int n, const double * x, double rx)
{
    double b[DFLOAT_MAX_N + 1], r[DFLOAT_MAX_N + 1];
    int i;
    for (i = 0; i < n; i++)
        b[i] = x[i];
    b[n] = rx;
    switch (n)
    {
        case 1: d1b_atan((d1b_ptr) r, (d1b_srcptr) b); break;
        case 2: d2b_atan((d2b_ptr) r, (d2b_srcptr) b); break;
        case 3: d3b_atan((d3b_ptr) r, (d3b_srcptr) b); break;
        default: d4b_atan((d4b_ptr) r, (d4b_srcptr) b); break;
    }
    for (i = 0; i < n; i++)
        res[i] = r[i];
    *rrad = r[n];
}

/* the plain functions via arb to a relative accuracy of 53 n + 40
   bits (the static bounds are relative); the special values follow
   the double functions (log: nan for x < 0 or nan, -inf for 0, inf
   for inf; atan: nan for nan, +-pi/2 for +-inf) */
static void
_dfloat_unary_slow(double * res, int n, const double * x, void (*f)(arb_t, const arb_t, slong), int is_log)
{
    arb_t t;
    slong prec;
    double err = 0.0;

    arb_init(t);
    _dfloat_get_arb(t, x, n, 0.0);
    if (!arb_is_finite(t) || (is_log && !arb_is_positive(t)))
    {
        double v = x[0];
        if (is_log)
            v = (v == 0.0 && _dfloat_abs_sum(x, n) == 0.0) ? -D_INF : ((v > 0.0 && v == D_INF) ? D_INF : D_NAN);
        else
            v = (v == D_INF) ? 1.5707963267948966 : (v == -D_INF) ? -1.5707963267948966 : D_NAN;
        res[0] = v;
        for (prec = 1; prec < n; prec++)
            res[prec] = 0.0;
        if (!is_log && (v == 1.5707963267948966 || v == -1.5707963267948966))
        {
            /* +-pi/2 to full precision */
            arb_const_pi(t, 53 * n + 64);
            arb_mul_2exp_si(t, t, -1);
            if (v < 0)
                arb_neg(t, t);
            _dfloat_set_arf(res, n, &err, arb_midref(t));
        }
    }
    else
    {
        for (prec = 53 * n + 64; ; prec *= 2)
        {
            f(t, t, prec);
            if (arb_is_zero(t) || arb_rel_accuracy_bits(t) >= 53 * n + 40)
                break;
            _dfloat_get_arb(t, x, n, 0.0);
        }
        _dfloat_set_arf(res, n, &err, arb_midref(t));
    }
    arb_clear(t);
}

void
_dfloat_log_slow(double * res, int n, const double * x)
{
    _dfloat_unary_slow(res, n, x, arb_log, 1);
}

void
_dfloat_atan_slow(double * res, int n, const double * x)
{
    _dfloat_unary_slow(res, n, x, arb_atan, 0);
}

int
_dfloat_ball_log_slow(double * res, double * rrad, int n, const double * x, double rx)
{
    arb_t t;
    int status;
    arb_init(t);
    _dfloat_get_arb(t, x, n, rx);
    if (arb_is_positive(t))
    {
        arb_log(t, t, 53 * n + 40);
        _dfloat_set_arb(res, rrad, n, t);
        status = GR_SUCCESS;
    }
    else
        status = GR_UNABLE;
    arb_clear(t);
    return status;
}

void
_dfloat_ball_atan_slow(double * res, double * rrad, int n, const double * x, double rx)
{
    arb_t t;
    arb_init(t);
    _dfloat_get_arb(t, x, n, rx);
    arb_atan(t, t, 53 * n + 40);
    _dfloat_set_arb(res, rrad, n, t);
    arb_clear(t);
}

/* ---- expm1 and log1p: dispatchers and arb fallbacks ---- */

void
_dfloat_expm1(double * res, int n, const double * x)
{
    switch (n)
    {
        case 1: d1_expm1((d1_ptr) res, (d1_srcptr) x); break;
        case 2: d2_expm1((d2_ptr) res, (d2_srcptr) x); break;
        case 3: d3_expm1((d3_ptr) res, (d3_srcptr) x); break;
        default: d4_expm1((d4_ptr) res, (d4_srcptr) x); break;
    }
}

void
_dfloat_log1p(double * res, int n, const double * x)
{
    switch (n)
    {
        case 1: d1_log1p((d1_ptr) res, (d1_srcptr) x); break;
        case 2: d2_log1p((d2_ptr) res, (d2_srcptr) x); break;
        case 3: d3_log1p((d3_ptr) res, (d3_srcptr) x); break;
        default: d4_log1p((d4_ptr) res, (d4_srcptr) x); break;
    }
}

void
_dfloat_ball_expm1(double * res, double * rrad, int n, const double * x, double rx)
{
    double b[DFLOAT_MAX_N + 1], r[DFLOAT_MAX_N + 1];
    int i;
    for (i = 0; i < n; i++)
        b[i] = x[i];
    b[n] = rx;
    switch (n)
    {
        case 1: d1b_expm1((d1b_ptr) r, (d1b_srcptr) b); break;
        case 2: d2b_expm1((d2b_ptr) r, (d2b_srcptr) b); break;
        case 3: d3b_expm1((d3b_ptr) r, (d3b_srcptr) b); break;
        default: d4b_expm1((d4b_ptr) r, (d4b_srcptr) b); break;
    }
    for (i = 0; i < n; i++)
        res[i] = r[i];
    *rrad = r[n];
}

int
_dfloat_ball_log1p(double * res, double * rrad, int n, const double * x, double rx)
{
    double b[DFLOAT_MAX_N + 1], r[DFLOAT_MAX_N + 1];
    int i, status;
    for (i = 0; i < n; i++)
        b[i] = x[i];
    b[n] = rx;
    switch (n)
    {
        case 1: status = d1b_log1p((d1b_ptr) r, (d1b_srcptr) b); break;
        case 2: status = d2b_log1p((d2b_ptr) r, (d2b_srcptr) b); break;
        case 3: status = d3b_log1p((d3b_ptr) r, (d3b_srcptr) b); break;
        default: status = d4b_log1p((d4b_ptr) r, (d4b_srcptr) b); break;
    }
    if (status == GR_SUCCESS)
    {
        for (i = 0; i < n; i++)
            res[i] = r[i];
        *rrad = r[n];
    }
    return status;
}

/* to a relative accuracy of 53 n + 40 bits; special values as for the
   double functions (expm1: nan, inf, -1; log1p: nan below -1, -inf at
   -1, inf at inf) */
void
_dfloat_expm1_slow(double * res, int n, const double * x)
{
    arb_t t;
    slong prec;
    double err = 0.0;

    arb_init(t);
    _dfloat_get_arb(t, x, n, 0.0);
    if (!arb_is_finite(t))
    {
        res[0] = expm1(x[0]);
        for (prec = 1; prec < n; prec++)
            res[prec] = 0.0;
    }
    else
    {
        for (prec = 53 * n + 64; ; prec *= 2)
        {
            arb_expm1(t, t, prec);
            if (arb_is_zero(t) || arb_rel_accuracy_bits(t) >= 53 * n + 40)
                break;
            _dfloat_get_arb(t, x, n, 0.0);
        }
        _dfloat_set_arf(res, n, &err, arb_midref(t));
    }
    arb_clear(t);
}

void
_dfloat_log1p_slow(double * res, int n, const double * x)
{
    arb_t t, s;
    slong prec;
    double err = 0.0;

    arb_init(t);
    arb_init(s);
    _dfloat_get_arb(t, x, n, 0.0);
    arb_add_ui(s, t, 1, ARF_PREC_EXACT);
    if (!arb_is_finite(t) || !arb_is_positive(s))
    {
        double v = x[0];
        v = (arb_is_zero(s)) ? -D_INF : ((v > 0.0 && v == D_INF) ? D_INF : D_NAN);
        res[0] = v;
        for (prec = 1; prec < n; prec++)
            res[prec] = 0.0;
    }
    else
    {
        for (prec = 53 * n + 64; ; prec *= 2)
        {
            arb_log1p(t, t, prec);
            if (arb_is_zero(t) || arb_rel_accuracy_bits(t) >= 53 * n + 40)
                break;
            _dfloat_get_arb(t, x, n, 0.0);
        }
        _dfloat_set_arf(res, n, &err, arb_midref(t));
    }
    arb_clear(t);
    arb_clear(s);
}

int
_dfloat_ball_log1p_slow(double * res, double * rrad, int n, const double * x, double rx)
{
    arb_t t, s;
    int status;
    arb_init(t);
    arb_init(s);
    _dfloat_get_arb(t, x, n, rx);
    arb_add_ui(s, t, 1, 53 * n + 40);
    if (arb_is_positive(s))
    {
        arb_log1p(t, t, 53 * n + 40);
        _dfloat_set_arb(res, rrad, n, t);
        status = GR_SUCCESS;
    }
    else
        status = GR_UNABLE;
    arb_clear(t);
    arb_clear(s);
    return status;
}
