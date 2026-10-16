/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Ball operations of any magnitude, done by scaling the inputs to
   heads in [1, 2), running the fixed-N kernels where no underflow can
   occur, and scaling back.

   Scaling a double by a power of two is exact unless the result is
   subnormal (rounding error at most 2^-1075) or overflows.  Tail
   components that would become subnormal, or that fall below the
   TwoProd underflow guard 2^-480 relative to the unit head, are flushed
   to zero and their magnitude is added to the radius.  A radius that
   underflows is bumped up.  Scaling the result back may overflow (the
   result is then the whole real line) or underflow (every rounded
   component contributes 2^-1075, bounded by a 2^-1070 bump). */

#include "fp_contract.h"
#include <float.h>
#include "dfloat.h"
#include "double_extras.h"

#define TINY_THRESHOLD 0x1p-480

/* ---------------------------------------------------- renormalization */

/* The general renormalization of len terms to n components, the same
   algorithm as DF(_renorm_terms) in template.inc (which the fixed-size
   kernels inline) with n and len at run time: it is only used by slow
   paths and conversions, so one copy serves every n.  The results are
   bitwise identical to the templated version. */

/* VecSumErrBranch over the residuals t[1..len-1] of a VecSum */
static void
_renorm_branchy(double * res, int n, double * err, double * t, int M)
{
    int i, j;
    double f, eps, eps_new, e;

    eps = t[0];
    j = 0;

    for (i = 1; i < M; i++)
    {
        DFLOAT_TWO_SUM(f, eps_new, eps, t[i]);
        if (eps_new != 0.0)
        {
            res[j] = f;
            j++;
            if (j == n)
            {
                e = 0.0;
                f = res[n - 1];
                DFLOAT_TWO_SUM(f, eps, f, eps_new);
                e += fabs(eps);
                for (i = i + 1; i < M; i++)
                {
                    DFLOAT_TWO_SUM(f, eps, f, t[i]);
                    e += fabs(eps);
                }
                res[n - 1] = f;
                if (err != NULL)
                    *err += e;
                return;
            }
            eps = eps_new;
        }
        else
            eps = f;
    }

    res[j] = eps;
    for (j = j + 1; j < n; j++)
        res[j] = 0.0;
}

void _dfloat_renorm(double * res, int n, double * err, double * t, int M)
{
    int i, bad;
    double g, e, r[DFLOAT_MAX_N], sfx[DFLOAT_MAX_N + 1];

    if (M <= 0)
    {
        for (i = 0; i < n; i++)
            res[i] = 0.0;
        return;
    }

    for (i = M - 2; i >= 0; i--)
        DFLOAT_TWO_SUM(t[i], t[i + 1], t[i], t[i + 1]);

    if (M == 1)
    {
        res[0] = t[0];
        for (i = 1; i < n; i++)
            res[i] = 0.0;
        return;
    }

    /* suffix sums of the magnitudes: sfx[i] = |t[i+2]| + ... + |t[M-1]| */
    sfx[n] = 0.0;
    for (i = M - 1; i >= n + 2; i--)
        sfx[n] += fabs(t[i]);
    for (i = n - 1; i >= 0; i--)
        sfx[i] = sfx[i + 1] + ((i + 2 < M) ? fabs(t[i + 2]) : 0.0);

    g = t[1];
    bad = 0;
    for (i = 0; i < n; i++)
    {
        if (i == 0)
        {
            DFLOAT_TWO_SUM(r[0], g, t[0], t[1]);
        }
        else if (i + 1 < M)
        {
            DFLOAT_TWO_SUM(r[i], g, g, t[i + 1]);
        }
        else
        {
            r[i] = g;
            g = 0.0;
        }
        bad |= (g == 0.0) & (sfx[i] != 0.0);
    }

    if (FLINT_UNLIKELY(bad || r[0] == 0.0))
    {
        _renorm_branchy(res, n, err, t, M);
        return;
    }

    e = fabs(g);
    if (FLINT_UNLIKELY(sfx[n - 1] > 0x1p-40 * fabs(r[n - 1])))
    {
        for (i = n + 1; i < M; i++)
        {
            DFLOAT_TWO_SUM(r[n - 1], g, r[n - 1], t[i]);
            e += fabs(g);
        }
    }
    else
        e += sfx[n - 1];

    for (i = 0; i < n; i++)
        res[i] = r[i];

    if (err != NULL)
        *err += e;
}

/* ---------------------------------------------------------- dispatch */

void _dfloat_canonicalise(double * res, int n, const double * x)
{
    switch (n)
    {
        case 1: d1_canonicalise((d1_ptr) res, (d1_srcptr) x); break;
        case 2: d2_canonicalise((d2_ptr) res, (d2_srcptr) x); break;
        case 3: d3_canonicalise((d3_ptr) res, (d3_srcptr) x); break;
        default: d4_canonicalise((d4_ptr) res, (d4_srcptr) x); break;
    }
}

void _dfloat_add_n(double * res, double * err, int n, const double * x, const double * y)
{
    switch (n)
    {
        case 1: _d1_add_tracked(res, err, x, y); break;
        case 2: _d2_add_tracked(res, err, x, y); break;
        case 3: _d3_add_tracked(res, err, x, y); break;
        default: _d4_add_tracked(res, err, x, y); break;
    }
}

void _dfloat_sub_n(double * res, double * err, int n, const double * x, const double * y)
{
    switch (n)
    {
        case 1: _d1_sub_tracked(res, err, x, y); break;
        case 2: _d2_sub_tracked(res, err, x, y); break;
        case 3: _d3_sub_tracked(res, err, x, y); break;
        default: _d4_sub_tracked(res, err, x, y); break;
    }
}

void _dfloat_mul_n(double * res, double * err, int n, const double * x, const double * y)
{
    switch (n)
    {
        case 1: _d1_mul_tracked(res, err, x, y); break;
        case 2: _d2_mul_tracked(res, err, x, y); break;
        case 3: _d3_mul_tracked(res, err, x, y); break;
        default: _d4_mul_tracked(res, err, x, y); break;
    }
}

void _dfloat_sqr_n(double * res, double * err, int n, const double * x)
{
    switch (n)
    {
        case 1: _d1_sqr_tracked(res, err, x); break;
        case 2: _d2_sqr_tracked(res, err, x); break;
        case 3: _d3_sqr_tracked(res, err, x); break;
        default: _d4_sqr_tracked(res, err, x); break;
    }
}

void _dfloat_div_approx_n(double * q, int n, const double * x, const double * y)
{
    switch (n)
    {
        case 1: _d1_div_approx(q, x, y); break;
        case 2: _d2_div_approx(q, x, y); break;
        case 3: _d3_div_approx(q, x, y); break;
        default: _d4_div_approx(q, x, y); break;
    }
}

void _dfloat_sqrt_approx_n(double * s, int n, const double * x)
{
    switch (n)
    {
        case 1: _d1_sqrt_approx(s, x); break;
        case 2: _d2_sqrt_approx(s, x); break;
        case 3: _d3_sqrt_approx(s, x); break;
        default: _d4_sqrt_approx(s, x); break;
    }
}

void _dfloat_rsqrt_approx_n(double * r, int n, const double * x)
{
    switch (n)
    {
        case 1: _d1_rsqrt_approx(r, x); break;
        case 2: _d2_rsqrt_approx(r, x); break;
        case 3: _d3_rsqrt_approx(r, x); break;
        default: _d4_rsqrt_approx(r, x); break;
    }
}

/* ----------------------------------------------------------- scaling */

/* e with |x| 2^-e in [1, 2), for finite nonzero x */
static int
_exponent(double x)
{
    int e;
    frexp(x, &e);
    return e - 1;
}

/* res = x 2^-e with the head in [1, 2); tails that would be inexact
   or below the guard are flushed into *err (in the scaled units).
   The radius is scaled too. */
static void
_scale_in(double * res, double * rres, const double * x, double rx, int n, int e)
{
    double err = 0.0, v;
    int i;

    res[0] = ldexp(x[0], -e);
    for (i = 1; i < n; i++)
    {
        v = ldexp(x[i], -e);
        if (x[i] == 0.0)
        {
            res[i] = 0.0;
        }
        else if (fabs(v) < TINY_THRESHOLD)
        {
            /* flushed; v may have been rounded (even to zero) if it
               is subnormal, by at most DBL_MIN */
            res[i] = 0.0;
            err += fabs(v) + DBL_MIN;
        }
        else
            res[i] = v;
    }

    v = ldexp(rx, -e);
    if (rx != 0.0 && v < DFLOAT_RAD_UNDERFLOW_THRESHOLD)
        v += DFLOAT_RAD_UNDERFLOW_BUMP;
    *rres = v + err;
}

/* res = q 2^e, rad = rq 2^e (+ bumps), or the whole line on overflow */
static void
_scale_out(double * res, double * rres, const double * q, double rq, int n, int e)
{
    double err = 0.0, v, r;
    int i;

    for (i = 0; i < n; i++)
    {
        v = ldexp(q[i], e);
        res[i] = v;
        if (q[i] != 0.0 && fabs(v) < DBL_MIN)
            err += DFLOAT_RAD_UNDERFLOW_BUMP;
    }

    r = ldexp(rq, e);
    if (rq != 0.0 && r < DFLOAT_RAD_UNDERFLOW_THRESHOLD)
        r += DFLOAT_RAD_UNDERFLOW_BUMP;
    r = _dfloat_rad_finish(r + err);

    if (!(_dfloat_abs_sum(res, n) + r <= DBL_MAX))
    {
        /* the whole line, as [0 +/- inf] */
        for (i = 0; i < n; i++)
            res[i] = 0.0;
        r = D_INF;
    }

    *rres = r;
}

static void
_set_zero(double * res, int n)
{
    int i;
    for (i = 0; i < n; i++)
        res[i] = 0.0;
}

/* -------------------------------------------------------------- mul */

void
_dfloat_ball_mul_scaled(double * res, double * rrad, int n,
    const double * x, double rx, const double * y, double ry)
{
    double xs[DFLOAT_MAX_N], ys[DFLOAT_MAX_N], zs[DFLOAT_MAX_N];
    double rxs, rys, rad, err = 0.0, ax, ay;
    int ex, ey;

    /* the whole line times anything but an exact zero */
    if (_dfloat_ball_is_whole(x, n, rx) || _dfloat_ball_is_whole(y, n, ry))
    {
        int zero = (!_dfloat_ball_is_whole(x, n, rx) && rx == 0.0 && _dfloat_abs_sum(x, n) == 0.0)
                || (!_dfloat_ball_is_whole(y, n, ry) && ry == 0.0 && _dfloat_abs_sum(y, n) == 0.0);
        _set_zero(res, n);
        *rrad = zero ? 0.0 : D_INF;
        return;
    }

    if (x[0] == 0.0 || y[0] == 0.0)
    {
        /* one midpoint is exactly zero */
        ax = _dfloat_abs_sum(x, n);
        ay = _dfloat_abs_sum(y, n);
        rad = _dfloat_rad_mul(ax, ry) + _dfloat_rad_mul(ay, rx);
        rad += _dfloat_rad_mul(rx, ry);
        _set_zero(res, n);
        *rrad = _dfloat_rad_finish(rad);
        if (!(*rrad <= DBL_MAX))
            *rrad = D_INF;
        return;
    }

    if (!(_dfloat_abs_sum(x, n) <= DBL_MAX) || !(_dfloat_abs_sum(y, n) <= DBL_MAX))
    {
        _set_zero(res, n);
        *rrad = D_INF;
        return;
    }

    ex = _exponent(x[0]);
    ey = _exponent(y[0]);
    _scale_in(xs, &rxs, x, rx, n, ex);
    _scale_in(ys, &rys, y, ry, n, ey);

    _dfloat_mul_n(zs, &err, n, xs, ys);
    ax = _dfloat_abs_sum(xs, n);
    ay = _dfloat_abs_sum(ys, n);
    rad = _dfloat_rad_mul(ax, rys) + _dfloat_rad_mul(ay, rxs);
    rad += _dfloat_rad_mul(rxs, rys);
    rad += err;

    /* ex + ey is at most 2046 in absolute value: ldexp handles it */
    _scale_out(res, rrad, zs, rad, n, ex + ey);
}

/* -------------------------------------------------------------- div */

/* Division: q = x / y with y not containing zero (checked by the
   caller).  After scaling, q is computed by long division, small
   tails of q are flushed (q is our choice; the residual accounts for
   it), and the residual x - q y is enclosed with the tracked kernels:
       x - q y  in  D +/- (Derr + Perr),   P = q y, D = x - P.
   Then for x' in x +/- rx and y' in y +/- ry,
       |x'/y' - q| = |x' - q y'| / |y'|
                  <= (|x - q y| + rx + |q| ry) / (|y| - ry). */
void
_dfloat_ball_div(double * res, double * rrad, int n,
    const double * x, double rx, const double * y, double ry)
{
    double xs[DFLOAT_MAX_N], ys[DFLOAT_MAX_N], q[DFLOAT_MAX_N];
    double P[DFLOAT_MAX_N], D[DFLOAT_MAX_N];
    double rxs, rys, Perr = 0.0, Derr = 0.0, E, ylo, aq, rad;
    int ex, ey, i;

    /* (y was checked by the caller) */
    if (_dfloat_ball_is_whole(x, n, rx))
    {
        _set_zero(res, n);
        *rrad = D_INF;
        return;
    }

    if (x[0] == 0.0 && rx == 0.0)
    {
        _set_zero(res, n);
        *rrad = 0.0;
        return;
    }

    if (!(_dfloat_abs_sum(x, n) <= DBL_MAX) || !(_dfloat_abs_sum(y, n) <= DBL_MAX)
        || y[0] == 0.0 || !(rx <= DBL_MAX))
    {
        _set_zero(res, n);
        *rrad = D_INF;
        return;
    }

    ey = _exponent(y[0]);
    _scale_in(ys, &rys, y, ry, n, ey);

    /* a lower bound on |y| after flushing must still exclude zero */
    ylo = _dfloat_lower(fabs(ys[0]) - _dfloat_rad_finish(_dfloat_abs_sum(ys + 1, n - 1) + rys));
    if (!(ylo > 0.0))
    {
        _set_zero(res, n);
        *rrad = D_INF;
        return;
    }

    if (x[0] == 0.0)
    {
        /* midpoint zero, nonzero radius: |x'/y'| <= rx / (|y| - ry) */
        _set_zero(res, n);
        rad = _dfloat_rad_div(rx, ylo);
        _scale_out(res, rrad, res, rad, n, -ey);
        return;
    }

    ex = _exponent(x[0]);
    _scale_in(xs, &rxs, x, rx, n, ex);

    _dfloat_div_approx_n(q, n, xs, ys);
    for (i = 1; i < n; i++)
        if (fabs(q[i]) < TINY_THRESHOLD)
            q[i] = 0.0;
    if (!(_dfloat_abs_sum(q, n) <= DBL_MAX))
    {
        _set_zero(res, n);
        *rrad = D_INF;
        return;
    }

    _dfloat_mul_n(P, &Perr, n, q, ys);
    _dfloat_sub_n(D, &Derr, n, xs, P);
    E = _dfloat_abs_sum(D, n) + Derr + Perr;

    aq = _dfloat_abs_sum(q, n);
    rad = E + rxs + _dfloat_rad_mul(aq, rys);
    rad = _dfloat_rad_div(rad, ylo);

    _scale_out(res, rrad, q, rad, n, ex - ey);
}

/* ------------------------------------------------------------- sqrt */

/* Square root.  x must not contain negative numbers: GR_DOMAIN if the
   ball is entirely negative, GR_UNABLE if it straddles zero (as arb).
   For x' in [x +/- rx], all >= 0:
       |sqrt(x') - s| <= |x' - s^2| / (sqrt(x') + s)
                      <= (|x - s^2| + rx) / s
   where s is the computed root.  Scaling is by an even power. */
int
_dfloat_ball_sqrt(double * res, double * rrad, int n,
    const double * x, double rx)
{
    double xs[DFLOAT_MAX_N], s[DFLOAT_MAX_N], P[DFLOAT_MAX_N], D[DFLOAT_MAX_N];
    double rxs, Perr = 0.0, Derr = 0.0, E, slo, rad, xlo, xhi;
    int ex, i;

    /* an infinite radius means the ball contains negative numbers */
    if (!(_dfloat_abs_sum(x, n) <= DBL_MAX) || !(rx <= DBL_MAX))
        return GR_UNABLE;

    /* sign analysis on the unscaled ball */
    xhi = x[0] + _dfloat_rad_finish(_dfloat_abs_sum(x + 1, n - 1) + rx);
    xlo = x[0] - _dfloat_rad_finish(_dfloat_abs_sum(x + 1, n - 1) + rx);

    if (xhi < 0.0)
        return GR_DOMAIN;
    if (x[0] == 0.0 && rx == 0.0)
    {
        _set_zero(res, n);
        *rrad = 0.0;
        return GR_SUCCESS;
    }
    if (!(xlo >= 0.0))
        return GR_UNABLE;   /* contains negative numbers */

    ex = _exponent(x[0]);
    if (ex & 1)
        ex--;
    _scale_in(xs, &rxs, x, rx, n, ex);

    _dfloat_sqrt_approx_n(s, n, xs);
    for (i = 1; i < n; i++)
        if (fabs(s[i]) < TINY_THRESHOLD)
            s[i] = 0.0;

    _dfloat_sqr_n(P, &Perr, n, s);
    _dfloat_sub_n(D, &Derr, n, xs, P);
    E = _dfloat_abs_sum(D, n) + Derr + Perr;

    slo = _dfloat_lower(fabs(s[0]) - _dfloat_rad_finish(_dfloat_abs_sum(s + 1, n - 1)));
    if (!(slo > 0.0))
    {
        _set_zero(res, n);
        *rrad = D_INF;
        return GR_SUCCESS;
    }

    rad = _dfloat_rad_div(E + rxs, slo);
    _scale_out(res, rrad, s, rad, n, ex / 2);
    return GR_SUCCESS;
}

/* Reciprocal square root.  x must be positive: GR_DOMAIN if the ball
   is entirely nonpositive, GR_UNABLE if it contains zero.  For x' in
   [x +/- rx], all > 0, and r the computed reciprocal root:
       1/sqrt(x') = r (1 - e')^(-1/2), e' = 1 - x' r^2,
       |1/sqrt(x') - r| <= |r| (|e'|/2) (1 + |e'|)   for |e'| <= 1/4,
       |e'| <= |1 - x r^2| + rx r^2.
   Scaling is by an even power. */
int
_dfloat_ball_rsqrt(double * res, double * rrad, int n,
    const double * x, double rx)
{
    double xs[DFLOAT_MAX_N], r[DFLOAT_MAX_N], P[DFLOAT_MAX_N], Q[DFLOAT_MAX_N], D[DFLOAT_MAX_N], one[DFLOAT_MAX_N];
    double rxs, Perr = 0.0, Qerr = 0.0, Derr = 0.0, E, rhi, xhi, rad, xlo;
    int ex, i;

    if (!(_dfloat_abs_sum(x, n) <= DBL_MAX) || !(rx <= DBL_MAX))
        return GR_UNABLE;

    xhi = x[0] + _dfloat_rad_finish(_dfloat_abs_sum(x + 1, n - 1) + rx);
    xlo = x[0] - _dfloat_rad_finish(_dfloat_abs_sum(x + 1, n - 1) + rx);

    if (xhi <= 0.0)
        return GR_DOMAIN;
    if (!(xlo > 0.0))
        return GR_UNABLE;

    ex = _exponent(x[0]);
    if (ex & 1)
        ex--;
    _scale_in(xs, &rxs, x, rx, n, ex);

    _dfloat_rsqrt_approx_n(r, n, xs);
    for (i = 1; i < n; i++)
        if (fabs(r[i]) < TINY_THRESHOLD)
            r[i] = 0.0;

    _dfloat_sqr_n(P, &Perr, n, r);
    _dfloat_mul_n(Q, &Qerr, n, xs, P);
    one[0] = 1.0;
    for (i = 1; i < n; i++)
        one[i] = 0.0;
    _dfloat_sub_n(D, &Derr, n, one, Q);
    rhi = _dfloat_rad_finish(_dfloat_abs_sum(r, n));
    xhi = xs[0] + _dfloat_rad_finish(_dfloat_abs_sum(xs + 1, n - 1) + rxs);
    E = _dfloat_abs_sum(D, n) + Derr + Qerr + _dfloat_rad_mul(xhi, Perr);
    E = _dfloat_rad_finish(E + _dfloat_rad_mul(rxs, _dfloat_rad_mul(rhi, rhi)));
    if (!(E <= 0.25))
    {
        _set_zero(res, n);
        *rrad = D_INF;
        return GR_SUCCESS;
    }
    rad = _dfloat_rad_mul(rhi, 0.5 * E * (1.0 + E));
    _scale_out(res, rrad, r, rad, n, -ex / 2);
    return GR_SUCCESS;
}
