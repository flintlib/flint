/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* The complex functions that are not performance-critical (division,
   square roots, the elementary functions, conversions and their gr
   methods), implemented once for all N through the real operations of
   the format (_dfloat_ops_struct). The arithmetic, the vector
   operations and the dot products are compiled per format
   (ctemplate.inc).

   Plain and ball formats share the code: the plain real operations
   never fail, and the checks that only concern balls (containment of
   zero, the whole line, the enclosures used where a formula breaks
   down) are made under ops->ball. Every ball function is rigorous by
   composition: each real operation returns an enclosure, and every
   formula is an identity at every point of the input balls (the
   exceptions, at branch cuts and singularities, are handled
   explicitly). Outputs may alias inputs. */

#include "fp_contract.h"
#include <float.h>
#include <math.h>
#include <string.h>
#include "dfloat.h"
#include "double_extras.h"
#include "fmpz.h"
#include "arb.h"
#include "acb.h"
#include "gr.h"
#include "gr_generic.h"

/* the doubles of a real element, and the parts of a complex element */
#define PSZ(ops) ((ops)->n + (ops)->ball)
#define IM(ops, x) ((x) + PSZ(ops))

/* temporaries: a real element, a complex element */
typedef double part_t[DFLOAT_MAX_N + 1];
typedef double cplx_t[2 * (DFLOAT_MAX_N + 1)];

#define AUTO DFLOAT_MODE_AUTO

/* everything here is compiled for size: the work is in the real
   functions called, and a call costs nothing in relative terms */
#if defined(__GNUC__)
#define COLD __attribute__((cold))
#else
#define COLD
#endif

static const _dfloat_ops_struct * const _dfloat_ops_tab[2][DFLOAT_MAX_N + 1] = {
    {NULL, &_d1_ops, &_d2_ops, &_d3_ops, &_d4_ops},
    {NULL, &_d1b_ops, &_d2b_ops, &_d3b_ops, &_d4b_ops},
};

_dfloat_ops_t
_dfloat_ops(int n, int ball)
{
    return _dfloat_ops_tab[ball != 0][n];
}

/* ---- real parts ---- */

static COLD void
part_set(_dfloat_ops_t ops, double * res, const double * x)
{
    int i;
    if (res != x)
        for (i = 0; i < PSZ(ops); i++)
            res[i] = x[i];
}

static COLD void
cplx_set(_dfloat_ops_t ops, double * res, const double * x)
{
    int i;
    if (res != x)
        for (i = 0; i < 2 * PSZ(ops); i++)
            res[i] = x[i];
}

/* an exact zero (plain: the head is zero only for zero) */
static COLD int
part_is_zero(_dfloat_ops_t ops, const double * x)
{
    if (ops->ball)
        return x[ops->n] == 0.0 && _dfloat_abs_sum(x, ops->n) == 0.0;
    return x[0] == 0.0;
}

/* (plain: an exact zero) */
static COLD int
part_contains_zero(_dfloat_ops_t ops, const double * x)
{
    return ops->ball ? ops->contains_zero(x) : (x[0] == 0.0);
}

/* a ball part that is the whole line */
static COLD int
part_is_W(_dfloat_ops_t ops, const double * x)
{
    return ops->ball && _dfloat_ball_is_whole(x, ops->n, x[ops->n]);
}

/* [0 +/- b] (b = inf: the whole line) */
static COLD void
part_zero_pm(_dfloat_ops_t ops, double * res, double b)
{
    int i;
    for (i = 0; i < ops->n; i++)
        res[i] = 0.0;
    res[ops->n] = b;
}

/* an upper bound on |x| (+inf for the whole line) */
static COLD double
part_upper(_dfloat_ops_t ops, const double * x)
{
    double v = _dfloat_rad_finish(_dfloat_abs_sum(x, ops->n) + (ops->ball ? x[ops->n] : 0.0));
    return (v <= DBL_MAX) ? v : D_INF;
}

/* a lower bound on |x| (0 if the ball contains zero) */
static COLD double
part_lower(_dfloat_ops_t ops, const double * x)
{
    if (part_is_W(ops, x))
        return 0.0;
    return _dfloat_lower(fabs(x[0]) - _dfloat_rad_finish(_dfloat_abs_sum(x + 1, ops->n - 1)
        + (ops->ball ? x[ops->n] : 0.0)));
}

/* the ball [lo, hi] for 0 <= lo <= hi (the whole line if hi = inf) */
static COLD void
part_interval(_dfloat_ops_t ops, double * res, double lo, double hi)
{
    double m;
    if (!(hi <= DBL_MAX))
    {
        part_zero_pm(ops, res, D_INF);
        return;
    }
    m = 0.5 * lo + 0.5 * hi;
    part_zero_pm(ops, res, _dfloat_rad_finish(FLINT_MAX(hi - m, m - lo)));
    res[0] = m;
}

/* signed bounds of a real element (+-inf for the whole line) */
static COLD double
part_lo(_dfloat_ops_t ops, const double * x)
{
    double r, v;
    if (part_is_W(ops, x))
        return -D_INF;
    r = _dfloat_rad_finish(_dfloat_abs_sum(x + 1, ops->n - 1) + (ops->ball ? x[ops->n] : 0.0));
    v = x[0] - r;
    return v - (fabs(v) * 0x1p-50 + 0x1p-1074);
}

static COLD double
part_hi(_dfloat_ops_t ops, const double * x)
{
    double r, v;
    if (part_is_W(ops, x))
        return D_INF;
    r = _dfloat_rad_finish(_dfloat_abs_sum(x + 1, ops->n - 1) + (ops->ball ? x[ops->n] : 0.0));
    v = x[0] + r;
    return v + (fabs(v) * 0x1p-50 + 0x1p-1074);
}

/* the ball [lo, hi] for any finite doubles lo <= hi */
static COLD void
part_interval_signed(_dfloat_ops_t ops, double * res, double lo, double hi)
{
    double m;
    if (!(lo >= -DBL_MAX && hi <= DBL_MAX))
    {
        part_zero_pm(ops, res, D_INF);
        return;
    }
    m = 0.5 * lo + 0.5 * hi;
    part_zero_pm(ops, res, _dfloat_rad_finish(FLINT_MAX(hi - m, m - lo)));
    res[0] = m;
}

/* a lower bound on |z| over the rectangle (0 if unknown) */
static COLD double
cplx_lower_abs(_dfloat_ops_t ops, const double * x)
{
    double la = part_lower(ops, x), lb = part_lower(ops, IM(ops, x));
    /* (max(la, lb) where the squares could overflow or underflow) */
    if (!(la >= 0x1p-500 && la <= 0x1p500 && lb >= 0x1p-500 && lb <= 0x1p500))
        return FLINT_MAX(la, lb);
    return _dfloat_lower(sqrt(_dfloat_lower(la * la + lb * lb)));
}

/* an upper bound on |z| over the rectangle (+inf if unbounded) */
static COLD double
cplx_upper_abs(_dfloat_ops_t ops, const double * x)
{
    double u = _dfloat_rad_finish(part_upper(ops, x) + part_upper(ops, IM(ops, x)));
    return (u <= DBL_MAX) ? u : D_INF;
}

/* both parts [0 +/- b] (the enclosure of a set of modulus <= b) */
static COLD void
cplx_zero_pm(_dfloat_ops_t ops, double * res, double b)
{
    if (!(b <= DBL_MAX))
        b = D_INF;
    part_zero_pm(ops, res, b);
    part_zero_pm(ops, IM(ops, res), b);
}

/* a constant from a 5-term table row (balls: the rest and the table
   residual in the radius) */
static COLD void
part_table(_dfloat_ops_t ops, double * res, const double * row)
{
    double r = _dfloat_exp_tab_err;
    int i;
    for (i = 0; i < ops->n; i++)
        res[i] = row[i];
    if (ops->ball)
    {
        for (i = ops->n; i < 5; i++)
            r += fabs(row[i]);
        res[ops->n] = _dfloat_rad_finish(r);
    }
}

#define ALL_GT(x, c) ops->cmp_all_d(x, c, 1, 1)
#define ALL_GE(x, c) ops->cmp_all_d(x, c, 1, 0)
#define ALL_LT(x, c) ops->cmp_all_d(x, c, -1, 1)
#define ALL_LE(x, c) ops->cmp_all_d(x, c, -1, 0)

static COLD int
cplx_is_zero(_dfloat_ops_t ops, const double * x)
{
    return part_is_zero(ops, x) && part_is_zero(ops, IM(ops, x));
}

static COLD int
cplx_contains_zero(_dfloat_ops_t ops, const double * x)
{
    return part_contains_zero(ops, x) && part_contains_zero(ops, IM(ops, x));
}

static COLD void
cplx_one(_dfloat_ops_t ops, double * res)
{
    ops->one(res);
    ops->zero(IM(ops, res));
}

/* the exponent e of a scaling by 2^-e that brings the larger head to
   [1/2, 1), or 0 if the heads are within [2^-500, 2^500] (or not
   finite), where no square overflows or underflows */
static COLD int
scale_exp(double a, double b)
{
    double m = FLINT_MAX(fabs(a), fabs(b));
    int e;
    if (m >= 0x1p-500 && m <= 0x1p500)
        return 0;
    if (!(m <= DBL_MAX) || m == 0.0)
        return 0;
    frexp(m, &e);
    return e;
}

/* |z|^2 = a^2 + b^2; for balls, the lower bound is raised to (min |a|)^2
   + (min |b|)^2 where the squares of wide balls give less, so that a
   ball z away from zero gives |z|^2 away from zero */
static COLD void
abs2(_dfloat_ops_t ops, double * res, const double * a, const double * b, int mode)
{
    part_t s, t;
    double la, lb, lo, hi, low;
    int n = ops->n;
    ops->sqr(s, a, mode);
    ops->sqr(t, b, mode);
    ops->add(res, s, t);
    if (ops->ball && res[n] > 0x1p-30 * fabs(res[0]) && !part_is_W(ops, res))
    {
        la = part_lower(ops, a);
        lb = part_lower(ops, b);
        lo = _dfloat_lower(la * la + lb * lb);
        low = res[0] - _dfloat_rad_finish(_dfloat_abs_sum(res + 1, n - 1) + res[n]);
        if (lo > low)
        {
            hi = part_upper(ops, res);
            if (hi >= lo)
                part_interval(ops, res, lo, hi);
        }
    }
}

/* the square root of a real element enclosing a nonnegative quantity
   (balls: [0, sqrt(upper bound)] where the ball square root fails) */
static COLD void
sqrt_nonneg(_dfloat_ops_t ops, double * res, const double * x, int mode)
{
    double hi;
    if (ops->sqrt(res, x, mode) == GR_SUCCESS || !ops->ball)
        return;
    hi = part_upper(ops, x);
    part_interval(ops, res, 0.0, (hi <= DBL_MAX) ? _dfloat_rad_finish(sqrt(hi)) : D_INF);
}

/* sinh and cosh of a real element from one exponential: with e =
   expm1 |x|, sinh |x| = e (e + 2) / (2 (e + 1)) and cosh x = 1 + e^2 /
   (2 (e + 1)) for |x| < 1 (no cancellation), and (E -+ 1/E) / 2, E =
   exp |x|, beyond. The formulas are identities for all real x, and
   sinh is odd, cosh even, so that the sign of the head can be taken
   out even for a ball containing zero. */
static COLD void
real_sinh_cosh(_dfloat_ops_t ops, double * sh, double * ch, const double * x)
{
    part_t a, e, t, u, v;
    int neg = (x[0] < 0.0), status;
    if (part_is_W(ops, x))
    {
        part_zero_pm(ops, sh, D_INF);
        part_zero_pm(ops, ch, D_INF);
        return;
    }
    ops->abs(a, x);
    if (a[0] < 1.0)
    {
        ops->expm1(e, a);
        ops->one(t);
        ops->add(u, e, t);                      /* e + 1 = exp(a) > 0 */
        ops->mul_2exp_si(u, u, 1);
        ops->set_d(t, 2.0);
        ops->add(v, e, t);
        ops->mul(v, v, e, AUTO);
        if (ops->div(sh, v, u, AUTO) != GR_SUCCESS)
            part_zero_pm(ops, sh, D_INF);
        ops->sqr(v, e, AUTO);
        if (ops->div(v, v, u, AUTO) != GR_SUCCESS)
            part_zero_pm(ops, ch, D_INF);
        else
        {
            ops->one(t);
            ops->add(ch, v, t);
        }
    }
    else
    {
        ops->exp(e, a);
        ops->one(t);
        status = ops->div(t, t, e, AUTO);
        if (status != GR_SUCCESS)
        {
            part_zero_pm(ops, sh, D_INF);
            part_zero_pm(ops, ch, D_INF);
            return;
        }
        ops->sub(u, e, t);
        ops->add(v, e, t);
        ops->mul_2exp_si(sh, u, -1);
        ops->mul_2exp_si(ch, v, -1);
    }
    if (neg)
        ops->neg(sh, sh);
}

/* ------------------------------------------------------------------ */
/* Arithmetic                                                          */
/* ------------------------------------------------------------------ */

/* where the formula fails for a ball y away from zero (a wide y, whose
   |y|^2 cannot be bounded away from zero as a ball, or a y with a part
   that is the whole line): |x / y| <= max |x| / min |y| */
static COLD int
div_crude(_dfloat_ops_t ops, double * res, const double * x, const double * y)
{
    double l = cplx_lower_abs(ops, y), u;
    if (!(l > 0.0))
        return GR_UNABLE;
    u = cplx_upper_abs(ops, x);
    cplx_zero_pm(ops, res, (u <= DBL_MAX) ? _dfloat_rad_finish(_dfloat_rad_div(u, l)) : D_INF);
    return GR_SUCCESS;
}

/* x / y: a real or imaginary y by two real divisions, otherwise x
   conj(y) / |y|^2 with y scaled by a power of two so that |y|^2 can
   neither overflow nor underflow. Balls: GR_DOMAIN for y = 0, GR_UNABLE
   for a y containing zero (or whose |y|^2 cannot be bounded away from
   zero). */
COLD int
_dfloat_complex_div(_dfloat_ops_t ops, double * res, const double * x, const double * y, int mode)
{
    part_t c, d, s, t, u, v;
    const double * a = x, * b = IM(ops, x);
    int e, st1, st2, status = GR_SUCCESS;

    if (cplx_is_zero(ops, y))
    {
        if (ops->ball)
            return GR_DOMAIN;
        status = GR_DOMAIN;
    }

    if (part_is_zero(ops, IM(ops, y)))
    {
        part_set(ops, c, y);
        st1 = ops->div(t, a, c, mode);
        st2 = ops->div(u, b, c, mode);
        if ((st1 | st2) != GR_SUCCESS)
            return div_crude(ops, res, x, y);
        part_set(ops, res, t);
        part_set(ops, IM(ops, res), u);
        return status;
    }

    if (part_is_zero(ops, y))
    {
        /* (a + b i) / (d i) = b / d - (a / d) i */
        part_set(ops, d, IM(ops, y));
        st1 = ops->div(t, b, d, mode);
        st2 = ops->div(u, a, d, mode);
        if ((st1 | st2) != GR_SUCCESS)
            return div_crude(ops, res, x, y);
        part_set(ops, res, t);
        ops->neg(IM(ops, res), u);
        return status;
    }

    if (ops->ball && cplx_contains_zero(ops, y))
        return GR_UNABLE;
    if (ops->ball && (part_is_W(ops, y) || part_is_W(ops, IM(ops, y))))
        return div_crude(ops, res, x, y);

    e = scale_exp(y[0], IM(ops, y)[0]);
    ops->mul_2exp_si(c, y, -e);
    ops->mul_2exp_si(d, IM(ops, y), -e);
    abs2(ops, s, c, d, mode);
    /* (a c + b d) + (b c - a d) i */
    ops->mul(t, a, c, mode);
    ops->mul(u, b, d, mode);
    ops->add(v, t, u);
    ops->mul(t, b, c, mode);
    ops->mul(u, a, d, mode);
    ops->sub(u, t, u);
    /* (two statements: the second division overwrites v, which the
       first reads, and the order of evaluation of the operands of | is
       unspecified; MSVC evaluated the second one first) */
    st1 = ops->div(t, v, s, mode);
    st2 = ops->div(v, u, s, mode);
    if ((st1 | st2) != GR_SUCCESS)
        return div_crude(ops, res, x, y);
    if (e != 0)
    {
        ops->mul_2exp_si(t, t, -e);
        ops->mul_2exp_si(v, v, -e);
    }
    part_set(ops, res, t);
    part_set(ops, IM(ops, res), v);
    return status;
}

COLD int
_dfloat_complex_inv(_dfloat_ops_t ops, double * res, const double * x, int mode)
{
    cplx_t one;
    cplx_one(ops, one);
    return _dfloat_complex_div(ops, res, one, x, mode);
}

/* |x| = sqrt(a^2 + b^2), scaled */
COLD void
_dfloat_complex_abs(_dfloat_ops_t ops, double * res, const double * x)
{
    part_t a, b, s;
    int e;
    if (part_is_zero(ops, IM(ops, x)))
    {
        ops->abs(res, x);
        return;
    }
    if (part_is_zero(ops, x))
    {
        ops->abs(res, IM(ops, x));
        return;
    }
    if (part_is_W(ops, x) || part_is_W(ops, IM(ops, x)))
    {
        part_zero_pm(ops, res, D_INF);
        return;
    }
    e = scale_exp(x[0], IM(ops, x)[0]);
    ops->mul_2exp_si(a, x, -e);
    ops->mul_2exp_si(b, IM(ops, x), -e);
    abs2(ops, s, a, b, AUTO);
    sqrt_nonneg(ops, res, s, AUTO);
    if (e != 0)
        ops->mul_2exp_si(res, res, e);
}

COLD void
_dfloat_complex_arg(_dfloat_ops_t ops, double * res, const double * x)
{
    ops->atan2(res, IM(ops, x), x);
}

/* an enclosure of sqrt over a ball: [0, s] + [-s, s] i with s =
   sqrt(an upper bound of |x|) */
static COLD void
sqrt_crude(_dfloat_ops_t ops, double * res, const double * x)
{
    double u = _dfloat_rad_finish(part_upper(ops, x) + part_upper(ops, IM(ops, x)));
    u = (u <= DBL_MAX) ? _dfloat_rad_finish(sqrt(u)) : D_INF;
    part_interval(ops, res, 0.0, u);
    part_zero_pm(ops, IM(ops, res), u);
}

/* The principal square root t + u i (t >= 0, and defined everywhere):
   with m = |x|, t = sqrt((m + a)/2) and u = b / (2 t) for a >= 0, and
   |u| = sqrt((m - a)/2) with the sign of b and t = b / (2 u) for a < 0
   (neither form cancels), scaled by an even power of two. For balls
   these hold at every point as long as the divisor does not vanish,
   which fails only for balls touching the branch cut (-inf, 0], where
   the crude enclosure is used; for a < 0 and b containing zero (a ball
   crossing the cut), u is [0 +/- |u|] and t = |b| / (2 |u|). */
COLD void
_dfloat_complex_sqrt(_dfloat_ops_t ops, double * res, const double * x, int mode)
{
    part_t a, b, s, m, t, u, w;
    int e;

    if (part_is_zero(ops, IM(ops, x)))
    {
        int pos, neg;
        if (ops->ball)
        {
            pos = ALL_GE(x, 0.0);
            neg = !pos && ALL_LE(x, 0.0);
        }
        else
        {
            pos = (x[0] >= 0.0 || x[0] != x[0]);
            neg = !pos;
        }
        if (pos && ops->sqrt(t, x, mode) == GR_SUCCESS)
        {
            part_set(ops, res, t);
            ops->zero(IM(ops, res));
            return;
        }
        if (neg)
        {
            ops->neg(u, x);
            if (ops->sqrt(t, u, mode) == GR_SUCCESS)
            {
                part_set(ops, IM(ops, res), t);
                ops->zero(res);
                return;
            }
        }
        sqrt_crude(ops, res, x);
        return;
    }

    if (ops->ball && (part_is_W(ops, x) || part_is_W(ops, IM(ops, x)) || cplx_contains_zero(ops, x)))
    {
        sqrt_crude(ops, res, x);
        return;
    }

    e = scale_exp(x[0], IM(ops, x)[0]);
    e += (e & 1);
    ops->mul_2exp_si(a, x, -e);
    ops->mul_2exp_si(b, IM(ops, x), -e);
    abs2(ops, s, a, b, mode);
    sqrt_nonneg(ops, m, s, mode);
    if (a[0] >= 0.0)
    {
        ops->add(w, m, a);
        ops->mul_2exp_si(w, w, -1);
        sqrt_nonneg(ops, t, w, mode);
        ops->mul_2exp_si(w, t, 1);
        if (ops->div(u, b, w, mode) != GR_SUCCESS)
        {
            sqrt_crude(ops, res, x);
            return;
        }
    }
    else
    {
        ops->sub(w, m, a);
        ops->mul_2exp_si(w, w, -1);
        sqrt_nonneg(ops, u, w, mode);
        ops->mul_2exp_si(w, u, 1);
        if (part_contains_zero(ops, b))
        {
            ops->abs(s, b);
            if (ops->div(t, s, w, mode) != GR_SUCCESS)
            {
                sqrt_crude(ops, res, x);
                return;
            }
            part_zero_pm(ops, u, part_upper(ops, u));
        }
        else
        {
            if (b[0] < 0.0)
            {
                ops->neg(u, u);
                ops->neg(w, w);
            }
            if (ops->div(t, b, w, mode) != GR_SUCCESS)
            {
                sqrt_crude(ops, res, x);
                return;
            }
        }
    }
    ops->mul_2exp_si(res, t, e / 2);
    ops->mul_2exp_si(IM(ops, res), u, e / 2);
}

COLD int
_dfloat_complex_rsqrt(_dfloat_ops_t ops, double * res, const double * x, int mode)
{
    cplx_t t;
    int status = GR_SUCCESS;
    if (cplx_is_zero(ops, x))
    {
        if (ops->ball)
            return GR_DOMAIN;
        status = GR_DOMAIN;
    }
    else if (ops->ball && cplx_contains_zero(ops, x))
        return GR_UNABLE;
    _dfloat_complex_sqrt(ops, t, x, mode);
    if (_dfloat_complex_inv(ops, res, t, mode) != GR_SUCCESS && ops->ball)
    {
        /* |1/sqrt(x)| <= 1/sqrt(min |x|) */
        double l = cplx_lower_abs(ops, x);
        if (!(l > 0.0))
            return GR_UNABLE;
        l = _dfloat_lower(sqrt(l));
        if (!(l > 0.0))
            return GR_UNABLE;
        cplx_zero_pm(ops, res, _dfloat_rad_finish(_dfloat_rad_div(1.0, l)));
    }
    return status;
}

/* ------------------------------------------------------------------ */
/* Elementary functions                                                */
/* ------------------------------------------------------------------ */

/* exp(a + b i) = e^a (cos b + i sin b) */
COLD void
_dfloat_complex_exp(_dfloat_ops_t ops, double * res, const double * x)
{
    part_t e, s, c;
    if (part_is_zero(ops, IM(ops, x)))
    {
        ops->exp(res, x);
        ops->zero(IM(ops, res));
        return;
    }
    ops->exp(e, x);
    ops->sin_cos(s, c, IM(ops, x));
    ops->mul(res, e, c, AUTO);
    ops->mul(IM(ops, res), e, s, AUTO);
}

/* where the formula fails for a (wide) ball away from zero: log|x| in
   [log min |x|, log max |x|], and the argument g */
static COLD int
log_crude(_dfloat_ops_t ops, double * res, const double * x, const double * g)
{
    part_t a, b;
    double l = cplx_lower_abs(ops, x), u = cplx_upper_abs(ops, x);
    if (!ops->ball || !(l > 0.0))
        return GR_UNABLE;
    if (u <= DBL_MAX)
    {
        ops->set_d(a, l);
        ops->set_d(b, u);
        if (ops->log(a, a) != GR_SUCCESS || ops->log(b, b) != GR_SUCCESS)
            return GR_UNABLE;
        part_interval_signed(ops, a, part_lo(ops, a), part_hi(ops, b));
    }
    else
        part_zero_pm(ops, a, D_INF);
    part_set(ops, IM(ops, res), g);
    part_set(ops, res, a);
    return GR_SUCCESS;
}

/* log|x| + arg(x) i, with log|x| = log1p((a - 1)(a + 1) + b^2) / 2 (a the
   larger part) for |x|^2 near 1, where it is relative, and log(|x
   2^-e|^2) / 2 + e log 2 otherwise; the imaginary part is atan2, which
   for a ball crossing the branch cut encloses both sides. Balls:
   GR_DOMAIN for 0, GR_UNABLE for a ball containing 0. */
COLD int
_dfloat_complex_log(_dfloat_ops_t ops, double * res, const double * x)
{
    part_t a, b, s, t, u, g;
    int e, status;

    if (ops->ball)
    {
        if (cplx_is_zero(ops, x))
            return GR_DOMAIN;
        if (cplx_contains_zero(ops, x))
            return GR_UNABLE;
    }

    ops->atan2(g, IM(ops, x), x);

    if (part_is_W(ops, x) || part_is_W(ops, IM(ops, x)))
    {
        /* |x| is unbounded (and bounded away from zero) */
        part_zero_pm(ops, res, D_INF);
        part_set(ops, IM(ops, res), g);
        return GR_SUCCESS;
    }

    if (part_is_zero(ops, IM(ops, x)))
    {
        ops->abs(t, x);
        if (ops->log(t, t) != GR_SUCCESS)
            return log_crude(ops, res, x, g);
        part_set(ops, res, t);
        part_set(ops, IM(ops, res), g);
        return GR_SUCCESS;
    }

    e = scale_exp(x[0], IM(ops, x)[0]);
    ops->mul_2exp_si(a, x, -e);
    ops->mul_2exp_si(b, IM(ops, x), -e);
    abs2(ops, s, a, b, AUTO);
    if (e == 0 && s[0] >= 0.5 && s[0] <= 2.0)
    {
        double * p = a, * q = b;
        if (fabs(a[0]) < fabs(b[0]))
            FLINT_SWAP(double *, p, q);
        ops->one(u);
        ops->sub(s, p, u);
        ops->add(t, p, u);
        ops->mul(s, s, t, AUTO);
        ops->sqr(t, q, AUTO);
        ops->add(s, s, t);
        status = ops->log1p(t, s);
    }
    else
    {
        status = ops->log(t, s);
        if (e != 0)
        {
            part_table(ops, u, _dfloat_const[DFLOAT_CONST_LOG2]);
            ops->mul_d(u, u, (double) (2 * e));
            ops->add(t, t, u);
        }
    }
    if (status != GR_SUCCESS)
        return log_crude(ops, res, x, g);
    ops->mul_2exp_si(res, t, -1);
    part_set(ops, IM(ops, res), g);
    return GR_SUCCESS;
}

/* sin(a + b i) = sin a cosh b + i cos a sinh b,
   cos(a + b i) = cos a cosh b - i sin a sinh b
   (sn or cs may be NULL) */
COLD void
_dfloat_complex_sin_cos(_dfloat_ops_t ops, double * sn, double * cs, const double * x)
{
    part_t s, c, sh, ch;
    ops->sin_cos(s, c, x);
    if (part_is_zero(ops, IM(ops, x)))
    {
        if (sn != NULL)
        {
            part_set(ops, sn, s);
            ops->zero(IM(ops, sn));
        }
        if (cs != NULL)
        {
            part_set(ops, cs, c);
            ops->zero(IM(ops, cs));
        }
        return;
    }
    real_sinh_cosh(ops, sh, ch, IM(ops, x));
    if (sn != NULL)
    {
        ops->mul(sn, s, ch, AUTO);
        ops->mul(IM(ops, sn), c, sh, AUTO);
    }
    if (cs != NULL)
    {
        ops->mul(cs, c, ch, AUTO);
        ops->mul(IM(ops, cs), s, sh, AUTO);
        ops->neg(IM(ops, cs), IM(ops, cs));
    }
}

/* sinh(a + b i) = sinh a cos b + i cosh a sin b,
   cosh(a + b i) = cosh a cos b + i sinh a sin b
   (sh or ch may be NULL) */
COLD void
_dfloat_complex_sinh_cosh(_dfloat_ops_t ops, double * sh, double * ch, const double * x)
{
    part_t s, c, h, k;
    ops->sin_cos(s, c, IM(ops, x));
    real_sinh_cosh(ops, h, k, x);
    if (sh != NULL)
    {
        ops->mul(sh, h, c, AUTO);
        ops->mul(IM(ops, sh), k, s, AUTO);
    }
    if (ch != NULL)
    {
        ops->mul(ch, k, c, AUTO);
        ops->mul(IM(ops, ch), h, s, AUTO);
    }
}

/* for b away from zero: |tan(a + b i)|^2 = (sin^2 a + sinh^2 b) /
   (cos^2 a + sinh^2 b) <= coth^2 b, and coth |b| < 1/|b| + 1 */
static COLD int
tan_crude(_dfloat_ops_t ops, double * res, const double * b)
{
    double l = part_lower(ops, b);
    if (!(l > 0.0))
        return GR_UNABLE;
    cplx_zero_pm(ops, res, _dfloat_rad_finish(_dfloat_rad_div(1.0, l) + 1.0));
    return GR_SUCCESS;
}

/* tan: the poles are the real numbers pi/2 + k pi, which are not
   representable, so that GR_DOMAIN never occurs; GR_UNABLE for a ball
   that may contain a pole. Off the real line: (t + u i) / (1 - t u i)
   with t = tan a, u = tanh b (the addition theorem; the denominator has
   real part 1), or, where tan a fails (a ball near a pole, or the whole
   line) and b is away from zero, tan x = i (1 - w) / (1 + w) with w =
   exp(2 i x) (for b > 0, where |w| = exp(-2 b) < 1; b < 0 by
   conjugation). */
COLD int
_dfloat_complex_tan(_dfloat_ops_t ops, double * res, const double * x)
{
    cplx_t p, q;
    part_t t;
    const double * a = x, * b = IM(ops, x);
    int status, neg;

    if (part_is_zero(ops, b))
    {
        status = ops->tan(t, a);
        if (status == GR_SUCCESS)
        {
            part_set(ops, res, t);
            ops->zero(IM(ops, res));
        }
        return status;
    }

    if (part_is_zero(ops, a))
    {
        ops->tanh(t, b);
        part_set(ops, IM(ops, res), t);
        ops->zero(res);
        return GR_SUCCESS;
    }

    if (ops->tan(p, a) == GR_SUCCESS)
    {
        ops->tanh(IM(ops, p), b);
        ops->one(q);
        ops->mul(IM(ops, q), p, IM(ops, p), AUTO);
        ops->neg(IM(ops, q), IM(ops, q));
        return (_dfloat_complex_div(ops, res, p, q, AUTO) == GR_SUCCESS) ? GR_SUCCESS : GR_UNABLE;
    }

    if (part_contains_zero(ops, b))
        return GR_UNABLE;
    if (part_is_W(ops, b))
        return tan_crude(ops, res, b);

    neg = (b[0] < 0.0);
    /* w = exp(2 i x') with x' = x or conj(x), whose imaginary part |b| is
       positive: 2 i x' = -2 |b| + 2 a i */
    ops->mul_2exp_si(q, b, 1);
    if (!neg)
        ops->neg(q, q);
    ops->mul_2exp_si(IM(ops, q), a, 1);
    _dfloat_complex_exp(ops, q, q);
    ops->one(t);
    ops->sub(p, t, q);                          /* 1 - w */
    ops->neg(IM(ops, p), IM(ops, q));
    ops->add(q, q, t);                          /* 1 + w */
    if (_dfloat_complex_div(ops, p, p, q, AUTO) != GR_SUCCESS)
        return tan_crude(ops, res, b);
    /* i (p0 + p1 i) = -p1 + p0 i, conjugated for b < 0 */
    part_set(ops, t, p);
    ops->neg(res, IM(ops, p));
    if (neg)
        ops->neg(IM(ops, res), t);
    else
        part_set(ops, IM(ops, res), t);
    return GR_SUCCESS;
}

/* tanh x = -i tan(i x) */
COLD int
_dfloat_complex_tanh(_dfloat_ops_t ops, double * res, const double * x)
{
    cplx_t t;
    part_t r;
    int status;
    /* i x = -b + a i */
    ops->neg(t, IM(ops, x));
    part_set(ops, IM(ops, t), x);
    status = _dfloat_complex_tan(ops, t, t);
    if (status != GR_SUCCESS)
        return status;
    part_set(ops, r, t);
    part_set(ops, res, IM(ops, t));
    ops->neg(IM(ops, res), r);
    return GR_SUCCESS;
}

/* whether y is an integer: GR_SUCCESS (with e) for an exact integer,
   GR_DOMAIN if the ball contains no integer, else GR_UNABLE */
static COLD int
cplx_get_fmpz(_dfloat_ops_t ops, fmpz_t e, const double * y)
{
    int status;
    if (part_is_zero(ops, IM(ops, y)))
        return ops->get_fmpz(e, y);
    if (!part_contains_zero(ops, IM(ops, y)))
        return GR_DOMAIN;
    status = ops->get_fmpz(e, y);
    return (status == GR_DOMAIN) ? GR_DOMAIN : GR_UNABLE;
}

/* x^e for |e| < 2^40 by binary powering */
static COLD int
cplx_pow_si(_dfloat_ops_t ops, double * res, const double * x, slong e)
{
    cplx_t b, r;
    ulong a = FLINT_UABS(e);
    if (e < 0)
    {
        if (_dfloat_complex_inv(ops, b, x, AUTO) != GR_SUCCESS && ops->ball)
            return GR_UNABLE;
    }
    else
        cplx_set(ops, b, x);
    cplx_one(ops, r);
    while (a)
    {
        if (a & 1)
            ops->cmul(r, r, b, AUTO);
        a >>= 1;
        if (a)
            ops->csqr(b, b, AUTO);
    }
    cplx_set(ops, res, r);
    return GR_SUCCESS;
}

/* x^y = exp(y log x) (the principal branch). Plain: an integer y (|y| <
   2^40) by binary powering, 0^y = 0 for Re(y) > 0. Balls: on the
   domain x != 0, or x = 0 and Re(y) > 0, or y an integer (x != 0 if y <
   0), GR_SUCCESS inside, GR_DOMAIN outside, GR_UNABLE otherwise; an
   exact integer y by binary powering for |y| < 2^40 and by exp(y log x)
   beyond, and for x containing zero and Re(y) > 0 the enclosure |x^y|
   <= max(1, |x|)^Re(y) exp(pi |Im(y)|) of both parts. */
COLD int
_dfloat_complex_pow(_dfloat_ops_t ops, double * res, const double * x, const double * y)
{
    cplx_t l;
    fmpz_t e;
    int status, ystat;

    if (!ops->ball)
    {
        double y0 = y[0];
        if (part_is_zero(ops, IM(ops, y)) && fabs(y0) < 0x1p40 && y0 == rint(y0)
            && _dfloat_abs_sum(y + 1, ops->n - 1) == 0.0)
            return cplx_pow_si(ops, res, x, (slong) y0);
        if (cplx_is_zero(ops, x) && y0 > 0.0)
        {
            ops->zero(res);
            ops->zero(IM(ops, res));
            return GR_SUCCESS;
        }
        _dfloat_complex_log(ops, l, x);
        ops->cmul(l, l, y, AUTO);
        _dfloat_complex_exp(ops, res, l);
        return GR_SUCCESS;
    }

    fmpz_init(e);
    ystat = cplx_get_fmpz(ops, e, y);
    if (ystat == GR_SUCCESS)
    {
        if (fmpz_sgn(e) < 0 && cplx_contains_zero(ops, x))
            status = cplx_is_zero(ops, x) ? GR_DOMAIN : GR_UNABLE;
        else if (fmpz_bits(e) <= 40)
            status = cplx_pow_si(ops, res, x, fmpz_get_si(e));
        else if (!cplx_contains_zero(ops, x))
        {
            status = _dfloat_complex_log(ops, l, x);
            if (status == GR_SUCCESS)
            {
                part_t f;
                ops->set_fmpz(f, e);
                ops->mul(l, l, f, AUTO);
                ops->mul(IM(ops, l), IM(ops, l), f, AUTO);
                _dfloat_complex_exp(ops, res, l);
            }
            else
                status = GR_UNABLE;
        }
        else
        {
            /* a positive e >= 2^40 and x containing zero: |x^e| <= 1 if
               |x| <= 1, else (probably) an overflow */
            double u = _dfloat_rad_finish(part_upper(ops, x) + part_upper(ops, IM(ops, x)));
            part_zero_pm(ops, res, (u <= 1.0) ? 1.0 : D_INF);
            part_zero_pm(ops, IM(ops, res), (u <= 1.0) ? 1.0 : D_INF);
            status = GR_SUCCESS;
        }
    }
    else if (!cplx_contains_zero(ops, x))
    {
        status = _dfloat_complex_log(ops, l, x);
        if (status == GR_SUCCESS)
        {
            ops->cmul(l, l, y, AUTO);
            _dfloat_complex_exp(ops, res, l);
        }
        else
            status = GR_UNABLE;
    }
    else if (ALL_GT(y, 0.0))
    {
        if (cplx_is_zero(ops, x))
        {
            ops->zero(res);
            ops->zero(IM(ops, res));
        }
        else
        {
            /* max(1, |x|)^Re(y) exp(pi |Im(y)|) <= exp(s log(max(1, u)) +
               pi v) with u >= |x|, s >= Re(y), v >= |Im(y)| */
            part_t t, w;
            double u, s, v;
            int n = ops->n;
            u = _dfloat_rad_finish(part_upper(ops, x) + part_upper(ops, IM(ops, x)));
            s = _dfloat_rad_finish(y[0] + _dfloat_rad_finish(_dfloat_abs_sum(y + 1, n - 1) + y[n]));
            v = part_upper(ops, IM(ops, y));
            if (u <= DBL_MAX && s <= DBL_MAX && v <= DBL_MAX)
            {
                ops->set_d(t, FLINT_MAX(u, 1.0));
                GR_IGNORE(ops->log(t, t));
                ops->mul_d(t, t, s);
                part_table(ops, w, _dfloat_const[DFLOAT_CONST_PI]);
                ops->mul_d(w, w, v);
                ops->add(t, t, w);
                ops->exp(t, t);
                u = part_upper(ops, t);
            }
            else
                u = D_INF;
            part_zero_pm(ops, res, u);
            part_zero_pm(ops, IM(ops, res), u);
        }
        status = GR_SUCCESS;
    }
    else if (cplx_is_zero(ops, x) && ALL_LE(y, 0.0) && !cplx_contains_zero(ops, y))
        status = GR_DOMAIN;          /* 0^y, Re(y) <= 0, y != 0 */
    else
        status = GR_UNABLE;
    fmpz_clear(e);
    return status;
}

/* ------------------------------------------------------------------ */
/* Conversions                                                         */
/* ------------------------------------------------------------------ */

COLD void
_dfloat_complex_get_acb(_dfloat_ops_t ops, acb_t res, const double * x)
{
    ops->get_arb(acb_realref(res), x);
    ops->get_arb(acb_imagref(res), IM(ops, x));
}

COLD void
_dfloat_complex_set_acb(_dfloat_ops_t ops, double * res, const acb_t x)
{
    ops->set_arb(res, acb_realref(x));
    ops->set_arb(IM(ops, res), acb_imagref(x));
}

COLD char *
_dfloat_complex_get_str_ops(_dfloat_ops_t ops, const double * x, slong digits)
{
    int n = ops->n;
    return _dfloat_complex_get_str(x, IM(ops, x), n, ops->ball ? x[n] : -1.0,
        ops->ball ? IM(ops, x)[n] : -1.0, digits);
}

COLD void
_dfloat_complex_randtest(_dfloat_ops_t ops, double * res, flint_rand_t state, int special)
{
    if (special)
        ops->randtest_special(res, state);
    else
        ops->randtest(res, state);
    if (n_randint(state, special ? 4 : 8) == 0)
        ops->zero(IM(ops, res));
    else if (special)
        ops->randtest_special(IM(ops, res), state);
    else
        ops->randtest(IM(ops, res), state);
    if (!special && n_randint(state, 16) == 0)
        ops->zero(res);
}

/* ------------------------------------------------------------------ */
/* Generic ring methods                                                */
/* ------------------------------------------------------------------ */

#define OPS(ctx) _dfloat_ops(DFLOAT_CTX_N(ctx), DFLOAT_CTX_BALL(ctx))
#define MODE(ctx) (DFLOAT_CTX_FAST(ctx) ? DFLOAT_MODE_FAST : DFLOAT_MODE_AUTO)

/* the real ring of the parts */
#define REAL_CTX(rctx, ctx) GR_MUST_SUCCEED(gr_ctx_init_dfloat(rctx, DFLOAT_CTX_N(ctx), DFLOAT_CTX_BALL(ctx) ? DFLOAT_BALL : 0))

/* canonical parts in the strong contexts */
static COLD int
finish(double * res, gr_ctx_t ctx, int status)
{
    if (DFLOAT_CTX_STRONG(ctx))
    {
        _dfloat_ops_t ops = OPS(ctx);
        ops->canonicalise(res);
        ops->canonicalise(IM(ops, res));
    }
    return status;
}

static COLD int
_gr_dfloat_complex_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    int status = GR_SUCCESS, n = DFLOAT_CTX_N(ctx);
    if (DFLOAT_CTX_BALL(ctx))
        status |= gr_stream_write(out, "Complex numbers as rectangular balls with ");
    else
        status |= gr_stream_write(out, "Complex floating-point numbers with ");
    status |= gr_stream_write_si(out, n);
    status |= gr_stream_write(out, DFLOAT_CTX_BALL(ctx) ? "-term double expansion midpoints (d"
                                                        : "-term double expansion parts (d");
    status |= gr_stream_write_si(out, n);
    status |= gr_stream_write(out, DFLOAT_CTX_BALL(ctx) ? "cb" : "c");
    status |= gr_stream_write(out, DFLOAT_CTX_STRONG(ctx) ? ", strong)" : ")");
    return status;
}

static COLD int
_gr_dfloat_complex_randtest(double * res, flint_rand_t state, gr_ctx_t ctx)
{
    _dfloat_complex_randtest(OPS(ctx), res, state, n_randint(state, 16) == 0);
    return GR_SUCCESS;
}

static COLD int
_gr_dfloat_complex_write(gr_stream_t out, const double * x, gr_ctx_t ctx)
{
    _dfloat_ops_t ops = OPS(ctx);
    if (part_is_zero(ops, IM(ops, x)))
    {
        /* a real number, as the real ring writes it */
        gr_ctx_t rctx;
        int status;
        REAL_CTX(rctx, ctx);
        status = gr_write(out, x, rctx);
        gr_ctx_clear(rctx);
        return status;
    }
    return gr_stream_write_free(out, _dfloat_complex_get_str_ops(ops, x, 0));
}

/* midpoints (the whole line as such) */
static COLD int
_gr_dfloat_complex_write_n(gr_stream_t out, const double * x, slong digits, gr_ctx_t ctx)
{
    _dfloat_ops_t ops = OPS(ctx);
    int n = ops->n;
    return gr_stream_write_free(out, _dfloat_complex_get_str(x, IM(ops, x), n,
        part_is_W(ops, x) ? x[n] : -1.0, part_is_W(ops, IM(ops, x)) ? IM(ops, x)[n] : -1.0,
        FLINT_MAX(digits, 1)));
}

static COLD int
_gr_dfloat_complex_pi(double * res, gr_ctx_t ctx)
{
    _dfloat_ops_t ops = OPS(ctx);
    part_table(ops, res, _dfloat_const[DFLOAT_CONST_PI]);
    ops->zero(IM(ops, res));
    return GR_SUCCESS;
}

static COLD int
_gr_dfloat_complex_pos_inf(double * res, gr_ctx_t ctx)
{
    _dfloat_ops_t ops = OPS(ctx);
    ops->set_d(res, D_INF);
    ops->zero(IM(ops, res));
    return GR_SUCCESS;
}

static COLD int
_gr_dfloat_complex_neg_inf(double * res, gr_ctx_t ctx)
{
    _dfloat_ops_t ops = OPS(ctx);
    ops->set_d(res, -D_INF);
    ops->zero(IM(ops, res));
    return GR_SUCCESS;
}

static COLD int
_gr_dfloat_complex_nan(double * res, gr_ctx_t ctx)
{
    _dfloat_ops_t ops = OPS(ctx);
    ops->set_d(res, D_NAN);
    ops->set_d(IM(ops, res), D_NAN);
    return GR_SUCCESS;
}

/* the conversions from real types, through the real ring of the parts */
static COLD int
_gr_dfloat_complex_set_si(double * res, slong x, gr_ctx_t ctx)
{
    gr_ctx_t rctx;
    int status;
    REAL_CTX(rctx, ctx);
    status = gr_set_si(res, x, rctx);
    status |= gr_zero(IM(OPS(ctx), res), rctx);
    gr_ctx_clear(rctx);
    return status;
}

static COLD int
_gr_dfloat_complex_set_ui(double * res, ulong x, gr_ctx_t ctx)
{
    gr_ctx_t rctx;
    int status;
    REAL_CTX(rctx, ctx);
    status = gr_set_ui(res, x, rctx);
    status |= gr_zero(IM(OPS(ctx), res), rctx);
    gr_ctx_clear(rctx);
    return status;
}

static COLD int
_gr_dfloat_complex_set_fmpz(double * res, const fmpz_t x, gr_ctx_t ctx)
{
    gr_ctx_t rctx;
    int status;
    REAL_CTX(rctx, ctx);
    status = gr_set_fmpz(res, x, rctx);
    status |= gr_zero(IM(OPS(ctx), res), rctx);
    gr_ctx_clear(rctx);
    return status;
}

static COLD int
_gr_dfloat_complex_set_fmpq(double * res, const fmpq_t x, gr_ctx_t ctx)
{
    gr_ctx_t rctx;
    int status;
    REAL_CTX(rctx, ctx);
    status = gr_set_fmpq(res, x, rctx);
    status |= gr_zero(IM(OPS(ctx), res), rctx);
    gr_ctx_clear(rctx);
    return status;
}

static COLD int
_gr_dfloat_complex_set_d(double * res, double x, gr_ctx_t ctx)
{
    gr_ctx_t rctx;
    int status;
    REAL_CTX(rctx, ctx);
    status = gr_set_d(res, x, rctx);
    status |= gr_zero(IM(OPS(ctx), res), rctx);
    gr_ctx_clear(rctx);
    return status;
}

static COLD int
_gr_dfloat_complex_set_other(double * res, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
{
    _dfloat_ops_t ops = OPS(ctx);
    gr_ctx_t rctx, xctx;
    int status;

    switch (x_ctx->which_ring)
    {
        case GR_CTX_FMPZ:
        case GR_CTX_FMPQ:
        case GR_CTX_REAL_FLOAT_ARF:
        case GR_CTX_RR_ARB:
        case GR_CTX_DFLOAT:
        case GR_CTX_DFLOAT_BALL:
            REAL_CTX(rctx, ctx);
            status = gr_set_other(res, x, x_ctx, rctx);
            status |= gr_zero(IM(ops, res), rctx);
            gr_ctx_clear(rctx);
            return status;
        case GR_CTX_CC_ACB:
            /* (nonfinite parts: the whole line, or nan) */
            _dfloat_complex_set_acb(ops, res, x);
            return GR_SUCCESS;
        case GR_CTX_DFLOAT_COMPLEX:
        case GR_CTX_DFLOAT_COMPLEX_BALL:
        {
            /* part by part, through the real rings */
            int nx = DFLOAT_CTX_N(x_ctx), bx = DFLOAT_CTX_BALL(x_ctx);
            REAL_CTX(rctx, ctx);
            GR_MUST_SUCCEED(gr_ctx_init_dfloat(xctx, nx, bx ? DFLOAT_BALL : 0));
            status = gr_set_other(res, x, xctx, rctx);
            status |= gr_set_other(IM(ops, res), (const double *) x + nx + bx, xctx, rctx);
            gr_ctx_clear(xctx);
            gr_ctx_clear(rctx);
            return status;
        }
        default:
        {
            acb_t t;
            gr_ctx_init_complex_acb(rctx, 53 * ops->n + 30);
            acb_init(t);
            status = gr_set_other(t, x, x_ctx, rctx);
            if (status == GR_SUCCESS)
                _dfloat_complex_set_acb(ops, res, t);
            acb_clear(t);
            gr_ctx_clear(rctx);
            return status;
        }
    }
}

/* the real-valued conversions: a zero imaginary part is required
   (GR_DOMAIN for one away from zero, GR_UNABLE for a ball containing
   zero, or a nan) */
static COLD int
real_check(_dfloat_ops_t ops, const double * x)
{
    const double * im = IM(ops, x);
    if (part_is_zero(ops, im))
        return GR_SUCCESS;
    if (!ops->ball)
        return (im[0] != im[0]) ? GR_UNABLE : GR_DOMAIN;
    return part_contains_zero(ops, im) ? GR_UNABLE : GR_DOMAIN;
}

static COLD int
_gr_dfloat_complex_get_fmpz(fmpz_t res, const double * x, gr_ctx_t ctx)
{
    _dfloat_ops_t ops = OPS(ctx);
    gr_ctx_t rctx;
    int status;
    if (ops->ball)
        return cplx_get_fmpz(ops, res, x);
    status = real_check(ops, x);
    if (status != GR_SUCCESS)
        return status;
    REAL_CTX(rctx, ctx);
    status = gr_get_fmpz(res, x, rctx);
    gr_ctx_clear(rctx);
    return status;
}

static COLD int
_gr_dfloat_complex_get_si(slong * res, const double * x, gr_ctx_t ctx)
{
    gr_ctx_t rctx;
    int status = real_check(OPS(ctx), x);
    if (status != GR_SUCCESS)
        return status;
    REAL_CTX(rctx, ctx);
    status = gr_get_si(res, x, rctx);
    gr_ctx_clear(rctx);
    return status;
}

static COLD int
_gr_dfloat_complex_get_ui(ulong * res, const double * x, gr_ctx_t ctx)
{
    gr_ctx_t rctx;
    int status = real_check(OPS(ctx), x);
    if (status != GR_SUCCESS)
        return status;
    REAL_CTX(rctx, ctx);
    status = gr_get_ui(res, x, rctx);
    gr_ctx_clear(rctx);
    return status;
}

static COLD int
_gr_dfloat_complex_get_d(double * res, const double * x, gr_ctx_t ctx)
{
    gr_ctx_t rctx;
    int status = real_check(OPS(ctx), x);
    if (status != GR_SUCCESS)
        return status;
    REAL_CTX(rctx, ctx);
    status = gr_get_d(res, x, rctx);
    gr_ctx_clear(rctx);
    return status;
}

static COLD truth_t
_gr_dfloat_complex_is_invertible(const double * x, gr_ctx_t ctx)
{
    _dfloat_ops_t ops = OPS(ctx);
    if (cplx_is_zero(ops, x))
        return T_FALSE;
    if (cplx_contains_zero(ops, x))
        return T_UNKNOWN;
    return T_TRUE;
}

/* the plain ring: GR_DOMAIN for an exact singularity; the ball ring:
   the statuses of the functions */
static COLD int
_gr_dfloat_complex_div(double * res, const double * x, const double * y, gr_ctx_t ctx)
{
    _dfloat_ops_t ops = OPS(ctx);
    if (!ops->ball && cplx_is_zero(ops, y))
        return GR_DOMAIN;
    return finish(res, ctx, _dfloat_complex_div(ops, res, x, y, MODE(ctx)));
}

static COLD int
_gr_dfloat_complex_inv(double * res, const double * x, gr_ctx_t ctx)
{
    _dfloat_ops_t ops = OPS(ctx);
    if (!ops->ball && cplx_is_zero(ops, x))
        return GR_DOMAIN;
    return finish(res, ctx, _dfloat_complex_inv(ops, res, x, MODE(ctx)));
}

static COLD int
_gr_dfloat_complex_sqrt(double * res, const double * x, gr_ctx_t ctx)
{
    _dfloat_complex_sqrt(OPS(ctx), res, x, MODE(ctx));
    return finish(res, ctx, GR_SUCCESS);
}

static COLD int
_gr_dfloat_complex_rsqrt(double * res, const double * x, gr_ctx_t ctx)
{
    _dfloat_ops_t ops = OPS(ctx);
    if (!ops->ball && cplx_is_zero(ops, x))
        return GR_DOMAIN;
    return finish(res, ctx, _dfloat_complex_rsqrt(ops, res, x, MODE(ctx)));
}

static COLD int
_gr_dfloat_complex_abs(double * res, const double * x, gr_ctx_t ctx)
{
    _dfloat_ops_t ops = OPS(ctx);
    part_t t;
    _dfloat_complex_abs(ops, t, x);
    part_set(ops, res, t);
    ops->zero(IM(ops, res));
    return finish(res, ctx, GR_SUCCESS);
}

static COLD int
_gr_dfloat_complex_arg(double * res, const double * x, gr_ctx_t ctx)
{
    _dfloat_ops_t ops = OPS(ctx);
    part_t t;
    _dfloat_complex_arg(ops, t, x);
    part_set(ops, res, t);
    ops->zero(IM(ops, res));
    return finish(res, ctx, GR_SUCCESS);
}

/* x / |x| (0 for 0; a ball containing zero gives [0 +/- 1] parts) */
static COLD int
_gr_dfloat_complex_sgn(double * res, const double * x, gr_ctx_t ctx)
{
    _dfloat_ops_t ops = OPS(ctx);
    part_t a, t;
    if (cplx_is_zero(ops, x))
    {
        ops->zero(res);
        ops->zero(IM(ops, res));
        return GR_SUCCESS;
    }
    if (!cplx_contains_zero(ops, x))
    {
        _dfloat_complex_abs(ops, a, x);
        ops->one(t);
        if (ops->div(t, t, a, AUTO) == GR_SUCCESS)
        {
            ops->mul(res, x, t, AUTO);
            ops->mul(IM(ops, res), IM(ops, x), t, AUTO);
            return finish(res, ctx, GR_SUCCESS);
        }
    }
    part_zero_pm(ops, res, 1.0);
    part_zero_pm(ops, IM(ops, res), 1.0);
    return GR_SUCCESS;
}

static COLD int
_gr_dfloat_complex_exp(double * res, const double * x, gr_ctx_t ctx)
{
    _dfloat_complex_exp(OPS(ctx), res, x);
    return finish(res, ctx, GR_SUCCESS);
}

/* (plain: log 0 = -inf, as in the real ring) */
static COLD int
_gr_dfloat_complex_log(double * res, const double * x, gr_ctx_t ctx)
{
    return finish(res, ctx, _dfloat_complex_log(OPS(ctx), res, x));
}

static COLD int
_gr_dfloat_complex_sin(double * res, const double * x, gr_ctx_t ctx)
{
    _dfloat_complex_sin_cos(OPS(ctx), res, NULL, x);
    return finish(res, ctx, GR_SUCCESS);
}

static COLD int
_gr_dfloat_complex_cos(double * res, const double * x, gr_ctx_t ctx)
{
    _dfloat_complex_sin_cos(OPS(ctx), NULL, res, x);
    return finish(res, ctx, GR_SUCCESS);
}

static COLD int
_gr_dfloat_complex_sin_cos(double * sn, double * cs, const double * x, gr_ctx_t ctx)
{
    _dfloat_complex_sin_cos(OPS(ctx), sn, cs, x);
    finish(sn, ctx, GR_SUCCESS);
    return finish(cs, ctx, GR_SUCCESS);
}

static COLD int
_gr_dfloat_complex_sinh(double * res, const double * x, gr_ctx_t ctx)
{
    _dfloat_complex_sinh_cosh(OPS(ctx), res, NULL, x);
    return finish(res, ctx, GR_SUCCESS);
}

static COLD int
_gr_dfloat_complex_cosh(double * res, const double * x, gr_ctx_t ctx)
{
    _dfloat_complex_sinh_cosh(OPS(ctx), NULL, res, x);
    return finish(res, ctx, GR_SUCCESS);
}

static COLD int
_gr_dfloat_complex_tan(double * res, const double * x, gr_ctx_t ctx)
{
    return finish(res, ctx, _dfloat_complex_tan(OPS(ctx), res, x));
}

static COLD int
_gr_dfloat_complex_tanh(double * res, const double * x, gr_ctx_t ctx)
{
    return finish(res, ctx, _dfloat_complex_tanh(OPS(ctx), res, x));
}

/* plain: 0^y is 0 for Re(y) > 0, 1 for y = 0 and undefined otherwise */
static COLD int
_gr_dfloat_complex_pow(double * res, const double * x, const double * y, gr_ctx_t ctx)
{
    _dfloat_ops_t ops = OPS(ctx);
    if (!ops->ball && cplx_is_zero(ops, x) && !(y[0] > 0.0) && !cplx_is_zero(ops, y))
        return (y[0] != y[0] || IM(ops, y)[0] != IM(ops, y)[0]) ? GR_UNABLE : GR_DOMAIN;
    return finish(res, ctx, _dfloat_complex_pow(ops, res, x, y));
}

/* the methods of both complex rings that are implemented here; the
   formats' own tables (the arithmetic and the vectors) extend these */
gr_method_tab_input _dfloat_complex_gr_methods_input[] =
{
    {GR_METHOD_CTX_WRITE,       (gr_funcptr) _gr_dfloat_complex_ctx_write},
    {GR_METHOD_RANDTEST,        (gr_funcptr) _gr_dfloat_complex_randtest},
    {GR_METHOD_WRITE,           (gr_funcptr) _gr_dfloat_complex_write},
    {GR_METHOD_WRITE_N,         (gr_funcptr) _gr_dfloat_complex_write_n},
    {GR_METHOD_PI,              (gr_funcptr) _gr_dfloat_complex_pi},
    {GR_METHOD_SET_SI,          (gr_funcptr) _gr_dfloat_complex_set_si},
    {GR_METHOD_SET_UI,          (gr_funcptr) _gr_dfloat_complex_set_ui},
    {GR_METHOD_SET_FMPZ,        (gr_funcptr) _gr_dfloat_complex_set_fmpz},
    {GR_METHOD_SET_FMPQ,        (gr_funcptr) _gr_dfloat_complex_set_fmpq},
    {GR_METHOD_SET_D,           (gr_funcptr) _gr_dfloat_complex_set_d},
    {GR_METHOD_SET_STR,         (gr_funcptr) gr_generic_set_str_ring_exponents},
    {GR_METHOD_SET_OTHER,       (gr_funcptr) _gr_dfloat_complex_set_other},
    {GR_METHOD_GET_SI,          (gr_funcptr) _gr_dfloat_complex_get_si},
    {GR_METHOD_GET_UI,          (gr_funcptr) _gr_dfloat_complex_get_ui},
    {GR_METHOD_GET_FMPZ,        (gr_funcptr) _gr_dfloat_complex_get_fmpz},
    {GR_METHOD_GET_D,           (gr_funcptr) _gr_dfloat_complex_get_d},
    {GR_METHOD_DIV,             (gr_funcptr) _gr_dfloat_complex_div},
    {GR_METHOD_INV,             (gr_funcptr) _gr_dfloat_complex_inv},
    {GR_METHOD_SQRT,            (gr_funcptr) _gr_dfloat_complex_sqrt},
    {GR_METHOD_RSQRT,           (gr_funcptr) _gr_dfloat_complex_rsqrt},
    {GR_METHOD_ABS,             (gr_funcptr) _gr_dfloat_complex_abs},
    {GR_METHOD_ARG,             (gr_funcptr) _gr_dfloat_complex_arg},
    {GR_METHOD_SGN,             (gr_funcptr) _gr_dfloat_complex_sgn},
    {GR_METHOD_EXP,             (gr_funcptr) _gr_dfloat_complex_exp},
    {GR_METHOD_LOG,             (gr_funcptr) _gr_dfloat_complex_log},
    {GR_METHOD_SIN,             (gr_funcptr) _gr_dfloat_complex_sin},
    {GR_METHOD_COS,             (gr_funcptr) _gr_dfloat_complex_cos},
    {GR_METHOD_SIN_COS,         (gr_funcptr) _gr_dfloat_complex_sin_cos},
    {GR_METHOD_TAN,             (gr_funcptr) _gr_dfloat_complex_tan},
    {GR_METHOD_SINH,            (gr_funcptr) _gr_dfloat_complex_sinh},
    {GR_METHOD_COSH,            (gr_funcptr) _gr_dfloat_complex_cosh},
    {GR_METHOD_TANH,            (gr_funcptr) _gr_dfloat_complex_tanh},
    {GR_METHOD_POW,             (gr_funcptr) _gr_dfloat_complex_pow},
    {0,                         (gr_funcptr) NULL},
};

/* the plain ring only */
gr_method_tab_input _dfloat_complex_plain_gr_methods_input[] =
{
    {GR_METHOD_POS_INF,         (gr_funcptr) _gr_dfloat_complex_pos_inf},
    {GR_METHOD_NEG_INF,         (gr_funcptr) _gr_dfloat_complex_neg_inf},
    {GR_METHOD_UINF,            (gr_funcptr) gr_not_in_domain},
    {GR_METHOD_UNDEFINED,       (gr_funcptr) _gr_dfloat_complex_nan},
    {GR_METHOD_UNKNOWN,         (gr_funcptr) _gr_dfloat_complex_nan},
    {0,                         (gr_funcptr) NULL},
};

/* the ball ring only */
gr_method_tab_input _dfloat_complex_ball_gr_methods_input[] =
{
    {GR_METHOD_IS_INVERTIBLE,   (gr_funcptr) _gr_dfloat_complex_is_invertible},
    {0,                         (gr_funcptr) NULL},
};
