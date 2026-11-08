/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "acb.h"
#include "acf.h"
#include "impl.h"

/*
    The less common complex elementary functions of acf (the inverse
    functions, expm1, log1p, sinc): rounded from the balls of the acb
    functions, evaluated at increasing working precisions.

    In the strict modes each part is rounded when both endpoints of its
    ball round to the same number; in the relaxed modes when the radii lie
    below 2^-(prec + 4) (ARF_RND_FAST) resp. 2^-(prec + 8) (ARF_RND_ACCURATE)
    times the larger part, ARF_RND_FAST truncating the endpoint farther
    from zero (within an ulp of the larger part below the value and a
    sliver above it) and ARF_RND_ACCURATE rounding the midpoint to nearest.

    The parts of these functions at dyadic arguments off the axes are
    transcendental, so that the strict loop terminates; its cost grows
    with the logarithm of the ratio of the parts of the argument when one
    is much smaller than the other (the size of the small part of the
    value, which acb computes with an error relative to the larger one).
    On the axes a part may be exactly zero or a multiple of pi/2 (of 1/2
    for the pi versions), which no ball could decide: the real arguments
    in the real domain of a function go to the real function of arf, and
    for the other arguments on the axes those parts are set directly (with
    the branch conventions of acb, the reciprocal functions being f(1/z)
    as in gr), the others coming from the loop.

    In the strict modes an argument z = z0 + eta very close to a point z0
    of an axis (|eta| < 2^-(prec + 64) |z0|) where f is analytic in a disk
    around z0 is evaluated from f(z0 + eta) = f(z0) + f'(z0) eta + R,
    |R| <= M t^2 / (1 - t) with t = |eta| / r and M a bound for |f| on the
    disk of radius r (the Cauchy estimates; acb bounds f on a box around
    z0), the parts of f(z0) and f'(z0) that vanish on the axis set to
    exact zeros: the small part of the value then comes out with its own
    relative accuracy at the usual working precision, rather than at a
    precision growing with log(|z0| / |eta|).
*/

enum
{
    F_EXPM1, F_LOG1P, F_SINC, F_SINC_PI,
    F_ASIN, F_ACOS, F_ATAN, F_ACOT, F_ASEC, F_ACSC,
    F_ASINH, F_ACOSH, F_ATANH, F_ACOTH, F_ASECH, F_ACSCH,
    F_ASIN_PI, F_ACOS_PI, F_ATAN_PI
};

/* the real domains */
enum { D_ALL, D_GT_M1, D_ABS_LE1, D_NONZERO, D_ABS_GE1, D_GE1, D_ABS_LT1, D_ABS_GT1, D_ASECH };

typedef int (*_real_func_t)(arf_ptr, arf_srcptr, slong, arf_rnd_t);

static const struct
{
    _real_func_t real;
    unsigned char domain;
}
_acb_funcs[] =
{
    { arf_expm1, D_ALL }, { arf_log1p, D_GT_M1 }, { arf_sinc, D_ALL }, { arf_sinc_pi, D_ALL },
    { arf_asin, D_ABS_LE1 }, { arf_acos, D_ABS_LE1 }, { arf_atan, D_ALL }, { arf_acot, D_NONZERO },
    { arf_asec, D_ABS_GE1 }, { arf_acsc, D_ABS_GE1 },
    { arf_asinh, D_ALL }, { arf_acosh, D_GE1 }, { arf_atanh, D_ABS_LT1 }, { arf_acoth, D_ABS_GT1 },
    { arf_asech, D_ASECH }, { arf_acsch, D_NONZERO },
    { arf_asin_pi, D_ABS_LE1 }, { arf_acos_pi, D_ABS_LE1 }, { arf_atan_pi, D_ALL },
};

static int
_in_domain(const arf_t x, int d)
{
    int c;

    switch (d)
    {
        case D_ALL: return 1;
        case D_GT_M1: return arf_cmp_si(x, -1) > 0;
        case D_NONZERO: return !arf_is_zero(x);
        case D_GE1: return arf_cmp_si(x, 1) >= 0;
        case D_ASECH: return arf_sgn(x) > 0 && arf_cmp_si(x, 1) <= 0;
        default:
            c = arf_cmpabs_ui(x, 1);
            return (d == D_ABS_LE1) ? (c <= 0) : (d == D_ABS_GE1) ? (c >= 0) : (d == D_ABS_LT1) ? (c < 0) : (c > 0);
    }
}

/* the value at the working precision wp */
static void
_acb_eval(acb_t r, int f, const acb_t z, slong wp)
{
    switch (f)
    {
        case F_EXPM1: acb_expm1(r, z, wp); break;
        case F_LOG1P: acb_log1p(r, z, wp); break;
        case F_SINC: acb_sinc(r, z, wp); break;
        case F_SINC_PI: acb_sinc_pi(r, z, wp); break;
        case F_ASIN: case F_ASIN_PI: acb_asin(r, z, wp); break;
        case F_ACOS: case F_ACOS_PI: acb_acos(r, z, wp); break;
        case F_ATAN: case F_ATAN_PI: acb_atan(r, z, wp); break;
        case F_ASINH: acb_asinh(r, z, wp); break;
        case F_ACOSH: acb_acosh(r, z, wp); break;
        case F_ATANH: acb_atanh(r, z, wp); break;
        default:
            /* the reciprocal functions: f(1/z) */
            acb_inv(r, z, wp);
            switch (f)
            {
                case F_ACOT: acb_atan(r, r, wp); break;
                case F_ASEC: acb_acos(r, r, wp); break;
                case F_ACSC: acb_asin(r, r, wp); break;
                case F_ACOTH: acb_atanh(r, r, wp); break;
                case F_ASECH: acb_acosh(r, r, wp); break;
                default: acb_asinh(r, r, wp); break;
            }
    }

    if (f >= F_ASIN_PI)
    {
        arb_t t;
        arb_init(t);
        arb_const_pi(t, wp);
        acb_div_arb(r, r, t, wp);
        arb_clear(t);
    }
}

/* the kinds of the parts on the axes */
enum { P_LOOP, P_ZERO, P_CONST, P_NAN };

/* the parts of f on the axes (x real outside the real domain of f, or
   imaginary): *kre, *kim the kinds, *cre, *cim the multiples of pi/2 (of
   1/2 for the pi versions) of the constant parts; 0 when neither part is
   exact (the loop computes both) */
static int
_acf_axis_parts(int * kre, int * cre, int * kim, int * cim, int f, const arf_t a, const arf_t b)
{
    int base = (f >= F_ASIN_PI) ? f - F_ASIN_PI + F_ASIN : f;
    int s, c;

    *kre = *kim = P_LOOP;
    *cre = *cim = 0;

    if (arf_is_zero(b))
    {
        /* a real, outside the real domain */
        s = arf_sgn(a);
        c = arf_cmpabs_ui(a, 1);

        switch (base)
        {
            case F_LOG1P:
                /* log(-1 - a) + pi i for a < -1 */
                if (c == 0)
                    *kre = P_NAN;
                else
                {
                    if (arf_cmp_si(a, -2) == 0)
                        *kre = P_ZERO;
                    *kim = P_CONST; *cim = 2;
                }
                return 1;
            case F_ASIN:
                *kre = P_CONST; *cre = s;
                return 1;
            case F_ACOS:
                if (s > 0) *kre = P_ZERO; else { *kre = P_CONST; *cre = 2; }
                return 1;
            case F_ACOT:
                *kre = P_NAN;
                return 1;
            case F_ASEC:
                if (s == 0) *kre = P_NAN;
                else if (s > 0) *kre = P_ZERO;
                else { *kre = P_CONST; *cre = 2; }
                return 1;
            case F_ACSC:
                if (s == 0) *kre = P_NAN;
                else { *kre = P_CONST; *cre = s; }
                return 1;
            case F_ACOSH:
                if (c < 0 || (s < 0 && c == 0)) *kre = P_ZERO;
                if (s < 0 && c >= 0) { *kim = P_CONST; *cim = 2; }
                return 1;
            case F_ATANH:
                if (c == 0) *kre = P_NAN;
                else { *kim = P_CONST; *cim = -s; }
                return 1;
            case F_ACOTH:
                if (s == 0 || c == 0) *kre = P_NAN;
                else { *kim = P_CONST; *cim = -s; }
                return 1;
            case F_ASECH:
                if (s == 0) *kre = P_NAN;
                else if (s > 0 || c > 0) *kre = P_ZERO;
                else
                {
                    if (c == 0) *kre = P_ZERO;
                    *kim = P_CONST; *cim = 2;
                }
                return 1;
            case F_ACSCH:
                *kre = P_NAN;
                return 1;
            default:
                return 0;
        }
    }
    else
    {
        /* a = 0, b nonzero */
        s = arf_sgn(b);
        c = arf_cmpabs_ui(b, 1);

        switch (base)
        {
            case F_SINC:
            case F_SINC_PI:
                *kim = P_ZERO;
                return 1;
            case F_ASIN:
            case F_ACSC:
            case F_ATANH:
            case F_ACOTH:
                *kre = P_ZERO;
                return 1;
            case F_ACOS:
            case F_ASEC:
                *kre = P_CONST; *cre = 1;
                return 1;
            case F_ATAN:
                if (c == 0) *kre = P_NAN;
                else if (c < 0) *kre = P_ZERO;
                else { *kre = P_CONST; *cre = s; }
                return 1;
            case F_ACOT:
                if (c == 0) *kre = P_NAN;
                else if (c > 0) *kre = P_ZERO;
                else { *kre = P_CONST; *cre = -s; }
                return 1;
            case F_ASINH:
                if (c <= 0) *kre = P_ZERO;
                else { *kim = P_CONST; *cim = s; }
                return 1;
            case F_ACOSH:
                *kim = P_CONST; *cim = s;
                return 1;
            case F_ASECH:
                *kim = P_CONST; *cim = -s;
                return 1;
            case F_ACSCH:
                if (c >= 0) *kre = P_ZERO;
                else { *kim = P_CONST; *cim = -s; }
                return 1;
            default:
                return 0;
        }
    }
}

/* res = c pi / 2 (c / 2 for the pi versions), c = +-1, +-2 */
static int
_acf_set_const(arf_t res, int c, int pi, slong prec, arf_rnd_t rnd)
{
    int s;

    if (pi)
    {
        /* c / 2: exact */
        arf_set_si_2exp_si(res, c, -1);
        return 0;
    }

    s = arf_const_pi(res, prec, (c < 0) ? _acf_rnd_flip(rnd) : rnd);
    if (FLINT_ABS(c) == 1)
        arf_mul_2exp_si(res, res, -1);
    if (c < 0)
        arf_neg(res, res);
    return s;
}

/* the rounding of the ball v in a strict mode: the inexact flag, or
   ARF_ELEM_UNCERTAIN when its endpoints (rounded outward to wp bits:
   rounding is monotone) round differently */
static int
_acf_round_ball_strict(arf_t res, const arb_t v, slong wp, slong prec, arf_rnd_t rnd)
{
    arf_t lo, hi;
    int ok;

    if (mag_is_zero(arb_radref(v)))
        return arf_set_round(res, arb_midref(v), prec, rnd);

    arf_init(lo);
    arf_init(hi);
    arb_get_lbound_arf(lo, v, wp);
    arb_get_ubound_arf(hi, v, wp);
    arf_set_round(lo, lo, prec, rnd);
    arf_set_round(hi, hi, prec, rnd);
    ok = arf_equal(lo, hi);
    if (ok)
        arf_swap(res, lo);
    arf_clear(lo);
    arf_clear(hi);
    return ok ? 1 : ARF_ELEM_UNCERTAIN;
}

/* the rounding of the ball v in a relaxed mode (its radius small enough
   against the larger part) */
static void
_acf_round_ball_relaxed(arf_t res, const arb_t v, slong prec, arf_rnd_t rnd)
{
    if (rnd == ARF_RND_FAST)
    {
        /* the endpoint farther from zero, truncated */
        if (arf_sgn(arb_midref(v)) >= 0)
            arb_get_ubound_arf(res, v, prec + 8);
        else
            arb_get_lbound_arf(res, v, prec + 8);
        arf_set_round(res, res, prec, ARF_RND_DOWN);
    }
    else
        arf_set_round(res, arb_midref(v), prec, ARF_RND_NEAR);
}

/* the spread of the exponents of the parts of x: |E(a) - E(b)| (both
   nonzero) plus the larger of |E(a)|, |E(b)|, at most ACF_ELEM_SPREAD_MAX */
static slong
_acf_exp_spread(const acf_t x)
{
    arf_srcptr a = acf_realref(x), b = acf_imagref(x);
    fmpz_t t, u;
    slong r;

    fmpz_init(t);
    fmpz_init(u);
    if (!arf_is_special(a) && !arf_is_special(b))
    {
        fmpz_sub(t, ARF_EXPREF(a), ARF_EXPREF(b));
        fmpz_abs(t, t);
    }
    if (!arf_is_special(a))
    {
        fmpz_abs(u, ARF_EXPREF(a));
        fmpz_add(t, t, u);
    }
    if (!arf_is_special(b))
    {
        fmpz_abs(u, ARF_EXPREF(b));
        fmpz_add(t, t, u);
    }
    r = (fmpz_cmp_si(t, ACF_ELEM_SPREAD_MAX) > 0) ? ACF_ELEM_SPREAD_MAX : fmpz_get_si(t);
    fmpz_clear(t);
    fmpz_clear(u);
    return r;
}

/* the parts of f(x) wanted (want bit 0 the real part, bit 1 the
   imaginary part) into res, the others left alone; the inexact flags */
static int
_acf_acb_loop(acf_t res, int f, const acf_t x, int want, slong prec, arf_rnd_t rnd)
{
    acb_t z, r;
    arf_t t1, t2;
    mag_t m, u;
    slong wp, acc, need, cap;
    int relaxed = (rnd >= ARF_RND_FAST), s1, s2, status;

    acb_init(z);
    acb_init(r);
    arf_init(t1);
    arf_init(t2);
    mag_init(m);
    mag_init(u);

    arb_set_arf(acb_realref(z), acf_realref(x));
    arb_set_arf(acb_imagref(z), acf_imagref(x));

    need = prec + (relaxed ? ((rnd == ARF_RND_FAST) ? 4 : 8) : 0);
    /* (the cap grows with the spread of the exponents of the parts of x,
       which bounds that of the result's parts where the functions are
       analytic and away from their cuts and singular points; up to
       ACF_ELEM_SPREAD_MAX bits) */
    cap = 4 * prec + ACF_ELEM_RETRY_LIMIT_BITS + 2 * _acf_exp_spread(x);
    wp = prec + 32;

    for (;;)
    {
        _acb_eval(r, f, z, wp);

        if (acb_is_finite(r))
        {
            if (relaxed)
            {
                /* the radii against the larger part (its lower bound) */
                arb_get_mag_lower(m, acb_realref(r));
                arb_get_mag_lower(u, acb_imagref(r));
                mag_max(m, m, u);
                mag_zero(u);
                if (want & 1) mag_max(u, u, arb_radref(acb_realref(r)));
                if (want & 2) mag_max(u, u, arb_radref(acb_imagref(r)));
                mag_mul_2exp_si(u, u, need);

                if (!mag_is_zero(m) && mag_cmp(u, m) <= 0)
                {
                    if (want & 1) _acf_round_ball_relaxed(t1, acb_realref(r), prec, rnd);
                    if (want & 2) _acf_round_ball_relaxed(t2, acb_imagref(r), prec, rnd);
                    status = want;
                    break;
                }

                acc = acb_rel_accuracy_bits(r);
            }
            else
            {
                s1 = (want & 1) ? _acf_round_ball_strict(t1, acb_realref(r), wp + 8, prec, rnd) : 0;
                s2 = (want & 2) ? _acf_round_ball_strict(t2, acb_imagref(r), wp + 8, prec, rnd) : 0;

                if (s1 != ARF_ELEM_UNCERTAIN && s2 != ARF_ELEM_UNCERTAIN)
                {
                    status = s1 | (s2 << 1);
                    break;
                }

                acc = WORD_MAX;
                if (s1 == ARF_ELEM_UNCERTAIN) acc = FLINT_MIN(acc, arb_rel_accuracy_bits(acb_realref(r)));
                if (s2 == ARF_ELEM_UNCERTAIN) acc = FLINT_MIN(acc, arb_rel_accuracy_bits(acb_imagref(r)));
            }
        }
        else
            acc = 0;

        if (wp > cap)
        {
            /* beyond the limit of the retries: the midpoints as they
               stand (nan when the balls are not finite) */
            if (!acb_is_finite(r))
            {
                arf_nan(t1);
                arf_nan(t2);
                status = 0;
            }
            else
            {
                arf_rnd_t r1 = relaxed ? arf_rnd_relaxed_to_strict(rnd) : rnd;
                status = 0;
                if (want & 1) status |= arf_set_round(t1, arb_midref(acb_realref(r)), prec, r1);
                if (want & 2) status |= arf_set_round(t2, arb_midref(acb_imagref(r)), prec, r1) << 1;
            }
            break;
        }

        /* the missing accuracy, or (near a rounding boundary) half as much
           again */
        /* (at least geometric growth: a part that vanishes to more bits
           than the accuracy shows is found in a logarithmic number of
           steps) */
        if (acc < need + 8)
            wp += FLINT_MAX(need + 16 - FLINT_MAX(acc, 0), wp / 2);
        else
            wp += wp / 2;
        wp = FLINT_MIN(wp, cap + 1);
    }

    if (want & 1) arf_swap(acf_realref(res), t1);
    if (want & 2) arf_swap(acf_imagref(res), t2);

    acb_clear(z);
    acb_clear(r);
    arf_clear(t1);
    arf_clear(t2);
    mag_clear(m);
    mag_clear(u);
    return status;
}

/* f'(z) (off the branch cuts) */
static void
_acb_deriv(acb_t r, int f, const acb_t z, slong wp)
{
    int base = (f >= F_ASIN_PI) ? f - F_ASIN_PI + F_ASIN : f;
    acb_t t, u;

    acb_init(t);
    acb_init(u);

    switch (base)
    {
        case F_EXPM1:
            acb_exp(r, z, wp);
            break;
        case F_LOG1P:
            acb_add_ui(t, z, 1, wp);
            acb_inv(r, t, wp);
            break;
        case F_SINC:
        case F_SINC_PI:
            /* (cos z - sinc z) / z, resp. (cos pi z - sinc_pi z) / z */
            if (base == F_SINC) { acb_cos(t, z, wp); acb_sinc(u, z, wp); }
            else { acb_cos_pi(t, z, wp); acb_sinc_pi(u, z, wp); }
            acb_sub(t, t, u, wp);
            acb_div(r, t, z, wp);
            break;
        case F_ASIN:
        case F_ACOS:
        case F_ASINH:
            /* -+1 / sqrt(1 -+ z^2) */
            acb_sqr(t, z, wp);
            if (base == F_ASINH) acb_add_ui(t, t, 1, wp);
            else { acb_neg(t, t); acb_add_ui(t, t, 1, wp); }
            acb_rsqrt(r, t, wp);
            if (base == F_ACOS) acb_neg(r, r);
            break;
        case F_ATAN:
        case F_ACOT:
        case F_ATANH:
        case F_ACOTH:
            /* -+1 / (1 +- z^2) */
            acb_sqr(t, z, wp);
            if (base == F_ATAN || base == F_ACOT) acb_add_ui(t, t, 1, wp);
            else { acb_neg(t, t); acb_add_ui(t, t, 1, wp); }
            acb_inv(r, t, wp);
            if (base == F_ACOT) acb_neg(r, r);
            break;
        case F_ACOSH:
            /* 1 / (sqrt(z - 1) sqrt(z + 1)) */
            acb_sub_ui(t, z, 1, wp);
            acb_add_ui(u, z, 1, wp);
            acb_sqrt(t, t, wp);
            acb_sqrt(u, u, wp);
            acb_mul(t, t, u, wp);
            acb_inv(r, t, wp);
            break;
        default:
            /* the reciprocal functions: g(1/z)' = -g'(1/z) / z^2 */
            {
                int g = (base == F_ASEC) ? F_ACOS : (base == F_ACSC) ? F_ASIN : (base == F_ASECH) ? F_ACOSH : F_ASINH;
                acb_inv(u, z, wp);
                _acb_deriv(t, g, u, wp);
                acb_sqr(u, z, wp);
                acb_div(r, t, u, wp);
                acb_neg(r, r);
            }
            break;
    }

    if (f >= F_ASIN_PI)
    {
        arb_t c;
        arb_init(c);
        arb_const_pi(c, wp);
        acb_div_arb(r, r, c, wp);
        arb_clear(c);
    }

    acb_clear(t);
    acb_clear(u);
}

/* the branch cuts (the singularities lie on them): segments of the real
   (axis 0) and imaginary (axis 1) axes */
static const struct { signed char n, axis[2]; double lo[2], hi[2]; } _acb_cuts[] =
{
    { 0, {0, 0}, {0, 0}, {0, 0} },                          /* expm1 */
    { 1, {0, 0}, {-HUGE_VAL, 0}, {-1, 0} },                  /* log1p */
    { 0, {0, 0}, {0, 0}, {0, 0} },                          /* sinc */
    { 0, {0, 0}, {0, 0}, {0, 0} },                          /* sinc_pi */
    { 2, {0, 0}, {-HUGE_VAL, 1}, {-1, HUGE_VAL} },           /* asin */
    { 2, {0, 0}, {-HUGE_VAL, 1}, {-1, HUGE_VAL} },           /* acos */
    { 2, {1, 1}, {-HUGE_VAL, 1}, {-1, HUGE_VAL} },           /* atan */
    { 1, {1, 0}, {-1, 0}, {1, 0} },                          /* acot */
    { 1, {0, 0}, {-1, 0}, {1, 0} },                          /* asec */
    { 1, {0, 0}, {-1, 0}, {1, 0} },                          /* acsc */
    { 2, {1, 1}, {-HUGE_VAL, 1}, {-1, HUGE_VAL} },           /* asinh */
    { 1, {0, 0}, {-HUGE_VAL, 0}, {1, 0} },                   /* acosh */
    { 2, {0, 0}, {-HUGE_VAL, 1}, {-1, HUGE_VAL} },           /* atanh */
    { 1, {0, 0}, {-1, 0}, {1, 0} },                          /* acoth */
    { 2, {0, 0}, {-HUGE_VAL, 1}, {0, HUGE_VAL} },            /* asech */
    { 1, {1, 0}, {-1, 0}, {1, 0} },                          /* acsch */
};

/* the distance from the point t (axis 0: t real, axis 1: t i) to the
   cuts of f (HUGE_VAL for none) */
static double
_acb_cut_dist(int f, int axis, double t)
{
    int base = (f >= F_ASIN_PI) ? f - F_ASIN_PI + F_ASIN : f;
    double d = HUGE_VAL, e, w;
    int i;

    for (i = 0; i < _acb_cuts[base].n; i++)
    {
        double lo = _acb_cuts[base].lo[i], hi = _acb_cuts[base].hi[i];

        if (_acb_cuts[base].axis[i] == axis)
            e = (t < lo) ? lo - t : (t > hi) ? t - hi : 0.0;
        else
        {
            /* the other axis: through the distance of 0 to the segment */
            w = (0 < lo) ? lo : (0 > hi) ? -hi : 0.0;
            e = sqrt(t * t + w * w);
        }
        d = FLINT_MIN(d, e);
    }

    return d;
}

/* f(z0 + eta) for z very close to the point z0 of an axis where f is
   analytic (the strict modes): the inexact flags, or ARF_ELEM_UNCERTAIN
   when the expansion does not apply or decide the rounding */
static int
_acf_acb_near_axis(acf_t res, int f, const acf_t x, slong prec, arf_rnd_t rnd)
{
    arf_srcptr a = acf_realref(x), b = acf_imagref(x);
    int axis, kre, cre, kim, cim, status = ARF_ELEM_UNCERTAIN, s1, s2;
    slong wp;
    double t0, dist, rad;
    acb_t z0, eta, c0, c1, box, v;
    mag_t M, t, err;
    arf_t r1, r2;

    /* z0 = a (axis 0) or b i (axis 1), eta the small part */
    axis = (arf_cmpabs(a, b) >= 0) ? 0 : 1;
    if (arf_is_zero(axis ? a : b))
        return ARF_ELEM_UNCERTAIN;

    /* (z0 in double, clamped to [2^-500, 2^500] in absolute value: the
       cuts lie at distances 0 and 1 from the origin, so that the clamped
       distance errs on the safe side by far less than the margin of
       rad = dist / 4) */
    if (arf_cmpabs_2exp_si(axis ? b : a, 500) > 0)
        t0 = arf_sgn(axis ? b : a) * ldexp(1.0, 500);
    else if (arf_cmpabs_2exp_si(axis ? b : a, -500) < 0)
        t0 = 0.0;
    else
        t0 = arf_get_d(axis ? b : a, ARF_RND_NEAR);
    dist = _acb_cut_dist(f, axis, t0);
    if (!(dist > 0))
        return ARF_ELEM_UNCERTAIN;

    /* the parts of f on the axis segment around z0 (constant kinds there:
       the segment avoids the cuts, where the kinds change) */
    if (axis == 0)
    {
        /* off the cuts on the real axis: the real domain, f real */
        kre = P_LOOP;
        kim = P_ZERO;
    }
    else
    {
        arf_t zero;
        arf_init(zero);
        if (!_acf_axis_parts(&kre, &cre, &kim, &cim, f, zero, b))
            kre = kim = P_LOOP;
        arf_clear(zero);
        if (kre == P_NAN || kim == P_NAN)
            return ARF_ELEM_UNCERTAIN;
    }

    rad = FLINT_MIN(dist / 4, 1.0);

    acb_init(z0);
    acb_init(eta);
    acb_init(c0);
    acb_init(c1);
    acb_init(box);
    acb_init(v);
    mag_init(M);
    mag_init(t);
    mag_init(err);
    arf_init(r1);
    arf_init(r2);

    if (axis == 0)
    {
        arb_set_arf(acb_realref(z0), a);
        arb_set_arf(acb_realref(eta), b);
        acb_mul_onei(eta, eta);
    }
    else
    {
        arb_set_arf(acb_imagref(z0), b);
        arb_set_arf(acb_realref(eta), a);
    }

    /* M: |f| on the box z0 + [-rad, rad] + [-rad, rad] i, containing the
       disk of radius rad */
    acb_set(box, z0);
    mag_set_d(arb_radref(acb_realref(box)), rad);
    mag_set_d(arb_radref(acb_imagref(box)), rad);
    _acb_eval(v, f, box, 64);
    acb_get_mag(M, v);

    /* t = |eta| / r' with r' = rad (1 - 2^-20) <= rad */
    acb_get_mag(t, eta);
    mag_set_d_lower(err, rad * (1 - ldexp(1, -20)));
    mag_div(t, t, err);

    if (mag_is_finite(M) && mag_cmp_2exp_si(t, -2) < 0)
    {
        /* err = M t^2 / (1 - t) <= M t^2 (1 + 2t) */
        mag_mul(err, t, t);
        mag_mul(err, err, M);
        mag_mul_2exp_si(t, t, 1);
        mag_add_ui(t, t, 1);
        mag_mul(err, err, t);

        for (wp = prec + 32; wp <= 8 * prec + 256; wp *= 2)
        {
            _acb_eval(c0, f, z0, wp);
            _acb_deriv(c1, f, z0, wp);

            /* the parts vanishing on the segment: f(z0) by the kinds; f'
               real on the real axis, and on the imaginary axis real when
               Re f is constant along it and imaginary when Im f is */
            if (kre == P_ZERO) arb_zero(acb_realref(c0));
            if (kim == P_ZERO) arb_zero(acb_imagref(c0));
            if (axis == 0 || kre != P_LOOP) arb_zero(acb_imagref(c1));
            if (axis == 1 && kim != P_LOOP) arb_zero(acb_realref(c1));

            acb_mul(v, c1, eta, wp);
            acb_add(v, v, c0, wp);
            arb_add_error_mag(acb_realref(v), err);
            arb_add_error_mag(acb_imagref(v), err);

            if (!acb_is_finite(v))
                break;

            s1 = _acf_round_ball_strict(r1, acb_realref(v), wp + 8, prec, rnd);
            s2 = (s1 == ARF_ELEM_UNCERTAIN) ? s1 : _acf_round_ball_strict(r2, acb_imagref(v), wp + 8, prec, rnd);

            if (s1 != ARF_ELEM_UNCERTAIN && s2 != ARF_ELEM_UNCERTAIN)
            {
                arf_swap(acf_realref(res), r1);
                arf_swap(acf_imagref(res), r2);
                status = s1 | (s2 << 1);
                break;
            }
        }
    }

    acb_clear(z0);
    acb_clear(eta);
    acb_clear(c0);
    acb_clear(c1);
    acb_clear(box);
    acb_clear(v);
    mag_clear(M);
    mag_clear(t);
    mag_clear(err);
    arf_clear(r1);
    arf_clear(r2);
    return status;
}

static int
_acf_acb_func(acf_t res, int f, const acf_t x, slong prec, arf_rnd_t rnd)
{
    arf_srcptr a = acf_realref(x), b = acf_imagref(x);
    int kre, cre, kim, cim, want, status;
    acf_t t;

    if ((arf_is_special(a) && !arf_is_zero(a)) || (arf_is_special(b) && !arf_is_zero(b)))
    {
        arf_nan(acf_realref(res));
        arf_nan(acf_imagref(res));
        return 0;
    }

    /* the real functions in their real domains */
    if (arf_is_zero(b) && _in_domain(a, _acb_funcs[f].domain))
    {
        status = _acb_funcs[f].real(acf_realref(res), a, prec, rnd);
        arf_zero(acf_imagref(res));
        return status;
    }

    /* (strict modes: a part much smaller than the other) */
    if (rnd < ARF_RND_FAST && !arf_is_zero(a) && !arf_is_zero(b))
    {
        slong d = _acf_part_gap(a, b) + _acf_part_gap(b, a);

        if (d > prec + 64)
        {
            acf_init(t);
            status = _acf_acb_near_axis(t, f, x, prec, rnd);
            if (status != ARF_ELEM_UNCERTAIN)
                acf_swap(res, t);
            acf_clear(t);
            if (status != ARF_ELEM_UNCERTAIN)
                return status;
        }
    }

    if (!((arf_is_zero(a) || arf_is_zero(b)) && _acf_axis_parts(&kre, &cre, &kim, &cim, f, a, b)))
        kre = kim = P_LOOP;

    if (kre == P_NAN || kim == P_NAN)
    {
        arf_nan(acf_realref(res));
        arf_nan(acf_imagref(res));
        return 0;
    }

    want = (kre == P_LOOP) | ((kim == P_LOOP) << 1);

    /* (through a temporary: res may alias x) */
    acf_init(t);
    status = (want != 0) ? _acf_acb_loop(t, f, x, want, prec, rnd) : 0;

    if (kre == P_ZERO)
        arf_zero(acf_realref(t));
    else if (kre == P_CONST)
        status |= _acf_set_const(acf_realref(t), cre, f >= F_ASIN_PI, prec, rnd);

    if (kim == P_ZERO)
        arf_zero(acf_imagref(t));
    else if (kim == P_CONST)
        status |= _acf_set_const(acf_imagref(t), cim, f >= F_ASIN_PI, prec, rnd) << 1;

    acf_swap(res, t);
    acf_clear(t);
    return status;
}

int
acf_expm1(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
{
    arf_srcptr a = acf_realref(x), b = acf_imagref(x);

    /* Re z <= -(prec + 3): e^z - 1 = -(1 + delta) + i e^a sin b with
       delta = -e^a cos b, |delta| < 2^(-prec - 1), whose sign decides the
       real part directly (the strict modes could not resolve it from
       balls); the imaginary part is that of e^z */
    if (!arf_is_special(a) && !arf_is_special(b) && arf_cmp_si(a, -prec - 3) <= 0)
    {
        acf_t t;
        arf_t c;
        int s, dsgn;

        acf_init(t);
        arf_init(c);
        s = acf_exp(t, x, prec, rnd) & 2;
        arf_cos(c, b, 16, ARF_RND_DOWN);
        dsgn = (arf_sgn(c) > 0);
        s |= _arf_one_plus_tiny(acf_realref(t), 1, dsgn, prec, rnd);
        acf_swap(res, t);
        acf_clear(t);
        arf_clear(c);
        return s;
    }

    return _acf_acb_func(res, F_EXPM1, x, prec, rnd);
}
int acf_log1p(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd) { return _acf_acb_func(res, F_LOG1P, x, prec, rnd); }
int acf_sinc(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd) { return _acf_acb_func(res, F_SINC, x, prec, rnd); }
int acf_sinc_pi(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd) { return _acf_acb_func(res, F_SINC_PI, x, prec, rnd); }
int acf_asin(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd) { return _acf_acb_func(res, F_ASIN, x, prec, rnd); }
int acf_acos(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd) { return _acf_acb_func(res, F_ACOS, x, prec, rnd); }
int acf_atan(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd) { return _acf_acb_func(res, F_ATAN, x, prec, rnd); }
int acf_acot(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd) { return _acf_acb_func(res, F_ACOT, x, prec, rnd); }
int acf_asec(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd) { return _acf_acb_func(res, F_ASEC, x, prec, rnd); }
int acf_acsc(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd) { return _acf_acb_func(res, F_ACSC, x, prec, rnd); }
int acf_asinh(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd) { return _acf_acb_func(res, F_ASINH, x, prec, rnd); }
int acf_acosh(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd) { return _acf_acb_func(res, F_ACOSH, x, prec, rnd); }
int acf_atanh(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd) { return _acf_acb_func(res, F_ATANH, x, prec, rnd); }
int acf_acoth(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd) { return _acf_acb_func(res, F_ACOTH, x, prec, rnd); }
int acf_asech(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd) { return _acf_acb_func(res, F_ASECH, x, prec, rnd); }
int acf_acsch(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd) { return _acf_acb_func(res, F_ACSCH, x, prec, rnd); }
int acf_asin_pi(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd) { return _acf_acb_func(res, F_ASIN_PI, x, prec, rnd); }
int acf_acos_pi(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd) { return _acf_acb_func(res, F_ACOS_PI, x, prec, rnd); }
int acf_atan_pi(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd) { return _acf_acb_func(res, F_ATAN_PI, x, prec, rnd); }
