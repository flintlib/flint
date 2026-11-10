/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Reference values for the tests of mp_real, computed with MPFR.  The
   elementary functions of arb are built on the kernels of mp_real, so
   they cannot serve as independent references.

   ref_<f>(res, x, prec) sets the ball res to f(x) for the ball x, with
   the signature of arb_<f>: f at the midpoint (taken exactly),
   correctly rounded by MPFR at prec + REF_GUARD bits (a radius of one
   ulp there, or zero when MPFR reports the result exact), plus, for an
   inexact x, the propagated error r C with C a bound for |f'| over the
   ball, from its endpoints (where |f'| is largest for these functions
   on the balls accepted: the tangents on a ball that keeps away from a
   pole).  Outside the domain of f the result is indeterminate.  The
   constants, logarithms of integers and the arctangents of fractions
   follow. */

#ifndef MP_REAL_TEST_REF_H
#define MP_REAL_TEST_REF_H

#include <mpfr.h>
#include "arb.h"

/* (arf.h declares its MPFR conversions only when mpfr.h was included first) */
int arf_get_mpfr(mpfr_t x, const arf_t y, mpfr_rnd_t rnd);
void arf_set_mpfr(arf_t x, const mpfr_t y);

#define REF_GUARD 30

/* the exponent range of MPFR opened as far as it goes (once) */
static void
ref_init(void)
{
    static int done = 0;
    if (!done)
    {
        mpfr_set_emin(mpfr_get_emin_min());
        mpfr_set_emax(mpfr_get_emax_max());
        done = 1;
    }
}

static int ref_cospi_fn(mpfr_ptr r, mpfr_srcptr x, mpfr_rnd_t rnd);

/* res = v within one ulp of v's precision (exact when tern = 0), or the
   special value */
static void
ref_set_mpfr(arb_t res, mpfr_srcptr v, int tern)
{
    if (mpfr_nan_p(v))
        arb_indeterminate(res);
    else if (mpfr_inf_p(v))
    {
        if (mpfr_sgn(v) > 0)
            arb_pos_inf(res);
        else
            arb_neg_inf(res);
    }
    else if (mpfr_zero_p(v) && tern != 0)
    {
        /* underflow (MPFR's exponent range is that of a long, 32 bits
           on some 64-bit systems): a nonzero value below 2^emin */
        arf_zero(arb_midref(res));
        mag_set_ui_2exp_si(arb_radref(res), 1, mpfr_get_emin());
    }
    else
    {
        arf_set_mpfr(arb_midref(res), v);
        mag_zero(arb_radref(res));
        if (tern != 0)
            arf_mag_set_ulp(arb_radref(res), arb_midref(res), mpfr_get_prec(v));
    }
}

/* v = m exactly (the precision of v adjusted to the mantissa) */
static void
ref_mpfr_set_arf(mpfr_t v, const arf_t m)
{
    slong b = FLINT_MAX(arf_bits(m), MPFR_PREC_MIN);
    if (mpfr_get_prec(v) < b)
        mpfr_set_prec(v, b);
    arf_get_mpfr(v, m, MPFR_RNDN);
}

/* the endpoints lo <= hi of the ball x (rounded outward at wp bits) */
static void
ref_endpoints(mpfr_t lo, mpfr_t hi, const arb_t x, slong wp)
{
    arf_t a, b;
    arf_init(a);
    arf_init(b);
    arb_get_interval_arf(a, b, x, wp);
    mpfr_set_prec(lo, wp);
    mpfr_set_prec(hi, wp);
    arf_get_mpfr(lo, a, MPFR_RNDD);
    arf_get_mpfr(hi, b, MPFR_RNDU);
    arf_clear(a);
    arf_clear(b);
}

typedef int (* ref_mpfr_fn)(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t);

typedef struct
{
    ref_mpfr_fn f;       /* the function */
    ref_mpfr_fn df;      /* an upper bound for |f'|, or NULL for the constant 1 */
    double dmul;         /* a constant factor on the bound (pi for the functions of pi x) */
    int domain;          /* 0 all reals; 1 x > 0; 2 |x| < 1; 3 x >= 1; 4 x > -1;
                            5 cos x != 0 on the ball; 6 cos(pi x) != 0;
                            7 |x| <= 1 (an infinite radius from the
                            derivative at the endpoints) */
}
ref_spec;

/* C = an upper bound for |f'| over [lo, hi]: df at the endpoints
   (rounded up), the larger one */
static void
ref_deriv_bound(mag_t C, const ref_spec * S, mpfr_srcptr lo, mpfr_srcptr hi)
{
    mpfr_t t, u;
    arf_t a;

    if (S->df == NULL)
    {
        mag_one(C);
    }
    else
    {
        mpfr_init2(t, 64);
        mpfr_init2(u, 64);
        S->df(t, lo, MPFR_RNDU);
        S->df(u, hi, MPFR_RNDU);
        mpfr_abs(t, t, MPFR_RNDU);
        mpfr_abs(u, u, MPFR_RNDU);
        if (mpfr_cmp(u, t) > 0)
            mpfr_set(t, u, MPFR_RNDU);

        if (!mpfr_number_p(t))
            mag_inf(C);
        else
        {
            arf_init(a);
            arf_set_mpfr(a, t);
            arf_get_mag(C, a);
            arf_clear(a);
        }
        mpfr_clear(t);
        mpfr_clear(u);
    }

    if (S->dmul != 1.0)
    {
        mag_t m;
        mag_init(m);
        mag_set_d(m, S->dmul * (1 + 1e-12));
        mag_mul(C, C, m);
        mag_clear(m);
    }
}

/* res = f(x) by the specification S */
static void
ref_eval_spec(arb_t res, const arb_t x, slong prec, const ref_spec * S)
{
    mpfr_t v, in, lo, hi;
    slong wp = prec + REF_GUARD;
    int ok = 1, tern;

    ref_init();

    if (!arb_is_finite(x))
    {
        arb_indeterminate(res);
        return;
    }

    mpfr_init2(v, wp);
    mpfr_init2(in, wp);
    mpfr_init2(lo, wp);
    mpfr_init2(hi, wp);

    ref_endpoints(lo, hi, x, wp);

    switch (S->domain)
    {
        case 1: ok = mpfr_sgn(lo) > 0; break;
        case 2: ok = mpfr_cmp_si(lo, -1) > 0 && mpfr_cmp_ui(hi, 1) < 0; break;
        case 3: ok = mpfr_cmp_ui(lo, 1) >= 0; break;
        case 4: ok = mpfr_cmp_si(lo, -1) > 0; break;
        case 7: ok = mpfr_cmp_si(lo, -1) >= 0 && mpfr_cmp_ui(hi, 1) <= 0; break;
        case 5:
        case 6:
            /* the zeros of the cosine are pi (resp. 1) apart: a ball
               of width below 2 (resp. 1/2) with the same sign of the
               cosine at both ends (an exact zero is a pole) */
            {
                mpfr_t s;
                int sl, sh;
                mpfr_init2(s, 64);
                if (S->domain == 5)
                {
                    mpfr_cos(s, lo, MPFR_RNDN); sl = mpfr_sgn(s);
                    mpfr_cos(s, hi, MPFR_RNDN); sh = mpfr_sgn(s);
                    ok = (sl * sh > 0) && mag_cmp_2exp_si(arb_radref(x), 0) < 0;
                }
                else
                {
                    ref_cospi_fn(s, lo, MPFR_RNDN); sl = mpfr_sgn(s);
                    ref_cospi_fn(s, hi, MPFR_RNDN); sh = mpfr_sgn(s);
                    ok = (sl * sh > 0) && mag_cmp_2exp_si(arb_radref(x), -2) < 0;
                }
                mpfr_clear(s);
            }
            break;
        default: break;
    }

    if (!ok)
    {
        arb_indeterminate(res);
    }
    else
    {
        mag_t r;
        mag_init(r);
        mag_set(r, arb_radref(x));      /* (res may alias x) */

        ref_mpfr_set_arf(in, arb_midref(x));
        tern = S->f(v, in, MPFR_RNDN);
        ref_set_mpfr(res, v, tern);

        if (!mag_is_zero(r) && arb_is_finite(res))
        {
            mag_t C;
            mag_init(C);
            ref_deriv_bound(C, S, lo, hi);
            mag_mul(C, C, r);
            mag_add(arb_radref(res), arb_radref(res), C);
            mag_clear(C);
        }
        mag_clear(r);
    }

    mpfr_clear(v);
    mpfr_clear(in);
    mpfr_clear(lo);
    mpfr_clear(hi);
}

/* bounds for |f'| as functions of one argument (rounded up when called
   with MPFR_RNDU) */
static int ref_d_recip(mpfr_ptr r, mpfr_srcptr x, mpfr_rnd_t rnd) { return mpfr_ui_div(r, 1, x, rnd); }
static int ref_d_recip1p(mpfr_ptr r, mpfr_srcptr x, mpfr_rnd_t rnd) { mpfr_add_ui(r, x, 1, MPFR_RNDD); return mpfr_ui_div(r, 1, r, rnd); }
static int ref_d_sec2(mpfr_ptr r, mpfr_srcptr x, mpfr_rnd_t rnd) { mpfr_cos(r, x, MPFR_RNDZ); mpfr_sqr(r, r, MPFR_RNDZ); return mpfr_ui_div(r, 1, r, rnd); }
static int ref_d_rsqrt_1m(mpfr_ptr r, mpfr_srcptr x, mpfr_rnd_t rnd) { mpfr_sqr(r, x, MPFR_RNDU); mpfr_ui_sub(r, 1, r, MPFR_RNDZ); return mpfr_rec_sqrt(r, r, rnd); }
static int ref_d_rsqrt_m1(mpfr_ptr r, mpfr_srcptr x, mpfr_rnd_t rnd) { mpfr_sqr(r, x, MPFR_RNDZ); mpfr_sub_ui(r, r, 1, MPFR_RNDZ); return mpfr_rec_sqrt(r, r, rnd); }
static int ref_d_recip_1m2(mpfr_ptr r, mpfr_srcptr x, mpfr_rnd_t rnd) { mpfr_sqr(r, x, MPFR_RNDU); mpfr_ui_sub(r, 1, r, MPFR_RNDZ); return mpfr_ui_div(r, 1, r, rnd); }

/* sin, cos, tan (which = 0, 1, 2) of pi x, correctly rounded, with the
   functions of MPFR 4.1 (mpfr_sinu and the like came with 4.2): x reduced
   exactly to y in [0, 1/4] by the period, the symmetries and the swap of
   sine and cosine about 1/4, so that the zeros and units (and tan = +-1
   at y = 1/4) come out exactly; then pi y at increasing working
   precisions until mpfr_can_round decides.  For y in (0, 1/4] the sine
   and cosine of t = pi y keep the relative error of t (t cot t <= 1,
   t tan t < 1): with pi and the product rounded to nearest, the values
   are within 2^(3 - wp) relative, the quotient within 2^(4 - wp). */
static int
ref_trigpi(mpfr_ptr r, mpfr_srcptr x, mpfr_rnd_t rnd, int which)
{
    mpfr_t y, t, p, s, c, v;
    mpfr_prec_t wp, py = mpfr_get_prec(x) + 4;
    int sgn_s = 1, sgn_c = 1, swap = 0, tern;

    if (!mpfr_number_p(x))
    {
        mpfr_set_nan(r);
        return 0;
    }

    /* y = |x| mod 2 in [0, 2), exactly (sin and tan being odd, cos even) */
    if (mpfr_sgn(x) < 0)
        sgn_s = -1;
    mpfr_init2(y, py);
    mpfr_abs(y, x, MPFR_RNDN);
    mpfr_div_2ui(y, y, 1, MPFR_RNDN);
    mpfr_frac(y, y, MPFR_RNDN);
    mpfr_mul_2ui(y, y, 1, MPFR_RNDN);

    if (mpfr_cmp_ui(y, 1) >= 0)
    {
        /* sin(pi (y - 1) + pi) = -sin, cos likewise */
        mpfr_sub_ui(y, y, 1, MPFR_RNDN);
        sgn_s = -sgn_s;
        sgn_c = -sgn_c;
    }
    if (mpfr_cmp_ui_2exp(y, 1, -1) > 0)
    {
        /* sin(pi y) = sin(pi (1 - y)), cos(pi y) = -cos(pi (1 - y)) */
        mpfr_ui_sub(y, 1, y, MPFR_RNDN);
        sgn_c = -sgn_c;
    }
    if (mpfr_cmp_ui_2exp(y, 1, -2) > 0)
    {
        /* sin(pi y) = cos(pi (1/2 - y)) and the other way around:
           y = (1 - 2y) / 2 */
        mpfr_mul_2ui(y, y, 1, MPFR_RNDN);
        mpfr_ui_sub(y, 1, y, MPFR_RNDN);
        mpfr_div_2ui(y, y, 1, MPFR_RNDN);
        swap = 1;
    }

    if (mpfr_zero_p(y))
    {
        /* sin = 0, cos = 1 before the swap: exact */
        int sv = swap ? sgn_s : 0, cv = swap ? 0 : sgn_c;
        if (which == 0)
            tern = mpfr_set_si(r, sv, rnd);
        else if (which == 1)
            tern = mpfr_set_si(r, cv, rnd);
        else if (cv == 0)
        {
            mpfr_set_inf(r, sv);
            tern = 0;
        }
        else
            tern = mpfr_set_si(r, sv * cv, rnd);
        mpfr_clear(y);
        return tern;
    }

    if (which == 2 && mpfr_cmp_ui_2exp(y, 1, -2) == 0)
    {
        /* tan(pi/4) = 1 */
        tern = mpfr_set_si(r, sgn_s * sgn_c, rnd);
        mpfr_clear(y);
        return tern;
    }

    for (wp = mpfr_get_prec(r) + 32; ; wp += wp / 2)
    {
        mpfr_inits2(wp, t, p, s, c, v, (mpfr_ptr) 0);
        mpfr_const_pi(p, MPFR_RNDN);
        mpfr_mul(t, p, y, MPFR_RNDN);
        mpfr_sin_cos(s, c, t, MPFR_RNDN);
        if (swap)
            mpfr_swap(s, c);
        if (sgn_s < 0)
            mpfr_neg(s, s, MPFR_RNDN);
        if (sgn_c < 0)
            mpfr_neg(c, c, MPFR_RNDN);
        if (which == 0)
            mpfr_set(v, s, MPFR_RNDN);
        else if (which == 1)
            mpfr_set(v, c, MPFR_RNDN);
        else
            mpfr_div(v, s, c, MPFR_RNDN);

        if (mpfr_cmpabs_ui(v, 1) == 0)
        {
            /* +-1 within the error (cos t near t = 0, or tan t near
               t = pi/4): the value lies strictly on a known side, below 1
               in absolute value except for the swapped quotient cot t > 1,
               and within 2^(6 - wp) < 2^-(prec + 2) of it, so that the
               neighbour of v at wp bits on that side rounds as it does */
            int below = !(which == 2 && swap);
            if ((mpfr_sgn(v) > 0) == below)
                mpfr_nextbelow(v);
            else
                mpfr_nextabove(v);
            tern = mpfr_set(r, v, rnd);
            mpfr_clears(t, p, s, c, v, (mpfr_ptr) 0);
            break;
        }

        if (mpfr_can_round(v, wp - 6, MPFR_RNDN, rnd, mpfr_get_prec(r) + (rnd == MPFR_RNDN)))
        {
            tern = mpfr_set(r, v, rnd);
            mpfr_clears(t, p, s, c, v, (mpfr_ptr) 0);
            break;
        }
        mpfr_clears(t, p, s, c, v, (mpfr_ptr) 0);
    }

    mpfr_clear(y);
    return tern;
}

/* the functions of pi x */
static int ref_sinpi_fn(mpfr_ptr r, mpfr_srcptr x, mpfr_rnd_t rnd) { return ref_trigpi(r, x, rnd, 0); }
static int ref_cospi_fn(mpfr_ptr r, mpfr_srcptr x, mpfr_rnd_t rnd) { return ref_trigpi(r, x, rnd, 1); }
static int ref_tanpi_fn(mpfr_ptr r, mpfr_srcptr x, mpfr_rnd_t rnd) { return ref_trigpi(r, x, rnd, 2); }
static int ref_d_secpi2(mpfr_ptr r, mpfr_srcptr x, mpfr_rnd_t rnd) { ref_cospi_fn(r, x, MPFR_RNDZ); mpfr_sqr(r, r, MPFR_RNDZ); return mpfr_ui_div(r, 1, r, rnd); }

#define REF_PI 3.141592653589793

/* the specifications: f, a bound for |f'| (NULL: 1), a factor, the domain */
static const ref_spec ref_spec_exp = {mpfr_exp, mpfr_exp, 1.0, 0};
static const ref_spec ref_spec_expm1 = {mpfr_expm1, mpfr_exp, 1.0, 0};
static const ref_spec ref_spec_log = {mpfr_log, ref_d_recip, 1.0, 1};
static const ref_spec ref_spec_log1p = {mpfr_log1p, ref_d_recip1p, 1.0, 4};
static const ref_spec ref_spec_sin = {mpfr_sin, NULL, 1.0, 0};
static const ref_spec ref_spec_cos = {mpfr_cos, NULL, 1.0, 0};
static const ref_spec ref_spec_tan = {mpfr_tan, ref_d_sec2, 1.0, 5};
static const ref_spec ref_spec_atan = {mpfr_atan, NULL, 1.0, 0};
static const ref_spec ref_spec_sinh = {mpfr_sinh, mpfr_cosh, 1.0, 0};
static const ref_spec ref_spec_cosh = {mpfr_cosh, mpfr_cosh, 1.0, 0};
static const ref_spec ref_spec_tanh = {mpfr_tanh, NULL, 1.0, 0};
static const ref_spec ref_spec_asin = {mpfr_asin, ref_d_rsqrt_1m, 1.0, 7};
static const ref_spec ref_spec_acos = {mpfr_acos, ref_d_rsqrt_1m, 1.0, 7};
static const ref_spec ref_spec_asinh = {mpfr_asinh, NULL, 1.0, 0};
static const ref_spec ref_spec_acosh = {mpfr_acosh, ref_d_rsqrt_m1, 1.0, 3};
static const ref_spec ref_spec_atanh = {mpfr_atanh, ref_d_recip_1m2, 1.0, 2};
static const ref_spec ref_spec_sin_pi = {ref_sinpi_fn, NULL, REF_PI, 0};
static const ref_spec ref_spec_cos_pi = {ref_cospi_fn, NULL, REF_PI, 0};
static const ref_spec ref_spec_tan_pi = {ref_tanpi_fn, ref_d_secpi2, REF_PI, 6};

#define REF_FN(name) \
    static void ref_##name(arb_t res, const arb_t x, slong prec) \
    { ref_eval_spec(res, x, prec, &ref_spec_##name); }

REF_FN(log) REF_FN(log1p) REF_FN(sin) REF_FN(cos) REF_FN(tan) REF_FN(atan)
REF_FN(tanh) REF_FN(asin) REF_FN(acos) REF_FN(asinh) REF_FN(acosh) REF_FN(atanh)
REF_FN(sin_pi) REF_FN(cos_pi) REF_FN(tan_pi)

static void
ref_sin_cos(arb_t s, arb_t c, const arb_t x, slong prec)
{
    if (s != NULL) ref_sin(s, x, prec);
    if (c != NULL) ref_cos(c, x, prec);
}

static void
ref_sin_cos_pi(arb_t s, arb_t c, const arb_t x, slong prec)
{
    if (s != NULL) ref_sin_pi(s, x, prec);
    if (c != NULL) ref_cos_pi(c, x, prec);
}

static void
ref_sin_cos_pi_fmpq(arb_t s, arb_t c, const fmpq_t x, slong prec)
{
    arb_t t;
    arb_init(t);
    arb_set_fmpq(t, x, prec + REF_GUARD + 10);
    ref_sin_cos_pi(s, c, t, prec);
    arb_clear(t);
}

/* constants */
static void
ref_const(arb_t res, int (* f)(mpfr_ptr, mpfr_rnd_t), slong prec)
{
    mpfr_t v;
    int tern;
    ref_init();
    mpfr_init2(v, prec + REF_GUARD);
    tern = f(v, MPFR_RNDN);
    ref_set_mpfr(res, v, tern);
    mpfr_clear(v);
}

static int ref_zeta3_fn(mpfr_ptr r, mpfr_rnd_t rnd) { return mpfr_zeta_ui(r, 3, rnd); }

static void ref_const_pi(arb_t res, slong prec) { ref_const(res, mpfr_const_pi, prec); }
static void ref_const_log2(arb_t res, slong prec) { ref_const(res, mpfr_const_log2, prec); }
/* exp(x), beyond 2^24 as 2^k exp(x - k log 2), within MPFR's exponent
   range also where that is the range of a 32-bit long */
static void
ref_exp(arb_t res, const arb_t x, slong prec)
{
    if (arb_is_finite(x) && arf_cmpabs_2exp_si(arb_midref(x), 24) >= 0)
    {
        arb_t L, r;
        fmpz_t k;
        slong wp = prec + 30 + arf_abs_bound_lt_2exp_si(arb_midref(x));

        arb_init(L);
        arb_init(r);
        fmpz_init(k);
        ref_const_log2(L, wp);
        arb_div(r, x, L, wp);
        arf_get_fmpz(k, arb_midref(r), ARF_RND_NEAR);
        arb_mul_fmpz(r, L, k, wp);
        arb_sub(r, x, r, wp);
        ref_eval_spec(res, r, prec + 10, &ref_spec_exp);
        arb_mul_2exp_fmpz(res, res, k);
        arb_clear(L);
        arb_clear(r);
        fmpz_clear(k);
    }
    else
        ref_eval_spec(res, x, prec, &ref_spec_exp);
}

/* expm1, sinh and cosh beyond 2^24 from exp */
static void
ref_expm1(arb_t res, const arb_t x, slong prec)
{
    if (arb_is_finite(x) && arf_cmpabs_2exp_si(arb_midref(x), 24) >= 0)
    {
        ref_exp(res, x, prec + 10);
        arb_sub_ui(res, res, 1, prec + 10);
    }
    else
        ref_eval_spec(res, x, prec, &ref_spec_expm1);
}

static void
ref_sinh_cosh_big(arb_t res, const arb_t x, slong prec, int cosh)
{
    arb_t t, u;
    arb_init(t);
    arb_init(u);
    ref_exp(t, x, prec + 10);
    arb_neg(u, x);
    ref_exp(u, u, prec + 10);
    if (cosh)
        arb_add(res, t, u, prec + 10);
    else
        arb_sub(res, t, u, prec + 10);
    arb_mul_2exp_si(res, res, -1);
    arb_clear(t);
    arb_clear(u);
}

static void
ref_sinh(arb_t res, const arb_t x, slong prec)
{
    if (arb_is_finite(x) && arf_cmpabs_2exp_si(arb_midref(x), 24) >= 0)
        ref_sinh_cosh_big(res, x, prec, 0);
    else
        ref_eval_spec(res, x, prec, &ref_spec_sinh);
}

static void
ref_cosh(arb_t res, const arb_t x, slong prec)
{
    if (arb_is_finite(x) && arf_cmpabs_2exp_si(arb_midref(x), 24) >= 0)
        ref_sinh_cosh_big(res, x, prec, 1);
    else
        ref_eval_spec(res, x, prec, &ref_spec_cosh);
}

static void ref_const_euler(arb_t res, slong prec) { ref_const(res, mpfr_const_euler, prec); }
static void ref_const_catalan(arb_t res, slong prec) { ref_const(res, mpfr_const_catalan, prec); }
static void ref_const_apery(arb_t res, slong prec) { ref_const(res, ref_zeta3_fn, prec); }

static void
ref_log_ui(arb_t res, ulong n, slong prec)
{
    arb_t t;
    arb_init(t);
    arb_set_ui(t, n);
    ref_log(res, t, prec);
    arb_clear(t);
}

/* atan(p/q) resp. atanh(p/q), with the signature of arb_atan_frac_bsplit:
   the quotient as a ball at prec + REF_GUARD + 10 bits, its error
   propagated by the derivative bound */
static void
ref_atan_frac(arb_t res, const fmpz_t p, const fmpz_t q, int hyperbolic, slong prec)
{
    arb_t t;
    arb_init(t);
    arb_fmpz_div_fmpz(t, p, q, prec + REF_GUARD + 10);
    if (hyperbolic)
        ref_atanh(res, t, prec);
    else
        ref_atan(res, t, prec);
    arb_clear(t);
}

#endif
