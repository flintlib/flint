/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb/elem.h"

/*
    sin, cos, tan of x and of pi x (arb/elem.h).

    The point evaluation _arb_sin_cos_arf handles the special values,
    the tiny arguments (sin x = x within |x|^3 / 6, cos x = 1 within
    x^2 / 2), the register fast paths at one and two limbs, and else
    the kernels _mp_real_sin_cos_fp, _mp_real_tan_fp and their variants
    for pi x (which carry the exact values at the multiples of 1/2 and
    the tiny products pi x); a tangent given as a reciprocal is inverted
    as a ball.

    The ball functions evaluate the midpoint at a precision bounded by
    the accuracy of the input and propagate the radius through the
    derivative: |sin(x + t) - sin x| <= |cos x| |t| + t^2 / 2 and
    likewise for the cosine, using the computed cosine resp. sine, with
    pi |t| in place of |t| for the functions of pi x.  Wide balls (a
    radius above 2^-24, or beyond the magnitude limit) go to
    _arb_sin_cos_wide, which bounds the range on the interval.
*/

#define MAGLIM(prec) MP_REAL_TRIG_MAGLIM(prec)

/* the ball of (-1)^neg / Y for the frame (y, w + 1) within err units */
static void
_arb_set_mpn_err_recip(arb_t res, nn_ptr y, slong w, ulong err, int neg, slong prec)
{
    _arb_set_mpn_err(res, y, w + 1, -FLINT_BITS * w, err, err, neg, prec + 4);
    arb_ui_div(res, 1, res, prec);
}

/* an upper bound for the frame (y, len) 2^e within err units (err far
   below 2^30) */
static void
_mag_set_frame(mag_t m, nn_srcptr y, slong len, slong e, ulong err)
{
    slong i;
    unsigned int b;
    ulong t;

    for (i = len - 1; i >= 0 && y[i] == 0; i--)
        ;
    if (i < 0)
    {
        mag_set_ui_2exp_si(m, err, e);
        return;
    }
    if (i == 0)
    {
        mag_set_ui_2exp_si(m, y[0], e);
        mag_add_ui_2exp_si(m, m, err + 1, e);
        return;
    }

    /* the top FLINT_BITS bits of (y, i + 1), truncated, then two units of
       its top 30 bits (the truncation and the error) */
    b = FLINT_BIT_COUNT(y[i]);
    t = (b == FLINT_BITS) ? y[i] : ((y[i] << (FLINT_BITS - b)) | (y[i - 1] >> b));
    mag_set_ui_2exp_si(m, (t >> (FLINT_BITS - MAG_BITS)) + 2, e + FLINT_BITS * i + b - MAG_BITS);
}

void
_arb_sin_cos_arf(arb_t rs, arb_t rc, arb_t rt, mag_t ms, mag_t mc, const arf_t x, int pi, slong prec, slong maglim)
{
    nn_srcptr d;
    slong xn, e, nk, w, k, eshift = 0;
    fmpz_t xe;
    int xsgn;

    if (ms != NULL) mag_one(ms);
    if (mc != NULL) mag_one(mc);

    if (arf_is_special(x))
    {
        if (arf_is_zero(x))
        {
            if (rs != NULL) arb_zero(rs);
            if (rc != NULL) arb_one(rc);
            if (rt != NULL) arb_zero(rt);
            if (ms != NULL) mag_zero(ms);
        }
        else if (arf_is_nan(x) || pi)
        {
            if (rs != NULL) arb_indeterminate(rs);
            if (rc != NULL) arb_indeterminate(rc);
            if (rt != NULL) arb_indeterminate(rt);
        }
        else
        {
            if (rs != NULL) arb_zero_pm_one(rs);
            if (rc != NULL) arb_zero_pm_one(rc);
            if (rt != NULL) arb_zero_pm_inf(rt);
        }
        return;
    }

    xsgn = ARF_SGNBIT(x);
    e = ARF_EXP(x);
    ARF_GET_MPN_READONLY(d, xn, x);

    if (FLINT_UNLIKELY(COEFF_IS_MPZ(e) || (e > maglim
            && !(pi && e <= FLINT_BITS * xn))))
    {
        if (pi && !(COEFF_IS_MPZ(e) && fmpz_sgn(ARF_EXPREF(x)) < 0))
        {
            /* (the reduction mod 2 is exact: no magnitude limit) beyond
               the mantissa an even integer, exact below */
            e = FLINT_BITS * xn + 1;
        }
        else if (COEFF_IS_MPZ(e) && fmpz_sgn(ARF_EXPREF(x)) < 0)
        {
            /* tiny with a bignum exponent: evaluated as if the exponent
               were e, the difference restored on the results of the pi
               functions (which carry e) */
            e = -prec - FLINT_BITS * (xn + 1) - 100;
            eshift = 1;
            fmpz_init_set(xe, ARF_EXPREF(x));
        }
        else
        {
            if (rs != NULL) arb_zero_pm_one(rs);
            if (rc != NULL) arb_zero_pm_one(rc);
            if (rt != NULL) arb_zero_pm_inf(rt);
            return;
        }
    }

    /* the functions of pi x at a multiple of 1/2: exact values (a pole
       of the tangent at a half-integer) */
    if (pi && e >= 0 && !eshift)
    {
        /* the bit of weight 1/2 is bit FLINT_BITS xn - e - 1 of the
           mantissa; an integer (half-integer) has nothing below it
           (below and including it) */
        slong b = FLINT_BITS * xn - e - 1;      /* (>= -1) */
        int half, k_odd;
        if (b < 0 || (flint_mpn_zero_p(d, b / FLINT_BITS) && (d[b / FLINT_BITS] & ((UWORD(1) << (b % FLINT_BITS)) - 1)) == 0))
        {
            half = (b >= 0) && ((d[b / FLINT_BITS] >> (b % FLINT_BITS)) & 1);
            k_odd = (b + 1 >= 0) && (b + 1 < FLINT_BITS * xn) && ((d[(b + 1) / FLINT_BITS] >> ((b + 1) % FLINT_BITS)) & 1);
            if (!half)
            {
                /* sin = 0, cos = (-1)^k, tan = 0 */
                if (rs != NULL) arb_zero(rs);
                if (rt != NULL) arb_zero(rt);
                if (rc != NULL) arb_set_si(rc, k_odd ? -1 : 1);
                if (ms != NULL) mag_zero(ms);
            }
            else
            {
                /* sin = (-1)^k (x > 0), cos = 0, tan a pole */
                if (rs != NULL) arb_set_si(rs, (k_odd ^ xsgn) ? -1 : 1);
                if (rc != NULL) arb_zero(rc);
                if (rt != NULL) arb_indeterminate(rt);
                if (mc != NULL) mag_zero(mc);
            }
            return;
        }
    }

    /* tiny: sin x = x within |x|^3 / 6 < 2^(3e - 2), tan x = x within
       |x|^3 / 2 < 2^(3e - 1), cos x = 1 within x^2 / 2 < 2^(2e - 1)
       (the functions of pi x: the kernels), with the true exponent of x
       (a bignum: the fake e only for the pi functions) */
    if (!pi && 2 * e <= -prec - 2)
    {
        fmpz_t t3, t2;
        fmpz_init(t3);
        fmpz_init(t2);

        /* (everything read from x before the outputs, which may alias it,
           are written) */
        fmpz_mul_ui(t3, ARF_EXPREF(x), 3);
        fmpz_sub_ui(t3, t3, 1);
        fmpz_mul_ui(t2, ARF_EXPREF(x), 2);
        fmpz_sub_ui(t2, t2, 1);
        if (ms != NULL)
            arf_get_mag(ms, x);

        if (rs != NULL || rt != NULL)
        {
            arb_ptr r = (rs != NULL) ? rs : rt;
            int inexact = arf_set_round(arb_midref(r), x, prec, ARB_RND);
            mag_one(arb_radref(r));
            mag_mul_2exp_fmpz(arb_radref(r), arb_radref(r), t3);
            if (inexact)
                arf_mag_add_ulp(arb_radref(r), arb_radref(r), arb_midref(r), prec);
            if (rs != NULL && rt != NULL)
                arb_set(rt, rs);
        }
        if (rc != NULL)
        {
            arf_one(arb_midref(rc));
            mag_one(arb_radref(rc));
            mag_mul_2exp_fmpz(arb_radref(rc), arb_radref(rc), t2);
        }

        fmpz_clear(t3);
        fmpz_clear(t2);
        if (eshift)
            fmpz_clear(xe);
        return;
    }

#if FLINT_BITS == 64
    if (xn <= 2 && prec <= 2 * FLINT_BITS && !eshift && e < FLINT_BITS)
    {
        /* the reductions and kernels of arf/elem.h in registers (a
           one-limb mantissa at more than FLINT_BITS bits padded to two) */
        int a = 0, sneg = 0, ok = 1;
        ulong v2 = 0, v1 = 0, v0 = 0, ev = 0;

        if (xn == 1 && prec <= FLINT_BITS)
        {
            if (pi)
                ok = (e >= -31 && _arf_small_trig_pi_reduce_1(&a, &sneg, &v1, &v0, d[0], e) == 0), ev = 2;
            else if (e < 0 || (e == 0 && rt == NULL))
            {
                ok = (e >= -32);
                v1 = (e == 0) ? d[0] : (d[0] >> (-e));
                v0 = (e == 0) ? 0 : (d[0] << (FLINT_BITS + e));
                ev = 0;
            }
            else
                ok = (_arf_small_trig_reduce_1(&a, &sneg, &v1, &v0, d[0], e) == 0), ev = 2;

            if (ok)
            {
                ulong s1, s0, g1, g0, c2, c1, c0;
                ulong err = MP_REAL_SMALL_SIN_COS_ERR_2 + ev;
                int nsv, ncv;

                _mp_real_small_sin_cos_g_2(&s1, &s0, &g1, &g0, v1, v0);
                sub_dddmmmsss(c2, c1, c0, UWORD(1), UWORD(0), UWORD(0), UWORD(0), g1, g0);
                nsv = (a & 1) ? (a == 3) : (sneg ^ (a == 2));
                ncv = (a & 1) ? (sneg ^ (a == 1)) : (a == 2);

                if (rs != NULL)
                {
                    if (!(a & 1))
                        _arb_small_set_1(rs, 0, s1, s0, 0, err, nsv ^ xsgn, prec);
                    else
                        _arb_small_set_1(rs, c2, c1, c0, 0, err, nsv ^ xsgn, prec);
                }
                if (rc != NULL)
                {
                    if (a & 1)
                        _arb_small_set_1(rc, 0, s1, s0, 0, err, ncv, prec);
                    else
                        _arb_small_set_1(rc, c2, c1, c0, 0, err, ncv, prec);
                }
                if (ms != NULL || mc != NULL)
                {
                    ulong t[3];
                    t[0] = s0; t[1] = s1; t[2] = 0;
                    if (((a & 1) == 0) ? (ms != NULL) : (mc != NULL))
                        _mag_set_frame((a & 1) ? mc : ms, t, 3, -2 * FLINT_BITS, err);
                    t[0] = c0; t[1] = c1; t[2] = c2;
                    if (((a & 1) == 0) ? (mc != NULL) : (ms != NULL))
                        _mag_set_frame((a & 1) ? ms : mc, t, 3, -2 * FLINT_BITS, err);
                }
                if (rt != NULL)
                {
                    /* the quotient in registers (as arf/sin_cos.c): tan v
                       = sin v / cos v < 1 for a even, else cot v */
                    ulong q[2];
                    if (!(a & 1) && (s1 < c1 || (s1 == c1 && s0 < c0)))
                    {
                        _mp_real_divq_4_2(q, s1, s0, UWORD(0), UWORD(0), c1, c0);
                        _arb_small_set_1(rt, 0, q[1], q[0], 0, 3 * err + 1, sneg ^ xsgn, prec);
                    }
                    else if (a & 1)
                    {
                        unsigned int sh = flint_clz(s1);
                        ulong t1 = s1, t0 = s0;
                        if (sh != 0)
                        {
                            t1 = (s1 << sh) | (s0 >> (FLINT_BITS - sh));
                            t0 = s0 << sh;
                        }
                        _mp_real_divq_4_2(q, c1 >> 1, (c1 << (FLINT_BITS - 1)) | (c0 >> 1),
                            c0 << (FLINT_BITS - 1), UWORD(0), t1, t0);
                        _arb_small_set_1(rt, 0, q[1], q[0], sh + 1,
                            err * ((UWORD(2) << sh) + 2) + 1, 1 ^ sneg ^ xsgn, prec);
                    }
                    else
                    {
                        arb_t u;
                        arb_init(u);
                        _arb_small_set_1(rt, 0, s1, s0, 0, err, sneg ^ xsgn, prec + 4);
                        _arb_small_set_1(u, c2, c1, c0, 0, err, 0, prec + 4);
                        arb_div(rt, rt, u, prec);
                        arb_clear(u);
                    }
                }
                return;
            }
        }
        else
        {
            ulong d1 = d[xn - 1], d0 = (xn == 2) ? d[0] : 0;

            if (pi)
                ok = (e >= -39 && _arf_small_trig_pi_reduce_2(&a, &sneg, &v2, &v1, &v0, d1, d0, e) == 0), ev = 8;
            else if (e < 0 || (e == 0 && rt == NULL))
            {
                ok = (e >= -40);
                if (e == 0)
                {
                    v2 = d1; v1 = d0; v0 = 0;
                }
                else
                {
                    unsigned int sh = -e;
                    v2 = d1 >> sh;
                    v1 = (d1 << (FLINT_BITS - sh)) | (d0 >> sh);
                    v0 = d0 << (FLINT_BITS - sh);
                }
                ev = 0;
            }
            else
                ok = (_arf_small_trig_reduce_2(&a, &sneg, &v2, &v1, &v0, d1, d0, e) == 0), ev = 2;

            if (ok)
            {
                ulong s2, s1, s0, g2, g1, g0, c3, c2, c1, c0;
                ulong err = MP_REAL_SMALL_SIN_COS_ERR_3 + ev;
                int nsv, ncv;

                _mp_real_small_sin_cos_g_3(&s2, &s1, &s0, &g2, &g1, &g0, v2, v1, v0);
                sub_ddddmmmmssss(c3, c2, c1, c0, UWORD(1), UWORD(0), UWORD(0), UWORD(0), UWORD(0), g2, g1, g0);
                nsv = (a & 1) ? (a == 3) : (sneg ^ (a == 2));
                ncv = (a & 1) ? (sneg ^ (a == 1)) : (a == 2);

                if (rs != NULL)
                {
                    if (!(a & 1))
                        _arb_small_set_2(rs, 0, s2, s1, s0, 0, err, nsv ^ xsgn, prec);
                    else
                        _arb_small_set_2(rs, c3, c2, c1, c0, 0, err, nsv ^ xsgn, prec);
                }
                if (rc != NULL)
                {
                    if (a & 1)
                        _arb_small_set_2(rc, 0, s2, s1, s0, 0, err, ncv, prec);
                    else
                        _arb_small_set_2(rc, c3, c2, c1, c0, 0, err, ncv, prec);
                }
                if (ms != NULL || mc != NULL)
                {
                    ulong t[4];
                    t[0] = s0; t[1] = s1; t[2] = s2; t[3] = 0;
                    if (((a & 1) == 0) ? (ms != NULL) : (mc != NULL))
                        _mag_set_frame((a & 1) ? mc : ms, t, 4, -3 * FLINT_BITS, err);
                    t[0] = c0; t[1] = c1; t[2] = c2; t[3] = c3;
                    if (((a & 1) == 0) ? (mc != NULL) : (ms != NULL))
                        _mag_set_frame((a & 1) ? ms : mc, t, 4, -3 * FLINT_BITS, err);
                }
                if (rt != NULL)
                {
                    /* the quotient in registers (as arf/sin_cos.c) */
                    ulong w2, w1, w0;
                    if (!(a & 1) && (s2 < c2 || (s2 == c2 && (s1 < c1 || (s1 == c1 && s0 < c0)))))
                    {
                        _mp_real_small_divq_3(&w2, &w1, &w0, s2, s1, s0, c2, c1, c0);
                        _arb_small_set_2(rt, 0, w2, w1, w0, 0, 3 * err + MP_REAL_SMALL_DIVQ_3_ERR, sneg ^ xsgn, prec);
                    }
                    else if (a & 1)
                    {
                        unsigned int sh = flint_clz(s2);
                        ulong t2 = s2, t1 = s1, t0 = s0;
                        if (sh != 0)
                        {
                            t2 = (s2 << sh) | (s1 >> (FLINT_BITS - sh));
                            t1 = (s1 << sh) | (s0 >> (FLINT_BITS - sh));
                            t0 = s0 << sh;
                        }
                        _mp_real_small_divq_3(&w2, &w1, &w0, c2 >> 1, (c2 << (FLINT_BITS - 1)) | (c1 >> 1),
                            (c1 << (FLINT_BITS - 1)) | (c0 >> 1), t2, t1, t0);
                        _arb_small_set_2(rt, 0, w2, w1, w0, sh + 1,
                            err * ((UWORD(2) << sh) + 2) + MP_REAL_SMALL_DIVQ_3_ERR + 1, 1 ^ sneg ^ xsgn, prec);
                    }
                    else
                    {
                        arb_t u;
                        arb_init(u);
                        _arb_small_set_2(rt, 0, s2, s1, s0, 0, err, sneg ^ xsgn, prec + 4);
                        _arb_small_set_2(u, c3, c2, c1, c0, 0, err, 0, prec + 4);
                        arb_div(rt, rt, u, prec);
                        arb_clear(u);
                    }
                }
                return;
            }
        }
    }
#endif

    /* small (4e <= -prec - 2): two terms (arb/elem.c), sin x = x - x^3/6
       within |x|^5 / 120 < 2^(5e - 6), tan x = x + x^3/3 within
       2 |x|^5 / 15 < 2^(5e - 2), cos x = 1 - x^2/2 within x^4 / 24 <
       2^(4e - 4) */
    if (!pi && e <= -((prec + 5) / 4))
    {
        if (ms != NULL)
            arf_get_mag(ms, x);

        /* (the outputs may alias x, whose mantissa is read for each:
           temporaries, swapped in at the end) */
        {
            arb_t ts, tt, tc;
            arb_init(ts);
            arb_init(tt);
            arb_init(tc);
            if (rs != NULL)
            {
                _arb_small_arg_series(ts, d, xn, e, xsgn, ARB_SMALL_SERIES_SIN, prec);
                arb_add_error_2exp_si(ts, 5 * e - 6);
            }
            if (rt != NULL)
            {
                _arb_small_arg_series(tt, d, xn, e, xsgn, ARB_SMALL_SERIES_TAN, prec);
                arb_add_error_2exp_si(tt, 5 * e - 2);
            }
            if (rc != NULL)
            {
                _arb_small_arg_series(tc, d, xn, e, 0, ARB_SMALL_SERIES_COS, prec);
                arb_add_error_2exp_si(tc, 4 * e - 4);
            }
            if (rs != NULL) arb_swap(rs, ts);
            if (rt != NULL) arb_swap(rt, tt);
            if (rc != NULL) arb_swap(rc, tc);
            arb_clear(ts);
            arb_clear(tt);
            arb_clear(tc);
        }
        if (eshift)
            fmpz_clear(xe);
        return;
    }


    nk = _arb_elem_limbs(prec);
    {
        nn_ptr ys, yc;
        ulong err;
        slong wmax, ylen, tlen;
        int code, recip;
        TMP_INIT;

        /* buffers for the relative accuracy near the zeros (the reduced
           argument of pi x has at most the bits of x) */
        wmax = nk + xn + prec / FLINT_BITS + 2;
        tlen = FLINT_MAX(xn, nk) + 3;
        ylen = FLINT_MAX(wmax + 3, tlen + 1);

        TMP_START;
        ys = TMP_ALLOC(2 * ylen * sizeof(ulong));
        yc = ys + ylen;

        if (rt == NULL)
        {
            if (!pi)
                code = _mp_real_sin_cos_fp(ys, yc, &w, &err, d, xn, e, _arb_elem_bits(prec), wmax, (rs != NULL) | ((rc != NULL) << 1)), k = 0;
            else
                code = _mp_real_sin_cos_pi_fp(ys, yc, &w, &k, &err, d, xn, e, _arb_elem_bits(prec), wmax, (rs != NULL) | ((rc != NULL) << 1));

            if (code & MP_REAL_TRIG_EXACT)
            {
                int a = code & 3;
                if (rs != NULL)
                {
                    if (a & 1)
                        arb_set_si(rs, ((a == 3) ^ xsgn) ? -1 : 1);
                    else
                        arb_zero(rs);
                }
                if (rc != NULL)
                {
                    if (a & 1)
                        arb_zero(rc);
                    else
                        arb_set_si(rc, (a == 2) ? -1 : 1);
                }
                if (ms != NULL && !(a & 1)) mag_zero(ms);
                if (mc != NULL && (a & 1)) mag_zero(mc);
            }
            else if (code & MP_REAL_TRIG_TINY)
            {
                /* the sine in (Y - 1, Y + 6) units, the cosine 1 from above */
                if (rs != NULL)
                {
                    _arb_set_mpn_err(rs, ys, tlen, 0, 1, 6, xsgn, prec);
                    if (eshift)
                    {
                        fmpz_t t;
                        fmpz_init(t);
                        fmpz_sub_si(t, xe, e);
                        fmpz_add_si(t, t, k);
                        arb_mul_2exp_fmpz(rs, rs, t);
                        fmpz_clear(t);
                    }
                    else
                        arb_mul_2exp_si(rs, rs, k);
                }
                if (rc != NULL)
                {
                    arf_one(arb_midref(rc));
                    mag_set_ui_2exp_si(arb_radref(rc), 1, -prec - 2);
                }
                if (ms != NULL)
                {
                    /* (the bound 2^k (Y + 6), with the true exponent) */
                    _mag_set_frame(ms, ys, tlen, 0, 6);
                    if (eshift)
                    {
                        fmpz_t t;
                        fmpz_init(t);
                        fmpz_sub_si(t, xe, e);
                        fmpz_add_si(t, t, k);
                        mag_mul_2exp_fmpz(ms, ms, t);
                        fmpz_clear(t);
                    }
                    else
                        mag_mul_2exp_si(ms, ms, k);
                }
            }
            else
            {
                if (rs != NULL)
                    _arb_set_mpn_err(rs, ys, w + 1, -FLINT_BITS * w, err, err, (code & 1) ^ xsgn, prec);
                if (rc != NULL)
                    _arb_set_mpn_err(rc, yc, w + 1, -FLINT_BITS * w, err, err, (code >> 1) & 1, prec);
                if (ms != NULL)
                    _mag_set_frame(ms, ys, w + 1, -FLINT_BITS * w, err);
                if (mc != NULL)
                    _mag_set_frame(mc, yc, w + 1, -FLINT_BITS * w, err);
            }
        }
        else
        {
            if (!pi)
                code = _mp_real_tan_fp(ys, &w, &err, &recip, d, xn, e, _arb_elem_bits(prec), wmax), k = 0;
            else
                code = _mp_real_tan_pi_fp(ys, &w, &k, &err, &recip, d, xn, e, _arb_elem_bits(prec), wmax);

            if (code & MP_REAL_TRIG_EXACT)
            {
                if (code & MP_REAL_TRIG_ONE)
                    arb_set_si(rt, ((code & 1) ^ xsgn) ? -1 : 1);
                else if (code & 1)
                    arb_indeterminate(rt);    /* a pole */
                else
                    arb_zero(rt);
            }
            else if (code & MP_REAL_TRIG_TINY)
            {
                /* the tangent in [Y, Y + 6) units */
                _arb_set_mpn_err(rt, ys, tlen, 0, 0, 6, xsgn, prec);
                if (eshift)
                {
                    fmpz_t t;
                    fmpz_init(t);
                    fmpz_sub_si(t, xe, e);
                    fmpz_add_si(t, t, k);
                    arb_mul_2exp_fmpz(rt, rt, t);
                    fmpz_clear(t);
                }
                else
                    arb_mul_2exp_si(rt, rt, k);
            }
            else if (!recip)
                _arb_set_mpn_err(rt, ys, w + 1, -FLINT_BITS * w, err, err, code ^ xsgn, prec);
            else
                _arb_set_mpn_err_recip(rt, ys, w, err, code ^ xsgn, prec);
        }

        TMP_END;
    }

    if (eshift)
        fmpz_clear(xe);
}

/* the propagated error d r + r^2 / 2 for the derivative bound d and the
   radius r, or pi r for the functions of pi x (pi < 2^2); r is finite
   and below 1, so the fast mag arithmetic applies once the exponents
   are in the lagom range */
static void
_mag_trig_propagate(mag_t res, const mag_t d, const mag_t r, int pi)
{
    mag_t t, u;

    mag_init(t);
    mag_init(u);

    if (pi)
        mag_mul_2exp_si(t, r, 2);
    else
        mag_set(t, r);

    if (MAG_IS_LAGOM(t) && MAG_IS_LAGOM(d) && MAG_IS_LAGOM(res))
    {
        mag_fast_mul(u, d, t);
        mag_add(res, res, u);
        mag_fast_mul(t, t, t);
        mag_fast_mul_2exp_si(t, t, -1);
        mag_add(res, res, t);
    }
    else
    {
        mag_mul(u, d, t);
        mag_mul(t, t, t);
        mag_mul_2exp_si(t, t, -1);
        mag_add(u, u, t);
        mag_add(res, res, u);
    }

    mag_clear(t);
    mag_clear(u);
}

/* sin, cos (either NULL) of the ball x = xmid +- xrad, or of pi x */
static void
_arb_sin_cos(arb_t zsin, arb_t zcos, const arf_t x, const mag_t xrad, int pi, slong prec)
{
    slong maglim = MAGLIM(prec), exp, radexp;

    if (mag_is_zero(xrad))
    {
        _arb_sin_cos_arf(zsin, zcos, NULL, NULL, NULL, x, pi, prec, maglim);
        return;
    }

    if (!arf_is_finite(x) || !mag_is_finite(xrad))
    {
        if (arf_is_nan(x) || (pi && !mag_is_finite(xrad)))
        {
            if (zsin != NULL) arb_indeterminate(zsin);
            if (zcos != NULL) arb_indeterminate(zcos);
        }
        else
        {
            if (zsin != NULL) arb_zero_pm_one(zsin);
            if (zcos != NULL) arb_zero_pm_one(zcos);
        }
        return;
    }

    /* the functions of pi x are 2-periodic: a midpoint |x| >= 2 is
       reduced exactly into [-1, 1] (an even integer, such as any x beyond
       its mantissa, to 0), so that neither the wide algorithm nor the
       magnitude limit sees more than the reduced value */
    if (pi && arf_cmpabs_2exp_si(x, 1) >= 0)
    {
        arf_t v;
        arf_init(v);
        if (!COEFF_IS_MPZ(ARF_EXP(x)) && ARF_EXP(x) <= (slong) arf_bits(x))
        {
            fmpz_t n;
            fmpz_init(n);
            arf_mul_2exp_si(v, x, -1);
            arf_get_fmpz(n, v, ARF_RND_NEAR);
            arf_sub_fmpz(v, v, n, ARF_PREC_EXACT, ARF_RND_DOWN);
            arf_mul_2exp_si(v, v, 1);
            fmpz_clear(n);
        }
        _arb_sin_cos(zsin, zcos, v, xrad, pi, prec);
        arf_clear(v);
        return;
    }

    radexp = MAG_EXP(xrad);
    if (COEFF_IS_MPZ(radexp))
        radexp = (fmpz_sgn(MAG_EXPREF(xrad)) < 0) ? -prec - 100 : 0;
    exp = COEFF_IS_MPZ(ARF_EXP(x)) ? ((fmpz_sgn(ARF_EXPREF(x)) < 0) ? -prec - 100 : maglim + 1) : ARF_EXP(x);

    /* wide balls (r >= 2^-24, resp. 2^-26 for pi x) and beyond the
       magnitude limit: the range on the interval */
    if (radexp >= (pi ? -25 : -23) || exp > maglim)
    {
        if (!pi)
        {
            _arb_sin_cos_wide(zsin, zcos, x, xrad, prec);
        }
        else
        {
            /* pi x as a ball (the radius pi r within a few bits), then
               the wide algorithm */
            arb_t t, u;
            slong wp = FLINT_MIN(prec, 64) + 10;
            arb_init(t);
            arb_init(u);
            arf_set(arb_midref(u), x);
            mag_set(arb_radref(u), xrad);
            arb_const_pi(t, wp);
            arb_mul(t, t, u, wp);
            _arb_sin_cos_wide(zsin, zcos, arb_midref(t), arb_radref(t), prec);
            arb_clear(t);
            arb_clear(u);
        }
        return;
    }

    /* the precision bounded by the accuracy of the input: r 2^-radexp
       <= 1 units below the point */
    /* near x = 0 the error of the cosine is quadratic */
    if (zcos != NULL && exp < -2 && !pi)
        prec = FLINT_MIN(prec, 20 - FLINT_MAX(exp, radexp) - radexp);
    else
        prec = FLINT_MIN(prec, 20 - radexp);
    prec = FLINT_MAX(prec, 2);

    /* the derivatives from the bounds on the other function (the outputs
       may alias x: the radius copied) */
    {
        mag_t r, ms, mc;
        mag_init_set(r, xrad);
        mag_init(ms);
        mag_init(mc);
        _arb_sin_cos_arf(zsin, zcos, NULL, (zcos != NULL) ? ms : NULL, (zsin != NULL) ? mc : NULL, x, pi, prec, maglim);

        if (zsin != NULL)
            _mag_trig_propagate(arb_radref(zsin), mc, r, pi);
        if (zcos != NULL)
            _mag_trig_propagate(arb_radref(zcos), ms, r, pi);

        mag_clear(r);
        mag_clear(ms);
        mag_clear(mc);
    }
}

void
arb_sin_cos(arb_t s, arb_t c, const arb_t x, slong prec)
{
    _arb_sin_cos(s, c, arb_midref(x), arb_radref(x), 0, prec);
}

void
arb_sin(arb_t s, const arb_t x, slong prec)
{
    _arb_sin_cos(s, NULL, arb_midref(x), arb_radref(x), 0, prec);
}

void
arb_cos(arb_t c, const arb_t x, slong prec)
{
    _arb_sin_cos(NULL, c, arb_midref(x), arb_radref(x), 0, prec);
}

void
arb_sin_cos_pi(arb_t s, arb_t c, const arb_t x, slong prec)
{
    _arb_sin_cos(s, c, arb_midref(x), arb_radref(x), 1, prec);
}

void
arb_sin_pi(arb_t s, const arb_t x, slong prec)
{
    _arb_sin_cos(s, NULL, arb_midref(x), arb_radref(x), 1, prec);
}

void
arb_cos_pi(arb_t c, const arb_t x, slong prec)
{
    _arb_sin_cos(NULL, c, arb_midref(x), arb_radref(x), 1, prec);
}

/* tan of the ball x, or of pi x: the point evaluation for an exact x,
   else sin / cos */
static void
_arb_tan(arb_t y, const arb_t x, int pi, slong prec)
{
    if (mag_is_zero(arb_radref(x)))
    {
        _arb_sin_cos_arf(NULL, NULL, y, NULL, NULL, arb_midref(x), pi, prec, MAGLIM(prec));
    }
    else
    {
        arb_t u;
        arb_init(u);
        _arb_sin_cos(y, u, arb_midref(x), arb_radref(x), pi, prec + 4);
        arb_div(y, y, u, prec);
        arb_clear(u);
    }
}

void
arb_tan(arb_t y, const arb_t x, slong prec)
{
    _arb_tan(y, x, 0, prec);
}

void
arb_tan_pi(arb_t y, const arb_t x, slong prec)
{
    _arb_tan(y, x, 1, prec);
}
