/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arf/elem.h"

/*
    arf_sin, arf_cos, arf_sin_cos, arf_tan and the functions of pi x
    (arf/elem.h explains the rounding).

    - tiny x: sin x = x (1 - delta) and tan x = x (1 + delta) with
      delta < x^2 / 3, cos x = 1 - delta with delta < x^2 / 2, on the known
      sides (_arf_set_tiny_rel_err, _arf_one_plus_tiny); for the functions
      of pi x the kernels return the product pi x with its bounds.
    - the register fast paths of arf/elem.h for mantissas of one and two
      limbs at prec <= FLINT_BITS resp. 2 FLINT_BITS (64-bit).
    - the kernels _mp_real_sin_cos_fp, _mp_real_tan_fp and the pi variants
      at increasing guard bits (and buffers, for the relative accuracy
      near the zeros) until the rounding is decided; a cotangent takes
      one approximate reciprocal of the normalised tangent
      (_arf_set_mpn_err_quot).
    - arguments beyond 2^ARF_TRIG_MAGLIM(prec) give nan.
*/

/* the fixed-point result (y, w + 1) within err ulps into res */
static int
_set(arf_t res, nn_ptr y, slong w, ulong err, int neg, slong prec, arf_rnd_t rnd)
{
    return _arf_set_mpn_err(res, y, w + 1, -FLINT_BITS * w, ARF_ELEM_EB(err), ARF_ELEM_EB(err), neg, prec, rnd);
}

#if FLINT_BITS == 64
/* the outputs from v at two limbs (v within ev ulps): sin, cos (either
   NULL) or tan */
static int
_small_out_1(arf_t rs, arf_t rc, arf_t rt, int a, int sneg, ulong v1, ulong v0, ulong ev, int xsgn, slong prec, arf_rnd_t rnd)
{
    ulong s1, s0, g1, g0, c2, c1, c0, q[2];
    ulong err = MP_REAL_SMALL_SIN_COS_ERR_2 + ev;
    int nsv, ncv, status;

    _mp_real_small_sin_cos_g_2(&s1, &s0, &g1, &g0, v1, v0);

    if (rt != NULL)
    {
        /* cos v = 1 - g < 1 (g > 0 for v >= 2^-32) */
        sub_ddmmss(c1, c0, UWORD(0), UWORD(0), g1, g0);

        if (!(a & 1))
        {
            /* tan v = sin v / cos v < 1 */
            if (s1 > c1 || (s1 == c1 && s0 >= c0))
                return ARF_ELEM_UNCERTAIN;
            _mp_real_divq_4_2(q, s1, s0, UWORD(0), UWORD(0), c1, c0);
            return _arf_small_set_1(rt, 0, q[1], q[0], 0, 3 * err + 1, sneg ^ xsgn, prec, rnd);
        }
        else
        {
            /* cot v = 2^(sh+1) cos v / (2 S), S = sin v 2^sh */
            unsigned int sh = flint_clz(s1);
            if (sh != 0)
            {
                s1 = (s1 << sh) | (s0 >> (FLINT_BITS - sh));
                s0 = s0 << sh;
            }
            _mp_real_divq_4_2(q, c1 >> 1, (c1 << (FLINT_BITS - 1)) | (c0 >> 1),
                c0 << (FLINT_BITS - 1), UWORD(0), s1, s0);
            return _arf_small_set_1(rt, 0, q[1], q[0], sh + 1,
                err * ((UWORD(2) << sh) + 2) + 1, 1 ^ sneg ^ xsgn, prec, rnd);
        }
    }

    /* the signs of sin |x| and cos |x|; sin |x| from sin v (a even) or
       cos v (a odd), and conversely */
    nsv = (a & 1) ? (a == 3) : (sneg ^ (a == 2));
    ncv = (a & 1) ? (sneg ^ (a == 1)) : (a == 2);
    status = 1;

    if (rs != NULL)
    {
        if (!(a & 1))
            status = _arf_small_set_1(rs, 0, s1, s0, 0, err, nsv ^ xsgn, prec, rnd);
        else
        {
            sub_dddmmmsss(c2, c1, c0, UWORD(1), UWORD(0), UWORD(0), UWORD(0), g1, g0);
            status = _arf_small_set_1(rs, c2, c1, c0, 0, err, nsv ^ xsgn, prec, rnd);
        }
        if (status == ARF_ELEM_UNCERTAIN)
            return status;
    }
    if (rc != NULL)
    {
        if (a & 1)
            status = _arf_small_set_1(rc, 0, s1, s0, 0, err, ncv, prec, rnd);
        else
        {
            sub_dddmmmsss(c2, c1, c0, UWORD(1), UWORD(0), UWORD(0), UWORD(0), g1, g0);
            status = _arf_small_set_1(rc, c2, c1, c0, 0, err, ncv, prec, rnd);
        }
    }

    return status;
}

static int
_small_out_2(arf_t rs, arf_t rc, arf_t rt, int a, int sneg, ulong v2, ulong v1, ulong v0, ulong ev, int xsgn, slong prec, arf_rnd_t rnd)
{
    ulong s2, s1, s0, g2, g1, g0, c3, c2, c1, c0, w2, w1, w0;
    ulong err = MP_REAL_SMALL_SIN_COS_ERR_3 + ev;
    int nsv, ncv, status;

    _mp_real_small_sin_cos_g_3(&s2, &s1, &s0, &g2, &g1, &g0, v2, v1, v0);

    if (rt != NULL)
    {
        sub_dddmmmsss(c2, c1, c0, UWORD(0), UWORD(0), UWORD(0), g2, g1, g0);

        if (!(a & 1))
        {
            if (s2 > c2 || (s2 == c2 && (s1 > c1 || (s1 == c1 && s0 >= c0))))
                return ARF_ELEM_UNCERTAIN;
            _mp_real_small_divq_3(&w2, &w1, &w0, s2, s1, s0, c2, c1, c0);
            return _arf_small_set_2(rt, 0, w2, w1, w0, 0,
                3 * err + MP_REAL_SMALL_DIVQ_3_ERR, sneg ^ xsgn, prec, rnd);
        }
        else
        {
            unsigned int sh = flint_clz(s2);
            if (sh != 0)
            {
                s2 = (s2 << sh) | (s1 >> (FLINT_BITS - sh));
                s1 = (s1 << sh) | (s0 >> (FLINT_BITS - sh));
                s0 = s0 << sh;
            }
            /* cos v / 2, the lowest bit dropped (half a unit) */
            _mp_real_small_divq_3(&w2, &w1, &w0, c2 >> 1, (c2 << (FLINT_BITS - 1)) | (c1 >> 1),
                (c1 << (FLINT_BITS - 1)) | (c0 >> 1), s2, s1, s0);
            return _arf_small_set_2(rt, 0, w2, w1, w0, sh + 1,
                err * ((UWORD(2) << sh) + 2) + MP_REAL_SMALL_DIVQ_3_ERR + 1, 1 ^ sneg ^ xsgn, prec, rnd);
        }
    }

    nsv = (a & 1) ? (a == 3) : (sneg ^ (a == 2));
    ncv = (a & 1) ? (sneg ^ (a == 1)) : (a == 2);
    status = 1;

    if (rs != NULL)
    {
        if (!(a & 1))
            status = _arf_small_set_2(rs, 0, s2, s1, s0, 0, err, nsv ^ xsgn, prec, rnd);
        else
        {
            sub_ddddmmmmssss(c3, c2, c1, c0, UWORD(1), UWORD(0), UWORD(0), UWORD(0), UWORD(0), g2, g1, g0);
            status = _arf_small_set_2(rs, c3, c2, c1, c0, 0, err, nsv ^ xsgn, prec, rnd);
        }
        if (status == ARF_ELEM_UNCERTAIN)
            return status;
    }
    if (rc != NULL)
    {
        if (a & 1)
            status = _arf_small_set_2(rc, 0, s2, s1, s0, 0, err, ncv, prec, rnd);
        else
        {
            sub_ddddmmmmssss(c3, c2, c1, c0, UWORD(1), UWORD(0), UWORD(0), UWORD(0), UWORD(0), g2, g1, g0);
            status = _arf_small_set_2(rc, c3, c2, c1, c0, 0, err, ncv, prec, rnd);
        }
    }

    return status;
}

/* sin, cos, tan of x (pi = 0) or of pi x (pi = 1) by the register
   paths, for a mantissa of xn <= 2 limbs with prec <= 2 FLINT_BITS (a
   one-limb mantissa at more than FLINT_BITS bits padded to two) */
static int
_small(arf_t rs, arf_t rc, arf_t rt, nn_srcptr d, slong xn, slong e, int xsgn, int pi, slong prec, arf_rnd_t rnd)
{
    int a = 0, sneg = 0;
    ulong v2, v1, v0, ev;

    if (e >= FLINT_BITS)
        return ARF_ELEM_UNCERTAIN;

    if (xn == 1 && prec <= FLINT_BITS)
    {
        if (pi)
        {
            if (e < -31 || _arf_small_trig_pi_reduce_1(&a, &sneg, &v1, &v0, d[0], e) != 0)
                return ARF_ELEM_UNCERTAIN;
            ev = 2;
        }
        else if (e < 0 || (e == 0 && rt == NULL))
        {
            /* (the tangent folds |x| > pi/4) */
            if (e < -32)
                return ARF_ELEM_UNCERTAIN;
            v1 = (e == 0) ? d[0] : (d[0] >> (-e));
            v0 = (e == 0) ? 0 : (d[0] << (FLINT_BITS + e));
            ev = 0;
        }
        else
        {
            if (_arf_small_trig_reduce_1(&a, &sneg, &v1, &v0, d[0], e) != 0)
                return ARF_ELEM_UNCERTAIN;
            ev = 2;
        }

        return _small_out_1(rs, rc, rt, a, sneg, v1, v0, ev, xsgn, prec, rnd);
    }
    else
    {
        ulong d1 = d[xn - 1], d0 = (xn == 2) ? d[0] : 0;

        if (pi)
        {
            if (e < -39 || _arf_small_trig_pi_reduce_2(&a, &sneg, &v2, &v1, &v0, d1, d0, e) != 0)
                return ARF_ELEM_UNCERTAIN;
            ev = 8;
        }
        else if (e < 0 || (e == 0 && rt == NULL))
        {
            if (e < -40)
                return ARF_ELEM_UNCERTAIN;
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
        {
            if (_arf_small_trig_reduce_2(&a, &sneg, &v2, &v1, &v0, d1, d0, e) != 0)
                return ARF_ELEM_UNCERTAIN;
            ev = 2;
        }

        return _small_out_2(rs, rc, rt, a, sneg, v2, v1, v0, ev, xsgn, prec, rnd);
    }
}
#endif

/* the exponent of a result evaluated with the fake exponent e in place
   of the bignum xe of x: the difference xe - e added (res scaling like
   x) or subtracted (res scaling like 1/x) */
static void
_arf_trig_restore_exp(arf_t res, const fmpz_t xe, slong e, int recip)
{
    fmpz_t t;
    fmpz_init(t);
    fmpz_sub_si(t, xe, e);
    if (recip)
        fmpz_sub(ARF_EXPREF(res), ARF_EXPREF(res), t);
    else
        fmpz_add(ARF_EXPREF(res), ARF_EXPREF(res), t);
    fmpz_clear(t);
}

/* the variants of the outputs: cot into rt, sec into rc, csc or sinc
   into rs (the reciprocals by one approximate division of the kernels'
   frames, sinc by one division by x) */
#define TRIG_COT 1
#define TRIG_SEC 2
#define TRIG_CSC 3
#define TRIG_SINC 4

/* sin, cos (either NULL) or tan of x (pi = 0) or pi x (pi = 1) for a
   finite nonzero x, or the variant mode; the outputs may not alias x */
static int
_arf_trig(arf_t rs, arf_t rc, arf_t rt, const arf_t x, int pi, int mode, slong prec, arf_rnd_t rnd)
{
    nn_srcptr d;
    slong xn, e, nk, extra, w, k, eshift = 0, pb;
    fmpz_t xe;
    int xsgn, status;

    if (prec == ARF_PREC_EXACT)
        flint_throw(FLINT_ERROR, "arf_sin_cos: exact precision not supported\n");

    xsgn = ARF_SGNBIT(x);
    e = ARF_EXP(x);
    ARF_GET_MPN_READONLY(d, xn, x);

    if (FLINT_UNLIKELY(COEFF_IS_MPZ(e) || (e > ARF_TRIG_MAGLIM(prec)
            && !(pi && e <= FLINT_BITS * xn))))
    {
        if (pi && !(COEFF_IS_MPZ(e) && fmpz_sgn(ARF_EXPREF(x)) < 0))
        {
            /* (the reduction mod 2 is exact: no magnitude limit) beyond
               the mantissa an even integer, which the kernels return as
               MP_REAL_TRIG_EXACT with a = 0 */
            e = FLINT_BITS * xn + 1;
        }
        else if (COEFF_IS_MPZ(e) && fmpz_sgn(ARF_EXPREF(x)) < 0)
        {
            /* tiny with a bignum exponent: evaluated as if the exponent
               were e (below every threshold), the difference restored
               on the result of the pi functions (which carries e) */
            e = -prec - FLINT_BITS * (xn + 1) - 100;
            eshift = 1;
            fmpz_init_set(xe, ARF_EXPREF(x));
        }
        else
        {
            if (rs != NULL) arf_nan(rs);
            if (rc != NULL) arf_nan(rc);
            if (rt != NULL) arf_nan(rt);
            return 0;
        }
    }

    /* tiny x (for the functions of pi x the kernels handle it): sin x =
       x (1 - delta), tan x = x (1 + delta), delta < x^2/3 < 2^(2e - 1);
       cos x = 1 - delta', delta' < x^2/2 < 2^(2e - 1) */
    if (!pi && 2 * e <= -prec - 2)
    {
        /* (the variants: cot x = (1/x)(1 - delta), csc x = (1/x)(1 + delta),
           sec x = 1 + delta', sinc x = 1 - delta'') */
        status = 1;
        if (rt != NULL)
        {
            if (mode == TRIG_COT)
                status = _arf_set_tiny_recip_rel_err(rt, x, -1, 1 - 2 * e, prec, rnd);
            else
                status = _arf_set_tiny_rel_err(rt, x, 1, 1 - 2 * e, prec, rnd);
        }
        if (status != ARF_ELEM_UNCERTAIN && rc != NULL)
            status = _arf_one_plus_tiny(rc, 0, mode != TRIG_SEC, prec, rnd);
        if (status != ARF_ELEM_UNCERTAIN && rs != NULL)
        {
            if (mode == TRIG_CSC)
                status = _arf_set_tiny_recip_rel_err(rs, x, 1, 1 - 2 * e, prec, rnd);
            else if (mode == TRIG_SINC)
                status = _arf_one_plus_tiny(rs, 0, 1, prec, rnd);
            else
                status = _arf_set_tiny_rel_err(rs, x, -1, 1 - 2 * e, prec, rnd);
        }
        if (status != ARF_ELEM_UNCERTAIN)
        {
            if (eshift)
                fmpz_clear(xe);
            return status;
        }
    }

#if FLINT_BITS == 64
    if (xn <= 2 && prec <= 2 * FLINT_BITS && !eshift && mode == 0)
    {
        status = _small(rs, rc, rt, d, xn, e, xsgn, pi, prec, rnd);
        if (status != ARF_ELEM_UNCERTAIN)
            return status;
    }
#endif

    /* small x (4e <= -prec - 2) by the series in fixed point (arf/elem.c):
       sin x = x - x^3/6 within |x|^5 / 120 < 2^(5e - 6), cos x = 1 -
       x^2/2 within x^4 / 24 < 2^(4e - 4), tan x = x + x^3/3 within
       2 |x|^5 / 15 < 2^(5e - 2) */
    if (!pi && !eshift && mode == 0 && e <= -((prec + 5) / 4))
    {
        status = 1;
        if (rt != NULL)
            status = _arf_small_arg_series(rt, d, xn, e, xsgn, ARF_SMALL_SERIES_TAN, -2, prec, rnd);
        if (status != ARF_ELEM_UNCERTAIN && rs != NULL)
            status = _arf_small_arg_series(rs, d, xn, e, xsgn, ARF_SMALL_SERIES_SIN, -6, prec, rnd);
        if (status != ARF_ELEM_UNCERTAIN && rc != NULL)
            status = _arf_small_arg_series(rc, d, xn, e, 0, ARF_SMALL_SERIES_COS, -4, prec, rnd);
        if (status != ARF_ELEM_UNCERTAIN)
            return status;
    }

    for (extra = 0; ; extra = (extra == 0) ? 32 : 2 * extra)
    {
        nn_ptr ys, yc;
        ulong err;
        slong wmax, ylen, tlen;
        int code, recip;
        TMP_INIT;

        pb = _arf_elem_bits(prec, rnd, extra);
        nk = _mp_real_elem_limbs(pb);
        /* buffers for the relative accuracy near the zeros (the
           reduced argument of pi x has at most the bits of x), grown
           with the retries */
        wmax = nk + xn + (prec + extra) / FLINT_BITS + 2;
        tlen = FLINT_MAX(xn, nk) + 3;
        ylen = FLINT_MAX(wmax + 3, tlen + 1);

        TMP_START;
        ys = TMP_ALLOC(2 * ylen * sizeof(ulong));
        yc = ys + ylen;

        if (rt == NULL)
        {
            if (!pi)
            {
                code = _mp_real_sin_cos_fp(ys, yc, &w, &err, d, xn, e, pb, wmax, (rs != NULL) | ((rc != NULL) << 1));
                k = 0;
            }
            else
                code = _mp_real_sin_cos_pi_fp(ys, yc, &w, &k, &err, d, xn, e, pb, wmax, (rs != NULL) | ((rc != NULL) << 1));

            if (code & MP_REAL_TRIG_EXACT)
            {
                /* (the reciprocals of the zeros are poles: nan) */
                int a = code & 3;
                if (rs != NULL)
                {
                    if (a & 1)
                    {
                        arf_one(rs);
                        if ((a == 3) ^ xsgn)
                            arf_neg(rs, rs);
                    }
                    else if (mode == TRIG_CSC)
                        arf_nan(rs);
                    else
                        arf_zero(rs);
                }
                if (rc != NULL)
                {
                    if (a & 1)
                    {
                        if (mode == TRIG_SEC)
                            arf_nan(rc);
                        else
                            arf_zero(rc);
                    }
                    else
                    {
                        arf_one(rc);
                        if (a == 2)
                            arf_neg(rc, rc);
                    }
                }
                TMP_END;
                if (eshift)
                    fmpz_clear(xe);
                return 0;
            }

            status = 1;
            if (code & MP_REAL_TRIG_TINY)
            {
                /* (sinc of a tiny argument never gets here: it is 1 - delta
                   from the exponent, before any kernel) */
                FLINT_ASSERT(mode != TRIG_SINC);
                /* the sine in (Y - 1, Y + 6) units, the cosine 1 from above
                   (its reciprocal 1 from below); the cosecant by one
                   division of the frame */
                if (rs != NULL)
                {
                    if (mode == TRIG_CSC)
                        status = _arf_set_mpn_err_quot(rs, ys, tlen, k, 6, NULL, 0, 0, 1, xsgn, nk, prec, rnd);
                    else
                        status = _arf_set_mpn_err(rs, ys, tlen, k, 0, 3, xsgn, prec, rnd);
                    if (status != ARF_ELEM_UNCERTAIN && eshift)
                        _arf_trig_restore_exp(rs, xe, e, mode == TRIG_CSC);
                }
                if (status != ARF_ELEM_UNCERTAIN && rc != NULL)
                    status = _arf_one_plus_tiny(rc, 0, mode != TRIG_SEC, prec, rnd);
            }
            else
            {
                if (rs != NULL)
                {
                    if (mode == TRIG_CSC)
                        status = _arf_set_mpn_err_quot(rs, ys, w + 1, -FLINT_BITS * w, err, NULL, 0, 0, 1, (code & 1) ^ xsgn, nk, prec, rnd);
                    else if (mode == TRIG_SINC)
                        status = _arf_set_mpn_err_quot(rs, ys, w + 1, -FLINT_BITS * w, err, d, xn, e, 0, code & 1, nk, prec, rnd);
                    else
                        status = _set(rs, ys, w, err, (code & 1) ^ xsgn, prec, rnd);
                }
                if (status != ARF_ELEM_UNCERTAIN && rc != NULL)
                {
                    if (mode == TRIG_SEC)
                        status = _arf_set_mpn_err_quot(rc, yc, w + 1, -FLINT_BITS * w, err, NULL, 0, 0, 1, (code >> 1) & 1, nk, prec, rnd);
                    else
                        status = _set(rc, yc, w, err, (code >> 1) & 1, prec, rnd);
                }
            }
        }
        else
        {
            if (!pi)
            {
                code = _mp_real_tan_fp(ys, &w, &err, &recip, d, xn, e, pb, wmax);
                k = 0;
            }
            else
                code = _mp_real_tan_pi_fp(ys, &w, &k, &err, &recip, d, xn, e, pb, wmax);

            if (code & MP_REAL_TRIG_EXACT)
            {
                if (code & MP_REAL_TRIG_ONE)
                {
                    arf_one(rt);
                    if ((code & 1) ^ xsgn)
                        arf_neg(rt, rt);
                }
                else if ((code & 1) ^ (mode == TRIG_COT))
                    arf_nan(rt);    /* a pole */
                else
                    arf_zero(rt);
                TMP_END;
                if (eshift)
                    fmpz_clear(xe);
                return 0;
            }

            if (code & MP_REAL_TRIG_TINY)
            {
                /* the tangent in [Y, Y + 6) units; its reciprocal by one
                   division of the frame */
                if (mode == TRIG_COT)
                    status = _arf_set_mpn_err_quot(rt, ys, tlen, k, 6, NULL, 0, 0, 1, xsgn, nk, prec, rnd);
                else
                    status = _arf_set_mpn_err(rt, ys, tlen, k, -1, 3, xsgn, prec, rnd);
                if (status != ARF_ELEM_UNCERTAIN && eshift)
                    _arf_trig_restore_exp(rt, xe, e, mode == TRIG_COT);
            }
            else if (recip == (mode == TRIG_COT))
                status = _set(rt, ys, w, err, code ^ xsgn, prec, rnd);
            else
                status = _arf_set_mpn_err_quot(rt, ys, w + 1, -FLINT_BITS * w, err, NULL, 0, 0, 1, code ^ xsgn, nk, prec, rnd);
        }

        TMP_END;

        if (status != ARF_ELEM_UNCERTAIN)
        {
            if (eshift)
                fmpz_clear(xe);
            return status;
        }
    }
}

/* the special values, and x copied when an output aliases it */
#define TRIG_SPECIAL(rs, rc, rt, x, sin0, cos0) \
    if (arf_is_special(x)) \
    { \
        if (arf_is_zero(x)) \
        { \
            if ((rs) != NULL) arf_set_si(rs, sin0); \
            if ((rc) != NULL) arf_set_si(rc, cos0); \
            if ((rt) != NULL) arf_zero(rt); \
        } \
        else \
        { \
            if ((rs) != NULL) arf_nan(rs); \
            if ((rc) != NULL) arf_nan(rc); \
            if ((rt) != NULL) arf_nan(rt); \
        } \
        return 0; \
    }

int
arf_sin_cos(arf_ptr res1, arf_ptr res2, arf_srcptr x, slong prec, arf_rnd_t rnd)
{
    TRIG_SPECIAL(res1, res2, NULL, x, 0, 1)

    if (res1 == x || res2 == x)
    {
        arf_t t;
        int status;
        arf_init(t);
        arf_set(t, x);
        status = _arf_trig(res1, res2, NULL, t, 0, 0, prec, rnd);
        arf_clear(t);
        return status;
    }

    return _arf_trig(res1, res2, NULL, x, 0, 0, prec, rnd);
}

int
arf_sin(arf_ptr res, arf_srcptr x, slong prec, arf_rnd_t rnd)
{
    TRIG_SPECIAL(res, NULL, NULL, x, 0, 1)
    return _arf_trig(res, NULL, NULL, x, 0, 0, prec, rnd);
}

int
arf_cos(arf_ptr res, arf_srcptr x, slong prec, arf_rnd_t rnd)
{
    TRIG_SPECIAL(NULL, res, NULL, x, 0, 1)
    return _arf_trig(NULL, res, NULL, x, 0, 0, prec, rnd);
}

int
arf_tan(arf_ptr res, arf_srcptr x, slong prec, arf_rnd_t rnd)
{
    TRIG_SPECIAL(NULL, NULL, res, x, 0, 1)
    return _arf_trig(NULL, NULL, res, x, 0, 0, prec, rnd);
}

int
arf_sin_cos_pi(arf_ptr res1, arf_ptr res2, arf_srcptr x, slong prec, arf_rnd_t rnd)
{
    TRIG_SPECIAL(res1, res2, NULL, x, 0, 1)

    if (res1 == x || res2 == x)
    {
        arf_t t;
        int status;
        arf_init(t);
        arf_set(t, x);
        status = _arf_trig(res1, res2, NULL, t, 1, 0, prec, rnd);
        arf_clear(t);
        return status;
    }

    return _arf_trig(res1, res2, NULL, x, 1, 0, prec, rnd);
}

int
arf_sin_pi(arf_ptr res, arf_srcptr x, slong prec, arf_rnd_t rnd)
{
    TRIG_SPECIAL(res, NULL, NULL, x, 0, 1)
    return _arf_trig(res, NULL, NULL, x, 1, 0, prec, rnd);
}

int
arf_cos_pi(arf_ptr res, arf_srcptr x, slong prec, arf_rnd_t rnd)
{
    TRIG_SPECIAL(NULL, res, NULL, x, 0, 1)
    return _arf_trig(NULL, res, NULL, x, 1, 0, prec, rnd);
}

int
arf_tan_pi(arf_ptr res, arf_srcptr x, slong prec, arf_rnd_t rnd)
{
    TRIG_SPECIAL(NULL, NULL, res, x, 0, 1)
    return _arf_trig(NULL, NULL, res, x, 1, 0, prec, rnd);
}

/* the variants: res may alias x (the special values first) */

int
arf_cot(arf_ptr res, arf_srcptr x, slong prec, arf_rnd_t rnd)
{
    if (arf_is_special(x))
    {
        arf_nan(res);
        return 0;
    }
    if (res == x)
    {
        arf_t t;
        int status;
        arf_init(t);
        arf_set(t, x);
        status = _arf_trig(NULL, NULL, res, t, 0, TRIG_COT, prec, rnd);
        arf_clear(t);
        return status;
    }
    return _arf_trig(NULL, NULL, res, x, 0, TRIG_COT, prec, rnd);
}

int
arf_sec(arf_ptr res, arf_srcptr x, slong prec, arf_rnd_t rnd)
{
    TRIG_SPECIAL(NULL, res, NULL, x, 0, 1)
    if (res == x)
    {
        arf_t t;
        int status;
        arf_init(t);
        arf_set(t, x);
        status = _arf_trig(NULL, res, NULL, t, 0, TRIG_SEC, prec, rnd);
        arf_clear(t);
        return status;
    }
    return _arf_trig(NULL, res, NULL, x, 0, TRIG_SEC, prec, rnd);
}

int
arf_csc(arf_ptr res, arf_srcptr x, slong prec, arf_rnd_t rnd)
{
    if (arf_is_special(x))
    {
        arf_nan(res);
        return 0;
    }
    if (res == x)
    {
        arf_t t;
        int status;
        arf_init(t);
        arf_set(t, x);
        status = _arf_trig(res, NULL, NULL, t, 0, TRIG_CSC, prec, rnd);
        arf_clear(t);
        return status;
    }
    return _arf_trig(res, NULL, NULL, x, 0, TRIG_CSC, prec, rnd);
}

int
arf_sinc(arf_ptr res, arf_srcptr x, slong prec, arf_rnd_t rnd)
{
    TRIG_SPECIAL(NULL, res, NULL, x, 0, 1)
    if (res == x)
    {
        arf_t t;
        int status;
        arf_init(t);
        arf_set(t, x);
        status = _arf_trig(res, NULL, NULL, t, 0, TRIG_SINC, prec, rnd);
        arf_clear(t);
        return status;
    }
    return _arf_trig(res, NULL, NULL, x, 0, TRIG_SINC, prec, rnd);
}

int
arf_cot_pi(arf_ptr res, arf_srcptr x, slong prec, arf_rnd_t rnd)
{
    if (arf_is_special(x))
    {
        arf_nan(res);
        return 0;
    }
    if (res == x)
    {
        arf_t t;
        int status;
        arf_init(t);
        arf_set(t, x);
        status = _arf_trig(NULL, NULL, res, t, 1, TRIG_COT, prec, rnd);
        arf_clear(t);
        return status;
    }
    return _arf_trig(NULL, NULL, res, x, 1, TRIG_COT, prec, rnd);
}

int
arf_sec_pi(arf_ptr res, arf_srcptr x, slong prec, arf_rnd_t rnd)
{
    TRIG_SPECIAL(NULL, res, NULL, x, 0, 1)
    if (res == x)
    {
        arf_t t;
        int status;
        arf_init(t);
        arf_set(t, x);
        status = _arf_trig(NULL, res, NULL, t, 1, TRIG_SEC, prec, rnd);
        arf_clear(t);
        return status;
    }
    return _arf_trig(NULL, res, NULL, x, 1, TRIG_SEC, prec, rnd);
}

int
arf_csc_pi(arf_ptr res, arf_srcptr x, slong prec, arf_rnd_t rnd)
{
    if (arf_is_special(x))
    {
        arf_nan(res);
        return 0;
    }
    if (res == x)
    {
        arf_t t;
        int status;
        arf_init(t);
        arf_set(t, x);
        status = _arf_trig(res, NULL, NULL, t, 1, TRIG_CSC, prec, rnd);
        arf_clear(t);
        return status;
    }
    return _arf_trig(res, NULL, NULL, x, 1, TRIG_CSC, prec, rnd);
}
