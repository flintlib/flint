/*
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nfloat.h"

/* Having some dummy code in the continue block makes the loop go faster. Don't ask me why. */
#if 1
#define DUMMY_OP add_ssaaaa(s1, s0, s1, s0, 0, 0)
#else
#define DUMMY_OP 0
#endif

/* Experimental: use fast but inaccurate product? */
#define FLINT_MPN_MUL_2X2H(r3, r2, r1, a1, a0, b1, b0)   \
    do {                                                                  \
        ulong __t1, __t2, __u1, __u2, __v3, __v2;                     \
        ulong __r3, __r2, __r1;                                       \
        ulong __a1 = (a1), __a0 = (a0), __b1 = (b1), __b0 = (b0);     \
        umul_ppmm(__t2, __t1, __a0, __b1);                                \
        umul_ppmm(__u2, __u1, __a1, __b0);                                \
        add_sssaaaaaa(__r3, __r2, __r1, 0, __t2, __t1, 0, __u2, __u1);    \
        umul_ppmm(__v3, __v2, __a1, __b1);                                \
        add_ssaaaa(__r3, __r2, __r3, __r2, __v3, __v2);                   \
        (r1) = __r1; (r2) = __r2; (r3) = __r3;                            \
    } while (0)

int
__nfloat_vec_dot(nfloat_ptr res, nfloat_srcptr initial, int subtract, nfloat_srcptr x, slong sizeof_xstep, nfloat_srcptr y, slong sizeof_ystep, slong len, gr_ctx_t ctx)
{
    if (NFLOAT_CTX_NLIMBS(ctx) == 1 && !(NFLOAT_CTX_FLAGS(ctx) & (NFLOAT_ALLOW_INF | NFLOAT_ALLOW_NAN)))
    {
        slong i, xexp, delta, sexp, norm, pad_bits;
        int xsgnbit;
        ulong t1, t0, s1, s0;
        nfloat_srcptr xi, yi;
        int have_zeros = 0;

        s1 = s0 = 0;
        sexp = WORD_MIN;

        if (initial != NULL && !NFLOAT_IS_ZERO(initial))
            sexp = NFLOAT_EXP(initial);

        /* Lookahead to determine a working exponent. This is about 25% slower
           than updating the exponent in a single loop, but easier, so
           let's go with it for now. */
        for (i = 0, xi = x, yi = y; i < len; i++, xi = (char *) xi + sizeof_xstep, yi = (char *) yi + sizeof_ystep)
        {
            if (NFLOAT_IS_ZERO(xi) || NFLOAT_IS_ZERO(yi))
            {
                /* DUMMY_OP; */
                have_zeros = 1;
            }
            else
            {
                xexp = NFLOAT_EXP(xi) + NFLOAT_EXP(yi);
                sexp = FLINT_MAX(sexp, xexp);
            }
        }

        if (sexp == WORD_MIN)
            return nfloat_zero(res, ctx);

        pad_bits = FLINT_BIT_COUNT(len + 1) + 1;
        sexp += pad_bits;

        if (initial != NULL && !NFLOAT_IS_ZERO(initial))
        {
            xexp = NFLOAT_EXP(initial);
            xsgnbit = NFLOAT_SGNBIT(initial) ^ subtract;
            delta = sexp - xexp;

            if (delta < FLINT_BITS)
            {
                t1 = NFLOAT_D(initial)[0];
                s1 = t1 >> delta;
                s0 = t1 << (FLINT_BITS - delta);
                if (xsgnbit)
                    sub_ddmmss(s1, s0, 0, 0, s1, s0);
            }
            else if (delta < 2 * FLINT_BITS)
            {
                t1 = NFLOAT_D(initial)[0];
                s0 = t1 >> (delta - FLINT_BITS);
                if (xsgnbit)
                    sub_ddmmss(s1, s0, 0, 0, s1, s0);
            }
        }

        if (have_zeros)
        {
            for (i = 0, xi = x, yi = y; i < len; i++, xi = (char *) xi + sizeof_xstep, yi = (char *) yi + sizeof_ystep)
            {
                if (NFLOAT_IS_ZERO(xi) || NFLOAT_IS_ZERO(yi))
                {
                    DUMMY_OP;
                    continue;
                }

                xexp = NFLOAT_EXP(xi) + NFLOAT_EXP(yi);
                delta = sexp - xexp;
                FLINT_ASSERT(delta != 0);

                if (delta < FLINT_BITS)
                {
                    xsgnbit = NFLOAT_SGNBIT(xi) ^ NFLOAT_SGNBIT(yi);
                    umul_ppmm(t1, t0, NFLOAT_D(xi)[0], NFLOAT_D(yi)[0]);
                    if (xsgnbit)
                        sub_ddmmss(s1, s0, s1, s0, t1 >> delta, (t0 >> delta) | (t1 << (FLINT_BITS - delta)));
                    else
                        add_ssaaaa(s1, s0, s1, s0, t1 >> delta, (t0 >> delta) | (t1 << (FLINT_BITS - delta)));
                }
                else if (delta < 2 * FLINT_BITS)
                {
                    xsgnbit = NFLOAT_SGNBIT(xi) ^ NFLOAT_SGNBIT(yi);
                    umul_ppmm(t1, t0, NFLOAT_D(xi)[0], NFLOAT_D(yi)[0]);
                    if (xsgnbit)
                        sub_ddmmss(s1, s0, s1, s0, 0, t1 >> (delta - FLINT_BITS));
                    else
                        add_ssaaaa(s1, s0, s1, s0, 0, t1 >> (delta - FLINT_BITS));
                }
            }
        }
        else
        {
            for (i = 0, xi = x, yi = y; i < len; i++, xi = (char *) xi + sizeof_xstep, yi = (char *) yi + sizeof_ystep)
            {
                xexp = NFLOAT_EXP(xi) + NFLOAT_EXP(yi);
                delta = sexp - xexp;

                if (delta < FLINT_BITS)
                {
                    xsgnbit = NFLOAT_SGNBIT(xi) ^ NFLOAT_SGNBIT(yi);
                    umul_ppmm(t1, t0, NFLOAT_D(xi)[0], NFLOAT_D(yi)[0]);
                    if (xsgnbit)
                        sub_ddmmss(s1, s0, s1, s0, t1 >> delta, (t0 >> delta) | (t1 << (FLINT_BITS - delta)));
                    else
                        add_ssaaaa(s1, s0, s1, s0, t1 >> delta, (t0 >> delta) | (t1 << (FLINT_BITS - delta)));
                }
                else if (delta < 2 * FLINT_BITS)
                {
                    xsgnbit = NFLOAT_SGNBIT(xi) ^ NFLOAT_SGNBIT(yi);
                    umul_ppmm(t1, t0, NFLOAT_D(xi)[0], NFLOAT_D(yi)[0]);
                    if (xsgnbit)
                        sub_ddmmss(s1, s0, s1, s0, 0, t1 >> (delta - FLINT_BITS));
                    else
                        add_ssaaaa(s1, s0, s1, s0, 0, t1 >> (delta - FLINT_BITS));
                }
            }
        }

        if (LIMB_MSB_IS_SET(s1))
        {
            sub_ddmmss(s1, s0, 0, 0, s1, s0);
            xsgnbit = 1;
        }
        else
        {
            xsgnbit = 0;
        }

        NFLOAT_SGNBIT(res) = xsgnbit ^ subtract;

        if (s1 != 0)
        {
            norm = flint_clz(s1);
            if (norm)
                NFLOAT_D(res)[0] = (s1 << norm) | (s0 >> (FLINT_BITS - norm));
            else
                NFLOAT_D(res)[0] = s1;
            NFLOAT_EXP(res) = sexp - norm;
        }
        else if (s0 != 0)
        {
            norm = flint_clz(s0);
            NFLOAT_D(res)[0] = s0 << norm;
            NFLOAT_EXP(res) = sexp - FLINT_BITS - norm;
        }
        else
        {
            return nfloat_zero(res, ctx);
        }

        NFLOAT_HANDLE_UNDERFLOW_OVERFLOW(res, ctx);
        return GR_SUCCESS;
    }

    if (NFLOAT_CTX_NLIMBS(ctx) == 2 && !(NFLOAT_CTX_FLAGS(ctx) & (NFLOAT_ALLOW_INF | NFLOAT_ALLOW_NAN)))
    {
        slong i, xexp, delta, sexp, norm, pad_bits;
        int xsgnbit;
        nfloat_srcptr xi, yi;
        ulong s0, s1, s2;
        ulong t0, t1, t2, t3;

        s0 = s1 = s2 = 0;
        sexp = WORD_MIN;

        if (initial != NULL && !NFLOAT_IS_ZERO(initial))
            sexp = NFLOAT_EXP(initial);

        for (i = 0, xi = x, yi = y; i < len; i++, xi = (char *) xi + sizeof_xstep, yi = (char *) yi + sizeof_ystep)
        {
            if (NFLOAT_IS_ZERO(xi) || NFLOAT_IS_ZERO(yi))
            {
                DUMMY_OP;
                continue;
            }

            xexp = NFLOAT_EXP(xi) + NFLOAT_EXP(yi);
            sexp = FLINT_MAX(sexp, xexp);
        }

        if (sexp == WORD_MIN)
            return nfloat_zero(res, ctx);

        pad_bits = FLINT_BIT_COUNT(len + 1) + 1;
        sexp += pad_bits;

        if (initial != NULL && !NFLOAT_IS_ZERO(initial))
        {
            xexp = NFLOAT_EXP(initial);
            xsgnbit = NFLOAT_SGNBIT(initial) ^ subtract;
            delta = sexp - xexp;
            FLINT_ASSERT(delta != 0);

            if (delta < 3 * FLINT_BITS)
            {
                if (delta < FLINT_BITS)
                {
                    s0 = 0;
                    s1 = NFLOAT_D(initial)[0];
                    s2 = NFLOAT_D(initial)[1];
                }
                else if (delta < 2 * FLINT_BITS)
                {
                    s0 = NFLOAT_D(initial)[0];
                    s1 = NFLOAT_D(initial)[1];
                    s2 = 0;

                    delta -= FLINT_BITS;
                }
                else
                {
                    s0 = NFLOAT_D(initial)[1];
                    s1 = 0;
                    s2 = 0;

                    delta -= 2 * FLINT_BITS;
                }

                if (delta != 0)
                {
                    s0 = (s0 >> delta) | s1 << (FLINT_BITS - delta);
                    s1 = (s1 >> delta) | (s2 << (FLINT_BITS - delta));
                    s2 = s2 >> delta;
                }

                if (xsgnbit)
                    sub_dddmmmsss(s2, s1, s0, 0, 0, 0, s2, s1, s0);
            }
        }

        for (i = 0, xi = x, yi = y; i < len; i++, xi = (char *) xi + sizeof_xstep, yi = (char *) yi + sizeof_ystep)
        {
            if (NFLOAT_IS_ZERO(xi) || NFLOAT_IS_ZERO(yi))
                continue;

            xexp = NFLOAT_EXP(xi) + NFLOAT_EXP(yi);
            xsgnbit = NFLOAT_SGNBIT(xi) ^ NFLOAT_SGNBIT(yi);
            delta = sexp - xexp;
            FLINT_ASSERT(delta != 0);

            if (delta < FLINT_BITS)
            {
                xsgnbit = NFLOAT_SGNBIT(xi) ^ NFLOAT_SGNBIT(yi);

#if 1
                FLINT_MPN_MUL_2X2H(t3, t2, t1, NFLOAT_D(xi)[1], NFLOAT_D(xi)[0], NFLOAT_D(yi)[1], NFLOAT_D(yi)[0]);
                (void) t0;
#else
                FLINT_MPN_MUL_2X2(t3, t2, t1, t0, NFLOAT_D(xi)[1], NFLOAT_D(xi)[0], NFLOAT_D(yi)[1], NFLOAT_D(yi)[0]);
#endif

                if (xsgnbit)
                    sub_dddmmmsss(s2, s1, s0, s2, s1, s0, t3 >> delta, (t2 >> delta) | (t3 << (FLINT_BITS - delta)), (t1 >> delta) | (t2 << (FLINT_BITS - delta)));
                else
                    add_sssaaaaaa(s2, s1, s0, s2, s1, s0, t3 >> delta, (t2 >> delta) | (t3 << (FLINT_BITS - delta)), (t1 >> delta) | (t2 << (FLINT_BITS - delta)));
            }
            else if (delta < 3 * FLINT_BITS)
            {
                if (delta < 2 * FLINT_BITS)
                {
                    FLINT_MPN_MUL_2X2H(t3, t2, t1, NFLOAT_D(xi)[1], NFLOAT_D(xi)[0], NFLOAT_D(yi)[1], NFLOAT_D(yi)[0]);

                    delta -= FLINT_BITS;

                    if (delta != 0)
                    {
                        t2 = (t2 >> delta) | (t3 << (FLINT_BITS - delta));
                        t3 = t3 >> delta;
                    }
                }
                else
                {
                    umul_ppmm(t3, t2, NFLOAT_D(xi)[1], NFLOAT_D(yi)[1]);

                    delta -= 2 * FLINT_BITS;

                    if (delta != 0)
                        t3 = t3 >> delta;

                    t2 = t3;
                    t3 = 0;
                }

                if (xsgnbit)
                    sub_dddmmmsss(s2, s1, s0, s2, s1, s0, 0, t3, t2);
                else
                    add_sssaaaaaa(s2, s1, s0, s2, s1, s0, 0, t3, t2);
            }
        }

        if (LIMB_MSB_IS_SET(s2))
        {
            sub_dddmmmsss(s2, s1, s0, 0, 0, 0, s2, s1, s0);
            xsgnbit = 1;
        }
        else
        {
            xsgnbit = 0;
        }

        NFLOAT_SGNBIT(res) = xsgnbit ^ subtract;

        if (s2 != 0)
        {
            norm = flint_clz(s2);
            if (norm)
            {
                NFLOAT_D(res)[0] = (s1 << norm) | (s0 >> (FLINT_BITS - norm));
                NFLOAT_D(res)[1] = (s2 << norm) | (s1 >> (FLINT_BITS - norm));
            }
            else
            {
                NFLOAT_D(res)[0] = s1;
                NFLOAT_D(res)[1] = s2;
            }
            NFLOAT_EXP(res) = sexp - norm;
        }
        else if (s1 != 0)
        {
            norm = flint_clz(s1);
            if (norm)
            {
                NFLOAT_D(res)[0] = (s0 << norm);
                NFLOAT_D(res)[1] = (s1 << norm) | (s0 >> (FLINT_BITS - norm));
            }
            else
            {
                NFLOAT_D(res)[0] = s0;
                NFLOAT_D(res)[1] = s1;
            }
            NFLOAT_EXP(res) = sexp - FLINT_BITS - norm;
        }
        else if (s0 != 0)
        {
            norm = flint_clz(s0);
            NFLOAT_D(res)[0] = 0;
            NFLOAT_D(res)[1] = s0 << norm;
            NFLOAT_EXP(res) = sexp - 2 * FLINT_BITS - norm;
        }
        else
        {
            return nfloat_zero(res, ctx);
        }

        NFLOAT_HANDLE_UNDERFLOW_OVERFLOW(res, ctx);
        return GR_SUCCESS;
    }

    if (!(NFLOAT_CTX_FLAGS(ctx) & (NFLOAT_ALLOW_INF | NFLOAT_ALLOW_NAN)))
    {
        slong i, xexp, delta, sexp, norm, pad_bits;
        int xsgnbit;
        ulong t[NFLOAT_MAX_LIMBS + 1];
        ulong s[NFLOAT_MAX_LIMBS + 1];
        nfloat_srcptr xi, yi;
        slong n;

        slong nlimbs = NFLOAT_CTX_NLIMBS(ctx);

        sexp = WORD_MIN;

        if (initial != NULL && !NFLOAT_IS_ZERO(initial))
            sexp = NFLOAT_EXP(initial);

        for (i = 0, xi = x, yi = y; i < len; i++, xi = (char *) xi + sizeof_xstep, yi = (char *) yi + sizeof_ystep)
        {
            if (NFLOAT_IS_ZERO(xi) || NFLOAT_IS_ZERO(yi))
                continue;

            xexp = NFLOAT_EXP(xi) + NFLOAT_EXP(yi);
            sexp = FLINT_MAX(sexp, xexp);
        }

        if (sexp == WORD_MIN)
            return nfloat_zero(res, ctx);

        flint_mpn_zero(s, nlimbs + 1);

        pad_bits = FLINT_BIT_COUNT(len + 1) + 1;
        sexp += pad_bits;

        if (initial != NULL && !NFLOAT_IS_ZERO(initial))
        {
            xexp = NFLOAT_EXP(initial);
            xsgnbit = NFLOAT_SGNBIT(initial) ^ subtract;
            delta = sexp - xexp;

            if (delta < FLINT_BITS)
            {
                FLINT_ASSERT(delta != 0);

                mpn_rshift(s + 1, NFLOAT_D(initial), nlimbs, delta);
                s[0] = NFLOAT_D(initial)[0] << (FLINT_BITS - delta);
                if (xsgnbit)
                    mpn_neg(s, s, nlimbs + 1);
            }
            else
            {
                slong delta_limbs, delta_bits;

                delta_limbs = delta / FLINT_BITS;
                delta_bits = delta % FLINT_BITS;

                if (delta_limbs < nlimbs + 1)
                {
                    if (delta_bits == 0)
                        flint_mpn_copyi(s, NFLOAT_D(initial) + delta_limbs - 1, nlimbs + 1 - delta_limbs);
                    else
                        mpn_rshift(s, NFLOAT_D(initial) + delta_limbs - 1, nlimbs + 1 - delta_limbs, delta_bits);

                    if (xsgnbit)
                        mpn_neg(s, s, nlimbs + 1);
                }
            }
        }

        for (i = 0, xi = x, yi = y; i < len; i++, xi = (char *) xi + sizeof_xstep, yi = (char *) yi + sizeof_ystep)
        {
            if (NFLOAT_IS_ZERO(xi) || NFLOAT_IS_ZERO(yi))
                continue;

            xexp = NFLOAT_EXP(xi) + NFLOAT_EXP(yi);
            xsgnbit = NFLOAT_SGNBIT(xi) ^ NFLOAT_SGNBIT(yi);
            delta = sexp - xexp;

            if (delta < FLINT_BITS)
            {
                t[0] = flint_mpn_mulhigh_n(t + 1, NFLOAT_D(xi), NFLOAT_D(yi), nlimbs);
                mpn_rshift(t, t, nlimbs + 1, delta);

                if (xsgnbit)
                    mpn_sub_n(s, s, t, nlimbs + 1);
                else
                    mpn_add_n(s, s, t, nlimbs + 1);
            }
            else
            {
                slong delta_limbs, delta_bits;

                delta_limbs = delta / FLINT_BITS;
                delta_bits = delta % FLINT_BITS;

                /* alternative criterion: if (delta < (FLINT_BITS * nlimbs) + 2 * pad_bits) */
                if (delta_limbs < nlimbs + 1)
                {
                    /* todo: squaring case */
                    flint_mpn_mulhigh_n(t, NFLOAT_D(xi) + delta_limbs - 1, NFLOAT_D(yi) + delta_limbs - 1, nlimbs + 1 - delta_limbs);

                    if (delta_bits != 0)
                        mpn_rshift(t, t, nlimbs + 1 - delta_limbs, delta_bits);

                    if (xsgnbit)
                        mpn_sub(s, s, nlimbs + 1, t, nlimbs + 1 - delta_limbs);
                    else
                        mpn_add(s, s, nlimbs + 1, t, nlimbs + 1 - delta_limbs);
                }
                else
                {
                    /* Skip term. */
                }
            }
        }

        if (LIMB_MSB_IS_SET(s[nlimbs]))
        {
            mpn_neg(s, s, nlimbs + 1);
            xsgnbit = 1;
        }
        else
        {
            xsgnbit = 0;
        }

        NFLOAT_SGNBIT(res) = xsgnbit ^ subtract;

        n = nlimbs + 1;
        MPN_NORM(s, n);

        xexp = sexp - (nlimbs + 1 - n) * FLINT_BITS;

        if (n == nlimbs + 1)
        {
            norm = flint_clz(s[n - 1]);
            if (norm)
            {
                mpn_lshift(NFLOAT_D(res), s + 1, nlimbs, norm);
                NFLOAT_D(res)[0] |= (s[0] >> (FLINT_BITS - norm));
            }
            else
                flint_mpn_copyi(NFLOAT_D(res), s + 1, nlimbs);
            xexp -= norm;
        }
        else
        {
            if (n == 0)
                return nfloat_zero(res, ctx);

            norm = flint_clz(s[n - 1]);
            if (norm)
                mpn_lshift(NFLOAT_D(res) + nlimbs - n, s, n, norm);
            else
                flint_mpn_copyi(NFLOAT_D(res) + nlimbs - n, s, n);
            flint_mpn_zero(NFLOAT_D(res), nlimbs - n);
            xexp -= norm;
        }

        NFLOAT_EXP(res) = xexp;
        NFLOAT_HANDLE_UNDERFLOW_OVERFLOW(res, ctx);
        return GR_SUCCESS;
    }

    return GR_UNABLE;
}

int
_nfloat_vec_dot(nfloat_ptr res, nfloat_srcptr initial, int subtract, nfloat_srcptr x, nfloat_srcptr y, slong len, gr_ctx_t ctx)
{
    return __nfloat_vec_dot(res, initial, subtract, x, ctx->sizeof_elem, y, ctx->sizeof_elem, len, ctx);
}

int
_nfloat_vec_dot_rev(nfloat_ptr res, nfloat_srcptr initial, int subtract, nfloat_srcptr x, nfloat_srcptr y, slong len, gr_ctx_t ctx)
{
    return __nfloat_vec_dot(res, initial, subtract, x, ctx->sizeof_elem, GR_ENTRY(y, len - 1, ctx->sizeof_elem), -ctx->sizeof_elem, len, ctx);
}
