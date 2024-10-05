/*
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nfloat.h"

/* TODO: don't use ctx->sizeof_elem, to allow real dot products to be
   called from complex context objects */

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

static void
_nfloat_vec_dot_set_initial(ulong * s, slong sexp, nfloat_srcptr x, int subtract, slong nlimbs)
{
    slong xexp, delta;
    slong delta_limbs, delta_bits;
    int xsgnbit;

    xexp = NFLOAT_EXP(x);
    xsgnbit = NFLOAT_SGNBIT(x) ^ subtract;
    delta = sexp - xexp;

    if (delta < FLINT_BITS)
    {
        FLINT_ASSERT(delta != 0);

        mpn_rshift(s + 1, NFLOAT_D(x), nlimbs, delta);
        s[0] = NFLOAT_D(x)[0] << (FLINT_BITS - delta);
    }
    else
    {
        delta_limbs = delta / FLINT_BITS;
        delta_bits = delta % FLINT_BITS;

        if (delta_limbs < nlimbs + 1)
        {
            if (delta_bits == 0)
                flint_mpn_copyi(s, NFLOAT_D(x) + delta_limbs - 1, nlimbs + 1 - delta_limbs);
            else
                mpn_rshift(s, NFLOAT_D(x) + delta_limbs - 1, nlimbs + 1 - delta_limbs, delta_bits);
        }
    }

    if (xsgnbit)
        mpn_neg(s, s, nlimbs + 1);
}

/* n is the number of limbs of s and x */
static void
_nfloat_vec_dot_add_limbs(ulong * s, slong sexp, mp_srcptr x, slong xexp, int subtract, slong n)
{
    slong delta, delta_limbs, delta_bits;
    ulong t[NFLOAT_MAX_LIMBS + 1];

    if (sexp < xexp)
    {
        delta_bits = xexp - sexp;
        delta_limbs = 0;

        FLINT_ASSERT(delta < FLINT_BITS);

        mpn_lshift(t, x, n, delta_bits);

        if (subtract)
            mpn_sub_n(s, s, t, n);
        else
            mpn_add_n(s, s, t, n);
    }
    else
    {
        delta = sexp - xexp;

        delta_limbs = delta / FLINT_BITS;
        delta_bits = delta % FLINT_BITS;

        if (delta_limbs < n)
        {
            if (delta_bits == 0)
            {
                if (subtract)
                    mpn_sub(s, s, n, x + delta_limbs, n - delta_limbs);
                else
                    mpn_add(s, s, n, x + delta_limbs, n - delta_limbs);
            }
            else
            {
                mpn_rshift(t, x + delta_limbs, n - delta_limbs, delta_bits);

                if (subtract)
                    mpn_sub(s, s, n, t, n - delta_limbs);
                else
                    mpn_add(s, s, n, t, n - delta_limbs);
            }
        }
    }
}

#include <math.h>

static void
_nfloat_vec_dot_addmul(ulong * s, slong sexp, nfloat_srcptr xi, nfloat_srcptr yi, int subtract, slong nlimbs)
{
    slong xexp, delta;
    slong delta_limbs, delta_bits;
    int xsgnbit;
    ulong t[NFLOAT_MAX_LIMBS + 1];

    xexp = NFLOAT_EXP(xi) + NFLOAT_EXP(yi);
    xsgnbit = NFLOAT_SGNBIT(xi) ^ NFLOAT_SGNBIT(yi) ^ subtract;
    delta = sexp - xexp;

    if (delta < FLINT_BITS)
    {
        t[0] = flint_mpn_mulhigh_n(t + 1, NFLOAT_D(xi), NFLOAT_D(yi), nlimbs);

        if (nlimbs == 3)
        {
            if (xsgnbit)
                sub_ddddmmmmssss(s[3], s[2], s[1], s[0], s[3], s[2], s[1], s[0],
                    t[3] >> delta,
                    (t[2] >> delta) | (t[3] << (FLINT_BITS - delta)),
                    (t[1] >> delta) | (t[2] << (FLINT_BITS - delta)),
                    (t[0] >> delta) | (t[1] << (FLINT_BITS - delta)));
            else
                add_ssssaaaaaaaa(s[3], s[2], s[1], s[0], s[3], s[2], s[1], s[0],
                    t[3] >> delta,
                    (t[2] >> delta) | (t[3] << (FLINT_BITS - delta)),
                    (t[1] >> delta) | (t[2] << (FLINT_BITS - delta)),
                    (t[0] >> delta) | (t[1] << (FLINT_BITS - delta)));
        }
        else if (nlimbs == 4)
        {
            if (xsgnbit)
                sub_dddddmmmmmsssss(s[4], s[3], s[2], s[1], s[0], s[4], s[3], s[2], s[1], s[0],
                    t[4] >> delta,
                    (t[3] >> delta) | (t[4] << (FLINT_BITS - delta)),
                    (t[2] >> delta) | (t[3] << (FLINT_BITS - delta)),
                    (t[1] >> delta) | (t[2] << (FLINT_BITS - delta)),
                    (t[0] >> delta) | (t[1] << (FLINT_BITS - delta)));
            else
                add_sssssaaaaaaaaaa(s[4], s[3], s[2], s[1], s[0], s[4], s[3], s[2], s[1], s[0],
                    t[4] >> delta,
                    (t[3] >> delta) | (t[4] << (FLINT_BITS - delta)),
                    (t[2] >> delta) | (t[3] << (FLINT_BITS - delta)),
                    (t[1] >> delta) | (t[2] << (FLINT_BITS - delta)),
                    (t[0] >> delta) | (t[1] << (FLINT_BITS - delta)));
        }
        else
        {
            mpn_rshift(t, t, nlimbs + 1, delta);

            if (xsgnbit)
                mpn_sub_n(s, s, t, nlimbs + 1);
            else
                mpn_add_n(s, s, t, nlimbs + 1);
        }
    }
    else
    {
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

static int
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

#if 0
                FLINT_MPN_MUL_2X2H(t3, t2, t1, NFLOAT_D(xi)[1], NFLOAT_D(xi)[0], NFLOAT_D(yi)[1], NFLOAT_D(yi)[0]);
                (void) t0;
#else
                FLINT_MPN_MUL_2X2(t3, t2, t1, t0, NFLOAT_D(xi)[1], NFLOAT_D(xi)[0], NFLOAT_D(yi)[1], NFLOAT_D(yi)[0]);
                (void) t0;
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
        slong i, xexp, sexp, pad_bits;
        int xsgnbit;
        ulong s[NFLOAT_MAX_LIMBS + 1];
        nfloat_srcptr xi, yi;

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
            _nfloat_vec_dot_set_initial(s, sexp, initial, subtract, nlimbs);

        for (i = 0, xi = x, yi = y; i < len; i++, xi = (char *) xi + sizeof_xstep, yi = (char *) yi + sizeof_ystep)
        {
            if (NFLOAT_IS_ZERO(xi) || NFLOAT_IS_ZERO(yi))
                continue;

            _nfloat_vec_dot_addmul(s, sexp, xi, yi, 0, nlimbs);
        }

        if (LIMB_MSB_IS_SET(s[nlimbs]))
        {
            mpn_neg(s, s, nlimbs + 1);
            xsgnbit = !subtract;
        }
        else
        {
            xsgnbit = subtract;
        }

        return nfloat_set_mpn_2exp(res, s, nlimbs + 1, sexp, xsgnbit, ctx);
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

/*
void
print_s(ulong * s, slong sexp, slong nlimbs)
{
    ulong x = s[nlimbs - 1];
    double c = 1.0;

    if (LIMB_MSB_IS_SET(x))
    {
        x = -x;
        c = -1.0;
    }

    flint_printf("SUM %.15e\n", (double) x * ldexp(1.0, sexp - FLINT_BITS) * c);
}
*/

FLINT_FORCE_INLINE void
_nfloat_vec_dot_addmul_1(ulong * _s1, ulong * _s0, slong sexp, nfloat_srcptr xi, nfloat_srcptr yi, int subtract)
{
    slong xexp, delta;
    ulong s0, s1, t0, t1;
    int xsgnbit;

    xexp = NFLOAT_EXP(xi) + NFLOAT_EXP(yi);
    delta = sexp - xexp;
    xsgnbit = NFLOAT_SGNBIT(xi) ^ NFLOAT_SGNBIT(yi) ^ subtract;

    s0 = *_s0;
    s1 = *_s1;

    if (delta < FLINT_BITS)
    {
        umul_ppmm(t1, t0, NFLOAT_D(xi)[0], NFLOAT_D(yi)[0]);
        if (xsgnbit)
            sub_ddmmss(s1, s0, s1, s0, t1 >> delta, (t0 >> delta) | (t1 << (FLINT_BITS - delta)));
        else
            add_ssaaaa(s1, s0, s1, s0, t1 >> delta, (t0 >> delta) | (t1 << (FLINT_BITS - delta)));
    }
    else if (delta < 2 * FLINT_BITS)
    {
        umul_ppmm(t1, t0, NFLOAT_D(xi)[0], NFLOAT_D(yi)[0]);
        if (xsgnbit)
            sub_ddmmss(s1, s0, s1, s0, 0, t1 >> (delta - FLINT_BITS));
        else
            add_ssaaaa(s1, s0, s1, s0, 0, t1 >> (delta - FLINT_BITS));
    }

    *_s0 = s0;
    *_s1 = s1;
}

FLINT_FORCE_INLINE void
_nfloat_vec_dot_addmul_2(ulong * _s2, ulong * _s1, ulong * _s0, slong sexp, nfloat_srcptr xi, nfloat_srcptr yi, int subtract)
{
    slong xexp, delta;
    ulong s0, s1, s2, t0, t1, t2, t3;
    int xsgnbit;

    s0 = *_s0;
    s1 = *_s1;
    s2 = *_s2;

    xexp = NFLOAT_EXP(xi) + NFLOAT_EXP(yi);
    xsgnbit = NFLOAT_SGNBIT(xi) ^ NFLOAT_SGNBIT(yi) ^ subtract;
    delta = sexp - xexp;
    FLINT_ASSERT(delta != 0);

    if (delta < FLINT_BITS)
    {
        FLINT_MPN_MUL_2X2(t3, t2, t1, t0, NFLOAT_D(xi)[1], NFLOAT_D(xi)[0], NFLOAT_D(yi)[1], NFLOAT_D(yi)[0]);
        (void) t0;

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

    *_s0 = s0;
    *_s1 = s1;
    *_s2 = s2;
}

static int
_flint_mpn_signed_add_n(nn_ptr res, nn_srcptr x, int xsgnbit, nn_srcptr y, int ysgnbit, mp_size_t n)
{
    if (xsgnbit == ysgnbit)
        mpn_add_n(res, x, y, n);
    else
    {
        if (mpn_cmp(x, y, n) >= 0)
            mpn_sub_n(res, x, y, n);
        else
        {
            mpn_sub_n(res, y, x, n);
            xsgnbit = !xsgnbit;
        }
    }

    return xsgnbit;
}

static int
__nfloat_complex_vec_dot(nfloat_complex_ptr res, nfloat_complex_srcptr initial, int subtract, nfloat_complex_srcptr x, slong xstep, nfloat_complex_srcptr y, slong ystep, slong len, gr_ctx_t ctx)
{
    slong i, pad_bits;
    slong re_exp, im_exp;
    ulong re_s[NFLOAT_MAX_LIMBS + 1];
    ulong im_s[NFLOAT_MAX_LIMBS + 1];
    nfloat_srcptr a, b, c, d;
    nfloat_srcptr re_initial, im_initial;
    nfloat_ptr re_res, im_res;
    slong aexp, bexp, cexp, dexp;
    nfloat_srcptr xi, yi;
    int status = GR_SUCCESS;
    int re_sgnbit, im_sgnbit;

    slong nlimbs = NFLOAT_CTX_NLIMBS(ctx);
    slong real_part_size = NFLOAT_CTX_DATA_NLIMBS(ctx);

    re_exp = WORD_MIN;
    im_exp = WORD_MIN;

    re_initial = NFLOAT_COMPLEX_RE(initial, ctx);
    im_initial = NFLOAT_COMPLEX_IM(initial, ctx);

    if (initial != NULL)
    {
        if (!NFLOAT_IS_ZERO(re_initial))
            re_exp = NFLOAT_EXP(re_initial);
        if (!NFLOAT_IS_ZERO(im_initial))
            im_exp = NFLOAT_EXP(im_initial);
    }

    for (i = 0, xi = x, yi = y; i < len; i++, xi = (ulong *) xi + xstep, yi = (ulong *) yi + ystep)
    {
        a = xi;
        b = (ulong *) xi + real_part_size;
        c = yi;
        d = (ulong *) yi + real_part_size;

        aexp = NFLOAT_EXP(a);
        bexp = NFLOAT_EXP(b);
        cexp = NFLOAT_EXP(c);
        dexp = NFLOAT_EXP(d);

        if (aexp != NFLOAT_EXP_ZERO && cexp != NFLOAT_EXP_ZERO)
            re_exp = FLINT_MAX(re_exp, aexp + cexp);
        if (bexp != NFLOAT_EXP_ZERO && dexp != NFLOAT_EXP_ZERO)
            re_exp = FLINT_MAX(re_exp, bexp + dexp);
        if (aexp != NFLOAT_EXP_ZERO && dexp != NFLOAT_EXP_ZERO)
            im_exp = FLINT_MAX(im_exp, aexp + dexp);
        if (bexp != NFLOAT_EXP_ZERO && cexp != NFLOAT_EXP_ZERO)
            im_exp = FLINT_MAX(im_exp, bexp + cexp);
    }

    if (re_exp == WORD_MIN && im_exp == WORD_MIN)
        return nfloat_complex_zero(res, ctx);

    pad_bits = FLINT_BIT_COUNT(len + 1) + 2;
    if (re_exp != WORD_MIN) re_exp += pad_bits;
    if (im_exp != WORD_MIN) im_exp += pad_bits;

    flint_mpn_zero(re_s, nlimbs + 1);
    flint_mpn_zero(im_s, nlimbs + 1);

    if (initial != NULL)
    {
        if (!NFLOAT_IS_ZERO(re_initial))
            _nfloat_vec_dot_set_initial(re_s, re_exp, re_initial, subtract, nlimbs);
        if (!NFLOAT_IS_ZERO(im_initial))
            _nfloat_vec_dot_set_initial(im_s, im_exp, im_initial, subtract, nlimbs);
    }

    if (nlimbs == 1)
    {
        ulong re0, re1;
        ulong im0, im1;

        re0 = re_s[0];
        re1 = re_s[1];
        im0 = im_s[0];
        im1 = im_s[1];

        for (i = 0, xi = x, yi = y; i < len; i++, xi = (ulong *) xi + xstep, yi = (ulong *) yi + ystep)
        {
            a = xi;
            b = (ulong *) xi + real_part_size;
            c = yi;
            d = (ulong *) yi + real_part_size;

            if (!NFLOAT_IS_ZERO(a))
            {
                if (!NFLOAT_IS_ZERO(c))
                    _nfloat_vec_dot_addmul_1(&re1, &re0, re_exp, a, c, 0);
                if (!NFLOAT_IS_ZERO(d))
                    _nfloat_vec_dot_addmul_1(&im1, &im0, im_exp, a, d, 0);
            }

            if (!NFLOAT_IS_ZERO(b))
            {
                if (!NFLOAT_IS_ZERO(d))
                    _nfloat_vec_dot_addmul_1(&re1, &re0, re_exp, b, d, 1);
                if (!NFLOAT_IS_ZERO(c))
                    _nfloat_vec_dot_addmul_1(&im1, &im0, im_exp, b, c, 0);
            }
        }

        re_s[0] = re0;
        re_s[1] = re1;
        im_s[0] = im0;
        im_s[1] = im1;
    }
    else if (nlimbs == 2)
    {
        ulong re0, re1, re2;
        ulong im0, im1, im2;

        re0 = re_s[0];
        re1 = re_s[1];
        re2 = re_s[2];
        im0 = im_s[0];
        im1 = im_s[1];
        im2 = im_s[2];

        for (i = 0, xi = x, yi = y; i < len; i++, xi = (ulong *) xi + xstep, yi = (ulong *) yi + ystep)
        {
            a = xi;
            b = (ulong *) xi + real_part_size;
            c = yi;
            d = (ulong *) yi + real_part_size;

            if (!NFLOAT_IS_ZERO(a))
            {
                if (!NFLOAT_IS_ZERO(c))
                    _nfloat_vec_dot_addmul_2(&re2, &re1, &re0, re_exp, a, c, 0);
                if (!NFLOAT_IS_ZERO(d))
                    _nfloat_vec_dot_addmul_2(&im2, &im1, &im0, im_exp, a, d, 0);
            }

            if (!NFLOAT_IS_ZERO(b))
            {
                if (!NFLOAT_IS_ZERO(d))
                    _nfloat_vec_dot_addmul_2(&re2, &re1, &re0, re_exp, b, d, 1);
                if (!NFLOAT_IS_ZERO(c))
                    _nfloat_vec_dot_addmul_2(&im2, &im1, &im0, im_exp, b, c, 0);
            }
        }

        re_s[0] = re0;
        re_s[1] = re1;
        re_s[2] = re2;

        im_s[0] = im0;
        im_s[1] = im1;
        im_s[2] = im2;
    }
    else
    {
        for (i = 0, xi = x, yi = y; i < len; i++, xi = (ulong *) xi + xstep, yi = (ulong *) yi + ystep)
        {
            a = xi;
            b = (ulong *) xi + real_part_size;
            c = yi;
            d = (ulong *) yi + real_part_size;

            /* Want Karatsuba? */
            if (nlimbs >= 12 && !NFLOAT_IS_ZERO(a) && !NFLOAT_IS_ZERO(b) &&
                    !NFLOAT_IS_ZERO(c) && !NFLOAT_IS_ZERO(d))
            {
                slong abexp, cdexp, abcdexp, adelta, bdelta, cdelta, ddelta;

                aexp = NFLOAT_EXP(a);
                bexp = NFLOAT_EXP(b);
                cexp = NFLOAT_EXP(c);
                dexp = NFLOAT_EXP(d);

                abexp = FLINT_MAX(aexp, bexp) + 2;
                cdexp = FLINT_MAX(cexp, dexp) + 2;
                abcdexp = abexp + cdexp;

                adelta = abexp - aexp;
                bdelta = abexp - bexp;
                cdelta = cdexp - cexp;
                ddelta = cdexp - dexp;

                if (adelta < FLINT_BITS && bdelta < FLINT_BITS &&
                    cdelta < FLINT_BITS && ddelta < FLINT_BITS &&
                    /* To do: check very carefully that the rounded sums
                       cannot overflow the guard bits in re_s and im_s. */
                    abcdexp < re_exp + FLINT_BITS - 8 &&
                    abcdexp < im_exp + FLINT_BITS - 8 &&
                    /* To do: use with truncation also when this term
                       is small. */
                    abcdexp > re_exp - 2 * FLINT_BITS &&
                    abcdexp > im_exp - 2 * FLINT_BITS)
                {
                    ulong aa[NFLOAT_MAX_LIMBS + 1];
                    ulong bb[NFLOAT_MAX_LIMBS + 1];
                    ulong cc[NFLOAT_MAX_LIMBS + 1];
                    ulong dd[NFLOAT_MAX_LIMBS + 1];
                    ulong s[NFLOAT_MAX_LIMBS + 1];
                    ulong t[NFLOAT_MAX_LIMBS + 1];
                    ulong u[NFLOAT_MAX_LIMBS + 1];
                    ulong v[NFLOAT_MAX_LIMBS + 1];

                    int asgnbit, bsgnbit, csgnbit, dsgnbit;
                    int ssgnbit, tsgnbit, usgnbit;

                    asgnbit = NFLOAT_SGNBIT(a);
                    bsgnbit = NFLOAT_SGNBIT(b);
                    csgnbit = NFLOAT_SGNBIT(c);
                    dsgnbit = NFLOAT_SGNBIT(d);

                    /*
                        s = c * (a + b)
                        t = a * (d - c)
                        u = b * (c + d)
                        re = s - u
                        im = s + t
                    */

                    aa[0] = mpn_rshift(aa + 1, NFLOAT_D(a), nlimbs, adelta);
                    bb[0] = mpn_rshift(bb + 1, NFLOAT_D(b), nlimbs, bdelta);
                    cc[0] = mpn_rshift(cc + 1, NFLOAT_D(c), nlimbs, cdelta);
                    dd[0] = mpn_rshift(dd + 1, NFLOAT_D(d), nlimbs, ddelta);

                    ssgnbit = csgnbit ^ _flint_mpn_signed_add_n(v, aa, asgnbit, bb, bsgnbit, nlimbs + 1);
                    flint_mpn_mulhigh_n(s, cc, v, nlimbs + 1);

                    tsgnbit = asgnbit ^ _flint_mpn_signed_add_n(v, dd, dsgnbit, cc, !csgnbit, nlimbs + 1);
                    flint_mpn_mulhigh_n(t, aa, v, nlimbs + 1);

                    usgnbit = bsgnbit ^ _flint_mpn_signed_add_n(v, cc, csgnbit, dd, dsgnbit, nlimbs + 1);
                    flint_mpn_mulhigh_n(u, bb, v, nlimbs + 1);

                    usgnbit = _flint_mpn_signed_add_n(u, s, ssgnbit, u, !usgnbit, nlimbs + 1);
                    tsgnbit = _flint_mpn_signed_add_n(t, s, ssgnbit, t, tsgnbit, nlimbs + 1);

                    /*
                    flint_printf("|u| = ");
                    print_s(u, abexp + cdexp, nlimbs + 1);
                    flint_printf("|t| = ");
                    print_s(t, abexp + cdexp, nlimbs + 1);
                    */

                    _nfloat_vec_dot_add_limbs(re_s, re_exp, u, abcdexp, usgnbit, nlimbs + 1);
                    _nfloat_vec_dot_add_limbs(im_s, im_exp, t, abcdexp, tsgnbit, nlimbs + 1);

                    continue;
                }
            }

            if (!NFLOAT_IS_ZERO(a))
            {
                if (!NFLOAT_IS_ZERO(c))
                    _nfloat_vec_dot_addmul(re_s, re_exp, a, c, 0, nlimbs);
                if (!NFLOAT_IS_ZERO(d))
                    _nfloat_vec_dot_addmul(im_s, im_exp, a, d, 0, nlimbs);
            }

            if (!NFLOAT_IS_ZERO(b))
            {
                if (!NFLOAT_IS_ZERO(d))
                    _nfloat_vec_dot_addmul(re_s, re_exp, b, d, 1, nlimbs);
                if (!NFLOAT_IS_ZERO(c))
                    _nfloat_vec_dot_addmul(im_s, im_exp, b, c, 0, nlimbs);
            }
        }
    }

    /* todo: maybe specialize for nlimbs=1, 2 */

    re_res = res;
    im_res = (nn_ptr) res + real_part_size;

    if (re_exp == WORD_MIN)
    {
        status |= nfloat_zero(re_res, ctx);
    }
    else
    {
        re_sgnbit = LIMB_MSB_IS_SET(re_s[nlimbs]);
        if (re_sgnbit)
            mpn_neg(re_s, re_s, nlimbs + 1);
        status |= nfloat_set_mpn_2exp(re_res, re_s, nlimbs + 1, re_exp, re_sgnbit ^ subtract, ctx);
    }

    if (im_exp == WORD_MIN)
    {
        status |= nfloat_zero(im_res, ctx);
    }
    else
    {
        im_sgnbit = LIMB_MSB_IS_SET(im_s[nlimbs]);
        if (im_sgnbit)
            mpn_neg(im_s, im_s, nlimbs + 1);
        status |= nfloat_set_mpn_2exp(im_res, im_s, nlimbs + 1, im_exp, im_sgnbit ^ subtract, ctx);
    }

    return status;
}

int
_nfloat_complex_vec_dot(nfloat_complex_ptr res, nfloat_complex_srcptr initial, int subtract, nfloat_complex_srcptr x, nfloat_complex_srcptr y, slong len, gr_ctx_t ctx)
{
    return __nfloat_complex_vec_dot(res, initial, subtract, x, NFLOAT_COMPLEX_CTX_DATA_NLIMBS(ctx), y, NFLOAT_COMPLEX_CTX_DATA_NLIMBS(ctx), len, ctx);
}

int
_nfloat_complex_vec_dot_rev(nfloat_complex_ptr res, nfloat_complex_srcptr initial, int subtract, nfloat_complex_srcptr x, nfloat_complex_srcptr y, slong len, gr_ctx_t ctx)
{
    return __nfloat_complex_vec_dot(res, initial, subtract, x, NFLOAT_COMPLEX_CTX_DATA_NLIMBS(ctx), GR_ENTRY(y, len - 1, ctx->sizeof_elem), -NFLOAT_COMPLEX_CTX_DATA_NLIMBS(ctx), len, ctx);
}
