/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "mpn_extras.h"

/* We need uint64_t instead of mp_limb_t on 32-bit systems for
   safe summation of 30-bit error bounds. */
#include <stdint.h>

void
_arb_dot_addmul_generic(mp_ptr sum, mp_ptr serr, mp_ptr tmp, mp_size_t sn,
    mp_srcptr xptr, mp_size_t xn, mp_srcptr yptr, mp_size_t yn,
    int negative, flint_bitcnt_t shift);

void
_arb_dot_add_generic(mp_ptr sum, mp_ptr serr, mp_ptr tmp, mp_size_t sn,
    mp_srcptr xptr, mp_size_t xn,
    int negative, flint_bitcnt_t shift);

void
arb_approx_dot_simple(arb_t res, const arb_t initial, int subtract,
    arb_srcptr x, slong xstep, arb_srcptr y, slong ystep, slong len, slong prec)
{
    slong i;

    if (len <= 0)
    {
        if (initial == NULL)
            arf_zero(arb_midref(res));
        else
            arf_set_round(arb_midref(res), arb_midref(initial), prec, ARB_RND);
        return;
    }

    if (initial == NULL)
    {
        arf_mul(arb_midref(res), arb_midref(x), arb_midref(y), prec, ARB_RND);
    }
    else
    {
        if (subtract)
            arf_neg(arb_midref(res), arb_midref(initial));
        else
            arf_set(arb_midref(res), arb_midref(initial));
        arf_addmul(arb_midref(res), arb_midref(x), arb_midref(y), prec, ARB_RND);
    }

    for (i = 1; i < len; i++)
        arf_addmul(arb_midref(res), arb_midref(x + i * xstep), arb_midref(y + i * ystep), prec, ARB_RND);

    if (subtract)
        arf_neg(arb_midref(res), arb_midref(res));
}

void
arb_approx_dot(arb_t res, const arb_t initial, int subtract, arb_srcptr x, slong xstep, arb_srcptr y, slong ystep, slong len, slong prec)
{
    slong i, j, nonzero, padding, extend;
    slong xexp, yexp, exp, max_exp, min_exp, sum_exp;
    int xnegative, ynegative;
    mp_size_t xn, yn, sn, alloc;
    flint_bitcnt_t shift;
    arb_srcptr xi, yi;
    arf_srcptr xm, ym;
    mp_limb_t serr;   /* Sum over arithmetic errors  - not used, but need dummy for calls */
    mp_ptr tmp, sum;  /* Workspace */
    ARF_ADD_TMP_DECL;

    /* todo: fast fma and fmma (len=2) code */
    if (len <= 1)
    {
        if (initial == NULL)
        {
            if (len <= 0)
                arf_zero(arb_midref(res));
            else
            {
                if (subtract)
                    arf_neg_mul(arb_midref(res), arb_midref(x), arb_midref(y), prec, ARB_RND);
                else
                    arf_mul(arb_midref(res), arb_midref(x), arb_midref(y), prec, ARB_RND);
            }
            return;
        }
        else if (len <= 0)
        {
            arf_set_round(arb_midref(res), arb_midref(initial), prec, ARB_RND);
            return;
        }
    }

    /* Number of nonzero midpoint terms in sum. */
    nonzero = 0;

    /* Terms are bounded by 2^max_exp (with WORD_MIN = -infty) */
    max_exp = WORD_MIN;

    /* Used to reduce the precision. */
    min_exp = WORD_MAX;

    /* Account for the initial term. */
    if (initial != NULL)
    {
        if (!ARF_IS_LAGOM(arb_midref(initial)))
        {
            arb_approx_dot_simple(res, initial, subtract, x, xstep, y, ystep, len, prec);
            return;
        }

        xm = arb_midref(initial);

        if (!arf_is_special(xm))
        {
            max_exp = ARF_EXP(xm);
            nonzero++;

            if (prec > 2 * FLINT_BITS)
                min_exp = ARF_EXP(xm) - ARF_SIZE(xm) * FLINT_BITS;
        }
    }

    /* Determine maximum exponents for the main sum and the radius sum. */
    for (i = 0; i < len; i++)
    {
        xi = x + i * xstep;
        yi = y + i * ystep;

        /* Fallback for huge exponents or non-finite values. */
        if (!ARF_IS_LAGOM(arb_midref(xi)) || !ARF_IS_LAGOM(arb_midref(yi)))
        {
            arb_approx_dot_simple(res, initial, subtract, x, xstep, y, ystep, len, prec);
            return;
        }

        xm = arb_midref(xi);
        ym = arb_midref(yi);

        if (!arf_is_special(xm))
        {
            xexp = ARF_EXP(xm);

            if (!arf_is_special(ym))
            {
                yexp = ARF_EXP(ym);

                max_exp = FLINT_MAX(max_exp, xexp + yexp);
                nonzero++;

                if (prec > 2 * FLINT_BITS)
                {
                    slong bot;
                    bot = (xexp + yexp) - (ARF_SIZE(xm) + ARF_SIZE(ym)) * FLINT_BITS;
                    min_exp = FLINT_MIN(min_exp, bot);
                }
            }
        }
    }

    /* The midpoint sum is zero. */
    if (max_exp == WORD_MIN)
    {
        arf_zero(arb_midref(res));
        return;
    }
    else
    {
        /* Reduce precision based on actual sizes. */
        if (min_exp != WORD_MAX)
            prec = FLINT_MIN(prec, max_exp - min_exp + MAG_BITS);

        prec = FLINT_MAX(prec, 2);
    }

    /* Extend sum so that we can use two's complement addition. */
    extend = FLINT_BIT_COUNT(nonzero) + 1;

    /* Extra bits to improve accuracy (optional). */
    padding = 4 + FLINT_BIT_COUNT(len);

    /* Number of limbs. */
    sn = (prec + extend + padding + FLINT_BITS - 1) / FLINT_BITS;

    /* Avoid having to make a special case for sn = 1. */
    sn = FLINT_MAX(sn, 2);

    /* Exponent for the main sum. */
    sum_exp = max_exp + extend;

    /* We need sn + 1 limb for the sum (sn limbs + 1 dummy limb
       for carry or borrow that avoids an extra branch). We need
       2 * (sn + 2) limbs to store the product of two numbers
       with up to (sn + 2) limbs, plus 1 extra limb for shifting
       the product. */
    alloc = (sn + 1) + 2 * (sn + 2) + 1;
    ARF_ADD_TMP_ALLOC(sum, alloc)
    tmp = sum + (sn + 1);

    /* Set sum to 0 */
    serr = 0;
    for (j = 0; j < sn + 1; j++)
        sum[j] = 0;

    if (initial != NULL)
    {
        xm = arb_midref(initial);

        if (!arf_is_special(xm))
        {
            mp_srcptr xptr;

            xexp = ARF_EXP(xm);
            xn = ARF_SIZE(xm);
            xnegative = ARF_SGNBIT(xm);

            shift = sum_exp - xexp;

            if (shift < sn * FLINT_BITS)
            {
                xptr = (xn <= ARF_NOPTR_LIMBS) ? ARF_NOPTR_D(xm) : ARF_PTR_D(xm);
                _arb_dot_add_generic(sum, &serr, tmp, sn, xptr, xn, xnegative ^ subtract, shift);
            }
        }
    }

    for (i = 0; i < len; i++)
    {
        xi = x + i * xstep;
        yi = y + i * ystep;
        xm = arb_midref(xi);
        ym = arb_midref(yi);

        /* The midpoints of x[i] and y[i] are both nonzero. */
        if (!arf_is_special(xm) && !arf_is_special(ym))
        {
            xexp = ARF_EXP(xm);
            xn = ARF_SIZE(xm);
            xnegative = ARF_SGNBIT(xm);

            yexp = ARF_EXP(ym);
            yn = ARF_SIZE(ym);
            ynegative = ARF_SGNBIT(ym);

            exp = xexp + yexp;
            shift = sum_exp - exp;

            if (shift >= sn * FLINT_BITS)
            {
                /* do nothing */
            }
#if 0
            else if (xn == 1 && yn == 1 && sn == 2 && shift < FLINT_BITS)  /* Fastest path. */
            {
                mp_limb_t hi, lo, x0, y0;

                x0 = ARF_NOPTR_D(xm)[0];
                y0 = ARF_NOPTR_D(ym)[0];

                umul_ppmm(hi, lo, x0, y0);

                lo = (lo >> shift) | (hi << (FLINT_BITS - shift));
                hi = (hi >> shift);

                if (xnegative ^ ynegative)
                    sub_ddmmss(sum[1], sum[0], sum[1], sum[0], hi, lo);
                else
                    add_ssaaaa(sum[1], sum[0], sum[1], sum[0], hi, lo);
            }
            else if (xn == 2 && yn == 2 && shift < FLINT_BITS && sn <= 3)
            {
                mp_limb_t x1, x0, y1, y0;
                mp_limb_t u3, u2, u1, u0;

                x0 = ARF_NOPTR_D(xm)[0];
                x1 = ARF_NOPTR_D(xm)[1];
                y0 = ARF_NOPTR_D(ym)[0];
                y1 = ARF_NOPTR_D(ym)[1];

                FLINT_MPN_MUL_2X2(u3, u2, u1, u0, x1, x0, y1, y0);

                u1 = (u1 >> shift) | (u2 << (FLINT_BITS - shift));
                u2 = (u2 >> shift) | (u3 << (FLINT_BITS - shift));
                u3 = (u3 >> shift);

                if (sn == 2)
                {
                    if (xnegative ^ ynegative)
                        sub_ddmmss(sum[1], sum[0], sum[1], sum[0], u3, u2);
                    else
                        add_ssaaaa(sum[1], sum[0], sum[1], sum[0], u3, u2);
                }
                else
                {
                    if (xnegative ^ ynegative)
                        sub_dddmmmsss(sum[2], sum[1], sum[0], sum[2], sum[1], sum[0], u3, u2, u1);
                    else
                        add_sssaaaaaa(sum[2], sum[1], sum[0], sum[2], sum[1], sum[0], u3, u2, u1);
                }
            }
#endif
            else if (xn <= 2 && yn <= 2 && sn <= 3)
            {
                mp_limb_t x1, x0, y1, y0;
                mp_limb_t u3, u2, u1, u0;

                if (xn == 1 && yn == 1)
                {
                    x0 = ARF_NOPTR_D(xm)[0];
                    y0 = ARF_NOPTR_D(ym)[0];
                    umul_ppmm(u3, u2, x0, y0);
                    u1 = u0 = 0;
                }
                else if (xn == 2 && yn == 2)
                {
                    x0 = ARF_NOPTR_D(xm)[0];
                    x1 = ARF_NOPTR_D(xm)[1];
                    y0 = ARF_NOPTR_D(ym)[0];
                    y1 = ARF_NOPTR_D(ym)[1];
                    FLINT_MPN_MUL_2X2(u3, u2, u1, u0, x1, x0, y1, y0);
                }
                else if (xn == 1)
                {
                    x0 = ARF_NOPTR_D(xm)[0];
                    y0 = ARF_NOPTR_D(ym)[0];
                    y1 = ARF_NOPTR_D(ym)[1];
                    FLINT_MPN_MUL_2X1(u3, u2, u1, y1, y0, x0);
                    u0 = 0;
                }
                else
                {
                    x0 = ARF_NOPTR_D(xm)[0];
                    x1 = ARF_NOPTR_D(xm)[1];
                    y0 = ARF_NOPTR_D(ym)[0];
                    FLINT_MPN_MUL_2X1(u3, u2, u1, x1, x0, y0);
                    u0 = 0;
                }

                if (sn == 2)
                {
                    if (shift < FLINT_BITS)
                    {
                        u2 = (u2 >> shift) | (u3 << (FLINT_BITS - shift));
                        u3 = (u3 >> shift);
                    }
                    else if (shift == FLINT_BITS)
                    {
                        u2 = u3;
                        u3 = 0;
                    }
                    else /* FLINT_BITS < shift < 2 * FLINT_BITS */
                    {
                        u2 = (u3 >> (shift - FLINT_BITS));
                        u3 = 0;
                    }

                    if (xnegative ^ ynegative)
                        sub_ddmmss(sum[1], sum[0], sum[1], sum[0], u3, u2);
                    else
                        add_ssaaaa(sum[1], sum[0], sum[1], sum[0], u3, u2);
                }
                else if (sn == 3)
                {
                    if (shift < FLINT_BITS)
                    {
                        u1 = (u1 >> shift) | (u2 << (FLINT_BITS - shift));
                        u2 = (u2 >> shift) | (u3 << (FLINT_BITS - shift));
                        u3 = (u3 >> shift);
                    }
                    else if (shift == FLINT_BITS)
                    {
                        u1 = u2;
                        u2 = u3;
                        u3 = 0;
                    }
                    else if (shift < 2 * FLINT_BITS)
                    {
                        u1 = (u3 << (2 * FLINT_BITS - shift)) | (u2 >> (shift - FLINT_BITS));
                        u2 = (u3 >> (shift - FLINT_BITS));
                        u3 = 0;
                    }
                    else if (shift == 2 * FLINT_BITS)
                    {
                        u1 = u3;
                        u2 = 0;
                        u3 = 0;
                    }
                    else  /* 2 * FLINT_BITS < shift < 3 * FLINT_BITS */
                    {
                        u1 = (u3 >> (shift - 2 * FLINT_BITS));
                        u2 = 0;
                        u3 = 0;
                    }

                    if (xnegative ^ ynegative)
                        sub_dddmmmsss(sum[2], sum[1], sum[0], sum[2], sum[1], sum[0], u3, u2, u1);
                    else
                        add_sssaaaaaa(sum[2], sum[1], sum[0], sum[2], sum[1], sum[0], u3, u2, u1);
                }
            }
            else
            {
                mp_srcptr xptr, yptr;

                xptr = (xn <= ARF_NOPTR_LIMBS) ? ARF_NOPTR_D(xm) : ARF_PTR_D(xm);
                yptr = (yn <= ARF_NOPTR_LIMBS) ? ARF_NOPTR_D(ym) : ARF_PTR_D(ym);

                _arb_dot_addmul_generic(sum, &serr, tmp, sn, xptr, xn, yptr, yn, xnegative ^ ynegative, shift);
            }
        }
    }

    xnegative = 0;
    if (sum[sn - 1] >= LIMB_TOP)
    {
        mpn_neg(sum, sum, sn);
        xnegative = 1;
    }

    if (sum[sn - 1] == 0)
    {
        slong sum_exp2;
        mp_size_t sn2;

        sn2 = sn;
        sum_exp2 = sum_exp;

        while (sn2 > 0 && sum[sn2 - 1] == 0)
        {
            sum_exp2 -= FLINT_BITS;
            sn2--;
        }

        if (sn2 == 0)
        {
            arf_zero(arb_midref(res));
        }
        else
        {
            _arf_set_round_mpn(arb_midref(res), &exp, sum, sn2, xnegative ^ subtract, prec, ARF_RND_DOWN);
            _fmpz_set_si_small(ARF_EXPREF(arb_midref(res)), exp + sum_exp2);
        }
    }
    else
    {

        if (sn == 2)
            _arf_set_round_uiui(arb_midref(res), &exp, sum[1], sum[0], xnegative ^ subtract, prec, ARF_RND_DOWN);
        else
            _arf_set_round_mpn(arb_midref(res), &exp, sum, sn, xnegative ^ subtract, prec, ARF_RND_DOWN);

        _fmpz_set_si_small(ARF_EXPREF(arb_midref(res)), exp + sum_exp);
    }

    ARF_ADD_TMP_FREE(sum, alloc);
}
