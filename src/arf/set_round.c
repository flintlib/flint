/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "arf.h"

int
arf_set_round(arf_t y, const arf_t x, slong prec, arf_rnd_t rnd)
{
    if (arf_is_special(x))
    {
        arf_set(y, x);
        return 0;
    }
    else
    {
        int inexact;
        slong fix;
        slong xn;
        nn_srcptr xptr;

        if (y == x)
        {
            nn_ptr xtmp;
            TMP_INIT;

            ARF_GET_MPN_READONLY(xptr, xn, x);

            /* exact */
            if (xn * FLINT_BITS <= prec)
                return 0;

            if ((xn - 1) * FLINT_BITS < prec)
            {
                /* exact */
                if ((xptr[0] << (prec - (xn-1) *  FLINT_BITS)) == 0)
                    return 0;
            }

            /* inexact */
            TMP_START;
            xtmp = TMP_ALLOC(xn * sizeof(ulong));
            flint_mpn_copyi(xtmp, xptr, xn);
            inexact = _arf_set_round_mpn(y, &fix, xtmp, xn, ARF_SGNBIT(x), prec, rnd);
            _fmpz_add_fast(ARF_EXPREF(y), ARF_EXPREF(x), fix);

            TMP_END;
            return inexact;
        }
        else
        {
            ARF_GET_MPN_READONLY(xptr, xn, x);
            inexact = _arf_set_round_mpn(y, &fix, xptr, xn,
                ARF_SGNBIT(x), prec, rnd);
            _fmpz_add_fast(ARF_EXPREF(y), ARF_EXPREF(x), fix);
            return inexact;
        }
    }
}

int
_arf_set_round_mpn(arf_t y, slong * exp_shift, nn_srcptr x, slong xn,
    int sgnbit, slong prec, arf_rnd_t rnd)
{
    unsigned int leading;
    flint_bitcnt_t exp, bc, val, val_bits;
    slong yn, val_limbs;
    nn_ptr yptr;
    ulong t;
    int increment, inexact;

    /* Compute the total bit length of x. */
    leading = flint_clz(x[xn - 1]);
    exp = xn * FLINT_BITS - leading;

    /* Set exponent. */
    *exp_shift = -(slong) leading;

    /* Find first nonzero bit. */
    val_limbs = 0;
    while (x[val_limbs] == 0)
        val_limbs++;
    val_bits = flint_ctz(x[val_limbs]);
    val = val_limbs * FLINT_BITS + val_bits;

    if (exp - val <= (ulong) prec)
    {
        inexact = 0;
        increment = 0;
    }
    else
    {
        inexact = 1;

        /* Limb and bit of the truncation point. */
        val_limbs = (exp - prec) / FLINT_BITS;
        val_bits = (exp - prec) % FLINT_BITS;

        if (rnd == ARF_RND_DOWN)
        {
            increment = 0;
        }
        else if (rnd == ARF_RND_NEAR)
        {
            /* If exactly one excess bit, there is a tie; the rounding
               direction is determined by the bit to the left of the
               truncation point. */
            if (exp - val - 1 == (ulong) prec)
            {
                increment = (x[val_limbs] >> val_bits) & 1;
            }
            else
            {
                /* The bit to the right of the truncation point determines
                   the rounding direction. */
                slong exc_limbs = (exp - prec - 1) / FLINT_BITS;
                flint_bitcnt_t exc_bits = (exp - prec - 1) % FLINT_BITS;

                increment = (x[exc_limbs] >> exc_bits) & 1;
            }
        }
        else
        {
            if (rnd == ARF_RND_UP)
                increment = 1;
            else if (rnd == ARF_RND_FLOOR)
                increment = sgnbit;
            else
                increment = !sgnbit;
        }

        if (!increment)
        {
            /* Find first nonzero bit from the truncation point. */
            t = x[val_limbs] & (LIMB_ONES << val_bits);
            while (t == 0)
            {
                val_limbs++;
                t = x[val_limbs];
            }

            val_bits = flint_ctz(t);
            val = val_limbs * FLINT_BITS + val_bits;
        }
        else
        {
            /* Find first zero bit from the truncation point */
            t = (~x[val_limbs]) & (LIMB_ONES << val_bits);
            while (t == 0)
            {
                val_limbs++;
                if (val_limbs < xn)
                    t = ~x[val_limbs];
                else  /* The array is all ones up to the highest limb. */
                {
                    val_bits = 0;
                    goto END_SCAN1;
                }
            }
            val_bits = flint_ctz(t);
            END_SCAN1:
            val = val_limbs * FLINT_BITS + val_bits;

            /* Overflow to next power of two (unlikely). */
            if (val == exp)
            {
                exp_shift[0]++;
                ARF_DEMOTE(y);
                ARF_NOPTR_D(y)[0] = LIMB_TOP;
                ARF_XSIZE(y) = ARF_MAKE_XSIZE(1, sgnbit);
                return 1;
            }
        }
    }

    /* Now copy the result to destination. */
    x += val_limbs;
    xn -= val_limbs;
    bc = exp - val;
    yn = (bc + FLINT_BITS - 1) / FLINT_BITS;

    ARF_GET_MPN_WRITE(yptr, yn, y);
    ARF_XSIZE(y) |= sgnbit;

    if (leading == 0)
    {
        flint_mpn_copyi(yptr, x, xn);
    }
    else if (xn == yn)
    {
        mpn_lshift(yptr, x, yn, leading);
    }
    else
    {
        mpn_lshift(yptr, x + 1, yn, leading);
        yptr[0] |= (x[0] >> (FLINT_BITS - leading));
    }

    if (increment)
    {
        /* Mask off bits from the last limb. */
        yptr[0] &= LIMB_ONES << (yn * FLINT_BITS - bc);

        /* Increment (no carry propagation). */
        yptr[0] += LIMB_ONE << (yn * FLINT_BITS - bc);
    }
    else if (inexact && prec < yn * FLINT_BITS)
    {
        /* Mask off bits from the last limb. */
        yptr[0] &= LIMB_ONES << (yn * FLINT_BITS - prec);
    }

    return inexact;
}

/* Top-aligns the single-limb integer value v and rounds it to prec bits.
   Writes inexact, v, exp. Warning: macro without parentheses. */
#define ARF_NORMALISE_ROUND_LIMB(inexact, exp, v, sgnbit, prec, rnd) \
    do { \
        exp = flint_clz(v); \
        v <<= exp; \
        exp = FLINT_BITS - exp; \
        if (prec >= exp) \
        { \
            inexact = 0; \
        } \
        else \
        { \
            ulong hi_mask, lo_mask, rndn_mask, __t, __u; \
            hi_mask = LIMB_ONES << (FLINT_BITS - prec); \
            __t = v & hi_mask; \
            inexact = (__t != v); \
            if (inexact && rnd != ARF_RND_DOWN) \
            { \
                if (rnd == ARF_RND_NEAR) \
                { \
                    lo_mask = LIMB_ONES >> prec; \
                    rndn_mask = LIMB_ONE << (FLINT_BITS - prec - 1); \
                    __u = v & lo_mask; \
                    if (__u > rndn_mask || (__u == rndn_mask && \
                            (__t << (prec - 1)))) \
                        __t += (LIMB_ONE << (FLINT_BITS - prec)); \
                } \
                else if (arf_rounds_up(rnd, sgnbit)) \
                { \
                    __t += (LIMB_ONE << (FLINT_BITS - prec)); \
                } \
                if (__t == 0) \
                { \
                    __t = LIMB_TOP; \
                    exp++; \
                } \
            } \
            v = __t; \
        } \
    } while (0)

int
_arf_set_round_ui(arf_t x, ulong v, int sgnbit, slong prec, arf_rnd_t rnd)
{
    _fmpz_demote(ARF_EXPREF(x));
    ARF_DEMOTE(x);

    if (v == 0)
    {
        ARF_EXP(x) = ARF_EXP_ZERO;
        ARF_XSIZE(x) = 0;
        return 0;
    }
    else
    {
        int exp, inexact;
        ARF_NORMALISE_ROUND_LIMB(inexact, exp, v, sgnbit, prec, rnd);
        ARF_EXP(x) = exp;
        ARF_XSIZE(x) = ARF_MAKE_XSIZE(1, sgnbit);
        ARF_NOPTR_D(x)[0] = v;
        return inexact;
    }
}

int
_arf_set_round_uiui(arf_t z, slong * fix, ulong hi, ulong lo, int sgnbit, slong prec, arf_rnd_t rnd)
{
    int leading, trailing, bc, inexact, zn, up, exp;

    if (hi == 0)
    {
        ARF_NORMALISE_ROUND_LIMB(inexact, exp, lo, sgnbit, prec, rnd);
        leading = 2 * FLINT_BITS - exp;
        zn = 1;
    }
    else if (lo == 0)
    {
        ARF_NORMALISE_ROUND_LIMB(inexact, exp, hi, sgnbit, prec, rnd);
        leading = FLINT_BITS - exp;
        lo = hi;
        zn = 1;
    }
    else
    {
        leading = flint_clz(hi);
        trailing = flint_ctz(lo);

        bc = 2 * FLINT_BITS - leading - trailing;

        if (bc <= prec)
        {
            inexact = 0;
            zn = 2;

            if (leading != 0)
            {
                if (bc <= FLINT_BITS)
                {
                    lo = (hi << leading) | (lo >> (FLINT_BITS - leading));
                    zn = 1;
                }
                else
                {
                    hi = (hi << leading) | (lo >> (FLINT_BITS - leading));
                    lo = (lo << leading);
                }
            }
        }
        else
        {
            inexact = 1;

            if (rnd == ARF_RND_DOWN)
            {
                up = 0;
            }
            else if (rnd == ARF_RND_NEAR)
            {
                if (bc == prec + 1)
                {
                    /* exactly one excess bit; check the parity bit which
                       must be either the lsb of hi or a bit in lo */
                    if (trailing == FLINT_BITS - 1)
                        up = hi & 1;
                    else
                        up = (lo >> (trailing + 1)) & 1;
                }
                else
                {
                    /* two or more excess bits; test the first excess bit */
                    flint_bitcnt_t pos = 2 * FLINT_BITS - leading - prec - 1;

                    if (pos < FLINT_BITS)
                        up = (lo >> pos) & 1;
                    else
                        up = (hi >> (pos - FLINT_BITS)) & 1;
                }
            }
            else
            {
                up = arf_rounds_up(rnd, sgnbit);
            }

            if (prec <= FLINT_BITS)
            {
                zn = 1;

                if (leading == 0)
                    lo = hi;
                else
                    lo = (hi << leading) | (lo >> (FLINT_BITS - leading));

                lo = lo & (LIMB_ONES << (FLINT_BITS - prec));

                if (up)
                {
                    ulong t, ovf;
                    t = lo + (LIMB_ONE << (FLINT_BITS - prec));
                    ovf = (t == 0);
                    leading -= ovf;
                    lo = (t >> ovf) | (ovf << (FLINT_BITS - 1));
                }
            }
            else
            {
                zn = 2;

                if (leading != 0)
                {
                    hi = (hi << leading) | (lo >> (FLINT_BITS - leading));
                    lo = (lo << leading);
                }

                lo = lo & (LIMB_ONES << (2 * FLINT_BITS - prec));

                if (up)
                {
                    add_ssaaaa(hi, lo, hi, lo, 0, (LIMB_ONE << (2 * FLINT_BITS - prec)));
                }

                if (lo == 0)
                {
                    if (hi == 0)
                    {
                        leading -= 1;
                        lo = LIMB_TOP;
                    }
                    else
                    {
                        lo = hi;
                    }

                    zn = 1;
                }
            }
        }
    }

    *fix = -leading;
    ARF_DEMOTE(z);
    ARF_XSIZE(z) = ARF_MAKE_XSIZE(zn, sgnbit);
    ARF_NOPTR_D(z)[0] = lo;
    ARF_NOPTR_D(z)[1] = hi;
    return inexact;
}
