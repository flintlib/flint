/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arf.h"

int
_arf_set_round_mpn(arf_t y, slong * exp_shift, mp_srcptr x, mp_size_t xn,
    int sgnbit, slong prec, arf_rnd_t rnd)
{
    unsigned int leading;
    flint_bitcnt_t exp, bc, val, val_bits;
    mp_size_t yn, val_limbs;
    mp_ptr yptr;
    mp_limb_t t;
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

    if (exp - val <= prec)
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
            if (exp - val - 1 == prec)
            {
                increment = (x[val_limbs] >> val_bits) & 1;
            }
            else
            {
                /* The bit to the right of the truncation point determines
                   the rounding direction. */
                mp_size_t exc_limbs = (exp - prec - 1) / FLINT_BITS;
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

