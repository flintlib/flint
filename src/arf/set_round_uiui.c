/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arf.h"

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
            mp_limb_t hi_mask, lo_mask, rndn_mask, __t, __u; \
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
_arf_set_round_uiui(arf_t z, slong * fix, mp_limb_t hi, mp_limb_t lo, int sgnbit, slong prec, arf_rnd_t rnd)
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
                    mp_limb_t t, ovf;
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

